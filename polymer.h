/* 
 * This code is part of the project "NUMA-aware Graph-structured Analytics"
 * 
 *
 * Copyright (C) 2014 Institute of Parallel And Distributed Systems (IPADS), Shanghai Jiao Tong University
 *     All rights reserved 
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 * 
 * For more about this software, visit:
 *
 *     http://ipads.se.sjtu.edu.cn/projects/polymer.html
 *
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <sys/mman.h>

#include "custom-barrier.h"
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "IO-numa.h"

#include <numa.h>
#include <pthread.h>

using namespace std;

#define PAGESIZE (4096)

#define MIN(x, y) ((x > y) ? (y) : (x))

//*****START FRAMEWORK*****

/* This is no longer ligra, 
the vertices set should be set by each node
*/

inline int SXCHG(char *ptr, char newv) {
    char ret = newv;
    __asm__ __volatile__ (
    "  xchgb %0,%1\n"
    : "=r" (ret)
    : "m" (*(volatile char *) ptr), "0" (newv)
    : "memory");
    return ret;
}

int roundUp(double x) {
    int ones = x / 1;
    double others = x - ones;
    if (others > 0) {
        ones++;
    }
    return ones;
}

struct Subworker_Partitioner {
    int tid;
    int subTid;
    int numOfSub;
    int dense_start;
    int dense_end;
    pthread_barrier_t *global_barr;
    pthread_barrier_t *local_barr;
    pthread_barrier_t *leader_barr;
    Custom_barrier local_custom;
    Custom_barrier subMaster_custom;

    Subworker_Partitioner(int nSub) : numOfSub(nSub) {}

    inline bool isMaster() { return (tid + subTid == 0); }

    inline bool isSubMaster() { return (subTid == 0); }

    inline intT getStartPos(intT m) { return subTid * (m / numOfSub); }

    inline intT getEndPos(intT m) { return (subTid == numOfSub - 1) ? m : ((subTid + 1) * (m / numOfSub)); }

    inline void localWait() {
        local_custom.wait();
    }

    inline void globalWait() {
        local_custom.wait();
        if (isSubMaster()) {
            subMaster_custom.wait();
        }
        local_custom.wait();
    }
};

struct Default_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;

    Default_Hash_F(int _n, int _shardNum) : n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum) {}

    inline int hashFunc(int index) {
        if (index >= shardNum * vertPerShard) {
            return index;
        }
        int idxOfShard = index % shardNum;
        int idxInShard = index / shardNum;
        return (idxOfShard * vertPerShard + idxInShard);
    }

    inline int hashBackFunc(int index) {
        if (index >= shardNum * vertPerShard) {
            return index;
        }
        int idxOfShard = index / vertPerShard;
        int idxInShard = index % vertPerShard;
        return (idxOfShard + idxInShard * shardNum);
    }
};

template<class vertex>
void partitionByDegree(graph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, bool useOutDegree = false) {
    const intT n = GA.n;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useOutDegree) {
        { parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree(); }
    } else {
        { parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree(); }
    }

    int accum[numOfShards];
    for (int i = 0; i < numOfShards; i++) {
        accum[i] = 0;
        sizeArr[i] = 0;
    }

    long totalDegree = 0;
    for (intT i = 0; i < n; i++) {
        totalDegree += degrees[i];
    }

    int averageDegree = totalDegree / numOfShards;
    cerr << "average is " + to_string(averageDegree) + "\n";
    int counter = 0;
    int tmpSizeCounter = 0;
    for (intT i = 0; i < n; i += PAGESIZE / sizeOfOneEle) {
        int localAccum = 0;
        int localSize = 0;
        for (intT j = 0; j < PAGESIZE / sizeOfOneEle; j++) {
            if (i + j >= n)
                break;
            //accum[counter] += degrees[i + j];
            //sizeArr[counter]++;
            localAccum += degrees[i + j];
            localSize++;
            tmpSizeCounter++;
        }
        accum[counter] += localAccum;
        sizeArr[counter] += localSize;
        if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
            int oldDiff = averageDegree - (accum[counter] - localAccum);
            int newDiff = accum[counter] - averageDegree;
            if (oldDiff < newDiff) {
                accum[counter] -= localAccum;
                sizeArr[counter] -= localSize;
            } else {
                localAccum = 0;
                localSize = 0;
            }
            counter++;
            accum[counter] += localAccum;
            sizeArr[counter] += localSize;
            //cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle) << endl;
            tmpSizeCounter = 0;
        }
    }

    free(degrees);
}

template<class vertex>
void subPartitionByDegree(graph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, bool useOutDegree = false,
                          bool useFakeDegree = false) {
    const intT n = GA.n;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useFakeDegree) {
        { parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getFakeDegree(); }
    } else {
        if (useOutDegree) {
            { parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree(); }
        } else {
            { parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree(); }
        }
    }

    int accum[numOfShards];
    for (int i = 0; i < numOfShards; i++) {
        accum[i] = 0;
        sizeArr[i] = 0;
    }

    long totalDegree = 0;
    for (intT i = 0; i < n; i++) {
        totalDegree += degrees[i];
    }

    int averageDegree = totalDegree / numOfShards;
    int counter = 0;
    int tmpSizeCounter = 0;
    for (intT i = 0; i < n; i++) {
        accum[counter] += degrees[i];
        sizeArr[counter]++;
        tmpSizeCounter++;
        if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
            counter++;
            //cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle) << endl;
            tmpSizeCounter = 0;
        }
    }

    free(degrees);
}

template<class vertex>
void subPartitionByDegree(graph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, int subStart, int subEnd,
                          bool useOutDegree = false, bool useFakeDegree = false) {
    const intT n = subEnd - subStart;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useFakeDegree) {
        { parallel_for (intT i = subStart; i < subEnd; i++) degrees[i - subStart] = GA.V[i].getFakeDegree(); }
    } else {
        if (useOutDegree) {
            { parallel_for (intT i = subStart; i < subEnd; i++) degrees[i - subStart] = GA.V[i].getOutDegree(); }
        } else {
            { parallel_for (intT i = subStart; i < subEnd; i++) degrees[i - subStart] = GA.V[i].getInDegree(); }
        }
    }

    int accum[numOfShards];
    for (int i = 0; i < numOfShards; i++) {
        accum[i] = 0;
        sizeArr[i] = 0;
    }

    long totalDegree = 0;
    for (intT i = 0; i < n; i++) {
        totalDegree += degrees[i];
    }

    int averageDegree = totalDegree / numOfShards;
    int counter = 0;
    int tmpSizeCounter = 0;
    for (intT i = 0; i < n; i++) {
        accum[counter] += degrees[i];
        sizeArr[counter]++;
        tmpSizeCounter++;
        if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
            counter++;
            //cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle) << endl;
            tmpSizeCounter = 0;
        }
    }

    free(degrees);
}

template<class vertex, class Hash_F>
void graphHasher(graph<vertex> &GA, Hash_F hash) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *) malloc(sizeof(vertex) * GA.n);

    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intT d = V[i].getOutDegree();
            V[i].setFakeDegree(d);
            intE *outEdges = V[i].getOutNeighborPtr();
            for (intT j = 0; j < d; j++) {
                outEdges[j] = hash.hashFunc(outEdges[j]);
            }
            newVertexSet[hash.hashFunc(i)] = V[i];
        }
    }
    GA.V = newVertexSet;
    free(V);
}

template<class vertex, class Hash_F>
void graphInEdgeHasher(graph<vertex> &GA, Hash_F hash) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *) malloc(sizeof(vertex) * GA.n);

    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intT d = V[i].getInDegree();
            //V[i].setFakeDegree(d);
            intE *inEdges = V[i].getInNeighborPtr();
            for (intT j = 0; j < d; j++) {
                inEdges[j] = hash.hashFunc(inEdges[j]);
            }
            newVertexSet[hash.hashFunc(i)] = V[i];
        }
    }
    GA.V = newVertexSet;
    free(V);
}

template<class vertex, class Hash_F>
void graphAllEdgeHasher(graph<vertex> &GA, Hash_F hash) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *) malloc(sizeof(vertex) * GA.n);

    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intT d = V[i].getOutDegree();
            V[i].setFakeDegree(d);
            intE *outEdges = V[i].getOutNeighborPtr();
            for (intT j = 0; j < d; j++) {
                outEdges[j] = hash.hashFunc(outEdges[j]);
            }
            d = V[i].getInDegree();
            intE *inEdges = V[i].getInNeighborPtr();
            for (intT j = 0; j < d; j++) {
                inEdges[j] = hash.hashFunc(inEdges[j]);
            }
            newVertexSet[hash.hashFunc(i)] = V[i];
        }
    }
    GA.V = newVertexSet;
    free(V);
}

template<class vertex>
graph<vertex> graphFilter(graph<vertex> &GA, int rangeLow, int rangeHi, bool useOutEdge = true) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *) numa_alloc_local(sizeof(vertex) * GA.n);
    int *counters = (int *) numa_alloc_local(sizeof(int) * GA.n);
    int *offsets = (int *) numa_alloc_local(sizeof(int) * GA.n);
    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intT d = (useOutEdge) ? (V[i].getOutDegree()) : (V[i].getInDegree());
            //V[i].setFakeDegree(d);
            newVertexSet[i].setOutDegree(V[i].getOutDegree());
            newVertexSet[i].setInDegree(V[i].getInDegree());

            counters[i] = 0;
            for (intT j = 0; j < d; j++) {
                intT ngh = (useOutEdge) ? (V[i].getOutNeighbor(j)) : (V[i].getInNeighbor(j));
                if (rangeLow <= ngh && ngh < rangeHi)
                    counters[i]++;
            }
            newVertexSet[i].setFakeDegree(counters[i]);
        }
    }

    intT totalSize = 0;
    for (intT i = 0; i < GA.n; i++) {
        offsets[i] = totalSize;
        totalSize += counters[i];
    }

    numa_free(counters, sizeof(int) * GA.n);

    intE *edges = (intE *) numa_alloc_local(sizeof(intE) * totalSize);

    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intE *localEdges = &edges[offsets[i]];
            intT counter = 0;
            intT d = (useOutEdge) ? (V[i].getOutDegree()) : (V[i].getInDegree());
            for (intT j = 0; j < d; j++) {
                intT ngh = (useOutEdge) ? (V[i].getOutNeighbor(j)) : (V[i].getInNeighbor(j));
                if (rangeLow <= ngh && ngh < rangeHi) {
                    localEdges[counter] = ngh;
                    counter++;
                }
            }
            if (counter != newVertexSet[i].getFakeDegree()) {
                printf("oops: %d %d\n", counter, newVertexSet[i].getFakeDegree());
            }
            if (i == 0) {
                cerr << "fake deg: " + to_string(newVertexSet[i].getFakeDegree()) + "\n";
            }
            if (useOutEdge)
                newVertexSet[i].setOutNeighbors(localEdges);
            else
                newVertexSet[i].setInNeighbors(localEdges);
        }
    }
    numa_free(offsets, sizeof(int) * GA.n);
    //printf("degree: %d\n", newVertexSet[0].getFakeDegree());
    return graph<vertex>(newVertexSet, GA.n, GA.m);
}

template<class vertex>
graph<vertex> graphFilter2Direction(graph<vertex> &GA, int rangeLow, int rangeHi) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *) numa_alloc_local(sizeof(vertex) * GA.n);
    int *counters = (int *) numa_alloc_local(sizeof(int) * GA.n);
    int *offsets = (int *) numa_alloc_local(sizeof(int) * GA.n);
    int *inCounters = (int *) numa_alloc_local(sizeof(int) * GA.n);
    int *inOffsets = (int *) numa_alloc_local(sizeof(int) * GA.n);
    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            newVertexSet[i].setOutDegree(V[i].getOutDegree());
            newVertexSet[i].setInDegree(V[i].getInDegree());

            intT d = V[i].getOutDegree();
            counters[i] = 0;
            for (intT j = 0; j < d; j++) {
                intT ngh = V[i].getOutNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi)
                    counters[i]++;
            }
            newVertexSet[i].setFakeDegree(counters[i]);

            d = V[i].getInDegree();
            inCounters[i] = 0;
            for (intT j = 0; j < d; j++) {
                intT ngh = V[i].getInNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi)
                    inCounters[i]++;
            }
            newVertexSet[i].setFakeDegree(counters[i]);
            newVertexSet[i].setFakeInDegree(inCounters[i]);
        }
    }

    intT totalSize = 0;
    intT totalInSize = 0;
    for (intT i = 0; i < GA.n; i++) {
        offsets[i] = totalSize;
        totalSize += counters[i];

        inOffsets[i] = totalInSize;
        totalInSize += inCounters[i];
    }

    numa_free(counters, sizeof(int) * GA.n);
    numa_free(inCounters, sizeof(int) * GA.n);

    intE *edges = (intE *) numa_alloc_local(sizeof(intE) * totalSize);
    intE *inEdges = (intE *) numa_alloc_local(sizeof(intE) * totalInSize);
    cerr << "totalInSize is " + to_string(totalInSize) + "\n";

    {
        parallel_for (intT i = 0; i < GA.n; i++) {
            intE *localEdges = &edges[offsets[i]];
            intT counter = 0;
            intT d = V[i].getOutDegree();
            for (intT j = 0; j < d; j++) {
                intT ngh = V[i].getOutNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi) {
                    localEdges[counter] = ngh;
                    counter++;
                }
            }
            if (counter != newVertexSet[i].getFakeDegree()) {
                printf("oops: %d %d\n", counter, newVertexSet[i].getFakeDegree());
            }

            intE *localInEdges = &inEdges[inOffsets[i]];
            counter = 0;
            d = V[i].getInDegree();
            for (intT j = 0; j < d; j++) {
                intT ngh = V[i].getInNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi) {
                    localInEdges[counter] = ngh;
                    counter++;
                }
            }
            if (counter != newVertexSet[i].getFakeInDegree()) {
                printf("oops: %d %d\n", counter, newVertexSet[i].getFakeInDegree());
            }

            if (i == 0) {
                cerr << "fake deg: " + to_string(newVertexSet[i].getFakeDegree()) + "\n";
            }

            newVertexSet[i].setOutNeighbors(localEdges);
            newVertexSet[i].setInNeighbors(localInEdges);
        }
    }
    numa_free(offsets, sizeof(int) * GA.n);
    numa_free(inOffsets, sizeof(int) * GA.n);
    //printf("degree: %d\n", newVertexSet[0].getFakeDegree());
    return graph<vertex>(newVertexSet, GA.n, GA.m);
}

void *mapDataArray(int numOfShards, int *sizeArr, int sizeOfOneEle) {
    int numOfPages = 0;
    for (int i = 0; i < numOfShards; i++) {
        numOfPages += sizeArr[i] / (double) (PAGESIZE / sizeOfOneEle);
    }
    numOfPages++;

    void *toBeReturned = mmap(NULL, numOfPages * PAGESIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (toBeReturned == NULL) {
        cout << "OOps" << endl;
    }

    int offset = 0;
    for (int i = 0; i < numOfShards; i++) {
        void *startPos = (void *) ((char *) toBeReturned + offset * sizeOfOneEle);
        //printf("start binding %d : %d\n", i, offset);
        numa_tonode_memory(startPos, sizeArr[i], i);
        offset = offset + sizeArr[i];
    }
    return toBeReturned;
}

struct AsyncChunk {
    int accessCounter;
    intT m;
    intT *s;
};

struct LocalFrontier {
    intT n;
    intT m;
    intT head;
    intT tail;
    intT insertTail;
    intT emptySignal;
    intT outEdgesCount;
    int startID;
    int endID;
    bool *b;
    intT *s;
    intT sparseCounter;
    intT **sparseChunks;
    intT *chunkSizes;
    intT *tmp;
    AsyncChunk **localQueue;
    bool isDense;

    LocalFrontier(bool *_b, int start, int end) : b(_b), startID(start), endID(end), n(end - start), m(0),
                                                  isDense(true), s(NULL), outEdgesCount(0), sparseChunks(NULL),
                                                  chunkSizes(NULL) {}

    bool inRange(int index) { return (startID <= index && index < endID); }

    inline void setBit(int index, bool val) { b[index - startID] = val; }

    inline bool getBit(int index) { return b[index - startID]; }

    void toSparse() {
        if (isDense) {
            if (s != NULL)
                free(s);
            _seq<intT> R = sequence::packIndex(b, n);
            s = R.A;
            m = R.n;
            { parallel_for (intT i = 0; i < m; i++) s[i] = s[i] + startID; }
            if (m == 0) {
                printf("%p\n", s);
            } else {
                printf("M is %d and first ele is %d\n", m, s[0]);
            }
        }
        isDense = false;
    }

    void toSparseAsync(int nextID, LocalFrontier *next) {
        if (isDense) {
            if (s != NULL)
                free(s);
            _seq<intT> R = sequence::packIndex(b, n);
            s = R.A;
            m = R.n;
            { parallel_for (intT i = 0; i < m; i++) s[i] = s[i] + startID; }
            if (m == 0) {
                printf("%p\n", s);
            } else {
                printf("M is %d and first ele is %d\n", m, s[0]);
            }
            AsyncChunk *myChunk = (AsyncChunk *) malloc(sizeof(AsyncChunk));
            myChunk->s = R.A;
            myChunk->m = R.n;
            myChunk->accessCounter = 0;
            next->localQueue[0] = myChunk;
            next->insertTail = 1;
            next->head = 0;
            next->tail = 1;
        }
        isDense = false;
    }

    void toDense() {
        if (!isDense) {
            { parallel_for (intT i = 0; i < n; i++) b[i] = false; }
            { parallel_for (intT i = 0; i < m; i++) b[s[i] - startID] = true; }
        }
        //printf("hehe\n");
        isDense = true;
    }

    void toDenseWithMerge(int numOfSub) {
        if (!isDense) {
            { parallel_for (intT i = 0; i < n; i++) b[i] = false; }
            for (int i = 0; i < numOfSub; i++) {
                intT *sparsePtr = sparseChunks[i];
                intT size = chunkSizes[i];
                { parallel_for (intT j = 0; j < size; j++) b[sparsePtr[j] - startID] = true; }
            }
        }
        //printf("hehe\n");
        isDense = true;
    }

    void setSparse(intT _m, intT *_s) {
        if (s != NULL) {
            free(s);
        }
        m = _m;
        s = _s;
        isDense = false;
    }

    bool *swapBitVector(bool *newB) {
        bool *tmp = b;
        b = newB;
        return tmp;
    }

    void clearFrontier() {
        if (s != NULL) {
            free(s);
        }
        s = NULL;
        m = 0;
        outEdgesCount = 0;
    }
};

//*****VERTEX OBJECT*****
struct vertices {
    intT n, m;
    int numOfNodes;
    intT numOfVertices;
    int *numOfVertexOnNode;
    int *offsets;
    int *numOfNonZero;
    bool **d;
    LocalFrontier **frontiers;
    LocalFrontier **nextFrontiers;
    bool isDense;
    bool firstSparse;
    AsyncChunk **asyncQueue;
    int asyncEndSignal;
    intT readerTail;
    intT insertTail;

    vertices(int _numOfNodes) {
        this->numOfNodes = _numOfNodes;
        d = (bool **) malloc(numOfNodes * sizeof(bool *));
        frontiers = (LocalFrontier **) malloc(numOfNodes * sizeof(LocalFrontier *));
        nextFrontiers = (LocalFrontier **) malloc(numOfNodes * sizeof(LocalFrontier *));
        numOfVertexOnNode = (int *) malloc(numOfNodes * sizeof(int));
        offsets = (int *) malloc((numOfNodes + 1) * sizeof(int));
        numOfNonZero = (int *) malloc(numOfNodes * sizeof(int));
        numOfVertices = 0;
        m = -1;
        firstSparse = false;
    }

    /*
    void registerArr(int nodeNum, bool *arr, int size) {
	d[nodeNum] = arr;
	numOfVertexOnNode[nodeNum] = size;
    }
    */
    void registerFrontier(int nodeNum, LocalFrontier *frontier) {
        frontiers[nodeNum] = frontier;
        numOfVertexOnNode[nodeNum] = frontier->n;
        isDense = frontier->isDense;
    }

    void toDense() {
        if (isDense)
            return;
        else
            isDense = true;
        for (int i = 0; i < numOfNodes; i++) {
            //printf("todense %d start m = %d n = %d\n", i, frontiers[i]->s[frontiers[i]->m - 1], frontiers[i]->n);
            frontiers[i]->toDense();
        }
    }

    void toDenseWithMerge(int numOfSub) {
        if (isDense)
            return;
        else
            isDense = true;
        for (int i = 0; i < numOfNodes; i++) {
            //printf("todense %d start m = %d n = %d\n", i, frontiers[i]->s[frontiers[i]->m - 1], frontiers[i]->n);
            frontiers[i]->toDenseWithMerge(numOfSub);
        }
    }

    void toSparse() {
        if (!isDense) {
            firstSparse = false;
            return;
        } else {
            isDense = false;
            firstSparse = true;
        }
        for (int i = 0; i < numOfNodes; i++) {
            if (frontiers[i]->isDense) {
                //printf("real convert\n");
            }
            frontiers[i]->toSparse();
        }
    }

    void toSparseAsync() {
        if (!isDense) {
            firstSparse = false;
            return;
        } else {
            isDense = false;
            firstSparse = true;
        }
        for (int i = 0; i < numOfNodes; i++) {
            if (frontiers[i]->isDense) {
                //printf("real convert\n");
            }
            int nextID = (i + 1) % numOfNodes;
            frontiers[i]->toSparseAsync(nextID, frontiers[nextID]);
        }
    }

    long long getEdgeStat() {
        long long sum = 0;
        for (int i = 0; i < numOfNodes; i++) {
            sum = sum + frontiers[i]->outEdgesCount;
        }
        return sum;
    }

    void calculateOffsets() {
        for (int i = 0; i < numOfNodes; i++) {
            numOfVertices += numOfVertexOnNode[i];
        }
        offsets[0] = 0;
        for (int i = 1; i < numOfNodes; i++) {
            offsets[i] = numOfVertexOnNode[i - 1] + offsets[i - 1];
            //printf("offset of %d: %d\n", i, offsets[i]);
        }
        offsets[numOfNodes] = numOfVertices;
    }

    int getSize(int nodeNum) {
        return numOfVertexOnNode[nodeNum];
    }

    int getSparseSize(int nodeNum) {
        return numOfNonZero[nodeNum];
    }

    void calculateNumOfNonZero(int nodeNum) {
        numOfNonZero[nodeNum] = 0;
        if (false && frontiers[nodeNum]->isDense) {
            numOfNonZero[nodeNum] = sequence::sum(frontiers[nodeNum]->b, numOfVertexOnNode[nodeNum]);
            frontiers[nodeNum]->m = numOfNonZero[nodeNum];
        } else {
            numOfNonZero[nodeNum] = frontiers[nodeNum]->m;
        }
        //printf("non zero count of %d: %d\n", nodeNum, frontiers[nodeNum]->m);
    }

    int numNonzeros() {
        if (m < 0) {
            intT sum = 0;
            for (int i = 0; i < numOfNodes; i++) {
                sum = sum + numOfNonZero[i];
            }
            //printf("num of non zero: %d\n", sum);
            m = sum;
        }
        return m;
    }

    bool isEmpty() {
        if (m < 0) {
            int sum = 0;
            for (int i = 0; i < numOfNodes; i++) {
                sum = sum + numOfNonZero[i];
            }
            m = sum;
        }
        return (m == 0);
    }

    int getNodeNumOfIndex(int index) {
        int result = 0;
        while (result < numOfNodes && offsets[result] <= index) {
            result++;
        }
        return result - 1;
    }

    int getNodeNumOfSparseIndex(int index) {
        int result = 0;
        int accum = 0;
        while (result < numOfNodes && accum <= index) {
            accum += numOfNonZero[result];
            result++;
        }
        return result - 1;
    }

    int getOffset(int nodeNum) {
        return offsets[nodeNum];
    }

    void setBit(int index, bool bit) {
        int accum = 0;
        int i = 0;
        while (index >= accum + numOfVertexOnNode[i]) {
            accum += numOfVertexOnNode[i];
            i++;
        }
        *(frontiers[i]->b + (index - accum)) = bit;
    }

    bool getBit(int index) {
        int accum = 0;
        int i = 0;
        while (index >= accum + numOfVertexOnNode[i]) {
            accum += numOfVertexOnNode[i];
            i++;
        }
        return *(frontiers[i]->b + (index - accum));
    }

    bool *getArr(int nodeNum) {
        return frontiers[nodeNum]->b;
    }

    bool *getNextArr(int nodeNum) {
        if (nextFrontiers[nodeNum] == NULL) return NULL;
        return nextFrontiers[nodeNum]->b;
    }

    intT *getSparseArr(int nodeNum) {
        return frontiers[nodeNum]->s;
    }

    LocalFrontier *getFrontier(int nodeNum) {
        return frontiers[nodeNum];
    }

    LocalFrontier *swapFrontier(int nodeNum, LocalFrontier *newOne) {
        if (nodeNum == 0) {
            isDense = newOne->isDense;
        }
        LocalFrontier *oldOne = frontiers[nodeNum];
        frontiers[nodeNum] = newOne;
        m = -1;
        return oldOne;
    }

    bool eq(vertices &b) {
        return false;
    }

    void print() {
    }

    void del() {
        free(offsets);
        free(numOfVertexOnNode);
        free(d);
    }
};

struct Default_worker_arg {
    void *GA;
    int maxIter;
    int tid;
    int numOfNode;
    int rangeLow;
    int rangeHi;
};

struct Default_subworker_arg {
    void *GA;
    void *Global_G;
    int maxIter;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    pthread_barrier_t *global_barr;
    pthread_barrier_t *node_barr;
    pthread_barrier_t *master_barr;
    LocalFrontier *localFrontier;

    volatile int *local_custom_counter;
    volatile int *local_custom_toggle;
};

struct nonNegF {
    bool operator()(intT a) { return (a >= 0); }
};

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options {
    DENSE, DENSE_PARALLEL, DENSE_FORWARD
};


Subworker_Partitioner dummyPartitioner(1);


//*****EDGE FUNCTIONS*****

template<class F, class vertex>
bool *edgeMapDense(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool parallel = 0,
                   Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    intT size = next->endID - next->startID;
    vertex *G = GA.V;

    if (subworker.isSubMaster()) {
        frontier->nextFrontiers[subworker.tid] = next;
    }

    subworker.globalWait();
    int localOffset = next->startID;
    bool *localBitVec = frontier->getArr(subworker.tid);
    int currNodeNum = 0;
    bool *currBitVector = frontier->getNextArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT startPos = subworker.dense_start;
    intT endPos = subworker.dense_end;

    while (startPos >= nextSwitchPoint) {
        currOffset += frontier->getSize(currNodeNum);
        nextSwitchPoint += frontier->getSize(currNodeNum + 1);
        currNodeNum++;
        currBitVector = frontier->getArr(currNodeNum);
    }

    for (intT i = startPos; i < endPos; i++) {
        //next->setBit(i, false);
        if (f.cond(i)) {
            intT d = G[i].getFakeInDegree();
            for (intT j = 0; j < d; j++) {
                intT ngh = G[i].getInNeighbor(j);
                if (localBitVec[ngh - localOffset] && f.updateAtomic(ngh, i)) {
                    currBitVector[i - currOffset] = true;
                }
                if (!f.cond(i)) break;
                //__builtin_prefetch(f.nextPrefetchAddr(G[i].getInNeighbor(j+3)), 1, 3);
            }
        }
    }
    return NULL;
}

template<class F, class vertex>
bool *
edgeMapDenseForward(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false, int start = 0,
                    int end = 0) {
    intT numVertices = GA.n;
    vertex *G = GA.V;

    int currNodeNum = 0;
    bool *currBitVector = frontier->getArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT m = 0;
    intT outEdgesCount = 0;
    bool *nextB = next->b;

    int startPos = 0;
    int endPos = numVertices;
    if (part) {
        startPos = start;
        endPos = end;
        currNodeNum = frontier->getNodeNumOfIndex(startPos);
        //printf("nodeNum: %d %d\n", currNodeNum, endPos);
        currBitVector = frontier->getArr(currNodeNum);
        nextSwitchPoint = frontier->getOffset(currNodeNum + 1);
        currOffset = frontier->getOffset(currNodeNum);
    }
    for (long i = startPos; i < endPos; i++) {
        if (i == nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            currBitVector = frontier->getArr(currNodeNum);
            //printf("OK\n");
        }
        //printf("edgemap: %p\n", currBitVector);
        m += G[i].getFakeDegree();
        if (currBitVector[i - currOffset]) {
            intT d = G[i].getFakeDegree();
            for (intT j = 0; j < d; j++) {
                uintT ngh = G[i].getOutNeighbor(j);
                if (/*next->inRange(ngh) &&*/ f.cond(ngh) && f.updateAtomic(i, ngh)) {
                    /*
                    if (!next->getBit(ngh)) {
                    m++;
                    outEdgesCount += G[ngh].getOutDegree();
                    }
                    */
                    //int idx = ngh - next->startID;
                    //m += 1 - SXCHG((char *)&(nextB[idx]), 1);
                    next->setBit(ngh, true);
                }
                //__builtin_prefetch(f.nextPrefetchAddr(G[i].getOutNeighbor(j+1)), 1, 0);
            }
        }
        //__builtin_prefetch(f.nextPrefetchAddr(i+1), 0, 3);
        //__builtin_prefetch(G[i+3].getOutNeighborPtr(), 0, 3);
    }
    //writeAdd(&(next->m), m);
    //writeAdd(&(next->outEdgesCount), outEdgesCount);
    //printf("edgeMap: %d %d\n", m, outEdgesCount);

    return NULL;
}

#define DYNAMIC_CHUNK_SIZE (64)

template<class F, class vertex>
bool *edgeMapDenseForwardDynamic(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next,
                                 Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    vertex *G = GA.V;
    if (subworker.isMaster()) {
        printf("we are here\n");
    }
    int currNodeNum = 0;
    bool *currBitVector = frontier->getArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT m = 0;
    intT outEdgesCount = 0;
    bool *nextB = next->b;

    intT *counterPtr = &(next->sparseCounter);

    intT oldStartPos = 0;

    intT startPos = __sync_fetch_and_add(counterPtr, DYNAMIC_CHUNK_SIZE);
    intT endPos = (startPos + DYNAMIC_CHUNK_SIZE > numVertices) ? (numVertices) : (startPos + DYNAMIC_CHUNK_SIZE);
    while (startPos < numVertices) {
        m = 0;
        while (startPos >= nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            currBitVector = frontier->getArr(currNodeNum);
        }
        for (long i = startPos; i < endPos; i++) {
            if (i == nextSwitchPoint) {
                currOffset += frontier->getSize(currNodeNum);
                nextSwitchPoint += frontier->getSize(currNodeNum + 1);
                currNodeNum++;
                currBitVector = frontier->getArr(currNodeNum);
            }
            m += G[i].getFakeDegree();
            if (currBitVector[i - currOffset]) {
                intT d = G[i].getFakeDegree();
                for (intT j = 0; j < d; j++) {
                    uintT ngh = G[i].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(i, ngh)) {
                        next->setBit(ngh, true);
                    }
                }
            }
        }

        if (subworker.tid == 7) {
            //printf("%d: %d to %d deg: %d\n", subworker.subTid, startPos, endPos, m);
        }

        oldStartPos = startPos;
        startPos = __sync_fetch_and_add(counterPtr, DYNAMIC_CHUNK_SIZE);
        endPos = (startPos + DYNAMIC_CHUNK_SIZE > numVertices) ? (numVertices) : (startPos + DYNAMIC_CHUNK_SIZE);
    }
    return NULL;
}

template<class F, class vertex>
bool *edgeMapDenseReduce(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool parallel = 0,
                         Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    intT size = next->endID - next->startID;
    vertex *G = GA.V;

    if (subworker.isSubMaster()) {
        frontier->nextFrontiers[subworker.tid] = next;
    }

    subworker.globalWait();

    struct timeval startT, endT;
    struct timezone tz = {0, 0};
    gettimeofday(&startT, &tz);

    int localOffset = next->startID;
    bool *localBitVec = frontier->getArr(subworker.tid);
    int currNodeNum = 0;
    bool *currBitVector = frontier->getNextArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT startPos = subworker.dense_start;
    intT endPos = subworker.dense_end;

    while (startPos >= nextSwitchPoint) {
        currOffset += frontier->getSize(currNodeNum);
        nextSwitchPoint += frontier->getSize(currNodeNum + 1);
        currNodeNum++;
        currBitVector = frontier->getNextArr(currNodeNum);
    }

    for (intT i = startPos; i < endPos; i++) {
        //next->setBit(i, false);
        if (i >= nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            currBitVector = frontier->getNextArr(currNodeNum);
        }
        if (f.cond(i)) {
            double data[2];
            intT d = G[i].getFakeInDegree();
            f.initFunc((void *) data, i);
            bool shouldActive = false;
            for (intT j = 0; j < d; j++) {
                intT ngh = G[i].getInNeighbor(j);
                if (localBitVec[ngh - localOffset] && f.reduceFunc((void *) data, ngh)) {
                    currBitVector[i - currOffset] = true;
                    //shouldActive = true;
                }
                if (!f.cond(i)) break;
                //__builtin_prefetch(f.nextPrefetchAddr(G[i].getInNeighbor(j+3)), 1, 3);
            }
            if (d > 0) {
                f.combineFunc((void *) data, i);
            }
        }
    }

    subworker.localWait();
    gettimeofday(&endT, &tz);
    if (subworker.isSubMaster()) {
        double time1 = ((double) startT.tv_sec) + ((double) startT.tv_usec) / 1000000.0;
        double time2 = ((double) endT.tv_sec) + ((double) endT.tv_usec) / 1000000.0;
        double duration = time2 - time1;
        //printf("time of %d: %lf\n", subworker.tid, duration);
    }
    return NULL;
}

template<class F, class vertex>
bool *edgeMapDenseDynamic(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next,
                          Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    vertex *G = GA.V;
    if (subworker.isMaster()) {
        //printf("we are here\n");
    }

    if (subworker.isSubMaster()) {
        frontier->nextFrontiers[subworker.tid] = next;
    }

    subworker.globalWait();
    int localOffset = next->startID;
    bool *localBitVec = frontier->getArr(subworker.tid);
    int currNodeNum = 0;
    bool *currBitVector = frontier->getNextArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT m = 0;
    intT outEdgesCount = 0;
    bool *nextB = next->b;

    intT chunkSize = numVertices / frontier->numOfNodes;
    intT rollingOffset = 0;//subworker.tid * chunkSize;

    intT *counterPtr = &(next->sparseCounter);

    intT oldStartPos = 0;

    intT startPos = __sync_fetch_and_add(counterPtr, DYNAMIC_CHUNK_SIZE);
    intT endPos = (startPos + DYNAMIC_CHUNK_SIZE > numVertices) ? (numVertices) : (startPos + DYNAMIC_CHUNK_SIZE);
    while (startPos < numVertices) {
        m = 0;
        intT actualSize = endPos - startPos;
        intT actualStart = startPos + rollingOffset;
        while (actualStart >= nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            if (currNodeNum >= frontier->numOfNodes) {
                currNodeNum = 0;
                currOffset = 0;
                nextSwitchPoint = frontier->getSize(0) + numVertices;
                //printf("switch back\n");
            } else {
                currBitVector = frontier->getNextArr(currNodeNum);
            }
        }
        for (long i = startPos; i < endPos; i++) {
            intT idx = (i + rollingOffset) % numVertices;
            if (i + rollingOffset == nextSwitchPoint) {
                currOffset += frontier->getSize(currNodeNum);
                nextSwitchPoint += frontier->getSize(currNodeNum + 1);
                currNodeNum++;
                if (currNodeNum >= frontier->numOfNodes) {
                    currNodeNum = 0;
                    currOffset = 0;
                    nextSwitchPoint = frontier->getSize(0) + numVertices;
                    //printf("switch back\n");
                } else {
                    currBitVector = frontier->getNextArr(currNodeNum);
                }
            }
            m += G[idx].getFakeInDegree();
            if (f.cond(idx)) {
                intT d = G[idx].getFakeInDegree();
                //printf("in deg of %d: %d\n", idx, d);
                for (intT j = 0; j < d; j++) {
                    uintT ngh = G[idx].getInNeighbor(j);
                    if (localBitVec[ngh - localOffset] && f.updateAtomic(ngh, idx)) {
                        currBitVector[idx - currOffset] = true;
                    }
                    if (!f.cond(idx)) {
                        break;
                    }
                }
            }
        }

        if (subworker.tid == 7) {
            //printf("%d: %d to %d deg: %d\n", subworker.subTid, startPos, endPos, m);
        }

        oldStartPos = startPos;
        startPos = __sync_fetch_and_add(counterPtr, DYNAMIC_CHUNK_SIZE);
        endPos = (startPos + DYNAMIC_CHUNK_SIZE > numVertices) ? (numVertices) : (startPos + DYNAMIC_CHUNK_SIZE);
    }
    return NULL;
}

template<class F, class vertex>
bool *edgeMapDenseBP(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false, int start = 0,
                     int end = 0) {
    intT numVertices = GA.n;
    vertex *G = GA.V;

    int currNodeNum = 0;
    bool *currBitVector = frontier->getArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;
    intT m = 0;
    intT outEdgesCount = 0;
    bool *nextB = next->b;
    int startPos = 0;
    int endPos = numVertices;
    if (part) {
        startPos = start;
        endPos = end;
        currNodeNum = frontier->getNodeNumOfIndex(startPos);
        currBitVector = frontier->getArr(currNodeNum);
        nextSwitchPoint = frontier->getOffset(currNodeNum + 1);
        currOffset = frontier->getOffset(currNodeNum);
    }
    for (long i = startPos; i < endPos; i++) {
        if (i == nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            currBitVector = frontier->getArr(currNodeNum);
        }
        m += G[i].getFakeDegree();
        if (currBitVector[i - currOffset]) {
            intT d = G[i].getFakeDegree();
            for (intT j = 0; j < d; j++) {
                uintT ngh = G[i].getOutNeighbor(j);
                if (/*next->inRange(ngh) &&*/ f.cond(ngh) && f.updateAtomic(i, ngh, j)) {
                    next->setBit(ngh, true);
                }
            }
        }
    }
    return NULL;
}

template<class F, class vertex>
bool *edgeMapDenseForwardGlobalWrite(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *nexts[],
                                     Subworker_Partitioner &subworker) {
    intT numVertices = GA.n;
    vertex *G = GA.V;

    bool *currBitVector = frontier->getArr(subworker.tid);
    int currOffset = frontier->getOffset(subworker.tid);
    int counter = 0;

    intT m = 0;
    intT outEdgesCount = 0;

    int startPos = subworker.dense_start;
    int endPos = subworker.dense_end;

    //printf("%d %d: start-end: %d %d\n", subworker.tid, subworker.subTid, startPos, endPos);

    int currNodeNum = frontier->getNodeNumOfIndex(startPos);
    bool *nextBitVector = nexts[currNodeNum]->b;
    intT nextSwitchPoint = frontier->getOffset(currNodeNum + 1);
    int offset = frontier->getOffset(currNodeNum);

    for (long i = startPos; i < endPos; i++) {
        if (i == nextSwitchPoint) {
            offset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            nextBitVector = nexts[currNodeNum]->b;
        }
        m += G[i].getFakeDegree();
        if (f.cond(i)) {
            intT d = G[i].getFakeDegree();
            for (intT j = 0; j < d; j++) {
                uintT ngh = G[i].getInNeighbor(j);
                if (currBitVector[ngh - currOffset] && f.updateAtomic(ngh, i)) {
                    nextBitVector[i - offset] = true;
                }
            }
        }
    }
    return NULL;
}

AsyncChunk *newChunk(int blockSize) {
    AsyncChunk *myChunk = (AsyncChunk *) malloc(sizeof(AsyncChunk));
    myChunk->s = (intT *) malloc(sizeof(intT) * blockSize);
    myChunk->m = 0;
    myChunk->accessCounter = 0;
    return myChunk;
}

template<class F, class vertex>
void edgeMapSparseAsync(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next,
                        Subworker_Partitioner &subworker = dummyPartitioner) {
    const int BLOCK_SIZE = 64;

    vertex *V = GA.V;

    int tid = subworker.tid;
    volatile intT *queueHead = &(frontier->frontiers[tid]->head);
    volatile intT *queueTail = &(frontier->readerTail);
    volatile intT *insertTail = &(frontier->insertTail);
    volatile intT *endSignal = &(frontier->asyncEndSignal);
    volatile intT *localSignal = &(frontier->frontiers[tid]->emptySignal);
    volatile intT *signals[frontier->numOfNodes];
    for (int i = 0; i < frontier->numOfNodes; i++) {
        signals[i] = &(frontier->frontiers[tid]->emptySignal);
    }
    *queueHead = 0;
    *insertTail = *queueTail;
    *localSignal = 0;
    *endSignal = 0;
    pthread_barrier_wait(subworker.local_barr);
    if (subworker.isSubMaster()) {
        printf("passed barrier\n");
    }
    int accumSize = 0;
    AsyncChunk *myChunk = newChunk(BLOCK_SIZE);
    bool shouldFinish = false;
    while (!shouldFinish) {
        //first fetch a block
        volatile intT currHead = 0;
        volatile intT currTail = 0;
        intT endPos = 0;
        AsyncChunk *currChunk = NULL;
        /*
        if (*endSignal == 1) {
            printf("spin: %d %d %d %d %d\n", subworker.tid, subworker.subTid, currHead, currTail, endPos);
        }
        */
        do {
            currHead = *queueHead;
            currTail = *queueTail;
            endPos = MIN(currHead + 1, currTail);
        } while (!__sync_bool_compare_and_swap((intT *) queueHead, currHead, endPos));

        int reallyGotOne = endPos - currHead;
        //printf("get: %d, %d\n", currHead, endPos);
        if (reallyGotOne > 0) {
            *localSignal = 0;

            if (*endSignal == 1) {
                //printf("possible? %d %d %d\n", currHead, endPos, currTail);
            }

            //process chunk
            currChunk = frontier->asyncQueue[currHead % GA.n];
            if (currChunk == NULL) {
                printf("oops: %p %p %d %d %d\n", currChunk, frontier->asyncQueue[currHead % GA.n], currHead, currTail,
                       endPos);
            }
            //printf("chunk pointer: %p\n", currChunk);
            int chunkSize = currChunk->m;
            for (intT i = 0; i < chunkSize; i++) {
                accumSize++;
                intT idx = currChunk->s[i];
                intT d = V[idx].getOutDegree();
                for (intT j = 0; j < d; j++) {
                    intT ngh = V[idx].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                        //add ngh into chunk
                        myChunk->s[myChunk->m] = ngh;
                        myChunk->m += 1;
                        if (myChunk->m >= BLOCK_SIZE) {
                            //if full, send it
                            intT insertPos = __sync_fetch_and_add(insertTail, 1);
                            /*
                            if (*endSignal == 1)
                            printf("before insert %d %d %d\n", *queueTail, insertPos, *insertTail);
                            */
                            frontier->asyncQueue[insertPos % GA.n] = myChunk;
                            //__asm__ __volatile__ ("mfence\n":::);
                            while (!__sync_bool_compare_and_swap((intT *) queueTail, insertPos, insertPos + 1)) {
                                if (*queueTail > insertPos) {
                                    break;
                                    printf("pending on insert %d %d %d\n", *queueTail, insertPos, *insertTail);
                                }

                            }
                            //printf("insert over: %d %d\n", insertPos, *queueTail);
                            myChunk = newChunk(BLOCK_SIZE);
                        }
                    }
                }
            }
            /*
            int oldCounter = __sync_fetch_and_add(&(currChunk->accessCounter), 1);
            oldCounter++;

            if (oldCounter >= frontier->numOfNodes) {
            frontier->asyncQueue[currHead % GA.n] = NULL;
            free(currChunk);
            }
            */
        } else {
            if (myChunk->m > 0) {
                //send it
                //printf("flush chunk with size: %d\n", myChunk->m);
                intT insertPos = __sync_fetch_and_add(insertTail, 1);
                frontier->asyncQueue[insertPos % GA.n] = myChunk;
                //__asm__ __volatile__ ("mfence\n":::);
                while (!__sync_bool_compare_and_swap((intT *) queueTail, insertPos, insertPos + 1)) {
                    if (*queueTail > insertPos) {
                        break;
                        printf("pending on insert %d %d %d\n", *queueTail, insertPos, *insertTail);
                    }
                }
                //printf("insert over: %d %d\n", insertPos, *queueTail);
                //__sync_fetch_and_add(queueTail, 1);
                myChunk = newChunk(BLOCK_SIZE);
                continue;
            }
            //end game part
            //pthread_barrier_wait(subworker.local_barr);
            //printf("%d %d here %d %d\n", subworker.tid, subworker.subTid, currHead, currTail);
            if (*endSignal == 1) {
                //printf("%d %d and here %d %d\n", subworker.tid, subworker.subTid, currHead, currTail);
                shouldFinish = true;
            }
            if (subworker.isMaster()) {
                //marker algorithm
                *localSignal = 1;
                int i = 0;
                int marker = 1;
                while (*localSignal == 1) {
                    //printf("checking\n");
                    i = (i + 1) % frontier->numOfNodes;
                    if (*(signals[i]) == 1) {
                        marker++;
                    } else {
                        marker = 0;
                    }
                    *(signals[i]) = 1;
                    if (marker > frontier->numOfNodes) {
                        *endSignal = 1;
                        printf("master out\n");
                        shouldFinish = true;
                        break;
                    }
                }
            }
            //pthread_barrier_wait(subworker.local_barr);
        }
    }
    //printf("end loop of %d %d: %d\n", subworker.tid, subworker.subTid, accumSize);
}

template<class F, class vertex>
void edgeMapSparseAsyncPipe(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next,
                            Subworker_Partitioner &subworker = dummyPartitioner) {
    const int BLOCK_SIZE = 64;
    vertex *V = GA.V;
    int tid = subworker.tid;

    frontier->frontiers[tid]->tmp = (intT *) malloc(sizeof(intT) * frontier->getSize(tid));
    pthread_barrier_wait(subworker.local_barr);

    AsyncChunk **localQueue = frontier->frontiers[tid]->localQueue;
    AsyncChunk **nextQueue = frontier->frontiers[(tid + 1) % frontier->numOfNodes]->localQueue;

    volatile intT *localHead = &(frontier->frontiers[tid]->head);
    volatile intT *localTail = &(frontier->frontiers[tid]->tail);
    volatile intT *nextTail = &(frontier->frontiers[(tid + 1) % frontier->numOfNodes]->tail);
    volatile intT *insertTail = &(frontier->frontiers[(tid + 1) % frontier->numOfNodes]->insertTail);
    volatile intT *localSignal = &(frontier->frontiers[tid]->emptySignal);
    volatile intT *endSignal = &(frontier->asyncEndSignal);
    volatile intT *endGameOnFly = &(frontier->readerTail);
    volatile intT *signals[frontier->numOfNodes];
    for (int i = 0; i < frontier->numOfNodes; i++) {
        signals[i] = &(frontier->frontiers[i]->emptySignal);
    }
    *localHead = 0;
    *insertTail = *nextTail;
    *localSignal = 0;
    *endSignal = 0;
    *endGameOnFly = 0;

    int offset = frontier->getOffset(tid);
    int *bitVec = frontier->frontiers[tid]->tmp;
    int size = frontier->getSize(tid);
    int subSize = size / subworker.numOfSub;
    int startPos = subSize * subworker.subTid;
    int endPos = subSize * (subworker.subTid + 1);
    if (subworker.subTid == subworker.numOfSub - 1) {
        endPos = size;
    }

    for (int i = startPos; i < endPos; i++) {
        bitVec[i] = 0;
    }

    pthread_barrier_wait(subworker.local_barr);

    int accumSize = 0;
    AsyncChunk *myChunk = newChunk(BLOCK_SIZE);
    bool shouldFinish = false;
    *localSignal = 0;
    while (!shouldFinish) {
        volatile intT currHead = 0;
        volatile intT currTail = 0;
        intT endPos = 0;
        AsyncChunk *currChunk = NULL;

        do {
            currHead = *localHead;
            currTail = *localTail;
            endPos = MIN(currHead + 1, currTail);
        } while (!__sync_bool_compare_and_swap((intT *) localHead, currHead, endPos));

        int reallyGotOne = endPos - currHead;
        if (reallyGotOne > 0) {
            *localSignal = 0;

            currChunk = localQueue[currHead % GA.n];
            int chunkSize = currChunk->m;
            if (chunkSize > 0) {
                for (intT i = 0; i < chunkSize; i++) {
                    accumSize++;
                    intT idx = currChunk->s[i];
                    intT d = V[idx].getOutDegree();
                    for (intT j = 0; j < d; j++) {
                        intT ngh = V[idx].getOutNeighbor(j);
                        if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                            //add ngh into chunk
                            int counter = __sync_fetch_and_add(&(bitVec[ngh - offset]), 1);
                            if (counter == 0) {
                                myChunk->s[myChunk->m] = ngh;
                                myChunk->m += 1;
                                if (myChunk->m >= BLOCK_SIZE) {
                                    //if full, send it
                                    //printf("%d send to %d\n", tid, (tid + 1) % frontier->numOfNodes);
                                    intT insertPos = __sync_fetch_and_add(insertTail, 1);
                                    nextQueue[insertPos % GA.n] = myChunk;
                                    while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                                        if (*nextTail > insertPos) {
                                            break;
                                            printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                                        }
                                    }
                                    myChunk = newChunk(BLOCK_SIZE);
                                }
                            }
                        }
                    }
                }
            } else {
                if (chunkSize < 0) {
                    // end game message.
                    if (currChunk->accessCounter >= 2 * frontier->numOfNodes && subworker.isMaster()) {
                        *endSignal = 1;
                        printf("master out\n");
                        shouldFinish = true;
                        continue;
                    }

                    if (*localHead >= *localTail && myChunk->m <= 0) {
                        //forward this to the next node
                        intT insertPos = __sync_fetch_and_add(insertTail, 1);
                        nextQueue[insertPos % GA.n] = currChunk;
                        currChunk->accessCounter++;
                        while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                            if (*nextTail > insertPos) {
                                break;
                                printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                            }
                        }
                    } else {
                        free(currChunk);
                        *endGameOnFly = 0;
                    }
                    continue;
                }
            }

            int oldCounter = __sync_fetch_and_add(&(currChunk->accessCounter), 1);
            oldCounter++;

            if (oldCounter >= frontier->numOfNodes) {
                localQueue[currHead % GA.n] = NULL;
                for (int i = 0; i < chunkSize; i++) {
                    intT idx = currChunk->s[i];
                    if (bitVec[idx - offset] <= 1) {
                        bitVec[idx - offset] = 0;
                    } else {
                        bitVec[idx - offset] = 1;
                        myChunk->s[myChunk->m] = idx;
                        myChunk->m += 1;
                        if (myChunk->m >= BLOCK_SIZE) {
                            //if full, send it
                            intT insertPos = __sync_fetch_and_add(insertTail, 1);
                            nextQueue[insertPos % GA.n] = myChunk;
                            while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                                if (*nextTail > insertPos) {
                                    break;
                                    printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                                }
                            }
                            myChunk = newChunk(BLOCK_SIZE);
                        }
                    }
                }
                free(currChunk);
            } else {
                //forward the chunk to next
                //printf("forward chunk from %d to %d\n", tid, (tid + 1) % frontier->numOfNodes);
                intT insertPos = __sync_fetch_and_add(insertTail, 1);
                nextQueue[insertPos % GA.n] = currChunk;
                while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                    if (*nextTail > insertPos) {
                        break;
                        printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                    }
                }
            }
        } else {
            //first send out what is left
            if (myChunk->m > 0) {
                //if full, send it
                intT insertPos = __sync_fetch_and_add(insertTail, 1);
                nextQueue[insertPos % GA.n] = myChunk;
                while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                    if (*nextTail > insertPos) {
                        break;
                        printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                    }

                }
                //printf("%d send leftover %p to %d at %d\n", tid, myChunk, (tid + 1) % frontier->numOfNodes, insertPos);
                myChunk = newChunk(BLOCK_SIZE);
            }

            if (*endSignal == 1) {
                shouldFinish = true;
            }

            //end game protocol
            if (subworker.isMaster()) {
                /*
                *localSignal = 1;
                int i = 0;
                int marker = 1;
                while (*localSignal == 1) {
                    //printf("checking\n");
                    i = (i + 1) % frontier->numOfNodes;
                    if (*(signals[i]) == 1) {
                    marker++;
                    } else {
                    marker = 0;
                    }
                    *(signals[i]) = 1;
                    if (marker > 3 * frontier->numOfNodes) {
                    *endSignal = 1;
                    //printf("master out\n");
                    shouldFinish = true;
                    break;
                    }
                }
                */

                // create end game chunk and send it.
                if (*endGameOnFly == 0) {
                    printf("sent end game\n");
                    AsyncChunk *endGameChunk = (AsyncChunk *) malloc(sizeof(AsyncChunk));
                    endGameChunk->accessCounter = 1;
                    endGameChunk->m = -1; //magic number for end game chunk.
                    intT insertPos = __sync_fetch_and_add(insertTail, 1);
                    nextQueue[insertPos % GA.n] = endGameChunk;
                    *endGameOnFly = 1;
                    while (!__sync_bool_compare_and_swap((intT *) nextTail, insertPos, insertPos + 1)) {
                        if (*nextTail > insertPos) {
                            break;
                            printf("pending on insert %d %d %d\n", *nextTail, insertPos, *insertTail);
                        }
                    }
                }
            } else {
                bool firstRound = true;
                while (*localHead >= *localTail && *endSignal != 1) {
                    if (firstRound) {
                        *localSignal = 1;
                        firstRound = false;
                    }
                }
                *localSignal = 0;
            }
        }
    }
}

void switchFrontier(int nodeNum, vertices *V, LocalFrontier *&next);

template<class F, class vertex>
void edgeMapSparseV5(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next,
                     Subworker_Partitioner &subworker = dummyPartitioner) {
    vertex *V = GA.V;
    intT currM = frontier->numNonzeros();
    int startPos = 0;//subworker.getStartPos(currM);
    int endPos = currM;//subworker.getEndPos(currM);
    if (!subworker.isSubMaster()) {
        startPos = 1;
        endPos = 0;
    }

    int counter = 0;

    while (startPos < endPos) { //only master get in
        counter++;
        next->m = 0;
        next->outEdgesCount = 0;
        int bufferLen = frontier->getEdgeStat();
        if (subworker.isSubMaster()) {
            next->s = (intT *) malloc(sizeof(intT) * bufferLen);
        }
        intT nextEdgesCount = 0;

        //pthread_barrier_wait(subworker.local_barr);
        intT *nextFrontier = next->s;
        int tmp = 0;
        int currNodeNum = frontier->getNodeNumOfSparseIndex(startPos);
        int offset = 0;
        for (int i = 0; i < currNodeNum; i++) {
            offset += frontier->getSparseSize(i);
        }
        intT *currActiveList = frontier->getSparseArr(currNodeNum);
        int lengthOfCurr = frontier->getSparseSize(currNodeNum) - (startPos - offset);
        for (int i = startPos; i < endPos; i++) {
            if (lengthOfCurr <= 0) {
                while (currNodeNum + 1 < frontier->numOfNodes && lengthOfCurr <= 0) {
                    offset += frontier->getSparseSize(currNodeNum);
                    currNodeNum++;
                    lengthOfCurr = frontier->getSparseSize(currNodeNum);
                }
                if (currNodeNum >= frontier->numOfNodes || lengthOfCurr <= 0) {
                    printf("oops\n");
                }
                currActiveList = frontier->getSparseArr(currNodeNum);
            }
            intT idx = currActiveList[i - offset];
            intT d = V[idx].getFakeDegree();
            for (intT j = 0; j < d; j++) {
                uintT ngh = V[idx].getOutNeighbor(j);
                if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                    next->s[tmp] = ngh;
                    tmp++;
                    nextEdgesCount += V[ngh].getOutDegree();
                }
            }
            lengthOfCurr--;
        }
        next->m = tmp;
        next->outEdgesCount = nextEdgesCount;
        pthread_barrier_wait(subworker.leader_barr);
        switchFrontier(subworker.tid, frontier, next);
        pthread_barrier_wait(subworker.leader_barr);
        frontier->calculateNumOfNonZero(subworker.tid);
        pthread_barrier_wait(subworker.leader_barr);
        currM = frontier->numNonzeros();
        startPos = 0;
        endPos = currM;
        if (frontier->isEmpty()) {
            if (subworker.tid == 0) {
                printf("Sparse ok: %d\n", counter);
            }
            break;
        }
    }
    //pthread_barrier_wait(subworker.local_barr);
}

template<class F, class vertex>
void edgeMapSparseV4(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool firstTime = false,
                     Subworker_Partitioner &subworker = dummyPartitioner) {
    // in V4, all thread has its own chunk.
    vertex *V = GA.V;

    intT nextM = 0;
    intT *nextChunk = (intT *) malloc(sizeof(intT) * next->n);
    intT nextEdgesCount = 0;
    if (firstTime) {
        intT currM = frontier->numNonzeros();
        int startPos = subworker.getStartPos(currM);
        int endPos = subworker.getEndPos(currM);

        next->outEdgesCount = 0;
        int bufferLen = frontier->getEdgeStat();

        //pthread_barrier_wait(subworker.local_barr);
        subworker.localWait();

        if (startPos < endPos) {
            int currNodeNum = frontier->getNodeNumOfSparseIndex(startPos);
            int offset = 0;
            for (int i = 0; i < currNodeNum; i++) {
                offset += frontier->getSparseSize(i);
            }
            intT *currActiveList = frontier->getSparseArr(currNodeNum);
            int lengthOfCurr = frontier->getSparseSize(currNodeNum) - (startPos - offset);
            for (int i = startPos; i < endPos; i++) {
                if (lengthOfCurr <= 0) {
                    while (currNodeNum + 1 < frontier->numOfNodes && lengthOfCurr <= 0) {
                        offset += frontier->getSparseSize(currNodeNum);
                        currNodeNum++;
                        lengthOfCurr = frontier->getSparseSize(currNodeNum);
                    }
                    if (currNodeNum >= frontier->numOfNodes || lengthOfCurr <= 0) {
                        printf("oops\n");
                    }
                    currActiveList = frontier->getSparseArr(currNodeNum);
                }
                intT idx = currActiveList[i - offset];
                intT d = V[idx].getFakeDegree();
                for (intT j = 0; j < d; j++) {
                    uintT ngh = V[idx].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                        nextChunk[nextM] = ngh;
                        nextM++;
                        nextEdgesCount += V[ngh].getOutDegree();
                    }
                }
                lengthOfCurr--;
            }
        }
        next->chunkSizes[subworker.subTid] = nextM;
        next->sparseChunks[subworker.subTid] = nextChunk;
        __sync_fetch_and_add(&(next->m), nextM);
        __sync_fetch_and_add(&(next->outEdgesCount), nextEdgesCount);
        //pthread_barrier_wait(subworker.global_barr);
        subworker.globalWait();
    } else {
        intT numOfChunks = subworker.numOfSub * frontier->numOfNodes;
        next->sparseCounter = 0;
        next->m = 0;
        //pthread_barrier_wait(subworker.local_barr);
        subworker.localWait();

        intT fetchedChunk = __sync_fetch_and_add(&(next->sparseCounter), 1);
        while (fetchedChunk < numOfChunks) {
            int nodeIdx = fetchedChunk / subworker.numOfSub;
            int subIdx = fetchedChunk % subworker.numOfSub;
            intT *chunk = frontier->frontiers[nodeIdx]->sparseChunks[subIdx];
            //printf("%d %d curr chunk is: %p\n", nodeIdx, subIdx, chunk);
            intT chunkSize = frontier->frontiers[nodeIdx]->chunkSizes[subIdx];
            for (int i = 0; i < chunkSize; i++) {
                intT idx = chunk[i];
                intT d = V[idx].getFakeDegree();
                for (intT j = 0; j < d; j++) {
                    uintT ngh = V[idx].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                        nextChunk[nextM] = ngh;
                        nextM++;
                        nextEdgesCount += V[ngh].getOutDegree();
                    }
                }
                __builtin_prefetch(&(chunk[i + 2]), 0, 3);
            }
            fetchedChunk = __sync_fetch_and_add(&(next->sparseCounter), 1);
        }
        next->chunkSizes[subworker.subTid] = nextM;
        next->sparseChunks[subworker.subTid] = nextChunk;
        __sync_fetch_and_add(&(next->m), nextM);
        __sync_fetch_and_add(&(next->outEdgesCount), nextEdgesCount);
        //pthread_barrier_wait(subworker.global_barr);
        subworker.globalWait();
    }
}

template<class F, class vertex>
void edgeMapSparseV3(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false,
                     Subworker_Partitioner &subworker = dummyPartitioner) {
    vertex *V = GA.V;
    if (part) {
        intT currM = frontier->numNonzeros();
        int startPos = subworker.getStartPos(currM);
        int endPos = subworker.getEndPos(currM);

        intT *mPtr = &(next->m);
        *mPtr = 0;
        next->outEdgesCount = 0;
        int bufferLen = frontier->getEdgeStat();
        if (subworker.isSubMaster())
            next->s = (intT *) malloc(sizeof(intT) * bufferLen);
        intT nextEdgesCount = 0;

        //pthread_barrier_wait(subworker.local_barr);
        subworker.localWait();
        intT *nextFrontier = next->s;

        if (startPos < endPos) {
            //printf("have ele: %d to %d %d, %p\n", startPos, endPos, subworker.tid, next);
            int currNodeNum = frontier->getNodeNumOfSparseIndex(startPos);
            int offset = 0;
            for (int i = 0; i < currNodeNum; i++) {
                offset += frontier->getSparseSize(i);
            }
            intT *currActiveList = frontier->getSparseArr(currNodeNum);
            int lengthOfCurr = frontier->getSparseSize(currNodeNum) - (startPos - offset);
            //printf("nodeNum of %d %d: %d from %d to %d\n", subworker.tid, subworker.subTid, currNodeNum, startPos, endPos);
            for (int i = startPos; i < endPos; i++) {
                if (lengthOfCurr <= 0) {
                    while (currNodeNum + 1 < frontier->numOfNodes && lengthOfCurr <= 0) {
                        offset += frontier->getSparseSize(currNodeNum);
                        currNodeNum++;
                        lengthOfCurr = frontier->getSparseSize(currNodeNum);
                        //printf("lengthOfCurr: %d\n", lengthOfCurr);
                    }
                    if (currNodeNum >= frontier->numOfNodes || lengthOfCurr <= 0) {
                        printf("oops\n");
                    }
                    currActiveList = frontier->getSparseArr(currNodeNum);
                    //printf("currList of %d %d: %p\n", subworker.tid, subworker.subTid, currActiveList);
                }
                //printf("s: %p %d\n", currActiveList, i-offset);
                intT idx = currActiveList[i - offset];
                //printf("vertex on %d %d: %d\n", subworker.tid, subworker.subTid, idx);
                intT d = V[idx].getFakeDegree();
                //printf("degree: %d\n", d);
                for (intT j = 0; j < d; j++) {
                    uintT ngh = V[idx].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                        //add to active list
                        //printf("out edge # %d: %d -> %d of %d %d\n", nextM, idx, ngh, subworker.tid, subworker.subTid);
                        /*
                        if (nextM >= bufferLen) {
                            printf("oops: %d %d\n", subworker.tid, subworker.subTid);
                        }
                        */
                        //printf("I am here\n");
                        int tmp = __sync_fetch_and_add(mPtr, 1);
                        if (tmp >= bufferLen)
                            printf("oops\n");
                        nextFrontier[tmp] = ngh;
                        nextEdgesCount += V[ngh].getOutDegree();
                    }
                }
                lengthOfCurr--;
                //printf("nextM: %d %d\n", nextM, nextEdgesCount);
            }
        }
        __sync_fetch_and_add(&(next->outEdgesCount), nextEdgesCount);
        //pthread_barrier_wait(subworker.local_barr);
        subworker.localWait();
    }
}

template<class F, class vertex>
void edgeMapSparseV2(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false,
                     Subworker_Partitioner &subworker = dummyPartitioner) {
    vertex *V = GA.V;
    if (part) {
        intT currM = frontier->numNonzeros();
        int startPos = subworker.getStartPos(currM);
        int endPos = subworker.getEndPos(currM);

        intT nextM = 0;
        intT nextEdgesCount = 0;
        intT *nextFrontier = NULL;
        pthread_barrier_wait(subworker.local_barr);

        if (startPos < endPos) {
            //printf("have ele: %d to %d %d, %p\n", startPos, endPos, subworker.tid, next);
            int bufferLen = frontier->getEdgeStat();
            nextFrontier = (intT *) malloc(sizeof(intT) * bufferLen);

            int currNodeNum = frontier->getNodeNumOfSparseIndex(startPos);
            int offset = 0;
            for (int i = 0; i < currNodeNum; i++) {
                offset += frontier->getSparseSize(i);
            }
            intT *currActiveList = frontier->getSparseArr(currNodeNum);
            int lengthOfCurr = frontier->getSparseSize(currNodeNum) - (startPos - offset);
            //printf("nodeNum of %d %d: %d from %d to %d\n", subworker.tid, subworker.subTid, currNodeNum, startPos, endPos);
            for (int i = startPos; i < endPos; i++) {
                if (lengthOfCurr <= 0) {
                    while (currNodeNum + 1 < frontier->numOfNodes && lengthOfCurr <= 0) {
                        offset += frontier->getSparseSize(currNodeNum);
                        currNodeNum++;
                        lengthOfCurr = frontier->getSparseSize(currNodeNum);
                        //printf("lengthOfCurr: %d\n", lengthOfCurr);
                    }
                    if (currNodeNum >= frontier->numOfNodes || lengthOfCurr <= 0) {
                        printf("oops\n");
                    }
                    currActiveList = frontier->getSparseArr(currNodeNum);
                    //printf("currList of %d %d: %p\n", subworker.tid, subworker.subTid, currActiveList);
                }
                //printf("s: %p %d\n", currActiveList, i-offset);
                intT idx = currActiveList[i - offset];
                //printf("vertex on %d %d: %d\n", subworker.tid, subworker.subTid, idx);
                intT d = V[idx].getFakeDegree();
                for (intT j = 0; j < d; j++) {
                    uintT ngh = V[idx].getOutNeighbor(j);
                    if (f.cond(ngh) && f.updateAtomic(idx, ngh)) {
                        //add to active list
                        //printf("out edge # %d: %d -> %d of %d %d\n", nextM, idx, ngh, subworker.tid, subworker.subTid);
                        if (nextM >= bufferLen) {
                            printf("oops: %d %d\n", subworker.tid, subworker.subTid);
                        }
                        nextFrontier[nextM] = ngh;
                        nextM++;
                        nextEdgesCount += V[ngh].getOutDegree();
                    }
                }
                lengthOfCurr--;
                //printf("nextM: %d %d\n", nextM, nextEdgesCount);
            }
        } else {
            nextM = 0;
            nextEdgesCount = 0;
        }
        struct timeval start1, end1, start2, end2, start3, end3, startT, endT;
        struct timezone tz = {0, 0};
        //gettimeofday(&startT, &tz);
        writeAdd(&(next->m), nextM);
        writeAdd(&(next->outEdgesCount), nextEdgesCount);
        if (subworker.isSubMaster()) {
            intT *offsets = (intT *) malloc(sizeof(intT) * subworker.numOfSub);
            next->tmp = offsets;
        }
        //gettimeofday(&start1, &tz);
        pthread_barrier_wait(subworker.local_barr);
        //gettimeofday(&end1, &tz);
        next->tmp[subworker.subTid] = nextM;
        if (subworker.isSubMaster()) {
            //printf("next of %d: %d %d\n", subworker.tid, next->m, nextM);
            if (next->m > 0) {
                intT *sparseArr = (intT *) malloc(sizeof(intT) * next->m);
                //printf("next sparse: %p\n", sparseArr);
                next->s = sparseArr;
            }
            next->isDense = false;
        }
        //gettimeofday(&start2, &tz);
        pthread_barrier_wait(subworker.local_barr);
        //gettimeofday(&end2, &tz);
        if (nextM > 0) {
            intT fillOffset = 0;
            for (int i = 0; i < subworker.subTid; i++) {
                fillOffset += next->tmp[i];
            }
            for (intT i = 0; i < nextM; i++) {
                next->s[i + fillOffset] = nextFrontier[i];
                if (i + fillOffset >= next->m) {
                    printf("oops\n");
                }
                //printf("filled to %d of %d: %d\n", i + fillOffset, subworker.tid, nextFrontier[i]);
            }
            if (startPos < endPos) {
                free(nextFrontier);
            }
        }
        //gettimeofday(&start3, &tz);
        pthread_barrier_wait(subworker.local_barr);
        //gettimeofday(&end3, &tz);
        //gettimeofday(&endT, &tz);
        /*
        if (subworker.isMaster()) {
            double timeStart = ((double)start1.tv_sec) + ((double)start1.tv_usec) / 1000000.0;
            double timeEnd = ((double)end1.tv_sec) + ((double)end1.tv_usec) / 1000000.0;
            double mapTime = timeEnd - timeStart;
            printf("#1:after map time: %lf\n", mapTime);

            timeStart = ((double)start2.tv_sec) + ((double)start2.tv_usec) / 1000000.0;
            timeEnd = ((double)end2.tv_sec) + ((double)end2.tv_usec) / 1000000.0;
            mapTime = timeEnd - timeStart;
            printf("#2:after map time: %lf\n", mapTime);

            timeStart = ((double)start3.tv_sec) + ((double)start3.tv_usec) / 1000000.0;
            timeEnd = ((double)end3.tv_sec) + ((double)end3.tv_usec) / 1000000.0;
            mapTime = timeEnd - timeStart;
            printf("#3:after map time: %lf\n", mapTime);

            timeStart = ((double)startT.tv_sec) + ((double)startT.tv_usec) / 1000000.0;
            timeEnd = ((double)endT.tv_sec) + ((double)endT.tv_usec) / 1000000.0;
            mapTime = timeEnd - timeStart;
            printf("Total:after map time: %lf\n", mapTime);
        }
        */
    }
}

template<class F, class Vert_F, class vertex>
void edgeMapSparse(graph<vertex> GA, vertices *frontier, F f, Vert_F vf, LocalFrontier *next,
                   Subworker_Partitioner &subworker = dummyPartitioner) {
    vertex *V = GA.V;
    intT sparseIter = 0;
    if (!subworker.isMaster()) {
        return;
    } else {
        frontier->toSparse();
    }
    //for first time put all frontiers together
    intT totM = frontier->numNonzeros();
    intT *sparseQueue = (intT *) malloc(sizeof(intT) * totM);
    intT frontierOffsets[frontier->numOfNodes];
    frontierOffsets[0] = 0;
    for (int i = 0; i < frontier->numOfNodes - 1; i++) {
        frontierOffsets[i + 1] = frontierOffsets[i] + frontier->frontiers[i]->m;
    }
    for (int i = 0; i < frontier->numOfNodes; i++) {
        intT o = frontierOffsets[i];
        intT *s = frontier->frontiers[i]->s;
        {
            parallel_for (intT j = 0; j < frontier->frontiers[i]->m; j++) {
                sparseQueue[o + j] = s[j];
            }
        }
    }
    while (true) {
        sparseIter++;
        intT *degrees = (intT *) malloc(sizeof(intT) * totM);
        {
            parallel_for (intT i = 0; i < totM; i++) {
                degrees[i] = V[sparseQueue[i]].getOutDegree();
            }
        }
        printf("sparse: %d\n", totM);
        //printf("in loop : %d %p\n", frontier->numOfNodes, frontierOffsets);
        uintT *offsets = (uintT *) degrees;
        uintT outEdgeCount = sequence::plusScan(offsets, (uintT *) degrees, (uintT) totM);
        intT newM = 0;
        intT *outEdges = (intT *) malloc(sizeof(intT) * outEdgeCount);
        {
            parallel_for (intT i = 0; i < totM; i++) {
                intT v = sparseQueue[i];
                intT o = offsets[i];
                vertex vert = V[sparseQueue[i]];
                intT d = vert.getOutDegree();
                if (d < 1000) {
                    for (intT j = 0; j < d; j++) {
                        intT ngh = vert.getOutNeighbor(j);
                        if (f.cond(ngh) && f.updateAtomic(v, ngh))
                            outEdges[o + j] = ngh;
                        else
                            outEdges[o + j] = -1;
                    }
                } else {
                    {
                        parallel_for (intT j = 0; j < d; j++) {
                            intT ngh = vert.getOutNeighbor(j);
                            if (f.cond(ngh) && f.updateAtomic(v, ngh))
                                outEdges[o + j] = ngh;
                            else
                                outEdges[o + j] = -1;
                        }
                    }
                }
            }
        }
        //printf("sparse mode %d: %d\n", next->startID, degreeCount);

        intT *nextIndices = (intT *) malloc(sizeof(intT) * outEdgeCount);
        newM = sequence::filter(outEdges, nextIndices, outEdgeCount, nonNegF());
        if (newM <= 0) {
            printf("sparseIter: %d %d\n", sparseIter, newM);
            break;
        } else {
            if (sparseQueue != NULL) {
                free(sparseQueue);
            }
            sparseQueue = nextIndices;
            totM = newM;
            {
                parallel_for (intT i = 0; i < totM; i++) {
                    vf(sparseQueue[i]);
                }
            }
        }
    }
}

static int edgesTraversed = 0;

void switchFrontier(int nodeNum, vertices *V, LocalFrontier *&next) {
    LocalFrontier *current = V->getFrontier(nodeNum);
    intT size = V->getSize(nodeNum);
    /*
    for (intT i = 0; i < size; i++) {
	(current->b)[i] = false;
    }
    */
    LocalFrontier *newF = V->swapFrontier(nodeNum, next);
    next = newF;
    next->clearFrontier();
    //V->registerArr(nodeNum, tmp, V->getSize(nodeNum));
}

void clearLocalFrontier(LocalFrontier *next, int nodeNum, int subNum, int totalSub);

// decides on sparse or dense base on number of nonzeros in the active vertices
template<class F, class vertex>
void edgeMap(graph<vertex> GA, vertices *V, F f, LocalFrontier *next, intT threshold = -1,
             char option = DENSE, bool remDups = false, bool part = false,
             Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    uintT numEdges = GA.m;
    vertex *G = GA.V;
    long long m = (long long) V->numNonzeros() + V->getEdgeStat();

    if (subworker.isMaster()) {
        //printf("%d %d\n", V->numNonzeros(), threshold);
    }

    /*
    if (subworker.isMaster())
	printf("%d\n", m);
    */
    int start = subworker.dense_start;
    int end = subworker.dense_end;

    if (m >= threshold) {
        //Dense part
        if (subworker.isMaster()) {
            printf("Dense: %d\n", m);
            V->toDense();
        }

        if (subworker.isSubMaster()) {
            next->sparseCounter = 0;
        }

        clearLocalFrontier(next, subworker.tid, subworker.subTid, subworker.numOfSub);

        //pthread_barrier_wait(subworker.global_barr);
        subworker.globalWait();

        bool *R = (option == DENSE_FORWARD) ?
                  edgeMapDenseForward(GA, V, f, next, part, start, end) :
                  //edgeMapDenseForwardDynamic(GA, V, f, next, subworker) :
                  //edgeMapDense(GA, V, f, next, option, subworker);
                  edgeMapDenseDynamic(GA, V, f, next, subworker);
        next->isDense = true;
    } else {
        //Sparse part
        if (subworker.isMaster()) {
            printf("Sparse: %d %d\n", V->numNonzeros(), m);
            V->toSparse();
        }
        /*
        if (subworker.isSubMaster()) {
            if (next->sparseChunks == NULL) {
            next->sparseChunks = (intT **)malloc(sizeof(intT *) * subworker.numOfSub);
            }
            if (next->chunkSizes == NULL) {
            next->chunkSizes = (intT *)malloc(sizeof(intT) * subworker.numOfSub);
            }
        }
        */
        //pthread_barrier_wait(subworker.global_barr);
        subworker.globalWait();
        if (V->firstSparse && subworker.isMaster()) {
            printf("my first sparse\n");
        }

        edgeMapSparseV3(GA, V, f, next, part, subworker);
        //edgeMapSparseV4(GA, V, f, next, V->firstSparse, subworker);
        //edgeMapSparseV5(GA, V, f, next, subworker);
        next->isDense = false;
    }
}

//*****VERTEX FUNCTIONS*****
template<class vertex>
void vertexCounter(graph<vertex> GA, LocalFrontier *frontier, int nodeNum, int subNum, int totalSub) {
    if (!frontier->isDense)
        return;

    int size = frontier->endID - frontier->startID;
    int offset = frontier->startID;
    bool *b = frontier->b;
    int subSize = size / totalSub;
    int startPos = subSize * subNum;
    int endPos = subSize * (subNum + 1);
    if (subNum == totalSub - 1) {
        endPos = size;
    }

    int m = 0;
    intT outEdges = 0;

    for (int i = startPos; i < endPos; i++) {
        if (b[i]) {
            outEdges += GA.V[i + offset].getOutDegree();
            m++;
        }
    }
    __sync_fetch_and_add(&(frontier->m), m);
    __sync_fetch_and_add(&(frontier->outEdgesCount), outEdges);
}

template<class F>
void vertexMap(vertices *V, F add, int nodeNum) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    for (int i = 0; i < size; i++) {
        if (b[i])
            add(i + offset);
    }
}

template<class F>
void vertexMap(vertices *V, F add, int nodeNum, int subNum, int totalSub) {
    if (V->isDense) {
        int size = V->getSize(nodeNum);
        int offset = V->getOffset(nodeNum);
        bool *b = V->getArr(nodeNum);
        int subSize = size / totalSub;
        int startPos = subSize * subNum;
        int endPos = subSize * (subNum + 1);
        if (subNum == totalSub - 1) {
            endPos = size;
        }

        for (int i = startPos; i < endPos; i++) {
            if (b[i])
                add(i + offset);
        }
    } else {
        int size = V->frontiers[nodeNum]->m;
        intT *s = V->frontiers[nodeNum]->s;
        int subSize = size / totalSub;
        int startPos = subSize * subNum;
        int endPos = subSize * (subNum + 1);
        if (subNum == totalSub - 1) {
            endPos = size;
        }
        for (int i = startPos; i < endPos; i++) {
            add(s[i]);
        }
    }
}

void clearLocalFrontier(LocalFrontier *next, int nodeNum, int subNum, int totalSub) {
    int size = next->endID - next->startID;
    //int offset = V->getOffset(nodeNum);
    bool *b = next->b;
    int subSize = size / totalSub;
    int startPos = subSize * subNum;
    int endPos = subSize * (subNum + 1);
    if (subNum == totalSub - 1) {
        endPos = size;
    }

    for (int i = startPos; i < endPos; i++) {
        b[i] = false;
    }
}

template<class F>
void vertexFilter(vertices *V, F filter, int nodeNum, bool *result) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    for (int i = 0; i < size; i++) {
        result[i] = false;
        if (b[i])
            result[i] = filter(i + offset);
    }
}

template<class F>
void vertexFilter(vertices *V, F filter, int nodeNum, int subNum, int totalSub, LocalFrontier *result) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    int subSize = size / totalSub;
    int startPos = subSize * subNum;
    int endPos = subSize * (subNum + 1);
    if (subNum == totalSub - 1) {
        endPos = size;
    }

    bool *dst = result->b;
    int m = 0;
    /*
    if (size != result->endID - result->startID || offset != result->startID)
	printf("oops\n");
    */
    for (int i = startPos; i < endPos; i++) {
        //result->setBit(i+offset, b[i] ? (filter(i+offset)) : (false));
        if (b[i]) {
            dst[i] = filter(i + offset);
            if (dst[i])
                m++;
        } else {
            dst[i] = false;
        }
    }
    //printf("filter over\n");
    writeAdd(&(result->m), m);
}

