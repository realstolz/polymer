// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <sys/mman.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "IO.h"

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
		: "m" (*(volatile char *)ptr), "0" (newv)
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
    Subworker_Partitioner(int nSub):numOfSub(nSub){}
    
    inline bool isMaster() {return (tid + subTid == 0);}
    inline bool isSubMaster() {return (subTid == 0);}
    inline intT getStartPos(intT m) {return subTid * (m / numOfSub);}
    inline intT getEndPos(intT m) {return (subTid == numOfSub - 1) ? m : ((subTid + 1) * (m / numOfSub));}
};

struct Default_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;
    Default_Hash_F(int _n, int _shardNum):n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum){}
    
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

template <class vertex>
void partitionByDegree(wghGraph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, bool useOutDegree=false) {
    const intT n = GA.n;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useOutDegree) {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();}
    } else {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree();}
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
    for (intT i = 0; i < n; i+=PAGESIZE/sizeOfOneEle) {
	for (intT j = 0; j < PAGESIZE / sizeOfOneEle; j++) {
	    if (i + j >= n)
		break;
	    accum[counter] += degrees[i + j];
	    sizeArr[counter]++;
	    tmpSizeCounter++;
	}
	if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
	    counter++;
	    //cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle) << endl;
	    tmpSizeCounter = 0;
	}
    }
    
    free(degrees);
}

template <class vertex>
void subPartitionByDegree(wghGraph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, bool useOutDegree=false, bool useFakeDegree=false) {
    const intT n = GA.n;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useFakeDegree) {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getFakeDegree();}
    } else {
	if (useOutDegree) {
	    {parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();}
	} else {
	    {parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree();}
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

template <class vertex>
void subPartitionByDegree(wghGraph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, int subStart, int subEnd, bool useOutDegree=false, bool useFakeDegree=false) {
    const intT n = subEnd - subStart;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useFakeDegree) {
	{parallel_for(intT i = subStart; i < subEnd; i++) degrees[i-subStart] = GA.V[i].getFakeDegree();}
    } else {
	if (useOutDegree) {
	    {parallel_for(intT i = subStart; i < subEnd; i++) degrees[i-subStart] = GA.V[i].getOutDegree();}
	} else {
	    {parallel_for(intT i = subStart; i < subEnd; i++) degrees[i-subStart] = GA.V[i].getInDegree();}
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

template <class vertex, class Hash_F>
void graphHasher(wghGraph<vertex> &GA, Hash_F hash) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *)malloc(sizeof(vertex) * GA.n);

    {parallel_for (intT i = 0; i < GA.n; i++) {
	    intT d = V[i].getOutDegree();
	    //V[i].setFakeDegree(d);
	    intE *outEdges = V[i].getOutNeighborPtr();
	    for (intT j = 0; j < d; j++) {
		outEdges[2*j] = hash.hashFunc(outEdges[j]);
	    }
	    newVertexSet[hash.hashFunc(i)] = V[i];	    
	}
    }
    GA.V = newVertexSet;
    free(V);
}

template <class vertex, class Hash_F>
void graphInEdgeHasher(wghGraph<vertex> &GA, Hash_F hash) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *)malloc(sizeof(vertex) * GA.n);

    {parallel_for (intT i = 0; i < GA.n; i++) {
	    intT d = V[i].getInDegree();
	    //V[i].setFakeDegree(d);
	    intE *inEdges = V[i].getInNeighborPtr();
	    for (intT j = 0; j < d; j++) {
		inEdges[2*j] = hash.hashFunc(inEdges[j]);
	    }
	    newVertexSet[hash.hashFunc(i)] = V[i];	    
	}
    }
    GA.V = newVertexSet;
    free(V);
}

template <class vertex>
wghGraph<vertex> graphFilter(wghGraph<vertex> &GA, int rangeLow, int rangeHi, bool useOutEdge=true) {
    vertex *V = GA.V;
    vertex *newVertexSet = (vertex *)numa_alloc_local(sizeof(vertex) * GA.n);
    int *counters = (int *)numa_alloc_local(sizeof(int) * GA.n);
    int *offsets = (int *)numa_alloc_local(sizeof(int) * GA.n);
    {parallel_for (intT i = 0; i < GA.n; i++) {
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

    intE *edges = (intE *)numa_alloc_local(sizeof(intE) * totalSize * 2);

    {parallel_for (intT i = 0; i < GA.n; i++) {
	    intE *localEdges = &edges[offsets[i]];
	    intT counter = 0;
	    intT d = (useOutEdge) ? (V[i].getOutDegree()) : (V[i].getInDegree());
	    for (intT j = 0; j < d; j++) {
		intT ngh = (useOutEdge) ? (V[i].getOutNeighbor(j)) : (V[i].getInNeighbor(j));
		intT wgh = (useOutEdge) ? (V[i].getOutWeight(j)) : (V[i].getInWeight(j));
		if (rangeLow <= ngh && ngh < rangeHi) {
		    localEdges[counter * 2] = ngh;
		    localEdges[counter * 2 + 1] = wgh;
		    counter++;
		}
	    }
	    if (counter != newVertexSet[i].getFakeDegree()) {
		printf("oops: %d %d\n", counter, newVertexSet[i].getFakeDegree());
	    }
	    if (i == 0) {
		printf("fake deg: %d\n", newVertexSet[i].getFakeDegree());
	    }
	    if (useOutEdge)
		newVertexSet[i].setOutNeighbors(localEdges);
	    else
		newVertexSet[i].setInNeighbors(localEdges);
	}
    }
    numa_free(offsets, sizeof(int) * GA.n);
    //printf("degree: %d\n", newVertexSet[0].getFakeDegree());
    return wghGraph<vertex>(newVertexSet, GA.n, GA.m);
}

void *mapDataArray(int numOfShards, int *sizeArr, int sizeOfOneEle) {
    int numOfPages = 0;
    for (int i = 0; i < numOfShards; i++) {	
        numOfPages += sizeArr[i] / (double)(PAGESIZE / sizeOfOneEle);
    }
    numOfPages++;

    void *toBeReturned = mmap(NULL, numOfPages * PAGESIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (toBeReturned == NULL) {
	cout << "OOps" << endl;
    }
    
    int offset = 0;
    for (int i = 0; i < numOfShards; i++) {
	void *startPos = (void *)((char *)toBeReturned + offset * sizeOfOneEle);
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
    intT emptySignal;
    intT outEdgesCount;
    int startID;
    int endID;
    bool *b;
    intT *s;
    intT *tmp;
    bool isDense;
    
    LocalFrontier(bool *_b, int start, int end):b(_b), startID(start), endID(end), n(end - start), m(0), isDense(true), s(NULL), outEdgesCount(0){}
    
    bool inRange(int index) { return (startID <= index && index < endID);}
    inline void setBit(int index, bool val) { b[index-startID] = val;}
    inline bool getBit(int index) { return b[index-startID];}

    void toSparse() {
	if (isDense) {
	    if (s != NULL)
		free(s);
	    _seq<intT> R = sequence::packIndex(b, n);
	    s = R.A;
	    m = R.n;
	    {parallel_for (intT i = 0; i < m; i++) s[i] = s[i] + startID;}
	    if (m == 0) {
		printf("%p\n", s);
	    } else {
		printf("M is %d and first ele is %d\n", m, s[0]);
	    }
	}
	isDense = false;
    }
    
    void toDense() {
	if (!isDense) {
	    {parallel_for(intT i=0;i<n;i++) b[i] = false;}
	    {parallel_for(intT i=0;i<m;i++) b[s[i] - startID] = true;}
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
    bool** d;
    LocalFrontier **frontiers;
    bool isDense;
    AsyncChunk **asyncQueue;
    int asyncEndSignal;
    intT readerTail;
    intT insertTail;
    
    vertices(int _numOfNodes) {
	this->numOfNodes = _numOfNodes;
	d = (bool **)malloc(numOfNodes * sizeof(bool*));
	frontiers = (LocalFrontier **)malloc(numOfNodes * sizeof(LocalFrontier*));
	numOfVertexOnNode = (int *)malloc(numOfNodes * sizeof(int));
	offsets = (int *)malloc((numOfNodes + 1) * sizeof(int));
	numOfNonZero = (int *)malloc(numOfNodes * sizeof(int));
	numOfVertices = 0;
	m = -1;
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
	for (int i = 0; i< numOfNodes; i++) {
	    //printf("todense %d start m = %d n = %d\n", i, frontiers[i]->s[frontiers[i]->m - 1], frontiers[i]->n);
	    frontiers[i]->toDense();
	}
    }

    void toSparse() {
	if (!isDense)
	    return;
	else
	    isDense = false;
	for (int i = 0; i< numOfNodes; i++) {
	    if (frontiers[i]->isDense) {
		//printf("real convert\n");
	    }
	    frontiers[i]->toSparse();
	}
    }

    intT getEdgeStat() {
	intT sum = 0;
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
	    offsets[i] = numOfVertexOnNode[i-1] + offsets[i-1];
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
	LocalFrontier *oldOne =  frontiers[nodeNum];
	frontiers[nodeNum] = newOne;
	m = -1;
	return oldOne;
    }

    bool eq (vertices& b) {
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
};

struct nonNegF{bool operator() (intT a) {return (a>=0);}};

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options {
  DENSE, DENSE_PARALLEL, DENSE_FORWARD
};




//*****EDGE FUNCTIONS*****

template <class F, class vertex>
bool* edgeMapDense(wghGraph<vertex> GA, vertices* frontier, F f, LocalFrontier *next, bool parallel = 0, Subworker_Partitioner &subworker = 0) {
    intT numVertices = GA.n;
    intT size = next->endID - next->startID;
    vertex *G = GA.V;
    intT startPos = subworker.dense_start;
    intT endPos = subworker.dense_end;
    bool *bitVec = frontier->getArr(subworker.tid);
    for (intT i = startPos; i < endPos; i++){
	next->setBit(i, false);
	if (f.cond(i) && bitVec[i - next->startID]) { 
	    intT d = G[i].getInDegree();
	    for(intT j=0; j<d; j++){
		intT ngh = G[i].getInNeighbor(j);
		if (f.update(ngh, i, G[i].getInWeight(j))) next->setBit(i, true);
		if(!f.cond(i)) break;
		__builtin_prefetch(f.nextPrefetchAddr(G[i].getInNeighbor(j+3)), 1, 3);
	    }
	}
    }
    return NULL;
}

template <class F, class vertex>
bool* edgeMapDenseForward(wghGraph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false, int start = 0, int end = 0) {
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
	nextSwitchPoint = frontier->getOffset(currNodeNum+1);
	currOffset = frontier->getOffset(currNodeNum);
    }
    for (long i=startPos; i<endPos; i++){
	if (i == nextSwitchPoint) {
	    currOffset += frontier->getSize(currNodeNum);
	    nextSwitchPoint += frontier->getSize(currNodeNum + 1);
	    currNodeNum++;
	    currBitVector = frontier->getArr(currNodeNum);
	    //printf("OK\n");
	}
	//printf("edgemap: %p\n", currBitVector);
	m += G[i].getFakeDegree();
	if (currBitVector[i-currOffset]) {
	    intT d = G[i].getFakeDegree();
	    for(intT j=0; j<d; j++){
		uintT ngh = G[i].getOutNeighbor(j);
		if (/*next->inRange(ngh) &&*/ f.cond(ngh) && f.updateAtomic(i, ngh, G[i].getOutWeight(j))) {
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

template <class F, class vertex>
bool* edgeMapDenseForwardGlobalWrite(wghGraph<vertex> GA, vertices *frontier, F f, LocalFrontier *nexts[], Subworker_Partitioner &subworker) {
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
    intT nextSwitchPoint = frontier->getOffset(currNodeNum+1);
    int offset = frontier->getOffset(currNodeNum);

    for (long i=startPos; i<endPos; i++) {
	if (i == nextSwitchPoint) {
	    offset += frontier->getSize(currNodeNum);
	    nextSwitchPoint += frontier->getSize(currNodeNum + 1);
	    currNodeNum++;
	    nextBitVector = nexts[currNodeNum]->b;
	}
	m += G[i].getFakeDegree();
	if (f.cond(i)) {
	    intT d = G[i].getFakeDegree();
	    for(intT j=0; j<d; j++) {
		uintT ngh = G[i].getInNeighbor(j);
		if (currBitVector[ngh-currOffset] && f.updateAtomic(ngh, i, G[i].getInWeight(j))) {
		    nextBitVector[i-offset] = true;
		}
	    }
	}
    }
    return NULL;
}

AsyncChunk *newChunk(int blockSize) {
    AsyncChunk *myChunk = (AsyncChunk *)malloc(sizeof(AsyncChunk));
    myChunk->s = (intT *)malloc(sizeof(intT) * blockSize);
    myChunk->m = 0;
    myChunk->accessCounter = 0;
    return myChunk;
}

template <class F, class vertex>
void edgeMapSparseAsync(wghGraph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, Subworker_Partitioner &subworker = NULL) {
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
	} while (!__sync_bool_compare_and_swap((intT *)queueHead, currHead, endPos));
	
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
		printf("oops: %p %p %d %d %d\n", currChunk, frontier->asyncQueue[currHead % GA.n], currHead, currTail, endPos);
	    }
	    //printf("chunk pointer: %p\n", currChunk);
	    int chunkSize = currChunk->m;
	    for (intT i = 0; i < chunkSize; i++) {
		accumSize++;
		intT idx = currChunk->s[i];
		intT d = V[idx].getOutDegree();
		for (intT j = 0; j < d; j++) {
		    intT ngh = V[idx].getOutNeighbor(j);
		    if (f.cond(ngh) && f.updateAtomic(idx, ngh, V[idx].getOutWeight(j))) {
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
			    while (!__sync_bool_compare_and_swap((intT *)queueTail, insertPos, insertPos+1)) {				
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
		while (!__sync_bool_compare_and_swap((intT *)queueTail, insertPos, insertPos+1)) {
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

template <class F, class vertex>
void edgeMapSparseV3(wghGraph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false, Subworker_Partitioner &subworker = NULL) {
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
	    next->s = (intT *)malloc(sizeof(intT) * bufferLen);
	intT nextEdgesCount = 0;
	
	pthread_barrier_wait(subworker.local_barr);
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
		for (intT j = 0; j < d; j++) {
		    uintT ngh = V[idx].getOutNeighbor(j);
		    if (f.cond(ngh) && f.updateAtomic(idx, ngh, V[idx].getOutWeight(j))) {
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
	pthread_barrier_wait(subworker.local_barr);
    }
}

static int edgesTraversed = 0;

void switchFrontier(int nodeNum, vertices *V, LocalFrontier* &next) {
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

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
void edgeMap(wghGraph<vertex> GA, vertices *V, F f, LocalFrontier *next, intT threshold = -1, 
	     char option=DENSE, bool remDups=false, bool part = false, Subworker_Partitioner &subworker = NULL) {
    intT numVertices = GA.n;
    uintT numEdges = GA.m;
    vertex *G = GA.V;    
    intT m = V->numNonzeros();// + V->getEdgeStat();
    
    /*
    if (subworker.isMaster())
	printf("%d\n", m);
    */
    int start = subworker.dense_start;
    int end = subworker.dense_end;

    if (m >= threshold) {
	//Dense part
	if (subworker.isMaster()) {
	    V->toDense();
	}

	pthread_barrier_wait(subworker.global_barr);
	
	bool* R = (option == DENSE_FORWARD) ? 
	    edgeMapDenseForward(GA, V, f, next, part, start, end) : 
	    edgeMapDense(GA, V, f, next, option, subworker);
	next->isDense = true;
    } else {
	//Sparse part
	if (subworker.isMaster()) {
	    V->toSparse();
	}

	pthread_barrier_wait(subworker.global_barr);
	
	edgeMapSparseV3(GA, V, f, next, part, subworker);
	next->isDense = false;
    }
}

//*****VERTEX FUNCTIONS*****
template<class vertex>
void vertexCounter(wghGraph<vertex> GA, LocalFrontier *frontier, int nodeNum, int subNum, int totalSub) {
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
	    outEdges += GA.V[i+offset].getOutDegree();
	    m++;
	}
    }
    writeAdd(&(frontier->m), m);
    writeAdd(&(frontier->outEdgesCount), outEdges);
}

template <class F>
void vertexMap(vertices *V, F add, int nodeNum) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    for (int i = 0; i < size; i++) {
	if (b[i])
	    add(i + offset);
    }
}

template <class F>
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

template <class F>
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

template <class F>
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

//*****WEIGHTED EDGE FUNCTIONS*****
/*
template <class F, class vertex>
  bool* edgeMapDense(wghGraph<vertex> GA, bool* vertices, F f, bool parallel = 0) {
  intT numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for (long i=0; i<numVertices; i++){
    next[i] = 0;
    if (f.cond(i)) {
      intT d = G[i].getInDegree();
      if(!parallel || d < 1000) {
	for(intT j=0; j<d; j++){
	  intT ngh = G[i].getInNeighbor(j);
	  if(vertices[ngh]){
	    intT edgeLen = G[i].getInWeight(j);
	    if (vertices[ngh] && f.update(ngh,i,edgeLen)) next[i] = 1;
	  }
	  if (!f.cond(i)) break;
	}
      }
      else {
	{parallel_for(intT j=0; j<d; j++){
	  intT ngh = G[i].getInNeighbor(j);
	  if(vertices[ngh]){
	    intT edgeLen = G[i].getInWeight(j);
	    if (vertices[ngh] && f.update(ngh,i,edgeLen)) next[i] = 1;
	  }
	  }}
      }
    }
    }
  }
  return next;
}

template <class F, class vertex>
bool* edgeMapDenseForward(wghGraph<vertex> GA, bool* vertices, F f) {
  intT numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for(long i=0;i<numVertices;i++) next[i] = 0;}
  {parallel_for (long i=0; i<numVertices; i++){
    if (vertices[i]) {
      intT d = G[i].getOutDegree();
      if(d < 1000) {
	for(intT j=0; j<d; j++){
	  intT ngh = G[i].getOutNeighbor(j);
	  intT edgeLen = G[i].getOutWeight(j);
	  if (f.cond(ngh) && f.updateAtomic(i,ngh,edgeLen)) next[ngh] = 1;
	}
      }
      else{
	parallel_for(intT j=0; j<d; j++){
	  intT ngh = G[i].getOutNeighbor(j);
	  intT edgeLen = G[i].getOutWeight(j);
	  if (f.cond(ngh) && f.updateAtomic(i,ngh,edgeLen)) next[ngh] = 1;
	}
      }
    }
  }}
  return next;
}

template <class F, class vertex>
  pair<uintT,intT*> edgeMapSparseW(vertex* frontierVertices, intT* indices, 
				   uintT* degrees, uintT m, F f, 
				   intT remDups=0, intT* flags=NULL) {
  uintT* offsets = degrees;
  uintT outEdgeCount = sequence::plusScan(offsets, degrees, m);

  intT* outEdges = newA(intT,outEdgeCount);

  {parallel_for (long i = 0; i < m; i++) {
    intT v = indices[i];
    uintT o = offsets[i];
    vertex vert = frontierVertices[i];
    intT d = vert.getOutDegree();
    if(d < 1000) {
      for (intT j=0; j < d; j++) {
	intT ngh = vert.getOutNeighbor(j);
	intT edgeLen = vert.getOutWeight(j);
	if (f.cond(ngh) && f.updateAtomic(v,ngh,edgeLen)) 
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = -1;
      }
    }
    else {
      {parallel_for (intT j=0; j < d; j++) {
	intT ngh = vert.getOutNeighbor(j);
	intT edgeLen = vert.getOutWeight(j);
	if (f.cond(ngh) && f.updateAtomic(v,ngh,edgeLen)) 
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = -1;
	
	}}
    }
    }}

  intT* nextIndices = newA(intT, outEdgeCount);
  if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);

  // Filter out the empty slots (marked with -1)
  uintT nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonNegF());
  free(outEdges);
  return pair<uintT,intT*>(nextM, nextIndices);
}

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
  vertices edgeMap(wghGraph<vertex> GA, vertices V, F f, intT threshold = -1, char option=DENSE, bool remDups=false) {
  intT numVertices = GA.n;
  uintT numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  intT m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }

  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(vertex,m);

  {parallel_for (long i=0; i < m; i++){
    vertex v = G[V.s[i]];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  uintT outDegrees = sequence::plusReduce(degrees, m);    
  if (outDegrees == 0) return vertices(numVertices);
  if (m + outDegrees > threshold) { 
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R = option == DENSE_FORWARD ? 
      edgeMapDenseForward(GA,V.d,f) : 
      edgeMapDense(GA, V.d, f,option);
    vertices v1 = vertices(numVertices, R);
    //cout << "size (D) = " << v1.m << endl;
    return  v1;
  } else { 
    pair<uintT,intT*> R = 
      remDups ? 
      edgeMapSparseW(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparseW(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    //cout << "size (S) = " << R.first << endl;
    free(degrees);
    free(frontierVertices);
    return vertices(numVertices, R.first, R.second);
  }
}
*/
