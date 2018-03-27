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

#include "polymer.h"
#include "gettime.h"
#include "math.h"

#include <pthread.h>
#include <sys/mman.h>
#include <numa.h>
using namespace std;

#define PAGE_SIZE (4096)

int CORES_PER_NODE = 6;

volatile int shouldStart = 0;

int vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

pthread_barrier_t barr;
pthread_barrier_t global_barr;
pthread_mutex_t mut;

volatile int global_counter = 0;
volatile int global_toggle = 0;

void *graph_full;

vertices *Frontier;

#define NSTATES 2

struct EdgeWeight {
    float potential[NSTATES][NSTATES];
};

struct EdgeData {
    float belief[NSTATES];
};

struct VertexInfo {
    float potential[NSTATES];
};

struct VertexData {
    float product[NSTATES];
};

template <class ET>
inline void writeDiv(ET *a, ET b) {
  volatile ET newV, oldV; 
  do {oldV = *a; newV = oldV / b;}
  while (!CAS(a, oldV, newV));
}

template <class ET>
inline void writeMult(ET *a, ET b) {
  volatile ET newV, oldV; 
  do {oldV = *a; newV = oldV * b;}
  while (!CAS(a, oldV, newV));
}

template <class vertex>
struct BP_F {
    EdgeWeight *edgeW;
    EdgeData *edgeD_curr;
    EdgeData *edgeD_next;
    VertexInfo *vertI;
    VertexData *vertD_curr;
    VertexData *vertD_next;
    intT *offsets;
    intT rangeLow;
    BP_F(EdgeWeight *_edgeW, EdgeData *_edgeD_curr, EdgeData *_edgeD_next, VertexInfo *_vertI, VertexData *_vertD_curr, VertexData *_vertD_next, intT *_offsets, intT _rangeLow=0) : 
	edgeW(_edgeW), edgeD_curr(_edgeD_curr), edgeD_next(_edgeD_next), vertI(_vertI), vertD_curr(_vertD_curr), vertD_next(_vertD_next), offsets(_offsets),rangeLow(_rangeLow) {}
    inline bool update(intT s, intT d, intT edgeIdx){
	intT dstIdx = offsets[s] + edgeIdx;
	for (int i = 0; i < NSTATES; i++) {
	    edgeD_next[dstIdx].belief[i] = 0.0;
	    for (int j = 0; j < NSTATES; j++) {
		edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
	    }
	    vertD_next[d].product[i] = vertD_next[d].product[i] * edgeD_next[dstIdx].belief[i];
	}
	return 1;
    }
    inline bool updateAtomic (intT s, intT d, intT edgeIdx) { //atomic Update
	//printf("we are here: %d\n", s);
	intT dstIdx = offsets[s-rangeLow] + edgeIdx;
	//printf("idx: %d\n", dstIdx);
	for (int i = 0; i < NSTATES; i++) {
	    edgeD_next[dstIdx].belief[i] = 0.0;
	    for (int j = 0; j < NSTATES; j++) {
		edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
	    }
	    //writeMult(&(vertD_next[d].product[i]), edgeD_next[dstIdx].belief[i]);
	    vertD_next[d].product[i] = vertD_next[d].product[i] * edgeD_next[dstIdx].belief[i];
	}
	return 1;
    }
    inline bool cond (intT d) {return 1; } //does nothing
};

//resets p
struct BP_Vertex_Reset {
    VertexData *vertD;
    BP_Vertex_Reset(VertexData *_vertD) :
	vertD(_vertD) {}
    inline bool operator () (intT i) {
	for (int i = 0; i < NSTATES; i++) {
	    vertD[i].product[i] = 1.0;
	}
	return 1;
    }
};

struct BP_worker_arg {
    void *GA;
    int maxIter;
    int tid;
    int numOfNode;
    int rangeLow;
    int rangeHi;
    
    VertexInfo *vertI;
    VertexData *vertD_curr;
    VertexData *vertD_next;
};

struct BP_subworker_arg {
    void *GA;
    int maxIter;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    pthread_barrier_t *node_barr;
    LocalFrontier *localFrontier;
    volatile int *barr_counter;
    volatile int *toggle;

    VertexInfo *vertI;
    VertexData *vertD_curr;
    VertexData *vertD_next;
    EdgeWeight *edgeW;
    EdgeData *edgeD_curr;
    EdgeData *edgeD_next;
    intT *localOffsets;
    intT *localOffsets2;
};

template <class F, class vertex>
bool* edgeMapDenseBPNoRep(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false, Subworker_Partitioner &subworker=dummyPartitioner) {
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
    
    int size = frontier->getSize(subworker.tid);
    int subSize = size / CORES_PER_NODE;
    intT startPos = subSize * subworker.subTid;
    intT endPos = subSize * (subworker.subTid + 1);
    if (subworker.subTid == CORES_PER_NODE - 1) {
	endPos = size;
    }

    intT low = frontier->getFrontier(subworker.tid)->startID;
    startPos += low;
    endPos += low;

    currNodeNum = frontier->getNodeNumOfIndex(startPos);
    currBitVector = frontier->getArr(currNodeNum);
    nextSwitchPoint = frontier->getOffset(currNodeNum+1);
    currOffset = frontier->getOffset(currNodeNum);

    for (long i=startPos; i<endPos; i++){
	if (i == nextSwitchPoint) {
	    currOffset += frontier->getSize(currNodeNum);
	    nextSwitchPoint += frontier->getSize(currNodeNum + 1);
	    currNodeNum++;
	    currBitVector = frontier->getArr(currNodeNum);
	}
	m += G[i].getOutDegree();
	if (currBitVector[i-currOffset]) {
	    intT d = G[i].getFakeDegree();
	    for(intT j=0; j<d; j++){
		uintT ngh = G[i].getOutNeighbor(j);
		if (/*next->inRange(ngh) &&*/ f.cond(ngh) && f.updateAtomic(i,ngh,j)) {
		    next->setBit(ngh, true);
		}
	    }
	}
    }
    return NULL;
}

template <class vertex>
void *BeliefPropagationSubWorker(void *arg) {
    BP_subworker_arg *my_arg = (BP_subworker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    const intT n = GA.n;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;
    pthread_barrier_t *local_barr = my_arg->node_barr;
    LocalFrontier *output = my_arg->localFrontier;

    int currIter = 0;
    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int start = my_arg->startPos;
    int end = my_arg->endPos;

    VertexInfo *vertI = my_arg->vertI;
    VertexData *vertD_curr = my_arg->vertD_curr;
    VertexData *vertD_next = my_arg->vertD_next;

    EdgeWeight *edgeW = my_arg->edgeW;
    EdgeData *edgeD_curr = my_arg->edgeD_curr;
    EdgeData *edgeD_next = my_arg->edgeD_next;

    intT *localOffsets = my_arg->localOffsets;
    intT *localOffsets2 = my_arg->localOffsets2;
    graph<vertex> &Global_G = *(graph<vertex> *)graph_full;

    Custom_barrier globalCustom(&global_counter, &global_toggle, Frontier->numOfNodes);
    Custom_barrier localCustom(my_arg->barr_counter, my_arg->toggle, CORES_PER_NODE);

    Subworker_Partitioner subworker(CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = &global_barr;
    subworker.local_custom = localCustom;
    subworker.subMaster_custom = globalCustom;

    if (subTid == 0) {
	Frontier->getFrontier(tid)->m = rangeHi - rangeLow;
    }

    pthread_barrier_wait(local_barr);
    pthread_barrier_wait(&global_barr);
    while(1) {
	if (maxIter > 0 && currIter >= maxIter)
            break;
        currIter++;
	if (subTid == 0)
	    Frontier->calculateNumOfNonZero(tid);
	if (subTid == 0) {
	    {parallel_for(long i=output->startID;i<output->endID;i++) output->setBit(i, false);}
	}
	
	pthread_barrier_wait(&global_barr);

	vertexMap(Frontier, BP_Vertex_Reset(vertD_next), tid, subTid, CORES_PER_NODE);
	output->m = 1;
	pthread_barrier_wait(&global_barr);	

        //edgeMapDenseBP(GA, Frontier, BP_F<vertex>(edgeW, edgeD_curr, edgeD_next, vertI, vertD_curr, vertD_next, localOffsets),output,true,start,end);
        edgeMapDenseBPNoRep(GA, Frontier, BP_F<vertex>(edgeW, edgeD_curr, edgeD_next, vertI, vertD_curr, vertD_next, localOffsets2,rangeLow),output,true,subworker);
	pthread_barrier_wait(&global_barr);

	pthread_barrier_wait(&global_barr);
	//pthread_barrier_wait(local_barr);

	swap(edgeD_curr, edgeD_next);
	swap(vertD_curr, vertD_next);
	
	if (subworker.isSubMaster()) {
	    pthread_barrier_wait(&global_barr);
	    switchFrontier(tid, Frontier, output);
	} else {
	    output = Frontier->getFrontier(tid);
	    pthread_barrier_wait(&global_barr);
	}
	
	pthread_barrier_wait(&global_barr);
	//pthread_barrier_wait(local_barr);
    }

    pthread_barrier_wait(local_barr);
    return NULL;
}

pthread_barrier_t timerBarr;

template <class vertex>
void *BeliefPropagationThread(void *arg) {
    BP_worker_arg *my_arg = (BP_worker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    graph<vertex> localGraph = graphFilter(GA, rangeLow, rangeHi);

    // create edge data

    intT *fakeDegrees = (intT *)numa_alloc_local(sizeof(intT) * localGraph.n);
    intT *localOffsets = (intT *)numa_alloc_local(sizeof(intT) * localGraph.n);

    {parallel_for (intT i = 0; i < localGraph.n; i++) {
	    fakeDegrees[i] = localGraph.V[i].getFakeDegree();
	}
    }

    localOffsets[0] = 0;
    for (intT i = 1; i < localGraph.n; i++) {
	localOffsets[i] = localOffsets[i-1] + fakeDegrees[i-1];
    }

    intT numLocalEdge = localOffsets[localGraph.n - 1] + fakeDegrees[localGraph.n - 1];

    intT *localDegrees = (intT *)numa_alloc_local(sizeof(intT) * localGraph.n);
    intT *localOffsets2 = (intT *)numa_alloc_local(sizeof(intT) * localGraph.n);

    {parallel_for (intT i = rangeLow; i < rangeHi; i++) {
	    localDegrees[i-rangeLow] = GA.V[i].getOutDegree();
	}
    }

    localOffsets2[0] = 0;
    for (intT i = 1; i < (rangeHi - rangeLow); i++) {
	localOffsets2[i] = localOffsets2[i-1] + localDegrees[i-1];
    }

    numLocalEdge = localOffsets2[rangeHi - rangeLow - 1] + localDegrees[rangeHi - rangeLow - 1];
    printf ("numLocalEdge of %d: %d\n", tid, numLocalEdge);

    EdgeWeight *edgeW = (EdgeWeight *)numa_alloc_local(sizeof(EdgeWeight) * numLocalEdge);
    
    EdgeData *edgeD_curr = (EdgeData *)numa_alloc_local(sizeof(EdgeData) * numLocalEdge);
    EdgeData *edgeD_next = (EdgeData *)numa_alloc_local(sizeof(EdgeData) * numLocalEdge);
    /*
    {parallel_for (intT i = 0; i < n; i++) {
	    intT o = offsets[i];
	    intT d = degrees[i];
	    vertex vert = GA.V[i];
	    for (intT j = 0; j < d; j++) {
		intT ngh = vert.getOutNeighbor(j);
		int idx = -1;
		intT ngh_d = GA.V[ngh].getOutDegree();
		vertex ngh_vert = GA.V[ngh];
		for (intT k = 0; k < ngh_d; k++) {
		    if (ngh_vert.getOutNeighbor(k) == i) {
			idx = k;
			break;
		    }
		}
		if (idx == -1) {
		    idx = 0;
		}
		edgeW[o + j].transposeOffset = offsets[ngh] + idx;
	    }
	}
    }
    */
    int sizeOfShards[CORES_PER_NODE];

    subPartitionByDegree(localGraph, CORES_PER_NODE, sizeOfShards, sizeof(VertexData), true, true);
    
    for (int i = 0; i < CORES_PER_NODE; i++) {
	//printf("subPartition: %d %d: %d\n", tid, i, sizeOfShards[i]);
    }

    while (shouldStart == 0) ;

    printf("over filtering\n");
    /*
    if (0 != __cilkrts_set_param("nworkers","1")) {
	printf("set failed: %d\n", tid);
    }
    */

    const intT n = GA.n;
    int numOfT = my_arg->numOfNode;

    int blockSize = rangeHi - rangeLow;

    //printf("blockSizeof %d: %d low: %d high: %d\n", tid, blockSize, rangeLow, rangeHi);

    bool* frontier = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    
    /*
    if (tid == 0)
	startTime();
    */
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    for(intT i=0;i<blockSize;i++) frontier[i] = true;
    if (tid == 0)
	Frontier = new vertices(numOfT);

    //printf("register %d: %p\n", tid, frontier);
    
    LocalFrontier *current = new LocalFrontier(frontier, rangeLow, rangeHi);

    bool* next = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    for(intT i=0;i<blockSize;i++) next[i] = false;
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);
    
    Frontier->registerFrontier(tid, current);

    pthread_barrier_wait(&barr);

    if (tid == 0)
	Frontier->calculateOffsets();

    pthread_barrier_t localBarr;
    pthread_barrier_init(&localBarr, NULL, CORES_PER_NODE+1);

    int startPos = 0;

    pthread_t subTids[CORES_PER_NODE];    

    volatile int local_custom_counter = 0;
    volatile int local_toggle = 0;

    for (int i = 0; i < CORES_PER_NODE; i++) {	
	BP_subworker_arg *arg = (BP_subworker_arg *)malloc(sizeof(BP_subworker_arg));
	arg->GA = (void *)(&localGraph);
	arg->maxIter = maxIter;
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;
	arg->node_barr = &localBarr;
	arg->localFrontier = output;

	arg->barr_counter = &local_custom_counter;
	arg->toggle = &local_toggle;
	
	arg->startPos = startPos;
	arg->endPos = startPos + sizeOfShards[i];

	arg->edgeW = edgeW;
	arg->edgeD_curr = edgeD_curr;
	arg->edgeD_next = edgeD_next;

	arg->vertI = my_arg->vertI;
	arg->vertD_curr = my_arg->vertD_curr;
	arg->vertD_next = my_arg->vertD_next;
	arg->localOffsets = localOffsets;
	arg->localOffsets2 = localOffsets2;
	startPos = arg->endPos;
        pthread_create(&subTids[i], NULL, BeliefPropagationSubWorker<vertex>, (void *)arg);
    }

    pthread_barrier_wait(&barr);
    pthread_barrier_wait(&timerBarr);
    printf("here we go\n");

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&barr);
    intT round = 0;
    return NULL;
}

struct BP_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;
    BP_Hash_F(int _n, int _shardNum):n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum){}
    
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
void BeliefPropagation(graph<vertex> &GA, int maxIter) {
    numOfNode = numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
    CORES_PER_NODE = numa_num_configured_cpus() / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&timerBarr, NULL, numOfNode+1);
    pthread_barrier_init(&global_barr, NULL, CORES_PER_NODE * numOfNode);
    pthread_mutex_init(&mut, NULL);
    int sizeArr[numOfNode];
    BP_Hash_F hasher(GA.n, numOfNode);
    graphHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(VertexData));
    graph_full = (void *)&GA;
    /*
    intT vertPerPage = PAGESIZE / sizeof(double);
    intT subShardSize = ((GA.n / numOfNode) / vertPerPage) * vertPerPage;
    for (int i = 0; i < numOfNode - 1; i++) {
	sizeArr[i] = subShardSize;
    }
    sizeArr[numOfNode - 1] = GA.n - subShardSize * (numOfNode - 1);
    */
    VertexInfo *vertI = (VertexInfo *)malloc(sizeof(VertexInfo) * GA.n);
    VertexData *vertD_curr = (VertexData *)mapDataArray(numOfNode, sizeArr, sizeof(VertexData));
    VertexData *vertD_next = (VertexData *)mapDataArray(numOfNode, sizeArr, sizeof(VertexData));

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    int prev = 0;
    for (int i = 0; i < numOfNode; i++) {
	BP_worker_arg *arg = (BP_worker_arg *)malloc(sizeof(BP_worker_arg));
	arg->GA = (void *)(&GA);
	arg->maxIter = maxIter;
	arg->tid = i;
	arg->numOfNode = numOfNode;
	arg->rangeLow = prev;
	arg->rangeHi = prev + sizeArr[i];
	
	arg->vertI = vertI;
	arg->vertD_curr = vertD_curr;
	arg->vertD_next = vertD_next;
	prev = prev + sizeArr[i];
	pthread_create(&tids[i], NULL, BeliefPropagationThread<vertex>, (void *)arg);
    }
    shouldStart = 1;
    pthread_barrier_wait(&timerBarr);
    nextTime("Graph Partition");
    startTime();
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("BeliefPropagation");
    if (needResult) {

    }
}

int parallel_main(int argc, char* argv[]) {  
    char* iFile;
    bool binary = false;
    bool symmetric = false;
    int maxIter = -1;
    needResult = false;
    if(argc > 1) iFile = argv[1];
    if(argc > 2) maxIter = atoi(argv[2]);
    if(argc > 3) if((string) argv[3] == (string) "-result") needResult = true;
    if(argc > 4) if((string) argv[4] == (string) "-s") symmetric = true;
    if(argc > 5) if((string) argv[5] == (string) "-b") binary = true;
    numa_set_interleave_mask(numa_all_nodes_ptr);
    if(symmetric) {
	graph<symmetricVertex> G = 
	    readGraph<symmetricVertex>(iFile,symmetric,binary);
	BeliefPropagation(G, maxIter);
	//G.del(); 
    } else {
	graph<asymmetricVertex> G = 
	    readGraph<asymmetricVertex>(iFile,symmetric,binary);
	BeliefPropagation(G, maxIter);
	//G.del();
    }
    return 0;
}
