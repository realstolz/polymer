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

#include <numa.h>
#include <pthread.h>
using namespace std;

pthread_barrier_t barr;
pthread_barrier_t subMasterBarr;
pthread_barrier_t global_barr;
pthread_barrier_t timerBarr;

volatile int shouldStart = 0;

volatile int global_counter = 0;
volatile int global_toggle = 0;

intT *IDs_global = NULL;
intT *PrevIDs_global = NULL;

int vPerNode = 0;
int numOfNode = 0;

int CORES_PER_NODE = 0;

bool needResult = false;

vertices *Frontier;

struct CC_F {
    intT* IDs;
    intT* prevIDs;
    CC_F(intT* _IDs, intT* _prevIDs) : 
	IDs(_IDs), prevIDs(_prevIDs) {}
    inline void *nextPrefetchAddr(intT index) {
	return NULL;
    }
    inline bool update(intT s, intT d){ //Update function writes min ID
	intT origID = IDs[d];
	if(IDs[s] < origID) {
	    IDs[d] = min(origID,IDs[s]);
	    if(origID == prevIDs[d]) return 1;
	}
	return 0;
    }
    inline bool updateAtomic (intT s, intT d) { //atomic Update
	intT origID = IDs[d];
	bool res = (writeMin(&IDs[d], IDs[s]) && origID == prevIDs[d]);
	return res;
    }

    inline void initFunc(void *dataPtr, intT d) {
	intT *tmp = (intT *)dataPtr;
	tmp[0] = IDs[d];
	tmp[1] = prevIDs[d];
	return;
    }

    inline bool reduceFunc(void *dataPtr, intT s) {
	intT *tmp = (intT *)dataPtr;
	intT origID = *(intT *)dataPtr;
	if(IDs[s] < origID) {
	    *(intT *)dataPtr = min(origID,IDs[s]);
	    if(origID == tmp[1]) return true;
	}
	return false;
    }

    inline bool combineFunc(void *dataPtr, intT d) {
	intT newVal = *(intT *)dataPtr;
	intT origID = IDs[d];
	bool res = (writeMin(&IDs[d], newVal) && origID == prevIDs[d]);
	return res;
    }
    
    inline void vertUpdate(intT v) {
	prevIDs[v] = IDs[v];
    }

    inline bool cond (intT d) { return 1; } //does nothing
};

//function used by vertex map to sync prevIDs with IDs
struct CC_Vertex_F {
  intT* IDs;
  intT* prevIDs;
  CC_Vertex_F(intT* _IDs, intT* _prevIDs) :
    IDs(_IDs), prevIDs(_prevIDs) {}
  inline bool operator () (intT i) {
    prevIDs[i] = IDs[i];
    //printf("update: %d: %d to %d\n", i, prevIDs[i], IDs[i]);
    return 1;
  }
};

template <class F, class vertex>
void edgeMapCustom(graph<vertex> GA, vertices *V, F f, LocalFrontier *next, intT threshold = -1, 
	     char option=DENSE, bool remDups=false, bool part = false, Subworker_Partitioner &subworker = dummyPartitioner) {
    intT numVertices = GA.n;
    uintT numEdges = GA.m;
    vertex *G = GA.V;    
    long long m = (long long)V->numNonzeros() + V->getEdgeStat();
    
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
	    //printf("Dense: %d\n", m);
	    V->toDense();
	}

	if (subworker.isSubMaster()) {
	    next->sparseCounter = 0;
	}

	clearLocalFrontier(next, subworker.tid, subworker.subTid, subworker.numOfSub);

	//pthread_barrier_wait(subworker.global_barr);
	subworker.globalWait();
	
	bool* R = (option == DENSE_FORWARD) ? 
	    edgeMapDenseForward(GA, V, f, next, part, start, end) :
	    //edgeMapDense(GA, V, f, next, option, subworker);
	    edgeMapDenseReduce(GA, V, f, next, option, subworker);
	next->isDense = true;
    } else {
	//Sparse part
	if (subworker.isMaster()) {
	    //printf("Sparse: %d %d\n", V->numNonzeros(), m);
	    V->toSparse();
	}
	subworker.globalWait();
	if (V->firstSparse && subworker.isMaster()) {
	    printf("my first sparse\n");
	}
	
	edgeMapSparseV3(GA, V, f, next, part, subworker);
	next->isDense = false;
    }
}

template <class vertex>
void *ComponentsSubWorker(void *args) {
    Default_subworker_arg *my_arg = (Default_subworker_arg *)args;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    graph<vertex> &GA_global = *(graph<vertex> *)my_arg->Global_G;
    const intT n = GA.n;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;
    pthread_barrier_t *local_barr = my_arg->node_barr;
    pthread_barrier_t *master_barr = my_arg->master_barr;
    pthread_barrier_t *global_barr = my_arg->global_barr;
    LocalFrontier *output = my_arg->localFrontier;
    
    int currIter = 0;
    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int start = my_arg->startPos;
    int end = my_arg->endPos;

    intT numVisited = 0;

    Custom_barrier local_custom(my_arg->local_custom_counter, my_arg->local_custom_toggle, CORES_PER_NODE);
    Custom_barrier global_custom(&global_counter, &global_toggle, Frontier->numOfNodes);

    Subworker_Partitioner subworker(CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = global_barr;
    subworker.local_barr = my_arg->node_barr;
    subworker.leader_barr = &subMasterBarr;
    subworker.local_custom = local_custom;
    subworker.subMaster_custom = global_custom;

    intT *IDs = IDs_global;
    intT *PrevIDs = PrevIDs_global;
    if (subworker.isMaster()) {
	pthread_barrier_init(&subMasterBarr, NULL, Frontier->numOfNodes);
    }
    pthread_barrier_wait(master_barr);

    if (subworker.isSubMaster()) {
	Frontier->calculateNumOfNonZero(tid);
    }

    pthread_barrier_wait(global_barr);
    
    intT switchThreshold = GA.m/8;

    struct timeval startT, endT;
    struct timezone tz = {0, 0};

    gettimeofday(&startT, &tz);

    while (!Frontier->isEmpty() || currIter == 0) {
	currIter++;
	intT currM = Frontier->numNonzeros();
	if (subworker.isMaster()) {
	    numVisited += currM;
	    //printf("non zeros: %d\n", currM);
	}
	
	//clearLocalFrontier(output, tid, subTid, CORES_PER_NODE);

	vertexMap(Frontier, CC_Vertex_F(IDs,PrevIDs), tid, subTid, CORES_PER_NODE);
	//pthread_barrier_wait(global_barr);
	subworker.globalWait();

	//edgeMap(GA, Frontier, CC_F(IDs,PrevIDs), output, switchThreshold, DENSE_FORWARD, false, true, subworker);
	edgeMapCustom(GA, Frontier, CC_F(IDs,PrevIDs), output, switchThreshold, DENSE_PARALLEL, false, true, subworker);
	/*
	if (currM >= switchThreshold) {
	    edgeMap(GA, Frontier, CC_F(IDs,PrevIDs), output, switchThreshold, DENSE_FORWARD, false, true, subworker);
	} else {
	    printf("here\n");	    
	    edgeMapSparse(GA_global, Frontier, CC_F(IDs, PrevIDs), CC_Vertex_F(IDs, PrevIDs), output, subworker);
	    pthread_barrier_wait(global_barr);
	    break;
	}
	*/
	//pthread_barrier_wait(global_barr);
	subworker.localWait();
	vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);

	if (subworker.isSubMaster()) {
	    //pthread_barrier_wait(global_barr);
	    subworker.globalWait();
	    switchFrontier(tid, Frontier, output); //set new frontier
	} else {
	    output = Frontier->getFrontier(tid);
	    //pthread_barrier_wait(global_barr);
	    subworker.globalWait();
	}

	if (subworker.isSubMaster()) {
	    Frontier->calculateNumOfNonZero(tid);	   	  	  	    
	}

	//pthread_barrier_wait(global_barr);
	subworker.globalWait();
	//swap(IDs, PrevIDs);
	//break;
    }

    if (subworker.isMaster()) {
	gettimeofday(&endT, &tz);
	double time1 = ((double)startT.tv_sec) + ((double)startT.tv_usec) / 1000000.0;
	double time2 = ((double)endT.tv_sec) + ((double)endT.tv_usec) / 1000000.0;
	printf("time used %lf\n", time2 - time1);
	cout << "Finished in " << currIter << " iterations.\n";
    }

    pthread_barrier_wait(master_barr);
    return NULL;
}

template <class vertex>
void *ComponentsWorker(void *args) {
    Default_worker_arg *my_arg = (Default_worker_arg *)args;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    int tid = my_arg->tid;
    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    graph<vertex> localGraph = graphFilter2Direction(GA, rangeLow, rangeHi);
    
    while (shouldStart == 0);
    
    const intT n = GA.n;
    int numOfT = my_arg->numOfNode;
    int blockSize = rangeHi - rangeLow;
    
    bool *frontier = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    intT outEdgesCount = 0;
    for(intT i=0;i<blockSize;i++) {
	frontier[i] = true;
	outEdgesCount += GA.V[i + rangeLow].getOutDegree();
    }

    LocalFrontier *current = new LocalFrontier(frontier, rangeLow, rangeHi);
    current->m = blockSize;
    current->outEdgesCount = outEdgesCount;

    if (tid == 0)
	Frontier = new vertices(numOfT);

    pthread_barrier_wait(&barr);
    Frontier->registerFrontier(tid, current);
    pthread_barrier_wait(&barr);

    if (tid == 0) {
	Frontier->calculateOffsets();
    }

    bool *next = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    for (intT i = 0; i < blockSize; i++) next[i] = false;
    
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);
    
    int sizeOfShards[CORES_PER_NODE];
    subPartitionByDegree(localGraph, CORES_PER_NODE, sizeOfShards, sizeof(intT), true, true);

    pthread_barrier_t masterBarr;
    pthread_barrier_init(&masterBarr, NULL, CORES_PER_NODE+1);

    pthread_barrier_t nodeBarr;
    pthread_barrier_init(&nodeBarr, NULL, CORES_PER_NODE);

    int startPos = 0;
    pthread_t subTids[CORES_PER_NODE];
    
    volatile int local_counter = 0;
    volatile int local_toggle = 0;

    for (int i = 0; i < CORES_PER_NODE; i++) {	
	Default_subworker_arg *arg = (Default_subworker_arg *)malloc(sizeof(Default_subworker_arg));
	arg->GA = (void *)(&localGraph);
	arg->Global_G = (void *)(&GA);
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;

	arg->node_barr = &nodeBarr;
	arg->master_barr = &masterBarr;
	arg->global_barr = &global_barr;
	arg->localFrontier = output;

	arg->local_custom_counter = &local_counter;
	arg->local_custom_toggle = &local_toggle;
	
	arg->startPos = startPos;
	arg->endPos = startPos + sizeOfShards[i];
	startPos = arg->endPos;
        pthread_create(&subTids[i], NULL, ComponentsSubWorker<vertex>, (void *)arg);
    }

    pthread_barrier_wait(&barr);
    if (tid == 0)
	GA.del();
    pthread_barrier_wait(&barr);
    pthread_barrier_wait(&timerBarr);

    intT *IDs = IDs_global;
    Default_Hash_F hasher(GA.n, numOfNode);
    for (intT i = rangeLow; i < rangeHi; i++) {
	IDs[hasher.hashFunc(i)] = i;
    }

    pthread_barrier_wait(&masterBarr);
    //computation of subworkers
    pthread_barrier_wait(&masterBarr);

    pthread_barrier_wait(&barr);

    return NULL;
}

template <class vertex>
void Components(graph<vertex> &GA) {
    numOfNode = numa_num_configured_nodes();
    int numOfCpu = numa_num_configured_cpus();
    CORES_PER_NODE = numOfCpu / numOfNode;
    printf("cores_per_node: %d\n", CORES_PER_NODE);
    vPerNode = GA.n / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&global_barr, NULL, numOfNode * CORES_PER_NODE);
    pthread_barrier_init(&timerBarr, NULL, numOfNode+1);
    int sizeArr[numOfNode];
    Default_Hash_F hasher(GA.n, numOfNode);
    graphAllEdgeHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(intT));
    /*
    intT vertPerPage = PAGESIZE / sizeof(double);
    intT subShardSize = ((GA.n / numOfNode) / vertPerPage) * vertPerPage;
    for (int i = 0; i < numOfNode - 1; i++) {
	sizeArr[i] = subShardSize;
    }
    sizeArr[numOfNode - 1] = GA.n - subShardSize * (numOfNode - 1);
    */
    IDs_global = (intT *)mapDataArray(numOfNode, sizeArr, sizeof(intT));
    PrevIDs_global = (intT *)mapDataArray(numOfNode, sizeArr, sizeof(intT));    

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    int prev = 0;
    for (int i = 0; i < numOfNode; i++) {
	Default_worker_arg *arg = (Default_worker_arg *)malloc(sizeof(Default_worker_arg));
	arg->GA = (void *)(&GA);
	arg->tid = i;
	arg->numOfNode = numOfNode;
	arg->rangeLow = prev;
	arg->rangeHi = prev + sizeArr[i];
	prev = prev + sizeArr[i];
	pthread_create(&tids[i], NULL, ComponentsWorker<vertex>, (void *)arg);
    }
    shouldStart = 1;
    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    startTime();
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("Components");


    if (needResult) {
	for (intT i = 0; i < GA.n; i++) {
	    printf("Result of %d : %d\n", i, IDs_global[hasher.hashFunc(i)]);
	}
    }
}

int parallel_main(int argc, char* argv[]) {  
  char* iFile;
  bool binary = false;
  bool symmetric = false;
  needResult = false;
  if(argc > 1) iFile = argv[1];
  if(argc > 2) if((string) argv[2] == (string) "-s") symmetric = true;
  if(argc > 3) if((string) argv[3] == (string) "-result") needResult = true;
  if(argc > 4) if((string) argv[4] == (string) "-b") binary = true;

  if(symmetric) {
    graph<symmetricVertex> G = 
      readGraph<symmetricVertex>(iFile,symmetric,binary);
    Components(G);
    //G.del(); 
  } else {
    graph<asymmetricVertex> G = 
      readGraph<asymmetricVertex>(iFile,symmetric,binary);
    Components(G);
    //G.del();
  }
}
