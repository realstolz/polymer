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
#include "ligra-rewrite.h"
#include "gettime.h"

#include <numa.h>
#include <pthread.h>
using namespace std;

pthread_barrier_t barr;
pthread_barrier_t global_barr;
pthread_barrier_t timerBarr;

volatile int shouldStart = 0;

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
	//printf("update of edge (%d)%d -> %d(%d / %d to %d), %s\n", IDs[s], s, d, origID, prevIDs[d], IDs[d], res ? "YES" : "NO");
	return res;
	/*
	  if (origID != prevIDs[d])
	  printf("diff ID: %d, %d, %d\n", d, origID, prevIDs[d]);
	*/
	/*
	  return (writeMin(&IDs[d],IDs[s]) 
	  && origID == prevIDs[d]);
	*/
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

template <class vertex>
void *ComponentsSubWorker(void *args) {
    Default_subworker_arg *my_arg = (Default_subworker_arg *)args;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
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

    Subworker_Partitioner subworker(CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = global_barr;
    subworker.local_barr = my_arg->node_barr;

    intT *IDs = IDs_global;
    intT *PrevIDs = PrevIDs_global;

    pthread_barrier_wait(master_barr);

    if (subworker.isSubMaster()) {
	Frontier->calculateNumOfNonZero(tid);
    }

    pthread_barrier_wait(global_barr);
    
    intT switchThreshold = GA.m/20;

    while (!Frontier->isEmpty() || currIter == 0) {
	currIter++;
	if (subworker.isMaster()) {
	    numVisited += Frontier->numNonzeros();
	    printf("non zeros: %d\n", Frontier->numNonzeros());
	}
	
	//clearLocalFrontier(output, tid, subTid, CORES_PER_NODE);

	vertexMap(Frontier, CC_Vertex_F(IDs,PrevIDs), tid, subTid, CORES_PER_NODE);
	pthread_barrier_wait(global_barr);
	edgeMap(GA, Frontier, CC_F(IDs,PrevIDs), output, switchThreshold, DENSE_FORWARD, false, true, subworker);
	pthread_barrier_wait(global_barr);

	vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);

	if (subworker.isSubMaster()) {
	    pthread_barrier_wait(global_barr);
	    switchFrontier(tid, Frontier, output); //set new frontier
	} else {
	    output = Frontier->getFrontier(tid);
	    pthread_barrier_wait(global_barr);
	}

	if (subworker.isSubMaster()) {
	    Frontier->calculateNumOfNonZero(tid);	   	  	  	    
	}

	pthread_barrier_wait(global_barr);
	//swap(IDs, PrevIDs);
	//break;
    }

    if (subworker.isMaster()) {
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

    graph<vertex> localGraph = graphFilter(GA, rangeLow, rangeHi);
    
    while (shouldStart == 0);
    pthread_barrier_wait(&timerBarr);
    
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
    for (int i = 0; i < CORES_PER_NODE; i++) {	
	Default_subworker_arg *arg = (Default_subworker_arg *)malloc(sizeof(Default_subworker_arg));
	arg->GA = (void *)(&localGraph);
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;

	arg->node_barr = &nodeBarr;
	arg->master_barr = &masterBarr;
	arg->global_barr = &global_barr;
	arg->localFrontier = output;
	
	arg->startPos = startPos;
	arg->endPos = startPos + sizeOfShards[i];
	startPos = arg->endPos;
        pthread_create(&subTids[i], NULL, ComponentsSubWorker<vertex>, (void *)arg);
    }

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
    graphHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(intT));
    
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
    G.del(); 
  } else {
    graph<asymmetricVertex> G = 
      readGraph<asymmetricVertex>(iFile,symmetric,binary);
    Components(G);
    G.del();
  }
}
