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
#include "ligra-rewrite-wgh.h"
#include "gettime.h"
#include "math.h"

#include <pthread.h>
#include <sys/mman.h>
#include <numa.h>
using namespace std;

#define PAGE_SIZE (4096)

int CORES_PER_NODE = 6;

volatile int shouldStart = 0;

int *ShortestPathLen_global = NULL;
int *Visited_global = NULL;

double *p_ans = NULL;
int vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

pthread_barrier_t barr;
pthread_barrier_t global_barr;
pthread_mutex_t mut;

volatile int global_counter = 0;
volatile int global_toggle = 0;

vertices *Frontier;

struct BF_F {
    int* ShortestPathLen;
    int* Visited;
    BF_F(int* _ShortestPathLen, int* _Visited) : 
	ShortestPathLen(_ShortestPathLen), Visited(_Visited) {}

    inline void *nextPrefetchAddr(intT index) {
	return &ShortestPathLen[index];
    }

    inline bool update (intT s, intT d, intT edgeLen) { //Update ShortestPathLen if found a shorter path
	intT newDist = ShortestPathLen[s] + edgeLen;
	if(ShortestPathLen[d] > newDist) {
	    ShortestPathLen[d] = newDist;
	    if(Visited[d] == 0) { Visited[d] = 1 ; return 1;}
	}
	return 0;
    }
    inline bool updateAtomic (intT s, intT d, int edgeLen){ //atomic Update
	int newDist = ShortestPathLen[s] + edgeLen;
	return (writeMin(&ShortestPathLen[d],newDist) &&
		CAS(&Visited[d],0,1));
    }
    inline bool cond (intT d) { return 1; } //does nothing
};

//reset visited vertices
struct BF_Vertex_F {
    int* Visited;
    BF_Vertex_F(int* _Visited) : Visited(_Visited) {}
    inline bool operator() (intT i){
	Visited[i] = 0;
	return 1;
    }
};

struct BF_worker_arg {
    void *GA;
    int tid;
    int numOfNode;
    intT start;
    int rangeLow;
    int rangeHi;
};

struct BF_subworker_arg {
    void *GA;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    int **ShortestPathLen_ptr;
    int **Visited_ptr;
    pthread_barrier_t *node_barr;
    pthread_barrier_t *node_barr2;
    LocalFrontier *localFrontier;

    volatile int *local_custom_counter;
    volatile int *local_custom_toggle;
};

template <class vertex>
void *BFSubWorker(void *arg) {
    BF_subworker_arg *my_arg = (BF_subworker_arg *)arg;
    wghGraph<vertex> &GA = *(wghGraph<vertex> *)my_arg->GA;
    const intT n = GA.n;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;
    pthread_barrier_t *local_barr = my_arg->node_barr;
    LocalFrontier *output = my_arg->localFrontier;

    int *ShortestPathLen = *(my_arg->ShortestPathLen_ptr);
    int *Visited = *(my_arg->Visited_ptr);
    
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
    subworker.global_barr = &global_barr;
    subworker.local_barr = my_arg->node_barr2;
    subworker.local_custom = local_custom;
    subworker.subMaster_custom = global_custom;

    pthread_barrier_wait(local_barr);
    if (subworker.isMaster())
	printf("started\n");
    if (subworker.isSubMaster()) {
	Frontier->calculateNumOfNonZero(tid);
    }
    pthread_barrier_wait(&global_barr);
    while(!Frontier->isEmpty() || currIter == 0){ //loop until frontier is empty
	currIter++;
	if (tid + subTid == 0) {
	    numVisited += Frontier->numNonzeros();
	    printf("Round %d: num of non zeros: %d\n", currIter, Frontier->numNonzeros());
	}

	if (subTid == 0) {
	    //{parallel_for(long i=output->startID;i<output->endID;i++) output->setBit(i, false);}
	}
	//pthread_barrier_wait(&global_barr);
	subworker.globalWait();
	//apply edgemap
	edgeMap(GA, Frontier, BF_F(ShortestPathLen, Visited), output, GA.n/10, DENSE_FORWARD, false, true, subworker);
	//pthread_barrier_wait(&global_barr);
	subworker.globalWait();
        vertexMap(Frontier, BF_Vertex_F(Visited), tid, subTid, CORES_PER_NODE);
	vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);
	//edgeMapSparseAsync(GA, Frontier, BF_F(parents), output, subworker);
	if (subTid == 0) {
	    //pthread_barrier_wait(&global_barr);
	    subworker.globalWait();
	    switchFrontier(tid, Frontier, output); //set new frontier
	} else {
	    output = Frontier->getFrontier(tid);
	    //pthread_barrier_wait(&global_barr);
	    subworker.globalWait();
	}

	if (subworker.isSubMaster()) {
	    Frontier->calculateNumOfNonZero(tid);	   	  	  	    
	}
	//pthread_barrier_wait(&global_barr);
	subworker.globalWait();
    }

    if (tid + subTid == 0) {
	cout << "Vertices visited = " << numVisited << "\n";
	cout << "Finished in " << currIter << " iterations\n";
    }
    pthread_barrier_wait(local_barr);
    return NULL;
}

pthread_barrier_t timerBarr;

template <class vertex>
void *BFThread(void *arg) {
    BF_worker_arg *my_arg = (BF_worker_arg *)arg;
    wghGraph<vertex> &GA = *(wghGraph<vertex> *)my_arg->GA;
    int tid = my_arg->tid;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    wghGraph<vertex> localGraph = graphFilter(GA, rangeLow, rangeHi);

    int sizeOfShards[CORES_PER_NODE];

    subPartitionByDegree(localGraph, CORES_PER_NODE, sizeOfShards, sizeof(int), true, true);
    
    for (int i = 0; i < CORES_PER_NODE; i++) {
	//printf("subPartition: %d %d: %d\n", tid, i, sizeOfShards[i]);
    }

    while (shouldStart == 0) ;
    pthread_barrier_wait(&timerBarr);
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

    double one_over_n = 1/(double)n;
    
    
    int *ShortestPathLen = ShortestPathLen_global;
    int *Visited = Visited_global;
    bool* frontier = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    
    /*
      double* p_curr = (double *)malloc(sizeof(double) * blockSize);
      double* p_next = (double *)malloc(sizeof(double) * blockSize);
      bool* frontier = (bool *)malloc(sizeof(bool) * blockSize);
    */
    
    /*
      if (tid == 0)
      startTime();
    */
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    for(intT i=rangeLow;i<rangeHi;i++) ShortestPathLen[i] = INT_MAX/2;
    for(intT i=rangeLow;i<rangeHi;i++) Visited[i] = 0;
    for(intT i=0;i<blockSize;i++) frontier[i] = false;
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

    if (tid == 0) {
	Frontier->calculateOffsets();
	Frontier->setBit(my_arg->start, true);
    }

    if (my_arg->start >= rangeLow && my_arg->start < rangeHi) {
	current->m = 1;
	current->outEdgesCount = GA.V[my_arg->start].getOutDegree();
	ShortestPathLen[my_arg->start] = 0;
    } else {
	current->m = 0;
	current->outEdgesCount = 0;
    }

    pthread_barrier_t localBarr;
    pthread_barrier_init(&localBarr, NULL, CORES_PER_NODE+1);

    pthread_barrier_t localBarr2;
    pthread_barrier_init(&localBarr2, NULL, CORES_PER_NODE);

    int startPos = 0;

    pthread_t subTids[CORES_PER_NODE];    
    
    volatile int local_counter = 0;
    volatile int local_toggle = 0;

    for (int i = 0; i < CORES_PER_NODE; i++) {	
	BF_subworker_arg *arg = (BF_subworker_arg *)malloc(sizeof(BF_subworker_arg));
	arg->GA = (void *)(&localGraph);
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;
	arg->ShortestPathLen_ptr = &ShortestPathLen;
	arg->Visited_ptr = &Visited;
	arg->node_barr = &localBarr;
	arg->node_barr2 = &localBarr2;
	arg->localFrontier = output;

	arg->local_custom_counter = &local_counter;
	arg->local_custom_toggle = &local_toggle;
	
	arg->startPos = startPos;
	arg->endPos = startPos + sizeOfShards[i];
	startPos = arg->endPos;
        pthread_create(&subTids[i], NULL, BFSubWorker<vertex>, (void *)arg);
    }

    pthread_barrier_wait(&barr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&barr);
    return NULL;
}

struct BF_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;
    BF_Hash_F(int _n, int _shardNum):n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum){}
    
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
void BF_main(wghGraph<vertex> &GA, intT start) {
    numOfNode = 1;//numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
    CORES_PER_NODE = 10;//numa_num_configured_cpus() / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&timerBarr, NULL, numOfNode+1);
    pthread_barrier_init(&global_barr, NULL, CORES_PER_NODE * numOfNode);
    pthread_mutex_init(&mut, NULL);
    int sizeArr[numOfNode];
    BF_Hash_F hasher(GA.n, numOfNode);
    graphHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(int));
    
    ShortestPathLen_global = (int *)mapDataArray(numOfNode, sizeArr, sizeof(int));
    Visited_global = (int *)mapDataArray(numOfNode, sizeArr, sizeof(int));

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    int prev = 0;
    for (int i = 0; i < numOfNode; i++) {
	BF_worker_arg *arg = (BF_worker_arg *)malloc(sizeof(BF_worker_arg));
	arg->GA = (void *)(&GA);
	arg->tid = i;
	arg->numOfNode = numOfNode;
	arg->start = hasher.hashFunc(start);
	arg->rangeLow = prev;
	arg->rangeHi = prev + sizeArr[i];
	prev = prev + sizeArr[i];
	pthread_create(&tids[i], NULL, BFThread<vertex>, (void *)arg);
    }
    shouldStart = 1;
    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    startTime();
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("BellmanFord");
    if (needResult) {
	for (intT i = 0; i < GA.n; i++) {
	    cout << i << "\t" << std::scientific << std::setprecision(9) << p_ans[hasher.hashFunc(i)] << "\n";
	}
    }
}

int parallel_main(int argc, char* argv[]) {  
    char* iFile;
    bool binary = false;
    bool symmetric = false;
    int startPos = 1;
    needResult = false;
    if(argc > 1) iFile = argv[1];
    if(argc > 2) startPos = atoi(argv[2]);
    if(argc > 3) if((string) argv[3] == (string) "-result") needResult = true;
    if(argc > 4) if((string) argv[4] == (string) "-s") symmetric = true;
    if(argc > 5) if((string) argv[5] == (string) "-b") binary = true;
    numa_set_interleave_mask(numa_all_nodes_ptr);
    if(symmetric) {
	wghGraph<symmetricWghVertex> WG = 
	    readWghGraph<symmetricWghVertex>(iFile,symmetric,binary);
	BF_main(WG, (intT)startPos);
	WG.del(); 
    } else {
	wghGraph<asymmetricWghVertex> WG = 
	    readWghGraph<asymmetricWghVertex>(iFile,symmetric,binary);
	BF_main(WG, (intT)startPos);
	WG.del();
    }
    return 0;
}
