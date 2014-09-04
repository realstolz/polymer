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
#include "parallel.h"
#include "custom-barrier.h"

#include <numa.h>
#include <sys/mman.h>
#include <pthread.h>

using namespace std;

int CORES_PER_NODE = 6;

volatile int shouldStart = 0;

vertices *Frontier;

intT *parents_global;

pthread_barrier_t barr;
pthread_barrier_t subMasterBarr;
pthread_barrier_t global_barr;
pthread_barrier_t timerBarr;

volatile int global_counter;
volatile int global_toggle;

int vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

struct BFS_F {
    intT* Parents;
    BFS_F(intT* _Parents) : Parents(_Parents) {}

    inline void *nextPrefetchAddr(intT index) {
	return NULL;
    }
    inline bool update (intT s, intT d) { //Update
	if(Parents[d] == -1) { Parents[d] = s; return 1; }
	else return 0;
    }
    inline bool updateAtomic (intT s, intT d){ //atomic version of Update
	return (CAS(&Parents[d],(intT)-1,s));
    }
    
    inline void vertUpdate(intT v) {
	return;
    }
    //cond function checks if vertex has been visited yet
    inline bool cond (intT d) { return (Parents[d] == -1); } 
};

struct BFS_Vert_F {
    inline bool operator () (intT i) {
	return 1;
    }
};

struct BFS_worker_arg {
    void *GA;
    int tid;
    int numOfNode;
    int rangeLow;
    int rangeHi;
    int start;
};

struct BFS_subworker_arg {
    void *GA;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    intT *parents_ptr;
    pthread_barrier_t *global_barr;
    pthread_barrier_t *node_barr;
    pthread_barrier_t *node_barr2;
    LocalFrontier *localFrontier;
    volatile int *barr_counter;
    volatile int *toggle;
};

template <class vertex>
void *BFSSubWorker(void *arg) {
    BFS_subworker_arg *my_arg = (BFS_subworker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    const intT n = GA.n;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;
    pthread_barrier_t *local_barr = my_arg->node_barr;
    pthread_barrier_t *local_barr2 = my_arg->node_barr2;
    pthread_barrier_t *global_barr = my_arg->global_barr;
    LocalFrontier *output = my_arg->localFrontier;

    intT *parents = my_arg->parents_ptr;
    
    int currIter = 0;
    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int start = my_arg->startPos;
    int end = my_arg->endPos;

    intT numVisited = 0;

    Custom_barrier localCustom(my_arg->barr_counter, my_arg->toggle, CORES_PER_NODE);
    Custom_barrier globalCustom(&global_counter, &global_toggle, Frontier->numOfNodes);

    Subworker_Partitioner subworker(CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = global_barr;
    subworker.local_barr = my_arg->node_barr2;
    subworker.leader_barr = &subMasterBarr;
    subworker.local_custom = localCustom;
    subworker.subMaster_custom = globalCustom;
    

    if (subworker.isMaster()) {
	pthread_barrier_init(&subMasterBarr, NULL, Frontier->numOfNodes);
    }

    pthread_barrier_wait(local_barr);

    if (subTid == 0) 
	Frontier->calculateNumOfNonZero(tid);

    pthread_barrier_wait(global_barr);

    struct timeval startT, endT;
    struct timezone tz = {0, 0};
    while(!Frontier->isEmpty() || currIter == 0){ //loop until frontier is empty
	currIter++;
	if (tid + subTid == 0) {
	    numVisited += Frontier->numNonzeros();
	    //printf("num of non zeros: %d\n", Frontier->numNonzeros());
	}

	//pthread_barrier_wait(global_barr);
	//apply edgemap
	gettimeofday(&startT, &tz);
	edgeMap(GA, Frontier, BFS_F(parents), output, GA.m/20, DENSE_FORWARD, false, true, subworker);
	subworker.localWait();
	vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);
	//edgeMapSparseAsync(GA, Frontier, BFS_F(parents), output, subworker);
	if (subTid == 0) {
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
	gettimeofday(&endT, &tz);
	double timeStart = ((double)startT.tv_sec) + ((double)startT.tv_usec) / 1000000.0;
	double timeEnd = ((double)endT.tv_sec) + ((double)endT.tv_usec) / 1000000.0;

	double mapTime = timeEnd - timeStart;
	if (tid + subTid == 0) {
	    printf("edge map time: %lf\n", mapTime);
	}
	//break;
    }

    if (tid + subTid == 0) {
	cout << "Vertices visited = " << numVisited << "\n";
	cout << "Finished in " << currIter << " iterations\n";
    }

    pthread_barrier_wait(local_barr);
    return NULL;
}

template <class vertex>
void *BFSWorker(void *arg) {
    BFS_worker_arg *my_arg = (BFS_worker_arg *)arg;
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
    const intT n = GA.n;
    int numOfT = my_arg->numOfNode;
    int blockSize = rangeHi - rangeLow;

    intT *parents = parents_global;

    for (intT i = rangeLow; i < rangeHi; i++) {
	parents[i] = -1;
    }
    
    bool *frontier = (bool *)numa_alloc_local(sizeof(bool) * blockSize);

    for(intT i=0;i<blockSize;i++) frontier[i] = false;

    LocalFrontier *current = new LocalFrontier(frontier, rangeLow, rangeHi);

    if (tid == 0)
	Frontier = new vertices(numOfT);

    pthread_barrier_wait(&barr);
    Frontier->registerFrontier(tid, current);
    pthread_barrier_wait(&barr);

    if (tid == 0) {
	Frontier->calculateOffsets();
	Frontier->setBit(my_arg->start, true);
	parents[my_arg->start] = my_arg->start;
    }

    if (my_arg->start >= rangeLow && my_arg->start < rangeHi) {
	current->m = 1;
	current->outEdgesCount = GA.V[my_arg->start].getOutDegree();
    }
    
    bool *next = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    for (intT i = 0; i < blockSize; i++) next[i] = false;
    
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);
    
    int sizeOfShards[CORES_PER_NODE];
    partitionByDegree(GA, CORES_PER_NODE, sizeOfShards, sizeof(intT), true);

    int startPos = 0;

    pthread_barrier_t localBarr;
    pthread_barrier_init(&localBarr, NULL, CORES_PER_NODE+1);

    pthread_barrier_t localBarr2;
    pthread_barrier_init(&localBarr2, NULL, CORES_PER_NODE);

    pthread_barrier_wait(&barr);
    /*
    if (tid == 0)
	Frontier->toSparse();
    */
    /*
    current->s = (intT *)malloc(sizeof(intT) * GA.n);
    current->s[0] = my_arg->start;
    current->m = 1;
    current->head = 0;
    current->tail = 1;
    */

    if (tid == 0) {
	Frontier->asyncQueue = (AsyncChunk **)malloc(sizeof(AsyncChunk *) * GA.n);
	{parallel_for(intT i = 0; i < GA.n; i++) Frontier->asyncQueue[i] = NULL;}
	AsyncChunk *firstChunk = newChunk(64);
	firstChunk->m = 1;
	firstChunk->s[0] = my_arg->start;
	Frontier->asyncQueue[0] = firstChunk;
	Frontier->readerTail = 1;
	Frontier->m = 1;
    }

    volatile int local_counter = 0;
    volatile int local_toggle = 0;

    pthread_t subTids[CORES_PER_NODE];
    for (int i = 0; i < CORES_PER_NODE; i++) {	
	BFS_subworker_arg *arg = (BFS_subworker_arg *)malloc(sizeof(BFS_subworker_arg));
	arg->GA = (void *)(&localGraph);
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;
	arg->parents_ptr = parents;

	arg->barr_counter = &local_counter;
	arg->toggle = &local_toggle;

	arg->node_barr = &localBarr;
	arg->node_barr2 = &localBarr2;
	arg->global_barr = &global_barr;
	arg->localFrontier = output;
	
	arg->startPos = startPos;
	arg->endPos = startPos + sizeOfShards[i];
	startPos = arg->endPos;
        pthread_create(&subTids[i], NULL, BFSSubWorker<vertex>, (void *)arg);
    }

    pthread_barrier_wait(&timerBarr);

    struct timeval timer_start, timer_end;
    struct timezone timer_zone = {0, 0};

    gettimeofday(&timer_start, &timer_zone);
    pthread_barrier_wait(&localBarr);
    //computation of subworkers
    pthread_barrier_wait(&localBarr);
    gettimeofday(&timer_end, &timer_zone);
    double timeStart = ((double)timer_start.tv_sec) + ((double)timer_start.tv_usec) / 1000000.0;
    double timeEnd = ((double)timer_end.tv_sec) + ((double)timer_end.tv_usec) / 1000000.0;
    
    double mapTime = timeEnd - timeStart;
    if (tid == 0) 
	printf("self counted time: %lf\n", mapTime);
    
    pthread_barrier_wait(&barr);
    return NULL;
}

struct PR_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;
    PR_Hash_F(int _n, int _shardNum):n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum){}
    
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
void BFS(intT start, graph<vertex> &GA) {
    numOfNode = numa_num_configured_nodes();
    int numOfCpu = numa_num_configured_cpus();
    CORES_PER_NODE = numOfCpu / numOfNode;
    vPerNode = GA.n / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&global_barr, NULL, numOfNode * CORES_PER_NODE);
    pthread_barrier_init(&timerBarr, NULL, numOfNode+1);
    int sizeArr[numOfNode];
    PR_Hash_F hasher(GA.n, numOfNode);
    graphHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(intT));
    /*
    intT vertPerPage = PAGESIZE / sizeof(double);
    intT subShardSize = ((GA.n / numOfNode) / vertPerPage) * vertPerPage;
    for (int i = 0; i < numOfNode - 1; i++) {
	sizeArr[i] = subShardSize;
    }
    sizeArr[numOfNode - 1] = GA.n - subShardSize * (numOfNode - 1);
    */
    parents_global = (intT *)mapDataArray(numOfNode, sizeArr, sizeof(intT));

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    int prev = 0;
    for (int i = 0; i < numOfNode; i++) {
	BFS_worker_arg *arg = (BFS_worker_arg *)malloc(sizeof(BFS_worker_arg));
	arg->GA = (void *)(&GA);
	arg->tid = i;
	arg->numOfNode = numOfNode;
	arg->rangeLow = prev;
	arg->rangeHi = prev + sizeArr[i];
	arg->start = hasher.hashFunc(start);	
	prev = prev + sizeArr[i];
	pthread_create(&tids[i], NULL, BFSWorker<vertex>, (void *)arg);
    }
    shouldStart = 1;
    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    startTime();
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("BFS");
    if (needResult) {
	int counter = 0;
	for (intT i = 0; i < GA.n; i++) {
	    if (parents_global[i] != -1)
		counter++;
	}
	printf("Vert visited: %d\n", counter);
    }
}

int parallel_main(int argc, char* argv[]) {  
    char* iFile;
    bool binary = false;
    bool symmetric = false;
    int start = 0;
    global_counter = 0;
    global_toggle = 0;
    if(argc > 1) iFile = argv[1];
    if(argc > 2) start = atoi(argv[2]);
    if(argc > 3) if((string) argv[3] == (string) "-result") needResult = true;
    //pass -s flag if graph is already symmetric
    if(argc > 4) if((string) argv[4] == (string) "-s") symmetric = true;
    //pass -b flag if using binary file (also need to pass 2nd arg for now)
    if(argc > 5) if((string) argv[5] == (string) "-b") binary = true;

    if(symmetric) {
	graph<symmetricVertex> G = 
	    readGraph<symmetricVertex>(iFile,symmetric,binary); //symmetric graph
	BFS((intT)start,G);
	G.del(); 
    } else {
	graph<asymmetricVertex> G = 
	    readGraph<asymmetricVertex>(iFile,symmetric,binary); //asymmetric graph
	BFS((intT)start,G);
	G.del();
    }
}
