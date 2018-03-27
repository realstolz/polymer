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

double *p_curr_global = NULL;
double *p_next_global = NULL;

double *p_ans = NULL;
intT vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

pthread_barrier_t barr;
pthread_barrier_t global_barr;
pthread_mutex_t mut;

vertices *Frontier;

template <class vertex>
struct PR_F {
    double* p_curr, *p_next;
    vertex* V;
    intT rangeLow;
    intT rangeHi;
    PR_F(double* _p_curr, double* _p_next, vertex* _V, intT _rangeLow, intT _rangeHi) : 
	p_curr(_p_curr), p_next(_p_next), V(_V), rangeLow(_rangeLow), rangeHi(_rangeHi) {}

    inline void *nextPrefetchAddr(intT index) {
	return &p_curr[index];
    }
    inline bool update(intT s, intT d){ //update function applies PageRank equation
	p_next[d] += p_curr[s]/V[s].getOutDegree();
	return 1;
    }
    inline bool updateAtomic (intT s, intT d) { //atomic Update
	writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
	return 1;
    }
    inline bool cond (intT d) { return (rangeLow <= d && d < rangeHi); } //does nothing
};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
    double damping;
    double addedConstant;
    double* p_curr;
    double* p_next;
    PR_Vertex_F(double* _p_curr, double* _p_next, double _damping, intT n) :
	p_curr(_p_curr), p_next(_p_next), 
	damping(_damping), addedConstant((1-_damping)*(1/(double)n)){}
    inline bool operator () (intT i) {
	p_next[i] = damping*p_next[i] + addedConstant;
	return 1;
    }
};

//resets p
struct PR_Vertex_Reset {
    double* p_curr;
    PR_Vertex_Reset(double* _p_curr) :
	p_curr(_p_curr) {}
    inline bool operator () (intT i) {
	p_curr[i] = 0.0;
	return 1;
    }
};

struct PR_worker_arg {
    void *GA;
    int maxIter;
    int tid;
    int numOfNode;
    intT rangeLow;
    intT rangeHi;
};

struct PR_subworker_arg {
    void *GA;
    int maxIter;
    int tid;
    int subTid;
    intT startPos;
    intT endPos;
    intT rangeLow;
    intT rangeHi;
    double **p_curr_ptr;
    double **p_next_ptr;
    double damping;
    pthread_barrier_t *node_barr;
    LocalFrontier *localFrontier;
};

template <class vertex>
void *PageRankSubWorker(void *arg) {
    PR_subworker_arg *my_arg = (PR_subworker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    const intT n = GA.n;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;
    pthread_barrier_t *local_barr = my_arg->node_barr;
    LocalFrontier *output = my_arg->localFrontier;

    double *p_curr = *(my_arg->p_curr_ptr);
    double *p_next = *(my_arg->p_next_ptr);
    
    double damping = my_arg->damping;
    int currIter = 0;
    intT rangeLow = my_arg->rangeLow;
    intT rangeHi = my_arg->rangeHi;

    intT start = my_arg->startPos;
    intT end = my_arg->endPos;

    Subworker_Partitioner subworker(CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = &global_barr;

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
	    {parallel_for(long long i=output->startID;i<output->endID;i++) output->setBit(i, false);}
	}
	
	pthread_barrier_wait(&global_barr);
	//pthread_barrier_wait(local_barr);

        edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,0,DENSE_PARALLEL, false, true, subworker);

	pthread_barrier_wait(&global_barr);
	//pthread_barrier_wait(local_barr);
	if (subTid == 0) {
	    //printf("next active: %d\n", output->m);
	}

        vertexMap(Frontier, PR_Vertex_F(p_curr, p_next, damping, n), tid, subTid, CORES_PER_NODE);
	//vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);
	output->m = 1;

	pthread_barrier_wait(&global_barr);	
	//pthread_barrier_wait(local_barr);

	vertexMap(Frontier,PR_Vertex_Reset(p_curr), tid, subTid, CORES_PER_NODE);
	pthread_barrier_wait(&global_barr);
	//pthread_barrier_wait(local_barr);
	swap(p_curr, p_next);
	if (subworker.isSubMaster()) {
	    pthread_barrier_wait(&global_barr);
	    switchFrontier(tid, Frontier, output);
	} else {
	    output = Frontier->getFrontier(tid);
	    pthread_barrier_wait(&global_barr);
	}
	//pthread_barrier_wait(local_barr);
    }
    if (subworker.isMaster()) {
	p_ans = p_curr;
    }
    pthread_barrier_wait(local_barr);
    return NULL;
}

pthread_barrier_t timerBarr;

template <class vertex>
void *PageRankThread(void *arg) {
    PR_worker_arg *my_arg = (PR_worker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    intT rangeLow = my_arg->rangeLow;
    intT rangeHi = my_arg->rangeHi;

    //graph<vertex> localGraph = graphFilter(GA, rangeLow, rangeHi);

    intT sizeOfShards[CORES_PER_NODE];

    subPartitionByDegree(GA, CORES_PER_NODE, sizeOfShards, sizeof(double), rangeLow, rangeHi);
    
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
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    int numOfT = my_arg->numOfNode;

    intT blockSize = rangeHi - rangeLow;

    //printf("blockSizeof %d: %d low: %d high: %d\n", tid, blockSize, rangeLow, rangeHi);

    double one_over_n = 1/(double)n;
    
    
    double* p_curr = p_curr_global;
    double* p_next = p_next_global;
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

    for(intT i=rangeLow;i<rangeHi;i++) p_curr[i] = one_over_n;
    for(intT i=rangeLow;i<rangeHi;i++) p_next[i] = 0; //0 if unchanged
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

    intT startPos = 0;

    pthread_t subTids[CORES_PER_NODE];    

    for (int i = 0; i < CORES_PER_NODE; i++) {	
	PR_subworker_arg *arg = (PR_subworker_arg *)malloc(sizeof(PR_subworker_arg));
	arg->GA = (void *)(&GA);
	arg->maxIter = maxIter;
	arg->tid = tid;
	arg->subTid = i;
	arg->rangeLow = rangeLow;
	arg->rangeHi = rangeHi;
	arg->p_curr_ptr = &p_curr;
	arg->p_next_ptr = &p_next;
	arg->damping = damping;
	arg->node_barr = &localBarr;
	arg->localFrontier = output;
	
	arg->startPos = startPos + rangeLow;
	arg->endPos = startPos + rangeLow + sizeOfShards[i];
	startPos = startPos + sizeOfShards[i];
        pthread_create(&subTids[i], NULL, PageRankSubWorker<vertex>, (void *)arg);
    }

    pthread_barrier_wait(&barr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&barr);
    intT round = 0;
    /*
    while(1){
	if (maxIter > 0 && round >= maxIter)
	    break;
	round++;

	pthread_barrier_wait(&localBarr);

	//edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,GA.m/20,DENSE_FORWARD);
	pthread_barrier_wait(&localBarr);
	
	//vertexMap(Frontier, PR_Vertex_F(p_curr, p_next, damping, n), tid);      

	pthread_barrier_wait(&barr);
	pthread_barrier_wait(&localBarr);
	//vertexMap(Frontier,PR_Vertex_Reset(p_curr), tid);
	pthread_barrier_wait(&localBarr);

	swap(p_curr,p_next);
	if (tid == 0) {
	    p_ans = p_curr;
	}

	switchFrontier(tid, Frontier, output);

	pthread_barrier_wait(&localBarr);	
	pthread_barrier_wait(&barr);
    }
    */
    return NULL;
}

struct PR_Hash_F {
    int shardNum;
    intT vertPerShard;
    intT n;
    PR_Hash_F(intT _n, int _shardNum):n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum){}
    
    inline intT hashFunc(intT index) {
	if (index >= shardNum * vertPerShard) {
	    return index;
	}
	intT idxOfShard = index % shardNum;
	intT idxInShard = index / shardNum;
	return (idxOfShard * vertPerShard + idxInShard);
    }

    inline intT hashBackFunc(intT index) {
	if (index >= shardNum * vertPerShard) {
	    return index;
	}
	intT idxOfShard = index / vertPerShard;
	intT idxInShard = index % vertPerShard;
	return (idxOfShard + idxInShard * shardNum);
    }
};

template <class vertex>
void PageRank(graph<vertex> &GA, int maxIter) {
    numOfNode = numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
    CORES_PER_NODE = numa_num_configured_cpus() / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&timerBarr, NULL, numOfNode+1);
    pthread_barrier_init(&global_barr, NULL, CORES_PER_NODE * numOfNode);
    pthread_mutex_init(&mut, NULL);
    intT sizeArr[numOfNode];
    PR_Hash_F hasher(GA.n, numOfNode);
    graphInEdgeHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(double));
    
    p_curr_global = (double *)mapDataArray(numOfNode, sizeArr, sizeof(double));
    p_next_global = (double *)mapDataArray(numOfNode, sizeArr, sizeof(double));

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    intT prev = 0;
    for (int i = 0; i < numOfNode; i++) {
	PR_worker_arg *arg = (PR_worker_arg *)malloc(sizeof(PR_worker_arg));
	arg->GA = (void *)(&GA);
	arg->maxIter = maxIter;
	arg->tid = i;
	arg->numOfNode = numOfNode;
	arg->rangeLow = prev;
	arg->rangeHi = prev + sizeArr[i];
	prev = prev + sizeArr[i];
	pthread_create(&tids[i], NULL, PageRankThread<vertex>, (void *)arg);
    }
    shouldStart = 1;
    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    startTime();
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("PageRank");
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
	PageRank(G, maxIter);
	//G.del(); 
    } else {
	graph<asymmetricVertex> G = 
	    readGraph<asymmetricVertex>(iFile,symmetric,binary);
	PageRank(G, maxIter);
	//G.del();
    }
    return 0;
}
