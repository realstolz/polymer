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
#include <sys/syscall.h>

#include <sched.h>
#include <string>
#include <chrono>

//#include <papi.h>
#define NUM_EVENTS 3

using namespace std;

#define PAGE_SIZE (4096)

int CORES_PER_NODE = 6;
int NODE_USED = -1;

volatile int shouldStart = 0;

double *p_curr_global = NULL;
double *p_next_global = NULL;

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

template<class vertex>
struct PR_F {
    double *p_curr, *p_next;
    vertex *V;
    int rangeLow;
    int rangeHi;

    PR_F(double *_p_curr, double *_p_next, vertex *_V, int _rangeLow, int _rangeHi) :
            p_curr(_p_curr), p_next(_p_next), V(_V), rangeLow(_rangeLow), rangeHi(_rangeHi) {}

    inline void *nextPrefetchAddr(intT index) {
        return &p_curr[index];
    }

    inline bool update(intT s, intT d) { //update function applies PageRank equation
        p_next[d] += p_curr[s] / V[s].getOutDegree();
        return 1;
    }

    inline double getCurrVal(intT i) {
        return p_curr[i];
    }

    inline bool updateValVer(intT s, double val, intT d) {
        writeAdd(&p_next[d], val / V[s].getOutDegree());
        return true;
    }

    inline bool updateAtomic(intT s, intT d) { //atomic Update
        writeAdd(&p_next[d], p_curr[s] / V[s].getOutDegree());
        //return (p_curr[s] / V[s].getOutDegree()) >= 0;
        /*
        if (d == 110101) {
            cout << "Update from " << s << "\t" << std::scientific << std::setprecision(9) << p_curr[s]/V[s].getOutDegree() << " -- " << p_next[d] << "\n";
        }
        */
        return 1;
    }

    inline void initFunc(void *dataPtr, intT d) {
        *(double *) dataPtr = 0.0;
    }

    inline bool reduceFunc(void *dataPtr, intT s, bool print_info = false) {
        *(double *) dataPtr += p_curr[s] / (double) V[s].getOutDegree();
        if (print_info) {
            //cout << "reduce: " << s << " " << std::scientific << std::setprecision(9) << p_curr[s] / (double)V[s].getOutDegree() << " " << p_curr[s] << " " << V[s].getOutDegree() << *(double *)dataPtr << "\n";
        }
        return true;
    }

    inline bool combineFunc(void *dataPtr, intT d) {
        double val = *(double *) dataPtr;
        writeAdd((double *) &p_next[d], val);
        /*
        if (d == 77) {
            cout << "combine result: " << std::scientific << std::setprecision(9) << val << " " << p_next[d] << "\n";
        }
        */
        return true;
    }

    inline bool cond(intT d) { return true; } //does nothing
};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
    double damping;
    double addedConstant;
    double *p_curr;
    double *p_next;

    PR_Vertex_F(double *_p_curr, double *_p_next, double _damping, intT n) :
            p_curr(_p_curr), p_next(_p_next),
            damping(_damping), addedConstant((1 - _damping) * (1 / (double) n)) {}

    inline bool operator()(intT i) {
        p_next[i] = damping * p_next[i] + addedConstant;
        return 1;
    }
};

//resets p
struct PR_Vertex_Reset {
    double *p_curr;

    PR_Vertex_Reset(double *_p_curr) :
            p_curr(_p_curr) {}

    inline bool operator()(intT i) {
        p_curr[i] = 0.0;
        return 1;
    }
};

struct PR_worker_arg {
    void *GA;
    int maxIter;
    int tid;
    int numOfNode;
    int rangeLow;
    int rangeHi;
};

struct PR_subworker_arg {
    void *GA;
    int maxIter;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    double **p_curr_ptr;
    double **p_next_ptr;
    double damping;
    pthread_barrier_t *node_barr;
    LocalFrontier *localFrontier;
    volatile int *barr_counter;
    volatile int *toggle;
};

template<class F, class vertex>
bool *edgeMapDenseForwardOTHER(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next, bool part = false,
                               int start = 0, int end = 0) {
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
        m += G[i].getFakeDegree();
        if (currBitVector[i - currOffset]) {
            intT d = G[i].getFakeDegree();
            double val = f.getCurrVal(i);
            uintT ngh = 0;
            for (intT j = 0; j < d; j++) {
                ngh += G[i].getOutNeighbor(j);
                if (/*next->inRange(ngh) &&*/ f.cond(ngh) && f.updateValVer(i, val, ngh)) {
                    /*
                    if (!next->getBit(ngh)) {
                    m++;
                    outEdgesCount += G[ngh].getOutDegree();
                    }
                    */
                    next->setBit(ngh, true);
                }
            }
        }
    }
    return NULL;
}

template<class vertex>
void *PageRankSubWorker(void *arg) {
    PR_subworker_arg *my_arg = (PR_subworker_arg *) arg;
    graph<vertex> &GA = *(graph<vertex> *) my_arg->GA;
    const intT n = GA.n;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    int P_CORES_PER_NODE = CORES_PER_NODE / 2;
    int offset = subTid < P_CORES_PER_NODE ? 0 : (numOfNode - 1) * P_CORES_PER_NODE;
    int core = tid * P_CORES_PER_NODE + subTid + offset;
    CPU_SET(core, &cpuset);
    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set_t), &cpuset);

    cerr << "On " + to_string(sched_getcpu()) + "\n";

    pthread_barrier_t *local_barr = my_arg->node_barr;
    LocalFrontier *output = my_arg->localFrontier;

    double *p_curr = *(my_arg->p_curr_ptr);
    double *p_next = *(my_arg->p_next_ptr);

    double damping = my_arg->damping;
    int currIter = 0;
    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int start = my_arg->startPos;
    int end = my_arg->endPos;

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

    while (1) {
        if (maxIter > 0 && currIter >= maxIter)
            break;
        currIter++;
        if (subTid == 0)
            Frontier->calculateNumOfNonZero(tid);
        if (subTid == 0) {
            //{parallel_for(long i=output->startID;i<output->endID;i++) output->setBit(i, false);}
        }

        pthread_barrier_wait(&global_barr);
        //pthread_barrier_wait(local_barr);

        //edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,0,DENSE_FORWARD, false, true, subworker);
        clearLocalFrontier(output, subworker.tid, subworker.subTid, subworker.numOfSub);
        output->sparseCounter = 0;
        subworker.globalWait();
        //edgeMapDenseReduce(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,false,subworker);

        subworker.localWait();
        struct timeval startT, endT;
        struct timezone tz = {0, 0};
        gettimeofday(&startT, &tz);
        //edgeMapDenseForward(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output, true, subworker.dense_start, subworker.dense_end);
        edgeMapDenseForwardOTHER(GA, Frontier, PR_F<vertex>(p_curr, p_next, GA.V, rangeLow, rangeHi), output, true,
                                 subworker.dense_start, subworker.dense_end);
        //edgeMapDenseForwardDynamic(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output, subworker);
        subworker.localWait();
        gettimeofday(&endT, &tz);
        if (subworker.isSubMaster()) {
            double time1 = ((double) startT.tv_sec) + ((double) startT.tv_usec) / 1000000.0;
            double time2 = ((double) endT.tv_sec) + ((double) endT.tv_usec) / 1000000.0;
            double duration = time2 - time1;
            //printf("time of %d: %lf\n", subworker.tid * CORES_PER_NODE + subworker.subTid, duration);
        }

        output->isDense = true;

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

        vertexMap(Frontier, PR_Vertex_Reset(p_curr), tid, subTid, CORES_PER_NODE);
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

template<class vertex>
void *PageRankThread(void *arg) {
    PR_worker_arg *my_arg = (PR_worker_arg *) arg;
    graph<vertex> &GA = *(graph<vertex> *) my_arg->GA;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    if (tid == 0) {
        //printf("average is: %lf\n", GA.m / (float) (my_arg->numOfNode));
    }
    pthread_barrier_wait(&barr);
    intT degreeSum = 0;
    for (intT i = rangeLow; i < rangeHi; i++) {
        degreeSum += GA.V[i].getInDegree();
    }
    cerr << to_string(tid) + " : degree count: " + to_string(degreeSum) + "\n";

    //graph<vertex> localGraph = graphFilter(GA, rangeLow, rangeHi);
    graph<vertex> localGraph = graphFilter2Direction(GA, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);
    if (tid == 0)
        GA.del();
    pthread_barrier_wait(&barr);

    int sizeOfShards[CORES_PER_NODE];

    subPartitionByDegree(localGraph, CORES_PER_NODE, sizeOfShards, sizeof(double), true, true);
    //intT localDegrees = (intT *)malloc(sizeof(intT) * localGraph.n);

    for (int i = 0; i < CORES_PER_NODE; i++) {
        //printf("subPartition: %d %d: %d\n", tid, i, sizeOfShards[i]);
    }

    while (shouldStart == 0);
    pthread_barrier_wait(&timerBarr);
    cerr << "over filtering\n";
    /*
    if (0 != __cilkrts_set_param("nworkers","1")) {
	printf("set failed: %d\n", tid);
    }
    */

    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    int numOfT = my_arg->numOfNode;

    int blockSize = rangeHi - rangeLow;

    //printf("blockSizeof %d: %d low: %d high: %d\n", tid, blockSize, rangeLow, rangeHi);

    double one_over_n = 1 / (double) n;


    double *p_curr = p_curr_global;
    double *p_next = p_next_global;
    bool *frontier = (bool *) numa_alloc_local(sizeof(bool) * blockSize);

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

    for (intT i = rangeLow; i < rangeHi; i++) p_curr[i] = one_over_n;
    for (intT i = rangeLow; i < rangeHi; i++) p_next[i] = 0; //0 if unchanged
    for (intT i = 0; i < blockSize; i++) frontier[i] = true;
    if (tid == 0)
        Frontier = new vertices(numOfT);

    //printf("register %d: %p\n", tid, frontier);

    LocalFrontier *current = new LocalFrontier(frontier, rangeLow, rangeHi);

    bool *next = (bool *) numa_alloc_local(sizeof(bool) * blockSize);
    for (intT i = 0; i < blockSize; i++) next[i] = false;
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);

    Frontier->registerFrontier(tid, current);

    pthread_barrier_wait(&barr);

    if (tid == 0)
        Frontier->calculateOffsets();

    pthread_barrier_t localBarr;
    pthread_barrier_init(&localBarr, NULL, CORES_PER_NODE + 1);

    int startPos = 0;

    pthread_t subTids[CORES_PER_NODE];

    volatile int local_custom_counter;
    volatile int local_toggle;

    for (int i = 0; i < CORES_PER_NODE; i++) {
        PR_subworker_arg *arg = (PR_subworker_arg *) malloc(sizeof(PR_subworker_arg));
        arg->GA = (void *) (&localGraph);
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

        arg->barr_counter = &local_custom_counter;
        arg->toggle = &local_toggle;

        arg->startPos = startPos;
        arg->endPos = startPos + sizeOfShards[i];
        startPos = arg->endPos;

        pthread_create(&subTids[i], NULL, PageRankSubWorker<vertex>, (void *) arg);
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
    int vertPerShard;
    int n;

    PR_Hash_F(int _n, int _shardNum) : n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum) {}

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
void PageRank(graph<vertex> &GA, int maxIter) {
//    numOfNode = numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
//    CORES_PER_NODE = numa_num_configured_cpus() / numOfNode;
    if (NODE_USED != -1)
        numOfNode = NODE_USED;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_barrier_init(&timerBarr, NULL, numOfNode + 1);
    pthread_barrier_init(&global_barr, NULL, CORES_PER_NODE * numOfNode);
    pthread_mutex_init(&mut, NULL);
    int sizeArr[numOfNode];
    PR_Hash_F hasher(GA.n, numOfNode);
    //graphHasher(GA, hasher);
    //graphAllEdgeHasher(GA, hasher);
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(double));
    /*
    intT vertPerPage = PAGESIZE / sizeof(double);
    intT subShardSize = ((GA.n / numOfNode) / vertPerPage) * vertPerPage;
    for (int i = 0; i < numOfNode - 1; i++) {
	sizeArr[i] = subShardSize;
    }
    sizeArr[numOfNode - 1] = GA.n - subShardSize * (numOfNode - 1);
    */
    int accum = 0;
    for (int i = 0; i < numOfNode; i++) {
        intT degreeSum = 0;
        for (intT j = accum; j < accum + sizeArr[i]; j++) {
            degreeSum += GA.V[j].getInDegree();
        }
        cerr << to_string(i) + ": degree sum: " + to_string(degreeSum) + "\n";
        accum += sizeArr[i];
    }
    //return;

    p_curr_global = (double *) mapDataArray(numOfNode, sizeArr, sizeof(double));
    p_next_global = (double *) mapDataArray(numOfNode, sizeArr, sizeof(double));

    cerr << "start create " + to_string(numOfNode) << "threads\n";
    pthread_t tids[numOfNode];
    int prev = 0;
    for (int i = 0; i < numOfNode; i++) {
        PR_worker_arg *arg = (PR_worker_arg *) malloc(sizeof(PR_worker_arg));
        arg->GA = (void *) (&GA);
        arg->maxIter = maxIter;
        arg->tid = i;
        arg->numOfNode = numOfNode;
        arg->rangeLow = prev;
        arg->rangeHi = prev + sizeArr[i];
        prev = prev + sizeArr[i];

        pthread_create(&tids[i], NULL, PageRankThread<vertex>, (void *) arg);
    }
    shouldStart = 1;

    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    nextTime("partition over");
    cerr << "all created\n";
    for (int i = 0; i < numOfNode; i++) {
        pthread_join(tids[i], NULL);
    }
    nextTime("PageRank");

    if (needResult) {
        for (intT i = 0; i < GA.n; i++) {
            cerr << i << "\t" << std::scientific << std::setprecision(9) << p_ans[hasher.hashFunc(i)] << "\n";
            //cout << i << "\t" << std::scientific << std::setprecision(9) << p_ans[i] << "\n";
        }
    }
}

int parallel_main(int argc, char *argv[]) {
    char *iFile;
    bool binary = false;
    bool symmetric = false;
    int maxIter = 20;
    needResult = false;

    iFile = argv[1];
    maxIter = atoi(argv[2]);
    numOfNode = atoi(argv[3]);
    CORES_PER_NODE = atoi(argv[4]);

//    if (argc > 3) NODE_USED = atoi(argv[3]);
//    if (argc > 4) if ((string) argv[4] == (string) "-result") needResult = true;
//    if (argc > 5) if ((string) argv[5] == (string) "-s") symmetric = true;
//    if (argc > 6) if ((string) argv[6] == (string) "-b") binary = true;
    numa_set_interleave_mask(numa_all_nodes_ptr);
    startTime();
    if (symmetric) {
        graph<symmetricVertex> G =
                readGraph<symmetricVertex>(iFile, symmetric, binary);
        PageRank(G, maxIter);
        //G.del();
    } else {
        const auto start = chrono::system_clock::now();
        graph<asymmetricVertex> G =
                readGraph<asymmetricVertex>(iFile, symmetric, binary);
        const auto mid = chrono::system_clock::now();

        PageRank(G, maxIter);
        //G.del();
        const auto end = chrono::system_clock::now();

        chrono::duration<double> prepro = mid - start;
        chrono::duration<double> pr = end - mid;
        cout << "preprocessing: " << prepro.count() << endl;
        cout << "PR: " << pr.count() << endl;
    }
    return 0;
}
