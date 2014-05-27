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
#include "math.h"

#include <pthread.h>
#include <sys/mman.h>
#include <numa.h>
using namespace std;

#define PAGE_SIZE (4096)

volatile int shouldStart = 0;

double *p_curr_global = NULL;
double *p_next_global = NULL;

double *p_ans = NULL;
int vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

template <class vertex>
struct PR_F {
    double* p_curr, *p_next;
    vertex* V;
    int rangeLow;
    int rangeHi;
    PR_F(double* _p_curr, double* _p_next, vertex* _V, int _rangeLow, int _rangeHi) : 
	p_curr(_p_curr), p_next(_p_next), V(_V), rangeLow(_rangeLow), rangeHi(_rangeHi) {}
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
    int rangeLow;
    int rangeHi;
};

pthread_barrier_t barr;
pthread_mutex_t mut;

vertices *Frontier;

template <class vertex>
void *PageRankThread(void *arg) {
    PR_worker_arg *my_arg = (PR_worker_arg *)arg;
    graph<vertex> &GA = *(graph<vertex> *)my_arg->GA;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;

    while (shouldStart == 0) ;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);
    /*
    if (0 != __cilkrts_set_param("nworkers","1")) {
	printf("set failed: %d\n", tid);
    }
    */

    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    int numOfT = my_arg->numOfNode;

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int blockSize = rangeHi - rangeLow;

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

    pthread_barrier_wait(&barr);
    
    Frontier->registerArr(tid, frontier, blockSize);

    pthread_barrier_wait(&barr);

    if (tid == 0)
	Frontier->calculateOffsets();

    bool *next = (bool *)numa_alloc_local(sizeof(bool) * blockSize);
    for (intT i = 0; i < blockSize; i++) next[i] = false;
    
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);

    intT round = 0;
    while(1){
	if (maxIter > 0 && round >= maxIter)
	    break;
	round++;

	//gettimeofday(&start, &tz);
	edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,GA.m/20,DENSE_FORWARD);

	//gettimeofday(&end, &tz);
	//double timeStart = ((double)start.tv_sec) + ((double)start.tv_usec) / 1000000.0;
	//double timeEnd = ((double)end.tv_sec) + ((double)end.tv_usec) / 1000000.0;
	//mapTime = timeEnd - timeStart;

	//output.del();

	vertexMap(Frontier, PR_Vertex_F(p_curr, p_next, damping, n), tid);

	/*
	double tmpConst = (1-damping)*(1/(double)n);

	for (intT i = rangeLow; i < rangeHi; i++) {
		p_next[i] = damping * p_next[i] + tmpConst;
	}
	*/
	//vertexMap(Frontier,PR_Vertex_F(p_curr,p_next,damping,n));
	/*
	//compute L1-norm between p_curr and p_next
	{parallel_for(intT i=0;i<n;i++) {
		p_curr[i] = fabs(p_curr[i]-p_next[i]);
	    }}
	double L1_norm = sequence::plusReduce(p_curr,n);
	//cout<<"Round "<<round<<", L1 norm = "<<L1_norm<<endl;
	if(L1_norm < epsilon) break;
	//reset p_curr
	*/

	/*
	for (intT i = rangeLow; i < rangeHi; i++) {
		p_curr[i] = p_next[i];
		p_next[i] = 0.0;
	}
	*/

	pthread_barrier_wait(&barr);
	vertexMap(Frontier,PR_Vertex_Reset(p_curr), tid);
	swap(p_curr,p_next);
	if (tid == 0) {
	    p_ans = p_curr;
	}
	
	//Frontier.del(); 
	//Frontier = output;
	switchFrontier(tid, Frontier, output);
	pthread_barrier_wait(&barr);

	/*
	  if (tid == 0)
	    printf("round %d: %lf\n", round, mapTime);
	*/
    }
    //cout<<"Finished in "<<round<<" iterations\n";
    //Frontier.del()
    
    /*
    if (tid == 0)
	nextTime("PageRank");
    */
    //printf("total init time: %lf\n", mapTime);
    //pthread_exit(NULL);    
    return NULL;
}

template <class vertex>
void PageRank(graph<vertex> GA, int maxIter) {
    numOfNode = numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_mutex_init(&mut, NULL);
    int sizeArr[numOfNode];
    partitionByDegree(GA, numOfNode, sizeArr, sizeof(double));
    
    p_curr_global = (double *)mapDataArray(numOfNode, sizeArr, sizeof(double));
    p_next_global = (double *)mapDataArray(numOfNode, sizeArr, sizeof(double));

    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    int prev = 0;
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
    startTime();
    shouldStart = 1;
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("PageRank");
    if (needResult) {
	for (intT i = 0; i < GA.n; i++) {
	    cout << i << "\t" << std::scientific << std::setprecision(9) << p_ans[i] << "\n";
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
    
    if(symmetric) {
	graph<symmetricVertex> G = 
	    readGraph<symmetricVertex>(iFile,symmetric,binary);
	PageRank(G, maxIter);
	G.del(); 
    } else {
	graph<asymmetricVertex> G = 
	    readGraph<asymmetricVertex>(iFile,symmetric,binary);
	PageRank(G, maxIter);
	G.del();
    }
    return 0;
}
