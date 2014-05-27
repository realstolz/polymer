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
#include "ligra-numa.h"
#include "gettime.h"
#include "math.h"

#include <pthread.h>
#include <numa.h>
using namespace std;

double **p_curr_list = NULL;
int vPerNode = 0;
int numOfNode = 0;

bool needResult = false;

double getCurrP(int index) {
    int nodeNo = index / vPerNode;
    int localInd = index % vPerNode;
    if (nodeNo >= numOfNode) {
	nodeNo = numOfNode - 1;
	localInd = localInd + vPerNode;
    }
    return *(*(p_curr_list+nodeNo) + localInd);
}

template <class vertex>
struct PR_F {
    double* p_curr, *p_next;
    vertex* V;
    int rangeLow;
    int rangeHi;
    PR_F(double* _p_curr, double* _p_next, vertex* _V, int _rangeLow, int _rangeHi) : 
	p_curr(_p_curr), p_next(_p_next), V(_V), rangeLow(_rangeLow), rangeHi(_rangeHi) {}
    inline bool update(intT s, intT d){ //update function applies PageRank equation
	p_next[d-rangeLow] += getCurrP(s)/V[s].getOutDegree();
	return 1;
    }
    inline bool updateAtomic (intT s, intT d) { //atomic Update
	writeAdd(&p_next[d-rangeLow],getCurrP(s)/V[s].getOutDegree());
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
};

pthread_barrier_t barr;
pthread_mutex_t mut;

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
    /*
    if (0 != __cilkrts_set_param("nworkers","1")) {
	printf("set failed: %d\n", tid);
    }
    */

    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    int numOfT = my_arg->numOfNode;
    
    int blockSize = n / numOfT;

    int rangeLow = tid * blockSize;
    int rangeHi = (tid + 1) * blockSize;
    if (tid == numOfT - 1) {
	rangeHi = rangeHi + n % numOfT;
	blockSize = blockSize + n % numOfT;
    }
    double one_over_n = 1/(double)n;

    
    double* p_curr = (double *)numa_alloc_local(sizeof(double) * blockSize);
    double* p_next = (double *)numa_alloc_local(sizeof(double) * blockSize);
    bool* frontier = (bool *)numa_alloc_local(sizeof(bool) * n);
    
    /*
    double* p_curr = (double *)malloc(sizeof(double) * blockSize);
    double* p_next = (double *)malloc(sizeof(double) * blockSize);
    bool* frontier = (bool *)malloc(sizeof(bool) * n);
    */
    

    p_curr_list[tid] = p_curr;
    /*
    if (tid == 0)
	startTime();
    */
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    {parallel_for(intT i=0;i<blockSize;i++) p_curr[i] = one_over_n;}
    {parallel_for(intT i=0;i<blockSize;i++) p_next[i] = 0;} //0 if unchanged
    {parallel_for(intT i=0;i<n;i++) frontier[i] = 1;}
    
    vertices Frontier(n,n,frontier);

    pthread_barrier_wait(&barr);

    intT round = 0;
    while(1){
	if (maxIter > 0 && round >= maxIter)
	    break;
	round++;

	//gettimeofday(&start, &tz);
	vertices output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),GA.m/20,DENSE_FORWARD);

	//gettimeofday(&end, &tz);
	//double timeStart = ((double)start.tv_sec) + ((double)start.tv_usec) / 1000000.0;
	//double timeEnd = ((double)end.tv_sec) + ((double)end.tv_usec) / 1000000.0;
	//mapTime = timeEnd - timeStart;

	//output.del();

	double tmpConst = (1-damping)*(1/(double)n);

	{parallel_for (intT i = 0; i < blockSize; i++) {
		p_next[i] = damping * p_next[i] + tmpConst;
	    }}
	
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
	//vertexMap(Frontier,PR_Vertex_Reset(p_curr));
	{parallel_for (intT i = 0; i < blockSize; i++) {
		p_curr[i] = p_next[i];
		p_next[i] = 0.0;
	    }}
	
	//swap(p_curr,p_next);
	//Frontier.del(); 
	//Frontier = output;
	pthread_barrier_wait(&barr);

	/*
	  if (tid == 0)
	    printf("round %d: %lf\n", round, mapTime);
	*/
    }
    //cout<<"Finished in "<<round<<" iterations\n";
    //Frontier.del();

    if (needResult) {
	for (intT i = 0; i < blockSize; i++) {
	    pthread_mutex_lock(&mut);
	    cout << i+rangeLow << "\t" << std::scientific << p_curr[i] << "\n";
	    pthread_mutex_unlock(&mut);
	}
    }
    

    numa_free(p_curr, blockSize); 
    numa_free(p_next, blockSize); 
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
    numOfNode = 6;//numa_num_configured_nodes();
    vPerNode = GA.n / numOfNode;
    pthread_barrier_init(&barr, NULL, numOfNode);
    pthread_mutex_init(&mut, NULL);
    printf("start create %d threads\n", numOfNode);
    pthread_t tids[numOfNode];
    p_curr_list = (double **)malloc(sizeof(double *) * numOfNode);
    startTime();
    for (int i = 0; i < numOfNode; i++) {
	PR_worker_arg *arg = (PR_worker_arg *)malloc(sizeof(PR_worker_arg));
	arg->GA = (void *)(&GA);
	arg->maxIter = maxIter;
	arg->tid = i;
	arg->numOfNode = numOfNode;
	pthread_create(&tids[i], NULL, PageRankThread<vertex>, (void *)arg);
    }
    printf("all created\n");
    for (int i = 0; i < numOfNode; i++) {
	pthread_join(tids[i], NULL);
    }
    nextTime("PageRank");
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
