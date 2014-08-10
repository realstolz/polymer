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
#include "ligra.h"
#include "gettime.h"
#include "math.h"
using namespace std;

bool needResult = false;

template <class ET>
inline void writeDiv(ET *a, ET b) {
  volatile ET newV, oldV; 
  do {oldV = *a; newV = oldV / b;}
  while (!CAS(a, oldV, newV));
}

template <class vertex>
struct BP_F {
    double** m_next;
    double* product;
    vertex* V;
    BP_F(double* _product, double** _p_next, vertex* _V) : 
	product(_product), m_next(_p_next), V(_V) {}
    inline bool update(intT s, intT d){ 
	(m_next[s])[d] = product[s];
	return 1;
    }
    inline bool updateAtomic (intT s, intT d) { //atomic Update
	(m_next[s])[d] = product[s];
	return 1;
    }
    inline bool cond (intT d) { return 1; } //does nothing
};

template <class vertex>
struct BP_F_2 {
    double** m_next;
    double** m_curr;
    vertex* V;
    BP_F_2(double** _p_curr, double** _p_next, vertex* _V) : 
	m_curr(_p_curr), m_next(_p_next), V(_V) {}
    inline bool update(intT s, intT d){ 
	(m_next[s])[d] /= (m_curr[d])[s];
	return 1;
    }
    inline bool updateAtomic (intT s, intT d) { //atomic Update
	writeDiv(&((m_next[s])[d]), (m_curr[d])[s]);
	return 1;
    }
    inline bool cond (intT d) { return 1; } //does nothing
};

//resets p
struct BP_Vertex_Reset {
    double* p_curr;
    BP_Vertex_Reset(double* _p_curr) :
    p_curr(_p_curr) {}
    inline bool operator () (intT i) {
	p_curr[i] = 0.0;
	return 1;
    }
};

template <class vertex>
graph<vertex> GenLineGraph(graph<vertex> &GA) {
    vertex *newVertexSet = (vertex *)malloc(sizeof(vertex) * GA.m);
    intT *newDegrees = (intT *)malloc(sizeof(intT) * GA.m);
    intT *newOffsets = (intT *)malloc(sizeof(intT) * GA.m);
    intT *degrees = (intT *)malloc(sizeof(intT) * GA.n);
    intT *offsets = (intT *)malloc(sizeof(intT) * GA.n);
    {parallel_for(intT i = 0; i < GA.n; i++) degrees[i] = GA.V[i].getOutDegree();}
    offsets[0] = 0;
    for (intT i = 1; i < GA.n; i++) {
	offsets[i] =  offsets[i-1] + degrees[i-1];
    }
    {parallel_for(intT i = 0; i < GA.n; i++) {
	    intT d = degrees[i];
	    intT offset = offsets[i];
	    for (intT j = 0; j < d; j++) {
		intT ngh = GA.V[i].getOutNeighbor(j);		
		intT subDeg = degrees[ngh];
		newVertexSet[offset + j].setOutDegree(subDeg);
		newDegrees[offset + j] = subDeg;
	    }
	}
    }
    printf("OK\n");
    intT totalDeg = newDegrees[0];//sequence::plusScan(newOffsets, newDegrees, GA.m);
    newOffsets[0] = 0;
    for (intT i = 1; i < GA.m; i++) {
	newOffsets[i] = newOffsets[i-1] + newDegrees[i-1];
	totalDeg += newDegrees[i];
    }
    printf("Total Edge in Line Graph: %d\n", totalDeg);
    
    intE *edges = (intE *)malloc(sizeof(intE) * totalDeg);

    {parallel_for(intT i = 0; i < GA.n; i++) {
	    intT d = degrees[i];
	    intT offset = offsets[i];
	    for (intT j = 0; j < d; j++) {
		intT ngh = GA.V[i].getOutNeighbor(j);		
		intT subDeg = degrees[ngh];
		intT nghOffset = offsets[ngh];
		intT edgeOffset = newOffsets[offset + j];
		intE *localEdges = &edges[edgeOffset];
		intT counter = 0;
		for (intT k = 0; k < subDeg; k++) {
		    localEdges[k] = nghOffset + k;
		}
		newVertexSet[offset + j].setOutNeighbors(localEdges);
	    }
	}
    }
    return graph<vertex>(newVertexSet, GA.m, totalDeg);
}

template <class vertex>
void BeliefPropagation(graph<vertex> GA, int maxIter = -1) {
    //GenLineGraph(GA);
    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    
    srand(time(NULL));
  
    double one_over_n = 1/(double)n;
    double **m_curr = (double **)malloc(sizeof(double *) * n);
    double **m_next = (double **)malloc(sizeof(double *) * n);
    double *product = (double *)malloc(sizeof(double) * n);
    printf("start init\n");
    for (intT i = 0; i < n; i++) {
	m_curr[i] = (double *)malloc(sizeof(double) * n);
	m_next[i] = (double *)malloc(sizeof(double) * n);
	{parallel_for (intT j = 0; j < n; j++) {
		(m_curr[i])[j] = 0.5;//(rand() % 1000) / 1000.0;
		(m_next[i])[j] = 0.5;//(rand() % 1000) / 1000.0;
	    }
	}
	product[i] = (rand() % 1000) / 1000.0;
    }
    printf("end init\n");
    bool* frontier = newA(bool,n);
  
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    
    {parallel_for(intT i=0;i<n;i++) frontier[i] = 1;}

    vertices Frontier(n,n,frontier);
    startTime();
  
    intT round = 0;
    while(1){
	if (maxIter > 0 && round >= maxIter)
	    break;
	round++;

	vertices output = edgeMap(GA, Frontier, BP_F<vertex>(product,m_next,GA.V),GA.m/20,DENSE_FORWARD);
	edgeMap(GA, Frontier, BP_F_2<vertex>(m_curr,m_next,GA.V),GA.m/20,DENSE_FORWARD);

	swap(m_curr, m_next);
    }
    cout<<"Finished in "<<round<<" iterations\n";
    Frontier.del();
    nextTime("BeliefPropagation");

    if (needResult) {
	for (intT i = 0; i < GA.n; i++) {
	    //cout << i << "\t" << std::scientific << std::setprecision(9)<< p_curr[i] << "\n";
	}
    }

    //free(p_curr); free(p_next); 
    printf("total init time: %lf\n", mapTime);
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
	BeliefPropagation(G, maxIter);
	G.del(); 
    } else {
	graph<asymmetricVertex> G = 
	    readGraph<asymmetricVertex>(iFile,symmetric,binary);
	BeliefPropagation(G, maxIter);
	G.del();
    }
}
