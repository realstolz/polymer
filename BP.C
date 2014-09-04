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
    BP_F(EdgeWeight *_edgeW, EdgeData *_edgeD_curr, EdgeData *_edgeD_next, VertexInfo *_vertI, VertexData *_vertD_curr, VertexData *_vertD_next, intT *_offsets) : 
	edgeW(_edgeW), edgeD_curr(_edgeD_curr), edgeD_next(_edgeD_next), vertI(_vertI), vertD_curr(_vertD_curr), vertD_next(_vertD_next), offsets(_offsets) {}
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
	intT dstIdx = offsets[s] + edgeIdx;
	//printf("idx: %d\n", dstIdx);
	for (int i = 0; i < NSTATES; i++) {
	    edgeD_next[dstIdx].belief[i] = 0.0;
	    for (int j = 0; j < NSTATES; j++) {
		edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
	    }
	    writeMult(&(vertD_next[d].product[i]), edgeD_next[dstIdx].belief[i]);
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
/*
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
*/

template <class vertex>
void BeliefPropagation(graph<vertex> GA, int maxIter = -1) {
    //GenLineGraph(GA);
    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    
    srand(time(NULL));

    printf("start init\n");

    intT *degrees = newA(intT, n);
    intT *offsets = newA(intT, n);

    {parallel_for (intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();}

    offsets[0] = 0;
    for (intT i = 1; i < n; i++) {
	offsets[i] = offsets[i-1] + degrees[i-1];
    }

    printf("%d %d %d\n", n, offsets[n-1] + degrees[n-1], GA.m);
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
		    printf("not symmetric\n");
		}
	    }
	}
    }
    */
    printf("check over\n");
  
    EdgeWeight *edgeW = (EdgeWeight *)malloc(sizeof(EdgeWeight) * GA.m);
    EdgeData *edgeD_curr = (EdgeData *)malloc(sizeof(EdgeData) * GA.m);
    EdgeData *edgeD_next = (EdgeData *)malloc(sizeof(EdgeData) * GA.m);
    {parallel_for (intT i = 0; i < GA.m; i++) {
	    for (int j = 0; j < NSTATES; j++) {
		for (int k = 0; k < NSTATES; k++) {
		    edgeW[i].potential[j][k] = rand() % 10000 / 20000.0;
		}
		edgeD_curr[i].belief[j] = 0.5;
	    }
	}
    }

    VertexInfo *vertI = (VertexInfo *)malloc(sizeof(VertexInfo) * n);
    VertexData *vertD_curr = (VertexData *)malloc(sizeof(VertexData) * n);
    VertexData *vertD_next = (VertexData *)malloc(sizeof(VertexData) * n);
    {parallel_for (intT i = 0; i < n; i++) {
	    for (int j = 0; j < NSTATES; j++) {
		vertI[i].potential[j] = rand() % 10000 / 10000.0;
		vertD_curr[i].product[j] = 1.0;
	    }
	}
    }
    printf("end init\n");

    bool* frontier = newA(bool,n);
  
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    
    {parallel_for(intT i=0;i<n;i++) frontier[i] = 1;}

    vertices Frontier(n,n,frontier);
    Frontier.toDense();
    startTime();
  
    intT round = 0;
    while(1){
	if (maxIter > 0 && round >= maxIter)
	    break;
	round++;

	//vertices output = edgeMap(GA, Frontier, BP_F<vertex>(product,m_next,GA.V),GA.m/20,DENSE_FORWARD);

	vertexMap(Frontier,BP_Vertex_Reset(vertD_next));
	edgeMapDenseBP(GA, Frontier.d, BP_F<vertex>(edgeW, edgeD_curr, edgeD_next, vertI, vertD_curr, vertD_next, offsets));
	swap(edgeD_curr, edgeD_next);
	swap(vertD_curr, vertD_next);
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
