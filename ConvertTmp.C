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
using namespace std;

struct BF_F {
    int* ShortestPathLen;
    int* Visited;
    BF_F(int* _ShortestPathLen, int* _Visited) : 
	ShortestPathLen(_ShortestPathLen), Visited(_Visited) {}
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

template <class vertex>
int* BellmanFord(intT start, wghGraph<vertex> GA) {
    for (intT i = 0; i < GA.n; i++) {
	if (GA.V[i].getOutDegree() + GA.V[i].getInDegree() <= 0)
	    continue;
	for (intT j = 0; j < GA.V[i].getOutDegree(); j++) {
	    printf("%d\t%d\t%d\n", i, GA.V[i].getOutNeighbor(j), GA.V[i].getOutWeight(j));
	}
    }
    return NULL;
}

int parallel_main(int argc, char* argv[]) {  
    char* iFile;
    bool binary = false;
    bool symmetric = false;
    intT startPos = 0;
    if(argc > 1) iFile = argv[1];
    if(argc > 2) startPos = atoi(argv[2]);
    if(argc > 3) if((string) argv[3] == (string) "-s") symmetric = true;
    if(argc > 4) if((string) argv[4] == (string) "-b") binary = true;
  
    if(symmetric) {
	wghGraph<symmetricWghVertex> WG = 
	    readWghGraph<symmetricWghVertex>(iFile,symmetric,binary);
	BellmanFord((intT)startPos,WG);
	WG.del(); 
    } else {
	wghGraph<asymmetricWghVertex> WG = 
	    readWghGraph<asymmetricWghVertex>(iFile,symmetric,binary);
	BellmanFord((intT)startPos,WG);
	WG.del();
    }
}
