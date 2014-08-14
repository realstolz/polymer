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

#define PAGESIZE (4096)

bool needOutDegree = false;

template <class vertex>
void countDegree(graph<vertex> GA, int numOfShards) {
    const intT n = GA.n;

    for (int i = 0; i < n; i++) {
	if (GA.V[i].getOutDegree() + GA.V[i].getInDegree() <= 0)
	    continue;
	//printf("%d", i);
	for (int j = 0; j < GA.V[i].getOutDegree(); j++) {
	    printf("%d\t%d\n", i, GA.V[i].getOutNeighbor(j));
	}
	//printf("\n");
    }
}

int parallel_main(int argc, char* argv[]) {  
  char* iFile;
  bool binary = false;
  bool symmetric = false;
  needOutDegree = false;
  int numOfShards = -1;
  if(argc > 1) iFile = argv[1];
  if(argc > 2) numOfShards = atoi(argv[2]);
  if(argc > 3) if((string) argv[3] == (string) "-outDeg") needOutDegree = true;
  if(argc > 4) if((string) argv[4] == (string) "-s") symmetric = true;
  if(argc > 5) if((string) argv[5] == (string) "-b") binary = true;
  
  if(symmetric) {
    graph<symmetricVertex> G = 
      readGraph<symmetricVertex>(iFile,symmetric,binary);
    countDegree(G, numOfShards);
    G.del(); 
  } else {
    graph<asymmetricVertex> G = 
      readGraph<asymmetricVertex>(iFile,symmetric,binary);
    countDegree(G, numOfShards);
    G.del();
  }
}
