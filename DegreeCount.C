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
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (needOutDegree) {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();}
    } else {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree();}
    }

    int accum[numOfShards];
    int partitionResult[numOfShards];
    for (int i = 0; i < numOfShards; i++) {
	accum[i] = 0;
	partitionResult[i] = 0;
    }

    long totalDegree = 0;
    for (intT i = 0; i < n; i++) {
	totalDegree += degrees[i];
    }

    int averageDegree = totalDegree / numOfShards;
    int counter = 0;
    int tmpSizeCounter = 0;
    for (intT i = 0; i < n; i+=PAGESIZE/sizeof(double)) {
	for (intT j = 0; j < PAGESIZE / sizeof(double); j++) {
	    if (i + j >= n)
		break;
	    accum[counter] += degrees[i + j];
	    partitionResult[counter] += 1;
	    tmpSizeCounter++;
	}
	if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
	    counter++;
	    cout << tmpSizeCounter / (double)(PAGESIZE / sizeof(double)) << endl;
	    tmpSizeCounter = 0;
	}
    }
    
    cout << tmpSizeCounter / (double)(PAGESIZE / sizeof(double));
    
    for (int i = 0; i < numOfShards; i++) {
	cout << accum[i] << endl;
    }

    int sum = 0;
    for (int i = 0; i < numOfShards; i++) {
	sum += partitionResult[i];
	cout << partitionResult[i] << " " << partitionResult[i]/(double)(PAGESIZE / sizeof(double)) << endl;
    }
    cout << sum << endl;
    cout << GA.n << endl;
    free(degrees);
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
