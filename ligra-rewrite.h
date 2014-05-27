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

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <sys/mman.h>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "IO.h"

#include <numa.h>

using namespace std;

#define PAGESIZE (4096)

//*****START FRAMEWORK*****

/* This is no longer ligra, 
the vertices set should be set by each node
*/

int roundUp(double x) {
    int ones = x / 1;
    double others = x - ones;
    if (others > 0) {
	ones++;
    }
    return ones;
}

template <class vertex>
void partitionByDegree(graph<vertex> GA, int numOfShards, int *sizeArr, int sizeOfOneEle, bool useOutDegree=false) {
    const intT n = GA.n;
    int *degrees = newA(int, n);

    int shardSize = n / numOfShards;

    if (useOutDegree) {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();}
    } else {
	{parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree();}
    }

    int accum[numOfShards];
    for (int i = 0; i < numOfShards; i++) {
	accum[i] = 0;
	sizeArr[i] = 0;
    }

    long totalDegree = 0;
    for (intT i = 0; i < n; i++) {
	totalDegree += degrees[i];
    }

    int averageDegree = totalDegree / numOfShards;
    int counter = 0;
    int tmpSizeCounter = 0;
    for (intT i = 0; i < n; i+=PAGESIZE/sizeOfOneEle) {
	for (intT j = 0; j < PAGESIZE / sizeOfOneEle; j++) {
	    if (i + j >= n)
		break;
	    accum[counter] += degrees[i + j];
	    sizeArr[counter]++;
	    tmpSizeCounter++;
	}
	if (accum[counter] >= averageDegree && counter < numOfShards - 1) {
	    counter++;
	    cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle) << endl;
	    tmpSizeCounter = 0;
	}
    }
    
    cout << tmpSizeCounter / (double)(PAGESIZE / sizeOfOneEle);
    cout << GA.n << endl;
    free(degrees);
}

void *mapDataArray(int numOfShards, int *sizeArr, int sizeOfOneEle) {
    int numOfPages = 0;
    for (int i = 0; i < numOfShards; i++) {
	cout << sizeArr[i] << endl;
        numOfPages += sizeArr[i] / (double)(PAGESIZE / sizeOfOneEle);
    }
    numOfPages++;

    void *toBeReturned = mmap(NULL, numOfPages * PAGESIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (toBeReturned == NULL) {
	cout << "OOps" << endl;
    }
    
    int offset = 0;
    for (int i = 0; i < numOfShards; i++) {
	void *startPos = (void *)((char *)toBeReturned + offset * sizeOfOneEle);
	//printf("start binding %d : %p\n", i, startPos);
	numa_tonode_memory(startPos, sizeArr[i], i);
	offset = offset + sizeArr[i];
    }
    return toBeReturned;
}


//*****VERTEX OBJECT*****
struct vertices {
    intT n, m;
    int numOfNodes;
    int *numOfVertexOnNode;
    int *offsets;
    bool** d;
    
    vertices(int _numOfNodes) {
	this->numOfNodes = _numOfNodes;
	d = (bool **)malloc(numOfNodes * sizeof(bool*));
	numOfVertexOnNode = (int *)malloc(numOfNodes * sizeof(int));
	offsets = (int *)malloc(numOfNodes * sizeof(int));
    }

    void registerArr(int nodeNum, bool* arr, int size) {
	d[nodeNum] = arr;
	numOfVertexOnNode[nodeNum] = size;
    }

    void calculateOffsets() {
	offsets[0] = 0;
	for (int i = 1; i < numOfNodes; i++) {
	    offsets[i] = numOfVertexOnNode[i-1] + offsets[i-1];
	}
    }

    int getSize(int nodeNum) {
	return numOfVertexOnNode[nodeNum];
    }

    int getOffset(int nodeNum) {
	return offsets[nodeNum];
    }

    bool getBit(int index) {
	int accum = 0;
	int i = 0;
        while (index >= accum + numOfVertexOnNode[i]) {
	    accum += numOfVertexOnNode[i];
	    i++;
	}
	return *(d[i] + (index - accum));
    }

    bool *getArr(int nodeNum) {
	return d[nodeNum];
    }

    bool eq (vertices& b) {
    }
    
    void print() {
    }

    void del() {
	free(offsets);
	free(numOfVertexOnNode);
	free(d);
    }
};

struct LocalFrontier {
    int startID;
    int endID;
    bool *b;
    
    LocalFrontier(bool *_b, int start, int end):b(_b), startID(start), endID(end){}
    
    bool inRange(int index) { return (startID <= index && index < endID);}
    void setBit(int index, bool val) { b[index-startID] = val;}

    bool *swapBitVector(bool *newB) {
	bool *tmp = b;
	b = newB;
	return tmp;
    }
};

struct nonNegF{bool operator() (intT a) {return (a>=0);}};

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options {
  DENSE, DENSE_PARALLEL, DENSE_FORWARD
};




//*****EDGE FUNCTIONS*****

template <class F, class vertex>
bool* edgeMapDense(graph<vertex> GA, vertices* frontier, F f, LocalFrontier *next, bool parallel = 0) {
    intT numVertices = GA.n;
    vertex *G = GA.V;
    for (intT i=next->startID; i<next->endID; i++){
	next->setBit(i, false);
	if (f.cond(i)) { 
	    intT d = G[i].getInDegree();
	    for(intT j=0; j<d; j++){
		intT ngh = G[i].getInNeighbor(j);
		if (frontier->getBit(ngh) && f.update(ngh,i)) next->setBit(i, true);
		if(!f.cond(i)) break;
	    }
	}
    }
    return NULL;
}

template <class F, class vertex>
bool* edgeMapDenseForward(graph<vertex> GA, vertices *frontier, F f, LocalFrontier *next) {
    intT numVertices = GA.n;
    vertex *G = GA.V;
    {parallel_for(long i=next->startID;i<next->endID;i++) next->setBit(i, false);}
    int currNodeNum = 0;
    bool *currBitVector = frontier->getArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    for (long i=0; i<numVertices; i++){
	if (i == nextSwitchPoint) {
	    currOffset += frontier->getSize(currNodeNum);
	    nextSwitchPoint += frontier->getSize(currNodeNum + 1);
	    currNodeNum++;
	    currBitVector = frontier->getArr(currNodeNum);
	}
	if (currBitVector[i-currOffset]) {
	    intT d = G[i].getOutDegree();
	    for(intT j=0; j<d; j++){
		uintT ngh = G[i].getOutNeighbor(j);
		if (next->inRange(ngh) && f.cond(ngh) && f.updateAtomic(i,ngh)) {
		    next->setBit(ngh, true);
		}
	    }
	}
    }
    return NULL;
}

static int edgesTraversed = 0;

void switchFrontier(int nodeNum, vertices *V, LocalFrontier *next) {
    bool *current = V->getArr(nodeNum);
    intT size = V->getSize(nodeNum);
    for (intT i = 0; i < size; i++) {
	current[i] = false;
    }
    
    bool *tmp = next->swapBitVector(V->getArr(nodeNum));
    V->registerArr(nodeNum, tmp, V->getSize(nodeNum));
}

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
void edgeMap(graph<vertex> GA, vertices *V, F f, LocalFrontier *next, intT threshold = -1, 
		 char option=DENSE, bool remDups=false) {
    intT numVertices = GA.n;
    uintT numEdges = GA.m;
    vertex *G = GA.V;
    //intT m = V.numNonzeros();
    
    /*
      if (numVertices != V.numRows()) {
      cout << "edgeMap: Sizes Don't match" << endl;
      abort();
      }
    */
    
    /*
    // used to generate nonzero indices to get degrees
    uintT* degrees = newA(uintT, m);
    vertex* frontierVertices;
    V.toSparse();
    frontierVertices = newA(vertex,m);
    {parallel_for (long i=0; i < m; i++){
    vertex v = G[V.s[i]];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
    uintT outDegrees = sequence::plusReduce(degrees, m);
    edgesTraversed += outDegrees;
    if (outDegrees == 0) return vertices(numVertices);
  */
    
    
    bool* R = (option == DENSE_FORWARD) ? 
	edgeMapDenseForward(GA, V, f, next) : 
	edgeMapDense(GA, V, f, next, option);
    //cout << "size (D) = " << v1.m << endl;    
}

//*****VERTEX FUNCTIONS*****

template <class F>
void vertexMap(vertices *V, F add, int nodeNum) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    for (int i = 0; i < size; i++) {
	if (b[i])
	    add(i + offset);
    }
}

template <class F>
void vertexFilter(vertices *V, F filter, int nodeNum, bool *result) {
    int size = V->getSize(nodeNum);
    int offset = V->getOffset(nodeNum);
    bool *b = V->getArr(nodeNum);
    for (int i = 0; i < size; i++) {
	result[i] = false;
	if (b[i])
	    result[i] = filter(i + offset);
    }
}

//*****WEIGHTED EDGE FUNCTIONS*****
/*
template <class F, class vertex>
  bool* edgeMapDense(wghGraph<vertex> GA, bool* vertices, F f, bool parallel = 0) {
  intT numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for (long i=0; i<numVertices; i++){
    next[i] = 0;
    if (f.cond(i)) {
      intT d = G[i].getInDegree();
      if(!parallel || d < 1000) {
	for(intT j=0; j<d; j++){
	  intT ngh = G[i].getInNeighbor(j);
	  if(vertices[ngh]){
	    intT edgeLen = G[i].getInWeight(j);
	    if (vertices[ngh] && f.update(ngh,i,edgeLen)) next[i] = 1;
	  }
	  if (!f.cond(i)) break;
	}
      }
      else {
	{parallel_for(intT j=0; j<d; j++){
	  intT ngh = G[i].getInNeighbor(j);
	  if(vertices[ngh]){
	    intT edgeLen = G[i].getInWeight(j);
	    if (vertices[ngh] && f.update(ngh,i,edgeLen)) next[i] = 1;
	  }
	  }}
      }
    }
    }
  }
  return next;
}

template <class F, class vertex>
bool* edgeMapDenseForward(wghGraph<vertex> GA, bool* vertices, F f) {
  intT numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for(long i=0;i<numVertices;i++) next[i] = 0;}
  {parallel_for (long i=0; i<numVertices; i++){
    if (vertices[i]) {
      intT d = G[i].getOutDegree();
      if(d < 1000) {
	for(intT j=0; j<d; j++){
	  intT ngh = G[i].getOutNeighbor(j);
	  intT edgeLen = G[i].getOutWeight(j);
	  if (f.cond(ngh) && f.updateAtomic(i,ngh,edgeLen)) next[ngh] = 1;
	}
      }
      else{
	parallel_for(intT j=0; j<d; j++){
	  intT ngh = G[i].getOutNeighbor(j);
	  intT edgeLen = G[i].getOutWeight(j);
	  if (f.cond(ngh) && f.updateAtomic(i,ngh,edgeLen)) next[ngh] = 1;
	}
      }
    }
  }}
  return next;
}

template <class F, class vertex>
  pair<uintT,intT*> edgeMapSparseW(vertex* frontierVertices, intT* indices, 
				   uintT* degrees, uintT m, F f, 
				   intT remDups=0, intT* flags=NULL) {
  uintT* offsets = degrees;
  uintT outEdgeCount = sequence::plusScan(offsets, degrees, m);

  intT* outEdges = newA(intT,outEdgeCount);

  {parallel_for (long i = 0; i < m; i++) {
    intT v = indices[i];
    uintT o = offsets[i];
    vertex vert = frontierVertices[i];
    intT d = vert.getOutDegree();
    if(d < 1000) {
      for (intT j=0; j < d; j++) {
	intT ngh = vert.getOutNeighbor(j);
	intT edgeLen = vert.getOutWeight(j);
	if (f.cond(ngh) && f.updateAtomic(v,ngh,edgeLen)) 
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = -1;
      }
    }
    else {
      {parallel_for (intT j=0; j < d; j++) {
	intT ngh = vert.getOutNeighbor(j);
	intT edgeLen = vert.getOutWeight(j);
	if (f.cond(ngh) && f.updateAtomic(v,ngh,edgeLen)) 
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = -1;
	
	}}
    }
    }}

  intT* nextIndices = newA(intT, outEdgeCount);
  if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);

  // Filter out the empty slots (marked with -1)
  uintT nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonNegF());
  free(outEdges);
  return pair<uintT,intT*>(nextM, nextIndices);
}

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
  vertices edgeMap(wghGraph<vertex> GA, vertices V, F f, intT threshold = -1, char option=DENSE, bool remDups=false) {
  intT numVertices = GA.n;
  uintT numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  intT m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }

  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(vertex,m);

  {parallel_for (long i=0; i < m; i++){
    vertex v = G[V.s[i]];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  uintT outDegrees = sequence::plusReduce(degrees, m);    
  if (outDegrees == 0) return vertices(numVertices);
  if (m + outDegrees > threshold) { 
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R = option == DENSE_FORWARD ? 
      edgeMapDenseForward(GA,V.d,f) : 
      edgeMapDense(GA, V.d, f,option);
    vertices v1 = vertices(numVertices, R);
    //cout << "size (D) = " << v1.m << endl;
    return  v1;
  } else { 
    pair<uintT,intT*> R = 
      remDups ? 
      edgeMapSparseW(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparseW(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    //cout << "size (S) = " << R.first << endl;
    free(degrees);
    free(frontierVertices);
    return vertices(numVertices, R.first, R.second);
  }
}
*/
