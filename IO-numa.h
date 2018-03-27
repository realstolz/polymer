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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
#include "quickSort.h"
using namespace std;

typedef pair<uintE,uintE> intPair;

typedef pair<uintE, pair<uintE,intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator() (pair<intE,E> a, pair<intE,E> b) {
    return a.first < b.first;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long long n; // total number of characters
  char* Chars;  // array storing all strings
  long long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long long nn, char** S, long long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};
 
inline bool isSpace(char c) {
  switch (c)  {
  case '\r': 
  case '\t': 
  case '\n': 
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file (fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long long end = file.tellg();
  file.seekg (0, ios::beg);
  long long n = end - file.tellg();
  char* bytes = newA(char,n+1);
  //printf("end is: %lu, n is: %lu, bytes is: %p\n", end, n, bytes);
  file.read (bytes,n);
  file.close();
  return _seq<char>(bytes,n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long long n) {
  {parallel_for (long long i=0; i < n; i++) 
      if (isSpace(Str[i])) Str[i] = 0; }

  // mark start of words
  bool *FL = newA(bool,n);
  FL[0] = Str[0];
  {parallel_for (long long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];}
    
  // offset for each start of word
  _seq<long long> Off = sequence::packIndex<long long>(FL, n);
  long long m = Off.n;
  long long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char*, m);
  {parallel_for (long long j=0; j < m; j++) SA[j] = Str+offsets[j];}

  free(offsets); free(FL);
  return words(Str,n,SA,m);
}

template <class vertex>
graph<vertex> readGraphFromFile(char* fname, bool isSymmetric) {
  _seq<char> S = readStringFromFile(fname);
  words W = stringToWords(S.A, S.n);
  if (W.Strings[0] != (string) "AdjacencyGraph") {
    cout << "Bad input file" << endl;
    abort();
  }

  long long len = W.m -1;
  long long n = atol(W.Strings[1]);
  long long m = atol(W.Strings[2]);
  if (len != n + m + 2) {
    cout << "Bad input file" << endl;
    abort();
  }

  intT* offsets = newA(intT,n);
  intE* edges = newA(intE,m);

  {parallel_for(long long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long long i=0; i<m; i++) edges[i] = atol(W.Strings[i+n+3]); }
  //W.del(); // to deal with performance bug in malloc
    
  vertex* v = newA(vertex,n);

  {parallel_for (uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l); 
    v[i].setOutNeighbors(edges+o);     
    }}

  if(!isSymmetric) {
    intT* tOffsets = newA(intT,n);
    {parallel_for(intT i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
    intE* inEdges = newA(intE,m);
    intPair* temp = newA(intPair,m);
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(intT j=0;j<v[i].getOutDegree();j++){	  
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
      }
      }}
    free(offsets);

    quickSort(temp,m,pairFirstCmp<intE>());
 
    tOffsets[0] = 0; inEdges[0] = temp[0].second;
    {parallel_for(intT i=1;i<m;i++) {
      inEdges[i] = temp[i].second;
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);

    uintT currOffset = m;
    for(intT i=n-1;i>=0;i--) {
      if(tOffsets[i] == INT_T_MAX) tOffsets[i] = currOffset;
      else currOffset = tOffsets[i];
    }

    {parallel_for(uintT i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors(inEdges+o);
      }}    

    free(tOffsets);
    return graph<vertex>(v,(intT)n,m,edges,inEdges);
  }

  else {
    free(offsets);
    return graph<vertex>(v,(intT)n,m,edges);
  }
}

template <class vertex>
wghGraph<vertex> readWghGraphFromFile(char* fname, bool isSymmetric) {
  _seq<char> S = readStringFromFile(fname);
  words W = stringToWords(S.A, S.n);
  //printf("convert over\n");
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
    cout << "Bad input file" << endl;
    abort();
  }

  long long len = W.m -1;
  long long n = atol(W.Strings[1]);
  long long m = atol(W.Strings[2]);
  if (len != n + 2*m + 2) {
    cout << "Bad input file" << endl;
    abort();
  }

  intT* offsets = newA(intT,n);
  intE* edgesAndWeights = newA(intE,2*m);

  {parallel_for(long long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long long i=0; i<m; i++) {
      edgesAndWeights[2*i] = atol(W.Strings[i+n+3]); 
      edgesAndWeights[2*i+1] = atol(W.Strings[i+n+m+3]);
    } }
  W.del(); // to deal with performance bug in malloc

  vertex *v = newA(vertex,n);

  {parallel_for (uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
    v[i].setOutNeighbors((intE*)(edgesAndWeights+2*o));
    }}

  if(!isSymmetric) {
    intT* tOffsets = newA(intT,n);
    {parallel_for(intT i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
    intE* inEdgesAndWghs = newA(intE,2*m);
    intTriple* temp = newA(intTriple,m);
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(intT j=0;j<v[i].getOutDegree();j++){	  
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
      }
      }}
    free(offsets);
    quickSort(temp,m,pairFirstCmp<intPair>());

    tOffsets[0] = 0; 
    inEdgesAndWghs[0] = temp[0].second.first;
    inEdgesAndWghs[1] = temp[0].second.second;
    {parallel_for(long long i=1;i<m;i++) {
      inEdgesAndWghs[2*i] = temp[i].second.first; 
      inEdgesAndWghs[2*i+1] = temp[i].second.second;
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    //printf("offset over\n");
    free(temp);

    uintT currOffset = m;
    for(intT i=n-1;i>=0;i--) {
      if(tOffsets[i] == INT_T_MAX) tOffsets[i] = currOffset;
      else currOffset = tOffsets[i];
    }
    
    {parallel_for(uintT i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors((intE*)(inEdgesAndWghs+2*o));
      }}

    free(tOffsets);
    return wghGraph<vertex>(v,(intT)n,m,edgesAndWeights, inEdgesAndWghs);
  }

  else {  
    free(offsets);
    return wghGraph<vertex>(v,(intT)n,m,edgesAndWeights); 
  }
}

template <class vertex>
graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+7];
  char adjFile[strlen(iFile)+4];
  char idxFile[strlen(iFile)+4];
  strcpy(configFile,iFile);
  strcpy(adjFile,iFile);
  strcpy(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  intT n;
  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long long size = in2.tellg();
  in2.seekg(0);
  uintT m = size/sizeof(uint);
  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();
  
  uintE* edges = (uintE*) s;
  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(long)) { cout << "File size wrong\n"; abort(); }

  char* t = (char *) malloc(size);
  in3.read(t,size);
  in3.close();
  intT* offsets = (intT*) t;

  vertex* v = newA(vertex,n);
  
  {parallel_for(long long i=0;i<n;i++) {
    uintT o = offsets[i];
    uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i];
      v[i].setOutDegree(l); 
      v[i].setOutNeighbors((intE*)edges+o); }}

  cout << "n = "<<n<<" m = "<<m<<endl;

  if(!isSymmetric) {
    intT* tOffsets = newA(intT,n);
    {parallel_for(intT i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
    uintE* inEdges = newA(uintE,m);
    intPair* temp = newA(intPair,m);
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(intT j=0;j<v[i].getOutDegree();j++){
	temp[o+j].first = v[i].getOutNeighbor(j);
	temp[o+j].second = i;
      }
      }}

    quickSort(temp,m,pairFirstCmp<intE>());

    tOffsets[0] = 0; inEdges[0] = temp[0].second;
    {parallel_for(intT i=1;i<m;i++) {
      inEdges[i] = temp[i].second;
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);

    uintT currOffset = m;
    for(intT i=n-1;i>=0;i--) {
      if(tOffsets[i] == INT_T_MAX) tOffsets[i] = currOffset;
      else currOffset = tOffsets[i];
    }

    {parallel_for(uintT i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors((intE*)inEdges+o);
      }}
    free(tOffsets);
    return graph<vertex>(v,(intT)n,m,(intE*)edges, (intE*)inEdges);
  }
  free(offsets);
  return graph<vertex>(v,n,m,(intE*)edges);
}

template <class vertex>
wghGraph<vertex> readWghGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+7];
  char adjFile[strlen(iFile)+4];
  char idxFile[strlen(iFile)+4];
  strcpy(configFile,iFile);
  strcpy(adjFile,iFile);
  strcpy(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  intT n;
  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long long size = in2.tellg();
  in2.seekg(0);
  uintT m = size/sizeof(uint);

  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();

  intE* edges = (intE*) s;
  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(long)) { cout << "File size wrong\n"; abort(); }

  char* t = (char *) malloc(size);
  in3.read(t,size);
  in3.close();
  intT* offsets = (intT*) t;

  vertex *V = newA(vertex, n);
  intE* edgesAndWeights = newA(intE,2*m);
  {parallel_for(long long i=0;i<m;i++) {
    edgesAndWeights[2*i] = edges[i];
    edgesAndWeights[2*i+1] = 1; //give them unit weight
    }}
  free(edges);

  {parallel_for(long long i=0;i<n;i++) {
    uintT o = offsets[i];
    uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i];
    V[i].setOutDegree(l);
    V[i].setOutNeighbors(edgesAndWeights+2*o);
    }}
  cout << "n = "<<n<<" m = "<<m<<endl;

  if(!isSymmetric) {
    intT* tOffsets = newA(intT,n);
    {parallel_for(intT i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
    intE* inEdgesAndWghs = newA(intE,2*m);
    intPair* temp = newA(intPair,m);
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(intT j=0;j<V[i].getOutDegree();j++){
	temp[o+j] = make_pair(V[i].getOutNeighbor(j),i);
      }
      }}

    quickSort(temp,m,pairFirstCmp<intE>());

    tOffsets[0] = 0; 
    inEdgesAndWghs[0] = temp[0].second;
    inEdgesAndWghs[1] = 1;
    {parallel_for(intT i=1;i<m;i++) {
      inEdgesAndWghs[2*i] = temp[i].second;
      inEdgesAndWghs[2*i+1] = 1;
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);
    uintT currOffset = m;
    for(intT i=n-1;i>=0;i--) {
      if(tOffsets[i] == INT_T_MAX) tOffsets[i] = currOffset;
      else currOffset = tOffsets[i];
    }

    {parallel_for(uintT i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      V[i].setInDegree(l);
      V[i].setInNeighbors((intE*)(inEdgesAndWghs+2*o));
      }}
    free(tOffsets);
    return wghGraph<vertex>(V,(intT)n,m,edges,inEdgesAndWghs);
  }
  free(offsets);
  return wghGraph<vertex>(V,n,m,edgesAndWeights);
}

template <class vertex>
graph<vertex> readGraph(char* iFile, bool symmetric, bool binary) {
  if(binary) return readGraphFromBinary<vertex>(iFile,symmetric); 
  else return readGraphFromFile<vertex>(iFile,symmetric);
}

template <class vertex>
wghGraph<vertex> readWghGraph(char* iFile, bool symmetric, bool binary) {
  if(binary) return readWghGraphFromBinary<vertex>(iFile,symmetric); 
  else return readWghGraphFromFile<vertex>(iFile,symmetric);
}


template <class vertex>
graph<vertex> loadGraphFromBin(char *fileName) {
    int fd = open(fileName, O_RDONLY, S_IREAD);
    long long totalSize = 0;
    read(fd, (void *)&totalSize, sizeof(long long));
    void *buf = (void *)malloc(totalSize);
    printf("totalSize is: %lld %p\n", totalSize, buf);
    lseek(fd, 0, SEEK_SET);

    long long readSize = 0;
    while (readSize < totalSize) {
	long long readOnce = read(fd, (void *)((char *)buf + readSize), totalSize - readSize);
	readSize += readOnce;
	printf("read %lld, accum %lld\n", readOnce, readSize);
    }

    char *ptr = (char *)buf;
    ptr += sizeof(long long);

    intT n = *(intT *)ptr;
    ptr += sizeof(intT);
    
    long long m = *(long long *)ptr;
    ptr += sizeof(long long);

    printf("n & m: %" PRIintT " %lld\n", n, m);

    vertex *vertices = newA(vertex, n);
    intE *edges = newA(intE, m);
    intE *inEdges = newA(intE, m);
    long long counter = 0;
    long long inCounter = 0;

    for (intT i = 0; i < n; i++) {
	if (*(intT *)ptr != -1) {
	    printf("oops: %" PRIintT "\n", i);
	    abort();
	}
	ptr += sizeof(intT);

	intT vertIdx = *(intT *)ptr;
	if (vertIdx != i) {
	    printf("oops\n");
	}
	ptr += sizeof(intT);
	
	intT outDeg = *(intT *)ptr;
	ptr += sizeof(intT);

	intT inDeg = *(intT *)ptr;
	ptr += sizeof(intT);
	
	vertices[i].setOutDegree(outDeg);
	vertices[i].setOutNeighbors(&edges[counter]);
	vertices[i].setInDegree(inDeg);
	vertices[i].setInNeighbors(&inEdges[inCounter]);

	for (intT j = 0; j < outDeg; j++) {
	    intE dst = *(intE *)ptr;
	    edges[counter++] = dst;
	    ptr += sizeof(intE);
	}
	
	for (intT j = 0; j < inDeg; j++) {
	    intE dst = *(intE *)ptr;
	    inEdges[inCounter++] = dst;
	    ptr += sizeof(intE);
	}
    }
    
    free(buf);
    return graph<vertex>(vertices, (intT)n, m, edges);
}

template <class vertex>
wghGraph<vertex> loadWghGraphFromBin(char *fileName) {
}


template <class vertex>
void dumpGraphToBin(graph<vertex> &graph, char *fileName) {
    const intT n = graph.n;
    
    long long totalSize = 0;
    for (intT i = 0; i < n; i++) {
	totalSize += 4 * sizeof(intT); // const -1, vertID, outDeg, inDeg
	totalSize += graph.V[i].getOutDegree() * sizeof(intE);
	totalSize += graph.V[i].getInDegree() * sizeof(intE);
    }
    totalSize += sizeof(intT) + sizeof(long long) * 2;
    
    printf("Total size is: %ld\n", totalSize);
    
    void *buf = (void *)malloc(totalSize);
    
    char *ptr = (char *)buf;
    
    *(long long *)ptr = totalSize;
    ptr += sizeof(long long);
    
    *(intT *)ptr = n;
    ptr += sizeof(intT);
    
    *(long long *)ptr = graph.m;
    ptr += sizeof(long long);
    
    printf("write n & m: %d %d\n", n, graph.m);
    
    for (intT i = 0; i < n; i++) {
	*(intT *)ptr = -1;
	ptr += sizeof(intT);
	
	*(intT *)ptr = i;
	ptr += sizeof(intT);
	
	*(intT *)ptr = graph.V[i].getOutDegree();
	ptr += sizeof(intT);
	
	*(intT *)ptr = graph.V[i].getInDegree();
	ptr += sizeof(intT);
	
	for (intT j = 0; j < graph.V[i].getOutDegree(); j++) {
	    *(intE *)ptr = graph.V[i].getOutNeighbor(j);
	    ptr += sizeof(intE);
	}
	
	for (intT j = 0; j < graph.V[i].getInDegree(); j++) {
	    *(intE *)ptr = graph.V[i].getInNeighbor(j);
	    ptr += sizeof(intE);
	}
    }
    
    int fd = open(fileName, O_RDWR | O_CREAT, S_IWRITE | S_IREAD);
    long long written = 0;
    while (written < totalSize) {
	long long sizeWritten = write(fd, (void *)((char *)buf + written), totalSize - written);
	if (sizeWritten < 0) {
	    printf("oops\n");
	}	    
	written += sizeWritten;
	printf("wrote: %ld, accum: %ld\n", sizeWritten, written);
    }
}

template<class vertex>
void dumpPartitionInfo(graph<vertex> &graph, char * fileName, intT *sizeArr, int numOfNodes) {
    int fd = open(fileName, O_RDWR | O_CREAT, S_IWRITE | S_IREAD);
    if ((write(fd, sizeArr, sizeof(int) * numOfNodes)) < 0)
	abort();

    //check size arr
    intT numOfVert = 0;
    for (int i = 0; i < numOfNodes; i++) {
	numOfVert += sizeArr[i];
    }

    if (numOfVert != graph.n) {
	printf("size check wrong!!\n");
	abort();
    }

    long long totalSize = 0;
    for (intT i = 0; i < graph.n; i++) {
	totalSize += 2 * sizeof(intT);
    }

    void *buf = (void *)malloc(totalSize);
    char *ptr = (char *)buf;
    
    for (intT i = 0; i < graph.n; i++) {
	*(intT *)ptr = graph.V[i].getOutDegree();
	ptr += sizeof(intT);

	*(intT *)ptr = graph.V[i].getInDegree();
	ptr += sizeof(intT);
    }

    long long written = 0;
    while (written < totalSize) {
	long long sizeWritten = write(fd, (void *)((char *)buf + written), totalSize - written);
	if (sizeWritten < 0) {
	    printf("oops\n");
	}	    
	written += sizeWritten;
    }
}

template <class vertex>
graph<vertex> loadPartitionFromFile(char *fileName, intT *sizeArr, intT numOfNodes) {
    int fd = open(fileName, O_RDONLY, S_IREAD);
    read(fd, (void *)sizeArr, sizeof(intT) * numOfNodes);

    long long totalSize = 0;
    for (int i = 0; i < numOfNodes; i++) {
	totalSize += sizeArr[i];
    }

    void *buf = (void *)malloc(totalSize * sizeof(intT) * 2);
    printf("number of vert is: %ld %p\n", totalSize, buf);

    long long readSize = 0;
    while (readSize < totalSize) {
	long long readOnce = read(fd, (void *)((char *)buf + readSize), totalSize - readSize);
	readSize += readOnce;
	printf("read %ld, accum %ld\n", readOnce, readSize);
    }

    const intT n = totalSize;
    char *ptr = (char *)buf;

    vertex *vertices = newA(vertex, n);
    long long counter = 0;
    long long inCounter = 0;

    for (intT i = 0; i < n; i++) {
	intT outDeg = *(intT *)ptr;
	ptr += sizeof(intT);

	intT inDeg = *(intT *)ptr;
	ptr += sizeof(intT);
	
	vertices[i].setOutDegree(outDeg);
	vertices[i].setInDegree(inDeg);
    }
    
    free(buf);
    return graph<vertex>(vertices, (intT)n, 0, NULL);
}

template <class vertex>
void dumpSubgraphToEdgeList(graph<vertex> &graph, char *fileName) {
    const intT n = graph.n;
    
    long long totalSize = 0;
    for (intT i = 0; i < n; i++) {
	//totalSize += 4 * sizeof(intT); // const -1, vertID, outDeg, inDeg
	totalSize += graph.V[i].getFakeDegree() * sizeof(intE) * 2; //save both src and dst
    }
    totalSize += sizeof(long long); //totalSize itself
    
    printf("subgraph total size is: %ld\n", totalSize);
    
    void *buf = (void *)malloc(totalSize);
    
    char *ptr = (char *)buf;
    
    *(long long *)ptr = totalSize;
    ptr += sizeof(long long);       
    
    for (intT i = 0; i < n; i++) {
	for (intT j = 0; j < graph.V[i].getFakeDegree(); j++) {
	    *(intE *)ptr = i;
	    ptr += sizeof(intE);
	    *(intE *)ptr = graph.V[i].getOutNeighbor(j);
	    ptr += sizeof(intE);
	}
    }
    
    int fd = open(fileName, O_RDWR | O_CREAT, S_IWRITE | S_IREAD);
    long long written = 0;
    while (written < totalSize) {
	long long sizeWritten = write(fd, (void *)((char *)buf + written), totalSize - written);
	if (sizeWritten < 0) {
	    printf("oops\n");
	}	    
	written += sizeWritten;
    }
}
