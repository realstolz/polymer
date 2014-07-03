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
#include "parallel.h"
#include "gettime.h"

using namespace std;

bool needResult = false;

struct BFS_F {
  intT* Parents;
  BFS_F(intT* _Parents) : Parents(_Parents) {}
  inline bool update (intT s, intT d) { //Update
    if(Parents[d] == -1) { Parents[d] = s; return 1; }
    else return 0;
  }
  inline bool updateAtomic (intT s, intT d){ //atomic version of Update
    return (CAS(&Parents[d],(intT)-1,s));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (intT d) { return (Parents[d] == -1); } 
};

template <class vertex>
void BFS(intT start, graph<vertex> GA) {
  startTime();
  intT n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  intT* Parents = newA(intT,GA.n);
  parallel_for(intT i=0;i<GA.n;i++) Parents[i] = -1;
  Parents[start] = start;

  vertices Frontier(n,start); //creates initial frontier

  intT round = 0;
  intT numVisited = 0;
  struct timeval startT, endT;
  struct timezone tz = {0, 0};
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    round++;
    numVisited+=Frontier.numNonzeros();
    //printf("num of non zeros: %d\n", Frontier.numNonzeros());
    gettimeofday(&startT, &tz);
    //apply edgemap
    vertices output = edgeMap(GA, Frontier, BFS_F(Parents),GA.m/20);    
    gettimeofday(&endT, &tz);

    double timeStart = ((double)startT.tv_sec) + ((double)startT.tv_usec) / 1000000.0;
    double timeEnd = ((double)endT.tv_sec) + ((double)endT.tv_usec) / 1000000.0;

    double mapTime = timeEnd - timeStart;
    
    //printf("edge map time: %lf\n", mapTime);

    Frontier.del();
    Frontier = output; //set new frontier
  } 

  Frontier.del();

  cout<<"Vertices visited = "<<numVisited<<endl;
  // cout<<"Edges traversed = "<<edgesTraversed<<endl;
  cout << "Finished in " << round << " iterations\n";
  nextTime("BFS");

  if (needResult) {
      for (intT i = 0; i < GA.n; i++) {
	  cout << i << "\t" << std::scientific << std::setprecision(9) << Parents[i] << "\n";
      }
  }
  free(Parents);
}

int parallel_main(int argc, char* argv[]) {  
  char* iFile;
  bool binary = false;
  bool symmetric = false;
  int start = 0;
  if(argc > 1) iFile = argv[1];
  if(argc > 2) start = atoi(argv[2]);
  if(argc > 3) if((string) argv[3] == (string) "-result") needResult = true;
  //pass -s flag if graph is already symmetric
  if(argc > 4) if((string) argv[4] == (string) "-s") symmetric = true;
  //pass -b flag if using binary file (also need to pass 2nd arg for now)
  if(argc > 5) if((string) argv[5] == (string) "-b") binary = true;

  if(symmetric) {
    graph<symmetricVertex> G = 
      readGraph<symmetricVertex>(iFile,symmetric,binary); //symmetric graph
    BFS((intT)start,G);
    G.del(); 
  } else {
    graph<asymmetricVertex> G = 
      readGraph<asymmetricVertex>(iFile,symmetric,binary); //asymmetric graph
    BFS((intT)start,G);
    G.del();
  }
}
