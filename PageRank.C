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

template <class vertex>
struct PR_F {
  double* p_curr, *p_next;
  vertex* V;
  PR_F(double* _p_curr, double* _p_next, vertex* _V) : 
    p_curr(_p_curr), p_next(_p_next), V(_V) {}
  inline bool update(intT s, intT d){ //update function applies PageRank equation
    p_next[d] += p_curr[s]/V[s].getOutDegree();
    return 1;
  }
  inline bool updateAtomic (intT s, intT d) { //atomic Update
    writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
    /*
    if (d == 0) {
	cout << "Update from " << s << "\t" << std::scientific << std::setprecision(9) << p_curr[s]/V[s].getOutDegree() << " -- " << p_next[d] << "\n";
    }
    */
    return 1;
  }
  inline bool cond (intT d) { return 1; } //does nothing
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

template <class vertex>
void PageRank(graph<vertex> GA, int maxIter = -1) {
    //startTime();
  const intT n = GA.n;
  const double damping = 0.85;
  const double epsilon = 0.0000001;
  
  double one_over_n = 1/(double)n;
  double* p_curr = newA(double,n);
  double* p_next = newA(double,n);
  bool* frontier = newA(bool,n);
  
  double mapTime = 0.0;
  struct timeval start, end;
  struct timezone tz = {0, 0};

  {parallel_for(intT i=0;i<n;i++) p_curr[i] = one_over_n;}
  {parallel_for(intT i=0;i<n;i++) p_next[i] = 0;} //0 if unchanged
  {parallel_for(intT i=0;i<n;i++) frontier[i] = 1;}

  vertices Frontier(n,n,frontier);
  nextTime("Init");
  intT round = 0;
  while(1){
    if (maxIter > 0 && round >= maxIter)
      break;
    round++;

    //gettimeofday(&start, &tz);

    vertices output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V),GA.m/20,DENSE_FORWARD);

    //gettimeofday(&end, &tz);
    //double timeStart = ((double)start.tv_sec) + ((double)start.tv_usec) / 1000000.0;
    //double timeEnd = ((double)end.tv_sec) + ((double)end.tv_usec) / 1000000.0;
    //mapTime = timeEnd - timeStart;

    vertexMap(Frontier,PR_Vertex_F(p_curr,p_next,damping,n));
    /*
    double tmpConst = (1-damping)*(1/(double)n);
    
    {parallel_for (intT i = 0; i < n; i++) {
	    p_next[i] = damping * p_next[i] + tmpConst;
	}}
    */
    /*
    //compute L1-norm between p_curr and p_next
    {parallel_for(intT i=0;i<n;i++) {
      p_curr[i] = fabs(p_curr[i]-p_next[i]);
      }}
    double L1_norm = sequence::plusReduce(p_curr,n);
    //cout<<"Round "<<round<<", L1 norm = "<<L1_norm<<endl;
    if(L1_norm < epsilon) break;
    */
    /*
    {parallel_for (intT i = 0; i < n; i++) {
	    p_curr[i] = p_next[i];
	    p_next[i] = 0.0;
	}}	   
    */
    //reset p_curr
    vertexMap(Frontier,PR_Vertex_Reset(p_curr));
    swap(p_curr,p_next);
    
    Frontier.del(); 
    Frontier = output;
    //printf("round %d: %lf\n", round, mapTime);
  }
  cout<<"Finished in "<<round<<" iterations\n";
  Frontier.del();
  nextTime("PageRank");

  if (needResult) {
      for (intT i = 0; i < GA.n; i++) {
	  cout << i << "\t" << std::scientific << std::setprecision(9)<< p_curr[i] << "\n";
      }
  }

  free(p_curr); free(p_next); 
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
  startTime();
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
}
