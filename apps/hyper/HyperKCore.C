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

// Parallel implementation of K-Core decomposition of a symmetric
// hypergraph.
#define HYPER 1
#include "hygra.h"

template <class vertex>
struct Init_Deg {
  long* Degrees;
  vertex* H;
  Init_Deg(long* _Degrees, vertex* _H) : Degrees(_Degrees), H(_H) {}
  inline bool update (uintE s, uintE d) { 
    xadd(&Degrees[s],(long)H[d].getOutDegree()-1);
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    Degrees[s]+=H[d].getOutDegree()-1;
    return 1;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

struct Update_Deg {
  long* Degrees;
  intE *Counts;
  Update_Deg(long* _Degrees, intE* _Counts) : Degrees(_Degrees), Counts(_Counts) {}
  inline bool update (uintE s, uintE d) {
    Degrees[d] -= Counts[s];
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    writeAdd(&Degrees[d],(long)-Counts[s]);
    //xadd(&Degrees[d],(long)-Counts[s]);
    return 1;
  }
  inline bool cond (uintE d) { return Degrees[d] > 0; }
};

//return true for neighbor the first time it's updated 
struct Count_Removed {
  intE* Counts;
  Count_Removed(intE* _Counts) : Counts(_Counts) {}
  inline bool update (uintE s, uintE d) {
    Counts[d]++;
    return Counts[d] == 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    volatile intE oldV, newV; 
    do { 
      oldV = Counts[d]; newV = oldV + 1;
    } while(!CAS(&Counts[d],oldV,newV));
    return oldV == 0.0;

    //return xadd(&Counts[d],1) == 0;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

template<class vertex>
struct Deg_LessThan_K {
  vertex* V;
  uintE* coreNumbers;
  long* Degrees;
  uintE k;
  Deg_LessThan_K(vertex* _V, long* _Degrees, uintE* _coreNumbers, uintE _k) : 
    V(_V), k(_k), Degrees(_Degrees), coreNumbers(_coreNumbers) {}
  inline bool operator () (uintE i) {
    if(Degrees[i] < k) { coreNumbers[i] = k-1; Degrees[i] = 0; return true; }
    else return false;
  }
};

template<class vertex>
struct Deg_AtLeast_K {
  vertex* V;
  long *Degrees;
  uintE k;
  Deg_AtLeast_K(vertex* _V, long* _Degrees, uintE _k) : 
    V(_V), k(_k), Degrees(_Degrees) {}
  inline bool operator () (uintE i) {
    return Degrees[i] >= k;
  }
};

//assumes symmetric hypergraph
// 1) iterate over all remaining active vertices
// 2) for each active vertex, remove if induced degree < k. Any vertex removed has
//    core-number (k-1) (part of (k-1)-core, but not k-core)
// 3) stop once no vertices are removed. Vertices remaining are in the k-core.
template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
  const long nv = GA.nv, nh = GA.nh;
  bool* active = newA(bool,nv);
  {parallel_for(long i=0;i<nv;i++) active[i] = 1;}
  vertexSubset Frontier(nv, nv, active);
  uintE* coreNumbers = newA(uintE,nv);
  long* Degrees = newA(long,nv);
  {parallel_for(long i=0;i<nv;i++) {
      coreNumbers[i] = 0;
      Degrees[i] = 0;
    }}
  edgeMap(GA,FROM_V,Frontier,Init_Deg<vertex>(Degrees,GA.H),INT_T_MAX,no_output);
  intE* Counts = newA(intE,nh);
  {parallel_for(long i=0;i<nh;i++) Counts[i] = 0;}
  long largestCore = -1;
  for (long k = 1;;k++) {
    while (true) {
      vertexSubset toRemove 
	= vertexFilter(Frontier,Deg_LessThan_K<vertex>(GA.V,Degrees,coreNumbers,k));
      //cout << "k = " << k <<  " to remove " << toRemove.numNonzeros() << " frontier size = " << Frontier.numNonzeros() << endl;
      vertexSubset remaining = vertexFilter(Frontier,Deg_AtLeast_K<vertex>(GA.V,Degrees,k));
      Frontier.del();
      Frontier = remaining;
      //cout << "num remaining vertices = " << remaining.numNonzeros() << endl;
      //long maxDegree = sequence::reduce(Degrees,nv,maxF<long>());
      //cout << "max remaining degree = " << maxDegree << endl;
      if (0 == toRemove.numNonzeros()) { // fixed point. found k-core
	toRemove.del();
        break;
      }
      else {
	vertexSubset FrontierH = edgeMap(GA,FROM_V,toRemove,Count_Removed(Counts));
	cout << "k="<<k-1<< " num active = " << toRemove.numNonzeros() << " frontierH = " << FrontierH.numNonzeros() << endl;
	edgeMap(GA,FROM_H,FrontierH,Update_Deg(Degrees,Counts),-1,no_output);
	auto reset_counts = [&] (const uintE& i) { Counts[i] = 0; };
	vertexMap(FrontierH,reset_counts);
	FrontierH.del();
	toRemove.del();
      }
    }
    if(Frontier.numNonzeros() == 0) { largestCore = k-1; break; }
  }
  cout << "largestCore was " << largestCore << endl;
  Frontier.del(); free(coreNumbers); free(Degrees); free(Counts);
}
