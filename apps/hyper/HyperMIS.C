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
//
// This is an implementation of the MIS algorithm on hyper graphs from
// "Parallel Search for Maximal Independence Given Minimal
// Dependence", Proceedings of the ACM-SIAM Symposium on Discrete
// Algorithms (SODA), 1990 by Paul Beame and Michael Luby. We choose a
// different sampling probability for better performance.
#define HYPER 1
#include "hygra.h"

#define CHECK 1

struct MIS_Check {
  uintT* flags;
  intT* Degrees;
  MIS_Check(uintT* _flags, intT* _Degrees) : flags(_flags), Degrees(_Degrees) {}
  inline bool update (uintE s, uintE d) {
    if(flags[d]>1) xadd(&Degrees[s],1);
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) {
    if(flags[d]>1) Degrees[s]++;
  }
  inline bool cond (uintE i) {return cond_true(i);}
};

struct MIS_Reset {
  uintT* flags;
  uintT round;
  MIS_Reset(uintT* _flags, uintT _round) : flags(_flags), round(_round) {}
  inline bool update (uintE s, uintE d) {
    flags[d] = 0;
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) {
    CAS(&flags[d],round,(uintT)0);
    return 1;
  }
  inline bool cond (uintE i) {return flags[i] == round;}
};

struct MIS_Select {
  uintT* flags;
  long offset, inverseProb, round;
  MIS_Select(uintT* _flags, long _round, long _offset, long _inverseProb) :
    flags(_flags), round(_round), offset(_offset), inverseProb(_inverseProb) {}
  inline void operator () (uintE i) {
    if(hashInt((ulong)(i+offset)) % inverseProb == 0) flags[i] = round;    
  }
};

struct MIS_Filter {
  uintT* flags;
  MIS_Filter(uintT* _flags) : flags(_flags) {}
  inline bool operator () (uintE i) {
    return flags[i] == 0;
  }
};

template <class vertex>
struct Edge_Filter {
  vertex* H;
  intT* Degree;
  Edge_Filter(vertex* _H, intT* _Degree) : H(_H), Degree(_Degree) {}
  inline bool operator () (uintE i) {
    return Degree[i] == H[i].getOutDegree();
  }
};

template <class vertex>
struct Check_H {
  vertex* H;
  uintT* flags;
  Check_H(vertex* _H, uintT* _flags) : H(_H), flags(_flags) {}
  inline bool operator () (uintE i) {
    if(H[i].getOutDegree() == 0) return 0;
    else if(H[i].getOutDegree() == 1) {
      if(flags[H[i].getOutNeighbor(0)] == 0) flags[H[i].getOutNeighbor(0)] = 1;
      return 0;
    } else return 1;
  }
};

struct Degree_Reset {
  intT* Degree;
  Degree_Reset(intT* _Degree) : Degree(_Degree) {}
  inline void operator () (uintE i) { Degree[i] = 0; }
};

//Takes a symmetric graph as input; priority of a vertex is its ID.
template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
  const intE nv = GA.nv, nh = GA.nh;
  uintT* flags = newA(uintT,nv); //undecided = 0, out = 1, in = anything else
  bool* frontier_data = newA(bool, nv);
  {parallel_for(long i=0;i<nv;i++) {
    flags[i] = 0;
    frontier_data[i] = 1;
  }}
  bool* frontierH = newA(bool, nh);
  intT* Degree = newA(intT,nh);
  {parallel_for(long i=0;i<nh;i++) {frontierH[i] = 1; }}
  long round = 1;
  long inverseProb = 2; //to be set
  vertexSubset FrontierV(nv, frontier_data);
  vertexSubset FrontierH(nh, frontierH);

  while (!FrontierV.isEmpty()) {
    round++;
    vertexMap(FrontierV,MIS_Select(flags,round,round*nv,inverseProb));
    vertexMap(FrontierH,Degree_Reset(Degree));
    cout << round << " " << FrontierV.numNonzeros() << " " << FrontierH.numNonzeros() << endl;

    edgeMap(GA, FROM_H, FrontierH, MIS_Check(flags,Degree), -1, no_output);
    vertexSubset fullEdges = vertexFilter(FrontierH, Edge_Filter<vertex>(GA.H,Degree));
    cout << "full edges = " << fullEdges.numNonzeros() << endl;

    edgeMap(GA, FROM_H, fullEdges, MIS_Reset(flags,round), -1, no_output);
    fullEdges.del();
    //pack edges
    auto pack_predicate = [&] (const uintE& u, const uintE& ngh) { return flags[ngh] == 0; };
    packEdges(GA.H, FrontierH, pack_predicate, no_output);

    vertexSubset remainingHyperedges = vertexFilter(FrontierH, Check_H<vertex>(GA.H,flags));

    FrontierH.del();
    FrontierH = remainingHyperedges;
    vertexSubset output = vertexFilter(FrontierV, MIS_Filter(flags));
    FrontierV.del();
    FrontierV = output;
  }

#ifdef CHECK
  {for(long i=0;i<nh;i++) {
      vertex h = GA.H[i];
      long inSet = 0;
      for(long j=0;j<h.getInDegree();j++) if(flags[h.getInNeighbor(j)] > 1) inSet++;
      if(inSet == h.getInDegree()) {cout << "incorrect answer " << inSet << " " << h.getInDegree() << " " << h.getOutDegree() << endl; exit(0);}
    }}
#endif
  
  free(flags); free(Degree);
  FrontierV.del(); FrontierH.del();
}
