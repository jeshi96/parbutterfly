#ifndef _BUTILS_
#define _BUTILS_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "hypergraphIO.h"
#include "parseCommandLine.h"
#include "vertex.h"
#include "sequence.h"

using namespace std;

// Takes the elements of a vertex array, and returns the out degree choose 2
template <class vertex, class E>
struct chooseV { 
  vertex* V;
  chooseV(vertex* _V) : V(_V) {}
  E operator() (const E& i) const {
    uintE v_deg = V[i].getOutDegree();
    return (E) ((v_deg * (v_deg-1)) / 2); 
  }
};

template <class vertex, class E>
struct seagullSumHelper { 
  vertex* V;
  vertex u;
  seagullSumHelper(vertex* _V,vertex _u) : V(_V),u(_u) {}
  E operator() (const E& i) const {
	return (E) (V[u.getOutNeighbor(i)].getOutDegree() - 1);
  }
};

template <class vertex, class E>
struct seagullSum { 
  vertex* V;
  vertex* U;
  uintE* active;
  seagullSum(vertex* _V,vertex* _U, uintE* _active) : V(_V),U(_U),active(_active) {}
  E operator() (const E& i) const {
	uintE u_idx = active[i];
	return sequence::reduce<E>((E) 0, (long) U[u_idx].getOutDegree(), addF<E>(),seagullSumHelper<vertex,E>(V,U[u_idx]));
  }
};

struct WedgeIntCons {
  long nu;
  WedgeIntCons(long _nu) : nu(_nu) {}
  tuple<uintE,uintE> operator() (uintE v1, uintE v2, uintE c) {
    return make_tuple(v1 <= v2 ? v1 * nu + v2 : v2 * nu + v1, c);
  }
};

struct Wedge {
  uintE v1;
  uintE v2;
  uintE u;
  Wedge(uintE _v1, uintE _v2, uintE _u) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1), u(_u) {}
};

struct WedgeCons { Wedge operator() (uintE v1, uintE v2, uintE c) { return Wedge(v1, v2, c); }};

struct WedgeCmp {
  long nv;
  WedgeCmp(long _nv) : nv(_nv) {}
  bool operator() (Wedge vs1, Wedge vs2) {
    return vs1.v1 * nv + vs1.v2 < vs2.v1 * nv + vs2.v2;
  }
};

struct WedgeEq { bool operator() (Wedge vs1, Wedge vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Represents a pair of vertices on one side of a bipartite graph (ordered, with least vertex first)
struct VertexPair {
  uintE v1;
  uintE v2;
  VertexPair(uintE _v1, uintE _v2) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1) {}
};

// Represents a pair of vertices on one side of a bipartite graph (unordered, stored based on constructor order)
struct UVertexPair {
  uintE v1;
  uintE v2;
  UVertexPair(uintE _v1, uintE _v2) : v1(_v1), v2(_v2) {}
};

// Comparer for VertexPair based on least vertex in pair and then greatest vertex in pair
struct VertexPairCmp {
  long nv;
  VertexPairCmp(long _nv) : nv(_nv) {}
  bool operator() (VertexPair vs1, VertexPair vs2) {
    return vs1.v1 * nv + vs1.v2 < vs2.v1 * nv + vs2.v2;
  }
};

// Comparer for VertexPair based on greatest vertex in pair and then least vertex in pair
struct VertexPairCmp2 {
  long nv;
  VertexPairCmp2(long _nv) : nv(_nv) {}
  bool operator() (VertexPair vs1, VertexPair vs2) {
    return vs1.v2 * nv + vs1.v1 < vs2.v2 * nv + vs2.v1;
  }
};

// Comparer for UVertexPair
struct UVertexPairCmp {
  long nv;
  UVertexPairCmp(long _nv) : nv(_nv) {}
  bool operator() (UVertexPair vs1, UVertexPair vs2) {
    return vs1.v2 * nv + vs1.v1 < vs2.v2 * nv + vs2.v1;
  }
};

// Equality for VertexPair and UVertexPair
struct VertexPairEq { bool operator() (VertexPair vs1, VertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };
struct UVertexPairEq { bool operator() (UVertexPair vs1, UVertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Constructs a VertexPair and UVertexPair
struct VertexPairCons { VertexPair operator() (uintE v1, uintE v2, uintE c) { return VertexPair(v1, v2); }};
struct UVertexPairCons { UVertexPair operator() (uintE v1, uintE v2) { return UVertexPair(v1, v2); }};

// Constructs a uintE form of a VertexPair and UVertexPair
struct VertexPairIntCons {
  long nu;
  VertexPairIntCons(long _nu) : nu(_nu) {}
  uintE operator() (uintE v1, uintE v2, uintE c) {
    return v1 <= v2 ? v1 * nu + v2 : v2 * nu + v1;
  }
};
struct UVertexPairIntCons {
  long nu;
  UVertexPairIntCons(long _nu) : nu(_nu) {}
  uintE operator() (uintE v1, uintE v2) {
    return v1 * nu + v2;
  }
};

// Comparer for indices for VertexPairs in nest, by v1 or v2 (based on use_v1)
struct NestedVPCmp {
  VertexPair* nest;
  bool use_v1;
  NestedVPCmp(VertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 < nest[idx2].v1;
    return nest[idx1].v2 < nest[idx2].v2;
  }
};
struct NestedVPEq {
  VertexPair* nest;
  bool use_v1;
  NestedVPEq(VertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 == nest[idx2].v1;
    return nest[idx1].v2 == nest[idx2].v2;
  }
};

struct nonZeroF{bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != 0);}};
struct greaterOneF{bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) > 1);}};

struct nonMaxTupleF{bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != UINT_E_MAX || get<0>(a) != UINT_E_MAX);}};

struct uintELt {bool operator () (uintE a, uintE b) {return a < b;};};
struct uintEEq {bool operator() (uintE a, uintE b) {return a == b;};};
struct uintETupleLt {bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<0>(a) < get<0>(b);} };
struct uintETupleEq {bool operator() (tuple<uintE,uintE> a,tuple<uintE,uintE> b) {return get<0>(a) == get<0>(b); }};
struct uintETupleAdd {
  tuple<uintE,uintE> operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) const {
    return make_tuple(get<0>(a), get<1>(a) + get<1>(b));
  };
};

template<class T> struct refl{T operator() (T obj) {return obj;}};
template<class T> struct reflCount{long operator() (T obj) {return 1;}};
struct UVertexPairV2{uintE operator() (UVertexPair obj) {return obj.v2;}};
struct choose2{uintE operator() (uintE obj) {return obj*(obj-1)/2;}};
struct uintETupleGet0{uintE operator() (tuple<uintE,uintE> obj) {return get<0>(obj);}};
struct uintECount{long operator() (tuple<uintE,uintE> obj) {return get<1>(obj);}};

template <class E>
struct writeAddArr {
  E* arr;
  writeAddArr(E* _arr) : arr(_arr) {}
  void operator() (long idx, E num) {
    writeAdd(&arr[idx], num);
  }
};

template <class E>
struct writeAddSet {
  sparseAdditiveSet<E> set;
  writeAddSet(sparseAdditiveSet<E> _set) : set(_set) {}
  void operator() (long idx, E num) {
    set.insert(pair<uintE, E>(idx, num));
  }
};

template<class E, class K>
E getAdd (E curr, tuple<K,E> v) {
  return curr + get<1>(v);
}

template<class E, class K>
tuple<K,E> getAddReduce (tuple<K,E> curr, tuple<K,E> v) {
  return make_tuple(get<0>(curr),get<1>(curr) + get<1>(v));
}


//***********************************************************************************************
//***********************************************************************************************
// Graph construction + reading

/*
 *  Constructs a complete bipartite graph, with nv vertices on one side and nu vertices on
 *  the otheer.
 * 
 *  nv: Number of vertices in one bipartition (V)
 *  nu: Number of vertices in the other bipartition (U)
 * 
 *  Returns: Complete bipartite graph
 */
template<class vertex>
bipartiteGraph<vertex> bpGraphComplete(long nv, long nu){
  // Construct vertex partitions
  vertex* v = newA(vertex,nv+1);
  vertex* u = newA(vertex,nu+1);
  // Add all neighbors to v
  parallel_for(int i=0;i<nv;++i) {
    uintE* neighbors_v = newA(uintE,nu);
    parallel_for(int j=0;j<nu;++j) {
      neighbors_v[j] = j;
    }
    v[i] = vertex(neighbors_v,nu);
  }
  // Add all neighbors to u
  parallel_for(int i=0;i<nu;++i) {
    uintE* neighbors_u = newA(uintE,nv);
    parallel_for(int j=0;j<nv;++j) {
      neighbors_u[j] = j;
    }
    u[i] = vertex(neighbors_u,nv);
  }

  Uncompressed_Membipartitegraph<vertex>* mem = 
    new Uncompressed_Membipartitegraph<vertex>(v,u,nv,nu);
  return bipartiteGraph<vertex>(v,u,nv,nu,mem);
}

/*
 *  Constructs a bipartite graph from a hypergraph (by taking the vertices of the
 *  hypergraph as one vertex partition, and the edges of the hypergraph as the
 *  other vertex partition; edges of the bipartite graph are given by inclusion of 
 *  vertices in edges of the hypergraph).
 *  The vertex class is assumed to be symmetricVertex.
 * 
 *  G: Hypergraph
 * 
 *  Returns: Bipartite graph from G
 */
template<class vertex>
bipartiteGraph<symmetricVertex> bpGraphFromHypergraph(hypergraph<vertex> G){
  long nv = G.nv;
  long nu = G.nh;
  symmetricVertex* v = newA(symmetricVertex,nv);
  symmetricVertex* u = newA(symmetricVertex,nu);

  // Copy all neighbors from inclusion in hypergraph
  parallel_for(int i=0;i<nv;++i) {
    uintE* neighbors_v = newA(uintE,G.V[i].getOutDegree());
    parallel_for(int j=0;j<G.V[i].getOutDegree();++j) {
      neighbors_v[j] = G.V[i].getOutNeighbor(j);
    }
    v[i] = symmetricVertex(neighbors_v, G.V[i].getOutDegree());
  }

  parallel_for(int i=0;i<nu;++i) {
    uintE* neighbors_u = newA(uintE,G.H[i].getInDegree());
    parallel_for(int j=0;j<G.H[i].getInDegree();++j) {
      neighbors_u[j] = G.H[i].getInNeighbor(j);
    }
    u[i] = symmetricVertex(neighbors_u, G.H[i].getInDegree());
  }

  Uncompressed_Membipartitegraph<symmetricVertex>* mem = 
    new Uncompressed_Membipartitegraph<symmetricVertex>(v,u,nv,nu);
  return bipartiteGraph<symmetricVertex>(v,u,nv,nu,mem);

}

//***********************************************************************************************
//***********************************************************************************************

/*
 *  Sort objects using cmp, and then retrieve indices of where consecutive objects are
 *  not equal (as given by eq) to give frequency counts of each object in objs.
 *  In other words, if we let arr be the returned array, arr[i+1] - arr[i] is the 
 *  frequency of objs[arr[i]] (objs[arr[i]] = ... = objs[arr[i+1]-1]).
 *  Iterating from 0 (inclusive) to the returned length - 1 (exclusive) will give
 *  all frequency counts.
 * 
 *  objs: Objects to count frequencies of
 *  num : Length of objs
 *  cmp : Comparator for T objects
 *  eq  : Equality comparator for T objects
 *  sort: If objs is already sorted, then this should be set to false so we don't resort;
 *        by default, we assume objs is not sorted.
 * 
 *  Returns: Array and length of array with frequency counts (as described above)
 */
template <class T, class Cmp, class Eq>
pair<uintE*, long> getFreqs(T* objs, long num, Cmp cmp, Eq eq, bool sort=true) {
  // Sort objects
  if (sort) sampleSort(objs, num, cmp);

  uintE* freqs = newA(uintE, num + 1);
  freqs[0] = 0;
  freqs[num] = num;

  // Retrieve indices where objects differ
  parallel_for(long i=1; i < num; ++i) {
    if (!eq(objs[i-1],objs[i])) freqs[i] = i;
    else freqs[i] = UINT_E_MAX;
  }
  uintE* freqs_f = newA(uintE, num+1);
  long num_freqs_f = sequence::filter(freqs, freqs_f, num+1, nonMaxF());
  free(freqs);
  return make_pair(freqs_f, num_freqs_f);
}

template <class S, class T, class Cmp, class Eq, class OpT, class OpuintE, class OpCount>
pair<tuple<S,uintE>*, long> getFreqs_seq(T* objs, long num, Cmp cmp, Eq eq, bool sort=true, OpT opt=refl<T>(),
  OpuintE opuinte=refl<uintE>(), OpCount opcount=reflCount<T>()) {
  if(sort) sampleSort(objs, num, cmp);

  using X = tuple<S,uintE>;
  X* freqs = newA(X, num);
  T prev = objs[0];
  T curr = objs[0];
  long idx = 0;
  long count = opcount(prev);
  for(long i=1; i < num; ++i) {
    curr = objs[i];
    if (!eq(prev, curr)) {
      freqs[idx] = make_tuple(opt(prev), opuinte(count));
      idx++;
      count = opcount(curr);
      prev = curr;
    }
    else {
      count += opcount(curr);
    }
  }
  freqs[idx] = make_tuple(opt(curr), opuinte(count));
  return make_pair(freqs, idx + 1);
}

//********************************************************************************************
//********************************************************************************************

/*
 *  Computes the total number of wedges on all vertices on one bipartition of our graph, as 
 *  specified by V (note that the vertices of V are considered the centers of the wedges).
 * 
 *  nv: Number of vertices in our bipartition V
 *  V : Bipartition of vertices that forms the center of our wedges
 * 
 *  Returns: Number of total wedges
 */
template<class vertex>
long countWedges(long nv, vertex* V) {
  return sequence::reduce<long>((long) 0, nv, addF<long>(), chooseV<vertex, long>(V));
}

/*
 *  Compares the total number of wedges on all vertices on one bipartition of our graph
 *  to the other bipartition. Returns the side with the least number of wedges, and
 *  the total number of wedges on that side.
 * 
 *  GA: Bipartite graph
 * 
 *  Returns: True if GA.V is the side with the least number of wedges total. False
 *           otherwise. Also, returns the total number of wedges on the corresponding
 *           side.
 */
template <class vertex>
pair<bool,long> cmpWedgeCounts(bipartiteGraph<vertex> GA) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = countWedges<vertex>(nv,GA.V);
  long num_wedges_u = countWedges<vertex>(nu, GA.U);
  return make_pair((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
}

/*
 *  Retrieves all wedges on the bipartition V of our graph, where the vertices V are 
 *  taken to be the center of the wedges. Stores these wedges using the constructor 
 *  cons on the two endpoints that make up the wedge.
 * 
 *  nv  : Number of vertices in our bipartition V
 *  V   : Bipartition of vertices that forms the center of our wedges
 *  cons: Constructor to create wedge objects, given two endpoints of a wedge + center
 * 
 *  Returns: Array of wedges on bipartition V
 */
template<class wedge, class vertex, class wedgeCons>
wedge* getWedges(const long nv, const vertex* V, wedgeCons cons) {
//timer t;
//t.start();
  // Retrieve the indices of each wedge associated with each vertex in V
  long* wedge_idxs = newA(long,nv+1);
  parallel_for(long i=0;i<nv+1;++i){
    wedge_idxs[i] = 0;
  }
  parallel_for (long i = 0; i < nv; ++i) {
    uintE v_deg = V[i].getOutDegree();
    if (v_deg >= 2)
      wedge_idxs[i] = v_deg * (v_deg - 1)/2;
  }
  sequence::plusScan(wedge_idxs, wedge_idxs, nv+1);
  long num_wedges = wedge_idxs[nv];

  // Retrieve each wedge associated with each vertex in V
  wedge* wedges = newA(wedge,num_wedges);
  parallel_for (long i = 0; i < nv; ++i) {
    const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    long wedge_idx = wedge_idxs[i];
    long idx = 0;
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
        wedges[wedge_idx + idx] = cons(v.getOutNeighbor(j), v.getOutNeighbor(k), i);
        ++idx;
      }
    }
  }
  free(wedge_idxs);
//t.stop();
//t.reportTotal("\tgetWedges:");

  return wedges;
}


template<class vertex, class wedgeCons>
void getWedgesHash(sparseAdditiveSet<uintE>& wedges, long nv, vertex* V, wedgeCons cons) {
//timer t;
//t.start();
// Count number of wedges by their key
  parallel_for (long i = 0; i < nv; ++i) {
    const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
        wedges.insert(pair<uintE,uintE>(cons(v.getOutNeighbor(j), v.getOutNeighbor(k), i), 1));
      }
    }
  }
//t.stop();
//t.reportTotal("\tgetWedgesHash:");
}

#endif