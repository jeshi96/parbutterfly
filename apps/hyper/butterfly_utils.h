#ifndef _BUTILS_
#define _BUTILS_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "hypergraphIO.h"
#include "parseCommandLine.h"
#include "vertex.h"
#include "sequence.h"
#include "binary_search.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }


using namespace std;

// Takes the elements of a vertex array, and returns the out degree choose 2
template <class vertex, class E>
struct chooseV { 
  vertex* V;
  chooseV(vertex* _V) : V(_V) {}
  inline E operator() (const E& i) const {
    uintE v_deg = V[i].getOutDegree();
    return (E) ((v_deg * (v_deg-1)) / 2); 
  }
};

template <class vertex, class E>
struct getV { 
  vertex* V;
  getV(vertex* _V) : V(_V) {}
  inline E operator() (const E& i) const {
    return (E) V[i].getOutDegree();
  }
};

template <class vertex, class E>
struct seagullSumHelper { 
  vertex* V;
  vertex u;
  seagullSumHelper(vertex* _V,vertex _u) : V(_V),u(_u) {}
  inline E operator() (const E& i) const {
	return (E) (V[u.getOutNeighbor(i)].getOutDegree() - 1);
  }
};

template <class vertex, class E>
struct seagullSum { 
  vertex* V;
  vertex* U;
  uintE* active;
  seagullSum(vertex* _V,vertex* _U, uintE* _active) : V(_V),U(_U),active(_active) {}
  inline E operator() (const E& i) const {
	const vertex u = U[active[i]];
  E ret=0;
  for (long k=0; k < u.getOutDegree(); ++k) {
    ret += V[u.getOutNeighbor(k)].getOutDegree() - 1;
  }
  return ret;
	//return sequence::reduce<E>((E) 0, (long) U[u_idx].getOutDegree(), addF<E>(),seagullSumHelper<vertex,E>(V,U[u_idx]));
  }
};

struct WedgeIntCons {
  long nu;
  WedgeIntCons(long _nu) : nu(_nu) {}
  inline tuple<uintE,uintE> operator() (uintE v1, uintE v2, uintE c) {
    return make_tuple(v1 <= v2 ? v1 * nu + v2 : v2 * nu + v1, c);
  }
};

struct UWedgeIntCons {
  long nu;
  UWedgeIntCons(long _nu) : nu(_nu) {}
  inline tuple<uintE,uintE> operator() (uintE v1, uintE v2, uintE c) {
    return make_tuple(v1 * nu + v2, c);
  }
};

struct UWedge {
  uintE v1;
  uintE v2;
  uintE u;
  UWedge(uintE _v1, uintE _v2, uintE _u) : v1(_v1), v2(_v2), u(_u) {}
};

struct UWedgeCons { inline UWedge operator() (uintE v1, uintE v2, uintE c) { return UWedge(v1, v2, c); }};

struct UWedgeCmp {
  inline bool operator() (UWedge vs1, UWedge vs2) {
  	if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
  	return vs1.v1 < vs2.v1;
  }
};

struct UWedgeEq { inline bool operator() (UWedge vs1, UWedge vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

struct Wedge {
  uintE v1;
  uintE v2;
  uintE u;
  Wedge(uintE _v1, uintE _v2, uintE _u) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1), u(_u) {}
};

struct WedgeCons { inline Wedge operator() (uintE v1, uintE v2, uintE c) { return Wedge(v1, v2, c); }};

struct WedgeCmp {
  long nv;
  WedgeCmp(long _nv) : nv(_nv) {}
  inline bool operator() (Wedge vs1, Wedge vs2) {
    return vs1.v1 * nv + vs1.v2 < vs2.v1 * nv + vs2.v2;
  }
};

struct WedgeEq { inline bool operator() (Wedge vs1, Wedge vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Represents a pair of vertices on one side of a bipartite graph (ordered, with least vertex first)
struct VertexPair {
  uintE v1;
  uintE v2;
  VertexPair(uintE _v1, uintE _v2) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1) {}
};

struct VertexPairHash {
  uintE v1;
  uintE v2;
  uintE hash;
  VertexPairHash(uintE _v1, uintE _v2, uintE nu) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1) {
    hash = v1 * nu + v2;
  }
  VertexPairHash() : v1(0), v2(0), hash(0) {};
};

// Represents a pair of vertices on one side of a bipartite graph (unordered, stored based on constructor order)
struct UVertexPair {
  uintE v1;
  uintE v2;
  UVertexPair(uintE _v1, uintE _v2) : v1(_v1), v2(_v2) {}
};

struct UVertexPairHash {
  uintE v1;
  uintE v2;
  uintE hash;
  UVertexPairHash(uintE _v1, uintE _v2, uintE nu) : v1(_v1), v2(_v2) { hash = v1 * nu + v2; }
  UVertexPairHash() : v1(0), v2(0), hash(0) {};
};

//TODO get rid of legacy _nv
// Comparer for VertexPair based on least vertex in pair and then greatest vertex in pair
struct VertexPairCmp {
  long nv;
  VertexPairCmp(long _nv) : nv(_nv) {}
  inline bool operator() (VertexPair vs1, VertexPair vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

struct VertexPairHashCmp {
  inline bool operator() (VertexPairHash vs1, VertexPairHash vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

struct UVertexPairHashCmp {
  inline bool operator() (UVertexPairHash vs1, UVertexPairHash vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

// Comparer for VertexPair based on greatest vertex in pair and then least vertex in pair
struct VertexPairCmp2 {
  long nv;
  VertexPairCmp2(long _nv) : nv(_nv) {}
  inline bool operator() (VertexPair vs1, VertexPair vs2) {
    if (vs1.v2 == vs2.v2) return vs1.v1 < vs2.v1;
    return vs1.v2 < vs2.v2;
  }
};

// Comparer for UVertexPair
struct UVertexPairCmp2{
  long nv;
  UVertexPairCmp2(long _nv) : nv(_nv) {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v2 == vs2.v2) return vs1.v1 < vs2.v1;
    return vs1.v2 < vs2.v2;
  }
};

struct UVertexPairCmp {
  long nv;
  UVertexPairCmp(long _nv) : nv(_nv) {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

// Equality for VertexPair and UVertexPair
struct VertexPairEq { inline bool operator() (VertexPair vs1, VertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };
struct UVertexPairEq { inline bool operator() (UVertexPair vs1, UVertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };
struct VertexPairHashEq { inline bool operator() (VertexPairHash vs1, VertexPairHash vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };
struct UVertexPairHashEq { inline bool operator() (UVertexPairHash vs1, UVertexPairHash vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Constructs a VertexPair and UVertexPair
struct VertexPairCons { inline VertexPair operator() (uintE v1, uintE v2, uintE c) { return VertexPair(v1, v2); }};
struct UVertexPairCons { inline UVertexPair operator() (uintE v1, uintE v2, uintE c) { return UVertexPair(v1, v2); }};
struct VertexPairHashCons {
  long nu;
  VertexPairHashCons(long _nu) : nu(_nu) {}
  inline VertexPairHash* operator() (uintE v1, uintE v2, uintE c) {
    return new VertexPairHash(v1, v2, nu);
  }
};
struct VertexPairHashConsN {
  long nu;
  VertexPairHashConsN(long _nu) : nu(_nu) {}
  inline VertexPairHash operator() (uintE v1, uintE v2, uintE c) {
    return VertexPairHash(v1, v2, nu);
  }
};
struct UVertexPairHashCons {
  long nu;
  UVertexPairHashCons(long _nu) : nu(_nu) {}
  inline UVertexPairHash* operator() (uintE v1, uintE v2, uintE c) {
    return new UVertexPairHash(v1, v2, nu);
  }
};
struct UVertexPairHashConsN {
  long nu;
  UVertexPairHashConsN(long _nu) : nu(_nu) {}
  inline UVertexPairHash operator() (uintE v1, uintE v2, uintE c) {
    return UVertexPairHash(v1, v2, nu);
  }
};

// Constructs a uintE form of a VertexPair and UVertexPair
struct VertexPairIntCons {
  long nu;
  VertexPairIntCons(long _nu) : nu(_nu) {}
  inline uintE operator() (uintE v1, uintE v2, uintE c) {
    return v1 <= v2 ? v1 * nu + v2 : v2 * nu + v1;
  }
};
struct UVertexPairIntCons {
  long nu;
  UVertexPairIntCons(long _nu) : nu(_nu) {}
  inline uintE operator() (uintE v1, uintE v2, uintE c) {
    return v1 * nu + v2;
  }
};

// Comparer for indices for VertexPairs in nest, by v1 or v2 (based on use_v1)
struct NestedVPCmp {
  VertexPair* nest;
  bool use_v1;
  NestedVPCmp(VertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 < nest[idx2].v1;
    return nest[idx1].v2 < nest[idx2].v2;
  }
};
struct NestedVPEq {
  VertexPair* nest;
  bool use_v1;
  NestedVPEq(VertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 == nest[idx2].v1;
    return nest[idx1].v2 == nest[idx2].v2;
  }
};

struct NestedUVPCmp {
  UVertexPair* nest;
  bool use_v1;
  NestedUVPCmp(UVertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 < nest[idx2].v1;
    return nest[idx1].v2 < nest[idx2].v2;
  }
};
struct NestedUVPEq {
  UVertexPair* nest;
  bool use_v1;
  NestedUVPEq(UVertexPair* _nest, bool _use_v1) : nest(_nest), use_v1(_use_v1) {}
  inline bool operator() (uintE idx1, uintE idx2) {
    if (use_v1) return nest[idx1].v1 == nest[idx2].v1;
    return nest[idx1].v2 == nest[idx2].v2;
  }
};

struct nonZeroF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != 0);}};
struct nonZeroPairF{inline bool operator() (pair<uintE,uintE> &a) {return (a.second != 0);}};
struct greaterOneF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) > 1);}};
template<class T> struct cmpF{inline bool operator() (T a, T b) {return a < b;}};

struct nonEmptyUVPF{inline bool operator() (UVertexPair &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX);}};
struct nonEmptyUVPHF{inline bool operator() (UVertexPairHash &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX);}};
struct nonEmptyUWF{inline bool operator() (UWedge &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX || a.u != UINT_E_MAX);}};

struct nonMaxTupleF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != UINT_E_MAX || get<0>(a) != UINT_E_MAX);}};

struct uintELt {inline bool operator () (uintE a, uintE b) {return a < b;};};
struct uintEEq {inline bool operator() (uintE a, uintE b) {return a == b;};};
struct uintETupleLt {inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<0>(a) < get<0>(b);} };
struct uintETupleEq {inline bool operator() (tuple<uintE,uintE> a,tuple<uintE,uintE> b) {return get<0>(a) == get<0>(b); }};
struct uintETupleAdd {
  inline tuple<uintE,uintE> operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) const {
    return make_tuple(get<0>(a), get<1>(a) + get<1>(b));
  };
};
struct uintETupleLtBoth {
  inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {
    if (get<0>(a) == get<0>(b)) return get<1>(a) < get<1>(b);
      return get<0>(a) < get<0>(b);
  }
};

template<class T> struct refl{inline T operator() (T obj) {return obj;}};
template<class T> struct reflCount{inline long operator() (T obj) {return 1;}};
struct UVertexPairV2{inline uintE operator() (UVertexPair obj) {return obj.v2;}};
struct choose2{inline uintE operator() (uintE obj) {return obj*(obj-1)/2;}};
struct uintETupleGet0{inline uintE operator() (tuple<uintE,uintE> obj) {return get<0>(obj);}};
struct uintECount{inline long operator() (tuple<uintE,uintE> obj) {return get<1>(obj);}};

template <class E>
struct writeAddArr {
  E* arr;
  writeAddArr(E* _arr) : arr(_arr) {}
  inline void operator() (long idx, E num) {
    writeAdd(&arr[idx], num);
  }
};

template <class E>
struct writeAddSet {
  sparseAdditiveSet<E> set;
  writeAddSet(sparseAdditiveSet<E> _set) : set(_set) {}
  inline void operator() (long idx, E num) {
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

template<class vertex>
bipartiteGraph<vertex> KONECTToBp(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    char* S2 = newA(char,S.n);
    //ignore starting lines with '#' and find where to start in file 
    long k=0;
    while(1) {
      if(S.A[k] == '%') {
	while(S.A[k++] != '\n') continue;
      }
      if(k >= S.n || S.A[k] != '%') break; 
    }
    
    parallel_for(long i=0;i<S.n-k;i++) S2[i] = S.A[k+i];
    S.del();

    long spaces = 0;
    bool prev_space = true;
    long t = 0;
    while(1) {
      if (!isSpace(S2[t]) && prev_space) {
        prev_space = false;
        spaces++;
      }
      else if(isspace(S2[t])) {
        prev_space = true;
      }
      if(S2[t] == '\n') break;
      t++;
    }
  
    words W = stringToWords(S2, S.n-k);
    long m = W.m / spaces;
    using T = pair<uintE, uintE>;
    T *edges = newA(T,m);

    parallel_for(long i=0; i < m; i++) {
	    edges[i] = make_pair(atol(W.Strings[spaces*i]), atol(W.Strings[spaces*i + 1]));
    }

    // Remove duplicates
    quickSort(edges, m, uintETupleLtBoth());
    T lastRead = make_pair(UINT_E_MAX, UINT_E_MAX);
    long offset = 0;
    for (long i=0; i < m; ++i) {
      if (!(edges[i].first == lastRead.first && edges[i].second == lastRead.second)) {
        edges[offset] = edges[i];
        offset++;
        lastRead = edges[i];
      }
    }
    m = offset;

    long maxV = 0, maxU = 0;
    for (long i=0; i < m; i++) {
      maxV = max<intT>(maxV, edges[i].first);
      maxU = max<intT>(maxU, edges[i].second);
    }
    maxV++; maxU++;
    long nv = maxV;
    long nu = maxU;

    uintE* degV = newA(uintE, nv);
    uintE* degU = newA(uintE, nu);
    parallel_for(long i=0; i < nv; ++i) {degV[i] = 0;}
    parallel_for(long i=0; i < nu; ++i) {degU[i] = 0;}

    uintE* idxV = newA(uintE, nv);
    uintE* idxU = newA(uintE, nu);
    parallel_for(long i=0; i < nv; ++i) {idxV[i] = 0;}
    parallel_for(long i=0; i < nu; ++i) {idxU[i] = 0;}

    for (long i=0; i < m; i++) {
      degV[edges[i].first]++;
      degU[edges[i].second]++;
    }

    vertex* v = newA(vertex,nv+1);
    vertex* u = newA(vertex,nu+1);
  
  // Add all neighbors to v
  parallel_for(long i=0;i<nv;++i) {
    uintE* neighbors_v = newA(uintE,degV[i]);
    v[i] = vertex(neighbors_v, degV[i]);
  }
  // Add all neighbors to u
  parallel_for(long i=0;i<nu;++i) {
    uintE* neighbors_u = newA(uintE,degU[i]);
    u[i] = vertex(neighbors_u, degU[i]);
  }

  for(long i=0;i<m;++i) {
    uintE v_idx = edges[i].first;
    uintE u_idx = edges[i].second;
    (v[v_idx].neighbors)[idxV[v_idx]] = u_idx;
    idxV[v_idx]++;
    (u[u_idx].neighbors)[idxU[u_idx]] = v_idx;
    idxU[u_idx]++;
  }

  free(edges);
  free(degV);
  free(degU);
  free(idxV);
  free(idxU);
    
    Uncompressed_Membipartitegraph<vertex>* mem = new Uncompressed_Membipartitegraph<vertex>(v,u,nv,nu);
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

template <class T>
tuple<long, T*> intersect_seq(T* a, T* b, size_t num_a, size_t num_b) {
  if (num_b < num_a) return intersect_seq(b, a, num_b, num_a);
  //sampleSort(a,num_a,cmpF<T>());
  uintE* ret = newA(uintE,num_a);
  long idx = 0;

  auto a_map = make_in_imap<uintT>(num_a, [&] (size_t i) { return a[i]; });
  auto lt = [] (const T& l, const T& r) { return l < r; };

  for(size_t i=0;i<num_b;++i) {
    size_t find_idx = pbbs::linear_search(a_map, b[i], lt);
    if(a[find_idx] == b[i]) ret[idx++] = b[i];
  }
  return make_tuple(idx, ret);
}

template <class T>
tuple<long, T*> intersect(T* a, T* b, size_t num_a, size_t num_b) {
  if (num_a < 1000 && num_b < 1000) return intersect_seq(a,b,num_a,num_b);
  if (num_b < num_a) return intersect(b, a, num_b, num_a);

  uintE* a_cpy = newA(uintE, num_a);
  parallel_for(size_t i=0; i < num_a; ++i) {a_cpy[i] = a[i];};

  sampleSort(a_cpy,num_a,cmpF<T>());
  uintE* ret = newA(uintE,num_a);
  long idx = 0;

  auto a_map = make_in_imap<uintT>(num_a, [&] (size_t i) { return a_cpy[i]; });
  auto lt = [] (const T& l, const T& r) { return l < r; };

  for(size_t i=0;i<num_b;++i) {
    size_t find_idx = pbbs::binary_search(a_map, b[i], lt);
    if(a[find_idx] == b[i]) ret[idx++] = b[i];
  }
  return make_tuple(idx, ret);
}

tuple<long, uintE*> intersect_hist(uintE* a, uintE* b, size_t num_a, size_t num_b, uintE max) {
  if (num_a < 1000 && num_b < 1000) return intersect_seq(a,b,num_a,num_b);
  using X = tuple<uintE,uintE>;
  size_t num_total = num_a + num_b;
  uintE* total = newA(uintE, num_total);

  parallel_for(long i=0; i < num_a; ++i) {total[i] = a[i]; }
  parallel_for(long i=0; i < num_b; ++i) {total[i+num_a] = b[i]; }

  pbbsa::sequence<uintE> seq_total = pbbsa::sequence<uintE>(total, num_total);
  tuple<size_t,X*> hist_total = pbbsa::sparse_histogram<uintE,uintE>(seq_total, max);

  X* arr = newA(X, get<0>(hist_total));
  long len = sequence::filter(get<1>(hist_total), arr, get<0>(hist_total), greaterOneF());

  uintE* ret = newA(uintE,len);
  parallel_for(long i=0; i < len; ++i) { ret[i] = get<0>(arr[i]); }

  free(total);
  free(get<1>(hist_total));
  free(arr);

  return make_tuple(len, ret);
}

tuple<long, uintE*> intersect_hash_seq(uintE* a, uintE* b, size_t num_a, size_t num_b) {
  if (num_b < num_a) return intersect_hash_seq(b, a, num_b, num_a);

  sparseAdditiveSet<uintE> set = sparseAdditiveSet<uintE>(num_a, 1, UINT_E_MAX);
  for(size_t i=0; i<num_a; ++i) { set.insert(make_pair(a[i],0)); }

  uintE* ret = newA(uintE,num_a);
  long idx = 0;
  for(size_t i=0; i<num_b; ++i) {
    if(set.find(b[i]).first != UINT_E_MAX) ret[idx++] = b[i];
  }

  set.del();
  return make_tuple(idx, ret);
}

tuple<long, uintE*> intersect_hash(uintE* a, uintE* b, size_t num_a, size_t num_b) {
  //if (num_b < num_a) return intersect_hash(b, a, num_b, num_a);
  if (num_a < 1000 && num_b < 1000) return intersect_hash_seq(a,b,num_a,num_b);
  using X=pair<uintE,uintE>;

  sparseAdditiveSet<uintE> set = sparseAdditiveSet<uintE>(num_a, 1, UINT_E_MAX);
  parallel_for(size_t i=0; i<num_a; ++i) { set.insert(make_pair(a[i],0)); }
  parallel_for(size_t i=0; i<num_b; ++i) { set.insert(make_pair(b[i],1)); }

  _seq<pair<uintE,uintE>> seq = set.entries();
  X* seq_f = newA(X, seq.n);
  long len = sequence::filter(seq.A, seq_f, seq.n, nonZeroPairF());

  uintE* ret = newA(uintE,len);
  parallel_for(long i=0; i < len; ++i) { ret[i] = seq_f[i].first; }

  free(seq_f);
  seq.del();
  set.del();
  return make_tuple(len, ret);
}

//********************************************************************************************
//********************************************************************************************

//TODO we don't use nbhrs or save nbhrs -- delete from everything
template<class vertex>
pair<uintE*, uintE**> countWedgesScan(long nu, vertex* V, vertex* U, bool half=false, bool save_nbhrs=false) {
  uintE* idxs = newA(uintE, nu + 1);
  idxs[nu] = 0;

  using T = uintE*;
  T* nbhd_idxs = newA(T, nu);

  parallel_for(long i=0; i < nu; ++i) {
    idxs[i] = 0;
    uintE u_deg = U[i].getOutDegree();

    nbhd_idxs[i] = newA(uintE, u_deg + 1);
    parallel_for(long j=0; j < u_deg+1; ++j) {(nbhd_idxs[i])[j] = 0;}

    parallel_for(long j=0; j < u_deg; ++j) { //TODO can parallelize this too technically
      if (!half) (nbhd_idxs[i])[j] = V[U[i].getOutNeighbor(j)].getOutDegree() - 1; //idxs[i] +=
      else {
        (nbhd_idxs[i])[j] = 0;
        vertex v = V[U[i].getOutNeighbor(j)];
        for (long k = 0; k < v.getOutDegree(); ++k) { //TODO can parallelize this too technically
          if (v.getOutNeighbor(k) < i) nbhd_idxs[i][j] ++;
          else break;
        }
      }
    }

    sequence::plusScan(nbhd_idxs[i], nbhd_idxs[i], u_deg + 1);
    // Set up indices associated with u
    idxs[i] = (nbhd_idxs[i])[u_deg];
    //if (!save_nbhrs)
    free(nbhd_idxs[i]);
  }
  //if (!save_nbhrs)
  free(nbhd_idxs);
  sequence::plusScan(idxs, idxs, nu + 1);

  //return save_nbhrs ? make_pair(idxs, nbhd_idxs) : make_pair(idxs, nullptr);
  return make_pair(idxs, nullptr);
}

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
  return sequence::reduce<long>((long) 0, nv, addF<long>(), chooseV<vertex, long>(V)); //TODO CHECK make sure 2* not needed
}

template<class vertex>
long countWedges_seq(long nv, vertex* V) {
  long num_wedges = 0;
  for(long i=0; i < nv; ++i) {
    uintE deg_v = V[i].getOutDegree();
    num_wedges += deg_v * (deg_v - 1) / 2;
  }
  return num_wedges; //TODO CHECK
}

template <class vertex>
pair<bool,long> cmpWedgeCounts_seq(bipartiteGraph<vertex> GA) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = countWedges_seq<vertex>(nv,GA.V);
  long num_wedges_u = countWedges_seq<vertex>(nu, GA.U);
  return make_pair((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
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

//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class vertex, class wedgeCons, class Sequence>
wedge* _getWedges_seq2(Sequence I, const long nu, const vertex* V, const vertex* U, wedgeCons cons, long num_wedges,
  long curr_idx, long next_idx) {
  wedge* seagulls = newA(wedge, num_wedges);
  long idx = 0;
  for(long i=curr_idx; i < next_idx; ++i) {
    uintE u_idx = I[i];
    vertex u = U[u_idx];
    for(long j=0; j < u.getOutDegree(); ++j) {
      vertex v = V[u.getOutNeighbor(j)];
      // Find neighbors (not equal to u) of v
      for (long k = 0; k < v.getOutDegree(); ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx != u_idx) {
          seagulls[idx] = cons(u_idx, u2_idx, u.getOutNeighbor(j));
          ++idx;
        }
      }
    }
  }
  return seagulls;
}

template<class wedge, class vertex, class wedgeCons, class Sequence, class F>
wedge* _getWedges2(Sequence I, const long nu, const vertex* V, const vertex* U, wedgeCons cons,
  long num_wedges, long curr_idx=0, long next_idx=-1) {
  if (next_idx == -1) next_idx = nu;
  if (num_wedges < 10000) return _getWedges_seq2<wedge>(I, nu, V, U, cons, num_wedges, curr_idx, next_idx);
  // First, we must retrive the indices at which to store each seagull (so that storage can be parallelized)
  // Array of indices associated with seagulls for each active vertex
  long* sg_idxs = newA(long, next_idx - curr_idx + 1);
  sg_idxs[next_idx-curr_idx] = 0;

  // 2D array of indices associated with seagulls for each neighbor of each active vertex
  using T = long*;
  T* nbhd_idxs = newA(T, next_idx-curr_idx); 

  parallel_for(long i=curr_idx; i < next_idx; ++i) {
    // Set up for each active vertex
    const uintE u_idx = I[i];
    const uintE u_deg = U[u_idx].getOutDegree();
    // Allocate space for indices associated with each neighbor of u
    nbhd_idxs[i-curr_idx] = newA(long, u_deg + 1);
    (nbhd_idxs[i-curr_idx])[u_deg] = 0;
    // Assign indices associated with each neighbor of u
    parallel_for(long j=0; j < u_deg; ++j) {
      (nbhd_idxs[i-curr_idx])[j] = V[U[u_idx].getOutNeighbor(j)].getOutDegree()-1;
    }
    sequence::plusScan(nbhd_idxs[i-curr_idx], nbhd_idxs[i-curr_idx], u_deg+1);
    // Set up indices associated with u
    sg_idxs[i-curr_idx] = (nbhd_idxs[i-curr_idx])[u_deg];
  }
  // Assign indices associated with each active vertex
  sequence::plusScan(sg_idxs, sg_idxs, next_idx-curr_idx + 1);

  // Allocate space for seagull storage
  long num_sg = sg_idxs[next_idx-curr_idx];
  wedge* seagulls = newA(wedge, num_sg);

  // Store seagulls in parallel
  parallel_for(long i=curr_idx; i < next_idx; ++i) {
    uintE u_idx = I[i];
    long sg_idx = sg_idxs[i-curr_idx];
    // Consider each neighbor v of active vertex u
    parallel_for(long j=0; j < U[u_idx].getOutDegree(); ++j) {
      const vertex v = V[U[u_idx].getOutNeighbor(j)];
      const uintE v_deg = v.getOutDegree();
      long nbhd_idx = (nbhd_idxs[i-curr_idx])[j];
      long idx = 0;
      // Find neighbors (not equal to u) of v
      for (long k = 0; k < v_deg; ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx != u_idx) {
          seagulls[sg_idx+nbhd_idx+idx] = cons(u_idx, u2_idx, U[u_idx].getOutNeighbor(j));
          ++idx;
        }
      }
    }
  }

  // Cleanup
  parallel_for(long i=curr_idx; i<next_idx;++i) { 
    free(nbhd_idxs[i]); 
  }
  free(nbhd_idxs);
  free(sg_idxs);

  return seagulls;
}

template<class vertex, class wedgeCons, class T, class Sequence>
void _getWedgesHash2(T& wedges, Sequence I, long nu, vertex* V, vertex* U, wedgeCons cons, long num_wedges, long curr_idx=0, long next_idx=-1) {
  if (next_idx == -1) next_idx = nu;
  wedges.resize(num_wedges);
  hashInsertTimer.start();
  parallel_for(long i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    const uintE u_idx = I[i];
    const uintE u_deg = U[u_idx].getOutDegree();

    parallel_for (long j=0; j < u_deg; ++j ) {
        const vertex v = V[U[u_idx].getOutNeighbor(j)];
        const uintE v_deg = v.getOutDegree();
        // Find all seagulls with center v and endpoint u
        for (long k=0; k < v_deg; ++k) { 
          const uintE u2_idx = v.getOutNeighbor(k);

          if (u2_idx != u_idx) wedges.insert(make_pair(cons(u_idx,u2_idx,U[u_idx].getOutNeighbor(j)),1)); //(u_idx * nu + u2_idx, 1)
	}
    }
  }
  hashInsertTimer.stop();
}

template<class vertex, class Sequence>
pair<long,long> getNextWedgeIdx_seq2(Sequence I, long nu, vertex* V, vertex* U, long max_wedges, long curr_idx) {
  long orig = max_wedges;
  for(long i=curr_idx; i < nu; ++i) {
    uintE u_idx = I[i];
    for(long j=0; j < U[u_idx].getOutDegree(); ++j) {
      uintE num = V[U[u_idx].getOutNeighbor(j)].getOutDegree() - 1;
      if (num > max_wedges) {
        if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
        return make_pair(i, orig - max_wedges);
      }
      else { max_wedges -= num; }
    }
  }
  return make_pair(nu, orig - max_wedges);
}

// TODO based on half can also optimize this stuff? but will be hard to parallelize later so maybe not
//JS: there looks to be some inefficiency here. we should have gotten the prefix sum of wedges for all vertices at the beginning, so we don't need to recompute it here. binary search can use a doubling search instead, starting at the last index until getting the right range, then binary search inside.
template<class vertex, class Sequence>
pair<long,long> getNextWedgeIdx2(Sequence I, long nu, vertex* V, vertex* U, long max_wedges, long curr_idx) {
  nextWedgeTimer.start();
  if (nu - curr_idx < 10000) return getNextWedgeIdx_seq2(I, nu, V, U, max_wedges, curr_idx);
  uintE* idxs = newA(uintE, nu - curr_idx + 1);
  idxs[nu-curr_idx] = 0;
  parallel_for(long i=curr_idx; i < nu; ++i) {
    idxs[i-curr_idx] = 0;
    uintE u_idx = I[i];
    for(long j=0; j < U[u_idx].getOutDegree(); ++j) {
      idxs[i-curr_idx] += V[U[u_idx].getOutNeighbor(j)].getOutDegree() - 1;
    }
  }
  sequence::plusScan(idxs, idxs, nu - curr_idx + 1);

  auto idx_map = make_in_imap<uintT>(nu - curr_idx, [&] (size_t i) { return idxs[i+1]; });
  auto lte = [] (const uintE& l, const uintE& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  long num = idxs[find_idx - curr_idx];
  free(idxs);
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
  nextWedgeTimer.stop();
  return make_pair(find_idx, num); //TODO make sure right
}

template<class wedge, class vertex, class wedgeCons, class Sequence, class F>
pair<pair<wedge*,long>, long> getWedges2(Sequence I, const long nu, const vertex* V, const vertex* U, wedgeCons cons,
  long max_wedges, long curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) return make_pair(make_pair(_getWedges2<wedge>(I, nu, V, U, cons, num_wedges), num_wedges), nu);

  pair<long, long> p = getNextWedgeIdx2(I, nu, V, U, max_wedges, curr_idx);
  long next_idx = p.first;
  num_wedges = p.second;
  return make_pair(make_pair(_getWedges2<wedge>(I, nu, V, U, cons, num_wedges, curr_idx, next_idx), num_wedges), next_idx);
}

template<class vertex, class wedgeCons, class T, class Sequence>
long getWedgesHash2(T& wedges, Sequence I, long nu, vertex* V, vertex* U, wedgeCons cons, long max_wedges, long curr_idx, long num_wedges) {
  if (max_wedges >= num_wedges) {
    _getWedgesHash2(wedges, I, nu, V, U, cons, num_wedges);
    return nu;
  }
  pair<long, long> p = getNextWedgeIdx2(I, nu, V, U, max_wedges, curr_idx);
  long next_idx = p.first;
  num_wedges = p.second;
  _getWedgesHash2(wedges, I, nu, V, U, cons, num_wedges, curr_idx, next_idx);
  return next_idx;
}

//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class vertex, class wedgeCons>
wedge* _getWedges_seq(const long nu, const vertex* V, const vertex* U, wedgeCons cons, long num_wedges,
  long curr_idx, long next_idx) {
  wedge* seagulls = newA(wedge, num_wedges);
  long idx = 0;
  for(long i=curr_idx; i < next_idx; ++i) {
    for(long j=0; j < U[i].getOutDegree(); ++j) {
      vertex v = V[U[i].getOutNeighbor(j)];
      // Find neighbors (not equal to u) of v
      for (long k = 0; k < v.getOutDegree(); ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx < i) {
          seagulls[idx] = cons(i, u2_idx, U[i].getOutNeighbor(j));
          ++idx;
        }
        else break; 
      }
    }
  }
  return seagulls;
}

template<class wedge, class vertex, class wedgeCons>
wedge* _getWedges(const long nu, const vertex* V, const vertex* U, wedgeCons cons,
  long num_wedges, uintE* wedge_idxs, uintE** nbhd_idxs, long curr_idx=0, long next_idx=-1) {
  if (next_idx == -1) next_idx = nu;
  if (num_wedges < 10000) return _getWedges_seq<wedge>(nu, V, U, cons, num_wedges, curr_idx, next_idx);
 
  // Allocate space for seagull storage
  long num_sg = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  wedge* seagulls = newA(wedge, num_sg);
  // Store seagulls in parallel
  parallel_for(long i=curr_idx; i < next_idx; ++i) {
    long sg_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    // Consider each neighbor v of active vertex u
    //parallel_for(long j=0; j < U[i].getOutDegree(); ++j) {
    long idx = 0;
    for(long j=0; j < U[i].getOutDegree(); ++j) {
      const vertex v = V[U[i].getOutNeighbor(j)];
      //long nbhd_idx = (nbhd_idxs[i])[j];
      //long idx = 0;
      // Find neighbors (not equal to u) of v
      for (long k = 0; k < v.getOutDegree(); ++k) {
        uintE u2_idx = v.getOutNeighbor(k);
        if (u2_idx < i) {
          seagulls[sg_idx+idx] = cons(i, u2_idx, U[i].getOutNeighbor(j));
          ++idx;
        }
        else break;
      }
    }
  }
  /*parallel_for(long i=curr_idx; i<next_idx;++i) { 
    free(nbhd_idxs[i]); 
  }*/

  return seagulls;
}

template<class vertex, class wedgeCons, class T>
void _getWedgesHash(T& wedges,long nu, vertex* V, vertex* U, wedgeCons cons, long num_wedges, long curr_idx=0, long next_idx=-1) {
  if (next_idx == -1) next_idx = nu;
  wedges.resize(num_wedges);
  hashInsertTimer.start();
  parallel_for(long i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    parallel_for (long j=0; j < U[i].getOutDegree(); ++j ) {
        const vertex v = V[U[i].getOutNeighbor(j)];
        // Find all seagulls with center v and endpoint u
        for (long k=0; k < v.getOutDegree(); ++k) { 
          const uintE u2_idx = v.getOutNeighbor(k);
          if (u2_idx < i) wedges.insert(make_pair(cons(i,u2_idx,U[i].getOutNeighbor(j)),1)); //(u_idx * nu + u2_idx, 1)
          else break;
        }
    }
  }
  hashInsertTimer.stop();
}

template<class vertex>
long getNextWedgeIdx_seq(long nu, vertex* V, vertex* U, long max_wedges, long curr_idx) {
  long orig = max_wedges;
  for(long i=curr_idx; i < nu; ++i) {
    for(long j=0; j < U[i].getOutDegree(); ++j) {
      vertex v = V[U[i].getOutNeighbor(j)];
      uintE num = 0;
      for (long k=0; k < v.getOutDegree(); ++k) {
        if (v.getOutNeighbor(k) < i) num ++;
        else break;
      }
      if (num > max_wedges) {
        if (i == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
        return i;
      }
      else { max_wedges -= num; }
    }
  }
  return nu;
}

// TODO doubling search
template<class vertex>
long getNextWedgeIdx(long nu, vertex* V, vertex* U, long max_wedges, long curr_idx, uintE* wedge_idxs) {
  nextWedgeTimer.start();
  if (nu - curr_idx < 2000) return getNextWedgeIdx_seq(nu, V, U, max_wedges, curr_idx);

  auto idx_map = make_in_imap<uintT>(nu - curr_idx, [&] (size_t i) { return wedge_idxs[curr_idx+i+1] - wedge_idxs[curr_idx]; });
  auto lte = [] (const uintE& l, const uintE& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
  nextWedgeTimer.stop();
  return find_idx; //TODO make sure right
}

template<class wedge, class vertex, class wedgeCons>
pair<pair<wedge*,long>, long> getWedges(const long nu, const vertex* V, const vertex* U, wedgeCons cons, long max_wedges, long curr_idx, long num_wedges, uintE* wedge_idxs, uintE** nbhd_idxs) {
  if (max_wedges >= num_wedges) return make_pair(make_pair(_getWedges<wedge>(nu, V, U, cons, num_wedges, wedge_idxs, nbhd_idxs), num_wedges), nu);
  long next_idx = getNextWedgeIdx(nu, V, U, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  wedge* wedges = _getWedges<wedge>(nu, V, U, cons, num_wedges, wedge_idxs, nbhd_idxs, curr_idx, next_idx);
  return make_pair(make_pair(wedges, num_wedges), next_idx);
}

template<class vertex, class wedgeCons, class T>
long getWedgesHash(T& wedges, long nu, vertex* V, vertex* U, wedgeCons cons, long max_wedges, long curr_idx, long num_wedges, uintE* wedge_idxs) {
if (max_wedges >= num_wedges) {
    _getWedgesHash(wedges, nu, V, U, cons, num_wedges);
    return nu;
  }
  long next_idx = getNextWedgeIdx(nu, V, U, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedgesHash(wedges, nu, V, U, cons, num_wedges, curr_idx, next_idx);
  return next_idx;
}

//***************************************************************************************************
//***************************************************************************************************

pair<tuple<uintE,uintE>*,long> updateBuckets(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE,uintE>;
  X* update = newA(X,num_updates);

    // Filter for bucket updates
  parallel_for(long i=0; i < num_updates; ++i) {
    const uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if (old_b > k) {
        uintE new_b = max(butterflies[u_idx],k);
        D.s[u_idx] = new_b;
        uintE new_bkt = b.get_bucket(old_b, new_b);
        update[i] = make_tuple(u_idx, new_bkt);
    }
    else {update[i] = make_tuple(UINT_E_MAX,UINT_E_MAX);}
  }

  X* update_filter = newA(X, num_updates);
  long num_updates_filter = sequence::filter(update,update_filter,num_updates, nonMaxTupleF());
  free(update);
  return make_pair(update_filter,  num_updates_filter);
}

pair<tuple<uintE,uintE>*,long> updateBuckets_seq(uintE* update_idxs, long num_updates, uintE* butterflies, 
  array_imap<uintE> D, buckets<array_imap<uintE>> b, uintE k) {
  using X = tuple<uintE, uintE>;
  X* update = newA(X, num_updates);
  long idx = 0;
  for(long i=0; i < num_updates; ++i) {
    uintE u_idx = update_idxs[i];
    uintE old_b = D.s[u_idx];

    if(old_b > k) {
      uintE new_b = max(butterflies[u_idx], k);
      D.s[u_idx] = new_b;
      uintE new_bkt = b.get_bucket(old_b, new_b);
      update[idx] = make_tuple(u_idx, new_bkt);
      ++idx;
    }
  }
  return make_pair(update, idx);
}

template <class vertex>
struct edgeToIdx { 
  long nv;
  long nu;
  vertex* V;
  vertex* U;
  bool overflow;
  long num_edges;
  long max_wedges;
  sparseAdditiveSet<uintE> edges;
  sparsePointerAdditiveSet<UVertexPairHash, uintE, UVertexPairHashEq> edges_overflow;

  edgeToIdx(bipartiteGraph<vertex> GA, bool use_v, long _max_wedges) : max_wedges(_max_wedges) {
    nv = use_v ? GA.nv : GA.nu;
    nu = use_v ? GA.nu : GA.nv;
    V = use_v ? GA.V : GA.U;
    U = use_v ? GA.U : GA.V;

    overflow = (nu > UINT_E_MAX / nu);
  
    if (nv < nu) num_edges = sequence::reduce<long>((long) 0, (long) nv, addF<long>(), getV<vertex, long>(V));
    else num_edges = sequence::reduce<long>((long) 0, (long) nu, addF<long>(), getV<vertex, long>(U));
  
    if (overflow)
      edges_overflow = sparsePointerAdditiveSet<UVertexPairHash, uintE, UVertexPairHashEq>(num_edges, 1, UINT_E_MAX, UVertexPairHashEq());
    else
      edges = sparseAdditiveSet<uintE>(num_edges, 1, UINT_E_MAX);
    
    if (!overflow && nu*nv < max_wedges) {
// go through edges in an array
      bool* edges_bool = newA(bool, nu*nv);
      parallel_for(long i=0; i < nu*nv; ++i) {edges_bool[i] = false;}
      parallel_for(long i=0; i < nv; ++i) {
        for (long j = 0; j < V[i].getOutDegree(); ++j) {
          edges_bool[nu*i + V[i].getOutNeighbor(j)] = true;
        }
      }
      auto f = [&] (size_t i) { return edges_bool[i]; };
      auto f_in = make_in_imap<bool>(nu*nv, f);
      auto out = pbbs::pack_index<uintE>(f_in);
      out.alloc = false;
      free(edges_bool);

      parallel_for(long i=0; i < out.size(); ++i) {
        edges.insert(make_pair(out.s[i], i));
      }

      free(out.s);
    }
    else {
// hash edges sequentially
// TODO maybe it's faster to hash all w/val 1 or something, then retrieve entries, then rehash w/entry idx??
      uintE idx = 0;
      UVertexPairHashCons overflow_cons = UVertexPairHashCons(nu);
      for (long i=0; i < nv; ++i) {
        for (long j=0; j < V[i].getOutDegree(); ++j) {
          if (overflow) edges_overflow.insert(make_pair(overflow_cons(i, V[i].getOutNeighbor(j), 0), idx));
          else edges.insert(make_pair(nu*i + V[i].getOutNeighbor(j), idx));
          idx++;
        }
      }
    }
  }

  void del() {
    if (overflow) edges_overflow.del();
    else edges.del();
  }

  // Note: always in the format (V, U)
  inline uintE operator() (const uintE& i, const uintE& j) {
    if (overflow) {
      UVertexPairHashConsN overflow_cons = UVertexPairHashConsN(nu);
      return (edges_overflow.find(overflow_cons(i, j, 0))).second;
    }
    return (edges.find((uintE) (nu*i + j))).second;
  }
};


#endif
