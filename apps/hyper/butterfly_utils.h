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

#define MAX_STEP_SIZE 1000

using namespace std;

struct UWedge {
  uintE v1;
  uintE v2;
  uintE u;
  intT j; intT k;
UWedge(uintE _v1, uintE _v2, uintE _u, intT _j, intT _k) : v1(_v1), v2(_v2), u(_u), j(_j), k(_k) {}
};

struct UWedgeCons { inline UWedge operator() (uintE v1, uintE v2, uintE c, intT j, intT k) { return UWedge(v1, v2, c, j, k); }};

struct UWedgeIntRankCons {
  long nu;
UWedgeIntRankCons(long _nu) : nu(_nu) {}
  inline long operator() (uintE v1, uintE v2, uintE c, intT j, intT k) {
    return ((((long) v1) * nu + ((long) v2 >> 1)) << 1) + (v2 & 0b1);
  }
};

struct UWedgeCmp {
  inline bool operator() (UWedge vs1, UWedge vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

struct UWedgeEq { inline bool operator() (UWedge vs1, UWedge vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Represents a pair of vertices on one side of a bipartite graph (ordered, with least vertex first)

// Represents a pair of vertices on one side of a bipartite graph (unordered, stored based on constructor order)
struct UVertexPair {
  uintE v1;
  uintE v2;
UVertexPair(uintE _v1, uintE _v2) : v1(_v1), v2(_v2) {}
};

//TODO get rid of legacy _nv
// Comparer for VertexPair based on least vertex in pair and then greatest vertex in pair

// Comparer for VertexPair based on greatest vertex in pair and then least vertex in pair

// Comparer for UVertexPair
struct UVertexPairCmp2{
  UVertexPairCmp2() {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v2 == vs2.v2) return vs1.v1 < vs2.v1;
    return vs1.v2 < vs2.v2;
  }
};

struct UVertexPairCmp {
  UVertexPairCmp() {}
  inline bool operator() (UVertexPair vs1, UVertexPair vs2) {
    if (vs1.v1 == vs2.v1) return vs1.v2 < vs2.v2;
    return vs1.v1 < vs2.v1;
  }
};

struct UVPFirst { uintE operator() (const UVertexPair& x) {return x.v1;}};
struct UVPSecond { uintE operator() (const UVertexPair& x) {return x.v2;}};
struct UWFirst { uintE operator() (const UWedge& x) {return x.v1;}};
struct UWSecond { uintE operator() (const UWedge& x) {return x.v2;}};
template<class T, class X>
  struct tupleFirst {T operator() (tuple<T,X> a) {return get<0>(a);} };

// Equality for VertexPair and UVertexPair
struct UVertexPairEq { inline bool operator() (UVertexPair vs1, UVertexPair vs2) { return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);} };

// Constructs a VertexPair and UVertexPair
struct UVertexPairCons { inline UVertexPair operator() (uintE v1, uintE v2, uintE c, intT j, intT k) { return UVertexPair(v1, v2); }};

struct UVertexPairRankCons {
  inline UVertexPair operator() (uintE v1, uintE v2, bool b) { return UVertexPair(v1, (v2 << 1) + b); }
};

// Constructs a uintE form of a VertexPair and UVertexPair
struct UVertexPairIntCons {
  long nu;
UVertexPairIntCons(long _nu) : nu(_nu) {}
  inline long operator() (uintE v1, uintE v2, uintE c, intT j, intT k) {
    return (long) v1 * nu + (long) v2;
  }
};
struct UVertexPairIntRankCons {
  long nu;
UVertexPairIntRankCons(long _nu) : nu(_nu) {}
  inline long operator() (uintE v1, uintE v2, bool b) {
    return ((((long) v1) * nu + (long) v2) << 1) + (b ? 0b1 : 0);
  }
};

// Comparer for indices for VertexPairs in nest, by v1 or v2 (based on use_v1)
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

struct NestedUIECmp{
  uintE* nest;
NestedUIECmp(uintE* _nest) : nest(_nest) {}
  inline bool operator() (uintE idx1, uintE idx2) { //
    return nest[idx1] < nest[idx2];
  }
};

struct nonZeroF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != 0);}};
struct nonZeroPairF{inline bool operator() (pair<uintE,uintE> &a) {return (a.second != 0);}};
struct greaterOneF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) > 1);}};
struct greaterOneLongF{inline bool operator() (tuple<long,uintE> &a) {return (get<1>(a) > 1);}};
template<class T> struct cmpF{inline bool operator() (T a, T b) {return a < b;}};

struct nonEmptyUVPF{inline bool operator() (UVertexPair &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX);}};
struct nonEmptyUWF{inline bool operator() (UWedge &a) {return (a.v1 != UINT_E_MAX || a.v2 != UINT_E_MAX || a.u != UINT_E_MAX);}};

struct nonMaxTupleF{inline bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != UINT_E_MAX || get<0>(a) != UINT_E_MAX);}};

template<class T, class X>
  struct tupleLt {inline bool operator() (tuple<T,X> a, tuple<T,X> b) {return get<0>(a) < get<0>(b);} };
template<class T, class X>
  struct tupleEq {inline bool operator() (tuple<T,X> a,tuple<T,X> b) {return get<0>(a) == get<0>(b); }};
template<class T, class X>
  struct tupleAdd {
    inline tuple<T,X> operator() (tuple<T,X> a, tuple<T,X> b) const {
      return make_tuple(get<0>(a), get<1>(a) + get<1>(b));
    };
  };

struct uintELt {inline bool operator () (uintE a, uintE b) {return a < b;};};
struct uintEEq {inline bool operator() (uintE a, uintE b) {return a == b;};};
struct uintETupleLt {inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<0>(a) < get<0>(b);} };
struct uintETupleGt {inline bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<1>(a) > get<1>(b);} };
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
struct choose2{inline long operator() (long obj) {return obj*(obj-1)/2;}};
struct uintETupleGet0{inline uintE operator() (tuple<uintE,long> obj) {return get<0>(obj);}};
struct uintECount{inline long operator() (tuple<uintE,long> obj) {return get<1>(obj);}};

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

template<class K, class E>
  E getAdd (E curr, tuple<K,E> v) {
  return curr + get<1>(v);
}

template<class K, class E>
  tuple<K,E> getAddReduce (tuple<K,E> curr, tuple<K,E> v) {
  return make_tuple(get<0>(curr),get<1>(curr) + get<1>(v));
}

template<class K, class E>
  pair<K,E> getAddReducePair (pair<K,E> curr, pair<K,E> v) {
  return make_pair(v.first, curr.second + v.second);
}

//***********************************************************************************************
//***********************************************************************************************


//symmetric compact bipartite
struct bipartiteCSR {
  uintT *offsetsV, *offsetsU;
  uintE *edgesV, *edgesU;
  long nv, nu, numEdges;

bipartiteCSR () {}
bipartiteCSR(uintT* _offsetsV, uintT* _offsetsU, uintE* _edgesV, uintE* _edgesU, long _nv, long _nu, long _ne) :
  offsetsV(_offsetsV), offsetsU(_offsetsU), edgesV(_edgesV), edgesU(_edgesU), nv(_nv), nu(_nu), numEdges(_ne)
  {}

  void del() {
    free(offsetsV); free(offsetsU); free(edgesV); free(edgesU);
  }
};

struct graphCSR {
  uintT *offsets;
  uintE *edges;
  long n, numEdges;

graphCSR(uintT* _offsets, uintE* _edges, long _n, long _ne) :
  offsets(_offsets), edges(_edges), n(_n), numEdges(_ne)
  {}

  void del() {
    free(offsets); free(edges);
  }
};

template <class E>
struct rankWedgeF { 
  uintT* offsets;
  uintE* edges;
  uintE* rank;
  uintE* orank;
rankWedgeF(uintT* _offsets, uintE* _edges, uintE* _rank, uintE* _orank) :
  offsets(_offsets), edges(_edges), rank(_rank), orank(_orank) {}
  inline E operator() (const uintT& i) const {
    intT v_offset = offsets[i];
    intT v_deg = offsets[i+1]-v_offset;
    intT in_deg = 0; intT out_deg = 0;
    for(intT j=0; j < v_deg; ++j) {
      uintE u = edges[v_offset+j];
      if (orank[u] < rank[i]) in_deg++;
      else out_deg++;
    }
    return (E) ( in_deg * out_deg + ((in_deg * (in_deg-1)) / 2)); 
  }
};

template<class F>
tuple<uintE*, uintE*, uintE*> getRanks(bipartiteCSR& G, F samplesort_f) {
  uintE* ranks = newA(uintE, G.nv + G.nu);
  uintE* rankV = newA(uintE, G.nv);
  uintE* rankU = newA(uintE, G.nu);
  
  parallel_for(long v=0; v < G.nv+G.nu; ++v) { ranks[v] = v; }

  sampleSort(ranks,G.nv+G.nu, samplesort_f); 

  parallel_for(long i=0; i < G.nv + G.nu; ++i) {
    if (ranks[i] >= G.nv) rankU[ranks[i] - G.nv] = i;
    else rankV[ranks[i]] = i;
  }

  return make_tuple(ranks, rankV, rankU);
}

tuple<uintE*, uintE*, uintE*> getApproxDegRanks(bipartiteCSR& G) {
  auto samplesort_f = [&] (const uintE a, const uintE b) -> const uintE {
    uintE deg_a = (a >= G.nv) ? G.offsetsU[a-G.nv+1]-G.offsetsU[a-G.nv] : G.offsetsV[a+1]-G.offsetsV[a];
    uintE deg_b = (b >= G.nv) ? G.offsetsU[b-G.nv+1]-G.offsetsU[b-G.nv] : G.offsetsV[b+1]-G.offsetsV[b];
    uintE log_deg_a = deg_a <= 0 ? 0 : (uintE) floor(log2(deg_a)) + 1;
    uintE log_deg_b = deg_b <= 0 ? 0 : (uintE) floor(log2(deg_b)) + 1;
    return log_deg_a > log_deg_b;
  };
  return getRanks(G, samplesort_f);
}

tuple<uintE*, uintE*, uintE*> getDegRanks(bipartiteCSR& G) {
  auto samplesort_f = [&] (const uintE a, const uintE b) -> const uintE {
    uintE deg_a = (a >= G.nv) ? G.offsetsU[a-G.nv+1]-G.offsetsU[a-G.nv] : G.offsetsV[a+1]-G.offsetsV[a];
    uintE deg_b = (b >= G.nv) ? G.offsetsU[b-G.nv+1]-G.offsetsU[b-G.nv] : G.offsetsV[b+1]-G.offsetsV[b];
    return deg_a > deg_b;
  };
  return getRanks(G, samplesort_f);
}

template <class E>
struct degF { 
  uintT* offsetsV;
  uintT* offsetsU;
  vertexSubset active;
  long nv;
degF(long _nv, uintT* _offsetsV, uintT* _offsetsU, vertexSubset& _active) : nv(_nv), offsetsV(_offsetsV), offsetsU(_offsetsU), active(_active) {}
  inline E operator() (const uintT& i) const {
    bool use_v = active.vtx(i) < nv;
    intT idx = use_v ? active.vtx(i) : active.vtx(i) - nv;
    intT deg = (use_v ? offsetsV[idx+1] - offsetsV[idx] : offsetsU[idx+1] - offsetsU[idx]);
    return (E) (deg); 
  }
};



////////////////////////////////////////////////////////////////////////

template <class E, class VS>
  vertexSubsetData<E> edgeMapInduced(bipartiteCSR& GA, VS& V) {
  long nTo = GA.nv + GA.nu;
  uintT m = V.size();
  V.toSparse();
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
      bool use_v = V.vtx(i) < GA.nv;
      intT idx = use_v ? V.vtx(i) : V.vtx(i) - GA.nv;
      intT deg = (use_v ? GA.offsetsV[idx+1] - GA.offsetsV[idx] : GA.offsetsU[idx+1] - GA.offsetsU[idx]);
      degrees[i] = deg;
    });
  long outEdgeCount = pbbs::scan_add(degrees, degrees);
  if (outEdgeCount == 0) {
    return vertexSubsetData<E>(nTo);
  }
  typedef tuple<uintE, E> VE;
  VE* outEdges = pbbs::new_array_no_init<VE>(outEdgeCount);

  parallel_for (size_t i = 0; i < m; i++) {
    uintT o = degrees[i];
    bool use_v = V.vtx(i) < GA.nv;
    intT idx = use_v ? V.vtx(i) : V.vtx(i) - GA.nv;
    intT offset  = use_v ? GA.offsetsV[idx] : GA.offsetsU[idx];
    intT deg = (use_v ? GA.offsetsV[idx+1] : GA.offsetsU[idx+1]) - offset;
    granular_for(j,0,deg,deg > 10000, {
	intT nbhr = use_v ? GA.edgesV[offset+j] + GA.nv : GA.edgesU[offset+j];
	// must decrement D.s[nbhr]
	outEdges[o+j] = make_tuple(nbhr, 1);
      });
  }
  auto vs = vertexSubsetData<E>(nTo, outEdgeCount, outEdges);
  return vs;
}

template <class Val>
struct BipartiteProp {
  using K = uintE; // keys are always uintE's (vertex-identifiers)
  using KV = tuple<K, Val>;
  bipartiteCSR GA;
  pbbs::hist_table<K, Val> ht;

BipartiteProp(bipartiteCSR& _GA, KV _empty, size_t ht_size=numeric_limits<size_t>::max()) : GA(_GA) {
  if (ht_size == numeric_limits<size_t>::max()) {
    ht_size = 1L << pbbs::log2_up(GA.numEdges/20);
  } else { ht_size = 1L << pbbs::log2_up(ht_size); }
  ht = pbbs::hist_table<K, Val>(_empty, ht_size);
}

  // map_f: (uintE v, uintE ngh) -> E
  // reduce_f: (E, tuple(uintE ngh, E ngh_val)) -> E
  // apply_f: (uintE ngh, E reduced_val) -> O
  template <class O, class M, class Reduce, class Apply, class VS>
    inline vertexSubsetData<O> edgeMapReduce(VS& vs, Reduce& reduce_f, Apply& apply_f) {
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(GA.nv+GA.nu);
    }

    auto oneHop = edgeMapInduced<M, VS>(GA, vs);
    oneHop.toSparse();

    auto get_elm = make_in_imap<tuple<K, M> >(oneHop.size(), [&] (size_t i) { return oneHop.vtxAndData(i); });
    auto get_key = make_in_imap<uintE>(oneHop.size(), [&] (size_t i) -> uintE { return oneHop.vtx(i); });

    auto q = [&] (sequentialHT<K, Val>& S, tuple<K, M> v) -> void { S.template insertF<M>(v, reduce_f); };
    auto res = pbbs::histogram_reduce<tuple<K, M>, tuple<K, O> >(get_elm, get_key, oneHop.size(), q, apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(GA.nv+GA.nu, res.first, res.second);
  }

  template <class O, class Apply, class VS>
    inline vertexSubsetData<O> bpedgePropCount(VS& vs, Apply& apply_f) {
    auto reduce_f = [&] (const uintE& cur, const tuple<uintE, uintE>& r) { return cur + 1; };
    return edgeMapReduce<O, uintE>(vs,reduce_f, apply_f);
  }

  ~BipartiteProp() {
    ht.del();
  }
};

// TODO serialize for small buckets, log buckets
tuple<uintE*,uintE*,uintE*> getApproxCoCoreRanks(bipartiteCSR& GA, size_t num_buckets=128) {
  using X = tuple<bool, uintE>;
  long n = GA.nv + GA.nu;
  //bool* active = newA(bool,n);
  //{parallel_for(long i=0;i<n;i++) active[i] = 1;}
  //vertexSubset Frontier(n, n, active);
  /*uintE* Degrees = newA(uintE,n);
    {parallel_for(long i=0;i<GA.nv;i++) {
    Degrees[i] = GA.offsetsV[i+1] - GA.offsetsV[i];
    }}
    {parallel_for(long i=0;i<GA.nu;i++) {
    Degrees[GA.nv+i] = GA.offsetsU[i+1] - GA.offsetsU[i];
    }}*/
  //bool* Flags = newA(bool,n);
  //{parallel_for(long i=0;i<n;i++) Flags[i] = 0;}
  auto D = array_imap<uintE>(n, [&] (size_t i) {
      if(i >= GA.nv) return GA.offsetsU[i-GA.nv+1] - GA.offsetsU[i-GA.nv] <= 0 ? 0 : (uintE) floor(log2(GA.offsetsU[i-GA.nv+1] - GA.offsetsU[i-GA.nv])) + 1;
      return GA.offsetsV[i+1] - GA.offsetsV[i] <= 0 ? 0 : (uintE) floor(log2(GA.offsetsV[i+1] - GA.offsetsV[i])) + 1;
    });

  auto D_act = array_imap<uintE>(n, [&] (size_t i) {
      if(i >= GA.nv) return GA.offsetsU[i-GA.nv+1] - GA.offsetsU[i-GA.nv];
      return GA.offsetsV[i+1] - GA.offsetsV[i];
    });

  auto hp = BipartiteProp<uintE>(GA, make_tuple(UINT_E_MAX, 0), (size_t)GA.numEdges/5);
  auto b = make_buckets(n, D, decreasing, num_buckets);

  size_t finished = 0;
  while (finished != n) {
    auto bkt = b.next_bucket();
    auto active2 = bkt.identifiers;
    uintE k = bkt.id;
    finished += active2.size();
    auto apply_f = [&] (const tuple<uintE, uintE>& p) -> const Maybe<tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      uintE deg = D_act.s[v];
      uintE deg_log = D.s[v];
      if (deg_log < k) {
        uintE new_deg = deg - edgesRemoved; //min(deg - edgesRemoved, k);
        uintE new_deg_log = deg - edgesRemoved <= 0 ? 0 : (uintE) floor(log2(new_deg)) + 1;
        D_act.s[v] = new_deg;
        D.s[v] = min(new_deg_log, k);
        uintE bkt = b.get_bucket(deg_log, new_deg_log);
        return wrap(v, bkt);
      }
      return Maybe<tuple<uintE, uintE> >();
    };

    vertexSubsetData<uintE> moved = hp.template bpedgePropCount<uintE>(active2, apply_f);
    b.update_buckets(moved.get_fn_repr(), moved.size());
    moved.del(); active2.del();
  }
  b.del();
  
  auto samplesort_f = [&] (const uintE a, const uintE b) -> const uintE {
    return D[a] > D[b];
  };
  //free(active);
  return getRanks(GA, samplesort_f);
}

tuple<uintE*,uintE*,uintE*> getCoCoreRanks(bipartiteCSR& GA, size_t num_buckets=128) {
  using X = tuple<bool, uintE>;
  long n = GA.nv + GA.nu;
  //bool* active = newA(bool,n);
  //{parallel_for(long i=0;i<n;i++) active[i] = 1;}
  //vertexSubset Frontier(n, n, active);
  /*uintE* Degrees = newA(uintE,n);
    {parallel_for(long i=0;i<GA.nv;i++) {
    Degrees[i] = GA.offsetsV[i+1] - GA.offsetsV[i];
    }}
    {parallel_for(long i=0;i<GA.nu;i++) {
    Degrees[GA.nv+i] = GA.offsetsU[i+1] - GA.offsetsU[i];
    }}*/

  auto D = array_imap<uintE>(n, [&] (size_t i) {
      if(i >= GA.nv) return GA.offsetsU[i-GA.nv+1] - GA.offsetsU[i-GA.nv];
      return GA.offsetsV[i+1] - GA.offsetsV[i];
    });

  auto hp = BipartiteProp<uintE>(GA, make_tuple(UINT_E_MAX, 0), (size_t)GA.numEdges/5);
  auto b = make_buckets(n, D, decreasing, num_buckets);

  size_t finished = 0;
  //long rounds = 0;
  //long small_rounds = 0;
  while (finished != n) {
    //rounds++;
    auto bkt = b.next_bucket();
    auto active2 = bkt.identifiers;
    uintE k = bkt.id;
    finished += active2.size();
    /*if (active.size() <= 10) {
      small_rounds++;
      for (intT i=0; i < active.size(); ++i) {
      bool use_v = active.vtx(i) < GA.nv;
      intT idx = use_v ? active.vtx(i) : active.vtx(i) - GA.nv;
      intT offset  = use_v ? GA.offsetsV[idx] : GA.offsetsU[idx];
      intT deg = (use_v ? GA.offsetsV[idx+1] : GA.offsetsU[idx+1]) - offset;
      X* updated = newA(X, deg);
      intT deg_idx = 0;
      granular_for(j,0,deg,deg > 10000, { 
      intT nbhr = use_v ? GA.edgesV[offset+j] + GA.nv : GA.edgesU[offset+j];
      // must decrement D.s[nbhr]
      //writeAdd(&update[nbhr], 1);
      uintE old_deg = D.s[nbhr];
      if (old_deg < k) {
      uintE new_deg = min(old_deg - 1, k);
      D.s[nbhr] = new_deg;
      uintE bkt = b.get_bucket(old_deg, new_deg);
      updated[deg_idx++] = make_tuple(nbhr, bkt);
      }
      });
      vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(GA.nu+GA.nv,deg_idx,updated);
      b.update_buckets(moved.get_fn_repr(), moved.size());
      moved.del();
      }
      active.del();
      }
      else {*/
    auto apply_f = [&] (const tuple<uintE, uintE>& p) -> const Maybe<tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      uintE deg = D.s[v];
      if (deg < k) {
        uintE new_deg = min(deg - edgesRemoved, k);
        D.s[v] = new_deg;
        uintE bkt = b.get_bucket(deg, new_deg);
        return wrap(v, bkt);
      }
      return Maybe<tuple<uintE, uintE> >();
    };

    //vertexSubset FrontierH = vertexProp(GA,active,Remove_BPedge(Flags));
    //cout << "k="<<k<< " num active = " << active.numNonzeros() << " frontierH = " << FrontierH.numNonzeros() << endl;
    vertexSubsetData<uintE> moved = hp.template bpedgePropCount<uintE>(active2, apply_f);
    b.update_buckets(moved.get_fn_repr(), moved.size());
    moved.del(); active2.del();
    //}
  }
  b.del();
  //free(active);

  auto samplesort_f = [&] (const uintE a, const uintE b) -> const uintE {
    return D[a] > D[b];
  };
  return getRanks(GA, samplesort_f);
}

tuple<uintE*, uintE*, uintE*> getCoreRanks(bipartiteCSR& G, size_t num_buckets=128) {
  using X = tuple<uintE, uintE>;
  const size_t n = G.nv + G.nu; const size_t m = G.numEdges;
  auto D = array_imap<uintE>(n, [&] (size_t i) {
      if(i >= G.nv) return G.offsetsU[i-G.nv+1] - G.offsetsU[i-G.nv];
      return G.offsetsV[i+1] - G.offsetsV[i];
    });

  //auto em = EdgeMap<uintE, vertex>(GA, make_tuple(UINT_E_MAX, 0), (size_t)G.numEdges/5);
  auto b = make_buckets(n, D, increasing, num_buckets);
  //intT* update = newA(intT, G.nu + G.nv);
  //parallel_for(long v=0; v < G.nv+G.nu; ++v) { update[v] = 0; }

  size_t finished = 0;
  while (finished != n) {
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();

    long num_hash = sequence::reduce<long>((long) 0, active.size(), addF<long>(), degF<long>(G.nv, G.offsetsV, G.offsetsU, active));
    sparseAdditiveSet<intT> update_hash = sparseAdditiveSet<intT>(num_hash, 1, INT_T_MAX);

    parallel_for(intT i=0; i < active.size(); ++i) {
      bool use_v = active.vtx(i) < G.nv;
      intT idx = use_v ? active.vtx(i) : active.vtx(i) - G.nv;
      intT offset  = use_v ? G.offsetsV[idx] : G.offsetsU[idx];
      intT deg = (use_v ? G.offsetsV[idx+1] : G.offsetsU[idx+1]) - offset;
      granular_for(j,0,deg,deg > 10000, { 
	  intT nbhr = use_v ? G.edgesV[offset+j] + G.nv : G.edgesU[offset+j];
	  // must decrement D.s[nbhr]
	  //writeAdd(&update[nbhr], 1);
	  update_hash.insert(make_pair(nbhr,1));
	});
    }

    auto update = update_hash.entries();
    X* update_b = newA(X, update.n);
    parallel_for(long i=0; i < update.n; ++i) { 
      uintE v = update.A[i].first;
      uintE deg = D.s[v];
      if (deg > k) {
        uintE new_deg = max(deg - update.A[i].second, k);
        D.s[v] = new_deg;
        uintE bkt = b.get_bucket(deg, new_deg);
        update_b[i] = make_tuple(v, bkt);
      }
      else update_b[i] = make_tuple(UINT_E_MAX, UINT_E_MAX);
    }
    update_hash.del();
    update.del();

    X* update_filter = newA(X, update.n);
    long num_updates_filter = sequence::filter(update_b ,update_filter, update.n, nonMaxTupleF());
    free(update_b);

    //vertexSubsetData<uintE> moved = em.template edgeMapCount<uintE>(active, apply_f);
    // second should be array of tuple<uintE,uintE> of idx, new bucket pairs; first is length of this array
    vertexSubsetData<uintE> moved = vertexSubsetData<uintE>(G.nu+G.nv, num_updates_filter, update_filter);
    b.update_buckets(moved.get_fn_repr(), moved.size());
    moved.del(); 
    active.del();
  }

  auto samplesort_f = [&] (const uintE a, const uintE b) -> const uintE {
    return D[a] > D[b];
  };
  
  return getRanks(G, samplesort_f);
}

graphCSR rankGraph(bipartiteCSR& G, bool use_vb, uintE* ranks, uintE* rankV, uintE* rankU) {
  using X = tuple<uintE,uintE>;
  // we put a 1 if nu and use_vb; 0 otherwise
  // store if 1, don't store if 0 (store if nu and use_v or if nv and not use_v)
  
  uintT* offsets = newA(uintT,G.nv+G.nu+1);
  offsets[G.nv+G.nu] = 0;

  parallel_for(long i=0; i < G.nv + G.nu; ++i) {
    if (ranks[i] >= G.nv) offsets[i] = G.offsetsU[ranks[i]-G.nv+1]-G.offsetsU[ranks[i]-G.nv];
    else offsets[i] = G.offsetsV[ranks[i]+1]-G.offsetsV[ranks[i]];
  }

  // Now we have to reformat the graph
  sequence::plusScan(offsets,offsets,G.nv+G.nu+1);

  uintE* edges = newA(uintE,offsets[G.nv+G.nu]);

  auto lt = [] (const uintE& l, const uintE& r) { return l > r; };
  parallel_for(long i=0; i < G.nv + G.nu; ++i) {
    // need to fill in offsets[i] to offsets[i+1] in edges array
    bool use_v = (ranks[i] < G.nv);
    intT idx = use_v ? ranks[i] : ranks[i] - G.nv;
    intT offset  = use_v ? G.offsetsV[idx] : G.offsetsU[idx];
    intT deg = (use_v ? G.offsetsV[idx+1] : G.offsetsU[idx+1])-offset;
    granular_for(j,0,deg,deg > 10000, { 
	intT nbhr = use_v ? G.edgesV[offset+j] : G.edgesU[offset+j];
	uintE r; 
	if (use_vb) r = use_v ? (rankU[nbhr] << 1) + 0b1  : (rankV[nbhr] << 1);
	else r = use_v ? (rankU[nbhr] << 1) : (rankV[nbhr] << 1) + 0b1;
	edges[offsets[i]+j] = r;
      });
    sampleSort(&edges[offsets[i]], deg, lt);
  }
  return graphCSR(offsets,edges,G.nv+G.nu,G.numEdges);
}

pair<graphCSR,tuple<uintE,uintE>*> rankGraphEdges(bipartiteCSR& G, bool use_vb, uintE* ranks, uintE* rankV, uintE* rankU) {
  using X = tuple<uintE,uintE>;
  // we put a 1 if nu and use_vb; 0 otherwise
  // store if 1, don't store if 0 (store if nu and use_v or if nv and not use_v)

  uintT* offsets = newA(uintT,G.nv+G.nu+1);
  offsets[G.nv+G.nu] = 0;

  parallel_for(long i=0; i < G.nv + G.nu; ++i) {
    if (ranks[i] >= G.nv) offsets[i] = G.offsetsU[ranks[i]-G.nv+1]-G.offsetsU[ranks[i]-G.nv];
    else offsets[i] = G.offsetsV[ranks[i]+1]-G.offsetsV[ranks[i]];
  }

  // Now we have to reformat the graph
  sequence::plusScan(offsets,offsets,G.nv+G.nu+1);

  uintE* edges = newA(uintE,offsets[G.nv+G.nu]);
  X* edges_convert = newA(X,offsets[G.nv+G.nu]);

  auto lt = [] (const uintE& l, const uintE& r) { return l > r; };
  auto ltx = [] (const X& l, const X& r) { return get<0>(l) > get<0>(r); };
  parallel_for(long i=0; i < G.nv + G.nu; ++i) {
    // need to fill in offsets[i] to offsets[i+1] in edges array
    bool use_v = (ranks[i] < G.nv);
    intT idx = use_v ? ranks[i] : ranks[i] - G.nv;
    intT offset  = use_v ? G.offsetsV[idx] : G.offsetsU[idx];
    intT deg = (use_v ? G.offsetsV[idx+1] : G.offsetsU[idx+1])-offset;
    granular_for(j,0,deg,deg > 10000, { 
	intT nbhr = use_v ? G.edgesV[offset+j] : G.edgesU[offset+j];
	uintE r; 
	if (use_vb) r = use_v ? (rankU[nbhr] << 1) + 0b1  : (rankV[nbhr] << 1);
	else r = use_v ? (rankU[nbhr] << 1) : (rankV[nbhr] << 1) + 0b1;
	//edges[offsets[i]+j] = r;
	edges_convert[offsets[i]+j] = make_tuple(r,offset+j);
      });
    sampleSort(&edges_convert[offsets[i]], deg, ltx);
    granular_for(j,0,deg,deg > 10000, {
	edges[offsets[i]+j] = get<0>(edges_convert[offsets[i]+j]); 
      });
  }
  return make_pair(graphCSR(offsets,edges,G.nv+G.nu,G.numEdges), edges_convert);
}

template <class E>
struct clrF { 
  uintE* edges; uintT offset; uintT color; uintT* colors;
  clrF(uintE* _edges, uintT _offset, uintT _color, uintT* _colors) : edges(_edges), offset(_offset), color(_color), colors(_colors) {}
  inline E operator() (const uintT& i) const {
    return (E) (colors[edges[offset+i]] == color ? 1 : 0); 
  }
};

template <class E>
struct eF { 
  uintT offset; uintE* colors;
  eF(uintT _offset, uintE* _colors) : offset(_offset), colors(_colors) {}
  inline E operator() (const uintT& i) const {
    return (E) (colors[offset+i] == 0);
  }
};

struct nonMaxUintTF{bool operator() (uintT &a) {return (a != UINT_T_MAX);}};

struct isSameColor{
  isSameColor(uintT mycolor, uintT *colors) : colors_(colors), me_(mycolor) {};
  bool operator () (uintT v) {return me_ == colors_[v];};
  uintT me_;
  uintT *colors_;
};

struct isZeroF{
  isZeroF(uintT offset, uintE *colors) : colors_(colors), offset_(offset) {};
  bool operator () (uintT v) {return colors_[offset_+v] == 0;};
  uintT offset_;
  uintT *colors_;
};

bipartiteCSR eSparseBipartite(bipartiteCSR& G, long denom, long seed) {
  uintE numColors = max<uintE>(1,denom);
  uintE* colorsV = newA(uintE,G.numEdges);
  uintE* colorsU = newA(uintE,G.numEdges);
  parallel_for(long i=0; i < G.nv; ++i) {
    uintT v_offset = G.offsetsV[i];
    uintT v_deg = G.offsetsV[i+1] - v_offset;
    parallel_for(long j=0; j < v_deg; ++j) {
      uintT u = G.edgesV[v_offset + j];
      long concat = i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      concat ^= u + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      //long concat = (i + u) * (i + u + 1) / 2 + i;
      colorsV[v_offset + j] = hashInt((ulong) seed + concat) % numColors;
    }
  }
  parallel_for(intT i=0; i < G.nu; ++i) {
    intT u_offset = G.offsetsU[i];
    intT u_deg = G.offsetsU[i+1] - u_offset;
    granular_for(j,0,u_deg,u_deg>1000,{
      uintE v = G.edgesU[u_offset + j];
      intT v_offset = G.offsetsV[v];
      intT v_deg = G.offsetsV[v+1] - v_offset;
      // find k such that edgesV[v_offset + k] = i
      auto idx_map = make_in_imap<uintE>(v_deg, [&] (size_t k) { return G.edgesV[v_offset + k]; });
      auto lte = [] (const uintE& l, const uintE& r) { return l < r; };
      size_t find_idx = pbbs::binary_search(idx_map, i, lte);
      colorsU[u_offset + j] = colorsV[v_offset + find_idx];
      });
  }
  /*parallel_for(long u=0; u < G.nu; ++u) {
    uintT u_offset = G.offsetsU[u];
    uintT u_deg = G.offsetsU[u+1] - u_offset;
    parallel_for(long j=0; j < u_deg; ++j) {
      uintT i = G.edgesU[u_offset + j];
      long concat = i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      concat ^= u + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      //long concat = (i + u) * (i + u + 1) / 2 + i;
      colorsU[v_offset + j] = hashInt((ulong) seed + concat) % numColors;
    }
  }*/
  uintT* offsetsV = newA(uintT,G.nv+1);
  uintT* offsetsU = newA(uintT,G.nu+1);
  offsetsV[G.nv] = 0; offsetsU[G.nu] = 0;
  parallel_for(long i=0; i < G.nv; ++i) {
    uintT v_offset = G.offsetsV[i];
    uintT v_deg = G.offsetsV[i+1] - v_offset;
    if (v_deg > 10000) offsetsV[i] = sequence::reduce<uintT>((uintT) 0, v_deg, addF<uintT>(), eF<uintT>(v_offset, colorsV));
    else {
      offsetsV[i] = 0;
      for(long j=0; j < v_deg; ++j) {
        if (colorsV[v_offset + j] == 0) offsetsV[i]++;
      }
    }
  }
  parallel_for(long i=0; i < G.nu; ++i) {
    uintT u_offset = G.offsetsU[i];
    uintT u_deg = G.offsetsU[i+1] - u_offset;
    if (u_deg > 10000) offsetsU[i] = sequence::reduce<uintT>((uintT) 0, u_deg, addF<uintT>(), eF<uintT>(u_offset, colorsU));
    else {
      offsetsU[i] = 0;
      for(long j=0; j < u_deg; ++j) {
        if (colorsU[u_offset + j] == 0) offsetsU[i]++;
      }
    }
  }
  long mv = sequence::plusScan(offsetsV,offsetsV,G.nv+1);
  long mu = sequence::plusScan(offsetsU,offsetsU,G.nu+1);
  uintE* edgesV = newA(uintE,mv);
  uintE* edgesU = newA(uintE,mu);
  parallel_for(long i=0; i<G.nv; ++i) {
    uintT v_offset = G.offsetsV[i];
    uintT v_deg = G.offsetsV[i+1] - v_offset;
    uintT v_clr_offset = offsetsV[i];
    if (v_deg > 10000) sequence::filter(&G.edgesV[v_offset],&edgesV[v_clr_offset],v_deg,isZeroF(v_offset, colorsV));
    else {
      long idx = 0;
      for(long j=0; j < v_deg; ++j) {
        uintE u = G.edgesU[v_offset + j];
        if (colorsV[v_offset + j] == 0) {edgesV[v_clr_offset + idx] = u; idx++;}
      }
    }
  }
  parallel_for(long i=0; i<G.nu; ++i) {
    uintT v_offset = G.offsetsU[i];
    uintT v_deg = G.offsetsU[i+1] - v_offset;
    uintT v_clr_offset = offsetsU[i];
    if (v_deg > 10000) sequence::filter(&G.edgesU[v_offset],&edgesU[v_clr_offset],v_deg,isZeroF(v_offset, colorsU));
    else {
      long idx = 0;
      for(long j=0; j < v_deg; ++j) {
        uintE u = G.edgesU[v_offset + j];
        if (colorsU[v_offset + j] == 0) {edgesU[v_clr_offset + idx] = u; idx++;}
      }
    }
  }
  free(colorsV); free(colorsU);
  uintT* offsetsV_f = newA(uintT,G.nv+1);
  uintT* offsetsU_f = newA(uintT,G.nu+1);
  parallel_for(long i=0; i < G.nv; ++i) {if (offsetsV[i] == offsetsV[i+1]) offsetsV_f[i] = UINT_T_MAX; else offsetsV_f[i] = i;}
  parallel_for(long i=0; i < G.nu; ++i) {if (offsetsU[i] == offsetsU[i+1]) offsetsU_f[i] = UINT_T_MAX; else offsetsU_f[i] = i;}
  uintT* offsetsV_ff = newA(uintT,G.nv+1);
  uintT* offsetsU_ff = newA(uintT,G.nu+1);
  long num_vff = sequence::filter(offsetsV_f,offsetsV_ff,G.nv,nonMaxUintTF());
  long num_uff = sequence::filter(offsetsU_f,offsetsU_ff,G.nu,nonMaxUintTF());
  parallel_for(long i=0; i < num_vff; ++i) {offsetsV_f[offsetsV_ff[i]] = i;}
  parallel_for(long i=0; i < num_uff; ++i) {offsetsU_f[offsetsU_ff[i]] = i;}
  parallel_for(long i=0; i < mv; ++i) { edgesV[i] = offsetsU_f[edgesV[i]]; }
  parallel_for(long i=0; i < mv; ++i) { edgesU[i] = offsetsV_f[edgesU[i]]; }

  parallel_for(long i=0; i < G.nv; ++i) {if (offsetsV[i] == offsetsV[i+1]) offsetsV_f[i] = UINT_T_MAX; else offsetsV_f[i] = offsetsV[i];}
  parallel_for(long i=0; i < G.nu; ++i) {if (offsetsU[i] == offsetsU[i+1]) offsetsU_f[i] = UINT_T_MAX; else offsetsU_f[i] = offsetsV[i];}
  offsetsV_f[G.nv] = offsetsV[G.nv]; offsetsU_f[G.nu] = offsetsU[G.nu];
  num_vff = sequence::filter(offsetsV_f,offsetsV_ff,G.nv+1,nonMaxUintTF());
  num_uff = sequence::filter(offsetsU_f,offsetsU_ff,G.nu+1,nonMaxUintTF());
  free(offsetsV_f); free(offsetsU_f); free(offsetsV); free(offsetsU);
  return bipartiteCSR(offsetsV_ff,offsetsU_ff,edgesV,edgesU,num_vff-1,num_uff-1,mv);
}

bipartiteCSR clrSparseBipartite(bipartiteCSR& G, long denom, long seed) {
  double p = 1/denom;
  uintT numColors = max<uintT>(1,denom);
  uintT* colorsV = newA(uintT,G.nv);
  uintT* colorsU = newA(uintT, G.nu);
  parallel_for(long i=0;i<G.nv;i++) colorsV[i] = hashInt((ulong) seed+i) % numColors;
  parallel_for(long i=0;i<G.nu;i++) colorsU[i] = hashInt((ulong) seed+i) % numColors;

  uintT* offsetsV = newA(uintT,G.nv+1);
  uintT* offsetsU = newA(uintT,G.nu+1);
  offsetsV[G.nv] = 0; offsetsU[G.nu] = 0;
  parallel_for(long i=0; i < G.nv; ++i) {
    uintT v_offset = G.offsetsV[i];
    uintT v_deg = G.offsetsV[i+1] - v_offset;
    if (v_deg > 10000) offsetsV[i] = sequence::reduce<uintT>((uintT) 0, v_deg, addF<uintT>(), clrF<uintT>(G.edgesV, v_offset, colorsV[i], colorsU));
    else {
      offsetsV[i] = 0;
      for(long j=0; j < v_deg; ++j) {
        if (colorsU[G.edgesV[v_offset + j]] == colorsV[i]) offsetsV[i]++;
      }
    }
  }
  parallel_for(long i=0; i < G.nu; ++i) {
    uintT v_offset = G.offsetsU[i];
    uintT v_deg = G.offsetsU[i+1] - v_offset;
    if (v_deg > 10000) offsetsU[i] = sequence::reduce<uintT>((uintT) 0, v_deg, addF<uintT>(), clrF<uintT>(G.edgesU, v_offset, colorsU[i], colorsV));
    else {
      offsetsU[i] = 0;
      for(long j=0; j < v_deg; ++j) {
        uintE u = G.edgesU[v_offset + j];
        if (colorsV[u] == colorsU[i]) offsetsU[i]++;
      }
    }
  }
  long mv = sequence::plusScan(offsetsV,offsetsV,G.nv+1);
  long mu = sequence::plusScan(offsetsU,offsetsU,G.nu+1);
  uintE* edgesV = newA(uintE,mv);
  uintE* edgesU = newA(uintE,mu);
  parallel_for(long i=0; i<G.nv; ++i) {
    uintT v_offset = G.offsetsV[i];
    uintT v_deg = G.offsetsV[i+1] - v_offset;
    uintT v_clr_offset = offsetsV[i];
    if (v_deg > 10000) sequence::filter(&G.edgesV[v_offset],&edgesV[v_clr_offset],v_deg,isSameColor(colorsV[i],colorsU));
    else {
      long idx = 0;
      for(long j=0; j < v_deg; ++j) {
        uintE u = G.edgesV[v_offset + j];
        if (colorsU[u] == colorsV[i]) {edgesV[v_clr_offset + idx] = u; idx++;}
      }
    }
  }
  parallel_for(long i=0; i<G.nu; ++i) {
    uintT v_offset = G.offsetsU[i];
    uintT v_deg = G.offsetsU[i+1] - v_offset;
    uintT v_clr_offset = offsetsU[i];
    if (v_deg > 10000) sequence::filter(&G.edgesU[v_offset],&edgesU[v_clr_offset],v_deg,isSameColor(colorsU[i],colorsV));
    else {
      long idx = 0;
      for(long j=0; j < v_deg; ++j) {
        uintE u = G.edgesU[v_offset + j];
        if (colorsV[u] == colorsU[i]) {edgesU[v_clr_offset + idx] = u; idx++;}
      }
    }
  }
  free(colorsV); free(colorsU);
  return bipartiteCSR(offsetsV,offsetsU,edgesV,edgesU,G.nv,G.nu,mv);  
}

bipartiteCSR readBipartite(char* fname) {
  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);

  if (W.Strings[0] != (string) "AdjacencyHypergraph") {
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long nv = atol(W.Strings[1]);
  long mv = atol(W.Strings[2]);
  long nu = atol(W.Strings[3]);
  long mu = atol(W.Strings[4]);

  if ((len != nv + mv + nu + mu + 4) | (mv != mu)) {
    cout << "Bad input file" << endl;
    abort();
  }

  uintT* offsetsV = newA(uintT,nv+1);
  uintT* offsetsU = newA(uintT,nu+1);
  uintE* edgesV = newA(uintE,mv);
  uintE* edgesU = newA(uintE,mu);

  {parallel_for(long i=0; i < nv; i++) offsetsV[i] = atol(W.Strings[i + 5]);}
  offsetsV[nv] = mv;
  
  {parallel_for(long i=0; i<mv; i++) {
      edgesV[i] = atol(W.Strings[i+nv+5]);
      if(edgesV[i] < 0 || edgesV[i] >= nu) { cout << "edgesV out of range: nu = " << nu << " edge = " << edgesV[i] << endl; exit(0); }
    }}

  {parallel_for(long i=0; i < nu; i++) offsetsU[i] = atol(W.Strings[i + nv + mv + 5]);}
  offsetsU[nu] = mu;
  
  {parallel_for(long i=0; i<mu; i++) {
      edgesU[i] = atol(W.Strings[i+nv+mv+nu+5]);
      if(edgesU[i] < 0 || edgesU[i] >= nv) { cout << "edgesU out of range: nv = " << nv << " edge = " << edgesU[i] << endl; exit(0); }
    }}

  S.del();
  free(W.Strings);
  //W.del(); // to deal with performance bug in malloc
  return bipartiteCSR(offsetsV,offsetsU,edgesV,edgesU,nv,nu,mv);  
}

// Takes the elements of a vertex array, and returns the out degree choose 2
template <class E>
struct wedgeF { 
  uintT* offsets;
wedgeF(uintT* _offsets) : offsets(_offsets) {}
  inline E operator() (const uintT& i) const {
    uintE v_deg = offsets[i+1]-offsets[i];
    return (E) ((v_deg * (v_deg-1)) / 2); 
  }
};

long* computeWorkPrefixSum(bipartiteCSR & GA, bool use_v) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long* workPrefixSum;
  workPrefixSum = newA(long,nu);
  parallel_for(intT i=0;i<nu;i++) {
    long u_work = 0;
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1]-offsetsU[i];
    for(intT j=0; j<u_deg;j++) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      // for(intT k=0; k<v_deg;k++) {
      //   uintE u2_idx = GA.edgesV[v_offset+k];
      //   if(u2_idx < i) u_work++;
      //   else break;
  
      // }
      u_work += offsetsV[v+1]-offsetsV[v];
    }
    workPrefixSum[i] = u_work;
  }
  sequence::plusScan(workPrefixSum,workPrefixSum,nu);
  return workPrefixSum;
}

tuple<bool,long> cmpWedgeCounts(bipartiteCSR & GA, long type=0) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = sequence::reduce<long>((long) 0, nv, addF<long>(), wedgeF<long>(GA.offsetsV));
  long num_wedges_u = sequence::reduce<long>((long) 0, nu, addF<long>(), wedgeF<long>(GA.offsetsU));

  return make_tuple((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
}

pair<bool,long> cmpWedgeCountsSeq(bipartiteCSR & GA) {
  const long nv = GA.nv, nu = GA.nu;

  long num_wedges_v = 0;
  for(long i=0; i < nv; ++i) {
    uintE deg_v = GA.offsetsV[i+1]-GA.offsetsV[i];
    num_wedges_v += deg_v * (deg_v - 1) / 2;
  }

  long num_wedges_u = 0;
  for(long i=0; i < nu; ++i) {
    uintE deg_u = GA.offsetsU[i+1]-GA.offsetsU[i];
    num_wedges_u += deg_u * (deg_u - 1) / 2;
  }
  return make_pair((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
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
// TODO turn cmp into extract
struct nonMaxLongF{bool operator() (long &a) {return (a != LONG_MAX);}};

template <class L, class T, class Cmp, class Eq, class F>
  pair<L*, long> getFreqs(T* objs, long num, Cmp cmp, Eq eq, L maxL, F nonF, bool sort=true) {
  // Sort objects
  //if (sort) parallelIntegerSort(objs, num, cmp); //sampleSort(objs, num, cmp);

  L* freqs = newA(L, num + 1);
  freqs[0] = 0;
  freqs[num] = num;

  // Retrieve indices where objects differ
  parallel_for(long i=1; i < num; ++i) {
    if (!eq(objs[i-1],objs[i])) freqs[i] = i;
    else freqs[i] = maxL;
  }
  L* freqs_f = newA(L, num+1);
  long num_freqs_f = sequence::filter(freqs, freqs_f, num+1, nonF);
  free(freqs);
  return make_pair(freqs_f, num_freqs_f);
}

template <class L, class S, class T, class Cmp, class Eq, class OpT, class OpuintE, class OpCount>
  pair<tuple<S,L>*, long> getFreqs_seq(T* objs, long num, Cmp cmp, Eq eq, bool sort=true, OpT opt=refl<T>(),
				       OpuintE opuinte=refl<L>(), OpCount opcount=reflCount<T>()) {
  //if(sort) parallelIntegerSort(objs, num, cmp); //sampleSort(objs, num, cmp);

  using X = tuple<S,L>;
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

long* countWedgesScan(graphCSR& G) {
  long* idxs = newA(long, G.n + 1);
  idxs[G.n] = 0;

  using T = long*;
  T* nbhd_idxs = newA(T, G.n);

  parallel_for(intT i=0; i < G.n; ++i) {
    idxs[i] = 0;
    intT offset = G.offsets[i];
    intT deg = G.offsets[i+1] - offset;

    nbhd_idxs[i] = newA(long, deg + 1);
    (nbhd_idxs[i])[deg] = 0;

    parallel_for(intT j=0; j < deg; ++j) {
      (nbhd_idxs[i])[j] = 0;
      uintE v = G.edges[offset+j] >> 1;
      intT v_offset = G.offsets[v];
      intT v_deg = G.offsets[v+1] - v_offset;
      if (v > i) {
	for (intT k = 0; k < v_deg; ++k) {
	  if ((G.edges[v_offset + k] >> 1) > i) (nbhd_idxs[i][j])++;
	  else break;
	}
      }
    }
    idxs[i] = sequence::plusReduce(nbhd_idxs[i], deg + 1);
    free(nbhd_idxs[i]);
  }
  free(nbhd_idxs);

  sequence::plusScan(idxs, idxs, G.n + 1);
  return idxs;
}

//TODO we don't use nbhrs or save nbhrs -- delete from everything
long* countWedgesScan(bipartiteCSR& GA, bool use_v, bool half=false) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long* idxs = newA(long, nu + 1);
  idxs[nu] = 0;

  using T = long*;
  T* nbhd_idxs = newA(T, nu);

  parallel_for(intT i=0; i < nu; ++i) {
    idxs[i] = 0;
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;

    nbhd_idxs[i] = newA(long, u_deg + 1);
    parallel_for(intT j=0; j < u_deg+1; ++j) {(nbhd_idxs[i])[j] = 0;} //JS: granular_for

    parallel_for(intT j=0; j < u_deg; ++j) { //TODO can parallelize this too technically //JS: granular_for
      if (!half) {
        uintE v = edgesU[u_offset + j];
        (nbhd_idxs[i])[j] = offsetsV[v+1] - offsetsV[v] - 1;//V[U[i].getOutNeighbor(j)].getOutDegree() - 1;
      }
      else {
        (nbhd_idxs[i])[j] = 0;
        uintE v = edgesU[u_offset + j];
        intT v_offset = offsetsV[v];
        intT v_deg = offsetsV[v+1] - v_offset;
        for (intT k = 0; k < v_deg; ++k) { //TODO can parallelize this too technically
          if (edgesV[v_offset + k] < i) nbhd_idxs[i][j] ++;
          else break;
        }
      }
    }

    idxs[i] = sequence::plusReduce(nbhd_idxs[i], u_deg + 1);
    free(nbhd_idxs[i]);
  }
  free(nbhd_idxs);
  sequence::plusScan(idxs, idxs, nu + 1);

  return idxs;
}

struct CountESpace {
  long type; long nu; bool rank;
  // for sort: 0, 1
  _seq<UWedge> wedges_seq_uw;
  // for sort: 1; for hist: 6
  _seq<tuple<uintE,long>> butterflies_seq_intt;
  // for hist: 6
  pbbsa::sequence<tuple<uintE, long>> tmp_uint;
  pbbsa::sequence<tuple<uintE, long>> out_uint;
  // for hist: 4
  pbbsa::sequence<tuple<long, uintE>> tmp;
  pbbsa::sequence<tuple<long, uintE>> out;
  _seq<long> wedges_seq_int;
  // for hist: 4, for hash: 2, 3
  sparseAdditiveSet<long, long> wedges_hash;
  // for hash: 3
  _seq<pair<long, long>> wedges_seq_intp;
  sparseAdditiveSet<long, long> butterflies_hash;


CountESpace(long _type, long _nu, bool _rank) : type(_type), nu(_nu), rank(_rank) {
  using T = pair<long,long>;
  using X = tuple<uintE,long>;
  using E = pair<long, uintE>;
  if (type == 2 || type == 3) {
    wedges_hash = sparseAdditiveSet<long, long>(nu,1,LONG_MAX, LONG_MAX);
    if (type == 3) {
      wedges_seq_intp = _seq<T>(newA(T, nu), nu);
      butterflies_hash = sparseAdditiveSet<long, long>(nu, 1, LONG_MAX, LONG_MAX);
    }
  }
  else if (type == 4 || type == 6) {
    tmp = pbbsa::sequence<tuple<long, uintE>>();
    out = pbbsa::sequence<tuple<long, uintE>>();
    wedges_seq_int = _seq<long>(newA(long, nu), nu);
    wedges_hash = sparseAdditiveSet<long, long>(nu,1,LONG_MAX, LONG_MAX);
    if (type == 6) {
      butterflies_seq_intt = _seq<X>(newA(X, 1), 1);
      tmp_uint = pbbsa::sequence<tuple<uintE, long>>();
      out_uint = pbbsa::sequence<tuple<uintE, long>>();
    }
  }
  else if (type == 0 || type == 1) {
    if (type == 1) butterflies_seq_intt = _seq<X>(newA(X, 1), 1);
    wedges_seq_uw = _seq<UWedge>(newA(UWedge, nu), nu);
  }
}

  void clear() {
    if (type == 2 || type == 3 || type == 4 || type == 6) wedges_hash.clear();
    if (type == 3) butterflies_hash.clear();
  }

  void del() {
    if (type == 2 || type == 3) {
      wedges_hash.del();
      if (type == 3) { wedges_seq_intp.del(); butterflies_hash.del(); }
    }
    else if (type == 4 || type == 6) {
      wedges_seq_int.del();
      wedges_hash.del();
      if (type == 6) butterflies_seq_intt.del();
    }
    else if (type == 0 || type == 1) {
      if (type == 1) butterflies_seq_intt.del();
      wedges_seq_uw.del();
    }
  }
};

struct CountSpace {
  long type; long nu; bool rank;
  // for hash: 2, 3
  sparseAdditiveSet<long, long> wedges_hash;
  _seq<pair<long,long>> wedges_seq_intp;
  // for hash: 3
  sparseAdditiveSet<long, long> butterflies_hash;
  _seq<pair<long,long>> butterflies_seq_intp;
  // for hist: 4, 6
  pbbsa::sequence<tuple<long, uintE>> tmp;
  pbbsa::sequence<tuple<long, uintE>> out;
  _seq<long> wedges_seq_int;
  // for hist: 6
  _seq<tuple<uintE,long>> butterflies_seq_intt;
  pbbsa::sequence<tuple<uintE, long>> tmp_uint;
  pbbsa::sequence<tuple<uintE, long>> out_uint;
  // for sort: 0, 1
  _seq<UVertexPair> wedges_seq_uvp;
  _seq<UWedge> wedges_seq_uw;

CountSpace(long _type, long _nu, bool _rank) : type(_type), nu(_nu), rank(_rank) {
  using T = pair<uintE,long>;
  using X = tuple<uintE,uintE>;
  using E = pair<long,long>;
  using L = tuple<uintE,long>;
  if (type == 2 || type == 3) {
    wedges_hash = sparseAdditiveSet<long, long>(nu,1,LONG_MAX, LONG_MAX);
    wedges_seq_intp = _seq<E>(newA(E, nu), nu);
    if (type == 3) {
      butterflies_hash = sparseAdditiveSet<long, long>(nu, 1, LONG_MAX, LONG_MAX);
      butterflies_seq_intp = _seq<E>(newA(E, nu), nu);
    }
  }
  else if (type == 4 || type == 6) {
    tmp = pbbsa::sequence<tuple<long, uintE>>();
    out = pbbsa::sequence<tuple<long, uintE>>();
    wedges_seq_int = _seq<long>(newA(long, nu), nu);
    if (rank) {
      wedges_hash = sparseAdditiveSet<long, long>(nu,1,LONG_MAX, LONG_MAX);
      wedges_seq_intp = _seq<E>(newA(E, nu), nu);
    }
    if (type == 6) {
      butterflies_seq_intt = _seq<L>(newA(L, 1), 1);
      tmp_uint = pbbsa::sequence<tuple<uintE, long>>();
      out_uint = pbbsa::sequence<tuple<uintE, long>>();
    }
  }
  else if (type == 0 || type == 1) {
    if (type == 1) butterflies_seq_intt = _seq<L>(newA(L, 1), 1);
    if (!rank) wedges_seq_uvp = _seq<UVertexPair>(newA(UVertexPair, nu), nu);
    else wedges_seq_uw = _seq<UWedge>(newA(UWedge, nu), nu);
  }
}

  void clear() {
    if (type == 2 || type == 3) {
      wedges_hash.clear();
      if (type == 3) butterflies_hash.clear();
    }
    else if ((type == 4 || type == 6) && rank) wedges_hash.clear();

  }

  void del() {
    if (type == 2 || type == 3) {
      wedges_hash.del(); wedges_seq_intp.del(); 
      if (type == 3) { butterflies_hash.del(); butterflies_seq_intp.del(); }
    }
    else if (type == 4 || type == 6) {
      wedges_seq_int.del();
      if (rank){
        wedges_hash.del(); wedges_seq_intp.del();
      }
      if (type == 6) butterflies_seq_intt.del();
    }
    else if (type == 0 || type == 1) {
      if (type == 1) butterflies_seq_intt.del(); 
      if (!rank) wedges_seq_uvp.del();
      else wedges_seq_uw.del();
    }
  }
};

//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class wedgeCons>
  void _getWedges_seq(wedge* wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, intT curr_idx, intT next_idx) {
  const long nu = use_v ? GA.nu : GA.nv;
  const long nv = use_v ? GA.nv : GA.nu;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long idx = 0;
  for(intT i=curr_idx; i < next_idx; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset+k];
        if (u2 < i) {
          wedges[idx] = cons(i, u2, v, j, k); //TODO here: also store j, k; so we can do offsetsV[v] + k and eti[offsetsU[i] + j]
          ++idx;
        }
        else break; 
      }
    }
  }
}

template<class wedge, class wedgeCons>
  void _getWedges(_seq<wedge>& wedges_seq, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, long* wedge_idxs, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  // Allocate space for seagull storage
  if (wedges_seq.n < num_wedges) {
    free(wedges_seq.A);
    wedges_seq.A = newA(wedge, num_wedges);
    wedges_seq.n = num_wedges;
  }

  if (next_idx == INT_T_MAX) next_idx = nu;
  if (num_wedges < 10000) return _getWedges_seq<wedge>(wedges_seq.A, GA, use_v, cons, num_wedges, curr_idx, next_idx);
 
  // Store seagulls in parallel
  parallel_for(intT i=curr_idx; i < next_idx; ++i) {
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    // Consider each neighbor v of active vertex u
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    long idx = 0;
    for(intT j=0; j < u_deg; ++j) { //JS: test granular_for
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = edgesV[v_offset+k];
        if (u2 < i) {
          wedges_seq.A[wedge_idx+idx] = cons(i, u2, v, j, k);
          ++idx;
        }
        else {break;}
      }
    }
  }
}

template<class wedgeCons, class T>
  void _getWedgesHash(T& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long num_wedges, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  if (next_idx == INT_T_MAX) next_idx = nu;

  wedges.resize(num_wedges);

  //hashInsertTimer.start();
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    parallel_for(intT j=0;j<u_deg;j++) { //JS: test granular_for
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1] - v_offset;
	// Find all seagulls with center v and endpoint u
	for (long k=0; k < v_deg; ++k) { 
	  uintE u2 = edgesV[v_offset+k];
	  if (u2 < i) wedges.insert(make_pair(cons(i, u2, v, j, k),1));
	  else break;
	}
    }
  }
  //hashInsertTimer.stop();
}

intT getNextWedgeIdx_seq(bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx) {
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long orig = max_wedges;
  for(intT i=curr_idx; i < nu; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      long num = 0;
      for (intT k=0; k < v_deg; ++k) {
        if (edgesV[v_offset+k] < i) num ++;
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
intT getNextWedgeIdx(bipartiteCSR& GA, bool use_v, long max_wedges, intT curr_idx, long* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  //nextWedgeTimer.start();
  if (nu - curr_idx < 2000) return getNextWedgeIdx_seq(GA, use_v, max_wedges, curr_idx);

  auto idx_map = make_in_imap<long>(nu - curr_idx, [&] (size_t i) { return wedge_idxs[curr_idx+i+1] - wedge_idxs[curr_idx]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }
  //nextWedgeTimer.stop();
  return find_idx; //TODO make sure right
}

//TODO 3 tuple instead of nested pairs
template<class wedge, class wedgeCons>
  pair<long, intT> getWedges(_seq<wedge>& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges, long* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  if (max_wedges >= num_wedges) {
    _getWedges<wedge>(wedges, GA, use_v, cons, num_wedges, wedge_idxs);
    return make_pair(num_wedges, nu);
  }
  long next_idx = getNextWedgeIdx(GA, use_v, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedges<wedge>(wedges, GA, use_v, cons, num_wedges, wedge_idxs, curr_idx, next_idx);
  return make_pair(num_wedges, next_idx);
}

template<class wedgeCons, class T>
  intT getWedgesHash(T& wedges, bipartiteCSR& GA, bool use_v, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges, long* wedge_idxs) {
  const long nu = use_v ? GA.nu : GA.nv;
  if (max_wedges >= num_wedges) {
    _getWedgesHash(wedges, GA, use_v, cons, num_wedges);
    return nu;
  }
  intT next_idx = getNextWedgeIdx(GA, use_v, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedgesHash(wedges, GA, use_v, cons, num_wedges, curr_idx, next_idx);
  return next_idx;
}

//***************************************************************************************************
//***************************************************************************************************

template<class wedge, class wedgeCons>
  void _getWedges_seq(wedge* wedges, graphCSR& GA, wedgeCons cons, long num_wedges, intT curr_idx, intT next_idx) {
  long idx = 0;
  for(intT i=curr_idx; i < next_idx; ++i) {
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    for(intT j=0; j < u_deg; ++j) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
      if (v <= i) break;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = GA.edges[v_offset+k] >> 1;
        if (u2 > i) {
          wedges[idx] = cons(i, GA.edges[v_offset+k], v, j, k);
          ++idx;
        }
        else break; 
      }
    }
  }
}

template<class wedge, class wedgeCons>
  void _getWedges(_seq<wedge>& wedges_seq, graphCSR& GA, wedgeCons cons, long num_wedges, 
		  long* wedge_idxs, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  // Allocate space for seagull storage
  if (wedges_seq.n < num_wedges) {
    free(wedges_seq.A);
    wedges_seq.A = newA(wedge, num_wedges);
    wedges_seq.n = num_wedges;
  }

  if (next_idx == INT_T_MAX) next_idx = GA.n;
  //if (num_wedges < 10000)
  return _getWedges_seq<wedge>(wedges_seq.A, GA, cons, num_wedges, curr_idx, next_idx);
 
  // Store seagulls in parallel
  parallel_for(intT i=curr_idx; i < next_idx; ++i) {
    long wedge_idx = wedge_idxs[i] - wedge_idxs[curr_idx];
    // Consider each neighbor v of active vertex u
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    long idx = 0;
    for(intT j=0; j < u_deg; ++j) { //JS: test granular_fr
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1] - v_offset;
      if (v <= i) break;
      // Find neighbors (not equal to u) of v
      for (intT k = 0; k < v_deg; ++k) {
        uintE u2 = GA.edges[v_offset+k] >> 1;
        if (u2 > i) {
          wedges_seq.A[wedge_idx+idx] = cons(i, GA.edges[v_offset+k], v, j, k);
          ++idx;
        }
        else break; //{if (i < next_idx-1) assertf(wedge_idx+idx == wedge_idxs[i+1] - wedge_idxs[curr_idx], "%d, %d", wedge_idx+idx,wedge_idxs[i+1] - wedge_idxs[curr_idx] ); break;}
      }
    }
  }
}

// Note: 1 in LSB if we store i/u2, 0 if we store v
template<class wedgeCons, class T>
  void _getWedgesHash(T& wedges, graphCSR& GA, wedgeCons cons, long num_wedges, intT curr_idx=0, intT next_idx=INT_T_MAX) {
  if (next_idx == INT_T_MAX) next_idx = GA.n;

  wedges.resize(num_wedges);
  parallel_for(intT i=curr_idx; i < next_idx; ++i){
    // Set up for each active vertex
    intT u_offset = GA.offsets[i];
    intT u_deg = GA.offsets[i+1] - u_offset;
    granular_for(j,0,u_deg,u_deg>1000,{
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1] - v_offset;
	if (v > i) {
	  // Find all seagulls with center v and endpoint u
	  for (intT k=0; k < v_deg; ++k) { 
	    uintE u2 = GA.edges[v_offset+k] >> 1;
	    if (u2 > i) wedges.insert(make_pair(cons(i, u2, (GA.edges[v_offset+k] & 0b1) ), 1));
	    else break;
	  }
	} 
      });
  }
}

intT getNextWedgeIdx(graphCSR& GA, long max_wedges, intT curr_idx, long* wedge_idxs) {
  //if (GA.n - curr_idx < 2000) return getNextWedgeIdx_seq(GA, use_v, max_wedges, curr_idx);
  auto idx_map = make_in_imap<long>(GA.n - curr_idx, [&] (size_t i) { return wedge_idxs[curr_idx+i+1] - wedge_idxs[curr_idx]; });
  auto lte = [] (const long& l, const long& r) { return l <= r; };
  size_t find_idx = pbbs::binary_search(idx_map, max_wedges, lte) + curr_idx; //this rets first # > searched num
  if (find_idx == curr_idx) {cout << "Space must accomodate seagulls originating from one vertex\n"; exit(0); }

  return find_idx;
}

template<class wedgeCons, class T>
  intT getWedgesHash(T& wedges, graphCSR& GA, wedgeCons cons, long max_wedges, 
		     intT curr_idx, long num_wedges, long* wedge_idxs) {
  if (max_wedges >= num_wedges) {
    _getWedgesHash(wedges, GA, cons, num_wedges);
    return GA.n;
  }
  intT next_idx = getNextWedgeIdx(GA, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedgesHash(wedges, GA, cons, num_wedges, curr_idx, next_idx);
  return next_idx;
}

template<class wedge, class wedgeCons>
  pair<long, intT> getWedges(_seq<wedge>& wedges, graphCSR& GA, wedgeCons cons, long max_wedges, intT curr_idx, long num_wedges, long* wedge_idxs) {
  if (max_wedges >= num_wedges) {
    _getWedges<wedge>(wedges, GA, cons, num_wedges, wedge_idxs);
    return make_pair(num_wedges, GA.n);
  }
  long next_idx = getNextWedgeIdx(GA, max_wedges, curr_idx, wedge_idxs);
  num_wedges = wedge_idxs[next_idx] - wedge_idxs[curr_idx];
  _getWedges<wedge>(wedges, GA, cons, num_wedges, wedge_idxs, curr_idx, next_idx);
  return make_pair(num_wedges, next_idx);
}

//***************************************************************************************************
//***************************************************************************************************


// eti[u_offset + j] = v_offset + i, where i refers to u and j refers to v
uintE* edgeToIdxs(bipartiteCSR& GA, bool use_v) {
  long nu = use_v ? GA.nu : GA.nv;
  long nv = use_v ? GA.nv : GA.nu;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  uintE* eti = newA(uintE, GA.numEdges);

  parallel_for(intT i=0; i < nu; ++i) {
    intT u_offset = offsetsU[i];
    intT u_deg = offsetsU[i+1] - u_offset;
    granular_for(j,0,u_deg,u_deg>1000,{ // (intT j = 0; j < u_deg; ++j) { //JS: test granular_for
      
      uintE v = edgesU[u_offset + j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1] - v_offset;
      // find k such that edgesV[v_offset + k] = i
      auto idx_map = make_in_imap<uintE>(v_deg, [&] (size_t k) { return edgesV[v_offset + k]; });
      auto lte = [] (const uintE& l, const uintE& r) { return l < r; };
      size_t find_idx = pbbs::binary_search(idx_map, i, lte); 
      eti[u_offset+j] = v_offset + find_idx;
      });
  }

  return eti;
}

uintE* idxsToEdge(bipartiteCSR& G, bool use_v) {
  return edgeToIdxs(G, !use_v);
}

/*pair<uintE, uintE> idxsToEdge(uintE idx, uintE* ite) {
  uintE u = edgesV[idx];
  uintE v = edgesU[ite[idx]];
  // we want to go through all neighbors of v
  // intersect each N(u') with N(u)
  // get all u, v, u2
  // intersection of N(u) and N(u2) - 1 is the num butterflies to subtract from v, u2
  }*/

#endif
