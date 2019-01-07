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
struct VertexPairCons { VertexPair operator() (uintE v1, uintE v2) { return VertexPair(v1, v2); }};
struct UVertexPairCons { UVertexPair operator() (uintE v1, uintE v2) { return UVertexPair(v1, v2); }};

// Constructs a uintE form of a VertexPair and UVertexPair
struct VertexPairIntCons {
  long nu;
  VertexPairIntCons(long _nu) : nu(_nu) {}
  uintE operator() (uintE v1, uintE v2) {
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

struct nonZeroF{bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) != 0);}};
struct greaterOneF{bool operator() (tuple<uintE,uintE> &a) {return (get<1>(a) > 1);}};

struct uintELt {bool operator () (uintE a, uintE b) {return a < b;};};
struct uintEEq {bool operator() (uintE a, uintE b) {return a == b;};};
struct uintETupleLt {bool operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) {return get<0>(a) < get<0>(b);} };
struct uintETupleEq {bool operator() (tuple<uintE,uintE> a,tuple<uintE,uintE> b) {return get<0>(a) == get<0>(b); }};
struct uintETupleAdd {
  tuple<uintE,uintE> operator() (tuple<uintE,uintE> a, tuple<uintE,uintE> b) const {
    return make_tuple(get<0>(a), get<1>(a) + get<1>(b));
  };
};

//***********************************************************************************************
// Graph construction + reading

// Returns a complete bipartite graph with nv vertices on one side and nu vertices on the other
template<class vertex>
bipartiteGraph<vertex> bpGraphComplete(long nv, long nu){
  	vertex* v = newA(vertex,nv+1);
  	vertex* u = newA(vertex,nu+1);
  	parallel_for(int i=0;i<nv;++i) {
  		uintE* neighbors_v = newA(uintE,nu);
  		parallel_for(int j=0;j<nu;++j) {
  			neighbors_v[j] = j;
  		}
  		v[i] = vertex(neighbors_v,nu);
  	}
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

// Reads a bipartite graph from file iFile (builds off of readHypergraphFromFile) -- not terribly efficient
template<class vertex>
bipartiteGraph<vertex> bpGraphFromFile(hypergraph<vertex> G){
	long nv = G.nv;
	long nu = G.nh;
	vertex* v = newA(vertex,nv);
	vertex* u = newA(vertex,nu);
	parallel_for(int i=0;i<nv;++i) {
		uintE* neighbors_v = newA(uintE,G.V[i].getOutDegree());
		parallel_for(int j=0;j<G.V[i].getOutDegree();++j) {
			neighbors_v[j] = G.V[i].getOutNeighbor(j);
		}
		v[i] = vertex(neighbors_v, G.V[i].getOutDegree());
	}
	parallel_for(int i=0;i<nu;++i) {
		uintE* neighbors_u = newA(uintE,G.H[i].getInDegree());
		parallel_for(int j=0;j<G.H[i].getInDegree();++j) {
			neighbors_u[j] = G.H[i].getInNeighbor(j);
		}
		u[i] = vertex(neighbors_u, G.H[i].getInDegree());
	}
	G.del();
	Uncompressed_Membipartitegraph<vertex>* mem = new Uncompressed_Membipartitegraph<vertex>(v,u,nv,nu);
	return bipartiteGraph<vertex>(v,u,nv,nu,mem);

}

//***********************************************************
// general utils

// Sort objects using cmp
// Then, retrieve indices of where consecutive objects are not equal (as given by eq)
// If we let arr be the returned array, arr[i+1] - arr[i] is the frequency of objs[arr[i]]
// Iterating from 0 (inclusive) to the returned length-1 (exclusive) will give all frequency counts
template <class T, class Cmp, class Eq>
pair<uintE*, long> getFreqs(T* objs, long num, Cmp cmp, Eq eq) {
  // Sort the wedges by the key
  sampleSort(objs, num, cmp);

  uintE* freqs = newA(uintE, num + 1);
  freqs[0] = 0;
  freqs[num] = num;
  parallel_for(long i=1; i < num; ++i) {
    if (!eq(objs[i-1],objs[i])) freqs[i] = i;
    else freqs[i] = UINT_E_MAX;
  }
  uintE* freqs_f = newA(uintE, num+1);
  long num_freqs_f = sequence::filter(freqs, freqs_f, num+1, nonMaxF());
  free(freqs);
  return make_pair(freqs_f, num_freqs_f);
}