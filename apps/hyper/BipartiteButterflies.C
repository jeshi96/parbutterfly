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
#include "sparseSet.h"
#include "sampleSort.h"
#include "../../lib/histogram.h"
#include "../../lib/gbbs-histogram.h"

#include <assert.h>

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }

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

// Represents a pair of vertices on one side of a bipartite graph
struct VertexPair {
	uintE v1;
	uintE v2;
	VertexPair(uintE _v1, uintE _v2) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1) {}
};

struct VertexPairCmp {
	long nv;
	VertexPairCmp(long _nv) : nv(_nv) {}
	bool operator() (VertexPair vs1, VertexPair vs2) {
		return vs1.v1 * nv + vs1.v2 < vs2.v1 * nv + vs2.v2;
	}
};

struct VertexPairCmp2 {
	long nv;
	VertexPairCmp2(long _nv) : nv(_nv) {}
	bool operator() (VertexPair vs1, VertexPair vs2) {
		return vs1.v2 * nv + vs1.v1 < vs2.v2 * nv + vs2.v1;
	}
};

// Constructs a VertexPair
struct VertexPairCons { VertexPair operator() (uintE v1, uintE v2) { return VertexPair(v1, v2); }};

// Constructs a uintE form of a VertexPair
struct VertexPairIntCons {
  long nu;
  VertexPairIntCons(long _nu) : nu(_nu) {}
  uintE operator() (uintE v1, uintE v2) {
    return v1 <= v2 ? v1 * nu + v2 : v2 * nu + v1;
  }
};

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
bipartiteGraph<vertex> bpGraphFromFile(char* iFile){
	hypergraph<vertex> G = readHypergraphFromFile<vertex>(iFile,1,0);
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

// Returns the number of wedges on all vertices on one side of our bipartite graph, as specified by V
template<class vertex>
long countWedges(long nv, vertex* V) {
  return sequence::reduce<long>((long) 0, nv, addF<long>(), chooseV<vertex, long>(V));
}

// Compares the wedge counts on one side of our bipartite graph versus the other; returns the side with the least number of wedges
template <class vertex>
pair<bool,long> cmpWedgeCounts(bipartiteGraph<vertex> GA) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = countWedges<vertex>(nv,GA.V);
  long num_wedges_u = countWedges<vertex>(nu, GA.U);
  return make_pair((num_wedges_v <= num_wedges_u), num_wedges_v <= num_wedges_u ? num_wedges_v : num_wedges_u);
}

template<class wedge, class vertex, class wedgeCons>
wedge* getWedges(const long nv, const vertex* V, wedgeCons cons) {
timer t;
t.start();
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
        wedges[wedge_idx + idx] = cons(v.getOutNeighbor(j), v.getOutNeighbor(k));
        ++idx;
      }
    }
  }
  free(wedge_idxs);
t.stop();
t.reportTotal("\tgetWedges:");

  return wedges;
}

//********************************************************************************************
//********************************************************************************************

// Not parallel
void storeButterfliesSortCE_seq(uintE* butterflies, VertexPair* wedges, uintE* wedge_freqs_f, 
                            uintE* par_idxs_f, long i, bool use_v1) {
timer t;
t.start();
  // Retrieve the frequency counts for each distinct key by taking the difference between
  // the wedge_freqs_f indices
  // Because of how we set up the par_idxs_f, we know that every element in this range has the same
  // first vertex
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  for (long idx = par_idxs_f[i]; idx < par_idxs_f[i+1]; ++idx) { //TODO this is without reduce -- do thresholding
      uintE num_butterflies = wedge_freqs_f[idx+1] - wedge_freqs_f[idx];
      long wedge_idx = wedge_freqs_f[idx];
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      
  	  VertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies[vs.v1] += num_butterflies;
      else butterflies[vs.v2] += num_butterflies;
    }
t.stop();
//t.reportTotal("\tstoreButterfliesSortCE_seq:");
}

// Parallel
void storeButterfliesSortCE(uintE* butterflies, long nu, VertexPair* wedges, uintE* wedge_freqs_f, 
                            uintE* par_idxs_f, long i, bool use_v1) {
timer t;
t.start();
  // Retrieve the frequency counts for each distinct key by taking the difference between
  // the wedge_freqs_f indices
  // Because of how we set up the par_idxs_f, we know that every element in this range has the same
  // first vertex, and a different second vertex
  // Store the number of butterflies on each key by the second vertex (no parallel issues b/c distinct)
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // The sum gives the number of butterflies on the first vertex
  uintE* butterflies_v = newA(uintE, nu);
  parallel_for(long i=0; i < nu; ++i) {butterflies_v[i] = 0; }
  parallel_for (long idx = par_idxs_f[i]; idx < par_idxs_f[i+1]; ++idx) {
      uintE num_butterflies = wedge_freqs_f[idx+1] - wedge_freqs_f[idx];
      long wedge_idx = wedge_freqs_f[idx];
      num_butterflies = num_butterflies * (num_butterflies - 1)/2;
      
  	  VertexPair vs = wedges[wedge_idx];

      if (use_v1) butterflies_v[vs.v2] = num_butterflies;
      else butterflies_v[vs.v1] = num_butterflies;
  }
  long sum = sequence::plusReduce(butterflies_v, nu);
  VertexPair start_vs = wedges[wedge_freqs_f[par_idxs_f[i]]];
  if (use_v1) butterflies[start_vs.v1] += sum;
  else butterflies[start_vs.v2] += sum;
t.stop();
//t.reportTotal("\tstoreButterfliesSortCE:");
}


// Retrieve frequency counts for all wedges with the same key, and for all wedges with the same first vertex
// (if use_v1; otherwise, second vertex -- in comments, we assume first vertex)
// First, retrieve a list of indices where consecutive wedges have different keys
pair<uintE*, long> getWedgeFreqs(long nv, long nu, VertexPair* wedges, long num_wedges, bool use_v1) {
timer t;
t.start();
  // Sort the wedges by the key (pair of two vertices in U)
  if (use_v1) sampleSort(wedges, num_wedges, VertexPairCmp(nu)); //TODO check if num_wedges > 2^32 to see if overflow; if it is, compile w/env var LONG/EDGELONG
  else sampleSort(wedges, num_wedges, VertexPairCmp2(nu));

  uintE* wedge_freqs = newA(uintE, num_wedges+1);
  wedge_freqs[0] = 0;
  wedge_freqs[num_wedges] = num_wedges;

  parallel_for (long i = 1; i < num_wedges; ++i) {
  	if((wedges[i-1].v1 != wedges[i].v1) || (wedges[i-1].v2 != wedges[i].v2))
  		wedge_freqs[i] = i;
  	else
  		wedge_freqs[i] = UINT_E_MAX;
  }
  uintE* wedge_freqs_f = newA(uintE, num_wedges+1);
  long num_wedge_freqs_f = sequence::filter(wedge_freqs, wedge_freqs_f, num_wedges+1, nonMaxF());
  free(wedge_freqs);
  return make_pair(wedge_freqs_f, num_wedge_freqs_f);
t.stop();
t.reportTotal("\tgetWedgeFreqs:");
}

void countButterfliesSortCE(uintE* butterflies, long nv, long nu, VertexPair* wedges, uintE* wedge_freqs_f, 
            long num_wedge_freqs_f, bool use_v1) {
timer t;
t.start();
  // Given this list wedge_freqs_f, retrieve a list of indices for wedge_freqs_f such that
  // consecutive wedges have different first vertices
  // That is to say, the list we construct par_idxs_f is a list of indices on wedge_freqs_f,
  // such that wedges[wedge_freqs_f[par_idxs_f[i]]] and wedges[wedge_freqs_f[par_idxs_f[i+1]]]
  // have different first vertices
  uintE* par_idxs = newA(uintE, num_wedge_freqs_f);
  par_idxs[0] = 0;
  par_idxs[num_wedge_freqs_f-1] = num_wedge_freqs_f-1;
  parallel_for (long i=1; i < num_wedge_freqs_f - 1; ++i) {
    uintE v_prev = use_v1 ? wedges[wedge_freqs_f[i-1]].v1 : wedges[wedge_freqs_f[i-1]].v2;
    uintE v_curr = use_v1 ? wedges[wedge_freqs_f[i]].v1 : wedges[wedge_freqs_f[i]].v2;
    if(v_prev != v_curr)
      par_idxs[i] = i;
    else
      par_idxs[i] = UINT_E_MAX;
  }
  uintE* par_idxs_f = newA(uintE, num_wedge_freqs_f);
  long num_par_idxs_f = sequence::filter(par_idxs, par_idxs_f, num_wedge_freqs_f, nonMaxF());

  // Use these two lists to retrieve the number of butterflies in a cache-efficient manner
  // Start by iterating through our par_idxs_f list
  parallel_for (long i = 0; i < num_par_idxs_f-1; ++i) {
    // The difference between consecutive elements tells us how many distinct values in wedge_freqs_f
    // have the same first vertex; use this to threshold whether to parallelize or not
    long range = par_idxs_f[i+1] - par_idxs_f[i];
    if (range > 10000) storeButterfliesSortCE(butterflies,nu,wedges,wedge_freqs_f,par_idxs_f,i,use_v1);
    else storeButterfliesSortCE_seq(butterflies,wedges,wedge_freqs_f,par_idxs_f,i,use_v1); //TODO set threshold var
  }

t.stop();
t.reportTotal("\tcountButterfliesSortCE:");
}

void countButterfliesSort(uintE* butterflies, VertexPair* wedges, uintE* wedge_freqs_f, long num_wedge_freqs_f) {
timer t;
t.start();
  // Then, retrieve the frequency counts for each distinct key by taking the difference between
  // these indices
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in U
  parallel_for (long i = 1; i < num_wedge_freqs_f; ++i) {
  	uintE num_butterflies = wedge_freqs_f[i] - wedge_freqs_f[i-1];
  	long wedge_idx = wedge_freqs_f[i-1];

  	num_butterflies = num_butterflies * (num_butterflies - 1)/2;

  	VertexPair vs = wedges[wedge_idx];

  	writeAdd(&butterflies[vs.v1],num_butterflies); 
  	writeAdd(&butterflies[vs.v2],num_butterflies);
  }
t.stop();
t.reportTotal("\tcountButterfliesSort:");
}

// This is the original compute function, without the more cache-efficient sorting method
template <class vertex>
void ComputeSort(bipartiteGraph<vertex> GA, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  long num_wedges = use_v.second;

  VertexPair* wedges = getWedges<VertexPair>(nv, V, VertexPairCons());

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
  	butterflies[i] = 0; 
  }

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges, true);

  countButterfliesSort(butterflies, wedges, freq_pair.first, freq_pair.second);

  free(freq_pair.first);

  /*for (long i=0; i < nu; ++i) {
  	cout << i << ", " << butterflies[i] << "\n";
  }*/
  free(butterflies);
  free(wedges);
}

// This is the new compute function, with cache efficient sorting
template <class vertex>
void ComputeSortCE(bipartiteGraph<vertex> GA, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  // TODO if 0? make sure doesn't break
  long num_wedges = use_v.second;

  VertexPair* wedges = getWedges<VertexPair>(nv, V, VertexPairCons());

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
  	butterflies[i] = 0; 
  }

  // Retrieve butterflies + store with first vertex
  pair<uintE*, long> freq_pair = getWedgeFreqs(nv, nu, wedges, num_wedges, true);
  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair.first, freq_pair.second , true);
  free(freq_pair.first);

  // Retrieve butterflies + store with second vertex
  pair<uintE*, long> freq_pair2 = getWedgeFreqs(nv, nu, wedges, num_wedges, false);
  countButterfliesSortCE(butterflies, nv, nu, wedges, freq_pair2.first, freq_pair2.second , false);
  free(freq_pair2.first);

  /*for (long i=0; i < nu; ++i) {
    cout << i << ", " << butterflies[i] << "\n";
  }*/
  free(butterflies);
  free(wedges);

}

//********************************************************************************************
//********************************************************************************************


template<class vertex>
void getWedgesHash(sparseAdditiveSet<uintE>& wedges, long nv, long nu,vertex* V) {
timer t;
t.start();
// Count number of wedges by their key
  parallel_for (long i = 0; i < nv; ++i) {
  	const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
      	VertexPair vs = VertexPair(v.getOutNeighbor(j), v.getOutNeighbor(k));
      	wedges.insert(pair<uintE,uintE>(vs.v1 * nu + vs.v2, 1));
      }
    }
  }
t.stop();
t.reportTotal("\tgetWedgesHash:");
}

template <class vertex>
void ComputeHash(bipartiteGraph<vertex> GA, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  long num_wedges = use_v.second;
  // TODO fix prob don't need that ceil stuff just divide w/float -- for mem on hash table
  float f = ((float)num_wedges)/((float) (nv*nu+nu));
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nv*nu + nu,f,
  	UINT_E_MAX);
  getWedgesHash(wedges,nv, nu,V);
  
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
  	butterflies[i] = 0;
  }
timer t;
t.start();
  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
  	pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
  	uintE num_butterflies = wedge_freq_pair.second;
  	uintE v2 = wedge_freq_pair.first % nu;
  	uintE v1 = wedge_freq_pair.first / nu;

  	num_butterflies = num_butterflies * (num_butterflies - 1)/2;

  	writeAdd(&butterflies[v1],num_butterflies);
  	writeAdd(&butterflies[v2],num_butterflies);
  }
t.stop();
t.reportTotal("\tcountButterfliesHash:");

  /*for (long i=0; i < nu; ++i) {
  	cout << i << ", " << butterflies[i] << "\n";
  }*/

  free(butterflies);
  wedge_freqs.del();
  wedges.del();
}

template <class vertex>
void ComputeHashCE(bipartiteGraph<vertex> GA, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  long num_wedges = use_v.second;
  // TODO fix prob don't need that ceil stuff just divide w/float -- for mem on hash table
  float f = ((float)num_wedges)/((float) (nv*nu+nu));
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nv*nu + nu,f,
  	UINT_E_MAX);
  getWedgesHash(wedges,nv, nu,V);
  
  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  sparseAdditiveSet<uintE> butterflies_set = sparseAdditiveSet<uintE>(nu,1,UINT_E_MAX);

timer t;
t.start();
  // Retrieve count on each key; that number choose 2 is the number of butterflies
  parallel_for (long i=0; i < wedge_freqs.n; ++i) {
  	pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
  	uintE num_butterflies = wedge_freq_pair.second;
  	uintE v2 = wedge_freq_pair.first % nu;
  	uintE v1 = wedge_freq_pair.first / nu;

  	num_butterflies = num_butterflies * (num_butterflies - 1)/2;
  	butterflies_set.insert(pair<uintE,uintE>(v1, num_butterflies));
  	butterflies_set.insert(pair<uintE,uintE>(v2, num_butterflies));
  }

  _seq<pair<uintE,uintE>> butterflies_seq = butterflies_set.entries();
  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
  parallel_for(long i=0; i < butterflies_seq.n; ++i) {
    pair<uintE, uintE> butterfly_pair = butterflies_seq.A[i];
    butterflies[butterfly_pair.first] = butterfly_pair.second;
  }
t.stop();
t.reportTotal("\tcountButterfliesHashCE:");
  
  /*for (long i=0; i < nu; ++i) {
  	cout << i << ", " << butterflies[i] << "\n";
  }*/

  free(butterflies);
  //butterflies.del();
  butterflies_set.del();
  wedge_freqs.del();
  wedges.del();
}

//********************************************************************************************
//********************************************************************************************


// TODO use more efficient hist from laxman
template <class vertex>
void ComputeHist(bipartiteGraph<vertex> GA, bool gbbs, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  long num_wedges = use_v.second;

  uintE* wedges_list = getWedges<uintE>(nv, V, VertexPairIntCons(nu));
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.
  tuple<uintE, uintE>* wedge_freqs;
  size_t wedge_freqs_n;

timer t2;
t2.start();
  if (gbbs) {
  	tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
      pbbsa::gbbs::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu);
  	wedge_freqs = get<1>(wedges_tuple);
 	  wedge_freqs_n = get<0>(wedges_tuple);
  }
  else {
  	tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
      pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu);
  	wedge_freqs = get<1>(wedges_tuple);
  	wedge_freqs_n = get<0>(wedges_tuple);
  }
t2.stop();
t2.reportTotal("\thistWedges:");

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }

timer t;
t.start();
  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    tuple<uintE,uintE> wedge_freq_pair = wedge_freqs[i];
    uintE num_butterflies = get<1>(wedge_freq_pair);
    uintE wedge_freq_pair_first = get<0>(wedge_freq_pair);

    uintE v2 = wedge_freq_pair_first % nu;
    uintE v1 = wedge_freq_pair_first / nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    writeAdd(&butterflies[v1],num_butterflies);
    writeAdd(&butterflies[v2],num_butterflies);
  }
t.stop();
t.reportTotal("\tcountButterfliesHist:");

  /*for (long i=0; i < nu; ++i) {
    cout << i << ", " << butterflies[i] << "\n";
  }*/

  free(butterflies);
  free(wedge_freqs);
  free(wedges_list);
}

template<class E, class K>
E chooseAdd (E curr, tuple<K,E> v) {
		return curr + get<1>(v);
}

template<class E, class K>
tuple<K,E> chooseAddReduce (tuple<K,E> curr, tuple<K,E> v) {
		return make_tuple(get<0>(curr),get<1>(curr) + get<1>(v));
}

// TODO use more efficient hist from laxman
template <class vertex>
void ComputeHistCE(bipartiteGraph<vertex> GA, commandLine P) {
  pair<bool,long> use_v = cmpWedgeCounts(GA);
  const long nv = use_v.first ? GA.nv : GA.nu;
  const long nu = use_v.first ? GA.nu : GA.nv;
  const vertex* V = use_v.first ? GA.V : GA.U;
  const vertex* U = use_v.first ? GA.U : GA.V;

  long num_wedges = use_v.second;

  uintE* wedges_list = getWedges<uintE>(nv, V, VertexPairIntCons(nu));
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges);
  //JS: one optimization is to fuse getWedgesInt into the histogram
  //code, instead of creating a new sequence, so that we don't need to
  //actually write out all the wedges.
timer t2;
t2.start();
  tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
    pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu);
  tuple<uintE, uintE>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
using X = tuple<uintE,uintE>;
  X* wedge_freqs_i = newA(X, 2*wedge_freqs_n);
  parallel_for(long i = 0; i < wedge_freqs_n; i++) {
    tuple<uintE,uintE> wedge_freq_pair = wedge_freqs[i];
    uintE num = get<1>(wedge_freq_pair);
    num = (num * (num-1)) / 2;
    uintE wedge_num = get<0>(wedge_freq_pair);
    wedge_freqs_i[2*i] = make_tuple(wedge_num % nu, num);
    wedge_freqs_i[2*i + 1] = make_tuple(wedge_num / nu, num);
  }
  pbbsa::sequence<tuple<uintE, uintE>> wedge_freqs_i_seq = pbbsa::sequence<tuple<uintE,uintE>>(wedge_freqs_i,2*wedge_freqs_n);
  tuple<size_t,tuple<uintE, uintE>*> butterflies_tuple = 
    pbbsa::sparse_histogram_f<uintE,uintE>(wedge_freqs_i_seq,nu,chooseAdd<uintE,uintE>, chooseAddReduce<uintE,uintE>);
  tuple<uintE,uintE>* butterflies_l = get<1>(butterflies_tuple);
  size_t butterflies_n = get<0>(butterflies_tuple);
t2.stop();
t2.reportTotal("\thistWedgesCE:");

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
  parallel_for (long i=0; i < butterflies_n; ++i) {
    butterflies[get<0>(butterflies_l[i])] = get<1>(butterflies_l[i]);
  }

  /*for (long i=0; i < nu; ++i) {
    cout << i << ", " << butterflies[i] << "\n";
  }*/

  free(butterflies);
  free(wedge_freqs);
  free(wedges_list);
}

int parallel_main(int argc, char* argv[]){
	commandLine P(argc,argv," [-i <ingraph>] [-t <type>] [-nv <numleftvertices>] [-nu <numrightvertices>] <inFile>");
  	std::string iFileConst = P.getOptionValue("-i", "");
    char* iFile = new char[iFileConst.length()+1];
    strcpy(iFile, iFileConst.c_str());
    long ty = P.getOptionLongValue("-t",0);
  	long nv = P.getOptionLongValue("-nv", 100);
  	long nu = P.getOptionLongValue("-nu", 100);
	bipartiteGraph<symmetricVertex> G = (iFileConst.length() != 0) ? 
		bpGraphFromFile<symmetricVertex>(iFile) : bpGraphComplete<symmetricVertex>(nv,nu);
  delete [] iFile;

  timer t;
  t.start();
  if (ty == 0) ComputeSort(G,P);
  else if (ty==1) ComputeSortCE(G,P);
  else if (ty==2) ComputeHash(G,P);
  else if (ty==3) ComputeHashCE(G,P);
  else if(ty==4) ComputeHist(G,false,P);
  else if(ty==5) ComputeHist(G,true,P);
  else ComputeHistCE(G,P);
	t.stop();

  if (ty==0) t.reportTotal("Sort:");
  else if (ty==1) t.reportTotal("SortCE:");
  else if (ty==2) t.reportTotal("Hash:");
  else if (ty==3) t.reportTotal("HashCE:");
  else if (ty==4) t.reportTotal("Hist:");
  else if (ty==5) t.reportTotal("HistNT:");
  else t.reportTotal("HistCE:");

	G.del();
}
