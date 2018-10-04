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
#include "../../lib/histogram.h"

#include <assert.h>

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }

using namespace std;

struct VertexIndSet {
	uintE v1;
	uintE v2;
	VertexIndSet(uintE _v1, uintE _v2) : v1(_v1<=_v2 ? _v1 : _v2), v2(_v1<=_v2 ? _v2 : _v1) {}
};

struct VertexIndSetCmp {
	long nv;
	VertexIndSetCmp(long _nv) : nv(_nv) {}
	bool operator() (VertexIndSet vs1, VertexIndSet vs2) {
		return vs1.v1 * nv + vs1.v2 < vs2.v1 * nv + vs2.v2;
	}
};

struct VertexIndSetEq {
	bool operator() (VertexIndSet vs1, VertexIndSet vs2) {
		return (vs1.v1 == vs2.v1) && (vs1.v2 == vs2.v2);
	}
};

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

template<class vertex>
long countWedges(long nv, vertex* V) {
  long num_wedges = 0;
  parallel_for (long i = 0; i < nv; ++i) { //TODO fix
    uintE v_deg = V[i].getOutDegree();
    if (v_deg >= 2)
      writeAdd(&num_wedges, (long) (v_deg * (v_deg - 1)/2));
  }
  return num_wedges;
}

template <class vertex>
bool cmpWedgeCounts(bipartiteGraph<vertex> GA) {
  const long nv = GA.nv, nu = GA.nu;
  long num_wedges_v = countWedges<vertex>(nv,GA.V);
  long num_wedges_u = countWedges<vertex>(nu, GA.U);
  return (num_wedges_v <= num_wedges_u);
}

template<class vertex>
VertexIndSet* getWedges(long nv, const vertex* V) {
  // Retrieve the indices of each wedge associated with each vertex in GA.V
  //long wedge_idxs[nv+1];
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

  // Retrieve each wedge associated with each vertex in GA.V
  VertexIndSet* wedges = newA(VertexIndSet,num_wedges);
  parallel_for (long i = 0; i < nv; ++i) {
    const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    long wedge_idx = wedge_idxs[i];
    long idx = 0;
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
        assert(wedge_idx + idx < num_wedges);
        //assertf(v.getOutNeighbor(j), "%li \n %x \n", j,v.getOutDegree())
        wedges[wedge_idx + idx] = VertexIndSet(v.getOutNeighbor(j), v.getOutNeighbor(k));
        ++idx;
      }
    }
  }
  free(wedge_idxs);

  return wedges;
}

template <class vertex>
void Compute(bipartiteGraph<vertex> GA, commandLine P) {
  bool use_v = cmpWedgeCounts(GA);
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  long num_wedges = countWedges(nv,V);

  VertexIndSet* wedges = getWedges<vertex>(nv, V);

  // Sort the wedges by the key (pair of two vertices in GA.H)
  quickSort(wedges,num_wedges,VertexIndSetCmp(nu));

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  uintE* butterflies = newA(uintE, nu);
  //uintE butterflies[nu];
  parallel_for(long i=0;i<nu;++i){
  	butterflies[i] = 0;
  }
  uintE* wedge_freqs = newA(uintE, num_wedges-1);
  //uintE wedge_freqs[num_wedges-1];
  VertexIndSetEq vertexIndSetEq = VertexIndSetEq();
  parallel_for (long i = 0; i < num_wedges - 1; ++i) {
  	if (!vertexIndSetEq(wedges[i],wedges[i+1]))
  		wedge_freqs[i] = i;
  	else
  		wedge_freqs[i] = UINT_E_MAX;
  }
  uintE* wedge_freqs_f = newA(uintE, num_wedges-1);
  //uintE wedge_freqs_f[num_wedges-1];
  long num_wedge_freqs_f = sequence::filter(wedge_freqs, wedge_freqs_f, num_wedges-1, nonMaxF());

  // Then, retrieve the frequency counts for each distinct key by taking the difference between
  // these indices
  // Take the frequency count choose 2 to receive the number of butterflies on that key
  // Store butterfly counts with each vertex in GA.H
  parallel_for (long i = 0; i < num_wedge_freqs_f + 1; ++i) { //TODO make add to butterflies ok
  	uintE num_butterflies = 0;
  	long wedge_idx;
  	if (i==0) {
  		num_butterflies = wedge_freqs_f[i] + 1;
  		wedge_idx = wedge_freqs_f[i];
  	}
  	else if (i == num_wedge_freqs_f){
  		num_butterflies = num_wedges - 1 - wedge_freqs_f[i-1];
  		wedge_idx = num_wedges - 1;
  	}
  	else{
  		num_butterflies = wedge_freqs_f[i] - wedge_freqs_f[i-1];
  		wedge_idx = wedge_freqs_f[i];
  	}

  	num_butterflies = num_butterflies * (num_butterflies - 1)/2;
  	assert(wedge_idx < num_wedges);
  	VertexIndSet vs = wedges[wedge_idx];

  	assert (vs.v1 < nu);
  	assert (vs.v2 < nu);
  	writeAdd(&butterflies[vs.v1],num_butterflies);
  	writeAdd(&butterflies[vs.v2],num_butterflies);
  }

  /*for (long i=0; i < nu; ++i) {
  	cout << i << ", " << butterflies[i] << "\n";
  }*/
  free(butterflies);
  free(wedge_freqs_f);
  free(wedge_freqs);
  free(wedges);
}

template <class vertex>
void ComputeHash(bipartiteGraph<vertex> GA, commandLine P) {
  bool use_v = cmpWedgeCounts(GA);
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  long num_wedges = countWedges(nv,V);
  float f = ceilf(((float)num_wedges)/((float) (nv*nu+nu))*1000000)/1000000;
  sparseAdditiveSet<uintE> wedges = sparseAdditiveSet<uintE>(nv*nu + nu,f,
  	UINT_E_MAX);

  parallel_for (long i = 0; i < nv; ++i) {
  	const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
      	VertexIndSet vs = VertexIndSet(v.getOutNeighbor(j), v.getOutNeighbor(k));
      	wedges.insert(pair<uintE,uintE>(vs.v1 * nu + vs.v2, 1));
      }
    }
  }

  _seq<pair<uintE,uintE>> wedge_freqs = wedges.entries();
  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
  	butterflies[i] = 0;
  }
  parallel_for (long i=0; i < wedge_freqs.n; ++i) { //TODO ISSUE parallel protect
  	pair<uintE,uintE> wedge_freq_pair = wedge_freqs.A[i];
  	uintE num_butterflies = wedge_freq_pair.second;
  	uintE v2 = wedge_freq_pair.first % nu;
  	uintE v1 = (wedge_freq_pair.first - v2) / nu;

  	num_butterflies = num_butterflies * (num_butterflies - 1)/2;

  	writeAdd(&butterflies[v1],num_butterflies);
  	writeAdd(&butterflies[v2],num_butterflies);
  }

  /*for (long i=0; i < nu; ++i) {
  	cout << i << ", " << butterflies[i] << "\n";
  }*/

  free(butterflies);
  wedge_freqs.del();
  wedges.del();
}

template<class vertex>
uintE* getWedgesInt(long nv, long nu, const vertex* V) {
  // Retrieve the indices of each wedge associated with each vertex in GA.V
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

  // Retrieve each wedge associated with each vertex in GA.V
  uintE* wedges = newA(uintE,num_wedges);
  parallel_for (long i = 0; i < nv; ++i) {
    const vertex v = V[i];
    const uintE v_deg = v.getOutDegree();
    long wedge_idx = wedge_idxs[i];
    long idx = 0;
    for (long j = 0; j < v_deg; ++j) {
      for (long k = j+1; k < v_deg; ++k) {
        assert(wedge_idx + idx < num_wedges);
        //assertf(v.getOutNeighbor(j), "%li \n %x \n", j,v.getOutDegree())
        VertexIndSet vs = VertexIndSet(v.getOutNeighbor(j), v.getOutNeighbor(k));
        wedges[wedge_idx + idx] = vs.v1 * nu + vs.v2;
        ++idx;
      }
    }
  }

  free(wedge_idxs);

  return wedges;
}

template <class vertex>
void ComputeHist(bipartiteGraph<vertex> GA, commandLine P) {
  bool use_v = cmpWedgeCounts(GA);
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const vertex* V = use_v ? GA.V : GA.U;
  const vertex* U = use_v ? GA.U : GA.V;

  long num_wedges = countWedges(nv,V);

  uintE* wedges_list = getWedgesInt<vertex>(nv, nu, V);
  pbbsa::sequence<uintE> wedges_seq = pbbsa::sequence<uintE>(wedges_list,num_wedges);

  tuple<size_t,tuple<uintE, uintE>*> wedges_tuple = 
    pbbsa::sparse_histogram<uintE, uintE>(wedges_seq, nv*nu + nu);
  tuple<uintE, uintE>* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);

  uintE* butterflies = newA(uintE, nu);
  parallel_for(long i=0;i<nu;++i){
    butterflies[i] = 0;
  }
  parallel_for (long i=0; i < wedge_freqs_n; ++i) { //TODO ISSUE parallel protect
    tuple<uintE,uintE> wedge_freq_pair = wedge_freqs[i];
    uintE num_butterflies = get<1>(wedge_freq_pair);
    uintE wedge_freq_pair_first = get<0>(wedge_freq_pair);

    uintE v2 = wedge_freq_pair_first % nu;
    uintE v1 = (wedge_freq_pair_first - v2) / nu;

    num_butterflies = num_butterflies * (num_butterflies - 1)/2;

    writeAdd(&butterflies[v1],num_butterflies);
    writeAdd(&butterflies[v2],num_butterflies);
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
  if (ty == 0) Compute(G,P);
  else if (ty==1) ComputeHash(G,P);
  else ComputeHist(G,P);
	t.stop();
  if (ty==0) t.reportTotal("time sort: ");
  else if (ty==1) t.reportTotal("time hash: ");
  else t.reportTotal("time hist: ");
	G.del();
	/*parallel_for(int i=0; i < nv; ++i) {
		v[i].del();
	}
	parallel_for(int i=0; i < nu; ++i) {
		u[i].del();
	}
	free(v);
	free(u);*/
}
