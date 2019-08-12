#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "vertex.h"
#include "sequence.h"
#include "sparseSet.h"
#include "sampleSort.h"
#include "../../lib/histogram.h"
#include "../../lib/sample_sort.h"
#include "../../radixsort/RadixSort/radixSort.h"

#include "butterfly_utils.h"

using namespace std;

pair<intT,long> CountSortTotal(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  pair<long, intT> wedges_pair  = getWedges<UVertexPair>(cs.wedges_seq_uvp, GA, use_v, UVertexPairCons(), max_wedges, curr_idx, num_wedges,
							 wedge_idxs);

  UVertexPair* wedges = cs.wedges_seq_uvp.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  pbbs::sample_sort(wedges, num_wedges_curr, UVertexPairCmp());
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());

  long* butterflies = newA(long, freq_pair.second-1);

  const size_t eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num = freq_pair.first[i+1] - freq_pair.first[i];
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }

  free(freq_pair.first);
  long num_butterflies = sequence::plusReduce(butterflies, freq_pair.second-1);
  free(butterflies);
  return make_pair(wedges_pair.second,num_butterflies);
}

pair<intT,long> CountSortTotal(CountSpace& cs, graphCSR& GA, long num_wedges, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  pair<long, intT> wedges_pair  = getWedges<UWedge>(cs.wedges_seq_uw, GA, UWedgeCons(), max_wedges, curr_idx, num_wedges,
						    wedge_idxs);
  UWedge* wedges = cs.wedges_seq_uw.A;
  long num_wedges_curr = wedges_pair.first;

  // Retrieve frequency counts for all wedges with the same key
  // First, retrieve a list of indices where consecutive wedges have different keys
  pbbs::sample_sort(wedges, num_wedges_curr, UWedgeCmp());
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  long* butterflies = newA(long, freq_pair.second-1);

  const size_t eltsPerCacheLine = 64/sizeof(long);
  parallel_for (long i = 0; i < freq_pair.second-1; ++i) {
    long num = freq_pair.first[i+1] - freq_pair.first[i];
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }

  free(freq_pair.first);
  long num_butterflies = sequence::plusReduce(butterflies, freq_pair.second-1);
  free(butterflies);
  return make_pair(wedges_pair.second,num_butterflies);
}

//************************************************************************************************************************
//************************************************************************************************************************

pair<intT,long> CountHashTotal(CountSpace& cs, graphCSR& GA, long num_wedges, long max_wedges, long* wedge_idxs, 
	       intT curr_idx=0) {
  const size_t eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, UVertexPairIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);

  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);
  long* butterflies = newA(long, num_wedges_seq);

  parallel_for(long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num = wedge_freq_pair.second;
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }
  long num_butterflies = sequence::plusReduce(butterflies, num_wedges_seq);
  free(butterflies);

  return make_pair(next_idx, num_butterflies);
}

pair<intT,long> CountHashTotal(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  intT next_idx = getWedgesHash(cs.wedges_hash, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);

  size_t num_wedges_seq = cs.wedges_hash.entries_no_init(cs.wedges_seq_intp);
  long* butterflies = newA(long, num_wedges_seq);

  // Retrieve count on each key; that number choose 2 is the number of butterflies  
  parallel_for (long i=0; i < num_wedges_seq; ++i) {
    auto wedge_freq_pair = cs.wedges_seq_intp.A[i];
    long num = wedge_freq_pair.second;
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }
  long num_butterflies = sequence::plusReduce(butterflies, num_wedges_seq);
  free(butterflies);

  return make_pair(next_idx, num_butterflies);
}

//************************************************************************************************************************
//************************************************************************************************************************

pair<intT,long> CountHistTotal(CountSpace& cs, graphCSR& GA, long num_wedges, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  using X = tuple<long,uintE>;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  pair<long, intT> wedges_pair  = getWedges<long>(cs.wedges_seq_int, GA, UWedgeIntRankCons(GA.n), max_wedges, curr_idx, num_wedges, wedge_idxs);
  intT next_idx = wedges_pair.second;

  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(cs.wedges_seq_int.A, wedges_pair.first);

  tuple<size_t,X*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, ((GA.n*(GA.n+1)) << 1), cs.tmp, cs.out);
  X* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
  long* butterflies = newA(long, wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    long num = get<1>(wedge_freqs[i]);
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }

  long num_butterflies = sequence::plusReduce(butterflies, wedge_freqs_n);
  free(butterflies);

  return make_pair(next_idx, num_butterflies);
}

pair<intT,long> CountHistTotal(CountSpace& cs, bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long* wedge_idxs, intT curr_idx=0) {
  const long nu = use_v ? GA.nu : GA.nv;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  using T = tuple<long,uintE>;

  pair<long, intT> wedges_list_pair = getWedges<long>(cs.wedges_seq_int, GA, use_v, UVertexPairIntCons(nu), max_wedges, curr_idx, num_wedges, wedge_idxs);
  long* wedges_list = cs.wedges_seq_int.A;
  long num_wedges_curr = wedges_list_pair.first;
  pbbsa::sequence<long> wedges_seq = pbbsa::sequence<long>(wedges_list,num_wedges_curr);

  tuple<size_t,T*> wedges_tuple = pbbsa::sparse_histogram<long, uintE>(wedges_seq, nu*nu, cs.tmp, cs.out);
  T* wedge_freqs = get<1>(wedges_tuple);
  size_t wedge_freqs_n = get<0>(wedges_tuple);
  long* butterflies = newA(long, wedge_freqs_n);

  parallel_for (long i=0; i < wedge_freqs_n; ++i) {
    long num = get<1>(wedge_freqs[i]);
    num = num * (num - 1)/2;
    butterflies[i] = num;
  }

  long num_butterflies = sequence::plusReduce(butterflies, wedge_freqs_n);
  free(butterflies);

  return make_pair(wedges_list_pair.second, num_butterflies);
}

//************************************************************************************************************************
//************************************************************************************************************************


long CountOrigCompactParallel_WedgeAwareTotal(bipartiteCSR& GA, bool use_v, long max_array_size, long* wedgesPrefixSum) {
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;
  const size_t eltsPerCacheLine = 64/sizeof(long);

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); //15 tunable parameter
//#ifdef VERBOSE
timer t1,t2,t3,t4;
  t1.start();
//#endif

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);
  long* results = newA(long, eltsPerCacheLine*stepSize);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[eltsPerCacheLine*i] = 0; });
  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");

  t2.start();
#endif

  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
	  intT used_idx = 0;
	  long shift = nu*(i-step*stepSize);
	  intT u_offset  = offsetsU[i];
	  intT u_deg = offsetsU[i+1]-u_offset;
	  for (intT j=0; j < u_deg; ++j ) {
	    uintE v = edgesU[u_offset+j];
	    intT v_offset = offsetsV[v];
	    intT v_deg = offsetsV[v+1]-offsetsV[v];
	    for (intT k=0; k < v_deg; ++k) {
	      uintE u2_idx = edgesV[v_offset+k];
	      if (u2_idx < i) {
		results[eltsPerCacheLine*(i % stepSize)] += wedges[shift+u2_idx];
		wedges[shift+u2_idx]++;
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	      }
	      else break;
	    }
	  }
	  for(intT j=0; j < used_idx; ++j) { wedges[used[shift+j]+shift] = 0; }
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,nu));
  }
  auto butterflies_extract_f = [&] (const long i) -> const long {
      return results[eltsPerCacheLine*i];
  };
  long num_butterflies = sequence::reduce<long>((long)0,(long)stepSize,addF<long>(),butterflies_extract_f);
#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
  free(results);

  return num_butterflies;
}

long CountOrigCompactParallel_WedgeAwareTotal(graphCSR& GA, long max_array_size, long* wedgesPrefixSum) {
  const size_t eltsPerCacheLine = 64/sizeof(long);

  long stepSize = min<long>(getWorkers() * 40, max_array_size/GA.n); //15 tunable parameter
//#ifdef VERBOSE
  timer t1,t2,t3,t4;
  t1.start();
//#endif

  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);
  long* results = newA(long, eltsPerCacheLine*stepSize);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[eltsPerCacheLine*i] = 0; });
  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");

  t2.start();
#endif
  
  for(intT step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    std::function<void(intT,intT)> recursive_lambda =
      [&]
      (intT start, intT end){
      if ((start == end-1) || (wedgesPrefixSum[end]-wedgesPrefixSum[start] < 1000)){ 
	for (intT i = start; i < end; i++){
      intT used_idx = 0;
      long shift = GA.n*(i-step*stepSize);
      intT u_offset  = GA.offsets[i];
      intT u_deg = GA.offsets[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
	if (v <= i) break; 
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
	  if (u2_idx > i) {
        results[eltsPerCacheLine*(i % stepSize)] += wedges[shift+u2_idx];
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	  }
	  else break;
	}
      }
      
      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        wedges[shift+u2_idx] = 0;
      }
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,GA.n));
  }
  auto butterflies_extract_f = [&] (const long i) -> const long {
      return results[eltsPerCacheLine*i];
  };
  long num_butterflies = sequence::reduce<long>((long)0,(long)stepSize,addF<long>(),butterflies_extract_f);
#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
  free(results);

  return num_butterflies;
}

//************************************************************************************************************************
//************************************************************************************************************************

long CountOrigCompactParallelTotal(bipartiteCSR& GA, bool use_v, long max_array_size) {
//#ifdef VERBOSE
  timer t1,t2,t3;
  t1.start();
//#endif
  const size_t eltsPerCacheLine = 64/sizeof(long);
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); //15 tunable parameter

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);
  long* results = newA(long, eltsPerCacheLine*stepSize);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[eltsPerCacheLine*i] = 0; });
  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");

  t2.start();
#endif

  for(intT step = 0; step < (nu+stepSize-1)/stepSize; step++) {
    parallel_for_1(intT i=step*stepSize; i < min((step+1)*stepSize,nu); ++i){
      intT used_idx = 0;
      long shift = nu*(i-step*stepSize);
      intT u_offset  = offsetsU[i];
      intT u_deg = offsetsU[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = edgesU[u_offset+j];
	intT v_offset = offsetsV[v];
	intT v_deg = offsetsV[v+1]-offsetsV[v];
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = edgesV[v_offset+k];
	  if (u2_idx < i) {
	    results[eltsPerCacheLine*(i % stepSize)] += wedges[shift+u2_idx];
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = u2_idx;
	  }
	  else break;
	}
      }
      for(intT j=0; j < used_idx; ++j) { wedges[used[shift+j]+shift] = 0; }
    }
  }
  auto butterflies_extract_f = [&] (const long i) -> const long {
      return results[eltsPerCacheLine*i];
  };
  long num_butterflies = sequence::reduce<long>((long)0,(long)stepSize,addF<long>(),butterflies_extract_f);
#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif

  free(wedges);
  free(used);
  free(results);

  return num_butterflies;
}

long CountWorkEfficientParallelTotal(graphCSR& GA, long max_array_size) {
//#ifdef VERBOSE
  timer t1,t2,t3;
  t1.start();
//#endif
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long stepSize = min<long>(getWorkers() * 60, max_array_size/GA.n); //15 tunable parameter
  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);
  long* results = newA(long, eltsPerCacheLine*stepSize);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[eltsPerCacheLine*i] = 0; });
  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");

  t2.start();
#endif
  for(long step = 0; step < (GA.n+stepSize-1)/stepSize; step++) {
    parallel_for_1(long i=step*stepSize; i < min((step+1)*stepSize,GA.n); ++i){
      intT used_idx = 0;
      long shift = GA.n*(i-step*stepSize);
      intT u_offset  = GA.offsets[i];
      intT u_deg = GA.offsets[i+1]-u_offset;
      for (intT j=0; j < u_deg; ++j ) {
	uintE v = GA.edges[u_offset+j] >> 1;
	intT v_offset = GA.offsets[v];
	intT v_deg = GA.offsets[v+1]-v_offset;
	if (v <= i) break; 
	for (intT k=0; k < v_deg; ++k) { 
	  uintE u2_idx = GA.edges[v_offset+k] >> 1;
	  if (u2_idx > i) {
        results[eltsPerCacheLine*(i % stepSize)] += wedges[shift+u2_idx];
	    wedges[shift+u2_idx]++;
	    if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = GA.edges[v_offset+k];
	  }
	  else break;
	}
      }
      
      for(long j=0; j < used_idx; ++j) {
        uintE u2_idx = used[shift+j] >> 1;
        wedges[shift+u2_idx] = 0;
      }
    }
  }
  auto butterflies_extract_f = [&] (const long i) -> const long {
      return results[eltsPerCacheLine*i];
  };
  long num_butterflies = sequence::reduce<long>((long)0,(long)stepSize,addF<long>(),butterflies_extract_f);
#ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
#endif
  
  free(wedges);
  free(used);
  free(results);

  return num_butterflies;
}


long CountWorkEfficientSerialTotal(graphCSR& GA) {
  //cout << "Work-efficient Serial (make sure running with CILK_NWORKERS=1)" << endl;
  timer t1,t2,t3;
  t1.start();
  long nb = 0;
  uintE* wedges = newA(uintE, GA.n);
  uintE* used = newA(uintE, GA.n);
  //long* butterflies = newA(long,GA.n);
  const size_t eltsPerCacheLine = 64/sizeof(long);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  for(long i=0;i<GA.n;i++) {
    wedges[i] = 0;
    //butterflies[i] = 0;
  }
  t3.stop();
#ifdef VERBOSE
  t3.reportTotal("preprocess (initialize)");
  t2.start();
#endif

  for(long i=0; i < GA.n; ++i){
    intT used_idx = 0;
    intT u_offset  = GA.offsets[i];
    intT u_deg = GA.offsets[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1]-v_offset;
      if (v <= i) break;
      for (intT k=0; k < v_deg; ++k) { 
	uintE u2_idx = GA.edges[v_offset+k] >> 1;
	if (u2_idx > i) { //TODO combine into one graph
	    nb += wedges[u2_idx];
	  //results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
	  wedges[u2_idx]++;
	  if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
	}
	else break;
      }
    }

    for(long j=0; j < used_idx; ++j) {
      wedges[used[j]] = 0;
    }
  }

  #ifdef VERBOSE
  t2.reportTotal("counting (main loop)");
  #endif
  
  free(wedges);
  free(used);

  return nb;
}

//************************************************************************************************************************
//************************************************************************************************************************

/*
 *  Computes total butterfly count.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges    : Number of wedges produced by the given ranking
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  ranks         : Sorted U and V vertices by rank (where U indices all come before V indices)
 *  rankV         : Array mapping V indices to their rank
 *  rankU         : Array mapping U indices to their rank
 * 
 *  Returns total butterfly count.
 *  Note: This function takes ownership of and frees ranks, rankV, and rankU
 */
long CountRankTotal(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, CountType type,
		                uintE* ranks, uintE* rankV, uintE* rankU) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  #ifdef VERBOSE
  timer t_rank;
  t_rank.start();
  #endif

  // Rank graph using ranks, rankV, and rankU
  // Vertices in g are indexed by rank
  auto g = rankGraph(GA, use_v, ranks, rankV, rankU);
  free(ranks); free(rankU); free(rankV);

  const size_t eltsPerCacheLine = 64/sizeof(long);

  #ifdef VERBOSE
  t_rank.reportTotal("ranking");
  timer t_time;
  t_time.start();
  #endif

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = (type == BATCHS || type == SERIAL) ? nullptr : countWedgesScan(g);
  long num_butterflies = 0;

  #ifdef VERBOSE
  if (type == BATCHWA) t_time.reportTotal("counting (scan)");
  #endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) num_butterflies = CountWorkEfficientParallelTotal(g, max_array_size);
  else if (type == SERIAL) num_butterflies = CountWorkEfficientSerialTotal(g);
  else if (type == BATCHWA) num_butterflies = CountOrigCompactParallel_WedgeAwareTotal(g, max_array_size, wedge_idxs);
  else {
    // Initialize structure that holds all space for counting algorithms, to be
    // reused with wedge batches
    CountSpace cs = CountSpace(type, g.n, false);
  
    if (max_wedges >= num_wedges) {
      if (type == SORT) num_butterflies = CountSortTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == HASH) num_butterflies = CountHashTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == HIST) num_butterflies = CountHistTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
    }
    else {
      intT curr_idx = 0;
      pair<intT,long> ret;
      while(curr_idx < g.n) {
	if (type == SORT) ret = CountSortTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == HASH) ret = CountHashTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == HIST) ret = CountHistTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
  curr_idx = ret.first;
  num_butterflies += ret.second;
	cs.clear();
      }
    }
    cs.del();
  }
  g.del();

  if (type != BATCHS && type != SERIAL) free(wedge_idxs);

  #ifdef VERBOSE
  if (type != BATCHS && type != BATCHWA && type != SERIAL) t_time.reportTotal("counting");
  #endif
  return num_butterflies;
}

/*
 *  Computes total butterfly count.
 * 
 *  GA            : Bipartite graph in CSR format
 *  use_v         : Denotes which bipartition produces the fewest wedges. If true, considering two-hop neighbors of U
 *                  produces the fewest wedges. If false, considering two-hop neighbors of V produces the fewest wedges.
 *  num_wedges    : Number of wedges produced by the bipartition given by use_v, assuming ranking by side
 *  max_wedges    : Max number of wedges (consisting of two endpoint indices and a center) that the system can hold in
 *                  memory
 *  max_array_size:
 *  type          : Wedge/butterfly aggregation type
 *  tw            : Vertex ordering option
 * 
 *  Returns total butterfly count.
 */
long CountTotal(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, CountType type,
                RankType tw) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;

  // Any ranking that is not by side
  if (tw != SIDE) {
    #ifdef VERBOSE
    timer t_rank;
    t_rank.start();
    #endif

    // Pick the correct ranking
    tuple<uintE*,uintE*,uintE*> rank_tup;
    if (tw == COCORE) rank_tup = getCoCoreRanks(GA);
    else if (tw == ACOCORE) rank_tup = getApproxCoCoreRanks(GA);
    else if (tw == DEG) rank_tup = getDegRanks(GA);
    else if (tw == ADEG) rank_tup = getApproxDegRanks(GA);

    // Compute the number of wedges that the ranking produces
    long num_ccwedges =
      sequence::reduce<long>((long) 0, GA.nu, addF<long>(), rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup),
                                                                             get<1>(rank_tup)));
    num_ccwedges +=
      sequence::reduce<long>((long) 0, GA.nv, addF<long>(), rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup),
                                                                             get<2>(rank_tup)));

    #ifdef VERBOSE
    t_rank.reportTotal("ranking");
    #endif

    // Compute total butterfly count using the specified ranking
    return CountRankTotal(GA, use_v, num_ccwedges, max_wedges, max_array_size, type, get<0>(rank_tup), get<1>(rank_tup),
                          get<2>(rank_tup));
  }

  // Ranking by side
  const size_t eltsPerCacheLine = 64/sizeof(long);
  #ifdef VERBOSE
  timer t_time;
  t_time.start();
  #endif

  // Compute wedge indices so that wedges can be stored in parallel (for any
  // aggregation type except simple batching)
  long* wedge_idxs = (type == BATCHS) ? nullptr : countWedgesScan(GA, use_v, true);

  // Initialize structure that holds all space for counting algorithms, to be
  // reused with wedge batches
  CountSpace cs = CountSpace(type, nu, false);
  long num_butterflies = 0;

  #ifdef VERBOSE
  if (type == BATCHWA) t_time.reportTotal("counting (scan)");
  #endif

  // Choose wedge/butterfly aggregation type
  if (type == BATCHS) num_butterflies = CountOrigCompactParallelTotal(GA, use_v, max_array_size);
  else if (type == BATCHWA) num_butterflies = CountOrigCompactParallel_WedgeAwareTotal(GA,use_v,max_array_size,wedge_idxs);
  else {
    if (max_wedges >= num_wedges) {
      if (type == SORT) num_butterflies = CountSortTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == HASH) num_butterflies = CountHashTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == HIST) num_butterflies = CountHistTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
    }
    else {
      intT curr_idx = 0;
      pair<intT,long> ret;
      while(curr_idx < nu) {
	if (type == SORT) ret = CountSortTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == HASH) ret = CountHashTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == HIST) ret = CountHistTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
  curr_idx = ret.first;
  num_butterflies += ret.second;
	cs.clear();
      }
    }
  }

  if (type != BATCHS) free(wedge_idxs);
  cs.del();

  #ifdef VERBOSE
  if (type != BATCHS && type != BATCHWA) t_time.reportTotal("counting");
  #endif

  return num_butterflies;
}