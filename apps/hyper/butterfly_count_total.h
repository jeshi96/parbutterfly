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
  radix::parallelIntegerSort<uintE>(wedges, num_wedges_curr, UVPFirst());
  radix::parallelIntegerSort<uintE>(wedges, num_wedges_curr, UVPSecond());
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UVertexPairCmp(), UVertexPairEq(), LONG_MAX, nonMaxLongF());

  long* butterflies = newA(long, freq_pair.second-1);

  const intT eltsPerCacheLine = 64/sizeof(long);
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
  radix::parallelIntegerSort<uintE>(wedges, num_wedges_curr, UWFirst());
  radix::parallelIntegerSort<uintE>(wedges, num_wedges_curr, UWSecond());
  auto freq_pair = getFreqs<long>(wedges, num_wedges_curr, UWedgeCmp(), UWedgeEq(), LONG_MAX, nonMaxLongF());
  long* butterflies = newA(long, freq_pair.second-1);

  const intT eltsPerCacheLine = 64/sizeof(long);
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
  const intT eltsPerCacheLine = 64/sizeof(long);

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
  const intT eltsPerCacheLine = 64/sizeof(long);

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
  const intT eltsPerCacheLine = 64/sizeof(long);

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
  const intT eltsPerCacheLine = 64/sizeof(long);

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
  const intT eltsPerCacheLine = 64/sizeof(long);

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); //15 tunable parameter
//#ifdef VERBOSE
timer t1,t2,t3,t4;
  t1.start();
//#endif

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);
  long* results = newA(long, stepSize*eltsPerCacheLine);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[i*eltsPerCacheLine] = 0; });
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
		results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
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
      return results[i*eltsPerCacheLine];
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
  const intT eltsPerCacheLine = 64/sizeof(long);

  long stepSize = min<long>(getWorkers() * 40, max_array_size/GA.n); //15 tunable parameter
//#ifdef VERBOSE
  timer t1,t2,t3,t4;
  t1.start();
//#endif

  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);
  long* results = newA(long, stepSize*eltsPerCacheLine);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[i*eltsPerCacheLine] = 0; });
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
        results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
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
      return results[i*eltsPerCacheLine];
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
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = min<long>(getWorkers() * 60, max_array_size/nu); //15 tunable parameter

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);
  long* results = newA(long, stepSize*eltsPerCacheLine);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[i*eltsPerCacheLine] = 0; });
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
	    results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
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
      return results[i*eltsPerCacheLine];
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
  const intT eltsPerCacheLine = 64/sizeof(long);
  long stepSize = min<long>(getWorkers() * 60, max_array_size/GA.n); //15 tunable parameter
  uintE* wedges = newA(uintE, GA.n*stepSize);
  uintE* used = newA(uintE, GA.n*stepSize);
  long* results = newA(long, stepSize*eltsPerCacheLine);
  t1.stop();
#ifdef VERBOSE
  t1.reportTotal("preprocess (malloc)");
#endif
  t3.start();
  granular_for(i,0,GA.n*stepSize,GA.n*stepSize > 10000, { wedges[i] = 0; });
  granular_for(i,0,stepSize,stepSize > 10000, { results[i*eltsPerCacheLine] = 0; });
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
        results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
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
      return results[i*eltsPerCacheLine];
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


long CountRankTotal(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, long type,
		uintE* ranks, uintE* rankV, uintE* rankU) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
#ifdef VERBOSE
  timer t_rank;
  t_rank.start();
#endif
  auto g = rankGraph(GA, use_v, ranks, rankV, rankU);
  free(ranks); free(rankU); free(rankV);
#ifdef VERBOSE
  t_rank.reportTotal("ranking");
#endif

  const intT eltsPerCacheLine = 64/sizeof(long);
#ifdef VERBOSE
  timer t_time;
  t_time.start();
#endif

  long* wedge_idxs = (type == 11) ? nullptr : countWedgesScan(g);

  long num_butterflies = 0;
#ifdef VERBOSE
  if (type == 8) t_time.reportTotal("counting (scan)");
#endif
  if (type == 11) num_butterflies = CountWorkEfficientParallelTotal(g, max_array_size);
  else if (type == 8) num_butterflies = CountOrigCompactParallel_WedgeAwareTotal(g, max_array_size, wedge_idxs);
  else {
    CountSpace cs = CountSpace(type, g.n, false);
    if (max_wedges >= num_wedges) {
      if (type == 0) num_butterflies = CountSortTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == 2) num_butterflies = CountHashTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == 4) num_butterflies = CountHistTotal(cs, g, num_wedges, max_wedges, wedge_idxs).second;
    }
    else {
      intT curr_idx = 0;
      pair<intT,long> ret;
      while(curr_idx < g.n) {
	if (type ==0) ret = CountSortTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type ==2) ret = CountHashTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == 4) ret = CountHistTotal(cs, g, num_wedges, max_wedges, wedge_idxs, curr_idx);
    curr_idx = ret.first; num_butterflies += ret.second;
	cs.clear();
      }
    }
    cs.del();
  }
  g.del();

  if (type != 11) free(wedge_idxs);
#ifdef VERBOSE
  if (type != 11 && type != 8) t_time.reportTotal("counting");
#endif
  return num_butterflies;
}

long CountTotal(bipartiteCSR& GA, bool use_v, long num_wedges, long max_wedges, long max_array_size, long type=0, long tw=0) {
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  // tw 0 is use side, tw 1 is use co core ranks, tw 2 is use approx co core ranks, tw 3 is use the better b/w approx co core and side (TODO put in deg for tw 3)
  if (tw != 0) {
#ifdef VERBOSE
    timer t_rank;
    t_rank.start();
#endif
    //auto rank_tup = getDegRanks(GA);
    // auto rank_tup = getCoreRanks(GA);
    tuple<uintE*,uintE*,uintE*> rank_tup;
    if (tw == 1) rank_tup = getCoCoreRanks(GA);
    else if (tw == 2) rank_tup = getApproxCoCoreRanks(GA);
    else if (tw == 3) rank_tup = getDegRanks(GA);
    else if (tw == 4) rank_tup = getApproxDegRanks(GA);

    //long num_rwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    //num_rwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));

    //long num_cwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup2), get<1>(rank_tup2)));
    //num_cwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
    //  rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup2), get<2>(rank_tup2)));

    long num_ccwedges = sequence::reduce<long>((long) 0, GA.nu, addF<long>(),
					       rankWedgeF<long>(GA.offsetsU, GA.edgesU, get<2>(rank_tup), get<1>(rank_tup)));
    num_ccwedges += sequence::reduce<long>((long) 0, GA.nv, addF<long>(),
					   rankWedgeF<long>(GA.offsetsV, GA.edgesV, get<1>(rank_tup), get<2>(rank_tup)));
#ifdef VERBOSE
    t_rank.reportTotal("ranking");
#endif
    //cout << "Rank wedges: " << num_rwedges << "\n";
    //cout << "Core wedges: " << num_cwedges << "\n";
    /*
    cout << "Side wedges: " << num_wedges << "\n";
    cout << "Co Core wedges: " << num_ccwedges << "\n"; */

    if (num_ccwedges < num_wedges + 1000 || tw == 1 || tw == 2 || tw == 3) return CountRankTotal(GA, use_v, num_ccwedges, max_wedges, max_array_size, type, get<0>(rank_tup), get<1>(rank_tup), get<2>(rank_tup));
    free(get<0>(rank_tup)); free(get<1>(rank_tup)); free(get<2>(rank_tup));
  }
  const intT eltsPerCacheLine = 64/sizeof(long);
#ifdef VERBOSE
  timer t_time;
  t_time.start();
#endif
  long* wedge_idxs = (type == 7 || type == 11) ? nullptr : countWedgesScan(GA, use_v, true);
  CountSpace cs = CountSpace(type, nu, false);
  long num_butterflies = 0;
#ifdef VERBOSE
  if (type == 8) t_time.reportTotal("counting (scan)");
#endif

  if (type == 7) num_butterflies = CountOrigCompactParallelTotal(GA, use_v, max_array_size);
  else if (type == 8) num_butterflies = CountOrigCompactParallel_WedgeAwareTotal(GA,use_v,max_array_size,wedge_idxs);
  else {

    if (max_wedges >= num_wedges) {
      if (type == 0) num_butterflies = CountSortTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == 2) num_butterflies = CountHashTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
      else if (type == 4) num_butterflies = CountHistTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs).second;
    }
    else {
      intT curr_idx = 0;
      pair<intT,long> ret;
      while(curr_idx < nu) {
	if (type == 0) ret = CountSortTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == 2) ret = CountHashTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
	else if (type == 4) ret = CountHistTotal(cs, GA, use_v, num_wedges, max_wedges, wedge_idxs, curr_idx);
    curr_idx = ret.first; num_butterflies += ret.second;
	cs.clear();
      }
    }
  }

  if (type != 7 && type != 11) free(wedge_idxs);
  cs.del();
#ifdef VERBOSE
  if (type != 7 && type != 11 && type != 8) t_time.reportTotal("counting");
#endif
  return num_butterflies;

}