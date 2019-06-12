#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

#define HYPER 1
//#define LONG 1
//#define EDGELONG 1
#define MCX16 1

#include "hypergraphIO.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"

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

#include "butterfly_count.h"
#include "butterfly_ecount.h"
#include "butterfly_peel.h"
#include "butterfly_epeel.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }

#include <vector>

using namespace std;

void CountOrigCompactSerial(bipartiteCSR& GA, bool use_v) {
  timer t1,t2;
  t1.start();
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  cout << "Original Serial (make sure running with CILK_NWORKERS=1)" << endl;  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long results = 0;
  uintE* wedges = newA(uintE, nu);
  uintE* used = newA(uintE, nu);

  for(long i=0; i < nu; ++i) { wedges[i] = 0; }

  long* butterflies = newA(long,nu);
  for(long i=0; i < nu; ++i) { butterflies[i] = 0; }

  t1.reportTotal("preprocess");
  t2.start();

  for(intT i=0; i < nu; ++i){
    intT used_idx = 0;
    intT u_offset  = offsetsU[i];
    intT u_deg = offsetsU[i+1]-u_offset;
    for (intT j=0; j < u_deg; ++j ) {
      uintE v = edgesU[u_offset+j];
      intT v_offset = offsetsV[v];
      intT v_deg = offsetsV[v+1]-offsetsV[v];
      for (intT k=0; k < v_deg; ++k) { 
        uintE u2_idx = edgesV[v_offset+k];
        if (u2_idx < i) {
          butterflies[i] += wedges[u2_idx];
          butterflies[u2_idx] += wedges[u2_idx];
          //results += wedges[u2_idx];
          wedges[u2_idx]++;
          if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
        }
        else break;
      }
    }
    for(intT j=0; j < used_idx; ++j) { wedges[used[j]] = 0; }
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);

  for(long i=0;i<nu;i++) results += butterflies[i];
  free(butterflies);
  cout << "num: " << results/2 << "\n";
}

void CountOrigCompactParallel_WedgeAware(bipartiteCSR& GA, bool use_v, long* wedgesPrefixSum) {
  timer t1,t2;
  t1.start();
  cout << "Original Parallel Wedge-Aware" << endl;
  //cout << GA.nv << " " << GA.nu << " " << GA.numEdges << endl;
  
  const long nv = use_v ? GA.nv : GA.nu;
  const long nu = use_v ? GA.nu : GA.nv;
  uintT* offsetsV = use_v ? GA.offsetsV : GA.offsetsU;
  uintT* offsetsU = use_v ? GA.offsetsU : GA.offsetsV;
  uintE* edgesV = use_v ? GA.edgesV : GA.edgesU;
  uintE* edgesU = use_v ? GA.edgesU : GA.edgesV;

  long stepSize = 1000; //tunable parameter

  uintE* wedges = newA(uintE, nu*stepSize);
  uintE* used = newA(uintE, nu*stepSize);

  
  granular_for(i,0,nu*stepSize,nu*stepSize > 10000, { wedges[i] = 0; });
  const intT eltsPerCacheLine = 64/sizeof(long);
  
  long* results = newA(long,eltsPerCacheLine*stepSize); //one entry per cache line
  granular_for(i,0,stepSize,stepSize > 10000, {results[eltsPerCacheLine*i] = 0;});
  
  uintE* butterflies = newA(uintE,nu);
  granular_for(i,0,nu,nu > 10000, { butterflies[i] = 0; });
  t1.reportTotal("preprocess");

  t2.start();
  
  //JS: try wedge-aware parallelism using wedge counts per vertex
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
		//butterflies[i] += wedges[u2_idx];
		//butterflies[u2_idx] += wedges[u2_idx];
		results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
		wedges[shift+u2_idx]++;
		if (wedges[shift+u2_idx] == 1) used[shift+used_idx++] = shift+u2_idx;
	      }
	      else break;
	    }
	  }
	  for(intT j=0; j < used_idx; ++j) { wedges[used[shift+j]] = 0; }
	}
      } else {
	cilk_spawn recursive_lambda(start, start + ((end-start) >> 1));
	recursive_lambda(start + ((end-start)>>1), end);
      }
    }; 
    recursive_lambda(step*stepSize,min((step+1)*stepSize,nu));
  }
  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);
  long total = 0;
  
  for(long i=0;i<nu;i++) total += butterflies[i];
  for(long i=0;i<stepSize;i++) total += results[i*eltsPerCacheLine];
  free(butterflies);
  free(results);
  cout << "num: " << total << "\n";
}

void CountWorkEfficientSerial(graphCSR& GA) {
  cout << "Work-efficient Serial (make sure running with CILK_NWORKERS=1)" << endl;
  timer t1,t2;
  t1.start();
  uintE* wedges = newA(uintE, GA.n);
  uintE* used = newA(uintE, GA.n);
  long* butterflies = newA(long,GA.n);
  for(long i=0;i<GA.n;i++) {
    wedges[i] = 0;
    butterflies[i] = 0;
  }
  t1.reportTotal("preprocess");
  t2.start();
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
	  if (GA.edges[v_offset+k] & 0b1) {
	    butterflies[i] += wedges[u2_idx];
	    butterflies[u2_idx] += wedges[u2_idx];
	  }
	  //results[(i % stepSize)*eltsPerCacheLine] += wedges[shift+u2_idx];
	  wedges[u2_idx]++;
	  if (wedges[u2_idx] == 1) used[used_idx++] = u2_idx;
	}
	else break;
      }
    }

    for (long j=0; j < u_deg; ++j ) {
      uintE v = GA.edges[u_offset+j] >> 1;
      intT v_offset = GA.offsets[v];
      intT v_deg = GA.offsets[v+1]-v_offset;
      if (v <= i) break;
      if (!(GA.edges[u_offset+j] & 0b1)) continue;
      for (long k=0; k < v_deg; ++k) { 
	uintE u2_idx = GA.edges[v_offset+k] >> 1;
	if (u2_idx > i) { //TODO combine into one graph
	  if (wedges[u2_idx] > 1) butterflies[v] += wedges[u2_idx]-1;
	}
	else break;
      }
    }

    for(long j=0; j < used_idx; ++j) {
      wedges[used[j]] = 0;
    }
  }

  t2.reportTotal("main loop");
  
  free(wedges);
  free(used);
  long total = 0;
  for(long i=0;i<GA.n;i++) total += butterflies[i];
  free(butterflies);
  cout << "num: " << total/2 << "\n";
}

// Note: must be invoked with symmetricVertex
void Compute(bipartiteCSR& GA, commandLine P) {
  // Method type for counting + peeling
  long ty = P.getOptionLongValue("-t",0);
  long tp = P.getOptionLongValue("-tp",0);
  long te = P.getOptionLongValue("-e",0);
  long tw = P.getOptionLongValue("-w",0);
  bool nopeel = P.getOptionValue("-nopeel");
  
  // # of max wedges
  long max_wedges = P.getOptionLongValue("-m",2577500000);
  long max_array_size = P.getOptionLongValue("-a",23090996160);

  cout << "count " << ty << ", " << "peel " << tp << ", " << "edge " << te << ", " << "rank " << tw << ", " << "peel " << nopeel << "\n";

  timer t1;
  t1.start();
  tuple<bool,long> use_v_tuple = cmpWedgeCounts(GA);
  bool use_v = get<0>(use_v_tuple);
  long num_wedges = get<1>(use_v_tuple);
  t1.reportTotal("preprocess (wedge counts)");
  
  //TODO seq code integrate w/count
  if (te == 0) {
  
    if (ty == 8) {
      long* workPrefixSum = computeWorkPrefixSum(GA,use_v);
      CountOrigCompactParallel_WedgeAware(GA,use_v,workPrefixSum);
      free(workPrefixSum);
    }
    else if (ty == 9) CountOrigCompactSerial(GA,use_v);
    else if (ty == 12) {
      timer t_rank;
      t_rank.start();
      auto rank_tup = getDegRanks(GA);
      auto g = rankGraph(GA, use_v, get<0>(rank_tup), get<1>(rank_tup), get<2>(rank_tup));
      free(get<0>(rank_tup)); free(get<1>(rank_tup)); free(get<2>(rank_tup));
      t_rank.reportTotal("ranking");
      CountWorkEfficientSerial(g);
    }

    if (ty == 8 || ty == 9 || ty == 12) return;
    const intT eltsPerCacheLine = 64/sizeof(long);

    timer t;
    t.start();
    long* butterflies = Count(GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t.stop();

    if (ty==0) t.reportTotal("Sort:");
    else if (ty==1) t.reportTotal("SortCE:");
    else if (ty==2) t.reportTotal("Hash:");
    else if (ty==3) t.reportTotal("HashCE:");
    else if (ty==4) t.reportTotal("Hist:");
    else if (ty==6) t.reportTotal("HistCE:");
    else t.reportTotal("Par");

  
    long num_idxs = use_v ? GA.nu : GA.nv;
    long b = 0;
    for (long i=0; i < num_idxs; ++i) {b += butterflies[eltsPerCacheLine*i];}
    b = b / 2;
    cout << "number of butterflies: " << b << "\n";
  
    //uintE* butterflies2 = Count(GA, use_v, num_wedges, max_wedges, 0, 0);
    //for (long i=0; i < num_idxs; ++i) { assertf(butterflies[eltsPerCacheLine*i] == butterflies2[eltsPerCacheLine*i], "%d, %d, %d", i, butterflies[eltsPerCacheLine*i], butterflies2[eltsPerCacheLine*i]); }
    if(!nopeel) {
      timer t2;
      t2.start();
      auto cores = Peel(GA, use_v, butterflies, max_wedges, tp);
      t2.stop();
      if (tp ==0) t2.reportTotal("Hash Peel:");
      else if (tp==1) t2.reportTotal("Sort Peel:");
      else if (tp==2) t2.reportTotal("Hist Peel:");
      else t2.reportTotal("Par Peel:");

      /*long mc = 0;
      for (size_t i=0; i < num_idxs; i++) { mc = std::max(mc, cores[i]); }
      cout << "### Max core: " << mc << endl;*/
    }
    free(butterflies);
  }
  else {
  
    timer t3;
 
    auto eti = edgeToIdxs(GA, use_v);

    t3.start();
    long* ebutterflies = CountE(eti, GA, use_v, num_wedges, max_wedges, max_array_size, ty, tw);
    t3.stop();
    if(ty==2) t3.reportTotal("E Hash:");
    else if (ty == 3) t3.reportTotal("E HashCE:");
    else if (ty == 0) t3.reportTotal("E Sort:");
    else if (ty==1) t3.reportTotal("E SortCE:");
    else if (ty==4) t3.reportTotal("E Hist:");
    else if (ty==6) t3.reportTotal("E HistCE:");
    else t3.reportTotal("E Par:");

    const intT eltsPerCacheLine = 64/sizeof(long);
    long b=0;
 
    for (long i=0; i < GA.numEdges; ++i) {b += ebutterflies[eltsPerCacheLine*i];}
    cout << "number of edge butterflies: " << b/4 << "\n";

    //uintE* butterflies2 = CountE(eti, GA, use_v, num_wedges, max_wedges, 0, 0);
    //for (long i=0; i < GA.numEdges; ++i) { assertf(ebutterflies[eltsPerCacheLine*i] == butterflies2[eltsPerCacheLine*i], "%d, %d, %d", i, ebutterflies[eltsPerCacheLine*i], butterflies2[eltsPerCacheLine*i]); }

    if(!nopeel) {
      timer t2;
      t2.start();
      auto ite = idxsToEdge(GA, use_v);
      auto cores = PeelE(eti, ite, GA, use_v, ebutterflies, max_wedges, tp);
      free(ite);	    
      t2.stop();
      if (tp ==0) t2.reportTotal("Hash E Peel:");
      else if (tp==1) t2.reportTotal("Sort E Peel:");
      else if (tp == 2) t2.reportTotal("Hist E Peel:");
      else t2.reportTotal("Par E Peel:");
    }
    free(eti);

    free(ebutterflies);
  }
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," <inFile>");
  char* iFile = P.getArgument(0);
  long rounds = P.getOptionLongValue("-r",3);

  bipartiteCSR G = readBipartite(iFile);

  Compute(G,P);
  for(int r=0;r<rounds;r++) {
    Compute(G,P);
  }
  G.del();
}

