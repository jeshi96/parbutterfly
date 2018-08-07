#define HYPER 1
#include "hygra.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"

template <class vertex>
struct Init_Deg {
  long* Degrees;
  vertex* H;
  Init_Deg(long* _Degrees, vertex* _H) : Degrees(_Degrees), H(_H) {}
  inline bool update (uintE s, uintE d) { 
    xadd(&Degrees[s],(long)H[d].getOutDegree()-1);
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    Degrees[s]+=H[d].getOutDegree()-1;
    return 1;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

//return true for neighbor the first time it's updated 
struct Count_Removed {
  intE* Counts;
  Count_Removed(intE* _Counts) : Counts(_Counts) {}
  inline bool update (uintE s, uintE d) { 
    return ++Counts[d] == 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    return xadd(&Counts[d],1) == 1;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

struct Update_Deg {
  long* Degrees;
  intE *Counts;
  Update_Deg(long* _Degrees, intE* _Counts) : Degrees(_Degrees), Counts(_Counts) {}
  inline bool update (uintE s, uintE d) {
    Degrees[d] -= Counts[s];
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d){
    xadd(&Degrees[d],(long)-Counts[s]);
    return 1;
  }
  inline bool cond (uintE d) { return Degrees[d] > 0; }
};

template <class vertex>
array_imap<long> KCore(hypergraph<vertex>& GA, size_t num_buckets=128) {
  const size_t nv = GA.nv, mv = GA.mv, nh = GA.nh;
  bool* active = newA(bool,nv);
  {parallel_for(long i=0;i<nv;i++) active[i] = 1;}
  vertexSubset Frontier(nv, nv, active);
  long* Degrees = newA(long,nv);
  {parallel_for(long i=0;i<nv;i++) {
      Degrees[i] = 0;
    }}
  intE* Counts = newA(intE,nh);
  {parallel_for(long i=0;i<nh;i++) Counts[i] = 0;}
  edgeMap(GA,FROM_V,Frontier,Init_Deg<vertex>(Degrees,GA.H),INT_T_MAX,no_output);
  auto D = array_imap<long>(Degrees,nv);
  //auto D = array_imap<uintE>(nv, [&] (size_t i) { return GA.V[i].getOutDegree(); });

  auto em = EdgeMapHypergraph<uintE, vertex>(GA, make_tuple(UINT_E_MAX, 0), (size_t)GA.mv/5);
  auto b = make_buckets(nv, D, increasing, num_buckets);

  size_t finished = 0;
  while (finished != nv) {
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    long k = bkt.id;
    finished += active.size();

    auto apply_f = [&] (const tuple<uintE, uintE>& p) -> const Maybe<tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      long deg = D.s[v];
      if (deg > k) {
        long new_deg = max(deg - edgesRemoved, k);
        D.s[v] = new_deg;
        uintE bkt = b.get_bucket(deg, new_deg);
        return wrap(v, bkt);
      }
      return Maybe<tuple<uintE, uintE> >();
    };

    vertexSubset FrontierH = edgeMap(GA,FROM_V,active,Count_Removed(Counts));
    cout << "num active = " << active.numNonzeros() << " frontierH = " << FrontierH.numNonzeros() << endl;
    //edgeMap(GA,FROM_H,FrontierH,Update_Deg(Degrees,Counts),-1,no_output);

    auto map_f = [&] (const uintE& i, const uintE& j) { return Counts[i]; };
    auto reduce_f = [&] (const uintE& cur, const tuple<uintE, intE>& r) { return cur + std::get<1>(r); };
    vertexSubsetData<uintE> moved = em.template edgeMapReduce<uintE, intE>(FrontierH, FROM_V, map_f, reduce_f, apply_f);

    auto reset_counts = [&] (const uintE& i) { Counts[i] = 0; };
    vertexMap(FrontierH,reset_counts);

    //vertexSubsetData<uintE> moved = em.template edgeMapCount<uintE>(active, FROM_V, apply_f);
    b.update_buckets(moved.get_fn_repr(), moved.size());
    moved.del(); active.del();
  }
  return D;
}

template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 128);
  if (num_buckets != (1 << pbbs::log2_up(num_buckets))) {
    cout << "Number of buckets must be a power of two." << endl;
    exit(-1);
  }
  cout << "### Application: k-core" << endl;
  cout << "### Graph: " << P.getArgument(0) << endl;
  cout << "### Workers: " << getWorkers() << endl;
  cout << "### Buckets: " << num_buckets << endl;

  auto cores = KCore(GA, num_buckets);
  long mc = 0;
  for (size_t i=0; i < GA.nv; i++) { mc = std::max(mc, cores[i]); }
  cout << "### Max core: " << mc << endl;
}
