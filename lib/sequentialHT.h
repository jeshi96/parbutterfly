#pragma once

#include "utilities.h"
using namespace std;

namespace pbbsa {

template <class K, class V>
class sequentialHT {
 typedef tuple<K,V> T;

 public:
  size_t m;
  size_t mask;
  T empty;
  K max_key;
  T* table;
  bool alloc;
  size_t n_elms;

  inline size_t toRange(size_t h) {return h & mask;}
  inline size_t firstIndex(K v) {return toRange(pbbs::hash64(v));}
  inline size_t incrementIndex(size_t h) {return toRange(h+1);}

  sequentialHT(T* _table, size_t size, float loadFactor, tuple<K, V> _empty) :
    m((size_t) 1 << pbbs::log2_up((size_t)(loadFactor*size))),
    mask(m-1),
    empty(_empty),
    table(_table), alloc(false), n_elms(0) {
      max_key = get<0>(empty);
      if (table == nullptr) {
        table = pbbs::new_array_no_init<T>(m);
        alloc = true;
        for (size_t i=0; i<m; i++) {
          table[i] = empty;
        }
      }
    }

  // m must be a power of two
  sequentialHT(T* _table, size_t _m, tuple<K, V> _empty) :
    m((size_t)_m), mask(m-1), table(_table), empty(_empty), alloc(false), n_elms(0) { max_key = get<0>(empty); }

  inline void del() {
    if (alloc) {
      alloc = false;
      free(table);
    }
  }

  template <class M, class F>
  inline void insertF(tuple<K, M>& v, F& f) {
    K vKey = get<0>(v);
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = get<0>(table[h]);
      if (k == max_key) {
        get<0>(table[h]) = vKey;
        V cur = get<1>(table[h]);
        get<1>(table[h]) = f(cur, v);
        n_elms++;
        return;
      } else if (k == vKey) {
        V cur = get<1>(table[h]);
        get<1>(table[h]) = f(cur, v);
        return;
      }
      h = incrementIndex(h);
    }
  }

  // V must support ++
  inline void insertAdd(K& vKey) {
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = get<0>(table[h]);
      if (k == max_key) {
        table[h] = make_tuple(vKey, 1);
        n_elms++;
        return;
      } else if (k == vKey) {
        get<1>(table[h])++;
        return;
      }
      h = incrementIndex(h);
    }
  }

  // V must support ++, T<1> must be numeric
  inline void insertAdd(T& v) {
    const K& vKey = get<0>(v);
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = get<0>(table[h]);
      if (k == max_key) {
        table[h] = make_tuple(vKey, 1);
        n_elms++;
        return;
      } else if (k == vKey) {
        get<1>(table[h]) += get<1>(v);
        return;
      }
      h = incrementIndex(h);
    }
  }

  inline T find(K& v) {
    size_t h = firstIndex(v);
    T c = table[h];
    while (1) {
      if (get<0>(c) == max_key) {
        return empty;
      } else if (get<0>(c) == v) {
      	return c;
      }
      h = incrementIndex(h);
      c = table[h];
    }
  }

  inline size_t compactInto(T* out) {
    size_t k = 0;
    for (size_t i=0; i<m; i++) {
      auto kv = table[i]; auto key = get<0>(kv);
      if (key != max_key) {
        table[i] = empty;
        out[k++] = kv;
      }
    }
    return k;
  }

  inline tuple<size_t, T*> compact() {
    T* out = pbbs::new_array_no_init<T>(m);
    size_t k = 0;
    for (size_t i=0; i<m; i++) {
      auto kv = table[i]; auto key = get<0>(kv);
      if (key != max_key) {
        table[i] = empty;
        out[k++] = kv;
      }
    }
    return make_tuple(k, out);
  }

};

} // namespace pbbs

