#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "vertex.h"
#include "compressedVertex.h"
#include "parallel.h"
using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

// Class that handles implementation specific freeing of memory
// owned by the graph
struct Deletable {
public:
  virtual void del() = 0;
};

template <class vertex>
struct Uncompressed_Mem : public Deletable {
public:
  vertex* V;
  long n;
  long m;
  void* allocatedInplace, * inEdges;

  Uncompressed_Mem(vertex* VV, long nn, long mm, void* ai, void* _inEdges = NULL)
  : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges) { }

  void del() {
    if (allocatedInplace == NULL)
      for (long i=0; i < n; i++) V[i].del();
    else free(allocatedInplace);
    free(V);
    if(inEdges != NULL) free(inEdges);
  }
};

template <class vertex>
struct Uncompressed_Mem_Hypergraph : public Deletable {
public:
  vertex* V;
  vertex* H;
  long n_v;
  long m_v;
  long n_h;
  long m_h;
  void* edgesV, *inEdgesV, *edgesH, *inEdgesH;

 Uncompressed_Mem_Hypergraph(vertex* VV, vertex* HH, long nn_v, long mm_v, long nn_h, long mm_h, void* ai, void* _edgesV, void* _edgesH, void* _inEdgesV = NULL, void* _inEdgesH = NULL)
   : V(VV), H(HH), n_v(nn_v), m_v(mm_v), n_h(nn_h), m_h(mm_h), edgesV(_edgesV), edgesH(_edgesH), inEdgesV(_inEdgesV), inEdgesH(_inEdgesH) { }

  void del() {
    free(edgesV);
    free(edgesH);
    free(V);
    free(H);
    if(inEdgesV != NULL) free(inEdgesV);
    if(inEdgesH != NULL) free(inEdgesH);
  }
};

template <class vertex>
struct Compressed_Mem : public Deletable {
public:
  vertex* V;
  char* s;

  Compressed_Mem(vertex* _V, char* _s) :
                 V(_V), s(_s) { }

  void del() {
    free(V);
    free(s);
  }
};

template <class vertex>
struct graph {
  vertex *V;
  long n;
  long m;
  bool transposed;
  uintE* flags;
  Deletable *D;

graph(vertex* _V, long _n, long _m, Deletable* _D) : V(_V), n(_n), m(_m),
  D(_D), flags(NULL), transposed(0) {}

graph(vertex* _V, long _n, long _m, Deletable* _D, uintE* _flags) : V(_V),
  n(_n), m(_m), D(_D), flags(_flags), transposed(0) {}

  void del() {
    if (flags != NULL) free(flags);
    D->del();
    free(D);
  }

  void transpose() {
    if ((sizeof(vertex) == sizeof(asymmetricVertex)) ||
        (sizeof(vertex) == sizeof(compressedAsymmetricVertex))) {
      parallel_for(long i=0;i<n;i++) {
        V[i].flipEdges();
      }
      transposed = !transposed;
    }
  }
};

template <class vertex>
struct hypergraph {
  vertex *V;
  vertex *H;
  long n_v;
  long m_v;
  long n_h;
  long m_h;
  bool transposed;
  uintE* flagsV;
  uintE* flagsH;
  Deletable *D;

hypergraph(vertex* _V, vertex* _H, long _n_v, long _m_v, long _n_h, long _m_h, Deletable* _D) : V(_V), H(_H), n_v(_n_v), m_v(_m_v), n_h(_n_h), m_h(_m_h),
    D(_D), flagsV(NULL), flagsH(NULL), transposed(0) {}

hypergraph(vertex* _V, vertex* _H, long _n_v, long _m_v, long _n_h, long _m_h,  Deletable* _D, uintE* _flagsV, uintE* _flagsH) : V(_V), H(_H),
    n_v(_n_v), m_v(_m_v), n_h(_n_h), m_h(_m_h), D(_D), flagsV(_flagsV), flagsH(_flagsH), transposed(0) {}

  void del() {
    if (flagsV != NULL) free(flagsV);
    if (flagsH != NULL) free(flagsH); 
    D->del();
    free(D);
  }

  void transpose() {
    if ((sizeof(vertex) == sizeof(asymmetricVertex)) ||
        (sizeof(vertex) == sizeof(compressedAsymmetricVertex))) {
      parallel_for(long i=0;i<n_v;i++) {
        V[i].flipEdges();
      }
      parallel_for(long i=0;i<n_h;i++) {
	H[i].flipEdges();
      }
      transposed = !transposed;
    }
  }
};
#endif
