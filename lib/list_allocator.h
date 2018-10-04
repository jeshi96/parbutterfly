// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2016 Guy Blelloch, Daniel Ferizovic, and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// A concurrent allocator for any fixed type T
// Keeps a local pool per processor
// Grabs list_size elements from a global pool if empty, and
// Returns list_size elements to the global pool when local pool=2*list_size
// Keeps track of number of allocated elements.
// Probably more efficient than a general purpose allocator

#pragma once

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <atomic>
#include "concurrent_stack.h"
#include "utilities.h"
#include "random_shuffle.h"

constexpr const size_t default_alloc_size = 1000000;
constexpr const size_t list_size = 1 << 16;
constexpr const size_t line_size = 64;

template <typename T>
class list_allocator {
  
  //union alignas(64) block {
  union block {
    T data;
    block* next;
  };

  struct thread_list {
    size_t sz;
    block* head;  
    block* mid;
    char cache_line[line_size];
  thread_list() : sz(0), head(NULL) {};
  };

  using block_p = block*;
    
  static block_p initialize_list(block_p);
  static block_p get_list();

 public:
  static bool initialized;
  static T* alloc();
  static void free(T*);
  static void init(size_t n = default_alloc_size,
		   bool randomize = 0,
		   size_t max_blocks = (((size_t) 1) << 32) - 1);
  static void reserve(size_t n = default_alloc_size, bool randomize=false);
  static void finish();
  static size_t block_size () {return _block_size;}
  static size_t num_allocated_blocks() {return blocks_allocated;}
  static size_t num_used_blocks();

 private:
  static void rand_shuffle();
  static concurrent_stack<block_p> pool_roots;
  static concurrent_stack<block_p> global_stack;
  static thread_list* local_lists;

  static int thread_count;
  static size_t list_length;
  static size_t max_blocks;
  static size_t _block_size;
  static std::atomic<size_t> blocks_allocated;
  static block_p allocate_blocks(size_t num_blocks);
};

template<typename T> concurrent_stack<typename list_allocator<T>::block_p>
list_allocator<T>::pool_roots;

template<typename T> concurrent_stack<typename list_allocator<T>::block_p>
list_allocator<T>::global_stack;

template<typename T> bool 
list_allocator<T>::initialized = false;

template<typename T> typename list_allocator<T>::thread_list*
list_allocator<T>::local_lists;

template<typename T> int 
list_allocator<T>::thread_count;

template<typename T> size_t 
list_allocator<T>::list_length = list_size;

template<typename T> size_t 
list_allocator<T>::max_blocks;

template<typename T> size_t 
list_allocator<T>::_block_size;

template<typename T> std::atomic<size_t>
list_allocator<T>::blocks_allocated;

// Allocate a new list of list_length elements
template<typename T>
auto list_allocator<T>::initialize_list(block_p start) -> block_p {
  block_p p = start;
  block_p end  = start + list_length - 1;

  while (p != end) {
    p->next = (p + 1);
    p++;
  }
  p->next = NULL;
  return start;
}

template<typename T>
size_t list_allocator<T>::num_used_blocks() {
  size_t free_blocks = global_stack.size()*list_length;
  for (int i = 0; i < thread_count; ++i) 
    free_blocks += local_lists[i].sz;
  return blocks_allocated - free_blocks;
}

template<typename T>
auto list_allocator<T>::allocate_blocks(size_t num_blocks) -> block_p { 
  block_p start = (block_p) aligned_alloc(line_size,
					  num_blocks * _block_size+ line_size);
  if (start == NULL) {
    fprintf(stderr, "Cannot allocate space in list_allocator");
    exit(1); }

  blocks_allocated += num_blocks; // atomic
  if (blocks_allocated > max_blocks) {
    fprintf(stderr, "Too many blocks in list_allocator, change max_blocks");
    exit(1);  }

  pool_roots.push(start); // keep track so can free later
  return start;
}

// Either grab a list from the global pool, or if there is none
// then allocate a new list
template<typename T>
auto list_allocator<T>::get_list() -> block_p {
    maybe<block_p> rem = global_stack.pop();
    if (rem) return *rem;
    block_p start = allocate_blocks(list_length);
    return initialize_list(start);
}

// Randomly orders the free blocks.  Only used for testing.
// Not safe if run concurrently with alloc and free
template<typename T>
void list_allocator<T>::rand_shuffle() {
  size_t num_free = num_allocated_blocks()-num_used_blocks();

  // pull all free blocks out
  T** P = new T*[num_free];
  cilk_for (int i=0; i < num_free; i ++)
    P[i] = alloc();

  // randomly shuffle them
  pbbs::random_shuffle(sequence<T*>(P,num_free));

  // put them back
  cilk_for (int i=0; i < num_free; i ++)
    free(P[i]);
  
  delete[] P; 
}

// Allocate n elements across however many lists are needed (rounded up)
template<typename T>
void list_allocator<T>::reserve(size_t n, bool randomize) {
  if (!initialized) init(n);
  else {
    size_t num_lists = thread_count + ceil(n / (double)list_length);
    block_p start = allocate_blocks(list_length*num_lists);
    cilk_for (int i = 0; i < num_lists; ++i) 
      global_stack.push(initialize_list(start + i*list_length));
    if (randomize) rand_shuffle();
  }
}

template<typename T>
void list_allocator<T>::init(size_t n, bool randomize, 
			     size_t _max_blocks) {
    if (initialized) return;
    initialized = true;
    blocks_allocated = 0;
    max_blocks = _max_blocks;
    thread_count = __cilkrts_get_nworkers();

    // Hack to account for possible allignment expansion
    // i.e. sizeof(T) might not work -- better way?
    block_p x;
    _block_size = (char*) (x+1) - (char*) x; 

    // reserve n blocks in the global pool
    reserve(n, randomize);

    // all local lists start out empty
    local_lists = new thread_list[thread_count];
}

template<typename T>
void list_allocator<T>::finish() {
    if (!initialized) return;

    delete[] local_lists;

    maybe<block_p> x;
    while (x = pool_roots.pop()) std::free(*x);
    pool_roots.clear();
    global_stack.clear();

    blocks_allocated = 0;
    initialized = false;
}

template<typename T>
void list_allocator<T>::free(T* node) {
    block_p new_node = (block_p) node;
    int id = __cilkrts_get_worker_number();

    if (local_lists[id].sz == list_length+1) {
      local_lists[id].mid = local_lists[id].head;
    } else if (local_lists[id].sz == 2*list_length) {
        global_stack.push(local_lists[id].mid->next);
        local_lists[id].mid->next = NULL;
        local_lists[id].sz = list_length;
    }

    new_node->next = local_lists[id].head;
    local_lists[id].head = new_node;
    local_lists[id].sz++;
}

template<typename T>
inline T* list_allocator<T>::alloc() {
    int id = __cilkrts_get_worker_number();

    if (!local_lists[id].sz)  {
      local_lists[id].head = get_list();
      local_lists[id].sz = list_length;
    }

    local_lists[id].sz--;
    block_p p = local_lists[id].head;
    local_lists[id].head = local_lists[id].head->next;

    return &p->data;
}

