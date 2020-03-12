#include "header.h"
#include <map>
#include <set>
using namespace std;

/* #define NO_OPT */

extern "C" {

static int is_init=0;

  static map<int,set<void*> > free_list_gpu, free_list_host;
  static map<void *,int> live_ptrs_gpu, live_ptrs_host;

  static void clearGpuFreeList() {
    for(map<int,set<void*> >::iterator it=free_list_gpu.begin(); 
	it!=free_list_gpu.end(); ++it) {
      for(set<void*>::iterator it2=it->second.begin();
	  it2!=it->second.end(); ++it2) {
	hipFree(*it2);
      }
    }
    free_list_gpu.clear();
  }
  
  static void clearHostFreeList() {
    for(map<int,set<void*> >::iterator it=free_list_host.begin(); 
	it!=free_list_host.end(); ++it) {
      for(set<void*>::iterator it2=it->second.begin();
	  it2!=it->second.end(); ++it2) {
	hipHostFree(*it2);
      }
    }
    free_list_host.clear();
  }

  static int num_resurrections=0, num_morecore=0;

  typedef hipError_t (*mallocfn_t)(void **ptr, size_t bytes);
  static void *morecore(mallocfn_t fn, size_t bytes) {
    void *ptr;
    CUDA_SAFE(fn((void **)&ptr, bytes));
    num_morecore += 1;
    if(ptr==NULL) {
      /*try one more time*/
      clearHostFreeList();
      clearGpuFreeList();
      fn((void **)&ptr, bytes);
    }
    assert(ptr!=NULL); /*We hopefully have a pointer*/
    return ptr;
  }

  static inline void *resurrect_from_free_list(map<int,set<void *> > &free_map,
					size_t bytes, map<void*,int>& liveset) {
    void *ptr;
    num_resurrections +=1 ;
    assert(free_map.find(bytes) != free_map.end());
/*     assert(free_map.find(bytes)->second.size() > 0); */
    set<void *> &st = free_map.find(bytes)->second;
    ptr = *st.begin();
    st.erase(ptr);
    if(st.size()==0)
      free_map.erase(bytes);
    liveset[ptr] = bytes;
    return ptr;
  }

  void initmemmodule_()  {
    is_init=1;
  }

void *getGpuMem(size_t bytes) {
  assert(is_init);
  void *ptr;
#ifdef NO_OPT
  CUDA_SAFE(hipMalloc((void **) &ptr, bytes));
#else
  if(free_list_gpu.find(bytes)!=free_list_gpu.end()) {
    set<void*> &lst = free_list_gpu.find(bytes)->second;
    if(lst.size()!=0) {
      ptr = resurrect_from_free_list(free_list_gpu, bytes, live_ptrs_gpu);
      return ptr;
    }
  }
  else {
    for(map<int,set<void *> >::iterator it=free_list_gpu.begin();
	it != free_list_gpu.end(); ++it) {
      if(it->first >= bytes && it->second.size()>0) {
	ptr = resurrect_from_free_list(free_list_gpu, it->first, live_ptrs_gpu);
	return ptr;
      }
    }
  }
  ptr = morecore(hipMalloc, bytes);
/*   cutilSafeCall(hipMalloc((void **) &ptr, bytes)); */
  live_ptrs_gpu[ptr] = bytes;
#endif
  return ptr;
}

void *getHostMem(size_t bytes) {
  assert(is_init);
  void *ptr;
#ifdef NO_OPT
  CUDA_SAFE(hipHostMalloc((void **) &ptr, bytes));
#else
  if(free_list_host.find(bytes)!=free_list_host.end()) {
    set<void*> &lst = free_list_host.find(bytes)->second;
    if(lst.size()!=0) {
      ptr = resurrect_from_free_list(free_list_host, bytes, live_ptrs_host);
/*       ptr = *lst.begin(); */
/*       lst.erase(lst.begin()); */
/*       live_ptrs_host[ptr] = bytes; */
      return ptr;
    }
  }
  else {
    for(map<int,set<void *> >::iterator it=free_list_host.begin();
	it != free_list_host.end(); ++it) {
      if(it->first >= bytes && it->second.size()>0) {
	ptr = resurrect_from_free_list(free_list_host, it->first, live_ptrs_host);
/* 	set<void*> &lst = it->second; */
/* 	ptr = *lst.begin(); */
/* 	lst.erase(lst.begin()); */
/* 	live_ptrs_gpu[ptr] = bytes; */
	return ptr;
      }
    }
  }
/*   cutilSafeCall(hipHostMalloc((void **) &ptr, bytes)); */
  ptr = morecore(hipMallocHost, bytes);
  live_ptrs_host[ptr] = bytes;
#endif
  return ptr;
}

void freeHostMem(void *p) {
  int bytes;
  assert(is_init);
#ifdef NO_OPT
  hipHostFree(p);
#else
  assert(live_ptrs_host.find(p) != live_ptrs_host.end());
  bytes = live_ptrs_host[p];
  live_ptrs_host.erase(p);
  free_list_host[bytes].insert(p);
#endif
}

void freeGpuMem(void *p) {
  int bytes;
  assert(is_init);
#ifdef NO_OPT
  hipFree(p);
#else
  assert(live_ptrs_gpu.find(p) != live_ptrs_gpu.end());
  bytes = live_ptrs_gpu[p];
  live_ptrs_gpu.erase(p);
  free_list_gpu[bytes].insert(p);
#endif
}

void finalizememmodule_() {
  assert(is_init);
  is_init = 0;
  
  /*there should be no live pointers*/
  assert(live_ptrs_gpu.size()==0);
  assert(live_ptrs_host.size()==0);

  /*release all freed pointers*/
  clearGpuFreeList();
  clearHostFreeList();
  //printf("num. resurrections=%d \t num. morecore=%d\n", num_resurrections, num_morecore);
}

}

