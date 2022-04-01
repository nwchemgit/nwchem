/*------------------------------------------hybrid execution------------*/
/* $Id$ */
#include <assert.h>
///#define NUM_DEVICES 1
#define OLD_CUDA 1
static long long device_id=-1;
#include <stdio.h>
#include <stdlib.h>
#ifdef TCE_HIP
#include <hip/hip_runtime_api.h>
#endif
#ifdef TCE_CUDA
#ifdef OLD_CUDA
#include <cuda_runtime_api.h>
#else
#include <cuda.h>
#endif
#endif
#include "ga.h"
#include "typesf2c.h"

extern void FATR util_getppn_(Integer *);
extern int util_my_smp_index();

int check_device_(long *icuda) {
  /* Check whether this process is associated with a GPU */
  if((util_my_smp_index())<*icuda) return 1;
  return 0;
}

//void device_init_(int *icuda) {
int device_init_(long *icuda,long *cuda_device_number ) {
  /* Set device_id */
  
  int dev_count_check=0;
  device_id = util_my_smp_index();
#ifdef TCE_CUDA
  cudaGetDeviceCount(&dev_count_check);
#endif
#ifdef TCE_HIP
  hipGetDeviceCount(&dev_count_check);
#endif
  if(dev_count_check < *icuda){
    printf("Warning: Please check whether you have %ld cuda devices per node\n",*icuda);
    fflush(stdout);
    *cuda_device_number = 30;
  }
  else {
#ifdef TCE_CUDA
    cudaSetDevice(device_id);
#endif
#ifdef TCE_HIP
    hipSetDevice(device_id);
#endif
  }
  return 1;
}

