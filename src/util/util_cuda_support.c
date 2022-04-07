#include <stdio.h>
#include <stdlib.h>
#define OLD_CUDA 1
#ifdef TCE_CUDA
#ifdef OLD_CUDA
#include <cuda_runtime_api.h>
#else
#include <cuda.h>
#endif
#endif
#ifdef TCE_HIP
#include <hip/hip_runtime_api.h>
#endif
#include "ga.h"
#include "typesf2c.h"
Integer FATR util_cuda_get_num_devices_(){
  int dev_count_check;
#ifdef TCE_CUDA
  cudaGetDeviceCount(&dev_count_check);
#endif
#ifdef TCE_HIP
  hipGetDeviceCount(&dev_count_check);
#endif
  return (Integer) dev_count_check;
}
/* $Id$ */
