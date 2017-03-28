#include <stdio.h>
#include <stdlib.h>
#ifdef OLD_CUDA
#include <cuda_runtime_api.h>
#else
#include <cuda.h>
#endif
#include "ga.h"
#include "typesf2c.h"
Integer FATR util_cuda_get_num_devices_(){
  int dev_count_check;
  cudaGetDeviceCount(&dev_count_check);
  return (Integer) dev_count_check;
}
/* $Id$ */
