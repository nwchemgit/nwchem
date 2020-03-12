#ifndef __header_h__
#define __header_h__

#ifdef TCE_HIP
#include <hip/hip_runtime_api.h>
#endif

//static int notset;
#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#ifdef TCE_CUDA
#ifdef OLD_CUDA
#include <cuda_runtime_api.h>
#else
#include <cuda.h>
#endif
#endif
#include <sys/types.h>
#include <sys/time.h>
#include <assert.h>
#include <time.h>
////#include "util.h"

#ifdef TCE_CUDA
#define CHECK_ERR(x) { \
    cudaError_t err = cudaGetLastError();\
    if (cudaSuccess != err) { \
        printf("%s\n",cudaGetErrorString(err)); \
        exit(1); \
    } } 

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) {\
    printf("CUDA CALL FAILED AT LINE %d OF FILE %s error %s\n", __LINE__, __FILE__, cudaGetErrorString(cudaGetLastError()) ); exit(1);}
#endif

#ifdef TCE_HIP
#define CHECK_ERR(x) { \
    hipError_t err = hipGetLastError();\
    if (hipSuccess != err) { \
        printf("%s\n",hipGetErrorString(err)); \
        exit(1); \
    } } 

#define CUDA_SAFE(x) if ( hipSuccess != (x) ) {\
    printf("HIP CALL FAILED AT LINE %d OF FILE %s error %s\n", __LINE__, __FILE__, hipGetErrorString(hipGetLastError()) ); exit(1);}
#endif
  
typedef long Integer;

#define DIV_UB(x,y) ((x)/(y)+((x)%(y)?1:0))
#define MIN(x,y) ((x)<(y)?(x):(y))

void initMemModule();
void *getGpuMem(size_t bytes);
void *getHostMem(size_t bytes);
void freeHostMem(void *p);
void freeGpuMem(void *p);
void finalizeMemModule();

#if defined(__cplusplus)
}
#endif

#endif /*__header_h__*/
/* $Id$ */
