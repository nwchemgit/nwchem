#ifndef __header_h__
#define __header_h__
__device__ double* t3_s_d;
__device__ double* t3_d;
//static int notset;
#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <sys/types.h>
#include <sys/time.h>
#include <assert.h>
#include <time.h>
#include "cuda.h"
////#include "util.h"

#define CHECK_ERR(x) { \
    cudaError_t err = cudaGetLastError();\
    if (cudaSuccess != err) { \
        printf("%s\n",cudaGetErrorString(err)); \
        exit(1); \
    } } 

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) {\
   	printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}


typedef int Integer;

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
