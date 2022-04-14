__device__ double* t3_s_d;
__device__ double* t3_d;

#include "header.h"
extern    "C" void set_dev_mem_d(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d,size_t p6d)
{
    size_t size_t3;
    size_t3 = h1d*h2d*h3d*p4d*p5d*p6d;
    t3_d = (double *) getGpuMem(size_t3*sizeof(double));
    cudaMemset(t3_d,0,size_t3*sizeof(double));
}
extern          "C" void
dev_mem_d_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d)
{
    set_dev_mem_d((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d);
}
extern "C" void
dev_release()
{

        freeGpuMem(t3_d);
        freeGpuMem(t3_s_d);
}
extern "C" void
dev_release_()
{
    dev_release();
}


/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p6,p5,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_1_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t p6ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p6_0*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p6_0*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p6_0*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+p6_0*p6ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p6_0*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p6_0*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p6_0*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p6_0*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p6_0*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p6_0*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p6_1*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p6_1*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p6_1*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+p6_1*p6ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p6_1*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p6_1*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p6_1*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p6_1*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p6_1*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p6_1*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p6_2*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p6_2*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p6_2*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+p6_2*p6ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p6_2*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p6_2*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p6_2*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p6_2*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p6_2*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p6_2*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p6_3*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p6_3*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p6_3*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+p6_3*p6ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p6_3*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p6_3*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p6_3*p6ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p6_3*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p6_3*p6ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p6_3*p6ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_1_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p6ld_triplesx,p5ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*p6d*p5d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_1_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  p6ld_v2sub=h3d;
  h7ld_v2sub=p6d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  p6ld_triplesx=h1d*h3d;
  p5ld_triplesx=p6d*h1d*h3d;
  p4ld_triplesx=p5d*p6d*h1d*h3d;
  size_t total_x = h3d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_1_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p6ld_triplesx,p5ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_1_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_1_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,h2,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_2_kernel(size_t h1d,size_t h2d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t h2ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t h2ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  h2_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  h2_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  h2_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  h2_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+h2_0*h2ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+h2_1*h2ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+h2_2*h2ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+h2_3*h2ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+h2_0*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+h2_1*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+h2_2*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+h2_3*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_2_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h2d=h2d*p6d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*h2d*p5d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*h2d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_2_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  //CUDA_SAFE(
    cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice); //);
  //CUDA_SAFE(  
    cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice); //);
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  h2ld_v2sub=h3d;
  h7ld_v2sub=h2d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  h2ld_triplesx=h1d*h3d;
  p5ld_triplesx=h2d*h1d*h3d;
  p4ld_triplesx=p5d*h2d*h1d*h3d;
  size_t total_x = h3d*h2d*1;
  size_t total_y = p4d*p5d*h1d;
//printf("Blocks %d %d\n", total_x, total_y); 
//fflush(stdout);    
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_2_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,h7d,p4d,p5d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_2_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_2_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_3_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t h7ld_v2sub,size_t h1ld_triplesx,size_t h3ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  h3_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  h3_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  h3_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  h3_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
      triplesx_d[h1_3*h1ld_triplesx+h3_0*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
      triplesx_d[h1_3*h1ld_triplesx+h3_1*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
      triplesx_d[h1_3*h1ld_triplesx+h3_2*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
      triplesx_d[h1_3*h1ld_triplesx+h3_3*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_3_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  h3d=h3d*p6d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h1d*h3d*p5d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_3_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  h7ld_v2sub=h3d;
  h1ld_triplesx=1;
  h3ld_triplesx=h1d;
  p5ld_triplesx=h3d*h1d;
  p4ld_triplesx=p5d*h3d*h1d;
  size_t total_x = h3d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_3_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_3_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_3_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p5,p4,p6] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_4_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,size_t p6ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_4_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*p5d*p4d*p6d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_4_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  p6ld_v2sub=h3d;
  h7ld_v2sub=p6d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  p5ld_triplesx=h1d*h3d;
  p4ld_triplesx=p5d*h1d*h3d;
  p6ld_triplesx=p4d*p5d*h1d*h3d;
  size_t total_x = h3d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_4_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_4_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_4_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p4,p6] += t2sub[h7,p4,p5,h1] * v2sub[h3,h2,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_5_kernel(size_t h1d,size_t h2d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t h2ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t h2ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,size_t p6ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_0=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_1=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_2=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_3=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+h2_0*h2ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+h2_1*h2ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+h2_2*h2ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+h2_3*h2ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+h2_0*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+h2_1*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+h2_2*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+h2_3*h2ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_5_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*h2d*p5d*p4d*p6d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*h2d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_5_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  h2ld_v2sub=h3d;
  p6ld_v2sub=h2d*h3d;
  h7ld_v2sub=p6d*h2d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  h2ld_triplesx=h1d*h3d;
  p5ld_triplesx=h2d*h1d*h3d;
  p4ld_triplesx=p5d*h2d*h1d*h3d;
  p6ld_triplesx=p4d*p5d*h2d*h1d*h3d;
  size_t total_x = h3d*h2d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_5_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_5_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_5_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p4,p6] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_6_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h1ld_triplesx,size_t h3ld_triplesx,size_t p5ld_triplesx,size_t p4ld_triplesx,size_t p6ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal3;
      triplesx_d[h1_3*h1ld_triplesx+h3_0*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_0*p6ld_triplesx]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal7;
      triplesx_d[h1_3*h1ld_triplesx+h3_1*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_1*p6ld_triplesx]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal11;
      triplesx_d[h1_3*h1ld_triplesx+h3_2*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_2*p6ld_triplesx]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal15;
      triplesx_d[h1_3*h1ld_triplesx+h3_3*h3ld_triplesx+p5_3*p5ld_triplesx+p4_3*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p4_2*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p4_1*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p4_0*p4ld_triplesx+p6_3*p6ld_triplesx]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_6_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h1d*h3d*p5d*p4d*p6d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_6_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  p6ld_v2sub=h3d;
  h7ld_v2sub=p6d*h3d;
  h1ld_triplesx=1;
  h3ld_triplesx=h1d;
  p5ld_triplesx=h3d*h1d;
  p4ld_triplesx=p5d*h3d*h1d;
  p6ld_triplesx=p4d*p5d*h3d*h1d;
  size_t total_x = h3d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_6_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p4ld_triplesx,p6ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_6_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_6_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p5,p6,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_7_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t p5ld_triplesx,size_t p6ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p6_0*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p6_1*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p6_2*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+p5_3*p5ld_triplesx+p6_3*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_7_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*p5d*p6d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_7_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  p6ld_v2sub=h3d;
  h7ld_v2sub=p6d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  p5ld_triplesx=h1d*h3d;
  p6ld_triplesx=p5d*h1d*h3d;
  p4ld_triplesx=p6d*p5d*h1d*h3d;
  size_t total_x = h3d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_7_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_7_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_7_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p6,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,h2,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_8_kernel(size_t h1d,size_t h2d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t h2ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h3ld_triplesx,size_t h1ld_triplesx,size_t h2ld_triplesx,size_t p5ld_triplesx,size_t p6ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_0=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_1=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_2=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h2_3=rest_x%h2d;
  rest_x=rest_x/h2d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+h2_0*h2ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+h2_1*h2ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+h2_2*h2ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+h2_3*h2ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
      triplesx_d[h3_0*h3ld_triplesx+h1_3*h1ld_triplesx+h2_0*h2ld_triplesx+p5_3*p5ld_triplesx+p6_0*p6ld_triplesx+p4_3*p4ld_triplesx]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
      triplesx_d[h3_0*h3ld_triplesx+h1_2*h1ld_triplesx+h2_0*h2ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
      triplesx_d[h3_0*h3ld_triplesx+h1_1*h1ld_triplesx+h2_0*h2ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_0*h3ld_triplesx+h1_0*h1ld_triplesx+h2_0*h2ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
      triplesx_d[h3_1*h3ld_triplesx+h1_3*h1ld_triplesx+h2_1*h2ld_triplesx+p5_3*p5ld_triplesx+p6_1*p6ld_triplesx+p4_3*p4ld_triplesx]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
      triplesx_d[h3_1*h3ld_triplesx+h1_2*h1ld_triplesx+h2_1*h2ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
      triplesx_d[h3_1*h3ld_triplesx+h1_1*h1ld_triplesx+h2_1*h2ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_1*h3ld_triplesx+h1_0*h1ld_triplesx+h2_1*h2ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
      triplesx_d[h3_2*h3ld_triplesx+h1_3*h1ld_triplesx+h2_2*h2ld_triplesx+p5_3*p5ld_triplesx+p6_2*p6ld_triplesx+p4_3*p4ld_triplesx]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
      triplesx_d[h3_2*h3ld_triplesx+h1_2*h1ld_triplesx+h2_2*h2ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
      triplesx_d[h3_2*h3ld_triplesx+h1_1*h1ld_triplesx+h2_2*h2ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_2*h3ld_triplesx+h1_0*h1ld_triplesx+h2_2*h2ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
      triplesx_d[h3_3*h3ld_triplesx+h1_3*h1ld_triplesx+h2_3*h2ld_triplesx+p5_3*p5ld_triplesx+p6_3*p6ld_triplesx+p4_3*p4ld_triplesx]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
      triplesx_d[h3_3*h3ld_triplesx+h1_2*h1ld_triplesx+h2_3*h2ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
      triplesx_d[h3_3*h3ld_triplesx+h1_1*h1ld_triplesx+h2_3*h2ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h3_3*h3ld_triplesx+h1_0*h1ld_triplesx+h2_3*h2ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_8_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h3d*h1d*h2d*p5d*p6d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*h2d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_8_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  h2ld_v2sub=h3d;
  p6ld_v2sub=h2d*h3d;
  h7ld_v2sub=p6d*h2d*h3d;
  h3ld_triplesx=1;
  h1ld_triplesx=h3d;
  h2ld_triplesx=h1d*h3d;
  p5ld_triplesx=h2d*h1d*h3d;
  p6ld_triplesx=p5d*h2d*h1d*h3d;
  p4ld_triplesx=p6d*p5d*h2d*h1d*h3d;
  size_t total_x = h3d*h2d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_8_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,h2ld_v2sub,p6ld_v2sub,h7ld_v2sub,h3ld_triplesx,h1ld_triplesx,h2ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_8_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_8_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p6,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d1_9_kernel(size_t h1d,size_t h3d,size_t h7d,size_t p4d,size_t p5d,size_t p6d,size_t h7ld_t2sub,size_t p4ld_t2sub,size_t p5ld_t2sub,size_t h1ld_t2sub,size_t h3ld_v2sub,size_t p6ld_v2sub,size_t h7ld_v2sub,size_t h1ld_triplesx,size_t h3ld_triplesx,size_t p5ld_triplesx,size_t p6ld_triplesx,size_t p4ld_triplesx,double *triplesx_d, double *t2sub_d, double *v2sub_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h3_0,h3_1,h3_2,h3_3,h7,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,h7l,h7T;
  __shared__ double t2sub_shm[4*T1][Tcomm];
  __shared__ double v2sub_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_0=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_1=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_2=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  p5_3=rest_y%p5d;
  rest_y=rest_y/p5d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2sub_d_off, v2sub_d_off;for(h7T=0;h7T<h7d;h7T+=Tcomm){size_t h7l_hi;
    h7l_hi = MIN(Tcomm+h7T,h7d)-h7T;
    t2sub_d_off=p4_0*p4ld_t2sub+p5_0*p5ld_t2sub+h1_0*h1ld_t2sub;
    v2sub_d_off=h3_0*h3ld_v2sub+p6_0*p6ld_v2sub;
    if(thread_y+T1*0<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*0][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*0<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*0] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_1*p4ld_t2sub+p5_1*p5ld_t2sub+h1_1*h1ld_t2sub;
    v2sub_d_off=h3_1*h3ld_v2sub+p6_1*p6ld_v2sub;
    if(thread_y+T1*1<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*1][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*1<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*1] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_2*p4ld_t2sub+p5_2*p5ld_t2sub+h1_2*h1ld_t2sub;
    v2sub_d_off=h3_2*h3ld_v2sub+p6_2*p6ld_v2sub;
    if(thread_y+T1*2<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*2][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*2<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*2] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    t2sub_d_off=p4_3*p4ld_t2sub+p5_3*p5ld_t2sub+h1_3*h1ld_t2sub;
    v2sub_d_off=h3_3*h3ld_v2sub+p6_3*p6ld_v2sub;
    if(thread_y+T1*3<total_y)for(h7l=threadIdx.x;h7l<h7l_hi;h7l+=blockDim.x){
  h7=h7l+h7T;
  t2sub_shm[in1_idxl+T1*3][h7l] = t2sub_d[t2sub_d_off+h7*h7ld_t2sub];
      }
    if(thread_x+T1*3<total_x)for(h7l=threadIdx.y;h7l<h7l_hi;h7l+=blockDim.y){
  h7=h7l+h7T;
  v2sub_shm[h7l][in2_idxl+T1*3] = v2sub_d[v2sub_d_off+h7*h7ld_v2sub];
      }
    __syncthreads();
    for(h7l=0;h7l<h7l_hi;++h7l){
      a1=t2sub_shm[in1_idxl+T1*0][h7l];
      a2=t2sub_shm[in1_idxl+T1*1][h7l];
      a3=t2sub_shm[in1_idxl+T1*2][h7l];
      a4=t2sub_shm[in1_idxl+T1*3][h7l];
      b1=v2sub_shm[h7l][in2_idxl+T2*0];
      b2=v2sub_shm[h7l][in2_idxl+T2*1];
      b3=v2sub_shm[h7l][in2_idxl+T2*2];
      b4=v2sub_shm[h7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
      triplesx_d[h1_3*h1ld_triplesx+h3_0*h3ld_triplesx+p5_3*p5ld_triplesx+p6_0*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
      triplesx_d[h1_2*h1ld_triplesx+h3_0*h3ld_triplesx+p5_2*p5ld_triplesx+p6_0*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
      triplesx_d[h1_1*h1ld_triplesx+h3_0*h3ld_triplesx+p5_1*p5ld_triplesx+p6_0*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_0*h3ld_triplesx+p5_0*p5ld_triplesx+p6_0*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
      triplesx_d[h1_3*h1ld_triplesx+h3_1*h3ld_triplesx+p5_3*p5ld_triplesx+p6_1*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
      triplesx_d[h1_2*h1ld_triplesx+h3_1*h3ld_triplesx+p5_2*p5ld_triplesx+p6_1*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
      triplesx_d[h1_1*h1ld_triplesx+h3_1*h3ld_triplesx+p5_1*p5ld_triplesx+p6_1*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_1*h3ld_triplesx+p5_0*p5ld_triplesx+p6_1*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
      triplesx_d[h1_3*h1ld_triplesx+h3_2*h3ld_triplesx+p5_3*p5ld_triplesx+p6_2*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
      triplesx_d[h1_2*h1ld_triplesx+h3_2*h3ld_triplesx+p5_2*p5ld_triplesx+p6_2*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
      triplesx_d[h1_1*h1ld_triplesx+h3_2*h3ld_triplesx+p5_1*p5ld_triplesx+p6_2*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_2*h3ld_triplesx+p5_0*p5ld_triplesx+p6_2*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
      triplesx_d[h1_3*h1ld_triplesx+h3_3*h3ld_triplesx+p5_3*p5ld_triplesx+p6_3*p6ld_triplesx+p4_3*p4ld_triplesx]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
      triplesx_d[h1_2*h1ld_triplesx+h3_3*h3ld_triplesx+p5_2*p5ld_triplesx+p6_3*p6ld_triplesx+p4_2*p4ld_triplesx]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
      triplesx_d[h1_1*h1ld_triplesx+h3_3*h3ld_triplesx+p5_1*p5ld_triplesx+p6_3*p6ld_triplesx+p4_1*p4ld_triplesx]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      triplesx_d[h1_0*h1ld_triplesx+h3_3*h3ld_triplesx+p5_0*p5ld_triplesx+p6_3*p6ld_triplesx+p4_0*p4ld_triplesx]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d1_9_cuda(size_t h1d, size_t h2d, size_t h3d, size_t h7d, size_t p4d, size_t p5d, size_t p6d, double *triplesx, double *t2sub, double *v2sub) {
  h3d=h3d*h2d;
  size_t h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx;
  size_t size_triplesx,size_block_triplesx,size_t2sub,size_v2sub;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2sub_d,*v2sub_d;
  size_triplesx=h1d*h3d*p5d*p6d*p4d*sizeof(double);
  size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
  size_v2sub=h3d*p6d*h7d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d1_9_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_triplesx=size_triplesx/nstreams;
  t2sub_d=(double*)getGpuMem(size_t2sub);
  v2sub_d=(double*)getGpuMem(size_v2sub);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2sub_d,t2sub,size_t2sub,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2sub_d,v2sub,size_v2sub,cudaMemcpyHostToDevice));
  h7ld_t2sub=1;
  p4ld_t2sub=h7d;
  p5ld_t2sub=p4d*h7d;
  h1ld_t2sub=p5d*p4d*h7d;
  h3ld_v2sub=1;
  p6ld_v2sub=h3d;
  h7ld_v2sub=p6d*h3d;
  h1ld_triplesx=1;
  h3ld_triplesx=h1d;
  p5ld_triplesx=h3d*h1d;
  p6ld_triplesx=p5d*h3d*h1d;
  p4ld_triplesx=p6d*p5d*h3d*h1d;
  size_t total_x = h3d*p6d*1;
  size_t total_y = p4d*p5d*h1d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d1_9_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h3d,h7d,p4d,p5d,p6d,h7ld_t2sub,p4ld_t2sub,p5ld_t2sub,h1ld_t2sub,h3ld_v2sub,p6ld_v2sub,h7ld_v2sub,h1ld_triplesx,h3ld_triplesx,p5ld_triplesx,p6ld_triplesx,p4ld_triplesx,t3_d,t2sub_d,v2sub_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  freeGpuMem(t2sub_d);
  freeGpuMem(v2sub_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d1_9_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_9_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*h7d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,triplesx,t2sub,v2sub);
}

/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_1_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]-=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3]-=tlocal3;
      t3d[h3_0*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_0*p6ld_t3+p4_3*p4ld_t3]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]-=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]-=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3]-=tlocal7;
      t3d[h3_1*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_1*p6ld_t3+p4_3*p4ld_t3]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]-=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]-=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3]-=tlocal11;
      t3d[h3_2*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_2*p6ld_t3+p4_3*p4ld_t3]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]-=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]-=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3]-=tlocal15;
      t3d[h3_3*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_3*p6ld_t3+p4_3*p4ld_t3]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]-=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_1_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  p6d=p6d*p5d;
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h3d*h2d*h1d*p6d*p4d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_1_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  h3ld_t3=1;
  h2ld_t3=h3d;
  h1ld_t3=h2d*h3d;
  p6ld_t3=h1d*h2d*h3d;
  p4ld_t3=p6d*h1d*h2d*h3d;
  size_t total_x = h3d*p6d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_1_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_1_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_1_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4] -= t2[p7,p4,h1,h2] * v2[p7,h3]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_2_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t h2ld_t3,size_t h1ld_t3,size_t h3ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_0=rest_y;
  h3_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_1=rest_y;
  h3_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_2=rest_y;
  h3_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_3=rest_y;
  h3_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3]-=tlocal3;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_0*h3ld_t3+p4_3*p4ld_t3]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3]-=tlocal7;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_1*h3ld_t3+p4_3*p4ld_t3]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3]-=tlocal11;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_2*h3ld_t3+p4_3*p4ld_t3]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3]-=tlocal15;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_3*h3ld_t3+p4_3*p4ld_t3]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_2_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  h3d=h3d*p6d;
  h3d=h3d*p5d;
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h1d*h3d*p4d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_2_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  h2ld_t3=1;
  h1ld_t3=h2d;
  h3ld_t3=h1d*h2d;
  p4ld_t3=h3d*h1d*h2d;
  size_t total_x = h3d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_2_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_2_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_2_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p6,p4] += t2[p7,p4,h1,h2] * v2[p7,h3,p6]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_3_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t h2ld_t3,size_t h3ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_0=rest_y;
  p6_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_1=rest_y;
  p6_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_2=rest_y;
  p6_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p4_3=rest_y;
  p6_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3]+=tlocal3;
      t3d[h2_3*h2ld_t3+h3_0*h3ld_t3+h1_3*h1ld_t3+p6_0*p6ld_t3+p4_3*p4ld_t3]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3]+=tlocal7;
      t3d[h2_3*h2ld_t3+h3_1*h3ld_t3+h1_3*h1ld_t3+p6_1*p6ld_t3+p4_3*p4ld_t3]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3]+=tlocal11;
      t3d[h2_3*h2ld_t3+h3_2*h3ld_t3+h1_3*h1ld_t3+p6_2*p6ld_t3+p4_3*p4ld_t3]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3]+=tlocal15;
      t3d[h2_3*h2ld_t3+h3_3*h3ld_t3+h1_3*h1ld_t3+p6_3*p6ld_t3+p4_3*p4ld_t3]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_3_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  p6d=p6d*p5d;
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p6ld_t3,p4ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h3d*h1d*p6d*p4d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_3_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  h2ld_t3=1;
  h3ld_t3=h2d;
  h1ld_t3=h3d*h2d;
  p6ld_t3=h1d*h3d*h2d;
  p4ld_t3=p6d*h1d*h3d*h2d;
  size_t total_x = h3d*p6d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_3_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_3_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_3_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_4_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_0=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_1=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_2=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_3=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]+=tlocal3;
      t3d[h3_0*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_0*p6ld_t3+p4_3*p4ld_t3+p5_0*p5ld_t3]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]+=tlocal7;
      t3d[h3_1*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_1*p6ld_t3+p4_3*p4ld_t3+p5_1*p5ld_t3]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]+=tlocal11;
      t3d[h3_2*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_2*p6ld_t3+p4_3*p4ld_t3+p5_2*p5ld_t3]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]+=tlocal15;
      t3d[h3_3*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p6_3*p6ld_t3+p4_3*p4ld_t3+p5_3*p5ld_t3]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_4_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h3d*h2d*h1d*p6d*p4d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_4_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  p5ld_v2=p6d*h3d*p7d;
  h3ld_t3=1;
  h2ld_t3=h3d;
  h1ld_t3=h2d*h3d;
  p6ld_t3=h1d*h2d*h3d;
  p4ld_t3=p6d*h1d*h2d*h3d;
  p5ld_t3=p4d*p6d*h1d*h2d*h3d;
  size_t total_x = h3d*p6d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_4_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_4_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_4_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_5_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p5ld_v2,size_t h2ld_t3,size_t h1ld_t3,size_t h3ld_t3,size_t p4ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]+=tlocal3;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_0*h3ld_t3+p4_3*p4ld_t3+p5_0*p5ld_t3]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]+=tlocal7;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_1*h3ld_t3+p4_3*p4ld_t3+p5_1*p5ld_t3]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]+=tlocal11;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_2*h3ld_t3+p4_3*p4ld_t3+p5_2*p5ld_t3]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]+=tlocal15;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_3*h3ld_t3+p4_3*p4ld_t3+p5_3*p5ld_t3]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_5_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  h3d=h3d*p6d;
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p5ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h1d*h3d*p4d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_5_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p5ld_v2=h3d*p7d;
  h2ld_t3=1;
  h1ld_t3=h2d;
  h3ld_t3=h1d*h2d;
  p4ld_t3=h3d*h1d*h2d;
  p5ld_t3=p4d*h3d*h1d*h2d;
  size_t total_x = h3d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_5_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p5ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_5_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_5_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p6,p4,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_6_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h2ld_t3,size_t h3ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_0=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_1=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_2=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_3=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]-=tlocal3;
      t3d[h2_3*h2ld_t3+h3_0*h3ld_t3+h1_3*h1ld_t3+p6_0*p6ld_t3+p4_3*p4ld_t3+p5_0*p5ld_t3]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p6_0*p6ld_t3+p4_2*p4ld_t3+p5_0*p5ld_t3]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p6_0*p6ld_t3+p4_1*p4ld_t3+p5_0*p5ld_t3]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p6_0*p6ld_t3+p4_0*p4ld_t3+p5_0*p5ld_t3]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]-=tlocal7;
      t3d[h2_3*h2ld_t3+h3_1*h3ld_t3+h1_3*h1ld_t3+p6_1*p6ld_t3+p4_3*p4ld_t3+p5_1*p5ld_t3]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p6_1*p6ld_t3+p4_2*p4ld_t3+p5_1*p5ld_t3]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p6_1*p6ld_t3+p4_1*p4ld_t3+p5_1*p5ld_t3]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p6_1*p6ld_t3+p4_0*p4ld_t3+p5_1*p5ld_t3]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]-=tlocal11;
      t3d[h2_3*h2ld_t3+h3_2*h3ld_t3+h1_3*h1ld_t3+p6_2*p6ld_t3+p4_3*p4ld_t3+p5_2*p5ld_t3]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p6_2*p6ld_t3+p4_2*p4ld_t3+p5_2*p5ld_t3]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p6_2*p6ld_t3+p4_1*p4ld_t3+p5_2*p5ld_t3]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p6_2*p6ld_t3+p4_0*p4ld_t3+p5_2*p5ld_t3]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]-=tlocal15;
      t3d[h2_3*h2ld_t3+h3_3*h3ld_t3+h1_3*h1ld_t3+p6_3*p6ld_t3+p4_3*p4ld_t3+p5_3*p5ld_t3]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p6_3*p6ld_t3+p4_2*p4ld_t3+p5_3*p5ld_t3]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p6_3*p6ld_t3+p4_1*p4ld_t3+p5_3*p5ld_t3]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p6_3*p6ld_t3+p4_0*p4ld_t3+p5_3*p5ld_t3]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_6_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h3d*h1d*p6d*p4d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_6_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  p5ld_v2=p6d*h3d*p7d;
  h2ld_t3=1;
  h3ld_t3=h2d;
  h1ld_t3=h3d*h2d;
  p6ld_t3=h1d*h3d*h2d;
  p4ld_t3=p6d*h1d*h3d*h2d;
  p5ld_t3=p4d*p6d*h1d*h3d*h2d;
  size_t total_x = h3d*p6d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_6_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR();
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_6_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_6_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p4,p6,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_7_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p4ld_t3,size_t p6ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_0=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_1=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_2=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_3=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal3;
      t3d[h3_0*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h3_0*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h3_0*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_0*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal7;
      t3d[h3_1*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h3_1*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h3_1*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_1*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal11;
      t3d[h3_2*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h3_2*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h3_2*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_2*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal15;
      t3d[h3_3*h3ld_t3+h2_3*h2ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h3_3*h3ld_t3+h2_2*h2ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h3_3*h3ld_t3+h2_1*h2ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h3_3*h3ld_t3+h2_0*h2ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_7_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p4ld_t3,p6ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h3d*h2d*h1d*p4d*p6d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_7_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  p5ld_v2=p6d*h3d*p7d;
  h3ld_t3=1;
  h2ld_t3=h3d;
  h1ld_t3=h2d*h3d;
  p4ld_t3=h1d*h2d*h3d;
  p6ld_t3=p4d*h1d*h2d*h3d;
  p5ld_t3=p6d*p4d*h1d*h2d*h3d;
  size_t total_x = h3d*p6d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_7_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p4ld_t3,p6ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_7_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_7_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4,p6,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_8_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h2ld_t3,size_t h1ld_t3,size_t h3ld_t3,size_t p4ld_t3,size_t p6ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  p6_0=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  p6_1=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  p6_2=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  p6_3=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal3;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_0*h3ld_t3+p4_3*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_0*h3ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_0*h3ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_0*h3ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]-=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal7;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_1*h3ld_t3+p4_3*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_1*h3ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_1*h3ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_1*h3ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]-=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal11;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_2*h3ld_t3+p4_3*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_2*h3ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_2*h3ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_2*h3ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]-=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal15;
      t3d[h2_3*h2ld_t3+h1_3*h1ld_t3+h3_3*h3ld_t3+p4_3*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
      t3d[h2_2*h2ld_t3+h1_2*h1ld_t3+h3_3*h3ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
      t3d[h2_1*h2ld_t3+h1_1*h1ld_t3+h3_3*h3ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h1_0*h1ld_t3+h3_3*h3ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]-=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_8_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3,p6ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h1d*h3d*p4d*p6d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_8_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  p5ld_v2=p6d*h3d*p7d;
  h2ld_t3=1;
  h1ld_t3=h2d;
  h3ld_t3=h1d*h2d;
  p4ld_t3=h3d*h1d*h2d;
  p6ld_t3=p4d*h3d*h1d*h2d;
  p5ld_t3=p6d*p4d*h3d*h1d*h2d;
  size_t total_x = h3d*p6d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_8_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h1ld_t3,h3ld_t3,p4ld_t3,p6ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern "C" void sd_t_d2_8_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_8_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p4,p6,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_d2_9_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p7d,size_t p7ld_t2,size_t p4ld_t2,size_t h1ld_t2,size_t h2ld_t2,size_t p7ld_v2,size_t h3ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h2ld_t3,size_t h3ld_t3,size_t h1ld_t3,size_t p4ld_t3,size_t p6ld_t3,size_t p5ld_t3,double *t3d, double *t2_d, double *v2_d,size_t unused_idx, size_t total_x, size_t total_y) {
  size_t h1_0,h1_1,h1_2,h1_3,h2_0,h2_1,h2_2,h2_3,h3_0,h3_1,h3_2,h3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3,p7;
  double a1,b1;
  double a2,b2;
  double a3,b3;
  double a4,b4;
  size_t in1_idxl,in2_idxl,p7l,p7T;
  __shared__ double t2_shm[4*T1][Tcomm];
  __shared__ double v2_shm[Tcomm][4*T2];
  size_t rest_x=blockIdx.x;
  size_t rest_y=blockIdx.y;
  size_t thread_x = T2*4 * rest_x + threadIdx.x;
  size_t thread_y = T1*4 * rest_y + threadIdx.y;
  in1_idxl=threadIdx.y;
  in2_idxl=threadIdx.x ;
  double tlocal1=0;
  double tlocal2=0;
  double tlocal3=0;
  double tlocal4=0;
  double tlocal5=0;
  double tlocal6=0;
  double tlocal7=0;
  double tlocal8=0;
  double tlocal9=0;
  double tlocal10=0;
  double tlocal11=0;
  double tlocal12=0;
  double tlocal13=0;
  double tlocal14=0;
  double tlocal15=0;
  double tlocal16=0;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*0;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*0;
  h2_0=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_0=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_0=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_0=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_0=rest_y;
  p5_0=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*1;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*1;
  h2_1=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_1=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_1=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_1=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_1=rest_y;
  p5_1=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*2;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*2;
  h2_2=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_2=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_2=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_2=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_2=rest_y;
  p5_2=rest_x;
  rest_x = T2 *4* blockIdx.x + threadIdx.x+T1*3;
  rest_y = T1 *4* blockIdx.y + threadIdx.y+T1*3;
  h2_3=rest_y%h2d;
  rest_y=rest_y/h2d;
  h3_3=rest_x%h3d;
  rest_x=rest_x/h3d;
  h1_3=rest_y%h1d;
  rest_y=rest_y/h1d;
  p6_3=rest_x%p6d;
  rest_x=rest_x/p6d;
  p4_3=rest_y;
  p5_3=rest_x;
  size_t t2_d_off, v2_d_off;for(p7T=0;p7T<p7d;p7T+=Tcomm){size_t p7l_hi;
    p7l_hi = MIN(Tcomm+p7T,p7d)-p7T;
    t2_d_off=p4_0*p4ld_t2+h1_0*h1ld_t2+h2_0*h2ld_t2;
    v2_d_off=h3_0*h3ld_v2+p6_0*p6ld_v2+p5_0*p5ld_v2;
    if(thread_y+T1*0<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*0][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*0<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*0] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_1*p4ld_t2+h1_1*h1ld_t2+h2_1*h2ld_t2;
    v2_d_off=h3_1*h3ld_v2+p6_1*p6ld_v2+p5_1*p5ld_v2;
    if(thread_y+T1*1<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*1][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*1<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*1] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_2*p4ld_t2+h1_2*h1ld_t2+h2_2*h2ld_t2;
    v2_d_off=h3_2*h3ld_v2+p6_2*p6ld_v2+p5_2*p5ld_v2;
    if(thread_y+T1*2<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*2][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*2<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*2] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    t2_d_off=p4_3*p4ld_t2+h1_3*h1ld_t2+h2_3*h2ld_t2;
    v2_d_off=h3_3*h3ld_v2+p6_3*p6ld_v2+p5_3*p5ld_v2;
    if(thread_y+T1*3<total_y)for(p7l=threadIdx.x;p7l<p7l_hi;p7l+=blockDim.x){
  p7=p7l+p7T;
  t2_shm[in1_idxl+T1*3][p7l] = t2_d[t2_d_off+p7*p7ld_t2];
      }
    if(thread_x+T1*3<total_x)for(p7l=threadIdx.y;p7l<p7l_hi;p7l+=blockDim.y){
  p7=p7l+p7T;
  v2_shm[p7l][in2_idxl+T1*3] = v2_d[v2_d_off+p7*p7ld_v2];
      }
    __syncthreads();
    for(p7l=0;p7l<p7l_hi;++p7l){
      a1=t2_shm[in1_idxl+T1*0][p7l];
      a2=t2_shm[in1_idxl+T1*1][p7l];
      a3=t2_shm[in1_idxl+T1*2][p7l];
      a4=t2_shm[in1_idxl+T1*3][p7l];
      b1=v2_shm[p7l][in2_idxl+T2*0];
      b2=v2_shm[p7l][in2_idxl+T2*1];
      b3=v2_shm[p7l][in2_idxl+T2*2];
      b4=v2_shm[p7l][in2_idxl+T2*3];
      tlocal1+=a1*b1;
      tlocal2+=a2*b1;
      tlocal3+=a3*b1;
      tlocal4+=a4*b1;
      tlocal5+=a1*b2;
      tlocal6+=a2*b2;
      tlocal7+=a3*b2;
      tlocal8+=a4*b2;
      tlocal9+=a1*b3;
      tlocal10+=a2*b3;
      tlocal11+=a3*b3;
      tlocal12+=a4*b3;
      tlocal13+=a1*b4;
      tlocal14+=a2*b4;
      tlocal15+=a3*b4;
      tlocal16+=a4*b4;
    }
    __syncthreads();
  }
  if(thread_x+T1*0<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal3;
      t3d[h2_3*h2ld_t3+h3_0*h3ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal4;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal2;
      t3d[h2_2*h2ld_t3+h3_0*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal3;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal1;
      t3d[h2_1*h2ld_t3+h3_0*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal2;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_0*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_0*p6ld_t3+p5_0*p5ld_t3]+=tlocal1;
    }
  }
  if(thread_x+T1*1<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal7;
      t3d[h2_3*h2ld_t3+h3_1*h3ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal8;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal6;
      t3d[h2_2*h2ld_t3+h3_1*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal7;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal5;
      t3d[h2_1*h2ld_t3+h3_1*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal6;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_1*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_1*p6ld_t3+p5_1*p5ld_t3]+=tlocal5;
    }
  }
  if(thread_x+T1*2<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal11;
      t3d[h2_3*h2ld_t3+h3_2*h3ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal12;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal10;
      t3d[h2_2*h2ld_t3+h3_2*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal11;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal9;
      t3d[h2_1*h2ld_t3+h3_2*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal10;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_2*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_2*p6ld_t3+p5_2*p5ld_t3]+=tlocal9;
    }
  }
  if(thread_x+T1*3<total_x){
    if(thread_y+T2*3<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal15;
      t3d[h2_3*h2ld_t3+h3_3*h3ld_t3+h1_3*h1ld_t3+p4_3*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal16;
    }
    else if(thread_y+T2*2<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal14;
      t3d[h2_2*h2ld_t3+h3_3*h3ld_t3+h1_2*h1ld_t3+p4_2*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal15;
    }
    else if(thread_y+T2*1<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal13;
      t3d[h2_1*h2ld_t3+h3_3*h3ld_t3+h1_1*h1ld_t3+p4_1*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal14;
    }
    else if(thread_y+T2*0<total_y) {
      t3d[h2_0*h2ld_t3+h3_3*h3ld_t3+h1_0*h1ld_t3+p4_0*p4ld_t3+p6_3*p6ld_t3+p5_3*p5ld_t3]+=tlocal13;
    }
  }
  __syncthreads();
}
extern "C" void sd_t_d2_9_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, size_t p7d, double *t3, double *t2, double *v2) {
  size_t p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p4ld_t3,p6ld_t3,p5ld_t3;
  size_t size_t3,size_block_t3,size_el_block_t3,size_t2,size_v2;
  cudaStream_t *streams;
  size_t nstreams,i;
  double *t2_d,*v2_d;
  size_t3=h2d*h3d*h1d*p4d*p6d*p5d*sizeof(double);
  size_t2=p7d*p4d*h1d*h2d*sizeof(double);
  size_v2=p7d*h3d*p6d*p5d*sizeof(double);
  cudaFuncSetCacheConfig(sd_t_d2_9_kernel, cudaFuncCachePreferShared);
  nstreams=1;
  size_block_t3=size_t3/nstreams;
  size_el_block_t3=size_block_t3/sizeof(double);
  //t3d=(double*)getGpuMem(size_t3);
  t2_d=(double*)getGpuMem(size_t2);
  v2_d=(double*)getGpuMem(size_v2);
  streams=(cudaStream_t*) malloc(nstreams*sizeof(cudaStream_t));
  assert(streams!= NULL);
  for(i=0;i<nstreams;++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i])) ;
  }
  CUDA_SAFE(cudaMemcpy(t2_d,t2,size_t2,cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d,v2,size_v2,cudaMemcpyHostToDevice));
  p7ld_t2=1;
  p4ld_t2=p7d;
  h1ld_t2=p4d*p7d;
  h2ld_t2=h1d*p4d*p7d;
  p7ld_v2=1;
  h3ld_v2=p7d;
  p6ld_v2=h3d*p7d;
  p5ld_v2=p6d*h3d*p7d;
  h2ld_t3=1;
  h3ld_t3=h2d;
  h1ld_t3=h3d*h2d;
  p4ld_t3=h1d*h3d*h2d;
  p6ld_t3=p4d*h1d*h3d*h2d;
  p5ld_t3=p6d*p4d*h1d*h3d*h2d;
  size_t total_x = h3d*p6d*p5d;
  size_t total_y = p4d*h1d*h2d;
  dim3 dimBlock(T2,T1);dim3 dimGrid(DIV_UB(total_x,(4*T2)), DIV_UB(total_y,(4*T1)));
  for(i=0;i<nstreams;++i){
    sd_t_d2_9_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p7d,p7ld_t2,p4ld_t2,h1ld_t2,h2ld_t2,p7ld_v2,h3ld_v2,p6ld_v2,p5ld_v2,h2ld_t3,h3ld_t3,h1ld_t3,p4ld_t3,p6ld_t3,p5ld_t3,t3_d,t2_d,v2_d,i,total_x,total_y);
    CHECK_ERR("Kernel execution failed");
  }
  cudaDeviceSynchronize();
  for(i=0;i<nstreams;++i){
    cudaStreamDestroy(streams[i]);}
  //freeGpuMem(t3d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  free(streams);
}

extern "C" void sd_t_d2_9_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_9_cuda((size_t)*h1d,(size_t)*h2d,(size_t)*h3d,(size_t)*p4d,(size_t)*p5d,(size_t)*p6d,(size_t)*p7d,t3,t2,v2);
}


#define MAX_h3 64
/* IMPORTANT!!!!
t3_d must be passed as parameter to kernel function. A __global__ function can't access the global variable directly*/

__global__ void compute_energy_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,double* eval1,double* eval2,double* eval3,double* eval4,double* eval5,double* eval6, double* energy, double factor, size_t total_size, double* t3d, double* t3_sd)
{
  size_t h1,h2,p6,p4,p5, h3,i=0;
  double e1,e2,e4,e5,e6;
//  __shared__ double t2_shm[MAX_h3];
  __shared__ double energy_s[T1];
  __shared__ double energy2_s[T1];
  double inner_fac;
  size_t limit;
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  if(threadIdx.x==0)
  {
        energy[blockIdx.x]=0;
        energy[blockIdx.x+gridDim.x]=0;
        energy_s[threadIdx.x] = 0.0;
        energy2_s[threadIdx.x] = 0.0;
  }

  for(size_t j =0; j<T2*T1;j++) {
    thread_x = T2*T1*blockIdx.x + j;  
    rest_x = thread_x;
    __syncthreads();
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    h1=rest_x%h1d;
    rest_x=rest_x/h1d;
    p6=rest_x%p6d;
    rest_x=rest_x/p6d;
    p5=rest_x%p5d;
    rest_x=rest_x/p5d;
    p4=rest_x%p4d;
    e1 = eval1[h1];
    e2 = eval2[h2];
    e4 = eval4[p4];
    e5 = eval5[p5];
    e6 = eval6[p6];
/*
  for(p4=0;p4<p4d;p4++) 
    for(p5 = 0;p5<p5d;p5++)
        for(p6=0;p6<p6d;p6++) 
            for(h1= 0;h1<h1d;h1++) 
                for(h2=0;h2<h2d;h2++) 
                    for(h3=0;h3<h3d;h3++) {
                        inner_fac = -eval4[p4]-eval5[p5]-eval6[p6]+eval1[h1]
                            +eval2[h2]+eval3[h3];
                        energy_s[0]+=factor*t3d[i]*t3d[i]/inner_fac;
                        energy2_s[0]+=factor*t3d[i]*(t3_sd[i]+t3d[i])/inner_fac;
                        i++;
                    }
*/
    if(thread_x<total_size)
    for(size_t i=0;i<h3d;i++)
    {
        inner_fac = -e4-e5-e6+e1+e2+eval3[i]; //t2_shm[i];
//ckbn avoid e1 in case we need just (T)
        energy_s[threadIdx.x] += factor* t3d[thread_x*h3d+i]*t3d[thread_x*h3d+i]/inner_fac;
        energy2_s[threadIdx.x] += factor* t3d[thread_x*h3d+i]*(t3_sd[thread_x*h3d+i]+t3d[thread_x*h3d+i])/inner_fac;
    }
    __syncthreads();
  }
  if(threadIdx.x==0)
  {
/*    limit = blockDim.x;
      if (blockIdx.x == (gridDim.x-1)) limit = total_size%blockDim.x;
      for(size_t i=0;i<limit;i++)
      {
        energy[blockIdx.x]+=energy_s[i];
        energy[blockIdx.x+gridDim.x]+=energy2_s[i];
      }
*/
    energy[blockIdx.x] = energy_s[0];
    energy[blockIdx.x+gridDim.x] = energy2_s[0];
   }
  __syncthreads();

}

extern   "C" void compute_energy(double factor, double* energy, double* eval1, double* eval2,double* eval3,double* eval4,double* eval5,double* eval6,size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d,size_t p6d, double* host1, double* host2)
//ckbn en_comment, double* total_d, double* total_s)
{
    double* energy_d, *energy_h;
    double* eval_d1,*eval_d2,*eval_d3,*eval_d4,*eval_d5,*eval_d6;
    size_t size_energy = 2*sizeof(double);
    size_t total_block = DIV_UB((h1d*h2d*p4d*p5d*p6d), (T2*T1));

//    size_t total_block = 1;
    size_t total_elements = h1d*h2d*p4d*p5d*p6d;

    energy_d = (double*)getGpuMem(size_energy*total_block*2);
    size_t i=0,in; 
    double* t3 = (double*)malloc(sizeof(double)*h3d*total_elements);
    double* ts3 = (double*)malloc(sizeof(double)*h3d*total_elements);

    energy_h = (double*)getHostMem(size_energy*2*total_block);
    eval_d1 = (double*)getGpuMem(h1d*sizeof(double));
    eval_d2 = (double*)getGpuMem(h2d*sizeof(double));
    eval_d3 = (double*)getGpuMem(h3d*sizeof(double));
    eval_d4 = (double*)getGpuMem(p4d*sizeof(double));
    eval_d5 = (double*)getGpuMem(p5d*sizeof(double));
    eval_d6 = (double*)getGpuMem(p6d*sizeof(double));

    CUDA_SAFE(cudaMemcpy(eval_d1, eval1, h1d*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE(cudaMemcpy(eval_d2, eval2, h2d*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE(cudaMemcpy(eval_d3, eval3, h3d*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE(cudaMemcpy(eval_d4, eval4, p4d*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE(cudaMemcpy(eval_d5, eval5, p5d*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE(cudaMemcpy(eval_d6, eval6, p6d*sizeof(double), cudaMemcpyHostToDevice));
/* for test only */
//printf("host 2 is %f %f\n", host2[0], host2[1]);
//    CUDA_SAFE(cudaMemcpy(t3_s_d, host2, total_elements*h3d*sizeof(double), cudaMemcpyHostToDevice));

    dim3 dimBlock(1); //T2*T1);
    dim3 dimGrid(total_block);
    compute_energy_kernel<<<dimGrid,dimBlock,0>>>(h1d,h2d,h3d,p4d,p5d,p6d, eval_d1,eval_d2,eval_d3,eval_d4,eval_d5,eval_d6,energy_d, factor, h1d*h2d*p4d*p5d*p6d, t3_d, t3_s_d);
  cudaDeviceSynchronize();
    //CHECK_ERR("Kernel execution failed");
    CUDA_SAFE(cudaMemcpy(((char *) energy_h) , ((char *) energy_d) , 
    size_energy*total_block*2, cudaMemcpyDeviceToHost));

    for(size_t i=1;i<dimGrid.x;i++)
      {
        energy_h[0]+=energy_h[i];
        energy_h[dimGrid.x]+=energy_h[i+dimGrid.x];
      }

     
//    printf("CUDA energy_h is %f %f %d %d %d %d %d %d\n", energy_h[0], energy_h[dimGrid.x]); //, total_size, h1d, h2d, p4d, p5d,p6d);
/*
    CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_d) , sizeof(double)*h3d*total_elements, cudaMemcpyDeviceToHost));
    CUDA_SAFE(cudaMemcpy(((char *) ts3) , ((char *) t3_s_d) , sizeof(double)*h3d*total_elements, cudaMemcpyDeviceToHost));
    total_s[0]=0.0, total_d[0]=0.0;
    for(size_t i=0;i<h3d*total_elements;i++) {
        total_s[0] += ts3[i];
        total_d[0] += t3[i];
    }
*/
//    printf("Total doubles and singles %f, %f\n", total_d, total_s);
    energy[0] = energy_h[0];
    energy[1] = energy_h[dimGrid.x];
    freeGpuMem(energy_d);
    freeGpuMem(eval_d1);
    freeGpuMem(eval_d2);
    freeGpuMem(eval_d3);
    freeGpuMem(eval_d4);
    freeGpuMem(eval_d5);
    freeGpuMem(eval_d6);
    freeHostMem(energy_h);
}
extern          "C" void
compute_en_(double * factor, double * energy, double * eval1,double* eval2,double* eval3,double* eval4,double* eval5,double* eval6, Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double* host1, double* host2)
//ckbn en_comment,double* total_d, double* total_s)
{
    compute_energy((double) *factor, energy, eval1,eval2, eval3, eval4, eval5, eval6,(size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d, host1, host2);
//ckbn en_comment    ,total_d, total_s);
}

//__device__ double* t3_d; 
extern    "C" void set_dev_mem_s(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d,size_t p6d)
{
    size_t size_t3;
    size_t3 = h1d*h2d*h3d*p4d*p5d*p6d;
    t3_s_d = (double *) getGpuMem(size_t3*sizeof(double));
    cudaMemset(t3_s_d,0,size_t3*sizeof(double));
}



extern          "C" void
dev_mem_s_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d)
{
    set_dev_mem_s((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d);
}

/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p5,p4] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_1_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3, double *t2_d, double *v2_d,size_t p4, size_t total_x, double* t3d) {
  size_t h1,h2,h3,p6;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p4*p4ld_t3]+=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2];
    }
  }
    __syncthreads();
}

extern          "C" void 
sd_t_s1_1_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
//CUDA_SAFE(cudaMalloc((void**) &t3_d, size_t3));
//CUDA_SAFE(cudaMalloc((void**) &t2_d, size_t2));
//CUDA_SAFE(cudaMalloc((void**) &v2_d, size_v2));
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h3ld_t3 = 1;
  h2ld_t3 = h3d;
  h1ld_t3 = h2d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_1_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t2_d,v2_d,i,total_x, t3_s_d);
    CHECK_ERR("Kernel execution failed");
  }
/*
    st = timer();
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }
*/
  cudaDeviceSynchronize();

//  CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));
  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
   //  cudaFree(t2_d);
   //  cudaFree(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_1_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_1_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d,  t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h1,h2,p6,p5,p4] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_2_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,double *t2_d, double *v2_d,size_t p4, size_t total_x, double* t3d) {
  size_t h1,h2,h3,p6;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p4*p4ld_t3]-=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2];
    }
  }
    __syncthreads();
}

extern          "C" void 
sd_t_s1_2_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d, double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
/*    if(first==1)
    {
    t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
  }*/
//CUDA_SAFE(cudaMalloc((void**) &t2_d, size_t2));
//CUDA_SAFE(cudaMalloc((void**) &v2_d, size_v2));
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
/*  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }*/
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d ;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h3ld_t3 = 1;
  h1ld_t3 = h3d;
  h2ld_t3 = h1d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
//  for(i=0;i<nstreams;++i){

    sd_t_s1_2_kernel<<<dimGrid,dimBlock,0>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t2_d,v2_d,i,total_x, t3_s_d);
    CHECK_ERR("Kernel execution failed");
//  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }*/
  cudaDeviceSynchronize();
//  CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));
/*
  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }*/
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
extern          "C" void 
sd_t_s1_2_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_2_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d,  t3, t2, v2);
}
extern          "C" void 
sd_t_s1_3_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d ;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h1ld_t3 = 1;
  h3ld_t3 = h1d;
  h2ld_t3 = h1d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_1_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t2_d,v2_d,i,total_x, t3_s_d);
    CHECK_ERR("Kernel execution failed");
  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }
*/  cudaDeviceSynchronize();
  //CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_3_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
  sd_t_s1_3_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d, t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_4_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p5ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t p4, size_t total_x) {
  size_t h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;
    rest_x=rest_x/p6d;
    p5=rest_x%p5d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p5*p5ld_t3+p4*p4ld_t3]-=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2+p5*p5ld_v2];
    }
  }
    __syncthreads();
}

extern          "C" void 
sd_t_s1_4_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
/*  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }*/
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
  p5ld_v2 = p6d * h3d * h2d;
  h3ld_t3 = 1;
  h2ld_t3 = h3d;
  h1ld_t3 = h2d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
  p4ld_t3 = p6d * h1d * h2d * h3d;
  p5ld_t3 = p4d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
   i=0;
 // for(i=0;i<nstreams;++i){
    sd_t_s1_4_kernel<<<dimGrid,dimBlock,0>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    //sd_t_s1_4_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
//  }


  cudaDeviceSynchronize();
  /*  CUDA_SAFE(cudaMemcpy(((char *) t3_p) , ((char *) t3_d) , size_block_t3, cudaMemcpyDeviceToHost));
  printf("Time for Async DeviceToHost %f\n", et-st);
  stream = 0;
//  while (stream < nstreams) {
//    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = t3_p; //[stream * size_el_block_t3];
    double         *dst = t3;  //[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] -= src[i];
    }
//    stream++;
//  }
*/
//  cudaDeviceSynchronize();
/*
  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }*/
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_4_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_4_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d,  t3, t2, v2);
}

/*----------------------------------------------------------------------*
 *t3[h3,h1,h2,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_5_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p5ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t p4, size_t total_x) {
  size_t h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;
    rest_x=rest_x/p6d;
    p5=rest_x%p5d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p5*p5ld_t3+p4*p4ld_t3]+=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2+p5*p5ld_v2];
    }
  }
    __syncthreads();
}

extern          "C" void 
sd_t_s1_5_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d ;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
  p5ld_v2 = p6d * h3d * h2d;
  h3ld_t3 = 1;
  h1ld_t3 = h3d;
  h2ld_t3 = h1d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
  p4ld_t3 = p6d * h1d * h2d * h3d;
  p5ld_t3 = p4d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_5_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }
*/
  cudaDeviceSynchronize();

  //CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));
  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_5_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
  sd_t_s1_5_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d,  t3, t2, v2);
}

/*----------------------------------------------------------------------*
 *t3[h1,h3,h2,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_6_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p5d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t p5ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p5ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t p4, size_t total_x) {
  size_t h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;
    rest_x=rest_x/p6d;
    p5=rest_x%p5d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p5*p5ld_t3+p4*p4ld_t3]-=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2+p5*p5ld_v2];
    }
  }
    __syncthreads();
}

extern          "C" void 
sd_t_s1_6_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double          *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
  p5ld_v2 = p6d * h3d * h2d;
  h1ld_t3 = 1;
  h3ld_t3 = h1d;
  h2ld_t3 = h1d * h3d;
  p6ld_t3 = h1d * h2d * h3d;
  p4ld_t3 = p6d * h1d * h2d * h3d;
  p5ld_t3 = p4d * p6d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_6_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
  }
/*  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }*/
  cudaDeviceSynchronize();
  //CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_6_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_6_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d, t3, t2, v2);
}









/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p4,p6,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_7_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t p4, size_t total_x) {
  size_t h1,h2,h3,p6;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p4*p4ld_t3]+=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2];
    }
  }
    __syncthreads();
}
extern          "C" void 
sd_t_s1_7_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double         *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h3ld_t3 = 1;
  h2ld_t3 = h3d;
  h1ld_t3 = h2d * h3d;
  p4ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p6ld_t3 = p4d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_7_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }

  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }*/
  cudaDeviceSynchronize();
  //CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_7_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_7_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d, t3, t2, v2);
}
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_8_kernel(size_t h1d,size_t h2d,size_t h3d,size_t p4d,size_t p6d,size_t p4ld_t2,size_t h1ld_t2,size_t h3ld_v2,size_t h2ld_v2,size_t p6ld_v2,size_t h3ld_t3,size_t h2ld_t3,size_t h1ld_t3,size_t p6ld_t3,size_t p4ld_t3,double *t3d, double *t2_d, double *v2_d,size_t p4, size_t total_x) {
  size_t h1,h2,h3,p6;
  __shared__ double t2_shm[T1*4*Tcomm];
  
  for(size_t i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  size_t rest_x=blockIdx.x;
  size_t thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(size_t i=0;i<total_x;i+=gridDim.x*blockDim.x)
  {
    rest_x += i;
    h3=rest_x%h3d;
    rest_x=rest_x/h3d;
    h2=rest_x%h2d;
    rest_x=rest_x/h2d;
    p6=rest_x%p6d;

    if((thread_x+i)<total_x)
    for(h1=0;h1<h1d;h1++)
    for(p4=0;p4<p4d;p4++)
    {
      t3d[h3*h3ld_t3+h2*h2ld_t3+h1*h1ld_t3+p6*p6ld_t3+p4*p4ld_t3]-=t2_shm[h1*p4d+p4]*v2_d[h3*h3ld_v2+h2*h2ld_v2+p6*p6ld_v2];
    }
  }
    __syncthreads();
}
/*----------------------------------------------------------------------*
 *t3[h3,h1,h2,p4,p6,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
extern          "C" void 
sd_t_s1_8_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double          *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h3ld_t3 = 1;
  h1ld_t3 = h3d;
  h2ld_t3 = h1d * h3d;
  p4ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p6ld_t3 = p4d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_8_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }
  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }*/
  cudaDeviceSynchronize();
//  CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
//  freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
extern          "C" void 
sd_t_s1_8_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
  sd_t_s1_8_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d, t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h1,h3,h2,p4,p6,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
extern          "C" void 
sd_t_s1_9_cuda(size_t h1d, size_t h2d, size_t h3d, size_t p4d, size_t p5d, size_t p6d,  double *t3, double *t2, double *v2)
{
    double st, et;
//ckbn    st = timer(); 
  size_t          p7ld_t2, p4ld_t2, h1ld_t2, h2ld_v2, p7ld_v2, h3ld_v2,
                  p6ld_v2, p5ld_v2, h3ld_t3, h2ld_t3, h1ld_t3, p6ld_t3,
                  p5ld_t3, p4ld_t3;
  size_t          size_t3, size_block_t3, size_el_block_t3, size_t2,
                  size_v2;
  cudaStream_t   *streams;
  size_t          nstreams, i;
  double          *t2_d, *v2_d, *t3_p;
  //size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
  size_t2 = p4d * h1d * sizeof(double);
  size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
  nstreams = 1;
  size_block_t3 = size_t3 / nstreams;
  size_el_block_t3 = size_block_t3 / sizeof(double);
  /*  if(first==1)
    {
        t3_d = (double *) getGpuMem(size_t3);
        cudaMemset(t3_d,0,size_t3*sizeof(double));
        first = 0;
    }
*/
//  t3_d = (double *) getGpuMem(size_t3);
  t2_d = (double *) getGpuMem(size_t2);
  v2_d = (double *) getGpuMem(size_v2);
  //t3_p = (double *) getHostMem(size_t3);
  streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  assert(streams != NULL);
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaStreamCreate(&streams[i]));
  }
  CUDA_SAFE(cudaMemcpy(t2_d, t2, size_t2, cudaMemcpyHostToDevice));
  CUDA_SAFE(cudaMemcpy(v2_d, v2, size_v2, cudaMemcpyHostToDevice));

  p4ld_t2 = 1;
  h1ld_t2 = p4d;

  h3ld_v2 = 1;
  h2ld_v2 = h3d;
  p6ld_v2 = h3d * h2d;
//  p5ld_v2 = p6d * h3d * p7d;
  h1ld_t3 = 1;
  h3ld_t3 = h1d;
  h2ld_t3 = h1d * h3d;
  p4ld_t3 = h1d * h2d * h3d;
//  p5ld_t3 = p6d * h1d * h2d * h3d;
  p6ld_t3 = p4d * h1d * h2d * h3d;
  size_t total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_7_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    CHECK_ERR("Kernel execution failed");
  }
/*
  for (i = 0; i < nstreams; ++i) {
    CUDA_SAFE(cudaMemcpyAsync(((char *) t3_p) + i * size_block_t3, ((char *) t3_s_d) + i * size_block_t3, size_block_t3, cudaMemcpyDeviceToHost, streams[i]));
  }
  stream = 0;
  while (stream < nstreams) {
    while (cudaStreamQuery(streams[stream]) != cudaSuccess);
    double         *src = &t3_p[stream * size_el_block_t3];
    double         *dst = &t3[stream * size_el_block_t3];
    for (i = 0; i < size_el_block_t3; ++i) {
      dst[i] = src[i];
    }
    stream++;
  }*/
  cudaDeviceSynchronize();
  //CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

//  printf("out is %lf\n", t3_p[0]);
  for (i = 0; i < nstreams; ++i) {
    cudaStreamDestroy(streams[i]);
  }
  //freeGpuMem(t3_d);
  freeGpuMem(t2_d);
  freeGpuMem(v2_d);
  //freeHostMem(t3_p);
  free(streams);
}
extern          "C" void 
sd_t_s1_9_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
  sd_t_s1_9_cuda((size_t) *h1d, (size_t) *h2d, (size_t) *h3d, (size_t) *p4d, (size_t) *p5d, (size_t) *p6d,  t3, t2, v2);
}
