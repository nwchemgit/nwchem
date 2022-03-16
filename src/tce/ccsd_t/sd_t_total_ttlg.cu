#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda.h>

__device__ double* t3_s_d;
__device__ double* t3_d;

#include "header.h"
#include "ourinclude.h"

#define T1 16
#define T2 16
#define Tcomm 16

cublasHandle_t handle;
double* output_d;
size_t current_i_size;
extern "C" void ttgt_init()
{
  cublasCreate(&handle);
  output_d = NULL;
  current_i_size = 0;
}
extern    "C" void set_dev_mem_d(int h1d, int h2d, int h3d, int p4d, int p5d,int p6d)
{
    int size_t3;
    size_t3 = h1d*h2d*h3d*p4d*p5d*p6d;
    t3_d = (double *) getGpuMem(size_t3*sizeof(double));
    cudaMemset(t3_d,0,size_t3*sizeof(double));
}
extern          "C" void
dev_mem_d_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d)
{
    set_dev_mem_d((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d);
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
extern "C" void sd_t_d1_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub, int id) {
double* output_d;
	static int count = 0;
	if(count == 0)
	{
		ttgt_init();
		count++;
	}

	size_t size_triplesx,size_block_triplesx,size_el_block_triplesx,size_t2sub,size_v2sub;
        size_t i;
        double *t2sub_d,*v2sub_d;
        size_triplesx= p4d * p5d * h1d * h3d * h2d * p6d *sizeof(double);
        size_t2sub=h7d*p4d*p5d*h1d*sizeof(double);
        size_v2sub=h3d*h2d*p6d*h7d*sizeof(double);
        int i1[4], i2[4], o[6];
        i1[0] = h7d;
        i1[1] = p4d;
        i1[2] = p5d;
        i1[3] = h1d;
        i2[0] = h3d;
        i2[1] = h2d;
        i2[2] = p6d;
        i2[3] = h7d;
        o[0] = p4d;
        o[1] = p5d;
        o[2] = h1d;
        o[3] = h3d;
        o[4] = h2d;
        o[5] = p6d;
        cublasOperation_t transa, transb;
        transa = CUBLAS_OP_T;
        transb = CUBLAS_OP_T;
        size_t m,n,k;
        m = p4d*p5d*h1d;
        k = h7d;
        n = h3d*h2d*p6d;
        double alpha, beta;
        alpha = 1;
        beta = 0;
        t2sub_d=(double*)getGpuMem(size_t2sub);
        v2sub_d=(double*)getGpuMem(size_v2sub);
	//if(size_triplesx > current_i_size)
	{
        	output_d=(double*)getGpuMem(size_triplesx);
		current_i_size = size_triplesx;
	}
	if(output_d == NULL) 
	{
		exit(0);
	}
        int perm[6];
	//double beta;
	switch(id)
	{
		case 1:
        		perm[0] = 3;
		        perm[1] = 4;
		        perm[2] = 2;
		        perm[3] = 5;
		        perm[4] = 1;
		        perm[5] = 0;
			beta = -1.0;
		break;
		case 2:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 4;
		        perm[3] = 5;
		        perm[4] = 1;
		        perm[5] = 0;
			beta = 1.0;
		break;
		case 3:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 4;
		        perm[3] = 5;
		        perm[4] = 1;
		        perm[5] = 0;
			beta = -1.0;
		break;
		case 4:
        		perm[0] = 3;
		        perm[1] = 4;
		        perm[2] = 2;
		        perm[3] = 1;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = -1.0;
		break;
		case 5:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 4;
		        perm[3] = 1;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = 1.0;
		break;
		case 6:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 4;
		        perm[3] = 1;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = -1.0;
		break;
		case 7:
        		perm[0] = 3;
		        perm[1] = 4;
		        perm[2] = 2;
		        perm[3] = 1;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = 1.0;
		break;
		case 8:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 4;
		        perm[3] = 1;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = -1.0;
		break;
		case 9:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 4;
		        perm[3] = 1;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = 1.0;
		break;
	}

        cublasDgemm(handle, transa, transb, m, n, k, &alpha, t2sub_d, h7d, v2sub_d, n, &beta, output_d, m);
       ttlg_transpose(6, o, perm, output_d, t3_d, 1, beta);
        cudaDeviceSynchronize();
        freeGpuMem(t2sub_d);
        freeGpuMem(v2sub_d);
  	freeGpuMem(output_d);
}

/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p6,p5,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_1_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 1);
}

extern "C" void sd_t_d1_1_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_1_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,h2,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_2_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 2);
}
extern "C" void sd_t_d1_2_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_2_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_3_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 3);
}
extern "C" void sd_t_d1_3_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_3_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p5,p4,p6] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_4_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 4);
}
extern "C" void sd_t_d1_4_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_4_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p4,p6] += t2sub[h7,p4,p5,h1] * v2sub[h3,h2,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_5_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 5);
}
extern "C" void sd_t_d1_5_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_5_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p4,p6] -= t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_6_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 6);
}
extern "C" void sd_t_d1_6_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_6_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,p5,p6,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_7_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 7);
}
extern "C" void sd_t_d1_7_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_7_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h3,h1,h2,p5,p6,p4] -= t2sub[h7,p4,p5,h1] * v2sub[h3,h2,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_8_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 8);
}
extern "C" void sd_t_d1_8_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_8_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}
/*----------------------------------------------------------------------*
 *triplesx[h1,h3,p5,p6,p4] += t2sub[h7,p4,p5,h1] * v2sub[h3,p6,h7]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d1_9_cuda(int h1d, int h2d, int h3d, int h7d, int p4d, int p5d, int p6d, double *triplesx, double *t2sub, double *v2sub) {
sd_t_d1_cuda( h1d,  h2d,  h3d,  h7d,  p4d,  p5d,  p6d, triplesx, t2sub, v2sub, 9);
}
extern "C" void sd_t_d1_9_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* h7d, Integer* p4d, Integer* p5d, Integer* p6d, double *triplesx, double *t2sub, double *v2sub) {
  sd_t_d1_9_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*h7d,(int)*p4d,(int)*p5d,(int)*p6d,triplesx,t2sub,v2sub);
}


extern "C" void sd_t_d2_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *triplesx, double *t2sub, double *v2sub, int id) {
 size_t size_triplesx,size_block_triplesx,size_el_block_triplesx,size_t2sub,size_v2sub;
double* output_d;
        size_t i;
        double *t2sub_d,*v2sub_d;
        size_triplesx= p4d * p5d * h1d * h3d * h2d * p6d *sizeof(double);
        size_t2sub=p7d*p4d*p5d*h1d*sizeof(double);
        size_v2sub=h3d*h2d*p6d*p7d*sizeof(double);
        int i1[4], i2[4], o[6];
        i1[0] = p7d;
        i1[1] = p4d;
        i1[2] = h1d;
        i1[3] = h2d;
        i2[0] = p7d;
        i2[1] = h3d;
        i2[2] = p6d;
        i2[3] = p5d;
        o[0] = p4d;
        o[1] = h1d;
        o[2] = h2d;
        o[3] = h3d;
        o[4] = p6d;
        o[5] = p5d;
        cublasOperation_t transa, transb;
        transa = CUBLAS_OP_T;
        transb = CUBLAS_OP_N;
        size_t m,n,k;
        m = p4d*h1d*h2d;
        k = p7d;
        n = h3d*p6d*p5d;
        double alpha, beta;
        alpha = 1;
        beta = 0;
        t2sub_d=(double*)getGpuMem(size_t2sub);
        v2sub_d=(double*)getGpuMem(size_v2sub);
	//if(size_triplesx > current_i_size)
	{
        	output_d=(double*)getGpuMem(size_triplesx);
		current_i_size = size_triplesx;
	}
	if(output_d == NULL) 
	{
		exit(0);
	}
        int perm[6];
	//double beta;
	switch(id)
	{
		case 1:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 1;
		        perm[3] = 4;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = -1.0;
		break;
		case 2:
        		perm[0] = 2;
		        perm[1] = 1;
		        perm[2] = 3;
		        perm[3] = 4;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = -1.0;
		break;
		case 3:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 1;
		        perm[3] = 4;
		        perm[4] = 5;
		        perm[5] = 0;
			beta = 1.0;
		break;
		case 4:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 1;
		        perm[3] = 4;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = 1.0;
		break;
		case 5:
        		perm[0] = 2;
		        perm[1] = 1;
		        perm[2] = 3;
		        perm[3] = 4;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = 1.0;
		break;
		case 6:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 1;
		        perm[3] = 4;
		        perm[4] = 0;
		        perm[5] = 5;
			beta = -1.0;
		break;
		case 7:
        		perm[0] = 3;
		        perm[1] = 2;
		        perm[2] = 1;
		        perm[3] = 0;
		        perm[4] = 4;
		        perm[5] = 5;
			beta = -1.0;
		break;
		case 8:
        		perm[0] = 2;
		        perm[1] = 1;
		        perm[2] = 3;
		        perm[3] = 0;
		        perm[4] = 4;
		        perm[5] = 5;
			beta = -1.0;
		break;
		case 9:
        		perm[0] = 2;
		        perm[1] = 3;
		        perm[2] = 1;
		        perm[3] = 0;
		        perm[4] = 4;
		        perm[5] = 5;
			beta = 1.0;
		break;
	}

        cublasDgemm(handle, transa, transb, m, n, k, &alpha, t2sub_d, p7d, v2sub_d, n, &beta, output_d, m);
     ttlg_transpose(6, o, perm, output_d, t3_d, 1, beta);
        cudaDeviceSynchronize();
        freeGpuMem(t2sub_d);
        freeGpuMem(v2sub_d);
        	freeGpuMem(output_d);
}

/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_1_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 1);
}
extern "C" void sd_t_d2_1_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_1_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4] -= t2[p7,p4,h1,h2] * v2[p7,h3]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_2_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 2);
}
extern "C" void sd_t_d2_2_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_2_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p6,p4] += t2[p7,p4,h1,h2] * v2[p7,h3,p6]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_3_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 3);
}
extern "C" void sd_t_d2_3_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_3_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/

extern "C" void sd_t_d2_4_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 4);
}
extern "C" void sd_t_d2_4_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_4_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p5]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_5_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 5);
}
extern "C" void sd_t_d2_5_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_5_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p6,p4,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_6_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 6);
}
extern "C" void sd_t_d2_6_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_6_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p4,p6,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_7_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 7);
}
extern "C" void sd_t_d2_7_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_7_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h1,h3,p4,p6,p5] -= t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_8_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 8);
}
extern "C" void sd_t_d2_8_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_8_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}
/*----------------------------------------------------------------------*
 *t3[h2,h3,h1,p4,p6,p5] += t2[p7,p4,h1,h2] * v2[p7,h3,p6,p5]
 *----------------------------------------------------------------------*/
extern "C" void sd_t_d2_9_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, int p7d, double *t3, double *t2, double *v2) {
	sd_t_d2_cuda(h1d, h2d, h3d, p4d,  p5d, p6d,  p7d, t3, t2, v2, 9);
}

extern "C" void sd_t_d2_9_cuda_(Integer *h1d, Integer* h2d, Integer* h3d, Integer* p4d, Integer* p5d, Integer* p6d, Integer* p7d, double *t3, double *t2, double *v2) {
  sd_t_d2_9_cuda((int)*h1d,(int)*h2d,(int)*h3d,(int)*p4d,(int)*p5d,(int)*p6d,(int)*p7d,t3,t2,v2);
}


#define MAX_h3 64
/* IMPORTANT!!!!
t3_d must be passed as parameter to kernel function. A __global__ function can't access the global variable directly*/

__global__ void compute_energy_kernel(int h1d,int h2d,int h3d,int p4d,int p5d,int p6d,double* eval1,double* eval2,double* eval3,double* eval4,double* eval5,double* eval6, double* energy, double factor, int total_size, double* t3d, double* t3_sd)
{
  int h1,h2,p6,p4,p5, h3,i=0;
  double e1,e2,e4,e5,e6;
//  __shared__ double t2_shm[MAX_h3];
  __shared__ double energy_s[T1];
  __shared__ double energy2_s[T1];
  double inner_fac;
  int limit;
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  if(threadIdx.x==0)
  {
        energy[blockIdx.x]=0;
        energy[blockIdx.x+gridDim.x]=0;
        energy_s[threadIdx.x] = 0.0;
        energy2_s[threadIdx.x] = 0.0;
  }

  for(int j =0; j<T2*T1;j++) {
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
    for(int i=0;i<h3d;i++)
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
/*	  limit = blockDim.x;
      if (blockIdx.x == (gridDim.x-1)) limit = total_size%blockDim.x;
      for(int i=0;i<limit;i++)
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

extern   "C" void compute_energy(double factor, double* energy, double* eval1, double* eval2,double* eval3,double* eval4,double* eval5,double* eval6,int h1d, int h2d, int h3d, int p4d, int p5d,int p6d, double* host1, double* host2)
//ckbn en_comment, double* total_d, double* total_s)
{
    double* energy_d, *energy_h;
    double* eval_d1,*eval_d2,*eval_d3,*eval_d4,*eval_d5,*eval_d6;
    int size_energy = 2*sizeof(double);
    int total_block = DIV_UB((h1d*h2d*p4d*p5d*p6d), (T2*T1));

//    int total_block = 1;
    int total_elements = h1d*h2d*p4d*p5d*p6d;

    energy_d = (double*)getGpuMem(size_energy*total_block*2);
    int i=0,in; 
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

    for(int i=1;i<dimGrid.x;i++)
      {
        energy_h[0]+=energy_h[i];
        energy_h[dimGrid.x]+=energy_h[i+dimGrid.x];
      }

     
//    printf("CUDA energy_h is %f %f %d %d %d %d %d %d\n", energy_h[0], energy_h[dimGrid.x]); //, total_size, h1d, h2d, p4d, p5d,p6d);
/*
    CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_d) , sizeof(double)*h3d*total_elements, cudaMemcpyDeviceToHost));
    CUDA_SAFE(cudaMemcpy(((char *) ts3) , ((char *) t3_s_d) , sizeof(double)*h3d*total_elements, cudaMemcpyDeviceToHost));
    total_s[0]=0.0, total_d[0]=0.0;
    for(int i=0;i<h3d*total_elements;i++) {
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
    compute_energy((double) *factor, energy, eval1,eval2, eval3, eval4, eval5, eval6,(int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d, host1, host2);
//ckbn en_comment    ,total_d, total_s);
}

//__device__ double* t3_d; 
extern    "C" void set_dev_mem_s(int h1d, int h2d, int h3d, int p4d, int p5d,int p6d)
{
    int size_t3;
    size_t3 = h1d*h2d*h3d*p4d*p5d*p6d;
    t3_s_d = (double *) getGpuMem(size_t3*sizeof(double));
    cudaMemset(t3_s_d,0,size_t3*sizeof(double));
}



extern          "C" void
dev_mem_s_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d)
{
    set_dev_mem_s((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d);
}

/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p5,p4] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_1_kernel(int h1d,int h2d,int h3d,int p4d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p4ld_t3, double *t2_d, double *v2_d,int p4, int total_x, double* t3d) {
  int h1,h2,h3,p6;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_1_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
	size_t2 = p4d * h1d * sizeof(double);
	size_v2 = h3d * h2d * p6d * p5d * sizeof(double);
	nstreams = 1;
	size_block_t3 = size_t3 / nstreams;
	size_el_block_t3 = size_block_t3 / sizeof(double);
//CUDA_SAFE(cudaMalloc((void**) &t3_d, size_t3));
//CUDA_SAFE(cudaMalloc((void**) &t2_d, size_t2));
//CUDA_SAFE(cudaMalloc((void**) &v2_d, size_v2));
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h3ld_t3 = 1;
	h2ld_t3 = h3d;
	h1ld_t3 = h2d * h3d;
	p6ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
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

//	CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));
	for (i = 0; i < nstreams; ++i) {
		cudaStreamDestroy(streams[i]);
	}
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
   //  cudaFree(t2_d);
   //  cudaFree(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_1_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_1_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d,  t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h1,h2,p6,p5,p4] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_2_kernel(int h1d,int h2d,int h3d,int p4d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p4ld_t3,double *t2_d, double *v2_d,int p4, int total_x, double* t3d) {
  int h1,h2,h3,p6;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_2_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d, double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
	t3_p = (double *) getHostMem(size_t3);
	streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
/*	assert(streams != NULL);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h3ld_t3 = 1;
	h1ld_t3 = h3d;
	h2ld_t3 = h1d * h3d;
	p6ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
//  for(i=0;i<nstreams;++i){

    sd_t_s1_2_kernel<<<dimGrid,dimBlock,0>>>(h1d,h2d,h3d,p4d,p5d*p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p4ld_t3,t2_d,v2_d,i,total_x, t3_s_d);
		CHECK_ERR("Kernel execution failed");
//	}
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
//	CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));
/*
	for (i = 0; i < nstreams; ++i) {
		cudaStreamDestroy(streams[i]);
	}*/
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
extern          "C" void 
sd_t_s1_2_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_2_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d,  t3, t2, v2);
}
extern          "C" void 
sd_t_s1_3_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
	t3_p = (double *) getHostMem(size_t3);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h1ld_t3 = 1;
	h3ld_t3 = h1d;
	h2ld_t3 = h1d * h3d;
	p6ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p4ld_t3 = p5d * p6d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
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
*/	cudaDeviceSynchronize();
	//CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

	for (i = 0; i < nstreams; ++i) {
		cudaStreamDestroy(streams[i]);
	}
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_3_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
	sd_t_s1_3_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d, t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_4_kernel(int h1d,int h2d,int h3d,int p4d,int p5d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int p5ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p5ld_t3,int p4ld_t3,double *t3d, double *t2_d, double *v2_d,int p4, int total_x) {
  int h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_4_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
	streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
/*	assert(streams != NULL);
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
  int total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
   i=0;
 // for(i=0;i<nstreams;++i){
    sd_t_s1_4_kernel<<<dimGrid,dimBlock,0>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
    //sd_t_s1_4_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_d,t2_d,v2_d,i,total_x);
		CHECK_ERR("Kernel execution failed");
//	}


	cudaDeviceSynchronize();
	/*	CUDA_SAFE(cudaMemcpy(((char *) t3_p) , ((char *) t3_d) , size_block_t3, cudaMemcpyDeviceToHost));
	printf("Time for Async DeviceToHost %f\n", et-st);
	stream = 0;
//	while (stream < nstreams) {
//		while (cudaStreamQuery(streams[stream]) != cudaSuccess);
		double         *src = t3_p; //[stream * size_el_block_t3];
		double         *dst = t3;  //[stream * size_el_block_t3];
		for (i = 0; i < size_el_block_t3; ++i) {
			dst[i] -= src[i];
		}
//		stream++;
//	}
*/
//	cudaDeviceSynchronize();
/*
	for (i = 0; i < nstreams; ++i) {
		cudaStreamDestroy(streams[i]);
	}*/
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_4_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_4_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d,  t3, t2, v2);
}

/*----------------------------------------------------------------------*
 *t3[h3,h1,h2,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_5_kernel(int h1d,int h2d,int h3d,int p4d,int p5d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int p5ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p5ld_t3,int p4ld_t3,double *t3d, double *t2_d, double *v2_d,int p4, int total_x) {
  int h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_5_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
  int total_x = h3d*h2d*p6d*p5d;
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
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_5_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
	sd_t_s1_5_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d,  t3, t2, v2);
}

/*----------------------------------------------------------------------*
 *t3[h1,h3,h2,p6,p4,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_6_kernel(int h1d,int h2d,int h3d,int p4d,int p5d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int p5ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p5ld_t3,int p4ld_t3,double *t3d, double *t2_d, double *v2_d,int p4, int total_x) {
  int h1,h2,h3,p6,p5;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_6_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
  int total_x = h3d*h2d*p6d*p5d;
  dim3 dimBlock(T2*T1);dim3 dimGrid(DIV_UB(total_x,T2*T1), 1);
  for(i=0;i<nstreams;++i){
    sd_t_s1_6_kernel<<<dimGrid,dimBlock,0,streams[i]>>>(h1d,h2d,h3d,p4d,p5d,p6d,p4ld_t2,h1ld_t2,h3ld_v2,h2ld_v2,p6ld_v2,p5ld_v2,h3ld_t3,h2ld_t3,h1ld_t3,p6ld_t3,p5ld_t3,p4ld_t3,t3_s_d,t2_d,v2_d,i,total_x);
		CHECK_ERR("Kernel execution failed");
	}
/*	for (i = 0; i < nstreams; ++i) {
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
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_6_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_6_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d, t3, t2, v2);
}









/*----------------------------------------------------------------------*
 *t3[h3,h2,h1,p4,p6,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_7_kernel(int h1d,int h2d,int h3d,int p4d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p4ld_t3,double *t3d, double *t2_d, double *v2_d,int p4, int total_x) {
  int h1,h2,h3,p6;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_7_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h3ld_t3 = 1;
	h2ld_t3 = h3d;
	h1ld_t3 = h2d * h3d;
	p4ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p6ld_t3 = p4d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
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
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
#undef T1
#undef T2
#undef Tcomm
extern          "C" void 
sd_t_s1_7_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_7_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d, t3, t2, v2);
}
#define T1 16
#define T2 16
#define Tcomm 16
__global__ void sd_t_s1_8_kernel(int h1d,int h2d,int h3d,int p4d,int p6d,int p4ld_t2,int h1ld_t2,int h3ld_v2,int h2ld_v2,int p6ld_v2,int h3ld_t3,int h2ld_t3,int h1ld_t3,int p6ld_t3,int p4ld_t3,double *t3d, double *t2_d, double *v2_d,int p4, int total_x) {
  int h1,h2,h3,p6;
  __shared__ double t2_shm[T1*2*Tcomm];
  
  for(int i=threadIdx.x;i<h1d*p4d;i+=blockDim.x)
  if(i<h1d*p4d)
  t2_shm[i] = t2_d[i];
  int rest_x=blockIdx.x;
  int thread_x = T2*T1 * rest_x + threadIdx.x;
  rest_x = thread_x;
    __syncthreads();
/* the following computation may need to happen inside the loop */
  for(int i=0;i<total_x;i+=gridDim.x*blockDim.x)
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
sd_t_s1_8_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h3ld_t3 = 1;
	h1ld_t3 = h3d;
	h2ld_t3 = h1d * h3d;
	p4ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p6ld_t3 = p4d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
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
//	CUDA_SAFE(cudaMemcpy(((char *) t3) , ((char *) t3_s_d) , size_t3, cudaMemcpyDeviceToHost));

	for (i = 0; i < nstreams; ++i) {
		cudaStreamDestroy(streams[i]);
	}
//	freeGpuMem(t3_d);
	freeGpuMem(t2_d);
	freeGpuMem(v2_d);
	freeHostMem(t3_p);
	free(streams);
}
extern          "C" void 
sd_t_s1_8_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d, double *t3, double *t2, double *v2)
{
	sd_t_s1_8_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d, t3, t2, v2);
}
/*----------------------------------------------------------------------*
 *t3[h1,h3,h2,p4,p6,p5] -= t2[p4,h1] * v2[h3,h2,p6,p5]
 *----------------------------------------------------------------------*/
extern          "C" void 
sd_t_s1_9_cuda(int h1d, int h2d, int h3d, int p4d, int p5d, int p6d,  double *t3, double *t2, double *v2)
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
	size_t3 = h3d * h2d * h1d * p6d * p5d * p4d * sizeof(double);
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
//	t3_d = (double *) getGpuMem(size_t3);
	t2_d = (double *) getGpuMem(size_t2);
	v2_d = (double *) getGpuMem(size_v2);
	t3_p = (double *) getHostMem(size_t3);
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
//	p5ld_v2 = p6d * h3d * p7d;
	h1ld_t3 = 1;
	h3ld_t3 = h1d;
	h2ld_t3 = h1d * h3d;
	p4ld_t3 = h1d * h2d * h3d;
//	p5ld_t3 = p6d * h1d * h2d * h3d;
	p6ld_t3 = p4d * h1d * h2d * h3d;
  int total_x = h3d*h2d*p6d*p5d;
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
	freeHostMem(t3_p);
	free(streams);
}
extern          "C" void 
sd_t_s1_9_cuda_(Integer * h1d, Integer * h2d, Integer * h3d, Integer * p4d, Integer * p5d, Integer * p6d,  double *t3, double *t2, double *v2)
{
	sd_t_s1_9_cuda((int) *h1d, (int) *h2d, (int) *h3d, (int) *p4d, (int) *p5d, (int) *p6d,  t3, t2, v2);
}
