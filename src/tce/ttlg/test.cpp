#include<iostream>
#include "ourproj.h"
#include "TransposeSpec.h"
#include "ParameterTuner.cpp"
#include <cstring>
#include <omp.h>
#include <fstream>
#include <time.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <cuComplex.h>
#include <complex.h>

#include <stdlib.h>
#include <cuda_runtime.h>
#include <complex.h>

#include "ourinclude.h"

#ifdef type
#undef type
#endif

#ifndef type
#define type double
#endif

#define MIN(a,b) (a < b?a:b)
int main(int argc, char* argv[])
{
	FILE * f;
	int dim[10];// = {15, 20, 20, 19, 16, 16};
	int permutation[10];// = {0, 5, 4, 1, 2, 3};
	int i;
	//#define type float
#define type double
	if(argc < 2)
	{
		fprintf(stderr, "Please input filename to read data\n");
		exit(0);
	}
	f = fopen(argv[1], "r");
	char line[255];
	char temp[50];
	cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
	int icase=-1;
	int ndim;
	while(1)
	{
		bool finish = false;
		size_t total = 1;
		fscanf(f, "%d", &ndim);
		if(feof(f)){
			break;}
		if(ndim == -1)
		{
			char line[255];
			fgets(line, 255, f);
			continue;
		}
		for(i = 0; i < ndim; i++)
		{
			fscanf(f, "%d", &dim[i]);
			if(feof(f)){
				finish = true;
				break;}
			printf("%d ", dim[i]);
			total *= dim[i];
		}
		printf("\t");
		for(i = 0; i < ndim; i++)
		{
			fscanf(f, "%d", &permutation[i]);
			printf("%d ", permutation[i]);
		}

		size_t totalsize = total * sizeof(type);

		type* A = (type*) malloc(totalsize);
		if(!A) {cout <<"Memory error "; exit(0);}
		type* B = (type*) malloc(totalsize);
		if(!B) {cout <<"Memory error "; exit(0);}
		type*	B_ref = (type *) malloc(totalsize);
		if(!B_ref) {cout <<"Memory error "; exit(0);}
		type *d_A, *d_B;
		for(i=0; i < total ; ++i){
			A[i] = (type)i;
			//B[i] = (type)-i;
		}

		cudaMalloc(&d_A, totalsize);
		//cout <<"A addr: " <<d_A;
		cudaMemcpy(d_A, A, totalsize, cudaMemcpyHostToDevice);
		cudaMalloc(&d_B, totalsize);
		//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
		//cout <<"\nB addr: " <<d_B;
		cudaMemcpy(d_B, B, totalsize, cudaMemcpyHostToDevice);

		type alpha=1, beta=0;
		TensorType  mytype = doubletype;
		//TensorType  mytype = floattype;
		double start, tmpTime;
		start = omp_get_wtime();

		ttlg_transpose(ndim,dim, permutation, d_A, d_B ,alpha, beta);
		cudaDeviceSynchronize();
		tmpTime = omp_get_wtime() - start;
		//cout << "\t"<<tmpTime<<"\t";
		fprintf(stdout, "\t%6.6lf\t", tmpTime);
		//
		//                              //dCuTranspose_03241_15x20x20x16x304_reference(A, B_ref, sizes1, sizes1, ldb);
		transpose_check(ndim,A, B_ref, alpha, beta, dim, permutation);
		//tranpose5_check(A, B_ref, sizes1, sizes1, sizes1, perm1);
		cudaMemcpy(B, d_B,totalsize, cudaMemcpyDeviceToHost);

#ifdef printd
		cout<<"\n";
		for(int i = 0; i < MIN(1000, total); i++)
			cout <<B[i]<<" ";
		cout<<"\n";
		for(int i = 0; i < MIN(1000, total); i++)
			cout <<A[i]<<" ";
		cout<<"\n"; 
#endif
#ifndef NOERRORC
		transpose_equal(B,B_ref,total);
		//	fprintf(stdout, "\t%f",tmpTime);
#endif
		fprintf(stdout, "\t%6.2lf\t %u\n",2*totalsize/(tmpTime*1000000000), totalsize);
		free(A);
		free(B);
		free(B_ref);
		cudaFree(d_A);
		cudaFree(d_B);

	}
#undef type
}
