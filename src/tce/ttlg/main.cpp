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

int getoff(int index,const int * dims,const int * stride, int n);


void transpose_check(int n, const type *A, type *B, const type alpha, const type beta,  const int *lda,  int *perm)
{
	int i = 0, j;
	int r_perm[100], lda_s[100], ldb_s[100], temp[100] ;
	lda_s[0] = 1;
	ldb_s[0] = 1;
	long size = 1;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			if(perm[j] == i)
			{
				r_perm[i] = j; 
				break;
			}
		}
		size *= lda[i];
	}	
	for(i = 1; i < n; i++)
	{
		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * lda[perm[i-1]];
	}
       for(i = 0; i < n; i++)
        {
                temp[i] = ldb_s[r_perm[i]];
        }
	
	for(i=0; i < size; i++)
	{
			B[getoff(i, lda, temp, n)] = A[i];
	}
}



void fuseIndices(TransposeSpec &spec);
using namespace std;

//template <typename tensortype>
//void ourproj_tranpose(int ndim, int *dims, int *perm, tensortype *input, tensortype *output, TensorType dataType = 1, tensortype alpha = 0.0, tensortype beta = 1.0){
//template <typename tensortype>
extern "C"
void ttlg_transpose(int ndim, int *dims, int *perm, double *input, double *output,double alpha = 1.0, double beta = 0.0){
#ifdef INPUT	
	cout <<"Dims: ";
		for(int i = 0; i < ndim; i++)
		{
			cout << dims[i] <<" ";
		}
	cout <<"\nPerm: ";
		for(int i = 0; i < ndim; i++)
		{
			cout << perm[i] <<" ";
		}
	cout <<"\n";
#endif

	TensorType dataType = (TensorType) 1;

	TransposeSpec *spec = new TransposeSpec(ndim, dims, perm,dataType );

	fuseIndices(*spec); //sets the fused indices to the spec
	spec->setDataType(dataType);
	spec->setBeta(beta);
	spec->setAlpha(alpha);


	ParameterTuner *parameterTuner = new ParameterTuner();
	Parameters & myParams = parameterTuner->tune(*spec);

	int * sizes = spec->getSizes();
	int n = spec->getNdim();
	if(n == 1 || n == 0)//simply copy out
	{
#ifdef printd
		cout <<"\n1 D\n";
#endif
		int ndim = 1;
		int lda[1];
		lda[0] = sizes[0];
#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif

		cudaMemcpy(output, input, sizes[0] * sizeof(input[0]), cudaMemcpyDeviceToDevice);
#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif

	}
	else
	{
#ifdef BW
cout << myParams . getBW()<<"\t"; 
#endif
#ifdef EFF
cout << myParams . getWarpEfficiency()<<"\t"; 
#endif
#ifdef TIME
cout << myParams. getTime(); 
#endif

		int params[12];
		int *perm_new = spec->getPermutation();
#ifdef printd
		cout <<endl<< myParams. getSharedMemSize1()<<"\t";
		cout<< myParams.getTileSize()<<"\t";
		cout<< myParams.getPaddingSize()<<"\t";
		cout<< myParams.getTbSize()<<endl;
		cout<< myParams.getNumElementsProcessedPerBlock()<<endl;
#endif
		params[0] = myParams.getTileSize();
		params[1] = myParams.getPaddingSize();
		params[2] = myParams.getTbSize();
		params[5] = myParams. getSharedMemSize1();
		int ldb[20];
		int rperm[20];
#ifdef printd
		cout <<"New dims: ";
#endif
		for(int i = 0; i < n; i++)
		{
#ifdef printd
			cout << sizes[i] <<" ";
#endif
			ldb[i] = sizes[perm_new[i]];
			for(int j = 0; j < n; j++)
			{
				if(perm_new[i] == j)
				{
					rperm[j] = i;
				}
			}
		}


#ifdef printd
		cout << "\n";
#endif
#ifdef printd
		cout << "Perm: ";
		for(int j = 0; j < n; j++)
		{
			cout <<perm_new[j]<<" ";
		}
		cout << "\n";
#endif
		params[3]=perm_new[1];
		params[4]= rperm[1];
		
		int tmp;
		BlockingCase caseId = parameterTuner->getCaseId();
		switch(caseId.getMode())
		{
			case BlockingCase::FVI_MATCH_AND_LESST32:
#ifdef printd
				cout << "...BlockingCase::FVI_MATCH_AND_LESST32 \n";
#endif
				tmp  = (sizes[1] + params[0] -1)/params[0];
				for(int j = 2; j < n; j++)
				{
					if(perm_new[1] == j && sizes[0] < 128)
					{
						tmp *=  (sizes[j] + params[0] -1)/params[0];
					}
					else 
					{
						tmp *= sizes[j];
					}
				}
				params[6] = tmp;
					fvimatchl32_transpose_kernel(n, input, output, sizes, ldb,params, rperm, alpha, beta);

				break;
			case BlockingCase::FVI_MATCH_AND_GREATERT32:
#ifdef printd
				cout << "...BlockingCase::FVI_MATCH_AND_GREATERT32 \n";
#endif
				if(sizes[0] < 128)
				{
					params[0] = 4;
					tmp  = (sizes[1] + params[0] -1)/params[0];
				}
				else//need fix
				{
					params[0] = 1;
					tmp = sizes[1];	
					//tmp  = (sizes[1] + params[0] -1)/params[0];
				}
				for(int j = 2; j < n; j++)
				{
					if(perm_new[1] == j && sizes[0] < 128)
					{
						tmp *=  (sizes[j] + params[0] -1)/params[0];
					}
					else 
					{
						tmp *= sizes[j];
					}

				}
				params[6] = tmp;
#ifndef notall
				if(sizes[0] >= 128)
					fvimatchg32_transpose_kernel(n, input, output, sizes, ldb,params, perm_new, alpha, beta);
				else
					fvimatchg32_blocking_transpose_kernel(n, input, output, sizes, ldb,params, perm_new, alpha, beta);

#endif
				break;
			case BlockingCase::FVI_NOMATCH_GENERAL:
				params[3] = myParams.getBlockAIndex();
				//params[1] = myParams.getPaddingSize();
				params[4] = myParams.getBlockBIndex();
				params[7] = myParams.getNumElementsProcessedPerBlock();
				params[8] = myParams.getNumElementsProcessedPerBlock1();
				params[11] = myParams.getTileSize1();

				//cout <<"\n...Here "<< params[11];
				if(params[0] == 1)
				{
					tmp = 1;
				}
				else
				{
					tmp  = (sizes[params[3]] + params[0] -1)/params[0];
				}
#ifdef printd
				cout <<"\nBa = "<< params[3]<<"\tBb = "<<params[4]<<"\ttile1 = "<<params[0]<<"\ttile2 = "<<params[11]<<"\t tmp = "<<tmp<<"\n";
#endif
				for(int j = params[3]+1; j < n; j++)
				{
					if(rperm[j] < params[4])
					{
						//tmp *=  (sizes[j] + params[0] -1)/params[0];
					}
					else if(rperm[j] == params[4])
					{
						if(params[11] == 1)
							tmp *= 1;
						else
							tmp *=  (sizes[j] + params[11] -1)/params[11];
					}
					else 
					{
						tmp *= sizes[j];
					}
#ifdef printd
					cout <<"j = "<<j<<" tmp = "<<tmp<<"\n";
#endif

				}
				//cout <<"\ntmp = "<<tmp<<"\n";
				params[6] = tmp;
				params[10] = myParams.getSharedMemSize2();
				fvinomatchgeneral_transpose_kernel(n, input, output, sizes, ldb,params, perm_new, rperm, alpha, beta);
				break;
			case BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP:
				params[3] = myParams.getBlockAIndex();
				params[4] = myParams.getBlockBIndex();
				params[7] = myParams.getNumElementsProcessedPerBlock();
				params[8] = myParams.getNumElementsProcessedPerBlock1();
				params[11] = myParams.getTileSize1();
				tmp = 1;
				if(params[0] != 1)
				{
					if(rperm[params[3]] == params[4])
					{
						tmp  = (sizes[params[3]] + params[0] -1)/params[0];
					}
					else
					{
						tmp  = (sizes[params[3]] + params[0] -1)/params[0];
					}
				}
				//else tmp = ;
#ifdef printd
				cout <<"\nBa = "<< params[3]<<"\tBb = "<<params[4]<<"\ttile1 = "<<params[0]<<"\ttile2 = "<<params[11]<<"\t tmp = "<<tmp<<"\n";
#endif
				for(int j = params[3]+1; j < n; j++)
				{
					if(rperm[j] < params[4])
					{
					}
					else if(rperm[j] == params[4])
					{
						if(params[11] == 1)
						{
						}
						else
						{
							tmp *=  (sizes[j] + params[11] -1)/params[11];
						}
					}
					else 
					{
						tmp *= sizes[j];
					}
#ifdef printd
					cout <<"j = "<<j<<" tmp = "<<tmp<<"\n";
#endif
				}
				//cout <<"\ntmp = "<<tmp<<"\n";
				params[6] = tmp;
				params[10] = myParams.getSharedMemSize2();
				fvigeneralolap_transpose_kernel(n, input, output, sizes, ldb,params, perm_new, rperm, alpha, beta);
				break;
/*			case BlockingCase::FVI_NOMATCH_AND_GREATERT32:
#ifdef printd
				cout << "...BlockingCase::FVI_NOMATCH_AND_GREATERT32 \n";
#endif
				tmp  = (sizes[0] + params[0] -1)/params[0];
				for(int j = 1; j < n; j++)
				{
					if(perm_new[0] == j)
					{
						tmp *=  (sizes[j] + params[0] -1)/params[0];
					}
					else 
					{
						tmp *= sizes[j];
					}

				}
				params[6] = tmp;
				params[10] = myParams. getSharedMemSize2();
				params[7]= perm_new[0];
				params[8]= rperm[0];
				params[9]= 0;//perm_new[0];
				params[1] =  myParams.getNumElementsProcessedPerBlock();
#ifndef notall
				fvinomatchg32_transpose_kernel(n, input, output, sizes, ldb,params, rperm, alpha, beta);
#endif
				break;
				*/
		}	}	
}

void fuseIndices(TransposeSpec &spec){
	int ndim = spec.getNdim();
	int *perm = spec.getPermutation();
	int *sizes = spec.getSizes();
	int rperm[100];
	for(int i = 0; i <  ndim; i++)
	{
		for(int j = i+1; j <  ndim; j++)
		{
		if(perm[i] == j)
		{
			rperm[j] = i;
		}
		}
	}
	for(int i = 0; i <  ndim; i++)
	{
		if(sizes[i] == 1)
		{
			for(int j = i+1; j <  ndim; j++)
			{
				sizes[j-1] = sizes[j];
			}
			int ind = -1;
			for(int j = 0; j < ndim; j++)
			{
				if(perm[j] > i)
					perm[j]--;
				else if (perm[j] == i)
					ind = j;
			}
			//for(int j = rperm[i] + 1; j < ndim; j++)
			for(int j = ind + 1; j < ndim; j++)
			{
				perm[j-1] = perm[j];
					
			}
			ndim--;
			i--;
		}
	}
	for(int i = 0; i <  ndim-1; i++)
	{
		//if((perm[i] == i) && (perm[i+1] == i+1))
		if(perm[i] +1  == perm[i+1])
		{
			sizes[perm[i]] *= sizes[perm[i+1]];
			for(int j = perm[i+1] + 1; j < ndim; j++)
			{
				sizes[j-1] = sizes[j];
			}
			for(int j = 0; j < ndim; j++)
			{
				if(perm[j] > perm[i])
					perm[j]--;
			}
			for(int j = i+2; j < ndim; j++)
			{
				perm[j-1] = perm[j];
			}
			ndim--;
			i--;
		}
	}
	spec.setNdim(ndim);
}
void transpose_init(type *A, int size)
{
	int i = 0;
	for(i = 0; i < size; i++)
	{
		A[i] = i;
	}
}
int transpose_equal(const type *A, const type*B, int total_size){ 
	int error = 0; 
	const type *Atmp= A; 
	const type *Btmp= B; 
	for(int i=0;i < total_size ; ++i){ 
		type Aabs = (Atmp[i] < 0) ? -Atmp[i] : Atmp[i]; 
		type Babs = (Btmp[i] < 0) ? -Btmp[i] : Btmp[i]; 
		type max =  (Aabs < Babs) ? Babs : Aabs; 
		type diff = (Aabs - Babs); 
		diff = (diff < 0) ? -diff : diff; 
		if(diff > 0){ 
			type relError = (diff/max); 
			if(relError > 1.000000e-10){ 
				//printf("i: %d relError: %.8e %e %e\n",i,relError,Atmp[i], Btmp[i]); 
				//		             exit(0); 
				error += 1; 
			} 
		} 
	} 
	printf("Number of errors: %d\t",error); 
	return (error > 0) ? 0 : 1; 
} 

