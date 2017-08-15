__global__ void CON(fvimg32_kernel_, ndim)  (const type * A, type * B,const int  size0,
		const int * __restrict__ lda_s, const int* __restrict__ ldb_s,const int * __restrict__ idx_s, type alpha, type beta)
{
	int idx;
	int tmp; 
	int i;

	idx = blockIdx.x;
	int aexpr = 0, bexpr = 0;
#pragma unroll
	//printf("%d\n ", idx);
	for(i = 1; i < ndim; i++)
	{
		tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		idx = tmp;
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
	}
	int B_base = bexpr;
	int A_base = aexpr;

	int i0=threadIdx.x;
	double regs[8];
	
	if(size0 < 1024)//this condition being repeated is just for better performance
	{
#pragma unroll 8
		for(int r=0; r<8; r++)
		{
			if(i0 < size0)
			{
				regs[r] = A[i0 + A_base];
			}
			else
				break;
			i0+=128;
		}

		i0=threadIdx.x;
#pragma unroll 8
		for(int r=0; r<8; r++)
		{
			if(i0 < size0)
			{
				B[i0 + B_base] = alpha*regs[r] + beta*B[i0 + B_base];
			}
			else 
				break;
			i0+=128;
		}
	}
	else
	{
		int i0=threadIdx.x;
		for(int iii=0; iii < (size0+1023)/1024; iii++)
		{
			i0=threadIdx.x + iii*1024;
#pragma unroll 8
			for(int r=0; r<8; r++)
			{
				if(i0 < size0)
				{
					regs[r] = A[i0 + A_base];
				}
				else
					break;
				i0+=128;
			}
			i0=threadIdx.x + iii*1024;

#pragma unroll 8
			for(int r=0; r<8; r++)
			{
				if(i0 < size0)
				{
					B[i0 + B_base] = alpha*regs[r] + beta*B[i0 + B_base];
					//B[i0 + B_base] = regs[r];
				}
				else 
					break;
				i0+=128;

			}
		}
	}
	return;
}
#undef ndim

