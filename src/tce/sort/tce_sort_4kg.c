#include <stdio.h>


/*
 * NOTE: 
 *        A and B must have nDim dimensions
 *        The size of A (dim0 * dim1 * ....  * dimN) must equal the 
 *           size of B
 *        aSize and perm must be single dimension arrys of size nDim
 *        aSize must contain the size of each dimension of A
 *        perm must contain the dimension permutation (i.e, let A has dimension sequence
 *         [0, 1, 2, 3] and size [2, 3, 4, 5]. The permutation array [2, 1, 3, 0] will
 *         cause nDTranspose to produce an array of dimension [4, 3, 5, 2]
 *
 */

#include "typesf2c.h"

void FATR tce_sort_4kg__(double * A, double * B, Integer * rowIdx, Integer * colIdx, 
			Integer * thirdIdx, Integer * fourthIdx, Integer * perm1, Integer * perm2,
			Integer * perm3, Integer * perm4, double * factor)
{
	Integer   i     = 0;
	Integer   j     = 0;
	Integer   k     = 0;
	Integer   l     = 0;

	Integer   aSize[4];
	Integer   bSize[4];
	Integer   perm[4];

	Integer   aJump[4];
	Integer   bJump[4];

	double * bPtr = B;
	double * aPtr[4];
	
	
	aSize[0] = *rowIdx;
	aSize[1] = *colIdx;
	aSize[2] = *thirdIdx;
	aSize[3] = *fourthIdx;

	aJump[0] = 1;
	for (i = 1; i < 4; i++)
	{
		aJump[i] = aJump[i-1] * aSize[i-1];
	}
  
	perm[0] = *perm1;
	perm[1] = *perm2;
	perm[2] = *perm3;
	perm[3] = *perm4;

	/* assuming column ordering, the -first stride of A will be 1, the second */
	for (i = 0; i < 4; i++)
	{
		bSize[i] = aSize[perm[i]-1];
		bJump[i] = aJump[perm[i]-1];
		aPtr[i] = A;
	}


	for (i = 0; i < bSize[3]; i++)
	{
		aPtr[2] = aPtr[3];
		for (j = 0; j < bSize[2]; j++)
		{
			aPtr[1] = aPtr[2];
			for (k = 0; k < bSize[1]; k++)
			{
				aPtr[0] = aPtr[1];
				for (l = 0; l < bSize[0]; l++)
				{
					*bPtr++ = *aPtr[0];
					aPtr[0] += bJump[0];
				}
				aPtr[1] += bJump[1];
			}
			aPtr[2] += bJump[2];
		}
		aPtr[3] += bJump[3];
	}

}


/* $Id$ */
