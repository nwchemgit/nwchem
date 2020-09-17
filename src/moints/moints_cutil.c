/*---------------------------------------------------------*\
$Id$
\*---------------------------------------------------------*/

#include <stdio.h>
#include "typesf2c.h"

void c_print_sparsemat( int, int, int *, int *, int, double * );



#if defined(EXT_INT)
long onbitmask_( long *len )
{
  unsigned long		mask;

  mask = ~((~0) << *len);
  return ((long)mask);
}
#else  
#if (defined(WIN32)) &&!defined(__MINGW32__)
int FATR ONBITMASK( int *len )
#else
int FATR onbitmask_( int *len )
#endif
{
  unsigned int		mask;

  mask = ~((~0) << *len);
  return ((int)mask);
}
#endif

  



#ifdef NOCOMPILE

void print_sparsemat_( int *nc, int *nr, int *cpi, int *ir, int *nnz, double *v )
{
    c_print_sparsemat( *nc, *nr, cpi, ir, *nnz, v );
}




void c_print_sparsemat( int nc, int nr, int *cpi, int *ir, int nnz, double *v )
{
    int		ic, colmax, cvlo, cvhi, ncv;
    int		npr, hasprint, ii2r, iiv, iir, mask16bit;
    int		ilab;

    colmax = nc < 6 ? nc : 6;
    mask16bit = ~((~0) << 16);
/*
    printf("\n");
    for (ic=0; ic<colmax; ++ic) {
	printf("   %3d %3d  ", cpi[2*ic], cpi[2*ic+1] );
    }
    printf("\n");
*/
    ii2r = sizeof(int)/2;             /* num of labels packed per int - 16 bits per label */
    npr = 1;
    hasprint = 1;
    while (hasprint) {
	hasprint = 0;
	for (ic=0; ic<colmax; ++ic) {
	    cvlo = cpi[2*ic];
	    cvhi = cpi[2*ic+1];
	    ncv = cvhi - cvlo + 1;
	    if ((cvlo>0)&&(ncv>=npr)) {
		iiv = cvlo + npr - 1;
		iir = iiv/ii2r + ((iiv%ii2r) ? 1 : 0);
		ilab = (ir[iir-1] >> (16*(iiv%ii2r))) & mask16bit;
		printf("   %2d %8.4f", ilab, v[iiv-1] );
		++hasprint;
	    }
	    else {
		printf("              ");
	    }
	}
	printf("\n");
	++npr;
    }
}
	







#endif
