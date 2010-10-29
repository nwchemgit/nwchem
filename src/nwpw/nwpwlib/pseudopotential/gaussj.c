/*
 $Id$
*/
#include        <stdlib.h>
#include	<math.h>
#include	<stdio.h>
#include	"gaussj.h"

#define	NMAX	50

/********************************
 *				*
 *          gaussj		*
 *				*
 ********************************/

/*
   Gauss-Jordan routine used by pseudt to find polynominal
 constants.  Taken from Numerical Recipes, page 28.

*/

void  gaussj(int n, double a[],
             int m, double b[])
{
    int i,j,k,l,ll;
    int icol,irow;
    int ipiv[NMAX], indxr[NMAX], indxc[NMAX];
    double big,dum,pivinv;

    for (j=0; j<n; ++j) ipiv[j] = 0;

    for (i=0; i<n; ++i)
    {
        big=0.0;
        for (j=0; j<n; ++j)
        {
            if (ipiv[j] != 1)
            {
                for (k=0; k<n; ++k)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j+k*n]) >= big)
                        {
                            big  = fabs(a[j+k*n]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                    {
                        printf("Error: singular matrix in gaussj\n");
                        exit(99);
                    }
                } /*for k*/
            } /*if ipiv!=1*/
        } /*for j */

        ++ipiv[icol];

        if (irow != icol)
        {
            for (l=0; l<n; ++l)
            {
                dum         = a[irow+l*n];
                a[irow+l*n] = a[icol+l*n];
                a[icol+l*n] = dum;
            }
            for (l=0; l<m; ++l)
            {
                dum         = b[irow+l*m];
                b[irow+l*m] = b[icol+l*m];
                b[icol+l*m] = dum;
            }
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol+n*icol] == 0.0)
        {
            printf("Error: singular matrix in gaussj\n");
            exit(99);
        }

        pivinv         = 1.0/a[icol+icol*n];
        a[icol+icol*n] = 1.0;
        for (l=0; l<n; ++l) a[icol+l*n] *= pivinv;
        for (l=0; l<m; ++l) b[icol+l*m] *= pivinv;

        for (ll=0; ll<n; ++ll)
        {
            if (ll != icol)
            {
                dum          = a[ll+icol*n];
                a[ll+icol*n] = 0.0;
                for (l=0; l<n; ++l) a[ll+l*n] -= a[icol+l*n]*dum;
                for (l=0; l<m; ++l) b[ll+l*m] -= b[icol+l*m]*dum;
            }
        }/*for ll*/

    } /*for i*/

    for (l=(n-1); l>=0; --l)
    {
        if (indxr[l] != indxc[l])
        {
            for (k=0; k<n; ++k)
            {
                dum             = a[k+indxr[l]*n];
                a[k+indxr[l]*n] = a[k+indxc[l]*n];
                a[k+indxc[l]*n] = dum;
            }
        }
    }/*for l*/

} /*end gaussj*/

