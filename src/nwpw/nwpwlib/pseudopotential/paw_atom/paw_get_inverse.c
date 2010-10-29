/*
   $Id$
*/

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>

#include        "paw_my_memory.h"


/***************************************************

        Calculates the inverse of matrix a using

        Gauss-Jordan elimination method with partial

        pivoting.

        For reference see Num. Rec. section2.1

        (that book sucks by the way )

***************************************************/



void    paw_get_inverse(double **a, int matrix_size)





{

    int     i, j, l;

    int     n;

    int     irow;

    double  big;

    double  temp;

    double  pivinv;

    int    *index;



    n = matrix_size;



    index = (int *) malloc(n * sizeof(int));

    /* Start the main loop over the rows */

    for (i = 0; i < n; ++i)

    {

        big = 0.0;



        /* Find the pivot in column i */

        for (j = i; j < n; ++j)

        {

            if (fabs(a[j][i]) >= big)

            {

                big = fabs(a[j][i]);

                irow = j;

            }

        }



        if (big == 0.0)

        {

            printf("Failed to invert matrix\n");

            exit(99);

        }

        index[i] = irow;



        /* Swap the rows to put pivot elemnt on the diag */

        if (irow != i)

        {

            for (j = 0; j < n; ++j)

            {

                temp = a[irow][j];

                a[irow][j] = a[i][j];

                a[i][j] = temp;

            }

        }



        pivinv = 1.0 / a[i][i];

        a[i][i] = 1.0;



        for (j = 0; j < n; ++j)

            a[i][j] *= pivinv;

        for (l = 0; l < n; ++l)

        {

            if (l != i)

            {

                temp = a[l][i];

                a[l][i] = 0.0;

                for (j = 0; j < n; ++j)

                    a[l][j] -= a[i][j] * temp;

            }

        }

    }





    for (i = n - 1; i >= 0; i--)

    {

        if (index[i] != i)

        {

            for (j = 0; j < n; ++j)

            {

                temp = a[j][index[i]];

                a[j][index[i]] = a[j][i];

                a[j][i] = temp;

            }



        }

    }

}


void paw_test_matrix_inverse()
{

    int i,j,k;
    int n;

    double** a;
    double** b;
    double** c;

    n = 2;

    a = paw_alloc_2d_array(n,n);
    b = paw_alloc_2d_array(n,n);
    c = paw_alloc_2d_array(n,n);

    a[0][0] = 1.0;
    a[0][1] = 2.0;
    a[1][0] = 2.0;
    a[1][1] = 3.0;


    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
        {

            b[i][j] = a[i][j];
            c[i][j] = 0.0;

        }

    paw_get_inverse(b,2);

    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
            for (k=0;k<n;k++)
            {

                c[i][j] = c[i][j] +  a[i][k]*b[k][j];

            }

    for (j=0;j<n;j++)
        printf("%f\t%f\n",c[0][j],c[1][j]);


}
