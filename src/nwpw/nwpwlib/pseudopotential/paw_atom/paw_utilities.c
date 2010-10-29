/*
   $Id$
*/

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include       "paw_utilities.h"
#include       "paw_my_memory.h"
#include       "paw_sdir.h"

static  char    word[50];


/****************************************
 Function name	  :   *paw_get_word
 Description	    :
 Return type		  : char
 Argument         : FILE *stream
 Author     		  : Eric Bylaska
 Date & Time		  : unknown
****************************************/

char    *paw_get_word(FILE *stream)
{
    if (fscanf(stream,"%s",word) != EOF)
        return word;
    else
        return NIL;
}


/****************************************
 Function name	  : paw_find_word
 Description	    : gives the location of the word
                    in the file;
                    if successfull returns 0 and places
                    the file pointer right after the word;
                    if fails returns 1
 Return type		  : int
 Argument         : char* my_word
 Argument         : FILE *fp
 Author     		  : Marat Valiev
 Date & Time		  : 1/7/99 3:58:25 PM
****************************************/
int paw_find_word(char* my_word, FILE *fp)
{

    char *w;

    rewind(fp);
    w = paw_get_word(fp);
    if (strcmp(my_word,w)!=0)
    {
        rewind(fp);
        w = paw_get_word(fp);
        while ((w!=NIL) && (strcmp(my_word,w)!=0))
            w = paw_get_word(fp);
    }


    if (w==NIL)
    {
        if (paw_debug()) printf("warning: %s section not found\n",my_word);
        return 1;
    }
    else
    {
        return 0;
    }

}


/****************************************
 Function name	  : *paw_spd_Name
 Description	    :
 Return type		  : char
 Argument         : int l
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:01:44 PM
****************************************/
char *paw_spd_Name(int l)
{
    char *s;
    if (l==0)
        s = "s";
    else if (l==1)
        s = "p";
    else if (l==2)
        s = "d";
    else if (l==3)
        s = "f";
    else if (l==4)
        s = "g";
    else
        s = "?";

    return s;
}

void paw_lu_decompose(int N, double **c, double **a, double **b)
{

    int i;
    int j;
    int k;
    double sum;

    for (j=0;j<=N-1;j++)
    {
        for (i=0;i<=j;i++)
        {

            a[i][j] = 0.0;
            b[i][j] = 0.0;

        }
    }

    for (i=0;i<=N-1;i++)
        a[i][i] = 1.0;

    for (j=0;j<=N-1;j++)
    {
        for (i=0;i<=j;i++)
        {
            sum = c[i][j];
            for (k=0;k<=i-1;k++)
            {

                sum = sum - a[i][k]*b[k][j];

            }

            b[i][j] = sum;
        }

        for (i=j+1;i<=N-1;i++)
        {
            sum = c[i][j];
            for (k=0;k<=j-1;k++)
            {

                sum = sum + a[i][k]*b[k][j];

            }

            if (b[j][j]==0)
            {
                printf("error, division by zero in lu decomposition\n");
                exit(1);
            }
            else
            {
                a[i][j] = sum/b[j][j];
            }

        }

    }

}

void paw_triang_matrix_inverse(char* matrix_type, int N, double **a, double **a_inv)
{

    int i;
    int j;
    int k;
    double sum;

    if (strcmp(matrix_type,"u")==0)
    {

        for (j=0;j<=N-1;j++)
        {
            for (i=0;i<=N-1;i++)
            {

                a_inv[i][j] = a[j][i];

            }
        }

    }
    else if (strcmp(matrix_type,"l")==0)
    {

        for (j=0;j<=N-1;j++)
        {
            for (i=0;i<=N-1;i++)
            {

                a_inv[i][j] = a[i][j];

            }
        }

    }
    else
    {

        printf("unknown matrix type in triang_matrix_inverse\n");
        exit(1);

    }

    for (i=0;i<=N-1;i++)
    {
        a_inv[i][i] = 1.0/a_inv[i][i];
        for (j=i+1;j<=N-1;j++)
        {
            sum = 0.0;
            for (k=i;k<=j-1;k++)
            {
                sum = sum - a_inv[j][k]*a_inv[k][i];
            }

            a_inv[j][i] = sum/a_inv[j][j];
        }
    }

    /*
      DO i=1,n
        a(i,i)=1.0_DP/a(i,i)
        DO j=i+1,n
          s=0
          DO k=i,j-1

            s = s - a(j,k)*a(k,i)

          END DO

          a(j,i) = s/a(j,j)

        END DO

      END DO */

    if (strcmp(matrix_type,"u")==0)
    {

        for (i=0;i<=N-1;i++)
        {
            for (j=i+1;j<=N-1;j++)
            {
                a_inv[i][j] = a_inv[j][i];
                a_inv[j][i] = 0.0;
            }
        }

    }

}
void paw_test_lu_decompose()
{

    int i,j,k;
    int n;

    double** a;
    double** b;
    double** c;
    double** c1;

    n = 3;

    a = paw_alloc_2d_array(n,n);
    b = paw_alloc_2d_array(n,n);
    c = paw_alloc_2d_array(n,n);
    c1 = paw_alloc_2d_array(n,n);
    /*
      c[0][0] = 1.0;
      c[0][1] = 2.0;
      c[1][0] = 2.0;
      c[1][1] = 3.0;

      */

    c[0][0] = 1.0; c[0][1] = 2.0; c[0][2] = 3.0;
    c[1][0] = 3.0; c[1][1] = 2.0; c[1][2] = 0.0;
    c[2][0] = 1.0; c[2][1] = 2.0; c[2][2] = 8.0;




    printf("*****lu decomposition test************\n");


    paw_lu_decompose(n, c, a, b);

    printf("*****original matrix************\n");

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",c[i][j]);

        }

        printf("\n");

    }

    printf("*****u matrix************\n");

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",b[i][j]);

        }

        printf("\n");

    }

    printf("*****l matrix************\n");

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",a[i][j]);

        }

        printf("\n");

    }

    printf("*****lu product************\n");

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            c1[i][j] = 0.0;
            for (k=0;k<n;k++)
            {

                c1[i][j] = c1[i][j] + a[i][k]*b[k][j];

            }

        }

        printf("\n");

    }

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",c1[i][j]);

        }

        printf("\n");

    }

}

void paw_test_triang_matrix_inverse()
{

    int i,j,k;
    int n;

    double** a;
    double** b;
    double** c;

    n = 3;

    a = paw_alloc_2d_array(n,n);
    b = paw_alloc_2d_array(n,n);
    c = paw_alloc_2d_array(n,n);

    a[0][0] = 1.0; a[0][1] = 2.0; a[0][2] = 3.0;
    a[1][0] = 0.0; a[1][1] = 3.0; a[1][2] = 6.3;
    a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 4.0;



    paw_triang_matrix_inverse("u",n,a,b);

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            c[i][j] = 0.0;
            for (k=0;k<n;k++)
            {

                c[i][j] = c[i][j] +  a[i][k]*b[k][j];

            }
        }
    }

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",c[i][j]);

        }

        printf("\n");

    }


}


void paw_square_matrix_product(int n, double **a,double **b,double **c)
{
    int i,j,k;

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            c[i][j] = 0.0;
            for (k=0;k<n;k++)
            {

                c[i][j] = c[i][j] +  a[i][k]*b[k][j];

            }
        }
    }


}


void paw_print_matrix(int n, double **a, char* title)
{
    int i,j;

    printf("\n\n %s \n",title);

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {

            printf("%f\t",a[i][j]);

        }

        printf("\n");

    }
}

