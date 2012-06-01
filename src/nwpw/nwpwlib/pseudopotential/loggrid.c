/*
 $Id$

   LogGrid.c - 6/9/95
   author - Eric Bylaska

   This file contains the data structure for handeling numerics
on a logarithmic grid.  The grid is defined from 0 to 45.0,
with the grid points defined by:

        r(i) = (a**i)*r0

        with r0 = 0.00625/Z
             a  = 1.0247
             i = 0,1,...,N; N = log(7200.0*Z)/AL

*/

#include        <stdio.h>
#include	"grids.h"
#include	"pred_cor.h"
#include        "loggrid.h"

/***********************/
/* LogGrid data structure */
/***********************/


/* Hamman definitions */
/*static double amesh = 1.0247;   */
/*static double r0Z   =  0.00625; */
/*static double Lmax  = 45.0;     */

/* My definitions */
/*static double Lmax  = 45.0;    */
/*static double amesh = 1.0050;  */
/*static double r0Z   = 0.00025; */

static double Lmax  = 45.0;
static double amesh = 1.0050;
static double r0Z   = 0.00025;



static int     N;
static double  r0;
static double  log_amesh;
static double  *rgrid;




/********************************
 *                              *
 *        init_LogGrid          *
 *                              *
 ********************************/

/* sets up the log grid data structure,
   and returns the number grid points, N.

   Entry - Z:    charge of the system
*/

void  init_LogGrid(Z)
double  Z;
{
    int  i;



    /* define r0 */
    r0 = r0Z/Z;

    /* define log(amesh) */
    log_amesh = log(amesh);

    /* find N */
    N = log(Lmax/r0)/log_amesh;

    /* Initialize the grids data structure */
    init_Grids(N);

    /* allocate rgrid and a tmp grids */
    rgrid = alloc_Grid();


    /* define rgrid */
    rgrid[0] = r0;
    for (i=1; i<N; ++i)
        rgrid[i] = amesh*rgrid[i-1];

} /* init_LogGrid */

/********************************
 *                              *
 *         end_LogGrid          *
 *                              *
 ********************************/
void  end_LogGrid()
{
     dealloc_Grid(rgrid);
     end_Grids();
}


/********************************
 *                              *
 *        alloc_LogGrid         *
 *                              *
 ********************************/

/* returns the pointer to a loggrid array.

*/
double  *alloc_LogGrid()
{
    double *tt;

    tt = alloc_Grid();
    Zero_LogGrid(tt);
    return tt;

} /* alloc_LogGrid */


/********************************
 *                              *
 *        dealloc_LogGrid       *
 *                              *
 ********************************/

/* deallocates the LogGrid array

*/
void  dealloc_LogGrid(grid)
double	*grid;
{

    dealloc_Grid(grid);

} /* dealloc_LogGrid */



/********************************
 *                              *
 *          r_LogGrid           *
 *                              *
 ********************************/

/* returns the pointer to the rgrid array.

*/
double  *r_LogGrid()
{
    return rgrid;

} /* r_LogGrid */

/********************************
 *                              *
 *      index_r_LogGrid         *
 *                              *
 ********************************/

/* returns the pointer to the rgrid array.

*/
int  index_r_LogGrid(r)
double r;
{
    int index;

    index = log(r/r0)/log_amesh;
    return index;

} /* index_r_LogGrid */



/********************************
 *                              *
 *          N_LogGrid           *
 *                              *
 ********************************/

/* returns the size of a log grid.

*/

int  N_LogGrid()
{
    return N;

} /* N_LogGrid */



/********************************
 *                              *
 *      log_amesh_LogGrid       *
 *                              *
 ********************************/

/* returns the value log(amesh).

*/

double log_amesh_LogGrid()
{
    return log_amesh;

} /* log_amesh_LogGrid */



/********************************
 *                              *
 *      amesh_LogGrid           *
 *                              *
 ********************************/

/* returns the value (amesh).

*/

double amesh_LogGrid()
{
    return(amesh);

} /* amesh_LogGrid */



/********************************
 *				*
 *      Integrate_LogGrid       *
 *				*
 ********************************/

double	Integrate_LogGrid(f)
double  f[];
{
    int    i;
    double sum;

    /*
       sum = (   9.0*f[0]*(rgrid[0]*rgrid[0]*rgrid[0])
              + 23.0*f[1]*(rgrid[1]*rgrid[1]*rgrid[1])
              + 28.0*f[2]*(rgrid[2]*rgrid[2]*rgrid[2])
             )/28.0;
    */
    sum = (   9.0*f[0]*(rgrid[0]*rgrid[0]*rgrid[0])
              + 28.0*f[1]*(rgrid[1]*rgrid[1]*rgrid[1])
              + 23.0*f[2]*(rgrid[2]*rgrid[2]*rgrid[2])
          )/24.0;
    for (i=3; i<N; ++i)
    {
        sum += f[i]*(rgrid[i]*rgrid[i]*rgrid[i]);
    }
    sum = log_amesh*sum + f[0]*(rgrid[0]*rgrid[0]*rgrid[0])/3.0;

    return sum;

}

/********************************
 *				*
 *  Integrate_LogGrid_na_nb     *
 *				*
 ********************************/

double	Integrate_LogGrid_na_nb(int na, int nb, double f[])
{
    int    i;
    double sum;

    sum = (   9.0*f[na  ]*(rgrid[na  ])
              + 28.0*f[na+1]*(rgrid[na+1])
              + 23.0*f[na+2]*(rgrid[na+2])
          )/24.0;
    for (i=na+3; i<=(nb-3); ++i)
    {
        sum += f[i]*(rgrid[i]);
    }
    sum += (  23.0*f[nb-2]*(rgrid[nb-2])
              + 28.0*f[nb-1]*(rgrid[nb-1])
              +  9.0*f[nb  ]*(rgrid[nb  ])
           )/24.0;

    sum = log_amesh*sum;

    return sum;

}

/********************************
 *				*
 *         Zero_LogGrid		*
 *				*
 ********************************/

void	Zero_LogGrid(grid)
double *grid;
{
    int i;

    for (i=0; i<N; ++i)
        grid[i] = 0.0;

} /* Zero_LogGrid */

/********************************
 *				*
 *       Copy_LogGrid		*
 *				*
 ********************************/

void Copy_LogGrid(gridnew,gridold)

double *gridnew,
*gridold;
{
    int i;
    for (i=0; i<N; ++i)
        gridnew[i] = gridold[i];
}


/**********************************
 *				  *
 *       Norm_LogGrid 	   	  *
 *				  *
 **********************************/

/* This routine calculates the Norm
   of a wavefunction assuming that
   the wavefunction decays like an
   exponential as r-> infinity
   and approaches zero as r**(2*gamma)

*/
double Norm_LogGrid(M,gamma,u)
int    M;
double gamma;
double u[];
{
    int	  i;
    double r0,sum;


    /* Find Integral(u**2) */
    r0 = rgrid[0]/sqrt(amesh);
    sum = pow(r0,(2.0*gamma+1.0))/(2.0*gamma+1.0);
    for (i=0; i<=(M-3); ++i)
        sum += log_amesh*rgrid[i]*(u[i]*u[i]);

    sum += log_amesh*(  23.0*rgrid[M-2]*(u[M-2]*u[M-2])
                        + 28.0*rgrid[M-1]*(u[M-1]*u[M-1])
                        +  9.0*rgrid[M]  *(u[M]  *u[M]))/24.0;

    return sum;

} /* Norm_LogGrid */


/**********************************
 *				  *
 *       Integrate2_LogGrid 	  *
 *				  *
 **********************************/

/* This routine calculates the integral fo function f,
   assuming that f approaches zero as r**nu, and decays
   like an exponential as r-> infinity

*/
double Integrate2_LogGrid(M,nu,f)
int    M;
double nu;
double f[];
{
    int	  i;
    double r0,sum;


    /* Find Integral of f */
    r0 = rgrid[0]/sqrt(amesh);
    sum = pow(r0,(nu+1.0))/(nu+1.0);
    for (i=0; i<=(M-3); ++i)
        sum += log_amesh*rgrid[i]*(f[i]);

    sum += log_amesh*(  23.0*rgrid[M-2]*(f[M-2])
                        + 28.0*rgrid[M-1]*(f[M-1])
                        +  9.0*rgrid[M]  *(f[M]  ))/24.0;

    return sum;

} /* Integrate2_LogGrid */

/********************************
 *				*
 *       Derivative_LogGrid	*
 *				*
 ********************************/

void	Derivative_LogGrid(f,df)

double	f[],
df[];
{
    int i;

    /* define dV/dr */
    df[0] = Derivative5_1(0,f)/(log_amesh*rgrid[0]);
    df[1] = Derivative5_2(1,f)/(log_amesh*rgrid[1]);
    for (i=2; i<N-2; ++i)
        df[i] = Derivative5_3(i,f)/(log_amesh*rgrid[i]);
    df[N-2] = Derivative5_4(N-2,f)/(log_amesh*rgrid[N-2]);
    df[N-1] = Derivative5_5(N-1,f)/(log_amesh*rgrid[N-1]);

} /* Derivative_LogGrid */

void	Plot_LogGrid(char *name, double *f)
{
    int i;
    FILE *fp;

    fp = fopen(name,"w+");
    for (i=0; i<N; ++i)
        fprintf(fp,"%le  %le\n",rgrid[i],f[i]);
    fclose(fp);
}
