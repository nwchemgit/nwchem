/*
   $Id$
*/

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>

#include        "paw_pred_cor.h"
#include        "paw_loggrid.h"
#include        "paw_utilities.h"
#include        "paw_my_memory.h"
#include        "paw_sdir.h"

/***********************/
/* LogGrid data structure */
/***********************/

/* Hamman definitions */
/*static double amesh = 1.0247*/
/*static double r0Z   =  0.00625 */
/*static double Lmax  = 45.0  */
/*static double Lmax  = 45.0; */


int     Ngrid;
double  log_amesh;
double  amesh;
double  *rgrid;
double  *rgrid2;
double  *rgrid3;
double  *scratch;


double Lmax;
static double r0Z   = 0.00025;
double  r0;

void  paw_end_LogGrid()
{
    paw_dealloc_LogGrid(rgrid);
    paw_dealloc_LogGrid(rgrid2);
    paw_dealloc_LogGrid(rgrid3);
    paw_dealloc_LogGrid(scratch);
}


/****************************************
Function name	  : paw_init_LogGrid
Description	    :
Return type		  : void
Argument        : double Z -> ion charge
Argument        : FILE *fp
Author     		  : Marat Valiev
Date & Time		  : 1/7/99 4:26:57 PM
****************************************/
void  paw_init_LogGrid_from_file( double Z, FILE *fp)
{
    int  i;
    char input[30];

    strcpy(input,"<grid>");
    if (paw_find_word(input,fp) != 0)
    {
        if (paw_debug()) printf("Using default parameters\n");
        Lmax=25.0;
        amesh=1.005;
        log_amesh = log(amesh);

        r0 = r0Z/Z;
        Ngrid = (int) floor(log(Lmax/r0)/log_amesh)+1;

        /* make sure Ngrid is odd */
        if ((Ngrid%2)==0) Ngrid += 1;

    }
    else
    {
        fscanf(fp,"%le",&Lmax);
        fscanf(fp,"%d", &Ngrid);
        fscanf(fp,"%le",&r0);

        log_amesh = log(Lmax/r0)/(Ngrid-1);
        amesh     = exp(log_amesh);
    }

    Lmax = r0*pow(amesh,Ngrid-1);
    rgrid   = paw_alloc_LogGrid();
    rgrid2  = paw_alloc_LogGrid();
    rgrid3  = paw_alloc_LogGrid();
    scratch = paw_alloc_LogGrid();

    /* define rgrid */
    rgrid[0] = r0;
    for (i=1; i <= Ngrid-1; ++i)
        rgrid[i] = amesh*rgrid[i-1];

    for (i=0; i <= Ngrid-1; ++i)
    {
        rgrid2[i] = rgrid[i]*rgrid[i];
        rgrid3[i] = rgrid[i]*rgrid[i]*rgrid[i];
    }

}



/****************************************
  Function name	  : paw_alloc_LogGrid
  Description	    : creates a loggrid array
  Return type		  : double*
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:31:15 PM
****************************************/
double* paw_alloc_LogGrid()
{
    double *tt;

    tt = paw_alloc_1d_array(Ngrid);

    return tt;

} /* paw_alloc_LogGrid */




/****************************************
  Function name	  : paw_dealloc_LogGrid
  Description	    : deallocates the LogGrid array
  Return type		  : void
  Argument         : grid
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:31:39 PM
****************************************/
void  paw_dealloc_LogGrid(double *grid)
{

    paw_dealloc_1d_array(grid);

} /* dealloc_LogGrid */





/****************************************
Function name	  : paw_r_LogGrid
Description	    : returns the pointer to the rgrid array
Return type		  : double*
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:35:15 PM
****************************************/
double* paw_r_LogGrid()
{
    return rgrid;

} /* r_LogGrid */


double* paw_r2_LogGrid()
{
    return rgrid2;

}

double* paw_r3_LogGrid()
{
    return rgrid3;

}
/****************************************
  Function name	  : paw_N_LogGrid
  Description	    :
  Return type		  : int
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:35:40 PM
****************************************/
int  paw_N_LogGrid()
{
    return Ngrid;

} /* paw_N_LogGrid */




/****************************************
Function name	  : paw_r0_LogGrid
Description	    : returns the first nonzero coordinate
of a log grid
Return type		  : double 
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:36:10 PM
****************************************/
double paw_r0_LogGrid()
{
    return r0;
}



/****************************************
Function name	  : paw_log_amesh_LogGrid
Description	    : returns the value log(amesh)
Return type		  : double
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:37:11 PM
****************************************/
double paw_log_amesh_LogGrid()
{
    return log_amesh;

} /* paw_log_amesh_LogGrid */




/****************************************
  Function name	  : paw_amesh_LogGrid
  Description	    : returns the value (amesh)
  Return type		  : double
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:37:43 PM
****************************************/
double paw_amesh_LogGrid()
{
    return(amesh);

} /* paw_amesh_LogGrid */





/****************************************
  Function name	  : paw_Def_Integr
  Description	    : calculates definite integral of
  the function f with the weight
  r**(rpow) from 0 to Nrange.
  Function f is assumed to behave
  as r**(fpow) near 0
  Return type		  : double
  Argument         : double fpow
  Argument         : double *f
  Argument         : double rpow
  Argument         : int Nrange
  Author     		  : Marat Valiev
  Date & Time		  : 1/7/99 4:38:43 PM
****************************************/
double  paw_Def_Integr(double fpow,double *f,double rpow,int Nrange)

{
    int i;
    double sum;

    sum = (   9.0*f[0]*pow(rgrid[0],rpow+1)
              + 23.0*f[1]*pow(rgrid[1],rpow+1)
              + 28.0*f[2]*pow(rgrid[2],rpow+1)
          )/28.0;

    for (i=3; i<Nrange; ++i)
    {
        sum += f[i]*pow(rgrid[i],rpow+1);
    }

    sum = log_amesh*sum + f[0]*pow(rgrid[0],rpow+1)/(rpow+fpow+1);

    return sum;

}


/****************************************
Function name	  : paw_Integrate_LogGrid
Description	    :  returns a definite integral of
of the given function f times
r squared
Return type		  : double
Argument         : double f[]
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:40:59 PM
****************************************/

double  paw_Integrate_LogGrid(double *f)
{
    int    i;
    double sum;

    sum = (   9.0*f[0]*(rgrid[0]*rgrid[0]*rgrid[0])
              + 23.0*f[1]*(rgrid[1]*rgrid[1]*rgrid[1])
              + 28.0*f[2]*(rgrid[2]*rgrid[2]*rgrid[2])
          )/28.0;

    for (i=3; i<=Ngrid-1; i++)
    {
        sum += f[i]*(rgrid[i]*rgrid[i]*rgrid[i]);
    }

    sum = log_amesh*sum + f[0]*(rgrid[0]*rgrid[0]*rgrid[0])/3.0;

    return sum;

}




/****************************************
Function name	  :   paw_Zero_LogGrid
Description	    :
Return type		  : void
Argument         : double* grid
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:42:43 PM
****************************************/
void    paw_Zero_LogGrid(double *grid)

{
    int i;

    for (i=0; i<Ngrid; ++i)
        grid[i] = 0.0;

} /* paw_Zero_LogGrid */


/****************************************
  Function name	  : paw_Copy_LogGrid
  Description	    :
  Return type		  : void
  Argument         : double *gridnew
  Argument         : double *gridold
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:44:22 PM
****************************************/
void paw_Copy_LogGrid(double *gridnew, double *gridold)

{
    int i;
    for (i=0; i<Ngrid; ++i)
        gridnew[i] = gridold[i];
}



/****************************************
Function name	  : paw_Norm_LogGrid
Description	    : This routine calculates the Norm
of a wavefunction assuming that
the wavefunction decays like an
exponential as r goes to  infinity
Return type		  : double
Argument         : int M -> endpoint
Argument         : double gamma -> power of u near the 0
Argument         : double *u
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:46:08 PM
****************************************/
double paw_Norm_LogGrid(int M, double gamma, double *u)
{
    int   i;
    double sum;


    sum = (   9.0*u[0]*u[0]*rgrid[0]
              + 23.0*u[1]*u[1]*rgrid[1]
              + 28.0*u[2]*u[2]*rgrid[2]
          )/28.0;

    for (i=3; i<=M; ++i)
    {
        sum += u[i]*u[i]*rgrid[i];
    }

    sum = log_amesh*sum + u[0]*u[0]*rgrid[0]/(2.0*gamma+1.0);

    return sum;


} /* paw_Norm_LogGrid */



/****************************************
  Function name	  :   paw_Derivative_LogGrid
  Description	    : calculates the derivative
  of the function defined on the
  loggrid array
  Return type		  : void
  Argument         : double *f  -> original function
  Argument         : double *df -> derivative
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:47:29 PM
****************************************/
void    paw_Derivative_LogGrid(double *f,double *df)
{
    int i;

    df[0] = paw_Derivative5_1(0,f)/(log_amesh*rgrid[0]);
    df[1] = paw_Derivative5_2(1,f)/(log_amesh*rgrid[1]);

    for (i=2; i<Ngrid-2; ++i)
        df[i] = paw_Derivative5_3(i,f)/(log_amesh*rgrid[i]);


    df[Ngrid-2] = paw_Derivative5_4(Ngrid-2,f)/(log_amesh*rgrid[Ngrid-2]);
    df[Ngrid-1] = paw_Derivative5_5(Ngrid-1,f)/(log_amesh*rgrid[Ngrid-1]);

} /* paw_Derivative_LogGrid */


/****************************************
  Function name	  : paw_dot_product
  Description	    :
  Return type		  : double
  Argument         : double *f
  Argument         : double *g
  Author     		  : Marat Valiev
  Date & Time		  : 2/7/99 4:53:24 PM
****************************************/
double paw_dot_product(double *f, double *g)
{

    int k;
    double norm;

    norm =0.0;

    /* Integrate from 0 to r0 */
    norm = 0.5 * f[0] * g[0] * rgrid[0];

    for (k = 0; k < Ngrid; ++k)
        norm += f[k] * g[k] * rgrid[k] * log_amesh;

    return norm;

}

double paw_dot_product1(int n, double *f, double *g)
{

    int k;
    double norm;

    norm =0.0;

    /* Integrate from 0 to r0 */
    norm = 0.5 * f[0] * g[0] * rgrid[0];

    for (k = 0; k < n; ++k)
        norm += f[k] * g[k] * rgrid[k] * log_amesh;

    return norm;

}
/****************************************
 Function name	  : paw_print_loggrid_information
 Description	    :
 Return type		  : void
 Argument         : FILE *fp
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:00:49 PM
****************************************/
void  paw_print_loggrid_information(FILE *fp)
{
    fprintf(fp,"\n");
    fprintf(fp," Logarithmic grid information ( r(i)=r0*pow(a,i) ):\n");
    fprintf(fp,"\n");

    fprintf(fp,"   a    = %le\n", paw_amesh_LogGrid());
    fprintf(fp,"   N    = %d\n", paw_N_LogGrid());
    fprintf(fp,"   r0   = %le\n",paw_r0_LogGrid());
    fprintf(fp,"   rmax = %le\n",paw_r_LogGrid()[paw_N_LogGrid()-1]);
    fprintf(fp,"\n");

}


int paw_get_grid_index(double r)
{
    int i;

    if (r > Lmax)
    {

        printf("grid point is out of range\n");
        exit(1);

    }
    else
    {
        i = (int) floor( log(r/rgrid[0])/log_amesh );
    }

    return i;
}

double* paw_scratch_LogGrid()
{

    return scratch;

}


