/*
   $Id$
*/

#include	<stdio.h>
#include  <math.h>
#include  "paw_loggrid.h"
#include	"paw_dirac_exchange.h"
#include  "paw_my_constants.h"

static	double onethird;

static	double	alpha=0.6666666666666666667;
static	double *ex_functional;
static	double *Vx;

/****************************************
 Function name	  : paw_init_dirac_exchange()
 Description	    :
****************************************/
void paw_init_dirac_exchange()
{

    Vx		        = paw_alloc_LogGrid();
    ex_functional = paw_alloc_LogGrid();

    /* define constants */
    onethird  = 1.0/3.0;

}

/****************************************
 Function name	  : paw_generate_exchange_pot_LDA(double *rho)
 Description	    :
****************************************/
void paw_generate_exchange_pot_LDA(double *rho)

{
    int	i;
    double n;
    double n_onethird;
    double ux_p;

    /* loggrid variables */
    int	   Ngrid;


    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();


    for (i=0; i<=Ngrid-1; i++)
    {

        n = rho[i]/(4.0*PI);
        n_onethird = pow((3.0*n/PI),onethird);

        ux_p = -(3.0/2.0)*alpha*n_onethird;

        Vx[i] = ux_p ;

    } /*for i*/


}
double* paw_get_exchange_potential()
{

    return Vx;

}

/****************************************
 Function name	  : paw_get_exchange_energy_LDA(double *rho)
 Description	    :
****************************************/
double paw_get_exchange_energy_LDA(double *rho)

{
    int	i;
    double n;
    double n_onethird;
    double ex_p;
    double Ex;
    double *tmp;

    /* loggrid variables */
    int	   Ngrid;

    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();
    tmp		    = paw_scratch_LogGrid();


    for (i=0; i<Ngrid; ++i)
    {

        n = rho[i]/(4.0*PI);
        n_onethird = pow((3.0*n/PI),onethird);


        ex_p = -(9.0/8.0)*alpha*n_onethird;

        ex_functional[i] = ex_p;

    } /*for i*/



    for (i=0; i<Ngrid; ++i)
        tmp[i] = rho[i]*ex_functional[i];
    Ex = paw_Integrate_LogGrid(tmp);

    return Ex;

}


