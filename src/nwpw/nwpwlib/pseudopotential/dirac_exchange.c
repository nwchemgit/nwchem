/*
 $Id$
   dirac_exchange.c - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the Dirac
   exchange potential and energy of a spin density rho, defined on
   a log grid.

*/

#include	<stdio.h>
#include        "loggrid.h"
#include	"dirac_exchange.h"

static	double	alpha=0.6666666666666666667;

/********************************
 *			 	*
 *        set_Dirac_alpha	*
 *				*
 ********************************/

void set_Dirac_alpha(aa)
double aa;
{
    alpha = aa;
}

/********************************
 *				*
 *         Dirac_alpha		*
 *				*
 ********************************/

double Dirac_alpha()
{

    return alpha;
}



/********************************
 *				*
 *         R_Dirac_Exchange	*
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Dirac exchange functional.

   Entry - rho[]: the density
   Exit  - Vx[]:  the exchange functional
	   Ex:	  the exchange energy
	   Px:  The variational exchange corrections for the eigenvalues.
*/

void  R_Dirac_Exchange(rho,Vx,Ex,Px)

double	rho[],
Vx[],
*Ex,
*Px;
{
    int	i;
    double onethird;
    double pi;
    double n;
    double n_onethird;
    double ex_p;
    double ux_p;

    /* loggrid variables */
    int	   Ngrid;

    /* tempory local grids */
    double *ex_functional;
    double *tmp;

    /* define constants */
    onethird  = 1.0/3.0;
    pi       = 4.0*atan(1.0);


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();

    /* allocate temporary memory */
    ex_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();

    for (i=0; i<Ngrid; ++i)
    {
        n     = rho[i]/(4.0*pi);
        n_onethird = pow((3.0*n/pi),onethird);

        ex_p = -(9.0/8.0)*alpha*n_onethird;
        ux_p = -(3.0/2.0)*alpha*n_onethird;

        ex_functional[i] = ex_p;
        Vx[i] = ux_p;
    } /*for i*/


    /* cacluate Ex, and Px */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ex = integrate(rho*ex_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = rho[i]*ex_functional[i];
    *Ex = Integrate_LogGrid(tmp);

    /* integrate px = integrate(rho*Vx) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*Vx[i];
    *Px = Integrate_LogGrid(tmp);


    /* deallocate temporary memory */
    dealloc_LogGrid(ex_functional);
    dealloc_LogGrid(tmp);

} /* Dirac_Exchange */
