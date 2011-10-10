/* pbe_exchange.c - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the PBE GGA
   exchange potential and energy of a spin density rho, defined on
   a log grid.

   THis GGA was developed by Perdew.

*/
#include	<stdio.h>
#include        "loggrid.h"
#include	"pbe_exchange.h"



/********************************
 *				*
 *         R_PBE96_Exchange	*
 *				*
 ********************************/

/* this routine calculates the PBE
   GGA exchange functional.

   Entry - rho[]: the density
   Exit  - Vx[]:  the exchange functional
	   Ex:	  the exchange energy
	   Px[]:  The variational exchange corrections for the eigenvalues.
*/

void  R_PBE96_Exchange(rho,Vx,Ex,Px)

double	rho[],
Vx[],
*Ex,
*Px;
{
    int	i;
    double mu,kappa;
    double onethird,fourthird,twothird;
    double pi;
    double n;
    double n_onethird;
    double ex_p;
    double ux_p;
    double s,F,Fs,Fss,kf,P0;
    double u,v;
    double agr,lap,delgr;

    /* loggrid variables */
    int	   Ngrid;

    /* tempory local grids */
    double *ex_functional;
    double *tmp;
    double *drho;
    double *ddrho;
    double *dadrho;
    double *rgrid;

    /* define constants */
    mu        = 0.2195149727645171;
    kappa     = 0.8040000000000000;
    pi        = 4.0*atan(1.0);
    onethird  = 1.0/3.0;
    fourthird = 4.0/3.0;
    twothird  = 2.0/3.0;


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();
    rgrid     = r_LogGrid();

    /* allocate temporary memory */
    ex_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();
    drho 	    = alloc_LogGrid();
    ddrho 	    = alloc_LogGrid();
    dadrho 	    = alloc_LogGrid();

    /* calculate drho,ddrho, */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = rho[i]/(4.0*pi);
    Derivative_LogGrid(tmp,drho);
    Derivative_LogGrid(drho,ddrho);

    /* calculate dadrho */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = fabs(drho[i]);
    Derivative_LogGrid(tmp,dadrho);



    for (i=0; i<Ngrid; ++i)
    {
        /* regular inputs to GGA */
        n     = rho[i]/(4.0*pi);
        if (n > 1.0e-18)
        {
            agr     = fabs(drho[i]);
            delgr   = drho[i]*dadrho[i];
            lap     = ddrho[i] + (2.0/rgrid[i])*drho[i];

            n_onethird = pow((3.0*n/pi),onethird);

            kf = pow( (3.0*pi*pi*n), onethird);
            s  = agr/(2.0*kf*n);
            u  = delgr/(n*n*(8.0*kf*kf*kf));
            v  = lap/(n*(4.0*kf*kf));
            P0 = 1.0 + (mu/kappa)*s*s;


            F   = (1.0 + kappa - kappa/P0);
            Fs  = 2.0*mu/(P0*P0);
            Fss = -4.0*(mu/kappa)*s*Fs/P0;

            ex_p = -(3.0/4.0)*n_onethird*F;
            ux_p = -(3.0/4.0)*n_onethird*(
                       fourthird*F
                       - v*Fs
                       - (u - fourthird*s*s*s)*Fss
                   );

            ex_functional[i] = ex_p;
            Vx[i]            = ux_p;
        }
        else
        {
            ex_functional[i] = 0.0;
            Vx[i]            = 0.0;
        }
    } /*for i*/


    /* cacluate Ex, and Px */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ex = integrate((rho_down+rho_down)*ex_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*ex_functional[i];
    *Ex = Integrate_LogGrid(tmp);

    /* integrate px_up = integrate(rho*Vx) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*Vx[i];
    *Px = Integrate_LogGrid(tmp);



    /* deallocate temporary memory */
    dealloc_LogGrid(ex_functional);
    dealloc_LogGrid(tmp);
    dealloc_LogGrid(drho);
    dealloc_LogGrid(ddrho);
    dealloc_LogGrid(dadrho);

} /* R_PBE_Exchange */
/* $Id$ */
