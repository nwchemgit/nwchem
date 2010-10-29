/*
 $Id$
   gga_exchange.c - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the GGA
   exchange potential and energy of a spin density rho, defined on
   a log grid.

   THis GGA was developed by Perdew. see Parr and Yangs' book
   page 197.

*/
#include	<stdio.h>
#include        "loggrid.h"
#include	"gga_exchange.h"



/********************************
 *				*
 *         R_GGA_Exchange	*
 *				*
 ********************************/

/* this routine calculates the spin
   polarized GGA exchange functional.

   Entry - rho[]: the spin density
   Exit  - Vx[]:  the spin dependent exchange functional
	   Ex:	  the exchange energy
	   Px[]:  The variational exchange corrections for the eigenvalues.
*/

void  R_GGA_Exchange(rho,Vx,Ex,Px)

double	*rho[],
*Vx[],
*Ex,
Px[];
{
    int	i;
    double onethird,fourthird,twothird;
    double two_to_onethird;
    double pi,Cx,rs_scale;
    double nup,ndown,xi;
    double n_onethird;
    double A,f,df;
    double ex_p, ex_f;
    double ux_p, ux_f;
    double s,Fs,kf;

    /* loggrid variables */
    int	   Ngrid;

    /* tempory local grids */
    double *ex_functional;
    double *tmp;
    double *drho;

    /* define constants */
    onethird  = 1.0/3.0;
    fourthird = 4.0/3.0;
    twothird  = 2.0/3.0;
    two_to_onethird = pow(2.0,onethird);
    A               = two_to_onethird - 1.0;

    pi       = 4.0*atan(1.0);


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();

    /* allocate temporary memory */
    ex_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();
    drho 	    = alloc_LogGrid();

    /* calculate drho */
    for (i=0; i<Ngrid; ++i)
    {
        nup   = rho[0][i]/(4.0*pi);
        ndown = rho[1][i]/(4.0*pi);
        tmp[i] = (nup + ndown);
    }
    Derivative_LogGrid(tmp,drho);


    for (i=0; i<Ngrid; ++i)
    {
        nup   = rho[0][i]/(4.0*pi);
        ndown = rho[1][i]/(4.0*pi);
        n     = (nup + ndown);
        n_onethird = pow((3.0*n/pi),onethird);

        if (n > 0.0)
            xi = (nup - ndown)/n;
        else
            xi = 1.0;

        kf = pow( (3.0*pi*pi*n), onethird);
        s  = drho[i]/(2.0*kf*n);
        Fs = (1.0 + 1.296*(s*s) + 14.0*(s*s*s*s) + 0.2*(s*s*s*s*s*s));
        Fs = pow(Fs, (1.0/15.0));

        ex_p = -(3.0/4.0)*n_onethird*Fs;
        ux_p = -(3.0/2.0)*alpha*n_onethird;


        ex_f = two_to_onethird*ex_p;
        ux_f = fourthird*ex_f;
        f =  (  pow((1.0+xi),fourthird)
                + pow((1.0-xi),fourthird)
                - 2.0)/(2.0*A);

        df = fourthird*(  pow((1.0+xi),onethird)
                          - pow((1.0-xi),onethird))
             /(2.0*A);

        ex_functional[i] = ex_p + f*(ex_f - ex_p);

        Vx[0][i] = ux_p + f*(ux_f-ux_p) + (+1.0-xi)*(df*(ex_f-ex_p));
        Vx[1][i] = ux_p + f*(ux_f-ux_p) + (-1.0-xi)*(df*(ex_f-ex_p));

    } /*for i*/


    /* cacluate Ex, and Px */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ex = integrate((rho_down+rho_down)*ex_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[0][i] + rho[1][i])*ex_functional[i];
    *Ex = Integrate_LogGrid(tmp);

    /* integrate px_up = integrate(rho_up*Vx_up) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[0][i])*Vx[0][i];
    Px[0] = Integrate_LogGrid(tmp);

    /* integrate px_down = integrate(rho_down*Vx_down) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[1][i])*Vx[1][i];
    Px[1] = Integrate_LogGrid(tmp);


    /* deallocate temporary memory */
    dealloc_LogGrid(ex_functional);
    dealloc_LogGrid(tmp);

} /* Dirac_Exchange */
