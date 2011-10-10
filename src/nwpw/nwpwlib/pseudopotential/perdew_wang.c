/* perdew_wang.c
   Author - Eric Bylaska
*/
#include	"loggrid.h"
#include	"perdew_wang.h"

/* PBE96 coefficients */

/* Perdew-Wang92 LDA correlation coefficients */
#define	A		0.0310907
#define	A1		0.213700
#define	B1	 	7.595700
#define	B2		3.587600
#define	B3		1.638200
#define B4		0.492940


#define small_number	1.0e-80

/********************************
 *				*
 *   R_Perdew_Wang              *
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Perdew and Wange correlation functional.
   This is a Ceperly and Alder parameterization

   Entry - rho[]: the density
   Exit  - Vc_out[]:  the dependent exchange functional
	   Ec_out:	  the exchange energy
	   Pc_out:  The variational exchange corrections for the eigenvalues.
*/

void  R_Perdew_Wang(rho,Vc_out,Ec_out,Pc_out)

double	rho[],
Vc_out[],
*Ec_out,
*Pc_out;
{
    int	i;
    double onethird;
    double onesixth;
    double sevensixth;
    double pi,rs_scale;
    double rs,rss,n;
    double Q0,Q1,Q2,Q3;
    double ec,ec_rs;
    double uc;

    /* loggrid variables */
    int	   Ngrid;

    /* temporary local grids */
    double *ec_functional;
    double *tmp;

    /* define constants */
    pi       = 4.0*atan(1.0);
    onethird  = 1.0/3.0;
    onesixth  = 1.0/6.0;
    sevensixth  = 7.0/6.0;

    pi       = 4.0*atan(1.0);
    rs_scale = pow( (0.75/pi),  (onethird));


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();

    /* allocate temporary memory */
    ec_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();

    for (i=0; i<Ngrid; ++i)
    {
        /* regular inputs to GGA */
        n     = rho[i]/(4.0*pi) + small_number;

        /* calculate rs */
        rs    = rs_scale/pow(n,onethird);
        rss   = sqrt(rs);


        /* unpolarized LDA correlation energy */
        /* ec_p = correlation energy
           ec_p_rs = dec_p/drs
           uc_p    = dec_p/dn
        */
        Q0 = -2*A*(1.0+A1*rs);
        Q1 =  2*A*rss*(B1+rss*(B2+rss*(B3+B4*rss)));
        Q2 = log(1.0+1.0/Q1);
        Q3 = A*(B1/rss + 2.0*B2 + rss*(3.0*B3+4.0*B4*rss));

        ec    = Q0*Q2;
        ec_rs = -2.0*A*A1*Q2-Q0*Q3/(Q1*(1+Q1));
        uc    = ec - rs*ec_rs/3.0;

        ec_functional[i] = ec;
        Vc_out[i] = uc;

    } /*for i*/


    /* cacluate Ec, and Pc */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ec = integrate((rho)*ec_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*ec_functional[i];
    *Ec_out = Integrate_LogGrid(tmp);

    /* integrate pc = integrate(rho*Vc_out) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*Vc_out[i];
    *Pc_out = Integrate_LogGrid(tmp);


    /* deallocate temporary memory */
    dealloc_LogGrid(ec_functional);
    dealloc_LogGrid(tmp);

} /* R_PBE_Correlation */
/* $Id$ */
