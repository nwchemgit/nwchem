/* revpbe_correlation.c
   Author - Eric Bylaska
   $Id$

*/
#include	"loggrid.h"
#include	"revpbe_correlation.h"

/* revPBE coefficients */

/* Perdew-Wang92 LDA correlation coefficients */
#define	A		0.0310907
#define	A1		0.213700
#define	B1	 	7.595700
#define	B2		3.587600
#define	B3		1.638200
#define B4		0.492940

/* PBE96  GGA correlation coefficients */
#define GAMMA	0.031090690869655
#define	BETA	0.066724550603149

#define small_number	1.0e-80

/********************************
 *				*
 *   R_revPBE_correlation          *
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Perdew and Zunger correlation functional.
   This is a Ceperly and Alder parameterization

   Entry - rho[]: the density
   Exit  - Vc_out[]:  the dependent exchange functional
	   Ec_out:	  the exchange energy
	   Pc_out:  The variational exchange corrections for the eigenvalues.
*/

void  R_revPBE_Correlation(rho,Vc_out,Ec_out,Pc_out)

double	rho[],
Vc_out[],
*Ec_out,
*Pc_out;
{
    int	i;
    double BOG;
    double onethird;
    double onesixth;
    double sevensixth;
    double pi,rs_scale;
    double rs,rss,n;
    double lap,agr,delgr;
    double t,t2,t4,t6;
    double Q0,Q1,Q2,Q3,Q4,Q5,Q8,Q9;
    double uu,vv;
    double ks,kf;
    double H,Hrs,Hrst,Ht,Htt,H_B,H_Bt,B_ec,B;
    double FACT0;
    double FACT1;
    double FACT2;
    double FACT3;
    double ec,ec_rs;
    double uc,duc;

    /* loggrid variables */
    int	   Ngrid;
    double *rgrid;

    /* temporary local grids */
    double *ec_functional;
    double *tmp;
    double *drho;
    double *ddrho;
    double *dadrho;

    /* define constants */
    pi       = 4.0*atan(1.0);
    BOG = BETA/GAMMA;
    onethird  = 1.0/3.0;
    onesixth  = 1.0/6.0;
    sevensixth  = 7.0/6.0;

    pi       = 4.0*atan(1.0);
    rs_scale = pow( (0.75/pi),  (onethird));


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();
    rgrid     = r_LogGrid();

    /* allocate temporary memory */
    ec_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();
    drho		    = alloc_LogGrid();
    ddrho	    = alloc_LogGrid();
    dadrho	    = alloc_LogGrid();

    /* calculate drho,ddrho */
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
        n     = rho[i]/(4.0*pi) + small_number;
        if (n > 1.0e-12)
        {
            agr   = fabs(drho[i]);
            delgr = drho[i]*dadrho[i];
            lap   = ddrho[i] + (2.0/rgrid[i])*drho[i];

            /* calculate rs */
            rs    = rs_scale/pow(n,onethird);
            rss   = sqrt(rs);

            /* calculate t, uu, and vv */
            kf = pow( (3.0*pi*pi*n), onethird);
            ks = sqrt(4.0*kf/pi);

            t  = agr/(2.0*ks*n);
            uu = delgr/(8*ks*ks*ks*n*n);
            vv = lap/(4*ks*ks*n);



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


            /* PBE96 correlation energy  corrections */
            t2 = t*t;
            t4 = t2*t2;
            B = -ec/GAMMA;
            B = BOG/(exp(B)-1.0);
            Q4 = 1.0 + B*t2;
            Q5 = 1.0 + B*t2 + B*B*t4;
            H = GAMMA*log(1.0 + BOG*Q4*t2/Q5);

            /* PBE96 correlation potential corrections */
            t6   = t4*t2;

            B_ec = (B/BETA)*(BOG+B);

            Q8   = Q5*Q5+BOG*Q4*Q5*t2;
            H_B  = -BETA*B*t6*(2.0+B*t2)/Q8;

            Hrs  = -(rs/3.0)*H_B*B_ec*ec_rs;

            Q9    = 1.0+2*B*t2;
            FACT0 = 2.0*BOG-6.0*B;
            FACT1 = Q5*Q9 + Q4*Q9*Q9;
            FACT2 = Q4*Q5 + B*t2*(Q4*Q9+Q5);
            FACT3 = 2.0*B*Q5*Q9 + BOG*FACT2;
            H_Bt  = 2.0*BETA*t4*((Q4*Q5*FACT0-BOG*FACT1)/Q8)/Q8;

            Hrst = (rs/3.0)*t2*H_Bt*B_ec*ec_rs;

            Ht  = 2.0*BETA*Q9/Q8;
            Htt = 3.0*BETA*t*(2.0*B/Q8-(Q9*FACT3/Q8)/Q8);

            duc = H + Hrs + Hrst + t2*Ht/6.0 + 7.0*t2*t*Htt/6.0
                  -  uu*Htt - vv*Ht;

            ec_functional[i] = ec + H;
            Vc_out[i] = uc + duc;
        }
        else
        {
            ec_functional[i] = 0.0;
            Vc_out[i] = 0.0;
        }

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
    dealloc_LogGrid(drho);
    dealloc_LogGrid(ddrho);
    dealloc_LogGrid(dadrho);

} /* R_revPBE_Correlation */
