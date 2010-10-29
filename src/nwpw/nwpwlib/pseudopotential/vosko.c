/*
 $Id$
   vosko.c
   Author - Eric Bylaska
*/

#include	"loggrid.h"
#include	"vosko.h"


/*---- parameters given by vosko et al -----------------*/
#define ap  3.109070e-02
#define af  1.554530e-02
#define x0p -1.049800e-01
#define x0f -3.250000e-01
#define bp  3.727440e+00
#define bf  7.060420e+00
#define cp  1.293520e+01
#define cf  1.805780e+01
/*------------------------------------------------------*/

/*     constants calculated from vosko's parameters */
#define xp   -4.581653e-01
#define xf   -5.772521e-01
#define qp    6.151991e+00
#define qf    4.730927e+00
#define xx0p  1.255491e+01
#define xx0f  1.586879e+01
#define cp1   3.109070e-02
#define cf1   1.554530e-02
#define cp2   9.690228e-04
#define cf2   2.247860e-03
#define cp3   1.049800e-01
#define cf3   3.250000e-01
#define cp4   3.878329e-02
#define cf4   5.249122e-02
#define cp5   3.075995e+00
#define cf5   2.365463e+00
#define cp6   1.863720e+00
#define cf6   3.530210e+00
#define dp1   6.218140e-02
#define df1   3.109060e-02
#define dp2   1.938045e-03
#define df2   4.495720e-03
#define dp3   1.049800e-01
#define df3   3.250000e-01
#define dp4  -3.205972e-02
#define df4  -1.779316e-02
#define dp5  -1.192972e-01
#define df5  -1.241661e-01
#define dp6   1.863720e+00
#define df6   3.530210e+00
#define dp7   9.461748e+00
#define df7   5.595417e+00
#define fc    1.923661e+00
#define fd    2.564881e+00
#define crs   7.876233e-01

#define	small_number	1.0e-80

/********************************
 *				*
 *   correlation_vosko et. al.  *
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Vosko et. al. correlation functional.
   This is a Ceperly and Alder parameterization

   Entry - rho: the density
   Exit  - Vc[]:  the exchange functional
	   Ec:	  the exchange energy
	   Pc:  The variational exchange corrections for the eigenvalues.
*/

void  R_Vosko(rho,Vc,Ec,Pc)

double	rho[],
Vc[],
*Ec,
*Pc;
{
    int	i;
    double onesixth,onethird,fourthird,twothird;
    double two_to_onethird;
    double pi,rs_scale;
    double rs,n;
    double x,xxp,dxxp;
    double ec_p;
    double uc_p;

    /* loggrid variables */
    int	   Ngrid;

    /* tempory local grids */
    double *ec_functional;
    double *tmp;

    /* define constants */
    onesixth  = 1.0/6.0;
    onethird  = 1.0/3.0;
    fourthird = 4.0/3.0;
    twothird  = 2.0/3.0;
    two_to_onethird = pow(2.0,onethird);

    pi       = 4.0*atan(1.0);
    rs_scale = pow( (0.75/pi), (1.0/3.0));


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();

    /* allocate temporary memory */
    ec_functional    = alloc_LogGrid();
    tmp		    = alloc_LogGrid();

    for (i=0; i<Ngrid; ++i)
    {
        n     = rho[i]/(4.0*pi) + small_number;
        rs    = rs_scale/pow(n,onethird);
        x     = sqrt(rs);

        xxp  = rs + bp*x + cp;
        dxxp = 2.0*x + bp;
        ec_p = cp1*log(rs/xxp) + cp2*log( (x+cp3)*(x+cp3)/xxp)
               + cp4*atan(cp5/(x+cp6));
        uc_p = ec_p
               - onesixth*x*(  dp1/x + dp2/(x+cp3) + dp4*dxxp/xxp
                               + dp5/( (x+dp6)*(x+dp6)+dp7)
                            );


        ec_functional[i] = ec_p;
        Vc[i] = uc_p;
    } /*for i*/


    /* cacluate Ec, and Pc */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ec = integrate((rhon)*ec_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*ec_functional[i];
    *Ec = Integrate_LogGrid(tmp);

    /* integrate pc = integrate(rhop*Vc) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*Vc[i];
    *Pc = Integrate_LogGrid(tmp);


    /* deallocate temporary memory */
    dealloc_LogGrid(ec_functional);
    dealloc_LogGrid(tmp);

} /* R_Vosko */
