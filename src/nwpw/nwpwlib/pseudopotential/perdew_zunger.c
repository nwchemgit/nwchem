/*
 $Id$
   perdew_zunger.c
   Author - Eric Bylaska
*/
#include	"loggrid.h"
#include	"perdew_zunger.h"

/* Perdew and Zunger coefficients */

/* low density limit: rs > 1.0 */
#define	nu_p 		-0.1423
#define nu_f   		-0.0843
#define	beta1_p 	1.0529
#define	beta1_f   	1.3981
#define	beta2_p 	0.3334
#define	beta2_f   	0.2611

/* high density limit: 0 <= rs <=1 */
#define	A_p		0.03110
#define	A_f		0.01555
#define	B_p		-0.0480
#define	B_f		-0.0269
#define	C_p		0.00200
#define C_f		0.00070
#define D_p 		-0.0116
#define D_f		-0.0048

#define small_number	1.0e-80

/********************************
 *				*
 *   correlation_perdew_zunger  *
 *				*
 ********************************/

/* this routine calculates the spin
   polarized Perdew and Zunger correlation functional.
   This is a Ceperly and Alder parameterization

   Entry - rho[]: the density
   Exit  - Vc[]:  the dependent exchange functional
	   Ec:	  the exchange energy
	   Pc:  The variational exchange corrections for the eigenvalues.
*/

void  R_Perdew_Zunger(rho,Vc,Ec,Pc)

double	rho[],
Vc[],
*Ec,
*Pc;
{
    int	i;
    double onethird,fourthird,twothird;
    double two_to_onethird;
    double pi,rs_scale;
    double rs,n;
    double log_rs,sqrt_rs;
    double denominator;
    double ec_p;
    double uc_p;

    /* loggrid variables */
    int	   Ngrid;

    /* tempory local grids */
    double *ec_functional;
    double *tmp;

    /* define constants */
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


        /* low density limit */
        if (rs >= 1.0)
        {
            sqrt_rs = sqrt(rs);

            denominator = (1.0 + beta1_p*sqrt_rs + beta2_p*rs);
            ec_p = nu_p/denominator;
            uc_p = ec_p
                   - onethird*nu_p*(0.5*beta1_p*sqrt_rs + beta2_p*rs)
                   /(denominator*denominator);
        }
        /* high denisity limit */
        else
        {
            log_rs  = log(rs);

            ec_p = A_p*log_rs + B_p + C_p*rs*log_rs + D_p*rs;
            uc_p = ec_p
                   - onethird*(A_p + C_p*(rs*log_rs + rs) + D_p*rs);
        }


        ec_functional[i] = ec_p;
        Vc[i] = uc_p;

    } /*for i*/


    /* cacluate Ec, and Pc */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */

    /* integrate Ec = integrate((rho)*ec_functional) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*ec_functional[i];
    *Ec = Integrate_LogGrid(tmp);

    /* integrate pc = integrate(rho*Vc) */
    for (i=0; i<Ngrid; ++i)
        tmp[i] = (rho[i])*Vc[i];
    *Pc = Integrate_LogGrid(tmp);


    /* deallocate temporary memory */
    dealloc_LogGrid(ec_functional);
    dealloc_LogGrid(tmp);

} /* R_Perdew_Zunger */
