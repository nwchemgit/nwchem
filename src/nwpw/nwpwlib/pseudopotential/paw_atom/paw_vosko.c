/*
   $Id$
*/

#include	<math.h>
#include	"paw_loggrid.h"
#include	"paw_vosko.h"
#include  "paw_my_constants.h"


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

static double onesixth,onethird;
static double rs_scale;

static double *Vc;
static double *ec_functional;

/****************************************
 Function name	  : paw_init_vosko()
 Description	    :
****************************************/
void paw_init_vosko()
{

    Vc		       = paw_alloc_LogGrid();
    ec_functional = paw_alloc_LogGrid();

    /* define constants */
    onesixth  = 1.0/6.0;
    onethird  = 1.0/3.0;
    rs_scale = pow( (0.75/PI), (1.0/3.0));


}
/****************************************
 Function name	  : get_corr_pot_LDA()
****************************************/
double* paw_get_corr_pot_LDA()
{
    return Vc;
}
/****************************************
 Function name	  : paw_generate_corr_pot
 Description	    : this routine calculates the spin
                    polarized Vosko et. al. correlation potential.
                    This is a Ceperly and Alder parameterization
 Return type		  : void
 Argument         : double **rho
 Argument         : double **Vc
 Argument         : double *Ec
 Author     		  : Eric Bylaska & Marat Valiev
 Date & Time		  : 1/8/99 1:19:44 PM
****************************************/
#define small_number    1.0e-80
void paw_generate_corr_pot_LDA(double *rho)
{
    int	i;
    double rs,n;
    double x,xxp,dxxp;
    double ec_p;
    double uc_p;

    /* loggrid variables */
    int	   Ngrid;


    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();


    for (i=0; i<= Ngrid-1; i++)
    {
        n     = rho[i]/(4.0*PI) + small_number;
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

        Vc[i] = uc_p;


    } /*for i*/



} /* R_Vosko */

/****************************************
 Function name	  : paw_get_correlation_energy_LDA(double *rho)
****************************************/
double paw_get_correlation_energy_LDA(double *rho)
{
    int	i;
    double rs,n;
    double x,xxp;
    double ec_p;
    double Ec;
    double *tmp;

    /* loggrid variables */
    int	   Ngrid;


    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();
    tmp		    = paw_scratch_LogGrid();

    /* allocate temporary memory */
    tmp		        = paw_alloc_LogGrid();
    ec_functional = paw_alloc_LogGrid();


    for (i=0; i<= Ngrid-1; i++)
    {
        n     = rho[i]/(4.0*PI) + small_number;
        rs    = rs_scale/pow(n,onethird);

        x     = sqrt(rs);

        xxp  = rs + bp*x + cp;
        ec_p = cp1*log(rs/xxp) + cp2*log( (x+cp3)*(x+cp3)/xxp)
               + cp4*atan(cp5/(x+cp6));

        ec_functional[i] = ec_p;


    } /*for i*/


    for (i=0; i<= Ngrid-1; i++)
        tmp[i] = rho[i]*ec_functional[i];

    Ec = paw_Integrate_LogGrid(tmp);

    return Ec;

}

