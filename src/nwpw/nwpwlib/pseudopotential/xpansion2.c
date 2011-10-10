/* xpansion2.c
   Author - Eric Bylaska
*/
#include	<stdio.h>
#include	"loggrid.h"
#include	"pred_cor.h"
#include	"schrodin.h"
#include	"gaussj.h"
#include	"xpansion2.h"

#define	SMALL	1.0e-9


static	int	l;
static	int	match;
static	double	poly[10];
static	double	c[10];
static	double	rc[13];

static	double	ldpsi;
static	double	el;
static	double	*Vall;
static	double	*ul;
static	double	*ul_prime;
static	double	*Vl;
static	double	*wl;
static	double	*wl_prime;

void	p_xpansion2(double p[])
{
    int    i;
    double r1,r2,r4,r6,r8,r10;
    double *r;

    r = r_LogGrid();

    /* generate p(r) */
    for (i=0; i<=match; ++i)
    {
        r1  = r[i];
        r2  = r1*r1;
        r4  = r2*r2;
        r6  = r4*r2;
        r8  = r6*r2;
        r10 = r8*r2;
        p[i] = c[0]
               + c[1]*r2
               + c[2]*r4
               + c[3]*r6
               + c[4]*r8
               + c[5]*r10;

    }
} /* p_xpansion */



void	dp_xpansion2(double dp[])
{
    int    i;
    double r1,r2,r3,r5,r7,r9,r11;
    double *r;

    r = r_LogGrid();

    /* generate (d/dr)p(r) */
    for (i=0; i<=match; ++i)
    {
        r1  = r[i];
        r2  = r1*r1;
        r3  = r1*r2;
        r5  = r3*r2;
        r7  = r5*r2;
        r9  = r7*r2;
        r11 = r9*r2;
        dp[i] = 2.0 *c[1]*r1
                + 4.0 *c[2]*r3
                + 6.0 *c[3]*r5
                + 8.0 *c[4]*r7
                + 10.0*c[5]*r9;
    }
} /* dp_xpansion2 */



void	ddp_xpansion2(double ddp[])
{
    int    i;
    double r1,r2,r4,r6,r8,r10;
    double *r;

    r = r_LogGrid();

    /* generate (d/dr)(d/dr) p(r) */
    for (i=0; i<=match; ++i)
    {
        r1  = r[i];
        r2  = r1*r1;
        r4  = r2*r2;
        r6  = r4*r2;
        r8  = r6*r2;
        r10 = r8*r2;
        ddp[i] = 2.0  *c[1]
                 + 12.0 *c[2]*r2
                 + 30.0 *c[3]*r4
                 + 56.0 *c[4]*r6
                 + 90.0 *c[5]*r8;
    }
} /* ddp_xpansion2 */


void	dddp_xpansion2(double dddp[])
{
    int    i;
    double r1,r2,r3,r5,r7,r9;
    double *r;

    r = r_LogGrid();

    /* generate (d/dr)**3 p(r) */
    for (i=0; i<=match; ++i)
    {
        r1  = r[i];
        r2  = r1*r1;
        r3  = r1*r2;
        r5  = r3*r2;
        r7  = r5*r2;
        r9  = r7*r2;
        dddp[i] =  24.0 *c[2]*r1
                   + 120.0 *c[3]*r3
                   + 336.0 *c[4]*r5
                   + 720.0 *c[5]*r7;
    }
} /* dddp_xpansion2 */

void	ddddp_xpansion2(double ddddp[])
{
    int    i;
    double r1,r2,r4,r6,r8;
    double *r;

    r = r_LogGrid();

    /* generate (d/dr)**4 p(r) */
    for (i=0; i<=match; ++i)
    {
        r1  = r[i];
        r2  = r1*r1;
        r4  = r2*r2;
        r6  = r4*r2;
        r8  = r6*r2;
        ddddp[i] =  24.0 *c[2]
                    + 360.0  *c[3]*r2
                    + 1680.0 *c[4]*r4
                    + 5040.0 *c[5]*r6;
    }
} /* ddddp_xpansion2 */


/********************************
 *				*
 *	 chi_xpansion2		*
 *				*
 ********************************/

void	chi_xpansion2()
{
    int    i,Ngrid;
    double *r,*dp,*ddp;


    Ngrid = N_LogGrid();
    r     = r_LogGrid();

    dp   = alloc_LogGrid();
    ddp  = alloc_LogGrid();

    dp_xpansion2(dp);
    ddp_xpansion2(ddp);

    for (i=0; i<=match; ++i)
    {
        Vl[i] = el + dp[i]*(l+1.0)/r[i]
                + 0.5*ddp[i] + 0.5*dp[i]*dp[i];
    }

    for (i=(match+1); i<Ngrid; ++i)
        Vl[i] = Vall[i];

    /* make it a projector rather than a potential */
    for (i=0; i<Ngrid; ++i)
        Vl[i] = Vl[i]*wl[i];

    dealloc_LogGrid(dp);
    dealloc_LogGrid(ddp);
} /* chi_xpansion2 */





/********************************
 *				*
 *	 psi_xpansion2    	*
 *				*
 ********************************/

void	psi_xpansion2()
{
    int    i,Ngrid;
    double *r,*p;

    Ngrid = N_LogGrid();
    r     = r_LogGrid();

    p = alloc_LogGrid();
    p_xpansion2(p);

    for (i=0; i<=match; ++i)
        wl[i] = pow(r[i],l+1.0)*exp(p[i]);

    for (i=(match+1); i<Ngrid; ++i)
        wl[i] = ul[i];

    dealloc_LogGrid(p);

} /* psi_xpansion2 */


/********************************
 *				*
 *	 dpsi_xpansion2    	*
 *				*
 ********************************/

void	dpsi_xpansion2()
{
    int    i,Ngrid;
    double al;
    double *r,*p,*dp;

    Ngrid = N_LogGrid();
    r     = r_LogGrid();
    al    = log_amesh_LogGrid();

    p  = alloc_LogGrid();
    dp = alloc_LogGrid();
    p_xpansion2(p);
    dp_xpansion2(dp);

    /* wl_prime = (d/di)wl */
    for (i=0; i<=match; ++i)
        wl_prime[i] = (al*r[i])*( (l+1.0) + r[i]*dp[i] )*pow(r[i],l)*exp(p[i]);

    for (i=(match+1); i<Ngrid; ++i)
        wl_prime[i] = ul_prime[i];

    dealloc_LogGrid(p);
    dealloc_LogGrid(dp);

} /* dpsi_xpansion2 */




/********************************
 *				*
 *	  get_c0_c10_xpansion2  *
 *				*
 ********************************/

int get_c0_c10_xpansion2()
{
    int    i;
    double delta,delta_old,tolerance;

    double  b[5],a[25],aa[25];
    double sum;

    /* define matrix a */
    a[0] = 1.0; a[5]= 1.0*rc[2]; a[10]=  1.0*rc[6]; a[15] =  1.0*rc[8]; a[20]=   1.0*rc[10];
    a[1] = 0.0; a[6]= 2.0*rc[1]; a[11]=  6.0*rc[5]; a[16] =  8.0*rc[7]; a[21]=  10.0*rc[9];
    a[2] = 0.0; a[7]= 2.0;       a[12]= 30.0*rc[4]; a[17]=  56.0*rc[6]; a[22]=  90.0*rc[8];
    a[3] = 0.0; a[8]= 0.0;       a[13]=120.0*rc[3]; a[18]= 336.0*rc[5]; a[23]= 720.0*rc[7];
    a[4] = 0.0; a[9]= 0.0;       a[14]=360.0*rc[2]; a[19]=1680.0*rc[4]; a[24]=5040.0*rc[6];


    /* iterate abount c4=c[2] */
    delta = c[2] = 0.0;
    tolerance    = 1.0;
    while (tolerance > SMALL)
    {
        b[0] = poly[0] - delta* 1.0 * rc[4];
        b[1] = poly[1] - delta* 4.0 * rc[3];
        b[2] = poly[2] - delta* 12.0* rc[2];
        b[3] = poly[3] - delta* 24.0* rc[1];
        b[4] = poly[4] - delta* 24.0;

        for (i=0; i<25; ++i) aa[i] = a[i];
        gaussj(5,aa,1,b);


        c[0]  = b[0]; /* c0 coefficient */
        c[1]  = b[1]; /* c2 coefficient */
        c[3]  = b[2]; /* c6 coefficient */
        c[4]  = b[3]; /* c8 coefficient */
        c[5]  = b[4]; /* c10 coefficient */



        /* can't optimize c[2], so set to zero  and return */
        sum = fabs(b[0]+b[1]+b[2]+b[3]+b[4]);

        if (sum>100.0)
        {
            delta = -0.01;
            b[0] = poly[0] - delta* 1.0 * rc[4];
            b[1] = poly[1] - delta* 4.0 * rc[3];
            b[2] = poly[2] - delta* 12.0* rc[2];
            b[3] = poly[3] - delta* 24.0* rc[1];
            b[4] = poly[4] - delta* 24.0;
            for (i=0; i<25; ++i) aa[i] = a[i];
            gaussj(5,aa,1,b);


            c[0]  = b[0]; /* c0 coefficient */
            c[1]  = b[1]; /* c2 coefficient */
            c[2]  = delta;
            c[3]  = b[2]; /* c6 coefficient */
            c[4]  = b[3]; /* c8 coefficient */
            c[5]  = b[4]; /* c10 coefficient */
            return 0;
        }

        delta_old = delta;
        delta     = -c[1]*c[1]/(2.0*l+5.0);
        tolerance = fabs(delta-delta_old);



    }

    c[2] = delta; /* c4 coefficient */

    return 1;

} /* get_c0_c10_xpansion2 */


/********************************
 *				*
 *	  init_xpansion2	*
 *				*
 ********************************/

int init_xpansion2(int    l_in,
                   int    match_in,
                   double el_in,
                   double *Vall_in,
                   double *ul_in,
                   double *ul_prime_in,
                   double *Vl_in,
                   double *wl_in,
                   double *wl_prime_in)
{
    int    i;
    double Vall_match, dVall_match,ddVall_match;
    double al;
    double *r;

    r  = r_LogGrid();
    al = log_amesh_LogGrid();

    l          = l_in;
    match      = match_in;
    el         = el_in;
    Vall       = Vall_in;
    ul         = ul_in;
    ul_prime   = ul_prime_in;
    Vl         = Vl_in;
    wl         = wl_in;
    wl_prime   = wl_prime_in;

    /* fix ul_prime = (d/dr)ul, rather then (d/di)ul */
    ldpsi        = ul_prime[match]/(al*r[match]*ul[match]);

    Vall_match   = Vall[match];
    dVall_match  = (1.0/(al*r[match]))*Derivative7_4(match,Vall);
    ddVall_match = (-1.0/(r[match]*r[match]*al))   *Derivative7_4(match,Vall)
                   + (+1.0/(r[match]*r[match]*al*al))*Laplacian7_4(match,Vall);

    rc[0]	       = 1.0;
    rc[1]	       = r[match];
    for (i=2; i<13; ++i)
        rc[i] = rc[1]*rc[i-1];


    /**************************************************************/
    /* define p(rcl), p'(rcl), p''(rcl), p'''(rcl) and p''''(rcl) */
    /**************************************************************/
    poly[0] = log(ul[match]/pow(rc[1],(l+1.0)));
    poly[1] = ldpsi - (l+1.0)/rc[1];
    poly[2] = 2.0*Vall_match
              - 2.0*el
              - (2.0*(l+1.0)/rc[1])*poly[1]
              - poly[1]*poly[1];
    poly[3] = 2.0*dVall_match
              + (2.0*(l+1.0)/rc[2])*poly[1]
              - (2.0*(l+1.0)/rc[1])*poly[2]
              - 2.0*poly[1]*poly[2];
    poly[4] = 2.0*ddVall_match
              - (4.0*(l+1.0)/rc[3])*poly[1]
              + (4.0*(l+1.0)/rc[2])*poly[2]
              - (2.0*(l+1.0)/rc[1])*poly[3]
              - 2.0*poly[2]*poly[2]
              - 2.0*poly[1]*poly[3];


    /* generate  trouulier-Martin expansion, if possible */
    return get_c0_c10_xpansion2();

} /*init_xpansion2 */
/* $Id$ */
