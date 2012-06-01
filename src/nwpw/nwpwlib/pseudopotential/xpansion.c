/*
 $Id$
   xpansion.c
   Author - Eric Bylaska
*/
#include        <stdlib.h>
#include	<stdio.h>
#include	"loggrid.h"
#include	"pred_cor.h"
#include	"schrodin.h"
#include	"gaussj.h"
#include	"xpansion.h"

#define	SMALL	1.0e-9


static	int	l;
static	int	match;
static	double  occupation;
static	double	poly[10];
static	double	c[10];
static	double	rc[13];

static	double	ae_core_charge;
static	double	ldpsi;
static	double	el;
static	double	*Vall;
static	double	*ul;
static	double	*ul_prime;
static	double	*Vl;
static	double	*wl;
static	double	*wl_prime;

void	p_xpansion(double p[])
{
    int    i;
    double r1,r2,r4,r6,r8,r10,r12;
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
        r12 = r10*r2;
        p[i] = c[0]
               + c[1]*r2
               + c[2]*r4
               + c[3]*r6
               + c[4]*r8
               + c[5]*r10
               + c[6]*r12;
    }
} /* p_xpansion */



void	dp_xpansion(double dp[])
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
                + 10.0*c[5]*r9
                + 12.0*c[6]*r11;
    }
} /* dp_xpansion */



void	ddp_xpansion(double ddp[])
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
                 + 90.0 *c[5]*r8
                 + 132.0*c[6]*r10;
    }
} /* ddp_xpansion */


void	dddp_xpansion(double dddp[])
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
                   + 720.0 *c[5]*r7
                   + 1320.0*c[6]*r9;
    }
} /* dddp_xpansion */

void	ddddp_xpansion(double ddddp[])
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
                    + 5040.0 *c[5]*r6
                    + 11880.0*c[6]*r8;
    }
} /* ddddp_xpansion */


/********************************
 *				*
 *	 psp_xpansion		*
 *				*
 ********************************/

void	psp_xpansion()
{
    int    i,Ngrid;
    double *r,*dp,*ddp;


    Ngrid = N_LogGrid();
    r     = r_LogGrid();

    dp   = alloc_LogGrid();
    ddp  = alloc_LogGrid();

    dp_xpansion(dp);
    ddp_xpansion(ddp);

    for (i=0; i<=match; ++i)
    {
        Vl[i] = el + dp[i]*(l+1.0)/r[i]
                + 0.5*ddp[i] + 0.5*dp[i]*dp[i];
    }

    for (i=(match+1); i<Ngrid; ++i)
        Vl[i] = Vall[i];

    dealloc_LogGrid(dp);
    dealloc_LogGrid(ddp);
} /* psp_xpansion */



/********************************
 *				*
 *	 psi_xpansion    	*
 *				*
 ********************************/

void	psi_xpansion()
{
    int    i,Ngrid;
    double *r,*p;

    Ngrid = N_LogGrid();
    r     = r_LogGrid();

    p = alloc_LogGrid();
    p_xpansion(p);

    for (i=0; i<=match; ++i)
        wl[i] = pow(r[i],l+1.0)*exp(p[i]);

    for (i=(match+1); i<Ngrid; ++i)
        wl[i] = ul[i];

    dealloc_LogGrid(p);

} /* psi_xpansion */


/********************************
 *				*
 *	 dpsi_xpansion    	*
 *				*
 ********************************/

void	dpsi_xpansion()
{
    int    i,Ngrid;
    double al;
    double *r,*p,*dp;

    Ngrid = N_LogGrid();
    r     = r_LogGrid();
    al    = log_amesh_LogGrid();

    p  = alloc_LogGrid();
    dp = alloc_LogGrid();
    p_xpansion(p);
    dp_xpansion(dp);

    /* wl_prime = (d/di)wl */
    for (i=0; i<=match; ++i)
        wl_prime[i] = (al*r[i])*( (l+1.0) + r[i]*dp[i] )*pow(r[i],l)*exp(p[i]);

    for (i=(match+1); i<Ngrid; ++i)
        wl_prime[i] = ul_prime[i];

    dealloc_LogGrid(p);
    dealloc_LogGrid(dp);

} /* dpsi_xpansion */




/********************************
 *				*
 *	  get_c0_c10		*
 *				*
 ********************************/

void	get_c0_c10()
{
    int    i,j,k,iteration;
    int    indx_rc,factor;
    double delta;
    double ddelta;
    double tolerance, fdnew,fdold;
    double  psp_core_charge;

    double  b[5],a[25];




    c[0]  = delta = 0.0;  /* c0 coefficient */
    c[6]  = 0.0;          /* c12 coefficient */
    fdold     = 0.0;
    fdnew     = 0.0;
    tolerance = 1.0;
    iteration = 0;
    while ((tolerance > SMALL) && ( iteration <= 50))
    {
        ++iteration;
        b[0] = poly[0] - delta;
        b[1] = poly[1];
        b[2] = poly[2];
        b[3] = poly[3];
        b[4] = poly[4];
        /* set up a matrix */
        for (i=0; i<5; ++i)
            for (j=0; j<5; ++j)
            {
                indx_rc = 2*(j+1)-i;
                factor  = 1;
                for (k=2*(j+1); k>(indx_rc); --k)
                    factor *= k;
                if (indx_rc >= 0)
                    a[i+j*5] = factor*rc[indx_rc];
                else
                    a[i+j*5] = 0.0;
            }

        gaussj(5,a,1,b);


        c[1]  = b[0]; /* c2 coefficient */
        c[2]  = b[1]; /* c4 coefficient */
        c[3]  = b[2]; /* c6 coefficient */
        c[4]  = b[3]; /* c8 coefficient */
        c[5]  = b[4]; /* c10 coefficient */

        /* generate pseudowavefunction neglecting exp(delta)=exp(c0) */
        psi_xpansion();

        /* Calculate the psp core charge */
        psp_core_charge = Norm_LogGrid(match,l+1.0,wl);
        if (psp_core_charge<0.0) psp_core_charge=fabs(psp_core_charge); /*debug*/

        /* calculate new delta */
        fdold = fdnew;
        fdnew = log(ae_core_charge/psp_core_charge) - 2.0*delta;
        tolerance = fabs(fdnew);
        if (iteration == 1)
            ddelta = -0.5;
        else
            ddelta = -fdnew*ddelta/(fdnew-fdold);
        delta = delta + ddelta;


    } /* while find delta */
    c[0] = delta;


    /* error - too many iterations */
    if (iteration > 50)
    {
        printf("Error in get_c0_c10: delta loop \n");
        exit(98);
    }

} /* get_c0_c10 */



/********************************
 *				*
 *	  get_c0_c12		*
 *				*
 ********************************/
/*

   entry - gamma = c[1]=c2,
   exit - returns the constraint value, and c vector
*/

double	get_c0_c12(double gamma)
{
    int    i,j,k,iteration;
    int    indx_rc,factor;
    double delta;
    double ddelta;
    double tolerance, fdnew,fdold;
    double  psp_core_charge;
    double constraint;

    double  b[5],a[25];


    c[0]  = delta = 0.0;    /* c0 coefficient */
    c[1]  = gamma;          /* c2 coefficient */
    fdold = 0.0;
    fdnew = 0.0;
    tolerance = 1.0;
    iteration = 0;
    while ((tolerance > SMALL) && ( iteration <= 50))
    {
        ++iteration;
        b[0] = poly[0] - 1.0*gamma*rc[2] - delta;
        b[1] = poly[1] - 2.0*gamma*rc[1];
        b[2] = poly[2] - 2.0*gamma;
        b[3] = poly[3];
        b[4] = poly[4];
        /* set up a matrix */
        for (i=0; i<5; ++i)
            for (j=0; j<5; ++j)
            {
                indx_rc = 2*(j+2)-i;
                factor  = 1;
                for (k=2*(j+2); k>(indx_rc); --k)
                    factor *= k;
                if (indx_rc >= 0)
                    a[i+j*5] = factor*rc[indx_rc];
                else
                    a[i+j*5] = 0.0;
            }

        gaussj(5,a,1,b);

        c[2] = b[0];  /* c4 coefficient */
        c[3] = b[1];  /* c6 coefficient */
        c[4] = b[2];  /* c8 coefficient */
        c[5] = b[3];  /* c10 coefficient */
        c[6] = b[4];  /* c12 coefficient */

        /* generate pseudowavefunction neglecting exp(delta)=exp(c0) */
        psi_xpansion();

        /* Calculate the psp core charge */
        psp_core_charge = Norm_LogGrid(match,l+1.0,wl);
        if (psp_core_charge<0.0) psp_core_charge=fabs(psp_core_charge);

        /* calculate new delta */
        fdold     = fdnew;
        fdnew     = log(ae_core_charge/psp_core_charge) - 2.0*delta;
        tolerance = fabs(fdnew);
        if (iteration == 1)
            ddelta = -0.5;
        else
            ddelta = -fdnew*ddelta/(fdnew-fdold);
        delta = delta + ddelta;


    } /* while find delta */
    c[0] = delta;


    /* error - too many iterations */
    if (iteration > 200)
    {
        printf("Error in get_c0_c12: delta loop \n");
        exit(98);
    }

    /* constraint = c2^2 + (2*l+5)*c4 = 0 */
    constraint = 8.0*(c[1]*c[1] + (2.0*l+5.0)*c[2]);
    return constraint;

} /* get_c0_c12 */


/********************************
 *				*
 *	  init_xpansion		*
 *				*
 ********************************/

void init_xpansion(int    l_in,
                   int    match_in,
                   double occupation_in,
                   double el_in,
                   double *Vall_in,
                   double *ul_in,
                   double *ul_prime_in,
                   double *Vl_in,
                   double *wl_in,
                   double *wl_prime_in)
{
    int    i,iteration;
    double gamma,gamma_mid,dgamma;
    double constraint1,constraint2,constraint_mid;
    double gamma1,gamma2;
    double Vall_match, dVall_match,ddVall_match;
    double al;
    double *r;

    r  = r_LogGrid();
    al = log_amesh_LogGrid();

    l          = l_in;
    match      = match_in;
    occupation = occupation_in;
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

    /******************************************/
    /* Calculate the all-electron core charge */
    /******************************************/
    ae_core_charge = Norm_LogGrid(match,(l+1.0),ul);


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


    /* get initial guess for gamma */
    get_c0_c10();


    /* Bracket gamma, so that constraint==0 can be found */
    gamma1 =  c[1];
    gamma2 = -c[1];

    constraint1 = get_c0_c12(gamma1);

    constraint2 = get_c0_c12(gamma2);

    iteration = 0;
    while ((constraint1*constraint2 > 0.0) && (iteration <=50))
    {
        ++iteration;
        if (fabs(constraint1) < fabs(constraint2))
        {
            gamma1 = gamma1 + 1.6*(gamma1-gamma2);
            constraint1 = get_c0_c12(gamma1);
        }
        else
        {
            gamma2 = gamma2 + 1.6*(gamma2-gamma1);
            constraint2 = get_c0_c12(gamma2);
        }
    }


    /* perform bisection of gamma1 and gamma2 until constraint==0 */
    constraint1 = get_c0_c12(gamma1);
    constraint2 = get_c0_c12(gamma2);
    if (constraint1 < 0.0)
    {
        gamma  = gamma1;
        dgamma = (gamma2-gamma1);
    }
    else
    {
        gamma  = gamma2;
        dgamma = (gamma1-gamma2);
    }

    iteration      = 0;
    constraint_mid = 1.0;
    while ((fabs(dgamma) > SMALL) && (constraint_mid != 0.0) && (iteration<=80))
    {
        ++iteration;

        dgamma    = 0.5*dgamma;
        gamma_mid = gamma+dgamma;

        constraint_mid = get_c0_c12(gamma_mid);
        if (constraint_mid < 0.0)
            gamma = gamma_mid;
    }

} /*init_xpansion */

