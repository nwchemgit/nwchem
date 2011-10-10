/* generate_rho_semicore.c
   Author - Eric Bylaska
*/
#include	<stdio.h>
#include	"loggrid.h"
#include	"generate_rho_semicore.h"
#include        "gaussj.h"


/********************************
 *				*
 *   generate_rho_semicore      *
 *				*
 ********************************/

/* this routine calculates  semicore density

   Entry - rho_core: the core density
         - r_semicore: the matching point
   Exit  - rho_semicore: the semicore density
*/

void  generate_rho_semicore(semicore_type,rho_core,r_semicore,rho_semicore)

int semicore_type;
double	rho_core[],
r_semicore,
rho_semicore[];
{

    int    k,index;
    int	   Ngrid;
    double *rgrid;
    double *drho_core,*ddrho_core;

    double c1,c2;
    double x1,x2,cs,ct,A,B,f,fp;
    double r,rho,drho,ddrho,dddrho;
    double c0,c3,c4,c5,c6;


    double bs[4],as[4*4];


    /* access the loggrid variables */
    Ngrid     = N_LogGrid();
    rgrid     = r_LogGrid();

    /* allocate temporary memory */
    drho_core  = alloc_LogGrid();
    ddrho_core = alloc_LogGrid();


    index = index_r_LogGrid(r_semicore);
    Derivative_LogGrid(rho_core,drho_core);
    Derivative_LogGrid(drho_core,ddrho_core);

    r     = rgrid[index];
    rho   = rho_core[index];
    drho  = drho_core[index];
    ddrho = ddrho_core[index];

    Derivative_LogGrid(ddrho_core,drho_core);
    dddrho = drho_core[index];

    /**************************/
    /*        r < rcore       */
    /**************************/

    /*********************************************************/
    /* find quadratic coefficients of rho_semicore expansion */
    /*********************************************************/
    if (semicore_type==0)
    {
        c1 = rho - (r*drho)/2.0;
        c2 = drho/(2.0*r);

        for (k=0; k<index; ++k)
            rho_semicore[k] = c1 + c2*rgrid[k]*rgrid[k];
    }


    /***********************************************************/
    /* find A and B coefficients of rho_semicore expansion     */
    /* S. Louie's semicore formula rho[r<rcore] = A*sin(B*r)/r */
    /***********************************************************/
    if (semicore_type==1)
    {
        x2 = 0.5*sqrt(-3.0*(r*drho/rho));
        x1 = x2-1.0;
        while (fabs(x2-x1) > 1.0e-6)
        {
            x1 = x2;
            ct = cos(x1)/sin(x1);
            cs = sin(x1)*sin(x1);
            cs = 1.0/cs;
            f  = x1*ct - (r*drho/rho+1.0);
            fp = ct - x1*cs;
            x2 = x1 - f/fp;
        }
        B = x2/r;
        A = (r*rho)/sin(x2);

        for (k=0; k<index; ++k)
            rho_semicore[k] = (A/rgrid[k])*sin(B*rgrid[k]);
    }

    /*************************************************************/
    /* find coefficients of rho_semicore expansion               */
    /* Fuchs and Scheffler's semicore formula:                   */
    /*   rho[r<rcore] = c0 + c3*r^3 + c4*r^4 + c5*r^5 + c6*r^6   */
    /*************************************************************/
    if (semicore_type==2)
    {
        /*
        printf(" + Utilizing Fuchs and Schefflers semicore formula\n");
        printf(" + ( rho[r<rcore] = c0 + c3*r^3 + c4*r^4 + c5*r^5 + c6*r^6 )\n");
        */

        c0 = rho - (r*drho)/2.0;

        bs[0] = rho-c0;
        bs[1] = drho;
        bs[2] = ddrho;
        bs[3] = dddrho;

        as[0] = 1.0*pow(r,3.0); as[4] = 1.0*pow(r,4.0); as[8]  =  1.0*pow(r,5.0);
        as[1] = 3.0*r*r;        as[5] = 4.0*pow(r,3.0); as[9]  =  5.0*pow(r,4.0);
        as[2] = 6.0*r;          as[6] = 12.0*r*r;       as[10] = 20.0*pow(r,3.0);
        as[3] = 6.0;            as[7] = 24.0*r;         as[11] = 60.0*r*r;

        as[12] =   1.0*pow(r,6.0);
        as[13] =   6.0*pow(r,5.0);
        as[14] =  30.0*pow(r,4.0);
        as[15] = 120.0*pow(r,3.0);

        gaussj(4,as,1,bs);
        c3 = bs[0]; c4 = bs[1]; c5 = bs[2]; c6 = bs[3];

        for (k=0; k<index; ++k)
            rho_semicore[k] = c0
                              + c3*rgrid[k]*rgrid[k]*rgrid[k]
                              + c4*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k]
                              + c5*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k]
                              + c6*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k];
    }

    /**************************/
    /*        r > rcore       */
    /**************************/
    for (k=index; k<Ngrid; ++k)
        rho_semicore[k] = rho_core[k];

    dealloc_LogGrid(drho_core);
    dealloc_LogGrid(ddrho_core);


} /* generate_rho_semicore */
/* $Id$ */
