
#include <stdio.h>
#include <math.h>
#include "loggrid.h"
#include "lyp_correlation.h"


/**R_BLYP_Correlation***/

/*calculates LYP correlation
  Entry: rho[]:  density
  Exit:  Vc[]: exchange potential
         Ec:   correlation energy
         Pc:   variational corrections for eigenvalues
*/

void  R_LYP_Correlation(rho,Vc,Ec,Pc)

double  rho[],
Vc[],
*Ec,
*Pc;
{
    /* local variables */
    int i;
    double pi;
    double a,b,c,d,Cf;
    double smallnumber;
    double F1,F12,F14,dF1,ddF1,dF1dn;
    double G1,dG1,ddG1;
    double n;
    double agr,agr2,lap;
    double e_c_rho;
    double ce,cp;


    /* loggrid variables */
    int Ngrid;
    double *rgrid;

    /* temporary local grids */
    double *ce_density;
    double *tmp;
    double *drho;
    double *ddrho;
    double *fn;
    double *fdn;
    double *dfdn;

    /* constants */
    pi = 4.0*atan(1.0);

    /* LYP parameters */
    a = 0.04918;
    b = 0.132;
    c = 0.2533;
    d = 0.349;

    /* more local variables */
    smallnumber = 1.0e-80;


    /* LDA parameter */
    Cf = (3.0/10.0)*pow((3.0*pi*pi),2.0/3.0);

    /* access the loggrid variables */
    Ngrid = N_LogGrid();
    rgrid = r_LogGrid();

    /* allocate temporary memeory */
    ce_density = alloc_LogGrid();
    tmp        = alloc_LogGrid();
    drho       = alloc_LogGrid();
    ddrho      = alloc_LogGrid();
    fn         = alloc_LogGrid();
    fdn        = alloc_LogGrid();
    dfdn        = alloc_LogGrid();
    /* calculate drho, ddrho */
    for (i=0;i<Ngrid;++i)
        tmp[i] = rho[i]/(4.0*pi);;
    Derivative_LogGrid(tmp,drho);
    Derivative_LogGrid(drho,ddrho);


    /* main loop */

    for (i=0;i<Ngrid;++i)
    {
        n = rho[i]/(4.0*pi);
        if ( n > 1.2e-18 )
        {
            /* define agr, agr2*/
            agr = fabs(drho[i]);
            agr2 = agr*agr;

            e_c_rho = exp(-1.0*c*pow(n,-1.0/3.0));

            F1 = 1.0/(1.0+d*pow(n,-1.0/3.0));

            F12 = F1*F1;

            dF1 = F12*pow(n,-4.0/3.0)*d/3.0;

            ddF1 = (d/3.0)*(2.0*F1*dF1*pow(n,-4.0/3.0) - (4.0/3.0)*F12*pow(n,-7.0/3.0));


            G1 = F1*pow(n,-5.0/3.0)*e_c_rho;

            dG1 = dF1*pow(n,-5.0/3.0)*e_c_rho - (5.0/3.0)*F1*pow(n,-8.0/3.0)*e_c_rho + (c/3.0)*F1*(1/(n*n*n))*e_c_rho;

            ddG1 = ddF1*pow(n,-5.0/3.0)*e_c_rho
                   + dF1*((2.0*c/3.0)*(1/(n*n*n))*e_c_rho
                          - (10.0/3.0)*pow(n,-8.0/3.0)*e_c_rho)
                   + F1*((40.0/9.0)*pow(n,-11.0/3.0)*e_c_rho
                         - (14.0*c/9.0)*(1/(n*n*n*n))*e_c_rho
                         + (c*c/9.0)*pow(n,-13.0/3.0)*e_c_rho);

            /****************************no lap in fn and fdn******/
            fn[i] = -1.0*a*(dF1*n + F1) - a*b*Cf*pow(n,5.0/3.0)*(dG1*n+(8.0/3.0)*G1) + (a*b/4.0)*(ddG1*n*agr2 + 3.0*dG1*agr2) + (a*b/72.0)*(3.0*ddG1*n*agr2 + 5.0*dG1*agr2);

            /************fdn is derivative of f w.r.t agr divided by agr*************/
            fdn[i] = (a*b/4.0)*(2.0*dG1*n + 4.0*G1) + (a*b/72.0)*(6.0*dG1*n + 4.0*G1);




            /***         expression for ce via integration by parts on the del square terms of correlation energy ****************/
            ce = -1.0*a*F1*n - a*b*Cf*G1*pow(n,8.0/3.0) + (a*b/4.0)*(dG1*n*agr2 + 2.0*G1*agr2) + (a*b/72.0)*(3.0*dG1*n*agr2 + 2.0*G1*agr2);


            ce_density[i] = ce/n;











        }/* endif */
        else
        {
            ce_density[i] = 0.0;

        }/* end if else */
    }/*for i*/


    /*calculate derivative fdn/agr--->dfdn*/
    Derivative_LogGrid(fdn,dfdn);

    /*calc correlation potential*******************************/
    for (i=0;i<Ngrid;++i)
    {
        n = rho[i]/(4.0*pi);
        if ( n > 1.2e-18 )
        {
            lap   = ddrho[i] + (2.0/rgrid[i])*drho[i];
            cp    = fn[i] - (dfdn[i]*drho[i] + fdn[i]*lap);
            Vc[i] = cp;
        }/* endif */
        else
        {
            Vc[i] = 0.0;
        }/* end if else*/
    }/*correlation potential*/


    /* calculate Ec and Pc */

    /* integrate rho*ec_density */
    for (i=0;i<Ngrid;++i)
        tmp[i] = (rho[i]) * ce_density[i];
    *Ec = Integrate_LogGrid(tmp);


    /* integrate rho*Vc */
    for (i=0;i<Ngrid;++i)
        tmp[i] = (rho[i])*Vc[i];
    *Pc = Integrate_LogGrid(tmp);




    /* deallocate temporary memory */
    dealloc_LogGrid(ce_density);
    dealloc_LogGrid(tmp);
    dealloc_LogGrid(drho);
    dealloc_LogGrid(ddrho);
    dealloc_LogGrid(fn);
    dealloc_LogGrid(fdn);
    dealloc_LogGrid(dfdn);


} /* R_BLYP_Correlation */

/* $Id$ */
