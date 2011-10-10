#include <stdio.h>
#include <math.h>
#include "loggrid.h"
#include "lyp_correlation.h"


/**R_LYP_Correlation***/

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
    double F1,F12,dF1,ddF1;
    double G1,dG1,ddG1,ddG1A,ddG1B,ddG1C;
    double T1,T2,T3,T4;
    double tw;
    double n;
    double agr,agr2,lap;
    double rho_onethird, rho_twothirds,rho_fourthirds,rho_fivethirds,rho_seventhirds,rho_eightthirds,rho3;
    double rho_eleventhirds,rho4,rho_thirteenthirds,e_c_rho,rho_fivethirds2;
    double ce,cp;
    double norm;

    /* loggrid variables */
    int Ngrid;
    double *rgrid;

    /* temporary local grids */
    double *ce_density;
    double *tmp;
    double *drho;
    double *ddrho;
    double *dadrho;

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
    dadrho     = alloc_LogGrid();

    /* calculate drho, ddrho */
    for (i=0;i<Ngrid;++i)
        tmp[i] = rho[i]/(4.0*pi);;
    Derivative_LogGrid(tmp,drho);
    Derivative_LogGrid(drho,ddrho);

    /* calculate dadrho */
    for (i=0;i<Ngrid;++i)
        tmp[i] = fabs(drho[i]);
    Derivative_LogGrid(tmp,dadrho);

    /* main loop */

    for (i=0;i<Ngrid;++i)
    {
        n = rho[i]/(4.0*pi);
        if ( n > 1.2e-18 )
        {
            /* define agr, agr2, lap */
            agr = fabs(drho[i]);
            agr2 = agr*agr;
            lap = ddrho[i] + (2.0/rgrid[i])*drho[i];

            rho_onethird = pow(n,-1.0/3.0);
            rho_twothirds = pow(n,-2.0/3.0);
            rho_fourthirds = pow(n,-4.0/3.0);
            rho_fivethirds = pow(n,-5.0/3.0);
            rho_seventhirds = pow(n,-7.0/3.0);
            rho_eightthirds = pow(n,-8.0/3.0);
            rho3 = 1/(n*n*n);
            rho_eleventhirds = pow(n,-11.0/3.0);
            rho4 = 1/(n*n*n*n);
            rho_thirteenthirds = pow(n,-13.0/3.0);

            rho_fivethirds2 = pow(n,5.0/3.0);



            e_c_rho = exp(-1.0*c*rho_onethird);






            F1 = 1.0/(1.0+d*rho_onethird);
            F12 = F1*F1;
            dF1 = F12*rho_fourthirds*d/3.0;
            ddF1 = (d/3.0)*(2.0*F1*dF1*rho_fourthirds - (4.0/3.0)*F12*rho_seventhirds);


            G1 = F1*rho_fivethirds*e_c_rho;
            dG1 = dF1*rho_fivethirds*e_c_rho - (5.0/3.0)*F1*rho_eightthirds*e_c_rho + (c/3.0)*F1*rho3*e_c_rho;              ddG1A = ddF1*rho_fivethirds*e_c_rho;
            ddG1B = dF1*((2.0*c/3.0)*rho3*e_c_rho - (10.0/3.0)*rho_eightthirds*e_c_rho);
            ddG1C = F1*((40.0/9.0)*rho_eleventhirds*e_c_rho - (14.0*c/9.0)*rho4*e_c_rho + (c*c/9.0)*rho_thirteenthirds*e_c_rho);
            ddG1 = ddG1A + ddG1B + ddG1C;

            /*debug set space derivatives to zero*/
            //agr = 0.0;
            //agr2 = 0.0;
            //lap = 0.0;

            T1 = -1.0*a*(dF1*n + F1);
            T2 = -1.0*a*b*Cf*rho_fivethirds2*(dG1*n + (8.0/3.0)*G1);
            T3 = -1.0*(a*b/4.0)*(ddG1*n*agr2 + dG1*(3.0*agr2 + 2*n*lap) + 4.0*G1*lap);
            T4 = -1.0*(a*b/72.0)*(3.0*ddG1*n*agr2 + dG1*(5.0*agr2 + 6.0*n*lap) + 4.0*G1*lap);

            cp = T1 + T2 + T3 + T4;

//             tw = (1.0/8.0)*((agr2/n) - lap);


//             ce = -1.0*(a*F1/n)*(n + b*rho_twothirds*(Cf*rho_fivethirds2-2.0*tw+(1.0/9.0)*tw+(1.0/18.0)*lap)*e_c_rho);
            tw = agr2/(8.0*n) - lap/8.0;
            ce = (-1.0*a/(1+d*rho_onethird))*(n + b*rho_twothirds*(Cf*rho_fivethirds2 - 2.0*tw + (1.0/9.0)*tw + (1.0/18.0)*lap)*e_c_rho);

            ce_density[i] = ce/n;
            Vc[i] = cp;


        }/* endif */
        else
        {
            ce_density[i] = 0.0;
            Vc[i] = 0.0;

        }/* end if else */
    }/*for i*/


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
    dealloc_LogGrid(dadrho);

//     printf("gerg\n");

} /* R_BLYP_Correlation */
/* $Id$ */
