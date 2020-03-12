#include <stdio.h>
#include <math.h>
#include "loggrid.h"
#include "becke_exchange.h"

/****restricted blyp exchange*******
entry: - rho[]: density
exit:  - Vx[]: exchange potential
       - Ex:   exchange energy
       - Px[]: variational exchange corrections for eigenvalues
*********************************/

#define tolrho 	2.0e-8
#define minden 	1.0e-10
#define minagr 	1.0e-10
#define minchi 	1.0e-6
#define maxchi 	10.0

void R_Becke_Exchange(rho,Vx,Ex,Px)

double rho[],
Vx[],
*Ex,
*Px;
{
    /*local variables*/
    int i;
    double lda_c;
    double beta;
    double pi;
    double X1, X2, X3, ex;
    double fdn_lda;
    double ux;
    double c;
    double norm;
    double lap,damp,ddamp;
    double z5,z4,z3,z2,z1;

    double rho_fourthirds,rho_onethird, two_onethird;


    int Ngrid;

    /* temporary local grids */
    double *rhoNRM;
    double *drho;
    double *ddrho;
    double *ex_density;
    double *tmp;
    double *rgrid;
    double *agr;
    double *chi;
    double *chidn;
    double *chiddn;
    double *H;
    double *F;
    double *Fdchi;
    double *G;
    double *G2;
    double *Gdr;
    double *G2dr;
    double *fdn;
    double *fddn;



    pi = 4.0*atan(1.0);
    c  = pow(2.0,1.0/3.0);


    /*Becke(beta) and lda (lda_c)  exchange parameter*/
    beta        = 0.0042;
    lda_c       = (3.0/2.0)*pow(3.0/(4.0*pi),1.0/3.0);

    /*access the loggrid variables*/
    Ngrid       = N_LogGrid();
    rgrid       = r_LogGrid();

    /*allocate temporary memory*/
    ex_density  = alloc_LogGrid();
    tmp         = alloc_LogGrid();
    drho        = alloc_LogGrid();
    ddrho        = alloc_LogGrid();
    rhoNRM      = alloc_LogGrid();
    agr         = alloc_LogGrid();
    chi         = alloc_LogGrid();
    chidn       = alloc_LogGrid();
    chiddn      = alloc_LogGrid();
    H           = alloc_LogGrid();
    F           = alloc_LogGrid();
    Fdchi       = alloc_LogGrid();
    G           = alloc_LogGrid();
    G2           = alloc_LogGrid();
    Gdr         = alloc_LogGrid();
    G2dr         = alloc_LogGrid();
    fdn         = alloc_LogGrid();
    fddn        = alloc_LogGrid();




    /*calculate derivatives******************************************/

    for (i=0; i<Ngrid; ++i)
       rhoNRM[i] = rho[i]/(4.0*pi) + minden;        /* normalize density     */


    Derivative_LogGrid(rhoNRM,drho);           /* drho                  */

    for (i=0;i<Ngrid;++i)
       agr[i] = fabs(drho[i]);             /* agr                   */

    Derivative_LogGrid(drho,ddrho);



    /* calculate chi chidn, chiddn ****************************/
    for (i=0;i<Ngrid;++i) {
       rho_fourthirds = pow(rhoNRM[i],4.0/3.0);
       chi[i] =  agr[i]/rho_fourthirds;

       chidn[i] = (-4.0/3.0)*chi[i]/rhoNRM[i];
       chiddn[i] = 1.0/rho_fourthirds;
    }


    /* calculate H         ******************************************/
    for (i=0;i<Ngrid;++i) {
       rho_fourthirds = pow(rhoNRM[i],4.0/3.0);
       H[i] = 1.0*beta*pow(2.0,1.0/3.0)*rho_fourthirds;
    }


    /* calculate F  and dF/dchi *************************************/
    for (i=0;i<Ngrid;++i) {
       if ((rhoNRM[i]>tolrho) && (agr[i]>minagr)) {
          F[i] = chi[i]*chi[i]/(1.0+6.0*beta*c*chi[i]*log(c*chi[i] + sqrt(1.0+c*c*chi[i]*chi[i])));

          Fdchi[i] = 2*F[i]/chi[i] - (F[i]*F[i]/(chi[i]*chi[i]))
                      *(6.0*beta*c*log(c*chi[i] + sqrt(1.0+c*c*chi[i]*chi[i]))
                       +6.0*beta*c*c*chi[i]/sqrt(1.0+c*c*chi[i]*chi[i]));
       }
    }

    /* calculate fdn and fddn******************************************/
    for (i=0;i<Ngrid;++i) {
       rho_onethird = pow(rhoNRM[i],1.0/3.0);
       fdn_lda = -lda_c/c*(4.0/3.0)*rho_onethird;
       fdn[i] = fdn_lda - (4.0/3.0)*H[i]*F[i]/rhoNRM[i] - H[i]*Fdchi[i]*chidn[i];
       fddn[i] = -1.0*H[i]*Fdchi[i]*chiddn[i];
    }


    /* calculate G ***************************************************/
    for (i=0;i<Ngrid;++i)
       if (agr[i] > minagr) 
          G[i] = fddn[i]/agr[i];
       else
          G[i] = 0.0;
    /* calculate dG/dr ***********************************************/

    Derivative_LogGrid(G,Gdr);


    /* calculate exchange potential and exchange energy density*******/
    for (i=0;i<Ngrid;++i) {
       /*exchange potential ****************************************/
       lap = ddrho[i] + (2.0/rgrid[i])*drho[i];
       ux  = fdn[i] - (Gdr[i]*drho[i] + G[i]*lap);

       /* exchange energy density ex ****/
       rho_onethird = pow(rhoNRM[i],1.0/3.0);
       X1 = -1.0*(lda_c/c)*rho_onethird;
       X2 = c*beta*rho_onethird*chi[i]*chi[i];
       X3 = 1.0 + 6.0*beta*c*chi[i]*log(c*chi[i] + sqrt(1.0+c*c*chi[i]*chi[i]));
       ex = X1 - X2/X3;

       Vx[i]         = ux;
       ex_density[i] = ex;            /*energy density*/

       if ((fdn[i]-ex)>0.0)
       {
            rho_onethird = pow(rhoNRM[i],1.0/3.0);
            ex_density[i] = -1.0*(lda_c/c)*rho_onethird;
            Vx[i]         = -lda_c/c*(4.0/3.0)*rho_onethird;
       }
    }/*for i*/


//    printf("------------------------------------\n");
//    for (i=0;i<Ngrid;++i)
//       printf("%le %le %le\n",rho[i],Vx[i],ex_density[i]);
//    printf("------------------------------------\n");



    /*calculate Ex, and Px */

    /*integrate rho*ex_density***/
    for (i=0;i<Ngrid;++i)
        tmp[i] = rho[i] * ex_density[i];
    *Ex = Integrate_LogGrid(tmp);

    /*integrate rho*Vx */
    for (i=0;i<Ngrid;++i)
        tmp[i] = rho[i] * Vx[i];
    *Px = Integrate_LogGrid(tmp);



    /*deallocate temporary memory */
    dealloc_LogGrid(ex_density);
    dealloc_LogGrid(tmp);
    dealloc_LogGrid(drho);
    dealloc_LogGrid(rhoNRM);
    dealloc_LogGrid(agr);
    dealloc_LogGrid(chi);
    dealloc_LogGrid(chidn);
    dealloc_LogGrid(chiddn);
    dealloc_LogGrid(H);
    dealloc_LogGrid(F);
    dealloc_LogGrid(Fdchi);
    dealloc_LogGrid(G);
    dealloc_LogGrid(Gdr);
    dealloc_LogGrid(fdn);
    dealloc_LogGrid(fddn);




}/* R_BLYP_Exchange */
/* $Id$ */
