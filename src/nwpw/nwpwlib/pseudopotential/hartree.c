/*
 $Id$
   Hartree.c - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the hartree
   potential and energy of a density n = rho_up+rho_down, defined on
   a log grid.

*/
#include	<stdio.h>
#include        "loggrid.h"
#include        "pred_cor.h"
#include        "hartree.h"


/********************************
 *                              *
 *           R_Hartree          *
 *                              *
 ********************************/

/* returns the hartree potential, and
   hartree energy of a density, n=rho_up+rho_down
   defined on the log grid.

   Entry - n[]
           charge
   Exit  - Vh[],

*/

double  R_Hartree(n,charge,Vh)
double 	n[],
charge,
Vh[];
{
    int i,Ngrid;
    double log_amesh;
    double *rgrid,*tmp;
    double r,r2,r3,tt;
    double E;

    /* access the size of grid, and log(amesh) */
    Ngrid     = N_LogGrid();
    log_amesh = log_amesh_LogGrid();

    /* get access to rgrid, and a tmp grid */
    rgrid = r_LogGrid();
    tmp   = alloc_LogGrid();

    for (i=0; i<Ngrid; ++i)
    {
        r  = rgrid[i];
        r3 = r*r*r;
        tmp[i] = (log_amesh*r3)*n[i];
    }

    Zero_LogGrid(Vh);

    /* define boundry at r->infinity */
    Vh[Ngrid-1] = charge; /* set the boundry condition */
    Vh[Ngrid-2] = charge;
    Vh[Ngrid-3] = charge;
    /*Vh[Ngrid-4] = charge; */

    /* Integrate tmp[i] to zero */
    for (i=(Ngrid-3); i>0; --i)
    {
        Vh[i-1] = Vh[i] + Corrector_In_F(i,tmp);
    }

    for (i=0; i<Ngrid; ++i)
    {
        r  = rgrid[i];
        r2 = r*r;
        tmp[i] = (log_amesh*r2)*n[i];
    }

    /* integrate tmp to zero */
    tt = 0.0;
    for (i=(Ngrid-3); i>0; --i)
    {
        tt        = tt + Corrector_In_F(i,tmp);
        Vh[i-1]   = Vh[i-1] - rgrid[i-1]*tt;
    }

    for (i=0; i<Ngrid; ++i)
        Vh[i] = Vh[i]/rgrid[i];


    /* calculate Eh */
    /* note that the integration is weird, because */
    /* we are integrating from 0 to infinity, and  */
    /* our log grid goes from r0 to 45.0           */
    for (i=0; i<Ngrid; ++i)
    {
        /*
         r  = rgrid[i];
         r3 = r*r*r;
         tmp[i] = r3*(n[i]*Vh[i]);
         */
        tmp[i] = n[i]*Vh[i];
    }
    /*
    E = (9.0*tmp[0] + 23.0*tmp[1] + 28.0*tmp[2])/28.0;
    for (i=3; i<Ngrid; ++i)
       E += tmp[i];
    E = log_amesh*E + tmp[0]/3.0;
    */
    E = Integrate_LogGrid(tmp);

    /* dealloc tmp memory */
    dealloc_LogGrid(tmp);

    return E;

} /* R_Hartree */

