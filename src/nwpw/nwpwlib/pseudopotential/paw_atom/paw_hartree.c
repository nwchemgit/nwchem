/*
   $Id$
*/

#include	<stdio.h>
#include  <math.h>
#include  "paw_loggrid.h"
#include  "paw_pred_cor.h"
#include  "paw_hartree.h"


static double *Vh;

/****************************************
 Function name	  : paw_init_hartree(double Z)
 Description	    :
****************************************/
void paw_init_hartree()
{

    Vh  = paw_alloc_LogGrid();

}
/****************************************
 Function name	  : paw_generate_hartree_pot
 Description	    :
 Return type		  : void
 Argument         : double *n
 Argument         : double *Vh
 Author     		  : Marat Valiev (modified from Erics code)
 Date & Time		  : 1/11/99 3:49:01 PM
****************************************/
void paw_generate_hartree_pot(double *n)
{

    int i;
    int Ngrid;
    double* r;
    double* r2;
    double* r3;
    double tt;
    double charge;
    double log_amesh;
    double *tmp;

    /* get access to rgrid, and a tmp grid */
    Ngrid     = paw_N_LogGrid();
    log_amesh = paw_log_amesh_LogGrid();
    r         = paw_r_LogGrid();
    r2        = paw_r2_LogGrid();
    r3        = paw_r3_LogGrid();
    tmp       = paw_scratch_LogGrid();

    paw_Zero_LogGrid(Vh);

    charge   = paw_Integrate_LogGrid(n);

    for (i=0; i <= Ngrid-1; i++)
        tmp[i] = (log_amesh*r3[i])*n[i];


    /* define boundry at r->infinity */

    Vh[Ngrid-1] = charge;
    Vh[Ngrid-2] = charge;
    Vh[Ngrid-3] = charge;

    /* Integrate tmp[i] to zero */
    for (i=(Ngrid-3); i>=1; i--)
        Vh[i-1] = Vh[i] + paw_Corrector_In_F(i,tmp);

    for (i=0; i<=Ngrid-1; i++)
        tmp[i] = (log_amesh*r2[i])*n[i];

    /* integrate tmp to zero */
    tt = 0.0;
    for (i=(Ngrid-3); i>=1; i--)
    {
        tt        = tt + paw_Corrector_In_F(i,tmp);
        Vh[i-1]   = Vh[i-1] - r[i-1]*tt;
    }

    for (i=0; i<=Ngrid-1; i++)
        Vh[i] = Vh[i]/r[i];

}


/****************************************
 Function name	  : paw_get_hartree_energy(double *n)
 Description	    :
****************************************/
double paw_get_hartree_energy(double *n)
{

    int i;
    int Ngrid;
    double Eh;
    double *tmp;

    Ngrid     = paw_N_LogGrid();
    tmp       = paw_scratch_LogGrid();

    paw_generate_hartree_pot(n);

    for (i=0; i<Ngrid; ++i)
        tmp[i] = n[i]*Vh[i];


    Eh = 0.5*paw_Integrate_LogGrid(tmp);

    return Eh;
}

/****************************************
 Function name	  : paw_get_hartree_pot
 Description	    :
****************************************/
double* paw_get_hartree_pot()
{

    return Vh;

}

