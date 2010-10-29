/*
   $Id$
*/

#include  <stdio.h>
#include  <math.h>
#include  "paw_loggrid.h"
#include  "paw_pred_cor.h"
#include  "paw_hartree.h"


static double Zion;
static double Eion;
static double *Vion;

/****************************************
 Function name	  : paw_init_ion(double Z)
 Description	    :
****************************************/
void paw_init_ion(double Z)
{
    int i;
    int Ngrid;
    double *rgrid;


    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    Zion = Z;

    Vion  = paw_alloc_LogGrid();
    for (i=0; i<Ngrid; ++i)
        Vion[i] = -Z/rgrid[i];

}

/****************************************
 Function name	  : paw_get_ion_energy
 Description	    :
 Return type		  : double
 Argument         : double *dn
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 2:12:25 PM
****************************************/
double paw_get_ion_energy(double *dn)
{
    int    k;
    int Ngrid;
    double *tmp;


    Ngrid = paw_N_LogGrid();
    tmp   = paw_scratch_LogGrid();

    for (k=0; k<Ngrid; ++k)
    {
        tmp[k] = Vion[k]*dn[k];
    }

    Eion = paw_Integrate_LogGrid(tmp);
    return Eion;
}

double* paw_get_ion_pot()
{
    return Vion;
}

double paw_get_ion_charge()
{

    return Zion;


}

