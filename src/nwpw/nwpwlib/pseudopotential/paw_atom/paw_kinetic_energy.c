/*
   $Id$
*/

#include   <stdio.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_loggrid.h"

double paw_get_kinetic_energy(int num_states, int *l, double *fill, double **psi, double **psi_prime)
{
    int    i;
    int    k;
    int Ngrid;
    double ekin;
    double ekin_total;
    double log_amesh;
    double *f;
    double *r;

    Ngrid     = paw_N_LogGrid();
    log_amesh = paw_log_amesh_LogGrid();
    r         = paw_r_LogGrid();

    f = paw_alloc_LogGrid();

    ekin_total = 0.0;

    for (i=0; i<=num_states-1; i++)
    {
        for (k=0; k<=Ngrid-1; k++)
        {
            f[k] = 0.5*psi_prime[i][k]/(r[k]*log_amesh)*
                   psi_prime[i][k]/(r[k]*log_amesh)
                   +0.5*l[i]*(l[i]+1)/(r[k]*r[k])*psi[i][k]*psi[i][k];
        }

        ekin  = paw_Def_Integr(0.0,f,0.0,Ngrid-1);
        ekin_total = ekin_total + ekin*fill[i];

    }

    paw_dealloc_LogGrid(f);

    return ekin_total;
}


