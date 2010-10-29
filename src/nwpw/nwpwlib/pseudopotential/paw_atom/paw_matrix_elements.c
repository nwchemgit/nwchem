/*
   $Id$
*/

#include   <stdio.h>
#include   <string.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_loggrid.h"
#include   "paw_potential.h"
#include   "paw_basis.h"
#include   "paw_my_memory.h"



static double **dH_ae;
static double **dO_ae;
static double **dO_ps;
static double **dH_ps;
static double **dT_ae;
static double **dT_ps;
static double **dV_ae;
static double **dV_ps;



/******************************


******************************/
void paw_generate_matrix_elements()
{

    int i;
    int j;
    int k;
    int Ngrid;
    int nbasis;
    double log_amesh;
    int    *orb_l;
    double **phi;
    double **phi_prime;
    double **phi_ps;
    double **phi_ps_prime;

    double *v;
    double *v_tilda;
    double *f1;
    double *f2;
    double *f3;
    double *rgrid;


    Ngrid     = paw_N_LogGrid();
    rgrid     = paw_r_LogGrid();
    log_amesh = paw_log_amesh_LogGrid();

    nbasis       = paw_get_nbasis();
    orb_l        = paw_get_pointer_paw_l_array();
    phi          = paw_get_pointer_paw_psi_array();
    phi_ps       = paw_get_pointer_paw_psi_ps_array();
    phi_prime    = paw_get_pointer_paw_psi_prime_array();
    phi_ps_prime = paw_get_pointer_paw_psi_ps_prime_array();

    f1 = paw_alloc_LogGrid();
    f2 = paw_alloc_LogGrid();
    f3 = paw_alloc_LogGrid();

    dH_ps = paw_alloc_2d_array(nbasis,nbasis);
    dO_ps = paw_alloc_2d_array(nbasis,nbasis);
    dH_ae = paw_alloc_2d_array(nbasis,nbasis);
    dO_ae = paw_alloc_2d_array(nbasis,nbasis);
    dT_ae = paw_alloc_2d_array(nbasis,nbasis);
    dT_ps = paw_alloc_2d_array(nbasis,nbasis);
    dV_ae = paw_alloc_2d_array(nbasis,nbasis);
    dV_ps = paw_alloc_2d_array(nbasis,nbasis);

    v_tilda = paw_get_ref_pot();
    v       = paw_get_kohn_sham_potential();

    for (i = 0; i <= nbasis-1; i++)
    {
        for (j = 0; j <= nbasis-1; j++)
        {
            if (orb_l[i] == orb_l[j])
            {
                for (k = 0; k <= Ngrid-1; k++)
                {

                    f1[k] = phi[i][k]*phi[j][k];

                    f2[k] = 0.5*(phi_prime[i][k]*phi_prime[j][k])/pow(rgrid[k]*log_amesh,2.0) +

                            0.5*orb_l[i]*(orb_l[i]+1)/pow(rgrid[k],2.0)*phi[i][k]*phi[j][k];

                    f3[k] = v[k]*phi[i][k]*phi[j][k];

                }

                dO_ae[i][j] = paw_Def_Integr(0.0,f1,0.0,Ngrid-1);
                dT_ae[i][j] = paw_Def_Integr(0.0,f2,0.0,Ngrid-1);
                dV_ae[i][j] = paw_Def_Integr(0.0,f3,0.0,Ngrid-1);

                for (k = 0; k < Ngrid; ++k)
                {

                    f1[k] = phi_ps[i][k]*phi_ps[j][k];

                    f2[k] = 0.5*(phi_ps_prime[i][k]*phi_ps_prime[j][k])/pow(rgrid[k]*log_amesh,2.0)+

                            0.5*orb_l[i]*(orb_l[i]+1)/pow(rgrid[k],2.0)*phi_ps[i][k]*phi_ps[j][k];


                    f3[k] = v_tilda[k]*phi_ps[i][k]*phi_ps[j][k];

                }

                dO_ps[i][j] = paw_Def_Integr(0.0,f1,0.0,Ngrid-1);
                dT_ps[i][j] = paw_Def_Integr(0.0,f2,0.0,Ngrid-1);
                dV_ps[i][j] = paw_Def_Integr(0.0,f3,0.0,Ngrid-1);

            }
            else
            {

                dO_ae[i][j] = 0.0;
                dT_ae[i][j] = 0.0;
                dV_ae[i][j] = 0.0;
                dO_ps[i][j] = 0.0;
                dT_ps[i][j] = 0.0;
                dV_ps[i][j] = 0.0;

            }

        }

    }

    for (i = 0; i < nbasis; ++i)
    {
        for (j = 0; j < nbasis; ++j)
        {

            dH_ae[i][j] = dT_ae[i][j]  + dV_ae[i][j];
            dH_ps[i][j] = dT_ps[i][j]  + dV_ps[i][j];

        }

    }

    paw_dealloc_LogGrid(f1);
    paw_dealloc_LogGrid(f2);
    paw_dealloc_LogGrid(f3);

}

double** paw_get_pointer_dH_ae()
{
    return dH_ae;
}

double** paw_get_pointer_dH_ps()
{
    return dH_ps;
}

double** paw_get_pointer_dO_ae()
{
    return dO_ae;
}

double** paw_get_pointer_dO_ps()
{
    return dO_ps;
}

