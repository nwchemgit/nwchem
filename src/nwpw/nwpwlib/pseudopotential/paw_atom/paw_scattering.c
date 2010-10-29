/*
   $Id$
*/

#include   <stdio.h>
#include   <string.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_loggrid.h"
#include   "paw_schrodin.h"
#include   "paw_my_constants.h"
#include   "paw_potential.h"
#include   "paw_core.h"
#include   "paw_utilities.h"
#include   "paw_orbitals.h"
#include   "paw_kinetic_energy.h"
#include   "paw_my_memory.h"
#include   "paw_matrix_elements.h"
#include   "paw_get_inverse.h"
#include   "paw_utilities.h"
#include   "paw_basis.h"
#include   "paw_scattering.h"

static int paw_scattering_initialized = False;

static double *u;
static double *u_prime;
static double **w;
static double **w_prime;
static double **A;
static double **B;
static double **D;
static double *g;
static double *c;

void paw_init_paw_scattering_set()
{
    paw_scattering_initialized = False;
}

void paw_init_paw_scattering()
{

    int nbasis;
    int Ngrid;

    Ngrid = paw_N_LogGrid();

    nbasis = paw_get_nbasis();

    u       = paw_alloc_LogGrid();
    u_prime = paw_alloc_LogGrid();

    w       = paw_alloc_2d_array(nbasis,Ngrid);
    w_prime = paw_alloc_2d_array(nbasis,Ngrid);

    g = paw_alloc_1d_array(nbasis);
    c = paw_alloc_1d_array(nbasis);

    A = paw_alloc_2d_array(nbasis,nbasis);
    B = paw_alloc_2d_array(nbasis,nbasis);
    D = paw_alloc_2d_array(nbasis,nbasis);

    paw_scattering_initialized = True;

}

void paw_end_paw_scattering()
{

    int nbasis;
    int Ngrid;

    if (paw_scattering_initialized)
    {
        Ngrid  = paw_N_LogGrid();
        nbasis = paw_get_nbasis();

        paw_dealloc_LogGrid(u);
        paw_dealloc_LogGrid(u_prime);

        paw_dealloc_2d_array(nbasis,Ngrid,w);
        paw_dealloc_2d_array(nbasis,Ngrid,w_prime);

        paw_dealloc_1d_array(g);
        paw_dealloc_1d_array(c);

        paw_dealloc_2d_array(nbasis,nbasis,A);
        paw_dealloc_2d_array(nbasis,nbasis,B);
        paw_dealloc_2d_array(nbasis,nbasis,D);
        paw_scattering_initialized = False;
    }

}


/****************************************

****************************************/
void paw_solve_paw_scattering(int l, double r, double e, double* psi,double* psi_prime)
{

    int i1;
    int j1;
    int s;
    int k;
    int i_end_point;
    int nbasis;
    int Ngrid;
    int *orb_l;
    double *V;
    double **prj_ps;

    double** dH_ae;
    double** dH_ps;
    double** dO_ae;
    double** dO_ps;

    if (!(paw_scattering_initialized))
        paw_init_paw_scattering();

    Ngrid = paw_N_LogGrid();

    nbasis = paw_get_nbasis();
    orb_l  = paw_get_pointer_paw_l_array();
    prj_ps = paw_get_pointer_paw_prj_ps_array();

    dH_ae = paw_get_pointer_dH_ae();
    dH_ps = paw_get_pointer_dH_ps();
    dO_ae = paw_get_pointer_dO_ae();
    dO_ps = paw_get_pointer_dO_ps();


    paw_Zero_LogGrid(u);
    paw_Zero_LogGrid(u_prime);

    i_end_point = paw_get_grid_index(r);

    V = paw_get_ref_pot();

    paw_Zero_LogGrid(u);
    paw_Zero_LogGrid(u_prime);


    paw_R_Schrodinger_Fixed_E(
        l,
        V,
        i_end_point,
        e,
        u,
        u_prime
    );



    for (i1 = 0; i1 <= nbasis-1; i1++)
    {
        paw_Zero_LogGrid(w[i1]);
        paw_Zero_LogGrid(w_prime[i1]);

        if (l==orb_l[i1])
        {
            paw_R_Schrodinger_Fixed_E1(
                l,
                V,
                prj_ps[i1],
                i_end_point,
                e,
                w[i1],
                w_prime[i1]
            );

        }

    }



    for (i1 = 0; i1 <= nbasis-1; i1++)
    {
        for (j1 = 0; j1 <= nbasis-1; j1++)
        {

            A[i1][j1] = 0.0;
            B[i1][j1] = 0.0;
            D[i1][j1] = 0.0;

        }

    }

    for (i1 = 0; i1 <= nbasis-1; i1++)
    {
        if (orb_l[i1]==l)
            g[i1] = paw_dot_product(prj_ps[i1],u);
        else
            g[i1] = 0.0;

    }

    for (i1 = 0; i1 <= nbasis-1; i1++)
    {
        for (j1 = 0; j1 <= nbasis-1; j1++)
        {

            A[i1][j1] = dH_ae[i1][j1]-e*dO_ae[i1][j1]
                        -dH_ps[i1][j1] + e*dO_ps[i1][j1];

        }

    }



    for (i1 = 0; i1 < nbasis; ++i1)
    {
        for (j1 = 0; j1 < nbasis; ++j1)
        {
            if (orb_l[i1]==orb_l[j1] && orb_l[i1]==l)
            {
                B[i1][j1] = paw_dot_product(prj_ps[i1],w[j1]);
            }
            else
            {

                B[i1][j1] = 0.0;

            }
        }
    }

    for (i1 = 0; i1 < nbasis; ++i1)
    {
        D[i1][i1] = 1.0;
        for (j1 = 0; j1 < nbasis; ++j1)
        {
            for (s=0;s<nbasis;++s)
                D[i1][j1] = D[i1][j1] + A[i1][s]*B[s][j1];

        }

    }

    paw_get_inverse(D, nbasis);

    for (i1 = 0; i1 < nbasis; ++i1)
    {
        c[i1] = 0.0;
        for (j1 = 0; j1 < nbasis; ++j1)
        {
            for (s=0;s<nbasis;++s)
            {

                c[i1] = c[i1] - D[i1][j1]*A[j1][s]*g[s];

            }

        }

    }


    for (k=0;k<=i_end_point;++k)
    {

        psi[k] = u[k];
        psi_prime[k] = u_prime[k];


    }

    for (k=0;k<=i_end_point;++k)
    {
        for (i1 = 0; i1 < nbasis; ++i1)
        {

            psi[k]       = psi[k] + c[i1]*w[i1][k];
            psi_prime[k] = psi_prime[k] + c[i1]*w_prime[i1][k];


        }

    }


}
