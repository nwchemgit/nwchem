/*
   $Id$
*/

/************************************
 REVISION LOG ENTRY
 Revision By: Marat Valiev
 Revised on 4/8/99 10:10:56 PM
 Comments: ...
************************************/
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
#include   "paw_atom.h"
#include   "paw_sdir.h"

static int nbasis;
static int  max_orb_l;
static int *prin_n;
static int *prin_n_ps;
static int *orb_l;
static int *orb_type;
static int *l_counter;

static char projector_method[20];
static char nodal_constraint[20];

static int projectors_done = False;
static int pseudo_orbitals_done = False;

static int max_i_r_orbital;

static int *i_r_orbital;

static double* r_orbital;

static double *log_deriv;
static double *e;
static double *e_ps;
static double *fill;
static double *log_deriv;

static double **phi0;
static double **phi0_prime;
static double **phi_ps0;
static double **phi_ps0_prime;

static double **phi;
static double **phi_ps;
static double **phi_prime;
static double **phi_ps_prime;

static double **prj_ps0;
static double **prj_ps;
static double *rho_ps;
static double *rho;

static double **psi_ps;
static double **psi_ps_prime;

static double ekin;
static int *delta_ekin;

static double *scaling_factor;
static double **tr_matrix;

static double Zvalence;

/****************************************
 Function name	  : paw_init_paw_orbitals
 Description	    :
 Return type		  : void
 Argument         : int a_nbasis
 Argument         : double *a_rc_orb
 Argument         : int *a_n
 Argument         : int *a_n_ps
 Argument         : int *a_l
 Argument         : int *a_s_z
 Argument         : double *a_e
 Argument         : double *a_fill

 Author     		  : Marat Valiev
 Date & Time		  : 4/4/99 5:11:29 PM
****************************************/
void paw_init_paw_basis(
    char* a_nodal_constraint,
    char* a_projector_method,
    int a_nbasis,
    int *a_n,
    int *a_l,
    double *r_match
)
{

    int i;
    int k;
    int i_match;
    int Ngrid;
    int index;
    double *rgrid;
    double *psi;
    double *psi_prime;


    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    strcpy(nodal_constraint,a_nodal_constraint);
    strcpy(projector_method,a_projector_method);

    nbasis = a_nbasis;
    /*  prin_n = a_n; */
    /*  orb_l  = a_l; */
    /* r_orbital = r_match; */



    /*maximum angular momentum in the basis*/
    max_orb_l = 0;
    for (i=0;i<=nbasis-1;i++)
    {
        if (max_orb_l < a_l[i])
            max_orb_l = a_l[i];
    }

    delta_ekin  = (int *) malloc( nbasis*      sizeof(int));

    orb_type    = (int *)    malloc( nbasis*      sizeof(int));
    prin_n_ps   = (int *)    malloc( nbasis*      sizeof(int));
    prin_n      = (int *)    malloc( nbasis*      sizeof(int));
    orb_l       = (int *)    malloc( nbasis*      sizeof(int));
    l_counter   = (int *)    malloc((max_orb_l+1)*sizeof(int));
    i_r_orbital = (int *)    malloc( nbasis*      sizeof(int));
    r_orbital   = (double *) malloc( nbasis*      sizeof(double));
    fill        = (double *) malloc( nbasis*      sizeof(double));
    e           = (double *) malloc( nbasis*      sizeof(double));
    e_ps        = (double *) malloc( nbasis*      sizeof(double));
    log_deriv   = (double *) malloc( nbasis*      sizeof(double));

    rho         = paw_alloc_LogGrid();
    rho_ps      = paw_alloc_LogGrid();


    phi          = (double **) malloc(nbasis * sizeof(double *));
    phi0         = (double **) malloc(nbasis * sizeof(double *));
    phi_prime    = (double **) malloc(nbasis * sizeof(double *));
    phi_ps0      = (double **) malloc(nbasis * sizeof(double *));
    phi_ps0_prime= (double **) malloc(nbasis * sizeof(double *));
    phi0_prime   = (double **) malloc(nbasis * sizeof(double *));
    phi_ps       = (double **) malloc(nbasis * sizeof(double *));
    phi_ps_prime = (double **) malloc(nbasis * sizeof(double *));
    prj_ps       = (double **) malloc(nbasis * sizeof(double *));
    prj_ps0      = (double **) malloc(nbasis * sizeof(double *));

    psi_ps       = (double **) malloc(nbasis * sizeof(double *));
    psi_ps_prime = (double **) malloc(nbasis * sizeof(double *));

    for (i = 0; i < nbasis; ++i)
    {
        prin_n[i] = a_n[i];
        orb_l[i]  = a_l[i];

        phi0[i]         = paw_alloc_LogGrid();
        phi[i]          = paw_alloc_LogGrid();
        phi_prime[i]    = paw_alloc_LogGrid();

        phi_ps0[i]      = paw_alloc_LogGrid();
        phi_ps[i]       = paw_alloc_LogGrid();
        phi_ps_prime[i] = paw_alloc_LogGrid();

        phi_ps0_prime[i] = paw_alloc_LogGrid();

        phi0_prime[i] = paw_alloc_LogGrid();

        prj_ps[i]       = paw_alloc_LogGrid();
        prj_ps0[i]      = paw_alloc_LogGrid();

        psi_ps[i]          = paw_alloc_LogGrid();
        psi_ps_prime[i] = paw_alloc_LogGrid();

    }


    Zvalence = 0.0;
    for (i=0;i<=nbasis-1;i++)
    {
        index = paw_get_orbital_index(prin_n[i],orb_l[i]);

        fill[i]      = paw_get_fill(index);
        Zvalence    += fill[i];
        e[i]         = paw_get_e(index);
        e_ps[i]      = e[i];
        orb_type[i]  = paw_get_orb_type(index);

        psi = paw_get_psi(index);
        for (k=0;k<=Ngrid-1;k++)
            phi[i][k] = psi[k];

        psi_prime = paw_get_psi_prime(index);
        for (k=0;k<=Ngrid-1;k++)
            phi_prime[i][k] = psi_prime[k];

    }


    /* set matching point for pseudoorbitals*/
    for (i=0;i<=nbasis-1;i++)
    {
        i_r_orbital[i] = paw_get_grid_index(r_match[i]);
        /*r_orbital[i] = rgrid[i_r_orbital[i]];  */
        r_orbital[i]   = r_match[i];
    }

    /*largest sphere in the basis*/
    max_i_r_orbital = 0;
    for (i=0;i<=nbasis-1;i++)
    {
        if (max_i_r_orbital < i_r_orbital[i])
            max_i_r_orbital = i_r_orbital[i];
    }


    /*check if prin_n array is monotonically increasing*/
    for (i=0;i<=nbasis-2;i++)
    {
        if (prin_n[i]>prin_n[i+1])
            printf("please order your states according to increasing n");

    }

    /* counter for number of orbitals per angular momentum*/
    for (i=0;i<(max_orb_l+1);++i)
        l_counter[i] = 0;

    if (strcmp(nodal_constraint,"off")==0)
    {

        for (i=0;i<nbasis;++i)
        {

            prin_n_ps[i] = orb_l[i] + 1;

        }


    }
    else if (strcmp(nodal_constraint,"on")==0)
    {

        for (i=0;i<nbasis;++i)
        {

            prin_n_ps[i] = orb_l[i] + 1 + l_counter[orb_l[i]];
            l_counter[orb_l[i]] = l_counter[orb_l[i]]+1;

        }

    }
    else
    {

        printf("unknown value for the nodal_constraint\n");
        exit(1);

    }



    for (i = 0; i < nbasis; ++i)
    {
        i_match = i_r_orbital[i];

        if (fabs(phi[i][i_match]) < SMALL)
        {
            printf("error, your orbital matching sphere is to close to the node of %d%d orbital \n",
                   prin_n[i],orb_l[i]);
            exit(1);
        }
        else
        {
            log_deriv[i] = phi_prime[i][i_match]/phi[i][i_match];
        }

    }

    /*set all paw orbitals to scattering ***/
    for (i = 0; i < nbasis; ++i)
    {
        orb_type[i]=scattering;
    }

    scaling_factor = paw_alloc_1d_array(nbasis);
    tr_matrix      = paw_alloc_2d_array(nbasis,nbasis);

}

void paw_end_paw_basis()
{
    int i;


    for (i = 0; i < nbasis; ++i)
    {
        paw_dealloc_LogGrid(phi0[i]);
        paw_dealloc_LogGrid(phi[i]);
        paw_dealloc_LogGrid(phi_prime[i]);

        paw_dealloc_LogGrid(phi_ps0[i]);
        paw_dealloc_LogGrid(phi_ps[i]);
        paw_dealloc_LogGrid(phi_ps_prime[i]);

        paw_dealloc_LogGrid(phi_ps0_prime[i]);

        paw_dealloc_LogGrid(phi0_prime[i]);

        paw_dealloc_LogGrid( prj_ps[i]);
        paw_dealloc_LogGrid(prj_ps0[i]);

        paw_dealloc_LogGrid(psi_ps[i]);
        paw_dealloc_LogGrid( psi_ps_prime[i]);
    }
    free(phi);
    free(phi0);
    free(phi_prime);
    free(phi_ps0);
    free(phi_ps0_prime);
    free(phi0_prime);
    free(phi_ps);
    free(phi_ps_prime);
    free(prj_ps);
    free(prj_ps0);
    free(psi_ps);
    free(psi_ps_prime);


    free(delta_ekin);
    free(orb_type);
    free(prin_n_ps);
    free(prin_n);
    free(orb_l);
    free(l_counter);
    free(i_r_orbital);
    free(r_orbital);
    free(fill);
    free(e);
    free(e_ps);
    free(log_deriv);

    paw_dealloc_LogGrid(rho);
    paw_dealloc_LogGrid(rho_ps);

}


/****************************************
 Function name	  : solve_pseudo_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 4/7/99 2:31:16 PM
****************************************/
void paw_generate_pseudo_orbitals()
{
    int     i;
    int     k;
    int Ngrid;
    int     max_iter;
    int     converged;
    int     iteration;
    int     match;
    int     status;
    double log_amesh;
    double f1,f2;
    double  *V;
    double  *rgrid;
    double  *f;
    double norm;

    char    output[300];
    FILE   *fp;

    Ngrid     = paw_N_LogGrid();
    rgrid     = paw_r_LogGrid();
    log_amesh = paw_log_amesh_LogGrid();

    /*set maximum number of iterations*/
    max_iter = 600;

    for (i = 0; i <= nbasis-1; i++)
    {
        iteration = 0;
        converged = False;
        status    = False;
        match = i_r_orbital[i];

        /*Start the selfconsistent loop over given ps state */

        while ((iteration <= max_iter) && (!converged))
        {
            ++iteration;

            /* get pseudopotential*/
            V = paw_get_paw_potential(i);

            /*
            if(orb_type[i]==bound || orb_type[i]==virt)
            {
              status = paw_R_Schrodinger(prin_n_ps[i], 
                orb_l[i], 
                V, 
                &e_ps[i], 
                phi_ps[i], 
                phi_ps_prime[i]);
              
            }
            else if(orb_type[i]==scattering)
            {
              status = paw_R_Schrodinger_Fixed_Logderiv(prin_n_ps[i], 
                orb_l[i], 
                V, 
                match,      
                log_deriv[i],
                &e_ps[i], 
                phi_ps[i], 
                phi_ps_prime[i]);
                  
            }
            else
            {
              printf("unknown orbital type\n");
              exit(1);
            }
            */

            status = paw_R_Schrodinger_Fixed_Logderiv(prin_n_ps[i],
                     orb_l[i],
                     V,
                     match,
                     log_deriv[i],
                     &e_ps[i],
                     phi_ps[i],
                     phi_ps_prime[i]);

            /*Update pseudopotential potential*/
            paw_update_paw_potential(&converged, i, e[i], e_ps[i],phi_ps[i]);

        }

        /*report on convergence status*/
        if (converged && (status))
        {
            if (paw_debug()) printf("\n%d%d pseudo orbital with the eigenvalue=%f has been found\n",
                                        prin_n[i], orb_l[i], e_ps[i]);
        }
        else
        {
            if (paw_debug()) printf("Unable to find %d%d orbital\n ",
                                        prin_n[i], orb_l[i]);
        }

        norm = phi[i][match]/phi_ps[i][match];

        /*scale ps orbital */
        for (k=0;k<=match;k++)
        {
            phi_ps[i][k]=phi_ps[i][k]*norm;
            phi_ps_prime[i][k]=phi_ps_prime[i][k]*norm;
        }


        for (k=match+1;k<=Ngrid-1;k++)
        {
            phi_ps[i][k]=phi[i][k];
            phi_ps_prime[i][k]=phi_prime[i][k];
        }


    }





    if (paw_debug())
    {
        for (i = 0; i <= nbasis-1; ++i)
        {
            sprintf(output, "%s%s%d%d", paw_sdir(),"test",prin_n[i],orb_l[i]);
            fp = fopen(output, "w+");
            for (k = 0; k < Ngrid; k++)
            {
                fprintf(fp, "%f\t%f\t%f\n", rgrid[k], phi[i][k],phi_ps[i][k]);

            }
            fclose(fp);
        }
    }

    pseudo_orbitals_done = True;

    /* save original basis functions */
    for (i=0;i<=nbasis-1;i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
        {
            phi0[i][k] = phi[i][k];
            phi0_prime[i][k] = phi_prime[i][k];
            phi_ps0[i][k] = phi_ps[i][k];
            phi_ps0_prime[i][k] = phi_ps_prime[i][k];
        }

    }

    /* calculate densities */
    paw_Zero_LogGrid(rho);
    paw_Zero_LogGrid(rho_ps);

    for (i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
        {
            rho[k] += fill[i]*pow((phi0[i][k]/rgrid[k]),2.0);
            rho_ps[k] += fill[i]*pow((phi_ps0[i][k]/rgrid[k]),2.0);
        }
    }

    /*calculate kinetic energy*/
    ekin    = paw_get_kinetic_energy(nbasis, orb_l, fill, phi, phi_prime);

    f = paw_alloc_LogGrid();

    for (i = 0; i <= nbasis-1; i++)
    {

        for (k=0; k<=Ngrid-1; k++)
        {
            f[k] = 0.5*phi_prime[i][k]/(rgrid[k]*log_amesh)*
                   phi_prime[i][k]/(rgrid[k]*log_amesh) -
                   0.5*phi_ps_prime[i][k]/(rgrid[k]*log_amesh)*
                   phi_ps_prime[i][k]/(rgrid[k]*log_amesh);
        }

        f1 = paw_Def_Integr(0.0,f,0.0,Ngrid-1);

        for (k=0; k<=Ngrid-1; k++)
        {
            f[k] = 0.5*phi_prime[i][k]/(rgrid[k]*log_amesh)*
                   phi_prime[i][k]/(rgrid[k]*log_amesh);
        }

        f2 = paw_Def_Integr(0.0,f,0.0,Ngrid-1);

        delta_ekin[i] = (int)(100*f1/f2);

        if (paw_debug())
            printf("kinetic energy of the %d%s orbital was reduced by %d%s\n",
                   prin_n[i],paw_spd_Name(orb_l[i]),delta_ekin[i],"%");
    }

    paw_dealloc_LogGrid(f);

}

/****************************************
 Function name	  :   generate_projectors
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 4/7/99 2:31:54 PM
****************************************/
void    paw_generate_projectors()
{
    int i;
    int j;
    double norm;
    double norm_error;

    if (strcmp(projector_method,"vanderbilt")==0)
    {
        paw_generate_projectors_vanderbilt();
    }
    else if (strcmp(projector_method,"blochl")==0)
    {
        paw_generate_projectors_blochl();
    }
    else
    {

        printf("error, unknown projector_method in generate_projectors \n");
        exit(1);

    }

    /*Check paw basis orthogonality*/
    norm_error = 0.0;
    for (i = 0; i <= nbasis-1; i++)
    {
        for (j = 0; j <= nbasis-1; j++)
        {

            if (orb_l[i] == orb_l[j])
            {
                norm = paw_dot_product(prj_ps[i],phi_ps[j]);

                if (i==j) norm = norm - 1.0;

                if (fabs(norm) > norm_error) norm_error=fabs(norm);


            }

        }
    }

    if (paw_debug())
        printf("paw basis is orthogonal to within %le \n",norm_error);

    projectors_done = True;
}

void    paw_generate_projectors_vanderbilt()
{
    int i;
    int j;
    int k;
    int Ngrid;
    double norm;
    double max_prj;
    double max_phi;
    double max_phi_ps;
    double tmp_prj;
    double *V_ref;
    double *V;
    double *rgrid;
    double **b;
    double **prj_ps0;
    char  output[300];
    FILE  *fp;


    /*make sure the pseudo orbitalas were found*/
    if (! pseudo_orbitals_done)
    {
        printf("cannot calculate projectors\n");
        printf("pseudo basis is not ready\n");
        exit(1);

    }

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    b = paw_alloc_2d_array(nbasis,nbasis);
    prj_ps0 = paw_alloc_2d_array(nbasis,Ngrid);

    /*obtain the ref potential*/
    V_ref = paw_get_ref_pot();

    /*form initial guess for projectors as in Blochl*/
    for (i = 0; i < nbasis; ++i)
    {
        V = paw_get_paw_potential(i);

        for (k = 0; k < Ngrid; ++k)
            prj_ps[i][k] = (V_ref[k] - V[k])*phi_ps[i][k];

        if (paw_debug())
        {
            sprintf(output, "%s%s_%d%d", paw_sdir(),"p",prin_n[i],orb_l[i]);
            fp = fopen(output, "w+");
            for (k = 0; k < Ngrid; ++k)
            {
                fprintf(fp, "%le\t  %le\t %le \n", rgrid[k], prj_ps[i][k],(V_ref[k] - V[k]));
            }
            fclose(fp);
        }


    }

    /*check for accidental null projectors*/
    for (i = 0; i < nbasis; ++i)
    {

        max_prj = 0.0;

        for (k = 0; k < Ngrid; ++k)
        {
            tmp_prj = fabs(prj_ps[i][k]);
            if (tmp_prj > max_prj)
                max_prj = tmp_prj;
        }

        if (max_prj <0.000001)
        {
            printf("found null projectors, aborting the program ... \n");
            exit(1);
        }
    }





    for (i=0;i<=nbasis-1;i++)
    {
        for (j=0;j<=nbasis-1;j++)
        {
            if (orb_l[i] == orb_l[j] )
            {
                b[i][j] = paw_dot_product(phi_ps[j],prj_ps[i]);

            }
            else
            {
                b[i][j] = 0.0;
            }



        }
    }

    paw_get_inverse(b,nbasis);

    for (i=0;i<=nbasis-1;i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            prj_ps0[i][k] = prj_ps[i][k];
    }

    for (i=0;i<=nbasis-1;i++)
    {
        paw_Zero_LogGrid(prj_ps[i]);
    }

    for (i=0;i<=nbasis-1;i++)
    {
        for (j=0;j<=nbasis-1;j++)
        {

            for (k = 0; k <= Ngrid-1; k++)
                prj_ps[i][k] = prj_ps[i][k] + b[i][j]*prj_ps0[j][k];

        }

    }

    for (i=0;i<=nbasis-1;i++)
    {

        norm = paw_dot_product(prj_ps[i],phi_ps[i]);
        if (fabs(norm) < SMALL)
        {
            printf("division by zero while normalizing prj_pss");
            exit(99);
        }

        if (paw_debug()) printf("prj_ps norm=%le\n",norm);

        for (k = 0; k < Ngrid; ++k)
            prj_ps[i][k] = prj_ps[i][k]/norm;

    }

    /*rescale basis*/
    for (i = 0; i <= nbasis-1; i++)
    {

        max_prj = 0.0;
        max_phi = 0.0;
        max_phi_ps = 0.0;

        for (k = 0; k < i_r_orbital[i]; k++)
        {
            if (max_prj < fabs(prj_ps[i][k]))
                max_prj =  fabs(prj_ps[i][k]);

            if (max_phi_ps < fabs(phi_ps[i][k]))
                max_phi_ps =  fabs(phi_ps[i][k]);

        }

        if (max_prj == 0.0 || max_phi_ps == 0.0)
        {
            printf("division by zero while rescaling basis\n");
            exit(99);
        }

        norm = 0.75*sqrt(max_phi_ps/max_prj);

        max_prj = 0.0;
        max_phi = 0.0;
        max_phi_ps = 0.0;

        for (k = 0; k < Ngrid; ++k)
        {
            prj_ps[i][k]       = prj_ps[i][k]*norm;
            phi[i][k]          = phi[i][k]/norm;
            phi_ps[i][k]       = phi_ps[i][k]/norm;
            phi_prime[i][k]    = phi_prime[i][k]/norm;
            phi_ps_prime[i][k] = phi_ps_prime[i][k]/norm;

        }

    }



}

/*************************





***************************/
void    paw_generate_projectors_blochl()
{
    int i;
    int j;
    int k;
    int Ngrid;
    double norm;
    double max_prj;
    double max_phi;
    double max_phi_ps;
    double tmp_prj;
    double *V_ref;
    double *V;
    double *rgrid;
    double **b;
    double **L;
    double **U;
    double **L_inv;
    double **U_inv;
    double **test_matrix;
    char  output[300];
    FILE  *fp;


    /*make sure the pseudo orbitalas were found*/
    if (! pseudo_orbitals_done)
    {
        printf("cannot calculate projectors\n");
        printf("pseudo basis is not ready\n");
        exit(1);

    }

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    b = paw_alloc_2d_array(nbasis,nbasis);
    L = paw_alloc_2d_array(nbasis,nbasis);
    L_inv = paw_alloc_2d_array(nbasis,nbasis);
    U = paw_alloc_2d_array(nbasis,nbasis);
    U_inv = paw_alloc_2d_array(nbasis,nbasis);
    test_matrix = paw_alloc_2d_array(nbasis,nbasis);

    /*obtain the ref potential*/
    V_ref = paw_get_ref_pot();

    /*form initial guess for projectors as in Blochl*/
    for (i = 0; i < nbasis; ++i)
    {
        V = paw_get_paw_potential(i);

        for (k = 0; k < Ngrid; ++k)
            prj_ps[i][k] = (V_ref[k] - V[k])*phi_ps[i][k];

        if (paw_debug())
        {
            sprintf(output, "%s%s_%d%d", paw_sdir(),"p",prin_n[i],orb_l[i]);
            fp = fopen(output, "w+");
            for (k = 0; k < Ngrid; ++k)
            {
                fprintf(fp, "%le\t  %le\t %le \n", rgrid[k], prj_ps[i][k],(V_ref[k] - V[k]));
            }
            fclose(fp);
        }

    }

    /*check for accidental null projectors*/
    for (i = 0; i < nbasis; ++i)
    {

        max_prj = 0.0;

        for (k = 0; k < Ngrid; ++k)
        {
            tmp_prj = fabs(prj_ps[i][k]);
            if (tmp_prj > max_prj)
                max_prj = tmp_prj;
        }

        if (max_prj <0.000001)
        {
            printf("found null projectors, aborting the program ... \n");
            exit(1);
        }
    }





    for (i=0;i<=nbasis-1;i++)
    {
        for (j=0;j<=nbasis-1;j++)
        {
            if (orb_l[i] == orb_l[j] )
            {
                b[i][j] = paw_dot_product(phi_ps[j],prj_ps[i]);

            }
            else
            {
                b[i][j] = 0.0;
            }

        }
    }

    paw_lu_decompose(nbasis, b, L, U);
    paw_triang_matrix_inverse("l",nbasis,L,L_inv);

    paw_triang_matrix_inverse("u",nbasis,U,U_inv);

    for (i=0;i<=nbasis-1;i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            prj_ps0[i][k] = prj_ps[i][k];
    }

    for (i=0;i<=nbasis-1;i++)
    {
        paw_Zero_LogGrid(prj_ps[i]);
        paw_Zero_LogGrid(phi[i]);
        paw_Zero_LogGrid(phi_ps[i]);
        paw_Zero_LogGrid(phi_prime[i]);
        paw_Zero_LogGrid(phi_ps_prime[i]);
    }

    for (i=0;i<=nbasis-1;i++)
    {
        for (j=0;j<=i;j++)
        {

            for (k = 0; k <= Ngrid-1; k++)
            {
                prj_ps[i][k] = prj_ps[i][k] + L_inv[i][j]*prj_ps0[j][k];

                phi[i][k]          = phi[i][k]          + U_inv[j][i]*phi0[j][k];
                phi_ps[i][k]       = phi_ps[i][k]       + U_inv[j][i]*phi_ps0[j][k];
                phi_prime[i][k]    = phi_prime[i][k]    + U_inv[j][i]*phi0_prime[j][k];
                phi_ps_prime[i][k] = phi_ps_prime[i][k] + U_inv[j][i]*phi_ps0_prime[j][k];

            }


        }

    }

    for (i=0;i<=nbasis-1;i++)
        paw_Zero_LogGrid(prj_ps0[i]);

    for (i=0;i<=nbasis-1;i++)
    {
        for (j=0;j<=nbasis-1;j++)
        {

            for (k = 0; k <= Ngrid-1; k++)
            {
                prj_ps0[i][k] = prj_ps0[i][k] + U_inv[i][j]*prj_ps[j][k];

            }


        }

    }


    /*rescale basis*/
    for (i = 0; i <= nbasis-1; i++)
    {

        max_prj = 0.0;
        max_phi = 0.0;
        max_phi_ps = 0.0;

        for (k = 0; k < i_r_orbital[i]; k++)
        {
            if (max_prj < fabs(prj_ps[i][k]))
                max_prj =  fabs(prj_ps[i][k]);

            if (max_phi_ps < fabs(phi_ps[i][k]))
                max_phi_ps =  fabs(phi_ps[i][k]);

        }

        if (max_prj == 0.0 || max_phi_ps == 0.0)
        {
            printf("division by zero while rescaling basis\n");
            exit(99);
        }

        norm = 0.75*sqrt(max_phi_ps/max_prj);

        max_prj = 0.0;
        max_phi = 0.0;
        max_phi_ps = 0.0;

        scaling_factor[i] = norm;

        for (k = 0; k < Ngrid; ++k)
        {
            prj_ps[i][k]       = prj_ps[i][k]*norm;
            phi[i][k]          = phi[i][k]/norm;
            phi_ps[i][k]       = phi_ps[i][k]/norm;
            phi_prime[i][k]    = phi_prime[i][k]/norm;
            phi_ps_prime[i][k] = phi_ps_prime[i][k]/norm;

        }

    }

    for (i = 0; i <= nbasis-1; i++)
    {
        for  (j = 0; j <= nbasis-1; j++)
        {
            tr_matrix[i][j] = U_inv[i][j];
        }
    }
    if (paw_debug())
        printf("derivative=%f\n",phi_ps_prime[0][0]/(rgrid[0]*paw_log_amesh_LogGrid()));

}

/****************************************
 Function name	  : paw_solve_pseudo_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 5/15/00
****************************************/
void paw_solve_pseudo_orbitals()
{

    int i;
    int k;
    int i_end_point;
    int Ngrid;
    double scale;
    double *rgrid;


    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();


    for (i=0;i<= nbasis-1;i++)
    {

        i_end_point = paw_get_grid_index(1.2*r_orbital[i]);


        paw_solve_paw_scattering(orb_l[i], 3*r_orbital[i], e_ps[i], psi_ps[i],psi_ps_prime[i]);


        if (paw_debug())
            printf("%d%s log derivative   %le  %le \n",prin_n[i],paw_spd_Name(orb_l[i]),
                   psi_ps_prime[i][i_end_point]/(psi_ps[i][i_end_point]*rgrid[i_end_point]*paw_log_amesh_LogGrid()),
                   phi0_prime[i][i_end_point]/(phi0[i][i_end_point]*rgrid[i_end_point]*paw_log_amesh_LogGrid()));

        /*rescale the orbitals*/
        scale = phi_ps0[i][i_r_orbital[i]]/psi_ps[i][i_r_orbital[i]];
        for (k=0;k<=i_end_point;++k)
        {
            psi_ps[i][k]       = psi_ps[i][k]*scale;
            psi_ps_prime[i][k] = psi_ps_prime[i][k]*scale;

        }

        for (k=i_end_point+1;k<=Ngrid-1;++k)
        {
            psi_ps[i][k]       = phi_ps0[i][k];
            psi_ps_prime[i][k] = phi_ps0_prime[i][k];

        }

    }


}



/****************************************
 Function name	  : paw_find_density
 Description	    :
 Return type		  : void
 Argument         : double **rho_ps
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 1:13:57 PM
****************************************/
double* paw_get_pointer_paw_ps_density()
{

    return rho_ps;


}

/****************************************
 Function name	  : paw_find_density
 Description	    :
 Return type		  : void
 Argument         : double **rho_ps
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 1:13:57 PM
****************************************/
double* paw_get_pointer_paw_density()
{

    return rho;


}

double paw_get_paw_kinetic_energy()
{

    return ekin;

}
/****************************************
 Function name	  : paw_print_basis_to_file
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:01:57 PM
****************************************/
void paw_print_basis_to_file(char* atom_name)
{
    int i;
    int j;
    int k;
    int Ngrid;
    double *rgrid;
    char data_filename[300];
    char script_filename[300];
    char nl_name[20];
    char title[20];
    FILE *fp;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    if (paw_debug())
    {
        sprintf(data_filename,"%s%s_paw.dat",paw_sdir(),atom_name);
        fp = fopen(data_filename,"w+");

        for (k=0; k<=Ngrid-1; k++)
        {
            fprintf(fp,"%le", rgrid[k]);

            for (i=0; i<=nbasis-1; i++)
                fprintf(fp,"\t%le \t%le \t%le ",phi[i][k],phi_ps[i][k],prj_ps[i][k]);

            fprintf(fp,"\n");
        }
        fclose(fp);


        /* individual orbitals gnu script file */
        for (i=0,j=2; i<=nbasis-1; i++,j=j+3)
        {

            sprintf(nl_name,"%d%s",prin_n[i],paw_spd_Name(orb_l[i]));
            sprintf(script_filename,"%s%s_%s_paw.plt",paw_sdir(),atom_name,nl_name);
            fp = fopen(script_filename,"w+");

            fprintf(fp,"set style data lines \n");
            fprintf(fp,"set nolabel \n");
            fprintf(fp,"set autoscale \n");
            fprintf(fp,"set xr[0:%f] \n",1.5*r_orbital[i]);
            fprintf(fp,"set grid \n");
            fprintf(fp,"set key \n");
            fprintf(fp,"set nolabel \n");

            fprintf(fp,"set xlabel \"r (a0)\" \n");
            fprintf(fp,"set title \" %s paw_basis for %s\\n (e=%f Hartree)\" \n",
                    nl_name,atom_name,e_ps[i]);

            sprintf(title,"%s (%f)",nl_name,e_ps[i]);

            fprintf(fp,"plot \"%s\" using 1:%d title \" all-electron orbital \" ,  ",data_filename,j);
            fprintf(fp,"\"%s\" using 1:%d title \" pseudo orbital \" ,  ",data_filename,j+1);
            fprintf(fp,"\"%s\" using 1:%d title \" projector \" \n",data_filename,j+2);

            fprintf(fp,"\n");
            fprintf(fp,"pause -1\n");
            fclose(fp);

        }
    }
}

/****************************************
 Function name	  : paw_print_basis_to_file
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:01:57 PM
****************************************/
void paw_print_basis_test_to_file(char* atom_name)
{
    int i;
    int j;
    int k;
    int Ngrid;
    double *rgrid;
    char data_filename[300];
    char script_filename[300];
    char nl_name[20];
    char title[20];
    FILE *fp;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    if (paw_debug())
    {
        sprintf(data_filename,"%s%s_test.dat",paw_sdir(),atom_name);
        printf("data_filename: %s\n",data_filename);
        fp = fopen(data_filename,"w+");

        for (k=0; k<=Ngrid-1; k++)
        {
            fprintf(fp,"%le", rgrid[k]);

            for (i=0; i<=nbasis-1; i++)
                fprintf(fp,"\t%le \t%le",phi_ps0[i][k],psi_ps[i][k]);



            fprintf(fp,"\n");

        }

        fclose(fp);


        /*gnu script file */
        for (i=0,j=2; i<=nbasis-1; i++,j=j+2)
        {

            sprintf(nl_name,"%d%s",prin_n[i],paw_spd_Name(orb_l[i]));

            sprintf(script_filename,"%s%s_%s_test.plt",paw_sdir(),atom_name,nl_name);

            printf("script_filename: %s \n",script_filename);
            fp = fopen(script_filename,"w+");

            fprintf(fp,"set style data lines \n");
            fprintf(fp,"set nolabel \n");
            fprintf(fp,"set autoscale \n");
            fprintf(fp,"set xr[0:%f] \n",1.5*r_orbital[i]);
            fprintf(fp,"set grid \n");
            fprintf(fp,"set key \n");
            fprintf(fp,"set nolabel \n");

            fprintf(fp,"set xlabel \"r (a0)\" \n");
            fprintf(fp,"set title \" %s paw_basis test for %s\\n (e=%f Hartree)\" \n",
                    nl_name,atom_name,e_ps[i]);

            sprintf(title,"%s (%f)",nl_name,e_ps[i]);

            fprintf(fp,"plot \"%s\" using 1:%d title \" original pseudo orbital \" ,  ",data_filename,j);
            fprintf(fp,"\"%s\" using 1:%d title \" solution of PAW Hamiltonian \"  with points\n",data_filename,j+1);

            fprintf(fp,"\n");
            fprintf(fp,"pause -1\n");
            fclose(fp);

        }
    }
}

double paw_get_Zvalence() { return Zvalence; }


int paw_get_nbasis()
{

    return nbasis;

}

int* paw_get_pointer_paw_l_array()
{

    return orb_l;
}

int* paw_get_pointer_paw_n_array()
{

    return prin_n;
}

int* paw_get_pointer_paw_n_ps_array()
{

    return prin_n_ps;
}

double* paw_get_pointer_paw_e_array()
{

    return e_ps;
}

double** paw_get_pointer_paw_psi_array()
{

    return phi;

}

double** paw_get_pointer_paw_psi_prime_array()
{

    return phi_prime;

}

double** paw_get_pointer_paw_psi_ps_array()
{

    return phi_ps;

}

double** paw_get_pointer_paw_psi_ps_prime_array()
{

    return phi_ps_prime;

}

double** paw_get_pointer_paw_prj_ps_array()
{

    return prj_ps;

}

double** paw_get_pointer_paw_prj_ps0_array()
{

    return prj_ps0;

}


int paw_get_max_i_r_orbital()
{

    return max_i_r_orbital;

}
double paw_get_r_orbital(int i)
{
    return r_orbital[i];
}


int paw_projectors_are_done()
{
    return projectors_done;
}

double** paw_get_pointer_paw_psi_ps_unscr_array()
{

    return phi_ps0;

}

void paw_print_basis_information(FILE *fp)
{
    int i;

    fprintf(fp,"\n");

    fprintf(fp," Basis information :\n");
    fprintf(fp,"\n");

    fprintf(fp,"   nodal constraint = %s\n",
            nodal_constraint);

    fprintf(fp,"   projector method = %s\n",
            projector_method);

    fprintf(fp,"\n");

    fprintf(fp,"   nl    e               fill      dEkin\n");
    for (i=0; i<nbasis; ++i)
    {
        fprintf(fp,"   %d%s   %le  %f   %d%s\n",prin_n[i],paw_spd_Name(orb_l[i]),
                e_ps[i],fill[i],delta_ekin[i],"%");
    }

    fprintf(fp,"\n\n");
}

/****************************************
 Function name	  : paw_solve_pseudo_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 5/15/00
****************************************/
void paw_scattering_test(double e1,double e2,int number_points ,int l, double r )
{

    int i;
    int k;
    int i_end_point;
    int Ngrid;
    double *rgrid;
    FILE   *fp;
    double *V_ks;


    double *psi1;
    double *psi1_prime;

    double *psi;
    double *psi_prime;

    double *log_grid_ae;
    double *log_grid_paw;

    double de;
    double* e3;
    double e_test;
    double log_amesh;

    /*char output[30];*/
    char data_filename[300];
    char script_filename[300];

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();
    log_amesh = paw_log_amesh_LogGrid();

    psi1       = paw_alloc_LogGrid();
    psi1_prime = paw_alloc_LogGrid();

    psi       = paw_alloc_LogGrid();
    psi_prime = paw_alloc_LogGrid();

    log_grid_ae = paw_alloc_1d_array(number_points);
    log_grid_paw = paw_alloc_1d_array(number_points);
    e3 = paw_alloc_1d_array(number_points);

    de = (e2-e1)/number_points;
    V_ks = paw_get_kohn_sham_potential();

    i_end_point = paw_get_grid_index(r);

    for (i=0;i<= number_points-1;i++)
    {
        e_test = e1+de*i;
        e3[i] = e_test;


        paw_solve_paw_scattering(l, r, e_test, psi,psi_prime);

        paw_R_Schrodinger_Fixed_E(
            l,
            V_ks,
            i_end_point,
            e_test,
            psi1,
            psi1_prime
        );




        log_grid_ae[i] =  psi1_prime[i_end_point-1]/(psi1[i_end_point-1]*rgrid[i_end_point-1]*log_amesh);
        log_grid_paw[i]  =  psi_prime[i_end_point-1]/(psi[i_end_point-1]*rgrid[i_end_point-1]*log_amesh);

    }

    if (paw_debug())
    {
        sprintf(data_filename,"%s%s_%s_scat_test.dat", paw_sdir(),paw_get_atom_name(),paw_spd_Name(l));
        fp = fopen(data_filename,"w+");

        for (k=0; k<=number_points-1; k++)
        {
            fprintf(fp,"%le\t%le\t%le\n", e3[k],log_grid_ae[k],log_grid_paw[k]);

        }
        fclose(fp);

        sprintf(script_filename,"%s%s_%s_scat_test.plt", paw_sdir(),paw_get_atom_name(),paw_spd_Name(l));
        printf("script_filename: %s\n",script_filename);
        fp = fopen(script_filename,"w+");

        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[%f:%f] \n",e1,e2);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"e (Hartree)\" \n");
        fprintf(fp,"set ylabel \"logarithmic derivative at r=%f\" \n",r);
        fprintf(fp,"set title \" %s %s channel scattering test\n",paw_get_atom_name(),paw_spd_Name(l));

        fprintf(fp,"plot \"%s\" using 1:2 title \"all electron\",",data_filename);
        fprintf(fp,"\"\" using 1:3 title \"paw\" \n");

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);
    }


}
