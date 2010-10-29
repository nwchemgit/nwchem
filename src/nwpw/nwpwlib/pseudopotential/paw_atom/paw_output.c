/*
   $Id$
*/

#include   <stdio.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_loggrid.h"
#include   "paw_atom.h"
#include   "paw_ion.h"
#include   "paw_basis.h"
#include   "paw_core.h"
#include   "paw_comp_charge.h"
#include   "paw_my_constants.h"
#include   "paw_potential.h"
#include   "paw_sdir.h"


void paw_generate_basis_file(char *outfile)
{
    char*   atom_name;
    FILE    *fp;
    int     i;
    int     k;
    double tmp;
    double *Vpseudo;
    double *rho_core;
    double *rho_core_ps;

    int nbasis;
    int Ngrid;
    double *rgrid;
    int* prin_n;
    int* prin_n_ps;
    int *l;
    double* e;
    double** psi;
    double** psi_ps0;
    double** psi_prime;
    double** psi_ps;
    double** psi_ps_prime;
    double** prj_ps;
    double** prj_ps0;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    atom_name = paw_get_atom_name();
    nbasis    = paw_get_nbasis();


    prin_n    = paw_get_pointer_paw_n_array();
    prin_n_ps = paw_get_pointer_paw_n_ps_array();
    l         = paw_get_pointer_paw_l_array();
    e         = paw_get_pointer_paw_e_array();

    psi          = paw_get_pointer_paw_psi_array();
    psi_ps       = paw_get_pointer_paw_psi_ps_array();
    psi_prime    = paw_get_pointer_paw_psi_prime_array();
    psi_ps_prime = paw_get_pointer_paw_psi_ps_prime_array();
    prj_ps       = paw_get_pointer_paw_prj_ps_array();
    prj_ps0       = paw_get_pointer_paw_prj_ps0_array();

    rho_core    = paw_get_pointer_core_density();
    rho_core_ps = paw_get_pointer_ps_core_density();

    Vpseudo = paw_get_pointer_pseudopotential();

    psi_ps0 =  paw_get_pointer_paw_psi_ps_unscr_array();

    /* output the basis file */
    /*sprintf(output, "%s_basis", atom_name);*/
    if (paw_debug()) printf("paw basis file generated: %s\n",outfile);
    fp = fopen(outfile, "w+");

    fprintf(fp,"4\n"); /*dummy tag*/    /* new*/
    fprintf(fp,"%s\n",atom_name);    /* new*/
    fprintf(fp,"%lf\n",paw_get_Zvalence());    /* new*/
    fprintf(fp,"%15.11e\n",rgrid[0]);

    fprintf(fp,"%15.11e\n",rgrid[Ngrid-1]);

    fprintf(fp,"%d\n",Ngrid);

    fprintf(fp,"%d\n",nbasis);

    /* printout cutoff radii */    /* new*/
    for (i=0; i<nbasis; ++i) fprintf(fp,"%le ",paw_get_r_orbital(i));    /* new*/
    fprintf(fp,"\n");    /* new*/

    fprintf(fp,"%d\n",paw_get_max_i_r_orbital());
    fprintf(fp,"%s\n",paw_get_comment());    /* new*/

    fprintf(fp,"%15.11e\n",paw_get_core_kinetic_energy());

    for ( i = 0; i < nbasis; ++i)
        fprintf(fp, "%d\t %15.11e\t %d\t %d\n",prin_n[i],e[i],prin_n_ps[i],l[i]);

    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            fprintf(fp, "%15.11e \n",psi[i][k]);
    }


    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
        {
            tmp = psi_prime[i][k]/(rgrid[k]*paw_log_amesh_LogGrid());
            fprintf(fp, "%15.11e \n",psi_prime[i][k]/(rgrid[k]*paw_log_amesh_LogGrid()));
        }
    }


    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            fprintf(fp, "%15.11e \n",psi_ps[i][k]);
    }


    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
        {
            tmp = psi_ps_prime[i][k]/(rgrid[k]*paw_log_amesh_LogGrid());
            fprintf(fp, "%15.11e \n",psi_ps_prime[i][k]/(rgrid[k]*paw_log_amesh_LogGrid()));
        }
    }


    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            fprintf(fp, "%15.11e \n",prj_ps[i][k]);
    }

    for (k = 0; k <= Ngrid-1; k++)
    {
        fprintf(fp, "%15.11e \n",rho_core[k]/(4.0*PI));
    }


    for (k = 0; k <= Ngrid-1; k++)
    {
        fprintf(fp, "%15.11e \n",rho_core_ps[k]/(4.0*PI));
    }


    for (k = 0; k <= Ngrid-1; k++)
    {
        fprintf(fp, "%15.11e \n",Vpseudo[k]);
    }


    fprintf(fp,"%15.11e\n",paw_get_sigma_comp());

    fprintf(fp,"%15.11e\n",paw_get_ion_charge());


    for ( i = 0; i <= nbasis-1; i++)
    {
        for (k = 0; k <= Ngrid-1; k++)
            fprintf(fp, "%15.11e \n",prj_ps0[i][k]);
    }

    fclose(fp);

}




