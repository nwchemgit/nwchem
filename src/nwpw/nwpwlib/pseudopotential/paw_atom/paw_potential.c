/*
   $Id$
*/

/************************************
  REVISION LOG ENTRY
  Revision By: ...
  Revised on 3/30/99 3:52:12 PM
  Comments: created file
************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <math.h>

#include "paw_loggrid.h"
#include "paw_potential.h"
#include "paw_vosko.h"
#include "paw_dirac_exchange.h"
#include "paw_hartree.h"
#include "paw_ion.h"
#include "paw_utilities.h"
#include "paw_sdir.h"

double    *Vi;

extern char     *atom_name;

/* global loggrid data */
extern int     Ngrid;
extern double  *rgrid;

extern   int   Hartree_Type ;
extern   int   Exchange_Type;
extern   int   Correlation_Type;


/****************************************
 Function name	  : paw_init_potential
 Description	    :
 Return type		  : void
 Argument         : double Z
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 5:05:57 PM
****************************************/
void paw_init_potential()
{

    Vi   = paw_alloc_LogGrid();


}



/****************************************
 Function name	  : paw_get_kohn_sham_potential
 Description	    :
 Return type		  : void
 Argument         : double *dn
 Argument         : double **rho
 Argument         : double **Vo
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 5:32:13 PM
****************************************/
void paw_find_kohn_sham_potential(double *rho,double *V_ks)
{
    int k;
    double    *Vion;
    double    *Vh;
    double    *Vx;
    double    *Vc;

    Vion = paw_get_ion_pot();

    paw_generate_hartree_pot(rho);
    Vh = paw_get_hartree_pot();

    paw_generate_exchange_pot_LDA(rho);
    Vx = paw_get_exchange_potential();

    paw_generate_corr_pot_LDA(rho);
    Vc = paw_get_corr_pot_LDA();

    for (k=0; k <= Ngrid-1; k++)
        V_ks[k] = Vion[k] + Vh[k] + Vx[k] + Vc[k];
}



/****************************************
 Function name	  : paw_set_Kohn_Sham_potential
 Description	    :
 Return type		  : void
 Argument         : double *rho
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 6:59:53 PM
****************************************/
void paw_set_kohn_sham_potential(double *rho)
{
    paw_find_kohn_sham_potential(rho,Vi);
}



/****************************************
 Function name	  : *paw_get_Kohn_Sham_potential
 Description	    :
 Return type		  : double
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 7:00:50 PM
****************************************/
double *paw_get_kohn_sham_potential()
{

    return Vi;

}


/****************************************
 Function name	  :   paw_Thomas_Fermi
 Description	    :
 Return type		  : void
 Argument         : double Z
 Argument         : double *V
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 3:41:34 PM
****************************************/
void    paw_Thomas_Fermi(double Z, double *V)
{
    int    i;
    double  x,t;

    for (i=0; i<Ngrid; ++i)
    {
        x = rgrid[i]* pow( (Z/0.69395656), (1.0/3.0));
        t = Z/(
                1.0+sqrt(x)*(0.02747 - x*(0.1486-0.007298*x))
                +           x*(1.243 + x*(0.2302+0.006944*x))
            );

        if (t < 1.0)
            t=1.0;

        V[i] = -t/rgrid[i];

    }

}


/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "paw_basis.h"
#include "paw_orbitals.h"
#include "paw_loggrid.h"
#include "paw_utilities.h"
#include "paw_my_constants.h"
#include "paw_bisect.h"
#include "paw_potential.h"
#include "paw_potential.h"
#include "paw_hartree.h"
#include "paw_comp_charge.h"
#include "paw_dirac_exchange.h"
#include "paw_vosko.h"
#include "paw_ion.h"
#include "paw_core.h"


static int nbasis;
static int* i_r_potential;
static double* rc_pot;
static double rc_ref;
static double lambda=6;
static double r_function_for_rc_pot;

static double *c;
static double pot_tolerance = 1.0e-10;
static double* r_potential;
static double r_ref;
static double **fcut;
static double **V_paw;
static double *V_ref;
static double *V_pseudo;

/****************************************
 Function name    : paw_function_for_rc_pot
 Description        :
 Return type              : double
 Argument         : double r
 Author                   : Marat Valiev
 Date & Time              : 4/9/99 12:50:27 PM
****************************************/
double paw_function_for_rc_pot( double r)
{

    double tmp;

    tmp = exp(-pow((r_function_for_rc_pot/r), lambda)) - pot_tolerance;

    return tmp;

}

double paw_find_rc_pot(double rcut_in)
{

    double tmp;
    double rc1;
    double rc2;

    r_function_for_rc_pot = rcut_in;

    rc1 = rcut_in/100;
    rc2 = 2*rcut_in;
    tmp = paw_bisection(paw_function_for_rc_pot, rc1, rc2, 0.00001*pot_tolerance);

    return tmp;

}


/****************************************
 Function name    : paw_init_paw_potential
 Description        :
 Return type              : void
 Argument         : int a_nbasis
 Argument         : double a_r_sphere
 Author                   : Marat Valiev
 Date & Time              : 4/11/99 4:11:34 PM
****************************************/
void paw_init_paw_potential(int a_nbasis,
                            double c0,
                            double a_r_ref,
                            double* a_r_potential,
                            double* V_ks)
{

    int   i;
    int   k;
    int    Ngrid;
    int ic;
    double rc;
    double a,d,b;
    double V_prime;
    double *rgrid;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    nbasis      = a_nbasis;

    r_ref  = a_r_ref;
    rc_ref = paw_find_rc_pot(r_ref);

    i_r_potential = (int *) malloc(nbasis * sizeof(int));
    r_potential   = (double *) malloc(nbasis * sizeof(double));
    rc_pot        = (double *) malloc(nbasis * sizeof(double));
    c             = (double *) malloc(nbasis * sizeof(double));
    fcut          = (double **) malloc(nbasis * sizeof(double *));
    V_pseudo      = paw_alloc_LogGrid();
    V_ref         = paw_alloc_LogGrid();
    V_paw         = (double **) malloc(nbasis * sizeof(double *));

    for (i = 0; i <= nbasis-1; ++i)
    {
        V_paw[i]  = paw_alloc_LogGrid();
        fcut[i]  = paw_alloc_LogGrid();
    }

    for (i=0; i <= nbasis-1;i++)
    {
        i_r_potential[i] = paw_get_grid_index(a_r_potential[i]);
        /*r_potential[i]   = rgrid[i_r_potential[i]];*/
        r_potential[i] = a_r_potential[i];
        rc_pot[i] = paw_find_rc_pot(r_potential[i]);
    }


    /* set ref potential
    for (k = 0; k <= Ngrid-1; ++k)
      V_ref[k] = V_ks[k]*(1 - exp(-pow((rgrid[k]/rc_ref),lambda)))+
                 c0*exp(-pow((rgrid[k]/rc_ref),lambda));

    */

    ic = paw_get_grid_index(a_r_ref);
    rc = rgrid[ic];

    V_prime = 0.5*(V_ks[ic+1] - V_ks[ic-1])/(rc*paw_log_amesh_LogGrid());


    a = c0;

    d = (0.5*V_prime*rc + a - V_ks[ic])/(rc*rc*rc*rc);

    b = (0.5*V_prime*rc - 2.0*d*rc*rc*rc*rc)/(rc*rc);

    for (k=0;k<ic;++k)
        V_ref[k] = a + b*rgrid[k]*rgrid[k] + d*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k];

    for (k=ic;k<Ngrid;++k)
        V_ref[k] = V_ks[k];

    /* Form cutoff function  for PS potential*/
    for (i = 0; i <= nbasis-1; i++)
        for (k = 0; k <= Ngrid-1; k++)
            fcut[i][k] = exp(-pow((rgrid[k]/rc_pot[i]),lambda));

    /*Initialize coefficients for the construction of PS potential*/
    for (i=0; i<nbasis; ++i)
        c[i] = V_ks[i_r_potential[i]];
    /*
      for (i = 0; i < nbasis; ++i)
        c[i] = V_ks[paw_get_grid_index(r_potential[i])];
    */


    /* set initial paw potential*/
    for (i = 0; i < nbasis; ++i)
        for (k = 0; k <= Ngrid-1; k++)
            V_paw[i][k] = V_ref[k] + c[i] * fcut[i][k];

}


/****************************************
 Function name    : paw_get_paw_potential
 Description        :
 Return type              : double*
 Argument         : int i
 Author                   : Marat Valiev
 Date & Time              : 1/11/99 11:11:16 AM
****************************************/
double* paw_get_paw_potential(int i)
{

    return V_paw[i];

}


double* paw_get_pointer_pseudopotential()
{

    return V_pseudo;

}

/****************************************
 Function name    : paw_update_paw_potential
 Description        :
 Return type              : void
 Argument         : int conv_status
 Argument         : int i
 Argument         : double *w
 Author                   : Marat Valiev
 Date & Time              : 1/10/99 6:27:44 PM
****************************************/
void paw_update_paw_potential(int *conv_status, int i, double eig, double eig_ps,double *w)
{

    int     k;
    int Ngrid;
    double  sv;
    double  dcl;
    double  eps;
    double  *tmp;

    eps = 1.0e-7;

    Ngrid = paw_N_LogGrid();
    tmp   = paw_scratch_LogGrid();

    for (k = 0; k <= Ngrid-1; k++)
    {
        tmp[k]= fcut[i][k] * w[k] * w[k];
    }

    sv = paw_Def_Integr(0.0,tmp,0.0,Ngrid-1);

    dcl = (eig - eig_ps) / sv;

    /*In case of the overflow rescale dcl*/

    if (fabs(dcl)>100)
        dcl = 10*fabs(dcl)/dcl;


    c[i] = c[i] + dcl;
    *conv_status = (fabs(dcl) <= eps);


    /*Form new paw potential if necessary*/
    if (!(*conv_status))
        for (k = 0; k <= Ngrid-1; k++)
            V_paw[i][k] = V_ref[k] + c[i] * fcut[i][k];

}

/****************************************
 Function name    : paw_get_ref_pot
 Description        :
 Return type              : double*
 Author                   : Marat Valiev
 Date & Time              : 1/25/99 11:29:24 AM
 Modifications    : 1/26/99
****************************************/
double* paw_get_ref_pot()
{
    return V_ref;
}


/****************************************
 Function name    : paw_generate_pseudopot
 Description        :
 Return type              : void
 Author                   : Marat Valiev
 Date & Time              : 4/10/99 7:40:00 PM
****************************************/
void paw_generate_pseudopot()
{

    int   k;
    int   Ngrid;
    double charge;
    double ps_charge;
    double Z;
    double *Vh;
    double *Vx;
    double *Vc;
    double *rho;
    double *rho_ps;
    double *rho_core;
    double *rho_core_ps;
    double *full_density;
    double *full_ps_density;
    double* V_comp;
    double *rgrid;
    FILE *fp;
    char data_filename[300];

    if ( !(paw_projectors_are_done()) )
    {
        printf("error, pseudopotential cannot be generated ");
        printf(" because projectors have not been yet \n");
        exit(1);
    }

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    Vh      = paw_alloc_LogGrid();

    Vx      = paw_alloc_LogGrid();
    Vc      = paw_alloc_LogGrid();

    full_density = paw_alloc_LogGrid();
    full_ps_density = paw_alloc_LogGrid();

    Z = paw_get_ion_charge();

    paw_set_core();

    /*get densities*/
    rho         = paw_get_pointer_paw_density();
    rho_ps      = paw_get_pointer_paw_ps_density();
    rho_core    = paw_get_pointer_core_density();
    rho_core_ps = paw_get_pointer_ps_core_density();

    paw_Zero_LogGrid(full_density);
    paw_Zero_LogGrid(full_ps_density);

    for (k=0;k<=Ngrid-1;k++)
    {
        full_density[k]    = rho[k] + rho_core[k];
        full_ps_density[k] = rho_ps[k] + rho_core_ps[k];

    }

    charge   = paw_Integrate_LogGrid(full_density);
    ps_charge      = paw_Integrate_LogGrid(full_ps_density);


    V_comp = paw_find_comp_charge_potential(Z,charge,ps_charge);

    paw_generate_hartree_pot(full_ps_density);
    Vh = paw_get_hartree_pot();

    paw_generate_exchange_pot_LDA(full_ps_density);
    Vx = paw_get_exchange_potential();

    paw_generate_corr_pot_LDA(full_ps_density);
    Vc = paw_get_corr_pot_LDA();


    /*form pseudopotential*/
    for (k=0;k<=Ngrid-1;k++)
    {

        V_pseudo[k] = V_ref[k] - Vh[k]- V_comp[k]- Vc[k] - Vx[k];



    }

    if (paw_debug())
    {
        sprintf(data_filename,"%sdensity",paw_sdir());
        fp = fopen(data_filename,"w");

        for (k=0;k<Ngrid;++k)
            fprintf(fp,"%f  %f   %f\n",rgrid[k],rho[k] - rho_ps[k],V_pseudo[k]);
        fclose(fp);
    }

    paw_dealloc_LogGrid(full_density);
    paw_dealloc_LogGrid(full_ps_density);

}

double paw_get_potential_matching_radius()
{
    return r_ref;
}

void  paw_print_paw_potential_information(FILE *fp)
{
    int i;
    int *n;
    int *l;

    fprintf(fp,"\n");

    fprintf(fp," Paw potential information :\n");
    fprintf(fp,"\n");

    fprintf(fp,"   reference potential matching radius    = %le\n",
            paw_get_potential_matching_radius());

    fprintf(fp,"   reference potential rcut parameter     = %le\n",
            rc_ref);

    fprintf(fp,"   lambda parameter                       = %le\n",
            lambda);

    fprintf(fp,"   potential tolerance                    = %le\n",
            pot_tolerance);

    fprintf(fp,"\n");

    n = paw_get_pointer_paw_n_array();
    //l = paw_get_pointer_l_array();
    l = paw_get_pointer_paw_l_array();

    fprintf(fp,"   nl     c[i]        rcut[i] \n");


    for (i=0;i<=nbasis-1;i++)
    {

        fprintf(fp,"   %d%s    %f    %f \n",n[i],paw_spd_Name(l[i]),c[i],
                rc_pot[i]);

    }

    fprintf(fp,"\n\n");

}

void paw_print_paw_potential_to_file(char* atom_name)
{
    int i;
    int j;
    int k;
    int Ngrid;
    int *prin_n;
    int * orb_l;
    double *rgrid;
    char data_filename[300];
    char script_filename[300];
    char nl_name[20];
    FILE *fp;

    if (paw_debug())
    {
        Ngrid = paw_N_LogGrid();
        rgrid = paw_r_LogGrid();

        prin_n = paw_get_pointer_paw_n_array();
        orb_l  = paw_get_pointer_paw_l_array();

        sprintf(data_filename,"%s%s_pot.dat",paw_sdir(),atom_name);
        fp = fopen(data_filename,"w+");

        for (k=0; k<=Ngrid-1; k++)
        {
            fprintf(fp,"%le\t%le\t%le", rgrid[k], V_ref[k], V_pseudo[k]);

            for (i=0; i<=nbasis-1; i++)
                fprintf(fp,"\t%le ",V_paw[i][k]);



            fprintf(fp,"\n");

        }

        fclose(fp);

        sprintf(script_filename,"%s%s_pot.plt",paw_sdir(),atom_name);
        fp = fopen(script_filename,"w+");
        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[0:%f] \n",2*r_ref);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set key \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");
        fprintf(fp,"set title \" %s potentials\n",atom_name);

        fprintf(fp,"plot \"%s\" using 1:2 title \"vref\"\n",data_filename);
        fprintf(fp,"replot \"%s\" using 1:3 title \"vlocal\"\n",data_filename);
        fprintf(fp,"replot \"%s\" using 1:($2+$3) title \"vref+vlocal\"\n",data_filename);
        for (i=0; i<=nbasis-1; i++)
           fprintf(fp,"replot \"%s\" using 1:%d title \"v%d%s\"\n",data_filename,i+4,prin_n[i],paw_spd_Name(orb_l[i]));


        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);

        fclose(fp);

        sprintf(script_filename,"%s%s_ref_pot.plt",paw_sdir(),atom_name);

        fp = fopen(script_filename,"w+");

        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[0:%f] \n",2*r_ref);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set nokey \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");
        fprintf(fp,"set title \" %s reference potential\n",atom_name);

        fprintf(fp,"plot \"%s\" using 1:2 \n",data_filename);

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);


        sprintf(script_filename,"%s%s_loc_pot.plt",paw_sdir(),atom_name);

        fp = fopen(script_filename,"w+");

        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[0:%f] \n",2*r_ref);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set nokey \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");
        fprintf(fp,"set title \" %s local potential\n",atom_name);

        fprintf(fp,"plot \"%s\" using 1:3 \n",data_filename);

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);


        /*gnu script file */
        for (i=0,j=4; i<=nbasis-1; i++,j=j+1)
        {

            sprintf(nl_name,"%d%s",prin_n[i],paw_spd_Name(orb_l[i]));

            sprintf(script_filename,"%s%s_%s_pot.plt",paw_sdir(),atom_name,nl_name);

            fp = fopen(script_filename,"w+");

            fprintf(fp,"set style data lines \n");
            fprintf(fp,"set nolabel \n");
            fprintf(fp,"set autoscale \n");
            fprintf(fp,"set xr[0:%f] \n",1.5*r_potential[i]);
            fprintf(fp,"set grid \n");
            fprintf(fp,"set nokey \n");
            fprintf(fp,"set nolabel \n");

            fprintf(fp,"set xlabel \"r (a0)\" \n");
            fprintf(fp,"set title \" %s paw potential for %s orbital\" \n",
                    atom_name,nl_name);

            fprintf(fp,"plot \"%s\" using 1:%d \n",data_filename,j);

            fprintf(fp,"\n");
            fprintf(fp,"pause -1\n");
            fclose(fp);

        }
    }
}









