/*
   $Id$
*/


#include   <stdio.h>
#include   <string.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_utilities.h"
#include   "paw_loggrid.h"

#include   "paw_atom.h"
#include   "paw_potential.h"
#include   "paw_orbitals.h"
#include   "paw_dirac_exchange.h"
#include   "paw_vosko.h"
#include   "paw_hartree.h"
#include   "paw_ion.h"
#include   "paw_kinetic_energy.h"
#include   "paw_sdir.h"

static char     comment[80];
static char     *atom_name;
static double   Zion;


static double   Total_E;
static double   E_kin;
static double   E_ion;
static double   E_Hartree;
static double   E_exchange;
static double   E_correlation;


/****************************************
 Function name	  :  paw_init_atom
 Description	    :
 Return type		  : void
 Argument         : char *aname
 Argument         : FILE   *fp
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 10:51:31 PM
****************************************/
void   paw_init_atom(char* atom, char *infile)
{

    int p,p1;
    char input[50];
    char *w,*tc;
    FILE *fp;


    if ( (fp  = fopen(infile, "r" )) == NULL )
    {
        printf("atom output file %s does not exist\n", infile);
        exit(1);
    }
    else
    {
        if (paw_debug())
            printf("atom input file %s was found\n", infile);
    }

    /* set the name of the atom */
    atom_name = atom;

    /* read ion charge */
    /*fscanf(fp,"%le", &Zion);*/
    strcpy(input,"<atom_charge>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("please specify atom charge\n");
        printf("Aborting the program\n");
        exit(1);

    }
    else
    {
        fscanf(fp,"%le", &Zion);
    }

    /*initialize log grid*/
    paw_init_LogGrid_from_file(Zion,fp);

    paw_init_ion(Zion);
    paw_init_hartree();
    paw_init_dirac_exchange();
    paw_init_vosko();

    /*initialize potential*/
    paw_init_potential(Zion);

    /*initialize orbitals*/
    paw_init_orbitals_from_file(fp);

    fclose(fp);

   strcpy(comment,"PAW Hamann pseudopotential");
   fp = fopen(infile,"r+");
   w = paw_get_word(fp);
   while ((w!=NIL) && (strcmp("<comment>",w)!=0))
      w = paw_get_word(fp);

   if(w!=NIL)
   {
      w = paw_get_word(fp);
      p  = 0;
      tc = comment;
      while ((w!=NIL)&&(strcmp("<end>",w) != 0))
      {
        p = (strlen(w));
        strcpy(tc, w);
        for(p1=0;p1<p; ++p1) ++tc;
        strcpy(tc, " ");
        ++tc;

        w = paw_get_word(fp);
      }
   }
   fclose(fp);


}




/****************************************
 Function name	  :  paw_solve_atom
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 1:36:54 PM
 Comments         : this is a modification
                    of Eric's code
****************************************/
void   paw_solve_atom()
{
    int nbasis;
    int *l;
    double *fill;
    double **psi;
    double **psi_prime;
    double *rho;


    /*find bound orbitals*/
    paw_solve_occupied_orbitals();

    /*find unoccupied orbitals*/
    paw_solve_unoccupied_orbitals();

    /*find scattering orbitals*/
    paw_solve_scattering_orbitals();

    paw_print_orbitals_to_file(atom_name);

    /* calculate total energy */
    nbasis = paw_get_Ntotal();
    l = paw_get_pointer_l_array();
    fill = paw_get_pointer_fill_array();
    psi  = paw_get_pointer_psi_array();
    psi_prime = paw_get_pointer_psi_prime_array();
    rho       = paw_get_pointer_density();

    E_kin         = paw_get_kinetic_energy(nbasis, l, fill, psi, psi_prime);
    E_ion         = paw_get_ion_energy(rho);
    E_Hartree     = paw_get_hartree_energy(rho);
    E_exchange    = paw_get_exchange_energy_LDA(rho);
    E_correlation = paw_get_correlation_energy_LDA(rho);

    Total_E = E_kin + E_ion + E_Hartree
              + E_exchange
              + E_correlation;

    if (paw_debug())
        printf("total energy = %f\n", Total_E);

    paw_print_atom();

}



/****************************************
 Function name	  :  paw_print_Atom
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 2:12:00 PM
****************************************/
void   paw_print_atom()
{

    char output[300];
    FILE *fp;


    if (paw_debug())
    {
        sprintf(output,"%s%s_out",paw_sdir(),atom_name);
        fp = fopen(output,"w+");

        fprintf(fp,"####################################################\n");
        fprintf(fp,"##                                                ##\n");
        fprintf(fp,"##        SINGLE ATOM GROUND STATE CALCULATION    ##\n");
        fprintf(fp,"##                                                ##\n");
        fprintf(fp,"####################################################\n");

        fprintf(fp,"atom name                   : %s\n",atom_name);
        fprintf(fp,"ion charge                  : %f\n",Zion);

        paw_print_loggrid_information(fp);

        paw_print_orbital_information(fp);

        fprintf(fp,"************ENERGY ************\n");

        fprintf(fp,"E_kinetic     = %12.6f\n",E_kin);
        fprintf(fp,"\n");

        fprintf(fp,"E_coulomb     = %12.6f\n",E_ion + E_Hartree);
        fprintf(fp,"\n");

        fprintf(fp,"E_exc         = %12.6f\n",E_exchange + E_correlation);
        fprintf(fp,"\n");
        fprintf(fp,"\nTotal E     = %12.6f\n",Total_E);
        fprintf(fp,"\n");



        fclose(fp);

    }


} /* print_Atom */

double paw_get_atom_kinetic_energy()
{

    return E_kin;

}

double paw_get_atom_total_energy()
{

    return Total_E;

}

char* paw_get_atom_name()
{

    return atom_name;

}
char* paw_get_comment()
{

    return comment;

}


/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

#include   <stdio.h>
#include   <string.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_utilities.h"
#include   "paw_loggrid.h"
#include   "paw_orbitals.h"
#include   "paw_basis.h"
#include   "paw_core.h"
#include   "paw_comp_charge.h"
#include   "paw_potential.h"
#include   "paw_matrix_elements.h"
#include   "paw_my_constants.h"
#include   "paw_ion.h"
#include   "paw_atom.h"


/****************************************
 Function name    : paw_init_paw_atom
 Description        :
 Return type              : void
 Argument         : char *aname
 Argument         : FILE   *fp
 Author                   : Marat Valiev
 Date & Time              : 3/30/99 10:51:31 PM
****************************************/
void   paw_init_paw_atom(char *infile)
{
    int i;
    int nbasis;
    char nodal_constraint[20];
    char projector_method[30];
    int *n;
    int *l;
    double c0;
    double r_pot;
    double r_comp;
    double *r_orb;
    char input[50];
    FILE *fp;
    double *rgrid;
    double *V_ks;


    rgrid = paw_r_LogGrid();

    /* open the atom output file */

    if ( (fp  = fopen(infile, "r" )) == NULL )
    {
        printf("atom output file %s does not exist\n", input);
        exit(1);
    }
    else
    {
        if (paw_debug())
            printf("reading paw basis parameters \n");
    }

    /*atom_name = aname;*/

    /* paw basis */
    strcpy(input,"<paw_basis>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("error, <paw_basis> section was not found \n");
        exit(1);
    }
    else
    {
        fscanf(fp, "%d", &nbasis);
        n     = (int *) malloc( nbasis*sizeof(int));
        l     = (int *) malloc( nbasis*sizeof(int));
        r_orb = (double *) malloc( nbasis*sizeof(double));

        for (i=0; i<=nbasis-1; i++)
            fscanf(fp,"%d %d %lf",&n[i],&l[i],&r_orb[i]);

    }

    /* potential matching radius */
    strcpy(input,"<ref_potential_matching_radius>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("Aborting the program\n");
        exit(1);
    }
    else
    {
        fscanf(fp, "%le", &r_pot);
    }

    /* potential at zero */
    strcpy(input,"<ref_potential_at_zero>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("Aborting the program\n");
        exit(1);
    }
    else
    {
        fscanf(fp, "%le", &c0);
    }

    /* compensation charge radius */
    strcpy(input,"<compensation_charge_radius>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("Aborting the program\n");
        exit(1);
    }
    else
    {
        fscanf(fp, "%le", &r_comp);
    }

    /* node constraint*/
    strcpy(input,"<nodal_constraint>");
    if (paw_find_word(input,fp) != 0)
    {
        if (paw_debug()) printf(" nodal_constraint will be set to \"on\" \n");
        strcpy(nodal_constraint,"off");
    }
    else
    {
        fscanf(fp, "%s",nodal_constraint);

    }

    /*projector method*/
    strcpy(input,"<projector_method>");
    if (paw_find_word(input,fp) != 0)
    {
        if (paw_debug()) printf(" <projector_method> will be set to Vanderbilt \n");
        strcpy(projector_method,"vanderbilt");
    }
    else
    {
        fscanf(fp, "%s",projector_method);

    }

    fclose(fp);


    paw_init_core();

    paw_init_comp_charge(r_comp);

    paw_init_paw_basis(nodal_constraint,projector_method,nbasis,n,l,r_orb);

    V_ks = paw_get_kohn_sham_potential();

    paw_init_paw_potential(nbasis,c0,r_pot,r_orb,V_ks);

    /* deallocate memory */
    free(r_orb);
    free(n);
    free(l);

}


/**************************************



***************************************/
void paw_solve_paw_atom(char *infile)
{
    char input[30];
    int i;
    int n;
    int *l;
    int *number_points;
    double *e1;
    double *e2;
    double *r;
    FILE *fp;

    paw_generate_pseudo_orbitals();
    paw_generate_projectors();
    paw_generate_pseudopot();


    paw_generate_matrix_elements();
    paw_solve_pseudo_orbitals();

    if (paw_debug())
    {
        paw_print_paw_atom();
        paw_print_basis_to_file(atom_name);
        paw_print_paw_potential_to_file(atom_name);
        paw_print_basis_test_to_file(atom_name);
    }

    if ( (fp  = fopen(infile, "r" )) == NULL )
    {
        printf("atom output file %s does not exist\n", input);
        exit(1);
    }

    strcpy(input,"<scattering_test>");
    if (paw_find_word(input,fp) != 0)
    {
        if (paw_debug()) printf(" scattering test will not be performed \n");
    }
    else
    {
        fscanf(fp, "%d",&n);
        l =  (int *) malloc( n*sizeof(int));
        number_points =  (int *) malloc( n*sizeof(int));
        e1 = (double *) malloc( n*sizeof(double));
        e2 = (double *) malloc( n*sizeof(double));
        r  = (double *) malloc( n*sizeof(double));

        for (i=0;i<=n-1;i++)
        {

            fscanf(fp,"%d %lf %lf %lf %d",&l[i],&e1[i],&e2[i],&r[i],&number_points[i]);

        }
        for (i=0;i<=n-1;i++)
        {
            paw_scattering_test(e1[i],e2[i],number_points[i],l[i],r[i]);
        }


    }
    fclose(fp);


}


/****************************************
 Function name    :  paw_print_paw_atom
 Description        :
 Return type              : void
 Author                   : Marat Valiev
****************************************/
void   paw_print_paw_atom()
{

    char output[300];
    FILE *fp;

    if (paw_debug())
    {
        sprintf(output,"%s%s_paw",paw_sdir(),atom_name);
        fp = fopen(output,"w+");

        fprintf(fp,"####################################################\n");
        fprintf(fp,"##                                                ##\n");
        fprintf(fp,"##        PAW basis generator                     ##\n");
        fprintf(fp,"##                                                ##\n");
        fprintf(fp,"##        Author: Marat Valiev                    ##\n");
        fprintf(fp,"####################################################\n");
        fprintf(fp,"\n");
        fprintf(fp,"\n");

        fprintf(fp," Atom information :\n");
        fprintf(fp,"\n");

        fprintf(fp,"   name                : %s\n",atom_name);
        fprintf(fp,"   ion charge          : %f\n",paw_get_ion_charge());
        fprintf(fp,"   total energy (LDA)  : %f\n",paw_get_atom_total_energy());
        fprintf(fp,"\n");
        fprintf(fp,"\n");

        paw_print_loggrid_information(fp);

        paw_print_comp_charge_information(fp);

        paw_print_paw_potential_information(fp);

        paw_print_basis_information(fp);

        paw_print_core_information(fp);

        fclose(fp);

    }

}
