/*
   $Id$
*/

/************************************
  REVISION LOG ENTRY
  Revision By: Marat Valiev
  Revised on 3/30/99 11:16:00 PM
  Comments: ...
 ************************************/

#include   <stdio.h>
#include   <string.h>
#include   <math.h>
#include   <stdlib.h>

#include   "paw_utilities.h"
#include   "paw_orbitals.h"
#include   "paw_loggrid.h"
#include   "paw_schrodin.h"
#include   "paw_pauli.h"
#include   "paw_potential.h"
#include   "paw_my_constants.h"
#include   "paw_dirac_exchange.h"
#include   "paw_vosko.h"
#include   "paw_hartree.h"
#include   "paw_kinetic_energy.h"
#include   "paw_ion.h"
#include   "paw_sdir.h"

static int   Nvalence;
static int   Nvirt;
static int   Nbound;
static int   Nscat;
static int   Ntotal;

static   int occupied_orbitals_done = False;

static   int      *n;
static   int      *l;
static   double   *fill;
static   int      *orb_type;


static  double    *eigenvalue;
static  double    **psi;
static  double    **psi_prime;
static  double    *psi_tmp;
static  double    *rho;

static int Solver_Type = Pauli;



/****************************************
 Function name	  : paw_init_orbitals
 Description	    :
 Return type		  : void
 Argument         : FILE *fp
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 11:47:04 PM
****************************************/
void paw_init_orbitals_from_file(FILE *fp)
{
    int     i;
    char   input[20];
    char *w;

    rewind(fp);

    strcpy(input,"<orbitals>");
    if (paw_find_word(input,fp) != 0)
    {
        printf("orbital section not found\n");
        printf("Aborting the program\n");
        exit(1);
    }
    else
    {
        /* read number of core, valence and virtual orbitals */
        fscanf(fp,"%d %d %d",&Nvalence,&Nvirt,&Nscat);

        Nbound   = Nvalence+Nvirt;
        Ntotal   = Nvalence+Nvirt+Nscat;

        /* set orbital indexes arrays */
        n          = (int *) malloc(Ntotal*sizeof(int));
        l          = (int *) malloc(Ntotal*sizeof(int));
        orb_type   = (int *) malloc(Ntotal*sizeof(int));
        fill       = (double *) malloc(Ntotal*sizeof(double));
        eigenvalue = (double *) malloc(Ntotal*sizeof(double));


        /*read orbital indexes arrays*/
        for (i=0; i<Nbound; ++i)
            fscanf(fp,"%d %d %lf",&n[i],&l[i],&fill[i]);

        for (i=Nbound; i<Ntotal; ++i)
            fscanf(fp,"%d %d %lf %lf",&n[i],&l[i],&fill[i],&eigenvalue[i]);

    }

    /*allocating additional memory*/
    psi       = (double **) malloc(Ntotal*sizeof(double*));
    psi_prime = (double **) malloc(Ntotal*sizeof(double*));
    for (i=0; i<Ntotal; i++)
    {
        psi[i]       = paw_alloc_LogGrid();
        psi_prime[i] = paw_alloc_LogGrid();
    }

    rho     = paw_alloc_LogGrid();

    psi_tmp   = paw_alloc_LogGrid();

    for (i=0; i<= Nvalence-1; i++)
        orb_type[i] = bound;

    for (i=Nvalence; i<= Nbound-1; i++)
        orb_type[i] = virt;

    for (i=Nbound; i<= Nbound+Nscat-1; i++)
        orb_type[i] = scattering;

    rewind(fp);
    if (paw_find_word("<solver>",fp) == 0)
    {
      w = paw_get_word(fp);
      if (strcmp ("schrodinger", w) == 0)
        Solver_Type = Schrodinger;
      if (strcmp ("pauli", w) == 0)
        Solver_Type = Pauli;
      if (strcmp ("dirac", w) == 0)
        Solver_Type = Dirac;
      if (strcmp ("zora", w) == 0)
        Solver_Type = ZORA;
    }
    else
        Solver_Type = Pauli;
}


/****************************************
 Function name	  : paw_guess_eigenvalues
 Description	    :
 Return type		  : void
 Argument         : double **V
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 11:46:56 PM
****************************************/
void paw_guess_eigenvalues(double Zion, double *V)
{
    int i;
    int Ngrid;
    double echarge;
    double Z;

    Ngrid = paw_N_LogGrid();

    echarge = 0.0;
    for (i=0; i<=Nbound-1; i++)
    {
        echarge += fill[i];
        Z = Zion - echarge + 1.0;
        eigenvalue[i] = -0.5*(Z*Z)/((double) (n[i]*n[i]));
        if (eigenvalue[i] > V[Ngrid-1])
            eigenvalue[i] = 2.0*V[Ngrid-1];
    }
}




/****************************************
 Function name	  :  paw_solve_occupied_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 6:13:22 PM
****************************************/
void   paw_solve_occupied_orbitals()
{
    int   i;
    int   j;
    int  k;
    int  it;
    int  converged;
    int   Ngrid;
    int   max_iter;
    double sum;
    double Etmp;
    double thl;
    double sn;
    double sd;
    double dr;
    double rl0;
    double rl1;
    double vn;
    double Zion;
    double *Vo;
    double *Vo1;
    double *Vi1;
    double *Vi;
    double *rgrid;

    max_iter = 100;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    /* allocate temporary grids */
    Vi      = paw_alloc_LogGrid();
    Vo1     = paw_alloc_LogGrid();
    Vi1     = paw_alloc_LogGrid();

    /* set initial guess for KS potential as Thomas-Fermi*/
    Zion = paw_get_ion_charge();
    paw_Thomas_Fermi(Zion,Vi);

    /*initial guess for the eigenvalues*/
    paw_guess_eigenvalues(Zion,Vi);

    it = 0;
    converged = False;
    while ((!converged) && (it < max_iter))
    {

        it        = it + 1;
        converged = True;

        /* solve for each of the eigenstates */
        for (i=0; i<= Nvalence-1; i++)
        {

            Etmp = eigenvalue[i];

            if (Solver_Type==Schrodinger)
            {
               paw_R_Schrodinger(n[i],l[i],Vi,
                                 &Etmp,psi[i],psi_prime[i]);
            }
            else
            {
               paw_R_Pauli(n[i],l[i],Zion,Vi,
                           &Etmp,psi[i],psi_prime[i]);
            }

            if (fabs(eigenvalue[i] - Etmp) >1.0e-10)
                converged = False;

            eigenvalue[i]=Etmp;


            /*orthogonalize to lower orbitals*/
            for (j=0;j<=i-1;j++)
            {
                if (l[i]==l[j])
                {
                    sum = paw_dot_product(psi[i],psi[j]);
                    for (k=0;k<Ngrid;++k)
                        psi[i][k] = psi[i][k] - sum*psi[j][k];
                }

            }

            /*normalize the orbital*/
            for (k=0; k<=Ngrid-1; k++)
                psi_tmp[k] = pow((psi[i][k]/rgrid[k]),2.0);

            sum = paw_Integrate_LogGrid(psi_tmp);

            for (k=0; k<=Ngrid-1; k++)
                psi[i][k]=psi[i][k]/pow(sum,0.5);

        }

        /*get new density*/
        paw_generate_density(rho);

        /*get new potential*/
        paw_set_kohn_sham_potential(rho);
        Vo = paw_get_kohn_sham_potential();


        /*****************************************/
        /* Generate the next iteration potential */
        /* using D.G. Anderson method            */
        /*****************************************/
        thl = 0.0;
        if (it > 1)
        {
            sn = 0.0;
            sd = 0.0;

            for (k=0; k<Ngrid; ++k)
            {
                rl0 = Vo[k]  - Vi[k];
                rl1 = Vo1[k] - Vi1[k];
                dr  = rl0 - rl1;
                sn = sn + rl0*dr*(rgrid[k]*rgrid[k]);
                sd = sd +  dr*dr*(rgrid[k]*rgrid[k]);
            }


            thl = sn/sd;

        }

        for (k=0; k<=Ngrid-1; k++)
        {

            vn = (1.0-0.5)*( (1.0-thl)*Vi[k] + thl*Vi1[k] )
                 + 0.5 *( (1.0-thl)*Vo[k] + thl*Vo1[k] );
            Vi1[k] = Vi[k];
            Vo1[k] = Vo[k];
            Vi[k]  = vn;

        }


    } /* end while */

    if (!converged)
    {
        printf("Not Converged\n");
        exit(1);
    }

    occupied_orbitals_done = True;
    paw_generate_density(rho);
    paw_set_kohn_sham_potential(rho);



    /* free up temporary memory */
    paw_dealloc_LogGrid(Vi);
    paw_dealloc_LogGrid(Vo1);
    paw_dealloc_LogGrid(Vi1);

}



/****************************************
 Function name	  : paw_solve_unoccupied_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 6:13:39 PM
****************************************/
void paw_solve_unoccupied_orbitals()
{
    int i;
    int j;
    int k;
    int state;
    int Ngrid;
    int status;
    double sum;
    double *V;
    double *rgrid;
    double Zion;


    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();
    Zion = paw_get_ion_charge();

    /*check if the occupied orbitals are done*/
    if (!(occupied_orbitals_done))
    {
        printf("cannot calculate unoccupied states\n");
        printf("calculate occupied states first\n");
    }

    /*get Kohn-Sham potential*/
    V = paw_get_kohn_sham_potential();

    for (i=Nvalence; i<= Nbound-1; i++)
    {

        state=paw_bound_state_test(l[i],V);
        if (state==0)
        {
            printf("This potential has no bound states with n=%d l=%d\n",n[i],l[i]);
            printf("please change your input file\n");
            exit(1);
        }
        else
        {
            if (Solver_Type==Schrodinger)
            {
               status = paw_R_Schrodinger(n[i],l[i],V,
                                          &eigenvalue[i],psi[i],psi_prime[i]);
            }
            else
            {
               status = paw_R_Pauli(n[i],l[i],Zion,V,
                                          &eigenvalue[i],psi[i],psi_prime[i]);
            }

        }

        if (!(status))
        {
            printf("This potential has no bound states with n=%d l=%d\n",n[i],l[i]);
            printf("please change your input file\n");
            exit(1);
        }


        /*orthogonalize to lower orbitals*/
        for (j=0;j<=i-1;j++)
        {
            if (l[i]==l[j])
            {
                sum = paw_dot_product(psi[i],psi[j]);
                for (k=0;k<=Ngrid-1;k++)
                    psi[i][k] = psi[i][k] - sum*psi[j][k];
            }

        }

        for (k=0;k<=Ngrid-1;k++)
            psi_tmp[k] = pow((psi[i][k]/rgrid[k]),2.0);

        sum = paw_Integrate_LogGrid(psi_tmp);

        for (k=0;k<=Ngrid-1;k++)
            psi[i][k]=psi[i][k]/pow(sum,0.5);
    }


}


/****************************************
 Function name	  : paw_solve_scattering_orbitals
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 4:58:27 PM
****************************************/
void paw_solve_scattering_orbitals()
{
    int i;
    int j;
    int k;
    int i_match;
    int Ngrid;
    double sum;
    double *V;
    double *rgrid;
    double r_sphere;
    double Zion;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();
    Zion = paw_get_ion_charge();

    /*normalization sphere radius*/
    r_sphere = 2.0;

    /*check if the occupied orbitals are done*/
    if (!(occupied_orbitals_done))
    {
        printf("cannot calculate scattering orbitals\n");
        printf("calculate occupied states first\n");
    }


    /*get Kohn-Sham potential*/
    V = paw_get_kohn_sham_potential();


    for (i=Nbound;i<Ntotal;i++)
    {

        /*set the end point*/
        i_match = Ngrid-1;

        if (Solver_Type==Schrodinger)
        {
           paw_R_Schrodinger_Fixed_E(l[i],V,i_match,eigenvalue[i],psi[i],psi_prime[i]);
        }
        else
        {
           paw_R_Pauli_Fixed_E(n[i],l[i],Zion,V,i_match,eigenvalue[i],psi[i],psi_prime[i]);
        }

        for (j=0;j<i-1;j++)
        {

            if (l[i]==l[j])
            {
                sum = paw_dot_product(psi[i],psi[j]);
                for (k=0;k<Ngrid;++k)
                    psi[i][k] = psi[i][k] - sum*psi[j][k];
            }

        }

        /*normalize*/
        sum = paw_dot_product1(paw_get_grid_index(r_sphere),psi[i],psi[i]);
        sum = 1.0/sqrt(sum);

        for (k=0;k<Ngrid;++k)
        {
            psi[i][k] = psi[i][k]*sum;
            psi_prime[i][k] = psi_prime[i][k]*sum;

        }


    }

    /*debug
    printf("orthogonality\n");
    for(i=0;i<=Ntotal-1;i++)
    {
        for(j=0;j<=Ntotal-1;j++)
        {
          if(l[i]==l[j])
          {
            printf("%d\t %d\t %f\n",i,j, paw_dot_product(psi[i],psi[j]));
          }
      
        }  
     }
        exit(1);
    end debug  */
}



/****************************************
 Function name	  : paw_bound_state_test
 Description	    :
 Return type		  : int
 Argument         : int n
 Argument         : int l
 Argument         : double *v
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:02:14 PM
****************************************/
int paw_bound_state_test(int l, double *v)
{
    int i;
    int Ngrid;
    double L2,r2;
    double Emin;
    double *r;


    Ngrid = paw_N_LogGrid();
    r      = (double *) paw_r_LogGrid();
    L2 = ((double) (l*(l+1)));
    Emin = 0;
    for (i=0; i<Ngrid; ++i)
    {
        r2 = r[i];
        r2 = r2*r2;
        Emin = Min(Emin, (v[i] + 0.5*L2/r2));
    }

    if (Emin>=0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}




/****************************************
 Function name	  : paw_print_orbital_information
 Description	    :
 Return type		  : void
 Argument         : FILE *fp
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:01:37 PM
****************************************/
void paw_print_orbital_information(FILE *fp)
{
    int i;

    if (Solver_Type==Schrodinger)
       fprintf(fp,"solver type : Schrodinger\n\n");
    else
       fprintf(fp,"solver type : Pauli\n\n");
    fprintf(fp,"\n************ ORBITALS ********\n");
    fprintf(fp,"number of valence    orbitals : %d\n",Nvalence);
    fprintf(fp,"number of virtual    orbitals : %d\n",Nvirt);
    fprintf(fp,"number of scattering orbitals : %d\n\n",Nscat);


    fprintf(fp,"n\tl\tpopulation\tEigenvalue\tDescription\n");
    for (i=0; i<Ntotal; ++i)
    {
        fprintf(fp,"%d\t%s\t%.2lf\t\t%f\t%s\n",n[i],paw_spd_Name(l[i]),
                fill[i],eigenvalue[i],paw_orbital_type_name(orb_type[i]));
    }

    fprintf(fp,"\n\n");
}






/****************************************
 Function name	  : paw_save_orbitals_to_file
 Description	    :
 Return type		  : void
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:01:57 PM
****************************************/
void paw_print_orbitals_to_file(char* atom_name)
{
    int i;
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
        sprintf(data_filename,"%s%s_orb.dat",paw_sdir(),atom_name);
        fp = fopen(data_filename,"w+");


        for (k=0; k<=Ngrid-1; k++)
        {
            fprintf(fp,"%le", rgrid[k]);

            for (i=0; i<=Ntotal-1; i++)
            {

                fprintf(fp,"\t%le",psi[i][k]);

            }
            fprintf(fp,"\n");
        }
        fclose(fp);

        /* all orbitals gnu script file */
        sprintf(script_filename,"%s%s_all_orb.plt",paw_sdir(),atom_name);

        fp = fopen(script_filename,"w+");

        fprintf(fp,
                "set title \" Kohn-Sham bound states for %s \" \n",
                atom_name);

        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set key right top box\n");

        fprintf(fp,"set grid \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");

        fprintf(fp,"set style data lines \n");

        fprintf(fp,"set autoscale \n");

        fprintf(fp,"set xr[0:5] \n");

        i=0;
        sprintf(nl_name,"%d%s",n[i],paw_spd_Name(l[i]));
        sprintf(title,"%s (%f)",nl_name,eigenvalue[i]);
        fprintf(fp,"plot \"%s\" using 1:2  title \"%s\"  ",data_filename,title);

        for (i=1; i<=Nbound-1; i++)
        {

            sprintf(nl_name,"%d%s",n[i],paw_spd_Name(l[i]));

            sprintf(title,"%s (%f)",nl_name,eigenvalue[i]);


            fprintf(fp,",\"\" using 1:%d   title \"%s\"  ",i+2,title);

        }

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");

        fclose(fp);

        /* individual orbitals gnu script file */
        for (i=0; i<=Ntotal-1; i++)
        {
            sprintf(nl_name,"%d%s",n[i],paw_spd_Name(l[i]));

            sprintf(script_filename,"%s%s_%s_orb.plt",paw_sdir(),atom_name,nl_name);

            fp = fopen(script_filename,"w+");

            fprintf(fp,"set style data lines \n");

            fprintf(fp,"set nolabel \n");

            fprintf(fp,"set autoscale \n");

            fprintf(fp,"set xr[0:%d] \n",n[i]+l[i]);

            fprintf(fp,"set xlabel \"r (a0)\" \n");

            fprintf(fp,"set grid \n");

            fprintf(fp,"set nokey \n");

            fprintf(fp,"set title \" %s orbital for %s\\n (e=%f Hartree)\" \n",
                    nl_name,atom_name,eigenvalue[i]);

            sprintf(title,"%s (%f)",nl_name,eigenvalue[i]);

            fprintf(fp,"plot \"%s\" using 1:%d   title \"%s\"  ",data_filename,i+2,title);


            fprintf(fp,"\n");
            fprintf(fp,"pause -1\n");
            fclose(fp);

        }
    }
}


int paw_get_l(int i)
{

    return l[i];

}

int paw_get_orb_type(int i)
{

    return orb_type[i];

}

int* paw_get_pointer_l_array()
{

    return l;

}

double paw_get_e(int i)
{

    return eigenvalue[i];

}

double paw_get_fill(int i)
{

    return fill[i];

}

double* paw_get_pointer_fill_array()
{

    return fill;

}

double* paw_get_psi(int i)
{

    return psi[i];

}

double** paw_get_pointer_psi_array()
{

    return psi;

}

double** paw_get_pointer_psi_prime_array()
{

    return psi_prime;

}

double* paw_get_psi_prime(int i)
{

    return psi_prime[i];

}



int paw_get_Ntotal()
{

    return Ntotal;

}



int paw_get_orbital_index(int prin_n, int orb_l)
{
    int i;
    int index;
    int orbital_found;

    index = -1;
    orbital_found = False;

    for (i=0; i <= Ntotal-1;i++)
    {
        if ( n[i]==prin_n && l[i]==orb_l)
        {
            index = i;
            orbital_found = True;
            break;
        }
    }

    if ( !(orbital_found))
    {
        printf("while generating paw basis \n");
        printf("the all electron orbital %d%s could not be located \n",prin_n,paw_spd_Name(orb_l));
        printf("please check your <orbitals> section \n");
        exit(1);
    }

    return index;

}

char* paw_orbital_type_name(int i)
{
    char *s;
    if (i==bound)
        s = "bound";
    else if (i==virt)
        s = "virtual";
    else if (i==scattering)
        s = "scattering";
    else
        s = "????";

    return s;
}


/****************************************
 Function name	  : paw_get_core_density
 Description	    :
 Return type		  : void
 Argument         : double** rho
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:35:46 PM
****************************************/
void paw_generate_density(double* rho)
{
    int i;
    int k;
    int Ngrid;
    double* rgrid;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    paw_Zero_LogGrid(rho);
    for (i=0; i<=Nvalence-1; ++i)
    {
        for (k=0; k<Ngrid; ++k)
            rho[k] += fill[i]*pow((psi[i][k]/rgrid[k]),2.0);
    }

}


/****************************************
 Function name	  : paw_get_density
 Description	    :
 Return type		  : void
 Argument         : double** rho
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:35:46 PM
****************************************/
double* paw_get_pointer_density()
{

    return rho;

}



