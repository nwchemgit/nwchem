/* atom.c -
   author - Eric Bylaska

*/

#include	<stdio.h>
#include	<string.h>
#include	<math.h>

#include	"name.h"
#include	"get_word.h"
#include	"loggrid.h"

#include	"schrodin.h"
#include	"pauli.h"

#include	"dft.h"
#include	"atom.h"

#define Max_Iterations	100
#define	False	0
#define	True	1
#define	Max(x,y)	((x>y) ? x : y)


/* atom structure variables */

static	int	Ncore,
		Nvalence,
		Ncv;
static	int	*n;
static	int	*l;
static	int	lmax;
static 	double	*fill;
static	int	*turning_point;
static  double	*peak;
static	double	Zion;
static	double	amass;
static	double	Total_E,
		E_Hartree, P_Hartree,
		E_exchange,
		P_exchange,
		E_correlation,
		P_correlation;	
static	double  *eigenvalue;
static	double	**r_psi;
static	double 	**r_psi_prime;
static  double  *rho;
static  double  *rho_core;
static  double  *rho_valence;
static	double	*Vion;
static	double	*Vall;
static	char	atom_name[10];
static int	solver_iterations;

/* solver parameters: this is the Kawai-Weare default */
static	int	Solver_Type      = Pauli;


/********************************
 *				*
 *	  init_Atom		*
 *				*
 ********************************/

void	init_Atom(char *filename)
{
   int	  i,Ngrid;
   double *rgrid;
   char	  *w;
   FILE	  *fp;

   /* open data file */
   fp = fopen(filename,"r+");

   /* move to <atom> section of data file */
   w = get_word(fp);
   while ((w!=NIL) && (strcmp("<atom>",w)!=0))
      w = get_word(fp);
   

   /* Error occured */
   if (w==NIL)
   {
      printf("Error: <atom> section not found\n");
      fclose(fp);
      exit(99);
   }

   /* set the name of the atom */
   fscanf(fp,"%s",atom_name);

   /* read in atom info. from the stream */
   fscanf(fp,"%le", &Zion);
   fscanf(fp,"%le", &amass);
   fscanf(fp,"%d %d",&Ncore,&Nvalence);
   Ncv	    = Ncore+Nvalence;


   /* allocate the necessary memory for eigenvalues */
   n    = (int *) malloc((Ncv+1)*sizeof(int));
   l    = (int *) malloc((Ncv+1)*sizeof(int));
   fill = (double *) malloc((Ncv+1)*sizeof(double));

   /* allocate memory for outer peak and turning_point positions */ 
   turning_point = (int *)    malloc((Ncv+1)*sizeof(int));
   peak          = (double *) malloc((Ncv+1)*sizeof(double));

   /* set eigenvalue arrays */
   lmax = 0;
   for (i=0; i<Ncv; ++i)
   {
      fscanf(fp,"%d %d %le",&n[i],&l[i],&fill[i]);
      lmax = Max(l[i],lmax);
   }

   /* set up logarithmic grid */
   init_LogGrid(Zion);
   Ngrid = N_LogGrid();
   rgrid = r_LogGrid();


   /* allocate the necessary memory */
   eigenvalue  = (double *)  malloc((Ncv+1)*sizeof(double));
   r_psi       = (double **) malloc((Ncv+1)*sizeof(double*));
   r_psi_prime = (double **) malloc((Ncv+1)*sizeof(double*));
   for (i=0; i<(Ncv+1); ++i)
   {
      r_psi[i]       = alloc_LogGrid();
      r_psi_prime[i] = alloc_LogGrid();
   }

   rho         = alloc_LogGrid();
   rho_core    = alloc_LogGrid();
   rho_valence = alloc_LogGrid();
   Vion        = alloc_LogGrid();
   Vall        = alloc_LogGrid();

   /* set Vion */
   for (i=0; i<Ngrid; ++i)
      Vion[i] = -Zion/rgrid[i];
  fclose(fp);

   /* set the solver type */
   fp = fopen(filename,"r+");
   w = get_word(fp);
   while ((w!=NIL) && (strcmp("<solver>",w)!=0))
      w = get_word(fp);
   if (w!=NIL)
   {
      w = get_word(fp);
      if (strcmp("schrodinger",w)==0) Solver_Type=Schrodinger;
      if (strcmp("pauli",w)==0)       Solver_Type=Pauli;
      if (strcmp("dirac",w)==0)       Solver_Type=Dirac;
   }
   fclose(fp);

   /* initialize DFT stuff */
   init_DFT(filename);

} /* init_Atom */


/********************************
 *				*
 *        Thomas_Fermi		*
 *				*
 ********************************/

void 	Thomas_Fermi(Z,Vtmp)

double	Z;
double 	Vtmp[];
{
   int    i,Ngrid;
   double *rgrid,
	  x,t;

   Ngrid = N_LogGrid();
   rgrid = r_LogGrid();
 
   for (i=0; i<Ngrid; ++i)
   {
      x = rgrid[i]* pow( (Z/0.69395656), (1.0/3.0));
      t = Z/( 
		1.0+sqrt(x)*(0.02747 - x*(0.1486-0.007298*x))
	      +           x*(1.243 + x*(0.2302+0.006944*x))
            );

      if (t < 1.0) 
	t=1.0;
      Vtmp[i] = -t/rgrid[i];
   }
} /* Thomas_Fermi */
   




/********************************
 *				*
 *	   solve_Atom		*
 *				*
 ********************************/

void	solve_Atom()
{
   int	i,k,Ngrid;
   int  it,converged;
   double Etmp,echarge,Z;
   double thl,sn,sd,dr,rl0,rl1,vn;
   double *rgrid;
   double *Vo,*Vo1,*Vi1;
   double *Vx,*Vc;
   double *Vh;




   /* get loggrid variables */
   Ngrid  = N_LogGrid();
   rgrid  = r_LogGrid();

   /* allocate temporary grids */
   Vh     = alloc_LogGrid();
   Vx     = alloc_LogGrid();
   Vc     = alloc_LogGrid();
   Vo     = alloc_LogGrid();
   Vo1    = alloc_LogGrid();
   Vi1    = alloc_LogGrid();


   /* using the TF potential and eigenvalues for the initial guess */
   /* get Thomas-Fermi potential */
   Thomas_Fermi(Zion,Vh);
   Copy_LogGrid(Vall, Vh);
   Zero_LogGrid(Vh);

   /* initial guess for eigenvalues */
   echarge = 0.0;
   for (i=0; i<Ncv; ++i)
   {
      echarge += fill[i];
      Z = Zion - echarge + 1.0;
      eigenvalue[i] = -0.5*(Z*Z)/((double) (n[i]*n[i]));
      if (eigenvalue[i] > Vall[Ngrid-1]) 
         eigenvalue[i] = 2.0*Vall[Ngrid-1];
   }

   E_Hartree     = 0.0;
   P_Hartree     = 0.0;
   E_exchange    = 0.0;
   E_correlation = 0.0;
   P_exchange    = 0.0;
   P_correlation = 0.0;


   it = 0;
   converged = False;
   while ((!converged) && (it < Max_Iterations))
   {
   
      it        = it + 1;
      solver_iterations = it;
      converged = True;

      /*************************************/
      /* solve for each of the eigenstates */
      /*************************************/
      for (i=0; i<Ncv; ++i)
      {
   
         Etmp = eigenvalue[i];

         /**********************************************************/
         /*          solve radial equation(s) for state i          */
         /**********************************************************/
         if (Solver_Type==Schrodinger)
         {
            R_Schrodinger(n[i],l[i],Vall,
			  &turning_point[i],&Etmp,r_psi[i],r_psi_prime[i]);
         }
         else if (Solver_Type==Pauli)
         {
            R_Pauli(n[i],l[i],Zion,Vall,
		    &turning_point[i],&Etmp,r_psi[i],r_psi_prime[i]);
         }
         /**********************************************************/

         if (eigenvalue[i] != Etmp)
            converged = False;
         eigenvalue[i] = Etmp;
      } /* solving eigenstates */

      /********************************/
      /* caluculate the total density */
      /********************************/
      Zero_LogGrid(rho);
      for (i=0; i<Ncv; ++i)
      {
         for (k=0; k<Ngrid; ++k)
            rho[k] += fill[i]*pow((r_psi[i][k]/rgrid[k]),2.0);
      } /*for i*/


      /***************************************/
      /* get the hartree potential an energy */
      /***************************************/
       P_Hartree = R_Hartree_DFT(rho,echarge,Vh);
       E_Hartree = 0.5*P_Hartree;
      
      /*****************************************/
      /* get the exchange potential and energy */
      /*****************************************/
      R_Exchange_DFT(rho,Vx,&E_exchange,&P_exchange);

      /********************************************/
      /* get the correlation potential and energy */
      /********************************************/
      R_Correlation_DFT(rho,Vc,&E_correlation,&P_correlation);
     
      
      /**********************************************************/
      /* add up all the scf potentials on to the -Z/r potential */
      /**********************************************************/
      for (k=0; k<Ngrid; ++k)
         Vo[k] = Vion[k] + Vh[k] + Vx[k] + Vc[k];


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
            rl0 = Vo[k]  - Vall[k];
            rl1 = Vo1[k] - Vi1[k];
            dr  = rl0 - rl1;
            sn = sn + rl0*dr*(rgrid[k]*rgrid[k]); 
            sd = sd +  dr*dr*(rgrid[k]*rgrid[k]); 
         } /* for k */
         thl = sn/sd;
      } /* if it>1 */
      for (k=0; k<Ngrid; ++k)
      {
         vn = (1.0-0.5)*( (1.0-thl)*Vall[k] + thl*Vi1[k] )
                 + 0.5 *( (1.0-thl)*Vo[k] + thl*Vo1[k] );
         Vi1[k] = Vall[k];
         Vo1[k] = Vo[k];
         Vall[k]  = vn;
      } /* for k */
      

   } /* end while */


   if (!converged)
      printf("Not Converged\n");


   /*************************/
   /* find the total energy */
   /*************************/
   Total_E = 0.0;
   for (i=0; i<Ncv; ++i)
   {
      Total_E += fill[i]*eigenvalue[i];
   }

   Total_E +=   E_Hartree     - P_Hartree
              + E_exchange    - P_exchange
              + E_correlation - P_correlation;


   /********************************/
   /* caluculate the core density  */
   /********************************/
   Zero_LogGrid(rho_core);
   for (i=0; i<Ncore; ++i)
   {
      for (k=0; k<Ngrid; ++k)
         rho_core[k] += fill[i]*pow((r_psi[i][k]/rgrid[k]),2.0);
   } /*for i*/

   /********************************/
   /* caluculate the valence density  */
   /********************************/
   Zero_LogGrid(rho_valence);
   for (i=Ncore; i<Ncv; ++i)
   {
      for (k=0; k<Ngrid; ++k)
         rho_valence[k] += fill[i]*pow((r_psi[i][k]/rgrid[k]),2.0);
   } /*for i*/
   

   /****************************/
   /* free up temporary memory */
   /****************************/
   dealloc_LogGrid(Vh);
   dealloc_LogGrid(Vx);
   dealloc_LogGrid(Vc);
   dealloc_LogGrid(Vo);
   dealloc_LogGrid(Vo1);
   dealloc_LogGrid(Vi1);

   /************************************************/
   /* find the outermost peak of the wavefunctions */
   /************************************************/
   for (i=0; i<Ncv; ++i)
   {
      k = Ngrid-2;
      while ((r_psi_prime[i][k]*r_psi_prime[i][k+1] >= 0.0) && (k>=0))
         --k;
      peak[i] = rgrid[k];
   }
      
} /* solve_Atom */


/********************************
 *				*
 * solve_Scattering_State_Atom  *
 *				*
 ********************************/

void	solve_Scattering_State_Atom(nt,lt,et,rmax)

int	nt,lt;
double	et;
double	rmax;
{
   double r0,al;

   r0   = r_LogGrid()[0];
   al = log_amesh_LogGrid();

   n[Ncv]	      = nt;
   l[Ncv]	      = lt;
   eigenvalue[Ncv]    = et;
   fill[Ncv]	      = 0.0;
   turning_point[Ncv] = rint(log(rmax/r0)/al);
   peak[Ncv]          = rmax;

   if (Solver_Type==Schrodinger)
   {
      R_Schrodinger_Fixed_E(nt,lt,Vall,
                            turning_point[Ncv],et,
			    r_psi[Ncv],r_psi_prime[Ncv]);
   }
   else if (Solver_Type==Pauli)
   {
      R_Pauli_Fixed_E(nt,lt,Zion,Vall,
                      turning_point[Ncv],et,
                      r_psi[Ncv],r_psi_prime[Ncv]);
   }

} /* solve_Scattering_State */


/********************************
 *				*
 *        print_Atom		*
 *				*
 ********************************/

char	*solver_Name_Atom()
{
   char *s;
   if      (Solver_Type==Schrodinger)
      s = "Schrodinger";
   else if (Solver_Type==Pauli)
      s = "Pauli";
   else
      s = "unknown?";

   return s;
}


void	print_Atom(fp)

FILE 	*fp;
{
   int i;
   
   fprintf(fp,"All electron atom solver\n\n");
   fprintf(fp,"Atom name: %s\n",atom_name);
   fprintf(fp,"Zcharge  : %le\n",Zion);
   fprintf(fp,"Amass    : %le\n",amass);
   fprintf(fp,"Ncore    : %d\n",Ncore);
   fprintf(fp,"Nvalence : %d\n",Nvalence);
   fprintf(fp,"         : restricted calculation\n");

   fprintf(fp,"\n------ Solver information ------\n");
   fprintf(fp,"solver type      : %s\n",solver_Name_Atom());
   fprintf(fp,"hartree type     : %s\n",hartree_Name_DFT());
   fprintf(fp,"exchange type    : %s\n",exchange_Name_DFT());
   if (strcmp(exchange_Name_DFT(),"Dirac")==0)
      fprintf(fp,"           alpha : %lf\n",Dirac_alpha());
   fprintf(fp,"correlation type : %s\n",correlation_Name_DFT());
   fprintf(fp,"Solver iterations: %d\n",solver_iterations);


   fprintf(fp,"\n------- Grid information -------\n");
   fprintf(fp,"Zcharge  : %le\n",Zion);
   fprintf(fp,"amesh    : %le\n", amesh_LogGrid());
   fprintf(fp,"R_max    : %le\n", r_LogGrid()[N_LogGrid()-1]);
   fprintf(fp,"Ngrid    : %d\n",N_LogGrid());
   fprintf(fp,"------------------------------------------------------------\n");
   fprintf(fp,"n\tl\tpopulation\tEigenvalue\tOuter Peak\n");
   for (i=0; i<Ncv; ++i)
   {
      fprintf(fp,"%d\t%s\t%.2lf\t\t%le\t%le\n",n[i],spd_Name(l[i]), 
				        fill[i],eigenvalue[i],peak[i]);
   }
   fprintf(fp,"------------------------------------------------------------\n");
   fprintf(fp,"charge      = %le\n",Integrate_LogGrid(rho));
   fprintf(fp,"core charge = %le\n",Integrate_LogGrid(rho_core));


   fprintf(fp,"\nTotal E       = %le\n",Total_E);
   fprintf(fp,"\n");
   

   fprintf(fp,"E_Hartree     = %le\n",E_Hartree);
   fprintf(fp,"<Vh>          = %le\n",P_Hartree);
   fprintf(fp,"\n");

   fprintf(fp,"E_exchange    = %le\n",E_exchange);
   fprintf(fp,"<Vx>          = %le\n",P_exchange);
   fprintf(fp,"\n");

   fprintf(fp,"E_correlation = %le\n",E_correlation);
   fprintf(fp,"<Vc>          = %le\n",P_correlation);
   
   
} /* print_Atom */

/********************************
 *				*
 * set_(solver parameters)_Atom	*
 *				*
 ********************************/

void set_Solver_Atom(solver)
int solver;
{
   Solver_Type = solver;
}


/********************************
 *				*
 *	      E_Atom 		*
 *				*
 ********************************/

double	E_Atom()
{
   return Total_E;
} /*E_Atom*/

double	eigenvalue_Atom(int i)
{
   return eigenvalue[i];
}

double	*rho_Atom()
{
   return rho;
}

double	*rho_core_Atom()
{
   return rho_core;
}

double	*rho_valence_Atom()
{
   return rho_valence;
}


double	*Vall_Atom()
{
   return Vall;
}


double	*r_psi_Atom(int i)
{
   return r_psi[i];
}

double	*r_psi_prime_Atom(int i)
{
   return r_psi_prime[i];
}

int	n_Atom(int i)
{
   return n[i];
}

int	l_Atom(int i)
{
   return l[i];
}

int	lmax_Atom()
{
   return lmax;
}

double	fill_Atom(int i)
{
   return fill[i];
}

int	Ncore_Atom()
{
   return Ncore;
}

int	Nvalence_Atom()
{
   return Nvalence;
}

double	peak_Atom(int i)
{
   return peak[i];
}

int	turning_point_Atom(int i)
{
   return turning_point[i];
}

double	Zion_Atom()
{
   return Zion;
}
double	Amass_Atom()
{
   return amass;
}

int	state_Atom(int nt, int lt)
{
   int i;
   i = 0;
   while ( ((nt != n[i]) || (lt != l[i])) && (i<=Ncv) )
     ++i;

   /* Error */
   if (i>Ncv)
      printf("Error: state_Atom\n");

   return i;
}

char	*name_Atom()
{
   return atom_name;
}
