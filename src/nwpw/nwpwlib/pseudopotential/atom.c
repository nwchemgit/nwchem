/* atom.c -
   author - Eric Bylaska and Patrick Nichols
   $Id$
*/

#include	<stdio.h>
#include        <stdlib.h>
#include	<string.h>

#include	"loggrid.h"
#include	"name.h"
#include	"get_word.h"

#include	"schrodin.h"
#include	"pauli.h"
#include        "dirac.h"
#include        "zora.h"
#include	"dft.h"
#include	"atom.h"

#define Max_Iterations	400
#define	False	0
#define	True	1
#define	Max(x,y)	((x>y) ? x : y)


/* atom structure variables */

static int Ncore, Nvalence, Ncv;
static int *n;
static int *l;
static int *s2;
static int lmax;
static double *fill;
static int *turning_point;
static double *peak;
static double Zion;
static double amass;
static double Total_E,
  E_Hartree, P_Hartree, E_exchange, P_exchange, E_correlation, P_correlation;
static double *eigenvalue;
static double **r_psi;
static double **r_psi_prime;
static double *rho;
static double *rho_core;
static double *rho_valence;
static double *Vion;
static double *Vall;
static char atom_name[10];
static int solver_iterations;

/* solver parameters: this is the Kawai-Weare default */
static int Solver_Type = Pauli;


/********************************
 *				*
 *	  init_Atom		*
 *				*
 ********************************/

void init_Atom (char *filename)
{
  int i, Ngrid, nx, lx, ncvh;
  double *rgrid, fillx;
  char *w;
  FILE *fp;

  /* set the solver type first */
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<solver>", w) != 0))
    w = get_word (fp);
  if (w != NIL)
    {
      w = get_word (fp);
      if (strcmp ("schrodinger", w) == 0)
	Solver_Type = Schrodinger;
      if (strcmp ("pauli", w) == 0)
	Solver_Type = Pauli;
      if (strcmp ("dirac", w) == 0)
	Solver_Type = Dirac;
      if (strcmp ("zora", w) == 0)
	Solver_Type = ZORA;
    }
  fclose (fp);
  /* open data file */
  fp = fopen (filename, "r+");

  /* move to <atom> section of data file */
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<atom>", w) != 0))
    w = get_word (fp);


  /* Error occured */
  if (w == NIL)
    {
      printf ("Error: <atom> section not found\n");
      fclose (fp);
      exit (99);
    }

  /* set the name of the atom */
  fscanf (fp, "%s", atom_name);

  /* read in atom info. from the stream */
  fscanf (fp, "%le", &Zion);
  fscanf (fp, "%le", &amass);
  fscanf (fp, "%d %d", &Ncore, &Nvalence);

  if (Solver_Type != Dirac && Solver_Type!=ZORA)
    {
      Ncv = Ncore + Nvalence;

      /* allocate the necessary memory for eigenvalues */
      n = (int *) malloc ((Ncv + 1) * sizeof (int));
      l = (int *) malloc ((Ncv + 1) * sizeof (int));
      s2 = (int *) malloc ((Ncv + 1) * sizeof (int));
      fill = (double *) malloc ((Ncv + 1) * sizeof (double));

      /* allocate memory for outer peak and turning_point positions */
      turning_point = (int *) malloc ((Ncv + 1) * sizeof (int));
      peak = (double *) malloc ((Ncv + 1) * sizeof (double));

      /* set eigenvalue arrays */
      lmax = 0;
      for (i = 0; i < Ncv; ++i)
	{
	  fscanf(fp, "%d %d %le", &n[i], &l[i], &fill[i]);
	  if (l[i] > lmax)
	    lmax = l[i];
	}

      /* set up logarithmic grid */
      init_LogGrid(Zion);
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();


      /* allocate the necessary memory */
      eigenvalue = (double *) malloc ((Ncv + 1) * sizeof (double));
      r_psi = (double **) malloc ((Ncv + 1) * sizeof (double *));
      r_psi_prime = (double **) malloc ((Ncv + 1) * sizeof (double *));
      for (i = 0; i < (Ncv + 1); ++i)
	{
	  r_psi[i] = alloc_LogGrid();
	  r_psi_prime[i] = alloc_LogGrid();
	}
    }
  else
    {

      Ncore = 2 * Ncore;
      Nvalence = 2 * Nvalence;
      Ncv = Ncore + Nvalence;

      /* allocate the necessary memory for eigenvalues */
      n = (int *) malloc((Ncv + 2) * sizeof (int));
      l = (int *) malloc((Ncv + 2) * sizeof (int));
      s2 = (int *) malloc((Ncv + 2) * sizeof (int));
      fill = (double *) malloc((Ncv + 2) * sizeof (double));

      /* allocate memory for outer peak and turning_point positions */
      turning_point = (int *) malloc((Ncv + 2) * sizeof (int));
      peak = (double *) malloc((Ncv + 2) * sizeof (double));

      /* set eigenvalue arrays */
      lmax = 0;
      ncvh = Ncv / 2;
      for (i = 0; i < ncvh; ++i)
	{
	  fscanf (fp, "%d %d %le", &nx, &lx, &fillx);
          lmax=(lmax>lx)?lmax:lx;
	  n[2 * i] = nx;
	  n[2 * i + 1] = nx;
	  l[2 * i] = lx;
	  l[2 * i + 1] = lx;
	  s2[2 * i] = -1;
	  s2[2 * i + 1] = 1;
/****************************************************************
 * fill in the j=l+0.5 and j=l-0.5 states with the proper
 * occupancy. A simple divide will give the wrong occupancies.
 ****************************************************************/
	  if (!lx)
	    {
	      fill[2 * i] = fillx * 0.5;
	      fill[2 * i + 1] = fillx * 0.5;
	    }
	  else
	    {
	      fill[2 * i] = (2.0 * lx) * (fillx * 0.5 / (2. * lx + 1.));
	      fill[2 * i + 1] = (2.0 * lx+ 2.0) * 
			(fillx * 0.5 / (2. * lx + 1.));
	    }
	}

      /* set up logarithmic grid */
      init_LogGrid(Zion);
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();


      /* allocate the necessary memory */
      eigenvalue = (double *) malloc((Ncv + 2) * sizeof (double));
      r_psi = (double **) malloc((Ncv + 2) * sizeof (double *));
      r_psi_prime = (double **) malloc((Ncv + 2) * sizeof (double *));
      for (i = 0; i < (Ncv + 2); ++i)
	{
	  r_psi[i] = alloc_LogGrid();
	  r_psi_prime[i] = alloc_LogGrid();
	}
    }
  rho = alloc_LogGrid();
  rho_core = alloc_LogGrid();
  rho_valence = alloc_LogGrid();
  Vion = alloc_LogGrid();
  Vall = alloc_LogGrid();

  /* set Vion */
  for (i = 0; i < Ngrid; ++i)
    Vion[i] = -Zion / rgrid[i];
  fclose(fp);


  /* initialize DFT stuff */
  init_DFT(filename);
}

/********************************
 *				*
 *	  end_Atom		*
 *				*
 ********************************/
void end_Atom()
{
   int i;

   dealloc_LogGrid(Vion);
   dealloc_LogGrid(Vall);
   dealloc_LogGrid(rho);
   dealloc_LogGrid(rho_core);
   dealloc_LogGrid(rho_valence);

  if (Solver_Type != Dirac && Solver_Type!=ZORA)
      for (i=0; i<(Ncv+1);++i)
      {
          dealloc_LogGrid(r_psi[i]);
          dealloc_LogGrid(r_psi_prime[i]);
      }
   else
      for (i=0; i<(Ncv+2);++i)
      {
          dealloc_LogGrid(r_psi[i]);
          dealloc_LogGrid(r_psi_prime[i]);
      }
   free(eigenvalue);
   free(r_psi);
   free(r_psi_prime);
   free(peak);
   free(turning_point);
   free(fill);
   free(s2);
   free(l);
   free(n);

   end_LogGrid();
}


/********************************
 *				*
 *        Thomas_Fermi		*
 *				*
 ********************************/

void
Thomas_Fermi (Z, Vtmp)
     double Z;
     double Vtmp[];
{
  int i, Ngrid;
  double *rgrid, x, t;

  Ngrid = N_LogGrid ();
  rgrid = r_LogGrid ();

  for (i = 0; i < Ngrid; ++i)
    {
      x = rgrid[i] * pow ((Z / 0.69395656), (1.0 / 3.0));
      t = Z / (1.0 + sqrt (x) * (0.02747 - x * (0.1486 - 0.007298 * x))
	       + x * (1.243 + x * (0.2302 + 0.006944 * x)));

      if (t < 1.0)
	t = 1.0;
      Vtmp[i] = -t / rgrid[i];
    }
}				/* Thomas_Fermi */





/********************************
 *				*
 *	   solve_Atom		*
 *				*
 ********************************/

void
solve_Atom ()
{
  int i, k, Ngrid;
  int it, converged, sz;
  double Etmp, echarge, Z;
  double thl, sn, sd, dr, rl0, rl1, vn;
  double *rgrid;
  double *Vo, *Vo1, *Vi1;
  double *Vx, *Vc;
  double *Vh;

  /* get loggrid variables */
  Ngrid = N_LogGrid ();
  rgrid = r_LogGrid ();

  /* allocate temporary grids */
  Vh = alloc_LogGrid ();
  Vx = alloc_LogGrid ();
  Vc = alloc_LogGrid ();
  Vo = alloc_LogGrid ();
  Vo1 = alloc_LogGrid ();
  Vi1 = alloc_LogGrid ();


  /* using the TF potential and eigenvalues for the initial guess */
  /* get Thomas-Fermi potential */
  Thomas_Fermi (Zion, Vh);
  Copy_LogGrid (Vall, Vh);
  Zero_LogGrid (Vh);

  /* initial guess for eigenvalues */
  echarge = 0.0;
  for (i = 0; i < Ncv; ++i)
    {
      echarge += fill[i];
      Z = Zion - echarge + 1.0;
      eigenvalue[i] = -0.5 * (Z * Z) / ((double) (n[i] * n[i]));
      if (eigenvalue[i] > Vall[Ngrid - 1])
	eigenvalue[i] = 2.0 * Vall[Ngrid - 1];
    }

  E_Hartree = 0.0;
  P_Hartree = 0.0;
  E_exchange = 0.0;
  E_correlation = 0.0;
  P_exchange = 0.0;
  P_correlation = 0.0;


  it = 0;
  converged = False;
  while ((!converged) && (it < Max_Iterations))
    {

      it = it + 1;
      solver_iterations = it;
      converged = True;

      /*************************************/
      /* solve for each of the eigenstates */
      /*************************************/
      for (i = 0; i < Ncv; ++i)
	{

	  Etmp = eigenvalue[i];
	 /**********************************************************/
	  /*          solve radial equation(s) for state i          */
	 /**********************************************************/
	  if (Solver_Type == Schrodinger)
	    {
	      R_Schrodinger(n[i], l[i], Vall,
		            &turning_point[i], &Etmp, r_psi[i],
		            r_psi_prime[i]);
	    }
	  else if (Solver_Type == Pauli)
	    {
	      R_Pauli(n[i], l[i], Zion, Vall,
		      &turning_point[i], &Etmp, r_psi[i], r_psi_prime[i]);
	    }
	  else if (Solver_Type == Dirac)
	    {
	      R_Dirac(n[i], l[i], s2[i], Zion, Vall,
		      &turning_point[i], &Etmp, r_psi[i], r_psi_prime[i]);
	    }
	  else if (Solver_Type == ZORA)
	    {
	      R_ZORA(n[i], l[i], s2[i], Zion, Vall,
		      &turning_point[i], &Etmp, r_psi[i], r_psi_prime[i]);
	    }
	 /**********************************************************/

	  if (fabs (eigenvalue[i] - Etmp) > 1.e-15)
	    converged = False;
	  eigenvalue[i] = Etmp;
	}			/* solving eigenstates */

      /********************************/
      /* caluculate the total density */
      /********************************/
      Zero_LogGrid (rho);
      for (i = 0; i < Ncv; ++i)
	{
	  for (k = 0; k < Ngrid; ++k)
	    rho[k] += fill[i] * pow ((r_psi[i][k] / rgrid[k]), 2.0);
	}


      /***************************************/
      /* get the hartree potential an energy */
      /***************************************/
      P_Hartree = R_Hartree_DFT(rho, echarge, Vh);
      E_Hartree = 0.5 * P_Hartree;

      /*****************************************/
      /* get the exchange potential and energy */
      /*****************************************/
      R_Exchange_DFT(rho, Vx, &E_exchange, &P_exchange);

      /********************************************/
      /* get the correlation potential and energy */
      /********************************************/
      R_Correlation_DFT(rho, Vc, &E_correlation, &P_correlation);


      /**********************************************************/
      /* add up all the scf potentials on to the -Z/r potential */
      /**********************************************************/
      for (k = 0; k < Ngrid; ++k)
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
	  for (k = 0; k < Ngrid; ++k)
	    {
	      rl0 = Vo[k] - Vall[k];
	      rl1 = Vo1[k] - Vi1[k];
	      dr = rl0 - rl1;
	      sn = sn + rl0 * dr * (rgrid[k] * rgrid[k]);
	      sd = sd + dr * dr * (rgrid[k] * rgrid[k]);
	    }			/* for k */
	  thl = sn / sd;
	}			/* if it>1 */
      for (k = 0; k < Ngrid; ++k)
	{
	  vn = (1.0 - 0.5) * ((1.0 - thl) * Vall[k] + thl * Vi1[k])
	    + 0.5 * ((1.0 - thl) * Vo[k] + thl * Vo1[k]);
	  Vi1[k] = Vall[k];
	  Vo1[k] = Vo[k];
	  Vall[k] = vn;
	}			/* for k */


    }				/* end while */


  if (!converged)
    printf ("Not Converged\n");


   /*************************/
  /* find the total energy */
   /*************************/
  Total_E = 0.0;
  for (i = 0; i < Ncv; ++i)
    {
      Total_E += fill[i] * eigenvalue[i];
    }

  Total_E += E_Hartree - P_Hartree
    + E_exchange - P_exchange + E_correlation - P_correlation;


   /********************************/
  /* caluculate the core density  */
   /********************************/
  Zero_LogGrid (rho_core);
  for (i = 0; i < Ncore; ++i)
    {
      for (k = 0; k < Ngrid; ++k)
	rho_core[k] += fill[i] * pow ((r_psi[i][k] / rgrid[k]), 2.0);
    }				/*for i */

   /********************************/
  /* caluculate the valence density  */
   /********************************/
  Zero_LogGrid (rho_valence);
  for (i = Ncore; i < Ncv; ++i)
    {
      for (k = 0; k < Ngrid; ++k)
	rho_valence[k] += fill[i] * pow ((r_psi[i][k] / rgrid[k]), 2.0);
    }				/*for i */


   /****************************/
  /* free up temporary memory */
   /****************************/
  dealloc_LogGrid (Vh);
  dealloc_LogGrid (Vx);
  dealloc_LogGrid (Vc);
  dealloc_LogGrid (Vo);
  dealloc_LogGrid (Vo1);
  dealloc_LogGrid (Vi1);

   /************************************************/
  /* find the outermost peak of the wavefunctions */
   /************************************************/
  for (i = 0; i < Ncv; ++i)
    {
      k = Ngrid - 2;
      while ((r_psi_prime[i][k] * r_psi_prime[i][k + 1] >= 0.0) && (k >= 0))
	--k;
      peak[i] = rgrid[k];
    }

}				/* solve_Atom */


/********************************
 *				*
 * solve_Scattering_State_Atom  *
 *				*
 ********************************/

void solve_Scattering_State_Atom (int nt, int lt, double et, double rmax)
{
  double r0, al;
  int st;

  r0 = r_LogGrid ()[0];
  al = log_amesh_LogGrid ();

  n[Ncv] = nt;
  l[Ncv] = lt;
  eigenvalue[Ncv] = et;
  fill[Ncv] = 0.0;
  turning_point[Ncv] = rint (log (rmax / r0) / al);
  peak[Ncv] = rmax;

  if (Solver_Type == Schrodinger)
    {
      R_Schrodinger_Fixed_E(nt, lt, Vall,
			     turning_point[Ncv], et,
			     r_psi[Ncv], r_psi_prime[Ncv]);
    }
  else if (Solver_Type == Pauli)
    {
      R_Pauli_Fixed_E(nt, lt, Zion, Vall,
		       turning_point[Ncv], et, r_psi[Ncv], r_psi_prime[Ncv]);
    }
  else if (Solver_Type == Dirac)
    {
      s2[Ncv] = 1;
      s2[Ncv + 1] = -1;
      n[Ncv]=n[Ncv+1]=nt;
      l[Ncv]=l[Ncv+1]=lt;
      eigenvalue[Ncv]=et;
      eigenvalue[Ncv + 1] = et;
      fill[Ncv]=fill[Ncv + 1] = 0.0;
      turning_point[Ncv]=turning_point[Ncv + 1] = rint (log (rmax / r0) / al);
      peak[Ncv]=peak[Ncv + 1] = rmax;
      R_Dirac_Fixed_E (nt, lt, 1, Zion, Vall,
		       turning_point[Ncv], et, r_psi[Ncv], r_psi_prime[Ncv]);
      R_Dirac_Fixed_E (nt, lt, -1, Zion, Vall,
		       turning_point[Ncv + 1], et, r_psi[Ncv + 1],
		       r_psi_prime[Ncv + 1]);
    }
  else if (Solver_Type == ZORA)
    {
      s2[Ncv] = 1;
      s2[Ncv + 1] = -1;
      n[Ncv]=n[Ncv+1]=nt;
      l[Ncv]=l[Ncv+1]=lt;
      eigenvalue[Ncv]=et;
      eigenvalue[Ncv + 1] = et;
      fill[Ncv]=fill[Ncv + 1] = 0.0;
      turning_point[Ncv]=turning_point[Ncv + 1] = rint (log (rmax / r0) / al);
      peak[Ncv]=peak[Ncv + 1] = rmax;
      R_ZORA_Fixed_E (nt, lt, 1, Zion, Vall,
		       turning_point[Ncv], et, r_psi[Ncv], r_psi_prime[Ncv]);
      R_ZORA_Fixed_E (nt, lt, -1, Zion, Vall,
		       turning_point[Ncv + 1], et, r_psi[Ncv + 1],
		       r_psi_prime[Ncv + 1]);
    }

}				/* solve_Scattering_State */


/********************************
 *				*
 *        print_Atom		*
 *				*
 ********************************/

char *
solver_Name_Atom ()
{
  char *s;
  if (Solver_Type == Schrodinger)
    s = "Schrodinger";
  else if (Solver_Type == Pauli)
    s = "Pauli";
  else if (Solver_Type == Dirac)
    s = "Dirac";
  else if (Solver_Type == ZORA)
    s = "ZORA";
  else
    s = "unknown?";

  return s;
}

int
solver_Type_Atom ()
{
  return Solver_Type;
}

void
print_Atom (FILE * fp)
{
  int i, st;

  fprintf (fp, "All electron atom solver\n\n");
  fprintf (fp, "Atom name: %s\n", atom_name);
  fprintf (fp, "Zcharge  : %le\n", Zion);
  fprintf (fp, "Amass    : %le\n", amass);
  fprintf (fp, "Ncore    : %d\n", Ncore);
  fprintf (fp, "Nvalence : %d\n", Nvalence);
  fprintf (fp, "         : restricted calculation\n");

  fprintf (fp, "\n------ Solver information ------\n");
  fprintf (fp, "solver type      : %s\n", solver_Name_Atom ());
  fprintf (fp, "hartree type     : %s\n", hartree_Name_DFT ());
  fprintf (fp, "exchange type    : %s\n", exchange_Name_DFT ());
  if (strcmp (exchange_Name_DFT (), "Dirac") == 0)
    fprintf (fp, "           alpha : %lf\n", Dirac_alpha ());
  fprintf (fp, "correlation type : %s\n", correlation_Name_DFT ());
  fprintf (fp, "Solver iterations: %d\n", solver_iterations);


  fprintf (fp, "\n------- Grid information -------\n");
  fprintf (fp, "Zcharge  : %le\n", Zion);
  fprintf (fp, "amesh    : %le\n", amesh_LogGrid ());
  fprintf (fp, "R_max    : %le\n", r_LogGrid ()[N_LogGrid () - 1]);
  fprintf (fp, "Ngrid    : %d\n", N_LogGrid ());
  fprintf (fp,
	   "------------------------------------------------------------\n");
  if (Solver_Type == Dirac || Solver_Type==ZORA)
    {
      fprintf (fp, "n\tl\tspin\tpopulation\tEigenvalue\tOuter Peak\n");
      for (i = 0; i < Ncv; ++i)
	{
	  fprintf (fp, "%d\t%s\t%.1lf\t%.2lf\t\t%le\t%le\n", n[i],
		   spd_Name (l[i]), 0.5 * s2[i], fill[i], eigenvalue[i],
		   peak[i]);
	}
    }
  else
    {
      fprintf (fp, "n\tl\tpopulation\tEigenvalue\tOuter Peak\n");
      for (i = 0; i < Ncv; ++i)
	{

	  fprintf (fp, "%d\t%s\t%.2lf\t\t%le\t%le\n", n[i], spd_Name (l[i]),
		   fill[i], eigenvalue[i], peak[i]);
	}
    }
  fprintf (fp,
	   "------------------------------------------------------------\n");

  fprintf (fp, "electronic charge      = %le\n", -Integrate_LogGrid (rho));
  fprintf (fp, "electronic core charge = %le\n",
	   -Integrate_LogGrid (rho_core));
  fprintf (fp, "total atom charge      = %le\n",
	   Zion - Integrate_LogGrid (rho));




  fprintf (fp, "\nTotal E       = %le\n", Total_E);
  fprintf (fp, "\n");


  fprintf (fp, "E_Hartree     = %le\n", E_Hartree);
  fprintf (fp, "<Vh>          = %le\n", P_Hartree);
  fprintf (fp, "\n");

  fprintf (fp, "E_exchange    = %le\n", E_exchange);
  fprintf (fp, "<Vx>          = %le\n", P_exchange);
  fprintf (fp, "\n");

  fprintf (fp, "E_correlation = %le\n", E_correlation);
  fprintf (fp, "<Vc>          = %le\n", P_correlation);

}				/* print_Atom */

/********************************
 *				*
 * set_(solver parameters)_Atom	*
 *				*
 ********************************/

void
set_Solver_Atom (solver)
     int solver;
{
  Solver_Type = solver;
}


/********************************
 *				*
 *	      E_Atom 		*
 *				*
 ********************************/

double E_Atom()
{
  return Total_E;
}				/*E_Atom */

double eigenvalue_Atom(int i)
{
  return eigenvalue[i];
}

double * rho_Atom()
{
  return rho;
}

double * rho_core_Atom()
{
  return rho_core;
}

double * rho_valence_Atom()
{
  return rho_valence;
}


double * Vall_Atom()
{
  return Vall;
}


double * r_psi_Atom(int i)
{
  return r_psi[i];
}

double * r_psi_prime_Atom(int i)
{
  return r_psi_prime[i];
}

int n_Atom(int i)
{
  return n[i];
}

int l_Atom(int i)
{
  return l[i];
}

int s_Atom(int i)
{
  return s2[i];
};

int lmax_Atom()
{
  return lmax;
}

double fill_Atom(int i)
{
  return fill[i];
}

int Ncore_Atom()
{
  return Ncore;
}

int Nvalence_Atom()
{
  return Nvalence;
}

double peak_Atom(int i)
{
  return peak[i];
}

int
turning_point_Atom(int i)
{
  return turning_point[i];
}

double
Zion_Atom ()
{
  return Zion;
}

double
Amass_Atom ()
{
  return amass;
}

int
state_Atom(int nt, int lt)
{
  int i;
  i = 0;
  while (((nt != n[i]) || (lt != l[i])) && (i <= Ncv))
    ++i;

  /* Error */
  if (i > Ncv)
    printf ("Error: state_Atom\n");

  return i;
}

int
state_RelAtom(int nt, int lt, int st)
{
  int i;
  i = 0;
  while (((nt != n[i]) || (lt != l[i]) || (st != s2[i])) && (i <= (Ncv + 1)))
    ++i;

  /* Error */
  if (i > Ncv)
    printf ("Error: state_RelAtom\n");

  return i;
}

char *
name_Atom ()
{
  return atom_name;
}

char *
spin_Name (int i)
{
  char *u = "U";
  char *d = "D";
  if (s2[i] > 0)
    return u;
  return d;
}

int
isRelativistic_Atom ()
{
  if (Solver_Type == Dirac || Solver_Type==ZORA)
    return 1;
  return 0;
}
