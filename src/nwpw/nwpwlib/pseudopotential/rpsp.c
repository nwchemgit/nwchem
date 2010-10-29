/* psp.c -
   author - Patrick Nichols
   $Id$
*/

#include	<stdio.h>
#include        <stdlib.h>
#include	<string.h>

#include	"name.h"
#include	"get_word.h"
#include	"loggrid.h"

#include	"dft.h"
#include	"atom.h"
#include        "rtroullier.h"
#include	"rhamann.h"
#include	"troullier.h"
#include        "vanderbilt.h"
#include	"generate_rho_semicore.h"
#include	"rpsp.h"
#include        "psp.h"
#define	False	0
#define	True	1
#define	Max(x,y)	((x>y) ? x : y)


/* atom structure variables */

static int Nvalence;
static int npsp_states;
static int *n;
static int *l;
static int *spin;
static int lmax;
static double *fill;
static double *rcut;
static double *peak;
static double Zion;
static double Total_E,
  E_Hartree, P_Hartree, E_exchange, P_exchange, E_correlation, P_correlation;
static double *eigenvalue;
static double **r_psi;
static double **r_psi_prime;
static double *rho;
static double *rho_semicore;
static double *drho_semicore;
static double r_semicore;
static double **V_psp;
static char comment[80];

/* extra Vanderbilt parameters */
static double rlocal, clocal;
static int ns[10], indx_il[4][10], indx_ijl[4][4][10];
static double *Vlocal;
static double **r_hard_psi;
static double *D0;
static double *q;

/* solver parameters: this is the Kawai-Weare default */
static int Solver_Type = Hamann;


/********************************
 *				*
 *	  init_RelPsp 		*
 *				*
 ********************************/

void
init_RelPsp (char *filename)
{
  int p, p1, p2;
  int ltmp, semicore_type;
  double rctmp;
  char *w, *tc;
  FILE *fp;


    /* find the psp type */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<pseudopotential>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        w = get_word(fp);
        if (strcmp("hamann",w)==0)            Solver_Type = Hamann;
        if (strcmp("troullier-martins",w)==0) Solver_Type = Troullier;
    }else{
        Solver_Type = Hamann;
    }
    fclose(fp);



    /* find the comment */  
    if (Solver_Type == Hamann)    strcpy(comment,"Hamann pseudopotential");
    if (Solver_Type == Troullier) strcpy(comment,"Troullier-Martins pseudopotential");

  /* set lmax  */
  lmax = lmax_Atom ();

  /* set the number psp projectors */
  npsp_states = 2 * lmax + 4;
  /* allocate memory for n,l,fill,rcut,peak, and eigenvalue */
  n = (int *) malloc ((npsp_states) * sizeof (int));
  l = (int *) malloc ((npsp_states) * sizeof (int));
  spin = (int *) malloc ((npsp_states) * sizeof (int));
  fill = (double *) malloc ((npsp_states) * sizeof (double));
  rcut = (double *) malloc ((npsp_states) * sizeof (double));
  peak = (double *) malloc ((npsp_states) * sizeof (double));
  eigenvalue = (double *) malloc ((npsp_states) * sizeof (double));


  /* allocate memory for r_psi, V_psp, and rho */
  r_psi = (double **) malloc ((npsp_states) * sizeof (double *));
  r_psi_prime = (double **) malloc ((npsp_states) * sizeof (double *));
  V_psp = (double **) malloc ((npsp_states) * sizeof (double *));
  for (p = 0; p < npsp_states; ++p)
    {
      r_psi[p] = alloc_LogGrid ();
      r_psi_prime[p] = alloc_LogGrid ();
      V_psp[p] = alloc_LogGrid ();
    }
  rho = alloc_LogGrid ();
  rho_semicore = alloc_LogGrid ();
  drho_semicore = alloc_LogGrid ();


  /* get the psp info */
  if (Solver_Type==Troullier) {
	  Suggested_Param_RelTroullier(&Nvalence, n, l, spin, eigenvalue, fill, rcut);
  }else{
	  Suggested_Param_RelHamann (&Nvalence, n, l, spin, eigenvalue, fill, rcut);
  }
  /* set the number psp projectors */
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<npsp-states>", w) != 0))
    w = get_word (fp);

  if (w != NIL)
    {
      w = get_word (fp);
      while ((w != NIL) && (strcmp ("<end>", w) != 0))
	{
	  sscanf (w, "%d", &ltmp);
	  w = get_word (fp);
	  Nvalence = 2 * ltmp + 2;
	}
    }
  fclose (fp);

  /* get rcut */
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<rcut>", w) != 0))
    w = get_word (fp);

  if (w != NIL)
    {
      w = get_word (fp);
      while ((w != NIL) && (strcmp ("<end>", w) != 0))
	{
	  sscanf (w, "%d", &ltmp);
	  w = get_word (fp);
	  sscanf (w, "%lf", &rctmp);
	  w = get_word (fp);
	  rcut[2 * ltmp] = rctmp;
	  rcut[2 * ltmp + 1] = rctmp;
	}
    }
  fclose (fp);

  /* get ecut */
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<ecut>", w) != 0))
    w = get_word (fp);

  if (w != NIL)
    {
      w = get_word (fp);
      while ((w != NIL) && (strcmp ("<end>", w) != 0))
	{
	  sscanf (w, "%d", &ltmp);
	  w = get_word (fp);
	  sscanf (w, "%lf", &rctmp);
	  w = get_word (fp);
	  eigenvalue[2 * ltmp] = rctmp;
	  eigenvalue[2 * ltmp + 1] = rctmp;
	}
    }
  fclose (fp);

  /* get r_semicore - if zero then no core corrections added */
  r_semicore = 0.0;
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<semicore>", w) != 0))
    w = get_word (fp);

  if (w != NIL)
    {
      w = get_word (fp);
      while ((w != NIL) && (strcmp ("<end>", w) != 0))
	{
	  sscanf (w, "%lf", &rctmp);
	  w = get_word (fp);
	  r_semicore = rctmp;
	}
    }
  fclose (fp);

  /* find the semicore_type */
  semicore_type = 2;
  fp = fopen (filename, "r+");
  w = get_word (fp);
  while ((w != NIL) && (strcmp ("<semicore_type>", w) != 0))
    w = get_word (fp);
  if (w != NIL)
    {
      w = get_word (fp);
      if (strcmp ("quadratic", w) == 0)
	semicore_type = 0;
      if (strcmp ("louie", w) == 0)
	semicore_type = 1;
      if (strcmp ("fuchs", w) == 0)
	semicore_type = 2;
    }
  fclose (fp);




  /* generate non-zero rho_semicore */
  if (r_semicore > 0.0)
    {
      /*
         printf("\n\n");
         printf("Generating non-zero semicore density\n");
       */
      generate_rho_semicore (semicore_type, rho_core_Atom (), r_semicore,
			     rho_semicore);
      Derivative_LogGrid (rho_semicore, drho_semicore);
    }

  /* define the ion charge */
  /*
     Zion=0.0;
     for (p=Ncore_Atom(); p<(Ncore_Atom()+Nvalence_Atom()); ++p)
     Zion += fill_Atom(p);
   */

  Zion = Zion_Atom ();
  for (p = 0; p < (Ncore_Atom ()); ++p)
    Zion -= fill_Atom (p);
  
}				/* init_RelPsp */



/********************************
 *				*
 *	   solve_RelPsp 		*
 *				*
 ********************************/

void
solve_RelPsp ()
{

  int p, k, Ngrid;
  double *rgrid;

  /* get loggrid variables */
  Ngrid = N_LogGrid ();
  rgrid = r_LogGrid ();

  if (Solver_Type==Troullier) {
  solve_RelTroullier(Nvalence, n, l, spin, eigenvalue, fill, rcut,
		   r_psi, r_psi_prime, rho, rho_semicore, V_psp,
		   &Total_E,
		   &E_Hartree, &P_Hartree,
		   &E_exchange, &P_exchange, &E_correlation, &P_correlation);
  }else{
  solve_RelHamann(Nvalence, n, l, spin, eigenvalue, fill, rcut,
		   r_psi, r_psi_prime, rho, rho_semicore, V_psp,
		   &Total_E,
		   &E_Hartree, &P_Hartree,
		   &E_exchange, &P_exchange, &E_correlation, &P_correlation);
  }
   /******************************************************/
  /* find the outermost peak of the pseudowavefunctions */
   /******************************************************/
  for (p = 0; p < Nvalence; ++p)
    {
      if (fill[p] != 0.0)
	{
	  k = Ngrid - 2;
	  while ((r_psi_prime[p][k] * r_psi_prime[p][k + 1] >= 0.0)
		 && (k >= 0))
	    --k;
	  peak[p] = rgrid[k];
	}
      else
	peak[p] = rgrid[Ngrid - 1];
    }
}				/* solve_RelPsp */



/********************************
 *				*
 *        print_RelPsp		*
 *				*
 ********************************/

char *
solver_Name_RelPsp ()
{
  char *s;
  if (Solver_Type==Troullier) {
	s="TroullierMartins";
  }else{
        s = "Hamann";
  }
  return s;
}


void
print_RelPsp (FILE * fp)
{
  int i;
  double dx;
  fprintf (fp, "\n\n");
  fprintf (fp, "PSP solver information\n\n");
  fprintf (fp, "Atom name: %s\n", name_Atom ());
  fprintf (fp, "Zcharge  : %le\n", Zion);
  fprintf (fp, "Nvalence : %d\n", Nvalence);
  fprintf (fp, "         : restricted calculation\n");

  fprintf (fp, "\n------ Solver information ------\n");
  fprintf (fp, "solver type      : %s\n", solver_Name_RelPsp ());
  fprintf (fp, "hartree type     : %s\n", hartree_Name_DFT ());
  fprintf (fp, "exchange type    : %s\n", exchange_Name_DFT ());
  if (strcmp (exchange_Name_DFT (), "Dirac") == 0)
    fprintf (fp, "           alpha : %lf\n", Dirac_alpha ());
  fprintf (fp, "correlation type : %s\n", correlation_Name_DFT ());

  fprintf (fp,
	   "----------------------------------------------------------------------------\n");
  fprintf (fp, "n\tl\ts\tpopulation\tEcut\t\tRcut\t\tOuter Peak\n");

  for (i = 0; i < Nvalence; ++i)
    {
      fprintf (fp, "%d\t%s\t%.1lf\t%.2lf\t\t%le\t%le\t%le\n", (l[i] + 1),
	       spd_Name (l[i]), 0.5 * spin[i], fill[i], eigenvalue[i],
	       rcut[i], peak[i]);
    }
  fprintf (fp,
	   "----------------------------------------------------------------------------\n");
  if (r_semicore > 0.0)
    {
      fprintf (fp, "SemiCore Corrections Added\n");
      fprintf (fp, "    rcore                      : %lf\n", r_semicore);
      fprintf (fp, "    Semicore Charge            : %lf\n",
	       Integrate_LogGrid (rho_semicore));
      fprintf (fp, "    Semicore Charge gradient   : %lf\n",
	       Integrate_LogGrid (drho_semicore));
    }

  fprintf (fp, "Pseudopotential ion charge       = %lf\n", Zion);
  dx= - Integrate_LogGrid(rho);
  fprintf (fp, "Pseudopotential electronic charge= %lf\n",
	   dx);
  fprintf (fp, "Pseudopotential atom charge      = %lf\n",
	   (Zion+dx));

  fprintf (fp, "\nTotal E       = %le\n", Total_E);
  fprintf (fp, "\n");


  fprintf (fp, "E_Hartree     = %le\n", E_Hartree);
  fprintf (fp, "<Vh>          = %le\n", P_Hartree);
  fprintf (fp, "\n");

  fprintf (fp, "E_exchange    = %le\n", E_exchange);
  fprintf (fp, "<Vx>          = %le\n", P_exchange);
  fprintf (fp, "\n");

  fprintf (fp, "E_correlation = %le\n", E_correlation);
  fprintf (fp, "<Vc>          = %le\n\n", P_correlation);


}				/* print_Atom */

/********************************
 *				*
 * set_(solver parameters)_Atom	*
 *				*
 ********************************/

void
set_Solver_RelPsp (solver)
     int solver;
{
  Solver_Type=solver; 
  fprintf (stdout," RelPsp:: Only Hamann or TM Solver is available\n");
  fprintf (stdout, "RelPsp::Cannot really set the solver yet!\n");
}


int
Vanderbilt_RelPsp ()
{
  int value;

  value = 0;
  if (Solver_Type == Vanderbilt)
    value = 1;
  return value;
}

int
NormConserving_RelPsp ()
{
  int value;

  value = 0;
  if (Solver_Type == Hamann)
    value = 1;
  if (Solver_Type == Troullier)
    value = 1;
  return value;
}

/********************************
 *				*
 *	      E_Atom 		*
 *				*
 ********************************/

double
E_RelPsp ()
{
  return Total_E;
}				/*E_Atom */


double
eigenvalue_RelPsp (int i)
{
  return eigenvalue[i];
}

double *
rho_RelPsp ()
{
  return rho;
}

double *
rho_semicore_RelPsp ()
{
  return rho_semicore;
}

double *
drho_semicore_RelPsp ()
{
  return drho_semicore;
}

double
r_semicore_RelPsp ()
{
  return r_semicore;
}


double *
Beta_RelPsp (int i, int l)
{
  return V_psp[indx_il[i][l]];
}

double *
r_psi_il_RelPsp (int i, int l)
{
  return r_psi[indx_il[i][l]];
}

double *
r_hard_psi_il_RelPsp (int i, int l)
{
  return r_hard_psi[indx_il[i][l]];
}

int
ns_RelPsp (int l)
{
  return ns[l];
}

double
D0_RelPsp (int i, int j, int l)
{
  return D0[indx_ijl[i][j][l]];
}

double
q_RelPsp (int i, int j, int l)
{
  return q[indx_ijl[i][j][l]];
}

double *
Vlocal_RelPsp ()
{
  return Vlocal;
}



double *
V_RelPsp (int i)
{
  return V_psp[i];
}


double *
r_psi_RelPsp (int i)
{
  return r_psi[i];
}

int
n_RelPsp (int i)
{
  return n[i];
}

int
l_RelPsp (int i)
{
  return l[i];
}

int
lmax_RelPsp ()
{
  return lmax;
}

double
fill_RelPsp (int i)
{
  return fill[i];
}


int
Nvalence_RelPsp ()
{
  return Nvalence;
}

double
peak_RelPsp (int i)
{
  return peak[i];
}

double
rcut_RelPsp (int l)
{
  return rcut[2*l];
}


double
rcut_il_RelPsp (int i, int l)
{
  return rcut[indx_il[i][l]];
}


double
Zion_RelPsp ()
{
  return Zion;
}

int
state_RelPsp (int nt, int lt, int st)
{
  int i;
  i = 0;
  while (((nt != n[i]) || (lt != l[i]) || st!=spin[i]) && (i < Nvalence))
    ++i;

  /* Error */
  if (i >= Nvalence)
    printf ("Error: state_RelPsp\n");

  return i;
}


char *
comment_RelPsp ()
{
  return comment;
}
