/* psp.c -
   author - Eric Bylaska
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
#include	"hamann.h"
#include	"troullier.h"
#include        "vanderbilt.h"
#include	"generate_rho_semicore.h"
#include	"psp.h"

#define	False	0
#define	True	1
#define	Max(x,y)	((x>y) ? x : y)


/* atom structure variables */

static	int	Nvalence;
static  int     npsp_states;
static	int	*n;
static	int	*l;
static	int	lmax;
static 	double	*fill;
static 	double	*rcut;
static  double	*peak;
static	double	Zion;
static	double	Total_E,
E_Hartree, P_Hartree,
E_exchange,
P_exchange,
E_correlation,
P_correlation;
static	double  *eigenvalue;
static	double	**r_psi;
static	double	**r_psi_prime;
static  double  *rho;
static  double  *rho_semicore;
static  double  r_semicore;
static	double	**V_psp;
static  char    comment[80];

/* extra Vanderbilt parameters */
static double rlocal,clocal;
static int    ns[10],indx_il[4][10],indx_ijl[4][4][10];
static double *Vlocal;
static double **r_hard_psi;
static double *D0;
static double *q;

/* solver parameters: this is the Kawai-Weare default */
static	int	Solver_Type      = Hamann;


/********************************
 *				*
 *	  init_Psp 		*
 *				*
 ********************************/

void	init_Psp(char *filename)
{
    int	  p,p1,p2;
    int	  ltmp;
    double rctmp;
    char   *w;
    FILE	  *fp;



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
        if (strcmp("vanderbilt",w)==0)        Solver_Type = Vanderbilt;
    }
    fclose(fp);



    /* find the comment */
    if (Solver_Type = Hamann)     comment = "Hamann pseudopotential";
    if (Solver_Type = Troullier)  comment = "Troullier-Martins pseudopotential";
    if (Solver_Type = Vanderbilt) comment = "Vanderbilt pseudopotential";
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<comment>",w)!=0))
        w = get_word(fp);
    if (w!=NIL)
    {
        fscanf(fp,"%s",comment);
    }
    fclose(fp);


    /* set lmax  */
    lmax = lmax_Atom() + 1;

    /* set the number psp projectors */
    npsp_states = lmax + 1;
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<npsp-states>",w)!=0))
        w = get_word(fp);

    if (w!=NIL)
    {
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<end>",w) != 0))
        {
            sscanf(w,"%d", &ltmp);
            w = get_word(fp);
            npsp_states = ltmp;
        }
    }
    fclose(fp);


    /* allocate memory for n,l,fill,rcut,peak, and eigenvalue */
    n    	 = (int *)    malloc((npsp_states)*sizeof(int));
    l    	 = (int *)    malloc((npsp_states)*sizeof(int));
    fill 	 = (double *) malloc((npsp_states)*sizeof(double));
    rcut 	 = (double *) malloc((npsp_states)*sizeof(double));
    peak 	 = (double *) malloc((npsp_states)*sizeof(double));
    eigenvalue = (double *) malloc((npsp_states)*sizeof(double));


    /* allocate memory for r_psi, V_psp, and rho */
    r_psi       = (double **) malloc((npsp_states)*sizeof(double*));
    r_psi_prime = (double **) malloc((npsp_states)*sizeof(double*));
    V_psp       = (double **) malloc((npsp_states)*sizeof(double*));
    for (p=0; p<npsp_states; ++p)
    {
        r_psi[p]       = alloc_LogGrid();
        r_psi_prime[p] = alloc_LogGrid();
        V_psp[p]       = alloc_LogGrid();
    }
    rho          = alloc_LogGrid();
    rho_semicore = alloc_LogGrid();

    /* allocate extra memory for Vanderbilt psp */
    if (Solver_Type==Vanderbilt)
    {
        r_hard_psi = (double **) malloc((npsp_states)*sizeof(double*));
        for (p=0; p<npsp_states; ++p)
        {
            r_hard_psi[p]       = alloc_LogGrid();
        }
        Vlocal = alloc_LogGrid();
        D0     = (double *) malloc((npsp_states*npsp_states)*sizeof(double));
        q      = (double *) malloc((npsp_states*npsp_states)*sizeof(double));

    } /* Solver_Type==Vanderbilt */


    /* get the psp info */
    if (Solver_Type==Hamann)
        Suggested_Param_Hamann(&Nvalence,n,l,eigenvalue,fill,rcut);
    else if (Solver_Type==Troullier)
        Suggested_Param_Troullier(&Nvalence,n,l,eigenvalue,fill,rcut);
    else if (Solver_Type==Vanderbilt)
        Suggested_Param_Vanderbilt(&Nvalence,n,l,eigenvalue,fill,rcut,
                                   &rlocal,&clocal);

    /* set the number psp projectors */
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<npsp-states>",w)!=0))
        w = get_word(fp);

    if (w!=NIL)
    {
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<end>",w) != 0))
        {
            sscanf(w,"%d", &ltmp);
            w = get_word(fp);
            Nvalence = ltmp;
        }
    }
    fclose(fp);

    /* Read in Vanderbilt parameters */
    if (Solver_Type==Vanderbilt)
    {

        /* get vanderbilt-local */
        fp = fopen(filename,"r+");
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<vanderbilt-local>",w)!=0))
            w = get_word(fp);

        if (w!=NIL)
        {
            w = get_word(fp);
            while ((w!=NIL) && (strcmp("<end>",w) != 0))
            {
                sscanf(w,"%lf",&rctmp);
                w = get_word(fp);
                rlocal = rctmp;

                sscanf(w,"%lf",&rctmp);
                w = get_word(fp);
                clocal = rctmp;

            }
        }
        fclose(fp);


        /* get vanderbilt-states */
        fp = fopen(filename,"r+");
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<vanderbilt-states>",w)!=0))
            w = get_word(fp);

        if (w!=NIL)
        {
            w = get_word(fp);
            while ((w!=NIL) && (strcmp("<end>",w) != 0))
            {
                sscanf(w,"%d", &ltmp);
                w = get_word(fp);
                sscanf(w,"%d", &p1);
                w = get_word(fp);
                sscanf(w,"%d",&p2);
                w = get_word(fp);

                n[ltmp] = p1;
                l[ltmp] = p2;
                lmax = Max(lmax,p2);

                sscanf(w,"%lf",&rctmp);
                w = get_word(fp);
                eigenvalue[ltmp] = rctmp;

                sscanf(w,"%lf",&rctmp);
                w = get_word(fp);
                rcut[ltmp] = rctmp;
                fill[ltmp] = 0.0;

            }
        }
        fclose(fp);

        /* define ns */
        for (p=0; p<=lmax; ++p)
        {
            ns[p] = 0;
            for (p1=0; p1<npsp_states; ++p1)
                if (l[p1]==p) ns[p] = ns[p] + 1;
        }
        /* define indx_il */
        p1 = 0;
        for (p=0; p<=lmax; ++p)
        {
            for (p2=0; p2 < ns[p]; ++p2)
            {
                if (p1 >= npsp_states)
                {
                    printf("Error in defining indx_il: need more npsp_states\n");
                    exit(99);
                }
                indx_il[p2][p] = p1;
                ++p1;
            }
        }

        /* define indx_ijl */
        ltmp = 0;
        for (p =0; p <=lmax; ++p)
            for (p1=0; p1<ns[p]; ++p1)
                for (p2=0; p2<ns[p]; ++p2)
                {
                    indx_ijl[p2][p1][p] = ltmp;
                    ++ltmp;
                }

    }

    /* Read in norm-conserving parameters */
    else
    {

        /* get rcut */
        fp = fopen(filename,"r+");
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<rcut>",w)!=0))
            w = get_word(fp);

        if (w!=NIL)
        {
            w = get_word(fp);
            while ((w!=NIL) && (strcmp("<end>",w) != 0))
            {
                sscanf(w,"%d", &ltmp);
                w = get_word(fp);
                sscanf(w,"%lf", &rctmp);
                w = get_word(fp);
                rcut[ltmp] = rctmp;
            }
        }
        fclose(fp);

        /* get ecut */
        fp = fopen(filename,"r+");
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<ecut>",w)!=0))
            w = get_word(fp);

        if (w!=NIL)
        {
            w = get_word(fp);
            while ((w!=NIL) && (strcmp("<end>",w) != 0))
            {
                sscanf(w,"%d", &ltmp);
                w = get_word(fp);
                sscanf(w,"%lf", &rctmp);
                w = get_word(fp);
                eigenvalue[ltmp] = rctmp;
            }
        }
        fclose(fp);

    }


    /* get r_semicore - if zero then no core corrections added */
    r_semicore = 0.0;
    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<semicore>",w)!=0))
        w = get_word(fp);

    if (w!=NIL)
    {
        w = get_word(fp);
        while ((w!=NIL) && (strcmp("<end>",w) != 0))
        {
            sscanf(w,"%lf", &rctmp);
            w = get_word(fp);
            r_semicore = rctmp;
        }
    }
    fclose(fp);

    /* generate non-zero rho_semicore */
    if (r_semicore > 0.0)
    {
        printf("\n\n");
        printf("Generating non-zero semicore density\n");
        generate_rho_semicore(rho_core_Atom(),r_semicore,rho_semicore);
    }

    /* define the ion charge */
    Zion=0.0;
    for (p=Ncore_Atom(); p<(Ncore_Atom()+Nvalence_Atom()); ++p)
        Zion += fill_Atom(p);


} /* init_Psp */



/********************************
 *				*
 *	   solve_Psp 		*
 *				*
 ********************************/

void	solve_Psp()
{

    int 		p,k,Ngrid;
    double	*rgrid;

    /* get loggrid variables */
    Ngrid = N_LogGrid();
    rgrid = r_LogGrid();

    if (Solver_Type == Hamann)
        solve_Hamann(Nvalence,n,l,eigenvalue,fill,rcut,
                     r_psi,r_psi_prime,rho,rho_semicore,V_psp,
                     &Total_E,
                     &E_Hartree,&P_Hartree,
                     &E_exchange,&P_exchange,
                     &E_correlation,&P_correlation);
    else if (Solver_Type == Troullier)
        solve_Troullier(Nvalence,n,l,eigenvalue,fill,rcut,
                        r_psi,r_psi_prime,rho,rho_semicore,V_psp,
                        &Total_E,
                        &E_Hartree,&P_Hartree,
                        &E_exchange,&P_exchange,
                        &E_correlation,&P_correlation);
    else if (Solver_Type == Vanderbilt)
        solve_Vanderbilt(Nvalence,n,l,eigenvalue,fill,rcut,rlocal,clocal,
                         ns,indx_il,indx_ijl,
                         r_hard_psi,r_psi,r_psi_prime,rho,rho_semicore,V_psp,
                         Vlocal,D0,q,
                         &Total_E,
                         &E_Hartree,&P_Hartree,
                         &E_exchange,&P_exchange,
                         &E_correlation,&P_correlation);


    /******************************************************/
    /* find the outermost peak of the pseudowavefunctions */
    /******************************************************/
    for (p=0; p<Nvalence; ++p)
    {
        if (fill[p] != 0.0)
        {
            k = Ngrid-2;
            while ((r_psi_prime[p][k]*r_psi_prime[p][k+1] >= 0.0) && (k>=0))
                --k;
            peak[p] = rgrid[k];
        }
        else
            peak[p] = rgrid[Ngrid-1];
    }


} /* solve_Psp */



/********************************
 *				*
 *        print_Psp		*
 *				*
 ********************************/

char	*solver_Name_Psp()
{
    char *s;
    if      (Solver_Type==Hamann)
        s = "Hamann";
    else if (Solver_Type==Troullier)
        s = "Troullier-Martins";
    else if (Solver_Type==Vanderbilt)
        s = "Vanderbilt";
    else
        s = "unknown?";

    return s;
}


void	print_Psp(fp)

FILE 	*fp;
{
    int i;

    fprintf(fp,"\n\n");
    fprintf(fp,"PSP solver information\n\n");
    fprintf(fp,"Atom name: %s\n",name_Atom());
    fprintf(fp,"Zcharge  : %le\n",Zion);
    fprintf(fp,"Nvalence : %d\n",Nvalence);
    fprintf(fp,"         : restricted calculation\n");

    fprintf(fp,"\n------ Solver information ------\n");
    fprintf(fp,"solver type      : %s\n",solver_Name_Psp());
    fprintf(fp,"hartree type     : %s\n",hartree_Name_DFT());
    fprintf(fp,"exchange type    : %s\n",exchange_Name_DFT());
    if (strcmp(exchange_Name_DFT(),"Dirac") == 0)
        fprintf(fp,"           alpha : %lf\n",Dirac_alpha());
    fprintf(fp,"correlation type : %s\n",correlation_Name_DFT());

    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"n\tl\tpopulation\tEcut\t\tRcut\t\tOuter Peak\n");

    for (i=0; i<Nvalence; ++i)
    {
        if (Solver_Type==Vanderbilt)
            fprintf(fp,"%d\t%s\t%s\t\t%le\t%le\t%le\n",l[i]+1,spd_Name(l[i]),
                    "  -",eigenvalue[i],rcut[i],peak[i]);
        else
            fprintf(fp,"%d\t%s\t%.2lf\t\t%le\t%le\t%le\n",l[i]+1,spd_Name(l[i]),
                    fill[i],eigenvalue[i],rcut[i],peak[i]);
    }
    fprintf(fp,"----------------------------------------------------------------------------\n");
    if (r_semicore > 0.0)
    {
        fprintf(fp,"SemiCore Corrections Added\n");
        fprintf(fp,"    rcore             : %lf\n",r_semicore);
        fprintf(fp,"    Semicore Charge   : %lf\n",
                Integrate_LogGrid(rho_semicore));
    }
    if (Solver_Type==Vanderbilt)
        fprintf(fp,"Using AE valence density for descreening\n");
    fprintf(fp,  "Pseudopotential Charge: %lf\n", Integrate_LogGrid(rho));
    fprintf(fp,"\nTotal E       = %le\n",Total_E);
    fprintf(fp,"\n");


    fprintf(fp,"E_Hartree     = %le\n",E_Hartree);
    fprintf(fp,"<Vh>          = %le\n",P_Hartree);
    fprintf(fp,"\n");

    fprintf(fp,"E_exchange    = %le\n",E_exchange);
    fprintf(fp,"<Vx>          = %le\n",P_exchange);
    fprintf(fp,"\n");

    fprintf(fp,"E_correlation = %le\n",E_correlation);
    fprintf(fp,"<Vc>          = %le\n\n",P_correlation);


} /* print_Atom */

/********************************
 *				*
 * set_(solver parameters)_Atom	*
 *				*
 ********************************/

void set_Solver_Psp(solver)
int solver;
{
    Solver_Type = solver;
}


int Vanderbilt_Psp()
{
    int value;

    value = 0;
    if (Solver_Type==Vanderbilt) value = 1;
    return value;
}

/********************************
 *				*
 *	      E_Atom 		*
 *				*
 ********************************/

double	E_Psp()
{
    return Total_E;
} /*E_Atom*/


double	eigenvalue_Psp(int i)
{
    return eigenvalue[i];
}

double	*rho_Psp()
{
    return rho;
}

double	*rho_semicore_Psp()
{
    return rho_semicore;
}

double	r_semicore_Psp()
{
    return r_semicore;
}


double *Beta_Psp(int i,int l)
{
    return V_psp[indx_il[i][l]];
}

double *r_psi_il_Psp(int i, int l)
{
    return r_psi[indx_il[i][l]];
}
double *r_hard_psi_il_Psp(int i, int l)
{
    return r_hard_psi[indx_il[i][l]];
}

int ns_Psp(int l)
{
    return ns[l];
}

double D0_Psp(int i, int j, int l)
{
    return D0[indx_ijl[i][j][l]];
}

double q_Psp(int i, int j, int l)
{
    return q[indx_ijl[i][j][l]];
}

double *Vlocal_Psp()
{
    return Vlocal;
}



double	*V_Psp(int i)
{
    return V_psp[i];
}


double	*r_psi_Psp(int i)
{
    return r_psi[i];
}

int	n_Psp(int i)
{
    return n[i];
}

int	l_Psp(int i)
{
    return l[i];
}

int	lmax_Psp()
{
    return lmax;
}

double	fill_Psp(int i)
{
    return fill[i];
}


int	Nvalence_Psp()
{
    return Nvalence;
}

double	peak_Psp(int i)
{
    return peak[i];
}

double	rcut_Psp(int i)
{
    return rcut[i];
}

double	rcut_il_Psp(int i, int l)
{
    return rcut[ indx_il[i][l] ];
}

double	Zion_Psp()
{
    return Zion;
}

int state_Psp(int nt, int lt)
{
    int i;
    i = 0;
    while ( ((nt != n[i]) || (lt != l[i])) && (i<Nvalence) )
        ++i;

    /* Error */
    if (i>=Nvalence)
        printf("Error: state_Psp\n");

    return i;
}



char    *comment_Psp()
{
    return comment;
}

