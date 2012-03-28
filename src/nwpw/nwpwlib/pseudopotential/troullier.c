/* troullier.c -
   Author - Eric Bylaska
   $Id$

*/

#include	<stdio.h>
#ifndef WIN32
#include	<strings.h>
#endif
#include	"pred_cor.h"
#include	"loggrid.h"
#include	"schrodin.h"
#include	"dft.h"
#include	"atom.h"
#include	"xpansion.h"
#include	"troullier.h"
#include	"debug.h"

#define	Max(x,y)	((x>y) ? x : y)
#define	True	1
#define	False	0
#define SMALL	1.0e-9

/********************************
 *				*
 *    Suggested_Param_Troullier *
 *				*
 ********************************/

/*  This routine returns suggested parameters
values for the Troullier prescription.

   Entry -
   Exit	 - num_states_psp
	   n_psp[],
	   l_psp[],
           e_psp[], an array of suggested psp eigenvalues.
           fill_psp[], an array of suggested psp fillings
  	   rcut_psp[], an array of rcut values

   Uses - Atom data structure

*/

void	Suggested_Param_Troullier(num_states_psp,n_psp,l_psp,
                               e_psp,fill_psp,rcut_psp)

int	*num_states_psp;
int	n_psp[];
int	l_psp[];
double	e_psp[];
double	fill_psp[];
double	rcut_psp[];
{
    int p,npsps;
    int i,l,lmax;
    int Nc,Nv,n;
    double rcmax,emax;

    Nc	= Ncore_Atom();
    Nv   = Nvalence_Atom();
    lmax = lmax_Atom();
    npsps = lmax+2;

    for (p=0; p<npsps; ++p)
    {
        rcut_psp[p] = 0.0;
        fill_psp[p] = 0.0;
        n_psp[p] = 0;
        l_psp[p] = p;
    }
    emax = 0.0;

    /*******************************************/
    /* iterate over core states                */
    /* - all core states are scattering states */
    /*******************************************/
    rcmax = 0.0;
    for (i=0; i<Nc; ++i)
    {
        n = n_Atom(i);
        l = l_Atom(i);

        /* lowest l state, i.e. 1s, 2p, 3d, ... */
        rcut_psp[l] = peak_Atom(i);
        rcmax       = Max(rcmax,rcut_psp[l]);
        fill_psp[l] = 0.0;  /* all core states are scattering states */
        n_psp[l]    = Max(n_psp[l],n);
        e_psp[l]    = 0.0;
    } /* core states */

    /***********************************/
    /* iterate over valence states     */
    /* - remove core scattering states */
    /***********************************/
    if (Nv > 0)
    {
        rcmax = 0.0;
        emax  = -100.0;
        for (i=Nc; i<(Nc+Nv); ++i)
        {
            n = n_Atom(i);
            l = l_Atom(i);

            /* lowest l state, i.e. 1s, 2p, 3d, ... */
            rcut_psp[l] = peak_Atom(i);
            n_psp[l]    = n;
            fill_psp[l] = fill_Atom(i);
            e_psp[l]    = eigenvalue_Atom(i);
            emax        = Max(emax,e_psp[l]);
            rcmax       = Max(rcmax,rcut_psp[l]);
        }

    } /* valence states */

    /* set n_psp for guarenteed scatttering state */
    n_psp[npsps-1] = l_psp[npsps-1]+1;

    /* set the rcut for the scattering states */
    for (p=0; p<npsps; ++p)
    {
        if (fill_psp[p] == 0.0)
        {
            rcut_psp[p]  = Max(rcut_psp[p],rcmax);
            e_psp[p]     = emax;
        }
    }


    *num_states_psp = npsps;

} /* Suggested_Params_Troullier */


/********************************
 *				*
 *         solve_Troullier      *
 *				*
 ********************************/

/*  This routine solves for the Troullier-Martins psp

   Entry - num_psp
	   n_psp[],
	   l_psp[],
           e_psp[], an array of suggested psp eigenvalues.
           fill_psp[], an array of suggested psp fillings
  	   rcut_psp[], an array of rcut values

   Uses - Atom data structure

*/

void	solve_Troullier(num_psp,n_psp,l_psp,e_psp,fill_psp,rcut_psp,
                     r_psi_psp,r_psi_prime_psp,rho_psp,rho_semicore,V_psp,
                     eall_psp,
                     eh_psp,ph_psp,
                     ex_psp,px_psp,
                     ec_psp,pc_psp,
                     kb_expansion,r_psi_extra,r_psi_prime_extra)


int	num_psp;
int	n_psp[];
int	l_psp[];
double	e_psp[];
double	fill_psp[];
double	rcut_psp[];
double	**r_psi_psp;
double	**r_psi_prime_psp;
double	*rho_psp;
double	*rho_semicore;
double	**V_psp;
double	*eall_psp;
double	*eh_psp;
double	*ph_psp;
double	*ex_psp;
double	*px_psp;
double	*ec_psp;
double	*pc_psp;
int     kb_expansion[];
double  **r_psi_extra;
double  **r_psi_prime_extra;

{
    int 		istate,
    i,l,k,
    match,mch,
    Ngrid;
    double	al,amesh,rmax,Zion;
    double	gamma,gpr,nu0;
    double	ldpsi_match;
    double	el,eeig;
    double	ph,px,pc,eh,ex,ec;

    double	*ul,*ul_prime;
    double	*wl,*wl_prime;
    double	*r,
    *Vh,*Vx,*Vc,*rho_total,
    *Vall,
    *Vl;


    /* Allocate Grids */
    Vall       = Vall_Atom();

    r     = r_LogGrid();
    Ngrid = N_LogGrid();
    al    = log_amesh_LogGrid();
    amesh = amesh_LogGrid();

   eeig   = 0.0;
    Zion  = 0.0;
    for (k=0; k<Ngrid; ++k)
        rho_psp[k] = 0.0;


    if (debug_print()){
        printf("\n\nTroullier pseudopotential check\n\n");
        printf("l\trcore     rmatch    E in       E psp      norm test slope test\n");
    }

    for (l=0; l<(num_psp); ++l)
    {
        wl        = r_psi_psp[l];
        wl_prime	= r_psi_prime_psp[l];
        Vl        = V_psp[l];


        /*******************************************/
        /* Solve for scattering state if necessary */
        /*******************************************/
        if (fill_psp[l] == 0.0)
        {
            rmax = 20.0;
            solve_Scattering_State_Atom(n_psp[l],l_psp[l],e_psp[l],rmax);

            /* scattering state saved at the end of the atom list */
            istate = Nvalence_Atom() + Ncore_Atom();

            mch       = rint(log(rmax/r[0])/al);
            ul        = r_psi_Atom(istate);
            ul_prime  = r_psi_prime_Atom(istate);
            nu0       = Norm_LogGrid(mch,(l+1.0),ul);
            nu0       = 1.0/sqrt(nu0);
            for (i=0; i<Ngrid; ++i)
            {
                ul[i]       = ul[i]*nu0;
                ul_prime[i] = ul_prime[i]*nu0;
            }
        }
        /*******************************************/
        /* find state of all-electron wavefunction */
        /*******************************************/
        else
            istate = state_Atom(n_psp[l],l_psp[l]);

        /*************************************/
        /* get the all-electron wavefunction */
        /*************************************/
        ul          = r_psi_Atom(istate);
        ul_prime    = r_psi_prime_Atom(istate);
        el          = e_psp[l];

        /*****************************/
        /* find matching point stuff */
        /*****************************/
        match            = rint(log(rcut_psp[l]/r[0])/al);
        ldpsi_match      = ul_prime[match]/ul[match];


        /* make sure that wavefunctions are non-negative at the matching point */
        if (ul[match]<0.0)
        {
            nu0 = -1.0;
            for (i=0; i<Ngrid; ++i)
            {
                ul[i] = ul[i]*nu0;
                ul_prime[i] = ul_prime[i]*nu0;
            }
        }

        /**************************************/
        /* generate troullier pseudopotential */
        /**************************************/
        init_xpansion(l,match,fill_psp[l],el,Vall,ul,ul_prime,Vl,wl,wl_prime);
        psi_xpansion();
        dpsi_xpansion();
        psp_xpansion();

        /******************/
        /* verify psp Vl */
        /******************/
        /*******************/
        /* psp bound state */
        /*******************/
        if (fill_psp[l] > 0.0)
        {

            R_Schrodinger(l_psp[l]+1,l_psp[l],Vl,&mch,&el,wl,wl_prime);

        }
        /* scattering state */
        else
        {
            R_Schrodinger_Fixed_Logderiv(l_psp[l]+1,l_psp[l],Vl,match,
                                         ldpsi_match,&el,wl,wl_prime);
            R_Schrodinger_Fixed_E(l_psp[l]+1,l_psp[l],Vl,
                                  Ngrid-1,el,wl,wl_prime);

            /* normalize the scattering state to mch = 20.0 a.u.*/
            rmax      = 20.0;
            mch       = rint(log(rmax/r[0])/al);
            nu0       = Norm_LogGrid(mch,(l+1.0),wl);
            nu0       = 1.0/sqrt(nu0);
            for (i=0; i<Ngrid; ++i)
            {
                wl[i]       = wl[i]*nu0;
                wl_prime[i] = wl_prime[i]*nu0;
            }
        }

        gamma=fabs(ul[match]/wl[match]);
        gpr  =fabs(ul_prime[match]/wl_prime[match]);
        if (debug_print()){
            /*printf("ul[match] wl[match]: %lf %lf\n",ul[match],wl[match]);*/
            printf("%d\t%lf  %lf  %lf  %lf  %lf  %lf\n",l_psp[l],
                   rcut_psp[l],r[match],
                   e_psp[l],el,
                   gamma,gpr);
        }

        /* Use the analytic form of pseudopotential */
        el = e_psp[l];
        psi_xpansion();
        dpsi_xpansion();
        psp_xpansion();

        if (fill_psp[l] > 0.0)
        {
            eeig += fill_psp[l]*el;

            /* accumulate charges */
            Zion += fill_psp[l];
            for (k=0; k<Ngrid; ++k)
                rho_psp[k] += fill_psp[l]*pow(wl[k]/r[k],2.0);
        }


    } /* for l */

   /* solve for other states */
   if (debug_print())
   {
      printf("\n\nTroullier pseudopotential extra states\n\n");
      printf("l\tn\tE psp\n");
   }
   i = 0;
   for (l=0; l<num_psp; ++l)
   {
      Vl = V_psp[l];
      el = e_psp[l]+0.1;
      for (k=0; k<(kb_expansion[l_psp[l]]-1); ++k)
      {
         wl       = r_psi_extra[i];
         wl_prime = r_psi_prime_extra[i];
         R_Schrodinger(l_psp[l]+2+k,l_psp[l],Vl,&mch,&el,wl,wl_prime);
         ++i;
         if (debug_print()) printf("%d\t%d\t%lf\n",l_psp[l],l_psp[l]+2+k,el);
      }
   }



    /***************************************/
    /* get the hartree potential an energy */
    /* get the exchange potential and energy */
    /* get the correlation potential and energy */
    /***************************************/
    Vh = alloc_LogGrid();
    Vx = alloc_LogGrid();
    Vc = alloc_LogGrid();
    ph = R_Hartree_DFT(rho_psp,Zion,Vh);
    eh = 0.5*ph;

    rho_total = alloc_LogGrid();
    for (k=0; k<Ngrid; ++k)
        rho_total[k] = rho_psp[k] + rho_semicore[k];

    R_Exchange_DFT(rho_total,Vx,&ex,&px);
    R_Correlation_DFT(rho_total,Vc,&ec,&pc);

    R_Screening_Cut(Vx);
    R_Screening_Cut(Vc);

    /* recalculate px and pc */
    for (k=0; k<Ngrid; ++k) rho_total[k] = (rho_psp[k])*Vx[k];
    px = Integrate_LogGrid(rho_total);
    for (k=0; k<Ngrid; ++k) rho_total[k] = (rho_psp[k])*Vc[k];
    pc = Integrate_LogGrid(rho_total);


    *eall_psp = eeig + eh + ex + ec - ph - px - pc;
    *eh_psp   = eh;
    *ph_psp   = ph;
    *ex_psp   = ex;
    *px_psp   = px;
    *ec_psp   = ec;
    *pc_psp   = pc;
    for (l=0; l<num_psp; ++l)
        for (k=0; k<Ngrid; ++k)
            V_psp[l][k] = V_psp[l][k] - Vh[k] - Vx[k] - Vc[k];

    /* deallocate memory */
    dealloc_LogGrid(Vh);
    dealloc_LogGrid(Vx);
    dealloc_LogGrid(Vc);
    dealloc_LogGrid(rho_total);

} /* solve_Troullier */
