/* hamann.c -
   Author - Eric Bylaska
   $Id$

*/

#include        <stdio.h>
#include	"loggrid.h"
#include	"schrodin.h"
#include	"dft.h"
#include	"atom.h"
#include	"hamann.h"
#include	"debug.h"

#define	Max(x,y)	((x>y) ? x : y)
#define	True	1
#define	False	0
#define	SMALL	1.0e-9

/********************************
 *				*
 *    Suggested_Param_Hamann	*
 *				*
 ********************************/

/*  This routine returns suggested parameters
values for the Hamann prescription.

   Entry -
   Exit	 - num_states_psp
	   n_psp[],
	   l_psp[],
           e_psp[], an array of suggested psp eigenvalues.
           fill_psp[], an array of suggested psp fillings
  	   rcut_psp[], an array of rcut values

   Uses - Atom data structure

*/
#define	VSC	0.6
#define	VSC1	0.4
#define	CSC	1.9
#define	CSC1	4.0

void	Suggested_Param_Hamann(num_states_psp,n_psp,l_psp,
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
        if (n == (l+1))
            rcut_psp[l] = CSC1*peak_Atom(i);
        else
            rcut_psp[l] = CSC*peak_Atom(i);

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
            if (n == (l+1))
                rcut_psp[l] = VSC1*peak_Atom(i);
            else
                rcut_psp[l] = VSC*peak_Atom(i);

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

} /* Suggested_Params_Hamann */


/********************************
 *				*
 *         solve_Hamann	        *
 *				*
 ********************************/

/*  This routine solves for the Hamann psp

   Entry - num_psp
	   n_psp[],
	   l_psp[],
           e_psp[], an array of suggested psp eigenvalues.
           fill_psp[], an array of suggested psp fillings
  	   rcut_psp[], an array of rcut values

   Uses - Atom data structure

*/
#define	ALAM	3.5

void	solve_Hamann(num_psp,n_psp,l_psp,e_psp,fill_psp,rcut_psp,
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
    int 		iteration,
    converged,
    p,
    i,k,
    match,mch,
    nrc,
    Ngrid;
    double	al,amesh,rmax,texp,Zion;
    double	gamma,gpr,del,nu0,nu1,nu2,r0;
    double	sv,sf,sx;
    double	dcl;
    double	cl[6];
    double	rpsi_match,rpsi_prime_match;
    double	ldpsi_match;
    double	e1l,e2l,eeig;
    double	ph,px,pc,eh,ex,ec;
    double	*w1l,*w1l_prime,
    *w2l,*w2l_prime;
    double	*r,
    *Vh,*Vx,*Vc,
    *Vall,
    *Vcut,
    *Fcut,
    *f,
    *V1l,
    *V2l;

    /* Allocate Grids */
    Vall       = Vall_Atom();
    Vcut       = alloc_LogGrid();
    Fcut       = alloc_LogGrid();
    f	      = alloc_LogGrid();

    r     = r_LogGrid();
    Ngrid = N_LogGrid();
    al    = log_amesh_LogGrid();
    amesh = amesh_LogGrid();

    eeig  = 0.0;
    Zion  = 0.0;
    for (k=0; k<Ngrid; ++k)
        rho_psp[k] = 0.0;

    if (debug_print()){
        printf("\n\nHamann pseudopotential check\n\n");
        printf("l\trcore     rmatch    E in       E psp      norm test slope test\n");
    }
    for (p=0; p<num_psp; ++p)
    {
        w1l       = r_psi_psp[p];
        w1l_prime	= r_psi_prime_psp[p];
        V1l       = V_psp[p];
        w2l       = w1l;
        w2l_prime = w1l_prime;
        V2l	= V1l;


        /* Solve for scattering state if necessary */
        if (fill_psp[p] == 0.0)
        {
            rmax = 2.5*rcut_psp[p];
            solve_Scattering_State_Atom(n_psp[p],l_psp[p],e_psp[p],rmax);

            /* scattering state saved at the end of the atom list */
            i = Nvalence_Atom() + Ncore_Atom();
        }
        /* find state of all-electron wavefunction */
        else
            i     = state_Atom(n_psp[p],l_psp[p]);

        /* find matching point stuff */
        nrc	       = rint(log(rcut_psp[p]/r[0])/al);
        match            = turning_point_Atom(i);
        rpsi_match       = r_psi_Atom(i)[match];
        rpsi_prime_match = r_psi_prime_Atom(i)[match];
        ldpsi_match      = rpsi_prime_match/rpsi_match;

        /* Form the cutoff potential */
        for (k=0; k<Ngrid; ++k)
        {
            texp = pow((r[k]/rcut_psp[p]),ALAM);

            if (texp < 700.0)
                Fcut[k] = exp(-texp);
            else
                Fcut[k] = 0.0;

            Vcut[k] = (1.0 - Fcut[k])*Vall[k];
        }

        iteration = 0;
        cl[p]     = Vall[nrc];
        converged = False;
        while ((iteration <= 30) && (!converged))
        {
            ++iteration;

            /* guess V1l */
            for (k=0; k<Ngrid; ++k)
                V1l[k] = Vcut[k] + cl[p]*Fcut[k];

            /* find w1l, and e1l */
            e1l = e_psp[p];

            /* psp bound state */
            if (fill_psp[p] > 0.0)
            {
                R_Schrodinger(l_psp[p]+1,l_psp[p],V1l,&mch,&e1l,w1l,w1l_prime);
                mch = Ngrid-1;
            }
            /* psp scattering state */
            else
            {
                mch = match;
                R_Schrodinger_Fixed_Logderiv(l_psp[p]+1,l_psp[p],V1l,mch,
                                             ldpsi_match,&e1l,w1l,w1l_prime);
            }


            /* perform integrals */
            nu0 = ((double) (2*l_psp[p]+3));
            nu1 = ((double) (l_psp[p]+1));
            nu2 = ((double) (l_psp[p]+2));
            r0 = r[0]/sqrt(amesh);
            sv = pow(r0,nu0) * (w1l[0]*Fcut[0]/pow(r[0],nu1)) / nu0;
            sf = pow(r0,nu0) * Fcut[0]/nu0;
            sx = sv;
            for (k=0; k<(mch-2); ++k)
            {
                sv += al*(r[k])          * (Fcut[k]*w1l[k]*w1l[k]);
                sf += al*(pow(r[k],nu0)) * (Fcut[k]*Fcut[k]);
                sx += al*(pow(r[k],nu2)) * (Fcut[k]*w1l[k]);
            }
            sv += al*( 23.0*(r[mch-2])   * (Fcut[mch-2]*w1l[mch-2]*w1l[mch-2])
                       + 28.0*(r[mch-1])   * (Fcut[mch-1]*w1l[mch-1]*w1l[mch-1])
                       +  9.0*(r[mch])     * (Fcut[mch]  *w1l[mch]  *w1l[mch])
                     )/24.0;

            sf += al*( 23.0*pow(r[mch-2],nu0)*(Fcut[mch-2]*Fcut[mch-2])
                       + 28.0*pow(r[mch-1],nu0)*(Fcut[mch-1]*Fcut[mch-1])
                       +  9.0*pow(r[mch],nu0)  *(Fcut[mch]  *Fcut[mch])
                     )/24.0;

            sx += al*( 23.0*pow(r[mch-2],nu2)*(Fcut[mch-2]*w1l[mch-2])
                       + 28.0*pow(r[mch-1],nu2)*(Fcut[mch-1]*w1l[mch-1])
                       +  9.0*pow(r[mch],nu2)  *(Fcut[mch]  *w1l[mch]  )
                     )/24.0;

            dcl       = (e_psp[p]-e1l)/sv;
            cl[p]     = cl[p] + dcl;
            converged = (fabs(dcl) <= SMALL);

        } /* while iteration */


        gamma = fabs(rpsi_match/w1l[match]);
        sx    = sx*gamma*gamma;
        sf    = sf*gamma*gamma;
        del   = (-sx
                 + (sx/fabs(sx))*sqrt(sx*sx - sf*(gamma*gamma-1.0))
                )/sf;

        /* construct final psp wavefunction w2l */
        nu1 = ((double) (l_psp[p]+1));
        for (k=0; k<Ngrid; ++k)
            w2l[k] = gamma*(w1l[k] + del*pow(r[k],nu1)*Fcut[k]);

        /* construct final psp V2l */
        for (k=0; k<Ngrid; ++k)
            if ((fabs(w2l[k]) > SMALL) || (r[k]<0.1))              /* hacking */
                V2l[k] = V1l[k]
                         + (gamma*del*pow(r[k],nu1)*Fcut[k]/(2.0*w2l[k]))
                         * (
                             (  ALAM*ALAM*pow(r[k]/rcut_psp[p],(2.0*ALAM))
                                - (2.0*ALAM*((double) l_psp[p]) + ALAM*(ALAM+1.0))
                                *pow(r[k]/rcut_psp[p],ALAM)
                             )/(r[k]*r[k])
                             + 2.0*e_psp[p] - 2.0*V1l[k]
                         );
        e2l = e_psp[p];

        /******************/
        /* verify psp V2l */
        /******************/
        /* psp bound state */
        if (fill_psp[p] > 0.0)
        {
            R_Schrodinger(l_psp[p]+1,l_psp[p],V2l,&mch,&e2l,w2l,w2l_prime);
        }
        /* scattering state */
        else
        {
            R_Schrodinger_Fixed_Logderiv(l_psp[p]+1,l_psp[p],V2l,mch,
                                         ldpsi_match,&e2l,w2l,w2l_prime);
        }

        eeig += fill_psp[p]*e2l;

        /* accumulate charges */
        Zion += fill_psp[p];
        for (k=0; k<Ngrid; ++k)
            rho_psp[k] += fill_psp[p]*pow(w2l[k]/r[k],2.0);

        gamma=fabs(rpsi_match/w2l[match]);
        gpr  =fabs(rpsi_prime_match/w2l_prime[match]);
        if (debug_print()){
            printf("%d\t%lf  %lf  %lf  %lf  %lf  %lf\n",l_psp[p],
                   rcut_psp[p],r[match],
                   e_psp[p],e2l,
                   gamma,gpr);
        }

        /* Extend scattering states to Ngrid */
        if (fill_psp[p] == 0.0)
        {
            R_Schrodinger_Fixed_E(l_psp[p]+1,l_psp[p],V2l,
                                  Ngrid-1,e2l,w2l,w2l_prime);
        }

    } /* for p */

   /* solve for other states */
   if (debug_print())
   {
      printf("\n\nHamann pseudopotential extra states\n\n");
      printf("l\tn\tE psp\n");
   }
   i = 0;
   for (p=0; p<num_psp; ++p)
   {
      V1l = V_psp[p];
      e1l = e_psp[p]+0.1;
      for (k=0; k<(kb_expansion[l_psp[p]]-1); ++k)
      {
         w1l       = r_psi_extra[i];
         w1l_prime = r_psi_prime_extra[i];
         R_Schrodinger(l_psp[p]+2+k,l_psp[p],V1l,&mch,&e1l,w1l,w1l_prime);
         ++i;
         if (debug_print()) printf("%d\t%d\t%lf\n",l_psp[p],l_psp[p]+2+k,e1l);
      }
   }



    /***************************************************/
    /* get the hartree potential an energy             */
    /* get the exchange potential and energy           */
    /* get the correlation potential and energy        */
    /* Semicore corrections added if rho_semicore != 0 */
    /***************************************************/
    Vh = Vcut;
    Vx = Fcut;
    Vc =  w1l_prime;
    ph = R_Hartree_DFT(rho_psp,Zion,Vh);
    eh = 0.5*ph;

    for (k=0;k<Ngrid; ++k)
        f[k] = rho_psp[k] + rho_semicore[k];
    R_Exchange_DFT(f,Vx,&ex,&px);
    R_Correlation_DFT(f,Vc,&ec,&pc);

    R_Screening_Cut(Vx);
    R_Screening_Cut(Vc);


    /* recalculate px and pc */
    for (k=0; k<Ngrid; ++k) f[k] = (rho_psp[k])*Vx[k];
    px = Integrate_LogGrid(f);
    for (k=0; k<Ngrid; ++k) f[k] = (rho_psp[k])*Vc[k];
    pc = Integrate_LogGrid(f);


    *eall_psp = eeig + eh + ex + ec - ph - px - pc;
    *eh_psp   = eh;
    *ph_psp   = ph;
    *ex_psp   = ex;
    *px_psp   = px;
    *ec_psp   = ec;
    *pc_psp   = pc;
    for (p=0; p<num_psp; ++p)
        for (k=0; k<Ngrid; ++k)
            V_psp[p][k] = V_psp[p][k] - Vh[k] - Vx[k] - Vc[k];

    /* deallocate memory */
    dealloc_LogGrid(Vcut);
    dealloc_LogGrid(Fcut);
    dealloc_LogGrid(f);

} /* solve_Hamann */
