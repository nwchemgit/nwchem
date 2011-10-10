/* vanderbilt.c -
   Author - Eric Bylaska
*/

#include        <stdio.h>
#include	"loggrid.h"
#include	"schrodin.h"
#include	"dft.h"
#include	"atom.h"
#include        "name.h"
#include	"vanderbilt.h"
#include        "xpansion2.h"
#include        "gaussj.h"

#define	Max(x,y)	((x>y) ? x : y)
#define	True	1
#define	False	0
#define	SMALL	1.0e-9
#define MAXIT   2

/********************************
 *				*
 *  Suggested_Param_Vanderbilt	*
 *				*
 ********************************/

/*  This routine returns suggested parameters
values for the Vanderbilt prescription.

   Entry -
   Exit	 - num_states_psp
	   n_psp[],
	   l_psp[],
           e_psp[], an array of suggested psp eigenvalues.
           fill_psp[], an array of suggested psp fillings
  	   rcut_psp[], an array of rcut values

   Uses - Atom data structure

*/
#define	VSC	1.5
#define	VSC1	1.5
#define	CSC	3.0
#define	CSC1	4.0

void	Suggested_Param_Vanderbilt(num_states_psp,n_psp,l_psp,
                                e_psp,fill_psp,rcut_psp,
                                rlocal,clocal)

int	*num_states_psp;
int	n_psp[];
int	l_psp[];
double	e_psp[];
double	fill_psp[];
double	rcut_psp[];
double  *rlocal,*clocal;
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

    *rlocal         = rcmax;
    *clocal         = emax;
    *num_states_psp = npsps;

} /* Suggested_Params_Vanderbilt */





/********************************
 *				*
 *         solve_Vanderbilt      *
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

void	solve_Vanderbilt(num_psp,n_psp,l_psp,e_psp,fill_psp,rcut_psp,
                      rlocal,clocal,
                      ns,indx_il,indx_ijl,
                      r_hard_psi_psp,
                      r_psi_psp, r_psi_prime_psp,
                      rho_psp,rho_semicore,
                      V_psp,
                      Vlocal,D0,q,
                      eall_psp,
                      eh_psp,ph_psp,
                      ex_psp,px_psp,
                      ec_psp,pc_psp)

int	num_psp;
int	n_psp[];
int	l_psp[];
double	e_psp[];
double	fill_psp[];
double	rcut_psp[];
double  rlocal,clocal;
int    ns[10],indx_il[4][10],indx_ijl[4][4][10];
double  **r_hard_psi_psp;
double	**r_psi_psp;
double	**r_psi_prime_psp;
double	*rho_psp;
double	*rho_semicore;
double	**V_psp;
double  *Vlocal;
double  *D0,*q;
double	*eall_psp;
double	*eh_psp;
double	*ph_psp;
double	*ex_psp;
double	*px_psp;
double	*ec_psp;
double	*pc_psp;
{
    int 		i,l,j,k,it,
    match,match_R,
    Ngrid;
    double	al,amesh,Zion;
    double	f;
    double	el,eeig,scale;
    double       sum1,sum2;
    double	ph,px,pc,eh,ex,ec;

    double	*ul,*ul_prime;
    double	*wl,*wl_prime;
    double	*r,
    *Vh,*Vx,*Vc,*rho_valence,
    *Vall,
    *Vl;
    double	*chi;


    int    lmax,Nc,Nv;
    double B[250],Binv[250],bb[250],aa[250];
    double D_shift;

    /* Allocate Grids */
    Vall       = Vall_Atom();

    r     = r_LogGrid();
    Ngrid = N_LogGrid();
    al    = log_amesh_LogGrid();
    amesh = amesh_LogGrid();

    Nc	= Ncore_Atom();
    Nv   = Nvalence_Atom();
    Zion        = 0.0;
    for (i=Nc; i<(Nc+Nv); ++i)
        Zion += fill_Atom(i);
    for (k=0; k<Ngrid; ++k)
        rho_psp[k] = 0.0;

    /* define the local pseudopotential Vlocal */
    for (k=0; k<Ngrid; ++k)
    {
        f         = (r[k]/rlocal);
        f         = pow(f,3.5);
        f         = exp(-f);
        Vlocal[k] = (1.0-f)*Vall[k] + clocal*f;

    }

    printf("\n\nVanderbilt pseudopotential check\n\n");
    printf("l\trcore     rmatch    E in       E psp      norm test slope test\n");


    /* compute ul, ul_prime, wl, wl_prime, and chi */
    lmax = 0;
    Vl = alloc_LogGrid();
    for (l=0; l<(num_psp); ++l)
    {
        match_R   = rint(log((10.0)/r[0])/al);
        match     = rint(log(rcut_psp[l]/r[0])/al);
        el        = e_psp[l];
        wl        = r_psi_psp[l];
        wl_prime	= r_psi_prime_psp[l];
        ul        = r_hard_psi_psp[l];
        ul_prime  = alloc_LogGrid();
        chi       = V_psp[l];
        if (l_psp[l] > lmax) lmax = l_psp[l];

        /* Solve for r_hard_psi */
        R_Schrodinger_Fixed_E(n_psp[l],l_psp[l],Vall,
                              match_R,el,ul,ul_prime);


        /* generate pseudo-wavefunction wl (and wl_prime) */
        match            = rint(log(rcut_psp[l]/r[0])/al);

        /* scale ul and ul_prime by fabs(1.0/ul[match]) */
        scale = fabs(1.0/ul[match]);
        for (it=1; it<=MAXIT; ++it)
        {
            for (k=0; k<Ngrid; ++k)
            {
                ul[k]       *= scale;
                ul_prime[k] *= scale;
            }

            /* ul is less than zero */
            if (ul[match] < 0.0)
            {
                for (k=0; k<Ngrid; ++k)
                {
                    ul[k]       = -ul[k];
                    ul_prime[k] = -ul_prime[k];
                }
            }

            /* Troullier-Martin expansion */
            if (init_xpansion2(l_psp[l],match,el,Vall,ul,ul_prime,chi,wl,wl_prime))
            {
                psi_xpansion2();
                dpsi_xpansion2();
                chi_xpansion2();
            }
            /* Troullier-Martin expansion with the curvature constraint removed */
            else
            {
                psi_xpansion2();
                dpsi_xpansion2();
                chi_xpansion2();
            }



            if (it==MAXIT)
                printf("%s\t%lf  %lf  %lf  %lf  %lf  %lf\n",spd_Name(l_psp[l]),
                       rcut_psp[l],r[match],
                       e_psp[l],el,
                       fabs(ul[match]/wl[match]),
                       fabs(ul_prime[match]/wl_prime[match]));


            /* generate chi */
            for (k=0; k<Ngrid; ++k)
                chi[k] = chi[k] - Vlocal[k]*wl[k];

            /* generate scaling factor */
            for (k=0; k<Ngrid; ++k)
                Vl[k] = wl[k]*chi[k]/(r[k]*r[k]);
            scale = fabs(Integrate_LogGrid(Vl));
            scale = 1.0/sqrt(scale);
        }

        dealloc_LogGrid(ul_prime);
    }

    /* compute B(i,j) = <wl(i)|chi(j)> and B^(-1) */
    for (l=0; l<=lmax; ++l)
    {
        /* calculate B */
        for (i=0; i<ns[l]; ++i)
        {
            wl = r_psi_psp[ indx_il[i][l] ];

            for (j=0; j<ns[l]; ++j)
            {
                chi = V_psp[ indx_il[j][l] ];
                for (k=0; k<Ngrid; ++k)
                    Vl[k] = wl[k]*chi[k]/(r[k]*r[k]);

                B[ indx_ijl[i][j][l] ] = Integrate_LogGrid(Vl);
            }
        }

        /* calculate inverse of B */
        for (i=0; i<(ns[l]*ns[l]); ++i) bb[i] =0.0;
        for (i=0; i<ns[l]; ++i)
        {
            bb[ i+i*ns[l] ] = 1.0;
            for (j=0; j<ns[l]; ++j)
                aa[i + j*ns[l] ] = B[ indx_ijl[i][j][l] ];
        }
        gaussj(ns[l],aa,ns[l],bb);
        for (i=0; i<ns[l]; ++i)
            for (j=0; j<ns[l]; ++j)
            {
                Binv[ indx_ijl[i][j][l] ] = bb[ i+j*ns[l] ];
            }

        /* calculated |beta> = B^(-1)*|chi> */
        for (k=0; k<Ngrid; ++k)
        {
            for (i=0; i<ns[l]; ++i)
            {
                bb[i] = 0.0;
                for (j=0; j<ns[l]; ++j)
                {
                    chi   = V_psp[ indx_il[j][l] ];
                    bb[i] = bb[i] + Binv[ indx_ijl[i][j][l] ]*chi[k];
                }/*for j*/
            }/*for i*/

            for (i=0; i<ns[l]; ++i)
            {
                chi    = V_psp[ indx_il[i][l] ];
                chi[k] = bb[i];
            }/*for i*/

        } /*for k*/

    }/*for l*/



    /* calculate Q(i,j) = <ul(i)|ul(j)> - <wl(i)|wl(j)> */
    /* calculate D(i,j) = B(i,j) + e(j)*Q(i,j)          */
    /* used to descreen D */
    rho_valence = rho_valence_Atom();
    for (k=0; k<Ngrid; ++k)
        Vl[k] = rho_valence[k]*Vlocal[k];
    D_shift = Integrate_LogGrid(Vl);

    for (l=0; l<=lmax; ++l)
    {
        /* calculate Q and D */
        for (i=0; i<ns[l]; ++i)
        {
            ul = r_hard_psi_psp[ indx_il[i][l] ];
            wl =      r_psi_psp[ indx_il[i][l] ];

            for (j=0; j<ns[l]; ++j)
            {
                ul_prime = r_hard_psi_psp[ indx_il[j][l] ];
                wl_prime =      r_psi_psp[ indx_il[j][l] ];
                for (k=0; k<Ngrid; ++k)
                    Vl[k] = ul[k]*ul_prime[k]/(r[k]*r[k]);
                sum1 = Integrate_LogGrid(Vl);
                for (k=0; k<Ngrid; ++k)
                    Vl[k] = wl[k]*wl_prime[k]/(r[k]*r[k]);
                sum2 = Integrate_LogGrid(Vl);

                q[ indx_ijl[i][j][l] ]  = sum1-sum2;
                D0[ indx_ijl[i][j][l] ] = B[ indx_ijl[i][j][l] ]
                                          + e_psp[ indx_il[j][l]]
                                          *q[ indx_ijl[i][j][l] ]
                                          - D_shift;
            }
        }
    }
    dealloc_LogGrid(Vl);

    printf("\n\nVanderbilt matrix elements\n");

    for (l=0; l<=lmax; ++l)
    {
        printf("\n");
        printf("B(i,j)   l=%s : ",spd_Name(l));
        for (i=0; i<ns[l]; ++i)
        {
            for (j=0; j<ns[l]; ++j)
                printf("%le ",B[ indx_ijl[i][j][l]]);
            printf("\n");
            printf("               ");
        }
    }

    for (l=0; l<=lmax; ++l)
    {
        printf("\n");
        printf("D0(i,j)  l=%s : ",spd_Name(l));
        for (i=0; i<ns[l]; ++i)
        {
            for (j=0; j<ns[l]; ++j)
                printf("%le ",D0[ indx_ijl[i][j][l]]);
            printf("\n");
            printf("               ");
        }
    }

    for (l=0; l<=lmax; ++l)
    {
        printf("\n");
        printf("q(i,j)  l=%s : ",spd_Name(l));
        for (i=0; i<ns[l]; ++i)
        {
            for (j=0; j<ns[l]; ++j)
                printf("%le ",q[ indx_ijl[i][j][l]]);
            printf("\n");
            printf("              ");
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

    ph = R_Hartree_DFT(rho_valence,Zion,Vh);
    eh = 0.5*ph;

    R_Exchange_DFT(rho_valence,Vx,&ex,&px);
    R_Correlation_DFT(rho_valence,Vc,&ec,&pc);

    /* recalculate px and pc */
    for (k=0; k<Ngrid; ++k) rho_psp[k] = (rho_valence[k])*Vx[k];
    px = Integrate_LogGrid(rho_psp);
    for (k=0; k<Ngrid; ++k) rho_psp[k] = (rho_valence[k])*Vc[k];
    pc = Integrate_LogGrid(rho_psp);

    /* addup the reference density's AE eigenvalues */
    eeig = 0.0;
    for (i=Nc; i<(Nc+Nv); ++i)
        eeig += fill_Atom(i)*eigenvalue_Atom(i);


    *eall_psp = eeig + eh + ex + ec - ph - px - pc;
    *eh_psp   = eh;
    *ph_psp   = ph;
    *ex_psp   = ex;
    *px_psp   = px;
    *ec_psp   = ec;
    *pc_psp   = pc;

    /* descreen Vlocal */
    for (k=0; k<Ngrid; ++k)
        Vlocal[k] = Vlocal[k] - Vh[k] - Vx[k] - Vc[k];

    /* make rho_psp = rho_valence */
    for (k=0; k<Ngrid; ++k)
        rho_psp[k] = rho_valence[k];


    /* deallocate memory */
    dealloc_LogGrid(Vh);
    dealloc_LogGrid(Vx);
    dealloc_LogGrid(Vc);



} /* solve_Vanderbilt */




/* $Id$ */
