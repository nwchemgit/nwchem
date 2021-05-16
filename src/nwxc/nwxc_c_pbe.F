#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_pbe.F
C> The PBE correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the PBE correlation functional
C>
C> The PBE correlation functional [1,2] can be seen as a simplification
C> of the PW91 functional [3,4]. The functional form is a 
C> straightforward GGA. This routine also provide the correlation
C> part of PBESOL [5,6].
C>
C> [1] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>77</b>, 3865-3865 (1996), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865">
C> 10.1103/PhysRevLett.77.3865</a>.
C>
C> [2] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Erratum: Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>78</b>, 1396-1396 (1997), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.78.1396">
C> 10.1103/PhysRevLett.78.1396</a>.
C>
C> [3] J.P. Perdew, Y. Wang, "Accurate and simple density functional
C> for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>33</b>, 8800-8802 (1986), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.33.8800">
C> 10.1103/PhysRevB.33.8800</a>.
C>
C> [4] J.P. Perdew, Y. Wang, "Erratum: Accurate and simple density
C> functional for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>40</b>, 3399-3399 (1989), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.40.3399">
C> 10.1103/PhysRevB.40.3399</a>.
C>
C> [5] J.P. Perdew, A. Ruzsinszky, G.I. Csonka, O.A. Vydrov,
C> G.E. Scuseria, L.A. Constantin, X. Zhou, K. Burke,
C> "Restoring the density-gradient expansion for exchange in solids
C> and surfaces", Phys. Rev. Lett. <b>100</b>, 136406-136409 (2008),
C> DOI: <a href="https://doi.org/10.1103/PhysRevLett.100.136406">
C> 10.1103/PhysRevLett.100.136406</a>.
C>
C> [6] J.P. Perdew, A. Ruzsinszky, G.I. Csonka, O.A. Vydrov,
C> G.E. Scuseria, L.A. Constantin, X. Zhou, K. Burke,
C> "Erratum: Restoring the density-gradient expansion for exchange in solids
C> and surfaces", Phys. Rev. Lett. <b>102</b>, 039902-039902 (2009),
C> DOI: <a href="https://doi.org/10.1103/PhysRevLett.102.039902">
C> 10.1103/PhysRevLett.102.039902</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_pbe_p(param,tol_rho, ipol, nq, wght, rho,
     &                        rgamma, ffunc)
#else
      Subroutine nwxc_c_pbe(param,tol_rho, ipol, nq, wght, rho, rgamma, 
     &                      ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_pbe_d2(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc)
#else
      Subroutine nwxc_c_pbe_d3(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc)
#endif
c
c$Id$
c
#include "nwad.fh"
c
      Implicit none
#include "nwxc_param.fh"
c
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)!< [Input] Parameters of functional
      type(nwad_dble)::beta
#else
      double precision param(*)!< [Input] Parameters of functional
      double precision beta
#endif
#else
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \beta \f$
      double precision beta
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::ffunc(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
c                                   !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
c     double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
c                                   !< and possibly rho
#endif
c
c     Constants in PBE functional
c
      double precision GAMMA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
c     parameter (BETA = 0.06672455060314922d0)
c     parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
      type(nwad_dble)::rho_t(1,2)
      type(nwad_dble)::rgamma_t(1,3)
      type(nwad_dble)::tau_t(1,2)
c
c     Intermediate derivative results, etc.
c
      integer n,i
      type(nwad_dble)::rhoval, gammaval
      type(nwad_dble)::nepsc(1),epsc,H0,zeta,arglog,fAt,fAtnum
      type(nwad_dble)::A,ks,argexp,expinA,Mz,Nzrt
      double precision dnepscdn(2)
      double precision depscdna, depscdnb
      double precision dH0dna, dH0dnb, dH0dg
      double precision dphidna, dphidnb, dphidzeta
      double precision dzetadna, dzetadnb
      double precision darglogdna, darglogdnb, darglogdg
      double precision dfAtdt, dfAtdA
      double precision dfAtnumdt, dfAtnumdA
      double precision dfAtdendt, dfAtdendA
      double precision dfAtdna, dfAtdnb, dfAtdg
      double precision dAdna, dAdnb
      double precision dtdna, dtdnb, dtdg
      double precision dksdna, dksdnb
      double precision dargexpdna, dargexpdnb
      double precision dMzdna, dMzdnb 
      double precision dNzrtdna, dNzrtdnb
      double precision dNzrtdg
      type(nwad_dble)::rrho
      type(nwad_dble)::opz, omz
      type(nwad_dble)::rks
      type(nwad_dble)::rphi, rphi3
      type(nwad_dble)::expmone, expmone2
c     type(nwad_dble)::t, t2, t3, t4
      type(nwad_dble)::t2, t4
      type(nwad_dble)::fAtden,fAtden2
      type(nwad_dble)::phi, phi2, phi3
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision d2nepscdn2(NCOL_AMAT2)
c     double precision d2epscdna2, d2epscdnadnb, d2epscdnb2
c     double precision d2H0dna2, d2H0dnadnb, d2H0dnb2
c     double precision d2H0dnadg, d2H0dnbdg, d2H0dg2
c     double precision d2phidzeta2, d2phidna2, d2phidnadnb, d2phidnb2
c     double precision d2zetadna2, d2zetadnadnb, d2zetadnb2
c     double precision d2arglogdna2, d2arglogdnb2, d2arglogdnadnb
c     double precision d2arglogdnadg, d2arglogdnbdg, d2arglogdg2
c     double precision d2fAtdt2, d2fAtdA2, d2fAtdtdA, d2fAtdg2
c     double precision d2fAtnumdt2, d2fAtnumdtdA, d2fAtnumdA2
c     double precision d2fAtdendt2, d2fAtdendtdA, d2fAtdendA2
c     double precision d2fAtdna2, d2fAtdnb2, d2fAtdnadnb
c     double precision d2fAtdnadg, d2fAtdnbdg
c     double precision d2Adna2, d2Adnadnb, d2Adnb2
c     double precision d2tdna2, d2tdnb2, d2tdnadnb
c     double precision d2tdg2, d2tdnadg, d2tdnbdg
c     double precision d2ksdna2, d2ksdnb2, d2ksdnadnb
c     double precision d2argexpdna2, d2argexpdnb2, d2argexpdnadnb
c     double precision d2Mzdna2, d2Mzdnadnb, d2Mzdnb2
c     double precision d2Nzrtdna2, d2Nzrtdnadnb, d2Nzrtdnb2
c     double precision d2Nzrtdnadg, d2Nzrtdnbdg, d2Nzrtdg2
c     double precision rrho3
c     double precision opz2, omz2
c     double precision rks2
c     double precision rphi5
c     double precision expmone3
c     double precision fAtden3
#endif
#if defined(THIRD_DERIV)
c     double precision d3nepscdn3(NCOL_AMAT3)
c     double precision d3epscdna3, d3epscdna2dnb, d3epscdnadnb2,
c    1                 d3epscdnb3
c     double precision d3ksdna3, d3ksdna2dnb, d3ksdnadnb2, d3ksdnb3
c     double precision d3zetadna3, d3zetadna2dnb, d3zetadnadnb2,
c    1                 d3zetadnb3
c     double precision d3phidzeta3, d3phidna3, d3phidna2dnb,
c    1                 d3phidnadnb2, d3phidnb3
c     double precision d3tdna3, d3tdna2dnb, d3tdnadnb2, d3tdnb3
c     double precision d3argexpdna3, d3argexpdna2dnb, d3argexpdnadnb2,
c    1                 d3argexpdnb3
c     double precision d3Adna3, d3Adna2dnb, d3Adnadnb2, d3Adnb3
c     double precision d3fAtnumdt3, d3fAtdendt3, d3fAtnumdt2dA,
c    1                 d3fAtdendt2dA, d3fAtnumdtdA2, d3fAtdendtdA2,
c    2                 d3fAtnumdA3, d3fAtdendA3   
c     double precision d3fAtdt3, d3fAtdt2dA, d3fAtdtdA2, d3fAtdA3   
c     double precision d3fAtdna3, d3fAtdna2dnb, d3fAtdnadnb2,
c    1                 d3fAtdnb3    
c     double precision d3arglogdna3, d3arglogdna2dnb, d3arglogdnadnb2,
c    1                 d3arglogdnb3
c     double precision d3H0dna3, d3H0dna2dnb, d3H0dnadnb2, d3H0dnb3
c     double precision d3tdg3, d3tdna2dg, d3tdnadnbdg, d3tdnbdnadg, 
c    1                 d3tdnb2dg, d3tdnadg2, d3tdnbdg2
c     double precision d3fAtdg3, d3fAtdna2dg, d3fAtdnadnbdg, 
c    1                 d3fAtdnbdnadg,
c    2                 d3fAtdnb2dg, d3fAtdnadg2, d3fAtdnbdg2
c     double precision d3arglogdg3, d3arglogdna2dg, d3arglogdnadnbdg, 
c    1                 d3arglogdnbdnadg,
c    2                 d3arglogdnb2dg, d3arglogdnadg2, d3arglogdnbdg2
c     double precision d3H0dg3, d3H0dna2dg, d3H0dnadnbdg, d3H0dnbdnadg,
c    1                 d3H0dnb2dg, d3H0dnadg2, d3H0dnbdg2
c     double precision d3Mzdna3, d3Mzdna2dnb, d3Mzdnadnb2, d3Mzdnb3
c     double precision d3Nzrtdna3, d3Nzrtdna2dnb, d3Nzrtdnadnb2, 
c    1                 d3Nzrtdnb3
c     double precision d3Nzrtdg3, d3Nzrtdnadg2, d3Nzrtdnbdg2,
c    1                 d3Nzrtdna2dg, d3Nzrtdnadnbdg, d3Nzrtdnb2dg
c     double precision rrho4
c     double precision opz3, omz3
c     double precision rks3
c     double precision rphi6
c     double precision expmone4
c     double precision fAtden4
#endif
c     
c References:
c [a] J. P. Perdew, K. Burke, and M. Ernzerhof,
c     {\it Generalized gradient approximation made simple},
c     Phys.\ Rev.\ Lett. {\bf 77,} 3865 (1996).
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\it Real-space cutoff
c     construction of a generalized gradient approximation: The PW91
c     density functional}, submitted to Phys.\ Rev.\ B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys.\ Rev.\ B {\bf 45}, 13244 (1992).
c     
c  E_c(PBE) = Int n (epsilon_c + H0) dxdydz
c
c  n*epsilon_c                <=== supplied by another subroutine
c  d(n*epsilon_c)/d(na)       <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(na)d(na) <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(na)d(nb) <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(nb)d(nb) <=== supplied by another subroutine
c
c  H0 = GAMMA * phi**3 * log{ 1 + BETA/GAMMA * t**2 * [ ... ]}
c
c  phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)]
c
c  zeta = (na - nb)/n
c
c  [ ... ] = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4)
c
c  A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
c  t = |Nabla n|/(2*phi*ks*n)
c
c  ks = 2 * (3 * PI**2 * n)**(1/6) / sqrt(PI)
c
c  |Nabla n| = sqrt(g_aa + g_bb + 2*g_ab)
c
c  Names of variables
c
c  E_c(PBE)                  : Ec
c  n (alpha+beta density)    : rhoval
c  na, nb                    : rho(*,2), rho(*,3)
c  epsilon_c                 : epsc
c  H0                        : H0
c  n*epsilon_c               : nepsc
c  phi                       : phi
c  zeta                      : zeta
c  { ... }                   : arglog
c  [ ... ]                   : fAt
c  (1 + A * t**2)            : fAtnum
c  (1 + A * t**2 + A**2 * t**4) : fAtden
c  A                         : A
c  t                         : t
c  |Nabla n|                 : gammaval
c  ks                        : ks
c  {-epsilon_c ... }         : argexp
c  g_aa, g_bb, g_ab          : g
c
c  Derivatives of these are named like d...dna, d2...dnadnb,
c  d2...dna2, etc.
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      beta = param(1)
      pi = acos(-1.0d0)
c
      do 20 n = 1, nq
c
c        n and zeta = (na - nb)/n
c
         if (ipol.eq.2) then
           rhoval   = 0.0d0
           rho_t(1,1) = 0.0d0
           rho_t(1,2) = 0.0d0
c          rgamma_t(1) = 0.0d0 ! dummy, not used
c          rgamma_t(2) = 0.0d0 ! dummy, not used
c          rgamma_t(3) = 0.0d0 ! dummy, not used
c          tau_t(1) = 0.0d0 ! dummy, not used
c          tau_t(2) = 0.0d0 ! dummy, not used
           gammaval = 0.0d0
           if (rho(n,R_A).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_A)
              rho_t(1,1) = rho(n,R_A)
              gammaval = gammaval + rgamma(n,G_AA)
           endif
           if (rho(n,R_B).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_B)
              rho_t(1,2) = rho(n,R_B)
              gammaval = gammaval + rgamma(n,G_BB)
              if (rho(n,R_A).gt.0.5d0*tol_rho) then
                gammaval = gammaval + 2.0d0*rgamma(n,G_AB)
              endif
           endif
         else
           rhoval = rho(n,R_T)
           rho_t(1,1) = rhoval
c          rgamma_t(1) = 0.0d0 ! dummy, not used
c          tau_t(1) = 0.0d0 ! dummy, not used
           gammaval = rgamma(n,G_TT)
         endif
         if (rhoval.le.tol_rho) goto 20
c Daniel (7-24-12): gammaval is gamma^2 in the correlation part of the
c calculation.
c        dsqgamma = sqrt(gammaval)
c Daniel (7-24-12): variable for storing the product of the density
c and the correlation energy per particle, and its derivatives.  This 
c is later used for storing the correlation energy from the LDA part
c of the calculation.
         nepsc(1) = 0.0d0
c        dnepscdn(D1_RA) = 0.0d0
c        if (ipol.eq.2) dnepscdn(D1_RB) = 0.0d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        d2nepscdn2(D2_RA_RA)=0.0d0
c        d2nepscdn2(D2_RA_RB)=0.0d0
c        if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
#if defined(THIRD_DERIV)
c        d3nepscdn3(D3_RA_RA_RA)=0.0d0
c        d3nepscdn3(D3_RA_RA_RB)=0.0d0
c        d3nepscdn3(D3_RA_RB_RB)=0.0d0
c        if (ipol.eq.2) d3nepscdn3(D3_RB_RB_RB)=0.0d0
#endif
c        ========
c        LDA part
c        ========
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_pw91lda",func_lda,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,nepsc)
         call nwxc_c_pw91lda_p(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                         nepsc) 
         call leave_subr("nwxc_c_pw91lda",func_lda,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,nepsc)
#else
         call nwxc_c_pw91lda(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                       nepsc) 
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_pw91lda_d2(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                          nepsc)
#else
         call nwxc_c_pw91lda_d3(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                          nepsc)
#endif
c
c        ==================
c        PBE non-local part
c        ==================
         if(abs(nepsc(1)).lt.tol_rho*tol_rho) goto 20
c
c        epsilon_c = n*epsilon_c / n
c
c Daniel (7-24-12): Regardless of the form of the HEG correlation used
c above, nepsc = eps*qwght*rhoval*cfac.  We need to gather the 
c contributions from different spins.  The derivatives given are
c just the derivative of the epsilon from the local (LDA) part.  Since
c epsilon and its derivatives are multiplied by the density in the
c routines above, we end up taking derivatives of n*epsilon_c / n.
c
         rrho = 1.0d0/rhoval
c        rrho2 = rrho*rrho
c
         epsc = nepsc(1)*rrho
c        if (ipol.eq.1) then
c           depscdna = dnepscdn(D1_RA)*rrho-nepsc*rrho2
c           depscdnb = depscdna
c        else
c           depscdna = dnepscdn(D1_RA)*rrho-nepsc*rrho2
c           depscdnb = dnepscdn(D1_RB)*rrho-nepsc*rrho2
c        endif
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        rrho3 = rrho2*rrho
c        if (ipol.eq.1) then
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)*rrho
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     +2.0d0*nepsc*rrho3
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)*rrho
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     +2.0d0*nepsc*rrho3
c           d2epscdnb2   = d2epscdna2
c        else
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)*rrho
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     +2.0d0*nepsc*rrho3
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)*rrho
c    &                     -dnepscdn(D1_RA)*rrho2
c    &                     -dnepscdn(D1_RB)*rrho2
c    &                     +2.0d0*nepsc*rrho3
c           d2epscdnb2   = d2nepscdn2(D2_RB_RB)*rrho
c    &                     -dnepscdn(D1_RB)*rrho2
c    &                     -dnepscdn(D1_RB)*rrho2
c    &                     +2.0d0*nepsc*rrho3
c        endif
#endif
#if defined(THIRD_DERIV)
c Subtle point here!  We don't actually calculate the first derivative
c of the LDA correlation energy with respect to the beta spin density
c for a restricted calculation, so we need to use that na = nb
c        rrho4 = rrho3*rrho
c        if (ipol.eq.1) then
c           d3epscdna3    = d3nepscdn3(D3_RA_RA_RA)*rrho
c    1                    - d2nepscdn2(D2_RA_RA)*rrho2
c    2                    - 2.0d0*d2nepscdn2(D2_RA_RA)*rrho2
c    3                    + 6.0d0*dnepscdn(D1_RA)*rrho3
c    4                    - 6.0d0*nepsc*rrho4
c           d3epscdna2dnb = d3nepscdn3(D3_RA_RA_RB)*rrho
c    1                    - d2nepscdn2(D2_RA_RA)*rrho2
c    2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
c    3                    + 4.0d0*dnepscdn(D1_RA)*rrho3
c    4                    + 2.0d0*dnepscdn(D1_RA)*rrho3
c    5                    - 6.0d0*nepsc*rrho4
c           d3epscdnadnb2 = d3epscdna2dnb
c           d3epscdnb3    = d3epscdna3
c        else
c           d3epscdna3    = d3nepscdn3(D3_RA_RA_RA)*rrho
c    1                    - 3.0d0*d2nepscdn2(D2_RA_RA)*rrho2
c    2                    + 6.0d0*dnepscdn(D1_RA)*rrho3
c    3                    - 6.0d0*nepsc*rrho4
c           d3epscdna2dnb = d3nepscdn3(D3_RA_RA_RB)*rrho
c    1                    - d2nepscdn2(D2_RA_RA)*rrho2
c    2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
c    3                    + 4.0d0*dnepscdn(D1_RA)*rrho3
c    4                    + 2.0d0*dnepscdn(D1_RB)*rrho3
c    5                    - 6.0d0*nepsc*rrho4
c           d3epscdnadnb2 = d3nepscdn3(D3_RA_RB_RB)*rrho
c    1                    - d2nepscdn2(D2_RB_RB)*rrho2
c    2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
c    3                    + 2.0d0*dnepscdn(D1_RA)*rrho3
c    4                    + 4.0d0*dnepscdn(D1_RB)*rrho3
c    5                    - 6.0d0*nepsc*rrho4
c           d3epscdnb3    = d3nepscdn3(D3_RB_RB_RB)*rrho
c    1                    - 3.0d0*d2nepscdn2(D2_RB_RB)*rrho2
c    2                    + 6.0d0*dnepscdn(D1_RB)*rrho3
c    3                    - 6.0d0*nepsc*rrho4
c        endif
#endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
c Daniel (7-24-12): Thomas-Fermi screening vector
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
c        dksdna = (1.0d0/6.0d0)*ks*rrho
c        dksdnb = dksdna
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        d2ksdna2   = (1.0d0/6.0d0)*dksdna*rrho
c    &              - (1.0d0/6.0d0)*ks*rrho2
c        d2ksdnadnb = d2ksdna2
c        d2ksdnb2   = d2ksdna2
#endif
#if defined(THIRD_DERIV)
c        d3ksdna3    = (1.0d0/6.0d0)*d2ksdna2*rrho
c    1               - (1.0d0/3.0d0)*dksdna*rrho2
c    2               + (1.0d0/3.0d0)*ks*rrho3
c        d3ksdna2dnb = d3ksdna3
c        d3ksdnadnb2 = d3ksdna3
c        d3ksdnb3    = d3ksdna3
#endif
c
c        zeta = (na-nb)/n and its derivs
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,R_A)-rho(n,R_B))/rhoval
         endif
         if(zeta.le.-1.0d0) zeta=-1.0d0
         if(zeta.ge. 1.0d0) zeta= 1.0d0
c        if (ipol.eq.1) then
c           dzetadna = 1.0d0*rrho
c           dzetadnb = -1.0d0*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2zetadna2   = -2.0d0*rrho2
c           d2zetadnadnb = 0.0d0
c           d2zetadnb2   = 2.0d0*rrho2
#endif
#if defined(THIRD_DERIV)
c           d3zetadna3    = 6.0d0*rrho3
c           d3zetadna2dnb = 2.0d0*rrho3
c           d3zetadnadnb2 = -2.0d0*rrho3
c           d3zetadnb3    = -6.0d0*rrho3
#endif
c        else
c           dzetadna =  2.0d0*rho(n,R_B)*rrho2
c           dzetadnb = -2.0d0*rho(n,R_A)*rrho2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2zetadna2   = -4.0d0*rho(n,R_B)*rrho3
c           d2zetadnadnb = 2.0d0*(rho(n,R_A)-rho(n,R_B))*rrho3
c           d2zetadnb2   = 4.0d0*rho(n,R_A)*rrho3
#endif
#if defined(THIRD_DERIV)
c           d3zetadna3    = 12.0d0*rho(n,R_B)*rrho4 
c           d3zetadna2dnb = -4.0d0*(rho(n,R_A)-2.0d0*rho(n,R_B))*
c    1                      rrho4
c           d3zetadnadnb2 = 4.0d0*(rho(n,R_A)-2.0d0*rho(n,R_B))*
c    1                      rrho4
c           d3zetadnb3    = -12.0d0*rho(n,R_A)*rrho4
#endif
c        endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
         opz  = (1.0d0+zeta)
         omz  = (1.0d0-zeta)
c
         phi = 0.0d0
         if (opz.gt.1.0d-10) then
           phi = phi + 0.5d0 * opz**(2.0d0/3.0d0)
         endif
         if (omz.gt.1.0d-10) then
           phi = phi + 0.5d0 * omz**(2.0d0/3.0d0)
         endif
c        phi = 0.5d0*(opz**(2.0d0/3.0d0)
c    &               +omz**(2.0d0/3.0d0))
c        if (omz.lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &             opz**(2.0d0/3.0d0)/opz)
c        else if (opz.lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &            -omz**(2.0d0/3.0d0)/omz)
c        else
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &         opz**(2.0d0/3.0d0)/opz
c    &        -omz**(2.0d0/3.0d0)/omz)
c        endif
c        dphidna = dphidzeta*dzetadna
c        dphidnb = dphidzeta*dzetadnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        opz2 = opz*opz
c        omz2 = omz*omz
c
c        if (omz.lt.tol_rho) then
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         opz**(2.0d0/3.0d0)/opz2)
c        else if (opz.lt.tol_rho) then
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         omz**(2.0d0/3.0d0)/omz2)
c        else
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         opz**(2.0d0/3.0d0)/opz2
c    &        +omz**(2.0d0/3.0d0)/omz2)
c        endif
c
c        d2phidna2   = d2phidzeta2*dzetadna*dzetadna
c    &               + dphidzeta*d2zetadna2 
c        d2phidnadnb = d2phidzeta2*dzetadna*dzetadnb
c    &               + dphidzeta*d2zetadnadnb
c        d2phidnb2   = d2phidzeta2*dzetadnb*dzetadnb
c    &               + dphidzeta*d2zetadnb2 
#endif
#if defined(THIRD_DERIV)
c Daniel (7-30-12): Testing is done here to prevent numerical problems
c        opz3 = opz2*opz
c        omz3 = omz2*omz
c
c        if (omz.lt.tol_rho) then
c           d3phidzeta3 = (4.0d0/27.0d0)*(
c    1         opz**(2.0d0/3.0d0)/opz3)
c        else if (opz.lt.tol_rho) then
c           d3phidzeta3 = (4.0d0/27.0d0)*(
c    1         -omz**(2.0d0/3.0d0)/omz3)
c        else
c           d3phidzeta3 = (4.0d0/27.0d0)*(
c    1         opz**(2.0d0/3.0d0)/opz3
c    2         -omz**(2.0d0/3.0d0)/omz3)
c        endif
c        d3phidna3    = d3phidzeta3*dzetadna*dzetadna*dzetadna
c    1                + 3.0d0*d2phidzeta2*d2zetadna2*dzetadna
c    2                + dphidzeta*d3zetadna3
c        d3phidna2dnb = d3phidzeta3*dzetadna*dzetadna*dzetadnb
c    1                + 2.0d0*d2phidzeta2*d2zetadnadnb*dzetadna
c    2                + d2phidzeta2*d2zetadna2*dzetadnb
c    3                + dphidzeta*d3zetadna2dnb
c        d3phidnadnb2 = d3phidzeta3*dzetadna*dzetadnb*dzetadnb
c    1                + 2.0d0*d2phidzeta2*d2zetadnadnb*dzetadnb
c    2                + d2phidzeta2*d2zetadnb2*dzetadna
c    3                + dphidzeta*d3zetadnadnb2
c        d3phidnb3    = d3phidzeta3*dzetadnb*dzetadnb*dzetadnb
c    1                + 3.0d0*d2phidzeta2*d2zetadnb2*dzetadnb
c    2                + dphidzeta*d3zetadnb3
#endif
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
         rks = 1.0d0/ks
         rphi = 1.0d0/phi
c        rphi2 = rphi*rphi
         rphi3 = rphi**3
c
c        t = dsqgamma*rphi*rks*rrho/2.0d0
c        dtdna = -t*rrho-t*rphi*dphidna-t*rks*dksdna
c        dtdnb = -t*rrho-t*rphi*dphidnb-t*rks*dksdnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        rks2 = rks*rks
c
c        d2tdna2 = - dtdna*rrho
c    &           + t*rrho2
c    &           - dtdna*rphi*dphidna
c    &           + t*rphi2*dphidna*dphidna
c    &           - t*rphi*d2phidna2         
c    &           - dtdna*rks*dksdna
c    &           + t*rks2*dksdna*dksdna   
c    &           - t*rks*d2ksdna2           
c
c        d2tdnadnb = - dtdnb*rrho
c    &           + t*rrho2
c    &           - dtdnb*rphi*dphidna
c    &           + t*rphi2*dphidna*dphidnb
c    &           - t*rphi*d2phidnadnb
c    &           - dtdnb*rks*dksdna
c    &           + t*rks2*dksdna*dksdnb
c    &           - t*rks*d2ksdnadnb
c
c        d2tdnb2 = - dtdnb*rrho
c    &           + t*rrho2
c    &           - dtdnb*rphi*dphidnb
c    &           + t*rphi2*dphidnb*dphidnb
c    &           - t*rphi*d2phidnb2
c    &           - dtdnb*rks*dksdnb
c    &           + t*rks2*dksdnb*dksdnb
c    &           - t*rks*d2ksdnb2
#endif
#if defined(THIRD_DERIV)
c        rks3 = rks2*rks
c        
c        d3tdna3    = -2.0d0*t*rphi3*dphidna*dphidna*dphidna
c    1              - t*rks*d3ksdna3
c    2              - 2.0d0*t*rrho3
c    3              + 3.0d0*t*rphi2*d2phidna2*dphidna
c    4              - t*rphi*d3phidna3
c    5              - 2.0d0*t*rks3*dksdna*dksdna*dksdna
c    6              + 3.0d0*t*rks2*d2ksdna2*dksdna
c    7              + (dphidna*dphidna*rphi2
c    8                - d2phidna2*rphi
c    9                + dksdna*dksdna*rks2
c    A                - d2ksdna2*rks
c    B                + rrho2)*dtdna*2.0d0
c    C              - (dphidna*rphi + dksdna*rks
c    D                + rrho)*d2tdna2
c        d3tdna2dnb = -2.0d0*t*rphi3*dphidna*dphidna*dphidnb
c    1              - t*rks*d3ksdna2dnb
c    2              - 2.0d0*t*rrho3
c    3              + 2.0d0*t*rphi2*d2phidnadnb*dphidna
c    4              + t*rphi2*d2phidna2*dphidnb
c    5              - t*rphi*d3phidna2dnb
c    6              - 2.0d0*t*rks3*dksdna*dksdna*dksdnb
c    7              + 2.0d0*t*rks2*d2ksdnadnb*dksdna
c    8              + t*rks2*d2ksdna2*dksdnb
c    9              + (dphidna*dphidna*rphi2
c    A                - d2phidna2*rphi
c    B                + dksdna*dksdna*rks2
c    C                - d2ksdna2*rks
c    D                + rrho2)*dtdnb
c    E              + (dphidna*dphidnb*rphi2
c    F                - d2phidnadnb*rphi
c    G                + dksdna*dksdnb*rks2
c    H                - d2ksdnadnb*rks
c    I                + rrho2)*dtdna
c    J              - (dphidna*rphi + dksdna*rks
c    K                + rrho)*d2tdnadnb
c        d3tdnadnb2 = -2.0d0*t*rphi3*dphidna*dphidnb*dphidnb
c    1              - t*rks*d3ksdnadnb2
c    2              - 2.0d0*t*rrho3
c    3              + 2.0d0*t*rphi2*d2phidnadnb*dphidnb
c    4              + t*rphi2*d2phidnb2*dphidna
c    5              - t*rphi*d3phidnadnb2
c    6              - 2.0d0*t*rks3*dksdna*dksdnb*dksdnb
c    7              + 2.0d0*t*rks2*d2ksdnadnb*dksdnb
c    8              + t*rks2*d2ksdnb2*dksdna
c    9              + (dphidna*dphidnb*rphi2
c    A                - d2phidnadnb*rphi
c    B                + dksdna*dksdnb*rks2
c    C                - d2ksdnadnb*rks
c    D                + rrho2)*dtdnb
c    E              + (dphidna*dphidnb*rphi2
c    F                - d2phidnadnb*rphi
c    G                + dksdna*dksdnb*rks2
c    H                - d2ksdnadnb*rks
c    I                + rrho2)*dtdnb
c    J              - (dphidna*rphi + dksdna*rks
c    K                + rrho)*d2tdnb2
c        d3tdnb3    = -2.0d0*t*rphi3*dphidnb*dphidnb*dphidnb
c    1              - t*rks*d3ksdnb3
c    2              - 2.0d0*t*rrho3
c    3              + 3.0d0*t*rphi2*d2phidnb2*dphidnb
c    4              - t*rphi*d3phidnb3
c    5              - 2.0d0*t*rks3*dksdnb*dksdnb*dksdnb
c    6              + 3.0d0*t*rks2*d2ksdnb2*dksdnb
c    7              + (dphidnb*dphidnb*rphi2
c    8                - d2phidnb2*rphi
c    9                + dksdnb*dksdnb*rks2
c    A                - d2ksdnb2*rks
c    B                + rrho2)*dtdnb*2.0d0
c    C              - (dphidnb*rphi + dksdnb*rks
c    D                + rrho)*d2tdnb2
#endif
c
c        { ... } in A (see below) and its derivs
c
c        rphi4 = rphi3*rphi
c
         argexp = -epsc*rphi3/GAMMA
c        dargexpdna = -depscdna*rphi3/GAMMA
c    &                +3.0d0*epsc*rphi4*dphidna/GAMMA
c        dargexpdnb = -depscdnb*rphi3/GAMMA
c    &                +3.0d0*epsc*rphi4*dphidnb/GAMMA
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        rphi5 = rphi4*rphi
c        d2argexpdna2 = -d2epscdna2*rphi3/GAMMA
c    &        +3.0d0*depscdna*rphi4*dphidna/GAMMA
c    &        +3.0d0*depscdna*rphi4*dphidna/GAMMA
c    &        -12.0d0*epsc*rphi5*dphidna*dphidna/GAMMA
c    &        +3.0d0*epsc*rphi4*d2phidna2/GAMMA
c        d2argexpdnadnb = -d2epscdnadnb*rphi3/GAMMA
c    &        +3.0d0*depscdna*rphi4*dphidnb/GAMMA
c    &        +3.0d0*depscdnb*rphi4*dphidna/GAMMA
c    &        -12.0d0*epsc*rphi5*dphidna*dphidnb/GAMMA
c    &        +3.0d0*epsc*rphi4*d2phidnadnb/GAMMA
c        d2argexpdnb2 = -d2epscdnb2*rphi3/GAMMA
c    &        +3.0d0*depscdnb*rphi4*dphidnb/GAMMA
c    &        +3.0d0*depscdnb*rphi4*dphidnb/GAMMA
c    &        -12.0d0*epsc*rphi5*dphidnb*dphidnb/GAMMA
c    &        +3.0d0*epsc*rphi4*d2phidnb2/GAMMA
#endif
c
#ifdef THIRD_DERIV
c        rphi6 = rphi5*rphi
c
c        d3argexpdna3    = -12.0d0*rphi5*
c    1                     ( 3.0d0*dphidna*dphidna*depscdna
c    2                     + 3.0d0*dphidna*d2phidna2*epsc)/GAMMA
c    3                   + 60.0d0*rphi6*
c    4                     dphidna*dphidna*dphidna*epsc/GAMMA
c    5                   + 3.0d0*rphi4*
c    6                     ( 3.0d0*d2phidna2*depscdna
c    7                     + 3.0d0*dphidna*d2epscdna2
c    8                     + d3phidna3*epsc)/GAMMA
c    9                   - d3epscdna3*rphi3/GAMMA
c        d3argexpdna2dnb = -12.0d0*rphi5*
c    1                     ( 2.0d0*dphidna*dphidnb*depscdna
c    2                     + dphidna*dphidna*depscdnb
c    3                     + 2.0d0*dphidna*d2phidnadnb*epsc
c    4                     + dphidnb*d2phidna2*epsc)/GAMMA
c    5                   + 60.0d0*rphi6*
c    6                     dphidna*dphidna*dphidnb*epsc/GAMMA
c    7                   + 3.0d0*rphi4*
c    8                     ( 2.0d0*d2phidnadnb*depscdna
c    9                     + d2phidna2*depscdnb
c    A                     + 2.0d0*dphidna*d2epscdnadnb
c    B                     + dphidnb*d2epscdna2
c    C                     + d3phidna2dnb*epsc)/GAMMA
c    D                   - d3epscdna2dnb*rphi3/GAMMA
c        d3argexpdnadnb2 = -12.0d0*rphi5*
c    1                     ( 2.0d0*dphidna*dphidnb*depscdnb
c    2                     + dphidnb*dphidnb*depscdna
c    3                     + 2.0d0*dphidnb*d2phidnadnb*epsc
c    4                     + dphidna*d2phidnb2*epsc)/GAMMA
c    5                   + 60.0d0*rphi6*
c    6                     dphidna*dphidnb*dphidnb*epsc/GAMMA
c    7                   + 3.0d0*rphi4*
c    8                     ( 2.0d0*d2phidnadnb*depscdnb
c    9                     + d2phidnb2*depscdna
c    A                     + 2.0d0*dphidnb*d2epscdnadnb
c    B                     + dphidna*d2epscdnb2
c    C                     + d3phidnadnb2*epsc)/GAMMA
c    D                   - d3epscdnadnb2*rphi3/GAMMA
c        d3argexpdnb3    = -12.0d0*rphi5*
c    1                     ( 3.0d0*dphidnb*dphidnb*depscdnb
c    2                     + 3.0d0*dphidnb*d2phidnb2*epsc)/GAMMA
c    3                   + 60.0d0*rphi6*
c    4                     dphidnb*dphidnb*dphidnb*epsc/GAMMA
c    5                   + 3.0d0*rphi4*
c    6                     ( 3.0d0*d2phidnb2*depscdnb
c    7                     + 3.0d0*dphidnb*d2epscdnb2
c    8                     + d3phidnb3*epsc)/GAMMA
c    9                   - d3epscdnb3*rphi3/GAMMA
#endif
c
c        A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
         if (abs(argexp).lt.EXPTOL) then
            expinA=exp(argexp)
         else
            expinA=0.0d0
         endif
         expmone = expinA - 1.0d0
         expmone2 = expmone*expmone
c
         A = BETA/GAMMA/expmone
c        dAdna = -BETA/GAMMA*dargexpdna*expinA/expmone2
c        dAdnb = -BETA/GAMMA*dargexpdnb*expinA/expmone2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        expmone3 = expmone2*expmone
c
c        d2Adna2   = -BETA/GAMMA*d2argexpdna2
c    &               *expinA/expmone2
c    &             - BETA/GAMMA*dargexpdna
c    &               *dargexpdna*expinA/expmone2
c    &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdna
c    &               *expinA*expinA/expmone3
c        d2Adnadnb  = -BETA/GAMMA*d2argexpdnadnb
c    &               *expinA/expmone2
c    &             - BETA/GAMMA*dargexpdna
c    &               *dargexpdnb*expinA/expmone2
c    &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdnb
c    &               *expinA*expinA/expmone3
c        d2Adnb2   = -BETA/GAMMA*d2argexpdnb2
c    &               *expinA/expmone2
c    &             - BETA/GAMMA*dargexpdnb
c    &               *dargexpdnb*expinA/expmone2
c    &             + 2.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb
c    &               *expinA*expinA/expmone3
#endif
#ifdef THIRD_DERIV
c        expmone4 = expmone3*expmone
c
c        d3Adna3    = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdna*
c    1                dargexpdna*expinA*expinA*expinA/
c    2                expmone4
c    3              + 6.0d0*BETA/GAMMA*
c    4                ( dargexpdna*d2argexpdna2
c    5                + dargexpdna*dargexpdna*dargexpdna)*
c    6                expinA*expinA/expmone3
c    7              - BETA/GAMMA*
c    8                ( d3argexpdna3 + 3.0d0*d2argexpdna2*dargexpdna
c    9                + dargexpdna*dargexpdna*dargexpdna)*
c    A                expinA/expmone2
c        d3Adna2dnb = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdna*
c    1                dargexpdnb*expinA*expinA*expinA/
c    2                expmone4
c    3              + 2.0d0*BETA/GAMMA*
c    4                ( 2.0d0*dargexpdna*d2argexpdnadnb
c    5                + 2.0d0*dargexpdna*dargexpdna*dargexpdnb
c    6                + dargexpdnb*d2argexpdna2
c    7                + dargexpdna*dargexpdna*dargexpdnb)*
c    8                expinA*expinA/expmone3
c    9              - BETA/GAMMA*
c    A                ( d3argexpdna2dnb
c    B                + 2.0d0*d2argexpdnadnb*dargexpdna
c    C                + d2argexpdna2*dargexpdnb
c    D                + dargexpdna*dargexpdna*dargexpdnb)*
c    E                expinA/expmone2
c        d3Adnadnb2 = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdnb*
c    1                dargexpdnb*expinA*expinA*expinA/
c    2                expmone4
c    3              + 2.0d0*BETA/GAMMA*
c    4                ( 2.0d0*dargexpdnb*d2argexpdnadnb
c    5                + 2.0d0*dargexpdna*dargexpdnb*dargexpdnb
c    6                + dargexpdna*d2argexpdnb2
c    7                + dargexpdna*dargexpdnb*dargexpdnb)*
c    8                expinA*expinA/expmone3
c    9              - BETA/GAMMA*
c    A                ( d3argexpdnadnb2
c    B                + 2.0d0*d2argexpdnadnb*dargexpdnb
c    C                + d2argexpdnb2*dargexpdna
c    D                + dargexpdna*dargexpdnb*dargexpdnb)*
c    E                expinA/expmone2
c        d3Adnb3    = -6.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb*
c    1                dargexpdnb*expinA*expinA*expinA/
c    2                expmone4
c    3              + 6.0d0*BETA/GAMMA*
c    4                ( dargexpdnb*d2argexpdnb2
c    5                + dargexpdnb*dargexpdnb*dargexpdnb)*
c    6                expinA*expinA/expmone3
c    7              - BETA/GAMMA*
c    8                ( d3argexpdnb3 + 3.0d0*d2argexpdnb2*dargexpdnb
c    9                + dargexpdnb*dargexpdnb*dargexpdnb)*
c    A                expinA/expmone2
#endif
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
         t2 = gammaval*(rphi*rks*rrho/2.0d0)**2.0d0
c        t2 = t*t
c        t3 = t2*t
         t4 = t2*t2
c
         fAtnum = 1.0d0+A*t2
         fAtden = 1.0d0+A*t2+A*A*t4
         fAtden2 = fAtden*fAtden
         fAt = fAtnum/fAtden
c        dfAtnumdt = 2.0d0*A*t
c        dfAtnumdA = t2
c        dfAtdendt = 2.0d0*A*t+4.0d0*A*A*t3
c        dfAtdendA = t2+2.0d0*A*t4
c        dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/fAtden2
c        dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/fAtden2
c        dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
c        dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        fAtden3 = fAtden2*fAtden
c
c        d2fAtnumdt2 = 2.0d0*A
c        d2fAtdendt2 = 2.0d0*A+12.0d0*A**2*t2
c        d2fAtnumdtdA = 2.0d0*t
c        d2fAtdendtdA = 2.0d0*t+8.0d0*A*t3
c        d2fAtnumdA2 = 0.0d0
c        d2fAtdendA2 = 2.0d0*t4
c        d2fAtdt2  = (d2fAtnumdt2*fAtden-fAtnum*d2fAtdendt2)
c    &               /(fAtden2)
c    &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
c    &               /(fAtden3)*dfAtdendt
c        d2fAtdtdA = (d2fAtnumdtdA*fAtden+dfAtnumdt*dfAtdendA
c    &               -dfAtnumdA*dfAtdendt-fAtnum*d2fAtdendtdA)
c    &               /(fAtden2)
c    &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
c    &               /(fAtden3)*dfAtdendA
c        d2fAtdA2  = (d2fAtnumdA2*fAtden-fAtnum*d2fAtdendA2)
c    &               /(fAtden2)
c    &               -2.0d0*(dfAtnumdA*fAtden-fAtnum*dfAtdendA)
c    &               /(fAtden3)*dfAtdendA
c        d2fAtdna2 = d2fAtdt2*dtdna*dtdna + d2fAtdtdA*dtdna*dAdna
c    &             + dfAtdt*d2tdna2 + d2fAtdtdA*dAdna*dtdna
c    &             + d2fAtdA2*dAdna*dAdna + dfAtdA*d2Adna2
c        d2fAtdnb2 = d2fAtdt2*dtdnb*dtdnb + d2fAtdtdA*dtdnb*dAdnb
c    &             + dfAtdt*d2tdnb2 + d2fAtdtdA*dAdnb*dtdnb
c    &             + d2fAtdA2*dAdnb*dAdnb + dfAtdA*d2Adnb2
c        d2fAtdnadnb = d2fAtdt2*dtdna*dtdnb + d2fAtdtdA*dtdna*dAdnb
c    &             + dfAtdt*d2tdnadnb + d2fAtdtdA*dAdna*dtdnb
c    &             + d2fAtdA2*dAdna*dAdnb + dfAtdA*d2Adnadnb
#endif
#ifdef THIRD_DERIV
c        fAtden4 = fAtden3*fAtden
c
c        d3fAtnumdt3   = 0.0d0
c        d3fAtdendt3   = 24.0d0*A**2*t
c        d3fAtnumdt2dA = 2.0d0
c        d3fAtdendt2dA = 2.0d0+24.0d0*A*t2
c        d3fAtnumdtdA2 = 0.0d0
c        d3fAtdendtdA2 = 8.0d0*t3
c        d3fAtnumdA3   = 0.0d0
c        d3fAtdendA3   = 0.0d0
c Here we pay attention to which derivatives are zero.
c        d3fAtdt3   = -( fAtnum*d3fAtdendt3
c    1                 + 3.0d0*dfAtnumdt*d2fAtdendt2
c    2                 + 3.0d0*d2fAtnumdt2*dfAtdendt)/
c    3                 (fAtden2)
c    4                + 6.0d0*( fAtnum*d2fAtdendt2*dfAtdendt
c    5                        + dfAtnumdt*dfAtdendt*dfAtdendt)/
c    6                 (fAtden3)
c    7                - 6.0d0*fAtnum*dfAtdendt*dfAtdendt*dfAtdendt/
c    8                 (fAtden4)
c        d3fAtdt2dA = -( fAtnum*d3fAtdendt2dA
c    1                 + 2.0d0*dfAtnumdt*d2fAtdendtdA
c    2                 + dfAtnumdA*d2fAtdendt2
c    3                 + 2.0d0*d2fAtnumdtdA*dfAtdendt
c    4                 + d2fAtnumdt2*dfAtdendA)/
c    5                 (fAtden2)
c    6                + 2.0d0*( 2.0d0*fAtnum*d2fAtdendtdA*dfAtdendt
c    7                        + fAtnum*d2fAtdendt2*dfAtdendA
c    8                        + 2.0d0*dfAtnumdt*dfAtdendA*dfAtdendt
c    9                        + dfAtnumdA*dfAtdendt*dfAtdendt)/
c    A                 (fAtden3)
c    B                - 6.0d0*fAtnum*dfAtdendt*dfAtdendt*dfAtdendA/
c    C                 (fAtden4)
c    D                + d3fAtnumdt2dA/fAtden
c        d3fAtdtdA2 = -( fAtnum*d3fAtdendtdA2
c    1                 + 2.0d0*dfAtnumdA*d2fAtdendtdA
c    2                 + dfAtnumdt*d2fAtdendA2
c    3                 + 2.0d0*d2fAtnumdtdA*dfAtdendA)/
c    4                 (fAtden2)
c    5                + 2.0d0*( 2.0d0*fAtnum*d2fAtdendtdA*dfAtdendA
c    6                        + fAtnum*d2fAtdendA2*dfAtdendt
c    7                        + 2.0d0*dfAtnumdA*dfAtdendA*dfAtdendt
c    8                        + dfAtnumdt*dfAtdendA*dfAtdendA)/
c    9                 (fAtden3)
c    A                - 6.0d0*fAtnum*dfAtdendt*dfAtdendA*dfAtdendA/
c    B                 (fAtden4)
c        d3fAtdA3   = -3.0d0*dfAtnumdA*d2fAtdendA2/
c    1                 (fAtden2)
c    2                + 6.0d0*( fAtnum*d2fAtdendA2*dfAtdendA
c    3                        + dfAtnumdA*dfAtdendA*dfAtdendA)/
c    3                 (fAtden3)
c    4                - 6.0d0*fAtnum*dfAtdendA*dfAtdendA*dfAtdendA/
c    5                 (fAtden4)
c
c        d3fAtdna3    = d3fAtdA3*dAdna*dAdna*dAdna
c    1                + 3.0d0*d3fAtdtdA2*dAdna*dAdna*dtdna
c    2                + 3.0d0*d2fAtdA2*d2Adna2*dAdna
c    3                + 3.0d0*d3fAtdt2dA*dAdna*dtdna*dtdna
c    4                + 3.0d0*d2fAtdtdA*(d2Adna2*dtdna + dAdna*d2tdna2)
c    5                + dfAtdA*d3Adna3 + d3fAtdt3*dtdna*dtdna*dtdna
c    6                + dfAtdt*d3tdna3
c    7                + 3.0d0*d2fAtdt2*d2tdna2*dtdna
c        d3fAtdna2dnb = d3fAtdA3*dAdna*dAdna*dAdnb
c    1                + 2.0d0*d3fAtdtdA2*dAdna*dAdnb*dtdna
c    2                + d3fAtdtdA2*dAdna*dAdna*dtdnb
c    3                + 2.0d0*d2fAtdA2*d2Adnadnb*dAdna
c    4                + d2fAtdA2*d2Adna2*dAdnb
c    5                + 2.0d0*d3fAtdt2dA*dAdna*dtdna*dtdnb
c    6                + d3fAtdt2dA*dAdnb*dtdna*dtdna
c    7                + 2.0d0*d2fAtdtdA*(d2Adnadnb*dtdna
c    8                                 + dAdna*d2tdnadnb)
c    9                + d2fAtdtdA*(d2Adna2*dtdnb + dAdnb*d2tdna2)
c    A                + dfAtdA*d3Adna2dnb + d3fAtdt3*dtdna*dtdna*dtdnb
c    B                + dfAtdt*d3tdna2dnb
c    C                + 2.0d0*d2fAtdt2*d2tdnadnb*dtdna
c    D                + d2fAtdt2*d2tdna2*dtdnb
c        d3fAtdnadnb2 = d3fAtdA3*dAdna*dAdnb*dAdnb
c    1                + 2.0d0*d3fAtdtdA2*dAdna*dAdnb*dtdnb
c    2                + d3fAtdtdA2*dAdnb*dAdnb*dtdna
c    3                + 2.0d0*d2fAtdA2*d2Adnadnb*dAdnb
c    4                + d2fAtdA2*d2Adnb2*dAdna
c    5                + 2.0d0*d3fAtdt2dA*dAdnb*dtdna*dtdnb
c    6                + d3fAtdt2dA*dAdna*dtdnb*dtdnb
c    7                + 2.0d0*d2fAtdtdA*(d2Adnadnb*dtdnb
c    8                                 + dAdnb*d2tdnadnb)
c    9                + d2fAtdtdA*(d2Adnb2*dtdna + dAdna*d2tdnb2)
c    A                + dfAtdA*d3Adnadnb2 + d3fAtdt3*dtdna*dtdnb*dtdnb
c    B                + dfAtdt*d3tdnadnb2
c    C                + 2.0d0*d2fAtdt2*d2tdnadnb*dtdnb
c    D                + d2fAtdt2*d2tdnb2*dtdna
c        d3fAtdnb3    = d3fAtdA3*dAdnb*dAdnb*dAdnb
c    1                + 3.0d0*d3fAtdtdA2*dAdnb*dAdnb*dtdnb
c    2                + 3.0d0*d2fAtdA2*d2Adnb2*dAdnb
c    3                + 3.0d0*d3fAtdt2dA*dAdnb*dtdnb*dtdnb
c    4                + 3.0d0*d2fAtdtdA*(d2Adnb2*dtdnb + dAdnb*d2tdnb2)
c    5                + dfAtdA*d3Adnb3 + d3fAtdt3*dtdnb*dtdnb*dtdnb
c    6                + dfAtdt*d3tdnb3
c    7                + 3.0d0*d2fAtdt2*d2tdnb2*dtdnb
#endif
c
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + BETA/GAMMA*t2*fAt
c        darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
c    &                            +t2*dfAtdna)
c        darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
c    &                            +t2*dfAtdnb)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        d2arglogdna2 = BETA/GAMMA*(2.0d0*dtdna*dtdna*fAt
c    &                             +2.0d0*t*d2tdna2*fAt
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +t2*d2fAtdna2)
c        d2arglogdnb2 = BETA/GAMMA*(2.0d0*dtdnb*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnb2*fAt
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +t2*d2fAtdnb2)
c        d2arglogdnadnb = BETA/GAMMA*(2.0d0*dtdna*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnadnb*fAt
c    &                             +2.0d0*t*dtdna*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdna
c    &                             +t2*d2fAtdnadnb)
#endif
#if defined(THIRD_DERIV)
c        d3arglogdna3    = BETA/GAMMA*(6.0d0*dfAtdna*dtdna*dtdna
c    1                                +6.0d0*dfAtdna*d2tdna2*t
c    2                                +6.0d0*d2fAtdna2*dtdna*t
c    3                                +6.0d0*fAt*dtdna*d2tdna2
c    4                                +d3fAtdna3*t2
c    5                                +2.0d0*fAt*d3tdna3*t)
c        d3arglogdna2dnb = BETA/GAMMA*(4.0d0*dfAtdna*dtdna*dtdnb
c    1                                +2.0d0*dfAtdnb*dtdna*dtdna
c    2                                +4.0d0*dfAtdna*d2tdnadnb*t
c    3                                +2.0d0*dfAtdnb*d2tdna2*t
c    3                                +4.0d0*d2fAtdnadnb*dtdna*t
c    4                                +2.0d0*d2fAtdna2*dtdnb*t
c    5                                +4.0d0*fAt*dtdna*d2tdnadnb
c    6                                +2.0d0*fAt*dtdnb*d2tdna2
c    7                                +d3fAtdna2dnb*t2
c    8                                +2.0d0*fAt*d3tdna2dnb*t)
c        d3arglogdnadnb2 = BETA/GAMMA*(4.0d0*dfAtdnb*dtdna*dtdnb
c    1                                +2.0d0*dfAtdna*dtdnb*dtdnb
c    2                                +4.0d0*dfAtdnb*d2tdnadnb*t
c    3                                +2.0d0*dfAtdna*d2tdnb2*t
c    3                                +4.0d0*d2fAtdnadnb*dtdnb*t
c    4                                +2.0d0*d2fAtdnb2*dtdna*t
c    5                                +4.0d0*fAt*dtdnb*d2tdnadnb
c    6                                +2.0d0*fAt*dtdna*d2tdnb2
c    7                                +d3fAtdnadnb2*t2
c    8                                +2.0d0*fAt*d3tdnadnb2*t)
c        d3arglogdnb3    = BETA/GAMMA*(6.0d0*dfAtdnb*dtdnb*dtdnb
c    1                                +6.0d0*dfAtdnb*d2tdnb2*t
c    2                                +6.0d0*d2fAtdnb2*dtdnb*t
c    3                                +6.0d0*fAt*dtdnb*d2tdnb2
c    4                                +d3fAtdnb3*t2
c    5                                +2.0d0*fAt*d3tdnb3*t)
#endif
c
c        H0 = GAMMA * phi**3 * log{arglog} and its derivs
c
c Daniel - I redefine the enhancement factor as a product of two
c functions:
c
c        Mz   = phi**3
c        Nzrt = log{arglog}
c        H0   = GAMMA * Mz * Nzrt
c
c This makes the third derivatives substantially easier to read.
         phi2 = phi*phi
         phi3 = phi2*phi
         Mz = phi3
c        dMzdna = 3.0d0*phi2*dphidna
c        dMzdnb = 3.0d0*phi2*dphidnb
c
         Nzrt = log(arglog)
c        dNzrtdna = darglogdna/arglog
c        dNzrtdnb = darglogdnb/arglog
c
         H0 = GAMMA*Mz*Nzrt
c        dH0dna = GAMMA*(dMzdna*Nzrt + Mz*dNzrtdna)
c        dH0dnb = GAMMA*(dMzdnb*Nzrt + Mz*dNzrtdnb)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        d2Mzdna2 = 3.0d0*phi2*d2phidna2 + 6.0d0*phi*dphidna*dphidna
c        d2Mzdnadnb = 3.0d0*phi2*d2phidnadnb
c    1              + 6.0d0*phi*dphidna*dphidnb
c        d2Mzdnb2 = 3.0d0*phi2*d2phidnb2 + 6.0d0*phi*dphidnb*dphidnb
c
c        d2Nzrtdna2 = d2arglogdna2/arglog -
c    1                darglogdna*darglogdna/arglog/arglog
c        d2Nzrtdnadnb = d2arglogdnadnb/arglog -
c    1                  darglogdna*darglogdnb/arglog/arglog
c        d2Nzrtdnb2 = d2arglogdnb2/arglog -
c    1                darglogdnb*darglogdnb/arglog/arglog
c
c        d2H0dna2 = GAMMA*( d2Mzdna2*Nzrt + 2.0d0*dMzdna*dNzrtdna
c    &                    + Mz*d2Nzrtdna2)
c        d2H0dnadnb = GAMMA*( d2Mzdnadnb*Nzrt + dMzdna*dNzrtdnb
c    &                      + dMzdnb*dNzrtdna + Mz*d2Nzrtdnadnb)
c        d2H0dnb2 = GAMMA*( d2Mzdnb2*Nzrt + 2.0d0*dMzdnb*dNzrtdnb
c    &                    + Mz*d2Nzrtdnb2)         
#endif
#ifdef THIRD_DERIV
c        d3Mzdna3 = 3.0d0*phi2*d3phidna3
c    1            + 18.0d0*phi*d2phidna2*dphidna
c    2            + 6.0d0*dphidna*dphidna*dphidna
c        d3Mzdna2dnb = 3.0d0*phi2*d3phidna2dnb
c    1               + 12.0d0*phi*d2phidnadnb*dphidna
c    2               + 6.0d0*phi*d2phidna2*dphidnb
c    3               + 6.0d0*dphidna*dphidna*dphidnb
c        d3Mzdnadnb2 = 3.0d0*phi2*d3phidnadnb2
c    1               + 12.0d0*phi*d2phidnadnb*dphidnb
c    2               + 6.0d0*phi*d2phidnb2*dphidna
c    3               + 6.0d0*dphidna*dphidnb*dphidnb
c        d3Mzdnb3 = 3.0d0*phi2*d3phidnb3
c    1            + 18.0d0*phi*d2phidnb2*dphidnb
c    2            + 6.0d0*dphidnb*dphidnb*dphidnb
c
c        d3Nzrtdna3 = d3arglogdna3/arglog
c    1              - 3.0d0*d2arglogdna2*darglogdna/arglog/arglog
c    2              + 2.0d0*darglogdna*darglogdna*darglogdna/
c    3                arglog/arglog/arglog
c        d3Nzrtdna2dnb = d3arglogdna2dnb/arglog
c    1                 - 2.0d0*d2arglogdnadnb*darglogdna/arglog/arglog
c    2                 - d2arglogdna2*darglogdnb/arglog/arglog
c    3                 + 2.0d0*darglogdna*darglogdna*darglogdnb/
c    4                   arglog/arglog/arglog
c        d3Nzrtdnadnb2 = d3arglogdnadnb2/arglog
c    1                 - 2.0d0*d2arglogdnadnb*darglogdnb/arglog/arglog
c    2                 - d2arglogdnb2*darglogdna/arglog/arglog
c    3                 + 2.0d0*darglogdna*darglogdnb*darglogdnb/
c    4                   arglog/arglog/arglog
c        d3Nzrtdnb3 = d3arglogdnb3/arglog
c    1              - 3.0d0*d2arglogdnb2*darglogdnb/arglog/arglog
c    2              + 2.0d0*darglogdnb*darglogdnb*darglogdnb/
c    3                arglog/arglog/arglog
c
c        d3H0dna3    = GAMMA*( d3Mzdna3*Nzrt
c    1                       + 3.0d0*d2Mzdna2*dNzrtdna
c    2                       + 3.0d0*dMzdna*d2Nzrtdna2
c    3                       + Mz*d3Nzrtdna3)
c        d3H0dna2dnb = GAMMA*( d3Mzdna2dnb*Nzrt
c    1                       + 2.0d0*d2Mzdnadnb*dNzrtdna
c    2                       + d2Mzdna2*dNzrtdnb
c    3                       + 2.0d0*dMzdna*d2Nzrtdnadnb
c    4                       + dMzdnb*d2Nzrtdna2
c    5                       + Mz*d3Nzrtdna2dnb)
c        d3H0dnadnb2 = GAMMA*( d3Mzdnadnb2*Nzrt
c    1                       + 2.0d0*d2Mzdnadnb*dNzrtdnb
c    2                       + d2Mzdnb2*dNzrtdna
c    3                       + 2.0d0*dMzdnb*d2Nzrtdnadnb
c    4                       + dMzdna*d2Nzrtdnb2
c    5                       + Mz*d3Nzrtdnadnb2)
c        d3H0dnb3    = GAMMA*( d3Mzdnb3*Nzrt
c    1                       + 3.0d0*d2Mzdnb2*dNzrtdnb
c    2                       + 3.0d0*dMzdnb*d2Nzrtdnb2
c    3                       + Mz*d3Nzrtdnb3)
#endif
c
c        Now we update Ec, Amat, and Amat2
c
c        if (lfac) then
c           if (ldew) ffunc(n)=ffunc(n)+nepsc*fac               
c        endif
c        if (lfac) then
c           Amat(n,1) = Amat(n,1) + dnepscdn(1)*fac
c           if (ipol.eq.2) Amat(n,2) = Amat(n,2) + dnepscdn(2)*fac
#ifdef SECOND_DERIV
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) 
c    &                        + d2nepscdn2(D2_RA_RA)*fac
c           Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) 
c    &                        + d2nepscdn2(D2_RA_RB)*fac
c           if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) 
c    &                        + d2nepscdn2(D2_RB_RB)*fac
#endif
c        endif
         ffunc(n)=ffunc(n)+wght*(nepsc(1) + H0*rhoval)
c        Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                 +  (dnepscdn(D1_RA) + H0 + rhoval*dH0dna)*wght
c        if (ipol.eq.2) then
c          Amat(n,D1_RB) = Amat(n,D1_RB)
c    &                   +  (dnepscdn(D1_RB) + H0 + rhoval*dH0dnb)*wght
c        endif
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c        Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &   + (d2nepscdn2(D2_RA_RA) + 2.d0*dH0dna + rhoval*d2H0dna2)*wght
c        Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &   + (d2nepscdn2(D2_RA_RB) + dH0dna + dH0dnb
c    &      + rhoval*d2H0dnadnb)*wght
c        if (ipol.eq.2) then
c          Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &     + (d2nepscdn2(D2_RB_RB) + 2.d0*dH0dnb + rhoval*d2H0dnb2)*wght
c        endif
#endif
#if defined(THIRD_DERIV)
c Daniel (7-31-12): XC third derivatives (LDA part)
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    1                           + d3nepscdn3(D3_RA_RA_RA)*wght        
c           Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)         
c    1                           + d3nepscdn3(D3_RA_RA_RB)*wght
c           Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
c    1                           + d3nepscdn3(D3_RA_RB_RB)*wght
c           if (ipol.eq.2)
c    1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)       
c    2                           + d3nepscdn3(D3_RB_RB_RB)*wght  
c Daniel (7-31-12): XC third derivatives (GGA part)
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    1        + (3.0d0*d2H0dna2 + rhoval*d3H0dna3)*wght
c           Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)
c    1        + (d2H0dna2 + 2.0d0*d2H0dnadnb + rhoval*d3H0dna2dnb)*wght
c           Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
c    1        + (d2H0dnb2 + 2.0d0*d2H0dnadnb + rhoval*d3H0dnadnb2)*wght
c           if (ipol.eq.2)
c    1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
c    2        + (3.0d0*d2H0dnb2 + rhoval*d3H0dnb3)*wght
#endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
c        if (dsqgamma.gt.TOLL)then
c           dtdg = 0.25d0*rphi*rks*rrho/dsqgamma
c           dfAtdg = dfAtdt*dtdg
c           darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t2*dfAtdg)
c           dNzrtdg = darglogdg/arglog
c           dH0dg = GAMMA*Mz*dNzrtdg
c           if (ipol.eq.1) then
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c           else
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*wght
c           endif
#ifdef SECOND_DERIV
c           d2tdg2 = -0.125d0*rphi*rks*rrho/(dsqgamma*dsqgamma*dsqgamma)
c           d2tdnadg = -dtdg*rrho-dtdg*rphi*dphidna
c    &                 -dtdg*rks*dksdna
c           d2tdnbdg = -dtdg*rrho-dtdg*rphi*dphidnb
c    &                 -dtdg*rks*dksdnb
c           d2fAtdg2 = d2fAtdt2*dtdg*dtdg+dfAtdt*d2tdg2
c           d2fAtdnadg = d2fAtdt2*dtdg*dtdna
c    &                  +d2fAtdtdA*dtdg*dAdna
c    &                  +dfAtdt*d2tdnadg 
c           d2fAtdnbdg = d2fAtdt2*dtdg*dtdnb
c    &                  +d2fAtdtdA*dtdg*dAdnb
c    &                  +dfAtdt*d2tdnbdg 
c           d2arglogdnadg = BETA/GAMMA*(2.0d0*dtdna*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnadg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdna
c    &                                 +2.0d0*t*dtdna*dfAtdg
c    &                                 +t2*d2fAtdnadg)
c           d2arglogdnbdg = BETA/GAMMA*(2.0d0*dtdnb*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnbdg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdnb
c    &                                 +2.0d0*t*dtdnb*dfAtdg
c    &                                 +t2*d2fAtdnbdg)
c           d2arglogdg2 = BETA/GAMMA*(2.0d0*dtdg*dtdg*fAt
c    &                               +2.0d0*t*d2tdg2*fAt
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +t2*d2fAtdg2)
c
c           d2Nzrtdg2 = d2arglogdg2/arglog -
c    1                  darglogdg*darglogdg/arglog/arglog
c           d2Nzrtdnadg = d2arglogdnadg/arglog -
c    1                    darglogdna*darglogdg/arglog/arglog
c           d2Nzrtdnbdg = d2arglogdnbdg/arglog -
c    1                    darglogdnb*darglogdg/arglog/arglog
c
c           d2H0dnadg = GAMMA*dMzdna*dNzrtdg
c    &                + GAMMA*Mz*d2Nzrtdnadg
c           d2H0dnbdg = GAMMA*dMzdnb*dNzrtdg
c    &                + GAMMA*Mz*d2Nzrtdnbdg
c           d2H0dg2 = GAMMA*Mz*d2Nzrtdg2
c
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
c    &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &             + d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
c    &             + 2.0d0*d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
c    &             + d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
c    &             + 4.0d0*d2H0dg2*rhoval*wght
c           if (ipol.eq.2) then
c             Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
c    &               + (dH0dg + d2H0dnbdg*rhoval)*wght
c             Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
c    &               + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*wght
c             Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &               + (dH0dg + d2H0dnbdg*rhoval)*wght
c             Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
c    &               + 2.0d0*d2H0dg2*rhoval*wght
c             Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &               + d2H0dg2*rhoval*wght
c           endif
#endif
#ifdef THIRD_DERIV
c Derivatives of t, the reduced density gradient
c           d3tdg3 = 0.1875d0*rphi*rks*rrho/dsqgamma/dsqgamma/
c    1               dsqgamma/dsqgamma/dsqgamma
c           d3tdna2dg = dtdg*( dphidna*dphidna*rphi2
c    1                       - d2phidna2*rphi
c    2                       + dksdna*dksdna*rks2
c    3                       - d2ksdna2*rks
c    4                       + rrho2)
c    5                - d2tdnadg*( dphidna*rphi
c    6                           + dksdna*rks
c    7                           + rrho)
c           d3tdnadnbdg = dtdg*( dphidna*dphidnb*rphi2
c    1                       - d2phidnadnb*rphi
c    2                       + dksdna*dksdnb*rks2
c    3                       - d2ksdnadnb*rks
c    4                       + rrho2)
c    5                  - d2tdnbdg*( dphidna*rphi
c    6                             + dksdna*rks
c    7                             + rrho)
c           d3tdnbdnadg = dtdg*( dphidna*dphidnb*rphi2
c    1                       - d2phidnadnb*rphi
c    2                       + dksdna*dksdnb*rks2
c    3                       - d2ksdnadnb*rks
c    4                       + rrho2)
c    5                  - d2tdnadg*( dphidnb*rphi
c    6                             + dksdnb*rks
c    7                             + rrho)
c           d3tdnb2dg = dtdg*( dphidnb*dphidnb*rphi2
c    1                       - d2phidnb2*rphi
c    2                       + dksdnb*dksdnb*rks2
c    3                       - d2ksdnb2*rks
c    4                       + rrho2)
c    5                - d2tdnbdg*( dphidnb*rphi
c    6                           + dksdnb*rks
c    7                           + rrho)
c           d3tdnadg2 = -d2tdg2*( dphidna*rphi + dksdna*rks
c    1                          + rrho)
c           d3tdnbdg2 = -d2tdg2*( dphidnb*rphi + dksdnb*rks
c    1                          + rrho)
c Derivatives of the rational function fAt
c           d3fAtdg3      = d3fAtdt3*dtdg*dtdg*dtdg
c    1                    + 3.0d0*d2fAtdt2*d2tdg2*dtdg
c    2                    + dfAtdt*d3tdg3
c           d3fAtdna2dg   = d3fAtdtdA2*dAdna*dAdna*dtdg
c    1                    + 2.0d0*d3fAtdt2dA*dAdna*dtdna*dtdg
c    2                    + d2fAtdtdA*( d2Adna2*dtdg
c    3                                + 2.0d0*dAdna*d2tdnadg)
c    4                    + d3fAtdt3*dtdna*dtdna*dtdg
c    5                    + d2fAtdt2*( d2tdna2*dtdg
c    6                               + 2.0d0*d2tdnadg*dtdna)
c    7                    + dfAtdt*d3tdna2dg
c           d3fAtdnadnbdg = d3fAtdtdA2*dAdna*dAdnb*dtdg
c    1                    + d3fAtdt2dA*( dAdna*dtdnb
c    2                                 + dAdnb*dtdna)*dtdg
c    3                    + d2fAtdtdA*( d2Adnadnb*dtdg
c    4                                + dAdna*d2tdnbdg
c    5                                + dAdnb*d2tdnadg)
c    6                    + d3fAtdt3*dtdna*dtdnb*dtdg
c    7                    + d2fAtdt2*( d2tdnadnb*dtdg
c    8                               + d2tdnadg*dtdnb
c    9                               + d2tdnbdg*dtdna)
c    A                    + dfAtdt*d3tdnadnbdg
c           d3fAtdnbdnadg = d3fAtdnadnbdg
c           d3fAtdnb2dg   = d3fAtdtdA2*dAdnb*dAdnb*dtdg
c    1                    + 2.0d0*d3fAtdt2dA*dAdnb*dtdnb*dtdg
c    2                    + d2fAtdtdA*( d2Adnb2*dtdg
c    3                                + 2.0d0*dAdnb*d2tdnbdg)
c    4                    + d3fAtdt3*dtdnb*dtdnb*dtdg
c    5                    + d2fAtdt2*( d2tdnb2*dtdg
c    6                               + 2.0d0*d2tdnbdg*dtdnb)
c    7                    + dfAtdt*d3tdnb2dg
c           d3fAtdnadg2   = d3fAtdt2dA*dAdna*dtdg*dtdg
c    1                    + d2fAtdtdA*dAdna*d2tdg2
c    2                    + d3fAtdt3*dtdna*dtdg*dtdg
c    3                    + d2fAtdt2*( 2.0d0*d2tdnadg*dtdg
c    4                               + dtdna*d2tdg2)
c    5                    + dfAtdt*d3tdnadg2
c           d3fAtdnbdg2   = d3fAtdt2dA*dAdnb*dtdg*dtdg
c    1                    + d2fAtdtdA*dAdnb*d2tdg2
c    2                    + d3fAtdt3*dtdnb*dtdg*dtdg
c    3                    + d2fAtdt2*( 2.0d0*d2tdnbdg*dtdg
c    4                               + dtdnb*d2tdg2)
c    5                    + dfAtdt*d3tdnbdg2
c Derivatives of the logarithm argument arglog
c           d3arglogdg3      = BETA/GAMMA*(d3fAtdg3*t2
c    1                       + 6.0d0*d2fAtdg2*dtdg*t
c    2                       + 6.0d0*dfAtdg*( dtdg*dtdg
c    3                                      + d2tdg2*t)
c    4                       + 6.0d0*fAt*dtdg*d2tdg2
c    5                       + 2.0d0*fAt*d3tdg3*t)
c           d3arglogdna2dg   = 2.0d0*BETA/GAMMA*(
c    1                         2.0d0*d2fAtdnadg*dtdna*t
c    2                       + 2.0d0*dfAtdna*d2tdnadg*t
c    3                       + 2.0d0*dfAtdna*dtdna*dtdg
c    4                       + dfAtdg*dtdna*dtdna
c    5                       + 2.0d0*fAt*d2tdnadg*dtdna
c    6                       + d2fAtdna2*dtdg*t
c    7                       + dfAtdg*d2tdna2*t
c    8                       + fAt*d3tdna2dg*t
c    9                       + fAt*d2tdna2*dtdg
c    A                       + 0.50d0*d3fAtdna2dg*t2)
c           d3arglogdnadnbdg = 2.0d0*BETA/GAMMA*(
c    1                         d2fAtdnadg*dtdnb*t
c    2                       + d2fAtdnbdg*dtdna*t
c    3                       + dfAtdna*d2tdnbdg*t
c    4                       + dfAtdnb*d2tdnadg*t
c    5                       + dfAtdna*dtdnb*dtdg
c    6                       + dfAtdnb*dtdna*dtdg
c    7                       + dfAtdg*dtdna*dtdnb
c    8                       + fAt*d2tdnadg*dtdnb
c    9                       + fAt*d2tdnbdg*dtdna
c    A                       + d2fAtdnadnb*dtdg*t
c    B                       + dfAtdg*d2tdnadnb*t
c    C                       + fAt*d3tdnadnbdg*t
c    D                       + fAt*d2tdnadnb*dtdg
c    E                       + 0.50d0*d3fAtdnadnbdg*t2)
c           d3arglogdnbdnadg = 2.0d0*BETA/GAMMA*(
c    1                         d2fAtdnbdg*dtdna*t
c    2                       + d2fAtdnadg*dtdnb*t
c    3                       + dfAtdnb*d2tdnadg*t
c    4                       + dfAtdna*d2tdnbdg*t
c    5                       + dfAtdnb*dtdna*dtdg
c    6                       + dfAtdna*dtdnb*dtdg
c    7                       + dfAtdg*dtdnb*dtdna
c    8                       + fAt*d2tdnbdg*dtdna
c    9                       + fAt*d2tdnadg*dtdnb
c    A                       + d2fAtdnadnb*dtdg*t
c    B                       + dfAtdg*d2tdnadnb*t
c    C                       + fAt*d3tdnbdnadg*t
c    D                       + fAt*d2tdnadnb*dtdg
c    E                       + 0.50d0*d3fAtdnbdnadg*t2)
c           d3arglogdnb2dg   = 2.0d0*BETA/GAMMA*(
c    1                         2.0d0*d2fAtdnbdg*dtdnb*t
c    2                       + 2.0d0*dfAtdnb*d2tdnbdg*t
c    3                       + 2.0d0*dfAtdnb*dtdnb*dtdg
c    4                       + dfAtdg*dtdnb*dtdnb
c    5                       + 2.0d0*fAt*d2tdnbdg*dtdnb
c    6                       + d2fAtdnb2*dtdg*t
c    7                       + dfAtdg*d2tdnb2*t
c    8                       + fAt*d3tdnb2dg*t
c    9                       + fAt*d2tdnb2*dtdg
c    A                       + 0.50d0*d3fAtdnb2dg*t2)
c           d3arglogdnadg2   = 2.0d0*BETA/GAMMA*(
c    1                         2.0d0*d2fAtdnadg*dtdg*t
c    2                       + dfAtdna*d2tdg2*t
c    3                       + dfAtdna*dtdg*dtdg
c    4                       + d2fAtdg2*dtdna*t
c    5                       + 2.0d0*dfAtdg*d2tdnadg*t
c    6                       + 2.0d0*dfAtdg*dtdna*dtdg
c    7                       + 2.0d0*fAt*d2tdnadg*dtdg
c    8                       + fAt*dtdna*d2tdg2
c    9                       + fAt*d3tdnadg2*t
c    A                       + 0.50d0*d3fAtdnadg2*t2)
c           d3arglogdnbdg2   = 2.0d0*BETA/GAMMA*(
c    1                         2.0d0*d2fAtdnbdg*dtdg*t
c    2                       + dfAtdnb*d2tdg2*t
c    3                       + dfAtdnb*dtdg*dtdg
c    4                       + d2fAtdg2*dtdnb*t
c    5                       + 2.0d0*dfAtdg*d2tdnbdg*t
c    6                       + 2.0d0*dfAtdg*dtdnb*dtdg
c    7                       + 2.0d0*fAt*d2tdnbdg*dtdg
c    8                       + fAt*dtdnb*d2tdg2
c    9                       + fAt*d3tdnbdg2*t
c    A                       + 0.50d0*d3fAtdnbdg2*t2)
c Derivatives of the correlation enhancement factor
c           d3Nzrtdg3 = d3arglogdg3/arglog
c    1                - 3.0d0*d2arglogdg2*darglogdg/arglog/arglog
c    2                + 2.0d0*darglogdg*darglogdg*darglogdg/arglog/
c    3                  arglog/arglog
c           d3Nzrtdnadg2 = d3arglogdnadg2/arglog
c    1                   - (2.0d0*d2arglogdnadg*darglogdg +
c    2                      darglogdna*d2arglogdg2)/arglog/arglog
c    3                   + 2.0d0*darglogdna*darglogdg*darglogdg/
c    4                     arglog/arglog/arglog
c           d3Nzrtdnbdg2 = d3arglogdnbdg2/arglog
c    1                   - (2.0d0*d2arglogdnbdg*darglogdg +
c    2                      darglogdnb*d2arglogdg2)/arglog/arglog
c    3                   + 2.0d0*darglogdnb*darglogdg*darglogdg/
c    4                     arglog/arglog/arglog
c           d3Nzrtdna2dg = d3arglogdna2dg/arglog
c    1                   - (2.0d0*d2arglogdnadg*darglogdna +
c    2                      d2arglogdna2*darglogdg)/arglog/arglog
c    3                   + 2.0d0*darglogdna*darglogdna*darglogdg/
c    4                     arglog/arglog/arglog
c           d3Nzrtdnadnbdg = d3arglogdnadnbdg/arglog
c    1                     - (d2arglogdnadg*darglogdnb +
c    2                        d2arglogdnbdg*darglogdna +
c    3                        d2arglogdnadnb*darglogdg)/arglog/arglog
c    4                     + 2.0d0*darglogdna*darglogdnb*darglogdg/
c    5                       arglog/arglog/arglog
c           d3Nzrtdnb2dg = d3arglogdnb2dg/arglog
c    1                   - (2.0d0*d2arglogdnbdg*darglogdnb +
c    2                      d2arglogdnb2*darglogdg)/arglog/arglog
c    3                   + 2.0d0*darglogdnb*darglogdnb*darglogdg/
c    4                     arglog/arglog/arglog
c
c           d3H0dg3      = GAMMA*Mz*d3Nzrtdg3
c           d3H0dna2dg   = GAMMA*( d2Mzdna2*dNzrtdg
c    1                           + 2.0d0*dMzdna*d2Nzrtdnadg
c    2                           + Mz*d3Nzrtdna2dg)
c           d3H0dnadnbdg = GAMMA*( d2Mzdnadnb*dNzrtdg 
c    1                           + dMzdna*d2Nzrtdnbdg
c    2                           + dMzdnb*d2Nzrtdnadg
c    3                           + Mz*d3Nzrtdnadnbdg)
c           d3H0dnbdnadg = d3H0dnadnbdg
c           d3H0dnb2dg   = GAMMA*( d2Mzdnb2*dNzrtdg
c    1                           + 2.0d0*dMzdnb*d2Nzrtdnbdg
c    2                           + Mz*d3Nzrtdnb2dg)
c           d3H0dnadg2   = GAMMA*( dMzdna*d2Nzrtdg2
c    1                           + Mz*d3Nzrtdnadg2)
c           d3H0dnbdg2   = GAMMA*( dMzdnb*d2Nzrtdg2
c    1                           + Mz*d3Nzrtdnbdg2)
c There are 31 unique 3rd order functional derivatives involving the 
c gradient of the electron density.  Note that many permutations of
c indices are identical so only one permutation is stored. 
c Keep in mind that gamma = gaa + 2*gab + gbb for a spin polarized
c functional, so dgamma/dgaa = 1 and dgamma/dgab = 2
c
c Notes so far (may reduce storage requirements, worry about this when
c we know the code works)
c 1. D3_RA_RA_GAA = D3_RA_RA_GBB
c 2. D3_RA_RB_GAA = D3_RA_RB_GBB
c 3. D3_RB_RB_GAA = D3_RB_RB_GBB 
c 4. D3_RA_GAA_GAA = D3_RA_GAA_GBB = D3_RA_GBB_GBB
c 5. D3_RA_GAA_GAB = D3_RA_GAB_GBB
c 6. D3_RB_GAA_GAA = D3_RB_GAA_GBB = D3_RB_GBB_GBB
c 7. D3_RB_GAA_GAB = D3_RB_GAB_GBB
c 8. D3_GAA_GAA_GAA = D3_GAA_GAA_GBB = D3_GAA_GBB_GBB = D3_GBB_GBB_GBB
c 9. D3_GAA_GAA_GAB = D3_GAA_GAB_GBB = D3_GAB_GBB_GBB
c 10. D3_GAA_GAB_GAB = D3_GAB_GAB_GBB
c 
c It looks like we can remove 15 of the 31 elements, which could be a
c huge benefit here. 
c
c Mixed derivatives dradradg
c           Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
c    1             + (2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*wght
c           Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB)
c    1             + 2.0d0*(2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*wght
c           Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB)
c    1             + (2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*wght
c Mixed derivatives dradrbdg
c           Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA)
c    1             + ( d2H0dnadg + d2H0dnbdg 
c    2               + d3H0dnadnbdg*rhoval)*wght
c           Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB)
c    1             + 2.0d0*( d2H0dnadg + d2H0dnbdg 
c    2                     + d3H0dnadnbdg*rhoval)*wght
c           Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB)
c    1             + ( d2H0dnadg + d2H0dnbdg 
c    2               + d3H0dnadnbdg*rhoval)*wght
c Mixed derivatives drbdrbdg
c           if (ipol.eq.2) then
c             Cmat3(n,D3_RB_RB_GAA) = Cmat3(n,D3_RB_RB_GAA)
c    1               + (2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*wght
c             Cmat3(n,D3_RB_RB_GAB) = Cmat3(n,D3_RB_RB_GAB)
c    1               + 2.0d0*(2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*wght
c             Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
c    1               + (2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*wght
c           endif
c Mixed derivatives dradgdg
c           Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
c    1             + (d2H0dg2 + d3H0dnadg2*rhoval)*wght
c           Cmat3(n,D3_RA_GAA_GAB) = Cmat3(n,D3_RA_GAA_GAB)
c    1             + 2.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*wght
c           Cmat3(n,D3_RA_GAA_GBB) = Cmat3(n,D3_RA_GAA_GBB)
c    1             + (d2H0dg2 + d3H0dnadg2*rhoval)*wght
c           Cmat3(n,D3_RA_GAB_GAB) = Cmat3(n,D3_RA_GAB_GAB)
c    1             + 4.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*wght
c           Cmat3(n,D3_RA_GAB_GBB) = Cmat3(n,D3_RA_GAB_GBB)
c    1             + 2.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*wght
c           Cmat3(n,D3_RA_GBB_GBB) = Cmat3(n,D3_RA_GBB_GBB)
c    1             + (d2H0dg2 + d3H0dnadg2*rhoval)*wght
c Mixed derivatives drbdgdg
c           if (ipol.eq.2) then
c             Cmat3(n,D3_RB_GAA_GAA) = Cmat3(n,D3_RB_GAA_GAA)
c    1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c             Cmat3(n,D3_RB_GAA_GAB) = Cmat3(n,D3_RB_GAA_GAB)
c    1               + 2.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c             Cmat3(n,D3_RB_GAA_GBB) = Cmat3(n,D3_RB_GAA_GBB)
c    1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c             Cmat3(n,D3_RB_GAB_GAB) = Cmat3(n,D3_RB_GAB_GAB)
c    1               + 4.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c             Cmat3(n,D3_RB_GAB_GBB) = Cmat3(n,D3_RB_GAB_GBB)
c    1               + 2.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c             Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB)
c    1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*wght
c           endif
c Derivatives dgaadgdg
c           Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
c    1             + (d3H0dg3*rhoval)*wght
c           Cmat3(n,D3_GAA_GAA_GAB) = Cmat3(n,D3_GAA_GAA_GAB)
c    1             + 2.0d0*(d3H0dg3*rhoval)*wght
c           Cmat3(n,D3_GAA_GAA_GBB) = Cmat3(n,D3_GAA_GAA_GBB)
c    1             + (d3H0dg3*rhoval)*wght
c           Cmat3(n,D3_GAA_GAB_GAB) = Cmat3(n,D3_GAA_GAB_GAB)
c    1             + 4.0d0*(d3H0dg3*rhoval)*wght
c           Cmat3(n,D3_GAA_GAB_GBB) = Cmat3(n,D3_GAA_GAB_GBB)
c    1             + 2.0d0*(d3H0dg3*rhoval)*wght
c           Cmat3(n,D3_GAA_GBB_GBB) = Cmat3(n,D3_GAA_GBB_GBB)
c    1             + (d3H0dg3*rhoval)*wght
c Derivatives dgabdgdg
c           Cmat3(n,D3_GAB_GAB_GAB) = Cmat3(n,D3_GAB_GAB_GAB)
c    1             + 8.0d0*(d3H0dg3*rhoval)*wght
c           if (ipol.eq.2) then
c             Cmat3(n,D3_GAB_GAB_GBB) = Cmat3(n,D3_GAB_GAB_GBB)
c    1               + 4.0d0*(d3H0dg3*rhoval)*wght
c             Cmat3(n,D3_GAB_GBB_GBB) = Cmat3(n,D3_GAB_GBB_GBB)
c    1               + 2.0d0*(d3H0dg3*rhoval)*wght
c Derivatives dgbbdgdg
c             Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
c    1               + (d3H0dg3*rhoval)*wght
c           endif
#endif
c        endif
   20 continue
c
      return
      end
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_pbe.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_pbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_pbe.F"
#endif
#undef NWAD_PRINT
C> @}
