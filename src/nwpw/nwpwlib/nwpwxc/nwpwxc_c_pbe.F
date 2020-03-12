#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_pbe.F
C> The PBE correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_pbe(param,tol_rho,ipol,nq,wght,rho,rgamma, 
     &                      ffunc, Amat, Cmat) 
#else
      Subroutine nwpwxc_c_pbe_d2(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc, Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id$
c
      Implicit none
#include "nwpwxc_param.fh"
c
c
c     Input and other parameters
c
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \beta \f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision ffunc(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
                                    !< and possibly rho
#endif
c
c     Constants in PBE functional
c
      double precision GAMMA, BETA, PI
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
      double precision rho_t(2)
      double precision dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n,i
      double precision rhoval, gammaval
      double precision nepsc, dnepscdn(2)
      double precision epsc, depscdna, depscdnb
      double precision H0, dH0dna, dH0dnb, dH0dg
      double precision phi, dphidna, dphidnb, dphidzeta
      double precision zeta, dzetadna, dzetadnb
      double precision arglog, darglogdna, darglogdnb, darglogdg
      double precision fAt, dfAtdt, dfAtdA
      double precision fAtnum, dfAtnumdt, dfAtnumdA
      double precision fAtden, dfAtdendt, dfAtdendA
      double precision dfAtdna, dfAtdnb, dfAtdg
      double precision A, dAdna, dAdnb
      double precision t, dtdna, dtdnb, dtdg
      double precision ks, dksdna, dksdnb
      double precision argexp, dargexpdna, dargexpdnb
      double precision expinA
#ifdef SECOND_DERIV
      double precision d2nepscdn2(NCOL_AMAT2)
      double precision d2epscdna2, d2epscdnadnb, d2epscdnb2
      double precision d2H0dna2, d2H0dnadnb, d2H0dnb2
      double precision d2H0dnadg, d2H0dnbdg, d2H0dg2
      double precision d2phidzeta2, d2phidna2, d2phidnadnb, d2phidnb2
      double precision d2zetadna2, d2zetadnadnb, d2zetadnb2
      double precision d2arglogdna2, d2arglogdnb2, d2arglogdnadnb
      double precision d2arglogdnadg, d2arglogdnbdg, d2arglogdg2
      double precision d2fAtdt2, d2fAtdA2, d2fAtdtdA, d2fAtdg2
      double precision d2fAtnumdt2, d2fAtnumdtdA, d2fAtnumdA2
      double precision d2fAtdendt2, d2fAtdendtdA, d2fAtdendA2
      double precision d2fAtdna2, d2fAtdnb2, d2fAtdnadnb
      double precision d2fAtdnadg, d2fAtdnbdg
      double precision d2Adna2, d2Adnadnb, d2Adnb2
      double precision d2tdna2, d2tdnb2, d2tdnadnb
      double precision d2tdg2, d2tdnadg, d2tdnbdg
      double precision d2ksdna2, d2ksdnb2, d2ksdnadnb
      double precision d2argexpdna2, d2argexpdnb2, d2argexpdnadnb
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
           rhoval = rho(n,R_A)+rho(n,R_B)
         else
           rhoval = rho(n,R_T)
         endif
         rho_t(1) = rhoval
         if (ipol.eq.2) then
            rho_t(1) = rho(n,R_A)
            rho_t(2) = rho(n,R_B)
         endif
         if (rhoval.le.tol_rho) goto 20
         if (ipol.eq.1) then
            gammaval = rgamma(n,G_TT)
         else
            gammaval = rgamma(n,G_AA)
     &               + rgamma(n,G_BB)
     &               + rgamma(n,G_AB)*2.0d0
         endif
         dsqgamma = max(dsqrt(gammaval),tol_rho)
         nepsc = 0.0d0
         dnepscdn(D1_RA) = 0.0d0
         if (ipol.eq.2) dnepscdn(D1_RB) = 0.0d0
#ifdef SECOND_DERIV
         d2nepscdn2(D2_RA_RA)=0.0d0
         d2nepscdn2(D2_RA_RB)=0.0d0
         if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
c        ========
c        LDA part
c        ========
c
#ifndef SECOND_DERIV
         call nwpwxc_c_pw91lda(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                       nepsc, dnepscdn) 
#else
         call nwpwxc_c_pw91lda_d2(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                          nepsc, dnepscdn, d2nepscdn2)
#endif
c
c        ==================
c        PBE non-local part
c        ==================
         if(abs(nepsc).lt.tol_rho*tol_rho) goto 20
c
c        epsilon_c = n*epsilon_c / n
c
         epsc = nepsc/rhoval
         if (ipol.eq.1) then
            depscdna = dnepscdn(D1_RA)/rhoval-nepsc/(rhoval**2)
            depscdnb = depscdna
         else
            depscdna = dnepscdn(D1_RA)/rhoval-nepsc/(rhoval**2)
            depscdnb = dnepscdn(D1_RB)/rhoval-nepsc/(rhoval**2)
         endif
#ifdef SECOND_DERIV
         if (ipol.eq.1) then
            d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnb2   = d2epscdna2
         else
            d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
     &                     -dnepscdn(D1_RA)/(rhoval**2)
     &                     -dnepscdn(D1_RB)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnb2   = d2nepscdn2(D2_RB_RB)/rhoval
     &                     -dnepscdn(D1_RB)/(rhoval**2)
     &                     -dnepscdn(D1_RB)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
         endif
#endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
         dksdna = (1.0d0/6.0d0)*ks/rhoval
         dksdnb = dksdna
#ifdef SECOND_DERIV
         d2ksdna2   = (1.0d0/6.0d0)*dksdna/rhoval
     &              - (1.0d0/6.0d0)*ks/(rhoval**2)
         d2ksdnadnb = d2ksdna2
         d2ksdnb2   = d2ksdna2
#endif
c
c        zeta = (na-nb)/n and its derivs
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,R_A)-rho(n,R_B))/rhoval
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
         if (ipol.eq.1) then
            dzetadna = 1.0d0/rhoval
            dzetadnb = -1.0d0/rhoval
#ifdef SECOND_DERIV
            d2zetadna2   = -2.0d0/(rhoval**2)
            d2zetadnadnb = 0.0d0
            d2zetadnb2   = 2.0d0/(rhoval**2)
#endif
         else
            dzetadna =  2.0d0*rho(n,R_B)/(rhoval**2)
            dzetadnb = -2.0d0*rho(n,R_A)/(rhoval**2)
#ifdef SECOND_DERIV
            d2zetadna2   = -4.0d0*rho(n,R_B)/(rhoval**3)
            d2zetadnadnb = 2.0d0*(rho(n,R_A)-rho(n,R_B))/(rhoval**3)
            d2zetadnb2   = 4.0d0*rho(n,R_A)/(rhoval**3)
#endif
         endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
         phi = 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0)
     &               +(1.0d0-zeta)**(2.0d0/3.0d0))
         if ((1.0d0-zeta).lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &             (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta))
         else if ((1.0d0+zeta).lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &            -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
         else
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &         (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta)
     &        -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
         endif
         dphidna = dphidzeta*dzetadna
         dphidnb = dphidzeta*dzetadnb
#ifdef SECOND_DERIV
         if ((1.0d0-zeta).lt.tol_rho) then
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         (1.0d0+zeta)**(2.0d0/3.0d0)/((1.0d0+zeta)**2))
         else if ((1.0d0+zeta).lt.tol_rho) then
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         (1.0d0-zeta)**(2.0d0/3.0d0)/((1.0d0-zeta)**2))
         else
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         (1.0d0+zeta)**(2.0d0/3.0d0)/((1.0d0+zeta)**2)
     &        +(1.0d0-zeta)**(2.0d0/3.0d0)/((1.0d0-zeta)**2))
         endif
         d2phidna2   = d2phidzeta2*dzetadna*dzetadna
     &               + dphidzeta*d2zetadna2 
         d2phidnadnb = d2phidzeta2*dzetadna*dzetadnb
     &               + dphidzeta*d2zetadnadnb
         d2phidnb2   = d2phidzeta2*dzetadnb*dzetadnb
     &               + dphidzeta*d2zetadnb2 
#endif
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
         t = dsqgamma/(2.0d0*phi*ks*rhoval)
         dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
         dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
#ifdef SECOND_DERIV
         d2tdna2 = - dtdna/rhoval
     &           + t/(rhoval**2)
     &           - dtdna/phi*dphidna
     &           + t/(phi**2)*(dphidna**2)
     &           - t/phi*d2phidna2
     &           - dtdna/ks*dksdna
     &           + t/(ks**2)*(dksdna**2)
     &           - t/ks*d2ksdna2
         d2tdnadnb = - dtdnb/rhoval
     &           + t/(rhoval**2)
     &           - dtdnb/phi*dphidna
     &           + t/(phi**2)*(dphidna*dphidnb)
     &           - t/phi*d2phidnadnb
     &           - dtdnb/ks*dksdna
     &           + t/(ks**2)*(dksdna*dksdnb)
     &           - t/ks*d2ksdnadnb
         d2tdnb2 = - dtdnb/rhoval
     &           + t/(rhoval**2)
     &           - dtdnb/phi*dphidnb
     &           + t/(phi**2)*(dphidnb**2)
     &           - t/phi*d2phidnb2
     &           - dtdnb/ks*dksdnb
     &           + t/(ks**2)*(dksdnb**2)
     &           - t/ks*d2ksdnb2
#endif
c
c        { ... } in A (see below) and its derivs
c
         argexp = -epsc/GAMMA/(phi**3)
         dargexpdna = -depscdna/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidna
         dargexpdnb = -depscdnb/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidnb
#ifdef SECOND_DERIV
         d2argexpdna2 = -d2epscdna2/GAMMA/(phi**3)
     &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidna
     &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidna
     &        -12.0d0*epsc/GAMMA/(phi**5)*dphidna**2
     &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidna2
         d2argexpdnadnb = -d2epscdnadnb/GAMMA/(phi**3)
     &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidnb
     &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidna
     &        -12.0d0*epsc/GAMMA/(phi**5)*dphidna*dphidnb
     &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidnadnb
         d2argexpdnb2 = -d2epscdnb2/GAMMA/(phi**3)
     &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidnb
     &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidnb
     &        -12.0d0*epsc/GAMMA/(phi**5)*dphidnb**2
     &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidnb2
#endif
c
c        A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
         if (dabs(argexp).lt.EXPTOL) then
            expinA=dexp(argexp)
         else
            expinA=0.0d0
         endif
         A = BETA/GAMMA/(expinA-1.0d0)
         dAdna = -BETA/GAMMA*dargexpdna*expinA/(expinA-1.0d0)**2
         dAdnb = -BETA/GAMMA*dargexpdnb*expinA/(expinA-1.0d0)**2
#ifdef SECOND_DERIV
         d2Adna2   = -BETA/GAMMA*d2argexpdna2
     &               *expinA/(expinA-1.0d0)**2 
     &             - BETA/GAMMA*dargexpdna
     &               *dargexpdna*expinA/(expinA-1.0d0)**2
     &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdna
     &               *expinA*expinA/(expinA-1.0d0)**3
         d2Adnadnb  = -BETA/GAMMA*d2argexpdnadnb
     &               *expinA/(expinA-1.0d0)**2 
     &             - BETA/GAMMA*dargexpdna
     &               *dargexpdnb*expinA/(expinA-1.0d0)**2
     &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdnb
     &               *expinA*expinA/(expinA-1.0d0)**3
         d2Adnb2   = -BETA/GAMMA*d2argexpdnb2
     &               *expinA/(expinA-1.0d0)**2 
     &             - BETA/GAMMA*dargexpdnb
     &               *dargexpdnb*expinA/(expinA-1.0d0)**2
     &             + 2.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb
     &               *expinA*expinA/(expinA-1.0d0)**3
#endif
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
         fAtnum = 1.0d0+A*t**2
         fAtden = 1.0d0+A*t**2+A**2*t**4
         fAt = fAtnum/fAtden
         dfAtnumdt = 2.0d0*A*t
         dfAtnumdA = t**2
         dfAtdendt = 2.0d0*A*t+4.0d0*A**2*t**3
         dfAtdendA = t**2+2.0d0*A*t**4
         dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/(fAtden**2)
         dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/(fAtden**2)
         dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
         dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
#ifdef SECOND_DERIV
         d2fAtnumdt2 = 2.0d0*A
         d2fAtdendt2 = 2.0d0*A+12.0d0*A**2*t**2
         d2fAtnumdtdA = 2.0d0*t
         d2fAtdendtdA = 2.0d0*t+8.0d0*A*t**3
         d2fAtnumdA2 = 0.0d0
         d2fAtdendA2 = 2.0d0*t**4
         d2fAtdt2  = (d2fAtnumdt2*fAtden-fAtnum*d2fAtdendt2)
     &               /(fAtden**2)
     &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
     &               /(fAtden**3)*dfAtdendt
         d2fAtdtdA = (d2fAtnumdtdA*fAtden+dfAtnumdt*dfAtdendA
     &               -dfAtnumdA*dfAtdendt-fAtnum*d2fAtdendtdA)
     &               /(fAtden**2)
     &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
     &               /(fAtden**3)*dfAtdendA
         d2fAtdA2  = (d2fAtnumdA2*fAtden-fAtnum*d2fAtdendA2)
     &               /(fAtden**2)
     &               -2.0d0*(dfAtnumdA*fAtden-fAtnum*dfAtdendA)
     &               /(fAtden**3)*dfAtdendA
         d2fAtdna2 = d2fAtdt2*dtdna*dtdna + d2fAtdtdA*dtdna*dAdna
     &             + dfAtdt*d2tdna2 + d2fAtdtdA*dAdna*dtdna
     &             + d2fAtdA2*dAdna*dAdna + dfAtdA*d2Adna2
         d2fAtdnb2 = d2fAtdt2*dtdnb*dtdnb + d2fAtdtdA*dtdnb*dAdnb
     &             + dfAtdt*d2tdnb2 + d2fAtdtdA*dAdnb*dtdnb
     &             + d2fAtdA2*dAdnb*dAdnb + dfAtdA*d2Adnb2
         d2fAtdnadnb = d2fAtdt2*dtdna*dtdnb + d2fAtdtdA*dtdna*dAdnb
     &             + dfAtdt*d2tdnadnb + d2fAtdtdA*dAdna*dtdnb
     &             + d2fAtdA2*dAdna*dAdnb + dfAtdA*d2Adnadnb
#endif
c
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + BETA/GAMMA*t**2*fAt
         darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
     &                            +t*t*dfAtdna)
         darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
     &                            +t*t*dfAtdnb)
#ifdef SECOND_DERIV
         d2arglogdna2 = BETA/GAMMA*(2.0d0*dtdna*dtdna*fAt
     &                             +2.0d0*t*d2tdna2*fAt
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +t*t*d2fAtdna2)
         d2arglogdnb2 = BETA/GAMMA*(2.0d0*dtdnb*dtdnb*fAt
     &                             +2.0d0*t*d2tdnb2*fAt
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +t*t*d2fAtdnb2)
         d2arglogdnadnb = BETA/GAMMA*(2.0d0*dtdna*dtdnb*fAt
     &                             +2.0d0*t*d2tdnadnb*fAt
     &                             +2.0d0*t*dtdna*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdna
     &                             +t*t*d2fAtdnadnb)
#endif
c
c        H0 = GAMMA * phi**3 * log{arglog} and its derivs
c
         H0 = GAMMA*(phi**3)*dlog(arglog)
         dH0dna = GAMMA*(3.0d0*(phi**2)*dphidna*dlog(arglog)
     &                  +(phi**3)*darglogdna/arglog)
         dH0dnb = GAMMA*(3.0d0*(phi**2)*dphidnb*dlog(arglog)
     &                  +(phi**3)*darglogdnb/arglog)
#ifdef SECOND_DERIV
         d2H0dna2 = GAMMA*(6.0d0*phi*dphidna*dphidna*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidna2*dlog(arglog)
     &                +6.0d0*(phi**2)*dphidna*darglogdna/arglog
     &                +(phi**3)*d2arglogdna2/arglog
     &                -(phi**3)*darglogdna*darglogdna/arglog/arglog)
         d2H0dnadnb = GAMMA*(6.0d0*phi*dphidna*dphidnb*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidnadnb*dlog(arglog)
     &                +3.0d0*(phi**2)*dphidna*darglogdnb/arglog
     &                +3.0d0*(phi**2)*dphidnb*darglogdna/arglog
     &                +(phi**3)*d2arglogdnadnb/arglog
     &                -(phi**3)*darglogdna*darglogdnb/arglog/arglog)
         d2H0dnb2 = GAMMA*(6.0d0*phi*dphidnb*dphidnb*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidnb2*dlog(arglog)
     &                +6.0d0*(phi**2)*dphidnb*darglogdnb/arglog
     &                +(phi**3)*d2arglogdnb2/arglog
     &                -(phi**3)*darglogdnb*darglogdnb/arglog/arglog)
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
         ffunc(n)=ffunc(n)+wght*(nepsc + H0*rhoval)
         Amat(n,D1_RA) = Amat(n,D1_RA)
     &                 +  (dnepscdn(D1_RA) + H0 + rhoval*dH0dna)*wght
         if (ipol.eq.2) then
           Amat(n,D1_RB) = Amat(n,D1_RB)
     &                   +  (dnepscdn(D1_RB) + H0 + rhoval*dH0dnb)*wght
         endif
#ifdef SECOND_DERIV
         Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &   + (d2nepscdn2(D2_RA_RA) + 2.d0*dH0dna + rhoval*d2H0dna2)*wght
         Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &   + (d2nepscdn2(D2_RA_RB) + dH0dna + dH0dnb
     &      + rhoval*d2H0dnadnb)*wght
         if (ipol.eq.2) then
           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &     + (d2nepscdn2(D2_RB_RB) + 2.d0*dH0dnb + rhoval*d2H0dnb2)*wght
         endif
#endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
         if (dsqgamma.gt.TOLL)then
            dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
            dfAtdg = dfAtdt*dtdg
            darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
            dH0dg = GAMMA*(phi**3)*darglogdg/arglog
            if (ipol.eq.1) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
            else
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*wght
            endif
#ifdef SECOND_DERIV
            d2tdg2 = -0.125d0/(phi*ks*rhoval)/(dsqgamma**3)
            d2tdnadg = -dtdg/rhoval-dtdg/phi*dphidna
     &                 -dtdg/ks*dksdna
            d2tdnbdg = -dtdg/rhoval-dtdg/phi*dphidnb
     &                 -dtdg/ks*dksdnb
            d2fAtdg2 = d2fAtdt2*(dtdg**2)+dfAtdt*d2tdg2
            d2fAtdnadg = d2fAtdt2*dtdg*dtdna
     &                  +d2fAtdtdA*dtdg*dAdna
     &                  +dfAtdt*d2tdnadg 
            d2fAtdnbdg = d2fAtdt2*dtdg*dtdnb
     &                  +d2fAtdtdA*dtdg*dAdnb
     &                  +dfAtdt*d2tdnbdg 
            d2arglogdnadg = BETA/GAMMA*(2.0d0*dtdna*dtdg*fAt
     &                                 +2.0d0*t*d2tdnadg*fAt
     &                                 +2.0d0*t*dtdg*dfAtdna
     &                                 +2.0d0*t*dtdna*dfAtdg
     &                                 +t*t*d2fAtdnadg)
            d2arglogdnbdg = BETA/GAMMA*(2.0d0*dtdnb*dtdg*fAt
     &                                 +2.0d0*t*d2tdnbdg*fAt
     &                                 +2.0d0*t*dtdg*dfAtdnb
     &                                 +2.0d0*t*dtdnb*dfAtdg
     &                                 +t*t*d2fAtdnbdg)
            d2arglogdg2 = BETA/GAMMA*(2.0d0*dtdg*dtdg*fAt
     &                               +2.0d0*t*d2tdg2*fAt
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +t*t*d2fAtdg2)
            d2H0dnadg = GAMMA*3.0d0*phi**2*dphidna*darglogdg/arglog 
     &                + GAMMA*phi**3*d2arglogdnadg/arglog
     &                - GAMMA*phi**3*darglogdg*darglogdna/arglog**2
            d2H0dnbdg = GAMMA*3.0d0*phi**2*dphidnb*darglogdg/arglog 
     &                + GAMMA*phi**3*d2arglogdnbdg/arglog
     &                - GAMMA*phi**3*darglogdg*darglogdnb/arglog**2
            d2H0dg2 = GAMMA*phi**3*d2arglogdg2/arglog
     &              - GAMMA*phi**3*darglogdg*darglogdg/arglog**2
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &             + (dH0dg + d2H0dnadg*rhoval)*wght
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
     &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*wght
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
     &             + (dH0dg + d2H0dnadg*rhoval)*wght
            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
     &             + (dH0dg + d2H0dnbdg*rhoval)*wght
            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
     &             + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &             + (dH0dg + d2H0dnbdg*rhoval)*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &             + d2H0dg2*rhoval*wght
            Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
     &             + 2.0d0*d2H0dg2*rhoval*wght
            Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
     &             + d2H0dg2*rhoval*wght
            Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
     &             + 4.0d0*d2H0dg2*rhoval*wght
            Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
     &             + 2.0d0*d2H0dg2*rhoval*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &             + d2H0dg2*rhoval*wght
#endif
         endif
   20 continue
c
      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_c_pbe.F"
#endif
