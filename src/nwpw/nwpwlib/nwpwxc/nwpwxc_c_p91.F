#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_p91.F
C> The Perdew 91 correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The Perdew 91 correlation functional
C>
C> The Perdew 91 correlation functional [1-3].
C>
C> [1] M. Rasolt, D.J.W. Geldart,
C>     "Exchange and correlation energy in a nonuniform fermion fluid",
C>     Phys. Rev. B <b>34</b>, 1325-1328 (1986), DOI:
C>     <A HREF="https://doi.org/10.1103/PhysRevB.34.1325">
C>     10.1103/PhysRevB.34.1325</A>.
C>
C> [2] J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
C>     M.R. Pederson, D.J. Singh, C. Fiolhais, 
C>     "Atoms, molecules, solids, and surfaces: Applications of the
C>     generalized gradient approximation for exchange and correlation",
C>     Phys. Rev. B <b>46</b>, 6671-6687 (1992), DOI:
C>     <A HREF="https://doi.org/10.1103/PhysRevB.46.6671">
C>     10.1103/PhysRevB.46.6671</A>.
C>
C> [3] J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
C>     M.R. Pederson, D.J. Singh, C. Fiolhais, 
C>     "Erratum: Atoms, molecules, solids, and surfaces: Applications
C>     of the generalized gradient approximation for exchange and
C>     correlation",
C>     Phys. Rev. B <b>48</b>, 4978-4978 (1993), DOI:
C>     <A HREF="https://doi.org/10.1103/PhysRevB.48.4978.2">
C>     10.1103/PhysRevB.48.4978.2</A>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_p91(tol_rho,ipol,nq,wght,rho,rgamma,ffunc,
     &                      Amat, Cmat)
#else
      Subroutine nwpwxc_c_p91_d2(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc, Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id$
c
      Implicit none
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
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
      double precision ALPHA, BETA, NU, CCZERO, CX, PI
      double precision cnoname, ca, cb, cc, cd
      parameter (ALPHA = 0.09d0)
      parameter (CCZERO = 0.004235d0)
      parameter (CX = -0.001667d0)
      parameter (PI = 3.1415926535897932385d0)
      parameter (cnoname = 2.568d0)
      parameter (ca = 23.266d0)
      parameter (cb = 7.389d-3)
      parameter (cc = 8.723d0)
      parameter (cd = 0.472d0)
c
      integer lnq
      double precision lqwght
c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
      double precision rho_t(3)
      double precision dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n
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
      double precision rs, drsdna, drsdnb
      double precision ccrs, dccrsdna, dccrsdnb
      double precision cnum, dcnumdna, dcnumdnb
      double precision cden, dcdendna, dcdendnb
      double precision kf, dkfdna, dkfdnb
      double precision H1argexp, dH1argexpdna, dH1argexpdnb
      double precision dH1argexpdg
      double precision H1prefac, dH1prefacdna, dH1prefacdnb
      double precision dH1prefacdg
      double precision expinH1
      double precision H1, dH1dna, dH1dnb, dH1dg
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
      double precision d2rsdna2, d2rsdnb2, d2rsdnadnb
      double precision d2ccrsdna2, d2ccrsdnb2, d2ccrsdnadnb
      double precision d2cnumdna2, d2cnumdnb2, d2cnumdnadnb
      double precision d2cdendna2, d2cdendnb2, d2cdendnadnb
      double precision d2kfdna2, d2kfdnb2, d2kfdnadnb
      double precision d2H1argexpdna2, d2H1argexpdnadnb, d2H1argexpdnb2
      double precision d2H1argexpdnadg, d2H1argexpdnbdg, d2H1argexpdg2
      double precision d2H1prefacdna2, d2H1prefacdnadnb, d2H1prefacdnb2
      double precision d2H1prefacdnadg, d2H1prefacdnbdg, d2H1prefacdg2
      double precision d2H1dna2, d2H1dnb2, d2H1dnadnb
      double precision d2H1dnadg, d2H1dnbdg, d2H1dg2
#endif
      logical lfac
      parameter (lfac = .false.)
c     
c References:
c [a] J. P. Perdew et al., Phys. Rev. B 46, 6671 (1992).
c [b] M. Rasolt and D. J. W. Geldart, Phys. Rev. B 34, 1325 (1986).
c     
c  E_c(PW91) = Int n (epsilon_c + H0 + H1) dxdydz
c
c  n*epsilon_c                <=== supplied by another subroutine
c  d(n*epsilon_c)/d(na)       <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(na)d(na) <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(na)d(nb) <=== supplied by another subroutine
c  d2(n*epsilon_c)/d(nb)d(nb) <=== supplied by another subroutine
c
c  H0 = BETA**2/2/ALPHA * phi**3 * log{ 1 + 2*ALPHA/BETA * t**2 * [ ... ]}
c
c  BETA = NU * CCZERO
c
c  NU = (16/PI)(3*PI**2)**(1/3)
c
c  phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)]
c
c  zeta = (na - nb)/n
c
c  [ ... ] = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4)
c
c  A = 2*ALPHA/BETA [exp{-2*ALPHA*epsilon_c/(BETA**2*phi**3)}-1]**(-1)
c
c  t = |Nabla n|/(2*phi*ks*n)
c
c  ks = 2 * (3 * PI**2 * n)**(1/6) / sqrt(PI)
c
c  |Nabla n| = sqrt(g_aa + g_bb + 2*g_ab)
c
c  H1 = NU * [CC(rs) - CCZERO - 3/7*CX] * Phi**3 * t**2 
c       * exp[-100 * phi**4 * (ks**2/kf**2) * t**2]
c
c  CC(rs) = 0.001 * (cnoname + ca*rs + cb*rs**2)/(1+cc*rs+cd*rs**2+10*cb*rs**3)
c         - CX
c
c  rs = (3/4 / PI / n)**(1/3)
c
c  Names of variables
c
c  E_c(PW91)                 : Ec
c  n (alpha+beta density)    : rhoval
c  na, nb                    : rho(*,2), rho(*,3)
c  epsilon_c                 : epsc
c  H0                        : H0
c  n*epsilon_c               : nepsc
c  phi (also called "g")     : phi
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
c  H1                        : H1
c  CC(rs)                    : ccrs
c  CCZERO                    : CCZERO
c  CX                        : CX
c  kf                        : kf
c  rs                        : rs
c
c  Derivatives of these are named like d...dna, d2...dnadnb,
c  d2...dna2, etc.
c
      NU = (16.0d0/PI)*(3.0d0*PI**2)**(1.0d0/3.0d0)
      BETA = NU * CCZERO
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      do 20 n = 1, nq
c
c        n and zeta = (na - nb)/n
c
         rhoval = rho(n,R_T)
         rho_t(1) = rho(n,R_T)
         if (ipol.eq.2) then
            rhoval = rho(n,R_A)+rho(n,R_B)
            rho_t(1) = rho(n,R_A)
            rho_t(2) = rho(n,R_B)
         endif
         if (rhoval.le.tol_rho) goto 20
         if (ipol.eq.1) then
c           gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1)
            gammaval = rgamma(n,G_TT)
         else
c           gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,1,2)*delrho(n,1,2) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,2,2)*delrho(n,2,2) +
c    &                 delrho(n,3,1)*delrho(n,3,1) +
c    &                 delrho(n,3,2)*delrho(n,3,2) +
c    &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &                 delrho(n,2,1)*delrho(n,2,2) +
c    &                 delrho(n,3,1)*delrho(n,3,2))
            gammaval = rgamma(n,G_AA)+rgamma(n,G_BB)+
     &                 rgamma(n,G_AB)*2
         endif
         dsqgamma = max(dsqrt(gammaval),tol_rho)
         nepsc = 0.0d0
         dnepscdn(1) = 0.0d0
         if (ipol.eq.2) dnepscdn(2) = 0.0d0
#ifdef SECOND_DERIV
         d2nepscdn2(D2_RA_RA)=0.0d0
         d2nepscdn2(D2_RA_RB)=0.0d0
         if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
c
c        call for LDA bit
c
         lnq = 1
         lqwght = 1.0d0
c
#ifndef SECOND_DERIV
         call nwpwxc_c_pw91lda(tol_rho,ipol,lnq,lqwght,rho_t,nepsc,
     &                       dnepscdn)
#else
         call nwpwxc_c_pw91lda_d2(tol_rho,ipol,lnq,lqwght,rho_t,nepsc,
     &                          dnepscdn,d2nepscdn2)
#endif
c
c        ============
c        PW91 H0 part
c        ============
         if(abs(nepsc).lt.tol_rho*tol_rho) goto 20
c
c        epsilon_c = n*epsilon_c / n
c
         epsc = nepsc/rhoval
         if (ipol.eq.1) then
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = depscdna
         else
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = dnepscdn(2)/rhoval-nepsc/(rhoval**2)
         endif
#ifdef SECOND_DERIV
         if (ipol.eq.1) then
            d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
     &                     -dnepscdn(1)/(rhoval**2)
     &                     -dnepscdn(1)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
     &                     -dnepscdn(1)/(rhoval**2)
     &                     -dnepscdn(1)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnb2   = d2epscdna2
         else
            d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
     &                     -dnepscdn(1)/(rhoval**2)
     &                     -dnepscdn(1)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
     &                     -dnepscdn(1)/(rhoval**2)
     &                     -dnepscdn(2)/(rhoval**2)
     &                     +2.0d0*nepsc/(rhoval**3)
            d2epscdnb2   = d2nepscdn2(D2_RB_RB)/rhoval
     &                     -dnepscdn(2)/(rhoval**2)
     &                     -dnepscdn(2)/(rhoval**2)
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
            d2zetadnb2   = -2.0d0/(rhoval**2)
#endif
         else
            dzetadna =  2.0d0*rho(n,R_B)/(rhoval**2)
            dzetadnb = -2.0d0*rho(n,R_A)/(rhoval**2)
#ifdef SECOND_DERIV
            d2zetadna2   = -4.0d0*rho(n,R_B)/(rhoval**3)
            d2zetadnadnb = 2.0d0*(rho(n,R_A)-rho(n,R_B))/(rhoval**3)
            d2zetadnb2   = -4.0d0*rho(n,R_A)/(rhoval**3)
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
         d2tdnb2 = - dtdna/rhoval
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
         argexp = -2.0d0*ALPHA*epsc/BETA**2/(phi**3)
         dargexpdna = -2.0d0*ALPHA*depscdna/BETA**2/(phi**3)
     &                +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*dphidna
         dargexpdnb = -2.0d0*ALPHA*depscdnb/BETA**2/(phi**3)
     &                +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*dphidnb
#ifdef SECOND_DERIV
         d2argexpdna2 = -2.0d0*ALPHA*d2epscdna2/BETA**2/(phi**3)
     &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidna
     &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidna
     &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidna**2
     &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidna2
         d2argexpdnadnb = -2.0d0*ALPHA*d2epscdnadnb/BETA**2/(phi**3)
     &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidnb
     &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidna
     &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidna*dphidnb
     &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidnadnb
         d2argexpdnb2 = -2.0d0*ALPHA*d2epscdnb2/BETA**2/(phi**3)
     &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidnb
     &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidnb
     &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidnb**2
     &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidnb2
#endif
c
c        A = 2*ALPHA/BETA [exp{-2*ALPHA*epsilon_c/(BETA**2*phi**3)}-1]**(-1)
c
         if (dabs(argexp).lt.EXPTOL) then
            expinA=dexp(argexp)
         else
            expinA=0.0d0
         endif
         A = 2.0d0*ALPHA/BETA/(expinA-1.0d0)
         dAdna = -2.0d0*ALPHA/BETA*dargexpdna*expinA/(expinA-1.0d0)**2
         dAdnb = -2.0d0*ALPHA/BETA*dargexpdnb*expinA/(expinA-1.0d0)**2
#ifdef SECOND_DERIV
         d2Adna2   = -2.0d0*ALPHA/BETA*d2argexpdna2
     &               *expinA/(expinA-1.0d0)**2 
     &             - 2.0d0*ALPHA/BETA*dargexpdna
     &               *dargexpdna*expinA/(expinA-1.0d0)**2
     &             + 4.0d0*ALPHA/BETA*dargexpdna*dargexpdna
     &               *expinA*expinA/(expinA-1.0d0)**3
         d2Adnadnb  = -2.0d0*ALPHA/BETA*d2argexpdnadnb
     &               *expinA/(expinA-1.0d0)**2 
     &             - 2.0d0*ALPHA/BETA*dargexpdna
     &               *dargexpdnb*expinA/(expinA-1.0d0)**2
     &             + 4.0d0*ALPHA/BETA*dargexpdna*dargexpdnb
     &               *expinA*expinA/(expinA-1.0d0)**3
         d2Adnb2   = -2.0d0*ALPHA/BETA*d2argexpdnb2
     &               *expinA/(expinA-1.0d0)**2 
     &             - 2.0d0*ALPHA/BETA*dargexpdnb
     &               *dargexpdnb*expinA/(expinA-1.0d0)**2
     &             + 4.0d0*ALPHA/BETA*dargexpdnb*dargexpdnb
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
c        arglog = 1 + 2*ALPHA/BETA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + 2.0d0*ALPHA/BETA*t**2*fAt
         darglogdna = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdna*fAt
     &                            +t*t*dfAtdna)
         darglogdnb = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdnb*fAt
     &                            +t*t*dfAtdnb)
#ifdef SECOND_DERIV
         d2arglogdna2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdna*fAt
     &                             +2.0d0*t*d2tdna2*fAt
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +t*t*d2fAtdna2)
         d2arglogdnb2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdnb*dtdnb*fAt
     &                             +2.0d0*t*d2tdnb2*fAt
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +t*t*d2fAtdnb2)
         d2arglogdnadnb = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdnb*fAt
     &                             +2.0d0*t*d2tdnadnb*fAt
     &                             +2.0d0*t*dtdna*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdna
     &                             +t*t*d2fAtdnadnb)
#endif
c
c        H0 = BETA**2/2/ALPHA * phi**3 * log{arglog} and its derivs
c
         H0 = 0.5d0*BETA**2/ALPHA*(phi**3)*dlog(arglog)
         dH0dna = 0.5d0*BETA**2/ALPHA*(3.0d0*(phi**2)
     &                  *dphidna*dlog(arglog)
     &                  +(phi**3)*darglogdna/arglog)
         dH0dnb = 0.5d0*BETA**2/ALPHA*(3.0d0*(phi**2)
     &                  *dphidnb*dlog(arglog)
     &                  +(phi**3)*darglogdnb/arglog)
#ifdef SECOND_DERIV
         d2H0dna2 = 0.5d0*BETA**2/ALPHA
     &                *(6.0d0*phi*dphidna*dphidna*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidna2*dlog(arglog)
     &                +6.0d0*(phi**2)*dphidna*darglogdna/arglog
     &                +(phi**3)*d2arglogdna2/arglog
     &                -(phi**3)*darglogdna*darglogdna/arglog/arglog)
         d2H0dnadnb = 0.5d0*BETA**2/ALPHA
     &                *(6.0d0*phi*dphidna*dphidnb*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidnadnb*dlog(arglog)
     &                +3.0d0*(phi**2)*dphidna*darglogdnb/arglog
     &                +3.0d0*(phi**2)*dphidnb*darglogdna/arglog
     &                +(phi**3)*d2arglogdnadnb/arglog
     &                -(phi**3)*darglogdna*darglogdnb/arglog/arglog)
         d2H0dnb2 = 0.5d0*BETA**2/ALPHA
     &                *(6.0d0*phi*dphidnb*dphidnb*dlog(arglog)
     &                +3.0d0*(phi**2)*d2phidnb2*dlog(arglog)
     &                +6.0d0*(phi**2)*dphidnb*darglogdnb/arglog
     &                +(phi**3)*d2arglogdnb2/arglog
     &                -(phi**3)*darglogdnb*darglogdnb/arglog/arglog)
#endif
c        ============
c        PW91 H1 part
c        ============
c
c        rs = (3/4 / PI / n)**(1/3) and its derivs
c        
         rs = (0.75d0 / PI / dabs(rhoval))**(1.0d0/3.0d0)
         drsdna = (-1.0d0/3.0d0)*rs/rhoval
         drsdnb = drsdna
#ifdef SECOND_DERIV
         d2rsdna2 = (-4.0d0/3.0d0)*drsdna/rhoval
         d2rsdnb2 = d2rsdna2
         d2rsdnadnb = d2rsdna2
#endif
c
c        CC(rs) and its derivs
c
         cnum = cnoname + ca*rs + cb*rs**2
         cden = 1.0d0 + cc*rs + cd*rs**2 + 10.0d0*cb*rs**3
         ccrs = 1.0d-3*cnum/cden - CX
         dcnumdna = (ca + 2.0d0*cb*rs)*drsdna
         dcnumdnb = (ca + 2.0d0*cb*rs)*drsdnb
         dcdendna = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*drsdna
         dcdendnb = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*drsdnb
         dccrsdna = 1.0d-3*(dcnumdna*cden-cnum*dcdendna)/cden**2
         dccrsdnb = 1.0d-3*(dcnumdnb*cden-cnum*dcdendnb)/cden**2
#ifdef SECOND_DERIV
         d2cnumdna2 = (ca + 2.0d0*cb*rs)*d2rsdna2
     &              + (2.0d0*cb)*drsdna*drsdna
         d2cnumdnb2 = (ca + 2.0d0*cb*rs)*d2rsdnb2
     &              + (2.0d0*cb)*drsdnb*drsdnb
         d2cnumdnadnb = (ca + 2.0d0*cb*rs)*d2rsdnadnb
     &              + (2.0d0*cb)*drsdna*drsdnb
         d2cdendna2 = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdna2
     &              + (2.0d0*cd+60.0d0*cb*rs)*drsdna*drsdna
         d2cdendnb2 = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdnb2
     &              + (2.0d0*cd+60.0d0*cb*rs)*drsdnb*drsdnb
         d2cdendnadnb = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdnadnb
     &              + (2.0d0*cd+60.0d0*cb*rs)*drsdna*drsdnb
         d2ccrsdna2 = 1.0d-3*(d2cnumdna2*cden-cnum*d2cdendna2)/cden**2
     &      - 1.0d-3*2.0d0*(dcnumdna*cden-cnum*dcdendna)
     &        /cden**3*dcdendna
         d2ccrsdnb2 = 1.0d-3*(d2cnumdnb2*cden-cnum*d2cdendnb2)/cden**2
     &      - 1.0d-3*2.0d0*(dcnumdnb*cden-cnum*dcdendnb)
     &        /cden**3*dcdendnb
         d2ccrsdnadnb 
     &      = 1.0d-3*(d2cnumdnadnb*cden-cnum*d2cdendnadnb)/cden**2
     &      + 1.0d-3*(dcnumdna*dcdendnb-dcnumdnb*dcdendna)/cden**2
     &      - 1.0d-3*2.0d0*(dcnumdna*cden-cnum*dcdendna)
     &        /cden**3*dcdendnb
#endif
c
c        kf = (3 * PI**2 * n)**(1/3)
c
         kf = (3.0d0*PI**2*rhoval)**(1.0d0/3.0d0)
         dkfdna = (1.0d0/3.0d0)*kf/rhoval
         dkfdnb = dkfdna
#ifdef SECOND_DERIV
         d2kfdna2 = (-2.0d0/3.0d0)*dkfdna/rhoval
         d2kfdnadnb = (-2.0d0/3.0d0)*dkfdna/rhoval
         d2kfdnb2 = (-2.0d0/3.0d0)*dkfdnb/rhoval
#endif
c
c        H1argexp = -100 * phi**4 * (ks**2/kf**2) * t**2 and its derivs
c
         H1argexp = -100.0d0*phi**4*(ks**2/kf**2)*t**2
         dH1argexpdna = 4.0d0*H1argexp/phi*dphidna
     &                 +2.0d0*H1argexp/ks*dksdna
     &                 -2.0d0*H1argexp/kf*dkfdna
     &                 +2.0d0*H1argexp/t*dtdna
         dH1argexpdnb = 4.0d0*H1argexp/phi*dphidnb
     &                 +2.0d0*H1argexp/ks*dksdnb
     &                 -2.0d0*H1argexp/kf*dkfdnb
     &                 +2.0d0*H1argexp/t*dtdnb
#ifdef SECOND_DERIV
         d2H1argexpdna2 = 4.0d0*dH1argexpdna/phi*dphidna
     &                   -4.0d0*H1argexp/phi**2*dphidna**2
     &                   +4.0d0*H1argexp/phi*d2phidna2
     &                   +2.0d0*dH1argexpdna/ks*dksdna
     &                   -2.0d0*H1argexp/ks**2*dksdna**2
     &                   +2.0d0*H1argexp/ks*d2ksdna2
     &                   -2.0d0*dH1argexpdna/kf*dkfdna
     &                   +2.0d0*H1argexp/kf**2*dkfdna**2
     &                   -2.0d0*H1argexp/kf*d2kfdna2
     &                   +2.0d0*dH1argexpdna/t*dtdna
     &                   -2.0d0*H1argexp/t**2*dtdna**2
     &                   +2.0d0*H1argexp/t*d2tdna2
         d2H1argexpdnb2 = 4.0d0*dH1argexpdnb/phi*dphidnb
     &                   -4.0d0*H1argexp/phi**2*dphidnb**2
     &                   +4.0d0*H1argexp/phi*d2phidnb2
     &                   +2.0d0*dH1argexpdnb/ks*dksdnb
     &                   -2.0d0*H1argexp/ks**2*dksdnb**2
     &                   +2.0d0*H1argexp/ks*d2ksdnb2
     &                   -2.0d0*dH1argexpdnb/kf*dkfdnb
     &                   +2.0d0*H1argexp/kf**2*dkfdnb**2
     &                   -2.0d0*H1argexp/kf*d2kfdnb2
     &                   +2.0d0*dH1argexpdnb/t*dtdnb
     &                   -2.0d0*H1argexp/t**2*dtdnb**2
     &                   +2.0d0*H1argexp/t*d2tdnb2
         d2H1argexpdnadnb = 4.0d0*dH1argexpdna/phi*dphidnb
     &                   -4.0d0*H1argexp/phi**2*dphidna*dphidnb
     &                   +4.0d0*H1argexp/phi*d2phidnadnb
     &                   +2.0d0*dH1argexpdna/ks*dksdnb
     &                   -2.0d0*H1argexp/ks**2*dksdna*dksdnb
     &                   +2.0d0*H1argexp/ks*d2ksdnadnb
     &                   -2.0d0*dH1argexpdna/kf*dkfdnb
     &                   +2.0d0*H1argexp/kf**2*dkfdna*dkfdnb
     &                   -2.0d0*H1argexp/kf*d2kfdnadnb
     &                   +2.0d0*dH1argexpdna/t*dtdnb
     &                   -2.0d0*H1argexp/t**2*dtdna*dtdnb
     &                   +2.0d0*H1argexp/t*d2tdnadnb
#endif
c
c        H1prefac = NU*[CC(rs) - CC(0)-3/7*CX] * g**3 * t**2
c
         H1prefac = NU*(ccrs - CCZERO - (3.0d0/7.0d0)*CX)
     &            * phi**3 * t**2
         dH1prefacdna = NU*dccrsdna * phi**3 * t**2
     &                + 3.0d0*H1prefac/phi*dphidna
     &                + 2.0d0*H1prefac/t*dtdna
         dH1prefacdnb = NU*dccrsdnb * phi**3 * t**2
     &                + 3.0d0*H1prefac/phi*dphidnb
     &                + 2.0d0*H1prefac/t*dtdnb
#ifdef SECOND_DERIV
         d2H1prefacdna2 = NU*d2ccrsdna2 * phi**3 * t**2
     &                  + 3.0d0* NU*dccrsdna * phi**2*dphidna * t**2
     &                  + 2.0d0* NU*dccrsdna * phi**3 * t*dtdna
     &                  + 3.0d0*dH1prefacdna/phi*dphidna
     &                  - 3.0d0*H1prefac/phi**2*dphidna**2
     &                  + 3.0d0*H1prefac/phi*d2phidna2
     &                  + 2.0d0*dH1prefacdna/t*dtdna
     &                  - 2.0d0*H1prefac/t**2*dtdna**2
     &                  + 2.0d0*H1prefac/t*d2tdna2
         d2H1prefacdnb2 = NU*d2ccrsdnb2 * phi**3 * t**2
     &                  + 3.0d0* NU*dccrsdnb * phi**2*dphidnb * t**2
     &                  + 2.0d0* NU*dccrsdnb * phi**3 * t*dtdnb
     &                  + 3.0d0*dH1prefacdnb/phi*dphidnb
     &                  - 3.0d0*H1prefac/phi**2*dphidnb**2
     &                  + 3.0d0*H1prefac/phi*d2phidnb2
     &                  + 2.0d0*dH1prefacdnb/t*dtdnb
     &                  - 2.0d0*H1prefac/t**2*dtdnb**2
     &                  + 2.0d0*H1prefac/t*d2tdnb2
         d2H1prefacdnadnb = NU*d2ccrsdnadnb * phi**3 * t**2
     &                  + 3.0d0* NU*dccrsdna * phi**2*dphidnb * t**2
     &                  + 2.0d0* NU*dccrsdna * phi**3 * t*dtdnb
     &                  + 3.0d0*dH1prefacdna/phi*dphidnb
     &                  - 3.0d0*H1prefac/phi**2*dphidna*dphidnb
     &                  + 3.0d0*H1prefac/phi*d2phidnadnb
     &                  + 2.0d0*dH1prefacdna/t*dtdnb
     &                  - 2.0d0*H1prefac/t**2*dtdna*dtdnb
     &                  + 2.0d0*H1prefac/t*d2tdnadnb
#endif
c
c        H1 = H1prefac * exp(H1argexp)
c 
         if (dabs(H1argexp).lt.EXPTOL) then
            expinH1 = dexp(H1argexp)
         else
            expinH1 = 0.0d0
         endif
         H1 = H1prefac * expinH1
         dH1dna = dH1prefacdna * expinH1
     &          + H1prefac * dH1argexpdna * expinH1
         dH1dnb = dH1prefacdnb * expinH1
     &          + H1prefac * dH1argexpdnb * expinH1
#ifdef SECOND_DERIV
         d2H1dna2 = d2H1prefacdna2 * expinH1
     &            + dH1prefacdna * dH1argexpdna * expinH1
     &            + dH1prefacdna * dH1argexpdna * expinH1
     &            + H1prefac * d2H1argexpdna2 * expinH1
     &            + H1prefac * dH1argexpdna * dH1argexpdna * expinH1
         d2H1dnb2 = d2H1prefacdnb2 * expinH1
     &            + dH1prefacdnb * dH1argexpdnb * expinH1
     &            + dH1prefacdnb * dH1argexpdnb * expinH1
     &            + H1prefac * d2H1argexpdnb2 * expinH1
     &            + H1prefac * dH1argexpdnb * dH1argexpdnb * expinH1
         d2H1dnadnb = d2H1prefacdnadnb * expinH1
     &            + dH1prefacdna * dH1argexpdnb * expinH1
     &            + dH1prefacdnb * dH1argexpdna * expinH1
     &            + H1prefac * d2H1argexpdnadnb * expinH1
     &            + H1prefac * dH1argexpdna * dH1argexpdnb * expinH1
#endif
c
c        Now we update Ec, Amat, and Amat2
c
c        if (lfac) then
c           ffunc(n)=ffunc(n)+nepsc*wght               
c        endif
c        if (nlfac) then
            ffunc(n)=ffunc(n)+(H0*rhoval+H1*rhoval)*wght
c        endif
c        if (lfac) then
c           Amat(n,1) = Amat(n,1) + dnepscdn(1)*wght
c           if (ipol.eq.2) Amat(n,2) = Amat(n,2) + dnepscdn(2)*wght
#ifdef SECOND_DERIV
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) 
c    &                        + d2nepscdn2(D2_RA_RA)*wght
c           Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) 
c    &                        + d2nepscdn2(D2_RA_RB)*wght
c           if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) 
c    &                        + d2nepscdn2(D2_RB_RB)*wght
#endif
c        endif
c        if (nlfac)then
            Amat(n,D1_RA) = Amat(n,D1_RA) + (H0 + H1 +
     &         rhoval*dH0dna + rhoval*dH1dna)*wght
            if (ipol.eq.2) Amat(n,D1_RB) = Amat(n,D1_RB) + (H0 + H1 +
     &         rhoval*dH0dnb + rhoval*dH1dnb)*wght
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &         + (2.0d0*dH0dna + rhoval*d2H0dna2
     &         +  2.0d0*dH1dna + rhoval*d2H1dna2)*wght
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &         + (dH0dna + dH0dnb + rhoval*d2H0dnadnb
     &         +  dH1dna + dH1dnb + rhoval*d2H1dnadnb)*wght
            if (ipol.eq.2)
     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &         + (2.0d0*dH0dnb + rhoval*d2H0dnb2
     &         +  2.0d0*dH1dnb + rhoval*d2H1dnb2)*wght
#endif
c        endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        goto 20
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (dsqgamma.gt.TOLL)then
c
c           H0 part
c
            dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
            dfAtdg = dfAtdt*dtdg
            darglogdg = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
            dH0dg = 0.5d0*BETA**2/ALPHA*(phi**3)*darglogdg/arglog
c
c           H1 part
c
            dH1argexpdg = 2.0d0*H1argexp/t*dtdg
            dH1prefacdg = 2.0d0*H1prefac/t*dtdg
            dH1dg = dH1prefacdg * expinH1
     &            + H1prefac * dH1argexpdg * expinH1
c
c           Now form Cmat
c
            if (ipol.eq.1) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
     &                                         + dH1dg*rhoval*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
     &                                         + dH1dg*rhoval*wght*2.0d0
            else
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
     &                                         + dH1dg*rhoval*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
     &                                         + dH1dg*rhoval*wght*2.0d0
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*wght
     &                                         + dH1dg*rhoval*wght
            endif
#ifdef SECOND_DERIV
c
c           H0 part
c
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
            d2arglogdnadg = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdg*fAt
     &                                 +2.0d0*t*d2tdnadg*fAt
     &                                 +2.0d0*t*dtdg*dfAtdna
     &                                 +2.0d0*t*dtdna*dfAtdg
     &                                 +t*t*d2fAtdnadg)
            d2arglogdnbdg = 2.0d0*ALPHA/BETA*(2.0d0*dtdnb*dtdg*fAt
     &                                 +2.0d0*t*d2tdnbdg*fAt
     &                                 +2.0d0*t*dtdg*dfAtdnb
     &                                 +2.0d0*t*dtdnb*dfAtdg
     &                                 +t*t*d2fAtdnbdg)
            d2arglogdg2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdg*dtdg*fAt
     &                               +2.0d0*t*d2tdg2*fAt
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +t*t*d2fAtdg2)
            d2H0dnadg = 0.5d0*BETA**2/ALPHA*3.0d0*phi**2
     &                  *dphidna*darglogdg/arglog
     &                + 0.5d0*BETA**2/ALPHA*phi**3
     &                  *d2arglogdnadg/arglog
     &                - 0.5d0*BETA**2/ALPHA*phi**3
     &                  *darglogdg*darglogdna/arglog**2
            d2H0dnbdg = 0.5d0*BETA**2/ALPHA*3.0d0
     &                  *phi**2*dphidnb*darglogdg/arglog 
     &                + 0.5d0*BETA**2/ALPHA*phi**3
     &                  *d2arglogdnbdg/arglog
     &                - 0.5d0*BETA**2/ALPHA*phi**3
     &                  *darglogdg*darglogdnb/arglog**2
            d2H0dg2 = 0.5d0*BETA**2/ALPHA*phi**3
     &                *d2arglogdg2/arglog
     &              - 0.5d0*BETA**2/ALPHA*phi**3
     &                *darglogdg*darglogdg/arglog**2
c
c           H1 part
c
            d2H1argexpdnadg = 2.0d0*dH1argexpdna/t*dtdg
     &                      - 2.0d0*H1argexp/t**2*dtdg*dtdna
     &                      + 2.0d0*H1argexp/t*d2tdnadg
            d2H1argexpdnbdg = 2.0d0*dH1argexpdnb/t*dtdg
     &                      - 2.0d0*H1argexp/t**2*dtdg*dtdnb
     &                      + 2.0d0*H1argexp/t*d2tdnbdg
            d2H1argexpdg2 = 2.0d0*dH1argexpdg/t*dtdg
     &                    - 2.0d0*H1argexp/t**2*dtdg*dtdg
     &                    + 2.0d0*H1argexp/t*d2tdg2
            d2H1prefacdnadg = 2.0d0*dH1prefacdna/t*dtdg
     &                  - 2.0d0*H1prefac/t**2*dtdna*dtdg
     &                  + 2.0d0*H1prefac/t*d2tdnadg
            d2H1prefacdnbdg = 2.0d0*dH1prefacdnb/t*dtdg
     &                  - 2.0d0*H1prefac/t**2*dtdnb*dtdg
     &                  + 2.0d0*H1prefac/t*d2tdnbdg
            d2H1prefacdg2 = 2.0d0*dH1prefacdg/t*dtdg
     &                  - 2.0d0*H1prefac/t**2*dtdg*dtdg
     &                  + 2.0d0*H1prefac/t*d2tdg2
            d2H1dnadg = d2H1prefacdnadg * expinH1
     &            + dH1prefacdna * dH1argexpdg * expinH1
     &            + dH1prefacdg * dH1argexpdna * expinH1
     &            + H1prefac * d2H1argexpdnadg * expinH1
     &            + H1prefac * dH1argexpdna * dH1argexpdg * expinH1
            d2H1dnbdg = d2H1prefacdnbdg * expinH1
     &            + dH1prefacdnb * dH1argexpdg * expinH1
     &            + dH1prefacdg * dH1argexpdnb * expinH1
     &            + H1prefac * d2H1argexpdnbdg * expinH1
     &            + H1prefac * dH1argexpdnb * dH1argexpdg * expinH1
            d2H1dg2 = d2H1prefacdg2 * expinH1
     &            + dH1prefacdg * dH1argexpdg * expinH1
     &            + dH1prefacdg * dH1argexpdg * expinH1
     &            + H1prefac * d2H1argexpdg2 * expinH1
     &            + H1prefac * dH1argexpdg * dH1argexpdg * expinH1
c
c           Now form Cmat2
c
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &             + (dH0dg + d2H0dnadg*rhoval)*wght
     &             + (dH1dg + d2H1dnadg*rhoval)*wght
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
     &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*wght
     &             + 2.0d0*(dH1dg + d2H1dnadg*rhoval)*wght
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
     &             + (dH0dg + d2H0dnadg*rhoval)*wght
     &             + (dH1dg + d2H1dnadg*rhoval)*wght
            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
     &             + (dH0dg + d2H0dnbdg*rhoval)*wght
     &             + (dH1dg + d2H1dnbdg*rhoval)*wght
            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
     &             + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*wght
     &             + 2.0d0*(dH1dg + d2H1dnbdg*rhoval)*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &             + (dH0dg + d2H0dnbdg*rhoval)*wght
     &             + (dH1dg + d2H1dnbdg*rhoval)*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &             + d2H0dg2*rhoval*wght
     &             + d2H1dg2*rhoval*wght
            Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
     &             + 2.0d0*d2H0dg2*rhoval*wght
     &             + 2.0d0*d2H1dg2*rhoval*wght
            Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
     &             + d2H0dg2*rhoval*wght
     &             + d2H1dg2*rhoval*wght
            Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
     &             + 4.0d0*d2H0dg2*rhoval*wght
     &             + 4.0d0*d2H1dg2*rhoval*wght
            Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
     &             + 2.0d0*d2H0dg2*rhoval*wght
     &             + 2.0d0*d2H1dg2*rhoval*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &             + d2H0dg2*rhoval*wght
     &             + d2H1dg2*rhoval*wght
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
#include "nwpwxc_c_p91.F"
#endif
C>
C> @}
