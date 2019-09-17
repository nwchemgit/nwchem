#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_p91.F
C> The Perdew 91 correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_p91_p(tol_rho,ipol,nq,wght,rho,rgamma,ffunc)
#else
      Subroutine nwxc_c_p91(tol_rho, ipol, nq, wght, rho, rgamma, ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_p91_d2(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc)
#else
      Subroutine nwxc_c_p91_d3(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, ffunc)
#endif
c
c$Id$
c
#include "nwad.fh"
c
      Implicit none
#include "intf_nwxc_c_pw91lda.fh"
c
#include "nwxc_param.fh"
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
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
c                                   !< and possibly rho
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
      type(nwad_dble)::rho_t(3)
c     type(nwad_dble)::dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n
      type(nwad_dble)::rhoval, gammaval
c     type(nwad_dble)::nepsc(1),epsc,kf,ks,zeta,t,rs,cnum,cden,ccrs
      type(nwad_dble)::nepsc(1),epsc,kf,ks,zeta,t2,rs,cnum,cden,ccrs
      type(nwad_dble)::fAt,fAtnum,fAtden,argexp,arglog,H0,H1
      type(nwad_dble)::H1argexp,H1prefac,expinH1
      type(nwad_dble)::phi,expinA,A
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
      double precision drsdna, drsdnb
      double precision dccrsdna, dccrsdnb
      double precision dcnumdna, dcnumdnb
      double precision dcdendna, dcdendnb
      double precision dkfdna, dkfdnb
      double precision dH1argexpdna, dH1argexpdnb
      double precision dH1argexpdg
      double precision dH1prefacdna, dH1prefacdnb
      double precision dH1prefacdg
      double precision dH1dna, dH1dnb, dH1dg
#ifdef SECOND_DERIV
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
c     double precision d2rsdna2, d2rsdnb2, d2rsdnadnb
c     double precision d2ccrsdna2, d2ccrsdnb2, d2ccrsdnadnb
c     double precision d2cnumdna2, d2cnumdnb2, d2cnumdnadnb
c     double precision d2cdendna2, d2cdendnb2, d2cdendnadnb
c     double precision d2kfdna2, d2kfdnb2, d2kfdnadnb
c     double precision d2H1argexpdna2, d2H1argexpdnadnb, d2H1argexpdnb2
c     double precision d2H1argexpdnadg, d2H1argexpdnbdg, d2H1argexpdg2
c     double precision d2H1prefacdna2, d2H1prefacdnadnb, d2H1prefacdnb2
c     double precision d2H1prefacdnadg, d2H1prefacdnbdg, d2H1prefacdg2
c     double precision d2H1dna2, d2H1dnb2, d2H1dnadnb
c     double precision d2H1dnadg, d2H1dnbdg, d2H1dg2
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
         rho_t(2) = 0.0d0
         gammaval = rgamma(n,G_TT)
         if (ipol.eq.2) then
            rhoval   = 0.0d0
            rho_t(1) = 0.0d0
            rho_t(2) = 0.0d0
            gammaval = 0.0d0
            if (rho(n,R_A).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_A)
              rho_t(1) = rho(n,R_A)
              gammaval = gammaval + rgamma(n,G_AA)
            endif
            if (rho(n,R_B).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_B)
              rho_t(2) = rho(n,R_B)
              gammaval = gammaval + rgamma(n,G_BB)
              if (rho(n,R_A).gt.0.5d0*tol_rho) then
                gammaval = gammaval + 2.0d0*rgamma(n,G_AB)
              endif
            endif
         endif
         if (rhoval.le.tol_rho) goto 20
         if (ipol.eq.1) then
c           gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1)
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
         endif
c        dsqgamma = max(sqrt(gammaval),tol_rho)
c        dsqgamma = sqrt(gammaval)
         nepsc(1) = 0.0d0
c        dnepscdn(1) = 0.0d0
c        if (ipol.eq.2) dnepscdn(2) = 0.0d0
#ifdef SECOND_DERIV
c        d2nepscdn2(D2_RA_RA)=0.0d0
c        d2nepscdn2(D2_RA_RB)=0.0d0
c        if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
c
c        call for LDA bit
c
         lnq = 1
         lqwght = 1.0d0
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_pw91lda_p(tol_rho,ipol,lnq,lqwght,rho_t,nepsc)
#else
         call nwxc_c_pw91lda(tol_rho,ipol,lnq,lqwght,rho_t,nepsc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_pw91lda_d2(tol_rho,ipol,lnq,lqwght,rho_t,nepsc)
#else
         call nwxc_c_pw91lda_d3(tol_rho,ipol,lnq,lqwght,rho_t,nepsc)
#endif
c
c        ============
c        PW91 H0 part
c        ============
c        if(abs(nepsc(1)).lt.tol_rho*tol_rho) goto 20
c
c        epsilon_c = n*epsilon_c / n
c
         epsc = nepsc(1)/rhoval
c        if (ipol.eq.1) then
c           depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
c           depscdnb = depscdna
c        else
c           depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
c           depscdnb = dnepscdn(2)/rhoval-nepsc/(rhoval**2)
c        endif
#ifdef SECOND_DERIV
c        if (ipol.eq.1) then
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnb2   = d2epscdna2
c        else
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
c    &                     -dnepscdn(1)/(rhoval**2)
c    &                     -dnepscdn(2)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnb2   = d2nepscdn2(D2_RB_RB)/rhoval
c    &                     -dnepscdn(2)/(rhoval**2)
c    &                     -dnepscdn(2)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c        endif
#endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
         ks = (2.0d0/sqrt(PI))*((3.0d0*PI*PI)*rhoval)**(1.0d0/6.0d0)
c        dksdna = (1.0d0/6.0d0)*ks/rhoval
c        dksdnb = dksdna
#ifdef SECOND_DERIV
c        d2ksdna2   = (1.0d0/6.0d0)*dksdna/rhoval
c    &              - (1.0d0/6.0d0)*ks/(rhoval**2)
c        d2ksdnadnb = d2ksdna2
c        d2ksdnb2   = d2ksdna2
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
c           dzetadna = 1.0d0/rhoval
c           dzetadnb = -1.0d0/rhoval
#ifdef SECOND_DERIV
c           d2zetadna2   = -2.0d0/(rhoval**2)
c           d2zetadnadnb = 0.0d0
c           d2zetadnb2   = -2.0d0/(rhoval**2)
#endif
c        else
c           dzetadna =  2.0d0*rho(n,R_B)/(rhoval**2)
c           dzetadnb = -2.0d0*rho(n,R_A)/(rhoval**2)
#ifdef SECOND_DERIV
c           d2zetadna2   = -4.0d0*rho(n,R_B)/(rhoval**3)
c           d2zetadnadnb = 2.0d0*(rho(n,R_A)-rho(n,R_B))/(rhoval**3)
c           d2zetadnb2   = -4.0d0*rho(n,R_A)/(rhoval**3)
#endif
c        endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
         phi = 0.0d0
         if (1.0d0+zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0))
         endif
         if (1.0d0-zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*((1.0d0-zeta)**(2.0d0/3.0d0))
         endif
c        phi = 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0)
c    &               +(1.0d0-zeta)**(2.0d0/3.0d0))
c        if ((1.0d0-zeta).lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &             (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta))
c        else if ((1.0d0+zeta).lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &            -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
c        else
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &         (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta)
c    &        -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
c        endif
c        dphidna = dphidzeta*dzetadna
c        dphidnb = dphidzeta*dzetadnb
#ifdef SECOND_DERIV
c        if ((1.0d0-zeta).lt.tol_rho) then
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         (1.0d0+zeta)**(2.0d0/3.0d0)/((1.0d0+zeta)**2))
c        else if ((1.0d0+zeta).lt.tol_rho) then
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         (1.0d0-zeta)**(2.0d0/3.0d0)/((1.0d0-zeta)**2))
c        else
c           d2phidzeta2 = -(1.0d0/9.0d0)*(
c    &         (1.0d0+zeta)**(2.0d0/3.0d0)/((1.0d0+zeta)**2)
c    &        +(1.0d0-zeta)**(2.0d0/3.0d0)/((1.0d0-zeta)**2))
c        endif
c        d2phidna2   = d2phidzeta2*dzetadna*dzetadna
c    &               + dphidzeta*d2zetadna2 
c        d2phidnadnb = d2phidzeta2*dzetadna*dzetadnb
c    &               + dphidzeta*d2zetadnadnb
c        d2phidnb2   = d2phidzeta2*dzetadnb*dzetadnb
c    &               + dphidzeta*d2zetadnb2 
#endif
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
c        t = dsqgamma/(2.0d0*phi*ks*rhoval)
         t2 = gammaval/((2.0d0*phi*ks*rhoval)**2.0d0)
c        dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
c        dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
#ifdef SECOND_DERIV
c        d2tdna2 = - dtdna/rhoval
c    &           + t/(rhoval**2)
c    &           - dtdna/phi*dphidna
c    &           + t/(phi**2)*(dphidna**2)
c    &           - t/phi*d2phidna2
c    &           - dtdna/ks*dksdna
c    &           + t/(ks**2)*(dksdna**2)
c    &           - t/ks*d2ksdna2
c        d2tdnadnb = - dtdnb/rhoval
c    &           + t/(rhoval**2)
c    &           - dtdnb/phi*dphidna
c    &           + t/(phi**2)*(dphidna*dphidnb)
c    &           - t/phi*d2phidnadnb
c    &           - dtdnb/ks*dksdna
c    &           + t/(ks**2)*(dksdna*dksdnb)
c    &           - t/ks*d2ksdnadnb
c        d2tdnb2 = - dtdna/rhoval
c    &           + t/(rhoval**2)
c    &           - dtdnb/phi*dphidnb
c    &           + t/(phi**2)*(dphidnb**2)
c    &           - t/phi*d2phidnb2
c    &           - dtdnb/ks*dksdnb
c    &           + t/(ks**2)*(dksdnb**2)
c    &           - t/ks*d2ksdnb2
#endif
c
c        { ... } in A (see below) and its derivs
c
         argexp = (-2.0d0*ALPHA/BETA**2)*epsc/(phi**3.0d0)
c        dargexpdna = -2.0d0*ALPHA*depscdna/BETA**2/(phi**3)
c    &                +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*dphidna
c        dargexpdnb = -2.0d0*ALPHA*depscdnb/BETA**2/(phi**3)
c    &                +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*dphidnb
#ifdef SECOND_DERIV
c        d2argexpdna2 = -2.0d0*ALPHA*d2epscdna2/BETA**2/(phi**3)
c    &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidna
c    &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidna
c    &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidna**2
c    &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidna2
c        d2argexpdnadnb = -2.0d0*ALPHA*d2epscdnadnb/BETA**2/(phi**3)
c    &        +6.0d0*ALPHA*depscdna/BETA**2/(phi**4)*dphidnb
c    &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidna
c    &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidna*dphidnb
c    &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidnadnb
c        d2argexpdnb2 = -2.0d0*ALPHA*d2epscdnb2/BETA**2/(phi**3)
c    &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidnb
c    &        +6.0d0*ALPHA*depscdnb/BETA**2/(phi**4)*dphidnb
c    &        -24.0d0*ALPHA*epsc/BETA**2/(phi**5)*dphidnb**2
c    &        +6.0d0*ALPHA*epsc/BETA**2/(phi**4)*d2phidnb2
#endif
c
c        A = 2*ALPHA/BETA [exp{-2*ALPHA*epsilon_c/(BETA**2*phi**3)}-1]**(-1)
c
         if (abs(argexp).lt.EXPTOL) then
            expinA=exp(argexp)
         else
            expinA=0.0d0
         endif
         A = (2.0d0*ALPHA/BETA)/(expinA-1.0d0)
c        dAdna = -2.0d0*ALPHA/BETA*dargexpdna*expinA/(expinA-1.0d0)**2
c        dAdnb = -2.0d0*ALPHA/BETA*dargexpdnb*expinA/(expinA-1.0d0)**2
#ifdef SECOND_DERIV
c        d2Adna2   = -2.0d0*ALPHA/BETA*d2argexpdna2
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - 2.0d0*ALPHA/BETA*dargexpdna
c    &               *dargexpdna*expinA/(expinA-1.0d0)**2
c    &             + 4.0d0*ALPHA/BETA*dargexpdna*dargexpdna
c    &               *expinA*expinA/(expinA-1.0d0)**3
c        d2Adnadnb  = -2.0d0*ALPHA/BETA*d2argexpdnadnb
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - 2.0d0*ALPHA/BETA*dargexpdna
c    &               *dargexpdnb*expinA/(expinA-1.0d0)**2
c    &             + 4.0d0*ALPHA/BETA*dargexpdna*dargexpdnb
c    &               *expinA*expinA/(expinA-1.0d0)**3
c        d2Adnb2   = -2.0d0*ALPHA/BETA*d2argexpdnb2
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - 2.0d0*ALPHA/BETA*dargexpdnb
c    &               *dargexpdnb*expinA/(expinA-1.0d0)**2
c    &             + 4.0d0*ALPHA/BETA*dargexpdnb*dargexpdnb
c    &               *expinA*expinA/(expinA-1.0d0)**3
#endif
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
c        fAtnum = 1.0d0+A*(t**2.0d0)
c        fAtden = 1.0d0+A*(t**2.0d0)+(A**2.0d0)*(t**4.0d0)
         fAtnum = 1.0d0+A*t2
         fAtden = 1.0d0+A*t2+(A**2.0d0)*(t2**2.0d0)
         fAt = fAtnum/fAtden
c        dfAtnumdt = 2.0d0*A*t
c        dfAtnumdA = t**2
c        dfAtdendt = 2.0d0*A*t+4.0d0*A**2*t**3
c        dfAtdendA = t**2+2.0d0*A*t**4
c        dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/(fAtden**2)
c        dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/(fAtden**2)
c        dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
c        dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
#ifdef SECOND_DERIV
c        d2fAtnumdt2 = 2.0d0*A
c        d2fAtdendt2 = 2.0d0*A+12.0d0*A**2*t**2
c        d2fAtnumdtdA = 2.0d0*t
c        d2fAtdendtdA = 2.0d0*t+8.0d0*A*t**3
c        d2fAtnumdA2 = 0.0d0
c        d2fAtdendA2 = 2.0d0*t**4
c        d2fAtdt2  = (d2fAtnumdt2*fAtden-fAtnum*d2fAtdendt2)
c    &               /(fAtden**2)
c    &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
c    &               /(fAtden**3)*dfAtdendt
c        d2fAtdtdA = (d2fAtnumdtdA*fAtden+dfAtnumdt*dfAtdendA
c    &               -dfAtnumdA*dfAtdendt-fAtnum*d2fAtdendtdA)
c    &               /(fAtden**2)
c    &               -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
c    &               /(fAtden**3)*dfAtdendA
c        d2fAtdA2  = (d2fAtnumdA2*fAtden-fAtnum*d2fAtdendA2)
c    &               /(fAtden**2)
c    &               -2.0d0*(dfAtnumdA*fAtden-fAtnum*dfAtdendA)
c    &               /(fAtden**3)*dfAtdendA
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
c
c        arglog = 1 + 2*ALPHA/BETA * t**2 * fAt and its derivs
c
c        arglog = 1.0d0 + (2.0d0*ALPHA/BETA)*(t**2.0d0)*fAt
         arglog = 1.0d0 + (2.0d0*ALPHA/BETA)*t2*fAt
c        darglogdna = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdna*fAt
c    &                            +t*t*dfAtdna)
c        darglogdnb = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdnb*fAt
c    &                            +t*t*dfAtdnb)
#ifdef SECOND_DERIV
c        d2arglogdna2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdna*fAt
c    &                             +2.0d0*t*d2tdna2*fAt
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +t*t*d2fAtdna2)
c        d2arglogdnb2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdnb*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnb2*fAt
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +t*t*d2fAtdnb2)
c        d2arglogdnadnb = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnadnb*fAt
c    &                             +2.0d0*t*dtdna*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdna
c    &                             +t*t*d2fAtdnadnb)
#endif
c
c        H0 = BETA**2/2/ALPHA * phi**3 * log{arglog} and its derivs
c
         H0 = (0.5d0*BETA**2/ALPHA)*(phi**3.0d0)*log(arglog)
c        dH0dna = 0.5d0*BETA**2/ALPHA*(3.0d0*(phi**2)
c    &                  *dphidna*dlog(arglog)
c    &                  +(phi**3)*darglogdna/arglog)
c        dH0dnb = 0.5d0*BETA**2/ALPHA*(3.0d0*(phi**2)
c    &                  *dphidnb*dlog(arglog)
c    &                  +(phi**3)*darglogdnb/arglog)
#ifdef SECOND_DERIV
c        d2H0dna2 = 0.5d0*BETA**2/ALPHA
c    &                *(6.0d0*phi*dphidna*dphidna*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidna2*dlog(arglog)
c    &                +6.0d0*(phi**2)*dphidna*darglogdna/arglog
c    &                +(phi**3)*d2arglogdna2/arglog
c    &                -(phi**3)*darglogdna*darglogdna/arglog/arglog)
c        d2H0dnadnb = 0.5d0*BETA**2/ALPHA
c    &                *(6.0d0*phi*dphidna*dphidnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidnadnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*dphidna*darglogdnb/arglog
c    &                +3.0d0*(phi**2)*dphidnb*darglogdna/arglog
c    &                +(phi**3)*d2arglogdnadnb/arglog
c    &                -(phi**3)*darglogdna*darglogdnb/arglog/arglog)
c        d2H0dnb2 = 0.5d0*BETA**2/ALPHA
c    &                *(6.0d0*phi*dphidnb*dphidnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidnb2*dlog(arglog)
c    &                +6.0d0*(phi**2)*dphidnb*darglogdnb/arglog
c    &                +(phi**3)*d2arglogdnb2/arglog
c    &                -(phi**3)*darglogdnb*darglogdnb/arglog/arglog)
#endif
c        ============
c        PW91 H1 part
c        ============
c
c        rs = (3/4 / PI / n)**(1/3) and its derivs
c        
         rs = ((0.75d0 / PI) / rhoval)**(1.0d0/3.0d0)
c        drsdna = (-1.0d0/3.0d0)*rs/rhoval
c        drsdnb = drsdna
#ifdef SECOND_DERIV
c        d2rsdna2 = (-4.0d0/3.0d0)*drsdna/rhoval
c        d2rsdnb2 = d2rsdna2
c        d2rsdnadnb = d2rsdna2
#endif
c
c        CC(rs) and its derivs
c
         cnum = cnoname + ca*rs + cb*(rs**2.0d0)
         cden = 1.0d0 + cc*rs + cd*(rs**2.0d0) + (10.0d0*cb)*(rs**3.0d0)
         ccrs = 1.0d-3*cnum/cden - CX
c        dcnumdna = (ca + 2.0d0*cb*rs)*drsdna
c        dcnumdnb = (ca + 2.0d0*cb*rs)*drsdnb
c        dcdendna = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*drsdna
c        dcdendnb = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*drsdnb
c        dccrsdna = 1.0d-3*(dcnumdna*cden-cnum*dcdendna)/cden**2
c        dccrsdnb = 1.0d-3*(dcnumdnb*cden-cnum*dcdendnb)/cden**2
#ifdef SECOND_DERIV
c        d2cnumdna2 = (ca + 2.0d0*cb*rs)*d2rsdna2
c    &              + (2.0d0*cb)*drsdna*drsdna
c        d2cnumdnb2 = (ca + 2.0d0*cb*rs)*d2rsdnb2
c    &              + (2.0d0*cb)*drsdnb*drsdnb
c        d2cnumdnadnb = (ca + 2.0d0*cb*rs)*d2rsdnadnb
c    &              + (2.0d0*cb)*drsdna*drsdnb
c        d2cdendna2 = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdna2
c    &              + (2.0d0*cd+60.0d0*cb*rs)*drsdna*drsdna
c        d2cdendnb2 = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdnb2
c    &              + (2.0d0*cd+60.0d0*cb*rs)*drsdnb*drsdnb
c        d2cdendnadnb = (cc + 2.0d0*cd*rs + 30.0d0*cb*rs**2)*d2rsdnadnb
c    &              + (2.0d0*cd+60.0d0*cb*rs)*drsdna*drsdnb
c        d2ccrsdna2 = 1.0d-3*(d2cnumdna2*cden-cnum*d2cdendna2)/cden**2
c    &      - 1.0d-3*2.0d0*(dcnumdna*cden-cnum*dcdendna)
c    &        /cden**3*dcdendna
c        d2ccrsdnb2 = 1.0d-3*(d2cnumdnb2*cden-cnum*d2cdendnb2)/cden**2
c    &      - 1.0d-3*2.0d0*(dcnumdnb*cden-cnum*dcdendnb)
c    &        /cden**3*dcdendnb
c        d2ccrsdnadnb 
c    &      = 1.0d-3*(d2cnumdnadnb*cden-cnum*d2cdendnadnb)/cden**2
c    &      + 1.0d-3*(dcnumdna*dcdendnb-dcnumdnb*dcdendna)/cden**2
c    &      - 1.0d-3*2.0d0*(dcnumdna*cden-cnum*dcdendna)
c    &        /cden**3*dcdendnb
#endif
c
c        kf = (3 * PI**2 * n)**(1/3)
c
         kf = ((3.0d0*PI**2)*rhoval)**(1.0d0/3.0d0)
c        dkfdna = (1.0d0/3.0d0)*kf/rhoval
c        dkfdnb = dkfdna
#ifdef SECOND_DERIV
c        d2kfdna2 = (-2.0d0/3.0d0)*dkfdna/rhoval
c        d2kfdnadnb = (-2.0d0/3.0d0)*dkfdna/rhoval
c        d2kfdnb2 = (-2.0d0/3.0d0)*dkfdnb/rhoval
#endif
c
c        H1argexp = -100 * phi**4 * (ks**2/kf**2) * t**2 and its derivs
c
c        H1argexp = -100.0d0*(phi**4.0d0)*((ks**2.0d0)
c    &            / (kf**2.0d0))*(t**2.0d0)
         H1argexp = -100.0d0*(phi**4.0d0)*((ks**2.0d0)
     &            / (kf**2.0d0))*t2
c        dH1argexpdna = 4.0d0*H1argexp/phi*dphidna
c    &                 +2.0d0*H1argexp/ks*dksdna
c    &                 -2.0d0*H1argexp/kf*dkfdna
c    &                 +2.0d0*H1argexp/t*dtdna
c        dH1argexpdnb = 4.0d0*H1argexp/phi*dphidnb
c    &                 +2.0d0*H1argexp/ks*dksdnb
c    &                 -2.0d0*H1argexp/kf*dkfdnb
c    &                 +2.0d0*H1argexp/t*dtdnb
#ifdef SECOND_DERIV
c        d2H1argexpdna2 = 4.0d0*dH1argexpdna/phi*dphidna
c    &                   -4.0d0*H1argexp/phi**2*dphidna**2
c    &                   +4.0d0*H1argexp/phi*d2phidna2
c    &                   +2.0d0*dH1argexpdna/ks*dksdna
c    &                   -2.0d0*H1argexp/ks**2*dksdna**2
c    &                   +2.0d0*H1argexp/ks*d2ksdna2
c    &                   -2.0d0*dH1argexpdna/kf*dkfdna
c    &                   +2.0d0*H1argexp/kf**2*dkfdna**2
c    &                   -2.0d0*H1argexp/kf*d2kfdna2
c    &                   +2.0d0*dH1argexpdna/t*dtdna
c    &                   -2.0d0*H1argexp/t**2*dtdna**2
c    &                   +2.0d0*H1argexp/t*d2tdna2
c        d2H1argexpdnb2 = 4.0d0*dH1argexpdnb/phi*dphidnb
c    &                   -4.0d0*H1argexp/phi**2*dphidnb**2
c    &                   +4.0d0*H1argexp/phi*d2phidnb2
c    &                   +2.0d0*dH1argexpdnb/ks*dksdnb
c    &                   -2.0d0*H1argexp/ks**2*dksdnb**2
c    &                   +2.0d0*H1argexp/ks*d2ksdnb2
c    &                   -2.0d0*dH1argexpdnb/kf*dkfdnb
c    &                   +2.0d0*H1argexp/kf**2*dkfdnb**2
c    &                   -2.0d0*H1argexp/kf*d2kfdnb2
c    &                   +2.0d0*dH1argexpdnb/t*dtdnb
c    &                   -2.0d0*H1argexp/t**2*dtdnb**2
c    &                   +2.0d0*H1argexp/t*d2tdnb2
c        d2H1argexpdnadnb = 4.0d0*dH1argexpdna/phi*dphidnb
c    &                   -4.0d0*H1argexp/phi**2*dphidna*dphidnb
c    &                   +4.0d0*H1argexp/phi*d2phidnadnb
c    &                   +2.0d0*dH1argexpdna/ks*dksdnb
c    &                   -2.0d0*H1argexp/ks**2*dksdna*dksdnb
c    &                   +2.0d0*H1argexp/ks*d2ksdnadnb
c    &                   -2.0d0*dH1argexpdna/kf*dkfdnb
c    &                   +2.0d0*H1argexp/kf**2*dkfdna*dkfdnb
c    &                   -2.0d0*H1argexp/kf*d2kfdnadnb
c    &                   +2.0d0*dH1argexpdna/t*dtdnb
c    &                   -2.0d0*H1argexp/t**2*dtdna*dtdnb
c    &                   +2.0d0*H1argexp/t*d2tdnadnb
#endif
c
c        H1prefac = NU*[CC(rs) - CC(0)-3/7*CX] * g**3 * t**2
c
c        H1prefac = NU*(ccrs - (CCZERO + (3.0d0/7.0d0)*CX))
c    &            * (phi**3.0d0) * (t**2.0d0)
         H1prefac = NU*(ccrs - (CCZERO + (3.0d0/7.0d0)*CX))
     &            * (phi**3.0d0) * t2
c        dH1prefacdna = NU*dccrsdna * phi**3 * t**2
c    &                + 3.0d0*H1prefac/phi*dphidna
c    &                + 2.0d0*H1prefac/t*dtdna
c        dH1prefacdnb = NU*dccrsdnb * phi**3 * t**2
c    &                + 3.0d0*H1prefac/phi*dphidnb
c    &                + 2.0d0*H1prefac/t*dtdnb
#ifdef SECOND_DERIV
c        d2H1prefacdna2 = NU*d2ccrsdna2 * phi**3 * t**2
c    &                  + 3.0d0* NU*dccrsdna * phi**2*dphidna * t**2
c    &                  + 2.0d0* NU*dccrsdna * phi**3 * t*dtdna
c    &                  + 3.0d0*dH1prefacdna/phi*dphidna
c    &                  - 3.0d0*H1prefac/phi**2*dphidna**2
c    &                  + 3.0d0*H1prefac/phi*d2phidna2
c    &                  + 2.0d0*dH1prefacdna/t*dtdna
c    &                  - 2.0d0*H1prefac/t**2*dtdna**2
c    &                  + 2.0d0*H1prefac/t*d2tdna2
c        d2H1prefacdnb2 = NU*d2ccrsdnb2 * phi**3 * t**2
c    &                  + 3.0d0* NU*dccrsdnb * phi**2*dphidnb * t**2
c    &                  + 2.0d0* NU*dccrsdnb * phi**3 * t*dtdnb
c    &                  + 3.0d0*dH1prefacdnb/phi*dphidnb
c    &                  - 3.0d0*H1prefac/phi**2*dphidnb**2
c    &                  + 3.0d0*H1prefac/phi*d2phidnb2
c    &                  + 2.0d0*dH1prefacdnb/t*dtdnb
c    &                  - 2.0d0*H1prefac/t**2*dtdnb**2
c    &                  + 2.0d0*H1prefac/t*d2tdnb2
c        d2H1prefacdnadnb = NU*d2ccrsdnadnb * phi**3 * t**2
c    &                  + 3.0d0* NU*dccrsdna * phi**2*dphidnb * t**2
c    &                  + 2.0d0* NU*dccrsdna * phi**3 * t*dtdnb
c    &                  + 3.0d0*dH1prefacdna/phi*dphidnb
c    &                  - 3.0d0*H1prefac/phi**2*dphidna*dphidnb
c    &                  + 3.0d0*H1prefac/phi*d2phidnadnb
c    &                  + 2.0d0*dH1prefacdna/t*dtdnb
c    &                  - 2.0d0*H1prefac/t**2*dtdna*dtdnb
c    &                  + 2.0d0*H1prefac/t*d2tdnadnb
#endif
c
c        H1 = H1prefac * exp(H1argexp)
c 
         if (abs(H1argexp).lt.EXPTOL) then
            expinH1 = exp(H1argexp)
         else
            expinH1 = 0.0d0
         endif
         H1 = H1prefac * expinH1
c        dH1dna = dH1prefacdna * expinH1
c    &          + H1prefac * dH1argexpdna * expinH1
c        dH1dnb = dH1prefacdnb * expinH1
c    &          + H1prefac * dH1argexpdnb * expinH1
#ifdef SECOND_DERIV
c        d2H1dna2 = d2H1prefacdna2 * expinH1
c    &            + dH1prefacdna * dH1argexpdna * expinH1
c    &            + dH1prefacdna * dH1argexpdna * expinH1
c    &            + H1prefac * d2H1argexpdna2 * expinH1
c    &            + H1prefac * dH1argexpdna * dH1argexpdna * expinH1
c        d2H1dnb2 = d2H1prefacdnb2 * expinH1
c    &            + dH1prefacdnb * dH1argexpdnb * expinH1
c    &            + dH1prefacdnb * dH1argexpdnb * expinH1
c    &            + H1prefac * d2H1argexpdnb2 * expinH1
c    &            + H1prefac * dH1argexpdnb * dH1argexpdnb * expinH1
c        d2H1dnadnb = d2H1prefacdnadnb * expinH1
c    &            + dH1prefacdna * dH1argexpdnb * expinH1
c    &            + dH1prefacdnb * dH1argexpdna * expinH1
c    &            + H1prefac * d2H1argexpdnadnb * expinH1
c    &            + H1prefac * dH1argexpdna * dH1argexpdnb * expinH1
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
c           Amat(n,D1_RA) = Amat(n,D1_RA) + (H0 + H1 +
c    &         rhoval*dH0dna + rhoval*dH1dna)*wght
c           if (ipol.eq.2) Amat(n,D1_RB) = Amat(n,D1_RB) + (H0 + H1 +
c    &         rhoval*dH0dnb + rhoval*dH1dnb)*wght
#ifdef SECOND_DERIV
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &         + (2.0d0*dH0dna + rhoval*d2H0dna2
c    &         +  2.0d0*dH1dna + rhoval*d2H1dna2)*wght
c           Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &         + (dH0dna + dH0dnb + rhoval*d2H0dnadnb
c    &         +  dH1dna + dH1dnb + rhoval*d2H1dnadnb)*wght
c           if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &         + (2.0d0*dH0dnb + rhoval*d2H0dnb2
c    &         +  2.0d0*dH1dnb + rhoval*d2H1dnb2)*wght
#endif
c        endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        goto 20
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        if (dsqgamma.gt.TOLL)then
c
c           H0 part
c
c           dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
c           dfAtdg = dfAtdt*dtdg
c           darglogdg = 2.0d0*ALPHA/BETA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
c           dH0dg = 0.5d0*BETA**2/ALPHA*(phi**3)*darglogdg/arglog
c
c           H1 part
c
c           dH1argexpdg = 2.0d0*H1argexp/t*dtdg
c           dH1prefacdg = 2.0d0*H1prefac/t*dtdg
c           dH1dg = dH1prefacdg * expinH1
c    &            + H1prefac * dH1argexpdg * expinH1
c
c           Now form Cmat
c
c           if (ipol.eq.1) then
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c    &                                         + dH1dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c    &                                         + dH1dg*rhoval*wght*2.0d0
c           else
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c    &                                         + dH1dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c    &                                         + dH1dg*rhoval*wght*2.0d0
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*wght
c    &                                         + dH1dg*rhoval*wght
c           endif
#ifdef SECOND_DERIV
c
c           H0 part
c
c           d2tdg2 = -0.125d0/(phi*ks*rhoval)/(dsqgamma**3)
c           d2tdnadg = -dtdg/rhoval-dtdg/phi*dphidna
c    &                 -dtdg/ks*dksdna
c           d2tdnbdg = -dtdg/rhoval-dtdg/phi*dphidnb
c    &                 -dtdg/ks*dksdnb
c           d2fAtdg2 = d2fAtdt2*(dtdg**2)+dfAtdt*d2tdg2
c           d2fAtdnadg = d2fAtdt2*dtdg*dtdna
c    &                  +d2fAtdtdA*dtdg*dAdna
c    &                  +dfAtdt*d2tdnadg 
c           d2fAtdnbdg = d2fAtdt2*dtdg*dtdnb
c    &                  +d2fAtdtdA*dtdg*dAdnb
c    &                  +dfAtdt*d2tdnbdg 
c           d2arglogdnadg = 2.0d0*ALPHA/BETA*(2.0d0*dtdna*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnadg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdna
c    &                                 +2.0d0*t*dtdna*dfAtdg
c    &                                 +t*t*d2fAtdnadg)
c           d2arglogdnbdg = 2.0d0*ALPHA/BETA*(2.0d0*dtdnb*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnbdg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdnb
c    &                                 +2.0d0*t*dtdnb*dfAtdg
c    &                                 +t*t*d2fAtdnbdg)
c           d2arglogdg2 = 2.0d0*ALPHA/BETA*(2.0d0*dtdg*dtdg*fAt
c    &                               +2.0d0*t*d2tdg2*fAt
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +t*t*d2fAtdg2)
c           d2H0dnadg = 0.5d0*BETA**2/ALPHA*3.0d0*phi**2
c    &                  *dphidna*darglogdg/arglog
c    &                + 0.5d0*BETA**2/ALPHA*phi**3
c    &                  *d2arglogdnadg/arglog
c    &                - 0.5d0*BETA**2/ALPHA*phi**3
c    &                  *darglogdg*darglogdna/arglog**2
c           d2H0dnbdg = 0.5d0*BETA**2/ALPHA*3.0d0
c    &                  *phi**2*dphidnb*darglogdg/arglog 
c    &                + 0.5d0*BETA**2/ALPHA*phi**3
c    &                  *d2arglogdnbdg/arglog
c    &                - 0.5d0*BETA**2/ALPHA*phi**3
c    &                  *darglogdg*darglogdnb/arglog**2
c           d2H0dg2 = 0.5d0*BETA**2/ALPHA*phi**3
c    &                *d2arglogdg2/arglog
c    &              - 0.5d0*BETA**2/ALPHA*phi**3
c    &                *darglogdg*darglogdg/arglog**2
c
c           H1 part
c
c           d2H1argexpdnadg = 2.0d0*dH1argexpdna/t*dtdg
c    &                      - 2.0d0*H1argexp/t**2*dtdg*dtdna
c    &                      + 2.0d0*H1argexp/t*d2tdnadg
c           d2H1argexpdnbdg = 2.0d0*dH1argexpdnb/t*dtdg
c    &                      - 2.0d0*H1argexp/t**2*dtdg*dtdnb
c    &                      + 2.0d0*H1argexp/t*d2tdnbdg
c           d2H1argexpdg2 = 2.0d0*dH1argexpdg/t*dtdg
c    &                    - 2.0d0*H1argexp/t**2*dtdg*dtdg
c    &                    + 2.0d0*H1argexp/t*d2tdg2
c           d2H1prefacdnadg = 2.0d0*dH1prefacdna/t*dtdg
c    &                  - 2.0d0*H1prefac/t**2*dtdna*dtdg
c    &                  + 2.0d0*H1prefac/t*d2tdnadg
c           d2H1prefacdnbdg = 2.0d0*dH1prefacdnb/t*dtdg
c    &                  - 2.0d0*H1prefac/t**2*dtdnb*dtdg
c    &                  + 2.0d0*H1prefac/t*d2tdnbdg
c           d2H1prefacdg2 = 2.0d0*dH1prefacdg/t*dtdg
c    &                  - 2.0d0*H1prefac/t**2*dtdg*dtdg
c    &                  + 2.0d0*H1prefac/t*d2tdg2
c           d2H1dnadg = d2H1prefacdnadg * expinH1
c    &            + dH1prefacdna * dH1argexpdg * expinH1
c    &            + dH1prefacdg * dH1argexpdna * expinH1
c    &            + H1prefac * d2H1argexpdnadg * expinH1
c    &            + H1prefac * dH1argexpdna * dH1argexpdg * expinH1
c           d2H1dnbdg = d2H1prefacdnbdg * expinH1
c    &            + dH1prefacdnb * dH1argexpdg * expinH1
c    &            + dH1prefacdg * dH1argexpdnb * expinH1
c    &            + H1prefac * d2H1argexpdnbdg * expinH1
c    &            + H1prefac * dH1argexpdnb * dH1argexpdg * expinH1
c           d2H1dg2 = d2H1prefacdg2 * expinH1
c    &            + dH1prefacdg * dH1argexpdg * expinH1
c    &            + dH1prefacdg * dH1argexpdg * expinH1
c    &            + H1prefac * d2H1argexpdg2 * expinH1
c    &            + H1prefac * dH1argexpdg * dH1argexpdg * expinH1
c
c           Now form Cmat2
c
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c    &             + (dH1dg + d2H1dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
c    &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*wght
c    &             + 2.0d0*(dH1dg + d2H1dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c    &             + (dH1dg + d2H1dnadg*rhoval)*wght
c           Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
c    &             + (dH0dg + d2H0dnbdg*rhoval)*wght
c    &             + (dH1dg + d2H1dnbdg*rhoval)*wght
c           Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
c    &             + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*wght
c    &             + 2.0d0*(dH1dg + d2H1dnbdg*rhoval)*wght
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &             + (dH0dg + d2H0dnbdg*rhoval)*wght
c    &             + (dH1dg + d2H1dnbdg*rhoval)*wght
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &             + d2H0dg2*rhoval*wght
c    &             + d2H1dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
c    &             + 2.0d0*d2H0dg2*rhoval*wght
c    &             + 2.0d0*d2H1dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
c    &             + d2H0dg2*rhoval*wght
c    &             + d2H1dg2*rhoval*wght
c           Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
c    &             + 4.0d0*d2H0dg2*rhoval*wght
c    &             + 4.0d0*d2H1dg2*rhoval*wght
c           Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
c    &             + 2.0d0*d2H0dg2*rhoval*wght
c    &             + 2.0d0*d2H1dg2*rhoval*wght
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &             + d2H0dg2*rhoval*wght
c    &             + d2H1dg2*rhoval*wght
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
c     Compile source again for Maxima
c
#include "nwxc_c_p91.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_p91.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_p91.F"
#endif
#undef NWAD_PRINT
C>
C> @}
