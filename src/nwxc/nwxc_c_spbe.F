#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_spbe.F
C> The simplified PBE correlation functional for SSB-D
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the SSB-D PBE like correlation functional
C>
C> Evaluates the SSB-D PBE like correlation functional [1].
C>
C> ### References ###
C>
C> [1] M. Swart, M. Sola, F.M. Bickelhaupt,
C>     "A new all-round density functional based on spin states and
C>     SN2 barriers", J. Chem. Phys. <b>131</b>, 094103 (2009), DOI:
C>     <a href="https://doi.org/10.1063/1.3213193">
C>     10.1063/1.3213193</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_spbe_p(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                         ffunc)
#else
      Subroutine nwxc_c_spbe(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                       ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_spbe_d2(tol_rho, ipol, nq, wght, rho, rgamma,
     &                          ffunc)
#else
      Subroutine nwxc_c_spbe_d3(tol_rho, ipol, nq, wght, rho, rgamma,
     &                          ffunc)
#endif
c
c$Id$
c
#include "nwad.fh"
c
      Implicit none
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
      double precision GAMMA, BETA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETA = 0.06672455060314922d0)
c     parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
      type(nwad_dble)::rho_t(2)
c     type(nwad_dble)::dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n
      type(nwad_dble)::rhoval, gammaval
      type(nwad_dble)::nepsc(1),epsc,H0,phi,zeta,arglog,
     +                 fAt,fAtnum,fAtden
      type(nwad_dble)::A,t2,ks,argexp,expinA
      double precision dnepscdn(2)
      double precision depscdna, depscdnb
      double precision dH0dna, dH0dnb, dH0dg
      double precision dphidna, dphidnb, dphidzeta
      double precision dzetadna, dzetadnb
      double precision darglogdna, darglogdnb, darglogdg
      double precision dfAtdt, dfAtdA
      double precision dfAtdendt, dfAtdendA
      double precision dfAtdna, dfAtdnb, dfAtdg
      double precision dAdna, dAdnb
      double precision dtdna, dtdnb, dtdg
      double precision dksdna, dksdnb
      double precision dargexpdna, dargexpdnb
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
c     double precision d2fAtdendt2, d2fAtdendtdA, d2fAtdendA2
c     double precision d2fAtdna2, d2fAtdnb2, d2fAtdnadnb
c     double precision d2fAtdnadg, d2fAtdnbdg
c     double precision d2Adna2, d2Adnadnb, d2Adnb2
c     double precision d2tdna2, d2tdnb2, d2tdnadnb
c     double precision d2tdg2, d2tdnadg, d2tdnbdg
c     double precision d2ksdna2, d2ksdnb2, d2ksdnadnb
c     double precision d2argexpdna2, d2argexpdnb2, d2argexpdnadnb
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
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      pi = acos(-1.0d0)
c
      do 20 n = 1, nq
c
c        n and zeta = (na - nb)/n
c
         if (ipol.eq.2) then
           rhoval = 0.0d0
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
         else
           rhoval = rho(n,R_T)
           rho_t(1) = rhoval
           gammaval = rgamma(n,G_TT)
         endif
         if (rhoval.le.tol_rho) goto 20
c        dsqgamma = sqrt(gammaval)
         nepsc(1) = 0.0d0
c        dnepscdn(D1_RA) = 0.0d0
c        if (ipol.eq.2) dnepscdn(D1_RB) = 0.0d0
#ifdef SECOND_DERIV
c        d2nepscdn2(D2_RA_RA)=0.0d0
c        d2nepscdn2(D2_RA_RB)=0.0d0
c        if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
c
c        ========
c        LDA part
c        ========
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_pw91lda_p(tol_rho, ipol, 1, 1.0d0, rho_t,
     &                         nepsc)
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
         epsc = nepsc(1)/rhoval
c        if (ipol.eq.1) then
c           depscdna = dnepscdn(D1_RA)/rhoval-nepsc/(rhoval**2)
c           depscdnb = depscdna
c        else
c           depscdna = dnepscdn(D1_RA)/rhoval-nepsc/(rhoval**2)
c           depscdnb = dnepscdn(D1_RB)/rhoval-nepsc/(rhoval**2)
c        endif
#ifdef SECOND_DERIV
c        if (ipol.eq.1) then
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnb2   = d2epscdna2
c        else
c           d2epscdna2   = d2nepscdn2(D2_RA_RA)/rhoval
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnadnb = d2nepscdn2(D2_RA_RB)/rhoval
c    &                     -dnepscdn(D1_RA)/(rhoval**2)
c    &                     -dnepscdn(D1_RB)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c           d2epscdnb2   = d2nepscdn2(D2_RB_RB)/rhoval
c    &                     -dnepscdn(D1_RB)/(rhoval**2)
c    &                     -dnepscdn(D1_RB)/(rhoval**2)
c    &                     +2.0d0*nepsc/(rhoval**3)
c        endif
#endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
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
c           d2zetadnb2   = 2.0d0/(rhoval**2)
#endif
c        else
c           dzetadna =  2.0d0*rho(n,R_B)/(rhoval**2)
c           dzetadnb = -2.0d0*rho(n,R_A)/(rhoval**2)
#ifdef SECOND_DERIV
c           d2zetadna2   = -4.0d0*rho(n,R_B)/(rhoval**3)
c           d2zetadnadnb = 2.0d0*(rho(n,R_A)-rho(n,R_B))/(rhoval**3)
c           d2zetadnb2   = 4.0d0*rho(n,R_A)/(rhoval**3)
#endif
c        endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
         phi = 0.0d0
         if (1.0d0+zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*(1.0d0+zeta)**(2.0d0/3.0d0)
         endif
         if (1.0d0-zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*(1.0d0-zeta)**(2.0d0/3.0d0)
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
c        d2tdnb2 = - dtdnb/rhoval
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
         argexp = -epsc/GAMMA/(phi**3.0d0)
c        dargexpdna = -depscdna/GAMMA/(phi**3)
c    &                +3.0d0*epsc/GAMMA/(phi**4)*dphidna
c        dargexpdnb = -depscdnb/GAMMA/(phi**3)
c    &                +3.0d0*epsc/GAMMA/(phi**4)*dphidnb
#ifdef SECOND_DERIV
c        d2argexpdna2 = -d2epscdna2/GAMMA/(phi**3)
c    &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidna
c    &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidna
c    &        -12.0d0*epsc/GAMMA/(phi**5)*dphidna**2
c    &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidna2
c        d2argexpdnadnb = -d2epscdnadnb/GAMMA/(phi**3)
c    &        +3.0d0*depscdna/GAMMA/(phi**4)*dphidnb
c    &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidna
c    &        -12.0d0*epsc/GAMMA/(phi**5)*dphidna*dphidnb
c    &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidnadnb
c        d2argexpdnb2 = -d2epscdnb2/GAMMA/(phi**3)
c    &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidnb
c    &        +3.0d0*depscdnb/GAMMA/(phi**4)*dphidnb
c    &        -12.0d0*epsc/GAMMA/(phi**5)*dphidnb**2
c    &        +3.0d0*epsc/GAMMA/(phi**4)*d2phidnb2
#endif
c
c        A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
         if (abs(argexp).lt.EXPTOL) then
            expinA=exp(argexp)
         else
            expinA=0.0d0
         endif
         A = BETA/GAMMA/(expinA-1.0d0)
c        dAdna = -BETA/GAMMA*dargexpdna*expinA/(expinA-1.0d0)**2
c        dAdnb = -BETA/GAMMA*dargexpdnb*expinA/(expinA-1.0d0)**2
#ifdef SECOND_DERIV
c        d2Adna2   = -BETA/GAMMA*d2argexpdna2
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - BETA/GAMMA*dargexpdna
c    &               *dargexpdna*expinA/(expinA-1.0d0)**2
c    &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdna
c    &               *expinA*expinA/(expinA-1.0d0)**3
c        d2Adnadnb  = -BETA/GAMMA*d2argexpdnadnb
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - BETA/GAMMA*dargexpdna
c    &               *dargexpdnb*expinA/(expinA-1.0d0)**2
c    &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdnb
c    &               *expinA*expinA/(expinA-1.0d0)**3
c        d2Adnb2   = -BETA/GAMMA*d2argexpdnb2
c    &               *expinA/(expinA-1.0d0)**2 
c    &             - BETA/GAMMA*dargexpdnb
c    &               *dargexpdnb*expinA/(expinA-1.0d0)**2
c    &             + 2.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb
c    &               *expinA*expinA/(expinA-1.0d0)**3
#endif
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
         fAtnum = 1.0d0
         fAtden = 1.0d0+A*t2
         fAt = fAtnum/fAtden
c        dfAtdendt = 2.0d0*A*t
c        dfAtdendA = t**2
c        dfAtdt = (-fAtnum*dfAtdendt)/(fAtden**2)
c        dfAtdA = (-fAtnum*dfAtdendA)/(fAtden**2)
c        dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
c        dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
#ifdef SECOND_DERIV
c        d2fAtdendt2 = 2.0d0*A
c        d2fAtdendtdA = 2.0d0*t
c        d2fAtdendA2 = 0.0d0
c        d2fAtdt2  = (-fAtnum*d2fAtdendt2)
c    &               /(fAtden**2)
c    &               -2.0d0*(-fAtnum*dfAtdendt)
c    &               /(fAtden**3)*dfAtdendt
c        d2fAtdtdA = (-fAtnum*d2fAtdendtdA)
c    &               /(fAtden**2)
c    &               -2.0d0*(-fAtnum*dfAtdendt)
c    &               /(fAtden**3)*dfAtdendA
c        d2fAtdA2  = (-fAtnum*d2fAtdendA2)
c    &               /(fAtden**2)
c    &               -2.0d0*(-fAtnum*dfAtdendA)
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
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + BETA/GAMMA*t2*fAt
c        darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
c    &                            +t*t*dfAtdna)
c        darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
c    &                            +t*t*dfAtdnb)
#ifdef SECOND_DERIV
c        d2arglogdna2 = BETA/GAMMA*(2.0d0*dtdna*dtdna*fAt
c    &                             +2.0d0*t*d2tdna2*fAt
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +2.0d0*t*dtdna*dfAtdna
c    &                             +t*t*d2fAtdna2)
c        d2arglogdnb2 = BETA/GAMMA*(2.0d0*dtdnb*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnb2*fAt
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdnb
c    &                             +t*t*d2fAtdnb2)
c        d2arglogdnadnb = BETA/GAMMA*(2.0d0*dtdna*dtdnb*fAt
c    &                             +2.0d0*t*d2tdnadnb*fAt
c    &                             +2.0d0*t*dtdna*dfAtdnb
c    &                             +2.0d0*t*dtdnb*dfAtdna
c    &                             +t*t*d2fAtdnadnb)
#endif
c
c        H0 = GAMMA * phi**3 * log{arglog} and its derivs
c
         H0 = GAMMA*(phi**3.0d0)*log(arglog)
c        dH0dna = GAMMA*(3.0d0*(phi**2)*dphidna*dlog(arglog)
c    &                  +(phi**3)*darglogdna/arglog)
c        dH0dnb = GAMMA*(3.0d0*(phi**2)*dphidnb*dlog(arglog)
c    &                  +(phi**3)*darglogdnb/arglog)
#ifdef SECOND_DERIV
c        d2H0dna2 = GAMMA*(6.0d0*phi*dphidna*dphidna*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidna2*dlog(arglog)
c    &                +6.0d0*(phi**2)*dphidna*darglogdna/arglog
c    &                +(phi**3)*d2arglogdna2/arglog
c    &                -(phi**3)*darglogdna*darglogdna/arglog/arglog)
c        d2H0dnadnb = GAMMA*(6.0d0*phi*dphidna*dphidnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidnadnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*dphidna*darglogdnb/arglog
c    &                +3.0d0*(phi**2)*dphidnb*darglogdna/arglog
c    &                +(phi**3)*d2arglogdnadnb/arglog
c    &                -(phi**3)*darglogdna*darglogdnb/arglog/arglog)
c        d2H0dnb2 = GAMMA*(6.0d0*phi*dphidnb*dphidnb*dlog(arglog)
c    &                +3.0d0*(phi**2)*d2phidnb2*dlog(arglog)
c    &                +6.0d0*(phi**2)*dphidnb*darglogdnb/arglog
c    &                +(phi**3)*d2arglogdnb2/arglog
c    &                -(phi**3)*darglogdnb*darglogdnb/arglog/arglog)
#endif
c
c        Now we update Ec, Amat, and Amat2
c
         ffunc(n)=ffunc(n)+nepsc(1)*wght               
         ffunc(n)=ffunc(n)+(H0*rhoval)*wght
c        Amat(n,D1_RA) = Amat(n,D1_RA) + dnepscdn(D1_RA)*wght
c        if (ipol.eq.2) Amat(n,D1_RB) = Amat(n,D1_RB) 
c    &                                + dnepscdn(D1_RB)*wght
#ifdef SECOND_DERIV
c        Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) 
c    &                     + d2nepscdn2(D2_RA_RA)*wght
c        Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) 
c    &                     + d2nepscdn2(D2_RA_RB)*wght
c        if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) 
c    &                        + d2nepscdn2(D2_RB_RB)*wght
#endif
c        Amat(n,D1_RA) = Amat(n,D1_RA) +  (H0 + rhoval*dH0dna)*wght
c        if (ipol.eq.2) Amat(n,D1_RB) = Amat(n,D1_RB) +  (H0 +
c    &                  rhoval*dH0dnb)*wght
#ifdef SECOND_DERIV
c        Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &                + (2.d0*dH0dna + rhoval*d2H0dna2)*wght
c        Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &                + (dH0dna + dH0dnb + rhoval*d2H0dnadnb)*wght
c        if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &                + (2.d0*dH0dnb + rhoval*d2H0dnb2)*wght
#endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
c        if (dsqgamma.gt.TOLL)then
c           dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
c           dfAtdg = dfAtdt*dtdg
c           darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
c           dH0dg = GAMMA*(phi**3)*darglogdg/arglog
c           if (ipol.eq.1) then
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c           else
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*wght*2.0d0
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*wght
c           endif
#ifdef SECOND_DERIV
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
c           d2arglogdnadg = BETA/GAMMA*(2.0d0*dtdna*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnadg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdna
c    &                                 +2.0d0*t*dtdna*dfAtdg
c    &                                 +t*t*d2fAtdnadg)
c           d2arglogdnbdg = BETA/GAMMA*(2.0d0*dtdnb*dtdg*fAt
c    &                                 +2.0d0*t*d2tdnbdg*fAt
c    &                                 +2.0d0*t*dtdg*dfAtdnb
c    &                                 +2.0d0*t*dtdnb*dfAtdg
c    &                                 +t*t*d2fAtdnbdg)
c           d2arglogdg2 = BETA/GAMMA*(2.0d0*dtdg*dtdg*fAt
c    &                               +2.0d0*t*d2tdg2*fAt
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +2.0d0*t*dtdg*dfAtdg
c    &                               +t*t*d2fAtdg2)
c           d2H0dnadg = GAMMA*3.0d0*phi**2*dphidna*darglogdg/arglog 
c    &                + GAMMA*phi**3*d2arglogdnadg/arglog
c    &                - GAMMA*phi**3*darglogdg*darglogdna/arglog**2
c           d2H0dnbdg = GAMMA*3.0d0*phi**2*dphidnb*darglogdg/arglog 
c    &                + GAMMA*phi**3*d2arglogdnbdg/arglog
c    &                - GAMMA*phi**3*darglogdg*darglogdnb/arglog**2
c           d2H0dg2 = GAMMA*phi**3*d2arglogdg2/arglog
c    &              - GAMMA*phi**3*darglogdg*darglogdg/arglog**2
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
c    &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
c    &             + (dH0dg + d2H0dnadg*rhoval)*wght
c           Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
c    &             + (dH0dg + d2H0dnbdg*rhoval)*wght
c           Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
c    &             + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*wght
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &             + (dH0dg + d2H0dnbdg*rhoval)*wght
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &             + d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
c    &             + 2.0d0*d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
c    &             + d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
c    &             + 4.0d0*d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
c    &             + 2.0d0*d2H0dg2*rhoval*wght
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &             + d2H0dg2*rhoval*wght
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
#include "nwxc_c_spbe.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_spbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_spbe.F"
#endif
#undef NWAD_PRINT
C>
C> @}
