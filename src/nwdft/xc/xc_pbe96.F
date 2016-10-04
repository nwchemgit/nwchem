#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine xc_cpbe96(tol_rho, cfac, lcfac, nlcfac, rho, delrho,
     &                    Amat, Cmat, nq, ipol, Ec, qwght,
     &                    ldew,ffunc)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
      Subroutine xc_cpbe96_d2(tol_rho, cfac, lcfac, nlcfac,
     &                    rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                    nq, ipol, Ec, qwght, ldew, ffunc)
#else 
      Subroutine xc_cpbe96_d3(tol_rho, cfac, lcfac, nlcfac,
     &                    rho, delrho, Amat, Amat2, Amat3, 
     &                    Cmat, Cmat2, Cmat3,
     &                    nq, ipol, Ec, qwght, ldew, ffunc)
#endif
c
      Implicit none
#include "dft2drv.fh"
#include "dft3drv.fh"
c
c
c     Input and other parameters
c
      integer lnq
      integer ipol, nq, id_cpbe
      double precision dummy(1)
      double precision cfac(*)
      logical lcfac(*), nlcfac(*), ldew
      logical lfac, nlfac, lfacl, nlfacl
      double precision ffunc(*)
      double precision fac, facl
      double precision lqwght
      double precision tol_rho
c
      double precision epsMS, p14f, p14a, p14b
      double precision rho13, rho23, rho43, rho53
      double precision BETA, dBETAdr, d2BETAdr2, d3BETAdr3
      parameter (epsMS=1.e-8)
      parameter (p14a=0.1d0)
      parameter (p14b=0.1778d0)
c
c     Constants in PBE functional
c
      double precision GAMMA, BETAzero, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETAzero = 0.06672455060314922d0)
      parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision TOLL
      double precision EXPTOL
      double precision EPS
      double precision eps_argexp
c
      parameter (TOLL = 1.0D-40)
      parameter (EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
      parameter (eps_argexp = 1.0e-14)

c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision rho_t(3)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
      double precision dsqgamma
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,3)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
#ifdef THIRD_DERIV
c 
c     Sampling Matrices for the XC third derivatives
c
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
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
      double precision Mz, dMzdna, dMzdnb 
      double precision Nzrt, dNzrtdna, dNzrtdnb
      double precision dNzrtdg
      double precision rrho, rrho2
      double precision opz, omz
      double precision rks
      double precision rphi, rphi2, rphi3, rphi4
      double precision expmone, expmone2
      double precision t2, t3, t4
      double precision fAtden2
      double precision phi2, phi3
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
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
      double precision d2Mzdna2, d2Mzdnadnb, d2Mzdnb2
      double precision d2Nzrtdna2, d2Nzrtdnadnb, d2Nzrtdnb2
      double precision d2Nzrtdnadg, d2Nzrtdnbdg, d2Nzrtdg2
      double precision rrho3
      double precision opz2, omz2
      double precision rks2
      double precision rphi5
      double precision expmone3
      double precision fAtden3
#endif
#ifdef THIRD_DERIV
      double precision d3nepscdn3(NCOL_AMAT3)
      double precision d3epscdna3, d3epscdna2dnb, d3epscdnadnb2,
     1                 d3epscdnb3
      double precision d3ksdna3, d3ksdna2dnb, d3ksdnadnb2, d3ksdnb3
      double precision d3zetadna3, d3zetadna2dnb, d3zetadnadnb2,
     1                 d3zetadnb3
      double precision d3phidzeta3, d3phidna3, d3phidna2dnb,
     1                 d3phidnadnb2, d3phidnb3
      double precision d3tdna3, d3tdna2dnb, d3tdnadnb2, d3tdnb3
      double precision d3argexpdna3, d3argexpdna2dnb, d3argexpdnadnb2,
     1                 d3argexpdnb3
      double precision d3Adna3, d3Adna2dnb, d3Adnadnb2, d3Adnb3
      double precision d3fAtnumdt3, d3fAtdendt3, d3fAtnumdt2dA,
     1                 d3fAtdendt2dA, d3fAtnumdtdA2, d3fAtdendtdA2,
     2                 d3fAtnumdA3, d3fAtdendA3   
      double precision d3fAtdt3, d3fAtdt2dA, d3fAtdtdA2, d3fAtdA3   
      double precision d3fAtdna3, d3fAtdna2dnb, d3fAtdnadnb2,
     1                 d3fAtdnb3    
      double precision d3arglogdna3, d3arglogdna2dnb, d3arglogdnadnb2,
     1                 d3arglogdnb3
      double precision d3H0dna3, d3H0dna2dnb, d3H0dnadnb2, d3H0dnb3
      double precision d3tdg3, d3tdna2dg, d3tdnadnbdg, d3tdnbdnadg, 
     1                 d3tdnb2dg, d3tdnadg2, d3tdnbdg2
      double precision d3fAtdg3, d3fAtdna2dg, d3fAtdnadnbdg, 
     1                 d3fAtdnbdnadg,
     2                 d3fAtdnb2dg, d3fAtdnadg2, d3fAtdnbdg2
      double precision d3arglogdg3, d3arglogdna2dg, d3arglogdnadnbdg, 
     1                 d3arglogdnbdnadg,
     2                 d3arglogdnb2dg, d3arglogdnadg2, d3arglogdnbdg2
      double precision d3H0dg3, d3H0dna2dg, d3H0dnadnbdg, d3H0dnbdnadg,
     1                 d3H0dnb2dg, d3H0dnadg2, d3H0dnbdg2
      double precision d3Mzdna3, d3Mzdna2dnb, d3Mzdnadnb2, d3Mzdnb3
      double precision d3Nzrtdna3, d3Nzrtdna2dnb, d3Nzrtdnadnb2, 
     1                 d3Nzrtdnb3
      double precision d3Nzrtdg3, d3Nzrtdnadg2, d3Nzrtdnbdg2,
     1                 d3Nzrtdna2dg, d3Nzrtdnadnbdg, d3Nzrtdnb2dg
      double precision rrho4
      double precision opz3, omz3
      double precision rks3
      double precision rphi6
      double precision expmone4
      double precision fAtden4
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
      p14f=(3d0/(4d0*acos(-1d0)))**(1d0/3d0)  ! gfortran error if this is defined as a parameter

      if (abs(cfac(12)).gt.epsMS) then
c       original PBE correlation functional
        id_cpbe = 1
        fac = cfac(12)
        lfac = lcfac(12)
        nlfac = nlcfac(12)
      elseif (abs(cfac(46)).gt.epsMS) then
c       simplified PBE correlation functional
        id_cpbe = 2
        fac = cfac(46)
        lfac = lcfac(46)
        nlfac = nlcfac(46)
      elseif (abs(cfac(64)).gt.epsMS) then
c       simplified PBE correlation functional
        id_cpbe = 3
        fac = cfac(64)
        lfac = lcfac(64)
        nlfac = nlcfac(64)
      endif
c
c     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
c
      do 20 n = 1, nq
c
c        n and zeta = (na - nb)/n
c
         rhoval = rho(n,1)
         rho_t(1) = rho(n,1)
         if (ipol.eq.2) then
            rho_t(2) = rho(n,2)
            rho_t(3) = rho(n,3)
         endif
         if (rhoval.le.tol_rho) goto 20
c Daniel (7-24-12): gammaval is gamma^2 in the exchange part of the
c calculation.
         if (ipol.eq.1) then
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
         else
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,1) +
     &                 delrho(n,3,2)*delrho(n,3,2) +
     &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,2))
         endif
c Daniel (7-26-12): Is this due to numerical instabilities?
         dsqgamma = max(dsqrt(gammaval),tol_rho)
c Daniel (7-24-12): variable for storing the product of the density
c and the correlation energy per particle, and its derivatives.  This 
c is later used for storing the correlation energy from the LDA part
c of the calculation.
         nepsc = 0.0d0
         dnepscdn(1) = 0.0d0
         if (ipol.eq.2) dnepscdn(2) = 0.0d0
c Daniel (7-29-12): Probably needed for XC third derivatives also
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         d2nepscdn2(D2_RA_RA)=0.0d0
         d2nepscdn2(D2_RA_RB)=0.0d0
         if (ipol.eq.2) d2nepscdn2(D2_RB_RB)=0.0d0
#endif
#ifdef THIRD_DERIV
         d3nepscdn3(D3_RA_RA_RA)=0.0d0
         d3nepscdn3(D3_RA_RA_RB)=0.0d0
         d3nepscdn3(D3_RA_RB_RB)=0.0d0
         if (ipol.eq.2) d3nepscdn3(D3_RB_RB_RB)=0.0d0
#endif
c
c        call for LDA bit
c
         lnq = 1
         lqwght = 1.0d0
c
         if (abs(cfac(1)).gt.EPS) then
            facl = cfac(1)
#if !defined SECOND_DERIV && !defined THIRD_DERIV 
            call xc_vwn_5(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,.false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_5_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_vwn_5_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(3)).gt.EPS) then
            facl = cfac(3)
            lfacl = lcfac(3)
            nlfacl = nlcfac(3)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_p81(tol_rho,facl,lfacl,nlfacl,rho_t,dnepscdn,
     &         lnq,ipol,nepsc,lqwght,.false.,dummy)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call xc_p81_d2(tol_rho,facl,lfacl,nlfacl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_p81_d3(tol_rho,facl,lfacl,nlfacl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
c Daniel (2-19-13): This is the default choice when PBE96 correlation
c is used.
         if (abs(cfac(6)).gt.EPS.or.lfac) then
            facl = cfac(6)
            if (abs(facl).lt.EPS) facl = 1.0d0
            lfacl = lcfac(6)
            nlfacl = nlcfac(6)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_pw91lda(tol_rho,facl,lfacl,nlfacl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_pw91lda_d2(tol_rho,facl,lfacl,nlfacl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_pw91lda_d3(tol_rho,facl,lfacl,nlfacl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(7)).gt.EPS) then
            facl = cfac(7)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_vwn_1_rpa(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_1_rpa_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_vwn_1_rpa_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(8)).gt.EPS) then
            facl = cfac(8)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_vwn_1(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_1_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_vwn_1_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(9)).gt.EPS) then
            facl = cfac(9)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_vwn_2(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_2_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_vwn_2_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(10)).gt.EPS) then
            facl = cfac(10)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_vwn_3(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_3_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else 
            call xc_vwn_3_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c
         if (abs(cfac(11)).gt.EPS) then
            facl = cfac(11)
#if !defined SECOND_DERIV && !defined THIRD_DERIV
            call xc_vwn_4(tol_rho,facl,rho_t,
     &         dnepscdn,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
            call xc_vwn_4_d2(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#else
            call xc_vwn_4_d3(tol_rho,facl,rho_t,
     &         dnepscdn,d2nepscdn2,d3nepscdn3,lnq,ipol,nepsc,lqwght,
     &         .false.,dummy)
#endif
         endif
c        ==================
c        PBE non-local part
c        ==================
         if(abs(nepsc).lt.tol_rho*tol_rho) goto 20
c
c        epsilon_c = n*epsilon_c / n
c
c Daniel (7-24-12): Regardless of the form of the HEG correlation used
c above, nepsc = eps*qwght*rhoval*cfac.  We need to gather the 
c contributions from different spins.  The derivatives given are
c just the derivative of the epsilon from the local (LDA) part.  Since
c epsilon and its derivatives are multiplied by the density in the
c routines above, we end up taking derivatives of n*epsilon_c / n.
         rrho = 1.0d0/rhoval
         rrho2 = rrho*rrho
c
         epsc = nepsc*rrho
         if (ipol.eq.1) then
            depscdna = dnepscdn(1)*rrho-nepsc*rrho2
            depscdnb = depscdna
         else
            depscdna = dnepscdn(1)*rrho-nepsc*rrho2
            depscdnb = dnepscdn(2)*rrho-nepsc*rrho2
         endif
c Daniel (7-29-12): Probably needed for XC third derivatives also
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         rrho3 = rrho2*rrho
         if (ipol.eq.1) then
            d2epscdna2   = d2nepscdn2(D2_RA_RA)*rrho
     &                     -dnepscdn(1)*rrho2
     &                     -dnepscdn(1)*rrho2
     &                     +2.0d0*nepsc*rrho3
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)*rrho
     &                     -dnepscdn(1)*rrho2
     &                     -dnepscdn(1)*rrho2
     &                     +2.0d0*nepsc*rrho3
            d2epscdnb2   = d2epscdna2
         else
            d2epscdna2   = d2nepscdn2(D2_RA_RA)*rrho
     &                     -dnepscdn(1)*rrho2
     &                     -dnepscdn(1)*rrho2
     &                     +2.0d0*nepsc*rrho3
            d2epscdnadnb = d2nepscdn2(D2_RA_RB)*rrho
     &                     -dnepscdn(1)*rrho2
     &                     -dnepscdn(2)*rrho2
     &                     +2.0d0*nepsc*rrho3
            d2epscdnb2   = d2nepscdn2(D2_RB_RB)*rrho
     &                     -dnepscdn(2)*rrho2
     &                     -dnepscdn(2)*rrho2
     &                     +2.0d0*nepsc*rrho3
         endif
#endif
c
#ifdef THIRD_DERIV
c Subtle point here!  We don't actually calculate the first derivative
c of the LDA correlation energy with respect to the beta spin density
c for a restricted calculation, so we need to use that na = nb
         rrho4 = rrho3*rrho
         if (ipol.eq.1) then
            d3epscdna3    = d3nepscdn3(D3_RA_RA_RA)*rrho
     1                    - d2nepscdn2(D2_RA_RA)*rrho2
     2                    - 2.0d0*d2nepscdn2(D2_RA_RA)*rrho2
     3                    + 6.0d0*dnepscdn(1)*rrho3
     4                    - 6.0d0*nepsc*rrho4
            d3epscdna2dnb = d3nepscdn3(D3_RA_RA_RB)*rrho
     1                    - d2nepscdn2(D2_RA_RA)*rrho2
     2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
     3                    + 4.0d0*dnepscdn(1)*rrho3
     4                    + 2.0d0*dnepscdn(1)*rrho3
     5                    - 6.0d0*nepsc*rrho4
            d3epscdnadnb2 = d3epscdna2dnb
            d3epscdnb3    = d3epscdna3
         else
            d3epscdna3    = d3nepscdn3(D3_RA_RA_RA)*rrho
     1                    - 3.0d0*d2nepscdn2(D2_RA_RA)*rrho2
     2                    + 6.0d0*dnepscdn(1)*rrho3
     3                    - 6.0d0*nepsc*rrho4
            d3epscdna2dnb = d3nepscdn3(D3_RA_RA_RB)*rrho
     1                    - d2nepscdn2(D2_RA_RA)*rrho2
     2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
     3                    + 4.0d0*dnepscdn(1)*rrho3
     4                    + 2.0d0*dnepscdn(2)*rrho3
     5                    - 6.0d0*nepsc*rrho4
            d3epscdnadnb2 = d3nepscdn3(D3_RA_RB_RB)*rrho
     1                    - d2nepscdn2(D2_RB_RB)*rrho2
     2                    - 2.0d0*d2nepscdn2(D2_RA_RB)*rrho2
     3                    + 2.0d0*dnepscdn(1)*rrho3
     4                    + 4.0d0*dnepscdn(2)*rrho3
     5                    - 6.0d0*nepsc*rrho4
            d3epscdnb3    = d3nepscdn3(D3_RB_RB_RB)*rrho
     1                    - 3.0d0*d2nepscdn2(D2_RB_RB)*rrho2
     2                    + 6.0d0*dnepscdn(2)*rrho3
     3                    - 6.0d0*nepsc*rrho4
         endif
#endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
c Daniel (7-24-12): Thomas-Fermi screening vector
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
         dksdna = (1.0d0/6.0d0)*ks*rrho
         dksdnb = dksdna
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         d2ksdna2   = (1.0d0/6.0d0)*dksdna*rrho
     &              - (1.0d0/6.0d0)*ks*rrho2
         d2ksdnadnb = d2ksdna2
         d2ksdnb2   = d2ksdna2
#endif
#ifdef THIRD_DERIV
         d3ksdna3    = (1.0d0/6.0d0)*d2ksdna2*rrho
     1               - (1.0d0/3.0d0)*dksdna*rrho2
     2               + (1.0d0/3.0d0)*ks*rrho3
         d3ksdna2dnb = d3ksdna3
         d3ksdnadnb2 = d3ksdna3
         d3ksdnb3    = d3ksdna3
#endif
c
c        zeta = (na-nb)/n and its derivs
c
c Daniel (7-24-12): Spin-polarization
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,2)-rho(n,3))*rrho
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
         if (ipol.eq.1) then
            dzetadna = 1.0d0*rrho
            dzetadnb = -1.0d0*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2zetadna2   = -2.0d0*rrho2
            d2zetadnadnb = 0.0d0
c Daniel (8-27-12): This should be a positive sign. 
c OLD
c            d2zetadnb2   = -2.0d0/(rhoval**2)
c OLD
            d2zetadnb2   = 2.0d0*rrho2
#endif
c Daniel (7-29-12): Third derivative stuff
#ifdef THIRD_DERIV
            d3zetadna3    = 6.0d0*rrho3
            d3zetadna2dnb = 2.0d0*rrho3
            d3zetadnadnb2 = -2.0d0*rrho3
            d3zetadnb3    = -6.0d0*rrho3
#endif
         else
            dzetadna =  2.0d0*rho(n,3)*rrho2
            dzetadnb = -2.0d0*rho(n,2)*rrho2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2zetadna2   = -4.0d0*rho(n,3)*rrho3
            d2zetadnadnb = 2.0d0*(rho(n,2)-rho(n,3))*rrho3
c Daniel (8-28-12): This should be a positive sign
c OLD
c            d2zetadnb2   = -4.0d0*rho(n,2)/(rhoval**3)
c OLD
            d2zetadnb2   = 4.0d0*rho(n,2)*rrho3
#endif
c Daniel (7-29-12): Third derivative stuff
#ifdef THIRD_DERIV
            d3zetadna3    = 12.0d0*rho(n,3)*rrho4 
            d3zetadna2dnb = -4.0d0*(rho(n,2)-2.0d0*rho(n,3))*
     1                      rrho4
            d3zetadnadnb2 = 4.0d0*(rho(n,2)-2.0d0*rho(n,3))*
     1                      rrho4
            d3zetadnb3    = -12.0d0*rho(n,2)*rrho4
#endif
         endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
c Daniel (7-24-12): Spin-scaling factor
         opz = 1.0d0 + zeta
         omz = 1.0d0 - zeta
c
         phi = 0.5d0*(opz**(2.0d0/3.0d0)
     &               +omz**(2.0d0/3.0d0))
         if (omz.lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &             opz**(2.0d0/3.0d0)/opz)
         else if (opz.lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &            -omz**(2.0d0/3.0d0)/omz)
         else
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &         opz**(2.0d0/3.0d0)/opz
     &        -omz**(2.0d0/3.0d0)/omz)
         endif
         dphidna = dphidzeta*dzetadna
         dphidnb = dphidzeta*dzetadnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         opz2 = opz*opz
         omz2 = omz*omz
c
         if (omz.lt.tol_rho) then
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         opz**(2.0d0/3.0d0)/opz2)
         else if (opz.lt.tol_rho) then
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         omz**(2.0d0/3.0d0)/omz2)
         else
            d2phidzeta2 = -(1.0d0/9.0d0)*(
     &         opz**(2.0d0/3.0d0)/opz2
     &        +omz**(2.0d0/3.0d0)/omz2)
         endif
         d2phidna2   = d2phidzeta2*dzetadna*dzetadna
     &               + dphidzeta*d2zetadna2 
         d2phidnadnb = d2phidzeta2*dzetadna*dzetadnb
     &               + dphidzeta*d2zetadnadnb
         d2phidnb2   = d2phidzeta2*dzetadnb*dzetadnb
     &               + dphidzeta*d2zetadnb2 
#endif
#ifdef THIRD_DERIV
c Daniel (7-30-12): Testing is done here to prevent numerical problems
         opz3 = opz2*opz
         omz3 = omz2*omz
c
         if (omz.lt.tol_rho) then
            d3phidzeta3 = (4.0d0/27.0d0)*(
     1         opz**(2.0d0/3.0d0)/opz3)
         else if (opz.lt.tol_rho) then
            d3phidzeta3 = (4.0d0/27.0d0)*(
     1         -omz**(2.0d0/3.0d0)/omz3)
         else
            d3phidzeta3 = (4.0d0/27.0d0)*(
     1         opz**(2.0d0/3.0d0)/opz3
     2         -omz**(2.0d0/3.0d0)/omz3)
         endif
         d3phidna3    = d3phidzeta3*dzetadna*dzetadna*dzetadna
     1                + 3.0d0*d2phidzeta2*d2zetadna2*dzetadna
     2                + dphidzeta*d3zetadna3
         d3phidna2dnb = d3phidzeta3*dzetadna*dzetadna*dzetadnb
     1                + 2.0d0*d2phidzeta2*d2zetadnadnb*dzetadna
     2                + d2phidzeta2*d2zetadna2*dzetadnb
     3                + dphidzeta*d3zetadna2dnb
         d3phidnadnb2 = d3phidzeta3*dzetadna*dzetadnb*dzetadnb
     1                + 2.0d0*d2phidzeta2*d2zetadnadnb*dzetadnb
     2                + d2phidzeta2*d2zetadnb2*dzetadna
     3                + dphidzeta*d3zetadnadnb2
         d3phidnb3    = d3phidzeta3*dzetadnb*dzetadnb*dzetadnb
     1                + 3.0d0*d2phidzeta2*d2zetadnb2*dzetadnb
     2                + dphidzeta*d3zetadnb3
#endif
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
         rks = 1.0d0/ks
         rphi = 1.0d0/phi
         rphi2 = rphi*rphi
         rphi3 = rphi2*rphi
c
         t = dsqgamma*rphi*rks*rrho/2.0d0
         dtdna = -t*rrho-t*rphi*dphidna-t*rks*dksdna
         dtdnb = -t*rrho-t*rphi*dphidnb-t*rks*dksdnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         rks2 = rks*rks
c
         d2tdna2 = - dtdna*rrho
     &           + t*rrho2
     &           - dtdna*rphi*dphidna
     &           + t*rphi2*dphidna*dphidna
     &           - t*rphi*d2phidna2         
     &           - dtdna*rks*dksdna
     &           + t*rks2*dksdna*dksdna   
     &           - t*rks*d2ksdna2           
         d2tdnadnb = - dtdnb*rrho
     &           + t*rrho2
     &           - dtdnb*rphi*dphidna
     &           + t*rphi2*dphidna*dphidnb
     &           - t*rphi*d2phidnadnb
     &           - dtdnb*rks*dksdna
     &           + t*rks2*dksdna*dksdnb
     &           - t*rks*d2ksdnadnb
c
         d2tdnb2 = - dtdnb*rrho
     &           + t*rrho2
     &           - dtdnb*rphi*dphidnb
     &           + t*rphi2*dphidnb*dphidnb
     &           - t*rphi*d2phidnb2
     &           - dtdnb*rks*dksdnb
     &           + t*rks2*dksdnb*dksdnb
     &           - t*rks*d2ksdnb2
#endif
c
#ifdef THIRD_DERIV
         rks3 = rks2*rks
c        
         d3tdna3    = -2.0d0*t*rphi3*dphidna*dphidna*dphidna
     1              - t*rks*d3ksdna3
     2              - 2.0d0*t*rrho3
     3              + 3.0d0*t*rphi2*d2phidna2*dphidna
     4              - t*rphi*d3phidna3
     5              - 2.0d0*t*rks3*dksdna*dksdna*dksdna
     6              + 3.0d0*t*rks2*d2ksdna2*dksdna
     7              + (dphidna*dphidna*rphi2
     8                - d2phidna2*rphi
     9                + dksdna*dksdna*rks2
     A                - d2ksdna2*rks
     B                + rrho2)*dtdna*2.0d0
     C              - (dphidna*rphi + dksdna*rks
     D                + rrho)*d2tdna2
         d3tdna2dnb = -2.0d0*t*rphi3*dphidna*dphidna*dphidnb
     1              - t*rks*d3ksdna2dnb
     2              - 2.0d0*t*rrho3
     3              + 2.0d0*t*rphi2*d2phidnadnb*dphidna
     4              + t*rphi2*d2phidna2*dphidnb
     5              - t*rphi*d3phidna2dnb
     6              - 2.0d0*t*rks3*dksdna*dksdna*dksdnb
     7              + 2.0d0*t*rks2*d2ksdnadnb*dksdna
     8              + t*rks2*d2ksdna2*dksdnb
     9              + (dphidna*dphidna*rphi2
     A                - d2phidna2*rphi
     B                + dksdna*dksdna*rks2
     C                - d2ksdna2*rks
     D                + rrho2)*dtdnb
     E              + (dphidna*dphidnb*rphi2
     F                - d2phidnadnb*rphi
     G                + dksdna*dksdnb*rks2
     H                - d2ksdnadnb*rks
     I                + rrho2)*dtdna
     J              - (dphidna*rphi + dksdna*rks
     K                + rrho)*d2tdnadnb
         d3tdnadnb2 = -2.0d0*t*rphi3*dphidna*dphidnb*dphidnb
     1              - t*rks*d3ksdnadnb2
     2              - 2.0d0*t*rrho3
     3              + 2.0d0*t*rphi2*d2phidnadnb*dphidnb
     4              + t*rphi2*d2phidnb2*dphidna
     5              - t*rphi*d3phidnadnb2
     6              - 2.0d0*t*rks3*dksdna*dksdnb*dksdnb
     7              + 2.0d0*t*rks2*d2ksdnadnb*dksdnb
     8              + t*rks2*d2ksdnb2*dksdna
     9              + (dphidna*dphidnb*rphi2
     A                - d2phidnadnb*rphi
     B                + dksdna*dksdnb*rks2
     C                - d2ksdnadnb*rks
     D                + rrho2)*dtdnb
     E              + (dphidna*dphidnb*rphi2
     F                - d2phidnadnb*rphi
     G                + dksdna*dksdnb*rks2
     H                - d2ksdnadnb*rks
     I                + rrho2)*dtdnb
     J              - (dphidna*rphi + dksdna*rks
     K                + rrho)*d2tdnb2
         d3tdnb3    = -2.0d0*t*rphi3*dphidnb*dphidnb*dphidnb
     1              - t*rks*d3ksdnb3
     2              - 2.0d0*t*rrho3
     3              + 3.0d0*t*rphi2*d2phidnb2*dphidnb
     4              - t*rphi*d3phidnb3
     5              - 2.0d0*t*rks3*dksdnb*dksdnb*dksdnb
     6              + 3.0d0*t*rks2*d2ksdnb2*dksdnb
     7              + (dphidnb*dphidnb*rphi2
     8                - d2phidnb2*rphi
     9                + dksdnb*dksdnb*rks2
     A                - d2ksdnb2*rks
     B                + rrho2)*dtdnb*2.0d0
     C              - (dphidnb*rphi + dksdnb*rks
     D                + rrho)*d2tdnb2
#endif
c
c        { ... } in A (see below) and its derivs
c
         rphi4 = rphi3*rphi
c
         argexp = -epsc*rphi3/GAMMA
         dargexpdna = -depscdna*rphi3/GAMMA
     &                +3.0d0*epsc*rphi4*dphidna/GAMMA
         dargexpdnb = -depscdnb*rphi3/GAMMA
     &                +3.0d0*epsc*rphi4*dphidnb/GAMMA
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         rphi5 = rphi4*rphi
         d2argexpdna2 = -d2epscdna2*rphi3/GAMMA
     &        +3.0d0*depscdna*rphi4*dphidna/GAMMA
     &        +3.0d0*depscdna*rphi4*dphidna/GAMMA
     &        -12.0d0*epsc*rphi5*dphidna*dphidna/GAMMA
     &        +3.0d0*epsc*rphi4*d2phidna2/GAMMA
         d2argexpdnadnb = -d2epscdnadnb*rphi3/GAMMA
     &        +3.0d0*depscdna*rphi4*dphidnb/GAMMA
     &        +3.0d0*depscdnb*rphi4*dphidna/GAMMA
     &        -12.0d0*epsc*rphi5*dphidna*dphidnb/GAMMA
     &        +3.0d0*epsc*rphi4*d2phidnadnb/GAMMA
         d2argexpdnb2 = -d2epscdnb2*rphi3/GAMMA
     &        +3.0d0*depscdnb*rphi4*dphidnb/GAMMA
     &        +3.0d0*depscdnb*rphi4*dphidnb/GAMMA
     &        -12.0d0*epsc*rphi5*dphidnb*dphidnb/GAMMA
     &        +3.0d0*epsc*rphi4*d2phidnb2/GAMMA
#endif
c
#ifdef THIRD_DERIV
         rphi6 = rphi5*rphi
c
         d3argexpdna3    = -12.0d0*rphi5*
     1                     ( 3.0d0*dphidna*dphidna*depscdna
     2                     + 3.0d0*dphidna*d2phidna2*epsc)/GAMMA
     3                   + 60.0d0*rphi6*
     4                     dphidna*dphidna*dphidna*epsc/GAMMA
     5                   + 3.0d0*rphi4*
     6                     ( 3.0d0*d2phidna2*depscdna
     7                     + 3.0d0*dphidna*d2epscdna2
     8                     + d3phidna3*epsc)/GAMMA
     9                   - d3epscdna3*rphi3/GAMMA
         d3argexpdna2dnb = -12.0d0*rphi5*
     1                     ( 2.0d0*dphidna*dphidnb*depscdna
     2                     + dphidna*dphidna*depscdnb
     3                     + 2.0d0*dphidna*d2phidnadnb*epsc
     4                     + dphidnb*d2phidna2*epsc)/GAMMA
     5                   + 60.0d0*rphi6*
     6                     dphidna*dphidna*dphidnb*epsc/GAMMA
     7                   + 3.0d0*rphi4*
     8                     ( 2.0d0*d2phidnadnb*depscdna
     9                     + d2phidna2*depscdnb
     A                     + 2.0d0*dphidna*d2epscdnadnb
     B                     + dphidnb*d2epscdna2
     C                     + d3phidna2dnb*epsc)/GAMMA
     D                   - d3epscdna2dnb*rphi3/GAMMA
         d3argexpdnadnb2 = -12.0d0*rphi5*
     1                     ( 2.0d0*dphidna*dphidnb*depscdnb
     2                     + dphidnb*dphidnb*depscdna
     3                     + 2.0d0*dphidnb*d2phidnadnb*epsc
     4                     + dphidna*d2phidnb2*epsc)/GAMMA
     5                   + 60.0d0*rphi6*
     6                     dphidna*dphidnb*dphidnb*epsc/GAMMA
     7                   + 3.0d0*rphi4*
     8                     ( 2.0d0*d2phidnadnb*depscdnb
     9                     + d2phidnb2*depscdna
     A                     + 2.0d0*dphidnb*d2epscdnadnb
     B                     + dphidna*d2epscdnb2
     C                     + d3phidnadnb2*epsc)/GAMMA
     D                   - d3epscdnadnb2*rphi3/GAMMA
         d3argexpdnb3    = -12.0d0*rphi5*
     1                     ( 3.0d0*dphidnb*dphidnb*depscdnb
     2                     + 3.0d0*dphidnb*d2phidnb2*epsc)/GAMMA
     3                   + 60.0d0*rphi6*
     4                     dphidnb*dphidnb*dphidnb*epsc/GAMMA
     5                   + 3.0d0*rphi4*
     6                     ( 3.0d0*d2phidnb2*depscdnb
     7                     + 3.0d0*dphidnb*d2epscdnb2
     8                     + d3phidnb3*epsc)/GAMMA
     9                   - d3epscdnb3*rphi3/GAMMA
#endif
c
c        A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
c Daniel (7-24-12): "A" function
c     this avoids expmone to get a zero value that will result in 1./0.
         if (dabs(argexp).lt.eps_argexp) then
            expinA=dexp(eps_argexp)
         elseif (dabs(argexp).lt.EXPTOL) then
            expinA=dexp(argexp)
         else
            expinA=0.0d0
         endif
         expmone = expinA - 1.0d0
         expmone2 = expmone*expmone
c
c        revTPSS variant of cPBE functional
c
         if (id_cpbe.eq.3) then
c
           rho13 = rhoval**(1d0/3d0)
           rho23 = rho13*rho13
           rho43 = rho23*rho23
           rho53 = rho43*rho13
c
           BETA = BETAzero*(rho13+p14a*p14f)/(rho13+p14b*p14f)
           dBETAdr = BETAzero*(1d0/3d0)*(p14b*p14f-p14a*p14f) /
     &                      (rho23*(rho13+p14b*p14f)**2)
           d2BETAdr2 = BETAzero*(-2d0/9d0)*(p14b*p14f-p14a*p14f)*
     &         (2d0/rho43 + p14b*p14f/rho53) / (rho13+p14b*p14f)**3
#ifdef THIRD_DERIV
           stop 'not yet implemented 3rd derivative of BETA vs r'
#endif
c
         else
           BETA = BETAzero
           dBETAdr = 0.0d0
           d2BETAdr2 = 0.0d0
           d3BETAdr3 = 0.0d0
         endif
c
         A = BETA/GAMMA/expmone
         dAdna = -BETA/GAMMA*dargexpdna*expinA/expmone2
         dAdnb = -BETA/GAMMA*dargexpdnb*expinA/expmone2
c
         if (id_cpbe.eq.3) then
           dAdna = dAdna + dBETAdr/GAMMA/expmone
           dAdnb = dAdnb + dBETAdr/GAMMA/expmone
         endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         expmone3 = expmone2*expmone
c
         d2Adna2   = -BETA/GAMMA*d2argexpdna2
     &               *expinA/expmone2
     &             - BETA/GAMMA*dargexpdna
     &               *dargexpdna*expinA/expmone2
     &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdna
     &               *expinA*expinA/expmone3
         d2Adnadnb  = -BETA/GAMMA*d2argexpdnadnb
     &               *expinA/expmone2
     &             - BETA/GAMMA*dargexpdna
     &               *dargexpdnb*expinA/expmone2
     &             + 2.0d0*BETA/GAMMA*dargexpdna*dargexpdnb
     &               *expinA*expinA/expmone3
         d2Adnb2   = -BETA/GAMMA*d2argexpdnb2
     &               *expinA/expmone2
     &             - BETA/GAMMA*dargexpdnb
     &               *dargexpdnb*expinA/expmone2
     &             + 2.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb
     &               *expinA*expinA/expmone3
c
         if (id_cpbe.eq.3) then
           d2Adna2 = d2Adna2 + d2BETAdr2/GAMMA/expmone
     &             - dBETAdr/GAMMA*dargexpdna*expinA/expmone2
           d2Adnadnb = d2Adnadnb + d2BETAdr2/GAMMA/expmone
     &             - dBETAdr/GAMMA*dargexpdna*expinA/expmone2
     &             - dBETAdr/GAMMA*dargexpdnb*expinA/expmone2
           d2Adnb2 = d2Adnb2 + d2BETAdr2/GAMMA/expmone
     &             - dBETAdr/GAMMA*dargexpdnb*expinA/expmone2
         endif
c
#endif
#ifdef THIRD_DERIV
         expmone4 = expmone3*expmone
c
         d3Adna3    = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdna*
     1                dargexpdna*expinA*expinA*expinA/
     2                expmone4
     3              + 6.0d0*BETA/GAMMA*
     4                ( dargexpdna*d2argexpdna2
     5                + dargexpdna*dargexpdna*dargexpdna)*
     6                expinA*expinA/expmone3
     7              - BETA/GAMMA*
     8                ( d3argexpdna3 + 3.0d0*d2argexpdna2*dargexpdna
     9                + dargexpdna*dargexpdna*dargexpdna)*
     A                expinA/expmone2
         d3Adna2dnb = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdna*
     1                dargexpdnb*expinA*expinA*expinA/
     2                expmone4
     3              + 2.0d0*BETA/GAMMA*
     4                ( 2.0d0*dargexpdna*d2argexpdnadnb
     5                + 2.0d0*dargexpdna*dargexpdna*dargexpdnb
     6                + dargexpdnb*d2argexpdna2
     7                + dargexpdna*dargexpdna*dargexpdnb)*
     8                expinA*expinA/expmone3
     9              - BETA/GAMMA*
     A                ( d3argexpdna2dnb
     B                + 2.0d0*d2argexpdnadnb*dargexpdna
     C                + d2argexpdna2*dargexpdnb
     D                + dargexpdna*dargexpdna*dargexpdnb)*
     E                expinA/expmone2
         d3Adnadnb2 = -6.0d0*BETA/GAMMA*dargexpdna*dargexpdnb*
     1                dargexpdnb*expinA*expinA*expinA/
     2                expmone4
     3              + 2.0d0*BETA/GAMMA*
     4                ( 2.0d0*dargexpdnb*d2argexpdnadnb
     5                + 2.0d0*dargexpdna*dargexpdnb*dargexpdnb
     6                + dargexpdna*d2argexpdnb2
     7                + dargexpdna*dargexpdnb*dargexpdnb)*
     8                expinA*expinA/expmone3
     9              - BETA/GAMMA*
     A                ( d3argexpdnadnb2
     B                + 2.0d0*d2argexpdnadnb*dargexpdnb
     C                + d2argexpdnb2*dargexpdna
     D                + dargexpdna*dargexpdnb*dargexpdnb)*
     E                expinA/expmone2
         d3Adnb3    = -6.0d0*BETA/GAMMA*dargexpdnb*dargexpdnb*
     1                dargexpdnb*expinA*expinA*expinA/
     2                expmone4
     3              + 6.0d0*BETA/GAMMA*
     4                ( dargexpdnb*d2argexpdnb2
     5                + dargexpdnb*dargexpdnb*dargexpdnb)*
     6                expinA*expinA/expmone3
     7              - BETA/GAMMA*
     8                ( d3argexpdnb3 + 3.0d0*d2argexpdnb2*dargexpdnb
     9                + dargexpdnb*dargexpdnb*dargexpdnb)*
     A                expinA/expmone2
#endif
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
         t2 = t*t
         t3 = t2*t
         t4 = t3*t
c
         if (id_cpbe.eq.1.or.id_cpbe.eq.3) then
           fAtnum = 1.0d0+A*t2
           fAtden = 1.0d0+A*t2+A*A*t4
           fAtden2 = fAtden*fAtden
           fAt = fAtnum/fAtden
           dfAtnumdt = 2.0d0*A*t
           dfAtnumdA = t2
           dfAtdendt = 2.0d0*A*t+4.0d0*A*A*t3
           dfAtdendA = t2+2.0d0*A*t4
           dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/fAtden2
           dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/fAtden2
         else
           fAtnum = 1.0d0
           fAtden = 1.0d0+A*t**2
           fAt = fAtnum/fAtden
           dfAtdendt = 2.0d0*A*t
           dfAtdendA = t**2
           dfAtdt = (-fAtnum*dfAtdendt)/(fAtden**2)
           dfAtdA = (-fAtnum*dfAtdendA)/(fAtden**2)
         endif
         dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
         dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         fAtden3 = fAtden2*fAtden
c
         if (id_cpbe.eq.1.or.id_cpbe.eq.3) then
           d2fAtnumdt2 = 2.0d0*A
           d2fAtdendt2 = 2.0d0*A+12.0d0*A**2*t2
           d2fAtnumdtdA = 2.0d0*t
           d2fAtdendtdA = 2.0d0*t+8.0d0*A*t3
           d2fAtnumdA2 = 0.0d0
           d2fAtdendA2 = 2.0d0*t4
           d2fAtdt2  = (d2fAtnumdt2*fAtden-fAtnum*d2fAtdendt2)
     &                 /(fAtden2)
     &                 -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
     &                 /(fAtden3)*dfAtdendt
           d2fAtdtdA = (d2fAtnumdtdA*fAtden+dfAtnumdt*dfAtdendA
     &                 -dfAtnumdA*dfAtdendt-fAtnum*d2fAtdendtdA)
     &                 /(fAtden2)
     &                 -2.0d0*(dfAtnumdt*fAtden-fAtnum*dfAtdendt)
     &                 /(fAtden3)*dfAtdendA
           d2fAtdA2  = (d2fAtnumdA2*fAtden-fAtnum*d2fAtdendA2)
     &                 /(fAtden2)
     &                 -2.0d0*(dfAtnumdA*fAtden-fAtnum*dfAtdendA)
     &                 /(fAtden3)*dfAtdendA
         else
           d2fAtdendt2 = 2.0d0*A
           d2fAtdendtdA = 2.0d0*t
           d2fAtdendA2 = 0.0d0
           d2fAtdt2  = (-fAtnum*d2fAtdendt2)
     &                 /(fAtden**2)
     &                 -2.0d0*(-fAtnum*dfAtdendt)
     &                 /(fAtden**3)*dfAtdendt
           d2fAtdtdA = (-fAtnum*d2fAtdendtdA)
     &                 /(fAtden**2)
     &                 -2.0d0*(-fAtnum*dfAtdendt)
     &                 /(fAtden**3)*dfAtdendA
           d2fAtdA2  = (-fAtnum*d2fAtdendA2)
     &                 /(fAtden**2)
     &                 -2.0d0*(-fAtnum*dfAtdendA)
     &                 /(fAtden**3)*dfAtdendA
         endif
c
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
#ifdef THIRD_DERIV
         fAtden4 = fAtden3*fAtden
c
         d3fAtnumdt3   = 0.0d0
         d3fAtdendt3   = 24.0d0*A**2*t
         d3fAtnumdt2dA = 2.0d0
         d3fAtdendt2dA = 2.0d0+24.0d0*A*t2
         d3fAtnumdtdA2 = 0.0d0
         d3fAtdendtdA2 = 8.0d0*t3
         d3fAtnumdA3   = 0.0d0
         d3fAtdendA3   = 0.0d0
c Here we pay attention to which derivatives are zero.
         d3fAtdt3   = -( fAtnum*d3fAtdendt3
     1                 + 3.0d0*dfAtnumdt*d2fAtdendt2
     2                 + 3.0d0*d2fAtnumdt2*dfAtdendt)/
     3                 (fAtden2)
     4                + 6.0d0*( fAtnum*d2fAtdendt2*dfAtdendt
     5                        + dfAtnumdt*dfAtdendt*dfAtdendt)/
     6                 (fAtden3)
     7                - 6.0d0*fAtnum*dfAtdendt*dfAtdendt*dfAtdendt/
     8                 (fAtden4)
         d3fAtdt2dA = -( fAtnum*d3fAtdendt2dA
     1                 + 2.0d0*dfAtnumdt*d2fAtdendtdA
     2                 + dfAtnumdA*d2fAtdendt2
     3                 + 2.0d0*d2fAtnumdtdA*dfAtdendt
     4                 + d2fAtnumdt2*dfAtdendA)/
     5                 (fAtden2)
     6                + 2.0d0*( 2.0d0*fAtnum*d2fAtdendtdA*dfAtdendt
     7                        + fAtnum*d2fAtdendt2*dfAtdendA
     8                        + 2.0d0*dfAtnumdt*dfAtdendA*dfAtdendt
     9                        + dfAtnumdA*dfAtdendt*dfAtdendt)/
     A                 (fAtden3)
     B                - 6.0d0*fAtnum*dfAtdendt*dfAtdendt*dfAtdendA/
     C                 (fAtden4)
     D                + d3fAtnumdt2dA/fAtden
         d3fAtdtdA2 = -( fAtnum*d3fAtdendtdA2
     1                 + 2.0d0*dfAtnumdA*d2fAtdendtdA
     2                 + dfAtnumdt*d2fAtdendA2
     3                 + 2.0d0*d2fAtnumdtdA*dfAtdendA)/
     4                 (fAtden2)
     5                + 2.0d0*( 2.0d0*fAtnum*d2fAtdendtdA*dfAtdendA
     6                        + fAtnum*d2fAtdendA2*dfAtdendt
     7                        + 2.0d0*dfAtnumdA*dfAtdendA*dfAtdendt
     8                        + dfAtnumdt*dfAtdendA*dfAtdendA)/
     9                 (fAtden3)
     A                - 6.0d0*fAtnum*dfAtdendt*dfAtdendA*dfAtdendA/
     B                 (fAtden4)
         d3fAtdA3   = -3.0d0*dfAtnumdA*d2fAtdendA2/
     1                 (fAtden2)
     2                + 6.0d0*( fAtnum*d2fAtdendA2*dfAtdendA
     3                        + dfAtnumdA*dfAtdendA*dfAtdendA)/
     3                 (fAtden3)
     4                - 6.0d0*fAtnum*dfAtdendA*dfAtdendA*dfAtdendA/
     5                 (fAtden4)
c
         d3fAtdna3    = d3fAtdA3*dAdna*dAdna*dAdna
     1                + 3.0d0*d3fAtdtdA2*dAdna*dAdna*dtdna
     2                + 3.0d0*d2fAtdA2*d2Adna2*dAdna
     3                + 3.0d0*d3fAtdt2dA*dAdna*dtdna*dtdna
     4                + 3.0d0*d2fAtdtdA*(d2Adna2*dtdna + dAdna*d2tdna2)
     5                + dfAtdA*d3Adna3 + d3fAtdt3*dtdna*dtdna*dtdna
     6                + dfAtdt*d3tdna3
     7                + 3.0d0*d2fAtdt2*d2tdna2*dtdna
         d3fAtdna2dnb = d3fAtdA3*dAdna*dAdna*dAdnb
     1                + 2.0d0*d3fAtdtdA2*dAdna*dAdnb*dtdna
     2                + d3fAtdtdA2*dAdna*dAdna*dtdnb
     3                + 2.0d0*d2fAtdA2*d2Adnadnb*dAdna
     4                + d2fAtdA2*d2Adna2*dAdnb
     5                + 2.0d0*d3fAtdt2dA*dAdna*dtdna*dtdnb
     6                + d3fAtdt2dA*dAdnb*dtdna*dtdna
     7                + 2.0d0*d2fAtdtdA*(d2Adnadnb*dtdna
     8                                 + dAdna*d2tdnadnb)
     9                + d2fAtdtdA*(d2Adna2*dtdnb + dAdnb*d2tdna2)
     A                + dfAtdA*d3Adna2dnb + d3fAtdt3*dtdna*dtdna*dtdnb
     B                + dfAtdt*d3tdna2dnb
     C                + 2.0d0*d2fAtdt2*d2tdnadnb*dtdna
     D                + d2fAtdt2*d2tdna2*dtdnb
         d3fAtdnadnb2 = d3fAtdA3*dAdna*dAdnb*dAdnb
     1                + 2.0d0*d3fAtdtdA2*dAdna*dAdnb*dtdnb
     2                + d3fAtdtdA2*dAdnb*dAdnb*dtdna
     3                + 2.0d0*d2fAtdA2*d2Adnadnb*dAdnb
     4                + d2fAtdA2*d2Adnb2*dAdna
     5                + 2.0d0*d3fAtdt2dA*dAdnb*dtdna*dtdnb
     6                + d3fAtdt2dA*dAdna*dtdnb*dtdnb
     7                + 2.0d0*d2fAtdtdA*(d2Adnadnb*dtdnb
     8                                 + dAdnb*d2tdnadnb)
     9                + d2fAtdtdA*(d2Adnb2*dtdna + dAdna*d2tdnb2)
     A                + dfAtdA*d3Adnadnb2 + d3fAtdt3*dtdna*dtdnb*dtdnb
     B                + dfAtdt*d3tdnadnb2
     C                + 2.0d0*d2fAtdt2*d2tdnadnb*dtdnb
     D                + d2fAtdt2*d2tdnb2*dtdna
         d3fAtdnb3    = d3fAtdA3*dAdnb*dAdnb*dAdnb
     1                + 3.0d0*d3fAtdtdA2*dAdnb*dAdnb*dtdnb
     2                + 3.0d0*d2fAtdA2*d2Adnb2*dAdnb
     3                + 3.0d0*d3fAtdt2dA*dAdnb*dtdnb*dtdnb
     4                + 3.0d0*d2fAtdtdA*(d2Adnb2*dtdnb + dAdnb*d2tdnb2)
     5                + dfAtdA*d3Adnb3 + d3fAtdt3*dtdnb*dtdnb*dtdnb
     6                + dfAtdt*d3tdnb3
     7                + 3.0d0*d2fAtdt2*d2tdnb2*dtdnb
#endif
c
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + BETA/GAMMA*t2*fAt
         darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
     &                            +t2*dfAtdna)
         darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
     &                            +t2*dfAtdnb)
c
         if (id_cpbe.eq.3) then
           darglogdna = darglogdna + dBETAdr/GAMMA*t**2*fAt
           darglogdnb = darglogdnb + dBETAdr/GAMMA*t**2*fAt
         endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         d2arglogdna2 = BETA/GAMMA*(2.0d0*dtdna*dtdna*fAt
     &                             +2.0d0*t*d2tdna2*fAt
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +2.0d0*t*dtdna*dfAtdna
     &                             +t2*d2fAtdna2)
         d2arglogdnb2 = BETA/GAMMA*(2.0d0*dtdnb*dtdnb*fAt
     &                             +2.0d0*t*d2tdnb2*fAt
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdnb
     &                             +t2*d2fAtdnb2)
         d2arglogdnadnb = BETA/GAMMA*(2.0d0*dtdna*dtdnb*fAt
     &                             +2.0d0*t*d2tdnadnb*fAt
     &                             +2.0d0*t*dtdna*dfAtdnb
     &                             +2.0d0*t*dtdnb*dfAtdna
     &                             +t2*d2fAtdnadnb)
c
         if (id_cpbe.eq.3) then
           d2arglogdna2 = d2arglogdna2
     &      + d2BETAdr2/GAMMA*t**2*fAt
     &      + dBETAdr/GAMMA*(2.0d0*t*dtdna*fAt+t*t*dfAtdna)
           d2arglogdnb2 = d2arglogdnb2
     &      + d2BETAdr2/GAMMA*t**2*fAt
     &      + dBETAdr/GAMMA*(2.0d0*t*dtdnb*fAt+t*t*dfAtdnb)
           d2arglogdnadnb = d2arglogdnadnb
     &      + d2BETAdr2/GAMMA*t**2*fAt
     &      + dBETAdr/GAMMA*(2.0d0*t*dtdna*fAt+t*t*dfAtdna)
     &      + dBETAdr/GAMMA*(2.0d0*t*dtdnb*fAt+t*t*dfAtdnb)
         endif
c
#endif
#ifdef THIRD_DERIV
         d3arglogdna3    = BETA/GAMMA*(6.0d0*dfAtdna*dtdna*dtdna
     1                                +6.0d0*dfAtdna*d2tdna2*t
     2                                +6.0d0*d2fAtdna2*dtdna*t
     3                                +6.0d0*fAt*dtdna*d2tdna2
     4                                +d3fAtdna3*t2
     5                                +2.0d0*fAt*d3tdna3*t)
         d3arglogdna2dnb = BETA/GAMMA*(4.0d0*dfAtdna*dtdna*dtdnb
     1                                +2.0d0*dfAtdnb*dtdna*dtdna
     2                                +4.0d0*dfAtdna*d2tdnadnb*t
     3                                +2.0d0*dfAtdnb*d2tdna2*t
     3                                +4.0d0*d2fAtdnadnb*dtdna*t
     4                                +2.0d0*d2fAtdna2*dtdnb*t
     5                                +4.0d0*fAt*dtdna*d2tdnadnb
     6                                +2.0d0*fAt*dtdnb*d2tdna2
     7                                +d3fAtdna2dnb*t2
     8                                +2.0d0*fAt*d3tdna2dnb*t)
         d3arglogdnadnb2 = BETA/GAMMA*(4.0d0*dfAtdnb*dtdna*dtdnb
     1                                +2.0d0*dfAtdna*dtdnb*dtdnb
     2                                +4.0d0*dfAtdnb*d2tdnadnb*t
     3                                +2.0d0*dfAtdna*d2tdnb2*t
     3                                +4.0d0*d2fAtdnadnb*dtdnb*t
     4                                +2.0d0*d2fAtdnb2*dtdna*t
     5                                +4.0d0*fAt*dtdnb*d2tdnadnb
     6                                +2.0d0*fAt*dtdna*d2tdnb2
     7                                +d3fAtdnadnb2*t2
     8                                +2.0d0*fAt*d3tdnadnb2*t)
         d3arglogdnb3    = BETA/GAMMA*(6.0d0*dfAtdnb*dtdnb*dtdnb
     1                                +6.0d0*dfAtdnb*d2tdnb2*t
     2                                +6.0d0*d2fAtdnb2*dtdnb*t
     3                                +6.0d0*fAt*dtdnb*d2tdnb2
     4                                +d3fAtdnb3*t2
     5                                +2.0d0*fAt*d3tdnb3*t)
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
         dMzdna = 3.0d0*phi2*dphidna
         dMzdnb = 3.0d0*phi2*dphidnb
c
         Nzrt = dlog(arglog)
         dNzrtdna = darglogdna/arglog
         dNzrtdnb = darglogdnb/arglog
c
         H0 = GAMMA*Mz*Nzrt
         dH0dna = GAMMA*(dMzdna*Nzrt + Mz*dNzrtdna)
         dH0dnb = GAMMA*(dMzdnb*Nzrt + Mz*dNzrtdnb)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         d2Mzdna2 = 3.0d0*phi2*d2phidna2 + 6.0d0*phi*dphidna*dphidna
         d2Mzdnadnb = 3.0d0*phi2*d2phidnadnb
     1              + 6.0d0*phi*dphidna*dphidnb
         d2Mzdnb2 = 3.0d0*phi2*d2phidnb2 + 6.0d0*phi*dphidnb*dphidnb
c
         d2Nzrtdna2 = d2arglogdna2/arglog -
     1                darglogdna*darglogdna/arglog/arglog
         d2Nzrtdnadnb = d2arglogdnadnb/arglog -
     1                  darglogdna*darglogdnb/arglog/arglog
         d2Nzrtdnb2 = d2arglogdnb2/arglog -
     1                darglogdnb*darglogdnb/arglog/arglog
c
         d2H0dna2 = GAMMA*( d2Mzdna2*Nzrt + 2.0d0*dMzdna*dNzrtdna
     &                    + Mz*d2Nzrtdna2)
         d2H0dnadnb = GAMMA*( d2Mzdnadnb*Nzrt + dMzdna*dNzrtdnb
     &                      + dMzdnb*dNzrtdna + Mz*d2Nzrtdnadnb)
         d2H0dnb2 = GAMMA*( d2Mzdnb2*Nzrt + 2.0d0*dMzdnb*dNzrtdnb
     &                    + Mz*d2Nzrtdnb2)         
#endif
#ifdef THIRD_DERIV
         d3Mzdna3 = 3.0d0*phi2*d3phidna3
     1            + 18.0d0*phi*d2phidna2*dphidna
     2            + 6.0d0*dphidna*dphidna*dphidna
         d3Mzdna2dnb = 3.0d0*phi2*d3phidna2dnb
     1               + 12.0d0*phi*d2phidnadnb*dphidna
     2               + 6.0d0*phi*d2phidna2*dphidnb
     3               + 6.0d0*dphidna*dphidna*dphidnb
         d3Mzdnadnb2 = 3.0d0*phi2*d3phidnadnb2
     1               + 12.0d0*phi*d2phidnadnb*dphidnb
     2               + 6.0d0*phi*d2phidnb2*dphidna
     3               + 6.0d0*dphidna*dphidnb*dphidnb
         d3Mzdnb3 = 3.0d0*phi2*d3phidnb3
     1            + 18.0d0*phi*d2phidnb2*dphidnb
     2            + 6.0d0*dphidnb*dphidnb*dphidnb
c
         d3Nzrtdna3 = d3arglogdna3/arglog
     1              - 3.0d0*d2arglogdna2*darglogdna/arglog/arglog
     2              + 2.0d0*darglogdna*darglogdna*darglogdna/
     3                arglog/arglog/arglog
         d3Nzrtdna2dnb = d3arglogdna2dnb/arglog
     1                 - 2.0d0*d2arglogdnadnb*darglogdna/arglog/arglog
     2                 - d2arglogdna2*darglogdnb/arglog/arglog
     3                 + 2.0d0*darglogdna*darglogdna*darglogdnb/
     4                   arglog/arglog/arglog
         d3Nzrtdnadnb2 = d3arglogdnadnb2/arglog
     1                 - 2.0d0*d2arglogdnadnb*darglogdnb/arglog/arglog
     2                 - d2arglogdnb2*darglogdna/arglog/arglog
     3                 + 2.0d0*darglogdna*darglogdnb*darglogdnb/
     4                   arglog/arglog/arglog
         d3Nzrtdnb3 = d3arglogdnb3/arglog
     1              - 3.0d0*d2arglogdnb2*darglogdnb/arglog/arglog
     2              + 2.0d0*darglogdnb*darglogdnb*darglogdnb/
     3                arglog/arglog/arglog
c
         d3H0dna3    = GAMMA*( d3Mzdna3*Nzrt
     1                       + 3.0d0*d2Mzdna2*dNzrtdna
     2                       + 3.0d0*dMzdna*d2Nzrtdna2
     3                       + Mz*d3Nzrtdna3)
         d3H0dna2dnb = GAMMA*( d3Mzdna2dnb*Nzrt
     1                       + 2.0d0*d2Mzdnadnb*dNzrtdna
     2                       + d2Mzdna2*dNzrtdnb
     3                       + 2.0d0*dMzdna*d2Nzrtdnadnb
     4                       + dMzdnb*d2Nzrtdna2
     5                       + Mz*d3Nzrtdna2dnb)
         d3H0dnadnb2 = GAMMA*( d3Mzdnadnb2*Nzrt
     1                       + 2.0d0*d2Mzdnadnb*dNzrtdnb
     2                       + d2Mzdnb2*dNzrtdna
     3                       + 2.0d0*dMzdnb*d2Nzrtdnadnb
     4                       + dMzdna*d2Nzrtdnb2
     5                       + Mz*d3Nzrtdnadnb2)
         d3H0dnb3    = GAMMA*( d3Mzdnb3*Nzrt
     1                       + 3.0d0*d2Mzdnb2*dNzrtdnb
     2                       + 3.0d0*dMzdnb*d2Nzrtdnb2
     3                       + Mz*d3Nzrtdnb3)
#endif
c
c        Now we update Ec, Amat, and Amat2
c
c Daniel (3-5-13): I believe that lfac is set to false for all 
c calculations that I do.
         if (lfac) then
            Ec = Ec+nepsc*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+nepsc*fac               
         endif
         if (nlfac) then
            Ec = Ec+(H0*rhoval)*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+(H0*rhoval)*fac
         endif
         if (lfac) then
            Amat(n,1) = Amat(n,1) + dnepscdn(1)*fac
            if (ipol.eq.2) Amat(n,2) = Amat(n,2) + dnepscdn(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) 
     &                        + d2nepscdn2(D2_RA_RA)*fac
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) 
     &                        + d2nepscdn2(D2_RA_RB)*fac
            if (ipol.eq.2)
     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) 
     &                        + d2nepscdn2(D2_RB_RB)*fac
#endif
c Daniel (7-31-12): XC third derivatives (LDA part)
#ifdef THIRD_DERIV
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1                           + d3nepscdn3(D3_RA_RA_RA)*fac        
            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)         
     1                           + d3nepscdn3(D3_RA_RA_RB)*fac
            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
     1                           + d3nepscdn3(D3_RA_RB_RB)*fac
            if (ipol.eq.2)
     1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)       
     2                           + d3nepscdn3(D3_RB_RB_RB)*fac  
#endif
         endif
         if (nlfac)then
            Amat(n,1) = Amat(n,1) +  (H0 +
     &                  rhoval*dH0dna)*fac
            if (ipol.eq.2) Amat(n,2) = Amat(n,2) +  (H0 +
     &                  rhoval*dH0dnb)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &                + (2.d0*dH0dna + rhoval*d2H0dna2)*fac
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &                + (dH0dna + dH0dnb + rhoval*d2H0dnadnb)*fac
            if (ipol.eq.2)
     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &                + (2.d0*dH0dnb + rhoval*d2H0dnb2)*fac
#endif
c Daniel (7-31-12): XC third derivatives (GGA part)
#ifdef THIRD_DERIV
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1        + (3.0d0*d2H0dna2 + rhoval*d3H0dna3)*fac
            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)
     1        + (d2H0dna2 + 2.0d0*d2H0dnadnb + rhoval*d3H0dna2dnb)*fac
            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
     1        + (d2H0dnb2 + 2.0d0*d2H0dnadnb + rhoval*d3H0dnadnb2)*fac
            if (ipol.eq.2)
     1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     2        + (3.0d0*d2H0dnb2 + rhoval*d3H0dnb3)*fac
#endif
         endif
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
         if (dsqgamma.gt.TOLL)then
            dtdg = 0.25d0*rphi*rks*rrho/dsqgamma
            dfAtdg = dfAtdt*dtdg
            darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t2*dfAtdg)
c
            dNzrtdg = darglogdg/arglog
c
            dH0dg = GAMMA*Mz*dNzrtdg
            if (ipol.eq.1) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*fac
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*fac*2.0d0
            else
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dH0dg*rhoval*fac
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dH0dg*rhoval*fac*2.0d0
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dH0dg*rhoval*fac
            endif
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2tdg2 = -0.125d0*rphi*rks*rrho/
     &               (dsqgamma*dsqgamma*dsqgamma)
            d2tdnadg = -dtdg*rrho-dtdg*rphi*dphidna
     &                 -dtdg*rks*dksdna
            d2tdnbdg = -dtdg*rrho-dtdg*rphi*dphidnb
     &                 -dtdg*rks*dksdnb
            d2fAtdg2 = d2fAtdt2*dtdg*dtdg+dfAtdt*d2tdg2
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
     &                                 +t2*d2fAtdnadg)
            d2arglogdnbdg = BETA/GAMMA*(2.0d0*dtdnb*dtdg*fAt
     &                                 +2.0d0*t*d2tdnbdg*fAt
     &                                 +2.0d0*t*dtdg*dfAtdnb
     &                                 +2.0d0*t*dtdnb*dfAtdg
     &                                 +t2*d2fAtdnbdg)
c
            if (id_cpbe.eq.3) then
              d2arglogdnadg = d2arglogdnadg
     &          + dBETAdr/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
              d2arglogdnbdg = d2arglogdnbdg
     &          + dBETAdr/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
            endif
c
            d2arglogdg2 = BETA/GAMMA*(2.0d0*dtdg*dtdg*fAt
     &                               +2.0d0*t*d2tdg2*fAt
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +2.0d0*t*dtdg*dfAtdg
     &                               +t2*d2fAtdg2)
c
            d2Nzrtdg2 = d2arglogdg2/arglog -
     1                  darglogdg*darglogdg/arglog/arglog
            d2Nzrtdnadg = d2arglogdnadg/arglog -
     1                    darglogdna*darglogdg/arglog/arglog
            d2Nzrtdnbdg = d2arglogdnbdg/arglog -
     1                    darglogdnb*darglogdg/arglog/arglog
c
            d2H0dnadg = GAMMA*dMzdna*dNzrtdg
     &                + GAMMA*Mz*d2Nzrtdnadg
            d2H0dnbdg = GAMMA*dMzdnb*dNzrtdg
     &                + GAMMA*Mz*d2Nzrtdnbdg
            d2H0dg2 = GAMMA*Mz*d2Nzrtdg2
c
c ORIGINAL CODE
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c     &             + (dH0dg + d2H0dnadg*rhoval)*fac
c            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
c     &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*fac
c            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
c     &             + (dH0dg + d2H0dnadg*rhoval)*fac
c            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
c     &             + (dH0dg + d2H0dnbdg*rhoval)*fac
c            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
c     &             + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*fac
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c     &             + (dH0dg + d2H0dnbdg*rhoval)*fac
c            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c     &             + d2H0dg2*rhoval*fac
c            Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
c     &             + 2.0d0*d2H0dg2*rhoval*fac
c            Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
c     &             + d2H0dg2*rhoval*fac
c            Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
c     &             + 4.0d0*d2H0dg2*rhoval*fac
c            Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
c     &             + 2.0d0*d2H0dg2*rhoval*fac
c            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c     &             + d2H0dg2*rhoval*fac
c ORIGINAL CODE
c
c Daniel (3-11-13): Only do beta terms if they are needed
c since the extra calculations are unnecessary for restricted
c calculations.
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &             + (dH0dg + d2H0dnadg*rhoval)*fac
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB)
     &             + 2.0d0*(dH0dg + d2H0dnadg*rhoval)*fac
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB)
     &             + (dH0dg + d2H0dnadg*rhoval)*fac
c
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &             + d2H0dg2*rhoval*fac
            Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB)
     &             + 2.0d0*d2H0dg2*rhoval*fac
            Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB)
     &             + d2H0dg2*rhoval*fac
            Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB)
     &             + 4.0d0*d2H0dg2*rhoval*fac
            if (ipol.eq.2) then
              Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA)
     &               + (dH0dg + d2H0dnbdg*rhoval)*fac
              Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB)
     &               + 2.0d0*(dH0dg + d2H0dnbdg*rhoval)*fac
              Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &               + (dH0dg + d2H0dnbdg*rhoval)*fac
c
              Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB)
     &               + 2.0d0*d2H0dg2*rhoval*fac
              Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &               + d2H0dg2*rhoval*fac
            endif
#endif
#ifdef THIRD_DERIV
c Derivatives of t, the reduced density gradient
            d3tdg3 = 0.1875d0*rphi*rks*rrho/dsqgamma/dsqgamma/
     1               dsqgamma/dsqgamma/dsqgamma
            d3tdna2dg = dtdg*( dphidna*dphidna*rphi2
     1                       - d2phidna2*rphi
     2                       + dksdna*dksdna*rks2
     3                       - d2ksdna2*rks
     4                       + rrho2)
     5                - d2tdnadg*( dphidna*rphi
     6                           + dksdna*rks
     7                           + rrho)
            d3tdnadnbdg = dtdg*( dphidna*dphidnb*rphi2
     1                       - d2phidnadnb*rphi
     2                       + dksdna*dksdnb*rks2
     3                       - d2ksdnadnb*rks
     4                       + rrho2)
     5                  - d2tdnbdg*( dphidna*rphi
     6                             + dksdna*rks
     7                             + rrho)
            d3tdnbdnadg = dtdg*( dphidna*dphidnb*rphi2
     1                       - d2phidnadnb*rphi
     2                       + dksdna*dksdnb*rks2
     3                       - d2ksdnadnb*rks
     4                       + rrho2)
     5                  - d2tdnadg*( dphidnb*rphi
     6                             + dksdnb*rks
     7                             + rrho)
            d3tdnb2dg = dtdg*( dphidnb*dphidnb*rphi2
     1                       - d2phidnb2*rphi
     2                       + dksdnb*dksdnb*rks2
     3                       - d2ksdnb2*rks
     4                       + rrho2)
     5                - d2tdnbdg*( dphidnb*rphi
     6                           + dksdnb*rks
     7                           + rrho)
            d3tdnadg2 = -d2tdg2*( dphidna*rphi + dksdna*rks
     1                          + rrho)
            d3tdnbdg2 = -d2tdg2*( dphidnb*rphi + dksdnb*rks
     1                          + rrho)
c Derivatives of the rational function fAt
            d3fAtdg3      = d3fAtdt3*dtdg*dtdg*dtdg
     1                    + 3.0d0*d2fAtdt2*d2tdg2*dtdg
     2                    + dfAtdt*d3tdg3
            d3fAtdna2dg   = d3fAtdtdA2*dAdna*dAdna*dtdg
     1                    + 2.0d0*d3fAtdt2dA*dAdna*dtdna*dtdg
     2                    + d2fAtdtdA*( d2Adna2*dtdg
     3                                + 2.0d0*dAdna*d2tdnadg)
     4                    + d3fAtdt3*dtdna*dtdna*dtdg
     5                    + d2fAtdt2*( d2tdna2*dtdg
     6                               + 2.0d0*d2tdnadg*dtdna)
     7                    + dfAtdt*d3tdna2dg
            d3fAtdnadnbdg = d3fAtdtdA2*dAdna*dAdnb*dtdg
     1                    + d3fAtdt2dA*( dAdna*dtdnb
     2                                 + dAdnb*dtdna)*dtdg
     3                    + d2fAtdtdA*( d2Adnadnb*dtdg
     4                                + dAdna*d2tdnbdg
     5                                + dAdnb*d2tdnadg)
     6                    + d3fAtdt3*dtdna*dtdnb*dtdg
     7                    + d2fAtdt2*( d2tdnadnb*dtdg
     8                               + d2tdnadg*dtdnb
     9                               + d2tdnbdg*dtdna)
     A                    + dfAtdt*d3tdnadnbdg
            d3fAtdnbdnadg = d3fAtdnadnbdg
            d3fAtdnb2dg   = d3fAtdtdA2*dAdnb*dAdnb*dtdg
     1                    + 2.0d0*d3fAtdt2dA*dAdnb*dtdnb*dtdg
     2                    + d2fAtdtdA*( d2Adnb2*dtdg
     3                                + 2.0d0*dAdnb*d2tdnbdg)
     4                    + d3fAtdt3*dtdnb*dtdnb*dtdg
     5                    + d2fAtdt2*( d2tdnb2*dtdg
     6                               + 2.0d0*d2tdnbdg*dtdnb)
     7                    + dfAtdt*d3tdnb2dg
            d3fAtdnadg2   = d3fAtdt2dA*dAdna*dtdg*dtdg
     1                    + d2fAtdtdA*dAdna*d2tdg2
     2                    + d3fAtdt3*dtdna*dtdg*dtdg
     3                    + d2fAtdt2*( 2.0d0*d2tdnadg*dtdg
     4                               + dtdna*d2tdg2)
     5                    + dfAtdt*d3tdnadg2
            d3fAtdnbdg2   = d3fAtdt2dA*dAdnb*dtdg*dtdg
     1                    + d2fAtdtdA*dAdnb*d2tdg2
     2                    + d3fAtdt3*dtdnb*dtdg*dtdg
     3                    + d2fAtdt2*( 2.0d0*d2tdnbdg*dtdg
     4                               + dtdnb*d2tdg2)
     5                    + dfAtdt*d3tdnbdg2
c Derivatives of the logarithm argument arglog
            d3arglogdg3      = BETA/GAMMA*(d3fAtdg3*t2
     1                       + 6.0d0*d2fAtdg2*dtdg*t
     2                       + 6.0d0*dfAtdg*( dtdg*dtdg
     3                                      + d2tdg2*t)
     4                       + 6.0d0*fAt*dtdg*d2tdg2
     5                       + 2.0d0*fAt*d3tdg3*t)
            d3arglogdna2dg   = 2.0d0*BETA/GAMMA*(
     1                         2.0d0*d2fAtdnadg*dtdna*t
     2                       + 2.0d0*dfAtdna*d2tdnadg*t
     3                       + 2.0d0*dfAtdna*dtdna*dtdg
     4                       + dfAtdg*dtdna*dtdna
     5                       + 2.0d0*fAt*d2tdnadg*dtdna
     6                       + d2fAtdna2*dtdg*t
     7                       + dfAtdg*d2tdna2*t
     8                       + fAt*d3tdna2dg*t
     9                       + fAt*d2tdna2*dtdg
     A                       + 0.50d0*d3fAtdna2dg*t2)
            d3arglogdnadnbdg = 2.0d0*BETA/GAMMA*(
     1                         d2fAtdnadg*dtdnb*t
     2                       + d2fAtdnbdg*dtdna*t
     3                       + dfAtdna*d2tdnbdg*t
     4                       + dfAtdnb*d2tdnadg*t
     5                       + dfAtdna*dtdnb*dtdg
     6                       + dfAtdnb*dtdna*dtdg
     7                       + dfAtdg*dtdna*dtdnb
     8                       + fAt*d2tdnadg*dtdnb
     9                       + fAt*d2tdnbdg*dtdna
     A                       + d2fAtdnadnb*dtdg*t
     B                       + dfAtdg*d2tdnadnb*t
     C                       + fAt*d3tdnadnbdg*t
     D                       + fAt*d2tdnadnb*dtdg
     E                       + 0.50d0*d3fAtdnadnbdg*t2)
            d3arglogdnbdnadg = 2.0d0*BETA/GAMMA*(
     1                         d2fAtdnbdg*dtdna*t
     2                       + d2fAtdnadg*dtdnb*t
     3                       + dfAtdnb*d2tdnadg*t
     4                       + dfAtdna*d2tdnbdg*t
     5                       + dfAtdnb*dtdna*dtdg
     6                       + dfAtdna*dtdnb*dtdg
     7                       + dfAtdg*dtdnb*dtdna
     8                       + fAt*d2tdnbdg*dtdna
     9                       + fAt*d2tdnadg*dtdnb
     A                       + d2fAtdnadnb*dtdg*t
     B                       + dfAtdg*d2tdnadnb*t
     C                       + fAt*d3tdnbdnadg*t
     D                       + fAt*d2tdnadnb*dtdg
     E                       + 0.50d0*d3fAtdnbdnadg*t2)
            d3arglogdnb2dg   = 2.0d0*BETA/GAMMA*(
     1                         2.0d0*d2fAtdnbdg*dtdnb*t
     2                       + 2.0d0*dfAtdnb*d2tdnbdg*t
     3                       + 2.0d0*dfAtdnb*dtdnb*dtdg
     4                       + dfAtdg*dtdnb*dtdnb
     5                       + 2.0d0*fAt*d2tdnbdg*dtdnb
     6                       + d2fAtdnb2*dtdg*t
     7                       + dfAtdg*d2tdnb2*t
     8                       + fAt*d3tdnb2dg*t
     9                       + fAt*d2tdnb2*dtdg
     A                       + 0.50d0*d3fAtdnb2dg*t2)
            d3arglogdnadg2   = 2.0d0*BETA/GAMMA*(
     1                         2.0d0*d2fAtdnadg*dtdg*t
     2                       + dfAtdna*d2tdg2*t
     3                       + dfAtdna*dtdg*dtdg
     4                       + d2fAtdg2*dtdna*t
     5                       + 2.0d0*dfAtdg*d2tdnadg*t
     6                       + 2.0d0*dfAtdg*dtdna*dtdg
     7                       + 2.0d0*fAt*d2tdnadg*dtdg
     8                       + fAt*dtdna*d2tdg2
     9                       + fAt*d3tdnadg2*t
     A                       + 0.50d0*d3fAtdnadg2*t2)
            d3arglogdnbdg2   = 2.0d0*BETA/GAMMA*(
     1                         2.0d0*d2fAtdnbdg*dtdg*t
     2                       + dfAtdnb*d2tdg2*t
     3                       + dfAtdnb*dtdg*dtdg
     4                       + d2fAtdg2*dtdnb*t
     5                       + 2.0d0*dfAtdg*d2tdnbdg*t
     6                       + 2.0d0*dfAtdg*dtdnb*dtdg
     7                       + 2.0d0*fAt*d2tdnbdg*dtdg
     8                       + fAt*dtdnb*d2tdg2
     9                       + fAt*d3tdnbdg2*t
     A                       + 0.50d0*d3fAtdnbdg2*t2)
c Derivatives of the correlation enhancement factor
            d3Nzrtdg3 = d3arglogdg3/arglog
     1                - 3.0d0*d2arglogdg2*darglogdg/arglog/arglog
     2                + 2.0d0*darglogdg*darglogdg*darglogdg/arglog/
     3                  arglog/arglog
            d3Nzrtdnadg2 = d3arglogdnadg2/arglog
     1                   - (2.0d0*d2arglogdnadg*darglogdg +
     2                      darglogdna*d2arglogdg2)/arglog/arglog
     3                   + 2.0d0*darglogdna*darglogdg*darglogdg/
     4                     arglog/arglog/arglog
            d3Nzrtdnbdg2 = d3arglogdnbdg2/arglog
     1                   - (2.0d0*d2arglogdnbdg*darglogdg +
     2                      darglogdnb*d2arglogdg2)/arglog/arglog
     3                   + 2.0d0*darglogdnb*darglogdg*darglogdg/
     4                     arglog/arglog/arglog
            d3Nzrtdna2dg = d3arglogdna2dg/arglog
     1                   - (2.0d0*d2arglogdnadg*darglogdna +
     2                      d2arglogdna2*darglogdg)/arglog/arglog
     3                   + 2.0d0*darglogdna*darglogdna*darglogdg/
     4                     arglog/arglog/arglog
            d3Nzrtdnadnbdg = d3arglogdnadnbdg/arglog
     1                     - (d2arglogdnadg*darglogdnb +
     2                        d2arglogdnbdg*darglogdna +
     3                        d2arglogdnadnb*darglogdg)/arglog/arglog
     4                     + 2.0d0*darglogdna*darglogdnb*darglogdg/
     5                       arglog/arglog/arglog
            d3Nzrtdnb2dg = d3arglogdnb2dg/arglog
     1                   - (2.0d0*d2arglogdnbdg*darglogdnb +
     2                      d2arglogdnb2*darglogdg)/arglog/arglog
     3                   + 2.0d0*darglogdnb*darglogdnb*darglogdg/
     4                     arglog/arglog/arglog
c
            d3H0dg3      = GAMMA*Mz*d3Nzrtdg3
            d3H0dna2dg   = GAMMA*( d2Mzdna2*dNzrtdg
     1                           + 2.0d0*dMzdna*d2Nzrtdnadg
     2                           + Mz*d3Nzrtdna2dg)
            d3H0dnadnbdg = GAMMA*( d2Mzdnadnb*dNzrtdg 
     1                           + dMzdna*d2Nzrtdnbdg
     2                           + dMzdnb*d2Nzrtdnadg
     3                           + Mz*d3Nzrtdnadnbdg)
            d3H0dnbdnadg = d3H0dnadnbdg
            d3H0dnb2dg   = GAMMA*( d2Mzdnb2*dNzrtdg
     1                           + 2.0d0*dMzdnb*d2Nzrtdnbdg
     2                           + Mz*d3Nzrtdnb2dg)
            d3H0dnadg2   = GAMMA*( dMzdna*d2Nzrtdg2
     1                           + Mz*d3Nzrtdnadg2)
            d3H0dnbdg2   = GAMMA*( dMzdnb*d2Nzrtdg2
     1                           + Mz*d3Nzrtdnbdg2)
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
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1             + (2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*fac
            Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB)
     1             + 2.0d0*(2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*fac
            Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB)
     1             + (2.0d0*d2H0dnadg + d3H0dna2dg*rhoval)*fac
c Mixed derivatives dradrbdg
            Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA)
     1             + ( d2H0dnadg + d2H0dnbdg 
     2               + d3H0dnadnbdg*rhoval)*fac
            Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB)
     1             + 2.0d0*( d2H0dnadg + d2H0dnbdg 
     2                     + d3H0dnadnbdg*rhoval)*fac
            Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB)
     1             + ( d2H0dnadg + d2H0dnbdg 
     2               + d3H0dnadnbdg*rhoval)*fac
c Mixed derivatives drbdrbdg
            if (ipol.eq.2) then
              Cmat3(n,D3_RB_RB_GAA) = Cmat3(n,D3_RB_RB_GAA)
     1               + (2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*fac
              Cmat3(n,D3_RB_RB_GAB) = Cmat3(n,D3_RB_RB_GAB)
     1               + 2.0d0*(2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*fac
              Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
     1               + (2.0d0*d2H0dnbdg + d3H0dnb2dg*rhoval)*fac
            endif
c Mixed derivatives dradgdg
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
     1             + (d2H0dg2 + d3H0dnadg2*rhoval)*fac
            Cmat3(n,D3_RA_GAA_GAB) = Cmat3(n,D3_RA_GAA_GAB)
     1             + 2.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*fac
            Cmat3(n,D3_RA_GAA_GBB) = Cmat3(n,D3_RA_GAA_GBB)
     1             + (d2H0dg2 + d3H0dnadg2*rhoval)*fac
            Cmat3(n,D3_RA_GAB_GAB) = Cmat3(n,D3_RA_GAB_GAB)
     1             + 4.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*fac
            Cmat3(n,D3_RA_GAB_GBB) = Cmat3(n,D3_RA_GAB_GBB)
     1             + 2.0d0*(d2H0dg2 + d3H0dnadg2*rhoval)*fac
            Cmat3(n,D3_RA_GBB_GBB) = Cmat3(n,D3_RA_GBB_GBB)
     1             + (d2H0dg2 + d3H0dnadg2*rhoval)*fac
c Mixed derivatives drbdgdg
            if (ipol.eq.2) then
              Cmat3(n,D3_RB_GAA_GAA) = Cmat3(n,D3_RB_GAA_GAA)
     1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*fac
              Cmat3(n,D3_RB_GAA_GAB) = Cmat3(n,D3_RB_GAA_GAB)
     1               + 2.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*fac
              Cmat3(n,D3_RB_GAA_GBB) = Cmat3(n,D3_RB_GAA_GBB)
     1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*fac
              Cmat3(n,D3_RB_GAB_GAB) = Cmat3(n,D3_RB_GAB_GAB)
     1               + 4.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*fac
              Cmat3(n,D3_RB_GAB_GBB) = Cmat3(n,D3_RB_GAB_GBB)
     1               + 2.0d0*(d2H0dg2 + d3H0dnbdg2*rhoval)*fac
              Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB)
     1               + (d2H0dg2 + d3H0dnbdg2*rhoval)*fac
            endif
c Derivatives dgaadgdg
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
     1             + (d3H0dg3*rhoval)*fac
            Cmat3(n,D3_GAA_GAA_GAB) = Cmat3(n,D3_GAA_GAA_GAB)
     1             + 2.0d0*(d3H0dg3*rhoval)*fac
            Cmat3(n,D3_GAA_GAA_GBB) = Cmat3(n,D3_GAA_GAA_GBB)
     1             + (d3H0dg3*rhoval)*fac
            Cmat3(n,D3_GAA_GAB_GAB) = Cmat3(n,D3_GAA_GAB_GAB)
     1             + 4.0d0*(d3H0dg3*rhoval)*fac
            Cmat3(n,D3_GAA_GAB_GBB) = Cmat3(n,D3_GAA_GAB_GBB)
     1             + 2.0d0*(d3H0dg3*rhoval)*fac
            Cmat3(n,D3_GAA_GBB_GBB) = Cmat3(n,D3_GAA_GBB_GBB)
     1             + (d3H0dg3*rhoval)*fac
c Derivatives dgabdgdg
            Cmat3(n,D3_GAB_GAB_GAB) = Cmat3(n,D3_GAB_GAB_GAB)
     1             + 8.0d0*(d3H0dg3*rhoval)*fac
            if (ipol.eq.2) then
              Cmat3(n,D3_GAB_GAB_GBB) = Cmat3(n,D3_GAB_GAB_GBB)
     1               + 4.0d0*(d3H0dg3*rhoval)*fac
              Cmat3(n,D3_GAB_GBB_GBB) = Cmat3(n,D3_GAB_GBB_GBB)
     1               + 2.0d0*(d3H0dg3*rhoval)*fac
c Derivatives dgbbdgdg
              Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
     1               + (d3H0dg3*rhoval)*fac
            endif
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
#include "xc_pbe96.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_pbe96.F"
#endif
c $Id$
