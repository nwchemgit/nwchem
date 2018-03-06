c note that cfac, lcfac, nlcfac are single numbers
c in the original cpbe96 file, they are arrays
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Utility routine to evaluate the TPSS and PKZB functionals
C>
      Subroutine nwpwxc_c_Mpbe(tol_rho, ipol, rho, rgamma,
     &                       Ec, Amat, Cmat)
c
c$Id$
c
      Implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
c
c     Charge Density 
c
      double precision rho(ipol) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(ipol+1) !< [Input] The norm of the density
                                      !< gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision Ec      !< [Output] The value of the functional
      double precision Amat(2) !< [Output] The derivative wrt rho
      double precision Cmat(3) !< [Output] The derivative wrt rgamma

      double precision dummy(1)
c
c     Constants in PBE functional
c
      double precision GAMMA, BETA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETA = 0.06672455060314922d0)
      parameter (PI = 3.1415926535897932385d0)
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
      double precision fac
      parameter (fac = 1.0d0)
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

c      write(0,*) 'upon arrival in cpbe  Ec=',Ec
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
c        n and zeta = (na - nb)/n
c
         rhoval = 0.0d0
         if (ipol.eq.1) then
            rhoval = rho(R_T)
            rho_t(1) = rho(R_T)
         else if (ipol.eq.2) then
            rho_t(1) = rho(R_A)
            rho_t(2) = rho(R_B)
            rhoval = rho_t(1)+rho_t(2)
         endif
         if (rhoval.le.tol_rho) goto 20
         if (ipol.eq.1) then
c           gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1)
            gammaval = rgamma(G_TT)
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
            gammaval = rgamma(G_AA)+rgamma(G_BB)+2*rgamma(G_AB)
         endif
         dsqgamma = max(dsqrt(gammaval),tol_rho)
         nepsc = 0.0d0
         dnepscdn(1) = 0.0d0
         if (ipol.eq.2) dnepscdn(2) = 0.0d0
c
c        call for LDA bit
c        this implementation temporarily assigns the pw91LDA for
c        use in the metaGGA local part
c
         call nwpwxc_c_pw91lda(tol_rho,ipol,1,1.0d0,rho_t,
     >                         nepsc,dnepscdn)
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
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = depscdna
         else
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = dnepscdn(2)/rhoval-nepsc/(rhoval**2)
         endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
         dksdna = (1.0d0/6.0d0)*ks/rhoval
         dksdnb = dksdna
c
c        zeta = (na-nb)/n and its derivs
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(R_A)-rho(R_B))/rhoval
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
         if (ipol.eq.1) then
            dzetadna = 1.0d0/rhoval
            dzetadnb = -1.0d0/rhoval
         else
            dzetadna =  2.0d0*rho(R_B)/(rhoval**2)
            dzetadnb = -2.0d0*rho(R_A)/(rhoval**2)
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
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
         t = dsqgamma/(2.0d0*phi*ks*rhoval)
         dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
         dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
c
c        { ... } in A (see below) and its derivs
c
         argexp = -epsc/GAMMA/(phi**3)
         dargexpdna = -depscdna/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidna
         dargexpdnb = -depscdnb/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidnb
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
c
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
         arglog = 1.0d0 + BETA/GAMMA*t**2*fAt
         darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
     &                            +t*t*dfAtdna)
         darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
     &                            +t*t*dfAtdnb)
c
c        H0 = GAMMA * phi**3 * log{arglog} and its derivs
c
         H0 = GAMMA*(phi**3)*dlog(arglog)
         dH0dna = GAMMA*(3.0d0*(phi**2)*dphidna*dlog(arglog)
     &                  +(phi**3)*darglogdna/arglog)
         dH0dnb = GAMMA*(3.0d0*(phi**2)*dphidnb*dlog(arglog)
     &                  +(phi**3)*darglogdnb/arglog)
c
c        Now we update Ec, Amat, and Amat2
c

c          NOTE:  this PBE does the LDA part of Ec in house
            Ec = Ec+epsc*fac
            Ec = Ec+H0*fac
            Amat(D1_RA) = Amat(D1_RA) + depscdna
            if (ipol.eq.2) Amat(D1_RB) = Amat(D1_RB) + depscdnb

            Amat(D1_RA) = Amat(D1_RA) +  dH0dna
            if (ipol.eq.2) Amat(D1_RB) = Amat(D1_RB) + dH0dnb
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
         if (dsqgamma.gt.TOLL)then
            dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
            dfAtdg = dfAtdt*dtdg
            darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
            dH0dg = GAMMA*(phi**3)*darglogdg/arglog

            Cmat(D1_GAA) = Cmat(D1_GAA) + dH0dg
            Cmat(D1_GAB) = Cmat(D1_GAB) + dH0dg*2.0d0
            if (ipol.eq.2) Cmat(D1_GBB) = Cmat(D1_GBB) + dH0dg
         endif
   20 continue
c

      return
      end
c
C> @}
