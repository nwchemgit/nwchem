#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_Mpbe.F
C> The Perdew, Kurth, Zupan, and Blaha modified PBE correlation
C> functional
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Modified version of the PBE correlation functional
C> to evaluate the TPSS and PKZB functionals
C>
C> The modified PBE correlation functional [1,2] that is used
C> in the TPSS and PKZB functionals. Although this functional is
C> not intended to be used stand-alone it is provided as a separate
C> subroutine like any other functional so that we can separate it
C> out in the symbolic algebra code generation.
C>
C> ### References ###
C>
C> [1] J.P. Perdew, S. Kurth, A. Zupan, P. Blaha,
C>     <i>"Accurate density functional with correct formal properties:
C>     A Step Beyond the Generalized Gradient Approximation"</i>,
C>     Phys. Rev. Lett. <b>82</b>, 2544-2547 (1999), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevLett.82.2544">
C>     10.1103/PhysRevLett.82.2544</a>.
C>
C> [2] J.P. Perdew, S. Kurth, A. Zupan, P. Blaha,
C>     <i>"Erratum: Accurate density functional with correct formal
C>     properties: A Step Beyond the Generalized Gradient
C>     Approximation"</i>,
C>     Phys. Rev. Lett. <b>82</b>, 5179-5179 (1999), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevLett.82.5179">
C>     10.1103/PhysRevLett.82.5179</a>.
C>
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_Mpbe_p(param,tol_rho,ipol,nq,wght,rho,rgamma,
     +                         Ec)
#else
      Subroutine nwxc_c_Mpbe(param,tol_rho,ipol,nq,wght,rho,rgamma,
     +                       Ec)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_Mpbe_d2(param,tol_rho,ipol,nq,wght,rho,rgamma,
     +                          Ec)
#else
      Subroutine nwxc_c_Mpbe_d3(param,tol_rho,ipol,nq,wght,rho,rgamma,
     +                          Ec)
#endif
c
c$Id$
c
#include "nwad.fh"
c
      Implicit none
c
#include "intf_nwxc_c_pw91lda.fh"
c
#include "nwxc_param.fh"
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
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the
                                    !< density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::Ec(nq)  !< [Output] The value of the functional
c     double precision Amat(2) !< [Output] The derivative wrt rho
c     double precision Cmat(3) !< [Output] The derivative wrt rgamma

      double precision dummy(1)
c
c     Constants in PBE functional
c
      double precision GAMMA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
c     parameter (BETA = 0.06672455060314922d0)
      parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
      type(nwad_dble)::rho_t(1,2),rgamma_t(1,3),tau_t(1,2)
c     type(nwad_dble)::dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n, iq
      type(nwad_dble)::rhoval, gammaval
      type(nwad_dble)::nepsc(1)
      double precision dnepscdn(2)
      type(nwad_dble)::epsc, depscdna, depscdnb
      type(nwad_dble)::H0, dH0dna, dH0dnb, dH0dg
      type(nwad_dble)::phi, dphidna, dphidnb, dphidzeta
      type(nwad_dble)::zeta, dzetadna, dzetadnb
      type(nwad_dble)::arglog, darglogdna, darglogdnb, darglogdg
      type(nwad_dble)::fAt, dfAtdt, dfAtdA
      type(nwad_dble)::fAtnum, dfAtnumdt, dfAtnumdA
      type(nwad_dble)::fAtden, dfAtdendt, dfAtdendA
      type(nwad_dble)::dfAtdna, dfAtdnb, dfAtdg
      type(nwad_dble)::A, dAdna, dAdnb
c     type(nwad_dble)::t, dtdna, dtdnb, dtdg
      type(nwad_dble)::t2, dtdna, dtdnb, dtdg
      type(nwad_dble)::ks, dksdna, dksdnb
      type(nwad_dble)::argexp, dargexpdna, dargexpdnb
      type(nwad_dble)::expinA
c     double precision fac
c     parameter (fac = 1.0d0)
c     
c      write(0,*) 'upon arrival in cpbe  Ec=',Ec
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
c        n and zeta = (na - nb)/n
c
      beta = param(1)
      do iq = 1, nq
         rhoval = 0.0d0
         if (ipol.eq.1) then
            rhoval = rho(iq,R_T)
            rho_t(1,1) = rho(iq,R_T)
            rho_t(1,2) = 0.0d0
         else if (ipol.eq.2) then
            rho_t(1,1) = rho(iq,R_A)
            rho_t(1,2) = rho(iq,R_B)
            rhoval = rho_t(1,1)+rho_t(1,2)
         endif
         if (rhoval.le.tol_rho) goto 20
         if (ipol.eq.1) then
c           gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1)
            gammaval = rgamma(iq,G_TT)
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
            gammaval = rgamma(iq,G_AA)+rgamma(iq,G_BB)
     +               + 2.0d0*rgamma(iq,G_AB)
         endif
c        dsqgamma = sqrt(gammaval)
         nepsc(1) = 0.0d0
c        dnepscdn(1) = 0.0d0
c        if (ipol.eq.2) dnepscdn(2) = 0.0d0
c
c        call for LDA bit
c        this implementation temporarily assigns the pw91LDA for
c        use in the metaGGA local part
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_pw91lda",func_lda,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,nepsc)
         call nwxc_c_pw91lda_p(tol_rho,ipol,1,1.0d0,rho_t,nepsc)
         call leave_subr("nwxc_c_pw91lda",func_lda,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,nepsc)
#else
         call nwxc_c_pw91lda(tol_rho,ipol,1,1.0d0,rho_t,nepsc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_pw91lda_d2(tol_rho,ipol,1,1.0d0,rho_t,nepsc)
#else
         call nwxc_c_pw91lda_d3(tol_rho,ipol,1,1.0d0,rho_t,nepsc)
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
c           depscdna = dnepscdn(1)/rhoval-nepsc(1)/(rhoval**2.0d0)
c           depscdnb = depscdna
c        else
c           depscdna = dnepscdn(1)/rhoval-nepsc(1)/(rhoval**2.0d0)
c           depscdnb = dnepscdn(2)/rhoval-nepsc(1)/(rhoval**2.0d0)
c        endif
c
c        ks = 2*(3*PI**2*n)**(1/6)/sqrt(PI) and its derivs
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
c        dksdna = (1.0d0/6.0d0)*ks/rhoval
c        dksdnb = dksdna
c
c        zeta = (na-nb)/n and its derivs
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(iq,R_A)-rho(iq,R_B))/rhoval
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
c        if (ipol.eq.1) then
c           dzetadna = 1.0d0/rhoval
c           dzetadnb = -1.0d0/rhoval
c        else
c           dzetadna =  2.0d0*rho(iq,R_B)/(rhoval**2.0d0)
c           dzetadnb = -2.0d0*rho(iq,R_A)/(rhoval**2.0d0)
c        endif
c
c        phi = (1/2)[(1+zeta)**(2/3)+(1-zeta)**(2/3)] and its derivs
c
c        phi = 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0)
c    &               +(1.0d0-zeta)**(2.0d0/3.0d0))
         phi = 0.0d0
         if ((1.0d0+zeta).gt.tol_rho) then
           phi = phi + 0.5d0*(1.0d0+zeta)**(2.0d0/3.0d0)
         endif
         if ((1.0d0-zeta).gt.tol_rho) then
           phi = phi + 0.5d0*(1.0d0-zeta)**(2.0d0/3.0d0)
         endif
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
c
c        t = |Nabla n|/(2*phi*ks*n) and its derivs
c
c        t = dsqgamma/(2.0d0*phi*ks*rhoval)
         t2 = gammaval/((2.0d0*phi*ks*rhoval)**2)
c        dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
c        dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
c
c        { ... } in A (see below) and its derivs
c
         argexp = -epsc/GAMMA/(phi**3)
c        dargexpdna = -depscdna/GAMMA/(phi**3.0d0)
c    &                +3.0d0*epsc/GAMMA/(phi**4.0d0)*dphidna
c        dargexpdnb = -depscdnb/GAMMA/(phi**3.0d0)
c    &                +3.0d0*epsc/GAMMA/(phi**4.0d0)*dphidnb
c
c        A = BETA/GAMMA [exp{-epsilon_c/(GAMMA*phi**3)}-1]**(-1)
c
c        if (abs(argexp).lt.EXPTOL) then
            expinA=exp(argexp)
c        else
c           expinA=0.0d0
c        endif
         A = BETA/GAMMA/(expinA-1.0d0)
c        dAdna = -BETA/GAMMA*dargexpdna*expinA/(expinA-1.0d0)**2.0d0
c        dAdnb = -BETA/GAMMA*dargexpdnb*expinA/(expinA-1.0d0)**2.0d0
c
c        fAt = (1 + A * t**2)/(1 + A * t**2 + A**2 * t**4) and its derivs
c
c        fAtnum = 1.0d0+A*t**2.0d0
c        fAtden = 1.0d0+A*t**2.0d0+A**2.0d0*t**4.0d0
         fAtnum = 1.0d0+A*t2
         fAtden = 1.0d0+A*t2+A**2*t2**2
         fAt = fAtnum/fAtden
c        dfAtnumdt = 2.0d0*A*t
c        dfAtnumdA = t**2.0d0
c        dfAtdendt = 2.0d0*A*t+4.0d0*A**2.0d0*t**3.0d0
c        dfAtdendA = t**2.0d0+2.0d0*A*t**4.0d0
c        dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/(fAtden**2.0d0)
c        dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/(fAtden**2.0d0)
c        dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
c        dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
c
c        arglog = 1 + BETA/GAMMA * t**2 * fAt and its derivs
c
c        arglog = 1.0d0 + BETA/GAMMA*t**2.0d0*fAt
         arglog = 1.0d0 + BETA/GAMMA*t2*fAt
c        darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
c    &                            +t*t*dfAtdna)
c        darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
c    &                            +t*t*dfAtdnb)
c
c        H0 = GAMMA * phi**3 * log{arglog} and its derivs
c
         H0 = GAMMA*(phi**3)*log(arglog)
c        dH0dna = GAMMA*(3.0d0*(phi**2.0d0)*dphidna*log(arglog)
c    &                  +(phi**3.0d0)*darglogdna/arglog)
c        dH0dnb = GAMMA*(3.0d0*(phi**2.0d0)*dphidnb*log(arglog)
c    &                  +(phi**3.0d0)*darglogdnb/arglog)
c
c        Now we update Ec, Amat, and Amat2
c

c          NOTE:  this PBE does the LDA part of Ec in house
            Ec(iq) = Ec(iq)+epsc*wght
            Ec(iq) = Ec(iq)+H0*wght
c           Amat(D1_RA) = Amat(D1_RA) + depscdna
c           if (ipol.eq.2) Amat(D1_RB) = Amat(D1_RB) + depscdnb

c           Amat(D1_RA) = Amat(D1_RA) +  dH0dna
c           if (ipol.eq.2) Amat(D1_RB) = Amat(D1_RB) + dH0dnb
c
c        Now we go into gradient-correction parts
c        Note that the functional depends on |Nabla n| through "t" only
c
c        if (dsqgamma.gt.TOLL)then
c           dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
c           dfAtdg = dfAtdt*dtdg
c           darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
c           dH0dg = GAMMA*(phi**3)*darglogdg/arglog

c           Cmat(D1_GAA) = Cmat(D1_GAA) + dH0dg
c           Cmat(D1_GAB) = Cmat(D1_GAB) + dH0dg*2.0d0
c           if (ipol.eq.2) Cmat(D1_GBB) = Cmat(D1_GBB) + dH0dg
c        endif
   20 continue
      enddo
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_Mpbe.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_Mpbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_Mpbe.F"
#endif
#undef NWAD_PRINT
C>
C> @}
