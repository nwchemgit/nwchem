#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_ssbD_1.F
C> The SSB-D exchange functional part 1
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the SSB-D exchange functional part 1
C>
C> Evaluates part 1 of the SSB-D exchange functional [1].
C>
C> ### References ###
C>
C> [1] M. Swart, M. Sola, F.M. Bickelhaupt,
C>     "A new all-round density functional based on spin states and
C>     SN2 barriers", J. Chem. Phys. <b>131</b>, 094103 (2009), DOI:
C>     <a href="https://doi.org/10.1063/1.3213193">
C>     10.1063/1.3213193</a>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_ssbD_1(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                         func, Amat, Cmat)
#else
      Subroutine nwpwxc_x_ssbD_1_d2(tol_rho,ipol,nq,wght,rho,rgamma, 
     &                            func, Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id: xc_ssbD_1.F 19696 2010-10-29 16:53:42Z d3y133 $
c
c     SSB-D exchange functional part 1
c         (the one that depends on s)
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] M. Swart, M. Sola, and F.M. Bickelhaupt, JCP 131, 094103 (2009).
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Inputs and Outputs
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision func(nq)   !< [Output] The value of the functional
      double precision amat(nq,*) !< [Output] The derivative wrt rho
      double precision cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
                                   !< and possibly rho
#endif
c
      double precision pi
      double precision rA, rB, rC, rD, rE, rU
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
      parameter (rA=1.079966d0, rB=0.197465d0, rC=0.272729d0)
      parameter (rE=5.873645d0, rU=-0.749940d0)
      parameter (rD=rB*(1.0d0-rU))
c
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
      parameter (F73=7.d0/3.d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2), facSSB
#ifdef SECOND_DERIV
      double precision rhom23, d2s(3), gpp, d2g(3), gssb2
#endif
      double precision gssb0,gssb1
      gssb0(s)= rB*s*s/(1d0+rC*s*s)
     +               - rD*s*s/(1d0+rE*s**4)
      gssb1(s)= 2d0*rB*s/(1d0+rC*s*s)**2 +
     +         (2d0*rD*rE*s**5 - 2d0*rD*s)/(1d0+rE*s**4)**2
#ifdef SECOND_DERIV
      gssb2(s)= 8d0*rB/(1d0+rC*s*s)**3 - 6d0*rB/(1d0+rC*s*s)**2
     +       + 36d0*rD/(1d0+rE*s**4)**2 - 32d0*rD/(1d0+rE*s**4)**3
     -        - 6d0*rD/(1d0+rE*s**4)
#endif
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
            rho43 = C*rho(n,R_T)**F43
            rrho = 1d0/rho(n,R_T)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
c
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = dsqrt(gamma)
            if (.not.(gam12.gt.tol_rho**2)) goto 10
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
            gp=gssb1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + (rho13*g+rho43*d1g(1))*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp=gssb2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*wght*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*wght*4d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*wght*8d0
#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 20 n = 1, nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
c
c     Alpha
c
            if (rho(n,R_A).lt.tol_rho) goto 25
            rho43 = C*(2d0*rho(n,R_A))**F43
            rrho = 0.5d0/rho(n,R_A)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            gamma = rgamma(n,G_AA)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(gam12.gt.tol_rho**2)) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
            gp=gssb1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght*0.5d0
            Amat(n,D1_RA) = Amat(n,D1_RA) + (rho13*g+rho43*d1g(1))*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp=gssb2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*wght*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*wght*0.5d0
#endif
c
c     Beta
c
 25         continue
            if (rho(n,R_B).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,R_B))**F43
            rrho = 0.5d0/rho(n,R_B)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            gamma = rgamma(n,G_BB)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(gam12.gt.tol_rho**2)) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
            gp=gssb1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght*0.5d0
            Amat(n,D1_RB) = Amat(n,D1_RB) + (rho13*g+rho43*d1g(1))*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp=gssb2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*wght*2d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           + rho43*d2g(3)*wght*0.5d0
#endif
c
 20      continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_x_ssbD_1.F"
#endif
C>
C> @}
