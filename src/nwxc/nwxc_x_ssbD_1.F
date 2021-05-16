#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_ssbD_1.F
C> The SSB-D exchange functional part 1
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_ssbD_1_p(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                           func)
#else
      Subroutine nwxc_x_ssbD_1(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                         func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_ssbD_1_d2(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                            func)
#else
      Subroutine nwxc_x_ssbD_1_d3(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                            func)
#endif
c
c$Id$
c
c     SSB-D exchange functional part 1
c         (the one that depends on s)
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] M. Swart, M. Sola, and F.M. Bickelhaupt, JCP 131, 094103 (2009).
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
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
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision amat(nq,*) !< [Output] The derivative wrt rho
c     double precision cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c     double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
c                                  !< and possibly rho
#endif
c
      double precision pi
      double precision rA, rB, rC, rD, rE, rU
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
c     double precision F73
#endif
      parameter (rA=1.079966d0, rB=0.197465d0, rC=0.272729d0)
      parameter (rE=5.873645d0, rU=-0.749940d0)
      parameter (rD=rB*(1.0d0-rU))
c
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
c     parameter (F73=7.d0/3.d0)
#endif
c
      integer n
      double precision d1s(2), d, gp, d1g(2), facSSB
      type(nwad_dble)::rrho, rho43, rho13, gamma, gam12, s, g
#ifdef SECOND_DERIV
c     double precision rhom23, d2s(3), gpp, d2g(3), gssb2
#endif
c     type(nwad_dble)::gssb0,gssb1
c     gssb0(s)= rB*s*s/(1d0+rC*s*s)
c    +               - rD*s*s/(1d0+rE*s**4.0d0)
c     gssb1(s)= 2d0*rB*s/(1d0+rC*s*s)**2 +
c    +         (2d0*rD*rE*s**5 - 2d0*rD*s)/(1d0+rE*s**4)**2
#ifdef SECOND_DERIV
c     gssb2(s)= 8d0*rB/(1d0+rC*s*s)**3 - 6d0*rB/(1d0+rC*s*s)**2
c    +       + 36d0*rD/(1d0+rE*s**4)**2 - 32d0*rD/(1d0+rE*s**4)**3
c    -        - 6d0*rD/(1d0+rE*s**4)
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
c           rhom23 = F13*rho13*rrho
#endif
c
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
c           gam12 = sqrt(gamma)
c           if (.not.(gam12.gt.tol_rho**2)) goto 10
c
c           s = Cs*gam12/rho43
            s = Cs*Cs*gamma/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
c           gp=gssb1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) + (rho13*g+rho43*d1g(1))*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c           gpp=gssb2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*d1g(1)
c    &           + rho43*d2g(1))*wght*2d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +(rho13*d1g(2)
c    &           + rho43*d2g(2))*wght*4d0
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           + rho43*d2g(3)*wght*8d0
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
c           rhom23 = F13*rho13*rrho
#endif
            gamma = rgamma(n,G_AA)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
c           gam12 = 2d0*sqrt(gamma)
c           if (.not.(gam12.gt.tol_rho**2)) goto 25
c
c           s = Cs*gam12/rho43
            s = 4.0d0*Cs*Cs*gamma/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
c           gp=gssb1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght*0.5d0
c           Amat(n,D1_RA) = Amat(n,D1_RA) + (rho13*g+rho43*d1g(1))*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c           gpp=gssb2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*d1g(1)
c    &           + rho43*d2g(1))*wght*2d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +(rho13*d1g(2)
c    &           + rho43*d2g(2))*wght
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           + rho43*d2g(3)*wght*0.5d0
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
c           rhom23 = F13*rho13*rrho
#endif
            gamma = rgamma(n,G_BB)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
c           gam12 = 2d0*sqrt(gamma)
c           if (.not.(gam12.gt.tol_rho**2)) goto 20
c
c           s = Cs*gam12/rho43
            s = 4.0d0*Cs*Cs*gamma/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=gssb0(s)
c           gp=gssb1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + rho43*g*wght*0.5d0
c           Amat(n,D1_RB) = Amat(n,D1_RB) + (rho13*g+rho43*d1g(1))*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1g(2)*wght
#ifdef SECOND_DERIV
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c           gpp=gssb2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*d1g(1)
c    &           + rho43*d2g(1))*wght*2d0
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &           +(rho13*d1g(2)
c    &           + rho43*d2g(2))*wght
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &           + rho43*d2g(3)*wght*0.5d0
#endif
c
 20      continue
      endif
c
      return
c
      contains
c
c     The combination of statement functions and derived types with
c     overloaded operators is not properly supported in the Fortran 
c     standard (apparently). Therefore the statement functions from
c     the original subroutine had to be transformed into contained 
c     functions.
c
c     WARNING: This code is EXTREMELY EVIL! Although there appears to be
c     only one function there are actually three with the same name,
c     each one returning results of a different data type. The trick is
c     that depending on the data type the the subroutine that contains
c     these functions changes its name and hence these different
c     functions of the same name do not lead to conflicts. The 
c     alternative would have been to add a forest of conditional 
c     compilation constructs (#ifdef's) to change the function names
c     in the declarations and all places where they are used. That 
c     would have been extremely ugly, so we are between a rock and a
c     hard place on this one.
c
      function gssb0(r) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble), intent(in) :: r
        type(nwad_dble)             :: s
        s = rB*r/(1d0+rC*r)
     +               - rD*r/(1d0+rE*r**2.0d0)
      end function
c
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_ssbD_1.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_ssbD_1.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_ssbD_1.F"
#endif
#undef NWAD_PRINT
C>
C> @}
