c     Perdew-Wang '91 exchange functional 
c
c     References:
c     ****[a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c
#ifndef SECOND_DERIV
      Subroutine xc_xpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_xpw91_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                       qwght,ldew,func)
#endif
c
c$Id: xc_xpw91.F,v 1.5 2000-07-27 16:22:28 bjohnson Exp $
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision amat(nq,ipol), cmat(nq,*)
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
      double precision tol_rho, pi,
     &     P19645, P77956, P2743, P1508, P100, P004
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
      parameter (P19645=0.19645d0, P77956=7.7956d0, P2743=0.2743d0,
     &           P1508=0.1508d0, P100=100d0, P004=0.004d0)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
      parameter (F73=7.d0/3.d0)
#endif
c
      double precision x, arcsinh, darcsinh
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
#ifdef SECOND_DERIV
      double precision d2arcsinh
      d2arcsinh(x)=-x*darcsinh(x)**3
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &     s2, as, d1as, t, d1t, expfac, d1expfac, Fn, d1Fn, Fd, d1Fd,
     &     F, g, gp, d1g(2)
#ifdef SECOND_DERIV
      double precision rhom23, d2s(3), d2as, d2t, d2expfac, d2Fn, d2Fd,
     &     gpp, d2g(3)
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
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 10
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            s2 = s*s
            as = arcsinh(P77956*s)
            d1as = P77956*darcsinh(P77956*s)
#ifdef SECOND_DERIV
            d2as = P77956*P77956*d2arcsinh(P77956*s)
#endif
            t = 1d0 + P19645*s*as
            d1t = P19645*(as+s*d1as)
#ifdef SECOND_DERIV
            d2t = P19645*(2d0*d1as+s*d2as)
#endif
            expfac = -P1508*dexp(-P100*s2)
            d1expfac = expfac*(-2d0*P100*s)
#ifdef SECOND_DERIV
            d2expfac = expfac*(-2d0*P100+(2d0*P100*s)**2)
#endif
            Fn = t + (P2743+expfac)*s2
            d1Fn = d1t + (P2743+expfac)*(2d0*s) + d1expfac*s2
#ifdef SECOND_DERIV
            d2Fn = d2t + (P2743+expfac)*2d0 + d1expfac*(4d0*s)
     &                 + d2expfac*s2
#endif
            Fd = t + P004*s2*s2
            d1Fd = d1t + 4d0*P004*s2*s
#ifdef SECOND_DERIV
            d2Fd = d2t + 12d0*P004*s2
#endif
            F = Fn/Fd
            g = F - 1d0
            gp = (d1Fn-F*d1Fd)/Fd
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp = (d2Fn-2d0*gp*d1Fd-F*d2Fd)/Fd
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac*4d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*fac*8d0
#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = C*(2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            s2 = s*s
            as = arcsinh(P77956*s)
            d1as = P77956*darcsinh(P77956*s)
#ifdef SECOND_DERIV
            d2as = P77956*P77956*d2arcsinh(P77956*s)
#endif
            t = 1d0 + P19645*s*as
            d1t = P19645*(as+s*d1as)
#ifdef SECOND_DERIV
            d2t = P19645*(2d0*d1as+s*d2as)
#endif
            expfac = -P1508*dexp(-P100*s2)
            d1expfac = expfac*(-2d0*P100*s)
#ifdef SECOND_DERIV
            d2expfac = expfac*(-2d0*P100+(2d0*P100*s)**2)
#endif
            Fn = t + (P2743+expfac)*s2
            d1Fn = d1t + (P2743+expfac)*(2d0*s) + d1expfac*s2
#ifdef SECOND_DERIV
            d2Fn = d2t + (P2743+expfac)*2d0 + d1expfac*(4d0*s)
     &                 + d2expfac*s2
#endif
            Fd = t + P004*s2*s2
            d1Fd = d1t + 4d0*P004*s2*s
#ifdef SECOND_DERIV
            d2Fd = d2t + 12d0*P004*s2
#endif
            F = Fn/Fd
            g = F - 1d0
            gp = (d1Fn-F*d1Fd)/Fd
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp = (d2Fn-2d0*gp*d1Fd-F*d2Fd)/Fd
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*fac*0.5d0
#endif
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,2) = Amat(n,2) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*rhom23*fac
#endif
            endif
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            s2 = s*s
            as = arcsinh(P77956*s)
            d1as = P77956*darcsinh(P77956*s)
#ifdef SECOND_DERIV
            d2as = P77956*P77956*d2arcsinh(P77956*s)
#endif
            t = 1d0 + P19645*s*as
            d1t = P19645*(as+s*d1as)
#ifdef SECOND_DERIV
            d2t = P19645*(2d0*d1as+s*d2as)
#endif
            expfac = -P1508*dexp(-P100*s2)
            d1expfac = expfac*(-2d0*P100*s)
#ifdef SECOND_DERIV
            d2expfac = expfac*(-2d0*P100+(2d0*P100*s)**2)
#endif
            Fn = t + (P2743+expfac)*s2
            d1Fn = d1t + (P2743+expfac)*(2d0*s) + d1expfac*s2
#ifdef SECOND_DERIV
            d2Fn = d2t + (P2743+expfac)*2d0 + d1expfac*(4d0*s)
     &                 + d2expfac*s2
#endif
            Fd = t + P004*s2*s2
            d1Fd = d1t + 4d0*P004*s2*s
#ifdef SECOND_DERIV
            d2Fd = d2t + 12d0*P004*s2
#endif
            F = Fn/Fd
            g = F - 1d0
            gp = (d1Fn-F*d1Fd)/Fd
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,2) = Amat(n,2) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            gpp = (d2Fn-2d0*gp*d1Fd-F*d2Fd)/Fd
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           + rho43*d2g(3)*fac*0.5d0
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
#include "xc_xpw91.F"
#endif
