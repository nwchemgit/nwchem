c Copyright 2019 (C) Orbital-free DFT group at University of Florida
c Licensed under the Educational Community License, Version 2.0 
c (the "License"); you may not use this file except in compliance with 
c the License. You may obtain a copy of the License at
c
c    http://www.osedu.org/licenses/ECL-2.0
c
c Unless required by applicable law or agreed to in writing,
c software distributed under the License is distributed on an "AS IS"
c BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
c or implied. See the License for the specific language governing
c permissions and limitations under the License.
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c     Nearly correct asymptotic potential (NCAP) functional
c     (Exchange part only)
c           GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c     QTP, Department of Physics, University of Florida
c
c     References:
c     J. Carmona-Espindola, J.L. Gazquez, A. Vela, S.B. Trickey
c     JCTC 15, 303 (2019).
c     DOI: 10.1021/acs.jctc.8b00998
c
c
#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine xc_xncap(tol_rho, fac, rho, delrho,
     &                    Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
      Subroutine xc_xncap_d2(tol_rho, fac, rho, delrho,
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                       qwght,ldew,func)
#else
      Subroutine xc_xncap_d3(tol_rho, fac, rho, delrho,
     &                       Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &                       nq, ipol, Ex, qwght,ldew,func)
#endif
c
      implicit none
c      
      double precision fac, Ex
      integer nq, ipol
      logical ldew
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
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
c
      double precision tol_rho, pi, mu, alpha, beta, zeta
      double precision ckf, Ax
      double precision F43, F13, F23, F49
c
      parameter(mu=0.2195149727645171d0, beta=0.01808569669d0)
      parameter(zeta=0.30412141859531383d0)

      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0, F23=2.0d0/3.0d0)
      parameter (F49=4d0/9d0)
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, s2, d1s
      double precision arcsinh, darcsinh, d2arcsinh, d3arcsinh
      double precision Fx, dFx, d2Fx, d3Fx
      double precision tansin, dtansin
      double precision factor, dfactor
      double precision asymp, dasymp
      double precision denom, ddenom

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2denom, d2asymp, d2tansin, d2factor, d2s, rhom23
#endif
#if defined(THIRD_DERIV)
      double precision d3denom, d3asymp, d3tansin, d3factor, d3s, rhom53
#endif
      double precision ops, omz, logops

      arcsinh(s)=dlog(s+dsqrt(1d0+s*s))
      darcsinh(s)=1d0/dsqrt(1d0+s*s)
      d2arcsinh(s) = -s/dsqrt(1d0+s*s)**3
      d3arcsinh(s) = (2d0*s*s - 1d0)/dsqrt(1d0+s*s)**5
c
      pi = acos(-1.d0)
      ckf = (3d0*pi*pi)**F13
      Ax = -3d0/(4d0*pi)*ckf
      alpha = 4d0*pi/3d0*beta/mu
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = rho43*rrho
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).le.tol_rho**2) goto 10

            s = dsqrt(gamma)/(2d0*ckf*rho43)
            s2 = s*s

            if (s.lt.1d-8) then
               Fx = 1d0 + mu*s2 + alpha*zeta*mu*s2*s + 
     &              mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s2
               dFx = 2d0*mu*s + 3d0*alpha*zeta*mu*s2 +
     &               4d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s
               d2Fx = 2d0*mu + 6d0*alpha*zeta*mu*s + 
     &                12d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2
               d3Fx = 6d0*alpha*zeta*mu + 
     &                24d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s
            else
              ops = 1d0 + s
              omz = 1d0 - zeta
              logops = dlog(ops)

              asymp = 1d0 + alpha*(omz*s*logops + zeta*s)
              dasymp = alpha*( (s + zeta)/ops + omz*logops )

              if (s.gt.30d0) then
                tansin = arcsinh(s)
                dtansin = darcsinh(s)
              else
                tansin = dtanh(s)*arcsinh(s)
                dtansin = arcsinh(s)/dcosh(s)**2 + dtanh(s)*darcsinh(s)
              endif

              denom = 1d0 + beta*tansin
              ddenom = beta*dtansin

              factor = tansin/denom
              dfactor = (dtansin - factor*ddenom)/denom

              Fx = 1d0 + mu*factor*asymp
              dFx = mu*(dfactor*asymp + factor*dasymp)

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
              if (s.gt.30d0) then
                d2tansin = d2arcsinh(s)
              else
                d2tansin = 2d0*darcsinh(s)/dcosh(s)**2 +
     &                     dtanh(s)*d2arcsinh(s) -
     &                     2d0*tansin/dcosh(s)**2
              endif
              d2asymp = alpha*omz*(2d0 + s)/ops**2
              d2denom = beta*d2tansin
              d2factor = (d2tansin - 2d0*dfactor*ddenom -
     &                    factor*d2denom)/denom
              d2Fx = mu*(d2factor*asymp + 2d0*dfactor*dasymp +
     &                   factor*d2asymp)

#endif
#if defined(THIRD_DERIV)
              if (s.gt.30d0) then
                d3tansin = d3arcsinh(s)
              else
                d3tansin = 3d0*d2arcsinh(s)/dcosh(s)**2 - 
     &                     4d0*darcsinh(s)*dtanh(s)/dcosh(s)**2 + 
     &                     dtanh(s)*d3arcsinh(s) -
     &                     2d0*dtansin/dcosh(s)**2 +
     &                     4d0*tansin*dtanh(s)/dcosh(s)**2
              endif
              d3asymp = -alpha*omz*(3d0 + s)/ops**3
              d3denom = beta*d3tansin
              d3factor = (d3tansin - 3d0*d2factor*ddenom -
     &                    3d0*dfactor*d2denom - factor*d3denom)/denom
              d3Fx = mu*(d3factor*asymp + 3d0*d2factor*dasymp +
     &                   3d0*dfactor*d2asymp + factor*d3asymp)
#endif
            endif

            Ex = Ex + Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + Ax*rho43*Fx*fac

            d1s = 0.5d0*s/gamma
            Amat(n,1) = Amat(n,1) + F43*Ax*rho13*(Fx - s*dFx)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*Ax*rho43*dFx*d1s*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)

            d2s = -0.5d0*d1s/gamma
            rhom23 = rho13*rrho
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*F49*Ax*rhom23*
     &                          (Fx - s*dFx + 4d0*s2*d2Fx)*fac

            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) -
     &                           (F43*Ax*rho13*d2Fx*d1s*s)*4d0*fac

            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + Ax*rho43*
     &                            (d2Fx*d1s**2 + dFx*d2s)*8d0*fac
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho
            d3s = -1.5d0*d2s/gamma

            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx - s*dFx + 18d0*s2*d2Fx +8d0*s2*s*d3Fx)

            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) +
     &      F49*Ax*rhom23*d1s*s*(7d0*d2Fx + 4d0*d3Fx*s)*8d0*fac

            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) -
     &      F43*Ax*rho13*(d2Fx*d1s**2 + d3Fx*d1s**2*s +
     &                     d2Fx*d2s*s)*16d0*fac

            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) +
     &      Ax*rho43*(d3Fx*d1s**3 + 3d0*d2Fx*d1s*d2s +
     &                dFx*d3s)*32d0*fac

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
            if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = (2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = rho43*rrho
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.le.tol_rho**2) goto 25
c
            s = gam12/(2d0*ckf*rho43)
            s2 = s*s

            d1s = 0.5d0*s/gamma

            if (s.lt.1d-8) then
               Fx = 1d0 + mu*s2 + alpha*zeta*mu*s2*s + 
     &              mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s2
               dFx = 2d0*mu*s + 3d0*alpha*zeta*mu*s2 +
     &               4d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s
               d2Fx = 2d0*mu + 6d0*alpha*zeta*mu*s + 
     &                12d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2
               d3Fx = 6d0*alpha*zeta*mu + 
     &                24d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s
            else
              ops = 1d0 + s
              omz = 1d0 - zeta
              logops = dlog(ops)

              asymp = 1d0 + alpha*(omz*s*logops + zeta*s)
              dasymp = alpha*( (s + zeta)/ops + omz*logops )

              if (s.gt.30d0) then
                tansin = arcsinh(s)
                dtansin = darcsinh(s)
              else
                tansin = dtanh(s)*arcsinh(s)
                dtansin = arcsinh(s)/dcosh(s)**2 + dtanh(s)*darcsinh(s)
              endif

              denom = 1d0 + beta*tansin
              ddenom = beta*dtansin

              factor = tansin/denom
              dfactor = (dtansin - factor*ddenom)/denom

              Fx = 1d0 + mu*factor*asymp
              dFx = mu*(dfactor*asymp + factor*dasymp)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
              if (s.gt.30d0) then
                d2tansin = d2arcsinh(s)
              else
                d2tansin = 2d0*darcsinh(s)/dcosh(s)**2 +
     &                     dtanh(s)*d2arcsinh(s) -
     &                     2d0*tansin/dcosh(s)**2
              endif
              d2asymp = alpha*omz*(2d0 + s)/ops**2
              d2denom = beta*d2tansin
              d2factor = (d2tansin - 2d0*dfactor*ddenom -
     &                    factor*d2denom)/denom
              d2Fx = mu*(d2factor*asymp + 2d0*dfactor*dasymp +
     &                   factor*d2asymp)
#endif
#if defined(THIRD_DERIV)
              if (s.gt.30d0) then
                d3tansin = d3arcsinh(s)
              else
                d3tansin = 3d0*d2arcsinh(s)/dcosh(s)**2 - 
     &                     4d0*darcsinh(s)*dtanh(s)/dcosh(s)**2 + 
     &                     dtanh(s)*d3arcsinh(s) -
     &                     2d0*dtansin/dcosh(s)**2 +
     &                     4d0*tansin*dtanh(s)/dcosh(s)**2
              endif
              d3asymp = -alpha*omz*(3d0 + s)/ops**3
              d3denom = beta*d3tansin
              d3factor = (d3tansin - 3d0*d2factor*ddenom -
     &                    3d0*dfactor*d2denom - factor*d3denom)/denom
              d3Fx = mu*(d3factor*asymp + 3d0*d2factor*dasymp +
     &                   3d0*dfactor*d2asymp + factor*d3asymp)
#endif
            endif

            Ex = Ex + 0.5d0*Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + 0.5d0*Ax*rho43*Fx*fac

            Amat(n,1) = Amat(n,1) + F43*Ax*rho13*(Fx - s*dFx)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +
     &                       0.5d0*Ax*rho43*dFx*d1s*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s = -0.5d0*d1s/gamma
            rhom23 = rho13*rrho
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*F49*Ax*rhom23*
     &                          (Fx - s*dFx + 4d0*s2*d2Fx)*fac

            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) -
     &                           (F43*Ax*rho13*d2Fx*d1s*s)*fac

            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + fac*Ax*rho43*
     &                            (d2Fx*d1s**2 + dFx*d2s)*0.5d0
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho
            d3s = -1.5d0*d2s/gamma

            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx - s*dFx + 18d0*s2*d2Fx +8d0*s2*s*d3Fx)

            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) +
     &      F49*Ax*rhom23*d1s*s*(7d0*d2Fx + 4d0*d3Fx*s)*2d0*fac

            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) -
     &      F43*Ax*rho13*(d2Fx*d1s**2 + d3Fx*d1s**2*s +
     &                     d2Fx*d2s*s)*fac

            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) +
     &      Ax*rho43*(d3Fx*d1s**3 + 3d0*d2Fx*d1s*d2s +
     &                dFx*d3s)*0.5d0*fac

#endif
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = (2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = rho43*rrho
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.le.tol_rho**2) goto 20
c
            s = gam12/(2d0*ckf*rho43)
            s2 = s*s

            d1s = 0.5d0*s/gamma

            if (s.lt.1d-8) then
               Fx = 1d0 + mu*s2 + alpha*zeta*mu*s2*s + 
     &              mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s2
               dFx = 2d0*mu*s + 3d0*alpha*zeta*mu*s2 +
     &               4d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2*s
               d2Fx = 2d0*mu + 6d0*alpha*zeta*mu*s + 
     &                12d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s2
               d3Fx = 6d0*alpha*zeta*mu + 
     &                24d0*mu*(alpha - beta - alpha*zeta - 0.5d0)*s
            else
              ops = 1d0 + s
              omz = 1d0 - zeta
              logops = dlog(ops)

              asymp = 1d0 + alpha*(omz*s*logops + zeta*s)
              dasymp = alpha*( (s + zeta)/ops + omz*logops )

              if (s.gt.30d0) then
                tansin = arcsinh(s)
                dtansin = darcsinh(s)
              else
                tansin = dtanh(s)*arcsinh(s)
                dtansin = arcsinh(s)/dcosh(s)**2 + dtanh(s)*darcsinh(s)
              endif

              denom = 1d0 + beta*tansin
              ddenom = beta*dtansin

              factor = tansin/denom
              dfactor = (dtansin - factor*ddenom)/denom

              Fx = 1d0 + mu*factor*asymp
              dFx = mu*(dfactor*asymp + factor*dasymp)

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
              if (s.gt.30d0) then
                d2tansin = d2arcsinh(s)
              else
                d2tansin = 2d0*darcsinh(s)/dcosh(s)**2 +
     &                     dtanh(s)*d2arcsinh(s) -
     &                     2d0*tansin/dcosh(s)**2
              endif
              d2asymp = alpha*omz*(2d0 + s)/ops**2
              d2denom = beta*d2tansin
              d2factor = (d2tansin - 2d0*dfactor*ddenom -
     &                    factor*d2denom)/denom
              d2Fx = mu*(d2factor*asymp + 2d0*dfactor*dasymp +
     &                   factor*d2asymp)
#endif
#if defined(THIRD_DERIV)

              if (s.gt.30d0) then
                d3tansin = d3arcsinh(s)
              else
                d3tansin = 3d0*d2arcsinh(s)/dcosh(s)**2 - 
     &                     4d0*darcsinh(s)*dtanh(s)/dcosh(s)**2 + 
     &                     dtanh(s)*d3arcsinh(s) -
     &                     2d0*dtansin/dcosh(s)**2 +
     &                     4d0*tansin*dtanh(s)/dcosh(s)**2
              endif
              d3asymp = -alpha*omz*(3d0 + s)/ops**3
              d3denom = beta*d3tansin
              d3factor = (d3tansin - 3d0*d2factor*ddenom -
     &                    3d0*dfactor*d2denom - factor*d3denom)/denom
              d3Fx = mu*(d3factor*asymp + 3d0*d2factor*dasymp +
     &                   3d0*dfactor*d2asymp + factor*d3asymp)
#endif
            endif

            Ex = Ex + 0.5d0*Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + 0.5d0*Ax*rho43*Fx*fac

            Amat(n,2) = Amat(n,2) + F43*Ax*rho13*(Fx - s*dFx)*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*Ax*rho43*dFx*d1s*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s = -0.5d0*d1s/gamma
            rhom23 = rho13*rrho
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*F49*Ax*rhom23*
     &                          (Fx - s*dFx + 4d0*s2*d2Fx)*fac

            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) -
     &                           (F43*Ax*rho13*d2Fx*d1s*s)*fac

            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + fac*Ax*rho43*
     &                            (d2Fx*d1s**2 + dFx*d2s)*0.5d0
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho
            d3s = -1.5d0*d2s/gamma

            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx - s*dFx + 18d0*s2*d2Fx +8d0*s2*s*d3Fx)

            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) +
     &      F49*Ax*rhom23*d1s*s*(7d0*d2Fx + 4d0*d3Fx*s)*2d0*fac

            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) -
     &      F43*Ax*rho13*(d2Fx*d1s**2 + d3Fx*d1s**2*s +
     &                     d2Fx*d2s*s)*fac

            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) +
     &      Ax*rho43*(d3Fx*d1s**3 + 3d0*d2Fx*d1s*d2s +
     &                dFx*d3s)*0.5d0*fac
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
#include "xc_xncap.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_xncap.F"
#endif
c $Id$
