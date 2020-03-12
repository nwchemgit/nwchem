#if !defined SECOND_DERIV 
      subroutine xc_lb94_e(tol_rho, xfac, rho, delrho, Amat, nq,
     I     ipol,  Ex, qwght,ldew,func,delrhoq)
#elif defined(SECOND_DERIV)
c     For locations of 2nd derivatives of functionals in array
      subroutine xc_lb94_e_d2(tol_rho, xfac, rho, delrho, Amat, 
     I     Amat2, nq, ipol,  Ex, qwght,ldew,func,delrhoq)
#include "dft2drv.fh"
      implicit none
#endif
c
      double precision tol_rho, xfac
      double precision Ex
      integer nq, ipol
      logical ldew
      double precision func(*),qwght(*)
      double precision rho(nq,*)
      double precision delrho(nq,3,*)
      double precision Amat(nq,*)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*)
#endif
      double precision delrhoq(nq,*) ! [in] 3*rho+ddot(delrho,qxyz)
c
      integer n,ii
      double precision ex_lda   ! ignored
c

      ex_lda=0d0
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)     
      call xc_dirac(tol_rho, 1d0, .true.,.false.,rho, 
     &     amat, nq, ipol,ex_lda,qwght,.false.,func)

      call xc_lb94(tol_rho, 0d0, rho, delrho, 
     &     amat, nq, ipol)
#endif      
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      call xc_dirac_d2(tol_rho, 1d0, .true., .false., rho, Amat,
     &                       Amat2, nq, ipol, ex_lda, qwght, .false., 
     &     func)

      call xc_lb94_d2(tol_rho, 0d0, rho, delrho, 
     &     amat, amat2, nq, ipol)
#endif

c
c     eq. 12 of DOI:10.1080/00268979600100011
c     or
c     eq. 35 of DOI:10.1103/PhysRevA.51.170
c     delrhoq is computed in xc_delrhodotr
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)     
      do ii=1,ipol
         do n=1,nq
            ex=ex+delrhoq(n,ii)*amat(n,ii)*qwght(n)
         enddo
      enddo
#endif      
      return
      end
#if !defined SECOND_DERIV 
      Subroutine xc_lb94(tol_rho, fac, rho, delrho, 
     &                      Amat, nq, ipol)
#elif defined(SECOND_DERIV)
      Subroutine xc_lb94_d2(tol_rho, fac, rho, delrho, 
     &     Amat, Amat2,nq, ipol)
#include "dft2drv.fh"
#endif
c
C$Id$
c
      implicit none
c      
c
      double precision tol_rho, fac
      integer nq, ipol
c
c     Charge Density
c
      double precision rho(nq,*)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,*)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*)
#endif
c
      double precision BETA
      Parameter (BETA = 0.05D0)
c
c References:
c
c    R. van Leeuwen & E. J. Baerends, Phys. Rev. A 49, 2421 (1994).
c
c***************************************************************************
c
      integer n
      double precision arcsinh,darcsinh
      double precision hrho
      double precision rho13, rho43, gamma, x, g, gdenom
c      double precisiobn rhom23
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision dg,dgdenom
#endif
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            hrho  = 0.5d0*rho(n,1)
            rho13 = hrho**(1.d0/3.d0)
            rho43 = rho13*hrho
c            rhom23=rho13/hrho
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,1) = Amat(n,1) + rho13*g*(1d0-fac)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            dgdenom = 3d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) /
     G           (gdenom*gdenom)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     R           BETA/(3d0*rho13*rho13)*(4d0*x*dg - g)
#endif
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,2).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13*rho(n,2)
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,1) = Amat(n,1) + rho13*g*(1d0-fac)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            dgdenom = 3d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) /
     G           (gdenom*gdenom)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     R              BETA/(3d0*rho13*rho13)*(4d0*x*dg - g)
#endif
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
            rho13 = rho(n,3)**(1.d0/3.d0)
            rho43 = rho13*rho(n,3)
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,2) = Amat(n,2) + rho13*g*(1d0-fac)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            dgdenom = 3d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) /
     G           (gdenom*gdenom)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
     R              BETA/(3d0*rho13*rho13)*(4d0*x*dg - g)
#endif
c
 20      continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_lb94.F"
#endif
