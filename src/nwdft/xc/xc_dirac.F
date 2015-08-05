c
c     Modified to handle second and third derivatives while reusing code
c
c     DWS - 6/12
c     BGJ - 8/98
c
#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine xc_dirac(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ex, qwght, ldew, func)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
c     For locations of 2nd derivatives of functionals in array
#include "dft2drv.fh"
      Subroutine xc_dirac_d2(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, nq, ipol, Ex, qwght, ldew, 
     &                       func)
#else
#include "dft2drv.fh"
#include "dft3drv.fh"
      Subroutine xc_dirac_d3(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, Amat3, nq, ipol, Ex, qwght, ldew,
     &                       func)
#endif
c
C$Id$
c
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c
      integer nq, ipol
      double precision fac, Ex
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial First Derivatives of the Exchange Energy Functional
c
      double precision Amat(nq,ipol)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*)
#endif
#ifdef THIRD_DERIV
c
c     Partial Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,*)
#endif
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
c
c Daniel (5-23-12): P1 and P2 are for restricted calculations, P3 
c and P4 are for unrestricted calculations.
c
c     P1 =       -(3/PI)**(1/3)
c     P2 = -(3/4)*(3/PI)**(1/3)
c     P3 =       -(6/PI)**(1/3)
c     P4 = -(3/4)*(6/PI)**(1/3)
c
      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      double precision rho13, rho32, rho33, one_third, term
      double precision two_ninth, rhom23
      Parameter (one_third = 1.d0/3.d0)
      Parameter (two_ninth = 2.d0/9.d0)
      integer n
c
      term = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef SECOND_DERIV
c         write(LuOut,*)'xc_dirac called: restricted'
#endif
         do 10 n = 1, nq
            if (rho(n,1).gt.tol_rho)then
               rho13=rho(n,1)**one_third
               Amat(n,1) = Amat(n,1) + rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     R              1d0/rho(n,1)**(2d0/3d0)
#endif
#ifdef THIRD_DERIV
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + 
     1              1d0/rho(n,1)**(5d0/3d0)
#endif
               term = term + rho(n,1)*rho13*qwght(n)
               if(ldew) func(n) = func(n) + rho(n,1)*rho13*fac*P2
            endif
   10    continue
         call dscal(nq, P1*fac, Amat(1,1), 1)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
         call dscal(nq, 2.0d0*one_third*P1*fac, Amat2(1,D2_RA_RA), 1)
#endif
c
#ifdef THIRD_DERIV
c Daniel (1-11-13): Negative sign comes from repeatedly taking the 
c derivative with respect to the electron density.
         call dscal(nq, -4.0d0*two_ninth*P1*fac, Amat3(1,D3_RA_RA_RA), 
     1     1)
#endif
c 
         Ex = Ex + P2*term*fac
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef SECOND_DERIV
c         write(LuOut,*)'xc_dirac called: unrestricted'
#endif
         do 20 n = 1,nq
            if (rho(n,1).gt.tol_rho)then
               rho32=0.0d0
               rho33=0.0d0
               if (rho(n,2).gt.tol_rho) rho32=rho(n,2)**one_third
               if (rho(n,3).gt.tol_rho) rho33=rho(n,3)**one_third
               Amat(n,1) = Amat(n,1) + P3*rho32*fac
               Amat(n,2) = Amat(n,2) + P3*rho33*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               if (rho(n,2).gt.tol_rho) then
                  Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &                 + one_third*P3*fac*(1d0/rho(n,2)**(2d0/3d0))
               endif
c              Guard against case of no beta electrons, e.g. H atom
               if (rho(n,3).gt.tol_rho) then
                  Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &                 + one_third*P3*fac*(1d0/rho(n,3)**(2d0/3d0))
               endif
#endif
#ifdef THIRD_DERIV
               if (rho(n,2).gt.tol_rho) then
                 Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1                         - two_ninth*P3*rho32/(rho(n,2)**2)*fac
               endif
c              Guard against case of no beta electrons, e.g. H atom
               if (rho(n,3).gt.tol_rho) then
                 Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     1                         - two_ninth*P3*rho33/(rho(n,3)**2)*fac
               endif
#endif
               if (ldew)func(n) = func(n) + ( rho32*rho(n,2) + 
     &                                        rho33*rho(n,3)   )*P4*fac
               term = term + (rho32*rho(n,2)+rho33*rho(n,3))*qwght(n)
c      write(LuOut,*)'from xc_dirac ; rho32, rho33, qwght(n): ',
c     &                           rho32, rho33, qwght(n)
            endif
   20    continue
c 
         Ex = Ex + P4*term*fac
c
      endif
c      write(LuOut,*)'from xc_dirac ; p4, exp, fac, Ex: ', 
c     &                           p4, exp, fac, Ex
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_dirac.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_dirac.F"
#endif
