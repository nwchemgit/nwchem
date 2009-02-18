c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#ifndef SECOND_DERIV
      Subroutine xc_dirac(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ex, qwght, ldew, func)
#else
c     For locations of 2nd derivatives of functionals in array
#include "dft2drv.fh"
      Subroutine xc_dirac_d2(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, nq, ipol, Ex, qwght, ldew, func)
#endif
c
C$Id: xc_dirac.F,v 1.20 2009-02-18 02:34:15 niri Exp $
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
#ifdef SECOND_DERIV
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*)
#endif
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
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
      Parameter (one_third = 1.d0/3.d0)
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
            if (abs(rho(n,1)).gt.tol_rho)then
               rho13=abs(rho(n,1))**one_third*sign(1d0,rho(n,1))
               Amat(n,1) = Amat(n,1) + rho13
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + rho13/rho(n,1)
#endif
               term = term + rho(n,1)*rho13*qwght(n)
               if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
            endif
   10    continue
         call dscal(nq, P1*fac, Amat(1,1), 1)
#ifdef SECOND_DERIV
         call dscal(nq, 2.0d0*one_third*P1*fac, Amat2(1,D2_RA_RA), 1)
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
            if (abs(rho(n,1)).gt.tol_rho)then
! we can have beta negative density eg he2+
!               if (rho(n,2).lt.0.d0.or.rho(n,3).lt.0.d0)then
!                  write(LuOut,*)
!     &                      ' n, nq, rho(n,1), rho(n,2), rho(n,3): ',
!     &                        n, nq, rho(n,1), rho(n,2), rho(n,3)
!                  call errquit(' negative density ',911)
!               endif
               rho32=abs(rho(n,2))**one_third*sign(1d0,rho(n,2))
               rho33=abs(rho(n,3))**one_third*sign(1d0,rho(n,3))
               Amat(n,1) = Amat(n,1) + P3*rho32*fac
               Amat(n,2) = Amat(n,2) + P3*rho33*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &                           + one_third*P3*rho32/rho(n,2)*fac
c              Guard against case of no beta electrons, e.g. H atom
               if (rho(n,3).gt.tol_rho) then
                  Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &                              + one_third*P3*rho33/rho(n,3)*fac
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
