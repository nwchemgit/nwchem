#ifndef SECOND_DERIV
      Subroutine xc_xpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)     
#else
      Subroutine xc_xpw91_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func)
#endif
c
C$Id$
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
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
      double precision DPOW
      double precision BETA,  tol_rho, CPW91

      Parameter (BETA = 0.0042D0, CPW91=1.6455D0,DPOW=4) ! pw91 paper
C
C

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
c References:
c
c
c***************************************************************************

      integer n
      double precision gamma
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            gamma=0.25d0*gamma
#ifndef SECOND_DERIV
            call xc_xpw91core(DPOW,BETA,n,1,
     &           rho(n,1),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Cmat, nq, ipol, Ex, ldew)     
#else 
               call xc_xpw91core_d2(DPOW,BETA,n,1,
     &           rho(n,1),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &           ldew)
#endif

   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
           if (rho(n,1).lt.tol_rho) goto 20
c
c           Spin alpha:
c
           if (rho(n,2).gt.tol_rho) then
              gamma =    delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,3,1)*delrho(n,3,1)  
#ifndef SECOND_DERIV
               call xc_xpw91core(DPOW,BETA,n,1,
     &           rho(n,2),gamma,qwght(n),func(n),
     &              tol_rho, fac, lfac, nlfac, 
     &              Amat, Cmat, nq, ipol, Ex, ldew)     
#else
               call xc_xpw91core_d2(DPOW,BETA,n,1,
     &           rho(n,2),gamma,qwght(n),func(n),
     &              tol_rho, fac, lfac, nlfac, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         ldew)
#endif
            endif
c
c           Spin beta:
c
            if (rho(n,3).gt.tol_rho) then

            gamma =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2)   
#ifndef SECOND_DERIV
            call xc_xpw91core(DPOW,BETA,n,2,
     &           rho(n,3),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Cmat, nq, ipol, Ex, ldew)     
#else 
               call xc_xpw91core_d2(DPOW,BETA,n,2,
     &           rho(n,3),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &           ldew)
#endif
            endif
c        
   20    continue
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
#include "xc_xpw91.F"
#endif
