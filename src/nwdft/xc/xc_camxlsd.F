c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#ifndef SECOND_DERIV
      Subroutine xc_camxlsd(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ex, qwght, ldew, func)
#else
c     For locations of 2nd derivatives of functionals in array
#include "dft2drv.fh"
      Subroutine xc_camxlsd_d2(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, nq, ipol, Ex, qwght, ldew, func)
#endif
c
C$Id: xc_camxlsd.F,v 1.2 2009-02-18 06:12:48 niri Exp $
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
      double precision Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp
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
      double precision rho13, rho32, rho33, one_third
      Parameter (one_third = 1.d0/3.d0)
      integer n
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (dabs(rho(n,1)).gt.tol_rho)then
             rho13=dabs(rho(n,1))**one_third*sign(1d0,rho(n,1))
             Etmp = rho(n,1)*rho13*P2*fac
             if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
             Atmp = rho13*P1*fac
             Ctmp = 0.d0
             call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp)
#ifdef SECOND_DERIV
             A2tmp = (rho13/rho(n,1))*2.0d0*one_third*P1*fac
             C2tmp = 0.d0
             C3tmp = 0.d0
             call xc_att_xc_d2(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#endif
             Ex = Ex + qwght(n)*Etmp
             Amat(n,1) = Amat(n,1) + Atmp
            endif
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (dabs(rho(n,1)).gt.tol_rho)then
             rho32=dabs(rho(n,2))**one_third*sign(1d0,rho(n,2))
             rho33=dabs(rho(n,3))**one_third*sign(1d0,rho(n,3))
c
             Etmp = rho32*rho(n,2)*P4*fac
             Atmp = P3*rho32*fac
             Ctmp = 0.d0
             call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
             Amat(n,1) = Amat(n,1) + Atmp
             Ex = Ex + qwght(n)*Etmp
#ifdef SECOND_DERIV
             A2tmp = 0.d0
             C2tmp = 0.d0
             C3tmp = 0.d0
             if (dabs(rho(n,2)).gt.tol_rho) then
               A2tmp = one_third*P3*rho32/rho(n,2)*fac
               call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
             end if
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#endif
c
             Etmp = rho33*rho(n,3)*P4*fac
             Atmp = P3*rho33*fac
             Ctmp = 0.d0
             call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
             Amat(n,2) = Amat(n,2) + Atmp
             Ex = Ex + qwght(n)*Etmp
#ifdef SECOND_DERIV
             A2tmp = 0.d0
             C2tmp = 0.d0
             C3tmp = 0.d0
             if (dabs(rho(n,3)).gt.tol_rho) then
               A2tmp = one_third*P3*rho33/rho(n,3)*fac
               call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
             end if
             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
#endif
c
             if (ldew)func(n) = func(n) + ( rho32*rho(n,2) + 
     &                                      rho33*rho(n,3)   )*P4*fac
            endif
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
#include "xc_camxlsd.F"
#endif
