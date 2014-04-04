c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_camxlsd(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ex, qwght, ldew, func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c     For locations of 2nd derivatives of functionals in array
#include "dft2drv.fh"
      Subroutine xc_camxlsd_d2(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, nq, ipol, Ex, qwght, ldew, func)
#else
c     For locations of 3rd derivatives of functionals in array
#include "dft3drv.fh"
      Subroutine xc_camxlsd_d3(tol_rho, fac, lfac, nlfac, rho, Amat,
     1                       Amat2, Amat3, nq, ipol, Ex, qwght, ldew, 
     2                       func)
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
      double precision Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
c      double precision Amat2(nq,*)
      double precision Amat2(nq,NCOL_AMAT2)
#endif
#ifdef THIRD_DERIV
c
c     Partial Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,NCOL_AMAT3)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
      double precision rhom23
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
      double precision two_ninth
      Parameter (two_ninth = 2.0d0/9.0d0)
      integer n
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).gt.tol_rho)then
             rho13=rho(n,1)**one_third
             Etmp = rho(n,1)*rho13*P2*fac
             if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
             Atmp = rho13*P1*fac
             Ctmp = 0.0d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
             A2tmp = (rho13/rho(n,1))*2.0d0*one_third*P1*fac
             C2tmp = 0.0d0
             C3tmp = 0.0d0
#endif
#ifdef THIRD_DERIV
             rhom23 = rho13/rho(n,1)
             A3tmp = (rhom23/rho(n,1))*(-4.0d0)*two_ninth*P1*fac
             C4tmp = 0.0d0
             C5tmp = 0.0d0
             C6tmp = 0.0d0
#endif
#ifdef THIRD_DERIV
             call xc_att_xc_d3(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c
             Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
#elif defined(SECOND_DERIV)
             call xc_att_xc_d2(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
c
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#else
             call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp)
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
           if (rho(n,1).gt.tol_rho)then
             rho32=0.0d0
             rho33=0.0d0
             if (rho(n,2).gt.tol_rho) rho32=rho(n,2)**one_third
             if (rho(n,3).gt.tol_rho) rho33=rho(n,3)**one_third
c ---- alpha ----
             Etmp = rho32*rho(n,2)*P4*fac
c
             Atmp = P3*rho32*fac
             Ctmp = 0.0d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
             A2tmp = 0.0d0
             C2tmp = 0.0d0
             C3tmp = 0.0d0
c
             if (rho(n,2).gt.tol_rho) then
               A2tmp = one_third*P3*rho32/rho(n,2)*fac
             endif
#endif
#ifdef THIRD_DERIV
             A3tmp = 0.0d0
             C4tmp = 0.0d0
             C5tmp = 0.0d0
             C6tmp = 0.0d0
c
             if (rho(n,2).gt.tol_rho) then
               A3tmp = -two_ninth*P3*rho32/(rho(n,2)**2)*fac 
             endif
#endif
#ifdef THIRD_DERIV
             if (rho(n,2).gt.tol_rho) then
               call xc_att_xc_d3(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
             endif
c
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c
             Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
#elif defined(SECOND_DERIV)
             if (rho(n,2).gt.tol_rho) then
               call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
             endif
c
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#else
             if (rho(n,2).gt.tol_rho) then
               call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
             endif
#endif
             Amat(n,1) = Amat(n,1) + Atmp
             Ex = Ex + qwght(n)*Etmp
c ---- beta ----
             Etmp = rho33*rho(n,3)*P4*fac
             Atmp = P3*rho33*fac
             Ctmp = 0.0d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
             A2tmp = 0.0d0
             C2tmp = 0.0d0
             C3tmp = 0.0d0
c
             if (rho(n,3).gt.tol_rho) then
               A2tmp = one_third*P3*rho33/rho(n,3)*fac
             end if
#endif
#ifdef THIRD_DERIV
             A3tmp = 0.0d0
             C4tmp = 0.0d0
             C5tmp = 0.0d0
             C6tmp = 0.0d0
c
             if (rho(n,3).gt.tol_rho) then
               A3tmp = -two_ninth*P3*rho33/(rho(n,3)**2)*fac
             endif
#endif
#ifdef THIRD_DERIV
             if (rho(n,3).gt.tol_rho) then
               call xc_att_xc_d3(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
             endif
c
             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c
             Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
#elif defined(SECOND_DERIV)
             if (rho(n,3).gt.tol_rho) then
               call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
             end if
c
             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
#else
             if (rho(n,3).gt.tol_rho) then
               call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
             endif
#endif
             Amat(n,2) = Amat(n,2) + Atmp
             Ex = Ex + qwght(n)*Etmp
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
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c     Compile source again for the 2nd derivative case
#include "xc_camxlsd.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c     Compile source again for the 3rd derivative case
#include "xc_camxlsd.F"
#endif
