c
c     LDA screened exchange functional
c
c     BGJ - 2/01
c
c     Reference:
c
c     D. M. Bylander and L. Kleinman, Phys. Rev. B 41, 7868 (1990).
c
#ifndef SECOND_DERIV
      Subroutine xc_dirac_screened(tol_rho, fac, lfac, nlfac, rho, 
     &     Amat, nq, ipol, Ex, Ks, qwght, ldew, func)
#else
#include "dft2drv.fh"
      Subroutine xc_dirac_screened_d2(tol_rho, fac, lfac, nlfac, rho,
     &     Amat, Amat2, nq, ipol, Ex, Ks, qwght, ldew, func)
#endif
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision Ks       ! Thomas-Fermi screening constant [input]
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol)
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2)
#endif
c
      double precision ONE3, FOUR3
      Parameter (ONE3 = 1d0/3d0, FOUR3 = 4d0/3d0)
      integer n
      double precision pi, C, TPP, rhoval, rho13, rho43, z, d1z,
     &     f, d1f
      double precision kf, d1kf, f1, f2, d1f1, d1f2
      double precision rz, Tworz, atan2rz, datan2rz, rdatan, lnrdatan,
     &     z243, Eightrz2, Eightrz3, dzdp
#ifdef SECOND_DERIV
      double precision d2z, d2f
c
      call errquit('Must implement 2nd derivs in xc_dirac_screened',0,
     &       CAPMIS_ERR)
#endif
c
      pi = acos(-1d0)
      C = -1.5d0*(0.75d0/pi)**ONE3
      TPP = (3.d0*pi*pi)**ONE3
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rhoval = 0.5d0*rho(n,1)
c
c           Spin alpha:
c
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c
            kf = TPP*rho13
            z = Ks/kf
            rz = 1d0/z
            Tworz = 2d0*rz
            Eightrz2 = 4d0*Tworz*rz
            Eightrz3 = Eightrz2*rz
            rdatan = 1d0 + Tworz*Tworz
            datan2rz = 1d0/rdatan
            lnrdatan = log(rdatan)
            atan2rz = atan(Tworz)
            z243 = z*z/4d0 + 3d0
            f1 = -z*z/6d0
            f2 = 1d0 - z243*lnrdatan + 8d0*rz*atan2rz
            d1f1 = -z/3d0
            d1f2 = -z*0.5d0*lnrdatan - Eightrz2*atan2rz
     &           - Eightrz3*datan2rz*(2d0 - z243)
            f = 1d0 + f1*f2
            dzdp = -ONE3*z/rhoval
            d1f = (d1f1*f2 + f1*d1f2)*dzdp
c
            if (lfac) then
               Ex = Ex + 2d0*rho43*C*f*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*C*f*fac
               Amat(n,1) = Amat(n,1) + (FOUR3*rho13*f+rho43*d1f)*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(FOUR3*rho13/rhoval)*C*fac
#endif
            endif
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
            rhoval = rho(n,2)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c
            kf = TPP*rho13
            z = Ks/kf
            rz = 1d0/z
            Tworz = 2d0*rz
            Eightrz2 = 4d0*Tworz*rz
            Eightrz3 = Eightrz2*rz
            rdatan = 1d0 + Tworz*Tworz
            datan2rz = 1d0/rdatan
            lnrdatan = log(rdatan)
            atan2rz = atan(Tworz)
            z243 = z*z/4d0 + 3d0
            f1 = -z*z/6d0
            f2 = 1d0 - z243*lnrdatan + 8d0*rz*atan2rz
            d1f1 = -z/3d0
            d1f2 = -z*0.5d0*lnrdatan - Eightrz2*atan2rz
     &           - Eightrz3*datan2rz*(2d0 - z243)
            f = 1d0 + f1*f2
            dzdp = -ONE3*z/rhoval
            d1f = (d1f1*f2 + f1*d1f2)*dzdp
c
            if (lfac) then
               Ex = Ex + rho43*C*f*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*C*f*fac
               Amat(n,1) = Amat(n,1) + (FOUR3*rho13*f+rho43*d1f)*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(FOUR3*rho13/rhoval)*C*fac
#endif
            endif
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
            rhoval = rho(n,3)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c
            kf = TPP*rho13
            z = Ks/kf
            rz = 1d0/z
            Tworz = 2d0*rz
            Eightrz2 = 4d0*Tworz*rz
            Eightrz3 = Eightrz2*rz
            rdatan = 1d0 + Tworz*Tworz
            datan2rz = 1d0/rdatan
            lnrdatan = log(rdatan)
            atan2rz = atan(Tworz)
            z243 = z*z/4d0 + 3d0
            f1 = -z*z/6d0
            f2 = 1d0 - z243*lnrdatan + 8d0*rz*atan2rz
c            f2 = 1d0 - z243*log(1d0+rz) + 8d0*rz*atan2rz
            d1f1 = -z/3d0
            d1f2 = -z*0.5d0*lnrdatan - Eightrz2*atan2rz
     &           - Eightrz3*datan2rz*(2d0 - z243)
            f = 1d0 + f1*f2
            dzdp = -ONE3*z/rhoval
            d1f = (d1f1*f2 + f1*d1f2)*dzdp
c
            if (lfac) then
               Ex = Ex + rho43*C*f*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*C*f*fac
               Amat(n,2) = Amat(n,2) + (FOUR3*rho13*f+rho43*d1f)*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + ONE3*(FOUR3*rho13/rhoval)*C*fac
#endif
            endif
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
#include "xc_dirac_screened.F"
#endif
