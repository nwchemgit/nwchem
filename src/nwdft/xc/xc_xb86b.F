c     Becke86b exchange functional 
c
c     References:
c     [a] A. D. Becke J. Chem. Phys. 85 (1986) 7184
c
#ifndef SECOND_DERIV
      Subroutine xc_xb86b(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_xb86b_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#endif
c
c$Id$
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision tol_rho ! zero-density crit. (input)
      double precision fac ! energy multiplier (input)
      logical lfac, nlfac ! local and non-local labels (input)
      logical ldew        ! compute funct(n) (input)
      double precision Ex ! exchange energy (output)
      integer nq  ! number of grid points (input)
      integer ipol ! spinpolarization flag (input)
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2) ! charge density (in)
      double precision delrho(nq,3,ipol) ! gradient of rho (in)
      double precision qwght(nq) ! quadrature weights (in)
      double precision amat(nq,ipol) ! dex/drho matrix
      double precision cmat(nq,*)  ! 2*dex/(dgrho)
#ifdef SECOND_DERIV
c     second derivatives? dkdc for now.
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2) 
#endif

      double precision pi
      double precision C, Cs
      double precision F13, F43, F45, FM45
      parameter (F13=1d0/3d0,F43=4d0/3d0,F45=4d0/5d0,FM45=-F45)
      integer n
      double precision rrho, rho43, rho73, exnl
      double precision crho43, crho13, crhom23, drho, drho2
      double precision dr2_r43, dr2_r73, dr2_r83, dr4_r5, fac1, fac45
      double precision gparam, bparam, gparamu, bparamu
      parameter (gparam=0.007d0, bparam=0.00375d0)
c     parameter (gparamu=2**(2d0/3d0)*gparam)
c     parameter (bparamu=2**(1d0/3d0)*bparam)
      gparamu=2.0d0**(2d0/3d0)*gparam
      bparamu=2.0d0**(1d0/3d0)*bparam

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
            if (rho(n,1).lt.tol_rho) goto 10

            rho43 = rho(n,1)**F43
            rho73 = rho43 * rho(n,1)
            crho43 = C*rho43
            rrho = 1d0/rho(n,1)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
            if (lfac) then
               Ex = Ex + crho43*qwght(n)*fac

               if(ldew) func(n) = func(n) + crho43*fac
               Amat(n,1) = Amat(n,1) + crho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*crhom23*fac
#endif
            endif

            drho2 = (delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))
            drho = dsqrt(drho2) 
            if (.not.(nlfac.and.drho.gt.tol_rho**2)) goto 10

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparamu * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparamu * dr2_r43 * fac45 * fac
            Ex = Ex + exnl * qwght(n) 
            if (ldew) func(n) = func(n) + exnl

c     GC: dex/drho
            Amat(n,1) = Amat(n,1) + bparamu * fac45 * (F43 * dr2_r73 - 
     &           32d0/15d0 * gparamu * dr4_r5 / fac1) * fac
c     GC: 2*dex/d(grho2)
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2 * exnl / drho2 * 
     &           (1 - F45 * gparamu * dr2_r83 / fac1)

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif
 10      continue
c
      else
c     
c     ======> SPIN-UNRESTRICTED <======
c     
#ifdef IFCV81
C     DEC$ NOSWP
#endif
c
c     alpha
c
         do 100 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 100

            if (rho(n,2).lt.tol_rho) goto 150
            rho43 = (rho(n,2))**F43
            rho73 = rho43 * (rho(n,2))
            crho43 = C*rho43
            rrho = 1.d0/rho(n,2)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
            if (lfac) then
               Ex = Ex + 0.5d0*2**(4d0/3d0)*crho43*qwght(n)*fac
               if(ldew) func(n) = func(n) + 0.5d0*2**(4d0/3d0)*crho43
     &              *fac 
               Amat(n,1) = Amat(n,1) + crho13*2**(1d0/3d0)*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
#endif
            endif

            drho2 = (delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)) 
            drho = dsqrt(drho2) 
            if (.not.(nlfac.and.drho.gt.tol_rho**2)) goto 150

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparam * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparam * dr2_r43 * fac45 * fac
            Ex = Ex + exnl * qwght(n)
            if (ldew) func(n) = func(n) + exnl

c     GC: dex/drho
            Amat(n,1) = Amat(n,1) + bparam * fac45 * 
     &           (F43 * dr2_r73 - 
     &           32d0/15d0 * gparam * dr4_r5 / fac1) * fac
c     GC: 2*dex/d(grho2)
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + exnl / drho2 * 
     &           (1 - F45 * gparam * dr2_r83 / fac1)

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif
c
c     beta
c
 150        continue
            if (rho(n,3).lt.tol_rho) goto 100
            rho43 = rho(n,3)**F43
            rho73 = rho43 * rho(n,3)
            crho43 = C*rho43
            rrho = 1d0/rho(n,3)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
            if (lfac) then
               Ex = Ex + 0.5d0*2**(4d0/3d0)*crho43*qwght(n)*fac
               if(ldew) func(n) = func(n) + 0.5d0*2**(4d0/3d0)*crho43
     &              *fac 
               Amat(n,2) = Amat(n,2) + crho13*2**(1d0/3d0)*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
#endif
            endif

            drho2 = (delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2))
            drho = dsqrt(drho2)  
            if (.not.(nlfac.and.drho.gt.tol_rho**2)) goto 100

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparam * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparam * dr2_r43 * fac45 * fac
            Ex = Ex + exnl * qwght(n)
            if (ldew) func(n) = func(n) + exnl

c     GC: dex/drho
            Amat(n,2) = Amat(n,2) + bparam * fac45 * 
     &           (F43 * dr2_r73 - 
     &           32d0/15d0 * gparam * dr4_r5 / fac1) * fac
c     GC: 2*dex/d(grho2)
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + exnl / drho2 * 
     &           (1 - F45 * gparam * dr2_r83 / fac1)

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif
         
 100     continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_xb86b.F"
#endif
