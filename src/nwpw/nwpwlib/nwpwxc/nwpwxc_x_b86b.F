#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_b86b.F
C> The B86b exchange functional
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the B86b exchange functional
C>
C> Evaluate the non-local part of B86b GGA [1]. This functional 
C> explorers the large gradient limit of the exchange energy.
C>
C> ### References ###
C>
C> [1] A.D. Becke,
C> "On the large gradient behavior of the density functional exchange
C> energy",
C> J. Chem. Phys. <b>85</b>, 7184-7187 (1986), DOI:
C> <a href="https://doi.org/10.1063/1.451353">
C> 10.1063/1.451353</a>.
C>
c
c     Becke86b exchange functional 
c
c     References:
c     [a] A. D. Becke J. Chem. Phys. 85 (1986) 7184
c
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_b86b(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                       func, Amat, Cmat)
#else
      Subroutine nwpwxc_x_b86b_d2(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                          func, Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id: xc_xb86b.F 19696 2010-10-29 16:53:42Z d3y133 $
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
                                    !< and possibly rho
#endif
      double precision pi
      double precision C, Cs,chi
      double precision F13, F43, F45, FM45
      parameter (F13=1d0/3d0,F43=4d0/3d0,F45=4d0/5d0,FM45=-F45)
      integer n
      double precision rrho, rho43, rho73, exnl
      double precision crho43, crho13, crhom23, drho, drho2
      double precision dr2_r43, dr2_r73, dr2_r83, dr4_r5, fac1, fac45
      double precision gparam, bparam, gparamu, bparamu
      parameter (gparam=0.007d0, bparam=0.00375d0)
c     parameter (gparamu=2.0d0**(2d0/3d0)*gparam)
c     parameter (bparamu=2.0d0**(1d0/3d0)*bparam)
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
            if (rho(n,R_T).lt.tol_rho) goto 10

            rho43 = rho(n,R_T)**F43
            rho73 = rho43 * rho(n,R_T)
            crho43 = C*rho43
            rrho = 1d0/rho(n,R_T)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
c           if (lfac) then
c              Ex = Ex + crho43*qwght(n)*fac

c              if(ldew) func(n) = func(n) + crho43*fac
c              Amat(n,1) = Amat(n,1) + crho13*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*crhom23*fac
#endif
c           endif

            drho2 = rgamma(n,G_TT)
c           drho2 = (delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))
            drho = dsqrt(drho2) 
            chi = drho/rho43
            if ((drho.le.tol_rho**2).or.(chi.gt.10.0d0)) goto 10

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparamu * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparamu * dr2_r43 * fac45 
            func(n) = func(n) + exnl * wght

c     GC: dex/drho
            Amat(n,D1_RA) = Amat(n,D1_RA) + bparamu * fac45 * 
     &           (F43 * dr2_r73 - 
     &           32d0/15d0 * gparamu * dr4_r5 / fac1) * wght
c     GC: 2*dex/d(grho2)
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2 * exnl / drho2 * 
     &           (1 - F45 * gparamu * dr2_r83 / fac1) * wght

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
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 100

            if (rho(n,R_A).lt.tol_rho) goto 150
            rho43 = (rho(n,R_A))**F43
            rho73 = rho43 * (rho(n,R_A))
            crho43 = C*rho43
            rrho = 1.d0/rho(n,R_A)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
c           if (lfac) then
c              Ex = Ex + 0.5d0*2**(4d0/3d0)*crho43*qwght(n)*fac
c              if(ldew) func(n) = func(n) + 0.5d0*2**(4d0/3d0)*crho43
c    &              *fac 
c              Amat(n,1) = Amat(n,1) + crho13*2**(1d0/3d0)*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
#endif
c           endif

            drho2 = rgamma(n,G_AA)
c           drho2 = (delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1)) 
            drho = dsqrt(drho2) 
            chi = drho/rho43
            if ((drho.le.tol_rho**2).or.(chi.gt.10.0d0)) goto 150

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparam * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparam * dr2_r43 * fac45 * wght
            func(n) = func(n) + exnl

c     GC: dex/drho
            Amat(n,D1_RA) = Amat(n,D1_RA) + bparam * fac45 * 
     &           (F43 * dr2_r73 - 
     &           32d0/15d0 * gparam * dr4_r5 / fac1) * wght
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
            if (rho(n,R_B).lt.tol_rho) goto 100
            rho43 = rho(n,R_B)**F43
            rho73 = rho43 * rho(n,R_B)
            crho43 = C*rho43
            rrho = 1d0/rho(n,R_B)
            crho13 = F43*crho43*rrho
#ifdef SECOND_DERIV
            crhom23 = F13*crho13*rrho
#endif
c           if (lfac) then
c              Ex = Ex + 0.5d0*2**(4d0/3d0)*crho43*qwght(n)*fac
c              if(ldew) func(n) = func(n) + 0.5d0*2**(4d0/3d0)*crho43
c    &              *fac 
c              Amat(n,2) = Amat(n,2) + crho13*2**(1d0/3d0)*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
#endif
c           endif

            drho2 = rgamma(n,G_BB)
c           drho2 = (delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &           delrho(n,3,2)*delrho(n,3,2))
            drho = dsqrt(drho2)  
            chi = drho/rho43
            if ((drho.le.tol_rho**2).or.(chi.gt.10.0d0)) goto 100

c     Gradient correction
            dr2_r43 = drho2 / rho43
            dr2_r73 = drho2 / rho73
            dr2_r83 = dr2_r43 / rho43
            dr4_r5 = dr2_r73 * dr2_r83 
            fac1 = (1 + gparam * dr2_r83)
            fac45 = fac1**FM45
c     GC: exchange energy and energy density
            exnl = - bparam * dr2_r43 * fac45 * wght
            func(n) = func(n) + exnl

c     GC: dex/drho
            Amat(n,D1_RB) = Amat(n,D1_RB) + bparam * fac45 * 
     &           (F43 * dr2_r73 - 
     &           32d0/15d0 * gparam * dr4_r5 / fac1) * wght
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
#include "nwpwxc_x_b86b.F"
#endif
C>
C> @}
