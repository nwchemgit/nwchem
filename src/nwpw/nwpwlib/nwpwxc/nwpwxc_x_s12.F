#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_s12.F
C> The S12 exchange functional
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the S12 exchange functional
C>
C> Evaluate the S12 exchange functional [1].
C>
C> ### References ###
C>
C> [1] M. Swart,
C> "A new family of hybrid density functionals",
C> Chem. Phys. Lett. <b>580</b> (2013) 166-171, DOI:
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045">
C> 10.1016/j.cplett.2013.06.045</a>.
C>
#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_s12(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func, Amat, Cmat)
#else
      Subroutine nwpwxc_x_s12_d2(param, tol_rho, ipol, nq, wght, rho,
     &                       rgamma, func, Amat, Amat2, Cmat, Cmat2)
#endif
c
C$Id$
c
      implicit none
c      
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision param(*) !< [Input] Parameters of the functional
      !< (see Table 2 of [1])
      !< - param(1): \f$ A \f$
      !< - param(2): \f$ K \f$ of Eq.(6) in [1] which is identical to
      !<             \f$ \kappa \f$ in the PBE exchange functional
      !< - param(3): \f$ C \f$
      !< - param(4): \f$ D \f$
      !< - param(5): \f$ E \f$
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
c
      double precision rB, rC, rD, rA, rK, rE, rH, rG2, rH2
      double precision ums, vms, dudx, dvdx, d2udx2, d2vdx2
c
c References:
c
c    Swart, Chem. Phys. Lett. (2013), DOI:10.1016/j.cplett.2013.06.045.
c
c***************************************************************************
c
      integer n
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg
      double precision dgdenom, t, x2, x3, x4, g1, g2
      double precision g1h1, g2h1, g1h2, g2h2, g1h3, g3h1
      double precision hdenom, dhdenom, d2hdenom, PI, rM
      double precision gc4, dgc4, d2gc4
      parameter (rM=60.770665d0)
      parameter (PI = 3.1415926535897932385d0)
#ifdef SECOND_DERIV
      double precision rhom23, d2g, d2gdenom, d2g1, d2g2, d2g3
#endif
c
      rA = param(1)
      rK = param(2)
      rC = param(3)
      rD = param(4)
      rE = param(5)

c     if (is12x.eq.1) then
c
cswar1      1.03842032     0.00403198     0.00104596     0.75700000     0.00000000
cswar2      0.00594635     1.17755954     0.84432515     0.00000000     0.00000000
cswar3      1.00000000     0.00000000     0.00000000     0.00000000     1.00000000     0.00000000
c
c       rA = 1.03842032d0
c       rK = 0.757d0
c       rC = 0.00403198d0
c       rD = 0.00104596d0
c       rE = 0.00594635d0
c     elseif (is12x.eq.2) then
c
cswar1      1.02543951     0.00761554     0.00211063     0.75700000     0.25000000
cswar2      0.00604672     1.07735222     0.37705816     0.00000000     0.00000000
cswar3      1.00000000     0.00000000     0.00000000     0.00000000     1.00000000     0.00000000
c
c       rA = 1.02543951d0
c       rK = 0.757d0
c       rC = 0.00761554d0
c       rD = 0.00211063d0
c       rE = 0.00604672d0
c     else
c       stop 'error in xc_s12x.F'
c     endif
      rB = 1d0 + rK - rA
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,R_T))**(1.d0/3.d0)
            rho43 = rho13**4 
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
               x2 = x*x
            else
               x = 0d0
               x2 = 0d0
            endif
c
            gdenom = 1d0 + rC*x2 + rD*x2*x2
            hdenom = 1d0 + rE*x2
            ums = 1d0 - 1d0 / gdenom
            vms = 1d0 - 1d0 / hdenom
            g = C*rB*ums*vms
c
            dudx = (2d0*rC*x + 4d0*rD*x2*x)/(gdenom**2)
            dvdx = 2d0*rE*x/(hdenom**2)
            dg = C*rB*(dudx*vms + ums*dvdx)
c
cmswart            if (lfac) then
cmswart               Ex = Ex + rA*2d0*rho43*C*qwght(n)*fac
cmswart               if(ldew)func(n) = func(n) + rA*2.d0*rho43*C*fac
cmswart               Amat(n,1) = Amat(n,1) + rA*(4d0/3d0)*rho13*C*fac
cmswart            endif
c
c           if (nlfac) then
               func(n) = func(n) + 2.d0*rho43*g*wght
               Amat(n,D1_RA) = Amat(n,D1_RA) 
     +                       + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
            endif
c
#ifdef SECOND_DERIV
cmswart            if(lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,R_T))

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*(rhom23**2/rho(n,R_T))*d2g*wght
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
            endif
#endif
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
            if (rho(n,R_A).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,R_A)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_A)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
               x2 = x*x
            else
               x = 0d0
               x2 = 0d0
            endif
c
            gdenom = 1d0 + rC*x2 + rD*x2*x2
            hdenom = 1d0 + rE*x2
            ums = 1d0 - 1d0 / gdenom
            vms = 1d0 - 1d0 / hdenom
            g = C*rB*ums*vms
c
            dudx = (2d0*rC*x + 4d0*rD*x2*x)/(gdenom**2)
            dvdx = 2d0*rE*x/(hdenom**2)
            dg = C*rB*(dudx*vms + ums*dvdx)
c
cmswart            if (lfac) then
cmswart               Ex = Ex + rA*rho43*C*qwght(n)*fac
cmswart               if (ldew)func(n) = func(n) + rA*rho43*C*fac
cmswart               Amat(n,1) = Amat(n,1) + rA*(4d0/3d0)*rho13*C*fac
cmswart            endif
c
c           if (nlfac) then
               func(n) = func(n) + rho43*g*wght
               Amat(n,D1_RA) = Amat(n,D1_RA)
     +                       + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
            endif
c
#ifdef SECOND_DERIV
cmswart            if (lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,R_A)

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (2d0/3d0)*(rhom23**2/rho(n,R_A))*d2g*wght
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
            endif
#endif
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.tol_rho) goto 20
c
            rho13 = rho(n,R_B)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_B)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
               x2 = x*x
            else
               x = 0d0
               x2 = 0d0
            endif
c
            gdenom = 1d0 + rC*x2 + rD*x2*x2
            hdenom = 1d0 + rE*x2
            ums = 1d0 - 1d0 / gdenom
            vms = 1d0 - 1d0 / hdenom
            g = C*rB*ums*vms
c
            dudx = (2d0*rC*x + 4d0*rD*x2*x)/(gdenom**2)
            dvdx = 2d0*rE*x/(hdenom**2)
            dg = C*rB*(dudx*vms + ums*dvdx)
c
cmswart            if (lfac) then
cmswart               Ex = Ex + rA*rho43*C*qwght(n)*fac
cmswart               if (ldew)func(n) = func(n) + rA*rho43*C*fac
cmswart               Amat(n,2) = Amat(n,2) + rA*(4d0/3d0)*rho13*C*fac
cmswart            endif
c
c           if (nlfac) then
               func(n) = func(n) +rho43*g*wght
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
            endif
c
#ifdef SECOND_DERIV
cmswart            if(lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,R_B)

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)

c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (2d0/3d0)*(rhom23**2/rho(n,R_B))*d2g*wght
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
            endif
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
#include "nwpwxc_x_s12.F"
#endif
C>
C> @}
