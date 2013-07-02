#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#ifndef SECOND_DERIV
      Subroutine xc_s12x(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &              Amat, Cmat, nq, ipol, Ex, qwght,ldew,func,is12x)
#else
      Subroutine xc_s12x_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &               Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func,is12x)
#endif
c
C$Id$
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, Ex
      integer nq, ipol, is12x
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
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
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
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

      if (is12x.eq.1) then
c
cswar1      1.03842032     0.00403198     0.00104596     0.75700000     0.00000000
cswar2      0.00594635     1.17755954     0.84432515     0.00000000     0.00000000
cswar3      1.00000000     0.00000000     0.00000000     0.00000000     1.00000000     0.00000000
c
        rA = 1.03842032d0
        rK = 0.757d0
        rC = 0.00403198d0
        rD = 0.00104596d0
        rE = 0.00594635d0
      elseif (is12x.eq.2) then
c
cswar1      1.02543951     0.00761554     0.00211063     0.75700000     0.25000000
cswar2      0.00604672     1.07735222     0.37705816     0.00000000     0.00000000
cswar3      1.00000000     0.00000000     0.00000000     0.00000000     1.00000000     0.00000000
c
        rA = 1.02543951d0
        rK = 0.757d0
        rC = 0.00761554d0
        rD = 0.00211063d0
        rE = 0.00604672d0
      else
        stop 'error in xc_s12x.F'
      endif
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
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
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
            if (nlfac) then
               Ex = Ex + 2d0*rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
            endif
c
#ifdef SECOND_DERIV
cmswart            if(lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
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
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
            endif
c
#ifdef SECOND_DERIV
cmswart            if (lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,2)

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
c
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
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) +rho43*g*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
            endif
c
#ifdef SECOND_DERIV
cmswart            if(lfac) g = g + rA*C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,3)

            d2udx2 = (2d0*rC-6d0*rC*rC*x2+12d0*rD*x2-18d0*rC*rD*x2*x2
     &                -20d0*rD*rD*x2*x2*x2)/(gdenom**3)
            d2vdx2 = (2d0*rE - 6d0*rE*rE*x2)/(hdenom**3)
            d2g = C*rB*(d2udx2*vms + 2d0*dudx*dvdx + ums*d2vdx2)

c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
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
#include "xc_s12x.F"
#endif
