c Copyright 2018 (C) Orbital-free DFT group at University of Florida
c Licensed under the Educational Community License, Version 2.0 
c (the "License"); you may not use this file except in compliance with 
c the License. You may obtain a copy of the License at
c
c    http://www.osedu.org/licenses/ECL-2.0
c
c Unless required by applicable law or agreed to in writing,
c software distributed under the License is distributed on an "AS IS"
c BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
c or implied. See the License for the specific language governing
c permissions and limitations under the License.
c
#include "dft2drv.fh"
c     Strongly constrained and appropriately normed (SCAN) 
c     functional (Exchange part only)
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c     QTP, Department of Physics, University of Florida
c
c     References:
c     D. Mejia-Rodriguez, S.B. Trickey
c     Phys. Rev. A 96, 052512 (2017)

      Subroutine xc_xscanl(tol_rho, fac,  rho, delrho, laprho,
     &                    Amat, Cmat, Lmat, nq, ipol, Ex, 
     &                    qwght, ldew, func)
      implicit none
c      
      double precision fac, Ex
      integer nq, ipol
      logical ldew
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
c     Charge Density Laplacian
c
      double precision laprho(nq,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*), Lmat(nq,ipol)
c
c     Laplacian of the density
c
      double precision tol_rho
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call xc_xscanl_cs(tol_rho, fac,  rho, delrho, laprho,
     &                    Amat, Cmat, Lmat, nq, Ex, 1d0,
     &                    qwght, ldew, func)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xscanl_cs(tol_rho, fac,  
     R           rho(1,ispin+1), delrho(1,1,ispin), laprho(1,ispin),
     &           Amat(1,ispin), Cmat(1,cmatpos), Lmat(1,ispin),
     &           nq, Ex, 2d0,
     &           qwght, ldew, func)
         enddo

      endif
      return
      end
      Subroutine xc_xscanl_cs(tol_rho, fac,  rho, delrho, laprho,
     &                     Amat, Cmat, Lmat, nq, Ex, facttwo,
     &                     qwght, ldew, func)
      implicit none
c      
      double precision fac, Ex
      integer nq
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(*)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3)
c
c     Charge Density Laplacian
c
      double precision laprho(nq)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq), Cmat(nq), Lmat(nq)
c
      double precision facttwo, afact2 ! 2 for o.s. 1 for c.s.
c
      integer n
      double precision tol_rho, pi
      double precision rhoval, rrho, rho13, rho43, rho53, rho83
      double precision g2, lapval
      double precision p, p14, a, q
c
      double precision thr1, thr2
      double precision F13, F23, F43, F53, F83 
      double precision Ax, Pconst
      double precision rH0, rK1, rA1, rC1, rC2, rD, rMu
      double precision rB1, rB2, rB3, rB4

      double precision oma, oma2
      double precision exp1, exp2, exp3, exp4, exp5
      double precision x1, x2, x
      double precision H, Hden, Hnum
      double precision G
      double precision Fa, Fs, Fx

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivl
      double precision dFada
      double precision dGdp, dHdp, dHdx, dHda
      double precision dxdp, dx1dp, dx2dp, dxda
      double precision dadp, dadq, dpdg, dpdr, dqdr, dqdl
      double precision dfsdp, dfsdq
      double precision dFxda, dFxdp, dFxdr, dFxdg, dFxdl
      
c     functional derivatives above FFFFFFFFFFFF

      parameter (F13=1.0d0/3.0d0)
      parameter (F23=2.0d0*F13)
      parameter (F43=1.0d0+F13)
      parameter (F53=1.0d0+F23)
      parameter (F83=1.0d0+F53)

      parameter (thr1=0.996d0, thr2=1.004d0)

      parameter (rH0=1.174d0)
      parameter (rK1=0.065d0)
      parameter (rA1=4.9479d0)
      parameter (rC1=0.667d0)
      parameter (rC2=0.8d0)
      parameter (rD=1.24d0)
      parameter (rMu=10.0d0/81.0d0)
      parameter (rB3=0.5d0)

      rB2=dsqrt(5913.d0/405000.d0)
      rB1=(511.d0/13500.d0)/(2.d0*rB2)
      rB4=rMu*rMu/rK1-1606.d0/18225.d0-rB1*rB1
c
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then
c
            call ts_pc(tol_rho, rho(n), delrho(n,1:3), laprho(n),
     &                 dfsdp, dfsdq, fs, facttwo)

            rhoval=rho(n)*facttwo
            rho43 = rhoval**F43  
            rrho  = 1d0/rhoval   
            rho13 = rho43*rrho 
            rho83 = rho43*rho43
            rho53 = rho43*rho13
      
            g2 = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)

            g2 = g2 *facttwo*facttwo

            lapval = laprho(n)*facttwo
c     
c     Evaluate the Fx
c     
            p   = g2/(4d0*Pconst*rho83)
            p14 = dsqrt(dsqrt(p))

            q = lapval/(4d0*Pconst*rho53)
c
            a=fs-F53*p
            if(a.lt.0d0)  a=0d0
            oma = 1d0 - a
            oma2 = oma*oma
            
            exp1 = dexp(-rB4/rMu*p)
            exp2 = dexp(-rB3*oma2)
            x1 = rMu*p*(1d0 + rB4/rMu*p*exp1)
            x2 = rB1*p + rB2*oma*exp2

            x = x1 + x2*x2

            Hden = rK1 + x
            Hnum = hden + rK1*x
            H = Hnum/Hden

            if (p14.lt.0.002d0) then
              exp3 = 0d0
            else
              exp3 = dexp(-rA1/p14)
            endif
            G = 1d0 - exp3

            if (a.ge.thr1) then
              exp4 = 0d0
            else
              exp4 = dexp(-rC1*a/oma)
            end if

            if (a.le.thr2) then
              exp5 = 0d0
            else
              exp5 = dexp(rC2/oma)
            end if

            Fa = exp4 - rD*exp5

            Fx = G*(H + Fa*(rH0 - H))

            Ex = Ex + Fx*Ax*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + Fx*Ax*rho43*fac*afact2

c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
            dpdr = -F83*p*rrho
            dpdg = 1d0/(4d0*Pconst*rho83)

            dqdr = -F53*q*rrho
            dqdl = 1.d0/(4d0*Pconst*rho53)

            dadp = dfsdp - F53
            dadq = dfsdq

            if (p14.lt.0.001d0) then
              dGdp = 0d0
            else
              dGdp = -0.25d0*rA1*exp3/(p*p14)
            end if

            dx1dp = rMu + rB4*p*exp1*(2d0 - p*rB4/rMu)
            dx2dp = rB1
            dxdp = dx1dp + 2d0*x2*dx2dp
            dxda = 2d0*rB2*exp2*x2*(2d0*rB3*oma2 - 1d0)

            dHdx = (rK1/Hden)**2
            dHdp = dHdx*dxdp
            dHda = dHdx*dxda

            if ((a.ge.thr1).and.(a.le.thr2)) then
              dFada = 0d0
            else
              dFada = -(rC1*exp4 + rD*exp5*rC2)/oma2
            end if

            dFxdp = dGdp*(H + Fa*(rH0 - H)) + G*dHdp*(1d0 - Fa)
            dFxda = G*(dHda + dFada*(rH0 - H) - Fa*dHda)

            dFxdr = dFxda*(dadp*dpdr + dadq*dqdr) + dFxdp*dpdr
            dFxdg = dFxda*dadp*dpdg + dFxdp*dpdg
            dFxdl = dFxda*dadq*dqdl

            derivr = F43*Ax*rho13*Fx + Ax*rho43*dFxdr
            derivg = Ax*rho43*dFxdg
            derivl = Ax*rho43*dFxdl

            Amat(n) = Amat(n) + derivr*fac
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
            Cmat(n)=  Cmat(n) + 2.0d0*derivg*fac
            Lmat(n)=  Lmat(n) + derivl*fac
         endif
      enddo
      return
      end

      Subroutine xc_xscanl_d2()
      call errquit(' xscanl: d2 not coded ',0,0)
      return
      end
