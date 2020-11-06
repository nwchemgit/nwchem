#include "dft2drv.fh"
c     Deorbitalized version of the regularized and restored 
c     strongly constrained and appropriately normed (r^2SCAN-L) 
c     functional (exchange part only)
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              laprho - Laplacian of density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c
c     References:
c     J.W. Furness, A.D. Kaplan, J. Ning, J.P. Perdew, J. Sun
c     JPCLett 11, 8208-8215 (2020)
c     DOI: 10.1021/acs.jpclett.0c02405
c
c     D. Mejia-Rodriguez, S.B. Trickey
c     PRB 102, 121109(R) (2020)
c     DOI: 10.1103/PhysRevB.102.121109

      Subroutine xc_xr2scanl(tol_rho, fac,  rho, delrho, laprho,
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
      double precision tol_rho
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call xc_xr2scanl_cs(tol_rho, fac,  rho, delrho, laprho,
     &                    Amat, Cmat, Lmat, nq, Ex, 1d0,
     &                    qwght, ldew, func)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xr2scanl_cs(tol_rho, fac,  
     R           rho(1,ispin+1), delrho(1,1,ispin), laprho(1,ispin),
     &           Amat(1,ispin), Cmat(1,cmatpos), Lmat(1,ispin),
     &           nq, Ex, 2d0,
     &           qwght, ldew, func)
         enddo

      endif
      return
      end
      Subroutine xc_xr2scanl_cs(tol_rho, fac,  rho, delrho, laprho,
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
      integer n, ifx
      double precision tol_rho, pi
      double precision rhoval, rrho, rho13, rho43, rho53, rho83
      double precision tauN, tauW, tauU
      double precision p, p14, a, g2

      double precision F13, F23, F43, F53, F83, F18
      double precision Ax, Pconst
      double precision rH0, rK1, rA1, rC1x, rC2x, rD, rMu
      double precision rCeta, rC2, eta, rDp2

      double precision oma, oma2
      double precision exp1, exp2, exp3, exp4, exp5
      double precision rtemp(0:7), rRegu(0:7), rRegu1(7)
      double precision regalpha, dregalpha
      double precision x1, x2, x
      double precision H, Hden, Hnum
      double precision G
      double precision Fx, Fa

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivl
      double precision dFada
      double precision dGdp, dHdp, dHdx, dHda
      double precision dxdp, dx1dp, dx2dp, dxda
      double precision dpdg, dpdr 
      double precision dFxda, dFxdp, dFxdr, dFxdg, dFxdl
      double precision dadp,dadg,dadt,dadr,alpha,dadq
      double precision fs, dfsdp, dfsdq, dqdr, dqdl, q
      double precision lapval
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0)
      parameter (F18=1.d0/8.d0, F53=5.d0/3.d0)

      parameter (rH0=1.174d0)
      parameter (rK1=0.065d0)
      parameter (rA1=4.9479d0)
      parameter (rC1x=0.667d0)
      parameter (rC2x=0.8d0)
      parameter (rD=1.24d0)
      parameter (rMu=10.0d0/81.0d0)
      parameter (eta=1d-3)
      parameter (rCeta=20d0/27d0+eta*f53)
      parameter (rDp2=0.361d0)

      parameter (rRegu = (/ 1d0, -0.667d0, -0.4445555d0, 
     &                    -6.63086601049291d-1, 1.45129704448975d0,
     &                    -8.87998041596655d-1, 2.34528941478571d-1,
     &                    -2.31858433223407d-2/) )
      parameter (rRegu1 = (/ -0.667d0, -2d0*0.4445555d0,
     &             -3d0*6.63086601049291d-1, 4d0*1.45129704448975d0,
     &             -5d0*8.87998041596655d-1, 6d0*2.34528941478571d-1,
     &             -7d0*2.31858433223407d-2/) )
c      parameter (rC2 = -sum(rRegu1)*(1d0-rH0))

c
      pi=acos(-1d0)
      rC2 = -sum(rRegu1)*(1d0-rH0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then

            call ts_pc(tol_rho, rho(n), delrho(n,1:3), laprho(n),
     &                 dfsdp, dfsdq, fs, facttwo)

            rhoval=rho(n)*facttwo
            rho43 = rhoval**F43  ! rho^4/3
            rrho  = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            rho83 = rho43*rho43
            rho53 = rho43*rho13
      
      
            g2 = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)

            g2 = g2 *facttwo*facttwo

            lapval = laprho(n)*facttwo

            tauW = F18*g2*rrho
            tauU = 0.3d0*Pconst*rho53
c     
c     Evaluate the Fx
c     
            p   = g2/(4d0*Pconst*rho83)
            p14 = dsqrt(dsqrt(p))
            q = lapval/(4d0*Pconst*rho53)
c
            alpha = (fs - f53*p)/(1d0 + eta*f53*p)

            oma = 1d0 - alpha
            oma2 = oma*oma
            
            exp1 = dexp(-p**2/rdp2**4)
            x = rCeta*rC2*exp1*p + rMu*p

            Hden = rK1 + x
            Hnum = hden + rK1*x
            H = Hnum/Hden

            if (p14.lt.0.002d0) then
              exp3 = 0d0
            else
              exp3 = dexp(-rA1/p14)
            endif
            G = 1d0 - exp3
c
c Switching function
c
            if (alpha.lt.0.0d0) then
              exp4 = dexp(-rC1x*alpha/oma)
              Fa = exp4
            else if (alpha.lt.2.5d0) then
              rtemp(0) = 1d0
              do ifx=1,7
                rtemp(ifx) = rtemp(ifx-1)*alpha
              enddo
              Fa = dot_product(rRegu,rtemp)
            else
              exp5 = dexp(rC2x/oma)
              Fa = -rD*exp5
            end if

            Fx = G*(H + Fa*(rH0 - H))

            Ex = Ex + Fx*Ax*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + Fx*Ax*rho43*fac*afact2

c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
            dpdr = -F83*p*rrho
            dpdg = 1d0/(4d0*Pconst*rho83)

            dqdr = -F53*q*rrho
            dqdl = 1.d0/(4d0*Pconst*rho53)

            dadp = (dfsdp - f53 - alpha*eta*f53)/(1d0 + eta*f53*p)
            dadq = dfsdq/(1d0 + eta*f53*p)

            if (p14.lt.0.001d0) then
              dGdp = 0d0
            else
              dGdp = -0.25d0*rA1*exp3/(p*p14)
            end if

            dxdp = rMu + rCeta*rC2*exp1*(1d0 - 2d0*p**2/rdp2**4)

            dHdx = (rK1/Hden)**2
            dHdp = dHdx*dxdp

c
c Switching function
c
            if (alpha.lt.0.0d0) then
               dFada = -rC1x*exp4/oma2
            else if (alpha.lt.2.5d0) then
               dFada = dot_product(rRegu1,rtemp(0:6))
            else
               dFada = -rD*exp5*rC2x/oma2
            endif

            dFxdp = dGdp*(H + Fa*(rH0 - H)) + G*dHdp*(1d0 - Fa)
            dFxda = G*(dFada*(rH0 - H))

            dFxdr = dFxda*(dadp*dpdr + dadq*dqdr) + dFxdp*dpdr
            dFxdg = dFxda*(dadp*dpdg) + dFxdp*dpdg
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

      Subroutine xc_xr2scanl_d2()
      call errquit(' xr2scanl: d2 not coded ',0,0)
      return
      end


