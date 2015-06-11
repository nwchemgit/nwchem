c   
c$Id: xc_xmvs15.F 21740 2012-01-11 00:25:15Z edo $
c
#include "dft2drv.fh"
c     Made Very Simple (MVS) functional (Exchange part only)
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density
c     References:
c     J. Sun, J.P. Perdew, A. Ruzsinszky
c     PNAS 2015, 112, 685-689
c     DOI: 10.1073/pnas.1423145112 

      Subroutine xc_xmvs15(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau,Mmat)
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
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol), Mmat(nq,*)
      double precision tol_rho
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call xc_xmvs15_cs(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, Ex, 1d0,
     &                     qwght, ldew, func, tau,Mmat)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xmvs15_cs(tol_rho, fac,  
     R           rho(1,ispin+1), delrho(1,1,ispin), 
     &           Amat(1,ispin), Cmat(1,cmatpos), 
     &           nq, Ex, 2d0,
     &           qwght, ldew, func, 
     T           tau(1,ispin),Mmat(1,ispin))
         enddo

      endif
      return
      end
      Subroutine xc_xmvs15_cs(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, Ex, facttwo,
     &                     qwght, ldew, func, tau,Mmat)
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
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq), Cmat(nq)
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,*), Mmat(nq)
c
      double precision facttwo ! 2 for o.s. 1 for c.s.
c
      integer n
      double precision tol_rho, pi
      double precision rrho, rho43, rho13, rho23, rho83
      double precision tauN, tauW, tauU, gamma
      double precision p, a, z, rz, g2, g, uge
      double precision F83, F23, F53, F43, F13, F18
      double precision afact2, Ax, rhoval, Pconst
      double precision rK0, rB, rE1, rC1

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivt
      double precision FaNum, BFaDen, FaDen, Fa
      double precision dFaNumda, dFaDenda, dBFaDenda, dFada
      double precision BFxDen, FxDen, FxNum
      double precision dFxNumda, dFxDendp
      double precision dadp, dadz, dpdg, dpdr, dzdr, dzdg, dzdt
      double precision Fx, dFxdp, dFxdz, dFxdr, dFxdg, dFxdt
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0)
      parameter (F18=1.d0/8.d0, F53=5.d0/3.d0)
      parameter (rK0=0.174d0, rB=0.0233d0,
     &           rE1=-1.66651d0, uge=10.d0/81.d0)
c
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      rC1=(20.d0*rK0/(27.d0*uge))**(4d0)-(1.d0+rE1)*(1.d0+rE1) ! gfortran error if this is defined as a parameter
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then

            rhoval=rho(n)*facttwo
            rho43 = rhoval**F43  ! rho^4/3
            rrho  = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            rho23 = rhoval**F23
            rho83 = rhoval**F83
            
            g2 = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)

            if (dsqrt(g2).gt.tol_rho)then
               g2 = g2 *facttwo*facttwo
               g = dsqrt(g2)
            else
               g  = tol_rho
               g2 = tol_rho*tol_rho
            endif

            tauN = tau(n,1)*facttwo
            tauW = F18*g2*rrho
            tauU = 0.3d0*Pconst*rhoval**F53

c     
c     Evaluate the Fx
c     
            p=g2/(4d0*Pconst*rho83)
            z=TauW/TauN
            rz=TauN/TauW
c
c            a=(tauN-tauW)/tauU
c
            a = F53*p*(rz - 1d0)
            if(a.lt.0d0)  a=0d0
            
            FaNum = (1d0 - a)
            BFaDen = (1d0 + rE1*a*a)**2d0 + rC1*a**4d0
            FaDen = BFaDen**0.25d0
            Fa = Fanum / FaDen

            BFxDen = 1d0 + rB*p*p
            FxDen = BFxDen**F18
            FxNum = 1.d0 + rK0*Fa
            Fx = FxNum / FxDen

            Ex = Ex + Fx*Ax*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + Fx*Ax*rho43*fac*afact2

c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
            dpdr = -F83*p*rrho
            dpdg = p/g2

            dzdr = -z*rrho
            dzdg =  z/g2
            dzdt = -z/TauN

            dadp = a/p
            dadz = F53*p*(-1d0*rz*rz)

            dFaNumda = -1d0
            dBFaDenda = 2d0*(1d0 + rE1*a*a)*2d0*rE1*a + 4d0*rC1*a**3d0
            dFaDenda = 0.25d0*(FaDen/BFaDen)*dBFaDenda
            dFada = (dFaNumda * FaDen - FaNum*dFaDenda)/(FaDen*FaDen)

            dFxNumda = rK0 * dFada
            dFxDendp = F18*(FxDen/BFxDen)*2d0*rB*p

            dFxdp = dFxNumda*dadp/FxDen - FxNum*dFxDendp/(FxDen*FxDen)
            dFxdz = dFxNumda*dadz/FxDen

            dFxdr = dFxdz * dzdr + dFxdp * dpdr
            dFxdg = dFxdz * dzdg + dFxdp * dpdg
            dFxdt = dFxdz * dzdt

            derivr = F43*Ax*rho13*Fx + Ax*rho43*dFxdr
            derivg = Ax*rho43*dFxdg
            derivt = Ax*rho43*dFxdt

            Amat(n) = Amat(n) + derivr*fac
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
            Cmat(n)=  Cmat(n) + 2.0d0*derivg*fac
            Mmat(n)=  Mmat(n) + 0.5d0*derivt*fac
         endif
      enddo
      return
      end

      Subroutine xc_xmvs15_d2()
      call errquit(' xmvs15: d2 not coded ',0,0)
      return
      end


