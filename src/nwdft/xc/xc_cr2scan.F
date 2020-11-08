c     Regularized and restored strongly constrained and 
c     appropriately normed (r^2SCAN) 
c     functional (Correlation part only)
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c
c     References:
c     J.W. Furness, A.D. Kaplan, J. Ning, J.P. Perdew, J. Sun
c     JPCLett 11, 8208-8215 (2020)
c     DOI: 10.1021/acs.jpclett.0c02405

      Subroutine xc_cr2scan(tol_rho, cfac, rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ec, qwght, ldew, func,
     &                    tau, Mmat)     
c
      implicit none
c
#include "errquit.fh"
#include "dft2drv.fh"
c
c     Input and other parameters
c
      integer ipol, nq
      double precision dummy(1)
      double precision cfac
      logical ldew
      double precision func(*)
      double precision fac
      double precision tol_rho
c
c     Threshold parameters
c
      double precision thr1,thr2
      parameter (thr1=0.996d0,thr2=1.004d0)
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol), gammaval, gam12
c
c     Kinetic Energy Density
c
      double precision tau(nq,ipol), tauN
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)
c
c     Intermediate derivatives, results, etc.
c
      integer n, ifc
      double precision ntot,n13,n83,tautot
      double precision dn2,p,rs,rs12,drsdn,dpdg,dpdn
      double precision epsc,depscdrs,depscdzeta
      double precision zeta,omz,opz
      double precision phi,phi2,phi3,dphidzeta
      double precision t2,dt2dp,dt2drs,dt2dzeta
      double precision BETA,dBETAdrs
      double precision A,dAdrs,dAdzeta
      double precision At2,Gaux,Gat2
      double precision w1fac,expw1,w1,dw1drs,dw1dzeta
      double precision arg1,darg1dp,darg1drs,darg1dzeta
      double precision H1,dH1dp,dH1drs,dH1dzeta
      double precision Ec1,dEc1drs,dEc1dp,dEc1dzeta,dEc1dn,dEc1dg
      double precision epsc0,depsc0drs,depsc0dn
      double precision dx,ddxdzeta
      double precision gc,dgcdzeta
      double precision ginf,dginfdp
      double precision w0fac,expw0,w0,dw0drs
      double precision arg0,darg0dp,darg0drs
      double precision H0,dH0dp,dH0drs
      double precision Ec0,dEc0dp,dEc0drs,dEc0dzeta,dEc0dn,dEc0dg
      double precision ds,ddsdzeta
      double precision tueg,tvw
      double precision alpha,dalphadzeta,dalphadn,dalphadg,dalphadt
      double precision oma,oma2
      double precision fca,dfcada,dfcadg,dfcadzeta,dfcadn
      double precision exp5,exp6
      double precision vcpol,vcn
      double precision eclda0, eclda1, declda0drs, d2eclda0drs2
      double precision declda1drs, declda1dz, d2eclda1drs2
      double precision d2eclda1drsdz
      double precision deltaY, deltaYp, ddeltaYpdrs, ddeltaYpdp
      double precision ddeltaYpdz, GAt2y, dGAt2ydrs, dGAt2ydp
      double precision dGAt2ydz,opz23,omz23,exp1

      double precision GAMMA,BETAzero,pi,p14a,p14b,p14f,ckf,ckf2
      parameter (BETAzero = 0.06672455060314922d0)
      parameter (p14a=0.1d0,p14b=0.1778d0)

      double precision b1c,b2c,b3c,c1c,c2c,dc,dxc,xi
      parameter (b1c=0.0285764d0,b2c=0.0889d0,b3c=0.125541d0)
      parameter (c1c=0.64d0,c2c=1.5d0,dc=0.7d0,dxc=2.3631d0)
      parameter (xi=0.12802585262625815d0)

      double precision F4,F13,F23,F43,F53,F83,F14,F8,F5,F16
      parameter (F4=4d0,F13=1d0/3d0,F23=F13+F13,F43=F13+1d0)
      parameter (F53=F23+1d0,F83=F53+1d0,F14=0.25d0)
      parameter (F8=8d0,F5=5d0,F16=1d0/6d0)

      double precision rRegu(0:7), rRegu1(7), rtemp(0:7)
      parameter ( rRegu = (/ 1d0, -0.64d0, -0.4352d0, 
     &            -1.535685604549d0, 3.061560252175d0,
     &            -1.915710236206d0, 5.16884468372d-1,
     &            -5.1848879792d-2 /) )
      parameter ( rRegu1 = (/ -0.64d0, -2d0*0.4352d0,
     &            -3d0*1.535685604549d0, 4d0*3.061560252175d0,
     &            -5d0*1.915710236206d0, 6d0*5.16884468372d-1,
     &            -7d0*5.1848879792d-2 /) )

      double precision eta, rDp2, rfc2
c      parameter (eta=1d-3, rDp2=0.361d0, rfc2 = sum(rRegu1))
      parameter (eta=1d-3, rDp2=0.361d0)
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      Pi = dacos(-1d0)
      rfc2 = sum(rRegu1)
      gamma = (1d0 - dlog(2d0))/pi**2
      p14f = (3d0/(4d0*Pi))**F13
      ckf = (3d0*Pi*Pi)**F13
      ckf2 = ckf*ckf

      do 20 n = 1, nq
c
         ntot=rho(n,1)
         if (ntot.le.tol_rho) goto 20

         n13=ntot**F13
         n83=ntot**F83
         rs=p14f/n13
         rs12=dsqrt(rs)
         drsdn=-F13*rs/ntot
c
c----------------------------------------------------------------------
c functions related to spin-polarization
c----------------------------------------------------------------------
c
         if (ipol.eq.1) then
           zeta = 0d0
           phi = 1d0
           ds = 1d0
           dx = 1d0
           gc = 1d0
         else
           zeta = (rho(n,2) - rho(n,3))/ntot      
           if (zeta.lt.-1d0) zeta=-1d0
           if (zeta.gt. 1d0) zeta= 1d0
           opz = 1d0 + zeta
           omz = 1d0 - zeta
           opz23 = opz**f23
           omz23 = omz**f23
           phi = 0.5d0*(opz23 + omz23)
           ds = 0.5d0*(opz23*opz + omz23*omz)
           dx = 0.5d0*(opz23*opz23 + omz23*omz23)
           gc = (1d0 - dxc*(dx - 1d0))*(1d0 - zeta**12)
         endif
         phi2 = phi*phi
         phi3 = phi2*phi
         tueg = 0.3d0*ckf2*ds*ntot**F53
c         
c----------------------------------------------------------------------
c functions related to delrho  
c----------------------------------------------------------------------
c
         if (ipol.eq.1) then
           dn2=delrho(n,1,1)**2 + delrho(n,2,1)**2 + delrho(n,3,1)**2
         else
           dn2=(delrho(n,1,1)+delrho(n,1,2))**2 +
     &         (delrho(n,2,1)+delrho(n,2,2))**2 +
     &         (delrho(n,3,1)+delrho(n,3,2))**2
         end if
c         
         p=dn2/(F4*ckf2*n83)
         dpdg = 1d0/(F4*ckf2*n83)
         dpdn = -F83*p/ntot
c         
         t2 = ckf2*p/(4d0*phi2*rs*2d0**F23)
         dt2dp = ckf2/(4d0*phi2*rs*2d0**F23)
         dt2drs = -t2/rs
c
         tvw = 0.125d0*dn2/ntot
c
c----------------------------------------------------------------------
c functions related to tau
c----------------------------------------------------------------------
c
         if (ipol.eq.1) then
           tautot=tau(n,1)
         else
           tautot=tau(n,1)+tau(n,2)
         endif

         alpha = (tautot - tvw)/(tueg + eta*tvw)
         oma = 1d0 - alpha
         oma2 = oma*oma

         if (alpha.lt.0d0) then
           exp5 = dexp(-c1c*alpha/oma)
           fca = exp5
           dfcada = -c1c*exp5/oma2
         elseif (alpha.lt.2.5d0) then
           rtemp(0) = 1d0
           do ifc=1,7
             rtemp(ifc) = rtemp(ifc-1)*alpha
           enddo
           fca = dot_product(rRegu,rtemp)
           dfcada = dot_product(rRegu1, rtemp(0:6))
         else
           exp6 = dexp(c2c/oma)
           fca = -dc*exp6
           dfcada = -dc*exp6*c2c/oma2
         endif
c
         dalphadt = 1d0/(tueg + eta*tvw)
         dalphadg = -0.125d0/(ntot*(tueg + eta*tvw))*(1d0 + alpha*eta)
         dalphadn = f53*tueg/(ntot*(tueg+eta*tvw))*
     &                (p/ds*(1d0+eta*alpha) - alpha)
c
c----------------------------------------------------------------------
c EC1
c----------------------------------------------------------------------
c
c eclda1 and derivatives
c
         call lsdac(tol_rho,rs,zeta,eclda1,declda1drs,declda1dz,
     &              d2eclda1drs2,d2eclda1drsdz,dummy)
c
c eclda0 and derivatives
c
         eclda0 = -b1c/(1d0 + b2c*rs12 + b3c*rs)
         declda0drs = -eclda0*(b3c + 0.5d0*b2c/rs12)/
     &                (1d0 + b2c*rs12 + b3c*rs)
         d2eclda0drs2 = eclda0/(1d0 + b2c*rs12 + b3c*rs) *
     &          ( 2d0*(b3c+0.5d0*b2c/rs12)**2/(1d0+b2c*rs12+b3c*rs) +
     &            0.25d0*b2c/rs/rs12 )
c
         BETA = BETAzero*(1d0 + p14a*rs)/(1d0 + p14b*rs)
         dBETAdrs = BETAzero*(p14a - p14b)/(1d0 + p14b*rs)**2
c
         w1fac = eclda1/(GAMMA*phi3)
         expw1 = dexp(-w1fac)
         w1 = expw1 - 1d0
         dw1drs = -expw1*declda1drs/(GAMMA*phi3)

         A = BETA/GAMMA/w1
         dAdrs = dBETAdrs/GAMMA/w1 - A*dw1drs/w1
c         
         At2 = A*t2
         exp1 = dexp(-p**2/rDp2**4)
         deltaY = rfc2/(27d0*gamma*phi3*w1*ds) *
     &              (20d0*rs*(declda0drs-declda1drs) - 
     &               45d0*eta*(eclda0-eclda1))
         deltaYp = deltaY*p*exp1
         ddeltaYpdrs = -deltaYp/w1*dw1drs + 
     &                rfc2/(27d0*gamma*phi3*w1*ds)*
     &                  (20d0*(declda0drs-declda1drs) +
     &                   20d0*rs*(d2eclda0drs2-d2eclda1drs2) - 
     &                   45d0*eta*(declda0drs-declda1drs))*p*exp1
         ddeltaYpdp = deltaY*exp1 - 2d0*p*deltaYp/rDp2**4

         Gaux = 1d0 + 4d0*(At2-deltaYp)
         GAt2y = 1d0/dsqrt(dsqrt(Gaux))
         dGat2ydrs = -Gat2y*(dAdrs*t2 + A*dt2drs - ddeltaYpdrs)/Gaux
         dGat2ydp = -Gat2y*(A*dt2dp - ddeltaYpdp)/Gaux

         arg1 = 1d0 + w1*(1d0 - GAt2y)
         darg1drs = dw1drs*(1d0 - GAt2y) - w1*dGat2ydrs
         darg1dp = -w1*dGAt2ydp

         H1 = GAMMA*phi3*dlog(arg1)
         dH1dp = GAMMA*phi3*darg1dp/arg1
         dH1drs = GAMMA*phi3*darg1drs/arg1

         Ec1 = eclda1 + H1
         dEc1drs = declda1drs + dH1drs
         dEc1dp = dH1dp
c
c--------------------------------------------------------------
c EC0
c--------------------------------------------------------------
c
         w0fac = eclda0/b1c
         expw0 = dexp(-w0fac)
         w0 = expw0 - 1d0
         dw0drs = -declda0drs*expw0/b1c

         ginf = (1d0/(1d0 + 4d0*xi*p))**F14
         dginfdp = -xi*ginf/(1d0 + 4d0*xi*p)

         arg0 = 1d0 + w0*(1d0 - ginf)
         darg0drs = dw0drs*(1d0 - ginf)
         darg0dp = -w0*dginfdp

         H0 = b1c*dlog(arg0)
         dH0dp = b1c*darg0dp/arg0
         dH0drs = b1c*darg0drs/arg0

         Ec0 = (eclda0 + H0)*gc
         dEc0drs = gc*(declda0drs + dH0drs)
         dEc0dp = gc*dH0dp
c
c        --------------------------------------------------------------
c
         Ec = Ec + cfac*ntot*(Ec1 + fca*(Ec0-Ec1))*qwght(n)
         if (ldew) func(n) = func(n) + cfac*ntot*(Ec1 + fca*(Ec0-Ec1))

         dEc1dn = dEc1drs*drsdn + dEc1dp*dpdn
         dEc0dn = dEc0drs*drsdn + dEc0dp*dpdn
         dfcadn = dfcada*dalphadn
         vcn = Ec1 + fca*(Ec0-Ec1) + 
     &         ntot*(dEc1dn + fca*(dEc0dn-dEc1dn) + dfcadn*(Ec0-Ec1))

         Amat(n,1) = Amat(n,1) + cfac*vcn

         dEc1dg = dEc1dp*dpdg
         dEc0dg = dEc0dp*dpdg
         dfcadg = dfcada*dalphadg         

         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + cfac*ntot*(dEc1dg + 
     &   fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1))
         Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + cfac*ntot*(dEc1dg + 
     &   fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1))*2d0

         Mmat(n,1) = Mmat(n,1) + 0.5d0*cfac*
     &                           ntot*dfcada*dalphadt*(Ec0-Ec1)

         if (ipol.eq.2) then
           Amat(n,2) = Amat(n,2) + cfac*vcn

           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + cfac*ntot*(dEc1dg +
     &     fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1))
           
           Mmat(n,2) = Mmat(n,2) + 0.5d0*cfac*
     &                             ntot*dfcada*dalphadt*(Ec0-Ec1)

           if (omz.lt.tol_rho) then
             dphidzeta = 0.5d0*F23*(opz23/opz)
           else if (opz.lt.tol_rho) then
             dphidzeta = -0.5d0*F23*(omz23/omz)
           else
             dphidzeta = 0.5d0*F23*(opz23/opz - omz23/omz)
           end if

           dt2dzeta = -2d0*t2*dphidzeta/phi
           dw1dzeta = (3d0*w1fac*dphidzeta/phi - 
     &                 declda1dz/(GAMMA*phi3))*expw1
           dAdzeta = -A*dw1dzeta/w1

           ddsdzeta = 0.5d0*F53*(opz23 - omz23)
           ddeltaYpdz = -deltaYp/w1*dw1dzeta - deltaYp/ds*ddsdzeta -
     &                  3d0*deltaYp*dphidzeta/phi +
     &                  rfc2/(27d0*gamma*phi3*w1*ds)*(
     &                     -20d0*rs*d2eclda1drsdz +
     &                     45d0*eta*declda1dz)*p*exp1

           dGAt2ydz = -GAt2y*(dAdzeta*t2 + A*dt2dzeta - ddeltaYpdz)/Gaux
           darg1dzeta = dw1dzeta*(1d0-GAt2y) - w1*dGAt2ydz
           dH1dzeta = 3d0*H1*dphidzeta/phi + GAMMA*phi3*darg1dzeta/arg1
           dEc1dzeta = declda1dz + dH1dzeta

           ddxdzeta = 0.5d0*F43*(opz**F13 - omz**F13)
           dgcdzeta = -dxc*ddxdzeta*(1d0 - zeta**12) -
     &                12d0*zeta**11*(1d0 - dxc*(dx - 1d0))
           dEc0dzeta = dgcdzeta*(eclda0 + H0)


           dalphadzeta = -alpha/(tueg+eta*tvw)*tueg*ddsdzeta/ds

           dfcadzeta = dfcada*dalphadzeta

           vcpol = dEc1dzeta + dfcadzeta*(Ec0-Ec1) + 
     &             fca*(dEc0dzeta-dEc1dzeta)

           Amat(n,1) = Amat(n,1) + cfac*omz*vcpol
           Amat(n,2) = Amat(n,2) - cfac*opz*vcpol
         
         end if

20    continue
      end

      Subroutine xc_cr2scan_d2()
      implicit none
      call errquit(' xc_cscan: d2 not coded ',0,0)
      return
      end
