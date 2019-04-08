c   
c$Id$
c
#include "dft2drv.fh"
c    Tao,Perdew, Staroverov, Scuseria exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density
c     References:
c     [a] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         PRL 91, 146401 (2003).
c     [b] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         JCP 120, 6898 (2004).
      Subroutine xc_xtpss03(tol_rho, fac,  rho, delrho, 
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
         call xc_xtpss03_cs(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, Ex, 1d0,
     &                     qwght, ldew, func, tau,Mmat)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xtpss03_cs(tol_rho, fac,  
     R           rho(1,ispin+1), delrho(1,1,ispin), 
     &           Amat(1,ispin), Cmat(1,cmatpos), 
     &           nq, Ex, 2d0,
     &           qwght, ldew, func, 
     T           tau(1,ispin),Mmat(1,ispin))
         enddo

      endif
      return
      end
      Subroutine xc_xtpss03_cs(tol_rho, fac,  rho, delrho, 
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
      double precision tol_rho, pi
      integer n
      double precision rrho, rho43, rho13, gamma
      double precision tauN, tauW, tauU

      double precision  p, qtil, x,  al, mt, z
      double precision   F83, F23, F53, F43, F13
      double precision   G920 
      double precision  b,c,e,es
      double precision    C1, C2, C3
      double precision    kap, mu
      double precision xb,xc,xd
      double precision x1,x2,x3,x4,x5,x6,x7
      double precision   P32, Ax
c     functional derivatives below FFFFFFFFFFFF
      double precision dzdn, dpdn, daldn, dqtdn
      double precision til1, til2
      double precision dtil2dn, dtil1dn
      double precision ax1, bx1, dx1dn
      double precision dx2dn
      double precision dxddn, dxcdn, dx3dn
      double precision dx4dn, dx5dn, dx6dn, dx7dn
      double precision  dxdn
      double precision xmany, dxmanydn
      double precision dmtdn, derivn
      
      double precision dzdg, dpdg, daldg, dqtdg
      double precision dtil2dg
      double precision dx1dg, dx2dg
      double precision dxcdg, dxddg,dx3dg
      double precision dx4dg, dx5dg, dx6dg, dx7dg
      double precision dxmanydg, dxdg, dmtdg, derivg
      
      double precision dzdt, daldt, dqtdt
      double precision dx1dt, dx2dt, dx3dt
      double precision dx5dt
      double precision dxmanydt, dxdt, dmtdt, derivt
      double precision afact2
      double precision rhoval
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter(kap=0.8040d0, mu=0.21951d0)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0, F53=5.d0/3.d0)
      parameter (G920  =9.d0/20.d0 )
      
      parameter(b=0.40d0, c=1.59096d0, e=1.537d0)
      parameter (C1  =  10.d0/81.d0, 
     &     C2  = 146.d0/2025.d0,
     &     C3  = -73.d0/405.d0 )
c     
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      P32 = (3.d0*pi**2)**F23
      es=dsqrt(e)
      afact2=1d0/facttwo
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then
         
c     rho43= n*e_x^unif=exchange energy per electron for uniform electron gas
c     = n* Ax*n^(1/3)   or n*C*n^(1/3) 
         
            rhoval=rho(n)*facttwo
            rho43 = Ax*rhoval**F43 ! Ax*n^4/3
            rrho = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            
C     Below we just sum up the LDA contribution to the functional
            Ex = Ex + rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + rho43*fac*afact2
            Amat(n) = Amat(n) + F43*rho13*fac

c     
            gamma = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)
            gamma=gamma*facttwo*facttwo
            tauN = tau(n,1)*facttwo
            tauW=0.125d0*gamma*rrho
            tauU=0.3d0*P32*rhoval**F53
c     
c     Evaluate the Fx, i.e. mt(x) = Fx - 1 (LDA bit already there)
c     
            p=gamma/(rhoval**F83*P32*4.d0)
            if(abs(p).gt.tol_rho) then
            z=tauW/tauN
            al=(tauN-tauW)/tauU
c     al=dabs(al)
            if(al.lt.0d0)  al=0d0
            
            qtil=(G920*(al-1.d0)/((1.d0+b*al*(al-1.d0))**.5d0)) +
     +           F23*p
            
            xb=(c*z**2)/( (1+z**2)**2 )
            x1=(C1 + xb)*p
            x2=C2*qtil*qtil
            xc=C3*qtil
            xd=(0.5d0*(.6d0*z)**2  + .5d0*p*p)**.5d0
            x3=xc*xd
            x4=C1*C1*p*p/kap
            x5=2.d0*es*C1*(.6d0*z)**2
            x6= e*mu*p*p*p
            x7 = (1.d0+es*p)**(-2.d0)
            
            x=(x1+x2+x3 +x4 +x5+x6)*x7
         
            if (abs(x).lt.tol_rho) write(0,*) ' x for fx ',x
         
c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
C     Derivatives wrt n, density   below
            dzdn=-z*rrho
            dpdn = -p*rrho*F83
            daldn=F53*( -p*dzdn/z**2 +dpdn*(-1.d0+1.d0/z) )
            
            til1=al-1.d0
            til2=(1.d0+b*al*(al-1.d0))**(-0.5d0)
            dtil1dn=daldn
            dtil2dn=b*daldn*(2.d0*al-1d0)*
     &           (-.5d0)*(til2**3)
            dqtdn = G920*(til2*dtil1dn+til1*dtil2dn)+F23*dpdn
            
            ax1=c*p*z*z
            bx1=(1+z*z)**(-2.d0)
            dx1dn=(x1/p)*dpdn + 2d0*c*p*z*dzdn/((1d0+z*z)**3)*(1d0-z*z)
            dx2dn=2.d0*C2*qtil*dqtdn
            
            dxddn=.5d0/xd*( (3d0/5d0)**2*z*dzdn +
     +           p*dpdn)
            dxcdn=C3*dqtdn
            dx3dn=xc*dxddn+xd*dxcdn
            
            dx4dn=(2.d0*x4/p)*dpdn
            dx5dn=(2.d0*x5/z)*dzdn
            dx6dn=(3.d0*x6/p)*dpdn
            dx7dn=-2.d0*es*dpdn/(1.d0+es*p)**3
            
            xmany=x1+x2+x3 +x4 +x5+x6
            dxmanydn= dx1dn+dx2dn+dx3dn+dx4dn+dx5dn+dx6dn
            dxdn=x7*dxmanydn+xmany*dx7dn
C     Derivatives wrt n, density   above
            
C     Derivatives wrt gamma,    below
         
            dpdg=p/gamma
            dzdg=z/gamma
            daldg=(al/p)*dpdg-F53*(p/(z*z))*dzdg
            
            dtil2dg=-0.5d0*daldg*b*(2.d0*al-1d0)*til2**3.d0               
            dqtdg=G920*(til1*dtil2dg+til2*daldg)+F23*dpdg
            dx1dg=(x1/p)*dpdg + 2d0*c*p*z*dzdg/((1d0+z*z)**3)*(1d0-z*z)
         
            dx2dg=C2*2.d0*qtil*dqtdg
         
            dxcdg=C3*dqtdg
            dxddg=.5d0/xd*( (3d0/5d0)**2*z*dzdg +
     +           p*dpdg)
            dx3dg=xc*dxddg+xd*dxcdg
         
            dx4dg=(2.d0*x4/p)*dpdg
            dx5dg=(2.d0*x5/z)*dzdg
            dx6dg=(3.d0*x6/p)*dpdg
         
            dx7dg=-2.d0*es*dpdg*(1.d0+p*es)**(-3.d0)
         
            dxmanydg= dx1dg+dx2dg+dx3dg+dx4dg+dx5dg+dx6dg
            dxdg=x7*dxmanydg+xmany*dx7dg

C     Derivatives wrt tau,    below
c     ttttttttttttttttttttttttttttttttttttttttttttttttt
            dzdt= -z/tauN
            daldt=1.d0/tauU
         
            dqtdt=g920*daldt*til2*(1d0-
     -           0.5d0*b*til1*til2*til2*(2d0*al-1d0))
            
            dx1dt=c*p*dzdt*2d0*z*(1d0-z*z)/((1.d0+z*z)**3)
            dx2dt=2*c2*qtil*dqtdt
            dx3dt=x3*(dqtdt/qtil +
     &           0.5d0*(3d0/5d0)**2*z*dzdt/(xd*xd))
            dx5dt=2d0*(x5/z)*dzdt
            
            dxmanydt= dx1dt+dx2dt+dx3dt+dx5dt
            dxdt=x7*dxmanydt
c     ttttttttttttttttttttttttttttttttttttttttttttttttttt

            mt = kap - kap/(1.d0 + x/kap)
            
            Ex = Ex + mt*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + mt*rho43*fac*afact2
            
            dmtdn=dxdn/(1.d0+x/kap)**2
            derivn=mt*F43*rho13+rho43*dmtdn
            
            dmtdg=dxdg/(1.d0+x/kap)**2
            derivg = rho43*dmtdg
         
            dmtdt=dxdt/(1.d0+x/kap)**2
            derivt = rho43*dmtdt
            Amat(n) = Amat(n) + derivn*fac
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
            Cmat(n)=  Cmat(n) + 2d0*derivg*fac
            Mmat(n)=  Mmat(n) +0.5d0*derivt*fac
         endif
         endif
      enddo
      return
      end

      Subroutine xc_xtpss03_d2()
      call errquit(' xtpss03: d2 not coded ',0,0)
      return
      end


