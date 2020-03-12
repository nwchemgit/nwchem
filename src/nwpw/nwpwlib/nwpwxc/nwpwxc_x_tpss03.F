C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_tpss03.F
C> The TPSS exchange functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the TPSS exchange functional
C>
C> Evaluate the TPSS meta-GGA exchange functional [1,2]. 
C>
C> ### References ###
C>
C> [1] J. Tao, J.P. Perdew, V.N. Staveroverov, G.E. Scuseria,
C> "Climbing the density functional ladder: Nonemperical meta-
C> generalized gradient approximation designed for molecules
C> and solids",
C> Phys. Rev. Lett. <b>91</b>, 146401-146404 (2003), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.91.146401">
C> 10.1103/PhysRevLett.91.146401</a>.
C>
C> [2] J.P. Perdew, J. Tao, V.N. Staveroverov, G.E. Scuseria,
C> "Meta-generalized gradient approximation: Explanation of a
C> realistic nonempirical density functional",
C> J. Chem. Phys. <b>120</b>, 6898-6911 (2004), DOI:
C> <a href="https://doi.org/10.1063/1.1665298">
C> 10.1103/1.1665298</a>.
C>
c   
c$Id$
c
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
      Subroutine nwpwxc_x_tpss03(tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func, Amat, Cmat, Mmat)
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c      
c     The functional
c
      double precision func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] The derivative wrt tau
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call nwpwxc_x_tpss03_cs(1.0d0, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func, Amat, Cmat, Mmat)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         call nwpwxc_x_tpss03_cs(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_A), rgamma(1,G_AA), tau(1,T_A),
     &                      func, Amat(1,D1_RA), Cmat(1,D1_GAA),
     &                      Mmat(1,D1_TA))
         call nwpwxc_x_tpss03_cs(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_B), rgamma(1,G_BB), tau(1,T_B),
     &                      func, Amat(1,D1_RB), Cmat(1,D1_GBB),
     &                      Mmat(1,D1_TB))
      endif
      return
      end
      Subroutine nwpwxc_x_tpss03_cs(facttwo, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func, Amat, Cmat, Mmat)
      implicit none
c
c     Input and other parameters
c
      double precision facttwo !< [Input] Scale factor
                               !< - 1 for closed shell calculations
                               !< - 2 for open shell calculations
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c      
c     Charge Density
c
      double precision rho(nq) !< [Input] The density
c   
c     Charge Density Gradient
c
      double precision rgamma(nq) !< [Input] The norm of the density gradients
c   
c     Kinetic Energy Density 
c   
      double precision tau(nq) !< [Input] The kinetic energy density
c      
c     The functional
c
      double precision func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq) !< [Output] The derivative wrt rho
      double precision Cmat(nq) !< [Output] The derivative wrt rgamma
      double precision Mmat(nq) !< [Output] The derivative wrt tau
c
      double precision pi
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
         rhoval=rho(n)*facttwo
         if (rhoval.ge.tol_rho) then
         
c     rho43= n*e_x^unif=exchange energy per electron for uniform electron gas
c     = n* Ax*n^(1/3)   or n*C*n^(1/3) 
         
            rho43 = Ax*rhoval**F43 ! Ax*n^4/3
            rrho = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            
C     Below we just sum up the LDA contribution to the functional
            func(n)= func(n) + rho43*wght*afact2
            Amat(n) = Amat(n) + F43*rho13*wght

c     
c           gamma = delrho(n,1)*delrho(n,1) +
c    &           delrho(n,2)*delrho(n,2) +
c    &           delrho(n,3)*delrho(n,3)
            gamma=rgamma(n)
            gamma=gamma*facttwo*facttwo
            tauN = tau(n)*facttwo
            tauW=0.125d0*gamma*rrho
            tauU=0.3d0*P32*rhoval**F53
c     
c     Evaluate the Fx, i.e. mt(x) = Fx - 1 (LDA bit already there)
c     
            p=gamma/(rhoval**F83*P32*4.d0)
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
            
            func(n)= func(n) + mt*rho43*wght*afact2
            
            dmtdn=dxdn/(1.d0+x/kap)**2
            derivn=mt*F43*rho13+rho43*dmtdn
            
            dmtdg=dxdg/(1.d0+x/kap)**2
            derivg = rho43*dmtdg
         
            dmtdt=dxdt/(1.d0+x/kap)**2
            derivt = rho43*dmtdt
            Amat(n) = Amat(n) + derivn*wght
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
            Cmat(n)=  Cmat(n) + 2d0*derivg*wght
            Mmat(n)=  Mmat(n) +0.5d0*derivt*wght
         endif
      enddo
      return
      end

      Subroutine nwpwxc_x_tpss03_d2()
      call errquit(' xtpss03: d2 not coded ',0,0)
      return
      end
C>
C> @}
