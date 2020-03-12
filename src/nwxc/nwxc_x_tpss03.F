#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_tpss03.F
C> The TPSS exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the TPSS exchange functional
C>
C> Evaluate the TPSS meta-GGA exchange functional [1,2]. 
C>
C> Due to the form of the meta-GGAs we need to screen on the kinetic
C> energy density to ensure that LDA will be obtained when the kinetic
C> energy density goes to zero [3].
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
C> [3] J. Gr&auml;fenstein, D. Izotov, D. Cremer,
C>     "Avoiding singularity problems associated with meta-GGA exchange
C>     and correlation functionals containing the kinetic energy
C>     density", J. Chem. Phys. <b>127</b>, 214103 (2007), DOI:
C>     <a href="https://doi.org/10.1063/1.2800011">
C>     10.1063/1.2800011</a>.
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_tpss03_p(tol_rho, ipol, nq, wght, 
     &                           rho, rgamma, tau, func)
#else
      Subroutine nwxc_x_tpss03(tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_tpss03_d2(tol_rho, ipol, nq, wght, 
     &                            rho, rgamma, tau, func)
#else
      Subroutine nwxc_x_tpss03_d3(tol_rho, ipol, nq, wght, 
     &                            rho, rgamma, tau, func)
#endif
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
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
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c      
c     The functional
c
      type(nwad_dble)::func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] The derivative wrt tau
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_x_tpss03_cs_p(1.0d0, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#else
         call nwxc_x_tpss03_cs(1.0d0, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_x_tpss03_cs_d2(1.0d0, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#else
         call nwxc_x_tpss03_cs_d3(1.0d0, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#endif
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_x_tpss03_cs_p(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_A), rgamma(1,G_AA), tau(1,T_A),
     &                      func)
         call nwxc_x_tpss03_cs_p(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_B), rgamma(1,G_BB), tau(1,T_B),
     &                      func)
#else
         call nwxc_x_tpss03_cs(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_A), rgamma(1,G_AA), tau(1,T_A),
     &                      func)
         call nwxc_x_tpss03_cs(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_B), rgamma(1,G_BB), tau(1,T_B),
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_x_tpss03_cs_d2(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_A), rgamma(1,G_AA), tau(1,T_A),
     &                      func)
         call nwxc_x_tpss03_cs_d2(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_B), rgamma(1,G_BB), tau(1,T_B),
     &                      func)
#else
         call nwxc_x_tpss03_cs_d3(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_A), rgamma(1,G_AA), tau(1,T_A),
     &                      func)
         call nwxc_x_tpss03_cs_d3(2.0d0, tol_rho, ipol, nq, wght,
     &                      rho(1,R_B), rgamma(1,G_BB), tau(1,T_B),
     &                      func)
#endif
      endif
      return
      end
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_tpss03_cs_p(facttwo, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#else
      Subroutine nwxc_x_tpss03_cs(facttwo, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_tpss03_cs_d2(facttwo, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#else
      Subroutine nwxc_x_tpss03_cs_d3(facttwo, tol_rho, ipol, nq, wght,
     &                      rho, rgamma, tau, func)
#endif
#include "nwad.fh"
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
      type(nwad_dble)::rho(nq) !< [Input] The density
c   
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq) !< [Input] The norm of the density gradients
c   
c     Kinetic Energy Density 
c   
      type(nwad_dble)::tau(nq) !< [Input] The kinetic energy density
c      
c     The functional
c
      type(nwad_dble)::func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
c     double precision Amat(nq) !< [Output] The derivative wrt rho
c     double precision Cmat(nq) !< [Output] The derivative wrt rgamma
c     double precision Mmat(nq) !< [Output] The derivative wrt tau
c
      double precision pi
      integer n
      type(nwad_dble)::rrho, rho43, gamma
      type(nwad_dble)::tauN, tauW, tauU

      type(nwad_dble):: p, qtil, x,  al, mt, z
      double precision   F83, F23, F53, F43, F13
      double precision   G920 
      double precision  b,c,e,es
      double precision    C1, C2, C3
      double precision    kap, mu
      type(nwad_dble)::xb,xc,xd
      type(nwad_dble)::x1,x2,x3,x4,x5,x6,x7
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
      type(nwad_dble)::rhoval
      
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
c           rho13 = rho43*rrho 
            
C     Below we just sum up the LDA contribution to the functional
            func(n)= func(n) + rho43*wght*afact2
c           Amat(n) = Amat(n) + F43*rho13*wght

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
              if (tauN.ge.tol_rho) then
                z=tauW/tauN
              else
                z=0.0d0
              endif
              al=(tauN-tauW)/tauU
c     al=dabs(al)
              if(al.lt.0d0)  al=0d0
            
              qtil=(G920*(al-1.d0)/((1.d0+b*al*(al-1.d0))**.5d0)) +
     +             F23*p
            
              xb=(c*z**2)/( (1.0d0+z**2)**2)
              x1=(C1 + xb)*p
              x2=C2*qtil*qtil
              xc=C3*qtil
              if (gamma.lt.tol_rho) then
                xd = 0.0d0
              else
                xd=(0.5d0*(.6d0*z)**2 + .5d0*p**2)**.5d0
              endif
              x3=xc*xd
              x4=C1*C1*p**2/kap
              x5=2.d0*es*C1*(.6d0*z)**2
              x6= e*mu*p**3
              x7 = (1.d0+es*p)**(-2)
              
              x=(x1+x2+x3 +x4 +x5+x6)*x7
         
              if (abs(x).lt.tol_rho) write(0,*) ' x for fx ',x
         
c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
C     Derivatives wrt n, density   below
c           dzdn=-z*rrho
c           dpdn = -p*rrho*F83
c           daldn=F53*( -p*dzdn/z**2 +dpdn*(-1.d0+1.d0/z) )
c           
c           til1=al-1.d0
c           til2=(1.d0+b*al*(al-1.d0))**(-0.5d0)
c           dtil1dn=daldn
c           dtil2dn=b*daldn*(2.d0*al-1d0)*
c    &           (-.5d0)*(til2**3)
c           dqtdn = G920*(til2*dtil1dn+til1*dtil2dn)+F23*dpdn
c           
c           ax1=c*p*z*z
c           bx1=(1+z*z)**(-2.d0)
c           dx1dn=(x1/p)*dpdn + 2d0*c*p*z*dzdn/((1d0+z*z)**3)*(1d0-z*z)
c           dx2dn=2.d0*C2*qtil*dqtdn
c           
c           dxddn=.5d0/xd*( (3d0/5d0)**2*z*dzdn +
c    +           p*dpdn)
c           dxcdn=C3*dqtdn
c           dx3dn=xc*dxddn+xd*dxcdn
c           
c           dx4dn=(2.d0*x4/p)*dpdn
c           dx5dn=(2.d0*x5/z)*dzdn
c           dx6dn=(3.d0*x6/p)*dpdn
c           dx7dn=-2.d0*es*dpdn/(1.d0+es*p)**3
c           
c           xmany=x1+x2+x3 +x4 +x5+x6
c           dxmanydn= dx1dn+dx2dn+dx3dn+dx4dn+dx5dn+dx6dn
c           dxdn=x7*dxmanydn+xmany*dx7dn
C     Derivatives wrt n, density   above
            
C     Derivatives wrt gamma,    below
         
c           dpdg=p/gamma
c           dzdg=z/gamma
c           daldg=(al/p)*dpdg-F53*(p/(z*z))*dzdg
c           
c           dtil2dg=-0.5d0*daldg*b*(2.d0*al-1d0)*til2**3.d0
c           dqtdg=G920*(til1*dtil2dg+til2*daldg)+F23*dpdg
c           dx1dg=(x1/p)*dpdg + 2d0*c*p*z*dzdg/((1d0+z*z)**3)*(1d0-z*z)
c        
c           dx2dg=C2*2.d0*qtil*dqtdg
c        
c           dxcdg=C3*dqtdg
c           dxddg=.5d0/xd*( (3d0/5d0)**2*z*dzdg +
c    +           p*dpdg)
c           dx3dg=xc*dxddg+xd*dxcdg
c        
c           dx4dg=(2.d0*x4/p)*dpdg
c           dx5dg=(2.d0*x5/z)*dzdg
c           dx6dg=(3.d0*x6/p)*dpdg
c        
c           dx7dg=-2.d0*es*dpdg*(1.d0+p*es)**(-3.d0)
c        
c           dxmanydg= dx1dg+dx2dg+dx3dg+dx4dg+dx5dg+dx6dg
c           dxdg=x7*dxmanydg+xmany*dx7dg

C     Derivatives wrt tau,    below
c     ttttttttttttttttttttttttttttttttttttttttttttttttt
c           dzdt= -z/tauN
c           daldt=1.d0/tauU
c        
c           dqtdt=g920*daldt*til2*(1d0-
c    -           0.5d0*b*til1*til2*til2*(2d0*al-1d0))
c           
c           dx1dt=c*p*dzdt*2d0*z*(1d0-z*z)/((1.d0+z*z)**3)
c           dx2dt=2*c2*qtil*dqtdt
c           dx3dt=x3*(dqtdt/qtil +
c    &           0.5d0*(3d0/5d0)**2*z*dzdt/(xd*xd))
c           dx5dt=2d0*(x5/z)*dzdt
c           
c           dxmanydt= dx1dt+dx2dt+dx3dt+dx5dt
c           dxdt=x7*dxmanydt
c     ttttttttttttttttttttttttttttttttttttttttttttttttttt

              mt = kap - kap/(1.d0 + x/kap)
            
              func(n)= func(n) + mt*rho43*wght*afact2
            
c           dmtdn=dxdn/(1.d0+x/kap)**2
c           derivn=mt*F43*rho13+rho43*dmtdn
            
c           dmtdg=dxdg/(1.d0+x/kap)**2
c           derivg = rho43*dmtdg
c        
c           dmtdt=dxdt/(1.d0+x/kap)**2
c           derivt = rho43*dmtdt
c           Amat(n) = Amat(n) + derivn*wght
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
c           Cmat(n)=  Cmat(n) + 2d0*derivg*wght
c           Mmat(n)=  Mmat(n) +0.5d0*derivt*wght
         endif
      enddo
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_tpss03.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_tpss03.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_tpss03.F"
#endif
#undef NWAD_PRINT
C>
C> @}
