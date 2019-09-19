#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_pw6.F
C> The PW6 exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the PW6 exchange functional
C>
C> The PW6 exchange functional [1] is a 6 parameter functional based
C> on the Perdew-Wang-91 exchange functional and the Becke-95 
C> correlation functional. This subroutine implements just the exchange
C> part but the parameters were co-optimized with the correlation
C> functional. Also this subroutine just implements the gradient 
C> correction term (and not the LDA term).
C>
C> ### References ###
C>
C> [1] Y. Zhao, D. G. Truhlar,
C>     "Design of density functionals that are broadly accurate for
C>     thermochemistry, thermochemical kinetics, and nonbonded
C>     interactions", J. Phys. Chem. A <b>109</b> 5656-5667 (2005),
C>     DOI: 
C>     <a href="https://doi.org/10.1021/jp050536c">
C>     10.1021/jp050536c</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_pw6_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, func)
#else
      Subroutine nwxc_x_pw6(param, tol_rho, ipol, nq, wght, rho, rgamma,
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_pw6_d2(param, tol_rho, ipol, nq, wght,
     &                         rho, rgamma, func)
#else
      Subroutine nwxc_x_pw6_d3(param, tol_rho, ipol, nq, wght,
     &                         rho, rgamma, func)
#endif
c
C$Id$
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)!< [Input] Parameters of functional
      type(nwad_dble)::DPOW, BETA, CPW91
#else
      double precision param(*)!< [Input] Parameters of functional
      double precision DPOW, BETA, CPW91
#endif
#else
      double precision param(*)!< [Input] Parameters of functional
                               !< (see [1] Eq.(5) and Table 2):
                               !< - param(1): \f$ b \f$
                               !< - param(2): \f$ c \f$
                               !< - param(3): \f$ d \f$
      double precision DPOW, BETA, CPW91
#endif
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
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision amat(nq,*) !< [Output] The derivative wrt rho
c     double precision cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c     double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
c                                  !< and possibly rho
#endif
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c

c     Data BETA/0.00538d0/, CPW91/1.7382d0/,DPOW/3.8901d0/! PW6B95

c
c***************************************************************************

      integer n
      type(nwad_dble)::gamma
c
c parameters for PWB6K
c
c     if (ijzy.eq.2) then
c       BETA = 0.00539d0
c       CPW91= 1.7077d0
c       DPOW= 4.0876d0
c     endif
c
      BETA  = param(1)
      CPW91 = param(2)
      DPOW  = param(3)
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1)
            gamma=0.25d0*gamma
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
            call nwxc_x_pw6core_p(DPOW,BETA,CPW91,n,1,
     &           rho(n,R_T),gamma,func(n),
     &           tol_rho, wght, 
     &           nq, ipol)     
#else
            call nwxc_x_pw6core(DPOW,BETA,CPW91,n,1,
     &           rho(n,R_T),gamma,func(n),
     &           tol_rho, wght, 
     &           nq, ipol)     
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxc_x_pw6core_d2(DPOW,BETA,CPW91,n,1,
     &           rho(n,R_T),gamma,func(n),
     &           tol_rho, wght, 
     &           nq, ipol)
#else 
            call nwxc_x_pw6core_d3(DPOW,BETA,CPW91,n,1,
     &           rho(n,R_T),gamma,func(n),
     &           tol_rho, wght, 
     &           nq, ipol)
#endif
c         write (*,*) "Ex",Ex,Amat(n,1),Cmat(n,1)
c         stop
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
           if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
c
c           Spin alpha:
c
           if (rho(n,R_A).gt.tol_rho) then
              gamma = rgamma(n,G_AA)
c             gamma =    delrho(n,1,1)*delrho(n,1,1) +
c    &                   delrho(n,2,1)*delrho(n,2,1) +
c    &                   delrho(n,3,1)*delrho(n,3,1)  
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
               call nwxc_x_pw6core_p(DPOW,BETA,CPW91,n,1,
     &              rho(n,R_A),gamma,func(n),
     &              tol_rho, wght, 
     &              nq, ipol)     
#else
               call nwxc_x_pw6core(DPOW,BETA,CPW91,n,1,
     &              rho(n,R_A),gamma,func(n),
     &              tol_rho, wght, 
     &              nq, ipol)     
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
               call nwxc_x_pw6core_d2(DPOW,BETA,CPW91,n,1,
     &              rho(n,R_A),gamma,func(n),
     &              tol_rho, wght,
     &              nq, ipol)
#else
               call nwxc_x_pw6core_d3(DPOW,BETA,CPW91,n,1,
     &              rho(n,R_A),gamma,func(n),
     &              tol_rho, wght,
     &              nq, ipol)
#endif
            endif
c
c           Spin beta:
c
            if (rho(n,R_B).gt.tol_rho) then
               gamma = rgamma(n,G_BB)
c              gamma =   delrho(n,1,2)*delrho(n,1,2) +
c    &                   delrho(n,2,2)*delrho(n,2,2) +
c    &                   delrho(n,3,2)*delrho(n,3,2)   
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
               call nwxc_x_pw6core_p(DPOW,BETA,CPW91,n,2,
     &              rho(n,R_B),gamma,func(n),
     &              tol_rho, wght, 
     &              nq, ipol)     
#else
               call nwxc_x_pw6core(DPOW,BETA,CPW91,n,2,
     &              rho(n,R_B),gamma,func(n),
     &              tol_rho, wght, 
     &              nq, ipol)     
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
               call nwxc_x_pw6core_d2(DPOW,BETA,CPW91,n,2,
     &              rho(n,R_B),gamma,func(n),
     &              tol_rho, wght,
     &              nq, ipol)
#else 
               call nwxc_x_pw6core_d3(DPOW,BETA,CPW91,n,2,
     &              rho(n,R_B),gamma,func(n),
     &              tol_rho, wght,
     &              nq, ipol)
#endif
            endif
c        
   20    continue
c
      endif
c
      return
      end




#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_pw6core_p(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,func, tol_rho, wght, nq, ipol)     
#else
      Subroutine nwxc_x_pw6core(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,func, tol_rho, wght, nq, ipol)     
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_pw6core_d2(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,func, tol_rho, wght, nq, ipol)
#else
      Subroutine nwxc_x_pw6core_d3(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,func, tol_rho, wght, nq, ipol)
#endif
c
C$Id$
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
c      
      double precision wght
      integer nq, ipol
      type(nwad_dble)::func  ! value of the functional [output]
      type(nwad_dble)::rho
      integer ispin ! alpha=1; beta=2
c
c     Sampling Matrices for the XC Potential & Energy
c
c     double precision Amat(nq,ipol), Cmat(nq,*)
c
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::DPOW, BETA, CPW91
#else
      double precision DPOW, BETA, CPW91
#endif
#else
      double precision DPOW, BETA, CPW91
#endif
      double precision tol_rho, AX, pi, BETAPW91,big

      Parameter (big=1d4)

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
c     double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
c     double precision rhom23, f2x, d2den,d2num
#endif
c
c References:
c
c
c***************************************************************************

      integer n
      type(nwad_dble)::x, sinhm1,dsinhm1
      type(nwad_dble)::rho13, rho43,  Xa,  Ha, denom, num,  
     &                  fprimex, d1num,d1den,
     &        gamma,fx,x2a,ten6xd,expo,bbx2
      integer D0R,D1G,D2RR,D2RG,D2GG

c
c     sinhm1(x)=log(x+sqrt(1d0+x*x))
c     dsinhm1(x)=1d0/dsqrt(1d0+x*x)
c
      pi=acos(-1.d0)
      BETAPW91=(pi*36.d0)**(-5.d0/3.d0)*5.d0
      AX=-(0.75d0/pi)**(1.d0/3.d0)*1.5d0
      if(ispin.eq.1) then
         D0R=D1_RA
         D1G=D1_GAA
         D2RR=D2_RA_RA
         D2RG=D2_RA_GAA
         D2GG=D2_GAA_GAA
      else
         D0R=D1_RB
         D1G=D1_GBB
         D2RR=D2_RB_RB
         D2RG=D2_RB_GBB
         D2GG=D2_GBB_GBB
      endif

      rho13 = (rho*(ipol/2d0))**(1.d0/3.d0)
      rho43 = rho13**4.0d0
      if (gamma.gt.tol_rho**2)then
         xa = sqrt(gamma)/rho43
         x2a=xa*xa
         ten6xd=Xa**DPOW*1.d-6
         expo=0d0
         if(CPW91*x2a.lt.big) expo=exp(-CPW91*x2a)
         bbx2=(BETA-BETAPW91)*x2a
         Ha = asinh(Xa)
         denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
         num = -BETA*x2a+bbx2*expo+ten6xd
         fx=num*denom
c        d1num=-2.d0*xa*(beta-bbx2*expo*(1d0/x2a-CPW91))+
c    +        ten6xd/xa*dpow
c        d1den=6.d0*beta*(ha + xa*dsinhm1(xa)) -
c    -        ten6xd/ax/xa*dpow
c        fprimex=(d1num - d1den*fx)*denom
      else
         gamma = 0.d0
         Xa = 0.d0
         fx=0.d0
         fprimex=0.d0
         denom=0d0
         d1den=0d0
         expo=0d0
         ten6xd=0d0
         x2a=0d0
         bbx2=0d0
      endif
c     
c     if (lfac)then
c        Ex = Ex + 2d0/ipol*rho43*AX*qwght*fac
c        func = func + 2.d0/ipol*rho43*AX*fac
c        Amat(n,D0R) = Amat(n,D0R) + (4.d0/3.d0)*rho13*AX*fac
c     endif
c     
c     if (nlfac)then
c        Ex = Ex + 2d0/ipol*rho43*fx*qwght*fac
         func = func + 2.d0/ipol*rho43*fx*wght
c        Amat(n,D0R) = Amat(n,D0R) +
c    +        (4.d0/3.d0)*rho13*(fx-xa*fprimex)*wght
c        if (xa.gt.tol_rho)  then
c           Cmat(n,D1G)=Cmat(n,D1G)+
c    +           .5d0*fprimex/sqrt(gamma)*wght
c        endif
c     endif
c     
#ifdef SECOND_DERIV
c     rhom23 = 2d0/ipol*rho13/rho
c     if (lfac)then
c        Amat2(n,D2RR) = Amat2(n,D2RR) +
c    &        (4d0/9d0)*rhom23*Ax*fac
c     endif
c     if (nlfac)then
c     if(gamma.gt.tol_rho**2)then
c        d2num=-2d0*beta +
c    +        2d0*bbx2*expo*
c    *        (1d0/x2a-5d0*CPW91+2d0*CPW91*CPW91*x2a)+
c    +        ten6xd/x2a*dpow*(dpow-1d0)
c        d2den=6.d0*beta*dsinhm1(xa)*(2d0-x2a/(1d0+x2a)) -
c    -        ten6xd/ax/x2a*dpow*(dpow-1d0)
c        f2x = denom*(d2num -fx*d2den - 2d0*d1den*fprimex)
c      else
c        f2x=0d0
c      endif
c        Amat2(n,D2RR) = Amat2(n,D2RR)
c    &        + (4d0/9d0)*rhom23*(fx-xa*fprimex+4d0*x2a*f2x)*wght
c        Cmat2(n,D2RG) = Cmat2(n,D2RG)
c    &        - (4d0/ipol/3d0)*(rhom23**2/rho)*f2x*wght
c        if (xa.gt.tol_rho) then
c           Cmat2(n,D2GG) = Cmat2(n,D2GG)
c    &           - 0.25d0*gamma**(-1.5d0)*(fprimex-xa*f2x)*wght
c        endif
c     endif
#endif
c
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_pw6.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_pw6.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_pw6.F"
#endif
#undef NWAD_PRINT
C>
C> @}
