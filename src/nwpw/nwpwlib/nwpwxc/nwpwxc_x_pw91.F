#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_pw91.F
C> The PW91 exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the PW91 exchange functional
C>
C> Evaluate the PW91 GGA exchange functional [1-3].
C>
C> ### References ###
C>
C> [1] J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
C>     M.R. Pederson, D.J. Singh, C. Fiolhais, 
C>     "Atoms, molecules, solids, and surfaces: Applications of the
C>     generalized gradient approximation for exchange and correlation",
C>     Phys. Rev. B <b>46</b>, 6671-6687 (1992), DOI:
C>     <A HREF="https://doi.org/10.1103/PhysRevB.46.6671">
C>     10.1103/PhysRevB.46.6671</A>.
C>
C> [2] J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
C>     M.R. Pederson, D.J. Singh, C. Fiolhais, 
C>     "Erratum: Atoms, molecules, solids, and surfaces: Applications
C>     of the generalized gradient approximation for exchange and
C>     correlation",
C>     Phys. Rev. B <b>48</b>, 4978-4978 (1993), DOI:
C>     <A HREF="https://doi.org/10.1103/PhysRevB.48.4978.2">
C>     10.1103/PhysRevB.48.4978.2</A>.
C>
C> [3] Y. Zhao, D.G. Truhlar,
C>     "Design of density functionals that are broadly accurate for
C>     thermochemistry, thermochemical kinetics, and nonbonded
C>     interactions", J. Phys. Chem. A <b>109</b>, 5656-5667 (2005)
C>     DOI:
C>     <a href="https://doi.org/10.1021/jp050536c">
C>     10.1021/jp050536c</a>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_pw91(param, tol_rho, ipol, nq, wght,
     &                       rho, rgamma, func, Amat, Cmat)
#else
      Subroutine nwpwxc_x_pw91_d2(param, tol_rho, ipol, nq, wght, rho,
     &                          rgamma, func, Amat, Amat2, Cmat, Cmat2)
#endif
c
C$Id: xc_xmpw91.F 19696 2010-10-29 16:53:42Z d3y133 $
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
c     Input and other parameters
c
      double precision param(*) !< [Input] Parameters of the functional
                                !< (see Eq.(5) of [3])
                                !< - param(1): \f$ d \f$ 
                                !< - param(2): \f$ b \f$ 
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
      double precision func(nq)     !< [Output] The value of the functional
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
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
      double precision DPOW
      double precision BETA

!      Parameter (BETA = 0.0042D0, CPW91=1.6455D0,DPOW=4) ! pw91 paper
!      Parameter (BETA = 0.0046D0, CPW91=1.6455D0,DPOW=3.73D0)!mpw91 paper
c     from Zhao Truhlar, J. Phys. Chem. A 109, 5656, 2005
!     Parameter (BETA = 0.00426D0, DPOW=3.72D0)!errata 
C
C     C. Adamo confirmed that there is a typo in the JCP paper
c     BETA is 0.00426 instead of 0.0046
C     adamo@ext.jussieu.fr
C
c***************************************************************************

      integer n
      double precision gamma
c
      DPOW = param(1)
      BETA = param(2)
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
#ifndef SECOND_DERIV
            call nwpwxc_x_pw91core(DPOW,BETA,n,1,
     &           rho(n,1),gamma,func(n),
     &           tol_rho, wght, 
     &           Amat, Cmat, nq, ipol)
#else 
            call nwpwxc_x_pw91core_d2(DPOW,BETA,n,1,
     &           rho(n,1),gamma,func(n),
     &           tol_rho, wght,
     &           Amat, Amat2, Cmat, Cmat2, nq, ipol)
#endif

   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
           if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
c
c          Spin alpha:
c
           if (rho(n,R_A).gt.tol_rho) then
             gamma = rgamma(n,G_AA)
c            gamma =    delrho(n,1,1)*delrho(n,1,1) +
c    &                  delrho(n,2,1)*delrho(n,2,1) +
c    &                  delrho(n,3,1)*delrho(n,3,1)  
#ifndef SECOND_DERIV
             call nwpwxc_x_pw91core(DPOW,BETA,n,1,
     &            rho(n,R_A),gamma,func(n),
     &            tol_rho, wght, Amat, Cmat, nq, ipol)
#else
             call nwpwxc_x_pw91core_d2(DPOW,BETA,n,1,
     &            rho(n,R_A),gamma,func(n),
     &            tol_rho, wght, Amat, Amat2, Cmat, Cmat2, nq, ipol)
#endif
           endif
c
c          Spin beta:
c
           if (rho(n,R_B).gt.tol_rho) then
             gamma = rgamma(n,G_BB)
c            gamma =   delrho(n,1,2)*delrho(n,1,2) +
c    &                 delrho(n,2,2)*delrho(n,2,2) +
c    &                 delrho(n,3,2)*delrho(n,3,2)   
#ifndef SECOND_DERIV
             call nwpwxc_x_pw91core(DPOW,BETA,n,2,
     &            rho(n,R_B),gamma,func(n),
     &            tol_rho, wght, Amat, Cmat, nq, ipol)     
#else 
             call nwpwxc_x_pw91core_d2(DPOW,BETA,n,2,
     &            rho(n,R_B),gamma,func(n),
     &            tol_rho, wght, Amat, Amat2, Cmat, Cmat2, nq, ipol)
#endif
           endif
c        
   20    continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_pw91core(DPOW,BETA,n,ispin,
     &                           rho,gamma,func,
     &                           tol_rho, wght,
     &                           Amat, Cmat, nq, ipol)
#else
      Subroutine nwpwxc_x_pw91core_d2(DPOW,BETA,n,ispin,
     &                              rho,gamma,func,
     &                              tol_rho, wght,
     &                              Amat, Amat2, Cmat, Cmat2, nq, ipol)
#endif
c
C$Id: xc_xmpw91.F 19696 2010-10-29 16:53:42Z d3y133 $
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
      double precision wght
      integer nq, ipol
      double precision func  ! value of the functional [output]
      double precision rho
      integer ispin ! alpha=1; beta=2
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
      double precision DPOW
      double precision BETA, tol_rho, AX, pi,
     &     CPW91,BETAPW91,big

      Parameter (CPW91=1.6455D0,big=1d4)

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*), Cmat2(nq,*)
      double precision rhom23, f2x, d2den,d2num
#endif
c
c References:
c
c
c***************************************************************************

      integer n
      double precision x, sinhm1,dsinhm1
      double precision rho13, rho43,  Xa,  Ha, denom, num,  
     &                  fprimex, d1num,d1den,
     &        gamma,fx,x2a,ten6xd,expo,bbx2
      integer D0R,D1G,D2RR,D2RG,D2GG

c
      sinhm1(x)=log(x+dsqrt(1d0+x*x))
      dsinhm1(x)=1d0/dsqrt(1d0+x*x)
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

      rho13 = (rho*ipol/2d0)**(1.d0/3.d0)
      rho43 = rho13**4 
      if (gamma.gt.tol_rho**2)then
         xa = sqrt(gamma)/rho43
         x2a=xa*xa
         ten6xd=Xa**DPOW*1.d-6
         expo=0d0
         if(CPW91*x2a.lt.big) expo=exp(-CPW91*x2a)
         bbx2=(BETA-BETAPW91)*x2a
         Ha = sinhm1(Xa)
         denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
         num = -BETA*x2a+bbx2*expo+ten6xd
         fx=num*denom
         d1num=-2.d0*xa*(beta-bbx2*expo*(1d0/x2a-CPW91))+
     +        ten6xd/xa*dpow
         d1den=6.d0*beta*(ha + xa*dsinhm1(xa)) -
     -        ten6xd/ax/xa*dpow
         fprimex=(d1num - d1den*fx)*denom
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
      func = func + 2.d0/ipol*rho43*AX*wght
      Amat(n,D0R) = Amat(n,D0R) + (4.d0/3.d0)*rho13*AX*wght
c     
      func = func + 2.d0/ipol*rho43*fx*wght
      Amat(n,D0R) = Amat(n,D0R) +
     +    (4.d0/3.d0)*rho13*(fx-xa*fprimex)*wght
      if (xa.gt.tol_rho)  then
            Cmat(n,D1G)=Cmat(n,D1G)+
     +           .5d0*fprimex/sqrt(gamma)*wght
      endif
c     
#ifdef SECOND_DERIV
      rhom23 = 2d0/ipol*rho13/rho
      Amat2(n,D2RR) = Amat2(n,D2RR) +
     &        (4d0/9d0)*rhom23*Ax*wght
      if(gamma.gt.tol_rho**2)then
         d2num=-2d0*beta +
     +        2d0*bbx2*expo*
     *        (1d0/x2a-5d0*CPW91+2d0*CPW91*CPW91*x2a)+
     +        ten6xd/x2a*dpow*(dpow-1d0)
         d2den=6.d0*beta*dsinhm1(xa)*(2d0-x2a/(1d0+x2a)) -
     -        ten6xd/ax/x2a*dpow*(dpow-1d0)
         f2x = denom*(d2num -fx*d2den - 2d0*d1den*fprimex)
      else
            f2x=0d0
      endif	
      Amat2(n,D2RR) = Amat2(n,D2RR)
     &              + (4d0/9d0)*rhom23*(fx-xa*fprimex+4d0*x2a*f2x)*wght
      Cmat2(n,D2RG) = Cmat2(n,D2RG)
     &              - (4d0/ipol/3d0)*(rhom23**2/rho)*f2x*wght
      if (xa.gt.tol_rho) then
         Cmat2(n,D2GG) = Cmat2(n,D2GG)
     &                 - 0.25d0*gamma**(-1.5d0)*(fprimex-xa*f2x)*wght
      endif
#endif
c
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_x_pw91.F"
#endif
C>
C> @}
