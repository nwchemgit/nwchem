#ifndef SECOND_DERIV
      Subroutine xc_xpw6(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                  Amat, Cmat, nq, ipol, Ex, qwght,ldew,func,ijzy)     
#else
      Subroutine xc_xpw6_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func,ijzy)
#endif
c
C$Id$
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol,ijzy
      logical lfac, nlfac,ldew
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
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
      double precision DPOW
      double precision BETA,  tol_rho, CPW91

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
c References:
c
c
c***************************************************************************

      integer n
      double precision gamma
c
      if (ijzy.eq.2) then
c
c       parameters for PWB6K
c
        BETA  = 0.00539d0
        CPW91 = 1.7077d0
        DPOW  = 4.0876d0
      else
c
c       parameters for PW6B95
c
        BETA  = 0.00538d0
        CPW91 = 1.7382d0
        DPOW  = 3.8901d0
      endif
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            gamma=0.25d0*gamma
#ifndef SECOND_DERIV
            call xc_xpw6core(DPOW,BETA,CPW91,n,1,
     &           rho(n,1),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Cmat, nq, ipol, Ex, ldew)     
#else 
               call xc_xpw6core_d2(DPOW,BETA,CPW91,n,1,
     &           rho(n,1),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &           ldew)
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
           if (rho(n,1).lt.tol_rho) goto 20
c
c           Spin alpha:
c
           if (rho(n,2).gt.tol_rho) then
              gamma =    delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,3,1)*delrho(n,3,1)  
#ifndef SECOND_DERIV
               call xc_xpw6core(DPOW,BETA,CPW91,n,1,
     &           rho(n,2),gamma,qwght(n),func(n),
     &              tol_rho, fac, lfac, nlfac, 
     &              Amat, Cmat, nq, ipol, Ex, ldew)     
#else
               call xc_xpw6core_d2(DPOW,BETA,CPW91,n,1,
     &           rho(n,2),gamma,qwght(n),func(n),
     &              tol_rho, fac, lfac, nlfac, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         ldew)
#endif
            endif
c
c           Spin beta:
c
            if (rho(n,3).gt.tol_rho) then

            gamma =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2)   
#ifndef SECOND_DERIV
            call xc_xpw6core(DPOW,BETA,CPW91,n,2,
     &           rho(n,3),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Cmat, nq, ipol, Ex, ldew)     
#else 
               call xc_xpw6core_d2(DPOW,BETA,CPW91,n,2,
     &           rho(n,3),gamma,qwght(n),func(n),
     &           tol_rho, fac, lfac, nlfac, 
     &           Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &           ldew)
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
      Subroutine xc_xpw6core(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,qwght,func,
     &     tol_rho, fac, lfac, nlfac,
     &                      Amat, Cmat, nq, ipol, Ex, ldew)     
#else
      Subroutine xc_xpw6core_d2(DPOW,BETA,CPW91,n,ispin,
     .     rho,gamma,qwght,func,
     &     tol_rho, fac, lfac, nlfac,
     &     Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &     ldew)
#endif
c
C$Id$
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func  ! value of the functional [output]
      double precision rho
      double precision qwght
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

      Parameter (big=1d4)

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
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
         D0R=1
         D1G=D1_GAA
         D2RR=D2_RA_RA
         D2RG=D2_RA_GAA
         D2GG=D2_GAA_GAA
      else
         D0R=2
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
      if (lfac)then
         Ex = Ex + 2d0/ipol*rho43*AX*qwght*fac
         if(ldew)func = func + 2.d0/ipol*rho43*AX*fac
         Amat(n,D0R) = Amat(n,D0R) + (4.d0/3.d0)*rho13*AX*fac
      endif
c     
      if (nlfac)then
         Ex = Ex + 2d0/ipol*rho43*fx*qwght*fac
         if(ldew)func = func + 2.d0/ipol*rho43*fx*fac
         Amat(n,D0R) = Amat(n,D0R) +
     +        (4.d0/3.d0)*rho13*(fx-xa*fprimex)*fac
         if (xa.gt.tol_rho)  then
            Cmat(n,D1G)=Cmat(n,D1G)+
     +           .5d0*fprimex/sqrt(gamma)*fac
         endif
      endif
c     
#ifdef SECOND_DERIV
      rhom23 = 2d0/ipol*rho13/rho
      if (lfac)then
         Amat2(n,D2RR) = Amat2(n,D2RR) +
     &        (4d0/9d0)*rhom23*Ax*fac
      endif
      if (nlfac)then
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
     &        + (4d0/9d0)*rhom23*(fx-xa*fprimex+4d0*x2a*f2x)*fac
         Cmat2(n,D2RG) = Cmat2(n,D2RG)
     &        - (4d0/ipol/3d0)*(rhom23**2/rho)*f2x*fac
         if (xa.gt.tol_rho) then
            Cmat2(n,D2GG) = Cmat2(n,D2GG)
     &           - 0.25d0*gamma**(-1.5d0)*(fprimex-xa*f2x)*fac
         endif
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
#include "xc_xpw6.F"
#endif
