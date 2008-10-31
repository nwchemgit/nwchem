c
c     -----------------------------------------------------------------------
c     Uniform electron gas exchange functional for the erfc(r)/r interaction
c     as implemented in the following paper:
c     "A well-tempered density functional theory of electrons in molecules"
c     Ester Livshits & Roi Baer, Phys. Chem. Chem. Phys., 9, 2932 (2007)
c     The other relevant publication is: 
c     R. Baer, D. Neuhauser, Phys. Rev. Lett., 94, 043002 (2005)
c     -----------------------------------------------------------------------
c
#ifndef SECOND_DERIV
      subroutine xc_rodaes_erf(tol_rho, fac, lfac, nlfac, rho, Amat, nq,
     &                    ipol, Ex, qwght, ldew, func)
#else
c     For locations of 2nd derivatives of functionals in array
#include "dft2drv.fh"
      subroutine xc_rodaes_erf_d2(tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, nq, ipol, Ex, qwght, ldew, func)
#endif
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
c
ccase...start
#include "case.fh"
ccase...end
c
      integer nq, ipol, n
      double precision fac, Ex, total
      logical ldew, lfac, nlfac
      double precision func(*) ! value of the functional [output]
      double precision tol_rho
      double precision rho(nq,(ipol*(ipol+1))/2) ! charge density
      double precision qwght(nq)                 ! quadrature weights
      double precision Amat(nq,ipol)             ! partial first derivatives  
      double precision F(nq),RA(nq),RB(nq)
      double precision rhoA, rhoB, rhoTotal, rhoA1, rhoB1
      double precision gamma
      double precision fA, fB, fpA, fpB, fppA, fppB
      double precision EpsX
      double precision EpsXprime
      double precision EpsTwoXprime
c
#ifdef SECOND_DERIV
      double precision Amat2(nq,*)               ! partial second derivatives
#endif
c
c     -----------------------------------------------------------------------
c     Preliminaries
c     -----------------------------------------------------------------------
c
      gamma = cam_omega

c      write(luout,*) "In xc_rodaes_erf"
c      write(luout,*) "cam_omega:",cam_omega
c      write(luout,*) "gamma:",gamma

      do n = 1,nq
         if (ipol.eq.1) then   ! spin-restricted
            rA(n) = rho(n,1)
            rB(n) = 0.d0
         else                  ! spin-unrestricted
            rA(n) = rho(n,2)
            rB(n) = rho(n,3)
         end if
      end do
c
c     -----------------------------------------------------------------------
c     Calculate the first and second derivatives
c     -----------------------------------------------------------------------
c
      total = 0.d0
      do n = 1,nq
         rhoA = rA(n)
         rhoB = rB(n)
         rhoTotal  = rhoA + rhoB   ! total density at point
         if (rhoTotal.gt.tol_rho) then

            if (ipol.eq.1) then    ! spin-restricted
              rhoA1 = rhoA
              rhoB1 = rhoB
            else                   ! spin-unrestricted
              rhoA1 = rhoA*2
              rhoB1 = rhoB*2
            end if

            fA   = EpsX(rhoA1,gamma)
            fB   = EpsX(rhoB1,gamma)
            fpA  = EpsXprime(rhoA1,gamma)
            fpB  = EpsXprime(rhoB1,gamma)

            f(n) = fA * rhoA + fB * rhoB
            Amat(n,1) = Amat(n,1) + (fpA*rhoA1+fA)*fac
            if (ipol.gt.1) then
              Amat(n,2) = Amat(n,2) + (fpB*rhoB1+fB)*fac
            end if

#ifdef SECOND_DERIV
            fppA = EpsTwoXprime(rhoA1,gamma)
            fppB = EpsTwoXprime(rhoB1,gamma)

            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     &        ((fppA*rhoA+fpA)*4)*fac
c
c           Guard against case of no beta electrons, e.g. H atom
            if (ipol.gt.1) then
             if (rho(n,3).gt.tol_rho) then
              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
     &          ((fppB*rhoB+fpB)*4)*fac
             end if
            end if
#endif
            if (ldew) func(n) = func(n) + f(n)*fac 
            total = total + f(n)*qwght(n)
         end if
      end do

      Ex = Ex + total*fac

      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
#include "xc_rodaes_erf.F"
c
c     ---------------------------------------------------------------------------------------
c     Utility functions
c     ---------------------------------------------------------------------------------------
c
c     ---------------------------------------------------------------------------------------
c     Return the value of pi
c     ---------------------------------------------------------------------------------------
c
      double precision function ValueOfPi()
c 
      implicit none
#include "xc_params.fh"
c      
      ValueOfPi = pi          

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Evaluates the actual function
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNL(q)

      implicit none
#include "xc_params.fh"

      double precision q,TwoSqrtPi,OneOverQ,q2,DERF

      OneOverQ = 1.0d0/q
      TwoSqrtPi = 2.0d0*dsqrt(pi) 
      q2 = q**2.0d0

      if (q .lt. 1D-15) then
         HqBNL=1.d0
         return
      end if

      if (q .lt. 0.1d0) then
         HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*DERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the local Fermi vector for the provided density
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiK(den)

      implicit none
#include "xc_params.fh"

      double precision F13, den

      F13 = 1.0D0 / 3.0D0
      FermiK = (3.d0*pi*pi*den)**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the function EpsX at the given density value and gamma
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsX(Rho,gamma)

      implicit none
#include "xc_params.fh"

      double precision  kF,RHO,gamma,Cs
      double precision HqBNL
      double precision FermiK

      if (RHO.le.0D0) then
         EpsX = 0.0D0
         return
      end if

      kF = FermiK(Rho)
      Cs = -3.0D0/(4.0d0*pi)
      EpsX = Cs * kF * HqBNL(gamma/kF)

      return
      end      
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the function
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNLPrime(q)

      implicit none
#include "xc_params.fh"

      double precision q,OneOverQ,q2,q3,DERF

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0
      q3 = q**3.0d0

      if (q .lt. 0.1d0) then
        HqBNLPrime = -4.0d0/3.0d0*(dsqrt(Pi)+2.0d0*q3-3.0d0*q)
        return
      end if

      HqBNLPrime = 4.0d0/3.0d0*(q*(exp(-OneOverQ*OneOverQ)*(2.0d0*q2
     $     -1.0d0)+(3.0d0-2.0d0*q2))-dsqrt(Pi)*DERF(OneOverQ))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the local Fermi vector (it depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiKPrime(den)

      implicit none
#include "xc_params.fh"
   
      double precision F23, den

      F23 = 2.0D0 / 3.0D0
      FermiKPrime = (Pi/(3.0d0*den))**F23

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of q (q=gamma/kf) (it implicitly depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function QPrime(gamma,kF)

      implicit none

      double precision  kF, FermiK2, gamma

      FermiK2 = kF**2.0d0
      QPrime = -gamma/FermiK2

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsXprime(Rho,gamma)

      implicit none
#include "xc_params.fh"

      double precision Rho,gamma
      double precision Cs,kF,CsPrime

      double precision HqBNL
      double precision HqBNLPrime
      double precision QPrime
      double precision FermiK
      double precision FermiKPrime

      kF = FermiK(Rho)
      CsPrime = -3.0D0/(4.0d0*Pi)
      Cs = CsPrime*kF

      if (Rho.le.0d0) then
         EpsXprime = 0.0d0
         return
      end if

      EpsXprime = FermiKPrime(Rho)*(CsPrime*HqBNL(gamma/kF)+
     $     QPrime(gamma,kF)*HqBNLPrime(gamma/kF)*Cs)

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the main function that consititutes the functional
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNLTwoPrime(q)

      implicit none
#include "xc_params.fh"

      double precision q,OneOverQ,q2

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         HqBNLTwoPrime = 4.0d0-8.0d0*q2
         return
      end if

      HqBNLTwoPrime = exp(-OneOverQ*OneOverQ)*(4.0d0+8.0d0*q2)
     $     -8.0d0*q2+4.0d0

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the local Fermi vector
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiKTwoPrime(den)

      implicit none
#include "xc_params.fh"

      double precision F13, den

      F13 = 1.0D0/3.0D0
      FermiKTwoPrime =  -(8.0d0*Pi**2.0d0/(243.0d0*den**5.0d0))**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of q    
c     ---------------------------------------------------------------------------------------
c
      double precision function QTwoPrime(gamma,kF)

      implicit none

      double precision gamma, kF, FermiK3

      FermiK3 = kF**3.0d0
      QTwoPrime = (2.0d0*gamma)/FermiK3

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsTwoXprime(Rho,gamma)

      implicit none
#include "xc_params.fh"

      double precision Rho,gamma
      double precision kF,kFPrim,kFprim2,kF2prim
      double precision q,qprim,qprim2,q2prim
      double precision g,gprim,g2prim
      double precision Cs,CsPrim

      double precision FermiK
      double precision FermiKPrime
      double precision FermiKTwoPrime
      double precision QPrime
      double precision QTwoPrime
      double precision HqBNL
      double precision HqBNLPrime
      double precision HqBNLTwoPrime

      if (Rho.le.0d0) then
         EpsTwoXprime = 0.0d0
         return
      end if

      kF = FermiK(Rho)
      kFPrim = FermiKPrime(Rho)
      kFPrim2=kFPrim**2.0d0
      kF2prim = FermiKTwoPrime(Rho)
      CsPrim = -3.0d0/(4.0d0*Pi)
      Cs = CsPrim * kF
      q = gamma / kF
      qprim = QPrime(gamma,kF)
      Qprim2=qprim**2.0d0
      q2prim = QTwoPrime(gamma,kF)
      g = HqBNL(q)
      gprim = HqBNLPrime(q)
      g2prim = HqBNLTwoPrime(q)

      EpsTwoXprime = 
     $     kFPrim2*(2.0d0*CsPrim*gprim*qprim
     $     +Cs*(QPrim2*g2prim+gprim*Q2Prim))
     $     +kF2Prim*(g*CsPrim+Cs*gprim*qprim)

      return
      end
c
#endif
