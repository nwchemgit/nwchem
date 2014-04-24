#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_bnl.F
C> The BNL LDA exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Implementation of the BNL LDA exchange functional
C>
C> Evaluates the BNL LDA exchange functional [1,2].
C>
C> ### References ###
C>
C> [1] E. Livshits, R. Baer,
C> "A well-tempered density functional theory of electrons in
C> molecules",
C> Phys. Chem. Chem. Phys. <b>9</b>, 2932-2941 (2007), DOI:
C> <a href="http://dx.doi.org/10.1039/B617919C">10.1039/B617919C</a>.
C>
C> [2] R. Baer, D. Neuhauser,
C> "Density functional theory with correct long-range asymptotic
C> behavior",
C> Phys. Rev. Lett. <b>94</b>, 043002 (2005), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.94.043002">
C> 10.1103/PhysRevLett.94.043002</a>.
C>
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_x_bnl(param,tol_rho, ipol, nq, wght, rho, 
     &                      func, Amat)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c     For locations of 2nd derivatives of functionals in array
      subroutine nwxc_x_bnl_d2(param,tol_rho, ipol, nq, wght, rho,
     &                         func, Amat, Amat2)
#else
      subroutine nwxc_x_bnl_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         func, Amat, Amat2, Amat3)
#endif
c
      implicit none
c
#include "nwxc_param.fh"
c
      double precision param(1) !< [Input] Parameters of the functional
                                !< - param(1): \f$ \omega_{CAM} \f$
      double precision tol_rho  !< [Input] The lower limit on the density
      integer nq                !< [Input] The number of points
      integer ipol              !< [Input] The number of spin channels
      double precision wght     !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     The Exchange Energy Functional
c
      double precision func(nq)  !< [Output] The value of the functional
c
c     Partial First Derivatives of the Exchange Energy Functional
c
      double precision Amat(nq,*) !< [Output] 1st order partial derivatives
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*) !< [Output] 2nd order partial derivatives
#endif
#if defined(THIRD_DERIV)
c
c     Partial Third Order Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,*) !< [Output] 3rd order partial derivatives
#endif
c
      double precision F(nq),RA(nq),RB(nq)
      double precision rhoA, rhoB, rhoTotal, rhoA1, rhoB1
      double precision gamma
      double precision fA, fB, fpA, fpB, fppA, fppB, fpppA, fpppB
      double precision nwxc_EpsX
      double precision nwxc_EpsXprime
      double precision nwxc_EpsTwoXprime
      double precision nwxc_EpsThreeXprime
      integer          n
c
c     -----------------------------------------------------------------------
c     Preliminaries
c     -----------------------------------------------------------------------
c
      gamma = param(1)
c
      do n = 1,nq
         if (ipol.eq.1) then   ! spin-restricted
            rA(n) = rho(n,R_T)
            rB(n) = 0.d0
         else                  ! spin-unrestricted
            rA(n) = rho(n,R_A)
            rB(n) = rho(n,R_B)
         end if
      end do
c
c     -----------------------------------------------------------------------
c     Calculate the first and second derivatives
c     -----------------------------------------------------------------------
c
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

            fA   = nwxc_EpsX(rhoA1,gamma)
            fB   = nwxc_EpsX(rhoB1,gamma)
            fpA  = nwxc_EpsXprime(rhoA1,gamma)
            fpB  = nwxc_EpsXprime(rhoB1,gamma)

            f(n) = fA * rhoA + fB * rhoB
            Amat(n,D1_RA) = Amat(n,D1_RA) + (fpA*rhoA1+fA)*wght
            if (ipol.gt.1) then
              Amat(n,D1_RB) = Amat(n,D1_RB) + (fpB*rhoB1+fB)*wght
            end if

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c           Guard against case of no beta electrons, e.g. H atom
            if (ipol.eq.1) then
              if (rhoA1.gt.tol_rho) then
                fppA = nwxc_EpsTwoXprime(rhoA1,gamma)
                Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     &            ((fppA*rhoA+2*fpA)*2)*wght
              endif
            else
              if (rhoA1.gt.tol_rho) then
                fppA = nwxc_EpsTwoXprime(rhoA1,gamma)
                Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
     &            ((fppA*rhoA+fpA)*4)*wght
              endif
              if (rhoB1.gt.tol_rho) then
                fppB = nwxc_EpsTwoXprime(rhoB1,gamma)
                Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
     &            ((fppB*rhoB+fpB)*4)*wght
              end if
            end if
#endif
#if defined(THIRD_DERIV)
c
c           Guard against case of no beta electrons, e.g. H atom
            if (ipol.eq.1) then
              if (rhoA1.gt.tol_rho) then
                fpppA = nwxc_EpsThreeXprime(rhoA1,gamma)
                Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     &               + ( fpppA*rhoA + 3.0d0*fppA )*wght*4.0d0
              endif
            else
              if (rhoA1.gt.tol_rho) then
                fpppA = nwxc_EpsThreeXprime(rhoA1,gamma)
                Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     &               + ( fpppA*rhoA + 1.5d0*fppA )*wght*8.0d0
              endif
              if (rhoB1.gt.tol_rho) then
                fpppB = nwxc_EpsThreeXprime(rhoB1,gamma)
                Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     &               + ( fpppB*rhoB + 1.5d0*fppB )*wght*8.0d0
              end if
            end if
#endif
            func(n) = func(n) + f(n)*wght 
         end if
      end do

      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
#include "nwxc_x_bnl.F"
c
c     ---------------------------------------------------------------------------------------
c     Utility functions
c     ---------------------------------------------------------------------------------------
c
c     ---------------------------------------------------------------------------------------
c     Return the value of pi
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_ValueOfPi()
c 
      implicit none
c      
      nwxc_ValueOfPi = acos(-1.0d0)          

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Evaluates the actual function
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_HqBNL(q)

      implicit none

      double precision q,TwoSqrtPi,OneOverQ,q2,DERF
      double precision nwxc_ValueOfPi

      OneOverQ = 1.0d0/q
      TwoSqrtPi = 2.0d0*dsqrt(nwxc_ValueOfPi()) 
      q2 = q**2.0d0

      if (q .lt. 1D-15) then
         nwxc_HqBNL=1.d0
         return
      end if

      if (q .lt. 0.1d0) then
         nwxc_HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      nwxc_HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*DERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the local Fermi vector for the provided density
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_FermiK(den)

      implicit none

      double precision F13, den, nwxc_ValueOfPi

      F13 = 1.0D0 / 3.0D0
      nwxc_FermiK = (3.d0*nwxc_ValueOfPi()*nwxc_ValueOfPi()*den)**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the function EpsX at the given density value and gamma
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_EpsX(Rho,gamma)

      implicit none

      double precision  kF,RHO,gamma,Cs
      double precision nwxc_HqBNL
      double precision nwxc_FermiK
      double precision nwxc_ValueOfPi

      if (RHO.le.0D0) then
         nwxc_EpsX = 0.0D0
         return
      end if

      kF = nwxc_FermiK(Rho)
      Cs = -3.0D0/(4.0d0*nwxc_ValueOfPi())
      nwxc_EpsX = Cs * kF * nwxc_HqBNL(gamma/kF)

      return
      end      
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the function
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_HqBNLPrime(q)

      implicit none

      double precision q,OneOverQ,q2,q3,DERF
      double precision nwxc_ValueOfPi

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0
      q3 = q**3.0d0

      if (q .lt. 0.1d0) then
        nwxc_HqBNLPrime = -4.0d0/3.0d0*
     &                  (dsqrt(nwxc_ValueOfPi())+2.0d0*q3-3.0d0*q)
        return
      end if

      nwxc_HqBNLPrime = 4.0d0/3.0d0*(q*(exp(-OneOverQ*OneOverQ)*
     $     (2.0d0*q2-1.0d0)+(3.0d0-2.0d0*q2))-dsqrt(nwxc_ValueOfPi())*
     $     DERF(OneOverQ))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the local Fermi vector (it depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_FermiKPrime(den)

      implicit none
   
      double precision F23, den, nwxc_ValueOfPi

      F23 = 2.0D0 / 3.0D0
      nwxc_FermiKPrime = (nwxc_ValueOfPi()/(3.0d0*den))**F23

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of q (q=gamma/kf) (it implicitly depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_QPrime(gamma,kF)

      implicit none

      double precision  kF, FermiK2, gamma

      FermiK2 = kF**2.0d0
      nwxc_QPrime = -gamma/FermiK2

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_EpsXprime(Rho,gamma)

      implicit none

      double precision Rho,gamma
      double precision Cs,kF,CsPrime

      double precision nwxc_HqBNL
      double precision nwxc_HqBNLPrime
      double precision nwxc_QPrime
      double precision nwxc_FermiK
      double precision nwxc_FermiKPrime
      double precision nwxc_ValueOfPi

      kF = nwxc_FermiK(Rho)
      CsPrime = -3.0D0/(4.0d0*nwxc_ValueOfPi())
      Cs = CsPrime*kF

      if (Rho.le.0d0) then
         nwxc_EpsXprime = 0.0d0
         return
      end if

      nwxc_EpsXprime = nwxc_FermiKPrime(Rho)*
     $     (CsPrime*nwxc_HqBNL(gamma/kF)+
     $     nwxc_QPrime(gamma,kF)*nwxc_HqBNLPrime(gamma/kF)*Cs)

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the main function that consititutes the functional
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_HqBNLTwoPrime(q)

      implicit none

      double precision q,OneOverQ,q2

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         nwxc_HqBNLTwoPrime = 4.0d0-8.0d0*q2
         return
      end if

      nwxc_HqBNLTwoPrime = exp(-OneOverQ*OneOverQ)*(4.0d0+8.0d0*q2)
     $     -8.0d0*q2+4.0d0

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the local Fermi vector
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_FermiKTwoPrime(den)

      implicit none

      double precision F13, den, nwxc_ValueOfPi

      F13 = 1.0D0/3.0D0
      nwxc_FermiKTwoPrime =  -(8.0d0*nwxc_ValueOfPi()**2.0d0/
     $                        (243.0d0*den**5.0d0))**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of q    
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_QTwoPrime(gamma,kF)

      implicit none

      double precision gamma, kF, FermiK3

      FermiK3 = kF**3.0d0
      nwxc_QTwoPrime = (2.0d0*gamma)/FermiK3

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_EpsTwoXprime(Rho,gamma)

      implicit none

      double precision Rho,gamma
      double precision kF,kFPrim,kFprim2,kF2prim
      double precision q,qprim,qprim2,q2prim
      double precision g,gprim,g2prim
      double precision Cs,CsPrim

      double precision nwxc_FermiK
      double precision nwxc_FermiKPrime
      double precision nwxc_FermiKTwoPrime
      double precision nwxc_QPrime
      double precision nwxc_QTwoPrime
      double precision nwxc_HqBNL
      double precision nwxc_HqBNLPrime
      double precision nwxc_HqBNLTwoPrime
      double precision nwxc_ValueOfPi

      if (Rho.le.0d0) then
         nwxc_EpsTwoXprime = 0.0d0
         return
      end if

      kF = nwxc_FermiK(Rho)
      kFPrim = nwxc_FermiKPrime(Rho)
      kFPrim2=kFPrim**2.0d0
      kF2prim = nwxc_FermiKTwoPrime(Rho)
      CsPrim = -3.0d0/(4.0d0*nwxc_ValueOfPi())
      Cs = CsPrim * kF
      q = gamma / kF
      qprim = nwxc_QPrime(gamma,kF)
      Qprim2=qprim**2.0d0
      q2prim = nwxc_QTwoPrime(gamma,kF)
      g = nwxc_HqBNL(q)
      gprim = nwxc_HqBNLPrime(q)
      g2prim = nwxc_HqBNLTwoPrime(q)

      nwxc_EpsTwoXprime = 
     $     kFPrim2*(2.0d0*CsPrim*gprim*qprim
     $     +Cs*(QPrim2*g2prim+gprim*Q2Prim))
     $     +kF2Prim*(g*CsPrim+Cs*gprim*qprim)

      return
      end
c
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
#include "nwxc_x_bnl.F"
c
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of the main function that consititutes the functional
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_HqBNLThreePrime(q)

      implicit none

      double precision q,OneOverQ
      double precision q2, q3, q4

      OneOverQ = 1.0d0/q
      q2 = q*q
      q3 = q2*q
      q4 = q3*q

      if (q .lt. 0.1d0) then
         nwxc_HqBNLThreePrime = -16.0d0*q
         return
      end if

      nwxc_HqBNLThreePrime = 8.0d0*( exp(-OneOverQ*OneOverQ)
     1                        + 2.0d0*q2*exp(-OneOverQ*OneOverQ)
     2                        - 2.0d0*q4
     3                        + 2.0d0*q4*exp(-OneOverQ*OneOverQ) ) / q3

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of the local Fermi vector
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_FermiKThreePrime(den)

      implicit none

      double precision F13, den
      double precision nwxc_ValueOfPi

      F13 = 1.0D0/3.0D0
      nwxc_FermiKThreePrime =  (10.0d0/9.0d0)*
     1                (nwxc_ValueOfPi()**2.0d0/(9.0d0*den**8.0d0))**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of q    
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_QThreePrime(gamma,kF)

      implicit none

      double precision gamma, kF, FermiK4

      FermiK4 = kF**4.0d0
      nwxc_QThreePrime = -(6.0d0*gamma)/FermiK4

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function nwxc_EpsThreeXprime(Rho,gamma)

      implicit none

      double precision Rho, gamma
c Fermi wavevector stuff
      double precision kF, kFPrim, kF2prim, kF3prim
      double precision kFPrim2, kFPrim3 
c q stuff
      double precision q, qprim, q2prim, q3prim
      double precision qprim2, qprim3 
c H(q) stuff
      double precision g, gprim, g2prim, g3prim
      double precision Cs

      double precision nwxc_FermiK
      double precision nwxc_FermiKPrime
      double precision nwxc_FermiKTwoPrime
      double precision nwxc_FermiKThreePrime
      double precision nwxc_QPrime
      double precision nwxc_QTwoPrime
      double precision nwxc_QThreePrime
      double precision nwxc_HqBNL
      double precision nwxc_HqBNLPrime
      double precision nwxc_HqBNLTwoPrime
      double precision nwxc_HqBNLThreePrime
      double precision nwxc_ValueOfPi

      if (Rho.le.0d0) then
         nwxc_EpsThreeXprime = 0.0d0
         return
      end if

      kF = nwxc_FermiK(Rho)
      kFPrim = nwxc_FermiKPrime(Rho)
      kF2Prim = nwxc_FermiKTwoPrime(Rho)
      kF3Prim = nwxc_FermiKThreePrime(Rho)
c
      kFPrim2 = kFPrim**2.0d0
      kFPrim3 = kFPrim**3.0d0
c
      Cs = -3.0d0/(4.0d0*nwxc_ValueOfPi())
c
      q = gamma / kF
      qprim = nwxc_QPrime(gamma,kF)
      q2prim = nwxc_QTwoPrime(gamma,kF)
      q3prim = nwxc_QThreePrime(gamma,kF)
c
      qprim2 = qprim**2.0d0
      qprim3 = qprim**3.0d0
c
      g = nwxc_HqBNL(q)
      gprim = nwxc_HqBNLPrime(q)
      g2prim = nwxc_HqBNLTwoPrime(q)
      g3prim = nwxc_HqBNLThreePrime(q)

      nwxc_EpsThreeXprime = Cs*kFPrim3*( 3.0d0*qprim2*g2prim
     1                            + 3.0d0*kF*qprim*g2prim*q2prim
     2                            + kF*qprim3*g3prim
     3                            + gprim*( 3.0d0*q2prim
     4                                    + kF*q3prim ) )
     5               + 3.0d0*Cs*kFPrim*kF2Prim*( kF*qprim2*g2prim
     6                                         + gprim*( 2.0d0*qprim
     7                                                 + kF*q2prim ) )
     8               + Cs*kF3Prim*( g + kF*gprim*qprim )

      return
      end
#endif
c $Id$
C> @}
