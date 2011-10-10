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
      subroutine bnl2007_x(n2ft3d,ispin,dn,x_parameter,gamma,xcp,xce)
      implicit none
      integer n2ft3d
      integer ispin
      real*8  dn(n2ft3d,2)
      real*8  x_parameter,gamma

      real*8  xcp(n2ft3d,2)
      real*8  xce(n2ft3d)

c
c     **** local variables ****
      double precision dncut
      parameter (dncut=1.0d-30)

      integer n
      double precision rhoA1,rhoB1,rhoTotal
      double precision fA, fB, fpA, fpB, fppA, fppB

      double precision TEpsX
      double precision TEpsXprime
      external         TEpsX
      external         TEpsXprime

c     -----------------------------------------------------------------------
c     Calculate the first derivatives
c     -----------------------------------------------------------------------
c
c     **** spin-restricted ****
      if (ispin.eq.1) then   
         do n = 1,n2ft3d
            rhoA1 = dn(n,1) + dn(n,ispin)
            if (rhoA1.gt.dncut) then
               fA    = TEpsX(rhoA1,gamma)
               fpA   = TEpsXprime(rhoA1,gamma)
               xcp(n,1) = xcp(n,1) + x_parameter*(fpA*rhoA1+fA)
               xce(n)   = xce(n)   + x_parameter*fA
            end if
         end do

c     **** spin-unrestricted ****
      else
         do n = 1,n2ft3d
            rhoTotal = dn(n,1) + dn(n,ispin)
            if (rhoTotal.gt.dncut) then
               rhoA1 = dn(n,1)     + dn(n,1)
               rhoB1 = dn(n,ispin) + dn(n,ispin)
               fA    = TEpsX(rhoA1,gamma)
               fB    = TEpsX(rhoB1,gamma)
               fpA   = TEpsXprime(rhoA1,gamma)
               fpB   = TEpsXprime(rhoB1,gamma)
               xcp(n,1) = xcp(n,1) + x_parameter*(fpA*rhoA1+fA)
               xcp(n,2) = xcp(n,2) + x_parameter*(fpB*rhoB1+fB)
               xce(n) = xce(n) 
     >                + x_parameter*( fA*dn(n,1) 
     >                              + fB*dn(n,ispin) )/(rhoTotal)
            end if
         end do
      end if

      return
      end 

c
c     ---------------------------------------------------------------------------------------
c     Evaluates the actual function
c     ---------------------------------------------------------------------------------------
c
      double precision function THqBNL(q)
      implicit none
      double precision q

      double precision Pi
      parameter (Pi = 3.141592653589793d0)
      double precision TwoSqrtPi,OneOverQ,q2
      double precision util_erf
      external util_erf

      OneOverQ = 1.0d0/q
      TwoSqrtPi = 2.0d0*dsqrt(Pi)
      q2 = q**2.0d0

      if (q .lt. 1.0d-15) then
         THqBNL=1.d0
         return
      end if

      if (q .lt. 0.1d0) then
         THqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      THqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*util_erf(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the function
c     ---------------------------------------------------------------------------------------
c
      double precision function THqBNLPrime(q)
      implicit none
      double precision q

      double precision Pi
      parameter (Pi = 3.141592653589793d0)
      double precision OneOverQ,q2,q3
      double precision util_erf
      external         util_erf

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0
      q3 = q**3.0d0

      if (q .lt. 0.1d0) then
        THqBNLPrime = -4.0d0/3.0d0*(dsqrt(Pi)+2.0d0*q3-3.0d0*q)
        return
      end if

      THqBNLPrime = 4.0d0/3.0d0*(q*(exp(-OneOverQ*OneOverQ)*(2.0d0*q2
     $     -1.0d0)+(3.0d0-2.0d0*q2))-dsqrt(Pi)*util_erf(OneOverQ))

      return
      end


c
c     ---------------------------------------------------------------------------------------
c     Calculate the local Fermi vector for the provided density
c     ---------------------------------------------------------------------------------------
c
      double precision function TFermiK(den)
      implicit none
      double precision den

      double precision Pi,F13
      parameter (Pi = 3.141592653589793d0)
      parameter (F13  = 1.0d0/3.0d0)

      TFermiK = (3.d0*Pi*Pi*den)**F13
      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the local Fermi vector (it depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function TFermiKPrime(den)
      implicit none
      double precision  den

      double precision Pi,F23
      parameter (Pi   = 3.141592653589793d0)
      parameter (F23  = 2.0d0/3.0d0)


      TFermiKPrime = (Pi/(3.0d0*den))**F23
      return
      end

c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of q (q=gamma/kf) (it implicitly depends on the density)
c     ---------------------------------------------------------------------------------------
      double precision function TQPrime(gamma,kF)
      implicit none
      double precision  kF, gamma


      TQPrime = -gamma/(kF*kF)
      return
      end


c
c     ---------------------------------------------------------------------------------------
c     Calculate the function EpsX at the given density value and gamma
c     ---------------------------------------------------------------------------------------
c
      double precision function TEpsX(Rho,gamma)
      implicit none
      double precision  Rho,gamma

      double precision Pi
      parameter (Pi = 3.141592653589793d0)
      double precision  kF,Cs
      double precision THqBNL
      double precision TFermiK
      external         THqBNL
      external         TFermiK

      if (Rho.le.0.0d0) then
         TEpsX = 0.0d0
         return
      end if

      kF = TFermiK(Rho)
      Cs = -3.0D0/(4.0d0*Pi)
      TEpsX = Cs * kF * THqBNL(gamma/kF)
      return
      end



c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of TEpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function TEpsXprime(Rho,gamma)
      implicit none
      double precision Rho,gamma

      double precision Pi
      parameter (Pi = 3.141592653589793d0)
      double precision Cs,kF,CsPrime

      double precision THqBNL
      double precision THqBNLPrime
      double precision TQPrime
      double precision TFermiK
      double precision TFermiKPrime
      external         THqBNL
      external         THqBNLPrime
      external         TQPrime
      external         TFermiK
      external         TFermiKPrime


      kF = TFermiK(Rho)
      CsPrime = -3.0D0/(4.0d0*Pi)
      Cs = CsPrime*kF

      if (Rho.le.0d0) then
         TEpsXprime = 0.0d0
         return
      end if

      TEpsXprime = TFermiKPrime(Rho)*(CsPrime*THqBNL(gamma/kF)+
     $     TQPrime(gamma,kF)*THqBNLPrime(gamma/kF)*Cs)
      return
      end

c $Id$
