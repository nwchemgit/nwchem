#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_bnl.F
C> The BNL LDA exchange functional
C>
C> @}
#endif
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
C> <a href="https://doi.org/10.1039/B617919C">10.1039/B617919C</a>.
C>
C> [2] R. Baer, D. Neuhauser,
C> "Density functional theory with correct long-range asymptotic
C> behavior",
C> Phys. Rev. Lett. <b>94</b>, 043002 (2005), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.94.043002">
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
#if defined(NWAD_PRINT)
      subroutine nwxc_x_bnl_p(param,tol_rho, ipol, nq, wght, rho, 
     &                        func)
#else
      subroutine nwxc_x_bnl(param,tol_rho, ipol, nq, wght, rho, 
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c     For locations of 2nd derivatives of functionals in array
      subroutine nwxc_x_bnl_d2(param,tol_rho, ipol, nq, wght, rho,
     &                         func)
#else
      subroutine nwxc_x_bnl_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         func)
#endif
c
#include "nwad.fh"
      implicit none
c
#include "nwxc_param.fh"
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(1) !< [Input] Parameters of the functional
      type(nwad_dble)::gamma
#else
      double precision param(1) !< [Input] Parameters of the functional
      double precision gamma
#endif
#else
      double precision param(1) !< [Input] Parameters of the functional
                                !< - param(1): \f$ \omega_{CAM} \f$
      double precision gamma
#endif
      double precision tol_rho  !< [Input] The lower limit on the density
      integer nq                !< [Input] The number of points
      integer ipol              !< [Input] The number of spin channels
      double precision wght     !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     The Exchange Energy Functional
c
      type(nwad_dble)::func(nq)  !< [Output] The value of the functional
c
c     Partial First Derivatives of the Exchange Energy Functional
c
c     double precision Amat(nq,*) !< [Output] 1st order partial derivatives
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
c     double precision Amat2(nq,*) !< [Output] 2nd order partial derivatives
#endif
#if defined(THIRD_DERIV)
c
c     Partial Third Order Derivatives of the Exchange Energy Functional
c
c     double precision Amat3(nq,*) !< [Output] 3rd order partial derivatives
#endif
c
      type(nwad_dble)::F(nq),RA(nq),RB(nq)
      type(nwad_dble)::rhoA, rhoB, rhoTotal, rhoA1, rhoB1
      type(nwad_dble)::fA, fB, fpA, fpB, fppA, fppB, fpppA, fpppB
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      type(nwad_dble)::nwxc_EpsX_p
#else
      type(nwad_dble)::nwxc_EpsX
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::nwxc_EpsX_d2
#else
      type(nwad_dble)::nwxc_EpsX_d3
#endif
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
            rA(n) = 0.0d0
            rB(n) = 0.0d0
            if (rho(n,R_A).gt.0.5d0*tol_rho) then
              rA(n) = rho(n,R_A)
            endif
            if (rho(n,R_B).gt.0.5d0*tol_rho) then
              rB(n) = rho(n,R_B)
            endif
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
              rhoA1 = rhoA*2.0d0
              rhoB1 = rhoB*2.0d0
            end if

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
            fA   = nwxc_EpsX_p(rhoA1,gamma)
            fB   = nwxc_EpsX_p(rhoB1,gamma)
#else
            fA   = nwxc_EpsX(rhoA1,gamma)
            fB   = nwxc_EpsX(rhoB1,gamma)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            fA   = nwxc_EpsX_d2(rhoA1,gamma)
            fB   = nwxc_EpsX_d2(rhoB1,gamma)
#else
            fA   = nwxc_EpsX_d3(rhoA1,gamma)
            fB   = nwxc_EpsX_d3(rhoB1,gamma)
#endif
c           fpA  = nwxc_EpsXprime(rhoA1,gamma)
c           fpB  = nwxc_EpsXprime(rhoB1,gamma)

            f(n) = fA * rhoA + fB * rhoB
c           Amat(n,D1_RA) = Amat(n,D1_RA) + (fpA*rhoA1+fA)*wght
c           if (ipol.gt.1) then
c             Amat(n,D1_RB) = Amat(n,D1_RB) + (fpB*rhoB1+fB)*wght
c           end if

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c           Guard against case of no beta electrons, e.g. H atom
c           if (ipol.eq.1) then
c             if (rhoA1.gt.tol_rho) then
c               fppA = nwxc_EpsTwoXprime(rhoA1,gamma)
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
c    &            ((fppA*rhoA+2*fpA)*2)*wght
c             endif
c           else
c             if (rhoA1.gt.tol_rho) then
c               fppA = nwxc_EpsTwoXprime(rhoA1,gamma)
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
c    &            ((fppA*rhoA+fpA)*4)*wght
c             endif
c             if (rhoB1.gt.tol_rho) then
c               fppB = nwxc_EpsTwoXprime(rhoB1,gamma)
c               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
c    &            ((fppB*rhoB+fpB)*4)*wght
c             end if
c           end if
#endif
#if defined(THIRD_DERIV)
c
c           Guard against case of no beta electrons, e.g. H atom
c           if (ipol.eq.1) then
c             if (rhoA1.gt.tol_rho) then
c               fpppA = nwxc_EpsThreeXprime(rhoA1,gamma)
c               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    &               + ( fpppA*rhoA + 3.0d0*fppA )*wght*4.0d0
c             endif
c           else
c             if (rhoA1.gt.tol_rho) then
c               fpppA = nwxc_EpsThreeXprime(rhoA1,gamma)
c               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    &               + ( fpppA*rhoA + 1.5d0*fppA )*wght*8.0d0
c             endif
c             if (rhoB1.gt.tol_rho) then
c               fpppB = nwxc_EpsThreeXprime(rhoB1,gamma)
c               Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
c    &               + ( fpppB*rhoB + 1.5d0*fppB )*wght*8.0d0
c             end if
c           end if
#endif
            func(n) = func(n) + f(n)*wght 
         end if
      end do

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Utility functions
c     ---------------------------------------------------------------------------------------
c
c     ---------------------------------------------------------------------------------------
c     Return the value of pi
c     ---------------------------------------------------------------------------------------
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
      double precision function nwxc_ValueOfPi()
c 
      implicit none
c      
      nwxc_ValueOfPi = acos(-1.0d0)          

      return
      end
#endif
#endif
c
c     ---------------------------------------------------------------------------------------
c     Evaluates the actual function
c     ---------------------------------------------------------------------------------------
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      function nwxc_HqBNL_p(q)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_HqBNL_p
      type(nwad_dble)::q,OneOverQ,q2
      double precision nwxc_ValueOfPi,TwoSqrtPi

      TwoSqrtPi = 2.0d0*dsqrt(nwxc_ValueOfPi()) 

      if (q .lt. 1D-15) then
         nwxc_HqBNL_p=1.d0
         return
      end if
      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

c     if (q .lt. 0.1d0) then
c        nwxc_HqBNL_p=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
c        return
c     end if

      nwxc_HqBNL_p=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*ERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
#else
      function nwxc_HqBNL(q)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_HqBNL
      type(nwad_dble)::q,OneOverQ,q2
      double precision nwxc_ValueOfPi,TwoSqrtPi

      TwoSqrtPi = 2.0d0*dsqrt(nwxc_ValueOfPi()) 

      if (q .lt. 1D-15) then
         nwxc_HqBNL=1.d0
         return
      end if
      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         nwxc_HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      nwxc_HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*ERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      function nwxc_HqBNL_d2(q)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_HqBNL_d2
      type(nwad_dble)::q,OneOverQ,q2
      double precision nwxc_ValueOfPi,TwoSqrtPi

      TwoSqrtPi = 2.0d0*dsqrt(nwxc_ValueOfPi()) 

      if (q .lt. 1D-15) then
         nwxc_HqBNL_d2=1.d0
         return
      end if
      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         nwxc_HqBNL_d2=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      nwxc_HqBNL_d2=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*ERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
#else
      function nwxc_HqBNL_d3(q)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_HqBNL_d3
      type(nwad_dble)::q,OneOverQ,q2
      double precision nwxc_ValueOfPi,TwoSqrtPi

      TwoSqrtPi = 2.0d0*dsqrt(nwxc_ValueOfPi()) 

      if (q .lt. 1D-15) then
         nwxc_HqBNL_d3=1.d0
         return
      end if
      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         nwxc_HqBNL_d3=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      nwxc_HqBNL_d3=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*ERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
#endif
c
c     ---------------------------------------------------------------------------------------
c     Calculate the local Fermi vector for the provided density
c     ---------------------------------------------------------------------------------------
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      function nwxc_FermiK_p(den)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_FermiK_p
      type(nwad_dble):: den
      double precision  F13, nwxc_ValueOfPi

      F13 = 1.0D0 / 3.0D0
      nwxc_FermiK_p = (3.d0*nwxc_ValueOfPi()*nwxc_ValueOfPi()*den)**F13

      return
      end
#else
      function nwxc_FermiK(den)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_FermiK
      type(nwad_dble):: den
      double precision  F13, nwxc_ValueOfPi

      F13 = 1.0D0 / 3.0D0
      nwxc_FermiK = (3.d0*nwxc_ValueOfPi()*nwxc_ValueOfPi()*den)**F13

      return
      end
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      function nwxc_FermiK_d2(den)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_FermiK_d2
      type(nwad_dble):: den
      double precision  F13, nwxc_ValueOfPi

      F13 = 1.0D0 / 3.0D0
      nwxc_FermiK_d2 = (3.d0*nwxc_ValueOfPi()*nwxc_ValueOfPi()*den)**F13

      return
      end
#else
      function nwxc_FermiK_d3(den)

#include "nwad.fh"

      implicit none

      type(nwad_dble):: nwxc_FermiK_d3
      type(nwad_dble):: den
      double precision  F13, nwxc_ValueOfPi

      F13 = 1.0D0 / 3.0D0
      nwxc_FermiK_d3 = (3.d0*nwxc_ValueOfPi()*nwxc_ValueOfPi()*den)**F13

      return
      end
#endif
c
c     ---------------------------------------------------------------------------------------
c     Calculate the function EpsX at the given density value and gamma
c     ---------------------------------------------------------------------------------------
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      function nwxc_EpsX_p(Rho,gamma)

#include "nwad.fh"

      implicit none

      type(nwad_dble)::nwxc_EpsX_p
      type(nwad_dble)::kF,RHO
      type(nwad_dble)::gamma
      type(nwad_dble)::nwxc_HqBNL_p
      type(nwad_dble)::nwxc_FermiK_p
      double precision Cs
      double precision nwxc_ValueOfPi

      if (RHO.le.0D0) then
         nwxc_EpsX_p = 0.0D0
         return
      end if

      kF = nwxc_FermiK_p(Rho)
      Cs = -3.0D0/(4.0d0*nwxc_ValueOfPi())
      nwxc_EpsX_p = Cs * kF * nwxc_HqBNL_p(gamma/kF)

      return
      end      
#else
      function nwxc_EpsX(Rho,gamma)

#include "nwad.fh"

      implicit none

      type(nwad_dble)::nwxc_EpsX
      type(nwad_dble)::kF,RHO
      double precision gamma,Cs
      type(nwad_dble)::nwxc_HqBNL
      type(nwad_dble)::nwxc_FermiK
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
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      function nwxc_EpsX_d2(Rho,gamma)

#include "nwad.fh"

      implicit none

      type(nwad_dble)::nwxc_EpsX_d2
      type(nwad_dble)::kF,RHO
      double precision gamma,Cs
      type(nwad_dble)::nwxc_HqBNL_d2
      type(nwad_dble)::nwxc_FermiK_d2
      double precision nwxc_ValueOfPi

      if (RHO.le.0D0) then
         nwxc_EpsX_d2 = 0.0D0
         return
      end if

      kF = nwxc_FermiK_d2(Rho)
      Cs = -3.0D0/(4.0d0*nwxc_ValueOfPi())
      nwxc_EpsX_d2 = Cs * kF * nwxc_HqBNL_d2(gamma/kF)

      return
      end      
#else
      function nwxc_EpsX_d3(Rho,gamma)

#include "nwad.fh"

      implicit none

      type(nwad_dble)::nwxc_EpsX_d3
      type(nwad_dble)::kF,RHO
      double precision gamma,Cs
      type(nwad_dble)::nwxc_HqBNL_d3
      type(nwad_dble)::nwxc_FermiK_d3
      double precision nwxc_ValueOfPi

      if (RHO.le.0D0) then
         nwxc_EpsX_d3 = 0.0D0
         return
      end if

      kF = nwxc_FermiK_d3(Rho)
      Cs = -3.0D0/(4.0d0*nwxc_ValueOfPi())
      nwxc_EpsX_d3 = Cs * kF * nwxc_HqBNL_d3(gamma/kF)

      return
      end      
#endif
c
c
c**********************************************************************
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
#include "nwxc_x_bnl.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
#include "nwxc_x_bnl.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
#include "nwxc_x_bnl.F"
#endif
#undef NWAD_PRINT
c $Id$
C> @}
