#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_camlsd.F
C> The CAM-LSD exchange functional
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the CAM-LSD exchange functional
C>
C> Evaluate the CAM-LSD functional [1,2]. This routine is
C> also used to implement CAM-B3LYP.
C>
C> ### References ###
C>
C> [1] T. Yanai, D.P. Tew, N.C. Handy,
C> "A new hybrid exchange-correlation functional using the Coulomb-attenuating
C> method (CAM-B3LYP)",
C> Chem. Phys. Lett. <b>393</b>, 51-57 (2004), DOI:
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011">
C> 10.1016/j.cplett.2004.06.011</a>.
C>
C> [2] A.D. Becke,
C> "Density-functional exchange-energy approximation with correct
C> asymptotic behavior",
C> Phys. Rev. A <b>38</b>, 3098-3100 (1998), DOI:
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098">
C> 10.1103/PhysRevA.38.3098</a>.
C>
c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_camlsd_p(param, tol_rho, ipol, nq, wght, rho,
     +                           func)
#else
      Subroutine nwxc_x_camlsd(param, tol_rho, ipol, nq, wght, rho,
     +                         func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c     For locations of 2nd derivatives of functionals in array
      Subroutine nwxc_x_camlsd_d2(param, tol_rho, ipol, nq, wght, rho,
     +                            func)
#else
      Subroutine nwxc_x_camlsd_d3(param, tol_rho, ipol, nq, wght, rho,
     +                            func)
#endif
c
C$Id$
c
#include "nwad.fh"
      Implicit none
c
#include "nwxc_param.fh"
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)!< [Input] Parameters of functional
#else
      double precision param(*)!< [Input] Parameters of functional
#endif
#else
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \alpha_{CAM} \f$
                               !< - param(2): \f$ \beta_{CAM} \f$
                               !< - param(3): \f$ \omega_{CAM} \f$
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
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision Atmp,Ctmp,A2tmp,C2tmp,C3tmp
      type(nwad_dble):: Etmp
      double precision A3tmp, C4tmp, C5tmp, C6tmp
      double precision rhom23
      double precision P1, P2, P3, P4
c
c     P1 =       -(3/PI)**(1/3)
c     P2 = -(3/4)*(3/PI)**(1/3)
c     P3 =       -(6/PI)**(1/3)
c     P4 = -(3/4)*(6/PI)**(1/3)
c
      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      double precision  one_third,two_ninth
      type(nwad_dble):: rho13, rho32, rho33
      Parameter (one_third = 1.d0/3.d0)
      Parameter (two_ninth = 2.d0/9.d0)
      integer n
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,R_T).gt.tol_rho)then
             rho13=rho(n,R_T)**one_third
             Etmp = rho(n,R_T)*rho13*P2*wght
c            Atmp = rho13*P1*wght
c            Ctmp = 0.d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = (rho13/rho(n,R_T))*2.0d0*one_third*P1*wght
c            C2tmp = 0.d0
c            C3tmp = 0.d0
#endif
#if defined(THIRD_DERIV)
c            rhom23 = rho13/rho(n,R_T)
c            A3tmp = (rhom23/rho(n,R_T))*(-4.0d0)*two_ninth*P1*wght
c            C4tmp = 0.0d0
c            C5tmp = 0.0d0
c            C6tmp = 0.0d0
#endif
#if defined(THIRD_DERIV)
             call nwxc_x_att_d3(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp)
c
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
#elif defined(SECOND_DERIV)

             call nwxc_x_att_d2(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp)
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#else
#if defined(NWAD_PRINT)
             call nwxc_x_att_p(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp)
#else
             call nwxc_x_att(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp)
#endif
#endif
             func(n) = func(n) + Etmp
c            Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
            endif
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               rho32=rho(n,R_A)**one_third
             else
               rho32=0.0d0
             endif
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               rho33=rho(n,R_B)**one_third
             else
               rho33=0.0d0
             endif
c
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               Etmp = rho32*rho(n,R_A)*P4*wght
             else
               Etmp = 0.0d0
             endif
c            Atmp = P3*rho32*wght
c            Ctmp = 0.d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = 0.d0
c            C2tmp = 0.d0
c            C3tmp = 0.d0
c            if (rho(n,R_A).gt.0.5d0*tol_rho) then
c              A2tmp = one_third*P3*rho32/rho(n,R_A)*wght
c            end if
#endif
#if defined(THIRD_DERIV)
c            A3tmp = 0.0d0
c            C4tmp = 0.0d0
c            C5tmp = 0.0d0
c            C6tmp = 0.0d0
c
c            if (rho(n,R_A).gt.0.5d0*tol_rho) then
c              A3tmp = -two_ninth*P3*rho32/(rho(n,R_A)**2)*wght
c            endif
#endif
#if defined(THIRD_DERIV)
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               call nwxc_x_att_d3(param,tol_rho,rho(n,R_A),ipol,
     &              Etmp)
             endif
c
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
#elif defined(SECOND_DERIV)
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               call nwxc_x_att_d2(param,tol_rho,rho(n,R_A),ipol,
     &              Etmp)
             end if
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#else
#if defined(NWAD_PRINT)
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               call nwxc_x_att_p(param,tol_rho,rho(n,R_A),ipol,
     &              Etmp)
             endif
#else
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               call nwxc_x_att(param,tol_rho,rho(n,R_A),ipol,
     &              Etmp)
             endif
#endif
#endif
             func(n) = func(n) + Etmp
c            Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
c
c            Beta spin channel
c
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               Etmp = rho33*rho(n,R_B)*P4*wght
             else
               Etmp = 0.0d0
             endif
c            Atmp = P3*rho33*wght
c            Ctmp = 0.d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = 0.d0
c            C2tmp = 0.d0
c            C3tmp = 0.d0
c            if (rho(n,R_B).gt.0.5d0*tol_rho) then
c              A2tmp = one_third*P3*rho33/rho(n,R_B)*wght
c            end if
#endif
#if defined(THIRD_DERIV)
c            A3tmp = 0.0d0
c            C4tmp = 0.0d0
c            C5tmp = 0.0d0
c            C6tmp = 0.0d0
c
c            if (rho(n,R_B).gt.0.5d0*tol_rho) then
c              A3tmp = -two_ninth*P3*rho33/(rho(n,R_B)**2)*wght
c            endif
#endif
#if defined(THIRD_DERIV)
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               call nwxc_x_att_d3(param,tol_rho,rho(n,R_B),ipol,
     &              Etmp)
             endif
c
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c
c            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
#elif defined(SECOND_DERIV)
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               call nwxc_x_att_d2(param,tol_rho,rho(n,R_B),ipol,
     &              Etmp)
             end if
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
#else
#if defined(NWAD_PRINT)
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               call nwxc_x_att_p(param,tol_rho,rho(n,R_B),ipol,
     &              Etmp)
             endif
#else
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               call nwxc_x_att(param,tol_rho,rho(n,R_B),ipol,
     &              Etmp)
             endif
#endif
#endif
             func(n) = func(n) + Etmp
c            Amat(n,D1_RB) = Amat(n,D1_RB) + Atmp
c
c            func(n) = func(n) + ( rho32*rho(n,R_A) + 
c    &                             rho33*rho(n,R_B)   )*P4*wght
   20    continue
c 
      endif
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_camlsd.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_camlsd.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_camlsd.F"
#endif
#undef NWAD_PRINT
C> @}
