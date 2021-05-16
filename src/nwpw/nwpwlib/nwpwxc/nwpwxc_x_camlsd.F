#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_camlsd.F
C> The CAM-LSD exchange functional
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_camlsd(param, tol_rho, ipol, nq, wght, rho,
     +                         func, Amat)
#else
c     For locations of 2nd derivatives of functionals in array
      Subroutine nwpwxc_x_camlsd_d2(param, tol_rho, ipol, nq, wght, rho,
     +                           func, Amat, Amat2)
#endif
c
C$Id$
c
      Implicit none
c
#include "nwpwxc_param.fh"
c
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \alpha_{CAM} \f$
                               !< - param(2): \f$ \beta_{CAM} \f$
                               !< - param(3): \f$ \omega_{CAM} \f$
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
#ifdef SECOND_DERIV
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,*) !< [Output] 2nd order partial derivatives
#endif
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp
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
      double precision rho13, rho32, rho33, one_third
      Parameter (one_third = 1.d0/3.d0)
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
             Atmp = rho13*P1*wght
             Ctmp = 0.d0
             call nwpwxc_x_att(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp,Atmp,Ctmp)
#ifdef SECOND_DERIV
             A2tmp = (rho13/rho(n,R_T))*2.0d0*one_third*P1*wght
             C2tmp = 0.d0
             C3tmp = 0.d0
             call nwpwxc_x_att_d2(param,tol_rho,rho(n,R_T),ipol,
     &            Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#endif
             func(n) = func(n) + Etmp
             Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
            endif
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
             rho32=max(rho(n,R_A),0.0d0)**one_third
             rho33=max(rho(n,R_B),0.0d0)**one_third
c
             Etmp = rho32*rho(n,R_A)*P4*wght
             Atmp = P3*rho32*wght
             Ctmp = 0.d0
             call nwpwxc_x_att(param,tol_rho,rho(n,R_A),ipol,
     &            Etmp,Atmp,Ctmp)
             func(n) = func(n) + Etmp
             Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
#ifdef SECOND_DERIV
             A2tmp = 0.d0
             C2tmp = 0.d0
             C3tmp = 0.d0
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               A2tmp = one_third*P3*rho32/rho(n,R_A)*wght
               call nwpwxc_x_att_d2(param,tol_rho,rho(n,R_A),ipol,
     &              Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
             end if
             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
#endif
c
             Etmp = rho33*rho(n,R_B)*P4*wght
             Atmp = P3*rho33*wght
             Ctmp = 0.d0
             call nwpwxc_x_att(param,tol_rho,rho(n,R_B),ipol,
     &            Etmp,Atmp,Ctmp)
             func(n) = func(n) + Etmp
             Amat(n,D1_RB) = Amat(n,D1_RB) + Atmp
#ifdef SECOND_DERIV
             A2tmp = 0.d0
             C2tmp = 0.d0
             C3tmp = 0.d0
             if (rho(n,R_B).gt.0.5d0*tol_rho) then
               A2tmp = one_third*P3*rho33/rho(n,R_B)*wght
               call nwpwxc_x_att_d2(param,tol_rho,rho(n,R_B),ipol,
     &              Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
             end if
             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
#endif
c
c            func(n) = func(n) + ( rho32*rho(n,R_A) + 
c    &                             rho33*rho(n,R_B)   )*P4*wght
   20    continue
c 
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_x_camlsd.F"
#endif
C> @}
