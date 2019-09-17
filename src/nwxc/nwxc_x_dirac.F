#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_dirac.F
C> The Dirac-Slater exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Implementation of the Dirac-Slater exchange functional
C>
C> This routine implements the Dirac-Slater exchange functional [1-3].
C> If \f$\sigma\f$ refers to the spin channels then this functional at a
C> particular density can be given as
C> \f{eqnarray*}{
C>   K &=& -c \sum_\sigma (\rho_\sigma)^{4/3} \\\\
C>   c &=& \frac{3}{8}4^{2/3}\sqrt[3]{\frac{3}{\pi}}
C> \f}
C> As this functional is simply additive singularities in case the density
C> in the other spin channel \f$\bar{\sigma}\f$ goes to zero can be avoided
C> using
C> \f{eqnarray*}{
C>    K &=& -c (\rho_\sigma)^{4/3}, \; \rho_{\bar{\sigma}} \rightarrow 0
C> \f}
C>
C> [1] F. Bloch, "Bemerkung zur Elektronentheorie des Ferromagnetismus und der elektrischen Leitfähigkeit",
C> Zeitschrift f&uuml;r Physik <b>57</b>, 545-555 (1929), DOI:
C> <a href="https://doi.org/10.1007/BF01340281">10.1007/BF01340281</a>.
C>
C> [2] P.A.M. Dirac, "Note on exchange phenomena in the Thomas atom",
C> Math. Proc. Cambridge Philos. Soc. <b>26</b>, 376-385 (1930), DOI:
C> <a href="https://doi.org/10.1017/S0305004100016108">10.1017/S0305004100016108</a>.
C>
C> [3] J.C. Slater, "A simplification of the Hartree-Fock method",
C> Phys. Rev. <b>81</b>, 385-390 (1951), DOI:
C> <a href="https://doi.org/10.1103/PhysRev.81.385">10.1103/PhysRev.81.385</a>.
c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_dirac_p(tol_rho, ipol, nq, wght, rho, func)
#else
      Subroutine nwxc_x_dirac(tol_rho, ipol, nq, wght, rho, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c     For locations of 2nd derivatives of functionals in array
      Subroutine nwxc_x_dirac_d2(tol_rho, ipol, nq, wght, rho, func)
#else
      Subroutine nwxc_x_dirac_d3(tol_rho, ipol, nq, wght, rho, func)
#endif
c
C$Id$
c
#include "nwad.fh"
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "nwxc_param.fh"
c
      double precision tol_rho  !< [Input] The lower limit on the density
      integer nq                !< [Input] The number of points
      integer ipol              !< [Input] The number of spin channels
      double precision wght     !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,ipol) !< [Input] The density
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
c     Partial Second Order Derivatives of the Exchange Energy Functional
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
      double precision P1, P2, P3, P4, Pi, one_third
c
c     P1 =       -(3/PI)**(1/3)
c     P2 = -(3/4)*(3/PI)**(1/3)
c     P3 =       -(6/PI)**(1/3)
c     P4 = -(3/4)*(6/PI)**(1/3)
c
c     Parameter (P1 = -(3.0d0/Pi)**one_third)
c     Parameter (P2 = (3.0d0/4.0d0)*P1)
c     Parameter (P3 = -(6.0d0/Pi)**one_third)
c     Parameter (P4 = (3.0d0/4.0d0)*P3)
c     Parameter (P1 = -0.9847450218426959D+00)
c     Parameter (P2 = -0.7385587663820219D+00)
c     Parameter (P3 = -0.1240700981798799D+01)
c     Parameter (P4 = -0.9305257363490993D+00)
c     double precision rho13, rho32, rho33, rhom23, rhom53, term
      type(nwad_dble)::rho13
      integer n
c
      Pi = acos(-1.0d0)
      one_third = 1.0d0/3.0d0
      P1 = -(3.0d0/Pi)**one_third
      P2 = (3.0d0/4.0d0)*P1
      P3 = (2.0d0**one_third)*P1
      P4 = (3.0d0/4.0d0)*P3
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef SECOND_DERIV
c         write(LuOut,*)'xc_dirac called: restricted'
#endif
         do n = 1, nq
            if (rho(n,R_T).gt.tol_rho)then
               rho13=rho(n,R_T)**one_third
c              Amat(n,D1_RA) = Amat(n,D1_RA) + wght*P1*rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c              rhom23 = rho13/rho(n,R_T)
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
c    +               wght*2.0d0*one_third*P1*rhom23
#endif
#if defined(THIRD_DERIV)
c              rhom53 = rhom23/rho(n,R_T)
c              Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + 
c    +               -wght*4.0d0*one_third*(2*one_third)*P1*rhom53
#endif
               func(n) = func(n) + wght*P2*rho(n,R_T)*rho13
            endif
         enddo
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef SECOND_DERIV
c         write(LuOut,*)'xc_dirac called: unrestricted'
#endif
         do n = 1,nq
            if (rho(n,R_A).gt.0.5d0*tol_rho)then
               rho13=rho(n,R_A)**one_third
c              Amat(n,D1_RA) = Amat(n,D1_RA) + wght*P3*rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c              rhom23 = rho13/rho(n,R_A)
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
c    +               wght*1.0d0*one_third*P3*rhom23
#endif
#if defined(THIRD_DERIV)
c              rhom53 = rhom23/rho(n,R_A)
c              Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + 
c    +              -wght*1.0d0*one_third*(2*one_third)*P3*rhom53
#endif
               func(n) = func(n) + wght*P4*rho(n,R_A)*rho13
            endif
            if (rho(n,R_B).gt.0.5d0*tol_rho)then
               rho13=rho(n,R_B)**one_third
c              Amat(n,D1_RB) = Amat(n,D1_RB) + wght*P3*rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c              rhom23 = rho13/rho(n,R_B)
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
c    +               wght*1.0d0*one_third*P3*rhom23
#endif
#if defined(THIRD_DERIV)
c              rhom53 = rhom23/rho(n,R_B)
c              Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + 
c    +              -wght*1.0d0*one_third*(2*one_third)*P3*rhom53
#endif
               func(n) = func(n) + wght*P4*rho(n,R_B)*rho13
            endif
         enddo
c 
      endif
c      write(LuOut,*)'from xc_dirac ; p4, exp, fac, Ex: ', 
c     &                           p4, exp, fac, Ex
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_dirac.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_dirac.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_dirac.F"
#endif
#undef NWAD_PRINT
C> @}
