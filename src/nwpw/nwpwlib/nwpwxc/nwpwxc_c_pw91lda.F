#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_pw91lda.F
C> The Perdew-Wang LDA correlation functional of 1991
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The homogeneous electron gas correlation functional of Perdew and
C> Wang of 1991
C>
C> This has the same form as VWN functional V, only with a different
C> form for the parameterized functionals of \f$r_s\f$. The VWN V code is
C> reused.
C>
C> [1] J.P. Perdew, Y. Wang,
C>     "Accuracte and simple analytic representation of the electron-gas
C>      correlation energy",
C>     Phys. Rev. B <b>45</b>, 13244-13249 (1992), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevB.45.13244">
C>     10.1103/PhysRevB.45.13244</a>.
*
* $Id$
*
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_pw91lda(tol_rho, ipol, nq, wght, rho,
     &                          func, Amat)
#else
      Subroutine nwpwxc_c_pw91lda_d2(tol_rho, ipol, nq, wght, rho,
     &                             func, Amat, Amat2)
#endif
      implicit none
#include "nwpwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     The Functional and the Partial Derivatives of the Correlation 
C     Energy Functional
c
      double precision func(*)     !< [Output] The value of the functional
      double precision Amat(nq,*)  !< [Output] The derivative wrt rho
#ifdef SECOND_DERIV
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
#endif
c
      double precision onethird, fourthirds, twothirds, pi
      Parameter (onethird = 1.D0/3.D0, fourthirds = 4.D0/3.D0)
      Parameter (twothirds = 2.D0/3.D0)
c     Parameter (pi = 3.1415926535898D0)
c
c     Functional Parameters
c
      double precision A(3), alp(3), b(4,3)
      save A, alp, b
c
      double precision e(3), d1e(3), rhoval, rs, d1rs, x, d1x,
     &     h1, d1h1, h2, d1h2,
     &     d1zeta(2), d1ersz(2), d1edrho(2), zeta, fz, d1fz, eps,
     &     dec_rs1, dec_rsz, d1dec_rs1, d1dec_rsz(2)
      double precision devwn_rsz, d1devwn_rsz(2), zeta2, zeta3, zeta4,
     &     d2fz0, beta_rs1, d1beta_rs1, t_vwn, d1t_vwn
#ifdef SECOND_DERIV
      double precision d2beta_rs1, d2t_vwn, d2devwn_rsz(3)
      double precision d2rs, d2x, d2h1, d2h2,
     &     d2e(3), d2zeta(3), d2dec_rs1, d2dec_rsz(3),
     &     d2ersz(3), d2edrho(3), d2fz, rrho2
#endif
      double precision p0, p1, p2, p3
c
      integer i, n, initial
      save initial
      data A   / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
      data alp / 0.21370d0, 0.20548d0, 0.11125d0 /
      data b   / 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0,
     &          14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0,
     &          10.357d0, 3.6231d0, 0.88026d0, 0.49671d0 /
      data initial /1/
c
c     Define miscellaneous parameters.
c
      pi = acos(-1.0d0)
      p0 = (1.0d0/(fourthirds*pi))**onethird
      p1 = 0.5D0/(2.d0**onethird - 1.d0)
      p2 = fourthirds*p1
      p3 = onethird*p2
      d2fz0 = 2.d0*p3
      if (initial.eq.1)then
         initial = 0
c        For convenience store -2A as A and multiply betas by 2A
         do i = 1, 3
            A(i) = -2d0*A(i)
            do n = 1, 4
               b(n,i) = -A(i)*b(n,i)
            enddo
         enddo
c        Finally, change the sign on A for spin stiffness since
c        the negative of that is fitted in the PW'91 paper.  We can't
c        just take the negative of A at the start since A also
c        contributes to the argument of the ln function.
         A(3) = -A(3)
      endif
c
c     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
c
      do 200 n = 1, nq
         if (ipol.eq.1) then 
            rhoval = rho(n,R_T)
         else 
            rhoval = rho(n,R_A)+rho(n,R_B)
         endif
         if (rhoval.lt.tol_rho)goto 200
c
         rs = p0*rhoval**(-onethird)
         d1rs = -onethird*rs/rhoval
         x = sqrt(rs)
         d1x = 0.5d0/x
#ifdef SECOND_DERIV
         d2rs = -fourthirds*d1rs/rhoval
         d2x = -0.5d0*d1x/rs
#endif
c
c        Evaluate the individual correlation energy formulas
c
c        Note that the Monte Carlo form (p = 1) is used for h2.
c
         do i = 1, 3
            h2 = x*(b(1,i) + x*(b(2,i) + x*(b(3,i) + x*b(4,i))))
            d1h2 = b(1,i)
     &           + x*(2d0*b(2,i) + x*(3d0*b(3,i) + 4d0*x*b(4,i)))
#ifdef SECOND_DERIV
            d2h2 = 2d0*b(2,i) + x*(6d0*b(3,i) + 12d0*x*b(4,i))
#endif
c
            h1 = DLOG(1d0+1d0/h2)
            d1h1 = -d1h2/(h2*(h2+1d0))
#ifdef SECOND_DERIV
            d2h1 = d1h1*d1h1*(2d0*h2+1d0) - d2h2/(h2*(h2+1d0))
#endif
c
            e(i) = A(i)*(1d0+alp(i)*rs)*h1
            d1e(i) = A(i)*(2d0*alp(i)*x*h1+(1d0+alp(i)*rs)*d1h1)
#ifdef SECOND_DERIV
            d2e(i) = A(i)*(2d0*alp(i)*h1+4d0*alp(i)*x*d1h1
     &                      +(1d0+alp(i)*rs)*d2h1)
#endif
c
c           Transform derivatives wrt x to derivatives wrt rs
c
#ifdef SECOND_DERIV
c           Do 2nd derivative first so the x first derivative in d1e
c           is not lost
            d2e(i) = d2e(i)*d1x*d1x + d1e(i)*d2x
#endif
            d1e(i) = d1e(i)*d1x
         enddo
c
c        Compute the polarization function and its derivatives
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,R_A) - rho(n,R_B))/rhoval
         endif
         if (zeta.gt.1.d0)then
            zeta = 1.d0
         elseif (zeta.lt.-1.d0)then
            zeta =-1.d0
         endif
         fz = ((1.d0+zeta)**fourthirds + 
     &         (1.d0-zeta)**fourthirds - 2.d0)*p1
         d1fz = ((1.d0+zeta)**onethird - 
     &           (1.d0-zeta)**onethird)*p2
         d1zeta(1) = (1.d0-zeta)/rhoval
         d1zeta(2) =-(1.d0+zeta)/rhoval
#ifdef SECOND_DERIV
         if(dabs(zeta).lt.tol_rho) then
            d2fz = d2fz0
         else
            if (1.0d0+zeta.le.tol_rho) then 
              d2fz = ((1.d0-zeta)**(-twothirds))*p3
            else if (1.0d0-zeta.le.tol_rho) then
              d2fz = ((1.d0+zeta)**(-twothirds))*p3
            else
              d2fz = ((1.d0+zeta)**(-twothirds) +
     &                (1.d0-zeta)**(-twothirds))*p3
            endif
         endif
         rrho2 = 2.d0/(rhoval*rhoval)
c        1 = aa, 2 = ab, 3 = bb
         d2zeta(1) =-rrho2*(1.d0-zeta)
         d2zeta(2) = rrho2*zeta
         d2zeta(3) = rrho2*(1.d0+zeta)
#endif
c
         dec_rs1 = e(2)-e(1)
         d1dec_rs1 = d1e(2)-d1e(1)
#ifdef SECOND_DERIV
         d2dec_rs1 = d2e(2)-d2e(1)
#endif
c
         beta_rs1 = e(2)-e(1)
         d1beta_rs1 = d1e(2)-d1e(1)
         zeta2 = zeta*zeta
         zeta3 = zeta2*zeta
         zeta4 = zeta3*zeta
         t_vwn = d2fz0*beta_rs1-e(3)
         d1t_vwn = d2fz0*d1beta_rs1-d1e(3)
         devwn_rsz = fz/d2fz0*(e(3)+t_vwn*zeta4)
         d1devwn_rsz(1) = fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
         d1devwn_rsz(2) = d1fz/d2fz0*(e(3)+t_vwn*zeta4)
     &        + fz/d2fz0*t_vwn*4.d0*zeta3
#ifdef SECOND_DERIV
         d2beta_rs1 = d2e(2)-d2e(1)
         d2t_vwn = d2fz0*d2beta_rs1-d2e(3)
         d2devwn_rsz(1) = fz/d2fz0*(d2e(3)+d2t_vwn*zeta4)
         d2devwn_rsz(2) = d1fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
     &        + fz/d2fz0*d1t_vwn*4.d0*zeta3
         d2devwn_rsz(3) = d2fz/d2fz0*(e(3)+t_vwn*zeta4)
     &        + d1fz/d2fz0*t_vwn*8.d0*zeta3
     &        + fz/d2fz0*t_vwn*12.d0*zeta2
#endif
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header
c     files.
c
         dec_rsz = devwn_rsz
         d1dec_rsz(1) = d1devwn_rsz(1)
         d1dec_rsz(2) = d1devwn_rsz(2)
#ifdef SECOND_DERIV
         d2dec_rsz(1) = d2devwn_rsz(1)
         d2dec_rsz(2) = d2devwn_rsz(2)
         d2dec_rsz(3) = d2devwn_rsz(3)
#endif
c
c     Finish off the unrestricted case:
c     Assemble the entire functional and its derivatives given the
c     parameterization-dependent part deltaEc(rs,zeta) and its derivatives
c
         eps = e(1) + dec_rsz
         d1ersz(1) = d1e(1) + d1dec_rsz(1)
         d1ersz(2) = d1dec_rsz(2)
         d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
         d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
         func(n) = func(n) + eps*rhoval*wght
         Amat(n,D1_RA) = Amat(n,D1_RA) + (eps+rhoval*d1edrho(1))*wght
         if (ipol.eq.2) then
            Amat(n,D1_RB) = Amat(n,D1_RB) + (eps+rhoval*d1edrho(2))*wght
         endif
#ifdef SECOND_DERIV
c        1 = rsrs, 2 = rsz, 3 = zz
         d2ersz(1) = d2e(1) + d2dec_rsz(1)
         d2ersz(2) = d2dec_rsz(2)
         d2ersz(3) = d2dec_rsz(3)
c        1 = aa, 2 = ab, 3 = bb
         d2edrho(1) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*d1zeta(1)*2.d0
     &              + d2ersz(3)*d1zeta(1)*d1zeta(1)
     &              + d1ersz(1)*d2rs
     &              + d1ersz(2)*d2zeta(1)
         d2edrho(2) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
     &              + d2ersz(3)*d1zeta(1)*d1zeta(2)
     &              + d1ersz(1)*d2rs
     &              + d1ersz(2)*d2zeta(2)
         d2edrho(3) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*d1zeta(2)*2.d0
     &              + d2ersz(3)*d1zeta(2)*d1zeta(2)
     &              + d1ersz(1)*d2rs
     &              + d1ersz(2)*d2zeta(3)
         Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &        + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*wght
         Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &        + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*wght
         if (ipol.eq.2) then
           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &          + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*wght
         endif
#endif
  200 continue
c
      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_c_pw91lda.F"
#endif
C> @}
