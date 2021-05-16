#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_perdew81.F
C> The Perdew Zunger correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the Perdew Zunger correlation functional
C>
C> Evaluates the Perdew Zunger local correlation density functional [1].
C>
C> ### References ###
C>
C> [1] J.P. Perdew, A. Zunger,
C>     "Self-interaction correction to density-functional approximations
C>     for many-electron systems", Phys. Rev. B <b>23</b>, 
C>     5048–5079 (1981), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevB.23.5048">
C>     10.1103/PhysRevB.23.5048</a>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_perdew81(tol_rho, ipol, nq, wght, rho, 
     &                           func, Amat)
#else
      Subroutine nwpwxc_c_perdew81_d2(tol_rho, ipol, nq, wght, rho,
     &                              func, Amat, Amat2)
#endif
c
c     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
c
      implicit none
#include "nwpwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin-channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     The functional values
c
      double precision func(*) !< [Output] The value of the functional
c
c     Partial Derivatives of the Correlation Energy Functional
c
      double precision Amat(nq,*)  !< [Output] The derivative wrt rho
#ifdef SECOND_DERIV
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
#endif
c
      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
     &                 pi, ONE3, FOUR3, TWO3
      save A, B, C, D, G, B1, B2
c     parameter (pi = 3.1415926535897932385d0)
      Parameter (ONE3 = 1.d0/3.d0, FOUR3 = 4.D0/3.D0)
      Parameter (TWO3 = 2.d0/3.d0)
      integer n, i
      double precision rhoval, rs, alnrs, d1rs, e(2), d1e(2), rden(2),
     &                 d1den(2), d1zeta(2), d1ersz(2), d1edrho(2), eps,
     &                 sqrtrs, fz, d1fz, zeta
#ifdef SECOND_DERIV
      double precision d2rs, d2e(2), d2den(2), d2zeta(3), d2ersz(3),
     &                 d2edrho(3), d2fzeta, d2fz, rrho2
#endif
      double precision x, fzeta, d1fzeta, rsfact
      fzeta(x) = ((1.d0+x)**FOUR3 + 
     &            (1.d0-x)**FOUR3 - 2.d0) / (2.d0**FOUR3-2.d0)
      d1fzeta(x) = FOUR3*((1.d0+x)**ONE3 - 
     &                    (1.d0-x)**ONE3) / (2.d0**FOUR3-2.d0)
#ifdef SECOND_DERIV
      d2fzeta(x) = ONE3*FOUR3*((1.d0+x)**(-TWO3) + 
     &                    (1.d0-x)**(-TWO3)) / (2.d0**FOUR3-2.d0)
#endif
      data A / 0.0311d0, 0.01555d0 /
      data B / -0.048d0, -0.0269d0 /
      data C / 0.0020d0, 0.0007d0 /
      data D / -0.0116d0, -0.0048d0 /
      data G / -.1423d0, -.0843d0 /
      data B1 / 1.0529d0, 1.3981d0 /
      data B2 / 0.3334d0, 0.2611d0 /
c
      pi = acos(-1.0d0)
      rsfact = (0.75d0/pi)**ONE3
c
c     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
c
      do n = 1, nq
         if (ipol.eq.1) then
            rhoval = rho(n,R_T)
         else
            rhoval = rho(n,R_A)+rho(n,R_B)
         endif
         if (rhoval.gt.tol_rho)then
            if (ipol.eq.1) then
               zeta = 0.0d0
               d1zeta(1) = 1.d0/rhoval
               d1zeta(2) =-1.d0/rhoval
               fz = 0d0
               d1fz = 0d0
            else
               zeta = (rho(n,R_A)-rho(n,R_B))/rhoval
               if(zeta.lt.-1d0) zeta=-1d0
               if(zeta.gt.1d0) zeta=1d0
               fz = fzeta(zeta)
               d1fz = d1fzeta(zeta)
               d1zeta(1) = (1.d0-zeta)/rhoval
               d1zeta(2) =-(1.d0+zeta)/rhoval
            endif
            rs = rsfact/abs(rhoval)**ONE3
            d1rs = -ONE3*rs/rhoval
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs/rhoval
            if ((1.d0-zeta).lt.tol_rho) then
              d2fz = (1.d0+zeta)**(-TWO3)
            else if ((1.d0+zeta).lt.tol_rho) then
              d2fz = (1.d0-zeta)**(-TWO3)
            else
              d2fz = (1.d0+zeta)**(-TWO3) + (1.d0-zeta)**(-TWO3)
            end if
            d2fz = d2fz*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
            rrho2 = 2.d0/(rhoval*rhoval)
c           1 = aa, 2 = ab, 3 = bb
            d2zeta(1) =-rrho2*(1.d0-zeta)
            d2zeta(2) = rrho2*zeta
            d2zeta(3) = rrho2*(1.d0+zeta)
#endif
            if (rs.lt.1.d0)then
               alnrs = log(rs)
               do i = 1, 2
                  e(i) = A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs
                  d1e(i) = A(i)/rs+C(i)*(alnrs+1d0)+D(i)
#ifdef SECOND_DERIV
                  d2e(i) = (C(i)-A(i)/rs)/rs
#endif
               enddo
            else
               sqrtrs = sqrt(rs)
               do i = 1, 2
                  rden(i) = 1.d0/(1.d0+B1(i)*sqrtrs+B2(i)*rs)
                  d1den(i) = B1(i)/(2.d0*sqrtrs)+B2(i)
                  e(i) = G(i)*rden(i)
                  d1e(i) = -G(i)*d1den(i)*rden(i)**2
#ifdef SECOND_DERIV
                  d2den(i) = -B1(i)/(4.d0*rs*sqrtrs)
                  d2e(i) = G(i)*rden(i)**2
     &                 *(2.d0*d1den(i)**2*rden(i)-d2den(i))
#endif
               enddo
            endif
            eps = e(1) + fz*(e(2)-e(1))
            d1ersz(1) = d1e(1) + fz*(d1e(2)-d1e(1))
            d1ersz(2) = d1fz*(e(2)-e(1))
            d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
            d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
            func(n) = func(n) + eps*rhoval*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + (eps+rhoval*d1edrho(1))*wght
            if (ipol.eq.2)
     &      Amat(n,D1_RB) = Amat(n,D1_RB) + (eps+rhoval*d1edrho(2))*wght
#ifdef SECOND_DERIV
c           1 = rsrs, 2 = rsz, 3 = zz
            d2ersz(1) = d2e(1) + fz*(d2e(2)-d2e(1))
            d2ersz(2) = d1fz*(d1e(2)-d1e(1))
            d2ersz(3) = d2fz*(e(2)-e(1))
c           1 = aa, 2 = ab, 3 = bb
            d2edrho(1) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*d1zeta(1)*2.d0
     &                 + d2ersz(3)*d1zeta(1)*d1zeta(1)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(1)
            d2edrho(2) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
     &                 + d2ersz(3)*d1zeta(1)*d1zeta(2)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(2)
            d2edrho(3) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*d1zeta(2)*2.d0
     &                 + d2ersz(3)*d1zeta(2)*d1zeta(2)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(3)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*wght
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &           + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*wght
            if (ipol.eq.2)
     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*wght
#endif
         endif
      enddo
      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_c_perdew81.F"
#endif
C>
C> @}
