#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_perdew81.F
C> The Perdew Zunger correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_perdew81_p(tol_rho, ipol, nq, wght, rho, 
     &                             func)
#else
      Subroutine nwxc_c_perdew81(tol_rho, ipol, nq, wght, rho, 
     &                           func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_perdew81_d2(tol_rho, ipol, nq, wght, rho,
     &                              func)
#else
      Subroutine nwxc_c_perdew81_d3(tol_rho, ipol, nq, wght, rho,
     &                              func)
#endif
c
c     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
c
#include "nwad.fh"
c
      implicit none
#include "nwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin-channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     The functional values
c
      type(nwad_dble)::func(*) !< [Output] The value of the functional
c
c     Partial Derivatives of the Correlation Energy Functional
c
c     double precision Amat(nq,*)  !< [Output] The derivative wrt rho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
#endif
#if defined(THIRD_DERIV)
c     double precision Amat3(nq,*) !< [Output] The 3rd derivative wrt rho
#endif
c
      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
     &                 pi, ONE3, FOUR3, TWO3
      double precision FIVE3, SEVEN3
      save A, B, C, D, G, B1, B2
c     parameter (pi = 3.1415926535897932385d0)
      Parameter (ONE3 = 1.d0/3.d0, FOUR3 = 4.D0/3.D0)
      Parameter (TWO3 = 2.d0/3.d0)
      Parameter (FIVE3 = 5.0d0/3.0d0, SEVEN3 = 7.0d0/3.0d0)
      integer n, i
      type(nwad_dble)::rhoval, rs, alnrs, e(2), eps, sqrtrs, fz, zeta
      type(nwad_dble)::rden(2)
      double precision d1e(2), 
     &                 d1den(2), d1zeta(2), d1ersz(2), d1edrho(2), 
     &                 d1fz
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision d2rs, d2e(2), d2den(2), d2zeta(3), d2ersz(3),
c    &                 d2edrho(3), d2fzeta, d2fz, rrho2
#endif
#if defined(THIRD_DERIV)
c     double precision d3rs, d3fz, rrho3, d3zeta(4), d3den(2), d3e(2),
c    1                 d3ersz(4), d3edrho(4)
#endif
      double precision x, fzeta, d1fzeta, rsfact
c     fzeta(x) = ((1.d0+x)**FOUR3 + 
c    &            (1.d0-x)**FOUR3 - 2.d0) / (2.d0**FOUR3-2.d0)
c     d1fzeta(x) = FOUR3*((1.d0+x)**ONE3 - 
c    &                    (1.d0-x)**ONE3) / (2.d0**FOUR3-2.d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     d2fzeta(x) = ONE3*FOUR3*((1.d0+x)**(-TWO3) + 
c    &                    (1.d0-x)**(-TWO3)) / (2.d0**FOUR3-2.d0)
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
            rhoval = 0.0d0
            if (rho(n,R_A).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_A)
            endif
            if (rho(n,R_B).gt.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_B)
            endif
         endif
         if (rhoval.gt.tol_rho)then
            if (ipol.eq.1) then
               zeta = 0.0d0
c              d1zeta(1) = 1.d0/rhoval
c              d1zeta(2) =-1.d0/rhoval
               fz = 0d0
               d1fz = 0d0
            else
               zeta = (rho(n,R_A)-rho(n,R_B))/rhoval
               if(zeta.le.-1d0) zeta=-1d0
               if(zeta.ge.1d0) zeta=1d0
               fz = -2.0d0
               if (1.0d0+zeta.gt.1.0d-10) then
                 fz = fz + (1.d0+zeta)**FOUR3  
               endif
               if (1.0d0-zeta.gt.1.0d-10) then
                 fz = fz + (1.d0-zeta)**FOUR3  
               endif
               fz = fz / (2.d0**FOUR3-2.d0)
c              fz = fzeta(zeta)
c              d1fz = d1fzeta(zeta)
c              d1zeta(1) = (1.d0-zeta)/rhoval
c              d1zeta(2) =-(1.d0+zeta)/rhoval
            endif
            rs = rsfact/rhoval**ONE3
c           d1rs = -ONE3*rs/rhoval
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2rs = -FOUR3*d1rs/rhoval
c           if ((1.d0-zeta).lt.tol_rho) then
c             d2fz = (1.d0+zeta)**(-TWO3)
c           else if ((1.d0+zeta).lt.tol_rho) then
c             d2fz = (1.d0-zeta)**(-TWO3)
c           else
c             d2fz = (1.d0+zeta)**(-TWO3) + (1.d0-zeta)**(-TWO3)
c           end if
c           d2fz = d2fz*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
c           rrho2 = 2.d0/(rhoval*rhoval)
c           1 = aa, 2 = ab, 3 = bb
c           d2zeta(1) =-rrho2*(1.d0-zeta)
c           d2zeta(2) = rrho2*zeta
c           d2zeta(3) = rrho2*(1.d0+zeta)
#endif
#if defined(THIRD_DERIV)
c           d3rs = -SEVEN3*d2rs/rhoval 
c           if ((1.d0-zeta).lt.tol_rho) then
c             d3fz = (1.d0+zeta)**(-FIVE3)
c           else if ((1.d0+zeta).lt.tol_rho) then
c             d3fz = (1.d0-zeta)**(-FIVE3)
c           else
c             d3fz = (1.d0+zeta)**(-FIVE3) + (1.d0-zeta)**(-FIVE3)
c           end if
c           d3fz = -d3fz*TWO3*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
c           rrho3 = 2.0d0/(rhoval*rhoval*rhoval)
c
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
c           d3zeta(1) = 3.0d0*rrho3*(1.0d0 - zeta)
c           d3zeta(2) = rrho3*(1.0d0 - 3.0d0*zeta)
c           d3zeta(3) = -rrho3*(1.0d0 + 3.0d0*zeta)
c           d3zeta(4) = -3.0d0*rrho3*(1.0d0 + zeta)
#endif
            alnrs = log(rs)
            sqrtrs = sqrt(rs)
            do i = 1, 2
              rden(i) = 1.d0/(1.d0+B1(i)*sqrtrs+B2(i)*rs)
              e(i) = heaviside(1.0d0-rs)*(
     +                   A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs)
     +             + heaviside(rs-1.0d0)*G(i)*rden(i)
            enddo
c           if (rs.lt.1.d0)then
c              alnrs = log(rs)
c              do i = 1, 2
c                 e(i) = A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs
c                 d1e(i) = A(i)/rs+C(i)*(alnrs+1d0)+D(i)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c                 d2e(i) = (C(i)-A(i)/rs)/rs
#endif
#if defined(THIRD_DERIV)
c                 d3e(i) = 2.0d0*A(i)/(rs*rs*rs)
c    1                   - C(i)/(rs*rs) 
#endif
c              enddo
c           else
c              sqrtrs = sqrt(rs)
c              do i = 1, 2
c                 rden(i) = 1.d0/(1.d0+B1(i)*sqrtrs+B2(i)*rs)
c                 d1den(i) = B1(i)/(2.d0*sqrtrs)+B2(i)
c                 e(i) = G(i)*rden(i)
c                 d1e(i) = -G(i)*d1den(i)*rden(i)**2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c                 d2den(i) = -B1(i)/(4.d0*rs*sqrtrs)
c                 d2e(i) = G(i)*rden(i)**2
c    &                 *(2.d0*d1den(i)**2*rden(i)-d2den(i))
#endif
#if defined(THIRD_DERIV)
c                 d3den(i) = 3.0d0*B1(i)/(8.0d0*rs*rs*sqrtrs)
c                 d3e(i) = G(i)*rden(i)*rden(i)*
c    1                   ( 6.0d0*( d1den(i)*d2den(i)*rden(i)
c    2                           - d1den(i)*d1den(i)*d1den(i)*
c    3                             rden(i)*rden(i) ) 
c    4                   - d3den(i) )
#endif
c              enddo
c           endif
            eps = e(1) + fz*(e(2)-e(1))
c           d1ersz(1) = d1e(1) + fz*(d1e(2)-d1e(1))
c           d1ersz(2) = d1fz*(e(2)-e(1))
c           d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
c           d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
            func(n) = func(n) + eps*rhoval*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) + (eps+rhoval*d1edrho(1))*wght
c           if (ipol.eq.2)
c    &      Amat(n,D1_RB) = Amat(n,D1_RB) + (eps+rhoval*d1edrho(2))*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           1 = rsrs, 2 = rsz, 3 = zz
c           d2ersz(1) = d2e(1) + fz*(d2e(2)-d2e(1))
c           d2ersz(2) = d1fz*(d1e(2)-d1e(1))
c           d2ersz(3) = d2fz*(e(2)-e(1))
c           1 = aa, 2 = ab, 3 = bb
c           d2edrho(1) = d2ersz(1)*d1rs*d1rs
c    &                 + d2ersz(2)*d1rs*d1zeta(1)*2.d0
c    &                 + d2ersz(3)*d1zeta(1)*d1zeta(1)
c    &                 + d1ersz(1)*d2rs
c    &                 + d1ersz(2)*d2zeta(1)
c           d2edrho(2) = d2ersz(1)*d1rs*d1rs
c    &                 + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
c    &                 + d2ersz(3)*d1zeta(1)*d1zeta(2)
c    &                 + d1ersz(1)*d2rs
c    &                 + d1ersz(2)*d2zeta(2)
c           d2edrho(3) = d2ersz(1)*d1rs*d1rs
c    &                 + d2ersz(2)*d1rs*d1zeta(2)*2.d0
c    &                 + d2ersz(3)*d1zeta(2)*d1zeta(2)
c    &                 + d1ersz(1)*d2rs
c    &                 + d1ersz(2)*d2zeta(3)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*wght
c           Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &           + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*wght
c           if (ipol.eq.2)
c    &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &           + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*wght
#endif
#if defined(THIRD_DERIV)
c 1 = rsrsrs, 2 = rsrsz, 3 = rszz, 4 = zzz
c           d3ersz(1) = d3e(1) + fz*(d3e(2)-d3e(1))
c           d3ersz(2) = d1fz*(d2e(2)-d2e(1))
c           d3ersz(3) = d2fz*(d1e(2)-d1e(1))
c           d3ersz(4) = d3fz*(e(2)-e(1))
c
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
c           d3edrho(1) = d3ersz(1)*d1rs*d1rs*d1rs
c    1                 + d2ersz(1)*d1rs*d2rs*3.0d0
c    2                 + d3ersz(3)*d1rs*d1zeta(1)*d1zeta(1)*3.0d0
c    3                 + d2ersz(2)*d1rs*d2zeta(1)*3.0d0
c    4                 + d1ersz(1)*d3rs
c    5                 + d2ersz(2)*d1zeta(1)*d2rs*3.0d0
c    6                 + d3ersz(2)*d1zeta(1)*d1rs*d1rs*3.0d0
c    7                 + d3ersz(4)*d1zeta(1)*d1zeta(1)*d1zeta(1)
c    8                 + d2ersz(3)*d1zeta(1)*d2zeta(1)*3.0d0
c    9                 + d1ersz(2)*d3zeta(1)
c           d3edrho(2) = d3ersz(1)*d1rs*d1rs*d1rs
c    1                 + d2ersz(1)*d1rs*d2rs*3.0d0
c    2                 + d3ersz(3)*d1rs*(d1zeta(1)*d1zeta(1)
c    3                                 + d1zeta(1)*d1zeta(2)*2.0d0)
c    4                 + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
c    5                                 + d2zeta(1))
c    6                 + d1ersz(1)*d3rs
c    7                 + d2ersz(2)*d2rs*(d1zeta(1)*2.0d0
c    8                                 + d1zeta(2))
c    9                 + d3ersz(2)*d1rs*d1rs*(d1zeta(2)
c    A                                      + d1zeta(1)*2.0d0)
c    B                 + d3ersz(4)*d1zeta(2)*d1zeta(1)*d1zeta(1)
c    C                 + d2ersz(3)*(d1zeta(1)*d2zeta(2)*2.0d0
c    D                            + d1zeta(2)*d2zeta(1))
c    E                 + d1ersz(2)*d3zeta(2)
c           d3edrho(3) = d3ersz(1)*d1rs*d1rs*d1rs
c    1                 + d2ersz(1)*d1rs*d2rs*3.0d0
c    2                 + d3ersz(3)*d1rs*(d1zeta(2)*d1zeta(2)
c    3                                 + d1zeta(2)*d1zeta(1)*2.0d0)
c    4                 + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
c    5                                 + d2zeta(3))
c    6                 + d1ersz(1)*d3rs
c    7                 + d2ersz(2)*d2rs*(d1zeta(2)*2.0d0
c    8                                 + d1zeta(1))
c    9                 + d3ersz(2)*d1rs*d1rs*(d1zeta(1)
c    A                                      + d1zeta(2)*2.0d0)
c    B                 + d3ersz(4)*d1zeta(1)*d1zeta(2)*d1zeta(2)
c    C                 + d2ersz(3)*(d1zeta(2)*d2zeta(2)*2.0d0
c    D                            + d1zeta(1)*d2zeta(3))
c    E                 + d1ersz(2)*d3zeta(3)
c           d3edrho(4) = d3ersz(1)*d1rs*d1rs*d1rs
c    1                 + d2ersz(1)*d1rs*d2rs*3.0d0
c    2                 + d3ersz(3)*d1rs*d1zeta(2)*d1zeta(2)*3.0d0
c    3                 + d2ersz(2)*d1rs*d2zeta(3)*3.0d0
c    4                 + d1ersz(1)*d3rs
c    5                 + d2ersz(2)*d1zeta(2)*d2rs*3.0d0
c    6                 + d3ersz(2)*d1zeta(2)*d1rs*d1rs*3.0d0
c    7                 + d3ersz(4)*d1zeta(2)*d1zeta(2)*d1zeta(2)
c    8                 + d2ersz(3)*d1zeta(2)*d2zeta(3)*3.0d0
c    9                 + d1ersz(2)*d3zeta(4)
c
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    1           + ( 3.0d0*d2edrho(1) + rhoval*d3edrho(1) )*wght
c           Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)
c    1           + ( d2edrho(1) + 2.0d0*d2edrho(2) 
c    2             + rhoval*d3edrho(2) )*wght
c           Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
c    1           + ( 2.0d0*d2edrho(2) + d2edrho(3) 
c    2             + rhoval*d3edrho(3) )*wght
c           if (ipol.eq.2)
c    1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
c    2           + ( 3.0d0*d2edrho(3) + rhoval*d3edrho(4) )*wght
#endif
         endif
      enddo
      return
      end
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_perdew81.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_perdew81.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_perdew81.F"
#endif
#undef NWAD_PRINT
C>
C> @}
c $Id$
