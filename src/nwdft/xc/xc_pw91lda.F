c     Uniform-gas correlation of Perdew and Wang 1991
c
c     This has the same form as VWN functional V, only with a different
c     form for the parameterized functionals of rs.  The VWN V code is
c     reused.
*
* $Id$
*
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_pw91lda(tol_rho, fac, lfac, nlfac, rho, Amat, nq,
     &                      ipol, Ec, qwght, ldew, func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#include "dft2drv.fh"
      Subroutine xc_pw91lda_d2(tol_rho, fac, lfac, nlfac, rho,
     &     Amat, Amat2, nq, ipol, Ec, qwght, ldew, func)
#else
#include "dft2drv.fh"
#include "dft3drv.fh"
      Subroutine xc_pw91lda_d3(tol_rho, fac, lfac, nlfac, rho,
     &     Amat, Amat2, Amat3, nq, ipol, Ec, qwght, ldew, func)
#endif
      implicit none
c
      integer nq, ipol
      double precision fac, Ec, tol_rho
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial Derivatives of the Correlation Energy Functional
c
      double precision Amat(nq,ipol)
#ifdef SECOND_DERIV
      double precision Amat2(nq,*)
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3)
#endif
c
      double precision onethird, fourthirds, twothirds, pi
      double precision fivethirds, seventhirds, threehalf
      Parameter (onethird = 1.D0/3.D0, fourthirds = 4.D0/3.D0)
      Parameter (twothirds = 2.D0/3.D0)
      Parameter (fivethirds = 5.0d0/3.0d0)
      Parameter (seventhirds = 7.0d0/3.0d0)
      Parameter (threehalf = 3.0d0/2.0d0)
      Parameter (pi = 3.1415926535898D0)
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
      double precision p4
c
#ifdef THIRD_DERIV
      double precision d3beta_rs1, d3t_vwn, d3devwn_rsz(4)
      double precision d3rs, d3x, d3h1, d3h2,
     &     d3e(3), d3zeta(4), d3dec_rs1, d3dec_rsz(4),
     &     d3ersz(4), d3edrho(4), d3fz, rrho3
#endif
c
      integer i, n, initial
      save initial
c Daniel (10-19-12): Parameters are taken from the paper:
c Phys. Rev. B 1992, 45, 13244.
      data A   / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
      data alp / 0.21370d0, 0.20548d0, 0.11125d0 /
      data b   / 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0,
     &          14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0,
     &          10.357d0, 3.6231d0, 0.88026d0, 0.49671d0 /
      data initial /1/
c
c     Define miscellaneous parameters.
c
      p0 = (1.0d0/(fourthirds*pi))**onethird
      p1 = 0.5D0/(2.d0**onethird - 1.d0)
      p2 = fourthirds*p1
      p3 = onethird*p2
c For XC-third derivative
      p4 = -twothirds*p3
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
c        just take the negative of A at the start since A also contributes
c        to the argument of the ln function.
         A(3) = -A(3)
      endif
c
c     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
c
      do 200 n = 1, nq
         if (rho(n,1).lt.tol_rho)goto 200
c
         rhoval = rho(n,1)
         rs = p0*rhoval**(-onethird)
         d1rs = -onethird*rs/rhoval
         x = sqrt(rs)
         d1x = 0.5d0/x
#ifdef SECOND_DERIV
         d2rs = -fourthirds*d1rs/rhoval
         d2x = -0.5d0*d1x/rs
#endif
#ifdef THIRD_DERIV
         d3rs = -seventhirds*d2rs/rhoval
         d3x = -threehalf*d2x/rs
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
#ifdef THIRD_DERIV
            d3h2 = 6.0d0*b(3,i) + 24.0d0*x*b(4,i)
#endif
c
            h1 = DLOG(1d0+1d0/h2)
            d1h1 = -d1h2/(h2*(h2+1d0))
#ifdef SECOND_DERIV
            d2h1 = d1h1*d1h1*(2d0*h2+1d0) - d2h2/(h2*(h2+1d0))
#endif
#ifdef THIRD_DERIV
            d3h1 = 2d0*d2h1*d1h1*(2d0*h2+1d0) 
     1           + 2d0*d1h1*d1h1*d1h2
     2           - d3h2/(h2*(h2+1d0))
     3           - d2h2*d1h1*(2d0*h2+1d0)/(h2*(h2+1d0))
#endif
c
            e(i) = A(i)*(1d0+alp(i)*rs)*h1
            d1e(i) = A(i)*(2d0*alp(i)*x*h1+(1d0+alp(i)*rs)*d1h1)
#ifdef SECOND_DERIV
            d2e(i) = A(i)*(2d0*alp(i)*h1+4d0*alp(i)*x*d1h1
     &                      +(1d0+alp(i)*rs)*d2h1)
#endif
#ifdef THIRD_DERIV
            d3e(i) = A(i)*( 6d0*alp(i)*d1h1 + 6d0*alp(i)*x*d2h1
     1                    + (1d0+alp(i)*rs)*d3h1 )
#endif
c
c           Transform derivatives wrt x to derivatives wrt rs
c
#ifdef THIRD_DERIV
            d3e(i) = d3e(i)*d1x*d1x*d1x + 3.0d0*d2e(i)*d1x*d2x
     &        + d1e(i)*d3x
#endif
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
            zeta = (rho(n,2) - rho(n,3))/rhoval
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
#ifdef THIRD_DERIV
         if (dabs(zeta+1.0d0).le.tol_rho) then
           d3fz = (-(1.0d0-zeta)**(-fivethirds))*p4
         else if (dabs(zeta-1.0d0).le.tol_rho) then
           d3fz = ((1.0d0+zeta)**(-fivethirds))*p4
         else
           d3fz = ((1.0d0+zeta)**(-fivethirds) -
     &             (1.0d0-zeta)**(-fivethirds))*p4
         end if
         rrho3 = 1.0d0/(rhoval*rhoval*rhoval)
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
         d3zeta(1) = 6.0d0*(1.0d0-zeta)*rrho3
         d3zeta(2) = 2.0d0*(1.0d0-3.0d0*zeta)*rrho3
         d3zeta(3) = -2.0d0*(1.0d0+3.0d0*zeta)*rrho3
         d3zeta(4) = -6.0d0*(1.0d0+3.0d0*zeta)*rrho3
#endif
c
         dec_rs1 = e(2)-e(1)
         d1dec_rs1 = d1e(2)-d1e(1)
#ifdef SECOND_DERIV
         d2dec_rs1 = d2e(2)-d2e(1)
#endif
#ifdef THIRD_DERIV
         d3dec_rs1 = d3e(2)-d3e(1)
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
#ifdef THIRD_DERIV
         d3beta_rs1 = d3e(2)-d3e(1)
         d3t_vwn = d2fz0*d3beta_rs1-d3e(3)
c Derivatives: 1 = drsdrsdrs, 2 = drsdrsdzeta, 3 = drsdzetadzeta,
c              4 = dzetadzetadzeta
         d3devwn_rsz(1) = fz/d2fz0*(d3e(3)+d3t_vwn*zeta4)
         d3devwn_rsz(2) = d1fz/d2fz0*(d2e(3)+d2t_vwn*zeta4)
     &        + fz/d2fz0*d2t_vwn*4.0d0*zeta3
         d3devwn_rsz(3) = d2fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
     &        + d1fz/d2fz0*d1t_vwn*8.0d0*zeta3
     &        + fz/d2fz0*d1t_vwn*12.0d0*zeta2
         d3devwn_rsz(4) = d3fz/d2fz0*(e(3)+t_vwn*zeta4)
     &        + d2fz/d2fz0*t_vwn*12.0d0*zeta3
     &        + d1fz/d2fz0*t_vwn*36.0d0*zeta2
     &        + fz/d2fz0*t_vwn*24.0d0*zeta
#endif
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
         dec_rsz = devwn_rsz
         d1dec_rsz(1) = d1devwn_rsz(1)
         d1dec_rsz(2) = d1devwn_rsz(2)
#ifdef SECOND_DERIV
         d2dec_rsz(1) = d2devwn_rsz(1)
         d2dec_rsz(2) = d2devwn_rsz(2)
         d2dec_rsz(3) = d2devwn_rsz(3)
#endif
#ifdef THIRD_DERIV
         d3dec_rsz(1) = d3devwn_rsz(1)
         d3dec_rsz(2) = d3devwn_rsz(2)
         d3dec_rsz(3) = d3devwn_rsz(3)
         d3dec_rsz(4) = d3devwn_rsz(4)
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
         Ec = Ec + eps*qwght(n)*rhoval*fac
         if (ldew) func(n) = func(n) + eps*rhoval*fac
         Amat(n,1) = Amat(n,1) + (eps + rhoval*d1edrho(1))*fac
         if (ipol.eq.2)
     &   Amat(n,2) = Amat(n,2) + (eps + rhoval*d1edrho(2))*fac
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
     &        + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
         Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &        + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*fac
         if (ipol.eq.2)
     &   Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &        + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*fac
#endif
#ifdef THIRD_DERIV
c 1 = rsrsrs, 2 = rsrsz, 3 = rszz, 4 = zzz
         d3ersz(1) = d3e(1) + d3dec_rsz(1)
         d3ersz(2) = d3dec_rsz(2)
         d3ersz(3) = d3dec_rsz(3)
         d3ersz(4) = d3dec_rsz(4)
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
         d3edrho(1) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*d1zeta(1)*d1zeta(1)*3.0d0
     &              + d2ersz(2)*d1rs*d2zeta(1)*3.0d0
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d1zeta(1)*d2rs*3.0d0
     &              + d3ersz(2)*d1zeta(1)*d1rs*d1rs*3.0d0
     &              + d3ersz(4)*d1zeta(1)*d1zeta(1)*d1zeta(1)
     &              + d2ersz(3)*d1zeta(1)*d2zeta(1)*3.0d0
     &              + d1ersz(2)*d3zeta(1)
         d3edrho(2) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*(d1zeta(1)*d1zeta(1)
     &                              + d1zeta(1)*d1zeta(2)*2.0d0)
     &              + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     &                              + d2zeta(1))
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d2rs*(d1zeta(1)*2.0d0
     &                              + d1zeta(2))
     &              + d3ersz(2)*d1rs*d1rs*(d1zeta(2)
     &                                   + d1zeta(1)*2.0d0)
     &              + d3ersz(4)*d1zeta(2)*d1zeta(1)*d1zeta(1)
     &              + d2ersz(3)*(d1zeta(1)*d2zeta(2)*2.0d0
     &                         + d1zeta(2)*d2zeta(1))
     &              + d1ersz(2)*d3zeta(2)
         d3edrho(3) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*(d1zeta(2)*d1zeta(2)
     &                              + d1zeta(2)*d1zeta(1)*2.0d0)
     &              + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     &                              + d2zeta(3))
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d2rs*(d1zeta(2)*2.0d0
     &                              + d1zeta(1))
     &              + d3ersz(2)*d1rs*d1rs*(d1zeta(1)
     &                                   + d1zeta(2)*2.0d0)
     &              + d3ersz(4)*d1zeta(1)*d1zeta(2)*d1zeta(2)
     &              + d2ersz(3)*(d1zeta(2)*d2zeta(2)*2.0d0
     &                         + d1zeta(1)*d2zeta(3))
     &              + d1ersz(2)*d3zeta(3)
         d3edrho(4) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*d1zeta(2)*d1zeta(2)*3.0d0
     &              + d2ersz(2)*d1rs*d2zeta(3)*3.0d0
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d1zeta(2)*d2rs*3.0d0
     &              + d3ersz(2)*d1zeta(2)*d1rs*d1rs*3.0d0
     &              + d3ersz(4)*d1zeta(2)*d1zeta(2)*d1zeta(2)
     &              + d2ersz(3)*d1zeta(2)*d2zeta(3)*3.0d0
     &              + d1ersz(2)*d3zeta(4)
c
         Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     &        + (3.0d0*d2edrho(1) + rhoval*d3edrho(1))*fac
         Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)
     &        + (d2edrho(1) + 2.0d0*d2edrho(2) + rhoval*d3edrho(2))*fac
         Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
     &        + (2.0d0*d2edrho(2) + d2edrho(3) + rhoval*d3edrho(3))*fac
         if (ipol.eq.2)
     &   Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     &        + (3.0d0*d2edrho(3) + rhoval*d3edrho(4))*fac
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
#include "xc_pw91lda.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_pw91lda.F"
#endif
