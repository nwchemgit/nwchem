#ifndef SECOND_DERIV
      Subroutine xc_perdew86(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Cmat, nq, ipol, Ec, qwght,
     ,     ldew, ffunc)
c
c$Id: xc_perdew86.F,v 1.18 2000-05-05 01:57:27 bjohnson Exp $
c
      implicit none
c
#include "dft2drv.fh"
c
      double precision tol_rho, fac ! [input]
      integer nq, ipol              ! [input]
      double precision Ec           ! [input/output]
      logical lfac, nlfac, ldew
      double precision ffunc(*)  ! value of the functional [output]
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONETHIRD, THREE, FOURTHIRDS, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONETHIRD = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOURTHIRDS = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision drsdrho, gammaval, dsqgamma, dphidrho
      double precision anum, aden, CC, CCPRIME, expfac, phi, func
      double precision AAadd, g, zeta, add, adp, dzetada, dzetadb
      double precision gaa, gab, gbb
      double precision dzeta, dprime
      double precision x
      dzeta(x) = 2.d0**ONETHIRD*dsqrt(((ONE + x)*.5d0)**FIVE3 + 
     &                                ((ONE - x)*.5d0)**FIVE3   )
      dprime(x) = (4.d0**(-ONETHIRD)*FIVE6)*(((ONE + x)*.5d0)**TWO3-
     -              ((ONE - x)*.5d0)**TWO3 ) /
     &              dsqrt(((ONE + x)*.5d0)**FIVE3 +
     &                    ((ONE - x)*.5d0)**FIVE3)
c
      rsfact = (0.75d0/pi)**ONETHIRD
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).lt.tol_rho) goto 10
            rho13  = rhoval**ONETHIRD
            rho43  = rhoval*rho13
            rho76  = rhoval**SEV6 
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONETHIRD*rs/rhoval
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
            dsqgamma = dsqrt(gammaval)
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
            CC = zzz+(anum*aden)
            CCPRIME = aden*((alpha+beta*(rs+rs)) -
     &                (anum*aden)*(pgamma+(rs+rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
            if (phi.lt.EXPTOL)expfac = exp(-phi)
c         
c           functional
c         
            func = expfac*CC*gammaval/rho43*fac
            dphidrho =  -phi * ( SEV6/rhoval  +   (CCPRIME/CC)) 
            AAadd = ((CCPRIME/CC) - (dphidrho + 
     &              FOURTHIRDS/rhoval))*func
            Amat(n,1) = Amat(n,1) + AAadd*fac
            if (dsqgamma.gt.0.d0)then
               g = func/gammaval*(2.d0-phi)*fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + g
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).lt.tol_rho) goto 20
            arho=1.d0/rhoval
            rho13  = rhoval**ONETHIRD
            rho43  = rhoval*rho13
            rho76  = rhoval**SEV6 
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONETHIRD*rs*arho
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1) +
     &                 delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,2)*delrho(n,3,2) +
     &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,2)   )
            dsqgamma = sqrt(gammaval)
            zeta = (rho(n,2) - rho(n,3))*arho
            add = 1.d0/dzeta(zeta)
            adp = dprime(zeta)
            dzetada =  (1.d0 - zeta)*arho*adp
            dzetadb = -(1.d0 + zeta)*arho*adp
c         
c           C(n)
c         
            anum = fff + alpha*rs + beta*rs2
            aden = 1.d0/(1.d0 + pgamma*rs + delta*rs2 + beta10*rs3)
            CC = zzz + (anum*aden)
            CCPRIME = aden*((alpha + beta*(rs + rs)) -
     &                (anum*aden)*(pgamma+(rs + rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
c         
c           functional
c         
            func = expfac*CC*gammaval/rho43*add
            dphidrho =  -phi*(SEV6/rhoval + (CCPRIME/CC)) 
            AAadd = ((CCPRIME/CC) -
     &               (dphidrho + FOURTHIRDS/rhoval))*func
            Amat(n,1) = Amat(n,1) + (AAadd - (func*add)*dzetada)*fac
            Amat(n,2) = Amat(n,2) + (AAadd - (func*add)*dzetadb)*fac
            if (dsqgamma.gt.TOLL)then
               gaa = func/gammaval*(2.d0-phi)*fac
               gbb = gaa
               gab = gaa
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + gaa*0.5D0
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + gab
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + gbb*0.5D0
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   20    continue
      endif
      return
      end
#endif

#ifndef SECOND_DERIV
      Subroutine xc_p81(tol_rho, fac, lfac, nlfac, rho, Amat, nq, ipol, 
     &                  Ec, qwght, ldew, func)
#else
#include "dft2drv.fh"
      Subroutine xc_p81_d2(tol_rho, fac, lfac, nlfac, rho, Amat, Amat2,
     &                     nq, ipol, Ec, qwght, ldew, func)
#endif
c
c     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
c
      implicit none
c
      integer nq, ipol
      logical lfac, nlfac, ldew
      double precision func(*)  ! value of the functional [output]
      double precision Ec, fac
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol)
#ifdef SECOND_DERIV
      double precision Amat2(nq,*)
#endif
      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
     &                 pi, tol_rho, ONE3, FOUR3, TWO3
      data A / 0.0311d0, 0.01555d0 /
      data B / -0.048d0, -0.0269d0 /
      data C / 0.0020d0, 0.0007d0 /
      data D / -0.0116d0, -0.0048d0 /
      data G / -.1423d0, -.0843d0 /
      data B1 / 1.0529d0, 1.3981d0 /
      data B2 / 0.3334d0, 0.2611d0 /
      save A, B, C, D, G, B1, B2
      parameter (pi = 3.1415926535897932385d0)
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
c
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do n = 1, nq
            if (abs(rho(n,1)).gt.tol_rho)then
               rhoval = rho(n,1)
               rs = rsfact/rhoval**ONE3
               d1rs = -ONE3*rs/rhoval
#ifdef SECOND_DERIV
               d2rs = -FOUR3*d1rs/rhoval
#endif
               if (rs.lt.1.d0)then
                  alnrs = log(rs)
                  e(1) = A(1)*alnrs+B(1)+C(1)*rs*alnrs+D(1)*rs
                  d1e(1) = A(1)/rs+C(1)*(alnrs+1d0)+D(1)
#ifdef SECOND_DERIV
                  d2e(1) = (C(1)-A(1)/rs)/rs
#endif
               else
                  sqrtrs = sqrt(rs)
                  rden(1) = 1.d0/(1.d0+B1(1)*sqrtrs+B2(1)*rs)
                  d1den(1) = B1(1)/(2.d0*sqrtrs)+B2(1)
                  e(1) = G(1)*rden(1)
                  d1e(1) = -G(1)*d1den(1)*rden(1)**2
#ifdef SECOND_DERIV
                  d2den(1) = -B1(1)/(4.d0*rs*sqrtrs)
                  d2e(1) = G(1)*rden(1)**2
     &                 *(2.d0*d1den(1)**2*rden(1)-d2den(1))
#endif
               endif
               d1edrho(1) = d1e(1)*d1rs
               Ec = Ec + e(1)*qwght(n)*rhoval*fac
               if (ldew) func(n) = func(n) + e(1)*rhoval*fac
               Amat(n,1) = Amat(n,1) + (e(1)+rhoval*d1edrho(1))*fac
#ifdef SECOND_DERIV
               d2edrho(1) = d2e(1)*d1rs*d1rs
     &                    + d1e(1)*d2rs
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
c     ab contribution is same as aa contribution for restricted
c     since in that case this functional depends on total rho only -
c     it's not that simple in general
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
#endif
            endif
         enddo
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do n = 1, nq
            if (abs(rho(n,1)).gt.tol_rho)then
               rhoval = rho(n,1)
               zeta = (rho(n,2)-rho(n,3))/rhoval
               rs = rsfact/rhoval**ONE3
               d1rs = -ONE3*rs/rhoval
               fz = fzeta(zeta)
               d1fz = d1fzeta(zeta)
               d1zeta(1) = (1.d0-zeta)/rhoval
               d1zeta(2) =-(1.d0+zeta)/rhoval
#ifdef SECOND_DERIV
               d2rs = -FOUR3*d1rs/rhoval
               d2fz = d2fzeta(zeta)
               rrho2 = 2.d0/(rhoval*rhoval)
c              1 = aa, 2 = ab, 3 = bb
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
     &                    *(2.d0*d1den(i)**2*rden(i)-d2den(i))
#endif
                  enddo
               endif
               eps = e(1) + fz*(e(2)-e(1))
               d1ersz(1) = d1e(1) + fz*(d1e(2)-d1e(1))
               d1ersz(2) = d1fz*(e(2)-e(1))
               d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
               d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
               Ec = Ec + eps*qwght(n)*rhoval*fac
               if (ldew) func(n) = func(n) + eps*rhoval*fac
               Amat(n,1) = Amat(n,1) + (eps + rhoval*d1edrho(1))*fac
               Amat(n,2) = Amat(n,2) + (eps + rhoval*d1edrho(2))*fac
#ifdef SECOND_DERIV
c              1 = rsrs, 2 = rsz, 3 = zz
               d2ersz(1) = d2e(1) + fz*(d2e(2)-d2e(1))
               d2ersz(2) = d1fz*(d1e(2)-d1e(1))
               d2ersz(3) = d2fz*(e(2)-e(1))
c              1 = aa, 2 = ab, 3 = bb
               d2edrho(1) = d2ersz(1)*d1rs*d1rs
     &                    + d2ersz(2)*d1rs*d1zeta(1)*2.d0
     &                    + d2ersz(3)*d1zeta(1)*d1zeta(1)
     &                    + d1ersz(1)*d2rs
     &                    + d1ersz(2)*d2zeta(1)
               d2edrho(2) = d2ersz(1)*d1rs*d1rs
     &                    + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
     &                    + d2ersz(3)*d1zeta(1)*d1zeta(2)
     &                    + d1ersz(1)*d2rs
     &                    + d1ersz(2)*d2zeta(2)
               d2edrho(3) = d2ersz(1)*d1rs*d1rs
     &                    + d2ersz(2)*d1rs*d1zeta(2)*2.d0
     &                    + d2ersz(3)*d1zeta(2)*d1zeta(2)
     &                    + d1ersz(1)*d2rs
     &                    + d1ersz(2)*d2zeta(3)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*fac
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*fac
#endif
            endif
         enddo
      endif
      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_perdew86.F"
#endif
