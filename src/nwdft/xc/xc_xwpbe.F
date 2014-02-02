c     wPBE exchange functional 
c
c     This functional is a long range corrected version of the PBE GGA
c     functional.  It is based on an exchange hole satisfying several
c     constraints of the exact functional.
c
c     References:
c     [a] E. Weintraub, T. M. Henderson, and G. E. Scuseria, J. Chem.
c           Theo. Comp. 5, 754 (2009).
c     [b] T. M. Henderson, B. G. Janesko, and G. E. Scuseria, J. Chem.
c            Phys. 128, 194105, (2008).
c
c     References (for PBE and its revisions):
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_xwpbe(tol_rho, fac, lfac, nlfac, rho, delrho,
     1                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_xwpbe_d2(tol_rho, fac, lfac, nlfac, rho, delrho,
     1                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     2                        qwght,ldew,func)
#else
      Subroutine xc_xwpbe_d3(tol_rho, fac, lfac, nlfac, rho, delrho,
     1                       Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     2                       nq, ipol, Ex, qwght, ldew, func)
#endif
c
c$Id$
c
      implicit none
c
#include "dft2drv.fh"
c Daniel (4-12-13): Header for XC-third derivatives
#include "dft3drv.fh"
#include "case.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision F73
#endif
c
c Daniel (4-12-13): XC-3rd derivative information for d3s
#ifdef THIRD_DERIV
      double precision F10d3, F23
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      parameter (F73=7.d0/3.d0)
#endif
c
#ifdef THIRD_DERIV
      parameter (F10d3=10.0d0/3.0d0, F23=2.d0/3.d0)
#endif
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     1      d, g, gp1, d1gs(2), gp2, d1gnu
      double precision nu, kf, d1nu
      double precision Fxhse,Fxhse10,Fxhse01
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision rhom23, gpp1, d2s(3),d2gs(3), gpp2, d2gnu
      double precision Fxhse20,Fxhse02,d2nu,gpp11,d2gsnu,Fxhse11
#endif
c
#ifdef THIRD_DERIV
      double precision rhom53, gppp1, gppp21, gppp12, gppp2
      double precision d3nu, d3s(4), d3gs(4)
      double precision Fxhse30, Fxhse21, Fxhse12, Fxhse03 
#endif
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
c Daniel (4-12-13): For information on the constants multiplied into 
c the partial derivatives, see xc_xpbe96.F
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13 
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c
c Daniel (4-12-13): Electron density for the XC-third derivatives
#ifdef THIRD_DERIV
            rhom53 = -F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2.0d0*rhom23*fac
#endif
c Daniel (4-13-13): XC-third derivative, coming from LDA.
#ifdef THIRD_DERIV
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
     1                              + 4.0d0*rhom53*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     1              delrho(n,2,1)*delrho(n,2,1) +
     2              delrho(n,3,1)*delrho(n,3,1)
            gam12 = dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 10
c

            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
#ifdef THIRD_DERIV
            call HSE08Fx_d3(ipol,rho(n,1),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11,fxhse30,fxhse21,
     2                      fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call HSE08Fx_d2(ipol,rho(n,1),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,1),s,fxhse,fxhse10,fxhse01)
#endif
c
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1)= Amat(n,1) + ( rho13*g 
     1                             + rho43*( d1gs(1) +d1gnu ) )*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2.0d0*rho43*d1gs(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*rrho*d1nu
c
            gpp1=fxhse20
            gpp2=fxhse02
            gpp11=fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     1                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2.0d0*gpp11*d1nu*d1s(1)

            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     1           + ( rhom23*g
     2             + 2.0d0*rho13*( d1gs(1) + d1gnu )
     3             + rho43*( d2gs(1) + d2gnu + d2gsnu ) )*fac*2.0d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     1           + ( rho13*d1gs(2)
     2             + rho43*d2gs(2) )*fac*4.0d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     1           + rho43*d2gs(3)*fac*8.0d0
#endif
c
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
      d3s(1) = -F10d3*d2s(1)*rrho
      d3s(2) = 0.5d0*d2s(1)/gamma
      d3s(3) = -F43*d2s(3)*rrho
      d3s(4) = -1.5d0*d2s(3)/gamma
c
      d3nu = -F73*d2nu*rrho
c
      d3gs(1) = d3nu*Fxhse01 
     1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 )  
     2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
     3        + d1nu*( d2nu*Fxhse02 
     4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
     5               + Fxhse11*d2s(1) 
     6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
     7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
     8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
     9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
     A        + Fxhse10*d3s(1)

      d3gs(2) = d2nu*d1s(2)*Fxhse11 
     1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
     2        + d1s(2)*Fxhse20*d2s(1) 
     3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
     4               + d1s(2)*d1s(1)*Fxhse21 ) 
     5        + Fxhse10*d3s(2) 
     6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
     7                 + d1s(2)*d1s(1)*Fxhse30 )

      d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
     1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 ) 
     2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

      d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2) 
     1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
      Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
     1              + ( rhom53*g
     2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
     3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
     4                + rho43*d3gs(1) )*fac*4.0d0
      Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1              + ( rhom23*d1gs(2)
     2                + 2.0d0*rho13*d2gs(2)
     3                + rho43*d3gs(2) )*fac*8.0d0
      Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) 
     1              + ( rho43*d3gs(3) + rho13*d2gs(3) )*fac*16.0d0
      Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
     1              + rho43*d3gs(4)*fac*32.0d0
#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = C*(2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-12-13): Derivative of the electron density
#ifdef THIRD_DERIV
            rhom53 = -F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,1) = Amat(n,1) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2.0d0*rhom23*fac
#endif
c
#ifdef THIRD_DERIV
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1                              + 4.0d0*rhom53*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     1              delrho(n,2,1)*delrho(n,2,1) +
     2              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 25
c
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
#ifdef THIRD_DERIV
            call HSE08Fx_d3(ipol,rho(n,2),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11,fxhse30,fxhse21,
     2                      fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call HSE08Fx_d2(ipol,rho(n,2),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,2),s,fxhse,fxhse10,fxhse01)
#endif
c
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01 
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + ( rho13*g 
     1                              + rho43*( d1gs(1) + d1gnu ) )*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1gs(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c     
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     1                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     1           + ( rhom23*g
     2             + 2.0d0*rho13*( d1gs(1) + d1gnu )
     3             + rho43*( d2gs(1) + d2gnu + d2gsnu ) )*fac*2.0d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     1           + ( rho13*d1gs(2)
     2             + rho43*d2gs(2) )*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     1           + rho43*d2gs(3)*fac*0.5d0
#endif
c
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
      d3s(1) = -F10d3*d2s(1)*rrho
      d3s(2) = 0.5d0*d2s(1)/gamma
      d3s(3) = -F43*d2s(3)*rrho
      d3s(4) = -1.5d0*d2s(3)/gamma
c
      d3nu = -F73*d2nu*rrho
c
      d3gs(1) = d3nu*Fxhse01 
     1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 ) 
     2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
     3        + d1nu*( d2nu*Fxhse02 
     4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
     5               + Fxhse11*d2s(1) 
     6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
     7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
     8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
     9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
     A        + Fxhse10*d3s(1)

      d3gs(2) = d2nu*d1s(2)*Fxhse11 
     1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
     2        + d1s(2)*Fxhse20*d2s(1) 
     3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
     4               + d1s(2)*d1s(1)*Fxhse21 ) 
     5        + Fxhse10*d3s(2) 
     6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
     7                 + d1s(2)*d1s(1)*Fxhse30 )

      d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
     1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 )   
     2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

      d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2)
     1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
      Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
     1              + ( rhom53*g
     2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
     3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
     4                + rho43*d3gs(1) )*fac*4.0d0
      Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1              + ( rhom23*d1gs(2)
     2                + 2.0d0*rho13*d2gs(2)
     3                + rho43*d3gs(2) )*fac*2.0d0
      Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) 
     1              + ( rho43*d3gs(3)
     2                + rho13*d2gs(3) )*fac
      Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
     1              + (rho43*d3gs(4))*fac*0.50d0
#endif
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c
#ifdef THIRD_DERIV
            rhom53 = -F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,2) = Amat(n,2) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2.0d0*rhom23*fac
#endif
c
#ifdef THIRD_DERIV
               Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     1                              + 4.0d0*rhom53*fac
#endif

            endif
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     1              delrho(n,2,2)*delrho(n,2,2) +
     2              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2.0d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 20
c
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
#ifdef THIRD_DERIV
            call HSE08Fx_d3(ipol,rho(n,3),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11,fxhse30,fxhse21,
     2                      fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call HSE08Fx_d2(ipol,rho(n,3),s,fxhse,fxhse10,fxhse01,
     1                      fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,3),s,fxhse,fxhse10,fxhse01)
#endif
c
            g = fxhse - 1.0d0
            gp1= fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,2) = Amat(n,2) + ( rho13*g 
     1                              + rho43*( d1gs(1) + d1gnu ) )*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1gs(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     1                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     1           + ( rhom23*g
     2             + 2.0d0*rho13*( d1gs(1) + d1gnu )
     3             + rho43*( d2gs(1) + d2gnu + d2gsnu ) )*fac*2.0d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     1           + ( rho13*d1gs(2)
     2             + rho43*d2gs(2) )*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     1           + rho43*d2gs(3)*fac*0.5d0
#endif
c
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
      d3s(1) = -F10d3*d2s(1)*rrho
      d3s(2) = 0.5d0*d2s(1)/gamma
      d3s(3) = -F43*d2s(3)*rrho
      d3s(4) = -1.5d0*d2s(3)/gamma
c
      d3nu = -F73*d2nu*rrho
c
      d3gs(1) = d3nu*Fxhse01 
     1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 ) 
     2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
     3        + d1nu*( d2nu*Fxhse02 
     4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
     5               + Fxhse11*d2s(1) 
     6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
     7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
     8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
     9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
     A        + Fxhse10*d3s(1)

      d3gs(2) = d2nu*d1s(2)*Fxhse11 
     1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
     2        + d1s(2)*Fxhse20*d2s(1) 
     3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
     4               + d1s(2)*d1s(1)*Fxhse21 ) 
     5        + Fxhse10*d3s(2) 
     6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
     7                 + d1s(2)*d1s(1)*Fxhse30 )

      d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
     1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 ) 
     2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

      d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2)
     1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) 
     1              + ( rhom53*g
     2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
     3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
     4                + rho43*d3gs(1) )*fac*4.0d0
      Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
     1              + ( rhom23*d1gs(2)
     2                + 2.0d0*rho13*d2gs(2)
     3                + rho43*d3gs(2) )*fac*2.0d0
      Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) 
     1              + ( rho43*d3gs(3) 
     2                + rho13*d2gs(3) )*fac
      Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) 
     1              + rho43*d3gs(4)*fac*0.50d0
#endif
c
 20      continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_xwpbe.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_xwpbe.F"
#endif
