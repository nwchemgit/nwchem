c     PBE exchange functional 
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine xc_xpbe96(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho,
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
      Subroutine xc_xpbe96_d2(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho,
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#else
      Subroutine xc_xpbe96_d3(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho,
     &                        Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &                        nq, ipol, Ex, qwght,ldew,func)
#endif
c
      implicit none
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c      
      character*4 whichf
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
      double precision amat(nq,ipol), cmat(nq,*)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13, F23
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision F73
#endif
#ifdef THIRD_DERIV
      double precision F10d3
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0, F23=2.0d0/3.0d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      parameter (F73=7.d0/3.d0)
#endif
#ifdef THIRD_DERIV
      parameter (F10d3=10.0d0/3.0d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision rhom23, d2s(3), gpp, d2g(3)
#endif
#ifdef THIRD_DERIV
      double precision d3s(4), d3g(4), rhom53, gppp
#endif
      double precision gpbe0, gpbe1, gpbe2, gpbe3
      double precision grpbe0, grpbe1, grpbe2, grpbe3
      double precision grevpbe0, grevpbe1, grevpbe2, grevpbe3
c Original PBE
      gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
      gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
      gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
      gpbe3(s)= 24.0d0*umk*um*s*
     1  (2.0d0*umk*s*s/(1.0d0+umk*s*s)-1.0d0)/(1.0d0+umk*s*s)**3
c revPBE by Zhang et al.
      grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
      grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
      grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
     /     (1d0+umkrev*s*s)**2
      grevpbe3(s)= 24.0d0*umkrev*um*s*
     1  (2.0d0*umkrev*s*s/(1.0d0+umkrev*s*s)-1.0d0)/
     2  (1.0d0+umkrev*s*s)**3
c RPBE by Hammer et al.
      grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
      grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
      grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
      grpbe3(s)= -4.0d0*umk*um*s*exp(-umk*s*s)*(3d0-2d0*umk*s*s)
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
c Daniel (9-28-12): There are somewhat mysterious coefficients involved
c in the evaluation of the functional and its derivatives.  We must
c recall that the exchange energy is always written based on the
c spin-scaling relationship for exchange:
c
c Ex[rho] = Ex[rho_a,rho_b] = 0.5*( Ex[2*rho_a] + Ex[2*rho_b] )
c 
c Thus, the electron density is always written:
c rho -> 2*rho_s
c gamma -> 4*gamma_ss
c
c Rationalization for the coefficients is mathematically justified below:
c
c ----------------------------
c Amat       -> 0.5*2 = 1
c Cmat       -> 0.5*4 = 2
c ----------------------------
c Amat2      -> 0.5*2*2 = 2
c Cmat2(rg)  -> 0.5*2*4 = 4
c Cmat2(gg)  -> 0.5*4*4 = 8
c ----------------------------
c Amat3      -> 0.5*2*2*2 = 4
c Cmat3(rrg) -> 0.5*2*2*4 = 8
c Cmat3(rgg) -> 0.5*2*4*4 = 16
c Cmat3(ggg) -> 0.5*4*4*4 = 32
c ----------------------------
c
c If, instead, the author of this code had decided to divide the total
c density (rho(n,1)) by 2 in constructing the density and gamma, those
c coefficients would be unnecessary.
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c#ifdef THIRD_DERIV
c            write(6,*) 'rho', rho(n,1)
c#endif
            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
#ifdef THIRD_DERIV
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
     1                              - 4.0d0*rhom53*fac
#endif 
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 10
c

            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c
c Daniel (7-27-12): gp is the derivative of the rational function,
c or whatever the function in the revision is.
c First derivatives of the enhancement factor
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1g(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
c Second derivatives of the enhancement factor
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac*4d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*fac*8d0
#endif
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1) 
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac*4.0d0
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac*8.0d0
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
     1           + (rho13*d2g(3) 
     2           +  rho43*d3g(3))*fac*16.0d0
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
     1           + (rho43*d3g(4))*fac*32.0d0
#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c Daniel (9-28-12): There are somewhat mysterious coefficients involved
c in the evaluation of the functional and its derivatives.  We must
c recall that the exchange energy is always written based on the
c spin-scaling relationship for exchange:
c
c Ex[rho] = Ex[rho_a,rho_b] = 0.5*( Ex[2*rho_a] + Ex[2*rho_b] )
c 
c Thus, the electron density is always written:
c rho -> 2*rho_s
c gamma -> 4*gamma_ss
c
c It seems like the derivatives should be correctly balanced by the
c following coefficients:
c
c -----------------------------
c Amat       -> 0.5*2 = 1
c Cmat       -> 0.5*1 = 0.5
c -----------------------------
c Amat2      -> 0.5*2*2 = 2
c Cmat2(rg)  -> 0.5*2*1 = 1
c Cmat2(gg)  -> 0.5*1*1 = 0.5
c -----------------------------
c Amat3      -> 0.5*2*2*2 = 4
c Cmat3(rrg) -> 0.5*2*2*1 = 2
c Cmat3(rgg) -> 0.5*2*1*1 = 1
c Cmat3(ggg) -> 0.5*1*1*1 = 0.5
c -----------------------------
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,1) = Amat(n,1) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
#ifdef THIRD_DERIV
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
     1                              - 4.0d0*rhom53*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c

            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c Daniel (9-28-12): Factors of 2 are inconsistent with the restricted
c calculations because a gam12 is doubled above and  
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1g(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2g(3)*fac*0.5d0
#endif
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac*4.0d0
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac*2.0d0
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
     1           + (rho13*d2g(3)
     2           +  rho43*d3g(3))*fac
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
     1           + (rho43*d3g(4))*fac*0.5d0
#endif
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,2) = Amat(n,2) + rho13*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*rhom23*fac
#endif
#ifdef THIRD_DERIV
               Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     1                              - 4.0d0*rhom53*fac
#endif
            endif
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,2) = Amat(n,2) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1g(2)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           + rho43*d2g(3)*fac*0.5d0
#endif
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac*4.0d0
            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac*2.0d0
            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB)
     1           + (rho13*d2g(3)
     2           +  rho43*d3g(3))*fac
            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
     1           + (rho43*d3g(4))*fac*0.5d0
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
#include "xc_xpbe96.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_xpbe96.F"
#endif
c $Id$
