c
c     Coulomb attenuated PBE exchange functional 
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_camxpbe96(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_camxpbe96_d2(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#else
      Subroutine xc_camxpbe96_d3(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho, 
     &     Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, nq, ipol, Ex,
     &     qwght,ldew,func)
#endif
c
c$Id$
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
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Atmp, Ctmp, Etmp
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
      double precision A2tmp, C2tmp, C3tmp
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
#endif
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision F73
#endif
#ifdef THIRD_DERIV
      double precision F10d3, F23
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      parameter (F73=7.d0/3.d0)
#endif
#ifdef THIRD_DERIV
      parameter (F10d3=10.0d0/3.0d0, F23=2.0d0/3.0d0)
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
      double precision gpbe0,gpbe1,gpbe2,gpbe3
      double precision grpbe0,grpbe1,grpbe2,grpbe3
      double precision grevpbe0,grevpbe1,grevpbe2,grevpbe3
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
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-11-13): Electron density for the XC-third derivatives
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0.0d0
#endif
c
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
#endif
            if (lfac) then
               Etmp = rho43
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp =  2d0*rhom23
#endif
c Daniel (7-27-12): XC-third derivative, coming from LDA.
#ifdef THIRD_DERIV
               A3tmp = -4.0d0*rhom53
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
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 2d0*rho43*d1g(2)
c
#if defined(THIRD_DERIV) 
c Second derivatives
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
            A2tmp = A2tmp
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*4d0
            C3tmp = rho43*d2g(3)*8d0
c Third derivatives
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
            A3tmp = A3tmp
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp = 
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*8.0d0
            C5tmp =
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))*16.0d0
            C6tmp =
     1             (rho43*d3g(4))*32.0d0
c
            call xc_att_xc_d3(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp*fac
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp*fac
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp*fac
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
     1                              + C6tmp*fac
#elif defined(SECOND_DERIV)
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
            A2tmp = A2tmp
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*4d0
            C3tmp = rho43*d2g(3)*8d0
c
            call xc_att_xc_d2(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
#else
            call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp*fac
            if (ldew) func(n) = func(n) + Etmp*fac
            Amat(n,1) = Amat(n,1) + Atmp*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp*fac
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-11-13): Derivative of the electron density
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0.0d0
#endif
c
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
#endif
            if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp = 2d0*rhom23
#endif
c
#ifdef THIRD_DERIV
               A3tmp = -4.0d0*rhom53
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
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#ifdef THIRD_DERIV
c Second derivatives
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
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp = rho43*d2g(3)*0.5d0
c Third derivatives
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
            A3tmp = A3tmp 
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp =
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*2.0d0
            C5tmp = 
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))
            C6tmp = 
     1             (rho43*d3g(4))*0.5d0
c
            call xc_att_xc_d3(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp*fac
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp*fac
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp*fac
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
     1                              + C6tmp*fac
#elif defined(SECOND_DERIV)
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
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp = rho43*d2g(3)*0.5d0
c
            call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
#else
            call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp*fac
            if (ldew) func(n) = func(n) + Etmp*fac
            Amat(n,1) = Amat(n,1) + Atmp*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp*fac
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
c Daniel (4-11-13): Derivative of the electron density
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0.0d0
#endif
c
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
#endif
            if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp= 2.d0*rhom23
#endif
c
#ifdef THIRD_DERIV
               A3tmp = -4.0d0*rhom53
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
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#ifdef THIRD_DERIV
c Second derivatives
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
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp =  rho43*d2g(3)*0.5d0
c Third derivatives
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
            A3tmp = A3tmp
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp = 
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*2.0d0
            C5tmp =
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))
            C6tmp = 
     1             (rho43*d3g(4))*0.5d0
c
            call xc_att_xc_d3(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*fac
c
            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp*fac
            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp*fac
            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp*fac
            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
     1                              + C6tmp*fac
#elif defined(SECOND_DERIV)
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
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp =  rho43*d2g(3)*0.5d0
c
            call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*fac
#else
            call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp*fac
            if (ldew) func(n) = func(n) + Etmp*fac
            Amat(n,2) = Amat(n,2) + Atmp*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp*fac
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
#include "xc_camxpbe96.F"
#endif
c Daniel (4-11-13): For 3rd derivatives
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_camxpbe96.F"
#endif
