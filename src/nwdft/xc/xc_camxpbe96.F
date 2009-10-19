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
#ifndef SECOND_DERIV
      Subroutine xc_camxpbe96(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_camxpbe96_d2(whichf,
     W     tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#endif
c
c$Id: xc_camxpbe96.F,v 1.3 2008/12/09 01:42:28 niri Exp $
c
      implicit none
c
#include "dft2drv.fh"
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
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
      double precision A2tmp, C2tmp, C3tmp
#endif
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
      parameter (F73=7.d0/3.d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
#ifdef SECOND_DERIV
      double precision rhom23, d2s(3), gpp, d2g(3)
#endif
      double precision gpbe0,gpbe1,gpbe2
      double precision grpbe0,grpbe1,grpbe2
      double precision grevpbe0,grevpbe1,grevpbe2
      gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
      gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
      gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
      grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
      grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
      grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
     /     (1d0+umkrev*s*s)**2
      grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
      grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
      grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = rho43*fac
               Atmp = rho13*fac
#ifdef SECOND_DERIV
               A2tmp =  2d0*rhom23*fac
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
            Etmp = Etmp + rho43*g*fac
            Atmp = Atmp + (rho13*g+rho43*d1g(1))*fac
            Ctmp = 2d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
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
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))*fac*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*fac*4d0
            C3tmp = rho43*d2g(3)*fac*8d0
c
            call xc_att_xc_d2(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
            call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if (ldew) func(n) = func(n) + Etmp
            Amat(n,1) = Amat(n,1) + Atmp
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = rho43*fac*0.5d0
               Atmp = rho13*fac
#ifdef SECOND_DERIV
               A2tmp = 2d0*rhom23*fac
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
            Etmp = Etmp + rho43*g*fac*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))*fac
            Ctmp = 0.5d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
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
     &           + rho43*d2g(1))*fac*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*fac
            C3tmp = rho43*d2g(3)*fac*0.5d0
c
            call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
            call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if (ldew) func(n) = func(n) + Etmp
            Amat(n,1) = Amat(n,1) + Atmp
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = rho43*fac*0.5d0
               Atmp = rho13*fac
#ifdef SECOND_DERIV
               A2tmp= 2.d0*rhom23*fac
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
            Etmp = Etmp + rho43*g*fac*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))*fac
            Ctmp = 0.5d0*rho43*d1g(2)*fac
#ifdef SECOND_DERIV
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
     &           + rho43*d2g(1))*fac*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*fac
            C3tmp =  rho43*d2g(3)*fac*0.5d0
c
            call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
#else
            call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if (ldew) func(n) = func(n) + Etmp
            Amat(n,2) = Amat(n,2) + Atmp
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp
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
