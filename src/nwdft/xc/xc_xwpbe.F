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
#ifndef SECOND_DERIV
      Subroutine xc_xwpbe(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_xwpbe_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#endif
c
      implicit none
c
#include "dft2drv.fh"
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
      double precision amat(nq,ipol), cmat(nq,*)
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
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
     &      d, g, gp1, d1gs(2), gp2, d1gnu
      double precision nu, kf, d1nu
      double precision Fxhse,Fxhse10,Fxhse01
#ifdef SECOND_DERIV
      double precision rhom23, gpp1, d2s(3),d2gs(3), gpp2, d2gnu
      double precision Fxhse20,Fxhse02,d2nu,gpp11,d2gsnu,Fxhse11
      double precision d1nud1rhonu,d1nud1rhos,d1sd1rhos,d1sd1rhonu
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
            kf = -1d0*pi*rho13 
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
               Amat(n,1) = Amat(n,1) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
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
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
#ifdef SECOND_DERIV
            call HSE08Fx_d2(ipol,rho(n,1),s,fxhse,fxhse10,fxhse01,
     &                      fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,1),s,fxhse,fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1)= Amat(n,1)+(rho13*g+rho43*(d1gs(1)+d1gnu))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1gs(2)*fac
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*rrho*d1nu
c
c            d1sd1rhos = -F43*rrho
c            d1sd1rhonu = -F13*2d0*cam_omega/gam12
c            d1nud1rhonu = -F13*rrho
c            d1nud1rhos = -F43*0.5d0*gam12*rrho**2d0/cam_omega
c
            gpp1=fxhse20
            gpp2=fxhse02
            gpp11=fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
c            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1) + 
c     &                gp1*d1sd1rhos*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
c            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu + 
c     &              gp2*d1nud1rhonu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c            d2gsnu = 2d0*gpp11*d1nu*d1s(1) + gp1*d1nud1rhos*d1nu +
c     &               gp2*d1sd1rhonu*d1s(1)

            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac*4d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2gs(3)*fac*8d0
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,1) = Amat(n,1) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
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
#ifdef SECOND_DERIV
            call HSE08Fx_d2(ipol,rho(n,2),s,fxhse,fxhse10,fxhse01,
     &                     fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,2),s,fxhse,fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01 
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*(d1gs(1)+d1gnu))*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1gs(2)*fac
#ifdef SECOND_DERIV
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
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2gs(3)*fac*0.5d0
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
               Amat(n,2) = Amat(n,2) + rho13*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*rhom23*fac
#endif
            endif
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
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
#ifdef SECOND_DERIV
            call HSE08Fx_d2(ipol,rho(n,3),s,fxhse,fxhse10,fxhse01,
     &                      fxhse20,fxhse02,fxhse11)
#else
            call HSE08Fx(ipol,rho(n,3),s,fxhse,fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1= fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,2) = Amat(n,2) + (rho13*g+rho43*(d1gs(1)+d1gnu))*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1gs(2)*fac
#ifdef SECOND_DERIV
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
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac*2d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           + rho43*d2gs(3)*fac*0.5d0
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
