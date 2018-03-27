c     Perdew-Wang 86 exchange functional 
c
c     References:
c     [a] J. P. Perdew and W. Yue (actually, Y. Wang)
c         Phys. Rev. B 33 (1986) 8800
c
#ifndef SECOND_DERIV
      Subroutine xc_xpw86(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_xpw86_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#endif
c
c$Id$
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision tol_rho ! zero-density crit. (input)
      double precision fac ! energy multiplier (input)
      logical lfac, nlfac ! local and non-local labels (input)
      logical ldew        ! compute funct(n) (input)
      double precision Ex ! exchange energy (output)
      integer nq  ! number of grid points (input)
      integer ipol ! spinpolarization flag (input)
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2) ! charge density (in)
      double precision delrho(nq,3,ipol) ! gradient of rho (in)
      double precision qwght(nq) ! quadrature weights (in)
      double precision amat(nq,ipol) ! dex/drho matrix
      double precision cmat(nq,*)  ! 2*dex/(dgrho)
#ifdef SECOND_DERIV
c     second derivatives? dkdc for now.
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2) 
#endif

c      constants for PW86
      double precision a, b, c, m1, Ax, F43, s_pfac
      parameter(a = 1.296d0, b = 14d0, c = 0.2d0) 
      parameter(m1 = 1d0/15d0)
      parameter(s_pfac = 6.18733545256027d0)
      parameter(Ax = -0.738558766382022d0)
      parameter(F43 = 4d0/3d0)

      double precision pi, axr43, rh
      double precision s, s2, s3, s4, s5, s6
      double precision fs, sx, drho, drho2, dfds,oneminus
      integer n
c


      pi = acos(-1.d0)

c     spin-restricted
      if (ipol.eq.1) then
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

            rh = rho(n,1)
            drho2 = (delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))
            drho = dsqrt(drho2)
            if (drho.gt.tol_rho) then
            s = drho / (s_pfac*rh**F43)
            s2 = s * s
            s3 = s2 * s
            s4 = s3 * s
            s5 = s4 * s
            s6 = s5 * s
            fs = (1 + a*s2 + b*s4 + c*s6)**m1
            dfds = (2*a*s + 4*b*s3 + 6*c*s5)/(15.d0*fs**14)
            oneminus=1d0-dfds/fs*s
            else
               s=0d0
               fs=1d0
               dfds=0d0
               oneminus=1d0
            endif

            axr43 = Ax * rh**F43
            sx = axr43 * fs
            Ex = Ex + sx * qwght(n) * fac
            if (ldew) func(n) = func(n) + sx * fac

c     GC: dex/drho
            Amat(n,1) = Amat(n,1) + F43*sx/rh * oneminus
            if(drho.gt.tol_rho) then
c     GC: 2*dex/d(grho2)
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sx/fs*dfds*s/drho2
            endif

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif
 10      continue
c
      else
c     spin unrestricted
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 100 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 100

c           alpha
            if (rho(n,2).lt.tol_rho) goto 150
            rh = rho(n,2) * 2
            drho2 = (delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)) * 4
            drho = dsqrt(drho2) 
            s = drho / (s_pfac*rh**F43)
            s2 = s * s
            s3 = s2 * s
            s4 = s3 * s
            s5 = s4 * s
            s6 = s5 * s

            axr43 = Ax * rh**F43
            fs = (1 + a*s2 + b*s4 + c*s6)**m1
            sx = axr43 * fs
            Ex = Ex + 0.5d0 * sx * qwght(n) * fac 
            if (ldew) func(n) = func(n) + 0.5d0 * sx * fac

            dfds = (2*a*s + 4*b*s3 + 6*c*s5)/(15.d0*fs**14)
c     GC: dex/drho
            Amat(n,1) = Amat(n,1) + F43*sx/rh * (1d0-dfds/fs*s)
c     GC: 2*dex/d(grho2)
            if(drho2.gt.tol_rho) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sx/fs*dfds*s/drho2
            endif

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif

 150        continue
c           beta
            if (rho(n,3).lt.tol_rho) goto 100
            rh = rho(n,3) * 2
            drho2 = (delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)) * 4
            drho = dsqrt(drho2) 
            s = drho / (s_pfac*rh**F43)
            s2 = s * s
            s3 = s2 * s
            s4 = s3 * s
            s5 = s4 * s
            s6 = s5 * s

            axr43 = Ax * rh**F43
            fs = (1 + a*s2 + b*s4 + c*s6)**m1
            sx = axr43 * fs
            Ex = Ex + 0.5d0 * sx * qwght(n) * fac 
            if (ldew) func(n) = func(n) + 0.5d0 * sx * fac

            dfds = (2*a*s + 4*b*s3 + 6*c*s5)/(15.d0*fs**14)
c     GC: dex/drho
            Amat(n,2) = Amat(n,2) + F43*sx/rh * (1d0-dfds/fs*s)
c     GC: 2*dex/d(grho2)
            if(drho2.gt.tol_rho) then
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sx/fs*dfds*s/drho2
            endif

c     xxxx missing
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + 0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + 0
#endif

 100      continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_xpw86.F"
#endif
