#ifndef SECOND_DERIV
      Subroutine xc_op(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ec, qwght,ldew,func)
#else
      Subroutine xc_op_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec,
     &                         qwght,ldew,func)
#endif
c
C$Id: xc_op.F,v 1.1 2006-01-06 01:58:46 edo Exp $
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, Ec
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
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
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
      double precision BETA,QABOP
      Parameter (BETA = 0.0042D0, QABOP=2.3670D0)
c
c References:
c    Tsuneda, Suzumura, Hirao, JCP 110, 10664 (1999)
c    Becke, Phys. Rev. A 38, 3098 (1988)
c
c***************************************************************************
c
      integer n
      double precision arcsinh, darcsinh
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg,
     &     dgdenom, t
      double precision kalpha,kbeta, rho13a, rho13b,rhoa,rhob
      double precision banb, hbab, hbabx
      double precision dhdab,dhdabx,dkadra,dkbdrb,dkadxa,dkbdxb,
     A     dbabdra,dbabdrb,dbabdga,dbabdgb,dkadga,dkbdgb
#ifdef SECOND_DERIV
      double precision rhom23, d2g, d2gdenom
#endif
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
      hbabx(x) = (1.5214d0*x + 0.5764d0)/
     /           (x**2*(x**2+1.1284d0*x+0.3183d0))
      dhdabx(x) = -(4.5642d0*x**3+5.7391d0*x**2+2.4355*x+0.3669d0)/
     /           (x**3*((x**2+1.1284d0*x+0.3183d0)**2))
c
c     Uniform electron gas constant
c
      C =  3d0*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            kalpha= C - 2d0*g
            
            banb = qabop * rho13 * kalpha *0.5d0

            hbab = hbabx(banb)

            Ec = Ec - 0.25d0*rho(n,1)**2*hbab*qwght(n)*fac
            if(ldew)func(n) = func(n) + 0.25d0*rho(n,1)**2*hbab*fac
            Amat(n,1) = Amat(n,1) + 0d0

c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0d0
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
            endif
c
#ifdef SECOND_DERIV
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (dabs(rho(n,1)).lt.tol_rho) goto 20
            if (dabs(rho(n,2)).ge.tol_rho)  then
c
c           Spin alpha:
c
               rhoa=rho(n,2)
               rho13a = abs(rhoa)**(1.d0/3.d0)*sign(1d0,rhoa)
               rho43 = rho13a**4 
               gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
               if (dsqrt(gamma).gt.tol_rho)then
                  x = sqrt(gamma) / rho43
                  gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
                  g = -BETA*x*x / gdenom
                  dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
                  dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
c
c     K \sigma from eq. 16 of JCP 1999
c
                  kalpha= C - 2d0*g
                  
                  dkadxa = (4d0*BETA*x + 12d0*BETA**2*x**2*
     *                 (arcsinh(x)-x*darcsinh(x)))/
     /                 (1+6d0*BETA*x*arcsinh(x))**2
                  
                  dkadra = -(4d0/3d0)*x*dkadxa/rhoa

                  dkadga = dkadxa*0.5d0/(rho43*sqrt(gamma))
               else
                  x = 0d0
               endif
            else
               rhoa=0d0
               rho13a=0d0
               rho43=0d0
            endif
            if(x.eq.0d0) then
               g=0d0
               kalpha= C
               dkadxa = 0d0
               dkadra = 0d0
               dkadga = 0d0
            endif
c     
c           Spin beta:
c
            if (dabs(rho(n,3)).ge.tol_rho) then
c
               rhob=rho(n,3)
               rho13b = abs(rhob)**(1.d0/3.d0)*sign(1d0,rhob)
               rho43 = rho13b**4 
               gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
               if (dsqrt(gamma).gt.tol_rho)then
                  x = sqrt(gamma) / rho43
                  gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
                  g = -BETA*x*x / gdenom
                  dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
                  dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c     
                  kbeta= C - 2d0*g
                  
                  dkbdxb = (4d0*BETA*x + 12d0*BETA**2*x**2*
     *                 (arcsinh(x)-x*darcsinh(x)))/
     /                 (1+6d0*BETA*x*arcsinh(x))**2
            
                  dkbdrb = -(4d0/3d0)*x*dkbdxb/rhob
            
                  dkbdgb = dkbdxb*0.5d0/(rho43*sqrt(gamma))
               else
                  x = 0d0
               endif
            else
               if(rho13a.eq.0) goto 20
               rhob=0d0
               rho13b=0d0
               rho43=0d0
               x=0d0
            endif
            if(x.eq.0d0) then
               kbeta=C
               dkbdxb = 0d0
               dkbdrb = 0d0
               dkbdgb=  0d0
            endif

            banb = qabop*(rho13a*kalpha*rho13b*kbeta)/
     /           (rho13a*kalpha+rho13b*kbeta)

            if(rhoa.ne.0) then
               dbabdra = qabop*
     *              (rho13b*kbeta/(rho13a*kalpha+rho13b*kbeta))**2 *
     *              (kalpha/(3d0*rho13a**2)+rho13a*dkadra)
            else
               dbabdra = 0d0
            endif

            if(rhob.ne.0) then
               dbabdrb = qabop*
     *              (rho13a*kalpha/(rho13a*kalpha+rho13b*kbeta))**2 *
     *              (kbeta/(3d0*rho13b**2)+rho13b*dkbdrb)
            else
               dbabdrb = 0d0
            endif

            dbabdga = qabop*rho13a*
     *           (rho13b*kbeta/(rho13a*kalpha+rho13b*kbeta))**2 *
     *           dkadga
            dbabdgb = qabop*rho13b*
     *           (rho13a*kalpha/(rho13a*kalpha+rho13b*kbeta))**2 *
     *           dkbdgb
            if(banb.ne.0) then
               hbab = hbabx(banb)
               dhdab = dhdabx(banb)
            else
               hbab = 0d0
               dhdab = 0d0
            endif

            Ec = Ec - rhoa*rhob*hbab*qwght(n)*fac
            if (ldew) func(n) = func(n) - rhoa*rhob*hbab*fac
            Amat(n,1) = Amat(n,1) - rhob*hbab - rhoa*rhob*dhdab*dbabdra
            Amat(n,2) = Amat(n,2) - rhoa*hbab - rhoa*rhob*dhdab*dbabdrb
c
c
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - rhoa*rhob*dhdab*dbabdga
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) - rhoa*rhob*dhdab*dbabdgb
               
c
#ifdef SECOND_DERIV
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,3)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
c
 20      continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_op.F"
#endif
