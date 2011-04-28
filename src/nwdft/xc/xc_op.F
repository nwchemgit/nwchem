#ifndef SECOND_DERIV
      Subroutine xc_op(tol_rho, whichf, 
     &     fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ec, qwght,ldew,func)
#else
      Subroutine xc_op_d2(tol_rho, whichf, 
     &     fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec,
     &                         qwght,ldew,func)
#endif
c
C$Id$
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, Ec
      character*4 whichf
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
      double precision QABOP,QAB88OP,QABPBOP
      Parameter (QAB88OP=2.3670D0,QABPBOP=2.3789D0)
c
c References:
c    Tsuneda, Suzumura, Hirao, JCP 110, 10664 (1999)
c    Tsuneda, Suzumura, Hirao, JCP 111, 5656 (1999)
c
c***************************************************************************
c
      integer n
      double precision rho13, rho43, gamma, x
      double precision kalpha,kbeta, rho13a, rho13b,rhoa,rhob
      double precision banb, hbab, hbabx
      double precision dhdab,dhdabx,dkadra,dkbdrb,dkadxa,dkbdxb,
     A     dbabdra,dbabdrb,dbabdga,dbabdgb,dkadga,dkbdgb,
     A     dbabdka,dbabdkb
c
      hbabx(x) = (1.5214d0*x + 0.5764d0)/
     /           (x**2*(x**2+1.1284d0*x+0.3183d0))
      dhdabx(x) = -(4.5642d0*x**4+5.7391d0*x**3+
     +     2.4355*x**2+0.3669d0*x)/
     /           ((x**4+1.1284d0*x**3+0.3183d0*x**2)**2)
c
      if(whichf.eq.'be88') then
         QABOP=QAB88OP
      endif
      if(whichf.eq.'pb96') then
         QABOP=QABPBOP
      endif
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rhoa=rho(n,1)*0.5d0
            rho13a = (rhoa)**(1.d0/3.d0)
            rho43 = rho13a**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gamma = 0.25d0 * gamma
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
               call xc_kop(tol_rho,whichf,x,
     &              kalpha, dkadxa)
               dkadra = -(4d0/3d0)*x*dkadxa/rhoa
               dkadga = (dkadxa/rho43)*0.5d0/dsqrt(gamma)
            else
               x=0d0
               call xc_kop(tol_rho,whichf,x,
     &              kalpha, dkadxa)
               dkadra = 0d0
               dkadga = 0d0
            endif
c
c
          
            banb = qabop * rho13a * kalpha *0.5d0

            if(banb.ne.0) then
               dbabdra = banb*0.5d0*
     /              (1d0/(3d0*rhoa)+dkadra/kalpha)

               dbabdga = banb/kalpha*dkadga*0.5d0
               
               hbab = hbabx(banb)
               dhdab = dhdabx(banb)
            else
               dbabdra =0d0
               dbabdga =0d0
               hbab = 0d0
               dhdab = 0d0
            endif

            Ec = Ec - rhoa**2*hbab*qwght(n)*fac
            if(ldew)func(n) = func(n) - rhoa**2*hbab*fac
            Amat(n,1) = Amat(n,1) - 
     -           (rhoa*hbab + rhoa**2*dhdab*dbabdra)*fac

c
            if (x.gt.tol_rho) then
                Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
     -              rhoa**2*dhdab*dbabdga*fac
             endif
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,2).ge.tol_rho*0.5d0)  then
c
c           Spin alpha:
c
               rhoa=rho(n,2)
               rho13a = rhoa**(1.d0/3.d0)
               rho43 = rho13a*rhoa
               gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
               if (dsqrt(gamma).gt.tol_rho)then
                  x = dsqrt(gamma) / rho43
                  call xc_kop(tol_rho,whichf,x,
     &                 kalpha, dkadxa)
                  
                  dkadra = -(4d0/3d0)*x*dkadxa/rhoa
                  dkadga = dkadxa*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               rhoa=0d0
               rho13a=0d0
               x = 0d0
            endif
            if(x.eq.0d0) then
               call xc_kop(tol_rho,whichf,x,
     &              kalpha, dkadxa)
               dkadra = 0d0
               dkadga = 0d0
            endif
c     
c           Spin beta:
c
            if (rho(n,3).ge.tol_rho*0.5d0) then
c
               rhob=rho(n,3)
               rho13b = rhob**(1.d0/3.d0)
               rho43 = rho13b*rhob
               gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
               if (dsqrt(gamma).gt.tol_rho)then
                  x = dsqrt(gamma) / rho43
                  call xc_kop(tol_rho,whichf,x,
     &                 kbeta, dkbdxb)
            
                  dkbdrb = -(4d0/3d0)*x*dkbdxb/rhob
                  dkbdgb = dkbdxb*0.5d0/(rho43*dsqrt(gamma))
               else
                  x = 0d0
               endif
            else
               if(rho13a.eq.0) goto 20
               rhob=0d0
               rho13b=0d0
               x=0d0
            endif
            if(x.eq.0d0) then
               call xc_kop(tol_rho,whichf,x,
     &              kbeta, dkbdxb)
               dkbdrb = 0d0
               dkbdgb=  0d0
            endif

            banb = qabop*(rho13a*kalpha*rho13b*kbeta)/
     /           (rho13a*kalpha+rho13b*kbeta)

            if(banb.ne.0) then
               dbabdra = banb*kbeta*rho13b/
     /              (rho13a*kalpha+rho13b*kbeta)*
     /              (1d0/(3d0*rhoa)+dkadra/kalpha)
               dbabdrb = banb*kalpha*rho13a/
     /              (rho13a*kalpha+rho13b*kbeta)*
     /              (1d0/(3d0*rhob)+dkbdrb/kbeta)

               dbabdga = banb*rho13b*kbeta/
     /              ((rho13a*kalpha+rho13b*kbeta)*kalpha)*
     *              dkadga
               dbabdgb = banb*rho13a*kalpha/
     /              ((rho13a*kalpha+rho13b*kbeta)*kbeta)*
     *              dkbdgb

               hbab = hbabx(banb)
               dhdab = dhdabx(banb)
            else
               dbabdra =0d0
               dbabdrb =0d0
               dbabdga =0d0
               dbabdgb =0d0
               hbab = 0d0
               dhdab = 0d0
            endif

            Ec = Ec - rhoa*rhob*hbab*qwght(n)*fac
            if (ldew) func(n) = func(n) - rhoa*rhob*hbab*fac
            Amat(n,1) = Amat(n,1) - 
     -           (rhob*hbab + rhoa*rhob*dhdab*dbabdra)*fac
            Amat(n,2) = Amat(n,2) - 
     -           (rhoa*hbab + rhoa*rhob*dhdab*dbabdrb)*fac
c
c
            if (x.gt.tol_rho) then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - 
     -              rhoa*rhob*dhdab*dbabdga*fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) - 
     -              rhoa*rhob*dhdab*dbabdgb*fac
            endif
               
c
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
