c   
c$Id: xc_ctpss03.F,v 1.8 2006-02-03 01:45:25 edo Exp $
c
#include "dft2drv.fh"
c   
c$Id: xc_ctpss03.F,v 1.8 2006-02-03 01:45:25 edo Exp $
c
#include "dft2drv.fh"

c    Tao,Perdew, Staroverov, Scuseria exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density


      Subroutine xc_ctpss03(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)     
c     References:
c     [a] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         PRL 91, 146401 (2003).
c     [b] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         JCP 120, 6898 (2004).

      Implicit none
c
c
c     Input and other parameters
c
      integer ipol, nq

      double precision cfac
      logical  ldew,lcfac,nlcfac
      double precision func(*)  ! value of the functional [output]

      double precision fac
      double precision tol_rho
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
 
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol), gammaval, gam12
      
c
c     Kinetic Energy Density
c
      double precision tau(nq,3,ipol), tauN
 
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)

      integer n
      double precision rhoval

c    call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision rho_t(3), delrho_t(3,2)
      double precision  tauNA,tauNB
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)
      double precision delrho_A(3,2), rho_A(3)

c     TPSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs

      double precision THRD, FACE, F73, F43
      double precision zeta, ccc
      double precision tauw
      double precision xx2,xx3
      double precision pbe, en
      double precision dd
      double precision rhoa,rhob
      double precision pbeup, delc, erevc, erevsic,
     A     delca,delcb

c     derivsssssssssssssssssssssssss

      double precision  decggadn, dtwt2dn, dtwt3dn
      double precision  dpbeupdn, drevdn
      double precision  atermn, btermn
      double precision  finaln
      double precision  drevdg, dpbeupdg
      double precision  dtwt2dg, dtwt3dg, decggadg
      double precision  atermg, btermg
      double precision  apartg, cpartg(2), finalg
      double precision  finalgbb
      double precision drevdt, apartt, bpartt(2),finalt
      double precision drevdta,drevdtb

      double precision dzetadna, dzetadnb
      double precision dCCCdna, dCCCdnb
      double precision drevdna, drevdnb

      double precision drevdgaa
      double precision drevdgbb
      double precision dtwt3dt
      double precision etildea,etildeb,
     N     detiladna,detiladnb,detilbdna,detilbdnb,
     D     detiladga,detiladgb,detilbdga,detilbdgb
      double precision pi,fabup,fabdown
      double precision czeta0,czeta1
      parameter (dd = 2.8d0)
      parameter ( THRD = 1d0/3d0)
      parameter ( F73 = 7d0/3d0)
      parameter ( F43 = 4d0/3d0)
c
      czeta0(zeta)=0.53d0 + 0.87d0*zeta**2+
     +     0.50d0*zeta**4 + 2.26d0*zeta**6
      czeta1(zeta)= 2d0*0.87d0*zeta*+
     +     4d0*0.50d0*zeta**3 + 6d0*2.26d0*zeta**5
c
       
      face=(3d0*pi*pi)**THRD

      fac = cfac


      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======

         do 22  n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 22

         rhoval = rho(n,1)

C   set up values to call PBE subroutine
         rho_t(1) = rho(n,1)
c do delrho
         delrho_t(1,1) = delrho(n,1,1)
         delrho_t(2,1) = delrho(n,2,1)
         delrho_t(3,1) = delrho(n,3,1)
         gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &        delrho(n,2,1)*delrho(n,2,1) +
     &        delrho(n,3,1)*delrho(n,3,1)
         gam12=dsqrt(gammaval)
c
c     get E_GGA[rho,gamma]
c
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE
         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE

         call xc_cMpbe96(tol_rho,
     &        rho_t, delrho_t,
     &        dneGGAdn,dneGGAdg, 
     &        1, ipol, neGGA)
         pbe = neGGA

c
c        epGGA = n*(epsilon_c^GGA) / n =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe

         tauN = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)
!         taun=5d0

         CCC = 0.53d0           !since zeta=0

         tauw = 0.125d0*gammaval/rhoval

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0

         en = pbe*(1.d0 + CCC*xx2)
c
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system
c

         rho_A(1) = (0.5d0)*rho(n,1)   ! total   equals (1/2)n_tot
         rho_A(2) = (0.5d0)*rho(n,1)   ! alpha   equals (1/2)n_tot
         rho_A(3) = 0.d0               ! beta  equals zero
         delrho_A(1,1) = (0.5d0)*delrho_t(1,1) ! nabla n_up x  
         delrho_A(2,1) = (0.5d0)*delrho_t(2,1) ! nabla n_up y
         delrho_A(3,1) = (0.5d0)*delrho_t(3,1) ! nabla n_up z
     
         delrho_A(1,2) = 0.d0   ! set beta gradient to zero
         delrho_A(2,2) = 0.d0   ! set beta gradient to zero
         delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

c
c     get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         call xc_cMpbe96(tol_rho, rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)

         pbeup = neFSP

c        functional deriv info below fffffffffffff

         dtwt2dn = -2.d0*xx2/rhoval
         dtwt3dn = -3.d0*xx3/rhoval
         decggadn= dneGGAdn(1)

         dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)

         dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)
         decggadg= dneGGAdg(1)

       if (pbeup.lt.pbe) then
          delc= xx2*pbe
C  eps-tilda is eps_c
c     erev = egga*(1-xx2)
C  functional deriv info below fffffffffffffffff

          drevdn= -pbe*dtwt2dn+decggadn*(1.d0 - xx2)
          drevdg= -pbe*dtwt2dg+decggadg*(1.d0 -xx2)
          drevdt= 2.d0*pbe*xx2/tauN
       else
          delc= xx2*pbeup

C  eps-tilda is eps^FSP
C  functional deriv info below fffffffffffffffff

        dpbeupdn = 0.5d0*dneFSPdn(1)
c  above note the .5's.  you are taking deriv wrt total density n
c                        not deriv wrt n_up
        dpbeupdg = 0.25d0*dneFSPdg(1)
c  note .25 above is because you want gamma=deln_tot*deln_tot


        atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*decggadn
        btermn=(1.d0+CCC)*(xx2*dpbeupdn +
     &             pbeup*dtwt2dn)
        drevdn=atermn - btermn

          atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*decggadg
          btermg=(1.d0+CCC)*(xx2*dpbeupdg+pbeup*dtwt2dg)
          drevdg=atermg-btermg

         drevdt=(CCC*pbe-(1.d0+CCC)*pbeup)*xx2*(-2.d0/tauN)

       endif

       delc = -(1.d0 + CCC)*delc
       erevc = en + delc
       erevsic = erevc*(1.d0 + dd*erevc*xx3)
         
       if(ldew) func(n) = func(n) + rhoval*erevsic*fac
       Ec = Ec + rhoval*erevsic*qwght(n)*fac

c     derivs wrt n
       finaln= rhoval*drevdn + erevc +
     &      dd*(erevc*erevc*xx3 +
     +      rhoval*(xx3*2.d0*erevc*drevdn +
     +      erevc*erevc*dtwt3dn))

       Amat(n,1)=Amat(n,1)+(finaln)*fac

c     derivs wrt g
         apartg=rhoval*drevdg
         cpartg(1)=erevc*erevc*dtwt3dg
         cpartg(2)=xx3*2.d0*erevc*drevdg
         finalg=apartg+rhoval*dd*( cpartg(1)+cpartg(2) )

         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2d0*finalg*fac

c     derivs wrt t
         apartt=rhoval*drevdt
         bpartt(1)=-erevc*erevc*xx3*3.d0/tauN
         bpartt(2)=xx3*2.d0*erevc*drevdt

         finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))

         Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac

22     continue

      else !ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
c            
            if (rho(n,1).lt.tol_rho) goto 20

         rhoval = rho(n,1)

         rho_t(1) = rho(n,1)
         rho_t(2) = rho(n,2)
         rho_t(3) = rho(n,3)
         delrho_t(1,1) = delrho(n,1,1)
         delrho_t(2,1) = delrho(n,2,1)
         delrho_t(3,1) = delrho(n,3,1)
         delrho_t(1,2) = delrho(n,1,2)
         delrho_t(2,2) = delrho(n,2,2)
         delrho_t(3,2) = delrho(n,3,2)
         
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE (n,1)
         dneGGAdn(2) = 0.0d0   !Amat in PBE (n,2)
         dneGGAdg(1) = 0.0d0  !Cmat in PBE--aa
         dneGGAdg(2) = 0.0d0  !Cmat in PBE--ab
         dneGGAdg(3) = 0.0d0  !Cmat in PBE--bb

c
c     get E_GGA[rho,gamma]
c
         call xc_cMpbe96(tol_rho,
     &        rho_t, delrho_t,
     &        dneGGAdn,dneGGAdg, 
     &        1, ipol, neGGA)
         pbe = neGGA

c
c        epGGA = (epsilon_c^GGA)  =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe

         gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &        delrho(n,1,2)*delrho(n,1,2) +
     &        delrho(n,2,1)*delrho(n,2,1) +
     &        delrho(n,2,2)*delrho(n,2,2) +
     &        delrho(n,3,1)*delrho(n,3,1) +
     &        delrho(n,3,2)*delrho(n,3,2) +
     &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &        delrho(n,2,1)*delrho(n,2,2) +
     &        delrho(n,3,1)*delrho(n,3,2))
         gam12=dsqrt(gammaval)
         tauNa = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)
         tauNb = tau(n,1,2) + tau(n,2,2) + tau(n,3,2)
         taun = tauna+taunb
         rhoa=rho(n,2)   
         rhob=rho(n,3)   
         zeta=(rhoa-rhob)/rhoval
         CCC=czeta0(zeta)
         dzetadna=2d0*rhob/(rhoval*rhoval)
         dzetadnb=-2d0*rhoa/(rhoval*rhoval)
         dcccdna=czeta1(zeta)*dzetadna
         dcccdnb=czeta1(zeta)*dzetadnb

         tauw = 0.125d0*gammaval/rhoval

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0
         dtwt2dn = -2.d0*xx2/rhoval
         dtwt3dn = -3.d0*xx3/rhoval
         dtwt3dt = -3.d0*xx3/taun
         dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)
         dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)

         en = pbe*(1.d0 + CCC*xx2)

c
c     Alpha bit
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system for alpha spin
c     to get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         rho_A(1) = rhoa
         rho_A(2) = rhoa
         rho_A(3) = 0.d0               ! beta  equals zero
         delrho_A(1,1) = delrho_t(1,1) ! nabla n_up x  
         delrho_A(2,1) = delrho_t(2,1) ! nabla n_up y
         delrho_A(3,1) = delrho_t(3,1) ! nabla n_up z
         delrho_A(1,2) = 0.d0   ! set beta gradient to zero
         delrho_A(2,2) = 0.d0   ! set beta gradient to zero
         delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
         call xc_cMpbe96(tol_rho, rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)

         pbeup = neFSP

c        functional deriv info below fffffffffffff
       if (pbeup.lt.pbe) then
          etildea=pbe
          detiladna=dneggadn(1)
          detiladnb=dneggadn(2)
          detiladga =dneggadg(D1_GAA)
          detiladgb =dneggadg(D1_GBB)
       else
          etildea= pbeup
          detiladna = dneFSPdn(1)
          detiladnb = 0d0
          detiladga = dneFSPdg(D1_GAA)*0.5d0
          detiladgb = 0d0
       endif

c     n_sigma/n_total factor
       fabup=rhoa/rhoval

       delc= xx2*etildea
       delca = -(1.d0 + CCC)*fabup*delc 
       erevc = en + delca

c
c     Beta bit
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system for beta spin
c     to get E_GGA[rho_beta,0,gamma_beta,0]
c
         rho_A(1) = rhob
         rho_A(2) = rhob
         rho_A(3) = 0.d0               ! beta  equals zero
         delrho_A(1,1) = delrho_t(1,2) ! nabla n_up x  
         delrho_A(2,1) = delrho_t(2,2) ! nabla n_up y
         delrho_A(3,1) = delrho_t(3,2) ! nabla n_up z
         delrho_A(1,2) = 0.d0   ! set beta gradient to zero
         delrho_A(2,2) = 0.d0   ! set beta gradient to zero
         delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
         call xc_cMpbe96(tol_rho, rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)

         pbeup = neFSP

c        functional deriv info below fffffffffffff

       if (pbeup.lt.pbe) then
          etildeb=pbe
          detilbdna=dneggadn(1)
          detilbdnb=dneggadn(2)
          detilbdga =dneggadg(D1_GAA)
          detilbdgb =dneggadg(D1_GBB)
       else
          etildeb= pbeup
          detilbdna=0d0
          detilbdnb = dneFSPdn(1)
          detilbdga=0d0
          detilbdgb = dneFSPdg(D1_GAA)*0.5d0
       endif

c     n_sigma/n_total factor
       fabdown=rhob/rhoval

       delc= xx2*etildeb
       delcb = -(1.d0 + CCC)*fabdown*delc
       erevc = erevc + delcb

       erevsic = erevc*(1.d0 + dd*erevc*xx3)
         
       if(ldew) func(n) = func(n) + rhoval*erevsic*fac
       Ec = Ec + rhoval*erevsic*qwght(n)*fac
cderivs
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(1)+
     Z      pbe*xx2*dcccdna
       btermn=(1.d0+CCC)*(
     &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
     +      xx2*( (etildea - etildeb)*rhob/(rhoval*rhoval) +
     +      fabup*detiladna + fabdown*detilbdna) ) +
     Z      xx2*(fabup*etildea+fabdown*etildeb)*dcccdna
       drevdna=atermn - btermn
       
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(2)+
     Z      pbe*xx2*dcccdnb
       btermn=(1.d0+CCC)*(
     &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
     +      xx2*( (etildeb-etildea)*rhoa/(rhoval*rhoval) +
     +      fabup*detilbdna+fabdown*detilbdnb) )+
     +      xx2*(fabup*etildea+fabdown*etildeb)*dcccdnb
       drevdnb=atermn - btermn

c gaa
       atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*dneggadg(D1_GAA)
       btermg=(1.d0+CCC)*
     *      (xx2*(fabup*detiladga+fabdown*detiladgb)+
     +      (etildea*fabup+etildeb*fabdown)*dtwt2dg)
       drevdgaa=atermg-btermg

c gbb
       atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*dneggadg(D1_GBB)
       btermg=(1.d0+CCC)*
     *      (xx2*(fabup*detilbdga+fabdown*detilbdgb)+
     +      (etildea*fabup+etildeb*fabdown)*dtwt2dg)
       drevdgbb=atermg-btermg

       drevdta=-2.d0*xx2/tauN*
     *(CCC*pbe-(1.d0+CCC)*(fabup*etildea+fabdown*etildeb))
       drevdtb=-2d0*xx2/tauN*
     *(CCC*pbe-(1.d0+CCC)*(fabup*etildea+fabdown*etildeb))

c     derivs wrt n

       finaln= rhoval*drevdna + erevc +
     &      dd*(erevc*erevc*xx3 +
     +      rhoval*(xx3*2.d0*erevc*drevdna +
     +      erevc*erevc*dtwt3dn))
       Amat(n,1)=Amat(n,1)+(finaln)*fac

       finaln= rhoval*drevdnb + erevc +
     &      dd*(erevc*erevc*xx3 +
     +      rhoval*(xx3*2.d0*erevc*drevdnb +
     +      erevc*erevc*dtwt3dn))
       Amat(n,2)=Amat(n,2)+(finaln)*fac

c     derivs wrt gaa
         apartg=rhoval*drevdgaa
         cpartg(1)=erevc*erevc*dtwt3dg
         cpartg(2)=xx3*2.d0*erevc*drevdgaa
         finalg=apartg+rhoval*dd*( cpartg(1)+cpartg(2) )

         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalg*fac

c     derivs wrt gbb
         apartg=rhoval*drevdgbb
         cpartg(1)=erevc*erevc*dtwt3dg
         cpartg(2)=xx3*2.d0*erevc*drevdgbb
         finalgbb=apartg+rhoval*dd*( cpartg(1)+cpartg(2) )

         Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+finalgbb*fac

         Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+(finalg+finalgbb)*fac

c     derivs wrt t

         apartt=rhoval*drevdta
         bpartt(1)=erevc*erevc*dtwt3dt
         bpartt(2)=xx3*2.d0*erevc*drevdta
         finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))
         Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac

         apartt=rhoval*drevdtb
         bpartt(1)=erevc*erevc*dtwt3dt
         bpartt(2)=xx3*2.d0*erevc*drevdtb
         finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))
         Mmat(n,2)=Mmat(n,2)+0.5d0*finalt*fac


20     continue

      endif  !end ipol=2 case

      return
      end

c

      Subroutine xc_ctpss03_d2()
      call errquit(' not coded ',0,0)
      return
      end



