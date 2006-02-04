#include "dft2drv.fh"
!#define PRINTA 1
      Subroutine xc_cpkzb99(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id: xc_cpkzb99.F,v 1.9 2006-02-04 19:31:45 edo Exp $
c

c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).

      Implicit none
c
c
c     Input and other parameters
c
      
      integer ipol, nq

      double precision cfac
      logical lcfac, nlcfac
      logical  ldew
      double precision func(*)

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
      double precision delrho(nq,3,ipol), gammaval
      
c
c     Kinetic Energy Density
c
      double precision tau(nq,3,ipol)
 
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
      double precision rhoval,rhoa,rhob

c    first sigma term
      double precision  taun
      double precision  CCC
      parameter (CCC = 0.53d0) !cpkzb empirical parameter

c   Second call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision rho_t(3), delrho_t(3,2)
      double precision  tauNA,tauNB
c    
      double precision gam12,pbe,tauw,xx2,en,
     ,     tauwa,tauwb,xx2a,xx2b,dtwat2dg,dtwat2dn,
     ,     dtwbt2dg,dtwbt2dn
      double precision pbeup,dtwt2dn,decggadn,dtwt2dg,
     ,     delc,decggadg,drevdn,drevdg,drevdt,
     ,     dpbeupdn,dpbeupdg,atermn,btermn,atermg,btermg,
     ,     erevc,finaln,apartg,finalg,apartt,finalt
c
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)
c
      double precision drevdna,drevdnb,drevdgaa,drevdgbb,
     A     drevdta,drevdtb,finalgbb
      double precision delca,delcb,
     A     detiladga,detiladgb,detilbdga,detilbdgb,
     A     detiladna,detiladnb,detilbdna,detilbdnb
      double precision etildea,etildeb,gaa,gbb
      double precision fabup,fabdown
      double precision delrho_A(3,2), rho_A(3)
c
      fac = cfac
      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======

         do 12  n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 12

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

         tauN = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)

         tauw = 0.125d0*gammaval/rhoval

         xx2 = (tauw/tauN)**2.d0

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
         decggadn= dneGGAdn(1)

         dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)

         decggadg= dneGGAdg(1)

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
#ifdef PRINTA
        write(0,'(A,4F18.6)') ' x d p dt',
     A       xx2,dpbeupdg,pbeup,dtwt2dg
        write(0,*) ' ag bg ',atermg,btermg
#endif
        drevdg=atermg-btermg
        
        drevdt=(CCC*pbe-(1.d0+CCC)*pbeup)*xx2*(-2.d0/tauN)

        
       delc = -(1.d0 + CCC)*delc
       erevc = en + delc
         
       if(ldew) func(n) = func(n) + rhoval*erevc*fac
       Ec = Ec + rhoval*erevc*qwght(n)*fac

c     derivs wrt n
       finaln= rhoval*drevdn + erevc 

       Amat(n,1)=Amat(n,1)+(finaln)*fac

c     derivs wrt g
         apartg=rhoval*drevdg
         finalg=apartg

         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2d0*finalg*fac

c     derivs wrt t
         apartt=rhoval*drevdt

         finalt=apartt

         Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac

12     continue
#ifdef PRINTA
            write(0,*) ' Amat ',Amat(1,1),Amat(2,1)
            write(0,*) ' Cmat ',Cmat(1,1),Cmat(2,1)
            write(0,*) ' Mmat ',Mmat(1,1),Mmat(2,1)
#endif

c
c     open-shell
c
       else   !ipol=2 and do alpha beta cases
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
         tauna = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)
         taunb = tau(n,1,2) + tau(n,2,2) + tau(n,3,2)
         taun = tauna+taunb
         rhoa=rho(n,2)   
         rhob=rho(n,3)   
         tauw = 0.125d0*gammaval/rhoval

         xx2 = (tauw/tauN)**2.d0
         dtwt2dn = -2.d0*xx2/rhoval
         dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)

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
          etildea= pbeup
          detiladna = dneFSPdn(1)
          detiladnb = 0d0
          detiladga = dneFSPdg(D1_GAA)*0.5d0
          detiladgb = 0d0

c     n_sigma/n_total factor
       fabup=rhoa/rhoval
       gaa = delrho(n,1,1)*delrho(n,1,1) +
     &      delrho(n,2,1)*delrho(n,2,1) +
     &      delrho(n,3,1)*delrho(n,3,1)
       tauwa = gaa/(8d0*rhoa) 
       xx2a = (tauwa/tauna)**2.d0
       dtwat2dn = -2.d0*xx2a/rhoa
       dtwat2dg = 2.d0*xx2a/gaa

       delc= xx2a*etildea
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
       rho_A(3) = 0.d0          ! beta  equals zero
       delrho_A(1,1) = delrho_t(1,2) ! nabla n_up x  
       delrho_A(2,1) = delrho_t(2,2) ! nabla n_up y
       delrho_A(3,1) = delrho_t(3,2) ! nabla n_up z
       delrho_A(1,2) = 0.d0     ! set beta gradient to zero
       delrho_A(2,2) = 0.d0     ! set beta gradient to zero
       delrho_A(3,2) = 0.d0     ! set beta gradient to zero
       
       neFSP = 0.0d0            !Ec in PBE
       dneFSPdn(1) = 0.0d0      !Amat in PBE
       dneFSPdn(2) = 0.0d0      !Amat in PBE
       dneFSPdg(1) = 0.0d0      !Cmat in PBE
       dneFSPdg(2) = 0.0d0      !Cmat in PBE
       dneFSPdg(3) = 0.0d0      !Cmat in PBE
c
         call xc_cMpbe96(tol_rho, rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)

         pbeup = neFSP

c        functional deriv info below fffffffffffff

         etildeb= pbeup
         detilbdna=0d0
         detilbdnb = dneFSPdn(1)
         detilbdga=0d0
         detilbdgb = dneFSPdg(D1_GAA)*0.5d0

c     n_sigma/n_total factor
       fabdown=rhob/rhoval
       gbb = delrho(n,1,2)*delrho(n,1,2) +
     &      delrho(n,2,2)*delrho(n,2,2) +
     &      delrho(n,3,2)*delrho(n,3,2)
       tauwb = gbb/(8d0*rhob)
       xx2b = (tauwb/taunb)**2.d0
       dtwbt2dn = -2.d0*xx2b/rhob
       dtwbt2dg = 2.d0*xx2b/gbb
       delcb = -(1.d0 + CCC)*fabdown*xx2b*etildeb
       erevc = erevc + delcb

       if(ldew) func(n) = func(n) + rhoval*erevc*fac
       Ec = Ec + rhoval*erevc*qwght(n)*fac
cderivs
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(1)
       btermn=(1.d0+CCC)*(
     &      dtwat2dn*fabup*etildea +
     +      ( (xx2a*etildea - xx2b*etildeb)*rhob/(rhoval*rhoval) +
     +      xx2a*fabup*detiladna  ))
       drevdna=atermn - btermn
       
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(2)
       btermn=(1.d0+CCC)*(
     &      dtwbt2dn*fabdown*etildeb +
     +      ( (xx2b*etildeb-xx2a*etildea)*rhoa/(rhoval*rhoval) +
     +      xx2b*fabdown*detilbdnb) )
       drevdnb=atermn - btermn

c gaa
       atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*dneggadg(D1_GAA)
       btermg=(1.d0+CCC)*fabup*(
     *      xx2a*detiladga + fabup*etildea*dtwat2dg)
#ifdef PRINTA
        write(0,'(A,5F18.6)') ' x d p dt',
     A       xx2a,detiladga,etildea,dtwat2dg,fabup
        write(0,*) ' ag bg ',atermg,btermg
#endif
       drevdgaa=atermg-btermg

c gbb
       atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*dneggadg(D1_GBB)
       btermg=(1.d0+CCC)*fabdown*
     *      (xx2b*detilbdgb + fabdown*etildeb*dtwbt2dg)

       drevdgbb=atermg-btermg

       drevdta=-2.d0*xx2/tauN*CCC*pbe+
     +      2d0*xx2a/tauna*(1.d0+CCC)*(fabup*etildea)
       drevdtb=-2d0*xx2/tauN*CCC*pbe+
     +      2d0*xx2b/taunb*(1.d0+CCC)*(fabdown*etildeb)

c     derivs wrt n

       finaln= rhoval*drevdna + erevc
       Amat(n,1)=Amat(n,1)+(finaln)*fac

       finaln= rhoval*drevdnb + erevc 
       Amat(n,2)=Amat(n,2)+(finaln)*fac

c     derivs wrt gaa
         finalg=rhoval*drevdgaa
         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalg*fac

c     derivs wrt gbb
         finalgbb=rhoval*drevdgbb

         Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+finalgbb*fac

         Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+(finalg+finalgbb)*fac

c     derivs wrt t

         apartt=rhoval*drevdta
         finalt=apartt
         Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac

         apartt=rhoval*drevdtb
         finalt=apartt
         Mmat(n,2)=Mmat(n,2)+0.5d0*finalt*fac


20     continue

#ifdef PRINTA
            write(0,*) ' Amat ',Amat(1,1),Amat(2,1)
            write(0,*) ' CmatA ',Cmat(1,1),
     2           0.5d0*(2d0*Cmat(1,1)+Cmat(1,2)),
     C           0.5d0*(2d0*Cmat(2,1)+Cmat(2,2))
            write(0,*) ' CmatB ',Cmat(1,3),
     2           0.5d0*(2d0*Cmat(1,3)+Cmat(1,2)),
     C           0.5d0*(2d0*Cmat(2,3)+Cmat(2,2))
            write(0,*) ' Mmat ',Mmat(1,1),Mmat(2,1)
#endif

          endif



      return
      end

c

      Subroutine xc_cpkzb99_d2()
      call errquit(' not coded ',0,0)
      return
      end




