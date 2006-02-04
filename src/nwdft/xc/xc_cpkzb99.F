#include "dft2drv.fh"
      Subroutine xc_cpkzb99(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id: xc_cpkzb99.F,v 1.7 2006-02-04 05:17:24 edo Exp $
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
      double precision dummy(1)

      double precision cfac
      logical lcfac, nlcfac
      logical  ldew
      double precision func(*)

      logical lfac, nlfac
      double precision fac
      double precision tol_rho

c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
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
      double precision  tauN,  Weiz
      double precision OE, CCC, mt
      parameter (OE = 0.125d0) !one eighth
      parameter (CCC = 0.53d0) !cpkzb empirical parameter

c   Second call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision  facG, Gqwght, epGGA
      logical  lfacG, nlfacG
      integer Gnq
      double precision rho_t(3), delrho_t(3,2)
      double precision  gammavalA, gammavalB
      double precision  tauNA,tauNB,  WeizA, WeizB
      double precision  sumTau, sumWeiz

c     spectral representation terms
      double precision aterma1, aterma2, aterma3
      double precision atermb1, atermb2, atermb3
      double precision aterm1, aterm2, aterm3
      double precision cterm1, cterm2, cterm3
      double precision mterm
      double precision ctermb1, ctermb2, ctermb3
      double precision gam12,pbe,tauw,xx2,xx3,en
      double precision pbeup,dtwt2dn,dtwt3dn,decggadn,dtwt2dg,
     ,     delc,dtwt3dg,decggadg,drevdn,drevdg,drevdt,
     ,     dpbeupdn,dpbeupdg,atermn,btermn,atermg,btermg,
     ,     erevc,erevsic,finaln,apartg,finalg,apartt,finalt
 

c  FULLY SPIN POLARIZED SECTION	
      double precision   rhovalF, epFSP
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)

      double precision  facF, Fqwght
      logical  lfacF, nlfacF
      double precision drevdna,drevdnb,drevdgaa,drevdgbb,
     A     drevdta,drevdtb,finalgbb
      double precision delca,delcb,
     A     detiladga,detiladgb,detilbdga,detilbdgb,
     A     detiladna,detiladnb,detilbdna,detilbdnb
      double precision dtwt3dt,etildea,etildeb
      double precision fabup,fabdown
      double precision delrho_A(3,2), rho_A(3)


      fac = cfac
      lfac = lcfac
      nlfac = nlcfac


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

        
       delc = -(1.d0 + CCC)*delc
       erevc = en + delc
       erevsic = erevc
         
       if(ldew) func(n) = func(n) + rhoval*erevsic*fac
       Ec = Ec + rhoval*erevsic*qwght(n)*fac

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
         tauNa = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)
         tauNb = tau(n,1,2) + tau(n,2,2) + tau(n,3,2)
         taun = tauna+taunb
         rhoa=rho(n,2)   
         rhob=rho(n,3)   
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
          etildea= pbeup
          detiladna = dneFSPdn(1)
          detiladnb = 0d0
          detiladga = dneFSPdg(D1_GAA)*0.5d0
          detiladgb = 0d0

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

          etildeb= pbeup
          detilbdna=0d0
          detilbdnb = dneFSPdn(1)
          detilbdga=0d0
          detilbdgb = dneFSPdg(D1_GAA)*0.5d0

c     n_sigma/n_total factor
       fabdown=rhob/rhoval

       delc= xx2*etildeb
       delcb = -(1.d0 + CCC)*fabdown*delc
       erevc = erevc + delcb

       erevsic = erevc
         
       if(ldew) func(n) = func(n) + rhoval*erevsic*fac
       Ec = Ec + rhoval*erevsic*qwght(n)*fac
cderivs
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(1)

       btermn=(1.d0+CCC)*(
     &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
     +      xx2*( (etildea - etildeb)*rhob/(rhoval*rhoval) +
     +      fabup*detiladna + fabdown*detilbdna) )

       drevdna=atermn - btermn
       
       atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*dneggadn(2)

       btermn=(1.d0+CCC)*(
     &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
     +      xx2*( (etildeb-etildea)*rhoa/(rhoval*rhoval) +
     +      fabup*detilbdna+fabdown*detilbdnb) )

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

       finaln= rhoval*drevdna + erevc
       Amat(n,1)=Amat(n,1)+(finaln)*fac

       finaln= rhoval*drevdnb + erevc 
       Amat(n,2)=Amat(n,2)+(finaln)*fac

c     derivs wrt gaa
         apartg=rhoval*drevdgaa
         finalg=apartg

         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalg*fac

c     derivs wrt gbb
         apartg=rhoval*drevdgbb
         finalgbb=apartg

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


          endif



      return
      end

c

      Subroutine xc_cpkzb99_d2()
      call errquit(' not coded ',0,0)
      return
      end




