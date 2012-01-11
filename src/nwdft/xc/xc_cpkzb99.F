#include "dft2drv.fh"
!#define PRINTA 1
      Subroutine xc_cpkzb99(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id$
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
      double precision tau(nq,ipol)
 
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
      double precision  ccc
      parameter (ccc = 0.53d0) !cpkzb empirical parameter

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
      double precision etildea,etildeb,gaa,gbb,gab
      double precision fabup,fabdown
      double precision delrho_A(3,2), rho_A(3)
c
      double precision xx1,xx1a,xx1b,pbedown
      double precision tauwplus,taunplus,rhoval2
      double precision dxx1dna,dxx1dnb
      double precision dxx1adna,dxx1bdnb
      double precision dxx1dgaa,dxx1dgbb
      double precision dxx1adgaa,dxx1bdgbb
      double precision drevdgab
      double precision dxx1dta,dxx1dtb
      double precision dxx1adta,dxx1bdtb
      double precision finalna,finalnb
      double precision finalgaa,finalgab
      double precision rhoa2,rhob2
      double precision detiladgaa,detiladgbb
      double precision detilbdgaa,detilbdgbb
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
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
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

         tauN = tau(n,1)
         tauw = 0.125d0*gammaval/rhoval
         xx2 = (tauw/tauN)**2.d0
         en = pbe*(1.d0 + ccc*xx2)
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
     &        dneFSPdn,dneFSPdg, 1, 2, neFSP)
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


        atermn=pbe*ccc*dtwt2dn+(1.d0+ccc*xx2)*decggadn
        btermn=(1.d0+ccc)*(xx2*dpbeupdn + pbeup*dtwt2dn)
        drevdn=atermn - btermn

        atermg=pbe*ccc*dtwt2dg+(1.d0+ccc*xx2)*decggadg
        btermg=(1.d0+ccc)*(xx2*dpbeupdg+pbeup*dtwt2dg)
        drevdg=atermg-btermg
        drevdt=(ccc*pbe-(1.d0+ccc)*pbeup)*xx2*(-2.d0/tauN)
        
       delc = -(1.d0 + ccc)*delc
       erevc = en + delc
         
       if(ldew) func(n) = func(n) + rhoval*erevc*fac
       Ec = Ec + rhoval*erevc*qwght(n)*fac

c     derivs wrt n
         finaln= rhoval*drevdn + erevc 
         Amat(n,1)=Amat(n,1)+(finaln)*fac

c     derivs wrt g
         apartg=rhoval*drevdg
         finalg=apartg
         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2.d0*finalg*fac

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
c
         rhoval = rho(n,1)
         rhoval2 = rhoval*rhoval
c
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
c
         gaa = delrho(n,1,1)*delrho(n,1,1) +
     &         delrho(n,2,1)*delrho(n,2,1) +
     &         delrho(n,3,1)*delrho(n,3,1)
         gbb = delrho(n,1,2)*delrho(n,1,2) +
     &         delrho(n,2,2)*delrho(n,2,2) +
     &         delrho(n,3,2)*delrho(n,3,2)
         gab = delrho(n,1,1)*delrho(n,1,2) +
     &         delrho(n,2,1)*delrho(n,2,2) +
     &         delrho(n,3,1)*delrho(n,3,2)
c
         rhoa=rho(n,2)   
         rhoa2 = rhoa*rhoa
         rhob=rho(n,3)   
         rhob2 = rhob*rhob
c
c        Check for small densities (H atom case as well)
c
         if ((rhoa.lt.tol_rho).or.
     &          (rhob.lt.tol_rho)) goto 20
c
         tauwa = 0.125d0*gaa/rhoa
         tauwb = 0.125d0*gbb/rhob
c
         tauna = tau(n,1)
         taunb = tau(n,2)
c
         tauw = tauwa+tauwb
         taun = tauna+taunb
c
         xx1 = tauw/taun
         xx2 = xx1*xx1
c
         xx1a = tauwa/tauna  
         xx2a = xx1a*xx1a 
c
         xx1b = tauwb/taunb  
         xx2b = xx1b*xx1b 
c
         en = pbe*(1.d0 + ccc*xx2)
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
     &        dneFSPdn,dneFSPdg, 1, 2, neFSP)
         pbeup = neFSP
c
c        functional deriv info below fffffffffffff
         etildea= pbeup
         detiladna = dneFSPdn(1)
         detiladnb = 0d0
         detiladgaa = dneFSPdg(D1_GAA)
         detiladgbb = 0d0
c
c     n_sigma/n_total factor
       fabup=rhoa/rhoval
       delca = -(1.d0 + ccc)*fabup*xx2a*etildea
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
     &        dneFSPdn,dneFSPdg, 1, 2, neFSP)
       pbedown = neFSP
c
c      functional deriv info below fffffffffffff
       etildeb= pbedown
       detilbdna=0d0
       detilbdnb = dneFSPdn(1)
       detilbdgaa=0d0
       detilbdgbb = dneFSPdg(D1_GAA)
c
c     n_sigma/n_total factor
       fabdown=rhob/rhoval
       delcb = -(1.d0 + ccc)*fabdown*xx2b*etildeb
       erevc = erevc + delcb
c
       if(ldew) func(n) = func(n) + rhoval*erevc*fac
       Ec = Ec + rhoval*erevc*qwght(n)*fac
c
c na
       dxx1dna = -0.125d0*gaa/(taun*rhoa2)
       dxx1adna = -0.125d0*gaa/(tauna*rhoa2)
       atermn=pbe*ccc*2.d0*xx1*dxx1dna + (1.d0+ccc*xx2)*dneggadn(1)
       btermn= (1.d0+ccc)*(2.d0*xx1a*dxx1adna*fabup*etildea + 
     &                     xx2a*etildea*fabdown/rhoval +
     &                     xx2a*fabup*detiladna - 
     &                     xx2b*etildeb*fabdown/rhoval)
       drevdna = atermn - btermn
c
c nb
       dxx1dnb = -0.125d0*gbb/(taun*rhob2)
       dxx1bdnb = -0.125d0*gbb/(taunb*rhob2)
       atermn=pbe*ccc*2.d0*xx1*dxx1dnb + (1.d0+ccc*xx2)*dneggadn(2)
       btermn= (1.d0+ccc)*(2.d0*xx1b*dxx1bdnb*fabdown*etildeb +
     &                     xx2b*etildeb*fabup/rhoval +
     &                     xx2b*fabdown*detilbdnb - 
     &                     xx2a*etildea*fabup/rhoval)
       drevdnb = atermn - btermn
c
c gaa
       dxx1dgaa = 0.125d0/(taun*rhoa)
       dxx1adgaa = 0.125d0/(tauna*rhoa)
       atermg=(1.d0+ccc*xx2)*dneggadg(D1_GAA)+ pbe*ccc*2.d0*xx1*dxx1dgaa
       btermg=(1.d0+ccc)*(2.d0*xx1a*dxx1adgaa*fabup*etildea + 
     &    xx2a*fabup*detiladgaa)
       drevdgaa = atermg - btermg
c
c gbb
       dxx1dgbb = 0.125d0/(taun*rhob)
       dxx1bdgbb = 0.125d0/(taunb*rhob)
       atermg=(1.d0+ccc*xx2)*dneggadg(D1_GBB)+ pbe*ccc*2.d0*xx1*dxx1dgbb
       btermg=(1.d0+ccc)*(2.d0*xx1b*dxx1bdgbb*fabdown*etildeb +
     &   xx2b*fabdown*detilbdgbb)
       drevdgbb = atermg - btermg
c
c gab
       atermg=(1.d0+ccc*xx2)*dneggadg(D1_GAB)
       drevdgab = atermg
c
c ta
       dxx1dta=-xx1/taun
       dxx1adta=-xx1a/tauna
       drevdta=pbe*2.d0*ccc*xx1*dxx1dta
     &        -(1.d0+ccc)*2.d0*xx1a*dxx1adta*fabup*etildea
c
c tb
       dxx1dtb=-xx1/taun
       dxx1bdtb=-xx1b/taunb
       drevdtb=pbe*2.d0*ccc*xx1*dxx1dtb
     &        -(1.d0+ccc)*2.d0*xx1b*dxx1bdtb*fabdown*etildeb
c
c derivs wrt na,nb
       finalna= rhoval*drevdna + erevc
       Amat(n,1)=Amat(n,1)+finalna*fac

       finalnb= rhoval*drevdnb + erevc 
       Amat(n,2)=Amat(n,2)+finalnb*fac
c
c     derivs wrt gaa
       finalgaa=rhoval*drevdgaa
       Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalgaa*fac
c
c     derivs wrt gbb
       finalgbb=rhoval*drevdgbb
       Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+ finalgbb*fac
c
c     derivs wrt gab
       finalgab=rhoval*drevdgab
       Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+ finalgab*fac
c
c     derivs wrt ta,tb
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
      call errquit(' cpkzb99: d2 not coded ',0,0)
      return
      end
