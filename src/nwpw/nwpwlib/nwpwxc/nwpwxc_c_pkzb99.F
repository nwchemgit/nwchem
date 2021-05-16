#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_pkzb99.F
C> The Perdew, Kurth, Zupan, and Blaha correlation functional
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The meta-GGA correlation functional by Perdew, Kurth, Zupan
C> and Blaha
C>
C> The meta-GGA correlation functional by Perdew, Kurth, Zupan and Blaha
C> that provides a self-interaction free energy [1,2].
C>
C> [1] J.P. Perdew, S. Kurth, A. Zupan, P. Blaha,
C>     <i>"Accurate density functional with correct formal properties:
C>     A Step Beyond the Generalized Gradient Approximation"</i>,
C>     Phys. Rev. Lett. <b>82</b>, 2544-2547 (1999), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevLett.82.2544">
C>     10.1103/PhysRevLett.82.2544</a>.
C>
C> [2] J.P. Perdew, S. Kurth, A. Zupan, P. Blaha,
C>     <i>"Erratum: Accurate density functional with correct formal
C>     properties: A Step Beyond the Generalized Gradient
C>     Approximation"</i>,
C>     Phys. Rev. Lett. <b>82</b>, 5179-5179 (1999), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevLett.82.5179">
C>     10.1103/PhysRevLett.82.5179</a>.
C>
      Subroutine nwpwxc_c_pkzb99(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func, Amat, Cmat, Mmat)
c
c$Id$
c

c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).
      Implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      double precision rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      double precision func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*) !< [Output] Derivative wrt density
      double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
      integer n
      double precision rhoval,rhoa,rhob

c    first sigma term
      double precision  taun
      double precision  ccc
      parameter (ccc = 0.53d0) !cpkzb empirical parameter

c   Second call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision rho_t(3), rgamma_t(3)
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
c     double precision delrho_A(3,2)
      double precision rho_A(3)
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
      double precision gammaval, rgamma_A(3)
c
      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======
         do 12  n = 1, nq
         if (rho(n,R_T).lt.tol_rho) goto 12

         rhoval = rho(n,R_T)

C   set up values to call PBE subroutine
         rho_t(1) = rho(n,R_T)
c do delrho
         rgamma_t(1) = rgamma(n,G_TT)
         rgamma_t(2) = 0.0d0
         rgamma_t(3) = 0.0d0
c        delrho_t(1,1) = delrho(n,1,1)
c        delrho_t(2,1) = delrho(n,2,1)
c        delrho_t(3,1) = delrho(n,3,1)
c        gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
         gammaval = rgamma(n,G_TT)
         gam12=dsqrt(gammaval)
c
c     get E_GGA[rho,gamma]
c
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE
         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE

         call nwpwxc_c_Mpbe(tol_rho, ipol,
     &        rho_t, rgamma_t, neGGA,
     &        dneGGAdn,dneGGAdg) 
         pbe = neGGA

         tauN = tau(n,T_T)
         tauw = 0.125d0*gammaval/rhoval
         xx2 = (tauw/tauN)**2.d0
         en = pbe*(1.d0 + ccc*xx2)
c
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system
c
         if (ipol.eq.1) then
           rho_A(1) = (0.5d0)*rho(n,R_T)   ! total   equals (1/2)n_tot
c          rho_A(2) = (0.5d0)*rho(n,R_T)   ! alpha   equals (1/2)n_tot
           rho_A(2) = 0.d0               ! beta  equals zero
           rho_A(3) = 0.d0               ! beta  equals zero
           rgamma_A(1) = 0.25d0*rgamma(n,G_TT)
           rgamma_A(2) = 0.0d0
           rgamma_A(3) = 0.0d0
         else
           rho_A(1) = rho(n,R_A)   ! total   equals (1/2)n_tot
c          rho_A(2) = rho(n,R_B)   ! alpha   equals (1/2)n_tot
           rho_A(2) = 0.d0         ! beta  equals zero
           rho_A(3) = 0.d0         ! beta  equals zero
           rgamma_A(1) = rgamma(n,G_AA)
           rgamma_A(2) = 0.0d0
           rgamma_A(3) = 0.0d0
         endif
c        delrho_A(1,1) = (0.5d0)*delrho_t(1,1) ! nabla n_up x  
c        delrho_A(2,1) = (0.5d0)*delrho_t(2,1) ! nabla n_up y
c        delrho_A(3,1) = (0.5d0)*delrho_t(3,1) ! nabla n_up z
     
c        delrho_A(1,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(2,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

c
c     get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         call nwpwxc_c_Mpbe(tol_rho, 2, rho_A, rgamma_A,
     &        neFSP, dneFSPdn,dneFSPdg)
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
         
         func(n) = func(n) + rhoval*erevc*wght
 
c     derivs wrt n
         finaln= rhoval*drevdn + erevc 
         Amat(n,D1_RA)=Amat(n,D1_RA)+(finaln)*wght

c     derivs wrt g
         apartg=rhoval*drevdg
         finalg=apartg
         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2.d0*finalg*wght

c     derivs wrt t
         apartt=rhoval*drevdt
         finalt=apartt
         Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght

12     continue
c
c     open-shell
c
       else   !ipol=2 and do alpha beta cases
         do 20 n = 1, nq
c            
         if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
c
         rhoval = rho(n,R_A)+rho(n,R_B)
         rhoval2 = rhoval*rhoval
c
c        rho_t(1) = rho(n,R_A)+rho(n,R_B)
         rho_t(1) = rho(n,R_A)
         rho_t(2) = rho(n,R_B)
c        delrho_t(1,1) = delrho(n,1,1)
c        delrho_t(2,1) = delrho(n,2,1)
c        delrho_t(3,1) = delrho(n,3,1)
c        delrho_t(1,2) = delrho(n,1,2)
c        delrho_t(2,2) = delrho(n,2,2)
c        delrho_t(3,2) = delrho(n,3,2)
         rgamma_t(1) = rgamma(n,G_AA)
         rgamma_t(2) = rgamma(n,G_AB)
         rgamma_t(3) = rgamma(n,G_BB)
         
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE (n,1)
         dneGGAdn(2) = 0.0d0   !Amat in PBE (n,2)
         dneGGAdg(1) = 0.0d0  !Cmat in PBE--aa
         dneGGAdg(2) = 0.0d0  !Cmat in PBE--ab
         dneGGAdg(3) = 0.0d0  !Cmat in PBE--bb
c
c     get E_GGA[rho,gamma]
c
         call nwpwxc_c_Mpbe(tol_rho, ipol,
     &        rho_t, rgamma_t, neGGA,
     &        dneGGAdn,dneGGAdg) 
         pbe = neGGA
c
c        epGGA = (epsilon_c^GGA)  =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe
c
c        gaa = delrho(n,1,1)*delrho(n,1,1) +
c    &         delrho(n,2,1)*delrho(n,2,1) +
c    &         delrho(n,3,1)*delrho(n,3,1)
c        gbb = delrho(n,1,2)*delrho(n,1,2) +
c    &         delrho(n,2,2)*delrho(n,2,2) +
c    &         delrho(n,3,2)*delrho(n,3,2)
c        gab = delrho(n,1,1)*delrho(n,1,2) +
c    &         delrho(n,2,1)*delrho(n,2,2) +
c    &         delrho(n,3,1)*delrho(n,3,2)
         gaa = rgamma(n,G_AA)
         gab = rgamma(n,G_AB)
         gbb = rgamma(n,G_BB)
c
         rhoa=rho(n,R_A)   
         rhoa2 = rhoa*rhoa
         rhob=rho(n,R_B)   
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
         tauna = tau(n,T_A)
         taunb = tau(n,T_B)
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
         rho_A(2) = 0.d0               ! beta  equals zero
c        delrho_A(1,1) = delrho_t(1,1) ! nabla n_up x  
c        delrho_A(2,1) = delrho_t(2,1) ! nabla n_up y
c        delrho_A(3,1) = delrho_t(3,1) ! nabla n_up z
c        delrho_A(1,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(2,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(3,2) = 0.d0   ! set beta gradient to zero
         rgamma_A(1) = rgamma(n,G_AA)
         rgamma_A(2) = 0.0d0
         rgamma_A(3) = 0.0d0

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
         call nwpwxc_c_Mpbe(tol_rho, 2, rho_A, rgamma_A,
     &        neFSP, dneFSPdn,dneFSPdg)
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
       rho_A(2) = 0.d0          ! beta  equals zero
c      delrho_A(1,1) = delrho_t(1,2) ! nabla n_up x  
c      delrho_A(2,1) = delrho_t(2,2) ! nabla n_up y
c      delrho_A(3,1) = delrho_t(3,2) ! nabla n_up z
c      delrho_A(1,2) = 0.d0     ! set beta gradient to zero
c      delrho_A(2,2) = 0.d0     ! set beta gradient to zero
c      delrho_A(3,2) = 0.d0     ! set beta gradient to zero
       rgamma_A(1) = rgamma(n,G_BB)
       rgamma_A(2) = 0.0d0
       rgamma_A(3) = 0.0d0
       
       neFSP = 0.0d0            !Ec in PBE
       dneFSPdn(1) = 0.0d0      !Amat in PBE
       dneFSPdn(2) = 0.0d0      !Amat in PBE
       dneFSPdg(1) = 0.0d0      !Cmat in PBE
       dneFSPdg(2) = 0.0d0      !Cmat in PBE
       dneFSPdg(3) = 0.0d0      !Cmat in PBE
c
       call nwpwxc_c_Mpbe(tol_rho, 2, rho_A, rgamma_A,
     &        neFSP, dneFSPdn,dneFSPdg)
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
       func(n) = func(n) + rhoval*erevc*wght
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
       Amat(n,D1_RA)=Amat(n,D1_RA)+finalna*wght

       finalnb= rhoval*drevdnb + erevc 
       Amat(n,D1_RB)=Amat(n,D1_RB)+finalnb*wght
c
c     derivs wrt gaa
       finalgaa=rhoval*drevdgaa
       Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalgaa*wght
c
c     derivs wrt gbb
       finalgbb=rhoval*drevdgbb
       Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+ finalgbb*wght
c
c     derivs wrt gab
       finalgab=rhoval*drevdgab
       Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+ finalgab*wght
c
c     derivs wrt ta,tb
       apartt=rhoval*drevdta
       finalt=apartt
       Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght

       apartt=rhoval*drevdtb
       finalt=apartt
       Mmat(n,D1_TB)=Mmat(n,D1_TB)+0.5d0*finalt*wght

20     continue

      endif

      return
      end

c
      Subroutine nwpwxc_c_pkzb99_d2()
      call errquit(' cpkzb99: d2 not coded ',0,0)
      return
      end
C>
C> @}
