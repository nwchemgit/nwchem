#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_pkzb99.F
C> The Perdew, Kurth, Zupan, and Blaha correlation functional
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The meta-GGA correlation functional by Perdew, Kurth, Zupan
C> and Blaha
C>
C> The meta-GGA correlation functional by Perdew, Kurth, Zupan and Blaha
C> that provides a self-interaction free energy [1,2].
C>
C> Due to the form of the meta-GGAs we need to screen on the kinetic
C> energy density to ensure that LDA will be obtained when the kinetic
C> energy density goes to zero [3].
C>
C> ### References ###
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
C> [3] J. Gr&auml;fenstein, D. Izotov, D. Cremer,
C>     "Avoiding singularity problems associated with meta-GGA exchange
C>     and correlation functionals containing the kinetic energy
C>     density", J. Chem. Phys. <b>127</b>, 214103 (2007), DOI:
C>     <a href="https://doi.org/10.1063/1.2800011">
C>     10.1063/1.2800011</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_pkzb99_p(param,tol_rho, ipol, nq, wght, rho,
     &                           rgamma, tau, func)
#else
      Subroutine nwxc_c_pkzb99(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_pkzb99_d2(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_c_pkzb99_d3(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif
c
c$Id$
c

c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).
c
#include "nwad.fh"
c
      Implicit none
c
#include "intf_nwxc_c_Mpbe.fh"
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      type(nwad_dble)::param(*)!< [Input] Parameters of functional
#else
      double precision param(*)!< [Input] Parameters of functional
#endif
#else
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \beta \f$
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      type(nwad_dble)::func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
c     double precision Amat(nq,*) !< [Output] Derivative wrt density
c     double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
      integer n
      type(nwad_dble)::rhoval,rhoa,rhob

c    first sigma term
      type(nwad_dble):: taun
      double precision  ccc
      parameter (ccc = 0.53d0) !cpkzb empirical parameter

c   Second call to the cPBE subroutine

      type(nwad_dble):: neGGA(1), dneGGAdn(2), dneGGAdg(3)
      type(nwad_dble):: rho_t(1,3), rgamma_t(1,3), tau_t(1,2)
      type(nwad_dble):: tauNA,tauNB
c    
      type(nwad_dble)::pbe,tauw,xx2,en,
     ,     tauwa,tauwb,xx2a,xx2b,dtwat2dg,dtwat2dn,
     ,     dtwbt2dg,dtwbt2dn
      type(nwad_dble)::pbeup,dtwt2dn,decggadn,dtwt2dg,
     ,     delc,decggadg,drevdn,drevdg,drevdt,
     ,     dpbeupdn,dpbeupdg,atermn,btermn,atermg,btermg,
     ,     erevc,finaln,apartg,finalg,apartt,finalt
c
      type(nwad_dble):: neFSP(1), dneFSPdn(2), dneFSPdg(3)
c
      double precision drevdna,drevdnb,drevdgaa,drevdgbb,
     A     drevdta,drevdtb,finalgbb
      type(nwad_dble)::delca,delcb,
     A     detiladga,detiladgb,detilbdga,detilbdgb,
     A     detiladna,detiladnb,detilbdna,detilbdnb
      type(nwad_dble)::etildea,etildeb,gaa,gbb,gab
      type(nwad_dble)::fabup,fabdown
c     double precision delrho_A(3,2)
      type(nwad_dble)::rho_A(1,3)
c
      type(nwad_dble)::xx1,xx1a,xx1b,pbedown
      type(nwad_dble)::tauwplus,taunplus,rhoval2
      double precision dxx1dna,dxx1dnb
      double precision dxx1adna,dxx1bdnb
      double precision dxx1dgaa,dxx1dgbb
      double precision dxx1adgaa,dxx1bdgbb
      double precision drevdgab
      double precision dxx1dta,dxx1dtb
      double precision dxx1adta,dxx1bdtb
      double precision finalna,finalnb
      double precision finalgaa,finalgab
      type(nwad_dble)::rhoa2,rhob2
      type(nwad_dble)::detiladgaa,detiladgbb
      type(nwad_dble)::detilbdgaa,detilbdgbb
      type(nwad_dble)::gammaval, rgamma_A(1,3)
c
      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======
         do 12  n = 1, nq
         if (rho(n,R_T).lt.tol_rho) goto 12

         rhoval = rho(n,R_T)

C   set up values to call PBE subroutine
         rho_t(1,1) = rho(n,R_T)
c do delrho
         rgamma_t(1,1) = rgamma(n,G_TT)
         rgamma_t(1,2) = 0.0d0
         rgamma_t(1,3) = 0.0d0
c        delrho_t(1,1) = delrho(n,1,1)
c        delrho_t(2,1) = delrho(n,2,1)
c        delrho_t(3,1) = delrho(n,3,1)
c        gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
         gammaval = rgamma(n,G_TT)
c        gam12=sqrt(gammaval)
c
c     get E_GGA[rho,gamma]
c
         neGGA(1) = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE
         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,neGGA)
         call nwxc_c_Mpbe_p(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                      neGGA)
         call leave_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,neGGA)
#else
         call nwxc_c_Mpbe(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                    neGGA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                       neGGA)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                       neGGA)
#endif
         pbe = neGGA(1)

         tauN = tau(n,T_T)
         if (tauN.lt.tol_rho) then
           tauw = 0.0d0
           xx2 = 0.0d0
         else
           tauw = 0.125d0*gammaval/rhoval
           xx2 = (tauw/tauN)**2.d0
         endif
         en = pbe*(1.d0 + ccc*xx2)
c
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system
c
         if (ipol.eq.1) then
           rho_A(1,1) = (0.5d0)*rho(n,R_T)   ! total   equals (1/2)n_tot
c          rho_A(2) = (0.5d0)*rho(n,R_T)   ! alpha   equals (1/2)n_tot
           rho_A(1,2) = 0.d0               ! beta  equals zero
           rho_A(1,3) = 0.d0               ! beta  equals zero
           rgamma_A(1,1) = 0.25d0*rgamma(n,G_TT)
           rgamma_A(1,2) = 0.0d0
           rgamma_A(1,3) = 0.0d0
         else
           rho_A(1,1) = rho(n,R_A)   ! total   equals (1/2)n_tot
c          rho_A(2) = rho(n,R_B)   ! alpha   equals (1/2)n_tot
           rho_A(1,2) = 0.d0         ! beta  equals zero
           rho_A(1,3) = 0.d0         ! beta  equals zero
           rgamma_A(1,1) = rgamma(n,G_AA)
           rgamma_A(1,2) = 0.0d0
           rgamma_A(1,3) = 0.0d0
         endif
c        delrho_A(1,1) = (0.5d0)*delrho_t(1,1) ! nabla n_up x  
c        delrho_A(2,1) = (0.5d0)*delrho_t(2,1) ! nabla n_up y
c        delrho_A(3,1) = (0.5d0)*delrho_t(3,1) ! nabla n_up z
     
c        delrho_A(1,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(2,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP(1) = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

c
c     get E_GGA[rho_alpha,0,gamma_alpha,0]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_mpbe",func_gga,param(1),2,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
         call nwxc_c_Mpbe_p(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                      neFSP)
         call leave_subr("nwxc_c_mpbe",func_gga,param(1),2,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
#else
         call nwxc_c_Mpbe(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                    neFSP)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#endif
         pbeup = neFSP(1)

c        functional deriv info below fffffffffffff
c        dtwt2dn = -2.d0*xx2/rhoval
c        decggadn= dneGGAdn(1)
c        dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2.0d0)
c        decggadg= dneGGAdg(1)
         delc= xx2*pbeup

C  eps-tilda is eps^FSP
C  functional deriv info below fffffffffffffffff

c        dpbeupdn = 0.5d0*dneFSPdn(1)
c  above note the .5's.  you are taking deriv wrt total density n
c                        not deriv wrt n_up
c        dpbeupdg = 0.25d0*dneFSPdg(1)
c  note .25 above is because you want gamma=deln_tot*deln_tot


c        atermn=pbe*ccc*dtwt2dn+(1.d0+ccc*xx2)*decggadn
c        btermn=(1.d0+ccc)*(xx2*dpbeupdn + pbeup*dtwt2dn)
c        drevdn=atermn - btermn

c        atermg=pbe*ccc*dtwt2dg+(1.d0+ccc*xx2)*decggadg
c        btermg=(1.d0+ccc)*(xx2*dpbeupdg+pbeup*dtwt2dg)
c        drevdg=atermg-btermg
c        drevdt=(ccc*pbe-(1.d0+ccc)*pbeup)*xx2*(-2.d0/tauN)
        
         delc = -(1.d0 + ccc)*delc
         erevc = en + delc
         
         func(n) = func(n) + rhoval*erevc*wght
 
c     derivs wrt n
c        finaln= rhoval*drevdn + erevc 
c        Amat(n,D1_RA)=Amat(n,D1_RA)+(finaln)*wght

c     derivs wrt g
c        apartg=rhoval*drevdg
c        finalg=apartg
c        Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2.d0*finalg*wght

c     derivs wrt t
c        apartt=rhoval*drevdt
c        finalt=apartt
c        Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght

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
         rho_t(1,1) = rho(n,R_A)
         rho_t(1,2) = rho(n,R_B)
c        delrho_t(1,1) = delrho(n,1,1)
c        delrho_t(2,1) = delrho(n,2,1)
c        delrho_t(3,1) = delrho(n,3,1)
c        delrho_t(1,2) = delrho(n,1,2)
c        delrho_t(2,2) = delrho(n,2,2)
c        delrho_t(3,2) = delrho(n,3,2)
         rgamma_t(1,1) = rgamma(n,G_AA)
         rgamma_t(1,2) = rgamma(n,G_AB)
         rgamma_t(1,3) = rgamma(n,G_BB)
         
         neGGA(1) = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE (n,1)
         dneGGAdn(2) = 0.0d0   !Amat in PBE (n,2)
         dneGGAdg(1) = 0.0d0  !Cmat in PBE--aa
         dneGGAdg(2) = 0.0d0  !Cmat in PBE--ab
         dneGGAdg(3) = 0.0d0  !Cmat in PBE--bb
c
c     get E_GGA[rho,gamma]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,neGGA)
         call nwxc_c_Mpbe_p(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                      neGGA)
         call leave_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_t,rgamma_t,tau_t,neGGA)
#else
         call nwxc_c_Mpbe(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                    neGGA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                       neGGA)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     +                       neGGA)
#endif
         pbe = neGGA(1)
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
         if ((rhoa.lt.tol_rho).or.(rhob.lt.tol_rho)) goto 20
c
         tauna = tau(n,T_A)
         taunb = tau(n,T_B)
         taun  = tauna+taunb
c
         if (rhoa.lt.tol_rho) then
           tauwa = 0.0d0
         else
           tauwa = 0.125d0*gaa/rhoa
         endif
         if (rhob.lt.tol_rho) then
           tauwb = 0.0d0
         else
           tauwb = 0.125d0*gbb/rhob
         endif
         tauw  = tauwa+tauwb
c
         if (taun.lt.tol_rho) then
           xx1 = 0.0d0
         else
           xx1 = tauw/taun
         endif
         xx2 = xx1*xx1
c
         if (tauna.lt.tol_rho) then
           xx1a = 0.0d0
         else
           xx1a = tauwa/tauna  
         endif
         xx2a = xx1a*xx1a 
c
         if (taunb.lt.tol_rho) then
           xx1b = 0.0d0
         else
           xx1b = tauwb/taunb  
         endif
         xx2b = xx1b*xx1b 
c
         en = pbe*(1.d0 + ccc*xx2)
c
c     Alpha bit
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system for alpha spin
c     to get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         rho_A(1,1) = rhoa
         rho_A(1,2) = 0.d0               ! beta  equals zero
c        delrho_A(1,1) = delrho_t(1,1) ! nabla n_up x  
c        delrho_A(2,1) = delrho_t(2,1) ! nabla n_up y
c        delrho_A(3,1) = delrho_t(3,1) ! nabla n_up z
c        delrho_A(1,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(2,2) = 0.d0   ! set beta gradient to zero
c        delrho_A(3,2) = 0.d0   ! set beta gradient to zero
         rgamma_A(1,1) = rgamma(n,G_AA)
         rgamma_A(1,2) = 0.0d0
         rgamma_A(1,3) = 0.0d0

         neFSP(1) = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
         call nwxc_c_Mpbe_p(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                      neFSP)
         call leave_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
#else
         call nwxc_c_Mpbe(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                    neFSP)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#endif
         pbeup = neFSP(1)
c
c        functional deriv info below fffffffffffff
         etildea= pbeup
c        detiladna = dneFSPdn(1)
c        detiladnb = 0d0
c        detiladgaa = dneFSPdg(D1_GAA)
c        detiladgbb = 0d0
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
         rho_A(1,1) = 0.0d0         ! alpha equals zero
         rho_A(1,2) = rhob
c        delrho_A(1,1) = delrho_t(1,2) ! nabla n_up x  
c        delrho_A(2,1) = delrho_t(2,2) ! nabla n_up y
c        delrho_A(3,1) = delrho_t(3,2) ! nabla n_up z
c        delrho_A(1,2) = 0.d0     ! set beta gradient to zero
c        delrho_A(2,2) = 0.d0     ! set beta gradient to zero
c        delrho_A(3,2) = 0.d0     ! set beta gradient to zero
         rgamma_A(1,1) = 0.0d0
         rgamma_A(1,2) = 0.0d0
         rgamma_A(1,3) = rgamma(n,G_BB)
       
         neFSP(1) = 0.0d0            !Ec in PBE
         dneFSPdn(1) = 0.0d0      !Amat in PBE
         dneFSPdn(2) = 0.0d0      !Amat in PBE
         dneFSPdg(1) = 0.0d0      !Cmat in PBE
         dneFSPdg(2) = 0.0d0      !Cmat in PBE
         dneFSPdg(3) = 0.0d0      !Cmat in PBE
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call enter_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
         call nwxc_c_Mpbe_p(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                      neFSP)
         call leave_subr("nwxc_c_mpbe",func_gga,param(1),ipol,1,
     +                   rho_A,rgamma_A,tau_t,neFSP)
#else
         call nwxc_c_Mpbe(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                    neFSP)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,2,1,1.0d0,rho_A,rgamma_A,
     +                       neFSP)
#endif
       pbedown = neFSP(1)
c
c      functional deriv info below fffffffffffff
       etildeb= pbedown
c      detilbdna=0d0
c      detilbdnb = dneFSPdn(1)
c      detilbdgaa=0d0
c      detilbdgbb = dneFSPdg(D1_GAA)
c
c     n_sigma/n_total factor
       fabdown=rhob/rhoval
       delcb = -(1.d0 + ccc)*fabdown*xx2b*etildeb
       erevc = erevc + delcb
c
       func(n) = func(n) + rhoval*erevc*wght
c
c na
c      dxx1dna = -0.125d0*gaa/(taun*rhoa2)
c      dxx1adna = -0.125d0*gaa/(tauna*rhoa2)
c      atermn=pbe*ccc*2.d0*xx1*dxx1dna + (1.d0+ccc*xx2)*dneggadn(1)
c      btermn= (1.d0+ccc)*(2.d0*xx1a*dxx1adna*fabup*etildea + 
c    &                     xx2a*etildea*fabdown/rhoval +
c    &                     xx2a*fabup*detiladna - 
c    &                     xx2b*etildeb*fabdown/rhoval)
c      drevdna = atermn - btermn
c
c nb
c      dxx1dnb = -0.125d0*gbb/(taun*rhob2)
c      dxx1bdnb = -0.125d0*gbb/(taunb*rhob2)
c      atermn=pbe*ccc*2.d0*xx1*dxx1dnb + (1.d0+ccc*xx2)*dneggadn(2)
c      btermn= (1.d0+ccc)*(2.d0*xx1b*dxx1bdnb*fabdown*etildeb +
c    &                     xx2b*etildeb*fabup/rhoval +
c    &                     xx2b*fabdown*detilbdnb - 
c    &                     xx2a*etildea*fabup/rhoval)
c      drevdnb = atermn - btermn
c
c gaa
c      dxx1dgaa = 0.125d0/(taun*rhoa)
c      dxx1adgaa = 0.125d0/(tauna*rhoa)
c      atermg=(1.d0+ccc*xx2)*dneggadg(D1_GAA)+ pbe*ccc*2.d0*xx1*dxx1dgaa
c      btermg=(1.d0+ccc)*(2.d0*xx1a*dxx1adgaa*fabup*etildea + 
c    &    xx2a*fabup*detiladgaa)
c      drevdgaa = atermg - btermg
c
c gbb
c      dxx1dgbb = 0.125d0/(taun*rhob)
c      dxx1bdgbb = 0.125d0/(taunb*rhob)
c      atermg=(1.d0+ccc*xx2)*dneggadg(D1_GBB)+ pbe*ccc*2.d0*xx1*dxx1dgbb
c      btermg=(1.d0+ccc)*(2.d0*xx1b*dxx1bdgbb*fabdown*etildeb +
c    &   xx2b*fabdown*detilbdgbb)
c      drevdgbb = atermg - btermg
c
c gab
c      atermg=(1.d0+ccc*xx2)*dneggadg(D1_GAB)
c      drevdgab = atermg
c
c ta
c      dxx1dta=-xx1/taun
c      dxx1adta=-xx1a/tauna
c      drevdta=pbe*2.d0*ccc*xx1*dxx1dta
c    &        -(1.d0+ccc)*2.d0*xx1a*dxx1adta*fabup*etildea
c
c tb
c      dxx1dtb=-xx1/taun
c      dxx1bdtb=-xx1b/taunb
c      drevdtb=pbe*2.d0*ccc*xx1*dxx1dtb
c    &        -(1.d0+ccc)*2.d0*xx1b*dxx1bdtb*fabdown*etildeb
c
c derivs wrt na,nb
c      finalna= rhoval*drevdna + erevc
c      Amat(n,D1_RA)=Amat(n,D1_RA)+finalna*wght

c      finalnb= rhoval*drevdnb + erevc 
c      Amat(n,D1_RB)=Amat(n,D1_RB)+finalnb*wght
c
c     derivs wrt gaa
c      finalgaa=rhoval*drevdgaa
c      Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ finalgaa*wght
c
c     derivs wrt gbb
c      finalgbb=rhoval*drevdgbb
c      Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+ finalgbb*wght
c
c     derivs wrt gab
c      finalgab=rhoval*drevdgab
c      Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+ finalgab*wght
c
c     derivs wrt ta,tb
c      apartt=rhoval*drevdta
c      finalt=apartt
c      Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght

c      apartt=rhoval*drevdtb
c      finalt=apartt
c      Mmat(n,D1_TB)=Mmat(n,D1_TB)+0.5d0*finalt*wght

20     continue

      endif

      return
      end
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_pkzb99.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_pkzb99.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_pkzb99.F"
#endif
#undef NWAD_PRINT
C>
C> @}
