#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_tpss03.F
C> The TPSS correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the TPSS correlation functional
C>
C> Evaluate the TPSS meta-GGA correlation functional [1,2]. 
C>
C> Due to the form of the meta-GGAs we need to screen on the kinetic
C> energy density to ensure that LDA will be obtained when the kinetic
C> energy density goes to zero [3].
C>
C> ### References ###
C>
C> [1] J. Tao, J.P. Perdew, V.N. Staveroverov, G.E. Scuseria,
C> "Climbing the density functional ladder: Nonemperical meta-
C> generalized gradient approximation designed for molecules
C> and solids",
C> Phys. Rev. Lett. <b>91</b>, 146401-146404 (2003), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.91.146401">
C> 10.1103/PhysRevLett.91.146401</a>.
C>
C> [2] J.P. Perdew, J. Tao, V.N. Staveroverov, G.E. Scuseria,
C> "Meta-generalized gradient approximation: Explanation of a
C> realistic nonempirical density functional",
C> J. Chem. Phys. <b>120</b>, 6898-6911 (2004), DOI:
C> <a href="https://doi.org/10.1063/1.1665298">
C> 10.1103/1.1665298</a>.
C>
C> [3] J. Gr&auml;fenstein, D. Izotov, D. Cremer,
C>     "Avoiding singularity problems associated with meta-GGA exchange
C>     and correlation functionals containing the kinetic energy
C>     density", J. Chem. Phys. <b>127</b>, 214103 (2007), DOI:
C>     <a href="https://doi.org/10.1063/1.2800011">
C>     10.1063/1.2800011</a>.
C>
c   
c$Id$
c
c
c    Tao,Perdew, Staroverov, Scuseria exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density


#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_tpss03_p(param,tol_rho, ipol, nq, wght,
     &                           rho, rgamma, tau, func)
#else
      Subroutine nwxc_c_tpss03(param,tol_rho, ipol, nq, wght,
     &                         rho, rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_tpss03_d2(param,tol_rho, ipol, nq, wght,
     &                            rho, rgamma, tau, func)
#else
      Subroutine nwxc_c_tpss03_d3(param,tol_rho, ipol, nq, wght,
     &                            rho, rgamma, tau, func)
#endif
c     References:
c     [a] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         PRL 91, 146401 (2003).
c     [b] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         JCP 120, 6898 (2004).
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
c     Local
c
      integer n
      type(nwad_dble)::rhoval

c    call to the cPBE subroutine

      type(nwad_dble):: neGGA(1)
      double precision  dneGGAdn(2), dneGGAdg(3)
      type(nwad_dble)::rho_t(1,3), rgamma_t(1,3)
      type(nwad_dble):: tauNA,tauNB
      type(nwad_dble):: neFSP(1)
      double precision  dneFSPdn(2), dneFSPdg(3)
      type(nwad_dble)::rgamma_A(1,3), rho_A(1,3)

c     TPSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs

      double precision THRD, F43, F73
      type(nwad_dble)::zeta, ccc
      type(nwad_dble)::tauw
      type(nwad_dble)::xx2,xx3
      type(nwad_dble)::pbe, en
      double precision dd
      type(nwad_dble)::rhoa,rhob
      type(nwad_dble)::pbeup, delc, erevc, erevsic,
     A     delca,delcb,pbedown

c     derivsssssssssssssssssssssssss

      double precision  decggadn, dtwt2dn, dtwt3dn
      double precision  dpbeupdn, drevdn
      double precision  atermn, btermn
      double precision  finaln
      double precision  drevdg, dpbeupdg
      double precision  dtwt2dg, dtwt3dg, decggadg
      double precision  atermg, btermg
      double precision  apartg, cpartg(2), finalg
      double precision  finalgaa,finalgbb,finalgab
      double precision drevdt, apartt, bpartt(2),finalt
      double precision drevdta,drevdtb

      double precision dzetadna, dzetadnb
      double precision dcccdna, dcccdnb
      double precision dcccdgaa,dcccdgbb,dcccdgab
      double precision drevdna, drevdnb

      double precision drevdgaa
      double precision drevdgbb
      double precision drevdgab
      double precision dtwt3dt
      type(nwad_dble)::etildea,etildeb
      double precision 
     N     detiladna,detiladnb,detilbdna,detilbdnb,
     D     detiladgaa,detiladgbb,detilbdgaa,detilbdgbb
      double precision detiladgab,detilbdgab
      double precision pi
      type(nwad_dble)::fabup,fabdown
c     type(nwad_dble)::czeta0
      double precision czeta1
      type(nwad_dble)::gaa,gbb,gab
      type(nwad_dble)::xi2
      type(nwad_dble)::delzeta2
      double precision dencxi2zeta,
     D     ddez2dna,ddez2dnb,dxi2dna,dxi2dnb,
     D     ddencxi2dna,ddencxi2dnb,
     D     ddencxi2dxi2,ddencxi2dzeta
      double precision dxi2dgaa
      double precision dxi2dgbb
      double precision dxi2dgab
c     type(nwad_dble)::denccc,onemzeta,onepzeta
      type(nwad_dble)::onemzeta,onepzeta
      type(nwad_dble)::rhoval2
      double precision rhoval3
      type(nwad_dble)::denxi22
      double precision term1,term2
      double precision ddelzeta2dna,ddelzeta2dnb
      double precision ddelzeta2dgaa,ddelzeta2dgbb,ddelzeta2dgab
      double precision ddencccdna,ddencccdnb
c     type(nwad_dble)::gam12,gammaval,tauN
      type(nwad_dble)::gammaval,tauN
      double precision fac
c
      parameter (dd = 2.8d0)
      parameter (THRD = 1d0/3d0)
      parameter (F43 = 4d0/3d0)
      parameter (F73 = 7d0/3d0)
c
c     czeta0(zeta)=
c    &  0.53d0+0.87d0*zeta**2.0d0+0.5d0*zeta**4.0d0+2.26d0*zeta**6.0d0
c     czeta1(zeta)=
c    &  2.d0*0.87d0*zeta+4d0*0.5d0*zeta**3+6d0*2.26d0*zeta**5
c     denccc(zeta,xi2)= 
c    &  1.d0+0.5d0*xi2*((1.d0+zeta)**(-F43)+(1.d0-zeta)**(-F43))
c
      pi=acos(-1d0)
      fac = wght
c
      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======

         do 22  n = 1, nq
         if (rho(n,R_T).lt.tol_rho) goto 22

         rhoval = rho(n,R_T)

C   set up values to call PBE subroutine
         rho_t(1,1) = rho(n,R_T)
c do delrho
         rgamma_t(1,1) = rgamma(n,G_TT)
         rgamma_t(1,2) = 0.0d0
         rgamma_t(1,3) = 0.0d0
c        gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &        delrho(n,2,1)*delrho(n,2,1) +
c    &        delrho(n,3,1)*delrho(n,3,1)
         gammaval = rgamma(n,G_TT)
c        gam12=sqrt(gammaval)
c
c     get E_GGA[rho,gamma]
c
         neGGA(1) = 0.0d0  !Ec in PBE
c        dneGGAdn(1) = 0.0d0   !Amat in PBE
c        dneGGAdg(1) = 0.0d0  !Cmat in PBE
c        dneGGAdg(2) = 0.0d0  !Cmat in PBE
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_Mpbe_p(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#else
         call nwxc_c_Mpbe(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho,ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#endif
         pbe = neGGA(1)

c
c        epGGA = n*(epsilon_c^GGA) / n =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe

         ccc = 0.53d0           !since zeta=0
         tauw = 0.125d0*gammaval/rhoval

         tauN = tau(n,T_T)
         if (tauN.ge.tol_rho) then

           xx2 = (tauw/tauN)**2
           xx3 = (tauw/tauN)**3

         else

           xx2 = 0.0d0
           xx3 = 0.0d0

         endif
         en = pbe*(1.d0 + ccc*xx2)
c
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system
c

         rho_A(1,1) = (0.5d0)*rho(n,R_T)   ! alpha   equals (1/2)n_tot
         rho_A(1,2) = 0.d0                 ! beta  equals zero
         rgamma_A(1,1) = 0.25d0*rgamma_t(1,1)
         rgamma_A(1,2) = 0.0d0
         rgamma_A(1,3) = 0.0d0

         neFSP(1) = 0.0d0  !Ec in PBE
c        dneFSPdn(1) = 0.0d0   !Amat in PBE
c        dneFSPdn(2) = 0.0d0   !Amat in PBE

c        dneFSPdg(1) = 0.0d0  !Cmat in PBE
c        dneFSPdg(2) = 0.0d0  !Cmat in PBE
c        dneFSPdg(3) = 0.0d0  !Cmat in PBE

c
c     get E_GGA[rho_alpha,0,gamma_alpha,0]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_Mpbe_p(param,tol_rho, 2,1,1.0d0, rho_A, rgamma_A,
     &        neFSP) 
#else
         call nwxc_c_Mpbe(param,tol_rho, 2,1,1.0d0, rho_A, rgamma_A,
     &        neFSP) 
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho, 2,1,1.0d0, rho_A, rgamma_A,
     &        neFSP) 
#else
         call nwxc_c_Mpbe_d3(param,tol_rho, 2,1,1.0d0, rho_A, rgamma_A,
     &        neFSP) 
#endif

         pbeup = neFSP(1)

c        functional deriv info below fffffffffffff

c        dtwt2dn = -2.d0*xx2/rhoval
c        dtwt3dn = -3.d0*xx3/rhoval
c        decggadn= dneGGAdn(1)
c        dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)
c        dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)
c        decggadg= dneGGAdg(1)

       if (pbeup.lt.pbe) then
          delc= xx2*pbe
C  eps-tilda is eps_c
c     erev = egga*(1-xx2)
C  functional deriv info below fffffffffffffffff

c         drevdn= -pbe*dtwt2dn+decggadn*(1.d0 - xx2)
c         drevdg= -pbe*dtwt2dg+decggadg*(1.d0 - xx2)
c         drevdt= 2.d0*pbe*xx2/tauN
       else
          delc= xx2*pbeup

C  eps-tilda is eps^FSP
C  functional deriv info below fffffffffffffffff

c         dpbeupdn = 0.5d0*dneFSPdn(1)
c  above note the .5's.  you are taking deriv wrt total density n
c                        not deriv wrt n_up
c         dpbeupdg = 0.25d0*dneFSPdg(1)
c  note .25 above is because you want gamma=deln_tot*deln_tot
c
c         atermn=pbe*ccc*dtwt2dn+(1.d0+ccc*xx2)*decggadn
c         btermn=(1.d0+ccc)*(xx2*dpbeupdn+pbeup*dtwt2dn)
c         drevdn=atermn-btermn
c
c         atermg=pbe*ccc*dtwt2dg+(1.d0+ccc*xx2)*decggadg
c         btermg=(1.d0+ccc)*(xx2*dpbeupdg+pbeup*dtwt2dg)
c         drevdg=atermg-btermg
c
c         drevdt=(ccc*pbe-(1.d0+ccc)*pbeup)*xx2*(-2.d0/tauN)
       endif

       delc = -(1.d0 + ccc)*delc
       erevc = en + delc
       erevsic = erevc*(1.d0 + dd*erevc*xx3)
         
       func(n) = func(n) + rhoval*erevsic*wght
c      Ec = Ec + rhoval*erevsic*qwght(n)*fac

c     derivs wrt n
c      finaln= rhoval*drevdn + erevc +
c    &      dd*(erevc*erevc*xx3 +
c    +      rhoval*(xx3*2.d0*erevc*drevdn +
c    +      erevc*erevc*dtwt3dn))
c        Amat(n,D1_RA)=Amat(n,D1_RA)+(finaln)*wght

c     derivs wrt g
c        apartg=rhoval*drevdg
c        cpartg(1)=erevc*erevc*dtwt3dg
c        cpartg(2)=xx3*2.d0*erevc*drevdg
c        finalg=apartg+rhoval*dd*( cpartg(1)+cpartg(2) )
c        Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2d0*finalg*wght

c     derivs wrt t
c        apartt=rhoval*drevdt
c        bpartt(1)=-erevc*erevc*xx3*3.d0/tauN
c        bpartt(2)=xx3*2.d0*erevc*drevdt
c        finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))
c        Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght

22     continue

      else !ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
c            
         if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20

         rhoval        = 0.0d0
         rho_t(1,1)    = 0.0d0
         rho_t(1,2)    = 0.0d0
         rgamma_t(1,1) = 0.0d0
         rgamma_t(1,2) = 0.0d0
         rgamma_t(1,3) = 0.0d0

         if (rho(n,R_A).ge.0.5d0*tol_rho) then
           rhoval           = rhoval + rho(n,R_A)
           rho_t(1,R_A)     = rho(n,R_A)
           rgamma_t(1,G_AA) = rgamma(n,G_AA)
         endif
         if (rho(n,R_B).ge.0.5d0*tol_rho) then
           rhoval           = rhoval + rho(n,R_B)
           rho_t(1,R_B)     = rho(n,R_B)
           rgamma_t(1,G_BB) = rgamma(n,G_BB)
           if (rho(n,R_A).ge.0.5d0*tol_rho) then
             rgamma_t(1,G_AB) = rgamma(n,G_AB)
           endif
         endif
c         
         neGGA(1) = 0.0d0  !Ec in PBE
c        dneGGAdn(1) = 0.0d0   !Amat in PBE (n,1)
c        dneGGAdn(2) = 0.0d0   !Amat in PBE (n,2)
c        dneGGAdg(1) = 0.0d0  !Cmat in PBE--aa
c        dneGGAdg(2) = 0.0d0  !Cmat in PBE--ab
c        dneGGAdg(3) = 0.0d0  !Cmat in PBE--bb
c
c     get E_GGA[rho,gamma]
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_Mpbe_p(param,tol_rho, ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#else
         call nwxc_c_Mpbe(param,tol_rho, ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho, ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho, ipol,1,1.0d0,rho_t,rgamma_t,
     &        neGGA)
#endif
         pbe = neGGA(1)
c
c        epGGA = (epsilon_c^GGA)  =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe
c
c        gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &        delrho(n,1,2)*delrho(n,1,2) +
c    &        delrho(n,2,1)*delrho(n,2,1) +
c    &        delrho(n,2,2)*delrho(n,2,2) +
c    &        delrho(n,3,1)*delrho(n,3,1) +
c    &        delrho(n,3,2)*delrho(n,3,2) +
c    &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &        delrho(n,2,1)*delrho(n,2,2) +
c    &        delrho(n,3,1)*delrho(n,3,2))
c        gammaval = rgamma(n,G_AA)+rgamma(n,G_BB)+2.0d0*rgamma(n,G_AB)
         gammaval = rgamma_t(1,G_AA)+rgamma_t(1,G_BB)
     &            + 2.0d0*rgamma_t(1,G_AB)
c        gam12=sqrt(gammaval)
         tauNa = tau(n,T_A)
         tauNb = tau(n,T_B)
         taun  = tauna+taunb
         rhoa  = rho_t(1,R_A)
         rhob  = rho_t(1,R_B)
c
c        Check for small densities (H atom case as well)
c
         if ((rhoa.lt.0.5d0*tol_rho).or.
     &       (rhob.lt.0.5d0*tol_rho)) goto 20
c
c        dcccdna=0.d0
c        dcccdnb=0.d0
c        dcccdgaa=0.d0
c        dcccdgbb=0.d0
c        dcccdgab=0.d0
c
         if(rhoval.gt.tol_rho) then
c
            zeta=(rhoa-rhob)/rhoval
            if (zeta.ge. 1.0d0) zeta =  1.0d0
            if (zeta.le.-1.0d0) zeta = -1.0d0
c
            onepzeta = 1.d0+zeta
            onemzeta = 1.d0-zeta
c
c           gaa = delrho(n,1,1)*delrho(n,1,1) +
c    &            delrho(n,2,1)*delrho(n,2,1) +
c    &            delrho(n,3,1)*delrho(n,3,1) 
c           gbb = delrho(n,1,2)*delrho(n,1,2) +
c    &            delrho(n,2,2)*delrho(n,2,2) +
c    &            delrho(n,3,2)*delrho(n,3,2) 
c           gab = delrho(n,1,1)*delrho(n,1,2) +
c    &            delrho(n,2,1)*delrho(n,2,2) +
c    &            delrho(n,3,1)*delrho(n,3,2)
            gaa = rgamma_t(1,G_AA)
            gbb = rgamma_t(1,G_BB)
            gab = rgamma_t(1,G_AB)
c
            rhoval2 = rhoval**2
c           rhoval3 = rhoval*rhoval*rhoval
c           denxi2 = 2.d0*(3.d0*pi*pi*rhoval)**(1.d0/3.d0)
c           denxi22 = denxi2**2
            denxi22 = 4.d0*(3.d0*pi*pi*rhoval)**(2.d0/3.d0)
            delzeta2 = (gaa*onemzeta**2 + gbb*onepzeta**2
     &                  -2.d0*gab*onemzeta*onepzeta)/rhoval2
c
c           dzetadna =  onemzeta/rhoval
c           dzetadnb = -onepzeta/rhoval
c
c           term1=-2.d0*gaa*onemzeta+2.d0*gbb*onepzeta+4.d0*gab*zeta
c           term1=(term1/rhoval2)*dzetadna
c           term2=-2.d0*delzeta2/rhoval
c           ddelzeta2dna = term1 + term2
c
c           term1=-2.d0*gaa*onemzeta+2.d0*gbb*onepzeta+4.d0*gab*zeta
c           term1=(term1/rhoval2)*dzetadnb
c           term2=-2.d0*delzeta2/rhoval
c           ddelzeta2dnb = term1 + term2
c
c           ddelzeta2dgaa =  onemzeta*onemzeta/rhoval2
c           ddelzeta2dgbb =  onepzeta*onepzeta/rhoval2
c           ddelzeta2dgab =  -2.d0*onepzeta*onemzeta/rhoval2
c
            xi2=delzeta2/denxi22
c           dxi2dna=(ddelzeta2dna -(2.d0/3d0)*delzeta2/rhoval)/denxi22
c           dxi2dnb=(ddelzeta2dnb -(2.d0/3d0)*delzeta2/rhoval)/denxi22
c
c           dxi2dgaa=onemzeta*onemzeta/rhoval2/denxi22
c           dxi2dgbb=onepzeta*onepzeta/rhoval2/denxi22
c           dxi2dgab=-2.d0*onepzeta*onemzeta/rhoval2/denxi22
c
            if (1.0d0-zeta.le.tol_rho) then
              ccc=0.0d0 ! limit established by Maxima
            else if (1.0d0+zeta.le.tol_rho) then
              ccc=0.0d0 ! limit established by Maxima
            else
              ccc=czeta0(zeta)/(denccc(zeta,xi2)**4)
            endif
c           ddencccdna=2.d0*(denccc(zeta,xi2)**3)*
c    &       (dxi2dna*(onepzeta**(-F43) + onemzeta**(-F43)) 
c    &        + xi2*F43*(onemzeta**(-F73) - onepzeta**(-F73))*dzetadna)
c           ddencccdnb=2.d0*(denccc(zeta,xi2)**3)*
c    &       (dxi2dnb*(onepzeta**(-F43) + onemzeta**(-F43)) 
c    &        + xi2*F43*(onemzeta**(-F73) - onepzeta**(-F73))*dzetadnb)
c
c           dcccdna=(czeta1(zeta)*dzetadna/denccc(zeta,xi2)**4) - 
c    &              (czeta0(zeta)*ddencccdna/(denccc(zeta,xi2)**8))
c           dcccdnb=(czeta1(zeta)*dzetadnb/denccc(zeta,xi2)**4) - 
c    &              (czeta0(zeta)*ddencccdnb/(denccc(zeta,xi2)**8))
c
c           dcccdgaa=-(czeta0(zeta)/(denccc(zeta,xi2)**8))*
c    &       2.d0*(denccc(zeta,xi2)**3)*
c    &       dxi2dgaa*(onepzeta**(-F43) + onemzeta**(-F43))
c           dcccdgbb=-(czeta0(zeta)/(denccc(zeta,xi2)**8))*
c    &       2.d0*(denccc(zeta,xi2)**3)*
c    &       dxi2dgbb*(onepzeta**(-F43) + onemzeta**(-F43))
c           dcccdgab=-(czeta0(zeta)/(denccc(zeta,xi2)**8))*
c    &       2.d0*(denccc(zeta,xi2)**3)*
c    &       dxi2dgab*(onepzeta**(-F43) + onemzeta**(-F43))
         endif
c
         tauw = 0.125d0*gammaval/rhoval
c
         if (tauN.ge.tol_rho) then
           xx2 = (tauw/tauN)**2
           xx3 = (tauw/tauN)**3
         else
           xx2 = 0.0d0
           xx3 = 0.0d0
         endif
c        dtwt2dn = -2.d0*xx2/rhoval
c        dtwt3dn = -3.d0*xx3/rhoval
c        dtwt3dt = -3.d0*xx3/taun
c        dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)
c        dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)
c
         en = pbe*(1.d0 + ccc*xx2)
c
c     Alpha bit
c     set up values to call PBE subroutine as 
c     Fully SpinPolarized system for alpha spin
c     to get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         rho_A(1,1) = rhoa
         rho_A(1,2) = 0.d0               ! beta  equals zero
         rgamma_A(1,1) = rgamma_t(1,G_AA) ! nabla n_up x  
         rgamma_A(1,2) = 0.0d0
         rgamma_A(1,3) = 0.0d0

         neFSP(1) = 0.0d0  !Ec in PBE
c        dneFSPdn(1) = 0.0d0   !Amat in PBE
c        dneFSPdn(2) = 0.0d0   !Amat in PBE

c        dneFSPdg(1) = 0.0d0  !Cmat in PBE
c        dneFSPdg(2) = 0.0d0  !Cmat in PBE
c        dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_Mpbe_p(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#else
         call nwxc_c_Mpbe(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#endif
         pbeup = neFSP(1)

c        functional deriv info below fffffffffffff
       if (pbeup.lt.pbe) then
          etildea    = pbe
c         detiladna  = dneggadn(1)
c         detiladnb  = dneggadn(2)
c         detiladgaa = dneggadg(D1_GAA)
c         detiladgbb = dneggadg(D1_GBB)
c         detiladgab = dneggadg(D1_GAB)
       else
          etildea    = pbeup
c         detiladna  = dneFSPdn(1)
c         detiladnb  = 0.d0
c         detiladgaa = dneFSPdg(D1_GAA)
c         detiladgbb = 0.d0
c         detiladgab = 0.d0
       endif

c     n_sigma/n_total factor
       fabup=rhoa/rhoval

       delc= xx2*etildea
       delca = -(1.d0 + ccc)*fabup*delc 
       erevc = en + delca
c
c     Beta bit
c     set up values to call PBE subroutine as 
c     Fully SpinPolarized system for beta spin
c     to get E_GGA[rho_beta,0,gamma_beta,0]
c
         rho_A(1,1) = rhob
         rho_A(1,2) = 0.d0                ! beta  equals zero
         rgamma_A(1,1) = rgamma_t(1,G_BB) ! nabla n_up x  
         rgamma_A(1,2) = 0.0d0
         rgamma_A(1,3) = 0.0d0

         neFSP(1) = 0.0d0  !Ec in PBE
c        dneFSPdn(1) = 0.0d0   !Amat in PBE
c        dneFSPdn(2) = 0.0d0   !Amat in PBE

c        dneFSPdg(1) = 0.0d0  !Cmat in PBE
c        dneFSPdg(2) = 0.0d0  !Cmat in PBE
c        dneFSPdg(3) = 0.0d0  !Cmat in PBE
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_Mpbe_p(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#else
         call nwxc_c_Mpbe(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_Mpbe_d2(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#else
         call nwxc_c_Mpbe_d3(param,tol_rho, 2,1,1.0d0,rho_A, rgamma_A,
     &        neFSP)
#endif
         pbedown = neFSP(1)

c        functional deriv info below fffffffffffff

       if (pbedown.lt.pbe) then
          etildeb=pbe
c         detilbdna  = dneggadn(1)
c         detilbdnb  = dneggadn(2)
c         detilbdgaa = dneggadg(D1_GAA)
c         detilbdgbb = dneggadg(D1_GBB)
c         detilbdgab = dneggadg(D1_GAB)
       else
          etildeb    = pbedown
c         detilbdna  = 0.d0
c         detilbdnb  = dneFSPdn(1)
c         detilbdgaa = 0.d0
c         detilbdgbb = dneFSPdg(D1_GAA)
c         detilbdgab = 0.d0
       endif
c
c     n_sigma/n_total factor
       fabdown=rhob/rhoval

       delc= xx2*etildeb
       delcb = -(1.d0 + ccc)*fabdown*delc
       erevc = erevc + delcb

       erevsic = erevc*(1.d0 + dd*erevc*xx3)
         
       func(n) = func(n) + rhoval*erevsic*wght
c
c na
c      atermn=pbe*ccc*dtwt2dn+(1.d0+ccc*xx2)*dneggadn(1)+
c    &      pbe*xx2*dcccdna
c      btermn=(1.d0+ccc)*(
c    &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
c    +      xx2*( (etildea - etildeb)*rhob/(rhoval*rhoval) +
c    +      fabup*detiladna + fabdown*detilbdna) ) +
c    Z      xx2*(fabup*etildea+fabdown*etildeb)*dcccdna
c      drevdna=atermn-btermn
c       
c nb
c      atermn=pbe*ccc*dtwt2dn+(1.d0+ccc*xx2)*dneggadn(2)+
c    &      pbe*xx2*dcccdnb
c      btermn=(1.d0+ccc)*(
c    &      dtwt2dn*(fabup*etildea+fabdown*etildeb) +
c    &      xx2*( (etildeb-etildea)*rhoa/(rhoval*rhoval) +
c    &      fabup*detilbdna+fabdown*detilbdnb) )+
c    &      xx2*(fabup*etildea+fabdown*etildeb)*dcccdnb
c      drevdnb=atermn-btermn
c
c gaa
c      atermg=pbe*ccc*dtwt2dg+(1.d0+ccc*xx2)*dneggadg(D1_GAA)+
c    &      pbe*xx2*dcccdgaa
c      btermg=(1.d0+ccc)*
c    &      (xx2*(fabup*detiladgaa+fabdown*detilbdgaa)+
c    &      (etildea*fabup+etildeb*fabdown)*dtwt2dg)+
c    &      xx2*(etildea*fabup+etildeb*fabdown)*dcccdgaa
c      drevdgaa=atermg-btermg
c
c gbb
c      atermg=pbe*ccc*dtwt2dg+(1.d0+ccc*xx2)*dneggadg(D1_GBB)+
c    &      pbe*xx2*dcccdgbb
c      btermg=(1.d0+ccc)*
c    &      (xx2*(fabup*detiladgbb+fabdown*detilbdgbb)+
c    &      (etildea*fabup+etildeb*fabdown)*dtwt2dg)+
c    &      xx2*(etildea*fabup+etildeb*fabdown)*dcccdgbb
c      drevdgbb=atermg-btermg
c
c gab
c      atermg=pbe*ccc*2.d0*dtwt2dg+(1.d0+ccc*xx2)*dneggadg(D1_GAB)+
c    &      pbe*xx2*dcccdgab
c      btermg=(1.d0+ccc)*
c    &      (xx2*(fabup*detiladgab+fabdown*detilbdgab)+
c    &      (etildea*fabup+etildeb*fabdown)*2.d0*dtwt2dg)+
c    &      xx2*(etildea*fabup+etildeb*fabdown)*dcccdgab
c      drevdgab=atermg-btermg
c
c ta,tb
c      drevdta=-2.d0*xx2/tauN*
c    *(ccc*pbe-(1.d0+ccc)*(fabup*etildea+fabdown*etildeb))
c      drevdtb=-2d0*xx2/tauN*
c    *(ccc*pbe-(1.d0+ccc)*(fabup*etildea+fabdown*etildeb))
c
c derivs wrt na
c        finaln= rhoval*drevdna + erevc +
c    &      dd*(erevc*erevc*xx3 + rhoval*(xx3*2.d0*erevc*drevdna 
c    &     + erevc*erevc*dtwt3dn))
c        Amat(n,D1_RA)=Amat(n,D1_RA)+finaln*wght
c
c derivs wrt nb
c        finaln= rhoval*drevdnb + erevc +
c    &      dd*(erevc*erevc*xx3 + rhoval*(xx3*2.d0*erevc*drevdnb 
c    &     + erevc*erevc*dtwt3dn))
c        Amat(n,D1_RB)=Amat(n,D1_RB)+finaln*wght
c
c derivs wrt gaa
c        apartg=rhoval*drevdgaa
c        cpartg(1)=erevc*erevc*dtwt3dg
c        cpartg(2)=xx3*2.d0*erevc*drevdgaa
c        finalgaa=apartg+rhoval*dd*(cpartg(1)+cpartg(2))
c        Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+finalgaa*wght
c
c derivs wrt gbb
c        apartg=rhoval*drevdgbb
c        cpartg(1)=erevc*erevc*dtwt3dg
c        cpartg(2)=xx3*2.d0*erevc*drevdgbb
c        finalgbb=apartg+rhoval*dd*(cpartg(1)+cpartg(2))
c        Cmat(n,D1_GBB)=Cmat(n,D1_GBB)+finalgbb*wght
c
c derivs wrt gab
c        apartg=rhoval*drevdgab
c        cpartg(1)=erevc*erevc*2.d0*dtwt3dg
c        cpartg(2)=xx3*2.d0*erevc*drevdgab
c        finalgab=apartg+rhoval*dd*(cpartg(1)+cpartg(2))
c        Cmat(n,D1_GAB)=Cmat(n,D1_GAB)+finalgab*wght
c
c derivs wrt ta
c        apartt=rhoval*drevdta
c        bpartt(1)=erevc*erevc*dtwt3dt
c        bpartt(2)=xx3*2.d0*erevc*drevdta
c        finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))
c        Mmat(n,D1_TA)=Mmat(n,D1_TA)+0.5d0*finalt*wght
c
c derivs wrt tb
c        apartt=rhoval*drevdtb
c        bpartt(1)=erevc*erevc*dtwt3dt
c        bpartt(2)=xx3*2.d0*erevc*drevdtb
c        finalt=apartt+dd*rhoval*(bpartt(1)+bpartt(2))
c        Mmat(n,D1_TB)=Mmat(n,D1_TB)+0.5d0*finalt*wght
c
20     continue

      endif  !end ipol=2 case

      return
c
      contains
c
c     The combination of statement functions and derived types with
c     overloaded operators is not properly supported in the Fortran 
c     standard (apparently). Therefore the statement functions from
c     the original subroutine had to be transformed into contained 
c     functions.
c
c     WARNING: This code is EXTREMELY EVIL! Although there appears to be
c     only one function there are actually three with the same name,
c     each one returning results of a different data type. The trick is
c     that depending on the data type the the subroutine that contains
c     these functions changes its name and hence these different
c     functions of the same name do not lead to conflicts. The 
c     alternative would have been to add a forest of conditional 
c     compilation constructs (#ifdef's) to change the function names
c     in the declarations and all places where they are used. That 
c     would have been extremely ugly, so we are between a rock and a
c     hard place on this one.
c
      function czeta0(zeta) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble), intent(in) :: zeta
        type(nwad_dble)             :: s
        s = 0.53d0 + 0.87d0*zeta**2
     &    + 0.5d0*zeta**4 + 2.26d0*zeta**6
      end function
c
      function denccc(zeta,xi2) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble), intent(in) :: zeta
        type(nwad_dble), intent(in) :: xi2
        type(nwad_dble)             :: s
        s = 1.d0 + 0.5d0*xi2*((1.d0+zeta)**(-F43) + (1.d0-zeta)**(-F43))
      end function
c
      end
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_tpss03.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_tpss03.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_tpss03.F"
#endif
#undef NWAD_PRINT
C> @}
