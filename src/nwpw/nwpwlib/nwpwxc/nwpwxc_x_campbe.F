#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_campbe.F
C> The CAM-PBE exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Coulomb attenuated PBE exchange functional
C>
C> Evaluates the Coulomb attenuated PBE exchange functional [1-8].
C>
C> ### References ###
C>
C> [1] O.A. Vydrov, G.E. Scuseria, "Assessment of a long range corrected
C> hybrid functional”, J. Chem. Phys. <b>125</b>, 234109 (2006),
C> DOI: <a href="https://doi.org/10.1063/1.2409292">
C> 10.1063/1.2409292</a>.
C>
C> [2] O.A. Vydrov, J. Heyd, A. Krukau, G.E. Scuseria, 
C> “Importance of short-range versus long-range Hartree-Fock exchange
C> for the performance of hybrid density functionals”, J. Chem. Phys.
C> <b>125</b>, 074106 (2006), DOI:
C> <a href="https://doi.org/10.1063/1.2244560">10.1063/1.2244560</a>.
C>
C> [3] O.A. Vydrov, G.E. Scuseria, J.P. Perdew, “Tests of
C> functionals for systems with fractional electron number”, 
C> J. Chem. Phys. <b>126</b>, 154109 (2007), DOI:
C> <a href="https://doi.org/10.1063/1.2723119">10.1063/1.2723119</a>.
C>
C> [4] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>77</b>, 3865-3865 (1996), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865">
C> 10.1103/PhysRevLett.77.3865</a>.
C>
C> [5] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Erratum: Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>78</b>, 1396-1396 (1997), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.78.1396">
C> 10.1103/PhysRevLett.78.1396</a>.
C>
C> [6] J.P. Perdew, Y. Wang, "Accurate and simple density functional
C> for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>33</b>, 8800-8802 (1986), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.33.8800">
C> 10.1103/PhysRevB.33.8800</a>.
C>
C> [7] J.P. Perdew, Y. Wang, "Erratum: Accurate and simple density
C> functional for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>40</b>, 3399-3399 (1989), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.40.3399">
C> 10.1103/PhysRevB.40.3399</a>.
C>
C> [8] Y. Zhang, W. Yang, "Comment on: Generalized gradient
C> approximation made simple", Phys. Rev. Lett. <b>80</b>, 890-890
C> (1998), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.80.890">
C> 10.1103/PhysRevLett.80.890</a>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_campbe(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &                     Amat, Cmat)
#else
      Subroutine nwpwxc_x_campbe_d2(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &                        Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \kappa \f$
                               !< - param(2): \f$ \mu \f$
                               !< - param(3): \f$ \alpha_{CAM} \f$
                               !< - param(4): \f$ \beta_{CAM} \f$
                               !< - param(5): \f$ \omega_{CAM} \f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision func(nq)   !< [Output] The value of the functional
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
                                   !< and possibly rho
#endif
      double precision Ex
      double precision Atmp, Ctmp, Etmp
      double precision A2tmp, C2tmp, C3tmp
c
      double precision pi, um, uk, umk
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
c     parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
c     parameter(um=0.2195149727645171d0)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
      parameter (F73=7.d0/3.d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
#ifdef SECOND_DERIV
      double precision rhom23, d2s(3), gpp, d2g(3)
#endif
      double precision gpbe0,gpbe1,gpbe2
c     double precision grpbe0,grpbe1,grpbe2
c     double precision grevpbe0,grevpbe1,grevpbe2
      gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
      gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
      gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
c     grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
c     grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
c     grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
c    /     (1d0+umkrev*s*s)**2
c     grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
c     grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
c     grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
c
      uk = param(1)
      um = param(2)
      umk = um/uk
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
            rho43 = C*rho(n,R_T)**F43
            rrho = 1d0/rho(n,R_T)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43
               Atmp = rho13
#ifdef SECOND_DERIV
               A2tmp =  2d0*rhom23
#endif
c           endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
            gam12 = dsqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 10
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
c           if(whichf.eq.'revp') then
c              g=grevpbe0(s)
c              gp=grevpbe1(s)
c           elseif(whichf.eq.'rpbe') then
c              g=grpbe0(s)
c              gp=grpbe1(s)
c           else
               g=gpbe0(s)
               gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 2d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
c              gpp=grpbe2(s)
c           else
               gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            A2tmp = A2tmp
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*4d0
            C3tmp = rho43*d2g(3)*8d0
c
            call nwpwxc_x_att_d2(param(3),tol_rho,rho(n,R_T),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
#else
            call nwpwxc_x_att(param(3),tol_rho,rho(n,R_T),ipol,
     &           Etmp,Atmp,Ctmp)
#endif
            func(n) = func(n) + Etmp*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp*wght
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 20 n = 1, nq
c
c     Alpha
c
            if (rho(n,R_A).lt.0.5d0*tol_rho) goto 25
            rho43 = C*(2d0*rho(n,R_A))**F43
            rrho = 0.5d0/rho(n,R_A)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#ifdef SECOND_DERIV
               A2tmp = 2d0*rhom23
#endif
c           endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
c           if(whichf.eq.'revp') then
c              g=grevpbe0(s)
c              gp=grevpbe1(s)
c           elseif(whichf.eq.'rpbe') then
c              g=grpbe0(s)
c              gp=grpbe1(s)
c           else
               g=gpbe0(s)
               gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
c              gpp=grpbe2(s)
c           else
               gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp = rho43*d2g(3)*0.5d0
c
            call nwpwxc_x_att_d2(param(3),tol_rho,rho(n,R_A),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
#else
            call nwpwxc_x_att(param(3),tol_rho,rho(n,R_A),ipol,
     &           Etmp,Atmp,Ctmp)
#endif
            func(n) = func(n) + Etmp*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp*wght
c
c     Beta
c
 25         continue
            if (rho(n,R_B).lt.0.5d0*tol_rho) goto 20
            rho43 = C*(2d0*rho(n,R_B))**F43
            rrho = 0.5d0/rho(n,R_B)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#ifdef SECOND_DERIV
               A2tmp= 2.d0*rhom23
#endif
c           endif
c
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
c           if(whichf.eq.'revp') then
c              g=grevpbe0(s)
c              gp=grevpbe1(s)
c           elseif(whichf.eq.'rpbe') then
c              g=grpbe0(s)
c              gp=grpbe1(s)
c           else
               g=gpbe0(s)
               gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
c              gpp=grpbe2(s)
c           else
               gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp =  rho43*d2g(3)*0.5d0
c
            call nwpwxc_x_att_d2(param(3),tol_rho,rho(n,R_B),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*wght
#else
            call nwpwxc_x_att(param(3),tol_rho,rho(n,R_B),ipol,
     &           Etmp,Atmp,Ctmp)
#endif
            func(n) = func(n) + Etmp*wght
            Amat(n,D1_RB) = Amat(n,D1_RB) + Atmp*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp*wght
c
 20      continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_x_campbe.F"
#endif
C> @}
