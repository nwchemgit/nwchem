#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_camrpbe.F
C> The CAM-RPBE exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
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
C> DOI: <a href="http://dx.doi.org/10.1063/1.2409292">
C> 10.1063/1.2409292</a>.
C>
C> [2] O.A. Vydrov, J. Heyd, A. Krukau, G.E. Scuseria, 
C> “Importance of short-range versus long-range Hartree-Fock exchange
C> for the performance of hybrid density functionals”, J. Chem. Phys.
C> <b>125</b>, 074106 (2006), DOI:
C> <a href="http://dx.doi.org/10.1063/1.2244560">10.1063/1.2244560</a>.
C>
C> [3] O.A. Vydrov, G.E. Scuseria, J.P. Perdew, “Tests of
C> functionals for systems with fractional electron number”, 
C> J. Chem. Phys. <b>126</b>, 154109 (2007), DOI:
C> <a href="http://dx.doi.org/10.1063/1.2723119">10.1063/1.2723119</a>.
C>
C> [4] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>77</b>, 3865-3865 (1996), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.77.3865">
C> 10.1103/PhysRevLett.77.3865</a>.
C>
C> [5] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Erratum: Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>78</b>, 1396-1396 (1997), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.78.1396">
C> 10.1103/PhysRevLett.78.1396</a>.
C>
C> [6] J.P. Perdew, Y. Wang, "Accurate and simple density functional
C> for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>33</b>, 8800-8802 (1986), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevB.33.8800">
C> 10.1103/PhysRevB.33.8800</a>.
C>
C> [7] J.P. Perdew, Y. Wang, "Erratum: Accurate and simple density
C> functional for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>40</b>, 3399-3399 (1989), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevB.40.3399">
C> 10.1103/PhysRevB.40.3399</a>.
C>
C> [8] Y. Zhang, W. Yang, "Comment on: Generalized gradient
C> approximation made simple", Phys. Rev. Lett. <b>80</b>, 890-890
C> (1998), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.80.890">
C> 10.1103/PhysRevLett.80.890</a>.
C>
C> [9] B. Hammer, L.B. Hansen, J.K. N&oslash;rskov,
C> "Improved adsorption energetics within density-functional theory
C> using revised Perdew-Burke-Ernzerhof functionals",
C> Phys. Rev. B <b>59</b>, 7413-7421 (1999), DOI:
C> <a href="http://dx.doi.org/10.1103/PhysRevB.59.7413">
C> 10.1103/PhysRevB.59.7413</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_camrpbe(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &                     Amat, Cmat)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_camrpbe_d2(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &                        Amat, Amat2, Cmat, Cmat2)
#else
      Subroutine nwxc_x_camrpbe_d3(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &     Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3)
#endif
c
c$Id: nwxc_x_campbe.F 23711 2013-03-08 22:18:42Z d3y133 $
c
      implicit none
c
#include "nwxc_param.fh"
c      
      double precision param(*)!< [Input] Parameters of functional
                               !< - param(1): \f$ \alpha_{CAM} \f$
                               !< - param(2): \f$ \beta_{CAM} \f$
                               !< - param(3): \f$ \omega_{CAM} \f$
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
                                   !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
      double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
                                    !< and possibly rho
#endif
      double precision Ex
      double precision Atmp, Ctmp, Etmp
c
      double precision pi, um, uk, umk
      double precision C, Cs
      double precision F43, F13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision F73
#endif
#if defined(THIRD_DERIV)
      double precision F10d3, F23
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
c     parameter(um=0.2195149727645171d0)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      parameter (F73=7.d0/3.d0)
#endif
#if defined(THIRD_DERIV)
      parameter (F10d3=10.0d0/3.0d0, F23=2.0d0/3.0d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision rhom23, d2s(3), gpp, d2g(3)
      double precision A2tmp, C2tmp, C3tmp
#endif
#if defined(THIRD_DERIV)
      double precision d3s(4), d3g(4), rhom53, gppp
      double precision A3tmp, C4tmp, C5tmp, C6tmp
#endif
c     double precision gpbe0,gpbe1,gpbe2,gpbe3
      double precision grpbe0,grpbe1,grpbe2,grpbe3
c     double precision grevpbe0,grevpbe1,grevpbe2
c     gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
c     gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
c     gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
c     gpbe3(s)= 24.0d0*umk*um*s*
c    1  (2.0d0*umk*s*s/(1.0d0+umk*s*s)-1.0d0)/(1.0d0+umk*s*s)**3
c     grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
c     grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
c     grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
c    /     (1d0+umkrev*s*s)**2
      grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
      grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
      grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
      grpbe3(s)= -4.0d0*umk*um*s*exp(-umk*s*s)*(3d0-2d0*umk*s*s)
c
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-11-13): Electron density for the XC-third derivatives
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0.0d0
#endif
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
#endif
c           if (lfac) then
               Etmp = rho43
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp =  2d0*rhom23
#endif
c Daniel (7-27-12): XC-third derivative, coming from LDA.
#ifdef THIRD_DERIV
               A3tmp = -4.0d0*rhom53
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
               g=grpbe0(s)
               gp=grpbe1(s)
c           else
c              g=gpbe0(s)
c              gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 2d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
c           else
c              gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            A2tmp = A2tmp
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))*4d0
            C3tmp = rho43*d2g(3)*8d0
#endif
#if defined(THIRD_DERIV)
c Third derivatives
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
c           if(whichf.eq.'revp') then
c              gppp = grevpbe3(s)
c           elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
c           else
c              gppp = gpbe3(s)
c           endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            A3tmp = A3tmp
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp = 
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*8.0d0
            C5tmp =
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))*16.0d0
            C6tmp =
     1             (rho43*d3g(4))*32.0d0
#endif
c
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_T),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp*wght
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp*wght
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp*wght
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
     &                              + C6tmp*wght

#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_T),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_T),ipol,
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-11-13): Derivative of the electron density
#ifdef THIRD_DERIV
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0.0d0
#endif
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
#endif
c           if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp = 2d0*rhom23
#endif
#if defined(THIRD_DERIV)
               A3tmp = -4.0d0*rhom53
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
               g=grpbe0(s)
               gp=grpbe1(s)
c           else
c              g=gpbe0(s)
c              gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
c           else
c              gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp = rho43*d2g(3)*0.5d0
#endif
#if defined(THIRD_DERIV)
c Third derivatives
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
c           if(whichf.eq.'revp') then
c              gppp = grevpbe3(s)
c           elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
c           else
c              gppp = gpbe3(s)
c           endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            A3tmp = A3tmp 
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp =
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*2.0d0
            C5tmp = 
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))
            C6tmp = 
     1             (rho43*d3g(4))*0.5d0
#endif
c
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_A),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp*wght
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp*wght
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp*wght
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
     1                              + C6tmp*wght
#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_A),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*wght
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_A),ipol,
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = F13*rho13*rrho
#endif
c Daniel (4-11-13): Derivative of the electron density
#if defined(THIRD_DERIV)
            rhom53 = F23*rhom23*rrho
#endif
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43*0.5d0
               Atmp = rho13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               A2tmp= 2.d0*rhom23
#endif
#if defined(THIRD_DERIV)
               A3tmp = -4.0d0*rhom53
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
               g=grpbe0(s)
               gp=grpbe1(s)
c           else
c              g=gpbe0(s)
c              gp=gpbe1(s)
c           endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g*0.5d0
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = 0.5d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c           if(whichf.eq.'revp') then
c              gpp=grevpbe2(s)
c           elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
c           else
c              gpp=gpbe2(s)
c           endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c
            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*2d0
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp =  rho43*d2g(3)*0.5d0
#endif
#if defined(THIRD_DERIV)
c Third derivatives
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
c           if(whichf.eq.'revp') then
c              gppp = grevpbe3(s)
c           elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
c           else
c              gppp = gpbe3(s)
c           endif
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            A3tmp = A3tmp
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*4.0d0
            C4tmp = 
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*2.0d0
            C5tmp =
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))
            C6tmp = 
     1             (rho43*d3g(4))*0.5d0
#endif
c
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_B),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*wght
c
            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp*wght
            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp*wght
            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp*wght
            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
     1                              + C6tmp*wght
#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_B),ipol,
     &           Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*wght
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_B),ipol,
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
#include "nwxc_x_camrpbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_camrpbe.F"
#endif
C> @}
