#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_rpbe.F
C> The RPBE exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the RPBE exchange functional
C>
C> The PBE exchange functional [1,2] can be seen as a simplification
C> of the PW91 functional [3,4]. The functional form is a 
C> straightforward GGA. The functional parameters were later revised
C> to achieve better atomization energies, and satisfy the Lieb-Oxford
C> criterion locally [5].
C>
C> [1] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>77</b>, 3865-3865 (1996), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865">
C> 10.1103/PhysRevLett.77.3865</a>.
C>
C> [2] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Erratum: Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>78</b>, 1396-1396 (1997), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.78.1396">
C> 10.1103/PhysRevLett.78.1396</a>.
C>
C> [3] J.P. Perdew, Y. Wang, "Accurate and simple density functional
C> for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>33</b>, 8800-8802 (1986), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.33.8800">
C> 10.1103/PhysRevB.33.8800</a>.
C>
C> [4] J.P. Perdew, Y. Wang, "Erratum: Accurate and simple density
C> functional for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>40</b>, 3399-3399 (1989), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.40.3399">
C> 10.1103/PhysRevB.40.3399</a>.
C>
C> [5] B. Hammer, L.B. Hansen, J.K. N&oslash;rskov,
C> "Improved adsorption energetics within density-functional theory
C> using revised Perdew-Burke-Ernzerhof functionals",
C> Phys. Rev. B <b>59</b>, 7413-7421 (1999), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.59.7413">
C> 10.1103/PhysRevB.59.7413</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_rpbe_p(
     W     tol_rho, ipol, nq, wght, rho, rgamma, func)
#else
      Subroutine nwxc_x_rpbe(
     W     tol_rho, ipol, nq, wght, rho, rgamma, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_rpbe_d2(
     W     tol_rho, ipol, nq, wght, rho, rgamma, func)
#else
      Subroutine nwxc_x_rpbe_d3(
     W     tol_rho, ipol, nq, wght, rho, rgamma, func)
#endif
c
c$Id$
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision amat(nq,*) !< [Output] The derivative wrt rho
c     double precision cmat(nq,*) !< [Output] The derivative wrt rgamma
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
c                                  !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
c     double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
c                                   !< and possibly rho
#endif
c
      double precision pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13, F23
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision F73
#endif
#if defined(THIRD_DERIV)
c     double precision F10d3
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0, F23=2.0d0/3.0d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     parameter (F73=7.d0/3.d0)
#endif
#ifdef THIRD_DERIV
c     parameter (F10d3=10.0d0/3.0d0)
#endif
c
      integer n
      type(nwad_dble)::rrho, rho43, rho13, gamma, gam12, s, g 
      double precision d1s(2), d, gp, d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision rhom23, d2s(3), gpp, d2g(3)
#endif
#if defined(THIRD_DERIV)
c     double precision d3s(4), d3g(4), rhom53, gppp
#endif
c     double precision gpbe0,gpbe1,gpbe2
c     type(nwad_dble)::grpbe0,grpbe1,grpbe2,grpbe3
c     double precision grevpbe0,grevpbe1,grevpbe2
c     gpbe0(s)= uk*(1.0d0 - 1.0d0/(1.0d0+umk*s*s))
c     gpbe1(s)= 2.0d0*um*s/(1.0d0+umk*s*s)**2
c     gpbe2(s)= 2.0d0*um*(1.0d0-4.0d0*umk*s*s/(1.0d0+umk*s*s))/
c    +     (1.0d0+umk*s*s)**2
c     grevpbe0(s)= ukrev*(1.0d0 - 1.0d0/(1.0d0+umkrev*s*s))
c     grevpbe1(s)= 2.0d0*um*s/(1.0d0+umkrev*s*s)**2
c     grevpbe2(s)= 2.0d0*um*(1.0d0-4.0d0*umkrev*s*s/(1.0d0+umkrev*s*s))/
c    /     (1.0d0+umkrev*s*s)**2
c     grpbe0(s)= uk*(1.0d0 - exp(-umk*s*s))
c     grpbe1(s)= 2.0d0*um*s*exp(-umk*s*s)
c     grpbe2(s)= 2.0d0*um*exp(-umk*s*s)*(1.0d0-2.0d0*umk*s*s)
c     grpbe3(s)= -4.0d0*umk*um*s*exp(-umk*s*s)*(3d0-2d0*umk*s*s)
c
      pi = acos(-1.0d0)
      C = -3.0d0/(4.0d0*pi)*(3.0d0*pi*pi)**F13
      Cs = 0.5d0/(3.0d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
c Daniel (9-28-12): There are somewhat mysterious coefficients involved
c in the evaluation of the functional and its derivatives.  We must
c recall that the exchange energy is always written based on the
c spin-scaling relationship for exchange:
c
c Ex[rho] = Ex[rho_a,rho_b] = 0.5*( Ex[2*rho_a] + Ex[2*rho_b] )
c 
c Thus, the electron density is always written:
c rho -> 2*rho_s
c gamma -> 4*gamma_ss
c
c Rationalization for the coefficients is mathematically justified below:
c
c ----------------------------
c Amat       -> 0.5*2 = 1
c Cmat       -> 0.5*4 = 2
c ----------------------------
c Amat2      -> 0.5*2*2 = 2
c Cmat2(rg)  -> 0.5*2*4 = 4
c Cmat2(gg)  -> 0.5*4*4 = 8
c ----------------------------
c Amat3      -> 0.5*2*2*2 = 4
c Cmat3(rrg) -> 0.5*2*2*4 = 8
c Cmat3(rgg) -> 0.5*2*4*4 = 16
c Cmat3(ggg) -> 0.5*4*4*4 = 32
c ----------------------------
c
c If, instead, the author of this code had decided to divide the total
c density (rho(n,1)) by 2 in constructing the density and gamma, those
c coefficients would be unnecessary.
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
            rho43 = C*rho(n,R_T)**F43
            rrho = 1.0d0/rho(n,R_T)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
#if defined(THIRD_DERIV)
c           rhom53 = F23*rhom23*rrho
#endif
c           if (lfac) then
c              func(n) = func(n) + wght*rho43*fac
c              Amat(n,1) = Amat(n,1) + wght*rho13*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + wght*2d0*rhom23*fac
#endif
c           endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
c           gam12 = sqrt(rgamma(n,G_TT))
c           if (gam12.lt.tol_rho**2) goto 10
c

c           s = Cs*gam12/rho43
            s = Cs*Cs*rgamma(n,G_TT)/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/rgamma(n,G_TT)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
c           gp=grpbe1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g
c           Amat(n,D1_RA)  = Amat(n,D1_RA) 
c    +                     + wght*(rho13*g+rho43*d1g(1))
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) 
c    +                     + wght*2.0d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_TT)
c           gpp=grpbe2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +wght*(rhom23*g
c    &           +      2.d0*rho13*d1g(1)
c    &           +      rho43*d2g(1))*2.0d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +wght*(rho13*d1g(2)
c    &           +      rho43*d2g(2))*4.0d0
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           +wght*rho43*d2g(3)*8.0d0
#endif
#if defined(THIRD_DERIV)
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c           d3s(1) = -F10d3*d2s(1)*rrho
c           d3s(2) = 0.5d0*d2s(1)/gamma
c           d3s(3) = -F43*d2s(3)*rrho
c           d3s(4) = -1.5d0*d2s(3)/gamma
c           gppp = grpbe3(s)
c Third derivatives of the enhancement factor
c           d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1) 
c    1             + gppp*d1s(1)*d1s(1)*d1s(1)
c           d3g(2) = gp*d3s(2)
c    1             + gpp*d1s(2)*d2s(1)
c    2             + 2.0d0*gpp*d1s(1)*d2s(2)
c    3             + gppp*d1s(1)*d1s(1)*d1s(2)
c           d3g(3) = gp*d3s(3)
c    1             + gpp*d1s(1)*d2s(3)
c    2             + 2.0d0*gpp*d1s(2)*d2s(2)
c    3             + gppp*d1s(1)*d1s(2)*d1s(2)
c           d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
c    1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    1           + (-rhom53*g
c    2           +  3.0d0*rhom23*d1g(1)
c    3           +  3.0d0*rho13*d2g(1)
c    4           +  rho43*d3g(1))*wght*4.0d0
c           Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
c    1           + (rhom23*d1g(2)
c    2           +  2.0d0*rho13*d2g(2)
c    3           +  rho43*d3g(2))*wght*8.0d0
c           Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
c    1           + (rho13*d2g(3) 
c    2           +  rho43*d3g(3))*wght*16.0d0
c           Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
c    1           + (rho43*d3g(4))*wght*32.0d0
#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c Daniel (9-28-12): There are somewhat mysterious coefficients involved
c in the evaluation of the functional and its derivatives.  We must
c recall that the exchange energy is always written based on the
c spin-scaling relationship for exchange:
c
c Ex[rho] = Ex[rho_a,rho_b] = 0.5*( Ex[2*rho_a] + Ex[2*rho_b] )
c 
c Thus, the electron density is always written:
c rho -> 2*rho_s
c gamma -> 4*gamma_ss
c
c It seems like the derivatives should be correctly balanced by the
c following coefficients:
c
c -----------------------------
c Amat       -> 0.5*2 = 1
c Cmat       -> 0.5*1 = 0.5
c -----------------------------
c Amat2      -> 0.5*2*2 = 2
c Cmat2(rg)  -> 0.5*2*1 = 1
c Cmat2(gg)  -> 0.5*1*1 = 0.5
c -----------------------------
c Amat3      -> 0.5*2*2*2 = 4
c Cmat3(rrg) -> 0.5*2*2*1 = 2
c Cmat3(rgg) -> 0.5*2*1*1 = 1
c Cmat3(ggg) -> 0.5*1*1*1 = 0.5
c -----------------------------
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 20 n = 1, nq
c           if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c     - factor 1/2 on tolerance needed for consistency with closed shell
c       case
c
            if (rho(n,R_A).lt.0.5d0*tol_rho) goto 25
            rho43 = C*(2d0*rho(n,R_A))**F43
            rrho = 0.5d0/rho(n,R_A)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
#if defined(THIRD_DERIV)
c           rhom53 = F23*rhom23*rrho
#endif
c           if (lfac) then
c              Ex = Ex + rho43*qwght(n)*fac*0.5d0
c              if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c              Amat(n,1) = Amat(n,1) + rho13*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
#endif
c           endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
c           gam12 = 2.0d0*sqrt(rgamma(n,G_AA))
c           if (gam12.lt.tol_rho**2) goto 25
c
c           s = Cs*gam12/rho43
            s = 4.0d0*Cs*Cs*rgamma(n,G_AA)/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/rgamma(n,G_AA)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
c           gp=grpbe1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g*0.5d0
c           Amat(n,D1_RA)  = Amat(n,D1_RA) 
c    +                     + wght*(rho13*g+rho43*d1g(1))
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
c    +                     + wght*0.5d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_AA)
c           gpp=grpbe2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +wght*(rhom23*g
c    &           +      2.d0*rho13*d1g(1)
c    &           +      rho43*d2g(1))*2.0d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +wght*(rho13*d1g(2)
c    &           +      rho43*d2g(2))
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           +wght*rho43*d2g(3)*0.5d0
#endif
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c           d3s(1) = -F10d3*d2s(1)*rrho
c           d3s(2) = 0.5d0*d2s(1)/gamma
c           d3s(3) = -F43*d2s(3)*rrho
c           d3s(4) = -1.5d0*d2s(3)/gamma
c           gppp = grpbe3(s)
c Third derivatives of the enhancement factor
c           d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
c    1             + gppp*d1s(1)*d1s(1)*d1s(1)
c           d3g(2) = gp*d3s(2)
c    1             + gpp*d1s(2)*d2s(1)
c    2             + 2.0d0*gpp*d1s(1)*d2s(2)
c    3             + gppp*d1s(1)*d1s(1)*d1s(2)
c           d3g(3) = gp*d3s(3)
c    1             + gpp*d1s(1)*d2s(3)
c    2             + 2.0d0*gpp*d1s(2)*d2s(2)
c    3             + gppp*d1s(1)*d1s(2)*d1s(2)
c           d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
c    1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c    1           + (-rhom53*g
c    2           +  3.0d0*rhom23*d1g(1)
c    3           +  3.0d0*rho13*d2g(1)
c    4           +  rho43*d3g(1))*wght*4.0d0
c           Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
c    1           + (rhom23*d1g(2)
c    2           +  2.0d0*rho13*d2g(2)
c    3           +  rho43*d3g(2))*wght*2.0d0
c           Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA)
c    1           + (rho13*d2g(3)
c    2           +  rho43*d3g(3))*wght
c           Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
c    1           + (rho43*d3g(4))*wght*0.5d0
#endif
c
c     Beta
c     - factor 1/2 on tolerance needed for consistency with closed shell
c       case
c
 25         continue
            if (rho(n,R_B).lt.0.5d0*tol_rho) goto 20
            rho43 = C*(2.0d0*rho(n,R_B))**F43
            rrho = 0.5d0/rho(n,R_B)
            rho13 = F43*rho43*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
#if defined(THIRD_DERIV)
c           rhom53 = F23*rhom23*rrho
#endif
c           if (lfac) then
c              Ex = Ex + rho43*qwght(n)*fac*0.5d0
c              if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c              Amat(n,2) = Amat(n,2) + rho13*fac
#ifdef SECOND_DERIV
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*rhom23*fac
#endif
c           endif
c
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
c           gam12 = 2.0d0*sqrt(rgamma(n,G_BB))
c           if (gam12.lt.tol_rho**2) goto 20
c
c           s = Cs*gam12/rho43
            s = 4.0d0*Cs*Cs*rgamma(n,G_BB)/(rho43*rho43)
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/rgamma(n,G_BB)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
c           gp=grpbe1(s)
c
c           d1g(1) = gp*d1s(1)
c           d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g*0.5d0
c           Amat(n,D1_RB)  = Amat(n,D1_RB)
c    +                     + wght*(rho13*g+rho43*d1g(1))
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
c    +                     + wght*0.5d0*rho43*d1g(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_BB)
c           gpp=grpbe2(s)
c           d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c           d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c           d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &           +wght*(rhom23*g
c    &           +      2.d0*rho13*d1g(1)
c    &           +      rho43*d2g(1))*2.0d0
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &           +wght*(rho13*d1g(2)
c    &           +      rho43*d2g(2))
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &           +wght*rho43*d2g(3)*0.5d0
#endif
#ifdef THIRD_DERIV
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c           d3s(1) = -F10d3*d2s(1)*rrho
c           d3s(2) = 0.5d0*d2s(1)/gamma
c           d3s(3) = -F43*d2s(3)*rrho
c           d3s(4) = -1.5d0*d2s(3)/gamma
c           gppp = grpbe3(s)
c Third derivatives of the enhancement factor
c           d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
c    1             + gppp*d1s(1)*d1s(1)*d1s(1)
c           d3g(2) = gp*d3s(2)
c    1             + gpp*d1s(2)*d2s(1)
c    2             + 2.0d0*gpp*d1s(1)*d2s(2)
c    3             + gppp*d1s(1)*d1s(1)*d1s(2)
c           d3g(3) = gp*d3s(3)
c    1             + gpp*d1s(1)*d2s(3)
c    2             + 2.0d0*gpp*d1s(2)*d2s(2)
c    3             + gppp*d1s(1)*d1s(2)*d1s(2)
c           d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
c    1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
c           Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
c    1           + (-rhom53*g
c    2           +  3.0d0*rhom23*d1g(1)
c    3           +  3.0d0*rho13*d2g(1)
c    4           +  rho43*d3g(1))*wght*4.0d0
c           Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
c    1           + (rhom23*d1g(2)
c    2           +  2.0d0*rho13*d2g(2)
c    3           +  rho43*d3g(2))*wght*2.0d0
c           Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB)
c    1           + (rho13*d2g(3)
c    2           +  rho43*d3g(3))*wght
c           Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
c    1           + (rho43*d3g(4))*wght*0.5d0
#endif
c
 20      continue
      endif
c
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
      function grpbe0(r) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble), intent(in) :: r
        type(nwad_dble)             :: s
        s = uk*(1d0 - exp(-umk*r))
      end function
c
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_rpbe.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_rpbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_rpbe.F"
#endif
#undef NWAD_PRINT
C> @}
