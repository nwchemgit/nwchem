#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_rpbe.F
C> The RPBE exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_rpbe(
     W     tol_rho, ipol, nq, wght, rho, rgamma, func,
     &                     Amat, Cmat)
#else
      Subroutine nwpwxc_x_rpbe_d2(
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
      double precision amat(nq,*) !< [Output] The derivative wrt rho
      double precision cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
      double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
                                   !< and possibly rho
#endif
c
      double precision pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
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
c     double precision gpbe0,gpbe1,gpbe2
      double precision grpbe0,grpbe1,grpbe2
c     double precision grevpbe0,grevpbe1,grevpbe2
c     gpbe0(s)= uk*(1.0d0 - 1.0d0/(1.0d0+umk*s*s))
c     gpbe1(s)= 2.0d0*um*s/(1.0d0+umk*s*s)**2
c     gpbe2(s)= 2.0d0*um*(1.0d0-4.0d0*umk*s*s/(1.0d0+umk*s*s))/
c    +     (1.0d0+umk*s*s)**2
c     grevpbe0(s)= ukrev*(1.0d0 - 1.0d0/(1.0d0+umkrev*s*s))
c     grevpbe1(s)= 2.0d0*um*s/(1.0d0+umkrev*s*s)**2
c     grevpbe2(s)= 2.0d0*um*(1.0d0-4.0d0*umkrev*s*s/(1.0d0+umkrev*s*s))/
c    /     (1.0d0+umkrev*s*s)**2
      grpbe0(s)= uk*(1.0d0 - exp(-umk*s*s))
      grpbe1(s)= 2.0d0*um*s*exp(-umk*s*s)
      grpbe2(s)= 2.0d0*um*exp(-umk*s*s)*(1.0d0-2.0d0*umk*s*s)
      
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
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
            rho43 = C*rho(n,R_T)**F43
            rrho = 1.0d0/rho(n,R_T)
            rho13 = F43*rho43*rrho
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
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
            gam12 = dsqrt(rgamma(n,G_TT))
            if (gam12.lt.tol_rho**2) goto 10
c

            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/rgamma(n,G_TT)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
            gp=grpbe1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g
            Amat(n,D1_RA)  = Amat(n,D1_RA) 
     +                     + wght*(rho13*g+rho43*d1g(1))
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) 
     +                     + wght*2.0d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_TT)
            gpp=grpbe2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +wght*(rhom23*g
     &           +      2.d0*rho13*d1g(1)
     &           +      rho43*d2g(1))*2.0d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +wght*(rho13*d1g(2)
     &           +      rho43*d2g(2))*4.0d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           +wght*rho43*d2g(3)*8.0d0
#endif
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
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
            gam12 = 2.0d0*dsqrt(rgamma(n,G_AA))
            if (gam12.lt.tol_rho**2) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/rgamma(n,G_AA)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
            gp=grpbe1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g*0.5d0
            Amat(n,D1_RA)  = Amat(n,D1_RA) 
     +                     + wght*(rho13*g+rho43*d1g(1))
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     +                     + wght*0.5d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_AA)
            gpp=grpbe2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +wght*(rhom23*g
     &           +      2.d0*rho13*d1g(1)
     &           +      rho43*d2g(1))*2.0d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +wght*(rho13*d1g(2)
     &           +      rho43*d2g(2))
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           +wght*rho43*d2g(3)*0.5d0
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
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
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
            gam12 = 2.0d0*dsqrt(rgamma(n,G_BB))
            if (gam12.lt.tol_rho**2) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/rgamma(n,G_BB)
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            g=grpbe0(s)
            gp=grpbe1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            func(n) = func(n) + wght*rho43*g*0.5d0
            Amat(n,D1_RB)  = Amat(n,D1_RB)
     +                     + wght*(rho13*g+rho43*d1g(1))
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
     +                     + wght*0.5d0*rho43*d1g(2)
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/rgamma(n,G_BB)
            gpp=grpbe2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +wght*(rhom23*g
     &           +      2.d0*rho13*d1g(1)
     &           +      rho43*d2g(1))*2.0d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +wght*(rho13*d1g(2)
     &           +      rho43*d2g(2))
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           +wght*rho43*d2g(3)*0.5d0
#endif
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
#include "nwpwxc_x_rpbe.F"
#endif
C> @}
