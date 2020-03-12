#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_wpbe.F
C> The long range corrected WPBE exchange functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the long range corrected WPBE exchange functional
C>
C> This functional is a long range corrected version [1,2] of the PBE
C> GGA functional [3-7]. It is based on an exchange hole satisfying
C> several constraints of the exact functional.
C>
C> ### References ###
C>
C> [1] E. Weintraub, T.M. Henderson, G.E. Scuseria,
C> "Long-range-corrected hybrids based on a new model exchange hole",
C> J. Chem. Theory. Comput. <b>5</b>, 754 (2009), DOI:
C> <a href="https://doi.org/10.1021/ct800530u">
C> 10.1021/ct800530u</a>.
C>
C> [2] T.M. Henderson, B.G. Janesko, G.E. Scuseria,
C> "Generalized gradient approximation model exchange holes for 
C> range-separated hybrids", J. Chem. Phys. <b>128</b>, 194105 (2008),
C> DOI: <a href="https://doi.org/10.1063/1.2921797">
C> 10.1063/1.2921797</a>.
C>
C> [3] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>77</b>, 3865-3865 (1996), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865">
C> 10.1103/PhysRevLett.77.3865</a>.
C>
C> [4] J.P. Perdew, K. Burke, M. Ernzerhof, 
C> "Erratum: Generalized Gradient Approximation Made Simple"
C> Phys. Rev. Lett. <b>78</b>, 1396-1396 (1997), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.78.1396">
C> 10.1103/PhysRevLett.78.1396</a>.
C>
C> [5] J.P. Perdew, Y. Wang, "Accurate and simple density functional
C> for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>33</b>, 8800-8802 (1986), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.33.8800">
C> 10.1103/PhysRevB.33.8800</a>.
C>
C> [6] J.P. Perdew, Y. Wang, "Erratum: Accurate and simple density
C> functional for the electronic exchange energy: Generalized gradient
C> approximation"
C> Phys. Rev. B <b>40</b>, 3399-3399 (1989), DOI:
C> <a href="https://doi.org/10.1103/PhysRevB.40.3399">
C> 10.1103/PhysRevB.40.3399</a>.
C>
C> [7] Y. Zhang, W. Yang, "Comment on: Generalized gradient
C> approximation made simple", Phys. Rev. Lett. <b>80</b>, 890-890
C> (1998), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.80.890">
C> 10.1103/PhysRevLett.80.890</a>.
C>
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_wpbe(param,
     &    tol_rho, ipol, nq, wght, rho, rgamma, func, Amat, Cmat) 
#else
      Subroutine nwpwxc_x_wpbe_d2(param,
     &    tol_rho, ipol, nq, wght, rho, rgamma, func,
     &    Amat, Amat2, Cmat, Cmat2) 
#endif
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ \omega_{CAM} \f$
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
c
      double precision pi
      double precision C, Cs
      double precision F43, F13
#ifdef SECOND_DERIV
      double precision F73
#endif
c     parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(um=0.2195149727645171d0)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#ifdef SECOND_DERIV
      parameter (F73=7.d0/3.d0)
#endif
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp1, d1gs(2), gp2, d1gnu
      double precision nu, kf, d1nu
      double precision Fxhse,Fxhse10,Fxhse01
      double precision cam_omega
#ifdef SECOND_DERIV
      double precision rhom23, gpp1, d2s(3),d2gs(3), gpp2, d2gnu
      double precision Fxhse20,Fxhse02,d2nu,gpp11,d2gsnu,Fxhse11
      double precision d1nud1rhonu,d1nud1rhos,d1sd1rhos,d1sd1rhonu
#endif
c     double precision gpbe0,gpbe1,gpbe2
c     double precision grpbe0,grpbe1,grpbe2
c     double precision grevpbe0,grevpbe1,grevpbe2
c     gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
c     gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
c     gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
c     grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
c     grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
c     grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
c    /     (1d0+umkrev*s*s)**2
c     grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
c     grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
c     grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
c
c     unpack parameters
c
      cam_omega = param(1)
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
            kf = -1d0*pi*rho13 
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
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
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
#ifdef SECOND_DERIV
            call nwpwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_T),s,fxhse,
     &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
            call nwpwxc_x_HSE08(cam_omega,ipol,rho(n,R_T),s,fxhse,
     &           fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght
            Amat(n,D1_RA)= Amat(n,D1_RA)
     +         +(rho13*g+rho43*(d1gs(1)+d1gnu))*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1gs(2)*wght
#ifdef SECOND_DERIV
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*rrho*d1nu
c
c            d1sd1rhos = -F43*rrho
c            d1sd1rhonu = -F13*2d0*cam_omega/gam12
c            d1nud1rhonu = -F13*rrho
c            d1nud1rhos = -F43*0.5d0*gam12*rrho**2d0/cam_omega
c
            gpp1=fxhse20
            gpp2=fxhse02
            gpp11=fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
c            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1) + 
c     &                gp1*d1sd1rhos*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
c            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu + 
c     &              gp2*d1nud1rhonu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c            d2gsnu = 2d0*gpp11*d1nu*d1s(1) + gp1*d1nud1rhos*d1nu +
c     &               gp2*d1sd1rhonu*d1s(1)

            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*wght*4d0
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2gs(3)*wght*8d0
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
c
c     Alpha
c
            if (rho(n,R_A).lt.0.5d0*tol_rho) goto 25
            rho43 = C*(2d0*rho(n,R_A))**F43
            rrho = 0.5d0/rho(n,R_A)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 25
c
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
#ifdef SECOND_DERIV
            call nwpwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_A),s,fxhse,
     &                     fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
            call nwpwxc_x_HSE08(cam_omega,ipol,rho(n,R_A),s,fxhse,
     &                     fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01 
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght*0.5d0
            Amat(n,D1_RA) = Amat(n,D1_RA)
     +          + (rho13*g+rho43*(d1gs(1)+d1gnu))*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1gs(2)*wght
#ifdef SECOND_DERIV
c
c
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c     
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*wght
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &           + rho43*d2gs(3)*wght*0.5d0
#endif
c
c     Beta
c
 25         continue
            if (rho(n,R_B).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,R_B))**F43
            rrho = 0.5d0/rho(n,R_B)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13
#ifdef SECOND_DERIV
            rhom23 = F13*rho13*rrho
#endif
c
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 20
c
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
#ifdef SECOND_DERIV
            call nwpwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_B),s,fxhse,
     &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
            call nwpwxc_x_HSE08(cam_omega,ipol,rho(n,R_B),s,fxhse,
     &           fxhse10,fxhse01)
#endif
            g = fxhse - 1.0d0
            gp1= fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght*0.5d0
            Amat(n,D1_RB) = Amat(n,D1_RB)
     +          + (rho13*g+rho43*(d1gs(1)+d1gnu))*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1gs(2)*wght
#ifdef SECOND_DERIV
c
c
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*wght
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &           + rho43*d2gs(3)*wght*0.5d0
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
#include "nwpwxc_x_wpbe.F"
#endif
c $Id$
C> @}
