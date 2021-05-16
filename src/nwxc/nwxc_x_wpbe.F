#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_wpbe.F
C> The long range corrected WPBE exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
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
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_wpbe_p(param,
     &    tol_rho, ipol, nq, wght, rho, rgamma, func)
#else
      Subroutine nwxc_x_wpbe(param,
     &    tol_rho, ipol, nq, wght, rho, rgamma, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_wpbe_d2(param,
     &    tol_rho, ipol, nq, wght, rho, rgamma, func)
c    &    Amat, Amat2, Cmat, Cmat2) 
#else
      Subroutine nwxc_x_wpbe_d3(param,
     W     tol_rho, ipol, nq, wght, rho, rgamma, func)
c    &     Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3)
#endif
c
#include "nwad.fh"
c      
      implicit none
c
#include "nwxc_param.fh"
c      
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of functional
      type(nwad_dble)::cam_omega
#else
      double precision param(*) !< [Input] Parameters of functional
      double precision cam_omega
#endif
#else
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ \omega_{CAM} \f$
      double precision cam_omega
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
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c     double precision Amat2(nq,*) !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*) !< [Output] The 2nd derivative wrt rgamma
                                   !< and possibly rho
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
      double precision pi
      double precision C, Cs
      double precision F43, F13
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision F73
#endif
c
c Daniel (4-12-13): XC-3rd derivative information for d3s
#if defined(THIRD_DERIV)
c     double precision F10d3, F23
#endif
c     parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
c     parameter(um=0.2195149727645171d0)
c     parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     parameter (F73=7.d0/3.d0)
#endif
#if defined(THIRD_DERIV)
c     parameter (F10d3=10.0d0/3.0d0, F23=2.d0/3.d0)
#endif
c
      integer n
      type(nwad_dble):: Fxhse,rrho,rho43,rho13,gamma,gam12,s,g,kf
      double precision d1s(2),
     &      d, gp1, d1gs(2), gp2, d1gnu
      double precision nu, d1nu
      double precision Fxhse10,Fxhse01
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision rhom23, gpp1, d2s(3),d2gs(3), gpp2, d2gnu
c     double precision Fxhse20,Fxhse02,d2nu,gpp11,d2gsnu,Fxhse11
c     double precision d1nud1rhonu,d1nud1rhos,d1sd1rhos,d1sd1rhonu
#endif
#if defined(THIRD_DERIV)
c     double precision rhom53, gppp1, gppp21, gppp12, gppp2
c     double precision d3nu, d3s(4), d3gs(4)
c     double precision Fxhse30, Fxhse21, Fxhse12, Fxhse03 
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
c Daniel (4-12-13): For information on the constants multiplied into 
c the partial derivatives, see nwxc_x_pbe.F
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
c Daniel (4-12-13): Electron density for the XC-third derivatives
#ifdef THIRD_DERIV
c           rhom53 = -F23*rhom23*rrho
#endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
            gam12 = sqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 10
c

            s = Cs*gam12/rho43
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c           nu = cam_omega/kf
c           d1nu = -F13*nu*rrho
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
#if defined(THIRD_DERIV)
            call nwxc_x_HSE08_d3(cam_omega,ipol,rho(n,R_T),s,fxhse)
c    &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11,
c    &           fxhse30,fxhse21,fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call nwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_T),s,fxhse)
c    &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
#if defined(NWAD_PRINT)
            call nwxc_x_HSE08_p(cam_omega,ipol,rho(n,R_T),s,fxhse)
#else
            call nwxc_x_HSE08(cam_omega,ipol,rho(n,R_T),s,fxhse)
c    &           fxhse10,fxhse01)
#endif
#endif
            g = fxhse - 1.0d0
c           gp1 = fxhse10
c           gp2 = fxhse01
c
c           d1gs(1) = gp1*d1s(1)
c           d1gs(2) = gp1*d1s(2)
c           d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght
c           Amat(n,D1_RA)= Amat(n,D1_RA)
c    +         +(rho13*g+rho43*(d1gs(1)+d1gnu))*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 2d0*rho43*d1gs(2)*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c
c           d2nu = -F43*rrho*d1nu
c
c           gpp1=fxhse20
c           gpp2=fxhse02
c           gpp11=fxhse11
c
c           d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
c           d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
c    &                gpp11*d1nu*d1s(2)
c           d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
c           d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
c           d2gsnu = 2d0*gpp11*d1nu*d1s(1)

c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*(d1gs(1)+d1gnu)
c    &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +(rho13*d1gs(2)
c    &           + rho43*d2gs(2))*wght*4d0
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           + rho43*d2gs(3)*wght*8d0
#endif
#if defined(THIRD_DERIV)
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c     d3s(1) = -F10d3*d2s(1)*rrho
c     d3s(2) = 0.5d0*d2s(1)/gamma
c     d3s(3) = -F43*d2s(3)*rrho
c     d3s(4) = -1.5d0*d2s(3)/gamma
c
c     d3nu = -F73*d2nu*rrho
c
c     d3gs(1) = d3nu*Fxhse01 
c    1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 )  
c    2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
c    3        + d1nu*( d2nu*Fxhse02 
c    4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
c    5               + Fxhse11*d2s(1) 
c    6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
c    7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
c    8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
c    9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
c    A        + Fxhse10*d3s(1)

c     d3gs(2) = d2nu*d1s(2)*Fxhse11 
c    1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
c    2        + d1s(2)*Fxhse20*d2s(1) 
c    3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
c    4               + d1s(2)*d1s(1)*Fxhse21 ) 
c    5        + Fxhse10*d3s(2) 
c    6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
c    7                 + d1s(2)*d1s(1)*Fxhse30 )

c     d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
c    1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 ) 
c    2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

c     d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2) 
c    1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
c     Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
c    1              + ( rhom53*g
c    2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
c    3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
c    4                + rho43*d3gs(1) )*wght*4.0d0
c     Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
c    1              + ( rhom23*d1gs(2)
c    2                + 2.0d0*rho13*d2gs(2)
c    3                + rho43*d3gs(2) )*wght*8.0d0
c     Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) 
c    1              + ( rho43*d3gs(3) + rho13*d2gs(3) )*wght*16.0d0
c     Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
c    1              + rho43*d3gs(4)*wght*32.0d0
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
c Daniel (4-12-13): Derivative of the electron density
#if defined(THIRD_DERIV)
c           rhom53 = -F23*rhom23*rrho
#endif
c
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            gam12 = 2d0*sqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 25
c
c
            s = Cs*gam12/rho43
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c
c           nu = cam_omega/kf
c           d1nu = -F13*nu*rrho
#if defined(THIRD_DERIV)
            call nwxc_x_HSE08_d3(cam_omega,ipol,rho(n,R_A),s,fxhse)
c    &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11,
c    &           fxhse30,fxhse21,fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call nwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_A),s,fxhse)
c    &                     fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
#if defined(NWAD_PRINT)
            call nwxc_x_HSE08_p(cam_omega,ipol,rho(n,R_A),s,fxhse)
#else
            call nwxc_x_HSE08(cam_omega,ipol,rho(n,R_A),s,fxhse)
c    &                     fxhse10,fxhse01)
#endif
#endif
            g = fxhse - 1.0d0
c           gp1 = fxhse10
c           gp2 = fxhse01 
c
c           d1gs(1) = gp1*d1s(1)
c           d1gs(2) = gp1*d1s(2)
c           d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght*0.5d0
c           Amat(n,D1_RA) = Amat(n,D1_RA)
c    +          + (rho13*g+rho43*(d1gs(1)+d1gnu))*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 0.5d0*rho43*d1gs(2)*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c
c           d2nu = -F43*d1nu*rrho
c     
c           gpp1 = fxhse20
c           gpp2 = fxhse02
c           gpp11 = fxhse11
c
c           d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
c           d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
c    &                gpp11*d1nu*d1s(2)
c           d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
c           d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
c           d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*(d1gs(1)+d1gnu)
c    &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           +(rho13*d1gs(2)
c    &           + rho43*d2gs(2))*wght
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &           + rho43*d2gs(3)*wght*0.5d0
#endif
#if defined(THIRD_DERIV)
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c     d3s(1) = -F10d3*d2s(1)*rrho
c     d3s(2) = 0.5d0*d2s(1)/gamma
c     d3s(3) = -F43*d2s(3)*rrho
c     d3s(4) = -1.5d0*d2s(3)/gamma
c
c     d3nu = -F73*d2nu*rrho
c
c     d3gs(1) = d3nu*Fxhse01 
c    1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 ) 
c    2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
c    3        + d1nu*( d2nu*Fxhse02 
c    4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
c    5               + Fxhse11*d2s(1) 
c    6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
c    7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
c    8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
c    9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
c    A        + Fxhse10*d3s(1)

c     d3gs(2) = d2nu*d1s(2)*Fxhse11 
c    1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
c    2        + d1s(2)*Fxhse20*d2s(1) 
c    3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
c    4               + d1s(2)*d1s(1)*Fxhse21 ) 
c    5        + Fxhse10*d3s(2) 
c    6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
c    7                 + d1s(2)*d1s(1)*Fxhse30 )

c     d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
c    1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 )   
c    2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

c     d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2)
c    1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
c     Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) 
c    1              + ( rhom53*g
c    2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
c    3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
c    4                + rho43*d3gs(1) )*wght*4.0d0
c     Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
c    1              + ( rhom23*d1gs(2)
c    2                + 2.0d0*rho13*d2gs(2)
c    3                + rho43*d3gs(2) )*wght*2.0d0
c     Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) 
c    1              + ( rho43*d3gs(3)
c    2                + rho13*d2gs(3) )*wght
c     Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) 
c    1              + (rho43*d3gs(4))*wght*0.50d0
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rhom23 = F13*rho13*rrho
#endif
#ifdef THIRD_DERIV
c           rhom53 = -F23*rhom23*rrho
#endif
c
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            gam12 = 2d0*sqrt(gamma)
            if (gam12.lt.tol_rho**2) goto 20
c
c
            s = Cs*gam12/rho43
c           d1s(1) = -F43*s*rrho
c           d1s(2) = 0.5d0*s/gamma
c
c
c           nu = cam_omega/kf
c           d1nu = -F13*nu*rrho
c
#if defined(THIRD_DERIV)
            call nwxc_x_HSE08_d3(cam_omega,ipol,rho(n,R_B),s,fxhse)
c    &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11,
c    &           fxhse30,fxhse21,fxhse12,fxhse03)
#elif defined(SECOND_DERIV)
            call nwxc_x_HSE08_d2(cam_omega,ipol,rho(n,R_B),s,fxhse)
c    &           fxhse10,fxhse01,fxhse20,fxhse02,fxhse11)
#else
#if defined(NWAD_PRINT)
            call nwxc_x_HSE08_p(cam_omega,ipol,rho(n,R_B),s,fxhse)
#else
            call nwxc_x_HSE08(cam_omega,ipol,rho(n,R_B),s,fxhse)
c    &           fxhse10,fxhse01)
#endif
#endif
            g = fxhse - 1.0d0
c           gp1= fxhse10
c           gp2 = fxhse01
c
c           d1gs(1) = gp1*d1s(1)
c           d1gs(2) = gp1*d1s(2)
c           d1gnu = gp2*d1nu
            func(n) = func(n) + rho43*g*wght*0.5d0
c           Amat(n,D1_RB) = Amat(n,D1_RB)
c    +          + (rho13*g+rho43*(d1gs(1)+d1gnu))*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + 0.5d0*rho43*d1gs(2)*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c
c           d2s(1) = -F73*d1s(1)*rrho
c           d2s(2) = -F43*d1s(2)*rrho
c           d2s(3) = -0.5d0*d1s(2)/gamma
c
c           d2nu = -F43*d1nu*rrho
c
c           gpp1 = fxhse20
c           gpp2 = fxhse02
c           gpp11 = fxhse11
c
c           d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
c           d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
c    &                gpp11*d1nu*d1s(2)
c           d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
c           d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
c           d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &           +(rhom23*g
c    &           + 2.d0*rho13*(d1gs(1)+d1gnu)
c    &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*wght*2d0
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &           +(rho13*d1gs(2)
c    &           + rho43*d2gs(2))*wght
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &           + rho43*d2gs(3)*wght*0.5d0
#endif
#if defined(THIRD_DERIV)
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c     d3s(1) = -F10d3*d2s(1)*rrho
c     d3s(2) = 0.5d0*d2s(1)/gamma
c     d3s(3) = -F43*d2s(3)*rrho
c     d3s(4) = -1.5d0*d2s(3)/gamma
c
c     d3nu = -F73*d2nu*rrho
c
c     d3gs(1) = d3nu*Fxhse01 
c    1        + 2.0d0*d2nu*( d1nu*Fxhse02 + d1s(1)*Fxhse11 ) 
c    2        + 2.0d0*( d1nu*Fxhse11 + d1s(1)*Fxhse20 )*d2s(1) 
c    3        + d1nu*( d2nu*Fxhse02 
c    4               + d1nu*( d1nu*Fxhse03 + d1s(1)*Fxhse12 ) 
c    5               + Fxhse11*d2s(1) 
c    6               + d1s(1)*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) ) 
c    7        + d1s(1)*( d2nu*Fxhse11 + Fxhse20*d2s(1) 
c    8                 + d1nu*( d1nu*Fxhse12 + d1s(1)*Fxhse21 ) 
c    9                 + d1s(1)*( d1nu*Fxhse21 + d1s(1)*Fxhse30 ) ) 
c    A        + Fxhse10*d3s(1)

c     d3gs(2) = d2nu*d1s(2)*Fxhse11 
c    1        + d2s(2)*( d1nu*Fxhse11 + d1s(1)*Fxhse20 ) 
c    2        + d1s(2)*Fxhse20*d2s(1) 
c    3        + d1nu*( Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 
c    4               + d1s(2)*d1s(1)*Fxhse21 ) 
c    5        + Fxhse10*d3s(2) 
c    6        + d1s(1)*( d2s(2)*Fxhse20 + d1nu*d1s(2)*Fxhse21 
c    7                 + d1s(2)*d1s(1)*Fxhse30 )

c     d3gs(3) = Fxhse10*d3s(3) + 2.0d0*d1s(2)*d2s(2)*Fxhse20 
c    1        + d1nu*( d2s(3)*Fxhse11 + d1s(2)*d1s(2)*Fxhse21 ) 
c    2        + d1s(1)*( d2s(3)*Fxhse20 + d1s(2)*d1s(2)*Fxhse30 )

c     d3gs(4) = Fxhse30*d1s(2)*d1s(2)*d1s(2)
c    1        + 3.0d0*Fxhse20*d1s(2)*d2s(3) + Fxhse10*d3s(4)
c
c     Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) 
c    1              + ( rhom53*g
c    2                + 3.0d0*( rhom23*( d1gs(1) + d1gnu ) )
c    3                + 3.0d0*( rho13*( d2gs(1) + d2gnu + d2gsnu ) )
c    4                + rho43*d3gs(1) )*wght*4.0d0
c     Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
c    1              + ( rhom23*d1gs(2)
c    2                + 2.0d0*rho13*d2gs(2)
c    3                + rho43*d3gs(2) )*wght*2.0d0
c     Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) 
c    1              + ( rho43*d3gs(3) 
c    2                + rho13*d2gs(3) )*wght
c     Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) 
c    1              + rho43*d3gs(4)*wght*0.50d0
#endif
c
 20      continue
      endif
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_wpbe.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_wpbe.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_wpbe.F"
#endif
#undef NWAD_PRINT
c $Id$
C> @}
