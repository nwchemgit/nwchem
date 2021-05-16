#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_vs98.F
C> Implementation of the VS98 exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The Voorhis and Scuseria exchange functional
C>
C> The VS98 functional [1,2] is a meta-GGA. This routine implements
C> the exchange component.
C>
C> ### References ###
C>
C> [1] T van Voorhis, GE Scuseria,
C>     "A novel form for the exchange-correlation energy functional",
C>     J.Chem.Phys. <b>109</b>, 400-410 (1998), DOI:
C>     <a href="https://doi.org/10.1063/1.476577">
C>     10.1063/1.476577</a>.
C>
C> [2] T van Voorhis, GE Scuseria,
C>     Erratum: "A novel form for the exchange-correlation energy
C>     functional",
C>     J.Chem.Phys. <b>129</b>, 219901-219901 (2008), DOI:
C>     <a href="https://doi.org/10.1063/1.3005348">
C>     10.1063/1.3005348</a>.
C>
c   VS98  exchange functional   
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauU - uniform-gas KE density
c                              ijzy - 1  VS98
c                              ijzy - 2  M06-L  
c                              ijzy - 3  M06-HF
c                              ijzy - 4  M06
c                              
c     References: 
c 
c     [a] T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998). 
c     [b] Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006).


#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_vs98_p(param,tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_x_vs98(param,tol_rho, ipol, nq, wght, rho, rgamma,
     &                       tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_vs98_d2(param,tol_rho, ipol, nq, wght, rho,
     &           rgamma, tau, func)
#else
      Subroutine nwxc_x_vs98_d3(param,tol_rho, ipol, nq, wght, rho,
     &           rgamma, tau, func)
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
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of the functional (see Table 1 of [1])
      type(nwad_dble)::r1, r2, r3, r4, r5, r6
#else
      double precision param(*) !< [Input] Parameters of the functional (see Table 1 of [1])
      double precision r1, r2, r3, r4, r5, r6
#endif
#else
      double precision param(*) !< [Input] Parameters of the functional (see Table 1 of [1])
                                !< - param(1): \f$ a \f$
                                !< - param(2): \f$ b \f$
                                !< - param(3): \f$ c \f$
                                !< - param(4): \f$ d \f$
                                !< - param(5): \f$ e \f$
                                !< - param(6): \f$ f \f$
      double precision r1, r2, r3, r4, r5, r6
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
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c      
c     The functional
c
      type(nwad_dble)::func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] The derivative wrt tau
c
c     Local
c
      integer n
      double precision pi
      type(nwad_dble):: rho43, rhoo, rho53, rho83
      type(nwad_dble):: Gamma
c
      type(nwad_dble)::tauN,tauu
      double precision DTol
      double precision Tiny, f13, f43, f53, f83, f113 
      type(nwad_dble)::gx, x, z, kx,xk,zk
      double precision gg
      double precision One, Two, Three, Four, Five, Six, Seven, Eight
      double precision Nine, F10, F11
      double precision cf, Axlsda

c     functional derivatives below FFFFFFFFFFFF

      double precision dxdr, dxdg, dzdr, dzdt, dgdx, dgdz

c     functional derivatives above FFFFFFFFFFFF


cedo       parameter( pi = 3.1415926535897932384626433832795d0 )
         
      parameter (cf = 9.115599720d0, Axlsda = -0.9305257363491d0 )
      parameter (gg  = 0.00186726d0) 
      parameter (f13=1.d0/3.d0,f43=4.0d0/3.0d0,f53=5.0d0/3.0d0)
      parameter (f83=8.d0/3.0d0, F113=11.0d0/3.d0)
      parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &           Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &           Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)
      pi=acos(-1.0d0)

      r1 = param(1)
      r2 = param(2)
      r3 = param(3)
      r4 = param(4)
      r5 = param(5)
      r6 = param(6)
c
c     Parameters for VS98
c
c     r1=  -9.800683d-01
c     r2=  -3.556788d-03
c     r3=   6.250326d-03
c     r4=  -2.354518d-05
c     r5=  -1.282732d-04
c     r6=   3.574822d-04
      
      DTol = tol_rho
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c
         do 10 n = 1, nq
            rhoo = rho(n,R_T)/Two   ! rho_sigma
            if (rhoo.lt.DTol) goto 10
            rho43 = rhoo**f43  
c           rrho = 1d0/rhoo       ! reciprocal of rho
c           rho13 = rho43*rrho
            rho53 = rhoo**f53
            rho83 = rho53*rhoo 
c
            tauN = tau(n,T_T)
c           if(taun.lt.dtol) goto 10
            tauu=tauN
            gamma = rgamma(n,G_TT)
            gamma = gamma/Four
            x = gamma/rho83
c           dxdr = -f83*x*rrho
c           dxdg = One/rho83
            z = tauu/rho53 - cf
c           dzdr = -f53 * tauu/rho83
c           dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
            call nwxc_gvt4_p(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     &                       r1,r2,r3,r4,r5,r6)
#else
            call nwxc_gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     &                     r1,r2,r3,r4,r5,r6)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxc_gvt4_d2(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     &                        r1,r2,r3,r4,r5,r6)
#else
            call nwxc_gvt4_d3(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     &                        r1,r2,r3,r4,r5,r6)
#endif
           
            func(n)=func(n)+ Two*rho43*gx*wght
c
c     functional derivatives 
c
c           Amat(n,D1_RA)  = Amat(n,D1_RA) + (f43*rho13*gx + 
c    &                       rho43*(dgdx*dxdr + dgdz*dzdr))*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + rho43*(dgdx*dxdg)*wght
c           Mmat(n,D1_TA)  = Mmat(n,D1_TA) + rho43*(dgdz*dzdt)*wght

10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======


         do 20 n = 1, nq
           if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
           if (rho(n,R_A).lt.DTol) goto 25
           rhoo = rho(n,R_A)
           rho43 = rhoo**f43
c          rrho = 1.0d0/rhoo       ! reciprocal of rho
c          rho13 = rho43*rrho
           rho53 = rhoo**f53
           rho83 = rho53*rhoo

c
           tauN = tau(n,T_A)
           tauu = tauN*Two

c          if (tauu).lt.DTol) goto 25

           Gamma = rgamma(n,G_AA)

           x = gamma/rho83
c          dxdr = -f83*x*rrho
c          dxdg = One/rho83
           z = tauu/rho53 - cf
c          dzdr = -f53 * tauu/rho83
c          dzdt = One/rho53
           kx = One + gg*x + gg*z
           xk = x/kx
           zk = z/kx
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_gvt4_p(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                      r1,r2,r3,r4,r5,r6)
#else
           call nwxc_gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_gvt4_d2(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                       r1,r2,r3,r4,r5,r6)
#else
           call nwxc_gvt4_d3(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                       r1,r2,r3,r4,r5,r6)
#endif
c
           func(n)=func(n)+ rho43*gx*wght
c
c     functional derivatives
c
c          Amat(n,D1_RA)  = Amat(n,D1_RA) + (f43*rho13*gx + 
c    &                      rho43*(dgdx*dxdr + dgdz*dzdr))*wght
c          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + rho43*(dgdx*dxdg)*wght
c          Mmat(n,D1_TA)  = Mmat(n,D1_TA) + rho43*(dgdz*dzdt)*wght
c
c     Beta               BETA           BETA
c
25         continue
c
c     Beta
c
           if (rho(n,R_B).lt.DTol) goto 20
           rhoo = rho(n,R_B)
           rho43 = rhoo**f43
c          rrho = 1.0d0/rhoo       ! reciprocal of rho
c          rho13 = rho43*rrho
           rho53 = rhoo**f53
           rho83 = rho53*rhoo
c
           tauN = tau(n,T_B)
           tauu = Two*tauN

c          if (tauu.lt.DTol) goto 20
c
           Gamma = rgamma(n,G_BB)
c
           x = gamma/rho83
c          dxdr = -f83*x*rrho
c          dxdg = One/rho83
           z = tauu/rho53 - cf
c          dzdr = -f53 * tauu/rho83
c          dzdt = One/rho53
           kx = One + gg*x + gg*z
           xk = x/kx
           zk = z/kx
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           call nwxc_gvt4_p(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                      r1,r2,r3,r4,r5,r6)
#else
           call nwxc_gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           call nwxc_gvt4_d2(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                       r1,r2,r3,r4,r5,r6)
#else
           call nwxc_gvt4_d3(gx,dgdx,dgdz,xk,zk,kx,gg,gg,
     +                       r1,r2,r3,r4,r5,r6)
#endif

           func(n)=func(n)+ rho43*gx*wght
c
c     functional derivatives
c
c          Amat(n,D1_RB)  = Amat(n,D1_RB) + (f43*rho13*gx +
c    &                      rho43*(dgdx*dxdr + dgdz*dzdr))*wght
c          Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + rho43*(dgdx*dxdg)*wght
c          Mmat(n,D1_TB)  = Mmat(n,D1_TB) + rho43*(dgdz*dzdt)*wght

c
20      continue
      endif
c
      return
      end
C>
C> \brief Evaluate the GVT4 form in the VS98 functional
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_gvt4_p(g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6)
#else
      Subroutine nwxc_gvt4(g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_gvt4_d2(g,dgdx,dgdz,xk,zk,k,c,ct,
     +                        r1,r2,r3,r4,r5,r6)
#else
      Subroutine nwxc_gvt4_d3(g,dgdx,dgdz,xk,zk,k,c,ct,
     +                        r1,r2,r3,r4,r5,r6)
#endif
#include "nwad.fh"
      Implicit none
c
c    some working variables
      type(nwad_dble)::g,xk,zk,k
      type(nwad_dble)::sxk,szk,sk
      double precision dgdx,dgdz,c,ct
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::r1,r2,r3,r4,r5,r6
      type(nwad_dble)::sr1,sr2,sr3,sr4,sr5,sr6
#else
      double precision r1,r2,r3,r4,r5,r6
      double precision sr1,sr2,sr3,sr4,sr5,sr6
#endif
#else
      double precision r1,r2,r3,r4,r5,r6
      double precision sr1,sr2,sr3,sr4,sr5,sr6
#endif
c     double precision sc,sct,sk2
      double precision sc,sct
      double precision One, Two, Three, Four, Six
      Data One/1.0d0/, Two/2.0d0/, Three/3.0d0/, Four/4.0d0/, Six/6.0d0/
C
      sxk = xk
      szk = zk
      sk = k
      sc = c
      sct = ct
      sr1 = r1
      sr2 = r2
      sr3 = r3
      sr4 = r4
      sr5 = r5
      sr6 = r6
c     sk2 = sk*sk
      g =  (sr1 + sr2*sxk + sr3*szk
     $  +sr4*sxk**2 + sr5*szk*sxk + sr6*szk**2)/sk
c     dgdx =   (-sr1*sc
c    $  +sr2*(One-Two*sc*sxk)
c    $  -Two*sr3*szk*sc
c    $  +sr4*(Two*sxk-Three*sxk*sxk*sc)
c    $  +sr5*(szk -Three*szk*sxk*sc)
c    $  -Three*sr6*szk*szk*sc )/sk2
c     dgdz =   (-sr1*sct
c    $  -Two*sr2*sxk*sct
c    $  +sr3*(One-Two*szk*sct)
c    $  -Three*sr4*sxk*sxk*sct
c    $  +sr5*(sxk-Three*sxk*szk*sct)
c    $  +sr6*(Two*szk-Three*szk*szk*sct))/sk2

      return
      end
#if !defined(NWAD_PRINT)
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_vs98.F"
#endif
#if !defined(SECOND_DERIV)
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_vs98.F"
#endif
#if !defined(THIRD_DERIV)
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_vs98.F"
#endif
#undef NWAD_PRINT
C> @}
