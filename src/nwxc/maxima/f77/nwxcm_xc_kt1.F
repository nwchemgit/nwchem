C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_xc_kt1.F 
C> The nwxcm_xc_kt1 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="https://doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_xc_kt1(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Cmat) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            fnc(iq) = fnc(iq)-1.2d-2*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**1.3333333333333333d+0+1.0d-1
            t4 = 1/t3
            fnc(iq) = (-6.0d-3*gammabb*t4-6.0d-3*gammaaa*t2)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.0d-3*gammabb*rhob**3.333333333333333d-1*w
     1         ght/t3**2+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t4*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            fnc(iq) = fnc(iq)-6.0d-3*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            fnc(iq) = fnc(iq)-6.0d-3*gammabb*t2*wght
            Amat(iq,D1_RB) = 8.0d-3*gammabb*rhob**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t2*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="https://doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_xc_kt1_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma 
                                    !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhoa**6.666666666666666d-1
            fnc(iq) = fnc(iq)-1.2d-2*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t3*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = -2.1333333333333332d-2*gammaaa*t5*wght/
     1         t1**3+2.6666666666666666d-3*gammaaa*t4*wght/t5+Amat2(iq,D
     2         2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**1.3333333333333333d+0+1.0d-1
            t4 = 1/t3
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t1**2
            t7 = rhob**3.333333333333333d-1
            t8 = 1/t3**2
            t9 = rhoa**6.666666666666666d-1
            t10 = rhob**6.666666666666666d-1
            fnc(iq) = (-6.0d-3*gammabb*t4-6.0d-3*gammaaa*t2)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t5*t6*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.0d-3*gammabb*t7*t8*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t4*wght
            Amat2(iq,D2_RA_RA) = -2.1333333333333332d-2*gammaaa*t9*wght/
     1         t1**3+2.6666666666666666d-3*gammaaa*t6*wght/t9+Amat2(iq,D
     2         2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 2.6666666666666666d-3*gammabb*t8*wght/t
     1         10-2.1333333333333332d-2*gammabb*t10*wght/t3**3+Amat2(iq,
     2         D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t5*t6*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 8.0d-3*t7*t8*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhoa**6.666666666666666d-1
            fnc(iq) = fnc(iq)-6.0d-3*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t3*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Amat2(iq,D2_RA_RA) = -2.1333333333333332d-2*gammaaa*t5*wght/
     1         t1**3+2.6666666666666666d-3*gammaaa*t4*wght/t5+Amat2(iq,D
     2         2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhob**6.666666666666666d-1
            fnc(iq) = fnc(iq)-6.0d-3*gammabb*t2*wght
            Amat(iq,D1_RB) = 8.0d-3*gammabb*t3*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t2*wght
            Amat2(iq,D2_RB_RB) = -2.1333333333333332d-2*gammabb*t5*wght/
     1         t1**3+2.6666666666666666d-3*gammabb*t4*wght/t5+Amat2(iq,D
     2         2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="https://doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_xc_kt1_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision rgamma(nq,*) !< [Input] The norm of the density 
                                    !< gradients 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma 
                                    !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho 
      double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma 
                                    !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhoa**6.666666666666666d-1
            t6 = 1/t1**3
            t7 = 1/t5
            fnc(iq) = fnc(iq)-1.2d-2*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t3*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 2.6666666666666666d-3*gammaaa*t4*t7*wgh
     1         t-2.1333333333333332d-2*gammaaa*t5*t6*wght+Amat2(iq,D2_RA
     2         _RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = -2.1333333333333332d-2*gammaaa*t6*wg
     1         ht/t3-1.7777777777777776d-3*gammaaa*t4*wght/rhoa**1.66666
     2         66666666669d+0+8.533333333333333d-2*gammaaa*rhoa*wght/t1*
     3         *4+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 2.6666666666666666d-3*t4*t7*wght-2.
     1         1333333333333332d-2*t5*t6*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**1.3333333333333333d+0+1.0d-1
            t4 = 1/t3
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t1**2
            t7 = rhob**3.333333333333333d-1
            t8 = 1/t3**2
            t9 = rhoa**6.666666666666666d-1
            t10 = 1/t1**3
            t11 = 1/t9
            t12 = rhob**6.666666666666666d-1
            t13 = 1/t3**3
            t14 = 1/t12
            fnc(iq) = (-6.0d-3*gammabb*t4-6.0d-3*gammaaa*t2)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t5*t6*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.0d-3*gammabb*t7*t8*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t4*wght
            Amat2(iq,D2_RA_RA) = -2.1333333333333332d-2*gammaaa*t10*t9*w
     1         ght+2.6666666666666666d-3*gammaaa*t11*t6*wght+Amat2(iq,D2
     2         _RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 2.6666666666666666d-3*gammabb*t14*t8*wg
     1         ht-2.1333333333333332d-2*gammabb*t12*t13*wght+Amat2(iq,D2
     2         _RB_RB)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t5*t6*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 8.0d-3*t7*t8*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = -1.7777777777777776d-3*gammaaa*t6*wg
     1         ht/rhoa**1.6666666666666669d+0-2.1333333333333332d-2*gamm
     2         aaa*t10*wght/t5+8.533333333333333d-2*gammaaa*rhoa*wght/t1
     3         **4+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = -1.7777777777777776d-3*gammabb*t8*wg
     1         ht/rhob**1.6666666666666669d+0-2.1333333333333332d-2*gamm
     2         abb*t13*wght/t7+8.533333333333333d-2*gammabb*rhob*wght/t3
     3         **4+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = -2.1333333333333332d-2*t10*t9*wght+
     1         2.6666666666666666d-3*t11*t6*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = 2.6666666666666666d-3*t14*t8*wght-2
     1         .1333333333333332d-2*t12*t13*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhoa**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhoa**6.666666666666666d-1
            t6 = 1/t1**3
            t7 = 1/t5
            fnc(iq) = fnc(iq)-6.0d-3*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*t3*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Amat2(iq,D2_RA_RA) = 2.6666666666666666d-3*gammaaa*t4*t7*wgh
     1         t-2.1333333333333332d-2*gammaaa*t5*t6*wght+Amat2(iq,D2_RA
     2         _RA)
            Cmat2(iq,D2_RA_GAA) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = -2.1333333333333332d-2*gammaaa*t6*wg
     1         ht/t3-1.7777777777777776d-3*gammaaa*t4*wght/rhoa**1.66666
     2         66666666669d+0+8.533333333333333d-2*gammaaa*rhoa*wght/t1*
     3         *4+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = 2.6666666666666666d-3*t4*t7*wght-2.
     1         1333333333333332d-2*t5*t6*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**3.333333333333333d-1
            t4 = 1/t1**2
            t5 = rhob**6.666666666666666d-1
            t6 = 1/t1**3
            t7 = 1/t5
            fnc(iq) = fnc(iq)-6.0d-3*gammabb*t2*wght
            Amat(iq,D1_RB) = 8.0d-3*gammabb*t3*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t2*wght
            Amat2(iq,D2_RB_RB) = 2.6666666666666666d-3*gammabb*t4*t7*wgh
     1         t-2.1333333333333332d-2*gammabb*t5*t6*wght+Amat2(iq,D2_RB
     2         _RB)
            Cmat2(iq,D2_RB_GBB) = 8.0d-3*t3*t4*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = -2.1333333333333332d-2*gammabb*t6*wg
     1         ht/t3-1.7777777777777776d-3*gammabb*t4*wght/rhob**1.66666
     2         66666666669d+0+8.533333333333333d-2*gammabb*rhob*wght/t1*
     3         *4+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = 2.6666666666666666d-3*t4*t7*wght-2.
     1         1333333333333332d-2*t5*t6*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
