C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_b88.F 
C> The nwxcm_x_b88 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252
C>    \,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)
C>    \,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}
C>    \,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}
C>    \,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}
C>    +1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b88(param,tol_rho,ipol,nq,wght, 
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
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 1/(gammaaa*t7+1)**5.0d-1
            t9 = 1/t4**2
            fnc(iq) = fnc(iq)-8.4d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = (4.2d-3*gammaaa*t2*(-3.3599999999999997d-2*
     1         gammaaa*t8/rhoa**3.6666666666666664d+0-3.3599999999999997
     2         d-2*t1*t3*t6)*t9+5.599999999999999d-3*gammaaa*t5*t6)*wght
     3         +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t2*(1.26d-2*t7*t8+1.26d-2*
     1         t2*t3/t1)*t9-4.2d-3*t2*t5)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = gammabb**5.0d-1
            t7 = 1/rhob**1.3333333333333333d+0
            t8 = asinh(t6*t7)
            t9 = 2.52d-2*t6*t7*t8+1.0d+0
            t10 = 1/t9
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/rhoa**2.6666666666666666d+0
            t13 = 1/(gammaaa*t12+1)**5.0d-1
            t14 = 1/t4**2
            t15 = 1/rhob**2.3333333333333334d+0
            t16 = 1/rhob**2.6666666666666666d+0
            t17 = 1/(gammabb*t16+1)**5.0d-1
            t18 = 1/t9**2
            fnc(iq) = (-4.2d-3*gammabb*t10*t7-4.2d-3*gammaaa*t2*t5)*wght
     1         +fnc(iq)
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t11*t5+4.2d-3
     1         *gammaaa*t14*t2*(-3.3599999999999997d-2*t1*t11*t3-3.35999
     2         99999999997d-2*gammaaa*t13/rhoa**3.6666666666666664d+0))*
     3         wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.2d-3*gammabb*t18*t7*(-3.3599999999999997
     1         d-2*t15*t6*t8-3.3599999999999997d-2*gammabb*t17/rhob**3.6
     2         666666666666664d+0)+5.599999999999999d-3*gammabb*t10*t15)
     3         *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t14*t2*(1.26d-2*t2*t3/t1+1
     1         .26d-2*t12*t13)-4.2d-3*t2*t5)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (4.2d-3*gammabb*t18*t7*(1.26d-2*t7*t8/t6+1
     1         .26d-2*t16*t17)-4.2d-3*t10*t7)*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 1/(gammaaa*t7+1)**5.0d-1
            t9 = 1/t4**2
            fnc(iq) = fnc(iq)-4.2d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = 4.2d-3*gammaaa*t2*(-3.3599999999999997d-2*g
     1         ammaaa*t8/rhoa**3.6666666666666664d+0-3.3599999999999997d
     2         -2*t1*t3*t6)*t9*wght+5.599999999999999d-3*gammaaa*t5*t6*w
     3         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 4.2d-3*gammaaa*t2*(1.26d-2*t7*t8+1.26d-2*t
     1         2*t3/t1)*t9*wght-4.2d-3*t2*t5*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = 1/rhob**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhob**2.3333333333333334d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 1/(gammabb*t7+1)**5.0d-1
            t9 = 1/t4**2
            fnc(iq) = fnc(iq)-4.2d-3*gammabb*t2*t5*wght
            Amat(iq,D1_RB) = 4.2d-3*gammabb*t2*(-3.3599999999999997d-2*g
     1         ammabb*t8/rhob**3.6666666666666664d+0-3.3599999999999997d
     2         -2*t1*t3*t6)*t9*wght+5.599999999999999d-3*gammabb*t5*t6*w
     3         ght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 4.2d-3*gammabb*t2*(1.26d-2*t7*t8+1.26d-2*t
     1         2*t3/t1)*t9*wght-4.2d-3*t2*t5*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252
C>    \,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)
C>    \,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}
C>    \,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}
C>    \,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}
C>    +1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b88_d2(param,tol_rho,ipol,nq,wght, 
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
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = (gammaaa*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammaaa*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhoa**3.3333333333333337d+0
            t16 = 1/t4**3
            t17 = 1/t8**3
            fnc(iq) = fnc(iq)-8.4d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t5*t6+4.2d-3*
     1         gammaaa*t11*t12*t2)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t12*t14*t2-4.2d-3*t2*t5)*w
     1         ght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (4.2d-3*gammaaa*t12*t2*(1.6799999999999
     1         998d-1*gammaaa*t9/rhoa**4.666666666666667d+0+7.84d-2*t1*t
     2         15*t3-4.4799999999999995d-2*gammaaa**2*t17/rhoa**7.333333
     3         333333333d+0)-1.1199999999999999d-2*gammaaa*t11*t12*t6-1.
     4         3066666666666665d-2*gammaaa*t15*t5-8.4d-3*gammaaa*t11**2*
     5         t16*t2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t12*t2*(4.2d-3*gammaaa*(-5.04d-2*t10*
     1         t9-1.6799999999999998d-2*t13*t3*t6+1.6799999999999998d-2*
     2         gammaaa*t17/rhoa**6.333333333333333d+0)+4.2d-3*t11)+5.599
     3         999999999999d-3*t5*t6-5.599999999999999d-3*gammaaa*t12*t1
     4         4*t6-8.4d-3*gammaaa*t11*t14*t16*t2)*wght+Cmat2(iq,D2_RA_G
     5         AA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t12*t2*(4.2d-3*gammaaa*(6.3d-3*t7*t9
     1         /gammaaa-6.3d-3*t2*t3/t1**3-6.3d-3*t17/rhoa**5.3333333333
     2         33333d+0)+8.4d-3*t14)-8.4d-3*gammaaa*t14**2*t16*t2)*wght+
     3         Cmat2(iq,D2_GAA_GAA)
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
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = gammabb**5.0d-1
            t7 = 1/rhob**1.3333333333333333d+0
            t8 = asinh(t6*t7)
            t9 = 2.52d-2*t6*t7*t8+1.0d+0
            t10 = 1/t9
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/rhoa**2.6666666666666666d+0
            t13 = (gammaaa*t12+1)**5.0d-1
            t14 = 1/t13
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = -3.3599999999999997d-2*t1*t11*t3-3.3599999999999997d-2
     1         *gammaaa*t14*t15
            t17 = 1/t4**2
            t18 = 1/rhob**2.3333333333333334d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = (gammabb*t19+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/rhob**3.6666666666666664d+0
            t23 = -3.3599999999999997d-2*t18*t6*t8-3.3599999999999997d-2
     1         *gammabb*t21*t22
            t24 = 1/t9**2
            t25 = 1/t1
            t26 = 1.26d-2*t2*t25*t3+1.26d-2*t12*t14
            t27 = 1/t6
            t28 = 1.26d-2*t27*t7*t8+1.26d-2*t19*t21
            t29 = 1/rhoa**3.3333333333333337d+0
            t30 = 1/t4**3
            t31 = 1/t13**3
            t32 = 1/rhob**3.3333333333333337d+0
            t33 = 1/t9**3
            t34 = 1/t20**3
            fnc(iq) = (-4.2d-3*gammabb*t10*t7-4.2d-3*gammaaa*t2*t5)*wght
     1         +fnc(iq)
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t11*t5+4.2d-3
     1         *gammaaa*t16*t17*t2)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.2d-3*gammabb*t23*t24*t7+5.59999999999999
     1         9d-3*gammabb*t10*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t17*t2*t26-4.2d-3*t2*t5)*w
     1         ght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (4.2d-3*gammabb*t24*t28*t7-4.2d-3*t10*t7)*
     1         wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (-1.3066666666666665d-2*gammaaa*t29*t5+
     1         4.2d-3*gammaaa*t17*t2*(-4.4799999999999995d-2*gammaaa**2*
     2         t31/rhoa**7.333333333333333d+0+7.84d-2*t1*t29*t3+1.679999
     3         9999999998d-1*gammaaa*t14/rhoa**4.666666666666667d+0)-8.4
     4         d-3*gammaaa*t16**2*t2*t30-1.1199999999999999d-2*gammaaa*t
     5         11*t16*t17)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (4.2d-3*gammabb*t24*t7*(7.84d-2*t32*t6*
     1         t8-4.4799999999999995d-2*gammabb**2*t34/rhob**7.333333333
     2         333333d+0+1.6799999999999998d-1*gammabb*t21/rhob**4.66666
     3         6666666667d+0)-8.4d-3*gammabb*t23**2*t33*t7-1.30666666666
     4         66665d-2*gammabb*t10*t32-1.1199999999999999d-2*gammabb*t1
     5         8*t23*t24)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (5.599999999999999d-3*t11*t5+t17*t2*(4
     1         .2d-3*gammaaa*(1.6799999999999998d-2*gammaaa*t31/rhoa**6.
     2         333333333333333d+0-1.6799999999999998d-2*t11*t25*t3-5.04d
     3         -2*t14*t15)+4.2d-3*t16)-8.4d-3*gammaaa*t16*t2*t26*t30-5.5
     4         99999999999999d-3*gammaaa*t11*t17*t26)*wght+Cmat2(iq,D2_R
     5         A_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t24*t7*(4.2d-3*gammabb*(-1.6799999999
     1         999998d-2*t18*t27*t8+1.6799999999999998d-2*gammabb*t34/rh
     2         ob**6.333333333333333d+0-5.04d-2*t21*t22)+4.2d-3*t23)-8.4
     3         d-3*gammabb*t23*t28*t33*t7-5.599999999999999d-3*gammabb*t
     4         18*t24*t28+5.599999999999999d-3*t10*t18)*wght+Cmat2(iq,D2
     5         _RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t17*t2*(4.2d-3*gammaaa*(-6.3d-3*t31/
     1         rhoa**5.333333333333333d+0-6.3d-3*t2*t3/t1**3+6.3d-3*t12*
     2         t14/gammaaa)+8.4d-3*t26)-8.4d-3*gammaaa*t2*t26**2*t30)*wg
     3         ht+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t24*t7*(4.2d-3*gammabb*(-6.3d-3*t7*t
     1         8/t6**3-6.3d-3*t34/rhob**5.333333333333333d+0+6.3d-3*t19*
     2         t21/gammabb)+8.4d-3*t28)-8.4d-3*gammabb*t28**2*t33*t7)*wg
     3         ht+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = (gammaaa*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammaaa*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhoa**3.3333333333333337d+0
            t16 = 1/t4**3
            t17 = 1/t8**3
            fnc(iq) = fnc(iq)-4.2d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = 5.599999999999999d-3*gammaaa*t5*t6*wght+4.2
     1         d-3*gammaaa*t11*t12*t2*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -4.2d-3*t2*t5*wght+4.2d-3*gammaaa*t12*t14*
     1         t2*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = 4.2d-3*gammaaa*t12*t2*(1.67999999999999
     1         98d-1*gammaaa*t9/rhoa**4.666666666666667d+0+7.84d-2*t1*t1
     2         5*t3-4.4799999999999995d-2*gammaaa**2*t17/rhoa**7.3333333
     3         33333333d+0)*wght-1.1199999999999999d-2*gammaaa*t11*t12*t
     4         6*wght-1.3066666666666665d-2*gammaaa*t15*t5*wght-8.4d-3*g
     5         ammaaa*t11**2*t16*t2*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t12*t2*(4.2d-3*gammaaa*(-5.04d-2*t10*t
     1         9-1.6799999999999998d-2*t13*t3*t6+1.6799999999999998d-2*g
     2         ammaaa*t17/rhoa**6.333333333333333d+0)*wght+4.2d-3*t11*wg
     3         ht)+5.599999999999999d-3*t5*t6*wght-5.599999999999999d-3*
     4         gammaaa*t12*t14*t6*wght-8.4d-3*gammaaa*t11*t14*t16*t2*wgh
     5         t+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t12*t2*(4.2d-3*gammaaa*(6.3d-3*t7*t9/
     1         gammaaa-6.3d-3*t2*t3/t1**3-6.3d-3*t17/rhoa**5.33333333333
     2         3333d+0)*wght+8.4d-3*t14*wght)-8.4d-3*gammaaa*t14**2*t16*
     3         t2*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = 1/rhob**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhob**2.3333333333333334d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = (gammabb*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhob**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammabb*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhob**3.3333333333333337d+0
            t16 = 1/t4**3
            t17 = 1/t8**3
            fnc(iq) = fnc(iq)-4.2d-3*gammabb*t2*t5*wght
            Amat(iq,D1_RB) = 5.599999999999999d-3*gammabb*t5*t6*wght+4.2
     1         d-3*gammabb*t11*t12*t2*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -4.2d-3*t2*t5*wght+4.2d-3*gammabb*t12*t14*
     1         t2*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = 4.2d-3*gammabb*t12*t2*(1.67999999999999
     1         98d-1*gammabb*t9/rhob**4.666666666666667d+0+7.84d-2*t1*t1
     2         5*t3-4.4799999999999995d-2*gammabb**2*t17/rhob**7.3333333
     3         33333333d+0)*wght-1.1199999999999999d-2*gammabb*t11*t12*t
     4         6*wght-1.3066666666666665d-2*gammabb*t15*t5*wght-8.4d-3*g
     5         ammabb*t11**2*t16*t2*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t12*t2*(4.2d-3*gammabb*(-5.04d-2*t10*t
     1         9-1.6799999999999998d-2*t13*t3*t6+1.6799999999999998d-2*g
     2         ammabb*t17/rhob**6.333333333333333d+0)*wght+4.2d-3*t11*wg
     3         ht)+5.599999999999999d-3*t5*t6*wght-5.599999999999999d-3*
     4         gammabb*t12*t14*t6*wght-8.4d-3*gammabb*t11*t14*t16*t2*wgh
     5         t+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t12*t2*(4.2d-3*gammabb*(6.3d-3*t7*t9/
     1         gammabb-6.3d-3*t2*t3/t1**3-6.3d-3*t17/rhob**5.33333333333
     2         3333d+0)*wght+8.4d-3*t14*wght)-8.4d-3*gammabb*t14**2*t16*
     3         t2*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252
C>    \,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)
C>    \,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}
C>    \,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}
C>    \,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}
C>    +1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b88_d3(param,tol_rho,ipol,nq,wght, 
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
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = (gammaaa*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammaaa*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhoa**3.3333333333333337d+0
            t16 = t11**2
            t17 = 1/t4**3
            t18 = gammaaa**2
            t19 = 1/t8**3
            t20 = 1/rhoa**7.333333333333333d+0
            t21 = 1/rhoa**4.666666666666667d+0
            t22 = 1.6799999999999998d-1*gammaaa*t21*t9+7.84d-2*t1*t15*t3
     1         -4.4799999999999995d-2*t18*t19*t20
            t23 = 1/rhoa**6.333333333333333d+0
            t24 = -5.04d-2*t10*t9-1.6799999999999998d-2*t13*t3*t6+1.6799
     1         999999999998d-2*gammaaa*t19*t23
            t25 = 4.2d-3*gammaaa*t24+4.2d-3*t11
            t26 = t14**2
            t27 = 1/rhoa**5.333333333333333d+0
            t28 = 1/gammaaa
            t29 = 1/t1**3
            t30 = 6.3d-3*t28*t7*t9-6.3d-3*t2*t29*t3-6.3d-3*t19*t27
            t31 = 4.2d-3*gammaaa*t30+8.4d-3*t14
            t32 = 1/rhoa**4.333333333333333d+0
            t33 = 1/t4**4
            t34 = 1/t8**5
            fnc(iq) = fnc(iq)-8.4d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t5*t6+4.2d-3*
     1         gammaaa*t11*t12*t2)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t12*t14*t2-4.2d-3*t2*t5)*w
     1         ght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-1.1199999999999999d-2*gammaaa*t11*t12
     1         *t6-1.3066666666666665d-2*gammaaa*t15*t5+4.2d-3*gammaaa*t
     2         12*t2*t22-8.4d-3*gammaaa*t16*t17*t2)*wght+Amat2(iq,D2_RA_
     3         RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (5.599999999999999d-3*t5*t6-5.59999999
     1         9999999d-3*gammaaa*t12*t14*t6+t12*t2*t25-8.4d-3*gammaaa*t
     2         11*t14*t17*t2)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t12*t2*t31-8.4d-3*gammaaa*t17*t2*t26
     1         )*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (4.2d-3*gammaaa*t12*t2*(-8.885333333
     1         333333d-1*gammaaa*t9/rhoa**5.666666666666667d+0-1.7919999
     2         999999997d-1*gammaaa**3*t34/rhoa**11-2.6133333333333336d-
     3         1*t1*t3*t32+5.525333333333333d-1*t18*t19/rhoa**8.33333333
     4         3333334d+0)-1.6799999999999995d-2*gammaaa*t12*t22*t6+3.35
     5         9999999999999d-2*gammaaa*t16*t17*t6+4.355555555555556d-2*
     6         gammaaa*t32*t5+2.52d-2*gammaaa*t11**3*t2*t33-2.52d-2*gamm
     7         aaa*t11*t17*t2*t22+3.92d-2*gammaaa*t11*t12*t15)*wght+Amat
     8         3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t12*t2*(4.2d-3*gammaaa*(2.072d-1*t
     1         21*t9+6.719999999999998d-2*t18*t34/rhoa**10+3.92d-2*t13*t
     2         15*t3-1.7359999999999998d-1*gammaaa*t19*t20)+4.2d-3*t22)+
     3         t12*(-1.1199999999999999d-2*gammaaa*t24-1.119999999999999
     4         9d-2*t11)*t6+2.2399999999999998d-2*gammaaa*t11*t14*t17*t6
     5         -1.3066666666666665d-2*t15*t5+2.52d-2*gammaaa*t14*t16*t2*
     6         t33+t17*t2*(-1.6799999999999998d-2*gammaaa*t11*t24-8.4d-3
     7         *gammaaa*t14*t22-8.4d-3*t16)+1.3066666666666665d-2*gammaa
     8         a*t12*t14*t15)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t12*t2*(4.2d-3*gammaaa*(-8.4d-3*t
     1         10*t28*t9+8.4d-3*t29*t3*t6-2.52d-2*gammaaa*t34/rhoa**9+4.
     2         1999999999999993d-2*t19*t23)+8.4d-3*t24)+t12*(-5.59999999
     3         9999999d-3*gammaaa*t30-1.1199999999999999d-2*t14)*t6+1.11
     4         99999999999999d-2*gammaaa*t17*t26*t6+2.52d-2*gammaaa*t11*
     5         t2*t26*t33+t17*t2*(-8.4d-3*gammaaa*t11*t30-2*t14*t25-8.4d
     6         -3*gammaaa*t14*t24-8.4d-3*t11*t14))*wght+Cmat3(iq,D3_RA_G
     7         AA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t12*t2*(4.2d-3*gammaaa*(-9.45000
     1         0000000001d-3*t7*t9/t18+9.450000000000001d-3*t34/rhoa**8+
     2         9.450000000000001d-3*t2*t3/t1**5-3.15d-3*t19*t27*t28)+1.2
     3         6d-2*t30)+2.52d-2*gammaaa*t14**3*t2*t33+t17*t2*(-2*t14*t3
     4         1-1.6799999999999998d-2*gammaaa*t14*t30-8.4d-3*t26))*wght
     5         +Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = gammabb**5.0d-1
            t7 = 1/rhob**1.3333333333333333d+0
            t8 = asinh(t6*t7)
            t9 = 2.52d-2*t6*t7*t8+1.0d+0
            t10 = 1/t9
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/rhoa**2.6666666666666666d+0
            t13 = (gammaaa*t12+1)**5.0d-1
            t14 = 1/t13
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = -3.3599999999999997d-2*t1*t11*t3-3.3599999999999997d-2
     1         *gammaaa*t14*t15
            t17 = 1/t4**2
            t18 = 1/rhob**2.3333333333333334d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = (gammabb*t19+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/rhob**3.6666666666666664d+0
            t23 = -3.3599999999999997d-2*t18*t6*t8-3.3599999999999997d-2
     1         *gammabb*t21*t22
            t24 = 1/t9**2
            t25 = 1/t1
            t26 = 1.26d-2*t2*t25*t3+1.26d-2*t12*t14
            t27 = 1/t6
            t28 = 1.26d-2*t27*t7*t8+1.26d-2*t19*t21
            t29 = 1/rhoa**3.3333333333333337d+0
            t30 = t16**2
            t31 = 1/t4**3
            t32 = gammaaa**2
            t33 = 1/t13**3
            t34 = 1/rhoa**7.333333333333333d+0
            t35 = 1/rhoa**4.666666666666667d+0
            t36 = 1.6799999999999998d-1*gammaaa*t14*t35-4.47999999999999
     1         95d-2*t32*t33*t34+7.84d-2*t1*t29*t3
            t37 = 1/rhob**3.3333333333333337d+0
            t38 = t23**2
            t39 = 1/t9**3
            t40 = gammabb**2
            t41 = 1/t20**3
            t42 = 1/rhob**7.333333333333333d+0
            t43 = 1/rhob**4.666666666666667d+0
            t44 = 7.84d-2*t37*t6*t8+1.6799999999999998d-1*gammabb*t21*t4
     1         3-4.4799999999999995d-2*t40*t41*t42
            t45 = 1/rhoa**6.333333333333333d+0
            t46 = 1.6799999999999998d-2*gammaaa*t33*t45-1.67999999999999
     1         98d-2*t11*t25*t3-5.04d-2*t14*t15
            t47 = 4.2d-3*gammaaa*t46+4.2d-3*t16
            t48 = 1/rhob**6.333333333333333d+0
            t49 = -1.6799999999999998d-2*t18*t27*t8+1.6799999999999998d-
     1         2*gammabb*t41*t48-5.04d-2*t21*t22
            t50 = 4.2d-3*gammabb*t49+4.2d-3*t23
            t51 = t26**2
            t52 = 1/rhoa**5.333333333333333d+0
            t53 = 1/gammaaa
            t54 = 1/t1**3
            t55 = -6.3d-3*t2*t3*t54+6.3d-3*t12*t14*t53-6.3d-3*t33*t52
            t56 = 4.2d-3*gammaaa*t55+8.4d-3*t26
            t57 = t28**2
            t58 = 1/rhob**5.333333333333333d+0
            t59 = 1/gammabb
            t60 = 1/t6**3
            t61 = -6.3d-3*t60*t7*t8+6.3d-3*t19*t21*t59-6.3d-3*t41*t58
            t62 = 4.2d-3*gammabb*t61+8.4d-3*t28
            t63 = 1/rhoa**4.333333333333333d+0
            t64 = 1/t4**4
            t65 = 1/t13**5
            t66 = 1/rhob**4.333333333333333d+0
            t67 = 1/t9**4
            t68 = 1/t20**5
            fnc(iq) = (-4.2d-3*gammabb*t10*t7-4.2d-3*gammaaa*t2*t5)*wght
     1         +fnc(iq)
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t11*t5+4.2d-3
     1         *gammaaa*t16*t17*t2)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.2d-3*gammabb*t23*t24*t7+5.59999999999999
     1         9d-3*gammabb*t10*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t17*t2*t26-4.2d-3*t2*t5)*w
     1         ght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (4.2d-3*gammabb*t24*t28*t7-4.2d-3*t10*t7)*
     1         wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (-1.3066666666666665d-2*gammaaa*t29*t5+
     1         4.2d-3*gammaaa*t17*t2*t36-8.4d-3*gammaaa*t2*t30*t31-1.119
     2         9999999999999d-2*gammaaa*t11*t16*t17)*wght+Amat2(iq,D2_RA
     3         _RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (4.2d-3*gammabb*t24*t44*t7-8.4d-3*gamma
     1         bb*t38*t39*t7-1.3066666666666665d-2*gammabb*t10*t37-1.119
     2         9999999999999d-2*gammabb*t18*t23*t24)*wght+Amat2(iq,D2_RB
     3         _RB)
            Cmat2(iq,D2_RA_GAA) = (5.599999999999999d-3*t11*t5+t17*t2*t4
     1         7-8.4d-3*gammaaa*t16*t2*t26*t31-5.599999999999999d-3*gamm
     2         aaa*t11*t17*t26)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t24*t50*t7-8.4d-3*gammabb*t23*t28*t39
     1         *t7-5.599999999999999d-3*gammabb*t18*t24*t28+5.5999999999
     2         99999d-3*t10*t18)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t17*t2*t56-8.4d-3*gammaaa*t2*t31*t51
     1         )*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t24*t62*t7-8.4d-3*gammabb*t39*t57*t7
     1         )*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (4.2d-3*gammaaa*t17*t2*(-1.791999999
     1         9999997d-1*gammaaa**3*t65/rhoa**11-2.6133333333333336d-1*
     2         t1*t3*t63+5.525333333333333d-1*t32*t33/rhoa**8.3333333333
     3         33334d+0-8.885333333333333d-1*gammaaa*t14/rhoa**5.6666666
     4         66666667d+0)+2.52d-2*gammaaa*t16**3*t2*t64+4.355555555555
     5         556d-2*gammaaa*t5*t63-2.52d-2*gammaaa*t16*t2*t31*t36-1.67
     6         99999999999995d-2*gammaaa*t11*t17*t36+3.359999999999999d-
     7         2*gammaaa*t11*t30*t31+3.92d-2*gammaaa*t16*t17*t29)*wght+A
     8         mat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (4.2d-3*gammabb*t24*t7*(-2.613333333
     1         3333336d-1*t6*t66*t8-1.7919999999999997d-1*gammabb**3*t68
     2         /rhob**11+5.525333333333333d-1*t40*t41/rhob**8.3333333333
     3         33334d+0-8.885333333333333d-1*gammabb*t21/rhob**5.6666666
     4         66666667d+0)+2.52d-2*gammabb*t23**3*t67*t7-2.52d-2*gammab
     5         b*t23*t39*t44*t7+4.355555555555556d-2*gammabb*t10*t66-1.6
     6         799999999999995d-2*gammabb*t18*t24*t44+3.359999999999999d
     7         -2*gammabb*t18*t38*t39+3.92d-2*gammabb*t23*t24*t37)*wght+
     8         Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t17*t2*(4.2d-3*gammaaa*(6.71999999
     1         9999998d-2*t32*t65/rhoa**10+2.072d-1*t14*t35-1.7359999999
     2         999998d-1*gammaaa*t33*t34+3.92d-2*t25*t29*t3)+4.2d-3*t36)
     3         +2.52d-2*gammaaa*t2*t26*t30*t64-1.3066666666666665d-2*t29
     4         *t5+t2*t31*(-1.6799999999999998d-2*gammaaa*t16*t46-8.4d-3
     5         *gammaaa*t26*t36-8.4d-3*t30)+t11*t17*(-1.1199999999999999
     6         d-2*gammaaa*t46-1.1199999999999999d-2*t16)+2.239999999999
     7         9998d-2*gammaaa*t11*t16*t26*t31+1.3066666666666665d-2*gam
     8         maaa*t17*t26*t29)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t24*t7*(4.2d-3*gammabb*(3.92d-2*t2
     1         7*t37*t8+6.719999999999998d-2*t40*t68/rhob**10+2.072d-1*t
     2         21*t43-1.7359999999999998d-1*gammabb*t41*t42)+4.2d-3*t44)
     3         +2.52d-2*gammabb*t28*t38*t67*t7+t39*(-1.6799999999999998d
     4         -2*gammabb*t23*t49-8.4d-3*gammabb*t28*t44-8.4d-3*t38)*t7+
     5         t18*t24*(-1.1199999999999999d-2*gammabb*t49-1.11999999999
     6         99999d-2*t23)+2.2399999999999998d-2*gammabb*t18*t23*t28*t
     7         39+1.3066666666666665d-2*gammabb*t24*t28*t37-1.3066666666
     8         666665d-2*t10*t37)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t17*t2*(4.2d-3*gammaaa*(-2.52d-2*
     1         gammaaa*t65/rhoa**9+8.4d-3*t11*t3*t54-8.4d-3*t14*t15*t53+
     2         4.1999999999999993d-2*t33*t45)+8.4d-3*t46)+2.52d-2*gammaa
     3         a*t16*t2*t51*t64+t2*t31*(-8.4d-3*gammaaa*t16*t55-2*t26*t4
     4         7-8.4d-3*gammaaa*t26*t46-8.4d-3*t16*t26)+t11*t17*(-5.5999
     5         99999999999d-3*gammaaa*t55-1.1199999999999999d-2*t26)+1.1
     6         199999999999999d-2*gammaaa*t11*t31*t51)*wght+Cmat3(iq,D3_
     7         RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t24*t7*(4.2d-3*gammabb*(8.4d-3*t1
     1         8*t60*t8-2.52d-2*gammabb*t68/rhob**9-8.4d-3*t21*t22*t59+4
     2         .1999999999999993d-2*t41*t48)+8.4d-3*t49)+2.52d-2*gammabb
     3         *t23*t57*t67*t7+t39*(-8.4d-3*gammabb*t23*t61-2*t28*t50-8.
     4         4d-3*gammabb*t28*t49-8.4d-3*t23*t28)*t7+t18*t24*(-5.59999
     5         9999999999d-3*gammabb*t61-1.1199999999999999d-2*t28)+1.11
     6         99999999999999d-2*gammabb*t18*t39*t57)*wght+Cmat3(iq,D3_R
     7         B_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t17*t2*(4.2d-3*gammaaa*(9.450000
     1         000000001d-3*t65/rhoa**8-3.15d-3*t33*t52*t53-9.4500000000
     2         00001d-3*t12*t14/t32+9.450000000000001d-3*t2*t3/t1**5)+1.
     3         26d-2*t55)+2.52d-2*gammaaa*t2*t26**3*t64+t2*t31*(-2*t26*t
     4         56-1.6799999999999998d-2*gammaaa*t26*t55-8.4d-3*t51))*wgh
     5         t+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (t24*t7*(4.2d-3*gammabb*(9.450000
     1         000000001d-3*t7*t8/t6**5+9.450000000000001d-3*t68/rhob**8
     2         -3.15d-3*t41*t58*t59-9.450000000000001d-3*t19*t21/t40)+1.
     3         26d-2*t61)+2.52d-2*gammabb*t28**3*t67*t7+t39*(-2*t28*t62-
     4         1.6799999999999998d-2*gammabb*t28*t61-8.4d-3*t57)*t7)*wgh
     5         t+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = (gammaaa*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammaaa*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhoa**3.3333333333333337d+0
            t16 = t11**2
            t17 = 1/t4**3
            t18 = gammaaa**2
            t19 = 1/t8**3
            t20 = 1/rhoa**7.333333333333333d+0
            t21 = 1/rhoa**4.666666666666667d+0
            t22 = 1.6799999999999998d-1*gammaaa*t21*t9+7.84d-2*t1*t15*t3
     1         -4.4799999999999995d-2*t18*t19*t20
            t23 = 1/rhoa**6.333333333333333d+0
            t24 = -5.04d-2*t10*t9-1.6799999999999998d-2*t13*t3*t6+1.6799
     1         999999999998d-2*gammaaa*t19*t23
            t25 = 4.2d-3*gammaaa*t24*wght+4.2d-3*t11*wght
            t26 = t14**2
            t27 = 1/rhoa**5.333333333333333d+0
            t28 = 1/gammaaa
            t29 = 1/t1**3
            t30 = 6.3d-3*t28*t7*t9-6.3d-3*t2*t29*t3-6.3d-3*t19*t27
            t31 = 4.2d-3*gammaaa*t30*wght+8.4d-3*t14*wght
            t32 = 1/rhoa**4.333333333333333d+0
            t33 = 1/t4**4
            t34 = 1/t8**5
            fnc(iq) = fnc(iq)-4.2d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = 5.599999999999999d-3*gammaaa*t5*t6*wght+4.2
     1         d-3*gammaaa*t11*t12*t2*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -4.2d-3*t2*t5*wght+4.2d-3*gammaaa*t12*t14*
     1         t2*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = -1.1199999999999999d-2*gammaaa*t11*t12*
     1         t6*wght-1.3066666666666665d-2*gammaaa*t15*t5*wght+4.2d-3*
     2         gammaaa*t12*t2*t22*wght-8.4d-3*gammaaa*t16*t17*t2*wght+Am
     3         at2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 5.599999999999999d-3*t5*t6*wght-5.5999
     1         99999999999d-3*gammaaa*t12*t14*t6*wght-8.4d-3*gammaaa*t11
     2         *t14*t17*t2*wght+t12*t2*t25+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -8.4d-3*gammaaa*t17*t2*t26*wght+t12*t
     1         2*t31+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = 4.2d-3*gammaaa*t12*t2*(-8.8853333333
     1         33333d-1*gammaaa*t9/rhoa**5.666666666666667d+0-1.79199999
     2         99999997d-1*gammaaa**3*t34/rhoa**11-2.6133333333333336d-1
     3         *t1*t3*t32+5.525333333333333d-1*t18*t19/rhoa**8.333333333
     4         333334d+0)*wght-1.6799999999999995d-2*gammaaa*t12*t22*t6*
     5         wght+3.359999999999999d-2*gammaaa*t16*t17*t6*wght+4.35555
     6         5555555556d-2*gammaaa*t32*t5*wght+2.52d-2*gammaaa*t11**3*
     7         t2*t33*wght-2.52d-2*gammaaa*t11*t17*t2*t22*wght+3.92d-2*g
     8         ammaaa*t11*t12*t15*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t12*t2*(4.2d-3*gammaaa*(2.072d-1*t2
     1         1*t9+6.719999999999998d-2*t18*t34/rhoa**10+3.92d-2*t13*t1
     2         5*t3-1.7359999999999998d-1*gammaaa*t19*t20)*wght+4.2d-3*t
     3         22*wght)+t17*t2*(-1.6799999999999998d-2*gammaaa*t11*t24*w
     4         ght-8.4d-3*gammaaa*t14*t22*wght-8.4d-3*t16*wght)+t12*t6*(
     5         -1.1199999999999999d-2*gammaaa*t24*wght-1.119999999999999
     6         9d-2*t11*wght)+2.2399999999999998d-2*gammaaa*t11*t14*t17*
     7         t6*wght-1.3066666666666665d-2*t15*t5*wght+2.52d-2*gammaaa
     8         *t14*t16*t2*t33*wght+1.3066666666666665d-2*gammaaa*t12*t1
     9         4*t15*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t12*t2*(4.2d-3*gammaaa*(-8.4d-3*t1
     1         0*t28*t9+8.4d-3*t29*t3*t6-2.52d-2*gammaaa*t34/rhoa**9+4.1
     2         999999999999993d-2*t19*t23)*wght+8.4d-3*t24*wght)+t17*t2*
     3         (-8.4d-3*gammaaa*t11*t30*wght-8.4d-3*gammaaa*t14*t24*wght
     4         -8.4d-3*t11*t14*wght-2*t14*t25)+t12*t6*(-5.59999999999999
     5         9d-3*gammaaa*t30*wght-1.1199999999999999d-2*t14*wght)+1.1
     6         199999999999999d-2*gammaaa*t17*t26*t6*wght+2.52d-2*gammaa
     7         a*t11*t2*t26*t33*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t12*t2*(4.2d-3*gammaaa*(-9.450000
     1         000000001d-3*t7*t9/t18+9.450000000000001d-3*t34/rhoa**8+9
     2         .450000000000001d-3*t2*t3/t1**5-3.15d-3*t19*t27*t28)*wght
     3         +1.26d-2*t30*wght)+t17*t2*(-1.6799999999999998d-2*gammaaa
     4         *t14*t30*wght-8.4d-3*t26*wght-2*t14*t31)+2.52d-2*gammaaa*
     5         t14**3*t2*t33*wght+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = 1/rhob**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhob**2.3333333333333334d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = (gammabb*t7+1)**5.0d-1
            t9 = 1/t8
            t10 = 1/rhob**3.6666666666666664d+0
            t11 = -3.3599999999999997d-2*gammabb*t10*t9-3.35999999999999
     1         97d-2*t1*t3*t6
            t12 = 1/t4**2
            t13 = 1/t1
            t14 = 1.26d-2*t7*t9+1.26d-2*t13*t2*t3
            t15 = 1/rhob**3.3333333333333337d+0
            t16 = t11**2
            t17 = 1/t4**3
            t18 = gammabb**2
            t19 = 1/t8**3
            t20 = 1/rhob**7.333333333333333d+0
            t21 = 1/rhob**4.666666666666667d+0
            t22 = 1.6799999999999998d-1*gammabb*t21*t9+7.84d-2*t1*t15*t3
     1         -4.4799999999999995d-2*t18*t19*t20
            t23 = 1/rhob**6.333333333333333d+0
            t24 = -5.04d-2*t10*t9-1.6799999999999998d-2*t13*t3*t6+1.6799
     1         999999999998d-2*gammabb*t19*t23
            t25 = 4.2d-3*gammabb*t24*wght+4.2d-3*t11*wght
            t26 = t14**2
            t27 = 1/rhob**5.333333333333333d+0
            t28 = 1/gammabb
            t29 = 1/t1**3
            t30 = 6.3d-3*t28*t7*t9-6.3d-3*t2*t29*t3-6.3d-3*t19*t27
            t31 = 4.2d-3*gammabb*t30*wght+8.4d-3*t14*wght
            t32 = 1/rhob**4.333333333333333d+0
            t33 = 1/t4**4
            t34 = 1/t8**5
            fnc(iq) = fnc(iq)-4.2d-3*gammabb*t2*t5*wght
            Amat(iq,D1_RB) = 5.599999999999999d-3*gammabb*t5*t6*wght+4.2
     1         d-3*gammabb*t11*t12*t2*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -4.2d-3*t2*t5*wght+4.2d-3*gammabb*t12*t14*
     1         t2*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = -1.1199999999999999d-2*gammabb*t11*t12*
     1         t6*wght-1.3066666666666665d-2*gammabb*t15*t5*wght+4.2d-3*
     2         gammabb*t12*t2*t22*wght-8.4d-3*gammabb*t16*t17*t2*wght+Am
     3         at2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 5.599999999999999d-3*t5*t6*wght-5.5999
     1         99999999999d-3*gammabb*t12*t14*t6*wght-8.4d-3*gammabb*t11
     2         *t14*t17*t2*wght+t12*t2*t25+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -8.4d-3*gammabb*t17*t2*t26*wght+t12*t
     1         2*t31+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = 4.2d-3*gammabb*t12*t2*(-8.8853333333
     1         33333d-1*gammabb*t9/rhob**5.666666666666667d+0-1.79199999
     2         99999997d-1*gammabb**3*t34/rhob**11-2.6133333333333336d-1
     3         *t1*t3*t32+5.525333333333333d-1*t18*t19/rhob**8.333333333
     4         333334d+0)*wght-1.6799999999999995d-2*gammabb*t12*t22*t6*
     5         wght+3.359999999999999d-2*gammabb*t16*t17*t6*wght+4.35555
     6         5555555556d-2*gammabb*t32*t5*wght+2.52d-2*gammabb*t11**3*
     7         t2*t33*wght-2.52d-2*gammabb*t11*t17*t2*t22*wght+3.92d-2*g
     8         ammabb*t11*t12*t15*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t12*t2*(4.2d-3*gammabb*(2.072d-1*t2
     1         1*t9+6.719999999999998d-2*t18*t34/rhob**10+3.92d-2*t13*t1
     2         5*t3-1.7359999999999998d-1*gammabb*t19*t20)*wght+4.2d-3*t
     3         22*wght)+t17*t2*(-1.6799999999999998d-2*gammabb*t11*t24*w
     4         ght-8.4d-3*gammabb*t14*t22*wght-8.4d-3*t16*wght)+t12*t6*(
     5         -1.1199999999999999d-2*gammabb*t24*wght-1.119999999999999
     6         9d-2*t11*wght)+2.2399999999999998d-2*gammabb*t11*t14*t17*
     7         t6*wght-1.3066666666666665d-2*t15*t5*wght+2.52d-2*gammabb
     8         *t14*t16*t2*t33*wght+1.3066666666666665d-2*gammabb*t12*t1
     9         4*t15*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t12*t2*(4.2d-3*gammabb*(-8.4d-3*t1
     1         0*t28*t9+8.4d-3*t29*t3*t6-2.52d-2*gammabb*t34/rhob**9+4.1
     2         999999999999993d-2*t19*t23)*wght+8.4d-3*t24*wght)+t17*t2*
     3         (-8.4d-3*gammabb*t11*t30*wght-8.4d-3*gammabb*t14*t24*wght
     4         -8.4d-3*t11*t14*wght-2*t14*t25)+t12*t6*(-5.59999999999999
     5         9d-3*gammabb*t30*wght-1.1199999999999999d-2*t14*wght)+1.1
     6         199999999999999d-2*gammabb*t17*t26*t6*wght+2.52d-2*gammab
     7         b*t11*t2*t26*t33*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t12*t2*(4.2d-3*gammabb*(-9.450000
     1         000000001d-3*t7*t9/t18+9.450000000000001d-3*t34/rhob**8+9
     2         .450000000000001d-3*t2*t3/t1**5-3.15d-3*t19*t27*t28)*wght
     3         +1.26d-2*t30*wght)+t17*t2*(-1.6799999999999998d-2*gammabb
     4         *t14*t30*wght-8.4d-3*t26*wght-2*t14*t31)+2.52d-2*gammabb*
     5         t14**3*t2*t33*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
