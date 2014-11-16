C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_tpss03.F 
C> The nwxcm_x_tpss03 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_tpss03 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sigma_{\alpha\alpha}^2\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& \left({{{\it t_6}\,\sigma_{\alpha\alpha}}\over{\tau_\alpha}}\right)^{2.0}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 2.0\,\tau_\alpha-0.25\,{\it t_6}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{10}} &=& .1097020523068037\,{\it t_8}\,{\it t_9}-1.0\\\\ 
C>   {\it t_{11}} &=& {{0.45\,{\it t_{10}}}\over{\left(0.0438808209227215\,{\it t_8}\,{\it t_9}\,{\it t_{10}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_3}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\beta^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \sigma_{\beta\beta}^2\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{16}} &=& \left({{{\it t_{15}}\,\sigma_{\beta\beta}}\over{\tau_\beta}}\right)^{2.0}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 2.0\,\tau_\beta-0.25\,{\it t_{15}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{19}} &=& .1097020523068037\,{\it t_{17}}\,{\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {{0.45\,{\it t_{19}}}\over{\left(0.0438808209227215\,{\it t_{17}}\,{\it t_{18}}\,{\it t_{19}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{12}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^2\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{26}} &=& \left({{{\it t_{25}}\,\sigma_{ss}}\over{\tau_s}}\right)^{2.0}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& 2.0\,\tau_s-0.25\,{\it t_{25}}\,\sigma_{ss}\\\\ 
C>   {\it t_{29}} &=& .1097020523068037\,{\it t_{27}}\,{\it t_{28}}-1.0\\\\ 
C>   {\it t_{30}} &=& {{0.45\,{\it t_{29}}}\over{\left(0.0438808209227215\,{\it t_{27}}\,{\it t_{28}}\,{\it t_{29}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{22}}\,\sigma_{ss}\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{20}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{16}}+1.35388578153653 \times 10^{-4}\,{\it t_{13}}\,{\it t_{14}}\right)^{0.5}\,{\it t_{20}}+.01645530784602056\,{\it t_{12}}\,\sigma_{\beta\beta}\,\left({{0.02485875\,{\it t_{16}}}\over{\left(0.015625\,{\it t_{16}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{16}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{13}}\,{\it t_{14}}\right)}\over{\left(.02040060021774676\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{11}}^2-.1802469135802469\,\left(0.0028125\,{\it t_7}+1.35388578153653 \times 10^{-4}\,{\it t_4}\,{\it t_5}\right)^{0.5}\,{\it t_{11}}+.01645530784602056\,{\it t_3}\,\sigma_{\alpha\alpha}\,\left({{0.02485875\,{\it t_7}}\over{\left(0.015625\,{\it t_7}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_7}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8}}+5.133173416322327 \times 10^{-6}\,{\it t_4}\,{\it t_5}\right)}\over{\left(.02040060021774676\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_{21}}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{30}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{26}}+1.35388578153653 \times 10^{-4}\,{\it t_{23}}\,{\it t_{24}}\right)^{0.5}\,{\it t_{30}}+.01645530784602056\,{\it t_{22}}\,\sigma_{ss}\,\left({{0.02485875\,{\it t_{26}}}\over{\left(0.015625\,{\it t_{26}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{26}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{ss}^3}\over{\rho_s^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}\right)}\over{\left(.02040060021774676\,{\it t_{22}}\,\sigma_{ss}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_{21}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] J Tao, JP Perdew, VN Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, 146401 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.91.146401 "> 
C> 10.1103/PhysRevLett.91.146401 </a> 
C> 
      subroutine nwxcm_x_tpss03(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 2.040060021774676d-2*gammaaa*t2+1.0d+0
            t4 = 1/t3**2.0d+0
            t5 = gammaaa**3
            t6 = 1/rhoa**8
            t7 = gammaaa**2
            t8 = 1/rhoa**5.333333333333333d+0
            t9 = 1/rhoa
            t10 = 1/taua
            t11 = gammaaa*t10*t9
            t12 = t11**2.0d+0
            t13 = 1.5625d-2*t12+1.0d+0
            t14 = 1/t13**2.0d+0
            t15 = 2.485875d-2*t12*t14+1.234567901234568d-1
            t16 = (1.35388578153653d-4*t7*t8+2.8125d-3*t12)**5.0d-1
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 2.0d+0*taua-2.5d-1*gammaaa*t9
            t19 = 1.0970205230680374d-1*t17*t18-1.0d+0
            t20 = 4.38808209227215d-2*t17*t18*t19+1.0d+0
            t21 = 1/t20**5.0d-1
            t22 = 4.5d-1*t19*t21+1.0970205230680376d-2*gammaaa*t2
            t23 = 5.133173416322327d-6*t7*t8+1.5033019185692253d-6*t5*t6
     1         +7.209876543209877d-2*t22**2-1.802469135802469d-1*t16*t22
     2         +1.6455307846020564d-2*gammaaa*t15*t2+1.7218861679299946d
     3         -3*t12
            t24 = 1.2437810945273631d+0*t23*t4+1.0d+0
            t25 = 8.04d-1-8.04d-1/t24
            t26 = rhoa**3.333333333333333d-1
            t27 = 1/t24**2
            t28 = 1/rhoa**6.333333333333333d+0
            t29 = 1/rhoa**3.6666666666666664d+0
            t30 = 1/rhoa**3
            t31 = 1/taua**2
            t32 = 1/rhoa**2
            t33 = 1/t13**3.0d+0
            t34 = t11**3.0d+0
            t35 = 2.7425513076700936d-2*gammaaa*t29-1.828367538446729d-1
     1         *t18*t2
            t36 = 1/t20**1.5d+0
            t37 = -2.25d-1*t19*(4.38808209227215d-2*t17*t18*t35+1.097020
     1         5230680374d-2*gammaaa*t19*t29-7.313470153786916d-2*t18*t1
     2         9*t2)*t36+4.5d-1*t21*t35-2.925388061514767d-2*gammaaa*t29
            t38 = 1/t16
            t39 = 1/t3**3.0d+0
            t40 = -2.25d-1*t19*(-1.0970205230680374d-2*t19*t2-1.20345402
     1         80324706d-3*t18/rhoa**4.333333333333333d+0)*t36-1.2341480
     2         884515421d-2*t2*t21+1.0970205230680376d-2*t2
            t41 = -1.4962853840493523d+0*t1*t27*(1.2437810945273631d+0*t
     1         4*(1.6455307846020564d-2*gammaaa*t2*(4.97175d-2*gammaaa*t
     2         14*t31*t32-1.553671875d-3*t10*t33*t34*t9)-9.0123456790123
     3         45d-2*t22*t38*(2.70777156307306d-4*gammaaa*t8+5.625d-3*ga
     4         mmaaa*t31*t32)+1.0266346832644656d-5*gammaaa*t8+4.5099057
     5         55707676d-6*t6*t7+1.4419753086419754d-1*t22*t40-1.8024691
     6         35802469d-1*t16*t40+3.443772335859989d-3*gammaaa*t31*t32+
     7         1.6455307846020564d-2*t15*t2)-5.074776173568845d-2*t2*t23
     8         *t39)*wght
            t42 = 1/taua**3
            t43 = 9.873184707612337d-2*t17*t21-2.25d-1*t19*(8.7761641845
     1         443d-2*t17*t19+9.627632224259765d-3*t18/rhoa**3.333333333
     2         3333337d+0)*t36
            fnc(iq) = (-1.8610514726981994d+0*t1*t25-1.8610514726981994d
     1         +0*t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.4962853840493523d+0*t1*t27*(1.243781094
     1         5273631d+0*t4*(1.6455307846020564d-2*gammaaa*t2*(1.553671
     2         875d-3*gammaaa*t10*t32*t33*t34-4.97175d-2*t14*t30*t31*t7)
     3         -9.012345679012345d-2*t22*t38*(-5.625d-3*t30*t31*t7-7.220
     4         724168194825d-4*t28*t7)-3.443772335859989d-3*t30*t31*t7-2
     5         .737692488705241d-5*t28*t7-1.2026415348553804d-5*t5/rhoa*
     6         *9+1.4419753086419754d-1*t22*t37-1.802469135802469d-1*t16
     7         *t37-4.3880820922721503d-2*gammaaa*t15*t29)+1.35327364628
     8         50253d-1*gammaaa*t23*t29*t39)-2.481401963597599d+0*t25*t2
     9         6-2.481401963597599d+0*t26)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t41+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t41+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = Amat(iq,D1_TA)-1.8610514726981994d+0*t1*t27
     1         *t4*(1.6455307846020564d-2*gammaaa*t2*(1.553671875d-3*gam
     2         maaa*t31*t33*t34*t9-4.97175d-2*t14*t32*t42*t7)+5.06944444
     3         4444443d-4*t22*t32*t38*t42*t7-3.443772335859989d-3*t32*t4
     4         2*t7+1.4419753086419754d-1*t22*t43-1.802469135802469d-1*t
     5         16*t43)*wght
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = rhob**1.3333333333333333d+0
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 2.040060021774676d-2*gammaaa*t3+1.0d+0
            t5 = 1/t4**2.0d+0
            t6 = gammaaa**3
            t7 = 1/rhoa**8
            t8 = gammaaa**2
            t9 = 1/rhoa**5.333333333333333d+0
            t10 = 1/rhoa
            t11 = 1/taua
            t12 = gammaaa*t10*t11
            t13 = t12**2.0d+0
            t14 = 1.5625d-2*t13+1.0d+0
            t15 = 1/t14**2.0d+0
            t16 = 2.485875d-2*t13*t15+1.234567901234568d-1
            t17 = (1.35388578153653d-4*t8*t9+2.8125d-3*t13)**5.0d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 2.0d+0*taua-2.5d-1*gammaaa*t10
            t20 = 1.0970205230680374d-1*t18*t19-1.0d+0
            t21 = 4.38808209227215d-2*t18*t19*t20+1.0d+0
            t22 = 1/t21**5.0d-1
            t23 = 1.0970205230680376d-2*gammaaa*t3+4.5d-1*t20*t22
            t24 = 5.133173416322327d-6*t8*t9+1.5033019185692253d-6*t6*t7
     1         +1.6455307846020564d-2*gammaaa*t16*t3+7.209876543209877d-
     2         2*t23**2-1.802469135802469d-1*t17*t23+1.7218861679299946d
     3         -3*t13
            t25 = 1.2437810945273631d+0*t24*t5+1.0d+0
            t26 = 8.04d-1-8.04d-1/t25
            t27 = 1/rhob**2.6666666666666666d+0
            t28 = 2.040060021774676d-2*gammabb*t27+1.0d+0
            t29 = 1/t28**2.0d+0
            t30 = gammabb**3
            t31 = 1/rhob**8
            t32 = gammabb**2
            t33 = 1/rhob**5.333333333333333d+0
            t34 = 1/rhob
            t35 = 1/taub
            t36 = gammabb*t34*t35
            t37 = t36**2.0d+0
            t38 = 1.5625d-2*t37+1.0d+0
            t39 = 1/t38**2.0d+0
            t40 = 2.485875d-2*t37*t39+1.234567901234568d-1
            t41 = (2.8125d-3*t37+1.35388578153653d-4*t32*t33)**5.0d-1
            t42 = 1/rhob**1.6666666666666669d+0
            t43 = 2.0d+0*taub-2.5d-1*gammabb*t34
            t44 = 1.0970205230680374d-1*t42*t43-1.0d+0
            t45 = 4.38808209227215d-2*t42*t43*t44+1.0d+0
            t46 = 1/t45**5.0d-1
            t47 = 4.5d-1*t44*t46+1.0970205230680376d-2*gammabb*t27
            t48 = 7.209876543209877d-2*t47**2-1.802469135802469d-1*t41*t
     1         47+1.6455307846020564d-2*gammabb*t27*t40+1.72188616792999
     2         46d-3*t37+5.133173416322327d-6*t32*t33+1.5033019185692253
     3         d-6*t30*t31
            t49 = 1.2437810945273631d+0*t29*t48+1.0d+0
            t50 = 8.04d-1-8.04d-1/t49
            t51 = rhoa**3.333333333333333d-1
            t52 = 1/t25**2
            t53 = 1/rhoa**6.333333333333333d+0
            t54 = 1/rhoa**3.6666666666666664d+0
            t55 = 1/rhoa**3
            t56 = 1/taua**2
            t57 = 1/rhoa**2
            t58 = 1/t14**3.0d+0
            t59 = t12**3.0d+0
            t60 = 2.7425513076700936d-2*gammaaa*t54-1.828367538446729d-1
     1         *t19*t3
            t61 = 1/t21**1.5d+0
            t62 = -2.25d-1*t20*(4.38808209227215d-2*t18*t19*t60+1.097020
     1         5230680374d-2*gammaaa*t20*t54-7.313470153786916d-2*t19*t2
     2         0*t3)*t61+4.5d-1*t22*t60-2.925388061514767d-2*gammaaa*t54
            t63 = 1/t17
            t64 = 1/t4**3.0d+0
            t65 = rhob**3.333333333333333d-1
            t66 = 1/t49**2
            t67 = 1/rhob**6.333333333333333d+0
            t68 = 1/rhob**3.6666666666666664d+0
            t69 = 1/rhob**3
            t70 = 1/taub**2
            t71 = 1/rhob**2
            t72 = 1/t38**3.0d+0
            t73 = t36**3.0d+0
            t74 = 2.7425513076700936d-2*gammabb*t68-1.828367538446729d-1
     1         *t27*t43
            t75 = 1/t45**1.5d+0
            t76 = -2.25d-1*t44*(4.38808209227215d-2*t42*t43*t74+1.097020
     1         5230680374d-2*gammabb*t44*t68-7.313470153786916d-2*t27*t4
     2         3*t44)*t75+4.5d-1*t46*t74-2.925388061514767d-2*gammabb*t6
     3         8
            t77 = 1/t41
            t78 = 1/t28**3.0d+0
            t79 = -2.25d-1*t20*(-1.0970205230680374d-2*t20*t3-1.20345402
     1         80324706d-3*t19/rhoa**4.333333333333333d+0)*t61-1.2341480
     2         884515421d-2*t22*t3+1.0970205230680376d-2*t3
            t80 = -2.25d-1*t44*(-1.0970205230680374d-2*t27*t44-1.2034540
     1         280324706d-3*t43/rhob**4.333333333333333d+0)*t75-1.234148
     2         0884515421d-2*t27*t46+1.0970205230680376d-2*t27
            t81 = 1/taua**3
            t82 = 9.873184707612337d-2*t18*t22-2.25d-1*t20*(8.7761641845
     1         443d-2*t18*t20+9.627632224259765d-3*t19/rhoa**3.333333333
     2         3333337d+0)*t61
            t83 = 1/taub**3
            t84 = 9.873184707612337d-2*t42*t46-2.25d-1*t44*(8.7761641845
     1         443d-2*t42*t44+9.627632224259765d-3*t43/rhob**3.333333333
     2         3333337d+0)*t75
            fnc(iq) = (-9.305257363490997d-1*t2*t50-9.305257363490997d-1
     1         *t1*t26-9.305257363490997d-1*t2-9.305257363490997d-1*t1)*
     2         wght+fnc(iq)
            Amat(iq,D1_RA) = (-7.481426920246762d-1*t1*t52*(1.2437810945
     1         273631d+0*t5*(1.6455307846020564d-2*gammaaa*t3*(1.5536718
     2         75d-3*gammaaa*t11*t57*t58*t59-4.97175d-2*t15*t55*t56*t8)-
     3         9.012345679012345d-2*t23*t63*(-5.625d-3*t55*t56*t8-7.2207
     4         24168194825d-4*t53*t8)-3.443772335859989d-3*t55*t56*t8-2.
     5         737692488705241d-5*t53*t8+1.4419753086419754d-1*t23*t62-1
     6         .802469135802469d-1*t17*t62-1.2026415348553804d-5*t6/rhoa
     7         **9-4.3880820922721503d-2*gammaaa*t16*t54)+1.353273646285
     8         0253d-1*gammaaa*t24*t54*t64)-1.2407009817987995d+0*t26*t5
     9         1-1.2407009817987995d+0*t51)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-7.481426920246762d-1*t2*t66*(1.3532736462
     1         850253d-1*gammabb*t48*t68*t78+1.2437810945273631d+0*t29*(
     2         -9.012345679012345d-2*t47*(-5.625d-3*t32*t69*t70-7.220724
     3         168194825d-4*t32*t67)*t77+1.4419753086419754d-1*t47*t76-1
     4         .802469135802469d-1*t41*t76+1.6455307846020564d-2*gammabb
     5         *t27*(1.553671875d-3*gammabb*t35*t71*t72*t73-4.97175d-2*t
     6         32*t39*t69*t70)-3.443772335859989d-3*t32*t69*t70-4.388082
     7         0922721503d-2*gammabb*t40*t68-2.737692488705241d-5*t32*t6
     8         7-1.2026415348553804d-5*t30/rhob**9))-1.2407009817987995d
     9         +0*t50*t65-1.2407009817987995d+0*t65)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.481426920246762d-1*t1*t5
     1         2*(1.2437810945273631d+0*t5*(-9.012345679012345d-2*t23*t6
     2         3*(2.70777156307306d-4*gammaaa*t9+5.625d-3*gammaaa*t56*t5
     3         7)+1.0266346832644656d-5*gammaaa*t9+4.509905755707676d-6*
     4         t7*t8+1.4419753086419754d-1*t23*t79-1.802469135802469d-1*
     5         t17*t79+1.6455307846020564d-2*gammaaa*t3*(4.97175d-2*gamm
     6         aaa*t15*t56*t57-1.553671875d-3*t10*t11*t58*t59)+3.4437723
     7         35859989d-3*gammaaa*t56*t57+1.6455307846020564d-2*t16*t3)
     8         -5.074776173568845d-2*t24*t3*t64)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.481426920246762d-1*t2*t6
     1         6*(1.2437810945273631d+0*t29*(1.4419753086419754d-1*t47*t
     2         80-1.802469135802469d-1*t41*t80-9.012345679012345d-2*t47*
     3         (5.625d-3*gammabb*t70*t71+2.70777156307306d-4*gammabb*t33
     4         )*t77+1.6455307846020564d-2*gammabb*t27*(4.97175d-2*gamma
     5         bb*t39*t70*t71-1.553671875d-3*t34*t35*t72*t73)+3.44377233
     6         5859989d-3*gammabb*t70*t71+1.6455307846020564d-2*t27*t40+
     7         1.0266346832644656d-5*gammabb*t33+4.509905755707676d-6*t3
     8         1*t32)-5.074776173568845d-2*t27*t48*t78)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t5*t
     1         52*(1.4419753086419754d-1*t23*t82-1.802469135802469d-1*t1
     2         7*t82+1.6455307846020564d-2*gammaaa*t3*(1.553671875d-3*ga
     3         mmaaa*t10*t56*t58*t59-4.97175d-2*t15*t57*t8*t81)+5.069444
     4         444444443d-4*t23*t57*t63*t8*t81-3.443772335859989d-3*t57*
     5         t8*t81)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t2*t29*
     1         t66*(1.4419753086419754d-1*t47*t84-1.802469135802469d-1*t
     2         41*t84+1.6455307846020564d-2*gammabb*t27*(1.553671875d-3*
     3         gammabb*t34*t70*t72*t73-4.97175d-2*t32*t39*t71*t83)+5.069
     4         444444444443d-4*t32*t47*t71*t77*t83-3.443772335859989d-3*
     5         t32*t71*t83)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 2.040060021774676d-2*gammaaa*t2+1.0d+0
            t4 = 1/t3**2.0d+0
            t5 = gammaaa**3
            t6 = 1/rhoa**8
            t7 = gammaaa**2
            t8 = 1/rhoa**5.333333333333333d+0
            t9 = 1/rhoa
            t10 = 1/taua
            t11 = gammaaa*t10*t9
            t12 = t11**2.0d+0
            t13 = 1.5625d-2*t12+1.0d+0
            t14 = 1/t13**2.0d+0
            t15 = 2.485875d-2*t12*t14+1.234567901234568d-1
            t16 = (1.35388578153653d-4*t7*t8+2.8125d-3*t12)**5.0d-1
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 2.0d+0*taua-2.5d-1*gammaaa*t9
            t19 = 1.0970205230680374d-1*t17*t18-1.0d+0
            t20 = 4.38808209227215d-2*t17*t18*t19+1.0d+0
            t21 = 1/t20**5.0d-1
            t22 = 4.5d-1*t19*t21+1.0970205230680376d-2*gammaaa*t2
            t23 = 5.133173416322327d-6*t7*t8+1.5033019185692253d-6*t5*t6
     1         +7.209876543209877d-2*t22**2-1.802469135802469d-1*t16*t22
     2         +1.6455307846020564d-2*gammaaa*t15*t2+1.7218861679299946d
     3         -3*t12
            t24 = 1.2437810945273631d+0*t23*t4+1.0d+0
            t25 = 8.04d-1-8.04d-1/t24
            t26 = rhoa**3.333333333333333d-1
            t27 = 1/t24**2
            t28 = 1/rhoa**6.333333333333333d+0
            t29 = 1/rhoa**3.6666666666666664d+0
            t30 = 1/rhoa**3
            t31 = 1/taua**2
            t32 = 1/rhoa**2
            t33 = 1/t13**3.0d+0
            t34 = t11**3.0d+0
            t35 = 2.7425513076700936d-2*gammaaa*t29-1.828367538446729d-1
     1         *t18*t2
            t36 = 1/t20**1.5d+0
            t37 = -2.25d-1*t19*(4.38808209227215d-2*t17*t18*t35+1.097020
     1         5230680374d-2*gammaaa*t19*t29-7.313470153786916d-2*t18*t1
     2         9*t2)*t36+4.5d-1*t21*t35-2.925388061514767d-2*gammaaa*t29
            t38 = 1/t16
            t39 = 1/t3**3.0d+0
            t40 = -2.25d-1*t19*(-1.0970205230680374d-2*t19*t2-1.20345402
     1         80324706d-3*t18/rhoa**4.333333333333333d+0)*t36-1.2341480
     2         884515421d-2*t2*t21+1.0970205230680376d-2*t2
            t41 = 1/taua**3
            t42 = 9.873184707612337d-2*t17*t21-2.25d-1*t19*(8.7761641845
     1         443d-2*t17*t19+9.627632224259765d-3*t18/rhoa**3.333333333
     2         3333337d+0)*t36
            fnc(iq) = (-9.305257363490997d-1*t1*t25-9.305257363490997d-1
     1         *t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-7.481426920246762d-1*t1*t27*(1.2437810945
     1         273631d+0*t4*(1.6455307846020564d-2*gammaaa*t2*(1.5536718
     2         75d-3*gammaaa*t10*t32*t33*t34-4.97175d-2*t14*t30*t31*t7)-
     3         9.012345679012345d-2*t22*t38*(-5.625d-3*t30*t31*t7-7.2207
     4         24168194825d-4*t28*t7)-3.443772335859989d-3*t30*t31*t7-2.
     5         737692488705241d-5*t28*t7-1.2026415348553804d-5*t5/rhoa**
     6         9+1.4419753086419754d-1*t22*t37-1.802469135802469d-1*t16*
     7         t37-4.3880820922721503d-2*gammaaa*t15*t29)+1.353273646285
     8         0253d-1*gammaaa*t23*t29*t39)-1.2407009817987995d+0*t25*t2
     9         6-1.2407009817987995d+0*t26)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.481426920246762d-1*t1*t2
     1         7*(1.2437810945273631d+0*t4*(1.6455307846020564d-2*gammaa
     2         a*t2*(4.97175d-2*gammaaa*t14*t31*t32-1.553671875d-3*t10*t
     3         33*t34*t9)-9.012345679012345d-2*t22*t38*(2.70777156307306
     4         d-4*gammaaa*t8+5.625d-3*gammaaa*t31*t32)+1.02663468326446
     5         56d-5*gammaaa*t8+4.509905755707676d-6*t6*t7+1.44197530864
     6         19754d-1*t22*t40-1.802469135802469d-1*t16*t40+3.443772335
     7         859989d-3*gammaaa*t31*t32+1.6455307846020564d-2*t15*t2)-5
     8         .074776173568845d-2*t2*t23*t39)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t27*
     1         t4*(1.6455307846020564d-2*gammaaa*t2*(1.553671875d-3*gamm
     2         aaa*t31*t33*t34*t9-4.97175d-2*t14*t32*t41*t7)+5.069444444
     3         444443d-4*t22*t32*t38*t41*t7-3.443772335859989d-3*t32*t41
     4         *t7+1.4419753086419754d-1*t22*t42-1.802469135802469d-1*t1
     5         6*t42)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = 2.040060021774676d-2*gammabb*t2+1.0d+0
            t4 = 1/t3**2.0d+0
            t5 = gammabb**3
            t6 = 1/rhob**8
            t7 = gammabb**2
            t8 = 1/rhob**5.333333333333333d+0
            t9 = 1/rhob
            t10 = 1/taub
            t11 = gammabb*t10*t9
            t12 = t11**2.0d+0
            t13 = 1.5625d-2*t12+1.0d+0
            t14 = 1/t13**2.0d+0
            t15 = 2.485875d-2*t12*t14+1.234567901234568d-1
            t16 = (1.35388578153653d-4*t7*t8+2.8125d-3*t12)**5.0d-1
            t17 = 1/rhob**1.6666666666666669d+0
            t18 = 2.0d+0*taub-2.5d-1*gammabb*t9
            t19 = 1.0970205230680374d-1*t17*t18-1.0d+0
            t20 = 4.38808209227215d-2*t17*t18*t19+1.0d+0
            t21 = 1/t20**5.0d-1
            t22 = 4.5d-1*t19*t21+1.0970205230680376d-2*gammabb*t2
            t23 = 5.133173416322327d-6*t7*t8+1.5033019185692253d-6*t5*t6
     1         +7.209876543209877d-2*t22**2-1.802469135802469d-1*t16*t22
     2         +1.6455307846020564d-2*gammabb*t15*t2+1.7218861679299946d
     3         -3*t12
            t24 = 1.2437810945273631d+0*t23*t4+1.0d+0
            t25 = 8.04d-1-8.04d-1/t24
            t26 = rhob**3.333333333333333d-1
            t27 = 1/t24**2
            t28 = 1/rhob**6.333333333333333d+0
            t29 = 1/rhob**3.6666666666666664d+0
            t30 = 1/rhob**3
            t31 = 1/taub**2
            t32 = 1/rhob**2
            t33 = 1/t13**3.0d+0
            t34 = t11**3.0d+0
            t35 = 2.7425513076700936d-2*gammabb*t29-1.828367538446729d-1
     1         *t18*t2
            t36 = 1/t20**1.5d+0
            t37 = -2.25d-1*t19*(4.38808209227215d-2*t17*t18*t35+1.097020
     1         5230680374d-2*gammabb*t19*t29-7.313470153786916d-2*t18*t1
     2         9*t2)*t36+4.5d-1*t21*t35-2.925388061514767d-2*gammabb*t29
            t38 = 1/t16
            t39 = 1/t3**3.0d+0
            t40 = -2.25d-1*t19*(-1.0970205230680374d-2*t19*t2-1.20345402
     1         80324706d-3*t18/rhob**4.333333333333333d+0)*t36-1.2341480
     2         884515421d-2*t2*t21+1.0970205230680376d-2*t2
            t41 = 1/taub**3
            t42 = 9.873184707612337d-2*t17*t21-2.25d-1*t19*(8.7761641845
     1         443d-2*t17*t19+9.627632224259765d-3*t18/rhob**3.333333333
     2         3333337d+0)*t36
            fnc(iq) = (-9.305257363490997d-1*t1*t25-9.305257363490997d-1
     1         *t1)*wght+fnc(iq)
            Amat(iq,D1_RB) = (-7.481426920246762d-1*t1*t27*(1.2437810945
     1         273631d+0*t4*(1.6455307846020564d-2*gammabb*t2*(1.5536718
     2         75d-3*gammabb*t10*t32*t33*t34-4.97175d-2*t14*t30*t31*t7)-
     3         9.012345679012345d-2*t22*t38*(-5.625d-3*t30*t31*t7-7.2207
     4         24168194825d-4*t28*t7)-3.443772335859989d-3*t30*t31*t7-2.
     5         737692488705241d-5*t28*t7-1.2026415348553804d-5*t5/rhob**
     6         9+1.4419753086419754d-1*t22*t37-1.802469135802469d-1*t16*
     7         t37-4.3880820922721503d-2*gammabb*t15*t29)+1.353273646285
     8         0253d-1*gammabb*t23*t29*t39)-1.2407009817987995d+0*t25*t2
     9         6-1.2407009817987995d+0*t26)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.481426920246762d-1*t1*t2
     1         7*(1.2437810945273631d+0*t4*(1.6455307846020564d-2*gammab
     2         b*t2*(4.97175d-2*gammabb*t14*t31*t32-1.553671875d-3*t10*t
     3         33*t34*t9)-9.012345679012345d-2*t22*t38*(2.70777156307306
     4         d-4*gammabb*t8+5.625d-3*gammabb*t31*t32)+1.02663468326446
     5         56d-5*gammabb*t8+4.509905755707676d-6*t6*t7+1.44197530864
     6         19754d-1*t22*t40-1.802469135802469d-1*t16*t40+3.443772335
     7         859989d-3*gammabb*t31*t32+1.6455307846020564d-2*t15*t2)-5
     8         .074776173568845d-2*t2*t23*t39)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t1*t27*
     1         t4*(1.6455307846020564d-2*gammabb*t2*(1.553671875d-3*gamm
     2         abb*t31*t33*t34*t9-4.97175d-2*t14*t32*t41*t7)+5.069444444
     3         444443d-4*t22*t32*t38*t41*t7-3.443772335859989d-3*t32*t41
     4         *t7+1.4419753086419754d-1*t22*t42-1.802469135802469d-1*t1
     5         6*t42)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_tpss03 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sigma_{\alpha\alpha}^2\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& \left({{{\it t_6}\,\sigma_{\alpha\alpha}}\over{\tau_\alpha}}\right)^{2.0}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 2.0\,\tau_\alpha-0.25\,{\it t_6}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{10}} &=& .1097020523068037\,{\it t_8}\,{\it t_9}-1.0\\\\ 
C>   {\it t_{11}} &=& {{0.45\,{\it t_{10}}}\over{\left(0.0438808209227215\,{\it t_8}\,{\it t_9}\,{\it t_{10}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_3}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\beta^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \sigma_{\beta\beta}^2\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{16}} &=& \left({{{\it t_{15}}\,\sigma_{\beta\beta}}\over{\tau_\beta}}\right)^{2.0}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 2.0\,\tau_\beta-0.25\,{\it t_{15}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{19}} &=& .1097020523068037\,{\it t_{17}}\,{\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {{0.45\,{\it t_{19}}}\over{\left(0.0438808209227215\,{\it t_{17}}\,{\it t_{18}}\,{\it t_{19}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{12}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^2\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{26}} &=& \left({{{\it t_{25}}\,\sigma_{ss}}\over{\tau_s}}\right)^{2.0}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& 2.0\,\tau_s-0.25\,{\it t_{25}}\,\sigma_{ss}\\\\ 
C>   {\it t_{29}} &=& .1097020523068037\,{\it t_{27}}\,{\it t_{28}}-1.0\\\\ 
C>   {\it t_{30}} &=& {{0.45\,{\it t_{29}}}\over{\left(0.0438808209227215\,{\it t_{27}}\,{\it t_{28}}\,{\it t_{29}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{22}}\,\sigma_{ss}\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{20}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{16}}+1.35388578153653 \times 10^{-4}\,{\it t_{13}}\,{\it t_{14}}\right)^{0.5}\,{\it t_{20}}+.01645530784602056\,{\it t_{12}}\,\sigma_{\beta\beta}\,\left({{0.02485875\,{\it t_{16}}}\over{\left(0.015625\,{\it t_{16}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{16}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{13}}\,{\it t_{14}}\right)}\over{\left(.02040060021774676\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{11}}^2-.1802469135802469\,\left(0.0028125\,{\it t_7}+1.35388578153653 \times 10^{-4}\,{\it t_4}\,{\it t_5}\right)^{0.5}\,{\it t_{11}}+.01645530784602056\,{\it t_3}\,\sigma_{\alpha\alpha}\,\left({{0.02485875\,{\it t_7}}\over{\left(0.015625\,{\it t_7}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_7}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8}}+5.133173416322327 \times 10^{-6}\,{\it t_4}\,{\it t_5}\right)}\over{\left(.02040060021774676\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_{21}}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{30}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{26}}+1.35388578153653 \times 10^{-4}\,{\it t_{23}}\,{\it t_{24}}\right)^{0.5}\,{\it t_{30}}+.01645530784602056\,{\it t_{22}}\,\sigma_{ss}\,\left({{0.02485875\,{\it t_{26}}}\over{\left(0.015625\,{\it t_{26}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{26}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{ss}^3}\over{\rho_s^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}\right)}\over{\left(.02040060021774676\,{\it t_{22}}\,\sigma_{ss}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_{21}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] J Tao, JP Perdew, VN Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, 146401 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.91.146401 "> 
C> 10.1103/PhysRevLett.91.146401 </a> 
C> 
      subroutine nwxcm_x_tpss03_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
                                             !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_tpss03 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sigma_{\alpha\alpha}^2\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& \left({{{\it t_6}\,\sigma_{\alpha\alpha}}\over{\tau_\alpha}}\right)^{2.0}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 2.0\,\tau_\alpha-0.25\,{\it t_6}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{10}} &=& .1097020523068037\,{\it t_8}\,{\it t_9}-1.0\\\\ 
C>   {\it t_{11}} &=& {{0.45\,{\it t_{10}}}\over{\left(0.0438808209227215\,{\it t_8}\,{\it t_9}\,{\it t_{10}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_3}\,\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\beta^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \sigma_{\beta\beta}^2\\\\ 
C>   {\it t_{15}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{16}} &=& \left({{{\it t_{15}}\,\sigma_{\beta\beta}}\over{\tau_\beta}}\right)^{2.0}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& 2.0\,\tau_\beta-0.25\,{\it t_{15}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{19}} &=& .1097020523068037\,{\it t_{17}}\,{\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {{0.45\,{\it t_{19}}}\over{\left(0.0438808209227215\,{\it t_{17}}\,{\it t_{18}}\,{\it t_{19}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{12}}\,\sigma_{\beta\beta}\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s^{{{16}\over{3}}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^2\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{26}} &=& \left({{{\it t_{25}}\,\sigma_{ss}}\over{\tau_s}}\right)^{2.0}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{\rho_s^{{{5}\over{3}}}}}\\\\ 
C>   {\it t_{28}} &=& 2.0\,\tau_s-0.25\,{\it t_{25}}\,\sigma_{ss}\\\\ 
C>   {\it t_{29}} &=& .1097020523068037\,{\it t_{27}}\,{\it t_{28}}-1.0\\\\ 
C>   {\it t_{30}} &=& {{0.45\,{\it t_{29}}}\over{\left(0.0438808209227215\,{\it t_{27}}\,{\it t_{28}}\,{\it t_{29}}+1.0\right)^{0.5}}}+.01097020523068038\,{\it t_{22}}\,\sigma_{ss}\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{20}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{16}}+1.35388578153653 \times 10^{-4}\,{\it t_{13}}\,{\it t_{14}}\right)^{0.5}\,{\it t_{20}}+.01645530784602056\,{\it t_{12}}\,\sigma_{\beta\beta}\,\left({{0.02485875\,{\it t_{16}}}\over{\left(0.015625\,{\it t_{16}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{16}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{13}}\,{\it t_{14}}\right)}\over{\left(.02040060021774676\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{11}}^2-.1802469135802469\,\left(0.0028125\,{\it t_7}+1.35388578153653 \times 10^{-4}\,{\it t_4}\,{\it t_5}\right)^{0.5}\,{\it t_{11}}+.01645530784602056\,{\it t_3}\,\sigma_{\alpha\alpha}\,\left({{0.02485875\,{\it t_7}}\over{\left(0.015625\,{\it t_7}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_7}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8}}+5.133173416322327 \times 10^{-6}\,{\it t_4}\,{\it t_5}\right)}\over{\left(.02040060021774676\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_{21}}\,\left(0.804-{{0.804}\over{{{1.243781094527363\,\left(.07209876543209877\,{\it t_{30}}^2-.1802469135802469\,\left(0.0028125\,{\it t_{26}}+1.35388578153653 \times 10^{-4}\,{\it t_{23}}\,{\it t_{24}}\right)^{0.5}\,{\it t_{30}}+.01645530784602056\,{\it t_{22}}\,\sigma_{ss}\,\left({{0.02485875\,{\it t_{26}}}\over{\left(0.015625\,{\it t_{26}}+1.0\right)^{2.0}}}+.1234567901234568\right)+.001721886167929995\,{\it t_{26}}+{{1.5033019185692253 \times 10^{-6}\,\sigma_{ss}^3}\over{\rho_s^8}}+5.133173416322327 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}\right)}\over{\left(.02040060021774676\,{\it t_{22}}\,\sigma_{ss}+1.0\right)^{2.0}}}+1.0}}\right)-.9305257363490997\,{\it t_{21}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] J Tao, JP Perdew, VN Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, 146401 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.91.146401 "> 
C> 10.1103/PhysRevLett.91.146401 </a> 
C> 
      subroutine nwxcm_x_tpss03_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
                                             !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat3(nq,NCOL_MMAT3)  !< [Output] The 3rd derivative wrt tau 
                                             !< and possibly rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
