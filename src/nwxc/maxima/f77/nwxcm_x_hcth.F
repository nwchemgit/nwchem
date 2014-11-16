C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_hcth.F 
C> The nwxcm_x_hcth functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_hcth functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{10}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{11}} &=& 0.004\,{\it t_{10}}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,\sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8\,{\it t_9}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{\it t_9}^4}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}\,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{\it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8\,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{\it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
      subroutine nwxcm_x_hcth(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Cmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = param(6)
            t2 = gammaaa**4
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**4
            t6 = 1/rhoa**1.0666666666666666d+1
            t7 = param(5)
            t8 = gammaaa**3
            t9 = 1/t4**3
            t10 = 1/rhoa**8
            t11 = param(4)
            t12 = gammaaa**2
            t13 = 1/t4**2
            t14 = 1/rhoa**5.333333333333333d+0
            t15 = param(3)
            t16 = 1/t4
            t17 = 6.4d-8*t10*t7*t8*t9+2.5600000000000006d-10*t1*t2*t5*t6
     1         +4.0d-3*gammaaa*t15*t16*t3+1.6d-5*t11*t12*t13*t14+param(2
     2         )
            t18 = rhoa**1.3333333333333333d+0
            t19 = 1/t4**5
            t20 = 1/rhoa**1.1666666666666665d+1
            t21 = 1/rhoa**9
            t22 = 1/rhoa**6.333333333333333d+0
            t23 = -1.8610514726981986d+0*t18*(t10*(1.9200000000000003d-7
     1         *t12*t7-1.28d-7*t11*t12)*t9+t5*t6*(1.0240000000000002d-9*
     2         t1*t8-7.68d-10*t7*t8)+4.0d-3*t15*t16*t3-4.096000000000001
     3         3d-12*t1*t19*t2/rhoa**1.3333333333333333d+1+t13*t14*(3.2d
     4         -5*gammaaa*t11-1.6d-5*gammaaa*t15))*wght
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t17*t18*wght
            Amat(iq,D1_RA) = -1.8610514726981986d+0*t18*(-5.12d-7*t21*t7
     1         *t8*t9+3.413333333333333d-7*t11*t21*t8*t9+2.0480000000000
     2         004d-9*t2*t20*t5*t7-2.7306666666666674d-9*t1*t2*t20*t5+4.
     3         266666666666667d-5*t12*t13*t15*t22-8.533333333333333d-5*t
     4         11*t12*t13*t22+1.0922666666666668d-11*gammaaa**5*t1*t19/r
     5         hoa**1.4333333333333334d+1-1.0666666666666666d-2*gammaaa*
     6         t15*t16/rhoa**3.6666666666666664d+0)*wght-2.4814019635975
     7         98d+0*rhoa**3.333333333333333d-1*t17*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t23+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t23+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)
            t2 = param(6)
            t3 = gammaaa**4
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 4.0d-3*gammaaa*t4+1.0d+0
            t6 = 1/t5**4
            t7 = 1/rhoa**1.0666666666666666d+1
            t8 = param(5)
            t9 = gammaaa**3
            t10 = 1/t5**3
            t11 = 1/rhoa**8
            t12 = param(4)
            t13 = gammaaa**2
            t14 = 1/t5**2
            t15 = 1/rhoa**5.333333333333333d+0
            t16 = param(3)
            t17 = 1/t5
            t18 = 6.4d-8*t10*t11*t8*t9+2.5600000000000006d-10*t2*t3*t6*t
     1         7+4.0d-3*gammaaa*t16*t17*t4+1.6d-5*t12*t13*t14*t15+t1
            t19 = rhoa**1.3333333333333333d+0
            t20 = gammabb**4
            t21 = 1/rhob**2.6666666666666666d+0
            t22 = 4.0d-3*gammabb*t21+1.0d+0
            t23 = 1/t22**4
            t24 = 1/rhob**1.0666666666666666d+1
            t25 = gammabb**3
            t26 = 1/t22**3
            t27 = 1/rhob**8
            t28 = gammabb**2
            t29 = 1/t22**2
            t30 = 1/rhob**5.333333333333333d+0
            t31 = 1/t22
            t32 = 6.4d-8*t25*t26*t27*t8+4.0d-3*gammabb*t16*t21*t31+1.6d-
     1         5*t12*t28*t29*t30+2.5600000000000006d-10*t2*t20*t23*t24+t
     2         1
            t33 = rhob**1.3333333333333333d+0
            t34 = 1/t5**5
            t35 = 1/rhoa**1.1666666666666665d+1
            t36 = 1/rhoa**9
            t37 = 1/rhoa**6.333333333333333d+0
            t38 = 1/t22**5
            t39 = 1/rhob**1.1666666666666665d+1
            t40 = 1/rhob**9
            t41 = 1/rhob**6.333333333333333d+0
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t32*t33-9.3052573634
     1         90993d-1*t18*t19)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t19*(-5.12d-7
     1         *t10*t36*t8*t9+3.413333333333333d-7*t10*t12*t36*t9+2.0480
     2         000000000004d-9*t3*t35*t6*t8-2.7306666666666674d-9*t2*t3*
     3         t35*t6+4.266666666666667d-5*t13*t14*t16*t37-8.53333333333
     4         3333d-5*t12*t13*t14*t37+1.0922666666666668d-11*gammaaa**5
     5         *t2*t34/rhoa**1.4333333333333334d+1-1.0666666666666666d-2
     6         *gammaaa*t16*t17/rhoa**3.6666666666666664d+0)-1.240700981
     7         798799d+0*rhoa**3.333333333333333d-1*t18)*wght+Amat(iq,D1
     8         _RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t33*(-5.12d-7
     1         *t25*t26*t40*t8+2.0480000000000004d-9*t20*t23*t39*t8+4.26
     2         6666666666667d-5*t16*t28*t29*t41-8.533333333333333d-5*t12
     3         *t28*t29*t41+3.413333333333333d-7*t12*t25*t26*t40-2.73066
     4         66666666674d-9*t2*t20*t23*t39+1.0922666666666668d-11*gamm
     5         abb**5*t2*t38/rhob**1.4333333333333334d+1-1.0666666666666
     6         666d-2*gammabb*t16*t31/rhob**3.6666666666666664d+0)-1.240
     7         700981798799d+0*rhob**3.333333333333333d-1*t32)*wght+Amat
     8         (iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t19*(
     1         t6*t7*(1.0240000000000002d-9*t2*t9-7.68d-10*t8*t9)+t10*t1
     2         1*(1.9200000000000003d-7*t13*t8-1.28d-7*t12*t13)+4.0d-3*t
     3         16*t17*t4-4.0960000000000013d-12*t2*t3*t34/rhoa**1.333333
     4         3333333333d+1+t14*t15*(3.2d-5*gammaaa*t12-1.6d-5*gammaaa*
     5         t16))*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t33*(
     1         t26*t27*(1.9200000000000003d-7*t28*t8-1.28d-7*t12*t28)+t2
     2         3*t24*(1.0240000000000002d-9*t2*t25-7.68d-10*t25*t8)-4.09
     3         60000000000013d-12*t2*t20*t38/rhob**1.3333333333333333d+1
     4         +4.0d-3*t16*t21*t31+(3.2d-5*gammabb*t12-1.6d-5*gammabb*t1
     5         6)*t29*t30)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(6)
            t2 = gammaaa**4
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**4
            t6 = 1/rhoa**1.0666666666666666d+1
            t7 = param(5)
            t8 = gammaaa**3
            t9 = 1/t4**3
            t10 = 1/rhoa**8
            t11 = param(4)
            t12 = gammaaa**2
            t13 = 1/t4**2
            t14 = 1/rhoa**5.333333333333333d+0
            t15 = param(3)
            t16 = 1/t4
            t17 = 6.4d-8*t10*t7*t8*t9+2.5600000000000006d-10*t1*t2*t5*t6
     1         +4.0d-3*gammaaa*t15*t16*t3+1.6d-5*t11*t12*t13*t14+param(2
     2         )
            t18 = rhoa**1.3333333333333333d+0
            t19 = 1/t4**5
            t20 = 1/rhoa**1.1666666666666665d+1
            t21 = 1/rhoa**9
            t22 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RA) = -9.305257363490993d-1*t18*(-5.12d-7*t21*t7*
     1         t8*t9+3.413333333333333d-7*t11*t21*t8*t9+2.04800000000000
     2         04d-9*t2*t20*t5*t7-2.7306666666666674d-9*t1*t2*t20*t5+4.2
     3         66666666666667d-5*t12*t13*t15*t22-8.533333333333333d-5*t1
     4         1*t12*t13*t22+1.0922666666666668d-11*gammaaa**5*t1*t19/rh
     5         oa**1.4333333333333334d+1-1.0666666666666666d-2*gammaaa*t
     6         15*t16/rhoa**3.6666666666666664d+0)*wght-1.24070098179879
     7         9d+0*rhoa**3.333333333333333d-1*t17*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*(1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12)*t9+t5*
     2         t6*(1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8)+4.0d-3*t1
     3         5*t16*t3-4.0960000000000013d-12*t1*t19*t2/rhoa**1.3333333
     4         333333333d+1+t13*t14*(3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t
     5         15))*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(6)
            t2 = gammabb**4
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 4.0d-3*gammabb*t3+1.0d+0
            t5 = 1/t4**4
            t6 = 1/rhob**1.0666666666666666d+1
            t7 = param(5)
            t8 = gammabb**3
            t9 = 1/t4**3
            t10 = 1/rhob**8
            t11 = param(4)
            t12 = gammabb**2
            t13 = 1/t4**2
            t14 = 1/rhob**5.333333333333333d+0
            t15 = param(3)
            t16 = 1/t4
            t17 = 6.4d-8*t10*t7*t8*t9+2.5600000000000006d-10*t1*t2*t5*t6
     1         +4.0d-3*gammabb*t15*t16*t3+1.6d-5*t11*t12*t13*t14+param(2
     2         )
            t18 = rhob**1.3333333333333333d+0
            t19 = 1/t4**5
            t20 = 1/rhob**1.1666666666666665d+1
            t21 = 1/rhob**9
            t22 = 1/rhob**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RB) = -9.305257363490993d-1*t18*(-5.12d-7*t21*t7*
     1         t8*t9+3.413333333333333d-7*t11*t21*t8*t9+2.04800000000000
     2         04d-9*t2*t20*t5*t7-2.7306666666666674d-9*t1*t2*t20*t5+4.2
     3         66666666666667d-5*t12*t13*t15*t22-8.533333333333333d-5*t1
     4         1*t12*t13*t22+1.0922666666666668d-11*gammabb**5*t1*t19/rh
     5         ob**1.4333333333333334d+1-1.0666666666666666d-2*gammabb*t
     6         15*t16/rhob**3.6666666666666664d+0)*wght-1.24070098179879
     7         9d+0*rhob**3.333333333333333d-1*t17*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t18*(
     1         t10*(1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12)*t9+t5*
     2         t6*(1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8)+4.0d-3*t1
     3         5*t16*t3-4.0960000000000013d-12*t1*t19*t2/rhob**1.3333333
     4         333333333d+1+t13*t14*(3.2d-5*gammabb*t11-1.6d-5*gammabb*t
     5         15))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_hcth functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{10}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{11}} &=& 0.004\,{\it t_{10}}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,\sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8\,{\it t_9}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{\it t_9}^4}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}\,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{\it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8\,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{\it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
      subroutine nwxcm_x_hcth_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_hcth functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{10}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{11}} &=& 0.004\,{\it t_{10}}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,\sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}\over{\rho_\beta^8\,{\it t_9}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{\it t_9}^4}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}\,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{\it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8\,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{\it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
      subroutine nwxcm_x_hcth_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
