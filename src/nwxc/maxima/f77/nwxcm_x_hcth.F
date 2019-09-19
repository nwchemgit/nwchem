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
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,
C>    \sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_6}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}
C>    \over{\rho_\beta^8\,{\it t_9}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{
C>    \it t_9}^4}}+{\it t_1}\right)-0.9305257363490993
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}
C>    +{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}
C>    \,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}
C>    \,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{
C>    \it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8
C>    \,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{
C>    \it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{
C>    \it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109
C>    , 6264 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_hcth(param,tol_rho,ipol,nq,wght, 
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
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t17*t18*wght
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t18*(-5.12d-7
     1         *t21*t7*t8*t9+3.413333333333333d-7*t11*t21*t8*t9+2.048000
     2         0000000004d-9*t2*t20*t5*t7-2.7306666666666674d-9*t1*t2*t2
     3         0*t5+4.266666666666667d-5*t12*t13*t15*t22-8.5333333333333
     4         33d-5*t11*t12*t13*t22+1.0922666666666668d-11*gammaaa**5*t
     5         1*t19/rhoa**1.4333333333333334d+1-1.0666666666666666d-2*g
     6         ammaaa*t15*t16/rhoa**3.6666666666666664d+0)-1.24070098179
     7         8799d+0*rhoa**3.333333333333333d-1*t17)*wght+Amat(iq,D1_R
     8         A)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*(1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12)*t9+t5*
     2         t6*(1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8)+4.0d-3*t1
     3         5*t16*t3-4.0960000000000013d-12*t1*t19*t2/rhoa**1.3333333
     4         333333333d+1+t13*t14*(3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t
     5         15))*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
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
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,
C>    \sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_6}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}
C>    \over{\rho_\beta^8\,{\it t_9}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{
C>    \it t_9}^4}}+{\it t_1}\right)-0.9305257363490993
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}
C>    +{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}
C>    \,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}
C>    \,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{
C>    \it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8
C>    \,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{
C>    \it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{
C>    \it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109
C>    , 6264 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_hcth_d2(param,tol_rho,ipol,nq,wght, 
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
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**5
            t21 = 1/t4**5
            t22 = 1/rhoa**1.4333333333333334d+1
            t23 = 1/rhoa**1.1666666666666665d+1
            t24 = 1/rhoa**9
            t25 = 1/rhoa**6.333333333333333d+0
            t26 = 1/rhoa**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammaaa*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhoa**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/t4**6
            t35 = 1/rhoa**1.5333333333333335d+1
            t36 = 1/rhoa**1.2666666666666665d+1
            t37 = 1/rhoa**10
            t38 = 1/rhoa**7.333333333333333d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t17*t18*wght
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t18*t27-1.240
     1         700981798799d+0*t17*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t18*(4.60
     1         80000000000004d-6*t37*t7*t8*t9+9.102222222222223d-7*t15*t
     2         37*t8*t9-4.8924444444444437d-6*t11*t37*t8*t9-4.0277333333
     3         333337d-8*t2*t36*t5*t7+8.738133333333334d-11*t20*t21*t35*
     4         t7+1.0922666666666667d-8*t11*t2*t36*t5+3.1857777777777785
     5         d-8*t1*t2*t36*t5-3.84d-4*t12*t13*t15*t38+5.40444444444444
     6         4d-4*t11*t12*t13*t38-2.730666666666667d-10*t1*t20*t21*t35
     7         +5.825422222222222d-13*gammaaa**6*t1*t34/rhoa**18+3.91111
     8         1111111111d-2*gammaaa*t15*t16/rhoa**4.666666666666667d+0)
     9         -2.481401963597598d+0*t19*t27-4.13566993932933d-1*t17/rho
     :         a**6.666666666666666d-1)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t18*(-1.
     1         5360000000000001d-6*t12*t24*t7*t9-3.413333333333334d-7*t1
     2         2*t15*t24*t9+1.7066666666666666d-6*t11*t12*t24*t9+1.43360
     3         00000000002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t11*t23
     4         *t5*t8-1.0922666666666668d-8*t1*t23*t5*t8-3.2768000000000
     5         01d-11*t2*t21*t22*t7-2.1845333333333336d-13*t1*t20*t34/rh
     6         oa**17-1.0666666666666666d-2*t15*t16*t26+1.28000000000000
     7         03d-4*gammaaa*t13*t15*t25-1.7066666666666666d-4*gammaaa*t
     8         11*t13*t25+9.830400000000002d-11*t1*t2*t21*t22)-1.2407009
     9         81798799d+0*t19*(1.9200000000000003d-7*t10*t12*t7*t9-1.28
     :         d-7*t10*t11*t12*t9-7.68d-10*t5*t6*t7*t8+1.024000000000000
     ;         2d-9*t1*t5*t6*t8+t33+t29-1.6d-5*gammaaa*t13*t14*t15+3.2d-
     <         5*gammaaa*t11*t13*t14))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t18*(t10*(3.8400000000000006d-7*gammaaa*t7-8.0d-3*t3
     2         2-2.56d-7*gammaaa*t11)*t9+t21*t28*(-1.6384000000000004d-1
     3         1*t1*t8-1.6d-2*t30)+t5*t6*(-2.3040000000000002d-9*t12*t7-
     4         1.2d-2*t31+3.0720000000000003d-9*t1*t12)+8.19200000000000
     5         2d-14*t1*t2*t34/rhoa**16+t13*t14*(3.2d-5*t11-3.2d-5*t15))
     6         *wght
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
            t34 = rhoa**3.333333333333333d-1
            t35 = gammaaa**5
            t36 = 1/t5**5
            t37 = 1/rhoa**1.4333333333333334d+1
            t38 = 1/rhoa**1.1666666666666665d+1
            t39 = 1/rhoa**9
            t40 = 1/rhoa**6.333333333333333d+0
            t41 = 1/rhoa**3.6666666666666664d+0
            t42 = -5.12d-7*t10*t39*t8*t9+3.413333333333333d-7*t10*t12*t3
     1         9*t9+2.0480000000000004d-9*t3*t38*t6*t8-2.730666666666667
     2         4d-9*t2*t3*t38*t6-1.0666666666666666d-2*gammaaa*t16*t17*t
     3         41+4.266666666666667d-5*t13*t14*t16*t40-8.533333333333333
     4         d-5*t12*t13*t14*t40+1.0922666666666668d-11*t2*t35*t36*t37
            t43 = rhob**3.333333333333333d-1
            t44 = gammabb**5
            t45 = 1/t22**5
            t46 = 1/rhob**1.4333333333333334d+1
            t47 = 1/rhob**1.1666666666666665d+1
            t48 = 1/rhob**9
            t49 = 1/rhob**6.333333333333333d+0
            t50 = 1/rhob**3.6666666666666664d+0
            t51 = -5.12d-7*t25*t26*t48*t8+2.0480000000000004d-9*t20*t23*
     1         t47*t8-1.0666666666666666d-2*gammabb*t16*t31*t50+4.266666
     2         666666667d-5*t16*t28*t29*t49-8.533333333333333d-5*t12*t28
     3         *t29*t49+3.413333333333333d-7*t12*t25*t26*t48-2.730666666
     4         6666674d-9*t2*t20*t23*t47+1.0922666666666668d-11*t2*t44*t
     5         45*t46
            t52 = 1/rhoa**1.3333333333333333d+1
            t53 = -4.0960000000000013d-12*t2*t3*t36*t52
            t54 = 1.0240000000000002d-9*t2*t9-7.68d-10*t8*t9
            t55 = 1.9200000000000003d-7*t13*t8-1.28d-7*t12*t13
            t56 = 3.2d-5*gammaaa*t12-1.6d-5*gammaaa*t16
            t57 = 4.0d-3*t16*t17*t4
            t58 = 1/rhob**1.3333333333333333d+1
            t59 = -4.0960000000000013d-12*t2*t20*t45*t58
            t60 = 1.0240000000000002d-9*t2*t25-7.68d-10*t25*t8
            t61 = 1.9200000000000003d-7*t28*t8-1.28d-7*t12*t28
            t62 = 3.2d-5*gammabb*t12-1.6d-5*gammabb*t16
            t63 = 4.0d-3*t16*t21*t31
            t64 = 1/t5**6
            t65 = 1/rhoa**1.5333333333333335d+1
            t66 = 1/rhoa**1.2666666666666665d+1
            t67 = 1/rhoa**10
            t68 = 1/rhoa**7.333333333333333d+0
            t69 = 1/t22**6
            t70 = 1/rhob**1.5333333333333335d+1
            t71 = 1/rhob**1.2666666666666665d+1
            t72 = 1/rhob**10
            t73 = 1/rhob**7.333333333333333d+0
            t74 = 3.2d-5*t12-3.2d-5*t16
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t32*t33-9.3052573634
     1         90993d-1*t18*t19)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t19*t42-1.240
     1         700981798799d+0*t18*t34)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t33*t51-1.240
     1         700981798799d+0*t32*t43)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t19*(
     1         t54*t6*t7+t57+t14*t15*t56+t10*t11*t55+t53)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t33*(
     1         t63+t29*t30*t62+t26*t27*t61+t23*t24*t60+t59)*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t19*(4.60
     1         80000000000004d-6*t10*t67*t8*t9+9.102222222222223d-7*t10*
     2         t16*t67*t9-4.8924444444444437d-6*t10*t12*t67*t9-4.0277333
     3         333333337d-8*t3*t6*t66*t8+8.738133333333334d-11*t35*t36*t
     4         65*t8-3.84d-4*t13*t14*t16*t68+5.404444444444444d-4*t12*t1
     5         3*t14*t68+3.1857777777777785d-8*t2*t3*t6*t66+1.0922666666
     6         666667d-8*t12*t3*t6*t66-2.730666666666667d-10*t2*t35*t36*
     7         t65+5.825422222222222d-13*gammaaa**6*t2*t64/rhoa**18+3.91
     8         1111111111111d-2*gammaaa*t16*t17/rhoa**4.666666666666667d
     9         +0)-2.481401963597598d+0*t34*t42-4.13566993932933d-1*t18/
     :         rhoa**6.666666666666666d-1)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363490993d-1*t33*(4.60
     1         80000000000004d-6*t25*t26*t72*t8-4.0277333333333337d-8*t2
     2         0*t23*t71*t8+8.738133333333334d-11*t44*t45*t70*t8-3.84d-4
     3         *t16*t28*t29*t73+5.404444444444444d-4*t12*t28*t29*t73+9.1
     4         02222222222223d-7*t16*t25*t26*t72-4.8924444444444437d-6*t
     5         12*t25*t26*t72+3.1857777777777785d-8*t2*t20*t23*t71+1.092
     6         2666666666667d-8*t12*t20*t23*t71-2.730666666666667d-10*t2
     7         *t44*t45*t70+5.825422222222222d-13*gammabb**6*t2*t69/rhob
     8         **18+3.911111111111111d-2*gammabb*t16*t31/rhob**4.6666666
     9         66666667d+0)-2.481401963597598d+0*t43*t51-4.1356699393293
     :         3d-1*t32/rhob**6.666666666666666d-1)*wght+Amat2(iq,D2_RB_
     ;         RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-1.240700981798799d+0*t34*(-7.
     1         68d-10*t6*t7*t8*t9+1.0240000000000002d-9*t2*t6*t7*t9+1.92
     2         00000000000003d-7*t10*t11*t13*t8+t57+t53-1.6d-5*gammaaa*t
     3         14*t15*t16+3.2d-5*gammaaa*t12*t14*t15-1.28d-7*t10*t11*t12
     4         *t13)-9.305257363490993d-1*t19*(1.4336000000000002d-8*t38
     5         *t6*t8*t9-1.0922666666666668d-8*t2*t38*t6*t9-4.0959999999
     6         999996d-9*t12*t38*t6*t9-1.5360000000000001d-6*t10*t13*t39
     7         *t8-3.276800000000001d-11*t3*t36*t37*t8-2.184533333333333
     8         6d-13*t2*t35*t64/rhoa**17-1.0666666666666666d-2*t16*t17*t
     9         41+1.2800000000000003d-4*gammaaa*t14*t16*t40-1.7066666666
     :         666666d-4*gammaaa*t12*t14*t40-3.413333333333334d-7*t10*t1
     ;         3*t16*t39+1.7066666666666666d-6*t10*t12*t13*t39+9.8304000
     <         00000002d-11*t2*t3*t36*t37))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-9.305257363490993d-1*t33*(-1.
     1         5360000000000001d-6*t26*t28*t48*t8+1.4336000000000002d-8*
     2         t23*t25*t47*t8-3.276800000000001d-11*t20*t45*t46*t8-2.184
     3         5333333333336d-13*t2*t44*t69/rhob**17-1.0666666666666666d
     4         -2*t16*t31*t50+1.2800000000000003d-4*gammabb*t16*t29*t49-
     5         1.7066666666666666d-4*gammabb*t12*t29*t49-3.4133333333333
     6         34d-7*t16*t26*t28*t48+1.7066666666666666d-6*t12*t26*t28*t
     7         48-1.0922666666666668d-8*t2*t23*t25*t47-4.095999999999999
     8         6d-9*t12*t23*t25*t47+9.830400000000002d-11*t2*t20*t45*t46
     9         )-1.240700981798799d+0*t43*(1.9200000000000003d-7*t26*t27
     :         *t28*t8-7.68d-10*t23*t24*t25*t8+t63+t59-1.6d-5*gammabb*t1
     ;         6*t29*t30+3.2d-5*gammabb*t12*t29*t30-1.28d-7*t12*t26*t27*
     <         t28+1.0240000000000002d-9*t2*t23*t24*t25))*wght+Cmat2(iq,
     =         D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t19*(t36*t52*(-1.6384000000000004d-11*t2*t9-1.6d-2*t
     2         54)+t6*t7*(-2.3040000000000002d-9*t13*t8-1.2d-2*t55+3.072
     3         0000000000003d-9*t13*t2)+t10*t11*(3.8400000000000006d-7*g
     4         ammaaa*t8-8.0d-3*t56-2.56d-7*gammaaa*t12)+t14*t15*t74+8.1
     5         92000000000002d-14*t2*t3*t64/rhoa**16)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t33*(t23*t24*(-2.3040000000000002d-9*t28*t8-1.2d-2*t
     2         61+3.0720000000000003d-9*t2*t28)+t26*t27*(3.8400000000000
     3         006d-7*gammabb*t8-8.0d-3*t62-2.56d-7*gammabb*t12)+t29*t30
     4         *t74+8.192000000000002d-14*t2*t20*t69/rhob**16+t45*t58*(-
     5         1.6d-2*t60-1.6384000000000004d-11*t2*t25))*wght
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
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**5
            t21 = 1/t4**5
            t22 = 1/rhoa**1.4333333333333334d+1
            t23 = 1/rhoa**1.1666666666666665d+1
            t24 = 1/rhoa**9
            t25 = 1/rhoa**6.333333333333333d+0
            t26 = 1/rhoa**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammaaa*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhoa**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/t4**6
            t35 = 1/rhoa**1.5333333333333335d+1
            t36 = 1/rhoa**1.2666666666666665d+1
            t37 = 1/rhoa**10
            t38 = 1/rhoa**7.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RA) = -9.305257363490993d-1*t18*t27*wght-1.240700
     1         981798799d+0*t17*t19*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Amat2(iq,D2_RA_RA) = -9.305257363490993d-1*t18*(4.6080000000
     1         000004d-6*t37*t7*t8*t9+9.102222222222223d-7*t15*t37*t8*t9
     2         -4.8924444444444437d-6*t11*t37*t8*t9-4.0277333333333337d-
     3         8*t2*t36*t5*t7+8.738133333333334d-11*t20*t21*t35*t7+1.092
     4         2666666666667d-8*t11*t2*t36*t5+3.1857777777777785d-8*t1*t
     5         2*t36*t5-3.84d-4*t12*t13*t15*t38+5.404444444444444d-4*t11
     6         *t12*t13*t38-2.730666666666667d-10*t1*t20*t21*t35+5.82542
     7         2222222222d-13*gammaaa**6*t1*t34/rhoa**18+3.9111111111111
     8         11d-2*gammaaa*t15*t16/rhoa**4.666666666666667d+0)*wght-2.
     9         481401963597598d+0*t19*t27*wght-4.13566993932933d-1*t17*w
     :         ght/rhoa**6.666666666666666d-1+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490993d-1*t18*(-1.53600000
     1         00000001d-6*t12*t24*t7*t9-3.413333333333334d-7*t12*t15*t2
     2         4*t9+1.7066666666666666d-6*t11*t12*t24*t9+1.4336000000000
     3         002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t11*t23*t5*t8-1
     4         .0922666666666668d-8*t1*t23*t5*t8-3.276800000000001d-11*t
     5         2*t21*t22*t7-2.1845333333333336d-13*t1*t20*t34/rhoa**17-1
     6         .0666666666666666d-2*t15*t16*t26+1.2800000000000003d-4*ga
     7         mmaaa*t13*t15*t25-1.7066666666666666d-4*gammaaa*t11*t13*t
     8         25+9.830400000000002d-11*t1*t2*t21*t22)*wght-1.2407009817
     9         98799d+0*t19*(1.9200000000000003d-7*t10*t12*t7*t9-1.28d-7
     :         *t10*t11*t12*t9-7.68d-10*t5*t6*t7*t8+1.0240000000000002d-
     ;         9*t1*t5*t6*t8+t33+t29-1.6d-5*gammaaa*t13*t14*t15+3.2d-5*g
     <         ammaaa*t11*t13*t14)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t18*(t10*(3.8400000000000006d-7*gammaaa*t7-8.0d-3*t3
     2         2-2.56d-7*gammaaa*t11)*t9+t21*t28*(-1.6384000000000004d-1
     3         1*t1*t8-1.6d-2*t30)+t5*t6*(-2.3040000000000002d-9*t12*t7-
     4         1.2d-2*t31+3.0720000000000003d-9*t1*t12)+8.19200000000000
     5         2d-14*t1*t2*t34/rhoa**16+t13*t14*(3.2d-5*t11-3.2d-5*t15))
     6         *wght
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
            t19 = rhob**3.333333333333333d-1
            t20 = gammabb**5
            t21 = 1/t4**5
            t22 = 1/rhob**1.4333333333333334d+1
            t23 = 1/rhob**1.1666666666666665d+1
            t24 = 1/rhob**9
            t25 = 1/rhob**6.333333333333333d+0
            t26 = 1/rhob**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammabb*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhob**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammabb*t11-1.6d-5*gammabb*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/t4**6
            t35 = 1/rhob**1.5333333333333335d+1
            t36 = 1/rhob**1.2666666666666665d+1
            t37 = 1/rhob**10
            t38 = 1/rhob**7.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RB) = -9.305257363490993d-1*t18*t27*wght-1.240700
     1         981798799d+0*t17*t19*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Amat2(iq,D2_RB_RB) = -9.305257363490993d-1*t18*(4.6080000000
     1         000004d-6*t37*t7*t8*t9+9.102222222222223d-7*t15*t37*t8*t9
     2         -4.8924444444444437d-6*t11*t37*t8*t9-4.0277333333333337d-
     3         8*t2*t36*t5*t7+8.738133333333334d-11*t20*t21*t35*t7+1.092
     4         2666666666667d-8*t11*t2*t36*t5+3.1857777777777785d-8*t1*t
     5         2*t36*t5-3.84d-4*t12*t13*t15*t38+5.404444444444444d-4*t11
     6         *t12*t13*t38-2.730666666666667d-10*t1*t20*t21*t35+5.82542
     7         2222222222d-13*gammabb**6*t1*t34/rhob**18+3.9111111111111
     8         11d-2*gammabb*t15*t16/rhob**4.666666666666667d+0)*wght-2.
     9         481401963597598d+0*t19*t27*wght-4.13566993932933d-1*t17*w
     :         ght/rhob**6.666666666666666d-1+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490993d-1*t18*(-1.53600000
     1         00000001d-6*t12*t24*t7*t9-3.413333333333334d-7*t12*t15*t2
     2         4*t9+1.7066666666666666d-6*t11*t12*t24*t9+1.4336000000000
     3         002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t11*t23*t5*t8-1
     4         .0922666666666668d-8*t1*t23*t5*t8-3.276800000000001d-11*t
     5         2*t21*t22*t7-2.1845333333333336d-13*t1*t20*t34/rhob**17-1
     6         .0666666666666666d-2*t15*t16*t26+1.2800000000000003d-4*ga
     7         mmabb*t13*t15*t25-1.7066666666666666d-4*gammabb*t11*t13*t
     8         25+9.830400000000002d-11*t1*t2*t21*t22)*wght-1.2407009817
     9         98799d+0*t19*(1.9200000000000003d-7*t10*t12*t7*t9-1.28d-7
     :         *t10*t11*t12*t9-7.68d-10*t5*t6*t7*t8+1.0240000000000002d-
     ;         9*t1*t5*t6*t8+t33+t29-1.6d-5*gammabb*t13*t14*t15+3.2d-5*g
     <         ammabb*t11*t13*t14)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t18*(t10*(3.8400000000000006d-7*gammabb*t7-8.0d-3*t3
     2         2-2.56d-7*gammabb*t11)*t9+t21*t28*(-1.6384000000000004d-1
     3         1*t1*t8-1.6d-2*t30)+t5*t6*(-2.3040000000000002d-9*t12*t7-
     4         1.2d-2*t31+3.0720000000000003d-9*t1*t12)+8.19200000000000
     5         2d-14*t1*t2*t34/rhob**16+t13*t14*(3.2d-5*t11-3.2d-5*t15))
     6         *wght
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
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_7}\,{\it t_8}\,
C>    \sigma_{\beta\beta}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_6}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_9}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{\beta\beta}^3}
C>    \over{\rho_\beta^8\,{\it t_9}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\beta\beta}^4}\over{\rho_\beta^{{{32}\over{3}}}\,{
C>    \it t_9}^4}}+{\it t_1}\right)-0.9305257363490993
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_7}\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}
C>    +{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}
C>    \,{\it t_4}^2}}+{{6.4 \times 10^{-8}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}^3}\over{\rho_\alpha^8\,{\it t_4}^3}}
C>    +{{2.560000000000001 \times 10^{-10}\,{\it t_2}\,
C>    \sigma_{\alpha\alpha}^4}\over{\rho_\alpha^{{{32}\over{3}}}
C>    \,{\it t_4}^4}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_7}\,{\it t_{10}}\,\sigma_{ss}}\over{{
C>    \it t_{11}}}}+{{1.6 \times 10^{-5}\,{\it t_6}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_{11}}^2}}
C>    +{{6.4 \times 10^{-8}\,{\it t_5}\,\sigma_{ss}^3}\over{\rho_s^8
C>    \,{\it t_{11}}^3}}+{{2.560000000000001 \times 10^{-10}\,{
C>    \it t_2}\,\sigma_{ss}^4}\over{\rho_s^{{{32}\over{3}}}\,{
C>    \it t_{11}}^4}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] FA Hamprecht, A Cohen, DJ Tozer, NC Handy, J.Chem.Phys. 109
C>    , 6264 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.477267 "> 
C> 10.1063/1.477267 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_hcth_d3(param,tol_rho,ipol,nq,wght, 
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
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**5
            t21 = 1/t4**5
            t22 = 1/rhoa**1.4333333333333334d+1
            t23 = 1/rhoa**1.1666666666666665d+1
            t24 = 1/rhoa**9
            t25 = 1/rhoa**6.333333333333333d+0
            t26 = 1/rhoa**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammaaa*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhoa**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/rhoa**6.666666666666666d-1
            t35 = gammaaa**6
            t36 = 1/t4**6
            t37 = 1/rhoa**18
            t38 = 1/rhoa**1.5333333333333335d+1
            t39 = 1/rhoa**1.2666666666666665d+1
            t40 = 1/rhoa**10
            t41 = 1/rhoa**7.333333333333333d+0
            t42 = 1/rhoa**4.666666666666667d+0
            t43 = 4.6080000000000004d-6*t40*t7*t8*t9+9.102222222222223d-
     1         7*t15*t40*t8*t9-4.8924444444444437d-6*t11*t40*t8*t9-4.027
     2         7333333333337d-8*t2*t39*t5*t7+8.738133333333334d-11*t20*t
     3         21*t38*t7+1.0922666666666667d-8*t11*t2*t39*t5+3.185777777
     4         7777785d-8*t1*t2*t39*t5+3.911111111111111d-2*gammaaa*t15*
     5         t16*t42-3.84d-4*t12*t13*t15*t41+5.404444444444444d-4*t11*
     6         t12*t13*t41-2.730666666666667d-10*t1*t20*t21*t38+5.825422
     7         222222222d-13*t1*t35*t36*t37
            t44 = 1.9200000000000003d-7*t10*t12*t7*t9-1.28d-7*t10*t11*t1
     1         2*t9-7.68d-10*t5*t6*t7*t8+1.0240000000000002d-9*t1*t5*t6*
     2         t8+t33+t29-1.6d-5*gammaaa*t13*t14*t15+3.2d-5*gammaaa*t11*
     3         t13*t14
            t45 = 1/rhoa**17
            t46 = -1.5360000000000001d-6*t12*t24*t7*t9-3.413333333333334
     1         d-7*t12*t15*t24*t9+1.7066666666666666d-6*t11*t12*t24*t9+1
     2         .4336000000000002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t
     3         11*t23*t5*t8-1.0922666666666668d-8*t1*t23*t5*t8-3.2768000
     4         00000001d-11*t2*t21*t22*t7-2.1845333333333336d-13*t1*t20*
     5         t36*t45-1.0666666666666666d-2*t15*t16*t26+1.2800000000000
     6         003d-4*gammaaa*t13*t15*t25-1.7066666666666666d-4*gammaaa*
     7         t11*t13*t25+9.830400000000002d-11*t1*t2*t21*t22
            t47 = 1/rhoa**16
            t48 = 8.192000000000002d-14*t1*t2*t36*t47
            t49 = -1.6384000000000004d-11*t1*t8-1.6d-2*t30
            t50 = -2.3040000000000002d-9*t12*t7
            t51 = 3.0720000000000003d-9*t1*t12
            t52 = t51+t50-1.2d-2*t31
            t53 = -2.56d-7*gammaaa*t11
            t54 = 3.8400000000000006d-7*gammaaa*t7
            t55 = t54+t53-8.0d-3*t32
            t56 = 3.2d-5*t11
            t57 = t56-3.2d-5*t15
            t58 = 1/t4**7
            t59 = 1/rhoa**19
            t60 = 1/rhoa**1.6333333333333333d+1
            t61 = 1/rhoa**1.3666666666666666d+1
            t62 = 1/rhoa**11
            t63 = 1/rhoa**8.333333333333334d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t17*t18*wght
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t18*t27-1.240
     1         700981798799d+0*t17*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t18*t43-4
     1         .13566993932933d-1*t17*t34-2.481401963597598d+0*t19*t27)*
     2         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t18*t46-
     1         1.240700981798799d+0*t19*t44)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t18*(t10*t55*t9+t5*t52*t6+t13*t14*t57+t21*t28*t49+t4
     2         8)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t18*(-
     1         4.6080000000000004d-5*t62*t7*t8*t9-1.7294222222222225d-5*
     2         t15*t62*t8*t9+6.045392592592592d-5*t11*t62*t8*t9+6.576355
     3         555555556d-7*t2*t5*t61*t7-3.058346666666667d-9*t20*t21*t6
     4         0*t7+4.6603377777777777d-12*t35*t36*t59*t7+3.233185185185
     5         185d-3*t12*t13*t15*t63-3.963259259259259d-3*t11*t12*t13*t
     6         63+2.9127111111111115d-8*t15*t2*t5*t61-2.94912d-7*t11*t2*
     7         t5*t61-4.035318518518519d-7*t1*t2*t5*t61+4.66033777777777
     8         8d-10*t11*t20*t21*t60+5.546287407407408d-9*t1*t20*t21*t60
     9         -2.504931555555555d-11*t1*t35*t36*t59+3.7282702222222225d
     :         -14*gammaaa**7*t1*t58/rhoa**2.1666666666666667d+1-1.82518
     ;         51851851852d-1*gammaaa*t15*t16/rhoa**5.666666666666667d+0
     <         )-3.722102945396397d+0*t19*t43-1.240700981798799d+0*t27*t
     =         34+2.7571132928862196d-1*t17/rhoa**1.6666666666666669d+0)
     >         *wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t18*(
     1         1.3824000000000003d-5*t12*t40*t7*t9+5.802666666666667d-6*
     2         t12*t15*t40*t9-1.9000888888888887d-5*t11*t12*t40*t9-2.164
     3         053333333334d-7*t39*t5*t7*t8-1.0922666666666667d-8*t15*t3
     4         9*t5*t8+1.0239999999999999d-7*t11*t39*t5*t8+1.27431111111
     5         11113d-7*t1*t39*t5*t8+1.0813440000000002d-9*t2*t21*t38*t7
     6         -1.7476266666666668d-12*t20*t36*t37*t7-1.3981013333333334
     7         d-14*t1*t35*t58/rhoa**2.0666666666666667d+1+3.91111111111
     8         1111d-2*t15*t16*t42-9.244444444444444d-4*gammaaa*t13*t15*
     9         t41+1.0808888888888889d-3*gammaaa*t11*t13*t41-1.747626666
     :         6666668d-10*t11*t2*t21*t38-1.875057777777778d-9*t1*t2*t21
     ;         *t38+8.956586666666668d-12*t1*t20*t36*t37)-2.481401963597
     <         598d+0*t19*t46-4.13566993932933d-1*t34*t44)*wght+Cmat3(iq
     =         ,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(-9.305257363490993d-1*t18*
     1         (-3.0720000000000003d-6*gammaaa*t24*t7*t9-1.7066666666666
     2         672d-6*gammaaa*t15*t24*t9+4.778666666666666d-6*gammaaa*t1
     3         1*t24*t9-3.6044800000000005d-10*t21*t22*t7*t8+6.5536d-11*
     4         t11*t21*t22*t8+5.679786666666667d-10*t1*t21*t22*t8+6.144d
     5         -8*t12*t23*t5*t7+6.553600000000002d-13*t2*t36*t45*t7+5.24
     6         2880000000001d-15*t1*t20*t58/rhoa**1.9666666666666668d+1+
     7         4.096000000000001d-9*t12*t15*t23*t5-3.2767999999999997d-8
     8         *t11*t12*t23*t5-3.276800000000001d-8*t1*t12*t23*t5-3.0583
     9         466666666675d-12*t1*t2*t36*t45+1.706666666666667d-4*t13*t
     :         15*t25-1.7066666666666666d-4*t11*t13*t25)-1.2407009817987
     ;         99d+0*t19*(3.8400000000000006d-7*gammaaa*t10*t7*t9+1.28d-
     <         7*gammaaa*t10*t15*t9-5.12d-7*gammaaa*t10*t11*t9+1.2288d-1
     =         1*t21*t28*t7*t8-3.276800000000001d-11*t1*t21*t28*t8-4.608
     >         0000000000004d-9*t12*t5*t6*t7+1.5360000000000001d-9*t11*t
     ?         12*t5*t6+3.0720000000000003d-9*t1*t12*t5*t6+t48-3.2d-5*t1
     @         3*t14*t15+3.2d-5*t11*t13*t14))*wght+Cmat3(iq,D3_RA_GAA_GA
     1         A)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t18*(t10*(3.8400000000000006d-7*t7-8.0d-3*t5
     2         7-8.0d-3*(t56-1.6d-5*t15)-2.56d-7*t11)*t9+t36*t47*(3.2768
     3         00000000001d-13*t1*t8-2.0d-2*t49)+t5*t6*(-4.6080000000000
     4         004d-9*gammaaa*t7-1.2d-2*t55-1.2d-2*(t54+t53)+6.144000000
     5         000001d-9*gammaaa*t1)-1.9660800000000006d-15*t1*t2*t58/rh
     6         oa**1.8666666666666668d+1+t21*t28*(-1.6d-2*t52-1.6d-2*(t5
     7         1+t50)-4.915200000000001d-11*t1*t12))*wght
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
            t34 = rhoa**3.333333333333333d-1
            t35 = gammaaa**5
            t36 = 1/t5**5
            t37 = 1/rhoa**1.4333333333333334d+1
            t38 = 1/rhoa**1.1666666666666665d+1
            t39 = 1/rhoa**9
            t40 = 1/rhoa**6.333333333333333d+0
            t41 = 1/rhoa**3.6666666666666664d+0
            t42 = -5.12d-7*t10*t39*t8*t9+3.413333333333333d-7*t10*t12*t3
     1         9*t9+2.0480000000000004d-9*t3*t38*t6*t8-2.730666666666667
     2         4d-9*t2*t3*t38*t6-1.0666666666666666d-2*gammaaa*t16*t17*t
     3         41+4.266666666666667d-5*t13*t14*t16*t40-8.533333333333333
     4         d-5*t12*t13*t14*t40+1.0922666666666668d-11*t2*t35*t36*t37
            t43 = rhob**3.333333333333333d-1
            t44 = gammabb**5
            t45 = 1/t22**5
            t46 = 1/rhob**1.4333333333333334d+1
            t47 = 1/rhob**1.1666666666666665d+1
            t48 = 1/rhob**9
            t49 = 1/rhob**6.333333333333333d+0
            t50 = 1/rhob**3.6666666666666664d+0
            t51 = -5.12d-7*t25*t26*t48*t8+2.0480000000000004d-9*t20*t23*
     1         t47*t8-1.0666666666666666d-2*gammabb*t16*t31*t50+4.266666
     2         666666667d-5*t16*t28*t29*t49-8.533333333333333d-5*t12*t28
     3         *t29*t49+3.413333333333333d-7*t12*t25*t26*t48-2.730666666
     4         6666674d-9*t2*t20*t23*t47+1.0922666666666668d-11*t2*t44*t
     5         45*t46
            t52 = 1/rhoa**1.3333333333333333d+1
            t53 = -4.0960000000000013d-12*t2*t3*t36*t52
            t54 = 1.0240000000000002d-9*t2*t9-7.68d-10*t8*t9
            t55 = 1.9200000000000003d-7*t13*t8-1.28d-7*t12*t13
            t56 = 3.2d-5*gammaaa*t12-1.6d-5*gammaaa*t16
            t57 = 4.0d-3*t16*t17*t4
            t58 = 1/rhob**1.3333333333333333d+1
            t59 = -4.0960000000000013d-12*t2*t20*t45*t58
            t60 = 1.0240000000000002d-9*t2*t25-7.68d-10*t25*t8
            t61 = 1.9200000000000003d-7*t28*t8-1.28d-7*t12*t28
            t62 = 3.2d-5*gammabb*t12-1.6d-5*gammabb*t16
            t63 = 4.0d-3*t16*t21*t31
            t64 = 1/rhoa**6.666666666666666d-1
            t65 = gammaaa**6
            t66 = 1/t5**6
            t67 = 1/rhoa**18
            t68 = 1/rhoa**1.5333333333333335d+1
            t69 = 1/rhoa**1.2666666666666665d+1
            t70 = 1/rhoa**10
            t71 = 1/rhoa**7.333333333333333d+0
            t72 = 1/rhoa**4.666666666666667d+0
            t73 = 4.6080000000000004d-6*t10*t70*t8*t9+9.102222222222223d
     1         -7*t10*t16*t70*t9-4.8924444444444437d-6*t10*t12*t70*t9-4.
     2         0277333333333337d-8*t3*t6*t69*t8+8.738133333333334d-11*t3
     3         5*t36*t68*t8+3.911111111111111d-2*gammaaa*t16*t17*t72-3.8
     4         4d-4*t13*t14*t16*t71+5.404444444444444d-4*t12*t13*t14*t71
     5         +3.1857777777777785d-8*t2*t3*t6*t69+1.0922666666666667d-8
     6         *t12*t3*t6*t69-2.730666666666667d-10*t2*t35*t36*t68+5.825
     7         422222222222d-13*t2*t65*t66*t67
            t74 = 1/rhob**6.666666666666666d-1
            t75 = gammabb**6
            t76 = 1/t22**6
            t77 = 1/rhob**18
            t78 = 1/rhob**1.5333333333333335d+1
            t79 = 1/rhob**1.2666666666666665d+1
            t80 = 1/rhob**10
            t81 = 1/rhob**7.333333333333333d+0
            t82 = 1/rhob**4.666666666666667d+0
            t83 = 3.911111111111111d-2*gammabb*t16*t31*t82-3.84d-4*t16*t
     1         28*t29*t81+5.404444444444444d-4*t12*t28*t29*t81+4.6080000
     2         000000004d-6*t25*t26*t8*t80+9.102222222222223d-7*t16*t25*
     3         t26*t80-4.8924444444444437d-6*t12*t25*t26*t80-4.027733333
     4         3333337d-8*t20*t23*t79*t8+8.738133333333334d-11*t44*t45*t
     5         78*t8+3.1857777777777785d-8*t2*t20*t23*t79+1.092266666666
     6         6667d-8*t12*t20*t23*t79-2.730666666666667d-10*t2*t44*t45*
     7         t78+5.825422222222222d-13*t2*t75*t76*t77
            t84 = -7.68d-10*t6*t7*t8*t9+1.0240000000000002d-9*t2*t6*t7*t
     1         9+1.9200000000000003d-7*t10*t11*t13*t8+t57+t53-1.6d-5*gam
     2         maaa*t14*t15*t16+3.2d-5*gammaaa*t12*t14*t15-1.28d-7*t10*t
     3         11*t12*t13
            t85 = 1/rhoa**17
            t86 = 1.4336000000000002d-8*t38*t6*t8*t9-1.0922666666666668d
     1         -8*t2*t38*t6*t9-4.0959999999999996d-9*t12*t38*t6*t9-2.184
     2         5333333333336d-13*t2*t35*t66*t85-1.5360000000000001d-6*t1
     3         0*t13*t39*t8-3.276800000000001d-11*t3*t36*t37*t8-1.066666
     4         6666666666d-2*t16*t17*t41+1.2800000000000003d-4*gammaaa*t
     5         14*t16*t40-1.7066666666666666d-4*gammaaa*t12*t14*t40-3.41
     6         3333333333334d-7*t10*t13*t16*t39+1.7066666666666666d-6*t1
     7         0*t12*t13*t39+9.830400000000002d-11*t2*t3*t36*t37
            t87 = 1.9200000000000003d-7*t26*t27*t28*t8-7.68d-10*t23*t24*
     1         t25*t8+t63+t59-1.6d-5*gammabb*t16*t29*t30+3.2d-5*gammabb*
     2         t12*t29*t30-1.28d-7*t12*t26*t27*t28+1.0240000000000002d-9
     3         *t2*t23*t24*t25
            t88 = 1/rhob**17
            t89 = -2.1845333333333336d-13*t2*t44*t76*t88-1.5360000000000
     1         001d-6*t26*t28*t48*t8+1.4336000000000002d-8*t23*t25*t47*t
     2         8-3.276800000000001d-11*t20*t45*t46*t8-1.0666666666666666
     3         d-2*t16*t31*t50+1.2800000000000003d-4*gammabb*t16*t29*t49
     4         -1.7066666666666666d-4*gammabb*t12*t29*t49-3.413333333333
     5         334d-7*t16*t26*t28*t48+1.7066666666666666d-6*t12*t26*t28*
     6         t48-1.0922666666666668d-8*t2*t23*t25*t47-4.09599999999999
     7         96d-9*t12*t23*t25*t47+9.830400000000002d-11*t2*t20*t45*t4
     8         6
            t90 = 1/rhoa**16
            t91 = 8.192000000000002d-14*t2*t3*t66*t90
            t92 = -1.6384000000000004d-11*t2*t9-1.6d-2*t54
            t93 = -2.3040000000000002d-9*t13*t8
            t94 = 3.0720000000000003d-9*t13*t2
            t95 = t94+t93-1.2d-2*t55
            t96 = -2.56d-7*gammaaa*t12
            t97 = 3.8400000000000006d-7*gammaaa*t8
            t98 = t97+t96-8.0d-3*t56
            t99 = 3.2d-5*t12
            t100 = t99-3.2d-5*t16
            t101 = 1/rhob**16
            t102 = 8.192000000000002d-14*t101*t2*t20*t76
            t103 = -1.6d-2*t60-1.6384000000000004d-11*t2*t25
            t104 = -2.3040000000000002d-9*t28*t8
            t105 = 3.0720000000000003d-9*t2*t28
            t106 = -1.2d-2*t61+t105+t104
            t107 = -2.56d-7*gammabb*t12
            t108 = 3.8400000000000006d-7*gammabb*t8
            t109 = -8.0d-3*t62+t108+t107
            t110 = 1/t5**7
            t111 = 1/rhoa**19
            t112 = 1/rhoa**1.6333333333333333d+1
            t113 = 1/rhoa**1.3666666666666666d+1
            t114 = 1/rhoa**11
            t115 = 1/rhoa**8.333333333333334d+0
            t116 = 1/t22**7
            t117 = 1/rhob**19
            t118 = 1/rhob**1.6333333333333333d+1
            t119 = 1/rhob**1.3666666666666666d+1
            t120 = 1/rhob**11
            t121 = 1/rhob**8.333333333333334d+0
            t122 = -8.0d-3*(t99-1.6d-5*t16)+3.8400000000000006d-7*t8-2.5
     1         6d-7*t12-8.0d-3*t100
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t32*t33-9.3052573634
     1         90993d-1*t18*t19)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t19*t42-1.240
     1         700981798799d+0*t18*t34)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t33*t51-1.240
     1         700981798799d+0*t32*t43)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t19*(
     1         t54*t6*t7+t57+t14*t15*t56+t10*t11*t55+t53)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t33*(
     1         t63+t29*t30*t62+t26*t27*t61+t23*t24*t60+t59)*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t19*t73-4
     1         .13566993932933d-1*t18*t64-2.481401963597598d+0*t34*t42)*
     2         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363490993d-1*t33*t83-4
     1         .13566993932933d-1*t32*t74-2.481401963597598d+0*t43*t51)*
     2         wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t19*t86-
     1         1.240700981798799d+0*t34*t84)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-9.305257363490993d-1*t33*t89-
     1         1.240700981798799d+0*t43*t87)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t19*(t10*t11*t98+t6*t7*t95+t36*t52*t92+t91+t100*t14*
     2         t15)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t33*(t103*t45*t58+t100*t29*t30+t109*t26*t27+t106*t23
     2         *t24+t102)*wght
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t19*(-
     1         4.6080000000000004d-5*t10*t114*t8*t9-1.7294222222222225d-
     2         5*t10*t114*t16*t9+6.045392592592592d-5*t10*t114*t12*t9+4.
     3         6603377777777777d-12*t111*t65*t66*t8+6.576355555555556d-7
     4         *t113*t3*t6*t8-3.058346666666667d-9*t112*t35*t36*t8-2.504
     5         931555555555d-11*t111*t2*t65*t66-4.035318518518519d-7*t11
     6         3*t2*t3*t6+2.9127111111111115d-8*t113*t16*t3*t6-2.94912d-
     7         7*t113*t12*t3*t6+5.546287407407408d-9*t112*t2*t35*t36+4.6
     8         60337777777778d-10*t112*t12*t35*t36+3.7282702222222225d-1
     9         4*gammaaa**7*t110*t2/rhoa**2.1666666666666667d+1-1.825185
     :         1851851852d-1*gammaaa*t16*t17/rhoa**5.666666666666667d+0+
     ;         3.233185185185185d-3*t115*t13*t14*t16-3.963259259259259d-
     <         3*t115*t12*t13*t14)-3.722102945396397d+0*t34*t73-1.240700
     =         981798799d+0*t42*t64+2.7571132928862196d-1*t18/rhoa**1.66
     >         66666666666669d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(-3.722102945396397d+0*t43*t8
     1         3-9.305257363490993d-1*t33*(4.6603377777777777d-12*t117*t
     2         75*t76*t8-3.058346666666667d-9*t118*t44*t45*t8-4.60800000
     3         00000004d-5*t120*t25*t26*t8+6.576355555555556d-7*t119*t20
     4         *t23*t8-2.504931555555555d-11*t117*t2*t75*t76+5.546287407
     5         407408d-9*t118*t2*t44*t45+4.660337777777778d-10*t118*t12*
     6         t44*t45-1.8251851851851852d-1*gammabb*t16*t31/rhob**5.666
     7         666666666667d+0+3.233185185185185d-3*t121*t16*t28*t29-3.9
     8         63259259259259d-3*t12*t121*t28*t29-1.7294222222222225d-5*
     9         t120*t16*t25*t26+6.045392592592592d-5*t12*t120*t25*t26-4.
     :         035318518518519d-7*t119*t2*t20*t23+2.9127111111111115d-8*
     ;         t119*t16*t20*t23-2.94912d-7*t119*t12*t20*t23+3.7282702222
     <         222225d-14*gammabb**7*t116*t2/rhob**2.1666666666666667d+1
     =         )-1.240700981798799d+0*t51*t74+2.7571132928862196d-1*t32/
     >         rhob**1.6666666666666669d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t19*(
     1         -2.164053333333334d-7*t6*t69*t8*t9+1.2743111111111113d-7*
     2         t2*t6*t69*t9-1.0922666666666667d-8*t16*t6*t69*t9+1.023999
     3         9999999999d-7*t12*t6*t69*t9+1.3824000000000003d-5*t10*t13
     4         *t70*t8+1.0813440000000002d-9*t3*t36*t68*t8-1.74762666666
     5         66668d-12*t35*t66*t67*t8+3.911111111111111d-2*t16*t17*t72
     6         -9.244444444444444d-4*gammaaa*t14*t16*t71+1.0808888888888
     7         889d-3*gammaaa*t12*t14*t71+5.802666666666667d-6*t10*t13*t
     8         16*t70-1.9000888888888887d-5*t10*t12*t13*t70-1.8750577777
     9         77778d-9*t2*t3*t36*t68-1.7476266666666668d-10*t12*t3*t36*
     :         t68+8.956586666666668d-12*t2*t35*t66*t67-1.39810133333333
     ;         34d-14*t110*t2*t65/rhoa**2.0666666666666667d+1)-2.4814019
     <         63597598d+0*t34*t86-4.13566993932933d-1*t64*t84)*wght+Cma
     =         t3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*(-2.481401963597598d+0*t43*t
     1         89-4.13566993932933d-1*t74*t87-9.305257363490993d-1*t33*(
     2         3.911111111111111d-2*t16*t31*t82-9.244444444444444d-4*gam
     3         mabb*t16*t29*t81+1.0808888888888889d-3*gammabb*t12*t29*t8
     4         1+1.3824000000000003d-5*t26*t28*t8*t80+5.802666666666667d
     5         -6*t16*t26*t28*t80-1.9000888888888887d-5*t12*t26*t28*t80-
     6         2.164053333333334d-7*t23*t25*t79*t8+1.0813440000000002d-9
     7         *t20*t45*t78*t8-1.7476266666666668d-12*t44*t76*t77*t8+1.2
     8         743111111111113d-7*t2*t23*t25*t79-1.0922666666666667d-8*t
     9         16*t23*t25*t79+1.0239999999999999d-7*t12*t23*t25*t79-1.87
     :         5057777777778d-9*t2*t20*t45*t78-1.7476266666666668d-10*t1
     ;         2*t20*t45*t78+8.956586666666668d-12*t2*t44*t76*t77-1.3981
     <         013333333334d-14*t116*t2*t75/rhob**2.0666666666666667d+1)
     =         )*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(-1.240700981798799d+0*t34*
     1         (t91+1.2288d-11*t36*t52*t8*t9-3.276800000000001d-11*t2*t3
     2         6*t52*t9-4.6080000000000004d-9*t13*t6*t7*t8+3.84000000000
     3         00006d-7*gammaaa*t10*t11*t8+3.0720000000000003d-9*t13*t2*
     4         t6*t7+1.5360000000000001d-9*t12*t13*t6*t7-3.2d-5*t14*t15*
     5         t16+1.28d-7*gammaaa*t10*t11*t16+3.2d-5*t12*t14*t15-5.12d-
     6         7*gammaaa*t10*t11*t12)-9.305257363490993d-1*t19*(-3.60448
     7         00000000005d-10*t36*t37*t8*t9+5.679786666666667d-10*t2*t3
     8         6*t37*t9+6.5536d-11*t12*t36*t37*t9+6.553600000000002d-13*
     9         t3*t66*t8*t85-3.0583466666666675d-12*t2*t3*t66*t85+6.144d
     :         -8*t13*t38*t6*t8-3.0720000000000003d-6*gammaaa*t10*t39*t8
     ;         -3.276800000000001d-8*t13*t2*t38*t6+4.096000000000001d-9*
     <         t13*t16*t38*t6-3.2767999999999997d-8*t12*t13*t38*t6+1.706
     =         666666666667d-4*t14*t16*t40-1.7066666666666666d-4*t12*t14
     >         *t40-1.7066666666666672d-6*gammaaa*t10*t16*t39+4.77866666
     ?         6666666d-6*gammaaa*t10*t12*t39+5.242880000000001d-15*t110
     @         *t2*t35/rhoa**1.9666666666666668d+1))*wght+Cmat3(iq,D3_RA
     1         _GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*(-9.305257363490993d-1*t33*
     1         (6.553600000000002d-13*t20*t76*t8*t88-3.0583466666666675d
     2         -12*t2*t20*t76*t88-3.0720000000000003d-6*gammabb*t26*t48*
     3         t8+6.144d-8*t23*t28*t47*t8-3.6044800000000005d-10*t25*t45
     4         *t46*t8+1.706666666666667d-4*t16*t29*t49-1.70666666666666
     5         66d-4*t12*t29*t49-1.7066666666666672d-6*gammabb*t16*t26*t
     6         48+4.778666666666666d-6*gammabb*t12*t26*t48-3.27680000000
     7         0001d-8*t2*t23*t28*t47+4.096000000000001d-9*t16*t23*t28*t
     8         47-3.2767999999999997d-8*t12*t23*t28*t47+5.67978666666666
     9         7d-10*t2*t25*t45*t46+6.5536d-11*t12*t25*t45*t46+5.2428800
     :         00000001d-15*t116*t2*t44/rhob**1.9666666666666668d+1)-1.2
     ;         40700981798799d+0*t43*(1.2288d-11*t25*t45*t58*t8-4.608000
     <         0000000004d-9*t23*t24*t28*t8+3.8400000000000006d-7*gammab
     =         b*t26*t27*t8-3.276800000000001d-11*t2*t25*t45*t58-3.2d-5*
     >         t16*t29*t30+3.2d-5*t12*t29*t30+3.0720000000000003d-9*t2*t
     ?         23*t24*t28+1.5360000000000001d-9*t12*t23*t24*t28+1.28d-7*
     @         gammabb*t16*t26*t27-5.12d-7*gammabb*t12*t26*t27+t102))*wg
     1         ht+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t19*(t6*t7*(-1.2d-2*t98-1.2d-2*(t97+t96)-4.6
     2         080000000000004d-9*gammaaa*t8+6.144000000000001d-9*gammaa
     3         a*t2)+t36*t52*(-1.6d-2*t95-1.6d-2*(t94+t93)-4.91520000000
     4         0001d-11*t13*t2)+t66*t90*(3.276800000000001d-13*t2*t9-2.0
     5         d-2*t92)-1.9660800000000006d-15*t110*t2*t3/rhoa**1.866666
     6         6666666668d+1+t10*t11*t122)*wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490993d-1*t33*(t23*t24*(-4.6080000000000004d-9*gammabb
     2         *t8+6.144000000000001d-9*gammabb*t2-1.2d-2*t109-1.2d-2*(t
     3         108+t107))+t101*(3.276800000000001d-13*t2*t25-2.0d-2*t103
     4         )*t76+(-4.915200000000001d-11*t2*t28-1.6d-2*t106-1.6d-2*(
     5         t105+t104))*t45*t58+t122*t26*t27-1.9660800000000006d-15*t
     6         116*t2*t20/rhob**1.8666666666666668d+1)*wght
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
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**5
            t21 = 1/t4**5
            t22 = 1/rhoa**1.4333333333333334d+1
            t23 = 1/rhoa**1.1666666666666665d+1
            t24 = 1/rhoa**9
            t25 = 1/rhoa**6.333333333333333d+0
            t26 = 1/rhoa**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammaaa*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhoa**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammaaa*t11-1.6d-5*gammaaa*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/rhoa**6.666666666666666d-1
            t35 = gammaaa**6
            t36 = 1/t4**6
            t37 = 1/rhoa**18
            t38 = 1/rhoa**1.5333333333333335d+1
            t39 = 1/rhoa**1.2666666666666665d+1
            t40 = 1/rhoa**10
            t41 = 1/rhoa**7.333333333333333d+0
            t42 = 1/rhoa**4.666666666666667d+0
            t43 = 4.6080000000000004d-6*t40*t7*t8*t9+9.102222222222223d-
     1         7*t15*t40*t8*t9-4.8924444444444437d-6*t11*t40*t8*t9-4.027
     2         7333333333337d-8*t2*t39*t5*t7+8.738133333333334d-11*t20*t
     3         21*t38*t7+1.0922666666666667d-8*t11*t2*t39*t5+3.185777777
     4         7777785d-8*t1*t2*t39*t5+3.911111111111111d-2*gammaaa*t15*
     5         t16*t42-3.84d-4*t12*t13*t15*t41+5.404444444444444d-4*t11*
     6         t12*t13*t41-2.730666666666667d-10*t1*t20*t21*t38+5.825422
     7         222222222d-13*t1*t35*t36*t37
            t44 = 1.9200000000000003d-7*t10*t12*t7*t9-1.28d-7*t10*t11*t1
     1         2*t9-7.68d-10*t5*t6*t7*t8+1.0240000000000002d-9*t1*t5*t6*
     2         t8+t33+t29-1.6d-5*gammaaa*t13*t14*t15+3.2d-5*gammaaa*t11*
     3         t13*t14
            t45 = 1/rhoa**17
            t46 = -1.5360000000000001d-6*t12*t24*t7*t9-3.413333333333334
     1         d-7*t12*t15*t24*t9+1.7066666666666666d-6*t11*t12*t24*t9+1
     2         .4336000000000002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t
     3         11*t23*t5*t8-1.0922666666666668d-8*t1*t23*t5*t8-3.2768000
     4         00000001d-11*t2*t21*t22*t7-2.1845333333333336d-13*t1*t20*
     5         t36*t45-1.0666666666666666d-2*t15*t16*t26+1.2800000000000
     6         003d-4*gammaaa*t13*t15*t25-1.7066666666666666d-4*gammaaa*
     7         t11*t13*t25+9.830400000000002d-11*t1*t2*t21*t22
            t47 = 1/rhoa**16
            t48 = 8.192000000000002d-14*t1*t2*t36*t47
            t49 = -1.6384000000000004d-11*t1*t8-1.6d-2*t30
            t50 = -2.3040000000000002d-9*t12*t7
            t51 = 3.0720000000000003d-9*t1*t12
            t52 = t51+t50-1.2d-2*t31
            t53 = -2.56d-7*gammaaa*t11
            t54 = 3.8400000000000006d-7*gammaaa*t7
            t55 = t54+t53-8.0d-3*t32
            t56 = 3.2d-5*t11
            t57 = t56-3.2d-5*t15
            t58 = 1/t4**7
            t59 = 1/rhoa**19
            t60 = 1/rhoa**1.6333333333333333d+1
            t61 = 1/rhoa**1.3666666666666666d+1
            t62 = 1/rhoa**11
            t63 = 1/rhoa**8.333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RA) = -9.305257363490993d-1*t18*t27*wght-1.240700
     1         981798799d+0*t17*t19*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Amat2(iq,D2_RA_RA) = -9.305257363490993d-1*t18*t43*wght-4.13
     1         566993932933d-1*t17*t34*wght-2.481401963597598d+0*t19*t27
     2         *wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490993d-1*t18*t46*wght-1.2
     1         40700981798799d+0*t19*t44*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t18*(t10*t55*t9+t5*t52*t6+t13*t14*t57+t21*t28*t49+t4
     2         8)*wght
            Amat3(iq,D3_RA_RA_RA) = -9.305257363490993d-1*t18*(-4.608000
     1         0000000004d-5*t62*t7*t8*t9-1.7294222222222225d-5*t15*t62*
     2         t8*t9+6.045392592592592d-5*t11*t62*t8*t9+6.57635555555555
     3         6d-7*t2*t5*t61*t7-3.058346666666667d-9*t20*t21*t60*t7+4.6
     4         603377777777777d-12*t35*t36*t59*t7+3.233185185185185d-3*t
     5         12*t13*t15*t63-3.963259259259259d-3*t11*t12*t13*t63+2.912
     6         7111111111115d-8*t15*t2*t5*t61-2.94912d-7*t11*t2*t5*t61-4
     7         .035318518518519d-7*t1*t2*t5*t61+4.660337777777778d-10*t1
     8         1*t20*t21*t60+5.546287407407408d-9*t1*t20*t21*t60-2.50493
     9         1555555555d-11*t1*t35*t36*t59+3.7282702222222225d-14*gamm
     :         aaa**7*t1*t58/rhoa**2.1666666666666667d+1-1.8251851851851
     ;         852d-1*gammaaa*t15*t16/rhoa**5.666666666666667d+0)*wght-3
     <         .722102945396397d+0*t19*t43*wght-1.240700981798799d+0*t27
     =         *t34*wght+2.7571132928862196d-1*t17*wght/rhoa**1.66666666
     >         66666669d+0+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363490993d-1*t18*(1.382400
     1         0000000003d-5*t12*t40*t7*t9+5.802666666666667d-6*t12*t15*
     2         t40*t9-1.9000888888888887d-5*t11*t12*t40*t9-2.16405333333
     3         3334d-7*t39*t5*t7*t8-1.0922666666666667d-8*t15*t39*t5*t8+
     4         1.0239999999999999d-7*t11*t39*t5*t8+1.2743111111111113d-7
     5         *t1*t39*t5*t8+1.0813440000000002d-9*t2*t21*t38*t7-1.74762
     6         66666666668d-12*t20*t36*t37*t7-1.3981013333333334d-14*t1*
     7         t35*t58/rhoa**2.0666666666666667d+1+3.911111111111111d-2*
     8         t15*t16*t42-9.244444444444444d-4*gammaaa*t13*t15*t41+1.08
     9         08888888888889d-3*gammaaa*t11*t13*t41-1.7476266666666668d
     :         -10*t11*t2*t21*t38-1.875057777777778d-9*t1*t2*t21*t38+8.9
     ;         56586666666668d-12*t1*t20*t36*t37)*wght-2.481401963597598
     <         d+0*t19*t46*wght-4.13566993932933d-1*t34*t44*wght+Cmat3(i
     =         q,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363490993d-1*t18*(-3.0720
     1         000000000003d-6*gammaaa*t24*t7*t9-1.7066666666666672d-6*g
     2         ammaaa*t15*t24*t9+4.778666666666666d-6*gammaaa*t11*t24*t9
     3         -3.6044800000000005d-10*t21*t22*t7*t8+6.5536d-11*t11*t21*
     4         t22*t8+5.679786666666667d-10*t1*t21*t22*t8+6.144d-8*t12*t
     5         23*t5*t7+6.553600000000002d-13*t2*t36*t45*t7+5.2428800000
     6         00001d-15*t1*t20*t58/rhoa**1.9666666666666668d+1+4.096000
     7         000000001d-9*t12*t15*t23*t5-3.2767999999999997d-8*t11*t12
     8         *t23*t5-3.276800000000001d-8*t1*t12*t23*t5-3.058346666666
     9         6675d-12*t1*t2*t36*t45+1.706666666666667d-4*t13*t15*t25-1
     :         .7066666666666666d-4*t11*t13*t25)*wght-1.240700981798799d
     ;         +0*t19*(3.8400000000000006d-7*gammaaa*t10*t7*t9+1.28d-7*g
     <         ammaaa*t10*t15*t9-5.12d-7*gammaaa*t10*t11*t9+1.2288d-11*t
     =         21*t28*t7*t8-3.276800000000001d-11*t1*t21*t28*t8-4.608000
     >         0000000004d-9*t12*t5*t6*t7+1.5360000000000001d-9*t11*t12*
     ?         t5*t6+3.0720000000000003d-9*t1*t12*t5*t6+t48-3.2d-5*t13*t
     @         14*t15+3.2d-5*t11*t13*t14)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t18*(t10*(3.8400000000000006d-7*t7-8.0d-3*t5
     2         7-8.0d-3*(t56-1.6d-5*t15)-2.56d-7*t11)*t9+t36*t47*(3.2768
     3         00000000001d-13*t1*t8-2.0d-2*t49)+t5*t6*(-4.6080000000000
     4         004d-9*gammaaa*t7-1.2d-2*t55-1.2d-2*(t54+t53)+6.144000000
     5         000001d-9*gammaaa*t1)-1.9660800000000006d-15*t1*t2*t58/rh
     6         oa**1.8666666666666668d+1+t21*t28*(-1.6d-2*t52-1.6d-2*(t5
     7         1+t50)-4.915200000000001d-11*t1*t12))*wght
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
            t19 = rhob**3.333333333333333d-1
            t20 = gammabb**5
            t21 = 1/t4**5
            t22 = 1/rhob**1.4333333333333334d+1
            t23 = 1/rhob**1.1666666666666665d+1
            t24 = 1/rhob**9
            t25 = 1/rhob**6.333333333333333d+0
            t26 = 1/rhob**3.6666666666666664d+0
            t27 = -5.12d-7*t24*t7*t8*t9+3.413333333333333d-7*t11*t24*t8*
     1         t9+2.0480000000000004d-9*t2*t23*t5*t7-2.7306666666666674d
     2         -9*t1*t2*t23*t5-1.0666666666666666d-2*gammabb*t15*t16*t26
     3         +4.266666666666667d-5*t12*t13*t15*t25-8.533333333333333d-
     4         5*t11*t12*t13*t25+1.0922666666666668d-11*t1*t20*t21*t22
            t28 = 1/rhob**1.3333333333333333d+1
            t29 = -4.0960000000000013d-12*t1*t2*t21*t28
            t30 = 1.0240000000000002d-9*t1*t8-7.68d-10*t7*t8
            t31 = 1.9200000000000003d-7*t12*t7-1.28d-7*t11*t12
            t32 = 3.2d-5*gammabb*t11-1.6d-5*gammabb*t15
            t33 = 4.0d-3*t15*t16*t3
            t34 = 1/rhob**6.666666666666666d-1
            t35 = gammabb**6
            t36 = 1/t4**6
            t37 = 1/rhob**18
            t38 = 1/rhob**1.5333333333333335d+1
            t39 = 1/rhob**1.2666666666666665d+1
            t40 = 1/rhob**10
            t41 = 1/rhob**7.333333333333333d+0
            t42 = 1/rhob**4.666666666666667d+0
            t43 = 4.6080000000000004d-6*t40*t7*t8*t9+9.102222222222223d-
     1         7*t15*t40*t8*t9-4.8924444444444437d-6*t11*t40*t8*t9-4.027
     2         7333333333337d-8*t2*t39*t5*t7+8.738133333333334d-11*t20*t
     3         21*t38*t7+1.0922666666666667d-8*t11*t2*t39*t5+3.185777777
     4         7777785d-8*t1*t2*t39*t5+3.911111111111111d-2*gammabb*t15*
     5         t16*t42-3.84d-4*t12*t13*t15*t41+5.404444444444444d-4*t11*
     6         t12*t13*t41-2.730666666666667d-10*t1*t20*t21*t38+5.825422
     7         222222222d-13*t1*t35*t36*t37
            t44 = 1.9200000000000003d-7*t10*t12*t7*t9-1.28d-7*t10*t11*t1
     1         2*t9-7.68d-10*t5*t6*t7*t8+1.0240000000000002d-9*t1*t5*t6*
     2         t8+t33+t29-1.6d-5*gammabb*t13*t14*t15+3.2d-5*gammabb*t11*
     3         t13*t14
            t45 = 1/rhob**17
            t46 = -1.5360000000000001d-6*t12*t24*t7*t9-3.413333333333334
     1         d-7*t12*t15*t24*t9+1.7066666666666666d-6*t11*t12*t24*t9+1
     2         .4336000000000002d-8*t23*t5*t7*t8-4.0959999999999996d-9*t
     3         11*t23*t5*t8-1.0922666666666668d-8*t1*t23*t5*t8-3.2768000
     4         00000001d-11*t2*t21*t22*t7-2.1845333333333336d-13*t1*t20*
     5         t36*t45-1.0666666666666666d-2*t15*t16*t26+1.2800000000000
     6         003d-4*gammabb*t13*t15*t25-1.7066666666666666d-4*gammabb*
     7         t11*t13*t25+9.830400000000002d-11*t1*t2*t21*t22
            t47 = 1/rhob**16
            t48 = 8.192000000000002d-14*t1*t2*t36*t47
            t49 = -1.6384000000000004d-11*t1*t8-1.6d-2*t30
            t50 = -2.3040000000000002d-9*t12*t7
            t51 = 3.0720000000000003d-9*t1*t12
            t52 = t51+t50-1.2d-2*t31
            t53 = -2.56d-7*gammabb*t11
            t54 = 3.8400000000000006d-7*gammabb*t7
            t55 = t54+t53-8.0d-3*t32
            t56 = 3.2d-5*t11
            t57 = t56-3.2d-5*t15
            t58 = 1/t4**7
            t59 = 1/rhob**19
            t60 = 1/rhob**1.6333333333333333d+1
            t61 = 1/rhob**1.3666666666666666d+1
            t62 = 1/rhob**11
            t63 = 1/rhob**8.333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t17*t18*wght
            Amat(iq,D1_RB) = -9.305257363490993d-1*t18*t27*wght-1.240700
     1         981798799d+0*t17*t19*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t18*(
     1         t10*t31*t9+t30*t5*t6+t33+t13*t14*t32+t29)*wght
            Amat2(iq,D2_RB_RB) = -9.305257363490993d-1*t18*t43*wght-4.13
     1         566993932933d-1*t17*t34*wght-2.481401963597598d+0*t19*t27
     2         *wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490993d-1*t18*t46*wght-1.2
     1         40700981798799d+0*t19*t44*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t18*(t10*t55*t9+t5*t52*t6+t13*t14*t57+t21*t28*t49+t4
     2         8)*wght
            Amat3(iq,D3_RB_RB_RB) = -9.305257363490993d-1*t18*(-4.608000
     1         0000000004d-5*t62*t7*t8*t9-1.7294222222222225d-5*t15*t62*
     2         t8*t9+6.045392592592592d-5*t11*t62*t8*t9+6.57635555555555
     3         6d-7*t2*t5*t61*t7-3.058346666666667d-9*t20*t21*t60*t7+4.6
     4         603377777777777d-12*t35*t36*t59*t7+3.233185185185185d-3*t
     5         12*t13*t15*t63-3.963259259259259d-3*t11*t12*t13*t63+2.912
     6         7111111111115d-8*t15*t2*t5*t61-2.94912d-7*t11*t2*t5*t61-4
     7         .035318518518519d-7*t1*t2*t5*t61+4.660337777777778d-10*t1
     8         1*t20*t21*t60+5.546287407407408d-9*t1*t20*t21*t60-2.50493
     9         1555555555d-11*t1*t35*t36*t59+3.7282702222222225d-14*gamm
     :         abb**7*t1*t58/rhob**2.1666666666666667d+1-1.8251851851851
     ;         852d-1*gammabb*t15*t16/rhob**5.666666666666667d+0)*wght-3
     <         .722102945396397d+0*t19*t43*wght-1.240700981798799d+0*t27
     =         *t34*wght+2.7571132928862196d-1*t17*wght/rhob**1.66666666
     >         66666669d+0+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363490993d-1*t18*(1.382400
     1         0000000003d-5*t12*t40*t7*t9+5.802666666666667d-6*t12*t15*
     2         t40*t9-1.9000888888888887d-5*t11*t12*t40*t9-2.16405333333
     3         3334d-7*t39*t5*t7*t8-1.0922666666666667d-8*t15*t39*t5*t8+
     4         1.0239999999999999d-7*t11*t39*t5*t8+1.2743111111111113d-7
     5         *t1*t39*t5*t8+1.0813440000000002d-9*t2*t21*t38*t7-1.74762
     6         66666666668d-12*t20*t36*t37*t7-1.3981013333333334d-14*t1*
     7         t35*t58/rhob**2.0666666666666667d+1+3.911111111111111d-2*
     8         t15*t16*t42-9.244444444444444d-4*gammabb*t13*t15*t41+1.08
     9         08888888888889d-3*gammabb*t11*t13*t41-1.7476266666666668d
     :         -10*t11*t2*t21*t38-1.875057777777778d-9*t1*t2*t21*t38+8.9
     ;         56586666666668d-12*t1*t20*t36*t37)*wght-2.481401963597598
     <         d+0*t19*t46*wght-4.13566993932933d-1*t34*t44*wght+Cmat3(i
     =         q,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363490993d-1*t18*(-3.0720
     1         000000000003d-6*gammabb*t24*t7*t9-1.7066666666666672d-6*g
     2         ammabb*t15*t24*t9+4.778666666666666d-6*gammabb*t11*t24*t9
     3         -3.6044800000000005d-10*t21*t22*t7*t8+6.5536d-11*t11*t21*
     4         t22*t8+5.679786666666667d-10*t1*t21*t22*t8+6.144d-8*t12*t
     5         23*t5*t7+6.553600000000002d-13*t2*t36*t45*t7+5.2428800000
     6         00001d-15*t1*t20*t58/rhob**1.9666666666666668d+1+4.096000
     7         000000001d-9*t12*t15*t23*t5-3.2767999999999997d-8*t11*t12
     8         *t23*t5-3.276800000000001d-8*t1*t12*t23*t5-3.058346666666
     9         6675d-12*t1*t2*t36*t45+1.706666666666667d-4*t13*t15*t25-1
     :         .7066666666666666d-4*t11*t13*t25)*wght-1.240700981798799d
     ;         +0*t19*(3.8400000000000006d-7*gammabb*t10*t7*t9+1.28d-7*g
     <         ammabb*t10*t15*t9-5.12d-7*gammabb*t10*t11*t9+1.2288d-11*t
     =         21*t28*t7*t8-3.276800000000001d-11*t1*t21*t28*t8-4.608000
     >         0000000004d-9*t12*t5*t6*t7+1.5360000000000001d-9*t11*t12*
     ?         t5*t6+3.0720000000000003d-9*t1*t12*t5*t6+t48-3.2d-5*t13*t
     @         14*t15+3.2d-5*t11*t13*t14)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490993d-1*t18*(t10*(3.8400000000000006d-7*t7-8.0d-3*t5
     2         7-8.0d-3*(t56-1.6d-5*t15)-2.56d-7*t11)*t9+t36*t47*(3.2768
     3         00000000001d-13*t1*t8-2.0d-2*t49)+t5*t6*(-4.6080000000000
     4         004d-9*gammabb*t7-1.2d-2*t55-1.2d-2*(t54+t53)+6.144000000
     5         000001d-9*gammabb*t1)-1.9660800000000006d-15*t1*t2*t58/rh
     6         ob**1.8666666666666668d+1+t21*t28*(-1.6d-2*t52-1.6d-2*(t5
     7         1+t50)-4.915200000000001d-11*t1*t12))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
