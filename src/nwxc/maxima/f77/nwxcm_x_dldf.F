C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_dldf.F 
C> The nwxcm_x_dldf functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_4} &=& 2.0\,\tau_\alpha+{\it t_2}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-2.0\,\tau_\beta\\\\ 
C>   {\it t_8} &=& 2.0\,\tau_\beta+{\it t_6}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-2.0\,\tau_s\\\\ 
C>   {\it t_{12}} &=& 2.0\,\tau_s+{\it t_{10}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_5}\,
C>    \sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(
C>    -{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{
C>    \it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}
C>    \over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{
C>    \it t_8}^4}}\right)+1.0\right)-0.9305257363490997
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678249
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0
C>    \,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028
C>    \,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{
C>    \it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}
C>    \over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_9}\,\sigma_{ss}}
C>    \over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{
C>    \it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{
C>    \it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{
C>    \it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{
C>    \it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz
C>    , Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_dldf(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhoa**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 1/t2**2
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t16*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(5.00501119582
     1         0146d-1*t17*t18*t8-3.2246206206851474d+0*t14*t17*t7-2.487
     2         9069649189497d+0*t15*t17-1.4754769541349888d+1*t12*t13*t1
     3         7+1.9966796007371967d+1*t10*t11*t17)-9.305257363490997d-1
     4         *t16*(1.823102830144513d-5*gammaaa**2*t19/rhoa**6.3333333
     5         33333333d+0-1.540711790047533d-2*gammaaa*t3/rhoa**3.66666
     6         66666666664d+0)*t5-1.2407009817987995d+0*rhoa**3.33333333
     7         3333333d-1*t16*t4)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*(
     1         5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaaa*t
     2         19/rhoa**5.333333333333333d+0)*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*t7+
     2         3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.7602526d+
     3         0*t10*t11)*wght
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = 1/rhob**2.6666666666666666d+0
            t18 = 1.1832860901832053d-3*gammabb*t17+1.0d+0
            t19 = 1/t18
            t20 = 5.777669212678249d-3*gammabb*t17*t19+1.0d+0
            t21 = rhob**1.3333333333333333d+0
            t22 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t23 = t22-2.0d+0*taub
            t24 = t23**4
            t25 = 2.0d+0*taub+t22
            t26 = 1/t25**4
            t27 = t23**3
            t28 = 1/t25**3
            t29 = t23**2
            t30 = 1/t25**2
            t31 = 1/t25
            t32 = 1.0d+0*(-1.637571d-1*t23*t31-1.880028d-1*t29*t30-4.490
     1         609d-1*t27*t28-8.2359d-3*t24*t26)+1.0d+0
            t33 = rhoa**6.666666666666666d-1
            t34 = 1/t9**5
            t35 = 1/t2**2
            t36 = rhob**6.666666666666666d-1
            t37 = 1/t25**5
            t38 = 1/t18**2
            fnc(iq) = (-9.305257363490997d-1*t16*t4*t5-9.305257363490997
     1         d-1*t20*t21*t32)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(5.00501119582
     1         0146d-1*t33*t34*t8-3.2246206206851474d+0*t14*t33*t7-2.487
     2         9069649189497d+0*t15*t33-1.4754769541349888d+1*t12*t13*t3
     3         3+1.9966796007371967d+1*t10*t11*t33)-9.305257363490997d-1
     4         *t16*(1.823102830144513d-5*gammaaa**2*t35/rhoa**6.3333333
     5         33333333d+0-1.540711790047533d-2*gammaaa*t3/rhoa**3.66666
     6         66666666664d+0)*t5-1.2407009817987995d+0*rhoa**3.33333333
     7         3333333d-1*t16*t4)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t21*t32*(1.823102830
     1         144513d-5*gammabb**2*t38/rhob**6.333333333333333d+0-1.540
     2         711790047533d-2*gammabb*t19/rhob**3.6666666666666664d+0)-
     3         9.305257363490997d-1*t20*t21*(5.005011195820146d-1*t24*t3
     4         6*t37-2.4879069649189497d+0*t31*t36-3.2246206206851474d+0
     5         *t23*t30*t36-1.4754769541349888d+1*t28*t29*t36+1.99667960
     6         07371967d+1*t26*t27*t36)-1.2407009817987995d+0*rhob**3.33
     7         3333333333333d-1*t20*t32)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*(
     1         5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaaa*t
     2         35/rhoa**5.333333333333333d+0)*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t21*t
     1         32*(5.777669212678249d-3*t17*t19-6.836635613041924d-6*gam
     2         mabb*t38/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         6.588720000000001d-2*t34*t8+1.0795254000000001d+0*t14*t7+
     2         3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.7602526d+
     3         0*t10*t11)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t20*t21
     1         *(6.588720000000001d-2*t24*t37+3.275142d-1*t31+1.07952540
     2         00000001d+0*t23*t30+3.4463765999999996d+0*t28*t29+2.76025
     3         26d+0*t26*t27)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhoa**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 1/t2**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t16*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t4*t5*(5.005011195820
     1         146d-1*t17*t18*t8-3.2246206206851474d+0*t14*t17*t7-2.4879
     2         069649189497d+0*t15*t17-1.4754769541349888d+1*t12*t13*t17
     3         +1.9966796007371967d+1*t10*t11*t17)*wght-9.30525736349099
     4         7d-1*t16*(1.823102830144513d-5*gammaaa**2*t19/rhoa**6.333
     5         333333333333d+0-1.540711790047533d-2*gammaaa*t3/rhoa**3.6
     6         666666666666664d+0)*t5*wght-1.2407009817987995d+0*rhoa**3
     7         .333333333333333d-1*t16*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*(
     1         5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaaa*t
     2         19/rhoa**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*t7+
     2         3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.7602526d+
     3         0*t10*t11)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t7 = t6-2.0d+0*taub
            t8 = t7**4
            t9 = 2.0d+0*taub+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhob**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 1/t2**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t16*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t4*t5*(5.005011195820
     1         146d-1*t17*t18*t8-3.2246206206851474d+0*t14*t17*t7-2.4879
     2         069649189497d+0*t15*t17-1.4754769541349888d+1*t12*t13*t17
     3         +1.9966796007371967d+1*t10*t11*t17)*wght-9.30525736349099
     4         7d-1*t16*(1.823102830144513d-5*gammabb**2*t19/rhob**6.333
     5         333333333333d+0-1.540711790047533d-2*gammabb*t3/rhob**3.6
     6         666666666666664d+0)*t5*wght-1.2407009817987995d+0*rhob**3
     7         .333333333333333d-1*t16*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t16*(
     1         5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammabb*t
     2         19/rhob**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*(
     1         6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*t7+
     2         3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.7602526d+
     3         0*t10*t11)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_4} &=& 2.0\,\tau_\alpha+{\it t_2}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-2.0\,\tau_\beta\\\\ 
C>   {\it t_8} &=& 2.0\,\tau_\beta+{\it t_6}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-2.0\,\tau_s\\\\ 
C>   {\it t_{12}} &=& 2.0\,\tau_s+{\it t_{10}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_5}\,
C>    \sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(
C>    -{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{
C>    \it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}
C>    \over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{
C>    \it t_8}^4}}\right)+1.0\right)-0.9305257363490997
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678249
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0
C>    \,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028
C>    \,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{
C>    \it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}
C>    \over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_9}\,\sigma_{ss}}
C>    \over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{
C>    \it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{
C>    \it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{
C>    \it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{
C>    \it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz
C>    , Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_dldf_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhoa**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 5.005011195820146d-1*t17*t18*t8-3.2246206206851474d+0*
     1         t14*t17*t7-2.4879069649189497d+0*t15*t17-1.47547695413498
     2         88d+1*t12*t13*t17+1.9966796007371967d+1*t10*t11*t17
            t20 = rhoa**3.333333333333333d-1
            t21 = gammaaa**2
            t22 = 1/t2**2
            t23 = 1/rhoa**6.333333333333333d+0
            t24 = 1/rhoa**3.6666666666666664d+0
            t25 = 1.823102830144513d-5*t21*t22*t23-1.540711790047533d-2*
     1         gammaaa*t24*t3
            t26 = 1/rhoa**5.333333333333333d+0
            t27 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t22*t26
            t28 = 6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*
     1         t7+3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.760252
     2         6d+0*t10*t11
            t29 = 1/t9**6
            t30 = 1/t20
            t31 = 1/t2**3
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t16*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t19*t4*t5-9.30525736
     1         3490997d-1*t16*t25*t5-1.2407009817987995d+0*t16*t20*t4)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*t
     1         27*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t28*t4*
     1         t5*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(-3.801973
     1         231566223d+1*t29*t5*t8+3.336674130546764d-1*t18*t30*t8-3.
     2         503474085248764d+2*t12*t5*t7-2.1497470804567648d+0*t14*t3
     3         0*t7-1.1829796847281436d+3*t11*t18*t5-1.1192644687354772d
     4         +1*t14*t5+1.5825394702560375d+3*t10*t13*t5-1.658604643279
     5         2997d+0*t15*t30-9.836513027566591d+0*t12*t13*t30+1.331119
     6         733824798d+1*t10*t11*t30)-9.305257363490997d-1*t16*(1.150
     7         5345172712729d-7*gammaaa**3*t31/rhoa**10+5.64927656350762
     8         1d-2*gammaaa*t3/rhoa**4.666666666666667d+0-1.640792547130
     9         0616d-4*t21*t22/rhoa**7.333333333333333d+0)*t5-1.86105147
     :         26981994d+0*t19*t25*t5-2.481401963597599d+0*t19*t20*t4-4.
     ;         1356699393293317d-1*t16*t4/t17-2.481401963597599d+0*t16*t
     <         20*t25)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t16*(-4.3145044
     1         39767274d-8*t21*t31/rhoa**9-1.540711790047533d-2*t24*t3+5
     2         .46930849043354d-5*gammaaa*t22*t23)*t5-9.305257363490997d
     3         -1*t19*t27*t5-1.2407009817987995d+0*t16*t20*t27)*wght+Cma
     4         t2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t16*(1.6179391649127278d-8*gammaaa*t31/rhoa**8-1.367
     2         3271226083847d-5*t22*t26)*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(-5.005011
     1         195820146d+0*t17*t29*t8+7.191756064814014d+1*t12*t17*t7-1
     2         .6373837701563188d+2*t11*t17*t18+1.1425055171208194d+1*t1
     3         4*t17-3.127215879613249d+1*t10*t13*t17)-9.305257363490997
     4         d-1*t25*t28*t5-1.2407009817987995d+0*t20*t28*t4)*wght+Mma
     5         t2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t27*t28*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-6.588720000000001d-1*t29*t8-1.8103608000000002d+
     2         1*t12*t7-2.26091184d+1*t11*t18-2.8140792000000003d+0*t14-
     3         3.7239775199999997d+1*t10*t13)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = 1/rhob**2.6666666666666666d+0
            t18 = 1.1832860901832053d-3*gammabb*t17+1.0d+0
            t19 = 1/t18
            t20 = 5.777669212678249d-3*gammabb*t17*t19+1.0d+0
            t21 = rhob**1.3333333333333333d+0
            t22 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t23 = t22-2.0d+0*taub
            t24 = t23**4
            t25 = 2.0d+0*taub+t22
            t26 = 1/t25**4
            t27 = t23**3
            t28 = 1/t25**3
            t29 = t23**2
            t30 = 1/t25**2
            t31 = 1/t25
            t32 = 1.0d+0*(-1.637571d-1*t23*t31-1.880028d-1*t29*t30-4.490
     1         609d-1*t27*t28-8.2359d-3*t24*t26)+1.0d+0
            t33 = rhoa**6.666666666666666d-1
            t34 = 1/t9**5
            t35 = 5.005011195820146d-1*t33*t34*t8-3.2246206206851474d+0*
     1         t14*t33*t7-2.4879069649189497d+0*t15*t33-1.47547695413498
     2         88d+1*t12*t13*t33+1.9966796007371967d+1*t10*t11*t33
            t36 = rhoa**3.333333333333333d-1
            t37 = gammaaa**2
            t38 = 1/t2**2
            t39 = 1/rhoa**6.333333333333333d+0
            t40 = 1/rhoa**3.6666666666666664d+0
            t41 = 1.823102830144513d-5*t37*t38*t39-1.540711790047533d-2*
     1         gammaaa*t3*t40
            t42 = rhob**6.666666666666666d-1
            t43 = 1/t25**5
            t44 = 5.005011195820146d-1*t24*t42*t43-2.4879069649189497d+0
     1         *t31*t42-3.2246206206851474d+0*t23*t30*t42-1.475476954134
     2         9888d+1*t28*t29*t42+1.9966796007371967d+1*t26*t27*t42
            t45 = rhob**3.333333333333333d-1
            t46 = gammabb**2
            t47 = 1/t18**2
            t48 = 1/rhob**6.333333333333333d+0
            t49 = 1/rhob**3.6666666666666664d+0
            t50 = 1.823102830144513d-5*t46*t47*t48-1.540711790047533d-2*
     1         gammabb*t19*t49
            t51 = 1/rhoa**5.333333333333333d+0
            t52 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t38*t51
            t53 = 1/rhob**5.333333333333333d+0
            t54 = 5.777669212678249d-3*t17*t19-6.836635613041924d-6*gamm
     1         abb*t47*t53
            t55 = 6.588720000000001d-2*t34*t8+1.0795254000000001d+0*t14*
     1         t7+3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.760252
     2         6d+0*t10*t11
            t56 = 6.588720000000001d-2*t24*t43+3.275142d-1*t31+1.0795254
     1         000000001d+0*t23*t30+3.4463765999999996d+0*t28*t29+2.7602
     2         526d+0*t26*t27
            t57 = 1/t9**6
            t58 = 1/t36
            t59 = 1/t2**3
            t60 = 1/t25**6
            t61 = 1/t45
            t62 = 1/t18**3
            fnc(iq) = (-9.305257363490997d-1*t16*t4*t5-9.305257363490997
     1         d-1*t20*t21*t32)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t16*t41*t5-9.3052573
     1         63490997d-1*t35*t4*t5-1.2407009817987995d+0*t16*t36*t4)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t21*t32*t50-1.240700
     1         9817987995d+0*t20*t32*t45-9.305257363490997d-1*t20*t21*t4
     2         4)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*t
     1         5*t52*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t21*t
     1         32*t54*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         55*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t20*t21
     1         *t56*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(3.3366741
     1         30546764d-1*t34*t58*t8-3.801973231566223d+1*t5*t57*t8-2.1
     2         497470804567648d+0*t14*t58*t7-3.503474085248764d+2*t12*t5
     3         *t7-1.6586046432792997d+0*t15*t58-9.836513027566591d+0*t1
     4         2*t13*t58+1.331119733824798d+1*t10*t11*t58-1.182979684728
     5         1436d+3*t11*t34*t5-1.1192644687354772d+1*t14*t5+1.5825394
     6         702560375d+3*t10*t13*t5)-9.305257363490997d-1*t16*t5*(1.1
     7         505345172712729d-7*gammaaa**3*t59/rhoa**10-1.640792547130
     8         0616d-4*t37*t38/rhoa**7.333333333333333d+0+5.649276563507
     9         621d-2*gammaaa*t3/rhoa**4.666666666666667d+0)-1.861051472
     :         6981994d+0*t35*t41*t5-2.481401963597599d+0*t16*t36*t41-2.
     ;         481401963597599d+0*t35*t36*t4-4.1356699393293317d-1*t16*t
     <         4/t33)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490997d-1*t21*t32*(1.15053
     1         45172712729d-7*gammabb**3*t62/rhob**10-1.6407925471300616
     2         d-4*t46*t47/rhob**7.333333333333333d+0+5.649276563507621d
     3         -2*gammabb*t19/rhob**4.666666666666667d+0)-9.305257363490
     4         997d-1*t20*t21*(3.336674130546764d-1*t24*t43*t61-1.658604
     5         6432792997d+0*t31*t61-2.1497470804567648d+0*t23*t30*t61-9
     6         .836513027566591d+0*t28*t29*t61+1.331119733824798d+1*t26*
     7         t27*t61-3.801973231566223d+1*t21*t24*t60-1.18297968472814
     8         36d+3*t21*t27*t43-1.1192644687354772d+1*t21*t30+1.5825394
     9         702560375d+3*t21*t26*t29-3.503474085248764d+2*t21*t23*t28
     :         )-2.481401963597599d+0*t32*t45*t50-1.8610514726981994d+0*
     ;         t21*t44*t50-2.481401963597599d+0*t20*t44*t45-4.1356699393
     <         293317d-1*t20*t32/t42)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t16*t5*(-4.3145
     1         04439767274d-8*t37*t59/rhoa**9-1.540711790047533d-2*t3*t4
     2         0+5.46930849043354d-5*gammaaa*t38*t39)-9.305257363490997d
     3         -1*t35*t5*t52-1.2407009817987995d+0*t16*t36*t52)*wght+Cma
     4         t2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363490997d-1*t21*t32*(-4.314
     1         504439767274d-8*t46*t62/rhob**9-1.540711790047533d-2*t19*
     2         t49+5.46930849043354d-5*gammabb*t47*t48)-1.24070098179879
     3         95d+0*t32*t45*t54-9.305257363490997d-1*t21*t44*t54)*wght+
     4         Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t16*t5*(1.6179391649127278d-8*gammaaa*t59/rhoa**8-1.
     2         3673271226083847d-5*t38*t51)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t21*t32*(1.6179391649127278d-8*gammabb*t62/rhob**8-1
     2         .3673271226083847d-5*t47*t53)*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(-5.005011
     1         195820146d+0*t33*t57*t8+7.191756064814014d+1*t12*t33*t7-1
     2         .6373837701563188d+2*t11*t33*t34+1.1425055171208194d+1*t1
     3         4*t33-3.127215879613249d+1*t10*t13*t33)-9.305257363490997
     4         d-1*t41*t5*t55-1.2407009817987995d+0*t36*t4*t55)*wght+Mma
     5         t2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t20*t21*(-5.0050
     1         11195820146d+0*t24*t42*t60-1.6373837701563188d+2*t27*t42*
     2         t43+1.1425055171208194d+1*t30*t42-3.127215879613249d+1*t2
     3         6*t29*t42+7.191756064814014d+1*t23*t28*t42)-9.30525736349
     4         0997d-1*t21*t50*t56-1.2407009817987995d+0*t20*t45*t56)*wg
     5         ht+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t52*t55*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t21*t54*t56*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-6.588720000000001d-1*t57*t8-1.8103608000000002d+
     2         1*t12*t7-2.26091184d+1*t11*t34-2.8140792000000003d+0*t14-
     3         3.7239775199999997d+1*t10*t13)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t20*t21*(-6.588720000000001d-1*t24*t60-2.26091184d+1*t27
     2         *t43-2.8140792000000003d+0*t30-3.7239775199999997d+1*t26*
     3         t29-1.8103608000000002d+1*t23*t28)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhoa**1.6666666666666669d+0
            t7 = t6-2.0d+0*taua
            t8 = t7**4
            t9 = 2.0d+0*taua+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhoa**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 5.005011195820146d-1*t17*t18*t8-3.2246206206851474d+0*
     1         t14*t17*t7-2.4879069649189497d+0*t15*t17-1.47547695413498
     2         88d+1*t12*t13*t17+1.9966796007371967d+1*t10*t11*t17
            t20 = rhoa**3.333333333333333d-1
            t21 = gammaaa**2
            t22 = 1/t2**2
            t23 = 1/rhoa**6.333333333333333d+0
            t24 = 1/rhoa**3.6666666666666664d+0
            t25 = 1.823102830144513d-5*t21*t22*t23-1.540711790047533d-2*
     1         gammaaa*t24*t3
            t26 = 1/rhoa**5.333333333333333d+0
            t27 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t22*t26
            t28 = 6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*
     1         t7+3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.760252
     2         6d+0*t10*t11
            t29 = 1/t9**6
            t30 = 1/t20
            t31 = 1/t2**3
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t16*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t19*t4*t5*wght-9.3052
     1         57363490997d-1*t16*t25*t5*wght-1.2407009817987995d+0*t16*
     2         t20*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t16*t
     1         27*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t28*t4*
     1         t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t4*t5*(-3.8019732
     1         31566223d+1*t29*t5*t8+3.336674130546764d-1*t18*t30*t8-3.5
     2         03474085248764d+2*t12*t5*t7-2.1497470804567648d+0*t14*t30
     3         *t7-1.1829796847281436d+3*t11*t18*t5-1.1192644687354772d+
     4         1*t14*t5+1.5825394702560375d+3*t10*t13*t5-1.6586046432792
     5         997d+0*t15*t30-9.836513027566591d+0*t12*t13*t30+1.3311197
     6         33824798d+1*t10*t11*t30)*wght-9.305257363490997d-1*t16*(1
     7         .1505345172712729d-7*gammaaa**3*t31/rhoa**10+5.6492765635
     8         07621d-2*gammaaa*t3/rhoa**4.666666666666667d+0-1.64079254
     9         71300616d-4*t21*t22/rhoa**7.333333333333333d+0)*t5*wght-1
     :         .8610514726981994d+0*t19*t25*t5*wght-2.481401963597599d+0
     ;         *t19*t20*t4*wght-4.1356699393293317d-1*t16*t4*wght/t17-2.
     <         481401963597599d+0*t16*t20*t25*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t16*(-4.31450443
     1         9767274d-8*t21*t31/rhoa**9-1.540711790047533d-2*t24*t3+5.
     2         46930849043354d-5*gammaaa*t22*t23)*t5*wght-9.305257363490
     3         997d-1*t19*t27*t5*wght-1.2407009817987995d+0*t16*t20*t27*
     4         wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t16*(1.6179391649127278d-8*gammaaa*t31/rhoa**8-1.367
     2         3271226083847d-5*t22*t26)*t5*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t4*t5*(-5.0050111
     1         95820146d+0*t17*t29*t8+7.191756064814014d+1*t12*t17*t7-1.
     2         6373837701563188d+2*t11*t17*t18+1.1425055171208194d+1*t14
     3         *t17-3.127215879613249d+1*t10*t13*t17)*wght-9.30525736349
     4         0997d-1*t25*t28*t5*wght-1.2407009817987995d+0*t20*t28*t4*
     5         wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t27*t28*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-6.588720000000001d-1*t29*t8-1.8103608000000002d+
     2         1*t12*t7-2.26091184d+1*t11*t18-2.8140792000000003d+0*t14-
     3         3.7239775199999997d+1*t10*t13)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = 9.115599744691192d+0*rhob**1.6666666666666669d+0
            t7 = t6-2.0d+0*taub
            t8 = t7**4
            t9 = 2.0d+0*taub+t6
            t10 = 1/t9**4
            t11 = t7**3
            t12 = 1/t9**3
            t13 = t7**2
            t14 = 1/t9**2
            t15 = 1/t9
            t16 = 1.0d+0*(-8.2359d-3*t10*t8-1.637571d-1*t15*t7-1.880028d
     1         -1*t13*t14-4.490609d-1*t11*t12)+1.0d+0
            t17 = rhob**6.666666666666666d-1
            t18 = 1/t9**5
            t19 = 5.005011195820146d-1*t17*t18*t8-3.2246206206851474d+0*
     1         t14*t17*t7-2.4879069649189497d+0*t15*t17-1.47547695413498
     2         88d+1*t12*t13*t17+1.9966796007371967d+1*t10*t11*t17
            t20 = rhob**3.333333333333333d-1
            t21 = gammabb**2
            t22 = 1/t2**2
            t23 = 1/rhob**6.333333333333333d+0
            t24 = 1/rhob**3.6666666666666664d+0
            t25 = 1.823102830144513d-5*t21*t22*t23-1.540711790047533d-2*
     1         gammabb*t24*t3
            t26 = 1/rhob**5.333333333333333d+0
            t27 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammab
     1         b*t22*t26
            t28 = 6.588720000000001d-2*t18*t8+1.0795254000000001d+0*t14*
     1         t7+3.275142d-1*t15+3.4463765999999996d+0*t12*t13+2.760252
     2         6d+0*t10*t11
            t29 = 1/t9**6
            t30 = 1/t20
            t31 = 1/t2**3
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t16*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t19*t4*t5*wght-9.3052
     1         57363490997d-1*t16*t25*t5*wght-1.2407009817987995d+0*t16*
     2         t20*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t16*t
     1         27*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t28*t4*
     1         t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t4*t5*(-3.8019732
     1         31566223d+1*t29*t5*t8+3.336674130546764d-1*t18*t30*t8-3.5
     2         03474085248764d+2*t12*t5*t7-2.1497470804567648d+0*t14*t30
     3         *t7-1.1829796847281436d+3*t11*t18*t5-1.1192644687354772d+
     4         1*t14*t5+1.5825394702560375d+3*t10*t13*t5-1.6586046432792
     5         997d+0*t15*t30-9.836513027566591d+0*t12*t13*t30+1.3311197
     6         33824798d+1*t10*t11*t30)*wght-9.305257363490997d-1*t16*(1
     7         .1505345172712729d-7*gammabb**3*t31/rhob**10+5.6492765635
     8         07621d-2*gammabb*t3/rhob**4.666666666666667d+0-1.64079254
     9         71300616d-4*t21*t22/rhob**7.333333333333333d+0)*t5*wght-1
     :         .8610514726981994d+0*t19*t25*t5*wght-2.481401963597599d+0
     ;         *t19*t20*t4*wght-4.1356699393293317d-1*t16*t4*wght/t17-2.
     <         481401963597599d+0*t16*t20*t25*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t16*(-4.31450443
     1         9767274d-8*t21*t31/rhob**9-1.540711790047533d-2*t24*t3+5.
     2         46930849043354d-5*gammabb*t22*t23)*t5*wght-9.305257363490
     3         997d-1*t19*t27*t5*wght-1.2407009817987995d+0*t16*t20*t27*
     4         wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t16*(1.6179391649127278d-8*gammabb*t31/rhob**8-1.367
     2         3271226083847d-5*t22*t26)*t5*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t4*t5*(-5.0050111
     1         95820146d+0*t17*t29*t8+7.191756064814014d+1*t12*t17*t7-1.
     2         6373837701563188d+2*t11*t17*t18+1.1425055171208194d+1*t14
     3         *t17-3.127215879613249d+1*t10*t13*t17)*wght-9.30525736349
     4         0997d-1*t25*t28*t5*wght-1.2407009817987995d+0*t20*t28*t4*
     5         wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t27*t28*t5*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t4*t5*(-6.588720000000001d-1*t29*t8-1.8103608000000002d+
     2         1*t12*t7-2.26091184d+1*t11*t18-2.8140792000000003d+0*t14-
     3         3.7239775199999997d+1*t10*t13)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& 9.115599744691192\,\rho_\alpha^{{{5}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-2.0\,\tau_\alpha\\\\ 
C>   {\it t_4} &=& 2.0\,\tau_\alpha+{\it t_2}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 9.115599744691192\,\rho_\beta^{{{5}\over{3}}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-2.0\,\tau_\beta\\\\ 
C>   {\it t_8} &=& 2.0\,\tau_\beta+{\it t_6}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& 9.115599744691192\,\rho_s^{{{5}\over{3}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-2.0\,\tau_s\\\\ 
C>   {\it t_{12}} &=& 2.0\,\tau_s+{\it t_{10}}\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_5}\,
C>    \sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}
C>    \,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(
C>    -{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{
C>    \it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}
C>    \over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{
C>    \it t_8}^4}}\right)+1.0\right)-0.9305257363490997
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678249
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205
C>    \,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0
C>    \,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028
C>    \,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{
C>    \it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}
C>    \over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.005777669212678249\,{\it t_9}\,\sigma_{ss}}
C>    \over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}
C>    +1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{
C>    \it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{
C>    \it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{
C>    \it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{
C>    \it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz
C>    , Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_dldf_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 9.115599744691192d+0*t6
            t8 = t7-2.0d+0*taua
            t9 = t8**4
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**4
            t12 = t8**3
            t13 = 1/t10**3
            t14 = t8**2
            t15 = 1/t10**2
            t16 = 1/t10
            t17 = 1.0d+0*(-8.2359d-3*t11*t9-1.637571d-1*t16*t8-1.880028d
     1         -1*t14*t15-4.490609d-1*t12*t13)+1.0d+0
            t18 = rhoa**6.666666666666666d-1
            t19 = 1/t10**5
            t20 = 5.005011195820146d-1*t18*t19*t9-3.2246206206851474d+0*
     1         t15*t18*t8-2.4879069649189497d+0*t16*t18-1.47547695413498
     2         88d+1*t13*t14*t18+1.9966796007371967d+1*t11*t12*t18
            t21 = rhoa**3.333333333333333d-1
            t22 = gammaaa**2
            t23 = 1/t2**2
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = 1.823102830144513d-5*t22*t23*t24-1.540711790047533d-2*
     1         gammaaa*t25*t3
            t27 = 1/rhoa**5.333333333333333d+0
            t28 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t23*t27
            t29 = 6.588720000000001d-2*t19*t9+1.0795254000000001d+0*t15*
     1         t8+3.275142d-1*t16+3.4463765999999996d+0*t13*t14+2.760252
     2         6d+0*t11*t12
            t30 = 1/t10**6
            t31 = 1/t21
            t32 = -3.801973231566223d+1*t30*t5*t9+3.336674130546764d-1*t
     1         19*t31*t9-3.503474085248764d+2*t13*t5*t8-2.14974708045676
     2         48d+0*t15*t31*t8-1.1829796847281436d+3*t12*t19*t5-1.11926
     3         44687354772d+1*t15*t5+1.5825394702560375d+3*t11*t14*t5-1.
     4         6586046432792997d+0*t16*t31-9.836513027566591d+0*t13*t14*
     5         t31+1.331119733824798d+1*t11*t12*t31
            t33 = 1/t18
            t34 = gammaaa**3
            t35 = 1/t2**3
            t36 = 1/rhoa**10
            t37 = 1/rhoa**7.333333333333333d+0
            t38 = 1/rhoa**4.666666666666667d+0
            t39 = 5.649276563507621d-2*gammaaa*t3*t38-1.6407925471300616
     1         d-4*t22*t23*t37+1.1505345172712729d-7*t34*t35*t36
            t40 = 1/rhoa**9
            t41 = -4.314504439767274d-8*t22*t35*t40-1.540711790047533d-2
     1         *t25*t3+5.46930849043354d-5*gammaaa*t23*t24
            t42 = 1/rhoa**8
            t43 = 1.6179391649127278d-8*gammaaa*t35*t42-1.36732712260838
     1         47d-5*t23*t27
            t44 = -5.005011195820146d+0*t18*t30*t9+7.191756064814014d+1*
     1         t13*t18*t8-1.6373837701563188d+2*t12*t18*t19+1.1425055171
     2         208194d+1*t15*t18-3.127215879613249d+1*t11*t14*t18
            t45 = rhoa**2
            t46 = 1/t10**7
            t47 = 1/t5
            t48 = 1/t2**4
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t17*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t20*t4*t5-9.30525736
     1         3490997d-1*t17*t26*t5-1.2407009817987995d+0*t17*t21*t4)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t17*t
     1         28*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t29*t4*
     1         t5*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t32*t4*t5-9.3052
     1         57363490997d-1*t17*t39*t5-1.8610514726981994d+0*t20*t26*t
     2         5-4.1356699393293317d-1*t17*t33*t4-2.481401963597599d+0*t
     3         20*t21*t4-2.481401963597599d+0*t17*t21*t26)*wght+Amat2(iq
     4         ,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t17*t41*t5-9.30
     1         5257363490997d-1*t20*t28*t5-1.2407009817987995d+0*t17*t21
     2         *t28)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t17*t43*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t44*t5-9.3052
     1         57363490997d-1*t26*t29*t5-1.2407009817987995d+0*t21*t29*t
     2         4)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t28*t29*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-6.588720000000001d-1*t30*t9-1.8103608000000002d+
     2         1*t13*t8-2.26091184d+1*t12*t19-2.8140792000000003d+0*t15-
     3         3.7239775199999997d+1*t11*t14)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t4*t5*(-1.112
     1         2247101822547d-1*t19*t47*t9+3.465726621898781d+3*t45*t46*
     2         t9-7.603946463132446d+1*t21*t30*t9+7.16582360152255d-1*t1
     3         5*t47*t8+6.405412170861182d+4*t11*t45*t8-7.00694817049752
     4         8d+2*t13*t21*t8+5.528682144264332d-1*t16*t47+3.2788376758
     5         5553d+0*t13*t14*t47-4.437065779415993d+0*t11*t12*t47+8.75
     6         5259318609024d+4*t12*t30*t45-1.5008982250061145d+5*t14*t1
     7         9*t45-4.9826190159894057d+3*t13*t45-2.3659593694562872d+3
     8         *t12*t19*t21-2.238528937470955d+1*t15*t21+3.1650789405120
     9         75d+3*t11*t14*t21)+2.757113292886221d-1*t17*t4/t6-9.30525
     :         7363490997d-1*t17*(1.089129192450197d-9*gammaaa**4*t48/rh
     ;         oa**1.3666666666666666d+1-2.1860155828154185d-6*t34*t35/r
     <         hoa**11-2.6363290629702235d-1*gammaaa*t3/rhoa**5.66666666
     =         6666667d+0+1.3815068112872864d-3*t22*t23/rhoa**8.33333333
     >         3333334d+0)*t5-2.791577209047299d+0*t20*t39*t5-2.79157720
     ?         9047299d+0*t26*t32*t5+t33*(-1.2407009817987995d+0*t20*t4-
     @         1.2407009817987995d+0*t17*t26)-3.7221029453963983d+0*t21*
     1         t32*t4-3.7221029453963983d+0*t17*t21*t39-7.44420589079279
     2         7d+0*t20*t21*t26)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363490997d-1*t17*(-4.0842
     1         344716882384d-10*t34*t48/rhoa**1.2666666666666665d+1+5.64
     2         9276563507621d-2*t3*t38-3.950056131979778d-4*gammaaa*t23*
     3         t37+7.334657547604365d-7*t22*t35*t36)*t5-1.86105147269819
     4         94d+0*t20*t41*t5-9.305257363490997d-1*t28*t32*t5-2.481401
     5         963597599d+0*t17*t21*t41-4.1356699393293317d-1*t17*t28*t3
     6         3-2.481401963597599d+0*t20*t21*t28)*wght+Cmat3(iq,D3_RA_R
     7         A_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t17*(1.5315
     1         879268830895d-10*t22*t48/rhoa**1.1666666666666665d+1-2.15
     2         7252219883637d-7*gammaaa*t35*t40+7.292411320578052d-5*t23
     3         *t24)*t5-9.305257363490997d-1*t20*t43*t5-1.24070098179879
     4         95d+0*t17*t21*t43)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t17*(4.853817494738183d-8*t35*t42-5.74345472
     2         5811585d-11*gammaaa*t48/rhoa**1.0666666666666666d+1)*t5*w
     3         ght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(4.5623
     1         67877879468d+2*t46*t5*t9-3.336674130546764d+0*t30*t31*t9-
     2         4.228073429874891d+3*t11*t5*t8+4.794504043209342d+1*t13*t
     3         31*t8+1.2133954705806733d+4*t12*t30*t5-5.562437653679438d
     4         +3*t14*t19*t5+7.454653957991719d+2*t13*t5-1.0915891801042
     5         124d+2*t12*t19*t31+7.61670344747213d+0*t15*t31-2.08481058
     6         64088326d+1*t11*t14*t31)-1.8610514726981994d+0*t26*t44*t5
     7         -9.305257363490997d-1*t29*t39*t5-2.481401963597599d+0*t21
     8         *t4*t44-4.1356699393293317d-1*t29*t33*t4-2.48140196359759
     9         9d+0*t21*t26*t29)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_GAA_TA) = (-9.305257363490997d-1*t28*t44*t5-9
     1         .305257363490997d-1*t29*t41*t5-1.2407009817987995d+0*t21*
     2         t28*t29)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t29*t43*t5*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(7.906464000000002d+0*t46*t9+2.575807488d+2*
     2         t11*t8+2.3136215999999998d+2*t12*t30+4.33572912d+2*t14*t1
     3         9+4.74635328d+1*t13)*wght
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 9.115599744691192d+0*t6
            t8 = t7-2.0d+0*taua
            t9 = t8**4
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**4
            t12 = t8**3
            t13 = 1/t10**3
            t14 = t8**2
            t15 = 1/t10**2
            t16 = 1/t10
            t17 = 1.0d+0*(-8.2359d-3*t11*t9-1.637571d-1*t16*t8-1.880028d
     1         -1*t14*t15-4.490609d-1*t12*t13)+1.0d+0
            t18 = 1/rhob**2.6666666666666666d+0
            t19 = 1.1832860901832053d-3*gammabb*t18+1.0d+0
            t20 = 1/t19
            t21 = 5.777669212678249d-3*gammabb*t18*t20+1.0d+0
            t22 = rhob**1.3333333333333333d+0
            t23 = rhob**1.6666666666666669d+0
            t24 = 9.115599744691192d+0*t23
            t25 = t24-2.0d+0*taub
            t26 = t25**4
            t27 = 2.0d+0*taub+t24
            t28 = 1/t27**4
            t29 = t25**3
            t30 = 1/t27**3
            t31 = t25**2
            t32 = 1/t27**2
            t33 = 1/t27
            t34 = 1.0d+0*(-1.637571d-1*t25*t33-1.880028d-1*t31*t32-4.490
     1         609d-1*t29*t30-8.2359d-3*t26*t28)+1.0d+0
            t35 = rhoa**6.666666666666666d-1
            t36 = 1/t10**5
            t37 = 5.005011195820146d-1*t35*t36*t9-3.2246206206851474d+0*
     1         t15*t35*t8-2.4879069649189497d+0*t16*t35-1.47547695413498
     2         88d+1*t13*t14*t35+1.9966796007371967d+1*t11*t12*t35
            t38 = rhoa**3.333333333333333d-1
            t39 = gammaaa**2
            t40 = 1/t2**2
            t41 = 1/rhoa**6.333333333333333d+0
            t42 = 1/rhoa**3.6666666666666664d+0
            t43 = 1.823102830144513d-5*t39*t40*t41-1.540711790047533d-2*
     1         gammaaa*t3*t42
            t44 = rhob**6.666666666666666d-1
            t45 = 1/t27**5
            t46 = 5.005011195820146d-1*t26*t44*t45-2.4879069649189497d+0
     1         *t33*t44-3.2246206206851474d+0*t25*t32*t44-1.475476954134
     2         9888d+1*t30*t31*t44+1.9966796007371967d+1*t28*t29*t44
            t47 = rhob**3.333333333333333d-1
            t48 = gammabb**2
            t49 = 1/t19**2
            t50 = 1/rhob**6.333333333333333d+0
            t51 = 1/rhob**3.6666666666666664d+0
            t52 = 1.823102830144513d-5*t48*t49*t50-1.540711790047533d-2*
     1         gammabb*t20*t51
            t53 = 1/rhoa**5.333333333333333d+0
            t54 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t40*t53
            t55 = 1/rhob**5.333333333333333d+0
            t56 = 5.777669212678249d-3*t18*t20-6.836635613041924d-6*gamm
     1         abb*t49*t55
            t57 = 6.588720000000001d-2*t36*t9+1.0795254000000001d+0*t15*
     1         t8+3.275142d-1*t16+3.4463765999999996d+0*t13*t14+2.760252
     2         6d+0*t11*t12
            t58 = 6.588720000000001d-2*t26*t45+3.275142d-1*t33+1.0795254
     1         000000001d+0*t25*t32+3.4463765999999996d+0*t30*t31+2.7602
     2         526d+0*t28*t29
            t59 = 1/t10**6
            t60 = 1/t38
            t61 = 3.336674130546764d-1*t36*t60*t9-3.801973231566223d+1*t
     1         5*t59*t9-2.1497470804567648d+0*t15*t60*t8-3.5034740852487
     2         64d+2*t13*t5*t8-1.6586046432792997d+0*t16*t60-9.836513027
     3         566591d+0*t13*t14*t60+1.331119733824798d+1*t11*t12*t60-1.
     4         1829796847281436d+3*t12*t36*t5-1.1192644687354772d+1*t15*
     5         t5+1.5825394702560375d+3*t11*t14*t5
            t62 = 1/t35
            t63 = gammaaa**3
            t64 = 1/t2**3
            t65 = 1/rhoa**10
            t66 = 1/rhoa**7.333333333333333d+0
            t67 = 1/rhoa**4.666666666666667d+0
            t68 = 5.649276563507621d-2*gammaaa*t3*t67-1.6407925471300616
     1         d-4*t39*t40*t66+1.1505345172712729d-7*t63*t64*t65
            t69 = 1/t27**6
            t70 = 1/t47
            t71 = 3.336674130546764d-1*t26*t45*t70-1.6586046432792997d+0
     1         *t33*t70-2.1497470804567648d+0*t25*t32*t70-9.836513027566
     2         591d+0*t30*t31*t70+1.331119733824798d+1*t28*t29*t70-3.801
     3         973231566223d+1*t22*t26*t69-1.1829796847281436d+3*t22*t29
     4         *t45-1.1192644687354772d+1*t22*t32+1.5825394702560375d+3*
     5         t22*t28*t31-3.503474085248764d+2*t22*t25*t30
            t72 = 1/t44
            t73 = gammabb**3
            t74 = 1/t19**3
            t75 = 1/rhob**10
            t76 = 1/rhob**7.333333333333333d+0
            t77 = 1/rhob**4.666666666666667d+0
            t78 = 5.649276563507621d-2*gammabb*t20*t77-1.640792547130061
     1         6d-4*t48*t49*t76+1.1505345172712729d-7*t73*t74*t75
            t79 = 1/rhoa**9
            t80 = -4.314504439767274d-8*t39*t64*t79-1.540711790047533d-2
     1         *t3*t42+5.46930849043354d-5*gammaaa*t40*t41
            t81 = 1/rhob**9
            t82 = -4.314504439767274d-8*t48*t74*t81-1.540711790047533d-2
     1         *t20*t51+5.46930849043354d-5*gammabb*t49*t50
            t83 = 1/rhoa**8
            t84 = 1.6179391649127278d-8*gammaaa*t64*t83-1.36732712260838
     1         47d-5*t40*t53
            t85 = 1/rhob**8
            t86 = 1.6179391649127278d-8*gammabb*t74*t85-1.36732712260838
     1         47d-5*t49*t55
            t87 = -5.005011195820146d+0*t35*t59*t9+7.191756064814014d+1*
     1         t13*t35*t8-1.6373837701563188d+2*t12*t35*t36+1.1425055171
     2         208194d+1*t15*t35-3.127215879613249d+1*t11*t14*t35
            t88 = -5.005011195820146d+0*t26*t44*t69-1.6373837701563188d+
     1         2*t29*t44*t45+1.1425055171208194d+1*t32*t44-3.12721587961
     2         3249d+1*t28*t31*t44+7.191756064814014d+1*t25*t30*t44
            t89 = -6.588720000000001d-1*t59*t9-1.8103608000000002d+1*t13
     1         *t8-2.26091184d+1*t12*t36-2.8140792000000003d+0*t15-3.723
     2         9775199999997d+1*t11*t14
            t90 = -6.588720000000001d-1*t26*t69-2.26091184d+1*t29*t45-2.
     1         8140792000000003d+0*t32-3.7239775199999997d+1*t28*t31-1.8
     2         103608000000002d+1*t25*t30
            t91 = rhoa**2
            t92 = 1/t10**7
            t93 = 1/t5
            t94 = 1/t2**4
            t95 = rhob**2
            t96 = 1/t27**7
            t97 = 1/t22
            t98 = 1/t19**4
            fnc(iq) = (-9.305257363490997d-1*t17*t4*t5-9.305257363490997
     1         d-1*t21*t22*t34)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t17*t43*t5-9.3052573
     1         63490997d-1*t37*t4*t5-1.2407009817987995d+0*t17*t38*t4)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t22*t34*t52-1.240700
     1         9817987995d+0*t21*t34*t47-9.305257363490997d-1*t21*t22*t4
     2         6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t17*t
     1         5*t54*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t22*t
     1         34*t56*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         57*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t21*t22
     1         *t58*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t17*t5*t68-4.135
     1         6699393293317d-1*t17*t4*t62-9.305257363490997d-1*t4*t5*t6
     2         1-1.8610514726981994d+0*t37*t43*t5-2.481401963597599d+0*t
     3         17*t38*t43-2.481401963597599d+0*t37*t38*t4)*wght+Amat2(iq
     4         ,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490997d-1*t22*t34*t78-4.13
     1         56699393293317d-1*t21*t34*t72-9.305257363490997d-1*t21*t2
     2         2*t71-2.481401963597599d+0*t34*t47*t52-1.8610514726981994
     3         d+0*t22*t46*t52-2.481401963597599d+0*t21*t46*t47)*wght+Am
     4         at2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t17*t5*t80-9.30
     1         5257363490997d-1*t37*t5*t54-1.2407009817987995d+0*t17*t38
     2         *t54)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363490997d-1*t22*t34*t82-1.2
     1         407009817987995d+0*t34*t47*t56-9.305257363490997d-1*t22*t
     2         46*t56)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t17*t5*t84*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t22*t34*t86*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*t87-9.3052
     1         57363490997d-1*t43*t5*t57-1.2407009817987995d+0*t38*t4*t5
     2         7)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t21*t22*t88-9.30
     1         5257363490997d-1*t22*t52*t58-1.2407009817987995d+0*t21*t4
     2         7*t58)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t54*t57*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t22*t56*t58*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*t89*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t21*t22*t90*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t17*t5*(1.089
     1         129192450197d-9*gammaaa**4*t94/rhoa**1.3666666666666666d+
     2         1-2.1860155828154185d-6*t63*t64/rhoa**11+1.38150681128728
     3         64d-3*t39*t40/rhoa**8.333333333333334d+0-2.63632906297022
     4         35d-1*gammaaa*t3/rhoa**5.666666666666667d+0)-9.3052573634
     5         90997d-1*t4*t5*(-1.1122247101822547d-1*t36*t9*t93+7.16582
     6         360152255d-1*t15*t8*t93+5.528682144264332d-1*t16*t93+3.27
     7         883767585553d+0*t13*t14*t93-4.437065779415993d+0*t11*t12*
     8         t93+3.465726621898781d+3*t9*t91*t92+6.405412170861182d+4*
     9         t11*t8*t91+8.755259318609024d+4*t12*t59*t91-1.50089822500
     :         61145d+5*t14*t36*t91-4.9826190159894057d+3*t13*t91-7.6039
     ;         46463132446d+1*t38*t59*t9-7.006948170497528d+2*t13*t38*t8
     <         -2.3659593694562872d+3*t12*t36*t38-2.238528937470955d+1*t
     =         15*t38+3.165078940512075d+3*t11*t14*t38)-2.79157720904729
     >         9d+0*t37*t5*t68-3.7221029453963983d+0*t17*t38*t68+(-1.240
     ?         7009817987995d+0*t17*t43-1.2407009817987995d+0*t37*t4)*t6
     @         2-2.791577209047299d+0*t43*t5*t61-3.7221029453963983d+0*t
     1         38*t4*t61+2.757113292886221d-1*t17*t4/t6-7.44420589079279
     2         7d+0*t37*t38*t43)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363490997d-1*t22*t34*(1.08
     1         9129192450197d-9*gammabb**4*t98/rhob**1.3666666666666666d
     2         +1-2.1860155828154185d-6*t73*t74/rhob**11+1.3815068112872
     3         864d-3*t48*t49/rhob**8.333333333333334d+0-2.6363290629702
     4         235d-1*gammabb*t20/rhob**5.666666666666667d+0)-9.30525736
     5         3490997d-1*t21*t22*(-1.1122247101822547d-1*t26*t45*t97+5.
     6         528682144264332d-1*t33*t97+7.16582360152255d-1*t25*t32*t9
     7         7+3.27883767585553d+0*t30*t31*t97-4.437065779415993d+0*t2
     8         8*t29*t97+3.465726621898781d+3*t26*t95*t96+8.755259318609
     9         024d+4*t29*t69*t95-1.5008982250061145d+5*t31*t45*t95-4.98
     :         26190159894057d+3*t30*t95+6.405412170861182d+4*t25*t28*t9
     ;         5-7.603946463132446d+1*t26*t47*t69-2.3659593694562872d+3*
     <         t29*t45*t47-2.238528937470955d+1*t32*t47+3.16507894051207
     =         5d+3*t28*t31*t47-7.006948170497528d+2*t25*t30*t47)-3.7221
     >         029453963983d+0*t34*t47*t78-2.791577209047299d+0*t22*t46*
     ?         t78+(-1.2407009817987995d+0*t34*t52-1.2407009817987995d+0
     @         *t21*t46)*t72-2.791577209047299d+0*t22*t52*t71-3.72210294
     1         53963983d+0*t21*t47*t71-7.444205890792797d+0*t46*t47*t52+
     2         2.757113292886221d-1*t21*t34/t23)*wght+Amat3(iq,D3_RB_RB_
     3         RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363490997d-1*t17*t5*(-4.0
     1         842344716882384d-10*t63*t94/rhoa**1.2666666666666665d+1+5
     2         .649276563507621d-2*t3*t67-3.950056131979778d-4*gammaaa*t
     3         40*t66+7.334657547604365d-7*t39*t64*t65)-1.86105147269819
     4         94d+0*t37*t5*t80-2.481401963597599d+0*t17*t38*t80-4.13566
     5         99393293317d-1*t17*t54*t62-9.305257363490997d-1*t5*t54*t6
     6         1-2.481401963597599d+0*t37*t38*t54)*wght+Cmat3(iq,D3_RA_R
     7         A_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-9.305257363490997d-1*t22*t34*(-4.
     1         0842344716882384d-10*t73*t98/rhob**1.2666666666666665d+1+
     2         5.649276563507621d-2*t20*t77-3.950056131979778d-4*gammabb
     3         *t49*t76+7.334657547604365d-7*t48*t74*t75)-2.481401963597
     4         599d+0*t34*t47*t82-1.8610514726981994d+0*t22*t46*t82-4.13
     5         56699393293317d-1*t34*t56*t72-9.305257363490997d-1*t22*t5
     6         6*t71-2.481401963597599d+0*t46*t47*t56)*wght+Cmat3(iq,D3_
     7         RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t17*t5*(1.5
     1         315879268830895d-10*t39*t94/rhoa**1.1666666666666665d+1-2
     2         .157252219883637d-7*gammaaa*t64*t79+7.292411320578052d-5*
     3         t40*t41)-9.305257363490997d-1*t37*t5*t84-1.24070098179879
     4         95d+0*t17*t38*t84)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363490997d-1*t22*t34*(1.
     1         5315879268830895d-10*t48*t98/rhob**1.1666666666666665d+1-
     2         2.157252219883637d-7*gammabb*t74*t81+7.292411320578052d-5
     3         *t49*t50)-1.2407009817987995d+0*t34*t47*t86-9.30525736349
     4         0997d-1*t22*t46*t86)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t17*t5*(4.853817494738183d-8*t64*t83-5.74345
     2         4725811585d-11*gammaaa*t94/rhoa**1.0666666666666666d+1)*w
     3         ght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*t22*t34*(4.853817494738183d-8*t74*t85-5.7434
     2         54725811585d-11*gammabb*t98/rhob**1.0666666666666666d+1)*
     3         wght
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(4.5623
     1         67877879468d+2*t5*t9*t92-3.336674130546764d+0*t59*t60*t9+
     2         4.794504043209342d+1*t13*t60*t8-4.228073429874891d+3*t11*
     3         t5*t8-1.0915891801042124d+2*t12*t36*t60+7.61670344747213d
     4         +0*t15*t60-2.0848105864088326d+1*t11*t14*t60+1.2133954705
     5         806733d+4*t12*t5*t59-5.562437653679438d+3*t14*t36*t5+7.45
     6         4653957991719d+2*t13*t5)-1.8610514726981994d+0*t43*t5*t87
     7         -2.481401963597599d+0*t38*t4*t87-9.305257363490997d-1*t5*
     8         t57*t68-4.1356699393293317d-1*t4*t57*t62-2.48140196359759
     9         9d+0*t38*t43*t57)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (-9.305257363490997d-1*t21*t22*(4.56
     1         2367877879468d+2*t22*t26*t96-3.336674130546764d+0*t26*t69
     2         *t70-1.0915891801042124d+2*t29*t45*t70+7.61670344747213d+
     3         0*t32*t70-2.0848105864088326d+1*t28*t31*t70+4.79450404320
     4         9342d+1*t25*t30*t70+1.2133954705806733d+4*t22*t29*t69-5.5
     5         62437653679438d+3*t22*t31*t45+7.454653957991719d+2*t22*t3
     6         0-4.228073429874891d+3*t22*t25*t28)-1.8610514726981994d+0
     7         *t22*t52*t88-2.481401963597599d+0*t21*t47*t88-9.305257363
     8         490997d-1*t22*t58*t78-4.1356699393293317d-1*t21*t58*t72-2
     9         .481401963597599d+0*t47*t52*t58)*wght+Mmat3(iq,D3_RB_RB_T
     :         B)
            Mmat3(iq,D3_RA_GAA_TA) = (-9.305257363490997d-1*t5*t54*t87-9
     1         .305257363490997d-1*t5*t57*t80-1.2407009817987995d+0*t38*
     2         t54*t57)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_RB_GAA_TA) = Mmat3(iq,D3_RB_GAA_TA)
            Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
            Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
            Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = (-9.305257363490997d-1*t22*t56*t88-
     1         9.305257363490997d-1*t22*t58*t82-1.2407009817987995d+0*t4
     2         7*t56*t58)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)-9.305257363490
     1         997d-1*t5*t57*t84*wght
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAA_GBB_TB) = Mmat3(iq,D3_GAA_GBB_TB)-9.30525736
     1         3490997d-1*t22*t58*t86*wght
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)-9.30525736
     1         3490997d-1*t5*t54*t89*wght
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)-9.3052573634
     1         90997d-1*t22*t56*t90*wght
            Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)-9.3052573634
     1         90997d-1*t4*t5*(7.906464000000002d+0*t9*t92+2.575807488d+
     2         2*t11*t8+2.3136215999999998d+2*t12*t59+4.33572912d+2*t14*
     3         t36+4.74635328d+1*t13)*wght
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)-9.3052573634
     1         90997d-1*t21*t22*(7.906464000000002d+0*t26*t96+2.31362159
     2         99999998d+2*t29*t69+4.33572912d+2*t31*t45+4.74635328d+1*t
     3         30+2.575807488d+2*t25*t28)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 9.115599744691192d+0*t6
            t8 = t7-2.0d+0*taua
            t9 = t8**4
            t10 = 2.0d+0*taua+t7
            t11 = 1/t10**4
            t12 = t8**3
            t13 = 1/t10**3
            t14 = t8**2
            t15 = 1/t10**2
            t16 = 1/t10
            t17 = 1.0d+0*(-8.2359d-3*t11*t9-1.637571d-1*t16*t8-1.880028d
     1         -1*t14*t15-4.490609d-1*t12*t13)+1.0d+0
            t18 = rhoa**6.666666666666666d-1
            t19 = 1/t10**5
            t20 = 5.005011195820146d-1*t18*t19*t9-3.2246206206851474d+0*
     1         t15*t18*t8-2.4879069649189497d+0*t16*t18-1.47547695413498
     2         88d+1*t13*t14*t18+1.9966796007371967d+1*t11*t12*t18
            t21 = rhoa**3.333333333333333d-1
            t22 = gammaaa**2
            t23 = 1/t2**2
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = 1.823102830144513d-5*t22*t23*t24-1.540711790047533d-2*
     1         gammaaa*t25*t3
            t27 = 1/rhoa**5.333333333333333d+0
            t28 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammaa
     1         a*t23*t27
            t29 = 6.588720000000001d-2*t19*t9+1.0795254000000001d+0*t15*
     1         t8+3.275142d-1*t16+3.4463765999999996d+0*t13*t14+2.760252
     2         6d+0*t11*t12
            t30 = 1/t10**6
            t31 = 1/t21
            t32 = -3.801973231566223d+1*t30*t5*t9+3.336674130546764d-1*t
     1         19*t31*t9-3.503474085248764d+2*t13*t5*t8-2.14974708045676
     2         48d+0*t15*t31*t8-1.1829796847281436d+3*t12*t19*t5-1.11926
     3         44687354772d+1*t15*t5+1.5825394702560375d+3*t11*t14*t5-1.
     4         6586046432792997d+0*t16*t31-9.836513027566591d+0*t13*t14*
     5         t31+1.331119733824798d+1*t11*t12*t31
            t33 = 1/t18
            t34 = gammaaa**3
            t35 = 1/t2**3
            t36 = 1/rhoa**10
            t37 = 1/rhoa**7.333333333333333d+0
            t38 = 1/rhoa**4.666666666666667d+0
            t39 = 5.649276563507621d-2*gammaaa*t3*t38-1.6407925471300616
     1         d-4*t22*t23*t37+1.1505345172712729d-7*t34*t35*t36
            t40 = 1/rhoa**9
            t41 = -4.314504439767274d-8*t22*t35*t40-1.540711790047533d-2
     1         *t25*t3+5.46930849043354d-5*gammaaa*t23*t24
            t42 = 1/rhoa**8
            t43 = 1.6179391649127278d-8*gammaaa*t35*t42-1.36732712260838
     1         47d-5*t23*t27
            t44 = -5.005011195820146d+0*t18*t30*t9+7.191756064814014d+1*
     1         t13*t18*t8-1.6373837701563188d+2*t12*t18*t19+1.1425055171
     2         208194d+1*t15*t18-3.127215879613249d+1*t11*t14*t18
            t45 = -6.588720000000001d-1*t30*t9-1.8103608000000002d+1*t13
     1         *t8-2.26091184d+1*t12*t19-2.8140792000000003d+0*t15-3.723
     2         9775199999997d+1*t11*t14
            t46 = rhoa**2
            t47 = 1/t10**7
            t48 = 1/t5
            t49 = 1/t2**4
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t17*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t20*t4*t5*wght-9.3052
     1         57363490997d-1*t17*t26*t5*wght-1.2407009817987995d+0*t17*
     2         t21*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t17*t
     1         28*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t29*t4*
     1         t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t32*t4*t5*wght-9.
     1         305257363490997d-1*t17*t39*t5*wght-1.8610514726981994d+0*
     2         t20*t26*t5*wght-4.1356699393293317d-1*t17*t33*t4*wght-2.4
     3         81401963597599d+0*t20*t21*t4*wght-2.481401963597599d+0*t1
     4         7*t21*t26*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t17*t41*t5*wght-
     1         9.305257363490997d-1*t20*t28*t5*wght-1.2407009817987995d+
     2         0*t17*t21*t28*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t17*t43*t5*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t4*t44*t5*wght-9.
     1         305257363490997d-1*t26*t29*t5*wght-1.2407009817987995d+0*
     2         t21*t29*t4*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t28*t29*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t45*t5*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = t33*(-1.2407009817987995d+0*t20*t4*w
     1         ght-1.2407009817987995d+0*t17*t26*wght)-9.305257363490997
     2         d-1*t4*t5*(-1.1122247101822547d-1*t19*t48*t9+3.4657266218
     3         98781d+3*t46*t47*t9-7.603946463132446d+1*t21*t30*t9+7.165
     4         82360152255d-1*t15*t48*t8+6.405412170861182d+4*t11*t46*t8
     5         -7.006948170497528d+2*t13*t21*t8+5.528682144264332d-1*t16
     6         *t48+3.27883767585553d+0*t13*t14*t48-4.437065779415993d+0
     7         *t11*t12*t48+8.755259318609024d+4*t12*t30*t46-1.500898225
     8         0061145d+5*t14*t19*t46-4.9826190159894057d+3*t13*t46-2.36
     9         59593694562872d+3*t12*t19*t21-2.238528937470955d+1*t15*t2
     :         1+3.165078940512075d+3*t11*t14*t21)*wght+2.75711329288622
     ;         1d-1*t17*t4*wght/t6-9.305257363490997d-1*t17*(1.089129192
     <         450197d-9*gammaaa**4*t49/rhoa**1.3666666666666666d+1-2.18
     =         60155828154185d-6*t34*t35/rhoa**11-2.6363290629702235d-1*
     >         gammaaa*t3/rhoa**5.666666666666667d+0+1.3815068112872864d
     ?         -3*t22*t23/rhoa**8.333333333333334d+0)*t5*wght-2.79157720
     @         9047299d+0*t20*t39*t5*wght-2.791577209047299d+0*t26*t32*t
     1         5*wght-3.7221029453963983d+0*t21*t32*t4*wght-3.7221029453
     2         963983d+0*t17*t21*t39*wght-7.444205890792797d+0*t20*t21*t
     3         26*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363490997d-1*t17*(-4.08423
     1         44716882384d-10*t34*t49/rhoa**1.2666666666666665d+1+5.649
     2         276563507621d-2*t3*t38-3.950056131979778d-4*gammaaa*t23*t
     3         37+7.334657547604365d-7*t22*t35*t36)*t5*wght-1.8610514726
     4         981994d+0*t20*t41*t5*wght-9.305257363490997d-1*t28*t32*t5
     5         *wght-2.481401963597599d+0*t17*t21*t41*wght-4.13566993932
     6         93317d-1*t17*t28*t33*wght-2.481401963597599d+0*t20*t21*t2
     7         8*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363490997d-1*t17*(1.53158
     1         79268830895d-10*t22*t49/rhoa**1.1666666666666665d+1-2.157
     2         252219883637d-7*gammaaa*t35*t40+7.292411320578052d-5*t23*
     3         t24)*t5*wght-9.305257363490997d-1*t20*t43*t5*wght-1.24070
     4         09817987995d+0*t17*t21*t43*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t17*(4.853817494738183d-8*t35*t42-5.74345472
     2         5811585d-11*gammaaa*t49/rhoa**1.0666666666666666d+1)*t5*w
     3         ght
            Mmat3(iq,D3_RA_RA_TA) = -9.305257363490997d-1*t4*t5*(4.56236
     1         7877879468d+2*t47*t5*t9-3.336674130546764d+0*t30*t31*t9-4
     2         .228073429874891d+3*t11*t5*t8+4.794504043209342d+1*t13*t3
     3         1*t8+1.2133954705806733d+4*t12*t30*t5-5.562437653679438d+
     4         3*t14*t19*t5+7.454653957991719d+2*t13*t5-1.09158918010421
     5         24d+2*t12*t19*t31+7.61670344747213d+0*t15*t31-2.084810586
     6         4088326d+1*t11*t14*t31)*wght-1.8610514726981994d+0*t26*t4
     7         4*t5*wght-9.305257363490997d-1*t29*t39*t5*wght-2.48140196
     8         3597599d+0*t21*t4*t44*wght-4.1356699393293317d-1*t29*t33*
     9         t4*wght-2.481401963597599d+0*t21*t26*t29*wght+Mmat3(iq,D3
     :         _RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_GAA_TA) = -9.305257363490997d-1*t28*t44*t5*wg
     1         ht-9.305257363490997d-1*t29*t41*t5*wght-1.240700981798799
     2         5d+0*t21*t28*t29*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t29*t43*t5*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_RA_TA_TA) = -9.305257363490997d-1*t4*t5*(6.00601
     1         3434984174d+1*t18*t47*t9-3.064167287043109d+2*t11*t18*t8+
     2         1.6774238597228797d+3*t12*t18*t30+1.2326075324628512d+3*t
     3         14*t18*t19-1.8953534198111305d+2*t13*t18)*wght-9.30525736
     4         3490997d-1*t26*t45*t5*wght-1.2407009817987995d+0*t21*t4*t
     5         45*wght+Mmat3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)-9.3052573634
     1         90997d-1*t28*t45*t5*wght
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(7.906464000000002d+0*t47*t9+2.575807488d+2*
     2         t11*t8+2.3136215999999998d+2*t12*t30+4.33572912d+2*t14*t1
     3         9+4.74635328d+1*t13)*wght
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832053d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678249d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = rhob**1.6666666666666669d+0
            t7 = 9.115599744691192d+0*t6
            t8 = t7-2.0d+0*taub
            t9 = t8**4
            t10 = 2.0d+0*taub+t7
            t11 = 1/t10**4
            t12 = t8**3
            t13 = 1/t10**3
            t14 = t8**2
            t15 = 1/t10**2
            t16 = 1/t10
            t17 = 1.0d+0*(-8.2359d-3*t11*t9-1.637571d-1*t16*t8-1.880028d
     1         -1*t14*t15-4.490609d-1*t12*t13)+1.0d+0
            t18 = rhob**6.666666666666666d-1
            t19 = 1/t10**5
            t20 = 5.005011195820146d-1*t18*t19*t9-3.2246206206851474d+0*
     1         t15*t18*t8-2.4879069649189497d+0*t16*t18-1.47547695413498
     2         88d+1*t13*t14*t18+1.9966796007371967d+1*t11*t12*t18
            t21 = rhob**3.333333333333333d-1
            t22 = gammabb**2
            t23 = 1/t2**2
            t24 = 1/rhob**6.333333333333333d+0
            t25 = 1/rhob**3.6666666666666664d+0
            t26 = 1.823102830144513d-5*t22*t23*t24-1.540711790047533d-2*
     1         gammabb*t25*t3
            t27 = 1/rhob**5.333333333333333d+0
            t28 = 5.777669212678249d-3*t1*t3-6.836635613041924d-6*gammab
     1         b*t23*t27
            t29 = 6.588720000000001d-2*t19*t9+1.0795254000000001d+0*t15*
     1         t8+3.275142d-1*t16+3.4463765999999996d+0*t13*t14+2.760252
     2         6d+0*t11*t12
            t30 = 1/t10**6
            t31 = 1/t21
            t32 = -3.801973231566223d+1*t30*t5*t9+3.336674130546764d-1*t
     1         19*t31*t9-3.503474085248764d+2*t13*t5*t8-2.14974708045676
     2         48d+0*t15*t31*t8-1.1829796847281436d+3*t12*t19*t5-1.11926
     3         44687354772d+1*t15*t5+1.5825394702560375d+3*t11*t14*t5-1.
     4         6586046432792997d+0*t16*t31-9.836513027566591d+0*t13*t14*
     5         t31+1.331119733824798d+1*t11*t12*t31
            t33 = 1/t18
            t34 = gammabb**3
            t35 = 1/t2**3
            t36 = 1/rhob**10
            t37 = 1/rhob**7.333333333333333d+0
            t38 = 1/rhob**4.666666666666667d+0
            t39 = 5.649276563507621d-2*gammabb*t3*t38-1.6407925471300616
     1         d-4*t22*t23*t37+1.1505345172712729d-7*t34*t35*t36
            t40 = 1/rhob**9
            t41 = -4.314504439767274d-8*t22*t35*t40-1.540711790047533d-2
     1         *t25*t3+5.46930849043354d-5*gammabb*t23*t24
            t42 = 1/rhob**8
            t43 = 1.6179391649127278d-8*gammabb*t35*t42-1.36732712260838
     1         47d-5*t23*t27
            t44 = -5.005011195820146d+0*t18*t30*t9+7.191756064814014d+1*
     1         t13*t18*t8-1.6373837701563188d+2*t12*t18*t19+1.1425055171
     2         208194d+1*t15*t18-3.127215879613249d+1*t11*t14*t18
            t45 = -6.588720000000001d-1*t30*t9-1.8103608000000002d+1*t13
     1         *t8-2.26091184d+1*t12*t19-2.8140792000000003d+0*t15-3.723
     2         9775199999997d+1*t11*t14
            t46 = rhob**2
            t47 = 1/t10**7
            t48 = 1/t5
            t49 = 1/t2**4
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t17*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t20*t4*t5*wght-9.3052
     1         57363490997d-1*t17*t26*t5*wght-1.2407009817987995d+0*t17*
     2         t21*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t17*t
     1         28*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t29*t4*
     1         t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t32*t4*t5*wght-9.
     1         305257363490997d-1*t17*t39*t5*wght-1.8610514726981994d+0*
     2         t20*t26*t5*wght-4.1356699393293317d-1*t17*t33*t4*wght-2.4
     3         81401963597599d+0*t20*t21*t4*wght-2.481401963597599d+0*t1
     4         7*t21*t26*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t17*t41*t5*wght-
     1         9.305257363490997d-1*t20*t28*t5*wght-1.2407009817987995d+
     2         0*t17*t21*t28*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t17*t43*t5*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t4*t44*t5*wght-9.
     1         305257363490997d-1*t26*t29*t5*wght-1.2407009817987995d+0*
     2         t21*t29*t4*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t28*t29*t5*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t4*t45*t5*wght
            Amat3(iq,D3_RB_RB_RB) = t33*(-1.2407009817987995d+0*t20*t4*w
     1         ght-1.2407009817987995d+0*t17*t26*wght)-9.305257363490997
     2         d-1*t4*t5*(-1.1122247101822547d-1*t19*t48*t9+3.4657266218
     3         98781d+3*t46*t47*t9-7.603946463132446d+1*t21*t30*t9+7.165
     4         82360152255d-1*t15*t48*t8+6.405412170861182d+4*t11*t46*t8
     5         -7.006948170497528d+2*t13*t21*t8+5.528682144264332d-1*t16
     6         *t48+3.27883767585553d+0*t13*t14*t48-4.437065779415993d+0
     7         *t11*t12*t48+8.755259318609024d+4*t12*t30*t46-1.500898225
     8         0061145d+5*t14*t19*t46-4.9826190159894057d+3*t13*t46-2.36
     9         59593694562872d+3*t12*t19*t21-2.238528937470955d+1*t15*t2
     :         1+3.165078940512075d+3*t11*t14*t21)*wght+2.75711329288622
     ;         1d-1*t17*t4*wght/t6-9.305257363490997d-1*t17*(1.089129192
     <         450197d-9*gammabb**4*t49/rhob**1.3666666666666666d+1-2.18
     =         60155828154185d-6*t34*t35/rhob**11-2.6363290629702235d-1*
     >         gammabb*t3/rhob**5.666666666666667d+0+1.3815068112872864d
     ?         -3*t22*t23/rhob**8.333333333333334d+0)*t5*wght-2.79157720
     @         9047299d+0*t20*t39*t5*wght-2.791577209047299d+0*t26*t32*t
     1         5*wght-3.7221029453963983d+0*t21*t32*t4*wght-3.7221029453
     2         963983d+0*t17*t21*t39*wght-7.444205890792797d+0*t20*t21*t
     3         26*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363490997d-1*t17*(-4.08423
     1         44716882384d-10*t34*t49/rhob**1.2666666666666665d+1+5.649
     2         276563507621d-2*t3*t38-3.950056131979778d-4*gammabb*t23*t
     3         37+7.334657547604365d-7*t22*t35*t36)*t5*wght-1.8610514726
     4         981994d+0*t20*t41*t5*wght-9.305257363490997d-1*t28*t32*t5
     5         *wght-2.481401963597599d+0*t17*t21*t41*wght-4.13566993932
     6         93317d-1*t17*t28*t33*wght-2.481401963597599d+0*t20*t21*t2
     7         8*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363490997d-1*t17*(1.53158
     1         79268830895d-10*t22*t49/rhob**1.1666666666666665d+1-2.157
     2         252219883637d-7*gammabb*t35*t40+7.292411320578052d-5*t23*
     3         t24)*t5*wght-9.305257363490997d-1*t20*t43*t5*wght-1.24070
     4         09817987995d+0*t17*t21*t43*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*t17*(4.853817494738183d-8*t35*t42-5.74345472
     2         5811585d-11*gammabb*t49/rhob**1.0666666666666666d+1)*t5*w
     3         ght
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = -9.305257363490997d-1*t4*t5*(4.56236
     1         7877879468d+2*t47*t5*t9-3.336674130546764d+0*t30*t31*t9-4
     2         .228073429874891d+3*t11*t5*t8+4.794504043209342d+1*t13*t3
     3         1*t8+1.2133954705806733d+4*t12*t30*t5-5.562437653679438d+
     4         3*t14*t19*t5+7.454653957991719d+2*t13*t5-1.09158918010421
     5         24d+2*t12*t19*t31+7.61670344747213d+0*t15*t31-2.084810586
     6         4088326d+1*t11*t14*t31)*wght-1.8610514726981994d+0*t26*t4
     7         4*t5*wght-9.305257363490997d-1*t29*t39*t5*wght-2.48140196
     8         3597599d+0*t21*t4*t44*wght-4.1356699393293317d-1*t29*t33*
     9         t4*wght-2.481401963597599d+0*t21*t26*t29*wght+Mmat3(iq,D3
     :         _RB_RB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = -9.305257363490997d-1*t28*t44*t5*wg
     1         ht-9.305257363490997d-1*t29*t41*t5*wght-1.240700981798799
     2         5d+0*t21*t28*t29*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)-9.30525736
     1         3490997d-1*t29*t43*t5*wght
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = -9.305257363490997d-1*t4*t5*(6.00601
     1         3434984174d+1*t18*t47*t9-3.064167287043109d+2*t11*t18*t8+
     2         1.6774238597228797d+3*t12*t18*t30+1.2326075324628512d+3*t
     3         14*t18*t19-1.8953534198111305d+2*t13*t18)*wght-9.30525736
     4         3490997d-1*t26*t45*t5*wght-1.2407009817987995d+0*t21*t4*t
     5         45*wght+Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
            Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)-9.3052573634
     1         90997d-1*t28*t45*t5*wght
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)-9.305257363490
     1         997d-1*t4*t5*(7.906464000000002d+0*t47*t9+2.575807488d+2*
     2         t11*t8+2.3136215999999998d+2*t12*t30+4.33572912d+2*t14*t1
     3         9+4.74635328d+1*t13)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
