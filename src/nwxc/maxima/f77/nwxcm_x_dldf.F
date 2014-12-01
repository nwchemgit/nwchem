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
C>   {\it t_2} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-1.0\\\\ 
C>   {\it t_4} &=& {\it t_2}+1.0\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-1.0\\\\ 
C>   {\it t_8} &=& {\it t_6}+1.0\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-1.0\\\\ 
C>   {\it t_{12}} &=& {\it t_{10}}+1.0\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_x_dldf(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/t2**2
            t20 = rhoa**6.666666666666666d-1
            t21 = 1/t11**5
            t22 = 1/taua**2
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t18*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-1.6123103103
     1         425737d+0*t16*t20*t7*t9+2.502505597910073d-1*t10*t20*t21*
     2         t7-1.2439534824594749d+0*t17*t20*t7-7.377384770674944d+0*
     3         t14*t15*t20*t7+9.983398003685984d+0*t12*t13*t20*t7)-9.305
     4         257363490997d-1*t18*(1.8231028301445137d-5*gammaaa**2*t19
     5         /rhoa**6.333333333333333d+0-1.5407117900475334d-2*gammaaa
     6         *t3/rhoa**3.6666666666666664d+0)*t5-1.2407009817987995d+0
     7         *rhoa**3.333333333333333d-1*t18*t4)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*(
     1         5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaaa*t
     2         19/rhoa**5.333333333333333d+0)*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         9.673861862055442d-1*t16*t22*t6*t9-1.5015033587460438d-1*
     2         t10*t21*t22*t6+7.463720894756849d-1*t17*t22*t6+4.42643086
     3         2404966d+0*t14*t15*t22*t6-5.99003880221159d+0*t12*t13*t22
     4         *t6)*wght
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
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = 1.1832860901832055d-3*gammabb*t19+1.0d+0
            t21 = 1/t20
            t22 = 5.777669212678251d-3*gammabb*t19*t21+1.0d+0
            t23 = rhob**1.3333333333333333d+0
            t24 = rhob**1.6666666666666669d+0
            t25 = 1/taub
            t26 = 4.557799872345596d+0*t24*t25
            t27 = t26-1.0d+0
            t28 = t27**4
            t29 = t26+1.0d+0
            t30 = 1/t29**4
            t31 = t27**3
            t32 = 1/t29**3
            t33 = t27**2
            t34 = 1/t29**2
            t35 = 1/t29
            t36 = 1.0d+0*(-1.637571d-1*t27*t35-1.880028d-1*t33*t34-4.490
     1         609d-1*t31*t32-8.2359d-3*t28*t30)+1.0d+0
            t37 = 1/t2**2
            t38 = rhoa**6.666666666666666d-1
            t39 = 1/t11**5
            t40 = 1/t20**2
            t41 = rhob**6.666666666666666d-1
            t42 = 1/t29**5
            t43 = 1/taua**2
            t44 = 1/taub**2
            fnc(iq) = (-9.305257363490997d-1*t18*t4*t5-9.305257363490997
     1         d-1*t22*t23*t36)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-1.6123103103
     1         425737d+0*t16*t38*t7*t9+2.502505597910073d-1*t10*t38*t39*
     2         t7-1.2439534824594749d+0*t17*t38*t7-7.377384770674944d+0*
     3         t14*t15*t38*t7+9.983398003685984d+0*t12*t13*t38*t7)-9.305
     4         257363490997d-1*t18*(1.8231028301445137d-5*gammaaa**2*t37
     5         /rhoa**6.333333333333333d+0-1.5407117900475334d-2*gammaaa
     6         *t3/rhoa**3.6666666666666664d+0)*t5-1.2407009817987995d+0
     7         *rhoa**3.333333333333333d-1*t18*t4)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t22*t23*(2.502505597
     1         910073d-1*t25*t28*t41*t42-1.2439534824594749d+0*t25*t35*t
     2         41-1.6123103103425737d+0*t25*t27*t34*t41-7.37738477067494
     3         4d+0*t25*t32*t33*t41+9.983398003685984d+0*t25*t30*t31*t41
     4         )-9.305257363490997d-1*t23*t36*(1.8231028301445137d-5*gam
     5         mabb**2*t40/rhob**6.333333333333333d+0-1.5407117900475334
     6         d-2*gammabb*t21/rhob**3.6666666666666664d+0)-1.2407009817
     7         987995d+0*rhob**3.333333333333333d-1*t22*t36)*wght+Amat(i
     8         q,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*(
     1         5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaaa*t
     2         37/rhoa**5.333333333333333d+0)*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t23*t
     1         36*(5.777669212678251d-3*t19*t21-6.836635613041927d-6*gam
     2         mabb*t40/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         9.673861862055442d-1*t16*t43*t6*t9-1.5015033587460438d-1*
     2         t10*t39*t43*t6+7.463720894756849d-1*t17*t43*t6+4.42643086
     3         2404966d+0*t14*t15*t43*t6-5.99003880221159d+0*t12*t13*t43
     4         *t6)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t22*t23
     1         *(-1.5015033587460438d-1*t24*t28*t42*t44+7.46372089475684
     2         9d-1*t24*t35*t44+9.673861862055442d-1*t24*t27*t34*t44+4.4
     3         26430862404966d+0*t24*t32*t33*t44-5.99003880221159d+0*t24
     4         *t30*t31*t44)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/t2**2
            t20 = rhoa**6.666666666666666d-1
            t21 = 1/t11**5
            t22 = 1/taua**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t4*t5*(-1.61231031034
     1         25737d+0*t16*t20*t7*t9+2.502505597910073d-1*t10*t20*t21*t
     2         7-1.2439534824594749d+0*t17*t20*t7-7.377384770674944d+0*t
     3         14*t15*t20*t7+9.983398003685984d+0*t12*t13*t20*t7)*wght-9
     4         .305257363490997d-1*t18*(1.8231028301445137d-5*gammaaa**2
     5         *t19/rhoa**6.333333333333333d+0-1.5407117900475334d-2*gam
     6         maaa*t3/rhoa**3.6666666666666664d+0)*t5*wght-1.2407009817
     7         987995d+0*rhoa**3.333333333333333d-1*t18*t4*wght+Amat(iq,
     8         D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*(
     1         5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaaa*t
     2         19/rhoa**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         9.673861862055442d-1*t16*t22*t6*t9-1.5015033587460438d-1*
     2         t10*t21*t22*t6+7.463720894756849d-1*t17*t22*t6+4.42643086
     3         2404966d+0*t14*t15*t22*t6-5.99003880221159d+0*t12*t13*t22
     4         *t6)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = rhob**1.6666666666666669d+0
            t7 = 1/taub
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/t2**2
            t20 = rhob**6.666666666666666d-1
            t21 = 1/t11**5
            t22 = 1/taub**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t4*t5*(-1.61231031034
     1         25737d+0*t16*t20*t7*t9+2.502505597910073d-1*t10*t20*t21*t
     2         7-1.2439534824594749d+0*t17*t20*t7-7.377384770674944d+0*t
     3         14*t15*t20*t7+9.983398003685984d+0*t12*t13*t20*t7)*wght-9
     4         .305257363490997d-1*t18*(1.8231028301445137d-5*gammabb**2
     5         *t19/rhob**6.333333333333333d+0-1.5407117900475334d-2*gam
     6         mabb*t3/rhob**3.6666666666666664d+0)*t5*wght-1.2407009817
     7         987995d+0*rhob**3.333333333333333d-1*t18*t4*wght+Amat(iq,
     8         D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t18*(
     1         5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammabb*t
     2         19/rhob**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*(
     1         9.673861862055442d-1*t16*t22*t6*t9-1.5015033587460438d-1*
     2         t10*t21*t22*t6+7.463720894756849d-1*t17*t22*t6+4.42643086
     3         2404966d+0*t14*t15*t22*t6-5.99003880221159d+0*t12*t13*t22
     4         *t6)*wght
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
C>   {\it t_2} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-1.0\\\\ 
C>   {\it t_4} &=& {\it t_2}+1.0\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-1.0\\\\ 
C>   {\it t_8} &=& {\it t_6}+1.0\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-1.0\\\\ 
C>   {\it t_{12}} &=& {\it t_{10}}+1.0\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_x_dldf_d2(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**2
            t21 = 1/t2**2
            t22 = 1/rhoa**6.333333333333333d+0
            t23 = 1/rhoa**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammaaa*t23*t3
            t25 = rhoa**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhoa**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t21*t28
            t30 = 1/taua**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t2**3
            t33 = 1/t11**6
            t34 = 1/t19
            t35 = rhoa**2.3333333333333334d+0
            t36 = 1/taua**3
            t37 = rhoa**3.3333333333333337d+0
            t38 = 1/taua**4
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t18*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t27*t4*t5-9.30525736
     1         3490997d-1*t18*t24*t5-1.2407009817987995d+0*t18*t19*t4)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(-1.074873
     1         5402283824d+0*t16*t34*t7*t9-8.758685213121911d+1*t14*t30*
     2         t5*t9+1.668337065273382d-1*t10*t26*t34*t7-8.2930232163964
     3         99d-1*t17*t34*t7-4.9182565137832956d+0*t14*t15*t34*t7+6.6
     4         5559866912399d+0*t12*t13*t34*t7-9.504933078915558d+0*t10*
     5         t30*t33*t5-2.9574492118203594d+2*t13*t26*t30*t5-2.7981611
     6         71838693d+0*t16*t30*t5+3.956348675640094d+2*t12*t15*t30*t
     7         5)-9.305257363490997d-1*t18*(1.1505345172712736d-7*gammaa
     8         a**3*t32/rhoa**10+5.649276563507623d-2*gammaaa*t3/rhoa**4
     9         .666666666666667d+0-1.6407925471300622d-4*t20*t21/rhoa**7
     :         .333333333333333d+0)*t5-1.8610514726981994d+0*t24*t27*t5-
     ;         2.481401963597599d+0*t19*t27*t4-4.1356699393293317d-1*t18
     <         *t4/t25-2.481401963597599d+0*t18*t19*t24)*wght+Amat2(iq,D
     =         2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t18*(-4.3145044
     1         397672766d-8*t20*t32/rhoa**9-1.5407117900475334d-2*t23*t3
     2         +5.469308490433541d-5*gammaaa*t21*t22)*t5-9.3052573634909
     3         97d-1*t27*t29*t5-1.2407009817987995d+0*t18*t19*t29)*wght+
     4         Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*(1.6179391649127287d-8*gammaaa*t32/rhoa**8-1.367
     2         3271226083852d-5*t21*t28)*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(5.2552111
     1         27873146d+1*t14*t35*t36*t9+1.6123103103425737d+0*t16*t25*
     2         t30*t9+5.702959847349336d+0*t10*t33*t35*t36+1.77446952709
     3         22156d+2*t13*t26*t35*t36+1.6788967031032173d+0*t16*t35*t3
     4         6-2.3738092053840557d+2*t12*t15*t35*t36-2.502505597910073
     5         d-1*t10*t25*t26*t30+1.2439534824594749d+0*t17*t25*t30+7.3
     6         77384770674944d+0*t14*t15*t25*t30-9.983398003685984d+0*t1
     7         2*t13*t25*t30)-9.305257363490997d-1*t24*t31*t5-1.24070098
     8         17987995d+0*t19*t31*t4)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-1.9347723724110882d+0*t16*t36*t6*t9-3.1531266767
     2         23888d+1*t14*t37*t38*t9+3.0030067174920877d-1*t10*t26*t36
     3         *t6-1.4927441789513698d+0*t17*t36*t6-8.852861724809932d+0
     4         *t14*t15*t36*t6+1.1980077604423181d+1*t12*t13*t36*t6-3.42
     5         17759084096017d+0*t10*t33*t37*t38-1.0646817162553295d+2*t
     6         13*t26*t37*t38-1.0073380218619295d+0*t16*t37*t38+1.424285
     7         5232304333d+2*t12*t15*t37*t38)*wght
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
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = 1.1832860901832055d-3*gammabb*t19+1.0d+0
            t21 = 1/t20
            t22 = 5.777669212678251d-3*gammabb*t19*t21+1.0d+0
            t23 = rhob**1.3333333333333333d+0
            t24 = rhob**1.6666666666666669d+0
            t25 = 1/taub
            t26 = 4.557799872345596d+0*t24*t25
            t27 = t26-1.0d+0
            t28 = t27**4
            t29 = t26+1.0d+0
            t30 = 1/t29**4
            t31 = t27**3
            t32 = 1/t29**3
            t33 = t27**2
            t34 = 1/t29**2
            t35 = 1/t29
            t36 = 1.0d+0*(-1.637571d-1*t27*t35-1.880028d-1*t33*t34-4.490
     1         609d-1*t31*t32-8.2359d-3*t28*t30)+1.0d+0
            t37 = rhoa**3.333333333333333d-1
            t38 = gammaaa**2
            t39 = 1/t2**2
            t40 = 1/rhoa**6.333333333333333d+0
            t41 = 1/rhoa**3.6666666666666664d+0
            t42 = 1.8231028301445137d-5*t38*t39*t40-1.5407117900475334d-
     1         2*gammaaa*t3*t41
            t43 = rhoa**6.666666666666666d-1
            t44 = 1/t11**5
            t45 = -1.6123103103425737d+0*t16*t43*t7*t9+2.502505597910073
     1         d-1*t10*t43*t44*t7-1.2439534824594749d+0*t17*t43*t7-7.377
     2         384770674944d+0*t14*t15*t43*t7+9.983398003685984d+0*t12*t
     3         13*t43*t7
            t46 = rhob**3.333333333333333d-1
            t47 = gammabb**2
            t48 = 1/t20**2
            t49 = 1/rhob**6.333333333333333d+0
            t50 = 1/rhob**3.6666666666666664d+0
            t51 = 1.8231028301445137d-5*t47*t48*t49-1.5407117900475334d-
     1         2*gammabb*t21*t50
            t52 = rhob**6.666666666666666d-1
            t53 = 1/t29**5
            t54 = 2.502505597910073d-1*t25*t28*t52*t53-1.243953482459474
     1         9d+0*t25*t35*t52-1.6123103103425737d+0*t25*t27*t34*t52-7.
     2         377384770674944d+0*t25*t32*t33*t52+9.983398003685984d+0*t
     3         25*t30*t31*t52
            t55 = 1/rhoa**5.333333333333333d+0
            t56 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t39*t55
            t57 = 1/rhob**5.333333333333333d+0
            t58 = 5.777669212678251d-3*t19*t21-6.836635613041927d-6*gamm
     1         abb*t48*t57
            t59 = 1/taua**2
            t60 = 9.673861862055442d-1*t16*t59*t6*t9-1.5015033587460438d
     1         -1*t10*t44*t59*t6+7.463720894756849d-1*t17*t59*t6+4.42643
     2         0862404966d+0*t14*t15*t59*t6-5.99003880221159d+0*t12*t13*
     3         t59*t6
            t61 = 1/taub**2
            t62 = -1.5015033587460438d-1*t24*t28*t53*t61+7.4637208947568
     1         49d-1*t24*t35*t61+9.673861862055442d-1*t24*t27*t34*t61+4.
     2         426430862404966d+0*t24*t32*t33*t61-5.99003880221159d+0*t2
     3         4*t30*t31*t61
            t63 = 1/t2**3
            t64 = 1/t11**6
            t65 = 1/t37
            t66 = 1/t20**3
            t67 = 1/t29**6
            t68 = 1/t46
            t69 = rhoa**2.3333333333333334d+0
            t70 = 1/taua**3
            t71 = rhob**2.3333333333333334d+0
            t72 = 1/taub**3
            t73 = rhoa**3.3333333333333337d+0
            t74 = 1/taua**4
            t75 = rhob**3.3333333333333337d+0
            t76 = 1/taub**4
            fnc(iq) = (-9.305257363490997d-1*t18*t4*t5-9.305257363490997
     1         d-1*t22*t23*t36)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t45*t5-9.30525736
     1         3490997d-1*t18*t42*t5-1.2407009817987995d+0*t18*t37*t4)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t22*t23*t54-9.305257
     1         363490997d-1*t23*t36*t51-1.2407009817987995d+0*t22*t36*t4
     2         6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         5*t56*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t23*t
     1         36*t58*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         60*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t22*t23
     1         *t62*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*(-1.074873
     1         5402283824d+0*t16*t65*t7*t9-8.758685213121911d+1*t14*t5*t
     2         59*t9+1.668337065273382d-1*t10*t44*t65*t7-8.2930232163964
     3         99d-1*t17*t65*t7-4.9182565137832956d+0*t14*t15*t65*t7+6.6
     4         5559866912399d+0*t12*t13*t65*t7-9.504933078915558d+0*t10*
     5         t5*t59*t64-2.9574492118203594d+2*t13*t44*t5*t59-2.7981611
     6         71838693d+0*t16*t5*t59+3.956348675640094d+2*t12*t15*t5*t5
     7         9)-9.305257363490997d-1*t18*t5*(1.1505345172712736d-7*gam
     8         maaa**3*t63/rhoa**10-1.6407925471300622d-4*t38*t39/rhoa**
     9         7.333333333333333d+0+5.649276563507623d-2*gammaaa*t3/rhoa
     :         **4.666666666666667d+0)-1.8610514726981994d+0*t42*t45*t5-
     ;         2.481401963597599d+0*t37*t4*t45-4.1356699393293317d-1*t18
     <         *t4/t43-2.481401963597599d+0*t18*t37*t42)*wght+Amat2(iq,D
     =         2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490997d-1*t22*t23*(1.66833
     1         7065273382d-1*t25*t28*t53*t68-8.293023216396499d-1*t25*t3
     2         5*t68-1.0748735402283824d+0*t25*t27*t34*t68-4.91825651378
     3         32956d+0*t25*t32*t33*t68+6.65559866912399d+0*t25*t30*t31*
     4         t68-9.504933078915558d+0*t23*t28*t61*t67-2.95744921182035
     5         94d+2*t23*t31*t53*t61-2.798161171838693d+0*t23*t34*t61+3.
     6         956348675640094d+2*t23*t30*t33*t61-8.758685213121911d+1*t
     7         23*t27*t32*t61)-9.305257363490997d-1*t23*t36*(1.150534517
     8         2712736d-7*gammabb**3*t66/rhob**10-1.6407925471300622d-4*
     9         t47*t48/rhob**7.333333333333333d+0+5.649276563507623d-2*g
     :         ammabb*t21/rhob**4.666666666666667d+0)-1.8610514726981994
     ;         d+0*t23*t51*t54-2.481401963597599d+0*t22*t46*t54-4.135669
     <         9393293317d-1*t22*t36/t52-2.481401963597599d+0*t36*t46*t5
     =         1)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t18*t5*(-4.3145
     1         044397672766d-8*t38*t63/rhoa**9-1.5407117900475334d-2*t3*
     2         t41+5.469308490433541d-5*gammaaa*t39*t40)-9.3052573634909
     3         97d-1*t45*t5*t56-1.2407009817987995d+0*t18*t37*t56)*wght+
     4         Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363490997d-1*t23*t36*(-4.314
     1         5044397672766d-8*t47*t66/rhob**9-1.5407117900475334d-2*t2
     2         1*t50+5.469308490433541d-5*gammabb*t48*t49)-9.30525736349
     3         0997d-1*t23*t54*t58-1.2407009817987995d+0*t36*t46*t58)*wg
     4         ht+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*t5*(1.6179391649127287d-8*gammaaa*t63/rhoa**8-1.
     2         3673271226083852d-5*t39*t55)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t23*t36*(1.6179391649127287d-8*gammabb*t66/rhob**8-1
     2         .3673271226083852d-5*t48*t57)*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*(5.2552111
     1         27873146d+1*t14*t69*t70*t9+1.6123103103425737d+0*t16*t43*
     2         t59*t9+5.702959847349336d+0*t10*t64*t69*t70+1.77446952709
     3         22156d+2*t13*t44*t69*t70+1.6788967031032173d+0*t16*t69*t7
     4         0-2.3738092053840557d+2*t12*t15*t69*t70-2.502505597910073
     5         d-1*t10*t43*t44*t59+1.2439534824594749d+0*t17*t43*t59+7.3
     6         77384770674944d+0*t14*t15*t43*t59-9.983398003685984d+0*t1
     7         2*t13*t43*t59)-9.305257363490997d-1*t42*t5*t60-1.24070098
     8         17987995d+0*t37*t4*t60)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t22*t23*(5.70295
     1         9847349336d+0*t28*t67*t71*t72+1.7744695270922156d+2*t31*t
     2         53*t71*t72+1.6788967031032173d+0*t34*t71*t72-2.3738092053
     3         840557d+2*t30*t33*t71*t72+5.255211127873146d+1*t27*t32*t7
     4         1*t72-2.502505597910073d-1*t28*t52*t53*t61+1.243953482459
     5         4749d+0*t35*t52*t61+1.6123103103425737d+0*t27*t34*t52*t61
     6         +7.377384770674944d+0*t32*t33*t52*t61-9.983398003685984d+
     7         0*t30*t31*t52*t61)-9.305257363490997d-1*t23*t51*t62-1.240
     8         7009817987995d+0*t22*t46*t62)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t56*t60*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t23*t58*t62*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-3.153126676723888d+1*t14*t73*t74*t9-1.9347723724
     2         110882d+0*t16*t6*t70*t9-3.4217759084096017d+0*t10*t64*t73
     3         *t74-1.0646817162553295d+2*t13*t44*t73*t74-1.007338021861
     4         9295d+0*t16*t73*t74+1.4242855232304333d+2*t12*t15*t73*t74
     5         +3.0030067174920877d-1*t10*t44*t6*t70-1.4927441789513698d
     6         +0*t17*t6*t70-8.852861724809932d+0*t14*t15*t6*t70+1.19800
     7         77604423181d+1*t12*t13*t6*t70)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t22*t23*(-3.4217759084096017d+0*t28*t67*t75*t76-1.064681
     2         7162553295d+2*t31*t53*t75*t76-1.0073380218619295d+0*t34*t
     3         75*t76+1.4242855232304333d+2*t30*t33*t75*t76-3.1531266767
     4         23888d+1*t27*t32*t75*t76+3.0030067174920877d-1*t24*t28*t5
     5         3*t72-1.4927441789513698d+0*t24*t35*t72-1.934772372411088
     6         2d+0*t24*t27*t34*t72-8.852861724809932d+0*t24*t32*t33*t72
     7         +1.1980077604423181d+1*t24*t30*t31*t72)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**2
            t21 = 1/t2**2
            t22 = 1/rhoa**6.333333333333333d+0
            t23 = 1/rhoa**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammaaa*t23*t3
            t25 = rhoa**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhoa**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t21*t28
            t30 = 1/taua**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t2**3
            t33 = 1/t11**6
            t34 = 1/t19
            t35 = rhoa**2.3333333333333334d+0
            t36 = 1/taua**3
            t37 = rhoa**3.3333333333333337d+0
            t38 = 1/taua**4
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t27*t4*t5*wght-9.3052
     1         57363490997d-1*t18*t24*t5*wght-1.2407009817987995d+0*t18*
     2         t19*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t4*t5*(-1.0748735
     1         402283824d+0*t16*t34*t7*t9-8.758685213121911d+1*t14*t30*t
     2         5*t9+1.668337065273382d-1*t10*t26*t34*t7-8.29302321639649
     3         9d-1*t17*t34*t7-4.9182565137832956d+0*t14*t15*t34*t7+6.65
     4         559866912399d+0*t12*t13*t34*t7-9.504933078915558d+0*t10*t
     5         30*t33*t5-2.9574492118203594d+2*t13*t26*t30*t5-2.79816117
     6         1838693d+0*t16*t30*t5+3.956348675640094d+2*t12*t15*t30*t5
     7         )*wght-9.305257363490997d-1*t18*(1.1505345172712736d-7*ga
     8         mmaaa**3*t32/rhoa**10+5.649276563507623d-2*gammaaa*t3/rho
     9         a**4.666666666666667d+0-1.6407925471300622d-4*t20*t21/rho
     :         a**7.333333333333333d+0)*t5*wght-1.8610514726981994d+0*t2
     ;         4*t27*t5*wght-2.481401963597599d+0*t19*t27*t4*wght-4.1356
     <         699393293317d-1*t18*t4*wght/t25-2.481401963597599d+0*t18*
     =         t19*t24*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t18*(-4.31450443
     1         97672766d-8*t20*t32/rhoa**9-1.5407117900475334d-2*t23*t3+
     2         5.469308490433541d-5*gammaaa*t21*t22)*t5*wght-9.305257363
     3         490997d-1*t27*t29*t5*wght-1.2407009817987995d+0*t18*t19*t
     4         29*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*(1.6179391649127287d-8*gammaaa*t32/rhoa**8-1.367
     2         3271226083852d-5*t21*t28)*t5*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t4*t5*(5.25521112
     1         7873146d+1*t14*t35*t36*t9+1.6123103103425737d+0*t16*t25*t
     2         30*t9+5.702959847349336d+0*t10*t33*t35*t36+1.774469527092
     3         2156d+2*t13*t26*t35*t36+1.6788967031032173d+0*t16*t35*t36
     4         -2.3738092053840557d+2*t12*t15*t35*t36-2.502505597910073d
     5         -1*t10*t25*t26*t30+1.2439534824594749d+0*t17*t25*t30+7.37
     6         7384770674944d+0*t14*t15*t25*t30-9.983398003685984d+0*t12
     7         *t13*t25*t30)*wght-9.305257363490997d-1*t24*t31*t5*wght-1
     8         .2407009817987995d+0*t19*t31*t4*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-1.9347723724110882d+0*t16*t36*t6*t9-3.1531266767
     2         23888d+1*t14*t37*t38*t9+3.0030067174920877d-1*t10*t26*t36
     3         *t6-1.4927441789513698d+0*t17*t36*t6-8.852861724809932d+0
     4         *t14*t15*t36*t6+1.1980077604423181d+1*t12*t13*t36*t6-3.42
     5         17759084096017d+0*t10*t33*t37*t38-1.0646817162553295d+2*t
     6         13*t26*t37*t38-1.0073380218619295d+0*t16*t37*t38+1.424285
     7         5232304333d+2*t12*t15*t37*t38)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = rhob**1.6666666666666669d+0
            t7 = 1/taub
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhob**3.333333333333333d-1
            t20 = gammabb**2
            t21 = 1/t2**2
            t22 = 1/rhob**6.333333333333333d+0
            t23 = 1/rhob**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammabb*t23*t3
            t25 = rhob**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhob**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammab
     1         b*t21*t28
            t30 = 1/taub**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t2**3
            t33 = 1/t11**6
            t34 = 1/t19
            t35 = rhob**2.3333333333333334d+0
            t36 = 1/taub**3
            t37 = rhob**3.3333333333333337d+0
            t38 = 1/taub**4
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t27*t4*t5*wght-9.3052
     1         57363490997d-1*t18*t24*t5*wght-1.2407009817987995d+0*t18*
     2         t19*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t4*t5*(-1.0748735
     1         402283824d+0*t16*t34*t7*t9-8.758685213121911d+1*t14*t30*t
     2         5*t9+1.668337065273382d-1*t10*t26*t34*t7-8.29302321639649
     3         9d-1*t17*t34*t7-4.9182565137832956d+0*t14*t15*t34*t7+6.65
     4         559866912399d+0*t12*t13*t34*t7-9.504933078915558d+0*t10*t
     5         30*t33*t5-2.9574492118203594d+2*t13*t26*t30*t5-2.79816117
     6         1838693d+0*t16*t30*t5+3.956348675640094d+2*t12*t15*t30*t5
     7         )*wght-9.305257363490997d-1*t18*(1.1505345172712736d-7*ga
     8         mmabb**3*t32/rhob**10+5.649276563507623d-2*gammabb*t3/rho
     9         b**4.666666666666667d+0-1.6407925471300622d-4*t20*t21/rho
     :         b**7.333333333333333d+0)*t5*wght-1.8610514726981994d+0*t2
     ;         4*t27*t5*wght-2.481401963597599d+0*t19*t27*t4*wght-4.1356
     <         699393293317d-1*t18*t4*wght/t25-2.481401963597599d+0*t18*
     =         t19*t24*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t18*(-4.31450443
     1         97672766d-8*t20*t32/rhob**9-1.5407117900475334d-2*t23*t3+
     2         5.469308490433541d-5*gammabb*t21*t22)*t5*wght-9.305257363
     3         490997d-1*t27*t29*t5*wght-1.2407009817987995d+0*t18*t19*t
     4         29*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t18*(1.6179391649127287d-8*gammabb*t32/rhob**8-1.367
     2         3271226083852d-5*t21*t28)*t5*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t4*t5*(5.25521112
     1         7873146d+1*t14*t35*t36*t9+1.6123103103425737d+0*t16*t25*t
     2         30*t9+5.702959847349336d+0*t10*t33*t35*t36+1.774469527092
     3         2156d+2*t13*t26*t35*t36+1.6788967031032173d+0*t16*t35*t36
     4         -2.3738092053840557d+2*t12*t15*t35*t36-2.502505597910073d
     5         -1*t10*t25*t26*t30+1.2439534824594749d+0*t17*t25*t30+7.37
     6         7384770674944d+0*t14*t15*t25*t30-9.983398003685984d+0*t12
     7         *t13*t25*t30)*wght-9.305257363490997d-1*t24*t31*t5*wght-1
     8         .2407009817987995d+0*t19*t31*t4*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t4*t5*(-1.9347723724110882d+0*t16*t36*t6*t9-3.1531266767
     2         23888d+1*t14*t37*t38*t9+3.0030067174920877d-1*t10*t26*t36
     3         *t6-1.4927441789513698d+0*t17*t36*t6-8.852861724809932d+0
     4         *t14*t15*t36*t6+1.1980077604423181d+1*t12*t13*t36*t6-3.42
     5         17759084096017d+0*t10*t33*t37*t38-1.0646817162553295d+2*t
     6         13*t26*t37*t38-1.0073380218619295d+0*t16*t37*t38+1.424285
     7         5232304333d+2*t12*t15*t37*t38)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_dldf functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_3} &=& {\it t_2}-1.0\\\\ 
C>   {\it t_4} &=& {\it t_2}+1.0\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_7} &=& {\it t_6}-1.0\\\\ 
C>   {\it t_8} &=& {\it t_6}+1.0\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{11}} &=& {\it t_{10}}-1.0\\\\ 
C>   {\it t_{12}} &=& {\it t_{10}}+1.0\\\\ 
C>   f &=& -0.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{0.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-0.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{0.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, 263201 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.103.263201 "> 
C> 10.1103/PhysRevLett.103.263201 </a> 
C> 
      subroutine nwxcm_x_dldf_d3(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**2
            t21 = 1/t2**2
            t22 = 1/rhoa**6.333333333333333d+0
            t23 = 1/rhoa**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammaaa*t23*t3
            t25 = rhoa**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhoa**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t21*t28
            t30 = 1/taua**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t25
            t33 = gammaaa**3
            t34 = 1/t2**3
            t35 = 1/rhoa**10
            t36 = 1/rhoa**7.333333333333333d+0
            t37 = 1/rhoa**4.666666666666667d+0
            t38 = 5.649276563507623d-2*gammaaa*t3*t37-1.6407925471300622
     1         d-4*t20*t21*t36+1.1505345172712736d-7*t33*t34*t35
            t39 = 1/t11**6
            t40 = 1/t19
            t41 = -1.0748735402283824d+0*t16*t40*t7*t9-8.758685213121911
     1         d+1*t14*t30*t5*t9+1.668337065273382d-1*t10*t26*t40*t7-8.2
     2         93023216396499d-1*t17*t40*t7-4.9182565137832956d+0*t14*t1
     3         5*t40*t7+6.65559866912399d+0*t12*t13*t40*t7-9.50493307891
     4         5558d+0*t10*t30*t39*t5-2.9574492118203594d+2*t13*t26*t30*
     5         t5-2.798161171838693d+0*t16*t30*t5+3.956348675640094d+2*t
     6         12*t15*t30*t5
            t42 = 1/rhoa**9
            t43 = -4.3145044397672766d-8*t20*t34*t42-1.5407117900475334d
     1         -2*t23*t3+5.469308490433541d-5*gammaaa*t21*t22
            t44 = 1/rhoa**8
            t45 = 1.6179391649127287d-8*gammaaa*t34*t44-1.36732712260838
     1         52d-5*t21*t28
            t46 = rhoa**2.3333333333333334d+0
            t47 = 1/taua**3
            t48 = 5.255211127873146d+1*t14*t46*t47*t9+1.6123103103425737
     1         d+0*t16*t25*t30*t9+5.702959847349336d+0*t10*t39*t46*t47+1
     2         .7744695270922156d+2*t13*t26*t46*t47+1.6788967031032173d+
     3         0*t16*t46*t47-2.3738092053840557d+2*t12*t15*t46*t47-2.502
     4         505597910073d-1*t10*t25*t26*t30+1.2439534824594749d+0*t17
     5         *t25*t30+7.377384770674944d+0*t14*t15*t25*t30-9.983398003
     6         685984d+0*t12*t13*t25*t30
            t49 = rhoa**3.3333333333333337d+0
            t50 = 1/taua**4
            t51 = 1/t2**4
            t52 = rhoa**2
            t53 = 1/t11**7
            t54 = 1/t5
            t55 = rhoa**3
            t56 = rhoa**5
            t57 = 1/taua**6
            t58 = 1/taua**5
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t18*t4*t5*wght
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t27*t4*t5-9.30525736
     1         3490997d-1*t18*t24*t5-1.2407009817987995d+0*t18*t19*t4)*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t41*t5-9.3052
     1         57363490997d-1*t18*t38*t5-1.8610514726981994d+0*t24*t27*t
     2         5-4.1356699393293317d-1*t18*t32*t4-2.481401963597599d+0*t
     3         19*t27*t4-2.481401963597599d+0*t18*t19*t24)*wght+Amat2(iq
     4         ,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t18*t43*t5-9.30
     1         5257363490997d-1*t27*t29*t5-1.2407009817987995d+0*t18*t19
     2         *t29)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*t45*t5*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t48*t5-9.3052
     1         57363490997d-1*t24*t31*t5-1.2407009817987995d+0*t19*t31*t
     2         4)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*(-1.9347723724110882d+0*t16*t47*t6*t9-3.1531266767
     2         23888d+1*t14*t49*t50*t9+3.0030067174920877d-1*t10*t26*t47
     3         *t6-1.4927441789513698d+0*t17*t47*t6-8.852861724809932d+0
     4         *t14*t15*t47*t6+1.1980077604423181d+1*t12*t13*t47*t6-3.42
     5         17759084096017d+0*t10*t39*t49*t50-1.0646817162553295d+2*t
     6         13*t26*t49*t50-1.0073380218619295d+0*t16*t49*t50+1.424285
     7         5232304333d+2*t12*t15*t49*t50)*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t4*t5*(3.5829
     1         11800761275d-1*t16*t54*t7*t9+8.006765213576478d+3*t12*t47
     2         *t52*t9-1.751737042624382d+2*t14*t19*t30*t9-5.56112355091
     3         1273d-2*t10*t26*t54*t7+2.764341072132166d-1*t17*t54*t7+1.
     4         639418837927765d+0*t14*t15*t54*t7-2.2185328897079964d+0*t
     5         12*t13*t54*t7+4.332158277373476d+2*t10*t47*t52*t53+1.0944
     6         07414826128d+4*t13*t39*t47*t52-1.876122781257643d+4*t15*t
     7         26*t47*t52-6.228273769986757d+2*t14*t47*t52-1.90098661578
     8         31115d+1*t10*t19*t30*t39-5.914898423640719d+2*t13*t19*t26
     9         *t30-5.596322343677388d+0*t16*t19*t30+7.912697351280188d+
     :         2*t12*t15*t19*t30)+2.757113292886221d-1*t18*t4/t6-9.30525
     ;         7363490997d-1*t18*t5*(1.0891291924501978d-9*gammaaa**4*t5
     <         1/rhoa**1.3666666666666666d+1-2.1860155828154199d-6*t33*t
     =         34/rhoa**11-2.636329062970224d-1*gammaaa*t3/rhoa**5.66666
     >         6666666667d+0+1.381506811287287d-3*t20*t21/rhoa**8.333333
     ?         333333334d+0)-2.791577209047299d+0*t24*t41*t5-2.791577209
     @         047299d+0*t27*t38*t5-3.7221029453963983d+0*t19*t4*t41+t32
     1         *(-1.2407009817987995d+0*t27*t4-1.2407009817987995d+0*t18
     2         *t24)-3.7221029453963983d+0*t18*t19*t38-7.444205890792797
     3         d+0*t19*t24*t27)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363490997d-1*t18*t5*(-4.0
     1         84234471688241d-10*t33*t51/rhoa**1.2666666666666665d+1+5.
     2         649276563507623d-2*t3*t37-3.9500561319797795d-4*gammaaa*t
     3         21*t36+7.334657547604368d-7*t20*t34*t35)-1.86105147269819
     4         94d+0*t27*t43*t5-9.305257363490997d-1*t29*t41*t5-2.481401
     5         963597599d+0*t18*t19*t43-4.1356699393293317d-1*t18*t29*t3
     6         2-2.481401963597599d+0*t19*t27*t29)*wght+Cmat3(iq,D3_RA_R
     7         A_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t18*t5*(1.5
     1         315879268830906d-10*t20*t51/rhoa**1.1666666666666665d+1-2
     2         .1572522198836383d-7*gammaaa*t34*t42+7.292411320578055d-5
     3         *t21*t22)-9.305257363490997d-1*t27*t45*t5-1.2407009817987
     4         995d+0*t18*t19*t45)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t18*t5*(4.8538174947381857d-8*t34*t44-5.7434
     2         5472581159d-11*gammaaa*t51/rhoa**1.0666666666666666d+1)*w
     3         ght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(-4.804
     1         0591281458866d+3*t12*t50*t55*t9+2.1020844511492584d+2*t14
     2         *t47*t5*t9+1.0748735402283824d+0*t16*t30*t40*t9-2.5992949
     3         66424086d+2*t10*t50*t53*t55-6.56644448895677d+3*t13*t39*t
     4         50*t55+1.125673668754586d+4*t15*t26*t50*t55+3.73696426199
     5         2054d+2*t14*t50*t55+2.2811839389397337d+1*t10*t39*t47*t5+
     6         7.097878108368861d+2*t13*t26*t47*t5+6.715586812412864d+0*
     7         t16*t47*t5-9.495236821536224d+2*t12*t15*t47*t5-1.66833706
     8         5273382d-1*t10*t26*t30*t40+8.293023216396499d-1*t17*t30*t
     9         40+4.9182565137832956d+0*t14*t15*t30*t40-6.65559866912399
     :         d+0*t12*t13*t30*t40)-1.8610514726981994d+0*t24*t48*t5-9.3
     ;         05257363490997d-1*t31*t38*t5-2.481401963597599d+0*t19*t4*
     <         t48-4.1356699393293317d-1*t31*t32*t4-2.481401963597599d+0
     =         *t19*t24*t31)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_GAA_TA) = (-9.305257363490997d-1*t29*t48*t5-9
     1         .305257363490997d-1*t31*t43*t5-1.2407009817987995d+0*t19*
     2         t29*t31)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t31*t45*t5*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(5.804317117233265d+0*t16*t50*t6*t9+1.891876
     2         006034333d+2*t14*t49*t58*t9-1.7294612861325187d+3*t12*t56
     3         *t57*t9-9.009020152476264d-1*t10*t26*t50*t6+4.47823253685
     4         411d+0*t17*t50*t6+2.6558585174429794d+1*t14*t15*t50*t6-3.
     5         5940232813269546d+1*t12*t13*t50*t6+2.0530655450457608d+1*
     6         t10*t39*t49*t58+6.388090297531976d+2*t13*t26*t49*t58+6.04
     7         4028131171577d+0*t16*t49*t58-8.5457131393826d+2*t12*t15*t
     8         49*t58-9.35746187912671d+1*t10*t53*t56*t57-2.363920016024
     9         4372d+3*t13*t39*t56*t57+4.0524252075165085d+3*t15*t26*t56
     :         *t57+1.3453071343171397d+2*t14*t56*t57)*wght
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
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = 1.1832860901832055d-3*gammabb*t19+1.0d+0
            t21 = 1/t20
            t22 = 5.777669212678251d-3*gammabb*t19*t21+1.0d+0
            t23 = rhob**1.3333333333333333d+0
            t24 = rhob**1.6666666666666669d+0
            t25 = 1/taub
            t26 = 4.557799872345596d+0*t24*t25
            t27 = t26-1.0d+0
            t28 = t27**4
            t29 = t26+1.0d+0
            t30 = 1/t29**4
            t31 = t27**3
            t32 = 1/t29**3
            t33 = t27**2
            t34 = 1/t29**2
            t35 = 1/t29
            t36 = 1.0d+0*(-1.637571d-1*t27*t35-1.880028d-1*t33*t34-4.490
     1         609d-1*t31*t32-8.2359d-3*t28*t30)+1.0d+0
            t37 = rhoa**3.333333333333333d-1
            t38 = gammaaa**2
            t39 = 1/t2**2
            t40 = 1/rhoa**6.333333333333333d+0
            t41 = 1/rhoa**3.6666666666666664d+0
            t42 = 1.8231028301445137d-5*t38*t39*t40-1.5407117900475334d-
     1         2*gammaaa*t3*t41
            t43 = rhoa**6.666666666666666d-1
            t44 = 1/t11**5
            t45 = -1.6123103103425737d+0*t16*t43*t7*t9+2.502505597910073
     1         d-1*t10*t43*t44*t7-1.2439534824594749d+0*t17*t43*t7-7.377
     2         384770674944d+0*t14*t15*t43*t7+9.983398003685984d+0*t12*t
     3         13*t43*t7
            t46 = rhob**3.333333333333333d-1
            t47 = gammabb**2
            t48 = 1/t20**2
            t49 = 1/rhob**6.333333333333333d+0
            t50 = 1/rhob**3.6666666666666664d+0
            t51 = 1.8231028301445137d-5*t47*t48*t49-1.5407117900475334d-
     1         2*gammabb*t21*t50
            t52 = rhob**6.666666666666666d-1
            t53 = 1/t29**5
            t54 = 2.502505597910073d-1*t25*t28*t52*t53-1.243953482459474
     1         9d+0*t25*t35*t52-1.6123103103425737d+0*t25*t27*t34*t52-7.
     2         377384770674944d+0*t25*t32*t33*t52+9.983398003685984d+0*t
     3         25*t30*t31*t52
            t55 = 1/rhoa**5.333333333333333d+0
            t56 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t39*t55
            t57 = 1/rhob**5.333333333333333d+0
            t58 = 5.777669212678251d-3*t19*t21-6.836635613041927d-6*gamm
     1         abb*t48*t57
            t59 = 1/taua**2
            t60 = 9.673861862055442d-1*t16*t59*t6*t9-1.5015033587460438d
     1         -1*t10*t44*t59*t6+7.463720894756849d-1*t17*t59*t6+4.42643
     2         0862404966d+0*t14*t15*t59*t6-5.99003880221159d+0*t12*t13*
     3         t59*t6
            t61 = 1/taub**2
            t62 = -1.5015033587460438d-1*t24*t28*t53*t61+7.4637208947568
     1         49d-1*t24*t35*t61+9.673861862055442d-1*t24*t27*t34*t61+4.
     2         426430862404966d+0*t24*t32*t33*t61-5.99003880221159d+0*t2
     3         4*t30*t31*t61
            t63 = 1/t43
            t64 = gammaaa**3
            t65 = 1/t2**3
            t66 = 1/rhoa**10
            t67 = 1/rhoa**7.333333333333333d+0
            t68 = 1/rhoa**4.666666666666667d+0
            t69 = 5.649276563507623d-2*gammaaa*t3*t68-1.6407925471300622
     1         d-4*t38*t39*t67+1.1505345172712736d-7*t64*t65*t66
            t70 = 1/t11**6
            t71 = 1/t37
            t72 = -1.0748735402283824d+0*t16*t7*t71*t9-8.758685213121911
     1         d+1*t14*t5*t59*t9+1.668337065273382d-1*t10*t44*t7*t71-8.2
     2         93023216396499d-1*t17*t7*t71-4.9182565137832956d+0*t14*t1
     3         5*t7*t71+6.65559866912399d+0*t12*t13*t7*t71-9.50493307891
     4         5558d+0*t10*t5*t59*t70-2.9574492118203594d+2*t13*t44*t5*t
     5         59-2.798161171838693d+0*t16*t5*t59+3.956348675640094d+2*t
     6         12*t15*t5*t59
            t73 = 1/t52
            t74 = gammabb**3
            t75 = 1/t20**3
            t76 = 1/rhob**10
            t77 = 1/rhob**7.333333333333333d+0
            t78 = 1/rhob**4.666666666666667d+0
            t79 = 5.649276563507623d-2*gammabb*t21*t78-1.640792547130062
     1         2d-4*t47*t48*t77+1.1505345172712736d-7*t74*t75*t76
            t80 = 1/t29**6
            t81 = 1/t46
            t82 = 1.668337065273382d-1*t25*t28*t53*t81-8.293023216396499
     1         d-1*t25*t35*t81-1.0748735402283824d+0*t25*t27*t34*t81-4.9
     2         182565137832956d+0*t25*t32*t33*t81+6.65559866912399d+0*t2
     3         5*t30*t31*t81-9.504933078915558d+0*t23*t28*t61*t80-2.9574
     4         492118203594d+2*t23*t31*t53*t61-2.798161171838693d+0*t23*
     5         t34*t61+3.956348675640094d+2*t23*t30*t33*t61-8.7586852131
     6         21911d+1*t23*t27*t32*t61
            t83 = 1/rhoa**9
            t84 = -4.3145044397672766d-8*t38*t65*t83-1.5407117900475334d
     1         -2*t3*t41+5.469308490433541d-5*gammaaa*t39*t40
            t85 = 1/rhob**9
            t86 = -4.3145044397672766d-8*t47*t75*t85-1.5407117900475334d
     1         -2*t21*t50+5.469308490433541d-5*gammabb*t48*t49
            t87 = 1/rhoa**8
            t88 = 1.6179391649127287d-8*gammaaa*t65*t87-1.36732712260838
     1         52d-5*t39*t55
            t89 = 1/rhob**8
            t90 = 1.6179391649127287d-8*gammabb*t75*t89-1.36732712260838
     1         52d-5*t48*t57
            t91 = rhoa**2.3333333333333334d+0
            t92 = 1/taua**3
            t93 = 5.255211127873146d+1*t14*t9*t91*t92+5.702959847349336d
     1         +0*t10*t70*t91*t92+1.7744695270922156d+2*t13*t44*t91*t92+
     2         1.6788967031032173d+0*t16*t91*t92-2.3738092053840557d+2*t
     3         12*t15*t91*t92+1.6123103103425737d+0*t16*t43*t59*t9-2.502
     4         505597910073d-1*t10*t43*t44*t59+1.2439534824594749d+0*t17
     5         *t43*t59+7.377384770674944d+0*t14*t15*t43*t59-9.983398003
     6         685984d+0*t12*t13*t43*t59
            t94 = rhob**2.3333333333333334d+0
            t95 = 1/taub**3
            t96 = 5.702959847349336d+0*t28*t80*t94*t95+1.774469527092215
     1         6d+2*t31*t53*t94*t95+1.6788967031032173d+0*t34*t94*t95-2.
     2         3738092053840557d+2*t30*t33*t94*t95+5.255211127873146d+1*
     3         t27*t32*t94*t95-2.502505597910073d-1*t28*t52*t53*t61+1.24
     4         39534824594749d+0*t35*t52*t61+1.6123103103425737d+0*t27*t
     5         34*t52*t61+7.377384770674944d+0*t32*t33*t52*t61-9.9833980
     6         03685984d+0*t30*t31*t52*t61
            t97 = rhoa**3.3333333333333337d+0
            t98 = 1/taua**4
            t99 = -3.153126676723888d+1*t14*t9*t97*t98-3.421775908409601
     1         7d+0*t10*t70*t97*t98-1.0646817162553295d+2*t13*t44*t97*t9
     2         8-1.0073380218619295d+0*t16*t97*t98+1.4242855232304333d+2
     3         *t12*t15*t97*t98-1.9347723724110882d+0*t16*t6*t9*t92+3.00
     4         30067174920877d-1*t10*t44*t6*t92-1.4927441789513698d+0*t1
     5         7*t6*t92-8.852861724809932d+0*t14*t15*t6*t92+1.1980077604
     6         423181d+1*t12*t13*t6*t92
            t100 = rhob**3.3333333333333337d+0
            t101 = 1/taub**4
            t102 = 3.0030067174920877d-1*t24*t28*t53*t95-1.4927441789513
     1         698d+0*t24*t35*t95-1.9347723724110882d+0*t24*t27*t34*t95-
     2         8.852861724809932d+0*t24*t32*t33*t95+1.1980077604423181d+
     3         1*t24*t30*t31*t95-3.4217759084096017d+0*t100*t101*t28*t80
     4         -1.0646817162553295d+2*t100*t101*t31*t53-1.00733802186192
     5         95d+0*t100*t101*t34+1.4242855232304333d+2*t100*t101*t30*t
     6         33-3.153126676723888d+1*t100*t101*t27*t32
            t103 = 1/t2**4
            t104 = rhoa**2
            t105 = 1/t11**7
            t106 = 1/t5
            t107 = 1/t20**4
            t108 = rhob**2
            t109 = 1/t29**7
            t110 = 1/t23
            t111 = rhoa**3
            t112 = rhob**3
            t113 = rhoa**5
            t114 = 1/taua**6
            t115 = 1/taua**5
            t116 = rhob**5
            t117 = 1/taub**6
            t118 = 1/taub**5
            fnc(iq) = (-9.305257363490997d-1*t18*t4*t5-9.305257363490997
     1         d-1*t22*t23*t36)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t45*t5-9.30525736
     1         3490997d-1*t18*t42*t5-1.2407009817987995d+0*t18*t37*t4)*w
     2         ght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t22*t23*t54-9.305257
     1         363490997d-1*t23*t36*t51-1.2407009817987995d+0*t22*t36*t4
     2         6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         5*t56*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t23*t
     1         36*t58*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*t
     1         60*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t22*t23
     1         *t62*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363490997d-1*t4*t5*t72-9.3052
     1         57363490997d-1*t18*t5*t69-4.1356699393293317d-1*t18*t4*t6
     2         3-1.8610514726981994d+0*t42*t45*t5-2.481401963597599d+0*t
     3         37*t4*t45-2.481401963597599d+0*t18*t37*t42)*wght+Amat2(iq
     4         ,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363490997d-1*t22*t23*t82-9.30
     1         5257363490997d-1*t23*t36*t79-4.1356699393293317d-1*t22*t3
     2         6*t73-1.8610514726981994d+0*t23*t51*t54-2.481401963597599
     3         d+0*t22*t46*t54-2.481401963597599d+0*t36*t46*t51)*wght+Am
     4         at2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363490997d-1*t18*t5*t84-9.30
     1         5257363490997d-1*t45*t5*t56-1.2407009817987995d+0*t18*t37
     2         *t56)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363490997d-1*t23*t36*t86-9.3
     1         05257363490997d-1*t23*t54*t58-1.2407009817987995d+0*t36*t
     2         46*t58)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*t5*t88*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t23*t36*t90*wght
            Mmat2(iq,D2_RA_TA) = (-9.305257363490997d-1*t4*t5*t93-9.3052
     1         57363490997d-1*t42*t5*t60-1.2407009817987995d+0*t37*t4*t6
     2         0)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (-9.305257363490997d-1*t22*t23*t96-9.30
     1         5257363490997d-1*t23*t51*t62-1.2407009817987995d+0*t22*t4
     2         6*t62)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t5*t56*t60*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t23*t58*t62*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*t99*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t102*t22*t23*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363490997d-1*t4*t5*(8.0067
     1         65213576478d+3*t104*t12*t9*t92+1.094407414826128d+4*t104*
     2         t13*t70*t92-1.876122781257643d+4*t104*t15*t44*t92-6.22827
     3         3769986757d+2*t104*t14*t92+4.332158277373476d+2*t10*t104*
     4         t105*t92+3.582911800761275d-1*t106*t16*t7*t9-1.7517370426
     5         24382d+2*t14*t37*t59*t9-1.9009866157831115d+1*t10*t37*t59
     6         *t70-5.561123550911273d-2*t10*t106*t44*t7+2.7643410721321
     7         66d-1*t106*t17*t7+1.639418837927765d+0*t106*t14*t15*t7-2.
     8         2185328897079964d+0*t106*t12*t13*t7-5.914898423640719d+2*
     9         t13*t37*t44*t59-5.596322343677388d+0*t16*t37*t59+7.912697
     :         351280188d+2*t12*t15*t37*t59)-2.791577209047299d+0*t42*t5
     ;         *t72-3.7221029453963983d+0*t37*t4*t72-2.791577209047299d+
     <         0*t45*t5*t69-3.7221029453963983d+0*t18*t37*t69-9.30525736
     =         3490997d-1*t18*t5*(-2.1860155828154199d-6*t64*t65/rhoa**1
     >         1+1.381506811287287d-3*t38*t39/rhoa**8.333333333333334d+0
     ?         -2.636329062970224d-1*gammaaa*t3/rhoa**5.666666666666667d
     @         +0+1.0891291924501978d-9*gammaaa**4*t103/rhoa**1.36666666
     1         66666666d+1)+(-1.2407009817987995d+0*t4*t45-1.24070098179
     2         87995d+0*t18*t42)*t63+2.757113292886221d-1*t18*t4/t6-7.44
     3         4205890792797d+0*t37*t42*t45)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363490997d-1*t22*t23*(1.09
     1         4407414826128d+4*t108*t31*t80*t95-1.876122781257643d+4*t1
     2         08*t33*t53*t95-6.228273769986757d+2*t108*t32*t95+8.006765
     3         213576478d+3*t108*t27*t30*t95+4.332158277373476d+2*t108*t
     4         109*t28*t95-1.9009866157831115d+1*t28*t46*t61*t80-5.91489
     5         8423640719d+2*t31*t46*t53*t61-5.596322343677388d+0*t34*t4
     6         6*t61+7.912697351280188d+2*t30*t33*t46*t61-1.751737042624
     7         382d+2*t27*t32*t46*t61-5.561123550911273d-2*t110*t25*t28*
     8         t53+2.764341072132166d-1*t110*t25*t35+3.582911800761275d-
     9         1*t110*t25*t27*t34+1.639418837927765d+0*t110*t25*t32*t33-
     :         2.2185328897079964d+0*t110*t25*t30*t31)-2.791577209047299
     ;         d+0*t23*t51*t82-3.7221029453963983d+0*t22*t46*t82-2.79157
     <         7209047299d+0*t23*t54*t79-3.7221029453963983d+0*t36*t46*t
     =         79-9.305257363490997d-1*t23*t36*(-2.1860155828154199d-6*t
     >         74*t75/rhob**11+1.381506811287287d-3*t47*t48/rhob**8.3333
     ?         33333333334d+0-2.636329062970224d-1*gammabb*t21/rhob**5.6
     @         66666666666667d+0+1.0891291924501978d-9*gammabb**4*t107/r
     1         hob**1.3666666666666666d+1)+(-1.2407009817987995d+0*t22*t
     2         54-1.2407009817987995d+0*t36*t51)*t73-7.444205890792797d+
     3         0*t46*t51*t54+2.757113292886221d-1*t22*t36/t24)*wght+Amat
     4         3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-1.8610514726981994d+0*t45*t5*t84-
     1         2.481401963597599d+0*t18*t37*t84-9.305257363490997d-1*t5*
     2         t56*t72-9.305257363490997d-1*t18*t5*(5.649276563507623d-2
     3         *t3*t68-3.9500561319797795d-4*gammaaa*t39*t67+7.334657547
     4         604368d-7*t38*t65*t66-4.084234471688241d-10*t103*t64/rhoa
     5         **1.2666666666666665d+1)-4.1356699393293317d-1*t18*t56*t6
     6         3-2.481401963597599d+0*t37*t45*t56)*wght+Cmat3(iq,D3_RA_R
     7         A_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-1.8610514726981994d+0*t23*t54*t86
     1         -2.481401963597599d+0*t36*t46*t86-9.305257363490997d-1*t2
     2         3*t58*t82-9.305257363490997d-1*t23*t36*(5.649276563507623
     3         d-2*t21*t78-3.9500561319797795d-4*gammabb*t48*t77+7.33465
     4         7547604368d-7*t47*t75*t76-4.084234471688241d-10*t107*t74/
     5         rhob**1.2666666666666665d+1)-4.1356699393293317d-1*t36*t5
     6         8*t73-2.481401963597599d+0*t46*t54*t58)*wght+Cmat3(iq,D3_
     7         RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363490997d-1*t45*t5*t88-
     1         1.2407009817987995d+0*t18*t37*t88-9.305257363490997d-1*t1
     2         8*t5*(-2.1572522198836383d-7*gammaaa*t65*t83+7.2924113205
     3         78055d-5*t39*t40+1.5315879268830906d-10*t103*t38/rhoa**1.
     4         1666666666666665d+1))*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363490997d-1*t23*t54*t90
     1         -1.2407009817987995d+0*t36*t46*t90-9.305257363490997d-1*t
     2         23*t36*(-2.1572522198836383d-7*gammabb*t75*t85+7.29241132
     3         0578055d-5*t48*t49+1.5315879268830906d-10*t107*t47/rhob**
     4         1.1666666666666665d+1))*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t18*t5*(4.8538174947381857d-8*t65*t87-5.7434
     2         5472581159d-11*gammaaa*t103/rhoa**1.0666666666666666d+1)*
     3         wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*t23*t36*(4.8538174947381857d-8*t75*t89-5.743
     2         45472581159d-11*gammabb*t107/rhob**1.0666666666666666d+1)
     3         *wght
            Mmat3(iq,D3_RA_RA_TA) = (-9.305257363490997d-1*t4*t5*(-4.804
     1         0591281458866d+3*t111*t12*t9*t98-6.56644448895677d+3*t111
     2         *t13*t70*t98+1.125673668754586d+4*t111*t15*t44*t98+3.7369
     3         64261992054d+2*t111*t14*t98-2.599294966424086d+2*t10*t105
     4         *t111*t98+2.1020844511492584d+2*t14*t5*t9*t92+2.281183938
     5         9397337d+1*t10*t5*t70*t92+7.097878108368861d+2*t13*t44*t5
     6         *t92+6.715586812412864d+0*t16*t5*t92-9.495236821536224d+2
     7         *t12*t15*t5*t92+1.0748735402283824d+0*t16*t59*t71*t9-1.66
     8         8337065273382d-1*t10*t44*t59*t71+8.293023216396499d-1*t17
     9         *t59*t71+4.9182565137832956d+0*t14*t15*t59*t71-6.65559866
     :         912399d+0*t12*t13*t59*t71)-1.8610514726981994d+0*t42*t5*t
     ;         93-2.481401963597599d+0*t37*t4*t93-9.305257363490997d-1*t
     <         5*t60*t69-4.1356699393293317d-1*t4*t60*t63-2.481401963597
     =         599d+0*t37*t42*t60)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (-1.8610514726981994d+0*t23*t51*t96-
     1         2.481401963597599d+0*t22*t46*t96-9.305257363490997d-1*t22
     2         *t23*(2.2811839389397337d+1*t23*t28*t80*t95+7.09787810836
     3         8861d+2*t23*t31*t53*t95+6.715586812412864d+0*t23*t34*t95-
     4         9.495236821536224d+2*t23*t30*t33*t95+2.1020844511492584d+
     5         2*t23*t27*t32*t95-1.668337065273382d-1*t28*t53*t61*t81+8.
     6         293023216396499d-1*t35*t61*t81+1.0748735402283824d+0*t27*
     7         t34*t61*t81+4.9182565137832956d+0*t32*t33*t61*t81-6.65559
     8         866912399d+0*t30*t31*t61*t81-6.56644448895677d+3*t101*t11
     9         2*t31*t80+1.125673668754586d+4*t101*t112*t33*t53+3.736964
     :         261992054d+2*t101*t112*t32-4.8040591281458866d+3*t101*t11
     ;         2*t27*t30-2.599294966424086d+2*t101*t109*t112*t28)-9.3052
     <         57363490997d-1*t23*t62*t79-4.1356699393293317d-1*t22*t62*
     =         t73-2.481401963597599d+0*t46*t51*t62)*wght+Mmat3(iq,D3_RB
     >         _RB_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (-9.305257363490997d-1*t5*t56*t93-9
     1         .305257363490997d-1*t5*t60*t84-1.2407009817987995d+0*t37*
     2         t56*t60)*wght+Mmat3(iq,D3_RA_GAA_TA)
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
            Mmat3(iq,D3_RB_GBB_TB) = (-9.305257363490997d-1*t23*t58*t96-
     1         9.305257363490997d-1*t23*t62*t86-1.2407009817987995d+0*t4
     2         6*t58*t62)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = Mmat3(iq,D3_RA_TA_TA)-9.305257363490
     1         997d-1*t5*t60*t88*wght
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
     1         3490997d-1*t23*t62*t90*wght
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)-9.30525736
     1         3490997d-1*t5*t56*t99*wght
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)-9.3052573634
     1         90997d-1*t102*t23*t58*wght
            Mmat3(iq,D3_GAB_TA_TA) = Mmat3(iq,D3_GAB_TA_TA)-9.3052573634
     1         90997d-1*t4*t5*(5.804317117233265d+0*t16*t6*t9*t98-9.0090
     2         20152476264d-1*t10*t44*t6*t98+4.47823253685411d+0*t17*t6*
     3         t98+2.6558585174429794d+1*t14*t15*t6*t98-3.59402328132695
     4         46d+1*t12*t13*t6*t98+1.891876006034333d+2*t115*t14*t9*t97
     5         +2.0530655450457608d+1*t10*t115*t70*t97+6.388090297531976
     6         d+2*t115*t13*t44*t97+6.044028131171577d+0*t115*t16*t97-8.
     7         5457131393826d+2*t115*t12*t15*t97-1.7294612861325187d+3*t
     8         113*t114*t12*t9-2.3639200160244372d+3*t113*t114*t13*t70+4
     9         .0524252075165085d+3*t113*t114*t15*t44+1.3453071343171397
     :         d+2*t113*t114*t14-9.35746187912671d+1*t10*t105*t113*t114)
     ;         *wght
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)-9.3052573634
     1         90997d-1*t22*t23*(-2.3639200160244372d+3*t116*t117*t31*t8
     2         0+2.0530655450457608d+1*t100*t118*t28*t80+4.0524252075165
     3         085d+3*t116*t117*t33*t53+6.388090297531976d+2*t100*t118*t
     4         31*t53-9.009020152476264d-1*t101*t24*t28*t53+4.4782325368
     5         5411d+0*t101*t24*t35+5.804317117233265d+0*t101*t24*t27*t3
     6         4+6.044028131171577d+0*t100*t118*t34+2.6558585174429794d+
     7         1*t101*t24*t32*t33-8.5457131393826d+2*t100*t118*t30*t33+1
     8         .891876006034333d+2*t100*t118*t27*t32+1.3453071343171397d
     9         +2*t116*t117*t32-3.5940232813269546d+1*t101*t24*t30*t31-1
     :         .7294612861325187d+3*t116*t117*t27*t30-9.35746187912671d+
     ;         1*t109*t116*t117*t28)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = rhoa**1.6666666666666669d+0
            t7 = 1/taua
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhoa**3.333333333333333d-1
            t20 = gammaaa**2
            t21 = 1/t2**2
            t22 = 1/rhoa**6.333333333333333d+0
            t23 = 1/rhoa**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammaaa*t23*t3
            t25 = rhoa**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhoa**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammaa
     1         a*t21*t28
            t30 = 1/taua**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t25
            t33 = gammaaa**3
            t34 = 1/t2**3
            t35 = 1/rhoa**10
            t36 = 1/rhoa**7.333333333333333d+0
            t37 = 1/rhoa**4.666666666666667d+0
            t38 = 5.649276563507623d-2*gammaaa*t3*t37-1.6407925471300622
     1         d-4*t20*t21*t36+1.1505345172712736d-7*t33*t34*t35
            t39 = 1/t11**6
            t40 = 1/t19
            t41 = -1.0748735402283824d+0*t16*t40*t7*t9-8.758685213121911
     1         d+1*t14*t30*t5*t9+1.668337065273382d-1*t10*t26*t40*t7-8.2
     2         93023216396499d-1*t17*t40*t7-4.9182565137832956d+0*t14*t1
     3         5*t40*t7+6.65559866912399d+0*t12*t13*t40*t7-9.50493307891
     4         5558d+0*t10*t30*t39*t5-2.9574492118203594d+2*t13*t26*t30*
     5         t5-2.798161171838693d+0*t16*t30*t5+3.956348675640094d+2*t
     6         12*t15*t30*t5
            t42 = 1/rhoa**9
            t43 = -4.3145044397672766d-8*t20*t34*t42-1.5407117900475334d
     1         -2*t23*t3+5.469308490433541d-5*gammaaa*t21*t22
            t44 = 1/rhoa**8
            t45 = 1.6179391649127287d-8*gammaaa*t34*t44-1.36732712260838
     1         52d-5*t21*t28
            t46 = rhoa**2.3333333333333334d+0
            t47 = 1/taua**3
            t48 = 5.255211127873146d+1*t14*t46*t47*t9+1.6123103103425737
     1         d+0*t16*t25*t30*t9+5.702959847349336d+0*t10*t39*t46*t47+1
     2         .7744695270922156d+2*t13*t26*t46*t47+1.6788967031032173d+
     3         0*t16*t46*t47-2.3738092053840557d+2*t12*t15*t46*t47-2.502
     4         505597910073d-1*t10*t25*t26*t30+1.2439534824594749d+0*t17
     5         *t25*t30+7.377384770674944d+0*t14*t15*t25*t30-9.983398003
     6         685984d+0*t12*t13*t25*t30
            t49 = rhoa**3.3333333333333337d+0
            t50 = 1/taua**4
            t51 = -1.9347723724110882d+0*t16*t47*t6*t9-3.153126676723888
     1         d+1*t14*t49*t50*t9+3.0030067174920877d-1*t10*t26*t47*t6-1
     2         .4927441789513698d+0*t17*t47*t6-8.852861724809932d+0*t14*
     3         t15*t47*t6+1.1980077604423181d+1*t12*t13*t47*t6-3.4217759
     4         084096017d+0*t10*t39*t49*t50-1.0646817162553295d+2*t13*t2
     5         6*t49*t50-1.0073380218619295d+0*t16*t49*t50+1.42428552323
     6         04333d+2*t12*t15*t49*t50
            t52 = 1/t2**4
            t53 = rhoa**2
            t54 = 1/t11**7
            t55 = 1/t5
            t56 = rhoa**3
            t57 = rhoa**4
            t58 = 1/taua**5
            t59 = rhoa**5
            t60 = 1/taua**6
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t27*t4*t5*wght-9.3052
     1         57363490997d-1*t18*t24*t5*wght-1.2407009817987995d+0*t18*
     2         t19*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = -9.305257363490997d-1*t4*t41*t5*wght-9.
     1         305257363490997d-1*t18*t38*t5*wght-1.8610514726981994d+0*
     2         t24*t27*t5*wght-4.1356699393293317d-1*t18*t32*t4*wght-2.4
     3         81401963597599d+0*t19*t27*t4*wght-2.481401963597599d+0*t1
     4         8*t19*t24*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490997d-1*t18*t43*t5*wght-
     1         9.305257363490997d-1*t27*t29*t5*wght-1.2407009817987995d+
     2         0*t18*t19*t29*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         7d-1*t18*t45*t5*wght
            Mmat2(iq,D2_RA_TA) = -9.305257363490997d-1*t4*t48*t5*wght-9.
     1         305257363490997d-1*t24*t31*t5*wght-1.2407009817987995d+0*
     2         t19*t31*t4*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = Mmat2(iq,D2_GAA_TA)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)-9.305257363490997d-1
     1         *t4*t5*t51*wght
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = t32*(-1.2407009817987995d+0*t27*t4*w
     1         ght-1.2407009817987995d+0*t18*t24*wght)-9.305257363490997
     2         d-1*t4*t5*(3.582911800761275d-1*t16*t55*t7*t9+8.006765213
     3         576478d+3*t12*t47*t53*t9-1.751737042624382d+2*t14*t19*t30
     4         *t9-5.561123550911273d-2*t10*t26*t55*t7+2.764341072132166
     5         d-1*t17*t55*t7+1.639418837927765d+0*t14*t15*t55*t7-2.2185
     6         328897079964d+0*t12*t13*t55*t7+4.332158277373476d+2*t10*t
     7         47*t53*t54+1.094407414826128d+4*t13*t39*t47*t53-1.8761227
     8         81257643d+4*t15*t26*t47*t53-6.228273769986757d+2*t14*t47*
     9         t53-1.9009866157831115d+1*t10*t19*t30*t39-5.9148984236407
     :         19d+2*t13*t19*t26*t30-5.596322343677388d+0*t16*t19*t30+7.
     ;         912697351280188d+2*t12*t15*t19*t30)*wght+2.75711329288622
     <         1d-1*t18*t4*wght/t6-9.305257363490997d-1*t18*t5*(1.089129
     =         1924501978d-9*gammaaa**4*t52/rhoa**1.3666666666666666d+1-
     >         2.1860155828154199d-6*t33*t34/rhoa**11-2.636329062970224d
     ?         -1*gammaaa*t3/rhoa**5.666666666666667d+0+1.38150681128728
     @         7d-3*t20*t21/rhoa**8.333333333333334d+0)*wght-2.791577209
     1         047299d+0*t24*t41*t5*wght-2.791577209047299d+0*t27*t38*t5
     2         *wght-3.7221029453963983d+0*t19*t4*t41*wght-3.72210294539
     3         63983d+0*t18*t19*t38*wght-7.444205890792797d+0*t19*t24*t2
     4         7*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363490997d-1*t18*t5*(-4.08
     1         4234471688241d-10*t33*t52/rhoa**1.2666666666666665d+1+5.6
     2         49276563507623d-2*t3*t37-3.9500561319797795d-4*gammaaa*t2
     3         1*t36+7.334657547604368d-7*t20*t34*t35)*wght-1.8610514726
     4         981994d+0*t27*t43*t5*wght-9.305257363490997d-1*t29*t41*t5
     5         *wght-2.481401963597599d+0*t18*t19*t43*wght-4.13566993932
     6         93317d-1*t18*t29*t32*wght-2.481401963597599d+0*t19*t27*t2
     7         9*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363490997d-1*t18*t5*(1.53
     1         15879268830906d-10*t20*t52/rhoa**1.1666666666666665d+1-2.
     2         1572522198836383d-7*gammaaa*t34*t42+7.292411320578055d-5*
     3         t21*t22)*wght-9.305257363490997d-1*t27*t45*t5*wght-1.2407
     4         009817987995d+0*t18*t19*t45*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490997d-1*t18*t5*(4.8538174947381857d-8*t34*t44-5.7434
     2         5472581159d-11*gammaaa*t52/rhoa**1.0666666666666666d+1)*w
     3         ght
            Mmat3(iq,D3_RA_RA_TA) = -9.305257363490997d-1*t4*t5*(-4.8040
     1         591281458866d+3*t12*t50*t56*t9+2.1020844511492584d+2*t14*
     2         t47*t5*t9+1.0748735402283824d+0*t16*t30*t40*t9-2.59929496
     3         6424086d+2*t10*t50*t54*t56-6.56644448895677d+3*t13*t39*t5
     4         0*t56+1.125673668754586d+4*t15*t26*t50*t56+3.736964261992
     5         054d+2*t14*t50*t56+2.2811839389397337d+1*t10*t39*t47*t5+7
     6         .097878108368861d+2*t13*t26*t47*t5+6.715586812412864d+0*t
     7         16*t47*t5-9.495236821536224d+2*t12*t15*t47*t5-1.668337065
     8         273382d-1*t10*t26*t30*t40+8.293023216396499d-1*t17*t30*t4
     9         0+4.9182565137832956d+0*t14*t15*t30*t40-6.65559866912399d
     :         +0*t12*t13*t30*t40)*wght-1.8610514726981994d+0*t24*t48*t5
     ;         *wght-9.305257363490997d-1*t31*t38*t5*wght-2.481401963597
     <         599d+0*t19*t4*t48*wght-4.1356699393293317d-1*t31*t32*t4*w
     =         ght-2.481401963597599d+0*t19*t24*t31*wght+Mmat3(iq,D3_RA_
     >         RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_GAA_TA) = -9.305257363490997d-1*t29*t48*t5*wg
     1         ht-9.305257363490997d-1*t31*t43*t5*wght-1.240700981798799
     2         5d+0*t19*t29*t31*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)-9.30525736
     1         3490997d-1*t31*t45*t5*wght
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_RA_TA_TA) = -9.305257363490997d-1*t4*t5*(2.88243
     1         54768875315d+3*t12*t57*t58*t9-2.1020844511492587d+2*t14*t
     2         46*t50*t9-3.2246206206851474d+0*t16*t25*t47*t9+1.55957697
     3         98544515d+2*t10*t54*t57*t58+3.9398666933740617d+3*t13*t39
     4         *t57*t58-6.754042012527514d+3*t15*t26*t57*t58-2.242178557
     5         1952324d+2*t14*t57*t58-2.2811839389397343d+1*t10*t39*t46*
     6         t50-7.097878108368861d+2*t13*t26*t46*t50-6.71558681241286
     7         9d+0*t16*t46*t50+9.495236821536223d+2*t12*t15*t46*t50+5.0
     8         05011195820146d-1*t10*t25*t26*t47-2.4879069649189497d+0*t
     9         17*t25*t47-1.4754769541349888d+1*t14*t15*t25*t47+1.996679
     :         6007371967d+1*t12*t13*t25*t47)*wght-9.305257363490997d-1*
     ;         t24*t5*t51*wght-1.2407009817987995d+0*t19*t4*t51*wght+Mma
     <         t3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)-9.3052573634
     1         90997d-1*t29*t5*t51*wght
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)-9.305257363490
     1         997d-1*t4*t5*(-1.7294612861325187d+3*t12*t59*t60*t9+5.804
     2         317117233265d+0*t16*t50*t6*t9+1.891876006034333d+2*t14*t4
     3         9*t58*t9-9.35746187912671d+1*t10*t54*t59*t60-2.3639200160
     4         244372d+3*t13*t39*t59*t60+4.0524252075165085d+3*t15*t26*t
     5         59*t60+1.3453071343171397d+2*t14*t59*t60-9.00902015247626
     6         4d-1*t10*t26*t50*t6+4.47823253685411d+0*t17*t50*t6+2.6558
     7         585174429794d+1*t14*t15*t50*t6-3.5940232813269546d+1*t12*
     8         t13*t50*t6+2.0530655450457608d+1*t10*t39*t49*t58+6.388090
     9         297531976d+2*t13*t26*t49*t58+6.044028131171577d+0*t16*t49
     :         *t58-8.5457131393826d+2*t12*t15*t49*t58)*wght
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 1.1832860901832055d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 5.777669212678251d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = rhob**1.6666666666666669d+0
            t7 = 1/taub
            t8 = 4.557799872345596d+0*t6*t7
            t9 = t8-1.0d+0
            t10 = t9**4
            t11 = t8+1.0d+0
            t12 = 1/t11**4
            t13 = t9**3
            t14 = 1/t11**3
            t15 = t9**2
            t16 = 1/t11**2
            t17 = 1/t11
            t18 = 1.0d+0*(-1.637571d-1*t17*t9-1.880028d-1*t15*t16-4.4906
     1         09d-1*t13*t14-8.2359d-3*t10*t12)+1.0d+0
            t19 = rhob**3.333333333333333d-1
            t20 = gammabb**2
            t21 = 1/t2**2
            t22 = 1/rhob**6.333333333333333d+0
            t23 = 1/rhob**3.6666666666666664d+0
            t24 = 1.8231028301445137d-5*t20*t21*t22-1.5407117900475334d-
     1         2*gammabb*t23*t3
            t25 = rhob**6.666666666666666d-1
            t26 = 1/t11**5
            t27 = -1.6123103103425737d+0*t16*t25*t7*t9+2.502505597910073
     1         d-1*t10*t25*t26*t7-1.2439534824594749d+0*t17*t25*t7-7.377
     2         384770674944d+0*t14*t15*t25*t7+9.983398003685984d+0*t12*t
     3         13*t25*t7
            t28 = 1/rhob**5.333333333333333d+0
            t29 = 5.777669212678251d-3*t1*t3-6.836635613041927d-6*gammab
     1         b*t21*t28
            t30 = 1/taub**2
            t31 = 9.673861862055442d-1*t16*t30*t6*t9-1.5015033587460438d
     1         -1*t10*t26*t30*t6+7.463720894756849d-1*t17*t30*t6+4.42643
     2         0862404966d+0*t14*t15*t30*t6-5.99003880221159d+0*t12*t13*
     3         t30*t6
            t32 = 1/t25
            t33 = gammabb**3
            t34 = 1/t2**3
            t35 = 1/rhob**10
            t36 = 1/rhob**7.333333333333333d+0
            t37 = 1/rhob**4.666666666666667d+0
            t38 = 5.649276563507623d-2*gammabb*t3*t37-1.6407925471300622
     1         d-4*t20*t21*t36+1.1505345172712736d-7*t33*t34*t35
            t39 = 1/t11**6
            t40 = 1/t19
            t41 = -1.0748735402283824d+0*t16*t40*t7*t9-8.758685213121911
     1         d+1*t14*t30*t5*t9+1.668337065273382d-1*t10*t26*t40*t7-8.2
     2         93023216396499d-1*t17*t40*t7-4.9182565137832956d+0*t14*t1
     3         5*t40*t7+6.65559866912399d+0*t12*t13*t40*t7-9.50493307891
     4         5558d+0*t10*t30*t39*t5-2.9574492118203594d+2*t13*t26*t30*
     5         t5-2.798161171838693d+0*t16*t30*t5+3.956348675640094d+2*t
     6         12*t15*t30*t5
            t42 = 1/rhob**9
            t43 = -4.3145044397672766d-8*t20*t34*t42-1.5407117900475334d
     1         -2*t23*t3+5.469308490433541d-5*gammabb*t21*t22
            t44 = 1/rhob**8
            t45 = 1.6179391649127287d-8*gammabb*t34*t44-1.36732712260838
     1         52d-5*t21*t28
            t46 = rhob**2.3333333333333334d+0
            t47 = 1/taub**3
            t48 = 5.255211127873146d+1*t14*t46*t47*t9+1.6123103103425737
     1         d+0*t16*t25*t30*t9+5.702959847349336d+0*t10*t39*t46*t47+1
     2         .7744695270922156d+2*t13*t26*t46*t47+1.6788967031032173d+
     3         0*t16*t46*t47-2.3738092053840557d+2*t12*t15*t46*t47-2.502
     4         505597910073d-1*t10*t25*t26*t30+1.2439534824594749d+0*t17
     5         *t25*t30+7.377384770674944d+0*t14*t15*t25*t30-9.983398003
     6         685984d+0*t12*t13*t25*t30
            t49 = rhob**3.3333333333333337d+0
            t50 = 1/taub**4
            t51 = -1.9347723724110882d+0*t16*t47*t6*t9-3.153126676723888
     1         d+1*t14*t49*t50*t9+3.0030067174920877d-1*t10*t26*t47*t6-1
     2         .4927441789513698d+0*t17*t47*t6-8.852861724809932d+0*t14*
     3         t15*t47*t6+1.1980077604423181d+1*t12*t13*t47*t6-3.4217759
     4         084096017d+0*t10*t39*t49*t50-1.0646817162553295d+2*t13*t2
     5         6*t49*t50-1.0073380218619295d+0*t16*t49*t50+1.42428552323
     6         04333d+2*t12*t15*t49*t50
            t52 = 1/t2**4
            t53 = rhob**2
            t54 = 1/t11**7
            t55 = 1/t5
            t56 = rhob**3
            t57 = rhob**4
            t58 = 1/taub**5
            t59 = rhob**5
            t60 = 1/taub**6
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t18*t4*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t27*t4*t5*wght-9.3052
     1         57363490997d-1*t18*t24*t5*wght-1.2407009817987995d+0*t18*
     2         t19*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t18*t
     1         29*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t31*t4*
     1         t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = -9.305257363490997d-1*t4*t41*t5*wght-9.
     1         305257363490997d-1*t18*t38*t5*wght-1.8610514726981994d+0*
     2         t24*t27*t5*wght-4.1356699393293317d-1*t18*t32*t4*wght-2.4
     3         81401963597599d+0*t19*t27*t4*wght-2.481401963597599d+0*t1
     4         8*t19*t24*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490997d-1*t18*t43*t5*wght-
     1         9.305257363490997d-1*t27*t29*t5*wght-1.2407009817987995d+
     2         0*t18*t19*t29*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         7d-1*t18*t45*t5*wght
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = -9.305257363490997d-1*t4*t48*t5*wght-9.
     1         305257363490997d-1*t24*t31*t5*wght-1.2407009817987995d+0*
     2         t19*t31*t4*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = Mmat2(iq,D2_GBB_TB)-9.305257363490997d
     1         -1*t29*t31*t5*wght
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)-9.305257363490997d-1
     1         *t4*t5*t51*wght
            Amat3(iq,D3_RB_RB_RB) = t32*(-1.2407009817987995d+0*t27*t4*w
     1         ght-1.2407009817987995d+0*t18*t24*wght)-9.305257363490997
     2         d-1*t4*t5*(3.582911800761275d-1*t16*t55*t7*t9+8.006765213
     3         576478d+3*t12*t47*t53*t9-1.751737042624382d+2*t14*t19*t30
     4         *t9-5.561123550911273d-2*t10*t26*t55*t7+2.764341072132166
     5         d-1*t17*t55*t7+1.639418837927765d+0*t14*t15*t55*t7-2.2185
     6         328897079964d+0*t12*t13*t55*t7+4.332158277373476d+2*t10*t
     7         47*t53*t54+1.094407414826128d+4*t13*t39*t47*t53-1.8761227
     8         81257643d+4*t15*t26*t47*t53-6.228273769986757d+2*t14*t47*
     9         t53-1.9009866157831115d+1*t10*t19*t30*t39-5.9148984236407
     :         19d+2*t13*t19*t26*t30-5.596322343677388d+0*t16*t19*t30+7.
     ;         912697351280188d+2*t12*t15*t19*t30)*wght+2.75711329288622
     <         1d-1*t18*t4*wght/t6-9.305257363490997d-1*t18*t5*(1.089129
     =         1924501978d-9*gammabb**4*t52/rhob**1.3666666666666666d+1-
     >         2.1860155828154199d-6*t33*t34/rhob**11-2.636329062970224d
     ?         -1*gammabb*t3/rhob**5.666666666666667d+0+1.38150681128728
     @         7d-3*t20*t21/rhob**8.333333333333334d+0)*wght-2.791577209
     1         047299d+0*t24*t41*t5*wght-2.791577209047299d+0*t27*t38*t5
     2         *wght-3.7221029453963983d+0*t19*t4*t41*wght-3.72210294539
     3         63983d+0*t18*t19*t38*wght-7.444205890792797d+0*t19*t24*t2
     4         7*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363490997d-1*t18*t5*(-4.08
     1         4234471688241d-10*t33*t52/rhob**1.2666666666666665d+1+5.6
     2         49276563507623d-2*t3*t37-3.9500561319797795d-4*gammabb*t2
     3         1*t36+7.334657547604368d-7*t20*t34*t35)*wght-1.8610514726
     4         981994d+0*t27*t43*t5*wght-9.305257363490997d-1*t29*t41*t5
     5         *wght-2.481401963597599d+0*t18*t19*t43*wght-4.13566993932
     6         93317d-1*t18*t29*t32*wght-2.481401963597599d+0*t19*t27*t2
     7         9*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363490997d-1*t18*t5*(1.53
     1         15879268830906d-10*t20*t52/rhob**1.1666666666666665d+1-2.
     2         1572522198836383d-7*gammabb*t34*t42+7.292411320578055d-5*
     3         t21*t22)*wght-9.305257363490997d-1*t27*t45*t5*wght-1.2407
     4         009817987995d+0*t18*t19*t45*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490997d-1*t18*t5*(4.8538174947381857d-8*t34*t44-5.7434
     2         5472581159d-11*gammabb*t52/rhob**1.0666666666666666d+1)*w
     3         ght
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = -9.305257363490997d-1*t4*t5*(-4.8040
     1         591281458866d+3*t12*t50*t56*t9+2.1020844511492584d+2*t14*
     2         t47*t5*t9+1.0748735402283824d+0*t16*t30*t40*t9-2.59929496
     3         6424086d+2*t10*t50*t54*t56-6.56644448895677d+3*t13*t39*t5
     4         0*t56+1.125673668754586d+4*t15*t26*t50*t56+3.736964261992
     5         054d+2*t14*t50*t56+2.2811839389397337d+1*t10*t39*t47*t5+7
     6         .097878108368861d+2*t13*t26*t47*t5+6.715586812412864d+0*t
     7         16*t47*t5-9.495236821536224d+2*t12*t15*t47*t5-1.668337065
     8         273382d-1*t10*t26*t30*t40+8.293023216396499d-1*t17*t30*t4
     9         0+4.9182565137832956d+0*t14*t15*t30*t40-6.65559866912399d
     :         +0*t12*t13*t30*t40)*wght-1.8610514726981994d+0*t24*t48*t5
     ;         *wght-9.305257363490997d-1*t31*t38*t5*wght-2.481401963597
     <         599d+0*t19*t4*t48*wght-4.1356699393293317d-1*t31*t32*t4*w
     =         ght-2.481401963597599d+0*t19*t24*t31*wght+Mmat3(iq,D3_RB_
     >         RB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = -9.305257363490997d-1*t29*t48*t5*wg
     1         ht-9.305257363490997d-1*t31*t43*t5*wght-1.240700981798799
     2         5d+0*t19*t29*t31*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBA_TB) = Mmat3(iq,D3_GBB_GBA_TB)-9.30525736
     1         3490997d-1*t31*t45*t5*wght
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = -9.305257363490997d-1*t4*t5*(2.88243
     1         54768875315d+3*t12*t57*t58*t9-2.1020844511492587d+2*t14*t
     2         46*t50*t9-3.2246206206851474d+0*t16*t25*t47*t9+1.55957697
     3         98544515d+2*t10*t54*t57*t58+3.9398666933740617d+3*t13*t39
     4         *t57*t58-6.754042012527514d+3*t15*t26*t57*t58-2.242178557
     5         1952324d+2*t14*t57*t58-2.2811839389397343d+1*t10*t39*t46*
     6         t50-7.097878108368861d+2*t13*t26*t46*t50-6.71558681241286
     7         9d+0*t16*t46*t50+9.495236821536223d+2*t12*t15*t46*t50+5.0
     8         05011195820146d-1*t10*t25*t26*t47-2.4879069649189497d+0*t
     9         17*t25*t47-1.4754769541349888d+1*t14*t15*t25*t47+1.996679
     :         6007371967d+1*t12*t13*t25*t47)*wght-9.305257363490997d-1*
     ;         t24*t5*t51*wght-1.2407009817987995d+0*t19*t4*t51*wght+Mma
     <         t3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
            Mmat3(iq,D3_GBB_TB_TB) = Mmat3(iq,D3_GBB_TB_TB)-9.3052573634
     1         90997d-1*t29*t5*t51*wght
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)-9.305257363490
     1         997d-1*t4*t5*(-1.7294612861325187d+3*t12*t59*t60*t9+5.804
     2         317117233265d+0*t16*t50*t6*t9+1.891876006034333d+2*t14*t4
     3         9*t58*t9-9.35746187912671d+1*t10*t54*t59*t60-2.3639200160
     4         244372d+3*t13*t39*t59*t60+4.0524252075165085d+3*t15*t26*t
     5         59*t60+1.3453071343171397d+2*t14*t59*t60-9.00902015247626
     6         4d-1*t10*t26*t50*t6+4.47823253685411d+0*t17*t50*t6+2.6558
     7         585174429794d+1*t14*t15*t50*t6-3.5940232813269546d+1*t12*
     8         t13*t50*t6+2.0530655450457608d+1*t10*t39*t49*t58+6.388090
     9         297531976d+2*t13*t26*t49*t58+6.044028131171577d+0*t16*t49
     :         *t58-8.5457131393826d+2*t12*t15*t49*t58)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
