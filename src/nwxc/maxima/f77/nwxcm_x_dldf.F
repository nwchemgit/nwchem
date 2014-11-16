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
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
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
            t22 = -1.8610514726981994d+0*t18*(5.777669212678251d-3*t1*t3
     1         -6.836635613041927d-6*gammaaa*t19/rhoa**5.333333333333333
     2         d+0)*t5*wght
            t23 = 1/taua**2
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t18*t4*t5*wght
            Amat(iq,D1_RA) = -1.8610514726981994d+0*t4*t5*(-1.6123103103
     1         425737d+0*t16*t20*t7*t9+2.502505597910073d-1*t10*t20*t21*
     2         t7-1.2439534824594749d+0*t17*t20*t7-7.377384770674944d+0*
     3         t14*t15*t20*t7+9.983398003685984d+0*t12*t13*t20*t7)*wght-
     4         1.8610514726981994d+0*t18*(1.8231028301445137d-5*gammaaa*
     5         *2*t19/rhoa**6.333333333333333d+0-1.5407117900475334d-2*g
     6         ammaaa*t3/rhoa**3.6666666666666664d+0)*t5*wght-2.48140196
     7         3597599d+0*rhoa**3.333333333333333d-1*t18*t4*wght+Amat(iq
     8         ,D1_RA)
            Cmat(iq,D1_GAA) = t22+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t22+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = Amat(iq,D1_TA)-1.8610514726981994d+0*t4*t5*
     1         (9.673861862055442d-1*t16*t23*t6*t9-1.5015033587460438d-1
     2         *t10*t21*t23*t6+7.463720894756849d-1*t17*t23*t6+4.4264308
     3         62404966d+0*t14*t15*t23*t6-5.99003880221159d+0*t12*t13*t2
     4         3*t6)*wght
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
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
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
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_5}\,\sigma_{\beta\beta}}\over{.001183286090183205\,{\it t_5}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_7}}\over{{\it t_8}}}-{{0.1880028\,{\it t_7}^2}\over{{\it t_8}^2}}-{{0.4490609\,{\it t_7}^3}\over{{\it t_8}^3}}-{{0.0082359\,{\it t_7}^4}\over{{\it t_8}^4}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.001183286090183205\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_3}}\over{{\it t_4}}}-{{0.1880028\,{\it t_3}^2}\over{{\it t_4}^2}}-{{0.4490609\,{\it t_3}^3}\over{{\it t_4}^3}}-{{0.0082359\,{\it t_3}^4}\over{{\it t_4}^4}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{.005777669212678251\,{\it t_9}\,\sigma_{ss}}\over{.001183286090183205\,{\it t_9}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left(-{{0.1637571\,{\it t_{11}}}\over{{\it t_{12}}}}-{{0.1880028\,{\it t_{11}}^2}\over{{\it t_{12}}^2}}-{{0.4490609\,{\it t_{11}}^3}\over{{\it t_{12}}^3}}-{{0.0082359\,{\it t_{11}}^4}\over{{\it t_{12}}^4}}\right)+1.0\right)\\\\ 
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
