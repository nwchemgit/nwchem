C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_m05.F 
C> The nwxcm_x_m05 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-1.0\\\\ 
C>   {\it t_5} &=& {\it t_3}+1.0\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-1.0\\\\ 
C>   {\it t_{19}} &=& {\it t_{17}}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-1.0\\\\ 
C>   {\it t_{23}} &=& {\it t_{21}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{16}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{16}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}\,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{\it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{\it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{\it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{\it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{\it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{\it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{\it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{\it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{\it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{\it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}\,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{\it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{\it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{\it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}\over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}+{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{\it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{\it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{20}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{20}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{\it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{\it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{\it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{\it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{\it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{\it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{\it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{\it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{\it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{\it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123, 161103 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
      subroutine nwxcm_x_m05(param,tol_rho,ipol,nq,wght, 
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
            t2 = 4.492769220958891d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509496d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 1/taua
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(10)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(9)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(8)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(7)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(6)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(5)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(4)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(3)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(2)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(1)
            t42 = 1/t12
            t43 = 1.0d+0*(t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37
     1         +t32*t33*t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*
     2         t22+t17*t18*t19+t14*t15*t16)+1.0d+0
            t44 = 1/t2**2
            t45 = rhoa**6.666666666666666d-1
            t46 = 1/t12**12
            t47 = -1.8610514726981994d+0*t43*(3.6121864536509496d-3*t1*t
     1         3-1.6228720119327633d-5*gammaaa*t44/rhoa**5.3333333333333
     2         33d+0)*t5*wght
            t48 = 1/taua**2
            fnc(iq) = fnc(iq)-1.8610514726981994d+0*t4*t43*t5*wght
            Amat(iq,D1_RA) = -1.8610514726981994d+0*t4*t5*(-8.3559664326
     1         33593d+1*t11*t45*t46*t6*t8+8.355966432633593d+1*t13*t15*t
     2         45*t6*t8+7.596333120575993d+0*t41*t42*t45*t8-7.5963331205
     3         75993d+0*t10*t40*t41*t45*t8+1.5192666241151986d+1*t10*t38
     4         *t40*t45*t8-1.5192666241151986d+1*t37*t38*t39*t45*t8+2.27
     5         88999361727982d+1*t35*t37*t39*t45*t8-2.2788999361727982d+
     6         1*t34*t35*t36*t45*t8+3.038533248230397d+1*t32*t34*t36*t45
     7         *t8-3.038533248230397d+1*t31*t32*t33*t45*t8+3.79816656028
     8         79966d+1*t29*t31*t33*t45*t8-3.7981665602879966d+1*t28*t29
     9         *t30*t45*t8+4.5577998723455965d+1*t26*t28*t30*t45*t8-4.55
     :         77998723455965d+1*t25*t26*t27*t45*t8+5.317433184403195d+1
     ;         *t23*t25*t27*t45*t8-5.317433184403195d+1*t22*t23*t24*t45*
     <         t8+6.077066496460795d+1*t20*t22*t24*t45*t8-6.077066496460
     =         795d+1*t19*t20*t21*t45*t8+6.836699808518394d+1*t17*t19*t2
     >         1*t45*t8-6.836699808518394d+1*t16*t17*t18*t45*t8+7.596333
     ?         120575993d+1*t14*t16*t18*t45*t8-7.596333120575993d+1*t13*
     @         t14*t15*t45*t8)*wght-1.8610514726981994d+0*t43*(4.3276586
     1         98487369d-5*gammaaa**2*t44/rhoa**6.333333333333333d+0-9.6
     2         32497209735864d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)
     3         *t5*wght-2.481401963597599d+0*rhoa**3.333333333333333d-1*
     4         t4*t43*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t47+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t47+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = Amat(iq,D1_TA)-1.8610514726981994d+0*t4*t5*
     1         (t13*t48*(4.5577998723455965d+1*t14*t15*t7-5.013579859580
     2         155d+1*t15*t6*t7)+t40*t48*(4.557799872345596d+0*t10*t41*t
     3         7-9.115599744691192d+0*t10*t38*t7)+t37*t48*(9.11559974469
     4         1192d+0*t38*t39*t7-1.3673399617036788d+1*t35*t39*t7)+t34*
     5         t48*(1.3673399617036788d+1*t35*t36*t7-1.8231199489382383d
     6         +1*t32*t36*t7)+t31*t48*(1.8231199489382383d+1*t32*t33*t7-
     7         2.2788999361727982d+1*t29*t33*t7)+t28*t48*(2.278899936172
     8         7982d+1*t29*t30*t7-2.7346799234073577d+1*t26*t30*t7)+t25*
     9         t48*(2.7346799234073577d+1*t26*t27*t7-3.190459910641917d+
     :         1*t23*t27*t7)+t22*t48*(3.190459910641917d+1*t23*t24*t7-3.
     ;         6462398978764765d+1*t20*t24*t7)+t19*t48*(3.64623989787647
     <         65d+1*t20*t21*t7-4.102019885111037d+1*t17*t21*t7)+t16*t48
     =         *(4.102019885111037d+1*t17*t18*t7-4.5577998723455965d+1*t
     >         14*t18*t7)+5.013579859580155d+1*t11*t46*t48*t6*t7-4.55779
     ?         9872345596d+0*t41*t42*t48*t7)*wght
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
            t2 = 4.492769220958891d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509496d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 1/taua
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(10)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(9)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(8)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(7)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(6)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(5)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(4)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(3)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(2)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(1)
            t42 = 1/t12
            t43 = 1.0d+0*(t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37
     1         +t32*t33*t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*
     2         t22+t17*t18*t19+t14*t15*t16)+1.0d+0
            t44 = 1/rhob**2.6666666666666666d+0
            t45 = 4.492769220958891d-3*gammabb*t44+1.0d+0
            t46 = 1/t45
            t47 = 3.6121864536509496d-3*gammabb*t44*t46+1.0d+0
            t48 = rhob**1.3333333333333333d+0
            t49 = rhob**1.6666666666666669d+0
            t50 = 1/taub
            t51 = 4.557799872345596d+0*t49*t50
            t52 = t51-1.0d+0
            t53 = t52**11
            t54 = t51+1.0d+0
            t55 = 1/t54**11
            t56 = t52**10
            t57 = 1/t54**10
            t58 = t52**9
            t59 = 1/t54**9
            t60 = t52**8
            t61 = 1/t54**8
            t62 = t52**7
            t63 = 1/t54**7
            t64 = t52**6
            t65 = 1/t54**6
            t66 = t52**5
            t67 = 1/t54**5
            t68 = t52**4
            t69 = 1/t54**4
            t70 = t52**3
            t71 = 1/t54**3
            t72 = t52**2
            t73 = 1/t54**2
            t74 = 1/t54
            t75 = 1.0d+0*(t41*t52*t74+t38*t72*t73+t35*t70*t71+t32*t68*t6
     1         9+t29*t66*t67+t26*t64*t65+t23*t62*t63+t20*t60*t61+t53*t55
     2         *t6+t17*t58*t59+t14*t56*t57)+1.0d+0
            t76 = 1/t2**2
            t77 = rhoa**6.666666666666666d-1
            t78 = 1/t12**12
            t79 = 1/t45**2
            t80 = rhob**6.666666666666666d-1
            t81 = 1/t54**12
            t82 = 1/taua**2
            t83 = 1/taub**2
            fnc(iq) = (-9.305257363490997d-1*t47*t48*t75-9.3052573634909
     1         97d-1*t4*t43*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-8.3559664326
     1         33593d+1*t11*t6*t77*t78*t8+8.355966432633593d+1*t13*t15*t
     2         6*t77*t8+7.596333120575993d+0*t41*t42*t77*t8-7.5963331205
     3         75993d+0*t10*t40*t41*t77*t8+1.5192666241151986d+1*t10*t38
     4         *t40*t77*t8-1.5192666241151986d+1*t37*t38*t39*t77*t8+2.27
     5         88999361727982d+1*t35*t37*t39*t77*t8-2.2788999361727982d+
     6         1*t34*t35*t36*t77*t8+3.038533248230397d+1*t32*t34*t36*t77
     7         *t8-3.038533248230397d+1*t31*t32*t33*t77*t8+3.79816656028
     8         79966d+1*t29*t31*t33*t77*t8-3.7981665602879966d+1*t28*t29
     9         *t30*t77*t8+4.5577998723455965d+1*t26*t28*t30*t77*t8-4.55
     :         77998723455965d+1*t25*t26*t27*t77*t8+5.317433184403195d+1
     ;         *t23*t25*t27*t77*t8-5.317433184403195d+1*t22*t23*t24*t77*
     <         t8+6.077066496460795d+1*t20*t22*t24*t77*t8-6.077066496460
     =         795d+1*t19*t20*t21*t77*t8+6.836699808518394d+1*t17*t19*t2
     >         1*t77*t8-6.836699808518394d+1*t16*t17*t18*t77*t8+7.596333
     ?         120575993d+1*t14*t16*t18*t77*t8-7.596333120575993d+1*t13*
     @         t14*t15*t77*t8)-9.305257363490997d-1*t43*t5*(4.3276586984
     1         87369d-5*gammaaa**2*t76/rhoa**6.333333333333333d+0-9.6324
     2         97209735864d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)-1.
     3         2407009817987995d+0*rhoa**3.333333333333333d-1*t4*t43)*wg
     4         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t47*t48*(-8.35596643
     1         2633593d+1*t50*t53*t6*t80*t81+7.596333120575993d+0*t41*t5
     2         0*t74*t80-7.596333120575993d+0*t41*t50*t52*t73*t80+1.5192
     3         666241151986d+1*t38*t50*t52*t73*t80-1.5192666241151986d+1
     4         *t38*t50*t71*t72*t80+2.2788999361727982d+1*t35*t50*t71*t7
     5         2*t80-2.2788999361727982d+1*t35*t50*t69*t70*t80+3.0385332
     6         48230397d+1*t32*t50*t69*t70*t80-3.038533248230397d+1*t32*
     7         t50*t67*t68*t80+3.7981665602879966d+1*t29*t50*t67*t68*t80
     8         -3.7981665602879966d+1*t29*t50*t65*t66*t80+4.557799872345
     9         5965d+1*t26*t50*t65*t66*t80-4.5577998723455965d+1*t26*t50
     :         *t63*t64*t80+5.317433184403195d+1*t23*t50*t63*t64*t80-5.3
     ;         17433184403195d+1*t23*t50*t61*t62*t80+6.077066496460795d+
     <         1*t20*t50*t61*t62*t80-6.077066496460795d+1*t20*t50*t59*t6
     =         0*t80+6.836699808518394d+1*t17*t50*t59*t60*t80+8.35596643
     >         2633593d+1*t50*t55*t56*t6*t80-6.836699808518394d+1*t17*t5
     ?         0*t57*t58*t80+7.596333120575993d+1*t14*t50*t57*t58*t80-7.
     @         596333120575993d+1*t14*t50*t55*t56*t80)-9.305257363490997
     1         d-1*t48*t75*(4.327658698487369d-5*gammabb**2*t79/rhob**6.
     2         333333333333333d+0-9.632497209735864d-3*gammabb*t46/rhob*
     3         *3.6666666666666664d+0)-1.2407009817987995d+0*rhob**3.333
     4         333333333333d-1*t47*t75)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t43*t
     1         5*(3.6121864536509496d-3*t1*t3-1.6228720119327633d-5*gamm
     2         aaa*t76/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t48*t
     1         75*(3.6121864536509496d-3*t44*t46-1.6228720119327633d-5*g
     2         ammabb*t79/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         5.013579859580155d+1*t11*t6*t7*t78*t82+t13*(4.55779987234
     2         55965d+1*t14*t15*t7-5.013579859580155d+1*t15*t6*t7)*t82+t
     3         40*(4.557799872345596d+0*t10*t41*t7-9.115599744691192d+0*
     4         t10*t38*t7)*t82+t37*(9.115599744691192d+0*t38*t39*t7-1.36
     5         73399617036788d+1*t35*t39*t7)*t82+t34*(1.3673399617036788
     6         d+1*t35*t36*t7-1.8231199489382383d+1*t32*t36*t7)*t82+t31*
     7         (1.8231199489382383d+1*t32*t33*t7-2.2788999361727982d+1*t
     8         29*t33*t7)*t82+t28*(2.2788999361727982d+1*t29*t30*t7-2.73
     9         46799234073577d+1*t26*t30*t7)*t82+t25*(2.7346799234073577
     :         d+1*t26*t27*t7-3.190459910641917d+1*t23*t27*t7)*t82+t22*(
     ;         3.190459910641917d+1*t23*t24*t7-3.6462398978764765d+1*t20
     <         *t24*t7)*t82+t19*(3.6462398978764765d+1*t20*t21*t7-4.1020
     =         19885111037d+1*t17*t21*t7)*t82+t16*(4.102019885111037d+1*
     >         t17*t18*t7-4.5577998723455965d+1*t14*t18*t7)*t82-4.557799
     ?         872345596d+0*t41*t42*t7*t82)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t47*t48
     1         *(5.013579859580155d+1*t49*t53*t6*t81*t83-4.5577998723455
     2         96d+0*t41*t49*t74*t83+(4.557799872345596d+0*t41*t49*t52-9
     3         .115599744691192d+0*t38*t49*t52)*t73*t83+t71*(9.115599744
     4         691192d+0*t38*t49*t72-1.3673399617036788d+1*t35*t49*t72)*
     5         t83+t69*(1.3673399617036788d+1*t35*t49*t70-1.823119948938
     6         2383d+1*t32*t49*t70)*t83+t67*(1.8231199489382383d+1*t32*t
     7         49*t68-2.2788999361727982d+1*t29*t49*t68)*t83+t65*(2.2788
     8         999361727982d+1*t29*t49*t66-2.7346799234073577d+1*t26*t49
     9         *t66)*t83+t63*(2.7346799234073577d+1*t26*t49*t64-3.190459
     :         910641917d+1*t23*t49*t64)*t83+t61*(3.190459910641917d+1*t
     ;         23*t49*t62-3.6462398978764765d+1*t20*t49*t62)*t83+t59*(3.
     <         6462398978764765d+1*t20*t49*t60-4.102019885111037d+1*t17*
     =         t49*t60)*t83+t55*(4.5577998723455965d+1*t14*t49*t56-5.013
     >         579859580155d+1*t49*t56*t6)*t83+t57*(4.102019885111037d+1
     ?         *t17*t49*t58-4.5577998723455965d+1*t14*t49*t58)*t83)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.492769220958891d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509496d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(11)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 1/taua
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(10)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(9)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(8)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(7)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(6)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(5)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(4)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(3)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(2)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(1)
            t42 = 1/t12
            t43 = 1.0d+0*(t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37
     1         +t32*t33*t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*
     2         t22+t17*t18*t19+t14*t15*t16)+1.0d+0
            t44 = 1/t2**2
            t45 = rhoa**6.666666666666666d-1
            t46 = 1/t12**12
            t47 = 1/taua**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t43*t5*wght
            Amat(iq,D1_RA) = -9.305257363490997d-1*t4*t5*(-8.35596643263
     1         3593d+1*t11*t45*t46*t6*t8+8.355966432633593d+1*t13*t15*t4
     2         5*t6*t8+7.596333120575993d+0*t41*t42*t45*t8-7.59633312057
     3         5993d+0*t10*t40*t41*t45*t8+1.5192666241151986d+1*t10*t38*
     4         t40*t45*t8-1.5192666241151986d+1*t37*t38*t39*t45*t8+2.278
     5         8999361727982d+1*t35*t37*t39*t45*t8-2.2788999361727982d+1
     6         *t34*t35*t36*t45*t8+3.038533248230397d+1*t32*t34*t36*t45*
     7         t8-3.038533248230397d+1*t31*t32*t33*t45*t8+3.798166560287
     8         9966d+1*t29*t31*t33*t45*t8-3.7981665602879966d+1*t28*t29*
     9         t30*t45*t8+4.5577998723455965d+1*t26*t28*t30*t45*t8-4.557
     :         7998723455965d+1*t25*t26*t27*t45*t8+5.317433184403195d+1*
     ;         t23*t25*t27*t45*t8-5.317433184403195d+1*t22*t23*t24*t45*t
     <         8+6.077066496460795d+1*t20*t22*t24*t45*t8-6.0770664964607
     =         95d+1*t19*t20*t21*t45*t8+6.836699808518394d+1*t17*t19*t21
     >         *t45*t8-6.836699808518394d+1*t16*t17*t18*t45*t8+7.5963331
     ?         20575993d+1*t14*t16*t18*t45*t8-7.596333120575993d+1*t13*t
     @         14*t15*t45*t8)*wght-9.305257363490997d-1*t43*(4.327658698
     1         487369d-5*gammaaa**2*t44/rhoa**6.333333333333333d+0-9.632
     2         497209735864d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)*t
     3         5*wght-1.2407009817987995d+0*rhoa**3.333333333333333d-1*t
     4         4*t43*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t43*(
     1         3.6121864536509496d-3*t1*t3-1.6228720119327633d-5*gammaaa
     2         *t44/rhoa**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         t13*t47*(4.5577998723455965d+1*t14*t15*t7-5.0135798595801
     2         55d+1*t15*t6*t7)+t40*t47*(4.557799872345596d+0*t10*t41*t7
     3         -9.115599744691192d+0*t10*t38*t7)+t37*t47*(9.115599744691
     4         192d+0*t38*t39*t7-1.3673399617036788d+1*t35*t39*t7)+t34*t
     5         47*(1.3673399617036788d+1*t35*t36*t7-1.8231199489382383d+
     6         1*t32*t36*t7)+t31*t47*(1.8231199489382383d+1*t32*t33*t7-2
     7         .2788999361727982d+1*t29*t33*t7)+t28*t47*(2.2788999361727
     8         982d+1*t29*t30*t7-2.7346799234073577d+1*t26*t30*t7)+t25*t
     9         47*(2.7346799234073577d+1*t26*t27*t7-3.190459910641917d+1
     :         *t23*t27*t7)+t22*t47*(3.190459910641917d+1*t23*t24*t7-3.6
     ;         462398978764765d+1*t20*t24*t7)+t19*t47*(3.646239897876476
     <         5d+1*t20*t21*t7-4.102019885111037d+1*t17*t21*t7)+t16*t47*
     =         (4.102019885111037d+1*t17*t18*t7-4.5577998723455965d+1*t1
     >         4*t18*t7)+5.013579859580155d+1*t11*t46*t47*t6*t7-4.557799
     ?         872345596d+0*t41*t42*t47*t7)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 4.492769220958891d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509496d-3*gammabb*t1*t3+1.0d+0
            t5 = rhob**1.3333333333333333d+0
            t6 = param(11)
            t7 = rhob**1.6666666666666669d+0
            t8 = 1/taub
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(10)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(9)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(8)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(7)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(6)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(5)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(4)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(3)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(2)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(1)
            t42 = 1/t12
            t43 = 1.0d+0*(t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37
     1         +t32*t33*t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*
     2         t22+t17*t18*t19+t14*t15*t16)+1.0d+0
            t44 = 1/t2**2
            t45 = rhob**6.666666666666666d-1
            t46 = 1/t12**12
            t47 = 1/taub**2
            fnc(iq) = fnc(iq)-9.305257363490997d-1*t4*t43*t5*wght
            Amat(iq,D1_RB) = -9.305257363490997d-1*t4*t5*(-8.35596643263
     1         3593d+1*t11*t45*t46*t6*t8+8.355966432633593d+1*t13*t15*t4
     2         5*t6*t8+7.596333120575993d+0*t41*t42*t45*t8-7.59633312057
     3         5993d+0*t10*t40*t41*t45*t8+1.5192666241151986d+1*t10*t38*
     4         t40*t45*t8-1.5192666241151986d+1*t37*t38*t39*t45*t8+2.278
     5         8999361727982d+1*t35*t37*t39*t45*t8-2.2788999361727982d+1
     6         *t34*t35*t36*t45*t8+3.038533248230397d+1*t32*t34*t36*t45*
     7         t8-3.038533248230397d+1*t31*t32*t33*t45*t8+3.798166560287
     8         9966d+1*t29*t31*t33*t45*t8-3.7981665602879966d+1*t28*t29*
     9         t30*t45*t8+4.5577998723455965d+1*t26*t28*t30*t45*t8-4.557
     :         7998723455965d+1*t25*t26*t27*t45*t8+5.317433184403195d+1*
     ;         t23*t25*t27*t45*t8-5.317433184403195d+1*t22*t23*t24*t45*t
     <         8+6.077066496460795d+1*t20*t22*t24*t45*t8-6.0770664964607
     =         95d+1*t19*t20*t21*t45*t8+6.836699808518394d+1*t17*t19*t21
     >         *t45*t8-6.836699808518394d+1*t16*t17*t18*t45*t8+7.5963331
     ?         20575993d+1*t14*t16*t18*t45*t8-7.596333120575993d+1*t13*t
     @         14*t15*t45*t8)*wght-9.305257363490997d-1*t43*(4.327658698
     1         487369d-5*gammabb**2*t44/rhob**6.333333333333333d+0-9.632
     2         497209735864d-3*gammabb*t3/rhob**3.6666666666666664d+0)*t
     3         5*wght-1.2407009817987995d+0*rhob**3.333333333333333d-1*t
     4         4*t43*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t43*(
     1         3.6121864536509496d-3*t1*t3-1.6228720119327633d-5*gammabb
     2         *t44/rhob**5.333333333333333d+0)*t5*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t4*t5*(
     1         t13*t47*(4.5577998723455965d+1*t14*t15*t7-5.0135798595801
     2         55d+1*t15*t6*t7)+t40*t47*(4.557799872345596d+0*t10*t41*t7
     3         -9.115599744691192d+0*t10*t38*t7)+t37*t47*(9.115599744691
     4         192d+0*t38*t39*t7-1.3673399617036788d+1*t35*t39*t7)+t34*t
     5         47*(1.3673399617036788d+1*t35*t36*t7-1.8231199489382383d+
     6         1*t32*t36*t7)+t31*t47*(1.8231199489382383d+1*t32*t33*t7-2
     7         .2788999361727982d+1*t29*t33*t7)+t28*t47*(2.2788999361727
     8         982d+1*t29*t30*t7-2.7346799234073577d+1*t26*t30*t7)+t25*t
     9         47*(2.7346799234073577d+1*t26*t27*t7-3.190459910641917d+1
     :         *t23*t27*t7)+t22*t47*(3.190459910641917d+1*t23*t24*t7-3.6
     ;         462398978764765d+1*t20*t24*t7)+t19*t47*(3.646239897876476
     <         5d+1*t20*t21*t7-4.102019885111037d+1*t17*t21*t7)+t16*t47*
     =         (4.102019885111037d+1*t17*t18*t7-4.5577998723455965d+1*t1
     >         4*t18*t7)+5.013579859580155d+1*t11*t46*t47*t6*t7-4.557799
     ?         872345596d+0*t41*t42*t47*t7)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-1.0\\\\ 
C>   {\it t_5} &=& {\it t_3}+1.0\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-1.0\\\\ 
C>   {\it t_{19}} &=& {\it t_{17}}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-1.0\\\\ 
C>   {\it t_{23}} &=& {\it t_{21}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{16}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{16}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}\,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{\it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{\it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{\it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{\it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{\it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{\it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{\it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{\it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{\it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{\it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}\,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{\it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{\it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{\it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}\over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}+{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{\it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{\it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{20}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{20}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{\it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{\it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{\it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{\it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{\it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{\it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{\it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{\it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{\it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{\it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123, 161103 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
      subroutine nwxcm_x_m05_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_m05 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_3} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_4} &=& {\it t_3}-1.0\\\\ 
C>   {\it t_5} &=& {\it t_3}+1.0\\\\ 
C>   {\it t_6} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}-1.0\\\\ 
C>   {\it t_{19}} &=& {\it t_{17}}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{21}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}-1.0\\\\ 
C>   {\it t_{23}} &=& {\it t_{21}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{16}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{16}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{18}}}\over{{\it t_{19}}}}+{{{\it t_{14}}\,{\it t_{18}}^2}\over{{\it t_{19}}^2}}+{{{\it t_{13}}\,{\it t_{18}}^3}\over{{\it t_{19}}^3}}+{{{\it t_{12}}\,{\it t_{18}}^4}\over{{\it t_{19}}^4}}+{{{\it t_{11}}\,{\it t_{18}}^5}\over{{\it t_{19}}^5}}+{{{\it t_{10}}\,{\it t_{18}}^6}\over{{\it t_{19}}^6}}+{{{\it t_9}\,{\it t_{18}}^7}\over{{\it t_{19}}^7}}+{{{\it t_8}\,{\it t_{18}}^8}\over{{\it t_{19}}^8}}+{{{\it t_7}\,{\it t_{18}}^9}\over{{\it t_{19}}^9}}+{{{\it t_6}\,{\it t_{18}}^{10}}\over{{\it t_{19}}^{10}}}+{{{\it t_2}\,{\it t_{18}}^{11}}\over{{\it t_{19}}^{11}}}\right)+1.0\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_4}}\over{{\it t_5}}}+{{{\it t_{14}}\,{\it t_4}^2}\over{{\it t_5}^2}}+{{{\it t_{13}}\,{\it t_4}^3}\over{{\it t_5}^3}}+{{{\it t_{12}}\,{\it t_4}^4}\over{{\it t_5}^4}}+{{{\it t_{11}}\,{\it t_4}^5}\over{{\it t_5}^5}}+{{{\it t_{10}}\,{\it t_4}^6}\over{{\it t_5}^6}}+{{{\it t_9}\,{\it t_4}^7}\over{{\it t_5}^7}}+{{{\it t_8}\,{\it t_4}^8}\over{{\it t_5}^8}}+{{{\it t_7}\,{\it t_4}^9}\over{{\it t_5}^9}}+{{{\it t_6}\,{\it t_4}^{10}}\over{{\it t_5}^{10}}}+{{{\it t_2}\,{\it t_4}^{11}}\over{{\it t_5}^{11}}}\right)+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{20}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{20}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left(1.0\,\left({{{\it t_{15}}\,{\it t_{22}}}\over{{\it t_{23}}}}+{{{\it t_{14}}\,{\it t_{22}}^2}\over{{\it t_{23}}^2}}+{{{\it t_{13}}\,{\it t_{22}}^3}\over{{\it t_{23}}^3}}+{{{\it t_{12}}\,{\it t_{22}}^4}\over{{\it t_{23}}^4}}+{{{\it t_{11}}\,{\it t_{22}}^5}\over{{\it t_{23}}^5}}+{{{\it t_{10}}\,{\it t_{22}}^6}\over{{\it t_{23}}^6}}+{{{\it t_9}\,{\it t_{22}}^7}\over{{\it t_{23}}^7}}+{{{\it t_8}\,{\it t_{22}}^8}\over{{\it t_{23}}^8}}+{{{\it t_7}\,{\it t_{22}}^9}\over{{\it t_{23}}^9}}+{{{\it t_6}\,{\it t_{22}}^{10}}\over{{\it t_{23}}^{10}}}+{{{\it t_2}\,{\it t_{22}}^{11}}\over{{\it t_{23}}^{11}}}\right)+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, NE Schultz, DG Truhlar, J.Chem.Phys. 123, 161103 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2126975 "> 
C> 10.1063/1.2126975 </a> 
C> 
      subroutine nwxcm_x_m05_d3(param,tol_rho,ipol,nq,wght, 
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
