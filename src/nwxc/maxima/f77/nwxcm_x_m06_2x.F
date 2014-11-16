C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_m06_2x.F 
C> The nwxcm_x_m06_2x functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-1.0\\\\ 
C>   {\it t_6} &=& {\it t_4}+1.0\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {\it t_{18}}+1.0\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-1.0\\\\ 
C>   {\it t_{24}} &=& {\it t_{22}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{17}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{17}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{\it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{\it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{\it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{\it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{\it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{\it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{\it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{\it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{\it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{\it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{\it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{\it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{\it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{\it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{\it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{\it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}\over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}+{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{\it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{\it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{21}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{21}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{\it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{\it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{\it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{\it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{\it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{\it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{\it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{\it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{\it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{\it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{\it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
      subroutine nwxcm_x_m06_2x(param,tol_rho,ipol,nq,wght, 
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
            t6 = param(12)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 1/taua
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(11)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(10)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(9)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(8)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(7)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(6)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(5)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(4)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(3)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(2)
            t42 = 1/t12
            t43 = t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37+t32*t33
     1         *t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*t22+t17*
     2         t18*t19+t14*t15*t16+param(1)
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
            t6 = param(1)
            t7 = param(12)
            t8 = rhoa**1.6666666666666669d+0
            t9 = 1/taua
            t10 = 4.557799872345596d+0*t8*t9
            t11 = t10-1.0d+0
            t12 = t11**11
            t13 = t10+1.0d+0
            t14 = 1/t13**11
            t15 = param(11)
            t16 = t11**10
            t17 = 1/t13**10
            t18 = param(10)
            t19 = t11**9
            t20 = 1/t13**9
            t21 = param(9)
            t22 = t11**8
            t23 = 1/t13**8
            t24 = param(8)
            t25 = t11**7
            t26 = 1/t13**7
            t27 = param(7)
            t28 = t11**6
            t29 = 1/t13**6
            t30 = param(6)
            t31 = t11**5
            t32 = 1/t13**5
            t33 = param(5)
            t34 = t11**4
            t35 = 1/t13**4
            t36 = param(4)
            t37 = t11**3
            t38 = 1/t13**3
            t39 = param(3)
            t40 = t11**2
            t41 = 1/t13**2
            t42 = param(2)
            t43 = 1/t13
            t44 = t12*t14*t7+t6+t11*t42*t43+t39*t40*t41+t36*t37*t38+t33*
     1         t34*t35+t30*t31*t32+t27*t28*t29+t24*t25*t26+t21*t22*t23+t
     2         18*t19*t20+t15*t16*t17
            t45 = 1/rhob**2.6666666666666666d+0
            t46 = 4.492769220958891d-3*gammabb*t45+1.0d+0
            t47 = 1/t46
            t48 = 3.6121864536509496d-3*gammabb*t45*t47+1.0d+0
            t49 = rhob**1.3333333333333333d+0
            t50 = rhob**1.6666666666666669d+0
            t51 = 1/taub
            t52 = 4.557799872345596d+0*t50*t51
            t53 = t52-1.0d+0
            t54 = t53**11
            t55 = t52+1.0d+0
            t56 = 1/t55**11
            t57 = t53**10
            t58 = 1/t55**10
            t59 = t53**9
            t60 = 1/t55**9
            t61 = t53**8
            t62 = 1/t55**8
            t63 = t53**7
            t64 = 1/t55**7
            t65 = t53**6
            t66 = 1/t55**6
            t67 = t53**5
            t68 = 1/t55**5
            t69 = t53**4
            t70 = 1/t55**4
            t71 = t53**3
            t72 = 1/t55**3
            t73 = t53**2
            t74 = 1/t55**2
            t75 = 1/t55
            t76 = t42*t53*t75+t39*t73*t74+t36*t71*t72+t33*t69*t70+t54*t5
     1         6*t7+t30*t67*t68+t27*t65*t66+t24*t63*t64+t21*t61*t62+t18*
     2         t59*t60+t6+t15*t57*t58
            t77 = 1/t2**2
            t78 = rhoa**6.666666666666666d-1
            t79 = 1/t13**12
            t80 = 1/t46**2
            t81 = rhob**6.666666666666666d-1
            t82 = 1/t55**12
            t83 = 1/taua**2
            t84 = 1/taub**2
            fnc(iq) = (-9.305257363490997d-1*t48*t49*t76-9.3052573634909
     1         97d-1*t4*t44*t5)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363490997d-1*t4*t5*(-8.3559664326
     1         33593d+1*t12*t7*t78*t79*t9+8.355966432633593d+1*t14*t16*t
     2         7*t78*t9+7.596333120575993d+0*t42*t43*t78*t9-7.5963331205
     3         75993d+0*t11*t41*t42*t78*t9+1.5192666241151986d+1*t11*t39
     4         *t41*t78*t9-1.5192666241151986d+1*t38*t39*t40*t78*t9+2.27
     5         88999361727982d+1*t36*t38*t40*t78*t9-2.2788999361727982d+
     6         1*t35*t36*t37*t78*t9+3.038533248230397d+1*t33*t35*t37*t78
     7         *t9-3.038533248230397d+1*t32*t33*t34*t78*t9+3.79816656028
     8         79966d+1*t30*t32*t34*t78*t9-3.7981665602879966d+1*t29*t30
     9         *t31*t78*t9+4.5577998723455965d+1*t27*t29*t31*t78*t9-4.55
     :         77998723455965d+1*t26*t27*t28*t78*t9+5.317433184403195d+1
     ;         *t24*t26*t28*t78*t9-5.317433184403195d+1*t23*t24*t25*t78*
     <         t9+6.077066496460795d+1*t21*t23*t25*t78*t9-6.077066496460
     =         795d+1*t20*t21*t22*t78*t9+6.836699808518394d+1*t18*t20*t2
     >         2*t78*t9-6.836699808518394d+1*t17*t18*t19*t78*t9+7.596333
     ?         120575993d+1*t15*t17*t19*t78*t9-7.596333120575993d+1*t14*
     @         t15*t16*t78*t9)-9.305257363490997d-1*t44*t5*(4.3276586984
     1         87369d-5*gammaaa**2*t77/rhoa**6.333333333333333d+0-9.6324
     2         97209735864d-3*gammaaa*t3/rhoa**3.6666666666666664d+0)-1.
     3         2407009817987995d+0*rhoa**3.333333333333333d-1*t4*t44)*wg
     4         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t48*t49*(-8.35596643
     1         2633593d+1*t51*t54*t7*t81*t82+7.596333120575993d+0*t42*t5
     2         1*t75*t81-7.596333120575993d+0*t42*t51*t53*t74*t81+1.5192
     3         666241151986d+1*t39*t51*t53*t74*t81-1.5192666241151986d+1
     4         *t39*t51*t72*t73*t81+2.2788999361727982d+1*t36*t51*t72*t7
     5         3*t81-2.2788999361727982d+1*t36*t51*t70*t71*t81+3.0385332
     6         48230397d+1*t33*t51*t70*t71*t81+8.355966432633593d+1*t51*
     7         t56*t57*t7*t81-3.038533248230397d+1*t33*t51*t68*t69*t81+3
     8         .7981665602879966d+1*t30*t51*t68*t69*t81-3.79816656028799
     9         66d+1*t30*t51*t66*t67*t81+4.5577998723455965d+1*t27*t51*t
     :         66*t67*t81-4.5577998723455965d+1*t27*t51*t64*t65*t81+5.31
     ;         7433184403195d+1*t24*t51*t64*t65*t81-5.317433184403195d+1
     <         *t24*t51*t62*t63*t81+6.077066496460795d+1*t21*t51*t62*t63
     =         *t81-6.077066496460795d+1*t21*t51*t60*t61*t81+6.836699808
     >         518394d+1*t18*t51*t60*t61*t81-6.836699808518394d+1*t18*t5
     ?         1*t58*t59*t81+7.596333120575993d+1*t15*t51*t58*t59*t81-7.
     @         596333120575993d+1*t15*t51*t56*t57*t81)-9.305257363490997
     1         d-1*t49*t76*(4.327658698487369d-5*gammabb**2*t80/rhob**6.
     2         333333333333333d+0-9.632497209735864d-3*gammabb*t47/rhob*
     3         *3.6666666666666664d+0)-1.2407009817987995d+0*rhob**3.333
     4         333333333333d-1*t48*t76)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t44*t
     1         5*(3.6121864536509496d-3*t1*t3-1.6228720119327633d-5*gamm
     2         aaa*t77/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t49*t
     1         76*(3.6121864536509496d-3*t45*t47-1.6228720119327633d-5*g
     2         ammabb*t80/rhob**5.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t4*t5*(
     1         t14*(4.5577998723455965d+1*t15*t16*t8-5.013579859580155d+
     2         1*t16*t7*t8)*t83+t41*(4.557799872345596d+0*t11*t42*t8-9.1
     3         15599744691192d+0*t11*t39*t8)*t83+t38*(9.115599744691192d
     4         +0*t39*t40*t8-1.3673399617036788d+1*t36*t40*t8)*t83+t35*(
     5         1.3673399617036788d+1*t36*t37*t8-1.8231199489382383d+1*t3
     6         3*t37*t8)*t83+t32*(1.8231199489382383d+1*t33*t34*t8-2.278
     7         8999361727982d+1*t30*t34*t8)*t83+t29*(2.2788999361727982d
     8         +1*t30*t31*t8-2.7346799234073577d+1*t27*t31*t8)*t83+t26*(
     9         2.7346799234073577d+1*t27*t28*t8-3.190459910641917d+1*t24
     :         *t28*t8)*t83+t23*(3.190459910641917d+1*t24*t25*t8-3.64623
     ;         98978764765d+1*t21*t25*t8)*t83+t20*(3.6462398978764765d+1
     <         *t21*t22*t8-4.102019885111037d+1*t18*t22*t8)*t83+t17*(4.1
     =         02019885111037d+1*t18*t19*t8-4.5577998723455965d+1*t15*t1
     >         9*t8)*t83+5.013579859580155d+1*t12*t7*t79*t8*t83-4.557799
     ?         872345596d+0*t42*t43*t8*t83)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t48*t49
     1         *(5.013579859580155d+1*t50*t54*t7*t82*t84-4.5577998723455
     2         96d+0*t42*t50*t75*t84+(4.557799872345596d+0*t42*t50*t53-9
     3         .115599744691192d+0*t39*t50*t53)*t74*t84+t72*(9.115599744
     4         691192d+0*t39*t50*t73-1.3673399617036788d+1*t36*t50*t73)*
     5         t84+t70*(1.3673399617036788d+1*t36*t50*t71-1.823119948938
     6         2383d+1*t33*t50*t71)*t84+t56*(4.5577998723455965d+1*t15*t
     7         50*t57-5.013579859580155d+1*t50*t57*t7)*t84+t68*(1.823119
     8         9489382383d+1*t33*t50*t69-2.2788999361727982d+1*t30*t50*t
     9         69)*t84+t66*(2.2788999361727982d+1*t30*t50*t67-2.73467992
     :         34073577d+1*t27*t50*t67)*t84+t64*(2.7346799234073577d+1*t
     ;         27*t50*t65-3.190459910641917d+1*t24*t50*t65)*t84+t62*(3.1
     <         90459910641917d+1*t24*t50*t63-3.6462398978764765d+1*t21*t
     =         50*t63)*t84+t60*(3.6462398978764765d+1*t21*t50*t61-4.1020
     >         19885111037d+1*t18*t50*t61)*t84+t58*(4.102019885111037d+1
     ?         *t18*t50*t59-4.5577998723455965d+1*t15*t50*t59)*t84)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 4.492769220958891d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = 3.6121864536509496d-3*gammaaa*t1*t3+1.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = param(12)
            t7 = rhoa**1.6666666666666669d+0
            t8 = 1/taua
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(11)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(10)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(9)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(8)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(7)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(6)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(5)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(4)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(3)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(2)
            t42 = 1/t12
            t43 = t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37+t32*t33
     1         *t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*t22+t17*
     2         t18*t19+t14*t15*t16+param(1)
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
            t6 = param(12)
            t7 = rhob**1.6666666666666669d+0
            t8 = 1/taub
            t9 = 4.557799872345596d+0*t7*t8
            t10 = t9-1.0d+0
            t11 = t10**11
            t12 = t9+1.0d+0
            t13 = 1/t12**11
            t14 = param(11)
            t15 = t10**10
            t16 = 1/t12**10
            t17 = param(10)
            t18 = t10**9
            t19 = 1/t12**9
            t20 = param(9)
            t21 = t10**8
            t22 = 1/t12**8
            t23 = param(8)
            t24 = t10**7
            t25 = 1/t12**7
            t26 = param(7)
            t27 = t10**6
            t28 = 1/t12**6
            t29 = param(6)
            t30 = t10**5
            t31 = 1/t12**5
            t32 = param(5)
            t33 = t10**4
            t34 = 1/t12**4
            t35 = param(4)
            t36 = t10**3
            t37 = 1/t12**3
            t38 = param(3)
            t39 = t10**2
            t40 = 1/t12**2
            t41 = param(2)
            t42 = 1/t12
            t43 = t11*t13*t6+t10*t41*t42+t38*t39*t40+t35*t36*t37+t32*t33
     1         *t34+t29*t30*t31+t26*t27*t28+t23*t24*t25+t20*t21*t22+t17*
     2         t18*t19+t14*t15*t16+param(1)
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
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-1.0\\\\ 
C>   {\it t_6} &=& {\it t_4}+1.0\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {\it t_{18}}+1.0\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-1.0\\\\ 
C>   {\it t_{24}} &=& {\it t_{22}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{17}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{17}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{\it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{\it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{\it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{\it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{\it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{\it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{\it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{\it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{\it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{\it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{\it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{\it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{\it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{\it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{\it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{\it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}\over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}+{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{\it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{\it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{21}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{21}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{\it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{\it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{\it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{\it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{\it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{\it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{\it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{\it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{\it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{\it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{\it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
      subroutine nwxcm_x_m06_2x_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_m06_2x functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_4} &=& {{4.557799872345596\,\rho_\alpha^{{{5}\over{3}}}}\over{\tau_\alpha}}\\\\ 
C>   {\it t_5} &=& {\it t_4}-1.0\\\\ 
C>   {\it t_6} &=& {\it t_4}+1.0\\\\ 
C>   {\it t_7} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{16}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{4.557799872345596\,\rho_\beta^{{{5}\over{3}}}}\over{\tau_\beta}}\\\\ 
C>   {\it t_{19}} &=& {\it t_{18}}-1.0\\\\ 
C>   {\it t_{20}} &=& {\it t_{18}}+1.0\\\\ 
C>   {\it t_{21}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{22}} &=& {{4.557799872345596\,\rho_s^{{{5}\over{3}}}}\over{\tau_s}}\\\\ 
C>   {\it t_{23}} &=& {\it t_{22}}-1.0\\\\ 
C>   {\it t_{24}} &=& {\it t_{22}}+1.0\\\\ 
C>   f &=& -.9305257363490997\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{17}}\,\sigma_{\beta\beta}}\over{.004492769220958891\,{\it t_{17}}\,\sigma_{\beta\beta}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{19}}}\over{{\it t_{20}}}}+{{{\it t_{15}}\,{\it t_{19}}^2}\over{{\it t_{20}}^2}}+{{{\it t_{14}}\,{\it t_{19}}^3}\over{{\it t_{20}}^3}}+{{{\it t_{13}}\,{\it t_{19}}^4}\over{{\it t_{20}}^4}}+{{{\it t_{12}}\,{\it t_{19}}^5}\over{{\it t_{20}}^5}}+{{{\it t_{11}}\,{\it t_{19}}^6}\over{{\it t_{20}}^6}}+{{{\it t_{10}}\,{\it t_{19}}^7}\over{{\it t_{20}}^7}}+{{{\it t_9}\,{\it t_{19}}^8}\over{{\it t_{20}}^8}}+{{{\it t_8}\,{\it t_{19}}^9}\over{{\it t_{20}}^9}}+{{{\it t_7}\,{\it t_{19}}^{10}}\over{{\it t_{20}}^{10}}}+{{{\it t_3}\,{\it t_{19}}^{11}}\over{{\it t_{20}}^{11}}}+{\it t_2}\right)-.9305257363490997\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004492769220958891\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_5}}\over{{\it t_6}}}+{{{\it t_{15}}\,{\it t_5}^2}\over{{\it t_6}^2}}+{{{\it t_{14}}\,{\it t_5}^3}\over{{\it t_6}^3}}+{{{\it t_{13}}\,{\it t_5}^4}\over{{\it t_6}^4}}+{{{\it t_{12}}\,{\it t_5}^5}\over{{\it t_6}^5}}+{{{\it t_{11}}\,{\it t_5}^6}\over{{\it t_6}^6}}+{{{\it t_{10}}\,{\it t_5}^7}\over{{\it t_6}^7}}+{{{\it t_9}\,{\it t_5}^8}\over{{\it t_6}^8}}+{{{\it t_8}\,{\it t_5}^9}\over{{\it t_6}^9}}+{{{\it t_7}\,{\it t_5}^{10}}\over{{\it t_6}^{10}}}+{{{\it t_3}\,{\it t_5}^{11}}\over{{\it t_6}^{11}}}+{\it t_2}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,\rho_s^{{{4}\over{3}}}\,\left({{0.00361218645365095\,{\it t_{21}}\,\sigma_{ss}}\over{.004492769220958891\,{\it t_{21}}\,\sigma_{ss}+1.0}}+1.0\right)\,\left({{{\it t_{16}}\,{\it t_{23}}}\over{{\it t_{24}}}}+{{{\it t_{15}}\,{\it t_{23}}^2}\over{{\it t_{24}}^2}}+{{{\it t_{14}}\,{\it t_{23}}^3}\over{{\it t_{24}}^3}}+{{{\it t_{13}}\,{\it t_{23}}^4}\over{{\it t_{24}}^4}}+{{{\it t_{12}}\,{\it t_{23}}^5}\over{{\it t_{24}}^5}}+{{{\it t_{11}}\,{\it t_{23}}^6}\over{{\it t_{24}}^6}}+{{{\it t_{10}}\,{\it t_{23}}^7}\over{{\it t_{24}}^7}}+{{{\it t_9}\,{\it t_{23}}^8}\over{{\it t_{24}}^8}}+{{{\it t_8}\,{\it t_{23}}^9}\over{{\it t_{24}}^9}}+{{{\it t_7}\,{\it t_{23}}^{10}}\over{{\it t_{24}}^{10}}}+{{{\it t_3}\,{\it t_{23}}^{11}}\over{{\it t_{24}}^{11}}}+{\it t_2}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 123, 194101 (2006)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2370993 "> 
C> 10.1063/1.2370993 </a> 
C> 
      subroutine nwxcm_x_m06_2x_d3(param,tol_rho,ipol,nq,wght, 
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
