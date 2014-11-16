C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_vwn5.F 
C> The nwxcm_c_vwn5 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_vwn5 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{6.15199081975908}\over{{\it t_4}+3.72744}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 2.935818660072219\,{\it t_3}+{\it t_8}+12.9352\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.10498\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& \arctan \left({{7.123108917818118}\over{{\it t_4}+1.13107}}\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\left| .8908571061768626\,{\it t_3}+{\it t_8}+13.0045\right| }}\\\\ 
C>   {\it t_{18}} &=& \log \left(1.0\,\left({\it t_6}+0.0047584\right)^2\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{19}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_3}+{\it t_8}+18.0578\right| }}\\\\ 
C>   {\it t_{21}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}^{.3333333333333333}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\sqrt{{\it t_{22}}}}}\\\\ 
C>   {\it t_{24}} &=& 1.575246635799487\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& \arctan \left({{6.15199081975908}\over{{\it t_{24}}+3.72744}}\right)\\\\ 
C>   {\it t_{26}} &=& .7876233178997433\,{\it t_{23}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{{\it t_{22}}}}\\\\ 
C>   {\it t_{28}} &=& .6203504908994001\,{\it t_{27}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 2.935818660072219\,{\it t_{23}}+{\it t_{28}}+12.9352\right| }}\\\\ 
C>   {\it t_{30}} &=& \log \left(1.0\,\left({\it t_{26}}+0.10498\right)^2\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{31}} &=& \left| {\it t_{27}}\right| \\\\ 
C>   {\it t_{32}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{33}} &=& \rho_s+{\it t_{13}}\\\\ 
C>   {\it t_{34}} &=& \arctan \left({{7.123108917818118}\over{{\it t_{24}}+1.13107}}\right)\\\\ 
C>   {\it t_{35}} &=& {{1}\over{\left| .8908571061768626\,{\it t_{23}}+{\it t_{28}}+13.0045\right| }}\\\\ 
C>   {\it t_{36}} &=& \log \left(1.0\,\left({\it t_{26}}+0.0047584\right)^2\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{37}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{38}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_{23}}+{\it t_{28}}+18.0578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.125\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{14}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{20}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_6}+0.325\right)^2\,{\it t_{20}}\right)-0.0310907\,{\it t_{12}}-9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_4}+7.06042}}\right)-.03878329487811301\,{\it t_5}\right)+.01688686394038963\,{\it t_{19}}+6.991732350764465 \times 10^{-6}\,{\it t_{18}}+.005365091848883603\,{\it t_{16}}\right)}\over{{\it t_1}^4}}-.01688686394038963\,{\it t_{19}}-6.991732350764465 \times 10^{-6}\,{\it t_{18}}-.005365091848883603\,{\it t_{16}}\right)+0.0310907\,{\it t_{12}}+9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.03878329487811301\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{21}}\,\left(1.125\,\left(\left({{{\it t_{33}}}\over{{\it t_{21}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{33}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{38}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_{26}}+0.325\right)^2\,{\it t_{38}}\right)-0.0310907\,{\it t_{32}}-9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_{24}}+7.06042}}\right)-.03878329487811301\,{\it t_{25}}\right)+.01688686394038963\,{\it t_{37}}+6.991732350764465 \times 10^{-6}\,{\it t_{36}}+.005365091848883603\,{\it t_{34}}\right)}\over{{\it t_{21}}^4}}-.01688686394038963\,{\it t_{37}}-6.991732350764465 \times 10^{-6}\,{\it t_{36}}-.005365091848883603\,{\it t_{34}}\right)+0.0310907\,{\it t_{32}}+9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.03878329487811301\,{\it t_{25}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn5(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = 1/rhoa**1.6666666666666666d-1
            t2 = 1.4033852085886447d+0*t1+3.72744d+0
            t3 = 7.016926042943223d-1*t1+1.0498d-1
            t4 = t3**2
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t5
            t7 = 4.9237251092134837d-1*t6
            t8 = t7+2.615517080950829d+0*t1+1.29352d+1
            t9 = 1/t8
            t10 = t7+2.615517080950829d+0/t5**5.0d-1+1.29352d+1
            t11 = 1/t10
            t12 = 9.690227711544375d-4*log(1.0d+0*t4*t9)+3.10907d-2*log(
     1         4.9237251092134837d-1*t11*t6)+3.8783294878113006d-2*atan(
     2         6.15199081975908d+0/t2)
            t13 = 1/t2**2
            t14 = 1/rhoa**1.1666666666666667d+0
            t15 = 1/rhoa**1.3333333333333333d+0
            t16 = -1.6412417030711612d-1*t15-4.359195134918048d-1*t14
            fnc(iq) = 2.0d+0*rhoa*t12*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*(9.690227711544375d-4*t8*(-2.33
     1         8975347647741d-1*t14*t3*t9-1.0d+0*t16*t4/t8**2)/t4+6.3144
     2         67057030003d-2*t10*t5*(-4.9237251092134837d-1*t16*t6/t10*
     3         *2-1.6412417030711612d-1*t11*t15)+5.580665928883042d-2*t1
     4         3*t14/(3.784699104639999d+1*t13+1))*wght+2.0d+0*t12*wght+
     5         Amat(iq,D1_RA)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob+rhoa
            t2 = t1**3.333333333333333d-1
            t3 = t2**5.0d-1
            t4 = 1/t3
            t5 = 1.5752466357994865d+0*t4
            t6 = t5+3.72744d+0
            t7 = atan(6.15199081975908d+0/t6)
            t8 = 7.876233178997433d-1*t4
            t9 = t8+1.0498d-1
            t10 = t9**2
            t11 = 1/t2
            t12 = 6.203504908994001d-1*t11
            t13 = 2.935818660072219d+0*t4+t12+1.29352d+1
            t14 = abs(t13)
            t15 = 1/t14
            t16 = log(1.0d+0*t10*t15)
            t17 = abs(t11)
            t18 = log(6.203504908994001d-1*t15*t17)
            t19 = rhoa-rhob
            t20 = 1/t1
            t21 = 1.0d+0-t19*t20
            t22 = t19*t20+1.0d+0
            t23 = t22**1.3333333333333333d+0+t21**1.3333333333333333d+0-
     1         2.0d+0
            t24 = t5+1.13107d+0
            t25 = atan(7.123108917818117d+0/t24)
            t26 = t8+4.7584000000000004d-3
            t27 = t26**2
            t28 = 8.908571061768626d-1*t4+t12+1.30045d+1
            t29 = abs(t28)
            t30 = 1/t29
            t31 = log(1.0d+0*t27*t30)
            t32 = log(6.203504908994001d-1*t17*t30)
            t33 = t19**4
            t34 = 1/t1**4
            t35 = t5+7.06042d+0
            t36 = t8+3.25d-1
            t37 = t36**2
            t38 = 5.560951426165704d+0*t4+t12+1.8057800000000002d+1
            t39 = abs(t38)
            t40 = 1/t39
            t41 = 1.709920934161365d+0*(-3.8783294878113006d-2*t7+2.2478
     1         670955426108d-3*log(1.0d+0*t37*t40)+1.554535d-2*log(6.203
     2         504908994001d-1*t17*t40)+5.249139316978092d-2*atan(4.7309
     3         269095601136d+0/t35)-3.10907d-2*t18-9.690227711544375d-4*
     4         t16)+1.688686394038963d-2*t32+6.991732350764465d-6*t31+5.
     5         365091848883603d-3*t25
            t42 = t33*t34*t41-1.688686394038963d-2*t32-6.991732350764465
     1         d-6*t31-5.365091848883603d-3*t25
            t43 = 3.8783294878113006d-2*t7+1.1250000000000002d+0*t23*t42
     1         +3.10907d-2*t18+9.690227711544375d-4*t16
            t44 = 1/t1**6.666666666666667d-1
            t45 = 1/t3**3
            t46 = 1/t6**2
            t47 = 1/(3.784699104639999d+1*t46+1)
            t48 = 6.264085709464389d-2*t44*t45*t46*t47
            t49 = 1/t10
            t50 = -2.067834969664667d-1/t1**1.3333333333333333d+0
            t51 = t50-4.8930311001203647d-1*t44*t45
            t52 = 1/t13
            t53 = -2.625411059665811d-1*t15*t44*t45*t9-1.0d+0*t10*t15*t5
     1         1*t52
            t54 = 9.690227711544375d-4*t14*t49*t53
            t55 = 1/t17
            t56 = 1/t1**1.6666666666666666d+0
            t57 = -2.067834969664667d-1*t15*t55*t56-6.203504908994001d-1
     1         *t15*t17*t51*t52
            t58 = 5.011795824473985d-2*t14*t55*t57
            t59 = 1/t24**2
            t60 = 1/(5.07386806551d+1*t59+1)
            t61 = -1.0033305979455584d-2*t44*t45*t59*t60
            t62 = 1/t27
            t63 = t50-1.4847618436281043d-1*t44*t45
            t64 = 1/t28
            t65 = -1.0d+0*t27*t30*t63*t64-2.625411059665811d-1*t26*t30*t
     1         44*t45
            t66 = -6.991732350764465d-6*t29*t62*t65
            t67 = -6.203504908994001d-1*t17*t30*t63*t64-2.06783496966466
     1         7d-1*t30*t55*t56
            t68 = -2.722148880047894d-2*t29*t55*t67
            t69 = 1/t35**2
            t70 = t50-9.26825237694284d-1*t44*t45
            t71 = 1/t38
            t72 = t33*t34*(1.709920934161365d+0*(2.2478670955426108d-3*t
     1         39*(-1.0d+0*t37*t40*t70*t71-2.625411059665811d-1*t36*t40*
     2         t44*t45)/t37+2.5058979122369923d-2*t39*t55*(-6.2035049089
     3         94001d-1*t17*t40*t70*t71-2.067834969664667d-1*t40*t55*t56
     4         )+6.519760588836059d-2*t44*t45*t69/(2.238166942360001d+1*
     5         t69+1)-5.011795824473985d-2*t14*t55*t57-9.690227711544375
     6         d-4*t14*t49*t53-6.264085709464389d-2*t44*t45*t46*t47)+2.7
     7         22148880047894d-2*t29*t55*t67+6.991732350764465d-6*t29*t6
     8         2*t65+1.0033305979455584d-2*t44*t45*t59*t60)
            t73 = -4*t33*t41/t1**5
            t74 = t19**3
            t75 = 1/t1**2
            t76 = t19*t75
            t77 = -t20
            t78 = t21**3.333333333333333d-1
            t79 = -t19*t75
            t80 = t22**3.333333333333333d-1
            t81 = 1.0d+0*t43*wght
            fnc(iq) = 1.0d+0*t1*t43*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.1250000000000002d+0*t42*(1.333
     1         3333333333333d+0*(t79+t20)*t80+1.3333333333333333d+0*(t77
     2         +t76)*t78)+1.1250000000000002d+0*t23*(4*t34*t41*t74+t73+t
     3         72+t68+t66+t61)+t58+t54+t48)*wght+t81+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*t1*(1.1250000000000002d+0*t42*(1.333
     1         3333333333333d+0*(t79+t77)*t80+1.3333333333333333d+0*(t76
     2         +t20)*t78)+1.1250000000000002d+0*t23*(-4*t34*t41*t74+t73+
     3         t72+t68+t66+t61)+t58+t54+t48)*wght+t81+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhob+rhoa
            t2 = t1**3.333333333333333d-1
            t3 = t2**5.0d-1
            t4 = 1/t3
            t5 = 1.5752466357994865d+0*t4
            t6 = t5+3.72744d+0
            t7 = atan(6.15199081975908d+0/t6)
            t8 = 7.876233178997433d-1*t4
            t9 = t8+1.0498d-1
            t10 = t9**2
            t11 = 1/t2
            t12 = 6.203504908994001d-1*t11
            t13 = 2.935818660072219d+0*t4+t12+1.29352d+1
            t14 = abs(t13)
            t15 = 1/t14
            t16 = log(1.0d+0*t10*t15)
            t17 = abs(t11)
            t18 = log(6.203504908994001d-1*t15*t17)
            t19 = rhoa-rhob
            t20 = 1/t1
            t21 = t19*t20+1.0d+0
            t22 = t21**1.3333333333333333d+0-2.0d+0
            t23 = t5+1.13107d+0
            t24 = atan(7.123108917818117d+0/t23)
            t25 = t8+4.7584000000000004d-3
            t26 = t25**2
            t27 = 8.908571061768626d-1*t4+t12+1.30045d+1
            t28 = abs(t27)
            t29 = 1/t28
            t30 = log(1.0d+0*t26*t29)
            t31 = log(6.203504908994001d-1*t17*t29)
            t32 = t19**4
            t33 = 1/t1**4
            t34 = t5+7.06042d+0
            t35 = t8+3.25d-1
            t36 = t35**2
            t37 = 5.560951426165704d+0*t4+t12+1.8057800000000002d+1
            t38 = abs(t37)
            t39 = 1/t38
            t40 = 1.709920934161365d+0*(-3.8783294878113006d-2*t7+2.2478
     1         670955426108d-3*log(1.0d+0*t36*t39)+1.554535d-2*log(6.203
     2         504908994001d-1*t17*t39)+5.249139316978092d-2*atan(4.7309
     3         269095601136d+0/t34)-3.10907d-2*t18-9.690227711544375d-4*
     4         t16)+1.688686394038963d-2*t31+6.991732350764465d-6*t30+5.
     5         365091848883603d-3*t24
            t41 = t32*t33*t40-1.688686394038963d-2*t31-6.991732350764465
     1         d-6*t30-5.365091848883603d-3*t24
            t42 = 3.8783294878113006d-2*t7+1.1250000000000002d+0*t22*t41
     1         +3.10907d-2*t18+9.690227711544375d-4*t16
            t43 = 1/t1**6.666666666666667d-1
            t44 = 1/t3**3
            t45 = 1/t6**2
            t46 = 1/(3.784699104639999d+1*t45+1)
            t47 = 1/t10
            t48 = -2.067834969664667d-1/t1**1.3333333333333333d+0
            t49 = t48-4.8930311001203647d-1*t43*t44
            t50 = 1/t13
            t51 = -2.625411059665811d-1*t15*t43*t44*t9-1.0d+0*t10*t15*t4
     1         9*t50
            t52 = 1/t17
            t53 = 1/t1**1.6666666666666666d+0
            t54 = -2.067834969664667d-1*t15*t52*t53-6.203504908994001d-1
     1         *t15*t17*t49*t50
            t55 = 1/t23**2
            t56 = 1/(5.07386806551d+1*t55+1)
            t57 = 1/t26
            t58 = t48-1.4847618436281043d-1*t43*t44
            t59 = 1/t27
            t60 = -1.0d+0*t26*t29*t58*t59-2.625411059665811d-1*t25*t29*t
     1         43*t44
            t61 = -6.203504908994001d-1*t17*t29*t58*t59-2.06783496966466
     1         7d-1*t29*t52*t53
            t62 = 1/t34**2
            t63 = t48-9.26825237694284d-1*t43*t44
            t64 = 1/t37
            fnc(iq) = 1.0d+0*t1*t42*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.1250000000000002d+0*t22*(t32*t
     1         33*(1.709920934161365d+0*(2.2478670955426108d-3*t38*(-1.0
     2         d+0*t36*t39*t63*t64-2.625411059665811d-1*t35*t39*t43*t44)
     3         /t36+2.5058979122369923d-2*t38*t52*(-6.203504908994001d-1
     4         *t17*t39*t63*t64-2.067834969664667d-1*t39*t52*t53)+6.5197
     5         60588836059d-2*t43*t44*t62/(2.238166942360001d+1*t62+1)-5
     6         .011795824473985d-2*t14*t52*t54-9.690227711544375d-4*t14*
     7         t47*t51-6.264085709464389d-2*t43*t44*t45*t46)+2.722148880
     8         047894d-2*t28*t52*t61+6.991732350764465d-6*t28*t57*t60+1.
     9         0033305979455584d-2*t43*t44*t55*t56)-2.722148880047894d-2
     :         *t28*t52*t61-6.991732350764465d-6*t28*t57*t60-1.003330597
     ;         9455584d-2*t43*t44*t55*t56+4*t19**3*t33*t40-4*t32*t40/t1*
     <         *5)+5.011795824473985d-2*t14*t52*t54+9.690227711544375d-4
     =         *t14*t47*t51+6.264085709464389d-2*t43*t44*t45*t46+1.50000
     >         00000000002d+0*(t20-t19/t1**2)*t21**3.333333333333333d-1*
     ?         t41)*wght+1.0d+0*t42*wght+Amat(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**1.6666666666666666d-1
            t2 = 1.4033852085886447d+0*t1
            t3 = t2+3.72744d+0
            t4 = 7.016926042943223d-1*t1
            t5 = t4+1.0498d-1
            t6 = t5**2
            t7 = rhob**3.333333333333333d-1
            t8 = 1/t7
            t9 = 4.9237251092134837d-1*t8
            t10 = t9+2.615517080950829d+0*t1+1.29352d+1
            t11 = 1/t10
            t12 = t9+2.615517080950829d+0/t7**5.0d-1+1.29352d+1
            t13 = 1/t12
            t14 = t4+4.7584000000000004d-3
            t15 = t14**2
            t16 = t9+7.936634539391791d-1*t1+1.30045d+1
            t17 = 1/t16
            t18 = t2+1.13107d+0
            t19 = -1.1250000000000002d+0*(-1.688686394038963d-2*log(4.92
     1         37251092134837d-1*t17*t8)-5.365091848883603d-3*atan(7.123
     2         108917818117d+0/t18)-6.991732350764465d-6*log(1.0d+0*t15*
     3         t17))+3.10907d-2*log(4.9237251092134837d-1*t13*t8)+9.6902
     4         27711544375d-4*log(1.0d+0*t11*t6)+3.8783294878113006d-2*a
     5         tan(6.15199081975908d+0/t3)
            t20 = 1/t3**2
            t21 = 1/rhob**1.1666666666666667d+0
            t22 = 1/rhob**1.3333333333333333d+0
            t23 = -1.6412417030711612d-1*t22
            t24 = t23-4.359195134918048d-1*t21
            t25 = 1/t18**2
            t26 = t23-1.3227724232319651d-1*t21
            t27 = 1/t16**2
            fnc(iq) = 2.0d+0*rhob*t19*wght+fnc(iq)
            Amat(iq,D1_RB) = 2.0d+0*rhob*(-1.1250000000000002d+0*(-3.429
     1         692674920095d-2*t16*t7*(-4.9237251092134837d-1*t26*t27*t8
     2         -1.6412417030711612d-1*t17*t22)-6.991732350764465d-6*t16*
     3         (-1.0d+0*t15*t26*t27-2.338975347647741d-1*t14*t17*t21)/t1
     4         5-8.938659435806782d-3*t21*t25/(5.07386806551d+1*t25+1))+
     5         6.314467057030003d-2*t12*t7*(-4.9237251092134837d-1*t24*t
     6         8/t12**2-1.6412417030711612d-1*t13*t22)+9.690227711544375
     7         d-4*t10*(-1.0d+0*t24*t6/t10**2-2.338975347647741d-1*t11*t
     8         21*t5)/t6+5.580665928883042d-2*t20*t21/(3.784699104639999
     9         d+1*t20+1))*wght+2.0d+0*t19*wght+Amat(iq,D1_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_vwn5 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{6.15199081975908}\over{{\it t_4}+3.72744}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 2.935818660072219\,{\it t_3}+{\it t_8}+12.9352\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.10498\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& \arctan \left({{7.123108917818118}\over{{\it t_4}+1.13107}}\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\left| .8908571061768626\,{\it t_3}+{\it t_8}+13.0045\right| }}\\\\ 
C>   {\it t_{18}} &=& \log \left(1.0\,\left({\it t_6}+0.0047584\right)^2\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{19}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_3}+{\it t_8}+18.0578\right| }}\\\\ 
C>   {\it t_{21}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}^{.3333333333333333}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\sqrt{{\it t_{22}}}}}\\\\ 
C>   {\it t_{24}} &=& 1.575246635799487\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& \arctan \left({{6.15199081975908}\over{{\it t_{24}}+3.72744}}\right)\\\\ 
C>   {\it t_{26}} &=& .7876233178997433\,{\it t_{23}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{{\it t_{22}}}}\\\\ 
C>   {\it t_{28}} &=& .6203504908994001\,{\it t_{27}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 2.935818660072219\,{\it t_{23}}+{\it t_{28}}+12.9352\right| }}\\\\ 
C>   {\it t_{30}} &=& \log \left(1.0\,\left({\it t_{26}}+0.10498\right)^2\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{31}} &=& \left| {\it t_{27}}\right| \\\\ 
C>   {\it t_{32}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{33}} &=& \rho_s+{\it t_{13}}\\\\ 
C>   {\it t_{34}} &=& \arctan \left({{7.123108917818118}\over{{\it t_{24}}+1.13107}}\right)\\\\ 
C>   {\it t_{35}} &=& {{1}\over{\left| .8908571061768626\,{\it t_{23}}+{\it t_{28}}+13.0045\right| }}\\\\ 
C>   {\it t_{36}} &=& \log \left(1.0\,\left({\it t_{26}}+0.0047584\right)^2\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{37}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{38}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_{23}}+{\it t_{28}}+18.0578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.125\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{14}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{20}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_6}+0.325\right)^2\,{\it t_{20}}\right)-0.0310907\,{\it t_{12}}-9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_4}+7.06042}}\right)-.03878329487811301\,{\it t_5}\right)+.01688686394038963\,{\it t_{19}}+6.991732350764465 \times 10^{-6}\,{\it t_{18}}+.005365091848883603\,{\it t_{16}}\right)}\over{{\it t_1}^4}}-.01688686394038963\,{\it t_{19}}-6.991732350764465 \times 10^{-6}\,{\it t_{18}}-.005365091848883603\,{\it t_{16}}\right)+0.0310907\,{\it t_{12}}+9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.03878329487811301\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{21}}\,\left(1.125\,\left(\left({{{\it t_{33}}}\over{{\it t_{21}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{33}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{38}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_{26}}+0.325\right)^2\,{\it t_{38}}\right)-0.0310907\,{\it t_{32}}-9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_{24}}+7.06042}}\right)-.03878329487811301\,{\it t_{25}}\right)+.01688686394038963\,{\it t_{37}}+6.991732350764465 \times 10^{-6}\,{\it t_{36}}+.005365091848883603\,{\it t_{34}}\right)}\over{{\it t_{21}}^4}}-.01688686394038963\,{\it t_{37}}-6.991732350764465 \times 10^{-6}\,{\it t_{36}}-.005365091848883603\,{\it t_{34}}\right)+0.0310907\,{\it t_{32}}+9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.03878329487811301\,{\it t_{25}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn5_d2(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_vwn5 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{6.15199081975908}\over{{\it t_4}+3.72744}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 2.935818660072219\,{\it t_3}+{\it t_8}+12.9352\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.10498\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& \arctan \left({{7.123108917818118}\over{{\it t_4}+1.13107}}\right)\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\left| .8908571061768626\,{\it t_3}+{\it t_8}+13.0045\right| }}\\\\ 
C>   {\it t_{18}} &=& \log \left(1.0\,\left({\it t_6}+0.0047584\right)^2\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{19}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{17}}\right)\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_3}+{\it t_8}+18.0578\right| }}\\\\ 
C>   {\it t_{21}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{22}} &=& {\it t_{21}}^{.3333333333333333}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\sqrt{{\it t_{22}}}}}\\\\ 
C>   {\it t_{24}} &=& 1.575246635799487\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& \arctan \left({{6.15199081975908}\over{{\it t_{24}}+3.72744}}\right)\\\\ 
C>   {\it t_{26}} &=& .7876233178997433\,{\it t_{23}}\\\\ 
C>   {\it t_{27}} &=& {{1}\over{{\it t_{22}}}}\\\\ 
C>   {\it t_{28}} &=& .6203504908994001\,{\it t_{27}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 2.935818660072219\,{\it t_{23}}+{\it t_{28}}+12.9352\right| }}\\\\ 
C>   {\it t_{30}} &=& \log \left(1.0\,\left({\it t_{26}}+0.10498\right)^2\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{31}} &=& \left| {\it t_{27}}\right| \\\\ 
C>   {\it t_{32}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{29}}\right)\\\\ 
C>   {\it t_{33}} &=& \rho_s+{\it t_{13}}\\\\ 
C>   {\it t_{34}} &=& \arctan \left({{7.123108917818118}\over{{\it t_{24}}+1.13107}}\right)\\\\ 
C>   {\it t_{35}} &=& {{1}\over{\left| .8908571061768626\,{\it t_{23}}+{\it t_{28}}+13.0045\right| }}\\\\ 
C>   {\it t_{36}} &=& \log \left(1.0\,\left({\it t_{26}}+0.0047584\right)^2\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{37}} &=& \log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{35}}\right)\\\\ 
C>   {\it t_{38}} &=& {{1}\over{\left| 5.560951426165705\,{\it t_{23}}+{\it t_{28}}+18.0578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.125\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{14}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{20}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_6}+0.325\right)^2\,{\it t_{20}}\right)-0.0310907\,{\it t_{12}}-9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_4}+7.06042}}\right)-.03878329487811301\,{\it t_5}\right)+.01688686394038963\,{\it t_{19}}+6.991732350764465 \times 10^{-6}\,{\it t_{18}}+.005365091848883603\,{\it t_{16}}\right)}\over{{\it t_1}^4}}-.01688686394038963\,{\it t_{19}}-6.991732350764465 \times 10^{-6}\,{\it t_{18}}-.005365091848883603\,{\it t_{16}}\right)+0.0310907\,{\it t_{12}}+9.690227711544375 \times 10^{-4}\,{\it t_{10}}+.03878329487811301\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{21}}\,\left(1.125\,\left(\left({{{\it t_{33}}}\over{{\it t_{21}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left({{{\it t_{33}}^4\,\left(1.709920934161365\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{31}}\,{\it t_{38}}\right)+.002247867095542611\,\log \left(1.0\,\left({\it t_{26}}+0.325\right)^2\,{\it t_{38}}\right)-0.0310907\,{\it t_{32}}-9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.05249139316978092\,\arctan \left({{4.730926909560114}\over{{\it t_{24}}+7.06042}}\right)-.03878329487811301\,{\it t_{25}}\right)+.01688686394038963\,{\it t_{37}}+6.991732350764465 \times 10^{-6}\,{\it t_{36}}+.005365091848883603\,{\it t_{34}}\right)}\over{{\it t_{21}}^4}}-.01688686394038963\,{\it t_{37}}-6.991732350764465 \times 10^{-6}\,{\it t_{36}}-.005365091848883603\,{\it t_{34}}\right)+0.0310907\,{\it t_{32}}+9.690227711544375 \times 10^{-4}\,{\it t_{30}}+.03878329487811301\,{\it t_{25}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn5_d3(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2,Amat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho 
c 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
