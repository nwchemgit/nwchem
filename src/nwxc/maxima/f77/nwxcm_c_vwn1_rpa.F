C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_vwn1_rpa.F 
C> The nwxcm_c_vwn1_rpa functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_vwn1_rpa functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{0.0448998886415768}\over{{\it t_4}+13.072}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 10.29581201158544\,{\it t_3}+{\it t_8}+42.7198\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.409286\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_3}+{\it t_8}+101.578\right| }}\\\\ 
C>   {\it t_{17}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}^{.3333333333333333}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\sqrt{{\it t_{18}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.575246635799487\,{\it t_{19}}\\\\ 
C>   {\it t_{21}} &=& \arctan \left({{0.0448998886415768}\over{{\it t_{20}}+13.072}}\right)\\\\ 
C>   {\it t_{22}} &=& .7876233178997433\,{\it t_{19}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{18}}}}\\\\ 
C>   {\it t_{24}} &=& .6203504908994001\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\left| 10.29581201158544\,{\it t_{19}}+{\it t_{24}}+42.7198\right| }}\\\\ 
C>   {\it t_{26}} &=& \log \left(1.0\,\left({\it t_{22}}+0.409286\right)^2\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{27}} &=& \left| {\it t_{23}}\right| \\\\ 
C>   {\it t_{28}} &=& \log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_{19}}+{\it t_{24}}+101.578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{16}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_6}+0.743294\right)^2\,{\it t_{16}}\right)-0.0310907\,{\it t_{12}}-.004431373767749538\,{\it t_{10}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_4}+20.1231}}\right)-20.52197293770518\,{\it t_5}\right)+0.0310907\,{\it t_{12}}+.004431373767749538\,{\it t_{10}}+20.52197293770518\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{17}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{13}}}\over{{\it t_{17}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{29}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_{22}}+0.743294\right)^2\,{\it t_{29}}\right)-0.0310907\,{\it t_{28}}-.004431373767749538\,{\it t_{26}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_{20}}+20.1231}}\right)-20.52197293770518\,{\it t_{21}}\right)+0.0310907\,{\it t_{28}}+.004431373767749538\,{\it t_{26}}+20.52197293770518\,{\it t_{21}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn1_rpa(param,tol_rho,ipol,nq,wght, 
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
            t2 = 1.4033852085886447d+0*t1+1.3072d+1
            t3 = 7.016926042943223d-1*t1+4.09286d-1
            t4 = t3**2
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/t5
            t7 = 4.9237251092134837d-1*t6
            t8 = t7+9.172525723335381d+0*t1+4.2719799999999997d+1
            t9 = 1/t8
            t10 = t7+9.172525723335381d+0/t5**5.0d-1+4.2719799999999997d
     1         +1
            t11 = 1/t10
            t12 = 4.431373767749538d-3*log(1.0d+0*t4*t9)+3.10907d-2*log(
     1         4.9237251092134837d-1*t11*t6)+2.052197293770518d+1*atan(4
     2         .4899888641576796d-2/t2)
            t13 = 1/t2**2
            t14 = 1/rhoa**1.1666666666666667d+0
            t15 = 1/rhoa**1.3333333333333333d+0
            t16 = -1.6412417030711612d-1*t15-1.5287542872225635d+0*t14
            fnc(iq) = 2.0d+0*rhoa*t12*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*rhoa*(4.431373767749538d-3*t8*(-2.33
     1         8975347647741d-1*t14*t3*t9-1.0d+0*t16*t4/t8**2)/t4+6.3144
     2         67057030003d-2*t10*t5*(-4.9237251092134837d-1*t16*t6/t10*
     3         *2-1.6412417030711612d-1*t11*t15)+2.1552121112611455d-1*t
     4         13*t14/(2.016000000025997d-3*t13+1))*wght+2.0d+0*t12*wght
     5         +Amat(iq,D1_RA)
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
            t6 = t5+1.3072d+1
            t7 = atan(4.4899888641576796d-2/t6)
            t8 = 7.876233178997433d-1*t4
            t9 = t8+4.09286d-1
            t10 = t9**2
            t11 = 1/t2
            t12 = 6.203504908994001d-1*t11
            t13 = 1.0295812011585444d+1*t4+t12+4.2719799999999997d+1
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
            t24 = t5+2.0123100000000002d+1
            t25 = t8+7.43294d-1
            t26 = t25**2
            t27 = 1.5849422788428325d+1*t4+t12+1.01578d+2
            t28 = abs(t27)
            t29 = 1/t28
            t30 = -2.052197293770518d+1*t7+2.667310007273315d-3*log(1.0d
     1         +0*t26*t29)+1.554535d-2*log(6.203504908994001d-1*t17*t29)
     2         +6.188180297906176d-1*atan(1.171685277708971d+0/t24)-3.10
     3         907d-2*t18-4.431373767749538d-3*t16
            t31 = 2.052197293770518d+1*t7+1.923661050931536d+0*t23*t30+3
     1         .10907d-2*t18+4.431373767749538d-3*t16
            t32 = 1/t1**6.666666666666667d-1
            t33 = 1/t3**3
            t34 = 1/t6**2
            t35 = 1/(2.016000000025997d-3*t34+1)
            t36 = 2.4191438009473537d-1*t32*t33*t34*t35
            t37 = 1/t10
            t38 = -2.067834969664667d-1/t1**1.3333333333333333d+0
            t39 = t38-1.715968668597574d+0*t32*t33
            t40 = 1/t13
            t41 = -2.625411059665811d-1*t15*t32*t33*t9-1.0d+0*t10*t15*t3
     1         9*t40
            t42 = 4.431373767749538d-3*t14*t37*t41
            t43 = 1/t17
            t44 = 1/t1**1.6666666666666666d+0
            t45 = -2.067834969664667d-1*t15*t43*t44-6.203504908994001d-1
     1         *t15*t17*t39*t40
            t46 = 5.011795824473985d-2*t14*t43*t45
            t47 = 1/t24**2
            t48 = t38-2.641570464738054d+0*t32*t33
            t49 = 1/t27
            t50 = 1.923661050931536d+0*t23*(2.667310007273315d-3*t28*(-1
     1         .0d+0*t26*t29*t48*t49-2.625411059665811d-1*t25*t29*t32*t3
     2         3)/t26+2.5058979122369923d-2*t28*t43*(-6.203504908994001d
     3         -1*t17*t29*t48*t49-2.067834969664667d-1*t29*t43*t44)+1.90
     4         35804775132142d-1*t32*t33*t47/(1.3728463899999485d+0*t47+
     5         1)-5.011795824473985d-2*t14*t43*t45-4.431373767749538d-3*
     6         t14*t37*t41-2.4191438009473537d-1*t32*t33*t34*t35)
            t51 = 1/t1**2
            t52 = t19*t51
            t53 = -t20
            t54 = t21**3.333333333333333d-1
            t55 = -t19*t51
            t56 = t22**3.333333333333333d-1
            t57 = 1.0d+0*t31*wght
            fnc(iq) = 1.0d+0*t1*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.923661050931536d+0*t30*(1.3333
     1         333333333333d+0*(t55+t20)*t56+1.3333333333333333d+0*(t53+
     2         t52)*t54)+t50+t46+t42+t36)*wght+t57+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*t1*(1.923661050931536d+0*t30*(1.3333
     1         333333333333d+0*(t55+t53)*t56+1.3333333333333333d+0*(t52+
     2         t20)*t54)+t50+t46+t42+t36)*wght+t57+Amat(iq,D1_RB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhob+rhoa
            t2 = t1**3.333333333333333d-1
            t3 = t2**5.0d-1
            t4 = 1/t3
            t5 = 1.5752466357994865d+0*t4
            t6 = t5+1.3072d+1
            t7 = atan(4.4899888641576796d-2/t6)
            t8 = 7.876233178997433d-1*t4
            t9 = t8+4.09286d-1
            t10 = t9**2
            t11 = 1/t2
            t12 = 6.203504908994001d-1*t11
            t13 = 1.0295812011585444d+1*t4+t12+4.2719799999999997d+1
            t14 = abs(t13)
            t15 = 1/t14
            t16 = log(1.0d+0*t10*t15)
            t17 = abs(t11)
            t18 = log(6.203504908994001d-1*t15*t17)
            t19 = rhoa-rhob
            t20 = 1/t1
            t21 = t19*t20+1.0d+0
            t22 = t21**1.3333333333333333d+0-2.0d+0
            t23 = t5+2.0123100000000002d+1
            t24 = t8+7.43294d-1
            t25 = t24**2
            t26 = 1.5849422788428325d+1*t4+t12+1.01578d+2
            t27 = abs(t26)
            t28 = 1/t27
            t29 = -2.052197293770518d+1*t7+2.667310007273315d-3*log(1.0d
     1         +0*t25*t28)+1.554535d-2*log(6.203504908994001d-1*t17*t28)
     2         +6.188180297906176d-1*atan(1.171685277708971d+0/t23)-3.10
     3         907d-2*t18-4.431373767749538d-3*t16
            t30 = 2.052197293770518d+1*t7+1.923661050931536d+0*t22*t29+3
     1         .10907d-2*t18+4.431373767749538d-3*t16
            t31 = 1/t1**6.666666666666667d-1
            t32 = 1/t3**3
            t33 = 1/t6**2
            t34 = 1/(2.016000000025997d-3*t33+1)
            t35 = 1/t10
            t36 = -2.067834969664667d-1/t1**1.3333333333333333d+0
            t37 = t36-1.715968668597574d+0*t31*t32
            t38 = 1/t13
            t39 = -2.625411059665811d-1*t15*t31*t32*t9-1.0d+0*t10*t15*t3
     1         7*t38
            t40 = 1/t17
            t41 = 1/t1**1.6666666666666666d+0
            t42 = -2.067834969664667d-1*t15*t40*t41-6.203504908994001d-1
     1         *t15*t17*t37*t38
            t43 = 1/t23**2
            t44 = t36-2.641570464738054d+0*t31*t32
            t45 = 1/t26
            fnc(iq) = 1.0d+0*t1*t30*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t1*(1.923661050931536d+0*t22*(2.6673
     1         10007273315d-3*t27*(-1.0d+0*t25*t28*t44*t45-2.62541105966
     2         5811d-1*t24*t28*t31*t32)/t25+2.5058979122369923d-2*t27*t4
     3         0*(-6.203504908994001d-1*t17*t28*t44*t45-2.06783496966466
     4         7d-1*t28*t40*t41)+1.9035804775132142d-1*t31*t32*t43/(1.37
     5         28463899999485d+0*t43+1)-5.011795824473985d-2*t14*t40*t42
     6         -4.431373767749538d-3*t14*t35*t39-2.4191438009473537d-1*t
     7         31*t32*t33*t34)+5.011795824473985d-2*t14*t40*t42+4.431373
     8         767749538d-3*t14*t35*t39+2.4191438009473537d-1*t31*t32*t3
     9         3*t34+2.564881401242048d+0*(t20-t19/t1**2)*t21**3.3333333
     :         33333333d-1*t29)*wght+1.0d+0*t30*wght+Amat(iq,D1_RA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**1.6666666666666666d-1
            t2 = 1.4033852085886447d+0*t1
            t3 = t2+1.3072d+1
            t4 = atan(4.4899888641576796d-2/t3)
            t5 = 7.016926042943223d-1*t1
            t6 = t5+4.09286d-1
            t7 = t6**2
            t8 = rhob**3.333333333333333d-1
            t9 = 1/t8
            t10 = 4.9237251092134837d-1*t9
            t11 = t10+9.172525723335381d+0*t1+4.2719799999999997d+1
            t12 = 1/t11
            t13 = log(1.0d+0*t12*t7)
            t14 = 9.172525723335381d+0/t8**5.0d-1+t10+4.2719799999999997
     1         d+1
            t15 = 1/t14
            t16 = t2+2.0123100000000002d+1
            t17 = t5+7.43294d-1
            t18 = t17**2
            t19 = t10+1.412023044547508d+1*t1+1.01578d+2
            t20 = 1/t19
            t21 = -1.923661050931536d+0*(1.554535d-2*log(4.9237251092134
     1         837d-1*t20*t9)-3.10907d-2*log(4.9237251092134837d-1*t12*t
     2         9)-2.052197293770518d+1*t4+2.667310007273315d-3*log(1.0d+
     3         0*t18*t20)+6.188180297906176d-1*atan(1.171685277708971d+0
     4         /t16)-4.431373767749538d-3*t13)+3.10907d-2*log(4.92372510
     5         92134837d-1*t15*t9)+2.052197293770518d+1*t4+4.43137376774
     6         9538d-3*t13
            t22 = 1/t3**2
            t23 = 1/(2.016000000025997d-3*t22+1)
            t24 = 1/rhob**1.1666666666666667d+0
            t25 = 1/rhob**1.3333333333333333d+0
            t26 = -1.6412417030711612d-1*t25
            t27 = t26-1.5287542872225635d+0*t24
            t28 = 1/t11**2
            t29 = -1.0d+0*t27*t28*t7-2.338975347647741d-1*t12*t24*t6
            t30 = 1/t7
            t31 = 1/t16**2
            t32 = t26-2.3533717409125127d+0*t24
            t33 = 1/t19**2
            fnc(iq) = 2.0d+0*rhob*t21*wght+fnc(iq)
            Amat(iq,D1_RB) = 2.0d+0*rhob*(-1.923661050931536d+0*(3.15723
     1         35285150016d-2*t19*t8*(-4.9237251092134837d-1*t32*t33*t9-
     2         1.6412417030711612d-1*t20*t25)-6.314467057030003d-2*t11*t
     3         8*(-4.9237251092134837d-1*t27*t28*t9-1.6412417030711612d-
     4         1*t12*t25)+2.667310007273315d-3*t19*(-1.0d+0*t18*t32*t33-
     5         2.338975347647741d-1*t17*t20*t24)/t18+1.695897407293498d-
     6         1*t24*t31/(1.3728463899999485d+0*t31+1)-4.431373767749538
     7         d-3*t11*t29*t30-2.1552121112611455d-1*t22*t23*t24)+6.3144
     8         67057030003d-2*t14*t8*(-4.9237251092134837d-1*t27*t9/t14*
     9         *2-1.6412417030711612d-1*t15*t25)+4.431373767749538d-3*t1
     :         1*t29*t30+2.1552121112611455d-1*t22*t23*t24)*wght+2.0d+0*
     ;         t21*wght+Amat(iq,D1_RB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_vwn1_rpa functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{0.0448998886415768}\over{{\it t_4}+13.072}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 10.29581201158544\,{\it t_3}+{\it t_8}+42.7198\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.409286\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_3}+{\it t_8}+101.578\right| }}\\\\ 
C>   {\it t_{17}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}^{.3333333333333333}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\sqrt{{\it t_{18}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.575246635799487\,{\it t_{19}}\\\\ 
C>   {\it t_{21}} &=& \arctan \left({{0.0448998886415768}\over{{\it t_{20}}+13.072}}\right)\\\\ 
C>   {\it t_{22}} &=& .7876233178997433\,{\it t_{19}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{18}}}}\\\\ 
C>   {\it t_{24}} &=& .6203504908994001\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\left| 10.29581201158544\,{\it t_{19}}+{\it t_{24}}+42.7198\right| }}\\\\ 
C>   {\it t_{26}} &=& \log \left(1.0\,\left({\it t_{22}}+0.409286\right)^2\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{27}} &=& \left| {\it t_{23}}\right| \\\\ 
C>   {\it t_{28}} &=& \log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_{19}}+{\it t_{24}}+101.578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{16}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_6}+0.743294\right)^2\,{\it t_{16}}\right)-0.0310907\,{\it t_{12}}-.004431373767749538\,{\it t_{10}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_4}+20.1231}}\right)-20.52197293770518\,{\it t_5}\right)+0.0310907\,{\it t_{12}}+.004431373767749538\,{\it t_{10}}+20.52197293770518\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{17}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{13}}}\over{{\it t_{17}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{29}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_{22}}+0.743294\right)^2\,{\it t_{29}}\right)-0.0310907\,{\it t_{28}}-.004431373767749538\,{\it t_{26}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_{20}}+20.1231}}\right)-20.52197293770518\,{\it t_{21}}\right)+0.0310907\,{\it t_{28}}+.004431373767749538\,{\it t_{26}}+20.52197293770518\,{\it t_{21}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn1_rpa_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_c_vwn1_rpa functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {\it t_1}^{.3333333333333333}\\\\ 
C>   {\it t_3} &=& {{1}\over{\sqrt{{\it t_2}}}}\\\\ 
C>   {\it t_4} &=& 1.575246635799487\,{\it t_3}\\\\ 
C>   {\it t_5} &=& \arctan \left({{0.0448998886415768}\over{{\it t_4}+13.072}}\right)\\\\ 
C>   {\it t_6} &=& .7876233178997433\,{\it t_3}\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_8} &=& .6203504908994001\,{\it t_7}\\\\ 
C>   {\it t_9} &=& {{1}\over{\left| 10.29581201158544\,{\it t_3}+{\it t_8}+42.7198\right| }}\\\\ 
C>   {\it t_{10}} &=& \log \left(1.0\,\left({\it t_6}+0.409286\right)^2\,{\it t_9}\right)\\\\ 
C>   {\it t_{11}} &=& \left| {\it t_7}\right| \\\\ 
C>   {\it t_{12}} &=& \log \left(.6203504908994001\,{\it t_{11}}\,{\it t_9}\right)\\\\ 
C>   {\it t_{13}} &=& -\rho_\beta\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}+\rho_\alpha\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_3}+{\it t_8}+101.578\right| }}\\\\ 
C>   {\it t_{17}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{18}} &=& {\it t_{17}}^{.3333333333333333}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\sqrt{{\it t_{18}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.575246635799487\,{\it t_{19}}\\\\ 
C>   {\it t_{21}} &=& \arctan \left({{0.0448998886415768}\over{{\it t_{20}}+13.072}}\right)\\\\ 
C>   {\it t_{22}} &=& .7876233178997433\,{\it t_{19}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{18}}}}\\\\ 
C>   {\it t_{24}} &=& .6203504908994001\,{\it t_{23}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\left| 10.29581201158544\,{\it t_{19}}+{\it t_{24}}+42.7198\right| }}\\\\ 
C>   {\it t_{26}} &=& \log \left(1.0\,\left({\it t_{22}}+0.409286\right)^2\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{27}} &=& \left| {\it t_{23}}\right| \\\\ 
C>   {\it t_{28}} &=& \log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{25}}\right)\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\left| 15.84942278842832\,{\it t_{19}}+{\it t_{24}}+101.578\right| }}\\\\ 
C>   f &=& 1.0\,{\it t_1}\,\left(1.923661050931536\,\left(\left({\it t_{14}}\,{\it t_{15}}+1.0\right)^{{{4}\over{3}}}+\left(1.0-{\it t_{14}}\,{\it t_{15}}\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{11}}\,{\it t_{16}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_6}+0.743294\right)^2\,{\it t_{16}}\right)-0.0310907\,{\it t_{12}}-.004431373767749538\,{\it t_{10}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_4}+20.1231}}\right)-20.52197293770518\,{\it t_5}\right)+0.0310907\,{\it t_{12}}+.004431373767749538\,{\it t_{10}}+20.52197293770518\,{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,{\it t_{17}}\,\left(1.923661050931536\,\left(\left({{\rho_s+{\it t_{13}}}\over{{\it t_{17}}}}+1.0\right)^{{{4}\over{3}}}-2.0\right)\,\left(0.01554535\,\log \left(.6203504908994001\,{\it t_{27}}\,{\it t_{29}}\right)+.002667310007273315\,\log \left(1.0\,\left({\it t_{22}}+0.743294\right)^2\,{\it t_{29}}\right)-0.0310907\,{\it t_{28}}-.004431373767749538\,{\it t_{26}}+.6188180297906176\,\arctan \left({{1.171685277708971}\over{{\it t_{20}}+20.1231}}\right)-20.52197293770518\,{\it t_{21}}\right)+0.0310907\,{\it t_{28}}+.004431373767749538\,{\it t_{26}}+20.52197293770518\,{\it t_{21}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] SH Vosko, L Wilk, M Nusair, Can.J.Phys. 58, 1200 (1980)  , DOI: 
C> <a href="http://dx.doi.org/10.1139/p80-159 "> 
C> 10.1139/p80-159 </a> 
C> 
      subroutine nwxcm_c_vwn1_rpa_d3(param,tol_rho,ipol,nq,wght, 
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
