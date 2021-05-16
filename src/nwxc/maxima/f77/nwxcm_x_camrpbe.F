C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_camrpbe.F 
C> The nwxcm_x_camrpbe functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_camrpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& 1.0-{{1}\over{e^{{{0.00449276922095889\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_3} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& 0.7481426920246764\,{\it t_1}\,{\it t_2}
C>    +0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_8} &=& \sqrt{{\it t_7}}\\\\ 
C>   {\it t_9} &=& {\it t_5}^2\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_9}}}\\\\ 
C>   {\it t_{12}} &=& e^ {- {{14.13716694115413\,{\it t_{10}}\,{
C>    \it t_{11}}}\over{{\it t_7}}} }\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_5}}}\\\\ 
C>   {\it t_{14}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{15}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{17}} &=& 0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    +0.9305257363491\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{{\it t_{17}}}\\\\ 
C>   {\it t_{19}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{20}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{19}}}\over{{\it t_{17}}}} }\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{24}} &=& 0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    +0.9305257363491\,{\it t_{21}}\\\\ 
C>   {\it t_{25}} &=& \sqrt{{\it t_{24}}}\\\\ 
C>   {\it t_{26}} &=& \rho_s^2\\\\ 
C>   {\it t_{27}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{26}}}\over{{\it t_{24}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    -0.9305257363491\,{\it t_{14}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{16}}\,{
C>    \it t_{18}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\beta}\over{{
C>    \it t_{18}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{16}}\,{\it t_{18}}\,\left({\it t_{20}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{17}}\,\left({
C>    \it t_{20}}-1.0\right)}\over{{\it t_{19}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)+1.0\,\left(-0.7481426920246764\,{
C>    \it t_1}\,{\it t_2}-0.9305257363491\,{\it t_1}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_6}\,{
C>    \it t_8}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\alpha}\over{{
C>    \it t_8}}}\right)+0.2659615202676212\,{\it t_5}\,{\it t_6}
C>    \,{\it t_8}\,\left({\it t_{12}}-{{0.03536776513153214\,{
C>    \it t_9}\,{\it t_7}\,\left({\it t_{12}}-1.0\right)}\over{{
C>    \it t_{10}}}}-1.5\right)\right)+{\it t_3}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    -0.9305257363491\,{\it t_{21}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{23}}\,{
C>    \it t_{25}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_s}\over{{
C>    \it t_{25}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{23}}\,{\it t_{25}}\,\left({\it t_{27}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{24}}\,\left({
C>    \it t_{27}}-1.0\right)}\over{{\it t_{26}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3,4,5]. 
C> 
C> ### References ### 
C> 
C> [1] unpublished , DOI: 
C> <a href="https://doi.org/unpublished "> 
C> unpublished </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27340 2015-08-13 
C> 
C> [4] rewrap.py revision 27339 2015-08-13 
C> 
C> [5] call_subroutine.py revision 27339 2015-08-13 
C> 
      subroutine nwxcm_x_camrpbe(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,fnc,Amat,Cmat) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = 5.641895835477563d-1
            t26 = 1/t10**3
            t27 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t28 = 1/t9**2
            t29 = 1.413716694115413d+1*t12*t14*t27*t28
     +        -2.827433388230826d+1*rhoa*t14*t15
            t30 = 1/t1
            t31 = exp(-4.49276922095889d-3*gammaaa*t2
     +        -1.413716694115413d+1*t12*t14*t15)
            fnc(iq) = 2.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*(-3.5461536035682834d-1*t10*t6
     +        *t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t29*t9
     +        +7.073553026306428d-2*t11*t17*t9/rhoa**3+t16*t29
     +        -3.536776513153214d-2*t11*t13*t17*t27)
     +        +1.329807601338106d-1*t18*t20*t27*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t25*(3.759942411946509d+0*t19
     +        *t20-1.8799712059732546d+0*rhoa*t19*t26*t27))
     +        -1.7730768017841417d-1*t20*t21*t27*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7)+1.0d+0*t22*(
     +        -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23))*wght+Amat
     +        (iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1*t10
     +        *t6*t7*t8*(2.659615202676212d-1*t10*(
     +        -1.680616229806896d-3*t15*t30*t31+4.751830440958629d-2
     +        *rhoa**6.666666666666666d-1*t14*t28*t31
     +        -1.1887928018410271d-4*t11*t17*t3/rhoa*
     +        *3.3333333333333337d+0)*t7*t8+4.4697924746587997d-4*t18
     +        *t20*t24*t3*t7-2.24003435190034d-2*t19*t25*t26*t31/t23)
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7)
     +        -3.361232459613792d-3*t22*t3*t30)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = -param(1)
            t7 = param(2)
            t8 = param(3)
            t9 = 1/rhoa
            t10 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t11 = t10**5.0d-1
            t12 = t8**2
            t13 = rhoa**2
            t14 = 1/t13
            t15 = 1/t12
            t16 = 1/t10
            t17 = exp(-1.413716694115413d+1*t13*t15*t16)
            t18 = t17-1.0d+0
            t19 = -3.536776513153214d-2*t10*t12*t14*t18+t17-1.5d+0
            t20 = 1/t8
            t21 = 1/t11
            t22 = 2.659615202676212d-1*t11*t19*t8*t9
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t20
     +        *t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9+t6+1.0d+0
            t24 = rhob**1.3333333333333333d+0
            t25 = 1/rhob**2.6666666666666666d+0
            t26 = exp(-4.49276922095889d-3*gammabb*t25)
            t27 = 1.0d+0-t26
            t28 = -7.481426920246764d-1*t24*t27-9.305257363491d-1*t24
            t29 = 1/rhob
            t30 = 7.481426920246764d-1*t24*t27+9.305257363491d-1*t24
            t31 = t30**5.0d-1
            t32 = rhob**2
            t33 = 1/t32
            t34 = 1/t30
            t35 = exp(-1.413716694115413d+1*t15*t32*t34)
            t36 = t35-1.0d+0
            t37 = -3.536776513153214d-2*t12*t30*t33*t36+t35-1.5d+0
            t38 = 1/t31
            t39 = 2.659615202676212d-1*t29*t31*t37*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t20
     +        *t38)
            t40 = -3.5461536035682834d-1*t29*t31*t39*t7*t8+t6+1.0d+0
            t41 = rhoa**3.333333333333333d-1
            t42 = 1/rhoa**2.3333333333333334d+0
            t43 = 5.641895835477563d-1
            t44 = 1/t11**3
            t45 = -8.963286558970113d-3*gammaaa*t3*t42
     +        +9.975235893662352d-1*t4*t41+1.2407009817987999d+0*t41
            t46 = 1/t10**2
            t47 = 1.413716694115413d+1*t13*t15*t45*t46
     +        -2.827433388230826d+1*rhoa*t15*t16
            t48 = rhob**3.333333333333333d-1
            t49 = 1/rhob**2.3333333333333334d+0
            t50 = 1/t31**3
            t51 = -8.963286558970113d-3*gammabb*t26*t49
     +        +9.975235893662352d-1*t27*t48+1.2407009817987999d+0*t48
            t52 = 1/t30**2
            t53 = 1.413716694115413d+1*t15*t32*t51*t52
     +        -2.827433388230826d+1*rhob*t15*t34
            t54 = 1/t1
            t55 = exp(-4.49276922095889d-3*gammaaa*t2
     +        -1.413716694115413d+1*t13*t15*t16)
            t56 = 1/t24
            t57 = exp(-1.413716694115413d+1*t15*t32*t34
     +        -4.49276922095889d-3*gammabb*t25)
            fnc(iq) = (1.0d+0*t23*t5+1.0d+0*t28*t40)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*(-3.5461536035682834d-1*t11*t7
     +        *t8*t9*(2.659615202676212d-1*t11*(-3.536776513153214d-2
     +        *t10*t12*t14*t17*t47+t17*t47-3.536776513153214d-2*t12
     +        *t14*t18*t45+7.073553026306428d-2*t10*t12*t18/rhoa**3)
     +        *t8*t9+1.329807601338106d-1*t19*t21*t45*t8*t9
     +        -2.659615202676212d-1*t11*t14*t19*t8
     +        +3.5449077018110403d+0*t17*t43*(3.759942411946509d+0*t20
     +        *t21-1.8799712059732546d+0*rhoa*t20*t44*t45))
     +        -1.7730768017841417d-1*t21*t22*t45*t7*t8*t9
     +        +3.5461536035682834d-1*t11*t14*t22*t7*t8)+1.0d+0*t23*
     +        (8.963286558970113d-3*gammaaa*t3*t42
     +        -9.975235893662352d-1*t4*t41-1.2407009817987999d+0*t41))
     +        *wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t28*(-3.5461536035682834d-1*t29
     +        *t31*t7*t8*(2.659615202676212d-1*t29*t31*(
     +        -3.536776513153214d-2*t12*t30*t33*t35*t53+t35*t53
     +        -3.536776513153214d-2*t12*t33*t36*t51
     +        +7.073553026306428d-2*t12*t30*t36/rhob**3)*t8
     +        +1.329807601338106d-1*t29*t37*t38*t51*t8
     +        -2.659615202676212d-1*t31*t33*t37*t8
     +        +3.5449077018110403d+0*t35*t43*(3.759942411946509d+0*t20
     +        *t38-1.8799712059732546d+0*rhob*t20*t50*t51))
     +        -1.7730768017841417d-1*t29*t38*t39*t51*t7*t8
     +        +3.5461536035682834d-1*t31*t33*t39*t7*t8)+1.0d+0*t40*
     +        (8.963286558970113d-3*gammabb*t26*t49
     +        -9.975235893662352d-1*t27*t48-1.2407009817987999d+0*t48)
     +        )*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1*t11
     +        *t7*t8*t9*(2.659615202676212d-1*t11*(
     +        -1.680616229806896d-3*t16*t54*t55+4.751830440958629d-2
     +        *rhoa**6.666666666666666d-1*t15*t46*t55
     +        -1.1887928018410271d-4*t12*t18*t3/rhoa*
     +        *3.3333333333333337d+0)*t8*t9+4.4697924746587997d-4*t19
     +        *t21*t3*t42*t8-2.24003435190034d-2*t20*t43*t44*t55/t41)
     +        -5.959723299545067d-4*t21*t22*t3*t42*t7*t8)
     +        -3.361232459613792d-3*t23*t3*t54)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t28*(-3.5461536035682834d-1*t29
     +        *t31*t7*t8*(2.659615202676212d-1*t29*t31*(
     +        -1.680616229806896d-3*t34*t56*t57+4.751830440958629d-2
     +        *rhob**6.666666666666666d-1*t15*t52*t57
     +        -1.1887928018410271d-4*t12*t26*t36/rhob*
     +        *3.3333333333333337d+0)*t8+4.4697924746587997d-4*t26*t37
     +        *t38*t49*t8-2.24003435190034d-2*t20*t43*t50*t57/t48)
     +        -5.959723299545067d-4*t26*t38*t39*t49*t7*t8)
     +        -3.361232459613792d-3*t26*t40*t56)*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = 5.641895835477563d-1
            t26 = 1/t10**3
            t27 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t28 = 1/t9**2
            t29 = 1.413716694115413d+1*t12*t14*t27*t28
     +        -2.827433388230826d+1*rhoa*t14*t15
            t30 = 1/t1
            t31 = exp(-4.49276922095889d-3*gammaaa*t2
     +        -1.413716694115413d+1*t12*t14*t15)
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t5*(-3.5461536035682834d-1*t10*t6
     +        *t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t29*t9
     +        +7.073553026306428d-2*t11*t17*t9/rhoa**3+t16*t29
     +        -3.536776513153214d-2*t11*t13*t17*t27)
     +        +1.329807601338106d-1*t18*t20*t27*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t25*(3.759942411946509d+0*t19
     +        *t20-1.8799712059732546d+0*rhoa*t19*t26*t27))
     +        -1.7730768017841417d-1*t20*t21*t27*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7)*wght+1.0d+0
     +        *t22*(-9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23)*wght+Amat(iq,
     +        D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t5*(-3.5461536035682834d-1*t10*t6
     +        *t7*t8*(2.659615202676212d-1*t10*(-1.680616229806896d-3
     +        *t15*t30*t31+4.751830440958629d-2*rhoa*
     +        *6.666666666666666d-1*t14*t28*t31-1.1887928018410271d-4
     +        *t11*t17*t3/rhoa**3.3333333333333337d+0)*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t25*t26*t31/t23)
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7)*wght
     +        -3.361232459613792d-3*t22*t3*t30*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammabb*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhob
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhob**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhob**3.333333333333333d-1
            t24 = 1/rhob**2.3333333333333334d+0
            t25 = 5.641895835477563d-1
            t26 = 1/t10**3
            t27 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammabb*t24*t3+1.2407009817987999d+0*t23
            t28 = 1/t9**2
            t29 = 1.413716694115413d+1*t12*t14*t27*t28
     +        -2.827433388230826d+1*rhob*t14*t15
            t30 = 1/t1
            t31 = exp(-4.49276922095889d-3*gammabb*t2
     +        -1.413716694115413d+1*t12*t14*t15)
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t5*(-3.5461536035682834d-1*t10*t6
     +        *t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t29*t9
     +        +7.073553026306428d-2*t11*t17*t9/rhob**3+t16*t29
     +        -3.536776513153214d-2*t11*t13*t17*t27)
     +        +1.329807601338106d-1*t18*t20*t27*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t25*(3.759942411946509d+0*t19
     +        *t20-1.8799712059732546d+0*rhob*t19*t26*t27))
     +        -1.7730768017841417d-1*t20*t21*t27*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7)*wght+1.0d+0
     +        *t22*(-9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammabb*t24*t3-1.2407009817987999d+0*t23)*wght+Amat(iq,
     +        D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t5*(-3.5461536035682834d-1*t10*t6
     +        *t7*t8*(2.659615202676212d-1*t10*(-1.680616229806896d-3
     +        *t15*t30*t31+4.751830440958629d-2*rhob*
     +        *6.666666666666666d-1*t14*t28*t31-1.1887928018410271d-4
     +        *t11*t17*t3/rhob**3.3333333333333337d+0)*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t25*t26*t31/t23)
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7)*wght
     +        -3.361232459613792d-3*t22*t3*t30*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_x_camrpbe 
C> 
C> \brief Evaluate the nwxcm_x_camrpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& 1.0-{{1}\over{e^{{{0.00449276922095889\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_3} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& 0.7481426920246764\,{\it t_1}\,{\it t_2}
C>    +0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_8} &=& \sqrt{{\it t_7}}\\\\ 
C>   {\it t_9} &=& {\it t_5}^2\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_9}}}\\\\ 
C>   {\it t_{12}} &=& e^ {- {{14.13716694115413\,{\it t_{10}}\,{
C>    \it t_{11}}}\over{{\it t_7}}} }\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_5}}}\\\\ 
C>   {\it t_{14}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{15}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{17}} &=& 0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    +0.9305257363491\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{{\it t_{17}}}\\\\ 
C>   {\it t_{19}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{20}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{19}}}\over{{\it t_{17}}}} }\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{24}} &=& 0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    +0.9305257363491\,{\it t_{21}}\\\\ 
C>   {\it t_{25}} &=& \sqrt{{\it t_{24}}}\\\\ 
C>   {\it t_{26}} &=& \rho_s^2\\\\ 
C>   {\it t_{27}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{26}}}\over{{\it t_{24}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    -0.9305257363491\,{\it t_{14}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{16}}\,{
C>    \it t_{18}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\beta}\over{{
C>    \it t_{18}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{16}}\,{\it t_{18}}\,\left({\it t_{20}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{17}}\,\left({
C>    \it t_{20}}-1.0\right)}\over{{\it t_{19}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)+1.0\,\left(-0.7481426920246764\,{
C>    \it t_1}\,{\it t_2}-0.9305257363491\,{\it t_1}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_6}\,{
C>    \it t_8}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\alpha}\over{{
C>    \it t_8}}}\right)+0.2659615202676212\,{\it t_5}\,{\it t_6}
C>    \,{\it t_8}\,\left({\it t_{12}}-{{0.03536776513153214\,{
C>    \it t_9}\,{\it t_7}\,\left({\it t_{12}}-1.0\right)}\over{{
C>    \it t_{10}}}}-1.5\right)\right)+{\it t_3}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    -0.9305257363491\,{\it t_{21}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{23}}\,{
C>    \it t_{25}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_s}\over{{
C>    \it t_{25}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{23}}\,{\it t_{25}}\,\left({\it t_{27}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{24}}\,\left({
C>    \it t_{27}}-1.0\right)}\over{{\it t_{26}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3,4,5]. 
C> 
C> ### References ### 
C> 
C> [1] unpublished , DOI: 
C> <a href="https://doi.org/unpublished "> 
C> unpublished </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27340 2015-08-13 
C> 
C> [4] rewrap.py revision 27339 2015-08-13 
C> 
C> [5] call_subroutine.py revision 27339 2015-08-13 
C> 
      subroutine nwxcm_x_camrpbe_d2(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhoa*t19*t27*t28
            t30 = 1/rhoa**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhoa*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammaaa*t2)
            t40 = 1/rhoa**3.3333333333333337d+0
            t41 = rhoa**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammaaa**2
            t47 = 1/rhoa**6
            t48 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammaaa*t3*t40
            t49 = 1/t10**5
            t50 = t28**2
            t51 = 1/rhoa**4
            t52 = 1/t9**3
            t53 = -2.827433388230826d+1*t12*t14*t50*t52
     +        +1.413716694115413d+1*t12*t14*t31*t48
     +        +5.654866776461652d+1*rhoa*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t54 = t32**2
            t55 = 1/rhoa**5
            t56 = 4.026997797077544d-5*gammaaa*t3*t55
     +        -4.4816432794850564d-3*t24*t3
            t57 = 1.413716694115413d+1*t12*t14*t31*t56
     +        -9.503660881917259d-2*t14*t28*t3*t41*t52
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t58 = 1/rhoa**3.6666666666666664d+0
            t59 = exp(-8.98553844191778d-3*gammaaa*t2)
            t60 = exp(t38-8.98553844191778d-3*gammaaa*t2)
            t61 = 4.751830440958629d-2*t14*t3*t31*t41
     +        -4.49276922095889d-3*t2
            fnc(iq) = 2.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t35*t5+1.0d+0*t22*t25)*wght+Amat
     +        (iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t44*t5-3.361232459613792d-3*t22
     +        *t3*t36)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t54*t9
     +        -3.536776513153214d-2*t11*t13*t16*t53*t9
     +        -2.1220659078919285d-1*t11*t17*t51*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t54+t16
     +        *t53-3.536776513153214d-2*t11*t13*t17*t48
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30)
     +        -6.649038006690531d-2*t18*t27*t50*t7*t8
     +        +1.329807601338106d-1*t18*t20*t48*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*(2.8199568089598814d+0
     +        *rhoa*t19*t49*t50-1.8799712059732546d+0*rhoa*t19*t27*t48
     +        -3.759942411946509d+0*t19*t27*t28)+3.5449077018110403d+0
     +        *t16*t26*t29*t32)+8.865384008920708d-2*t21*t27*t50*t6*t7
     +        *t8-1.7730768017841417d-1*t20*t21*t48*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7)+1.0d+0*t22
     +        *(1.0738660792206783d-4*t3*t46*t47-3.3250786312207836d-1
     +        *t4*t45-4.135669939329333d-1*t45-8.963286558970114d-3
     +        *gammaaa*t3*t40)+2.0d+0*t25*t35)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t57*t9+t16*t57
     +        -3.536776513153214d-2*t11*t13*t17*t56
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39+2.3775856036820542d-4
     +        *t11*t17*t3/rhoa**4.333333333333333d+0)
     +        +1.329807601338106d-1*t18*t20*t56*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t56+9.478530360984885d-3*t19*t28*t3*t37
     +        *t49-6.31902024065659d-3*t19*t27*t3*t36)
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t56*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7)+1.0d+0
     +        *t22*(4.4816432794850564d-3*t24*t3-4.026997797077544d-5
     +        *gammaaa*t3*t55)+1.0d+0*t25*t44-3.361232459613792d-3*t3
     +        *t35*t36)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*
     +        (4.751830440958629d-2*t14*t31*t39*t41*t61
     +        -1.680616229806896d-3*t15*t36*t39*t61
     +        -3.194401344146213d-4*t14*t45*t52*t60
     +        +5.340971710208847d-7*t11*t17*t3*t47)*t7*t8
     +        -7.512005776780308d-7*t18*t27*t58*t59*t7
     +        -2.0081746054220725d-6*t18*t20*t3*t55*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t61
     +        +1.1293914261386051d-4*t19*t26*t49*t60/rhoa*
     +        *1.6666666666666669d+0)+1.0016007702373746d-6*t21*t27
     +        *t58*t59*t6*t7+2.677566140562764d-6*t20*t21*t3*t55*t6*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7)
     +        +1.510124173904079d-5*t22*t3*t51-6.722464919227584d-3*t3
     +        *t36*t44)*wght+Cmat2(iq,D2_GAA_GAA)
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = -param(1)
            t7 = param(2)
            t8 = param(3)
            t9 = 1/rhoa
            t10 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t11 = t10**5.0d-1
            t12 = t8**2
            t13 = rhoa**2
            t14 = 1/t13
            t15 = 1/t12
            t16 = 1/t10
            t17 = exp(-1.413716694115413d+1*t13*t15*t16)
            t18 = t17-1.0d+0
            t19 = -3.536776513153214d-2*t10*t12*t14*t18+t17-1.5d+0
            t20 = 1/t8
            t21 = 1/t11
            t22 = 2.659615202676212d-1*t11*t19*t8*t9
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t20
     +        *t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9+t6+1.0d+0
            t24 = rhob**1.3333333333333333d+0
            t25 = 1/rhob**2.6666666666666666d+0
            t26 = exp(-4.49276922095889d-3*gammabb*t25)
            t27 = 1.0d+0-t26
            t28 = -7.481426920246764d-1*t24*t27-9.305257363491d-1*t24
            t29 = 1/rhob
            t30 = 7.481426920246764d-1*t24*t27+9.305257363491d-1*t24
            t31 = t30**5.0d-1
            t32 = rhob**2
            t33 = 1/t32
            t34 = 1/t30
            t35 = exp(-1.413716694115413d+1*t15*t32*t34)
            t36 = t35-1.0d+0
            t37 = -3.536776513153214d-2*t12*t30*t33*t36+t35-1.5d+0
            t38 = 1/t31
            t39 = 2.659615202676212d-1*t29*t31*t37*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t20
     +        *t38)
            t40 = -3.5461536035682834d-1*t29*t31*t39*t7*t8+t6+1.0d+0
            t41 = rhoa**3.333333333333333d-1
            t42 = 1/rhoa**2.3333333333333334d+0
            t43 = 8.963286558970113d-3*gammaaa*t3*t42
     +        -9.975235893662352d-1*t4*t41-1.2407009817987999d+0*t41
            t44 = 5.641895835477563d-1
            t45 = 1/t11**3
            t46 = -8.963286558970113d-3*gammaaa*t3*t42
     +        +9.975235893662352d-1*t4*t41+1.2407009817987999d+0*t41
            t47 = 3.759942411946509d+0*t20*t21-1.8799712059732546d+0
     +        *rhoa*t20*t45*t46
            t48 = 1/rhoa**3
            t49 = 1/t10**2
            t50 = 1.413716694115413d+1*t13*t15*t46*t49
     +        -2.827433388230826d+1*rhoa*t15*t16
            t51 = -3.536776513153214d-2*t10*t12*t14*t17*t50+t17*t50
     +        +7.073553026306428d-2*t10*t12*t18*t48
     +        -3.536776513153214d-2*t12*t14*t18*t46
            t52 = 2.659615202676212d-1*t11*t51*t8*t9
     +        +1.329807601338106d-1*t19*t21*t46*t8*t9
     +        -2.659615202676212d-1*t11*t14*t19*t8
     +        +3.5449077018110403d+0*t17*t44*t47
            t53 = -3.5461536035682834d-1*t11*t52*t7*t8*t9
     +        -1.7730768017841417d-1*t21*t22*t46*t7*t8*t9
     +        +3.5461536035682834d-1*t11*t14*t22*t7*t8
            t54 = rhob**3.333333333333333d-1
            t55 = 1/rhob**2.3333333333333334d+0
            t56 = 8.963286558970113d-3*gammabb*t26*t55
     +        -9.975235893662352d-1*t27*t54-1.2407009817987999d+0*t54
            t57 = 1/t31**3
            t58 = -8.963286558970113d-3*gammabb*t26*t55
     +        +9.975235893662352d-1*t27*t54+1.2407009817987999d+0*t54
            t59 = 3.759942411946509d+0*t20*t38-1.8799712059732546d+0
     +        *rhob*t20*t57*t58
            t60 = 1/rhob**3
            t61 = 1/t30**2
            t62 = 1.413716694115413d+1*t15*t32*t58*t61
     +        -2.827433388230826d+1*rhob*t15*t34
            t63 = -3.536776513153214d-2*t12*t30*t33*t35*t62+t35*t62
     +        +7.073553026306428d-2*t12*t30*t36*t60
     +        -3.536776513153214d-2*t12*t33*t36*t58
            t64 = 2.659615202676212d-1*t29*t31*t63*t8
     +        +1.329807601338106d-1*t29*t37*t38*t58*t8
     +        -2.659615202676212d-1*t31*t33*t37*t8
     +        +3.5449077018110403d+0*t35*t44*t59
            t65 = -3.5461536035682834d-1*t29*t31*t64*t7*t8
     +        -1.7730768017841417d-1*t29*t38*t39*t58*t7*t8
     +        +3.5461536035682834d-1*t31*t33*t39*t7*t8
            t66 = 1/t1
            t67 = 1/t41
            t68 = -1.413716694115413d+1*t13*t15*t16
            t69 = exp(t68-4.49276922095889d-3*gammaaa*t2)
            t70 = 1/rhoa**3.3333333333333337d+0
            t71 = rhoa**6.666666666666666d-1
            t72 = 4.751830440958629d-2*t15*t49*t69*t71
     +        -1.1887928018410271d-4*t12*t18*t3*t70
     +        -1.680616229806896d-3*t16*t66*t69
            t73 = 2.659615202676212d-1*t11*t72*t8*t9
     +        +4.4697924746587997d-4*t19*t21*t3*t42*t8
     +        -2.24003435190034d-2*t20*t44*t45*t67*t69
            t74 = -3.5461536035682834d-1*t11*t7*t73*t8*t9
     +        -5.959723299545067d-4*t21*t22*t3*t42*t7*t8
            t75 = 1/t24
            t76 = 1/t54
            t77 = -1.413716694115413d+1*t15*t32*t34
            t78 = exp(t77-4.49276922095889d-3*gammabb*t25)
            t79 = 1/rhob**3.3333333333333337d+0
            t80 = rhob**6.666666666666666d-1
            t81 = 4.751830440958629d-2*t15*t61*t78*t80
     +        -1.1887928018410271d-4*t12*t26*t36*t79
     +        -1.680616229806896d-3*t34*t75*t78
            t82 = 2.659615202676212d-1*t29*t31*t8*t81
     +        +4.4697924746587997d-4*t26*t37*t38*t55*t8
     +        -2.24003435190034d-2*t20*t44*t57*t76*t78
            t83 = -3.5461536035682834d-1*t29*t31*t7*t8*t82
     +        -5.959723299545067d-4*t26*t38*t39*t55*t7*t8
            t84 = 1/t71
            t85 = gammaaa**2
            t86 = 1/rhoa**6
            t87 = -1.0738660792206783d-4*t3*t85*t86
     +        +3.3250786312207836d-1*t4*t84+4.135669939329333d-1*t84
     +        +8.963286558970114d-3*gammaaa*t3*t70
            t88 = 1/t11**5
            t89 = t46**2
            t90 = 1/rhoa**4
            t91 = 1/t10**3
            t92 = -2.827433388230826d+1*t13*t15*t89*t91
     +        +1.413716694115413d+1*t13*t15*t49*t87
     +        +5.654866776461652d+1*rhoa*t15*t46*t49
     +        -2.827433388230826d+1*t15*t16
            t93 = t50**2
            t94 = 1/t80
            t95 = gammabb**2
            t96 = 1/rhob**6
            t97 = -1.0738660792206783d-4*t26*t95*t96
     +        +3.3250786312207836d-1*t27*t94+4.135669939329333d-1*t94
     +        +8.963286558970114d-3*gammabb*t26*t79
            t98 = 1/t31**5
            t99 = t58**2
            t100 = 1/rhob**4
            t101 = 1/t30**3
            t102 = -2.827433388230826d+1*t101*t15*t32*t99
     +        +1.413716694115413d+1*t15*t32*t61*t97
     +        +5.654866776461652d+1*rhob*t15*t58*t61
     +        -2.827433388230826d+1*t15*t34
            t103 = t62**2
            t104 = 1/rhoa**5
            t105 = 4.026997797077544d-5*gammaaa*t104*t3
     +        -4.4816432794850564d-3*t3*t42
            t106 = -9.503660881917259d-2*t15*t3*t46*t71*t91
     +        +9.503660881917259d-2*t15*t3*t49*t67
     +        +1.413716694115413d+1*t105*t13*t15*t49
            t107 = 1/rhob**5
            t108 = 4.026997797077544d-5*gammabb*t107*t26
     +        -4.4816432794850564d-3*t26*t55
            t109 = -9.503660881917259d-2*t101*t15*t26*t58*t80
     +        +9.503660881917259d-2*t15*t26*t61*t76
     +        +1.413716694115413d+1*t108*t15*t32*t61
            t110 = 1/rhoa**3.6666666666666664d+0
            t111 = exp(-8.98553844191778d-3*gammaaa*t2)
            t112 = exp(t68-8.98553844191778d-3*gammaaa*t2)
            t113 = 4.751830440958629d-2*t15*t3*t49*t71
     +        -4.49276922095889d-3*t2
            t114 = 1/rhob**3.6666666666666664d+0
            t115 = exp(-8.98553844191778d-3*gammabb*t25)
            t116 = exp(t77-8.98553844191778d-3*gammabb*t25)
            t117 = 4.751830440958629d-2*t15*t26*t61*t80
     +        -4.49276922095889d-3*t25
            fnc(iq) = (1.0d+0*t23*t5+1.0d+0*t28*t40)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*t53+1.0d+0*t23*t43)*wght+Amat
     +        (iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t28*t65+1.0d+0*t40*t56)*wght+Amat
     +        (iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t5*t74-3.361232459613792d-3*t23
     +        *t3*t66)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t28*t83-3.361232459613792d-3*t26
     +        *t40*t75)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t11*t7*t8*t9*(2.659615202676212d-1*t11*t8*t9*(
     +        -3.536776513153214d-2*t10*t12*t14*t17*t93+t17*t93
     +        -3.536776513153214d-2*t10*t12*t14*t17*t92+t17*t92
     +        -2.1220659078919285d-1*t10*t12*t18*t90
     +        -3.536776513153214d-2*t12*t14*t18*t87
     +        +1.4147106052612857d-1*t10*t12*t17*t48*t50
     +        -7.073553026306428d-2*t12*t14*t17*t46*t50
     +        +1.4147106052612857d-1*t12*t18*t46*t48)
     +        -6.649038006690531d-2*t19*t45*t8*t89*t9
     +        +1.329807601338106d-1*t19*t21*t8*t87*t9
     +        +2.659615202676212d-1*t21*t46*t51*t8*t9
     +        +3.5449077018110403d+0*t17*t44*(2.8199568089598814d+0
     +        *rhoa*t20*t88*t89-1.8799712059732546d+0*rhoa*t20*t45*t87
     +        -3.759942411946509d+0*t20*t45*t46)-5.319230405352424d-1
     +        *t11*t14*t51*t8+5.319230405352424d-1*t11*t19*t48*t8
     +        -2.659615202676212d-1*t14*t19*t21*t46*t8
     +        +3.5449077018110403d+0*t17*t44*t47*t50)
     +        +8.865384008920708d-2*t22*t45*t7*t8*t89*t9
     +        -1.7730768017841417d-1*t21*t22*t7*t8*t87*t9
     +        -3.5461536035682834d-1*t21*t46*t52*t7*t8*t9
     +        +7.092307207136567d-1*t11*t14*t52*t7*t8
     +        -7.092307207136567d-1*t11*t22*t48*t7*t8
     +        +3.5461536035682834d-1*t14*t21*t22*t46*t7*t8)+1.0d+0*t23
     +        *(1.0738660792206783d-4*t3*t85*t86-3.3250786312207836d-1
     +        *t4*t84-4.135669939329333d-1*t84-8.963286558970114d-3
     +        *gammaaa*t3*t70)+2.0d+0*t43*t53)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (1.0d+0*t28*(-3.5461536035682834d-1
     +        *t29*t31*t7*t8*(3.5449077018110403d+0*t35*t44*
     +        (2.8199568089598814d+0*rhob*t20*t98*t99
     +        -1.8799712059732546d+0*rhob*t20*t57*t97
     +        -3.759942411946509d+0*t20*t57*t58)-6.649038006690531d-2
     +        *t29*t37*t57*t8*t99+2.659615202676212d-1*t29*t31*t8*(
     +        -3.536776513153214d-2*t12*t33*t36*t97
     +        +1.4147106052612857d-1*t12*t30*t35*t60*t62
     +        -7.073553026306428d-2*t12*t33*t35*t58*t62
     +        +1.4147106052612857d-1*t12*t36*t58*t60
     +        -2.1220659078919285d-1*t100*t12*t30*t36
     +        -3.536776513153214d-2*t103*t12*t30*t33*t35
     +        -3.536776513153214d-2*t102*t12*t30*t33*t35+t103*t35+t102
     +        *t35)+1.329807601338106d-1*t29*t37*t38*t8*t97
     +        +2.659615202676212d-1*t29*t38*t58*t63*t8
     +        -5.319230405352424d-1*t31*t33*t63*t8
     +        +5.319230405352424d-1*t31*t37*t60*t8
     +        -2.659615202676212d-1*t33*t37*t38*t58*t8
     +        +3.5449077018110403d+0*t35*t44*t59*t62)
     +        +8.865384008920708d-2*t29*t39*t57*t7*t8*t99
     +        -1.7730768017841417d-1*t29*t38*t39*t7*t8*t97
     +        -3.5461536035682834d-1*t29*t38*t58*t64*t7*t8
     +        +7.092307207136567d-1*t31*t33*t64*t7*t8
     +        -7.092307207136567d-1*t31*t39*t60*t7*t8
     +        +3.5461536035682834d-1*t33*t38*t39*t58*t7*t8)+1.0d+0*t40
     +        *(1.0738660792206783d-4*t26*t95*t96
     +        -3.3250786312207836d-1*t27*t94-4.135669939329333d-1*t94
     +        -8.963286558970114d-3*gammabb*t26*t79)+2.0d+0*t56*t65)
     +        *wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t11*t7*t8*t9*(1.329807601338106d-1*t21*t46*t72*t8*t9
     +        +2.659615202676212d-1*t11*(4.751830440958629d-2*t15*t49
     +        *t50*t69*t71-1.1887928018410271d-4*t12*t50*t69*t70
     +        -1.680616229806896d-3*t16*t50*t66*t69
     +        -1.680616229806896d-3*t46*t49*t66*t69
     +        +3.361232459613792d-3*t16*t42*t69+2.3775856036820542d-4
     +        *t12*t18*t3/rhoa**4.333333333333333d+0
     +        -3.536776513153214d-2*t105*t12*t14*t18
     +        -3.536776513153214d-2*t10*t106*t12*t14*t17+t106*t17)*t8
     +        *t9+1.329807601338106d-1*t105*t19*t21*t8*t9
     +        +3.5449077018110403d+0*t17*t44*(9.478530360984885d-3*t20
     +        *t3*t46*t67*t88-6.31902024065659d-3*t20*t3*t45*t66
     +        -1.8799712059732546d+0*rhoa*t105*t20*t45)
     +        -2.659615202676212d-1*t11*t14*t72*t8
     +        -4.4697924746587997d-4*t19*t21*t3*t70*t8
     +        +4.4697924746587997d-4*t21*t3*t42*t51*t8
     +        -2.2348962373293998d-4*t19*t3*t42*t45*t46*t8
     +        +1.6844800327854395d-1*t15*t44*t47*t49*t69*t71)
     +        -1.7730768017841417d-1*t21*t46*t7*t73*t8*t9
     +        -1.7730768017841417d-1*t105*t21*t22*t7*t8*t9
     +        +3.5461536035682834d-1*t11*t14*t7*t73*t8
     +        +5.959723299545067d-4*t21*t22*t3*t7*t70*t8
     +        -5.959723299545067d-4*t21*t3*t42*t52*t7*t8
     +        +2.9798616497725333d-4*t22*t3*t42*t45*t46*t7*t8)+1.0d+0
     +        *t43*t74-3.361232459613792d-3*t3*t53*t66+1.0d+0*t23*
     +        (4.4816432794850564d-3*t3*t42-4.026997797077544d-5
     +        *gammaaa*t104*t3))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (1.0d+0*t28*(-3.5461536035682834d-1
     +        *t29*t31*t7*t8*(3.5449077018110403d+0*t35*t44*
     +        (9.478530360984885d-3*t20*t26*t58*t76*t98
     +        -6.31902024065659d-3*t20*t26*t57*t75
     +        -1.8799712059732546d+0*rhob*t108*t20*t57)
     +        +1.329807601338106d-1*t29*t38*t58*t8*t81
     +        -2.659615202676212d-1*t31*t33*t8*t81
     +        +2.659615202676212d-1*t29*t31*t8*(4.751830440958629d-2
     +        *t15*t61*t62*t78*t80-1.1887928018410271d-4*t12*t62*t78
     +        *t79-1.680616229806896d-3*t34*t62*t75*t78
     +        -1.680616229806896d-3*t58*t61*t75*t78
     +        +3.361232459613792d-3*t34*t55*t78-3.536776513153214d-2
     +        *t108*t12*t33*t36+2.3775856036820542d-4*t12*t26*t36/rhob
     +        **4.333333333333333d+0-3.536776513153214d-2*t109*t12*t30
     +        *t33*t35+t109*t35)+1.6844800327854395d-1*t15*t44*t59*t61
     +        *t78*t80-4.4697924746587997d-4*t26*t37*t38*t79*t8
     +        +4.4697924746587997d-4*t26*t38*t55*t63*t8
     +        -2.2348962373293998d-4*t26*t37*t55*t57*t58*t8
     +        +1.329807601338106d-1*t108*t29*t37*t38*t8)
     +        -1.7730768017841417d-1*t29*t38*t58*t7*t8*t82
     +        +3.5461536035682834d-1*t31*t33*t7*t8*t82
     +        +5.959723299545067d-4*t26*t38*t39*t7*t79*t8
     +        -5.959723299545067d-4*t26*t38*t55*t64*t7*t8
     +        +2.9798616497725333d-4*t26*t39*t55*t57*t58*t7*t8
     +        -1.7730768017841417d-1*t108*t29*t38*t39*t7*t8)+1.0d+0
     +        *t56*t83-3.361232459613792d-3*t26*t65*t75+1.0d+0*t40*
     +        (4.4816432794850564d-3*t26*t55-4.026997797077544d-5
     +        *gammabb*t107*t26))*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t11*t7*t8*t9*(2.659615202676212d-1*t11*t8*t9*(
     +        -3.194401344146213d-4*t112*t15*t84*t91
     +        +5.340971710208847d-7*t12*t18*t3*t86
     +        +4.751830440958629d-2*t113*t15*t49*t69*t71
     +        -1.680616229806896d-3*t113*t16*t66*t69)
     +        +1.1293914261386051d-4*t112*t20*t44*t88/rhoa*
     +        *1.6666666666666669d+0+8.939584949317599d-4*t21*t3*t42
     +        *t72*t8-7.512005776780308d-7*t110*t111*t19*t45*t8
     +        -2.0081746054220725d-6*t104*t19*t21*t3*t8
     +        -2.24003435190034d-2*t113*t20*t44*t45*t67*t69)
     +        -1.1919446599090135d-3*t21*t3*t42*t7*t73*t8
     +        +1.0016007702373746d-6*t110*t111*t22*t45*t7*t8
     +        +2.677566140562764d-6*t104*t21*t22*t3*t7*t8)
     +        +1.510124173904079d-5*t23*t3*t90-6.722464919227584d-3*t3
     +        *t66*t74)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (1.0d+0*t28*(-3.5461536035682834d-1
     +        *t29*t31*t7*t8*(1.1293914261386051d-4*t116*t20*t44*t98
     +        /rhob**1.6666666666666669d+0+2.659615202676212d-1*t29
     +        *t31*t8*(5.340971710208847d-7*t12*t26*t36*t96
     +        -3.194401344146213d-4*t101*t116*t15*t94
     +        +4.751830440958629d-2*t117*t15*t61*t78*t80
     +        -1.680616229806896d-3*t117*t34*t75*t78)
     +        +8.939584949317599d-4*t26*t38*t55*t8*t81
     +        -7.512005776780308d-7*t114*t115*t37*t57*t8
     +        -2.0081746054220725d-6*t107*t26*t37*t38*t8
     +        -2.24003435190034d-2*t117*t20*t44*t57*t76*t78)
     +        -1.1919446599090135d-3*t26*t38*t55*t7*t8*t82
     +        +1.0016007702373746d-6*t114*t115*t39*t57*t7*t8
     +        +2.677566140562764d-6*t107*t26*t38*t39*t7*t8)
     +        -6.722464919227584d-3*t26*t75*t83+1.510124173904079d-5
     +        *t100*t26*t40)*wght+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhoa*t19*t27*t28
            t30 = 1/rhoa**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhoa*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammaaa*t2)
            t40 = 1/rhoa**3.3333333333333337d+0
            t41 = rhoa**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammaaa**2
            t47 = 1/rhoa**6
            t48 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammaaa*t3*t40
            t49 = 1/t10**5
            t50 = t28**2
            t51 = 1/rhoa**4
            t52 = 1/t9**3
            t53 = -2.827433388230826d+1*t12*t14*t50*t52
     +        +1.413716694115413d+1*t12*t14*t31*t48
     +        +5.654866776461652d+1*rhoa*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t54 = t32**2
            t55 = 1/rhoa**5
            t56 = 4.026997797077544d-5*gammaaa*t3*t55
     +        -4.4816432794850564d-3*t24*t3
            t57 = 1.413716694115413d+1*t12*t14*t31*t56
     +        -9.503660881917259d-2*t14*t28*t3*t41*t52
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t58 = 1/rhoa**3.6666666666666664d+0
            t59 = exp(-8.98553844191778d-3*gammaaa*t2)
            t60 = exp(t38-8.98553844191778d-3*gammaaa*t2)
            t61 = 4.751830440958629d-2*t14*t3*t31*t41
     +        -4.49276922095889d-3*t2
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t35*t5*wght+1.0d+0*t22*t25*wght
     +        +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t44*t5*wght-3.361232459613792d-3
     +        *t22*t3*t36*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*t5*(-3.5461536035682834d-1*t10
     +        *t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t54*t9
     +        -3.536776513153214d-2*t11*t13*t16*t53*t9
     +        -2.1220659078919285d-1*t11*t17*t51*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t54+t16
     +        *t53-3.536776513153214d-2*t11*t13*t17*t48
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30)
     +        -6.649038006690531d-2*t18*t27*t50*t7*t8
     +        +1.329807601338106d-1*t18*t20*t48*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*(2.8199568089598814d+0
     +        *rhoa*t19*t49*t50-1.8799712059732546d+0*rhoa*t19*t27*t48
     +        -3.759942411946509d+0*t19*t27*t28)+3.5449077018110403d+0
     +        *t16*t26*t29*t32)+8.865384008920708d-2*t21*t27*t50*t6*t7
     +        *t8-1.7730768017841417d-1*t20*t21*t48*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7)*wght
     +        +1.0d+0*t22*(1.0738660792206783d-4*t3*t46*t47
     +        -3.3250786312207836d-1*t4*t45-4.135669939329333d-1*t45
     +        -8.963286558970114d-3*gammaaa*t3*t40)*wght+2.0d+0*t25
     +        *t35*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t57*t9+t16*t57
     +        -3.536776513153214d-2*t11*t13*t17*t56
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39+2.3775856036820542d-4
     +        *t11*t17*t3/rhoa**4.333333333333333d+0)
     +        +1.329807601338106d-1*t18*t20*t56*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t56+9.478530360984885d-3*t19*t28*t3*t37
     +        *t49-6.31902024065659d-3*t19*t27*t3*t36)
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t56*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7)*wght
     +        +1.0d+0*t22*(4.4816432794850564d-3*t24*t3
     +        -4.026997797077544d-5*gammaaa*t3*t55)*wght+1.0d+0*t25
     +        *t44*wght-3.361232459613792d-3*t3*t35*t36*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*
     +        (4.751830440958629d-2*t14*t31*t39*t41*t61
     +        -1.680616229806896d-3*t15*t36*t39*t61
     +        -3.194401344146213d-4*t14*t45*t52*t60
     +        +5.340971710208847d-7*t11*t17*t3*t47)*t7*t8
     +        -7.512005776780308d-7*t18*t27*t58*t59*t7
     +        -2.0081746054220725d-6*t18*t20*t3*t55*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t61
     +        +1.1293914261386051d-4*t19*t26*t49*t60/rhoa*
     +        *1.6666666666666669d+0)+1.0016007702373746d-6*t21*t27
     +        *t58*t59*t6*t7+2.677566140562764d-6*t20*t21*t3*t55*t6*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7)*wght
     +        +1.510124173904079d-5*t22*t3*t51*wght
     +        -6.722464919227584d-3*t3*t36*t44*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammabb*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhob
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhob**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhob**3.333333333333333d-1
            t24 = 1/rhob**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammabb*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammabb*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhob*t19*t27*t28
            t30 = 1/rhob**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhob*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammabb*t2)
            t40 = 1/rhob**3.3333333333333337d+0
            t41 = rhob**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammabb**2
            t47 = 1/rhob**6
            t48 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammabb*t3*t40
            t49 = 1/t10**5
            t50 = t28**2
            t51 = 1/rhob**4
            t52 = 1/t9**3
            t53 = -2.827433388230826d+1*t12*t14*t50*t52
     +        +1.413716694115413d+1*t12*t14*t31*t48
     +        +5.654866776461652d+1*rhob*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t54 = t32**2
            t55 = 1/rhob**5
            t56 = 4.026997797077544d-5*gammabb*t3*t55
     +        -4.4816432794850564d-3*t24*t3
            t57 = 1.413716694115413d+1*t12*t14*t31*t56
     +        -9.503660881917259d-2*t14*t28*t3*t41*t52
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t58 = 1/rhob**3.6666666666666664d+0
            t59 = exp(-8.98553844191778d-3*gammabb*t2)
            t60 = exp(t38-8.98553844191778d-3*gammabb*t2)
            t61 = 4.751830440958629d-2*t14*t3*t31*t41
     +        -4.49276922095889d-3*t2
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t35*t5*wght+1.0d+0*t22*t25*wght
     +        +Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t44*t5*wght-3.361232459613792d-3
     +        *t22*t3*t36*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*t5*(-3.5461536035682834d-1*t10
     +        *t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t54*t9
     +        -3.536776513153214d-2*t11*t13*t16*t53*t9
     +        -2.1220659078919285d-1*t11*t17*t51*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t54+t16
     +        *t53-3.536776513153214d-2*t11*t13*t17*t48
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30)
     +        -6.649038006690531d-2*t18*t27*t50*t7*t8
     +        +1.329807601338106d-1*t18*t20*t48*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*(2.8199568089598814d+0
     +        *rhob*t19*t49*t50-1.8799712059732546d+0*rhob*t19*t27*t48
     +        -3.759942411946509d+0*t19*t27*t28)+3.5449077018110403d+0
     +        *t16*t26*t29*t32)+8.865384008920708d-2*t21*t27*t50*t6*t7
     +        *t8-1.7730768017841417d-1*t20*t21*t48*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7)*wght
     +        +1.0d+0*t22*(1.0738660792206783d-4*t3*t46*t47
     +        -3.3250786312207836d-1*t4*t45-4.135669939329333d-1*t45
     +        -8.963286558970114d-3*gammabb*t3*t40)*wght+2.0d+0*t25
     +        *t35*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t57*t9+t16*t57
     +        -3.536776513153214d-2*t11*t13*t17*t56
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39+2.3775856036820542d-4
     +        *t11*t17*t3/rhob**4.333333333333333d+0)
     +        +1.329807601338106d-1*t18*t20*t56*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhob*t19*t27*t56+9.478530360984885d-3*t19*t28*t3*t37
     +        *t49-6.31902024065659d-3*t19*t27*t3*t36)
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t56*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7)*wght
     +        +1.0d+0*t22*(4.4816432794850564d-3*t24*t3
     +        -4.026997797077544d-5*gammabb*t3*t55)*wght+1.0d+0*t25
     +        *t44*wght-3.361232459613792d-3*t3*t35*t36*wght+Cmat2(iq,
     +        D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*
     +        (4.751830440958629d-2*t14*t31*t39*t41*t61
     +        -1.680616229806896d-3*t15*t36*t39*t61
     +        -3.194401344146213d-4*t14*t45*t52*t60
     +        +5.340971710208847d-7*t11*t17*t3*t47)*t7*t8
     +        -7.512005776780308d-7*t18*t27*t58*t59*t7
     +        -2.0081746054220725d-6*t18*t20*t3*t55*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t61
     +        +1.1293914261386051d-4*t19*t26*t49*t60/rhob*
     +        *1.6666666666666669d+0)+1.0016007702373746d-6*t21*t27
     +        *t58*t59*t6*t7+2.677566140562764d-6*t20*t21*t3*t55*t6*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7)*wght
     +        +1.510124173904079d-5*t22*t3*t51*wght
     +        -6.722464919227584d-3*t3*t36*t44*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_x_camrpbe_d2 
C> 
C> \brief Evaluate the nwxcm_x_camrpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& 1.0-{{1}\over{e^{{{0.00449276922095889\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_3} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_7} &=& 0.7481426920246764\,{\it t_1}\,{\it t_2}
C>    +0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_8} &=& \sqrt{{\it t_7}}\\\\ 
C>   {\it t_9} &=& {\it t_5}^2\\\\ 
C>   {\it t_{10}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_9}}}\\\\ 
C>   {\it t_{12}} &=& e^ {- {{14.13716694115413\,{\it t_{10}}\,{
C>    \it t_{11}}}\over{{\it t_7}}} }\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_5}}}\\\\ 
C>   {\it t_{14}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{15}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{17}} &=& 0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    +0.9305257363491\,{\it t_{14}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{{\it t_{17}}}\\\\ 
C>   {\it t_{19}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{20}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{19}}}\over{{\it t_{17}}}} }\\\\ 
C>   {\it t_{21}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& 1.0-{{1}\over{e^{{{0.00449276922095889
C>    \,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}}}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{24}} &=& 0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    +0.9305257363491\,{\it t_{21}}\\\\ 
C>   {\it t_{25}} &=& \sqrt{{\it t_{24}}}\\\\ 
C>   {\it t_{26}} &=& \rho_s^2\\\\ 
C>   {\it t_{27}} &=& e^ {- {{14.13716694115413\,{\it t_{11}}\,{
C>    \it t_{26}}}\over{{\it t_{24}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.7481426920246764\,{\it t_{14}}\,{\it t_{15}}
C>    -0.9305257363491\,{\it t_{14}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{16}}\,{
C>    \it t_{18}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\beta}\over{{
C>    \it t_{18}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{16}}\,{\it t_{18}}\,\left({\it t_{20}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{17}}\,\left({
C>    \it t_{20}}-1.0\right)}\over{{\it t_{19}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)+1.0\,\left(-0.7481426920246764\,{
C>    \it t_1}\,{\it t_2}-0.9305257363491\,{\it t_1}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_6}\,{
C>    \it t_8}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_\alpha}\over{{
C>    \it t_8}}}\right)+0.2659615202676212\,{\it t_5}\,{\it t_6}
C>    \,{\it t_8}\,\left({\it t_{12}}-{{0.03536776513153214\,{
C>    \it t_9}\,{\it t_7}\,\left({\it t_{12}}-1.0\right)}\over{{
C>    \it t_{10}}}}-1.5\right)\right)+{\it t_3}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.7481426920246764\,{\it t_{21}}\,{\it t_{22}}
C>    -0.9305257363491\,{\it t_{21}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_4}\,{\it t_5}\,{\it t_{23}}\,{
C>    \it t_{25}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{13}}\,\rho_s}\over{{
C>    \it t_{25}}}}\right)+0.2659615202676212\,{\it t_5}\,{
C>    \it t_{23}}\,{\it t_{25}}\,\left({\it t_{27}}
C>    -{{0.03536776513153214\,{\it t_9}\,{\it t_{24}}\,\left({
C>    \it t_{27}}-1.0\right)}\over{{\it t_{26}}}}-1.5\right)\right)
C>    +{\it t_3}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3,4,5]. 
C> 
C> ### References ### 
C> 
C> [1] unpublished , DOI: 
C> <a href="https://doi.org/unpublished "> 
C> unpublished </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27340 2015-08-13 
C> 
C> [4] rewrap.py revision 27339 2015-08-13 
C> 
C> [5] call_subroutine.py revision 27339 2015-08-13 
C> 
      subroutine nwxcm_x_camrpbe_d3(param,tol_rho,ipol,nq,wght, rho,
     +   rgamma,fnc,Amat,Amat2,Amat3, Cmat,Cmat2,Cmat3) 
c $Id: $ 
#ifdef NWXC_QUAD_PREC 
      implicit real(kind=selected_real_kind(30))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(30) 
#else 
      implicit real(kind=selected_real_kind(15))(a-h,o-z),integer(i-n) 
      integer, parameter :: rk=selected_real_kind(15) 
#endif 
#include "nwxc_param.fh" 
      double precision sr(NCOL_RHO)
      double precision sg(NCOL_GAMMA)
      double precision st(NCOL_TAU)
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the 
                                    !< density 
      integer ipol                  !< [Input] The number of spin 
                                    !< channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the 
                                    !< functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the 
                                             !< density gradients 
      double precision fnc(nq
     +   )      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative 
                                            !< wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative 
                                            !< wrt rgamma 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd 
                                             !< derivative wrt rgamma 
                                             !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd 
                                             !< derivative wrt rgamma 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhoa*t19*t27*t28
            t30 = 1/rhoa**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhoa*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammaaa*t2)
            t40 = 1/rhoa**3.3333333333333337d+0
            t41 = rhoa**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammaaa**2
            t47 = 1/rhoa**6
            t48 = 1.0738660792206783d-4*t3*t46*t47
     +        -3.3250786312207836d-1*t4*t45-4.135669939329333d-1*t45
     +        -8.963286558970114d-3*gammaaa*t3*t40
            t49 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammaaa*t3*t40
            t50 = 1/t10**5
            t51 = t28**2
            t52 = 2.8199568089598814d+0*rhoa*t19*t50*t51
     +        -1.8799712059732546d+0*rhoa*t19*t27*t49
     +        -3.759942411946509d+0*t19*t27*t28
            t53 = 1/rhoa**4
            t54 = 1/t9**3
            t55 = -2.827433388230826d+1*t12*t14*t51*t54
     +        +1.413716694115413d+1*t12*t14*t31*t49
     +        +5.654866776461652d+1*rhoa*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t56 = t32**2
            t57 = -3.536776513153214d-2*t11*t13*t16*t56*t9
     +        -3.536776513153214d-2*t11*t13*t16*t55*t9
     +        -2.1220659078919285d-1*t11*t17*t53*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t56+t16
     +        *t55-3.536776513153214d-2*t11*t13*t17*t49
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30
            t58 = 2.659615202676212d-1*t10*t57*t7*t8
     +        -6.649038006690531d-2*t18*t27*t51*t7*t8
     +        +1.329807601338106d-1*t18*t20*t49*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*t52+3.5449077018110403d+0
     +        *t16*t26*t29*t32
            t59 = -3.5461536035682834d-1*t10*t58*t6*t7*t8
     +        +8.865384008920708d-2*t21*t27*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t49*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7
            t60 = 1/rhoa**5
            t61 = 4.4816432794850564d-3*t24*t3-4.026997797077544d-5
     +        *gammaaa*t3*t60
            t62 = 4.026997797077544d-5*gammaaa*t3*t60
     +        -4.4816432794850564d-3*t24*t3
            t63 = -1.8799712059732546d+0*rhoa*t19*t27*t62
     +        +9.478530360984885d-3*t19*t28*t3*t37*t50
     +        -6.31902024065659d-3*t19*t27*t3*t36
            t64 = 1/rhoa**4.333333333333333d+0
            t65 = 1.413716694115413d+1*t12*t14*t31*t62
     +        -9.503660881917259d-2*t14*t28*t3*t41*t54
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t66 = -3.536776513153214d-2*t11*t13*t16*t65*t9+t16*t65
     +        +2.3775856036820542d-4*t11*t17*t3*t64
     +        -3.536776513153214d-2*t11*t13*t17*t62
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39
            t67 = 2.659615202676212d-1*t10*t66*t7*t8
     +        +1.329807601338106d-1*t18*t20*t62*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*t63+1.6844800327854395d-1
     +        *t14*t26*t29*t31*t39*t41
            t68 = -3.5461536035682834d-1*t10*t6*t67*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t6*t62*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7
            t69 = 1/rhoa**3.6666666666666664d+0
            t70 = exp(-8.98553844191778d-3*gammaaa*t2)
            t71 = 1/rhoa**1.6666666666666669d+0
            t72 = exp(t38-8.98553844191778d-3*gammaaa*t2)
            t73 = 4.751830440958629d-2*t14*t3*t31*t41
            t74 = t73-4.49276922095889d-3*t2
            t75 = 4.751830440958629d-2*t14*t31*t39*t41*t74
     +        -1.680616229806896d-3*t15*t36*t39*t74
     +        -3.194401344146213d-4*t14*t45*t54*t72
     +        +5.340971710208847d-7*t11*t17*t3*t47
            t76 = 2.659615202676212d-1*t10*t7*t75*t8
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t74
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72
     +        -7.512005776780308d-7*t18*t27*t69*t7*t70
     +        -2.0081746054220725d-6*t18*t20*t3*t60*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
            t77 = -3.5461536035682834d-1*t10*t6*t7*t76*t8
     +        +1.0016007702373746d-6*t21*t27*t6*t69*t7*t70
     +        +2.677566140562764d-6*t20*t21*t3*t6*t60*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7
            t78 = gammaaa**3
            t79 = 1/rhoa**9.666666666666666d+0
            t80 = 1/rhoa**7
            t81 = 7.517062554544749d-4*t3*t46*t80-1.286568658174524d-6
     +        *t3*t78*t79-2.2167190874805223d-1*t4*t71
     +        -2.757113292886222d-1*t71-3.386130477833154d-2*gammaaa
     +        *t3*t64
            t82 = 1/t10**7
            t83 = t28**3
            t84 = 1/t9**4
            t85 = 8.482300164692478d+1*t12*t14*t83*t84
     +        +1.413716694115413d+1*t12*t14*t31*t81
     +        -1.6964600329384957d+2*rhoa*t14*t51*t54
     +        -8.482300164692478d+1*t12*t14*t28*t49*t54
     +        +8.482300164692478d+1*rhoa*t14*t31*t49
     +        +8.482300164692478d+1*t14*t28*t31
            t86 = t32**3
            t87 = 1/rhoa**8.666666666666666d+0
            t88 = 4.8246324681544644d-7*t3*t46*t87
     +        -2.550431938149111d-4*gammaaa*t3*t47
     +        +1.0457167652131799d-2*t3*t40
            t89 = 1/rhoa**5.333333333333333d+0
            t90 = 1.413716694115413d+1*t12*t14*t31*t88
     +        +2.8510982645751776d-1*t14*t3*t41*t51*t84
     +        -5.654866776461652d+1*t12*t14*t28*t54*t62
     +        +5.654866776461652d+1*rhoa*t14*t31*t62
     +        -9.503660881917259d-2*t14*t3*t41*t49*t54
     +        -3.8014643527669034d-1*t14*t28*t3*t37*t54
     +        +9.503660881917259d-2*t14*t3*t31*t36
            t91 = 1/rhoa**7.666666666666667d+0
            t92 = 6.040496695616315d-5*t3*t60-1.8092371755579242d-7
     +        *gammaaa*t3*t91
            t93 = 1/rhoa**4.666666666666667d+0
            t94 = 1.413716694115413d+1*t12*t14*t31*t92
     +        +9.583204032438639d-4*t14*t28*t45*t70*t84
     +        -6.388802688292426d-4*t14*t54*t70*t71
     +        -1.9007321763834517d-1*t14*t3*t41*t54*t62
     +        +4.269775509670888d-4*t13*t14*t28*t3*t54
     +        -4.269775509670888d-4*t14*t3*t30*t31
            t95 = exp(-1.3478307662876668d-2*gammaaa*t2)
            t96 = 1/rhoa**6.333333333333333d+0
            t97 = exp(t38-1.3478307662876668d-2*gammaaa*t2)
            t98 = t73-8.98553844191778d-3*t2
            t99 = -3.194401344146213d-4*t14*t45*t54*t70
     +        -2.134887754835444d-4*t13*t14*t3*t31
            t100 = t74**2
            fnc(iq) = 2.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t35*t5+1.0d+0*t22*t25)*wght+Amat
     +        (iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t44*t5-3.361232459613792d-3*t22
     +        *t3*t36)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (1.0d+0*t5*t59+1.0d+0*t22*t48+2.0d+0
     +        *t25*t35)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (1.0d+0*t5*t68+1.0d+0*t22*t61+1.0d+0
     +        *t25*t44-3.361232459613792d-3*t3*t35*t36)*wght+Cmat2(iq,
     +        D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t5*t77+1.510124173904079d-5
     +        *t22*t3*t53-6.722464919227584d-3*t3*t36*t44)*wght+Cmat2
     +        (iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t86*t9
     +        -3.536776513153214d-2*t11*t13*t16*t85*t9
     +        +8.488263631567714d-1*t11*t17*t60*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t56*t9
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t55*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t55*t9
     +        -6.366197723675786d-1*t11*t16*t32*t53*t9+t16*t86+t16*t85
     +        -3.536776513153214d-2*t11*t13*t17*t81
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t56+3*t16*t32*t55
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t55
     +        -6.366197723675786d-1*t11*t17*t28*t53
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t49
     +        +2.1220659078919285d-1*t11*t17*t30*t49
     +        +4.244131815783857d-1*t11*t16*t28*t30*t32)
     +        +3.5449077018110403d+0*t16*t26*(-7.049892022399704d+0
     +        *rhoa*t19*t82*t83-1.8799712059732546d+0*rhoa*t19*t27*t81
     +        +8.459870426879645d+0*t19*t50*t51+8.459870426879645d+0
     +        *rhoa*t19*t28*t49*t50-5.639913617919763d+0*t19*t27*t49)
     +        +9.973557010035795d-2*t18*t50*t7*t8*t83
     +        +1.329807601338106d-1*t18*t20*t7*t8*t81
     +        +3.989422804014318d-1*t20*t28*t57*t7*t8
     +        -1.994711402007159d-1*t27*t33*t51*t7*t8
     +        +3.989422804014318d-1*t20*t33*t49*t7*t8
     +        -1.994711402007159d-1*t18*t27*t28*t49*t7*t8
     +        -7.978845608028636d-1*t10*t13*t57*t7
     +        -1.5957691216057274d+0*t10*t18*t53*t7
     +        +1.994711402007159d-1*t13*t18*t27*t51*t7
     +        -3.989422804014318d-1*t13*t18*t20*t49*t7
     +        +1.5957691216057274d+0*t10*t30*t33*t7
     +        -7.978845608028636d-1*t13*t20*t28*t33*t7
     +        +7.978845608028636d-1*t18*t20*t28*t30*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t56
     +        +3.5449077018110403d+0*t16*t26*t29*t55
     +        +7.089815403622081d+0*t16*t26*t32*t52)
     +        -1.3298076013381063d-1*t21*t50*t6*t7*t8*t83
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t81
     +        -5.319230405352425d-1*t20*t28*t58*t6*t7*t8
     +        +2.6596152026762127d-1*t27*t34*t51*t6*t7*t8
     +        -5.319230405352425d-1*t20*t34*t49*t6*t7*t8
     +        +2.6596152026762127d-1*t21*t27*t28*t49*t6*t7*t8
     +        +1.063846081070485d+0*t10*t13*t58*t6*t7
     +        +2.12769216214097d+0*t10*t21*t53*t6*t7
     +        -2.6596152026762127d-1*t13*t21*t27*t51*t6*t7
     +        +5.319230405352425d-1*t13*t20*t21*t49*t6*t7
     +        -2.12769216214097d+0*t10*t30*t34*t6*t7
     +        +1.063846081070485d+0*t13*t20*t28*t34*t6*t7
     +        -1.063846081070485d+0*t20*t21*t28*t30*t6*t7)+1.0d+0*t22*
     +        (-7.517062554544749d-4*t3*t46*t80+1.286568658174524d-6
     +        *t3*t78*t79+2.2167190874805223d-1*t4*t71
     +        +2.757113292886222d-1*t71+3.386130477833154d-2*gammaaa
     +        *t3*t64)+3.0d+0*t25*t59+3.0d+0*t35*t48)*wght+Amat3(iq,
     +        D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(-3.536776513153214d-2
     +        *t11*t13*t16*t9*t90+t16*t90-7.073553026306428d-2*t11*t13
     +        *t16*t32*t65*t9+1.4147106052612857d-1*t11*t16*t30*t65*t9
     +        -7.132756811046163d-4*t11*t17*t3*t89
     +        -3.536776513153214d-2*t11*t13*t17*t88+2*t16*t32*t65
     +        -7.073553026306428d-2*t11*t13*t16*t28*t65
     +        +4.7551712073641084d-4*t11*t32*t39*t64
     +        -7.073553026306428d-2*t11*t13*t16*t32*t62
     +        +1.4147106052612857d-1*t11*t17*t30*t62
     +        +4.751830440958629d-2*t14*t31*t39*t41*t56
     +        -1.1887928018410271d-4*t11*t39*t40*t56
     +        -1.680616229806896d-3*t15*t36*t39*t56
     +        +4.751830440958629d-2*t14*t31*t39*t41*t55
     +        -1.1887928018410271d-4*t11*t39*t40*t55
     +        -1.680616229806896d-3*t15*t36*t39*t55
     +        -1.680616229806896d-3*t31*t36*t39*t49
     +        -1.0083697378841376d-2*t15*t39*t40-3.361232459613792d-3
     +        *t28*t31*t32*t36*t39+6.722464919227584d-3*t15*t24*t32
     +        *t39+6.722464919227584d-3*t24*t28*t31*t39)
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t88-2.369632590246221d-2*t19*t3*t37*t51
     +        *t82+5.639913617919763d+0*rhoa*t19*t28*t50*t62
     +        -3.759942411946509d+0*t19*t27*t62+9.478530360984885d-3
     +        *t19*t3*t37*t49*t50+1.895706072196977d-2*t19*t28*t3*t36
     +        *t50)+1.329807601338106d-1*t18*t20*t7*t8*t88
     +        +2.659615202676212d-1*t20*t28*t66*t7*t8
     +        +2.659615202676212d-1*t20*t33*t62*t7*t8
     +        -1.329807601338106d-1*t18*t27*t28*t62*t7*t8
     +        -6.649038006690531d-2*t27*t42*t51*t7*t8
     +        +1.329807601338106d-1*t20*t42*t49*t7*t8
     +        -5.319230405352424d-1*t10*t13*t66*t7
     +        +8.939584949317599d-4*t18*t20*t3*t64*t7
     +        -2.659615202676212d-1*t13*t18*t20*t62*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t57*t7
     +        +3.3523443559941d-4*t18*t24*t3*t50*t51*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t3*t49*t7
     +        +5.319230405352424d-1*t10*t30*t42*t7
     +        -2.659615202676212d-1*t13*t20*t28*t42*t7
     +        -8.939584949317599d-4*t20*t3*t33*t40*t7
     +        +4.4697924746587997d-4*t18*t27*t28*t3*t40*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t33*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t65
     +        +3.5449077018110403d+0*t16*t26*t32*t63
     +        +1.6844800327854395d-1*t14*t26*t31*t39*t41*t52
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t32*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t88
     +        -3.5461536035682834d-1*t20*t28*t6*t67*t7*t8
     +        -3.5461536035682834d-1*t20*t34*t6*t62*t7*t8
     +        +1.7730768017841417d-1*t21*t27*t28*t6*t62*t7*t8
     +        +8.865384008920708d-2*t27*t43*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t43*t49*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t6*t67*t7
     +        -1.1919446599090135d-3*t20*t21*t3*t6*t64*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t6*t62*t7
     +        -5.959723299545067d-4*t20*t24*t3*t58*t6*t7
     +        -4.469792474658801d-4*t21*t24*t3*t50*t51*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t3*t49*t6*t7
     +        -7.092307207136567d-1*t10*t30*t43*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t28*t43*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t34*t40*t6*t7
     +        -5.959723299545067d-4*t21*t27*t28*t3*t40*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t34*t6*t7)+1.0d+0
     +        *t22*(-4.8246324681544644d-7*t3*t46*t87
     +        +2.550431938149111d-4*gammaaa*t3*t47
     +        -1.0457167652131799d-2*t3*t40)+2.0d+0*t25*t68+2.0d+0*t35
     +        *t61-3.361232459613792d-3*t3*t36*t59+1.0d+0*t44*t48)
     +        *wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(-3.536776513153214d-2
     +        *t11*t13*t16*t9*t94+t16*t94-3.536776513153214d-2*t11*t13
     +        *t17*t92-1.0681943420417696d-6*t11*t17*t3*t80
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41*t74
     +        -1.1887928018410271d-4*t11*t32*t39*t40*t74
     +        -1.680616229806896d-3*t15*t32*t36*t39*t74
     +        -1.680616229806896d-3*t28*t31*t36*t39*t74
     +        +3.361232459613792d-3*t15*t24*t39*t74
     +        -3.194401344146213d-4*t14*t32*t45*t54*t72
     +        +1.1297883647561385d-5*t2*t28*t54*t72
     +        +5.648941823780692d-6*t2*t31*t32*t72
     +        +9.503660881917259d-2*t14*t31*t39*t41*t65
     +        -2.3775856036820542d-4*t11*t39*t40*t65
     +        -3.361232459613792d-3*t15*t36*t39*t65
     +        -3.361232459613792d-3*t31*t36*t39*t62)
     +        +7.512005776780308d-7*t18*t27*t7*t70*t93
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t92-7.964885979694308d-5*t19*t28*t70*t71
     +        *t82+3.1859543918777233d-5*t19*t2*t50*t70
     +        +1.895706072196977d-2*t19*t3*t37*t50*t62
     +        +2.838989964383816d-5*t19*t27*t3*t53
     +        -4.2584849465757246d-5*t19*t28*t3*t30*t50)
     +        +1.329807601338106d-1*t18*t20*t7*t8*t92
     +        +1.329807601338106d-1*t20*t28*t7*t75*t8
     +        +2.659615202676212d-1*t20*t42*t62*t7*t8
     +        -2.659615202676212d-1*t10*t13*t7*t75
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41*t74
     +        -1.1323857927539448d-3*t14*t26*t29*t45*t54*t72
     +        +1.1268008665170462d-6*t18*t28*t50*t69*t7*t70
     +        -7.512005776780308d-7*t27*t33*t69*t7*t70
     +        +8.939584949317599d-4*t20*t24*t3*t66*t7
     +        -4.4697924746587997d-4*t18*t24*t27*t3*t62*t7
     +        -2.0081746054220725d-6*t20*t3*t33*t60*t7
     +        +1.0040873027110363d-6*t18*t27*t28*t3*t60*t7
     +        +2.0081746054220725d-6*t18*t20*t3*t47*t7
     +        -8.939584949317599d-4*t20*t3*t40*t42*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t42*t7
     +        +3.368960065570879d-1*t14*t26*t31*t39*t41*t63)
     +        -1.0016007702373746d-6*t21*t27*t6*t7*t70*t93
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t92
     +        -1.7730768017841417d-1*t20*t28*t6*t7*t76*t8
     +        -3.5461536035682834d-1*t20*t43*t6*t62*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t6*t7*t76
     +        -1.502401155356062d-6*t21*t28*t50*t6*t69*t7*t70
     +        +1.0016007702373746d-6*t27*t34*t6*t69*t7*t70
     +        -1.1919446599090135d-3*t20*t24*t3*t6*t67*t7
     +        +5.959723299545067d-4*t21*t24*t27*t3*t6*t62*t7
     +        +2.677566140562764d-6*t20*t3*t34*t6*t60*t7
     +        -1.338783070281382d-6*t21*t27*t28*t3*t6*t60*t7
     +        -2.677566140562764d-6*t20*t21*t3*t47*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t40*t43*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t43*t6*t7)+1.0d+0
     +        *t22*(1.8092371755579242d-7*gammaaa*t3*t91
     +        -6.040496695616315d-5*t3*t60)+1.0d+0*t25*t77
     +        -6.722464919227584d-3*t3*t36*t68+2.0d+0*t44*t61
     +        +1.510124173904079d-5*t3*t35*t53)*wght+Cmat3(iq,
     +        D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(4.751830440958629d-2
     +        *t14*t31*t39*t41*t99-1.680616229806896d-3*t15*t36*t39
     +        *t99-3.194401344146213d-4*t14*t45*t54*t72*t98
     +        +3.221137646093454d-6*t13*t14*t84*t97
     +        +2.5379391956869274d-8*t31*t72*t89-2.3995753309638476d-9
     +        *t11*t17*t3*t87-3.194401344146213d-4*t14*t45*t54*t72*t74
     +        +5.648941823780692d-6*t2*t31*t72*t74
     +        +4.751830440958629d-2*t100*t14*t31*t39*t41
     +        -1.680616229806896d-3*t100*t15*t36*t39)
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t99
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t98
     +        -9.490367802866481d-7*t19*t26*t30*t82*t97
     +        +1.0124912502475182d-8*t18*t27*t7*t70*t96
     +        +3.7874396480580436d-9*t18*t50*t60*t7*t95
     +        +9.022265057551551d-9*t18*t20*t3*t7*t91
     +        +1.3409377423976399d-3*t20*t24*t3*t7*t75
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t74
     +        -2.2536017330340924d-6*t27*t42*t69*t7*t70
     +        -6.024523816266217d-6*t20*t3*t42*t60*t7
     +        -2.24003435190034d-2*t100*t19*t26*t27*t37*t39)
     +        -1.3499883336633578d-8*t21*t27*t6*t7*t70*t96
     +        -5.049919530744059d-9*t21*t50*t6*t60*t7*t95
     +        -1.202968674340207d-8*t20*t21*t3*t6*t7*t91
     +        -1.7879169898635203d-3*t20*t24*t3*t6*t7*t76
     +        +3.0048023107121236d-6*t27*t43*t6*t69*t7*t70
     +        +8.032698421688291d-6*t20*t3*t43*t6*t60*t7)
     +        -1.0083697378841378d-2*t3*t36*t77+4.5303725217122376d-5
     +        *t3*t44*t53-6.784639408342217d-8*t22*t3/rhoa*
     +        *6.666666666666667d+0)*wght+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = -param(1)
            t7 = param(2)
            t8 = param(3)
            t9 = 1/rhoa
            t10 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t11 = t10**5.0d-1
            t12 = t8**2
            t13 = rhoa**2
            t14 = 1/t13
            t15 = 1/t12
            t16 = 1/t10
            t17 = exp(-1.413716694115413d+1*t13*t15*t16)
            t18 = t17-1.0d+0
            t19 = -3.536776513153214d-2*t10*t12*t14*t18+t17-1.5d+0
            t20 = 1/t8
            t21 = 1/t11
            t22 = 2.659615202676212d-1*t11*t19*t8*t9
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t20
     +        *t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9+t6+1.0d+0
            t24 = rhob**1.3333333333333333d+0
            t25 = 1/rhob**2.6666666666666666d+0
            t26 = exp(-4.49276922095889d-3*gammabb*t25)
            t27 = 1.0d+0-t26
            t28 = -7.481426920246764d-1*t24*t27-9.305257363491d-1*t24
            t29 = 1/rhob
            t30 = 7.481426920246764d-1*t24*t27+9.305257363491d-1*t24
            t31 = t30**5.0d-1
            t32 = rhob**2
            t33 = 1/t32
            t34 = 1/t30
            t35 = exp(-1.413716694115413d+1*t15*t32*t34)
            t36 = t35-1.0d+0
            t37 = -3.536776513153214d-2*t12*t30*t33*t36+t35-1.5d+0
            t38 = 1/t31
            t39 = 2.659615202676212d-1*t29*t31*t37*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t20
     +        *t38)
            t40 = -3.5461536035682834d-1*t29*t31*t39*t7*t8+t6+1.0d+0
            t41 = rhoa**3.333333333333333d-1
            t42 = 1/rhoa**2.3333333333333334d+0
            t43 = 8.963286558970113d-3*gammaaa*t3*t42
     +        -9.975235893662352d-1*t4*t41-1.2407009817987999d+0*t41
            t44 = 5.641895835477563d-1
            t45 = 1/t11**3
            t46 = -8.963286558970113d-3*gammaaa*t3*t42
     +        +9.975235893662352d-1*t4*t41+1.2407009817987999d+0*t41
            t47 = 3.759942411946509d+0*t20*t21-1.8799712059732546d+0
     +        *rhoa*t20*t45*t46
            t48 = 1/rhoa**3
            t49 = 1/t10**2
            t50 = 1.413716694115413d+1*t13*t15*t46*t49
     +        -2.827433388230826d+1*rhoa*t15*t16
            t51 = -3.536776513153214d-2*t10*t12*t14*t17*t50+t17*t50
     +        +7.073553026306428d-2*t10*t12*t18*t48
     +        -3.536776513153214d-2*t12*t14*t18*t46
            t52 = 2.659615202676212d-1*t11*t51*t8*t9
     +        +1.329807601338106d-1*t19*t21*t46*t8*t9
     +        -2.659615202676212d-1*t11*t14*t19*t8
     +        +3.5449077018110403d+0*t17*t44*t47
            t53 = -3.5461536035682834d-1*t11*t52*t7*t8*t9
     +        -1.7730768017841417d-1*t21*t22*t46*t7*t8*t9
     +        +3.5461536035682834d-1*t11*t14*t22*t7*t8
            t54 = rhob**3.333333333333333d-1
            t55 = 1/rhob**2.3333333333333334d+0
            t56 = 8.963286558970113d-3*gammabb*t26*t55
     +        -9.975235893662352d-1*t27*t54-1.2407009817987999d+0*t54
            t57 = 1/t31**3
            t58 = -8.963286558970113d-3*gammabb*t26*t55
     +        +9.975235893662352d-1*t27*t54+1.2407009817987999d+0*t54
            t59 = 3.759942411946509d+0*t20*t38-1.8799712059732546d+0
     +        *rhob*t20*t57*t58
            t60 = 1/rhob**3
            t61 = 1/t30**2
            t62 = 1.413716694115413d+1*t15*t32*t58*t61
     +        -2.827433388230826d+1*rhob*t15*t34
            t63 = -3.536776513153214d-2*t12*t30*t33*t35*t62+t35*t62
     +        +7.073553026306428d-2*t12*t30*t36*t60
     +        -3.536776513153214d-2*t12*t33*t36*t58
            t64 = 2.659615202676212d-1*t29*t31*t63*t8
     +        +1.329807601338106d-1*t29*t37*t38*t58*t8
     +        -2.659615202676212d-1*t31*t33*t37*t8
     +        +3.5449077018110403d+0*t35*t44*t59
            t65 = -3.5461536035682834d-1*t29*t31*t64*t7*t8
     +        -1.7730768017841417d-1*t29*t38*t39*t58*t7*t8
     +        +3.5461536035682834d-1*t31*t33*t39*t7*t8
            t66 = 1/t1
            t67 = 1/t41
            t68 = -1.413716694115413d+1*t13*t15*t16
            t69 = exp(t68-4.49276922095889d-3*gammaaa*t2)
            t70 = 1/rhoa**3.3333333333333337d+0
            t71 = rhoa**6.666666666666666d-1
            t72 = 4.751830440958629d-2*t15*t49*t69*t71
     +        -1.1887928018410271d-4*t12*t18*t3*t70
     +        -1.680616229806896d-3*t16*t66*t69
            t73 = 2.659615202676212d-1*t11*t72*t8*t9
     +        +4.4697924746587997d-4*t19*t21*t3*t42*t8
     +        -2.24003435190034d-2*t20*t44*t45*t67*t69
            t74 = -3.5461536035682834d-1*t11*t7*t73*t8*t9
     +        -5.959723299545067d-4*t21*t22*t3*t42*t7*t8
            t75 = 1/t24
            t76 = 1/t54
            t77 = -1.413716694115413d+1*t15*t32*t34
            t78 = exp(t77-4.49276922095889d-3*gammabb*t25)
            t79 = 1/rhob**3.3333333333333337d+0
            t80 = rhob**6.666666666666666d-1
            t81 = 4.751830440958629d-2*t15*t61*t78*t80
     +        -1.1887928018410271d-4*t12*t26*t36*t79
     +        -1.680616229806896d-3*t34*t75*t78
            t82 = 2.659615202676212d-1*t29*t31*t8*t81
     +        +4.4697924746587997d-4*t26*t37*t38*t55*t8
     +        -2.24003435190034d-2*t20*t44*t57*t76*t78
            t83 = -3.5461536035682834d-1*t29*t31*t7*t8*t82
     +        -5.959723299545067d-4*t26*t38*t39*t55*t7*t8
            t84 = 1/t71
            t85 = gammaaa**2
            t86 = 1/rhoa**6
            t87 = 1.0738660792206783d-4*t3*t85*t86
     +        -3.3250786312207836d-1*t4*t84-4.135669939329333d-1*t84
     +        -8.963286558970114d-3*gammaaa*t3*t70
            t88 = -1.0738660792206783d-4*t3*t85*t86
     +        +3.3250786312207836d-1*t4*t84+4.135669939329333d-1*t84
     +        +8.963286558970114d-3*gammaaa*t3*t70
            t89 = 1/t11**5
            t90 = t46**2
            t91 = 2.8199568089598814d+0*rhoa*t20*t89*t90
     +        -1.8799712059732546d+0*rhoa*t20*t45*t88
     +        -3.759942411946509d+0*t20*t45*t46
            t92 = 1/rhoa**4
            t93 = 1/t10**3
            t94 = -2.827433388230826d+1*t13*t15*t90*t93
     +        +1.413716694115413d+1*t13*t15*t49*t88
     +        +5.654866776461652d+1*rhoa*t15*t46*t49
     +        -2.827433388230826d+1*t15*t16
            t95 = t50**2
            t96 = -3.536776513153214d-2*t10*t12*t14*t17*t95+t17*t95
     +        -3.536776513153214d-2*t10*t12*t14*t17*t94+t17*t94
     +        -2.1220659078919285d-1*t10*t12*t18*t92
     +        -3.536776513153214d-2*t12*t14*t18*t88
     +        +1.4147106052612857d-1*t10*t12*t17*t48*t50
     +        -7.073553026306428d-2*t12*t14*t17*t46*t50
     +        +1.4147106052612857d-1*t12*t18*t46*t48
            t97 = 2.659615202676212d-1*t11*t8*t9*t96
     +        +3.5449077018110403d+0*t17*t44*t91-6.649038006690531d-2
     +        *t19*t45*t8*t9*t90+1.329807601338106d-1*t19*t21*t8*t88
     +        *t9+2.659615202676212d-1*t21*t46*t51*t8*t9
     +        -5.319230405352424d-1*t11*t14*t51*t8
     +        +5.319230405352424d-1*t11*t19*t48*t8
     +        -2.659615202676212d-1*t14*t19*t21*t46*t8
     +        +3.5449077018110403d+0*t17*t44*t47*t50
            t98 = -3.5461536035682834d-1*t11*t7*t8*t9*t97
     +        +8.865384008920708d-2*t22*t45*t7*t8*t9*t90
     +        -1.7730768017841417d-1*t21*t22*t7*t8*t88*t9
     +        -3.5461536035682834d-1*t21*t46*t52*t7*t8*t9
     +        +7.092307207136567d-1*t11*t14*t52*t7*t8
     +        -7.092307207136567d-1*t11*t22*t48*t7*t8
     +        +3.5461536035682834d-1*t14*t21*t22*t46*t7*t8
            t99 = 1/t80
            t100 = gammabb**2
            t101 = 1/rhob**6
            t102 = -3.3250786312207836d-1*t27*t99-4.135669939329333d-1
     +        *t99-8.963286558970114d-3*gammabb*t26*t79
     +        +1.0738660792206783d-4*t100*t101*t26
            t103 = 3.3250786312207836d-1*t27*t99+4.135669939329333d-1
     +        *t99+8.963286558970114d-3*gammabb*t26*t79
     +        -1.0738660792206783d-4*t100*t101*t26
            t104 = 1/t31**5
            t105 = t58**2
            t106 = -3.759942411946509d+0*t20*t57*t58
     +        -1.8799712059732546d+0*rhob*t103*t20*t57
     +        +2.8199568089598814d+0*rhob*t104*t105*t20
            t107 = 1/rhob**4
            t108 = 1/t30**3
            t109 = 5.654866776461652d+1*rhob*t15*t58*t61
     +        +1.413716694115413d+1*t103*t15*t32*t61
     +        -2.827433388230826d+1*t15*t34-2.827433388230826d+1*t105
     +        *t108*t15*t32
            t110 = t62**2
            t111 = 1.4147106052612857d-1*t12*t30*t35*t60*t62
     +        -7.073553026306428d-2*t12*t33*t35*t58*t62
     +        +1.4147106052612857d-1*t12*t36*t58*t60
     +        -3.536776513153214d-2*t103*t12*t33*t36
     +        -2.1220659078919285d-1*t107*t12*t30*t36
     +        -3.536776513153214d-2*t110*t12*t30*t33*t35
     +        -3.536776513153214d-2*t109*t12*t30*t33*t35+t110*t35+t109
     +        *t35
            t112 = 2.659615202676212d-1*t29*t38*t58*t63*t8
     +        -5.319230405352424d-1*t31*t33*t63*t8
     +        +5.319230405352424d-1*t31*t37*t60*t8
     +        -2.659615202676212d-1*t33*t37*t38*t58*t8
     +        -6.649038006690531d-2*t105*t29*t37*t57*t8
     +        +1.329807601338106d-1*t103*t29*t37*t38*t8
     +        +2.659615202676212d-1*t111*t29*t31*t8
     +        +3.5449077018110403d+0*t35*t44*t59*t62
     +        +3.5449077018110403d+0*t106*t35*t44
            t113 = -3.5461536035682834d-1*t29*t38*t58*t64*t7*t8
     +        +7.092307207136567d-1*t31*t33*t64*t7*t8
     +        -7.092307207136567d-1*t31*t39*t60*t7*t8
     +        +3.5461536035682834d-1*t33*t38*t39*t58*t7*t8
     +        +8.865384008920708d-2*t105*t29*t39*t57*t7*t8
     +        -1.7730768017841417d-1*t103*t29*t38*t39*t7*t8
     +        -3.5461536035682834d-1*t112*t29*t31*t7*t8
            t114 = 1/rhoa**5
            t115 = 4.4816432794850564d-3*t3*t42-4.026997797077544d-5
     +        *gammaaa*t114*t3
            t116 = 4.026997797077544d-5*gammaaa*t114*t3
     +        -4.4816432794850564d-3*t3*t42
            t117 = 9.478530360984885d-3*t20*t3*t46*t67*t89
     +        -6.31902024065659d-3*t20*t3*t45*t66
     +        -1.8799712059732546d+0*rhoa*t116*t20*t45
            t118 = 1/rhoa**4.333333333333333d+0
            t119 = -9.503660881917259d-2*t15*t3*t46*t71*t93
     +        +9.503660881917259d-2*t15*t3*t49*t67
     +        +1.413716694115413d+1*t116*t13*t15*t49
            t120 = 4.751830440958629d-2*t15*t49*t50*t69*t71
     +        -1.1887928018410271d-4*t12*t50*t69*t70
     +        -1.680616229806896d-3*t16*t50*t66*t69
     +        -1.680616229806896d-3*t46*t49*t66*t69
     +        +3.361232459613792d-3*t16*t42*t69+2.3775856036820542d-4
     +        *t118*t12*t18*t3-3.536776513153214d-2*t116*t12*t14*t18
     +        -3.536776513153214d-2*t10*t119*t12*t14*t17+t119*t17
            t121 = 1.329807601338106d-1*t21*t46*t72*t8*t9
     +        +1.329807601338106d-1*t116*t19*t21*t8*t9
     +        +2.659615202676212d-1*t11*t120*t8*t9
     +        -2.659615202676212d-1*t11*t14*t72*t8
     +        -4.4697924746587997d-4*t19*t21*t3*t70*t8
     +        +4.4697924746587997d-4*t21*t3*t42*t51*t8
     +        -2.2348962373293998d-4*t19*t3*t42*t45*t46*t8
     +        +1.6844800327854395d-1*t15*t44*t47*t49*t69*t71
     +        +3.5449077018110403d+0*t117*t17*t44
            t122 = -1.7730768017841417d-1*t21*t46*t7*t73*t8*t9
     +        -1.7730768017841417d-1*t116*t21*t22*t7*t8*t9
     +        -3.5461536035682834d-1*t11*t121*t7*t8*t9
     +        +3.5461536035682834d-1*t11*t14*t7*t73*t8
     +        +5.959723299545067d-4*t21*t22*t3*t7*t70*t8
     +        -5.959723299545067d-4*t21*t3*t42*t52*t7*t8
     +        +2.9798616497725333d-4*t22*t3*t42*t45*t46*t7*t8
            t123 = 1/rhob**5
            t124 = 4.4816432794850564d-3*t26*t55-4.026997797077544d-5
     +        *gammabb*t123*t26
            t125 = 4.026997797077544d-5*gammabb*t123*t26
     +        -4.4816432794850564d-3*t26*t55
            t126 = 9.478530360984885d-3*t104*t20*t26*t58*t76
     +        -6.31902024065659d-3*t20*t26*t57*t75
     +        -1.8799712059732546d+0*rhob*t125*t20*t57
            t127 = 1/rhob**4.333333333333333d+0
            t128 = -9.503660881917259d-2*t108*t15*t26*t58*t80
     +        +9.503660881917259d-2*t15*t26*t61*t76
     +        +1.413716694115413d+1*t125*t15*t32*t61
            t129 = 4.751830440958629d-2*t15*t61*t62*t78*t80
     +        -1.1887928018410271d-4*t12*t62*t78*t79
     +        -1.680616229806896d-3*t34*t62*t75*t78
     +        -1.680616229806896d-3*t58*t61*t75*t78
     +        +3.361232459613792d-3*t34*t55*t78-3.536776513153214d-2
     +        *t12*t125*t33*t36+2.3775856036820542d-4*t12*t127*t26*t36
     +        -3.536776513153214d-2*t12*t128*t30*t33*t35+t128*t35
            t130 = 1.329807601338106d-1*t29*t38*t58*t8*t81
     +        -2.659615202676212d-1*t31*t33*t8*t81
     +        +1.6844800327854395d-1*t15*t44*t59*t61*t78*t80
     +        -4.4697924746587997d-4*t26*t37*t38*t79*t8
     +        +4.4697924746587997d-4*t26*t38*t55*t63*t8
     +        -2.2348962373293998d-4*t26*t37*t55*t57*t58*t8
     +        +1.329807601338106d-1*t125*t29*t37*t38*t8
     +        +2.659615202676212d-1*t129*t29*t31*t8
     +        +3.5449077018110403d+0*t126*t35*t44
            t131 = -1.7730768017841417d-1*t29*t38*t58*t7*t8*t82
     +        +3.5461536035682834d-1*t31*t33*t7*t8*t82
     +        +5.959723299545067d-4*t26*t38*t39*t7*t79*t8
     +        -5.959723299545067d-4*t26*t38*t55*t64*t7*t8
     +        +2.9798616497725333d-4*t26*t39*t55*t57*t58*t7*t8
     +        -1.7730768017841417d-1*t125*t29*t38*t39*t7*t8
     +        -3.5461536035682834d-1*t130*t29*t31*t7*t8
            t132 = 1/rhoa**3.6666666666666664d+0
            t133 = exp(-8.98553844191778d-3*gammaaa*t2)
            t134 = 1/rhoa**1.6666666666666669d+0
            t135 = exp(t68-8.98553844191778d-3*gammaaa*t2)
            t136 = 4.751830440958629d-2*t15*t3*t49*t71
            t137 = t136-4.49276922095889d-3*t2
            t138 = -3.194401344146213d-4*t135*t15*t84*t93
     +        +5.340971710208847d-7*t12*t18*t3*t86
     +        +4.751830440958629d-2*t137*t15*t49*t69*t71
     +        -1.680616229806896d-3*t137*t16*t66*t69
            t139 = 2.659615202676212d-1*t11*t138*t8*t9
     +        +1.1293914261386051d-4*t134*t135*t20*t44*t89
     +        +8.939584949317599d-4*t21*t3*t42*t72*t8
     +        -7.512005776780308d-7*t132*t133*t19*t45*t8
     +        -2.0081746054220725d-6*t114*t19*t21*t3*t8
     +        -2.24003435190034d-2*t137*t20*t44*t45*t67*t69
            t140 = -3.5461536035682834d-1*t11*t139*t7*t8*t9
     +        -1.1919446599090135d-3*t21*t3*t42*t7*t73*t8
     +        +1.0016007702373746d-6*t132*t133*t22*t45*t7*t8
     +        +2.677566140562764d-6*t114*t21*t22*t3*t7*t8
            t141 = 1/rhob**3.6666666666666664d+0
            t142 = exp(-8.98553844191778d-3*gammabb*t25)
            t143 = 1/rhob**1.6666666666666669d+0
            t144 = exp(t77-8.98553844191778d-3*gammabb*t25)
            t145 = 4.751830440958629d-2*t15*t26*t61*t80
            t146 = t145-4.49276922095889d-3*t25
            t147 = -3.194401344146213d-4*t108*t144*t15*t99
     +        +4.751830440958629d-2*t146*t15*t61*t78*t80
     +        -1.680616229806896d-3*t146*t34*t75*t78
     +        +5.340971710208847d-7*t101*t12*t26*t36
            t148 = 8.939584949317599d-4*t26*t38*t55*t8*t81
     +        -7.512005776780308d-7*t141*t142*t37*t57*t8
     +        -2.0081746054220725d-6*t123*t26*t37*t38*t8
     +        +2.659615202676212d-1*t147*t29*t31*t8
     +        -2.24003435190034d-2*t146*t20*t44*t57*t76*t78
     +        +1.1293914261386051d-4*t104*t143*t144*t20*t44
            t149 = -1.1919446599090135d-3*t26*t38*t55*t7*t8*t82
     +        +1.0016007702373746d-6*t141*t142*t39*t57*t7*t8
     +        +2.677566140562764d-6*t123*t26*t38*t39*t7*t8
     +        -3.5461536035682834d-1*t148*t29*t31*t7*t8
            t150 = gammaaa**3
            t151 = 1/rhoa**9.666666666666666d+0
            t152 = 1/rhoa**7
            t153 = 7.517062554544749d-4*t152*t3*t85
     +        -2.2167190874805223d-1*t134*t4-1.286568658174524d-6*t150
     +        *t151*t3-3.386130477833154d-2*gammaaa*t118*t3
     +        -2.757113292886222d-1*t134
            t154 = 1/t11**7
            t155 = t46**3
            t156 = 1/t10**4
            t157 = -1.6964600329384957d+2*rhoa*t15*t90*t93
     +        -8.482300164692478d+1*t13*t15*t46*t88*t93
     +        +8.482300164692478d+1*rhoa*t15*t49*t88
     +        +8.482300164692478d+1*t15*t46*t49+1.413716694115413d+1
     +        *t13*t15*t153*t49+8.482300164692478d+1*t13*t15*t155*t156
            t158 = t50**3
            t159 = gammabb**3
            t160 = 1/rhob**9.666666666666666d+0
            t161 = 1/rhob**7
            t162 = -2.2167190874805223d-1*t143*t27
     +        +7.517062554544749d-4*t100*t161*t26-1.286568658174524d-6
     +        *t159*t160*t26-3.386130477833154d-2*gammabb*t127*t26
     +        -2.757113292886222d-1*t143
            t163 = 1/t31**7
            t164 = t58**3
            t165 = 1/t30**4
            t166 = 8.482300164692478d+1*t15*t58*t61
     +        +1.413716694115413d+1*t15*t162*t32*t61
     +        +8.482300164692478d+1*rhob*t103*t15*t61
     +        -8.482300164692478d+1*t103*t108*t15*t32*t58
     +        +8.482300164692478d+1*t15*t164*t165*t32
     +        -1.6964600329384957d+2*rhob*t105*t108*t15
            t167 = t62**3
            t168 = 1/rhoa**8.666666666666666d+0
            t169 = -2.550431938149111d-4*gammaaa*t3*t86
     +        +4.8246324681544644d-7*t168*t3*t85+1.0457167652131799d-2
     +        *t3*t70
            t170 = 1/rhoa**5.333333333333333d+0
            t171 = -9.503660881917259d-2*t15*t3*t71*t88*t93
     +        -3.8014643527669034d-1*t15*t3*t46*t67*t93
     +        -5.654866776461652d+1*t116*t13*t15*t46*t93
     +        +2.8510982645751776d-1*t15*t156*t3*t71*t90
     +        +9.503660881917259d-2*t15*t3*t49*t66
     +        +1.413716694115413d+1*t13*t15*t169*t49
     +        +5.654866776461652d+1*rhoa*t116*t15*t49
            t172 = 1/rhob**8.666666666666666d+0
            t173 = 1.0457167652131799d-2*t26*t79+4.8246324681544644d-7
     +        *t100*t172*t26-2.550431938149111d-4*gammabb*t101*t26
            t174 = 1/rhob**5.333333333333333d+0
            t175 = 2.8510982645751776d-1*t105*t15*t165*t26*t80
     +        -9.503660881917259d-2*t103*t108*t15*t26*t80
     +        -3.8014643527669034d-1*t108*t15*t26*t58*t76
     +        +9.503660881917259d-2*t15*t26*t61*t75
     +        +1.413716694115413d+1*t15*t173*t32*t61
     +        +5.654866776461652d+1*rhob*t125*t15*t61
     +        -5.654866776461652d+1*t108*t125*t15*t32*t58
            t176 = 1/rhoa**7.666666666666667d+0
            t177 = 6.040496695616315d-5*t114*t3-1.8092371755579242d-7
     +        *gammaaa*t176*t3
            t178 = 1/rhoa**4.666666666666667d+0
            t179 = -1.9007321763834517d-1*t116*t15*t3*t71*t93
     +        +4.269775509670888d-4*t14*t15*t3*t46*t93
     +        -6.388802688292426d-4*t133*t134*t15*t93
     +        +9.583204032438639d-4*t133*t15*t156*t46*t84
     +        -4.269775509670888d-4*t15*t3*t48*t49
     +        +1.413716694115413d+1*t13*t15*t177*t49
            t180 = 1/rhob**7.666666666666667d+0
            t181 = 6.040496695616315d-5*t123*t26-1.8092371755579242d-7
     +        *gammabb*t180*t26
            t182 = 1/rhob**4.666666666666667d+0
            t183 = 9.583204032438639d-4*t142*t15*t165*t58*t99
     +        -1.9007321763834517d-1*t108*t125*t15*t26*t80
     +        -4.269775509670888d-4*t15*t26*t60*t61
     +        +1.413716694115413d+1*t15*t181*t32*t61
     +        +4.269775509670888d-4*t108*t15*t26*t33*t58
     +        -6.388802688292426d-4*t108*t142*t143*t15
            t184 = exp(-1.3478307662876668d-2*gammaaa*t2)
            t185 = 1/rhoa**6.333333333333333d+0
            t186 = exp(t68-1.3478307662876668d-2*gammaaa*t2)
            t187 = t136-8.98553844191778d-3*t2
            t188 = -3.194401344146213d-4*t133*t15*t84*t93
     +        -2.134887754835444d-4*t14*t15*t3*t49
            t189 = t137**2
            t190 = exp(-1.3478307662876668d-2*gammabb*t25)
            t191 = 1/rhob**6.333333333333333d+0
            t192 = exp(t77-1.3478307662876668d-2*gammabb*t25)
            t193 = t145-8.98553844191778d-3*t25
            t194 = -3.194401344146213d-4*t108*t142*t15*t99
     +        -2.134887754835444d-4*t15*t26*t33*t61
            t195 = t146**2
            fnc(iq) = (1.0d+0*t23*t5+1.0d+0*t28*t40)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t5*t53+1.0d+0*t23*t43)*wght+Amat
     +        (iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t28*t65+1.0d+0*t40*t56)*wght+Amat
     +        (iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t5*t74-3.361232459613792d-3*t23
     +        *t3*t66)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t28*t83-3.361232459613792d-3*t26
     +        *t40*t75)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (1.0d+0*t5*t98+1.0d+0*t23*t87+2.0d+0
     +        *t43*t53)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (2.0d+0*t56*t65+1.0d+0*t102*t40
     +        +1.0d+0*t113*t28)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (1.0d+0*t43*t74-3.361232459613792d-3
     +        *t3*t53*t66+1.0d+0*t122*t5+1.0d+0*t115*t23)*wght+Cmat2
     +        (iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (1.0d+0*t56*t83-3.361232459613792d-3
     +        *t26*t65*t75+1.0d+0*t124*t40+1.0d+0*t131*t28)*wght+Cmat2
     +        (iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.510124173904079d-5*t23*t3*t92
     +        -6.722464919227584d-3*t3*t66*t74+1.0d+0*t140*t5)*wght
     +        +Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (-6.722464919227584d-3*t26*t75*t83
     +        +1.510124173904079d-5*t107*t26*t40+1.0d+0*t149*t28)*wght
     +        +Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (3.0d+0*t43*t98+1.0d+0*t5*(
     +        -5.319230405352425d-1*t21*t46*t7*t8*t9*t97
     +        +1.063846081070485d+0*t11*t14*t7*t8*t97
     +        -3.5461536035682834d-1*t11*t7*t8*t9*
     +        (3.989422804014318d-1*t21*t46*t8*t9*t96
     +        -7.978845608028636d-1*t11*t14*t8*t96
     +        +2.659615202676212d-1*t11*t8*t9*(2.1220659078919285d-1
     +        *t10*t12*t17*t48*t95-1.0610329539459643d-1*t12*t14*t17
     +        *t46*t95-1.0610329539459643d-1*t10*t12*t14*t17*t50*t94+3
     +        *t17*t50*t94+2.1220659078919285d-1*t10*t12*t17*t48*t94
     +        -1.0610329539459643d-1*t12*t14*t17*t46*t94
     +        -6.366197723675786d-1*t10*t12*t17*t50*t92
     +        -6.366197723675786d-1*t12*t18*t46*t92
     +        -1.0610329539459643d-1*t12*t14*t17*t50*t88
     +        +2.1220659078919285d-1*t12*t18*t48*t88
     +        +4.244131815783857d-1*t12*t17*t46*t48*t50
     +        -3.536776513153214d-2*t12*t14*t153*t18
     +        +8.488263631567714d-1*t10*t114*t12*t18
     +        -3.536776513153214d-2*t10*t12*t14*t158*t17+t158*t17
     +        -3.536776513153214d-2*t10*t12*t14*t157*t17+t157*t17)
     +        +3.5449077018110403d+0*t17*t44*t47*t95
     +        +3.5449077018110403d+0*t17*t44*t47*t94
     +        -1.5957691216057274d+0*t11*t19*t8*t92
     +        +7.089815403622081d+0*t17*t44*t50*t91
     +        +3.5449077018110403d+0*t17*t44*(8.459870426879645d+0*t20
     +        *t89*t90+8.459870426879645d+0*rhoa*t20*t46*t88*t89
     +        -5.639913617919763d+0*t20*t45*t88-1.8799712059732546d+0
     +        *rhoa*t153*t20*t45-7.049892022399704d+0*rhoa*t154*t155
     +        *t20)-1.994711402007159d-1*t45*t51*t8*t9*t90
     +        +1.994711402007159d-1*t14*t19*t45*t8*t90
     +        +9.973557010035795d-2*t155*t19*t8*t89*t9
     +        +3.989422804014318d-1*t21*t51*t8*t88*t9
     +        -1.994711402007159d-1*t19*t45*t46*t8*t88*t9
     +        +1.329807601338106d-1*t153*t19*t21*t8*t9
     +        -3.989422804014318d-1*t14*t19*t21*t8*t88
     +        +1.5957691216057274d+0*t11*t48*t51*t8
     +        -7.978845608028636d-1*t14*t21*t46*t51*t8
     +        +7.978845608028636d-1*t19*t21*t46*t48*t8)
     +        +2.12769216214097d+0*t11*t22*t7*t8*t92
     +        +2.6596152026762127d-1*t45*t52*t7*t8*t9*t90
     +        -2.6596152026762127d-1*t14*t22*t45*t7*t8*t90
     +        -1.3298076013381063d-1*t155*t22*t7*t8*t89*t9
     +        -5.319230405352425d-1*t21*t52*t7*t8*t88*t9
     +        +2.6596152026762127d-1*t22*t45*t46*t7*t8*t88*t9
     +        -1.7730768017841417d-1*t153*t21*t22*t7*t8*t9
     +        +5.319230405352425d-1*t14*t21*t22*t7*t8*t88
     +        -2.12769216214097d+0*t11*t48*t52*t7*t8
     +        +1.063846081070485d+0*t14*t21*t46*t52*t7*t8
     +        -1.063846081070485d+0*t21*t22*t46*t48*t7*t8)+3.0d+0*t53
     +        *t87+1.0d+0*t23*(-7.517062554544749d-4*t152*t3*t85
     +        +2.2167190874805223d-1*t134*t4+1.286568658174524d-6*t150
     +        *t151*t3+3.386130477833154d-2*gammaaa*t118*t3
     +        +2.757113292886222d-1*t134))*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (1.0d+0*t28*(
     +        -3.5461536035682834d-1*t29*t31*t7*t8*
     +        (1.5957691216057274d+0*t31*t60*t63*t8
     +        -7.978845608028636d-1*t33*t38*t58*t63*t8
     +        -1.994711402007159d-1*t105*t29*t57*t63*t8
     +        +3.989422804014318d-1*t103*t29*t38*t63*t8
     +        +2.659615202676212d-1*t29*t31*(4.244131815783857d-1*t12
     +        *t35*t58*t60*t62-1.0610329539459643d-1*t109*t12*t30*t33
     +        *t35*t62-1.0610329539459643d-1*t103*t12*t33*t35*t62
     +        -6.366197723675786d-1*t107*t12*t30*t35*t62+3*t109*t35
     +        *t62+2.1220659078919285d-1*t103*t12*t36*t60
     +        +2.1220659078919285d-1*t110*t12*t30*t35*t60
     +        +2.1220659078919285d-1*t109*t12*t30*t35*t60
     +        -6.366197723675786d-1*t107*t12*t36*t58
     +        -1.0610329539459643d-1*t110*t12*t33*t35*t58
     +        -1.0610329539459643d-1*t109*t12*t33*t35*t58
     +        -3.536776513153214d-2*t12*t162*t33*t36
     +        +8.488263631567714d-1*t12*t123*t30*t36
     +        -3.536776513153214d-2*t12*t167*t30*t33*t35
     +        -3.536776513153214d-2*t12*t166*t30*t33*t35+t167*t35+t166
     +        *t35)*t8+7.978845608028636d-1*t37*t38*t58*t60*t8
     +        -1.994711402007159d-1*t103*t29*t37*t57*t58*t8
     +        +3.989422804014318d-1*t111*t29*t38*t58*t8
     +        +1.994711402007159d-1*t105*t33*t37*t57*t8
     +        -3.989422804014318d-1*t103*t33*t37*t38*t8
     +        +1.329807601338106d-1*t162*t29*t37*t38*t8
     +        -1.5957691216057274d+0*t107*t31*t37*t8
     +        +9.973557010035795d-2*t104*t164*t29*t37*t8
     +        -7.978845608028636d-1*t111*t31*t33*t8
     +        +7.089815403622081d+0*t106*t35*t44*t62
     +        +3.5449077018110403d+0*t110*t35*t44*t59
     +        +3.5449077018110403d+0*t109*t35*t44*t59
     +        +3.5449077018110403d+0*t35*t44*(8.459870426879645d+0
     +        *rhob*t103*t104*t20*t58-1.8799712059732546d+0*rhob*t162
     +        *t20*t57-5.639913617919763d+0*t103*t20*t57
     +        -7.049892022399704d+0*rhob*t163*t164*t20
     +        +8.459870426879645d+0*t104*t105*t20))
     +        -2.12769216214097d+0*t31*t60*t64*t7*t8
     +        +1.063846081070485d+0*t33*t38*t58*t64*t7*t8
     +        +2.6596152026762127d-1*t105*t29*t57*t64*t7*t8
     +        -5.319230405352425d-1*t103*t29*t38*t64*t7*t8
     +        -1.063846081070485d+0*t38*t39*t58*t60*t7*t8
     +        +2.6596152026762127d-1*t103*t29*t39*t57*t58*t7*t8
     +        -5.319230405352425d-1*t112*t29*t38*t58*t7*t8
     +        -2.6596152026762127d-1*t105*t33*t39*t57*t7*t8
     +        +5.319230405352425d-1*t103*t33*t38*t39*t7*t8
     +        -1.7730768017841417d-1*t162*t29*t38*t39*t7*t8
     +        +2.12769216214097d+0*t107*t31*t39*t7*t8
     +        -1.3298076013381063d-1*t104*t164*t29*t39*t7*t8
     +        +1.063846081070485d+0*t112*t31*t33*t7*t8)+3.0d+0*t102
     +        *t65+3.0d+0*t113*t56+1.0d+0*(2.2167190874805223d-1*t143
     +        *t27-7.517062554544749d-4*t100*t161*t26
     +        +1.286568658174524d-6*t159*t160*t26+3.386130477833154d-2
     +        *gammabb*t127*t26+2.757113292886222d-1*t143)*t40)*wght
     +        +Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-3.361232459613792d-3*t3*t66*t98
     +        +1.0d+0*t5*(-5.959723299545067d-4*t21*t3*t42*t7*t8*t97
     +        -3.5461536035682834d-1*t11*t7*t8*t9*
     +        (4.4697924746587997d-4*t21*t3*t42*t8*t96
     +        +2.659615202676212d-1*t11*t8*t9*(4.751830440958629d-2
     +        *t15*t49*t69*t71*t95-1.1887928018410271d-4*t12*t69*t70
     +        *t95-1.680616229806896d-3*t16*t66*t69*t95
     +        +4.751830440958629d-2*t15*t49*t69*t71*t94
     +        -1.1887928018410271d-4*t12*t69*t70*t94
     +        -1.680616229806896d-3*t16*t66*t69*t94
     +        -1.680616229806896d-3*t49*t66*t69*t88
     +        -1.0083697378841376d-2*t16*t69*t70-3.361232459613792d-3
     +        *t46*t49*t50*t66*t69+6.722464919227584d-3*t16*t42*t50
     +        *t69+4.7551712073641084d-4*t118*t12*t50*t69
     +        +6.722464919227584d-3*t42*t46*t49*t69
     +        -7.073553026306428d-2*t10*t119*t12*t14*t17*t50
     +        -7.073553026306428d-2*t116*t12*t14*t17*t50+2*t119*t17
     +        *t50+1.4147106052612857d-1*t116*t12*t18*t48
     +        +1.4147106052612857d-1*t10*t119*t12*t17*t48
     +        -7.073553026306428d-2*t119*t12*t14*t17*t46
     +        -7.132756811046163d-4*t12*t170*t18*t3
     +        -3.536776513153214d-2*t12*t14*t169*t18
     +        -3.536776513153214d-2*t10*t12*t14*t17*t171+t17*t171)
     +        +1.6844800327854395d-1*t15*t44*t49*t69*t71*t91
     +        +3.5449077018110403d+0*t17*t44*(-2.369632590246221d-2
     +        *t154*t20*t3*t67*t90+9.478530360984885d-3*t20*t3*t67*t88
     +        *t89+1.895706072196977d-2*t20*t3*t46*t66*t89
     +        +5.639913617919763d+0*rhoa*t116*t20*t46*t89
     +        -1.8799712059732546d+0*rhoa*t169*t20*t45
     +        -3.759942411946509d+0*t116*t20*t45)-6.649038006690531d-2
     +        *t45*t72*t8*t9*t90+3.3523443559941d-4*t19*t3*t42*t8*t89
     +        *t90+1.329807601338106d-1*t21*t72*t8*t88*t9
     +        +2.659615202676212d-1*t116*t21*t51*t8*t9
     +        -1.329807601338106d-1*t116*t19*t45*t46*t8*t9
     +        +2.659615202676212d-1*t120*t21*t46*t8*t9
     +        +1.329807601338106d-1*t169*t19*t21*t8*t9
     +        -2.2348962373293998d-4*t19*t3*t42*t45*t8*t88
     +        +5.319230405352424d-1*t11*t48*t72*t8
     +        -2.659615202676212d-1*t14*t21*t46*t72*t8
     +        -8.939584949317599d-4*t21*t3*t51*t70*t8
     +        +4.4697924746587997d-4*t19*t3*t45*t46*t70*t8
     +        -4.4697924746587997d-4*t3*t42*t45*t46*t51*t8
     +        +8.939584949317599d-4*t118*t19*t21*t3*t8
     +        -2.659615202676212d-1*t116*t14*t19*t21*t8
     +        -5.319230405352424d-1*t11*t120*t14*t8
     +        +1.6844800327854395d-1*t15*t44*t47*t49*t50*t69*t71
     +        +3.5449077018110403d+0*t117*t17*t44*t50
     +        +3.5449077018110403d+0*t119*t17*t44*t47)
     +        +8.865384008920708d-2*t45*t7*t73*t8*t9*t90
     +        -4.469792474658801d-4*t22*t3*t42*t7*t8*t89*t90
     +        -1.7730768017841417d-1*t21*t7*t73*t8*t88*t9
     +        -3.5461536035682834d-1*t116*t21*t52*t7*t8*t9
     +        +1.7730768017841417d-1*t116*t22*t45*t46*t7*t8*t9
     +        -3.5461536035682834d-1*t121*t21*t46*t7*t8*t9
     +        -1.7730768017841417d-1*t169*t21*t22*t7*t8*t9
     +        +2.9798616497725333d-4*t22*t3*t42*t45*t7*t8*t88
     +        -7.092307207136567d-1*t11*t48*t7*t73*t8
     +        +3.5461536035682834d-1*t14*t21*t46*t7*t73*t8
     +        +1.1919446599090135d-3*t21*t3*t52*t7*t70*t8
     +        -5.959723299545067d-4*t22*t3*t45*t46*t7*t70*t8
     +        +5.959723299545067d-4*t3*t42*t45*t46*t52*t7*t8
     +        -1.1919446599090135d-3*t118*t21*t22*t3*t7*t8
     +        +3.5461536035682834d-1*t116*t14*t21*t22*t7*t8
     +        +7.092307207136567d-1*t11*t121*t14*t7*t8)+1.0d+0*t74*t87
     +        +1.0d+0*t23*(2.550431938149111d-4*gammaaa*t3*t86
     +        -4.8246324681544644d-7*t168*t3*t85-1.0457167652131799d-2
     +        *t3*t70)+2.0d+0*t115*t53+2.0d+0*t122*t43)*wght+Cmat3(iq,
     +        D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (1.0d+0*t102*t83+1.0d+0*t28*(
     +        -7.092307207136567d-1*t31*t60*t7*t8*t82
     +        +3.5461536035682834d-1*t33*t38*t58*t7*t8*t82
     +        +8.865384008920708d-2*t105*t29*t57*t7*t8*t82
     +        -1.7730768017841417d-1*t103*t29*t38*t7*t8*t82
     +        -3.5461536035682834d-1*t29*t31*t7*t8*
     +        (5.319230405352424d-1*t31*t60*t8*t81
     +        -2.659615202676212d-1*t33*t38*t58*t8*t81
     +        -6.649038006690531d-2*t105*t29*t57*t8*t81
     +        +1.329807601338106d-1*t103*t29*t38*t8*t81
     +        +2.659615202676212d-1*t29*t31*t8*(4.751830440958629d-2
     +        *t110*t15*t61*t78*t80+4.751830440958629d-2*t109*t15*t61
     +        *t78*t80-1.0083697378841376d-2*t34*t78*t79
     +        -1.1887928018410271d-4*t110*t12*t78*t79
     +        -1.1887928018410271d-4*t109*t12*t78*t79
     +        -3.361232459613792d-3*t58*t61*t62*t75*t78
     +        -1.680616229806896d-3*t103*t61*t75*t78
     +        -1.680616229806896d-3*t110*t34*t75*t78
     +        -1.680616229806896d-3*t109*t34*t75*t78
     +        +6.722464919227584d-3*t34*t55*t62*t78
     +        +4.7551712073641084d-4*t12*t127*t62*t78
     +        +6.722464919227584d-3*t55*t58*t61*t78
     +        -7.073553026306428d-2*t12*t128*t30*t33*t35*t62
     +        -7.073553026306428d-2*t12*t125*t33*t35*t62+2*t128*t35
     +        *t62+1.4147106052612857d-1*t12*t125*t36*t60
     +        +1.4147106052612857d-1*t12*t128*t30*t35*t60
     +        -7.073553026306428d-2*t12*t128*t33*t35*t58
     +        -3.536776513153214d-2*t12*t173*t33*t36
     +        -7.132756811046163d-4*t12*t174*t26*t36
     +        -3.536776513153214d-2*t12*t175*t30*t33*t35+t175*t35)
     +        +1.6844800327854395d-1*t15*t44*t59*t61*t62*t78*t80
     +        +1.6844800327854395d-1*t106*t15*t44*t61*t78*t80
     +        -8.939584949317599d-4*t26*t38*t63*t79*t8
     +        +4.4697924746587997d-4*t26*t37*t57*t58*t79*t8
     +        -4.4697924746587997d-4*t26*t55*t57*t58*t63*t8
     +        +2.659615202676212d-1*t125*t29*t38*t63*t8
     +        -1.329807601338106d-1*t125*t29*t37*t57*t58*t8
     +        +2.659615202676212d-1*t129*t29*t38*t58*t8
     +        -2.2348962373293998d-4*t103*t26*t37*t55*t57*t8
     +        +4.4697924746587997d-4*t111*t26*t38*t55*t8
     +        +3.3523443559941d-4*t104*t105*t26*t37*t55*t8
     +        -2.659615202676212d-1*t125*t33*t37*t38*t8
     +        +1.329807601338106d-1*t173*t29*t37*t38*t8
     +        +8.939584949317599d-4*t127*t26*t37*t38*t8
     +        -5.319230405352424d-1*t129*t31*t33*t8
     +        +3.5449077018110403d+0*t35*t44*(-2.369632590246221d-2
     +        *t105*t163*t20*t26*t76+9.478530360984885d-3*t103*t104
     +        *t20*t26*t76+1.895706072196977d-2*t104*t20*t26*t58*t75
     +        +5.639913617919763d+0*rhob*t104*t125*t20*t58
     +        -1.8799712059732546d+0*rhob*t173*t20*t57
     +        -3.759942411946509d+0*t125*t20*t57)
     +        +3.5449077018110403d+0*t126*t35*t44*t62
     +        +3.5449077018110403d+0*t128*t35*t44*t59)
     +        +1.1919446599090135d-3*t26*t38*t64*t7*t79*t8
     +        -5.959723299545067d-4*t26*t39*t57*t58*t7*t79*t8
     +        +5.959723299545067d-4*t26*t55*t57*t58*t64*t7*t8
     +        -3.5461536035682834d-1*t125*t29*t38*t64*t7*t8
     +        +1.7730768017841417d-1*t125*t29*t39*t57*t58*t7*t8
     +        -3.5461536035682834d-1*t130*t29*t38*t58*t7*t8
     +        +2.9798616497725333d-4*t103*t26*t39*t55*t57*t7*t8
     +        -4.469792474658801d-4*t104*t105*t26*t39*t55*t7*t8
     +        -5.959723299545067d-4*t112*t26*t38*t55*t7*t8
     +        +3.5461536035682834d-1*t125*t33*t38*t39*t7*t8
     +        -1.7730768017841417d-1*t173*t29*t38*t39*t7*t8
     +        -1.1919446599090135d-3*t127*t26*t38*t39*t7*t8
     +        +7.092307207136567d-1*t130*t31*t33*t7*t8)+1.0d+0*t40*(
     +        -1.0457167652131799d-2*t26*t79-4.8246324681544644d-7
     +        *t100*t172*t26+2.550431938149111d-4*gammabb*t101*t26)
     +        -3.361232459613792d-3*t113*t26*t75+2.0d+0*t124*t65
     +        +2.0d+0*t131*t56)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (1.0d+0*t5*(
     +        -3.5461536035682834d-1*t11*t7*t8*t9*
     +        (2.659615202676212d-1*t11*t8*t9*(-3.194401344146213d-4
     +        *t135*t15*t50*t84*t93+1.1297883647561385d-5*t135*t2*t46
     +        *t93+4.751830440958629d-2*t137*t15*t49*t50*t69*t71
     +        +9.503660881917259d-2*t119*t15*t49*t69*t71
     +        -1.1887928018410271d-4*t12*t137*t50*t69*t70
     +        -2.3775856036820542d-4*t119*t12*t69*t70
     +        -1.680616229806896d-3*t137*t16*t50*t66*t69
     +        -1.680616229806896d-3*t137*t46*t49*t66*t69
     +        -3.361232459613792d-3*t116*t49*t66*t69
     +        -3.361232459613792d-3*t119*t16*t66*t69
     +        +3.361232459613792d-3*t137*t16*t42*t69
     +        +5.648941823780692d-6*t135*t2*t49*t50
     +        -1.0681943420417696d-6*t12*t152*t18*t3
     +        -3.536776513153214d-2*t12*t14*t177*t18
     +        -3.536776513153214d-2*t10*t12*t14*t17*t179+t17*t179)
     +        -1.1323857927539448d-3*t135*t15*t44*t47*t84*t93
     +        +3.5449077018110403d+0*t17*t44*(2.838989964383816d-5*t20
     +        *t3*t45*t92+1.895706072196977d-2*t116*t20*t3*t67*t89
     +        -4.2584849465757246d-5*t20*t3*t46*t48*t89
     +        +3.1859543918777233d-5*t133*t2*t20*t89
     +        -7.964885979694308d-5*t133*t134*t154*t20*t46
     +        -1.8799712059732546d+0*rhoa*t177*t20*t45)
     +        +2.659615202676212d-1*t116*t21*t72*t8*t9
     +        +1.329807601338106d-1*t138*t21*t46*t8*t9
     +        +1.329807601338106d-1*t177*t19*t21*t8*t9
     +        +1.1268008665170462d-6*t132*t133*t19*t46*t8*t89
     +        +2.0081746054220725d-6*t19*t21*t3*t8*t86
     +        -8.939584949317599d-4*t21*t3*t70*t72*t8
     +        -4.4697924746587997d-4*t3*t42*t45*t46*t72*t8
     +        -7.512005776780308d-7*t132*t133*t45*t51*t8
     +        -2.0081746054220725d-6*t114*t21*t3*t51*t8
     +        +1.0040873027110363d-6*t114*t19*t3*t45*t46*t8
     +        -4.4697924746587997d-4*t116*t19*t3*t42*t45*t8
     +        +7.512005776780308d-7*t133*t178*t19*t45*t8
     +        +8.939584949317599d-4*t120*t21*t3*t42*t8
     +        -2.659615202676212d-1*t11*t138*t14*t8
     +        +1.6844800327854395d-1*t137*t15*t44*t47*t49*t69*t71
     +        +3.368960065570879d-1*t117*t15*t44*t49*t69*t71)
     +        -3.5461536035682834d-1*t116*t21*t7*t73*t8*t9
     +        -1.7730768017841417d-1*t139*t21*t46*t7*t8*t9
     +        -1.7730768017841417d-1*t177*t21*t22*t7*t8*t9
     +        -1.502401155356062d-6*t132*t133*t22*t46*t7*t8*t89
     +        -2.677566140562764d-6*t21*t22*t3*t7*t8*t86
     +        +1.1919446599090135d-3*t21*t3*t7*t70*t73*t8
     +        +5.959723299545067d-4*t3*t42*t45*t46*t7*t73*t8
     +        +1.0016007702373746d-6*t132*t133*t45*t52*t7*t8
     +        +2.677566140562764d-6*t114*t21*t3*t52*t7*t8
     +        -1.338783070281382d-6*t114*t22*t3*t45*t46*t7*t8
     +        +5.959723299545067d-4*t116*t22*t3*t42*t45*t7*t8
     +        -1.0016007702373746d-6*t133*t178*t22*t45*t7*t8
     +        -1.1919446599090135d-3*t121*t21*t3*t42*t7*t8
     +        +3.5461536035682834d-1*t11*t139*t14*t7*t8)
     +        +1.510124173904079d-5*t3*t53*t92+2.0d+0*t115*t74
     +        -6.722464919227584d-3*t122*t3*t66+1.0d+0*t140*t43+1.0d+0
     +        *t23*(1.8092371755579242d-7*gammaaa*t176*t3
     +        -6.040496695616315d-5*t114*t3))*wght+Cmat3(iq,
     +        D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (1.0d+0*t28*(
     +        -3.5461536035682834d-1*t29*t31*t7*t8*
     +        (2.659615202676212d-1*t29*t31*t8*(-3.194401344146213d-4
     +        *t108*t144*t15*t62*t99+4.751830440958629d-2*t146*t15*t61
     +        *t62*t78*t80+9.503660881917259d-2*t128*t15*t61*t78*t80
     +        -1.1887928018410271d-4*t12*t146*t62*t78*t79
     +        -2.3775856036820542d-4*t12*t128*t78*t79
     +        -1.680616229806896d-3*t146*t34*t62*t75*t78
     +        -1.680616229806896d-3*t146*t58*t61*t75*t78
     +        -3.361232459613792d-3*t125*t61*t75*t78
     +        -3.361232459613792d-3*t128*t34*t75*t78
     +        +3.361232459613792d-3*t146*t34*t55*t78
     +        +5.648941823780692d-6*t144*t25*t61*t62
     +        +1.1297883647561385d-5*t108*t144*t25*t58
     +        -3.536776513153214d-2*t12*t181*t33*t36
     +        -1.0681943420417696d-6*t12*t161*t26*t36
     +        -3.536776513153214d-2*t12*t183*t30*t33*t35+t183*t35)
     +        -1.1323857927539448d-3*t108*t144*t15*t44*t59*t99
     +        -8.939584949317599d-4*t26*t38*t79*t8*t81
     +        -4.4697924746587997d-4*t26*t55*t57*t58*t8*t81
     +        +2.659615202676212d-1*t125*t29*t38*t8*t81
     +        +1.6844800327854395d-1*t146*t15*t44*t59*t61*t78*t80
     +        +3.368960065570879d-1*t126*t15*t44*t61*t78*t80
     +        -7.512005776780308d-7*t141*t142*t57*t63*t8
     +        -2.0081746054220725d-6*t123*t26*t38*t63*t8
     +        +1.0040873027110363d-6*t123*t26*t37*t57*t58*t8
     +        +1.329807601338106d-1*t147*t29*t38*t58*t8
     +        +1.1268008665170462d-6*t104*t141*t142*t37*t58*t8
     +        -4.4697924746587997d-4*t125*t26*t37*t55*t57*t8
     +        +7.512005776780308d-7*t142*t182*t37*t57*t8
     +        +8.939584949317599d-4*t129*t26*t38*t55*t8
     +        +1.329807601338106d-1*t181*t29*t37*t38*t8
     +        +2.0081746054220725d-6*t101*t26*t37*t38*t8
     +        -2.659615202676212d-1*t147*t31*t33*t8
     +        +3.5449077018110403d+0*t35*t44*(1.895706072196977d-2
     +        *t104*t125*t20*t26*t76-4.2584849465757246d-5*t104*t20
     +        *t26*t58*t60-7.964885979694308d-5*t142*t143*t163*t20*t58
     +        +2.838989964383816d-5*t107*t20*t26*t57
     +        -1.8799712059732546d+0*rhob*t181*t20*t57
     +        +3.1859543918777233d-5*t104*t142*t20*t25))
     +        +1.1919446599090135d-3*t26*t38*t7*t79*t8*t82
     +        +5.959723299545067d-4*t26*t55*t57*t58*t7*t8*t82
     +        -3.5461536035682834d-1*t125*t29*t38*t7*t8*t82
     +        +1.0016007702373746d-6*t141*t142*t57*t64*t7*t8
     +        +2.677566140562764d-6*t123*t26*t38*t64*t7*t8
     +        -1.338783070281382d-6*t123*t26*t39*t57*t58*t7*t8
     +        -1.502401155356062d-6*t104*t141*t142*t39*t58*t7*t8
     +        -1.7730768017841417d-1*t148*t29*t38*t58*t7*t8
     +        +5.959723299545067d-4*t125*t26*t39*t55*t57*t7*t8
     +        -1.0016007702373746d-6*t142*t182*t39*t57*t7*t8
     +        -1.1919446599090135d-3*t130*t26*t38*t55*t7*t8
     +        -1.7730768017841417d-1*t181*t29*t38*t39*t7*t8
     +        -2.677566140562764d-6*t101*t26*t38*t39*t7*t8
     +        +3.5461536035682834d-1*t148*t31*t33*t7*t8)+2.0d+0*t124
     +        *t83-6.722464919227584d-3*t131*t26*t75
     +        +1.510124173904079d-5*t107*t26*t65+1.0d+0*t149*t56
     +        +1.0d+0*(1.8092371755579242d-7*gammabb*t180*t26
     +        -6.040496695616315d-5*t123*t26)*t40)*wght+Cmat3(iq,
     +        D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (1.0d+0*t5*(
     +        -3.5461536035682834d-1*t11*t7*t8*t9*
     +        (2.659615202676212d-1*t11*t8*t9*(-3.194401344146213d-4
     +        *t135*t15*t187*t84*t93-3.194401344146213d-4*t135*t137
     +        *t15*t84*t93+4.751830440958629d-2*t15*t189*t49*t69*t71
     +        +4.751830440958629d-2*t15*t188*t49*t69*t71
     +        -1.680616229806896d-3*t16*t189*t66*t69
     +        -1.680616229806896d-3*t16*t188*t66*t69
     +        +5.648941823780692d-6*t135*t137*t2*t49
     +        +2.5379391956869274d-8*t135*t170*t49
     +        -2.3995753309638476d-9*t12*t168*t18*t3
     +        +3.221137646093454d-6*t14*t15*t156*t186)
     +        +3.7874396480580436d-9*t114*t184*t19*t8*t89
     +        +1.1293914261386051d-4*t134*t135*t187*t20*t44*t89
     +        +1.1293914261386051d-4*t134*t135*t137*t20*t44*t89
     +        -2.2536017330340924d-6*t132*t133*t45*t72*t8
     +        -6.024523816266217d-6*t114*t21*t3*t72*t8
     +        +1.0124912502475182d-8*t133*t185*t19*t45*t8
     +        +1.3409377423976399d-3*t138*t21*t3*t42*t8
     +        +9.022265057551551d-9*t176*t19*t21*t3*t8
     +        -2.24003435190034d-2*t189*t20*t44*t45*t67*t69
     +        -2.24003435190034d-2*t188*t20*t44*t45*t67*t69
     +        -9.490367802866481d-7*t154*t186*t20*t44*t48)
     +        -5.049919530744059d-9*t114*t184*t22*t7*t8*t89
     +        +3.0048023107121236d-6*t132*t133*t45*t7*t73*t8
     +        +8.032698421688291d-6*t114*t21*t3*t7*t73*t8
     +        -1.3499883336633578d-8*t133*t185*t22*t45*t7*t8
     +        -1.7879169898635203d-3*t139*t21*t3*t42*t7*t8
     +        -1.202968674340207d-8*t176*t21*t22*t3*t7*t8)
     +        +4.5303725217122376d-5*t3*t74*t92-1.0083697378841378d-2
     +        *t140*t3*t66-6.784639408342217d-8*t23*t3/rhoa*
     +        *6.666666666666667d+0)*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (1.0d+0*t28*(
     +        -3.5461536035682834d-1*t29*t31*t7*t8*
     +        (2.659615202676212d-1*t29*t31*t8*(-3.194401344146213d-4
     +        *t108*t144*t15*t193*t99-3.194401344146213d-4*t108*t144
     +        *t146*t15*t99+4.751830440958629d-2*t15*t195*t61*t78*t80
     +        +4.751830440958629d-2*t15*t194*t61*t78*t80
     +        -1.680616229806896d-3*t195*t34*t75*t78
     +        -1.680616229806896d-3*t194*t34*t75*t78
     +        +5.648941823780692d-6*t144*t146*t25*t61
     +        +2.5379391956869274d-8*t144*t174*t61
     +        -2.3995753309638476d-9*t12*t172*t26*t36
     +        +3.221137646093454d-6*t15*t165*t192*t33)
     +        -2.2536017330340924d-6*t141*t142*t57*t8*t81
     +        -6.024523816266217d-6*t123*t26*t38*t8*t81
     +        +1.0124912502475182d-8*t142*t191*t37*t57*t8
     +        +1.3409377423976399d-3*t147*t26*t38*t55*t8
     +        +9.022265057551551d-9*t180*t26*t37*t38*t8
     +        +3.7874396480580436d-9*t104*t123*t190*t37*t8
     +        -2.24003435190034d-2*t195*t20*t44*t57*t76*t78
     +        -2.24003435190034d-2*t194*t20*t44*t57*t76*t78
     +        -9.490367802866481d-7*t163*t192*t20*t44*t60
     +        +1.1293914261386051d-4*t104*t143*t144*t193*t20*t44
     +        +1.1293914261386051d-4*t104*t143*t144*t146*t20*t44)
     +        +3.0048023107121236d-6*t141*t142*t57*t7*t8*t82
     +        +8.032698421688291d-6*t123*t26*t38*t7*t8*t82
     +        -1.3499883336633578d-8*t142*t191*t39*t57*t7*t8
     +        -1.7879169898635203d-3*t148*t26*t38*t55*t7*t8
     +        -1.202968674340207d-8*t180*t26*t38*t39*t7*t8
     +        -5.049919530744059d-9*t104*t123*t190*t39*t7*t8)
     +        +4.5303725217122376d-5*t107*t26*t83
     +        -1.0083697378841378d-2*t149*t26*t75-6.784639408342217d-8
     +        *t26*t40/rhob**6.666666666666667d+0)*wght+Cmat3(iq,
     +        D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammaaa*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhoa
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhoa**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammaaa*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammaaa*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhoa*t19*t27*t28
            t30 = 1/rhoa**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhoa*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammaaa*t2)
            t40 = 1/rhoa**3.3333333333333337d+0
            t41 = rhoa**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammaaa**2
            t47 = 1/rhoa**6
            t48 = 1.0738660792206783d-4*t3*t46*t47
     +        -3.3250786312207836d-1*t4*t45-4.135669939329333d-1*t45
     +        -8.963286558970114d-3*gammaaa*t3*t40
            t49 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammaaa*t3*t40
            t50 = 1/t10**5
            t51 = t28**2
            t52 = 2.8199568089598814d+0*rhoa*t19*t50*t51
     +        -1.8799712059732546d+0*rhoa*t19*t27*t49
     +        -3.759942411946509d+0*t19*t27*t28
            t53 = 1/rhoa**4
            t54 = 1/t9**3
            t55 = -2.827433388230826d+1*t12*t14*t51*t54
     +        +1.413716694115413d+1*t12*t14*t31*t49
     +        +5.654866776461652d+1*rhoa*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t56 = t32**2
            t57 = -3.536776513153214d-2*t11*t13*t16*t56*t9
     +        -3.536776513153214d-2*t11*t13*t16*t55*t9
     +        -2.1220659078919285d-1*t11*t17*t53*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t56+t16
     +        *t55-3.536776513153214d-2*t11*t13*t17*t49
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30
            t58 = 2.659615202676212d-1*t10*t57*t7*t8
     +        -6.649038006690531d-2*t18*t27*t51*t7*t8
     +        +1.329807601338106d-1*t18*t20*t49*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*t52+3.5449077018110403d+0
     +        *t16*t26*t29*t32
            t59 = -3.5461536035682834d-1*t10*t58*t6*t7*t8
     +        +8.865384008920708d-2*t21*t27*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t49*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7
            t60 = 1/rhoa**5
            t61 = 4.4816432794850564d-3*t24*t3-4.026997797077544d-5
     +        *gammaaa*t3*t60
            t62 = 4.026997797077544d-5*gammaaa*t3*t60
     +        -4.4816432794850564d-3*t24*t3
            t63 = -1.8799712059732546d+0*rhoa*t19*t27*t62
     +        +9.478530360984885d-3*t19*t28*t3*t37*t50
     +        -6.31902024065659d-3*t19*t27*t3*t36
            t64 = 1/rhoa**4.333333333333333d+0
            t65 = 1.413716694115413d+1*t12*t14*t31*t62
     +        -9.503660881917259d-2*t14*t28*t3*t41*t54
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t66 = -3.536776513153214d-2*t11*t13*t16*t65*t9+t16*t65
     +        +2.3775856036820542d-4*t11*t17*t3*t64
     +        -3.536776513153214d-2*t11*t13*t17*t62
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39
            t67 = 2.659615202676212d-1*t10*t66*t7*t8
     +        +1.329807601338106d-1*t18*t20*t62*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*t63+1.6844800327854395d-1
     +        *t14*t26*t29*t31*t39*t41
            t68 = -3.5461536035682834d-1*t10*t6*t67*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t6*t62*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7
            t69 = 1/rhoa**3.6666666666666664d+0
            t70 = exp(-8.98553844191778d-3*gammaaa*t2)
            t71 = 1/rhoa**1.6666666666666669d+0
            t72 = exp(t38-8.98553844191778d-3*gammaaa*t2)
            t73 = 4.751830440958629d-2*t14*t3*t31*t41
            t74 = t73-4.49276922095889d-3*t2
            t75 = 4.751830440958629d-2*t14*t31*t39*t41*t74
     +        -1.680616229806896d-3*t15*t36*t39*t74
     +        -3.194401344146213d-4*t14*t45*t54*t72
     +        +5.340971710208847d-7*t11*t17*t3*t47
            t76 = 2.659615202676212d-1*t10*t7*t75*t8
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t74
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72
     +        -7.512005776780308d-7*t18*t27*t69*t7*t70
     +        -2.0081746054220725d-6*t18*t20*t3*t60*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
            t77 = -3.5461536035682834d-1*t10*t6*t7*t76*t8
     +        +1.0016007702373746d-6*t21*t27*t6*t69*t7*t70
     +        +2.677566140562764d-6*t20*t21*t3*t6*t60*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7
            t78 = gammaaa**3
            t79 = 1/rhoa**9.666666666666666d+0
            t80 = 1/rhoa**7
            t81 = 7.517062554544749d-4*t3*t46*t80-1.286568658174524d-6
     +        *t3*t78*t79-2.2167190874805223d-1*t4*t71
     +        -2.757113292886222d-1*t71-3.386130477833154d-2*gammaaa
     +        *t3*t64
            t82 = 1/t10**7
            t83 = t28**3
            t84 = 1/t9**4
            t85 = 8.482300164692478d+1*t12*t14*t83*t84
     +        +1.413716694115413d+1*t12*t14*t31*t81
     +        -1.6964600329384957d+2*rhoa*t14*t51*t54
     +        -8.482300164692478d+1*t12*t14*t28*t49*t54
     +        +8.482300164692478d+1*rhoa*t14*t31*t49
     +        +8.482300164692478d+1*t14*t28*t31
            t86 = t32**3
            t87 = 1/rhoa**8.666666666666666d+0
            t88 = 4.8246324681544644d-7*t3*t46*t87
     +        -2.550431938149111d-4*gammaaa*t3*t47
     +        +1.0457167652131799d-2*t3*t40
            t89 = 1/rhoa**5.333333333333333d+0
            t90 = 1.413716694115413d+1*t12*t14*t31*t88
     +        +2.8510982645751776d-1*t14*t3*t41*t51*t84
     +        -5.654866776461652d+1*t12*t14*t28*t54*t62
     +        +5.654866776461652d+1*rhoa*t14*t31*t62
     +        -9.503660881917259d-2*t14*t3*t41*t49*t54
     +        -3.8014643527669034d-1*t14*t28*t3*t37*t54
     +        +9.503660881917259d-2*t14*t3*t31*t36
            t91 = 1/rhoa**7.666666666666667d+0
            t92 = 6.040496695616315d-5*t3*t60-1.8092371755579242d-7
     +        *gammaaa*t3*t91
            t93 = 1/rhoa**4.666666666666667d+0
            t94 = 1.413716694115413d+1*t12*t14*t31*t92
     +        +9.583204032438639d-4*t14*t28*t45*t70*t84
     +        -6.388802688292426d-4*t14*t54*t70*t71
     +        -1.9007321763834517d-1*t14*t3*t41*t54*t62
     +        +4.269775509670888d-4*t13*t14*t28*t3*t54
     +        -4.269775509670888d-4*t14*t3*t30*t31
            t95 = exp(-1.3478307662876668d-2*gammaaa*t2)
            t96 = 1/rhoa**6.333333333333333d+0
            t97 = exp(t38-1.3478307662876668d-2*gammaaa*t2)
            t98 = t73-8.98553844191778d-3*t2
            t99 = -3.194401344146213d-4*t14*t45*t54*t70
     +        -2.134887754835444d-4*t13*t14*t3*t31
            t100 = t74**2
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t35*t5*wght+1.0d+0*t22*t25*wght
     +        +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t44*t5*wght-3.361232459613792d-3
     +        *t22*t3*t36*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = 1.0d+0*t5*t59*wght+1.0d+0*t22*t48
     +        *wght+2.0d+0*t25*t35*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*t5*t68*wght+1.0d+0*t22*t61
     +        *wght+1.0d+0*t25*t44*wght-3.361232459613792d-3*t3*t35
     +        *t36*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t5*t77*wght
     +        +1.510124173904079d-5*t22*t3*t53*wght
     +        -6.722464919227584d-3*t3*t36*t44*wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t86*t9
     +        -3.536776513153214d-2*t11*t13*t16*t85*t9
     +        +8.488263631567714d-1*t11*t17*t60*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t56*t9
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t55*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t55*t9
     +        -6.366197723675786d-1*t11*t16*t32*t53*t9+t16*t86+t16*t85
     +        -3.536776513153214d-2*t11*t13*t17*t81
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t56+3*t16*t32*t55
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t55
     +        -6.366197723675786d-1*t11*t17*t28*t53
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t49
     +        +2.1220659078919285d-1*t11*t17*t30*t49
     +        +4.244131815783857d-1*t11*t16*t28*t30*t32)
     +        +3.5449077018110403d+0*t16*t26*(-7.049892022399704d+0
     +        *rhoa*t19*t82*t83-1.8799712059732546d+0*rhoa*t19*t27*t81
     +        +8.459870426879645d+0*t19*t50*t51+8.459870426879645d+0
     +        *rhoa*t19*t28*t49*t50-5.639913617919763d+0*t19*t27*t49)
     +        +9.973557010035795d-2*t18*t50*t7*t8*t83
     +        +1.329807601338106d-1*t18*t20*t7*t8*t81
     +        +3.989422804014318d-1*t20*t28*t57*t7*t8
     +        -1.994711402007159d-1*t27*t33*t51*t7*t8
     +        +3.989422804014318d-1*t20*t33*t49*t7*t8
     +        -1.994711402007159d-1*t18*t27*t28*t49*t7*t8
     +        -7.978845608028636d-1*t10*t13*t57*t7
     +        -1.5957691216057274d+0*t10*t18*t53*t7
     +        +1.994711402007159d-1*t13*t18*t27*t51*t7
     +        -3.989422804014318d-1*t13*t18*t20*t49*t7
     +        +1.5957691216057274d+0*t10*t30*t33*t7
     +        -7.978845608028636d-1*t13*t20*t28*t33*t7
     +        +7.978845608028636d-1*t18*t20*t28*t30*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t56
     +        +3.5449077018110403d+0*t16*t26*t29*t55
     +        +7.089815403622081d+0*t16*t26*t32*t52)
     +        -1.3298076013381063d-1*t21*t50*t6*t7*t8*t83
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t81
     +        -5.319230405352425d-1*t20*t28*t58*t6*t7*t8
     +        +2.6596152026762127d-1*t27*t34*t51*t6*t7*t8
     +        -5.319230405352425d-1*t20*t34*t49*t6*t7*t8
     +        +2.6596152026762127d-1*t21*t27*t28*t49*t6*t7*t8
     +        +1.063846081070485d+0*t10*t13*t58*t6*t7
     +        +2.12769216214097d+0*t10*t21*t53*t6*t7
     +        -2.6596152026762127d-1*t13*t21*t27*t51*t6*t7
     +        +5.319230405352425d-1*t13*t20*t21*t49*t6*t7
     +        -2.12769216214097d+0*t10*t30*t34*t6*t7
     +        +1.063846081070485d+0*t13*t20*t28*t34*t6*t7
     +        -1.063846081070485d+0*t20*t21*t28*t30*t6*t7)*wght+1.0d+0
     +        *t22*(-7.517062554544749d-4*t3*t46*t80
     +        +1.286568658174524d-6*t3*t78*t79+2.2167190874805223d-1
     +        *t4*t71+2.757113292886222d-1*t71+3.386130477833154d-2
     +        *gammaaa*t3*t64)*wght+3.0d+0*t25*t59*wght+3.0d+0*t35*t48
     +        *wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t9*t90+t16*t90
     +        -7.073553026306428d-2*t11*t13*t16*t32*t65*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t65*t9
     +        -7.132756811046163d-4*t11*t17*t3*t89
     +        -3.536776513153214d-2*t11*t13*t17*t88+2*t16*t32*t65
     +        -7.073553026306428d-2*t11*t13*t16*t28*t65
     +        +4.7551712073641084d-4*t11*t32*t39*t64
     +        -7.073553026306428d-2*t11*t13*t16*t32*t62
     +        +1.4147106052612857d-1*t11*t17*t30*t62
     +        +4.751830440958629d-2*t14*t31*t39*t41*t56
     +        -1.1887928018410271d-4*t11*t39*t40*t56
     +        -1.680616229806896d-3*t15*t36*t39*t56
     +        +4.751830440958629d-2*t14*t31*t39*t41*t55
     +        -1.1887928018410271d-4*t11*t39*t40*t55
     +        -1.680616229806896d-3*t15*t36*t39*t55
     +        -1.680616229806896d-3*t31*t36*t39*t49
     +        -1.0083697378841376d-2*t15*t39*t40-3.361232459613792d-3
     +        *t28*t31*t32*t36*t39+6.722464919227584d-3*t15*t24*t32
     +        *t39+6.722464919227584d-3*t24*t28*t31*t39)
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t88-2.369632590246221d-2*t19*t3*t37*t51
     +        *t82+5.639913617919763d+0*rhoa*t19*t28*t50*t62
     +        -3.759942411946509d+0*t19*t27*t62+9.478530360984885d-3
     +        *t19*t3*t37*t49*t50+1.895706072196977d-2*t19*t28*t3*t36
     +        *t50)+1.329807601338106d-1*t18*t20*t7*t8*t88
     +        +2.659615202676212d-1*t20*t28*t66*t7*t8
     +        +2.659615202676212d-1*t20*t33*t62*t7*t8
     +        -1.329807601338106d-1*t18*t27*t28*t62*t7*t8
     +        -6.649038006690531d-2*t27*t42*t51*t7*t8
     +        +1.329807601338106d-1*t20*t42*t49*t7*t8
     +        -5.319230405352424d-1*t10*t13*t66*t7
     +        +8.939584949317599d-4*t18*t20*t3*t64*t7
     +        -2.659615202676212d-1*t13*t18*t20*t62*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t57*t7
     +        +3.3523443559941d-4*t18*t24*t3*t50*t51*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t3*t49*t7
     +        +5.319230405352424d-1*t10*t30*t42*t7
     +        -2.659615202676212d-1*t13*t20*t28*t42*t7
     +        -8.939584949317599d-4*t20*t3*t33*t40*t7
     +        +4.4697924746587997d-4*t18*t27*t28*t3*t40*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t33*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t65
     +        +3.5449077018110403d+0*t16*t26*t32*t63
     +        +1.6844800327854395d-1*t14*t26*t31*t39*t41*t52
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t32*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t88
     +        -3.5461536035682834d-1*t20*t28*t6*t67*t7*t8
     +        -3.5461536035682834d-1*t20*t34*t6*t62*t7*t8
     +        +1.7730768017841417d-1*t21*t27*t28*t6*t62*t7*t8
     +        +8.865384008920708d-2*t27*t43*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t43*t49*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t6*t67*t7
     +        -1.1919446599090135d-3*t20*t21*t3*t6*t64*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t6*t62*t7
     +        -5.959723299545067d-4*t20*t24*t3*t58*t6*t7
     +        -4.469792474658801d-4*t21*t24*t3*t50*t51*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t3*t49*t6*t7
     +        -7.092307207136567d-1*t10*t30*t43*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t28*t43*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t34*t40*t6*t7
     +        -5.959723299545067d-4*t21*t27*t28*t3*t40*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t34*t6*t7)*wght
     +        +1.0d+0*t22*(-4.8246324681544644d-7*t3*t46*t87
     +        +2.550431938149111d-4*gammaaa*t3*t47
     +        -1.0457167652131799d-2*t3*t40)*wght+2.0d+0*t25*t68*wght
     +        +2.0d+0*t35*t61*wght-3.361232459613792d-3*t3*t36*t59
     +        *wght+1.0d+0*t44*t48*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(-3.536776513153214d-2
     +        *t11*t13*t16*t9*t94+t16*t94-3.536776513153214d-2*t11*t13
     +        *t17*t92-1.0681943420417696d-6*t11*t17*t3*t80
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41*t74
     +        -1.1887928018410271d-4*t11*t32*t39*t40*t74
     +        -1.680616229806896d-3*t15*t32*t36*t39*t74
     +        -1.680616229806896d-3*t28*t31*t36*t39*t74
     +        +3.361232459613792d-3*t15*t24*t39*t74
     +        -3.194401344146213d-4*t14*t32*t45*t54*t72
     +        +1.1297883647561385d-5*t2*t28*t54*t72
     +        +5.648941823780692d-6*t2*t31*t32*t72
     +        +9.503660881917259d-2*t14*t31*t39*t41*t65
     +        -2.3775856036820542d-4*t11*t39*t40*t65
     +        -3.361232459613792d-3*t15*t36*t39*t65
     +        -3.361232459613792d-3*t31*t36*t39*t62)
     +        +7.512005776780308d-7*t18*t27*t7*t70*t93
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhoa*t19*t27*t92-7.964885979694308d-5*t19*t28*t70*t71
     +        *t82+3.1859543918777233d-5*t19*t2*t50*t70
     +        +1.895706072196977d-2*t19*t3*t37*t50*t62
     +        +2.838989964383816d-5*t19*t27*t3*t53
     +        -4.2584849465757246d-5*t19*t28*t3*t30*t50)
     +        +1.329807601338106d-1*t18*t20*t7*t8*t92
     +        +1.329807601338106d-1*t20*t28*t7*t75*t8
     +        +2.659615202676212d-1*t20*t42*t62*t7*t8
     +        -2.659615202676212d-1*t10*t13*t7*t75
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41*t74
     +        -1.1323857927539448d-3*t14*t26*t29*t45*t54*t72
     +        +1.1268008665170462d-6*t18*t28*t50*t69*t7*t70
     +        -7.512005776780308d-7*t27*t33*t69*t7*t70
     +        +8.939584949317599d-4*t20*t24*t3*t66*t7
     +        -4.4697924746587997d-4*t18*t24*t27*t3*t62*t7
     +        -2.0081746054220725d-6*t20*t3*t33*t60*t7
     +        +1.0040873027110363d-6*t18*t27*t28*t3*t60*t7
     +        +2.0081746054220725d-6*t18*t20*t3*t47*t7
     +        -8.939584949317599d-4*t20*t3*t40*t42*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t42*t7
     +        +3.368960065570879d-1*t14*t26*t31*t39*t41*t63)
     +        -1.0016007702373746d-6*t21*t27*t6*t7*t70*t93
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t92
     +        -1.7730768017841417d-1*t20*t28*t6*t7*t76*t8
     +        -3.5461536035682834d-1*t20*t43*t6*t62*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t6*t7*t76
     +        -1.502401155356062d-6*t21*t28*t50*t6*t69*t7*t70
     +        +1.0016007702373746d-6*t27*t34*t6*t69*t7*t70
     +        -1.1919446599090135d-3*t20*t24*t3*t6*t67*t7
     +        +5.959723299545067d-4*t21*t24*t27*t3*t6*t62*t7
     +        +2.677566140562764d-6*t20*t3*t34*t6*t60*t7
     +        -1.338783070281382d-6*t21*t27*t28*t3*t6*t60*t7
     +        -2.677566140562764d-6*t20*t21*t3*t47*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t40*t43*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t43*t6*t7)*wght
     +        +1.0d+0*t22*(1.8092371755579242d-7*gammaaa*t3*t91
     +        -6.040496695616315d-5*t3*t60)*wght+1.0d+0*t25*t77*wght
     +        -6.722464919227584d-3*t3*t36*t68*wght+2.0d+0*t44*t61
     +        *wght+1.510124173904079d-5*t3*t35*t53*wght+Cmat3(iq,
     +        D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(4.751830440958629d-2
     +        *t14*t31*t39*t41*t99-1.680616229806896d-3*t15*t36*t39
     +        *t99-3.194401344146213d-4*t14*t45*t54*t72*t98
     +        +3.221137646093454d-6*t13*t14*t84*t97
     +        +2.5379391956869274d-8*t31*t72*t89-2.3995753309638476d-9
     +        *t11*t17*t3*t87-3.194401344146213d-4*t14*t45*t54*t72*t74
     +        +5.648941823780692d-6*t2*t31*t72*t74
     +        +4.751830440958629d-2*t100*t14*t31*t39*t41
     +        -1.680616229806896d-3*t100*t15*t36*t39)
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t99
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t98
     +        -9.490367802866481d-7*t19*t26*t30*t82*t97
     +        +1.0124912502475182d-8*t18*t27*t7*t70*t96
     +        +3.7874396480580436d-9*t18*t50*t60*t7*t95
     +        +9.022265057551551d-9*t18*t20*t3*t7*t91
     +        +1.3409377423976399d-3*t20*t24*t3*t7*t75
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t74
     +        -2.2536017330340924d-6*t27*t42*t69*t7*t70
     +        -6.024523816266217d-6*t20*t3*t42*t60*t7
     +        -2.24003435190034d-2*t100*t19*t26*t27*t37*t39)
     +        -1.3499883336633578d-8*t21*t27*t6*t7*t70*t96
     +        -5.049919530744059d-9*t21*t50*t6*t60*t7*t95
     +        -1.202968674340207d-8*t20*t21*t3*t6*t7*t91
     +        -1.7879169898635203d-3*t20*t24*t3*t6*t7*t76
     +        +3.0048023107121236d-6*t27*t43*t6*t69*t7*t70
     +        +8.032698421688291d-6*t20*t3*t43*t6*t60*t7)*wght
     +        -1.0083697378841378d-2*t3*t36*t77*wght
     +        +4.5303725217122376d-5*t3*t44*t53*wght
     +        -6.784639408342217d-8*t22*t3*wght/rhoa*
     +        *6.666666666666667d+0+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = exp(-4.49276922095889d-3*gammabb*t2)
            t4 = 1.0d+0-t3
            t5 = -7.481426920246764d-1*t1*t4-9.305257363491d-1*t1
            t6 = param(2)
            t7 = param(3)
            t8 = 1/rhob
            t9 = 7.481426920246764d-1*t1*t4+9.305257363491d-1*t1
            t10 = t9**5.0d-1
            t11 = t7**2
            t12 = rhob**2
            t13 = 1/t12
            t14 = 1/t11
            t15 = 1/t9
            t16 = exp(-1.413716694115413d+1*t12*t14*t15)
            t17 = t16-1.0d+0
            t18 = -3.536776513153214d-2*t11*t13*t17*t9+t16-1.5d+0
            t19 = 1/t7
            t20 = 1/t10
            t21 = 2.659615202676212d-1*t10*t18*t7*t8
     +        +1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t19
     +        *t20)
            t22 = -3.5461536035682834d-1*t10*t21*t6*t7*t8-param(1)
     +        +1.0d+0
            t23 = rhob**3.333333333333333d-1
            t24 = 1/rhob**2.3333333333333334d+0
            t25 = -9.975235893662352d-1*t23*t4+8.963286558970113d-3
     +        *gammabb*t24*t3-1.2407009817987999d+0*t23
            t26 = 5.641895835477563d-1
            t27 = 1/t10**3
            t28 = 9.975235893662352d-1*t23*t4-8.963286558970113d-3
     +        *gammabb*t24*t3+1.2407009817987999d+0*t23
            t29 = 3.759942411946509d+0*t19*t20-1.8799712059732546d+0
     +        *rhob*t19*t27*t28
            t30 = 1/rhob**3
            t31 = 1/t9**2
            t32 = 1.413716694115413d+1*t12*t14*t28*t31
     +        -2.827433388230826d+1*rhob*t14*t15
            t33 = -3.536776513153214d-2*t11*t13*t16*t32*t9
     +        +7.073553026306428d-2*t11*t17*t30*t9+t16*t32
     +        -3.536776513153214d-2*t11*t13*t17*t28
            t34 = 2.659615202676212d-1*t10*t33*t7*t8
     +        +1.329807601338106d-1*t18*t20*t28*t7*t8
     +        -2.659615202676212d-1*t10*t13*t18*t7
     +        +3.5449077018110403d+0*t16*t26*t29
            t35 = -3.5461536035682834d-1*t10*t34*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t28*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t21*t6*t7
            t36 = 1/t1
            t37 = 1/t23
            t38 = -1.413716694115413d+1*t12*t14*t15
            t39 = exp(t38-4.49276922095889d-3*gammabb*t2)
            t40 = 1/rhob**3.3333333333333337d+0
            t41 = rhob**6.666666666666666d-1
            t42 = 4.751830440958629d-2*t14*t31*t39*t41
     +        -1.1887928018410271d-4*t11*t17*t3*t40
     +        -1.680616229806896d-3*t15*t36*t39
            t43 = 2.659615202676212d-1*t10*t42*t7*t8
     +        +4.4697924746587997d-4*t18*t20*t24*t3*t7
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39
            t44 = -3.5461536035682834d-1*t10*t43*t6*t7*t8
     +        -5.959723299545067d-4*t20*t21*t24*t3*t6*t7
            t45 = 1/t41
            t46 = gammabb**2
            t47 = 1/rhob**6
            t48 = 1.0738660792206783d-4*t3*t46*t47
     +        -3.3250786312207836d-1*t4*t45-4.135669939329333d-1*t45
     +        -8.963286558970114d-3*gammabb*t3*t40
            t49 = -1.0738660792206783d-4*t3*t46*t47
     +        +3.3250786312207836d-1*t4*t45+4.135669939329333d-1*t45
     +        +8.963286558970114d-3*gammabb*t3*t40
            t50 = 1/t10**5
            t51 = t28**2
            t52 = 2.8199568089598814d+0*rhob*t19*t50*t51
     +        -1.8799712059732546d+0*rhob*t19*t27*t49
     +        -3.759942411946509d+0*t19*t27*t28
            t53 = 1/rhob**4
            t54 = 1/t9**3
            t55 = -2.827433388230826d+1*t12*t14*t51*t54
     +        +1.413716694115413d+1*t12*t14*t31*t49
     +        +5.654866776461652d+1*rhob*t14*t28*t31
     +        -2.827433388230826d+1*t14*t15
            t56 = t32**2
            t57 = -3.536776513153214d-2*t11*t13*t16*t56*t9
     +        -3.536776513153214d-2*t11*t13*t16*t55*t9
     +        -2.1220659078919285d-1*t11*t17*t53*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t32*t9+t16*t56+t16
     +        *t55-3.536776513153214d-2*t11*t13*t17*t49
     +        -7.073553026306428d-2*t11*t13*t16*t28*t32
     +        +1.4147106052612857d-1*t11*t17*t28*t30
            t58 = 2.659615202676212d-1*t10*t57*t7*t8
     +        -6.649038006690531d-2*t18*t27*t51*t7*t8
     +        +1.329807601338106d-1*t18*t20*t49*t7*t8
     +        +2.659615202676212d-1*t20*t28*t33*t7*t8
     +        -5.319230405352424d-1*t10*t13*t33*t7
     +        +5.319230405352424d-1*t10*t18*t30*t7
     +        -2.659615202676212d-1*t13*t18*t20*t28*t7
     +        +3.5449077018110403d+0*t16*t26*t52+3.5449077018110403d+0
     +        *t16*t26*t29*t32
            t59 = -3.5461536035682834d-1*t10*t58*t6*t7*t8
     +        +8.865384008920708d-2*t21*t27*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t49*t6*t7*t8
     +        -3.5461536035682834d-1*t20*t28*t34*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t34*t6*t7
     +        -7.092307207136567d-1*t10*t21*t30*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t28*t6*t7
            t60 = 1/rhob**5
            t61 = 4.4816432794850564d-3*t24*t3-4.026997797077544d-5
     +        *gammabb*t3*t60
            t62 = 4.026997797077544d-5*gammabb*t3*t60
     +        -4.4816432794850564d-3*t24*t3
            t63 = -1.8799712059732546d+0*rhob*t19*t27*t62
     +        +9.478530360984885d-3*t19*t28*t3*t37*t50
     +        -6.31902024065659d-3*t19*t27*t3*t36
            t64 = 1/rhob**4.333333333333333d+0
            t65 = 1.413716694115413d+1*t12*t14*t31*t62
     +        -9.503660881917259d-2*t14*t28*t3*t41*t54
     +        +9.503660881917259d-2*t14*t3*t31*t37
            t66 = -3.536776513153214d-2*t11*t13*t16*t65*t9+t16*t65
     +        +2.3775856036820542d-4*t11*t17*t3*t64
     +        -3.536776513153214d-2*t11*t13*t17*t62
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41
     +        -1.1887928018410271d-4*t11*t32*t39*t40
     +        -1.680616229806896d-3*t15*t32*t36*t39
     +        -1.680616229806896d-3*t28*t31*t36*t39
     +        +3.361232459613792d-3*t15*t24*t39
            t67 = 2.659615202676212d-1*t10*t66*t7*t8
     +        +1.329807601338106d-1*t18*t20*t62*t7*t8
     +        +1.329807601338106d-1*t20*t28*t42*t7*t8
     +        -2.659615202676212d-1*t10*t13*t42*t7
     +        -4.4697924746587997d-4*t18*t20*t3*t40*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t33*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t28*t3*t7
     +        +3.5449077018110403d+0*t16*t26*t63+1.6844800327854395d-1
     +        *t14*t26*t29*t31*t39*t41
            t68 = -3.5461536035682834d-1*t10*t6*t67*t7*t8
     +        -1.7730768017841417d-1*t20*t21*t6*t62*t7*t8
     +        -1.7730768017841417d-1*t20*t28*t43*t6*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t43*t6*t7
     +        +5.959723299545067d-4*t20*t21*t3*t40*t6*t7
     +        -5.959723299545067d-4*t20*t24*t3*t34*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t28*t3*t6*t7
            t69 = 1/rhob**3.6666666666666664d+0
            t70 = exp(-8.98553844191778d-3*gammabb*t2)
            t71 = 1/rhob**1.6666666666666669d+0
            t72 = exp(t38-8.98553844191778d-3*gammabb*t2)
            t73 = 4.751830440958629d-2*t14*t3*t31*t41
            t74 = t73-4.49276922095889d-3*t2
            t75 = 4.751830440958629d-2*t14*t31*t39*t41*t74
     +        -1.680616229806896d-3*t15*t36*t39*t74
     +        -3.194401344146213d-4*t14*t45*t54*t72
     +        +5.340971710208847d-7*t11*t17*t3*t47
            t76 = 2.659615202676212d-1*t10*t7*t75*t8
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t74
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72
     +        -7.512005776780308d-7*t18*t27*t69*t7*t70
     +        -2.0081746054220725d-6*t18*t20*t3*t60*t7
     +        +8.939584949317599d-4*t20*t24*t3*t42*t7
            t77 = -3.5461536035682834d-1*t10*t6*t7*t76*t8
     +        +1.0016007702373746d-6*t21*t27*t6*t69*t7*t70
     +        +2.677566140562764d-6*t20*t21*t3*t6*t60*t7
     +        -1.1919446599090135d-3*t20*t24*t3*t43*t6*t7
            t78 = gammabb**3
            t79 = 1/rhob**9.666666666666666d+0
            t80 = 1/rhob**7
            t81 = 7.517062554544749d-4*t3*t46*t80-1.286568658174524d-6
     +        *t3*t78*t79-2.2167190874805223d-1*t4*t71
     +        -2.757113292886222d-1*t71-3.386130477833154d-2*gammabb
     +        *t3*t64
            t82 = 1/t10**7
            t83 = t28**3
            t84 = 1/t9**4
            t85 = 8.482300164692478d+1*t12*t14*t83*t84
     +        +1.413716694115413d+1*t12*t14*t31*t81
     +        -1.6964600329384957d+2*rhob*t14*t51*t54
     +        -8.482300164692478d+1*t12*t14*t28*t49*t54
     +        +8.482300164692478d+1*rhob*t14*t31*t49
     +        +8.482300164692478d+1*t14*t28*t31
            t86 = t32**3
            t87 = 1/rhob**8.666666666666666d+0
            t88 = 4.8246324681544644d-7*t3*t46*t87
     +        -2.550431938149111d-4*gammabb*t3*t47
     +        +1.0457167652131799d-2*t3*t40
            t89 = 1/rhob**5.333333333333333d+0
            t90 = 1.413716694115413d+1*t12*t14*t31*t88
     +        +2.8510982645751776d-1*t14*t3*t41*t51*t84
     +        -5.654866776461652d+1*t12*t14*t28*t54*t62
     +        +5.654866776461652d+1*rhob*t14*t31*t62
     +        -9.503660881917259d-2*t14*t3*t41*t49*t54
     +        -3.8014643527669034d-1*t14*t28*t3*t37*t54
     +        +9.503660881917259d-2*t14*t3*t31*t36
            t91 = 1/rhob**7.666666666666667d+0
            t92 = 6.040496695616315d-5*t3*t60-1.8092371755579242d-7
     +        *gammabb*t3*t91
            t93 = 1/rhob**4.666666666666667d+0
            t94 = 1.413716694115413d+1*t12*t14*t31*t92
     +        +9.583204032438639d-4*t14*t28*t45*t70*t84
     +        -6.388802688292426d-4*t14*t54*t70*t71
     +        -1.9007321763834517d-1*t14*t3*t41*t54*t62
     +        +4.269775509670888d-4*t13*t14*t28*t3*t54
     +        -4.269775509670888d-4*t14*t3*t30*t31
            t95 = exp(-1.3478307662876668d-2*gammabb*t2)
            t96 = 1/rhob**6.333333333333333d+0
            t97 = exp(t38-1.3478307662876668d-2*gammabb*t2)
            t98 = t73-8.98553844191778d-3*t2
            t99 = -3.194401344146213d-4*t14*t45*t54*t70
     +        -2.134887754835444d-4*t13*t14*t3*t31
            t100 = t74**2
            fnc(iq) = 1.0d+0*t22*t5*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t35*t5*wght+1.0d+0*t22*t25*wght
     +        +Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t44*t5*wght-3.361232459613792d-3
     +        *t22*t3*t36*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*t5*t59*wght+1.0d+0*t22*t48
     +        *wght+2.0d+0*t25*t35*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*t5*t68*wght+1.0d+0*t22*t61
     +        *wght+1.0d+0*t25*t44*wght-3.361232459613792d-3*t3*t35
     +        *t36*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t5*t77*wght
     +        +1.510124173904079d-5*t22*t3*t53*wght
     +        -6.722464919227584d-3*t3*t36*t44*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t86*t9
     +        -3.536776513153214d-2*t11*t13*t16*t85*t9
     +        +8.488263631567714d-1*t11*t17*t60*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t56*t9
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t55*t9
     +        +2.1220659078919285d-1*t11*t16*t30*t55*t9
     +        -6.366197723675786d-1*t11*t16*t32*t53*t9+t16*t86+t16*t85
     +        -3.536776513153214d-2*t11*t13*t17*t81
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t56+3*t16*t32*t55
     +        -1.0610329539459643d-1*t11*t13*t16*t28*t55
     +        -6.366197723675786d-1*t11*t17*t28*t53
     +        -1.0610329539459643d-1*t11*t13*t16*t32*t49
     +        +2.1220659078919285d-1*t11*t17*t30*t49
     +        +4.244131815783857d-1*t11*t16*t28*t30*t32)
     +        +3.5449077018110403d+0*t16*t26*(-7.049892022399704d+0
     +        *rhob*t19*t82*t83-1.8799712059732546d+0*rhob*t19*t27*t81
     +        +8.459870426879645d+0*t19*t50*t51+8.459870426879645d+0
     +        *rhob*t19*t28*t49*t50-5.639913617919763d+0*t19*t27*t49)
     +        +9.973557010035795d-2*t18*t50*t7*t8*t83
     +        +1.329807601338106d-1*t18*t20*t7*t8*t81
     +        +3.989422804014318d-1*t20*t28*t57*t7*t8
     +        -1.994711402007159d-1*t27*t33*t51*t7*t8
     +        +3.989422804014318d-1*t20*t33*t49*t7*t8
     +        -1.994711402007159d-1*t18*t27*t28*t49*t7*t8
     +        -7.978845608028636d-1*t10*t13*t57*t7
     +        -1.5957691216057274d+0*t10*t18*t53*t7
     +        +1.994711402007159d-1*t13*t18*t27*t51*t7
     +        -3.989422804014318d-1*t13*t18*t20*t49*t7
     +        +1.5957691216057274d+0*t10*t30*t33*t7
     +        -7.978845608028636d-1*t13*t20*t28*t33*t7
     +        +7.978845608028636d-1*t18*t20*t28*t30*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t56
     +        +3.5449077018110403d+0*t16*t26*t29*t55
     +        +7.089815403622081d+0*t16*t26*t32*t52)
     +        -1.3298076013381063d-1*t21*t50*t6*t7*t8*t83
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t81
     +        -5.319230405352425d-1*t20*t28*t58*t6*t7*t8
     +        +2.6596152026762127d-1*t27*t34*t51*t6*t7*t8
     +        -5.319230405352425d-1*t20*t34*t49*t6*t7*t8
     +        +2.6596152026762127d-1*t21*t27*t28*t49*t6*t7*t8
     +        +1.063846081070485d+0*t10*t13*t58*t6*t7
     +        +2.12769216214097d+0*t10*t21*t53*t6*t7
     +        -2.6596152026762127d-1*t13*t21*t27*t51*t6*t7
     +        +5.319230405352425d-1*t13*t20*t21*t49*t6*t7
     +        -2.12769216214097d+0*t10*t30*t34*t6*t7
     +        +1.063846081070485d+0*t13*t20*t28*t34*t6*t7
     +        -1.063846081070485d+0*t20*t21*t28*t30*t6*t7)*wght+1.0d+0
     +        *t22*(-7.517062554544749d-4*t3*t46*t80
     +        +1.286568658174524d-6*t3*t78*t79+2.2167190874805223d-1
     +        *t4*t71+2.757113292886222d-1*t71+3.386130477833154d-2
     +        *gammabb*t3*t64)*wght+3.0d+0*t25*t59*wght+3.0d+0*t35*t48
     +        *wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*t5*(-3.5461536035682834d-1
     +        *t10*t6*t7*t8*(2.659615202676212d-1*t10*t7*t8*(
     +        -3.536776513153214d-2*t11*t13*t16*t9*t90+t16*t90
     +        -7.073553026306428d-2*t11*t13*t16*t32*t65*t9
     +        +1.4147106052612857d-1*t11*t16*t30*t65*t9
     +        -7.132756811046163d-4*t11*t17*t3*t89
     +        -3.536776513153214d-2*t11*t13*t17*t88+2*t16*t32*t65
     +        -7.073553026306428d-2*t11*t13*t16*t28*t65
     +        +4.7551712073641084d-4*t11*t32*t39*t64
     +        -7.073553026306428d-2*t11*t13*t16*t32*t62
     +        +1.4147106052612857d-1*t11*t17*t30*t62
     +        +4.751830440958629d-2*t14*t31*t39*t41*t56
     +        -1.1887928018410271d-4*t11*t39*t40*t56
     +        -1.680616229806896d-3*t15*t36*t39*t56
     +        +4.751830440958629d-2*t14*t31*t39*t41*t55
     +        -1.1887928018410271d-4*t11*t39*t40*t55
     +        -1.680616229806896d-3*t15*t36*t39*t55
     +        -1.680616229806896d-3*t31*t36*t39*t49
     +        -1.0083697378841376d-2*t15*t39*t40-3.361232459613792d-3
     +        *t28*t31*t32*t36*t39+6.722464919227584d-3*t15*t24*t32
     +        *t39+6.722464919227584d-3*t24*t28*t31*t39)
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhob*t19*t27*t88-2.369632590246221d-2*t19*t3*t37*t51
     +        *t82+5.639913617919763d+0*rhob*t19*t28*t50*t62
     +        -3.759942411946509d+0*t19*t27*t62+9.478530360984885d-3
     +        *t19*t3*t37*t49*t50+1.895706072196977d-2*t19*t28*t3*t36
     +        *t50)+1.329807601338106d-1*t18*t20*t7*t8*t88
     +        +2.659615202676212d-1*t20*t28*t66*t7*t8
     +        +2.659615202676212d-1*t20*t33*t62*t7*t8
     +        -1.329807601338106d-1*t18*t27*t28*t62*t7*t8
     +        -6.649038006690531d-2*t27*t42*t51*t7*t8
     +        +1.329807601338106d-1*t20*t42*t49*t7*t8
     +        -5.319230405352424d-1*t10*t13*t66*t7
     +        +8.939584949317599d-4*t18*t20*t3*t64*t7
     +        -2.659615202676212d-1*t13*t18*t20*t62*t7
     +        +4.4697924746587997d-4*t20*t24*t3*t57*t7
     +        +3.3523443559941d-4*t18*t24*t3*t50*t51*t7
     +        -2.2348962373293998d-4*t18*t24*t27*t3*t49*t7
     +        +5.319230405352424d-1*t10*t30*t42*t7
     +        -2.659615202676212d-1*t13*t20*t28*t42*t7
     +        -8.939584949317599d-4*t20*t3*t33*t40*t7
     +        +4.4697924746587997d-4*t18*t27*t28*t3*t40*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t33*t7
     +        +3.5449077018110403d+0*t16*t26*t29*t65
     +        +3.5449077018110403d+0*t16*t26*t32*t63
     +        +1.6844800327854395d-1*t14*t26*t31*t39*t41*t52
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t32*t39*t41)
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t88
     +        -3.5461536035682834d-1*t20*t28*t6*t67*t7*t8
     +        -3.5461536035682834d-1*t20*t34*t6*t62*t7*t8
     +        +1.7730768017841417d-1*t21*t27*t28*t6*t62*t7*t8
     +        +8.865384008920708d-2*t27*t43*t51*t6*t7*t8
     +        -1.7730768017841417d-1*t20*t43*t49*t6*t7*t8
     +        +7.092307207136567d-1*t10*t13*t6*t67*t7
     +        -1.1919446599090135d-3*t20*t21*t3*t6*t64*t7
     +        +3.5461536035682834d-1*t13*t20*t21*t6*t62*t7
     +        -5.959723299545067d-4*t20*t24*t3*t58*t6*t7
     +        -4.469792474658801d-4*t21*t24*t3*t50*t51*t6*t7
     +        +2.9798616497725333d-4*t21*t24*t27*t3*t49*t6*t7
     +        -7.092307207136567d-1*t10*t30*t43*t6*t7
     +        +3.5461536035682834d-1*t13*t20*t28*t43*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t34*t40*t6*t7
     +        -5.959723299545067d-4*t21*t27*t28*t3*t40*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t34*t6*t7)*wght
     +        +1.0d+0*t22*(-4.8246324681544644d-7*t3*t46*t87
     +        +2.550431938149111d-4*gammabb*t3*t47
     +        -1.0457167652131799d-2*t3*t40)*wght+2.0d+0*t25*t68*wght
     +        +2.0d+0*t35*t61*wght-3.361232459613792d-3*t3*t36*t59
     +        *wght+1.0d+0*t44*t48*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(-3.536776513153214d-2
     +        *t11*t13*t16*t9*t94+t16*t94-3.536776513153214d-2*t11*t13
     +        *t17*t92-1.0681943420417696d-6*t11*t17*t3*t80
     +        +4.751830440958629d-2*t14*t31*t32*t39*t41*t74
     +        -1.1887928018410271d-4*t11*t32*t39*t40*t74
     +        -1.680616229806896d-3*t15*t32*t36*t39*t74
     +        -1.680616229806896d-3*t28*t31*t36*t39*t74
     +        +3.361232459613792d-3*t15*t24*t39*t74
     +        -3.194401344146213d-4*t14*t32*t45*t54*t72
     +        +1.1297883647561385d-5*t2*t28*t54*t72
     +        +5.648941823780692d-6*t2*t31*t32*t72
     +        +9.503660881917259d-2*t14*t31*t39*t41*t65
     +        -2.3775856036820542d-4*t11*t39*t40*t65
     +        -3.361232459613792d-3*t15*t36*t39*t65
     +        -3.361232459613792d-3*t31*t36*t39*t62)
     +        +7.512005776780308d-7*t18*t27*t7*t70*t93
     +        +3.5449077018110403d+0*t16*t26*(-1.8799712059732546d+0
     +        *rhob*t19*t27*t92-7.964885979694308d-5*t19*t28*t70*t71
     +        *t82+3.1859543918777233d-5*t19*t2*t50*t70
     +        +1.895706072196977d-2*t19*t3*t37*t50*t62
     +        +2.838989964383816d-5*t19*t27*t3*t53
     +        -4.2584849465757246d-5*t19*t28*t3*t30*t50)
     +        +1.329807601338106d-1*t18*t20*t7*t8*t92
     +        +1.329807601338106d-1*t20*t28*t7*t75*t8
     +        +2.659615202676212d-1*t20*t42*t62*t7*t8
     +        -2.659615202676212d-1*t10*t13*t7*t75
     +        +1.6844800327854395d-1*t14*t26*t29*t31*t39*t41*t74
     +        -1.1323857927539448d-3*t14*t26*t29*t45*t54*t72
     +        +1.1268008665170462d-6*t18*t28*t50*t69*t7*t70
     +        -7.512005776780308d-7*t27*t33*t69*t7*t70
     +        +8.939584949317599d-4*t20*t24*t3*t66*t7
     +        -4.4697924746587997d-4*t18*t24*t27*t3*t62*t7
     +        -2.0081746054220725d-6*t20*t3*t33*t60*t7
     +        +1.0040873027110363d-6*t18*t27*t28*t3*t60*t7
     +        +2.0081746054220725d-6*t18*t20*t3*t47*t7
     +        -8.939584949317599d-4*t20*t3*t40*t42*t7
     +        -4.4697924746587997d-4*t24*t27*t28*t3*t42*t7
     +        +3.368960065570879d-1*t14*t26*t31*t39*t41*t63)
     +        -1.0016007702373746d-6*t21*t27*t6*t7*t70*t93
     +        -1.7730768017841417d-1*t20*t21*t6*t7*t8*t92
     +        -1.7730768017841417d-1*t20*t28*t6*t7*t76*t8
     +        -3.5461536035682834d-1*t20*t43*t6*t62*t7*t8
     +        +3.5461536035682834d-1*t10*t13*t6*t7*t76
     +        -1.502401155356062d-6*t21*t28*t50*t6*t69*t7*t70
     +        +1.0016007702373746d-6*t27*t34*t6*t69*t7*t70
     +        -1.1919446599090135d-3*t20*t24*t3*t6*t67*t7
     +        +5.959723299545067d-4*t21*t24*t27*t3*t6*t62*t7
     +        +2.677566140562764d-6*t20*t3*t34*t6*t60*t7
     +        -1.338783070281382d-6*t21*t27*t28*t3*t6*t60*t7
     +        -2.677566140562764d-6*t20*t21*t3*t47*t6*t7
     +        +1.1919446599090135d-3*t20*t3*t40*t43*t6*t7
     +        +5.959723299545067d-4*t24*t27*t28*t3*t43*t6*t7)*wght
     +        +1.0d+0*t22*(1.8092371755579242d-7*gammabb*t3*t91
     +        -6.040496695616315d-5*t3*t60)*wght+1.0d+0*t25*t77*wght
     +        -6.722464919227584d-3*t3*t36*t68*wght+2.0d+0*t44*t61
     +        *wght+1.510124173904079d-5*t3*t35*t53*wght+Cmat3(iq,
     +        D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 1.0d+0*t5*(
     +        -3.5461536035682834d-1*t10*t6*t7*t8*
     +        (2.659615202676212d-1*t10*t7*t8*(4.751830440958629d-2
     +        *t14*t31*t39*t41*t99-1.680616229806896d-3*t15*t36*t39
     +        *t99-3.194401344146213d-4*t14*t45*t54*t72*t98
     +        +3.221137646093454d-6*t13*t14*t84*t97
     +        +2.5379391956869274d-8*t31*t72*t89-2.3995753309638476d-9
     +        *t11*t17*t3*t87-3.194401344146213d-4*t14*t45*t54*t72*t74
     +        +5.648941823780692d-6*t2*t31*t72*t74
     +        +4.751830440958629d-2*t100*t14*t31*t39*t41
     +        -1.680616229806896d-3*t100*t15*t36*t39)
     +        -2.24003435190034d-2*t19*t26*t27*t37*t39*t99
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t98
     +        -9.490367802866481d-7*t19*t26*t30*t82*t97
     +        +1.0124912502475182d-8*t18*t27*t7*t70*t96
     +        +3.7874396480580436d-9*t18*t50*t60*t7*t95
     +        +9.022265057551551d-9*t18*t20*t3*t7*t91
     +        +1.3409377423976399d-3*t20*t24*t3*t7*t75
     +        +1.1293914261386051d-4*t19*t26*t50*t71*t72*t74
     +        -2.2536017330340924d-6*t27*t42*t69*t7*t70
     +        -6.024523816266217d-6*t20*t3*t42*t60*t7
     +        -2.24003435190034d-2*t100*t19*t26*t27*t37*t39)
     +        -1.3499883336633578d-8*t21*t27*t6*t7*t70*t96
     +        -5.049919530744059d-9*t21*t50*t6*t60*t7*t95
     +        -1.202968674340207d-8*t20*t21*t3*t6*t7*t91
     +        -1.7879169898635203d-3*t20*t24*t3*t6*t7*t76
     +        +3.0048023107121236d-6*t27*t43*t6*t69*t7*t70
     +        +8.032698421688291d-6*t20*t3*t43*t6*t60*t7)*wght
     +        -1.0083697378841378d-2*t3*t36*t77*wght
     +        +4.5303725217122376d-5*t3*t44*t53*wght
     +        -6.784639408342217d-8*t22*t3*wght/rhob*
     +        *6.666666666666667d+0+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end subroutine nwxcm_x_camrpbe_d3 
C> @} 
