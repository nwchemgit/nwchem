C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_camb88.F 
C> The nwxcm_x_camb88 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_camb88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; 
C>    \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}
C>    +0.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{
C>    \it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}
C>    \,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{
C>    \it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}+0.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}
C>    \,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{
C>    \it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    +0.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}-0.9305257363491002\,{\it t_{16}}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{
C>    \it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552
C>    \,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}
C>    \,\rho_\beta}\over{{\it t_{22}}}}\right)+0.2659615202676212
C>    \,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({
C>    \it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}
C>    \,\sigma_{\alpha\alpha}-0.9305257363491002\,{\it t_1}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{
C>    \it t_{10}}}}\right)+0.2659615202676212\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left({\it t_{14}}-{{0.03536776513153214\,{
C>    \it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}
C>    \over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    -0.9305257363491002\,{\it t_{25}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{
C>    \it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{
C>    \it t_{31}}}}\right)+0.2659615202676212\,{\it t_7}\,{
C>    \it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({
C>    \it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camb88(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = 1/(gammaaa*t27+1)**5.0d-1
            t29 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t28/rhoa**3.6666666666666664d+0
            t30 = 1/t5**2
            t31 = rhoa**3.333333333333333d-1
            t32 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t31-4.2d-3*gammaaa*t29*t3*t30
            t33 = 1/t12**3
            t34 = 1/t11**2
            t35 = 1.413716694115413d+1*t14*t16*t32*t34-2.827433388230826
     1         d+1*rhoa*t16*t17
            t36 = 1.26d-2*t3*t4/t1+1.26d-2*t27*t28
            t37 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t30*t36
            fnc(iq) = 2*t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = (t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(
     1         2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t11*t
     2         13*t15*t18*t35+t18*t35-3.536776513153214d-2*t13*t15*t19*t
     3         32+7.073553026306428d-2*t11*t13*t19/rhoa**3)*t9+1.3298076
     4         01338106d-1*t10*t20*t22*t32*t9-2.659615202676212d-1*t12*t
     5         15*t20*t9+3.5449077018110403d+0*t18*t25*(3.75994241194650
     6         9d+0*t21*t22-1.8799712059732546d+0*rhoa*t21*t32*t33))-1.7
     7         730768017841417d-1*t10*t22*t23*t32*t8*t9+3.54615360356828
     8         34d-1*t12*t15*t23*t8*t9)+t24*(5.599999999999999d-3*gammaa
     9         a*t26*t6-1.2407009817988002d+0*t31+4.2d-3*gammaaa*t29*t3*
     :         t30))*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t7*(-3.5461536035682834d-1*t10*t12*t8*t9*
     1         (2.659615202676212d-1*t10*t12*(1.413716694115413d+1*t14*t
     2         16*t18*t34*t37-3.536776513153214d-2*t13*t15*t19*t37-5.0d-
     3         1*t17*t18*t37)*t9+1.329807601338106d-1*t10*t20*t22*t37*t9
     4         -6.66432440723758d+0*rhoa*t18*t21*t25*t33*t37)-1.77307680
     5         17841417d-1*t10*t22*t23*t37*t8*t9)+t24*(4.2d-3*gammaaa*t3
     6         *t30*t36-4.2d-3*t3*t6))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = -param(1)
            t9 = param(2)
            t10 = param(3)
            t11 = 1/rhoa
            t12 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t13 = t12**5.0d-1
            t14 = t10**2
            t15 = rhoa**2
            t16 = 1/t15
            t17 = 1/t14
            t18 = 1/t12
            t19 = exp(-1.413716694115413d+1*t15*t17*t18)
            t20 = t19-1.0d+0
            t21 = -3.536776513153214d-2*t12*t14*t16*t20+t19-1.5d+0
            t22 = 1/t10
            t23 = 1/t13
            t24 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         2*t23)+2.659615202676212d-1*t10*t11*t13*t21
            t25 = -3.5461536035682834d-1*t10*t11*t13*t24*t9+t8+1.0d+0
            t26 = gammabb**5.0d-1
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/t27
            t29 = asinh(t26*t28)
            t30 = 2.52d-2*t26*t28*t29+1.0d+0
            t31 = 1/t30
            t32 = -4.2d-3*gammabb*t28*t31-9.305257363491002d-1*t27
            t33 = 1/rhob
            t34 = 4.2d-3*gammabb*t28*t31+9.305257363491002d-1*t27
            t35 = t34**5.0d-1
            t36 = rhob**2
            t37 = 1/t36
            t38 = 1/t34
            t39 = exp(-1.413716694115413d+1*t17*t36*t38)
            t40 = t39-1.0d+0
            t41 = -3.536776513153214d-2*t14*t34*t37*t40+t39-1.5d+0
            t42 = 1/t35
            t43 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         2*t42)+2.659615202676212d-1*t10*t33*t35*t41
            t44 = -3.5461536035682834d-1*t10*t33*t35*t43*t9+t8+1.0d+0
            t45 = 5.641895835477563d-1
            t46 = 1/rhoa**2.3333333333333334d+0
            t47 = 1/rhoa**2.6666666666666666d+0
            t48 = 1/(gammaaa*t47+1)**5.0d-1
            t49 = -3.3599999999999997d-2*gammaaa*t48/rhoa**3.66666666666
     1         66664d+0-3.3599999999999997d-2*t1*t4*t46
            t50 = 1/t5**2
            t51 = rhoa**3.333333333333333d-1
            t52 = -5.599999999999999d-3*gammaaa*t46*t6+1.240700981798800
     1         2d+0*t51-4.2d-3*gammaaa*t3*t49*t50
            t53 = 1/t13**3
            t54 = 1/t12**2
            t55 = 1.413716694115413d+1*t15*t17*t52*t54-2.827433388230826
     1         d+1*rhoa*t17*t18
            t56 = 1/rhob**2.3333333333333334d+0
            t57 = 1/rhob**2.6666666666666666d+0
            t58 = 1/(gammabb*t57+1)**5.0d-1
            t59 = -3.3599999999999997d-2*gammabb*t58/rhob**3.66666666666
     1         66664d+0-3.3599999999999997d-2*t26*t29*t56
            t60 = 1/t30**2
            t61 = rhob**3.333333333333333d-1
            t62 = 1.2407009817988002d+0*t61-4.2d-3*gammabb*t28*t59*t60-5
     1         .599999999999999d-3*gammabb*t31*t56
            t63 = 1/t35**3
            t64 = 1/t34**2
            t65 = 1.413716694115413d+1*t17*t36*t62*t64-2.827433388230826
     1         d+1*rhob*t17*t38
            t66 = 1.26d-2*t47*t48+1.26d-2*t3*t4/t1
            t67 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t50*t66
            t68 = 1.26d-2*t57*t58+1.26d-2*t28*t29/t26
            t69 = 4.2d-3*t28*t31-4.2d-3*gammabb*t28*t60*t68
            fnc(iq) = (t25*t7+t32*t44)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t7*(-3.5461536035682834d-1*t10*t11*t13*(2.
     1         659615202676212d-1*t10*t11*t13*(-3.536776513153214d-2*t12
     2         *t14*t16*t19*t55+t19*t55-3.536776513153214d-2*t14*t16*t20
     3         *t52+7.073553026306428d-2*t12*t14*t20/rhoa**3)+3.54490770
     4         18110403d+0*t19*t45*(3.759942411946509d+0*t22*t23-1.87997
     5         12059732546d+0*rhoa*t22*t52*t53)+1.329807601338106d-1*t10
     6         *t11*t21*t23*t52-2.659615202676212d-1*t10*t13*t16*t21)*t9
     7         -1.7730768017841417d-1*t10*t11*t23*t24*t52*t9+3.546153603
     8         5682834d-1*t10*t13*t16*t24*t9)+t25*(5.599999999999999d-3*
     9         gammaaa*t46*t6-1.2407009817988002d+0*t51+4.2d-3*gammaaa*t
     :         3*t49*t50))*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t32*(-3.5461536035682834d-1*t10*t33*t35*(2
     1         .659615202676212d-1*t10*t33*t35*(-3.536776513153214d-2*t1
     2         4*t34*t37*t39*t65+t39*t65-3.536776513153214d-2*t14*t37*t4
     3         0*t62+7.073553026306428d-2*t14*t34*t40/rhob**3)+3.5449077
     4         018110403d+0*t39*t45*(3.759942411946509d+0*t22*t42-1.8799
     5         712059732546d+0*rhob*t22*t62*t63)+1.329807601338106d-1*t1
     6         0*t33*t41*t42*t62-2.659615202676212d-1*t10*t35*t37*t41)*t
     7         9-1.7730768017841417d-1*t10*t33*t42*t43*t62*t9+3.54615360
     8         35682834d-1*t10*t35*t37*t43*t9)+t44*(-1.2407009817988002d
     9         +0*t61+4.2d-3*gammabb*t28*t59*t60+5.599999999999999d-3*ga
     :         mmabb*t31*t56))*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t7*(-3.5461536035682834d-1*t10*t11*t13*(2
     1         .659615202676212d-1*t10*t11*t13*(1.413716694115413d+1*t15
     2         *t17*t19*t54*t67-3.536776513153214d-2*t14*t16*t20*t67-5.0
     3         d-1*t18*t19*t67)-6.66432440723758d+0*rhoa*t19*t22*t45*t53
     4         *t67+1.329807601338106d-1*t10*t11*t21*t23*t67)*t9-1.77307
     5         68017841417d-1*t10*t11*t23*t24*t67*t9)+t25*(4.2d-3*gammaa
     6         a*t3*t50*t66-4.2d-3*t3*t6))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t32*(-3.5461536035682834d-1*t10*t33*t35*(
     1         2.659615202676212d-1*t10*t33*t35*(1.413716694115413d+1*t1
     2         7*t36*t39*t64*t69-3.536776513153214d-2*t14*t37*t40*t69-5.
     3         0d-1*t38*t39*t69)-6.66432440723758d+0*rhob*t22*t39*t45*t6
     4         3*t69+1.329807601338106d-1*t10*t33*t41*t42*t69)*t9-1.7730
     5         768017841417d-1*t10*t33*t42*t43*t69*t9)+t44*(4.2d-3*gamma
     6         bb*t28*t60*t68-4.2d-3*t28*t31))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = 1/(gammaaa*t27+1)**5.0d-1
            t29 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t28/rhoa**3.6666666666666664d+0
            t30 = 1/t5**2
            t31 = rhoa**3.333333333333333d-1
            t32 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t31-4.2d-3*gammaaa*t29*t3*t30
            t33 = 1/t12**3
            t34 = 1/t11**2
            t35 = 1.413716694115413d+1*t14*t16*t32*t34-2.827433388230826
     1         d+1*rhoa*t16*t17
            t36 = 1.26d-2*t3*t4/t1+1.26d-2*t27*t28
            t37 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t30*t36
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(2
     1         .659615202676212d-1*t10*t12*(-3.536776513153214d-2*t11*t1
     2         3*t15*t18*t35+t18*t35-3.536776513153214d-2*t13*t15*t19*t3
     3         2+7.073553026306428d-2*t11*t13*t19/rhoa**3)*t9+1.32980760
     4         1338106d-1*t10*t20*t22*t32*t9-2.659615202676212d-1*t12*t1
     5         5*t20*t9+3.5449077018110403d+0*t18*t25*(3.759942411946509
     6         d+0*t21*t22-1.8799712059732546d+0*rhoa*t21*t32*t33))-1.77
     7         30768017841417d-1*t10*t22*t23*t32*t8*t9+3.546153603568283
     8         4d-1*t12*t15*t23*t8*t9)*wght+t24*(5.599999999999999d-3*ga
     9         mmaaa*t26*t6-1.2407009817988002d+0*t31+4.2d-3*gammaaa*t29
     :         *t3*t30)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(
     1         2.659615202676212d-1*t10*t12*(1.413716694115413d+1*t14*t1
     2         6*t18*t34*t37-3.536776513153214d-2*t13*t15*t19*t37-5.0d-1
     3         *t17*t18*t37)*t9+1.329807601338106d-1*t10*t20*t22*t37*t9-
     4         6.66432440723758d+0*rhoa*t18*t21*t25*t33*t37)-1.773076801
     5         7841417d-1*t10*t22*t23*t37*t8*t9)*wght+t24*(4.2d-3*gammaa
     6         a*t3*t30*t36-4.2d-3*t3*t6)*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = rhob**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammabb*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhob
            t11 = 4.2d-3*gammabb*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhob**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhob*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhob**2.3333333333333334d+0
            t27 = 1/rhob**2.6666666666666666d+0
            t28 = 1/(gammabb*t27+1)**5.0d-1
            t29 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammabb*t28/rhob**3.6666666666666664d+0
            t30 = 1/t5**2
            t31 = rhob**3.333333333333333d-1
            t32 = -5.599999999999999d-3*gammabb*t26*t6+1.240700981798800
     1         2d+0*t31-4.2d-3*gammabb*t29*t3*t30
            t33 = 1/t12**3
            t34 = 1/t11**2
            t35 = 1.413716694115413d+1*t14*t16*t32*t34-2.827433388230826
     1         d+1*rhob*t16*t17
            t36 = 1.26d-2*t3*t4/t1+1.26d-2*t27*t28
            t37 = 4.2d-3*t3*t6-4.2d-3*gammabb*t3*t30*t36
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RB) = t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(2
     1         .659615202676212d-1*t10*t12*(-3.536776513153214d-2*t11*t1
     2         3*t15*t18*t35+t18*t35-3.536776513153214d-2*t13*t15*t19*t3
     3         2+7.073553026306428d-2*t11*t13*t19/rhob**3)*t9+1.32980760
     4         1338106d-1*t10*t20*t22*t32*t9-2.659615202676212d-1*t12*t1
     5         5*t20*t9+3.5449077018110403d+0*t18*t25*(3.759942411946509
     6         d+0*t21*t22-1.8799712059732546d+0*rhob*t21*t32*t33))-1.77
     7         30768017841417d-1*t10*t22*t23*t32*t8*t9+3.546153603568283
     8         4d-1*t12*t15*t23*t8*t9)*wght+t24*(5.599999999999999d-3*ga
     9         mmabb*t26*t6-1.2407009817988002d+0*t31+4.2d-3*gammabb*t29
     :         *t3*t30)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(
     1         2.659615202676212d-1*t10*t12*(1.413716694115413d+1*t14*t1
     2         6*t18*t34*t37-3.536776513153214d-2*t13*t15*t19*t37-5.0d-1
     3         *t17*t18*t37)*t9+1.329807601338106d-1*t10*t20*t22*t37*t9-
     4         6.66432440723758d+0*rhob*t18*t21*t25*t33*t37)-1.773076801
     5         7841417d-1*t10*t22*t23*t37*t8*t9)*wght+t24*(4.2d-3*gammab
     6         b*t3*t30*t36-4.2d-3*t3*t6)*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_camb88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; 
C>    \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}
C>    +0.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{
C>    \it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}
C>    \,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{
C>    \it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}+0.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}
C>    \,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{
C>    \it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    +0.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}-0.9305257363491002\,{\it t_{16}}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{
C>    \it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552
C>    \,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}
C>    \,\rho_\beta}\over{{\it t_{22}}}}\right)+0.2659615202676212
C>    \,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({
C>    \it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}
C>    \,\sigma_{\alpha\alpha}-0.9305257363491002\,{\it t_1}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{
C>    \it t_{10}}}}\right)+0.2659615202676212\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left({\it t_{14}}-{{0.03536776513153214\,{
C>    \it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}
C>    \over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    -0.9305257363491002\,{\it t_{25}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{
C>    \it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{
C>    \it t_{31}}}}\right)+0.2659615202676212\,{\it t_7}\,{
C>    \it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({
C>    \it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camb88_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = (gammaaa*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhoa**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t29*t30
            t32 = 1/t5**2
            t33 = rhoa**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammaaa*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         a*t21*t34*t35
            t37 = rhoa**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhoa*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammaaa*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammaaa*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammaaa*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = 1/t28**3
            t58 = -4.4799999999999995d-2*gammaaa**2*t57/rhoa**7.33333333
     1         3333333d+0+7.84d-2*t1*t4*t54+1.6799999999999998d-1*gammaa
     2         a*t29/rhoa**4.666666666666667d+0
            t59 = 1/rhoa**6.666666666666666d-1
            t60 = 1.3066666666666665d-2*gammaaa*t54*t6+4.135669939329334
     1         d-1*t59-4.2d-3*gammaaa*t3*t32*t58+8.4d-3*gammaaa*t3*t55*t
     2         56+1.1199999999999999d-2*gammaaa*t26*t31*t32
            t61 = rhoa**4
            t62 = 1/t11**3
            t63 = -2.827433388230826d+1*t14*t16*t52*t62+1.41371669411541
     1         3d+1*t14*t16*t39*t60+5.654866776461652d+1*rhoa*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t64 = t40**2
            t65 = 1.6799999999999998d-2*gammaaa*t57/rhoa**6.333333333333
     1         333d+0-1.6799999999999998d-2*t26*t4*t45-5.04d-2*t29*t30
            t66 = -4.2d-3*gammaaa*t3*t32*t65-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammaaa*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         aa*t26*t32*t46-4.2d-3*t3*t31*t32
            t67 = 1.413716694115413d+1*t14*t16*t39*t66-2.827433388230826
     1         d+1*t14*t16*t34*t47*t62+2.827433388230826d+1*rhoa*t16*t39
     2         *t47
            t68 = t47**2
            t69 = t46**2
            t70 = -6.3d-3*t57/rhoa**5.333333333333333d+0-6.3d-3*t3*t4/t1
     1         **3+6.3d-3*t27*t29/gammaaa
            t71 = -4.2d-3*gammaaa*t3*t32*t70+8.4d-3*gammaaa*t3*t56*t69-8
     1         .4d-3*t3*t32*t46
            fnc(iq) = 2*t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = (t43*t7+t24*t44)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t50*t7+t24*t51)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t7*(-3.5461536035682834d-1*t10*t12*t8*
     1         t9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t
     2         11*t13*t15*t18*t64+t18*t64-3.536776513153214d-2*t11*t13*t
     3         15*t18*t63+t18*t63-2.1220659078919285d-1*t11*t13*t19/t61-
     4         3.536776513153214d-2*t13*t15*t19*t60+1.4147106052612857d-
     5         1*t11*t13*t18*t38*t40-7.073553026306428d-2*t13*t15*t18*t3
     6         4*t40+1.4147106052612857d-1*t13*t19*t34*t38)*t9+1.3298076
     7         01338106d-1*t10*t20*t22*t60*t9-6.649038006690531d-2*t10*t
     8         20*t35*t52*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.3
     9         19230405352424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12
     :         *t20*t38*t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.544
     ;         9077018110403d+0*t18*t25*(-1.8799712059732546d+0*rhoa*t21
     <         *t35*t60+2.8199568089598814d+0*rhoa*t21*t52*t53-3.7599424
     =         11946509d+0*t21*t34*t35)+3.5449077018110403d+0*t18*t25*t3
     >         6*t40)-1.7730768017841417d-1*t10*t22*t23*t60*t8*t9+8.8653
     ?         84008920708d-2*t10*t23*t35*t52*t8*t9-3.5461536035682834d-
     @         1*t10*t22*t34*t42*t8*t9+7.092307207136567d-1*t12*t15*t42*
     1         t8*t9-7.092307207136567d-1*t12*t23*t38*t8*t9+3.5461536035
     2         682834d-1*t15*t22*t23*t34*t8*t9)+t24*(-1.3066666666666665
     3         d-2*gammaaa*t54*t6-4.135669939329334d-1*t59+4.2d-3*gammaa
     4         a*t3*t32*t58-8.4d-3*gammaaa*t3*t55*t56-1.1199999999999999
     5         d-2*gammaaa*t26*t31*t32)+2*t43*t44)*wght+Amat2(iq,D2_RA_R
     6         A)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t7*(-3.5461536035682834d-1*t10*t12*t8
     1         *t9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*
     2         t11*t13*t15*t18*t67+t18*t67-3.536776513153214d-2*t13*t15*
     3         t19*t66+1.413716694115413d+1*t14*t16*t18*t39*t40*t47-5.0d
     4         -1*t17*t18*t40*t47-3.536776513153214d-2*t13*t15*t18*t40*t
     5         47-5.0d-1*t18*t34*t39*t47+7.073553026306428d-2*t13*t19*t3
     6         8*t47+1.0d+0*t10*t17*t18*t47)*t9+1.329807601338106d-1*t10
     7         *t20*t22*t66*t9+1.329807601338106d-1*t10*t22*t34*t48*t9-2
     8         .659615202676212d-1*t12*t15*t48*t9+1.329807601338106d-1*t
     9         10*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*t35*t4
     :         7*t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.5449077018
     ;         110403d+0*t18*t25*(-1.8799712059732546d+0*rhoa*t21*t35*t6
     <         6+2.8199568089598814d+0*rhoa*t21*t34*t47*t53-1.8799712059
     =         732546d+0*t21*t35*t47)+5.01149519714857d+1*t14*t16*t18*t2
     >         5*t36*t39*t47)-1.7730768017841417d-1*t10*t22*t23*t66*t8*t
     ?         9-1.7730768017841417d-1*t10*t22*t34*t49*t8*t9+3.546153603
     @         5682834d-1*t12*t15*t49*t8*t9-1.7730768017841417d-1*t10*t2
     1         2*t42*t47*t8*t9+8.865384008920708d-2*t10*t23*t34*t35*t47*
     2         t8*t9+1.7730768017841417d-1*t15*t22*t23*t47*t8*t9)+t24*(4
     3         .2d-3*gammaaa*t3*t32*t65+5.599999999999999d-3*t26*t6-8.4d
     4         -3*gammaaa*t3*t31*t46*t56-5.599999999999999d-3*gammaaa*t2
     5         6*t32*t46+4.2d-3*t3*t31*t32)+t43*t51+t44*t50)*wght+Cmat2(
     6         iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t7*(-3.5461536035682834d-1*t10*t12*t
     1         8*t9*(2.659615202676212d-1*t10*t12*(1.9985948912206125d+2
     2         *t18*t61*t68/(t11**4*t9**4)+1.413716694115413d+1*t14*t16*
     3         t18*t39*t71-3.536776513153214d-2*t13*t15*t19*t71-5.0d-1*t
     4         17*t18*t71-3.534291735288533d+1*t14*t16*t18*t62*t68)*t9+1
     5         .329807601338106d-1*t10*t20*t22*t71*t9-6.649038006690531d
     6         -2*t10*t20*t35*t68*t9+2.659615202676212d-1*t10*t22*t47*t4
     7         8*t9-9.421466669512572d+1*t18*t25*t37*t68/(t12**7*t9**3)-
     8         6.66432440723758d+0*rhoa*t18*t21*t25*t35*t71+9.9964866108
     9         56369d+0*rhoa*t18*t21*t25*t53*t68)-1.7730768017841417d-1*
     :         t10*t22*t23*t71*t8*t9+8.865384008920708d-2*t10*t23*t35*t6
     ;         8*t8*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8*t9)+t24*
     <         (4.2d-3*gammaaa*t3*t32*t70-8.4d-3*gammaaa*t3*t56*t69+8.4d
     =         -3*t3*t32*t46)+2*t50*t51)*wght+Cmat2(iq,D2_GAA_GAA)
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
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = -param(1)
            t9 = param(2)
            t10 = param(3)
            t11 = 1/rhoa
            t12 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t13 = t12**5.0d-1
            t14 = t10**2
            t15 = rhoa**2
            t16 = 1/t15
            t17 = 1/t14
            t18 = 1/t12
            t19 = exp(-1.413716694115413d+1*t15*t17*t18)
            t20 = t19-1.0d+0
            t21 = -3.536776513153214d-2*t12*t14*t16*t20+t19-1.5d+0
            t22 = 1/t10
            t23 = 1/t13
            t24 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         2*t23)+2.659615202676212d-1*t10*t11*t13*t21
            t25 = -3.5461536035682834d-1*t10*t11*t13*t24*t9+t8+1.0d+0
            t26 = gammabb**5.0d-1
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/t27
            t29 = asinh(t26*t28)
            t30 = 2.52d-2*t26*t28*t29+1.0d+0
            t31 = 1/t30
            t32 = -4.2d-3*gammabb*t28*t31-9.305257363491002d-1*t27
            t33 = 1/rhob
            t34 = 4.2d-3*gammabb*t28*t31+9.305257363491002d-1*t27
            t35 = t34**5.0d-1
            t36 = rhob**2
            t37 = 1/t36
            t38 = 1/t34
            t39 = exp(-1.413716694115413d+1*t17*t36*t38)
            t40 = t39-1.0d+0
            t41 = -3.536776513153214d-2*t14*t34*t37*t40+t39-1.5d+0
            t42 = 1/t35
            t43 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         2*t42)+2.659615202676212d-1*t10*t33*t35*t41
            t44 = -3.5461536035682834d-1*t10*t33*t35*t43*t9+t8+1.0d+0
            t45 = 5.641895835477563d-1
            t46 = 1/rhoa**2.3333333333333334d+0
            t47 = 1/rhoa**2.6666666666666666d+0
            t48 = (gammaaa*t47+1)**5.0d-1
            t49 = 1/t48
            t50 = 1/rhoa**3.6666666666666664d+0
            t51 = -3.3599999999999997d-2*gammaaa*t49*t50-3.3599999999999
     1         997d-2*t1*t4*t46
            t52 = 1/t5**2
            t53 = rhoa**3.333333333333333d-1
            t54 = -5.599999999999999d-3*gammaaa*t46*t6+1.240700981798800
     1         2d+0*t53-4.2d-3*gammaaa*t3*t51*t52
            t55 = 1/t13**3
            t56 = 3.759942411946509d+0*t22*t23-1.8799712059732546d+0*rho
     1         a*t22*t54*t55
            t57 = rhoa**3
            t58 = 1/t57
            t59 = 1/t12**2
            t60 = 1.413716694115413d+1*t15*t17*t54*t59-2.827433388230826
     1         d+1*rhoa*t17*t18
            t61 = -3.536776513153214d-2*t12*t14*t16*t19*t60+t19*t60+7.07
     1         3553026306428d-2*t12*t14*t20*t58-3.536776513153214d-2*t14
     2         *t16*t20*t54
            t62 = 2.659615202676212d-1*t10*t11*t13*t61+3.544907701811040
     1         3d+0*t19*t45*t56+1.329807601338106d-1*t10*t11*t21*t23*t54
     2         -2.659615202676212d-1*t10*t13*t16*t21
            t63 = -3.5461536035682834d-1*t10*t11*t13*t62*t9-1.7730768017
     1         841417d-1*t10*t11*t23*t24*t54*t9+3.5461536035682834d-1*t1
     2         0*t13*t16*t24*t9
            t64 = 5.599999999999999d-3*gammaaa*t46*t6-1.2407009817988002
     1         d+0*t53+4.2d-3*gammaaa*t3*t51*t52
            t65 = 1/rhob**2.3333333333333334d+0
            t66 = 1/rhob**2.6666666666666666d+0
            t67 = (gammabb*t66+1)**5.0d-1
            t68 = 1/t67
            t69 = 1/rhob**3.6666666666666664d+0
            t70 = -3.3599999999999997d-2*gammabb*t68*t69-3.3599999999999
     1         997d-2*t26*t29*t65
            t71 = 1/t30**2
            t72 = rhob**3.333333333333333d-1
            t73 = 1.2407009817988002d+0*t72-4.2d-3*gammabb*t28*t70*t71-5
     1         .599999999999999d-3*gammabb*t31*t65
            t74 = 1/t35**3
            t75 = 3.759942411946509d+0*t22*t42-1.8799712059732546d+0*rho
     1         b*t22*t73*t74
            t76 = rhob**3
            t77 = 1/t76
            t78 = 1/t34**2
            t79 = 1.413716694115413d+1*t17*t36*t73*t78-2.827433388230826
     1         d+1*rhob*t17*t38
            t80 = -3.536776513153214d-2*t14*t34*t37*t39*t79+t39*t79+7.07
     1         3553026306428d-2*t14*t34*t40*t77-3.536776513153214d-2*t14
     2         *t37*t40*t73
            t81 = 2.659615202676212d-1*t10*t33*t35*t80+3.544907701811040
     1         3d+0*t39*t45*t75+1.329807601338106d-1*t10*t33*t41*t42*t73
     2         -2.659615202676212d-1*t10*t35*t37*t41
            t82 = -3.5461536035682834d-1*t10*t33*t35*t81*t9-1.7730768017
     1         841417d-1*t10*t33*t42*t43*t73*t9+3.5461536035682834d-1*t1
     2         0*t35*t37*t43*t9
            t83 = -1.2407009817988002d+0*t72+4.2d-3*gammabb*t28*t70*t71+
     1         5.599999999999999d-3*gammabb*t31*t65
            t84 = 1/t1
            t85 = 1.26d-2*t3*t4*t84+1.26d-2*t47*t49
            t86 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t52*t85
            t87 = 1.413716694115413d+1*t15*t17*t19*t59*t86-3.53677651315
     1         3214d-2*t14*t16*t20*t86-5.0d-1*t18*t19*t86
            t88 = 2.659615202676212d-1*t10*t11*t13*t87-6.66432440723758d
     1         +0*rhoa*t19*t22*t45*t55*t86+1.329807601338106d-1*t10*t11*
     2         t21*t23*t86
            t89 = -3.5461536035682834d-1*t10*t11*t13*t88*t9-1.7730768017
     1         841417d-1*t10*t11*t23*t24*t86*t9
            t90 = 4.2d-3*gammaaa*t3*t52*t85-4.2d-3*t3*t6
            t91 = 1/t26
            t92 = 1.26d-2*t28*t29*t91+1.26d-2*t66*t68
            t93 = 4.2d-3*t28*t31-4.2d-3*gammabb*t28*t71*t92
            t94 = 1.413716694115413d+1*t17*t36*t39*t78*t93-3.53677651315
     1         3214d-2*t14*t37*t40*t93-5.0d-1*t38*t39*t93
            t95 = 2.659615202676212d-1*t10*t33*t35*t94-6.66432440723758d
     1         +0*rhob*t22*t39*t45*t74*t93+1.329807601338106d-1*t10*t33*
     2         t41*t42*t93
            t96 = -3.5461536035682834d-1*t10*t33*t35*t9*t95-1.7730768017
     1         841417d-1*t10*t33*t42*t43*t9*t93
            t97 = 4.2d-3*gammabb*t28*t71*t92-4.2d-3*t28*t31
            t98 = t54**2
            t99 = 1/t13**5
            t100 = 1/rhoa**3.3333333333333337d+0
            t101 = t51**2
            t102 = 1/t5**3
            t103 = 1/t48**3
            t104 = 1.6799999999999998d-1*gammaaa*t49/rhoa**4.66666666666
     1         6667d+0+7.84d-2*t1*t100*t4-4.4799999999999995d-2*gammaaa*
     2         *2*t103/rhoa**7.333333333333333d+0
            t105 = 1/rhoa**6.666666666666666d-1
            t106 = 1.3066666666666665d-2*gammaaa*t100*t6+1.1199999999999
     1         999d-2*gammaaa*t46*t51*t52-4.2d-3*gammaaa*t104*t3*t52+8.4
     2         d-3*gammaaa*t101*t102*t3+4.135669939329334d-1*t105
            t107 = rhoa**4
            t108 = 1/t12**3
            t109 = -2.827433388230826d+1*t108*t15*t17*t98+5.654866776461
     1         652d+1*rhoa*t17*t54*t59+1.413716694115413d+1*t106*t15*t17
     2         *t59-2.827433388230826d+1*t17*t18
            t110 = t60**2
            t111 = t73**2
            t112 = 1/t35**5
            t113 = 1/rhob**3.3333333333333337d+0
            t114 = t70**2
            t115 = 1/t30**3
            t116 = 1/t67**3
            t117 = 1.6799999999999998d-1*gammabb*t68/rhob**4.66666666666
     1         6667d+0+7.84d-2*t113*t26*t29-4.4799999999999995d-2*gammab
     2         b**2*t116/rhob**7.333333333333333d+0
            t118 = 1/rhob**6.666666666666666d-1
            t119 = 1.1199999999999999d-2*gammabb*t65*t70*t71-4.2d-3*gamm
     1         abb*t117*t28*t71+1.3066666666666665d-2*gammabb*t113*t31+8
     2         .4d-3*gammabb*t114*t115*t28+4.135669939329334d-1*t118
            t120 = rhob**4
            t121 = 1/t34**3
            t122 = 5.654866776461652d+1*rhob*t17*t73*t78+1.4137166941154
     1         13d+1*t119*t17*t36*t78-2.827433388230826d+1*t17*t38-2.827
     2         433388230826d+1*t111*t121*t17*t36
            t123 = t79**2
            t124 = -1.6799999999999998d-2*t4*t46*t84-5.04d-2*t49*t50+1.6
     1         799999999999998d-2*gammaaa*t103/rhoa**6.333333333333333d+
     2         0
            t125 = 5.599999999999999d-3*gammaaa*t46*t52*t85+8.4d-3*gamma
     1         aa*t102*t3*t51*t85-5.599999999999999d-3*t46*t6-4.2d-3*t3*
     2         t51*t52-4.2d-3*gammaaa*t124*t3*t52
            t126 = 2.827433388230826d+1*rhoa*t17*t59*t86-2.8274333882308
     1         26d+1*t108*t15*t17*t54*t86+1.413716694115413d+1*t125*t15*
     2         t17*t59
            t127 = -1.6799999999999998d-2*t29*t65*t91-5.04d-2*t68*t69+1.
     1         6799999999999998d-2*gammabb*t116/rhob**6.333333333333333d
     2         +0
            t128 = 5.599999999999999d-3*gammabb*t65*t71*t92+8.4d-3*gamma
     1         bb*t115*t28*t70*t92-4.2d-3*t28*t70*t71-4.2d-3*gammabb*t12
     2         7*t28*t71-5.599999999999999d-3*t31*t65
            t129 = 2.827433388230826d+1*rhob*t17*t78*t93-2.8274333882308
     1         26d+1*t121*t17*t36*t73*t93+1.413716694115413d+1*t128*t17*
     2         t36*t78
            t130 = 1/t10**3
            t131 = t86**2
            t132 = t85**2
            t133 = 6.3d-3*t47*t49/gammaaa-6.3d-3*t3*t4/t1**3-6.3d-3*t103
     1         /rhoa**5.333333333333333d+0
            t134 = -8.4d-3*t3*t52*t85-4.2d-3*gammaaa*t133*t3*t52+8.4d-3*
     1         gammaaa*t102*t132*t3
            t135 = 1/t10**4
            t136 = t93**2
            t137 = t92**2
            t138 = 6.3d-3*t66*t68/gammabb-6.3d-3*t28*t29/t26**3-6.3d-3*t
     1         116/rhob**5.333333333333333d+0
            t139 = -8.4d-3*t28*t71*t92-4.2d-3*gammabb*t138*t28*t71+8.4d-
     1         3*gammabb*t115*t137*t28
            fnc(iq) = (t25*t7+t32*t44)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t63*t7+t25*t64)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t44*t83+t32*t82)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t25*t90+t7*t89)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t44*t97+t32*t96)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t7*(-3.5461536035682834d-1*t10*t11*t13
     1         *t9*(3.5449077018110403d+0*t19*t45*(2.8199568089598814d+0
     2         *rhoa*t22*t98*t99-3.759942411946509d+0*t22*t54*t55-1.8799
     3         712059732546d+0*rhoa*t106*t22*t55)-6.649038006690531d-2*t
     4         10*t11*t21*t55*t98+2.659615202676212d-1*t10*t11*t23*t54*t
     5         61-5.319230405352424d-1*t10*t13*t16*t61+2.659615202676212
     6         d-1*t10*t11*t13*(1.4147106052612857d-1*t12*t14*t19*t58*t6
     7         0-7.073553026306428d-2*t14*t16*t19*t54*t60+1.414710605261
     8         2857d-1*t14*t20*t54*t58-3.536776513153214d-2*t106*t14*t16
     9         *t20-2.1220659078919285d-1*t12*t14*t20/t107-3.53677651315
     :         3214d-2*t110*t12*t14*t16*t19-3.536776513153214d-2*t109*t1
     ;         2*t14*t16*t19+t110*t19+t109*t19)+3.5449077018110403d+0*t1
     <         9*t45*t56*t60+5.319230405352424d-1*t10*t13*t21*t58-2.6596
     =         15202676212d-1*t10*t16*t21*t23*t54+1.329807601338106d-1*t
     >         10*t106*t11*t21*t23)+8.865384008920708d-2*t10*t11*t24*t55
     ?         *t9*t98-3.5461536035682834d-1*t10*t11*t23*t54*t62*t9+7.09
     @         2307207136567d-1*t10*t13*t16*t62*t9-7.092307207136567d-1*
     1         t10*t13*t24*t58*t9+3.5461536035682834d-1*t10*t16*t23*t24*
     2         t54*t9-1.7730768017841417d-1*t10*t106*t11*t23*t24*t9)+2*t
     3         63*t64+t25*(-1.3066666666666665d-2*gammaaa*t100*t6-1.1199
     4         999999999999d-2*gammaaa*t46*t51*t52+4.2d-3*gammaaa*t104*t
     5         3*t52-8.4d-3*gammaaa*t101*t102*t3-4.135669939329334d-1*t1
     6         05))*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t32*(-3.5461536035682834d-1*t10*t33*t4
     1         2*t73*t81*t9+7.092307207136567d-1*t10*t35*t37*t81*t9-3.54
     2         61536035682834d-1*t10*t33*t35*(2.659615202676212d-1*t10*t
     3         33*t42*t73*t80-5.319230405352424d-1*t10*t35*t37*t80+2.659
     4         615202676212d-1*t10*t33*t35*(1.4147106052612857d-1*t14*t3
     5         4*t39*t77*t79-7.073553026306428d-2*t14*t37*t39*t73*t79+1.
     6         4147106052612857d-1*t14*t40*t73*t77-3.536776513153214d-2*
     7         t119*t14*t37*t40-2.1220659078919285d-1*t14*t34*t40/t120-3
     8         .536776513153214d-2*t123*t14*t34*t37*t39-3.53677651315321
     9         4d-2*t122*t14*t34*t37*t39+t123*t39+t122*t39)+3.5449077018
     :         110403d+0*t39*t45*t75*t79+5.319230405352424d-1*t10*t35*t4
     ;         1*t77+3.5449077018110403d+0*t39*t45*(-3.759942411946509d+
     <         0*t22*t73*t74-1.8799712059732546d+0*rhob*t119*t22*t74+2.8
     =         199568089598814d+0*rhob*t111*t112*t22)-6.649038006690531d
     >         -2*t10*t111*t33*t41*t74-2.659615202676212d-1*t10*t37*t41*
     ?         t42*t73+1.329807601338106d-1*t10*t119*t33*t41*t42)*t9-7.0
     @         92307207136567d-1*t10*t35*t43*t77*t9+8.865384008920708d-2
     1         *t10*t111*t33*t43*t74*t9+3.5461536035682834d-1*t10*t37*t4
     2         2*t43*t73*t9-1.7730768017841417d-1*t10*t119*t33*t42*t43*t
     3         9)+2*t82*t83+t44*(-1.1199999999999999d-2*gammabb*t65*t70*
     4         t71+4.2d-3*gammabb*t117*t28*t71-1.3066666666666665d-2*gam
     5         mabb*t113*t31-8.4d-3*gammabb*t114*t115*t28-4.135669939329
     6         334d-1*t118))*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t7*(-3.5461536035682834d-1*t10*t11*t1
     1         3*t9*(3.5449077018110403d+0*t19*t45*(2.8199568089598814d+
     2         0*rhoa*t22*t54*t86*t99-1.8799712059732546d+0*t22*t55*t86-
     3         1.8799712059732546d+0*rhoa*t125*t22*t55)+1.32980760133810
     4         6d-1*t10*t11*t23*t54*t87-2.659615202676212d-1*t10*t13*t16
     5         *t87+2.659615202676212d-1*t10*t11*t13*(1.413716694115413d
     6         +1*t15*t17*t19*t59*t60*t86-5.0d-1*t18*t19*t60*t86-3.53677
     7         6513153214d-2*t14*t16*t19*t60*t86-5.0d-1*t19*t54*t59*t86+
     8         7.073553026306428d-2*t14*t20*t58*t86+1.0d+0*t11*t18*t19*t
     9         86-3.536776513153214d-2*t125*t14*t16*t20-3.53677651315321
     :         4d-2*t12*t126*t14*t16*t19+t126*t19)+1.329807601338106d-1*
     ;         t10*t11*t23*t61*t86+5.01149519714857d+1*t15*t17*t19*t45*t
     <         56*t59*t86-6.649038006690531d-2*t10*t11*t21*t54*t55*t86-1
     =         .329807601338106d-1*t10*t16*t21*t23*t86+1.329807601338106
     >         d-1*t10*t11*t125*t21*t23)-1.7730768017841417d-1*t10*t11*t
     ?         23*t54*t88*t9+3.5461536035682834d-1*t10*t13*t16*t88*t9-1.
     @         7730768017841417d-1*t10*t11*t23*t62*t86*t9+8.865384008920
     1         708d-2*t10*t11*t24*t54*t55*t86*t9+1.7730768017841417d-1*t
     2         10*t16*t23*t24*t86*t9-1.7730768017841417d-1*t10*t11*t125*
     3         t23*t24*t9)+t63*t90+t64*t89+t25*(-5.599999999999999d-3*ga
     4         mmaaa*t46*t52*t85-8.4d-3*gammaaa*t102*t3*t51*t85+5.599999
     5         999999999d-3*t46*t6+4.2d-3*t3*t51*t52+4.2d-3*gammaaa*t124
     6         *t3*t52))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t82*t97+t83*t96+t32*(-1.7730768017841
     1         417d-1*t10*t33*t42*t73*t9*t95+3.5461536035682834d-1*t10*t
     2         35*t37*t9*t95-3.5461536035682834d-1*t10*t33*t35*t9*(1.329
     3         807601338106d-1*t10*t33*t42*t73*t94-2.659615202676212d-1*
     4         t10*t35*t37*t94+2.659615202676212d-1*t10*t33*t35*(1.41371
     5         6694115413d+1*t17*t36*t39*t78*t79*t93-5.0d-1*t38*t39*t79*
     6         t93-3.536776513153214d-2*t14*t37*t39*t79*t93-5.0d-1*t39*t
     7         73*t78*t93+7.073553026306428d-2*t14*t40*t77*t93+1.0d+0*t3
     8         3*t38*t39*t93-3.536776513153214d-2*t128*t14*t37*t40-3.536
     9         776513153214d-2*t129*t14*t34*t37*t39+t129*t39)+3.54490770
     :         18110403d+0*t39*t45*(-1.8799712059732546d+0*t22*t74*t93+2
     ;         .8199568089598814d+0*rhob*t112*t22*t73*t93-1.879971205973
     <         2546d+0*rhob*t128*t22*t74)+1.329807601338106d-1*t10*t33*t
     =         42*t80*t93+5.01149519714857d+1*t17*t36*t39*t45*t75*t78*t9
     >         3-6.649038006690531d-2*t10*t33*t41*t73*t74*t93-1.32980760
     ?         1338106d-1*t10*t37*t41*t42*t93+1.329807601338106d-1*t10*t
     @         128*t33*t41*t42)-1.7730768017841417d-1*t10*t33*t42*t81*t9
     1         *t93+8.865384008920708d-2*t10*t33*t43*t73*t74*t9*t93+1.77
     2         30768017841417d-1*t10*t37*t42*t43*t9*t93-1.77307680178414
     3         17d-1*t10*t128*t33*t42*t43*t9)+t44*(-5.599999999999999d-3
     4         *gammabb*t65*t71*t92-8.4d-3*gammabb*t115*t28*t70*t92+4.2d
     5         -3*t28*t70*t71+4.2d-3*gammabb*t127*t28*t71+5.599999999999
     6         999d-3*t31*t65))*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t7*(-3.5461536035682834d-1*t10*t11*t
     1         13*t9*(9.996486610856369d+0*rhoa*t131*t19*t22*t45*t99+2.6
     2         59615202676212d-1*t10*t11*t23*t86*t87+2.659615202676212d-
     3         1*t10*t11*t13*(1.413716694115413d+1*t134*t15*t17*t19*t59-
     4         3.536776513153214d-2*t134*t14*t16*t20-5.0d-1*t134*t18*t19
     5         -3.534291735288533d+1*t108*t131*t15*t17*t19+1.99859489122
     6         06125d+2*t107*t131*t135*t19/t12**4)-9.421466669512572d+1*
     7         t130*t131*t19*t45*t57/t13**7-6.66432440723758d+0*rhoa*t13
     8         4*t19*t22*t45*t55-6.649038006690531d-2*t10*t11*t131*t21*t
     9         55+1.329807601338106d-1*t10*t11*t134*t21*t23)-3.546153603
     :         5682834d-1*t10*t11*t23*t86*t88*t9+8.865384008920708d-2*t1
     ;         0*t11*t131*t24*t55*t9-1.7730768017841417d-1*t10*t11*t134*
     <         t23*t24*t9)+2*t89*t90+t25*(8.4d-3*t3*t52*t85+4.2d-3*gamma
     =         aa*t133*t3*t52-8.4d-3*gammaaa*t102*t132*t3))*wght+Cmat2(i
     >         q,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (2*t96*t97+t32*(-3.5461536035682834d-
     1         1*t10*t33*t42*t9*t93*t95-3.5461536035682834d-1*t10*t33*t3
     2         5*t9*(2.659615202676212d-1*t10*t33*t42*t93*t94+2.65961520
     3         2676212d-1*t10*t33*t35*(1.413716694115413d+1*t139*t17*t36
     4         *t39*t78-3.536776513153214d-2*t139*t14*t37*t40-5.0d-1*t13
     5         9*t38*t39-3.534291735288533d+1*t121*t136*t17*t36*t39+1.99
     6         85948912206125d+2*t120*t135*t136*t39/t34**4)-9.4214666695
     7         12572d+1*t130*t136*t39*t45*t76/t35**7-6.66432440723758d+0
     8         *rhob*t139*t22*t39*t45*t74-6.649038006690531d-2*t10*t136*
     9         t33*t41*t74+9.996486610856369d+0*rhob*t112*t136*t22*t39*t
     :         45+1.329807601338106d-1*t10*t139*t33*t41*t42)+8.865384008
     ;         920708d-2*t10*t136*t33*t43*t74*t9-1.7730768017841417d-1*t
     <         10*t139*t33*t42*t43*t9)+t44*(8.4d-3*t28*t71*t92+4.2d-3*ga
     =         mmabb*t138*t28*t71-8.4d-3*gammabb*t115*t137*t28))*wght+Cm
     >         at2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = (gammaaa*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhoa**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t29*t30
            t32 = 1/t5**2
            t33 = rhoa**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammaaa*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         a*t21*t34*t35
            t37 = rhoa**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhoa*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammaaa*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammaaa*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammaaa*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = 1/t28**3
            t58 = -4.4799999999999995d-2*gammaaa**2*t57/rhoa**7.33333333
     1         3333333d+0+7.84d-2*t1*t4*t54+1.6799999999999998d-1*gammaa
     2         a*t29/rhoa**4.666666666666667d+0
            t59 = 1/rhoa**6.666666666666666d-1
            t60 = 1.3066666666666665d-2*gammaaa*t54*t6+4.135669939329334
     1         d-1*t59-4.2d-3*gammaaa*t3*t32*t58+8.4d-3*gammaaa*t3*t55*t
     2         56+1.1199999999999999d-2*gammaaa*t26*t31*t32
            t61 = rhoa**4
            t62 = 1/t11**3
            t63 = -2.827433388230826d+1*t14*t16*t52*t62+1.41371669411541
     1         3d+1*t14*t16*t39*t60+5.654866776461652d+1*rhoa*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t64 = t40**2
            t65 = 1.6799999999999998d-2*gammaaa*t57/rhoa**6.333333333333
     1         333d+0-1.6799999999999998d-2*t26*t4*t45-5.04d-2*t29*t30
            t66 = -4.2d-3*gammaaa*t3*t32*t65-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammaaa*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         aa*t26*t32*t46-4.2d-3*t3*t31*t32
            t67 = 1.413716694115413d+1*t14*t16*t39*t66-2.827433388230826
     1         d+1*t14*t16*t34*t47*t62+2.827433388230826d+1*rhoa*t16*t39
     2         *t47
            t68 = t47**2
            t69 = t46**2
            t70 = -6.3d-3*t57/rhoa**5.333333333333333d+0-6.3d-3*t3*t4/t1
     1         **3+6.3d-3*t27*t29/gammaaa
            t71 = -4.2d-3*gammaaa*t3*t32*t70+8.4d-3*gammaaa*t3*t56*t69-8
     1         .4d-3*t3*t32*t46
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = t43*t7*wght+t24*t44*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t50*t7*wght+t24*t51*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t7*(-3.5461536035682834d-1*t10*t12*t8*t
     1         9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t1
     2         1*t13*t15*t18*t64+t18*t64-3.536776513153214d-2*t11*t13*t1
     3         5*t18*t63+t18*t63-2.1220659078919285d-1*t11*t13*t19/t61-3
     4         .536776513153214d-2*t13*t15*t19*t60+1.4147106052612857d-1
     5         *t11*t13*t18*t38*t40-7.073553026306428d-2*t13*t15*t18*t34
     6         *t40+1.4147106052612857d-1*t13*t19*t34*t38)*t9+1.32980760
     7         1338106d-1*t10*t20*t22*t60*t9-6.649038006690531d-2*t10*t2
     8         0*t35*t52*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.31
     9         9230405352424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12*
     :         t20*t38*t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.5449
     ;         077018110403d+0*t18*t25*(-1.8799712059732546d+0*rhoa*t21*
     <         t35*t60+2.8199568089598814d+0*rhoa*t21*t52*t53-3.75994241
     =         1946509d+0*t21*t34*t35)+3.5449077018110403d+0*t18*t25*t36
     >         *t40)-1.7730768017841417d-1*t10*t22*t23*t60*t8*t9+8.86538
     ?         4008920708d-2*t10*t23*t35*t52*t8*t9-3.5461536035682834d-1
     @         *t10*t22*t34*t42*t8*t9+7.092307207136567d-1*t12*t15*t42*t
     1         8*t9-7.092307207136567d-1*t12*t23*t38*t8*t9+3.54615360356
     2         82834d-1*t15*t22*t23*t34*t8*t9)*wght+t24*(-1.306666666666
     3         6665d-2*gammaaa*t54*t6-4.135669939329334d-1*t59+4.2d-3*ga
     4         mmaaa*t3*t32*t58-8.4d-3*gammaaa*t3*t55*t56-1.119999999999
     5         9999d-2*gammaaa*t26*t31*t32)*wght+2*t43*t44*wght+Amat2(iq
     6         ,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t7*(-3.5461536035682834d-1*t10*t12*t8*
     1         t9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t
     2         11*t13*t15*t18*t67+t18*t67-3.536776513153214d-2*t13*t15*t
     3         19*t66+1.413716694115413d+1*t14*t16*t18*t39*t40*t47-5.0d-
     4         1*t17*t18*t40*t47-3.536776513153214d-2*t13*t15*t18*t40*t4
     5         7-5.0d-1*t18*t34*t39*t47+7.073553026306428d-2*t13*t19*t38
     6         *t47+1.0d+0*t10*t17*t18*t47)*t9+1.329807601338106d-1*t10*
     7         t20*t22*t66*t9+1.329807601338106d-1*t10*t22*t34*t48*t9-2.
     8         659615202676212d-1*t12*t15*t48*t9+1.329807601338106d-1*t1
     9         0*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*t35*t47
     :         *t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.54490770181
     ;         10403d+0*t18*t25*(-1.8799712059732546d+0*rhoa*t21*t35*t66
     <         +2.8199568089598814d+0*rhoa*t21*t34*t47*t53-1.87997120597
     =         32546d+0*t21*t35*t47)+5.01149519714857d+1*t14*t16*t18*t25
     >         *t36*t39*t47)-1.7730768017841417d-1*t10*t22*t23*t66*t8*t9
     ?         -1.7730768017841417d-1*t10*t22*t34*t49*t8*t9+3.5461536035
     @         682834d-1*t12*t15*t49*t8*t9-1.7730768017841417d-1*t10*t22
     1         *t42*t47*t8*t9+8.865384008920708d-2*t10*t23*t34*t35*t47*t
     2         8*t9+1.7730768017841417d-1*t15*t22*t23*t47*t8*t9)*wght+t2
     3         4*(4.2d-3*gammaaa*t3*t32*t65+5.599999999999999d-3*t26*t6-
     4         8.4d-3*gammaaa*t3*t31*t46*t56-5.599999999999999d-3*gammaa
     5         a*t26*t32*t46+4.2d-3*t3*t31*t32)*wght+t43*t51*wght+t44*t5
     6         0*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t7*(-3.5461536035682834d-1*t10*t12*t8
     1         *t9*(2.659615202676212d-1*t10*t12*(1.9985948912206125d+2*
     2         t18*t61*t68/(t11**4*t9**4)+1.413716694115413d+1*t14*t16*t
     3         18*t39*t71-3.536776513153214d-2*t13*t15*t19*t71-5.0d-1*t1
     4         7*t18*t71-3.534291735288533d+1*t14*t16*t18*t62*t68)*t9+1.
     5         329807601338106d-1*t10*t20*t22*t71*t9-6.649038006690531d-
     6         2*t10*t20*t35*t68*t9+2.659615202676212d-1*t10*t22*t47*t48
     7         *t9-9.421466669512572d+1*t18*t25*t37*t68/(t12**7*t9**3)-6
     8         .66432440723758d+0*rhoa*t18*t21*t25*t35*t71+9.99648661085
     9         6369d+0*rhoa*t18*t21*t25*t53*t68)-1.7730768017841417d-1*t
     :         10*t22*t23*t71*t8*t9+8.865384008920708d-2*t10*t23*t35*t68
     ;         *t8*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8*t9)*wght+
     <         t24*(4.2d-3*gammaaa*t3*t32*t70-8.4d-3*gammaaa*t3*t56*t69+
     =         8.4d-3*t3*t32*t46)*wght+2*t50*t51*wght+Cmat2(iq,D2_GAA_GA
     >         A)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = rhob**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammabb*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhob
            t11 = 4.2d-3*gammabb*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhob**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhob*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhob**2.3333333333333334d+0
            t27 = 1/rhob**2.6666666666666666d+0
            t28 = (gammabb*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhob**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammabb*t29*t30
            t32 = 1/t5**2
            t33 = rhob**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammabb*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammabb*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         b*t21*t34*t35
            t37 = rhob**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhob*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammabb*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammabb*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammabb*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhob*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammabb*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhob**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = 1/t28**3
            t58 = -4.4799999999999995d-2*gammabb**2*t57/rhob**7.33333333
     1         3333333d+0+7.84d-2*t1*t4*t54+1.6799999999999998d-1*gammab
     2         b*t29/rhob**4.666666666666667d+0
            t59 = 1/rhob**6.666666666666666d-1
            t60 = 1.3066666666666665d-2*gammabb*t54*t6+4.135669939329334
     1         d-1*t59-4.2d-3*gammabb*t3*t32*t58+8.4d-3*gammabb*t3*t55*t
     2         56+1.1199999999999999d-2*gammabb*t26*t31*t32
            t61 = rhob**4
            t62 = 1/t11**3
            t63 = -2.827433388230826d+1*t14*t16*t52*t62+1.41371669411541
     1         3d+1*t14*t16*t39*t60+5.654866776461652d+1*rhob*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t64 = t40**2
            t65 = 1.6799999999999998d-2*gammabb*t57/rhob**6.333333333333
     1         333d+0-1.6799999999999998d-2*t26*t4*t45-5.04d-2*t29*t30
            t66 = -4.2d-3*gammabb*t3*t32*t65-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammabb*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         bb*t26*t32*t46-4.2d-3*t3*t31*t32
            t67 = 1.413716694115413d+1*t14*t16*t39*t66-2.827433388230826
     1         d+1*t14*t16*t34*t47*t62+2.827433388230826d+1*rhob*t16*t39
     2         *t47
            t68 = t47**2
            t69 = t46**2
            t70 = -6.3d-3*t57/rhob**5.333333333333333d+0-6.3d-3*t3*t4/t1
     1         **3+6.3d-3*t27*t29/gammabb
            t71 = -4.2d-3*gammabb*t3*t32*t70+8.4d-3*gammabb*t3*t56*t69-8
     1         .4d-3*t3*t32*t46
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RB) = t43*t7*wght+t24*t44*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t50*t7*wght+t24*t51*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t7*(-3.5461536035682834d-1*t10*t12*t8*t
     1         9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t1
     2         1*t13*t15*t18*t64+t18*t64-3.536776513153214d-2*t11*t13*t1
     3         5*t18*t63+t18*t63-2.1220659078919285d-1*t11*t13*t19/t61-3
     4         .536776513153214d-2*t13*t15*t19*t60+1.4147106052612857d-1
     5         *t11*t13*t18*t38*t40-7.073553026306428d-2*t13*t15*t18*t34
     6         *t40+1.4147106052612857d-1*t13*t19*t34*t38)*t9+1.32980760
     7         1338106d-1*t10*t20*t22*t60*t9-6.649038006690531d-2*t10*t2
     8         0*t35*t52*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.31
     9         9230405352424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12*
     :         t20*t38*t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.5449
     ;         077018110403d+0*t18*t25*(-1.8799712059732546d+0*rhob*t21*
     <         t35*t60+2.8199568089598814d+0*rhob*t21*t52*t53-3.75994241
     =         1946509d+0*t21*t34*t35)+3.5449077018110403d+0*t18*t25*t36
     >         *t40)-1.7730768017841417d-1*t10*t22*t23*t60*t8*t9+8.86538
     ?         4008920708d-2*t10*t23*t35*t52*t8*t9-3.5461536035682834d-1
     @         *t10*t22*t34*t42*t8*t9+7.092307207136567d-1*t12*t15*t42*t
     1         8*t9-7.092307207136567d-1*t12*t23*t38*t8*t9+3.54615360356
     2         82834d-1*t15*t22*t23*t34*t8*t9)*wght+t24*(-1.306666666666
     3         6665d-2*gammabb*t54*t6-4.135669939329334d-1*t59+4.2d-3*ga
     4         mmabb*t3*t32*t58-8.4d-3*gammabb*t3*t55*t56-1.119999999999
     5         9999d-2*gammabb*t26*t31*t32)*wght+2*t43*t44*wght+Amat2(iq
     6         ,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t7*(-3.5461536035682834d-1*t10*t12*t8*
     1         t9*(2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t
     2         11*t13*t15*t18*t67+t18*t67-3.536776513153214d-2*t13*t15*t
     3         19*t66+1.413716694115413d+1*t14*t16*t18*t39*t40*t47-5.0d-
     4         1*t17*t18*t40*t47-3.536776513153214d-2*t13*t15*t18*t40*t4
     5         7-5.0d-1*t18*t34*t39*t47+7.073553026306428d-2*t13*t19*t38
     6         *t47+1.0d+0*t10*t17*t18*t47)*t9+1.329807601338106d-1*t10*
     7         t20*t22*t66*t9+1.329807601338106d-1*t10*t22*t34*t48*t9-2.
     8         659615202676212d-1*t12*t15*t48*t9+1.329807601338106d-1*t1
     9         0*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*t35*t47
     :         *t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.54490770181
     ;         10403d+0*t18*t25*(-1.8799712059732546d+0*rhob*t21*t35*t66
     <         +2.8199568089598814d+0*rhob*t21*t34*t47*t53-1.87997120597
     =         32546d+0*t21*t35*t47)+5.01149519714857d+1*t14*t16*t18*t25
     >         *t36*t39*t47)-1.7730768017841417d-1*t10*t22*t23*t66*t8*t9
     ?         -1.7730768017841417d-1*t10*t22*t34*t49*t8*t9+3.5461536035
     @         682834d-1*t12*t15*t49*t8*t9-1.7730768017841417d-1*t10*t22
     1         *t42*t47*t8*t9+8.865384008920708d-2*t10*t23*t34*t35*t47*t
     2         8*t9+1.7730768017841417d-1*t15*t22*t23*t47*t8*t9)*wght+t2
     3         4*(4.2d-3*gammabb*t3*t32*t65+5.599999999999999d-3*t26*t6-
     4         8.4d-3*gammabb*t3*t31*t46*t56-5.599999999999999d-3*gammab
     5         b*t26*t32*t46+4.2d-3*t3*t31*t32)*wght+t43*t51*wght+t44*t5
     6         0*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t7*(-3.5461536035682834d-1*t10*t12*t8
     1         *t9*(2.659615202676212d-1*t10*t12*(1.9985948912206125d+2*
     2         t18*t61*t68/(t11**4*t9**4)+1.413716694115413d+1*t14*t16*t
     3         18*t39*t71-3.536776513153214d-2*t13*t15*t19*t71-5.0d-1*t1
     4         7*t18*t71-3.534291735288533d+1*t14*t16*t18*t62*t68)*t9+1.
     5         329807601338106d-1*t10*t20*t22*t71*t9-6.649038006690531d-
     6         2*t10*t20*t35*t68*t9+2.659615202676212d-1*t10*t22*t47*t48
     7         *t9-9.421466669512572d+1*t18*t25*t37*t68/(t12**7*t9**3)-6
     8         .66432440723758d+0*rhob*t18*t21*t25*t35*t71+9.99648661085
     9         6369d+0*rhob*t18*t21*t25*t53*t68)-1.7730768017841417d-1*t
     :         10*t22*t23*t71*t8*t9+8.865384008920708d-2*t10*t23*t35*t68
     ;         *t8*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8*t9)*wght+
     <         t24*(4.2d-3*gammabb*t3*t32*t70-8.4d-3*gammabb*t3*t56*t69+
     =         8.4d-3*t3*t32*t46)*wght+2*t50*t51*wght+Cmat2(iq,D2_GBB_GB
     >         B)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_camb88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; 
C>    \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}
C>    +0.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{
C>    \it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}
C>    \,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{
C>    \it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}+0.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}
C>    \,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{
C>    \it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    +0.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{
C>    \it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}-0.9305257363491002\,{\it t_{16}}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{
C>    \it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552
C>    \,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}
C>    \,\rho_\beta}\over{{\it t_{22}}}}\right)+0.2659615202676212
C>    \,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({
C>    \it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}
C>    \,\sigma_{\alpha\alpha}-0.9305257363491002\,{\it t_1}\right)
C>    \,\left(-0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{
C>    \it t_{10}}}}\right)+0.2659615202676212\,{\it t_7}\,{\it t_8}
C>    \,{\it t_{10}}\,\left({\it t_{14}}-{{0.03536776513153214\,{
C>    \it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}
C>    \over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}
C>    -0.9305257363491002\,{\it t_{25}}\right)\,\left(
C>    -0.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{
C>    \it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{
C>    \it t_{31}}}}\right)+0.2659615202676212\,{\it t_7}\,{
C>    \it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}
C>    -{{0.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({
C>    \it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)
C>    +{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  
C>    , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_camb88_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = (gammaaa*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhoa**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t29*t30
            t32 = 1/t5**2
            t33 = rhoa**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammaaa*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         a*t21*t34*t35
            t37 = rhoa**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhoa*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammaaa*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammaaa*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammaaa*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = gammaaa**2
            t58 = 1/t28**3
            t59 = 1/rhoa**7.333333333333333d+0
            t60 = 1/rhoa**4.666666666666667d+0
            t61 = 1.6799999999999998d-1*gammaaa*t29*t60-4.47999999999999
     1         95d-2*t57*t58*t59+7.84d-2*t1*t4*t54
            t62 = 1/rhoa**6.666666666666666d-1
            t63 = 4.135669939329334d-1*t62-4.2d-3*gammaaa*t3*t32*t61+1.3
     1         066666666666665d-2*gammaaa*t54*t6+8.4d-3*gammaaa*t3*t55*t
     2         56+1.1199999999999999d-2*gammaaa*t26*t31*t32
            t64 = -1.8799712059732546d+0*rhoa*t21*t35*t63+2.819956808959
     1         8814d+0*rhoa*t21*t52*t53-3.759942411946509d+0*t21*t34*t35
            t65 = rhoa**4
            t66 = 1/t65
            t67 = 1/t11**3
            t68 = -2.827433388230826d+1*t14*t16*t52*t67+1.41371669411541
     1         3d+1*t14*t16*t39*t63+5.654866776461652d+1*rhoa*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t69 = t40**2
            t70 = -3.536776513153214d-2*t11*t13*t15*t18*t69+t18*t69-3.53
     1         6776513153214d-2*t11*t13*t15*t18*t68+t18*t68-2.1220659078
     2         919285d-1*t11*t13*t19*t66-3.536776513153214d-2*t13*t15*t1
     3         9*t63+1.4147106052612857d-1*t11*t13*t18*t38*t40-7.0735530
     4         26306428d-2*t13*t15*t18*t34*t40+1.4147106052612857d-1*t13
     5         *t19*t34*t38
            t71 = 2.659615202676212d-1*t10*t12*t70*t9+1.329807601338106d
     1         -1*t10*t20*t22*t63*t9-6.649038006690531d-2*t10*t20*t35*t5
     2         2*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.3192304053
     3         52424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12*t20*t38*
     4         t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.544907701811
     5         0403d+0*t18*t25*t64+3.5449077018110403d+0*t18*t25*t36*t40
            t72 = -3.5461536035682834d-1*t10*t12*t71*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t63*t8*t9+8.865384008920708d-2*t10*t
     2         23*t35*t52*t8*t9-3.5461536035682834d-1*t10*t22*t34*t42*t8
     3         *t9+7.092307207136567d-1*t12*t15*t42*t8*t9-7.092307207136
     4         567d-1*t12*t23*t38*t8*t9+3.5461536035682834d-1*t15*t22*t2
     5         3*t34*t8*t9
            t73 = -4.135669939329334d-1*t62+4.2d-3*gammaaa*t3*t32*t61-1.
     1         3066666666666665d-2*gammaaa*t54*t6-8.4d-3*gammaaa*t3*t55*
     2         t56-1.1199999999999999d-2*gammaaa*t26*t31*t32
            t74 = 1/rhoa**6.333333333333333d+0
            t75 = 1.6799999999999998d-2*gammaaa*t58*t74-1.67999999999999
     1         98d-2*t26*t4*t45-5.04d-2*t29*t30
            t76 = -4.2d-3*gammaaa*t3*t32*t75-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammaaa*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         aa*t26*t32*t46-4.2d-3*t3*t31*t32
            t77 = -1.8799712059732546d+0*rhoa*t21*t35*t76+2.819956808959
     1         8814d+0*rhoa*t21*t34*t47*t53-1.8799712059732546d+0*t21*t3
     2         5*t47
            t78 = 1.413716694115413d+1*t14*t16*t39*t76-2.827433388230826
     1         d+1*t14*t16*t34*t47*t67+2.827433388230826d+1*rhoa*t16*t39
     2         *t47
            t79 = -3.536776513153214d-2*t11*t13*t15*t18*t78+t18*t78-3.53
     1         6776513153214d-2*t13*t15*t19*t76+1.413716694115413d+1*t14
     2         *t16*t18*t39*t40*t47-5.0d-1*t17*t18*t40*t47-3.53677651315
     3         3214d-2*t13*t15*t18*t40*t47-5.0d-1*t18*t34*t39*t47+7.0735
     4         53026306428d-2*t13*t19*t38*t47+1.0d+0*t10*t17*t18*t47
            t80 = 2.659615202676212d-1*t10*t12*t79*t9+1.329807601338106d
     1         -1*t10*t20*t22*t76*t9+1.329807601338106d-1*t10*t22*t34*t4
     2         8*t9-2.659615202676212d-1*t12*t15*t48*t9+1.32980760133810
     3         6d-1*t10*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*
     4         t35*t47*t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.5449
     5         077018110403d+0*t18*t25*t77+5.01149519714857d+1*t14*t16*t
     6         18*t25*t36*t39*t47
            t81 = -3.5461536035682834d-1*t10*t12*t8*t80*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*
     2         t22*t34*t49*t8*t9+3.5461536035682834d-1*t12*t15*t49*t8*t9
     3         -1.7730768017841417d-1*t10*t22*t42*t47*t8*t9+8.8653840089
     4         20708d-2*t10*t23*t34*t35*t47*t8*t9+1.7730768017841417d-1*
     5         t15*t22*t23*t47*t8*t9
            t82 = 4.2d-3*gammaaa*t3*t32*t75+5.599999999999999d-3*t26*t6-
     1         8.4d-3*gammaaa*t3*t31*t46*t56-5.599999999999999d-3*gammaa
     2         a*t26*t32*t46+4.2d-3*t3*t31*t32
            t83 = 1/t9**3
            t84 = t47**2
            t85 = 1/t12**7
            t86 = t46**2
            t87 = 1/rhoa**5.333333333333333d+0
            t88 = 1/gammaaa
            t89 = 1/t1**3
            t90 = -6.3d-3*t3*t4*t89+6.3d-3*t27*t29*t88-6.3d-3*t58*t87
            t91 = -4.2d-3*gammaaa*t3*t32*t90+8.4d-3*gammaaa*t3*t56*t86-8
     1         .4d-3*t3*t32*t46
            t92 = 1/t9**4
            t93 = 1/t11**4
            t94 = 1.9985948912206125d+2*t18*t65*t84*t92*t93+1.4137166941
     1         15413d+1*t14*t16*t18*t39*t91-3.536776513153214d-2*t13*t15
     2         *t19*t91-5.0d-1*t17*t18*t91-3.534291735288533d+1*t14*t16*
     3         t18*t67*t84
            t95 = 2.659615202676212d-1*t10*t12*t9*t94+1.329807601338106d
     1         -1*t10*t20*t22*t9*t91-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t91-6.649038006690531d-2*t10*t20*t35*t84*t9+2.65961
     3         5202676212d-1*t10*t22*t47*t48*t9-9.421466669512572d+1*t18
     4         *t25*t37*t83*t84*t85+9.996486610856369d+0*rhoa*t18*t21*t2
     5         5*t53*t84
            t96 = -3.5461536035682834d-1*t10*t12*t8*t9*t95-1.77307680178
     1         41417d-1*t10*t22*t23*t8*t9*t91+8.865384008920708d-2*t10*t
     2         23*t35*t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8
     3         *t9
            t97 = 4.2d-3*gammaaa*t3*t32*t90-8.4d-3*gammaaa*t3*t56*t86+8.
     1         4d-3*t3*t32*t46
            t98 = t34**3
            t99 = 1/rhoa**4.333333333333333d+0
            t100 = 1/rhoa**1.6666666666666669d+0
            t101 = t31**3
            t102 = 1/t5**4
            t103 = 1/t28**5
            t104 = -2.6133333333333336d-1*t1*t4*t99+5.525333333333333d-1
     1         *t57*t58/rhoa**8.333333333333334d+0-8.885333333333333d-1*
     2         gammaaa*t29/rhoa**5.666666666666667d+0-1.7919999999999997
     3         d-1*gammaaa**3*t103/rhoa**11
            t105 = -4.355555555555556d-2*gammaaa*t6*t99+2.52d-2*gammaaa*
     1         t3*t31*t56*t61+1.6799999999999995d-2*gammaaa*t26*t32*t61-
     2         3.359999999999999d-2*gammaaa*t26*t55*t56-3.92d-2*gammaaa*
     3         t31*t32*t54-4.2d-3*gammaaa*t104*t3*t32-2.52d-2*gammaaa*t1
     4         01*t102*t3-2.7571132928862224d-1*t100
            t106 = rhoa**5
            t107 = 8.482300164692478d+1*t14*t16*t93*t98-8.48230016469247
     1         8d+1*t14*t16*t34*t63*t67-1.6964600329384957d+2*rhoa*t16*t
     2         52*t67+8.482300164692478d+1*rhoa*t16*t39*t63+8.4823001646
     3         92478d+1*t16*t34*t39+1.413716694115413d+1*t105*t14*t16*t3
     4         9
            t108 = t40**3
            t109 = 2.072d-1*t29*t60-1.7359999999999998d-1*gammaaa*t58*t5
     1         9+6.719999999999998d-2*t103*t57/rhoa**10+3.92d-2*t4*t45*t
     2         54
            t110 = 1.6799999999999998d-2*gammaaa*t3*t31*t56*t75+1.119999
     1         9999999999d-2*gammaaa*t26*t32*t75+8.4d-3*gammaaa*t3*t46*t
     2         56*t61-4.2d-3*t3*t32*t61+1.3066666666666665d-2*t54*t6+8.4
     3         d-3*t3*t55*t56-2.2399999999999998d-2*gammaaa*t26*t31*t46*
     4         t56-2.52d-2*gammaaa*t102*t3*t46*t55-1.3066666666666665d-2
     5         *gammaaa*t32*t46*t54+1.1199999999999999d-2*t26*t31*t32-4.
     6         2d-3*gammaaa*t109*t3*t32
            t111 = 8.482300164692478d+1*t14*t16*t47*t52*t93-5.6548667764
     1         61652d+1*t14*t16*t34*t67*t76+5.654866776461652d+1*rhoa*t1
     2         6*t39*t76-2.827433388230826d+1*t14*t16*t47*t63*t67-1.1309
     3         733552923305d+2*rhoa*t16*t34*t47*t67+2.827433388230826d+1
     4         *t16*t39*t47+1.413716694115413d+1*t110*t14*t16*t39
            t112 = 8.4d-3*t26*t4*t89-8.4d-3*t29*t30*t88+4.19999999999999
     1         93d-2*t58*t74-2.52d-2*gammaaa*t103/rhoa**9
            t113 = 8.4d-3*gammaaa*t3*t31*t56*t90+5.599999999999999d-3*ga
     1         mmaaa*t26*t32*t90-1.1199999999999999d-2*gammaaa*t26*t56*t
     2         86-2.52d-2*gammaaa*t102*t3*t31*t86+1.6799999999999998d-2*
     3         gammaaa*t3*t46*t56*t75-8.4d-3*t3*t32*t75+1.67999999999999
     4         98d-2*t3*t31*t46*t56+1.1199999999999999d-2*t26*t32*t46-4.
     5         2d-3*gammaaa*t112*t3*t32
            t114 = 8.482300164692478d+1*t14*t16*t34*t84*t93-2.8274333882
     1         30826d+1*t14*t16*t34*t67*t91+2.827433388230826d+1*rhoa*t1
     2         6*t39*t91-5.654866776461652d+1*rhoa*t16*t67*t84-5.6548667
     3         76461652d+1*t14*t16*t47*t67*t76+1.413716694115413d+1*t113
     4         *t14*t16*t39
            t115 = t47**3
            t116 = t46**3
            t117 = -3.15d-3*t58*t87*t88-9.450000000000001d-3*t27*t29/t57
     1         +9.450000000000001d-3*t3*t4/t1**5+9.450000000000001d-3*t1
     2         03/rhoa**8
            t118 = 2.52d-2*gammaaa*t3*t46*t56*t90-1.26d-2*t3*t32*t90+2.5
     1         2d-2*t3*t56*t86-4.2d-3*gammaaa*t117*t3*t32-2.52d-2*gammaa
     2         a*t102*t116*t3
            fnc(iq) = 2*t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = (t43*t7+t24*t44)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t50*t7+t24*t51)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t24*t73+t7*t72+2*t43*t44)*wght+Amat2(i
     1         q,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t24*t82+t7*t81+t43*t51+t44*t50)*wght+
     1         Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t24*t97+t7*t96+2*t50*t51)*wght+Cmat2
     1         (iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (t24*(4.355555555555556d-2*gammaaa*t
     1         6*t99-2.52d-2*gammaaa*t3*t31*t56*t61-1.6799999999999995d-
     2         2*gammaaa*t26*t32*t61+3.359999999999999d-2*gammaaa*t26*t5
     3         5*t56+3.92d-2*gammaaa*t31*t32*t54+4.2d-3*gammaaa*t104*t3*
     4         t32+2.52d-2*gammaaa*t101*t102*t3+2.7571132928862224d-1*t1
     5         00)+t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(3.544907701
     6         8110403d+0*t18*t25*(-7.049892022399704d+0*rhoa*t21*t85*t9
     7         8+8.459870426879645d+0*rhoa*t21*t34*t53*t63-5.63991361791
     8         9763d+0*t21*t35*t63+8.459870426879645d+0*t21*t52*t53-1.87
     9         99712059732546d+0*rhoa*t105*t21*t35)+9.973557010035795d-2
     :         *t10*t20*t53*t9*t98+3.989422804014318d-1*t10*t22*t34*t70*
     ;         t9-7.978845608028636d-1*t12*t15*t70*t9+2.659615202676212d
     <         -1*t10*t12*(2.1220659078919285d-1*t11*t13*t18*t38*t69-1.0
     =         610329539459643d-1*t13*t15*t18*t34*t69-1.0610329539459643
     >         d-1*t11*t13*t15*t18*t40*t68+3*t18*t40*t68+2.1220659078919
     ?         285d-1*t11*t13*t18*t38*t68-1.0610329539459643d-1*t13*t15*
     @         t18*t34*t68-6.366197723675786d-1*t11*t13*t18*t40*t66-6.36
     1         6197723675786d-1*t13*t19*t34*t66-1.0610329539459643d-1*t1
     2         3*t15*t18*t40*t63+2.1220659078919285d-1*t13*t19*t38*t63+4
     3         .244131815783857d-1*t13*t18*t34*t38*t40-3.536776513153214
     4         d-2*t105*t13*t15*t19+8.488263631567714d-1*t11*t13*t19/t10
     5         6-3.536776513153214d-2*t108*t11*t13*t15*t18-3.53677651315
     6         3214d-2*t107*t11*t13*t15*t18+t108*t18+t107*t18)*t9-1.5957
     7         691216057274d+0*t12*t20*t66*t9+3.989422804014318d-1*t10*t
     8         22*t41*t63*t9-1.994711402007159d-1*t10*t20*t34*t35*t63*t9
     9         -3.989422804014318d-1*t15*t20*t22*t63*t9-1.99471140200715
     :         9d-1*t10*t35*t41*t52*t9+1.994711402007159d-1*t15*t20*t35*
     ;         t52*t9+1.5957691216057274d+0*t12*t38*t41*t9-7.97884560802
     <         8636d-1*t15*t22*t34*t41*t9+7.978845608028636d-1*t20*t22*t
     =         34*t38*t9+1.329807601338106d-1*t10*t105*t20*t22*t9+3.5449
     >         077018110403d+0*t18*t25*t36*t69+3.5449077018110403d+0*t18
     ?         *t25*t36*t68+7.089815403622081d+0*t18*t25*t40*t64)-1.3298
     @         076013381063d-1*t10*t23*t53*t8*t9*t98-5.319230405352425d-
     1         1*t10*t22*t34*t71*t8*t9+1.063846081070485d+0*t12*t15*t71*
     2         t8*t9+2.12769216214097d+0*t12*t23*t66*t8*t9-5.31923040535
     3         2425d-1*t10*t22*t42*t63*t8*t9+2.6596152026762127d-1*t10*t
     4         23*t34*t35*t63*t8*t9+5.319230405352425d-1*t15*t22*t23*t63
     5         *t8*t9+2.6596152026762127d-1*t10*t35*t42*t52*t8*t9-2.6596
     6         152026762127d-1*t15*t23*t35*t52*t8*t9-2.12769216214097d+0
     7         *t12*t38*t42*t8*t9+1.063846081070485d+0*t15*t22*t34*t42*t
     8         8*t9-1.063846081070485d+0*t22*t23*t34*t38*t8*t9-1.7730768
     9         017841417d-1*t10*t105*t22*t23*t8*t9)+3*t43*t73+3*t44*t72)
     :         *wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t7*(-3.5461536035682834d-1*t10*t12
     1         *t8*t9*(2.659615202676212d-1*t10*t22*t34*t79*t9-5.3192304
     2         05352424d-1*t12*t15*t79*t9+2.659615202676212d-1*t10*t12*(
     3         -7.073553026306428d-2*t11*t13*t15*t18*t40*t78+2*t18*t40*t
     4         78+1.4147106052612857d-1*t11*t13*t18*t38*t78-7.0735530263
     5         06428d-2*t13*t15*t18*t34*t78-7.073553026306428d-2*t13*t15
     6         *t18*t40*t76+1.4147106052612857d-1*t13*t19*t38*t76+1.4137
     7         16694115413d+1*t14*t16*t18*t39*t47*t69-5.0d-1*t17*t18*t47
     8         *t69-3.536776513153214d-2*t13*t15*t18*t47*t69+1.413716694
     9         115413d+1*t14*t16*t18*t39*t47*t68-5.0d-1*t17*t18*t47*t68-
     :         3.536776513153214d-2*t13*t15*t18*t47*t68-2.12206590789192
     ;         85d-1*t13*t19*t47*t66-5.0d-1*t18*t39*t47*t63-1.0d+0*t18*t
     <         34*t39*t40*t47+1.4147106052612857d-1*t13*t18*t38*t40*t47+
     =         2.0d+0*t10*t17*t18*t40*t47+2.0d+0*t10*t18*t34*t39*t47-3.0
     >         d+0*t15*t17*t18*t47-3.536776513153214d-2*t110*t13*t15*t19
     ?         -3.536776513153214d-2*t11*t111*t13*t15*t18+t111*t18)*t9+2
     @         .659615202676212d-1*t10*t22*t41*t76*t9-1.329807601338106d
     1         -1*t10*t20*t34*t35*t76*t9-2.659615202676212d-1*t15*t20*t2
     2         2*t76*t9+1.329807601338106d-1*t10*t22*t47*t70*t9+1.329807
     3         601338106d-1*t10*t22*t48*t63*t9-6.649038006690531d-2*t10*
     4         t20*t35*t47*t63*t9+9.973557010035795d-2*t10*t20*t47*t52*t
     5         53*t9-6.649038006690531d-2*t10*t35*t48*t52*t9+5.319230405
     6         352424d-1*t12*t38*t48*t9-2.659615202676212d-1*t15*t22*t34
     7         *t48*t9-1.329807601338106d-1*t10*t34*t35*t41*t47*t9-2.659
     8         615202676212d-1*t15*t22*t41*t47*t9+2.659615202676212d-1*t
     9         20*t22*t38*t47*t9+1.329807601338106d-1*t15*t20*t34*t35*t4
     :         7*t9+1.329807601338106d-1*t10*t110*t20*t22*t9+3.544907701
     ;         8110403d+0*t18*t25*(-7.049892022399704d+0*rhoa*t21*t47*t5
     <         2*t85+5.639913617919763d+0*rhoa*t21*t34*t53*t76-3.7599424
     =         11946509d+0*t21*t35*t76+2.8199568089598814d+0*rhoa*t21*t4
     >         7*t53*t63+5.639913617919763d+0*t21*t34*t47*t53-1.87997120
     ?         59732546d+0*rhoa*t110*t21*t35)+3.5449077018110403d+0*t18*
     @         t25*t36*t78+3.5449077018110403d+0*t18*t25*t40*t77+5.01149
     1         519714857d+1*t14*t16*t18*t25*t39*t47*t64+5.01149519714857
     2         d+1*t14*t16*t18*t25*t36*t39*t40*t47)-3.5461536035682834d-
     3         1*t10*t22*t34*t8*t80*t9+7.092307207136567d-1*t12*t15*t8*t
     4         80*t9-3.5461536035682834d-1*t10*t22*t42*t76*t8*t9+1.77307
     5         68017841417d-1*t10*t23*t34*t35*t76*t8*t9+3.54615360356828
     6         34d-1*t15*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*t22
     7         *t47*t71*t8*t9-1.7730768017841417d-1*t10*t22*t49*t63*t8*t
     8         9+8.865384008920708d-2*t10*t23*t35*t47*t63*t8*t9-1.329807
     9         6013381063d-1*t10*t23*t47*t52*t53*t8*t9+8.865384008920708
     :         d-2*t10*t35*t49*t52*t8*t9-7.092307207136567d-1*t12*t38*t4
     ;         9*t8*t9+3.5461536035682834d-1*t15*t22*t34*t49*t8*t9+1.773
     <         0768017841417d-1*t10*t34*t35*t42*t47*t8*t9+3.546153603568
     =         2834d-1*t15*t22*t42*t47*t8*t9-3.5461536035682834d-1*t22*t
     >         23*t38*t47*t8*t9-1.7730768017841417d-1*t15*t23*t34*t35*t4
     ?         7*t8*t9-1.7730768017841417d-1*t10*t110*t22*t23*t8*t9)+2*t
     @         43*t82+2*t44*t81+t24*(-1.6799999999999998d-2*gammaaa*t3*t
     1         31*t56*t75-1.1199999999999999d-2*gammaaa*t26*t32*t75-8.4d
     2         -3*gammaaa*t3*t46*t56*t61+4.2d-3*t3*t32*t61-1.30666666666
     3         66665d-2*t54*t6-8.4d-3*t3*t55*t56+2.2399999999999998d-2*g
     4         ammaaa*t26*t31*t46*t56+2.52d-2*gammaaa*t102*t3*t46*t55+1.
     5         3066666666666665d-2*gammaaa*t32*t46*t54-1.119999999999999
     6         9d-2*t26*t31*t32+4.2d-3*gammaaa*t109*t3*t32)+t50*t73+t51*
     7         t72)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t43*t97+t44*t96+t7*(-1.7730768017
     1         841417d-1*t10*t22*t34*t8*t9*t95+3.5461536035682834d-1*t12
     2         *t15*t8*t9*t95-3.5461536035682834d-1*t10*t12*t8*t9*(1.329
     3         807601338106d-1*t10*t22*t34*t9*t94-2.659615202676212d-1*t
     4         12*t15*t9*t94+2.659615202676212d-1*t10*t12*t9*(1.99859489
     5         12206125d+2*t18*t40*t65*t84*t92*t93-7.068583470577066d+0*
     6         t14*t16*t18*t34*t84*t93+1.413716694115413d+1*t14*t16*t18*
     7         t39*t40*t91-5.0d-1*t17*t18*t40*t91-3.536776513153214d-2*t
     8         13*t15*t18*t40*t91-5.0d-1*t18*t34*t39*t91+7.0735530263064
     9         28d-2*t13*t19*t38*t91+1.0d+0*t10*t17*t18*t91-3.5342917352
     :         88533d+1*t14*t16*t18*t40*t67*t84+1.0d+0*t18*t34*t67*t84+1
     ;         .413716694115413d+1*rhoa*t16*t18*t67*t84+2.82743338823082
     <         6d+1*t14*t16*t18*t39*t47*t78-1.0d+0*t17*t18*t47*t78-7.073
     =         553026306428d-2*t13*t15*t18*t47*t78-1.0d+0*t18*t39*t47*t7
     >         6-3.536776513153214d-2*t113*t13*t15*t19-3.536776513153214
     ?         d-2*t11*t114*t13*t15*t18+t114*t18)+7.084834422688148d+2*t
     @         18*t25*t36*t65*t84*t92*t93+3.5449077018110403d+0*t18*t25*
     1         (2.8199568089598814d+0*rhoa*t21*t34*t53*t91-1.87997120597
     2         32546d+0*t21*t35*t91-7.049892022399704d+0*rhoa*t21*t34*t8
     3         4*t85+2.8199568089598814d+0*t21*t53*t84+5.639913617919763
     4         d+0*rhoa*t21*t47*t53*t76-1.8799712059732546d+0*rhoa*t113*
     5         t21*t35)+1.329807601338106d-1*t10*t22*t41*t9*t91-6.649038
     6         006690531d-2*t10*t20*t34*t35*t9*t91-1.329807601338106d-1*
     7         t15*t20*t22*t9*t91+5.01149519714857d+1*t14*t16*t18*t25*t3
     8         6*t39*t91+9.973557010035795d-2*t10*t20*t34*t53*t84*t9-6.6
     9         49038006690531d-2*t10*t35*t41*t84*t9+6.649038006690531d-2
     :         *t15*t20*t35*t84*t9+2.659615202676212d-1*t10*t22*t47*t79*
     ;         t9+2.659615202676212d-1*t10*t22*t48*t76*t9-1.329807601338
     <         106d-1*t10*t20*t35*t47*t76*t9-1.329807601338106d-1*t10*t3
     =         4*t35*t47*t48*t9-2.659615202676212d-1*t15*t22*t47*t48*t9+
     >         1.329807601338106d-1*t10*t113*t20*t22*t9-1.00229903942971
     ?         4d+2*t14*t16*t18*t25*t36*t67*t84+1.002299039429714d+2*t14
     @         *t16*t18*t25*t39*t47*t77)-1.7730768017841417d-1*t10*t22*t
     1         42*t8*t9*t91+8.865384008920708d-2*t10*t23*t34*t35*t8*t9*t
     2         91+1.7730768017841417d-1*t15*t22*t23*t8*t9*t91-1.32980760
     3         13381063d-1*t10*t23*t34*t53*t8*t84*t9+8.865384008920708d-
     4         2*t10*t35*t42*t8*t84*t9-8.865384008920708d-2*t15*t23*t35*
     5         t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t8*t80*t9-3.5
     6         461536035682834d-1*t10*t22*t49*t76*t8*t9+1.77307680178414
     7         17d-1*t10*t23*t35*t47*t76*t8*t9+1.7730768017841417d-1*t10
     8         *t34*t35*t47*t49*t8*t9+3.5461536035682834d-1*t15*t22*t47*
     9         t49*t8*t9-1.7730768017841417d-1*t10*t113*t22*t23*t8*t9)+t
     :         24*(-8.4d-3*gammaaa*t3*t31*t56*t90-5.599999999999999d-3*g
     ;         ammaaa*t26*t32*t90+1.1199999999999999d-2*gammaaa*t26*t56*
     <         t86+2.52d-2*gammaaa*t102*t3*t31*t86-1.6799999999999998d-2
     =         *gammaaa*t3*t46*t56*t75+8.4d-3*t3*t32*t75-1.6799999999999
     >         998d-2*t3*t31*t46*t56-1.1199999999999999d-2*t26*t32*t46+4
     ?         .2d-3*gammaaa*t112*t3*t32)+2*t50*t82+2*t51*t81)*wght+Cmat
     @         3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (3*t50*t97+3*t51*t96+t7*(-5.31923
     1         0405352425d-1*t10*t22*t47*t8*t9*t95-3.5461536035682834d-1
     2         *t10*t12*t8*t9*(3.989422804014318d-1*t10*t22*t47*t9*t94+2
     3         .659615202676212d-1*t10*t12*t9*(5.995784673661838d+2*t18*
     4         t47*t65*t91*t92*t93+1.0602875205865599d+2*t115*t14*t16*t1
     5         8*t93-1.2990866792933983d+3*t115*t18*t65*t92/t11**5-1.060
     6         2875205865599d+2*t14*t16*t18*t47*t67*t91+2.82544696249235
     7         84d+3*rhoa**6*t115*t18/(t11**6*t9**6)+1.413716694115413d+
     8         1*t118*t14*t16*t18*t39-3.536776513153214d-2*t118*t13*t15*
     9         t19-5.0d-1*t118*t17*t18)+3.989422804014318d-1*t10*t22*t48
     :         *t9*t91-1.994711402007159d-1*t10*t20*t35*t47*t9*t91-2.826
     ;         440000853771d+2*t18*t25*t37*t47*t83*t85*t91+2.99894598325
     <         6911d+1*rhoa*t18*t21*t25*t47*t53*t91-1.994711402007159d-1
     =         *t10*t35*t48*t84*t9+9.973557010035795d-2*t10*t115*t20*t53
     >         *t9+1.329807601338106d-1*t10*t118*t20*t22*t9-1.3319284713
     ?         741864d+3*t106*t115*t18*t25/(t12**11*t9**5)-2.49912165271
     @         40922d+1*rhoa*t115*t18*t21*t25*t85+4.710733334756285d+2*t
     1         115*t18*t25*t37*t83/t12**9-6.66432440723758d+0*rhoa*t118*
     2         t18*t21*t25*t35)-5.319230405352425d-1*t10*t22*t49*t8*t9*t
     3         91+2.6596152026762127d-1*t10*t23*t35*t47*t8*t9*t91+2.6596
     4         152026762127d-1*t10*t35*t49*t8*t84*t9-1.3298076013381063d
     5         -1*t10*t115*t23*t53*t8*t9-1.7730768017841417d-1*t10*t118*
     6         t22*t23*t8*t9)+t24*(-2.52d-2*gammaaa*t3*t46*t56*t90+1.26d
     7         -2*t3*t32*t90-2.52d-2*t3*t56*t86+4.2d-3*gammaaa*t117*t3*t
     8         32+2.52d-2*gammaaa*t102*t116*t3))*wght+Cmat3(iq,D3_GAA_GA
     9         A_GAA)
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
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = -param(1)
            t9 = param(2)
            t10 = param(3)
            t11 = 1/rhoa
            t12 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t13 = t12**5.0d-1
            t14 = t10**2
            t15 = rhoa**2
            t16 = 1/t15
            t17 = 1/t14
            t18 = 1/t12
            t19 = exp(-1.413716694115413d+1*t15*t17*t18)
            t20 = t19-1.0d+0
            t21 = -3.536776513153214d-2*t12*t14*t16*t20+t19-1.5d+0
            t22 = 1/t10
            t23 = 1/t13
            t24 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         2*t23)+2.659615202676212d-1*t10*t11*t13*t21
            t25 = -3.5461536035682834d-1*t10*t11*t13*t24*t9+t8+1.0d+0
            t26 = gammabb**5.0d-1
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/t27
            t29 = asinh(t26*t28)
            t30 = 2.52d-2*t26*t28*t29+1.0d+0
            t31 = 1/t30
            t32 = -4.2d-3*gammabb*t28*t31-9.305257363491002d-1*t27
            t33 = 1/rhob
            t34 = 4.2d-3*gammabb*t28*t31+9.305257363491002d-1*t27
            t35 = t34**5.0d-1
            t36 = rhob**2
            t37 = 1/t36
            t38 = 1/t34
            t39 = exp(-1.413716694115413d+1*t17*t36*t38)
            t40 = t39-1.0d+0
            t41 = -3.536776513153214d-2*t14*t34*t37*t40+t39-1.5d+0
            t42 = 1/t35
            t43 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         2*t42)+2.659615202676212d-1*t10*t33*t35*t41
            t44 = -3.5461536035682834d-1*t10*t33*t35*t43*t9+t8+1.0d+0
            t45 = 5.641895835477563d-1
            t46 = 1/rhoa**2.3333333333333334d+0
            t47 = 1/rhoa**2.6666666666666666d+0
            t48 = (gammaaa*t47+1)**5.0d-1
            t49 = 1/t48
            t50 = 1/rhoa**3.6666666666666664d+0
            t51 = -3.3599999999999997d-2*gammaaa*t49*t50-3.3599999999999
     1         997d-2*t1*t4*t46
            t52 = 1/t5**2
            t53 = rhoa**3.333333333333333d-1
            t54 = -5.599999999999999d-3*gammaaa*t46*t6+1.240700981798800
     1         2d+0*t53-4.2d-3*gammaaa*t3*t51*t52
            t55 = 1/t13**3
            t56 = 3.759942411946509d+0*t22*t23-1.8799712059732546d+0*rho
     1         a*t22*t54*t55
            t57 = rhoa**3
            t58 = 1/t57
            t59 = 1/t12**2
            t60 = 1.413716694115413d+1*t15*t17*t54*t59-2.827433388230826
     1         d+1*rhoa*t17*t18
            t61 = -3.536776513153214d-2*t12*t14*t16*t19*t60+t19*t60+7.07
     1         3553026306428d-2*t12*t14*t20*t58-3.536776513153214d-2*t14
     2         *t16*t20*t54
            t62 = 2.659615202676212d-1*t10*t11*t13*t61+3.544907701811040
     1         3d+0*t19*t45*t56+1.329807601338106d-1*t10*t11*t21*t23*t54
     2         -2.659615202676212d-1*t10*t13*t16*t21
            t63 = -3.5461536035682834d-1*t10*t11*t13*t62*t9-1.7730768017
     1         841417d-1*t10*t11*t23*t24*t54*t9+3.5461536035682834d-1*t1
     2         0*t13*t16*t24*t9
            t64 = 5.599999999999999d-3*gammaaa*t46*t6-1.2407009817988002
     1         d+0*t53+4.2d-3*gammaaa*t3*t51*t52
            t65 = 1/rhob**2.3333333333333334d+0
            t66 = 1/rhob**2.6666666666666666d+0
            t67 = (gammabb*t66+1)**5.0d-1
            t68 = 1/t67
            t69 = 1/rhob**3.6666666666666664d+0
            t70 = -3.3599999999999997d-2*gammabb*t68*t69-3.3599999999999
     1         997d-2*t26*t29*t65
            t71 = 1/t30**2
            t72 = rhob**3.333333333333333d-1
            t73 = 1.2407009817988002d+0*t72-4.2d-3*gammabb*t28*t70*t71-5
     1         .599999999999999d-3*gammabb*t31*t65
            t74 = 1/t35**3
            t75 = 3.759942411946509d+0*t22*t42-1.8799712059732546d+0*rho
     1         b*t22*t73*t74
            t76 = rhob**3
            t77 = 1/t76
            t78 = 1/t34**2
            t79 = 1.413716694115413d+1*t17*t36*t73*t78-2.827433388230826
     1         d+1*rhob*t17*t38
            t80 = -3.536776513153214d-2*t14*t34*t37*t39*t79+t39*t79+7.07
     1         3553026306428d-2*t14*t34*t40*t77-3.536776513153214d-2*t14
     2         *t37*t40*t73
            t81 = 2.659615202676212d-1*t10*t33*t35*t80+3.544907701811040
     1         3d+0*t39*t45*t75+1.329807601338106d-1*t10*t33*t41*t42*t73
     2         -2.659615202676212d-1*t10*t35*t37*t41
            t82 = -3.5461536035682834d-1*t10*t33*t35*t81*t9-1.7730768017
     1         841417d-1*t10*t33*t42*t43*t73*t9+3.5461536035682834d-1*t1
     2         0*t35*t37*t43*t9
            t83 = -1.2407009817988002d+0*t72+4.2d-3*gammabb*t28*t70*t71+
     1         5.599999999999999d-3*gammabb*t31*t65
            t84 = 1/t1
            t85 = 1.26d-2*t3*t4*t84+1.26d-2*t47*t49
            t86 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t52*t85
            t87 = 1.413716694115413d+1*t15*t17*t19*t59*t86-3.53677651315
     1         3214d-2*t14*t16*t20*t86-5.0d-1*t18*t19*t86
            t88 = 2.659615202676212d-1*t10*t11*t13*t87-6.66432440723758d
     1         +0*rhoa*t19*t22*t45*t55*t86+1.329807601338106d-1*t10*t11*
     2         t21*t23*t86
            t89 = -3.5461536035682834d-1*t10*t11*t13*t88*t9-1.7730768017
     1         841417d-1*t10*t11*t23*t24*t86*t9
            t90 = 4.2d-3*gammaaa*t3*t52*t85-4.2d-3*t3*t6
            t91 = 1/t26
            t92 = 1.26d-2*t28*t29*t91+1.26d-2*t66*t68
            t93 = 4.2d-3*t28*t31-4.2d-3*gammabb*t28*t71*t92
            t94 = 1.413716694115413d+1*t17*t36*t39*t78*t93-3.53677651315
     1         3214d-2*t14*t37*t40*t93-5.0d-1*t38*t39*t93
            t95 = 2.659615202676212d-1*t10*t33*t35*t94-6.66432440723758d
     1         +0*rhob*t22*t39*t45*t74*t93+1.329807601338106d-1*t10*t33*
     2         t41*t42*t93
            t96 = -3.5461536035682834d-1*t10*t33*t35*t9*t95-1.7730768017
     1         841417d-1*t10*t33*t42*t43*t9*t93
            t97 = 4.2d-3*gammabb*t28*t71*t92-4.2d-3*t28*t31
            t98 = t54**2
            t99 = 1/t13**5
            t100 = 1/rhoa**3.3333333333333337d+0
            t101 = t51**2
            t102 = 1/t5**3
            t103 = gammaaa**2
            t104 = 1/t48**3
            t105 = 1/rhoa**7.333333333333333d+0
            t106 = 1/rhoa**4.666666666666667d+0
            t107 = 1.6799999999999998d-1*gammaaa*t106*t49+7.84d-2*t1*t10
     1         0*t4-4.4799999999999995d-2*t103*t104*t105
            t108 = 1/rhoa**6.666666666666666d-1
            t109 = 1.3066666666666665d-2*gammaaa*t100*t6+1.1199999999999
     1         999d-2*gammaaa*t46*t51*t52-4.2d-3*gammaaa*t107*t3*t52+8.4
     2         d-3*gammaaa*t101*t102*t3+4.135669939329334d-1*t108
            t110 = 2.8199568089598814d+0*rhoa*t22*t98*t99-3.759942411946
     1         509d+0*t22*t54*t55-1.8799712059732546d+0*rhoa*t109*t22*t5
     2         5
            t111 = rhoa**4
            t112 = 1/t111
            t113 = 1/t12**3
            t114 = -2.827433388230826d+1*t113*t15*t17*t98+5.654866776461
     1         652d+1*rhoa*t17*t54*t59+1.413716694115413d+1*t109*t15*t17
     2         *t59-2.827433388230826d+1*t17*t18
            t115 = t60**2
            t116 = 1.4147106052612857d-1*t12*t14*t19*t58*t60-7.073553026
     1         306428d-2*t14*t16*t19*t54*t60+1.4147106052612857d-1*t14*t
     2         20*t54*t58-3.536776513153214d-2*t109*t14*t16*t20-2.122065
     3         9078919285d-1*t112*t12*t14*t20-3.536776513153214d-2*t115*
     4         t12*t14*t16*t19-3.536776513153214d-2*t114*t12*t14*t16*t19
     5         +t115*t19+t114*t19
            t117 = -6.649038006690531d-2*t10*t11*t21*t55*t98+2.659615202
     1         676212d-1*t10*t11*t23*t54*t61-5.319230405352424d-1*t10*t1
     2         3*t16*t61+3.5449077018110403d+0*t19*t45*t56*t60+5.3192304
     3         05352424d-1*t10*t13*t21*t58-2.659615202676212d-1*t10*t16*
     4         t21*t23*t54+3.5449077018110403d+0*t110*t19*t45+1.32980760
     5         1338106d-1*t10*t109*t11*t21*t23+2.659615202676212d-1*t10*
     6         t11*t116*t13
            t118 = 8.865384008920708d-2*t10*t11*t24*t55*t9*t98-3.5461536
     1         035682834d-1*t10*t11*t23*t54*t62*t9+7.092307207136567d-1*
     2         t10*t13*t16*t62*t9-7.092307207136567d-1*t10*t13*t24*t58*t
     3         9+3.5461536035682834d-1*t10*t16*t23*t24*t54*t9-1.77307680
     4         17841417d-1*t10*t109*t11*t23*t24*t9-3.5461536035682834d-1
     5         *t10*t11*t117*t13*t9
            t119 = -1.3066666666666665d-2*gammaaa*t100*t6-1.119999999999
     1         9999d-2*gammaaa*t46*t51*t52+4.2d-3*gammaaa*t107*t3*t52-8.
     2         4d-3*gammaaa*t101*t102*t3-4.135669939329334d-1*t108
            t120 = t73**2
            t121 = 1/t35**5
            t122 = 1/rhob**3.3333333333333337d+0
            t123 = t70**2
            t124 = 1/t30**3
            t125 = gammabb**2
            t126 = 1/t67**3
            t127 = 1/rhob**7.333333333333333d+0
            t128 = 1/rhob**4.666666666666667d+0
            t129 = 1.6799999999999998d-1*gammabb*t128*t68+7.84d-2*t122*t
     1         26*t29-4.4799999999999995d-2*t125*t126*t127
            t130 = 1/rhob**6.666666666666666d-1
            t131 = 1.1199999999999999d-2*gammabb*t65*t70*t71-4.2d-3*gamm
     1         abb*t129*t28*t71+1.3066666666666665d-2*gammabb*t122*t31+8
     2         .4d-3*gammabb*t123*t124*t28+4.135669939329334d-1*t130
            t132 = -3.759942411946509d+0*t22*t73*t74-1.8799712059732546d
     1         +0*rhob*t131*t22*t74+2.8199568089598814d+0*rhob*t120*t121
     2         *t22
            t133 = rhob**4
            t134 = 1/t133
            t135 = 1/t34**3
            t136 = 5.654866776461652d+1*rhob*t17*t73*t78+1.4137166941154
     1         13d+1*t131*t17*t36*t78-2.827433388230826d+1*t17*t38-2.827
     2         433388230826d+1*t120*t135*t17*t36
            t137 = t79**2
            t138 = 1.4147106052612857d-1*t14*t34*t39*t77*t79-7.073553026
     1         306428d-2*t14*t37*t39*t73*t79+1.4147106052612857d-1*t14*t
     2         40*t73*t77-3.536776513153214d-2*t131*t14*t37*t40-2.122065
     3         9078919285d-1*t134*t14*t34*t40-3.536776513153214d-2*t137*
     4         t14*t34*t37*t39-3.536776513153214d-2*t136*t14*t34*t37*t39
     5         +t137*t39+t136*t39
            t139 = 2.659615202676212d-1*t10*t33*t42*t73*t80-5.3192304053
     1         52424d-1*t10*t35*t37*t80+3.5449077018110403d+0*t39*t45*t7
     2         5*t79+5.319230405352424d-1*t10*t35*t41*t77-6.649038006690
     3         531d-2*t10*t120*t33*t41*t74-2.659615202676212d-1*t10*t37*
     4         t41*t42*t73+3.5449077018110403d+0*t132*t39*t45+1.32980760
     5         1338106d-1*t10*t131*t33*t41*t42+2.659615202676212d-1*t10*
     6         t138*t33*t35
            t140 = -3.5461536035682834d-1*t10*t33*t42*t73*t81*t9+7.09230
     1         7207136567d-1*t10*t35*t37*t81*t9-7.092307207136567d-1*t10
     2         *t35*t43*t77*t9+8.865384008920708d-2*t10*t120*t33*t43*t74
     3         *t9+3.5461536035682834d-1*t10*t37*t42*t43*t73*t9-1.773076
     4         8017841417d-1*t10*t131*t33*t42*t43*t9-3.5461536035682834d
     5         -1*t10*t139*t33*t35*t9
            t141 = -1.1199999999999999d-2*gammabb*t65*t70*t71+4.2d-3*gam
     1         mabb*t129*t28*t71-1.3066666666666665d-2*gammabb*t122*t31-
     2         8.4d-3*gammabb*t123*t124*t28-4.135669939329334d-1*t130
            t142 = 1/rhoa**6.333333333333333d+0
            t143 = -1.6799999999999998d-2*t4*t46*t84-5.04d-2*t49*t50+1.6
     1         799999999999998d-2*gammaaa*t104*t142
            t144 = 5.599999999999999d-3*gammaaa*t46*t52*t85+8.4d-3*gamma
     1         aa*t102*t3*t51*t85-5.599999999999999d-3*t46*t6-4.2d-3*t3*
     2         t51*t52-4.2d-3*gammaaa*t143*t3*t52
            t145 = 2.8199568089598814d+0*rhoa*t22*t54*t86*t99-1.87997120
     1         59732546d+0*t22*t55*t86-1.8799712059732546d+0*rhoa*t144*t
     2         22*t55
            t146 = 2.827433388230826d+1*rhoa*t17*t59*t86-2.8274333882308
     1         26d+1*t113*t15*t17*t54*t86+1.413716694115413d+1*t144*t15*
     2         t17*t59
            t147 = 1.413716694115413d+1*t15*t17*t19*t59*t60*t86-5.0d-1*t
     1         18*t19*t60*t86-3.536776513153214d-2*t14*t16*t19*t60*t86-5
     2         .0d-1*t19*t54*t59*t86+7.073553026306428d-2*t14*t20*t58*t8
     3         6+1.0d+0*t11*t18*t19*t86-3.536776513153214d-2*t14*t144*t1
     4         6*t20-3.536776513153214d-2*t12*t14*t146*t16*t19+t146*t19
            t148 = 1.329807601338106d-1*t10*t11*t23*t54*t87-2.6596152026
     1         76212d-1*t10*t13*t16*t87+1.329807601338106d-1*t10*t11*t23
     2         *t61*t86+5.01149519714857d+1*t15*t17*t19*t45*t56*t59*t86-
     3         6.649038006690531d-2*t10*t11*t21*t54*t55*t86-1.3298076013
     4         38106d-1*t10*t16*t21*t23*t86+3.5449077018110403d+0*t145*t
     5         19*t45+1.329807601338106d-1*t10*t11*t144*t21*t23+2.659615
     6         202676212d-1*t10*t11*t13*t147
            t149 = -1.7730768017841417d-1*t10*t11*t23*t54*t88*t9+3.54615
     1         36035682834d-1*t10*t13*t16*t88*t9-1.7730768017841417d-1*t
     2         10*t11*t23*t62*t86*t9+8.865384008920708d-2*t10*t11*t24*t5
     3         4*t55*t86*t9+1.7730768017841417d-1*t10*t16*t23*t24*t86*t9
     4         -1.7730768017841417d-1*t10*t11*t144*t23*t24*t9-3.54615360
     5         35682834d-1*t10*t11*t13*t148*t9
            t150 = -5.599999999999999d-3*gammaaa*t46*t52*t85-8.4d-3*gamm
     1         aaa*t102*t3*t51*t85+5.599999999999999d-3*t46*t6+4.2d-3*t3
     2         *t51*t52+4.2d-3*gammaaa*t143*t3*t52
            t151 = 1/rhob**6.333333333333333d+0
            t152 = -1.6799999999999998d-2*t29*t65*t91-5.04d-2*t68*t69+1.
     1         6799999999999998d-2*gammabb*t126*t151
            t153 = 5.599999999999999d-3*gammabb*t65*t71*t92+8.4d-3*gamma
     1         bb*t124*t28*t70*t92-4.2d-3*t28*t70*t71-4.2d-3*gammabb*t15
     2         2*t28*t71-5.599999999999999d-3*t31*t65
            t154 = -1.8799712059732546d+0*t22*t74*t93+2.8199568089598814
     1         d+0*rhob*t121*t22*t73*t93-1.8799712059732546d+0*rhob*t153
     2         *t22*t74
            t155 = 2.827433388230826d+1*rhob*t17*t78*t93-2.8274333882308
     1         26d+1*t135*t17*t36*t73*t93+1.413716694115413d+1*t153*t17*
     2         t36*t78
            t156 = 1.413716694115413d+1*t17*t36*t39*t78*t79*t93-5.0d-1*t
     1         38*t39*t79*t93-3.536776513153214d-2*t14*t37*t39*t79*t93-5
     2         .0d-1*t39*t73*t78*t93+7.073553026306428d-2*t14*t40*t77*t9
     3         3+1.0d+0*t33*t38*t39*t93-3.536776513153214d-2*t14*t153*t3
     4         7*t40-3.536776513153214d-2*t14*t155*t34*t37*t39+t155*t39
            t157 = 1.329807601338106d-1*t10*t33*t42*t73*t94-2.6596152026
     1         76212d-1*t10*t35*t37*t94+1.329807601338106d-1*t10*t33*t42
     2         *t80*t93+5.01149519714857d+1*t17*t36*t39*t45*t75*t78*t93-
     3         6.649038006690531d-2*t10*t33*t41*t73*t74*t93-1.3298076013
     4         38106d-1*t10*t37*t41*t42*t93+3.5449077018110403d+0*t154*t
     5         39*t45+1.329807601338106d-1*t10*t153*t33*t41*t42+2.659615
     6         202676212d-1*t10*t156*t33*t35
            t158 = -1.7730768017841417d-1*t10*t33*t42*t73*t9*t95+3.54615
     1         36035682834d-1*t10*t35*t37*t9*t95-1.7730768017841417d-1*t
     2         10*t33*t42*t81*t9*t93+8.865384008920708d-2*t10*t33*t43*t7
     3         3*t74*t9*t93+1.7730768017841417d-1*t10*t37*t42*t43*t9*t93
     4         -1.7730768017841417d-1*t10*t153*t33*t42*t43*t9-3.54615360
     5         35682834d-1*t10*t157*t33*t35*t9
            t159 = -5.599999999999999d-3*gammabb*t65*t71*t92-8.4d-3*gamm
     1         abb*t124*t28*t70*t92+4.2d-3*t28*t70*t71+4.2d-3*gammabb*t1
     2         52*t28*t71+5.599999999999999d-3*t31*t65
            t160 = 1/t10**3
            t161 = t86**2
            t162 = 1/t13**7
            t163 = t85**2
            t164 = 1/rhoa**5.333333333333333d+0
            t165 = 1/gammaaa
            t166 = 1/t1**3
            t167 = 6.3d-3*t165*t47*t49-6.3d-3*t166*t3*t4-6.3d-3*t104*t16
     1         4
            t168 = -8.4d-3*t3*t52*t85-4.2d-3*gammaaa*t167*t3*t52+8.4d-3*
     1         gammaaa*t102*t163*t3
            t169 = 1/t10**4
            t170 = 1/t12**4
            t171 = 1.413716694115413d+1*t15*t168*t17*t19*t59-3.536776513
     1         153214d-2*t14*t16*t168*t20-5.0d-1*t168*t18*t19+1.99859489
     2         12206125d+2*t111*t161*t169*t170*t19-3.534291735288533d+1*
     3         t113*t15*t161*t17*t19
            t172 = 9.996486610856369d+0*rhoa*t161*t19*t22*t45*t99+2.6596
     1         15202676212d-1*t10*t11*t23*t86*t87-9.421466669512572d+1*t
     2         160*t161*t162*t19*t45*t57-6.66432440723758d+0*rhoa*t168*t
     3         19*t22*t45*t55-6.649038006690531d-2*t10*t11*t161*t21*t55+
     4         1.329807601338106d-1*t10*t11*t168*t21*t23+2.6596152026762
     5         12d-1*t10*t11*t13*t171
            t173 = -3.5461536035682834d-1*t10*t11*t23*t86*t88*t9+8.86538
     1         4008920708d-2*t10*t11*t161*t24*t55*t9-1.7730768017841417d
     2         -1*t10*t11*t168*t23*t24*t9-3.5461536035682834d-1*t10*t11*
     3         t13*t172*t9
            t174 = 8.4d-3*t3*t52*t85+4.2d-3*gammaaa*t167*t3*t52-8.4d-3*g
     1         ammaaa*t102*t163*t3
            t175 = t93**2
            t176 = 1/t35**7
            t177 = t92**2
            t178 = 1/rhob**5.333333333333333d+0
            t179 = 1/gammabb
            t180 = 1/t26**3
            t181 = 6.3d-3*t179*t66*t68-6.3d-3*t180*t28*t29-6.3d-3*t126*t
     1         178
            t182 = -8.4d-3*t28*t71*t92-4.2d-3*gammabb*t181*t28*t71+8.4d-
     1         3*gammabb*t124*t177*t28
            t183 = 1/t34**4
            t184 = 1.413716694115413d+1*t17*t182*t36*t39*t78-3.536776513
     1         153214d-2*t14*t182*t37*t40-5.0d-1*t182*t38*t39-3.53429173
     2         5288533d+1*t135*t17*t175*t36*t39+1.9985948912206125d+2*t1
     3         33*t169*t175*t183*t39
            t185 = 2.659615202676212d-1*t10*t33*t42*t93*t94-9.4214666695
     1         12572d+1*t160*t175*t176*t39*t45*t76-6.66432440723758d+0*r
     2         hob*t182*t22*t39*t45*t74-6.649038006690531d-2*t10*t175*t3
     3         3*t41*t74+9.996486610856369d+0*rhob*t121*t175*t22*t39*t45
     4         +1.329807601338106d-1*t10*t182*t33*t41*t42+2.659615202676
     5         212d-1*t10*t184*t33*t35
            t186 = -3.5461536035682834d-1*t10*t33*t42*t9*t93*t95+8.86538
     1         4008920708d-2*t10*t175*t33*t43*t74*t9-1.7730768017841417d
     2         -1*t10*t182*t33*t42*t43*t9-3.5461536035682834d-1*t10*t185
     3         *t33*t35*t9
            t187 = 8.4d-3*t28*t71*t92+4.2d-3*gammabb*t181*t28*t71-8.4d-3
     1         *gammabb*t124*t177*t28
            t188 = t54**3
            t189 = 1/rhoa**4.333333333333333d+0
            t190 = 1/rhoa**1.6666666666666669d+0
            t191 = t51**3
            t192 = 1/t5**4
            t193 = 1/t48**5
            t194 = -8.885333333333333d-1*gammaaa*t49/rhoa**5.66666666666
     1         6667d+0-2.6133333333333336d-1*t1*t189*t4-1.79199999999999
     2         97d-1*gammaaa**3*t193/rhoa**11+5.525333333333333d-1*t103*
     3         t104/rhoa**8.333333333333334d+0
            t195 = -4.355555555555556d-2*gammaaa*t189*t6-3.92d-2*gammaaa
     1         *t100*t51*t52+1.6799999999999995d-2*gammaaa*t107*t46*t52-
     2         4.2d-3*gammaaa*t194*t3*t52+2.52d-2*gammaaa*t102*t107*t3*t
     3         51-3.359999999999999d-2*gammaaa*t101*t102*t46-2.52d-2*gam
     4         maaa*t191*t192*t3-2.7571132928862224d-1*t190
            t196 = rhoa**5
            t197 = -1.6964600329384957d+2*rhoa*t113*t17*t98+8.4823001646
     1         92478d+1*t17*t54*t59+1.413716694115413d+1*t15*t17*t195*t5
     2         9+8.482300164692478d+1*rhoa*t109*t17*t59-8.48230016469247
     3         8d+1*t109*t113*t15*t17*t54+8.482300164692478d+1*t15*t17*t
     4         170*t188
            t198 = t60**3
            t199 = t73**3
            t200 = 1/rhob**4.333333333333333d+0
            t201 = 1/rhob**1.6666666666666669d+0
            t202 = t70**3
            t203 = 1/t30**4
            t204 = 1/t67**5
            t205 = -8.885333333333333d-1*gammabb*t68/rhob**5.66666666666
     1         6667d+0-2.6133333333333336d-1*t200*t26*t29-1.791999999999
     2         9997d-1*gammabb**3*t204/rhob**11+5.525333333333333d-1*t12
     3         5*t126/rhob**8.333333333333334d+0
            t206 = -3.92d-2*gammabb*t122*t70*t71+1.6799999999999995d-2*g
     1         ammabb*t129*t65*t71-4.2d-3*gammabb*t205*t28*t71+2.52d-2*g
     2         ammabb*t124*t129*t28*t70-3.359999999999999d-2*gammabb*t12
     3         3*t124*t65-4.355555555555556d-2*gammabb*t200*t31-2.52d-2*
     4         gammabb*t202*t203*t28-2.7571132928862224d-1*t201
            t207 = rhob**5
            t208 = 8.482300164692478d+1*t17*t73*t78+1.413716694115413d+1
     1         *t17*t206*t36*t78+8.482300164692478d+1*rhob*t131*t17*t78-
     2         8.482300164692478d+1*t131*t135*t17*t36*t73+8.482300164692
     3         478d+1*t17*t183*t199*t36-1.6964600329384957d+2*rhob*t120*
     4         t135*t17
            t209 = t79**3
            t210 = 3.92d-2*t100*t4*t84+2.072d-1*t106*t49+6.7199999999999
     1         98d-2*t103*t193/rhoa**10-1.7359999999999998d-1*gammaaa*t1
     2         04*t105
            t211 = -1.3066666666666665d-2*gammaaa*t100*t52*t85-2.2399999
     1         999999998d-2*gammaaa*t102*t46*t51*t85-2.52d-2*gammaaa*t10
     2         1*t192*t3*t85+8.4d-3*gammaaa*t102*t107*t3*t85+1.306666666
     3         6666665d-2*t100*t6+1.1199999999999999d-2*t46*t51*t52+1.11
     4         99999999999999d-2*gammaaa*t143*t46*t52-4.2d-3*gammaaa*t21
     5         0*t3*t52-4.2d-3*t107*t3*t52+1.6799999999999998d-2*gammaaa
     6         *t102*t143*t3*t51+8.4d-3*t101*t102*t3
            t212 = 8.482300164692478d+1*t15*t17*t170*t86*t98+2.827433388
     1         230826d+1*t17*t59*t86-1.1309733552923305d+2*rhoa*t113*t17
     2         *t54*t86-2.827433388230826d+1*t109*t113*t15*t17*t86+1.413
     3         716694115413d+1*t15*t17*t211*t59+5.654866776461652d+1*rho
     4         a*t144*t17*t59-5.654866776461652d+1*t113*t144*t15*t17*t54
            t213 = 3.92d-2*t122*t29*t91+2.072d-1*t128*t68+6.719999999999
     1         998d-2*t125*t204/rhob**10-1.7359999999999998d-1*gammabb*t
     2         126*t127
            t214 = -1.3066666666666665d-2*gammabb*t122*t71*t92-2.2399999
     1         999999998d-2*gammabb*t124*t65*t70*t92-2.52d-2*gammabb*t12
     2         3*t203*t28*t92+8.4d-3*gammabb*t124*t129*t28*t92+1.1199999
     3         999999999d-2*t65*t70*t71+1.1199999999999999d-2*gammabb*t1
     4         52*t65*t71-4.2d-3*gammabb*t213*t28*t71-4.2d-3*t129*t28*t7
     5         1+1.6799999999999998d-2*gammabb*t124*t152*t28*t70+1.30666
     6         66666666665d-2*t122*t31+8.4d-3*t123*t124*t28
            t215 = 2.827433388230826d+1*t17*t78*t93-1.1309733552923305d+
     1         2*rhob*t135*t17*t73*t93+8.482300164692478d+1*t120*t17*t18
     2         3*t36*t93-2.827433388230826d+1*t131*t135*t17*t36*t93+1.41
     3         3716694115413d+1*t17*t214*t36*t78+5.654866776461652d+1*rh
     4         ob*t153*t17*t78-5.654866776461652d+1*t135*t153*t17*t36*t7
     5         3
            t216 = -8.4d-3*t165*t49*t50+8.4d-3*t166*t4*t46-2.52d-2*gamma
     1         aa*t193/rhoa**9+4.1999999999999993d-2*t104*t142
            t217 = 1.1199999999999999d-2*t46*t52*t85+1.6799999999999998d
     1         -2*t102*t3*t51*t85+1.6799999999999998d-2*gammaaa*t102*t14
     2         3*t3*t85+5.599999999999999d-3*gammaaa*t167*t46*t52-4.2d-3
     3         *gammaaa*t216*t3*t52-8.4d-3*t143*t3*t52-2.52d-2*gammaaa*t
     4         163*t192*t3*t51+8.4d-3*gammaaa*t102*t167*t3*t51-1.1199999
     5         999999999d-2*gammaaa*t102*t163*t46
            t218 = -5.654866776461652d+1*t113*t144*t15*t17*t86+1.4137166
     1         94115413d+1*t15*t17*t217*t59+2.827433388230826d+1*rhoa*t1
     2         68*t17*t59+8.482300164692478d+1*t15*t161*t17*t170*t54-2.8
     3         27433388230826d+1*t113*t15*t168*t17*t54-5.654866776461652
     4         d+1*rhoa*t113*t161*t17
            t219 = -8.4d-3*t179*t68*t69+8.4d-3*t180*t29*t65-2.52d-2*gamm
     1         abb*t204/rhob**9+4.1999999999999993d-2*t126*t151
            t220 = 1.1199999999999999d-2*t65*t71*t92+1.6799999999999998d
     1         -2*t124*t28*t70*t92+1.6799999999999998d-2*gammabb*t124*t1
     2         52*t28*t92+5.599999999999999d-3*gammabb*t181*t65*t71-4.2d
     3         -3*gammabb*t219*t28*t71-8.4d-3*t152*t28*t71-2.52d-2*gamma
     4         bb*t177*t203*t28*t70+8.4d-3*gammabb*t124*t181*t28*t70-1.1
     5         199999999999999d-2*gammabb*t124*t177*t65
            t221 = -5.654866776461652d+1*t135*t153*t17*t36*t93+1.4137166
     1         94115413d+1*t17*t220*t36*t78+2.827433388230826d+1*rhob*t1
     2         7*t182*t78+8.482300164692478d+1*t17*t175*t183*t36*t73-2.8
     3         27433388230826d+1*t135*t17*t182*t36*t73-5.654866776461652
     4         d+1*rhob*t135*t17*t175
            t222 = 1/t10**5
            t223 = t86**3
            t224 = t85**3
            t225 = -9.450000000000001d-3*t47*t49/t103+9.450000000000001d
     1         -3*t3*t4/t1**5+9.450000000000001d-3*t193/rhoa**8-3.15d-3*
     2         t104*t164*t165
            t226 = 2.52d-2*gammaaa*t102*t167*t3*t85-4.2d-3*gammaaa*t225*
     1         t3*t52-1.26d-2*t167*t3*t52-2.52d-2*gammaaa*t192*t224*t3+2
     2         .52d-2*t102*t163*t3
            t227 = 1/t10**6
            t228 = t93**3
            t229 = t92**3
            t230 = -9.450000000000001d-3*t66*t68/t125+9.450000000000001d
     1         -3*t28*t29/t26**5+9.450000000000001d-3*t204/rhob**8-3.15d
     2         -3*t126*t178*t179
            t231 = 2.52d-2*gammabb*t124*t181*t28*t92-4.2d-3*gammabb*t230
     1         *t28*t71-1.26d-2*t181*t28*t71-2.52d-2*gammabb*t203*t229*t
     2         28+2.52d-2*t124*t177*t28
            fnc(iq) = (t25*t7+t32*t44)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t63*t7+t25*t64)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t44*t83+t32*t82)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t25*t90+t7*t89)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t44*t97+t32*t96)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t118*t7+2*t63*t64+t119*t25)*wght+Amat2
     1         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (2*t82*t83+t141*t44+t140*t32)*wght+Amat
     1         2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t63*t90+t64*t89+t149*t7+t150*t25)*wgh
     1         t+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t82*t97+t83*t96+t159*t44+t158*t32)*wg
     1         ht+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2*t89*t90+t173*t7+t174*t25)*wght+Cma
     1         t2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (2*t96*t97+t187*t44+t186*t32)*wght+Cm
     1         at2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (t7*(-3.5461536035682834d-1*t10*t11*
     1         t13*t9*(3.5449077018110403d+0*t19*t45*(8.459870426879645d
     2         +0*t22*t98*t99+8.459870426879645d+0*rhoa*t109*t22*t54*t99
     3         -1.8799712059732546d+0*rhoa*t195*t22*t55-5.63991361791976
     4         3d+0*t109*t22*t55-7.049892022399704d+0*rhoa*t162*t188*t22
     5         )+9.973557010035795d-2*t10*t11*t188*t21*t99-1.99471140200
     6         7159d-1*t10*t11*t55*t61*t98+1.994711402007159d-1*t10*t16*
     7         t21*t55*t98+1.5957691216057274d+0*t10*t13*t58*t61-7.97884
     8         5608028636d-1*t10*t16*t23*t54*t61+3.989422804014318d-1*t1
     9         0*t109*t11*t23*t61+2.659615202676212d-1*t10*t11*t13*(4.24
     :         4131815783857d-1*t14*t19*t54*t58*t60-1.0610329539459643d-
     ;         1*t114*t12*t14*t16*t19*t60-1.0610329539459643d-1*t109*t14
     <         *t16*t19*t60-6.366197723675786d-1*t112*t12*t14*t19*t60+3*
     =         t114*t19*t60+2.1220659078919285d-1*t109*t14*t20*t58+2.122
     >         0659078919285d-1*t115*t12*t14*t19*t58+2.1220659078919285d
     ?         -1*t114*t12*t14*t19*t58-6.366197723675786d-1*t112*t14*t20
     @         *t54-1.0610329539459643d-1*t115*t14*t16*t19*t54-1.0610329
     1         539459643d-1*t114*t14*t16*t19*t54+8.488263631567714d-1*t1
     2         2*t14*t20/t196-3.536776513153214d-2*t14*t16*t195*t20-3.53
     3         6776513153214d-2*t12*t14*t16*t19*t198+t19*t198-3.53677651
     4         3153214d-2*t12*t14*t16*t19*t197+t19*t197)+7.0898154036220
     5         81d+0*t110*t19*t45*t60+7.978845608028636d-1*t10*t21*t23*t
     6         54*t58+3.5449077018110403d+0*t115*t19*t45*t56+3.544907701
     7         8110403d+0*t114*t19*t45*t56-1.994711402007159d-1*t10*t109
     8         *t11*t21*t54*t55+3.989422804014318d-1*t10*t11*t116*t23*t5
     9         4+1.329807601338106d-1*t10*t11*t195*t21*t23-3.98942280401
     :         4318d-1*t10*t109*t16*t21*t23-1.5957691216057274d+0*t10*t1
     ;         12*t13*t21-7.978845608028636d-1*t10*t116*t13*t16)-1.32980
     <         76013381063d-1*t10*t11*t188*t24*t9*t99+2.6596152026762127
     =         d-1*t10*t11*t55*t62*t9*t98-2.6596152026762127d-1*t10*t16*
     >         t24*t55*t9*t98-2.12769216214097d+0*t10*t13*t58*t62*t9+1.0
     ?         63846081070485d+0*t10*t16*t23*t54*t62*t9-5.31923040535242
     @         5d-1*t10*t109*t11*t23*t62*t9-1.063846081070485d+0*t10*t23
     1         *t24*t54*t58*t9+2.6596152026762127d-1*t10*t109*t11*t24*t5
     2         4*t55*t9-5.319230405352425d-1*t10*t11*t117*t23*t54*t9-1.7
     3         730768017841417d-1*t10*t11*t195*t23*t24*t9+5.319230405352
     4         425d-1*t10*t109*t16*t23*t24*t9+2.12769216214097d+0*t10*t1
     5         12*t13*t24*t9+1.063846081070485d+0*t10*t117*t13*t16*t9)+3
     6         *t118*t64+3*t119*t63+t25*(4.355555555555556d-2*gammaaa*t1
     7         89*t6+3.92d-2*gammaaa*t100*t51*t52-1.6799999999999995d-2*
     8         gammaaa*t107*t46*t52+4.2d-3*gammaaa*t194*t3*t52-2.52d-2*g
     9         ammaaa*t102*t107*t3*t51+3.359999999999999d-2*gammaaa*t101
     :         *t102*t46+2.52d-2*gammaaa*t191*t192*t3+2.7571132928862224
     ;         d-1*t190))*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t32*(-2.12769216214097d+0*t10*t35*t
     1         77*t81*t9+2.6596152026762127d-1*t10*t120*t33*t74*t81*t9+1
     2         .063846081070485d+0*t10*t37*t42*t73*t81*t9-5.319230405352
     3         425d-1*t10*t131*t33*t42*t81*t9-3.5461536035682834d-1*t10*
     4         t33*t35*(1.5957691216057274d+0*t10*t35*t77*t80-1.99471140
     5         2007159d-1*t10*t120*t33*t74*t80-7.978845608028636d-1*t10*
     6         t37*t42*t73*t80+3.989422804014318d-1*t10*t131*t33*t42*t80
     7         +2.659615202676212d-1*t10*t33*t35*(4.244131815783857d-1*t
     8         14*t39*t73*t77*t79-1.0610329539459643d-1*t136*t14*t34*t37
     9         *t39*t79-1.0610329539459643d-1*t131*t14*t37*t39*t79-6.366
     :         197723675786d-1*t134*t14*t34*t39*t79+3*t136*t39*t79+2.122
     ;         0659078919285d-1*t131*t14*t40*t77+2.1220659078919285d-1*t
     <         137*t14*t34*t39*t77+2.1220659078919285d-1*t136*t14*t34*t3
     =         9*t77-6.366197723675786d-1*t134*t14*t40*t73-1.06103295394
     >         59643d-1*t137*t14*t37*t39*t73-1.0610329539459643d-1*t136*
     ?         t14*t37*t39*t73-3.536776513153214d-2*t14*t206*t37*t40+8.4
     @         88263631567714d-1*t14*t34*t40/t207-3.536776513153214d-2*t
     1         14*t209*t34*t37*t39-3.536776513153214d-2*t14*t208*t34*t37
     2         *t39+t209*t39+t208*t39)+7.089815403622081d+0*t132*t39*t45
     3         *t79+7.978845608028636d-1*t10*t41*t42*t73*t77+3.544907701
     4         8110403d+0*t137*t39*t45*t75+3.5449077018110403d+0*t136*t3
     5         9*t45*t75+3.5449077018110403d+0*t39*t45*(-1.8799712059732
     6         546d+0*rhob*t206*t22*t74-5.639913617919763d+0*t131*t22*t7
     7         4+8.459870426879645d+0*rhob*t121*t131*t22*t73-7.049892022
     8         399704d+0*rhob*t176*t199*t22+8.459870426879645d+0*t120*t1
     9         21*t22)-1.994711402007159d-1*t10*t131*t33*t41*t73*t74+1.9
     :         94711402007159d-1*t10*t120*t37*t41*t74+3.989422804014318d
     ;         -1*t10*t138*t33*t42*t73-3.989422804014318d-1*t10*t131*t37
     <         *t41*t42+1.329807601338106d-1*t10*t206*t33*t41*t42-1.5957
     =         691216057274d+0*t10*t134*t35*t41+9.973557010035795d-2*t10
     >         *t121*t199*t33*t41-7.978845608028636d-1*t10*t138*t35*t37)
     ?         *t9-1.063846081070485d+0*t10*t42*t43*t73*t77*t9+2.6596152
     @         026762127d-1*t10*t131*t33*t43*t73*t74*t9-2.65961520267621
     1         27d-1*t10*t120*t37*t43*t74*t9-5.319230405352425d-1*t10*t1
     2         39*t33*t42*t73*t9+5.319230405352425d-1*t10*t131*t37*t42*t
     3         43*t9-1.7730768017841417d-1*t10*t206*t33*t42*t43*t9+2.127
     4         69216214097d+0*t10*t134*t35*t43*t9-1.3298076013381063d-1*
     5         t10*t121*t199*t33*t43*t9+1.063846081070485d+0*t10*t139*t3
     6         5*t37*t9)+3*t140*t83+3*t141*t82+t44*(3.92d-2*gammabb*t122
     7         *t70*t71-1.6799999999999995d-2*gammabb*t129*t65*t71+4.2d-
     8         3*gammabb*t205*t28*t71-2.52d-2*gammabb*t124*t129*t28*t70+
     9         3.359999999999999d-2*gammabb*t123*t124*t65+4.355555555555
     :         556d-2*gammabb*t200*t31+2.52d-2*gammabb*t202*t203*t28+2.7
     ;         571132928862224d-1*t201))*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t7*(-3.5461536035682834d-1*t10*t11
     1         *t13*t9*(3.5449077018110403d+0*t19*t45*(5.639913617919763
     2         d+0*t22*t54*t86*t99+2.8199568089598814d+0*rhoa*t109*t22*t
     3         86*t99+5.639913617919763d+0*rhoa*t144*t22*t54*t99-7.04989
     4         2022399704d+0*rhoa*t162*t22*t86*t98-1.8799712059732546d+0
     5         *rhoa*t211*t22*t55-3.759942411946509d+0*t144*t22*t55)+9.9
     6         73557010035795d-2*t10*t11*t21*t86*t98*t99-6.6490380066905
     7         31d-2*t10*t11*t55*t87*t98+5.319230405352424d-1*t10*t13*t5
     8         8*t87-2.659615202676212d-1*t10*t16*t23*t54*t87+1.32980760
     9         1338106d-1*t10*t109*t11*t23*t87+2.659615202676212d-1*t10*
     :         t11*t13*(-1.0d+0*t19*t54*t59*t60*t86+1.4147106052612857d-
     ;         1*t14*t19*t58*t60*t86+2.0d+0*t11*t18*t19*t60*t86+2.0d+0*t
     <         11*t19*t54*t59*t86+1.413716694115413d+1*t115*t15*t17*t19*
     =         t59*t86+1.413716694115413d+1*t114*t15*t17*t19*t59*t86-5.0
     >         d-1*t109*t19*t59*t86-2.1220659078919285d-1*t112*t14*t20*t
     ?         86-3.0d+0*t16*t18*t19*t86-5.0d-1*t115*t18*t19*t86-5.0d-1*
     @         t114*t18*t19*t86-3.536776513153214d-2*t115*t14*t16*t19*t8
     1         6-3.536776513153214d-2*t114*t14*t16*t19*t86-7.07355302630
     2         6428d-2*t12*t14*t146*t16*t19*t60-7.073553026306428d-2*t14
     3         *t144*t16*t19*t60+2*t146*t19*t60+1.4147106052612857d-1*t1
     4         4*t144*t20*t58+1.4147106052612857d-1*t12*t14*t146*t19*t58
     5         -7.073553026306428d-2*t14*t146*t16*t19*t54-3.536776513153
     6         214d-2*t12*t14*t16*t19*t212+t19*t212-3.536776513153214d-2
     7         *t14*t16*t20*t211)-1.329807601338106d-1*t10*t11*t54*t55*t
     8         61*t86-2.659615202676212d-1*t10*t16*t23*t61*t86+5.0114951
     9         9714857d+1*t15*t17*t19*t45*t56*t59*t60*t86+5.011495197148
     :         57d+1*t110*t15*t17*t19*t45*t59*t86+2.659615202676212d-1*t
     ;         10*t21*t23*t58*t86+1.329807601338106d-1*t10*t16*t21*t54*t
     <         55*t86-6.649038006690531d-2*t10*t109*t11*t21*t55*t86+1.32
     =         9807601338106d-1*t10*t11*t116*t23*t86+2.659615202676212d-
     >         1*t10*t11*t144*t23*t61+3.5449077018110403d+0*t145*t19*t45
     ?         *t60+3.5449077018110403d+0*t146*t19*t45*t56-1.32980760133
     @         8106d-1*t10*t11*t144*t21*t54*t55+2.659615202676212d-1*t10
     1         *t11*t147*t23*t54+1.329807601338106d-1*t10*t11*t21*t211*t
     2         23-2.659615202676212d-1*t10*t144*t16*t21*t23-5.3192304053
     3         52424d-1*t10*t13*t147*t16)-1.3298076013381063d-1*t10*t11*
     4         t24*t86*t9*t98*t99+8.865384008920708d-2*t10*t11*t55*t88*t
     5         9*t98-7.092307207136567d-1*t10*t13*t58*t88*t9+3.546153603
     6         5682834d-1*t10*t16*t23*t54*t88*t9-1.7730768017841417d-1*t
     7         10*t109*t11*t23*t88*t9+1.7730768017841417d-1*t10*t11*t54*
     8         t55*t62*t86*t9+3.5461536035682834d-1*t10*t16*t23*t62*t86*
     9         t9-3.5461536035682834d-1*t10*t23*t24*t58*t86*t9-1.7730768
     :         017841417d-1*t10*t16*t24*t54*t55*t86*t9+8.865384008920708
     ;         d-2*t10*t109*t11*t24*t55*t86*t9-1.7730768017841417d-1*t10
     <         *t11*t117*t23*t86*t9-3.5461536035682834d-1*t10*t11*t144*t
     =         23*t62*t9+1.7730768017841417d-1*t10*t11*t144*t24*t54*t55*
     >         t9-3.5461536035682834d-1*t10*t11*t148*t23*t54*t9-1.773076
     ?         8017841417d-1*t10*t11*t211*t23*t24*t9+3.5461536035682834d
     @         -1*t10*t144*t16*t23*t24*t9+7.092307207136567d-1*t10*t13*t
     1         148*t16*t9)+t118*t90+t119*t89+t25*(1.3066666666666665d-2*
     2         gammaaa*t100*t52*t85+2.2399999999999998d-2*gammaaa*t102*t
     3         46*t51*t85+2.52d-2*gammaaa*t101*t192*t3*t85-8.4d-3*gammaa
     4         a*t102*t107*t3*t85-1.3066666666666665d-2*t100*t6-1.119999
     5         9999999999d-2*t46*t51*t52-1.1199999999999999d-2*gammaaa*t
     6         143*t46*t52+4.2d-3*gammaaa*t210*t3*t52+4.2d-3*t107*t3*t52
     7         -1.6799999999999998d-2*gammaaa*t102*t143*t3*t51-8.4d-3*t1
     8         01*t102*t3)+2*t149*t64+2*t150*t63)*wght+Cmat3(iq,D3_RA_RA
     9         _GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t140*t97+t141*t96+t32*(-7.09230720
     1         7136567d-1*t10*t35*t77*t9*t95+8.865384008920708d-2*t10*t1
     2         20*t33*t74*t9*t95+3.5461536035682834d-1*t10*t37*t42*t73*t
     3         9*t95-1.7730768017841417d-1*t10*t131*t33*t42*t9*t95-3.546
     4         1536035682834d-1*t10*t33*t35*t9*(5.319230405352424d-1*t10
     5         *t35*t77*t94-6.649038006690531d-2*t10*t120*t33*t74*t94-2.
     6         659615202676212d-1*t10*t37*t42*t73*t94+1.329807601338106d
     7         -1*t10*t131*t33*t42*t94+2.659615202676212d-1*t10*t33*t35*
     8         (-1.0d+0*t39*t73*t78*t79*t93+1.4147106052612857d-1*t14*t3
     9         9*t77*t79*t93+2.0d+0*t33*t38*t39*t79*t93+2.0d+0*t33*t39*t
     :         73*t78*t93+1.413716694115413d+1*t137*t17*t36*t39*t78*t93+
     ;         1.413716694115413d+1*t136*t17*t36*t39*t78*t93-5.0d-1*t131
     <         *t39*t78*t93-2.1220659078919285d-1*t134*t14*t40*t93-3.0d+
     =         0*t37*t38*t39*t93-5.0d-1*t137*t38*t39*t93-5.0d-1*t136*t38
     >         *t39*t93-3.536776513153214d-2*t137*t14*t37*t39*t93-3.5367
     ?         76513153214d-2*t136*t14*t37*t39*t93-7.073553026306428d-2*
     @         t14*t155*t34*t37*t39*t79-7.073553026306428d-2*t14*t153*t3
     1         7*t39*t79+2*t155*t39*t79+1.4147106052612857d-1*t14*t153*t
     2         40*t77+1.4147106052612857d-1*t14*t155*t34*t39*t77-7.07355
     3         3026306428d-2*t14*t155*t37*t39*t73-3.536776513153214d-2*t
     4         14*t214*t37*t40-3.536776513153214d-2*t14*t215*t34*t37*t39
     5         +t215*t39)+3.5449077018110403d+0*t39*t45*(5.6399136179197
     6         63d+0*t121*t22*t73*t93-7.049892022399704d+0*rhob*t120*t17
     7         6*t22*t93+2.8199568089598814d+0*rhob*t121*t131*t22*t93-1.
     8         8799712059732546d+0*rhob*t214*t22*t74-3.759942411946509d+
     9         0*t153*t22*t74+5.639913617919763d+0*rhob*t121*t153*t22*t7
     :         3)-1.329807601338106d-1*t10*t33*t73*t74*t80*t93-2.6596152
     ;         02676212d-1*t10*t37*t42*t80*t93+5.01149519714857d+1*t17*t
     <         36*t39*t45*t75*t78*t79*t93+5.01149519714857d+1*t132*t17*t
     =         36*t39*t45*t78*t93+2.659615202676212d-1*t10*t41*t42*t77*t
     >         93+1.329807601338106d-1*t10*t37*t41*t73*t74*t93-6.6490380
     ?         06690531d-2*t10*t131*t33*t41*t74*t93+1.329807601338106d-1
     @         *t10*t138*t33*t42*t93+9.973557010035795d-2*t10*t120*t121*
     1         t33*t41*t93+2.659615202676212d-1*t10*t153*t33*t42*t80+3.5
     2         449077018110403d+0*t154*t39*t45*t79+3.5449077018110403d+0
     3         *t155*t39*t45*t75-1.329807601338106d-1*t10*t153*t33*t41*t
     4         73*t74+2.659615202676212d-1*t10*t156*t33*t42*t73-2.659615
     5         202676212d-1*t10*t153*t37*t41*t42+1.329807601338106d-1*t1
     6         0*t214*t33*t41*t42-5.319230405352424d-1*t10*t156*t35*t37)
     7         +1.7730768017841417d-1*t10*t33*t73*t74*t81*t9*t93+3.54615
     8         36035682834d-1*t10*t37*t42*t81*t9*t93-3.5461536035682834d
     9         -1*t10*t42*t43*t77*t9*t93-1.7730768017841417d-1*t10*t37*t
     :         43*t73*t74*t9*t93+8.865384008920708d-2*t10*t131*t33*t43*t
     ;         74*t9*t93-1.3298076013381063d-1*t10*t120*t121*t33*t43*t9*
     <         t93-1.7730768017841417d-1*t10*t139*t33*t42*t9*t93-3.54615
     =         36035682834d-1*t10*t153*t33*t42*t81*t9+1.7730768017841417
     >         d-1*t10*t153*t33*t43*t73*t74*t9-3.5461536035682834d-1*t10
     ?         *t157*t33*t42*t73*t9+3.5461536035682834d-1*t10*t153*t37*t
     @         42*t43*t9-1.7730768017841417d-1*t10*t214*t33*t42*t43*t9+7
     1         .092307207136567d-1*t10*t157*t35*t37*t9)+t44*(1.306666666
     2         6666665d-2*gammabb*t122*t71*t92+2.2399999999999998d-2*gam
     3         mabb*t124*t65*t70*t92+2.52d-2*gammabb*t123*t203*t28*t92-8
     4         .4d-3*gammabb*t124*t129*t28*t92-1.1199999999999999d-2*t65
     5         *t70*t71-1.1199999999999999d-2*gammabb*t152*t65*t71+4.2d-
     6         3*gammabb*t213*t28*t71+4.2d-3*t129*t28*t71-1.679999999999
     7         9998d-2*gammabb*t124*t152*t28*t70-1.3066666666666665d-2*t
     8         122*t31-8.4d-3*t123*t124*t28)+2*t158*t83+2*t159*t82)*wght
     9         +Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t7*(-3.5461536035682834d-1*t10*t1
     1         1*t13*t9*(3.5449077018110403d+0*t19*t45*(5.63991361791976
     2         3d+0*rhoa*t144*t22*t86*t99+2.8199568089598814d+0*rhoa*t16
     3         8*t22*t54*t99+2.8199568089598814d+0*t161*t22*t99-1.879971
     4         2059732546d+0*rhoa*t217*t22*t55-1.8799712059732546d+0*t16
     5         8*t22*t55-7.049892022399704d+0*rhoa*t161*t162*t22*t54)+9.
     6         973557010035795d-2*t10*t11*t161*t21*t54*t99-1.32980760133
     7         8106d-1*t10*t11*t54*t55*t86*t87-2.659615202676212d-1*t10*
     8         t16*t23*t86*t87+2.659615202676212d-1*t10*t11*t144*t23*t87
     9         +2.659615202676212d-1*t10*t11*t13*(2.827433388230826d+1*t
     :         146*t15*t17*t19*t59*t86-1.0d+0*t144*t19*t59*t86-1.0d+0*t1
     ;         46*t18*t19*t86-7.073553026306428d-2*t14*t146*t16*t19*t86+
     <         1.413716694115413d+1*t15*t168*t17*t19*t59*t60-5.0d-1*t168
     =         *t18*t19*t60+1.9985948912206125d+2*t111*t161*t169*t170*t1
     >         9*t60-3.534291735288533d+1*t113*t15*t161*t17*t19*t60-3.53
     ?         6776513153214d-2*t14*t16*t168*t19*t60-5.0d-1*t168*t19*t54
     @         *t59+7.073553026306428d-2*t14*t168*t20*t58-7.068583470577
     1         066d+0*t15*t161*t17*t170*t19*t54+1.0d+0*t113*t161*t19*t54
     2         -3.536776513153214d-2*t12*t14*t16*t19*t218+t19*t218-3.536
     3         776513153214d-2*t14*t16*t20*t217+1.0d+0*t11*t168*t18*t19+
     4         1.413716694115413d+1*rhoa*t113*t161*t17*t19)+1.0022990394
     5         29714d+2*t145*t15*t17*t19*t45*t59*t86-1.329807601338106d-
     6         1*t10*t11*t144*t21*t55*t86+2.659615202676212d-1*t10*t11*t
     7         147*t23*t86-6.649038006690531d-2*t10*t11*t161*t55*t61+1.3
     8         29807601338106d-1*t10*t11*t168*t23*t61+5.01149519714857d+
     9         1*t15*t168*t17*t19*t45*t56*t59+7.084834422688148d+2*t111*
     :         t161*t169*t170*t19*t45*t56-1.002299039429714d+2*t113*t15*
     ;         t161*t17*t19*t45*t56-6.649038006690531d-2*t10*t11*t168*t2
     <         1*t54*t55+6.649038006690531d-2*t10*t16*t161*t21*t55+1.329
     =         807601338106d-1*t10*t11*t171*t23*t54+1.329807601338106d-1
     >         *t10*t11*t21*t217*t23-1.329807601338106d-1*t10*t16*t168*t
     ?         21*t23-2.659615202676212d-1*t10*t13*t16*t171)-1.329807601
     @         3381063d-1*t10*t11*t161*t24*t54*t9*t99+1.7730768017841417
     1         d-1*t10*t11*t54*t55*t86*t88*t9+3.5461536035682834d-1*t10*
     2         t16*t23*t86*t88*t9-3.5461536035682834d-1*t10*t11*t144*t23
     3         *t88*t9+1.7730768017841417d-1*t10*t11*t144*t24*t55*t86*t9
     4         -3.5461536035682834d-1*t10*t11*t148*t23*t86*t9+8.86538400
     5         8920708d-2*t10*t11*t161*t55*t62*t9-1.7730768017841417d-1*
     6         t10*t11*t168*t23*t62*t9+8.865384008920708d-2*t10*t11*t168
     7         *t24*t54*t55*t9-8.865384008920708d-2*t10*t16*t161*t24*t55
     8         *t9-1.7730768017841417d-1*t10*t11*t172*t23*t54*t9-1.77307
     9         68017841417d-1*t10*t11*t217*t23*t24*t9+1.7730768017841417
     :         d-1*t10*t16*t168*t23*t24*t9+3.5461536035682834d-1*t10*t13
     ;         *t16*t172*t9)+2*t149*t90+2*t150*t89+t25*(-1.1199999999999
     <         999d-2*t46*t52*t85-1.6799999999999998d-2*t102*t3*t51*t85-
     =         1.6799999999999998d-2*gammaaa*t102*t143*t3*t85-5.59999999
     >         9999999d-3*gammaaa*t167*t46*t52+4.2d-3*gammaaa*t216*t3*t5
     ?         2+8.4d-3*t143*t3*t52+2.52d-2*gammaaa*t163*t192*t3*t51-8.4
     @         d-3*gammaaa*t102*t167*t3*t51+1.1199999999999999d-2*gammaa
     1         a*t102*t163*t46)+t173*t64+t174*t63)*wght+Cmat3(iq,D3_RA_G
     2         AA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (2*t158*t97+2*t159*t96+t32*(1.7730
     1         768017841417d-1*t10*t33*t73*t74*t9*t93*t95+3.546153603568
     2         2834d-1*t10*t37*t42*t9*t93*t95-3.5461536035682834d-1*t10*
     3         t153*t33*t42*t9*t95-3.5461536035682834d-1*t10*t33*t35*t9*
     4         (-1.329807601338106d-1*t10*t33*t73*t74*t93*t94-2.65961520
     5         2676212d-1*t10*t37*t42*t93*t94+2.659615202676212d-1*t10*t
     6         153*t33*t42*t94+2.659615202676212d-1*t10*t33*t35*(2.82743
     7         3388230826d+1*t155*t17*t36*t39*t78*t93-1.0d+0*t153*t39*t7
     8         8*t93-1.0d+0*t155*t38*t39*t93-7.073553026306428d-2*t14*t1
     9         55*t37*t39*t93+1.413716694115413d+1*t17*t182*t36*t39*t78*
     :         t79-5.0d-1*t182*t38*t39*t79-3.536776513153214d-2*t14*t182
     ;         *t37*t39*t79-3.534291735288533d+1*t135*t17*t175*t36*t39*t
     <         79+1.9985948912206125d+2*t133*t169*t175*t183*t39*t79-5.0d
     =         -1*t182*t39*t73*t78+7.073553026306428d-2*t14*t182*t40*t77
     >         -7.068583470577066d+0*t17*t175*t183*t36*t39*t73+1.0d+0*t1
     ?         35*t175*t39*t73-3.536776513153214d-2*t14*t220*t37*t40+1.0
     @         d+0*t182*t33*t38*t39-3.536776513153214d-2*t14*t221*t34*t3
     1         7*t39+t221*t39+1.413716694115413d+1*rhob*t135*t17*t175*t3
     2         9)+3.5449077018110403d+0*t39*t45*(5.639913617919763d+0*rh
     3         ob*t121*t153*t22*t93-1.8799712059732546d+0*rhob*t22*t220*
     4         t74-1.8799712059732546d+0*t182*t22*t74+2.8199568089598814
     5         d+0*rhob*t121*t182*t22*t73-7.049892022399704d+0*rhob*t175
     6         *t176*t22*t73+2.8199568089598814d+0*t121*t175*t22)+1.0022
     7         99039429714d+2*t154*t17*t36*t39*t45*t78*t93-1.32980760133
     8         8106d-1*t10*t153*t33*t41*t74*t93+2.659615202676212d-1*t10
     9         *t156*t33*t42*t93-6.649038006690531d-2*t10*t175*t33*t74*t
     :         80+1.329807601338106d-1*t10*t182*t33*t42*t80+5.0114951971
     ;         4857d+1*t17*t182*t36*t39*t45*t75*t78-1.002299039429714d+2
     <         *t135*t17*t175*t36*t39*t45*t75+7.084834422688148d+2*t133*
     =         t169*t175*t183*t39*t45*t75-6.649038006690531d-2*t10*t182*
     >         t33*t41*t73*t74+6.649038006690531d-2*t10*t175*t37*t41*t74
     ?         +1.329807601338106d-1*t10*t184*t33*t42*t73+9.973557010035
     @         795d-2*t10*t121*t175*t33*t41*t73-1.329807601338106d-1*t10
     1         *t182*t37*t41*t42+1.329807601338106d-1*t10*t220*t33*t41*t
     2         42-2.659615202676212d-1*t10*t184*t35*t37)+1.7730768017841
     3         417d-1*t10*t153*t33*t43*t74*t9*t93-3.5461536035682834d-1*
     4         t10*t157*t33*t42*t9*t93+8.865384008920708d-2*t10*t175*t33
     5         *t74*t81*t9-1.7730768017841417d-1*t10*t182*t33*t42*t81*t9
     6         +8.865384008920708d-2*t10*t182*t33*t43*t73*t74*t9-8.86538
     7         4008920708d-2*t10*t175*t37*t43*t74*t9-1.3298076013381063d
     8         -1*t10*t121*t175*t33*t43*t73*t9-1.7730768017841417d-1*t10
     9         *t185*t33*t42*t73*t9+1.7730768017841417d-1*t10*t182*t37*t
     :         42*t43*t9-1.7730768017841417d-1*t10*t220*t33*t42*t43*t9+3
     ;         .5461536035682834d-1*t10*t185*t35*t37*t9)+t44*(-1.1199999
     <         999999999d-2*t65*t71*t92-1.6799999999999998d-2*t124*t28*t
     =         70*t92-1.6799999999999998d-2*gammabb*t124*t152*t28*t92-5.
     >         599999999999999d-3*gammabb*t181*t65*t71+4.2d-3*gammabb*t2
     ?         19*t28*t71+8.4d-3*t152*t28*t71+2.52d-2*gammabb*t177*t203*
     @         t28*t70-8.4d-3*gammabb*t124*t181*t28*t70+1.11999999999999
     1         99d-2*gammabb*t124*t177*t65)+t186*t83+t187*t82)*wght+Cmat
     2         3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t7*(-3.5461536035682834d-1*t10*t
     1         11*t13*t9*(2.998945983256911d+1*rhoa*t168*t19*t22*t45*t86
     2         *t99+9.973557010035795d-2*t10*t11*t21*t223*t99-1.99471140
     3         2007159d-1*t10*t11*t161*t55*t87+3.989422804014318d-1*t10*
     4         t11*t168*t23*t87+2.659615202676212d-1*t10*t11*t13*(5.9957
     5         84673661838d+2*t111*t168*t169*t170*t19*t86-1.060287520586
     6         5599d+2*t113*t15*t168*t17*t19*t86+1.413716694115413d+1*t1
     7         5*t17*t19*t226*t59+2.8254469624923584d+3*rhoa**6*t19*t223
     8         *t227/t12**6-3.536776513153214d-2*t14*t16*t20*t226-5.0d-1
     9         *t18*t19*t226+1.0602875205865599d+2*t15*t17*t170*t19*t223
     :         -1.2990866792933983d+3*t111*t169*t19*t223/t12**5)-2.82644
     ;         0000853771d+2*t160*t162*t168*t19*t45*t57*t86-1.9947114020
     <         07159d-1*t10*t11*t168*t21*t55*t86+3.989422804014318d-1*t1
     =         0*t11*t171*t23*t86+4.710733334756285d+2*t160*t19*t223*t45
     >         *t57/t13**9-6.66432440723758d+0*rhoa*t19*t22*t226*t45*t55
     ?         -1.3319284713741864d+3*t19*t196*t222*t223*t45/t13**11-2.4
     @         991216527140922d+1*rhoa*t162*t19*t22*t223*t45+1.329807601
     1         338106d-1*t10*t11*t21*t226*t23)-1.3298076013381063d-1*t10
     2         *t11*t223*t24*t9*t99+2.6596152026762127d-1*t10*t11*t161*t
     3         55*t88*t9-5.319230405352425d-1*t10*t11*t168*t23*t88*t9+2.
     4         6596152026762127d-1*t10*t11*t168*t24*t55*t86*t9-5.3192304
     5         05352425d-1*t10*t11*t172*t23*t86*t9-1.7730768017841417d-1
     6         *t10*t11*t226*t23*t24*t9)+3*t173*t90+3*t174*t89+t25*(-2.5
     7         2d-2*gammaaa*t102*t167*t3*t85+4.2d-3*gammaaa*t225*t3*t52+
     8         1.26d-2*t167*t3*t52+2.52d-2*gammaaa*t192*t224*t3-2.52d-2*
     9         t102*t163*t3))*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (3*t186*t97+3*t187*t96+t32*(2.659
     1         6152026762127d-1*t10*t175*t33*t74*t9*t95-5.31923040535242
     2         5d-1*t10*t182*t33*t42*t9*t95-3.5461536035682834d-1*t10*t3
     3         3*t35*t9*(-1.994711402007159d-1*t10*t175*t33*t74*t94+3.98
     4         9422804014318d-1*t10*t182*t33*t42*t94+2.659615202676212d-
     5         1*t10*t33*t35*(-1.0602875205865599d+2*t135*t17*t182*t36*t
     6         39*t93+5.995784673661838d+2*t133*t169*t182*t183*t39*t93+1
     7         .413716694115413d+1*t17*t231*t36*t39*t78-3.53677651315321
     8         4d-2*t14*t231*t37*t40-5.0d-1*t231*t38*t39+1.0602875205865
     9         599d+2*t17*t183*t228*t36*t39-1.2990866792933983d+3*t133*t
     :         169*t228*t39/t34**5+2.8254469624923584d+3*rhob**6*t227*t2
     ;         28*t39/t34**6)-2.826440000853771d+2*t160*t176*t182*t39*t4
     <         5*t76*t93-1.994711402007159d-1*t10*t182*t33*t41*t74*t93+2
     =         .998945983256911d+1*rhob*t121*t182*t22*t39*t45*t93+3.9894
     >         22804014318d-1*t10*t184*t33*t42*t93+4.710733334756285d+2*
     ?         t160*t228*t39*t45*t76/t35**9-6.66432440723758d+0*rhob*t22
     @         *t231*t39*t45*t74-1.3319284713741864d+3*t207*t222*t228*t3
     1         9*t45/t35**11-2.4991216527140922d+1*rhob*t176*t22*t228*t3
     2         9*t45+1.329807601338106d-1*t10*t231*t33*t41*t42+9.9735570
     3         10035795d-2*t10*t121*t228*t33*t41)+2.6596152026762127d-1*
     4         t10*t182*t33*t43*t74*t9*t93-5.319230405352425d-1*t10*t185
     5         *t33*t42*t9*t93-1.7730768017841417d-1*t10*t231*t33*t42*t4
     6         3*t9-1.3298076013381063d-1*t10*t121*t228*t33*t43*t9)+t44*
     7         (-2.52d-2*gammabb*t124*t181*t28*t92+4.2d-3*gammabb*t230*t
     8         28*t71+1.26d-2*t181*t28*t71+2.52d-2*gammabb*t203*t229*t28
     9         -2.52d-2*t124*t177*t28))*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = rhoa**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammaaa*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhoa
            t11 = 4.2d-3*gammaaa*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhoa**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhoa*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = (gammaaa*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhoa**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammaaa*t29*t30
            t32 = 1/t5**2
            t33 = rhoa**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammaaa*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammaaa*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         a*t21*t34*t35
            t37 = rhoa**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhoa*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammaaa*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammaaa*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammaaa*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammaaa*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = gammaaa**2
            t58 = 1/t28**3
            t59 = 1/rhoa**7.333333333333333d+0
            t60 = 1/rhoa**4.666666666666667d+0
            t61 = 1.6799999999999998d-1*gammaaa*t29*t60-4.47999999999999
     1         95d-2*t57*t58*t59+7.84d-2*t1*t4*t54
            t62 = 1/rhoa**6.666666666666666d-1
            t63 = 4.135669939329334d-1*t62-4.2d-3*gammaaa*t3*t32*t61+1.3
     1         066666666666665d-2*gammaaa*t54*t6+8.4d-3*gammaaa*t3*t55*t
     2         56+1.1199999999999999d-2*gammaaa*t26*t31*t32
            t64 = -1.8799712059732546d+0*rhoa*t21*t35*t63+2.819956808959
     1         8814d+0*rhoa*t21*t52*t53-3.759942411946509d+0*t21*t34*t35
            t65 = rhoa**4
            t66 = 1/t65
            t67 = 1/t11**3
            t68 = -2.827433388230826d+1*t14*t16*t52*t67+1.41371669411541
     1         3d+1*t14*t16*t39*t63+5.654866776461652d+1*rhoa*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t69 = t40**2
            t70 = -3.536776513153214d-2*t11*t13*t15*t18*t69+t18*t69-3.53
     1         6776513153214d-2*t11*t13*t15*t18*t68+t18*t68-2.1220659078
     2         919285d-1*t11*t13*t19*t66-3.536776513153214d-2*t13*t15*t1
     3         9*t63+1.4147106052612857d-1*t11*t13*t18*t38*t40-7.0735530
     4         26306428d-2*t13*t15*t18*t34*t40+1.4147106052612857d-1*t13
     5         *t19*t34*t38
            t71 = 2.659615202676212d-1*t10*t12*t70*t9+1.329807601338106d
     1         -1*t10*t20*t22*t63*t9-6.649038006690531d-2*t10*t20*t35*t5
     2         2*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.3192304053
     3         52424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12*t20*t38*
     4         t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.544907701811
     5         0403d+0*t18*t25*t64+3.5449077018110403d+0*t18*t25*t36*t40
            t72 = -3.5461536035682834d-1*t10*t12*t71*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t63*t8*t9+8.865384008920708d-2*t10*t
     2         23*t35*t52*t8*t9-3.5461536035682834d-1*t10*t22*t34*t42*t8
     3         *t9+7.092307207136567d-1*t12*t15*t42*t8*t9-7.092307207136
     4         567d-1*t12*t23*t38*t8*t9+3.5461536035682834d-1*t15*t22*t2
     5         3*t34*t8*t9
            t73 = -4.135669939329334d-1*t62+4.2d-3*gammaaa*t3*t32*t61-1.
     1         3066666666666665d-2*gammaaa*t54*t6-8.4d-3*gammaaa*t3*t55*
     2         t56-1.1199999999999999d-2*gammaaa*t26*t31*t32
            t74 = 1/rhoa**6.333333333333333d+0
            t75 = 1.6799999999999998d-2*gammaaa*t58*t74-1.67999999999999
     1         98d-2*t26*t4*t45-5.04d-2*t29*t30
            t76 = -4.2d-3*gammaaa*t3*t32*t75-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammaaa*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         aa*t26*t32*t46-4.2d-3*t3*t31*t32
            t77 = -1.8799712059732546d+0*rhoa*t21*t35*t76+2.819956808959
     1         8814d+0*rhoa*t21*t34*t47*t53-1.8799712059732546d+0*t21*t3
     2         5*t47
            t78 = 1.413716694115413d+1*t14*t16*t39*t76-2.827433388230826
     1         d+1*t14*t16*t34*t47*t67+2.827433388230826d+1*rhoa*t16*t39
     2         *t47
            t79 = -3.536776513153214d-2*t11*t13*t15*t18*t78+t18*t78-3.53
     1         6776513153214d-2*t13*t15*t19*t76+1.413716694115413d+1*t14
     2         *t16*t18*t39*t40*t47-5.0d-1*t17*t18*t40*t47-3.53677651315
     3         3214d-2*t13*t15*t18*t40*t47-5.0d-1*t18*t34*t39*t47+7.0735
     4         53026306428d-2*t13*t19*t38*t47+1.0d+0*t10*t17*t18*t47
            t80 = 2.659615202676212d-1*t10*t12*t79*t9+1.329807601338106d
     1         -1*t10*t20*t22*t76*t9+1.329807601338106d-1*t10*t22*t34*t4
     2         8*t9-2.659615202676212d-1*t12*t15*t48*t9+1.32980760133810
     3         6d-1*t10*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*
     4         t35*t47*t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.5449
     5         077018110403d+0*t18*t25*t77+5.01149519714857d+1*t14*t16*t
     6         18*t25*t36*t39*t47
            t81 = -3.5461536035682834d-1*t10*t12*t8*t80*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*
     2         t22*t34*t49*t8*t9+3.5461536035682834d-1*t12*t15*t49*t8*t9
     3         -1.7730768017841417d-1*t10*t22*t42*t47*t8*t9+8.8653840089
     4         20708d-2*t10*t23*t34*t35*t47*t8*t9+1.7730768017841417d-1*
     5         t15*t22*t23*t47*t8*t9
            t82 = 4.2d-3*gammaaa*t3*t32*t75+5.599999999999999d-3*t26*t6-
     1         8.4d-3*gammaaa*t3*t31*t46*t56-5.599999999999999d-3*gammaa
     2         a*t26*t32*t46+4.2d-3*t3*t31*t32
            t83 = 1/t9**3
            t84 = t47**2
            t85 = 1/t12**7
            t86 = t46**2
            t87 = 1/rhoa**5.333333333333333d+0
            t88 = 1/gammaaa
            t89 = 1/t1**3
            t90 = -6.3d-3*t3*t4*t89+6.3d-3*t27*t29*t88-6.3d-3*t58*t87
            t91 = -4.2d-3*gammaaa*t3*t32*t90+8.4d-3*gammaaa*t3*t56*t86-8
     1         .4d-3*t3*t32*t46
            t92 = 1/t9**4
            t93 = 1/t11**4
            t94 = 1.9985948912206125d+2*t18*t65*t84*t92*t93+1.4137166941
     1         15413d+1*t14*t16*t18*t39*t91-3.536776513153214d-2*t13*t15
     2         *t19*t91-5.0d-1*t17*t18*t91-3.534291735288533d+1*t14*t16*
     3         t18*t67*t84
            t95 = 2.659615202676212d-1*t10*t12*t9*t94+1.329807601338106d
     1         -1*t10*t20*t22*t9*t91-6.66432440723758d+0*rhoa*t18*t21*t2
     2         5*t35*t91-6.649038006690531d-2*t10*t20*t35*t84*t9+2.65961
     3         5202676212d-1*t10*t22*t47*t48*t9-9.421466669512572d+1*t18
     4         *t25*t37*t83*t84*t85+9.996486610856369d+0*rhoa*t18*t21*t2
     5         5*t53*t84
            t96 = -3.5461536035682834d-1*t10*t12*t8*t9*t95-1.77307680178
     1         41417d-1*t10*t22*t23*t8*t9*t91+8.865384008920708d-2*t10*t
     2         23*t35*t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8
     3         *t9
            t97 = 4.2d-3*gammaaa*t3*t32*t90-8.4d-3*gammaaa*t3*t56*t86+8.
     1         4d-3*t3*t32*t46
            t98 = t34**3
            t99 = 1/rhoa**4.333333333333333d+0
            t100 = 1/rhoa**1.6666666666666669d+0
            t101 = t31**3
            t102 = 1/t5**4
            t103 = 1/t28**5
            t104 = -2.6133333333333336d-1*t1*t4*t99+5.525333333333333d-1
     1         *t57*t58/rhoa**8.333333333333334d+0-8.885333333333333d-1*
     2         gammaaa*t29/rhoa**5.666666666666667d+0-1.7919999999999997
     3         d-1*gammaaa**3*t103/rhoa**11
            t105 = -4.355555555555556d-2*gammaaa*t6*t99+2.52d-2*gammaaa*
     1         t3*t31*t56*t61+1.6799999999999995d-2*gammaaa*t26*t32*t61-
     2         3.359999999999999d-2*gammaaa*t26*t55*t56-3.92d-2*gammaaa*
     3         t31*t32*t54-4.2d-3*gammaaa*t104*t3*t32-2.52d-2*gammaaa*t1
     4         01*t102*t3-2.7571132928862224d-1*t100
            t106 = rhoa**5
            t107 = 8.482300164692478d+1*t14*t16*t93*t98-8.48230016469247
     1         8d+1*t14*t16*t34*t63*t67-1.6964600329384957d+2*rhoa*t16*t
     2         52*t67+8.482300164692478d+1*rhoa*t16*t39*t63+8.4823001646
     3         92478d+1*t16*t34*t39+1.413716694115413d+1*t105*t14*t16*t3
     4         9
            t108 = t40**3
            t109 = 2.072d-1*t29*t60-1.7359999999999998d-1*gammaaa*t58*t5
     1         9+6.719999999999998d-2*t103*t57/rhoa**10+3.92d-2*t4*t45*t
     2         54
            t110 = 1.6799999999999998d-2*gammaaa*t3*t31*t56*t75+1.119999
     1         9999999999d-2*gammaaa*t26*t32*t75+8.4d-3*gammaaa*t3*t46*t
     2         56*t61-4.2d-3*t3*t32*t61+1.3066666666666665d-2*t54*t6+8.4
     3         d-3*t3*t55*t56-2.2399999999999998d-2*gammaaa*t26*t31*t46*
     4         t56-2.52d-2*gammaaa*t102*t3*t46*t55-1.3066666666666665d-2
     5         *gammaaa*t32*t46*t54+1.1199999999999999d-2*t26*t31*t32-4.
     6         2d-3*gammaaa*t109*t3*t32
            t111 = 8.482300164692478d+1*t14*t16*t47*t52*t93-5.6548667764
     1         61652d+1*t14*t16*t34*t67*t76+5.654866776461652d+1*rhoa*t1
     2         6*t39*t76-2.827433388230826d+1*t14*t16*t47*t63*t67-1.1309
     3         733552923305d+2*rhoa*t16*t34*t47*t67+2.827433388230826d+1
     4         *t16*t39*t47+1.413716694115413d+1*t110*t14*t16*t39
            t112 = 8.4d-3*t26*t4*t89-8.4d-3*t29*t30*t88+4.19999999999999
     1         93d-2*t58*t74-2.52d-2*gammaaa*t103/rhoa**9
            t113 = 8.4d-3*gammaaa*t3*t31*t56*t90+5.599999999999999d-3*ga
     1         mmaaa*t26*t32*t90-1.1199999999999999d-2*gammaaa*t26*t56*t
     2         86-2.52d-2*gammaaa*t102*t3*t31*t86+1.6799999999999998d-2*
     3         gammaaa*t3*t46*t56*t75-8.4d-3*t3*t32*t75+1.67999999999999
     4         98d-2*t3*t31*t46*t56+1.1199999999999999d-2*t26*t32*t46-4.
     5         2d-3*gammaaa*t112*t3*t32
            t114 = 8.482300164692478d+1*t14*t16*t34*t84*t93-2.8274333882
     1         30826d+1*t14*t16*t34*t67*t91+2.827433388230826d+1*rhoa*t1
     2         6*t39*t91-5.654866776461652d+1*rhoa*t16*t67*t84-5.6548667
     3         76461652d+1*t14*t16*t47*t67*t76+1.413716694115413d+1*t113
     4         *t14*t16*t39
            t115 = t47**3
            t116 = t46**3
            t117 = -3.15d-3*t58*t87*t88-9.450000000000001d-3*t27*t29/t57
     1         +9.450000000000001d-3*t3*t4/t1**5+9.450000000000001d-3*t1
     2         03/rhoa**8
            t118 = 2.52d-2*gammaaa*t3*t46*t56*t90-1.26d-2*t3*t32*t90+2.5
     1         2d-2*t3*t56*t86-4.2d-3*gammaaa*t117*t3*t32-2.52d-2*gammaa
     2         a*t102*t116*t3
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = t43*t7*wght+t24*t44*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t50*t7*wght+t24*t51*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t24*t73*wght+t7*t72*wght+2*t43*t44*wght
     1         +Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t24*t82*wght+t7*t81*wght+t43*t51*wght+
     1         t44*t50*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t24*t97*wght+t7*t96*wght+2*t50*t51*wg
     1         ht+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t24*(4.355555555555556d-2*gammaaa*t6
     1         *t99-2.52d-2*gammaaa*t3*t31*t56*t61-1.6799999999999995d-2
     2         *gammaaa*t26*t32*t61+3.359999999999999d-2*gammaaa*t26*t55
     3         *t56+3.92d-2*gammaaa*t31*t32*t54+4.2d-3*gammaaa*t104*t3*t
     4         32+2.52d-2*gammaaa*t101*t102*t3+2.7571132928862224d-1*t10
     5         0)*wght+t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(3.54490
     6         77018110403d+0*t18*t25*(-7.049892022399704d+0*rhoa*t21*t8
     7         5*t98+8.459870426879645d+0*rhoa*t21*t34*t53*t63-5.6399136
     8         17919763d+0*t21*t35*t63+8.459870426879645d+0*t21*t52*t53-
     9         1.8799712059732546d+0*rhoa*t105*t21*t35)+9.97355701003579
     :         5d-2*t10*t20*t53*t9*t98+3.989422804014318d-1*t10*t22*t34*
     ;         t70*t9-7.978845608028636d-1*t12*t15*t70*t9+2.659615202676
     <         212d-1*t10*t12*(2.1220659078919285d-1*t11*t13*t18*t38*t69
     =         -1.0610329539459643d-1*t13*t15*t18*t34*t69-1.061032953945
     >         9643d-1*t11*t13*t15*t18*t40*t68+3*t18*t40*t68+2.122065907
     ?         8919285d-1*t11*t13*t18*t38*t68-1.0610329539459643d-1*t13*
     @         t15*t18*t34*t68-6.366197723675786d-1*t11*t13*t18*t40*t66-
     1         6.366197723675786d-1*t13*t19*t34*t66-1.0610329539459643d-
     2         1*t13*t15*t18*t40*t63+2.1220659078919285d-1*t13*t19*t38*t
     3         63+4.244131815783857d-1*t13*t18*t34*t38*t40-3.53677651315
     4         3214d-2*t105*t13*t15*t19+8.488263631567714d-1*t11*t13*t19
     5         /t106-3.536776513153214d-2*t108*t11*t13*t15*t18-3.5367765
     6         13153214d-2*t107*t11*t13*t15*t18+t108*t18+t107*t18)*t9-1.
     7         5957691216057274d+0*t12*t20*t66*t9+3.989422804014318d-1*t
     8         10*t22*t41*t63*t9-1.994711402007159d-1*t10*t20*t34*t35*t6
     9         3*t9-3.989422804014318d-1*t15*t20*t22*t63*t9-1.9947114020
     :         07159d-1*t10*t35*t41*t52*t9+1.994711402007159d-1*t15*t20*
     ;         t35*t52*t9+1.5957691216057274d+0*t12*t38*t41*t9-7.9788456
     <         08028636d-1*t15*t22*t34*t41*t9+7.978845608028636d-1*t20*t
     =         22*t34*t38*t9+1.329807601338106d-1*t10*t105*t20*t22*t9+3.
     >         5449077018110403d+0*t18*t25*t36*t69+3.5449077018110403d+0
     ?         *t18*t25*t36*t68+7.089815403622081d+0*t18*t25*t40*t64)-1.
     @         3298076013381063d-1*t10*t23*t53*t8*t9*t98-5.3192304053524
     1         25d-1*t10*t22*t34*t71*t8*t9+1.063846081070485d+0*t12*t15*
     2         t71*t8*t9+2.12769216214097d+0*t12*t23*t66*t8*t9-5.3192304
     3         05352425d-1*t10*t22*t42*t63*t8*t9+2.6596152026762127d-1*t
     4         10*t23*t34*t35*t63*t8*t9+5.319230405352425d-1*t15*t22*t23
     5         *t63*t8*t9+2.6596152026762127d-1*t10*t35*t42*t52*t8*t9-2.
     6         6596152026762127d-1*t15*t23*t35*t52*t8*t9-2.1276921621409
     7         7d+0*t12*t38*t42*t8*t9+1.063846081070485d+0*t15*t22*t34*t
     8         42*t8*t9-1.063846081070485d+0*t22*t23*t34*t38*t8*t9-1.773
     9         0768017841417d-1*t10*t105*t22*t23*t8*t9)*wght+3*t43*t73*w
     :         ght+3*t44*t72*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t7*(-3.5461536035682834d-1*t10*t12*
     1         t8*t9*(2.659615202676212d-1*t10*t22*t34*t79*t9-5.31923040
     2         5352424d-1*t12*t15*t79*t9+2.659615202676212d-1*t10*t12*(-
     3         7.073553026306428d-2*t11*t13*t15*t18*t40*t78+2*t18*t40*t7
     4         8+1.4147106052612857d-1*t11*t13*t18*t38*t78-7.07355302630
     5         6428d-2*t13*t15*t18*t34*t78-7.073553026306428d-2*t13*t15*
     6         t18*t40*t76+1.4147106052612857d-1*t13*t19*t38*t76+1.41371
     7         6694115413d+1*t14*t16*t18*t39*t47*t69-5.0d-1*t17*t18*t47*
     8         t69-3.536776513153214d-2*t13*t15*t18*t47*t69+1.4137166941
     9         15413d+1*t14*t16*t18*t39*t47*t68-5.0d-1*t17*t18*t47*t68-3
     :         .536776513153214d-2*t13*t15*t18*t47*t68-2.122065907891928
     ;         5d-1*t13*t19*t47*t66-5.0d-1*t18*t39*t47*t63-1.0d+0*t18*t3
     <         4*t39*t40*t47+1.4147106052612857d-1*t13*t18*t38*t40*t47+2
     =         .0d+0*t10*t17*t18*t40*t47+2.0d+0*t10*t18*t34*t39*t47-3.0d
     >         +0*t15*t17*t18*t47-3.536776513153214d-2*t110*t13*t15*t19-
     ?         3.536776513153214d-2*t11*t111*t13*t15*t18+t111*t18)*t9+2.
     @         659615202676212d-1*t10*t22*t41*t76*t9-1.329807601338106d-
     1         1*t10*t20*t34*t35*t76*t9-2.659615202676212d-1*t15*t20*t22
     2         *t76*t9+1.329807601338106d-1*t10*t22*t47*t70*t9+1.3298076
     3         01338106d-1*t10*t22*t48*t63*t9-6.649038006690531d-2*t10*t
     4         20*t35*t47*t63*t9+9.973557010035795d-2*t10*t20*t47*t52*t5
     5         3*t9-6.649038006690531d-2*t10*t35*t48*t52*t9+5.3192304053
     6         52424d-1*t12*t38*t48*t9-2.659615202676212d-1*t15*t22*t34*
     7         t48*t9-1.329807601338106d-1*t10*t34*t35*t41*t47*t9-2.6596
     8         15202676212d-1*t15*t22*t41*t47*t9+2.659615202676212d-1*t2
     9         0*t22*t38*t47*t9+1.329807601338106d-1*t15*t20*t34*t35*t47
     :         *t9+1.329807601338106d-1*t10*t110*t20*t22*t9+3.5449077018
     ;         110403d+0*t18*t25*(-7.049892022399704d+0*rhoa*t21*t47*t52
     <         *t85+5.639913617919763d+0*rhoa*t21*t34*t53*t76-3.75994241
     =         1946509d+0*t21*t35*t76+2.8199568089598814d+0*rhoa*t21*t47
     >         *t53*t63+5.639913617919763d+0*t21*t34*t47*t53-1.879971205
     ?         9732546d+0*rhoa*t110*t21*t35)+3.5449077018110403d+0*t18*t
     @         25*t36*t78+3.5449077018110403d+0*t18*t25*t40*t77+5.011495
     1         19714857d+1*t14*t16*t18*t25*t39*t47*t64+5.01149519714857d
     2         +1*t14*t16*t18*t25*t36*t39*t40*t47)-3.5461536035682834d-1
     3         *t10*t22*t34*t8*t80*t9+7.092307207136567d-1*t12*t15*t8*t8
     4         0*t9-3.5461536035682834d-1*t10*t22*t42*t76*t8*t9+1.773076
     5         8017841417d-1*t10*t23*t34*t35*t76*t8*t9+3.546153603568283
     6         4d-1*t15*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*t22*
     7         t47*t71*t8*t9-1.7730768017841417d-1*t10*t22*t49*t63*t8*t9
     8         +8.865384008920708d-2*t10*t23*t35*t47*t63*t8*t9-1.3298076
     9         013381063d-1*t10*t23*t47*t52*t53*t8*t9+8.865384008920708d
     :         -2*t10*t35*t49*t52*t8*t9-7.092307207136567d-1*t12*t38*t49
     ;         *t8*t9+3.5461536035682834d-1*t15*t22*t34*t49*t8*t9+1.7730
     <         768017841417d-1*t10*t34*t35*t42*t47*t8*t9+3.5461536035682
     =         834d-1*t15*t22*t42*t47*t8*t9-3.5461536035682834d-1*t22*t2
     >         3*t38*t47*t8*t9-1.7730768017841417d-1*t15*t23*t34*t35*t47
     ?         *t8*t9-1.7730768017841417d-1*t10*t110*t22*t23*t8*t9)*wght
     @         +2*t43*t82*wght+2*t44*t81*wght+t24*(-1.6799999999999998d-
     1         2*gammaaa*t3*t31*t56*t75-1.1199999999999999d-2*gammaaa*t2
     2         6*t32*t75-8.4d-3*gammaaa*t3*t46*t56*t61+4.2d-3*t3*t32*t61
     3         -1.3066666666666665d-2*t54*t6-8.4d-3*t3*t55*t56+2.2399999
     4         999999998d-2*gammaaa*t26*t31*t46*t56+2.52d-2*gammaaa*t102
     5         *t3*t46*t55+1.3066666666666665d-2*gammaaa*t32*t46*t54-1.1
     6         199999999999999d-2*t26*t31*t32+4.2d-3*gammaaa*t109*t3*t32
     7         )*wght+t50*t73*wght+t51*t72*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t43*t97*wght+t44*t96*wght+t7*(-1.7
     1         730768017841417d-1*t10*t22*t34*t8*t9*t95+3.54615360356828
     2         34d-1*t12*t15*t8*t9*t95-3.5461536035682834d-1*t10*t12*t8*
     3         t9*(1.329807601338106d-1*t10*t22*t34*t9*t94-2.65961520267
     4         6212d-1*t12*t15*t9*t94+2.659615202676212d-1*t10*t12*t9*(1
     5         .9985948912206125d+2*t18*t40*t65*t84*t92*t93-7.0685834705
     6         77066d+0*t14*t16*t18*t34*t84*t93+1.413716694115413d+1*t14
     7         *t16*t18*t39*t40*t91-5.0d-1*t17*t18*t40*t91-3.53677651315
     8         3214d-2*t13*t15*t18*t40*t91-5.0d-1*t18*t34*t39*t91+7.0735
     9         53026306428d-2*t13*t19*t38*t91+1.0d+0*t10*t17*t18*t91-3.5
     :         34291735288533d+1*t14*t16*t18*t40*t67*t84+1.0d+0*t18*t34*
     ;         t67*t84+1.413716694115413d+1*rhoa*t16*t18*t67*t84+2.82743
     <         3388230826d+1*t14*t16*t18*t39*t47*t78-1.0d+0*t17*t18*t47*
     =         t78-7.073553026306428d-2*t13*t15*t18*t47*t78-1.0d+0*t18*t
     >         39*t47*t76-3.536776513153214d-2*t113*t13*t15*t19-3.536776
     ?         513153214d-2*t11*t114*t13*t15*t18+t114*t18)+7.08483442268
     @         8148d+2*t18*t25*t36*t65*t84*t92*t93+3.5449077018110403d+0
     1         *t18*t25*(2.8199568089598814d+0*rhoa*t21*t34*t53*t91-1.87
     2         99712059732546d+0*t21*t35*t91-7.049892022399704d+0*rhoa*t
     3         21*t34*t84*t85+2.8199568089598814d+0*t21*t53*t84+5.639913
     4         617919763d+0*rhoa*t21*t47*t53*t76-1.8799712059732546d+0*r
     5         hoa*t113*t21*t35)+1.329807601338106d-1*t10*t22*t41*t9*t91
     6         -6.649038006690531d-2*t10*t20*t34*t35*t9*t91-1.3298076013
     7         38106d-1*t15*t20*t22*t9*t91+5.01149519714857d+1*t14*t16*t
     8         18*t25*t36*t39*t91+9.973557010035795d-2*t10*t20*t34*t53*t
     9         84*t9-6.649038006690531d-2*t10*t35*t41*t84*t9+6.649038006
     :         690531d-2*t15*t20*t35*t84*t9+2.659615202676212d-1*t10*t22
     ;         *t47*t79*t9+2.659615202676212d-1*t10*t22*t48*t76*t9-1.329
     <         807601338106d-1*t10*t20*t35*t47*t76*t9-1.329807601338106d
     =         -1*t10*t34*t35*t47*t48*t9-2.659615202676212d-1*t15*t22*t4
     >         7*t48*t9+1.329807601338106d-1*t10*t113*t20*t22*t9-1.00229
     ?         9039429714d+2*t14*t16*t18*t25*t36*t67*t84+1.0022990394297
     @         14d+2*t14*t16*t18*t25*t39*t47*t77)-1.7730768017841417d-1*
     1         t10*t22*t42*t8*t9*t91+8.865384008920708d-2*t10*t23*t34*t3
     2         5*t8*t9*t91+1.7730768017841417d-1*t15*t22*t23*t8*t9*t91-1
     3         .3298076013381063d-1*t10*t23*t34*t53*t8*t84*t9+8.86538400
     4         8920708d-2*t10*t35*t42*t8*t84*t9-8.865384008920708d-2*t15
     5         *t23*t35*t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t8*t
     6         80*t9-3.5461536035682834d-1*t10*t22*t49*t76*t8*t9+1.77307
     7         68017841417d-1*t10*t23*t35*t47*t76*t8*t9+1.77307680178414
     8         17d-1*t10*t34*t35*t47*t49*t8*t9+3.5461536035682834d-1*t15
     9         *t22*t47*t49*t8*t9-1.7730768017841417d-1*t10*t113*t22*t23
     :         *t8*t9)*wght+t24*(-8.4d-3*gammaaa*t3*t31*t56*t90-5.599999
     ;         999999999d-3*gammaaa*t26*t32*t90+1.1199999999999999d-2*ga
     <         mmaaa*t26*t56*t86+2.52d-2*gammaaa*t102*t3*t31*t86-1.67999
     =         99999999998d-2*gammaaa*t3*t46*t56*t75+8.4d-3*t3*t32*t75-1
     >         .6799999999999998d-2*t3*t31*t46*t56-1.1199999999999999d-2
     ?         *t26*t32*t46+4.2d-3*gammaaa*t112*t3*t32)*wght+2*t50*t82*w
     @         ght+2*t51*t81*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 3*t50*t97*wght+3*t51*t96*wght+t7*
     1         (-5.319230405352425d-1*t10*t22*t47*t8*t9*t95-3.5461536035
     2         682834d-1*t10*t12*t8*t9*(3.989422804014318d-1*t10*t22*t47
     3         *t9*t94+2.659615202676212d-1*t10*t12*t9*(5.99578467366183
     4         8d+2*t18*t47*t65*t91*t92*t93+1.0602875205865599d+2*t115*t
     5         14*t16*t18*t93-1.2990866792933983d+3*t115*t18*t65*t92/t11
     6         **5-1.0602875205865599d+2*t14*t16*t18*t47*t67*t91+2.82544
     7         69624923584d+3*rhoa**6*t115*t18/(t11**6*t9**6)+1.41371669
     8         4115413d+1*t118*t14*t16*t18*t39-3.536776513153214d-2*t118
     9         *t13*t15*t19-5.0d-1*t118*t17*t18)+3.989422804014318d-1*t1
     :         0*t22*t48*t9*t91-1.994711402007159d-1*t10*t20*t35*t47*t9*
     ;         t91-2.826440000853771d+2*t18*t25*t37*t47*t83*t85*t91+2.99
     <         8945983256911d+1*rhoa*t18*t21*t25*t47*t53*t91-1.994711402
     =         007159d-1*t10*t35*t48*t84*t9+9.973557010035795d-2*t10*t11
     >         5*t20*t53*t9+1.329807601338106d-1*t10*t118*t20*t22*t9-1.3
     ?         319284713741864d+3*t106*t115*t18*t25/(t12**11*t9**5)-2.49
     @         91216527140922d+1*rhoa*t115*t18*t21*t25*t85+4.71073333475
     1         6285d+2*t115*t18*t25*t37*t83/t12**9-6.66432440723758d+0*r
     2         hoa*t118*t18*t21*t25*t35)-5.319230405352425d-1*t10*t22*t4
     3         9*t8*t9*t91+2.6596152026762127d-1*t10*t23*t35*t47*t8*t9*t
     4         91+2.6596152026762127d-1*t10*t35*t49*t8*t84*t9-1.32980760
     5         13381063d-1*t10*t115*t23*t53*t8*t9-1.7730768017841417d-1*
     6         t10*t118*t22*t23*t8*t9)*wght+t24*(-2.52d-2*gammaaa*t3*t46
     7         *t56*t90+1.26d-2*t3*t32*t90-2.52d-2*t3*t56*t86+4.2d-3*gam
     8         maaa*t117*t3*t32+2.52d-2*gammaaa*t102*t116*t3)*wght+Cmat3
     9         (iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = rhob**1.3333333333333333d+0
            t3 = 1/t2
            t4 = asinh(t1*t3)
            t5 = 2.52d-2*t1*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = -4.2d-3*gammabb*t3*t6-9.305257363491002d-1*t2
            t8 = param(2)
            t9 = param(3)
            t10 = 1/rhob
            t11 = 4.2d-3*gammabb*t3*t6+9.305257363491002d-1*t2
            t12 = t11**5.0d-1
            t13 = t9**2
            t14 = rhob**2
            t15 = 1/t14
            t16 = 1/t13
            t17 = 1/t11
            t18 = exp(-1.413716694115413d+1*t14*t16*t17)
            t19 = t18-1.0d+0
            t20 = -3.536776513153214d-2*t11*t13*t15*t19+t18-1.5d+0
            t21 = 1/t9
            t22 = 1/t12
            t23 = 2.659615202676212d-1*t10*t12*t20*t9+1.7724538509055202
     1         d+0*erf(3.759942411946509d+0*rhob*t21*t22)
            t24 = -3.5461536035682834d-1*t10*t12*t23*t8*t9-param(1)+1.0d
     1         +0
            t25 = 5.641895835477563d-1
            t26 = 1/rhob**2.3333333333333334d+0
            t27 = 1/rhob**2.6666666666666666d+0
            t28 = (gammabb*t27+1)**5.0d-1
            t29 = 1/t28
            t30 = 1/rhob**3.6666666666666664d+0
            t31 = -3.3599999999999997d-2*t1*t26*t4-3.3599999999999997d-2
     1         *gammabb*t29*t30
            t32 = 1/t5**2
            t33 = rhob**3.333333333333333d-1
            t34 = -5.599999999999999d-3*gammabb*t26*t6+1.240700981798800
     1         2d+0*t33-4.2d-3*gammabb*t3*t31*t32
            t35 = 1/t12**3
            t36 = 3.759942411946509d+0*t21*t22-1.8799712059732546d+0*rho
     1         b*t21*t34*t35
            t37 = rhob**3
            t38 = 1/t37
            t39 = 1/t11**2
            t40 = 1.413716694115413d+1*t14*t16*t34*t39-2.827433388230826
     1         d+1*rhob*t16*t17
            t41 = -3.536776513153214d-2*t11*t13*t15*t18*t40+t18*t40+7.07
     1         3553026306428d-2*t11*t13*t19*t38-3.536776513153214d-2*t13
     2         *t15*t19*t34
            t42 = 2.659615202676212d-1*t10*t12*t41*t9+1.329807601338106d
     1         -1*t10*t20*t22*t34*t9-2.659615202676212d-1*t12*t15*t20*t9
     2         +3.5449077018110403d+0*t18*t25*t36
            t43 = -3.5461536035682834d-1*t10*t12*t42*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t34*t8*t9+3.5461536035682834d-1*t12*
     2         t15*t23*t8*t9
            t44 = 5.599999999999999d-3*gammabb*t26*t6-1.2407009817988002
     1         d+0*t33+4.2d-3*gammabb*t3*t31*t32
            t45 = 1/t1
            t46 = 1.26d-2*t3*t4*t45+1.26d-2*t27*t29
            t47 = 4.2d-3*t3*t6-4.2d-3*gammabb*t3*t32*t46
            t48 = 1.413716694115413d+1*t14*t16*t18*t39*t47-3.53677651315
     1         3214d-2*t13*t15*t19*t47-5.0d-1*t17*t18*t47
            t49 = 2.659615202676212d-1*t10*t12*t48*t9+1.329807601338106d
     1         -1*t10*t20*t22*t47*t9-6.66432440723758d+0*rhob*t18*t21*t2
     2         5*t35*t47
            t50 = -3.5461536035682834d-1*t10*t12*t49*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t47*t8*t9
            t51 = 4.2d-3*gammabb*t3*t32*t46-4.2d-3*t3*t6
            t52 = t34**2
            t53 = 1/t12**5
            t54 = 1/rhob**3.3333333333333337d+0
            t55 = t31**2
            t56 = 1/t5**3
            t57 = gammabb**2
            t58 = 1/t28**3
            t59 = 1/rhob**7.333333333333333d+0
            t60 = 1/rhob**4.666666666666667d+0
            t61 = 1.6799999999999998d-1*gammabb*t29*t60-4.47999999999999
     1         95d-2*t57*t58*t59+7.84d-2*t1*t4*t54
            t62 = 1/rhob**6.666666666666666d-1
            t63 = 4.135669939329334d-1*t62-4.2d-3*gammabb*t3*t32*t61+1.3
     1         066666666666665d-2*gammabb*t54*t6+8.4d-3*gammabb*t3*t55*t
     2         56+1.1199999999999999d-2*gammabb*t26*t31*t32
            t64 = -1.8799712059732546d+0*rhob*t21*t35*t63+2.819956808959
     1         8814d+0*rhob*t21*t52*t53-3.759942411946509d+0*t21*t34*t35
            t65 = rhob**4
            t66 = 1/t65
            t67 = 1/t11**3
            t68 = -2.827433388230826d+1*t14*t16*t52*t67+1.41371669411541
     1         3d+1*t14*t16*t39*t63+5.654866776461652d+1*rhob*t16*t34*t3
     2         9-2.827433388230826d+1*t16*t17
            t69 = t40**2
            t70 = -3.536776513153214d-2*t11*t13*t15*t18*t69+t18*t69-3.53
     1         6776513153214d-2*t11*t13*t15*t18*t68+t18*t68-2.1220659078
     2         919285d-1*t11*t13*t19*t66-3.536776513153214d-2*t13*t15*t1
     3         9*t63+1.4147106052612857d-1*t11*t13*t18*t38*t40-7.0735530
     4         26306428d-2*t13*t15*t18*t34*t40+1.4147106052612857d-1*t13
     5         *t19*t34*t38
            t71 = 2.659615202676212d-1*t10*t12*t70*t9+1.329807601338106d
     1         -1*t10*t20*t22*t63*t9-6.649038006690531d-2*t10*t20*t35*t5
     2         2*t9+2.659615202676212d-1*t10*t22*t34*t41*t9-5.3192304053
     3         52424d-1*t12*t15*t41*t9+5.319230405352424d-1*t12*t20*t38*
     4         t9-2.659615202676212d-1*t15*t20*t22*t34*t9+3.544907701811
     5         0403d+0*t18*t25*t64+3.5449077018110403d+0*t18*t25*t36*t40
            t72 = -3.5461536035682834d-1*t10*t12*t71*t8*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t63*t8*t9+8.865384008920708d-2*t10*t
     2         23*t35*t52*t8*t9-3.5461536035682834d-1*t10*t22*t34*t42*t8
     3         *t9+7.092307207136567d-1*t12*t15*t42*t8*t9-7.092307207136
     4         567d-1*t12*t23*t38*t8*t9+3.5461536035682834d-1*t15*t22*t2
     5         3*t34*t8*t9
            t73 = -4.135669939329334d-1*t62+4.2d-3*gammabb*t3*t32*t61-1.
     1         3066666666666665d-2*gammabb*t54*t6-8.4d-3*gammabb*t3*t55*
     2         t56-1.1199999999999999d-2*gammabb*t26*t31*t32
            t74 = 1/rhob**6.333333333333333d+0
            t75 = 1.6799999999999998d-2*gammabb*t58*t74-1.67999999999999
     1         98d-2*t26*t4*t45-5.04d-2*t29*t30
            t76 = -4.2d-3*gammabb*t3*t32*t75-5.599999999999999d-3*t26*t6
     1         +8.4d-3*gammabb*t3*t31*t46*t56+5.599999999999999d-3*gamma
     2         bb*t26*t32*t46-4.2d-3*t3*t31*t32
            t77 = -1.8799712059732546d+0*rhob*t21*t35*t76+2.819956808959
     1         8814d+0*rhob*t21*t34*t47*t53-1.8799712059732546d+0*t21*t3
     2         5*t47
            t78 = 1.413716694115413d+1*t14*t16*t39*t76-2.827433388230826
     1         d+1*t14*t16*t34*t47*t67+2.827433388230826d+1*rhob*t16*t39
     2         *t47
            t79 = -3.536776513153214d-2*t11*t13*t15*t18*t78+t18*t78-3.53
     1         6776513153214d-2*t13*t15*t19*t76+1.413716694115413d+1*t14
     2         *t16*t18*t39*t40*t47-5.0d-1*t17*t18*t40*t47-3.53677651315
     3         3214d-2*t13*t15*t18*t40*t47-5.0d-1*t18*t34*t39*t47+7.0735
     4         53026306428d-2*t13*t19*t38*t47+1.0d+0*t10*t17*t18*t47
            t80 = 2.659615202676212d-1*t10*t12*t79*t9+1.329807601338106d
     1         -1*t10*t20*t22*t76*t9+1.329807601338106d-1*t10*t22*t34*t4
     2         8*t9-2.659615202676212d-1*t12*t15*t48*t9+1.32980760133810
     3         6d-1*t10*t22*t41*t47*t9-6.649038006690531d-2*t10*t20*t34*
     4         t35*t47*t9-1.329807601338106d-1*t15*t20*t22*t47*t9+3.5449
     5         077018110403d+0*t18*t25*t77+5.01149519714857d+1*t14*t16*t
     6         18*t25*t36*t39*t47
            t81 = -3.5461536035682834d-1*t10*t12*t8*t80*t9-1.77307680178
     1         41417d-1*t10*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*
     2         t22*t34*t49*t8*t9+3.5461536035682834d-1*t12*t15*t49*t8*t9
     3         -1.7730768017841417d-1*t10*t22*t42*t47*t8*t9+8.8653840089
     4         20708d-2*t10*t23*t34*t35*t47*t8*t9+1.7730768017841417d-1*
     5         t15*t22*t23*t47*t8*t9
            t82 = 4.2d-3*gammabb*t3*t32*t75+5.599999999999999d-3*t26*t6-
     1         8.4d-3*gammabb*t3*t31*t46*t56-5.599999999999999d-3*gammab
     2         b*t26*t32*t46+4.2d-3*t3*t31*t32
            t83 = 1/t9**3
            t84 = t47**2
            t85 = 1/t12**7
            t86 = t46**2
            t87 = 1/rhob**5.333333333333333d+0
            t88 = 1/gammabb
            t89 = 1/t1**3
            t90 = -6.3d-3*t3*t4*t89+6.3d-3*t27*t29*t88-6.3d-3*t58*t87
            t91 = -4.2d-3*gammabb*t3*t32*t90+8.4d-3*gammabb*t3*t56*t86-8
     1         .4d-3*t3*t32*t46
            t92 = 1/t9**4
            t93 = 1/t11**4
            t94 = 1.9985948912206125d+2*t18*t65*t84*t92*t93+1.4137166941
     1         15413d+1*t14*t16*t18*t39*t91-3.536776513153214d-2*t13*t15
     2         *t19*t91-5.0d-1*t17*t18*t91-3.534291735288533d+1*t14*t16*
     3         t18*t67*t84
            t95 = 2.659615202676212d-1*t10*t12*t9*t94+1.329807601338106d
     1         -1*t10*t20*t22*t9*t91-6.66432440723758d+0*rhob*t18*t21*t2
     2         5*t35*t91-6.649038006690531d-2*t10*t20*t35*t84*t9+2.65961
     3         5202676212d-1*t10*t22*t47*t48*t9-9.421466669512572d+1*t18
     4         *t25*t37*t83*t84*t85+9.996486610856369d+0*rhob*t18*t21*t2
     5         5*t53*t84
            t96 = -3.5461536035682834d-1*t10*t12*t8*t9*t95-1.77307680178
     1         41417d-1*t10*t22*t23*t8*t9*t91+8.865384008920708d-2*t10*t
     2         23*t35*t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t49*t8
     3         *t9
            t97 = 4.2d-3*gammabb*t3*t32*t90-8.4d-3*gammabb*t3*t56*t86+8.
     1         4d-3*t3*t32*t46
            t98 = t34**3
            t99 = 1/rhob**4.333333333333333d+0
            t100 = 1/rhob**1.6666666666666669d+0
            t101 = t31**3
            t102 = 1/t5**4
            t103 = 1/t28**5
            t104 = -2.6133333333333336d-1*t1*t4*t99+5.525333333333333d-1
     1         *t57*t58/rhob**8.333333333333334d+0-8.885333333333333d-1*
     2         gammabb*t29/rhob**5.666666666666667d+0-1.7919999999999997
     3         d-1*gammabb**3*t103/rhob**11
            t105 = -4.355555555555556d-2*gammabb*t6*t99+2.52d-2*gammabb*
     1         t3*t31*t56*t61+1.6799999999999995d-2*gammabb*t26*t32*t61-
     2         3.359999999999999d-2*gammabb*t26*t55*t56-3.92d-2*gammabb*
     3         t31*t32*t54-4.2d-3*gammabb*t104*t3*t32-2.52d-2*gammabb*t1
     4         01*t102*t3-2.7571132928862224d-1*t100
            t106 = rhob**5
            t107 = 8.482300164692478d+1*t14*t16*t93*t98-8.48230016469247
     1         8d+1*t14*t16*t34*t63*t67-1.6964600329384957d+2*rhob*t16*t
     2         52*t67+8.482300164692478d+1*rhob*t16*t39*t63+8.4823001646
     3         92478d+1*t16*t34*t39+1.413716694115413d+1*t105*t14*t16*t3
     4         9
            t108 = t40**3
            t109 = 2.072d-1*t29*t60-1.7359999999999998d-1*gammabb*t58*t5
     1         9+6.719999999999998d-2*t103*t57/rhob**10+3.92d-2*t4*t45*t
     2         54
            t110 = 1.6799999999999998d-2*gammabb*t3*t31*t56*t75+1.119999
     1         9999999999d-2*gammabb*t26*t32*t75+8.4d-3*gammabb*t3*t46*t
     2         56*t61-4.2d-3*t3*t32*t61+1.3066666666666665d-2*t54*t6+8.4
     3         d-3*t3*t55*t56-2.2399999999999998d-2*gammabb*t26*t31*t46*
     4         t56-2.52d-2*gammabb*t102*t3*t46*t55-1.3066666666666665d-2
     5         *gammabb*t32*t46*t54+1.1199999999999999d-2*t26*t31*t32-4.
     6         2d-3*gammabb*t109*t3*t32
            t111 = 8.482300164692478d+1*t14*t16*t47*t52*t93-5.6548667764
     1         61652d+1*t14*t16*t34*t67*t76+5.654866776461652d+1*rhob*t1
     2         6*t39*t76-2.827433388230826d+1*t14*t16*t47*t63*t67-1.1309
     3         733552923305d+2*rhob*t16*t34*t47*t67+2.827433388230826d+1
     4         *t16*t39*t47+1.413716694115413d+1*t110*t14*t16*t39
            t112 = 8.4d-3*t26*t4*t89-8.4d-3*t29*t30*t88+4.19999999999999
     1         93d-2*t58*t74-2.52d-2*gammabb*t103/rhob**9
            t113 = 8.4d-3*gammabb*t3*t31*t56*t90+5.599999999999999d-3*ga
     1         mmabb*t26*t32*t90-1.1199999999999999d-2*gammabb*t26*t56*t
     2         86-2.52d-2*gammabb*t102*t3*t31*t86+1.6799999999999998d-2*
     3         gammabb*t3*t46*t56*t75-8.4d-3*t3*t32*t75+1.67999999999999
     4         98d-2*t3*t31*t46*t56+1.1199999999999999d-2*t26*t32*t46-4.
     5         2d-3*gammabb*t112*t3*t32
            t114 = 8.482300164692478d+1*t14*t16*t34*t84*t93-2.8274333882
     1         30826d+1*t14*t16*t34*t67*t91+2.827433388230826d+1*rhob*t1
     2         6*t39*t91-5.654866776461652d+1*rhob*t16*t67*t84-5.6548667
     3         76461652d+1*t14*t16*t47*t67*t76+1.413716694115413d+1*t113
     4         *t14*t16*t39
            t115 = t47**3
            t116 = t46**3
            t117 = -3.15d-3*t58*t87*t88-9.450000000000001d-3*t27*t29/t57
     1         +9.450000000000001d-3*t3*t4/t1**5+9.450000000000001d-3*t1
     2         03/rhob**8
            t118 = 2.52d-2*gammabb*t3*t46*t56*t90-1.26d-2*t3*t32*t90+2.5
     1         2d-2*t3*t56*t86-4.2d-3*gammabb*t117*t3*t32-2.52d-2*gammab
     2         b*t102*t116*t3
            fnc(iq) = t24*t7*wght+fnc(iq)
            Amat(iq,D1_RB) = t43*t7*wght+t24*t44*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t50*t7*wght+t24*t51*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t24*t73*wght+t7*t72*wght+2*t43*t44*wght
     1         +Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t24*t82*wght+t7*t81*wght+t43*t51*wght+
     1         t44*t50*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t24*t97*wght+t7*t96*wght+2*t50*t51*wg
     1         ht+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t24*(4.355555555555556d-2*gammabb*t6
     1         *t99-2.52d-2*gammabb*t3*t31*t56*t61-1.6799999999999995d-2
     2         *gammabb*t26*t32*t61+3.359999999999999d-2*gammabb*t26*t55
     3         *t56+3.92d-2*gammabb*t31*t32*t54+4.2d-3*gammabb*t104*t3*t
     4         32+2.52d-2*gammabb*t101*t102*t3+2.7571132928862224d-1*t10
     5         0)*wght+t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(3.54490
     6         77018110403d+0*t18*t25*(-7.049892022399704d+0*rhob*t21*t8
     7         5*t98+8.459870426879645d+0*rhob*t21*t34*t53*t63-5.6399136
     8         17919763d+0*t21*t35*t63+8.459870426879645d+0*t21*t52*t53-
     9         1.8799712059732546d+0*rhob*t105*t21*t35)+9.97355701003579
     :         5d-2*t10*t20*t53*t9*t98+3.989422804014318d-1*t10*t22*t34*
     ;         t70*t9-7.978845608028636d-1*t12*t15*t70*t9+2.659615202676
     <         212d-1*t10*t12*(2.1220659078919285d-1*t11*t13*t18*t38*t69
     =         -1.0610329539459643d-1*t13*t15*t18*t34*t69-1.061032953945
     >         9643d-1*t11*t13*t15*t18*t40*t68+3*t18*t40*t68+2.122065907
     ?         8919285d-1*t11*t13*t18*t38*t68-1.0610329539459643d-1*t13*
     @         t15*t18*t34*t68-6.366197723675786d-1*t11*t13*t18*t40*t66-
     1         6.366197723675786d-1*t13*t19*t34*t66-1.0610329539459643d-
     2         1*t13*t15*t18*t40*t63+2.1220659078919285d-1*t13*t19*t38*t
     3         63+4.244131815783857d-1*t13*t18*t34*t38*t40-3.53677651315
     4         3214d-2*t105*t13*t15*t19+8.488263631567714d-1*t11*t13*t19
     5         /t106-3.536776513153214d-2*t108*t11*t13*t15*t18-3.5367765
     6         13153214d-2*t107*t11*t13*t15*t18+t108*t18+t107*t18)*t9-1.
     7         5957691216057274d+0*t12*t20*t66*t9+3.989422804014318d-1*t
     8         10*t22*t41*t63*t9-1.994711402007159d-1*t10*t20*t34*t35*t6
     9         3*t9-3.989422804014318d-1*t15*t20*t22*t63*t9-1.9947114020
     :         07159d-1*t10*t35*t41*t52*t9+1.994711402007159d-1*t15*t20*
     ;         t35*t52*t9+1.5957691216057274d+0*t12*t38*t41*t9-7.9788456
     <         08028636d-1*t15*t22*t34*t41*t9+7.978845608028636d-1*t20*t
     =         22*t34*t38*t9+1.329807601338106d-1*t10*t105*t20*t22*t9+3.
     >         5449077018110403d+0*t18*t25*t36*t69+3.5449077018110403d+0
     ?         *t18*t25*t36*t68+7.089815403622081d+0*t18*t25*t40*t64)-1.
     @         3298076013381063d-1*t10*t23*t53*t8*t9*t98-5.3192304053524
     1         25d-1*t10*t22*t34*t71*t8*t9+1.063846081070485d+0*t12*t15*
     2         t71*t8*t9+2.12769216214097d+0*t12*t23*t66*t8*t9-5.3192304
     3         05352425d-1*t10*t22*t42*t63*t8*t9+2.6596152026762127d-1*t
     4         10*t23*t34*t35*t63*t8*t9+5.319230405352425d-1*t15*t22*t23
     5         *t63*t8*t9+2.6596152026762127d-1*t10*t35*t42*t52*t8*t9-2.
     6         6596152026762127d-1*t15*t23*t35*t52*t8*t9-2.1276921621409
     7         7d+0*t12*t38*t42*t8*t9+1.063846081070485d+0*t15*t22*t34*t
     8         42*t8*t9-1.063846081070485d+0*t22*t23*t34*t38*t8*t9-1.773
     9         0768017841417d-1*t10*t105*t22*t23*t8*t9)*wght+3*t43*t73*w
     :         ght+3*t44*t72*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t7*(-3.5461536035682834d-1*t10*t12*
     1         t8*t9*(2.659615202676212d-1*t10*t22*t34*t79*t9-5.31923040
     2         5352424d-1*t12*t15*t79*t9+2.659615202676212d-1*t10*t12*(-
     3         7.073553026306428d-2*t11*t13*t15*t18*t40*t78+2*t18*t40*t7
     4         8+1.4147106052612857d-1*t11*t13*t18*t38*t78-7.07355302630
     5         6428d-2*t13*t15*t18*t34*t78-7.073553026306428d-2*t13*t15*
     6         t18*t40*t76+1.4147106052612857d-1*t13*t19*t38*t76+1.41371
     7         6694115413d+1*t14*t16*t18*t39*t47*t69-5.0d-1*t17*t18*t47*
     8         t69-3.536776513153214d-2*t13*t15*t18*t47*t69+1.4137166941
     9         15413d+1*t14*t16*t18*t39*t47*t68-5.0d-1*t17*t18*t47*t68-3
     :         .536776513153214d-2*t13*t15*t18*t47*t68-2.122065907891928
     ;         5d-1*t13*t19*t47*t66-5.0d-1*t18*t39*t47*t63-1.0d+0*t18*t3
     <         4*t39*t40*t47+1.4147106052612857d-1*t13*t18*t38*t40*t47+2
     =         .0d+0*t10*t17*t18*t40*t47+2.0d+0*t10*t18*t34*t39*t47-3.0d
     >         +0*t15*t17*t18*t47-3.536776513153214d-2*t110*t13*t15*t19-
     ?         3.536776513153214d-2*t11*t111*t13*t15*t18+t111*t18)*t9+2.
     @         659615202676212d-1*t10*t22*t41*t76*t9-1.329807601338106d-
     1         1*t10*t20*t34*t35*t76*t9-2.659615202676212d-1*t15*t20*t22
     2         *t76*t9+1.329807601338106d-1*t10*t22*t47*t70*t9+1.3298076
     3         01338106d-1*t10*t22*t48*t63*t9-6.649038006690531d-2*t10*t
     4         20*t35*t47*t63*t9+9.973557010035795d-2*t10*t20*t47*t52*t5
     5         3*t9-6.649038006690531d-2*t10*t35*t48*t52*t9+5.3192304053
     6         52424d-1*t12*t38*t48*t9-2.659615202676212d-1*t15*t22*t34*
     7         t48*t9-1.329807601338106d-1*t10*t34*t35*t41*t47*t9-2.6596
     8         15202676212d-1*t15*t22*t41*t47*t9+2.659615202676212d-1*t2
     9         0*t22*t38*t47*t9+1.329807601338106d-1*t15*t20*t34*t35*t47
     :         *t9+1.329807601338106d-1*t10*t110*t20*t22*t9+3.5449077018
     ;         110403d+0*t18*t25*(-7.049892022399704d+0*rhob*t21*t47*t52
     <         *t85+5.639913617919763d+0*rhob*t21*t34*t53*t76-3.75994241
     =         1946509d+0*t21*t35*t76+2.8199568089598814d+0*rhob*t21*t47
     >         *t53*t63+5.639913617919763d+0*t21*t34*t47*t53-1.879971205
     ?         9732546d+0*rhob*t110*t21*t35)+3.5449077018110403d+0*t18*t
     @         25*t36*t78+3.5449077018110403d+0*t18*t25*t40*t77+5.011495
     1         19714857d+1*t14*t16*t18*t25*t39*t47*t64+5.01149519714857d
     2         +1*t14*t16*t18*t25*t36*t39*t40*t47)-3.5461536035682834d-1
     3         *t10*t22*t34*t8*t80*t9+7.092307207136567d-1*t12*t15*t8*t8
     4         0*t9-3.5461536035682834d-1*t10*t22*t42*t76*t8*t9+1.773076
     5         8017841417d-1*t10*t23*t34*t35*t76*t8*t9+3.546153603568283
     6         4d-1*t15*t22*t23*t76*t8*t9-1.7730768017841417d-1*t10*t22*
     7         t47*t71*t8*t9-1.7730768017841417d-1*t10*t22*t49*t63*t8*t9
     8         +8.865384008920708d-2*t10*t23*t35*t47*t63*t8*t9-1.3298076
     9         013381063d-1*t10*t23*t47*t52*t53*t8*t9+8.865384008920708d
     :         -2*t10*t35*t49*t52*t8*t9-7.092307207136567d-1*t12*t38*t49
     ;         *t8*t9+3.5461536035682834d-1*t15*t22*t34*t49*t8*t9+1.7730
     <         768017841417d-1*t10*t34*t35*t42*t47*t8*t9+3.5461536035682
     =         834d-1*t15*t22*t42*t47*t8*t9-3.5461536035682834d-1*t22*t2
     >         3*t38*t47*t8*t9-1.7730768017841417d-1*t15*t23*t34*t35*t47
     ?         *t8*t9-1.7730768017841417d-1*t10*t110*t22*t23*t8*t9)*wght
     @         +2*t43*t82*wght+2*t44*t81*wght+t24*(-1.6799999999999998d-
     1         2*gammabb*t3*t31*t56*t75-1.1199999999999999d-2*gammabb*t2
     2         6*t32*t75-8.4d-3*gammabb*t3*t46*t56*t61+4.2d-3*t3*t32*t61
     3         -1.3066666666666665d-2*t54*t6-8.4d-3*t3*t55*t56+2.2399999
     4         999999998d-2*gammabb*t26*t31*t46*t56+2.52d-2*gammabb*t102
     5         *t3*t46*t55+1.3066666666666665d-2*gammabb*t32*t46*t54-1.1
     6         199999999999999d-2*t26*t31*t32+4.2d-3*gammabb*t109*t3*t32
     7         )*wght+t50*t73*wght+t51*t72*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t43*t97*wght+t44*t96*wght+t7*(-1.7
     1         730768017841417d-1*t10*t22*t34*t8*t9*t95+3.54615360356828
     2         34d-1*t12*t15*t8*t9*t95-3.5461536035682834d-1*t10*t12*t8*
     3         t9*(1.329807601338106d-1*t10*t22*t34*t9*t94-2.65961520267
     4         6212d-1*t12*t15*t9*t94+2.659615202676212d-1*t10*t12*t9*(1
     5         .9985948912206125d+2*t18*t40*t65*t84*t92*t93-7.0685834705
     6         77066d+0*t14*t16*t18*t34*t84*t93+1.413716694115413d+1*t14
     7         *t16*t18*t39*t40*t91-5.0d-1*t17*t18*t40*t91-3.53677651315
     8         3214d-2*t13*t15*t18*t40*t91-5.0d-1*t18*t34*t39*t91+7.0735
     9         53026306428d-2*t13*t19*t38*t91+1.0d+0*t10*t17*t18*t91-3.5
     :         34291735288533d+1*t14*t16*t18*t40*t67*t84+1.0d+0*t18*t34*
     ;         t67*t84+1.413716694115413d+1*rhob*t16*t18*t67*t84+2.82743
     <         3388230826d+1*t14*t16*t18*t39*t47*t78-1.0d+0*t17*t18*t47*
     =         t78-7.073553026306428d-2*t13*t15*t18*t47*t78-1.0d+0*t18*t
     >         39*t47*t76-3.536776513153214d-2*t113*t13*t15*t19-3.536776
     ?         513153214d-2*t11*t114*t13*t15*t18+t114*t18)+7.08483442268
     @         8148d+2*t18*t25*t36*t65*t84*t92*t93+3.5449077018110403d+0
     1         *t18*t25*(2.8199568089598814d+0*rhob*t21*t34*t53*t91-1.87
     2         99712059732546d+0*t21*t35*t91-7.049892022399704d+0*rhob*t
     3         21*t34*t84*t85+2.8199568089598814d+0*t21*t53*t84+5.639913
     4         617919763d+0*rhob*t21*t47*t53*t76-1.8799712059732546d+0*r
     5         hob*t113*t21*t35)+1.329807601338106d-1*t10*t22*t41*t9*t91
     6         -6.649038006690531d-2*t10*t20*t34*t35*t9*t91-1.3298076013
     7         38106d-1*t15*t20*t22*t9*t91+5.01149519714857d+1*t14*t16*t
     8         18*t25*t36*t39*t91+9.973557010035795d-2*t10*t20*t34*t53*t
     9         84*t9-6.649038006690531d-2*t10*t35*t41*t84*t9+6.649038006
     :         690531d-2*t15*t20*t35*t84*t9+2.659615202676212d-1*t10*t22
     ;         *t47*t79*t9+2.659615202676212d-1*t10*t22*t48*t76*t9-1.329
     <         807601338106d-1*t10*t20*t35*t47*t76*t9-1.329807601338106d
     =         -1*t10*t34*t35*t47*t48*t9-2.659615202676212d-1*t15*t22*t4
     >         7*t48*t9+1.329807601338106d-1*t10*t113*t20*t22*t9-1.00229
     ?         9039429714d+2*t14*t16*t18*t25*t36*t67*t84+1.0022990394297
     @         14d+2*t14*t16*t18*t25*t39*t47*t77)-1.7730768017841417d-1*
     1         t10*t22*t42*t8*t9*t91+8.865384008920708d-2*t10*t23*t34*t3
     2         5*t8*t9*t91+1.7730768017841417d-1*t15*t22*t23*t8*t9*t91-1
     3         .3298076013381063d-1*t10*t23*t34*t53*t8*t84*t9+8.86538400
     4         8920708d-2*t10*t35*t42*t8*t84*t9-8.865384008920708d-2*t15
     5         *t23*t35*t8*t84*t9-3.5461536035682834d-1*t10*t22*t47*t8*t
     6         80*t9-3.5461536035682834d-1*t10*t22*t49*t76*t8*t9+1.77307
     7         68017841417d-1*t10*t23*t35*t47*t76*t8*t9+1.77307680178414
     8         17d-1*t10*t34*t35*t47*t49*t8*t9+3.5461536035682834d-1*t15
     9         *t22*t47*t49*t8*t9-1.7730768017841417d-1*t10*t113*t22*t23
     :         *t8*t9)*wght+t24*(-8.4d-3*gammabb*t3*t31*t56*t90-5.599999
     ;         999999999d-3*gammabb*t26*t32*t90+1.1199999999999999d-2*ga
     <         mmabb*t26*t56*t86+2.52d-2*gammabb*t102*t3*t31*t86-1.67999
     =         99999999998d-2*gammabb*t3*t46*t56*t75+8.4d-3*t3*t32*t75-1
     >         .6799999999999998d-2*t3*t31*t46*t56-1.1199999999999999d-2
     ?         *t26*t32*t46+4.2d-3*gammabb*t112*t3*t32)*wght+2*t50*t82*w
     @         ght+2*t51*t81*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 3*t50*t97*wght+3*t51*t96*wght+t7*
     1         (-5.319230405352425d-1*t10*t22*t47*t8*t9*t95-3.5461536035
     2         682834d-1*t10*t12*t8*t9*(3.989422804014318d-1*t10*t22*t47
     3         *t9*t94+2.659615202676212d-1*t10*t12*t9*(5.99578467366183
     4         8d+2*t18*t47*t65*t91*t92*t93+1.0602875205865599d+2*t115*t
     5         14*t16*t18*t93-1.2990866792933983d+3*t115*t18*t65*t92/t11
     6         **5-1.0602875205865599d+2*t14*t16*t18*t47*t67*t91+2.82544
     7         69624923584d+3*rhob**6*t115*t18/(t11**6*t9**6)+1.41371669
     8         4115413d+1*t118*t14*t16*t18*t39-3.536776513153214d-2*t118
     9         *t13*t15*t19-5.0d-1*t118*t17*t18)+3.989422804014318d-1*t1
     :         0*t22*t48*t9*t91-1.994711402007159d-1*t10*t20*t35*t47*t9*
     ;         t91-2.826440000853771d+2*t18*t25*t37*t47*t83*t85*t91+2.99
     <         8945983256911d+1*rhob*t18*t21*t25*t47*t53*t91-1.994711402
     =         007159d-1*t10*t35*t48*t84*t9+9.973557010035795d-2*t10*t11
     >         5*t20*t53*t9+1.329807601338106d-1*t10*t118*t20*t22*t9-1.3
     ?         319284713741864d+3*t106*t115*t18*t25/(t12**11*t9**5)-2.49
     @         91216527140922d+1*rhob*t115*t18*t21*t25*t85+4.71073333475
     1         6285d+2*t115*t18*t25*t37*t83/t12**9-6.66432440723758d+0*r
     2         hob*t118*t18*t21*t25*t35)-5.319230405352425d-1*t10*t22*t4
     3         9*t8*t9*t91+2.6596152026762127d-1*t10*t23*t35*t47*t8*t9*t
     4         91+2.6596152026762127d-1*t10*t35*t49*t8*t84*t9-1.32980760
     5         13381063d-1*t10*t115*t23*t53*t8*t9-1.7730768017841417d-1*
     6         t10*t118*t22*t23*t8*t9)*wght+t24*(-2.52d-2*gammabb*t3*t46
     7         *t56*t90+1.26d-2*t3*t32*t90-2.52d-2*t3*t56*t86+4.2d-3*gam
     8         mabb*t117*t3*t32+2.52d-2*gammabb*t102*t116*t3)*wght+Cmat3
     9         (iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
