C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_cams12.F 
C> The nwxcm_x_cams12 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_cams12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)-{\it t_1}+1.0\\\\ 
C>   {\it t_4} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_2}\,{
C>    \it t_6}\,{\it t_9}+0.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{
C>    \it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}
C>    \,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{
C>    \it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{21}}
C>    \,{\it t_{23}}\,{\it t_{24}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}
C>    \,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{
C>    \it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{30}}
C>    \,{\it t_{32}}\,{\it t_{33}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{
C>    \it t_{23}}\,{\it t_{24}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{21}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}
C>    \,\left({\it t_{29}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(
C>    -0.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}
C>    -0.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(
C>    -0.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}
C>    \,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{
C>    \it t_{15}}}}\right)+0.2659615202676212\,{\it t_{12}}\,{
C>    \it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}
C>    -{{0.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({
C>    \it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)
C>    +{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{
C>    \it t_{32}}\,{\it t_{33}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{30}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}
C>    \,\left({\it t_{38}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_cams12(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = -2.6666666666666666d+0*gammaaa*t7/rhoa**3.666666666666
     1         6664d+0-5.333333333333333d+0*t4*t5/rhoa**6.33333333333333
     2         3d+0
            t37 = 1/t9**2
            t38 = 9.305257363491002d-1*t13*t2*t3*t36*t37+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t39 = 1/t19**3
            t40 = 1/t18**2
            t41 = 1.413716694115413d+1*t21*t23*t38*t40-2.827433388230826
     1         d+1*rhoa*t23*t24
            t42 = 1/t2
            t43 = t7*t8+2*gammaaa*t4*t6
            t44 = 9.305257363491002d-1*t13*t2*t3*t37*t43+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t42
            fnc(iq) = 2*t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = (t14*(-3.5461536035682834d-1*t15*t16*t17*t1
     1         9*(2.659615202676212d-1*t16*t17*t19*(-3.536776513153214d-
     2         2*t18*t20*t22*t25*t41+t25*t41-3.536776513153214d-2*t20*t2
     3         2*t26*t38+7.073553026306428d-2*t18*t20*t26/rhoa**3)+3.544
     4         9077018110403d+0*t25*t32*(3.759942411946509d+0*t28*t29-1.
     5         8799712059732546d+0*rhoa*t28*t38*t39)+1.329807601338106d-
     6         1*t16*t17*t27*t29*t38-2.659615202676212d-1*t16*t19*t22*t2
     7         7)-1.7730768017841417d-1*t15*t16*t17*t29*t30*t38+3.546153
     8         6035682834d-1*t15*t16*t19*t22*t30)+t31*(-9.30525736349100
     9         2d-1*t13*t2*t3*t36*t37-1.2407009817988002d+0*t10*t13*t3*t
     :         35-1.2407009817988002d+0*t1*t35+2.4814019635976003d+0*gam
     ;         maaa*t10*t11*t3*t33*t34))*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t14*(-3.5461536035682834d-1*t15*t16*t17*t
     1         19*(2.659615202676212d-1*t16*t17*t19*(1.413716694115413d+
     2         1*t21*t23*t25*t40*t44-3.536776513153214d-2*t20*t22*t26*t4
     3         4-5.0d-1*t24*t25*t44)-6.66432440723758d+0*rhoa*t25*t28*t3
     4         2*t39*t44+1.329807601338106d-1*t16*t17*t27*t29*t44)-1.773
     5         0768017841417d-1*t15*t16*t17*t29*t30*t44)+t31*(-9.3052573
     6         63491002d-1*t13*t2*t3*t37*t43-9.305257363491002d-1*t10*t1
     7         1*t3*t33*t42))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = -param(6)
            t16 = param(7)
            t17 = param(8)
            t18 = 1/rhoa
            t19 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t20 = t19**5.0d-1
            t21 = t17**2
            t22 = rhoa**2
            t23 = 1/t22
            t24 = 1/t21
            t25 = 1/t19
            t26 = exp(-1.413716694115413d+1*t22*t24*t25)
            t27 = t26-1.0d+0
            t28 = -3.536776513153214d-2*t19*t21*t23*t27+t26-1.5d+0
            t29 = 1/t17
            t30 = 1/t20
            t31 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         9*t30)+2.659615202676212d-1*t17*t18*t20*t28
            t32 = -3.5461536035682834d-1*t16*t17*t18*t20*t31+t15+1.0d+0
            t33 = rhob**1.3333333333333333d+0
            t34 = gammabb**2
            t35 = 1/rhob**5.333333333333333d+0
            t36 = 1/rhob**2.6666666666666666d+0
            t37 = gammabb*t36*t7+t34*t35*t4+1.0d+0
            t38 = 1.0d+0-1.0d+0/t37
            t39 = gammabb*t11*t36+1.0d+0
            t40 = 1.0d+0-1.0d+0/t39
            t41 = -9.305257363491002d-1*t3*t33*t38*t40-9.305257363491002
     1         d-1*t1*t33
            t42 = 1/rhob
            t43 = 9.305257363491002d-1*t3*t33*t38*t40+9.305257363491002d
     1         -1*t1*t33
            t44 = t43**5.0d-1
            t45 = rhob**2
            t46 = 1/t45
            t47 = 1/t43
            t48 = exp(-1.413716694115413d+1*t24*t45*t47)
            t49 = t48-1.0d+0
            t50 = -3.536776513153214d-2*t21*t43*t46*t49+t48-1.5d+0
            t51 = 1/t44
            t52 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         9*t51)+2.659615202676212d-1*t17*t42*t44*t50
            t53 = -3.5461536035682834d-1*t16*t17*t42*t44*t52+t15+1.0d+0
            t54 = 5.641895835477563d-1
            t55 = 1/t12**2
            t56 = 1/rhoa**2.3333333333333334d+0
            t57 = rhoa**3.333333333333333d-1
            t58 = -2.6666666666666666d+0*gammaaa*t7/rhoa**3.666666666666
     1         6664d+0-5.333333333333333d+0*t4*t5/rhoa**6.33333333333333
     2         3d+0
            t59 = 1/t9**2
            t60 = 9.305257363491002d-1*t13*t2*t3*t58*t59+1.2407009817988
     1         002d+0*t10*t13*t3*t57+1.2407009817988002d+0*t1*t57-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t55*t56
            t61 = 1/t20**3
            t62 = 1/t19**2
            t63 = 1.413716694115413d+1*t22*t24*t60*t62-2.827433388230826
     1         d+1*rhoa*t24*t25
            t64 = 1/t39**2
            t65 = 1/rhob**2.3333333333333334d+0
            t66 = rhob**3.333333333333333d-1
            t67 = -2.6666666666666666d+0*gammabb*t7/rhob**3.666666666666
     1         6664d+0-5.333333333333333d+0*t34*t4/rhob**6.3333333333333
     2         33d+0
            t68 = 1/t37**2
            t69 = 9.305257363491002d-1*t3*t33*t40*t67*t68+1.240700981798
     1         8002d+0*t3*t38*t40*t66+1.2407009817988002d+0*t1*t66-2.481
     2         4019635976003d+0*gammabb*t11*t3*t38*t64*t65
            t70 = 1/t44**3
            t71 = 1/t43**2
            t72 = 1.413716694115413d+1*t24*t45*t69*t71-2.827433388230826
     1         d+1*rhob*t24*t47
            t73 = 1/t2
            t74 = t7*t8+2*gammaaa*t4*t6
            t75 = 9.305257363491002d-1*t13*t2*t3*t59*t74+9.3052573634910
     1         02d-1*t10*t11*t3*t55*t73
            t76 = 1/t33
            t77 = t36*t7+2*gammabb*t35*t4
            t78 = 9.305257363491002d-1*t3*t33*t40*t68*t77+9.305257363491
     1         002d-1*t11*t3*t38*t64*t76
            fnc(iq) = (t41*t53+t14*t32)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t14*(-3.5461536035682834d-1*t16*t17*t18*t2
     1         0*(2.659615202676212d-1*t17*t18*t20*(-3.536776513153214d-
     2         2*t19*t21*t23*t26*t63+t26*t63-3.536776513153214d-2*t21*t2
     3         3*t27*t60+7.073553026306428d-2*t19*t21*t27/rhoa**3)+3.544
     4         9077018110403d+0*t26*t54*(3.759942411946509d+0*t29*t30-1.
     5         8799712059732546d+0*rhoa*t29*t60*t61)+1.329807601338106d-
     6         1*t17*t18*t28*t30*t60-2.659615202676212d-1*t17*t20*t23*t2
     7         8)-1.7730768017841417d-1*t16*t17*t18*t30*t31*t60+3.546153
     8         6035682834d-1*t16*t17*t20*t23*t31)+t32*(-9.30525736349100
     9         2d-1*t13*t2*t3*t58*t59-1.2407009817988002d+0*t10*t13*t3*t
     :         57-1.2407009817988002d+0*t1*t57+2.4814019635976003d+0*gam
     ;         maaa*t10*t11*t3*t55*t56))*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t41*(-3.5461536035682834d-1*t16*t17*t42*t4
     1         4*(2.659615202676212d-1*t17*t42*t44*(-3.536776513153214d-
     2         2*t21*t43*t46*t48*t72+t48*t72-3.536776513153214d-2*t21*t4
     3         6*t49*t69+7.073553026306428d-2*t21*t43*t49/rhob**3)+3.544
     4         9077018110403d+0*t48*t54*(3.759942411946509d+0*t29*t51-1.
     5         8799712059732546d+0*rhob*t29*t69*t70)+1.329807601338106d-
     6         1*t17*t42*t50*t51*t69-2.659615202676212d-1*t17*t44*t46*t5
     7         0)-1.7730768017841417d-1*t16*t17*t42*t51*t52*t69+3.546153
     8         6035682834d-1*t16*t17*t44*t46*t52)+t53*(-9.30525736349100
     9         2d-1*t3*t33*t40*t67*t68-1.2407009817988002d+0*t3*t38*t40*
     :         t66-1.2407009817988002d+0*t1*t66+2.4814019635976003d+0*ga
     ;         mmabb*t11*t3*t38*t64*t65))*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t14*(-3.5461536035682834d-1*t16*t17*t18*t
     1         20*(2.659615202676212d-1*t17*t18*t20*(1.413716694115413d+
     2         1*t22*t24*t26*t62*t75-3.536776513153214d-2*t21*t23*t27*t7
     3         5-5.0d-1*t25*t26*t75)-6.66432440723758d+0*rhoa*t26*t29*t5
     4         4*t61*t75+1.329807601338106d-1*t17*t18*t28*t30*t75)-1.773
     5         0768017841417d-1*t16*t17*t18*t30*t31*t75)+t32*(-9.3052573
     6         63491002d-1*t13*t2*t3*t59*t74-9.305257363491002d-1*t10*t1
     7         1*t3*t55*t73))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t41*(-3.5461536035682834d-1*t16*t17*t42*t
     1         44*(2.659615202676212d-1*t17*t42*t44*(1.413716694115413d+
     2         1*t24*t45*t48*t71*t78-3.536776513153214d-2*t21*t46*t49*t7
     3         8-5.0d-1*t47*t48*t78)-6.66432440723758d+0*rhob*t29*t48*t5
     4         4*t70*t78+1.329807601338106d-1*t17*t42*t50*t51*t78)-1.773
     5         0768017841417d-1*t16*t17*t42*t51*t52*t78)+t53*(-9.3052573
     6         63491002d-1*t3*t33*t40*t68*t77-9.305257363491002d-1*t11*t
     7         3*t38*t64*t76))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = -2.6666666666666666d+0*gammaaa*t7/rhoa**3.666666666666
     1         6664d+0-5.333333333333333d+0*t4*t5/rhoa**6.33333333333333
     2         3d+0
            t37 = 1/t9**2
            t38 = 9.305257363491002d-1*t13*t2*t3*t36*t37+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t39 = 1/t19**3
            t40 = 1/t18**2
            t41 = 1.413716694115413d+1*t21*t23*t38*t40-2.827433388230826
     1         d+1*rhoa*t23*t24
            t42 = 1/t2
            t43 = t7*t8+2*gammaaa*t4*t6
            t44 = 9.305257363491002d-1*t13*t2*t3*t37*t43+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t42
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = t14*(-3.5461536035682834d-1*t15*t16*t17*t19
     1         *(2.659615202676212d-1*t16*t17*t19*(-3.536776513153214d-2
     2         *t18*t20*t22*t25*t41+t25*t41-3.536776513153214d-2*t20*t22
     3         *t26*t38+7.073553026306428d-2*t18*t20*t26/rhoa**3)+3.5449
     4         077018110403d+0*t25*t32*(3.759942411946509d+0*t28*t29-1.8
     5         799712059732546d+0*rhoa*t28*t38*t39)+1.329807601338106d-1
     6         *t16*t17*t27*t29*t38-2.659615202676212d-1*t16*t19*t22*t27
     7         )-1.7730768017841417d-1*t15*t16*t17*t29*t30*t38+3.5461536
     8         035682834d-1*t15*t16*t19*t22*t30)*wght+t31*(-9.3052573634
     9         91002d-1*t13*t2*t3*t36*t37-1.2407009817988002d+0*t10*t13*
     :         t3*t35-1.2407009817988002d+0*t1*t35+2.4814019635976003d+0
     ;         *gammaaa*t10*t11*t3*t33*t34)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t14*(-3.5461536035682834d-1*t15*t16*t17*t1
     1         9*(2.659615202676212d-1*t16*t17*t19*(1.413716694115413d+1
     2         *t21*t23*t25*t40*t44-3.536776513153214d-2*t20*t22*t26*t44
     3         -5.0d-1*t24*t25*t44)-6.66432440723758d+0*rhoa*t25*t28*t32
     4         *t39*t44+1.329807601338106d-1*t16*t17*t27*t29*t44)-1.7730
     5         768017841417d-1*t15*t16*t17*t29*t30*t44)*wght+t31*(-9.305
     6         257363491002d-1*t13*t2*t3*t37*t43-9.305257363491002d-1*t1
     7         0*t11*t3*t33*t42)*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhob**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammabb**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhob**2.6666666666666666d+0
            t9 = gammabb*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammabb*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhob
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhob**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhob**2.3333333333333334d+0
            t35 = rhob**3.333333333333333d-1
            t36 = -2.6666666666666666d+0*gammabb*t7/rhob**3.666666666666
     1         6664d+0-5.333333333333333d+0*t4*t5/rhob**6.33333333333333
     2         3d+0
            t37 = 1/t9**2
            t38 = 9.305257363491002d-1*t13*t2*t3*t36*t37+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammabb*t10*t11*t3*t33*t34
            t39 = 1/t19**3
            t40 = 1/t18**2
            t41 = 1.413716694115413d+1*t21*t23*t38*t40-2.827433388230826
     1         d+1*rhob*t23*t24
            t42 = 1/t2
            t43 = t7*t8+2*gammabb*t4*t6
            t44 = 9.305257363491002d-1*t13*t2*t3*t37*t43+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t42
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RB) = t14*(-3.5461536035682834d-1*t15*t16*t17*t19
     1         *(2.659615202676212d-1*t16*t17*t19*(-3.536776513153214d-2
     2         *t18*t20*t22*t25*t41+t25*t41-3.536776513153214d-2*t20*t22
     3         *t26*t38+7.073553026306428d-2*t18*t20*t26/rhob**3)+3.5449
     4         077018110403d+0*t25*t32*(3.759942411946509d+0*t28*t29-1.8
     5         799712059732546d+0*rhob*t28*t38*t39)+1.329807601338106d-1
     6         *t16*t17*t27*t29*t38-2.659615202676212d-1*t16*t19*t22*t27
     7         )-1.7730768017841417d-1*t15*t16*t17*t29*t30*t38+3.5461536
     8         035682834d-1*t15*t16*t19*t22*t30)*wght+t31*(-9.3052573634
     9         91002d-1*t13*t2*t3*t36*t37-1.2407009817988002d+0*t10*t13*
     :         t3*t35-1.2407009817988002d+0*t1*t35+2.4814019635976003d+0
     ;         *gammabb*t10*t11*t3*t33*t34)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t14*(-3.5461536035682834d-1*t15*t16*t17*t1
     1         9*(2.659615202676212d-1*t16*t17*t19*(1.413716694115413d+1
     2         *t21*t23*t25*t40*t44-3.536776513153214d-2*t20*t22*t26*t44
     3         -5.0d-1*t24*t25*t44)-6.66432440723758d+0*rhob*t25*t28*t32
     4         *t39*t44+1.329807601338106d-1*t16*t17*t27*t29*t44)-1.7730
     5         768017841417d-1*t15*t16*t17*t29*t30*t44)*wght+t31*(-9.305
     6         257363491002d-1*t13*t2*t3*t37*t43-9.305257363491002d-1*t1
     7         0*t11*t3*t33*t42)*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_cams12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)-{\it t_1}+1.0\\\\ 
C>   {\it t_4} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_2}\,{
C>    \it t_6}\,{\it t_9}+0.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{
C>    \it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}
C>    \,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{
C>    \it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{21}}
C>    \,{\it t_{23}}\,{\it t_{24}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}
C>    \,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{
C>    \it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{30}}
C>    \,{\it t_{32}}\,{\it t_{33}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{
C>    \it t_{23}}\,{\it t_{24}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{21}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}
C>    \,\left({\it t_{29}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(
C>    -0.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}
C>    -0.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(
C>    -0.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}
C>    \,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{
C>    \it t_{15}}}}\right)+0.2659615202676212\,{\it t_{12}}\,{
C>    \it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}
C>    -{{0.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({
C>    \it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)
C>    +{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{
C>    \it t_{32}}\,{\it t_{33}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{30}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}
C>    \,\left({\it t_{38}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_cams12_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = 1/rhoa**6.333333333333333d+0
            t37 = 1/rhoa**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammaaa*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         a*t28*t40*t41
            t43 = rhoa**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhoa*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammaaa*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhoa*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = 1/rhoa**6
            t63 = 1/rhoa**3.3333333333333337d+0
            t64 = 1/rhoa**6.666666666666666d-1
            t65 = t38**2
            t66 = 1/t9**3
            t67 = 9.777777777777779d+0*gammaaa*t7/rhoa**4.66666666666666
     1         7d+0+3.377777777777778d+1*t4*t5/rhoa**7.333333333333333d+
     2         0
            t68 = 9.305257363491002d-1*t13*t2*t3*t39*t67-1.8610514726982
     1         003d+0*t13*t2*t3*t65*t66+4.135669939329334d-1*t10*t13*t3*
     2         t64+4.135669939329334d-1*t1*t64+2.481401963597601d+0*gamm
     3         aaa*t10*t11*t3*t33*t63-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t62+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t69 = rhoa**4
            t70 = 1/t69
            t71 = 1/t18**3
            t72 = -2.827433388230826d+1*t21*t23*t58*t71+1.41371669411541
     1         3d+1*t21*t23*t45*t68+5.654866776461652d+1*rhoa*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t73 = t46**2
            t74 = 1/rhoa**5
            t75 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmaaa*t36*t4
            t76 = 9.305257363491002d-1*t13*t2*t3*t39*t75+4.9628039271952
     1         007d+0*gammaaa*t10*t3*t60*t61*t74-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t66+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t77 = 1.413716694115413d+1*t21*t23*t45*t76-2.827433388230826
     1         d+1*t21*t23*t40*t53*t71+2.827433388230826d+1*rhoa*t23*t45
     2         *t53
            t78 = t53**2
            t79 = t52**2
            t80 = -1.8610514726982003d+0*t13*t2*t3*t66*t79-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t70+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t70+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            fnc(iq) = 2*t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = (t31*t50+t14*t49)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t31*t57+t14*t56)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t14*(-3.5461536035682834d-1*t15*t16*t1
     1         7*t19*(2.659615202676212d-1*t16*t17*t19*(-3.5367765131532
     2         14d-2*t18*t20*t22*t25*t73+t25*t73-3.536776513153214d-2*t1
     3         8*t20*t22*t25*t72+t25*t72-2.1220659078919285d-1*t18*t20*t
     4         26*t70-3.536776513153214d-2*t20*t22*t26*t68+1.41471060526
     5         12857d-1*t18*t20*t25*t44*t46-7.073553026306428d-2*t20*t22
     6         *t25*t40*t46+1.4147106052612857d-1*t20*t26*t40*t44)+3.544
     7         9077018110403d+0*t25*t32*(-1.8799712059732546d+0*rhoa*t28
     8         *t41*t68+2.8199568089598814d+0*rhoa*t28*t58*t59-3.7599424
     9         11946509d+0*t28*t40*t41)+1.329807601338106d-1*t16*t17*t27
     :         *t29*t68-6.649038006690531d-2*t16*t17*t27*t41*t58+2.65961
     ;         5202676212d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t1
     <         6*t19*t22*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.319
     =         230405352424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*
     >         t22*t27*t29*t40)-1.7730768017841417d-1*t15*t16*t17*t29*t3
     ?         0*t68+8.865384008920708d-2*t15*t16*t17*t30*t41*t58-3.5461
     @         536035682834d-1*t15*t16*t17*t29*t40*t48+7.092307207136567
     1         d-1*t15*t16*t19*t22*t48-7.092307207136567d-1*t15*t16*t19*
     2         t30*t44+3.5461536035682834d-1*t15*t16*t22*t29*t30*t40)+t3
     3         1*(-9.305257363491002d-1*t13*t2*t3*t39*t67+1.861051472698
     4         2003d+0*t13*t2*t3*t65*t66-4.135669939329334d-1*t10*t13*t3
     5         *t64-4.135669939329334d-1*t1*t64-2.481401963597601d+0*gam
     6         maaa*t10*t11*t3*t33*t63+1.3234143805853868d+1*t10*t3*t5*t
     7         60*t61*t62-2.4814019635976003d+0*t13*t3*t35*t38*t39+4.962
     8         8039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39)+2*t49*t5
     9         0)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t14*(-3.5461536035682834d-1*t15*t16*t
     1         17*t19*(2.659615202676212d-1*t16*t17*t19*(-3.536776513153
     2         214d-2*t18*t20*t22*t25*t77+t25*t77-3.536776513153214d-2*t
     3         20*t22*t26*t76+1.413716694115413d+1*t21*t23*t25*t45*t46*t
     4         53-5.0d-1*t24*t25*t46*t53-3.536776513153214d-2*t20*t22*t2
     5         5*t46*t53-5.0d-1*t25*t40*t45*t53+7.073553026306428d-2*t20
     6         *t26*t44*t53+1.0d+0*t17*t24*t25*t53)+3.5449077018110403d+
     7         0*t25*t32*(-1.8799712059732546d+0*rhoa*t28*t41*t76+2.8199
     8         568089598814d+0*rhoa*t28*t40*t53*t59-1.8799712059732546d+
     9         0*t28*t41*t53)+1.329807601338106d-1*t16*t17*t27*t29*t76+1
     :         .329807601338106d-1*t16*t17*t29*t40*t54-2.659615202676212
     ;         d-1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t47*
     <         t53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.649
     =         038006690531d-2*t16*t17*t27*t40*t41*t53-1.329807601338106
     >         d-1*t16*t22*t27*t29*t53)-1.7730768017841417d-1*t15*t16*t1
     ?         7*t29*t30*t76-1.7730768017841417d-1*t15*t16*t17*t29*t40*t
     @         55+3.5461536035682834d-1*t15*t16*t19*t22*t55-1.7730768017
     1         841417d-1*t15*t16*t17*t29*t48*t53+8.865384008920708d-2*t1
     2         5*t16*t17*t30*t40*t41*t53+1.7730768017841417d-1*t15*t16*t
     3         22*t29*t30*t53)+t31*(-9.305257363491002d-1*t13*t2*t3*t39*
     4         t75-4.9628039271952007d+0*gammaaa*t10*t3*t60*t61*t74+1.86
     5         10514726982003d+0*t13*t2*t3*t38*t52*t66-1.240700981798800
     6         2d+0*t13*t3*t35*t39*t52+2.4814019635976003d+0*gammaaa*t11
     7         *t3*t33*t34*t39*t52-9.305257363491002d-1*t11*t3*t33*t38*t
     8         39*t51+1.2407009817988002d+0*t10*t11*t3*t33*t34)+t49*t57+
     9         t50*t56)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t14*(-3.5461536035682834d-1*t15*t16*
     1         t17*t19*(2.659615202676212d-1*t16*t17*t19*(1.413716694115
     2         413d+1*t21*t23*t25*t45*t80-3.536776513153214d-2*t20*t22*t
     3         26*t80-5.0d-1*t24*t25*t80-3.534291735288533d+1*t21*t23*t2
     4         5*t71*t78+1.9985948912206125d+2*t25*t69*t78/(t16**4*t18**
     5         4))-6.66432440723758d+0*rhoa*t25*t28*t32*t41*t80+1.329807
     6         601338106d-1*t16*t17*t27*t29*t80+9.996486610856369d+0*rho
     7         a*t25*t28*t32*t59*t78-9.421466669512572d+1*t25*t32*t43*t7
     8         8/(t16**3*t19**7)-6.649038006690531d-2*t16*t17*t27*t41*t7
     9         8+2.659615202676212d-1*t16*t17*t29*t53*t54)-1.77307680178
     :         41417d-1*t15*t16*t17*t29*t30*t80+8.865384008920708d-2*t15
     ;         *t16*t17*t30*t41*t78-3.5461536035682834d-1*t15*t16*t17*t2
     <         9*t53*t55)+t31*(1.8610514726982003d+0*t13*t2*t3*t66*t79+1
     =         .8610514726982003d+0*t10*t3*t60*t61*t70-1.861051472698200
     >         3d+0*t13*t3*t39*t4*t70-1.8610514726982003d+0*t11*t3*t33*t
     ?         39*t51*t52)+2*t56*t57)*wght+Cmat2(iq,D2_GAA_GAA)
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
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = -param(6)
            t16 = param(7)
            t17 = param(8)
            t18 = 1/rhoa
            t19 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t20 = t19**5.0d-1
            t21 = t17**2
            t22 = rhoa**2
            t23 = 1/t22
            t24 = 1/t21
            t25 = 1/t19
            t26 = exp(-1.413716694115413d+1*t22*t24*t25)
            t27 = t26-1.0d+0
            t28 = -3.536776513153214d-2*t19*t21*t23*t27+t26-1.5d+0
            t29 = 1/t17
            t30 = 1/t20
            t31 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         9*t30)+2.659615202676212d-1*t17*t18*t20*t28
            t32 = -3.5461536035682834d-1*t16*t17*t18*t20*t31+t15+1.0d+0
            t33 = rhob**1.3333333333333333d+0
            t34 = gammabb**2
            t35 = 1/rhob**5.333333333333333d+0
            t36 = 1/rhob**2.6666666666666666d+0
            t37 = gammabb*t36*t7+t34*t35*t4+1.0d+0
            t38 = 1.0d+0-1.0d+0/t37
            t39 = gammabb*t11*t36+1.0d+0
            t40 = 1.0d+0-1.0d+0/t39
            t41 = -9.305257363491002d-1*t3*t33*t38*t40-9.305257363491002
     1         d-1*t1*t33
            t42 = 1/rhob
            t43 = 9.305257363491002d-1*t3*t33*t38*t40+9.305257363491002d
     1         -1*t1*t33
            t44 = t43**5.0d-1
            t45 = rhob**2
            t46 = 1/t45
            t47 = 1/t43
            t48 = exp(-1.413716694115413d+1*t24*t45*t47)
            t49 = t48-1.0d+0
            t50 = -3.536776513153214d-2*t21*t43*t46*t49+t48-1.5d+0
            t51 = 1/t44
            t52 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         9*t51)+2.659615202676212d-1*t17*t42*t44*t50
            t53 = -3.5461536035682834d-1*t16*t17*t42*t44*t52+t15+1.0d+0
            t54 = 5.641895835477563d-1
            t55 = 1/t12**2
            t56 = 1/rhoa**2.3333333333333334d+0
            t57 = rhoa**3.333333333333333d-1
            t58 = 1/rhoa**6.333333333333333d+0
            t59 = 1/rhoa**3.6666666666666664d+0
            t60 = -2.6666666666666666d+0*gammaaa*t59*t7-5.33333333333333
     1         3d+0*t4*t5*t58
            t61 = 1/t9**2
            t62 = 9.305257363491002d-1*t13*t2*t3*t60*t61+1.2407009817988
     1         002d+0*t10*t13*t3*t57+1.2407009817988002d+0*t1*t57-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t55*t56
            t63 = 1/t20**3
            t64 = 3.759942411946509d+0*t29*t30-1.8799712059732546d+0*rho
     1         a*t29*t62*t63
            t65 = rhoa**3
            t66 = 1/t65
            t67 = 1/t19**2
            t68 = 1.413716694115413d+1*t22*t24*t62*t67-2.827433388230826
     1         d+1*rhoa*t24*t25
            t69 = -3.536776513153214d-2*t19*t21*t23*t26*t68+t26*t68+7.07
     1         3553026306428d-2*t19*t21*t27*t66-3.536776513153214d-2*t21
     2         *t23*t27*t62
            t70 = 2.659615202676212d-1*t17*t18*t20*t69+3.544907701811040
     1         3d+0*t26*t54*t64+1.329807601338106d-1*t17*t18*t28*t30*t62
     2         -2.659615202676212d-1*t17*t20*t23*t28
            t71 = -3.5461536035682834d-1*t16*t17*t18*t20*t70-1.773076801
     1         7841417d-1*t16*t17*t18*t30*t31*t62+3.5461536035682834d-1*
     2         t16*t17*t20*t23*t31
            t72 = -9.305257363491002d-1*t13*t2*t3*t60*t61-1.240700981798
     1         8002d+0*t10*t13*t3*t57-1.2407009817988002d+0*t1*t57+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t55*t56
            t73 = 1/t39**2
            t74 = 1/rhob**2.3333333333333334d+0
            t75 = rhob**3.333333333333333d-1
            t76 = 1/rhob**6.333333333333333d+0
            t77 = 1/rhob**3.6666666666666664d+0
            t78 = -2.6666666666666666d+0*gammabb*t7*t77-5.33333333333333
     1         3d+0*t34*t4*t76
            t79 = 1/t37**2
            t80 = 9.305257363491002d-1*t3*t33*t40*t78*t79+1.240700981798
     1         8002d+0*t3*t38*t40*t75+1.2407009817988002d+0*t1*t75-2.481
     2         4019635976003d+0*gammabb*t11*t3*t38*t73*t74
            t81 = 1/t44**3
            t82 = 3.759942411946509d+0*t29*t51-1.8799712059732546d+0*rho
     1         b*t29*t80*t81
            t83 = rhob**3
            t84 = 1/t83
            t85 = 1/t43**2
            t86 = 1.413716694115413d+1*t24*t45*t80*t85-2.827433388230826
     1         d+1*rhob*t24*t47
            t87 = -3.536776513153214d-2*t21*t43*t46*t48*t86+t48*t86+7.07
     1         3553026306428d-2*t21*t43*t49*t84-3.536776513153214d-2*t21
     2         *t46*t49*t80
            t88 = 2.659615202676212d-1*t17*t42*t44*t87+3.544907701811040
     1         3d+0*t48*t54*t82+1.329807601338106d-1*t17*t42*t50*t51*t80
     2         -2.659615202676212d-1*t17*t44*t46*t50
            t89 = -3.5461536035682834d-1*t16*t17*t42*t44*t88-1.773076801
     1         7841417d-1*t16*t17*t42*t51*t52*t80+3.5461536035682834d-1*
     2         t16*t17*t44*t46*t52
            t90 = -9.305257363491002d-1*t3*t33*t40*t78*t79-1.24070098179
     1         88002d+0*t3*t38*t40*t75-1.2407009817988002d+0*t1*t75+2.48
     2         14019635976003d+0*gammabb*t11*t3*t38*t73*t74
            t91 = 1/t2
            t92 = t7*t8+2*gammaaa*t4*t6
            t93 = 9.305257363491002d-1*t13*t2*t3*t61*t92+9.3052573634910
     1         02d-1*t10*t11*t3*t55*t91
            t94 = 1.413716694115413d+1*t22*t24*t26*t67*t93-3.53677651315
     1         3214d-2*t21*t23*t27*t93-5.0d-1*t25*t26*t93
            t95 = 2.659615202676212d-1*t17*t18*t20*t94-6.66432440723758d
     1         +0*rhoa*t26*t29*t54*t63*t93+1.329807601338106d-1*t17*t18*
     2         t28*t30*t93
            t96 = -3.5461536035682834d-1*t16*t17*t18*t20*t95-1.773076801
     1         7841417d-1*t16*t17*t18*t30*t31*t93
            t97 = -9.305257363491002d-1*t13*t2*t3*t61*t92-9.305257363491
     1         002d-1*t10*t11*t3*t55*t91
            t98 = 1/t33
            t99 = t36*t7+2*gammabb*t35*t4
            t100 = 9.305257363491002d-1*t3*t33*t40*t79*t99+9.30525736349
     1         1002d-1*t11*t3*t38*t73*t98
            t101 = 1.413716694115413d+1*t100*t24*t45*t48*t85-3.536776513
     1         153214d-2*t100*t21*t46*t49-5.0d-1*t100*t47*t48
            t102 = -6.66432440723758d+0*rhob*t100*t29*t48*t54*t81+1.3298
     1         07601338106d-1*t100*t17*t42*t50*t51+2.659615202676212d-1*
     2         t101*t17*t42*t44
            t103 = -1.7730768017841417d-1*t100*t16*t17*t42*t51*t52-3.546
     1         1536035682834d-1*t102*t16*t17*t42*t44
            t104 = -9.305257363491002d-1*t3*t33*t40*t79*t99-9.3052573634
     1         91002d-1*t11*t3*t38*t73*t98
            t105 = t62**2
            t106 = 1/t20**5
            t107 = t11**2
            t108 = 1/t12**3
            t109 = 1/rhoa**6
            t110 = 1/rhoa**3.3333333333333337d+0
            t111 = 1/rhoa**6.666666666666666d-1
            t112 = t60**2
            t113 = 1/t9**3
            t114 = 9.777777777777779d+0*gammaaa*t7/rhoa**4.6666666666666
     1         67d+0+3.377777777777778d+1*t4*t5/rhoa**7.333333333333333d
     2         +0
            t115 = 2.4814019635976003d+0*t13*t3*t57*t60*t61-4.9628039271
     1         952007d+0*gammaaa*t11*t3*t55*t56*t60*t61+9.30525736349100
     2         2d-1*t114*t13*t2*t3*t61+2.481401963597601d+0*gammaaa*t10*
     3         t11*t110*t3*t55-1.3234143805853868d+1*t10*t107*t108*t109*
     4         t3*t5-1.8610514726982003d+0*t112*t113*t13*t2*t3+4.1356699
     5         39329334d-1*t10*t111*t13*t3+4.135669939329334d-1*t1*t111
            t116 = rhoa**4
            t117 = 1/t116
            t118 = 1/t19**3
            t119 = 5.654866776461652d+1*rhoa*t24*t62*t67+1.4137166941154
     1         13d+1*t115*t22*t24*t67-2.827433388230826d+1*t24*t25-2.827
     2         433388230826d+1*t105*t118*t22*t24
            t120 = t68**2
            t121 = t80**2
            t122 = 1/t44**5
            t123 = 1/t39**3
            t124 = 1/rhob**6
            t125 = 1/rhob**3.3333333333333337d+0
            t126 = 1/rhob**6.666666666666666d-1
            t127 = t78**2
            t128 = 1/t37**3
            t129 = 9.777777777777779d+0*gammabb*t7/rhob**4.6666666666666
     1         67d+0+3.377777777777778d+1*t34*t4/rhob**7.333333333333333
     2         d+0
            t130 = 2.4814019635976003d+0*t3*t40*t75*t78*t79-4.9628039271
     1         952007d+0*gammabb*t11*t3*t73*t74*t78*t79+9.30525736349100
     2         2d-1*t129*t3*t33*t40*t79+2.481401963597601d+0*gammabb*t11
     3         *t125*t3*t38*t73+4.135669939329334d-1*t126*t3*t38*t40-1.8
     4         610514726982003d+0*t127*t128*t3*t33*t40-1.323414380585386
     5         8d+1*t107*t123*t124*t3*t34*t38+4.135669939329334d-1*t1*t1
     6         26
            t131 = rhob**4
            t132 = 1/t131
            t133 = 1/t43**3
            t134 = 5.654866776461652d+1*rhob*t24*t80*t85+1.4137166941154
     1         13d+1*t130*t24*t45*t85-2.827433388230826d+1*t24*t47-2.827
     2         433388230826d+1*t121*t133*t24*t45
            t135 = t86**2
            t136 = 1/rhoa**5
            t137 = -2.6666666666666666d+0*t59*t7-1.0666666666666666d+1*g
     1         ammaaa*t4*t58
            t138 = 1.2407009817988002d+0*t13*t3*t57*t61*t92-2.4814019635
     1         976003d+0*gammaaa*t11*t3*t55*t56*t61*t92-1.86105147269820
     2         03d+0*t113*t13*t2*t3*t60*t92+9.305257363491002d-1*t11*t3*
     3         t55*t60*t61*t91+9.305257363491002d-1*t13*t137*t2*t3*t61-1
     4         .2407009817988002d+0*t10*t11*t3*t55*t56+4.962803927195200
     5         7d+0*gammaaa*t10*t107*t108*t136*t3
            t139 = 2.827433388230826d+1*rhoa*t24*t67*t93-2.8274333882308
     1         26d+1*t118*t22*t24*t62*t93+1.413716694115413d+1*t138*t22*
     2         t24*t67
            t140 = 1/rhob**5
            t141 = -2.6666666666666666d+0*t7*t77-1.0666666666666666d+1*g
     1         ammabb*t4*t76
            t142 = 1.2407009817988002d+0*t3*t40*t75*t79*t99-2.4814019635
     1         976003d+0*gammabb*t11*t3*t73*t74*t79*t99-1.86105147269820
     2         03d+0*t128*t3*t33*t40*t78*t99+9.305257363491002d-1*t11*t3
     3         *t73*t78*t79*t98+9.305257363491002d-1*t141*t3*t33*t40*t79
     4         -1.2407009817988002d+0*t11*t3*t38*t73*t74+4.9628039271952
     5         007d+0*gammabb*t107*t123*t140*t3*t38
            t143 = 1.413716694115413d+1*t142*t24*t45*t85+2.8274333882308
     1         26d+1*rhob*t100*t24*t85-2.827433388230826d+1*t100*t133*t2
     2         4*t45*t80
            t144 = 1/t17**3
            t145 = t93**2
            t146 = t92**2
            t147 = 1.8610514726982003d+0*t11*t3*t55*t61*t91*t92+1.861051
     1         4726982003d+0*t117*t13*t3*t4*t61-1.8610514726982003d+0*t1
     2         13*t13*t146*t2*t3-1.8610514726982003d+0*t10*t107*t108*t11
     3         7*t3
            t148 = 1/t17**4
            t149 = t100**2
            t150 = t99**2
            t151 = 1.8610514726982003d+0*t11*t3*t73*t79*t98*t99+1.861051
     1         4726982003d+0*t132*t3*t4*t40*t79-1.8610514726982003d+0*t1
     2         28*t150*t3*t33*t40-1.8610514726982003d+0*t107*t123*t132*t
     3         3*t38
            fnc(iq) = (t41*t53+t14*t32)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t32*t72+t14*t71)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t53*t90+t41*t89)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t32*t97+t14*t96)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t104*t53+t103*t41)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (2*t71*t72+t14*(-3.5461536035682834d-1*
     1         t16*t17*t18*t30*t62*t70+7.092307207136567d-1*t16*t17*t20*
     2         t23*t70-3.5461536035682834d-1*t16*t17*t18*t20*(2.65961520
     3         2676212d-1*t17*t18*t30*t62*t69-5.319230405352424d-1*t17*t
     4         20*t23*t69+2.659615202676212d-1*t17*t18*t20*(1.4147106052
     5         612857d-1*t19*t21*t26*t66*t68-7.073553026306428d-2*t21*t2
     6         3*t26*t62*t68+1.4147106052612857d-1*t21*t27*t62*t66-3.536
     7         776513153214d-2*t115*t21*t23*t27-2.1220659078919285d-1*t1
     8         17*t19*t21*t27-3.536776513153214d-2*t120*t19*t21*t23*t26-
     9         3.536776513153214d-2*t119*t19*t21*t23*t26+t120*t26+t119*t
     :         26)+3.5449077018110403d+0*t26*t54*t64*t68+5.3192304053524
     ;         24d-1*t17*t20*t28*t66+3.5449077018110403d+0*t26*t54*(-3.7
     <         59942411946509d+0*t29*t62*t63-1.8799712059732546d+0*rhoa*
     =         t115*t29*t63+2.8199568089598814d+0*rhoa*t105*t106*t29)-6.
     >         649038006690531d-2*t105*t17*t18*t28*t63-2.659615202676212
     ?         d-1*t17*t23*t28*t30*t62+1.329807601338106d-1*t115*t17*t18
     @         *t28*t30)-7.092307207136567d-1*t16*t17*t20*t31*t66+8.8653
     1         84008920708d-2*t105*t16*t17*t18*t31*t63+3.546153603568283
     2         4d-1*t16*t17*t23*t30*t31*t62-1.7730768017841417d-1*t115*t
     3         16*t17*t18*t30*t31)+t32*(-2.4814019635976003d+0*t13*t3*t5
     4         7*t60*t61+4.9628039271952007d+0*gammaaa*t11*t3*t55*t56*t6
     5         0*t61-9.305257363491002d-1*t114*t13*t2*t3*t61-2.481401963
     6         597601d+0*gammaaa*t10*t11*t110*t3*t55+1.3234143805853868d
     7         +1*t10*t107*t108*t109*t3*t5+1.8610514726982003d+0*t112*t1
     8         13*t13*t2*t3-4.135669939329334d-1*t10*t111*t13*t3-4.13566
     9         9939329334d-1*t1*t111))*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (2*t89*t90+t41*(-3.5461536035682834d-1*
     1         t16*t17*t42*t51*t80*t88+7.092307207136567d-1*t16*t17*t44*
     2         t46*t88-3.5461536035682834d-1*t16*t17*t42*t44*(2.65961520
     3         2676212d-1*t17*t42*t51*t80*t87-5.319230405352424d-1*t17*t
     4         44*t46*t87+2.659615202676212d-1*t17*t42*t44*(1.4147106052
     5         612857d-1*t21*t43*t48*t84*t86-7.073553026306428d-2*t21*t4
     6         6*t48*t80*t86+1.4147106052612857d-1*t21*t49*t80*t84-3.536
     7         776513153214d-2*t130*t21*t46*t49-2.1220659078919285d-1*t1
     8         32*t21*t43*t49-3.536776513153214d-2*t135*t21*t43*t46*t48-
     9         3.536776513153214d-2*t134*t21*t43*t46*t48+t135*t48+t134*t
     :         48)+3.5449077018110403d+0*t48*t54*t82*t86+5.3192304053524
     ;         24d-1*t17*t44*t50*t84+3.5449077018110403d+0*t48*t54*(-3.7
     <         59942411946509d+0*t29*t80*t81-1.8799712059732546d+0*rhob*
     =         t130*t29*t81+2.8199568089598814d+0*rhob*t121*t122*t29)-6.
     >         649038006690531d-2*t121*t17*t42*t50*t81-2.659615202676212
     ?         d-1*t17*t46*t50*t51*t80+1.329807601338106d-1*t130*t17*t42
     @         *t50*t51)-7.092307207136567d-1*t16*t17*t44*t52*t84+8.8653
     1         84008920708d-2*t121*t16*t17*t42*t52*t81+3.546153603568283
     2         4d-1*t16*t17*t46*t51*t52*t80-1.7730768017841417d-1*t130*t
     3         16*t17*t42*t51*t52)+t53*(-2.4814019635976003d+0*t3*t40*t7
     4         5*t78*t79+4.9628039271952007d+0*gammabb*t11*t3*t73*t74*t7
     5         8*t79-9.305257363491002d-1*t129*t3*t33*t40*t79-2.48140196
     6         3597601d+0*gammabb*t11*t125*t3*t38*t73-4.135669939329334d
     7         -1*t126*t3*t38*t40+1.8610514726982003d+0*t127*t128*t3*t33
     8         *t40+1.3234143805853868d+1*t107*t123*t124*t3*t34*t38-4.13
     9         5669939329334d-1*t1*t126))*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t71*t97+t72*t96+t14*(-1.7730768017841
     1         417d-1*t16*t17*t18*t30*t62*t95+3.5461536035682834d-1*t16*
     2         t17*t20*t23*t95-3.5461536035682834d-1*t16*t17*t18*t20*(1.
     3         329807601338106d-1*t17*t18*t30*t62*t94-2.659615202676212d
     4         -1*t17*t20*t23*t94+2.659615202676212d-1*t17*t18*t20*(1.41
     5         3716694115413d+1*t22*t24*t26*t67*t68*t93-5.0d-1*t25*t26*t
     6         68*t93-3.536776513153214d-2*t21*t23*t26*t68*t93-5.0d-1*t2
     7         6*t62*t67*t93+7.073553026306428d-2*t21*t27*t66*t93+1.0d+0
     8         *t18*t25*t26*t93-3.536776513153214d-2*t138*t21*t23*t27-3.
     9         536776513153214d-2*t139*t19*t21*t23*t26+t139*t26)+3.54490
     :         77018110403d+0*t26*t54*(-1.8799712059732546d+0*t29*t63*t9
     ;         3+2.8199568089598814d+0*rhoa*t106*t29*t62*t93-1.879971205
     <         9732546d+0*rhoa*t138*t29*t63)+1.329807601338106d-1*t17*t1
     =         8*t30*t69*t93+5.01149519714857d+1*t22*t24*t26*t54*t64*t67
     >         *t93-6.649038006690531d-2*t17*t18*t28*t62*t63*t93-1.32980
     ?         7601338106d-1*t17*t23*t28*t30*t93+1.329807601338106d-1*t1
     @         38*t17*t18*t28*t30)-1.7730768017841417d-1*t16*t17*t18*t30
     1         *t70*t93+8.865384008920708d-2*t16*t17*t18*t31*t62*t63*t93
     2         +1.7730768017841417d-1*t16*t17*t23*t30*t31*t93-1.77307680
     3         17841417d-1*t138*t16*t17*t18*t30*t31)+t32*(-1.24070098179
     4         88002d+0*t13*t3*t57*t61*t92+2.4814019635976003d+0*gammaaa
     5         *t11*t3*t55*t56*t61*t92+1.8610514726982003d+0*t113*t13*t2
     6         *t3*t60*t92-9.305257363491002d-1*t11*t3*t55*t60*t61*t91-9
     7         .305257363491002d-1*t13*t137*t2*t3*t61+1.2407009817988002
     8         d+0*t10*t11*t3*t55*t56-4.9628039271952007d+0*gammaaa*t10*
     9         t107*t108*t136*t3))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t53*(-1.2407009817988002d+0*t3*t40*t7
     1         5*t79*t99+2.4814019635976003d+0*gammabb*t11*t3*t73*t74*t7
     2         9*t99+1.8610514726982003d+0*t128*t3*t33*t40*t78*t99-9.305
     3         257363491002d-1*t11*t3*t73*t78*t79*t98-9.305257363491002d
     4         -1*t141*t3*t33*t40*t79+1.2407009817988002d+0*t11*t3*t38*t
     5         73*t74-4.9628039271952007d+0*gammabb*t107*t123*t140*t3*t3
     6         8)+t103*t90+t104*t89+t41*(-1.7730768017841417d-1*t100*t16
     7         *t17*t42*t51*t88-3.5461536035682834d-1*t16*t17*t42*t44*(1
     8         .329807601338106d-1*t100*t17*t42*t51*t87+2.65961520267621
     9         2d-1*t17*t42*t44*(1.413716694115413d+1*t100*t24*t45*t48*t
     :         85*t86-5.0d-1*t100*t47*t48*t86-3.536776513153214d-2*t100*
     ;         t21*t46*t48*t86-5.0d-1*t100*t48*t80*t85+7.073553026306428
     <         d-2*t100*t21*t49*t84-3.536776513153214d-2*t142*t21*t46*t4
     =         9+1.0d+0*t100*t42*t47*t48-3.536776513153214d-2*t143*t21*t
     >         43*t46*t48+t143*t48)+5.01149519714857d+1*t100*t24*t45*t48
     ?         *t54*t82*t85+3.5449077018110403d+0*t48*t54*(-1.8799712059
     @         732546d+0*rhob*t142*t29*t81-1.8799712059732546d+0*t100*t2
     1         9*t81+2.8199568089598814d+0*rhob*t100*t122*t29*t80)-6.649
     2         038006690531d-2*t100*t17*t42*t50*t80*t81+1.32980760133810
     3         6d-1*t101*t17*t42*t51*t80-1.329807601338106d-1*t100*t17*t
     4         46*t50*t51+1.329807601338106d-1*t142*t17*t42*t50*t51-2.65
     5         9615202676212d-1*t101*t17*t44*t46)+8.865384008920708d-2*t
     6         100*t16*t17*t42*t52*t80*t81-1.7730768017841417d-1*t102*t1
     7         6*t17*t42*t51*t80+1.7730768017841417d-1*t100*t16*t17*t46*
     8         t51*t52-1.7730768017841417d-1*t142*t16*t17*t42*t51*t52+3.
     9         5461536035682834d-1*t102*t16*t17*t44*t46))*wght+Cmat2(iq,
     :         D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2*t96*t97+t14*(-3.5461536035682834d-
     1         1*t16*t17*t18*t30*t93*t95-3.5461536035682834d-1*t16*t17*t
     2         18*t20*(2.659615202676212d-1*t17*t18*t30*t93*t94+2.659615
     3         202676212d-1*t17*t18*t20*(1.413716694115413d+1*t147*t22*t
     4         24*t26*t67-3.536776513153214d-2*t147*t21*t23*t27-5.0d-1*t
     5         147*t25*t26-3.534291735288533d+1*t118*t145*t22*t24*t26+1.
     6         9985948912206125d+2*t116*t145*t148*t26/t19**4)-9.42146666
     7         9512572d+1*t144*t145*t26*t54*t65/t20**7-6.66432440723758d
     8         +0*rhoa*t147*t26*t29*t54*t63-6.649038006690531d-2*t145*t1
     9         7*t18*t28*t63+9.996486610856369d+0*rhoa*t106*t145*t26*t29
     :         *t54+1.329807601338106d-1*t147*t17*t18*t28*t30)+8.8653840
     ;         08920708d-2*t145*t16*t17*t18*t31*t63-1.7730768017841417d-
     <         1*t147*t16*t17*t18*t30*t31)+t32*(-1.8610514726982003d+0*t
     =         11*t3*t55*t61*t91*t92-1.8610514726982003d+0*t117*t13*t3*t
     >         4*t61+1.8610514726982003d+0*t113*t13*t146*t2*t3+1.8610514
     ?         726982003d+0*t10*t107*t108*t117*t3))*wght+Cmat2(iq,D2_GAA
     @         _GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t53*(-1.8610514726982003d+0*t11*t3*t
     1         73*t79*t98*t99-1.8610514726982003d+0*t132*t3*t4*t40*t79+1
     2         .8610514726982003d+0*t128*t150*t3*t33*t40+1.8610514726982
     3         003d+0*t107*t123*t132*t3*t38)+t41*(-3.5461536035682834d-1
     4         *t16*t17*t42*t44*(2.659615202676212d-1*t17*t42*t44*(1.413
     5         716694115413d+1*t151*t24*t45*t48*t85-3.536776513153214d-2
     6         *t151*t21*t46*t49-5.0d-1*t151*t47*t48-3.534291735288533d+
     7         1*t133*t149*t24*t45*t48+1.9985948912206125d+2*t131*t148*t
     8         149*t48/t43**4)-9.421466669512572d+1*t144*t149*t48*t54*t8
     9         3/t44**7-6.66432440723758d+0*rhob*t151*t29*t48*t54*t81-6.
     :         649038006690531d-2*t149*t17*t42*t50*t81+9.996486610856369
     ;         d+0*rhob*t122*t149*t29*t48*t54+1.329807601338106d-1*t151*
     <         t17*t42*t50*t51+2.659615202676212d-1*t100*t101*t17*t42*t5
     =         1)+8.865384008920708d-2*t149*t16*t17*t42*t52*t81-1.773076
     >         8017841417d-1*t151*t16*t17*t42*t51*t52-3.5461536035682834
     ?         d-1*t100*t102*t16*t17*t42*t51)+2*t103*t104)*wght+Cmat2(iq
     @         ,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = 1/rhoa**6.333333333333333d+0
            t37 = 1/rhoa**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammaaa*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         a*t28*t40*t41
            t43 = rhoa**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhoa*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammaaa*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhoa*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = 1/rhoa**6
            t63 = 1/rhoa**3.3333333333333337d+0
            t64 = 1/rhoa**6.666666666666666d-1
            t65 = t38**2
            t66 = 1/t9**3
            t67 = 9.777777777777779d+0*gammaaa*t7/rhoa**4.66666666666666
     1         7d+0+3.377777777777778d+1*t4*t5/rhoa**7.333333333333333d+
     2         0
            t68 = 9.305257363491002d-1*t13*t2*t3*t39*t67-1.8610514726982
     1         003d+0*t13*t2*t3*t65*t66+4.135669939329334d-1*t10*t13*t3*
     2         t64+4.135669939329334d-1*t1*t64+2.481401963597601d+0*gamm
     3         aaa*t10*t11*t3*t33*t63-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t62+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t69 = rhoa**4
            t70 = 1/t69
            t71 = 1/t18**3
            t72 = -2.827433388230826d+1*t21*t23*t58*t71+1.41371669411541
     1         3d+1*t21*t23*t45*t68+5.654866776461652d+1*rhoa*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t73 = t46**2
            t74 = 1/rhoa**5
            t75 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmaaa*t36*t4
            t76 = 9.305257363491002d-1*t13*t2*t3*t39*t75+4.9628039271952
     1         007d+0*gammaaa*t10*t3*t60*t61*t74-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t66+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t77 = 1.413716694115413d+1*t21*t23*t45*t76-2.827433388230826
     1         d+1*t21*t23*t40*t53*t71+2.827433388230826d+1*rhoa*t23*t45
     2         *t53
            t78 = t53**2
            t79 = t52**2
            t80 = -1.8610514726982003d+0*t13*t2*t3*t66*t79-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t70+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t70+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = t31*t50*wght+t14*t49*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t31*t57*wght+t14*t56*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t14*(-3.5461536035682834d-1*t15*t16*t17
     1         *t19*(2.659615202676212d-1*t16*t17*t19*(-3.53677651315321
     2         4d-2*t18*t20*t22*t25*t73+t25*t73-3.536776513153214d-2*t18
     3         *t20*t22*t25*t72+t25*t72-2.1220659078919285d-1*t18*t20*t2
     4         6*t70-3.536776513153214d-2*t20*t22*t26*t68+1.414710605261
     5         2857d-1*t18*t20*t25*t44*t46-7.073553026306428d-2*t20*t22*
     6         t25*t40*t46+1.4147106052612857d-1*t20*t26*t40*t44)+3.5449
     7         077018110403d+0*t25*t32*(-1.8799712059732546d+0*rhoa*t28*
     8         t41*t68+2.8199568089598814d+0*rhoa*t28*t58*t59-3.75994241
     9         1946509d+0*t28*t40*t41)+1.329807601338106d-1*t16*t17*t27*
     :         t29*t68-6.649038006690531d-2*t16*t17*t27*t41*t58+2.659615
     ;         202676212d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t16
     <         *t19*t22*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.3192
     =         30405352424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*t
     >         22*t27*t29*t40)-1.7730768017841417d-1*t15*t16*t17*t29*t30
     ?         *t68+8.865384008920708d-2*t15*t16*t17*t30*t41*t58-3.54615
     @         36035682834d-1*t15*t16*t17*t29*t40*t48+7.092307207136567d
     1         -1*t15*t16*t19*t22*t48-7.092307207136567d-1*t15*t16*t19*t
     2         30*t44+3.5461536035682834d-1*t15*t16*t22*t29*t30*t40)*wgh
     3         t+t31*(-9.305257363491002d-1*t13*t2*t3*t39*t67+1.86105147
     4         26982003d+0*t13*t2*t3*t65*t66-4.135669939329334d-1*t10*t1
     5         3*t3*t64-4.135669939329334d-1*t1*t64-2.481401963597601d+0
     6         *gammaaa*t10*t11*t3*t33*t63+1.3234143805853868d+1*t10*t3*
     7         t5*t60*t61*t62-2.4814019635976003d+0*t13*t3*t35*t38*t39+4
     8         .9628039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39)*wght
     9         +2*t49*t50*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t14*(-3.5461536035682834d-1*t15*t16*t1
     1         7*t19*(2.659615202676212d-1*t16*t17*t19*(-3.5367765131532
     2         14d-2*t18*t20*t22*t25*t77+t25*t77-3.536776513153214d-2*t2
     3         0*t22*t26*t76+1.413716694115413d+1*t21*t23*t25*t45*t46*t5
     4         3-5.0d-1*t24*t25*t46*t53-3.536776513153214d-2*t20*t22*t25
     5         *t46*t53-5.0d-1*t25*t40*t45*t53+7.073553026306428d-2*t20*
     6         t26*t44*t53+1.0d+0*t17*t24*t25*t53)+3.5449077018110403d+0
     7         *t25*t32*(-1.8799712059732546d+0*rhoa*t28*t41*t76+2.81995
     8         68089598814d+0*rhoa*t28*t40*t53*t59-1.8799712059732546d+0
     9         *t28*t41*t53)+1.329807601338106d-1*t16*t17*t27*t29*t76+1.
     :         329807601338106d-1*t16*t17*t29*t40*t54-2.659615202676212d
     ;         -1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t47*t
     <         53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.6490
     =         38006690531d-2*t16*t17*t27*t40*t41*t53-1.329807601338106d
     >         -1*t16*t22*t27*t29*t53)-1.7730768017841417d-1*t15*t16*t17
     ?         *t29*t30*t76-1.7730768017841417d-1*t15*t16*t17*t29*t40*t5
     @         5+3.5461536035682834d-1*t15*t16*t19*t22*t55-1.77307680178
     1         41417d-1*t15*t16*t17*t29*t48*t53+8.865384008920708d-2*t15
     2         *t16*t17*t30*t40*t41*t53+1.7730768017841417d-1*t15*t16*t2
     3         2*t29*t30*t53)*wght+t31*(-9.305257363491002d-1*t13*t2*t3*
     4         t39*t75-4.9628039271952007d+0*gammaaa*t10*t3*t60*t61*t74+
     5         1.8610514726982003d+0*t13*t2*t3*t38*t52*t66-1.24070098179
     6         88002d+0*t13*t3*t35*t39*t52+2.4814019635976003d+0*gammaaa
     7         *t11*t3*t33*t34*t39*t52-9.305257363491002d-1*t11*t3*t33*t
     8         38*t39*t51+1.2407009817988002d+0*t10*t11*t3*t33*t34)*wght
     9         +t49*t57*wght+t50*t56*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t14*(-3.5461536035682834d-1*t15*t16*t
     1         17*t19*(2.659615202676212d-1*t16*t17*t19*(1.4137166941154
     2         13d+1*t21*t23*t25*t45*t80-3.536776513153214d-2*t20*t22*t2
     3         6*t80-5.0d-1*t24*t25*t80-3.534291735288533d+1*t21*t23*t25
     4         *t71*t78+1.9985948912206125d+2*t25*t69*t78/(t16**4*t18**4
     5         ))-6.66432440723758d+0*rhoa*t25*t28*t32*t41*t80+1.3298076
     6         01338106d-1*t16*t17*t27*t29*t80+9.996486610856369d+0*rhoa
     7         *t25*t28*t32*t59*t78-9.421466669512572d+1*t25*t32*t43*t78
     8         /(t16**3*t19**7)-6.649038006690531d-2*t16*t17*t27*t41*t78
     9         +2.659615202676212d-1*t16*t17*t29*t53*t54)-1.773076801784
     :         1417d-1*t15*t16*t17*t29*t30*t80+8.865384008920708d-2*t15*
     ;         t16*t17*t30*t41*t78-3.5461536035682834d-1*t15*t16*t17*t29
     <         *t53*t55)*wght+t31*(1.8610514726982003d+0*t13*t2*t3*t66*t
     =         79+1.8610514726982003d+0*t10*t3*t60*t61*t70-1.86105147269
     >         82003d+0*t13*t3*t39*t4*t70-1.8610514726982003d+0*t11*t3*t
     ?         33*t39*t51*t52)*wght+2*t56*t57*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhob**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammabb**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhob**2.6666666666666666d+0
            t9 = gammabb*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammabb*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhob
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhob**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhob**2.3333333333333334d+0
            t35 = rhob**3.333333333333333d-1
            t36 = 1/rhob**6.333333333333333d+0
            t37 = 1/rhob**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammabb*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammabb*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         b*t28*t40*t41
            t43 = rhob**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhob*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammabb*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammabb*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhob*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = 1/rhob**6
            t63 = 1/rhob**3.3333333333333337d+0
            t64 = 1/rhob**6.666666666666666d-1
            t65 = t38**2
            t66 = 1/t9**3
            t67 = 9.777777777777779d+0*gammabb*t7/rhob**4.66666666666666
     1         7d+0+3.377777777777778d+1*t4*t5/rhob**7.333333333333333d+
     2         0
            t68 = 9.305257363491002d-1*t13*t2*t3*t39*t67-1.8610514726982
     1         003d+0*t13*t2*t3*t65*t66+4.135669939329334d-1*t10*t13*t3*
     2         t64+4.135669939329334d-1*t1*t64+2.481401963597601d+0*gamm
     3         abb*t10*t11*t3*t33*t63-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t62+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammabb*t11*t3*t33*t34*t38*t39
            t69 = rhob**4
            t70 = 1/t69
            t71 = 1/t18**3
            t72 = -2.827433388230826d+1*t21*t23*t58*t71+1.41371669411541
     1         3d+1*t21*t23*t45*t68+5.654866776461652d+1*rhob*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t73 = t46**2
            t74 = 1/rhob**5
            t75 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmabb*t36*t4
            t76 = 9.305257363491002d-1*t13*t2*t3*t39*t75+4.9628039271952
     1         007d+0*gammabb*t10*t3*t60*t61*t74-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t66+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammabb*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t77 = 1.413716694115413d+1*t21*t23*t45*t76-2.827433388230826
     1         d+1*t21*t23*t40*t53*t71+2.827433388230826d+1*rhob*t23*t45
     2         *t53
            t78 = t53**2
            t79 = t52**2
            t80 = -1.8610514726982003d+0*t13*t2*t3*t66*t79-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t70+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t70+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RB) = t31*t50*wght+t14*t49*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t31*t57*wght+t14*t56*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t14*(-3.5461536035682834d-1*t15*t16*t17
     1         *t19*(2.659615202676212d-1*t16*t17*t19*(-3.53677651315321
     2         4d-2*t18*t20*t22*t25*t73+t25*t73-3.536776513153214d-2*t18
     3         *t20*t22*t25*t72+t25*t72-2.1220659078919285d-1*t18*t20*t2
     4         6*t70-3.536776513153214d-2*t20*t22*t26*t68+1.414710605261
     5         2857d-1*t18*t20*t25*t44*t46-7.073553026306428d-2*t20*t22*
     6         t25*t40*t46+1.4147106052612857d-1*t20*t26*t40*t44)+3.5449
     7         077018110403d+0*t25*t32*(-1.8799712059732546d+0*rhob*t28*
     8         t41*t68+2.8199568089598814d+0*rhob*t28*t58*t59-3.75994241
     9         1946509d+0*t28*t40*t41)+1.329807601338106d-1*t16*t17*t27*
     :         t29*t68-6.649038006690531d-2*t16*t17*t27*t41*t58+2.659615
     ;         202676212d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t16
     <         *t19*t22*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.3192
     =         30405352424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*t
     >         22*t27*t29*t40)-1.7730768017841417d-1*t15*t16*t17*t29*t30
     ?         *t68+8.865384008920708d-2*t15*t16*t17*t30*t41*t58-3.54615
     @         36035682834d-1*t15*t16*t17*t29*t40*t48+7.092307207136567d
     1         -1*t15*t16*t19*t22*t48-7.092307207136567d-1*t15*t16*t19*t
     2         30*t44+3.5461536035682834d-1*t15*t16*t22*t29*t30*t40)*wgh
     3         t+t31*(-9.305257363491002d-1*t13*t2*t3*t39*t67+1.86105147
     4         26982003d+0*t13*t2*t3*t65*t66-4.135669939329334d-1*t10*t1
     5         3*t3*t64-4.135669939329334d-1*t1*t64-2.481401963597601d+0
     6         *gammabb*t10*t11*t3*t33*t63+1.3234143805853868d+1*t10*t3*
     7         t5*t60*t61*t62-2.4814019635976003d+0*t13*t3*t35*t38*t39+4
     8         .9628039271952007d+0*gammabb*t11*t3*t33*t34*t38*t39)*wght
     9         +2*t49*t50*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t14*(-3.5461536035682834d-1*t15*t16*t1
     1         7*t19*(2.659615202676212d-1*t16*t17*t19*(-3.5367765131532
     2         14d-2*t18*t20*t22*t25*t77+t25*t77-3.536776513153214d-2*t2
     3         0*t22*t26*t76+1.413716694115413d+1*t21*t23*t25*t45*t46*t5
     4         3-5.0d-1*t24*t25*t46*t53-3.536776513153214d-2*t20*t22*t25
     5         *t46*t53-5.0d-1*t25*t40*t45*t53+7.073553026306428d-2*t20*
     6         t26*t44*t53+1.0d+0*t17*t24*t25*t53)+3.5449077018110403d+0
     7         *t25*t32*(-1.8799712059732546d+0*rhob*t28*t41*t76+2.81995
     8         68089598814d+0*rhob*t28*t40*t53*t59-1.8799712059732546d+0
     9         *t28*t41*t53)+1.329807601338106d-1*t16*t17*t27*t29*t76+1.
     :         329807601338106d-1*t16*t17*t29*t40*t54-2.659615202676212d
     ;         -1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t47*t
     <         53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.6490
     =         38006690531d-2*t16*t17*t27*t40*t41*t53-1.329807601338106d
     >         -1*t16*t22*t27*t29*t53)-1.7730768017841417d-1*t15*t16*t17
     ?         *t29*t30*t76-1.7730768017841417d-1*t15*t16*t17*t29*t40*t5
     @         5+3.5461536035682834d-1*t15*t16*t19*t22*t55-1.77307680178
     1         41417d-1*t15*t16*t17*t29*t48*t53+8.865384008920708d-2*t15
     2         *t16*t17*t30*t40*t41*t53+1.7730768017841417d-1*t15*t16*t2
     3         2*t29*t30*t53)*wght+t31*(-9.305257363491002d-1*t13*t2*t3*
     4         t39*t75-4.9628039271952007d+0*gammabb*t10*t3*t60*t61*t74+
     5         1.8610514726982003d+0*t13*t2*t3*t38*t52*t66-1.24070098179
     6         88002d+0*t13*t3*t35*t39*t52+2.4814019635976003d+0*gammabb
     7         *t11*t3*t33*t34*t39*t52-9.305257363491002d-1*t11*t3*t33*t
     8         38*t39*t51+1.2407009817988002d+0*t10*t11*t3*t33*t34)*wght
     9         +t49*t57*wght+t50*t56*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t14*(-3.5461536035682834d-1*t15*t16*t
     1         17*t19*(2.659615202676212d-1*t16*t17*t19*(1.4137166941154
     2         13d+1*t21*t23*t25*t45*t80-3.536776513153214d-2*t20*t22*t2
     3         6*t80-5.0d-1*t24*t25*t80-3.534291735288533d+1*t21*t23*t25
     4         *t71*t78+1.9985948912206125d+2*t25*t69*t78/(t16**4*t18**4
     5         ))-6.66432440723758d+0*rhob*t25*t28*t32*t41*t80+1.3298076
     6         01338106d-1*t16*t17*t27*t29*t80+9.996486610856369d+0*rhob
     7         *t25*t28*t32*t59*t78-9.421466669512572d+1*t25*t32*t43*t78
     8         /(t16**3*t19**7)-6.649038006690531d-2*t16*t17*t27*t41*t78
     9         +2.659615202676212d-1*t16*t17*t29*t53*t54)-1.773076801784
     :         1417d-1*t15*t16*t17*t29*t30*t80+8.865384008920708d-2*t15*
     ;         t16*t17*t30*t41*t78-3.5461536035682834d-1*t15*t16*t17*t29
     <         *t53*t55)*wght+t31*(1.8610514726982003d+0*t13*t2*t3*t66*t
     =         79+1.8610514726982003d+0*t10*t3*t60*t61*t70-1.86105147269
     >         82003d+0*t13*t3*t39*t4*t70-1.8610514726982003d+0*t11*t3*t
     ?         33*t39*t51*t52)*wght+2*t56*t57*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_cams12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)-{\it t_1}+1.0\\\\ 
C>   {\it t_4} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,
C>    \sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_2}\,{
C>    \it t_6}\,{\it t_9}+0.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{
C>    \it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}
C>    \,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{
C>    \it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{21}}
C>    \,{\it t_{23}}\,{\it t_{24}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}
C>    \,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{
C>    \it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& 0.9305257363491002\,{\it t_3}\,{\it t_{30}}
C>    \,{\it t_{32}}\,{\it t_{33}}+0.9305257363491002\,{\it t_1}
C>    \,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{
C>    \it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{
C>    \it t_{23}}\,{\it t_{24}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{21}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}
C>    \,\left({\it t_{29}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(
C>    -0.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}
C>    -0.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(
C>    -0.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}
C>    \,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}
C>    \left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{
C>    \it t_{15}}}}\right)+0.2659615202676212\,{\it t_{12}}\,{
C>    \it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}
C>    -{{0.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({
C>    \it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)
C>    +{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{
C>    \it t_{32}}\,{\it t_{33}}-0.9305257363491002\,{\it t_1}\,{
C>    \it t_{30}}\right)\,\left(-0.3546153603568283\,{\it t_{11}}
C>    \,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,
C>    \left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509
C>    \,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)
C>    +0.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}
C>    \,\left({\it t_{38}}-{{0.03536776513153214\,{\it t_{16}}\,{
C>    \it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}
C>    -1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_cams12_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = 1/rhoa**6.333333333333333d+0
            t37 = 1/rhoa**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammaaa*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         a*t28*t40*t41
            t43 = rhoa**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhoa*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammaaa*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhoa*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = rhoa**6
            t63 = 1/t62
            t64 = 1/rhoa**3.3333333333333337d+0
            t65 = 1/rhoa**6.666666666666666d-1
            t66 = t38**2
            t67 = 1/t9**3
            t68 = 1/rhoa**7.333333333333333d+0
            t69 = 1/rhoa**4.666666666666667d+0
            t70 = 9.777777777777779d+0*gammaaa*t69*t7+3.377777777777778d
     1         +1*t4*t5*t68
            t71 = 9.305257363491002d-1*t13*t2*t3*t39*t70-1.8610514726982
     1         003d+0*t13*t2*t3*t66*t67+4.135669939329334d-1*t10*t13*t3*
     2         t65+4.135669939329334d-1*t1*t65+2.481401963597601d+0*gamm
     3         aaa*t10*t11*t3*t33*t64-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t63+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t72 = -1.8799712059732546d+0*rhoa*t28*t41*t71+2.819956808959
     1         8814d+0*rhoa*t28*t58*t59-3.759942411946509d+0*t28*t40*t41
            t73 = rhoa**4
            t74 = 1/t73
            t75 = 1/t18**3
            t76 = -2.827433388230826d+1*t21*t23*t58*t75+1.41371669411541
     1         3d+1*t21*t23*t45*t71+5.654866776461652d+1*rhoa*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t77 = t46**2
            t78 = -3.536776513153214d-2*t18*t20*t22*t25*t77+t25*t77-3.53
     1         6776513153214d-2*t18*t20*t22*t25*t76+t25*t76-2.1220659078
     2         919285d-1*t18*t20*t26*t74-3.536776513153214d-2*t20*t22*t2
     3         6*t71+1.4147106052612857d-1*t18*t20*t25*t44*t46-7.0735530
     4         26306428d-2*t20*t22*t25*t40*t46+1.4147106052612857d-1*t20
     5         *t26*t40*t44
            t79 = 2.659615202676212d-1*t16*t17*t19*t78+3.544907701811040
     1         3d+0*t25*t32*t72+1.329807601338106d-1*t16*t17*t27*t29*t71
     2         -6.649038006690531d-2*t16*t17*t27*t41*t58+2.6596152026762
     3         12d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t16*t19*t2
     4         2*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.31923040535
     5         2424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*t22*t27*
     6         t29*t40
            t80 = -3.5461536035682834d-1*t15*t16*t17*t19*t79-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t71+8.865384008920708d-2*t
     2         15*t16*t17*t30*t41*t58-3.5461536035682834d-1*t15*t16*t17*
     3         t29*t40*t48+7.092307207136567d-1*t15*t16*t19*t22*t48-7.09
     4         2307207136567d-1*t15*t16*t19*t30*t44+3.5461536035682834d-
     5         1*t15*t16*t22*t29*t30*t40
            t81 = -9.305257363491002d-1*t13*t2*t3*t39*t70+1.861051472698
     1         2003d+0*t13*t2*t3*t66*t67-4.135669939329334d-1*t10*t13*t3
     2         *t65-4.135669939329334d-1*t1*t65-2.481401963597601d+0*gam
     3         maaa*t10*t11*t3*t33*t64+1.3234143805853868d+1*t10*t3*t5*t
     4         60*t61*t63-2.4814019635976003d+0*t13*t3*t35*t38*t39+4.962
     5         8039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t82 = rhoa**5
            t83 = 1/t82
            t84 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmaaa*t36*t4
            t85 = 9.305257363491002d-1*t13*t2*t3*t39*t84+4.9628039271952
     1         007d+0*gammaaa*t10*t3*t60*t61*t83-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t67+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t86 = -1.8799712059732546d+0*rhoa*t28*t41*t85+2.819956808959
     1         8814d+0*rhoa*t28*t40*t53*t59-1.8799712059732546d+0*t28*t4
     2         1*t53
            t87 = 1.413716694115413d+1*t21*t23*t45*t85-2.827433388230826
     1         d+1*t21*t23*t40*t53*t75+2.827433388230826d+1*rhoa*t23*t45
     2         *t53
            t88 = -3.536776513153214d-2*t18*t20*t22*t25*t87+t25*t87-3.53
     1         6776513153214d-2*t20*t22*t26*t85+1.413716694115413d+1*t21
     2         *t23*t25*t45*t46*t53-5.0d-1*t24*t25*t46*t53-3.53677651315
     3         3214d-2*t20*t22*t25*t46*t53-5.0d-1*t25*t40*t45*t53+7.0735
     4         53026306428d-2*t20*t26*t44*t53+1.0d+0*t17*t24*t25*t53
            t89 = 2.659615202676212d-1*t16*t17*t19*t88+3.544907701811040
     1         3d+0*t25*t32*t86+1.329807601338106d-1*t16*t17*t27*t29*t85
     2         +1.329807601338106d-1*t16*t17*t29*t40*t54-2.6596152026762
     3         12d-1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t4
     4         7*t53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.6
     5         49038006690531d-2*t16*t17*t27*t40*t41*t53-1.3298076013381
     6         06d-1*t16*t22*t27*t29*t53
            t90 = -3.5461536035682834d-1*t15*t16*t17*t19*t89-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t85-1.7730768017841417d-1*
     2         t15*t16*t17*t29*t40*t55+3.5461536035682834d-1*t15*t16*t19
     3         *t22*t55-1.7730768017841417d-1*t15*t16*t17*t29*t48*t53+8.
     4         865384008920708d-2*t15*t16*t17*t30*t40*t41*t53+1.77307680
     5         17841417d-1*t15*t16*t22*t29*t30*t53
            t91 = -9.305257363491002d-1*t13*t2*t3*t39*t84-4.962803927195
     1         2007d+0*gammaaa*t10*t3*t60*t61*t83+1.8610514726982003d+0*
     2         t13*t2*t3*t38*t52*t67-1.2407009817988002d+0*t13*t3*t35*t3
     3         9*t52+2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t5
     4         2-9.305257363491002d-1*t11*t3*t33*t38*t39*t51+1.240700981
     5         7988002d+0*t10*t11*t3*t33*t34
            t92 = 1/t16**3
            t93 = t53**2
            t94 = 1/t19**7
            t95 = t52**2
            t96 = -1.8610514726982003d+0*t13*t2*t3*t67*t95-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t97 = 1/t16**4
            t98 = 1/t18**4
            t99 = 1.9985948912206125d+2*t25*t73*t93*t97*t98+1.4137166941
     1         15413d+1*t21*t23*t25*t45*t96-3.536776513153214d-2*t20*t22
     2         *t26*t96-5.0d-1*t24*t25*t96-3.534291735288533d+1*t21*t23*
     3         t25*t75*t93
            t100 = 2.659615202676212d-1*t16*t17*t19*t99-6.66432440723758
     1         d+0*rhoa*t25*t28*t32*t41*t96+1.329807601338106d-1*t16*t17
     2         *t27*t29*t96-9.421466669512572d+1*t25*t32*t43*t92*t93*t94
     3         +9.996486610856369d+0*rhoa*t25*t28*t32*t59*t93-6.64903800
     4         6690531d-2*t16*t17*t27*t41*t93+2.659615202676212d-1*t16*t
     5         17*t29*t53*t54
            t101 = -1.7730768017841417d-1*t15*t16*t17*t29*t30*t96+8.8653
     1         84008920708d-2*t15*t16*t17*t30*t41*t93-3.5461536035682834
     2         d-1*t15*t16*t17*t29*t53*t55-3.5461536035682834d-1*t100*t1
     3         5*t16*t17*t19
            t102 = 1.8610514726982003d+0*t13*t2*t3*t67*t95+1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74-1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74-1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t103 = t40**3
            t104 = t11**3
            t105 = gammaaa**3
            t106 = 1/t12**4
            t107 = 1/rhoa**9.666666666666666d+0
            t108 = 1/rhoa**7
            t109 = 1/rhoa**4.333333333333333d+0
            t110 = 1/rhoa**1.6666666666666669d+0
            t111 = t38**3
            t112 = 1/t9**4
            t113 = -4.5629629629629626d+1*gammaaa*t7/rhoa**5.66666666666
     1         6667d+0-2.477037037037037d+2*t4*t5/rhoa**8.33333333333333
     2         4d+0
            t114 = -5.583154418094601d+0*t13*t2*t3*t38*t67*t70+3.7221029
     1         453964005d+0*t13*t3*t35*t39*t70-7.444205890792801d+0*gamm
     2         aaa*t11*t3*t33*t34*t39*t70-7.444205890792801d+0*t13*t3*t3
     3         5*t66*t67+1.4888411781585603d+1*gammaaa*t11*t3*t33*t34*t6
     4         6*t67+1.2407009817988002d+0*t13*t3*t38*t39*t65+7.44420589
     5         0792801d+0*gammaaa*t11*t3*t33*t38*t39*t64-3.9702431417561
     6         607d+1*t3*t38*t39*t5*t60*t61*t63+9.263900664097708d+1*t10
     7         *t108*t3*t5*t60*t61+9.305257363491002d-1*t113*t13*t2*t3*t
     8         39-9.374185195813161d+0*gammaaa*t10*t109*t11*t3*t33+5.583
     9         154418094601d+0*t111*t112*t13*t2*t3-2.7571132928862224d-1
     :         *t10*t110*t13*t3-1.0587315044683095d+2*t10*t104*t105*t106
     ;         *t107*t3-2.7571132928862224d-1*t1*t110
            t115 = 8.482300164692478d+1*t103*t21*t23*t98-8.4823001646924
     1         78d+1*t21*t23*t40*t71*t75-1.6964600329384957d+2*rhoa*t23*
     2         t58*t75+8.482300164692478d+1*rhoa*t23*t45*t71+8.482300164
     3         692478d+1*t23*t40*t45+1.413716694115413d+1*t114*t21*t23*t
     4         45
            t116 = t46**3
            t117 = 1/rhoa**8.666666666666666d+0
            t118 = 9.777777777777779d+0*t69*t7+6.755555555555556d+1*gamm
     1         aaa*t4*t68
            t119 = -3.7221029453964005d+0*t13*t2*t3*t38*t67*t84+2.481401
     1         9635976003d+0*t13*t3*t35*t39*t84-4.9628039271952007d+0*ga
     2         mmaaa*t11*t3*t33*t34*t39*t84+9.925607854390401d+0*gammaaa
     3         *t3*t38*t39*t60*t61*t83-1.8610514726982003d+0*t13*t2*t3*t
     4         52*t67*t70+9.305257363491002d-1*t11*t3*t33*t39*t51*t70-1.
     5         8610514726982003d+0*t11*t3*t33*t51*t66*t67-4.962803927195
     6         2007d+0*t13*t3*t35*t38*t52*t67+9.925607854390401d+0*gamma
     7         aa*t11*t3*t33*t34*t38*t52*t67+5.583154418094601d+0*t112*t
     8         13*t2*t3*t52*t66+4.135669939329334d-1*t13*t3*t39*t52*t65+
     9         2.481401963597601d+0*gammaaa*t11*t3*t33*t39*t52*t64+2.894
     :         9689575305343d+0*t10*t11*t3*t33*t64-1.3234143805853868d+1
     ;         *t3*t39*t5*t52*t60*t61*t63-3.143109153890294d+1*gammaaa*t
     <         10*t3*t60*t61*t63+3.9702431417561607d+1*t10*t104*t106*t11
     =         7*t3*t5-2.4814019635976003d+0*t11*t3*t33*t34*t38*t39+9.30
     >         5257363491002d-1*t118*t13*t2*t3*t39
            t120 = 8.482300164692478d+1*t21*t23*t53*t58*t98-5.6548667764
     1         61652d+1*t21*t23*t40*t75*t85+5.654866776461652d+1*rhoa*t2
     2         3*t45*t85-2.827433388230826d+1*t21*t23*t53*t71*t75-1.1309
     3         733552923305d+2*rhoa*t23*t40*t53*t75+2.827433388230826d+1
     4         *t23*t45*t53+1.413716694115413d+1*t119*t21*t23*t45
            t121 = 1/rhoa**7.666666666666667d+0
            t122 = -2.4814019635976003d+0*t13*t3*t35*t67*t95+4.962803927
     1         1952007d+0*gammaaa*t11*t3*t33*t34*t67*t95+5.5831544180946
     2         01d+0*t112*t13*t2*t3*t38*t95-3.7221029453964005d+0*t13*t2
     3         *t3*t52*t67*t84+1.8610514726982003d+0*t11*t3*t33*t39*t51*
     4         t84+9.925607854390401d+0*gammaaa*t3*t39*t52*t60*t61*t83+7
     5         .444205890792801d+0*t10*t3*t60*t61*t83-7.444205890792801d
     6         +0*t13*t3*t39*t4*t83-3.7221029453964005d+0*t13*t3*t38*t4*
     7         t67*t74-1.8610514726982003d+0*t3*t38*t39*t60*t61*t74-3.72
     8         21029453964005d+0*t11*t3*t33*t38*t51*t52*t67-2.4814019635
     9         97601d+0*t11*t3*t33*t34*t39*t52-4.9628039271952007d+0*gam
     :         maaa*t11*t121*t3*t33*t39*t4-1.4888411781585603d+1*gammaaa
     ;         *t10*t104*t106*t121*t3
            t123 = 8.482300164692478d+1*t21*t23*t40*t93*t98-2.8274333882
     1         30826d+1*t21*t23*t40*t75*t96+2.827433388230826d+1*rhoa*t2
     2         3*t45*t96-5.654866776461652d+1*rhoa*t23*t75*t93-5.6548667
     3         76461652d+1*t21*t23*t53*t75*t85+1.413716694115413d+1*t122
     4         *t21*t23*t45
            t124 = t53**3
            t125 = 1/rhoa**6.666666666666667d+0
            t126 = t52**3
            t127 = -5.583154418094601d+0*t11*t3*t33*t51*t67*t95-1.116630
     1         8836189202d+1*t13*t3*t4*t52*t67*t74-5.583154418094601d+0*
     2         t3*t39*t52*t60*t61*t74+5.583154418094601d+0*t11*t125*t3*t
     3         33*t39*t4+5.583154418094601d+0*t112*t126*t13*t2*t3+5.5831
     4         54418094601d+0*t10*t104*t106*t125*t3
            fnc(iq) = 2*t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = (t31*t50+t14*t49)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (t31*t57+t14*t56)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t31*t81+t14*t80+2*t49*t50)*wght+Amat2(
     1         iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t31*t91+t14*t90+t49*t57+t50*t56)*wght
     1         +Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2*t56*t57+t102*t31+t101*t14)*wght+Cm
     1         at2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (t14*(-3.5461536035682834d-1*t15*t16
     1         *t17*t19*(3.5449077018110403d+0*t25*t32*(-7.0498920223997
     2         04d+0*rhoa*t103*t28*t94+8.459870426879645d+0*rhoa*t28*t40
     3         *t59*t71-5.639913617919763d+0*t28*t41*t71+8.4598704268796
     4         45d+0*t28*t58*t59-1.8799712059732546d+0*rhoa*t114*t28*t41
     5         )+2.659615202676212d-1*t16*t17*t19*(8.488263631567714d-1*
     6         t18*t20*t26*t83+2.1220659078919285d-1*t18*t20*t25*t44*t77
     7         -1.0610329539459643d-1*t20*t22*t25*t40*t77-1.061032953945
     8         9643d-1*t18*t20*t22*t25*t46*t76+3*t25*t46*t76+2.122065907
     9         8919285d-1*t18*t20*t25*t44*t76-1.0610329539459643d-1*t20*
     :         t22*t25*t40*t76-6.366197723675786d-1*t18*t20*t25*t46*t74-
     ;         6.366197723675786d-1*t20*t26*t40*t74-1.0610329539459643d-
     <         1*t20*t22*t25*t46*t71+2.1220659078919285d-1*t20*t26*t44*t
     =         71+4.244131815783857d-1*t20*t25*t40*t44*t46-3.53677651315
     >         3214d-2*t114*t20*t22*t26-3.536776513153214d-2*t116*t18*t2
     ?         0*t22*t25-3.536776513153214d-2*t115*t18*t20*t22*t25+t116*
     @         t25+t115*t25)+3.989422804014318d-1*t16*t17*t29*t40*t78-7.
     1         978845608028636d-1*t16*t19*t22*t78+3.5449077018110403d+0*
     2         t25*t32*t42*t77+3.5449077018110403d+0*t25*t32*t42*t76-1.5
     3         957691216057274d+0*t16*t19*t27*t74+7.089815403622081d+0*t
     4         25*t32*t46*t72+3.989422804014318d-1*t16*t17*t29*t47*t71-1
     5         .994711402007159d-1*t16*t17*t27*t40*t41*t71-3.98942280401
     6         4318d-1*t16*t22*t27*t29*t71+9.973557010035795d-2*t103*t16
     7         *t17*t27*t59-1.994711402007159d-1*t16*t17*t41*t47*t58+1.9
     8         94711402007159d-1*t16*t22*t27*t41*t58+1.5957691216057274d
     9         +0*t16*t19*t44*t47-7.978845608028636d-1*t16*t22*t29*t40*t
     :         47+7.978845608028636d-1*t16*t27*t29*t40*t44+1.32980760133
     ;         8106d-1*t114*t16*t17*t27*t29)-5.319230405352425d-1*t15*t1
     <         6*t17*t29*t40*t79+1.063846081070485d+0*t15*t16*t19*t22*t7
     =         9+2.12769216214097d+0*t15*t16*t19*t30*t74-5.3192304053524
     >         25d-1*t15*t16*t17*t29*t48*t71+2.6596152026762127d-1*t15*t
     ?         16*t17*t30*t40*t41*t71+5.319230405352425d-1*t15*t16*t22*t
     @         29*t30*t71-1.3298076013381063d-1*t103*t15*t16*t17*t30*t59
     1         +2.6596152026762127d-1*t15*t16*t17*t41*t48*t58-2.65961520
     2         26762127d-1*t15*t16*t22*t30*t41*t58-2.12769216214097d+0*t
     3         15*t16*t19*t44*t48+1.063846081070485d+0*t15*t16*t22*t29*t
     4         40*t48-1.063846081070485d+0*t15*t16*t29*t30*t40*t44-1.773
     5         0768017841417d-1*t114*t15*t16*t17*t29*t30)+3*t49*t81+3*t5
     6         0*t80+t31*(5.583154418094601d+0*t13*t2*t3*t38*t67*t70-3.7
     7         221029453964005d+0*t13*t3*t35*t39*t70+7.444205890792801d+
     8         0*gammaaa*t11*t3*t33*t34*t39*t70+7.444205890792801d+0*t13
     9         *t3*t35*t66*t67-1.4888411781585603d+1*gammaaa*t11*t3*t33*
     :         t34*t66*t67-1.2407009817988002d+0*t13*t3*t38*t39*t65-7.44
     ;         4205890792801d+0*gammaaa*t11*t3*t33*t38*t39*t64+3.9702431
     <         417561607d+1*t3*t38*t39*t5*t60*t61*t63-9.263900664097708d
     =         +1*t10*t108*t3*t5*t60*t61-9.305257363491002d-1*t113*t13*t
     >         2*t3*t39+9.374185195813161d+0*gammaaa*t10*t109*t11*t3*t33
     ?         -5.583154418094601d+0*t111*t112*t13*t2*t3+2.7571132928862
     @         224d-1*t10*t110*t13*t3+1.0587315044683095d+2*t10*t104*t10
     1         5*t106*t107*t3+2.7571132928862224d-1*t1*t110))*wght+Amat3
     2         (iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t14*(-3.5461536035682834d-1*t15*t1
     1         6*t17*t19*(3.5449077018110403d+0*t25*t32*(-7.049892022399
     2         704d+0*rhoa*t28*t53*t58*t94+5.639913617919763d+0*rhoa*t28
     3         *t40*t59*t85-3.759942411946509d+0*t28*t41*t85+2.819956808
     4         9598814d+0*rhoa*t28*t53*t59*t71+5.639913617919763d+0*t28*
     5         t40*t53*t59-1.8799712059732546d+0*rhoa*t119*t28*t41)+2.65
     6         9615202676212d-1*t16*t17*t29*t40*t88-5.319230405352424d-1
     7         *t16*t19*t22*t88+2.659615202676212d-1*t16*t17*t19*(-7.073
     8         553026306428d-2*t18*t20*t22*t25*t46*t87+2*t25*t46*t87+1.4
     9         147106052612857d-1*t18*t20*t25*t44*t87-7.073553026306428d
     :         -2*t20*t22*t25*t40*t87-7.073553026306428d-2*t20*t22*t25*t
     ;         46*t85+1.4147106052612857d-1*t20*t26*t44*t85+1.4137166941
     <         15413d+1*t21*t23*t25*t45*t53*t77-5.0d-1*t24*t25*t53*t77-3
     =         .536776513153214d-2*t20*t22*t25*t53*t77+1.413716694115413
     >         d+1*t21*t23*t25*t45*t53*t76-5.0d-1*t24*t25*t53*t76-3.5367
     ?         76513153214d-2*t20*t22*t25*t53*t76-2.1220659078919285d-1*
     @         t20*t26*t53*t74-5.0d-1*t25*t45*t53*t71-1.0d+0*t25*t40*t45
     1         *t46*t53+1.4147106052612857d-1*t20*t25*t44*t46*t53+2.0d+0
     2         *t17*t24*t25*t46*t53+2.0d+0*t17*t25*t40*t45*t53-3.0d+0*t2
     3         2*t24*t25*t53-3.536776513153214d-2*t119*t20*t22*t26-3.536
     4         776513153214d-2*t120*t18*t20*t22*t25+t120*t25)+3.54490770
     5         18110403d+0*t25*t32*t42*t87+3.5449077018110403d+0*t25*t32
     6         *t46*t86+2.659615202676212d-1*t16*t17*t29*t47*t85-1.32980
     7         7601338106d-1*t16*t17*t27*t40*t41*t85-2.659615202676212d-
     8         1*t16*t22*t27*t29*t85+1.329807601338106d-1*t16*t17*t29*t5
     9         3*t78+5.01149519714857d+1*t21*t23*t25*t32*t45*t53*t72+1.3
     :         29807601338106d-1*t16*t17*t29*t54*t71-6.649038006690531d-
     ;         2*t16*t17*t27*t41*t53*t71+9.973557010035795d-2*t16*t17*t2
     <         7*t53*t58*t59-6.649038006690531d-2*t16*t17*t41*t54*t58+5.
     =         319230405352424d-1*t16*t19*t44*t54-2.659615202676212d-1*t
     >         16*t22*t29*t40*t54-1.329807601338106d-1*t16*t17*t40*t41*t
     ?         47*t53-2.659615202676212d-1*t16*t22*t29*t47*t53+5.0114951
     @         9714857d+1*t21*t23*t25*t32*t42*t45*t46*t53+2.659615202676
     1         212d-1*t16*t27*t29*t44*t53+1.329807601338106d-1*t16*t22*t
     2         27*t40*t41*t53+1.329807601338106d-1*t119*t16*t17*t27*t29)
     3         -3.5461536035682834d-1*t15*t16*t17*t29*t40*t89+7.09230720
     4         7136567d-1*t15*t16*t19*t22*t89-3.5461536035682834d-1*t15*
     5         t16*t17*t29*t48*t85+1.7730768017841417d-1*t15*t16*t17*t30
     6         *t40*t41*t85+3.5461536035682834d-1*t15*t16*t22*t29*t30*t8
     7         5-1.7730768017841417d-1*t15*t16*t17*t29*t53*t79-1.7730768
     8         017841417d-1*t15*t16*t17*t29*t55*t71+8.865384008920708d-2
     9         *t15*t16*t17*t30*t41*t53*t71-1.3298076013381063d-1*t15*t1
     :         6*t17*t30*t53*t58*t59+8.865384008920708d-2*t15*t16*t17*t4
     ;         1*t55*t58-7.092307207136567d-1*t15*t16*t19*t44*t55+3.5461
     <         536035682834d-1*t15*t16*t22*t29*t40*t55+1.773076801784141
     =         7d-1*t15*t16*t17*t40*t41*t48*t53+3.5461536035682834d-1*t1
     >         5*t16*t22*t29*t48*t53-3.5461536035682834d-1*t15*t16*t29*t
     ?         30*t44*t53-1.7730768017841417d-1*t15*t16*t22*t30*t40*t41*
     @         t53-1.7730768017841417d-1*t119*t15*t16*t17*t29*t30)+2*t49
     1         *t91+2*t50*t90+t31*(3.7221029453964005d+0*t13*t2*t3*t38*t
     2         67*t84-2.4814019635976003d+0*t13*t3*t35*t39*t84+4.9628039
     3         271952007d+0*gammaaa*t11*t3*t33*t34*t39*t84-9.92560785439
     4         0401d+0*gammaaa*t3*t38*t39*t60*t61*t83+1.8610514726982003
     5         d+0*t13*t2*t3*t52*t67*t70-9.305257363491002d-1*t11*t3*t33
     6         *t39*t51*t70+1.8610514726982003d+0*t11*t3*t33*t51*t66*t67
     7         +4.9628039271952007d+0*t13*t3*t35*t38*t52*t67-9.925607854
     8         390401d+0*gammaaa*t11*t3*t33*t34*t38*t52*t67-5.5831544180
     9         94601d+0*t112*t13*t2*t3*t52*t66-4.135669939329334d-1*t13*
     :         t3*t39*t52*t65-2.481401963597601d+0*gammaaa*t11*t3*t33*t3
     ;         9*t52*t64-2.8949689575305343d+0*t10*t11*t3*t33*t64+1.3234
     <         143805853868d+1*t3*t39*t5*t52*t60*t61*t63+3.1431091538902
     =         94d+1*gammaaa*t10*t3*t60*t61*t63-3.9702431417561607d+1*t1
     >         0*t104*t106*t117*t3*t5+2.4814019635976003d+0*t11*t3*t33*t
     ?         34*t38*t39-9.305257363491002d-1*t118*t13*t2*t3*t39)+t56*t
     @         81+t57*t80)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t14*(-3.5461536035682834d-1*t15*t
     1         16*t17*t19*(1.329807601338106d-1*t16*t17*t29*t40*t99-2.65
     2         9615202676212d-1*t16*t19*t22*t99+2.659615202676212d-1*t16
     3         *t17*t19*(1.9985948912206125d+2*t25*t46*t73*t93*t97*t98-7
     4         .068583470577066d+0*t21*t23*t25*t40*t93*t98+1.41371669411
     5         5413d+1*t21*t23*t25*t45*t46*t96-5.0d-1*t24*t25*t46*t96-3.
     6         536776513153214d-2*t20*t22*t25*t46*t96-5.0d-1*t25*t40*t45
     7         *t96+7.073553026306428d-2*t20*t26*t44*t96+1.0d+0*t17*t24*
     8         t25*t96-3.534291735288533d+1*t21*t23*t25*t46*t75*t93+1.0d
     9         +0*t25*t40*t75*t93+1.413716694115413d+1*rhoa*t23*t25*t75*
     :         t93+2.827433388230826d+1*t21*t23*t25*t45*t53*t87-1.0d+0*t
     ;         24*t25*t53*t87-7.073553026306428d-2*t20*t22*t25*t53*t87-1
     <         .0d+0*t25*t45*t53*t85-3.536776513153214d-2*t122*t20*t22*t
     =         26-3.536776513153214d-2*t123*t18*t20*t22*t25+t123*t25)+7.
     >         084834422688148d+2*t25*t32*t42*t73*t93*t97*t98+3.54490770
     ?         18110403d+0*t25*t32*(2.8199568089598814d+0*rhoa*t28*t40*t
     @         59*t96-1.8799712059732546d+0*t28*t41*t96-7.04989202239970
     1         4d+0*rhoa*t28*t40*t93*t94+2.8199568089598814d+0*t28*t59*t
     2         93+5.639913617919763d+0*rhoa*t28*t53*t59*t85-1.8799712059
     3         732546d+0*rhoa*t122*t28*t41)+1.329807601338106d-1*t16*t17
     4         *t29*t47*t96+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*
     5         t96-6.649038006690531d-2*t16*t17*t27*t40*t41*t96-1.329807
     6         601338106d-1*t16*t22*t27*t29*t96-1.002299039429714d+2*t21
     7         *t23*t25*t32*t42*t75*t93+9.973557010035795d-2*t16*t17*t27
     8         *t40*t59*t93-6.649038006690531d-2*t16*t17*t41*t47*t93+6.6
     9         49038006690531d-2*t16*t22*t27*t41*t93+2.659615202676212d-
     :         1*t16*t17*t29*t53*t88+1.002299039429714d+2*t21*t23*t25*t3
     ;         2*t45*t53*t86+2.659615202676212d-1*t16*t17*t29*t54*t85-1.
     <         329807601338106d-1*t16*t17*t27*t41*t53*t85-1.329807601338
     =         106d-1*t16*t17*t40*t41*t53*t54-2.659615202676212d-1*t16*t
     >         22*t29*t53*t54+1.329807601338106d-1*t122*t16*t17*t27*t29)
     ?         -1.7730768017841417d-1*t15*t16*t17*t29*t48*t96+8.86538400
     @         8920708d-2*t15*t16*t17*t30*t40*t41*t96+1.7730768017841417
     1         d-1*t15*t16*t22*t29*t30*t96-1.3298076013381063d-1*t15*t16
     2         *t17*t30*t40*t59*t93+8.865384008920708d-2*t15*t16*t17*t41
     3         *t48*t93-8.865384008920708d-2*t15*t16*t22*t30*t41*t93-3.5
     4         461536035682834d-1*t15*t16*t17*t29*t53*t89-3.546153603568
     5         2834d-1*t15*t16*t17*t29*t55*t85+1.7730768017841417d-1*t15
     6         *t16*t17*t30*t41*t53*t85+1.7730768017841417d-1*t15*t16*t1
     7         7*t40*t41*t53*t55+3.5461536035682834d-1*t15*t16*t22*t29*t
     8         53*t55-1.7730768017841417d-1*t100*t15*t16*t17*t29*t40-1.7
     9         730768017841417d-1*t122*t15*t16*t17*t29*t30+3.54615360356
     :         82834d-1*t100*t15*t16*t19*t22)+t31*(2.4814019635976003d+0
     ;         *t13*t3*t35*t67*t95-4.9628039271952007d+0*gammaaa*t11*t3*
     <         t33*t34*t67*t95-5.583154418094601d+0*t112*t13*t2*t3*t38*t
     =         95+3.7221029453964005d+0*t13*t2*t3*t52*t67*t84-1.86105147
     >         26982003d+0*t11*t3*t33*t39*t51*t84-9.925607854390401d+0*g
     ?         ammaaa*t3*t39*t52*t60*t61*t83-7.444205890792801d+0*t10*t3
     @         *t60*t61*t83+7.444205890792801d+0*t13*t3*t39*t4*t83+3.722
     1         1029453964005d+0*t13*t3*t38*t4*t67*t74+1.8610514726982003
     2         d+0*t3*t38*t39*t60*t61*t74+3.7221029453964005d+0*t11*t3*t
     3         33*t38*t51*t52*t67+2.481401963597601d+0*t11*t3*t33*t34*t3
     4         9*t52+4.9628039271952007d+0*gammaaa*t11*t121*t3*t33*t39*t
     5         4+1.4888411781585603d+1*gammaaa*t10*t104*t106*t121*t3)+2*
     6         t56*t91+2*t57*t90+t101*t50+t102*t49)*wght+Cmat3(iq,D3_RA_
     7         GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t14*(-3.5461536035682834d-1*t15*
     1         t16*t17*t19*(3.989422804014318d-1*t16*t17*t29*t53*t99+2.6
     2         59615202676212d-1*t16*t17*t19*(5.995784673661838d+2*t25*t
     3         53*t73*t96*t97*t98+1.0602875205865599d+2*t124*t21*t23*t25
     4         *t98-1.2990866792933983d+3*t124*t25*t73*t97/t18**5-1.0602
     5         875205865599d+2*t21*t23*t25*t53*t75*t96+2.825446962492358
     6         4d+3*t124*t25*t62/(t16**6*t18**6)+1.413716694115413d+1*t1
     7         27*t21*t23*t25*t45-3.536776513153214d-2*t127*t20*t22*t26-
     8         5.0d-1*t127*t24*t25)-2.826440000853771d+2*t25*t32*t43*t53
     9         *t92*t94*t96+2.998945983256911d+1*rhoa*t25*t28*t32*t53*t5
     :         9*t96+3.989422804014318d-1*t16*t17*t29*t54*t96-1.99471140
     ;         2007159d-1*t16*t17*t27*t41*t53*t96-2.4991216527140922d+1*
     <         rhoa*t124*t25*t28*t32*t94-1.994711402007159d-1*t16*t17*t4
     =         1*t54*t93+4.710733334756285d+2*t124*t25*t32*t43*t92/t19**
     >         9-1.3319284713741864d+3*t124*t25*t32*t82/(t16**5*t19**11)
     ?         +9.973557010035795d-2*t124*t16*t17*t27*t59-6.664324407237
     @         58d+0*rhoa*t127*t25*t28*t32*t41+1.329807601338106d-1*t127
     1         *t16*t17*t27*t29)-5.319230405352425d-1*t15*t16*t17*t29*t5
     2         5*t96+2.6596152026762127d-1*t15*t16*t17*t30*t41*t53*t96+2
     3         .6596152026762127d-1*t15*t16*t17*t41*t55*t93-1.3298076013
     4         381063d-1*t124*t15*t16*t17*t30*t59-5.319230405352425d-1*t
     5         100*t15*t16*t17*t29*t53-1.7730768017841417d-1*t127*t15*t1
     6         6*t17*t29*t30)+t31*(5.583154418094601d+0*t11*t3*t33*t51*t
     7         67*t95+1.1166308836189202d+1*t13*t3*t4*t52*t67*t74+5.5831
     8         54418094601d+0*t3*t39*t52*t60*t61*t74-5.583154418094601d+
     9         0*t11*t125*t3*t33*t39*t4-5.583154418094601d+0*t112*t126*t
     :         13*t2*t3-5.583154418094601d+0*t10*t104*t106*t125*t3)+3*t1
     ;         01*t57+3*t102*t56)*wght+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = -param(6)
            t16 = param(7)
            t17 = param(8)
            t18 = 1/rhoa
            t19 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t20 = t19**5.0d-1
            t21 = t17**2
            t22 = rhoa**2
            t23 = 1/t22
            t24 = 1/t21
            t25 = 1/t19
            t26 = exp(-1.413716694115413d+1*t22*t24*t25)
            t27 = t26-1.0d+0
            t28 = -3.536776513153214d-2*t19*t21*t23*t27+t26-1.5d+0
            t29 = 1/t17
            t30 = 1/t20
            t31 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         9*t30)+2.659615202676212d-1*t17*t18*t20*t28
            t32 = -3.5461536035682834d-1*t16*t17*t18*t20*t31+t15+1.0d+0
            t33 = rhob**1.3333333333333333d+0
            t34 = gammabb**2
            t35 = 1/rhob**5.333333333333333d+0
            t36 = 1/rhob**2.6666666666666666d+0
            t37 = gammabb*t36*t7+t34*t35*t4+1.0d+0
            t38 = 1.0d+0-1.0d+0/t37
            t39 = gammabb*t11*t36+1.0d+0
            t40 = 1.0d+0-1.0d+0/t39
            t41 = -9.305257363491002d-1*t3*t33*t38*t40-9.305257363491002
     1         d-1*t1*t33
            t42 = 1/rhob
            t43 = 9.305257363491002d-1*t3*t33*t38*t40+9.305257363491002d
     1         -1*t1*t33
            t44 = t43**5.0d-1
            t45 = rhob**2
            t46 = 1/t45
            t47 = 1/t43
            t48 = exp(-1.413716694115413d+1*t24*t45*t47)
            t49 = t48-1.0d+0
            t50 = -3.536776513153214d-2*t21*t43*t46*t49+t48-1.5d+0
            t51 = 1/t44
            t52 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         9*t51)+2.659615202676212d-1*t17*t42*t44*t50
            t53 = -3.5461536035682834d-1*t16*t17*t42*t44*t52+t15+1.0d+0
            t54 = 5.641895835477563d-1
            t55 = 1/t12**2
            t56 = 1/rhoa**2.3333333333333334d+0
            t57 = rhoa**3.333333333333333d-1
            t58 = 1/rhoa**6.333333333333333d+0
            t59 = 1/rhoa**3.6666666666666664d+0
            t60 = -2.6666666666666666d+0*gammaaa*t59*t7-5.33333333333333
     1         3d+0*t4*t5*t58
            t61 = 1/t9**2
            t62 = 9.305257363491002d-1*t13*t2*t3*t60*t61+1.2407009817988
     1         002d+0*t10*t13*t3*t57+1.2407009817988002d+0*t1*t57-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t55*t56
            t63 = 1/t20**3
            t64 = 3.759942411946509d+0*t29*t30-1.8799712059732546d+0*rho
     1         a*t29*t62*t63
            t65 = rhoa**3
            t66 = 1/t65
            t67 = 1/t19**2
            t68 = 1.413716694115413d+1*t22*t24*t62*t67-2.827433388230826
     1         d+1*rhoa*t24*t25
            t69 = -3.536776513153214d-2*t19*t21*t23*t26*t68+t26*t68+7.07
     1         3553026306428d-2*t19*t21*t27*t66-3.536776513153214d-2*t21
     2         *t23*t27*t62
            t70 = 2.659615202676212d-1*t17*t18*t20*t69+3.544907701811040
     1         3d+0*t26*t54*t64+1.329807601338106d-1*t17*t18*t28*t30*t62
     2         -2.659615202676212d-1*t17*t20*t23*t28
            t71 = -3.5461536035682834d-1*t16*t17*t18*t20*t70-1.773076801
     1         7841417d-1*t16*t17*t18*t30*t31*t62+3.5461536035682834d-1*
     2         t16*t17*t20*t23*t31
            t72 = -9.305257363491002d-1*t13*t2*t3*t60*t61-1.240700981798
     1         8002d+0*t10*t13*t3*t57-1.2407009817988002d+0*t1*t57+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t55*t56
            t73 = 1/t39**2
            t74 = 1/rhob**2.3333333333333334d+0
            t75 = rhob**3.333333333333333d-1
            t76 = 1/rhob**6.333333333333333d+0
            t77 = 1/rhob**3.6666666666666664d+0
            t78 = -2.6666666666666666d+0*gammabb*t7*t77-5.33333333333333
     1         3d+0*t34*t4*t76
            t79 = 1/t37**2
            t80 = 9.305257363491002d-1*t3*t33*t40*t78*t79+1.240700981798
     1         8002d+0*t3*t38*t40*t75+1.2407009817988002d+0*t1*t75-2.481
     2         4019635976003d+0*gammabb*t11*t3*t38*t73*t74
            t81 = 1/t44**3
            t82 = 3.759942411946509d+0*t29*t51-1.8799712059732546d+0*rho
     1         b*t29*t80*t81
            t83 = rhob**3
            t84 = 1/t83
            t85 = 1/t43**2
            t86 = 1.413716694115413d+1*t24*t45*t80*t85-2.827433388230826
     1         d+1*rhob*t24*t47
            t87 = -3.536776513153214d-2*t21*t43*t46*t48*t86+t48*t86+7.07
     1         3553026306428d-2*t21*t43*t49*t84-3.536776513153214d-2*t21
     2         *t46*t49*t80
            t88 = 2.659615202676212d-1*t17*t42*t44*t87+3.544907701811040
     1         3d+0*t48*t54*t82+1.329807601338106d-1*t17*t42*t50*t51*t80
     2         -2.659615202676212d-1*t17*t44*t46*t50
            t89 = -3.5461536035682834d-1*t16*t17*t42*t44*t88-1.773076801
     1         7841417d-1*t16*t17*t42*t51*t52*t80+3.5461536035682834d-1*
     2         t16*t17*t44*t46*t52
            t90 = -9.305257363491002d-1*t3*t33*t40*t78*t79-1.24070098179
     1         88002d+0*t3*t38*t40*t75-1.2407009817988002d+0*t1*t75+2.48
     2         14019635976003d+0*gammabb*t11*t3*t38*t73*t74
            t91 = 1/t2
            t92 = t7*t8+2*gammaaa*t4*t6
            t93 = 9.305257363491002d-1*t13*t2*t3*t61*t92+9.3052573634910
     1         02d-1*t10*t11*t3*t55*t91
            t94 = 1.413716694115413d+1*t22*t24*t26*t67*t93-3.53677651315
     1         3214d-2*t21*t23*t27*t93-5.0d-1*t25*t26*t93
            t95 = 2.659615202676212d-1*t17*t18*t20*t94-6.66432440723758d
     1         +0*rhoa*t26*t29*t54*t63*t93+1.329807601338106d-1*t17*t18*
     2         t28*t30*t93
            t96 = -3.5461536035682834d-1*t16*t17*t18*t20*t95-1.773076801
     1         7841417d-1*t16*t17*t18*t30*t31*t93
            t97 = -9.305257363491002d-1*t13*t2*t3*t61*t92-9.305257363491
     1         002d-1*t10*t11*t3*t55*t91
            t98 = 1/t33
            t99 = t36*t7+2*gammabb*t35*t4
            t100 = 9.305257363491002d-1*t3*t33*t40*t79*t99+9.30525736349
     1         1002d-1*t11*t3*t38*t73*t98
            t101 = 1.413716694115413d+1*t100*t24*t45*t48*t85-3.536776513
     1         153214d-2*t100*t21*t46*t49-5.0d-1*t100*t47*t48
            t102 = -6.66432440723758d+0*rhob*t100*t29*t48*t54*t81+1.3298
     1         07601338106d-1*t100*t17*t42*t50*t51+2.659615202676212d-1*
     2         t101*t17*t42*t44
            t103 = -1.7730768017841417d-1*t100*t16*t17*t42*t51*t52-3.546
     1         1536035682834d-1*t102*t16*t17*t42*t44
            t104 = -9.305257363491002d-1*t3*t33*t40*t79*t99-9.3052573634
     1         91002d-1*t11*t3*t38*t73*t98
            t105 = t62**2
            t106 = 1/t20**5
            t107 = t11**2
            t108 = 1/t12**3
            t109 = rhoa**6
            t110 = 1/t109
            t111 = 1/rhoa**3.3333333333333337d+0
            t112 = 1/rhoa**6.666666666666666d-1
            t113 = t60**2
            t114 = 1/t9**3
            t115 = 1/rhoa**7.333333333333333d+0
            t116 = 1/rhoa**4.666666666666667d+0
            t117 = 9.777777777777779d+0*gammaaa*t116*t7+3.37777777777777
     1         8d+1*t115*t4*t5
            t118 = 2.4814019635976003d+0*t13*t3*t57*t60*t61-4.9628039271
     1         952007d+0*gammaaa*t11*t3*t55*t56*t60*t61+9.30525736349100
     2         2d-1*t117*t13*t2*t3*t61+2.481401963597601d+0*gammaaa*t10*
     3         t11*t111*t3*t55-1.3234143805853868d+1*t10*t107*t108*t110*
     4         t3*t5-1.8610514726982003d+0*t113*t114*t13*t2*t3+4.1356699
     5         39329334d-1*t10*t112*t13*t3+4.135669939329334d-1*t1*t112
            t119 = -3.759942411946509d+0*t29*t62*t63-1.8799712059732546d
     1         +0*rhoa*t118*t29*t63+2.8199568089598814d+0*rhoa*t105*t106
     2         *t29
            t120 = rhoa**4
            t121 = 1/t120
            t122 = 1/t19**3
            t123 = 5.654866776461652d+1*rhoa*t24*t62*t67+1.4137166941154
     1         13d+1*t118*t22*t24*t67-2.827433388230826d+1*t24*t25-2.827
     2         433388230826d+1*t105*t122*t22*t24
            t124 = t68**2
            t125 = 1.4147106052612857d-1*t19*t21*t26*t66*t68-7.073553026
     1         306428d-2*t21*t23*t26*t62*t68+1.4147106052612857d-1*t21*t
     2         27*t62*t66-3.536776513153214d-2*t118*t21*t23*t27-2.122065
     3         9078919285d-1*t121*t19*t21*t27-3.536776513153214d-2*t124*
     4         t19*t21*t23*t26-3.536776513153214d-2*t123*t19*t21*t23*t26
     5         +t124*t26+t123*t26
            t126 = 2.659615202676212d-1*t17*t18*t30*t62*t69-5.3192304053
     1         52424d-1*t17*t20*t23*t69+3.5449077018110403d+0*t26*t54*t6
     2         4*t68+5.319230405352424d-1*t17*t20*t28*t66-6.649038006690
     3         531d-2*t105*t17*t18*t28*t63-2.659615202676212d-1*t17*t23*
     4         t28*t30*t62+3.5449077018110403d+0*t119*t26*t54+1.32980760
     5         1338106d-1*t118*t17*t18*t28*t30+2.659615202676212d-1*t125
     6         *t17*t18*t20
            t127 = -3.5461536035682834d-1*t16*t17*t18*t30*t62*t70+7.0923
     1         07207136567d-1*t16*t17*t20*t23*t70-7.092307207136567d-1*t
     2         16*t17*t20*t31*t66+8.865384008920708d-2*t105*t16*t17*t18*
     3         t31*t63+3.5461536035682834d-1*t16*t17*t23*t30*t31*t62-1.7
     4         730768017841417d-1*t118*t16*t17*t18*t30*t31-3.54615360356
     5         82834d-1*t126*t16*t17*t18*t20
            t128 = -2.4814019635976003d+0*t13*t3*t57*t60*t61+4.962803927
     1         1952007d+0*gammaaa*t11*t3*t55*t56*t60*t61-9.3052573634910
     2         02d-1*t117*t13*t2*t3*t61-2.481401963597601d+0*gammaaa*t10
     3         *t11*t111*t3*t55+1.3234143805853868d+1*t10*t107*t108*t110
     4         *t3*t5+1.8610514726982003d+0*t113*t114*t13*t2*t3-4.135669
     5         939329334d-1*t10*t112*t13*t3-4.135669939329334d-1*t1*t112
            t129 = t80**2
            t130 = 1/t44**5
            t131 = 1/t39**3
            t132 = rhob**6
            t133 = 1/t132
            t134 = 1/rhob**3.3333333333333337d+0
            t135 = 1/rhob**6.666666666666666d-1
            t136 = t78**2
            t137 = 1/t37**3
            t138 = 1/rhob**7.333333333333333d+0
            t139 = 1/rhob**4.666666666666667d+0
            t140 = 9.777777777777779d+0*gammabb*t139*t7+3.37777777777777
     1         8d+1*t138*t34*t4
            t141 = 2.4814019635976003d+0*t3*t40*t75*t78*t79-4.9628039271
     1         952007d+0*gammabb*t11*t3*t73*t74*t78*t79+9.30525736349100
     2         2d-1*t140*t3*t33*t40*t79+2.481401963597601d+0*gammabb*t11
     3         *t134*t3*t38*t73+4.135669939329334d-1*t135*t3*t38*t40-1.8
     4         610514726982003d+0*t136*t137*t3*t33*t40-1.323414380585386
     5         8d+1*t107*t131*t133*t3*t34*t38+4.135669939329334d-1*t1*t1
     6         35
            t142 = -3.759942411946509d+0*t29*t80*t81-1.8799712059732546d
     1         +0*rhob*t141*t29*t81+2.8199568089598814d+0*rhob*t129*t130
     2         *t29
            t143 = rhob**4
            t144 = 1/t143
            t145 = 1/t43**3
            t146 = 5.654866776461652d+1*rhob*t24*t80*t85+1.4137166941154
     1         13d+1*t141*t24*t45*t85-2.827433388230826d+1*t24*t47-2.827
     2         433388230826d+1*t129*t145*t24*t45
            t147 = t86**2
            t148 = 1.4147106052612857d-1*t21*t43*t48*t84*t86-7.073553026
     1         306428d-2*t21*t46*t48*t80*t86+1.4147106052612857d-1*t21*t
     2         49*t80*t84-3.536776513153214d-2*t141*t21*t46*t49-2.122065
     3         9078919285d-1*t144*t21*t43*t49-3.536776513153214d-2*t147*
     4         t21*t43*t46*t48-3.536776513153214d-2*t146*t21*t43*t46*t48
     5         +t147*t48+t146*t48
            t149 = 2.659615202676212d-1*t17*t42*t51*t80*t87-5.3192304053
     1         52424d-1*t17*t44*t46*t87+3.5449077018110403d+0*t48*t54*t8
     2         2*t86+5.319230405352424d-1*t17*t44*t50*t84-6.649038006690
     3         531d-2*t129*t17*t42*t50*t81-2.659615202676212d-1*t17*t46*
     4         t50*t51*t80+3.5449077018110403d+0*t142*t48*t54+1.32980760
     5         1338106d-1*t141*t17*t42*t50*t51+2.659615202676212d-1*t148
     6         *t17*t42*t44
            t150 = -3.5461536035682834d-1*t16*t17*t42*t51*t80*t88+7.0923
     1         07207136567d-1*t16*t17*t44*t46*t88-7.092307207136567d-1*t
     2         16*t17*t44*t52*t84+8.865384008920708d-2*t129*t16*t17*t42*
     3         t52*t81+3.5461536035682834d-1*t16*t17*t46*t51*t52*t80-1.7
     4         730768017841417d-1*t141*t16*t17*t42*t51*t52-3.54615360356
     5         82834d-1*t149*t16*t17*t42*t44
            t151 = -2.4814019635976003d+0*t3*t40*t75*t78*t79+4.962803927
     1         1952007d+0*gammabb*t11*t3*t73*t74*t78*t79-9.3052573634910
     2         02d-1*t140*t3*t33*t40*t79-2.481401963597601d+0*gammabb*t1
     3         1*t134*t3*t38*t73-4.135669939329334d-1*t135*t3*t38*t40+1.
     4         8610514726982003d+0*t136*t137*t3*t33*t40+1.32341438058538
     5         68d+1*t107*t131*t133*t3*t34*t38-4.135669939329334d-1*t1*t
     6         135
            t152 = rhoa**5
            t153 = 1/t152
            t154 = -2.6666666666666666d+0*t59*t7-1.0666666666666666d+1*g
     1         ammaaa*t4*t58
            t155 = 1.2407009817988002d+0*t13*t3*t57*t61*t92-2.4814019635
     1         976003d+0*gammaaa*t11*t3*t55*t56*t61*t92-1.86105147269820
     2         03d+0*t114*t13*t2*t3*t60*t92+9.305257363491002d-1*t11*t3*
     3         t55*t60*t61*t91+9.305257363491002d-1*t13*t154*t2*t3*t61-1
     4         .2407009817988002d+0*t10*t11*t3*t55*t56+4.962803927195200
     5         7d+0*gammaaa*t10*t107*t108*t153*t3
            t156 = -1.8799712059732546d+0*t29*t63*t93+2.8199568089598814
     1         d+0*rhoa*t106*t29*t62*t93-1.8799712059732546d+0*rhoa*t155
     2         *t29*t63
            t157 = 2.827433388230826d+1*rhoa*t24*t67*t93-2.8274333882308
     1         26d+1*t122*t22*t24*t62*t93+1.413716694115413d+1*t155*t22*
     2         t24*t67
            t158 = 1.413716694115413d+1*t22*t24*t26*t67*t68*t93-5.0d-1*t
     1         25*t26*t68*t93-3.536776513153214d-2*t21*t23*t26*t68*t93-5
     2         .0d-1*t26*t62*t67*t93+7.073553026306428d-2*t21*t27*t66*t9
     3         3+1.0d+0*t18*t25*t26*t93-3.536776513153214d-2*t155*t21*t2
     4         3*t27-3.536776513153214d-2*t157*t19*t21*t23*t26+t157*t26
            t159 = 1.329807601338106d-1*t17*t18*t30*t62*t94-2.6596152026
     1         76212d-1*t17*t20*t23*t94+1.329807601338106d-1*t17*t18*t30
     2         *t69*t93+5.01149519714857d+1*t22*t24*t26*t54*t64*t67*t93-
     3         6.649038006690531d-2*t17*t18*t28*t62*t63*t93-1.3298076013
     4         38106d-1*t17*t23*t28*t30*t93+3.5449077018110403d+0*t156*t
     5         26*t54+1.329807601338106d-1*t155*t17*t18*t28*t30+2.659615
     6         202676212d-1*t158*t17*t18*t20
            t160 = -1.7730768017841417d-1*t16*t17*t18*t30*t62*t95+3.5461
     1         536035682834d-1*t16*t17*t20*t23*t95-1.7730768017841417d-1
     2         *t16*t17*t18*t30*t70*t93+8.865384008920708d-2*t16*t17*t18
     3         *t31*t62*t63*t93+1.7730768017841417d-1*t16*t17*t23*t30*t3
     4         1*t93-1.7730768017841417d-1*t155*t16*t17*t18*t30*t31-3.54
     5         61536035682834d-1*t159*t16*t17*t18*t20
            t161 = -1.2407009817988002d+0*t13*t3*t57*t61*t92+2.481401963
     1         5976003d+0*gammaaa*t11*t3*t55*t56*t61*t92+1.8610514726982
     2         003d+0*t114*t13*t2*t3*t60*t92-9.305257363491002d-1*t11*t3
     3         *t55*t60*t61*t91-9.305257363491002d-1*t13*t154*t2*t3*t61+
     4         1.2407009817988002d+0*t10*t11*t3*t55*t56-4.96280392719520
     5         07d+0*gammaaa*t10*t107*t108*t153*t3
            t162 = rhob**5
            t163 = 1/t162
            t164 = -2.6666666666666666d+0*t7*t77-1.0666666666666666d+1*g
     1         ammabb*t4*t76
            t165 = 1.2407009817988002d+0*t3*t40*t75*t79*t99-2.4814019635
     1         976003d+0*gammabb*t11*t3*t73*t74*t79*t99-1.86105147269820
     2         03d+0*t137*t3*t33*t40*t78*t99+9.305257363491002d-1*t11*t3
     3         *t73*t78*t79*t98+9.305257363491002d-1*t164*t3*t33*t40*t79
     4         -1.2407009817988002d+0*t11*t3*t38*t73*t74+4.9628039271952
     5         007d+0*gammabb*t107*t131*t163*t3*t38
            t166 = -1.8799712059732546d+0*rhob*t165*t29*t81-1.8799712059
     1         732546d+0*t100*t29*t81+2.8199568089598814d+0*rhob*t100*t1
     2         30*t29*t80
            t167 = 1.413716694115413d+1*t165*t24*t45*t85+2.8274333882308
     1         26d+1*rhob*t100*t24*t85-2.827433388230826d+1*t100*t145*t2
     2         4*t45*t80
            t168 = 1.413716694115413d+1*t100*t24*t45*t48*t85*t86-5.0d-1*
     1         t100*t47*t48*t86-3.536776513153214d-2*t100*t21*t46*t48*t8
     2         6-5.0d-1*t100*t48*t80*t85+7.073553026306428d-2*t100*t21*t
     3         49*t84-3.536776513153214d-2*t165*t21*t46*t49+1.0d+0*t100*
     4         t42*t47*t48-3.536776513153214d-2*t167*t21*t43*t46*t48+t16
     5         7*t48
            t169 = 1.329807601338106d-1*t100*t17*t42*t51*t87+5.011495197
     1         14857d+1*t100*t24*t45*t48*t54*t82*t85-6.649038006690531d-
     2         2*t100*t17*t42*t50*t80*t81+1.329807601338106d-1*t101*t17*
     3         t42*t51*t80+3.5449077018110403d+0*t166*t48*t54-1.32980760
     4         1338106d-1*t100*t17*t46*t50*t51+1.329807601338106d-1*t165
     5         *t17*t42*t50*t51-2.659615202676212d-1*t101*t17*t44*t46+2.
     6         659615202676212d-1*t168*t17*t42*t44
            t170 = -1.7730768017841417d-1*t100*t16*t17*t42*t51*t88+8.865
     1         384008920708d-2*t100*t16*t17*t42*t52*t80*t81-1.7730768017
     2         841417d-1*t102*t16*t17*t42*t51*t80+1.7730768017841417d-1*
     3         t100*t16*t17*t46*t51*t52-1.7730768017841417d-1*t16*t165*t
     4         17*t42*t51*t52+3.5461536035682834d-1*t102*t16*t17*t44*t46
     5         -3.5461536035682834d-1*t16*t169*t17*t42*t44
            t171 = -1.2407009817988002d+0*t3*t40*t75*t79*t99+2.481401963
     1         5976003d+0*gammabb*t11*t3*t73*t74*t79*t99+1.8610514726982
     2         003d+0*t137*t3*t33*t40*t78*t99-9.305257363491002d-1*t11*t
     3         3*t73*t78*t79*t98-9.305257363491002d-1*t164*t3*t33*t40*t7
     4         9+1.2407009817988002d+0*t11*t3*t38*t73*t74-4.962803927195
     5         2007d+0*gammabb*t107*t131*t163*t3*t38
            t172 = 1/t17**3
            t173 = t93**2
            t174 = 1/t20**7
            t175 = t92**2
            t176 = 1.8610514726982003d+0*t11*t3*t55*t61*t91*t92+1.861051
     1         4726982003d+0*t121*t13*t3*t4*t61-1.8610514726982003d+0*t1
     2         14*t13*t175*t2*t3-1.8610514726982003d+0*t10*t107*t108*t12
     3         1*t3
            t177 = 1/t17**4
            t178 = 1/t19**4
            t179 = 1.413716694115413d+1*t176*t22*t24*t26*t67-3.536776513
     1         153214d-2*t176*t21*t23*t27-5.0d-1*t176*t25*t26-3.53429173
     2         5288533d+1*t122*t173*t22*t24*t26+1.9985948912206125d+2*t1
     3         20*t173*t177*t178*t26
            t180 = 2.659615202676212d-1*t17*t18*t30*t93*t94-9.4214666695
     1         12572d+1*t172*t173*t174*t26*t54*t65-6.66432440723758d+0*r
     2         hoa*t176*t26*t29*t54*t63-6.649038006690531d-2*t17*t173*t1
     3         8*t28*t63+9.996486610856369d+0*rhoa*t106*t173*t26*t29*t54
     4         +1.329807601338106d-1*t17*t176*t18*t28*t30+2.659615202676
     5         212d-1*t17*t179*t18*t20
            t181 = -3.5461536035682834d-1*t16*t17*t18*t30*t93*t95+8.8653
     1         84008920708d-2*t16*t17*t173*t18*t31*t63-1.773076801784141
     2         7d-1*t16*t17*t176*t18*t30*t31-3.5461536035682834d-1*t16*t
     3         17*t18*t180*t20
            t182 = -1.8610514726982003d+0*t11*t3*t55*t61*t91*t92-1.86105
     1         14726982003d+0*t121*t13*t3*t4*t61+1.8610514726982003d+0*t
     2         114*t13*t175*t2*t3+1.8610514726982003d+0*t10*t107*t108*t1
     3         21*t3
            t183 = t100**2
            t184 = 1/t44**7
            t185 = t99**2
            t186 = 1.8610514726982003d+0*t11*t3*t73*t79*t98*t99+1.861051
     1         4726982003d+0*t144*t3*t4*t40*t79-1.8610514726982003d+0*t1
     2         37*t185*t3*t33*t40-1.8610514726982003d+0*t107*t131*t144*t
     3         3*t38
            t187 = 1/t43**4
            t188 = 1.413716694115413d+1*t186*t24*t45*t48*t85-3.536776513
     1         153214d-2*t186*t21*t46*t49-5.0d-1*t186*t47*t48-3.53429173
     2         5288533d+1*t145*t183*t24*t45*t48+1.9985948912206125d+2*t1
     3         43*t177*t183*t187*t48
            t189 = -9.421466669512572d+1*t172*t183*t184*t48*t54*t83-6.66
     1         432440723758d+0*rhob*t186*t29*t48*t54*t81-6.6490380066905
     2         31d-2*t17*t183*t42*t50*t81+9.996486610856369d+0*rhob*t130
     3         *t183*t29*t48*t54+1.329807601338106d-1*t17*t186*t42*t50*t
     4         51+2.659615202676212d-1*t100*t101*t17*t42*t51+2.659615202
     5         676212d-1*t17*t188*t42*t44
            t190 = 8.865384008920708d-2*t16*t17*t183*t42*t52*t81-1.77307
     1         68017841417d-1*t16*t17*t186*t42*t51*t52-3.546153603568283
     2         4d-1*t100*t102*t16*t17*t42*t51-3.5461536035682834d-1*t16*
     3         t17*t189*t42*t44
            t191 = -1.8610514726982003d+0*t11*t3*t73*t79*t98*t99-1.86105
     1         14726982003d+0*t144*t3*t4*t40*t79+1.8610514726982003d+0*t
     2         137*t185*t3*t33*t40+1.8610514726982003d+0*t107*t131*t144*
     3         t3*t38
            t192 = t62**3
            t193 = t11**3
            t194 = gammaaa**3
            t195 = 1/t12**4
            t196 = 1/rhoa**9.666666666666666d+0
            t197 = 1/rhoa**7
            t198 = 1/rhoa**4.333333333333333d+0
            t199 = 1/rhoa**1.6666666666666669d+0
            t200 = t60**3
            t201 = 1/t9**4
            t202 = -4.5629629629629626d+1*gammaaa*t7/rhoa**5.66666666666
     1         6667d+0-2.477037037037037d+2*t4*t5/rhoa**8.33333333333333
     2         4d+0
            t203 = 7.444205890792801d+0*gammaaa*t11*t111*t3*t55*t60*t61-
     1         3.9702431417561607d+1*t107*t108*t110*t3*t5*t60*t61+1.2407
     2         009817988002d+0*t112*t13*t3*t60*t61+3.7221029453964005d+0
     3         *t117*t13*t3*t57*t61-7.444205890792801d+0*gammaaa*t11*t11
     4         7*t3*t55*t56*t61+9.305257363491002d-1*t13*t2*t202*t3*t61-
     5         5.583154418094601d+0*t114*t117*t13*t2*t3*t60-7.4442058907
     6         92801d+0*t113*t114*t13*t3*t57+1.4888411781585603d+1*gamma
     7         aa*t11*t113*t114*t3*t55*t56-9.374185195813161d+0*gammaaa*
     8         t10*t11*t198*t3*t55+9.263900664097708d+1*t10*t107*t108*t1
     9         97*t3*t5+5.583154418094601d+0*t13*t2*t200*t201*t3-2.75711
     :         32928862224d-1*t10*t13*t199*t3-1.0587315044683095d+2*t10*
     ;         t193*t194*t195*t196*t3-2.7571132928862224d-1*t1*t199
            t204 = 8.482300164692478d+1*t24*t62*t67+1.413716694115413d+1
     1         *t203*t22*t24*t67+8.482300164692478d+1*rhoa*t118*t24*t67-
     2         8.482300164692478d+1*t118*t122*t22*t24*t62+8.482300164692
     3         478d+1*t178*t192*t22*t24-1.6964600329384957d+2*rhoa*t105*
     4         t122*t24
            t205 = t68**3
            t206 = t80**3
            t207 = gammabb**3
            t208 = 1/t39**4
            t209 = 1/rhob**9.666666666666666d+0
            t210 = 1/rhob**7
            t211 = 1/rhob**4.333333333333333d+0
            t212 = 1/rhob**1.6666666666666669d+0
            t213 = t78**3
            t214 = 1/t37**4
            t215 = -4.5629629629629626d+1*gammabb*t7/rhob**5.66666666666
     1         6667d+0-2.477037037037037d+2*t34*t4/rhob**8.3333333333333
     2         34d+0
            t216 = 7.444205890792801d+0*gammabb*t11*t134*t3*t73*t78*t79+
     1         1.2407009817988002d+0*t135*t3*t40*t78*t79-3.9702431417561
     2         607d+1*t107*t131*t133*t3*t34*t78*t79+3.7221029453964005d+
     3         0*t140*t3*t40*t75*t79-7.444205890792801d+0*gammabb*t11*t1
     4         40*t3*t73*t74*t79+9.305257363491002d-1*t215*t3*t33*t40*t7
     5         9-5.583154418094601d+0*t137*t140*t3*t33*t40*t78-7.4442058
     6         90792801d+0*t136*t137*t3*t40*t75+1.4888411781585603d+1*ga
     7         mmabb*t11*t136*t137*t3*t73*t74-9.374185195813161d+0*gamma
     8         bb*t11*t211*t3*t38*t73-2.7571132928862224d-1*t212*t3*t38*
     9         t40+5.583154418094601d+0*t213*t214*t3*t33*t40+9.263900664
     :         097708d+1*t107*t131*t210*t3*t34*t38-1.0587315044683095d+2
     ;         *t193*t207*t208*t209*t3*t38-2.7571132928862224d-1*t1*t212
            t217 = 8.482300164692478d+1*t24*t80*t85+1.413716694115413d+1
     1         *t216*t24*t45*t85+8.482300164692478d+1*rhob*t141*t24*t85-
     2         8.482300164692478d+1*t141*t145*t24*t45*t80+8.482300164692
     3         478d+1*t187*t206*t24*t45-1.6964600329384957d+2*rhob*t129*
     4         t145*t24
            t218 = t86**3
            t219 = 1/rhoa**8.666666666666666d+0
            t220 = 9.777777777777779d+0*t116*t7+6.755555555555556d+1*gam
     1         maaa*t115*t4
            t221 = 2.481401963597601d+0*gammaaa*t11*t111*t3*t55*t61*t92-
     1         1.3234143805853868d+1*t107*t108*t110*t3*t5*t61*t92+4.1356
     2         69939329334d-1*t112*t13*t3*t61*t92-4.9628039271952007d+0*
     3         t114*t13*t3*t57*t60*t92+9.925607854390401d+0*gammaaa*t11*
     4         t114*t3*t55*t56*t60*t92+5.583154418094601d+0*t113*t13*t2*
     5         t201*t3*t92-1.8610514726982003d+0*t114*t117*t13*t2*t3*t92
     6         +9.305257363491002d-1*t11*t117*t3*t55*t61*t91-1.861051472
     7         6982003d+0*t11*t113*t114*t3*t55*t91-2.4814019635976003d+0
     8         *t11*t3*t55*t56*t60*t61+9.925607854390401d+0*gammaaa*t107
     9         *t108*t153*t3*t60*t61+2.4814019635976003d+0*t13*t154*t3*t
     :         57*t61-4.9628039271952007d+0*gammaaa*t11*t154*t3*t55*t56*
     ;         t61+9.305257363491002d-1*t13*t2*t220*t3*t61-3.72210294539
     <         64005d+0*t114*t13*t154*t2*t3*t60+2.8949689575305343d+0*t1
     =         0*t11*t111*t3*t55+3.9702431417561607d+1*t10*t193*t195*t21
     >         9*t3*t5-3.143109153890294d+1*gammaaa*t10*t107*t108*t110*t
     ?         3
            t222 = 2.827433388230826d+1*t24*t67*t93-1.1309733552923305d+
     1         2*rhoa*t122*t24*t62*t93+8.482300164692478d+1*t105*t178*t2
     2         2*t24*t93-2.827433388230826d+1*t118*t122*t22*t24*t93+1.41
     3         3716694115413d+1*t22*t221*t24*t67+5.654866776461652d+1*rh
     4         oa*t155*t24*t67-5.654866776461652d+1*t122*t155*t22*t24*t6
     5         2
            t223 = 1/rhob**8.666666666666666d+0
            t224 = 9.777777777777779d+0*t139*t7+6.755555555555556d+1*gam
     1         mabb*t138*t4
            t225 = 2.481401963597601d+0*gammabb*t11*t134*t3*t73*t79*t99+
     1         4.135669939329334d-1*t135*t3*t40*t79*t99-1.32341438058538
     2         68d+1*t107*t131*t133*t3*t34*t79*t99-4.9628039271952007d+0
     3         *t137*t3*t40*t75*t78*t99+9.925607854390401d+0*gammabb*t11
     4         *t137*t3*t73*t74*t78*t99+5.583154418094601d+0*t136*t214*t
     5         3*t33*t40*t99-1.8610514726982003d+0*t137*t140*t3*t33*t40*
     6         t99+9.305257363491002d-1*t11*t140*t3*t73*t79*t98-1.861051
     7         4726982003d+0*t11*t136*t137*t3*t73*t98-2.4814019635976003
     8         d+0*t11*t3*t73*t74*t78*t79+9.925607854390401d+0*gammabb*t
     9         107*t131*t163*t3*t78*t79+2.4814019635976003d+0*t164*t3*t4
     :         0*t75*t79-4.9628039271952007d+0*gammabb*t11*t164*t3*t73*t
     ;         74*t79+9.305257363491002d-1*t224*t3*t33*t40*t79-3.7221029
     <         453964005d+0*t137*t164*t3*t33*t40*t78+2.8949689575305343d
     =         +0*t11*t134*t3*t38*t73+3.9702431417561607d+1*t193*t208*t2
     >         23*t3*t34*t38-3.143109153890294d+1*gammabb*t107*t131*t133
     ?         *t3*t38
            t226 = 1.413716694115413d+1*t225*t24*t45*t85+5.6548667764616
     1         52d+1*rhob*t165*t24*t85+2.827433388230826d+1*t100*t24*t85
     2         -5.654866776461652d+1*t145*t165*t24*t45*t80-1.13097335529
     3         23305d+2*rhob*t100*t145*t24*t80+8.482300164692478d+1*t100
     4         *t129*t187*t24*t45-2.827433388230826d+1*t100*t141*t145*t2
     5         4*t45
            t227 = 1/rhoa**7.666666666666667d+0
            t228 = -3.7221029453964005d+0*t11*t114*t3*t55*t60*t91*t92-2.
     1         481401963597601d+0*t11*t3*t55*t56*t61*t92+9.9256078543904
     2         01d+0*gammaaa*t107*t108*t153*t3*t61*t92-3.722102945396400
     3         5d+0*t114*t13*t154*t2*t3*t92+1.8610514726982003d+0*t11*t1
     4         54*t3*t55*t61*t91-1.8610514726982003d+0*t107*t108*t121*t3
     5         *t60*t61-4.9628039271952007d+0*gammaaa*t11*t227*t3*t4*t55
     6         *t61-7.444205890792801d+0*t13*t153*t3*t4*t61-3.7221029453
     7         964005d+0*t114*t121*t13*t3*t4*t60+5.583154418094601d+0*t1
     8         3*t175*t2*t201*t3*t60-2.4814019635976003d+0*t114*t13*t175
     9         *t3*t57+4.9628039271952007d+0*gammaaa*t11*t114*t175*t3*t5
     :         5*t56-1.4888411781585603d+1*gammaaa*t10*t193*t195*t227*t3
     ;         +7.444205890792801d+0*t10*t107*t108*t153*t3
            t229 = -5.654866776461652d+1*t122*t155*t22*t24*t93+1.4137166
     1         94115413d+1*t22*t228*t24*t67+2.827433388230826d+1*rhoa*t1
     2         76*t24*t67+8.482300164692478d+1*t173*t178*t22*t24*t62-2.8
     3         27433388230826d+1*t122*t176*t22*t24*t62-5.654866776461652
     4         d+1*rhoa*t122*t173*t24
            t230 = 1/rhob**7.666666666666667d+0
            t231 = -3.7221029453964005d+0*t11*t137*t3*t73*t78*t98*t99-2.
     1         481401963597601d+0*t11*t3*t73*t74*t79*t99+9.9256078543904
     2         01d+0*gammabb*t107*t131*t163*t3*t79*t99-3.722102945396400
     3         5d+0*t137*t164*t3*t33*t40*t99+1.8610514726982003d+0*t11*t
     4         164*t3*t73*t79*t98-1.8610514726982003d+0*t107*t131*t144*t
     5         3*t78*t79-4.9628039271952007d+0*gammabb*t11*t230*t3*t4*t7
     6         3*t79-7.444205890792801d+0*t163*t3*t4*t40*t79-3.722102945
     7         3964005d+0*t137*t144*t3*t4*t40*t78+5.583154418094601d+0*t
     8         185*t214*t3*t33*t40*t78-2.4814019635976003d+0*t137*t185*t
     9         3*t40*t75+4.9628039271952007d+0*gammabb*t11*t137*t185*t3*
     :         t73*t74-1.4888411781585603d+1*gammabb*t193*t208*t230*t3*t
     ;         38+7.444205890792801d+0*t107*t131*t163*t3*t38
            t232 = 1.413716694115413d+1*t231*t24*t45*t85+2.8274333882308
     1         26d+1*rhob*t186*t24*t85+8.482300164692478d+1*t183*t187*t2
     2         4*t45*t80-2.827433388230826d+1*t145*t186*t24*t45*t80-5.65
     3         4866776461652d+1*t100*t145*t165*t24*t45-5.654866776461652
     4         d+1*rhob*t145*t183*t24
            t233 = 1/t17**5
            t234 = t93**3
            t235 = 1/rhoa**6.666666666666667d+0
            t236 = t92**3
            t237 = -5.583154418094601d+0*t107*t108*t121*t3*t61*t92-1.116
     1         6308836189202d+1*t114*t121*t13*t3*t4*t92-5.58315441809460
     2         1d+0*t11*t114*t175*t3*t55*t91+5.583154418094601d+0*t11*t2
     3         35*t3*t4*t55*t61+5.583154418094601d+0*t13*t2*t201*t236*t3
     4         +5.583154418094601d+0*t10*t193*t195*t235*t3
            t238 = 1/t17**6
            t239 = t100**3
            t240 = 1/rhob**6.666666666666667d+0
            t241 = t99**3
            t242 = -5.583154418094601d+0*t107*t131*t144*t3*t79*t99-1.116
     1         6308836189202d+1*t137*t144*t3*t4*t40*t99-5.58315441809460
     2         1d+0*t11*t137*t185*t3*t73*t98+5.583154418094601d+0*t11*t2
     3         40*t3*t4*t73*t79+5.583154418094601d+0*t214*t241*t3*t33*t4
     4         0+5.583154418094601d+0*t193*t208*t240*t3*t38
            fnc(iq) = (t41*t53+t14*t32)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t32*t72+t14*t71)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t53*t90+t41*t89)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (t32*t97+t14*t96)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (t104*t53+t103*t41)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (2*t71*t72+t128*t32+t127*t14)*wght+Amat
     1         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (2*t89*t90+t151*t53+t150*t41)*wght+Amat
     1         2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t71*t97+t72*t96+t161*t32+t14*t160)*wg
     1         ht+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t103*t90+t104*t89+t171*t53+t170*t41)*
     1         wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2*t96*t97+t182*t32+t14*t181)*wght+Cm
     1         at2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t191*t53+t190*t41+2*t103*t104)*wght+
     1         Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (3*t127*t72+3*t128*t71+t14*(-2.12769
     1         216214097d+0*t16*t17*t20*t66*t70+2.6596152026762127d-1*t1
     2         05*t16*t17*t18*t63*t70+1.063846081070485d+0*t16*t17*t23*t
     3         30*t62*t70-5.319230405352425d-1*t118*t16*t17*t18*t30*t70-
     4         3.5461536035682834d-1*t16*t17*t18*t20*(1.5957691216057274
     5         d+0*t17*t20*t66*t69-1.994711402007159d-1*t105*t17*t18*t63
     6         *t69-7.978845608028636d-1*t17*t23*t30*t62*t69+3.989422804
     7         014318d-1*t118*t17*t18*t30*t69+2.659615202676212d-1*t17*t
     8         18*t20*(4.244131815783857d-1*t21*t26*t62*t66*t68-1.061032
     9         9539459643d-1*t123*t19*t21*t23*t26*t68-1.0610329539459643
     :         d-1*t118*t21*t23*t26*t68-6.366197723675786d-1*t121*t19*t2
     ;         1*t26*t68+3*t123*t26*t68+2.1220659078919285d-1*t118*t21*t
     <         27*t66+2.1220659078919285d-1*t124*t19*t21*t26*t66+2.12206
     =         59078919285d-1*t123*t19*t21*t26*t66-6.366197723675786d-1*
     >         t121*t21*t27*t62-1.0610329539459643d-1*t124*t21*t23*t26*t
     ?         62-1.0610329539459643d-1*t123*t21*t23*t26*t62-3.536776513
     @         153214d-2*t203*t21*t23*t27+8.488263631567714d-1*t153*t19*
     1         t21*t27-3.536776513153214d-2*t19*t205*t21*t23*t26-3.53677
     2         6513153214d-2*t19*t204*t21*t23*t26+t205*t26+t204*t26)+7.0
     3         89815403622081d+0*t119*t26*t54*t68+7.978845608028636d-1*t
     4         17*t28*t30*t62*t66+3.5449077018110403d+0*t124*t26*t54*t64
     5         +3.5449077018110403d+0*t123*t26*t54*t64+3.544907701811040
     6         3d+0*t26*t54*(-1.8799712059732546d+0*rhoa*t203*t29*t63-5.
     7         639913617919763d+0*t118*t29*t63+8.459870426879645d+0*rhoa
     8         *t106*t118*t29*t62-7.049892022399704d+0*rhoa*t174*t192*t2
     9         9+8.459870426879645d+0*t105*t106*t29)-1.994711402007159d-
     :         1*t118*t17*t18*t28*t62*t63+1.994711402007159d-1*t105*t17*
     ;         t23*t28*t63+3.989422804014318d-1*t125*t17*t18*t30*t62-3.9
     <         89422804014318d-1*t118*t17*t23*t28*t30+1.329807601338106d
     =         -1*t17*t18*t203*t28*t30-1.5957691216057274d+0*t121*t17*t2
     >         0*t28+9.973557010035795d-2*t106*t17*t18*t192*t28-7.978845
     ?         608028636d-1*t125*t17*t20*t23)-1.063846081070485d+0*t16*t
     @         17*t30*t31*t62*t66+2.6596152026762127d-1*t118*t16*t17*t18
     1         *t31*t62*t63-2.6596152026762127d-1*t105*t16*t17*t23*t31*t
     2         63-5.319230405352425d-1*t126*t16*t17*t18*t30*t62+5.319230
     3         405352425d-1*t118*t16*t17*t23*t30*t31-1.7730768017841417d
     4         -1*t16*t17*t18*t203*t30*t31+2.12769216214097d+0*t121*t16*
     5         t17*t20*t31-1.3298076013381063d-1*t106*t16*t17*t18*t192*t
     6         31+1.063846081070485d+0*t126*t16*t17*t20*t23)+t32*(-7.444
     7         205890792801d+0*gammaaa*t11*t111*t3*t55*t60*t61+3.9702431
     8         417561607d+1*t107*t108*t110*t3*t5*t60*t61-1.2407009817988
     9         002d+0*t112*t13*t3*t60*t61-3.7221029453964005d+0*t117*t13
     :         *t3*t57*t61+7.444205890792801d+0*gammaaa*t11*t117*t3*t55*
     ;         t56*t61-9.305257363491002d-1*t13*t2*t202*t3*t61+5.5831544
     <         18094601d+0*t114*t117*t13*t2*t3*t60+7.444205890792801d+0*
     =         t113*t114*t13*t3*t57-1.4888411781585603d+1*gammaaa*t11*t1
     >         13*t114*t3*t55*t56+9.374185195813161d+0*gammaaa*t10*t11*t
     ?         198*t3*t55-9.263900664097708d+1*t10*t107*t108*t197*t3*t5-
     @         5.583154418094601d+0*t13*t2*t200*t201*t3+2.75711329288622
     1         24d-1*t10*t13*t199*t3+1.0587315044683095d+2*t10*t193*t194
     2         *t195*t196*t3+2.7571132928862224d-1*t1*t199))*wght+Amat3(
     3         iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (3*t150*t90+3*t151*t89+t41*(-2.12769
     1         216214097d+0*t16*t17*t44*t84*t88+2.6596152026762127d-1*t1
     2         29*t16*t17*t42*t81*t88+1.063846081070485d+0*t16*t17*t46*t
     3         51*t80*t88-5.319230405352425d-1*t141*t16*t17*t42*t51*t88-
     4         3.5461536035682834d-1*t16*t17*t42*t44*(1.5957691216057274
     5         d+0*t17*t44*t84*t87-1.994711402007159d-1*t129*t17*t42*t81
     6         *t87-7.978845608028636d-1*t17*t46*t51*t80*t87+3.989422804
     7         014318d-1*t141*t17*t42*t51*t87+2.659615202676212d-1*t17*t
     8         42*t44*(4.244131815783857d-1*t21*t48*t80*t84*t86-1.061032
     9         9539459643d-1*t146*t21*t43*t46*t48*t86-1.0610329539459643
     :         d-1*t141*t21*t46*t48*t86-6.366197723675786d-1*t144*t21*t4
     ;         3*t48*t86+3*t146*t48*t86+2.1220659078919285d-1*t141*t21*t
     <         49*t84+2.1220659078919285d-1*t147*t21*t43*t48*t84+2.12206
     =         59078919285d-1*t146*t21*t43*t48*t84-6.366197723675786d-1*
     >         t144*t21*t49*t80-1.0610329539459643d-1*t147*t21*t46*t48*t
     ?         80-1.0610329539459643d-1*t146*t21*t46*t48*t80-3.536776513
     @         153214d-2*t21*t216*t46*t49+8.488263631567714d-1*t163*t21*
     1         t43*t49-3.536776513153214d-2*t21*t218*t43*t46*t48-3.53677
     2         6513153214d-2*t21*t217*t43*t46*t48+t218*t48+t217*t48)+7.0
     3         89815403622081d+0*t142*t48*t54*t86+7.978845608028636d-1*t
     4         17*t50*t51*t80*t84+3.5449077018110403d+0*t147*t48*t54*t82
     5         +3.5449077018110403d+0*t146*t48*t54*t82+3.544907701811040
     6         3d+0*t48*t54*(-1.8799712059732546d+0*rhob*t216*t29*t81-5.
     7         639913617919763d+0*t141*t29*t81+8.459870426879645d+0*rhob
     8         *t130*t141*t29*t80-7.049892022399704d+0*rhob*t184*t206*t2
     9         9+8.459870426879645d+0*t129*t130*t29)-1.994711402007159d-
     :         1*t141*t17*t42*t50*t80*t81+1.994711402007159d-1*t129*t17*
     ;         t46*t50*t81+3.989422804014318d-1*t148*t17*t42*t51*t80-3.9
     <         89422804014318d-1*t141*t17*t46*t50*t51+1.329807601338106d
     =         -1*t17*t216*t42*t50*t51-1.5957691216057274d+0*t144*t17*t4
     >         4*t50+9.973557010035795d-2*t130*t17*t206*t42*t50-7.978845
     ?         608028636d-1*t148*t17*t44*t46)-1.063846081070485d+0*t16*t
     @         17*t51*t52*t80*t84+2.6596152026762127d-1*t141*t16*t17*t42
     1         *t52*t80*t81-2.6596152026762127d-1*t129*t16*t17*t46*t52*t
     2         81-5.319230405352425d-1*t149*t16*t17*t42*t51*t80+5.319230
     3         405352425d-1*t141*t16*t17*t46*t51*t52-1.7730768017841417d
     4         -1*t16*t17*t216*t42*t51*t52+2.12769216214097d+0*t144*t16*
     5         t17*t44*t52-1.3298076013381063d-1*t130*t16*t17*t206*t42*t
     6         52+1.063846081070485d+0*t149*t16*t17*t44*t46)+t53*(-7.444
     7         205890792801d+0*gammabb*t11*t134*t3*t73*t78*t79-1.2407009
     8         817988002d+0*t135*t3*t40*t78*t79+3.9702431417561607d+1*t1
     9         07*t131*t133*t3*t34*t78*t79-3.7221029453964005d+0*t140*t3
     :         *t40*t75*t79+7.444205890792801d+0*gammabb*t11*t140*t3*t73
     ;         *t74*t79-9.305257363491002d-1*t215*t3*t33*t40*t79+5.58315
     <         4418094601d+0*t137*t140*t3*t33*t40*t78+7.444205890792801d
     =         +0*t136*t137*t3*t40*t75-1.4888411781585603d+1*gammabb*t11
     >         *t136*t137*t3*t73*t74+9.374185195813161d+0*gammabb*t11*t2
     ?         11*t3*t38*t73+2.7571132928862224d-1*t212*t3*t38*t40-5.583
     @         154418094601d+0*t213*t214*t3*t33*t40-9.263900664097708d+1
     1         *t107*t131*t210*t3*t34*t38+1.0587315044683095d+2*t193*t20
     2         7*t208*t209*t3*t38+2.7571132928862224d-1*t1*t212))*wght+A
     3         mat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t127*t97+t128*t96+t14*(-7.09230720
     1         7136567d-1*t16*t17*t20*t66*t95+8.865384008920708d-2*t105*
     2         t16*t17*t18*t63*t95+3.5461536035682834d-1*t16*t17*t23*t30
     3         *t62*t95-1.7730768017841417d-1*t118*t16*t17*t18*t30*t95-3
     4         .5461536035682834d-1*t16*t17*t18*t20*(5.319230405352424d-
     5         1*t17*t20*t66*t94-6.649038006690531d-2*t105*t17*t18*t63*t
     6         94-2.659615202676212d-1*t17*t23*t30*t62*t94+1.32980760133
     7         8106d-1*t118*t17*t18*t30*t94+2.659615202676212d-1*t17*t18
     8         *t20*(-1.0d+0*t26*t62*t67*t68*t93+1.4147106052612857d-1*t
     9         21*t26*t66*t68*t93+2.0d+0*t18*t25*t26*t68*t93+2.0d+0*t18*
     :         t26*t62*t67*t93+1.413716694115413d+1*t124*t22*t24*t26*t67
     ;         *t93+1.413716694115413d+1*t123*t22*t24*t26*t67*t93-5.0d-1
     <         *t118*t26*t67*t93-2.1220659078919285d-1*t121*t21*t27*t93-
     =         3.0d+0*t23*t25*t26*t93-5.0d-1*t124*t25*t26*t93-5.0d-1*t12
     >         3*t25*t26*t93-3.536776513153214d-2*t124*t21*t23*t26*t93-3
     ?         .536776513153214d-2*t123*t21*t23*t26*t93-7.07355302630642
     @         8d-2*t157*t19*t21*t23*t26*t68-7.073553026306428d-2*t155*t
     1         21*t23*t26*t68+2*t157*t26*t68+1.4147106052612857d-1*t155*
     2         t21*t27*t66+1.4147106052612857d-1*t157*t19*t21*t26*t66-7.
     3         073553026306428d-2*t157*t21*t23*t26*t62-3.536776513153214
     4         d-2*t21*t221*t23*t27-3.536776513153214d-2*t19*t21*t222*t2
     5         3*t26+t222*t26)+3.5449077018110403d+0*t26*t54*(5.63991361
     6         7919763d+0*t106*t29*t62*t93-7.049892022399704d+0*rhoa*t10
     7         5*t174*t29*t93+2.8199568089598814d+0*rhoa*t106*t118*t29*t
     8         93-1.8799712059732546d+0*rhoa*t221*t29*t63-3.759942411946
     9         509d+0*t155*t29*t63+5.639913617919763d+0*rhoa*t106*t155*t
     :         29*t62)-1.329807601338106d-1*t17*t18*t62*t63*t69*t93-2.65
     ;         9615202676212d-1*t17*t23*t30*t69*t93+5.01149519714857d+1*
     <         t22*t24*t26*t54*t64*t67*t68*t93+5.01149519714857d+1*t119*
     =         t22*t24*t26*t54*t67*t93+2.659615202676212d-1*t17*t28*t30*
     >         t66*t93+1.329807601338106d-1*t17*t23*t28*t62*t63*t93-6.64
     ?         9038006690531d-2*t118*t17*t18*t28*t63*t93+1.3298076013381
     @         06d-1*t125*t17*t18*t30*t93+9.973557010035795d-2*t105*t106
     1         *t17*t18*t28*t93+2.659615202676212d-1*t155*t17*t18*t30*t6
     2         9+3.5449077018110403d+0*t156*t26*t54*t68+3.54490770181104
     3         03d+0*t157*t26*t54*t64-1.329807601338106d-1*t155*t17*t18*
     4         t28*t62*t63+2.659615202676212d-1*t158*t17*t18*t30*t62-2.6
     5         59615202676212d-1*t155*t17*t23*t28*t30+1.329807601338106d
     6         -1*t17*t18*t221*t28*t30-5.319230405352424d-1*t158*t17*t20
     7         *t23)+1.7730768017841417d-1*t16*t17*t18*t62*t63*t70*t93+3
     8         .5461536035682834d-1*t16*t17*t23*t30*t70*t93-3.5461536035
     9         682834d-1*t16*t17*t30*t31*t66*t93-1.7730768017841417d-1*t
     :         16*t17*t23*t31*t62*t63*t93+8.865384008920708d-2*t118*t16*
     ;         t17*t18*t31*t63*t93-1.3298076013381063d-1*t105*t106*t16*t
     <         17*t18*t31*t93-1.7730768017841417d-1*t126*t16*t17*t18*t30
     =         *t93-3.5461536035682834d-1*t155*t16*t17*t18*t30*t70+1.773
     >         0768017841417d-1*t155*t16*t17*t18*t31*t62*t63-3.546153603
     ?         5682834d-1*t159*t16*t17*t18*t30*t62+3.5461536035682834d-1
     @         *t155*t16*t17*t23*t30*t31-1.7730768017841417d-1*t16*t17*t
     1         18*t221*t30*t31+7.092307207136567d-1*t159*t16*t17*t20*t23
     2         )+t32*(-2.481401963597601d+0*gammaaa*t11*t111*t3*t55*t61*
     3         t92+1.3234143805853868d+1*t107*t108*t110*t3*t5*t61*t92-4.
     4         135669939329334d-1*t112*t13*t3*t61*t92+4.9628039271952007
     5         d+0*t114*t13*t3*t57*t60*t92-9.925607854390401d+0*gammaaa*
     6         t11*t114*t3*t55*t56*t60*t92-5.583154418094601d+0*t113*t13
     7         *t2*t201*t3*t92+1.8610514726982003d+0*t114*t117*t13*t2*t3
     8         *t92-9.305257363491002d-1*t11*t117*t3*t55*t61*t91+1.86105
     9         14726982003d+0*t11*t113*t114*t3*t55*t91+2.481401963597600
     :         3d+0*t11*t3*t55*t56*t60*t61-9.925607854390401d+0*gammaaa*
     ;         t107*t108*t153*t3*t60*t61-2.4814019635976003d+0*t13*t154*
     <         t3*t57*t61+4.9628039271952007d+0*gammaaa*t11*t154*t3*t55*
     =         t56*t61-9.305257363491002d-1*t13*t2*t220*t3*t61+3.7221029
     >         453964005d+0*t114*t13*t154*t2*t3*t60-2.8949689575305343d+
     ?         0*t10*t11*t111*t3*t55-3.9702431417561607d+1*t10*t193*t195
     @         *t219*t3*t5+3.143109153890294d+1*gammaaa*t10*t107*t108*t1
     1         10*t3)+2*t160*t72+2*t161*t71)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t53*(-2.481401963597601d+0*gammabb
     1         *t11*t134*t3*t73*t79*t99-4.135669939329334d-1*t135*t3*t40
     2         *t79*t99+1.3234143805853868d+1*t107*t131*t133*t3*t34*t79*
     3         t99+4.9628039271952007d+0*t137*t3*t40*t75*t78*t99-9.92560
     4         7854390401d+0*gammabb*t11*t137*t3*t73*t74*t78*t99-5.58315
     5         4418094601d+0*t136*t214*t3*t33*t40*t99+1.8610514726982003
     6         d+0*t137*t140*t3*t33*t40*t99-9.305257363491002d-1*t11*t14
     7         0*t3*t73*t79*t98+1.8610514726982003d+0*t11*t136*t137*t3*t
     8         73*t98+2.4814019635976003d+0*t11*t3*t73*t74*t78*t79-9.925
     9         607854390401d+0*gammabb*t107*t131*t163*t3*t78*t79-2.48140
     :         19635976003d+0*t164*t3*t40*t75*t79+4.9628039271952007d+0*
     ;         gammabb*t11*t164*t3*t73*t74*t79-9.305257363491002d-1*t224
     <         *t3*t33*t40*t79+3.7221029453964005d+0*t137*t164*t3*t33*t4
     =         0*t78-2.8949689575305343d+0*t11*t134*t3*t38*t73-3.9702431
     >         417561607d+1*t193*t208*t223*t3*t34*t38+3.143109153890294d
     ?         +1*gammabb*t107*t131*t133*t3*t38)+2*t170*t90+2*t171*t89+t
     @         41*(1.7730768017841417d-1*t100*t16*t17*t42*t80*t81*t88+3.
     1         5461536035682834d-1*t100*t16*t17*t46*t51*t88-3.5461536035
     2         682834d-1*t16*t165*t17*t42*t51*t88-3.5461536035682834d-1*
     3         t16*t17*t42*t44*(-1.329807601338106d-1*t100*t17*t42*t80*t
     4         81*t87-2.659615202676212d-1*t100*t17*t46*t51*t87+2.659615
     5         202676212d-1*t165*t17*t42*t51*t87+2.659615202676212d-1*t1
     6         7*t42*t44*(-1.0d+0*t100*t48*t80*t85*t86+1.414710605261285
     7         7d-1*t100*t21*t48*t84*t86+2.0d+0*t100*t42*t47*t48*t86-7.0
     8         73553026306428d-2*t167*t21*t43*t46*t48*t86-7.073553026306
     9         428d-2*t165*t21*t46*t48*t86+2*t167*t48*t86+2.0d+0*t100*t4
     :         2*t48*t80*t85+1.413716694115413d+1*t100*t147*t24*t45*t48*
     ;         t85+1.413716694115413d+1*t100*t146*t24*t45*t48*t85-5.0d-1
     <         *t100*t141*t48*t85+1.4147106052612857d-1*t165*t21*t49*t84
     =         +1.4147106052612857d-1*t167*t21*t43*t48*t84-7.07355302630
     >         6428d-2*t167*t21*t46*t48*t80-3.536776513153214d-2*t21*t22
     ?         5*t46*t49-2.1220659078919285d-1*t100*t144*t21*t49-3.0d+0*
     @         t100*t46*t47*t48-5.0d-1*t100*t147*t47*t48-5.0d-1*t100*t14
     1         6*t47*t48-3.536776513153214d-2*t21*t226*t43*t46*t48-3.536
     2         776513153214d-2*t100*t147*t21*t46*t48-3.536776513153214d-
     3         2*t100*t146*t21*t46*t48+t226*t48)+5.01149519714857d+1*t10
     4         0*t24*t45*t48*t54*t82*t85*t86+3.5449077018110403d+0*t166*
     5         t48*t54*t86+5.01149519714857d+1*t100*t142*t24*t45*t48*t54
     6         *t85+2.659615202676212d-1*t100*t17*t50*t51*t84+5.31923040
     7         5352424d-1*t101*t17*t44*t84+3.5449077018110403d+0*t167*t4
     8         8*t54*t82+3.5449077018110403d+0*t48*t54*(-1.8799712059732
     9         546d+0*rhob*t225*t29*t81-3.759942411946509d+0*t165*t29*t8
     :         1+5.639913617919763d+0*rhob*t130*t165*t29*t80+5.639913617
     ;         919763d+0*t100*t130*t29*t80-7.049892022399704d+0*rhob*t10
     <         0*t129*t184*t29+2.8199568089598814d+0*rhob*t100*t130*t141
     =         *t29)+1.329807601338106d-1*t100*t17*t46*t50*t80*t81-1.329
     >         807601338106d-1*t165*t17*t42*t50*t80*t81-6.64903800669053
     ?         1d-2*t100*t141*t17*t42*t50*t81-6.649038006690531d-2*t101*
     @         t129*t17*t42*t81-2.659615202676212d-1*t101*t17*t46*t51*t8
     1         0+2.659615202676212d-1*t168*t17*t42*t51*t80-2.65961520267
     2         6212d-1*t165*t17*t46*t50*t51+1.329807601338106d-1*t17*t22
     3         5*t42*t50*t51+1.329807601338106d-1*t100*t148*t17*t42*t51+
     4         1.329807601338106d-1*t101*t141*t17*t42*t51+9.973557010035
     5         795d-2*t100*t129*t130*t17*t42*t50-5.319230405352424d-1*t1
     6         68*t17*t44*t46)-3.5461536035682834d-1*t100*t16*t17*t51*t5
     7         2*t84-7.092307207136567d-1*t102*t16*t17*t44*t84-1.7730768
     8         017841417d-1*t100*t16*t17*t46*t52*t80*t81+1.7730768017841
     9         417d-1*t16*t165*t17*t42*t52*t80*t81+8.865384008920708d-2*
     :         t100*t141*t16*t17*t42*t52*t81+8.865384008920708d-2*t102*t
     ;         129*t16*t17*t42*t81+3.5461536035682834d-1*t102*t16*t17*t4
     <         6*t51*t80-3.5461536035682834d-1*t16*t169*t17*t42*t51*t80+
     =         3.5461536035682834d-1*t16*t165*t17*t46*t51*t52-1.77307680
     >         17841417d-1*t16*t17*t225*t42*t51*t52-1.3298076013381063d-
     ?         1*t100*t129*t130*t16*t17*t42*t52-1.7730768017841417d-1*t1
     @         00*t149*t16*t17*t42*t51-1.7730768017841417d-1*t102*t141*t
     1         16*t17*t42*t51+7.092307207136567d-1*t16*t169*t17*t44*t46)
     2         +t103*t151+t104*t150)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (2*t160*t97+2*t161*t96+t14*(1.7730
     1         768017841417d-1*t16*t17*t18*t62*t63*t93*t95+3.54615360356
     2         82834d-1*t16*t17*t23*t30*t93*t95-3.5461536035682834d-1*t1
     3         55*t16*t17*t18*t30*t95-3.5461536035682834d-1*t16*t17*t18*
     4         t20*(-1.329807601338106d-1*t17*t18*t62*t63*t93*t94-2.6596
     5         15202676212d-1*t17*t23*t30*t93*t94+2.659615202676212d-1*t
     6         155*t17*t18*t30*t94+2.659615202676212d-1*t17*t18*t20*(2.8
     7         27433388230826d+1*t157*t22*t24*t26*t67*t93-1.0d+0*t155*t2
     8         6*t67*t93-1.0d+0*t157*t25*t26*t93-7.073553026306428d-2*t1
     9         57*t21*t23*t26*t93+1.413716694115413d+1*t176*t22*t24*t26*
     :         t67*t68-5.0d-1*t176*t25*t26*t68-3.534291735288533d+1*t122
     ;         *t173*t22*t24*t26*t68-3.536776513153214d-2*t176*t21*t23*t
     <         26*t68+1.9985948912206125d+2*t120*t173*t177*t178*t26*t68-
     =         5.0d-1*t176*t26*t62*t67+7.073553026306428d-2*t176*t21*t27
     >         *t66-7.068583470577066d+0*t173*t178*t22*t24*t26*t62+1.0d+
     ?         0*t122*t173*t26*t62-3.536776513153214d-2*t21*t228*t23*t27
     @         +1.0d+0*t176*t18*t25*t26+1.413716694115413d+1*rhoa*t122*t
     1         173*t24*t26-3.536776513153214d-2*t19*t21*t229*t23*t26+t22
     2         9*t26)+3.5449077018110403d+0*t26*t54*(5.639913617919763d+
     3         0*rhoa*t106*t155*t29*t93-1.8799712059732546d+0*rhoa*t228*
     4         t29*t63-1.8799712059732546d+0*t176*t29*t63+2.819956808959
     5         8814d+0*rhoa*t106*t176*t29*t62-7.049892022399704d+0*rhoa*
     6         t173*t174*t29*t62+2.8199568089598814d+0*t106*t173*t29)+1.
     7         002299039429714d+2*t156*t22*t24*t26*t54*t67*t93-1.3298076
     8         01338106d-1*t155*t17*t18*t28*t63*t93+2.659615202676212d-1
     9         *t158*t17*t18*t30*t93-6.649038006690531d-2*t17*t173*t18*t
     :         63*t69+1.329807601338106d-1*t17*t176*t18*t30*t69+5.011495
     ;         19714857d+1*t176*t22*t24*t26*t54*t64*t67-1.00229903942971
     <         4d+2*t122*t173*t22*t24*t26*t54*t64+7.084834422688148d+2*t
     =         120*t173*t177*t178*t26*t54*t64-6.649038006690531d-2*t17*t
     >         176*t18*t28*t62*t63+6.649038006690531d-2*t17*t173*t23*t28
     ?         *t63+1.329807601338106d-1*t17*t179*t18*t30*t62+9.97355701
     @         0035795d-2*t106*t17*t173*t18*t28*t62-1.329807601338106d-1
     1         *t17*t176*t23*t28*t30+1.329807601338106d-1*t17*t18*t228*t
     2         28*t30-2.659615202676212d-1*t17*t179*t20*t23)+1.773076801
     3         7841417d-1*t155*t16*t17*t18*t31*t63*t93-3.546153603568283
     4         4d-1*t159*t16*t17*t18*t30*t93+8.865384008920708d-2*t16*t1
     5         7*t173*t18*t63*t70-1.7730768017841417d-1*t16*t17*t176*t18
     6         *t30*t70+8.865384008920708d-2*t16*t17*t176*t18*t31*t62*t6
     7         3-8.865384008920708d-2*t16*t17*t173*t23*t31*t63-1.3298076
     8         013381063d-1*t106*t16*t17*t173*t18*t31*t62-1.773076801784
     9         1417d-1*t16*t17*t18*t180*t30*t62+1.7730768017841417d-1*t1
     :         6*t17*t176*t23*t30*t31-1.7730768017841417d-1*t16*t17*t18*
     ;         t228*t30*t31+3.5461536035682834d-1*t16*t17*t180*t20*t23)+
     <         t32*(3.7221029453964005d+0*t11*t114*t3*t55*t60*t91*t92+2.
     =         481401963597601d+0*t11*t3*t55*t56*t61*t92-9.9256078543904
     >         01d+0*gammaaa*t107*t108*t153*t3*t61*t92+3.722102945396400
     ?         5d+0*t114*t13*t154*t2*t3*t92-1.8610514726982003d+0*t11*t1
     @         54*t3*t55*t61*t91+1.8610514726982003d+0*t107*t108*t121*t3
     1         *t60*t61+4.9628039271952007d+0*gammaaa*t11*t227*t3*t4*t55
     2         *t61+7.444205890792801d+0*t13*t153*t3*t4*t61+3.7221029453
     3         964005d+0*t114*t121*t13*t3*t4*t60-5.583154418094601d+0*t1
     4         3*t175*t2*t201*t3*t60+2.4814019635976003d+0*t114*t13*t175
     5         *t3*t57-4.9628039271952007d+0*gammaaa*t11*t114*t175*t3*t5
     6         5*t56+1.4888411781585603d+1*gammaaa*t10*t193*t195*t227*t3
     7         -7.444205890792801d+0*t10*t107*t108*t153*t3)+t181*t72+t18
     8         2*t71)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t53*(3.7221029453964005d+0*t11*t1
     1         37*t3*t73*t78*t98*t99+2.481401963597601d+0*t11*t3*t73*t74
     2         *t79*t99-9.925607854390401d+0*gammabb*t107*t131*t163*t3*t
     3         79*t99+3.7221029453964005d+0*t137*t164*t3*t33*t40*t99-1.8
     4         610514726982003d+0*t11*t164*t3*t73*t79*t98+1.861051472698
     5         2003d+0*t107*t131*t144*t3*t78*t79+4.9628039271952007d+0*g
     6         ammabb*t11*t230*t3*t4*t73*t79+7.444205890792801d+0*t163*t
     7         3*t4*t40*t79+3.7221029453964005d+0*t137*t144*t3*t4*t40*t7
     8         8-5.583154418094601d+0*t185*t214*t3*t33*t40*t78+2.4814019
     9         635976003d+0*t137*t185*t3*t40*t75-4.9628039271952007d+0*g
     :         ammabb*t11*t137*t185*t3*t73*t74+1.4888411781585603d+1*gam
     ;         mabb*t193*t208*t230*t3*t38-7.444205890792801d+0*t107*t131
     <         *t163*t3*t38)+t190*t90+t191*t89+t41*(8.865384008920708d-2
     =         *t16*t17*t183*t42*t81*t88-1.7730768017841417d-1*t16*t17*t
     >         186*t42*t51*t88-3.5461536035682834d-1*t16*t17*t42*t44*(-6
     ?         .649038006690531d-2*t17*t183*t42*t81*t87+1.32980760133810
     @         6d-1*t17*t186*t42*t51*t87+2.659615202676212d-1*t17*t42*t4
     1         4*(1.413716694115413d+1*t186*t24*t45*t48*t85*t86-5.0d-1*t
     2         186*t47*t48*t86-3.536776513153214d-2*t186*t21*t46*t48*t86
     3         -3.534291735288533d+1*t145*t183*t24*t45*t48*t86+1.9985948
     4         912206125d+2*t143*t177*t183*t187*t48*t86-5.0d-1*t186*t48*
     5         t80*t85+2.827433388230826d+1*t100*t167*t24*t45*t48*t85-1.
     6         0d+0*t100*t165*t48*t85+7.073553026306428d-2*t186*t21*t49*
     7         t84-7.068583470577066d+0*t183*t187*t24*t45*t48*t80+1.0d+0
     8         *t145*t183*t48*t80-3.536776513153214d-2*t21*t231*t46*t49+
     9         1.0d+0*t186*t42*t47*t48-1.0d+0*t100*t167*t47*t48-3.536776
     :         513153214d-2*t21*t232*t43*t46*t48-7.073553026306428d-2*t1
     ;         00*t167*t21*t46*t48+1.413716694115413d+1*rhob*t145*t183*t
     <         24*t48+t232*t48)+5.01149519714857d+1*t186*t24*t45*t48*t54
     =         *t82*t85+1.002299039429714d+2*t100*t166*t24*t45*t48*t54*t
     >         85-1.002299039429714d+2*t145*t183*t24*t45*t48*t54*t82+7.0
     ?         84834422688148d+2*t143*t177*t183*t187*t48*t54*t82+3.54490
     @         77018110403d+0*t48*t54*(-1.8799712059732546d+0*rhob*t231*
     1         t29*t81-1.8799712059732546d+0*t186*t29*t81+2.819956808959
     2         8814d+0*rhob*t130*t186*t29*t80-7.049892022399704d+0*rhob*
     3         t183*t184*t29*t80+2.8199568089598814d+0*t130*t183*t29+5.6
     4         39913617919763d+0*rhob*t100*t130*t165*t29)-6.649038006690
     5         531d-2*t17*t186*t42*t50*t80*t81-1.329807601338106d-1*t100
     6         *t101*t17*t42*t80*t81+6.649038006690531d-2*t17*t183*t46*t
     7         50*t81-1.329807601338106d-1*t100*t165*t17*t42*t50*t81+1.3
     8         29807601338106d-1*t17*t188*t42*t51*t80+9.973557010035795d
     9         -2*t130*t17*t183*t42*t50*t80-1.329807601338106d-1*t17*t18
     :         6*t46*t50*t51+1.329807601338106d-1*t17*t231*t42*t50*t51-2
     ;         .659615202676212d-1*t100*t101*t17*t46*t51+2.6596152026762
     <         12d-1*t100*t168*t17*t42*t51+2.659615202676212d-1*t101*t16
     =         5*t17*t42*t51-2.659615202676212d-1*t17*t188*t44*t46)+8.86
     >         5384008920708d-2*t16*t17*t186*t42*t52*t80*t81+1.773076801
     ?         7841417d-1*t100*t102*t16*t17*t42*t80*t81-8.86538400892070
     @         8d-2*t16*t17*t183*t46*t52*t81+1.7730768017841417d-1*t100*
     1         t16*t165*t17*t42*t52*t81-1.3298076013381063d-1*t130*t16*t
     2         17*t183*t42*t52*t80-1.7730768017841417d-1*t16*t17*t189*t4
     3         2*t51*t80+1.7730768017841417d-1*t16*t17*t186*t46*t51*t52-
     4         1.7730768017841417d-1*t16*t17*t231*t42*t51*t52+3.54615360
     5         35682834d-1*t100*t102*t16*t17*t46*t51-3.5461536035682834d
     6         -1*t100*t16*t169*t17*t42*t51-3.5461536035682834d-1*t102*t
     7         16*t165*t17*t42*t51+3.5461536035682834d-1*t16*t17*t189*t4
     8         4*t46)+2*t103*t171+2*t104*t170)*wght+Cmat3(iq,D3_RB_GBB_G
     9         BB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (3*t181*t97+3*t182*t96+t14*(2.659
     1         6152026762127d-1*t16*t17*t173*t18*t63*t95-5.3192304053524
     2         25d-1*t16*t17*t176*t18*t30*t95-3.5461536035682834d-1*t16*
     3         t17*t18*t20*(-1.994711402007159d-1*t17*t173*t18*t63*t94+3
     4         .989422804014318d-1*t17*t176*t18*t30*t94+2.65961520267621
     5         2d-1*t17*t18*t20*(-1.0602875205865599d+2*t122*t176*t22*t2
     6         4*t26*t93+5.995784673661838d+2*t120*t176*t177*t178*t26*t9
     7         3+1.413716694115413d+1*t22*t237*t24*t26*t67-3.53677651315
     8         3214d-2*t21*t23*t237*t27-5.0d-1*t237*t25*t26+1.0602875205
     9         865599d+2*t178*t22*t234*t24*t26+2.8254469624923584d+3*t10
     :         9*t234*t238*t26/t19**6-1.2990866792933983d+3*t120*t177*t2
     ;         34*t26/t19**5)-2.826440000853771d+2*t172*t174*t176*t26*t5
     <         4*t65*t93-1.994711402007159d-1*t17*t176*t18*t28*t63*t93+2
     =         .998945983256911d+1*rhoa*t106*t176*t26*t29*t54*t93+3.9894
     >         22804014318d-1*t17*t179*t18*t30*t93+4.710733334756285d+2*
     ?         t172*t234*t26*t54*t65/t20**9-6.66432440723758d+0*rhoa*t23
     @         7*t26*t29*t54*t63-2.4991216527140922d+1*rhoa*t174*t234*t2
     1         6*t29*t54-1.3319284713741864d+3*t152*t233*t234*t26*t54/t2
     2         0**11+1.329807601338106d-1*t17*t18*t237*t28*t30+9.9735570
     3         10035795d-2*t106*t17*t18*t234*t28)+2.6596152026762127d-1*
     4         t16*t17*t176*t18*t31*t63*t93-5.319230405352425d-1*t16*t17
     5         *t18*t180*t30*t93-1.7730768017841417d-1*t16*t17*t18*t237*
     6         t30*t31-1.3298076013381063d-1*t106*t16*t17*t18*t234*t31)+
     7         t32*(5.583154418094601d+0*t107*t108*t121*t3*t61*t92+1.116
     8         6308836189202d+1*t114*t121*t13*t3*t4*t92+5.58315441809460
     9         1d+0*t11*t114*t175*t3*t55*t91-5.583154418094601d+0*t11*t2
     :         35*t3*t4*t55*t61-5.583154418094601d+0*t13*t2*t201*t236*t3
     ;         -5.583154418094601d+0*t10*t193*t195*t235*t3))*wght+Cmat3(
     <         iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (t53*(5.583154418094601d+0*t107*t
     1         131*t144*t3*t79*t99+1.1166308836189202d+1*t137*t144*t3*t4
     2         *t40*t99+5.583154418094601d+0*t11*t137*t185*t3*t73*t98-5.
     3         583154418094601d+0*t11*t240*t3*t4*t73*t79-5.5831544180946
     4         01d+0*t214*t241*t3*t33*t40-5.583154418094601d+0*t193*t208
     5         *t240*t3*t38)+t41*(-3.5461536035682834d-1*t16*t17*t42*t44
     6         *(2.659615202676212d-1*t17*t42*t44*(1.413716694115413d+1*
     7         t24*t242*t45*t48*t85-3.536776513153214d-2*t21*t242*t46*t4
     8         9-5.0d-1*t242*t47*t48+1.0602875205865599d+2*t187*t239*t24
     9         *t45*t48-1.0602875205865599d+2*t100*t145*t186*t24*t45*t48
     :         -1.2990866792933983d+3*t143*t177*t239*t48/t43**5+2.825446
     ;         9624923584d+3*t132*t238*t239*t48/t43**6+5.995784673661838
     <         d+2*t100*t143*t177*t186*t187*t48)+4.710733334756285d+2*t1
     =         72*t239*t48*t54*t83/t44**9-2.826440000853771d+2*t100*t172
     >         *t184*t186*t48*t54*t83-6.66432440723758d+0*rhob*t242*t29*
     ?         t48*t54*t81-1.994711402007159d-1*t100*t17*t186*t42*t50*t8
     @         1-1.994711402007159d-1*t101*t17*t183*t42*t81-1.3319284713
     1         741864d+3*t162*t233*t239*t48*t54/t44**11-2.49912165271409
     2         22d+1*rhob*t184*t239*t29*t48*t54+2.998945983256911d+1*rho
     3         b*t100*t130*t186*t29*t48*t54+1.329807601338106d-1*t17*t24
     4         2*t42*t50*t51+3.989422804014318d-1*t100*t17*t188*t42*t51+
     5         3.989422804014318d-1*t101*t17*t186*t42*t51+9.973557010035
     6         795d-2*t130*t17*t239*t42*t50)+2.6596152026762127d-1*t100*
     7         t16*t17*t186*t42*t52*t81+2.6596152026762127d-1*t102*t16*t
     8         17*t183*t42*t81-1.7730768017841417d-1*t16*t17*t242*t42*t5
     9         1*t52-1.3298076013381063d-1*t130*t16*t17*t239*t42*t52-5.3
     :         19230405352425d-1*t100*t16*t17*t189*t42*t51-5.31923040535
     ;         2425d-1*t102*t16*t17*t186*t42*t51)+3*t103*t191+3*t104*t19
     <         0)*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = rhoa**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammaaa**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = gammaaa*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammaaa*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhoa
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhoa**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhoa*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = rhoa**3.333333333333333d-1
            t36 = 1/rhoa**6.333333333333333d+0
            t37 = 1/rhoa**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammaaa*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         a*t28*t40*t41
            t43 = rhoa**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhoa*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammaaa*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammaaa*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhoa*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = rhoa**6
            t63 = 1/t62
            t64 = 1/rhoa**3.3333333333333337d+0
            t65 = 1/rhoa**6.666666666666666d-1
            t66 = t38**2
            t67 = 1/t9**3
            t68 = 1/rhoa**7.333333333333333d+0
            t69 = 1/rhoa**4.666666666666667d+0
            t70 = 9.777777777777779d+0*gammaaa*t69*t7+3.377777777777778d
     1         +1*t4*t5*t68
            t71 = 9.305257363491002d-1*t13*t2*t3*t39*t70-1.8610514726982
     1         003d+0*t13*t2*t3*t66*t67+4.135669939329334d-1*t10*t13*t3*
     2         t65+4.135669939329334d-1*t1*t65+2.481401963597601d+0*gamm
     3         aaa*t10*t11*t3*t33*t64-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t63+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t72 = -1.8799712059732546d+0*rhoa*t28*t41*t71+2.819956808959
     1         8814d+0*rhoa*t28*t58*t59-3.759942411946509d+0*t28*t40*t41
            t73 = rhoa**4
            t74 = 1/t73
            t75 = 1/t18**3
            t76 = -2.827433388230826d+1*t21*t23*t58*t75+1.41371669411541
     1         3d+1*t21*t23*t45*t71+5.654866776461652d+1*rhoa*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t77 = t46**2
            t78 = -3.536776513153214d-2*t18*t20*t22*t25*t77+t25*t77-3.53
     1         6776513153214d-2*t18*t20*t22*t25*t76+t25*t76-2.1220659078
     2         919285d-1*t18*t20*t26*t74-3.536776513153214d-2*t20*t22*t2
     3         6*t71+1.4147106052612857d-1*t18*t20*t25*t44*t46-7.0735530
     4         26306428d-2*t20*t22*t25*t40*t46+1.4147106052612857d-1*t20
     5         *t26*t40*t44
            t79 = 2.659615202676212d-1*t16*t17*t19*t78+3.544907701811040
     1         3d+0*t25*t32*t72+1.329807601338106d-1*t16*t17*t27*t29*t71
     2         -6.649038006690531d-2*t16*t17*t27*t41*t58+2.6596152026762
     3         12d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t16*t19*t2
     4         2*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.31923040535
     5         2424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*t22*t27*
     6         t29*t40
            t80 = -3.5461536035682834d-1*t15*t16*t17*t19*t79-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t71+8.865384008920708d-2*t
     2         15*t16*t17*t30*t41*t58-3.5461536035682834d-1*t15*t16*t17*
     3         t29*t40*t48+7.092307207136567d-1*t15*t16*t19*t22*t48-7.09
     4         2307207136567d-1*t15*t16*t19*t30*t44+3.5461536035682834d-
     5         1*t15*t16*t22*t29*t30*t40
            t81 = -9.305257363491002d-1*t13*t2*t3*t39*t70+1.861051472698
     1         2003d+0*t13*t2*t3*t66*t67-4.135669939329334d-1*t10*t13*t3
     2         *t65-4.135669939329334d-1*t1*t65-2.481401963597601d+0*gam
     3         maaa*t10*t11*t3*t33*t64+1.3234143805853868d+1*t10*t3*t5*t
     4         60*t61*t63-2.4814019635976003d+0*t13*t3*t35*t38*t39+4.962
     5         8039271952007d+0*gammaaa*t11*t3*t33*t34*t38*t39
            t82 = rhoa**5
            t83 = 1/t82
            t84 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmaaa*t36*t4
            t85 = 9.305257363491002d-1*t13*t2*t3*t39*t84+4.9628039271952
     1         007d+0*gammaaa*t10*t3*t60*t61*t83-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t67+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t86 = -1.8799712059732546d+0*rhoa*t28*t41*t85+2.819956808959
     1         8814d+0*rhoa*t28*t40*t53*t59-1.8799712059732546d+0*t28*t4
     2         1*t53
            t87 = 1.413716694115413d+1*t21*t23*t45*t85-2.827433388230826
     1         d+1*t21*t23*t40*t53*t75+2.827433388230826d+1*rhoa*t23*t45
     2         *t53
            t88 = -3.536776513153214d-2*t18*t20*t22*t25*t87+t25*t87-3.53
     1         6776513153214d-2*t20*t22*t26*t85+1.413716694115413d+1*t21
     2         *t23*t25*t45*t46*t53-5.0d-1*t24*t25*t46*t53-3.53677651315
     3         3214d-2*t20*t22*t25*t46*t53-5.0d-1*t25*t40*t45*t53+7.0735
     4         53026306428d-2*t20*t26*t44*t53+1.0d+0*t17*t24*t25*t53
            t89 = 2.659615202676212d-1*t16*t17*t19*t88+3.544907701811040
     1         3d+0*t25*t32*t86+1.329807601338106d-1*t16*t17*t27*t29*t85
     2         +1.329807601338106d-1*t16*t17*t29*t40*t54-2.6596152026762
     3         12d-1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t4
     4         7*t53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.6
     5         49038006690531d-2*t16*t17*t27*t40*t41*t53-1.3298076013381
     6         06d-1*t16*t22*t27*t29*t53
            t90 = -3.5461536035682834d-1*t15*t16*t17*t19*t89-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t85-1.7730768017841417d-1*
     2         t15*t16*t17*t29*t40*t55+3.5461536035682834d-1*t15*t16*t19
     3         *t22*t55-1.7730768017841417d-1*t15*t16*t17*t29*t48*t53+8.
     4         865384008920708d-2*t15*t16*t17*t30*t40*t41*t53+1.77307680
     5         17841417d-1*t15*t16*t22*t29*t30*t53
            t91 = -9.305257363491002d-1*t13*t2*t3*t39*t84-4.962803927195
     1         2007d+0*gammaaa*t10*t3*t60*t61*t83+1.8610514726982003d+0*
     2         t13*t2*t3*t38*t52*t67-1.2407009817988002d+0*t13*t3*t35*t3
     3         9*t52+2.4814019635976003d+0*gammaaa*t11*t3*t33*t34*t39*t5
     4         2-9.305257363491002d-1*t11*t3*t33*t38*t39*t51+1.240700981
     5         7988002d+0*t10*t11*t3*t33*t34
            t92 = 1/t16**3
            t93 = t53**2
            t94 = 1/t19**7
            t95 = t52**2
            t96 = -1.8610514726982003d+0*t13*t2*t3*t67*t95-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t97 = 1/t16**4
            t98 = 1/t18**4
            t99 = 1.9985948912206125d+2*t25*t73*t93*t97*t98+1.4137166941
     1         15413d+1*t21*t23*t25*t45*t96-3.536776513153214d-2*t20*t22
     2         *t26*t96-5.0d-1*t24*t25*t96-3.534291735288533d+1*t21*t23*
     3         t25*t75*t93
            t100 = 2.659615202676212d-1*t16*t17*t19*t99-6.66432440723758
     1         d+0*rhoa*t25*t28*t32*t41*t96+1.329807601338106d-1*t16*t17
     2         *t27*t29*t96-9.421466669512572d+1*t25*t32*t43*t92*t93*t94
     3         +9.996486610856369d+0*rhoa*t25*t28*t32*t59*t93-6.64903800
     4         6690531d-2*t16*t17*t27*t41*t93+2.659615202676212d-1*t16*t
     5         17*t29*t53*t54
            t101 = -1.7730768017841417d-1*t15*t16*t17*t29*t30*t96+8.8653
     1         84008920708d-2*t15*t16*t17*t30*t41*t93-3.5461536035682834
     2         d-1*t15*t16*t17*t29*t53*t55-3.5461536035682834d-1*t100*t1
     3         5*t16*t17*t19
            t102 = 1.8610514726982003d+0*t13*t2*t3*t67*t95+1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74-1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74-1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t103 = t40**3
            t104 = t11**3
            t105 = gammaaa**3
            t106 = 1/t12**4
            t107 = 1/rhoa**9.666666666666666d+0
            t108 = 1/rhoa**7
            t109 = 1/rhoa**4.333333333333333d+0
            t110 = 1/rhoa**1.6666666666666669d+0
            t111 = t38**3
            t112 = 1/t9**4
            t113 = -4.5629629629629626d+1*gammaaa*t7/rhoa**5.66666666666
     1         6667d+0-2.477037037037037d+2*t4*t5/rhoa**8.33333333333333
     2         4d+0
            t114 = -5.583154418094601d+0*t13*t2*t3*t38*t67*t70+3.7221029
     1         453964005d+0*t13*t3*t35*t39*t70-7.444205890792801d+0*gamm
     2         aaa*t11*t3*t33*t34*t39*t70-7.444205890792801d+0*t13*t3*t3
     3         5*t66*t67+1.4888411781585603d+1*gammaaa*t11*t3*t33*t34*t6
     4         6*t67+1.2407009817988002d+0*t13*t3*t38*t39*t65+7.44420589
     5         0792801d+0*gammaaa*t11*t3*t33*t38*t39*t64-3.9702431417561
     6         607d+1*t3*t38*t39*t5*t60*t61*t63+9.263900664097708d+1*t10
     7         *t108*t3*t5*t60*t61+9.305257363491002d-1*t113*t13*t2*t3*t
     8         39-9.374185195813161d+0*gammaaa*t10*t109*t11*t3*t33+5.583
     9         154418094601d+0*t111*t112*t13*t2*t3-2.7571132928862224d-1
     :         *t10*t110*t13*t3-1.0587315044683095d+2*t10*t104*t105*t106
     ;         *t107*t3-2.7571132928862224d-1*t1*t110
            t115 = 8.482300164692478d+1*t103*t21*t23*t98-8.4823001646924
     1         78d+1*t21*t23*t40*t71*t75-1.6964600329384957d+2*rhoa*t23*
     2         t58*t75+8.482300164692478d+1*rhoa*t23*t45*t71+8.482300164
     3         692478d+1*t23*t40*t45+1.413716694115413d+1*t114*t21*t23*t
     4         45
            t116 = t46**3
            t117 = 1/rhoa**8.666666666666666d+0
            t118 = 9.777777777777779d+0*t69*t7+6.755555555555556d+1*gamm
     1         aaa*t4*t68
            t119 = -3.7221029453964005d+0*t13*t2*t3*t38*t67*t84+2.481401
     1         9635976003d+0*t13*t3*t35*t39*t84-4.9628039271952007d+0*ga
     2         mmaaa*t11*t3*t33*t34*t39*t84+9.925607854390401d+0*gammaaa
     3         *t3*t38*t39*t60*t61*t83-1.8610514726982003d+0*t13*t2*t3*t
     4         52*t67*t70+9.305257363491002d-1*t11*t3*t33*t39*t51*t70-1.
     5         8610514726982003d+0*t11*t3*t33*t51*t66*t67-4.962803927195
     6         2007d+0*t13*t3*t35*t38*t52*t67+9.925607854390401d+0*gamma
     7         aa*t11*t3*t33*t34*t38*t52*t67+5.583154418094601d+0*t112*t
     8         13*t2*t3*t52*t66+4.135669939329334d-1*t13*t3*t39*t52*t65+
     9         2.481401963597601d+0*gammaaa*t11*t3*t33*t39*t52*t64+2.894
     :         9689575305343d+0*t10*t11*t3*t33*t64-1.3234143805853868d+1
     ;         *t3*t39*t5*t52*t60*t61*t63-3.143109153890294d+1*gammaaa*t
     <         10*t3*t60*t61*t63+3.9702431417561607d+1*t10*t104*t106*t11
     =         7*t3*t5-2.4814019635976003d+0*t11*t3*t33*t34*t38*t39+9.30
     >         5257363491002d-1*t118*t13*t2*t3*t39
            t120 = 8.482300164692478d+1*t21*t23*t53*t58*t98-5.6548667764
     1         61652d+1*t21*t23*t40*t75*t85+5.654866776461652d+1*rhoa*t2
     2         3*t45*t85-2.827433388230826d+1*t21*t23*t53*t71*t75-1.1309
     3         733552923305d+2*rhoa*t23*t40*t53*t75+2.827433388230826d+1
     4         *t23*t45*t53+1.413716694115413d+1*t119*t21*t23*t45
            t121 = 1/rhoa**7.666666666666667d+0
            t122 = -2.4814019635976003d+0*t13*t3*t35*t67*t95+4.962803927
     1         1952007d+0*gammaaa*t11*t3*t33*t34*t67*t95+5.5831544180946
     2         01d+0*t112*t13*t2*t3*t38*t95-3.7221029453964005d+0*t13*t2
     3         *t3*t52*t67*t84+1.8610514726982003d+0*t11*t3*t33*t39*t51*
     4         t84+9.925607854390401d+0*gammaaa*t3*t39*t52*t60*t61*t83+7
     5         .444205890792801d+0*t10*t3*t60*t61*t83-7.444205890792801d
     6         +0*t13*t3*t39*t4*t83-3.7221029453964005d+0*t13*t3*t38*t4*
     7         t67*t74-1.8610514726982003d+0*t3*t38*t39*t60*t61*t74-3.72
     8         21029453964005d+0*t11*t3*t33*t38*t51*t52*t67-2.4814019635
     9         97601d+0*t11*t3*t33*t34*t39*t52-4.9628039271952007d+0*gam
     :         maaa*t11*t121*t3*t33*t39*t4-1.4888411781585603d+1*gammaaa
     ;         *t10*t104*t106*t121*t3
            t123 = 8.482300164692478d+1*t21*t23*t40*t93*t98-2.8274333882
     1         30826d+1*t21*t23*t40*t75*t96+2.827433388230826d+1*rhoa*t2
     2         3*t45*t96-5.654866776461652d+1*rhoa*t23*t75*t93-5.6548667
     3         76461652d+1*t21*t23*t53*t75*t85+1.413716694115413d+1*t122
     4         *t21*t23*t45
            t124 = t53**3
            t125 = 1/rhoa**6.666666666666667d+0
            t126 = t52**3
            t127 = -5.583154418094601d+0*t11*t3*t33*t51*t67*t95-1.116630
     1         8836189202d+1*t13*t3*t4*t52*t67*t74-5.583154418094601d+0*
     2         t3*t39*t52*t60*t61*t74+5.583154418094601d+0*t11*t125*t3*t
     3         33*t39*t4+5.583154418094601d+0*t112*t126*t13*t2*t3+5.5831
     4         54418094601d+0*t10*t104*t106*t125*t3
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = t31*t50*wght+t14*t49*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t31*t57*wght+t14*t56*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t31*t81*wght+t14*t80*wght+2*t49*t50*wgh
     1         t+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t31*t91*wght+t14*t90*wght+t49*t57*wght
     1         +t50*t56*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 2*t56*t57*wght+t102*t31*wght+t101*t14
     1         *wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t14*(-3.5461536035682834d-1*t15*t16*
     1         t17*t19*(3.5449077018110403d+0*t25*t32*(-7.04989202239970
     2         4d+0*rhoa*t103*t28*t94+8.459870426879645d+0*rhoa*t28*t40*
     3         t59*t71-5.639913617919763d+0*t28*t41*t71+8.45987042687964
     4         5d+0*t28*t58*t59-1.8799712059732546d+0*rhoa*t114*t28*t41)
     5         +2.659615202676212d-1*t16*t17*t19*(8.488263631567714d-1*t
     6         18*t20*t26*t83+2.1220659078919285d-1*t18*t20*t25*t44*t77-
     7         1.0610329539459643d-1*t20*t22*t25*t40*t77-1.0610329539459
     8         643d-1*t18*t20*t22*t25*t46*t76+3*t25*t46*t76+2.1220659078
     9         919285d-1*t18*t20*t25*t44*t76-1.0610329539459643d-1*t20*t
     :         22*t25*t40*t76-6.366197723675786d-1*t18*t20*t25*t46*t74-6
     ;         .366197723675786d-1*t20*t26*t40*t74-1.0610329539459643d-1
     <         *t20*t22*t25*t46*t71+2.1220659078919285d-1*t20*t26*t44*t7
     =         1+4.244131815783857d-1*t20*t25*t40*t44*t46-3.536776513153
     >         214d-2*t114*t20*t22*t26-3.536776513153214d-2*t116*t18*t20
     ?         *t22*t25-3.536776513153214d-2*t115*t18*t20*t22*t25+t116*t
     @         25+t115*t25)+3.989422804014318d-1*t16*t17*t29*t40*t78-7.9
     1         78845608028636d-1*t16*t19*t22*t78+3.5449077018110403d+0*t
     2         25*t32*t42*t77+3.5449077018110403d+0*t25*t32*t42*t76-1.59
     3         57691216057274d+0*t16*t19*t27*t74+7.089815403622081d+0*t2
     4         5*t32*t46*t72+3.989422804014318d-1*t16*t17*t29*t47*t71-1.
     5         994711402007159d-1*t16*t17*t27*t40*t41*t71-3.989422804014
     6         318d-1*t16*t22*t27*t29*t71+9.973557010035795d-2*t103*t16*
     7         t17*t27*t59-1.994711402007159d-1*t16*t17*t41*t47*t58+1.99
     8         4711402007159d-1*t16*t22*t27*t41*t58+1.5957691216057274d+
     9         0*t16*t19*t44*t47-7.978845608028636d-1*t16*t22*t29*t40*t4
     :         7+7.978845608028636d-1*t16*t27*t29*t40*t44+1.329807601338
     ;         106d-1*t114*t16*t17*t27*t29)-5.319230405352425d-1*t15*t16
     <         *t17*t29*t40*t79+1.063846081070485d+0*t15*t16*t19*t22*t79
     =         +2.12769216214097d+0*t15*t16*t19*t30*t74-5.31923040535242
     >         5d-1*t15*t16*t17*t29*t48*t71+2.6596152026762127d-1*t15*t1
     ?         6*t17*t30*t40*t41*t71+5.319230405352425d-1*t15*t16*t22*t2
     @         9*t30*t71-1.3298076013381063d-1*t103*t15*t16*t17*t30*t59+
     1         2.6596152026762127d-1*t15*t16*t17*t41*t48*t58-2.659615202
     2         6762127d-1*t15*t16*t22*t30*t41*t58-2.12769216214097d+0*t1
     3         5*t16*t19*t44*t48+1.063846081070485d+0*t15*t16*t22*t29*t4
     4         0*t48-1.063846081070485d+0*t15*t16*t29*t30*t40*t44-1.7730
     5         768017841417d-1*t114*t15*t16*t17*t29*t30)*wght+3*t49*t81*
     6         wght+3*t50*t80*wght+t31*(5.583154418094601d+0*t13*t2*t3*t
     7         38*t67*t70-3.7221029453964005d+0*t13*t3*t35*t39*t70+7.444
     8         205890792801d+0*gammaaa*t11*t3*t33*t34*t39*t70+7.44420589
     9         0792801d+0*t13*t3*t35*t66*t67-1.4888411781585603d+1*gamma
     :         aa*t11*t3*t33*t34*t66*t67-1.2407009817988002d+0*t13*t3*t3
     ;         8*t39*t65-7.444205890792801d+0*gammaaa*t11*t3*t33*t38*t39
     <         *t64+3.9702431417561607d+1*t3*t38*t39*t5*t60*t61*t63-9.26
     =         3900664097708d+1*t10*t108*t3*t5*t60*t61-9.305257363491002
     >         d-1*t113*t13*t2*t3*t39+9.374185195813161d+0*gammaaa*t10*t
     ?         109*t11*t3*t33-5.583154418094601d+0*t111*t112*t13*t2*t3+2
     @         .7571132928862224d-1*t10*t110*t13*t3+1.0587315044683095d+
     1         2*t10*t104*t105*t106*t107*t3+2.7571132928862224d-1*t1*t11
     2         0)*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t14*(-3.5461536035682834d-1*t15*t16
     1         *t17*t19*(3.5449077018110403d+0*t25*t32*(-7.0498920223997
     2         04d+0*rhoa*t28*t53*t58*t94+5.639913617919763d+0*rhoa*t28*
     3         t40*t59*t85-3.759942411946509d+0*t28*t41*t85+2.8199568089
     4         598814d+0*rhoa*t28*t53*t59*t71+5.639913617919763d+0*t28*t
     5         40*t53*t59-1.8799712059732546d+0*rhoa*t119*t28*t41)+2.659
     6         615202676212d-1*t16*t17*t29*t40*t88-5.319230405352424d-1*
     7         t16*t19*t22*t88+2.659615202676212d-1*t16*t17*t19*(-7.0735
     8         53026306428d-2*t18*t20*t22*t25*t46*t87+2*t25*t46*t87+1.41
     9         47106052612857d-1*t18*t20*t25*t44*t87-7.073553026306428d-
     :         2*t20*t22*t25*t40*t87-7.073553026306428d-2*t20*t22*t25*t4
     ;         6*t85+1.4147106052612857d-1*t20*t26*t44*t85+1.41371669411
     <         5413d+1*t21*t23*t25*t45*t53*t77-5.0d-1*t24*t25*t53*t77-3.
     =         536776513153214d-2*t20*t22*t25*t53*t77+1.413716694115413d
     >         +1*t21*t23*t25*t45*t53*t76-5.0d-1*t24*t25*t53*t76-3.53677
     ?         6513153214d-2*t20*t22*t25*t53*t76-2.1220659078919285d-1*t
     @         20*t26*t53*t74-5.0d-1*t25*t45*t53*t71-1.0d+0*t25*t40*t45*
     1         t46*t53+1.4147106052612857d-1*t20*t25*t44*t46*t53+2.0d+0*
     2         t17*t24*t25*t46*t53+2.0d+0*t17*t25*t40*t45*t53-3.0d+0*t22
     3         *t24*t25*t53-3.536776513153214d-2*t119*t20*t22*t26-3.5367
     4         76513153214d-2*t120*t18*t20*t22*t25+t120*t25)+3.544907701
     5         8110403d+0*t25*t32*t42*t87+3.5449077018110403d+0*t25*t32*
     6         t46*t86+2.659615202676212d-1*t16*t17*t29*t47*t85-1.329807
     7         601338106d-1*t16*t17*t27*t40*t41*t85-2.659615202676212d-1
     8         *t16*t22*t27*t29*t85+1.329807601338106d-1*t16*t17*t29*t53
     9         *t78+5.01149519714857d+1*t21*t23*t25*t32*t45*t53*t72+1.32
     :         9807601338106d-1*t16*t17*t29*t54*t71-6.649038006690531d-2
     ;         *t16*t17*t27*t41*t53*t71+9.973557010035795d-2*t16*t17*t27
     <         *t53*t58*t59-6.649038006690531d-2*t16*t17*t41*t54*t58+5.3
     =         19230405352424d-1*t16*t19*t44*t54-2.659615202676212d-1*t1
     >         6*t22*t29*t40*t54-1.329807601338106d-1*t16*t17*t40*t41*t4
     ?         7*t53-2.659615202676212d-1*t16*t22*t29*t47*t53+5.01149519
     @         714857d+1*t21*t23*t25*t32*t42*t45*t46*t53+2.6596152026762
     1         12d-1*t16*t27*t29*t44*t53+1.329807601338106d-1*t16*t22*t2
     2         7*t40*t41*t53+1.329807601338106d-1*t119*t16*t17*t27*t29)-
     3         3.5461536035682834d-1*t15*t16*t17*t29*t40*t89+7.092307207
     4         136567d-1*t15*t16*t19*t22*t89-3.5461536035682834d-1*t15*t
     5         16*t17*t29*t48*t85+1.7730768017841417d-1*t15*t16*t17*t30*
     6         t40*t41*t85+3.5461536035682834d-1*t15*t16*t22*t29*t30*t85
     7         -1.7730768017841417d-1*t15*t16*t17*t29*t53*t79-1.77307680
     8         17841417d-1*t15*t16*t17*t29*t55*t71+8.865384008920708d-2*
     9         t15*t16*t17*t30*t41*t53*t71-1.3298076013381063d-1*t15*t16
     :         *t17*t30*t53*t58*t59+8.865384008920708d-2*t15*t16*t17*t41
     ;         *t55*t58-7.092307207136567d-1*t15*t16*t19*t44*t55+3.54615
     <         36035682834d-1*t15*t16*t22*t29*t40*t55+1.7730768017841417
     =         d-1*t15*t16*t17*t40*t41*t48*t53+3.5461536035682834d-1*t15
     >         *t16*t22*t29*t48*t53-3.5461536035682834d-1*t15*t16*t29*t3
     ?         0*t44*t53-1.7730768017841417d-1*t15*t16*t22*t30*t40*t41*t
     @         53-1.7730768017841417d-1*t119*t15*t16*t17*t29*t30)*wght+2
     1         *t49*t91*wght+2*t50*t90*wght+t31*(3.7221029453964005d+0*t
     2         13*t2*t3*t38*t67*t84-2.4814019635976003d+0*t13*t3*t35*t39
     3         *t84+4.9628039271952007d+0*gammaaa*t11*t3*t33*t34*t39*t84
     4         -9.925607854390401d+0*gammaaa*t3*t38*t39*t60*t61*t83+1.86
     5         10514726982003d+0*t13*t2*t3*t52*t67*t70-9.305257363491002
     6         d-1*t11*t3*t33*t39*t51*t70+1.8610514726982003d+0*t11*t3*t
     7         33*t51*t66*t67+4.9628039271952007d+0*t13*t3*t35*t38*t52*t
     8         67-9.925607854390401d+0*gammaaa*t11*t3*t33*t34*t38*t52*t6
     9         7-5.583154418094601d+0*t112*t13*t2*t3*t52*t66-4.135669939
     :         329334d-1*t13*t3*t39*t52*t65-2.481401963597601d+0*gammaaa
     ;         *t11*t3*t33*t39*t52*t64-2.8949689575305343d+0*t10*t11*t3*
     <         t33*t64+1.3234143805853868d+1*t3*t39*t5*t52*t60*t61*t63+3
     =         .143109153890294d+1*gammaaa*t10*t3*t60*t61*t63-3.97024314
     >         17561607d+1*t10*t104*t106*t117*t3*t5+2.4814019635976003d+
     ?         0*t11*t3*t33*t34*t38*t39-9.305257363491002d-1*t118*t13*t2
     @         *t3*t39)*wght+t56*t81*wght+t57*t80*wght+Cmat3(iq,D3_RA_RA
     1         _GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t14*(-3.5461536035682834d-1*t15*t1
     1         6*t17*t19*(1.329807601338106d-1*t16*t17*t29*t40*t99-2.659
     2         615202676212d-1*t16*t19*t22*t99+2.659615202676212d-1*t16*
     3         t17*t19*(1.9985948912206125d+2*t25*t46*t73*t93*t97*t98-7.
     4         068583470577066d+0*t21*t23*t25*t40*t93*t98+1.413716694115
     5         413d+1*t21*t23*t25*t45*t46*t96-5.0d-1*t24*t25*t46*t96-3.5
     6         36776513153214d-2*t20*t22*t25*t46*t96-5.0d-1*t25*t40*t45*
     7         t96+7.073553026306428d-2*t20*t26*t44*t96+1.0d+0*t17*t24*t
     8         25*t96-3.534291735288533d+1*t21*t23*t25*t46*t75*t93+1.0d+
     9         0*t25*t40*t75*t93+1.413716694115413d+1*rhoa*t23*t25*t75*t
     :         93+2.827433388230826d+1*t21*t23*t25*t45*t53*t87-1.0d+0*t2
     ;         4*t25*t53*t87-7.073553026306428d-2*t20*t22*t25*t53*t87-1.
     <         0d+0*t25*t45*t53*t85-3.536776513153214d-2*t122*t20*t22*t2
     =         6-3.536776513153214d-2*t123*t18*t20*t22*t25+t123*t25)+7.0
     >         84834422688148d+2*t25*t32*t42*t73*t93*t97*t98+3.544907701
     ?         8110403d+0*t25*t32*(2.8199568089598814d+0*rhoa*t28*t40*t5
     @         9*t96-1.8799712059732546d+0*t28*t41*t96-7.049892022399704
     1         d+0*rhoa*t28*t40*t93*t94+2.8199568089598814d+0*t28*t59*t9
     2         3+5.639913617919763d+0*rhoa*t28*t53*t59*t85-1.87997120597
     3         32546d+0*rhoa*t122*t28*t41)+1.329807601338106d-1*t16*t17*
     4         t29*t47*t96+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t
     5         96-6.649038006690531d-2*t16*t17*t27*t40*t41*t96-1.3298076
     6         01338106d-1*t16*t22*t27*t29*t96-1.002299039429714d+2*t21*
     7         t23*t25*t32*t42*t75*t93+9.973557010035795d-2*t16*t17*t27*
     8         t40*t59*t93-6.649038006690531d-2*t16*t17*t41*t47*t93+6.64
     9         9038006690531d-2*t16*t22*t27*t41*t93+2.659615202676212d-1
     :         *t16*t17*t29*t53*t88+1.002299039429714d+2*t21*t23*t25*t32
     ;         *t45*t53*t86+2.659615202676212d-1*t16*t17*t29*t54*t85-1.3
     <         29807601338106d-1*t16*t17*t27*t41*t53*t85-1.3298076013381
     =         06d-1*t16*t17*t40*t41*t53*t54-2.659615202676212d-1*t16*t2
     >         2*t29*t53*t54+1.329807601338106d-1*t122*t16*t17*t27*t29)-
     ?         1.7730768017841417d-1*t15*t16*t17*t29*t48*t96+8.865384008
     @         920708d-2*t15*t16*t17*t30*t40*t41*t96+1.7730768017841417d
     1         -1*t15*t16*t22*t29*t30*t96-1.3298076013381063d-1*t15*t16*
     2         t17*t30*t40*t59*t93+8.865384008920708d-2*t15*t16*t17*t41*
     3         t48*t93-8.865384008920708d-2*t15*t16*t22*t30*t41*t93-3.54
     4         61536035682834d-1*t15*t16*t17*t29*t53*t89-3.5461536035682
     5         834d-1*t15*t16*t17*t29*t55*t85+1.7730768017841417d-1*t15*
     6         t16*t17*t30*t41*t53*t85+1.7730768017841417d-1*t15*t16*t17
     7         *t40*t41*t53*t55+3.5461536035682834d-1*t15*t16*t22*t29*t5
     8         3*t55-1.7730768017841417d-1*t100*t15*t16*t17*t29*t40-1.77
     9         30768017841417d-1*t122*t15*t16*t17*t29*t30+3.546153603568
     :         2834d-1*t100*t15*t16*t19*t22)*wght+t31*(2.481401963597600
     ;         3d+0*t13*t3*t35*t67*t95-4.9628039271952007d+0*gammaaa*t11
     <         *t3*t33*t34*t67*t95-5.583154418094601d+0*t112*t13*t2*t3*t
     =         38*t95+3.7221029453964005d+0*t13*t2*t3*t52*t67*t84-1.8610
     >         514726982003d+0*t11*t3*t33*t39*t51*t84-9.925607854390401d
     ?         +0*gammaaa*t3*t39*t52*t60*t61*t83-7.444205890792801d+0*t1
     @         0*t3*t60*t61*t83+7.444205890792801d+0*t13*t3*t39*t4*t83+3
     1         .7221029453964005d+0*t13*t3*t38*t4*t67*t74+1.861051472698
     2         2003d+0*t3*t38*t39*t60*t61*t74+3.7221029453964005d+0*t11*
     3         t3*t33*t38*t51*t52*t67+2.481401963597601d+0*t11*t3*t33*t3
     4         4*t39*t52+4.9628039271952007d+0*gammaaa*t11*t121*t3*t33*t
     5         39*t4+1.4888411781585603d+1*gammaaa*t10*t104*t106*t121*t3
     6         )*wght+2*t56*t91*wght+2*t57*t90*wght+t101*t50*wght+t102*t
     7         49*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t14*(-3.5461536035682834d-1*t15*t
     1         16*t17*t19*(3.989422804014318d-1*t16*t17*t29*t53*t99+2.65
     2         9615202676212d-1*t16*t17*t19*(5.995784673661838d+2*t25*t5
     3         3*t73*t96*t97*t98+1.0602875205865599d+2*t124*t21*t23*t25*
     4         t98-1.2990866792933983d+3*t124*t25*t73*t97/t18**5-1.06028
     5         75205865599d+2*t21*t23*t25*t53*t75*t96+2.8254469624923584
     6         d+3*t124*t25*t62/(t16**6*t18**6)+1.413716694115413d+1*t12
     7         7*t21*t23*t25*t45-3.536776513153214d-2*t127*t20*t22*t26-5
     8         .0d-1*t127*t24*t25)-2.826440000853771d+2*t25*t32*t43*t53*
     9         t92*t94*t96+2.998945983256911d+1*rhoa*t25*t28*t32*t53*t59
     :         *t96+3.989422804014318d-1*t16*t17*t29*t54*t96-1.994711402
     ;         007159d-1*t16*t17*t27*t41*t53*t96-2.4991216527140922d+1*r
     <         hoa*t124*t25*t28*t32*t94-1.994711402007159d-1*t16*t17*t41
     =         *t54*t93+4.710733334756285d+2*t124*t25*t32*t43*t92/t19**9
     >         -1.3319284713741864d+3*t124*t25*t32*t82/(t16**5*t19**11)+
     ?         9.973557010035795d-2*t124*t16*t17*t27*t59-6.6643244072375
     @         8d+0*rhoa*t127*t25*t28*t32*t41+1.329807601338106d-1*t127*
     1         t16*t17*t27*t29)-5.319230405352425d-1*t15*t16*t17*t29*t55
     2         *t96+2.6596152026762127d-1*t15*t16*t17*t30*t41*t53*t96+2.
     3         6596152026762127d-1*t15*t16*t17*t41*t55*t93-1.32980760133
     4         81063d-1*t124*t15*t16*t17*t30*t59-5.319230405352425d-1*t1
     5         00*t15*t16*t17*t29*t53-1.7730768017841417d-1*t127*t15*t16
     6         *t17*t29*t30)*wght+t31*(5.583154418094601d+0*t11*t3*t33*t
     7         51*t67*t95+1.1166308836189202d+1*t13*t3*t4*t52*t67*t74+5.
     8         583154418094601d+0*t3*t39*t52*t60*t61*t74-5.5831544180946
     9         01d+0*t11*t125*t3*t33*t39*t4-5.583154418094601d+0*t112*t1
     :         26*t13*t2*t3-5.583154418094601d+0*t10*t104*t106*t125*t3)*
     ;         wght+3*t101*t57*wght+3*t102*t56*wght+Cmat3(iq,D3_GAA_GAA_
     <         GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = rhob**1.3333333333333333d+0
            t3 = -t1+param(2)+1.0d+0
            t4 = param(4)
            t5 = gammabb**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/rhob**2.6666666666666666d+0
            t9 = gammabb*t7*t8+t4*t5*t6+1.0d+0
            t10 = 1.0d+0-1.0d+0/t9
            t11 = param(5)
            t12 = gammabb*t11*t8+1.0d+0
            t13 = 1.0d+0-1.0d+0/t12
            t14 = -9.305257363491002d-1*t10*t13*t2*t3-9.305257363491002d
     1         -1*t1*t2
            t15 = param(7)
            t16 = param(8)
            t17 = 1/rhob
            t18 = 9.305257363491002d-1*t10*t13*t2*t3+9.305257363491002d-
     1         1*t1*t2
            t19 = t18**5.0d-1
            t20 = t16**2
            t21 = rhob**2
            t22 = 1/t21
            t23 = 1/t20
            t24 = 1/t18
            t25 = exp(-1.413716694115413d+1*t21*t23*t24)
            t26 = t25-1.0d+0
            t27 = -3.536776513153214d-2*t18*t20*t22*t26+t25-1.5d+0
            t28 = 1/t16
            t29 = 1/t19
            t30 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         8*t29)+2.659615202676212d-1*t16*t17*t19*t27
            t31 = -3.5461536035682834d-1*t15*t16*t17*t19*t30-param(6)+1.
     1         0d+0
            t32 = 5.641895835477563d-1
            t33 = 1/t12**2
            t34 = 1/rhob**2.3333333333333334d+0
            t35 = rhob**3.333333333333333d-1
            t36 = 1/rhob**6.333333333333333d+0
            t37 = 1/rhob**3.6666666666666664d+0
            t38 = -2.6666666666666666d+0*gammabb*t37*t7-5.33333333333333
     1         3d+0*t36*t4*t5
            t39 = 1/t9**2
            t40 = 9.305257363491002d-1*t13*t2*t3*t38*t39+1.2407009817988
     1         002d+0*t10*t13*t3*t35+1.2407009817988002d+0*t1*t35-2.4814
     2         019635976003d+0*gammabb*t10*t11*t3*t33*t34
            t41 = 1/t19**3
            t42 = 3.759942411946509d+0*t28*t29-1.8799712059732546d+0*rho
     1         b*t28*t40*t41
            t43 = rhob**3
            t44 = 1/t43
            t45 = 1/t18**2
            t46 = 1.413716694115413d+1*t21*t23*t40*t45-2.827433388230826
     1         d+1*rhob*t23*t24
            t47 = -3.536776513153214d-2*t18*t20*t22*t25*t46+t25*t46+7.07
     1         3553026306428d-2*t18*t20*t26*t44-3.536776513153214d-2*t20
     2         *t22*t26*t40
            t48 = 2.659615202676212d-1*t16*t17*t19*t47+3.544907701811040
     1         3d+0*t25*t32*t42+1.329807601338106d-1*t16*t17*t27*t29*t40
     2         -2.659615202676212d-1*t16*t19*t22*t27
            t49 = -3.5461536035682834d-1*t15*t16*t17*t19*t48-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t40+3.5461536035682834d-1*
     2         t15*t16*t19*t22*t30
            t50 = -9.305257363491002d-1*t13*t2*t3*t38*t39-1.240700981798
     1         8002d+0*t10*t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481
     2         4019635976003d+0*gammabb*t10*t11*t3*t33*t34
            t51 = 1/t2
            t52 = t7*t8+2*gammabb*t4*t6
            t53 = 9.305257363491002d-1*t13*t2*t3*t39*t52+9.3052573634910
     1         02d-1*t10*t11*t3*t33*t51
            t54 = 1.413716694115413d+1*t21*t23*t25*t45*t53-3.53677651315
     1         3214d-2*t20*t22*t26*t53-5.0d-1*t24*t25*t53
            t55 = 2.659615202676212d-1*t16*t17*t19*t54-6.66432440723758d
     1         +0*rhob*t25*t28*t32*t41*t53+1.329807601338106d-1*t16*t17*
     2         t27*t29*t53
            t56 = -3.5461536035682834d-1*t15*t16*t17*t19*t55-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t53
            t57 = -9.305257363491002d-1*t13*t2*t3*t39*t52-9.305257363491
     1         002d-1*t10*t11*t3*t33*t51
            t58 = t40**2
            t59 = 1/t19**5
            t60 = t11**2
            t61 = 1/t12**3
            t62 = rhob**6
            t63 = 1/t62
            t64 = 1/rhob**3.3333333333333337d+0
            t65 = 1/rhob**6.666666666666666d-1
            t66 = t38**2
            t67 = 1/t9**3
            t68 = 1/rhob**7.333333333333333d+0
            t69 = 1/rhob**4.666666666666667d+0
            t70 = 9.777777777777779d+0*gammabb*t69*t7+3.377777777777778d
     1         +1*t4*t5*t68
            t71 = 9.305257363491002d-1*t13*t2*t3*t39*t70-1.8610514726982
     1         003d+0*t13*t2*t3*t66*t67+4.135669939329334d-1*t10*t13*t3*
     2         t65+4.135669939329334d-1*t1*t65+2.481401963597601d+0*gamm
     3         abb*t10*t11*t3*t33*t64-1.3234143805853868d+1*t10*t3*t5*t6
     4         0*t61*t63+2.4814019635976003d+0*t13*t3*t35*t38*t39-4.9628
     5         039271952007d+0*gammabb*t11*t3*t33*t34*t38*t39
            t72 = -1.8799712059732546d+0*rhob*t28*t41*t71+2.819956808959
     1         8814d+0*rhob*t28*t58*t59-3.759942411946509d+0*t28*t40*t41
            t73 = rhob**4
            t74 = 1/t73
            t75 = 1/t18**3
            t76 = -2.827433388230826d+1*t21*t23*t58*t75+1.41371669411541
     1         3d+1*t21*t23*t45*t71+5.654866776461652d+1*rhob*t23*t40*t4
     2         5-2.827433388230826d+1*t23*t24
            t77 = t46**2
            t78 = -3.536776513153214d-2*t18*t20*t22*t25*t77+t25*t77-3.53
     1         6776513153214d-2*t18*t20*t22*t25*t76+t25*t76-2.1220659078
     2         919285d-1*t18*t20*t26*t74-3.536776513153214d-2*t20*t22*t2
     3         6*t71+1.4147106052612857d-1*t18*t20*t25*t44*t46-7.0735530
     4         26306428d-2*t20*t22*t25*t40*t46+1.4147106052612857d-1*t20
     5         *t26*t40*t44
            t79 = 2.659615202676212d-1*t16*t17*t19*t78+3.544907701811040
     1         3d+0*t25*t32*t72+1.329807601338106d-1*t16*t17*t27*t29*t71
     2         -6.649038006690531d-2*t16*t17*t27*t41*t58+2.6596152026762
     3         12d-1*t16*t17*t29*t40*t47-5.319230405352424d-1*t16*t19*t2
     4         2*t47+3.5449077018110403d+0*t25*t32*t42*t46+5.31923040535
     5         2424d-1*t16*t19*t27*t44-2.659615202676212d-1*t16*t22*t27*
     6         t29*t40
            t80 = -3.5461536035682834d-1*t15*t16*t17*t19*t79-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t71+8.865384008920708d-2*t
     2         15*t16*t17*t30*t41*t58-3.5461536035682834d-1*t15*t16*t17*
     3         t29*t40*t48+7.092307207136567d-1*t15*t16*t19*t22*t48-7.09
     4         2307207136567d-1*t15*t16*t19*t30*t44+3.5461536035682834d-
     5         1*t15*t16*t22*t29*t30*t40
            t81 = -9.305257363491002d-1*t13*t2*t3*t39*t70+1.861051472698
     1         2003d+0*t13*t2*t3*t66*t67-4.135669939329334d-1*t10*t13*t3
     2         *t65-4.135669939329334d-1*t1*t65-2.481401963597601d+0*gam
     3         mabb*t10*t11*t3*t33*t64+1.3234143805853868d+1*t10*t3*t5*t
     4         60*t61*t63-2.4814019635976003d+0*t13*t3*t35*t38*t39+4.962
     5         8039271952007d+0*gammabb*t11*t3*t33*t34*t38*t39
            t82 = rhob**5
            t83 = 1/t82
            t84 = -2.6666666666666666d+0*t37*t7-1.0666666666666666d+1*ga
     1         mmabb*t36*t4
            t85 = 9.305257363491002d-1*t13*t2*t3*t39*t84+4.9628039271952
     1         007d+0*gammabb*t10*t3*t60*t61*t83-1.8610514726982003d+0*t
     2         13*t2*t3*t38*t52*t67+1.2407009817988002d+0*t13*t3*t35*t39
     3         *t52-2.4814019635976003d+0*gammabb*t11*t3*t33*t34*t39*t52
     4         +9.305257363491002d-1*t11*t3*t33*t38*t39*t51-1.2407009817
     5         988002d+0*t10*t11*t3*t33*t34
            t86 = -1.8799712059732546d+0*rhob*t28*t41*t85+2.819956808959
     1         8814d+0*rhob*t28*t40*t53*t59-1.8799712059732546d+0*t28*t4
     2         1*t53
            t87 = 1.413716694115413d+1*t21*t23*t45*t85-2.827433388230826
     1         d+1*t21*t23*t40*t53*t75+2.827433388230826d+1*rhob*t23*t45
     2         *t53
            t88 = -3.536776513153214d-2*t18*t20*t22*t25*t87+t25*t87-3.53
     1         6776513153214d-2*t20*t22*t26*t85+1.413716694115413d+1*t21
     2         *t23*t25*t45*t46*t53-5.0d-1*t24*t25*t46*t53-3.53677651315
     3         3214d-2*t20*t22*t25*t46*t53-5.0d-1*t25*t40*t45*t53+7.0735
     4         53026306428d-2*t20*t26*t44*t53+1.0d+0*t17*t24*t25*t53
            t89 = 2.659615202676212d-1*t16*t17*t19*t88+3.544907701811040
     1         3d+0*t25*t32*t86+1.329807601338106d-1*t16*t17*t27*t29*t85
     2         +1.329807601338106d-1*t16*t17*t29*t40*t54-2.6596152026762
     3         12d-1*t16*t19*t22*t54+1.329807601338106d-1*t16*t17*t29*t4
     4         7*t53+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t53-6.6
     5         49038006690531d-2*t16*t17*t27*t40*t41*t53-1.3298076013381
     6         06d-1*t16*t22*t27*t29*t53
            t90 = -3.5461536035682834d-1*t15*t16*t17*t19*t89-1.773076801
     1         7841417d-1*t15*t16*t17*t29*t30*t85-1.7730768017841417d-1*
     2         t15*t16*t17*t29*t40*t55+3.5461536035682834d-1*t15*t16*t19
     3         *t22*t55-1.7730768017841417d-1*t15*t16*t17*t29*t48*t53+8.
     4         865384008920708d-2*t15*t16*t17*t30*t40*t41*t53+1.77307680
     5         17841417d-1*t15*t16*t22*t29*t30*t53
            t91 = -9.305257363491002d-1*t13*t2*t3*t39*t84-4.962803927195
     1         2007d+0*gammabb*t10*t3*t60*t61*t83+1.8610514726982003d+0*
     2         t13*t2*t3*t38*t52*t67-1.2407009817988002d+0*t13*t3*t35*t3
     3         9*t52+2.4814019635976003d+0*gammabb*t11*t3*t33*t34*t39*t5
     4         2-9.305257363491002d-1*t11*t3*t33*t38*t39*t51+1.240700981
     5         7988002d+0*t10*t11*t3*t33*t34
            t92 = 1/t16**3
            t93 = t53**2
            t94 = 1/t19**7
            t95 = t52**2
            t96 = -1.8610514726982003d+0*t13*t2*t3*t67*t95-1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74+1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74+1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t97 = 1/t16**4
            t98 = 1/t18**4
            t99 = 1.9985948912206125d+2*t25*t73*t93*t97*t98+1.4137166941
     1         15413d+1*t21*t23*t25*t45*t96-3.536776513153214d-2*t20*t22
     2         *t26*t96-5.0d-1*t24*t25*t96-3.534291735288533d+1*t21*t23*
     3         t25*t75*t93
            t100 = 2.659615202676212d-1*t16*t17*t19*t99-6.66432440723758
     1         d+0*rhob*t25*t28*t32*t41*t96+1.329807601338106d-1*t16*t17
     2         *t27*t29*t96-9.421466669512572d+1*t25*t32*t43*t92*t93*t94
     3         +9.996486610856369d+0*rhob*t25*t28*t32*t59*t93-6.64903800
     4         6690531d-2*t16*t17*t27*t41*t93+2.659615202676212d-1*t16*t
     5         17*t29*t53*t54
            t101 = -1.7730768017841417d-1*t15*t16*t17*t29*t30*t96+8.8653
     1         84008920708d-2*t15*t16*t17*t30*t41*t93-3.5461536035682834
     2         d-1*t15*t16*t17*t29*t53*t55-3.5461536035682834d-1*t100*t1
     3         5*t16*t17*t19
            t102 = 1.8610514726982003d+0*t13*t2*t3*t67*t95+1.86105147269
     1         82003d+0*t10*t3*t60*t61*t74-1.8610514726982003d+0*t13*t3*
     2         t39*t4*t74-1.8610514726982003d+0*t11*t3*t33*t39*t51*t52
            t103 = t40**3
            t104 = t11**3
            t105 = gammabb**3
            t106 = 1/t12**4
            t107 = 1/rhob**9.666666666666666d+0
            t108 = 1/rhob**7
            t109 = 1/rhob**4.333333333333333d+0
            t110 = 1/rhob**1.6666666666666669d+0
            t111 = t38**3
            t112 = 1/t9**4
            t113 = -4.5629629629629626d+1*gammabb*t7/rhob**5.66666666666
     1         6667d+0-2.477037037037037d+2*t4*t5/rhob**8.33333333333333
     2         4d+0
            t114 = -5.583154418094601d+0*t13*t2*t3*t38*t67*t70+3.7221029
     1         453964005d+0*t13*t3*t35*t39*t70-7.444205890792801d+0*gamm
     2         abb*t11*t3*t33*t34*t39*t70-7.444205890792801d+0*t13*t3*t3
     3         5*t66*t67+1.4888411781585603d+1*gammabb*t11*t3*t33*t34*t6
     4         6*t67+1.2407009817988002d+0*t13*t3*t38*t39*t65+7.44420589
     5         0792801d+0*gammabb*t11*t3*t33*t38*t39*t64-3.9702431417561
     6         607d+1*t3*t38*t39*t5*t60*t61*t63+9.263900664097708d+1*t10
     7         *t108*t3*t5*t60*t61+9.305257363491002d-1*t113*t13*t2*t3*t
     8         39-9.374185195813161d+0*gammabb*t10*t109*t11*t3*t33+5.583
     9         154418094601d+0*t111*t112*t13*t2*t3-2.7571132928862224d-1
     :         *t10*t110*t13*t3-1.0587315044683095d+2*t10*t104*t105*t106
     ;         *t107*t3-2.7571132928862224d-1*t1*t110
            t115 = 8.482300164692478d+1*t103*t21*t23*t98-8.4823001646924
     1         78d+1*t21*t23*t40*t71*t75-1.6964600329384957d+2*rhob*t23*
     2         t58*t75+8.482300164692478d+1*rhob*t23*t45*t71+8.482300164
     3         692478d+1*t23*t40*t45+1.413716694115413d+1*t114*t21*t23*t
     4         45
            t116 = t46**3
            t117 = 1/rhob**8.666666666666666d+0
            t118 = 9.777777777777779d+0*t69*t7+6.755555555555556d+1*gamm
     1         abb*t4*t68
            t119 = -3.7221029453964005d+0*t13*t2*t3*t38*t67*t84+2.481401
     1         9635976003d+0*t13*t3*t35*t39*t84-4.9628039271952007d+0*ga
     2         mmabb*t11*t3*t33*t34*t39*t84+9.925607854390401d+0*gammabb
     3         *t3*t38*t39*t60*t61*t83-1.8610514726982003d+0*t13*t2*t3*t
     4         52*t67*t70+9.305257363491002d-1*t11*t3*t33*t39*t51*t70-1.
     5         8610514726982003d+0*t11*t3*t33*t51*t66*t67-4.962803927195
     6         2007d+0*t13*t3*t35*t38*t52*t67+9.925607854390401d+0*gamma
     7         bb*t11*t3*t33*t34*t38*t52*t67+5.583154418094601d+0*t112*t
     8         13*t2*t3*t52*t66+4.135669939329334d-1*t13*t3*t39*t52*t65+
     9         2.481401963597601d+0*gammabb*t11*t3*t33*t39*t52*t64+2.894
     :         9689575305343d+0*t10*t11*t3*t33*t64-1.3234143805853868d+1
     ;         *t3*t39*t5*t52*t60*t61*t63-3.143109153890294d+1*gammabb*t
     <         10*t3*t60*t61*t63+3.9702431417561607d+1*t10*t104*t106*t11
     =         7*t3*t5-2.4814019635976003d+0*t11*t3*t33*t34*t38*t39+9.30
     >         5257363491002d-1*t118*t13*t2*t3*t39
            t120 = 8.482300164692478d+1*t21*t23*t53*t58*t98-5.6548667764
     1         61652d+1*t21*t23*t40*t75*t85+5.654866776461652d+1*rhob*t2
     2         3*t45*t85-2.827433388230826d+1*t21*t23*t53*t71*t75-1.1309
     3         733552923305d+2*rhob*t23*t40*t53*t75+2.827433388230826d+1
     4         *t23*t45*t53+1.413716694115413d+1*t119*t21*t23*t45
            t121 = 1/rhob**7.666666666666667d+0
            t122 = -2.4814019635976003d+0*t13*t3*t35*t67*t95+4.962803927
     1         1952007d+0*gammabb*t11*t3*t33*t34*t67*t95+5.5831544180946
     2         01d+0*t112*t13*t2*t3*t38*t95-3.7221029453964005d+0*t13*t2
     3         *t3*t52*t67*t84+1.8610514726982003d+0*t11*t3*t33*t39*t51*
     4         t84+9.925607854390401d+0*gammabb*t3*t39*t52*t60*t61*t83+7
     5         .444205890792801d+0*t10*t3*t60*t61*t83-7.444205890792801d
     6         +0*t13*t3*t39*t4*t83-3.7221029453964005d+0*t13*t3*t38*t4*
     7         t67*t74-1.8610514726982003d+0*t3*t38*t39*t60*t61*t74-3.72
     8         21029453964005d+0*t11*t3*t33*t38*t51*t52*t67-2.4814019635
     9         97601d+0*t11*t3*t33*t34*t39*t52-4.9628039271952007d+0*gam
     :         mabb*t11*t121*t3*t33*t39*t4-1.4888411781585603d+1*gammabb
     ;         *t10*t104*t106*t121*t3
            t123 = 8.482300164692478d+1*t21*t23*t40*t93*t98-2.8274333882
     1         30826d+1*t21*t23*t40*t75*t96+2.827433388230826d+1*rhob*t2
     2         3*t45*t96-5.654866776461652d+1*rhob*t23*t75*t93-5.6548667
     3         76461652d+1*t21*t23*t53*t75*t85+1.413716694115413d+1*t122
     4         *t21*t23*t45
            t124 = t53**3
            t125 = 1/rhob**6.666666666666667d+0
            t126 = t52**3
            t127 = -5.583154418094601d+0*t11*t3*t33*t51*t67*t95-1.116630
     1         8836189202d+1*t13*t3*t4*t52*t67*t74-5.583154418094601d+0*
     2         t3*t39*t52*t60*t61*t74+5.583154418094601d+0*t11*t125*t3*t
     3         33*t39*t4+5.583154418094601d+0*t112*t126*t13*t2*t3+5.5831
     4         54418094601d+0*t10*t104*t106*t125*t3
            fnc(iq) = t14*t31*wght+fnc(iq)
            Amat(iq,D1_RB) = t31*t50*wght+t14*t49*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t31*t57*wght+t14*t56*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t31*t81*wght+t14*t80*wght+2*t49*t50*wgh
     1         t+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t31*t91*wght+t14*t90*wght+t49*t57*wght
     1         +t50*t56*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 2*t56*t57*wght+t102*t31*wght+t101*t14
     1         *wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t14*(-3.5461536035682834d-1*t15*t16*
     1         t17*t19*(3.5449077018110403d+0*t25*t32*(-7.04989202239970
     2         4d+0*rhob*t103*t28*t94+8.459870426879645d+0*rhob*t28*t40*
     3         t59*t71-5.639913617919763d+0*t28*t41*t71+8.45987042687964
     4         5d+0*t28*t58*t59-1.8799712059732546d+0*rhob*t114*t28*t41)
     5         +2.659615202676212d-1*t16*t17*t19*(8.488263631567714d-1*t
     6         18*t20*t26*t83+2.1220659078919285d-1*t18*t20*t25*t44*t77-
     7         1.0610329539459643d-1*t20*t22*t25*t40*t77-1.0610329539459
     8         643d-1*t18*t20*t22*t25*t46*t76+3*t25*t46*t76+2.1220659078
     9         919285d-1*t18*t20*t25*t44*t76-1.0610329539459643d-1*t20*t
     :         22*t25*t40*t76-6.366197723675786d-1*t18*t20*t25*t46*t74-6
     ;         .366197723675786d-1*t20*t26*t40*t74-1.0610329539459643d-1
     <         *t20*t22*t25*t46*t71+2.1220659078919285d-1*t20*t26*t44*t7
     =         1+4.244131815783857d-1*t20*t25*t40*t44*t46-3.536776513153
     >         214d-2*t114*t20*t22*t26-3.536776513153214d-2*t116*t18*t20
     ?         *t22*t25-3.536776513153214d-2*t115*t18*t20*t22*t25+t116*t
     @         25+t115*t25)+3.989422804014318d-1*t16*t17*t29*t40*t78-7.9
     1         78845608028636d-1*t16*t19*t22*t78+3.5449077018110403d+0*t
     2         25*t32*t42*t77+3.5449077018110403d+0*t25*t32*t42*t76-1.59
     3         57691216057274d+0*t16*t19*t27*t74+7.089815403622081d+0*t2
     4         5*t32*t46*t72+3.989422804014318d-1*t16*t17*t29*t47*t71-1.
     5         994711402007159d-1*t16*t17*t27*t40*t41*t71-3.989422804014
     6         318d-1*t16*t22*t27*t29*t71+9.973557010035795d-2*t103*t16*
     7         t17*t27*t59-1.994711402007159d-1*t16*t17*t41*t47*t58+1.99
     8         4711402007159d-1*t16*t22*t27*t41*t58+1.5957691216057274d+
     9         0*t16*t19*t44*t47-7.978845608028636d-1*t16*t22*t29*t40*t4
     :         7+7.978845608028636d-1*t16*t27*t29*t40*t44+1.329807601338
     ;         106d-1*t114*t16*t17*t27*t29)-5.319230405352425d-1*t15*t16
     <         *t17*t29*t40*t79+1.063846081070485d+0*t15*t16*t19*t22*t79
     =         +2.12769216214097d+0*t15*t16*t19*t30*t74-5.31923040535242
     >         5d-1*t15*t16*t17*t29*t48*t71+2.6596152026762127d-1*t15*t1
     ?         6*t17*t30*t40*t41*t71+5.319230405352425d-1*t15*t16*t22*t2
     @         9*t30*t71-1.3298076013381063d-1*t103*t15*t16*t17*t30*t59+
     1         2.6596152026762127d-1*t15*t16*t17*t41*t48*t58-2.659615202
     2         6762127d-1*t15*t16*t22*t30*t41*t58-2.12769216214097d+0*t1
     3         5*t16*t19*t44*t48+1.063846081070485d+0*t15*t16*t22*t29*t4
     4         0*t48-1.063846081070485d+0*t15*t16*t29*t30*t40*t44-1.7730
     5         768017841417d-1*t114*t15*t16*t17*t29*t30)*wght+3*t49*t81*
     6         wght+3*t50*t80*wght+t31*(5.583154418094601d+0*t13*t2*t3*t
     7         38*t67*t70-3.7221029453964005d+0*t13*t3*t35*t39*t70+7.444
     8         205890792801d+0*gammabb*t11*t3*t33*t34*t39*t70+7.44420589
     9         0792801d+0*t13*t3*t35*t66*t67-1.4888411781585603d+1*gamma
     :         bb*t11*t3*t33*t34*t66*t67-1.2407009817988002d+0*t13*t3*t3
     ;         8*t39*t65-7.444205890792801d+0*gammabb*t11*t3*t33*t38*t39
     <         *t64+3.9702431417561607d+1*t3*t38*t39*t5*t60*t61*t63-9.26
     =         3900664097708d+1*t10*t108*t3*t5*t60*t61-9.305257363491002
     >         d-1*t113*t13*t2*t3*t39+9.374185195813161d+0*gammabb*t10*t
     ?         109*t11*t3*t33-5.583154418094601d+0*t111*t112*t13*t2*t3+2
     @         .7571132928862224d-1*t10*t110*t13*t3+1.0587315044683095d+
     1         2*t10*t104*t105*t106*t107*t3+2.7571132928862224d-1*t1*t11
     2         0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t14*(-3.5461536035682834d-1*t15*t16
     1         *t17*t19*(3.5449077018110403d+0*t25*t32*(-7.0498920223997
     2         04d+0*rhob*t28*t53*t58*t94+5.639913617919763d+0*rhob*t28*
     3         t40*t59*t85-3.759942411946509d+0*t28*t41*t85+2.8199568089
     4         598814d+0*rhob*t28*t53*t59*t71+5.639913617919763d+0*t28*t
     5         40*t53*t59-1.8799712059732546d+0*rhob*t119*t28*t41)+2.659
     6         615202676212d-1*t16*t17*t29*t40*t88-5.319230405352424d-1*
     7         t16*t19*t22*t88+2.659615202676212d-1*t16*t17*t19*(-7.0735
     8         53026306428d-2*t18*t20*t22*t25*t46*t87+2*t25*t46*t87+1.41
     9         47106052612857d-1*t18*t20*t25*t44*t87-7.073553026306428d-
     :         2*t20*t22*t25*t40*t87-7.073553026306428d-2*t20*t22*t25*t4
     ;         6*t85+1.4147106052612857d-1*t20*t26*t44*t85+1.41371669411
     <         5413d+1*t21*t23*t25*t45*t53*t77-5.0d-1*t24*t25*t53*t77-3.
     =         536776513153214d-2*t20*t22*t25*t53*t77+1.413716694115413d
     >         +1*t21*t23*t25*t45*t53*t76-5.0d-1*t24*t25*t53*t76-3.53677
     ?         6513153214d-2*t20*t22*t25*t53*t76-2.1220659078919285d-1*t
     @         20*t26*t53*t74-5.0d-1*t25*t45*t53*t71-1.0d+0*t25*t40*t45*
     1         t46*t53+1.4147106052612857d-1*t20*t25*t44*t46*t53+2.0d+0*
     2         t17*t24*t25*t46*t53+2.0d+0*t17*t25*t40*t45*t53-3.0d+0*t22
     3         *t24*t25*t53-3.536776513153214d-2*t119*t20*t22*t26-3.5367
     4         76513153214d-2*t120*t18*t20*t22*t25+t120*t25)+3.544907701
     5         8110403d+0*t25*t32*t42*t87+3.5449077018110403d+0*t25*t32*
     6         t46*t86+2.659615202676212d-1*t16*t17*t29*t47*t85-1.329807
     7         601338106d-1*t16*t17*t27*t40*t41*t85-2.659615202676212d-1
     8         *t16*t22*t27*t29*t85+1.329807601338106d-1*t16*t17*t29*t53
     9         *t78+5.01149519714857d+1*t21*t23*t25*t32*t45*t53*t72+1.32
     :         9807601338106d-1*t16*t17*t29*t54*t71-6.649038006690531d-2
     ;         *t16*t17*t27*t41*t53*t71+9.973557010035795d-2*t16*t17*t27
     <         *t53*t58*t59-6.649038006690531d-2*t16*t17*t41*t54*t58+5.3
     =         19230405352424d-1*t16*t19*t44*t54-2.659615202676212d-1*t1
     >         6*t22*t29*t40*t54-1.329807601338106d-1*t16*t17*t40*t41*t4
     ?         7*t53-2.659615202676212d-1*t16*t22*t29*t47*t53+5.01149519
     @         714857d+1*t21*t23*t25*t32*t42*t45*t46*t53+2.6596152026762
     1         12d-1*t16*t27*t29*t44*t53+1.329807601338106d-1*t16*t22*t2
     2         7*t40*t41*t53+1.329807601338106d-1*t119*t16*t17*t27*t29)-
     3         3.5461536035682834d-1*t15*t16*t17*t29*t40*t89+7.092307207
     4         136567d-1*t15*t16*t19*t22*t89-3.5461536035682834d-1*t15*t
     5         16*t17*t29*t48*t85+1.7730768017841417d-1*t15*t16*t17*t30*
     6         t40*t41*t85+3.5461536035682834d-1*t15*t16*t22*t29*t30*t85
     7         -1.7730768017841417d-1*t15*t16*t17*t29*t53*t79-1.77307680
     8         17841417d-1*t15*t16*t17*t29*t55*t71+8.865384008920708d-2*
     9         t15*t16*t17*t30*t41*t53*t71-1.3298076013381063d-1*t15*t16
     :         *t17*t30*t53*t58*t59+8.865384008920708d-2*t15*t16*t17*t41
     ;         *t55*t58-7.092307207136567d-1*t15*t16*t19*t44*t55+3.54615
     <         36035682834d-1*t15*t16*t22*t29*t40*t55+1.7730768017841417
     =         d-1*t15*t16*t17*t40*t41*t48*t53+3.5461536035682834d-1*t15
     >         *t16*t22*t29*t48*t53-3.5461536035682834d-1*t15*t16*t29*t3
     ?         0*t44*t53-1.7730768017841417d-1*t15*t16*t22*t30*t40*t41*t
     @         53-1.7730768017841417d-1*t119*t15*t16*t17*t29*t30)*wght+2
     1         *t49*t91*wght+2*t50*t90*wght+t31*(3.7221029453964005d+0*t
     2         13*t2*t3*t38*t67*t84-2.4814019635976003d+0*t13*t3*t35*t39
     3         *t84+4.9628039271952007d+0*gammabb*t11*t3*t33*t34*t39*t84
     4         -9.925607854390401d+0*gammabb*t3*t38*t39*t60*t61*t83+1.86
     5         10514726982003d+0*t13*t2*t3*t52*t67*t70-9.305257363491002
     6         d-1*t11*t3*t33*t39*t51*t70+1.8610514726982003d+0*t11*t3*t
     7         33*t51*t66*t67+4.9628039271952007d+0*t13*t3*t35*t38*t52*t
     8         67-9.925607854390401d+0*gammabb*t11*t3*t33*t34*t38*t52*t6
     9         7-5.583154418094601d+0*t112*t13*t2*t3*t52*t66-4.135669939
     :         329334d-1*t13*t3*t39*t52*t65-2.481401963597601d+0*gammabb
     ;         *t11*t3*t33*t39*t52*t64-2.8949689575305343d+0*t10*t11*t3*
     <         t33*t64+1.3234143805853868d+1*t3*t39*t5*t52*t60*t61*t63+3
     =         .143109153890294d+1*gammabb*t10*t3*t60*t61*t63-3.97024314
     >         17561607d+1*t10*t104*t106*t117*t3*t5+2.4814019635976003d+
     ?         0*t11*t3*t33*t34*t38*t39-9.305257363491002d-1*t118*t13*t2
     @         *t3*t39)*wght+t56*t81*wght+t57*t80*wght+Cmat3(iq,D3_RB_RB
     1         _GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t14*(-3.5461536035682834d-1*t15*t1
     1         6*t17*t19*(1.329807601338106d-1*t16*t17*t29*t40*t99-2.659
     2         615202676212d-1*t16*t19*t22*t99+2.659615202676212d-1*t16*
     3         t17*t19*(1.9985948912206125d+2*t25*t46*t73*t93*t97*t98-7.
     4         068583470577066d+0*t21*t23*t25*t40*t93*t98+1.413716694115
     5         413d+1*t21*t23*t25*t45*t46*t96-5.0d-1*t24*t25*t46*t96-3.5
     6         36776513153214d-2*t20*t22*t25*t46*t96-5.0d-1*t25*t40*t45*
     7         t96+7.073553026306428d-2*t20*t26*t44*t96+1.0d+0*t17*t24*t
     8         25*t96-3.534291735288533d+1*t21*t23*t25*t46*t75*t93+1.0d+
     9         0*t25*t40*t75*t93+1.413716694115413d+1*rhob*t23*t25*t75*t
     :         93+2.827433388230826d+1*t21*t23*t25*t45*t53*t87-1.0d+0*t2
     ;         4*t25*t53*t87-7.073553026306428d-2*t20*t22*t25*t53*t87-1.
     <         0d+0*t25*t45*t53*t85-3.536776513153214d-2*t122*t20*t22*t2
     =         6-3.536776513153214d-2*t123*t18*t20*t22*t25+t123*t25)+7.0
     >         84834422688148d+2*t25*t32*t42*t73*t93*t97*t98+3.544907701
     ?         8110403d+0*t25*t32*(2.8199568089598814d+0*rhob*t28*t40*t5
     @         9*t96-1.8799712059732546d+0*t28*t41*t96-7.049892022399704
     1         d+0*rhob*t28*t40*t93*t94+2.8199568089598814d+0*t28*t59*t9
     2         3+5.639913617919763d+0*rhob*t28*t53*t59*t85-1.87997120597
     3         32546d+0*rhob*t122*t28*t41)+1.329807601338106d-1*t16*t17*
     4         t29*t47*t96+5.01149519714857d+1*t21*t23*t25*t32*t42*t45*t
     5         96-6.649038006690531d-2*t16*t17*t27*t40*t41*t96-1.3298076
     6         01338106d-1*t16*t22*t27*t29*t96-1.002299039429714d+2*t21*
     7         t23*t25*t32*t42*t75*t93+9.973557010035795d-2*t16*t17*t27*
     8         t40*t59*t93-6.649038006690531d-2*t16*t17*t41*t47*t93+6.64
     9         9038006690531d-2*t16*t22*t27*t41*t93+2.659615202676212d-1
     :         *t16*t17*t29*t53*t88+1.002299039429714d+2*t21*t23*t25*t32
     ;         *t45*t53*t86+2.659615202676212d-1*t16*t17*t29*t54*t85-1.3
     <         29807601338106d-1*t16*t17*t27*t41*t53*t85-1.3298076013381
     =         06d-1*t16*t17*t40*t41*t53*t54-2.659615202676212d-1*t16*t2
     >         2*t29*t53*t54+1.329807601338106d-1*t122*t16*t17*t27*t29)-
     ?         1.7730768017841417d-1*t15*t16*t17*t29*t48*t96+8.865384008
     @         920708d-2*t15*t16*t17*t30*t40*t41*t96+1.7730768017841417d
     1         -1*t15*t16*t22*t29*t30*t96-1.3298076013381063d-1*t15*t16*
     2         t17*t30*t40*t59*t93+8.865384008920708d-2*t15*t16*t17*t41*
     3         t48*t93-8.865384008920708d-2*t15*t16*t22*t30*t41*t93-3.54
     4         61536035682834d-1*t15*t16*t17*t29*t53*t89-3.5461536035682
     5         834d-1*t15*t16*t17*t29*t55*t85+1.7730768017841417d-1*t15*
     6         t16*t17*t30*t41*t53*t85+1.7730768017841417d-1*t15*t16*t17
     7         *t40*t41*t53*t55+3.5461536035682834d-1*t15*t16*t22*t29*t5
     8         3*t55-1.7730768017841417d-1*t100*t15*t16*t17*t29*t40-1.77
     9         30768017841417d-1*t122*t15*t16*t17*t29*t30+3.546153603568
     :         2834d-1*t100*t15*t16*t19*t22)*wght+t31*(2.481401963597600
     ;         3d+0*t13*t3*t35*t67*t95-4.9628039271952007d+0*gammabb*t11
     <         *t3*t33*t34*t67*t95-5.583154418094601d+0*t112*t13*t2*t3*t
     =         38*t95+3.7221029453964005d+0*t13*t2*t3*t52*t67*t84-1.8610
     >         514726982003d+0*t11*t3*t33*t39*t51*t84-9.925607854390401d
     ?         +0*gammabb*t3*t39*t52*t60*t61*t83-7.444205890792801d+0*t1
     @         0*t3*t60*t61*t83+7.444205890792801d+0*t13*t3*t39*t4*t83+3
     1         .7221029453964005d+0*t13*t3*t38*t4*t67*t74+1.861051472698
     2         2003d+0*t3*t38*t39*t60*t61*t74+3.7221029453964005d+0*t11*
     3         t3*t33*t38*t51*t52*t67+2.481401963597601d+0*t11*t3*t33*t3
     4         4*t39*t52+4.9628039271952007d+0*gammabb*t11*t121*t3*t33*t
     5         39*t4+1.4888411781585603d+1*gammabb*t10*t104*t106*t121*t3
     6         )*wght+2*t56*t91*wght+2*t57*t90*wght+t101*t50*wght+t102*t
     7         49*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t14*(-3.5461536035682834d-1*t15*t
     1         16*t17*t19*(3.989422804014318d-1*t16*t17*t29*t53*t99+2.65
     2         9615202676212d-1*t16*t17*t19*(5.995784673661838d+2*t25*t5
     3         3*t73*t96*t97*t98+1.0602875205865599d+2*t124*t21*t23*t25*
     4         t98-1.2990866792933983d+3*t124*t25*t73*t97/t18**5-1.06028
     5         75205865599d+2*t21*t23*t25*t53*t75*t96+2.8254469624923584
     6         d+3*t124*t25*t62/(t16**6*t18**6)+1.413716694115413d+1*t12
     7         7*t21*t23*t25*t45-3.536776513153214d-2*t127*t20*t22*t26-5
     8         .0d-1*t127*t24*t25)-2.826440000853771d+2*t25*t32*t43*t53*
     9         t92*t94*t96+2.998945983256911d+1*rhob*t25*t28*t32*t53*t59
     :         *t96+3.989422804014318d-1*t16*t17*t29*t54*t96-1.994711402
     ;         007159d-1*t16*t17*t27*t41*t53*t96-2.4991216527140922d+1*r
     <         hob*t124*t25*t28*t32*t94-1.994711402007159d-1*t16*t17*t41
     =         *t54*t93+4.710733334756285d+2*t124*t25*t32*t43*t92/t19**9
     >         -1.3319284713741864d+3*t124*t25*t32*t82/(t16**5*t19**11)+
     ?         9.973557010035795d-2*t124*t16*t17*t27*t59-6.6643244072375
     @         8d+0*rhob*t127*t25*t28*t32*t41+1.329807601338106d-1*t127*
     1         t16*t17*t27*t29)-5.319230405352425d-1*t15*t16*t17*t29*t55
     2         *t96+2.6596152026762127d-1*t15*t16*t17*t30*t41*t53*t96+2.
     3         6596152026762127d-1*t15*t16*t17*t41*t55*t93-1.32980760133
     4         81063d-1*t124*t15*t16*t17*t30*t59-5.319230405352425d-1*t1
     5         00*t15*t16*t17*t29*t53-1.7730768017841417d-1*t127*t15*t16
     6         *t17*t29*t30)*wght+t31*(5.583154418094601d+0*t11*t3*t33*t
     7         51*t67*t95+1.1166308836189202d+1*t13*t3*t4*t52*t67*t74+5.
     8         583154418094601d+0*t3*t39*t52*t60*t61*t74-5.5831544180946
     9         01d+0*t11*t125*t3*t33*t39*t4-5.583154418094601d+0*t112*t1
     :         26*t13*t2*t3-5.583154418094601d+0*t10*t104*t106*t125*t3)*
     ;         wght+3*t101*t57*wght+3*t102*t56*wght+Cmat3(iq,D3_GBB_GBB_
     <         GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
