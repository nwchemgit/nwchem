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
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}+.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{\it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}\,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{\it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}+.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}\,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{\it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}+.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}-.9305257363491002\,{\it t_{16}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\beta}\over{{\it t_{22}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({\it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)+{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}-.9305257363491002\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{\it t_{10}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left({\it t_{14}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}-.9305257363491002\,{\it t_{25}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{\it t_{31}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({\it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camb88(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Cmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
            t38 = 2*t7*(-3.5461536035682834d-1*t10*t12*t8*t9*(2.65961520
     1         2676212d-1*t10*t12*(1.413716694115413d+1*t14*t16*t18*t34*
     2         t37-3.536776513153214d-2*t13*t15*t19*t37-5.0d-1*t17*t18*t
     3         37)*t9+1.329807601338106d-1*t10*t20*t22*t37*t9-6.66432440
     4         723758d+0*rhoa*t18*t21*t25*t33*t37)-1.7730768017841417d-1
     5         *t10*t22*t23*t37*t8*t9)*wght+2*t24*(4.2d-3*gammaaa*t3*t30
     6         *t36-4.2d-3*t3*t6)*wght
            fnc(iq) = 2*t24*t7*wght+fnc(iq)
            Amat(iq,D1_RA) = 2*t7*(-3.5461536035682834d-1*t10*t12*t8*t9*
     1         (2.659615202676212d-1*t10*t12*(-3.536776513153214d-2*t11*
     2         t13*t15*t18*t35+t18*t35-3.536776513153214d-2*t13*t15*t19*
     3         t32+7.073553026306428d-2*t11*t13*t19/rhoa**3)*t9+1.329807
     4         601338106d-1*t10*t20*t22*t32*t9-2.659615202676212d-1*t12*
     5         t15*t20*t9+3.5449077018110403d+0*t18*t25*(3.7599424119465
     6         09d+0*t21*t22-1.8799712059732546d+0*rhoa*t21*t32*t33))-1.
     7         7730768017841417d-1*t10*t22*t23*t32*t8*t9+3.5461536035682
     8         834d-1*t12*t15*t23*t8*t9)*wght+2*t24*(5.599999999999999d-
     9         3*gammaaa*t26*t6-1.2407009817988002d+0*t31+4.2d-3*gammaaa
     :         *t29*t3*t30)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t38+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t38+Cmat(iq,D1_GAB)
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
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}+.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{\it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}\,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{\it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}+.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}\,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{\it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}+.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}-.9305257363491002\,{\it t_{16}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\beta}\over{{\it t_{22}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({\it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)+{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}-.9305257363491002\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{\it t_{10}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left({\it t_{14}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}-.9305257363491002\,{\it t_{25}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{\it t_{31}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({\it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camb88_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Cmat,Cmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
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
C>   {\it t_4} &=& {{1}\over{0.0252\,{\it t_2}\,{\rm asinh}\; \left({\it t_2}\,{\it t_3}\right)\,{\it t_3}+1.0}}\\\\ 
C>   {\it t_5} &=& -{\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_9} &=& 0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}+.9305257363491002\,{\it t_1}\\\\ 
C>   {\it t_{10}} &=& \sqrt{{\it t_9}}\\\\ 
C>   {\it t_{11}} &=& {\it t_7}^2\\\\ 
C>   {\it t_{12}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{14}} &=& e^ {- {{14.13716694115413\,{\it t_{12}}\,{\it t_{13}}}\over{{\it t_9}}} }\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_7}}}\\\\ 
C>   {\it t_{16}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{18}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{0.0252\,{\it t_{17}}\,{\rm asinh}\; \left({\it t_{17}}\,{\it t_{18}}\right)\,{\it t_{18}}+1.0}}\\\\ 
C>   {\it t_{20}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{21}} &=& 0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}+.9305257363491002\,{\it t_{16}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{{\it t_{21}}}\\\\ 
C>   {\it t_{23}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{24}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{23}}}\over{{\it t_{21}}}} }\\\\ 
C>   {\it t_{25}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{25}}}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{28}} &=& {{1}\over{0.0252\,{\it t_{26}}\,{\rm asinh}\; \left({\it t_{26}}\,{\it t_{27}}\right)\,{\it t_{27}}+1.0}}\\\\ 
C>   {\it t_{29}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{30}} &=& 0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}+.9305257363491002\,{\it t_{25}}\\\\ 
C>   {\it t_{31}} &=& \sqrt{{\it t_{30}}}\\\\ 
C>   {\it t_{32}} &=& \rho_s^2\\\\ 
C>   {\it t_{33}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{32}}}\over{{\it t_{30}}}} }\\\\ 
C>   f &=& \left(-0.0042\,{\it t_{17}}\,{\it t_{19}}\,\sigma_{\beta\beta}-.9305257363491002\,{\it t_{16}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\beta}\over{{\it t_{22}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{20}}\,{\it t_{22}}\,\left({\it t_{24}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{21}}\,\left({\it t_{24}}-1.0\right)}\over{{\it t_{23}}}}-1.5\right)\right)+{\it t_5}+1.0\right)+\left(-0.0042\,{\it t_2}\,{\it t_4}\,\sigma_{\alpha\alpha}-.9305257363491002\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_\alpha}\over{{\it t_{10}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_8}\,{\it t_{10}}\,\left({\it t_{14}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_9}\,\left({\it t_{14}}-1.0\right)}\over{{\it t_{12}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-0.0042\,{\it t_{26}}\,{\it t_{28}}\,\sigma_{ss}-.9305257363491002\,{\it t_{25}}\right)\,\left(-.3546153603568283\,{\it t_6}\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{15}}\,\rho_s}\over{{\it t_{31}}}}\right)+.2659615202676212\,{\it t_7}\,{\it t_{29}}\,{\it t_{31}}\,\left({\it t_{33}}-{{.03536776513153214\,{\it t_{11}}\,{\it t_{30}}\,\left({\it t_{33}}-1.0\right)}\over{{\it t_{32}}}}-1.5\right)\right)+{\it t_5}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Yanai, DP Tew, NC Handy, Chem.Phys.Lett. 393, 51 (2004)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2004.06.011 "> 
C> 10.1016/j.cplett.2004.06.011 </a> 
C> 
      subroutine nwxcm_x_camb88_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
