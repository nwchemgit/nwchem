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
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& .9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}+.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{\it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& .9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}+.9305257363491002\,{\it t_1}\,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& .9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}+.9305257363491002\,{\it t_1}\,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}-.9305257363491002\,{\it t_1}\,{\it t_{21}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left({\it t_{29}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(-.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}-.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{\it t_{15}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({\it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}-.9305257363491002\,{\it t_1}\,{\it t_{30}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left({\it t_{38}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_cams12(param,tol_rho,ipol,nq,wght, 
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
            t45 = 2*t14*(-3.5461536035682834d-1*t15*t16*t17*t19*(2.65961
     1         5202676212d-1*t16*t17*t19*(1.413716694115413d+1*t21*t23*t
     2         25*t40*t44-3.536776513153214d-2*t20*t22*t26*t44-5.0d-1*t2
     3         4*t25*t44)-6.66432440723758d+0*rhoa*t25*t28*t32*t39*t44+1
     4         .329807601338106d-1*t16*t17*t27*t29*t44)-1.77307680178414
     5         17d-1*t15*t16*t17*t29*t30*t44)*wght+2*t31*(-9.30525736349
     6         1002d-1*t13*t2*t3*t37*t43-9.305257363491002d-1*t10*t11*t3
     7         *t33*t42)*wght
            fnc(iq) = 2*t14*t31*wght+fnc(iq)
            Amat(iq,D1_RA) = 2*t14*(-3.5461536035682834d-1*t15*t16*t17*t
     1         19*(2.659615202676212d-1*t16*t17*t19*(-3.536776513153214d
     2         -2*t18*t20*t22*t25*t41+t25*t41-3.536776513153214d-2*t20*t
     3         22*t26*t38+7.073553026306428d-2*t18*t20*t26/rhoa**3)+3.54
     4         49077018110403d+0*t25*t32*(3.759942411946509d+0*t28*t29-1
     5         .8799712059732546d+0*rhoa*t28*t38*t39)+1.329807601338106d
     6         -1*t16*t17*t27*t29*t38-2.659615202676212d-1*t16*t19*t22*t
     7         27)-1.7730768017841417d-1*t15*t16*t17*t29*t30*t38+3.54615
     8         36035682834d-1*t15*t16*t19*t22*t30)*wght+2*t31*(-9.305257
     9         363491002d-1*t13*t2*t3*t36*t37-1.2407009817988002d+0*t10*
     :         t13*t3*t35-1.2407009817988002d+0*t1*t35+2.481401963597600
     ;         3d+0*gammaaa*t10*t11*t3*t33*t34)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t45+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t45+Cmat(iq,D1_GAB)
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
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& .9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}+.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{\it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& .9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}+.9305257363491002\,{\it t_1}\,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& .9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}+.9305257363491002\,{\it t_1}\,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}-.9305257363491002\,{\it t_1}\,{\it t_{21}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left({\it t_{29}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(-.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}-.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{\it t_{15}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({\it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}-.9305257363491002\,{\it t_1}\,{\it t_{30}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left({\it t_{38}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_cams12_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_cams12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)-{\it t_1}+1.0\\\\ 
C>   {\it t_4} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_7} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_9} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_5}\,\sigma_{\alpha\alpha}+1.0}}\\\\ 
C>   {\it t_{10}} &=& -{\it param}\left(6\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{13}} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{14}} &=& .9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}+.9305257363491002\,{\it t_1}\,{\it t_2}\\\\ 
C>   {\it t_{15}} &=& \sqrt{{\it t_{14}}}\\\\ 
C>   {\it t_{16}} &=& {\it t_{12}}^2\\\\ 
C>   {\it t_{17}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_{16}}}}\\\\ 
C>   {\it t_{19}} &=& e^ {- {{14.13716694115413\,{\it t_{17}}\,{\it t_{18}}}\over{{\it t_{14}}}} }\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{21}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{23}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{24}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{22}}\,\sigma_{\beta\beta}+1.0}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{26}} &=& .9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}+.9305257363491002\,{\it t_1}\,{\it t_{21}}\\\\ 
C>   {\it t_{27}} &=& \sqrt{{\it t_{26}}}\\\\ 
C>   {\it t_{28}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{29}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{28}}}\over{{\it t_{26}}}} }\\\\ 
C>   {\it t_{30}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{31}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{32}} &=& 1.0-{{1.0}\over{{\it t_4}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{33}} &=& 1.0-{{1.0}\over{{{{\it t_8}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_7}\,{\it t_{31}}\,\sigma_{ss}+1.0}}\\\\ 
C>   {\it t_{34}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{35}} &=& .9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}+.9305257363491002\,{\it t_1}\,{\it t_{30}}\\\\ 
C>   {\it t_{36}} &=& \sqrt{{\it t_{35}}}\\\\ 
C>   {\it t_{37}} &=& \rho_s^2\\\\ 
C>   {\it t_{38}} &=& e^ {- {{14.13716694115413\,{\it t_{18}}\,{\it t_{37}}}\over{{\it t_{35}}}} }\\\\ 
C>   f &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{21}}\,{\it t_{23}}\,{\it t_{24}}-.9305257363491002\,{\it t_1}\,{\it t_{21}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\beta}\over{{\it t_{27}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{25}}\,{\it t_{27}}\,\left({\it t_{29}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{26}}\,\left({\it t_{29}}-1.0\right)}\over{{\it t_{28}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)+\left(-.9305257363491002\,{\it t_3}\,{\it t_2}\,{\it t_6}\,{\it t_9}-.9305257363491002\,{\it t_1}\,{\it t_2}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_\alpha}\over{{\it t_{15}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{13}}\,{\it t_{15}}\,\left({\it t_{19}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{14}}\,\left({\it t_{19}}-1.0\right)}\over{{\it t_{17}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& \left(-.9305257363491002\,{\it t_3}\,{\it t_{30}}\,{\it t_{32}}\,{\it t_{33}}-.9305257363491002\,{\it t_1}\,{\it t_{30}}\right)\,\left(-.3546153603568283\,{\it t_{11}}\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{20}}\,\rho_s}\over{{\it t_{36}}}}\right)+.2659615202676212\,{\it t_{12}}\,{\it t_{34}}\,{\it t_{36}}\,\left({\it t_{38}}-{{.03536776513153214\,{\it t_{16}}\,{\it t_{35}}\,\left({\it t_{38}}-1.0\right)}\over{{\it t_{37}}}}-1.5\right)\right)+{\it t_{10}}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_cams12_d3(param,tol_rho,ipol,nq,wght, 
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
