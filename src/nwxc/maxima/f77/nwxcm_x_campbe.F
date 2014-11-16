C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_campbe.F 
C> The nwxcm_x_campbe functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_campbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& -{\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{10}} &=& 0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}+0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_{11}} &=& \sqrt{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {\it t_8}^2\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{15}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{14}}}\over{{\it t_{10}}}} }\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   {\it t_{17}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{20}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}+0.9305257363491\,{\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& \sqrt{{\it t_{20}}}\\\\ 
C>   {\it t_{22}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{23}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{22}}}\over{{\it t_{20}}}} }\\\\ 
C>   {\it t_{24}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{25}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{27}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}+0.9305257363491\,{\it t_{24}}\\\\ 
C>   {\it t_{28}} &=& \sqrt{{\it t_{27}}}\\\\ 
C>   {\it t_{29}} &=& \rho_s^2\\\\ 
C>   {\it t_{30}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{29}}}\over{{\it t_{27}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}-0.9305257363491\,{\it t_{17}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\beta}\over{{\it t_{21}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left({\it t_{23}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{20}}\,\left({\it t_{23}}-1.0\right)}\over{{\it t_{22}}}}-1.5\right)\right)+{\it t_6}+1.0\right)+1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}-0.9305257363491\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\alpha}\over{{\it t_{11}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left({\it t_{15}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{10}}\,\left({\it t_{15}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}-0.9305257363491\,{\it t_{24}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_s}\over{{\it t_{28}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left({\it t_{30}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{27}}\,\left({\it t_{30}}-1.0\right)}\over{{\it t_{29}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] OA Vydrov, GE Scuseria, JP Perdew, J.Chem.Phys. 126, 154109 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2723119 "> 
C> 10.1063/1.2723119 </a> 
C> 
      subroutine nwxcm_x_campbe(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = param(2)
            t4 = 1.6455307846020562d-2*gammaaa*t3/(rhoa**2.6666666666666
     1         666d+0*t2)+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = -9.305257363491d-1*t1*t2*t5-9.305257363491d-1*t1
            t7 = param(4)
            t8 = param(5)
            t9 = 1/rhoa
            t10 = 9.305257363491d-1*t1*t2*t5+9.305257363491d-1*t1
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
            t22 = 2.659615202676212d-1*t11*t19*t8*t9+1.7724538509055202d
     1         +0*erf(3.759942411946509d+0*rhoa*t20*t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9-param(3)+1.0d+
     1         0
            t24 = 5.641895835477563d-1
            t25 = 1/t4**2
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = rhoa**3.333333333333333d-1
            t28 = 1.2407009817987999d+0*t2*t27*t5-4.0832233200718415d-2*
     1         gammaaa*t25*t26*t3+1.2407009817987999d+0*t27
            t29 = 1/t11**3
            t30 = 1/t10**2
            t31 = 1.413716694115413d+1*t13*t15*t28*t30-2.827433388230826
     1         d+1*rhoa*t15*t16
            t32 = 1/t1
            t33 = 2.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t8*t9*(2.6596
     1         15202676212d-1*t11*(-7.656043725134702d-3*t16*t17*t25*t3*
     2         t32+2.164695365020097d-1*rhoa**6.666666666666666d-1*t15*t
     3         17*t25*t3*t30-5.415543126146092d-4*t12*t18*t25*t3/rhoa**3
     4         .3333333333333337d+0)*t8*t9+2.0362130283722077d-3*t19*t21
     5         *t25*t26*t3*t8-1.0204471812058663d-1*t17*t20*t24*t25*t29*
     6         t3/t27)-2.7149507044962773d-3*t21*t22*t25*t26*t3*t7*t8)*w
     7         ght-3.062417490053881d-2*t23*t25*t3*t32*wght
            fnc(iq) = 2.0d+0*t23*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = 2.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t8
     1         *t9*(2.659615202676212d-1*t11*(-3.536776513153214d-2*t10*
     2         t12*t14*t17*t31+t17*t31-3.536776513153214d-2*t12*t14*t18*
     3         t28+7.073553026306428d-2*t10*t12*t18/rhoa**3)*t8*t9+1.329
     4         807601338106d-1*t19*t21*t28*t8*t9-2.659615202676212d-1*t1
     5         1*t14*t19*t8+3.5449077018110403d+0*t17*t24*(3.75994241194
     6         6509d+0*t20*t21-1.8799712059732546d+0*rhoa*t20*t28*t29))-
     7         1.7730768017841417d-1*t21*t22*t28*t7*t8*t9+3.546153603568
     8         2834d-1*t11*t14*t22*t7*t8)*wght+2.0d+0*t23*(-1.2407009817
     9         987999d+0*t2*t27*t5+4.0832233200718415d-2*gammaaa*t25*t26
     :         *t3-1.2407009817987999d+0*t27)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t33+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t33+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = 1/t2
            t4 = param(2)
            t5 = 1.6455307846020562d-2*gammaaa*t3*t4/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = -9.305257363491d-1*t1*t2*t6-9.305257363491d-1*t1
            t8 = -param(3)
            t9 = param(4)
            t10 = param(5)
            t11 = 1/rhoa
            t12 = 9.305257363491d-1*t1*t2*t6+9.305257363491d-1*t1
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
            t26 = rhob**1.3333333333333333d+0
            t27 = 1.6455307846020562d-2*gammabb*t3*t4/rhob**2.6666666666
     1         666666d+0+1.0d+0
            t28 = 1.0d+0-1.0d+0/t27
            t29 = -9.305257363491d-1*t2*t26*t28-9.305257363491d-1*t26
            t30 = 1/rhob
            t31 = 9.305257363491d-1*t2*t26*t28+9.305257363491d-1*t26
            t32 = t31**5.0d-1
            t33 = rhob**2
            t34 = 1/t33
            t35 = 1/t31
            t36 = exp(-1.413716694115413d+1*t17*t33*t35)
            t37 = t36-1.0d+0
            t38 = -3.536776513153214d-2*t14*t31*t34*t37+t36-1.5d+0
            t39 = 1/t32
            t40 = 1.7724538509055202d+0*erf(3.759942411946509d+0*rhob*t2
     1         2*t39)+2.659615202676212d-1*t10*t30*t32*t38
            t41 = -3.5461536035682834d-1*t10*t30*t32*t40*t9+t8+1.0d+0
            t42 = 5.641895835477563d-1
            t43 = 1/t5**2
            t44 = 1/rhoa**2.3333333333333334d+0
            t45 = rhoa**3.333333333333333d-1
            t46 = 1.2407009817987999d+0*t2*t45*t6+1.2407009817987999d+0*
     1         t45-4.0832233200718415d-2*gammaaa*t4*t43*t44
            t47 = 1/t13**3
            t48 = 1/t12**2
            t49 = 1.413716694115413d+1*t15*t17*t46*t48-2.827433388230826
     1         d+1*rhoa*t17*t18
            t50 = 1/t27**2
            t51 = 1/rhob**2.3333333333333334d+0
            t52 = rhob**3.333333333333333d-1
            t53 = 1.2407009817987999d+0*t2*t28*t52+1.2407009817987999d+0
     1         *t52-4.0832233200718415d-2*gammabb*t4*t50*t51
            t54 = 1/t32**3
            t55 = 1/t31**2
            t56 = 1.413716694115413d+1*t17*t33*t53*t55-2.827433388230826
     1         d+1*rhob*t17*t35
            t57 = 1/t1
            t58 = 1/t26
            fnc(iq) = (1.0d+0*t25*t7+1.0d+0*t29*t41)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t7*(-3.5461536035682834d-1*t10*t11*
     1         t13*(2.659615202676212d-1*t10*t11*t13*(-3.536776513153214
     2         d-2*t12*t14*t16*t19*t49+t19*t49-3.536776513153214d-2*t14*
     3         t16*t20*t46+7.073553026306428d-2*t12*t14*t20/rhoa**3)+3.5
     4         449077018110403d+0*t19*t42*(3.759942411946509d+0*t22*t23-
     5         1.8799712059732546d+0*rhoa*t22*t46*t47)+1.329807601338106
     6         d-1*t10*t11*t21*t23*t46-2.659615202676212d-1*t10*t13*t16*
     7         t21)*t9-1.7730768017841417d-1*t10*t11*t23*t24*t46*t9+3.54
     8         61536035682834d-1*t10*t13*t16*t24*t9)+1.0d+0*t25*(-1.2407
     9         009817987999d+0*t2*t45*t6-1.2407009817987999d+0*t45+4.083
     :         2233200718415d-2*gammaaa*t4*t43*t44))*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t29*(-3.5461536035682834d-1*t10*t30
     1         *t32*(2.659615202676212d-1*t10*t30*t32*(-3.53677651315321
     2         4d-2*t14*t31*t34*t36*t56+t36*t56-3.536776513153214d-2*t14
     3         *t34*t37*t53+7.073553026306428d-2*t14*t31*t37/rhob**3)+3.
     4         5449077018110403d+0*t36*t42*(3.759942411946509d+0*t22*t39
     5         -1.8799712059732546d+0*rhob*t22*t53*t54)+1.32980760133810
     6         6d-1*t10*t30*t38*t39*t53-2.659615202676212d-1*t10*t32*t34
     7         *t38)*t9-1.7730768017841417d-1*t10*t30*t39*t40*t53*t9+3.5
     8         461536035682834d-1*t10*t32*t34*t40*t9)+1.0d+0*t41*(-1.240
     9         7009817987999d+0*t2*t28*t52-1.2407009817987999d+0*t52+4.0
     :         832233200718415d-2*gammabb*t4*t50*t51))*wght+Amat(iq,D1_R
     ;         B)
            Cmat(iq,D1_GAA) = (1.0d+0*t7*(-3.5461536035682834d-1*t10*t11
     1         *t13*(2.659615202676212d-1*t10*t11*t13*(-7.65604372513470
     2         2d-3*t18*t19*t4*t43*t57+2.164695365020097d-1*rhoa**6.6666
     3         66666666666d-1*t17*t19*t4*t43*t48-5.415543126146092d-4*t1
     4         4*t20*t4*t43/rhoa**3.3333333333333337d+0)-1.0204471812058
     5         663d-1*t19*t22*t4*t42*t43*t47/t45+2.0362130283722077d-3*t
     6         10*t21*t23*t4*t43*t44)*t9-2.7149507044962773d-3*t10*t23*t
     7         24*t4*t43*t44*t9)-1.5312087450269404d-2*t25*t4*t43*t57)*w
     8         ght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t29*(-3.5461536035682834d-1*t10*t3
     1         0*t32*(2.659615202676212d-1*t10*t30*t32*(-7.6560437251347
     2         02d-3*t35*t36*t4*t50*t58+2.164695365020097d-1*rhob**6.666
     3         666666666666d-1*t17*t36*t4*t50*t55-5.415543126146092d-4*t
     4         14*t37*t4*t50/rhob**3.3333333333333337d+0)-1.020447181205
     5         8663d-1*t22*t36*t4*t42*t50*t54/t52+2.0362130283722077d-3*
     6         t10*t38*t39*t4*t50*t51)*t9-2.7149507044962773d-3*t10*t39*
     7         t4*t40*t50*t51*t9)-1.5312087450269404d-2*t4*t41*t50*t58)*
     8         wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(1)
            t3 = param(2)
            t4 = 1.6455307846020562d-2*gammaaa*t3/(rhoa**2.6666666666666
     1         666d+0*t2)+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = -9.305257363491d-1*t1*t2*t5-9.305257363491d-1*t1
            t7 = param(4)
            t8 = param(5)
            t9 = 1/rhoa
            t10 = 9.305257363491d-1*t1*t2*t5+9.305257363491d-1*t1
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
            t22 = 2.659615202676212d-1*t11*t19*t8*t9+1.7724538509055202d
     1         +0*erf(3.759942411946509d+0*rhoa*t20*t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9-param(3)+1.0d+
     1         0
            t24 = 5.641895835477563d-1
            t25 = 1/t4**2
            t26 = 1/rhoa**2.3333333333333334d+0
            t27 = rhoa**3.333333333333333d-1
            t28 = 1.2407009817987999d+0*t2*t27*t5-4.0832233200718415d-2*
     1         gammaaa*t25*t26*t3+1.2407009817987999d+0*t27
            t29 = 1/t11**3
            t30 = 1/t10**2
            t31 = 1.413716694115413d+1*t13*t15*t28*t30-2.827433388230826
     1         d+1*rhoa*t15*t16
            t32 = 1/t1
            fnc(iq) = 1.0d+0*t23*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t8
     1         *t9*(2.659615202676212d-1*t11*(-3.536776513153214d-2*t10*
     2         t12*t14*t17*t31+t17*t31-3.536776513153214d-2*t12*t14*t18*
     3         t28+7.073553026306428d-2*t10*t12*t18/rhoa**3)*t8*t9+1.329
     4         807601338106d-1*t19*t21*t28*t8*t9-2.659615202676212d-1*t1
     5         1*t14*t19*t8+3.5449077018110403d+0*t17*t24*(3.75994241194
     6         6509d+0*t20*t21-1.8799712059732546d+0*rhoa*t20*t28*t29))-
     7         1.7730768017841417d-1*t21*t22*t28*t7*t8*t9+3.546153603568
     8         2834d-1*t11*t14*t22*t7*t8)*wght+1.0d+0*t23*(-1.2407009817
     9         987999d+0*t2*t27*t5+4.0832233200718415d-2*gammaaa*t25*t26
     :         *t3-1.2407009817987999d+0*t27)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t
     1         8*t9*(2.659615202676212d-1*t11*(-7.656043725134702d-3*t16
     2         *t17*t25*t3*t32+2.164695365020097d-1*rhoa**6.666666666666
     3         666d-1*t15*t17*t25*t3*t30-5.415543126146092d-4*t12*t18*t2
     4         5*t3/rhoa**3.3333333333333337d+0)*t8*t9+2.036213028372207
     5         7d-3*t19*t21*t25*t26*t3*t8-1.0204471812058663d-1*t17*t20*
     6         t24*t25*t29*t3/t27)-2.7149507044962773d-3*t21*t22*t25*t26
     7         *t3*t7*t8)*wght-1.5312087450269404d-2*t23*t25*t3*t32*wght
     8         +Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(1)
            t3 = param(2)
            t4 = 1.6455307846020562d-2*gammabb*t3/(rhob**2.6666666666666
     1         666d+0*t2)+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = -9.305257363491d-1*t1*t2*t5-9.305257363491d-1*t1
            t7 = param(4)
            t8 = param(5)
            t9 = 1/rhob
            t10 = 9.305257363491d-1*t1*t2*t5+9.305257363491d-1*t1
            t11 = t10**5.0d-1
            t12 = t8**2
            t13 = rhob**2
            t14 = 1/t13
            t15 = 1/t12
            t16 = 1/t10
            t17 = exp(-1.413716694115413d+1*t13*t15*t16)
            t18 = t17-1.0d+0
            t19 = -3.536776513153214d-2*t10*t12*t14*t18+t17-1.5d+0
            t20 = 1/t8
            t21 = 1/t11
            t22 = 2.659615202676212d-1*t11*t19*t8*t9+1.7724538509055202d
     1         +0*erf(3.759942411946509d+0*rhob*t20*t21)
            t23 = -3.5461536035682834d-1*t11*t22*t7*t8*t9-param(3)+1.0d+
     1         0
            t24 = 5.641895835477563d-1
            t25 = 1/t4**2
            t26 = 1/rhob**2.3333333333333334d+0
            t27 = rhob**3.333333333333333d-1
            t28 = 1.2407009817987999d+0*t2*t27*t5-4.0832233200718415d-2*
     1         gammabb*t25*t26*t3+1.2407009817987999d+0*t27
            t29 = 1/t11**3
            t30 = 1/t10**2
            t31 = 1.413716694115413d+1*t13*t15*t28*t30-2.827433388230826
     1         d+1*rhob*t15*t16
            t32 = 1/t1
            fnc(iq) = 1.0d+0*t23*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t8
     1         *t9*(2.659615202676212d-1*t11*(-3.536776513153214d-2*t10*
     2         t12*t14*t17*t31+t17*t31-3.536776513153214d-2*t12*t14*t18*
     3         t28+7.073553026306428d-2*t10*t12*t18/rhob**3)*t8*t9+1.329
     4         807601338106d-1*t19*t21*t28*t8*t9-2.659615202676212d-1*t1
     5         1*t14*t19*t8+3.5449077018110403d+0*t17*t24*(3.75994241194
     6         6509d+0*t20*t21-1.8799712059732546d+0*rhob*t20*t28*t29))-
     7         1.7730768017841417d-1*t21*t22*t28*t7*t8*t9+3.546153603568
     8         2834d-1*t11*t14*t22*t7*t8)*wght+1.0d+0*t23*(-1.2407009817
     9         987999d+0*t2*t27*t5+4.0832233200718415d-2*gammabb*t25*t26
     :         *t3-1.2407009817987999d+0*t27)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t6*(-3.5461536035682834d-1*t11*t7*t
     1         8*t9*(2.659615202676212d-1*t11*(-7.656043725134702d-3*t16
     2         *t17*t25*t3*t32+2.164695365020097d-1*rhob**6.666666666666
     3         666d-1*t15*t17*t25*t3*t30-5.415543126146092d-4*t12*t18*t2
     4         5*t3/rhob**3.3333333333333337d+0)*t8*t9+2.036213028372207
     5         7d-3*t19*t21*t25*t26*t3*t8-1.0204471812058663d-1*t17*t20*
     6         t24*t25*t29*t3/t27)-2.7149507044962773d-3*t21*t22*t25*t26
     7         *t3*t7*t8)*wght-1.5312087450269404d-2*t23*t25*t3*t32*wght
     8         +Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_campbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& -{\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{10}} &=& 0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}+0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_{11}} &=& \sqrt{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {\it t_8}^2\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{15}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{14}}}\over{{\it t_{10}}}} }\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   {\it t_{17}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{20}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}+0.9305257363491\,{\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& \sqrt{{\it t_{20}}}\\\\ 
C>   {\it t_{22}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{23}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{22}}}\over{{\it t_{20}}}} }\\\\ 
C>   {\it t_{24}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{25}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{27}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}+0.9305257363491\,{\it t_{24}}\\\\ 
C>   {\it t_{28}} &=& \sqrt{{\it t_{27}}}\\\\ 
C>   {\it t_{29}} &=& \rho_s^2\\\\ 
C>   {\it t_{30}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{29}}}\over{{\it t_{27}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}-0.9305257363491\,{\it t_{17}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\beta}\over{{\it t_{21}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left({\it t_{23}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{20}}\,\left({\it t_{23}}-1.0\right)}\over{{\it t_{22}}}}-1.5\right)\right)+{\it t_6}+1.0\right)+1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}-0.9305257363491\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\alpha}\over{{\it t_{11}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left({\it t_{15}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{10}}\,\left({\it t_{15}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}-0.9305257363491\,{\it t_{24}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_s}\over{{\it t_{28}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left({\it t_{30}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{27}}\,\left({\it t_{30}}-1.0\right)}\over{{\it t_{29}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] OA Vydrov, GE Scuseria, JP Perdew, J.Chem.Phys. 126, 154109 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2723119 "> 
C> 10.1063/1.2723119 </a> 
C> 
      subroutine nwxcm_x_campbe_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_campbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& -{\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_\alpha}}\\\\ 
C>   {\it t_{10}} &=& 0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}+0.9305257363491\,{\it t_1}\\\\ 
C>   {\it t_{11}} &=& \sqrt{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {\it t_8}^2\\\\ 
C>   {\it t_{13}} &=& \rho_\alpha^2\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_{12}}}}\\\\ 
C>   {\it t_{15}} &=& e^ {- {{14.13716694115413\,{\it t_{13}}\,{\it t_{14}}}\over{{\it t_{10}}}} }\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_8}}}\\\\ 
C>   {\it t_{17}} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{\rho_\beta}}\\\\ 
C>   {\it t_{20}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}+0.9305257363491\,{\it t_{17}}\\\\ 
C>   {\it t_{21}} &=& \sqrt{{\it t_{20}}}\\\\ 
C>   {\it t_{22}} &=& \rho_\beta^2\\\\ 
C>   {\it t_{23}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{22}}}\over{{\it t_{20}}}} }\\\\ 
C>   {\it t_{24}} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_{25}} &=& 1.0-{{1.0}\over{{{.01645530784602056\,{\it t_3}\,{\it t_4}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{\rho_s}}\\\\ 
C>   {\it t_{27}} &=& 0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}+0.9305257363491\,{\it t_{24}}\\\\ 
C>   {\it t_{28}} &=& \sqrt{{\it t_{27}}}\\\\ 
C>   {\it t_{29}} &=& \rho_s^2\\\\ 
C>   {\it t_{30}} &=& e^ {- {{14.13716694115413\,{\it t_{14}}\,{\it t_{29}}}\over{{\it t_{27}}}} }\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{17}}\,{\it t_{18}}-0.9305257363491\,{\it t_{17}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\beta}\over{{\it t_{21}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{19}}\,{\it t_{21}}\,\left({\it t_{23}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{20}}\,\left({\it t_{23}}-1.0\right)}\over{{\it t_{22}}}}-1.5\right)\right)+{\it t_6}+1.0\right)+1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_1}\,{\it t_5}-0.9305257363491\,{\it t_1}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_\alpha}\over{{\it t_{11}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_9}\,{\it t_{11}}\,\left({\it t_{15}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{10}}\,\left({\it t_{15}}-1.0\right)}\over{{\it t_{13}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\left(-0.9305257363491\,{\it t_2}\,{\it t_{24}}\,{\it t_{25}}-0.9305257363491\,{\it t_{24}}\right)\,\left(-.3546153603568283\,{\it t_7}\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left(1.77245385090552\,\mathrm{erf}\left({{3.759942411946509\,{\it t_{16}}\,\rho_s}\over{{\it t_{28}}}}\right)+.2659615202676212\,{\it t_8}\,{\it t_{26}}\,{\it t_{28}}\,\left({\it t_{30}}-{{.03536776513153214\,{\it t_{12}}\,{\it t_{27}}\,\left({\it t_{30}}-1.0\right)}\over{{\it t_{29}}}}-1.5\right)\right)+{\it t_6}+1.0\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] OA Vydrov, GE Scuseria, JP Perdew, J.Chem.Phys. 126, 154109 (2007)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2723119 "> 
C> 10.1063/1.2723119 </a> 
C> 
      subroutine nwxcm_x_campbe_d3(param,tol_rho,ipol,nq,wght, 
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
