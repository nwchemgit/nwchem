C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_sogga.F 
C> The nwxcm_x_sogga functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_sogga functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_4} &=& {{0.003680288926083876\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{\it t_4}+1.0}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{11}} &=& 1.0-{{1}\over{e^{{\it t_4}}}}\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{16}} &=& {{0.003680288926083876\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}
C>    \over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{
C>    \it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{
C>    \it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}
C>    +{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}
C>    +{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}
C>    +{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}
C>    +{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{
C>    \it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{
C>    \it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}
C>    +{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}
C>    +{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}
C>    +{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_sogga(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = 1/t5**2
            t26 = 1/rhoa**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t24*wght
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*(-3.9256415211561
     1         34d-2*gammaaa*t10*t25*t26*t9-2.9442311408671007d-2*gammaa
     2         a*t25*t26*t7*t8-1.962820760578067d-2*gammaaa*t2*t25*t26*t
     3         6-9.814103802890335d-3*gammaaa*t13*t25*t26-4.907051901445
     4         1675d-2*gammaaa*t11*t12*t25*t26-9.814103802890335d-3*gamm
     5         aaa*t15*t23*t26-4.9070519014451675d-2*gammaaa*t15*t21*t22
     6         *t26-3.925641521156134d-2*gammaaa*t15*t19*t20*t26-2.94423
     7         11408671007d-2*gammaaa*t15*t17*t18*t26-1.962820760578067d
     8         -2*gammaaa*t14*t15*t16*t26)-1.2407009817988002d+0*rhoa**3
     9         .333333333333333d-1*t24)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         25*t3*(1.4721155704335503d-2*t10*t9+1.1040866778251626d-2
     2         *t7*t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t
     3         13+1.8401444630419378d-2*t11*t12)+(3.6802889260838756d-3*
     4         t15*t23+1.8401444630419378d-2*t15*t21*t22+1.4721155704335
     5         503d-2*t15*t19*t20+1.1040866778251626d-2*t15*t17*t18+7.36
     6         0577852167751d-3*t14*t15*t16)*t3)*wght
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
            t2 = param(1)
            t3 = param(7)
            t4 = param(3)
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 3.6802889260838756d-3*gammaaa*t5
            t7 = t6+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = t8**2.0d+0
            t10 = param(4)
            t11 = t8**3.0d+0
            t12 = param(5)
            t13 = t8**4.0d+0
            t14 = param(6)
            t15 = param(2)
            t16 = param(9)
            t17 = exp(-t6)
            t18 = 1.0d+0-t17
            t19 = t18**2.0d+0
            t20 = param(10)
            t21 = t18**3.0d+0
            t22 = param(11)
            t23 = t18**4.0d+0
            t24 = param(12)
            t25 = param(8)
            t26 = t4*t9+t14*t8**5.0d+0+t15*t8+t3+t18*t25+t18**5.0d+0*t24
     1         +t22*t23+t20*t21+t2+t16*t19+t12*t13+t10*t11
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/rhob**2.6666666666666666d+0
            t29 = 3.6802889260838756d-3*gammabb*t28
            t30 = t29+1.0d+0
            t31 = 1.0d+0-1.0d+0/t30
            t32 = t31**2.0d+0
            t33 = t31**3.0d+0
            t34 = t31**4.0d+0
            t35 = exp(-t29)
            t36 = 1.0d+0-t35
            t37 = t36**2.0d+0
            t38 = t36**3.0d+0
            t39 = t36**4.0d+0
            t40 = t32*t4+t22*t39+t20*t38+t16*t37+t24*t36**5.0d+0+t25*t36
     1         +t12*t34+t10*t33+t14*t31**5.0d+0+t15*t31+t3+t2
            t41 = 1/t7**2
            t42 = 1/rhoa**3.6666666666666664d+0
            t43 = 1/t30**2
            t44 = 1/rhob**3.6666666666666664d+0
            t45 = 3.6802889260838756d-3*t15
            fnc(iq) = (-9.305257363491002d-1*t27*t40-9.305257363491002d-
     1         1*t1*t26)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*(-2.9442311408671
     1         007d-2*gammaaa*t10*t41*t42*t9-1.962820760578067d-2*gammaa
     2         a*t4*t41*t42*t8-9.814103802890335d-3*gammaaa*t15*t41*t42-
     3         4.9070519014451675d-2*gammaaa*t13*t14*t41*t42-3.925641521
     4         156134d-2*gammaaa*t11*t12*t41*t42-9.814103802890335d-3*ga
     5         mmaaa*t17*t25*t42-4.9070519014451675d-2*gammaaa*t17*t23*t
     6         24*t42-3.925641521156134d-2*gammaaa*t17*t21*t22*t42-2.944
     7         2311408671007d-2*gammaaa*t17*t19*t20*t42-1.96282076057806
     8         7d-2*gammaaa*t16*t17*t18*t42)-1.2407009817988002d+0*rhoa*
     9         *3.333333333333333d-1*t26)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491002d-1*t27*(-1.962820760578
     1         067d-2*gammabb*t31*t4*t43*t44-4.9070519014451675d-2*gamma
     2         bb*t14*t34*t43*t44-3.925641521156134d-2*gammabb*t12*t33*t
     3         43*t44-2.9442311408671007d-2*gammabb*t10*t32*t43*t44-9.81
     4         4103802890335d-3*gammabb*t15*t43*t44-4.9070519014451675d-
     5         2*gammabb*t24*t35*t39*t44-3.925641521156134d-2*gammabb*t2
     6         2*t35*t38*t44-2.9442311408671007d-2*gammabb*t20*t35*t37*t
     7         44-1.962820760578067d-2*gammabb*t16*t35*t36*t44-9.8141038
     8         02890335d-3*gammabb*t25*t35*t44)-1.2407009817988002d+0*rh
     9         ob**3.333333333333333d-1*t40)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         41*t5*(1.1040866778251626d-2*t10*t9+7.360577852167751d-3*
     2         t4*t8+t45+1.8401444630419378d-2*t13*t14+1.472115570433550
     3         3d-2*t11*t12)+(3.6802889260838756d-3*t17*t25+1.8401444630
     4         419378d-2*t17*t23*t24+1.4721155704335503d-2*t17*t21*t22+1
     5         .1040866778251626d-2*t17*t19*t20+7.360577852167751d-3*t16
     6         *t17*t18)*t5)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t27*(
     1         t28*t43*(t45+7.360577852167751d-3*t31*t4+1.84014446304193
     2         78d-2*t14*t34+1.4721155704335503d-2*t12*t33+1.10408667782
     3         51626d-2*t10*t32)+t28*(1.8401444630419378d-2*t24*t35*t39+
     4         1.4721155704335503d-2*t22*t35*t38+1.1040866778251626d-2*t
     5         20*t35*t37+7.360577852167751d-3*t16*t35*t36+3.68028892608
     6         38756d-3*t25*t35))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = 1/t5**2
            t26 = 1/rhoa**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RA) = -9.305257363491002d-1*t1*(-3.92564152115613
     1         4d-2*gammaaa*t10*t25*t26*t9-2.9442311408671007d-2*gammaaa
     2         *t25*t26*t7*t8-1.962820760578067d-2*gammaaa*t2*t25*t26*t6
     3         -9.814103802890335d-3*gammaaa*t13*t25*t26-4.9070519014451
     4         675d-2*gammaaa*t11*t12*t25*t26-9.814103802890335d-3*gamma
     5         aa*t15*t23*t26-4.9070519014451675d-2*gammaaa*t15*t21*t22*
     6         t26-3.925641521156134d-2*gammaaa*t15*t19*t20*t26-2.944231
     7         1408671007d-2*gammaaa*t15*t17*t18*t26-1.962820760578067d-
     8         2*gammaaa*t14*t15*t16*t26)*wght-1.2407009817988002d+0*rho
     9         a**3.333333333333333d-1*t24*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         25*t3*(1.4721155704335503d-2*t10*t9+1.1040866778251626d-2
     2         *t7*t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t
     3         13+1.8401444630419378d-2*t11*t12)+(3.6802889260838756d-3*
     4         t15*t23+1.8401444630419378d-2*t15*t21*t22+1.4721155704335
     5         503d-2*t15*t19*t20+1.1040866778251626d-2*t15*t17*t18+7.36
     6         0577852167751d-3*t14*t15*t16)*t3)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammabb*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = 1/t5**2
            t26 = 1/rhob**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RB) = -9.305257363491002d-1*t1*(-3.92564152115613
     1         4d-2*gammabb*t10*t25*t26*t9-2.9442311408671007d-2*gammabb
     2         *t25*t26*t7*t8-1.962820760578067d-2*gammabb*t2*t25*t26*t6
     3         -9.814103802890335d-3*gammabb*t13*t25*t26-4.9070519014451
     4         675d-2*gammabb*t11*t12*t25*t26-9.814103802890335d-3*gamma
     5         bb*t15*t23*t26-4.9070519014451675d-2*gammabb*t15*t21*t22*
     6         t26-3.925641521156134d-2*gammabb*t15*t19*t20*t26-2.944231
     7         1408671007d-2*gammabb*t15*t17*t18*t26-1.962820760578067d-
     8         2*gammabb*t14*t15*t16*t26)*wght-1.2407009817988002d+0*rho
     9         b**3.333333333333333d-1*t24*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t1*(t
     1         25*t3*(1.4721155704335503d-2*t10*t9+1.1040866778251626d-2
     2         *t7*t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t
     3         13+1.8401444630419378d-2*t11*t12)+(3.6802889260838756d-3*
     4         t15*t23+1.8401444630419378d-2*t15*t21*t22+1.4721155704335
     5         503d-2*t15*t19*t20+1.1040866778251626d-2*t15*t17*t18+7.36
     6         0577852167751d-3*t14*t15*t16)*t3)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_sogga functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_4} &=& {{0.003680288926083876\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{\it t_4}+1.0}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{11}} &=& 1.0-{{1}\over{e^{{\it t_4}}}}\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{16}} &=& {{0.003680288926083876\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}
C>    \over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{
C>    \it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{
C>    \it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}
C>    +{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}
C>    +{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}
C>    +{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}
C>    +{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{
C>    \it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{
C>    \it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}
C>    +{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}
C>    +{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}
C>    +{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_sogga_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhoa**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhoa**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammaaa*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammaaa*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammaaa*t2*t26*t27*t6-9.814103802890335d-3*gammaaa*t13*t26
     3         *t27-4.9070519014451675d-2*gammaaa*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammaaa*t15*t23*t27-4.9070519014451675d-2*
     5         gammaaa*t15*t21*t22*t27-3.925641521156134d-2*gammaaa*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammaaa*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammaaa*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = gammaaa**2
            t31 = 1/t5**4
            t32 = 1/rhoa**7.333333333333333d+0
            t33 = 1/t5**3
            t34 = 1/rhoa**4.666666666666667d+0
            t35 = exp(-7.360577852167751d-3*gammaaa*t3)
            t36 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t24*wght
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*t28-1.24070098179
     1         88002d+0*t24*t25)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-9.305257363491002d-1*t1*(1.4394018910
     1         905823d-1*gammaaa*t10*t26*t34*t9-7.705330676312523d-4*t10
     2         *t30*t32*t33*t9+1.9263326690781307d-3*t12*t30*t31*t32*t9+
     3         1.0795514183179368d-1*gammaaa*t26*t34*t7*t8-5.77899800723
     4         4393d-4*t30*t32*t33*t7*t8+5.778998007234393d-4*t30*t31*t3
     5         2*t6*t8+1.1557996014468785d-3*t10*t30*t31*t32*t7+7.197009
     6         455452912d-2*gammaaa*t2*t26*t34*t6-3.8526653381562614d-4*
     7         t2*t30*t32*t33*t6+1.9263326690781307d-3*t19*t22*t30*t32*t
     8         35+1.1557996014468785d-3*t17*t20*t30*t32*t35+5.7789980072
     9         34393d-4*t16*t18*t30*t32*t35+1.9263326690781307d-4*t14*t3
     :         0*t32*t35+3.598504727726456d-2*gammaaa*t13*t26*t34+1.7992
     ;         52363863228d-1*gammaaa*t11*t12*t26*t34+3.598504727726456d
     <         -2*gammaaa*t15*t23*t34+1.799252363863228d-1*gammaaa*t15*t
     =         21*t22*t34+1.4394018910905823d-1*gammaaa*t15*t19*t20*t34+
     >         1.0795514183179368d-1*gammaaa*t15*t17*t18*t34+7.197009455
     ?         452912d-2*gammaaa*t14*t15*t16*t34-1.9263326690781307d-4*t
     @         13*t30*t32*t33-9.631663345390654d-4*t11*t12*t30*t32*t33+1
     1         .9263326690781307d-4*t2*t30*t31*t32-9.631663345390653d-5*
     2         t15*t23*t30*t32-4.815831672695327d-4*t15*t21*t22*t30*t32-
     3         3.8526653381562614d-4*t15*t19*t20*t30*t32-2.8894990036171
     4         964d-4*t15*t17*t18*t30*t32-1.9263326690781307d-4*t14*t15*
     5         t16*t30*t32)-2.4814019635976003d+0*t25*t28-4.135669939329
     6         334d-1*t24/rhoa**6.666666666666666d-1)*wght+Amat2(iq,D2_R
     7         A_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491002d-1*t1*(2.889499003
     1         6171964d-4*gammaaa*t10*t33*t36*t9-7.22374750904299d-4*gam
     2         maaa*t12*t31*t36*t9-3.925641521156134d-2*t10*t26*t27*t9+2
     3         .167124252712897d-4*gammaaa*t33*t36*t7*t8-2.9442311408671
     4         007d-2*t26*t27*t7*t8-2.167124252712897d-4*gammaaa*t31*t36
     5         *t6*t8-4.334248505425794d-4*gammaaa*t10*t31*t36*t7+1.4447
     6         495018085982d-4*gammaaa*t2*t33*t36*t6-1.962820760578067d-
     7         2*t2*t26*t27*t6-7.22374750904299d-4*gammaaa*t19*t22*t35*t
     8         36-4.334248505425794d-4*gammaaa*t17*t20*t35*t36-2.1671242
     9         52712897d-4*gammaaa*t16*t18*t35*t36-7.22374750904299d-5*g
     :         ammaaa*t14*t35*t36+7.22374750904299d-5*gammaaa*t13*t33*t3
     ;         6+3.611873754521495d-4*gammaaa*t11*t12*t33*t36-7.22374750
     <         904299d-5*gammaaa*t2*t31*t36+3.611873754521495d-5*gammaaa
     =         *t15*t23*t36+1.8059368772607476d-4*gammaaa*t15*t21*t22*t3
     >         6+1.4447495018085982d-4*gammaaa*t15*t19*t20*t36+1.0835621
     ?         263564485d-4*gammaaa*t15*t17*t18*t36+7.22374750904299d-5*
     @         gammaaa*t14*t15*t16*t36-9.814103802890335d-3*t13*t26*t27-
     1         4.9070519014451675d-2*t11*t12*t26*t27-9.814103802890335d-
     2         3*t15*t23*t27-4.9070519014451675d-2*t15*t21*t22*t27-3.925
     3         641521156134d-2*t15*t19*t20*t27-2.9442311408671007d-2*t15
     4         *t17*t18*t27-1.962820760578067d-2*t14*t15*t16*t27)-1.2407
     5         009817988002d+0*t25*(1.4721155704335503d-2*t10*t26*t3*t9+
     6         1.1040866778251626d-2*t26*t3*t7*t8+7.360577852167751d-3*t
     7         2*t26*t3*t6+3.6802889260838756d-3*t13*t26*t3+1.8401444630
     8         419378d-2*t11*t12*t26*t3+3.6802889260838756d-3*t15*t23*t3
     9         +1.8401444630419378d-2*t15*t21*t22*t3+1.4721155704335503d
     :         -2*t15*t19*t20*t3+1.1040866778251626d-2*t15*t17*t18*t3+7.
     ;         360577852167751d-3*t14*t15*t16*t3))*wght+Cmat2(iq,D2_RA_G
     <         AA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t26*t3*(2.708905315891122d-4*t12*t26*t3*t9+8.126
     2         715947673364d-5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26
     3         *t3*t7+2.7089053158911214d-5*t2*t26*t3)+t3*(2.70890531589
     4         1122d-4*t19*t22*t3*t35+1.625343189534673d-4*t17*t20*t3*t3
     5         5+8.126715947673364d-5*t16*t18*t3*t35+2.7089053158911214d
     6         -5*t14*t3*t35-1.3544526579455607d-5*t15*t23*t3-6.77226328
     7         9727805d-5*t15*t21*t22*t3-5.417810631782243d-5*t15*t19*t2
     8         0*t3-4.063357973836682d-5*t15*t17*t18*t3-2.70890531589112
     9         14d-5*t14*t15*t16*t3)-7.360577852167751d-3*t29*t33/rhoa**
     :         5.333333333333333d+0)*wght
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
            t2 = param(1)
            t3 = param(7)
            t4 = param(3)
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 3.6802889260838756d-3*gammaaa*t5
            t7 = t6+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = t8**2.0d+0
            t10 = param(4)
            t11 = t8**3.0d+0
            t12 = param(5)
            t13 = t8**4.0d+0
            t14 = param(6)
            t15 = param(2)
            t16 = param(9)
            t17 = exp(-t6)
            t18 = 1.0d+0-t17
            t19 = t18**2.0d+0
            t20 = param(10)
            t21 = t18**3.0d+0
            t22 = param(11)
            t23 = t18**4.0d+0
            t24 = param(12)
            t25 = param(8)
            t26 = t4*t9+t14*t8**5.0d+0+t15*t8+t3+t18*t25+t18**5.0d+0*t24
     1         +t22*t23+t20*t21+t2+t16*t19+t12*t13+t10*t11
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/rhob**2.6666666666666666d+0
            t29 = 3.6802889260838756d-3*gammabb*t28
            t30 = t29+1.0d+0
            t31 = 1.0d+0-1.0d+0/t30
            t32 = t31**2.0d+0
            t33 = t31**3.0d+0
            t34 = t31**4.0d+0
            t35 = exp(-t29)
            t36 = 1.0d+0-t35
            t37 = t36**2.0d+0
            t38 = t36**3.0d+0
            t39 = t36**4.0d+0
            t40 = t32*t4+t22*t39+t20*t38+t16*t37+t24*t36**5.0d+0+t25*t36
     1         +t12*t34+t10*t33+t14*t31**5.0d+0+t15*t31+t3+t2
            t41 = rhoa**3.333333333333333d-1
            t42 = 1/t7**2
            t43 = 1/rhoa**3.6666666666666664d+0
            t44 = -2.9442311408671007d-2*gammaaa*t10*t42*t43*t9-1.962820
     1         760578067d-2*gammaaa*t4*t42*t43*t8-9.814103802890335d-3*g
     2         ammaaa*t15*t42*t43-4.9070519014451675d-2*gammaaa*t13*t14*
     3         t42*t43-3.925641521156134d-2*gammaaa*t11*t12*t42*t43-9.81
     4         4103802890335d-3*gammaaa*t17*t25*t43-4.9070519014451675d-
     5         2*gammaaa*t17*t23*t24*t43-3.925641521156134d-2*gammaaa*t1
     6         7*t21*t22*t43-2.9442311408671007d-2*gammaaa*t17*t19*t20*t
     7         43-1.962820760578067d-2*gammaaa*t16*t17*t18*t43
            t45 = rhob**3.333333333333333d-1
            t46 = 1/t30**2
            t47 = 1/rhob**3.6666666666666664d+0
            t48 = -1.962820760578067d-2*gammabb*t31*t4*t46*t47-4.9070519
     1         014451675d-2*gammabb*t14*t34*t46*t47-3.925641521156134d-2
     2         *gammabb*t12*t33*t46*t47-2.9442311408671007d-2*gammabb*t1
     3         0*t32*t46*t47-9.814103802890335d-3*gammabb*t15*t46*t47-4.
     4         9070519014451675d-2*gammabb*t24*t35*t39*t47-3.92564152115
     5         6134d-2*gammabb*t22*t35*t38*t47-2.9442311408671007d-2*gam
     6         mabb*t20*t35*t37*t47-1.962820760578067d-2*gammabb*t16*t35
     7         *t36*t47-9.814103802890335d-3*gammabb*t25*t35*t47
            t49 = 3.6802889260838756d-3*t15
            t50 = 1.1040866778251626d-2*t10*t9+7.360577852167751d-3*t4*t
     1         8+t49+1.8401444630419378d-2*t13*t14+1.4721155704335503d-2
     2         *t11*t12
            t51 = t49+7.360577852167751d-3*t31*t4+1.8401444630419378d-2*
     1         t14*t34+1.4721155704335503d-2*t12*t33+1.1040866778251626d
     2         -2*t10*t32
            t52 = gammaaa**2
            t53 = 1/t7**4
            t54 = 1/rhoa**7.333333333333333d+0
            t55 = 1/t7**3
            t56 = 1/rhoa**4.666666666666667d+0
            t57 = exp(-7.360577852167751d-3*gammaaa*t5)
            t58 = gammabb**2
            t59 = 1/t30**4
            t60 = 1/rhob**7.333333333333333d+0
            t61 = 1/t30**3
            t62 = 1/rhob**4.666666666666667d+0
            t63 = exp(-7.360577852167751d-3*gammabb*t28)
            t64 = 1/rhoa**6.333333333333333d+0
            t65 = 1/rhob**6.333333333333333d+0
            fnc(iq) = (-9.305257363491002d-1*t27*t40-9.305257363491002d-
     1         1*t1*t26)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*t44-1.24070098179
     1         88002d+0*t26*t41)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491002d-1*t27*t48-1.2407009817
     1         988002d+0*t40*t45)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         42*t5*t50+(3.6802889260838756d-3*t17*t25+1.84014446304193
     2         78d-2*t17*t23*t24+1.4721155704335503d-2*t17*t21*t22+1.104
     3         0866778251626d-2*t17*t19*t20+7.360577852167751d-3*t16*t17
     4         *t18)*t5)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t27*(
     1         t28*t46*t51+t28*(1.8401444630419378d-2*t24*t35*t39+1.4721
     2         155704335503d-2*t22*t35*t38+1.1040866778251626d-2*t20*t35
     3         *t37+7.360577852167751d-3*t16*t35*t36+3.6802889260838756d
     4         -3*t25*t35))*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491002d-1*t1*(1.0795514183
     1         179368d-1*gammaaa*t10*t42*t56*t9-5.778998007234393d-4*t10
     2         *t52*t54*t55*t9+1.1557996014468785d-3*t12*t52*t53*t54*t9+
     3         7.197009455452912d-2*gammaaa*t4*t42*t56*t8-3.852665338156
     4         2614d-4*t4*t52*t54*t55*t8+5.778998007234393d-4*t10*t52*t5
     5         3*t54*t8+1.9263326690781307d-3*t21*t24*t52*t54*t57+1.1557
     6         996014468785d-3*t19*t22*t52*t54*t57+5.778998007234393d-4*
     7         t18*t20*t52*t54*t57+1.9263326690781307d-4*t16*t52*t54*t57
     8         +3.598504727726456d-2*gammaaa*t15*t42*t56+1.7992523638632
     9         28d-1*gammaaa*t13*t14*t42*t56+1.4394018910905823d-1*gamma
     :         aa*t11*t12*t42*t56+3.598504727726456d-2*gammaaa*t17*t25*t
     ;         56+1.799252363863228d-1*gammaaa*t17*t23*t24*t56+1.4394018
     <         910905823d-1*gammaaa*t17*t21*t22*t56+1.0795514183179368d-
     =         1*gammaaa*t17*t19*t20*t56+7.197009455452912d-2*gammaaa*t1
     >         6*t17*t18*t56-1.9263326690781307d-4*t15*t52*t54*t55-9.631
     ?         663345390654d-4*t13*t14*t52*t54*t55-7.705330676312523d-4*
     @         t11*t12*t52*t54*t55+1.9263326690781307d-4*t4*t52*t53*t54+
     1         1.9263326690781307d-3*t11*t14*t52*t53*t54-9.6316633453906
     2         53d-5*t17*t25*t52*t54-4.815831672695327d-4*t17*t23*t24*t5
     3         2*t54-3.8526653381562614d-4*t17*t21*t22*t52*t54-2.8894990
     4         036171964d-4*t17*t19*t20*t52*t54-1.9263326690781307d-4*t1
     5         6*t17*t18*t52*t54)-2.4814019635976003d+0*t41*t44-4.135669
     6         939329334d-1*t26/rhoa**6.666666666666666d-1)*wght+Amat2(i
     7         q,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491002d-1*t27*(1.926332669
     1         0781307d-3*t24*t38*t58*t60*t63+1.1557996014468785d-3*t22*
     2         t37*t58*t60*t63+5.778998007234393d-4*t20*t36*t58*t60*t63+
     3         1.9263326690781307d-4*t16*t58*t60*t63+7.197009455452912d-
     4         2*gammabb*t31*t4*t46*t62+1.799252363863228d-1*gammabb*t14
     5         *t34*t46*t62+1.4394018910905823d-1*gammabb*t12*t33*t46*t6
     6         2+1.0795514183179368d-1*gammabb*t10*t32*t46*t62+3.5985047
     7         27726456d-2*gammabb*t15*t46*t62+1.799252363863228d-1*gamm
     8         abb*t24*t35*t39*t62+1.4394018910905823d-1*gammabb*t22*t35
     9         *t38*t62+1.0795514183179368d-1*gammabb*t20*t35*t37*t62+7.
     :         197009455452912d-2*gammabb*t16*t35*t36*t62+3.598504727726
     ;         456d-2*gammabb*t25*t35*t62-3.8526653381562614d-4*t31*t4*t
     <         58*t60*t61-9.631663345390654d-4*t14*t34*t58*t60*t61-7.705
     =         330676312523d-4*t12*t33*t58*t60*t61-5.778998007234393d-4*
     >         t10*t32*t58*t60*t61-1.9263326690781307d-4*t15*t58*t60*t61
     ?         +1.9263326690781307d-4*t4*t58*t59*t60+1.9263326690781307d
     @         -3*t14*t33*t58*t59*t60+1.1557996014468785d-3*t12*t32*t58*
     1         t59*t60+5.778998007234393d-4*t10*t31*t58*t59*t60-4.815831
     2         672695327d-4*t24*t35*t39*t58*t60-3.8526653381562614d-4*t2
     3         2*t35*t38*t58*t60-2.8894990036171964d-4*t20*t35*t37*t58*t
     4         60-1.9263326690781307d-4*t16*t35*t36*t58*t60-9.6316633453
     5         90653d-5*t25*t35*t58*t60)-2.4814019635976003d+0*t45*t48-4
     6         .135669939329334d-1*t40/rhob**6.666666666666666d-1)*wght+
     7         Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491002d-1*t1*(2.167124252
     1         712897d-4*gammaaa*t10*t55*t64*t9-4.334248505425794d-4*gam
     2         maaa*t12*t53*t64*t9-2.9442311408671007d-2*t10*t42*t43*t9+
     3         1.4447495018085982d-4*gammaaa*t4*t55*t64*t8-2.16712425271
     4         2897d-4*gammaaa*t10*t53*t64*t8-1.962820760578067d-2*t4*t4
     5         2*t43*t8-7.22374750904299d-4*gammaaa*t21*t24*t57*t64-4.33
     6         4248505425794d-4*gammaaa*t19*t22*t57*t64-2.16712425271289
     7         7d-4*gammaaa*t18*t20*t57*t64-7.22374750904299d-5*gammaaa*
     8         t16*t57*t64+7.22374750904299d-5*gammaaa*t15*t55*t64+3.611
     9         873754521495d-4*gammaaa*t13*t14*t55*t64+2.889499003617196
     :         4d-4*gammaaa*t11*t12*t55*t64-7.22374750904299d-5*gammaaa*
     ;         t4*t53*t64-7.22374750904299d-4*gammaaa*t11*t14*t53*t64+3.
     <         611873754521495d-5*gammaaa*t17*t25*t64+1.8059368772607476
     =         d-4*gammaaa*t17*t23*t24*t64+1.4447495018085982d-4*gammaaa
     >         *t17*t21*t22*t64+1.0835621263564485d-4*gammaaa*t17*t19*t2
     ?         0*t64+7.22374750904299d-5*gammaaa*t16*t17*t18*t64-9.81410
     @         3802890335d-3*t15*t42*t43-4.9070519014451675d-2*t13*t14*t
     1         42*t43-3.925641521156134d-2*t11*t12*t42*t43-9.81410380289
     2         0335d-3*t17*t25*t43-4.9070519014451675d-2*t17*t23*t24*t43
     3         -3.925641521156134d-2*t17*t21*t22*t43-2.9442311408671007d
     4         -2*t17*t19*t20*t43-1.962820760578067d-2*t16*t17*t18*t43)-
     5         1.2407009817988002d+0*t41*(1.1040866778251626d-2*t10*t42*
     6         t5*t9+7.360577852167751d-3*t4*t42*t5*t8+3.680288926083875
     7         6d-3*t15*t42*t5+1.8401444630419378d-2*t13*t14*t42*t5+1.47
     8         21155704335503d-2*t11*t12*t42*t5+3.6802889260838756d-3*t1
     9         7*t25*t5+1.8401444630419378d-2*t17*t23*t24*t5+1.472115570
     :         4335503d-2*t17*t21*t22*t5+1.1040866778251626d-2*t17*t19*t
     ;         20*t5+7.360577852167751d-3*t16*t17*t18*t5))*wght+Cmat2(iq
     <         ,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363491002d-1*t27*(-7.2237475
     1         0904299d-4*gammabb*t24*t38*t63*t65-4.334248505425794d-4*g
     2         ammabb*t22*t37*t63*t65-2.167124252712897d-4*gammabb*t20*t
     3         36*t63*t65-7.22374750904299d-5*gammabb*t16*t63*t65+1.4447
     4         495018085982d-4*gammabb*t31*t4*t61*t65+3.611873754521495d
     5         -4*gammabb*t14*t34*t61*t65+2.8894990036171964d-4*gammabb*
     6         t12*t33*t61*t65+2.167124252712897d-4*gammabb*t10*t32*t61*
     7         t65+7.22374750904299d-5*gammabb*t15*t61*t65-7.22374750904
     8         299d-5*gammabb*t4*t59*t65-7.22374750904299d-4*gammabb*t14
     9         *t33*t59*t65-4.334248505425794d-4*gammabb*t12*t32*t59*t65
     :         -2.167124252712897d-4*gammabb*t10*t31*t59*t65+1.805936877
     ;         2607476d-4*gammabb*t24*t35*t39*t65+1.4447495018085982d-4*
     <         gammabb*t22*t35*t38*t65+1.0835621263564485d-4*gammabb*t20
     =         *t35*t37*t65+7.22374750904299d-5*gammabb*t16*t35*t36*t65+
     >         3.611873754521495d-5*gammabb*t25*t35*t65-1.96282076057806
     ?         7d-2*t31*t4*t46*t47-4.9070519014451675d-2*t14*t34*t46*t47
     @         -3.925641521156134d-2*t12*t33*t46*t47-2.9442311408671007d
     1         -2*t10*t32*t46*t47-9.814103802890335d-3*t15*t46*t47-4.907
     2         0519014451675d-2*t24*t35*t39*t47-3.925641521156134d-2*t22
     3         *t35*t38*t47-2.9442311408671007d-2*t20*t35*t37*t47-1.9628
     4         20760578067d-2*t16*t35*t36*t47-9.814103802890335d-3*t25*t
     5         35*t47)-1.2407009817988002d+0*t45*(7.360577852167751d-3*t
     6         28*t31*t4*t46+1.8401444630419378d-2*t14*t28*t34*t46+1.472
     7         1155704335503d-2*t12*t28*t33*t46+1.1040866778251626d-2*t1
     8         0*t28*t32*t46+3.6802889260838756d-3*t15*t28*t46+1.8401444
     9         630419378d-2*t24*t28*t35*t39+1.4721155704335503d-2*t22*t2
     :         8*t35*t38+1.1040866778251626d-2*t20*t28*t35*t37+7.3605778
     ;         52167751d-3*t16*t28*t35*t36+3.6802889260838756d-3*t25*t28
     <         *t35))*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t42*t5*(1.625343189534673d-4*t12*t42*t5*t9+8.126
     2         715947673364d-5*t10*t42*t5*t8+2.7089053158911214d-5*t4*t4
     3         2*t5+2.708905315891122d-4*t11*t14*t42*t5)+t5*(2.708905315
     4         891122d-4*t21*t24*t5*t57+1.625343189534673d-4*t19*t22*t5*
     5         t57+8.126715947673364d-5*t18*t20*t5*t57+2.708905315891121
     6         4d-5*t16*t5*t57-1.3544526579455607d-5*t17*t25*t5-6.772263
     7         289727805d-5*t17*t23*t24*t5-5.417810631782243d-5*t17*t21*
     8         t22*t5-4.063357973836682d-5*t17*t19*t20*t5-2.708905315891
     9         1214d-5*t16*t17*t18*t5)-7.360577852167751d-3*t50*t55/rhoa
     :         **5.333333333333333d+0)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349100
     1         2d-1*t27*(t28*(2.708905315891122d-4*t24*t28*t38*t63+1.625
     2         343189534673d-4*t22*t28*t37*t63+8.126715947673364d-5*t20*
     3         t28*t36*t63+2.7089053158911214d-5*t16*t28*t63-6.772263289
     4         727805d-5*t24*t28*t35*t39-5.417810631782243d-5*t22*t28*t3
     5         5*t38-4.063357973836682d-5*t20*t28*t35*t37-2.708905315891
     6         1214d-5*t16*t28*t35*t36-1.3544526579455607d-5*t25*t28*t35
     7         )-7.360577852167751d-3*t51*t61/rhob**5.333333333333333d+0
     8         +t28*t46*(2.7089053158911214d-5*t28*t4*t46+2.708905315891
     9         122d-4*t14*t28*t33*t46+1.625343189534673d-4*t12*t28*t32*t
     :         46+8.126715947673364d-5*t10*t28*t31*t46))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhoa**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhoa**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammaaa*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammaaa*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammaaa*t2*t26*t27*t6-9.814103802890335d-3*gammaaa*t13*t26
     3         *t27-4.9070519014451675d-2*gammaaa*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammaaa*t15*t23*t27-4.9070519014451675d-2*
     5         gammaaa*t15*t21*t22*t27-3.925641521156134d-2*gammaaa*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammaaa*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammaaa*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = gammaaa**2
            t31 = 1/t5**4
            t32 = 1/rhoa**7.333333333333333d+0
            t33 = 1/t5**3
            t34 = 1/rhoa**4.666666666666667d+0
            t35 = exp(-7.360577852167751d-3*gammaaa*t3)
            t36 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RA) = -9.305257363491002d-1*t1*t28*wght-1.2407009
     1         817988002d+0*t24*t25*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Amat2(iq,D2_RA_RA) = -9.305257363491002d-1*t1*(1.43940189109
     1         05823d-1*gammaaa*t10*t26*t34*t9-7.705330676312523d-4*t10*
     2         t30*t32*t33*t9+1.9263326690781307d-3*t12*t30*t31*t32*t9+1
     3         .0795514183179368d-1*gammaaa*t26*t34*t7*t8-5.778998007234
     4         393d-4*t30*t32*t33*t7*t8+5.778998007234393d-4*t30*t31*t32
     5         *t6*t8+1.1557996014468785d-3*t10*t30*t31*t32*t7+7.1970094
     6         55452912d-2*gammaaa*t2*t26*t34*t6-3.8526653381562614d-4*t
     7         2*t30*t32*t33*t6+1.9263326690781307d-3*t19*t22*t30*t32*t3
     8         5+1.1557996014468785d-3*t17*t20*t30*t32*t35+5.77899800723
     9         4393d-4*t16*t18*t30*t32*t35+1.9263326690781307d-4*t14*t30
     :         *t32*t35+3.598504727726456d-2*gammaaa*t13*t26*t34+1.79925
     ;         2363863228d-1*gammaaa*t11*t12*t26*t34+3.598504727726456d-
     <         2*gammaaa*t15*t23*t34+1.799252363863228d-1*gammaaa*t15*t2
     =         1*t22*t34+1.4394018910905823d-1*gammaaa*t15*t19*t20*t34+1
     >         .0795514183179368d-1*gammaaa*t15*t17*t18*t34+7.1970094554
     ?         52912d-2*gammaaa*t14*t15*t16*t34-1.9263326690781307d-4*t1
     @         3*t30*t32*t33-9.631663345390654d-4*t11*t12*t30*t32*t33+1.
     1         9263326690781307d-4*t2*t30*t31*t32-9.631663345390653d-5*t
     2         15*t23*t30*t32-4.815831672695327d-4*t15*t21*t22*t30*t32-3
     3         .8526653381562614d-4*t15*t19*t20*t30*t32-2.88949900361719
     4         64d-4*t15*t17*t18*t30*t32-1.9263326690781307d-4*t14*t15*t
     5         16*t30*t32)*wght-2.4814019635976003d+0*t25*t28*wght-4.135
     6         669939329334d-1*t24*wght/rhoa**6.666666666666666d-1+Amat2
     7         (iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363491002d-1*t1*(2.8894990036
     1         171964d-4*gammaaa*t10*t33*t36*t9-7.22374750904299d-4*gamm
     2         aaa*t12*t31*t36*t9-3.925641521156134d-2*t10*t26*t27*t9+2.
     3         167124252712897d-4*gammaaa*t33*t36*t7*t8-2.94423114086710
     4         07d-2*t26*t27*t7*t8-2.167124252712897d-4*gammaaa*t31*t36*
     5         t6*t8-4.334248505425794d-4*gammaaa*t10*t31*t36*t7+1.44474
     6         95018085982d-4*gammaaa*t2*t33*t36*t6-1.962820760578067d-2
     7         *t2*t26*t27*t6-7.22374750904299d-4*gammaaa*t19*t22*t35*t3
     8         6-4.334248505425794d-4*gammaaa*t17*t20*t35*t36-2.16712425
     9         2712897d-4*gammaaa*t16*t18*t35*t36-7.22374750904299d-5*ga
     :         mmaaa*t14*t35*t36+7.22374750904299d-5*gammaaa*t13*t33*t36
     ;         +3.611873754521495d-4*gammaaa*t11*t12*t33*t36-7.223747509
     <         04299d-5*gammaaa*t2*t31*t36+3.611873754521495d-5*gammaaa*
     =         t15*t23*t36+1.8059368772607476d-4*gammaaa*t15*t21*t22*t36
     >         +1.4447495018085982d-4*gammaaa*t15*t19*t20*t36+1.08356212
     ?         63564485d-4*gammaaa*t15*t17*t18*t36+7.22374750904299d-5*g
     @         ammaaa*t14*t15*t16*t36-9.814103802890335d-3*t13*t26*t27-4
     1         .9070519014451675d-2*t11*t12*t26*t27-9.814103802890335d-3
     2         *t15*t23*t27-4.9070519014451675d-2*t15*t21*t22*t27-3.9256
     3         41521156134d-2*t15*t19*t20*t27-2.9442311408671007d-2*t15*
     4         t17*t18*t27-1.962820760578067d-2*t14*t15*t16*t27)*wght-1.
     5         2407009817988002d+0*t25*(1.4721155704335503d-2*t10*t26*t3
     6         *t9+1.1040866778251626d-2*t26*t3*t7*t8+7.360577852167751d
     7         -3*t2*t26*t3*t6+3.6802889260838756d-3*t13*t26*t3+1.840144
     8         4630419378d-2*t11*t12*t26*t3+3.6802889260838756d-3*t15*t2
     9         3*t3+1.8401444630419378d-2*t15*t21*t22*t3+1.4721155704335
     :         503d-2*t15*t19*t20*t3+1.1040866778251626d-2*t15*t17*t18*t
     ;         3+7.360577852167751d-3*t14*t15*t16*t3)*wght+Cmat2(iq,D2_R
     <         A_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t26*t3*(2.708905315891122d-4*t12*t26*t3*t9+8.126
     2         715947673364d-5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26
     3         *t3*t7+2.7089053158911214d-5*t2*t26*t3)+t3*(2.70890531589
     4         1122d-4*t19*t22*t3*t35+1.625343189534673d-4*t17*t20*t3*t3
     5         5+8.126715947673364d-5*t16*t18*t3*t35+2.7089053158911214d
     6         -5*t14*t3*t35-1.3544526579455607d-5*t15*t23*t3-6.77226328
     7         9727805d-5*t15*t21*t22*t3-5.417810631782243d-5*t15*t19*t2
     8         0*t3-4.063357973836682d-5*t15*t17*t18*t3-2.70890531589112
     9         14d-5*t14*t15*t16*t3)-7.360577852167751d-3*t29*t33/rhoa**
     :         5.333333333333333d+0)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammabb*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhob**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhob**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammabb*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammabb*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammabb*t2*t26*t27*t6-9.814103802890335d-3*gammabb*t13*t26
     3         *t27-4.9070519014451675d-2*gammabb*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammabb*t15*t23*t27-4.9070519014451675d-2*
     5         gammabb*t15*t21*t22*t27-3.925641521156134d-2*gammabb*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammabb*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammabb*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = gammabb**2
            t31 = 1/t5**4
            t32 = 1/rhob**7.333333333333333d+0
            t33 = 1/t5**3
            t34 = 1/rhob**4.666666666666667d+0
            t35 = exp(-7.360577852167751d-3*gammabb*t3)
            t36 = 1/rhob**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RB) = -9.305257363491002d-1*t1*t28*wght-1.2407009
     1         817988002d+0*t24*t25*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Amat2(iq,D2_RB_RB) = -9.305257363491002d-1*t1*(1.43940189109
     1         05823d-1*gammabb*t10*t26*t34*t9-7.705330676312523d-4*t10*
     2         t30*t32*t33*t9+1.9263326690781307d-3*t12*t30*t31*t32*t9+1
     3         .0795514183179368d-1*gammabb*t26*t34*t7*t8-5.778998007234
     4         393d-4*t30*t32*t33*t7*t8+5.778998007234393d-4*t30*t31*t32
     5         *t6*t8+1.1557996014468785d-3*t10*t30*t31*t32*t7+7.1970094
     6         55452912d-2*gammabb*t2*t26*t34*t6-3.8526653381562614d-4*t
     7         2*t30*t32*t33*t6+1.9263326690781307d-3*t19*t22*t30*t32*t3
     8         5+1.1557996014468785d-3*t17*t20*t30*t32*t35+5.77899800723
     9         4393d-4*t16*t18*t30*t32*t35+1.9263326690781307d-4*t14*t30
     :         *t32*t35+3.598504727726456d-2*gammabb*t13*t26*t34+1.79925
     ;         2363863228d-1*gammabb*t11*t12*t26*t34+3.598504727726456d-
     <         2*gammabb*t15*t23*t34+1.799252363863228d-1*gammabb*t15*t2
     =         1*t22*t34+1.4394018910905823d-1*gammabb*t15*t19*t20*t34+1
     >         .0795514183179368d-1*gammabb*t15*t17*t18*t34+7.1970094554
     ?         52912d-2*gammabb*t14*t15*t16*t34-1.9263326690781307d-4*t1
     @         3*t30*t32*t33-9.631663345390654d-4*t11*t12*t30*t32*t33+1.
     1         9263326690781307d-4*t2*t30*t31*t32-9.631663345390653d-5*t
     2         15*t23*t30*t32-4.815831672695327d-4*t15*t21*t22*t30*t32-3
     3         .8526653381562614d-4*t15*t19*t20*t30*t32-2.88949900361719
     4         64d-4*t15*t17*t18*t30*t32-1.9263326690781307d-4*t14*t15*t
     5         16*t30*t32)*wght-2.4814019635976003d+0*t25*t28*wght-4.135
     6         669939329334d-1*t24*wght/rhob**6.666666666666666d-1+Amat2
     7         (iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363491002d-1*t1*(2.8894990036
     1         171964d-4*gammabb*t10*t33*t36*t9-7.22374750904299d-4*gamm
     2         abb*t12*t31*t36*t9-3.925641521156134d-2*t10*t26*t27*t9+2.
     3         167124252712897d-4*gammabb*t33*t36*t7*t8-2.94423114086710
     4         07d-2*t26*t27*t7*t8-2.167124252712897d-4*gammabb*t31*t36*
     5         t6*t8-4.334248505425794d-4*gammabb*t10*t31*t36*t7+1.44474
     6         95018085982d-4*gammabb*t2*t33*t36*t6-1.962820760578067d-2
     7         *t2*t26*t27*t6-7.22374750904299d-4*gammabb*t19*t22*t35*t3
     8         6-4.334248505425794d-4*gammabb*t17*t20*t35*t36-2.16712425
     9         2712897d-4*gammabb*t16*t18*t35*t36-7.22374750904299d-5*ga
     :         mmabb*t14*t35*t36+7.22374750904299d-5*gammabb*t13*t33*t36
     ;         +3.611873754521495d-4*gammabb*t11*t12*t33*t36-7.223747509
     <         04299d-5*gammabb*t2*t31*t36+3.611873754521495d-5*gammabb*
     =         t15*t23*t36+1.8059368772607476d-4*gammabb*t15*t21*t22*t36
     >         +1.4447495018085982d-4*gammabb*t15*t19*t20*t36+1.08356212
     ?         63564485d-4*gammabb*t15*t17*t18*t36+7.22374750904299d-5*g
     @         ammabb*t14*t15*t16*t36-9.814103802890335d-3*t13*t26*t27-4
     1         .9070519014451675d-2*t11*t12*t26*t27-9.814103802890335d-3
     2         *t15*t23*t27-4.9070519014451675d-2*t15*t21*t22*t27-3.9256
     3         41521156134d-2*t15*t19*t20*t27-2.9442311408671007d-2*t15*
     4         t17*t18*t27-1.962820760578067d-2*t14*t15*t16*t27)*wght-1.
     5         2407009817988002d+0*t25*(1.4721155704335503d-2*t10*t26*t3
     6         *t9+1.1040866778251626d-2*t26*t3*t7*t8+7.360577852167751d
     7         -3*t2*t26*t3*t6+3.6802889260838756d-3*t13*t26*t3+1.840144
     8         4630419378d-2*t11*t12*t26*t3+3.6802889260838756d-3*t15*t2
     9         3*t3+1.8401444630419378d-2*t15*t21*t22*t3+1.4721155704335
     :         503d-2*t15*t19*t20*t3+1.1040866778251626d-2*t15*t17*t18*t
     ;         3+7.360577852167751d-3*t14*t15*t16*t3)*wght+Cmat2(iq,D2_R
     <         B_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349100
     1         2d-1*t1*(t26*t3*(2.708905315891122d-4*t12*t26*t3*t9+8.126
     2         715947673364d-5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26
     3         *t3*t7+2.7089053158911214d-5*t2*t26*t3)+t3*(2.70890531589
     4         1122d-4*t19*t22*t3*t35+1.625343189534673d-4*t17*t20*t3*t3
     5         5+8.126715947673364d-5*t16*t18*t3*t35+2.7089053158911214d
     6         -5*t14*t3*t35-1.3544526579455607d-5*t15*t23*t3-6.77226328
     7         9727805d-5*t15*t21*t22*t3-5.417810631782243d-5*t15*t19*t2
     8         0*t3-4.063357973836682d-5*t15*t17*t18*t3-2.70890531589112
     9         14d-5*t14*t15*t16*t3)-7.360577852167751d-3*t29*t33/rhob**
     :         5.333333333333333d+0)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_sogga functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_4} &=& {{0.003680288926083876\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& 1.0-{{1.0}\over{{\it t_4}+1.0}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(9\right)\\\\ 
C>   {\it t_{11}} &=& 1.0-{{1}\over{e^{{\it t_4}}}}\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(10\right)\\\\ 
C>   {\it t_{13}} &=& {\it param}\left(11\right)\\\\ 
C>   {\it t_{14}} &=& {\it param}\left(12\right)\\\\ 
C>   {\it t_{15}} &=& {\it param}\left(8\right)\\\\ 
C>   {\it t_{16}} &=& {{0.003680288926083876\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}
C>    \over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{
C>    \it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{
C>    \it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}
C>    +{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}
C>    +{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}
C>    +{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}
C>    +{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{
C>    \it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{
C>    \it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({
C>    \it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{
C>    \it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}
C>    +{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}
C>    +{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}
C>    +{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}
C>    +{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="https://doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_sogga_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhoa**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhoa**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammaaa*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammaaa*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammaaa*t2*t26*t27*t6-9.814103802890335d-3*gammaaa*t13*t26
     3         *t27-4.9070519014451675d-2*gammaaa*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammaaa*t15*t23*t27-4.9070519014451675d-2*
     5         gammaaa*t15*t21*t22*t27-3.925641521156134d-2*gammaaa*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammaaa*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammaaa*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = 1/rhoa**6.666666666666666d-1
            t31 = gammaaa**2
            t32 = 1/t5**4
            t33 = 1/rhoa**7.333333333333333d+0
            t34 = 1/t5**3
            t35 = 1/rhoa**4.666666666666667d+0
            t36 = exp(-7.360577852167751d-3*gammaaa*t3)
            t37 = 1.4394018910905823d-1*gammaaa*t10*t26*t35*t9-7.7053306
     1         76312523d-4*t10*t31*t33*t34*t9+1.9263326690781307d-3*t12*
     2         t31*t32*t33*t9+1.0795514183179368d-1*gammaaa*t26*t35*t7*t
     3         8-5.778998007234393d-4*t31*t33*t34*t7*t8+5.77899800723439
     4         3d-4*t31*t32*t33*t6*t8+1.1557996014468785d-3*t10*t31*t32*
     5         t33*t7+7.197009455452912d-2*gammaaa*t2*t26*t35*t6-3.85266
     6         53381562614d-4*t2*t31*t33*t34*t6+1.9263326690781307d-3*t1
     7         9*t22*t31*t33*t36+1.1557996014468785d-3*t17*t20*t31*t33*t
     8         36+5.778998007234393d-4*t16*t18*t31*t33*t36+1.92633266907
     9         81307d-4*t14*t31*t33*t36+3.598504727726456d-2*gammaaa*t13
     :         *t26*t35+1.799252363863228d-1*gammaaa*t11*t12*t26*t35+3.5
     ;         98504727726456d-2*gammaaa*t15*t23*t35+1.799252363863228d-
     <         1*gammaaa*t15*t21*t22*t35+1.4394018910905823d-1*gammaaa*t
     =         15*t19*t20*t35+1.0795514183179368d-1*gammaaa*t15*t17*t18*
     >         t35+7.197009455452912d-2*gammaaa*t14*t15*t16*t35-1.926332
     ?         6690781307d-4*t13*t31*t33*t34-9.631663345390654d-4*t11*t1
     @         2*t31*t33*t34+1.9263326690781307d-4*t2*t31*t32*t33-9.6316
     1         63345390653d-5*t15*t23*t31*t33-4.815831672695327d-4*t15*t
     2         21*t22*t31*t33-3.8526653381562614d-4*t15*t19*t20*t31*t33-
     3         2.8894990036171964d-4*t15*t17*t18*t31*t33-1.9263326690781
     4         307d-4*t14*t15*t16*t31*t33
            t38 = 1/rhoa**6.333333333333333d+0
            t39 = 2.8894990036171964d-4*gammaaa*t10*t34*t38*t9-7.2237475
     1         0904299d-4*gammaaa*t12*t32*t38*t9-3.925641521156134d-2*t1
     2         0*t26*t27*t9+2.167124252712897d-4*gammaaa*t34*t38*t7*t8-2
     3         .9442311408671007d-2*t26*t27*t7*t8-2.167124252712897d-4*g
     4         ammaaa*t32*t38*t6*t8-4.334248505425794d-4*gammaaa*t10*t32
     5         *t38*t7+1.4447495018085982d-4*gammaaa*t2*t34*t38*t6-1.962
     6         820760578067d-2*t2*t26*t27*t6-7.22374750904299d-4*gammaaa
     7         *t19*t22*t36*t38-4.334248505425794d-4*gammaaa*t17*t20*t36
     8         *t38-2.167124252712897d-4*gammaaa*t16*t18*t36*t38-7.22374
     9         750904299d-5*gammaaa*t14*t36*t38+7.22374750904299d-5*gamm
     :         aaa*t13*t34*t38+3.611873754521495d-4*gammaaa*t11*t12*t34*
     ;         t38-7.22374750904299d-5*gammaaa*t2*t32*t38+3.611873754521
     <         495d-5*gammaaa*t15*t23*t38+1.8059368772607476d-4*gammaaa*
     =         t15*t21*t22*t38+1.4447495018085982d-4*gammaaa*t15*t19*t20
     >         *t38+1.0835621263564485d-4*gammaaa*t15*t17*t18*t38+7.2237
     ?         4750904299d-5*gammaaa*t14*t15*t16*t38-9.814103802890335d-
     @         3*t13*t26*t27-4.9070519014451675d-2*t11*t12*t26*t27-9.814
     1         103802890335d-3*t15*t23*t27-4.9070519014451675d-2*t15*t21
     2         *t22*t27-3.925641521156134d-2*t15*t19*t20*t27-2.944231140
     3         8671007d-2*t15*t17*t18*t27-1.962820760578067d-2*t14*t15*t
     4         16*t27
            t40 = 1.4721155704335503d-2*t10*t26*t3*t9+1.1040866778251626
     1         d-2*t26*t3*t7*t8+7.360577852167751d-3*t2*t26*t3*t6+3.6802
     2         889260838756d-3*t13*t26*t3+1.8401444630419378d-2*t11*t12*
     3         t26*t3+3.6802889260838756d-3*t15*t23*t3+1.840144463041937
     4         8d-2*t15*t21*t22*t3+1.4721155704335503d-2*t15*t19*t20*t3+
     5         1.1040866778251626d-2*t15*t17*t18*t3+7.360577852167751d-3
     6         *t14*t15*t16*t3
            t41 = 1/rhoa**5.333333333333333d+0
            t42 = 2.708905315891122d-4*t12*t26*t3*t9+8.126715947673364d-
     1         5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26*t3*t7+2.70890
     2         53158911214d-5*t2*t26*t3
            t43 = gammaaa**3
            t44 = 1/t5**6
            t45 = 1/rhoa**11
            t46 = 1/t5**5
            t47 = 1/rhoa**8.333333333333334d+0
            t48 = 1/rhoa**5.666666666666667d+0
            t49 = exp(-1.1040866778251626d-2*gammaaa*t3)
            t50 = 1/rhoa**10
            t51 = 1/rhoa**9
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t24*wght
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*t28-1.24070098179
     1         88002d+0*t24*t25)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-9.305257363491002d-1*t1*t37-4.1356699
     1         39329334d-1*t24*t30-2.4814019635976003d+0*t25*t28)*wght+A
     2         mat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817988002d+0*t25*t40-9.3052
     1         57363491002d-1*t1*t39)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t26*t3*t42-7.360577852167751d-3*t29*t34*t41+t3*(
     2         2.708905315891122d-4*t19*t22*t3*t36+1.625343189534673d-4*
     3         t17*t20*t3*t36+8.126715947673364d-5*t16*t18*t3*t36+2.7089
     4         053158911214d-5*t14*t3*t36-1.3544526579455607d-5*t15*t23*
     5         t3-6.772263289727805d-5*t15*t21*t22*t3-5.417810631782243d
     6         -5*t15*t19*t20*t3-4.063357973836682d-5*t15*t17*t18*t3-2.7
     7         089053158911214d-5*t14*t15*t16*t3))*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491002d-1*t1*(-6.717208
     1         825089385d-1*gammaaa*t10*t26*t48*t9+8.475863743943775d-3*
     2         t10*t31*t34*t47*t9-2.1189659359859436d-2*t12*t31*t32*t47*
     3         t9+1.1343137263938945d-4*t12*t43*t45*t46*t9-2.26862745278
     4         77887d-5*t10*t32*t43*t45*t9-5.037906618817038d-1*gammaaa*
     5         t26*t48*t7*t8+6.356897807957831d-3*t31*t34*t47*t7*t8-1.70
     6         14705895908416d-5*t32*t43*t45*t7*t8-6.356897807957831d-3*
     7         t31*t32*t47*t6*t8+3.4029411791816827d-5*t43*t45*t46*t6*t8
     8         -5.671568631969471d-6*t43*t44*t45*t8-1.2713795615915663d-
     9         2*t10*t31*t32*t47*t7+6.805882358363365d-5*t10*t43*t45*t46
     :         *t7-5.671568631969472d-5*t12*t43*t44*t45*t7-3.35860441254
     ;         46923d-1*gammaaa*t2*t26*t48*t6+4.237931871971887d-3*t2*t3
     <         1*t34*t47*t6-2.2686274527877887d-5*t10*t43*t44*t45*t6-1.1
     =         343137263938943d-5*t2*t32*t43*t45*t6-5.671568631969472d-5
     >         *t17*t22*t43*t45*t49-2.2686274527877887d-5*t16*t20*t43*t4
     ?         5*t49-5.671568631969471d-6*t18*t43*t45*t49-1.679302206272
     @         3462d-1*gammaaa*t13*t26*t48-8.396511031361731d-1*gammaaa*
     1         t11*t12*t26*t48-1.6793022062723462d-1*gammaaa*t15*t23*t48
     2         -8.396511031361731d-1*gammaaa*t15*t21*t22*t48-6.717208825
     3         089385d-1*gammaaa*t15*t19*t20*t48-5.037906618817038d-1*ga
     4         mmaaa*t15*t17*t18*t48-3.3586044125446923d-1*gammaaa*t14*t
     5         15*t16*t48-2.1189659359859436d-2*t19*t22*t31*t36*t47-1.27
     6         13795615915663d-2*t17*t20*t31*t36*t47-6.356897807957831d-
     7         3*t16*t18*t31*t36*t47-2.1189659359859436d-3*t14*t31*t36*t
     8         47+2.1189659359859436d-3*t13*t31*t34*t47+1.05948296799297
     9         18d-2*t11*t12*t31*t34*t47-2.1189659359859436d-3*t2*t31*t3
     :         2*t47+1.0594829679929718d-3*t15*t23*t31*t47+5.29741483996
     ;         4859d-3*t15*t21*t22*t31*t47+4.237931871971887d-3*t15*t19*
     <         t20*t31*t47+3.1784489039789154d-3*t15*t17*t18*t31*t47+2.1
     =         189659359859436d-3*t14*t15*t16*t31*t47+1.1343137263938943
     >         d-5*t2*t43*t45*t46+5.671568631969472d-5*t19*t22*t36*t43*t
     ?         45+3.4029411791816827d-5*t17*t20*t36*t43*t45+1.7014705895
     @         908414d-5*t16*t18*t36*t43*t45+5.671568631969471d-6*t14*t3
     1         6*t43*t45-5.671568631969471d-6*t13*t32*t43*t45-2.83578431
     2         5984736d-5*t11*t12*t32*t43*t45-9.452614386615786d-7*t15*t
     3         23*t43*t45-4.726307193307893d-6*t15*t21*t22*t43*t45-3.781
     4         0457546463144d-6*t15*t19*t20*t43*t45-2.8357843159847357d-
     5         6*t15*t17*t18*t43*t45-1.8905228773231572d-6*t14*t15*t16*t
     6         43*t45)-3.7221029453964005d+0*t25*t37-1.2407009817988002d
     7         +0*t28*t30+2.7571132928862224d-1*t24/rhoa**1.666666666666
     8         6669d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491002d-1*t1*(-4.25367
     1         64739771044d-5*t12*t31*t46*t50*t9+8.507352947954206d-6*t1
     2         0*t31*t32*t50*t9+1.4394018910905823d-1*t10*t26*t35*t9-2.6
     3         00549103255476d-3*gammaaa*t10*t33*t34*t9+6.50137275813869
     4         d-3*gammaaa*t12*t32*t33*t9+6.380514710965656d-6*t31*t32*t
     5         50*t7*t8+1.0795514183179368d-1*t26*t35*t7*t8-1.9504118274
     6         416074d-3*gammaaa*t33*t34*t7*t8-1.2761029421931314d-5*t31
     7         *t46*t50*t6*t8+1.9504118274416074d-3*gammaaa*t32*t33*t6*t
     8         8+2.1268382369885522d-6*t31*t44*t50*t8-2.552205884386262d
     9         -5*t10*t31*t46*t50*t7+2.1268382369885516d-5*t12*t31*t44*t
     :         50*t7+3.900823654883215d-3*gammaaa*t10*t32*t33*t7+8.50735
     ;         2947954209d-6*t10*t31*t44*t50*t6+4.253676473977103d-6*t2*
     <         t31*t32*t50*t6+7.197009455452912d-2*t2*t26*t35*t6-1.30027
     =         4551627738d-3*gammaaa*t2*t33*t34*t6+2.1268382369885516d-5
     >         *t17*t22*t31*t49*t50+8.507352947954209d-6*t16*t20*t31*t49
     ?         *t50+2.1268382369885522d-6*t18*t31*t49*t50-4.253676473977
     @         103d-6*t2*t31*t46*t50-2.1268382369885522d-5*t19*t22*t31*t
     1         36*t50-1.276102942193131d-5*t17*t20*t31*t36*t50-6.3805147
     2         10965657d-6*t16*t18*t31*t36*t50-2.1268382369885516d-6*t14
     3         *t31*t36*t50+2.1268382369885516d-6*t13*t31*t32*t50+1.0634
     4         191184942758d-5*t11*t12*t31*t32*t50+3.54473039498092d-7*t
     5         15*t23*t31*t50+1.77236519749046d-6*t15*t21*t22*t31*t50+1.
     6         4178921579923678d-6*t15*t19*t20*t31*t50+1.063419118494276
     7         1d-6*t15*t17*t18*t31*t50+7.08946078996184d-7*t14*t15*t16*
     8         t31*t50+6.50137275813869d-3*gammaaa*t19*t22*t33*t36+3.900
     9         823654883215d-3*gammaaa*t17*t20*t33*t36+1.950411827441607
     :         4d-3*gammaaa*t16*t18*t33*t36+6.501372758138692d-4*gammaaa
     ;         *t14*t33*t36+3.598504727726456d-2*t13*t26*t35+1.799252363
     <         863228d-1*t11*t12*t26*t35+3.598504727726456d-2*t15*t23*t3
     =         5+1.799252363863228d-1*t15*t21*t22*t35+1.4394018910905823
     >         d-1*t15*t19*t20*t35+1.0795514183179368d-1*t15*t17*t18*t35
     ?         +7.197009455452912d-2*t14*t15*t16*t35-6.501372758138692d-
     @         4*gammaaa*t13*t33*t34-3.250686379069345d-3*gammaaa*t11*t1
     1         2*t33*t34+6.501372758138692d-4*gammaaa*t2*t32*t33-3.25068
     2         6379069346d-4*gammaaa*t15*t23*t33-1.6253431895346726d-3*g
     3         ammaaa*t15*t21*t22*t33-1.300274551627738d-3*gammaaa*t15*t
     4         19*t20*t33-9.752059137208037d-4*gammaaa*t15*t17*t18*t33-6
     5         .501372758138692d-4*gammaaa*t14*t15*t16*t33)-4.1356699393
     6         29334d-1*t30*t40-2.4814019635976003d+0*t25*t39)*wght+Cmat
     7         3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363491002d-1*t1*(1.59512
     1         8677741414d-5*gammaaa*t12*t46*t51*t9-3.190257355482828d-6
     2         *gammaaa*t10*t32*t51*t9+5.778998007234393d-4*t10*t34*t38*
     3         t9-1.4447495018085982d-3*t12*t32*t38*t9-2.392693016612120
     4         7d-6*gammaaa*t32*t51*t7*t8+4.334248505425794d-4*t34*t38*t
     5         7*t8+4.785386033224242d-6*gammaaa*t46*t51*t6*t8-4.3342485
     6         05425794d-4*t32*t38*t6*t8-7.97564338870707d-7*gammaaa*t44
     7         *t51*t8+9.570772066448484d-6*gammaaa*t10*t46*t51*t7-7.975
     8         64338870707d-6*gammaaa*t12*t44*t51*t7-8.668497010851588d-
     9         4*t10*t32*t38*t7-3.190257355482828d-6*gammaaa*t10*t44*t51
     :         *t6-1.595128677741414d-6*gammaaa*t2*t32*t51*t6+2.88949900
     ;         36171964d-4*t2*t34*t38*t6-7.97564338870707d-6*gammaaa*t17
     <         *t22*t49*t51-3.190257355482828d-6*gammaaa*t16*t20*t49*t51
     =         -7.97564338870707d-7*gammaaa*t18*t49*t51+1.59512867774141
     >         42d-6*gammaaa*t2*t46*t51+7.97564338870707d-6*gammaaa*t19*
     ?         t22*t36*t51+4.785386033224242d-6*gammaaa*t17*t20*t36*t51+
     @         2.392693016612121d-6*gammaaa*t16*t18*t36*t51+7.9756433887
     1         07071d-7*gammaaa*t14*t36*t51-7.97564338870707d-7*gammaaa*
     2         t13*t32*t51-3.987821694353535d-6*gammaaa*t11*t12*t32*t51-
     3         1.329273898117845d-7*gammaaa*t15*t23*t51-6.64636949058922
     4         6d-7*gammaaa*t15*t21*t22*t51-5.31709559247138d-7*gammaaa*
     5         t15*t19*t20*t51-3.987821694353535d-7*gammaaa*t15*t17*t18*
     6         t51-2.65854779623569d-7*gammaaa*t14*t15*t16*t51-1.4447495
     7         018085982d-3*t19*t22*t36*t38-8.668497010851588d-4*t17*t20
     8         *t36*t38-4.334248505425794d-4*t16*t18*t36*t38-1.444749501
     9         8085982d-4*t14*t36*t38+1.4447495018085982d-4*t13*t34*t38+
     :         7.22374750904299d-4*t11*t12*t34*t38-1.4447495018085982d-4
     ;         *t2*t32*t38+7.22374750904299d-5*t15*t23*t38+3.61187375452
     <         1495d-4*t15*t21*t22*t38+2.8894990036171964d-4*t15*t19*t20
     =         *t38+2.167124252712897d-4*t15*t17*t18*t38+1.4447495018085
     >         982d-4*t14*t15*t16*t38)-1.2407009817988002d+0*t25*(-1.083
     ?         5621263564485d-4*t10*t34*t41*t9+2.708905315891122d-4*t12*
     @         t32*t41*t9-8.126715947673364d-5*t34*t41*t7*t8+8.126715947
     1         673364d-5*t32*t41*t6*t8+1.625343189534673d-4*t10*t32*t41*
     2         t7-5.417810631782243d-5*t2*t34*t41*t6+2.708905315891122d-
     3         4*t19*t22*t36*t41+1.625343189534673d-4*t17*t20*t36*t41+8.
     4         126715947673364d-5*t16*t18*t36*t41+2.7089053158911214d-5*
     5         t14*t36*t41-2.7089053158911214d-5*t13*t34*t41-1.354452657
     6         945561d-4*t11*t12*t34*t41+2.7089053158911214d-5*t2*t32*t4
     7         1-1.3544526579455607d-5*t15*t23*t41-6.772263289727805d-5*
     8         t15*t21*t22*t41-5.417810631782243d-5*t15*t19*t20*t41-4.06
     9         3357973836682d-5*t15*t17*t18*t41-2.7089053158911214d-5*t1
     :         4*t15*t16*t41))*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491002d-1*t1*(t26*t3*(-1.9939108471767675d-6*t12*t34*t
     2         41*t9-5.981732541530301d-7*t34*t41*t6*t8+2.99086627076515
     3         07d-7*t32*t41*t8-1.1963465083060604d-6*t10*t34*t41*t7+2.9
     4         90866270765151d-6*t12*t32*t41*t7+1.1963465083060604d-6*t1
     5         0*t32*t41*t6-1.9939108471767675d-7*t2*t34*t41)+t3*(2.9908
     6         66270765151d-6*t17*t22*t41*t49+1.1963465083060604d-6*t16*
     7         t20*t41*t49+2.9908662707651507d-7*t18*t41*t49-2.990866270
     8         765151d-6*t19*t22*t36*t41-1.7945197624590903d-6*t17*t20*t
     9         36*t41-8.972598812295453d-7*t16*t18*t36*t41-2.99086627076
     :         5151d-7*t14*t36*t41+4.9847771179419187d-8*t15*t23*t41+2.4
     ;         923885589709593d-7*t15*t21*t22*t41+1.9939108471767675d-7*
     <         t15*t19*t20*t41+1.4954331353825753d-7*t15*t17*t18*t41+9.9
     =         69554235883837d-8*t14*t15*t16*t41)-1.4721155704335503d-2*
     >         t34*t41*t42+8.126715947673364d-5*t29*t32/rhoa**8)*wght
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
            t2 = param(1)
            t3 = param(7)
            t4 = param(3)
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 3.6802889260838756d-3*gammaaa*t5
            t7 = t6+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = t8**2.0d+0
            t10 = param(4)
            t11 = t8**3.0d+0
            t12 = param(5)
            t13 = t8**4.0d+0
            t14 = param(6)
            t15 = param(2)
            t16 = param(9)
            t17 = exp(-t6)
            t18 = 1.0d+0-t17
            t19 = t18**2.0d+0
            t20 = param(10)
            t21 = t18**3.0d+0
            t22 = param(11)
            t23 = t18**4.0d+0
            t24 = param(12)
            t25 = param(8)
            t26 = t4*t9+t14*t8**5.0d+0+t15*t8+t3+t18*t25+t18**5.0d+0*t24
     1         +t22*t23+t20*t21+t2+t16*t19+t12*t13+t10*t11
            t27 = rhob**1.3333333333333333d+0
            t28 = 1/rhob**2.6666666666666666d+0
            t29 = 3.6802889260838756d-3*gammabb*t28
            t30 = t29+1.0d+0
            t31 = 1.0d+0-1.0d+0/t30
            t32 = t31**2.0d+0
            t33 = t31**3.0d+0
            t34 = t31**4.0d+0
            t35 = exp(-t29)
            t36 = 1.0d+0-t35
            t37 = t36**2.0d+0
            t38 = t36**3.0d+0
            t39 = t36**4.0d+0
            t40 = t32*t4+t22*t39+t20*t38+t16*t37+t24*t36**5.0d+0+t25*t36
     1         +t12*t34+t10*t33+t14*t31**5.0d+0+t15*t31+t3+t2
            t41 = rhoa**3.333333333333333d-1
            t42 = 1/t7**2
            t43 = 1/rhoa**3.6666666666666664d+0
            t44 = -2.9442311408671007d-2*gammaaa*t10*t42*t43*t9-1.962820
     1         760578067d-2*gammaaa*t4*t42*t43*t8-9.814103802890335d-3*g
     2         ammaaa*t15*t42*t43-4.9070519014451675d-2*gammaaa*t13*t14*
     3         t42*t43-3.925641521156134d-2*gammaaa*t11*t12*t42*t43-9.81
     4         4103802890335d-3*gammaaa*t17*t25*t43-4.9070519014451675d-
     5         2*gammaaa*t17*t23*t24*t43-3.925641521156134d-2*gammaaa*t1
     6         7*t21*t22*t43-2.9442311408671007d-2*gammaaa*t17*t19*t20*t
     7         43-1.962820760578067d-2*gammaaa*t16*t17*t18*t43
            t45 = rhob**3.333333333333333d-1
            t46 = 1/t30**2
            t47 = 1/rhob**3.6666666666666664d+0
            t48 = -1.962820760578067d-2*gammabb*t31*t4*t46*t47-4.9070519
     1         014451675d-2*gammabb*t14*t34*t46*t47-3.925641521156134d-2
     2         *gammabb*t12*t33*t46*t47-2.9442311408671007d-2*gammabb*t1
     3         0*t32*t46*t47-9.814103802890335d-3*gammabb*t15*t46*t47-4.
     4         9070519014451675d-2*gammabb*t24*t35*t39*t47-3.92564152115
     5         6134d-2*gammabb*t22*t35*t38*t47-2.9442311408671007d-2*gam
     6         mabb*t20*t35*t37*t47-1.962820760578067d-2*gammabb*t16*t35
     7         *t36*t47-9.814103802890335d-3*gammabb*t25*t35*t47
            t49 = 3.6802889260838756d-3*t15
            t50 = 1.1040866778251626d-2*t10*t9+7.360577852167751d-3*t4*t
     1         8+t49+1.8401444630419378d-2*t13*t14+1.4721155704335503d-2
     2         *t11*t12
            t51 = t49+7.360577852167751d-3*t31*t4+1.8401444630419378d-2*
     1         t14*t34+1.4721155704335503d-2*t12*t33+1.1040866778251626d
     2         -2*t10*t32
            t52 = 1/rhoa**6.666666666666666d-1
            t53 = gammaaa**2
            t54 = 1/t7**4
            t55 = 1/rhoa**7.333333333333333d+0
            t56 = 1/t7**3
            t57 = 1/rhoa**4.666666666666667d+0
            t58 = exp(-7.360577852167751d-3*gammaaa*t5)
            t59 = 1.0795514183179368d-1*gammaaa*t10*t42*t57*t9-5.7789980
     1         07234393d-4*t10*t53*t55*t56*t9+1.1557996014468785d-3*t12*
     2         t53*t54*t55*t9+7.197009455452912d-2*gammaaa*t4*t42*t57*t8
     3         -3.8526653381562614d-4*t4*t53*t55*t56*t8+5.77899800723439
     4         3d-4*t10*t53*t54*t55*t8+1.9263326690781307d-3*t21*t24*t53
     5         *t55*t58+1.1557996014468785d-3*t19*t22*t53*t55*t58+5.7789
     6         98007234393d-4*t18*t20*t53*t55*t58+1.9263326690781307d-4*
     7         t16*t53*t55*t58+3.598504727726456d-2*gammaaa*t15*t42*t57+
     8         1.799252363863228d-1*gammaaa*t13*t14*t42*t57+1.4394018910
     9         905823d-1*gammaaa*t11*t12*t42*t57+3.598504727726456d-2*ga
     :         mmaaa*t17*t25*t57+1.799252363863228d-1*gammaaa*t17*t23*t2
     ;         4*t57+1.4394018910905823d-1*gammaaa*t17*t21*t22*t57+1.079
     <         5514183179368d-1*gammaaa*t17*t19*t20*t57+7.19700945545291
     =         2d-2*gammaaa*t16*t17*t18*t57-1.9263326690781307d-4*t15*t5
     >         3*t55*t56-9.631663345390654d-4*t13*t14*t53*t55*t56-7.7053
     ?         30676312523d-4*t11*t12*t53*t55*t56+1.9263326690781307d-4*
     @         t4*t53*t54*t55+1.9263326690781307d-3*t11*t14*t53*t54*t55-
     1         9.631663345390653d-5*t17*t25*t53*t55-4.815831672695327d-4
     2         *t17*t23*t24*t53*t55-3.8526653381562614d-4*t17*t21*t22*t5
     3         3*t55-2.8894990036171964d-4*t17*t19*t20*t53*t55-1.9263326
     4         690781307d-4*t16*t17*t18*t53*t55
            t60 = 1/rhob**6.666666666666666d-1
            t61 = gammabb**2
            t62 = 1/t30**4
            t63 = 1/rhob**7.333333333333333d+0
            t64 = 1/t30**3
            t65 = 1/rhob**4.666666666666667d+0
            t66 = exp(-7.360577852167751d-3*gammabb*t28)
            t67 = 1.9263326690781307d-3*t24*t38*t61*t63*t66+1.1557996014
     1         468785d-3*t22*t37*t61*t63*t66+5.778998007234393d-4*t20*t3
     2         6*t61*t63*t66+1.9263326690781307d-4*t16*t61*t63*t66+7.197
     3         009455452912d-2*gammabb*t31*t4*t46*t65+1.799252363863228d
     4         -1*gammabb*t14*t34*t46*t65+1.4394018910905823d-1*gammabb*
     5         t12*t33*t46*t65+1.0795514183179368d-1*gammabb*t10*t32*t46
     6         *t65+3.598504727726456d-2*gammabb*t15*t46*t65+1.799252363
     7         863228d-1*gammabb*t24*t35*t39*t65+1.4394018910905823d-1*g
     8         ammabb*t22*t35*t38*t65+1.0795514183179368d-1*gammabb*t20*
     9         t35*t37*t65+7.197009455452912d-2*gammabb*t16*t35*t36*t65+
     :         3.598504727726456d-2*gammabb*t25*t35*t65-3.85266533815626
     ;         14d-4*t31*t4*t61*t63*t64-9.631663345390654d-4*t14*t34*t61
     <         *t63*t64-7.705330676312523d-4*t12*t33*t61*t63*t64-5.77899
     =         8007234393d-4*t10*t32*t61*t63*t64-1.9263326690781307d-4*t
     >         15*t61*t63*t64+1.9263326690781307d-4*t4*t61*t62*t63+1.926
     ?         3326690781307d-3*t14*t33*t61*t62*t63+1.1557996014468785d-
     @         3*t12*t32*t61*t62*t63+5.778998007234393d-4*t10*t31*t61*t6
     1         2*t63-4.815831672695327d-4*t24*t35*t39*t61*t63-3.85266533
     2         81562614d-4*t22*t35*t38*t61*t63-2.8894990036171964d-4*t20
     3         *t35*t37*t61*t63-1.9263326690781307d-4*t16*t35*t36*t61*t6
     4         3-9.631663345390653d-5*t25*t35*t61*t63
            t68 = 1/rhoa**6.333333333333333d+0
            t69 = 2.167124252712897d-4*gammaaa*t10*t56*t68*t9-4.33424850
     1         5425794d-4*gammaaa*t12*t54*t68*t9-2.9442311408671007d-2*t
     2         10*t42*t43*t9+1.4447495018085982d-4*gammaaa*t4*t56*t68*t8
     3         -2.167124252712897d-4*gammaaa*t10*t54*t68*t8-1.9628207605
     4         78067d-2*t4*t42*t43*t8-7.22374750904299d-4*gammaaa*t21*t2
     5         4*t58*t68-4.334248505425794d-4*gammaaa*t19*t22*t58*t68-2.
     6         167124252712897d-4*gammaaa*t18*t20*t58*t68-7.223747509042
     7         99d-5*gammaaa*t16*t58*t68+7.22374750904299d-5*gammaaa*t15
     8         *t56*t68+3.611873754521495d-4*gammaaa*t13*t14*t56*t68+2.8
     9         894990036171964d-4*gammaaa*t11*t12*t56*t68-7.223747509042
     :         99d-5*gammaaa*t4*t54*t68-7.22374750904299d-4*gammaaa*t11*
     ;         t14*t54*t68+3.611873754521495d-5*gammaaa*t17*t25*t68+1.80
     <         59368772607476d-4*gammaaa*t17*t23*t24*t68+1.4447495018085
     =         982d-4*gammaaa*t17*t21*t22*t68+1.0835621263564485d-4*gamm
     >         aaa*t17*t19*t20*t68+7.22374750904299d-5*gammaaa*t16*t17*t
     ?         18*t68-9.814103802890335d-3*t15*t42*t43-4.907051901445167
     @         5d-2*t13*t14*t42*t43-3.925641521156134d-2*t11*t12*t42*t43
     1         -9.814103802890335d-3*t17*t25*t43-4.9070519014451675d-2*t
     2         17*t23*t24*t43-3.925641521156134d-2*t17*t21*t22*t43-2.944
     3         2311408671007d-2*t17*t19*t20*t43-1.962820760578067d-2*t16
     4         *t17*t18*t43
            t70 = 1.1040866778251626d-2*t10*t42*t5*t9+7.360577852167751d
     1         -3*t4*t42*t5*t8+3.6802889260838756d-3*t15*t42*t5+1.840144
     2         4630419378d-2*t13*t14*t42*t5+1.4721155704335503d-2*t11*t1
     3         2*t42*t5+3.6802889260838756d-3*t17*t25*t5+1.8401444630419
     4         378d-2*t17*t23*t24*t5+1.4721155704335503d-2*t17*t21*t22*t
     5         5+1.1040866778251626d-2*t17*t19*t20*t5+7.360577852167751d
     6         -3*t16*t17*t18*t5
            t71 = 1/rhob**6.333333333333333d+0
            t72 = -7.22374750904299d-4*gammabb*t24*t38*t66*t71-4.3342485
     1         05425794d-4*gammabb*t22*t37*t66*t71-2.167124252712897d-4*
     2         gammabb*t20*t36*t66*t71-7.22374750904299d-5*gammabb*t16*t
     3         66*t71+1.4447495018085982d-4*gammabb*t31*t4*t64*t71+3.611
     4         873754521495d-4*gammabb*t14*t34*t64*t71+2.889499003617196
     5         4d-4*gammabb*t12*t33*t64*t71+2.167124252712897d-4*gammabb
     6         *t10*t32*t64*t71+7.22374750904299d-5*gammabb*t15*t64*t71-
     7         7.22374750904299d-5*gammabb*t4*t62*t71-7.22374750904299d-
     8         4*gammabb*t14*t33*t62*t71-4.334248505425794d-4*gammabb*t1
     9         2*t32*t62*t71-2.167124252712897d-4*gammabb*t10*t31*t62*t7
     :         1+1.8059368772607476d-4*gammabb*t24*t35*t39*t71+1.4447495
     ;         018085982d-4*gammabb*t22*t35*t38*t71+1.0835621263564485d-
     <         4*gammabb*t20*t35*t37*t71+7.22374750904299d-5*gammabb*t16
     =         *t35*t36*t71+3.611873754521495d-5*gammabb*t25*t35*t71-1.9
     >         62820760578067d-2*t31*t4*t46*t47-4.9070519014451675d-2*t1
     ?         4*t34*t46*t47-3.925641521156134d-2*t12*t33*t46*t47-2.9442
     @         311408671007d-2*t10*t32*t46*t47-9.814103802890335d-3*t15*
     1         t46*t47-4.9070519014451675d-2*t24*t35*t39*t47-3.925641521
     2         156134d-2*t22*t35*t38*t47-2.9442311408671007d-2*t20*t35*t
     3         37*t47-1.962820760578067d-2*t16*t35*t36*t47-9.81410380289
     4         0335d-3*t25*t35*t47
            t73 = 7.360577852167751d-3*t28*t31*t4*t46+1.8401444630419378
     1         d-2*t14*t28*t34*t46+1.4721155704335503d-2*t12*t28*t33*t46
     2         +1.1040866778251626d-2*t10*t28*t32*t46+3.6802889260838756
     3         d-3*t15*t28*t46+1.8401444630419378d-2*t24*t28*t35*t39+1.4
     4         721155704335503d-2*t22*t28*t35*t38+1.1040866778251626d-2*
     5         t20*t28*t35*t37+7.360577852167751d-3*t16*t28*t35*t36+3.68
     6         02889260838756d-3*t25*t28*t35
            t74 = 1/rhoa**5.333333333333333d+0
            t75 = 1.625343189534673d-4*t12*t42*t5*t9+8.126715947673364d-
     1         5*t10*t42*t5*t8+2.7089053158911214d-5*t4*t42*t5+2.7089053
     2         15891122d-4*t11*t14*t42*t5
            t76 = 1/rhob**5.333333333333333d+0
            t77 = 2.7089053158911214d-5*t28*t4*t46+2.708905315891122d-4*
     1         t14*t28*t33*t46+1.625343189534673d-4*t12*t28*t32*t46+8.12
     2         6715947673364d-5*t10*t28*t31*t46
            t78 = gammaaa**3
            t79 = 1/t7**6
            t80 = 1/rhoa**11
            t81 = 1/t7**5
            t82 = 1/rhoa**8.333333333333334d+0
            t83 = 1/rhoa**5.666666666666667d+0
            t84 = exp(-1.1040866778251626d-2*gammaaa*t5)
            t85 = gammabb**3
            t86 = 1/t30**6
            t87 = 1/rhob**11
            t88 = 1/t30**5
            t89 = 1/rhob**8.333333333333334d+0
            t90 = 1/rhob**5.666666666666667d+0
            t91 = exp(-1.1040866778251626d-2*gammabb*t28)
            t92 = 1/rhoa**10
            t93 = 1/rhob**10
            t94 = 1/rhoa**9
            t95 = 1/rhob**9
            fnc(iq) = (-9.305257363491002d-1*t27*t40-9.305257363491002d-
     1         1*t1*t26)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491002d-1*t1*t44-1.24070098179
     1         88002d+0*t26*t41)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491002d-1*t27*t48-1.2407009817
     1         988002d+0*t40*t45)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         42*t5*t50+(3.6802889260838756d-3*t17*t25+1.84014446304193
     2         78d-2*t17*t23*t24+1.4721155704335503d-2*t17*t21*t22+1.104
     3         0866778251626d-2*t17*t19*t20+7.360577852167751d-3*t16*t17
     4         *t18)*t5)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t27*(
     1         t28*t46*t51+t28*(1.8401444630419378d-2*t24*t35*t39+1.4721
     2         155704335503d-2*t22*t35*t38+1.1040866778251626d-2*t20*t35
     3         *t37+7.360577852167751d-3*t16*t35*t36+3.6802889260838756d
     4         -3*t25*t35))*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491002d-1*t1*t59-4.1356699
     1         39329334d-1*t26*t52-2.4814019635976003d+0*t41*t44)*wght+A
     2         mat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491002d-1*t27*t67-4.135669
     1         939329334d-1*t40*t60-2.4814019635976003d+0*t45*t48)*wght+
     2         Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817988002d+0*t41*t70-9.3052
     1         57363491002d-1*t1*t69)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.2407009817988002d+0*t45*t73-9.3052
     1         57363491002d-1*t27*t72)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t42*t5*t75-7.360577852167751d-3*t50*t56*t74+t5*(
     2         2.708905315891122d-4*t21*t24*t5*t58+1.625343189534673d-4*
     3         t19*t22*t5*t58+8.126715947673364d-5*t18*t20*t5*t58+2.7089
     4         053158911214d-5*t16*t5*t58-1.3544526579455607d-5*t17*t25*
     5         t5-6.772263289727805d-5*t17*t23*t24*t5-5.417810631782243d
     6         -5*t17*t21*t22*t5-4.063357973836682d-5*t17*t19*t20*t5-2.7
     7         089053158911214d-5*t16*t17*t18*t5))*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349100
     1         2d-1*t27*(t28*t46*t77-7.360577852167751d-3*t51*t64*t76+t2
     2         8*(2.708905315891122d-4*t24*t28*t38*t66+1.625343189534673
     3         d-4*t22*t28*t37*t66+8.126715947673364d-5*t20*t28*t36*t66+
     4         2.7089053158911214d-5*t16*t28*t66-6.772263289727805d-5*t2
     5         4*t28*t35*t39-5.417810631782243d-5*t22*t28*t35*t38-4.0633
     6         57973836682d-5*t20*t28*t35*t37-2.7089053158911214d-5*t16*
     7         t28*t35*t36-1.3544526579455607d-5*t25*t28*t35))*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491002d-1*t1*(-5.037906
     1         618817038d-1*gammaaa*t10*t42*t83*t9+6.356897807957831d-3*
     2         t10*t53*t56*t82*t9-1.2713795615915663d-2*t12*t53*t54*t82*
     3         t9+6.805882358363365d-5*t12*t78*t80*t81*t9-5.671568631969
     4         472d-5*t14*t78*t79*t80*t9-1.7014705895908416d-5*t10*t54*t
     5         78*t80*t9-5.671568631969472d-5*t19*t24*t78*t80*t84-2.2686
     6         274527877887d-5*t18*t22*t78*t80*t84-5.671568631969471d-6*
     7         t20*t78*t80*t84-3.3586044125446923d-1*gammaaa*t4*t42*t8*t
     8         83-1.6793022062723462d-1*gammaaa*t15*t42*t83-8.3965110313
     9         61731d-1*gammaaa*t13*t14*t42*t83-6.717208825089385d-1*gam
     :         maaa*t11*t12*t42*t83-1.6793022062723462d-1*gammaaa*t17*t2
     ;         5*t83-8.396511031361731d-1*gammaaa*t17*t23*t24*t83-6.7172
     <         08825089385d-1*gammaaa*t17*t21*t22*t83-5.037906618817038d
     =         -1*gammaaa*t17*t19*t20*t83-3.3586044125446923d-1*gammaaa*
     >         t16*t17*t18*t83+4.237931871971887d-3*t4*t53*t56*t8*t82-6.
     ?         356897807957831d-3*t10*t53*t54*t8*t82-2.1189659359859436d
     @         -2*t21*t24*t53*t58*t82-1.2713795615915663d-2*t19*t22*t53*
     1         t58*t82-6.356897807957831d-3*t18*t20*t53*t58*t82-2.118965
     2         9359859436d-3*t16*t53*t58*t82+2.1189659359859436d-3*t15*t
     3         53*t56*t82+1.0594829679929718d-2*t13*t14*t53*t56*t82+8.47
     4         5863743943775d-3*t11*t12*t53*t56*t82-2.1189659359859436d-
     5         3*t4*t53*t54*t82-2.1189659359859436d-2*t11*t14*t53*t54*t8
     6         2+1.0594829679929718d-3*t17*t25*t53*t82+5.297414839964859
     7         d-3*t17*t23*t24*t53*t82+4.237931871971887d-3*t17*t21*t22*
     8         t53*t82+3.1784489039789154d-3*t17*t19*t20*t53*t82+2.11896
     9         59359859436d-3*t16*t17*t18*t53*t82+3.4029411791816827d-5*
     :         t10*t78*t8*t80*t81+1.1343137263938943d-5*t4*t78*t80*t81+1
     ;         .1343137263938945d-4*t11*t14*t78*t80*t81-2.26862745278778
     <         87d-5*t12*t78*t79*t8*t80-1.1343137263938943d-5*t4*t54*t78
     =         *t8*t80-5.671568631969471d-6*t10*t78*t79*t80+5.6715686319
     >         69472d-5*t21*t24*t58*t78*t80+3.4029411791816827d-5*t19*t2
     ?         2*t58*t78*t80+1.7014705895908414d-5*t18*t20*t58*t78*t80+5
     @         .671568631969471d-6*t16*t58*t78*t80-5.671568631969471d-6*
     1         t15*t54*t78*t80-2.835784315984736d-5*t13*t14*t54*t78*t80-
     2         2.2686274527877887d-5*t11*t12*t54*t78*t80-9.4526143866157
     3         86d-7*t17*t25*t78*t80-4.726307193307893d-6*t17*t23*t24*t7
     4         8*t80-3.7810457546463144d-6*t17*t21*t22*t78*t80-2.8357843
     5         159847357d-6*t17*t19*t20*t78*t80-1.8905228773231572d-6*t1
     6         6*t17*t18*t78*t80)-3.7221029453964005d+0*t41*t59-1.240700
     7         9817988002d+0*t44*t52+2.7571132928862224d-1*t26/rhoa**1.6
     8         666666666666669d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363491002d-1*t27*(-5.67156
     1         8631969472d-5*t24*t37*t85*t87*t91-2.2686274527877887d-5*t
     2         22*t36*t85*t87*t91-5.671568631969471d-6*t20*t85*t87*t91-3
     3         .3586044125446923d-1*gammabb*t31*t4*t46*t90-8.39651103136
     4         1731d-1*gammabb*t14*t34*t46*t90-6.717208825089385d-1*gamm
     5         abb*t12*t33*t46*t90-5.037906618817038d-1*gammabb*t10*t32*
     6         t46*t90-1.6793022062723462d-1*gammabb*t15*t46*t90-8.39651
     7         1031361731d-1*gammabb*t24*t35*t39*t90-6.717208825089385d-
     8         1*gammabb*t22*t35*t38*t90-5.037906618817038d-1*gammabb*t2
     9         0*t35*t37*t90-3.3586044125446923d-1*gammabb*t16*t35*t36*t
     :         90-1.6793022062723462d-1*gammabb*t25*t35*t90-2.1189659359
     ;         859436d-2*t24*t38*t61*t66*t89-1.2713795615915663d-2*t22*t
     <         37*t61*t66*t89-6.356897807957831d-3*t20*t36*t61*t66*t89-2
     =         .1189659359859436d-3*t16*t61*t66*t89+4.237931871971887d-3
     >         *t31*t4*t61*t64*t89+1.0594829679929718d-2*t14*t34*t61*t64
     ?         *t89+8.475863743943775d-3*t12*t33*t61*t64*t89+6.356897807
     @         957831d-3*t10*t32*t61*t64*t89+2.1189659359859436d-3*t15*t
     1         61*t64*t89-2.1189659359859436d-3*t4*t61*t62*t89-2.1189659
     2         359859436d-2*t14*t33*t61*t62*t89-1.2713795615915663d-2*t1
     3         2*t32*t61*t62*t89-6.356897807957831d-3*t10*t31*t61*t62*t8
     4         9+5.297414839964859d-3*t24*t35*t39*t61*t89+4.237931871971
     5         887d-3*t22*t35*t38*t61*t89+3.1784489039789154d-3*t20*t35*
     6         t37*t61*t89+2.1189659359859436d-3*t16*t35*t36*t61*t89+1.0
     7         594829679929718d-3*t25*t35*t61*t89+1.1343137263938943d-5*
     8         t4*t85*t87*t88+1.1343137263938945d-4*t14*t33*t85*t87*t88+
     9         6.805882358363365d-5*t12*t32*t85*t87*t88+3.40294117918168
     :         27d-5*t10*t31*t85*t87*t88-5.671568631969472d-5*t14*t32*t8
     ;         5*t86*t87-2.2686274527877887d-5*t12*t31*t85*t86*t87-5.671
     <         568631969471d-6*t10*t85*t86*t87+5.671568631969472d-5*t24*
     =         t38*t66*t85*t87+3.4029411791816827d-5*t22*t37*t66*t85*t87
     >         +1.7014705895908414d-5*t20*t36*t66*t85*t87+5.671568631969
     ?         471d-6*t16*t66*t85*t87-1.1343137263938943d-5*t31*t4*t62*t
     @         85*t87-2.835784315984736d-5*t14*t34*t62*t85*t87-2.2686274
     1         527877887d-5*t12*t33*t62*t85*t87-1.7014705895908416d-5*t1
     2         0*t32*t62*t85*t87-5.671568631969471d-6*t15*t62*t85*t87-4.
     3         726307193307893d-6*t24*t35*t39*t85*t87-3.7810457546463144
     4         d-6*t22*t35*t38*t85*t87-2.8357843159847357d-6*t20*t35*t37
     5         *t85*t87-1.8905228773231572d-6*t16*t35*t36*t85*t87-9.4526
     6         14386615786d-7*t25*t35*t85*t87)-3.7221029453964005d+0*t45
     7         *t67-1.2407009817988002d+0*t48*t60+2.7571132928862224d-1*
     8         t40/rhob**1.6666666666666669d+0)*wght+Amat3(iq,D3_RB_RB_R
     9         B)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491002d-1*t1*(-2.55220
     1         5884386262d-5*t12*t53*t81*t9*t92+2.1268382369885516d-5*t1
     2         4*t53*t79*t9*t92+6.380514710965656d-6*t10*t53*t54*t9*t92+
     3         2.1268382369885516d-5*t19*t24*t53*t84*t92+8.5073529479542
     4         09d-6*t18*t22*t53*t84*t92+2.1268382369885522d-6*t20*t53*t
     5         84*t92-1.2761029421931314d-5*t10*t53*t8*t81*t92-4.2536764
     6         73977103d-6*t4*t53*t81*t92-4.2536764739771044d-5*t11*t14*
     7         t53*t81*t92+8.507352947954209d-6*t12*t53*t79*t8*t92+4.253
     8         676473977103d-6*t4*t53*t54*t8*t92+2.1268382369885522d-6*t
     9         10*t53*t79*t92-2.1268382369885522d-5*t21*t24*t53*t58*t92-
     :         1.276102942193131d-5*t19*t22*t53*t58*t92-6.38051471096565
     ;         7d-6*t18*t20*t53*t58*t92-2.1268382369885516d-6*t16*t53*t5
     <         8*t92+2.1268382369885516d-6*t15*t53*t54*t92+1.06341911849
     =         42758d-5*t13*t14*t53*t54*t92+8.507352947954206d-6*t11*t12
     >         *t53*t54*t92+3.54473039498092d-7*t17*t25*t53*t92+1.772365
     ?         19749046d-6*t17*t23*t24*t53*t92+1.4178921579923678d-6*t17
     @         *t21*t22*t53*t92+1.0634191184942761d-6*t17*t19*t20*t53*t9
     1         2+7.08946078996184d-7*t16*t17*t18*t53*t92+1.0795514183179
     2         368d-1*t10*t42*t57*t9-1.9504118274416074d-3*gammaaa*t10*t
     3         55*t56*t9+3.900823654883215d-3*gammaaa*t12*t54*t55*t9+7.1
     4         97009455452912d-2*t4*t42*t57*t8-1.300274551627738d-3*gamm
     5         aaa*t4*t55*t56*t8+1.9504118274416074d-3*gammaaa*t10*t54*t
     6         55*t8+6.50137275813869d-3*gammaaa*t21*t24*t55*t58+3.90082
     7         3654883215d-3*gammaaa*t19*t22*t55*t58+1.9504118274416074d
     8         -3*gammaaa*t18*t20*t55*t58+6.501372758138692d-4*gammaaa*t
     9         16*t55*t58+3.598504727726456d-2*t15*t42*t57+1.79925236386
     :         3228d-1*t13*t14*t42*t57+1.4394018910905823d-1*t11*t12*t42
     ;         *t57+3.598504727726456d-2*t17*t25*t57+1.799252363863228d-
     <         1*t17*t23*t24*t57+1.4394018910905823d-1*t17*t21*t22*t57+1
     =         .0795514183179368d-1*t17*t19*t20*t57+7.197009455452912d-2
     >         *t16*t17*t18*t57-6.501372758138692d-4*gammaaa*t15*t55*t56
     ?         -3.250686379069345d-3*gammaaa*t13*t14*t55*t56-2.600549103
     @         255476d-3*gammaaa*t11*t12*t55*t56+6.501372758138692d-4*ga
     1         mmaaa*t4*t54*t55+6.50137275813869d-3*gammaaa*t11*t14*t54*
     2         t55-3.250686379069346d-4*gammaaa*t17*t25*t55-1.6253431895
     3         346726d-3*gammaaa*t17*t23*t24*t55-1.300274551627738d-3*ga
     4         mmaaa*t17*t21*t22*t55-9.752059137208037d-4*gammaaa*t17*t1
     5         9*t20*t55-6.501372758138692d-4*gammaaa*t16*t17*t18*t55)-4
     6         .135669939329334d-1*t52*t70-2.4814019635976003d+0*t41*t69
     7         )*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-9.305257363491002d-1*t27*(2.12683
     1         82369885516d-5*t24*t37*t61*t91*t93+8.507352947954209d-6*t
     2         22*t36*t61*t91*t93+2.1268382369885522d-6*t20*t61*t91*t93-
     3         4.253676473977103d-6*t4*t61*t88*t93-4.2536764739771044d-5
     4         *t14*t33*t61*t88*t93-2.552205884386262d-5*t12*t32*t61*t88
     5         *t93-1.2761029421931314d-5*t10*t31*t61*t88*t93+2.12683823
     6         69885516d-5*t14*t32*t61*t86*t93+8.507352947954209d-6*t12*
     7         t31*t61*t86*t93+2.1268382369885522d-6*t10*t61*t86*t93-2.1
     8         268382369885522d-5*t24*t38*t61*t66*t93-1.276102942193131d
     9         -5*t22*t37*t61*t66*t93-6.380514710965657d-6*t20*t36*t61*t
     :         66*t93-2.1268382369885516d-6*t16*t61*t66*t93+4.2536764739
     ;         77103d-6*t31*t4*t61*t62*t93+1.0634191184942758d-5*t14*t34
     <         *t61*t62*t93+8.507352947954206d-6*t12*t33*t61*t62*t93+6.3
     =         80514710965656d-6*t10*t32*t61*t62*t93+2.1268382369885516d
     >         -6*t15*t61*t62*t93+1.77236519749046d-6*t24*t35*t39*t61*t9
     ?         3+1.4178921579923678d-6*t22*t35*t38*t61*t93+1.06341911849
     @         42761d-6*t20*t35*t37*t61*t93+7.08946078996184d-7*t16*t35*
     1         t36*t61*t93+3.54473039498092d-7*t25*t35*t61*t93+6.5013727
     2         5813869d-3*gammabb*t24*t38*t63*t66+3.900823654883215d-3*g
     3         ammabb*t22*t37*t63*t66+1.9504118274416074d-3*gammabb*t20*
     4         t36*t63*t66+6.501372758138692d-4*gammabb*t16*t63*t66+7.19
     5         7009455452912d-2*t31*t4*t46*t65+1.799252363863228d-1*t14*
     6         t34*t46*t65+1.4394018910905823d-1*t12*t33*t46*t65+1.07955
     7         14183179368d-1*t10*t32*t46*t65+3.598504727726456d-2*t15*t
     8         46*t65+1.799252363863228d-1*t24*t35*t39*t65+1.43940189109
     9         05823d-1*t22*t35*t38*t65+1.0795514183179368d-1*t20*t35*t3
     :         7*t65+7.197009455452912d-2*t16*t35*t36*t65+3.598504727726
     ;         456d-2*t25*t35*t65-1.300274551627738d-3*gammabb*t31*t4*t6
     <         3*t64-3.250686379069345d-3*gammabb*t14*t34*t63*t64-2.6005
     =         49103255476d-3*gammabb*t12*t33*t63*t64-1.9504118274416074
     >         d-3*gammabb*t10*t32*t63*t64-6.501372758138692d-4*gammabb*
     ?         t15*t63*t64+6.501372758138692d-4*gammabb*t4*t62*t63+6.501
     @         37275813869d-3*gammabb*t14*t33*t62*t63+3.900823654883215d
     1         -3*gammabb*t12*t32*t62*t63+1.9504118274416074d-3*gammabb*
     2         t10*t31*t62*t63-1.6253431895346726d-3*gammabb*t24*t35*t39
     3         *t63-1.300274551627738d-3*gammabb*t22*t35*t38*t63-9.75205
     4         9137208037d-4*gammabb*t20*t35*t37*t63-6.501372758138692d-
     5         4*gammabb*t16*t35*t36*t63-3.250686379069346d-4*gammabb*t2
     6         5*t35*t63)-4.135669939329334d-1*t60*t73-2.481401963597600
     7         3d+0*t45*t72)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363491002d-1*t1*(9.57077
     1         2066448484d-6*gammaaa*t12*t81*t9*t94-7.97564338870707d-6*
     2         gammaaa*t14*t79*t9*t94-2.3926930166121207d-6*gammaaa*t10*
     3         t54*t9*t94-7.97564338870707d-6*gammaaa*t19*t24*t84*t94-3.
     4         190257355482828d-6*gammaaa*t18*t22*t84*t94-7.975643388707
     5         07d-7*gammaaa*t20*t84*t94+4.785386033224242d-6*gammaaa*t1
     6         0*t8*t81*t94+1.5951286777414142d-6*gammaaa*t4*t81*t94+1.5
     7         95128677741414d-5*gammaaa*t11*t14*t81*t94-3.1902573554828
     8         28d-6*gammaaa*t12*t79*t8*t94-1.595128677741414d-6*gammaaa
     9         *t4*t54*t8*t94-7.97564338870707d-7*gammaaa*t10*t79*t94+7.
     :         97564338870707d-6*gammaaa*t21*t24*t58*t94+4.7853860332242
     ;         42d-6*gammaaa*t19*t22*t58*t94+2.392693016612121d-6*gammaa
     <         a*t18*t20*t58*t94+7.975643388707071d-7*gammaaa*t16*t58*t9
     =         4-7.97564338870707d-7*gammaaa*t15*t54*t94-3.9878216943535
     >         35d-6*gammaaa*t13*t14*t54*t94-3.190257355482828d-6*gammaa
     ?         a*t11*t12*t54*t94-1.329273898117845d-7*gammaaa*t17*t25*t9
     @         4-6.646369490589226d-7*gammaaa*t17*t23*t24*t94-5.31709559
     1         247138d-7*gammaaa*t17*t21*t22*t94-3.987821694353535d-7*ga
     2         mmaaa*t17*t19*t20*t94-2.65854779623569d-7*gammaaa*t16*t17
     3         *t18*t94+4.334248505425794d-4*t10*t56*t68*t9-8.6684970108
     4         51588d-4*t12*t54*t68*t9+2.8894990036171964d-4*t4*t56*t68*
     5         t8-4.334248505425794d-4*t10*t54*t68*t8-1.4447495018085982
     6         d-3*t21*t24*t58*t68-8.668497010851588d-4*t19*t22*t58*t68-
     7         4.334248505425794d-4*t18*t20*t58*t68-1.4447495018085982d-
     8         4*t16*t58*t68+1.4447495018085982d-4*t15*t56*t68+7.2237475
     9         0904299d-4*t13*t14*t56*t68+5.778998007234393d-4*t11*t12*t
     :         56*t68-1.4447495018085982d-4*t4*t54*t68-1.444749501808598
     ;         2d-3*t11*t14*t54*t68+7.22374750904299d-5*t17*t25*t68+3.61
     <         1873754521495d-4*t17*t23*t24*t68+2.8894990036171964d-4*t1
     =         7*t21*t22*t68+2.167124252712897d-4*t17*t19*t20*t68+1.4447
     >         495018085982d-4*t16*t17*t18*t68)-1.2407009817988002d+0*t4
     ?         1*(-8.126715947673364d-5*t10*t56*t74*t9+1.625343189534673
     @         d-4*t12*t54*t74*t9-5.417810631782243d-5*t4*t56*t74*t8+8.1
     1         26715947673364d-5*t10*t54*t74*t8+2.708905315891122d-4*t21
     2         *t24*t58*t74+1.625343189534673d-4*t19*t22*t58*t74+8.12671
     3         5947673364d-5*t18*t20*t58*t74+2.7089053158911214d-5*t16*t
     4         58*t74-2.7089053158911214d-5*t15*t56*t74-1.35445265794556
     5         1d-4*t13*t14*t56*t74-1.0835621263564485d-4*t11*t12*t56*t7
     6         4+2.7089053158911214d-5*t4*t54*t74+2.708905315891122d-4*t
     7         11*t14*t54*t74-1.3544526579455607d-5*t17*t25*t74-6.772263
     8         289727805d-5*t17*t23*t24*t74-5.417810631782243d-5*t17*t21
     9         *t22*t74-4.063357973836682d-5*t17*t19*t20*t74-2.708905315
     :         8911214d-5*t16*t17*t18*t74))*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363491002d-1*t27*(-7.975
     1         64338870707d-6*gammabb*t24*t37*t91*t95-3.190257355482828d
     2         -6*gammabb*t22*t36*t91*t95-7.97564338870707d-7*gammabb*t2
     3         0*t91*t95+1.5951286777414142d-6*gammabb*t4*t88*t95+1.5951
     4         28677741414d-5*gammabb*t14*t33*t88*t95+9.570772066448484d
     5         -6*gammabb*t12*t32*t88*t95+4.785386033224242d-6*gammabb*t
     6         10*t31*t88*t95-7.97564338870707d-6*gammabb*t14*t32*t86*t9
     7         5-3.190257355482828d-6*gammabb*t12*t31*t86*t95-7.97564338
     8         870707d-7*gammabb*t10*t86*t95+7.97564338870707d-6*gammabb
     9         *t24*t38*t66*t95+4.785386033224242d-6*gammabb*t22*t37*t66
     :         *t95+2.392693016612121d-6*gammabb*t20*t36*t66*t95+7.97564
     ;         3388707071d-7*gammabb*t16*t66*t95-1.595128677741414d-6*ga
     <         mmabb*t31*t4*t62*t95-3.987821694353535d-6*gammabb*t14*t34
     =         *t62*t95-3.190257355482828d-6*gammabb*t12*t33*t62*t95-2.3
     >         926930166121207d-6*gammabb*t10*t32*t62*t95-7.975643388707
     ?         07d-7*gammabb*t15*t62*t95-6.646369490589226d-7*gammabb*t2
     @         4*t35*t39*t95-5.31709559247138d-7*gammabb*t22*t35*t38*t95
     1         -3.987821694353535d-7*gammabb*t20*t35*t37*t95-2.658547796
     2         23569d-7*gammabb*t16*t35*t36*t95-1.329273898117845d-7*gam
     3         mabb*t25*t35*t95-1.4447495018085982d-3*t24*t38*t66*t71-8.
     4         668497010851588d-4*t22*t37*t66*t71-4.334248505425794d-4*t
     5         20*t36*t66*t71-1.4447495018085982d-4*t16*t66*t71+2.889499
     6         0036171964d-4*t31*t4*t64*t71+7.22374750904299d-4*t14*t34*
     7         t64*t71+5.778998007234393d-4*t12*t33*t64*t71+4.3342485054
     8         25794d-4*t10*t32*t64*t71+1.4447495018085982d-4*t15*t64*t7
     9         1-1.4447495018085982d-4*t4*t62*t71-1.4447495018085982d-3*
     :         t14*t33*t62*t71-8.668497010851588d-4*t12*t32*t62*t71-4.33
     ;         4248505425794d-4*t10*t31*t62*t71+3.611873754521495d-4*t24
     <         *t35*t39*t71+2.8894990036171964d-4*t22*t35*t38*t71+2.1671
     =         24252712897d-4*t20*t35*t37*t71+1.4447495018085982d-4*t16*
     >         t35*t36*t71+7.22374750904299d-5*t25*t35*t71)-1.2407009817
     ?         988002d+0*t45*(2.708905315891122d-4*t24*t38*t66*t76+1.625
     @         343189534673d-4*t22*t37*t66*t76+8.126715947673364d-5*t20*
     1         t36*t66*t76+2.7089053158911214d-5*t16*t66*t76-5.417810631
     2         782243d-5*t31*t4*t64*t76-1.354452657945561d-4*t14*t34*t64
     3         *t76-1.0835621263564485d-4*t12*t33*t64*t76-8.126715947673
     4         364d-5*t10*t32*t64*t76-2.7089053158911214d-5*t15*t64*t76+
     5         2.7089053158911214d-5*t4*t62*t76+2.708905315891122d-4*t14
     6         *t33*t62*t76+1.625343189534673d-4*t12*t32*t62*t76+8.12671
     7         5947673364d-5*t10*t31*t62*t76-6.772263289727805d-5*t24*t3
     8         5*t39*t76-5.417810631782243d-5*t22*t35*t38*t76-4.06335797
     9         3836682d-5*t20*t35*t37*t76-2.7089053158911214d-5*t16*t35*
     :         t36*t76-1.3544526579455607d-5*t25*t35*t76))*wght+Cmat3(iq
     ;         ,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491002d-1*t1*(t42*t5*(-1.1963465083060604d-6*t12*t56*t
     2         74*t9+2.990866270765151d-6*t14*t54*t74*t9-5.9817325415303
     3         01d-7*t10*t56*t74*t8+1.1963465083060604d-6*t12*t54*t74*t8
     4         -1.9939108471767675d-7*t4*t56*t74-1.9939108471767675d-6*t
     5         11*t14*t56*t74+2.9908662707651507d-7*t10*t54*t74)+t5*(2.9
     6         90866270765151d-6*t19*t24*t74*t84+1.1963465083060604d-6*t
     7         18*t22*t74*t84+2.9908662707651507d-7*t20*t74*t84-2.990866
     8         270765151d-6*t21*t24*t58*t74-1.7945197624590903d-6*t19*t2
     9         2*t58*t74-8.972598812295453d-7*t18*t20*t58*t74-2.99086627
     :         0765151d-7*t16*t58*t74+4.9847771179419187d-8*t17*t25*t74+
     ;         2.4923885589709593d-7*t17*t23*t24*t74+1.9939108471767675d
     <         -7*t17*t21*t22*t74+1.4954331353825753d-7*t17*t19*t20*t74+
     =         9.969554235883837d-8*t16*t17*t18*t74)-1.4721155704335503d
     >         -2*t56*t74*t75+8.126715947673364d-5*t50*t54/rhoa**8)*wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491002d-1*t27*(t28*(2.990866270765151d-6*t24*t37*t76*t
     2         91+1.1963465083060604d-6*t22*t36*t76*t91+2.99086627076515
     3         07d-7*t20*t76*t91-2.990866270765151d-6*t24*t38*t66*t76-1.
     4         7945197624590903d-6*t22*t37*t66*t76-8.972598812295453d-7*
     5         t20*t36*t66*t76-2.990866270765151d-7*t16*t66*t76+2.492388
     6         5589709593d-7*t24*t35*t39*t76+1.9939108471767675d-7*t22*t
     7         35*t38*t76+1.4954331353825753d-7*t20*t35*t37*t76+9.969554
     8         235883837d-8*t16*t35*t36*t76+4.9847771179419187d-8*t25*t3
     9         5*t76)-1.4721155704335503d-2*t64*t76*t77+t28*t46*(-1.9939
     :         108471767675d-7*t4*t64*t76-1.9939108471767675d-6*t14*t33*
     ;         t64*t76-1.1963465083060604d-6*t12*t32*t64*t76-5.981732541
     <         530301d-7*t10*t31*t64*t76+2.990866270765151d-6*t14*t32*t6
     =         2*t76+1.1963465083060604d-6*t12*t31*t62*t76+2.99086627076
     >         51507d-7*t10*t62*t76)+8.126715947673364d-5*t51*t62/rhob**
     ?         8)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammaaa*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhoa**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhoa**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammaaa*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammaaa*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammaaa*t2*t26*t27*t6-9.814103802890335d-3*gammaaa*t13*t26
     3         *t27-4.9070519014451675d-2*gammaaa*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammaaa*t15*t23*t27-4.9070519014451675d-2*
     5         gammaaa*t15*t21*t22*t27-3.925641521156134d-2*gammaaa*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammaaa*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammaaa*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = 1/rhoa**6.666666666666666d-1
            t31 = gammaaa**2
            t32 = 1/t5**4
            t33 = 1/rhoa**7.333333333333333d+0
            t34 = 1/t5**3
            t35 = 1/rhoa**4.666666666666667d+0
            t36 = exp(-7.360577852167751d-3*gammaaa*t3)
            t37 = 1.4394018910905823d-1*gammaaa*t10*t26*t35*t9-7.7053306
     1         76312523d-4*t10*t31*t33*t34*t9+1.9263326690781307d-3*t12*
     2         t31*t32*t33*t9+1.0795514183179368d-1*gammaaa*t26*t35*t7*t
     3         8-5.778998007234393d-4*t31*t33*t34*t7*t8+5.77899800723439
     4         3d-4*t31*t32*t33*t6*t8+1.1557996014468785d-3*t10*t31*t32*
     5         t33*t7+7.197009455452912d-2*gammaaa*t2*t26*t35*t6-3.85266
     6         53381562614d-4*t2*t31*t33*t34*t6+1.9263326690781307d-3*t1
     7         9*t22*t31*t33*t36+1.1557996014468785d-3*t17*t20*t31*t33*t
     8         36+5.778998007234393d-4*t16*t18*t31*t33*t36+1.92633266907
     9         81307d-4*t14*t31*t33*t36+3.598504727726456d-2*gammaaa*t13
     :         *t26*t35+1.799252363863228d-1*gammaaa*t11*t12*t26*t35+3.5
     ;         98504727726456d-2*gammaaa*t15*t23*t35+1.799252363863228d-
     <         1*gammaaa*t15*t21*t22*t35+1.4394018910905823d-1*gammaaa*t
     =         15*t19*t20*t35+1.0795514183179368d-1*gammaaa*t15*t17*t18*
     >         t35+7.197009455452912d-2*gammaaa*t14*t15*t16*t35-1.926332
     ?         6690781307d-4*t13*t31*t33*t34-9.631663345390654d-4*t11*t1
     @         2*t31*t33*t34+1.9263326690781307d-4*t2*t31*t32*t33-9.6316
     1         63345390653d-5*t15*t23*t31*t33-4.815831672695327d-4*t15*t
     2         21*t22*t31*t33-3.8526653381562614d-4*t15*t19*t20*t31*t33-
     3         2.8894990036171964d-4*t15*t17*t18*t31*t33-1.9263326690781
     4         307d-4*t14*t15*t16*t31*t33
            t38 = 1/rhoa**6.333333333333333d+0
            t39 = 2.8894990036171964d-4*gammaaa*t10*t34*t38*t9-7.2237475
     1         0904299d-4*gammaaa*t12*t32*t38*t9-3.925641521156134d-2*t1
     2         0*t26*t27*t9+2.167124252712897d-4*gammaaa*t34*t38*t7*t8-2
     3         .9442311408671007d-2*t26*t27*t7*t8-2.167124252712897d-4*g
     4         ammaaa*t32*t38*t6*t8-4.334248505425794d-4*gammaaa*t10*t32
     5         *t38*t7+1.4447495018085982d-4*gammaaa*t2*t34*t38*t6-1.962
     6         820760578067d-2*t2*t26*t27*t6-7.22374750904299d-4*gammaaa
     7         *t19*t22*t36*t38-4.334248505425794d-4*gammaaa*t17*t20*t36
     8         *t38-2.167124252712897d-4*gammaaa*t16*t18*t36*t38-7.22374
     9         750904299d-5*gammaaa*t14*t36*t38+7.22374750904299d-5*gamm
     :         aaa*t13*t34*t38+3.611873754521495d-4*gammaaa*t11*t12*t34*
     ;         t38-7.22374750904299d-5*gammaaa*t2*t32*t38+3.611873754521
     <         495d-5*gammaaa*t15*t23*t38+1.8059368772607476d-4*gammaaa*
     =         t15*t21*t22*t38+1.4447495018085982d-4*gammaaa*t15*t19*t20
     >         *t38+1.0835621263564485d-4*gammaaa*t15*t17*t18*t38+7.2237
     ?         4750904299d-5*gammaaa*t14*t15*t16*t38-9.814103802890335d-
     @         3*t13*t26*t27-4.9070519014451675d-2*t11*t12*t26*t27-9.814
     1         103802890335d-3*t15*t23*t27-4.9070519014451675d-2*t15*t21
     2         *t22*t27-3.925641521156134d-2*t15*t19*t20*t27-2.944231140
     3         8671007d-2*t15*t17*t18*t27-1.962820760578067d-2*t14*t15*t
     4         16*t27
            t40 = 1.4721155704335503d-2*t10*t26*t3*t9+1.1040866778251626
     1         d-2*t26*t3*t7*t8+7.360577852167751d-3*t2*t26*t3*t6+3.6802
     2         889260838756d-3*t13*t26*t3+1.8401444630419378d-2*t11*t12*
     3         t26*t3+3.6802889260838756d-3*t15*t23*t3+1.840144463041937
     4         8d-2*t15*t21*t22*t3+1.4721155704335503d-2*t15*t19*t20*t3+
     5         1.1040866778251626d-2*t15*t17*t18*t3+7.360577852167751d-3
     6         *t14*t15*t16*t3
            t41 = 1/rhoa**5.333333333333333d+0
            t42 = 2.708905315891122d-4*t12*t26*t3*t9+8.126715947673364d-
     1         5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26*t3*t7+2.70890
     2         53158911214d-5*t2*t26*t3
            t43 = gammaaa**3
            t44 = 1/t5**6
            t45 = 1/rhoa**11
            t46 = 1/t5**5
            t47 = 1/rhoa**8.333333333333334d+0
            t48 = 1/rhoa**5.666666666666667d+0
            t49 = exp(-1.1040866778251626d-2*gammaaa*t3)
            t50 = 1/rhoa**10
            t51 = 1/rhoa**9
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RA) = -9.305257363491002d-1*t1*t28*wght-1.2407009
     1         817988002d+0*t24*t25*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Amat2(iq,D2_RA_RA) = -9.305257363491002d-1*t1*t37*wght-4.135
     1         669939329334d-1*t24*t30*wght-2.4814019635976003d+0*t25*t2
     2         8*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -1.2407009817988002d+0*t25*t40*wght-9.
     1         305257363491002d-1*t1*t39*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349100
     1         2d-1*t1*(t26*t3*t42-7.360577852167751d-3*t29*t34*t41+t3*(
     2         2.708905315891122d-4*t19*t22*t3*t36+1.625343189534673d-4*
     3         t17*t20*t3*t36+8.126715947673364d-5*t16*t18*t3*t36+2.7089
     4         053158911214d-5*t14*t3*t36-1.3544526579455607d-5*t15*t23*
     5         t3-6.772263289727805d-5*t15*t21*t22*t3-5.417810631782243d
     6         -5*t15*t19*t20*t3-4.063357973836682d-5*t15*t17*t18*t3-2.7
     7         089053158911214d-5*t14*t15*t16*t3))*wght
            Amat3(iq,D3_RA_RA_RA) = -9.305257363491002d-1*t1*(-6.7172088
     1         25089385d-1*gammaaa*t10*t26*t48*t9+8.475863743943775d-3*t
     2         10*t31*t34*t47*t9-2.1189659359859436d-2*t12*t31*t32*t47*t
     3         9+1.1343137263938945d-4*t12*t43*t45*t46*t9-2.268627452787
     4         7887d-5*t10*t32*t43*t45*t9-5.037906618817038d-1*gammaaa*t
     5         26*t48*t7*t8+6.356897807957831d-3*t31*t34*t47*t7*t8-1.701
     6         4705895908416d-5*t32*t43*t45*t7*t8-6.356897807957831d-3*t
     7         31*t32*t47*t6*t8+3.4029411791816827d-5*t43*t45*t46*t6*t8-
     8         5.671568631969471d-6*t43*t44*t45*t8-1.2713795615915663d-2
     9         *t10*t31*t32*t47*t7+6.805882358363365d-5*t10*t43*t45*t46*
     :         t7-5.671568631969472d-5*t12*t43*t44*t45*t7-3.358604412544
     ;         6923d-1*gammaaa*t2*t26*t48*t6+4.237931871971887d-3*t2*t31
     <         *t34*t47*t6-2.2686274527877887d-5*t10*t43*t44*t45*t6-1.13
     =         43137263938943d-5*t2*t32*t43*t45*t6-5.671568631969472d-5*
     >         t17*t22*t43*t45*t49-2.2686274527877887d-5*t16*t20*t43*t45
     ?         *t49-5.671568631969471d-6*t18*t43*t45*t49-1.6793022062723
     @         462d-1*gammaaa*t13*t26*t48-8.396511031361731d-1*gammaaa*t
     1         11*t12*t26*t48-1.6793022062723462d-1*gammaaa*t15*t23*t48-
     2         8.396511031361731d-1*gammaaa*t15*t21*t22*t48-6.7172088250
     3         89385d-1*gammaaa*t15*t19*t20*t48-5.037906618817038d-1*gam
     4         maaa*t15*t17*t18*t48-3.3586044125446923d-1*gammaaa*t14*t1
     5         5*t16*t48-2.1189659359859436d-2*t19*t22*t31*t36*t47-1.271
     6         3795615915663d-2*t17*t20*t31*t36*t47-6.356897807957831d-3
     7         *t16*t18*t31*t36*t47-2.1189659359859436d-3*t14*t31*t36*t4
     8         7+2.1189659359859436d-3*t13*t31*t34*t47+1.059482967992971
     9         8d-2*t11*t12*t31*t34*t47-2.1189659359859436d-3*t2*t31*t32
     :         *t47+1.0594829679929718d-3*t15*t23*t31*t47+5.297414839964
     ;         859d-3*t15*t21*t22*t31*t47+4.237931871971887d-3*t15*t19*t
     <         20*t31*t47+3.1784489039789154d-3*t15*t17*t18*t31*t47+2.11
     =         89659359859436d-3*t14*t15*t16*t31*t47+1.1343137263938943d
     >         -5*t2*t43*t45*t46+5.671568631969472d-5*t19*t22*t36*t43*t4
     ?         5+3.4029411791816827d-5*t17*t20*t36*t43*t45+1.70147058959
     @         08414d-5*t16*t18*t36*t43*t45+5.671568631969471d-6*t14*t36
     1         *t43*t45-5.671568631969471d-6*t13*t32*t43*t45-2.835784315
     2         984736d-5*t11*t12*t32*t43*t45-9.452614386615786d-7*t15*t2
     3         3*t43*t45-4.726307193307893d-6*t15*t21*t22*t43*t45-3.7810
     4         457546463144d-6*t15*t19*t20*t43*t45-2.8357843159847357d-6
     5         *t15*t17*t18*t43*t45-1.8905228773231572d-6*t14*t15*t16*t4
     6         3*t45)*wght-3.7221029453964005d+0*t25*t37*wght-1.24070098
     7         17988002d+0*t28*t30*wght+2.7571132928862224d-1*t24*wght/r
     8         hoa**1.6666666666666669d+0+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363491002d-1*t1*(-4.253676
     1         4739771044d-5*t12*t31*t46*t50*t9+8.507352947954206d-6*t10
     2         *t31*t32*t50*t9+1.4394018910905823d-1*t10*t26*t35*t9-2.60
     3         0549103255476d-3*gammaaa*t10*t33*t34*t9+6.50137275813869d
     4         -3*gammaaa*t12*t32*t33*t9+6.380514710965656d-6*t31*t32*t5
     5         0*t7*t8+1.0795514183179368d-1*t26*t35*t7*t8-1.95041182744
     6         16074d-3*gammaaa*t33*t34*t7*t8-1.2761029421931314d-5*t31*
     7         t46*t50*t6*t8+1.9504118274416074d-3*gammaaa*t32*t33*t6*t8
     8         +2.1268382369885522d-6*t31*t44*t50*t8-2.552205884386262d-
     9         5*t10*t31*t46*t50*t7+2.1268382369885516d-5*t12*t31*t44*t5
     :         0*t7+3.900823654883215d-3*gammaaa*t10*t32*t33*t7+8.507352
     ;         947954209d-6*t10*t31*t44*t50*t6+4.253676473977103d-6*t2*t
     <         31*t32*t50*t6+7.197009455452912d-2*t2*t26*t35*t6-1.300274
     =         551627738d-3*gammaaa*t2*t33*t34*t6+2.1268382369885516d-5*
     >         t17*t22*t31*t49*t50+8.507352947954209d-6*t16*t20*t31*t49*
     ?         t50+2.1268382369885522d-6*t18*t31*t49*t50-4.2536764739771
     @         03d-6*t2*t31*t46*t50-2.1268382369885522d-5*t19*t22*t31*t3
     1         6*t50-1.276102942193131d-5*t17*t20*t31*t36*t50-6.38051471
     2         0965657d-6*t16*t18*t31*t36*t50-2.1268382369885516d-6*t14*
     3         t31*t36*t50+2.1268382369885516d-6*t13*t31*t32*t50+1.06341
     4         91184942758d-5*t11*t12*t31*t32*t50+3.54473039498092d-7*t1
     5         5*t23*t31*t50+1.77236519749046d-6*t15*t21*t22*t31*t50+1.4
     6         178921579923678d-6*t15*t19*t20*t31*t50+1.0634191184942761
     7         d-6*t15*t17*t18*t31*t50+7.08946078996184d-7*t14*t15*t16*t
     8         31*t50+6.50137275813869d-3*gammaaa*t19*t22*t33*t36+3.9008
     9         23654883215d-3*gammaaa*t17*t20*t33*t36+1.9504118274416074
     :         d-3*gammaaa*t16*t18*t33*t36+6.501372758138692d-4*gammaaa*
     ;         t14*t33*t36+3.598504727726456d-2*t13*t26*t35+1.7992523638
     <         63228d-1*t11*t12*t26*t35+3.598504727726456d-2*t15*t23*t35
     =         +1.799252363863228d-1*t15*t21*t22*t35+1.4394018910905823d
     >         -1*t15*t19*t20*t35+1.0795514183179368d-1*t15*t17*t18*t35+
     ?         7.197009455452912d-2*t14*t15*t16*t35-6.501372758138692d-4
     @         *gammaaa*t13*t33*t34-3.250686379069345d-3*gammaaa*t11*t12
     1         *t33*t34+6.501372758138692d-4*gammaaa*t2*t32*t33-3.250686
     2         379069346d-4*gammaaa*t15*t23*t33-1.6253431895346726d-3*ga
     3         mmaaa*t15*t21*t22*t33-1.300274551627738d-3*gammaaa*t15*t1
     4         9*t20*t33-9.752059137208037d-4*gammaaa*t15*t17*t18*t33-6.
     5         501372758138692d-4*gammaaa*t14*t15*t16*t33)*wght-4.135669
     6         939329334d-1*t30*t40*wght-2.4814019635976003d+0*t25*t39*w
     7         ght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363491002d-1*t1*(1.595128
     1         677741414d-5*gammaaa*t12*t46*t51*t9-3.190257355482828d-6*
     2         gammaaa*t10*t32*t51*t9+5.778998007234393d-4*t10*t34*t38*t
     3         9-1.4447495018085982d-3*t12*t32*t38*t9-2.3926930166121207
     4         d-6*gammaaa*t32*t51*t7*t8+4.334248505425794d-4*t34*t38*t7
     5         *t8+4.785386033224242d-6*gammaaa*t46*t51*t6*t8-4.33424850
     6         5425794d-4*t32*t38*t6*t8-7.97564338870707d-7*gammaaa*t44*
     7         t51*t8+9.570772066448484d-6*gammaaa*t10*t46*t51*t7-7.9756
     8         4338870707d-6*gammaaa*t12*t44*t51*t7-8.668497010851588d-4
     9         *t10*t32*t38*t7-3.190257355482828d-6*gammaaa*t10*t44*t51*
     :         t6-1.595128677741414d-6*gammaaa*t2*t32*t51*t6+2.889499003
     ;         6171964d-4*t2*t34*t38*t6-7.97564338870707d-6*gammaaa*t17*
     <         t22*t49*t51-3.190257355482828d-6*gammaaa*t16*t20*t49*t51-
     =         7.97564338870707d-7*gammaaa*t18*t49*t51+1.595128677741414
     >         2d-6*gammaaa*t2*t46*t51+7.97564338870707d-6*gammaaa*t19*t
     ?         22*t36*t51+4.785386033224242d-6*gammaaa*t17*t20*t36*t51+2
     @         .392693016612121d-6*gammaaa*t16*t18*t36*t51+7.97564338870
     1         7071d-7*gammaaa*t14*t36*t51-7.97564338870707d-7*gammaaa*t
     2         13*t32*t51-3.987821694353535d-6*gammaaa*t11*t12*t32*t51-1
     3         .329273898117845d-7*gammaaa*t15*t23*t51-6.646369490589226
     4         d-7*gammaaa*t15*t21*t22*t51-5.31709559247138d-7*gammaaa*t
     5         15*t19*t20*t51-3.987821694353535d-7*gammaaa*t15*t17*t18*t
     6         51-2.65854779623569d-7*gammaaa*t14*t15*t16*t51-1.44474950
     7         18085982d-3*t19*t22*t36*t38-8.668497010851588d-4*t17*t20*
     8         t36*t38-4.334248505425794d-4*t16*t18*t36*t38-1.4447495018
     9         085982d-4*t14*t36*t38+1.4447495018085982d-4*t13*t34*t38+7
     :         .22374750904299d-4*t11*t12*t34*t38-1.4447495018085982d-4*
     ;         t2*t32*t38+7.22374750904299d-5*t15*t23*t38+3.611873754521
     <         495d-4*t15*t21*t22*t38+2.8894990036171964d-4*t15*t19*t20*
     =         t38+2.167124252712897d-4*t15*t17*t18*t38+1.44474950180859
     >         82d-4*t14*t15*t16*t38)*wght-1.2407009817988002d+0*t25*(-1
     ?         .0835621263564485d-4*t10*t34*t41*t9+2.708905315891122d-4*
     @         t12*t32*t41*t9-8.126715947673364d-5*t34*t41*t7*t8+8.12671
     1         5947673364d-5*t32*t41*t6*t8+1.625343189534673d-4*t10*t32*
     2         t41*t7-5.417810631782243d-5*t2*t34*t41*t6+2.7089053158911
     3         22d-4*t19*t22*t36*t41+1.625343189534673d-4*t17*t20*t36*t4
     4         1+8.126715947673364d-5*t16*t18*t36*t41+2.7089053158911214
     5         d-5*t14*t36*t41-2.7089053158911214d-5*t13*t34*t41-1.35445
     6         2657945561d-4*t11*t12*t34*t41+2.7089053158911214d-5*t2*t3
     7         2*t41-1.3544526579455607d-5*t15*t23*t41-6.772263289727805
     8         d-5*t15*t21*t22*t41-5.417810631782243d-5*t15*t19*t20*t41-
     9         4.063357973836682d-5*t15*t17*t18*t41-2.7089053158911214d-
     :         5*t14*t15*t16*t41)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491002d-1*t1*(t26*t3*(-1.9939108471767675d-6*t12*t34*t
     2         41*t9-5.981732541530301d-7*t34*t41*t6*t8+2.99086627076515
     3         07d-7*t32*t41*t8-1.1963465083060604d-6*t10*t34*t41*t7+2.9
     4         90866270765151d-6*t12*t32*t41*t7+1.1963465083060604d-6*t1
     5         0*t32*t41*t6-1.9939108471767675d-7*t2*t34*t41)+t3*(2.9908
     6         66270765151d-6*t17*t22*t41*t49+1.1963465083060604d-6*t16*
     7         t20*t41*t49+2.9908662707651507d-7*t18*t41*t49-2.990866270
     8         765151d-6*t19*t22*t36*t41-1.7945197624590903d-6*t17*t20*t
     9         36*t41-8.972598812295453d-7*t16*t18*t36*t41-2.99086627076
     :         5151d-7*t14*t36*t41+4.9847771179419187d-8*t15*t23*t41+2.4
     ;         923885589709593d-7*t15*t21*t22*t41+1.9939108471767675d-7*
     <         t15*t19*t20*t41+1.4954331353825753d-7*t15*t17*t18*t41+9.9
     =         69554235883837d-8*t14*t15*t16*t41)-1.4721155704335503d-2*
     >         t34*t41*t42+8.126715947673364d-5*t29*t32/rhoa**8)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = param(3)
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 3.6802889260838756d-3*gammabb*t3
            t5 = t4+1.0d+0
            t6 = 1.0d+0-1.0d+0/t5
            t7 = t6**2.0d+0
            t8 = param(4)
            t9 = t6**3.0d+0
            t10 = param(5)
            t11 = t6**4.0d+0
            t12 = param(6)
            t13 = param(2)
            t14 = param(9)
            t15 = exp(-t4)
            t16 = 1.0d+0-t15
            t17 = t16**2.0d+0
            t18 = param(10)
            t19 = t16**3.0d+0
            t20 = param(11)
            t21 = t16**4.0d+0
            t22 = param(12)
            t23 = param(8)
            t24 = t8*t9+t2*t7+t12*t6**5.0d+0+t13*t6+t16*t23+t16**5.0d+0*
     1         t22+t20*t21+t18*t19+t14*t17+t10*t11+param(7)+param(1)
            t25 = rhob**3.333333333333333d-1
            t26 = 1/t5**2
            t27 = 1/rhob**3.6666666666666664d+0
            t28 = -3.925641521156134d-2*gammabb*t10*t26*t27*t9-2.9442311
     1         408671007d-2*gammabb*t26*t27*t7*t8-1.962820760578067d-2*g
     2         ammabb*t2*t26*t27*t6-9.814103802890335d-3*gammabb*t13*t26
     3         *t27-4.9070519014451675d-2*gammabb*t11*t12*t26*t27-9.8141
     4         03802890335d-3*gammabb*t15*t23*t27-4.9070519014451675d-2*
     5         gammabb*t15*t21*t22*t27-3.925641521156134d-2*gammabb*t15*
     6         t19*t20*t27-2.9442311408671007d-2*gammabb*t15*t17*t18*t27
     7         -1.962820760578067d-2*gammabb*t14*t15*t16*t27
            t29 = 1.4721155704335503d-2*t10*t9+1.1040866778251626d-2*t7*
     1         t8+7.360577852167751d-3*t2*t6+3.6802889260838756d-3*t13+1
     2         .8401444630419378d-2*t11*t12
            t30 = 1/rhob**6.666666666666666d-1
            t31 = gammabb**2
            t32 = 1/t5**4
            t33 = 1/rhob**7.333333333333333d+0
            t34 = 1/t5**3
            t35 = 1/rhob**4.666666666666667d+0
            t36 = exp(-7.360577852167751d-3*gammabb*t3)
            t37 = 1.4394018910905823d-1*gammabb*t10*t26*t35*t9-7.7053306
     1         76312523d-4*t10*t31*t33*t34*t9+1.9263326690781307d-3*t12*
     2         t31*t32*t33*t9+1.0795514183179368d-1*gammabb*t26*t35*t7*t
     3         8-5.778998007234393d-4*t31*t33*t34*t7*t8+5.77899800723439
     4         3d-4*t31*t32*t33*t6*t8+1.1557996014468785d-3*t10*t31*t32*
     5         t33*t7+7.197009455452912d-2*gammabb*t2*t26*t35*t6-3.85266
     6         53381562614d-4*t2*t31*t33*t34*t6+1.9263326690781307d-3*t1
     7         9*t22*t31*t33*t36+1.1557996014468785d-3*t17*t20*t31*t33*t
     8         36+5.778998007234393d-4*t16*t18*t31*t33*t36+1.92633266907
     9         81307d-4*t14*t31*t33*t36+3.598504727726456d-2*gammabb*t13
     :         *t26*t35+1.799252363863228d-1*gammabb*t11*t12*t26*t35+3.5
     ;         98504727726456d-2*gammabb*t15*t23*t35+1.799252363863228d-
     <         1*gammabb*t15*t21*t22*t35+1.4394018910905823d-1*gammabb*t
     =         15*t19*t20*t35+1.0795514183179368d-1*gammabb*t15*t17*t18*
     >         t35+7.197009455452912d-2*gammabb*t14*t15*t16*t35-1.926332
     ?         6690781307d-4*t13*t31*t33*t34-9.631663345390654d-4*t11*t1
     @         2*t31*t33*t34+1.9263326690781307d-4*t2*t31*t32*t33-9.6316
     1         63345390653d-5*t15*t23*t31*t33-4.815831672695327d-4*t15*t
     2         21*t22*t31*t33-3.8526653381562614d-4*t15*t19*t20*t31*t33-
     3         2.8894990036171964d-4*t15*t17*t18*t31*t33-1.9263326690781
     4         307d-4*t14*t15*t16*t31*t33
            t38 = 1/rhob**6.333333333333333d+0
            t39 = 2.8894990036171964d-4*gammabb*t10*t34*t38*t9-7.2237475
     1         0904299d-4*gammabb*t12*t32*t38*t9-3.925641521156134d-2*t1
     2         0*t26*t27*t9+2.167124252712897d-4*gammabb*t34*t38*t7*t8-2
     3         .9442311408671007d-2*t26*t27*t7*t8-2.167124252712897d-4*g
     4         ammabb*t32*t38*t6*t8-4.334248505425794d-4*gammabb*t10*t32
     5         *t38*t7+1.4447495018085982d-4*gammabb*t2*t34*t38*t6-1.962
     6         820760578067d-2*t2*t26*t27*t6-7.22374750904299d-4*gammabb
     7         *t19*t22*t36*t38-4.334248505425794d-4*gammabb*t17*t20*t36
     8         *t38-2.167124252712897d-4*gammabb*t16*t18*t36*t38-7.22374
     9         750904299d-5*gammabb*t14*t36*t38+7.22374750904299d-5*gamm
     :         abb*t13*t34*t38+3.611873754521495d-4*gammabb*t11*t12*t34*
     ;         t38-7.22374750904299d-5*gammabb*t2*t32*t38+3.611873754521
     <         495d-5*gammabb*t15*t23*t38+1.8059368772607476d-4*gammabb*
     =         t15*t21*t22*t38+1.4447495018085982d-4*gammabb*t15*t19*t20
     >         *t38+1.0835621263564485d-4*gammabb*t15*t17*t18*t38+7.2237
     ?         4750904299d-5*gammabb*t14*t15*t16*t38-9.814103802890335d-
     @         3*t13*t26*t27-4.9070519014451675d-2*t11*t12*t26*t27-9.814
     1         103802890335d-3*t15*t23*t27-4.9070519014451675d-2*t15*t21
     2         *t22*t27-3.925641521156134d-2*t15*t19*t20*t27-2.944231140
     3         8671007d-2*t15*t17*t18*t27-1.962820760578067d-2*t14*t15*t
     4         16*t27
            t40 = 1.4721155704335503d-2*t10*t26*t3*t9+1.1040866778251626
     1         d-2*t26*t3*t7*t8+7.360577852167751d-3*t2*t26*t3*t6+3.6802
     2         889260838756d-3*t13*t26*t3+1.8401444630419378d-2*t11*t12*
     3         t26*t3+3.6802889260838756d-3*t15*t23*t3+1.840144463041937
     4         8d-2*t15*t21*t22*t3+1.4721155704335503d-2*t15*t19*t20*t3+
     5         1.1040866778251626d-2*t15*t17*t18*t3+7.360577852167751d-3
     6         *t14*t15*t16*t3
            t41 = 1/rhob**5.333333333333333d+0
            t42 = 2.708905315891122d-4*t12*t26*t3*t9+8.126715947673364d-
     1         5*t26*t3*t6*t8+1.625343189534673d-4*t10*t26*t3*t7+2.70890
     2         53158911214d-5*t2*t26*t3
            t43 = gammabb**3
            t44 = 1/t5**6
            t45 = 1/rhob**11
            t46 = 1/t5**5
            t47 = 1/rhob**8.333333333333334d+0
            t48 = 1/rhob**5.666666666666667d+0
            t49 = exp(-1.1040866778251626d-2*gammabb*t3)
            t50 = 1/rhob**10
            t51 = 1/rhob**9
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t24*wght
            Amat(iq,D1_RB) = -9.305257363491002d-1*t1*t28*wght-1.2407009
     1         817988002d+0*t24*t25*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491002d-1*t1*(t
     1         26*t29*t3+(3.6802889260838756d-3*t15*t23+1.84014446304193
     2         78d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.104
     3         0866778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15
     4         *t16)*t3)*wght
            Amat2(iq,D2_RB_RB) = -9.305257363491002d-1*t1*t37*wght-4.135
     1         669939329334d-1*t24*t30*wght-2.4814019635976003d+0*t25*t2
     2         8*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -1.2407009817988002d+0*t25*t40*wght-9.
     1         305257363491002d-1*t1*t39*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349100
     1         2d-1*t1*(t26*t3*t42-7.360577852167751d-3*t29*t34*t41+t3*(
     2         2.708905315891122d-4*t19*t22*t3*t36+1.625343189534673d-4*
     3         t17*t20*t3*t36+8.126715947673364d-5*t16*t18*t3*t36+2.7089
     4         053158911214d-5*t14*t3*t36-1.3544526579455607d-5*t15*t23*
     5         t3-6.772263289727805d-5*t15*t21*t22*t3-5.417810631782243d
     6         -5*t15*t19*t20*t3-4.063357973836682d-5*t15*t17*t18*t3-2.7
     7         089053158911214d-5*t14*t15*t16*t3))*wght
            Amat3(iq,D3_RB_RB_RB) = -9.305257363491002d-1*t1*(-6.7172088
     1         25089385d-1*gammabb*t10*t26*t48*t9+8.475863743943775d-3*t
     2         10*t31*t34*t47*t9-2.1189659359859436d-2*t12*t31*t32*t47*t
     3         9+1.1343137263938945d-4*t12*t43*t45*t46*t9-2.268627452787
     4         7887d-5*t10*t32*t43*t45*t9-5.037906618817038d-1*gammabb*t
     5         26*t48*t7*t8+6.356897807957831d-3*t31*t34*t47*t7*t8-1.701
     6         4705895908416d-5*t32*t43*t45*t7*t8-6.356897807957831d-3*t
     7         31*t32*t47*t6*t8+3.4029411791816827d-5*t43*t45*t46*t6*t8-
     8         5.671568631969471d-6*t43*t44*t45*t8-1.2713795615915663d-2
     9         *t10*t31*t32*t47*t7+6.805882358363365d-5*t10*t43*t45*t46*
     :         t7-5.671568631969472d-5*t12*t43*t44*t45*t7-3.358604412544
     ;         6923d-1*gammabb*t2*t26*t48*t6+4.237931871971887d-3*t2*t31
     <         *t34*t47*t6-2.2686274527877887d-5*t10*t43*t44*t45*t6-1.13
     =         43137263938943d-5*t2*t32*t43*t45*t6-5.671568631969472d-5*
     >         t17*t22*t43*t45*t49-2.2686274527877887d-5*t16*t20*t43*t45
     ?         *t49-5.671568631969471d-6*t18*t43*t45*t49-1.6793022062723
     @         462d-1*gammabb*t13*t26*t48-8.396511031361731d-1*gammabb*t
     1         11*t12*t26*t48-1.6793022062723462d-1*gammabb*t15*t23*t48-
     2         8.396511031361731d-1*gammabb*t15*t21*t22*t48-6.7172088250
     3         89385d-1*gammabb*t15*t19*t20*t48-5.037906618817038d-1*gam
     4         mabb*t15*t17*t18*t48-3.3586044125446923d-1*gammabb*t14*t1
     5         5*t16*t48-2.1189659359859436d-2*t19*t22*t31*t36*t47-1.271
     6         3795615915663d-2*t17*t20*t31*t36*t47-6.356897807957831d-3
     7         *t16*t18*t31*t36*t47-2.1189659359859436d-3*t14*t31*t36*t4
     8         7+2.1189659359859436d-3*t13*t31*t34*t47+1.059482967992971
     9         8d-2*t11*t12*t31*t34*t47-2.1189659359859436d-3*t2*t31*t32
     :         *t47+1.0594829679929718d-3*t15*t23*t31*t47+5.297414839964
     ;         859d-3*t15*t21*t22*t31*t47+4.237931871971887d-3*t15*t19*t
     <         20*t31*t47+3.1784489039789154d-3*t15*t17*t18*t31*t47+2.11
     =         89659359859436d-3*t14*t15*t16*t31*t47+1.1343137263938943d
     >         -5*t2*t43*t45*t46+5.671568631969472d-5*t19*t22*t36*t43*t4
     ?         5+3.4029411791816827d-5*t17*t20*t36*t43*t45+1.70147058959
     @         08414d-5*t16*t18*t36*t43*t45+5.671568631969471d-6*t14*t36
     1         *t43*t45-5.671568631969471d-6*t13*t32*t43*t45-2.835784315
     2         984736d-5*t11*t12*t32*t43*t45-9.452614386615786d-7*t15*t2
     3         3*t43*t45-4.726307193307893d-6*t15*t21*t22*t43*t45-3.7810
     4         457546463144d-6*t15*t19*t20*t43*t45-2.8357843159847357d-6
     5         *t15*t17*t18*t43*t45-1.8905228773231572d-6*t14*t15*t16*t4
     6         3*t45)*wght-3.7221029453964005d+0*t25*t37*wght-1.24070098
     7         17988002d+0*t28*t30*wght+2.7571132928862224d-1*t24*wght/r
     8         hob**1.6666666666666669d+0+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363491002d-1*t1*(-4.253676
     1         4739771044d-5*t12*t31*t46*t50*t9+8.507352947954206d-6*t10
     2         *t31*t32*t50*t9+1.4394018910905823d-1*t10*t26*t35*t9-2.60
     3         0549103255476d-3*gammabb*t10*t33*t34*t9+6.50137275813869d
     4         -3*gammabb*t12*t32*t33*t9+6.380514710965656d-6*t31*t32*t5
     5         0*t7*t8+1.0795514183179368d-1*t26*t35*t7*t8-1.95041182744
     6         16074d-3*gammabb*t33*t34*t7*t8-1.2761029421931314d-5*t31*
     7         t46*t50*t6*t8+1.9504118274416074d-3*gammabb*t32*t33*t6*t8
     8         +2.1268382369885522d-6*t31*t44*t50*t8-2.552205884386262d-
     9         5*t10*t31*t46*t50*t7+2.1268382369885516d-5*t12*t31*t44*t5
     :         0*t7+3.900823654883215d-3*gammabb*t10*t32*t33*t7+8.507352
     ;         947954209d-6*t10*t31*t44*t50*t6+4.253676473977103d-6*t2*t
     <         31*t32*t50*t6+7.197009455452912d-2*t2*t26*t35*t6-1.300274
     =         551627738d-3*gammabb*t2*t33*t34*t6+2.1268382369885516d-5*
     >         t17*t22*t31*t49*t50+8.507352947954209d-6*t16*t20*t31*t49*
     ?         t50+2.1268382369885522d-6*t18*t31*t49*t50-4.2536764739771
     @         03d-6*t2*t31*t46*t50-2.1268382369885522d-5*t19*t22*t31*t3
     1         6*t50-1.276102942193131d-5*t17*t20*t31*t36*t50-6.38051471
     2         0965657d-6*t16*t18*t31*t36*t50-2.1268382369885516d-6*t14*
     3         t31*t36*t50+2.1268382369885516d-6*t13*t31*t32*t50+1.06341
     4         91184942758d-5*t11*t12*t31*t32*t50+3.54473039498092d-7*t1
     5         5*t23*t31*t50+1.77236519749046d-6*t15*t21*t22*t31*t50+1.4
     6         178921579923678d-6*t15*t19*t20*t31*t50+1.0634191184942761
     7         d-6*t15*t17*t18*t31*t50+7.08946078996184d-7*t14*t15*t16*t
     8         31*t50+6.50137275813869d-3*gammabb*t19*t22*t33*t36+3.9008
     9         23654883215d-3*gammabb*t17*t20*t33*t36+1.9504118274416074
     :         d-3*gammabb*t16*t18*t33*t36+6.501372758138692d-4*gammabb*
     ;         t14*t33*t36+3.598504727726456d-2*t13*t26*t35+1.7992523638
     <         63228d-1*t11*t12*t26*t35+3.598504727726456d-2*t15*t23*t35
     =         +1.799252363863228d-1*t15*t21*t22*t35+1.4394018910905823d
     >         -1*t15*t19*t20*t35+1.0795514183179368d-1*t15*t17*t18*t35+
     ?         7.197009455452912d-2*t14*t15*t16*t35-6.501372758138692d-4
     @         *gammabb*t13*t33*t34-3.250686379069345d-3*gammabb*t11*t12
     1         *t33*t34+6.501372758138692d-4*gammabb*t2*t32*t33-3.250686
     2         379069346d-4*gammabb*t15*t23*t33-1.6253431895346726d-3*ga
     3         mmabb*t15*t21*t22*t33-1.300274551627738d-3*gammabb*t15*t1
     4         9*t20*t33-9.752059137208037d-4*gammabb*t15*t17*t18*t33-6.
     5         501372758138692d-4*gammabb*t14*t15*t16*t33)*wght-4.135669
     6         939329334d-1*t30*t40*wght-2.4814019635976003d+0*t25*t39*w
     7         ght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363491002d-1*t1*(1.595128
     1         677741414d-5*gammabb*t12*t46*t51*t9-3.190257355482828d-6*
     2         gammabb*t10*t32*t51*t9+5.778998007234393d-4*t10*t34*t38*t
     3         9-1.4447495018085982d-3*t12*t32*t38*t9-2.3926930166121207
     4         d-6*gammabb*t32*t51*t7*t8+4.334248505425794d-4*t34*t38*t7
     5         *t8+4.785386033224242d-6*gammabb*t46*t51*t6*t8-4.33424850
     6         5425794d-4*t32*t38*t6*t8-7.97564338870707d-7*gammabb*t44*
     7         t51*t8+9.570772066448484d-6*gammabb*t10*t46*t51*t7-7.9756
     8         4338870707d-6*gammabb*t12*t44*t51*t7-8.668497010851588d-4
     9         *t10*t32*t38*t7-3.190257355482828d-6*gammabb*t10*t44*t51*
     :         t6-1.595128677741414d-6*gammabb*t2*t32*t51*t6+2.889499003
     ;         6171964d-4*t2*t34*t38*t6-7.97564338870707d-6*gammabb*t17*
     <         t22*t49*t51-3.190257355482828d-6*gammabb*t16*t20*t49*t51-
     =         7.97564338870707d-7*gammabb*t18*t49*t51+1.595128677741414
     >         2d-6*gammabb*t2*t46*t51+7.97564338870707d-6*gammabb*t19*t
     ?         22*t36*t51+4.785386033224242d-6*gammabb*t17*t20*t36*t51+2
     @         .392693016612121d-6*gammabb*t16*t18*t36*t51+7.97564338870
     1         7071d-7*gammabb*t14*t36*t51-7.97564338870707d-7*gammabb*t
     2         13*t32*t51-3.987821694353535d-6*gammabb*t11*t12*t32*t51-1
     3         .329273898117845d-7*gammabb*t15*t23*t51-6.646369490589226
     4         d-7*gammabb*t15*t21*t22*t51-5.31709559247138d-7*gammabb*t
     5         15*t19*t20*t51-3.987821694353535d-7*gammabb*t15*t17*t18*t
     6         51-2.65854779623569d-7*gammabb*t14*t15*t16*t51-1.44474950
     7         18085982d-3*t19*t22*t36*t38-8.668497010851588d-4*t17*t20*
     8         t36*t38-4.334248505425794d-4*t16*t18*t36*t38-1.4447495018
     9         085982d-4*t14*t36*t38+1.4447495018085982d-4*t13*t34*t38+7
     :         .22374750904299d-4*t11*t12*t34*t38-1.4447495018085982d-4*
     ;         t2*t32*t38+7.22374750904299d-5*t15*t23*t38+3.611873754521
     <         495d-4*t15*t21*t22*t38+2.8894990036171964d-4*t15*t19*t20*
     =         t38+2.167124252712897d-4*t15*t17*t18*t38+1.44474950180859
     >         82d-4*t14*t15*t16*t38)*wght-1.2407009817988002d+0*t25*(-1
     ?         .0835621263564485d-4*t10*t34*t41*t9+2.708905315891122d-4*
     @         t12*t32*t41*t9-8.126715947673364d-5*t34*t41*t7*t8+8.12671
     1         5947673364d-5*t32*t41*t6*t8+1.625343189534673d-4*t10*t32*
     2         t41*t7-5.417810631782243d-5*t2*t34*t41*t6+2.7089053158911
     3         22d-4*t19*t22*t36*t41+1.625343189534673d-4*t17*t20*t36*t4
     4         1+8.126715947673364d-5*t16*t18*t36*t41+2.7089053158911214
     5         d-5*t14*t36*t41-2.7089053158911214d-5*t13*t34*t41-1.35445
     6         2657945561d-4*t11*t12*t34*t41+2.7089053158911214d-5*t2*t3
     7         2*t41-1.3544526579455607d-5*t15*t23*t41-6.772263289727805
     8         d-5*t15*t21*t22*t41-5.417810631782243d-5*t15*t19*t20*t41-
     9         4.063357973836682d-5*t15*t17*t18*t41-2.7089053158911214d-
     :         5*t14*t15*t16*t41)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491002d-1*t1*(t26*t3*(-1.9939108471767675d-6*t12*t34*t
     2         41*t9-5.981732541530301d-7*t34*t41*t6*t8+2.99086627076515
     3         07d-7*t32*t41*t8-1.1963465083060604d-6*t10*t34*t41*t7+2.9
     4         90866270765151d-6*t12*t32*t41*t7+1.1963465083060604d-6*t1
     5         0*t32*t41*t6-1.9939108471767675d-7*t2*t34*t41)+t3*(2.9908
     6         66270765151d-6*t17*t22*t41*t49+1.1963465083060604d-6*t16*
     7         t20*t41*t49+2.9908662707651507d-7*t18*t41*t49-2.990866270
     8         765151d-6*t19*t22*t36*t41-1.7945197624590903d-6*t17*t20*t
     9         36*t41-8.972598812295453d-7*t16*t18*t36*t41-2.99086627076
     :         5151d-7*t14*t36*t41+4.9847771179419187d-8*t15*t23*t41+2.4
     ;         923885589709593d-7*t15*t21*t22*t41+1.9939108471767675d-7*
     <         t15*t19*t20*t41+1.4954331353825753d-7*t15*t17*t18*t41+9.9
     =         69554235883837d-8*t14*t15*t16*t41)-1.4721155704335503d-2*
     >         t34*t41*t42+8.126715947673364d-5*t29*t32/rhob**8)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
