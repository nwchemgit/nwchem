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
C>   {\it t_4} &=& {{0.003680288926083876\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
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
C>   {\it t_{16}} &=& {{0.003680288926083876\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{\it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{\it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}+{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}+{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}+{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)-0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{\it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}+{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}+{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{\it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{\it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{\it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}+{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}+{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}+{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
      subroutine nwxcm_x_sogga(param,tol_rho,ipol,nq,wght, 
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
            t27 = -1.8610514726982003d+0*t1*(t25*t3*(1.4721155704335503d
     1         -2*t10*t9+1.1040866778251626d-2*t7*t8+7.360577852167751d-
     2         3*t2*t6+3.6802889260838756d-3*t13+1.8401444630419378d-2*t
     3         11*t12)+(3.6802889260838756d-3*t15*t23+1.8401444630419378
     4         d-2*t15*t21*t22+1.4721155704335503d-2*t15*t19*t20+1.10408
     5         66778251626d-2*t15*t17*t18+7.360577852167751d-3*t14*t15*t
     6         16)*t3)*wght
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t24*wght
            Amat(iq,D1_RA) = -1.8610514726982003d+0*t1*(-3.9256415211561
     1         34d-2*gammaaa*t10*t25*t26*t9-2.9442311408671007d-2*gammaa
     2         a*t25*t26*t7*t8-1.962820760578067d-2*gammaaa*t2*t25*t26*t
     3         6-9.814103802890335d-3*gammaaa*t13*t25*t26-4.907051901445
     4         1675d-2*gammaaa*t11*t12*t25*t26-9.814103802890335d-3*gamm
     5         aaa*t15*t23*t26-4.9070519014451675d-2*gammaaa*t15*t21*t22
     6         *t26-3.925641521156134d-2*gammaaa*t15*t19*t20*t26-2.94423
     7         11408671007d-2*gammaaa*t15*t17*t18*t26-1.962820760578067d
     8         -2*gammaaa*t14*t15*t16*t26)*wght-2.4814019635976003d+0*rh
     9         oa**3.333333333333333d-1*t24*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t27+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t27+Cmat(iq,D1_GAB)
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
C>   {\it t_4} &=& {{0.003680288926083876\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
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
C>   {\it t_{16}} &=& {{0.003680288926083876\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{\it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{\it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}+{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}+{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}+{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)-0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{\it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}+{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}+{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{\it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{\it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{\it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}+{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}+{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}+{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
      subroutine nwxcm_x_sogga_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_sogga functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(7\right)\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_4} &=& {{0.003680288926083876\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
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
C>   {\it t_{16}} &=& {{0.003680288926083876\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{17}} &=& 1.0-{{1.0}\over{{\it t_{16}}+1.0}}\\\\ 
C>   {\it t_{18}} &=& 1.0-{{1}\over{e^{{\it t_{16}}}}}\\\\ 
C>   {\it t_{19}} &=& {{0.003680288926083876\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{20}} &=& 1.0-{{1.0}\over{{\it t_{19}}+1.0}}\\\\ 
C>   {\it t_{21}} &=& 1.0-{{1}\over{e^{{\it t_{19}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{18}}+{\it t_{14}}\,{\it t_{18}}^{5.0}+{\it t_{13}}\,{\it t_{18}}^{4.0}+{\it t_{12}}\,{\it t_{18}}^{3.0}+{\it t_{10}}\,{\it t_{18}}^{2.0}+{\it t_9}\,{\it t_{17}}^{5.0}+{\it t_8}\,{\it t_{17}}^{4.0}+{\it t_7}\,{\it t_{17}}^{3.0}+{\it t_6}\,{\it t_{17}}^{2.0}+{\it t_3}\,{\it t_{17}}+{\it t_2}+{\it t_1}\right)-0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{11}}+{\it t_{14}}\,{\it t_{11}}^{5.0}+{\it t_{13}}\,{\it t_{11}}^{4.0}+{\it t_{12}}\,{\it t_{11}}^{3.0}+{\it t_{10}}\,{\it t_{11}}^{2.0}+{\it t_9}\,{\it t_5}^{5.0}+{\it t_8}\,{\it t_5}^{4.0}+{\it t_7}\,{\it t_5}^{3.0}+{\it t_6}\,{\it t_5}^{2.0}+{\it t_3}\,{\it t_5}+{\it t_2}+{\it t_1}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,\rho_s^{{{4}\over{3}}}\,\left({\it t_{15}}\,{\it t_{21}}+{\it t_{14}}\,{\it t_{21}}^{5.0}+{\it t_{13}}\,{\it t_{21}}^{4.0}+{\it t_{12}}\,{\it t_{21}}^{3.0}+{\it t_{10}}\,{\it t_{21}}^{2.0}+{\it t_9}\,{\it t_{20}}^{5.0}+{\it t_8}\,{\it t_{20}}^{4.0}+{\it t_7}\,{\it t_{20}}^{3.0}+{\it t_6}\,{\it t_{20}}^{2.0}+{\it t_3}\,{\it t_{20}}+{\it t_2}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Chem.Phys. 128, 184109 (2008)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.2912068 "> 
C> 10.1063/1.2912068 </a> 
C> 
      subroutine nwxcm_x_sogga_d3(param,tol_rho,ipol,nq,wght, 
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
