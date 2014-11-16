C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_pbeop.F 
C> The nwxcm_c_pbeop functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_pbeop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}+1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}\,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_pbeop(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = param(1)
            t2 = 4.492667444378077d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t3 = 1.804d+0-8.04d-1/t2
            t4 = rhoa**3.333333333333333d-1
            t5 = t1*t3*t4
            t6 = t5**2.0d+0
            t7 = 1/t6
            t8 = 8.658781460080349d-1*t6+1.0500052408963245d+0*t1*t3*t4+
     1         3.183d-1
            t9 = 1/t8
            t10 = 1.415701855281521d+0*t1*t3*t4+5.764d-1
            t11 = rhoa**2
            t12 = 1/t2**2
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**6.666666666666666d-1
            t15 = 3.333333333333333d-1*t1*t14*t3-9.632279000746596d-3*ga
     1         mmaaa*t1*t12*t13
            t16 = 1/t5**3.0d+0
            t17 = 1/t8**2
            t18 = 1/t4
            t19 = -5.905753878944037d-3*t1*t12*t18*t7*t9*wght+8.34321698
     1         0201866d-3*t1*t10*t12*t16*t18*t9*wght+1.1548969154728159d
     2         +0*t10*t11*t17*(6.255284912248943d-3*t1**2*t12*t3/t11+3.7
     3         92728787209827d-3*t1*t12/rhoa**2.3333333333333334d+0)*t7*
     4         wght
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t10*t11*t7*t9*wght
            Amat(iq,D1_RA) = t7*t9*(-1.1548969154728159d+0*t11*(4.719006
     1         184271737d-1*t1*t14*t3-1.363643525194619d-2*gammaaa*t1*t1
     2         2*t13)*wght-2.3097938309456317d+0*rhoa*t10*wght)+2.309793
     3         8309456317d+0*t10*t11*t15*t16*t9*wght+1.1548969154728159d
     4         +0*t10*t11*t17*(1.7317562920160698d+0*t1*t15*t3*t4+3.5000
     5         17469654415d-1*t1*t14*t3-1.0113943432559537d-2*gammaaa*t1
     6         *t12*t13)*t7*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t19+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t19+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = 4.492667444378077d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t3 = 1.804d+0-8.04d-1/t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 4.492667444378077d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t6 = 1.804d+0-8.04d-1/t5
            t7 = rhob**3.333333333333333d-1
            t8 = 1.8610514726982d+0*t6*t7+1.8610514726982d+0*t3*t4
            t9 = 1/t8
            t10 = t1*t3*t4*t6*t7*t9
            t11 = t10**2.0d+0
            t12 = 1/t11
            t13 = 3.908227599821866d+0*t1*t3*t4*t6*t7*t9+1.1995919419748
     1         984d+1*t11+3.183d-1
            t14 = 1/t13
            t15 = 5.269388045346497d+0*t1*t3*t4*t6*t7*t9+5.764d-1
            t16 = 1/t10**3.0d+0
            t17 = 1/t2**2
            t18 = 1/rhoa**3.3333333333333337d+0
            t19 = 1/rhoa**6.666666666666666d-1
            t20 = 6.203504908994d-1*t19*t3-1.7926167019779396d-2*gammaaa
     1         *t17*t18
            t21 = 1/t8**2
            t22 = 3.333333333333333d-1*t1*t19*t3*t6*t7*t9-9.632279000746
     1         596d-3*gammaaa*t1*t17*t18*t6*t7*t9-t1*t20*t21*t3*t4*t6*t7
            t23 = 1/t13**2
            t24 = 1/t5**2
            t25 = 1/rhob**3.3333333333333337d+0
            t26 = 1/rhob**6.666666666666666d-1
            t27 = 6.203504908994d-1*t26*t6-1.7926167019779396d-2*gammabb
     1         *t24*t25
            t28 = 3.333333333333333d-1*t1*t26*t3*t4*t6*t9-9.632279000746
     1         596d-3*gammabb*t1*t24*t25*t3*t4*t9-t1*t21*t27*t3*t4*t6*t7
            t29 = 1/rhoa**2
            t30 = 1/rhoa**2.3333333333333334d+0
            t31 = 3.612104625279974d-3*t1*t17*t30*t6*t7*t9-6.72231263241
     1         7276d-3*t1*t17*t21*t29*t3*t6*t7
            t32 = 1/rhob**2.3333333333333334d+0
            t33 = 1/rhob**2
            t34 = 3.612104625279974d-3*t1*t24*t3*t32*t4*t9-6.72231263241
     1         7276d-3*t1*t21*t24*t3*t33*t4*t6
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t12*t14*t15
     1         *wght
            Amat(iq,D1_RA) = t12*t14*(-8.336168033553906d-2*rhoa*rhob*(1
     1         .7564626817821655d+0*t1*t19*t3*t6*t7*t9-5.075621581597621
     2         d-2*gammaaa*t1*t17*t18*t6*t7*t9-5.269388045346497d+0*t1*t
     3         20*t21*t3*t4*t6*t7)*wght-8.336168033553906d-2*rhob*t15*wg
     4         ht)+8.336168033553906d-2*rhoa*rhob*t12*t15*t23*(2.3991838
     5         839497967d+1*t1*t22*t3*t4*t6*t7*t9+1.302742533273955d+0*t
     6         1*t19*t3*t6*t7*t9-3.7645138639902426d-2*gammaaa*t1*t17*t1
     7         8*t6*t7*t9-3.908227599821866d+0*t1*t20*t21*t3*t4*t6*t7)*w
     8         ght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t22*wght+
     9         Amat(iq,D1_RA)
            Amat(iq,D1_RB) = t12*t14*(-8.336168033553906d-2*rhoa*rhob*(1
     1         .7564626817821655d+0*t1*t26*t3*t4*t6*t9-5.075621581597621
     2         d-2*gammabb*t1*t24*t25*t3*t4*t9-5.269388045346497d+0*t1*t
     3         21*t27*t3*t4*t6*t7)*wght-8.336168033553906d-2*rhoa*t15*wg
     4         ht)+8.336168033553906d-2*rhoa*rhob*t12*t15*t23*(2.3991838
     5         839497967d+1*t1*t28*t3*t4*t6*t7*t9+1.302742533273955d+0*t
     6         1*t26*t3*t4*t6*t9-3.7645138639902426d-2*gammabb*t1*t24*t2
     7         5*t3*t4*t9-3.908227599821866d+0*t1*t21*t27*t3*t4*t6*t7)*w
     8         ght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t28*wght+
     9         Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23
     1         *(2.3991838839497967d+1*t1*t3*t31*t4*t6*t7*t9+1.411692698
     2         9963413d-2*t1*t17*t30*t6*t7*t9-2.6272327764644376d-2*t1*t
     3         17*t21*t29*t3*t6*t7)*wght-8.336168033553906d-2*rhoa*rhob*
     4         t12*t14*(1.9033580930991081d-2*t1*t17*t30*t6*t7*t9-3.5422
     5         473822341327d-2*t1*t17*t21*t29*t3*t6*t7)*wght+1.667233606
     6         7107812d-1*rhoa*rhob*t14*t15*t16*t31*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23
     1         *(2.3991838839497967d+1*t1*t3*t34*t4*t6*t7*t9+1.411692698
     2         9963413d-2*t1*t24*t3*t32*t4*t9-2.6272327764644376d-2*t1*t
     3         21*t24*t3*t33*t4*t6)*wght-8.336168033553906d-2*rhoa*rhob*
     4         t12*t14*(1.9033580930991081d-2*t1*t24*t3*t32*t4*t9-3.5422
     5         473822341327d-2*t1*t21*t24*t3*t33*t4*t6)*wght+1.667233606
     6         7107812d-1*rhoa*rhob*t14*t15*t16*t34*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_pbeop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}+1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}\,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_pbeop_d2(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
C> \brief Evaluate the nwxcm_c_pbeop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}+1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,\left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}\,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111, 5656 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
      subroutine nwxcm_c_pbeop_d3(param,tol_rho,ipol,nq,wght, 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
