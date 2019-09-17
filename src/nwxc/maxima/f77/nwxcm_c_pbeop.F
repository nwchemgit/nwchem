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
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}
C>    +1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}
C>    \,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}
C>    \,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{
C>    \it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}
C>    +0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_pbeop(param,tol_rho,ipol,nq,wght, 
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
            t2 = 4.492667444378077d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t3 = 1.804d+0-8.04d-1/t2
            t4 = rhoa**3.333333333333333d-1
            t5 = t1*t3*t4
            t6 = t5**2.0d+0
            t7 = 1/t6
            t8 = 1.0500052408963245d+0*t1*t3*t4
            t9 = 1.415701855281521d+0*t1*t3*t4+5.764d-1
            t10 = rhoa**2
            t11 = 1/t2**2
            t12 = 1/rhoa**3.3333333333333337d+0
            t13 = 1/rhoa**6.666666666666666d-1
            t14 = 6.203504908994d-1*t13*t3-1.7926167019779396d-2*gammaaa
     1         *t11*t12
            t15 = 8.955510855647054d-2*t1*t13*t3-7.218105721705101d-2*t1
     1         *t14-2.587859374674219d-3*gammaaa*t1*t11*t12
            t16 = 1/t5**3.0d+0
            t17 = t8+8.658781460080348d-1*t6+3.183d-1
            t18 = 1/t17
            t19 = 1/t17**2
            t20 = 1/t4
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t10*t7*t9*wght/(t8+8
     1         .658781460080349d-1*t6+3.183d-1)
            Amat(iq,D1_RA) = 1.385405033612855d+1*t18*t7*(-8.33616803355
     1         3906d-2*rhoa*t9*wght-8.336168033553906d-2*t10*(4.71900618
     2         4271736d-1*t1*t13*t3-3.8035d-1*t1*t14-1.363643525194619d-
     3         2*gammaaa*t1*t11*t12)*wght)+1.1548969154728159d+0*t10*t19
     4         *(6.44577519522176d+0*t1*t15*t3*t4+3.5000174696544145d-1*
     5         t1*t13*t3-2.821d-1*t1*t14-1.0113943432559537d-2*gammaaa*t
     6         1*t11*t12)*t7*t9*wght+8.597290421421169d+0*t10*t15*t16*t1
     7         8*t9*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.1548969154728159d+0*t10*t19*(3.127642456
     1         1244703d-3*t1**2*t11*t3/t10+1.8963643936049138d-3*t1*t11/
     2         rhoa**2.3333333333333334d+0)*t7*t9*wght+4.171608490100933
     3         d-3*t1*t11*t16*t18*t20*t9*wght-2.9528769394720183d-3*t1*t
     4         11*t18*t20*t7*wght+Cmat(iq,D1_GAA)
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
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}
C>    +1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}
C>    \,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}
C>    \,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{
C>    \it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}
C>    +0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_pbeop_d2(param,tol_rho,ipol,nq,wght, 
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
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.492667444378077d-3*gammaaa*t2+1.0d+0
            t4 = 1.804d+0-8.04d-1/t3
            t5 = rhoa**3.333333333333333d-1
            t6 = t1*t4*t5
            t7 = t6**2.0d+0
            t8 = 1/t7
            t9 = 1.0500052408963245d+0*t1*t4*t5
            t10 = 1.415701855281521d+0*t1*t4*t5+5.764d-1
            t11 = rhoa**2
            t12 = 1/t3**2
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**6.666666666666666d-1
            t15 = 6.203504908994d-1*t14*t4-1.7926167019779396d-2*gammaaa
     1         *t12*t13
            t16 = 8.955510855647054d-2*t1*t14*t4-7.218105721705101d-2*t1
     1         *t15-2.587859374674219d-3*gammaaa*t1*t12*t13
            t17 = 1/t6**3.0d+0
            t18 = t9+8.658781460080348d-1*t7+3.183d-1
            t19 = 1/t18
            t20 = 1/t18**2
            t21 = 6.44577519522176d+0*t1*t16*t4*t5+3.5000174696544145d-1
     1         *t1*t14*t4-2.821d-1*t1*t15-1.0113943432559537d-2*gammaaa*
     2         t1*t12*t13
            t22 = 4.719006184271736d-1*t1*t14*t4-3.8035d-1*t1*t15-1.3636
     1         43525194619d-2*gammaaa*t1*t12*t13
            t23 = -8.336168033553906d-2*t11*t22*wght-8.336168033553906d-
     1         2*rhoa*t10*wght
            t24 = 1/t5
            t25 = 1/rhoa**2.3333333333333334d+0
            t26 = 1.8963643936049138d-3*t1*t12*t25
            t27 = t1**2
            t28 = 1/t11
            t29 = 3.1276424561244703d-3*t12*t27*t28*t4+t26
            t30 = 1/t18**3
            t31 = -2.3097938309456317d+0*t10*t11*t21**2*t30*t8*wght
            t32 = 1/t6**4.0d+0
            t33 = -9.599999999999997d+1*t10*t11*t16**2*t19*t32*wght
            t34 = -1.7194580842842339d+1*t10*t11*t16*t17*t20*t21*wght
            t35 = -1.6672336067107812d-1*rhoa*t22*wght
            t36 = gammaaa**2
            t37 = 1/t3**3
            t38 = 1/rhoa**7
            t39 = 1/rhoa**4.333333333333333d+0
            t40 = 1/t4
            t41 = 1/rhoa**3.6666666666666664d+0
            t42 = 7.327274635867936d-3*gammaaa*t1*t12*t15*t40*t41
            t43 = 1/rhoa**1.6666666666666669d+0
            t44 = -4.135669939329333d-1*t4*t43+5.377850105933819d-2*gamm
     1         aaa*t12*t39-4.295269705186499d-4*t36*t37*t38
            t45 = 1/rhoa
            t46 = -2.5356666666666666d-1*t1*t15*t45
            t47 = t15**2
            t48 = 2.0437371323672143d-1*t1*t24*t40*t47
            t49 = 1.6672336067107812d-1*rhoa*t10*t16*wght
            t50 = 1.6672336067107812d-1*t11*t16*t22*wght
            t51 = 1.390536163366978d-3*gammaaa*t1*t12*t15*t40*t41
            t52 = -4.812070481136733d-2*t1*t15*t45
            t53 = 3.8785094488762883d-2*t1*t24*t40*t47
            t54 = t53+t52+t51-7.218105721705101d-2*t1*t44-5.970340570431
     1         37d-2*t1*t4*t43+7.76357812402266d-3*gammaaa*t1*t12*t39-6.
     2         20074216705465d-5*t1*t36*t37*t38
            t55 = -2.0d+0*t16*t23
            t56 = 8.336168033553906d-2*rhoa*t10*t21*wght
            t57 = 8.336168033553906d-2*t11*t21*t22*wght
            t58 = 5.434531812221229d-3*gammaaa*t1*t12*t15*t40*t41
            t59 = -6.208750505646785d-2*gammaaa*t1*t12*t13*t16
            t60 = -1.8806666666666666d-1*t1*t15*t45
            t61 = 2.1485917317405864d+0*t1*t14*t16*t4
            t62 = -1.7317562920160695d+0*t1*t15*t16
            t63 = 1.5158097674268203d-1*t1*t24*t40*t47
            t64 = -t21*t23
            t65 = 1/t3**4
            t66 = 2.4926983511559697d-5*t1*t36*t38*t40*t65+t53+t52+t51+2
     1         .985170285215685d-2*t1*t4*t43-1.7252395831161463d-3*gamma
     2         aa*t1*t12*t39
            t67 = -4.6581468744135934d-2*t1*t10*t12*t16*t19*t24*t32*wght
            t68 = -2.1314177932294023d-4*t1*t12*wght/rhoa**1.33333333333
     1         33333d+0
            t69 = 1/rhoa**6
            t70 = 1.610726139444937d-4*gammaaa*t37*t69-1.568539614230697
     1         3d-2*t12*t13
            t71 = -8.08981147293451d-5*t1*t10*t12*t21*t24*wght
            t72 = 4.2628355864588047d-4*t1*t12*t16*t24*wght
            t73 = -7.218105721705101d-2*t1*t70+4.673809408417444d-6*gamm
     1         aaa*t1*t40*t65*t69+2.3252783126454943d-5*gammaaa*t1*t37*t
     2         69-2.4261181637570808d-3*t1*t12*t13
            t74 = -9.704472655028322d-4*t1*t12*t23*t25
            t75 = 2.1314177932294023d-4*t1*t12*t21*t24*wght
            t76 = -2.168404344971009d-19*t1*t12*t15*t2*t40
            t77 = 1.1641407198087722d-2*t1*t12*t16*t25
            t78 = 3.127642456124471d-3*t12*t27*t28*t4+t26
            t79 = 1.617412109171387d-4*t1*t12*t13-4.673809408417444d-6*g
     1         ammaaa*t1*t40*t65*t69
            t80 = 1/rhoa**5
            t81 = 4.1368571381671093d-7*t2*t27*t65*wght
            t82 = -1.7526785281565413d-6*t1*t40*t65*t80-4.35989683621030
     1         2d-6*t1*t37*t80
            t83 = 1/rhoa**4.666666666666667d+0
            t84 = 1/rhoa**3
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t10*t11*t8*wght/(t9+
     1         8.658781460080349d-1*t7+3.183d-1)
            Amat(iq,D1_RA) = 1.1548969154728159d+0*t10*t11*t20*t21*t8*wg
     1         ht+8.597290421421169d+0*t10*t11*t16*t17*t19*wght+1.385405
     2         033612855d+1*t19*t23*t8+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.1548969154728159d+0*t10*t11*t20*t29*t8*w
     1         ght-2.9528769394720183d-3*t1*t12*t19*t24*t8*wght+4.171608
     2         490100933d-3*t1*t10*t12*t17*t19*t24*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.385405033612855d+1*t20*t8*(8.33616803
     1         3553906d-2*t10*t11*(t63+t62+t61+t60+t59+t58+6.44577519522
     2         176d+0*t1*t4*t5*t54-2.821d-1*t1*t44-2.3333449797696096d-1
     3         *t1*t4*t43+3.034183029767861d-2*gammaaa*t1*t12*t39-2.4233
     4         911676662234d-4*t1*t36*t37*t38)*wght+t64+t57+t56)+5.15662
     5         0156177406d+1*t17*t19*(1.6672336067107812d-1*t10*t11*t54*
     6         wght+t55+t50+t49)+1.385405033612855d+1*t19*t8*(t35-8.3361
     7         68033553906d-2*t11*(t48+t46-3.8035d-1*t1*t44-3.1460041228
     8         478236d-1*t1*t4*t43+t42+4.090930575583857d-2*gammaaa*t1*t
     9         12*t39-3.2674116647353707d-4*t1*t36*t37*t38)*wght)+t34+t3
     :         3+t31+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.385405033612855d+1*t20*t8*(8.33616803
     1         3553906d-2*t10*t11*(6.44577519522176d+0*t1*t4*t5*t66+9.74
     2         2032494018218d-5*t1*t36*t38*t40*t65+t63+t62+t61+t60+t59+t
     3         58+1.1666724898848048d-1*t1*t4*t43-6.742628955039691d-3*g
     4         ammaaa*t1*t12*t39)*wght+t64+t57+t56)+5.156620156177406d+1
     5         *t17*t19*(1.6672336067107812d-1*t10*t11*t66*wght+t55+t50+
     6         t49)+1.385405033612855d+1*t19*t8*(-8.336168033553906d-2*t
     7         11*(1.313499489223619d-4*t1*t36*t38*t40*t65+t48+t46+1.573
     8         0020614239118d-1*t1*t4*t43+t42-9.090956834630795d-3*gamma
     9         aa*t1*t12*t39)*wght-8.336168033553906d-2*t10*wght+t35)+t3
     :         4+t33+t31+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.385405033612855d+1*t20*t8*(8.3361680
     1         33553906d-2*t10*t11*(t77+t76+6.44577519522176d+0*t1*t4*t5
     2         *t73-2.821d-1*t1*t70+1.8266310926284157d-5*gammaaa*t1*t40
     3         *t65*t69+9.087716878748338d-5*gammaaa*t1*t37*t69-9.481821
     4         968024564d-3*t1*t12*t13)*wght+t75-t23*t29)+5.156620156177
     5         406d+1*t17*t19*(1.6672336067107812d-1*t10*t11*t73*wght+t7
     6         4+t72)+1.385405033612855d+1*t19*t8*(t68-8.336168033553906
     7         d-2*t11*(-3.8035d-1*t1*t70+2.4628115422942856d-5*gammaaa*
     8         t1*t40*t65*t69+1.2252793742757642d-4*gammaaa*t1*t37*t69-1
     9         .2784158048699554d-2*t1*t12*t13)*wght)+5.156620156177406d
     :         +1*t17*t20*(t71-1.6672336067107812d-1*t10*t11*t16*t29*wgh
     ;         t)-2.3097938309456317d+0*t10*t11*t21*t29*t30*t8*wght+t67+
     <         Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 1.385405033612855d+1*t20*t8*(8.3361680
     1         33553906d-2*t10*t11*(6.44577519522176d+0*t1*t4*t5*t79+t77
     2         +t76-1.8266310926284166d-5*gammaaa*t1*t40*t65*t69+6.32121
     3         4645349711d-4*t1*t12*t13)*wght-t23*t78+t75)+5.15662015617
     4         7406d+1*t17*t19*(1.6672336067107812d-1*t10*t11*t79*wght+t
     5         74+t72)+5.156620156177406d+1*t17*t20*(t71-1.6672336067107
     6         812d-1*t10*t11*t16*t78*wght)+1.385405033612855d+1*t19*t8*
     7         (t68-8.336168033553906d-2*t11*(8.522772032466369d-4*t1*t1
     8         2*t13-2.4628115422942862d-5*gammaaa*t1*t40*t65*t69)*wght)
     9         -2.3097938309456317d+0*t10*t11*t21*t30*t78*t8*wght+t67+Cm
     :         at2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.385405033612855d+1*t20*t8*(8.336168
     1         033553906d-2*t10*t11*(5.648685890994609d-6*t27*t65*t83+6.
     2         44577519522176d+0*t1*t4*t5*t82-6.84986659735656d-6*t1*t40
     3         *t65*t80-1.7039469147653133d-5*t1*t37*t80)*wght+4.2628355
     4         864588047d-4*t1*t12*t24*t29*wght)+5.156620156177406d+1*t1
     5         7*t19*(1.6672336067107812d-1*t10*t11*t82*wght+t81)-1.1548
     6         969154728159d+0*t11*t19*t8*(-9.235543283603571d-6*t1*t40*
     7         t65*t80-2.2973988267670578d-5*t1*t37*t80)*wght-2.30979383
     8         09456317d+0*t10*t11*t29**2*t30*t8*wght-2.260242948292618d
     9         -5*t10*t19*t2*t27*t32*t65*wght-8.343216980201866d-3*t1*t1
     :         0*t12*t17*t20*t24*t29*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = 5.156620156177406d+1*t17*t19*(2.92212
     1         4543902974d-7*t1*t10*t40*t65*t84*wght+t81)+1.385405033612
     2         855d+1*t20*t8*(8.336168033553906d-2*t10*t11*(1.6946057672
     3         983827d-5*t27*t65*t83+6.849866597356562d-6*t1*t40*t65*t80
     4         )*wght+2.1314177932294023d-4*t1*t12*t24*t78*wght+2.131417
     5         7932294023d-4*t1*t12*t24*t29*wght)+5.156620156177406d+1*t
     6         17*t20*(-8.08981147293451d-5*t1*t10*t12*t24*t78*wght-8.08
     7         981147293451d-5*t1*t10*t12*t24*t29*wght)-1.06661004509494
     8         46d-5*t1*t19*t40*t65*t8*t84*wght-2.3097938309456317d+0*t1
     9         0*t11*t29*t30*t78*t8*wght-2.2602429482926184d-5*t10*t19*t
     :         2*t27*t32*t65*wght+Cmat2(iq,D2_GAA_GBB)
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
            t24 = 2.3991838839497967d+1*t1*t22*t3*t4*t6*t7*t9+1.30274253
     1         3273955d+0*t1*t19*t3*t6*t7*t9-3.7645138639902426d-2*gamma
     2         aa*t1*t17*t18*t6*t7*t9-3.908227599821866d+0*t1*t20*t21*t3
     3         *t4*t6*t7
            t25 = 1.7564626817821655d+0*t1*t19*t3*t6*t7*t9-5.07562158159
     1         7621d-2*gammaaa*t1*t17*t18*t6*t7*t9-5.269388045346497d+0*
     2         t1*t20*t21*t3*t4*t6*t7
            t26 = -8.336168033553906d-2*rhoa*rhob*t25*wght-8.33616803355
     1         3906d-2*rhob*t15*wght
            t27 = 1/t5**2
            t28 = 1/rhob**3.3333333333333337d+0
            t29 = 1/rhob**6.666666666666666d-1
            t30 = 6.203504908994d-1*t29*t6-1.7926167019779396d-2*gammabb
     1         *t27*t28
            t31 = 3.333333333333333d-1*t1*t29*t3*t4*t6*t9-9.632279000746
     1         596d-3*gammabb*t1*t27*t28*t3*t4*t9-t1*t21*t3*t30*t4*t6*t7
            t32 = 2.3991838839497967d+1*t1*t3*t31*t4*t6*t7*t9+1.30274253
     1         3273955d+0*t1*t29*t3*t4*t6*t9-3.7645138639902426d-2*gamma
     2         bb*t1*t27*t28*t3*t4*t9-3.908227599821866d+0*t1*t21*t3*t30
     3         *t4*t6*t7
            t33 = 1.7564626817821655d+0*t1*t29*t3*t4*t6*t9-5.07562158159
     1         7621d-2*gammabb*t1*t27*t28*t3*t4*t9-5.269388045346497d+0*
     2         t1*t21*t3*t30*t4*t6*t7
            t34 = -8.336168033553906d-2*rhoa*rhob*t33*wght-8.33616803355
     1         3906d-2*rhoa*t15*wght
            t35 = 1/rhoa**2
            t36 = 1/rhoa**2.3333333333333334d+0
            t37 = 1.9033580930991081d-2*t1*t17*t36*t6*t7*t9-3.5422473822
     1         341327d-2*t1*t17*t21*t3*t35*t6*t7
            t38 = 3.612104625279974d-3*t1*t17*t36*t6*t7*t9-6.72231263241
     1         7276d-3*t1*t17*t21*t3*t35*t6*t7
            t39 = 2.3991838839497967d+1*t1*t3*t38*t4*t6*t7*t9+1.41169269
     1         89963413d-2*t1*t17*t36*t6*t7*t9-2.6272327764644376d-2*t1*
     2         t17*t21*t3*t35*t6*t7
            t40 = 1/rhob**2.3333333333333334d+0
            t41 = 1/rhob**2
            t42 = 1.9033580930991081d-2*t1*t27*t3*t4*t40*t9-3.5422473822
     1         341327d-2*t1*t21*t27*t3*t4*t41*t6
            t43 = 3.612104625279974d-3*t1*t27*t3*t4*t40*t9-6.72231263241
     1         7276d-3*t1*t21*t27*t3*t4*t41*t6
            t44 = 2.3991838839497967d+1*t1*t3*t4*t43*t6*t7*t9+1.41169269
     1         89963413d-2*t1*t27*t3*t4*t40*t9-2.6272327764644376d-2*t1*
     2         t21*t27*t3*t4*t41*t6
            t45 = 1/t13**3
            t46 = 1/t10**4.0d+0
            t47 = t20**2
            t48 = 1/t8**3
            t49 = gammaaa**2
            t50 = 1/t2**3
            t51 = 1/rhoa**7
            t52 = 1/rhoa**4.333333333333333d+0
            t53 = 1/rhoa**1.6666666666666669d+0
            t54 = -4.135669939329333d-1*t3*t53+5.377850105933819d-2*gamm
     1         aaa*t17*t52-4.295269705186499d-4*t49*t50*t51
            t55 = -2.222222222222222d-1*t1*t3*t53*t6*t7*t9+2.88968370022
     1         39795d-2*gammaaa*t1*t17*t52*t6*t7*t9-2.3079800683637772d-
     2         4*t1*t49*t50*t51*t6*t7*t9-t1*t21*t3*t4*t54*t6*t7+2*t1*t3*
     3         t4*t47*t48*t6*t7-6.666666666666666d-1*t1*t19*t20*t21*t3*t
     4         6*t7+1.9264558001493193d-2*gammaaa*t1*t17*t18*t20*t21*t6*
     5         t7
            t56 = 1.111111111111111d-1*t1*t19*t29*t3*t6*t9-3.21075966691
     1         5532d-3*gammaaa*t1*t17*t18*t29*t6*t9-3.210759666915532d-3
     2         *gammabb*t1*t19*t27*t28*t3*t9+9.278079874822385d-5*gammaa
     3         a*gammabb*t1*t17*t18*t27*t28*t9+2*t1*t20*t3*t30*t4*t48*t6
     4         *t7-3.333333333333333d-1*t1*t19*t21*t3*t30*t6*t7+9.632279
     5         000746596d-3*gammaaa*t1*t17*t18*t21*t30*t6*t7-3.333333333
     6         333333d-1*t1*t20*t21*t29*t3*t4*t6+9.632279000746596d-3*ga
     7         mmabb*t1*t20*t21*t27*t28*t3*t4
            t57 = gammabb**2
            t58 = 1/t5**3
            t59 = 1/rhob**7
            t60 = 1/rhob**4.333333333333333d+0
            t61 = 1/rhob**1.6666666666666669d+0
            t62 = t30**2
            t63 = -4.135669939329333d-1*t6*t61+5.377850105933819d-2*gamm
     1         abb*t27*t60-4.295269705186499d-4*t57*t58*t59
            t64 = -2.222222222222222d-1*t1*t3*t4*t6*t61*t9+2.88968370022
     1         39795d-2*gammabb*t1*t27*t3*t4*t60*t9-2.3079800683637772d-
     2         4*t1*t3*t4*t57*t58*t59*t9-t1*t21*t3*t4*t6*t63*t7+2*t1*t3*
     3         t4*t48*t6*t62*t7-6.666666666666666d-1*t1*t21*t29*t3*t30*t
     4         4*t6+1.9264558001493193d-2*gammabb*t1*t21*t27*t28*t3*t30*
     5         t4
            t65 = 1/t2**4
            t66 = 1/rhoa**5.666666666666667d+0
            t67 = 1/rhoa**3
            t68 = 1/rhoa**6
            t69 = 1.610726139444937d-4*gammaaa*t50*t68-1.568539614230697
     1         3d-2*t17*t18
            t70 = 8.654925256364164d-5*gammaaa*t1*t50*t6*t68*t7*t9-8.428
     1         244125653273d-3*t1*t17*t18*t6*t7*t9-t1*t21*t3*t4*t6*t69*t
     2         7-2.240770877472425d-3*t1*t17*t21*t3*t6*t67*t7+6.47511908
     3         056865d-5*gammaaa*t1*t21*t6*t65*t66*t7+1.344462526483455d
     4         -2*t1*t17*t20*t3*t35*t48*t6*t7-3.612104625279974d-3*t1*t1
     5         7*t20*t21*t36*t6*t7
            t71 = 1.2040348750933245d-3*t1*t19*t27*t3*t40*t9-3.479279953
     1         0583945d-5*gammaaa*t1*t17*t18*t27*t40*t9+1.34446252648345
     2         5d-2*t1*t20*t27*t3*t4*t41*t48*t6-2.240770877472425d-3*t1*
     3         t19*t21*t27*t3*t41*t6+6.47511908056865d-5*gammaaa*t1*t17*
     4         t18*t21*t27*t41*t6-3.612104625279974d-3*t1*t20*t21*t27*t3
     5         *t4*t40
            t72 = 1.2040348750933245d-3*t1*t17*t29*t36*t6*t9-3.479279953
     1         0583945d-5*gammabb*t1*t17*t27*t28*t36*t9+1.34446252648345
     2         5d-2*t1*t17*t3*t30*t35*t48*t6*t7-3.612104625279974d-3*t1*
     3         t17*t21*t30*t36*t6*t7-2.240770877472425d-3*t1*t17*t21*t29
     4         *t3*t35*t6+6.47511908056865d-5*gammabb*t1*t17*t21*t27*t28
     5         *t3*t35
            t73 = 1/rhob**6
            t74 = 1/t5**4
            t75 = 1/rhob**5.666666666666667d+0
            t76 = 1/rhob**3
            t77 = 1.610726139444937d-4*gammabb*t58*t73-1.568539614230697
     1         3d-2*t27*t28
            t78 = 8.654925256364164d-5*gammabb*t1*t3*t4*t58*t73*t9-8.428
     1         244125653273d-3*t1*t27*t28*t3*t4*t9-t1*t21*t3*t4*t6*t7*t7
     2         7-2.240770877472425d-3*t1*t21*t27*t3*t4*t6*t76+6.47511908
     3         056865d-5*gammabb*t1*t21*t3*t4*t74*t75+1.344462526483455d
     4         -2*t1*t27*t3*t30*t4*t41*t48*t6-3.612104625279974d-3*t1*t2
     5         1*t27*t3*t30*t4*t40
            t79 = 1/rhoa**4.666666666666667d+0
            t80 = 1/rhoa**5
            t81 = -3.245596971136562d-5*t1*t50*t6*t7*t80*t9-4.8563393104
     1         26487d-5*t1*t21*t6*t65*t7*t79+6.040223022918516d-5*t1*t21
     2         *t3*t50*t6*t7*t79+9.037897425591375d-5*t1*t3*t48*t52*t6*t
     3         65*t7
            t82 = 1.304729982396898d-5*t1*t17*t27*t36*t40*t9+9.037897425
     1         591375d-5*t1*t17*t27*t3*t35*t41*t48*t6-2.4281696552132434
     2         d-5*t1*t17*t21*t27*t36*t41*t6-2.4281696552132434d-5*t1*t1
     3         7*t21*t27*t3*t35*t40
            t83 = 1/rhob**5
            t84 = 1/rhob**4.666666666666667d+0
            t85 = -3.245596971136562d-5*t1*t3*t4*t58*t83*t9-4.8563393104
     1         26487d-5*t1*t21*t3*t4*t74*t84+6.040223022918516d-5*t1*t21
     2         *t3*t4*t58*t6*t84+9.037897425591375d-5*t1*t3*t4*t48*t6*t6
     3         0*t74
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t12*t14*t15
     1         *wght
            Amat(iq,D1_RA) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23*
     1         t24*wght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t22*
     2         wght+t12*t14*t26+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23*
     1         t32*wght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t31*
     2         wght+t12*t14*t34+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23
     1         *t39*wght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t38
     2         *wght-8.336168033553906d-2*rhoa*rhob*t12*t14*t37*wght+Cma
     3         t(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 8.336168033553906d-2*rhoa*rhob*t12*t15*t23
     1         *t44*wght+1.6672336067107812d-1*rhoa*rhob*t14*t15*t16*t43
     2         *wght-8.336168033553906d-2*rhoa*rhob*t12*t14*t42*wght+Cma
     3         t(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = t12*t23*(8.336168033553906d-2*rhoa*rhob
     1         *t15*(2.3991838839497967d+1*t1*t3*t4*t55*t6*t7*t9-8.68495
     2         0221826366d-1*t1*t3*t53*t6*t7*t9+1.1293541591970728d-1*ga
     3         mmaaa*t1*t17*t52*t6*t7*t9-9.02011140301807d-4*t1*t49*t50*
     4         t51*t6*t7*t9+7.997279613165988d+0*t1*t19*t22*t3*t6*t7*t9-
     5         2.3109608544299284d-1*gammaaa*t1*t17*t18*t22*t6*t7*t9-3.9
     6         08227599821866d+0*t1*t21*t3*t4*t54*t6*t7+7.81645519964373
     7         1d+0*t1*t3*t4*t47*t48*t6*t7-2.3991838839497967d+1*t1*t20*
     8         t21*t22*t3*t4*t6*t7-2.60548506654791d+0*t1*t19*t20*t21*t3
     9         *t6*t7+7.529027727980485d-2*gammaaa*t1*t17*t18*t20*t21*t6
     :         *t7)*wght+8.336168033553906d-2*rhoa*rhob*t24*t25*wght+8.3
     ;         36168033553906d-2*rhob*t15*t24*wght-t24*t26)+t12*t14*(-8.
     <         336168033553906d-2*rhoa*rhob*(-1.17097512118811d+0*t1*t3*
     =         t53*t6*t7*t9+1.5226864744792862d-1*gammaaa*t1*t17*t52*t6*
     >         t7*t9-1.2161642581134076d-3*t1*t49*t50*t51*t6*t7*t9-5.269
     ?         388045346497d+0*t1*t21*t3*t4*t54*t6*t7+1.0538776090692993
     @         d+1*t1*t3*t4*t47*t48*t6*t7-3.512925363564331d+0*t1*t19*t2
     1         0*t21*t3*t6*t7+1.0151243163195242d-1*gammaaa*t1*t17*t18*t
     2         20*t21*t6*t7)*wght-1.6672336067107812d-1*rhob*t25*wght)+t
     3         14*t16*(1.6672336067107812d-1*rhoa*rhob*t15*t55*wght+1.66
     4         72336067107812d-1*rhoa*rhob*t22*t25*wght+1.66723360671078
     5         12d-1*rhob*t15*t22*wght-2.0d+0*t22*t26)-5.001700820132344
     6         d-1*rhoa*rhob*t14*t15*t22**2*t46*wght-1.6672336067107812d
     7         -1*rhoa*rhob*t12*t15*t24**2*t45*wght-3.3344672134215625d-
     8         1*rhoa*rhob*t15*t16*t22*t23*t24*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t12*t23*(8.336168033553906d-2*rhoa*rhob
     1         *t15*(2.3991838839497967d+1*t1*t3*t4*t56*t6*t7*t9+7.99727
     2         9613165988d+0*t1*t22*t29*t3*t4*t6*t9+4.342475110913183d-1
     3         *t1*t19*t29*t3*t6*t9-1.2548379546634142d-2*gammaaa*t1*t17
     4         *t18*t29*t6*t9-2.3109608544299284d-1*gammabb*t1*t22*t27*t
     5         28*t3*t4*t9-1.2548379546634142d-2*gammabb*t1*t19*t27*t28*
     6         t3*t9+3.6260847840132643d-4*gammaaa*gammabb*t1*t17*t18*t2
     7         7*t28*t9+7.816455199643731d+0*t1*t20*t3*t30*t4*t48*t6*t7-
     8         2.3991838839497967d+1*t1*t21*t22*t3*t30*t4*t6*t7-1.302742
     9         533273955d+0*t1*t19*t21*t3*t30*t6*t7+3.7645138639902426d-
     :         2*gammaaa*t1*t17*t18*t21*t30*t6*t7-1.302742533273955d+0*t
     ;         1*t20*t21*t29*t3*t4*t6+3.7645138639902426d-2*gammabb*t1*t
     <         20*t21*t27*t28*t3*t4)*wght+8.336168033553906d-2*rhoa*rhob
     =         *t24*t33*wght+8.336168033553906d-2*rhoa*t15*t24*wght-t26*
     >         t32)+t12*t14*(-8.336168033553906d-2*rhoa*rhob*(5.85487560
     ?         594055d-1*t1*t19*t29*t3*t6*t9-1.6918738605325404d-2*gamma
     @         aa*t1*t17*t18*t29*t6*t9-1.6918738605325404d-2*gammabb*t1*
     1         t19*t27*t28*t3*t9+4.8889803176158986d-4*gammaaa*gammabb*t
     2         1*t17*t18*t27*t28*t9+1.0538776090692993d+1*t1*t20*t3*t30*
     3         t4*t48*t6*t7-1.7564626817821655d+0*t1*t19*t21*t3*t30*t6*t
     4         7+5.075621581597621d-2*gammaaa*t1*t17*t18*t21*t30*t6*t7-1
     5         .7564626817821655d+0*t1*t20*t21*t29*t3*t4*t6+5.0756215815
     6         97621d-2*gammabb*t1*t20*t21*t27*t28*t3*t4)*wght-8.3361680
     7         33553906d-2*rhob*t33*wght-8.336168033553906d-2*rhoa*t25*w
     8         ght-8.336168033553906d-2*t15*wght)+t14*t16*(1.66723360671
     9         07812d-1*rhoa*rhob*t15*t56*wght+1.6672336067107812d-1*rho
     :         a*rhob*t22*t33*wght+1.6672336067107812d-1*rhoa*t15*t22*wg
     ;         ht-2.0d+0*t26*t31)+t16*t23*(-1.6672336067107812d-1*rhoa*r
     <         hob*t15*t22*t32*wght-1.6672336067107812d-1*rhoa*rhob*t15*
     =         t24*t31*wght)-5.001700820132344d-1*rhoa*rhob*t14*t15*t22*
     >         t31*t46*wght-1.6672336067107812d-1*rhoa*rhob*t12*t15*t24*
     ?         t32*t45*wght+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = t12*t23*(8.336168033553906d-2*rhoa*rhob
     1         *t15*(2.3991838839497967d+1*t1*t3*t4*t6*t64*t7*t9-8.68495
     2         0221826366d-1*t1*t3*t4*t6*t61*t9+1.1293541591970728d-1*ga
     3         mmabb*t1*t27*t3*t4*t60*t9+7.997279613165988d+0*t1*t29*t3*
     4         t31*t4*t6*t9-9.02011140301807d-4*t1*t3*t4*t57*t58*t59*t9-
     5         2.3109608544299284d-1*gammabb*t1*t27*t28*t3*t31*t4*t9-3.9
     6         08227599821866d+0*t1*t21*t3*t4*t6*t63*t7+7.81645519964373
     7         1d+0*t1*t3*t4*t48*t6*t62*t7-2.3991838839497967d+1*t1*t21*
     8         t3*t30*t31*t4*t6*t7-2.60548506654791d+0*t1*t21*t29*t3*t30
     9         *t4*t6+7.529027727980485d-2*gammabb*t1*t21*t27*t28*t3*t30
     :         *t4)*wght+8.336168033553906d-2*rhoa*rhob*t32*t33*wght+8.3
     ;         36168033553906d-2*rhoa*t15*t32*wght-t32*t34)+t12*t14*(-8.
     <         336168033553906d-2*rhoa*rhob*(-1.17097512118811d+0*t1*t3*
     =         t4*t6*t61*t9+1.5226864744792862d-1*gammabb*t1*t27*t3*t4*t
     >         60*t9-1.2161642581134076d-3*t1*t3*t4*t57*t58*t59*t9-5.269
     ?         388045346497d+0*t1*t21*t3*t4*t6*t63*t7+1.0538776090692993
     @         d+1*t1*t3*t4*t48*t6*t62*t7-3.512925363564331d+0*t1*t21*t2
     1         9*t3*t30*t4*t6+1.0151243163195242d-1*gammabb*t1*t21*t27*t
     2         28*t3*t30*t4)*wght-1.6672336067107812d-1*rhoa*t33*wght)+t
     3         14*t16*(1.6672336067107812d-1*rhoa*rhob*t15*t64*wght+1.66
     4         72336067107812d-1*rhoa*rhob*t31*t33*wght+1.66723360671078
     5         12d-1*rhoa*t15*t31*wght-2.0d+0*t31*t34)-5.001700820132344
     6         d-1*rhoa*rhob*t14*t15*t31**2*t46*wght-1.6672336067107812d
     7         -1*rhoa*rhob*t12*t15*t32**2*t45*wght-3.3344672134215625d-
     8         1*rhoa*rhob*t15*t16*t23*t31*t32*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t12*t23*(8.336168033553906d-2*rhoa*rho
     1         b*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t70*t9+3.3825
     2         41776131776d-4*gammaaa*t1*t50*t6*t68*t7*t9+8.666103204112
     3         231d-2*t1*t17*t22*t36*t6*t7*t9-3.293949630991462d-2*t1*t1
     4         7*t18*t6*t7*t9-3.908227599821866d+0*t1*t21*t3*t4*t6*t69*t
     5         7-8.75744258821479d-3*t1*t17*t21*t3*t6*t67*t7+2.530623910
     6         2811573d-4*gammaaa*t1*t21*t6*t65*t66*t7+5.254465552928875
     7         d-2*t1*t17*t20*t3*t35*t48*t6*t7-1.4116926989963413d-2*t1*
     8         t17*t20*t21*t36*t6*t7-1.612806413056766d-1*t1*t17*t21*t22
     9         *t3*t35*t6*t7)*wght+8.336168033553906d-2*rhoa*rhob*t24*t3
     :         7*wght-t26*t39)+t12*t14*(-8.336168033553906d-2*rhoa*rhob*
     ;         (4.560615967925279d-4*gammaaa*t1*t50*t6*t68*t7*t9-4.44116
     <         8883897918d-2*t1*t17*t18*t6*t7*t9-5.269388045346497d+0*t1
     =         *t21*t3*t4*t6*t69*t7-1.1807491274113775d-2*t1*t17*t21*t3*
     >         t6*t67*t7+3.411991507534344d-4*gammaaa*t1*t21*t6*t65*t66*
     ?         t7+7.084494764468265d-2*t1*t17*t20*t3*t35*t48*t6*t7-1.903
     @         3580930991081d-2*t1*t17*t20*t21*t36*t6*t7)*wght-8.3361680
     1         33553906d-2*rhob*t37*wght)+t14*t16*(1.6672336067107812d-1
     2         *rhoa*rhob*t15*t70*wght+1.6672336067107812d-1*rhoa*rhob*t
     3         22*t37*wght-2.0d+0*t26*t38)+t16*t23*(-1.6672336067107812d
     4         -1*rhoa*rhob*t15*t22*t39*wght-1.6672336067107812d-1*rhoa*
     5         rhob*t15*t24*t38*wght)-5.001700820132344d-1*rhoa*rhob*t14
     6         *t15*t22*t38*t46*wght-1.6672336067107812d-1*rhoa*rhob*t12
     7         *t15*t24*t39*t45*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t12*t23*(8.336168033553906d-2*rhoa*rho
     1         b*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t71*t9+8.6661
     2         03204112231d-2*t1*t22*t27*t3*t4*t40*t9+4.705642329987803d
     3         -3*t1*t19*t27*t3*t40*t9-1.3597817940049742d-4*gammaaa*t1*
     4         t17*t18*t27*t40*t9+5.254465552928875d-2*t1*t20*t27*t3*t4*
     5         t41*t48*t6-1.612806413056766d-1*t1*t21*t22*t27*t3*t4*t41*
     6         t6-8.75744258821479d-3*t1*t19*t21*t27*t3*t41*t6+2.5306239
     7         102811573d-4*gammaaa*t1*t17*t18*t21*t27*t41*t6-1.41169269
     8         89963413d-2*t1*t20*t21*t27*t3*t4*t40)*wght+8.336168033553
     9         906d-2*rhoa*rhob*t24*t42*wght-t26*t44)+t12*t14*(-8.336168
     :         033553906d-2*rhoa*rhob*(6.344526976997026d-3*t1*t19*t27*t
     ;         3*t40*t9-1.8333676191059622d-4*gammaaa*t1*t17*t18*t27*t40
     <         *t9+7.084494764468265d-2*t1*t20*t27*t3*t4*t41*t48*t6-1.18
     =         07491274113775d-2*t1*t19*t21*t27*t3*t41*t6+3.411991507534
     >         344d-4*gammaaa*t1*t17*t18*t21*t27*t41*t6-1.90335809309910
     ?         81d-2*t1*t20*t21*t27*t3*t4*t40)*wght-8.336168033553906d-2
     @         *rhob*t42*wght)+t14*t16*(1.6672336067107812d-1*rhoa*rhob*
     1         t15*t71*wght+1.6672336067107812d-1*rhoa*rhob*t22*t42*wght
     2         -2.0d+0*t26*t43)+t16*t23*(-1.6672336067107812d-1*rhoa*rho
     3         b*t15*t22*t44*wght-1.6672336067107812d-1*rhoa*rhob*t15*t2
     4         4*t43*wght)-5.001700820132344d-1*rhoa*rhob*t14*t15*t22*t4
     5         3*t46*wght-1.6672336067107812d-1*rhoa*rhob*t12*t15*t24*t4
     6         4*t45*wght+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t12*t23*(8.336168033553906d-2*rhoa*rho
     1         b*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t72*t9+8.6661
     2         03204112231d-2*t1*t17*t31*t36*t6*t7*t9+4.705642329987803d
     3         -3*t1*t17*t29*t36*t6*t9-1.3597817940049742d-4*gammabb*t1*
     4         t17*t27*t28*t36*t9+5.254465552928875d-2*t1*t17*t3*t30*t35
     5         *t48*t6*t7-1.4116926989963413d-2*t1*t17*t21*t30*t36*t6*t7
     6         -1.612806413056766d-1*t1*t17*t21*t3*t31*t35*t6*t7-8.75744
     7         258821479d-3*t1*t17*t21*t29*t3*t35*t6+2.5306239102811573d
     8         -4*gammabb*t1*t17*t21*t27*t28*t3*t35)*wght+8.336168033553
     9         906d-2*rhoa*rhob*t32*t37*wght-t34*t39)+t12*t14*(-8.336168
     :         033553906d-2*rhoa*rhob*(6.344526976997026d-3*t1*t17*t29*t
     ;         36*t6*t9-1.8333676191059622d-4*gammabb*t1*t17*t27*t28*t36
     <         *t9+7.084494764468265d-2*t1*t17*t3*t30*t35*t48*t6*t7-1.90
     =         33580930991081d-2*t1*t17*t21*t30*t36*t6*t7-1.180749127411
     >         3775d-2*t1*t17*t21*t29*t3*t35*t6+3.411991507534344d-4*gam
     ?         mabb*t1*t17*t21*t27*t28*t3*t35)*wght-8.336168033553906d-2
     @         *rhoa*t37*wght)+t14*t16*(1.6672336067107812d-1*rhoa*rhob*
     1         t15*t72*wght+1.6672336067107812d-1*rhoa*rhob*t31*t37*wght
     2         -2.0d+0*t34*t38)+t16*t23*(-1.6672336067107812d-1*rhoa*rho
     3         b*t15*t31*t39*wght-1.6672336067107812d-1*rhoa*rhob*t15*t3
     4         2*t38*wght)-5.001700820132344d-1*rhoa*rhob*t14*t15*t31*t3
     5         8*t46*wght-1.6672336067107812d-1*rhoa*rhob*t12*t15*t32*t3
     6         9*t45*wght+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t12*t23*(8.336168033553906d-2*rhoa*rho
     1         b*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t78*t9+3.3825
     2         41776131776d-4*gammabb*t1*t3*t4*t58*t73*t9+8.666103204112
     3         231d-2*t1*t27*t3*t31*t4*t40*t9-3.293949630991462d-2*t1*t2
     4         7*t28*t3*t4*t9-3.908227599821866d+0*t1*t21*t3*t4*t6*t7*t7
     5         7-8.75744258821479d-3*t1*t21*t27*t3*t4*t6*t76+2.530623910
     6         2811573d-4*gammabb*t1*t21*t3*t4*t74*t75+5.254465552928875
     7         d-2*t1*t27*t3*t30*t4*t41*t48*t6-1.612806413056766d-1*t1*t
     8         21*t27*t3*t31*t4*t41*t6-1.4116926989963413d-2*t1*t21*t27*
     9         t3*t30*t4*t40)*wght+8.336168033553906d-2*rhoa*rhob*t32*t4
     :         2*wght-t34*t44)+t12*t14*(-8.336168033553906d-2*rhoa*rhob*
     ;         (4.560615967925279d-4*gammabb*t1*t3*t4*t58*t73*t9-4.44116
     <         8883897918d-2*t1*t27*t28*t3*t4*t9-5.269388045346497d+0*t1
     =         *t21*t3*t4*t6*t7*t77-1.1807491274113775d-2*t1*t21*t27*t3*
     >         t4*t6*t76+3.411991507534344d-4*gammabb*t1*t21*t3*t4*t74*t
     ?         75+7.084494764468265d-2*t1*t27*t3*t30*t4*t41*t48*t6-1.903
     @         3580930991081d-2*t1*t21*t27*t3*t30*t4*t40)*wght-8.3361680
     1         33553906d-2*rhoa*t42*wght)+t14*t16*(1.6672336067107812d-1
     2         *rhoa*rhob*t15*t78*wght+1.6672336067107812d-1*rhoa*rhob*t
     3         31*t42*wght-2.0d+0*t34*t43)+t16*t23*(-1.6672336067107812d
     4         -1*rhoa*rhob*t15*t31*t44*wght-1.6672336067107812d-1*rhoa*
     5         rhob*t15*t32*t43*wght)-5.001700820132344d-1*rhoa*rhob*t14
     6         *t15*t31*t43*t46*wght-1.6672336067107812d-1*rhoa*rhob*t12
     7         *t15*t32*t44*t45*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t12*t23*(8.336168033553906d-2*rhoa*rh
     1         ob*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t81*t9-1.268
     2         453166049416d-4*t1*t50*t6*t7*t80*t9+8.666103204112231d-2*
     3         t1*t17*t36*t38*t6*t7*t9-1.8979679327108684d-4*t1*t21*t6*t
     4         65*t7*t79+2.3606566327249603d-4*t1*t21*t3*t50*t6*t7*t79+3
     5         .53221601630552d-4*t1*t3*t48*t52*t6*t65*t7-1.612806413056
     6         766d-1*t1*t17*t21*t3*t35*t38*t6*t7)*wght+1.66723360671078
     7         12d-1*rhoa*rhob*t37*t39*wght)+t14*t16*(1.6672336067107812
     8         d-1*rhoa*rhob*t15*t81*wght+3.3344672134215625d-1*rhoa*rho
     9         b*t37*t38*wght)-8.336168033553906d-2*rhoa*rhob*t12*t14*(-
     :         1.71023098797198d-4*t1*t50*t6*t7*t80*t9-2.558993630650758
     ;         d-4*t1*t21*t6*t65*t7*t79+3.1828278988193504d-4*t1*t21*t3*
     <         t50*t6*t7*t79+4.762418864947907d-4*t1*t3*t48*t52*t6*t65*t
     =         7)*wght-5.001700820132344d-1*rhoa*rhob*t14*t15*t38**2*t46
     >         *wght-1.6672336067107812d-1*rhoa*rhob*t12*t15*t39**2*t45*
     ?         wght-3.3344672134215625d-1*rhoa*rhob*t15*t16*t23*t38*t39*
     @         wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t12*t23*(8.336168033553906d-2*rhoa*rh
     1         ob*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t82*t9+8.666
     2         103204112231d-2*t1*t27*t3*t38*t4*t40*t9+5.099181727518654
     3         d-5*t1*t17*t27*t36*t40*t9+3.53221601630552d-4*t1*t17*t27*
     4         t3*t35*t41*t48*t6-1.612806413056766d-1*t1*t21*t27*t3*t38*
     5         t4*t41*t6-9.489839663554341d-5*t1*t17*t21*t27*t36*t41*t6-
     6         9.489839663554341d-5*t1*t17*t21*t27*t3*t35*t40)*wght+8.33
     7         6168033553906d-2*rhoa*rhob*t37*t44*wght+8.336168033553906
     8         d-2*rhoa*rhob*t39*t42*wght)+t14*t16*(1.6672336067107812d-
     9         1*rhoa*rhob*t15*t82*wght+1.6672336067107812d-1*rhoa*rhob*
     :         t37*t43*wght+1.6672336067107812d-1*rhoa*rhob*t38*t42*wght
     ;         )+t16*t23*(-1.6672336067107812d-1*rhoa*rhob*t15*t38*t44*w
     <         ght-1.6672336067107812d-1*rhoa*rhob*t15*t39*t43*wght)-8.3
     =         36168033553906d-2*rhoa*rhob*t12*t14*(6.875128571647359d-5
     >         *t1*t17*t27*t36*t40*t9+4.762418864947907d-4*t1*t17*t27*t3
     ?         *t35*t41*t48*t6-1.279496815325379d-4*t1*t17*t21*t27*t36*t
     @         41*t6-1.279496815325379d-4*t1*t17*t21*t27*t3*t35*t40)*wgh
     1         t-5.001700820132344d-1*rhoa*rhob*t14*t15*t38*t43*t46*wght
     2         -1.6672336067107812d-1*rhoa*rhob*t12*t15*t39*t44*t45*wght
     3         +Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t12*t23*(8.336168033553906d-2*rhoa*rh
     1         ob*t15*(2.3991838839497967d+1*t1*t3*t4*t6*t7*t85*t9-1.268
     2         453166049416d-4*t1*t3*t4*t58*t83*t9+8.666103204112231d-2*
     3         t1*t27*t3*t4*t40*t43*t9-1.8979679327108684d-4*t1*t21*t3*t
     4         4*t74*t84+2.3606566327249603d-4*t1*t21*t3*t4*t58*t6*t84+3
     5         .53221601630552d-4*t1*t3*t4*t48*t6*t60*t74-1.612806413056
     6         766d-1*t1*t21*t27*t3*t4*t41*t43*t6)*wght+1.66723360671078
     7         12d-1*rhoa*rhob*t42*t44*wght)+t14*t16*(1.6672336067107812
     8         d-1*rhoa*rhob*t15*t85*wght+3.3344672134215625d-1*rhoa*rho
     9         b*t42*t43*wght)-8.336168033553906d-2*rhoa*rhob*t12*t14*(-
     :         1.71023098797198d-4*t1*t3*t4*t58*t83*t9-2.558993630650758
     ;         d-4*t1*t21*t3*t4*t74*t84+3.1828278988193504d-4*t1*t21*t3*
     <         t4*t58*t6*t84+4.762418864947907d-4*t1*t3*t4*t48*t6*t60*t7
     =         4)*wght-5.001700820132344d-1*rhoa*rhob*t14*t15*t43**2*t46
     >         *wght-1.6672336067107812d-1*rhoa*rhob*t12*t15*t44**2*t45*
     ?         wght-3.3344672134215625d-1*rhoa*rhob*t15*t16*t23*t43*t44*
     @         wght+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
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
C>   {\it t_4} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_5} &=& 1.804-{{0.804}\over{{{0.004492667444378077
C>    \,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\\\\ 
C>   {\it t_6} &=& {{1}\over{1.8610514726982\,{\it t_3}\,{\it t_5}
C>    +1.8610514726982\,{\it t_2}\,{\it t_4}}}\\\\ 
C>   {\it t_7} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}
C>    \,{\it t_6}\,{\it t_5}\right)^{2.0}\\\\ 
C>   f &=& -{{0.08336168033553906\,\rho_\alpha\,\rho_\beta\,
C>    \left(5.269388045346497\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_4}\,{\it t_6}\,{\it t_5}+0.5764\right)}\over{{\it t_7}
C>    \,\left(11.99591941974898\,{\it t_7}+3.908227599821866\,{
C>    \it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_4}\,{\it t_6}\,{\it t_5}
C>    +0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 111
C>    , 5656 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479954 "> 
C> 10.1063/1.479954 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_pbeop_d3(param,tol_rho,ipol,nq,wght, 
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
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.492667444378077d-3*gammaaa*t2+1.0d+0
            t4 = 1.804d+0-8.04d-1/t3
            t5 = rhoa**3.333333333333333d-1
            t6 = t1*t4*t5
            t7 = t6**2.0d+0
            t8 = 1/t7
            t9 = 1.0500052408963245d+0*t1*t4*t5
            t10 = 1.415701855281521d+0*t1*t4*t5+5.764d-1
            t11 = rhoa**2
            t12 = 1/t3**2
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**6.666666666666666d-1
            t15 = 6.203504908994d-1*t14*t4-1.7926167019779396d-2*gammaaa
     1         *t12*t13
            t16 = 8.955510855647054d-2*t1*t14*t4-7.218105721705101d-2*t1
     1         *t15-2.587859374674219d-3*gammaaa*t1*t12*t13
            t17 = 1/t6**3.0d+0
            t18 = t9+8.658781460080348d-1*t7+3.183d-1
            t19 = 1/t18
            t20 = 1/t18**2
            t21 = 6.44577519522176d+0*t1*t16*t4*t5+3.5000174696544145d-1
     1         *t1*t14*t4-2.821d-1*t1*t15-1.0113943432559537d-2*gammaaa*
     2         t1*t12*t13
            t22 = 4.719006184271736d-1*t1*t14*t4-3.8035d-1*t1*t15-1.3636
     1         43525194619d-2*gammaaa*t1*t12*t13
            t23 = -8.336168033553906d-2*t11*t22*wght-8.336168033553906d-
     1         2*rhoa*t10*wght
            t24 = 1/t5
            t25 = 1/rhoa**2.3333333333333334d+0
            t26 = 1.8963643936049138d-3*t1*t12*t25
            t27 = t1**2
            t28 = 1/t11
            t29 = 3.1276424561244703d-3*t12*t27*t28*t4+t26
            t30 = 1/t18**3
            t31 = t21**2
            t32 = -2.3097938309456317d+0*t10*t11*t30*t31*t8*wght
            t33 = t16**2
            t34 = 1/t6**4.0d+0
            t35 = -9.599999999999997d+1*t10*t11*t19*t33*t34*wght
            t36 = -1.7194580842842339d+1*t10*t11*t16*t17*t20*t21*wght
            t37 = -1.6672336067107812d-1*rhoa*t22*wght
            t38 = gammaaa**2
            t39 = 1/t3**3
            t40 = 1/rhoa**7
            t41 = 1/rhoa**4.333333333333333d+0
            t42 = 1/t4
            t43 = 1/rhoa**3.6666666666666664d+0
            t44 = 7.327274635867936d-3*gammaaa*t1*t12*t15*t42*t43
            t45 = 1/rhoa**1.6666666666666669d+0
            t46 = -4.135669939329333d-1*t4*t45+5.377850105933819d-2*gamm
     1         aaa*t12*t41-4.295269705186499d-4*t38*t39*t40
            t47 = 1/rhoa
            t48 = -2.5356666666666666d-1*t1*t15*t47
            t49 = t15**2
            t50 = 2.0437371323672143d-1*t1*t24*t42*t49
            t51 = t50+t48-3.8035d-1*t1*t46-3.1460041228478236d-1*t1*t4*t
     1         45+t44+4.090930575583857d-2*gammaaa*t1*t12*t41-3.26741166
     2         47353707d-4*t1*t38*t39*t40
            t52 = t37-8.336168033553906d-2*t11*t51*wght
            t53 = 1.6672336067107812d-1*rhoa*t10*t16*wght
            t54 = 1.6672336067107812d-1*t11*t16*t22*wght
            t55 = 1.390536163366978d-3*gammaaa*t1*t12*t15*t42*t43
            t56 = -4.812070481136733d-2*t1*t15*t47
            t57 = 3.8785094488762883d-2*t1*t24*t42*t49
            t58 = t57+t56+t55-7.218105721705101d-2*t1*t46-5.970340570431
     1         37d-2*t1*t4*t45+7.76357812402266d-3*gammaaa*t1*t12*t41-6.
     2         20074216705465d-5*t1*t38*t39*t40
            t59 = -2.0d+0*t16*t23
            t60 = 1.6672336067107812d-1*t10*t11*t58*wght+t59+t54+t53
            t61 = 8.336168033553906d-2*rhoa*t10*t21*wght
            t62 = 8.336168033553906d-2*t11*t21*t22*wght
            t63 = 5.434531812221229d-3*gammaaa*t1*t12*t15*t42*t43
            t64 = -6.208750505646785d-2*gammaaa*t1*t12*t13*t16
            t65 = -1.8806666666666666d-1*t1*t15*t47
            t66 = 2.1485917317405864d+0*t1*t14*t16*t4
            t67 = -1.7317562920160695d+0*t1*t15*t16
            t68 = 1.5158097674268203d-1*t1*t24*t42*t49
            t69 = t68+t67+t66+t65+t64+t63+6.44577519522176d+0*t1*t4*t5*t
     1         58-2.821d-1*t1*t46-2.3333449797696096d-1*t1*t4*t45+3.0341
     2         83029767861d-2*gammaaa*t1*t12*t41-2.4233911676662234d-4*t
     3         1*t38*t39*t40
            t70 = -t21*t23
            t71 = 8.336168033553906d-2*t10*t11*t69*wght+t70+t62+t61
            t72 = 1/t3**4
            t73 = 1.313499489223619d-4*t1*t38*t40*t42*t72+t50+t48+1.5730
     1         020614239118d-1*t1*t4*t45+t44-9.090956834630795d-3*gammaa
     2         a*t1*t12*t41
            t74 = -8.336168033553906d-2*t11*t73*wght-8.336168033553906d-
     1         2*t10*wght+t37
            t75 = 2.4926983511559697d-5*t1*t38*t40*t42*t72+t57+t56+t55+2
     1         .985170285215685d-2*t1*t4*t45-1.7252395831161463d-3*gamma
     2         aa*t1*t12*t41
            t76 = 1.6672336067107812d-1*t10*t11*t75*wght+t59+t54+t53
            t77 = 6.44577519522176d+0*t1*t4*t5*t75+9.742032494018218d-5*
     1         t1*t38*t40*t42*t72+t68+t67+t66+t65+t64+t63+1.166672489884
     2         8048d-1*t1*t4*t45-6.742628955039691d-3*gammaaa*t1*t12*t41
            t78 = 8.336168033553906d-2*t10*t11*t77*wght+t70+t62+t61
            t79 = -4.6581468744135934d-2*t1*t10*t12*t16*t19*t24*t34*wght
            t80 = 1/rhoa**1.3333333333333333d+0
            t81 = -2.1314177932294023d-4*t1*t12*t80*wght
            t82 = 1/rhoa**6
            t83 = 1.610726139444937d-4*gammaaa*t39*t82-1.568539614230697
     1         3d-2*t12*t13
            t84 = -3.8035d-1*t1*t83+2.4628115422942856d-5*gammaaa*t1*t42
     1         *t72*t82+1.2252793742757642d-4*gammaaa*t1*t39*t82-1.27841
     2         58048699554d-2*t1*t12*t13
            t85 = t81-8.336168033553906d-2*t11*t84*wght
            t86 = -8.08981147293451d-5*t1*t10*t12*t21*t24*wght
            t87 = t86-1.6672336067107812d-1*t10*t11*t16*t29*wght
            t88 = 4.2628355864588047d-4*t1*t12*t16*t24*wght
            t89 = -7.218105721705101d-2*t1*t83+4.673809408417444d-6*gamm
     1         aaa*t1*t42*t72*t82+2.3252783126454943d-5*gammaaa*t1*t39*t
     2         82-2.4261181637570808d-3*t1*t12*t13
            t90 = -9.704472655028322d-4*t1*t12*t23*t25
            t91 = 1.6672336067107812d-1*t10*t11*t89*wght+t90+t88
            t92 = 2.1314177932294023d-4*t1*t12*t21*t24*wght
            t93 = -2.168404344971009d-19*t1*t12*t15*t2*t42
            t94 = 1.1641407198087722d-2*t1*t12*t16*t25
            t95 = t94+t93+6.44577519522176d+0*t1*t4*t5*t89-2.821d-1*t1*t
     1         83+1.8266310926284157d-5*gammaaa*t1*t42*t72*t82+9.0877168
     2         78748338d-5*gammaaa*t1*t39*t82-9.481821968024564d-3*t1*t1
     3         2*t13
            t96 = 8.336168033553906d-2*t10*t11*t95*wght+t92-t23*t29
            t97 = 3.127642456124471d-3*t12*t27*t28*t4+t26
            t98 = 8.522772032466369d-4*t1*t12*t13-2.4628115422942862d-5*
     1         gammaaa*t1*t42*t72*t82
            t99 = t81-8.336168033553906d-2*t11*t98*wght
            t100 = t86-1.6672336067107812d-1*t10*t11*t16*t97*wght
            t101 = 1.617412109171387d-4*t1*t12*t13-4.673809408417444d-6*
     1         gammaaa*t1*t42*t72*t82
            t102 = 1.6672336067107812d-1*t10*t101*t11*wght+t90+t88
            t103 = t94+t93-1.8266310926284166d-5*gammaaa*t1*t42*t72*t82+
     1         6.44577519522176d+0*t1*t101*t4*t5+6.321214645349711d-4*t1
     2         *t12*t13
            t104 = 8.336168033553906d-2*t10*t103*t11*wght-t23*t97+t92
            t105 = 1/rhoa**5
            t106 = -9.235543283603571d-6*t1*t105*t42*t72-2.2973988267670
     1         578d-5*t1*t105*t39
            t107 = t29**2
            t108 = 4.1368571381671093d-7*t2*t27*t72*wght
            t109 = -1.7526785281565413d-6*t1*t105*t42*t72-4.359896836210
     1         302d-6*t1*t105*t39
            t110 = 1.6672336067107812d-1*t10*t109*t11*wght+t108
            t111 = 1/rhoa**4.666666666666667d+0
            t112 = -6.84986659735656d-6*t1*t105*t42*t72+5.64868589099460
     1         9d-6*t111*t27*t72+6.44577519522176d+0*t1*t109*t4*t5-1.703
     2         9469147653133d-5*t1*t105*t39
            t113 = 4.2628355864588047d-4*t1*t12*t24*t29*wght+8.336168033
     1         553906d-2*t10*t11*t112*wght
            t114 = 1/rhoa**3
            t115 = 6.849866597356562d-6*t1*t105*t42*t72+1.69460576729838
     1         27d-5*t111*t27*t72
            t116 = 5.158374252852702d+1*t10*t11*t16*t17*t30*t31*wght
            t117 = 1/t18**4
            t118 = 6.929381492836896d+0*t10*t11*t117*t21**3*t8*wght
            t119 = 1/t6**5.0d+0
            t120 = 1.4292875310322173d+3*t10*t11*t119*t16**3*t19*wght
            t121 = 2.879999999999999d+2*t10*t11*t20*t21*t33*t34*wght
            t122 = gammaaa**3
            t123 = 1/rhoa**1.0666666666666666d+1
            t124 = 1/rhoa**8
            t125 = 1/rhoa**7.333333333333333d+0
            t126 = 1/rhoa**5.333333333333333d+0
            t127 = 1/t4**2
            t128 = 1/rhoa**4
            t129 = -1.5437774695451842d-5*t122*t123*t72+6.89278323221555
     1         5d-1*t2*t4+4.2952697051864996d-3*t124*t38*t39-2.290565785
     2         860701d-1*gammaaa*t12*t126
            t130 = t15**3
            t131 = -1.6472439067503214d-1*t1*t127*t130*t14
            t132 = 3.3344672134215625d-1*rhoa*t16*t22*wght
            t133 = 3.8785094488762883d-2*t1*t42*t49*t80
            t134 = -3.126063012582715d-2*t1*t127*t130*t14
            t135 = -2.2286297971927035d-6*t1*t122*t123*t72-1.12076655355
     1         82496d-3*gammaaa*t1*t12*t127*t128*t49-7.218105721705101d-
     2         2*t1*t46*t47+2.0858042450504666d-3*gammaaa*t1*t12*t42*t43
     3         *t46+1.1635528346628865d-1*t1*t15*t24*t42*t46+4.997773241
     4         1113723d-5*t1*t125*t15*t38*t39*t42-6.257412735151401d-3*g
     5         ammaaa*t1*t111*t12*t15*t42+9.950567617385617d-2*t1*t2*t4+
     6         6.200742167054651d-4*t1*t124*t38*t39+4.812070481136733d-2
     7         *t1*t15*t28+t134+t133-7.218105721705101d-2*t1*t129-3.3067
     8         09200972614d-2*gammaaa*t1*t12*t126
            t136 = 1.6672336067107812d-1*rhoa*t21*t22*wght
            t137 = -4.3802107775735827d-3*gammaaa*t1*t12*t127*t128*t49
            t138 = 3.336151953199435d-2*gammaaa*t1*t12*t15*t16*t42*t43
            t139 = 1.51580976742682d-1*t1*t42*t49*t80
            t140 = -1.1545041946773796d+0*t1*t15*t16*t47
            t141 = -1.2217365744558055d-1*t1*t127*t130*t14
            t142 = 9.305257363490999d-1*t1*t16*t24*t42*t49
            t143 = -5.001700820132344d-1*rhoa*t10*t33*wght
            t144 = -5.001700820132344d-1*t11*t22*t33*wght
            t145 = -3.0d+0*t16*t60
            t146 = -3.3344672134215625d-1*rhoa*t10*t16*t21*wght
            t147 = -3.3344672134215625d-1*t11*t16*t21*t22*wght
            t148 = -t21*t60
            t149 = -2.0d+0*t16*t71
            t150 = -1.6672336067107812d-1*rhoa*t10*t31*wght
            t151 = -1.6672336067107812d-1*t11*t22*t31*wght
            t152 = -2*t21*t71
            t153 = 1/t3**5
            t154 = -7.478095053467911d-5*t1*t124*t38*t42*t72-1.339403228
     1         6178478d-5*t1*t125*t127*t15*t38*t72-1.1207665535582498d-3
     2         *gammaaa*t1*t12*t127*t128*t49-2.4060352405683666d-2*t1*t4
     3         6*t47+6.95268081683489d-4*gammaaa*t1*t12*t42*t43*t46+3.87
     4         85094488762883d-2*t1*t15*t24*t42*t46+1.6659244137037907d-
     5         5*t1*t125*t15*t38*t39*t42+5.972727856476445d-7*t1*t122*t1
     6         23*t153*t42-1.158780136139149d-3*gammaaa*t1*t111*t12*t15*
     7         t42-1.990113523477123d-2*t1*t2*t4-2.0669140556848836d-5*t
     8         1*t124*t38*t39+t134+t133+3.162939235712935d-3*gammaaa*t1*
     9         t12*t126
            t155 = 6.935240880537548d-1*t1*t10*t119*t12*t19*t24*t33*wght
            t156 = 6.929381492836896d+0*t10*t11*t117*t29*t31*t8*wght
            t157 = 1/rhoa**9.666666666666666d+0
            t158 = 1/rhoa**6.333333333333333d+0
            t159 = 5.789165510794442d-6*t157*t38*t72+5.228465380768991d-
     1         2*t12*t41-1.342271782870781d-3*gammaaa*t39*t40
            t160 = -3.6910961743400234d-4*t1*t114*t12*t127*t49
            t161 = 4.853886883760707d-4*t1*t10*t12*t16*t21*t24*wght
            t162 = 5.001700820132344d-1*t10*t11*t29*t33*wght
            t163 = 1.617962294586902d-4*t1*t10*t12*t24*t31*wght
            t164 = 5.156620156177406d+1*t17*t30*(6.668934426843125d-1*t1
     1         0*t11*t16*t21*t29*wght+t163)
            t165 = 4.2628355864588047d-4*t1*t12*t16*t80*wght
            t166 = 4.2628355864588047d-4*t1*t12*t24*t58*wght
            t167 = 1.6672336067107812d-1*rhoa*t10*t89*wght
            t168 = 1.6672336067107812d-1*t11*t22*t89*wght
            t169 = -7.004790959739062d-5*t1*t114*t12*t127*t49
            t170 = -4.812070481136733d-2*t1*t47*t83+1.390536163366978d-3
     1         *gammaaa*t1*t12*t42*t43*t83+7.757018897752577d-2*t1*t15*t
     2         24*t42*t83-1.4021428225252333d-5*gammaaa*t1*t40*t42*t72+8
     3         .357361739472638d-7*t1*t157*t38*t72-5.02276210731693d-6*g
     4         ammaaa*t1*t127*t15*t158*t72+1.390536163366978d-3*t1*t12*t
     5         15*t42*t43-1.2494433102778432d-5*gammaaa*t1*t15*t158*t39*
     6         t42+1.1198864730893335d-7*t1*t153*t157*t38*t42+7.65575065
     7         0077899d-3*t1*t12*t41-1.9377319272045787d-4*gammaaa*t1*t3
     8         9*t40+t169-7.218105721705101d-2*t1*t159
            t171 = -2.0d+0*t16*t85
            t172 = -9.704472655028322d-4*t1*t12*t25*t52
            t173 = 2.1314177932294023d-4*t1*t12*t21*t80*wght
            t174 = 2.1314177932294023d-4*t1*t12*t24*t69*wght
            t175 = 8.336168033553906d-2*rhoa*t10*t95*wght
            t176 = 8.336168033553906d-2*t11*t22*t95*wght
            t177 = -6.208750505646785d-2*gammaaa*t1*t12*t13*t89
            t178 = -2.73763173598349d-4*t1*t114*t12*t127*t49
            t179 = -2.168404344971009d-19*t1*t12*t2*t42*t46
            t180 = 8.673617379884035d-19*t1*t12*t15*t16*t2*t42
            t181 = 1.1641407198087722d-2*t1*t12*t25*t58
            t182 = 2.1485917317405864d+0*t1*t14*t4*t89
            t183 = -1.7317562920160695d+0*t1*t15*t89
            t184 = -t21*t85
            t185 = -1.2788506759376417d-3*t1*t12*t24*t33*wght
            t186 = -1.4556708982542485d-3*t1*t12*t25*t60
            t187 = -8.525671172917609d-4*t1*t12*t16*t21*t24*wght
            t188 = -9.704472655028322d-4*t1*t12*t25*t71
            t189 = -4.2628355864588047d-4*t1*t12*t24*t31*wght
            t190 = 3.3881317890172014d-21*gammaaa*t1*t127*t15*t158*t72
            t191 = 1.6672336067107812d-1*t11*t16*t98*wght
            t192 = 1.4021428225252333d-5*gammaaa*t1*t40*t42*t72-1.119886
     1         4730893338d-7*t1*t153*t157*t38*t42-1.0782747394475914d-4*
     2         t1*t12*t41+t169
            t193 = -2.0d+0*t101*t23
            t194 = -1.1213328209328066d-4*gammaaa*t1*t16*t42*t72*t82
            t195 = 3.880469066029241d-3*t1*t12*t13*t16
            t196 = 8.336168033553906d-2*t11*t21*t98*wght
            t197 = -t103*t23
            t198 = 1/t3**6
            t199 = -2.4060352405683666d-2*t1*t47*t83+6.95268081683489d-4
     1         *gammaaa*t1*t12*t42*t43*t83+3.8785094488762883d-2*t1*t15*
     2         t24*t42*t83+2.4926983511559705d-5*gammaaa*t1*t40*t42*t72-
     3         2.511381053658465d-6*gammaaa*t1*t127*t15*t158*t72+6.95268
     4         081683489d-4*t1*t12*t15*t42*t43-6.247216551389216d-6*gamm
     5         aaa*t1*t15*t158*t39*t42-2.2397729461786672d-7*t1*t153*t15
     6         7*t38*t42-8.087060545856936d-4*t1*t12*t41+7.7509277088183
     7         15d-6*gammaaa*t1*t39*t40-4.5019436218191217d-8*t1*t127*t1
     8         57*t198*t38+t169
            t200 = 3.365142774060559d-4*t10*t119*t16*t19*t2*t27*t72*wght
            t201 = -8.336168033553906d-2*rhoa*t106*wght
            t202 = 1/rhoa**8.666666666666666d+0
            t203 = 3.0201115114592575d-4*t39*t82-2.1709370665479152d-6*g
     1         ammaaa*t202*t72
            t204 = 2.4812702440232998d-6*t1*t126*t127*t15*t72
            t205 = 1.6672336067107812d-1*t106*t11*t16*wght
            t206 = 4.708839475609621d-7*t1*t126*t127*t15*t72
            t207 = 8.763392640782708d-6*t1*t42*t72*t82+4.505226730750645
     1         5d-5*t1*t39*t82-3.13401065230224d-7*gammaaa*t1*t202*t72-1
     2         .2598722822255004d-7*gammaaa*t1*t153*t202*t42+t206-7.2181
     3         05721705101d-2*t1*t203-1.6882288581821706d-8*gammaaa*t1*t
     4         127*t198*t202
            t208 = -2.0d+0*t109*t23
            t209 = 1.840321640170825d-6*t1*t126*t127*t15*t72
            t210 = -4.204998078498024d-5*t1*t105*t16*t42*t72
            t211 = -1.046019422511946d-4*t1*t105*t16*t39
            t212 = 8.336168033553906d-2*t106*t11*t21*wght
            t213 = -t112*t23
            t214 = -1.4556708982542485d-3*t1*t12*t25*t87
            t215 = 1.617962294586902d-4*t1*t10*t12*t21*t24*t29*wght
            t216 = -6.205285707250665d-7*t16*t2*t27*t72*wght
            t217 = -5.001700820132344d-1*t10*t109*t11*t16*wght
            t218 = -1.4556708982542485d-3*t1*t12*t25*t91
            t219 = -2.0684285690835547d-7*t2*t21*t27*t72*wght
            t220 = -4.2628355864588047d-4*t1*t12*t16*t24*t29*wght
            t221 = -1.6672336067107812d-1*t10*t11*t112*t16*wght
            t222 = -1.6672336067107812d-1*t10*t109*t11*t21*wght
            t223 = -9.704472655028322d-4*t1*t12*t25*t96
            t224 = -4.2628355864588047d-4*t1*t12*t21*t24*t29*wght
            t225 = -1.6672336067107812d-1*t10*t11*t112*t21*wght
            t226 = 2.2129186475408202d-7*gammaaa*t1*t153*t202*t42
            t227 = -4.0895832323652637d-6*t1*t42*t72*t82-4.7088394756096
     1         23d-7*t1*t126*t127*t15*t72+4.199574274085001d-8*gammaaa*t
     2         1*t153*t202*t42+3.3087224502121104d-24*gammaaa*t1*t127*t1
     3         98*t202
            t228 = 1.6412892085480882d-7*gammaaa*t1*t153*t202*t42
            t229 = 2.4269434418803534d-4*t1*t10*t12*t16*t24*t97*wght
            t230 = -2.4269434418803534d-4*t1*t10*t101*t12*t24*wght
            t231 = -8.08981147293451d-5*t1*t10*t103*t12*t24*wght
            t232 = -5.842261760521806d-7*t1*t42*t72*t82-1.45329894540343
     1         4d-6*t1*t39*t82+4.1995742740850023d-8*gammaaa*t1*t153*t20
     2         2*t42+t206+1.688228858182171d-8*gammaaa*t1*t127*t198*t202
            t233 = t1**3
            t234 = 1.632846801556842d-7*t10*t105*t119*t19*t198*t233*wght
            t235 = 1/rhoa**7.666666666666667d+0
            t236 = 5.876269973066699d-8*t1*t235*t72+4.7245210583456254d-
     1         8*t1*t153*t235*t42+9.49628732727471d-9*t1*t127*t198*t235
            t237 = -1.4556708982542485d-3*t1*t110*t12*t25
            t238 = -4.1368571381671093d-7*t2*t27*t29*t72*wght
            t239 = -9.704472655028322d-4*t1*t113*t12*t25
            t240 = -4.2628355864588047d-4*t1*t107*t12*t24*wght
            t241 = -1.5748403527818752d-8*t1*t153*t235*t42-3.16542910909
     1         1572d-9*t1*t127*t198*t235
            fnc(iq) = fnc(iq)-1.1548969154728159d+0*t10*t11*t8*wght/(t9+
     1         8.658781460080349d-1*t7+3.183d-1)
            Amat(iq,D1_RA) = 1.1548969154728159d+0*t10*t11*t20*t21*t8*wg
     1         ht+8.597290421421169d+0*t10*t11*t16*t17*t19*wght+1.385405
     2         033612855d+1*t19*t23*t8+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.1548969154728159d+0*t10*t11*t20*t29*t8*w
     1         ght-2.9528769394720183d-3*t1*t12*t19*t24*t8*wght+4.171608
     2         490100933d-3*t1*t10*t12*t17*t19*t24*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.385405033612855d+1*t20*t71*t8+1.38540
     1         5033612855d+1*t19*t52*t8+5.156620156177406d+1*t17*t19*t60
     2         +t36+t35+t32+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.385405033612855d+1*t20*t78*t8+1.38540
     1         5033612855d+1*t19*t74*t8+5.156620156177406d+1*t17*t19*t76
     2         +t36+t35+t32+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = -2.3097938309456317d+0*t10*t11*t21*t29
     1         *t30*t8*wght+1.385405033612855d+1*t20*t8*t96+5.1566201561
     2         77406d+1*t17*t19*t91+5.156620156177406d+1*t17*t20*t87+1.3
     3         85405033612855d+1*t19*t8*t85+t79+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = -2.3097938309456317d+0*t10*t11*t21*t30
     1         *t8*t97*wght+1.385405033612855d+1*t19*t8*t99+1.3854050336
     2         12855d+1*t104*t20*t8+t79+5.156620156177406d+1*t100*t17*t2
     3         0+5.156620156177406d+1*t102*t17*t19+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = -2.3097938309456317d+0*t10*t107*t11*t
     1         30*t8*wght-1.1548969154728159d+0*t106*t11*t19*t8*wght-2.2
     2         60242948292618d-5*t10*t19*t2*t27*t34*t72*wght-8.343216980
     3         201866d-3*t1*t10*t12*t17*t20*t24*t29*wght+1.3854050336128
     4         55d+1*t113*t20*t8+5.156620156177406d+1*t110*t17*t19+Cmat2
     5         (iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = 5.156620156177406d+1*t17*t20*(-8.0898
     1         1147293451d-5*t1*t10*t12*t24*t97*wght-8.08981147293451d-5
     2         *t1*t10*t12*t24*t29*wght)+1.385405033612855d+1*t20*t8*(2.
     3         1314177932294023d-4*t1*t12*t24*t97*wght+2.131417793229402
     4         3d-4*t1*t12*t24*t29*wght+8.336168033553906d-2*t10*t11*t11
     5         5*wght)+5.156620156177406d+1*t17*t19*(2.922124543902974d-
     6         7*t1*t10*t114*t42*t72*wght+t108)-2.3097938309456317d+0*t1
     7         0*t11*t29*t30*t8*t97*wght-1.0666100450949446d-5*t1*t114*t
     8         19*t42*t72*t8*wght-2.2602429482926184d-5*t10*t19*t2*t27*t
     9         34*t72*wght+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.385405033612855d+1*t19*t8*(-8.3361
     1         68033553906d-2*t11*(2.0437371323672138d-1*t1*t42*t49*t80-
     2         1.1743515210830219d-5*t1*t122*t123*t72-5.905753878944034d
     3         -3*gammaaa*t1*t12*t127*t128*t49-3.8034999999999997d-1*t1*
     4         t46*t47+1.0990911953801903d-2*gammaaa*t1*t12*t42*t43*t46+
     5         6.131211397101644d-1*t1*t15*t24*t42*t46+2.633520657006488
     6         4d-4*t1*t125*t15*t38*t39*t42-3.297273586140571d-2*gammaaa
     7         *t1*t111*t12*t15*t42+5.243340204746373d-1*t1*t2*t4+3.2674
     8         116647353707d-3*t1*t124*t38*t39+2.5356666666666666d-1*t1*
     9         t15*t28+t131-3.8035d-1*t1*t129-1.7424333933042352d-1*gamm
     :         aaa*t1*t12*t126)*wght-2.500850410066172d-1*rhoa*t51*wght)
     ;         +1.385405033612855d+1*t20*t8*(8.336168033553906d-2*t10*t1
     <         1*(-8.709992483173931d-6*t1*t122*t123*t72+4.2971834634811
     =         73d+0*t1*t14*t4*t58-3.463512584032139d+0*t1*t15*t58-1.241
     >         750101129357d-1*gammaaa*t1*t12*t13*t58+6.44577519522176d+
     ?         0*t1*t135*t4*t5-2.8209999999999996d-1*t1*t46*t47+8.151797
     @         718331845d-3*gammaaa*t1*t12*t42*t43*t46+4.547429302280461
     1         d-1*t1*t15*t24*t42*t46-1.7317562920160695d+0*t1*t16*t46-1
     2         .4323944878270575d+0*t1*t16*t4*t45+1.9532435318562644d-4*
     3         t1*t125*t15*t38*t39*t42-2.4455393154995533d-2*gammaaa*t1*
     4         t111*t12*t15*t42+1.862625151694036d-1*gammaaa*t1*t12*t16*
     5         t41-1.4876720675725458d-3*t1*t16*t38*t39*t40+3.8889082996
     6         16016d-1*t1*t2*t4+2.4233911676662234d-3*t1*t124*t38*t39+1
     7         .8806666666666666d-1*t1*t15*t28+t142+t141+t140+t139+t138+
     8         t137-2.821d-1*t1*t129-1.2923372163826075d-1*gammaaa*t1*t1
     9         2*t126)*wght+1.6672336067107812d-1*t11*t22*t69*wght+1.667
     :         2336067107812d-1*rhoa*t10*t69*wght+8.336168033553906d-2*t
     ;         11*t21*t51*wght-t23*t69-2*t21*t52+t136)+1.385405033612855
     <         d+1*t30*t8*(-3.3344672134215625d-1*t10*t11*t21*t69*wght+t
     =         152+t151+t150)+5.156620156177406d+1*t17*t20*(-3.334467213
     >         4215625d-1*t10*t11*t16*t69*wght-3.3344672134215625d-1*t10
     ?         *t11*t21*t58*wght+t149+t148+t147+t146)+5.156620156177406d
     @         +1*t17*t19*(3.3344672134215625d-1*t11*t22*t58*wght+3.3344
     1         672134215625d-1*rhoa*t10*t58*wght+1.6672336067107812d-1*t
     2         11*t16*t51*wght+1.6672336067107812d-1*t10*t11*t135*wght-2
     3         .0d+0*t23*t58-4.0d+0*t16*t52+t132)+1.9193471071598367d+2*
     4         t19*t34*(-1.0003401640264689d+0*t10*t11*t16*t58*wght+t145
     5         +t144+t143)+t121+t120+t118+t116+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 1.385405033612855d+1*t19*t8*(-8.3361
     1         68033553906d-2*t11*(2.0437371323672143d-1*t1*t42*t49*t80-
     2         3.940498467670857d-4*t1*t124*t38*t42*t72-7.05783536077738
     3         8d-5*t1*t125*t127*t15*t38*t72-5.905753878944036d-3*gammaa
     4         a*t1*t12*t127*t128*t49-1.2678333333333333d-1*t1*t46*t47+3
     5         .663637317933968d-3*gammaaa*t1*t12*t42*t43*t46+2.04373713
     6         23672143d-1*t1*t15*t24*t42*t46+8.778402190021627d-5*t1*t1
     7         25*t15*t38*t39*t42+3.147262076502499d-6*t1*t122*t123*t153
     8         *t42-6.106062196556612d-3*gammaaa*t1*t111*t12*t15*t42-1.0
     9         486680409492745d-1*t1*t2*t4-1.089137221578457d-4*t1*t124*
     :         t38*t39+t131+1.666675419682312d-2*gammaaa*t1*t12*t126)*wg
     ;         ht-1.6672336067107812d-1*rhoa*t73*wght-8.336168033553906d
     <         -2*rhoa*t51*wght-1.6672336067107812d-1*t22*wght)+1.385405
     =         033612855d+1*t20*t8*(8.336168033553906d-2*t11*t22*t77*wgh
     >         t+8.336168033553906d-2*rhoa*t10*t77*wght+8.33616803355390
     ?         6d-2*t10*t11*(2.1485917317405864d+0*t1*t14*t4*t75-1.73175
     @         62920160695d+0*t1*t15*t75-6.208750505646785d-2*gammaaa*t1
     1         *t12*t13*t75+5.980441711641634d-4*t1*t16*t38*t40*t42*t72-
     2         2.922609748205466d-4*t1*t124*t38*t42*t72-5.23469266537478
     3         9d-5*t1*t125*t127*t15*t38*t72+2.1485917317405864d+0*t1*t1
     4         4*t4*t58-1.7317562920160695d+0*t1*t15*t58-6.2087505056467
     5         85d-2*gammaaa*t1*t12*t13*t58+6.44577519522176d+0*t1*t154*
     6         t4*t5-9.403333333333332d-2*t1*t46*t47+2.7172659061106147d
     7         -3*gammaaa*t1*t12*t42*t43*t46+1.5158097674268203d-1*t1*t1
     8         5*t24*t42*t46+7.161972439135288d-1*t1*t16*t4*t45+6.510811
     9         772854215d-5*t1*t125*t15*t38*t39*t42+2.3342779854906137d-
     :         6*t1*t122*t123*t153*t42-4.528776510184358d-3*gammaaa*t1*t
     ;         111*t12*t15*t42-4.1391670037645234d-2*gammaaa*t1*t12*t16*
     <         t41-7.777816599232033d-2*t1*t2*t4-8.077970558887411d-5*t1
     =         *t124*t38*t39+t142+t141+t140+t139+t138+t137+1.23614864175
     >         72767d-2*gammaaa*t1*t12*t126)*wght+8.336168033553906d-2*t
     ?         11*t21*t73*wght+8.336168033553906d-2*t11*t22*t69*wght+8.3
     @         36168033553906d-2*rhoa*t10*t69*wght+8.336168033553906d-2*
     1         t10*t21*wght-t23*t77-t21*t74-t21*t52+t136)+1.385405033612
     2         855d+1*t30*t8*(-3.3344672134215625d-1*t10*t11*t21*t77*wgh
     3         t+t152+t151+t150)+5.156620156177406d+1*t17*t20*(-3.334467
     4         2134215625d-1*t10*t11*t16*t77*wght-3.3344672134215625d-1*
     5         t10*t11*t21*t75*wght+t149+t148+t147+t146)+5.1566201561774
     6         06d+1*t17*t19*(1.6672336067107812d-1*t11*t22*t75*wght+1.6
     7         672336067107812d-1*rhoa*t10*t75*wght+1.6672336067107812d-
     8         1*t11*t16*t73*wght+1.6672336067107812d-1*t11*t22*t58*wght
     9         +1.6672336067107812d-1*rhoa*t10*t58*wght+1.66723360671078
     :         12d-1*t10*t16*wght+1.6672336067107812d-1*t10*t11*t154*wgh
     ;         t-2.0d+0*t23*t75-2.0d+0*t16*t74-2.0d+0*t16*t52+t132)+1.91
     <         93471071598367d+2*t19*t34*(-1.0003401640264689d+0*t10*t11
     =         *t16*t75*wght+t145+t144+t143)+t121+t120+t118+t116+Amat3(i
     >         q,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.385405033612855d+1*t30*t8*(-3.334
     1         4672134215625d-1*t10*t11*t21*t95*wght-2*t29*t71+t189)+5.1
     2         56620156177406d+1*t17*t20*(-3.3344672134215625d-1*t10*t11
     3         *t16*t95*wght-3.3344672134215625d-1*t10*t11*t21*t89*wght-
     4         t29*t60+t188+t187)+1.9193471071598367d+2*t19*t34*(-1.0003
     5         401640264689d+0*t10*t11*t16*t89*wght+t186+t185)+1.3854050
     6         33612855d+1*t20*t8*(8.336168033553906d-2*t11*t21*t84*wght
     7         +8.336168033553906d-2*t10*t11*(-1.8806666666666666d-1*t1*
     8         t47*t83+5.434531812221229d-3*gammaaa*t1*t12*t42*t43*t83+3
     9         .0316195348536407d-1*t1*t15*t24*t42*t83-1.731756292016069
     :         5d+0*t1*t16*t83+1.1213328209328063d-4*gammaaa*t1*t16*t42*
     ;         t72*t82+5.578770253397046d-4*gammaaa*t1*t16*t39*t82-5.479
     <         893277885247d-5*gammaaa*t1*t40*t42*t72+3.2662471811902244
     =         d-6*t1*t157*t38*t72-1.963009749515546d-5*gammaaa*t1*t127*
     >         t15*t158*t72+6.44577519522176d+0*t1*t170*t4*t5+5.43453181
     ?         2221229d-3*t1*t12*t15*t42*t43-4.883108829640661d-5*gammaa
     @         a*t1*t15*t158*t39*t42+4.3767712227949007d-7*t1*t153*t157*
     1         t38*t42+2.9920415987988636d-2*t1*t12*t41-7.57309739895694
     2         7d-4*gammaaa*t1*t39*t40+t183+t182+t181+t180+t179+t178+t17
     3         7-5.820703599043862d-2*t1*t12*t13*t16-2.821d-1*t1*t159)*w
     4         ght-t23*t95-t29*t52+t184+t176+t175+t174+t173)+5.156620156
     5         177406d+1*t17*t19*(1.6672336067107812d-1*t11*t16*t84*wght
     6         +1.6672336067107812d-1*t10*t11*t170*wght-2.0d+0*t23*t89+t
     7         172+t171+t168+t167+t166+t165)+1.385405033612855d+1*t19*t8
     8         *(-1.6672336067107812d-1*rhoa*t84*wght-8.336168033553906d
     9         -2*t11*(-2.5356666666666666d-1*t1*t47*t83+7.3272746358679
     :         36d-3*gammaaa*t1*t12*t42*t43*t83+4.0874742647344287d-1*t1
     ;         *t15*t24*t42*t83-7.388434626882856d-5*gammaaa*t1*t40*t42*
     <         t72+4.403818204061332d-6*t1*t157*t38*t72-2.64668826029152
     =         04d-5*gammaaa*t1*t127*t15*t158*t72+7.327274635867936d-3*t
     >         1*t12*t15*t42*t43-6.58380164251622d-5*gammaaa*t1*t15*t158
     ?         *t39*t42+5.901116393442186d-7*t1*t153*t157*t38*t42+4.0341
     @         12095367414d-2*t1*t12*t41-1.0210661452298034d-3*gammaaa*t
     1         1*t39*t40+t160-3.8035d-1*t1*t159)*wght)+1.919347107159836
     2         7d+2*(t162+t161)*t20*t34+t164+t156+t155+Cmat3(iq,D3_RA_RA
     3         _GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = 1.385405033612855d+1*t19*t8*(-1.667
     1         2336067107812d-1*rhoa*t98*wght-8.336168033553906d-2*t11*(
     2         7.388434626882859d-5*gammaaa*t1*t40*t42*t72-5.90111639344
     3         2186d-7*t1*t153*t157*t38*t42-5.681848021644245d-4*t1*t12*
     4         t41+t190+t160)*wght)+1.9193471071598367d+2*t20*t34*(5.001
     5         700820132344d-1*t10*t11*t33*t97*wght+t161)+5.156620156177
     6         406d+1*t17*t30*(6.668934426843125d-1*t10*t11*t16*t21*t97*
     7         wght+t163)+1.385405033612855d+1*t20*t8*(8.336168033553906
     8         d-2*t10*t11*(5.479893277885249d-5*gammaaa*t1*t40*t42*t72+
     9         6.44577519522176d+0*t1*t192*t4*t5-1.0842021724855044d-19*
     :         t1*t12*t15*t42*t43-4.3767712227948996d-7*t1*t153*t157*t38
     ;         *t42-4.2141430968998067d-4*t1*t12*t41+2.1485917317405864d
     <         +0*t1*t101*t14*t4+t195+t194+t190+t181+t180+t179+t178-1.73
     =         17562920160695d+0*t1*t101*t15-6.208750505646785d-2*gammaa
     >         a*t1*t101*t12*t13)*wght+8.336168033553906d-2*t103*t11*t22
     ?         *wght+8.336168033553906d-2*rhoa*t10*t103*wght-t21*t99-t52
     @         *t97+t197+t196+t174+t173)+5.156620156177406d+1*t17*t19*(1
     1         .6672336067107812d-1*t101*t11*t22*wght+1.6672336067107812
     2         d-1*t10*t11*t192*wght+1.6672336067107812d-1*rhoa*t10*t101
     3         *wght-2.0d+0*t16*t99+t193+t191+t172+t166+t165)+1.38540503
     4         3612855d+1*t30*t8*(-3.3344672134215625d-1*t10*t103*t11*t2
     5         1*wght-2*t71*t97+t189)+5.156620156177406d+1*t17*t20*(-3.3
     6         344672134215625d-1*t10*t101*t11*t21*wght-3.33446721342156
     7         25d-1*t10*t103*t11*t16*wght-t60*t97+t188+t187)+1.91934710
     8         71598367d+2*t19*t34*(-1.0003401640264689d+0*t10*t101*t11*
     9         t16*wght+t186+t185)+6.929381492836896d+0*t10*t11*t117*t31
     :         *t8*t97*wght+t155+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 1.385405033612855d+1*t19*t8*(-8.336
     1         168033553906d-2*rhoa*t98*wght-8.336168033553906d-2*rhoa*t
     2         84*wght-8.336168033553906d-2*t11*(-1.2678333333333333d-1*
     3         t1*t47*t83+3.663637317933968d-3*gammaaa*t1*t12*t42*t43*t8
     4         3+2.0437371323672143d-1*t1*t15*t24*t42*t83+1.313499489223
     5         619d-4*gammaaa*t1*t40*t42*t72-1.3233441301457602d-5*gamma
     6         aa*t1*t127*t15*t158*t72+3.663637317933968d-3*t1*t12*t15*t
     7         42*t43-3.29190082125811d-5*gammaaa*t1*t15*t158*t39*t42-1.
     8         1802232786884372d-6*t1*t153*t157*t38*t42-4.26138601623318
     9         46d-3*t1*t12*t41+4.0842645809192135d-5*gammaaa*t1*t39*t40
     :         -2.3722487901637587d-7*t1*t127*t157*t198*t38+t160)*wght-2
     ;         .1314177932294023d-4*t1*t12*t25*wght)+1.385405033612855d+
     <         1*t30*t8*(-1.6672336067107812d-1*t10*t11*t21*t95*wght-1.6
     =         672336067107812d-1*t10*t103*t11*t21*wght-2*t29*t78+t189)+
     >         5.156620156177406d+1*t17*t20*(-1.6672336067107812d-1*t10*
     ?         t11*t16*t95*wght-1.6672336067107812d-1*t10*t11*t21*t89*wg
     @         ht-1.6672336067107812d-1*t10*t101*t11*t21*wght-1.66723360
     1         67107812d-1*t10*t103*t11*t16*wght-9.704472655028322d-4*t1
     2         *t12*t25*t78-t29*t76+t187)+1.9193471071598367d+2*t19*t34*
     3         (-5.001700820132344d-1*t10*t11*t16*t89*wght-5.00170082013
     4         2344d-1*t10*t101*t11*t16*wght-1.4556708982542485d-3*t1*t1
     5         2*t25*t76+t185)+1.385405033612855d+1*t20*t8*(8.3361680335
     6         53906d-2*t10*t11*(-9.403333333333332d-2*t1*t47*t83+2.7172
     7         659061106147d-3*gammaaa*t1*t12*t42*t43*t83+1.515809767426
     8         8203d-1*t1*t15*t24*t42*t83+1.1641407198087722d-2*t1*t12*t
     9         25*t75+9.742032494018218d-5*gammaaa*t1*t40*t42*t72-9.8150
     :         48747577728d-6*gammaaa*t1*t127*t15*t158*t72+6.44577519522
     ;         176d+0*t1*t199*t4*t5+2.7172659061106147d-3*t1*t12*t15*t42
     <         *t43-2.4415544148203305d-5*gammaaa*t1*t15*t158*t39*t42-8.
     =         7535424455898d-7*t1*t153*t157*t38*t42-3.1606073226748554d
     >         -3*t1*t12*t41+3.029238959582779d-5*gammaaa*t1*t39*t40-1.7
     ?         5946203156355d-7*t1*t127*t157*t198*t38+t195+t194+t183+t18
     @         2+t180+t178+t177)*wght+2.1314177932294023d-4*t1*t12*t24*t
     1         77*wght-t29*t74+t197+t196+t184+t176+t175+t173)+5.15662015
     2         6177406d+1*t17*t19*(4.2628355864588047d-4*t1*t12*t24*t75*
     3         wght+1.6672336067107812d-1*t10*t11*t199*wght-9.7044726550
     4         28322d-4*t1*t12*t25*t74+t193+t191+t171+t168+t167+t165)+1.
     5         9193471071598367d+2*t20*t34*(4.8538868837607063d-4*t1*t10
     6         *t12*t16*t21*t24*wght+t162)+t164+t156+t155+Cmat3(iq,D3_RA
     7         _RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.385405033612855d+1*t30*t8*(-1.66
     1         72336067107812d-1*t10*t11*t29*t95*wght-2*t29*t96+t225+t22
     2         4)+5.156620156177406d+1*t17*t20*(-8.08981147293451d-5*t1*
     3         t10*t12*t24*t95*wght-1.6672336067107812d-1*t10*t11*t29*t8
     4         9*wght-t29*t91+t223+t222+t221+t220+t219)+1.38540503361285
     5         5d+1*t20*t8*(4.2628355864588047d-4*t1*t12*t24*t95*wght+8.
     6         336168033553906d-2*t10*t11*(2.3282814396175444d-2*t1*t12*
     7         t25*t89-4.336808689942018d-19*t1*t12*t2*t42*t83+3.4249332
     8         986782793d-5*t1*t42*t72*t82+1.7607451452574904d-4*t1*t39*
     9         t82-1.2248426929463342d-6*gammaaa*t1*t202*t72+6.445775195
     :         22176d+0*t1*t207*t4*t5-4.923867625644263d-7*gammaaa*t1*t1
     ;         53*t202*t42+t211+t210+t209-2.821d-1*t1*t203-6.59798261836
     <         3312d-8*gammaaa*t1*t127*t198*t202)*wght-2*t29*t85+t213+t2
     =         12)+1.9193471071598367d+2*t19*t34*(-2.4269434418803534d-4
     >         *t1*t10*t12*t24*t89*wght+t218+t217+t216)+5.15662015617740
     ?         6d+1*t17*t19*(8.525671172917609d-4*t1*t12*t24*t89*wght+1.
     @         6672336067107812d-1*t10*t11*t207*wght-1.9408945310056644d
     1         -3*t1*t12*t25*t85+t208+t205)+1.385405033612855d+1*t19*t8*
     2         (t201-8.336168033553906d-2*t11*(4.6177716418017856d-5*t1*
     3         t42*t72*t82+2.3739787876592927d-4*t1*t39*t82-1.6514318265
     4         23d-6*gammaaa*t1*t202*t72-6.63875594262246d-7*gammaaa*t1*
     5         t153*t202*t42+t204-3.8035d-1*t1*t203-8.895932963114095d-8
     6         *gammaaa*t1*t127*t198*t202)*wght)+1.9193471071598367d+2*t
     7         20*t34*(2.4269434418803534d-4*t1*t10*t12*t16*t24*t29*wght
     8         +t214)+6.929381492836896d+0*t10*t107*t11*t117*t21*t8*wght
     9         +5.156620156177406d+1*t17*t30*(t215-2*t29*t87)+t200+Cmat3
     :         (iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = 1.385405033612855d+1*t20*t8*(2.131
     1         4177932294023d-4*t1*t12*t24*t95*wght+8.336168033553906d-2
     2         *t10*t11*(1.1641407198087722d-2*t1*t12*t25*t89-2.16840434
     3         4971009d-19*t1*t12*t2*t42*t83-1.598302206049864d-5*t1*t42
     4         *t72*t82+4.2049980784980245d-5*t1*t105*t16*t42*t72-1.8403
     5         21640170825d-6*t1*t126*t127*t15*t72+6.44577519522176d+0*t
     6         1*t227*t4*t5+1.1641407198087722d-2*t1*t101*t12*t25+t228+1
     7         .3234889800848443d-23*gammaaa*t1*t127*t198*t202)*wght+7.6
     8         98904069327958d-7*t1*t114*t21*t42*t72*wght+2.131417793229
     9         4023d-4*t1*t103*t12*t24*wght-t29*t99-t85*t97-t115*t23)+5.
     :         156620156177406d+1*t17*t19*(4.2628355864588047d-4*t1*t12*
     ;         t24*t89*wght+1.5397808138655913d-6*t1*t114*t16*t42*t72*wg
     <         ht+4.2628355864588047d-4*t1*t101*t12*t24*wght+1.667233606
     =         7107812d-1*t10*t11*t227*wght-9.704472655028322d-4*t1*t12*
     >         t25*t99-9.704472655028322d-4*t1*t12*t25*t85-3.50535705631
     ?         30826d-6*t1*t105*t23*t42*t72)+1.385405033612855d+1*t19*t8
     @         *(-8.336168033553906d-2*t11*(-2.1549600995075005d-5*t1*t4
     1         2*t72*t82-2.4812702440233023d-6*t1*t126*t127*t15*t72+t226
     2         )*wght-7.698904069327958d-7*t1*t128*t42*t72*wght)+5.15662
     3         0156177406d+1*t17*t20*(-2.922124543902974d-7*t1*t10*t114*
     4         t21*t42*t72*wght-1.6672336067107812d-1*t10*t101*t11*t29*w
     5         ght-1.6672336067107812d-1*t10*t11*t115*t16*wght-t91*t97+t
     6         231+t223+t220+t219)+1.9193471071598367d+2*t19*t34*(-8.766
     7         373631708922d-7*t1*t10*t114*t16*t42*t72*wght+t230+t218+t2
     8         16)+1.385405033612855d+1*t30*t8*(-1.6672336067107812d-1*t
     9         10*t103*t11*t29*wght-1.6672336067107812d-1*t10*t11*t115*t
     :         21*wght-2*t96*t97+t224)+6.929381492836896d+0*t10*t11*t117
     ;         *t21*t29*t8*t97*wght+5.156620156177406d+1*t17*t30*(t215-2
     <         *t87*t97)+1.9193471071598367d+2*t20*(t229+t214)*t34+t200+
     =         Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = 5.156620156177406d+1*t17*t30*(1.61
     1         7962294586902d-4*t1*t10*t12*t21*t24*t97*wght-2*t100*t97)+
     2         1.385405033612855d+1*t30*t8*(-4.2628355864588047d-4*t1*t1
     3         2*t21*t24*t97*wght-1.6672336067107812d-1*t10*t103*t11*t97
     4         *wght-2*t104*t97+t225)+5.156620156177406d+1*t17*t20*(-4.2
     5         628355864588047d-4*t1*t12*t16*t24*t97*wght-1.667233606710
     6         7812d-1*t10*t101*t11*t97*wght-t102*t97-9.704472655028322d
     7         -4*t1*t104*t12*t25+t231+t222+t221+t219)+1.385405033612855
     8         d+1*t20*t8*(8.336168033553906d-2*t10*t11*(-2.283288865785
     9         5197d-6*t1*t42*t72*t82-5.67982304921771d-6*t1*t39*t82+6.4
     :         4577519522176d+0*t1*t232*t4*t5+2.3282814396175444d-2*t1*t
     ;         101*t12*t25+t228+t211+t210+t209+6.597982618363315d-8*gamm
     <         aaa*t1*t127*t198*t202)*wght+4.2628355864588047d-4*t1*t103
     =         *t12*t24*wght-2*t97*t99+t213+t212)+1.385405033612855d+1*t
     >         19*t8*(t201-8.336168033553906d-2*t11*(-3.078514427867857d
     ?         -6*t1*t42*t72*t82-7.657996089223525d-6*t1*t39*t82+t226+t2
     @         04+8.895932963114095d-8*gammaaa*t1*t127*t198*t202)*wght)+
     1         5.156620156177406d+1*t17*t19*(8.525671172917609d-4*t1*t10
     2         1*t12*t24*wght+1.6672336067107812d-1*t10*t11*t232*wght-1.
     3         9408945310056644d-3*t1*t12*t25*t99+t208+t205)+6.929381492
     4         836896d+0*t10*t11*t117*t21*t8*t97**2*wght+1.9193471071598
     5         367d+2*t19*(-1.4556708982542485d-3*t1*t102*t12*t25+t230+t
     6         217+t216)*t34+1.9193471071598367d+2*t20*(t229-1.455670898
     7         2542485d-3*t1*t100*t12*t25)*t34+t200+Cmat3(iq,D3_RA_GBB_G
     8         BB)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.385405033612855d+1*t20*t8*(8.33
     1         6168033553906d-2*t10*t11*(2.2965800492743765d-7*t1*t235*t
     2         72+6.44577519522176d+0*t1*t236*t4*t5-2.0403644433615356d-
     3         8*t125*t198*t27*t42+1.8464503596165988d-7*t1*t153*t235*t4
     4         2-5.075533441197848d-8*t125*t153*t27+2.3282814396175444d-
     5         2*t1*t109*t12*t25+3.711365222829364d-8*t1*t127*t198*t235)
     6         *wght+2.500850410066172d-1*t106*t11*t29*wght+6.3942533796
     7         88208d-4*t1*t112*t12*t24*wght)+1.385405033612855d+1*t30*t
     8         8*(-3.3344672134215625d-1*t10*t11*t112*t29*wght-2*t113*t2
     9         9+t240)+5.156620156177406d+1*t17*t20*(-3.3344672134215625
     :         d-1*t10*t109*t11*t29*wght-1.617962294586902d-4*t1*t10*t11
     ;         2*t12*t24*wght-t110*t29+t239+t238)+1.9193471071598367d+2*
     <         t19*t34*(-4.853886883760707d-4*t1*t10*t109*t12*t24*wght-3
     =         .010951273132608d-10*t105*t198*t233*wght+t237)+5.15662015
     >         6177406d+1*t17*t19*(1.2788506759376417d-3*t1*t109*t12*t24
     ?         *wght+2.4269434418803534d-4*t1*t106*t12*t24*wght+1.667233
     @         6067107812d-1*t10*t11*t236*wght)-1.1548969154728159d+0*t1
     1         1*t19*(3.096434674730625d-7*t1*t235*t72+2.489533478483423
     2         3d-7*t1*t153*t235*t42+5.003962291751678d-8*t1*t127*t198*t
     3         235)*t8*wght+6.929381492836896d+0*t10*t11*t117*t29**3*t8*
     4         wght+6.780728844877855d-5*t10*t2*t20*t27*t29*t34*t72*wght
     5         +2.50296509406056d-2*t1*t10*t107*t12*t17*t24*t30*wght+t23
     6         4+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = 1.9193471071598367d+2*t20*t34*(1.
     1         1776103133514104d-7*t10*t2*t27*t72*t97*wght+2.35522062670
     2         28208d-7*t10*t2*t27*t29*t72*wght)+5.156620156177406d+1*t1
     3         7*t30*(3.235924589173804d-4*t1*t10*t12*t24*t29*t97*wght+1
     4         .617962294586902d-4*t1*t10*t107*t12*t24*wght)+1.385405033
     5         612855d+1*t20*t8*(8.336168033553906d-2*t106*t11*t97*wght+
     6         1.5397808138655913d-6*t1*t114*t29*t42*t72*wght+8.33616803
     7         3553906d-2*t10*t11*(6.44577519522176d+0*t1*t241*t4*t5+4.0
     8         80728886723071d-8*t125*t198*t27*t42-6.154834532055332d-8*
     9         t1*t153*t235*t42+1.1641407198087722d-2*t1*t109*t12*t25-1.
     :         2371217409431227d-8*t1*t127*t198*t235)*wght+4.26283558645
     ;         88047d-4*t1*t115*t12*t24*wght+2.1314177932294023d-4*t1*t1
     <         12*t12*t24*wght)+5.156620156177406d+1*t17*t20*(-5.8442490
     =         87805948d-7*t1*t10*t114*t29*t42*t72*wght-1.61796229458690
     >         2d-4*t1*t10*t115*t12*t24*wght-t110*t97+t239+t238)+1.91934
     ?         71071598367d+2*t19*t34*(-8.507303319268057d-10*t10*t126*t
     @         198*t27*t42*wght-3.0109512731326077d-10*t105*t198*t233*wg
     1         ht+t237)+5.156620156177406d+1*t17*t19*(2.988552160579178d
     2         -9*t126*t198*t27*t42*wght+1.6672336067107812d-1*t10*t11*t
     3         241*wght+4.2628355864588047d-4*t1*t109*t12*t24*wght+8.089
     4         81147293451d-5*t1*t106*t12*t24*wght)+1.385405033612855d+1
     5         *t30*t8*(-3.3344672134215625d-1*t10*t11*t115*t29*wght-2*t
     6         113*t97+t240)+6.929381492836896d+0*t10*t107*t11*t117*t8*t
     7         97*wght-1.1548969154728159d+0*t11*t19*(-8.298444928278073
     8         d-8*t1*t153*t235*t42-1.6679874305838935d-8*t1*t127*t198*t
     9         235)*t8*wght+t234+Cmat3(iq,D3_GAA_GAA_GBB)
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
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.492667444378077d-3*gammaaa*t2+1.0d+0
            t4 = 1.804d+0-8.04d-1/t3
            t5 = rhoa**3.333333333333333d-1
            t6 = 1/rhob**2.6666666666666666d+0
            t7 = 4.492667444378077d-3*gammabb*t6+1.0d+0
            t8 = 1.804d+0-8.04d-1/t7
            t9 = rhob**3.333333333333333d-1
            t10 = 1.8610514726982d+0*t8*t9+1.8610514726982d+0*t4*t5
            t11 = 1/t10
            t12 = t1*t11*t4*t5*t8*t9
            t13 = t12**2.0d+0
            t14 = 1/t13
            t15 = 3.908227599821866d+0*t1*t11*t4*t5*t8*t9+1.199591941974
     1         8984d+1*t13+3.183d-1
            t16 = 1/t15
            t17 = 5.269388045346497d+0*t1*t11*t4*t5*t8*t9+5.764d-1
            t18 = 1/t12**3.0d+0
            t19 = 1/t3**2
            t20 = 1/rhoa**3.3333333333333337d+0
            t21 = 1/rhoa**6.666666666666666d-1
            t22 = 6.203504908994d-1*t21*t4-1.7926167019779396d-2*gammaaa
     1         *t19*t20
            t23 = 1/t10**2
            t24 = -t1*t22*t23*t4*t5*t8*t9+3.333333333333333d-1*t1*t11*t2
     1         1*t4*t8*t9-9.632279000746596d-3*gammaaa*t1*t11*t19*t20*t8
     2         *t9
            t25 = 1/t15**2
            t26 = 2.3991838839497967d+1*t1*t11*t24*t4*t5*t8*t9-3.9082275
     1         99821866d+0*t1*t22*t23*t4*t5*t8*t9+1.302742533273955d+0*t
     2         1*t11*t21*t4*t8*t9-3.7645138639902426d-2*gammaaa*t1*t11*t
     3         19*t20*t8*t9
            t27 = -5.269388045346497d+0*t1*t22*t23*t4*t5*t8*t9+1.7564626
     1         817821655d+0*t1*t11*t21*t4*t8*t9-5.075621581597621d-2*gam
     2         maaa*t1*t11*t19*t20*t8*t9
            t28 = -8.336168033553906d-2*rhoa*rhob*t27*wght-8.33616803355
     1         3906d-2*rhob*t17*wght
            t29 = 1/t7**2
            t30 = 1/rhob**3.3333333333333337d+0
            t31 = 1/rhob**6.666666666666666d-1
            t32 = 6.203504908994d-1*t31*t8-1.7926167019779396d-2*gammabb
     1         *t29*t30
            t33 = -t1*t23*t32*t4*t5*t8*t9+3.333333333333333d-1*t1*t11*t3
     1         1*t4*t5*t8-9.632279000746596d-3*gammabb*t1*t11*t29*t30*t4
     2         *t5
            t34 = 2.3991838839497967d+1*t1*t11*t33*t4*t5*t8*t9-3.9082275
     1         99821866d+0*t1*t23*t32*t4*t5*t8*t9+1.302742533273955d+0*t
     2         1*t11*t31*t4*t5*t8-3.7645138639902426d-2*gammabb*t1*t11*t
     3         29*t30*t4*t5
            t35 = -5.269388045346497d+0*t1*t23*t32*t4*t5*t8*t9+1.7564626
     1         817821655d+0*t1*t11*t31*t4*t5*t8-5.075621581597621d-2*gam
     2         mabb*t1*t11*t29*t30*t4*t5
            t36 = -8.336168033553906d-2*rhoa*rhob*t35*wght-8.33616803355
     1         3906d-2*rhoa*t17*wght
            t37 = 1/rhoa**2
            t38 = 1/rhoa**2.3333333333333334d+0
            t39 = 1.9033580930991081d-2*t1*t11*t19*t38*t8*t9-3.542247382
     1         2341327d-2*t1*t19*t23*t37*t4*t8*t9
            t40 = 3.612104625279974d-3*t1*t11*t19*t38*t8*t9-6.7223126324
     1         17276d-3*t1*t19*t23*t37*t4*t8*t9
            t41 = 2.3991838839497967d+1*t1*t11*t4*t40*t5*t8*t9-2.6272327
     1         764644376d-2*t1*t19*t23*t37*t4*t8*t9+1.4116926989963413d-
     2         2*t1*t11*t19*t38*t8*t9
            t42 = 1/rhob**2.3333333333333334d+0
            t43 = 1/rhob**2
            t44 = 1.9033580930991081d-2*t1*t11*t29*t4*t42*t5-3.542247382
     1         2341327d-2*t1*t23*t29*t4*t43*t5*t8
            t45 = 3.612104625279974d-3*t1*t11*t29*t4*t42*t5-6.7223126324
     1         17276d-3*t1*t23*t29*t4*t43*t5*t8
            t46 = 2.3991838839497967d+1*t1*t11*t4*t45*t5*t8*t9-2.6272327
     1         764644376d-2*t1*t23*t29*t4*t43*t5*t8+1.4116926989963413d-
     2         2*t1*t11*t29*t4*t42*t5
            t47 = 1/t15**3
            t48 = t26**2
            t49 = 1/t12**4.0d+0
            t50 = t24**2
            t51 = t22**2
            t52 = 1/t10**3
            t53 = gammaaa**2
            t54 = 1/t3**3
            t55 = 1/rhoa**7
            t56 = 1/rhoa**4.333333333333333d+0
            t57 = 1/rhoa**1.6666666666666669d+0
            t58 = -4.135669939329333d-1*t4*t57+5.377850105933819d-2*gamm
     1         aaa*t19*t56-4.295269705186499d-4*t53*t54*t55
            t59 = -5.269388045346497d+0*t1*t23*t4*t5*t58*t8*t9-1.1709751
     1         2118811d+0*t1*t11*t4*t57*t8*t9+1.5226864744792862d-1*gamm
     2         aaa*t1*t11*t19*t56*t8*t9-1.2161642581134076d-3*t1*t11*t53
     3         *t54*t55*t8*t9+1.0538776090692993d+1*t1*t4*t5*t51*t52*t8*
     4         t9-3.512925363564331d+0*t1*t21*t22*t23*t4*t8*t9+1.0151243
     5         163195242d-1*gammaaa*t1*t19*t20*t22*t23*t8*t9
            t60 = -8.336168033553906d-2*rhoa*rhob*t59*wght-1.66723360671
     1         07812d-1*rhob*t27*wght
            t61 = -t1*t23*t4*t5*t58*t8*t9-2.222222222222222d-1*t1*t11*t4
     1         *t57*t8*t9+2.8896837002239795d-2*gammaaa*t1*t11*t19*t56*t
     2         8*t9-2.3079800683637772d-4*t1*t11*t53*t54*t55*t8*t9+2*t1*
     3         t4*t5*t51*t52*t8*t9-6.666666666666666d-1*t1*t21*t22*t23*t
     4         4*t8*t9+1.9264558001493193d-2*gammaaa*t1*t19*t20*t22*t23*
     5         t8*t9
            t62 = 1.6672336067107812d-1*rhoa*rhob*t17*t61*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t24*t27*wght+1.6672336067107812d-1
     2         *rhob*t17*t24*wght-2.0d+0*t24*t28
            t63 = 2.3991838839497967d+1*t1*t11*t4*t5*t61*t8*t9-3.9082275
     1         99821866d+0*t1*t23*t4*t5*t58*t8*t9-8.684950221826366d-1*t
     2         1*t11*t4*t57*t8*t9+1.1293541591970728d-1*gammaaa*t1*t11*t
     3         19*t56*t8*t9-9.02011140301807d-4*t1*t11*t53*t54*t55*t8*t9
     4         +7.816455199643731d+0*t1*t4*t5*t51*t52*t8*t9-2.3991838839
     5         497967d+1*t1*t22*t23*t24*t4*t5*t8*t9+7.997279613165988d+0
     6         *t1*t11*t21*t24*t4*t8*t9-2.60548506654791d+0*t1*t21*t22*t
     7         23*t4*t8*t9-2.3109608544299284d-1*gammaaa*t1*t11*t19*t20*
     8         t24*t8*t9+7.529027727980485d-2*gammaaa*t1*t19*t20*t22*t23
     9         *t8*t9
            t64 = 8.336168033553906d-2*rhoa*rhob*t17*t63*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t26*t27*wght+8.336168033553906d-2*rh
     2         ob*t17*t26*wght-t26*t28
            t65 = 1.0538776090692993d+1*t1*t22*t32*t4*t5*t52*t8*t9-1.756
     1         4626817821655d+0*t1*t21*t23*t32*t4*t8*t9+5.07562158159762
     2         1d-2*gammaaa*t1*t19*t20*t23*t32*t8*t9-1.7564626817821655d
     3         +0*t1*t22*t23*t31*t4*t5*t8+5.85487560594055d-1*t1*t11*t21
     4         *t31*t4*t8-1.6918738605325404d-2*gammaaa*t1*t11*t19*t20*t
     5         31*t8+5.075621581597621d-2*gammabb*t1*t22*t23*t29*t30*t4*
     6         t5-1.6918738605325404d-2*gammabb*t1*t11*t21*t29*t30*t4+4.
     7         8889803176158986d-4*gammaaa*gammabb*t1*t11*t19*t20*t29*t3
     8         0
            t66 = -8.336168033553906d-2*rhoa*rhob*t65*wght-8.33616803355
     1         3906d-2*rhob*t35*wght-8.336168033553906d-2*rhoa*t27*wght-
     2         8.336168033553906d-2*t17*wght
            t67 = -1.6672336067107812d-1*rhoa*rhob*t17*t24*t34*wght-1.66
     1         72336067107812d-1*rhoa*rhob*t17*t26*t33*wght
            t68 = 2*t1*t22*t32*t4*t5*t52*t8*t9-3.333333333333333d-1*t1*t
     1         21*t23*t32*t4*t8*t9+9.632279000746596d-3*gammaaa*t1*t19*t
     2         20*t23*t32*t8*t9-3.333333333333333d-1*t1*t22*t23*t31*t4*t
     3         5*t8+1.111111111111111d-1*t1*t11*t21*t31*t4*t8-3.21075966
     4         6915532d-3*gammaaa*t1*t11*t19*t20*t31*t8+9.63227900074659
     5         6d-3*gammabb*t1*t22*t23*t29*t30*t4*t5-3.210759666915532d-
     6         3*gammabb*t1*t11*t21*t29*t30*t4+9.278079874822385d-5*gamm
     7         aaa*gammabb*t1*t11*t19*t20*t29*t30
            t69 = 1.6672336067107812d-1*rhoa*rhob*t17*t68*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t24*t35*wght+1.6672336067107812d-1
     2         *rhoa*t17*t24*wght-2.0d+0*t28*t33
            t70 = 2.3991838839497967d+1*t1*t11*t4*t5*t68*t8*t9+7.8164551
     1         99643731d+0*t1*t22*t32*t4*t5*t52*t8*t9-2.3991838839497967
     2         d+1*t1*t23*t24*t32*t4*t5*t8*t9-1.302742533273955d+0*t1*t2
     3         1*t23*t32*t4*t8*t9+3.7645138639902426d-2*gammaaa*t1*t19*t
     4         20*t23*t32*t8*t9+7.997279613165988d+0*t1*t11*t24*t31*t4*t
     5         5*t8-1.302742533273955d+0*t1*t22*t23*t31*t4*t5*t8+4.34247
     6         5110913183d-1*t1*t11*t21*t31*t4*t8-1.2548379546634142d-2*
     7         gammaaa*t1*t11*t19*t20*t31*t8-2.3109608544299284d-1*gamma
     8         bb*t1*t11*t24*t29*t30*t4*t5+3.7645138639902426d-2*gammabb
     9         *t1*t22*t23*t29*t30*t4*t5-1.2548379546634142d-2*gammabb*t
     :         1*t11*t21*t29*t30*t4+3.6260847840132643d-4*gammaaa*gammab
     ;         b*t1*t11*t19*t20*t29*t30
            t71 = 8.336168033553906d-2*rhoa*rhob*t17*t70*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t26*t35*wght+8.336168033553906d-2*rh
     2         oa*t17*t26*wght-t28*t34
            t72 = t34**2
            t73 = t33**2
            t74 = gammabb**2
            t75 = 1/t7**3
            t76 = 1/rhob**7
            t77 = 1/rhob**4.333333333333333d+0
            t78 = 1/rhob**1.6666666666666669d+0
            t79 = t32**2
            t80 = -4.135669939329333d-1*t78*t8+5.377850105933819d-2*gamm
     1         abb*t29*t77-4.295269705186499d-4*t74*t75*t76
            t81 = -5.269388045346497d+0*t1*t23*t4*t5*t8*t80*t9+1.0538776
     1         090692993d+1*t1*t4*t5*t52*t79*t8*t9-1.17097512118811d+0*t
     2         1*t11*t4*t5*t78*t8-3.512925363564331d+0*t1*t23*t31*t32*t4
     3         *t5*t8+1.5226864744792862d-1*gammabb*t1*t11*t29*t4*t5*t77
     4         -1.2161642581134076d-3*t1*t11*t4*t5*t74*t75*t76+1.0151243
     5         163195242d-1*gammabb*t1*t23*t29*t30*t32*t4*t5
            t82 = -8.336168033553906d-2*rhoa*rhob*t81*wght-1.66723360671
     1         07812d-1*rhoa*t35*wght
            t83 = -t1*t23*t4*t5*t8*t80*t9+2*t1*t4*t5*t52*t79*t8*t9-2.222
     1         222222222222d-1*t1*t11*t4*t5*t78*t8-6.666666666666666d-1*
     2         t1*t23*t31*t32*t4*t5*t8+2.8896837002239795d-2*gammabb*t1*
     3         t11*t29*t4*t5*t77-2.3079800683637772d-4*t1*t11*t4*t5*t74*
     4         t75*t76+1.9264558001493193d-2*gammabb*t1*t23*t29*t30*t32*
     5         t4*t5
            t84 = 1.6672336067107812d-1*rhoa*rhob*t17*t83*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t33*t35*wght+1.6672336067107812d-1
     2         *rhoa*t17*t33*wght-2.0d+0*t33*t36
            t85 = 2.3991838839497967d+1*t1*t11*t4*t5*t8*t83*t9-3.9082275
     1         99821866d+0*t1*t23*t4*t5*t8*t80*t9+7.816455199643731d+0*t
     2         1*t4*t5*t52*t79*t8*t9-2.3991838839497967d+1*t1*t23*t32*t3
     3         3*t4*t5*t8*t9-8.684950221826366d-1*t1*t11*t4*t5*t78*t8+7.
     4         997279613165988d+0*t1*t11*t31*t33*t4*t5*t8-2.605485066547
     5         91d+0*t1*t23*t31*t32*t4*t5*t8+1.1293541591970728d-1*gamma
     6         bb*t1*t11*t29*t4*t5*t77-9.02011140301807d-4*t1*t11*t4*t5*
     7         t74*t75*t76-2.3109608544299284d-1*gammabb*t1*t11*t29*t30*
     8         t33*t4*t5+7.529027727980485d-2*gammabb*t1*t23*t29*t30*t32
     9         *t4*t5
            t86 = 8.336168033553906d-2*rhoa*rhob*t17*t85*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t34*t35*wght+8.336168033553906d-2*rh
     2         oa*t17*t34*wght-t34*t36
            t87 = 1/t3**4
            t88 = 1/rhoa**5.666666666666667d+0
            t89 = 1/rhoa**3
            t90 = 1/rhoa**6
            t91 = 1.610726139444937d-4*gammaaa*t54*t90-1.568539614230697
     1         3d-2*t19*t20
            t92 = -5.269388045346497d+0*t1*t23*t4*t5*t8*t9*t91+4.5606159
     1         67925279d-4*gammaaa*t1*t11*t54*t8*t9*t90-1.18074912741137
     2         75d-2*t1*t19*t23*t4*t8*t89*t9+3.411991507534344d-4*gammaa
     3         a*t1*t23*t8*t87*t88*t9+7.084494764468265d-2*t1*t19*t22*t3
     4         7*t4*t52*t8*t9-1.9033580930991081d-2*t1*t19*t22*t23*t38*t
     5         8*t9-4.441168883897918d-2*t1*t11*t19*t20*t8*t9
            t93 = -8.336168033553906d-2*rhoa*rhob*t92*wght-8.33616803355
     1         3906d-2*rhob*t39*wght
            t94 = -1.6672336067107812d-1*rhoa*rhob*t17*t24*t41*wght-1.66
     1         72336067107812d-1*rhoa*rhob*t17*t26*t40*wght
            t95 = -t1*t23*t4*t5*t8*t9*t91+8.654925256364164d-5*gammaaa*t
     1         1*t11*t54*t8*t9*t90-2.240770877472425d-3*t1*t19*t23*t4*t8
     2         *t89*t9+6.47511908056865d-5*gammaaa*t1*t23*t8*t87*t88*t9+
     3         1.344462526483455d-2*t1*t19*t22*t37*t4*t52*t8*t9-3.612104
     4         625279974d-3*t1*t19*t22*t23*t38*t8*t9-8.428244125653273d-
     5         3*t1*t11*t19*t20*t8*t9
            t96 = 1.6672336067107812d-1*rhoa*rhob*t17*t95*wght+1.6672336
     1         067107812d-1*rhoa*rhob*t24*t39*wght-2.0d+0*t28*t40
            t97 = 2.3991838839497967d+1*t1*t11*t4*t5*t8*t9*t95-3.9082275
     1         99821866d+0*t1*t23*t4*t5*t8*t9*t91+3.382541776131776d-4*g
     2         ammaaa*t1*t11*t54*t8*t9*t90-8.75744258821479d-3*t1*t19*t2
     3         3*t4*t8*t89*t9+2.5306239102811573d-4*gammaaa*t1*t23*t8*t8
     4         7*t88*t9+5.254465552928875d-2*t1*t19*t22*t37*t4*t52*t8*t9
     5         -1.612806413056766d-1*t1*t19*t23*t24*t37*t4*t8*t9+8.66610
     6         3204112231d-2*t1*t11*t19*t24*t38*t8*t9-1.4116926989963413
     7         d-2*t1*t19*t22*t23*t38*t8*t9-3.293949630991462d-2*t1*t11*
     8         t19*t20*t8*t9
            t98 = 8.336168033553906d-2*rhoa*rhob*t17*t97*wght+8.33616803
     1         3553906d-2*rhoa*rhob*t26*t39*wght-t28*t41
            t99 = 7.084494764468265d-2*t1*t22*t29*t4*t43*t5*t52*t8-1.180
     1         7491274113775d-2*t1*t21*t23*t29*t4*t43*t8+3.4119915075343
     2         44d-4*gammaaa*t1*t19*t20*t23*t29*t43*t8-1.903358093099108
     3         1d-2*t1*t22*t23*t29*t4*t42*t5+6.344526976997026d-3*t1*t11
     4         *t21*t29*t4*t42-1.8333676191059622d-4*gammaaa*t1*t11*t19*
     5         t20*t29*t42
            t100 = -8.336168033553906d-2*rhoa*rhob*t99*wght-8.3361680335
     1         53906d-2*rhob*t44*wght
            t101 = -1.6672336067107812d-1*rhoa*rhob*t17*t24*t46*wght-1.6
     1         672336067107812d-1*rhoa*rhob*t17*t26*t45*wght
            t102 = 1.344462526483455d-2*t1*t22*t29*t4*t43*t5*t52*t8-2.24
     1         0770877472425d-3*t1*t21*t23*t29*t4*t43*t8+6.4751190805686
     2         5d-5*gammaaa*t1*t19*t20*t23*t29*t43*t8-3.612104625279974d
     3         -3*t1*t22*t23*t29*t4*t42*t5+1.2040348750933245d-3*t1*t11*
     4         t21*t29*t4*t42-3.4792799530583945d-5*gammaaa*t1*t11*t19*t
     5         20*t29*t42
            t103 = 1.6672336067107812d-1*rhoa*rhob*t24*t44*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t102*t17*wght-2.0d+0*t28*t45
            t104 = 2.3991838839497967d+1*t1*t102*t11*t4*t5*t8*t9+5.25446
     1         5552928875d-2*t1*t22*t29*t4*t43*t5*t52*t8-1.6128064130567
     2         66d-1*t1*t23*t24*t29*t4*t43*t5*t8-8.75744258821479d-3*t1*
     3         t21*t23*t29*t4*t43*t8+2.5306239102811573d-4*gammaaa*t1*t1
     4         9*t20*t23*t29*t43*t8+8.666103204112231d-2*t1*t11*t24*t29*
     5         t4*t42*t5-1.4116926989963413d-2*t1*t22*t23*t29*t4*t42*t5+
     6         4.705642329987803d-3*t1*t11*t21*t29*t4*t42-1.359781794004
     7         9742d-4*gammaaa*t1*t11*t19*t20*t29*t42
            t105 = 8.336168033553906d-2*rhoa*rhob*t26*t44*wght+8.3361680
     1         33553906d-2*rhoa*rhob*t104*t17*wght-t28*t46
            t106 = 7.084494764468265d-2*t1*t19*t32*t37*t4*t52*t8*t9-1.90
     1         33580930991081d-2*t1*t19*t23*t32*t38*t8*t9-1.180749127411
     2         3775d-2*t1*t19*t23*t31*t37*t4*t8+6.344526976997026d-3*t1*
     3         t11*t19*t31*t38*t8+3.411991507534344d-4*gammabb*t1*t19*t2
     4         3*t29*t30*t37*t4-1.8333676191059622d-4*gammabb*t1*t11*t19
     5         *t29*t30*t38
            t107 = -8.336168033553906d-2*rhoa*t39*wght-8.336168033553906
     1         d-2*rhoa*rhob*t106*wght
            t108 = -1.6672336067107812d-1*rhoa*rhob*t17*t33*t41*wght-1.6
     1         672336067107812d-1*rhoa*rhob*t17*t34*t40*wght
            t109 = 1.344462526483455d-2*t1*t19*t32*t37*t4*t52*t8*t9-3.61
     1         2104625279974d-3*t1*t19*t23*t32*t38*t8*t9-2.2407708774724
     2         25d-3*t1*t19*t23*t31*t37*t4*t8+1.2040348750933245d-3*t1*t
     3         11*t19*t31*t38*t8+6.47511908056865d-5*gammabb*t1*t19*t23*
     4         t29*t30*t37*t4-3.4792799530583945d-5*gammabb*t1*t11*t19*t
     5         29*t30*t38
            t110 = 1.6672336067107812d-1*rhoa*rhob*t33*t39*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t109*t17*wght-2.0d+0*t36*t40
            t111 = 5.254465552928875d-2*t1*t19*t32*t37*t4*t52*t8*t9+2.39
     1         91838839497967d+1*t1*t109*t11*t4*t5*t8*t9-1.6128064130567
     2         66d-1*t1*t19*t23*t33*t37*t4*t8*t9+8.666103204112231d-2*t1
     3         *t11*t19*t33*t38*t8*t9-1.4116926989963413d-2*t1*t19*t23*t
     4         32*t38*t8*t9-8.75744258821479d-3*t1*t19*t23*t31*t37*t4*t8
     5         +4.705642329987803d-3*t1*t11*t19*t31*t38*t8+2.53062391028
     6         11573d-4*gammabb*t1*t19*t23*t29*t30*t37*t4-1.359781794004
     7         9742d-4*gammabb*t1*t11*t19*t29*t30*t38
            t112 = 8.336168033553906d-2*rhoa*rhob*t34*t39*wght+8.3361680
     1         33553906d-2*rhoa*rhob*t111*t17*wght-t36*t41
            t113 = 1/rhob**6
            t114 = 1/t7**4
            t115 = 1/rhob**5.666666666666667d+0
            t116 = 1/rhob**3
            t117 = 1.610726139444937d-4*gammabb*t113*t75-1.5685396142306
     1         973d-2*t29*t30
            t118 = -5.269388045346497d+0*t1*t117*t23*t4*t5*t8*t9+7.08449
     1         4764468265d-2*t1*t29*t32*t4*t43*t5*t52*t8-1.1807491274113
     2         775d-2*t1*t116*t23*t29*t4*t5*t8+4.560615967925279d-4*gamm
     3         abb*t1*t11*t113*t4*t5*t75-1.9033580930991081d-2*t1*t23*t2
     4         9*t32*t4*t42*t5-4.441168883897918d-2*t1*t11*t29*t30*t4*t5
     5         +3.411991507534344d-4*gammabb*t1*t114*t115*t23*t4*t5
            t119 = -8.336168033553906d-2*rhoa*t44*wght-8.336168033553906
     1         d-2*rhoa*rhob*t118*wght
            t120 = -1.6672336067107812d-1*rhoa*rhob*t17*t33*t46*wght-1.6
     1         672336067107812d-1*rhoa*rhob*t17*t34*t45*wght
            t121 = -t1*t117*t23*t4*t5*t8*t9+1.344462526483455d-2*t1*t29*
     1         t32*t4*t43*t5*t52*t8-2.240770877472425d-3*t1*t116*t23*t29
     2         *t4*t5*t8+8.654925256364164d-5*gammabb*t1*t11*t113*t4*t5*
     3         t75-3.612104625279974d-3*t1*t23*t29*t32*t4*t42*t5-8.42824
     4         4125653273d-3*t1*t11*t29*t30*t4*t5+6.47511908056865d-5*ga
     5         mmabb*t1*t114*t115*t23*t4*t5
            t122 = 1.6672336067107812d-1*rhoa*rhob*t33*t44*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t121*t17*wght-2.0d+0*t36*t45
            t123 = -3.908227599821866d+0*t1*t117*t23*t4*t5*t8*t9+2.39918
     1         38839497967d+1*t1*t11*t121*t4*t5*t8*t9+5.254465552928875d
     2         -2*t1*t29*t32*t4*t43*t5*t52*t8-1.612806413056766d-1*t1*t2
     3         3*t29*t33*t4*t43*t5*t8-8.75744258821479d-3*t1*t116*t23*t2
     4         9*t4*t5*t8+3.382541776131776d-4*gammabb*t1*t11*t113*t4*t5
     5         *t75+8.666103204112231d-2*t1*t11*t29*t33*t4*t42*t5-1.4116
     6         926989963413d-2*t1*t23*t29*t32*t4*t42*t5-3.29394963099146
     7         2d-2*t1*t11*t29*t30*t4*t5+2.5306239102811573d-4*gammabb*t
     8         1*t114*t115*t23*t4*t5
            t124 = 8.336168033553906d-2*rhoa*rhob*t34*t44*wght+8.3361680
     1         33553906d-2*rhoa*rhob*t123*t17*wght-t36*t46
            t125 = 1/rhoa**4.666666666666667d+0
            t126 = 1/rhoa**5
            t127 = 4.762418864947907d-4*t1*t4*t52*t56*t8*t87*t9-2.558993
     1         630650758d-4*t1*t125*t23*t8*t87*t9+3.1828278988193504d-4*
     2         t1*t125*t23*t4*t54*t8*t9-1.71023098797198d-4*t1*t11*t126*
     3         t54*t8*t9
            t128 = t41**2
            t129 = t40**2
            t130 = 9.037897425591375d-5*t1*t4*t52*t56*t8*t87*t9-4.856339
     1         310426487d-5*t1*t125*t23*t8*t87*t9+6.040223022918516d-5*t
     2         1*t125*t23*t4*t54*t8*t9-3.245596971136562d-5*t1*t11*t126*
     3         t54*t8*t9
            t131 = 3.3344672134215625d-1*rhoa*rhob*t39*t40*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t130*t17*wght
            t132 = 3.53221601630552d-4*t1*t4*t52*t56*t8*t87*t9-1.8979679
     1         327108684d-4*t1*t125*t23*t8*t87*t9+2.3606566327249603d-4*
     2         t1*t125*t23*t4*t54*t8*t9-1.268453166049416d-4*t1*t11*t126
     3         *t54*t8*t9+2.3991838839497967d+1*t1*t11*t130*t4*t5*t8*t9-
     4         1.612806413056766d-1*t1*t19*t23*t37*t4*t40*t8*t9+8.666103
     5         204112231d-2*t1*t11*t19*t38*t40*t8*t9
            t133 = 1.6672336067107812d-1*rhoa*rhob*t39*t41*wght+8.336168
     1         033553906d-2*rhoa*rhob*t132*t17*wght
            t134 = 4.762418864947907d-4*t1*t19*t29*t37*t4*t43*t52*t8-1.2
     1         79496815325379d-4*t1*t19*t23*t29*t38*t43*t8-1.27949681532
     2         5379d-4*t1*t19*t23*t29*t37*t4*t42+6.875128571647359d-5*t1
     3         *t11*t19*t29*t38*t42
            t135 = 9.037897425591375d-5*t1*t19*t29*t37*t4*t43*t52*t8-2.4
     1         281696552132434d-5*t1*t19*t23*t29*t38*t43*t8-2.4281696552
     2         132434d-5*t1*t19*t23*t29*t37*t4*t42+1.304729982396898d-5*
     3         t1*t11*t19*t29*t38*t42
            t136 = 1.6672336067107812d-1*rhoa*rhob*t39*t45*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t40*t44*wght+1.6672336067107812d-
     2         1*rhoa*rhob*t135*t17*wght
            t137 = 2.3991838839497967d+1*t1*t11*t135*t4*t5*t8*t9+3.53221
     1         601630552d-4*t1*t19*t29*t37*t4*t43*t52*t8-1.6128064130567
     2         66d-1*t1*t23*t29*t4*t40*t43*t5*t8-9.489839663554341d-5*t1
     3         *t19*t23*t29*t38*t43*t8+8.666103204112231d-2*t1*t11*t29*t
     4         4*t40*t42*t5-9.489839663554341d-5*t1*t19*t23*t29*t37*t4*t
     5         42+5.099181727518654d-5*t1*t11*t19*t29*t38*t42
            t138 = 8.336168033553906d-2*rhoa*rhob*t39*t46*wght+8.3361680
     1         33553906d-2*rhoa*rhob*t41*t44*wght+8.336168033553906d-2*r
     2         hoa*rhob*t137*t17*wght
            t139 = -1.6672336067107812d-1*rhoa*rhob*t17*t40*t46*wght-1.6
     1         672336067107812d-1*rhoa*rhob*t17*t41*t45*wght
            t140 = 1/rhob**5
            t141 = 1/rhob**4.666666666666667d+0
            t142 = 4.762418864947907d-4*t1*t114*t4*t5*t52*t77*t8+3.18282
     1         78988193504d-4*t1*t141*t23*t4*t5*t75*t8-1.71023098797198d
     2         -4*t1*t11*t140*t4*t5*t75-2.558993630650758d-4*t1*t114*t14
     3         1*t23*t4*t5
            t143 = t46**2
            t144 = t45**2
            t145 = 9.037897425591375d-5*t1*t114*t4*t5*t52*t77*t8+6.04022
     1         3022918516d-5*t1*t141*t23*t4*t5*t75*t8-3.245596971136562d
     2         -5*t1*t11*t140*t4*t5*t75-4.856339310426487d-5*t1*t114*t14
     3         1*t23*t4*t5
            t146 = 3.3344672134215625d-1*rhoa*rhob*t44*t45*wght+1.667233
     1         6067107812d-1*rhoa*rhob*t145*t17*wght
            t147 = 2.3991838839497967d+1*t1*t11*t145*t4*t5*t8*t9+3.53221
     1         601630552d-4*t1*t114*t4*t5*t52*t77*t8+2.3606566327249603d
     2         -4*t1*t141*t23*t4*t5*t75*t8-1.612806413056766d-1*t1*t23*t
     3         29*t4*t43*t45*t5*t8-1.268453166049416d-4*t1*t11*t140*t4*t
     4         5*t75+8.666103204112231d-2*t1*t11*t29*t4*t42*t45*t5-1.897
     5         9679327108684d-4*t1*t114*t141*t23*t4*t5
            t148 = 1.6672336067107812d-1*rhoa*rhob*t44*t46*wght+8.336168
     1         033553906d-2*rhoa*rhob*t147*t17*wght
            t149 = 1/t15**4
            t150 = 1/t12**5.0d+0
            t151 = t22**3
            t152 = 1/t10**4
            t153 = gammaaa**3
            t154 = 1/rhoa**1.0666666666666666d+1
            t155 = 1/rhoa**8
            t156 = 1/rhoa**5.333333333333333d+0
            t157 = -1.5437774695451842d-5*t153*t154*t87+4.29526970518649
     1         96d-3*t155*t53*t54+6.892783232215555d-1*t2*t4-2.290565785
     2         860701d-1*gammaaa*t156*t19
            t158 = -8.295189532329142d-6*t1*t11*t153*t154*t8*t87*t9+6*t1
     1         *t22*t4*t5*t52*t58*t8*t9-t1*t21*t23*t4*t58*t8*t9+2.889683
     2         7002239784d-2*gammaaa*t1*t19*t20*t23*t58*t8*t9+6.66666666
     3         6666666d-1*t1*t22*t23*t4*t57*t8*t9-8.669051100671937d-2*g
     4         ammaaa*t1*t19*t22*t23*t56*t8*t9+6.92394020509133d-4*t1*t2
     5         2*t23*t53*t54*t55*t8*t9+2.3079800683637772d-3*t1*t11*t155
     6         *t53*t54*t8*t9+2*t1*t21*t4*t51*t52*t8*t9-5.77936740044795
     7         7d-2*gammaaa*t1*t19*t20*t51*t52*t8*t9-t1*t157*t23*t4*t5*t
     8         8*t9-6*t1*t151*t152*t4*t5*t8*t9+3.7037037037037035d-1*t1*
     9         t11*t2*t4*t8*t9-1.2307912056509542d-1*gammaaa*t1*t11*t156
     :         *t19*t8*t9
            t159 = 2*t1*t32*t4*t5*t52*t58*t8*t9+2.222222222222222d-1*t1*
     1         t23*t32*t4*t57*t8*t9-2.8896837002239795d-2*gammaaa*t1*t19
     2         *t23*t32*t56*t8*t9+2.3079800683637772d-4*t1*t23*t32*t53*t
     3         54*t55*t8*t9+1.3333333333333333d+0*t1*t21*t22*t32*t4*t52*
     4         t8*t9-3.8529116002986386d-2*gammaaa*t1*t19*t20*t22*t32*t5
     5         2*t8*t9-6*t1*t152*t32*t4*t5*t51*t8*t9-3.333333333333333d-
     6         1*t1*t23*t31*t4*t5*t58*t8-7.407407407407407d-2*t1*t11*t31
     7         *t4*t57*t8+9.632279000746597d-3*gammaaa*t1*t11*t19*t31*t5
     8         6*t8-7.693266894545924d-5*t1*t11*t31*t53*t54*t55*t8+6.666
     9         666666666666d-1*t1*t31*t4*t5*t51*t52*t8-2.222222222222222
     :         d-1*t1*t21*t22*t23*t31*t4*t8+6.421519333831064d-3*gammaaa
     ;         *t1*t19*t20*t22*t23*t31*t8+9.632279000746596d-3*gammabb*t
     <         1*t23*t29*t30*t4*t5*t58+2.1405064446103544d-3*gammabb*t1*
     =         t11*t29*t30*t4*t57-2.7834239624467155d-4*gammaaa*gammabb*
     >         t1*t11*t19*t29*t30*t56+2.22311079466421d-6*gammabb*t1*t11
     ?         *t29*t30*t53*t54*t55-1.9264558001493193d-2*gammabb*t1*t29
     @         *t30*t4*t5*t51*t52+6.421519333831064d-3*gammabb*t1*t21*t2
     1         2*t23*t29*t30*t4-1.8556159749644768d-4*gammaaa*gammabb*t1
     2         *t19*t20*t22*t23*t29*t30
            t160 = 2*t1*t22*t4*t5*t52*t8*t80*t9-3.333333333333333d-1*t1*
     1         t21*t23*t4*t8*t80*t9+9.632279000746596d-3*gammaaa*t1*t19*
     2         t20*t23*t8*t80*t9+6.666666666666666d-1*t1*t21*t4*t52*t79*
     3         t8*t9-1.9264558001493193d-2*gammaaa*t1*t19*t20*t52*t79*t8
     4         *t9-6*t1*t152*t22*t4*t5*t79*t8*t9+2.222222222222222d-1*t1
     5         *t22*t23*t4*t5*t78*t8-7.407407407407407d-2*t1*t11*t21*t4*
     6         t78*t8+2.1405064446103544d-3*gammaaa*t1*t11*t19*t20*t78*t
     7         8+1.3333333333333333d+0*t1*t22*t31*t32*t4*t5*t52*t8-2.222
     8         222222222222d-1*t1*t21*t23*t31*t32*t4*t8+6.42151933383106
     9         4d-3*gammaaa*t1*t19*t20*t23*t31*t32*t8-2.8896837002239795
     :         d-2*gammabb*t1*t22*t23*t29*t4*t5*t77+9.632279000746596d-3
     ;         *gammabb*t1*t11*t21*t29*t4*t77-2.7834239624467155d-4*gamm
     <         aaa*gammabb*t1*t11*t19*t20*t29*t77+2.3079800683637772d-4*
     =         t1*t22*t23*t4*t5*t74*t75*t76-7.693266894545924d-5*t1*t11*
     >         t21*t4*t74*t75*t76+2.22311079466421d-6*gammaaa*t1*t11*t19
     ?         *t20*t74*t75*t76-3.8529116002986386d-2*gammabb*t1*t22*t29
     @         *t30*t32*t4*t5*t52+6.421519333831064d-3*gammabb*t1*t21*t2
     1         3*t29*t30*t32*t4-1.8556159749644768d-4*gammaaa*gammabb*t1
     2         *t19*t20*t23*t29*t30*t32
            t161 = gammabb**3
            t162 = 1/rhob**1.0666666666666666d+1
            t163 = 1/rhob**8
            t164 = 1/rhob**5.333333333333333d+0
            t165 = t32**3
            t166 = 6.892783232215555d-1*t6*t8+4.2952697051864996d-3*t163
     1         *t74*t75-2.290565785860701d-1*gammabb*t164*t29-1.54377746
     2         95451842d-5*t114*t161*t162
            t167 = 6*t1*t32*t4*t5*t52*t8*t80*t9-t1*t166*t23*t4*t5*t8*t9-
     1         6*t1*t152*t165*t4*t5*t8*t9-t1*t23*t31*t4*t5*t8*t80+2.8896
     2         837002239784d-2*gammabb*t1*t23*t29*t30*t4*t5*t80+2*t1*t31
     3         *t4*t5*t52*t79*t8+6.666666666666666d-1*t1*t23*t32*t4*t5*t
     4         78*t8+3.7037037037037035d-1*t1*t11*t4*t5*t6*t8-5.77936740
     5         0447957d-2*gammabb*t1*t29*t30*t4*t5*t52*t79-8.66905110067
     6         1937d-2*gammabb*t1*t23*t29*t32*t4*t5*t77+6.92394020509133
     7         d-4*t1*t23*t32*t4*t5*t74*t75*t76+2.3079800683637772d-3*t1
     8         *t11*t163*t4*t5*t74*t75-1.2307912056509542d-1*gammabb*t1*
     9         t11*t164*t29*t4*t5-8.295189532329142d-6*t1*t11*t114*t161*
     :         t162*t4*t5
            t168 = 1/t3**5
            t169 = 1/rhoa**9.333333333333333d+0
            t170 = 1/rhoa**6.666666666666667d+0
            t171 = 1/rhoa**4
            t172 = 1/rhoa**9.666666666666666d+0
            t173 = 5.789165510794442d-6*t172*t53*t87+5.228465380768991d-
     1         2*t19*t56-1.342271782870781d-3*gammaaa*t54*t55
            t174 = 4*t1*t22*t4*t5*t52*t8*t9*t91-6.666666666666666d-1*t1*
     1         t21*t23*t4*t8*t9*t91+1.9264558001493193d-2*gammaaa*t1*t19
     2         *t20*t23*t8*t9*t91-1.7309850512728328d-4*gammaaa*t1*t22*t
     3         23*t54*t8*t9*t90+8.9630835098897d-3*t1*t19*t22*t4*t52*t8*
     4         t89*t9-2.5900476322274596d-4*gammaaa*t1*t22*t52*t8*t87*t8
     5         8*t9+3.1106960746234286d-6*t1*t11*t172*t53*t8*t87*t9-1.94
     6         25357241705948d-4*gammaaa*t1*t170*t23*t8*t87*t9+1.3444625
     7         26483455d-2*t1*t19*t37*t4*t52*t58*t8*t9-3.612104625279974
     8         d-3*t1*t19*t23*t38*t58*t8*t9+2.809414708551091d-2*t1*t11*
     9         t19*t56*t8*t9-7.212437713636803d-4*gammaaa*t1*t11*t54*t55
     :         *t8*t9+1.5514963568929105d-6*t1*t168*t169*t23*t53*t8*t9+7
     ;         .224209250559948d-3*t1*t19*t38*t51*t52*t8*t9-4.0333875794
     <         50365d-2*t1*t152*t19*t37*t4*t51*t8*t9-t1*t173*t23*t4*t5*t
     =         8*t9+1.4938472516482834d-3*t1*t171*t19*t23*t4*t8*t9+1.685
     >         6488251306545d-2*t1*t19*t20*t22*t23*t8*t9
            t175 = 1.344462526483455d-2*t1*t29*t4*t43*t5*t52*t58*t8+1.49
     1         38472516482834d-3*t1*t23*t29*t4*t43*t57*t8-1.942535724170
     2         5948d-4*gammaaa*t1*t19*t23*t29*t43*t56*t8+1.5514963568929
     3         105d-6*t1*t23*t29*t43*t53*t54*t55*t8+8.9630835098897d-3*t
     4         1*t21*t22*t29*t4*t43*t52*t8-2.5900476322274596d-4*gammaaa
     5         *t1*t19*t20*t22*t29*t43*t52*t8-4.033387579450365d-2*t1*t1
     6         52*t29*t4*t43*t5*t51*t8-3.612104625279974d-3*t1*t23*t29*t
     7         4*t42*t5*t58-8.02689916728883d-4*t1*t11*t29*t4*t42*t57+1.
     8         0437839859175184d-4*gammaaa*t1*t11*t19*t29*t42*t56-8.3366
     9         65479990789d-7*t1*t11*t29*t42*t53*t54*t55+7.2242092505599
     :         48d-3*t1*t29*t4*t42*t5*t51*t52-2.408069750186649d-3*t1*t2
     ;         1*t22*t23*t29*t4*t42+6.958559906116789d-5*gammaaa*t1*t19*
     <         t20*t22*t23*t29*t42
            t176 = 2*t1*t32*t4*t5*t52*t8*t9*t91-3.333333333333333d-1*t1*
     1         t23*t31*t4*t5*t8*t91+9.632279000746596d-3*gammabb*t1*t23*
     2         t29*t30*t4*t5*t91-8.654925256364164d-5*gammaaa*t1*t23*t32
     3         *t54*t8*t9*t90+2.884975085454722d-5*gammaaa*t1*t11*t31*t5
     4         4*t8*t90-8.336665479990789d-7*gammaaa*gammabb*t1*t11*t29*
     5         t30*t54*t90+4.48154175494485d-3*t1*t19*t32*t4*t52*t8*t89*
     6         t9-1.2950238161137298d-4*gammaaa*t1*t32*t52*t8*t87*t88*t9
     7         +7.224209250559948d-3*t1*t19*t22*t32*t38*t52*t8*t9-4.0333
     8         87579450365d-2*t1*t152*t19*t22*t32*t37*t4*t8*t9+8.4282441
     9         25653273d-3*t1*t19*t20*t23*t32*t8*t9-7.469236258241416d-4
     :         *t1*t19*t23*t31*t4*t8*t89+2.1583730268562162d-5*gammabb*t
     ;         1*t19*t23*t29*t30*t4*t89+2.1583730268562162d-5*gammaaa*t1
     <         *t23*t31*t8*t87*t88-6.2370153547095d-7*gammaaa*gammabb*t1
     =         *t23*t29*t30*t87*t88+4.48154175494485d-3*t1*t19*t22*t31*t
     >         37*t4*t52*t8-1.2040348750933245d-3*t1*t19*t22*t23*t31*t38
     ?         *t8-2.8094147085510907d-3*t1*t11*t19*t20*t31*t8-1.2950238
     @         161137298d-4*gammabb*t1*t19*t22*t29*t30*t37*t4*t52+3.4792
     1         799530583945d-5*gammabb*t1*t19*t22*t23*t29*t30*t38+8.1183
     2         19890469586d-5*gammabb*t1*t11*t19*t20*t29*t30
            t177 = 2*t1*t117*t22*t4*t5*t52*t8*t9-3.333333333333333d-1*t1
     1         *t117*t21*t23*t4*t8*t9+9.632279000746596d-3*gammaaa*t1*t1
     2         17*t19*t20*t23*t8*t9+4.48154175494485d-3*t1*t116*t22*t29*
     3         t4*t5*t52*t8+4.48154175494485d-3*t1*t21*t29*t32*t4*t43*t5
     4         2*t8-1.2950238161137298d-4*gammaaa*t1*t19*t20*t29*t32*t43
     5         *t52*t8-4.033387579450365d-2*t1*t152*t22*t29*t32*t4*t43*t
     6         5*t8-7.469236258241416d-4*t1*t116*t21*t23*t29*t4*t8+2.158
     7         3730268562162d-5*gammaaa*t1*t116*t19*t20*t23*t29*t8-8.654
     8         925256364164d-5*gammabb*t1*t113*t22*t23*t4*t5*t75+2.88497
     9         5085454722d-5*gammabb*t1*t11*t113*t21*t4*t75-8.3366654799
     :         90789d-7*gammaaa*gammabb*t1*t11*t113*t19*t20*t75+7.224209
     ;         250559948d-3*t1*t22*t29*t32*t4*t42*t5*t52-1.2950238161137
     <         298d-4*gammabb*t1*t114*t115*t22*t4*t5*t52+8.4282441256532
     =         73d-3*t1*t22*t23*t29*t30*t4*t5-1.2040348750933245d-3*t1*t
     >         21*t23*t29*t32*t4*t42+3.4792799530583945d-5*gammaaa*t1*t1
     ?         9*t20*t23*t29*t32*t42-2.8094147085510907d-3*t1*t11*t21*t2
     @         9*t30*t4+2.1583730268562162d-5*gammabb*t1*t114*t115*t21*t
     1         23*t4+8.118319890469586d-5*gammaaa*t1*t11*t19*t20*t29*t30
     2         -6.2370153547095d-7*gammaaa*gammabb*t1*t114*t115*t19*t20*
     3         t23
            t178 = 1.344462526483455d-2*t1*t19*t37*t4*t52*t8*t80*t9-3.61
     1         2104625279974d-3*t1*t19*t23*t38*t8*t80*t9+7.2242092505599
     2         48d-3*t1*t19*t38*t52*t79*t8*t9-4.033387579450365d-2*t1*t1
     3         52*t19*t37*t4*t79*t8*t9+1.4938472516482834d-3*t1*t19*t23*
     4         t37*t4*t78*t8-8.02689916728883d-4*t1*t11*t19*t38*t78*t8+8
     5         .9630835098897d-3*t1*t19*t31*t32*t37*t4*t52*t8-2.40806975
     6         0186649d-3*t1*t19*t23*t31*t32*t38*t8-1.9425357241705948d-
     7         4*gammabb*t1*t19*t23*t29*t37*t4*t77+1.0437839859175184d-4
     8         *gammabb*t1*t11*t19*t29*t38*t77+1.5514963568929105d-6*t1*
     9         t19*t23*t37*t4*t74*t75*t76-8.336665479990789d-7*t1*t11*t1
     :         9*t38*t74*t75*t76-2.5900476322274596d-4*gammabb*t1*t19*t2
     ;         9*t30*t32*t37*t4*t52+6.958559906116789d-5*gammabb*t1*t19*
     <         t23*t29*t30*t32*t38
            t179 = 1/rhob**9.666666666666666d+0
            t180 = 1/t7**5
            t181 = 1/rhob**9.333333333333333d+0
            t182 = 1/rhob**6.666666666666667d+0
            t183 = 1/rhob**4
            t184 = 5.228465380768991d-2*t29*t77-1.342271782870781d-3*gam
     1         mabb*t75*t76+5.789165510794442d-6*t114*t179*t74
            t185 = 4*t1*t117*t32*t4*t5*t52*t8*t9-t1*t184*t23*t4*t5*t8*t9
     1         +1.344462526483455d-2*t1*t29*t4*t43*t5*t52*t8*t80-3.61210
     2         4625279974d-3*t1*t23*t29*t4*t42*t5*t80-4.033387579450365d
     3         -2*t1*t152*t29*t4*t43*t5*t79*t8+8.9630835098897d-3*t1*t11
     4         6*t29*t32*t4*t5*t52*t8-6.666666666666666d-1*t1*t117*t23*t
     5         31*t4*t5*t8+1.4938472516482834d-3*t1*t183*t23*t29*t4*t5*t
     6         8+7.224209250559948d-3*t1*t29*t4*t42*t5*t52*t79+2.8094147
     7         08551091d-2*t1*t11*t29*t4*t5*t77-7.212437713636803d-4*gam
     8         mabb*t1*t11*t4*t5*t75*t76-1.7309850512728328d-4*gammabb*t
     9         1*t113*t23*t32*t4*t5*t75+1.5514963568929105d-6*t1*t180*t1
     :         81*t23*t4*t5*t74+3.1106960746234286d-6*t1*t11*t114*t179*t
     ;         4*t5*t74-2.5900476322274596d-4*gammabb*t1*t114*t115*t32*t
     <         4*t5*t52+1.6856488251306545d-2*t1*t23*t29*t30*t32*t4*t5+1
     =         .9264558001493193d-2*gammabb*t1*t117*t23*t29*t30*t4*t5-1.
     >         9425357241705948d-4*gammabb*t1*t114*t182*t23*t4*t5
            t186 = 1/t3**6
            t187 = 1/rhoa**8.333333333333334d+0
            t188 = 1/rhoa**8.666666666666666d+0
            t189 = 3.0201115114592575d-4*t54*t90-2.1709370665479152d-6*g
     1         ammaaa*t188*t87
            t190 = 2.68892505296691d-2*t1*t19*t37*t4*t52*t8*t9*t91-7.224
     1         209250559948d-3*t1*t19*t23*t38*t8*t9*t91+1.62279848556828
     2         1d-4*t1*t11*t54*t8*t9*t90+1.1331458390995136d-4*t1*t23*t8
     3         *t87*t88*t9+2.0134076743061718d-5*t1*t23*t4*t54*t8*t88*t9
     4         -2.711369227677412d-4*t1*t152*t22*t4*t56*t8*t87*t9+3.0126
     5         32475197125d-5*t1*t156*t4*t52*t8*t87*t9+9.712678620852974
     6         d-5*t1*t125*t22*t52*t8*t87*t9-1.1665110279837858d-6*gamma
     7         aa*t1*t11*t188*t8*t87*t9-1.2080446045837032d-4*t1*t125*t2
     8         2*t4*t52*t54*t8*t9+3.245596971136562d-5*t1*t126*t22*t23*t
     9         54*t8*t9-8.705554958342553d-7*gammaaa*t1*t155*t186*t52*t8
     :         *t9-t1*t189*t23*t4*t5*t8*t9-1.7454334015045245d-6*gammaaa
     ;         *t1*t168*t187*t23*t8*t9
            t191 = 1.344462526483455d-2*t1*t29*t4*t43*t5*t52*t8*t91-3.61
     1         2104625279974d-3*t1*t23*t29*t4*t42*t5*t91-5.8181113383484
     2         14d-7*gammaaa*t1*t23*t29*t43*t54*t8*t90+3.126249554996546
     3         d-7*gammaaa*t1*t11*t29*t42*t54*t90+3.012632475197125d-5*t
     4         1*t19*t29*t4*t43*t52*t8*t89-8.093898850710812d-6*t1*t19*t
     5         23*t29*t4*t42*t89-8.705554958342553d-7*gammaaa*t1*t29*t43
     6         *t52*t8*t87*t88+2.338880758016063d-7*gammaaa*t1*t23*t29*t
     7         42*t87*t88+4.856339310426487d-5*t1*t19*t22*t29*t38*t43*t5
     8         2*t8-2.711369227677412d-4*t1*t152*t19*t22*t29*t37*t4*t43*
     9         t8+5.665729195497569d-5*t1*t19*t20*t23*t29*t43*t8+4.85633
     :         9310426487d-5*t1*t19*t22*t29*t37*t4*t42*t52-1.30472998239
     ;         6898d-5*t1*t19*t22*t23*t29*t38*t42-3.0443699589260953d-5*
     <         t1*t11*t19*t20*t29*t42
            t192 = 3.012632475197125d-5*t1*t114*t21*t4*t52*t77*t8-8.7055
     1         54958342553d-7*gammaaa*t1*t114*t19*t20*t52*t77*t8-2.71136
     2         9227677412d-4*t1*t114*t152*t22*t4*t5*t77*t8-1.20804460458
     3         37032d-4*t1*t141*t22*t4*t5*t52*t75*t8+2.0134076743061718d
     4         -5*t1*t141*t21*t23*t4*t75*t8-5.818111338348415d-7*gammaaa
     5         *t1*t141*t19*t20*t23*t75*t8+3.245596971136562d-5*t1*t140*
     6         t22*t23*t4*t5*t75-1.0818656570455205d-5*t1*t11*t140*t21*t
     7         4*t75+3.1262495549965463d-7*gammaaa*t1*t11*t140*t19*t20*t
     8         75+9.712678620852974d-5*t1*t114*t141*t22*t4*t5*t52-1.6187
     9         797701421625d-5*t1*t114*t141*t21*t23*t4+4.677761516032126
     :         d-7*gammaaa*t1*t114*t141*t19*t20*t23
            t193 = -2.711369227677412d-4*t1*t152*t32*t4*t56*t8*t87*t9+9.
     1         712678620852974d-5*t1*t125*t32*t52*t8*t87*t9-1.2080446045
     2         837032d-4*t1*t125*t32*t4*t52*t54*t8*t9+3.245596971136562d
     3         -5*t1*t126*t23*t32*t54*t8*t9+3.012632475197125d-5*t1*t31*
     4         t4*t52*t56*t8*t87-1.6187797701421625d-5*t1*t125*t23*t31*t
     5         8*t87-8.705554958342553d-7*gammabb*t1*t29*t30*t4*t52*t56*
     6         t87+4.677761516032126d-7*gammabb*t1*t125*t23*t29*t30*t87+
     7         2.0134076743061718d-5*t1*t125*t23*t31*t4*t54*t8-1.0818656
     8         570455205d-5*t1*t11*t126*t31*t54*t8-5.818111338348415d-7*
     9         gammabb*t1*t125*t23*t29*t30*t4*t54+3.1262495549965463d-7*
     :         gammabb*t1*t11*t126*t29*t30*t54
            t194 = 1.344462526483455d-2*t1*t117*t19*t37*t4*t52*t8*t9-3.6
     1         12104625279974d-3*t1*t117*t19*t23*t38*t8*t9+4.85633931042
     2         6487d-5*t1*t19*t29*t32*t38*t43*t52*t8+3.012632475197125d-
     3         5*t1*t116*t19*t29*t37*t4*t52*t8-2.711369227677412d-4*t1*t
     4         152*t19*t29*t32*t37*t4*t43*t8-8.093898850710812d-6*t1*t11
     5         6*t19*t23*t29*t38*t8-5.818111338348415d-7*gammabb*t1*t113
     6         *t19*t23*t37*t4*t75+3.1262495549965463d-7*gammabb*t1*t11*
     7         t113*t19*t38*t75+4.856339310426487d-5*t1*t19*t29*t32*t37*
     8         t4*t42*t52-8.705554958342553d-7*gammabb*t1*t114*t115*t19*
     9         t37*t4*t52-1.304729982396898d-5*t1*t19*t23*t29*t32*t38*t4
     :         2+5.665729195497569d-5*t1*t19*t23*t29*t30*t37*t4-3.044369
     ;         9589260953d-5*t1*t11*t19*t29*t30*t38+2.338880758016063d-7
     <         *gammabb*t1*t114*t115*t19*t23*t38
            t195 = 1/rhob**8.666666666666666d+0
            t196 = 1/rhob**8.333333333333334d+0
            t197 = 1/t7**6
            t198 = 3.0201115114592575d-4*t113*t75-2.1709370665479152d-6*
     1         gammabb*t114*t195
            t199 = -t1*t198*t23*t4*t5*t8*t9-2.711369227677412d-4*t1*t114
     1         *t152*t32*t4*t5*t77*t8-1.2080446045837032d-4*t1*t141*t32*
     2         t4*t5*t52*t75*t8+2.0134076743061718d-5*t1*t115*t23*t4*t5*
     3         t75*t8+2.68892505296691d-2*t1*t117*t29*t4*t43*t5*t52*t8+3
     4         .012632475197125d-5*t1*t114*t164*t4*t5*t52*t8+3.245596971
     5         136562d-5*t1*t140*t23*t32*t4*t5*t75+1.622798485568281d-4*
     6         t1*t11*t113*t4*t5*t75+9.712678620852974d-5*t1*t114*t141*t
     7         32*t4*t5*t52-8.705554958342553d-7*gammabb*t1*t163*t197*t4
     8         *t5*t52-7.224209250559948d-3*t1*t117*t23*t29*t4*t42*t5-1.
     9         7454334015045245d-6*gammabb*t1*t180*t196*t23*t4*t5+1.1331
     :         458390995136d-4*t1*t114*t115*t23*t4*t5-1.1665110279837858
     ;         d-6*gammabb*t1*t11*t114*t195*t4*t5
            t200 = 1/rhoa**7.333333333333333d+0
            t201 = 1/rhoa**7.666666666666667d+0
            t202 = -8.141013999554685d-7*t1*t200*t23*t4*t8*t87*t9+4.3744
     1         16354939197d-7*t1*t11*t201*t8*t87*t9-2.436256051774968d-6
     2         *t1*t168*t4*t52*t55*t8*t9+9.79374932813537d-7*t1*t186*t52
     3         *t55*t8*t9-1.8226671610363337d-6*t1*t152*t170*t186*t4*t8*
     4         t9+1.3090750511283933d-6*t1*t168*t200*t23*t8*t9
            t203 = -1.8226671610363337d-6*t1*t152*t29*t4*t43*t56*t8*t87+
     1         6.529166218756914d-7*t1*t125*t29*t43*t52*t8*t87+3.2645831
     2         09378457d-7*t1*t29*t4*t42*t52*t56*t87-1.7541605685120473d
     3         -7*t1*t125*t23*t29*t42*t87-8.12085350591656d-7*t1*t125*t2
     4         9*t4*t43*t52*t54*t8+2.1817917518806554d-7*t1*t126*t23*t29
     5         *t43*t54*t8+2.1817917518806554d-7*t1*t125*t23*t29*t4*t42*
     6         t54-1.1723435831237047d-7*t1*t11*t126*t29*t42*t54
            t204 = 3.264583109378457d-7*t1*t114*t19*t38*t52*t77*t8-1.822
     1         6671610363337d-6*t1*t114*t152*t19*t37*t4*t77*t8-8.1208535
     2         0591656d-7*t1*t141*t19*t37*t4*t52*t75*t8+2.18179175188065
     3         54d-7*t1*t141*t19*t23*t38*t75*t8+2.1817917518806554d-7*t1
     4         *t140*t19*t23*t37*t4*t75-1.1723435831237049d-7*t1*t11*t14
     5         0*t19*t38*t75+6.529166218756914d-7*t1*t114*t141*t19*t37*t
     6         4*t52-1.7541605685120473d-7*t1*t114*t141*t19*t23*t38
            t205 = 1/rhob**7.666666666666667d+0
            t206 = 1/rhob**7.333333333333333d+0
            t207 = -2.436256051774968d-6*t1*t180*t4*t5*t52*t76*t8-8.1410
     1         13999554685d-7*t1*t114*t206*t23*t4*t5*t8-1.82266716103633
     2         37d-6*t1*t152*t182*t197*t4*t5*t8+9.79374932813537d-7*t1*t
     3         197*t4*t5*t52*t76+1.3090750511283933d-6*t1*t180*t206*t23*
     4         t4*t5+4.374416354939197d-7*t1*t11*t114*t205*t4*t5
            fnc(iq) = fnc(iq)-8.336168033553906d-2*rhoa*rhob*t14*t16*t17
     1         *wght
            Amat(iq,D1_RA) = 8.336168033553906d-2*rhoa*rhob*t14*t17*t25*
     1         t26*wght+1.6672336067107812d-1*rhoa*rhob*t16*t17*t18*t24*
     2         wght+t14*t16*t28+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.336168033553906d-2*rhoa*rhob*t14*t17*t25*
     1         t34*wght+1.6672336067107812d-1*rhoa*rhob*t16*t17*t18*t33*
     2         wght+t14*t16*t36+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 8.336168033553906d-2*rhoa*rhob*t14*t17*t25
     1         *t41*wght+1.6672336067107812d-1*rhoa*rhob*t16*t17*t18*t40
     2         *wght-8.336168033553906d-2*rhoa*rhob*t14*t16*t39*wght+Cma
     3         t(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 8.336168033553906d-2*rhoa*rhob*t14*t17*t25
     1         *t46*wght+1.6672336067107812d-1*rhoa*rhob*t16*t17*t18*t45
     2         *wght-8.336168033553906d-2*rhoa*rhob*t14*t16*t44*wght+Cma
     3         t(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = -5.001700820132344d-1*rhoa*rhob*t16*t17
     1         *t49*t50*wght-1.6672336067107812d-1*rhoa*rhob*t14*t17*t47
     2         *t48*wght-3.3344672134215625d-1*rhoa*rhob*t17*t18*t24*t25
     3         *t26*wght+t14*t25*t64+t16*t18*t62+t14*t16*t60+Amat2(iq,D2
     4         _RA_RA)
            Amat2(iq,D2_RA_RB) = -5.001700820132344d-1*rhoa*rhob*t16*t17
     1         *t24*t33*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t17
     2         *t26*t34*t47*wght+t14*t25*t71+t16*t18*t69+t18*t25*t67+t14
     3         *t16*t66+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = -5.001700820132344d-1*rhoa*rhob*t16*t17
     1         *t49*t73*wght-1.6672336067107812d-1*rhoa*rhob*t14*t17*t47
     2         *t72*wght-3.3344672134215625d-1*rhoa*rhob*t17*t18*t25*t33
     3         *t34*wght+t14*t25*t86+t16*t18*t84+t14*t16*t82+Amat2(iq,D2
     4         _RB_RB)
            Cmat2(iq,D2_RA_GAA) = -5.001700820132344d-1*rhoa*rhob*t16*t1
     1         7*t24*t40*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t1
     2         7*t26*t41*t47*wght+t14*t25*t98+t16*t18*t96+t18*t25*t94+t1
     3         4*t16*t93+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = -5.001700820132344d-1*rhoa*rhob*t16*t1
     1         7*t24*t45*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t1
     2         7*t26*t46*t47*wght+t101*t18*t25+t105*t14*t25+t103*t16*t18
     3         +t100*t14*t16+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = -5.001700820132344d-1*rhoa*rhob*t16*t1
     1         7*t33*t40*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t1
     2         7*t34*t41*t47*wght+t108*t18*t25+t112*t14*t25+t110*t16*t18
     3         +t107*t14*t16+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = -5.001700820132344d-1*rhoa*rhob*t16*t1
     1         7*t33*t45*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t1
     2         7*t34*t46*t47*wght+t120*t18*t25+t124*t14*t25+t122*t16*t18
     3         +t119*t14*t16+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = -5.001700820132344d-1*rhoa*rhob*t129*
     1         t16*t17*t49*wght-1.6672336067107812d-1*rhoa*rhob*t128*t14
     2         *t17*t47*wght-3.3344672134215625d-1*rhoa*rhob*t17*t18*t25
     3         *t40*t41*wght-8.336168033553906d-2*rhoa*rhob*t127*t14*t16
     4         *wght+t133*t14*t25+t131*t16*t18+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = -5.001700820132344d-1*rhoa*rhob*t16*t
     1         17*t40*t45*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t
     2         17*t41*t46*t47*wght-8.336168033553906d-2*rhoa*rhob*t134*t
     3         14*t16*wght+t139*t18*t25+t138*t14*t25+t136*t16*t18+Cmat2(
     4         iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -5.001700820132344d-1*rhoa*rhob*t144*
     1         t16*t17*t49*wght-1.6672336067107812d-1*rhoa*rhob*t14*t143
     2         *t17*t47*wght-3.3344672134215625d-1*rhoa*rhob*t17*t18*t25
     3         *t45*t46*wght-8.336168033553906d-2*rhoa*rhob*t14*t142*t16
     4         *wght+t14*t148*t25+t146*t16*t18+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = t14*t25*(8.336168033553906d-2*rhoa*r
     1         hob*t17*(-3.2419488676002184d-5*t1*t11*t153*t154*t8*t87*t
     2         9-4.7983677678995934d+1*t1*t22*t23*t4*t5*t61*t8*t9+1.5994
     3         559226331975d+1*t1*t11*t21*t4*t61*t8*t9-4.621921708859857
     4         d-1*gammaaa*t1*t11*t19*t20*t61*t8*t9+2.3449365598931193d+
     5         1*t1*t22*t4*t5*t52*t58*t8*t9-2.3991838839497967d+1*t1*t23
     6         *t24*t4*t5*t58*t8*t9-3.908227599821865d+0*t1*t21*t23*t4*t
     7         58*t8*t9+1.1293541591970727d-1*gammaaa*t1*t19*t20*t23*t58
     8         *t8*t9-5.331519742110659d+0*t1*t11*t24*t4*t57*t8*t9+2.605
     9         48506654791d+0*t1*t22*t23*t4*t57*t8*t9+6.932882563289786d
     :         -1*gammaaa*t1*t11*t19*t24*t56*t8*t9-3.388062477591218d-1*
     ;         gammaaa*t1*t19*t22*t23*t56*t8*t9-5.537268584495724d-3*t1*
     <         t11*t24*t53*t54*t55*t8*t9+2.706033420905421d-3*t1*t22*t23
     =         *t53*t54*t55*t8*t9+9.02011140301807d-3*t1*t11*t155*t53*t5
     >         4*t8*t9+4.7983677678995934d+1*t1*t24*t4*t5*t51*t52*t8*t9+
     ?         7.81645519964373d+0*t1*t21*t4*t51*t52*t8*t9-2.25870831839
     @         41453d-1*gammaaa*t1*t19*t20*t51*t52*t8*t9-3.9082275998218
     1         66d+0*t1*t157*t23*t4*t5*t8*t9+2.3991838839497967d+1*t1*t1
     2         1*t158*t4*t5*t8*t9-2.3449365598931193d+1*t1*t151*t152*t4*
     3         t5*t8*t9-1.5994559226331975d+1*t1*t21*t22*t23*t24*t4*t8*t
     4         9+1.4474917036377277d+0*t1*t11*t2*t4*t8*t9+4.621921708859
     5         857d-1*gammaaa*t1*t19*t20*t22*t23*t24*t8*t9-4.81021215954
     6         30875d-1*gammaaa*t1*t11*t156*t19*t8*t9)*wght+1.6672336067
     7         107812d-1*rhoa*rhob*t27*t63*wght+1.6672336067107812d-1*rh
     8         ob*t17*t63*wght+8.336168033553906d-2*rhoa*rhob*t26*t59*wg
     9         ht+1.6672336067107812d-1*rhob*t26*t27*wght-t28*t63-2*t26*
     :         t60)+t14*t16*(-8.336168033553906d-2*rhoa*rhob*(-4.3710572
     ;         55553858d-5*t1*t11*t153*t154*t8*t87*t9+3.161632827207898d
     <         +1*t1*t22*t4*t5*t52*t58*t8*t9-5.269388045346496d+0*t1*t21
     =         *t23*t4*t58*t8*t9+1.5226864744792862d-1*gammaaa*t1*t19*t2
     >         0*t23*t58*t8*t9+3.5129253635643304d+0*t1*t22*t23*t4*t57*t
     ?         8*t9-4.5680594234378585d-1*gammaaa*t1*t19*t22*t23*t56*t8*
     @         t9+3.648492774340223d-3*t1*t22*t23*t53*t54*t55*t8*t9+1.21
     1         61642581134074d-2*t1*t11*t155*t53*t54*t8*t9+1.05387760906
     2         92992d+1*t1*t21*t4*t51*t52*t8*t9-3.0453729489585724d-1*ga
     3         mmaaa*t1*t19*t20*t51*t52*t8*t9-5.269388045346497d+0*t1*t1
     4         57*t23*t4*t5*t8*t9-3.161632827207898d+1*t1*t151*t152*t4*t
     5         5*t8*t9+1.9516252019801836d+0*t1*t11*t2*t4*t8*t9-6.485516
     6         465374737d-1*gammaaa*t1*t11*t156*t19*t8*t9)*wght-2.500850
     7         410066172d-1*rhob*t59*wght)+t14*t47*(-3.3344672134215625d
     8         -1*rhoa*rhob*t17*t26*t63*wght-1.6672336067107812d-1*rhoa*
     9         rhob*t27*t48*wght-1.6672336067107812d-1*rhob*t17*t48*wght
     :         -2*t26*t64)+t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t17
     ;         *t24*t63*wght-3.3344672134215625d-1*rhoa*rhob*t17*t26*t61
     <         *wght-3.3344672134215625d-1*rhoa*rhob*t24*t26*t27*wght-3.
     =         3344672134215625d-1*rhob*t17*t24*t26*wght-2.0d+0*t24*t64-
     >         t26*t62)+t16*t18*(3.3344672134215625d-1*rhoa*rhob*t27*t61
     ?         *wght+3.3344672134215625d-1*rhob*t17*t61*wght+1.667233606
     @         7107812d-1*rhoa*rhob*t24*t59*wght+3.3344672134215625d-1*r
     1         hob*t24*t27*wght+1.6672336067107812d-1*rhoa*rhob*t158*t17
     2         *wght-2.0d+0*t28*t61-4.0d+0*t24*t60)+t16*t49*(-1.00034016
     3         40264689d+0*rhoa*rhob*t17*t24*t61*wght-5.001700820132344d
     4         -1*rhoa*rhob*t27*t50*wght-5.001700820132344d-1*rhob*t17*t
     5         50*wght-3.0d+0*t24*t62)+1.5005102460397032d+0*rhoa*rhob*t
     6         17*t25*t26*t49*t50*wght+1.0003401640264689d+0*rhoa*rhob*t
     7         17*t18*t24*t47*t48*wght+5.001700820132344d-1*rhoa*rhob*t1
     8         4*t149*t17*t26**3*wght+2.0006803280529378d+0*rhoa*rhob*t1
     9         50*t16*t17*t24**3*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t14*t25*(8.336168033553906d-2*rhoa*r
     1         hob*t17*(-2.3991838839497967d+1*t1*t22*t23*t4*t5*t68*t8*t
     2         9+7.997279613165988d+0*t1*t11*t21*t4*t68*t8*t9-2.31096085
     3         44299284d-1*gammaaa*t1*t11*t19*t20*t68*t8*t9-2.3991838839
     4         497967d+1*t1*t23*t32*t4*t5*t61*t8*t9+7.816455199643731d+0
     5         *t1*t32*t4*t5*t52*t58*t8*t9+8.684950221826366d-1*t1*t23*t
     6         32*t4*t57*t8*t9-1.1293541591970728d-1*gammaaa*t1*t19*t23*
     7         t32*t56*t8*t9+9.02011140301807d-4*t1*t23*t32*t53*t54*t55*
     8         t8*t9+4.7983677678995934d+1*t1*t22*t24*t32*t4*t5*t52*t8*t
     9         9+5.21097013309582d+0*t1*t21*t22*t32*t4*t52*t8*t9-1.50580
     :         5545596097d-1*gammaaa*t1*t19*t20*t22*t32*t52*t8*t9-2.3449
     ;         365598931193d+1*t1*t152*t32*t4*t5*t51*t8*t9+2.39918388394
     <         97967d+1*t1*t11*t159*t4*t5*t8*t9-7.997279613165988d+0*t1*
     =         t21*t23*t24*t32*t4*t8*t9+2.3109608544299284d-1*gammaaa*t1
     >         *t19*t20*t23*t24*t32*t8*t9+7.997279613165988d+0*t1*t11*t3
     ?         1*t4*t5*t61*t8-1.302742533273955d+0*t1*t23*t31*t4*t5*t58*
     @         t8-2.8949834072754554d-1*t1*t11*t31*t4*t57*t8+3.764513863
     1         9902426d-2*gammaaa*t1*t11*t19*t31*t56*t8-3.00670380100602
     2         3d-4*t1*t11*t31*t53*t54*t55*t8+2.60548506654791d+0*t1*t31
     3         *t4*t5*t51*t52*t8-7.997279613165988d+0*t1*t22*t23*t24*t31
     4         *t4*t5*t8+2.6657598710553293d+0*t1*t11*t21*t24*t31*t4*t8-
     5         8.684950221826366d-1*t1*t21*t22*t23*t31*t4*t8-7.703202848
     6         09976d-2*gammaaa*t1*t11*t19*t20*t24*t31*t8+2.509675909326
     7         8284d-2*gammaaa*t1*t19*t20*t22*t23*t31*t8-2.3109608544299
     8         284d-1*gammabb*t1*t11*t29*t30*t4*t5*t61+3.764513863990242
     9         6d-2*gammabb*t1*t23*t29*t30*t4*t5*t58+8.36558636442276d-3
     :         *gammabb*t1*t11*t29*t30*t4*t57-1.0878254352039793d-3*gamm
     ;         aaa*gammabb*t1*t11*t19*t29*t30*t56+8.688422965168586d-6*g
     <         ammabb*t1*t11*t29*t30*t53*t54*t55-7.529027727980485d-2*ga
     =         mmabb*t1*t29*t30*t4*t5*t51*t52+2.3109608544299284d-1*gamm
     >         abb*t1*t22*t23*t24*t29*t30*t4*t5-7.70320284809976d-2*gamm
     ?         abb*t1*t11*t21*t24*t29*t30*t4+2.5096759093268284d-2*gamma
     @         bb*t1*t21*t22*t23*t29*t30*t4+2.2259819709672812d-3*gammaa
     1         a*gammabb*t1*t11*t19*t20*t24*t29*t30-7.252169568026529d-4
     2         *gammaaa*gammabb*t1*t19*t20*t22*t23*t29*t30)*wght+8.33616
     3         8033553906d-2*rhoa*rhob*t27*t70*wght+8.336168033553906d-2
     4         *rhob*t17*t70*wght+8.336168033553906d-2*rhoa*rhob*t26*t65
     5         *wght+8.336168033553906d-2*rhoa*rhob*t35*t63*wght+8.33616
     6         8033553906d-2*rhoa*t17*t63*wght+8.336168033553906d-2*rhob
     7         *t26*t35*wght+8.336168033553906d-2*rhoa*t26*t27*wght+8.33
     8         6168033553906d-2*t17*t26*wght-t28*t70-t26*t66-t34*t60)+t1
     9         4*t16*(-8.336168033553906d-2*rhoa*rhob*(1.053877609069299
     :         3d+1*t1*t32*t4*t5*t52*t58*t8*t9+1.17097512118811d+0*t1*t2
     ;         3*t32*t4*t57*t8*t9-1.5226864744792862d-1*gammaaa*t1*t19*t
     <         23*t32*t56*t8*t9+1.2161642581134076d-3*t1*t23*t32*t53*t54
     =         *t55*t8*t9+7.025850727128662d+0*t1*t21*t22*t32*t4*t52*t8*
     >         t9-2.0302486326390484d-1*gammaaa*t1*t19*t20*t22*t32*t52*t
     ?         8*t9-3.161632827207898d+1*t1*t152*t32*t4*t5*t51*t8*t9-1.7
     @         564626817821655d+0*t1*t23*t31*t4*t5*t58*t8-3.903250403960
     1         3666d-1*t1*t11*t31*t4*t57*t8+5.07562158159762d-2*gammaaa*
     2         t1*t11*t19*t31*t56*t8-4.053880860378026d-4*t1*t11*t31*t53
     3         *t54*t55*t8+3.512925363564331d+0*t1*t31*t4*t5*t51*t52*t8-
     4         1.17097512118811d+0*t1*t21*t22*t23*t31*t4*t8+3.3837477210
     5         650807d-2*gammaaa*t1*t19*t20*t22*t23*t31*t8+5.07562158159
     6         7621d-2*gammabb*t1*t23*t29*t30*t4*t5*t58+1.12791590702169
     7         33d-2*gammabb*t1*t11*t29*t30*t4*t57-1.4666940952847696d-3
     8         *gammaaa*gammabb*t1*t11*t19*t29*t30*t56+1.171443344488434
     9         d-5*gammabb*t1*t11*t29*t30*t53*t54*t55-1.0151243163195242
     :         d-1*gammabb*t1*t29*t30*t4*t5*t51*t52+3.3837477210650807d-
     ;         2*gammabb*t1*t21*t22*t23*t29*t30*t4-9.777960635231797d-4*
     <         gammaaa*gammabb*t1*t19*t20*t22*t23*t29*t30)*wght-1.667233
     =         6067107812d-1*rhob*t65*wght-8.336168033553906d-2*rhoa*t59
     >         *wght-1.6672336067107812d-1*t27*wght)+t14*t47*(-3.3344672
     ?         134215625d-1*rhoa*rhob*t17*t26*t70*wght-1.667233606710781
     @         2d-1*rhoa*rhob*t35*t48*wght-1.6672336067107812d-1*rhoa*t1
     1         7*t48*wght-2*t34*t64)+t18*t25*(-3.3344672134215625d-1*rho
     2         a*rhob*t17*t24*t70*wght-3.3344672134215625d-1*rhoa*rhob*t
     3         17*t26*t68*wght-3.3344672134215625d-1*rhoa*rhob*t24*t26*t
     4         35*wght-3.3344672134215625d-1*rhoa*t17*t24*t26*wght-2.0d+
     5         0*t33*t64-t34*t62)+t16*t18*(1.6672336067107812d-1*rhoa*rh
     6         ob*t27*t68*wght+1.6672336067107812d-1*rhob*t17*t68*wght+1
     7         .6672336067107812d-1*rhoa*rhob*t24*t65*wght+1.66723360671
     8         07812d-1*rhoa*rhob*t35*t61*wght+1.6672336067107812d-1*rho
     9         a*t17*t61*wght+1.6672336067107812d-1*rhob*t24*t35*wght+1.
     :         6672336067107812d-1*rhoa*t24*t27*wght+1.6672336067107812d
     ;         -1*t17*t24*wght+1.6672336067107812d-1*rhoa*rhob*t159*t17*
     <         wght-2.0d+0*t28*t68-2.0d+0*t24*t66-2.0d+0*t33*t60)+t16*t4
     =         9*(-1.0003401640264689d+0*rhoa*rhob*t17*t24*t68*wght-5.00
     >         1700820132344d-1*rhoa*rhob*t35*t50*wght-5.001700820132344
     ?         d-1*rhoa*t17*t50*wght-3.0d+0*t33*t62)+t25*t49*(5.00170082
     @         0132344d-1*rhoa*rhob*t17*t34*t50*wght+1.0003401640264689d
     1         +0*rhoa*rhob*t17*t24*t26*t33*wght)+t18*t47*(3.33446721342
     2         15625d-1*rhoa*rhob*t17*t33*t48*wght+6.668934426843125d-1*
     3         rhoa*rhob*t17*t24*t26*t34*wght)+2.0006803280529378d+0*rho
     4         a*rhob*t150*t16*t17*t33*t50*wght+5.001700820132344d-1*rho
     5         a*rhob*t14*t149*t17*t34*t48*wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = t14*t16*(-8.336168033553906d-2*rhoa*
     1         rhob*(1.0538776090692993d+1*t1*t22*t4*t5*t52*t8*t80*t9-1.
     2         7564626817821655d+0*t1*t21*t23*t4*t8*t80*t9+5.07562158159
     3         7621d-2*gammaaa*t1*t19*t20*t23*t8*t80*t9+3.51292536356433
     4         1d+0*t1*t21*t4*t52*t79*t8*t9-1.0151243163195242d-1*gammaa
     5         a*t1*t19*t20*t52*t79*t8*t9-3.161632827207898d+1*t1*t152*t
     6         22*t4*t5*t79*t8*t9+1.17097512118811d+0*t1*t22*t23*t4*t5*t
     7         78*t8-3.9032504039603666d-1*t1*t11*t21*t4*t78*t8+1.127915
     8         9070216935d-2*gammaaa*t1*t11*t19*t20*t78*t8+7.02585072712
     9         8662d+0*t1*t22*t31*t32*t4*t5*t52*t8-1.17097512118811d+0*t
     :         1*t21*t23*t31*t32*t4*t8+3.3837477210650807d-2*gammaaa*t1*
     ;         t19*t20*t23*t31*t32*t8-1.5226864744792862d-1*gammabb*t1*t
     <         22*t23*t29*t4*t5*t77+5.075621581597621d-2*gammabb*t1*t11*
     =         t21*t29*t4*t77-1.4666940952847698d-3*gammaaa*gammabb*t1*t
     >         11*t19*t20*t29*t77+1.2161642581134076d-3*t1*t22*t23*t4*t5
     ?         *t74*t75*t76-4.053880860378026d-4*t1*t11*t21*t4*t74*t75*t
     @         76+1.171443344488434d-5*gammaaa*t1*t11*t19*t20*t74*t75*t7
     1         6-2.0302486326390484d-1*gammabb*t1*t22*t29*t30*t32*t4*t5*
     2         t52+3.3837477210650807d-2*gammabb*t1*t21*t23*t29*t30*t32*
     3         t4-9.777960635231797d-4*gammaaa*gammabb*t1*t19*t20*t23*t2
     4         9*t30*t32)*wght-8.336168033553906d-2*rhob*t81*wght-1.6672
     5         336067107812d-1*rhoa*t65*wght-1.6672336067107812d-1*t35*w
     6         ght)+t14*t25*(8.336168033553906d-2*rhoa*rhob*t17*(7.81645
     7         5199643731d+0*t1*t22*t4*t5*t52*t8*t80*t9-2.39918388394979
     8         67d+1*t1*t23*t24*t4*t5*t8*t80*t9-1.302742533273955d+0*t1*
     9         t21*t23*t4*t8*t80*t9+3.7645138639902426d-2*gammaaa*t1*t19
     :         *t20*t23*t8*t80*t9+4.7983677678995934d+1*t1*t24*t4*t5*t52
     ;         *t79*t8*t9+2.60548506654791d+0*t1*t21*t4*t52*t79*t8*t9-7.
     <         529027727980485d-2*gammaaa*t1*t19*t20*t52*t79*t8*t9-2.344
     =         9365598931193d+1*t1*t152*t22*t4*t5*t79*t8*t9-4.7983677678
     >         995934d+1*t1*t23*t32*t4*t5*t68*t8*t9+2.3991838839497967d+
     ?         1*t1*t11*t160*t4*t5*t8*t9-5.331519742110659d+0*t1*t11*t24
     @         *t4*t5*t78*t8+8.684950221826366d-1*t1*t22*t23*t4*t5*t78*t
     1         8-2.8949834072754554d-1*t1*t11*t21*t4*t78*t8+8.3655863644
     2         2276d-3*gammaaa*t1*t11*t19*t20*t78*t8+1.5994559226331975d
     3         +1*t1*t11*t31*t4*t5*t68*t8+5.21097013309582d+0*t1*t22*t31
     4         *t32*t4*t5*t52*t8-1.5994559226331975d+1*t1*t23*t24*t31*t3
     5         2*t4*t5*t8-8.684950221826366d-1*t1*t21*t23*t31*t32*t4*t8+
     6         2.5096759093268284d-2*gammaaa*t1*t19*t20*t23*t31*t32*t8+6
     7         .932882563289786d-1*gammabb*t1*t11*t24*t29*t4*t5*t77-1.12
     8         93541591970728d-1*gammabb*t1*t22*t23*t29*t4*t5*t77+3.7645
     9         138639902426d-2*gammabb*t1*t11*t21*t29*t4*t77-1.087825435
     :         2039793d-3*gammaaa*gammabb*t1*t11*t19*t20*t29*t77-5.53726
     ;         8584495724d-3*t1*t11*t24*t4*t5*t74*t75*t76+9.020111403018
     <         07d-4*t1*t22*t23*t4*t5*t74*t75*t76-3.006703801006023d-4*t
     =         1*t11*t21*t4*t74*t75*t76+8.688422965168586d-6*gammaaa*t1*
     >         t11*t19*t20*t74*t75*t76-4.621921708859857d-1*gammabb*t1*t
     ?         11*t29*t30*t4*t5*t68-1.505805545596097d-1*gammabb*t1*t22*
     @         t29*t30*t32*t4*t5*t52+4.621921708859857d-1*gammabb*t1*t23
     1         *t24*t29*t30*t32*t4*t5+2.5096759093268284d-2*gammabb*t1*t
     2         21*t23*t29*t30*t32*t4-7.252169568026529d-4*gammaaa*gammab
     3         b*t1*t19*t20*t23*t29*t30*t32)*wght+8.336168033553906d-2*r
     4         hoa*rhob*t26*t81*wght+1.6672336067107812d-1*rhoa*rhob*t35
     5         *t70*wght+1.6672336067107812d-1*rhoa*t17*t70*wght+1.66723
     6         36067107812d-1*rhoa*t26*t35*wght-t28*t85-2*t34*t66)+t14*t
     7         47*(-1.6672336067107812d-1*rhoa*rhob*t17*t26*t85*wght-1.6
     8         672336067107812d-1*rhoa*rhob*t17*t34*t70*wght-1.667233606
     9         7107812d-1*rhoa*rhob*t26*t34*t35*wght-1.6672336067107812d
     :         -1*rhoa*t17*t26*t34*wght-2*t34*t71)+t18*t25*(-1.667233606
     ;         7107812d-1*rhoa*rhob*t17*t24*t85*wght-1.6672336067107812d
     <         -1*rhoa*rhob*t17*t26*t83*wght-1.6672336067107812d-1*rhoa*
     =         rhob*t17*t33*t70*wght-1.6672336067107812d-1*rhoa*rhob*t17
     >         *t34*t68*wght-1.6672336067107812d-1*rhoa*rhob*t24*t34*t35
     ?         *wght-1.6672336067107812d-1*rhoa*rhob*t26*t33*t35*wght-1.
     @         6672336067107812d-1*rhoa*t17*t24*t34*wght-1.6672336067107
     1         812d-1*rhoa*t17*t26*t33*wght-2.0d+0*t33*t71-t34*t69)+t16*
     2         t49*(-5.001700820132344d-1*rhoa*rhob*t17*t24*t83*wght-5.0
     3         01700820132344d-1*rhoa*rhob*t17*t33*t68*wght-5.0017008201
     4         32344d-1*rhoa*rhob*t24*t33*t35*wght-5.001700820132344d-1*
     5         rhoa*t17*t24*t33*wght-3.0d+0*t33*t69)+t16*t18*(1.66723360
     6         67107812d-1*rhoa*rhob*t24*t81*wght+3.3344672134215625d-1*
     7         rhoa*rhob*t35*t68*wght+3.3344672134215625d-1*rhoa*t17*t68
     8         *wght+3.3344672134215625d-1*rhoa*t24*t35*wght+1.667233606
     9         7107812d-1*rhoa*rhob*t160*t17*wght-2.0d+0*t28*t83-4.0d+0*
     :         t33*t66)+t18*t47*(3.3344672134215625d-1*rhoa*rhob*t17*t26
     ;         *t33*t34*wght-2*t34*t67)+t25*t49*(5.001700820132344d-1*rh
     <         oa*rhob*t17*t24*t33*t34*wght-3.0d+0*t33*t67)+2.0006803280
     =         529378d+0*rhoa*rhob*t150*t16*t17*t24*t73*wght+5.001700820
     >         132344d-1*rhoa*rhob*t14*t149*t17*t26*t72*wght+Amat3(iq,D3
     ?         _RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = t14*t25*(8.336168033553906d-2*rhoa*r
     1         hob*t17*(-4.7983677678995934d+1*t1*t23*t32*t4*t5*t8*t83*t
     2         9+2.3449365598931193d+1*t1*t32*t4*t5*t52*t8*t80*t9-2.3991
     3         838839497967d+1*t1*t23*t33*t4*t5*t8*t80*t9+4.798367767899
     4         5934d+1*t1*t33*t4*t5*t52*t79*t8*t9-3.908227599821866d+0*t
     5         1*t166*t23*t4*t5*t8*t9+2.3991838839497967d+1*t1*t11*t167*
     6         t4*t5*t8*t9-2.3449365598931193d+1*t1*t152*t165*t4*t5*t8*t
     7         9+1.5994559226331975d+1*t1*t11*t31*t4*t5*t8*t83-4.6219217
     8         08859857d-1*gammabb*t1*t11*t29*t30*t4*t5*t83-3.9082275998
     9         21865d+0*t1*t23*t31*t4*t5*t8*t80+1.1293541591970727d-1*ga
     :         mmabb*t1*t23*t29*t30*t4*t5*t80+7.81645519964373d+0*t1*t31
     ;         *t4*t5*t52*t79*t8-5.331519742110659d+0*t1*t11*t33*t4*t5*t
     <         78*t8+2.60548506654791d+0*t1*t23*t32*t4*t5*t78*t8+1.44749
     =         17036377277d+0*t1*t11*t4*t5*t6*t8-1.5994559226331975d+1*t
     >         1*t23*t31*t32*t33*t4*t5*t8-2.2587083183941453d-1*gammabb*
     ?         t1*t29*t30*t4*t5*t52*t79+6.932882563289786d-1*gammabb*t1*
     @         t11*t29*t33*t4*t5*t77-3.388062477591218d-1*gammabb*t1*t23
     1         *t29*t32*t4*t5*t77-5.537268584495724d-3*t1*t11*t33*t4*t5*
     2         t74*t75*t76+2.706033420905421d-3*t1*t23*t32*t4*t5*t74*t75
     3         *t76+9.02011140301807d-3*t1*t11*t163*t4*t5*t74*t75+4.6219
     4         21708859857d-1*gammabb*t1*t23*t29*t30*t32*t33*t4*t5-4.810
     5         2121595430875d-1*gammabb*t1*t11*t164*t29*t4*t5-3.24194886
     6         76002184d-5*t1*t11*t114*t161*t162*t4*t5)*wght+1.667233606
     7         7107812d-1*rhoa*rhob*t35*t85*wght+1.6672336067107812d-1*r
     8         hoa*t17*t85*wght+8.336168033553906d-2*rhoa*rhob*t34*t81*w
     9         ght+1.6672336067107812d-1*rhoa*t34*t35*wght-t36*t85-2*t34
     :         *t82)+t14*t16*(-8.336168033553906d-2*rhoa*rhob*(3.1616328
     ;         27207898d+1*t1*t32*t4*t5*t52*t8*t80*t9-5.269388045346497d
     <         +0*t1*t166*t23*t4*t5*t8*t9-3.161632827207898d+1*t1*t152*t
     =         165*t4*t5*t8*t9-5.269388045346496d+0*t1*t23*t31*t4*t5*t8*
     >         t80+1.5226864744792862d-1*gammabb*t1*t23*t29*t30*t4*t5*t8
     ?         0+1.0538776090692992d+1*t1*t31*t4*t5*t52*t79*t8+3.5129253
     @         635643304d+0*t1*t23*t32*t4*t5*t78*t8+1.9516252019801836d+
     1         0*t1*t11*t4*t5*t6*t8-3.0453729489585724d-1*gammabb*t1*t29
     2         *t30*t4*t5*t52*t79-4.568059423437859d-1*gammabb*t1*t23*t2
     3         9*t32*t4*t5*t77+3.648492774340223d-3*t1*t23*t32*t4*t5*t74
     4         *t75*t76+1.2161642581134074d-2*t1*t11*t163*t4*t5*t74*t75-
     5         6.485516465374737d-1*gammabb*t1*t11*t164*t29*t4*t5-4.3710
     6         57255553858d-5*t1*t11*t114*t161*t162*t4*t5)*wght-2.500850
     7         410066172d-1*rhoa*t81*wght)+t14*t47*(-3.3344672134215625d
     8         -1*rhoa*rhob*t17*t34*t85*wght-1.6672336067107812d-1*rhoa*
     9         rhob*t35*t72*wght-1.6672336067107812d-1*rhoa*t17*t72*wght
     :         -2*t34*t86)+t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t17
     ;         *t33*t85*wght-3.3344672134215625d-1*rhoa*rhob*t17*t34*t83
     <         *wght-3.3344672134215625d-1*rhoa*rhob*t33*t34*t35*wght-3.
     =         3344672134215625d-1*rhoa*t17*t33*t34*wght-2.0d+0*t33*t86-
     >         t34*t84)+t16*t18*(3.3344672134215625d-1*rhoa*rhob*t35*t83
     ?         *wght+3.3344672134215625d-1*rhoa*t17*t83*wght+1.667233606
     @         7107812d-1*rhoa*rhob*t33*t81*wght+3.3344672134215625d-1*r
     1         hoa*t33*t35*wght+1.6672336067107812d-1*rhoa*rhob*t167*t17
     2         *wght-2.0d+0*t36*t83-4.0d+0*t33*t82)+t16*t49*(-1.00034016
     3         40264689d+0*rhoa*rhob*t17*t33*t83*wght-5.001700820132344d
     4         -1*rhoa*rhob*t35*t73*wght-5.001700820132344d-1*rhoa*t17*t
     5         73*wght-3.0d+0*t33*t84)+1.5005102460397032d+0*rhoa*rhob*t
     6         17*t25*t34*t49*t73*wght+1.0003401640264689d+0*rhoa*rhob*t
     7         17*t18*t33*t47*t72*wght+5.001700820132344d-1*rhoa*rhob*t1
     8         4*t149*t17*t34**3*wght+2.0006803280529378d+0*rhoa*rhob*t1
     9         50*t16*t17*t33**3*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t14*t25*(8.336168033553906d-2*rhoa*
     1         rhob*t27*t97*wght+8.336168033553906d-2*rhob*t17*t97*wght+
     2         8.336168033553906d-2*rhoa*rhob*t17*(-2.3991838839497967d+
     3         1*t1*t22*t23*t4*t5*t8*t9*t95+7.997279613165988d+0*t1*t11*
     4         t21*t4*t8*t9*t95-2.3109608544299284d-1*gammaaa*t1*t11*t19
     5         *t20*t8*t9*t95+1.563291039928746d+1*t1*t22*t4*t5*t52*t8*t
     6         9*t91-2.3991838839497967d+1*t1*t23*t24*t4*t5*t8*t9*t91-2.
     7         60548506654791d+0*t1*t21*t23*t4*t8*t9*t91+7.5290277279804
     8         85d-2*gammaaa*t1*t19*t20*t23*t8*t9*t91+2.0764757191858965
     9         d-3*gammaaa*t1*t11*t24*t54*t8*t9*t90-6.765083552263552d-4
     :         *gammaaa*t1*t22*t23*t54*t8*t9*t90+3.502977035285916d-2*t1
     ;         *t19*t22*t4*t52*t8*t89*t9-5.376021376855886d-2*t1*t19*t23
     <         *t24*t4*t8*t89*t9-1.0122495641124629d-3*gammaaa*t1*t22*t5
     =         2*t8*t87*t88*t9+1.5535001344756127d-3*gammaaa*t1*t23*t24*
     >         t8*t87*t88*t9+1.2157308253500822d-5*t1*t11*t172*t53*t8*t8
     ?         7*t9-7.591871730843472d-4*gammaaa*t1*t170*t23*t8*t87*t9-1
     @         .612806413056766d-1*t1*t19*t23*t37*t4*t61*t8*t9+8.6661032
     1         04112231d-2*t1*t11*t19*t38*t61*t8*t9+5.254465552928875d-2
     2         *t1*t19*t37*t4*t52*t58*t8*t9-1.4116926989963413d-2*t1*t19
     3         *t23*t38*t58*t8*t9+1.0979832103304875d-1*t1*t11*t19*t56*t
     4         8*t9-2.8187848134431465d-3*gammaaa*t1*t11*t54*t55*t8*t9+6
     5         .063600883031948d-6*t1*t168*t169*t23*t53*t8*t9+2.82338539
     6         79926826d-2*t1*t19*t38*t51*t52*t8*t9+3.225612826113532d-1
     7         *t1*t19*t22*t24*t37*t4*t52*t8*t9-1.5763396658786624d-1*t1
     8         *t152*t19*t37*t4*t51*t8*t9-3.908227599821866d+0*t1*t173*t
     9         23*t4*t5*t8*t9+2.3991838839497967d+1*t1*t11*t174*t4*t5*t8
     :         *t9+5.83829505880986d-3*t1*t171*t19*t23*t4*t8*t9-8.666103
     ;         204112231d-2*t1*t19*t22*t23*t24*t38*t8*t9-2.0220907476261
     <         873d-1*t1*t11*t19*t20*t24*t8*t9+6.587899261982924d-2*t1*t
     =         19*t20*t22*t23*t8*t9)*wght+8.336168033553906d-2*rhoa*rhob
     >         *t26*t92*wght+8.336168033553906d-2*rhoa*rhob*t39*t63*wght
     ?         +8.336168033553906d-2*rhob*t26*t39*wght-t28*t97-t26*t93-t
     @         41*t60)+t14*t47*(-3.3344672134215625d-1*rhoa*rhob*t17*t26
     1         *t97*wght-1.6672336067107812d-1*rhoa*rhob*t39*t48*wght-2*
     2         t41*t64)+t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t17*t2
     3         4*t97*wght-3.3344672134215625d-1*rhoa*rhob*t17*t26*t95*wg
     4         ht-3.3344672134215625d-1*rhoa*rhob*t24*t26*t39*wght-2.0d+
     5         0*t40*t64-t41*t62)+t16*t18*(1.6672336067107812d-1*rhoa*rh
     6         ob*t27*t95*wght+1.6672336067107812d-1*rhob*t17*t95*wght+1
     7         .6672336067107812d-1*rhoa*rhob*t24*t92*wght+1.66723360671
     8         07812d-1*rhoa*rhob*t39*t61*wght+1.6672336067107812d-1*rho
     9         b*t24*t39*wght+1.6672336067107812d-1*rhoa*rhob*t17*t174*w
     :         ght-2.0d+0*t28*t95-2.0d+0*t24*t93-2.0d+0*t40*t60)+t16*t49
     ;         *(-1.0003401640264689d+0*rhoa*rhob*t17*t24*t95*wght-5.001
     <         700820132344d-1*rhoa*rhob*t39*t50*wght-3.0d+0*t40*t62)+t1
     =         4*t16*(-1.6672336067107812d-1*rhob*t92*wght-8.33616803355
     >         3906d-2*rhoa*rhob*(2.1077552181385986d+1*t1*t22*t4*t5*t52
     ?         *t8*t9*t91-3.512925363564331d+0*t1*t21*t23*t4*t8*t9*t91+1
     @         .0151243163195242d-1*gammaaa*t1*t19*t20*t23*t8*t9*t91-9.1
     1         21231935850558d-4*gammaaa*t1*t22*t23*t54*t8*t9*t90+4.7229
     2         9650964551d-2*t1*t19*t22*t4*t52*t8*t89*t9-1.3647966030137
     3         373d-3*gammaaa*t1*t22*t52*t8*t87*t88*t9+1.639146470832696
     4         7d-5*t1*t11*t172*t53*t8*t87*t9-1.023597452260303d-3*gamma
     5         aa*t1*t170*t23*t8*t87*t9+7.084494764468265d-2*t1*t19*t37*
     6         t4*t52*t58*t8*t9-1.9033580930991081d-2*t1*t19*t23*t38*t58
     7         *t8*t9+1.4803896279659726d-1*t1*t11*t19*t56*t8*t9-3.80051
     8         3306604399d-3*gammaaa*t1*t11*t54*t55*t8*t9+8.175436355410
     9         144d-6*t1*t168*t169*t23*t53*t8*t9+3.8067161861982163d-2*t
     :         1*t19*t38*t51*t52*t8*t9-2.1253484293404795d-1*t1*t152*t19
     ;         *t37*t4*t51*t8*t9-5.269388045346497d+0*t1*t173*t23*t4*t5*
     <         t8*t9+7.871660849409182d-3*t1*t171*t19*t23*t4*t8*t9+8.882
     =         337767795836d-2*t1*t19*t20*t22*t23*t8*t9)*wght)+t25*t49*(
     >         5.001700820132344d-1*rhoa*rhob*t17*t41*t50*wght+1.0003401
     ?         640264689d+0*rhoa*rhob*t17*t24*t26*t40*wght)+t18*t47*(3.3
     @         344672134215625d-1*rhoa*rhob*t17*t40*t48*wght+6.668934426
     1         843125d-1*rhoa*rhob*t17*t24*t26*t41*wght)+2.0006803280529
     2         378d+0*rhoa*rhob*t150*t16*t17*t40*t50*wght+5.001700820132
     3         344d-1*rhoa*rhob*t14*t149*t17*t41*t48*wght+Cmat3(iq,D3_RA
     4         _RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t14*t25*(8.336168033553906d-2*rhoa*
     1         rhob*t26*t99*wght+8.336168033553906d-2*rhoa*rhob*t17*(-2.
     2         3991838839497967d+1*t1*t102*t22*t23*t4*t5*t8*t9+2.3991838
     3         839497967d+1*t1*t11*t175*t4*t5*t8*t9+7.997279613165988d+0
     4         *t1*t102*t11*t21*t4*t8*t9-2.3109608544299284d-1*gammaaa*t
     5         1*t102*t11*t19*t20*t8*t9-1.612806413056766d-1*t1*t23*t29*
     6         t4*t43*t5*t61*t8+5.254465552928875d-2*t1*t29*t4*t43*t5*t5
     7         2*t58*t8+5.83829505880986d-3*t1*t23*t29*t4*t43*t57*t8-7.5
     8         91871730843472d-4*gammaaa*t1*t19*t23*t29*t43*t56*t8+6.063
     9         600883031948d-6*t1*t23*t29*t43*t53*t54*t55*t8+3.225612826
     :         113532d-1*t1*t22*t24*t29*t4*t43*t5*t52*t8+3.5029770352859
     ;         16d-2*t1*t21*t22*t29*t4*t43*t52*t8-1.0122495641124629d-3*
     <         gammaaa*t1*t19*t20*t22*t29*t43*t52*t8-1.5763396658786624d
     =         -1*t1*t152*t29*t4*t43*t5*t51*t8-5.376021376855886d-2*t1*t
     >         21*t23*t24*t29*t4*t43*t8+1.5535001344756127d-3*gammaaa*t1
     ?         *t19*t20*t23*t24*t29*t43*t8+8.666103204112231d-2*t1*t11*t
     @         29*t4*t42*t5*t61-1.4116926989963413d-2*t1*t23*t29*t4*t42*
     1         t5*t58-3.137094886658535d-3*t1*t11*t29*t4*t42*t57+4.07934
     2         5382014923d-4*gammaaa*t1*t11*t19*t29*t42*t56-3.2581586119
     3         3822d-6*t1*t11*t29*t42*t53*t54*t55+2.8233853979926826d-2*
     4         t1*t29*t4*t42*t5*t51*t52-8.666103204112231d-2*t1*t22*t23*
     5         t24*t29*t4*t42*t5+2.8887010680374103d-2*t1*t11*t21*t24*t2
     6         9*t4*t42-9.411284659975606d-3*t1*t21*t22*t23*t29*t4*t42-8
     7         .347432391127305d-4*gammaaa*t1*t11*t19*t20*t24*t29*t42+2.
     8         7195635880099484d-4*gammaaa*t1*t19*t20*t22*t23*t29*t42)*w
     9         ght+8.336168033553906d-2*rhoa*rhob*t44*t63*wght+8.3361680
     :         33553906d-2*rhob*t26*t44*wght+8.336168033553906d-2*rhoa*r
     ;         hob*t104*t27*wght+8.336168033553906d-2*rhob*t104*t17*wght
     <         -t46*t60-t104*t28-t100*t26)+t16*t18*(1.6672336067107812d-
     =         1*rhoa*rhob*t24*t99*wght+1.6672336067107812d-1*rhoa*rhob*
     >         t44*t61*wght+1.6672336067107812d-1*rhob*t24*t44*wght+1.66
     ?         72336067107812d-1*rhoa*rhob*t102*t27*wght+1.6672336067107
     @         812d-1*rhoa*rhob*t17*t175*wght+1.6672336067107812d-1*rhob
     1         *t102*t17*wght-2.0d+0*t45*t60-2.0d+0*t102*t28-2.0d+0*t100
     2         *t24)+t14*t16*(-1.6672336067107812d-1*rhob*t99*wght-8.336
     3         168033553906d-2*rhoa*rhob*(7.084494764468265d-2*t1*t29*t4
     4         *t43*t5*t52*t58*t8+7.871660849409182d-3*t1*t23*t29*t4*t43
     5         *t57*t8-1.023597452260303d-3*gammaaa*t1*t19*t23*t29*t43*t
     6         56*t8+8.175436355410144d-6*t1*t23*t29*t43*t53*t54*t55*t8+
     7         4.72299650964551d-2*t1*t21*t22*t29*t4*t43*t52*t8-1.364796
     8         6030137373d-3*gammaaa*t1*t19*t20*t22*t29*t43*t52*t8-2.125
     9         3484293404795d-1*t1*t152*t29*t4*t43*t5*t51*t8-1.903358093
     :         0991081d-2*t1*t23*t29*t4*t42*t5*t58-4.22968465133135d-3*t
     ;         1*t11*t29*t4*t42*t57+5.500102857317887d-4*gammaaa*t1*t11*
     <         t19*t29*t42*t56-4.392912541831627d-6*t1*t11*t29*t42*t53*t
     =         54*t55+3.8067161861982163d-2*t1*t29*t4*t42*t5*t51*t52-1.2
     >         689053953994053d-2*t1*t21*t22*t23*t29*t4*t42+3.6667352382
     ?         119245d-4*gammaaa*t1*t19*t20*t22*t23*t29*t42)*wght)+t25*t
     @         49*(5.001700820132344d-1*rhoa*rhob*t17*t46*t50*wght+1.000
     1         3401640264689d+0*rhoa*rhob*t17*t24*t26*t45*wght)+t16*t49*
     2         (-5.001700820132344d-1*rhoa*rhob*t44*t50*wght-1.000340164
     3         0264689d+0*rhoa*rhob*t102*t17*t24*wght-3.0d+0*t45*t62)+t1
     4         8*t47*(3.3344672134215625d-1*rhoa*rhob*t17*t45*t48*wght+6
     5         .668934426843125d-1*rhoa*rhob*t17*t24*t26*t46*wght)+t14*t
     6         47*(-1.6672336067107812d-1*rhoa*rhob*t44*t48*wght-3.33446
     7         72134215625d-1*rhoa*rhob*t104*t17*t26*wght-2*t46*t64)+t18
     8         *t25*(-3.3344672134215625d-1*rhoa*rhob*t24*t26*t44*wght-3
     9         .3344672134215625d-1*rhoa*rhob*t102*t17*t26*wght-3.334467
     :         2134215625d-1*rhoa*rhob*t104*t17*t24*wght-2.0d+0*t45*t64-
     ;         t46*t62)+2.0006803280529378d+0*rhoa*rhob*t150*t16*t17*t45
     <         *t50*wght+5.001700820132344d-1*rhoa*rhob*t14*t149*t17*t46
     =         *t48*wght+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t14*t25*(8.336168033553906d-2*rhoa*
     1         rhob*t35*t97*wght+8.336168033553906d-2*rhoa*t17*t97*wght+
     2         8.336168033553906d-2*rhoa*rhob*t17*(-2.3991838839497967d+
     3         1*t1*t23*t32*t4*t5*t8*t9*t95+7.997279613165988d+0*t1*t11*
     4         t31*t4*t5*t8*t95-2.3109608544299284d-1*gammabb*t1*t11*t29
     5         *t30*t4*t5*t95+7.816455199643731d+0*t1*t32*t4*t5*t52*t8*t
     6         9*t91-1.302742533273955d+0*t1*t23*t31*t4*t5*t8*t91+3.7645
     7         138639902426d-2*gammabb*t1*t23*t29*t30*t4*t5*t91-3.382541
     8         776131776d-4*gammaaa*t1*t23*t32*t54*t8*t9*t90+1.127513925
     9         3772587d-4*gammaaa*t1*t11*t31*t54*t8*t90-3.25815861193822
     :         d-6*gammaaa*gammabb*t1*t11*t29*t30*t54*t90+1.751488517642
     ;         958d-2*t1*t19*t32*t4*t52*t8*t89*t9-5.061247820562315d-4*g
     <         ammaaa*t1*t32*t52*t8*t87*t88*t9-1.612806413056766d-1*t1*t
     =         19*t23*t37*t4*t68*t8*t9+8.666103204112231d-2*t1*t11*t19*t
     >         38*t68*t8*t9+3.225612826113532d-1*t1*t19*t24*t32*t37*t4*t
     ?         52*t8*t9+2.8233853979926826d-2*t1*t19*t22*t32*t38*t52*t8*
     @         t9+2.3991838839497967d+1*t1*t11*t176*t4*t5*t8*t9-1.576339
     1         6658786624d-1*t1*t152*t19*t22*t32*t37*t4*t8*t9-8.66610320
     2         4112231d-2*t1*t19*t23*t24*t32*t38*t8*t9+3.293949630991462
     3         d-2*t1*t19*t20*t23*t32*t8*t9-2.91914752940493d-3*t1*t19*t
     4         23*t31*t4*t8*t89+8.435413034270525d-5*gammabb*t1*t19*t23*
     5         t29*t30*t4*t89+8.435413034270525d-5*gammaaa*t1*t23*t31*t8
     6         *t87*t88-2.437567554978843d-6*gammaaa*gammabb*t1*t23*t29*
     7         t30*t87*t88+1.751488517642958d-2*t1*t19*t22*t31*t37*t4*t5
     8         2*t8-5.376021376855886d-2*t1*t19*t23*t24*t31*t37*t4*t8+2.
     9         8887010680374103d-2*t1*t11*t19*t24*t31*t38*t8-4.705642329
     :         987803d-3*t1*t19*t22*t23*t31*t38*t8-1.0979832103304873d-2
     ;         *t1*t11*t19*t20*t31*t8-5.061247820562315d-4*gammabb*t1*t1
     <         9*t22*t29*t30*t37*t4*t52+1.5535001344756127d-3*gammabb*t1
     =         *t19*t23*t24*t29*t30*t37*t4-8.347432391127305d-4*gammabb*
     >         t1*t11*t19*t24*t29*t30*t38+1.3597817940049742d-4*gammabb*
     ?         t1*t19*t22*t23*t29*t30*t38+3.1728241860116063d-4*gammabb*
     @         t1*t11*t19*t20*t29*t30)*wght+8.336168033553906d-2*rhoa*rh
     1         ob*t39*t70*wght+8.336168033553906d-2*rhoa*t26*t39*wght+8.
     2         336168033553906d-2*rhoa*rhob*t106*t26*wght-t34*t93-t41*t6
     3         6-t111*t28)+t14*t47*(-1.6672336067107812d-1*rhoa*rhob*t17
     4         *t34*t97*wght-1.6672336067107812d-1*rhoa*rhob*t26*t34*t39
     5         *wght-1.6672336067107812d-1*rhoa*rhob*t111*t17*t26*wght-2
     6         *t41*t71)+t18*t25*(-1.6672336067107812d-1*rhoa*rhob*t17*t
     7         33*t97*wght-1.6672336067107812d-1*rhoa*rhob*t17*t34*t95*w
     8         ght-1.6672336067107812d-1*rhoa*rhob*t24*t34*t39*wght-1.66
     9         72336067107812d-1*rhoa*rhob*t26*t33*t39*wght-1.6672336067
     :         107812d-1*rhoa*rhob*t109*t17*t26*wght-1.6672336067107812d
     ;         -1*rhoa*rhob*t111*t17*t24*wght-2.0d+0*t40*t71-t41*t69)+t1
     <         6*t18*(1.6672336067107812d-1*rhoa*rhob*t35*t95*wght+1.667
     =         2336067107812d-1*rhoa*t17*t95*wght+1.6672336067107812d-1*
     >         rhoa*rhob*t39*t68*wght+1.6672336067107812d-1*rhoa*t24*t39
     ?         *wght+1.6672336067107812d-1*rhoa*rhob*t106*t24*wght+1.667
     @         2336067107812d-1*rhoa*rhob*t17*t176*wght-2.0d+0*t33*t93-2
     1         .0d+0*t40*t66-2.0d+0*t109*t28)+t16*t49*(-5.00170082013234
     2         4d-1*rhoa*rhob*t17*t33*t95*wght-5.001700820132344d-1*rhoa
     3         *rhob*t24*t33*t39*wght-5.001700820132344d-1*rhoa*rhob*t10
     4         9*t17*t24*wght-3.0d+0*t40*t69)+t14*t16*(-8.33616803355390
     5         6d-2*rhoa*t92*wght-8.336168033553906d-2*rhoa*rhob*(1.0538
     6         776090692993d+1*t1*t32*t4*t5*t52*t8*t9*t91-1.756462681782
     7         1655d+0*t1*t23*t31*t4*t5*t8*t91+5.075621581597621d-2*gamm
     8         abb*t1*t23*t29*t30*t4*t5*t91-4.560615967925279d-4*gammaaa
     9         *t1*t23*t32*t54*t8*t9*t90+1.5202053226417594d-4*gammaaa*t
     :         1*t11*t31*t54*t8*t90-4.392912541831627d-6*gammaaa*gammabb
     ;         *t1*t11*t29*t30*t54*t90+2.361498254822755d-2*t1*t19*t32*t
     <         4*t52*t8*t89*t9-6.823983015068688d-4*gammaaa*t1*t32*t52*t
     =         8*t87*t88*t9+3.8067161861982163d-2*t1*t19*t22*t32*t38*t52
     >         *t8*t9-2.1253484293404795d-1*t1*t152*t19*t22*t32*t37*t4*t
     ?         8*t9+4.441168883897918d-2*t1*t19*t20*t23*t32*t8*t9-3.9358
     @         30424704591d-3*t1*t19*t23*t31*t4*t8*t89+1.137330502511447
     1         9d-4*gammabb*t1*t19*t23*t29*t30*t4*t89+1.1373305025114479
     2         d-4*gammaaa*t1*t23*t31*t8*t87*t88-3.286525414874878d-6*ga
     3         mmaaa*gammabb*t1*t23*t29*t30*t87*t88+2.361498254822755d-2
     4         *t1*t19*t22*t31*t37*t4*t52*t8-6.344526976997026d-3*t1*t19
     5         *t22*t23*t31*t38*t8-1.4803896279659728d-2*t1*t11*t19*t20*
     6         t31*t8-6.823983015068688d-4*gammabb*t1*t19*t22*t29*t30*t3
     7         7*t4*t52+1.8333676191059622d-4*gammabb*t1*t19*t22*t23*t29
     8         *t30*t38+4.277857777913911d-4*gammabb*t1*t11*t19*t20*t29*
     9         t30)*wght-8.336168033553906d-2*t39*wght-8.336168033553906
     :         d-2*rhob*t106*wght)+t25*t49*(5.001700820132344d-1*rhoa*rh
     ;         ob*t17*t24*t33*t41*wght-3.0d+0*t40*t67)+t18*t47*(3.334467
     <         2134215625d-1*rhoa*rhob*t17*t26*t34*t40*wght-2*t41*t67)+5
     =         .001700820132344d-1*rhoa*rhob*t14*t149*t17*t26*t34*t41*wg
     >         ht+2.0006803280529378d+0*rhoa*rhob*t150*t16*t17*t24*t33*t
     ?         40*wght+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = t14*t16*(-8.336168033553906d-2*rhoa
     1         *t99*wght-8.336168033553906d-2*rhoa*rhob*(1.0538776090692
     2         993d+1*t1*t117*t22*t4*t5*t52*t8*t9-1.7564626817821655d+0*
     3         t1*t117*t21*t23*t4*t8*t9+5.075621581597621d-2*gammaaa*t1*
     4         t117*t19*t20*t23*t8*t9+2.361498254822755d-2*t1*t116*t22*t
     5         29*t4*t5*t52*t8+2.361498254822755d-2*t1*t21*t29*t32*t4*t4
     6         3*t52*t8-6.823983015068688d-4*gammaaa*t1*t19*t20*t29*t32*
     7         t43*t52*t8-2.1253484293404795d-1*t1*t152*t22*t29*t32*t4*t
     8         43*t5*t8-3.935830424704591d-3*t1*t116*t21*t23*t29*t4*t8+1
     9         .1373305025114479d-4*gammaaa*t1*t116*t19*t20*t23*t29*t8-4
     :         .560615967925279d-4*gammabb*t1*t113*t22*t23*t4*t5*t75+1.5
     ;         202053226417594d-4*gammabb*t1*t11*t113*t21*t4*t75-4.39291
     <         2541831627d-6*gammaaa*gammabb*t1*t11*t113*t19*t20*t75+3.8
     =         067161861982163d-2*t1*t22*t29*t32*t4*t42*t5*t52-6.8239830
     >         15068688d-4*gammabb*t1*t114*t115*t22*t4*t5*t52+4.44116888
     ?         3897918d-2*t1*t22*t23*t29*t30*t4*t5-6.344526976997026d-3*
     @         t1*t21*t23*t29*t32*t4*t42+1.8333676191059622d-4*gammaaa*t
     1         1*t19*t20*t23*t29*t32*t42-1.4803896279659728d-2*t1*t11*t2
     2         1*t29*t30*t4+1.1373305025114479d-4*gammabb*t1*t114*t115*t
     3         21*t23*t4+4.277857777913911d-4*gammaaa*t1*t11*t19*t20*t29
     4         *t30-3.286525414874878d-6*gammaaa*gammabb*t1*t114*t115*t1
     5         9*t20*t23)*wght-8.336168033553906d-2*t44*wght-8.336168033
     6         553906d-2*rhob*t118*wght)+t14*t25*(8.336168033553906d-2*r
     7         hoa*rhob*t17*(7.816455199643731d+0*t1*t117*t22*t4*t5*t52*
     8         t8*t9-2.3991838839497967d+1*t1*t102*t23*t32*t4*t5*t8*t9-2
     9         .3991838839497967d+1*t1*t117*t23*t24*t4*t5*t8*t9+2.399183
     :         8839497967d+1*t1*t11*t177*t4*t5*t8*t9-1.302742533273955d+
     ;         0*t1*t117*t21*t23*t4*t8*t9+3.7645138639902426d-2*gammaaa*
     <         t1*t117*t19*t20*t23*t8*t9-1.612806413056766d-1*t1*t23*t29
     =         *t4*t43*t5*t68*t8+3.225612826113532d-1*t1*t24*t29*t32*t4*
     >         t43*t5*t52*t8+1.751488517642958d-2*t1*t116*t22*t29*t4*t5*
     ?         t52*t8+1.751488517642958d-2*t1*t21*t29*t32*t4*t43*t52*t8-
     @         5.061247820562315d-4*gammaaa*t1*t19*t20*t29*t32*t43*t52*t
     1         8-1.5763396658786624d-1*t1*t152*t22*t29*t32*t4*t43*t5*t8+
     2         7.997279613165988d+0*t1*t102*t11*t31*t4*t5*t8-5.376021376
     3         855886d-2*t1*t116*t23*t24*t29*t4*t5*t8-2.91914752940493d-
     4         3*t1*t116*t21*t23*t29*t4*t8+8.435413034270525d-5*gammaaa*
     5         t1*t116*t19*t20*t23*t29*t8+2.0764757191858965d-3*gammabb*
     6         t1*t11*t113*t24*t4*t5*t75-3.382541776131776d-4*gammabb*t1
     7         *t113*t22*t23*t4*t5*t75+1.1275139253772587d-4*gammabb*t1*
     8         t11*t113*t21*t4*t75-3.25815861193822d-6*gammaaa*gammabb*t
     9         1*t11*t113*t19*t20*t75+8.666103204112231d-2*t1*t11*t29*t4
     :         *t42*t5*t68+2.8233853979926826d-2*t1*t22*t29*t32*t4*t42*t
     ;         5*t52-5.061247820562315d-4*gammabb*t1*t114*t115*t22*t4*t5
     <         *t52-8.666103204112231d-2*t1*t23*t24*t29*t32*t4*t42*t5-2.
     =         0220907476261873d-1*t1*t11*t24*t29*t30*t4*t5+3.2939496309
     >         91462d-2*t1*t22*t23*t29*t30*t4*t5-2.3109608544299284d-1*g
     ?         ammabb*t1*t102*t11*t29*t30*t4*t5+1.5535001344756127d-3*ga
     @         mmabb*t1*t114*t115*t23*t24*t4*t5-4.705642329987803d-3*t1*
     1         t21*t23*t29*t32*t4*t42+1.3597817940049742d-4*gammaaa*t1*t
     2         19*t20*t23*t29*t32*t42-1.0979832103304873d-2*t1*t11*t21*t
     3         29*t30*t4+8.435413034270525d-5*gammabb*t1*t114*t115*t21*t
     4         23*t4+3.1728241860116063d-4*gammaaa*t1*t11*t19*t20*t29*t3
     5         0-2.437567554978843d-6*gammaaa*gammabb*t1*t114*t115*t19*t
     6         20*t23)*wght+8.336168033553906d-2*rhoa*rhob*t44*t70*wght+
     7         8.336168033553906d-2*rhoa*t26*t44*wght+8.336168033553906d
     8         -2*rhoa*rhob*t104*t35*wght+8.336168033553906d-2*rhoa*rhob
     9         *t118*t26*wght+8.336168033553906d-2*rhoa*t104*t17*wght-t4
     :         6*t66-t100*t34-t123*t28)+t16*t18*(1.6672336067107812d-1*r
     ;         hoa*rhob*t44*t68*wght+1.6672336067107812d-1*rhoa*t24*t44*
     <         wght+1.6672336067107812d-1*rhoa*rhob*t102*t35*wght+1.6672
     =         336067107812d-1*rhoa*rhob*t118*t24*wght+1.667233606710781
     >         2d-1*rhoa*rhob*t17*t177*wght+1.6672336067107812d-1*rhoa*t
     ?         102*t17*wght-2.0d+0*t45*t66-2.0d+0*t100*t33-2.0d+0*t121*t
     @         28)+t25*t49*(5.001700820132344d-1*rhoa*rhob*t17*t24*t33*t
     1         46*wght-3.0d+0*t45*t67)+t18*t47*(3.3344672134215625d-1*rh
     2         oa*rhob*t17*t26*t34*t45*wght-2*t46*t67)+t14*t47*(-1.66723
     3         36067107812d-1*rhoa*rhob*t26*t34*t44*wght-1.6672336067107
     4         812d-1*rhoa*rhob*t104*t17*t34*wght-1.6672336067107812d-1*
     5         rhoa*rhob*t123*t17*t26*wght-2*t46*t71)+t18*t25*(-1.667233
     6         6067107812d-1*rhoa*rhob*t24*t34*t44*wght-1.66723360671078
     7         12d-1*rhoa*rhob*t26*t33*t44*wght-1.6672336067107812d-1*rh
     8         oa*rhob*t102*t17*t34*wght-1.6672336067107812d-1*rhoa*rhob
     9         *t104*t17*t33*wght-1.6672336067107812d-1*rhoa*rhob*t121*t
     :         17*t26*wght-1.6672336067107812d-1*rhoa*rhob*t123*t17*t24*
     ;         wght-2.0d+0*t45*t71-t46*t69)+t16*t49*(-5.001700820132344d
     <         -1*rhoa*rhob*t24*t33*t44*wght-5.001700820132344d-1*rhoa*r
     =         hob*t102*t17*t33*wght-5.001700820132344d-1*rhoa*rhob*t121
     >         *t17*t24*wght-3.0d+0*t45*t69)+5.001700820132344d-1*rhoa*r
     ?         hob*t14*t149*t17*t26*t34*t46*wght+2.0006803280529378d+0*r
     @         hoa*rhob*t150*t16*t17*t24*t33*t45*wght+Cmat3(iq,D3_RA_RB_
     1         GBB)
            Cmat3(iq,D3_RB_RB_GAA) = t14*t25*(8.336168033553906d-2*rhoa*
     1         rhob*t17*(-1.612806413056766d-1*t1*t19*t23*t37*t4*t8*t83*
     2         t9+8.666103204112231d-2*t1*t11*t19*t38*t8*t83*t9+5.254465
     3         552928875d-2*t1*t19*t37*t4*t52*t8*t80*t9-1.41169269899634
     4         13d-2*t1*t19*t23*t38*t8*t80*t9+2.8233853979926826d-2*t1*t
     5         19*t38*t52*t79*t8*t9-1.5763396658786624d-1*t1*t152*t19*t3
     6         7*t4*t79*t8*t9+3.225612826113532d-1*t1*t19*t32*t33*t37*t4
     7         *t52*t8*t9-2.3991838839497967d+1*t1*t109*t23*t32*t4*t5*t8
     8         *t9+2.3991838839497967d+1*t1*t11*t178*t4*t5*t8*t9-8.66610
     9         3204112231d-2*t1*t19*t23*t32*t33*t38*t8*t9+5.838295058809
     :         86d-3*t1*t19*t23*t37*t4*t78*t8-3.137094886658535d-3*t1*t1
     ;         1*t19*t38*t78*t8+3.502977035285916d-2*t1*t19*t31*t32*t37*
     <         t4*t52*t8+7.997279613165988d+0*t1*t109*t11*t31*t4*t5*t8-5
     =         .376021376855886d-2*t1*t19*t23*t31*t33*t37*t4*t8+2.888701
     >         0680374103d-2*t1*t11*t19*t31*t33*t38*t8-9.411284659975606
     ?         d-3*t1*t19*t23*t31*t32*t38*t8-7.591871730843472d-4*gammab
     @         b*t1*t19*t23*t29*t37*t4*t77+4.079345382014923d-4*gammabb*
     1         t1*t11*t19*t29*t38*t77+6.063600883031948d-6*t1*t19*t23*t3
     2         7*t4*t74*t75*t76-3.25815861193822d-6*t1*t11*t19*t38*t74*t
     3         75*t76-1.0122495641124629d-3*gammabb*t1*t19*t29*t30*t32*t
     4         37*t4*t52-2.3109608544299284d-1*gammabb*t1*t109*t11*t29*t
     5         30*t4*t5+1.5535001344756127d-3*gammabb*t1*t19*t23*t29*t30
     6         *t33*t37*t4-8.347432391127305d-4*gammabb*t1*t11*t19*t29*t
     7         30*t33*t38+2.7195635880099484d-4*gammabb*t1*t19*t23*t29*t
     8         30*t32*t38)*wght+8.336168033553906d-2*rhoa*rhob*t39*t85*w
     9         ght+8.336168033553906d-2*rhoa*t34*t39*wght+8.336168033553
     :         906d-2*rhoa*rhob*t111*t35*wght+8.336168033553906d-2*rhoa*
     ;         rhob*t106*t34*wght+8.336168033553906d-2*rhoa*t111*t17*wgh
     <         t-t41*t82-t111*t36-t107*t34)+t14*t16*(-8.336168033553906d
     =         -2*rhoa*rhob*(7.084494764468265d-2*t1*t19*t37*t4*t52*t8*t
     >         80*t9-1.9033580930991081d-2*t1*t19*t23*t38*t8*t80*t9+3.80
     ?         67161861982163d-2*t1*t19*t38*t52*t79*t8*t9-2.125348429340
     @         4795d-1*t1*t152*t19*t37*t4*t79*t8*t9+7.871660849409182d-3
     1         *t1*t19*t23*t37*t4*t78*t8-4.22968465133135d-3*t1*t11*t19*
     2         t38*t78*t8+4.72299650964551d-2*t1*t19*t31*t32*t37*t4*t52*
     3         t8-1.2689053953994053d-2*t1*t19*t23*t31*t32*t38*t8-1.0235
     4         97452260303d-3*gammabb*t1*t19*t23*t29*t37*t4*t77+5.500102
     5         857317887d-4*gammabb*t1*t11*t19*t29*t38*t77+8.17543635541
     6         0144d-6*t1*t19*t23*t37*t4*t74*t75*t76-4.392912541831627d-
     7         6*t1*t11*t19*t38*t74*t75*t76-1.3647966030137373d-3*gammab
     8         b*t1*t19*t29*t30*t32*t37*t4*t52+3.6667352382119245d-4*gam
     9         mabb*t1*t19*t23*t29*t30*t32*t38)*wght-1.6672336067107812d
     :         -1*rhoa*t106*wght)+t16*t18*(1.6672336067107812d-1*rhoa*rh
     ;         ob*t39*t83*wght+1.6672336067107812d-1*rhoa*t33*t39*wght+1
     <         .6672336067107812d-1*rhoa*rhob*t109*t35*wght+1.6672336067
     =         107812d-1*rhoa*rhob*t106*t33*wght+1.6672336067107812d-1*r
     >         hoa*rhob*t17*t178*wght+1.6672336067107812d-1*rhoa*t109*t1
     ?         7*wght-2.0d+0*t40*t82-2.0d+0*t109*t36-2.0d+0*t107*t33)+t2
     @         5*t49*(5.001700820132344d-1*rhoa*rhob*t17*t41*t73*wght+1.
     1         0003401640264689d+0*rhoa*rhob*t17*t33*t34*t40*wght)+t16*t
     2         49*(-5.001700820132344d-1*rhoa*rhob*t39*t73*wght-1.000340
     3         1640264689d+0*rhoa*rhob*t109*t17*t33*wght-3.0d+0*t40*t84)
     4         +t18*t47*(3.3344672134215625d-1*rhoa*rhob*t17*t40*t72*wgh
     5         t+6.668934426843125d-1*rhoa*rhob*t17*t33*t34*t41*wght)+t1
     6         4*t47*(-1.6672336067107812d-1*rhoa*rhob*t39*t72*wght-3.33
     7         44672134215625d-1*rhoa*rhob*t111*t17*t34*wght-2*t41*t86)+
     8         t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t33*t34*t39*wgh
     9         t-3.3344672134215625d-1*rhoa*rhob*t109*t17*t34*wght-3.334
     :         4672134215625d-1*rhoa*rhob*t111*t17*t33*wght-2.0d+0*t40*t
     ;         86-t41*t84)+2.0006803280529378d+0*rhoa*rhob*t150*t16*t17*
     <         t40*t73*wght+5.001700820132344d-1*rhoa*rhob*t14*t149*t17*
     =         t41*t72*wght+Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = t14*t16*(-8.336168033553906d-2*rhoa
     1         *rhob*(2.1077552181385986d+1*t1*t117*t32*t4*t5*t52*t8*t9-
     2         5.269388045346497d+0*t1*t184*t23*t4*t5*t8*t9+7.0844947644
     3         68265d-2*t1*t29*t4*t43*t5*t52*t8*t80-1.9033580930991081d-
     4         2*t1*t23*t29*t4*t42*t5*t80-2.1253484293404795d-1*t1*t152*
     5         t29*t4*t43*t5*t79*t8+4.72299650964551d-2*t1*t116*t29*t32*
     6         t4*t5*t52*t8-3.512925363564331d+0*t1*t117*t23*t31*t4*t5*t
     7         8+7.871660849409182d-3*t1*t183*t23*t29*t4*t5*t8+3.8067161
     8         861982163d-2*t1*t29*t4*t42*t5*t52*t79+1.4803896279659726d
     9         -1*t1*t11*t29*t4*t5*t77-3.800513306604399d-3*gammabb*t1*t
     :         11*t4*t5*t75*t76-9.121231935850558d-4*gammabb*t1*t113*t23
     ;         *t32*t4*t5*t75+8.175436355410144d-6*t1*t180*t181*t23*t4*t
     <         5*t74+1.6391464708326967d-5*t1*t11*t114*t179*t4*t5*t74-1.
     =         3647966030137373d-3*gammabb*t1*t114*t115*t32*t4*t5*t52+8.
     >         882337767795836d-2*t1*t23*t29*t30*t32*t4*t5+1.01512431631
     ?         95242d-1*gammabb*t1*t117*t23*t29*t30*t4*t5-1.023597452260
     @         303d-3*gammabb*t1*t114*t182*t23*t4*t5)*wght-1.66723360671
     1         07812d-1*rhoa*t118*wght)+t14*t25*(8.336168033553906d-2*rh
     2         oa*rhob*t17*(1.563291039928746d+1*t1*t117*t32*t4*t5*t52*t
     3         8*t9-2.3991838839497967d+1*t1*t117*t23*t33*t4*t5*t8*t9-2.
     4         3991838839497967d+1*t1*t121*t23*t32*t4*t5*t8*t9-3.9082275
     5         99821866d+0*t1*t184*t23*t4*t5*t8*t9+2.3991838839497967d+1
     6         *t1*t11*t185*t4*t5*t8*t9-1.612806413056766d-1*t1*t23*t29*
     7         t4*t43*t5*t8*t83+8.666103204112231d-2*t1*t11*t29*t4*t42*t
     8         5*t83+5.254465552928875d-2*t1*t29*t4*t43*t5*t52*t8*t80-1.
     9         4116926989963413d-2*t1*t23*t29*t4*t42*t5*t80-1.5763396658
     :         786624d-1*t1*t152*t29*t4*t43*t5*t79*t8+3.225612826113532d
     ;         -1*t1*t29*t32*t33*t4*t43*t5*t52*t8+3.502977035285916d-2*t
     <         1*t116*t29*t32*t4*t5*t52*t8-5.376021376855886d-2*t1*t116*
     =         t23*t29*t33*t4*t5*t8-2.60548506654791d+0*t1*t117*t23*t31*
     >         t4*t5*t8+7.997279613165988d+0*t1*t11*t121*t31*t4*t5*t8+5.
     ?         83829505880986d-3*t1*t183*t23*t29*t4*t5*t8+2.823385397992
     @         6826d-2*t1*t29*t4*t42*t5*t52*t79+1.0979832103304875d-1*t1
     1         *t11*t29*t4*t5*t77-2.8187848134431465d-3*gammabb*t1*t11*t
     2         4*t5*t75*t76+2.0764757191858965d-3*gammabb*t1*t11*t113*t3
     3         3*t4*t5*t75-6.765083552263552d-4*gammabb*t1*t113*t23*t32*
     4         t4*t5*t75+6.063600883031948d-6*t1*t180*t181*t23*t4*t5*t74
     5         +1.2157308253500822d-5*t1*t11*t114*t179*t4*t5*t74-1.01224
     6         95641124629d-3*gammabb*t1*t114*t115*t32*t4*t5*t52-8.66610
     7         3204112231d-2*t1*t23*t29*t32*t33*t4*t42*t5-2.022090747626
     8         1873d-1*t1*t11*t29*t30*t33*t4*t5+1.5535001344756127d-3*ga
     9         mmabb*t1*t114*t115*t23*t33*t4*t5+6.587899261982924d-2*t1*
     :         t23*t29*t30*t32*t4*t5+7.529027727980485d-2*gammabb*t1*t11
     ;         7*t23*t29*t30*t4*t5-2.3109608544299284d-1*gammabb*t1*t11*
     <         t121*t29*t30*t4*t5-7.591871730843472d-4*gammabb*t1*t114*t
     =         182*t23*t4*t5)*wght+8.336168033553906d-2*rhoa*rhob*t44*t8
     >         5*wght+8.336168033553906d-2*rhoa*t34*t44*wght+8.336168033
     ?         553906d-2*rhoa*rhob*t123*t35*wght+8.336168033553906d-2*rh
     @         oa*rhob*t118*t34*wght+8.336168033553906d-2*rhoa*t123*t17*
     1         wght-t46*t82-t123*t36-t119*t34)+t16*t18*(1.66723360671078
     2         12d-1*rhoa*rhob*t44*t83*wght+1.6672336067107812d-1*rhoa*t
     3         33*t44*wght+1.6672336067107812d-1*rhoa*rhob*t121*t35*wght
     4         +1.6672336067107812d-1*rhoa*rhob*t118*t33*wght+1.66723360
     5         67107812d-1*rhoa*rhob*t17*t185*wght+1.6672336067107812d-1
     6         *rhoa*t121*t17*wght-2.0d+0*t45*t82-2.0d+0*t121*t36-2.0d+0
     7         *t119*t33)+t25*t49*(5.001700820132344d-1*rhoa*rhob*t17*t4
     8         6*t73*wght+1.0003401640264689d+0*rhoa*rhob*t17*t33*t34*t4
     9         5*wght)+t16*t49*(-5.001700820132344d-1*rhoa*rhob*t44*t73*
     :         wght-1.0003401640264689d+0*rhoa*rhob*t121*t17*t33*wght-3.
     ;         0d+0*t45*t84)+t18*t47*(3.3344672134215625d-1*rhoa*rhob*t1
     <         7*t45*t72*wght+6.668934426843125d-1*rhoa*rhob*t17*t33*t34
     =         *t46*wght)+t14*t47*(-1.6672336067107812d-1*rhoa*rhob*t44*
     >         t72*wght-3.3344672134215625d-1*rhoa*rhob*t123*t17*t34*wgh
     ?         t-2*t46*t86)+t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t3
     @         3*t34*t44*wght-3.3344672134215625d-1*rhoa*rhob*t121*t17*t
     1         34*wght-3.3344672134215625d-1*rhoa*rhob*t123*t17*t33*wght
     2         -2.0d+0*t45*t86-t46*t84)+2.0006803280529378d+0*rhoa*rhob*
     3         t150*t16*t17*t45*t73*wght+5.001700820132344d-1*rhoa*rhob*
     4         t14*t149*t17*t46*t72*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = t14*t47*(-1.6672336067107812d-1*rh
     1         oa*rhob*t17*t41*t97*wght-1.6672336067107812d-1*rhoa*rhob*
     2         t26*t39*t41*wght-1.6672336067107812d-1*rhoa*rhob*t132*t17
     3         *t26*wght-2*t41*t98)+t18*t25*(-1.6672336067107812d-1*rhoa
     4         *rhob*t17*t40*t97*wght-1.6672336067107812d-1*rhoa*rhob*t1
     5         7*t41*t95*wght-1.6672336067107812d-1*rhoa*rhob*t24*t39*t4
     6         1*wght-1.6672336067107812d-1*rhoa*rhob*t26*t39*t40*wght-1
     7         .6672336067107812d-1*rhoa*rhob*t130*t17*t26*wght-1.667233
     8         6067107812d-1*rhoa*rhob*t132*t17*t24*wght-2.0d+0*t40*t98-
     9         t41*t96)+t14*t25*(1.6672336067107812d-1*rhoa*rhob*t39*t97
     :         *wght+8.336168033553906d-2*rhoa*rhob*t17*(-3.225612826113
     ;         532d-1*t1*t19*t23*t37*t4*t8*t9*t95+1.7332206408224463d-1*
     <         t1*t11*t19*t38*t8*t9*t95+1.050893110585775d-1*t1*t19*t37*
     =         t4*t52*t8*t9*t91-2.8233853979926826d-2*t1*t19*t23*t38*t8*
     >         t9*t91+6.34226583024708d-4*t1*t11*t54*t8*t9*t90+4.4285918
     ?         42992025d-4*t1*t23*t8*t87*t88*t9+7.868855442416534d-5*t1*
     @         t23*t4*t54*t8*t88*t9+2.168357784827018d-3*t1*t24*t4*t52*t
     1         56*t8*t87*t9-1.0596648048916558d-3*t1*t152*t22*t4*t56*t8*
     2         t87*t9+1.1774053387685064d-4*t1*t156*t4*t52*t8*t87*t9+3.7
     3         95935865421737d-4*t1*t125*t22*t52*t8*t87*t9-1.16512510085
     4         67096d-3*t1*t125*t23*t24*t8*t87*t9-4.558990595062808d-6*g
     5         ammaaa*t1*t11*t188*t8*t87*t9-4.7213132654499207d-4*t1*t12
     6         5*t22*t4*t52*t54*t8*t9+1.4491605732048626d-3*t1*t125*t23*
     7         t24*t4*t54*t8*t9-7.786783946947111d-4*t1*t11*t126*t24*t54
     8         *t8*t9+1.268453166049416d-4*t1*t126*t22*t23*t54*t8*t9-3.4
     9         02329015996045d-6*gammaaa*t1*t155*t186*t52*t8*t9-3.908227
     :         599821866d+0*t1*t189*t23*t4*t5*t8*t9+2.3991838839497967d+
     ;         1*t1*t11*t190*t4*t5*t8*t9-6.821550993410941d-6*gammaaa*t1
     <         *t168*t187*t23*t8*t9)*wght+8.336168033553906d-2*rhoa*rhob
     =         *t127*t26*wght-2*t41*t93-t132*t28)+t16*t49*(-5.0017008201
     >         32344d-1*rhoa*rhob*t17*t40*t95*wght-5.001700820132344d-1*
     ?         rhoa*rhob*t24*t39*t40*wght-5.001700820132344d-1*rhoa*rhob
     @         *t130*t17*t24*wght-3.0d+0*t40*t96)+t16*t18*(3.33446721342
     1         15625d-1*rhoa*rhob*t39*t95*wght+1.6672336067107812d-1*rho
     2         a*rhob*t127*t24*wght+1.6672336067107812d-1*rhoa*rhob*t17*
     3         t190*wght-4.0d+0*t40*t93-2.0d+0*t130*t28)+t14*t16*(-8.336
     4         168033553906d-2*rhoa*rhob*(1.4168989528936532d-1*t1*t19*t
     5         37*t4*t52*t8*t9*t91-3.8067161861982163d-2*t1*t19*t23*t38*
     6         t8*t9*t91+8.551154939859897d-4*t1*t11*t54*t8*t9*t90+5.970
     7         985138185101d-4*t1*t23*t8*t87*t88*t9+1.0609426329397834d-
     8         4*t1*t23*t4*t54*t8*t88*t9-1.4287256594843722d-3*t1*t152*t
     9         22*t4*t56*t8*t87*t9+1.5874729549826355d-4*t1*t156*t4*t52*
     :         t8*t87*t9+5.117987261301516d-4*t1*t125*t22*t52*t8*t87*t9-
     ;         6.146799265622613d-6*gammaaa*t1*t11*t188*t8*t87*t9-6.3656
     <         55797638701d-4*t1*t125*t22*t4*t52*t54*t8*t9+1.71023098797
     =         198d-4*t1*t126*t22*t23*t54*t8*t9-4.587294722559716d-6*gam
     >         maaa*t1*t155*t186*t52*t8*t9-5.269388045346497d+0*t1*t189*
     ?         t23*t4*t5*t8*t9-9.197365899836413d-6*gammaaa*t1*t168*t187
     @         *t23*t8*t9)*wght-8.336168033553906d-2*rhob*t127*wght)+t18
     1         *t47*(3.3344672134215625d-1*rhoa*rhob*t17*t26*t40*t41*wgh
     2         t-2*t41*t94)+t25*t49*(5.001700820132344d-1*rhoa*rhob*t17*
     3         t24*t40*t41*wght-3.0d+0*t40*t94)+5.001700820132344d-1*rho
     4         a*rhob*t128*t14*t149*t17*t26*wght+2.0006803280529378d+0*r
     5         hoa*rhob*t129*t150*t16*t17*t24*wght+Cmat3(iq,D3_RA_GAA_GA
     6         A)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t14*t25*(8.336168033553906d-2*rhoa
     1         *rhob*t44*t97*wght+8.336168033553906d-2*rhoa*rhob*t17*(-1
     2         .612806413056766d-1*t1*t23*t29*t4*t43*t5*t8*t95+8.6661032
     3         04112231d-2*t1*t11*t29*t4*t42*t5*t95+5.254465552928875d-2
     4         *t1*t29*t4*t43*t5*t52*t8*t91-1.4116926989963413d-2*t1*t23
     5         *t29*t4*t42*t5*t91-2.2738503311369804d-6*gammaaa*t1*t23*t
     6         29*t43*t54*t8*t90+1.2218094794768326d-6*gammaaa*t1*t11*t2
     7         9*t42*t54*t90+2.3991838839497967d+1*t1*t11*t191*t4*t5*t8*
     8         t9-1.612806413056766d-1*t1*t102*t19*t23*t37*t4*t8*t9+8.66
     9         6103204112231d-2*t1*t102*t11*t19*t38*t8*t9+1.177405338768
     :         5064d-4*t1*t19*t29*t4*t43*t52*t8*t89-3.163279887851447d-5
     ;         *t1*t19*t23*t29*t4*t42*t89-3.402329015996045d-6*gammaaa*t
     <         1*t29*t43*t52*t8*t87*t88+9.140878331170662d-7*gammaaa*t1*
     =         t23*t29*t42*t87*t88+2.168357784827018d-3*t1*t19*t24*t29*t
     >         37*t4*t43*t52*t8+1.8979679327108684d-4*t1*t19*t22*t29*t38
     ?         *t43*t52*t8-1.0596648048916558d-3*t1*t152*t19*t22*t29*t37
     @         *t4*t43*t8-5.825625504283547d-4*t1*t19*t23*t24*t29*t38*t4
     1         3*t8+2.2142959214960128d-4*t1*t19*t20*t23*t29*t43*t8+1.89
     2         79679327108684d-4*t1*t19*t22*t29*t37*t4*t42*t52-5.8256255
     3         04283548d-4*t1*t19*t23*t24*t29*t37*t4*t42+3.1302871466727
     4         39d-4*t1*t11*t19*t24*t29*t38*t42-5.099181727518654d-5*t1*
     5         t19*t22*t23*t29*t38*t42-1.1898090697543523d-4*t1*t11*t19*
     6         t20*t29*t42)*wght+8.336168033553906d-2*rhoa*rhob*t104*t39
     7         *wght+8.336168033553906d-2*rhoa*rhob*t134*t26*wght-t46*t9
     8         3-t100*t41-t137*t28)+t16*t18*(1.6672336067107812d-1*rhoa*
     9         rhob*t44*t95*wght+1.6672336067107812d-1*rhoa*rhob*t102*t3
     :         9*wght+1.6672336067107812d-1*rhoa*rhob*t134*t24*wght+1.66
     ;         72336067107812d-1*rhoa*rhob*t17*t191*wght-2.0d+0*t45*t93-
     <         2.0d+0*t100*t40-2.0d+0*t135*t28)+t14*t16*(-8.336168033553
     =         906d-2*rhoa*rhob*(7.084494764468265d-2*t1*t29*t4*t43*t5*t
     >         52*t8*t91-1.9033580930991081d-2*t1*t23*t29*t4*t42*t5*t91-
     ?         3.0657886332788037d-6*gammaaa*t1*t23*t29*t43*t54*t8*t90+1
     @         .6473422031868606d-6*gammaaa*t1*t11*t29*t42*t54*t90+1.587
     1         4729549826355d-4*t1*t19*t29*t4*t43*t52*t8*t89-4.264989384
     2         4179293d-5*t1*t19*t23*t29*t4*t42*t89-4.587294722559716d-6
     3         *gammaaa*t1*t29*t43*t52*t8*t87*t88+1.2324470305780792d-6*
     4         gammaaa*t1*t23*t29*t42*t87*t88+2.558993630650758d-4*t1*t1
     5         9*t22*t29*t38*t43*t52*t8-1.4287256594843722d-3*t1*t152*t1
     6         9*t22*t29*t37*t4*t43*t8+2.9854925690925505d-4*t1*t19*t20*
     7         t23*t29*t43*t8+2.558993630650758d-4*t1*t19*t22*t29*t37*t4
     8         *t42*t52-6.875128571647359d-5*t1*t19*t22*t23*t29*t38*t42-
     9         1.6041966667177168d-4*t1*t11*t19*t20*t29*t42)*wght-8.3361
     :         68033553906d-2*rhob*t134*wght)+t25*t49*(5.001700820132344
     ;         d-1*rhoa*rhob*t17*t24*t40*t46*wght-3.0d+0*t45*t94)+t18*t4
     <         7*(3.3344672134215625d-1*rhoa*rhob*t17*t26*t41*t45*wght-2
     =         *t46*t94)+t14*t47*(-1.6672336067107812d-1*rhoa*rhob*t26*t
     >         41*t44*wght-1.6672336067107812d-1*rhoa*rhob*t104*t17*t41*
     ?         wght-1.6672336067107812d-1*rhoa*rhob*t137*t17*t26*wght-2*
     @         t46*t98)+t18*t25*(-1.6672336067107812d-1*rhoa*rhob*t24*t4
     1         1*t44*wght-1.6672336067107812d-1*rhoa*rhob*t26*t40*t44*wg
     2         ht-1.6672336067107812d-1*rhoa*rhob*t102*t17*t41*wght-1.66
     3         72336067107812d-1*rhoa*rhob*t104*t17*t40*wght-1.667233606
     4         7107812d-1*rhoa*rhob*t135*t17*t26*wght-1.6672336067107812
     5         d-1*rhoa*rhob*t137*t17*t24*wght-2.0d+0*t45*t98-t46*t96)+t
     6         16*t49*(-5.001700820132344d-1*rhoa*rhob*t24*t40*t44*wght-
     7         5.001700820132344d-1*rhoa*rhob*t102*t17*t40*wght-5.001700
     8         820132344d-1*rhoa*rhob*t135*t17*t24*wght-3.0d+0*t45*t96)+
     9         5.001700820132344d-1*rhoa*rhob*t14*t149*t17*t26*t41*t46*w
     :         ght+2.0006803280529378d+0*rhoa*rhob*t150*t16*t17*t24*t40*
     ;         t45*wght+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t14*t25*(8.336168033553906d-2*rhoa
     1         *rhob*t17*(2.3991838839497967d+1*t1*t11*t192*t4*t5*t8*t9+
     2         2.168357784827018d-3*t1*t114*t24*t4*t5*t52*t77*t8+1.17740
     3         53387685064d-4*t1*t114*t21*t4*t52*t77*t8-3.40232901599604
     4         5d-6*gammaaa*t1*t114*t19*t20*t52*t77*t8-1.059664804891655
     5         8d-3*t1*t114*t152*t22*t4*t5*t77*t8-4.7213132654499207d-4*
     6         t1*t141*t22*t4*t5*t52*t75*t8+1.4491605732048626d-3*t1*t14
     7         1*t23*t24*t4*t5*t75*t8+7.868855442416534d-5*t1*t141*t21*t
     8         23*t4*t75*t8-2.2738503311369804d-6*gammaaa*t1*t141*t19*t2
     9         0*t23*t75*t8-3.225612826113532d-1*t1*t102*t23*t29*t4*t43*
     :         t5*t8-7.786783946947111d-4*t1*t11*t140*t24*t4*t5*t75+1.26
     ;         8453166049416d-4*t1*t140*t22*t23*t4*t5*t75-4.228177220164
     <         72d-5*t1*t11*t140*t21*t4*t75+1.2218094794768326d-6*gammaa
     =         a*t1*t11*t140*t19*t20*t75+3.795935865421737d-4*t1*t114*t1
     >         41*t22*t4*t5*t52+1.7332206408224463d-1*t1*t102*t11*t29*t4
     ?         *t42*t5-1.1651251008567096d-3*t1*t114*t141*t23*t24*t4*t5-
     @         6.326559775702894d-5*t1*t114*t141*t21*t23*t4+1.8281756662
     1         341322d-6*gammaaa*t1*t114*t141*t19*t20*t23)*wght+1.667233
     2         6067107812d-1*rhoa*rhob*t104*t44*wght+8.336168033553906d-
     3         2*rhoa*rhob*t142*t26*wght-2*t100*t46-t147*t28)+t14*t16*(-
     4         8.336168033553906d-2*rhoa*rhob*(1.5874729549826355d-4*t1*
     5         t114*t21*t4*t52*t77*t8-4.587294722559716d-6*gammaaa*t1*t1
     6         14*t19*t20*t52*t77*t8-1.4287256594843722d-3*t1*t114*t152*
     7         t22*t4*t5*t77*t8-6.365655797638701d-4*t1*t141*t22*t4*t5*t
     8         52*t75*t8+1.0609426329397834d-4*t1*t141*t21*t23*t4*t75*t8
     9         -3.0657886332788037d-6*gammaaa*t1*t141*t19*t20*t23*t75*t8
     :         +1.71023098797198d-4*t1*t140*t22*t23*t4*t5*t75-5.70076995
     ;         9906599d-5*t1*t11*t140*t21*t4*t75+1.6473422031868606d-6*g
     <         ammaaa*t1*t11*t140*t19*t20*t75+5.117987261301516d-4*t1*t1
     =         14*t141*t22*t4*t5*t52-8.529978768835859d-5*t1*t114*t141*t
     >         21*t23*t4+2.4648940611561584d-6*gammaaa*t1*t114*t141*t19*
     ?         t20*t23)*wght-8.336168033553906d-2*rhob*t142*wght)+t18*t4
     @         7*(3.3344672134215625d-1*rhoa*rhob*t17*t26*t45*t46*wght-2
     1         *t101*t46)+t25*t49*(5.001700820132344d-1*rhoa*rhob*t17*t2
     2         4*t45*t46*wght-3.0d+0*t101*t45)+t14*t47*(-1.6672336067107
     3         812d-1*rhoa*rhob*t26*t44*t46*wght-1.6672336067107812d-1*r
     4         hoa*rhob*t104*t17*t46*wght-1.6672336067107812d-1*rhoa*rho
     5         b*t147*t17*t26*wght-2*t105*t46)+t18*t25*(-1.6672336067107
     6         812d-1*rhoa*rhob*t24*t44*t46*wght-1.6672336067107812d-1*r
     7         hoa*rhob*t102*t17*t46*wght-1.6672336067107812d-1*rhoa*rho
     8         b*t26*t44*t45*wght-1.6672336067107812d-1*rhoa*rhob*t104*t
     9         17*t45*wght-1.6672336067107812d-1*rhoa*rhob*t145*t17*t26*
     :         wght-1.6672336067107812d-1*rhoa*rhob*t147*t17*t24*wght-t1
     ;         03*t46-2.0d+0*t105*t45)+t16*t49*(-5.001700820132344d-1*rh
     <         oa*rhob*t24*t44*t45*wght-5.001700820132344d-1*rhoa*rhob*t
     =         102*t17*t45*wght-5.001700820132344d-1*rhoa*rhob*t145*t17*
     >         t24*wght-3.0d+0*t103*t45)+t16*t18*(3.3344672134215625d-1*
     ?         rhoa*rhob*t102*t44*wght+1.6672336067107812d-1*rhoa*rhob*t
     @         142*t24*wght+1.6672336067107812d-1*rhoa*rhob*t17*t192*wgh
     1         t-4.0d+0*t100*t45-2.0d+0*t145*t28)+5.001700820132344d-1*r
     2         hoa*rhob*t14*t143*t149*t17*t26*wght+2.0006803280529378d+0
     3         *rhoa*rhob*t144*t150*t16*t17*t24*wght+Cmat3(iq,D3_RA_GBB_
     4         GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = t14*t25*(8.336168033553906d-2*rhoa
     1         *rhob*t17*(2.168357784827018d-3*t1*t33*t4*t52*t56*t8*t87*
     2         t9-1.0596648048916558d-3*t1*t152*t32*t4*t56*t8*t87*t9+3.7
     3         95935865421737d-4*t1*t125*t32*t52*t8*t87*t9-1.16512510085
     4         67096d-3*t1*t125*t23*t33*t8*t87*t9-4.7213132654499207d-4*
     5         t1*t125*t32*t4*t52*t54*t8*t9+1.4491605732048626d-3*t1*t12
     6         5*t23*t33*t4*t54*t8*t9-7.786783946947111d-4*t1*t11*t126*t
     7         33*t54*t8*t9+1.268453166049416d-4*t1*t126*t23*t32*t54*t8*
     8         t9+2.3991838839497967d+1*t1*t11*t193*t4*t5*t8*t9-3.225612
     9         826113532d-1*t1*t109*t19*t23*t37*t4*t8*t9+1.7332206408224
     :         463d-1*t1*t109*t11*t19*t38*t8*t9+1.1774053387685064d-4*t1
     ;         *t31*t4*t52*t56*t8*t87-6.326559775702894d-5*t1*t125*t23*t
     <         31*t8*t87-3.402329015996045d-6*gammabb*t1*t29*t30*t4*t52*
     =         t56*t87+1.8281756662341322d-6*gammabb*t1*t125*t23*t29*t30
     >         *t87+7.868855442416534d-5*t1*t125*t23*t31*t4*t54*t8-4.228
     ?         17722016472d-5*t1*t11*t126*t31*t54*t8-2.2738503311369804d
     @         -6*gammabb*t1*t125*t23*t29*t30*t4*t54+1.2218094794768326d
     1         -6*gammabb*t1*t11*t126*t29*t30*t54)*wght+1.66723360671078
     2         12d-1*rhoa*rhob*t111*t39*wght+8.336168033553906d-2*rhoa*r
     3         hob*t127*t34*wght-2*t107*t41-t132*t36)+t14*t16*(-8.336168
     4         033553906d-2*rhoa*rhob*(-1.4287256594843722d-3*t1*t152*t3
     5         2*t4*t56*t8*t87*t9+5.117987261301516d-4*t1*t125*t32*t52*t
     6         8*t87*t9-6.365655797638701d-4*t1*t125*t32*t4*t52*t54*t8*t
     7         9+1.71023098797198d-4*t1*t126*t23*t32*t54*t8*t9+1.5874729
     8         549826355d-4*t1*t31*t4*t52*t56*t8*t87-8.529978768835859d-
     9         5*t1*t125*t23*t31*t8*t87-4.587294722559716d-6*gammabb*t1*
     :         t29*t30*t4*t52*t56*t87+2.4648940611561584d-6*gammabb*t1*t
     ;         125*t23*t29*t30*t87+1.0609426329397834d-4*t1*t125*t23*t31
     <         *t4*t54*t8-5.700769959906599d-5*t1*t11*t126*t31*t54*t8-3.
     =         0657886332788037d-6*gammabb*t1*t125*t23*t29*t30*t4*t54+1.
     >         6473422031868606d-6*gammabb*t1*t11*t126*t29*t30*t54)*wght
     ?         -8.336168033553906d-2*rhoa*t127*wght)+t18*t47*(3.33446721
     @         34215625d-1*rhoa*rhob*t17*t34*t40*t41*wght-2*t108*t41)+t2
     1         5*t49*(5.001700820132344d-1*rhoa*rhob*t17*t33*t40*t41*wgh
     2         t-3.0d+0*t108*t40)+t14*t47*(-1.6672336067107812d-1*rhoa*r
     3         hob*t34*t39*t41*wght-1.6672336067107812d-1*rhoa*rhob*t111
     4         *t17*t41*wght-1.6672336067107812d-1*rhoa*rhob*t132*t17*t3
     5         4*wght-2*t112*t41)+t18*t25*(-1.6672336067107812d-1*rhoa*r
     6         hob*t33*t39*t41*wght-1.6672336067107812d-1*rhoa*rhob*t109
     7         *t17*t41*wght-1.6672336067107812d-1*rhoa*rhob*t34*t39*t40
     8         *wght-1.6672336067107812d-1*rhoa*rhob*t111*t17*t40*wght-1
     9         .6672336067107812d-1*rhoa*rhob*t130*t17*t34*wght-1.667233
     :         6067107812d-1*rhoa*rhob*t132*t17*t33*wght-t110*t41-2.0d+0
     ;         *t112*t40)+t16*t49*(-5.001700820132344d-1*rhoa*rhob*t33*t
     <         39*t40*wght-5.001700820132344d-1*rhoa*rhob*t109*t17*t40*w
     =         ght-5.001700820132344d-1*rhoa*rhob*t130*t17*t33*wght-3.0d
     >         +0*t110*t40)+t16*t18*(3.3344672134215625d-1*rhoa*rhob*t10
     ?         9*t39*wght+1.6672336067107812d-1*rhoa*rhob*t127*t33*wght+
     @         1.6672336067107812d-1*rhoa*rhob*t17*t193*wght-4.0d+0*t107
     1         *t40-2.0d+0*t130*t36)+5.001700820132344d-1*rhoa*rhob*t128
     2         *t14*t149*t17*t34*wght+2.0006803280529378d+0*rhoa*rhob*t1
     3         29*t150*t16*t17*t33*wght+Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = t14*t16*(-8.336168033553906d-2*rho
     1         a*rhob*(7.084494764468265d-2*t1*t117*t19*t37*t4*t52*t8*t9
     2         -1.9033580930991081d-2*t1*t117*t19*t23*t38*t8*t9+2.558993
     3         630650758d-4*t1*t19*t29*t32*t38*t43*t52*t8+1.587472954982
     4         6355d-4*t1*t116*t19*t29*t37*t4*t52*t8-1.4287256594843722d
     5         -3*t1*t152*t19*t29*t32*t37*t4*t43*t8-4.2649893844179293d-
     6         5*t1*t116*t19*t23*t29*t38*t8-3.0657886332788037d-6*gammab
     7         b*t1*t113*t19*t23*t37*t4*t75+1.6473422031868606d-6*gammab
     8         b*t1*t11*t113*t19*t38*t75+2.558993630650758d-4*t1*t19*t29
     9         *t32*t37*t4*t42*t52-4.587294722559716d-6*gammabb*t1*t114*
     :         t115*t19*t37*t4*t52-6.875128571647359d-5*t1*t19*t23*t29*t
     ;         32*t38*t42+2.9854925690925505d-4*t1*t19*t23*t29*t30*t37*t
     <         4-1.6041966667177168d-4*t1*t11*t19*t29*t30*t38+1.23244703
     =         05780794d-6*gammabb*t1*t114*t115*t19*t23*t38)*wght-8.3361
     >         68033553906d-2*rhoa*t134*wght)+t14*t25*(8.336168033553906
     ?         d-2*rhoa*rhob*t17*(5.254465552928875d-2*t1*t117*t19*t37*t
     @         4*t52*t8*t9+2.3991838839497967d+1*t1*t11*t194*t4*t5*t8*t9
     1         -1.612806413056766d-1*t1*t121*t19*t23*t37*t4*t8*t9-1.4116
     2         926989963413d-2*t1*t117*t19*t23*t38*t8*t9+8.6661032041122
     3         31d-2*t1*t11*t121*t19*t38*t8*t9+2.168357784827018d-3*t1*t
     4         19*t29*t33*t37*t4*t43*t52*t8+1.8979679327108684d-4*t1*t19
     5         *t29*t32*t38*t43*t52*t8+1.1774053387685064d-4*t1*t116*t19
     6         *t29*t37*t4*t52*t8-1.612806413056766d-1*t1*t109*t23*t29*t
     7         4*t43*t5*t8-1.0596648048916558d-3*t1*t152*t19*t29*t32*t37
     8         *t4*t43*t8-5.825625504283547d-4*t1*t19*t23*t29*t33*t38*t4
     9         3*t8-3.163279887851447d-5*t1*t116*t19*t23*t29*t38*t8-2.27
     :         38503311369804d-6*gammabb*t1*t113*t19*t23*t37*t4*t75+1.22
     ;         18094794768326d-6*gammabb*t1*t11*t113*t19*t38*t75+1.89796
     <         79327108684d-4*t1*t19*t29*t32*t37*t4*t42*t52-3.4023290159
     =         96045d-6*gammabb*t1*t114*t115*t19*t37*t4*t52+8.6661032041
     >         12231d-2*t1*t109*t11*t29*t4*t42*t5-5.825625504283548d-4*t
     ?         1*t19*t23*t29*t33*t37*t4*t42+3.130287146672739d-4*t1*t11*
     @         t19*t29*t33*t38*t42-5.099181727518654d-5*t1*t19*t23*t29*t
     1         32*t38*t42+2.2142959214960126d-4*t1*t19*t23*t29*t30*t37*t
     2         4-1.1898090697543524d-4*t1*t11*t19*t29*t30*t38+9.14087833
     3         1170663d-7*gammabb*t1*t114*t115*t19*t23*t38)*wght+8.33616
     4         8033553906d-2*rhoa*rhob*t111*t44*wght+8.336168033553906d-
     5         2*rhoa*rhob*t123*t39*wght+8.336168033553906d-2*rhoa*rhob*
     6         t134*t34*wght-t107*t46-t119*t41-t137*t36)+t25*t49*(5.0017
     7         00820132344d-1*rhoa*rhob*t17*t33*t40*t46*wght-3.0d+0*t108
     8         *t45)+t18*t47*(3.3344672134215625d-1*rhoa*rhob*t17*t34*t4
     9         1*t45*wght-2*t108*t46)+t14*t47*(-1.6672336067107812d-1*rh
     :         oa*rhob*t34*t41*t44*wght-1.6672336067107812d-1*rhoa*rhob*
     ;         t123*t17*t41*wght-1.6672336067107812d-1*rhoa*rhob*t137*t1
     <         7*t34*wght-2*t112*t46)+t18*t25*(-1.6672336067107812d-1*rh
     =         oa*rhob*t33*t41*t44*wght-1.6672336067107812d-1*rhoa*rhob*
     >         t34*t40*t44*wght-1.6672336067107812d-1*rhoa*rhob*t121*t17
     ?         *t41*wght-1.6672336067107812d-1*rhoa*rhob*t123*t17*t40*wg
     @         ht-1.6672336067107812d-1*rhoa*rhob*t135*t17*t34*wght-1.66
     1         72336067107812d-1*rhoa*rhob*t137*t17*t33*wght-t110*t46-2.
     2         0d+0*t112*t45)+t16*t49*(-5.001700820132344d-1*rhoa*rhob*t
     3         33*t40*t44*wght-5.001700820132344d-1*rhoa*rhob*t121*t17*t
     4         40*wght-5.001700820132344d-1*rhoa*rhob*t135*t17*t33*wght-
     5         3.0d+0*t110*t45)+t16*t18*(1.6672336067107812d-1*rhoa*rhob
     6         *t109*t44*wght+1.6672336067107812d-1*rhoa*rhob*t121*t39*w
     7         ght+1.6672336067107812d-1*rhoa*rhob*t134*t33*wght+1.66723
     8         36067107812d-1*rhoa*rhob*t17*t194*wght-2.0d+0*t107*t45-2.
     9         0d+0*t119*t40-2.0d+0*t135*t36)+5.001700820132344d-1*rhoa*
     :         rhob*t14*t149*t17*t34*t41*t46*wght+2.0006803280529378d+0*
     ;         rhoa*rhob*t150*t16*t17*t33*t40*t45*wght+Cmat3(iq,D3_RB_GA
     <         A_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t14*t25*(8.336168033553906d-2*rhoa
     1         *rhob*t17*(-3.908227599821866d+0*t1*t198*t23*t4*t5*t8*t9+
     2         2.3991838839497967d+1*t1*t11*t199*t4*t5*t8*t9+2.168357784
     3         827018d-3*t1*t114*t33*t4*t5*t52*t77*t8-1.0596648048916558
     4         d-3*t1*t114*t152*t32*t4*t5*t77*t8-4.7213132654499207d-4*t
     5         1*t141*t32*t4*t5*t52*t75*t8+1.4491605732048626d-3*t1*t141
     6         *t23*t33*t4*t5*t75*t8+7.868855442416534d-5*t1*t115*t23*t4
     7         *t5*t75*t8+1.050893110585775d-1*t1*t117*t29*t4*t43*t5*t52
     8         *t8+1.1774053387685064d-4*t1*t114*t164*t4*t5*t52*t8-3.225
     9         612826113532d-1*t1*t121*t23*t29*t4*t43*t5*t8-7.7867839469
     :         47111d-4*t1*t11*t140*t33*t4*t5*t75+1.268453166049416d-4*t
     ;         1*t140*t23*t32*t4*t5*t75+6.34226583024708d-4*t1*t11*t113*
     <         t4*t5*t75+3.795935865421737d-4*t1*t114*t141*t32*t4*t5*t52
     =         -3.402329015996045d-6*gammabb*t1*t163*t197*t4*t5*t52-2.82
     >         33853979926826d-2*t1*t117*t23*t29*t4*t42*t5+1.73322064082
     ?         24463d-1*t1*t11*t121*t29*t4*t42*t5-1.1651251008567096d-3*
     @         t1*t114*t141*t23*t33*t4*t5-6.821550993410941d-6*gammabb*t
     1         1*t180*t196*t23*t4*t5+4.4285918429920257d-4*t1*t114*t115*
     2         t23*t4*t5-4.558990595062808d-6*gammabb*t1*t11*t114*t195*t
     3         4*t5)*wght+1.6672336067107812d-1*rhoa*rhob*t123*t44*wght+
     4         8.336168033553906d-2*rhoa*rhob*t142*t34*wght-2*t119*t46-t
     5         147*t36)+t14*t16*(-8.336168033553906d-2*rhoa*rhob*(-5.269
     6         388045346497d+0*t1*t198*t23*t4*t5*t8*t9-1.428725659484372
     7         2d-3*t1*t114*t152*t32*t4*t5*t77*t8-6.365655797638701d-4*t
     8         1*t141*t32*t4*t5*t52*t75*t8+1.0609426329397834d-4*t1*t115
     9         *t23*t4*t5*t75*t8+1.4168989528936532d-1*t1*t117*t29*t4*t4
     :         3*t5*t52*t8+1.5874729549826355d-4*t1*t114*t164*t4*t5*t52*
     ;         t8+1.71023098797198d-4*t1*t140*t23*t32*t4*t5*t75+8.551154
     <         939859897d-4*t1*t11*t113*t4*t5*t75+5.117987261301516d-4*t
     =         1*t114*t141*t32*t4*t5*t52-4.587294722559716d-6*gammabb*t1
     >         *t163*t197*t4*t5*t52-3.8067161861982163d-2*t1*t117*t23*t2
     ?         9*t4*t42*t5-9.197365899836413d-6*gammabb*t1*t180*t196*t23
     @         *t4*t5+5.970985138185101d-4*t1*t114*t115*t23*t4*t5-6.1467
     1         99265622613d-6*gammabb*t1*t11*t114*t195*t4*t5)*wght-8.336
     2         168033553906d-2*rhoa*t142*wght)+t18*t47*(3.33446721342156
     3         25d-1*rhoa*rhob*t17*t34*t45*t46*wght-2*t120*t46)+t25*t49*
     4         (5.001700820132344d-1*rhoa*rhob*t17*t33*t45*t46*wght-3.0d
     5         +0*t120*t45)+t14*t47*(-1.6672336067107812d-1*rhoa*rhob*t3
     6         4*t44*t46*wght-1.6672336067107812d-1*rhoa*rhob*t123*t17*t
     7         46*wght-1.6672336067107812d-1*rhoa*rhob*t147*t17*t34*wght
     8         -2*t124*t46)+t18*t25*(-1.6672336067107812d-1*rhoa*rhob*t3
     9         3*t44*t46*wght-1.6672336067107812d-1*rhoa*rhob*t121*t17*t
     :         46*wght-1.6672336067107812d-1*rhoa*rhob*t34*t44*t45*wght-
     ;         1.6672336067107812d-1*rhoa*rhob*t123*t17*t45*wght-1.66723
     <         36067107812d-1*rhoa*rhob*t145*t17*t34*wght-1.667233606710
     =         7812d-1*rhoa*rhob*t147*t17*t33*wght-t122*t46-2.0d+0*t124*
     >         t45)+t16*t49*(-5.001700820132344d-1*rhoa*rhob*t33*t44*t45
     ?         *wght-5.001700820132344d-1*rhoa*rhob*t121*t17*t45*wght-5.
     @         001700820132344d-1*rhoa*rhob*t145*t17*t33*wght-3.0d+0*t12
     1         2*t45)+t16*t18*(3.3344672134215625d-1*rhoa*rhob*t121*t44*
     2         wght+1.6672336067107812d-1*rhoa*rhob*t142*t33*wght+1.6672
     3         336067107812d-1*rhoa*rhob*t17*t199*wght-4.0d+0*t119*t45-2
     4         .0d+0*t145*t36)+5.001700820132344d-1*rhoa*rhob*t14*t143*t
     5         149*t17*t34*wght+2.0006803280529378d+0*rhoa*rhob*t144*t15
     6         0*t16*t17*t33*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t14*t25*(8.336168033553906d-2*rho
     1         a*rhob*t17*(2.168357784827018d-3*t1*t4*t40*t52*t56*t8*t87
     2         *t9-1.1651251008567096d-3*t1*t125*t23*t40*t8*t87*t9-3.181
     3         693560359581d-6*t1*t200*t23*t4*t8*t87*t9+1.70962147314855
     4         32d-6*t1*t11*t201*t8*t87*t9-9.521443141779977d-6*t1*t168*
     5         t4*t52*t55*t8*t9+3.827620142995552d-6*t1*t186*t52*t55*t8*
     6         t9+1.4491605732048626d-3*t1*t125*t23*t4*t40*t54*t8*t9-7.7
     7         86783946947111d-4*t1*t11*t126*t40*t54*t8*t9+2.39918388394
     8         97967d+1*t1*t11*t202*t4*t5*t8*t9-3.225612826113532d-1*t1*
     9         t130*t19*t23*t37*t4*t8*t9-7.123398104051166d-6*t1*t152*t1
     :         70*t186*t4*t8*t9+1.7332206408224463d-1*t1*t11*t130*t19*t3
     ;         8*t8*t9+5.116163245058206d-6*t1*t168*t200*t23*t8*t9)*wght
     <         +2.500850410066172d-1*rhoa*rhob*t127*t41*wght+2.500850410
     =         066172d-1*rhoa*rhob*t132*t39*wght)+t18*t25*(-3.3344672134
     >         215625d-1*rhoa*rhob*t39*t40*t41*wght-3.3344672134215625d-
     ?         1*rhoa*rhob*t130*t17*t41*wght-3.3344672134215625d-1*rhoa*
     @         rhob*t132*t17*t40*wght-t131*t41-2.0d+0*t133*t40)+t14*t47*
     1         (-3.3344672134215625d-1*rhoa*rhob*t132*t17*t41*wght-1.667
     2         2336067107812d-1*rhoa*rhob*t128*t39*wght-2*t133*t41)+t16*
     3         t49*(-1.0003401640264689d+0*rhoa*rhob*t130*t17*t40*wght-5
     4         .001700820132344d-1*rhoa*rhob*t129*t39*wght-3.0d+0*t131*t
     5         40)+t16*t18*(5.001700820132344d-1*rhoa*rhob*t127*t40*wght
     6         +5.001700820132344d-1*rhoa*rhob*t130*t39*wght+1.667233606
     7         7107812d-1*rhoa*rhob*t17*t202*wght)-8.336168033553906d-2*
     8         rhoa*rhob*t14*t16*(-4.2898161846251925d-6*t1*t200*t23*t4*
     9         t8*t87*t9+2.3050497246084803d-6*t1*t11*t201*t8*t87*t9-1.2
     :         837578514626072d-5*t1*t168*t4*t52*t55*t8*t9+5.16070656287
     ;         9681d-6*t1*t186*t52*t55*t8*t9-9.604340549010496d-6*t1*t15
     <         2*t170*t186*t4*t8*t9+6.898024424877311d-6*t1*t168*t200*t2
     =         3*t8*t9)*wght+1.5005102460397032d+0*rhoa*rhob*t129*t17*t2
     >         5*t41*t49*wght+1.0003401640264689d+0*rhoa*rhob*t128*t17*t
     ?         18*t40*t47*wght+5.001700820132344d-1*rhoa*rhob*t14*t149*t
     @         17*t41**3*wght+2.0006803280529378d+0*rhoa*rhob*t150*t16*t
     1         17*t40**3*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t14*t25*(8.336168033553906d-2*rho
     1         a*rhob*t17*(2.3991838839497967d+1*t1*t11*t203*t4*t5*t8*t9
     2         -1.612806413056766d-1*t1*t135*t19*t23*t37*t4*t8*t9+8.6661
     3         03204112231d-2*t1*t11*t135*t19*t38*t8*t9-7.12339810405116
     4         6d-6*t1*t152*t29*t4*t43*t56*t8*t87+2.551746761997034d-6*t
     5         1*t125*t29*t43*t52*t8*t87+1.275873380998517d-6*t1*t29*t4*
     6         t42*t52*t56*t87-6.855658748377997d-7*t1*t125*t23*t29*t42*
     7         t87-3.173814380593326d-6*t1*t125*t29*t4*t43*t52*t54*t8+8.
     8         526938741763678d-7*t1*t126*t23*t29*t43*t54*t8+2.168357784
     9         827018d-3*t1*t19*t29*t37*t4*t40*t43*t52*t8-1.612806413056
     :         766d-1*t1*t130*t23*t29*t4*t43*t5*t8-5.825625504283547d-4*
     ;         t1*t19*t23*t29*t38*t40*t43*t8+8.526938741763678d-7*t1*t12
     <         5*t23*t29*t4*t42*t54-4.581785548038123d-7*t1*t11*t126*t29
     =         *t42*t54+8.666103204112231d-2*t1*t11*t130*t29*t4*t42*t5-5
     >         .825625504283548d-4*t1*t19*t23*t29*t37*t4*t40*t42+3.13028
     ?         7146672739d-4*t1*t11*t19*t29*t38*t40*t42)*wght+8.33616803
     @         3553906d-2*rhoa*rhob*t127*t46*wght+8.336168033553906d-2*r
     1         hoa*rhob*t132*t44*wght+1.6672336067107812d-1*rhoa*rhob*t1
     2         34*t41*wght+1.6672336067107812d-1*rhoa*rhob*t137*t39*wght
     3         )+t18*t47*(6.668934426843125d-1*rhoa*rhob*t17*t40*t41*t46
     4         *wght+3.3344672134215625d-1*rhoa*rhob*t128*t17*t45*wght)+
     5         t25*t49*(5.001700820132344d-1*rhoa*rhob*t129*t17*t46*wght
     6         +1.0003401640264689d+0*rhoa*rhob*t17*t40*t41*t45*wght)+t1
     7         6*t18*(1.6672336067107812d-1*rhoa*rhob*t127*t45*wght+1.66
     8         72336067107812d-1*rhoa*rhob*t130*t44*wght+3.3344672134215
     9         625d-1*rhoa*rhob*t134*t40*wght+3.3344672134215625d-1*rhoa
     :         *rhob*t135*t39*wght+1.6672336067107812d-1*rhoa*rhob*t17*t
     ;         203*wght)+t18*t25*(-3.3344672134215625d-1*rhoa*rhob*t40*t
     <         41*t44*wght-3.3344672134215625d-1*rhoa*rhob*t135*t17*t41*
     =         wght-3.3344672134215625d-1*rhoa*rhob*t137*t17*t40*wght-t1
     >         31*t46-2.0d+0*t133*t45)+t16*t49*(-5.001700820132344d-1*rh
     ?         oa*rhob*t129*t44*wght-1.0003401640264689d+0*rhoa*rhob*t13
     @         5*t17*t40*wght-3.0d+0*t131*t45)+t14*t47*(-1.6672336067107
     1         812d-1*rhoa*rhob*t128*t44*wght-3.3344672134215625d-1*rhoa
     2         *rhob*t137*t17*t41*wght-2*t133*t46)-8.336168033553906d-2*
     3         rhoa*rhob*t14*t16*(-9.604340549010496d-6*t1*t152*t29*t4*t
     4         43*t56*t8*t87+3.4404710419197876d-6*t1*t125*t29*t43*t52*t
     5         8*t87+1.7202355209598938d-6*t1*t29*t4*t42*t52*t56*t87-9.2
     6         43352729335597d-7*t1*t125*t23*t29*t42*t87-4.2791928382086
     7         9d-6*t1*t125*t29*t4*t43*t52*t54*t8+1.1496707374795517d-6*
     8         t1*t126*t23*t29*t43*t54*t8+1.1496707374795517d-6*t1*t125*
     9         t23*t29*t4*t42*t54-6.177533261950727d-7*t1*t11*t126*t29*t
     :         42*t54)*wght+5.001700820132344d-1*rhoa*rhob*t128*t14*t149
     ;         *t17*t46*wght+2.0006803280529378d+0*rhoa*rhob*t129*t150*t
     <         16*t17*t45*wght+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = t14*t25*(8.336168033553906d-2*rho
     1         a*rhob*t17*(2.3991838839497967d+1*t1*t11*t204*t4*t5*t8*t9
     2         +2.168357784827018d-3*t1*t114*t4*t40*t5*t52*t77*t8+1.2758
     3         73380998517d-6*t1*t114*t19*t38*t52*t77*t8-7.1233981040511
     4         66d-6*t1*t114*t152*t19*t37*t4*t77*t8-3.173814380593326d-6
     5         *t1*t141*t19*t37*t4*t52*t75*t8+1.4491605732048626d-3*t1*t
     6         141*t23*t4*t40*t5*t75*t8+8.526938741763679d-7*t1*t141*t19
     7         *t23*t38*t75*t8-3.225612826113532d-1*t1*t135*t23*t29*t4*t
     8         43*t5*t8-7.786783946947111d-4*t1*t11*t140*t4*t40*t5*t75+8
     9         .526938741763679d-7*t1*t140*t19*t23*t37*t4*t75-4.58178554
     :         8038123d-7*t1*t11*t140*t19*t38*t75+2.551746761997034d-6*t
     ;         1*t114*t141*t19*t37*t4*t52+1.7332206408224463d-1*t1*t11*t
     <         135*t29*t4*t42*t5-1.1651251008567096d-3*t1*t114*t141*t23*
     =         t4*t40*t5-6.855658748377997d-7*t1*t114*t141*t19*t23*t38)*
     >         wght+1.6672336067107812d-1*rhoa*rhob*t134*t46*wght+1.6672
     ?         336067107812d-1*rhoa*rhob*t137*t44*wght+8.336168033553906
     @         d-2*rhoa*rhob*t142*t41*wght+8.336168033553906d-2*rhoa*rho
     1         b*t147*t39*wght)+t18*t47*(3.3344672134215625d-1*rhoa*rhob
     2         *t17*t41*t45*t46*wght-2*t139*t46)+t25*t49*(5.001700820132
     3         344d-1*rhoa*rhob*t17*t40*t45*t46*wght-3.0d+0*t139*t45)+t1
     4         4*t47*(-1.6672336067107812d-1*rhoa*rhob*t41*t44*t46*wght-
     5         1.6672336067107812d-1*rhoa*rhob*t137*t17*t46*wght-1.66723
     6         36067107812d-1*rhoa*rhob*t147*t17*t41*wght-2*t138*t46)+t1
     7         8*t25*(-1.6672336067107812d-1*rhoa*rhob*t40*t44*t46*wght-
     8         1.6672336067107812d-1*rhoa*rhob*t135*t17*t46*wght-1.66723
     9         36067107812d-1*rhoa*rhob*t41*t44*t45*wght-1.6672336067107
     :         812d-1*rhoa*rhob*t137*t17*t45*wght-1.6672336067107812d-1*
     ;         rhoa*rhob*t145*t17*t41*wght-1.6672336067107812d-1*rhoa*rh
     <         ob*t147*t17*t40*wght-t136*t46-2.0d+0*t138*t45)+t16*t49*(-
     =         5.001700820132344d-1*rhoa*rhob*t40*t44*t45*wght-5.0017008
     >         20132344d-1*rhoa*rhob*t135*t17*t45*wght-5.001700820132344
     ?         d-1*rhoa*rhob*t145*t17*t40*wght-3.0d+0*t136*t45)+t16*t18*
     @         (3.3344672134215625d-1*rhoa*rhob*t134*t45*wght+3.33446721
     1         34215625d-1*rhoa*rhob*t135*t44*wght+1.6672336067107812d-1
     2         *rhoa*rhob*t142*t40*wght+1.6672336067107812d-1*rhoa*rhob*
     3         t145*t39*wght+1.6672336067107812d-1*rhoa*rhob*t17*t204*wg
     4         ht)-8.336168033553906d-2*rhoa*rhob*t14*t16*(1.72023552095
     5         98938d-6*t1*t114*t19*t38*t52*t77*t8-9.604340549010496d-6*
     6         t1*t114*t152*t19*t37*t4*t77*t8-4.2791928382086913d-6*t1*t
     7         141*t19*t37*t4*t52*t75*t8+1.1496707374795517d-6*t1*t141*t
     8         19*t23*t38*t75*t8+1.1496707374795517d-6*t1*t140*t19*t23*t
     9         37*t4*t75-6.177533261950726d-7*t1*t11*t140*t19*t38*t75+3.
     :         4404710419197876d-6*t1*t114*t141*t19*t37*t4*t52-9.2433527
     ;         29335595d-7*t1*t114*t141*t19*t23*t38)*wght+5.001700820132
     <         344d-1*rhoa*rhob*t14*t143*t149*t17*t41*wght+2.00068032805
     =         29378d+0*rhoa*rhob*t144*t150*t16*t17*t40*wght+Cmat3(iq,D3
     >         _GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t14*t25*(8.336168033553906d-2*rho
     1         a*rhob*t17*(2.3991838839497967d+1*t1*t11*t207*t4*t5*t8*t9
     2         +2.168357784827018d-3*t1*t114*t4*t45*t5*t52*t77*t8-9.5214
     3         43141779977d-6*t1*t180*t4*t5*t52*t76*t8+1.449160573204862
     4         6d-3*t1*t141*t23*t4*t45*t5*t75*t8-3.225612826113532d-1*t1
     5         *t145*t23*t29*t4*t43*t5*t8-3.181693560359581d-6*t1*t114*t
     6         206*t23*t4*t5*t8-7.123398104051166d-6*t1*t152*t182*t197*t
     7         4*t5*t8+3.827620142995552d-6*t1*t197*t4*t5*t52*t76-7.7867
     8         83946947111d-4*t1*t11*t140*t4*t45*t5*t75-1.16512510085670
     9         96d-3*t1*t114*t141*t23*t4*t45*t5+1.7332206408224463d-1*t1
     :         *t11*t145*t29*t4*t42*t5+5.116163245058206d-6*t1*t180*t206
     ;         *t23*t4*t5+1.7096214731485532d-6*t1*t11*t114*t205*t4*t5)*
     <         wght+2.500850410066172d-1*rhoa*rhob*t142*t46*wght+2.50085
     =         0410066172d-1*rhoa*rhob*t147*t44*wght)+t18*t25*(-3.334467
     >         2134215625d-1*rhoa*rhob*t44*t45*t46*wght-3.33446721342156
     ?         25d-1*rhoa*rhob*t145*t17*t46*wght-3.3344672134215625d-1*r
     @         hoa*rhob*t147*t17*t45*wght-t146*t46-2.0d+0*t148*t45)+t14*
     1         t47*(-3.3344672134215625d-1*rhoa*rhob*t147*t17*t46*wght-1
     2         .6672336067107812d-1*rhoa*rhob*t143*t44*wght-2*t148*t46)+
     3         t16*t49*(-1.0003401640264689d+0*rhoa*rhob*t145*t17*t45*wg
     4         ht-5.001700820132344d-1*rhoa*rhob*t144*t44*wght-3.0d+0*t1
     5         46*t45)+t16*t18*(5.001700820132344d-1*rhoa*rhob*t142*t45*
     6         wght+5.001700820132344d-1*rhoa*rhob*t145*t44*wght+1.66723
     7         36067107812d-1*rhoa*rhob*t17*t207*wght)-8.336168033553906
     8         d-2*rhoa*rhob*t14*t16*(-1.2837578514626072d-5*t1*t180*t4*
     9         t5*t52*t76*t8-4.2898161846251925d-6*t1*t114*t206*t23*t4*t
     :         5*t8-9.604340549010496d-6*t1*t152*t182*t197*t4*t5*t8+5.16
     ;         0706562879681d-6*t1*t197*t4*t5*t52*t76+6.898024424877311d
     <         -6*t1*t180*t206*t23*t4*t5+2.3050497246084803d-6*t1*t11*t1
     =         14*t205*t4*t5)*wght+1.5005102460397032d+0*rhoa*rhob*t144*
     >         t17*t25*t46*t49*wght+1.0003401640264689d+0*rhoa*rhob*t143
     ?         *t17*t18*t45*t47*wght+5.001700820132344d-1*rhoa*rhob*t14*
     @         t149*t17*t46**3*wght+2.0006803280529378d+0*rhoa*rhob*t150
     1         *t16*t17*t45**3*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
