C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_ssbd.F 
C> The nwxcm_x_ssbd functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_ssbd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{0.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{0.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{0.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{0.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{0.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{0.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{0.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{0.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131, 094103 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
      subroutine nwxcm_x_ssbd(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = abs(rhoa)
            t6 = 1/t5**5.333333333333333d+0
            t7 = 1.5904488902586256d-3*gammaaa**2.0d+0*t6+1.0d+0
            t8 = 1/t7
            t9 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3*
     1         gammaaa*t2*t8
            t10 = gammaaa**2
            t11 = 1/t3**2
            t12 = 1/rhoa**3.6666666666666664d+0
            t13 = 1/t7**2
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t9*wght
            Amat(iq,D1_RA) = (-1.2407009817987999d+0*rhoa**3.33333333333
     1         3333d-1*t9-9.305257363491d-1*t1*(1.5163101135555887d-2*ga
     2         mmaaa*t12*t8-4.823227474784833d-5*gammaaa**3.0d+0*t13/(rh
     3         oa**1.6666666666666669d+0*t5**7.333333333333333d+0)-8.664
     4         9263035052d-3*gammaaa*t12*t4+3.8886799859849375d-5*t10*t1
     5         1/rhoa**6.333333333333333d+0))*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*(-5.6
     1         86162925833458d-3*t2*t8+1.8087103030443127d-5*t10*t13*t2*
     2         t6+3.24934736381445d-3*t2*t4-1.4582549947443516d-5*gammaa
     3         a*t11/rhoa**5.333333333333333d+0)*wght
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
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = abs(rhoa)
            t6 = 1/t5**5.333333333333333d+0
            t7 = 1.5904488902586256d-3*gammaaa**2.0d+0*t6+1.0d+0
            t8 = 1/t7
            t9 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3*
     1         gammaaa*t2*t8
            t10 = rhob**1.3333333333333333d+0
            t11 = 1/rhob**2.6666666666666666d+0
            t12 = 4.4878396535373416d-3*gammabb*t11+1.0d+0
            t13 = 1/t12
            t14 = abs(rhob)
            t15 = 1/t14**5.333333333333333d+0
            t16 = 1.5904488902586256d-3*gammabb**2.0d+0*t15+1.0d+0
            t17 = 1/t16
            t18 = 3.24934736381445d-3*gammabb*t11*t13-5.686162925833458d
     1         -3*gammabb*t11*t17
            t19 = gammaaa**2
            t20 = 1/t3**2
            t21 = 1/rhoa**3.6666666666666664d+0
            t22 = 1/t7**2
            t23 = gammabb**2
            t24 = 1/t12**2
            t25 = 1/rhob**3.6666666666666664d+0
            t26 = 1/t16**2
            fnc(iq) = (-9.305257363491d-1*t1*t9-9.305257363491d-1*t10*t1
     1         8)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987999d+0*rhoa**3.33333333333
     1         3333d-1*t9-9.305257363491d-1*t1*(1.5163101135555887d-2*ga
     2         mmaaa*t21*t8-4.823227474784833d-5*gammaaa**3.0d+0*t22/(rh
     3         oa**1.6666666666666669d+0*t5**7.333333333333333d+0)-8.664
     4         9263035052d-3*gammaaa*t21*t4+3.8886799859849375d-5*t19*t2
     5         0/rhoa**6.333333333333333d+0))*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t10*(-4.823227474784833
     1         d-5*gammabb**3.0d+0*t26/(rhob**1.6666666666666669d+0*t14*
     2         *7.333333333333333d+0)+1.5163101135555887d-2*gammabb*t17*
     3         t25-8.6649263035052d-3*gammabb*t13*t25+3.8886799859849375
     4         d-5*t23*t24/rhob**6.333333333333333d+0)-1.240700981798799
     5         9d+0*rhob**3.333333333333333d-1*t18)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*(-5.6
     1         86162925833458d-3*t2*t8+1.8087103030443127d-5*t19*t2*t22*
     2         t6+3.24934736381445d-3*t2*t4-1.4582549947443516d-5*gammaa
     3         a*t20/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t10*(1.8
     1         087103030443127d-5*t11*t15*t23*t26-1.4582549947443516d-5*
     2         gammabb*t24/rhob**5.333333333333333d+0-5.686162925833458d
     3         -3*t11*t17+3.24934736381445d-3*t11*t13)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = abs(rhoa)
            t6 = 1/t5**5.333333333333333d+0
            t7 = 1.5904488902586256d-3*gammaaa**2.0d+0*t6+1.0d+0
            t8 = 1/t7
            t9 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3*
     1         gammaaa*t2*t8
            t10 = gammaaa**2
            t11 = 1/t3**2
            t12 = 1/rhoa**3.6666666666666664d+0
            t13 = 1/t7**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t9*wght
            Amat(iq,D1_RA) = -1.2407009817987999d+0*rhoa**3.333333333333
     1         333d-1*t9*wght-9.305257363491d-1*t1*(1.5163101135555887d-
     2         2*gammaaa*t12*t8-4.823227474784833d-5*gammaaa**3.0d+0*t13
     3         /(rhoa**1.6666666666666669d+0*t5**7.333333333333333d+0)-8
     4         .6649263035052d-3*gammaaa*t12*t4+3.8886799859849375d-5*t1
     5         0*t11/rhoa**6.333333333333333d+0)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*(-5.6
     1         86162925833458d-3*t2*t8+1.8087103030443127d-5*t10*t13*t2*
     2         t6+3.24934736381445d-3*t2*t4-1.4582549947443516d-5*gammaa
     3         a*t11/rhoa**5.333333333333333d+0)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammabb*t2+1.0d+0
            t4 = 1/t3
            t5 = abs(rhob)
            t6 = 1/t5**5.333333333333333d+0
            t7 = 1.5904488902586256d-3*gammabb**2.0d+0*t6+1.0d+0
            t8 = 1/t7
            t9 = 3.24934736381445d-3*gammabb*t2*t4-5.686162925833458d-3*
     1         gammabb*t2*t8
            t10 = gammabb**2
            t11 = 1/t3**2
            t12 = 1/rhob**3.6666666666666664d+0
            t13 = 1/t7**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t9*wght
            Amat(iq,D1_RB) = -1.2407009817987999d+0*rhob**3.333333333333
     1         333d-1*t9*wght-9.305257363491d-1*t1*(1.5163101135555887d-
     2         2*gammabb*t12*t8-4.823227474784833d-5*gammabb**3.0d+0*t13
     3         /(rhob**1.6666666666666669d+0*t5**7.333333333333333d+0)-8
     4         .6649263035052d-3*gammabb*t12*t4+3.8886799859849375d-5*t1
     5         0*t11/rhob**6.333333333333333d+0)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t1*(-5.6
     1         86162925833458d-3*t2*t8+1.8087103030443127d-5*t10*t13*t2*
     2         t6+3.24934736381445d-3*t2*t4-1.4582549947443516d-5*gammab
     3         b*t11/rhob**5.333333333333333d+0)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_ssbd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{0.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{0.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{0.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{0.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{0.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{0.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{0.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{0.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131, 094103 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
      subroutine nwxcm_x_ssbd_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**2
            t13 = 1/t3**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = gammaaa**3.0d+0
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammaaa*t15*t9-8.6649263035052d-
     1         3*gammaaa*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhoa**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammaaa*t13*t21
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = gammaaa**3
            t25 = 1/t3**3
            t26 = 1/rhoa**4.666666666666667d+0
            t27 = 1/t8**3
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t10*wght
            Amat(iq,D1_RA) = (-9.305257363491d-1*t1*t20-1.24070098179879
     1         99d+0*t10*t11)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t22*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t1*(-5.559803749703
     1         826d-2*gammaaa*t26*t9+3.537033481508877d-4*t16*t18*t23/t6
     2         **9.333333333333332d+0-8.182503237052212d-7*gammaaa**5.0d
     3         +0*t23*t27/t6**1.4666666666666667d+1+3.1771396446185735d-
     4         2*gammaaa*t26*t4+9.307611861877458d-7*t24*t25/rhoa**10+2.
     5         0900652390734276d-4*t16*t18*t19*t2-3.499811987386443d-4*t
     6         12*t13/rhoa**7.333333333333333d+0)-4.135669939329333d-1*t
     7         10*t23-2.4814019635975998d+0*t11*t20)*wght+Amat2(iq,D2_RA
     8         _RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*t1*(1.516310113555
     1         5887d-2*t15*t9-4.823227474784833d-5*t12*t15*t18*t7+3.0684
     2         3871389458d-7*gammaaa**4.0d+0*t17*t27/t6**1.2666666666666
     3         665d+1-1.44696824243545d-4*t17*t18*t19*t5-8.6649263035052
     4         d-3*t15*t4-3.490354448204047d-7*t12*t25/rhoa**9+1.1666039
     5         957954812d-4*gammaaa*t13*t14)-1.2407009817987999d+0*t11*t
     6         22)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*(5.426130909132938d-5*gammaaa*t18*t2*t7-1.1506645177
     2         104675d-7*t2*t24*t27/t6**1.0666666666666666d+1+1.30888291
     3         8076518d-7*gammaaa*t25/rhoa**8-2.916509989488703d-5*t13*t
     4         21)*wght
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
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhob**1.3333333333333333d+0
            t12 = 1/rhob**2.6666666666666666d+0
            t13 = 4.4878396535373416d-3*gammabb*t12+1.0d+0
            t14 = 1/t13
            t15 = gammabb**2.0d+0
            t16 = abs(rhob)
            t17 = 1/t16**5.333333333333333d+0
            t18 = 1.5904488902586256d-3*t15*t17+1.0d+0
            t19 = 1/t18
            t20 = 3.24934736381445d-3*gammabb*t12*t14-5.686162925833458d
     1         -3*gammabb*t12*t19
            t21 = rhoa**3.333333333333333d-1
            t22 = gammaaa**2
            t23 = 1/t3**2
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = gammaaa**3.0d+0
            t27 = 1/rhoa**1.6666666666666669d+0
            t28 = 1/t8**2
            t29 = 1/t6**7.333333333333333d+0
            t30 = 1.5163101135555887d-2*gammaaa*t25*t9-8.6649263035052d-
     1         3*gammaaa*t25*t4-4.823227474784833d-5*t26*t27*t28*t29+3.8
     2         886799859849375d-5*t22*t23*t24
            t31 = rhob**3.333333333333333d-1
            t32 = gammabb**2
            t33 = 1/t13**2
            t34 = 1/rhob**6.333333333333333d+0
            t35 = 1/rhob**3.6666666666666664d+0
            t36 = gammabb**3.0d+0
            t37 = 1/rhob**1.6666666666666669d+0
            t38 = 1/t18**2
            t39 = 1/t16**7.333333333333333d+0
            t40 = -4.823227474784833d-5*t36*t37*t38*t39+1.51631011355558
     1         87d-2*gammabb*t19*t35-8.6649263035052d-3*gammabb*t14*t35+
     2         3.8886799859849375d-5*t32*t33*t34
            t41 = 1/rhoa**5.333333333333333d+0
            t42 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t2*t
     1         22*t28*t7-1.4582549947443516d-5*gammaaa*t23*t41+3.2493473
     2         6381445d-3*t2*t4
            t43 = 1/rhob**5.333333333333333d+0
            t44 = -1.4582549947443516d-5*gammabb*t33*t43+1.8087103030443
     1         127d-5*t12*t17*t32*t38-5.686162925833458d-3*t12*t19+3.249
     2         34736381445d-3*t12*t14
            t45 = 1/rhoa**6.666666666666666d-1
            t46 = gammaaa**3
            t47 = 1/t3**3
            t48 = 1/rhoa**4.666666666666667d+0
            t49 = 1/t8**3
            t50 = 1/rhob**6.666666666666666d-1
            t51 = gammabb**3
            t52 = 1/t13**3
            t53 = 1/rhob**4.666666666666667d+0
            t54 = 1/t18**3
            fnc(iq) = (-9.305257363491d-1*t11*t20-9.305257363491d-1*t1*t
     1         10)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491d-1*t1*t30-1.24070098179879
     1         99d+0*t10*t21)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t11*t40-1.2407009817987
     1         999d+0*t20*t31)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t42*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t11*t44*
     1         wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t1*(-5.559803749703
     1         826d-2*gammaaa*t48*t9+3.537033481508877d-4*t26*t28*t45/t6
     2         **9.333333333333332d+0-8.182503237052212d-7*gammaaa**5.0d
     3         +0*t45*t49/t6**1.4666666666666667d+1+3.1771396446185735d-
     4         2*gammaaa*t4*t48+9.307611861877458d-7*t46*t47/rhoa**10+2.
     5         0900652390734276d-4*t2*t26*t28*t29-3.499811987386443d-4*t
     6         22*t23/rhoa**7.333333333333333d+0)-4.135669939329333d-1*t
     7         10*t45-2.4814019635975998d+0*t21*t30)*wght+Amat2(iq,D2_RA
     8         _RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491d-1*t11*(-8.18250323705
     1         2212d-7*gammabb**5.0d+0*t50*t54/t16**1.4666666666666667d+
     2         1-5.559803749703826d-2*gammabb*t19*t53+3.1771396446185735
     3         d-2*gammabb*t14*t53+9.307611861877458d-7*t51*t52/rhob**10
     4         +3.537033481508877d-4*t36*t38*t50/t16**9.333333333333332d
     5         +0+2.0900652390734276d-4*t12*t36*t38*t39-3.49981198738644
     6         3d-4*t32*t33/rhob**7.333333333333333d+0)-4.13566993932933
     7         3d-1*t20*t50-2.4814019635975998d+0*t31*t40)*wght+Amat2(iq
     8         ,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*t1*(1.516310113555
     1         5887d-2*t25*t9-4.823227474784833d-5*t22*t25*t28*t7+3.0684
     2         3871389458d-7*gammaaa**4.0d+0*t27*t49/t6**1.2666666666666
     3         665d+1-1.44696824243545d-4*t27*t28*t29*t5-3.4903544482040
     4         47d-7*t22*t47/rhoa**9-8.6649263035052d-3*t25*t4+1.1666039
     5         957954812d-4*gammaaa*t23*t24)-1.2407009817987999d+0*t21*t
     6         42)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363491d-1*t11*(3.06843871389
     1         458d-7*gammabb**4.0d+0*t37*t54/t16**1.2666666666666665d+1
     2         -3.490354448204047d-7*t32*t52/rhob**9-1.44696824243545d-4
     3         *t15*t37*t38*t39-4.823227474784833d-5*t17*t32*t35*t38+1.5
     4         163101135555887d-2*t19*t35-8.6649263035052d-3*t14*t35+1.1
     5         666039957954812d-4*gammabb*t33*t34)-1.2407009817987999d+0
     6         *t31*t44)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*(5.426130909132938d-5*gammaaa*t2*t28*t7-1.1506645177
     2         104675d-7*t2*t46*t49/t6**1.0666666666666666d+1+1.30888291
     3         8076518d-7*gammaaa*t47/rhoa**8-2.916509989488703d-5*t23*t
     4         41)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t11*(-1.1506645177104675d-7*t12*t51*t54/t16**1.06666666
     2         66666666d+1+1.308882918076518d-7*gammabb*t52/rhob**8-2.91
     3         6509989488703d-5*t33*t43+5.426130909132938d-5*gammabb*t12
     4         *t17*t38)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**2
            t13 = 1/t3**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = gammaaa**3.0d+0
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammaaa*t15*t9-8.6649263035052d-
     1         3*gammaaa*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhoa**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammaaa*t13*t21
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = gammaaa**3
            t25 = 1/t3**3
            t26 = 1/rhoa**4.666666666666667d+0
            t27 = 1/t8**3
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t10*wght
            Amat(iq,D1_RA) = -9.305257363491d-1*t1*t20*wght-1.2407009817
     1         987999d+0*t10*t11*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t22*w
     1         ght
            Amat2(iq,D2_RA_RA) = -9.305257363491d-1*t1*(-5.5598037497038
     1         26d-2*gammaaa*t26*t9+3.537033481508877d-4*t16*t18*t23/t6*
     2         *9.333333333333332d+0-8.182503237052212d-7*gammaaa**5.0d+
     3         0*t23*t27/t6**1.4666666666666667d+1+3.1771396446185735d-2
     4         *gammaaa*t26*t4+9.307611861877458d-7*t24*t25/rhoa**10+2.0
     5         900652390734276d-4*t16*t18*t19*t2-3.499811987386443d-4*t1
     6         2*t13/rhoa**7.333333333333333d+0)*wght-4.135669939329333d
     7         -1*t10*t23*wght-2.4814019635975998d+0*t11*t20*wght+Amat2(
     8         iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363491d-1*t1*(1.5163101135555
     1         887d-2*t15*t9-4.823227474784833d-5*t12*t15*t18*t7+3.06843
     2         871389458d-7*gammaaa**4.0d+0*t17*t27/t6**1.26666666666666
     3         65d+1-1.44696824243545d-4*t17*t18*t19*t5-8.6649263035052d
     4         -3*t15*t4-3.490354448204047d-7*t12*t25/rhoa**9+1.16660399
     5         57954812d-4*gammaaa*t13*t14)*wght-1.2407009817987999d+0*t
     6         11*t22*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*(5.426130909132938d-5*gammaaa*t18*t2*t7-1.1506645177
     2         104675d-7*t2*t24*t27/t6**1.0666666666666666d+1+1.30888291
     3         8076518d-7*gammaaa*t25/rhoa**8-2.916509989488703d-5*t13*t
     4         21)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammabb*t2+1.0d+0
            t4 = 1/t3
            t5 = gammabb**2.0d+0
            t6 = abs(rhob)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammabb*t2*t4-5.686162925833458d-3
     1         *gammabb*t2*t9
            t11 = rhob**3.333333333333333d-1
            t12 = gammabb**2
            t13 = 1/t3**2
            t14 = 1/rhob**6.333333333333333d+0
            t15 = 1/rhob**3.6666666666666664d+0
            t16 = gammabb**3.0d+0
            t17 = 1/rhob**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammabb*t15*t9-8.6649263035052d-
     1         3*gammabb*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhob**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammabb*t13*t21
            t23 = 1/rhob**6.666666666666666d-1
            t24 = gammabb**3
            t25 = 1/t3**3
            t26 = 1/rhob**4.666666666666667d+0
            t27 = 1/t8**3
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t10*wght
            Amat(iq,D1_RB) = -9.305257363491d-1*t1*t20*wght-1.2407009817
     1         987999d+0*t10*t11*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t1*t22*w
     1         ght
            Amat2(iq,D2_RB_RB) = -9.305257363491d-1*t1*(-5.5598037497038
     1         26d-2*gammabb*t26*t9+3.537033481508877d-4*t16*t18*t23/t6*
     2         *9.333333333333332d+0-8.182503237052212d-7*gammabb**5.0d+
     3         0*t23*t27/t6**1.4666666666666667d+1+3.1771396446185735d-2
     4         *gammabb*t26*t4+9.307611861877458d-7*t24*t25/rhob**10+2.0
     5         900652390734276d-4*t16*t18*t19*t2-3.499811987386443d-4*t1
     6         2*t13/rhob**7.333333333333333d+0)*wght-4.135669939329333d
     7         -1*t10*t23*wght-2.4814019635975998d+0*t11*t20*wght+Amat2(
     8         iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363491d-1*t1*(1.5163101135555
     1         887d-2*t15*t9-4.823227474784833d-5*t12*t15*t18*t7+3.06843
     2         871389458d-7*gammabb**4.0d+0*t17*t27/t6**1.26666666666666
     3         65d+1-1.44696824243545d-4*t17*t18*t19*t5-8.6649263035052d
     4         -3*t15*t4-3.490354448204047d-7*t12*t25/rhob**9+1.16660399
     5         57954812d-4*gammabb*t13*t14)*wght-1.2407009817987999d+0*t
     6         11*t22*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t1*(5.426130909132938d-5*gammabb*t18*t2*t7-1.1506645177
     2         104675d-7*t2*t24*t27/t6**1.0666666666666666d+1+1.30888291
     3         8076518d-7*gammabb*t25/rhob**8-2.916509989488703d-5*t13*t
     4         21)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_ssbd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{0.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{0.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{0.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{0.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{0.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{0.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{0.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{0.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131, 094103 (2009)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
      subroutine nwxcm_x_ssbd_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**2
            t13 = 1/t3**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = gammaaa**3.0d+0
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammaaa*t15*t9-8.6649263035052d-
     1         3*gammaaa*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhoa**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammaaa*t13*t21
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = gammaaa**3
            t25 = 1/t3**3
            t26 = 1/rhoa**10
            t27 = 1/rhoa**7.333333333333333d+0
            t28 = 1/rhoa**4.666666666666667d+0
            t29 = gammaaa**5.0d+0
            t30 = 1/t8**3
            t31 = 1/t6**1.4666666666666667d+1
            t32 = 1/t6**9.333333333333332d+0
            t33 = -5.559803749703826d-2*gammaaa*t28*t9+3.177139644618573
     1         5d-2*gammaaa*t28*t4+3.537033481508877d-4*t16*t18*t23*t32-
     2         8.182503237052212d-7*t23*t29*t30*t31-3.499811987386443d-4
     3         *t12*t13*t27+9.307611861877458d-7*t24*t25*t26+2.090065239
     4         0734276d-4*t16*t18*t19*t2
            t34 = 1/rhoa**9
            t35 = gammaaa**4.0d+0
            t36 = 1/t6**1.2666666666666665d+1
            t37 = 1.5163101135555887d-2*t15*t9-4.823227474784833d-5*t12*
     1         t15*t18*t7-1.44696824243545d-4*t17*t18*t19*t5-8.664926303
     2         5052d-3*t15*t4+3.06843871389458d-7*t17*t30*t35*t36-3.4903
     3         54448204047d-7*t12*t25*t34+1.1666039957954812d-4*gammaaa*
     4         t13*t14
            t38 = 1/rhoa**8
            t39 = 1/t6**1.0666666666666666d+1
            t40 = 5.426130909132938d-5*gammaaa*t18*t2*t7-1.1506645177104
     1         675d-7*t2*t24*t30*t39+1.308882918076518d-7*gammaaa*t25*t3
     2         8-2.916509989488703d-5*t13*t21
            t41 = gammaaa**4
            t42 = 1/t3**4
            t43 = 1/rhoa**5.666666666666667d+0
            t44 = 1/t8**4
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t10*wght
            Amat(iq,D1_RA) = (-9.305257363491d-1*t1*t20-1.24070098179879
     1         99d+0*t10*t11)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t22*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t1*t33-4.1356699393
     1         29333d-1*t10*t23-2.4814019635975998d+0*t11*t20)*wght+Amat
     2         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*t1*t37-1.240700981
     1         7987999d+0*t11*t22)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*t40*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491d-1*t1*(2.5945750831
     1         951186d-1*gammaaa*t43*t9-3.3012312494082846d-3*t11*t16*t1
     2         8/t6**1.1333333333333333d+1+1.8001507121514865d-5*t11*t29
     3         *t30/t6**1.6666666666666663d+1-2.0822165108651677d-8*gamm
     4         aaa**7.0d+0*t11*t44/t6**2.2d+1-1.4826651674886676d-1*gamm
     5         aaa*t4*t43+3.341685567477454d-8*t41*t42/rhoa**1.366666666
     6         6666666d+1-1.7685167407544386d-3*t16*t17*t18*t32+4.091251
     7         618526106d-6*t17*t29*t30*t31-1.768446253756717d-5*t24*t25
     8         /rhoa**11-1.0289551946207642d-3*t15*t16*t18*t19+2.9467552
     9         78268585d-3*t12*t13/rhoa**8.333333333333334d+0)-3.7221029
     :         453964d+0*t11*t33-1.2407009817987999d+0*t20*t23+2.7571132
     ;         92886222d-1*t10*t17)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491d-1*t1*(-5.55980374
     1         9703826d-2*t28*t9+1.768516740754439d-4*t12*t18*t28*t7-2.2
     2         501883901893585d-6*t23*t30*t35/t6**1.4666666666666664d+1+
     3         7.808311915744381d-9*gammaaa**6.0d+0*t23*t44/t6**2.0d+1+1
     4         .061110044452663d-3*t18*t23*t32*t5+6.270195717220283d-4*t
     5         18*t19*t2*t5-1.2531320878040456d-8*t24*t42/rhoa**1.266666
     6         6666666665d+1+3.1771396446185735d-2*t28*t4-1.329656776020
     7         9846d-6*t2*t30*t35*t36-4.091251618526106d-6*t23*t30*t31*t
     8         35-8.425473302967363d-4*gammaaa*t13*t27+5.93360256194688d
     9         -6*t12*t25*t26)-2.4814019635975998d+0*t11*t37-4.135669939
     :         329333d-1*t22*t23)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363491d-1*t1*(-1.4469682
     1         4243545d-4*gammaaa*t15*t18*t7-2.9281169684041436d-9*t17*t
     2         29*t44/t6**1.8d+1+4.699245329265171d-9*t12*t42/rhoa**1.16
     3         66666666666665d+1+3.06843871389458d-7*t15*t24*t30*t39+2.1
     4         479070997262062d-6*t16*t17*t30*t36-1.7451772241020236d-6*
     5         gammaaa*t25*t34-2.8939364848709d-4*gammaaa*t17*t18*t19+1.
     6         555471994393975d-4*t13*t14)-1.2407009817987999d+0*t11*t40
     7         )*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*t1*(5.426130909132938d-5*t18*t2*t7+1.0980438631
     2         515539d-9*t2*t41*t44/t6**1.6d+1-1.7622169984744393d-9*gam
     3         maaa*t42/rhoa**1.0666666666666666d+1-6.903987106262807d-7
     4         *t12*t2*t30*t39+3.9266487542295536d-7*t25*t38)*wght
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
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhob**1.3333333333333333d+0
            t12 = 1/rhob**2.6666666666666666d+0
            t13 = 4.4878396535373416d-3*gammabb*t12+1.0d+0
            t14 = 1/t13
            t15 = gammabb**2.0d+0
            t16 = abs(rhob)
            t17 = 1/t16**5.333333333333333d+0
            t18 = 1.5904488902586256d-3*t15*t17+1.0d+0
            t19 = 1/t18
            t20 = 3.24934736381445d-3*gammabb*t12*t14-5.686162925833458d
     1         -3*gammabb*t12*t19
            t21 = rhoa**3.333333333333333d-1
            t22 = gammaaa**2
            t23 = 1/t3**2
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = gammaaa**3.0d+0
            t27 = 1/rhoa**1.6666666666666669d+0
            t28 = 1/t8**2
            t29 = 1/t6**7.333333333333333d+0
            t30 = 1.5163101135555887d-2*gammaaa*t25*t9-8.6649263035052d-
     1         3*gammaaa*t25*t4-4.823227474784833d-5*t26*t27*t28*t29+3.8
     2         886799859849375d-5*t22*t23*t24
            t31 = rhob**3.333333333333333d-1
            t32 = gammabb**2
            t33 = 1/t13**2
            t34 = 1/rhob**6.333333333333333d+0
            t35 = 1/rhob**3.6666666666666664d+0
            t36 = gammabb**3.0d+0
            t37 = 1/rhob**1.6666666666666669d+0
            t38 = 1/t18**2
            t39 = 1/t16**7.333333333333333d+0
            t40 = -4.823227474784833d-5*t36*t37*t38*t39+1.51631011355558
     1         87d-2*gammabb*t19*t35-8.6649263035052d-3*gammabb*t14*t35+
     2         3.8886799859849375d-5*t32*t33*t34
            t41 = 1/rhoa**5.333333333333333d+0
            t42 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t2*t
     1         22*t28*t7-1.4582549947443516d-5*gammaaa*t23*t41+3.2493473
     2         6381445d-3*t2*t4
            t43 = 1/rhob**5.333333333333333d+0
            t44 = -1.4582549947443516d-5*gammabb*t33*t43+1.8087103030443
     1         127d-5*t12*t17*t32*t38-5.686162925833458d-3*t12*t19+3.249
     2         34736381445d-3*t12*t14
            t45 = 1/rhoa**6.666666666666666d-1
            t46 = gammaaa**3
            t47 = 1/t3**3
            t48 = 1/rhoa**10
            t49 = 1/rhoa**7.333333333333333d+0
            t50 = 1/rhoa**4.666666666666667d+0
            t51 = gammaaa**5.0d+0
            t52 = 1/t8**3
            t53 = 1/t6**1.4666666666666667d+1
            t54 = 1/t6**9.333333333333332d+0
            t55 = -5.559803749703826d-2*gammaaa*t50*t9+3.537033481508877
     1         d-4*t26*t28*t45*t54-8.182503237052212d-7*t45*t51*t52*t53+
     2         3.1771396446185735d-2*gammaaa*t4*t50-3.499811987386443d-4
     3         *t22*t23*t49+9.307611861877458d-7*t46*t47*t48+2.090065239
     4         0734276d-4*t2*t26*t28*t29
            t56 = 1/rhob**6.666666666666666d-1
            t57 = gammabb**3
            t58 = 1/t13**3
            t59 = 1/rhob**10
            t60 = 1/rhob**7.333333333333333d+0
            t61 = 1/rhob**4.666666666666667d+0
            t62 = gammabb**5.0d+0
            t63 = 1/t18**3
            t64 = 1/t16**1.4666666666666667d+1
            t65 = 1/t16**9.333333333333332d+0
            t66 = 3.537033481508877d-4*t36*t38*t56*t65-8.182503237052212
     1         d-7*t56*t62*t63*t64-5.559803749703826d-2*gammabb*t19*t61+
     2         3.1771396446185735d-2*gammabb*t14*t61-3.499811987386443d-
     3         4*t32*t33*t60+9.307611861877458d-7*t57*t58*t59+2.09006523
     4         90734276d-4*t12*t36*t38*t39
            t67 = 1/rhoa**9
            t68 = gammaaa**4.0d+0
            t69 = 1/t6**1.2666666666666665d+1
            t70 = 1.5163101135555887d-2*t25*t9-4.823227474784833d-5*t22*
     1         t25*t28*t7+3.06843871389458d-7*t27*t52*t68*t69-3.49035444
     2         8204047d-7*t22*t47*t67-1.44696824243545d-4*t27*t28*t29*t5
     3         -8.6649263035052d-3*t25*t4+1.1666039957954812d-4*gammaaa*
     4         t23*t24
            t71 = 1/rhob**9
            t72 = gammabb**4.0d+0
            t73 = 1/t16**1.2666666666666665d+1
            t74 = 3.06843871389458d-7*t37*t63*t72*t73-3.490354448204047d
     1         -7*t32*t58*t71-1.44696824243545d-4*t15*t37*t38*t39-4.8232
     2         27474784833d-5*t17*t32*t35*t38+1.5163101135555887d-2*t19*
     3         t35-8.6649263035052d-3*t14*t35+1.1666039957954812d-4*gamm
     4         abb*t33*t34
            t75 = 1/rhoa**8
            t76 = 1/t6**1.0666666666666666d+1
            t77 = -1.1506645177104675d-7*t2*t46*t52*t76+1.30888291807651
     1         8d-7*gammaaa*t47*t75+5.426130909132938d-5*gammaaa*t2*t28*
     2         t7-2.916509989488703d-5*t23*t41
            t78 = 1/rhob**8
            t79 = 1/t16**1.0666666666666666d+1
            t80 = -1.1506645177104675d-7*t12*t57*t63*t79+1.3088829180765
     1         18d-7*gammabb*t58*t78-2.916509989488703d-5*t33*t43+5.4261
     2         30909132938d-5*gammabb*t12*t17*t38
            t81 = gammaaa**4
            t82 = 1/t3**4
            t83 = 1/rhoa**5.666666666666667d+0
            t84 = 1/t8**4
            t85 = gammabb**4
            t86 = 1/t13**4
            t87 = 1/rhob**5.666666666666667d+0
            t88 = 1/t18**4
            fnc(iq) = (-9.305257363491d-1*t11*t20-9.305257363491d-1*t1*t
     1         10)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491d-1*t1*t30-1.24070098179879
     1         99d+0*t10*t21)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t11*t40-1.2407009817987
     1         999d+0*t20*t31)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t42*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t11*t44*
     1         wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t1*t55-4.1356699393
     1         29333d-1*t10*t45-2.4814019635975998d+0*t21*t30)*wght+Amat
     2         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491d-1*t11*t66-4.135669939
     1         329333d-1*t20*t56-2.4814019635975998d+0*t31*t40)*wght+Ama
     2         t2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*t1*t70-1.240700981
     1         7987999d+0*t21*t42)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363491d-1*t11*t74-1.24070098
     1         17987999d+0*t31*t44)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*t77*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t11*t80*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491d-1*t1*(2.5945750831
     1         951186d-1*gammaaa*t83*t9-2.0822165108651677d-8*gammaaa**7
     2         .0d+0*t21*t84/t6**2.2d+1-1.4826651674886676d-1*gammaaa*t4
     3         *t83+3.341685567477454d-8*t81*t82/rhoa**1.366666666666666
     4         6d+1-3.3012312494082846d-3*t21*t26*t28/t6**1.133333333333
     5         3333d+1+1.8001507121514865d-5*t21*t51*t52/t6**1.666666666
     6         6666663d+1-1.7685167407544386d-3*t26*t27*t28*t54+4.091251
     7         618526106d-6*t27*t51*t52*t53-1.768446253756717d-5*t46*t47
     8         /rhoa**11-1.0289551946207642d-3*t25*t26*t28*t29+2.9467552
     9         78268585d-3*t22*t23/rhoa**8.333333333333334d+0)-3.7221029
     :         453964d+0*t21*t55-1.2407009817987999d+0*t30*t45+2.7571132
     ;         92886222d-1*t10*t27)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363491d-1*t11*(-2.08221651
     1         08651677d-8*gammabb**7.0d+0*t31*t88/t16**2.2d+1+2.5945750
     2         831951186d-1*gammabb*t19*t87-1.4826651674886676d-1*gammab
     3         b*t14*t87+3.341685567477454d-8*t85*t86/rhob**1.3666666666
     4         666666d+1-1.7685167407544386d-3*t36*t37*t38*t65+4.0912516
     5         18526106d-6*t37*t62*t63*t64+1.8001507121514865d-5*t31*t62
     6         *t63/t16**1.6666666666666663d+1-1.768446253756717d-5*t57*
     7         t58/rhob**11-1.0289551946207642d-3*t35*t36*t38*t39-3.3012
     8         312494082846d-3*t31*t36*t38/t16**1.1333333333333333d+1+2.
     9         946755278268585d-3*t32*t33/rhob**8.333333333333334d+0)-3.
     :         7221029453964d+0*t31*t66-1.2407009817987999d+0*t40*t56+2.
     ;         757113292886222d-1*t20*t37)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491d-1*t1*(-5.55980374
     1         9703826d-2*t50*t9+7.808311915744381d-9*gammaaa**6.0d+0*t4
     2         5*t84/t6**2.0d+1-1.2531320878040456d-8*t46*t82/rhoa**1.26
     3         66666666666665d+1+1.768516740754439d-4*t22*t28*t50*t7-1.3
     4         296567760209846d-6*t2*t52*t68*t69-2.2501883901893585d-6*t
     5         45*t52*t68/t6**1.4666666666666664d+1-4.091251618526106d-6
     6         *t45*t52*t53*t68+1.061110044452663d-3*t28*t45*t5*t54+3.17
     7         71396446185735d-2*t4*t50+6.270195717220283d-4*t2*t28*t29*
     8         t5-8.425473302967363d-4*gammaaa*t23*t49+5.93360256194688d
     9         -6*t22*t47*t48)-2.4814019635975998d+0*t21*t70-4.135669939
     :         329333d-1*t42*t45)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-9.305257363491d-1*t11*(7.80831191
     1         5744381d-9*gammabb**6.0d+0*t56*t88/t16**2.0d+1-1.25313208
     2         78040456d-8*t57*t86/rhob**1.2666666666666665d+1-1.3296567
     3         760209846d-6*t12*t63*t72*t73-4.091251618526106d-6*t56*t63
     4         *t64*t72-2.2501883901893585d-6*t56*t63*t72/t16**1.4666666
     5         666666664d+1+1.061110044452663d-3*t15*t38*t56*t65+1.76851
     6         6740754439d-4*t17*t32*t38*t61-5.559803749703826d-2*t19*t6
     7         1+3.1771396446185735d-2*t14*t61-8.425473302967363d-4*gamm
     8         abb*t33*t60+5.93360256194688d-6*t32*t58*t59+6.27019571722
     9         0283d-4*t12*t15*t38*t39)-2.4814019635975998d+0*t31*t74-4.
     :         135669939329333d-1*t44*t56)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363491d-1*t1*(-2.9281169
     1         684041436d-9*t27*t51*t84/t6**1.8d+1+4.699245329265171d-9*
     2         t22*t82/rhoa**1.1666666666666665d+1+3.06843871389458d-7*t
     3         25*t46*t52*t76-1.44696824243545d-4*gammaaa*t25*t28*t7+2.1
     4         479070997262062d-6*t26*t27*t52*t69-1.7451772241020236d-6*
     5         gammaaa*t47*t67-2.8939364848709d-4*gammaaa*t27*t28*t29+1.
     6         555471994393975d-4*t23*t24)-1.2407009817987999d+0*t21*t77
     7         )*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-9.305257363491d-1*t11*(-2.928116
     1         9684041436d-9*t37*t62*t88/t16**1.8d+1+4.699245329265171d-
     2         9*t32*t86/rhob**1.1666666666666665d+1+3.06843871389458d-7
     3         *t35*t57*t63*t79+2.1479070997262062d-6*t36*t37*t63*t73-1.
     4         7451772241020236d-6*gammabb*t58*t71-2.8939364848709d-4*ga
     5         mmabb*t37*t38*t39-1.44696824243545d-4*gammabb*t17*t35*t38
     6         +1.555471994393975d-4*t33*t34)-1.2407009817987999d+0*t31*
     7         t80)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*t1*(1.0980438631515539d-9*t2*t81*t84/t6**1.6d+1
     2         -1.7622169984744393d-9*gammaaa*t82/rhoa**1.06666666666666
     3         66d+1-6.903987106262807d-7*t2*t22*t52*t76+3.9266487542295
     4         536d-7*t47*t75+5.426130909132938d-5*t2*t28*t7)*wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491d-1*t11*(1.0980438631515539d-9*t12*t85*t88/t16**1.6
     2         d+1-1.7622169984744393d-9*gammabb*t86/rhob**1.06666666666
     3         66666d+1-6.903987106262807d-7*t12*t32*t63*t79+3.926648754
     4         2295536d-7*t58*t78+5.426130909132938d-5*t12*t17*t38)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammaaa*t2+1.0d+0
            t4 = 1/t3
            t5 = gammaaa**2.0d+0
            t6 = abs(rhoa)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammaaa*t2*t4-5.686162925833458d-3
     1         *gammaaa*t2*t9
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**2
            t13 = 1/t3**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = gammaaa**3.0d+0
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammaaa*t15*t9-8.6649263035052d-
     1         3*gammaaa*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhoa**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammaaa*t13*t21
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = gammaaa**3
            t25 = 1/t3**3
            t26 = 1/rhoa**10
            t27 = 1/rhoa**7.333333333333333d+0
            t28 = 1/rhoa**4.666666666666667d+0
            t29 = gammaaa**5.0d+0
            t30 = 1/t8**3
            t31 = 1/t6**1.4666666666666667d+1
            t32 = 1/t6**9.333333333333332d+0
            t33 = -5.559803749703826d-2*gammaaa*t28*t9+3.177139644618573
     1         5d-2*gammaaa*t28*t4+3.537033481508877d-4*t16*t18*t23*t32-
     2         8.182503237052212d-7*t23*t29*t30*t31-3.499811987386443d-4
     3         *t12*t13*t27+9.307611861877458d-7*t24*t25*t26+2.090065239
     4         0734276d-4*t16*t18*t19*t2
            t34 = 1/rhoa**9
            t35 = gammaaa**4.0d+0
            t36 = 1/t6**1.2666666666666665d+1
            t37 = 1.5163101135555887d-2*t15*t9-4.823227474784833d-5*t12*
     1         t15*t18*t7-1.44696824243545d-4*t17*t18*t19*t5-8.664926303
     2         5052d-3*t15*t4+3.06843871389458d-7*t17*t30*t35*t36-3.4903
     3         54448204047d-7*t12*t25*t34+1.1666039957954812d-4*gammaaa*
     4         t13*t14
            t38 = 1/rhoa**8
            t39 = 1/t6**1.0666666666666666d+1
            t40 = 5.426130909132938d-5*gammaaa*t18*t2*t7-1.1506645177104
     1         675d-7*t2*t24*t30*t39+1.308882918076518d-7*gammaaa*t25*t3
     2         8-2.916509989488703d-5*t13*t21
            t41 = gammaaa**4
            t42 = 1/t3**4
            t43 = 1/rhoa**5.666666666666667d+0
            t44 = 1/t8**4
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t10*wght
            Amat(iq,D1_RA) = -9.305257363491d-1*t1*t20*wght-1.2407009817
     1         987999d+0*t10*t11*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t1*t22*w
     1         ght
            Amat2(iq,D2_RA_RA) = -9.305257363491d-1*t1*t33*wght-4.135669
     1         939329333d-1*t10*t23*wght-2.4814019635975998d+0*t11*t20*w
     2         ght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363491d-1*t1*t37*wght-1.24070
     1         09817987999d+0*t11*t22*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t1*t40*wght
            Amat3(iq,D3_RA_RA_RA) = -9.305257363491d-1*t1*(2.59457508319
     1         51186d-1*gammaaa*t43*t9-3.3012312494082846d-3*t11*t16*t18
     2         /t6**1.1333333333333333d+1+1.8001507121514865d-5*t11*t29*
     3         t30/t6**1.6666666666666663d+1-2.0822165108651677d-8*gamma
     4         aa**7.0d+0*t11*t44/t6**2.2d+1-1.4826651674886676d-1*gamma
     5         aa*t4*t43+3.341685567477454d-8*t41*t42/rhoa**1.3666666666
     6         666666d+1-1.7685167407544386d-3*t16*t17*t18*t32+4.0912516
     7         18526106d-6*t17*t29*t30*t31-1.768446253756717d-5*t24*t25/
     8         rhoa**11-1.0289551946207642d-3*t15*t16*t18*t19+2.94675527
     9         8268585d-3*t12*t13/rhoa**8.333333333333334d+0)*wght-3.722
     :         1029453964d+0*t11*t33*wght-1.2407009817987999d+0*t20*t23*
     ;         wght+2.757113292886222d-1*t10*t17*wght+Amat3(iq,D3_RA_RA_
     <         RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363491d-1*t1*(-5.559803749
     1         703826d-2*t28*t9+1.768516740754439d-4*t12*t18*t28*t7-2.25
     2         01883901893585d-6*t23*t30*t35/t6**1.4666666666666664d+1+7
     3         .808311915744381d-9*gammaaa**6.0d+0*t23*t44/t6**2.0d+1+1.
     4         061110044452663d-3*t18*t23*t32*t5+6.270195717220283d-4*t1
     5         8*t19*t2*t5-1.2531320878040456d-8*t24*t42/rhoa**1.2666666
     6         666666665d+1+3.1771396446185735d-2*t28*t4-1.3296567760209
     7         846d-6*t2*t30*t35*t36-4.091251618526106d-6*t23*t30*t31*t3
     8         5-8.425473302967363d-4*gammaaa*t13*t27+5.93360256194688d-
     9         6*t12*t25*t26)*wght-2.4814019635975998d+0*t11*t37*wght-4.
     :         135669939329333d-1*t22*t23*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -9.305257363491d-1*t1*(-1.44696824
     1         243545d-4*gammaaa*t15*t18*t7-2.9281169684041436d-9*t17*t2
     2         9*t44/t6**1.8d+1+4.699245329265171d-9*t12*t42/rhoa**1.166
     3         6666666666665d+1+3.06843871389458d-7*t15*t24*t30*t39+2.14
     4         79070997262062d-6*t16*t17*t30*t36-1.7451772241020236d-6*g
     5         ammaaa*t25*t34-2.8939364848709d-4*gammaaa*t17*t18*t19+1.5
     6         55471994393975d-4*t13*t14)*wght-1.2407009817987999d+0*t11
     7         *t40*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*t1*(5.426130909132938d-5*t18*t2*t7+1.0980438631
     2         515539d-9*t2*t41*t44/t6**1.6d+1-1.7622169984744393d-9*gam
     3         maaa*t42/rhoa**1.0666666666666666d+1-6.903987106262807d-7
     4         *t12*t2*t30*t39+3.9266487542295536d-7*t25*t38)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = 4.4878396535373416d-3*gammabb*t2+1.0d+0
            t4 = 1/t3
            t5 = gammabb**2.0d+0
            t6 = abs(rhob)
            t7 = 1/t6**5.333333333333333d+0
            t8 = 1.5904488902586256d-3*t5*t7+1.0d+0
            t9 = 1/t8
            t10 = 3.24934736381445d-3*gammabb*t2*t4-5.686162925833458d-3
     1         *gammabb*t2*t9
            t11 = rhob**3.333333333333333d-1
            t12 = gammabb**2
            t13 = 1/t3**2
            t14 = 1/rhob**6.333333333333333d+0
            t15 = 1/rhob**3.6666666666666664d+0
            t16 = gammabb**3.0d+0
            t17 = 1/rhob**1.6666666666666669d+0
            t18 = 1/t8**2
            t19 = 1/t6**7.333333333333333d+0
            t20 = 1.5163101135555887d-2*gammabb*t15*t9-8.6649263035052d-
     1         3*gammabb*t15*t4-4.823227474784833d-5*t16*t17*t18*t19+3.8
     2         886799859849375d-5*t12*t13*t14
            t21 = 1/rhob**5.333333333333333d+0
            t22 = -5.686162925833458d-3*t2*t9+1.8087103030443127d-5*t12*
     1         t18*t2*t7+3.24934736381445d-3*t2*t4-1.4582549947443516d-5
     2         *gammabb*t13*t21
            t23 = 1/rhob**6.666666666666666d-1
            t24 = gammabb**3
            t25 = 1/t3**3
            t26 = 1/rhob**10
            t27 = 1/rhob**7.333333333333333d+0
            t28 = 1/rhob**4.666666666666667d+0
            t29 = gammabb**5.0d+0
            t30 = 1/t8**3
            t31 = 1/t6**1.4666666666666667d+1
            t32 = 1/t6**9.333333333333332d+0
            t33 = -5.559803749703826d-2*gammabb*t28*t9+3.177139644618573
     1         5d-2*gammabb*t28*t4+3.537033481508877d-4*t16*t18*t23*t32-
     2         8.182503237052212d-7*t23*t29*t30*t31-3.499811987386443d-4
     3         *t12*t13*t27+9.307611861877458d-7*t24*t25*t26+2.090065239
     4         0734276d-4*t16*t18*t19*t2
            t34 = 1/rhob**9
            t35 = gammabb**4.0d+0
            t36 = 1/t6**1.2666666666666665d+1
            t37 = 1.5163101135555887d-2*t15*t9-4.823227474784833d-5*t12*
     1         t15*t18*t7-1.44696824243545d-4*t17*t18*t19*t5-8.664926303
     2         5052d-3*t15*t4+3.06843871389458d-7*t17*t30*t35*t36-3.4903
     3         54448204047d-7*t12*t25*t34+1.1666039957954812d-4*gammabb*
     4         t13*t14
            t38 = 1/rhob**8
            t39 = 1/t6**1.0666666666666666d+1
            t40 = 5.426130909132938d-5*gammabb*t18*t2*t7-1.1506645177104
     1         675d-7*t2*t24*t30*t39+1.308882918076518d-7*gammabb*t25*t3
     2         8-2.916509989488703d-5*t13*t21
            t41 = gammabb**4
            t42 = 1/t3**4
            t43 = 1/rhob**5.666666666666667d+0
            t44 = 1/t8**4
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t10*wght
            Amat(iq,D1_RB) = -9.305257363491d-1*t1*t20*wght-1.2407009817
     1         987999d+0*t10*t11*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t1*t22*w
     1         ght
            Amat2(iq,D2_RB_RB) = -9.305257363491d-1*t1*t33*wght-4.135669
     1         939329333d-1*t10*t23*wght-2.4814019635975998d+0*t11*t20*w
     2         ght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363491d-1*t1*t37*wght-1.24070
     1         09817987999d+0*t11*t22*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t1*t40*wght
            Amat3(iq,D3_RB_RB_RB) = -9.305257363491d-1*t1*(2.59457508319
     1         51186d-1*gammabb*t43*t9-3.3012312494082846d-3*t11*t16*t18
     2         /t6**1.1333333333333333d+1+1.8001507121514865d-5*t11*t29*
     3         t30/t6**1.6666666666666663d+1-2.0822165108651677d-8*gamma
     4         bb**7.0d+0*t11*t44/t6**2.2d+1-1.4826651674886676d-1*gamma
     5         bb*t4*t43+3.341685567477454d-8*t41*t42/rhob**1.3666666666
     6         666666d+1-1.7685167407544386d-3*t16*t17*t18*t32+4.0912516
     7         18526106d-6*t17*t29*t30*t31-1.768446253756717d-5*t24*t25/
     8         rhob**11-1.0289551946207642d-3*t15*t16*t18*t19+2.94675527
     9         8268585d-3*t12*t13/rhob**8.333333333333334d+0)*wght-3.722
     :         1029453964d+0*t11*t33*wght-1.2407009817987999d+0*t20*t23*
     ;         wght+2.757113292886222d-1*t10*t17*wght+Amat3(iq,D3_RB_RB_
     <         RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363491d-1*t1*(-5.559803749
     1         703826d-2*t28*t9+1.768516740754439d-4*t12*t18*t28*t7-2.25
     2         01883901893585d-6*t23*t30*t35/t6**1.4666666666666664d+1+7
     3         .808311915744381d-9*gammabb**6.0d+0*t23*t44/t6**2.0d+1+1.
     4         061110044452663d-3*t18*t23*t32*t5+6.270195717220283d-4*t1
     5         8*t19*t2*t5-1.2531320878040456d-8*t24*t42/rhob**1.2666666
     6         666666665d+1+3.1771396446185735d-2*t28*t4-1.3296567760209
     7         846d-6*t2*t30*t35*t36-4.091251618526106d-6*t23*t30*t31*t3
     8         5-8.425473302967363d-4*gammabb*t13*t27+5.93360256194688d-
     9         6*t12*t25*t26)*wght-2.4814019635975998d+0*t11*t37*wght-4.
     :         135669939329333d-1*t22*t23*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -9.305257363491d-1*t1*(-1.44696824
     1         243545d-4*gammabb*t15*t18*t7-2.9281169684041436d-9*t17*t2
     2         9*t44/t6**1.8d+1+4.699245329265171d-9*t12*t42/rhob**1.166
     3         6666666666665d+1+3.06843871389458d-7*t15*t24*t30*t39+2.14
     4         79070997262062d-6*t16*t17*t30*t36-1.7451772241020236d-6*g
     5         ammabb*t25*t34-2.8939364848709d-4*gammabb*t17*t18*t19+1.5
     6         55471994393975d-4*t13*t14)*wght-1.2407009817987999d+0*t11
     7         *t40*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491d-1*t1*(5.426130909132938d-5*t18*t2*t7+1.0980438631
     2         515539d-9*t2*t41*t44/t6**1.6d+1-1.7622169984744393d-9*gam
     3         mabb*t42/rhob**1.0666666666666666d+1-6.903987106262807d-7
     4         *t12*t2*t30*t39+3.9266487542295536d-7*t25*t38)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
