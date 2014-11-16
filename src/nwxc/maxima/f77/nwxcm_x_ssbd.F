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
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
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
            t14 = -1.8610514726982d+0*t1*(-5.686162925833458d-3*t2*t8+1.
     1         8087103030443127d-5*t10*t13*t2*t6+3.24934736381445d-3*t2*
     2         t4-1.4582549947443516d-5*gammaaa*t11/rhoa**5.333333333333
     3         333d+0)*wght
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t9*wght
            Amat(iq,D1_RA) = -2.4814019635975998d+0*rhoa**3.333333333333
     1         333d-1*t9*wght-1.8610514726982d+0*t1*(1.5163101135555887d
     2         -2*gammaaa*t12*t8-4.823227474784833d-5*gammaaa**3.0d+0*t1
     3         3/(rhoa**1.6666666666666669d+0*t5**7.333333333333333d+0)-
     4         8.6649263035052d-3*gammaaa*t12*t4+3.8886799859849375d-5*t
     5         10*t11/rhoa**6.333333333333333d+0)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t14+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t14+Cmat(iq,D1_GAB)
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
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
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
C> \brief Evaluate the nwxcm_x_ssbd functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_2}\,\sigma_{\beta\beta}}\over{.004487839653537342\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}-{{.005686162925833458\,{\it t_2}\,\sigma_{\beta\beta}}\over{{{.001590448890258625\,\sigma_{\beta\beta}^{2.0}}\over{\left| \rho_\beta\right| ^{5.333333333333333}}}+1.0}}\right)-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{.004487839653537342\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}-{{.005686162925833458\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{{{.001590448890258625\,\sigma_{\alpha\alpha}^{2.0}}\over{\left| \rho_\alpha\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,\left({{0.00324934736381445\,{\it t_3}\,\sigma_{ss}}\over{.004487839653537342\,{\it t_3}\,\sigma_{ss}+1.0}}-{{.005686162925833458\,{\it t_3}\,\sigma_{ss}}\over{{{.001590448890258625\,\sigma_{ss}^{2.0}}\over{\left| \rho_s\right| ^{5.333333333333333}}}+1.0}}\right)\\\\ 
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
