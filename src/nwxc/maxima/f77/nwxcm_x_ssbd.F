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
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_2}\,
C>    \sigma_{\beta\beta}}\over{0.004487839653537343\,{\it t_2}
C>    \,\sigma_{\beta\beta}+1.0}}-{{0.00568616292583346\,{\it t_2}
C>    \,\sigma_{\beta\beta}}\over{0.001590448890258626\,\left({
C>    \it t_2}\,\sigma_{\beta\beta}\right)^{2.0}+1.0}}\right)
C>    -0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.004487839653537343\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}-{{0.00568616292583346\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}}\over{0.001590448890258626\,\left({
C>    \it t_1}\,\sigma_{\alpha\alpha}\right)^{2.0}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.004487839653537343\,{\it t_3}\,\sigma_{ss}+1.0}}
C>    -{{0.00568616292583346\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.001590448890258626\,\left({\it t_3}\,
C>    \sigma_{ss}\right)^{2.0}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131
C>    , 094103 (2009)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ssbd(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.590448890258626d-3*(gammaaa*t1)**2.0d+0+1.0d+0
            t3 = 1/t2
            t4 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t5 = 1/t4
            t6 = 3.249347363814451d-3*gammaaa*t1*t5-5.68616292583346d-3*
     1         gammaaa*t1*t3
            t7 = rhoa**1.3333333333333333d+0
            t8 = 1/t2**2
            t9 = gammaaa**2
            t10 = 1/t4**2
            t11 = 1/rhoa**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t6*t7*wght
            Amat(iq,D1_RA) = (-9.305257363491d-1*t7*(3.8886799859849397d
     1         -5*t10*t9/rhoa**6.333333333333333d+0-4.8232274747848364d-
     2         5*gammaaa**3*t8/rhoa**9-8.664926303505202d-3*gammaaa*t11*
     3         t5+1.5163101135555893d-2*gammaaa*t11*t3)-1.24070098179879
     4         99d+0*rhoa**3.333333333333333d-1*t6)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t7*(1.80
     1         87103030443136d-5*t8*t9/rhoa**8+3.249347363814451d-3*t1*t
     2         5-5.68616292583346d-3*t1*t3-1.4582549947443524d-5*gammaaa
     3         *t10/rhoa**5.333333333333333d+0)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = 1/rhob**2.6666666666666666d+0
            t10 = gammabb*t9
            t11 = 1.590448890258626d-3*t10**2.0d+0+1.0d+0
            t12 = 1/t11
            t13 = 4.4878396535373427d-3*gammabb*t9+1.0d+0
            t14 = 1/t13
            t15 = 3.249347363814451d-3*gammabb*t14*t9-5.68616292583346d-
     1         3*gammabb*t12*t9
            t16 = rhob**1.3333333333333333d+0
            t17 = gammaaa**2
            t18 = 1/t5**2
            t19 = 1/rhoa**6.333333333333333d+0
            t20 = 1/t3**2
            t21 = 1/rhoa**3.6666666666666664d+0
            t22 = gammabb**2
            t23 = 1/t13**2
            t24 = 1/rhob**6.333333333333333d+0
            t25 = 1/t11**2
            t26 = 1/rhob**3.6666666666666664d+0
            fnc(iq) = (-9.305257363491d-1*t7*t8-9.305257363491d-1*t15*t1
     1         6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491d-1*(-8.664926303505202d-3*
     1         gammaaa*t21*t6+1.5163101135555893d-2*gammaaa*t21*t4-4.823
     2         2274747848364d-5*t17*t19*t2*t20+3.8886799859849397d-5*t17
     3         *t18*t19)*t8-1.2407009817987999d+0*rhoa**3.33333333333333
     4         3d-1*t7)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t16*(-8.664926303505202
     1         d-3*gammabb*t14*t26+1.5163101135555893d-2*gammabb*t12*t26
     2         -4.8232274747848364d-5*t10*t22*t24*t25+3.8886799859849397
     3         d-5*t22*t23*t24)-1.2407009817987999d+0*rhob**3.3333333333
     4         33333d-1*t15)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*(3.24934
     1         7363814451d-3*t1*t6-5.68616292583346d-3*t1*t4+1.808710303
     2         0443136d-5*t17*t20/rhoa**8-1.4582549947443524d-5*gammaaa*
     3         t18/rhoa**5.333333333333333d+0)*t8*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t16*(3.2
     1         49347363814451d-3*t14*t9-5.68616292583346d-3*t12*t9+1.808
     2         7103030443136d-5*t22*t25/rhob**8-1.4582549947443524d-5*ga
     3         mmabb*t23/rhob**5.333333333333333d+0)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = gammaaa**2
            t10 = 1/t5**2
            t11 = 1/rhoa**6.333333333333333d+0
            t12 = 1/t3**2
            t13 = 1/rhoa**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RA) = -9.305257363491d-1*t8*(-4.8232274747848364d
     1         -5*t11*t12*t2*t9+3.8886799859849397d-5*t10*t11*t9-8.66492
     2         6303505202d-3*gammaaa*t13*t6+1.5163101135555893d-2*gammaa
     3         a*t13*t4)*wght-1.2407009817987999d+0*rhoa**3.333333333333
     4         333d-1*t7*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t8*(1.80
     1         87103030443136d-5*t12*t9/rhoa**8+3.249347363814451d-3*t1*
     2         t6-5.68616292583346d-3*t1*t4-1.4582549947443524d-5*gammaa
     3         a*t10/rhoa**5.333333333333333d+0)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = gammabb*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammabb*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammabb*t1*t6-5.68616292583346d-3*
     1         gammabb*t1*t4
            t8 = rhob**1.3333333333333333d+0
            t9 = gammabb**2
            t10 = 1/t5**2
            t11 = 1/rhob**6.333333333333333d+0
            t12 = 1/t3**2
            t13 = 1/rhob**3.6666666666666664d+0
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RB) = -9.305257363491d-1*t8*(-4.8232274747848364d
     1         -5*t11*t12*t2*t9+3.8886799859849397d-5*t10*t11*t9-8.66492
     2         6303505202d-3*gammabb*t13*t6+1.5163101135555893d-2*gammab
     3         b*t13*t4)*wght-1.2407009817987999d+0*rhob**3.333333333333
     4         333d-1*t7*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t8*(1.80
     1         87103030443136d-5*t12*t9/rhob**8+3.249347363814451d-3*t1*
     2         t6-5.68616292583346d-3*t1*t4-1.4582549947443524d-5*gammab
     3         b*t10/rhob**5.333333333333333d+0)*wght
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
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_2}\,
C>    \sigma_{\beta\beta}}\over{0.004487839653537343\,{\it t_2}
C>    \,\sigma_{\beta\beta}+1.0}}-{{0.00568616292583346\,{\it t_2}
C>    \,\sigma_{\beta\beta}}\over{0.001590448890258626\,\left({
C>    \it t_2}\,\sigma_{\beta\beta}\right)^{2.0}+1.0}}\right)
C>    -0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.004487839653537343\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}-{{0.00568616292583346\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}}\over{0.001590448890258626\,\left({
C>    \it t_1}\,\sigma_{\alpha\alpha}\right)^{2.0}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.004487839653537343\,{\it t_3}\,\sigma_{ss}+1.0}}
C>    -{{0.00568616292583346\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.001590448890258626\,\left({\it t_3}\,
C>    \sigma_{ss}\right)^{2.0}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131
C>    , 094103 (2009)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ssbd_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.590448890258626d-3*(gammaaa*t1)**2.0d+0+1.0d+0
            t3 = 1/t2
            t4 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t5 = 1/t4
            t6 = 3.249347363814451d-3*gammaaa*t1*t5-5.68616292583346d-3*
     1         gammaaa*t1*t3
            t7 = rhoa**1.3333333333333333d+0
            t8 = rhoa**3.333333333333333d-1
            t9 = gammaaa**3
            t10 = 1/t2**2
            t11 = 1/rhoa**9
            t12 = gammaaa**2
            t13 = 1/t4**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = -4.8232274747848364d-5*t10*t11*t9-8.664926303505202d-3
     1         *gammaaa*t15*t5+1.5163101135555893d-2*gammaaa*t15*t3+3.88
     2         86799859849397d-5*t12*t13*t14
            t17 = 1/rhoa**8
            t18 = 1/rhoa**5.333333333333333d+0
            t19 = 3.249347363814451d-3*t1*t5-5.68616292583346d-3*t1*t3-1
     1         .4582549947443524d-5*gammaaa*t13*t18+1.8087103030443136d-
     2         5*t10*t12*t17
            t20 = 1/t2**3
            t21 = 1/rhoa**10
            t22 = 1/t4**3
            t23 = 1/rhoa**4.666666666666667d+0
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t6*t7*wght
            Amat(iq,D1_RA) = (-1.2407009817987999d+0*t6*t8-9.30525736349
     1         1d-1*t16*t7)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t19*t7*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t7*(9.3076118618774
     1         66d-7*t21*t22*t9+5.627098720582309d-4*t10*t21*t9+3.177139
     2         644618574d-2*gammaaa*t23*t5-5.559803749703828d-2*gammaaa*
     3         t23*t3-8.182503237052221d-7*gammaaa**5*t20/rhoa**1.533333
     4         3333333335d+1-3.4998119873864453d-4*t12*t13/rhoa**7.33333
     5         3333333333d+0)-2.4814019635975998d+0*t16*t8-4.13566993932
     6         9333d-1*t6/rhoa**6.666666666666666d-1)*wght+Amat2(iq,D2_R
     7         A_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987999d+0*t19*t8-9.30525
     1         7363491d-1*(-8.664926303505202d-3*t15*t5+1.51631011355558
     2         93d-2*t15*t3-3.49035444820405d-7*t11*t12*t22+3.0684387138
     3         94583d-7*gammaaa**4*t20/rhoa**1.4333333333333334d+1+1.166
     4         6039957954819d-4*gammaaa*t13*t14-1.9292909899139346d-4*t1
     5         0*t11*t12)*t7)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t7*(-1.1506645177104688d-7*t20*t9/rhoa**1.3333333333333
     2         333d+1+1.3088829180765188d-7*gammaaa*t17*t22-2.9165099894
     3         88705d-5*t13*t18+5.426130909132941d-5*gammaaa*t10*t17)*wg
     4         ht
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = 1/rhob**2.6666666666666666d+0
            t10 = gammabb*t9
            t11 = 1.590448890258626d-3*t10**2.0d+0+1.0d+0
            t12 = 1/t11
            t13 = 4.4878396535373427d-3*gammabb*t9+1.0d+0
            t14 = 1/t13
            t15 = 3.249347363814451d-3*gammabb*t14*t9-5.68616292583346d-
     1         3*gammabb*t12*t9
            t16 = rhob**1.3333333333333333d+0
            t17 = rhoa**3.333333333333333d-1
            t18 = gammaaa**2
            t19 = 1/t5**2
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = 1/t3**2
            t22 = 1/rhoa**3.6666666666666664d+0
            t23 = -8.664926303505202d-3*gammaaa*t22*t6+1.516310113555589
     1         3d-2*gammaaa*t22*t4-4.8232274747848364d-5*t18*t2*t20*t21+
     2         3.8886799859849397d-5*t18*t19*t20
            t24 = rhob**3.333333333333333d-1
            t25 = gammabb**2
            t26 = 1/t13**2
            t27 = 1/rhob**6.333333333333333d+0
            t28 = 1/t11**2
            t29 = 1/rhob**3.6666666666666664d+0
            t30 = -8.664926303505202d-3*gammabb*t14*t29+1.51631011355558
     1         93d-2*gammabb*t12*t29-4.8232274747848364d-5*t10*t25*t27*t
     2         28+3.8886799859849397d-5*t25*t26*t27
            t31 = 1/rhoa**8
            t32 = 1/rhoa**5.333333333333333d+0
            t33 = -1.4582549947443524d-5*gammaaa*t19*t32
            t34 = -5.68616292583346d-3*t1*t4
            t35 = 3.249347363814451d-3*t1*t6
            t36 = 1/rhob**8
            t37 = 1/rhob**5.333333333333333d+0
            t38 = -1.4582549947443524d-5*gammabb*t26*t37
            t39 = -5.68616292583346d-3*t12*t9
            t40 = 3.249347363814451d-3*t14*t9
            t41 = 1/t3**3
            t42 = gammaaa**3
            t43 = 1/rhoa**10
            t44 = 1/t5**3
            t45 = 1/rhoa**7.333333333333333d+0
            t46 = 1/rhoa**4.666666666666667d+0
            t47 = 1/t11**3
            t48 = gammabb**3
            t49 = 1/rhob**10
            t50 = 1/t13**3
            t51 = 1/rhob**7.333333333333333d+0
            t52 = 1/rhob**4.666666666666667d+0
            t53 = 1/rhoa**9
            t54 = 1/rhob**9
            fnc(iq) = (-9.305257363491d-1*t7*t8-9.305257363491d-1*t15*t1
     1         6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491d-1*t23*t8-1.24070098179879
     1         99d+0*t17*t7)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t16*t30-1.2407009817987
     1         999d+0*t15*t24)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*(t35+t34
     1         +t33+1.8087103030443136d-5*t18*t21*t31)*t8*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t16*(t40
     1         +t39+t38+1.8087103030443136d-5*t25*t28*t36)*wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*(3.177139644618574d
     1         -2*gammaaa*t46*t6-5.559803749703828d-2*gammaaa*t4*t46+1.2
     2         861939932759564d-4*t18*t2*t21*t45-3.4998119873864453d-4*t
     3         18*t19*t45+9.307611861877466d-7*t42*t43*t44+4.34090472730
     4         6352d-4*t21*t42*t43-8.182503237052221d-7*gammaaa**4*t2*t4
     5         1/rhoa**1.2666666666666665d+1)*t8-4.135669939329333d-1*t7
     6         /rhoa**6.666666666666666d-1-2.4814019635975998d+0*t17*t23
     7         )*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491d-1*t16*(3.177139644618
     1         574d-2*gammabb*t14*t52-5.559803749703828d-2*gammabb*t12*t
     2         52+1.2861939932759564d-4*t10*t25*t28*t51-3.49981198738644
     3         53d-4*t25*t26*t51+9.307611861877466d-7*t48*t49*t50+4.3409
     4         04727306352d-4*t28*t48*t49-8.182503237052221d-7*gammabb**
     5         4*t10*t47/rhob**1.2666666666666665d+1)-2.4814019635975998
     6         d+0*t24*t30-4.135669939329333d-1*t15/rhob**6.666666666666
     7         666d-1)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*(-8.66492630350520
     1         2d-3*t22*t6-3.49035444820405d-7*t18*t44*t53-1.44696824243
     2         5451d-4*t18*t21*t53+3.068438713894583d-7*t2*t41*t42/rhoa*
     3         *1.1666666666666665d+1+1.5163101135555893d-2*t22*t4-4.823
     4         227474784837d-5*gammaaa*t2*t20*t21+1.1666039957954819d-4*
     5         gammaaa*t19*t20)*t8-1.2407009817987999d+0*t17*(t35+t34+t3
     6         3+1.8087103030443136d-5*gammaaa*t2*t21*t32))*wght+Cmat2(i
     7         q,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363491d-1*t16*(-3.4903544482
     1         0405d-7*t25*t50*t54-1.446968242435451d-4*t25*t28*t54+3.06
     2         8438713894583d-7*t10*t47*t48/rhob**1.1666666666666665d+1-
     3         8.664926303505202d-3*t14*t29+1.5163101135555893d-2*t12*t2
     4         9-4.823227474784837d-5*gammabb*t10*t27*t28+1.166603995795
     5         4819d-4*gammabb*t26*t27)-1.2407009817987999d+0*t24*(t40+t
     6         39+t38+1.8087103030443136d-5*gammabb*t10*t28*t37))*wght+C
     7         mat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*(1.3088829180765188d-7*gammaaa*t31*t44-1.15066451771046
     2         88d-7*t41*t42/rhoa**1.3333333333333333d+1-2.9165099894887
     3         05d-5*t19*t32+5.426130909132941d-5*gammaaa*t21*t31)*t8*wg
     4         ht
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t16*(1.3088829180765188d-7*gammabb*t36*t50-1.1506645177
     2         104688d-7*t47*t48/rhob**1.3333333333333333d+1-2.916509989
     3         488705d-5*t26*t37+5.426130909132941d-5*gammabb*t28*t36)*w
     4         ght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = rhoa**3.333333333333333d-1
            t10 = gammaaa**2
            t11 = 1/t5**2
            t12 = 1/rhoa**6.333333333333333d+0
            t13 = 1/t3**2
            t14 = 1/rhoa**3.6666666666666664d+0
            t15 = -8.664926303505202d-3*gammaaa*t14*t6+1.516310113555589
     1         3d-2*gammaaa*t14*t4-4.8232274747848364d-5*t10*t12*t13*t2+
     2         3.8886799859849397d-5*t10*t11*t12
            t16 = 1/rhoa**8
            t17 = 1/rhoa**5.333333333333333d+0
            t18 = -1.4582549947443524d-5*gammaaa*t11*t17
            t19 = -5.68616292583346d-3*t1*t4
            t20 = 3.249347363814451d-3*t1*t6
            t21 = 1/t3**3
            t22 = gammaaa**3
            t23 = 1/rhoa**10
            t24 = 1/t5**3
            t25 = 1/rhoa**7.333333333333333d+0
            t26 = 1/rhoa**4.666666666666667d+0
            t27 = 1/rhoa**9
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RA) = -1.2407009817987999d+0*t7*t9*wght-9.3052573
     1         63491d-1*t15*t8*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*(t20+t19
     1         +t18+1.8087103030443136d-5*t10*t13*t16)*t8*wght
            Amat2(iq,D2_RA_RA) = -2.4814019635975998d+0*t15*t9*wght-9.30
     1         5257363491d-1*(3.177139644618574d-2*gammaaa*t26*t6-5.5598
     2         03749703828d-2*gammaaa*t26*t4+1.2861939932759564d-4*t10*t
     3         13*t2*t25-3.4998119873864453d-4*t10*t11*t25+9.30761186187
     4         7466d-7*t22*t23*t24+4.340904727306352d-4*t13*t22*t23-8.18
     5         2503237052221d-7*gammaaa**4*t2*t21/rhoa**1.26666666666666
     6         65d+1)*t8*wght-4.135669939329333d-1*t7*wght/rhoa**6.66666
     7         6666666666d-1+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -1.2407009817987999d+0*(t20+1.80871030
     1         30443136d-5*gammaaa*t13*t17*t2+t19+t18)*t9*wght-9.3052573
     2         63491d-1*(-8.664926303505202d-3*t14*t6+1.5163101135555893
     3         d-2*t14*t4-3.49035444820405d-7*t10*t24*t27-1.446968242435
     4         451d-4*t10*t13*t27+3.068438713894583d-7*t2*t21*t22/rhoa**
     5         1.1666666666666665d+1-4.823227474784837d-5*gammaaa*t12*t1
     6         3*t2+1.1666039957954819d-4*gammaaa*t11*t12)*t8*wght+Cmat2
     7         (iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*(1.3088829180765188d-7*gammaaa*t16*t24-1.15066451771046
     2         88d-7*t21*t22/rhoa**1.3333333333333333d+1-2.9165099894887
     3         05d-5*t11*t17+5.426130909132941d-5*gammaaa*t13*t16)*t8*wg
     4         ht
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = gammabb*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammabb*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammabb*t1*t6-5.68616292583346d-3*
     1         gammabb*t1*t4
            t8 = rhob**1.3333333333333333d+0
            t9 = rhob**3.333333333333333d-1
            t10 = gammabb**2
            t11 = 1/t5**2
            t12 = 1/rhob**6.333333333333333d+0
            t13 = 1/t3**2
            t14 = 1/rhob**3.6666666666666664d+0
            t15 = -8.664926303505202d-3*gammabb*t14*t6+1.516310113555589
     1         3d-2*gammabb*t14*t4-4.8232274747848364d-5*t10*t12*t13*t2+
     2         3.8886799859849397d-5*t10*t11*t12
            t16 = 1/rhob**8
            t17 = 1/rhob**5.333333333333333d+0
            t18 = -1.4582549947443524d-5*gammabb*t11*t17
            t19 = -5.68616292583346d-3*t1*t4
            t20 = 3.249347363814451d-3*t1*t6
            t21 = 1/t3**3
            t22 = gammabb**3
            t23 = 1/rhob**10
            t24 = 1/t5**3
            t25 = 1/rhob**7.333333333333333d+0
            t26 = 1/rhob**4.666666666666667d+0
            t27 = 1/rhob**9
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RB) = -1.2407009817987999d+0*t7*t9*wght-9.3052573
     1         63491d-1*t15*t8*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*(t20+t19
     1         +t18+1.8087103030443136d-5*t10*t13*t16)*t8*wght
            Amat2(iq,D2_RB_RB) = -2.4814019635975998d+0*t15*t9*wght-9.30
     1         5257363491d-1*(3.177139644618574d-2*gammabb*t26*t6-5.5598
     2         03749703828d-2*gammabb*t26*t4+1.2861939932759564d-4*t10*t
     3         13*t2*t25-3.4998119873864453d-4*t10*t11*t25+9.30761186187
     4         7466d-7*t22*t23*t24+4.340904727306352d-4*t13*t22*t23-8.18
     5         2503237052221d-7*gammabb**4*t2*t21/rhob**1.26666666666666
     6         65d+1)*t8*wght-4.135669939329333d-1*t7*wght/rhob**6.66666
     7         6666666666d-1+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -1.2407009817987999d+0*(t20+1.80871030
     1         30443136d-5*gammabb*t13*t17*t2+t19+t18)*t9*wght-9.3052573
     2         63491d-1*(-8.664926303505202d-3*t14*t6+1.5163101135555893
     3         d-2*t14*t4-3.49035444820405d-7*t10*t24*t27-1.446968242435
     4         451d-4*t10*t13*t27+3.068438713894583d-7*t2*t21*t22/rhob**
     5         1.1666666666666665d+1-4.823227474784837d-5*gammabb*t12*t1
     6         3*t2+1.1666039957954819d-4*gammabb*t11*t12)*t8*wght+Cmat2
     7         (iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*(1.3088829180765188d-7*gammabb*t16*t24-1.15066451771046
     2         88d-7*t21*t22/rhob**1.3333333333333333d+1-2.9165099894887
     3         05d-5*t11*t17+5.426130909132941d-5*gammabb*t13*t16)*t8*wg
     4         ht
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
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_2}\,
C>    \sigma_{\beta\beta}}\over{0.004487839653537343\,{\it t_2}
C>    \,\sigma_{\beta\beta}+1.0}}-{{0.00568616292583346\,{\it t_2}
C>    \,\sigma_{\beta\beta}}\over{0.001590448890258626\,\left({
C>    \it t_2}\,\sigma_{\beta\beta}\right)^{2.0}+1.0}}\right)
C>    -0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}}\over{0.004487839653537343\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0}}-{{0.00568616292583346\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}}\over{0.001590448890258626\,\left({
C>    \it t_1}\,\sigma_{\alpha\alpha}\right)^{2.0}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.003249347363814451\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.004487839653537343\,{\it t_3}\,\sigma_{ss}+1.0}}
C>    -{{0.00568616292583346\,{\it t_3}\,\sigma_{ss}}
C>    \over{0.001590448890258626\,\left({\it t_3}\,
C>    \sigma_{ss}\right)^{2.0}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, M Sola, FM Bickelhaupt, J.Chem.Phys. 131
C>    , 094103 (2009)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3213193 "> 
C> 10.1063/1.3213193 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ssbd_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 1.590448890258626d-3*(gammaaa*t1)**2.0d+0+1.0d+0
            t3 = 1/t2
            t4 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t5 = 1/t4
            t6 = 3.249347363814451d-3*gammaaa*t1*t5-5.68616292583346d-3*
     1         gammaaa*t1*t3
            t7 = rhoa**1.3333333333333333d+0
            t8 = rhoa**3.333333333333333d-1
            t9 = gammaaa**3
            t10 = 1/t2**2
            t11 = 1/rhoa**9
            t12 = gammaaa**2
            t13 = 1/t4**2
            t14 = 1/rhoa**6.333333333333333d+0
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = -4.8232274747848364d-5*t10*t11*t9-8.664926303505202d-3
     1         *gammaaa*t15*t5+1.5163101135555893d-2*gammaaa*t15*t3+3.88
     2         86799859849397d-5*t12*t13*t14
            t17 = 1/rhoa**8
            t18 = 1/rhoa**5.333333333333333d+0
            t19 = 3.249347363814451d-3*t1*t5-5.68616292583346d-3*t1*t3-1
     1         .4582549947443524d-5*gammaaa*t13*t18+1.8087103030443136d-
     2         5*t10*t12*t17
            t20 = 1/rhoa**6.666666666666666d-1
            t21 = gammaaa**5
            t22 = 1/t2**3
            t23 = 1/rhoa**1.5333333333333335d+1
            t24 = 1/rhoa**10
            t25 = 1/t4**3
            t26 = 1/rhoa**7.333333333333333d+0
            t27 = 1/rhoa**4.666666666666667d+0
            t28 = 9.307611861877466d-7*t24*t25*t9+5.627098720582309d-4*t
     1         10*t24*t9+3.177139644618574d-2*gammaaa*t27*t5-5.559803749
     2         703828d-2*gammaaa*t27*t3-3.4998119873864453d-4*t12*t13*t2
     3         6-8.182503237052221d-7*t21*t22*t23
            t29 = gammaaa**4
            t30 = 1/rhoa**1.4333333333333334d+1
            t31 = -8.664926303505202d-3*t15*t5+3.068438713894583d-7*t22*
     1         t29*t30+1.5163101135555893d-2*t15*t3-3.49035444820405d-7*
     2         t11*t12*t25+1.1666039957954819d-4*gammaaa*t13*t14-1.92929
     3         09899139346d-4*t10*t11*t12
            t32 = 1/rhoa**1.3333333333333333d+1
            t33 = -1.1506645177104688d-7*t22*t32*t9+1.3088829180765188d-
     1         7*gammaaa*t17*t25-2.916509989488705d-5*t13*t18+5.42613090
     2         9132941d-5*gammaaa*t10*t17
            t34 = 1/t2**4
            t35 = 1/t4**4
            t36 = 1/rhoa**11
            t37 = 1/rhoa**5.666666666666667d+0
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t6*t7*wght
            Amat(iq,D1_RA) = (-1.2407009817987999d+0*t6*t8-9.30525736349
     1         1d-1*t16*t7)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t19*t7*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-2.4814019635975998d+0*t16*t8-9.305257
     1         363491d-1*t28*t7-4.135669939329333d-1*t20*t6)*wght+Amat2(
     2         iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.2407009817987999d+0*t19*t8-9.30525
     1         7363491d-1*t31*t7)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*t33*t7*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491d-1*t7*(-1.768446253
     1         7567186d-5*t25*t36*t9-6.098703184783492d-3*t10*t36*t9-1.4
     2         826651674886682d-1*gammaaa*t37*t5+2.59457508319512d-1*gam
     3         maaa*t3*t37+3.341685567477458d-8*t29*t35/rhoa**1.36666666
     4         66666666d+1-2.082216510865171d-8*gammaaa**7*t34/rhoa**2.1
     5         666666666666667d+1+2.2092758740040996d-5*t21*t22/rhoa**1.
     6         6333333333333333d+1+2.946755278268587d-3*t12*t13/rhoa**8.
     7         333333333333334d+0)-3.7221029453964d+0*t28*t8+2.757113292
     8         886222d-1*t6/rhoa**1.6666666666666669d+0-1.24070098179879
     9         99d+0*t16*t20)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491d-1*t7*(-1.25313208
     1         7804047d-8*t35*t9/rhoa**1.2666666666666665d+1+3.177139644
     2         618574d-2*t27*t5+7.808311915744394d-9*gammaaa**6*t34/rhoa
     3         **2.0666666666666667d+1-5.559803749703828d-2*t27*t3-7.671
     4         096784736456d-6*t22*t23*t29-8.425473302967369d-4*gammaaa*
     5         t13*t26+5.933602561946885d-6*t12*t24*t25+1.86498129025013
     6         64d-3*t10*t12*t24)-2.4814019635975998d+0*t31*t8-4.1356699
     7         39329333d-1*t19*t20)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-9.305257363491d-1*t7*(2.45475097
     1         11156665d-6*t22*t30*t9+4.699245329265176d-9*t12*t35/rhoa*
     2         *1.1666666666666665d+1-2.9281169684041475d-9*t21*t34/rhoa
     3         **1.9666666666666668d+1-1.745177224102025d-6*gammaaa*t11*
     4         t25+1.5554719943939757d-4*t13*t14-4.3409047273063534d-4*g
     5         ammaaa*t10*t11)-1.2407009817987999d+0*t33*t8)*wght+Cmat3(
     6         iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*(-1.762216998474441d-9*gammaaa*t35/rhoa**1.0666
     2         666666666666d+1+1.0980438631515553d-9*t29*t34/rhoa**1.866
     3         6666666666668d+1-6.903987106262813d-7*t12*t22*t32+3.92664
     4         8754229556d-7*t17*t25+5.426130909132941d-5*t10*t17)*t7*wg
     5         ht
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = 1/rhob**2.6666666666666666d+0
            t10 = gammabb*t9
            t11 = 1.590448890258626d-3*t10**2.0d+0+1.0d+0
            t12 = 1/t11
            t13 = 4.4878396535373427d-3*gammabb*t9+1.0d+0
            t14 = 1/t13
            t15 = 3.249347363814451d-3*gammabb*t14*t9-5.68616292583346d-
     1         3*gammabb*t12*t9
            t16 = rhob**1.3333333333333333d+0
            t17 = rhoa**3.333333333333333d-1
            t18 = gammaaa**2
            t19 = 1/t5**2
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = 3.8886799859849397d-5*t18*t19*t20
            t22 = 1/t3**2
            t23 = 1/rhoa**3.6666666666666664d+0
            t24 = 1.5163101135555893d-2*gammaaa*t23*t4
            t25 = -8.664926303505202d-3*gammaaa*t23*t6
            t26 = t25+t24-4.8232274747848364d-5*t18*t2*t20*t22+t21
            t27 = rhob**3.333333333333333d-1
            t28 = gammabb**2
            t29 = 1/t13**2
            t30 = 1/rhob**6.333333333333333d+0
            t31 = 3.8886799859849397d-5*t28*t29*t30
            t32 = 1/t11**2
            t33 = 1/rhob**3.6666666666666664d+0
            t34 = 1.5163101135555893d-2*gammabb*t12*t33
            t35 = -8.664926303505202d-3*gammabb*t14*t33
            t36 = t35+t34-4.8232274747848364d-5*t10*t28*t30*t32+t31
            t37 = 1/rhoa**8
            t38 = 1/rhoa**5.333333333333333d+0
            t39 = -1.4582549947443524d-5*gammaaa*t19*t38
            t40 = -5.68616292583346d-3*t1*t4
            t41 = 3.249347363814451d-3*t1*t6
            t42 = t41+t40+t39+1.8087103030443136d-5*t18*t22*t37
            t43 = 1/rhob**8
            t44 = 1/rhob**5.333333333333333d+0
            t45 = -1.4582549947443524d-5*gammabb*t29*t44
            t46 = -5.68616292583346d-3*t12*t9
            t47 = 3.249347363814451d-3*t14*t9
            t48 = t47+t46+t45+1.8087103030443136d-5*t28*t32*t43
            t49 = 1/rhoa**6.666666666666666d-1
            t50 = gammaaa**4
            t51 = 1/t3**3
            t52 = 1/rhoa**1.2666666666666665d+1
            t53 = gammaaa**3
            t54 = 1/rhoa**10
            t55 = 1/t5**3
            t56 = 1/rhoa**7.333333333333333d+0
            t57 = 1/rhoa**4.666666666666667d+0
            t58 = 3.177139644618574d-2*gammaaa*t57*t6-5.559803749703828d
     1         -2*gammaaa*t4*t57+1.2861939932759564d-4*t18*t2*t22*t56-3.
     2         4998119873864453d-4*t18*t19*t56+9.307611861877466d-7*t53*
     3         t54*t55+4.340904727306352d-4*t22*t53*t54-8.18250323705222
     4         1d-7*t2*t50*t51*t52
            t59 = 1/rhob**6.666666666666666d-1
            t60 = gammabb**4
            t61 = 1/t11**3
            t62 = 1/rhob**1.2666666666666665d+1
            t63 = gammabb**3
            t64 = 1/rhob**10
            t65 = 1/t13**3
            t66 = 1/rhob**7.333333333333333d+0
            t67 = 1/rhob**4.666666666666667d+0
            t68 = 3.177139644618574d-2*gammabb*t14*t67-5.559803749703828
     1         d-2*gammabb*t12*t67+1.2861939932759564d-4*t10*t28*t32*t66
     2         -3.4998119873864453d-4*t28*t29*t66+9.307611861877466d-7*t
     3         63*t64*t65+4.340904727306352d-4*t32*t63*t64-8.18250323705
     4         2221d-7*t10*t60*t61*t62
            t69 = 1/rhoa**1.1666666666666665d+1
            t70 = 1/rhoa**9
            t71 = -3.49035444820405d-7*t18*t55*t70-1.446968242435451d-4*
     1         t18*t22*t70+3.068438713894583d-7*t2*t51*t53*t69-8.6649263
     2         03505202d-3*t23*t6+1.5163101135555893d-2*t23*t4-4.8232274
     3         74784837d-5*gammaaa*t2*t20*t22+1.1666039957954819d-4*gamm
     4         aaa*t19*t20
            t72 = 1/rhob**1.1666666666666665d+1
            t73 = 1/rhob**9
            t74 = -3.49035444820405d-7*t28*t65*t73-1.446968242435451d-4*
     1         t28*t32*t73+3.068438713894583d-7*t10*t61*t63*t72-8.664926
     2         303505202d-3*t14*t33+1.5163101135555893d-2*t12*t33-4.8232
     3         27474784837d-5*gammabb*t10*t30*t32+1.1666039957954819d-4*
     4         gammabb*t29*t30
            t75 = 1/rhoa**1.3333333333333333d+1
            t76 = 1.3088829180765188d-7*gammaaa*t37*t55
            t77 = -2.916509989488705d-5*t19*t38
            t78 = 1/rhob**1.3333333333333333d+1
            t79 = 1.3088829180765188d-7*gammabb*t43*t65
            t80 = -2.916509989488705d-5*t29*t44
            t81 = 1/t3**4
            t82 = gammaaa**5
            t83 = 1/t5**4
            t84 = 1/rhoa**11
            t85 = 1/rhoa**8.333333333333334d+0
            t86 = 1/rhoa**5.666666666666667d+0
            t87 = 1/t11**4
            t88 = gammabb**5
            t89 = 1/t13**4
            t90 = 1/rhob**11
            t91 = 1/rhob**8.333333333333334d+0
            t92 = 1/rhob**5.666666666666667d+0
            t93 = 1/rhoa**1.0666666666666666d+1
            t94 = 1/rhob**1.0666666666666666d+1
            fnc(iq) = (-9.305257363491d-1*t7*t8-9.305257363491d-1*t15*t1
     1         6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-9.305257363491d-1*t26*t8-1.24070098179879
     1         99d+0*t17*t7)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363491d-1*t16*t36-1.2407009817987
     1         999d+0*t15*t27)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t42*t8*w
     1         ght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t16*t48*
     1         wght
            Amat2(iq,D2_RA_RA) = (-9.305257363491d-1*t58*t8-4.1356699393
     1         29333d-1*t49*t7-2.4814019635975998d+0*t17*t26)*wght+Amat2
     2         (iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-9.305257363491d-1*t16*t68-4.135669939
     1         329333d-1*t15*t59-2.4814019635975998d+0*t27*t36)*wght+Ama
     2         t2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-9.305257363491d-1*t71*t8-1.240700981
     1         7987999d+0*t17*(t41+t40+t39+1.8087103030443136d-5*gammaaa
     2         *t2*t22*t38))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-9.305257363491d-1*t16*t74-1.24070098
     1         17987999d+0*t27*(t47+t46+t45+1.8087103030443136d-5*gammab
     2         b*t10*t32*t44))*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*(t77+t76-1.1506645177104688d-7*t51*t53*t75+5.4261309091
     2         32941d-5*gammaaa*t22*t37)*t8*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*t16*(t80+t79-1.1506645177104688d-7*t61*t63*t78+5.426130
     2         909132941d-5*gammabb*t32*t43)*wght
            Amat3(iq,D3_RA_RA_RA) = (-9.305257363491d-1*t8*(-1.482665167
     1         4886682d-1*gammaaa*t6*t86+2.59457508319512d-1*gammaaa*t4*
     2         t86-4.7160446420118396d-4*t18*t2*t22*t85+2.94675527826858
     3         7d-3*t18*t19*t85-1.7684462537567186d-5*t53*t55*t84-5.6270
     4         98720582308d-3*t22*t53*t84+3.341685567477458d-8*t50*t83/r
     5         hoa**1.3666666666666666d+1+2.2092758740040996d-5*t51*t82/
     6         rhoa**1.6333333333333333d+1-2.082216510865171d-8*gammaaa*
     7         *6*t2*t81/rhoa**19)-1.2407009817987999d+0*t49*(-4.8232274
     8         747848364d-5*t22*t53*t70+t25+t24+t21)+2.757113292886222d-
     9         1*t7/rhoa**1.6666666666666669d+0-3.7221029453964d+0*t17*t
     :         58)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-9.305257363491d-1*t16*(-1.48266516
     1         74886682d-1*gammabb*t14*t92+2.59457508319512d-1*gammabb*t
     2         12*t92-4.7160446420118396d-4*t10*t28*t32*t91+2.9467552782
     3         68587d-3*t28*t29*t91-1.7684462537567186d-5*t63*t65*t90-5.
     4         627098720582308d-3*t32*t63*t90+3.341685567477458d-8*t60*t
     5         89/rhob**1.3666666666666666d+1+2.2092758740040996d-5*t61*
     6         t88/rhob**1.6333333333333333d+1-2.082216510865171d-8*gamm
     7         abb**6*t10*t87/rhob**19)-1.2407009817987999d+0*t59*(-4.82
     8         32274747848364d-5*t32*t63*t73+t35+t34+t31)-3.722102945396
     9         4d+0*t27*t68+2.757113292886222d-1*t15/rhob**1.66666666666
     :         66669d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-9.305257363491d-1*t8*(-1.25313208
     1         7804047d-8*t52*t53*t83+7.808311915744394d-9*t2*t81*t82/rh
     2         oa**18+3.177139644618574d-2*t57*t6-5.559803749703828d-2*t
     3         4*t57+1.7685167407544403d-4*gammaaa*t2*t22*t56-8.42547330
     4         2967369d-4*gammaaa*t19*t56+5.933602561946885d-6*t18*t54*t
     5         55+1.6881296161746925d-3*t18*t22*t54-7.671096784736456d-6
     6         *t50*t51/rhoa**1.5333333333333335d+1)-2.4814019635975998d
     7         +0*t17*t71-4.135669939329333d-1*t42*t49)*wght+Cmat3(iq,D3
     8         _RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-9.305257363491d-1*t16*(-1.2531320
     1         87804047d-8*t62*t63*t89+7.808311915744394d-9*t10*t87*t88/
     2         rhob**18+3.177139644618574d-2*t14*t67-5.559803749703828d-
     3         2*t12*t67+1.7685167407544403d-4*gammabb*t10*t32*t66-8.425
     4         473302967369d-4*gammabb*t29*t66+5.933602561946885d-6*t28*
     5         t64*t65+1.6881296161746925d-3*t28*t32*t64-7.6710967847364
     6         56d-6*t60*t61/rhob**1.5333333333333335d+1)-2.481401963597
     7         5998d+0*t27*t74-4.135669939329333d-1*t48*t59)*wght+Cmat3(
     8         iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-1.2407009817987999d+0*t17*(-1.15
     1         06645177104688d-7*t18*t2*t51*t93+t77+t76+1.80871030304431
     2         36d-5*t2*t22*t38+3.617420606088627d-5*gammaaa*t22*t37)-9.
     3         305257363491d-1*t8*(4.699245329265176d-9*t18*t69*t83-2.92
     4         81169684041475d-9*t2*t50*t81/rhoa**17-1.745177224102025d-
     5         6*gammaaa*t55*t70-3.858581979827869d-4*gammaaa*t22*t70+2.
     6         4547509711156665d-6*t51*t53/rhoa**1.4333333333333334d+1-4
     7         .823227474784837d-5*t2*t20*t22+1.5554719943939757d-4*t19*
     8         t20))*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-1.2407009817987999d+0*t27*(-1.15
     1         06645177104688d-7*t10*t28*t61*t94+t80+t79+1.8087103030443
     2         136d-5*t10*t32*t44+3.617420606088627d-5*gammabb*t32*t43)-
     3         9.305257363491d-1*t16*(4.699245329265176d-9*t28*t72*t89-2
     4         .9281169684041475d-9*t10*t60*t87/rhob**17-1.7451772241020
     5         25d-6*gammabb*t65*t73-3.858581979827869d-4*gammabb*t32*t7
     6         3+2.4547509711156665d-6*t61*t63/rhob**1.4333333333333334d
     7         +1-4.823227474784837d-5*t10*t30*t32+1.5554719943939757d-4
     8         *t29*t30))*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*t8*(-1.762216998474441d-9*gammaaa*t83*t93+1.098
     2         0438631515553d-9*t50*t81/rhoa**1.8666666666666668d+1-6.90
     3         3987106262813d-7*t18*t51*t75+3.926648754229556d-7*t37*t55
     4         +5.426130909132941d-5*t22*t37)*wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491d-1*t16*(-1.762216998474441d-9*gammabb*t89*t94+1.09
     2         80438631515553d-9*t60*t87/rhob**1.8666666666666668d+1-6.9
     3         03987106262813d-7*t28*t61*t78+3.926648754229556d-7*t43*t6
     4         5+5.426130909132941d-5*t32*t43)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = gammaaa*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammaaa*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammaaa*t1*t6-5.68616292583346d-3*
     1         gammaaa*t1*t4
            t8 = rhoa**1.3333333333333333d+0
            t9 = rhoa**3.333333333333333d-1
            t10 = gammaaa**2
            t11 = 1/t5**2
            t12 = 1/rhoa**6.333333333333333d+0
            t13 = 3.8886799859849397d-5*t10*t11*t12
            t14 = 1/t3**2
            t15 = 1/rhoa**3.6666666666666664d+0
            t16 = 1.5163101135555893d-2*gammaaa*t15*t4
            t17 = -8.664926303505202d-3*gammaaa*t15*t6
            t18 = -4.8232274747848364d-5*t10*t12*t14*t2+t17+t16+t13
            t19 = 1/rhoa**8
            t20 = 1/rhoa**5.333333333333333d+0
            t21 = -1.4582549947443524d-5*gammaaa*t11*t20
            t22 = -5.68616292583346d-3*t1*t4
            t23 = 3.249347363814451d-3*t1*t6
            t24 = t23+t22+t21+1.8087103030443136d-5*t10*t14*t19
            t25 = 1/rhoa**6.666666666666666d-1
            t26 = gammaaa**4
            t27 = 1/t3**3
            t28 = 1/rhoa**1.2666666666666665d+1
            t29 = gammaaa**3
            t30 = 1/rhoa**10
            t31 = 1/t5**3
            t32 = 1/rhoa**7.333333333333333d+0
            t33 = 1/rhoa**4.666666666666667d+0
            t34 = 3.177139644618574d-2*gammaaa*t33*t6-5.559803749703828d
     1         -2*gammaaa*t33*t4+1.2861939932759564d-4*t10*t14*t2*t32-3.
     2         4998119873864453d-4*t10*t11*t32+9.307611861877466d-7*t29*
     3         t30*t31+4.340904727306352d-4*t14*t29*t30-8.18250323705222
     4         1d-7*t2*t26*t27*t28
            t35 = 1/rhoa**1.1666666666666665d+1
            t36 = 1/rhoa**9
            t37 = -8.664926303505202d-3*t15*t6+1.5163101135555893d-2*t15
     1         *t4-3.49035444820405d-7*t10*t31*t36-1.446968242435451d-4*
     2         t10*t14*t36+3.068438713894583d-7*t2*t27*t29*t35-4.8232274
     3         74784837d-5*gammaaa*t12*t14*t2+1.1666039957954819d-4*gamm
     4         aaa*t11*t12
            t38 = 1/rhoa**1.3333333333333333d+1
            t39 = 1.3088829180765188d-7*gammaaa*t19*t31
            t40 = -2.916509989488705d-5*t11*t20
            t41 = 1/t3**4
            t42 = gammaaa**5
            t43 = 1/t5**4
            t44 = 1/rhoa**11
            t45 = 1/rhoa**8.333333333333334d+0
            t46 = 1/rhoa**5.666666666666667d+0
            t47 = 1/rhoa**1.0666666666666666d+1
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RA) = -1.2407009817987999d+0*t7*t9*wght-9.3052573
     1         63491d-1*t18*t8*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363491d-1*t24*t8*w
     1         ght
            Amat2(iq,D2_RA_RA) = -2.4814019635975998d+0*t18*t9*wght-9.30
     1         5257363491d-1*t34*t8*wght-4.135669939329333d-1*t25*t7*wgh
     2         t+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -1.2407009817987999d+0*(t23+t22+t21+1.
     1         8087103030443136d-5*gammaaa*t14*t2*t20)*t9*wght-9.3052573
     2         63491d-1*t37*t8*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.305257363491d-
     1         1*(t40+t39-1.1506645177104688d-7*t27*t29*t38+5.4261309091
     2         32941d-5*gammaaa*t14*t19)*t8*wght
            Amat3(iq,D3_RA_RA_RA) = -3.7221029453964d+0*t34*t9*wght-9.30
     1         5257363491d-1*(-1.4826651674886682d-1*gammaaa*t46*t6+2.59
     2         457508319512d-1*gammaaa*t4*t46-4.7160446420118396d-4*t10*
     3         t14*t2*t45+2.946755278268587d-3*t10*t11*t45-1.76844625375
     4         67186d-5*t29*t31*t44-5.627098720582308d-3*t14*t29*t44+3.3
     5         41685567477458d-8*t26*t43/rhoa**1.3666666666666666d+1+2.2
     6         092758740040996d-5*t27*t42/rhoa**1.6333333333333333d+1-2.
     7         082216510865171d-8*gammaaa**6*t2*t41/rhoa**19)*t8*wght+2.
     8         757113292886222d-1*t7*wght/rhoa**1.6666666666666669d+0-1.
     9         2407009817987999d+0*t25*(-4.8232274747848364d-5*t14*t29*t
     :         36+t17+t16+t13)*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -2.4814019635975998d+0*t37*t9*wght-
     1         9.305257363491d-1*(3.177139644618574d-2*t33*t6-1.25313208
     2         7804047d-8*t28*t29*t43+7.808311915744394d-9*t2*t41*t42/rh
     3         oa**18-5.559803749703828d-2*t33*t4+1.7685167407544403d-4*
     4         gammaaa*t14*t2*t32-8.425473302967369d-4*gammaaa*t11*t32+5
     5         .933602561946885d-6*t10*t30*t31+1.6881296161746925d-3*t10
     6         *t14*t30-7.671096784736456d-6*t26*t27/rhoa**1.53333333333
     7         33335d+1)*t8*wght-4.135669939329333d-1*t24*t25*wght+Cmat3
     8         (iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -1.2407009817987999d+0*(-1.1506645
     1         177104688d-7*t10*t2*t27*t47+t40+t39+1.8087103030443136d-5
     2         *t14*t2*t20+3.617420606088627d-5*gammaaa*t14*t19)*t9*wght
     3         -9.305257363491d-1*(4.699245329265176d-9*t10*t35*t43-2.92
     4         81169684041475d-9*t2*t26*t41/rhoa**17-1.745177224102025d-
     5         6*gammaaa*t31*t36-3.858581979827869d-4*gammaaa*t14*t36+2.
     6         4547509711156665d-6*t27*t29/rhoa**1.4333333333333334d+1-4
     7         .823227474784837d-5*t12*t14*t2+1.5554719943939757d-4*t11*
     8         t12)*t8*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363491d-1*(-1.762216998474441d-9*gammaaa*t43*t47+1.098043
     2         8631515553d-9*t26*t41/rhoa**1.8666666666666668d+1-6.90398
     3         7106262813d-7*t10*t27*t38+3.926648754229556d-7*t19*t31+5.
     4         426130909132941d-5*t14*t19)*t8*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = gammabb*t1
            t3 = 1.590448890258626d-3*t2**2.0d+0+1.0d+0
            t4 = 1/t3
            t5 = 4.4878396535373427d-3*gammabb*t1+1.0d+0
            t6 = 1/t5
            t7 = 3.249347363814451d-3*gammabb*t1*t6-5.68616292583346d-3*
     1         gammabb*t1*t4
            t8 = rhob**1.3333333333333333d+0
            t9 = rhob**3.333333333333333d-1
            t10 = gammabb**2
            t11 = 1/t5**2
            t12 = 1/rhob**6.333333333333333d+0
            t13 = 3.8886799859849397d-5*t10*t11*t12
            t14 = 1/t3**2
            t15 = 1/rhob**3.6666666666666664d+0
            t16 = 1.5163101135555893d-2*gammabb*t15*t4
            t17 = -8.664926303505202d-3*gammabb*t15*t6
            t18 = -4.8232274747848364d-5*t10*t12*t14*t2+t17+t16+t13
            t19 = 1/rhob**8
            t20 = 1/rhob**5.333333333333333d+0
            t21 = -1.4582549947443524d-5*gammabb*t11*t20
            t22 = -5.68616292583346d-3*t1*t4
            t23 = 3.249347363814451d-3*t1*t6
            t24 = t23+t22+t21+1.8087103030443136d-5*t10*t14*t19
            t25 = 1/rhob**6.666666666666666d-1
            t26 = gammabb**4
            t27 = 1/t3**3
            t28 = 1/rhob**1.2666666666666665d+1
            t29 = gammabb**3
            t30 = 1/rhob**10
            t31 = 1/t5**3
            t32 = 1/rhob**7.333333333333333d+0
            t33 = 1/rhob**4.666666666666667d+0
            t34 = 3.177139644618574d-2*gammabb*t33*t6-5.559803749703828d
     1         -2*gammabb*t33*t4+1.2861939932759564d-4*t10*t14*t2*t32-3.
     2         4998119873864453d-4*t10*t11*t32+9.307611861877466d-7*t29*
     3         t30*t31+4.340904727306352d-4*t14*t29*t30-8.18250323705222
     4         1d-7*t2*t26*t27*t28
            t35 = 1/rhob**1.1666666666666665d+1
            t36 = 1/rhob**9
            t37 = -8.664926303505202d-3*t15*t6+1.5163101135555893d-2*t15
     1         *t4-3.49035444820405d-7*t10*t31*t36-1.446968242435451d-4*
     2         t10*t14*t36+3.068438713894583d-7*t2*t27*t29*t35-4.8232274
     3         74784837d-5*gammabb*t12*t14*t2+1.1666039957954819d-4*gamm
     4         abb*t11*t12
            t38 = 1/rhob**1.3333333333333333d+1
            t39 = 1.3088829180765188d-7*gammabb*t19*t31
            t40 = -2.916509989488705d-5*t11*t20
            t41 = 1/t3**4
            t42 = gammabb**5
            t43 = 1/t5**4
            t44 = 1/rhob**11
            t45 = 1/rhob**8.333333333333334d+0
            t46 = 1/rhob**5.666666666666667d+0
            t47 = 1/rhob**1.0666666666666666d+1
            fnc(iq) = fnc(iq)-9.305257363491d-1*t7*t8*wght
            Amat(iq,D1_RB) = -1.2407009817987999d+0*t7*t9*wght-9.3052573
     1         63491d-1*t18*t8*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363491d-1*t24*t8*w
     1         ght
            Amat2(iq,D2_RB_RB) = -2.4814019635975998d+0*t18*t9*wght-9.30
     1         5257363491d-1*t34*t8*wght-4.135669939329333d-1*t25*t7*wgh
     2         t+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -1.2407009817987999d+0*(t23+t22+t21+1.
     1         8087103030443136d-5*gammabb*t14*t2*t20)*t9*wght-9.3052573
     2         63491d-1*t37*t8*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.305257363491d-
     1         1*(t40+t39-1.1506645177104688d-7*t27*t29*t38+5.4261309091
     2         32941d-5*gammabb*t14*t19)*t8*wght
            Amat3(iq,D3_RB_RB_RB) = -3.7221029453964d+0*t34*t9*wght-9.30
     1         5257363491d-1*(-1.4826651674886682d-1*gammabb*t46*t6+2.59
     2         457508319512d-1*gammabb*t4*t46-4.7160446420118396d-4*t10*
     3         t14*t2*t45+2.946755278268587d-3*t10*t11*t45-1.76844625375
     4         67186d-5*t29*t31*t44-5.627098720582308d-3*t14*t29*t44+3.3
     5         41685567477458d-8*t26*t43/rhob**1.3666666666666666d+1+2.2
     6         092758740040996d-5*t27*t42/rhob**1.6333333333333333d+1-2.
     7         082216510865171d-8*gammabb**6*t2*t41/rhob**19)*t8*wght+2.
     8         757113292886222d-1*t7*wght/rhob**1.6666666666666669d+0-1.
     9         2407009817987999d+0*t25*(-4.8232274747848364d-5*t14*t29*t
     :         36+t17+t16+t13)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -2.4814019635975998d+0*t37*t9*wght-
     1         9.305257363491d-1*(3.177139644618574d-2*t33*t6-1.25313208
     2         7804047d-8*t28*t29*t43+7.808311915744394d-9*t2*t41*t42/rh
     3         ob**18-5.559803749703828d-2*t33*t4+1.7685167407544403d-4*
     4         gammabb*t14*t2*t32-8.425473302967369d-4*gammabb*t11*t32+5
     5         .933602561946885d-6*t10*t30*t31+1.6881296161746925d-3*t10
     6         *t14*t30-7.671096784736456d-6*t26*t27/rhob**1.53333333333
     7         33335d+1)*t8*wght-4.135669939329333d-1*t24*t25*wght+Cmat3
     8         (iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -1.2407009817987999d+0*(-1.1506645
     1         177104688d-7*t10*t2*t27*t47+t40+t39+1.8087103030443136d-5
     2         *t14*t2*t20+3.617420606088627d-5*gammabb*t14*t19)*t9*wght
     3         -9.305257363491d-1*(4.699245329265176d-9*t10*t35*t43-2.92
     4         81169684041475d-9*t2*t26*t41/rhob**17-1.745177224102025d-
     5         6*gammabb*t31*t36-3.858581979827869d-4*gammabb*t14*t36+2.
     6         4547509711156665d-6*t27*t29/rhob**1.4333333333333334d+1-4
     7         .823227474784837d-5*t12*t14*t2+1.5554719943939757d-4*t11*
     8         t12)*t8*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363491d-1*(-1.762216998474441d-9*gammabb*t43*t47+1.098043
     2         8631515553d-9*t26*t41/rhob**1.8666666666666668d+1-6.90398
     3         7106262813d-7*t10*t27*t38+3.926648754229556d-7*t19*t31+5.
     4         426130909132941d-5*t14*t19)*t8*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
