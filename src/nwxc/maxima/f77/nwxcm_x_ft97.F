C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_ft97.F 
C> The nwxcm_x_ft97 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_ft97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \left| \rho_\alpha\right| \\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{2.666666666666667}}}\\\\ 
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,\sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,\sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}\over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{\it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}\over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}\over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}+1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}\,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{\it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}\over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,\sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
      subroutine nwxcm_x_ft97(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa+6.255746320200999d+6
            t2 = 1/t1
            t3 = 9.474169d-4*gammaaa*t2+2.9136439999999997d-3
            t4 = rhoa**1.3333333333333333d+0
            t5 = abs(rhoa)
            t6 = 1/t5**2.6666666666666666d+0
            t7 = t3**2.0d+0
            t8 = gammaaa**2
            t9 = 1/t5**5.333333333333333d+0
            t10 = (t8*t9+1.0d+0)**5.0d-1
            t11 = gammaaa*t6+t10
            t12 = log(t11)
            t13 = t12**2.0d+0
            t14 = (9.0d+0*gammaaa*t13*t6*t7+1.0d+0)**5.0d-1
            t15 = 1/t14
            t16 = 1/t10
            t17 = 1/t5**4.666666666666666d+0
            t18 = 1/t11
            t19 = 1/t14**3
            t20 = 9.474169d-4*t2-9.474169d-4*gammaaa/t1**2
            t21 = t15*t6*(-2.0d+0*t3*t4*wght-2.0d+0*gammaaa*t20*t4*wght)
     1         +1.0d+0*gammaaa*t19*t3*t4*t6*(1.8d+1*gammaaa*t12*t18*t6*t
     2         7*(gammaaa*t16*t9+t6)+9.0d+0*t13*t6*t7+1.8d+1*gammaaa*t13
     3         *t20*t3*t6)*wght
            fnc(iq) = fnc(iq)-2.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 1.0d+0*gammaaa*t19*t3*t4*t6*(1.8d+1*gammaaa
     1         *t12*t18*t6*t7*(-2.6666666666666666d+0*rhoa*t16*t8/t5**7.
     2         333333333333333d+0-2.6666666666666666d+0*gammaaa*rhoa*t17
     3         )-2.4d+1*gammaaa*rhoa*t13*t17*t7)*wght-2.6666666666666666
     4         d+0*gammaaa*rhoa**3.333333333333333d-1*t15*t3*t6*wght+5.3
     5         33333333333333d+0*gammaaa*rhoa**2.3333333333333334d+0*t15
     6         *t17*t3*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t21+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t21+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammaaa+6.255746320200999d+6
            t2 = 1/t1
            t3 = 9.474169d-4*gammaaa*t2+2.9136439999999997d-3
            t4 = rhoa**1.3333333333333333d+0
            t5 = abs(rhoa)
            t6 = 1/t5**2.6666666666666666d+0
            t7 = t3**2.0d+0
            t8 = gammaaa**2
            t9 = 1/t5**5.333333333333333d+0
            t10 = (t8*t9+1.0d+0)**5.0d-1
            t11 = gammaaa*t6+t10
            t12 = log(t11)
            t13 = t12**2.0d+0
            t14 = (9.0d+0*gammaaa*t13*t6*t7+1.0d+0)**5.0d-1
            t15 = 1/t14
            t16 = gammabb+6.255746320200999d+6
            t17 = 1/t16
            t18 = 9.474169d-4*gammabb*t17+2.9136439999999997d-3
            t19 = rhob**1.3333333333333333d+0
            t20 = abs(rhob)
            t21 = 1/t20**2.6666666666666666d+0
            t22 = t18**2.0d+0
            t23 = gammabb**2
            t24 = 1/t20**5.333333333333333d+0
            t25 = (t23*t24+1.0d+0)**5.0d-1
            t26 = t25+gammabb*t21
            t27 = log(t26)
            t28 = t27**2.0d+0
            t29 = (9.0d+0*gammabb*t21*t22*t28+1.0d+0)**5.0d-1
            t30 = 1/t29
            t31 = 1/t10
            t32 = 1/t5**4.666666666666666d+0
            t33 = 1/t11
            t34 = 1/t14**3
            t35 = 1/t25
            t36 = 1/t20**4.666666666666666d+0
            t37 = 1/t26
            t38 = 1/t29**3
            t39 = 9.474169d-4*t2-9.474169d-4*gammaaa/t1**2
            t40 = 9.474169d-4*t17-9.474169d-4*gammabb/t16**2
            fnc(iq) = 1.0d+0*(-gammaaa*t15*t3*t4*t6-gammabb*t18*t19*t21*
     1         t30)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t3*t34*t4*t6*(1.8d+1
     1         *gammaaa*t12*t33*t6*t7*(-2.6666666666666666d+0*rhoa*t31*t
     2         8/t5**7.333333333333333d+0-2.6666666666666666d+0*gammaaa*
     3         rhoa*t32)-2.4d+1*gammaaa*rhoa*t13*t32*t7)-1.3333333333333
     4         333d+0*gammaaa*rhoa**3.333333333333333d-1*t15*t3*t6+2.666
     5         6666666666666d+0*gammaaa*rhoa**2.3333333333333334d+0*t15*
     6         t3*t32)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(5.0d-1*gammabb*t18*t19*t21*(1.8d+1*
     1         gammabb*t21*t22*t27*(-2.6666666666666666d+0*gammabb*rhob*
     2         t36-2.6666666666666666d+0*rhob*t23*t35/t20**7.33333333333
     3         3333d+0)*t37-2.4d+1*gammabb*rhob*t22*t28*t36)*t38+2.66666
     4         66666666666d+0*gammabb*rhob**2.3333333333333334d+0*t18*t3
     5         0*t36-1.3333333333333333d+0*gammabb*rhob**3.3333333333333
     6         33d-1*t18*t21*t30)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 1.0d+0*(5.0d-1*gammaaa*t3*t34*t4*t6*(1.8d+
     1         1*gammaaa*t12*t33*t6*t7*(gammaaa*t31*t9+t6)+9.0d+0*t13*t6
     2         *t7+1.8d+1*gammaaa*t13*t3*t39*t6)+t15*(-gammaaa*t39*t4-t3
     3         *t4)*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*(5.0d-1*gammabb*t18*t19*t21*t38*(1.
     1         8d+1*gammabb*t18*t21*t28*t40+1.8d+1*gammabb*t21*t22*t27*(
     2         gammabb*t24*t35+t21)*t37+9.0d+0*t21*t22*t28)+t21*t30*(-ga
     3         mmabb*t19*t40-t18*t19))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa+6.255746320200999d+6
            t2 = 1/t1
            t3 = 9.474169d-4*gammaaa*t2+2.9136439999999997d-3
            t4 = rhoa**1.3333333333333333d+0
            t5 = abs(rhoa)
            t6 = 1/t5**2.6666666666666666d+0
            t7 = t3**2.0d+0
            t8 = gammaaa**2
            t9 = 1/t5**5.333333333333333d+0
            t10 = (t8*t9+1.0d+0)**5.0d-1
            t11 = gammaaa*t6+t10
            t12 = log(t11)
            t13 = t12**2.0d+0
            t14 = (9.0d+0*gammaaa*t13*t6*t7+1.0d+0)**5.0d-1
            t15 = 1/t14
            t16 = 1/t10
            t17 = 1/t5**4.666666666666666d+0
            t18 = 1/t11
            t19 = 1/t14**3
            t20 = 9.474169d-4*t2-9.474169d-4*gammaaa/t1**2
            fnc(iq) = fnc(iq)-1.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 5.0d-1*gammaaa*t19*t3*t4*t6*(1.8d+1*gammaaa
     1         *t12*t18*t6*t7*(-2.6666666666666666d+0*rhoa*t16*t8/t5**7.
     2         333333333333333d+0-2.6666666666666666d+0*gammaaa*rhoa*t17
     3         )-2.4d+1*gammaaa*rhoa*t13*t17*t7)*wght-1.3333333333333333
     4         d+0*gammaaa*rhoa**3.333333333333333d-1*t15*t3*t6*wght+2.6
     5         666666666666666d+0*gammaaa*rhoa**2.3333333333333334d+0*t1
     6         5*t17*t3*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t15*t6*(-1.0d+0*t3*t4*wght-1.0d+0*gammaaa*
     1         t20*t4*wght)+5.0d-1*gammaaa*t19*t3*t4*t6*(1.8d+1*gammaaa*
     2         t12*t18*t6*t7*(gammaaa*t16*t9+t6)+9.0d+0*t13*t6*t7+1.8d+1
     3         *gammaaa*t13*t20*t3*t6)*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb+6.255746320200999d+6
            t2 = 1/t1
            t3 = 9.474169d-4*gammabb*t2+2.9136439999999997d-3
            t4 = rhob**1.3333333333333333d+0
            t5 = abs(rhob)
            t6 = 1/t5**2.6666666666666666d+0
            t7 = t3**2.0d+0
            t8 = gammabb**2
            t9 = 1/t5**5.333333333333333d+0
            t10 = (t8*t9+1.0d+0)**5.0d-1
            t11 = gammabb*t6+t10
            t12 = log(t11)
            t13 = t12**2.0d+0
            t14 = (9.0d+0*gammabb*t13*t6*t7+1.0d+0)**5.0d-1
            t15 = 1/t14
            t16 = 1/t10
            t17 = 1/t5**4.666666666666666d+0
            t18 = 1/t11
            t19 = 1/t14**3
            t20 = 9.474169d-4*t2-9.474169d-4*gammabb/t1**2
            fnc(iq) = fnc(iq)-1.0d+0*gammabb*t15*t3*t4*t6*wght
            Amat(iq,D1_RB) = 5.0d-1*gammabb*t19*t3*t4*t6*(1.8d+1*gammabb
     1         *t12*t18*t6*t7*(-2.6666666666666666d+0*rhob*t16*t8/t5**7.
     2         333333333333333d+0-2.6666666666666666d+0*gammabb*rhob*t17
     3         )-2.4d+1*gammabb*rhob*t13*t17*t7)*wght-1.3333333333333333
     4         d+0*gammabb*rhob**3.333333333333333d-1*t15*t3*t6*wght+2.6
     5         666666666666666d+0*gammabb*rhob**2.3333333333333334d+0*t1
     6         5*t17*t3*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = t15*t6*(-1.0d+0*t3*t4*wght-1.0d+0*gammabb*
     1         t20*t4*wght)+5.0d-1*gammabb*t19*t3*t4*t6*(1.8d+1*gammabb*
     2         t12*t18*t6*t7*(gammabb*t16*t9+t6)+9.0d+0*t13*t6*t7+1.8d+1
     3         *gammabb*t13*t20*t3*t6)*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_ft97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \left| \rho_\alpha\right| \\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{2.666666666666667}}}\\\\ 
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,\sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,\sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}\over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{\it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}\over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}\over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}+1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}\,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{\it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}\over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,\sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
      subroutine nwxcm_x_ft97_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_ft97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \left| \rho_\alpha\right| \\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}^{2.666666666666667}}}\\\\ 
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,\sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,\sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}+0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}\over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{\it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}\over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}\over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}+1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}\,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{\it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}\over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,\sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
      subroutine nwxcm_x_ft97_d3(param,tol_rho,ipol,nq,wght, 
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
