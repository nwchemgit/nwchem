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
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}
C>    +6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}
C>    +0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}
C>    \over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}
C>    \left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{
C>    \it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,
C>    \sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}
C>    \over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}
C>    \over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{
C>    \it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{
C>    \sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}
C>    +1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{
C>    \it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}
C>    \over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,
C>    \sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="https://doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ft97(param,tol_rho,ipol,nq,wght, 
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
            fnc(iq) = fnc(iq)-2.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t19*t3*t4*t6*(1.8d+1
     1         *gammaaa*t12*t18*t6*t7*(-2.6666666666666666d+0*rhoa*t16*t
     2         8/t5**7.333333333333333d+0-2.6666666666666666d+0*gammaaa*
     3         rhoa*t17)-2.4d+1*gammaaa*rhoa*t13*t17*t7)-1.3333333333333
     4         333d+0*gammaaa*rhoa**3.333333333333333d-1*t15*t3*t6+2.666
     5         6666666666666d+0*gammaaa*rhoa**2.3333333333333334d+0*t15*
     6         t17*t3)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*(5.0d-1*gammaaa*t19*t3*t4*t6*(1.8d+
     1         1*gammaaa*t12*t18*t6*t7*(gammaaa*t16*t9+t6)+9.0d+0*t13*t6
     2         *t7+1.8d+1*gammaaa*t13*t20*t3*t6)+t15*(-t3*t4-gammaaa*t20
     3         *t4)*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
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
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}
C>    +6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}
C>    +0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}
C>    \over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}
C>    \left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{
C>    \it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,
C>    \sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}
C>    \over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}
C>    \over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{
C>    \it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{
C>    \sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}
C>    +1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{
C>    \it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}
C>    \over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,
C>    \sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="https://doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ft97_d2(param,tol_rho,ipol,nq,wght, 
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhoa*t16*t17*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammaaa*t12*t19*t20*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhoa**2.3333333333333334d+0
            t24 = rhoa**3.333333333333333d-1
            t25 = gammaaa*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammaaa*t26
            t28 = 1.8d+1*gammaaa*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t27*t3*t6
            t29 = -t3*t4-gammaaa*t27*t4
            t30 = 1/t14**5
            t31 = t19**2
            t32 = 1/t11**2
            t33 = rhoa**2
            t34 = 1/t10**3
            t35 = 1/t5**6.666666666666666d+0
            t36 = t25**2
            t37 = 1.8948337999999998d-3*gammaaa/t1**3-1.8948337999999998
     1         d-3*t26
            fnc(iq) = fnc(iq)-2.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t21*t22*t3*t4*t6-1.3
     1         333333333333333d+0*gammaaa*t15*t24*t3*t6+2.66666666666666
     2         66d+0*gammaaa*t15*t18*t23*t3)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*(5.0d-1*gammaaa*t22*t28*t3*t4*t6+t1
     1         5*t29*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(5.0d-1*gammaaa*t22*t3*t4*t6*(1.
     1         8d+1*gammaaa*t12*t20*t6*t7*(1.9555555555555554d+1*t16*t33
     2         *t8/t5**9.333333333333332d+0-2.6666666666666666d+0*t16*t1
     3         7*t8-7.11111111111111d+0*gammaaa**4*t33*t34/t5**1.4666666
     4         666666667d+1+1.2444444444444443d+1*gammaaa*t33*t35-2.6666
     5         666666666666d+0*gammaaa*t18)-1.8d+1*gammaaa*t12*t31*t32*t
     6         6*t7+1.8d+1*gammaaa*t31*t32*t6*t7+1.1199999999999999d+2*g
     7         ammaaa*t13*t33*t35*t7-9.6d+1*gammaaa*rhoa*t12*t18*t19*t20
     8         *t7-2.4d+1*gammaaa*t13*t18*t7)-7.5d-1*gammaaa*t21**2*t3*t
     9         30*t4*t6+1.3333333333333333d+0*gammaaa*t21*t22*t24*t3*t6-
     :         4.444444444444444d-1*gammaaa*t15*t3*t6/rhoa**6.6666666666
     ;         66666d-1+9.777777777777779d+0*gammaaa*t15*t18*t3*t4-1.244
     <         4444444444443d+1*gammaaa*rhoa**3.3333333333333337d+0*t15*
     =         t3*t35-2.6666666666666666d+0*gammaaa*t18*t21*t22*t23*t3)*
     >         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(5.0d-1*t22*t6*(gammaaa*t3*t4*(
     1         1.8d+1*gammaaa*t12*t20*(2.6666666666666666d+0*gammaaa**3*
     2         rhoa*t34/t5**1.2666666666666665d+1-2.6666666666666666d+0*
     3         rhoa*t18-5.333333333333333d+0*gammaaa*rhoa*t16*t17)*t6*t7
     4         -1.8d+1*gammaaa*t12*t19*t25*t32*t6*t7+1.8d+1*gammaaa*t19*
     5         t25*t32*t6*t7+1.8d+1*t12*t19*t20*t6*t7-4.8d+1*gammaaa*rho
     6         a*t12*t18*t20*t25*t7-2.4d+1*rhoa*t13*t18*t7+3.6d+1*gammaa
     7         a*t12*t19*t20*t27*t3*t6-4.8d+1*gammaaa*rhoa*t13*t18*t27*t
     8         3)+t21*t3*t4+gammaaa*t21*t27*t4)-7.5d-1*gammaaa*t21*t28*t
     9         3*t30*t4*t6+3.333333333333333d-1*t15*(-4*t24*t3-4*gammaaa
     :         *t24*t27)*t6+6.666666666666666d-1*gammaaa*t22*t24*t28*t3*
     ;         t6+t15*t18*(2.6666666666666666d+0*t23*t3+2.66666666666666
     <         66d+0*gammaaa*t23*t27)-1.3333333333333333d+0*gammaaa*t18*
     =         t22*t23*t28*t3)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*(5.0d-1*t22*t6*(gammaaa*t3*t4*
     1         (1.8d+1*gammaaa*t12*t20*t6*t7*(t16*t9-t34*t8/t5**1.066666
     2         6666666666d+1)-1.8d+1*gammaaa*t12*t32*t36*t6*t7+1.8d+1*ga
     3         mmaaa*t32*t36*t6*t7+3.6d+1*t12*t20*t25*t6*t7+1.8d+1*gamma
     4         aa*t13*t3*t37*t6+7.2d+1*gammaaa*t12*t20*t25*t27*t3*t6+3.6
     5         d+1*t13*t27*t3*t6+1.8d+1*gammaaa*t13*t27**2*t6)+t28*t3*t4
     6         +gammaaa*t27*t28*t4-t28*t29)+t15*(-gammaaa*t37*t4-2*t27*t
     7         4)*t6-7.5d-1*gammaaa*t28**2*t3*t30*t4*t6)*wght+Cmat2(iq,D
     8         2_GAA_GAA)
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
            t32 = 1/t5**7.333333333333333d+0
            t33 = 1/t5**4.666666666666666d+0
            t34 = -2.6666666666666666d+0*rhoa*t31*t32*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t33
            t35 = 1/t11
            t36 = 1.8d+1*gammaaa*t12*t34*t35*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t33*t7
            t37 = 1/t14**3
            t38 = rhoa**2.3333333333333334d+0
            t39 = rhoa**3.333333333333333d-1
            t40 = 1/t25
            t41 = 1/t20**7.333333333333333d+0
            t42 = 1/t20**4.666666666666666d+0
            t43 = -2.6666666666666666d+0*gammabb*rhob*t42-2.666666666666
     1         6666d+0*rhob*t23*t40*t41
            t44 = 1/t26
            t45 = 1.8d+1*gammabb*t21*t22*t27*t43*t44-2.4d+1*gammabb*rhob
     1         *t22*t28*t42
            t46 = 1/t29**3
            t47 = rhob**2.3333333333333334d+0
            t48 = rhob**3.333333333333333d-1
            t49 = gammaaa*t31*t9+t6
            t50 = 1/t1**2
            t51 = 9.474169d-4*t2-9.474169d-4*gammaaa*t50
            t52 = 1.8d+1*gammaaa*t12*t35*t49*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t3*t51*t6
            t53 = -gammaaa*t4*t51-t3*t4
            t54 = gammabb*t24*t40+t21
            t55 = 1/t16**2
            t56 = 9.474169d-4*t17-9.474169d-4*gammabb*t55
            t57 = 1.8d+1*gammabb*t18*t21*t28*t56+1.8d+1*gammabb*t21*t22*
     1         t27*t44*t54+9.0d+0*t21*t22*t28
            t58 = -gammabb*t19*t56-t18*t19
            t59 = 1/t14**5
            t60 = t34**2
            t61 = 1/t11**2
            t62 = rhoa**2
            t63 = 1/t10**3
            t64 = 1/t5**6.666666666666666d+0
            t65 = 1/t29**5
            t66 = t43**2
            t67 = 1/t26**2
            t68 = rhob**2
            t69 = 1/t25**3
            t70 = 1/t20**6.666666666666666d+0
            t71 = t49**2
            t72 = 1.8948337999999998d-3*gammaaa/t1**3-1.8948337999999998
     1         d-3*t50
            t73 = t54**2
            t74 = 1.8948337999999998d-3*gammabb/t16**3-1.894833799999999
     1         8d-3*t55
            fnc(iq) = 1.0d+0*(-gammaaa*t15*t3*t4*t6-gammabb*t18*t19*t21*
     1         t30)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t3*t36*t37*t4*t6-1.3
     1         333333333333333d+0*gammaaa*t15*t3*t39*t6+2.66666666666666
     2         66d+0*gammaaa*t15*t3*t33*t38)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-1.3333333333333333d+0*gammabb*t18*
     1         t21*t30*t48+2.6666666666666666d+0*gammabb*t18*t30*t42*t47
     2         +5.0d-1*gammabb*t18*t19*t21*t45*t46)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 1.0d+0*(t15*t53*t6+5.0d-1*gammaaa*t3*t37*t
     1         4*t52*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*(t21*t30*t58+5.0d-1*gammabb*t18*t19
     1         *t21*t46*t57)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(5.0d-1*gammaaa*t3*t37*t4*t6*(1.
     1         8d+1*gammaaa*t12*t35*t6*t7*(1.9555555555555554d+1*t31*t62
     2         *t8/t5**9.333333333333332d+0-2.6666666666666666d+0*t31*t3
     3         2*t8+1.2444444444444443d+1*gammaaa*t62*t64-7.111111111111
     4         11d+0*gammaaa**4*t62*t63/t5**1.4666666666666667d+1-2.6666
     5         666666666666d+0*gammaaa*t33)+1.1199999999999999d+2*gammaa
     6         a*t13*t62*t64*t7-1.8d+1*gammaaa*t12*t6*t60*t61*t7+1.8d+1*
     7         gammaaa*t6*t60*t61*t7-9.6d+1*gammaaa*rhoa*t12*t33*t34*t35
     8         *t7-2.4d+1*gammaaa*t13*t33*t7)-1.2444444444444443d+1*gamm
     9         aaa*rhoa**3.3333333333333337d+0*t15*t3*t64-7.5d-1*gammaaa
     :         *t3*t36**2*t4*t59*t6+1.3333333333333333d+0*gammaaa*t3*t36
     ;         *t37*t39*t6-4.444444444444444d-1*gammaaa*t15*t3*t6/rhoa**
     <         6.666666666666666d-1+9.777777777777779d+0*gammaaa*t15*t3*
     =         t33*t4-2.6666666666666666d+0*gammaaa*t3*t33*t36*t37*t38)*
     >         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(5.0d-1*gammabb*t18*t19*t21*t46*
     1         (1.8d+1*gammabb*t21*t22*t27*t44*(1.2444444444444443d+1*ga
     2         mmabb*t68*t70-7.11111111111111d+0*gammabb**4*t68*t69/t20*
     3         *1.4666666666666667d+1+1.9555555555555554d+1*t23*t40*t68/
     4         t20**9.333333333333332d+0-2.6666666666666666d+0*gammabb*t
     5         42-2.6666666666666666d+0*t23*t40*t41)+1.1199999999999999d
     6         +2*gammabb*t22*t28*t68*t70-1.8d+1*gammabb*t21*t22*t27*t66
     7         *t67+1.8d+1*gammabb*t21*t22*t66*t67-9.6d+1*gammabb*rhob*t
     8         22*t27*t42*t43*t44-2.4d+1*gammabb*t22*t28*t42)-1.24444444
     9         44444443d+1*gammabb*rhob**3.3333333333333337d+0*t18*t30*t
     :         70-7.5d-1*gammabb*t18*t19*t21*t45**2*t65+1.33333333333333
     ;         33d+0*gammabb*t18*t21*t45*t46*t48-2.6666666666666666d+0*g
     <         ammabb*t18*t42*t45*t46*t47+9.777777777777779d+0*gammabb*t
     =         18*t19*t30*t42-4.444444444444444d-1*gammabb*t18*t21*t30/r
     >         hob**6.666666666666666d-1)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(5.0d-1*t37*t6*(gammaaa*t3*t4*(
     1         1.8d+1*gammaaa*t12*t35*t6*(2.6666666666666666d+0*gammaaa*
     2         *3*rhoa*t63/t5**1.2666666666666665d+1-2.6666666666666666d
     3         +0*rhoa*t33-5.333333333333333d+0*gammaaa*rhoa*t31*t32)*t7
     4         -1.8d+1*gammaaa*t12*t34*t49*t6*t61*t7+1.8d+1*gammaaa*t34*
     5         t49*t6*t61*t7+1.8d+1*t12*t34*t35*t6*t7-4.8d+1*gammaaa*rho
     6         a*t12*t33*t35*t49*t7-2.4d+1*rhoa*t13*t33*t7+3.6d+1*gammaa
     7         a*t12*t3*t34*t35*t51*t6-4.8d+1*gammaaa*rhoa*t13*t3*t33*t5
     8         1)+gammaaa*t36*t4*t51+t3*t36*t4)-7.5d-1*gammaaa*t3*t36*t4
     9         *t52*t59*t6+6.666666666666666d-1*gammaaa*t3*t37*t39*t52*t
     :         6+3.333333333333333d-1*t15*(-4*gammaaa*t39*t51-4*t3*t39)*
     ;         t6-1.3333333333333333d+0*gammaaa*t3*t33*t37*t38*t52+t15*t
     <         33*(2.6666666666666666d+0*gammaaa*t38*t51+2.6666666666666
     =         666d+0*t3*t38))*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(5.0d-1*t21*t46*(gammabb*t18*t1
     1         9*(1.8d+1*gammabb*t21*t22*t27*t44*(2.6666666666666666d+0*
     2         gammabb**3*rhob*t69/t20**1.2666666666666665d+1-2.66666666
     3         66666666d+0*rhob*t42-5.333333333333333d+0*gammabb*rhob*t4
     4         0*t41)-1.8d+1*gammabb*t21*t22*t27*t43*t54*t67+1.8d+1*gamm
     5         abb*t21*t22*t43*t54*t67+3.6d+1*gammabb*t18*t21*t27*t43*t4
     6         4*t56-4.8d+1*gammabb*rhob*t18*t28*t42*t56-4.8d+1*gammabb*
     7         rhob*t22*t27*t42*t44*t54+1.8d+1*t21*t22*t27*t43*t44-2.4d+
     8         1*rhob*t22*t28*t42)+gammabb*t19*t45*t56+t18*t19*t45)-7.5d
     9         -1*gammabb*t18*t19*t21*t45*t57*t65+6.666666666666666d-1*g
     :         ammabb*t18*t21*t46*t48*t57-1.3333333333333333d+0*gammabb*
     ;         t18*t42*t46*t47*t57+3.333333333333333d-1*t21*t30*(-4*gamm
     <         abb*t48*t56-4*t18*t48)+t30*t42*(2.6666666666666666d+0*gam
     =         mabb*t47*t56+2.6666666666666666d+0*t18*t47))*wght+Cmat2(i
     >         q,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*(5.0d-1*t37*t6*(gammaaa*t3*t4*
     1         (1.8d+1*gammaaa*t12*t35*t6*t7*(t31*t9-t63*t8/t5**1.066666
     2         6666666666d+1)+1.8d+1*gammaaa*t13*t3*t6*t72-1.8d+1*gammaa
     3         a*t12*t6*t61*t7*t71+1.8d+1*gammaaa*t6*t61*t7*t71+3.6d+1*t
     4         12*t35*t49*t6*t7+1.8d+1*gammaaa*t13*t51**2*t6+7.2d+1*gamm
     5         aaa*t12*t3*t35*t49*t51*t6+3.6d+1*t13*t3*t51*t6)-t52*t53+g
     6         ammaaa*t4*t51*t52+t3*t4*t52)+t15*t6*(-gammaaa*t4*t72-2*t4
     7         *t51)-7.5d-1*gammaaa*t3*t4*t52**2*t59*t6)*wght+Cmat2(iq,D
     8         2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*(5.0d-1*t21*t46*(gammabb*t18*t
     1         19*(1.8d+1*gammabb*t18*t21*t28*t74-1.8d+1*gammabb*t21*t22
     2         *t27*t67*t73+1.8d+1*gammabb*t21*t22*t67*t73+1.8d+1*gammab
     3         b*t21*t22*t27*t44*(t24*t40-t23*t69/t20**1.066666666666666
     4         6d+1)+1.8d+1*gammabb*t21*t28*t56**2+7.2d+1*gammabb*t18*t2
     5         1*t27*t44*t54*t56+3.6d+1*t18*t21*t28*t56+3.6d+1*t21*t22*t
     6         27*t44*t54)-t57*t58+gammabb*t19*t56*t57+t18*t19*t57)+t21*
     7         t30*(-gammabb*t19*t74-2*t19*t56)-7.5d-1*gammabb*t18*t19*t
     8         21*t57**2*t65)*wght+Cmat2(iq,D2_GBB_GBB)
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhoa*t16*t17*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammaaa*t12*t19*t20*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhoa**2.3333333333333334d+0
            t24 = rhoa**3.333333333333333d-1
            t25 = gammaaa*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammaaa*t26
            t28 = 1.8d+1*gammaaa*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t27*t3*t6
            t29 = -1.0d+0*t3*t4*wght-1.0d+0*gammaaa*t27*t4*wght
            t30 = 1/t14**5
            t31 = 1/t5**6.666666666666666d+0
            t32 = t19**2
            t33 = 1/t11**2
            t34 = rhoa**2
            t35 = 1/t10**3
            t36 = 1.8948337999999998d-3*gammaaa/t1**3-1.8948337999999998
     1         d-3*t26
            t37 = t25**2
            fnc(iq) = fnc(iq)-1.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 5.0d-1*gammaaa*t21*t22*t3*t4*t6*wght-1.3333
     1         333333333333d+0*gammaaa*t15*t24*t3*t6*wght+2.666666666666
     2         6666d+0*gammaaa*t15*t18*t23*t3*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 5.0d-1*gammaaa*t22*t28*t3*t4*t6*wght+t15*t
     1         29*t6+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d+1*ga
     1         mmaaa*t12*t20*t6*t7*(1.9555555555555554d+1*t16*t34*t8/t5*
     2         *9.333333333333332d+0-2.6666666666666666d+0*t16*t17*t8-7.
     3         11111111111111d+0*gammaaa**4*t34*t35/t5**1.46666666666666
     4         67d+1+1.2444444444444443d+1*gammaaa*t31*t34-2.66666666666
     5         66666d+0*gammaaa*t18)-1.8d+1*gammaaa*t12*t32*t33*t6*t7+1.
     6         8d+1*gammaaa*t32*t33*t6*t7+1.1199999999999999d+2*gammaaa*
     7         t13*t31*t34*t7-9.6d+1*gammaaa*rhoa*t12*t18*t19*t20*t7-2.4
     8         d+1*gammaaa*t13*t18*t7)*wght+1.3333333333333333d+0*gammaa
     9         a*t21*t24*t3*wght)-7.5d-1*gammaaa*t21**2*t3*t30*t4*t6*wgh
     :         t-4.444444444444444d-1*gammaaa*t15*t3*t6*wght/rhoa**6.666
     ;         666666666666d-1+9.777777777777779d+0*gammaaa*t15*t18*t3*t
     <         4*wght-1.2444444444444443d+1*gammaaa*rhoa**3.333333333333
     =         3337d+0*t15*t3*t31*wght-2.6666666666666666d+0*gammaaa*t18
     >         *t21*t22*t23*t3*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d+1*g
     1         ammaaa*t12*t20*(2.6666666666666666d+0*gammaaa**3*rhoa*t35
     2         /t5**1.2666666666666665d+1-2.6666666666666666d+0*rhoa*t18
     3         -5.333333333333333d+0*gammaaa*rhoa*t16*t17)*t6*t7-1.8d+1*
     4         gammaaa*t12*t19*t25*t33*t6*t7+1.8d+1*gammaaa*t19*t25*t33*
     5         t6*t7+1.8d+1*t12*t19*t20*t6*t7-4.8d+1*gammaaa*rhoa*t12*t1
     6         8*t20*t25*t7-2.4d+1*rhoa*t13*t18*t7+3.6d+1*gammaaa*t12*t1
     7         9*t20*t27*t3*t6-4.8d+1*gammaaa*rhoa*t13*t18*t27*t3)*wght+
     8         5.0d-1*t21*t3*t4*wght+5.0d-1*gammaaa*t21*t27*t4*wght+6.66
     9         6666666666666d-1*gammaaa*t24*t28*t3*wght)+t15*t6*(-1.3333
     :         333333333333d+0*t24*t3*wght-1.3333333333333333d+0*gammaaa
     ;         *t24*t27*wght)+t15*t18*(2.6666666666666666d+0*t23*t3*wght
     <         +2.6666666666666666d+0*gammaaa*t23*t27*wght)-7.5d-1*gamma
     =         aa*t21*t28*t3*t30*t4*t6*wght-1.3333333333333333d+0*gammaa
     >         a*t18*t22*t23*t28*t3*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d+1*
     1         gammaaa*t12*t20*t6*t7*(t16*t9-t35*t8/t5**1.06666666666666
     2         66d+1)-1.8d+1*gammaaa*t12*t33*t37*t6*t7+1.8d+1*gammaaa*t3
     3         3*t37*t6*t7+3.6d+1*t12*t20*t25*t6*t7+1.8d+1*gammaaa*t13*t
     4         3*t36*t6+7.2d+1*gammaaa*t12*t20*t25*t27*t3*t6+3.6d+1*t13*
     5         t27*t3*t6+1.8d+1*gammaaa*t13*t27**2*t6)*wght+5.0d-1*t28*t
     6         3*t4*wght+5.0d-1*gammaaa*t27*t28*t4*wght)+t15*t6*(-1.0d+0
     7         *gammaaa*t36*t4*wght-2.0d+0*t27*t4*wght)-7.5d-1*gammaaa*t
     8         28**2*t3*t30*t4*t6*wght-5.0d-1*t22*t28*t29*t6+Cmat2(iq,D2
     9         _GAA_GAA)
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhob*t16*t17*t8-2.6666666666666
     1         666d+0*gammabb*rhob*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammabb*t12*t19*t20*t6*t7-2.4d+1*gammabb*rhob*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhob**2.3333333333333334d+0
            t24 = rhob**3.333333333333333d-1
            t25 = gammabb*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammabb*t26
            t28 = 1.8d+1*gammabb*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammabb*t13*t27*t3*t6
            t29 = -1.0d+0*t3*t4*wght-1.0d+0*gammabb*t27*t4*wght
            t30 = 1/t14**5
            t31 = 1/t5**6.666666666666666d+0
            t32 = t19**2
            t33 = 1/t11**2
            t34 = rhob**2
            t35 = 1/t10**3
            t36 = 1.8948337999999998d-3*gammabb/t1**3-1.8948337999999998
     1         d-3*t26
            t37 = t25**2
            fnc(iq) = fnc(iq)-1.0d+0*gammabb*t15*t3*t4*t6*wght
            Amat(iq,D1_RB) = 5.0d-1*gammabb*t21*t22*t3*t4*t6*wght-1.3333
     1         333333333333d+0*gammabb*t15*t24*t3*t6*wght+2.666666666666
     2         6666d+0*gammabb*t15*t18*t23*t3*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 5.0d-1*gammabb*t22*t28*t3*t4*t6*wght+t15*t
     1         29*t6+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d+1*ga
     1         mmabb*t12*t20*t6*t7*(1.9555555555555554d+1*t16*t34*t8/t5*
     2         *9.333333333333332d+0-2.6666666666666666d+0*t16*t17*t8-7.
     3         11111111111111d+0*gammabb**4*t34*t35/t5**1.46666666666666
     4         67d+1+1.2444444444444443d+1*gammabb*t31*t34-2.66666666666
     5         66666d+0*gammabb*t18)-1.8d+1*gammabb*t12*t32*t33*t6*t7+1.
     6         8d+1*gammabb*t32*t33*t6*t7+1.1199999999999999d+2*gammabb*
     7         t13*t31*t34*t7-9.6d+1*gammabb*rhob*t12*t18*t19*t20*t7-2.4
     8         d+1*gammabb*t13*t18*t7)*wght+1.3333333333333333d+0*gammab
     9         b*t21*t24*t3*wght)-7.5d-1*gammabb*t21**2*t3*t30*t4*t6*wgh
     :         t-4.444444444444444d-1*gammabb*t15*t3*t6*wght/rhob**6.666
     ;         666666666666d-1+9.777777777777779d+0*gammabb*t15*t18*t3*t
     <         4*wght-1.2444444444444443d+1*gammabb*rhob**3.333333333333
     =         3337d+0*t15*t3*t31*wght-2.6666666666666666d+0*gammabb*t18
     >         *t21*t22*t23*t3*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d+1*g
     1         ammabb*t12*t20*(2.6666666666666666d+0*gammabb**3*rhob*t35
     2         /t5**1.2666666666666665d+1-2.6666666666666666d+0*rhob*t18
     3         -5.333333333333333d+0*gammabb*rhob*t16*t17)*t6*t7-1.8d+1*
     4         gammabb*t12*t19*t25*t33*t6*t7+1.8d+1*gammabb*t19*t25*t33*
     5         t6*t7+1.8d+1*t12*t19*t20*t6*t7-4.8d+1*gammabb*rhob*t12*t1
     6         8*t20*t25*t7-2.4d+1*rhob*t13*t18*t7+3.6d+1*gammabb*t12*t1
     7         9*t20*t27*t3*t6-4.8d+1*gammabb*rhob*t13*t18*t27*t3)*wght+
     8         5.0d-1*t21*t3*t4*wght+5.0d-1*gammabb*t21*t27*t4*wght+6.66
     9         6666666666666d-1*gammabb*t24*t28*t3*wght)+t15*t6*(-1.3333
     :         333333333333d+0*t24*t3*wght-1.3333333333333333d+0*gammabb
     ;         *t24*t27*wght)+t15*t18*(2.6666666666666666d+0*t23*t3*wght
     <         +2.6666666666666666d+0*gammabb*t23*t27*wght)-7.5d-1*gamma
     =         bb*t21*t28*t3*t30*t4*t6*wght-1.3333333333333333d+0*gammab
     >         b*t18*t22*t23*t28*t3*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d+1*
     1         gammabb*t12*t20*t6*t7*(t16*t9-t35*t8/t5**1.06666666666666
     2         66d+1)-1.8d+1*gammabb*t12*t33*t37*t6*t7+1.8d+1*gammabb*t3
     3         3*t37*t6*t7+3.6d+1*t12*t20*t25*t6*t7+1.8d+1*gammabb*t13*t
     4         3*t36*t6+7.2d+1*gammabb*t12*t20*t25*t27*t3*t6+3.6d+1*t13*
     5         t27*t3*t6+1.8d+1*gammabb*t13*t27**2*t6)*wght+5.0d-1*t28*t
     6         3*t4*wght+5.0d-1*gammabb*t27*t28*t4*wght)+t15*t6*(-1.0d+0
     7         *gammabb*t36*t4*wght-2.0d+0*t27*t4*wght)-7.5d-1*gammabb*t
     8         28**2*t3*t30*t4*t6*wght-5.0d-1*t22*t28*t29*t6+Cmat2(iq,D2
     9         _GBB_GBB)
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
C>   {\it t_3} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\alpha\alpha}}\over{\sigma_{\alpha\alpha}
C>    +6255746.320201}}+0.002913644\\\\ 
C>   {\it t_4} &=& \left| \rho_\beta\right| \\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_4}^{2.666666666666667}}}\\\\ 
C>   {\it t_6} &=& {{9.474169 \times 10^{-4}\,
C>    \sigma_{\beta\beta}}\over{\sigma_{\beta\beta}+6255746.320201}}
C>    +0.002913644\\\\ 
C>   {\it t_7} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_8} &=& {{1}\over{{\it t_7}^{2.666666666666667}}}\\\\ 
C>   {\it t_9} &=& {{9.474169 \times 10^{-4}\,\sigma_{ss}}
C>    \over{\sigma_{ss}+6255746.320201}}+0.002913644\\\\ 
C>   f &=& 1.0\,\left(-{{\rho_\beta^{{{4}\over{3}}}\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}}\over{\sqrt{9.0\,{\it t_5}
C>    \,\sigma_{\beta\beta}\,{\it t_6}^{2.0}\,\log ^{2.0}
C>    \left(\sqrt{{{\sigma_{\beta\beta}^2}\over{{
C>    \it t_4}^{5.333333333333333}}}+1.0}+{\it t_5}\,
C>    \sigma_{\beta\beta}\right)+1.0}}}-{{\rho_\alpha^{{{4}
C>    \over{3}}}\,{\it t_2}\,\sigma_{\alpha\alpha}\,{\it t_3}}
C>    \over{\sqrt{9.0\,{\it t_2}\,\sigma_{\alpha\alpha}\,{
C>    \it t_3}^{2.0}\,\log ^{2.0}\left(\sqrt{{{
C>    \sigma_{\alpha\alpha}^2}\over{{\it t_1}^{5.333333333333333}}}
C>    +1.0}+{\it t_2}\,\sigma_{\alpha\alpha}\right)+1.0}}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}}\over{\sqrt{9.0\,{\it t_8}\,\sigma_{ss}\,{
C>    \it t_9}^{2.0}\,\log ^{2.0}\left(\sqrt{{{\sigma_{ss}^2}
C>    \over{{\it t_7}^{5.333333333333333}}}+1.0}+{\it t_8}\,
C>    \sigma_{ss}\right)+1.0}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Filatov, W Thiel, Mol.Phys. 91, 847 (1997)  , DOI: 
C> <a href="https://doi.org/10.1080/002689797170950 "> 
C> 10.1080/002689797170950 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_ft97_d3(param,tol_rho,ipol,nq,wght, 
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhoa*t16*t17*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammaaa*t12*t19*t20*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhoa**2.3333333333333334d+0
            t24 = rhoa**3.333333333333333d-1
            t25 = gammaaa*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammaaa*t26
            t28 = 1.8d+1*gammaaa*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t27*t3*t6
            t29 = -t3*t4-gammaaa*t27*t4
            t30 = t21**2
            t31 = 1/t14**5
            t32 = t19**2
            t33 = 1/t11**2
            t34 = gammaaa**4
            t35 = rhoa**2
            t36 = 1/t10**3
            t37 = 1/t5**1.4666666666666667d+1
            t38 = 1/t5**9.333333333333332d+0
            t39 = 1/t5**6.666666666666666d+0
            t40 = 1.9555555555555554d+1*t16*t35*t38*t8-2.666666666666666
     1         6d+0*t16*t17*t8+1.2444444444444443d+1*gammaaa*t35*t39-7.1
     2         1111111111111d+0*t34*t35*t36*t37-2.6666666666666666d+0*ga
     3         mmaaa*t18
            t41 = 1.8d+1*gammaaa*t12*t20*t40*t6*t7-1.8d+1*gammaaa*t12*t3
     1         2*t33*t6*t7+1.8d+1*gammaaa*t32*t33*t6*t7+1.11999999999999
     2         99d+2*gammaaa*t13*t35*t39*t7-9.6d+1*gammaaa*rhoa*t12*t18*
     3         t19*t20*t7-2.4d+1*gammaaa*t13*t18*t7
            t42 = rhoa**3.3333333333333337d+0
            t43 = 1/rhoa**6.666666666666666d-1
            t44 = gammaaa**3
            t45 = 1/t5**1.2666666666666665d+1
            t46 = 2.6666666666666666d+0*rhoa*t36*t44*t45-2.6666666666666
     1         666d+0*rhoa*t18-5.333333333333333d+0*gammaaa*rhoa*t16*t17
            t47 = 1.8d+1*gammaaa*t12*t20*t46*t6*t7-1.8d+1*gammaaa*t12*t1
     1         9*t25*t33*t6*t7+1.8d+1*gammaaa*t19*t25*t33*t6*t7+1.8d+1*t
     2         12*t19*t20*t6*t7-4.8d+1*gammaaa*rhoa*t12*t18*t20*t25*t7-2
     3         .4d+1*rhoa*t13*t18*t7+3.6d+1*gammaaa*t12*t19*t20*t27*t3*t
     4         6-4.8d+1*gammaaa*rhoa*t13*t18*t27*t3
            t48 = gammaaa*t3*t4*t47+t21*t3*t4+gammaaa*t21*t27*t4
            t49 = 2.6666666666666666d+0*t23*t3+2.6666666666666666d+0*gam
     1         maaa*t23*t27
            t50 = -4*t24*t3-4*gammaaa*t24*t27
            t51 = t25**2
            t52 = 1/t5**1.0666666666666666d+1
            t53 = t16*t9-t36*t52*t8
            t54 = 1/t1**3
            t55 = 1.8948337999999998d-3*gammaaa*t54-1.8948337999999998d-
     1         3*t26
            t56 = t27**2
            t57 = 1.8d+1*gammaaa*t12*t20*t53*t6*t7-1.8d+1*gammaaa*t12*t3
     1         3*t51*t6*t7+1.8d+1*gammaaa*t33*t51*t6*t7+3.6d+1*t12*t20*t
     2         25*t6*t7+1.8d+1*gammaaa*t13*t56*t6+1.8d+1*gammaaa*t13*t3*
     3         t55*t6+7.2d+1*gammaaa*t12*t20*t25*t27*t3*t6+3.6d+1*t13*t2
     4         7*t3*t6
            t58 = gammaaa*t3*t4*t57+t28*t3*t4+gammaaa*t27*t28*t4-t28*t29
            t59 = t28**2
            t60 = -gammaaa*t4*t55-2*t27*t4
            t61 = 1/t14**7
            t62 = t19**3
            t63 = 1/t11**3
            t64 = rhoa**3
            t65 = 1/t10**5
            t66 = 1/t5**8.666666666666666d+0
            t67 = t25**3
            t68 = 5.684501399999999d-3*t54-5.684501399999999d-3*gammaaa/
     1         t1**4
            fnc(iq) = fnc(iq)-2.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t21*t22*t3*t4*t6-1.3
     1         333333333333333d+0*gammaaa*t15*t24*t3*t6+2.66666666666666
     2         66d+0*gammaaa*t15*t18*t23*t3)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*(5.0d-1*gammaaa*t22*t28*t3*t4*t6+t1
     1         5*t29*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-4.444444444444444d-1*gammaaa*t
     1         15*t3*t43*t6+5.0d-1*gammaaa*t22*t3*t4*t41*t6-7.5d-1*gamma
     2         aa*t3*t30*t31*t4*t6+1.3333333333333333d+0*gammaaa*t21*t22
     3         *t24*t3*t6-1.2444444444444443d+1*gammaaa*t15*t3*t39*t42+9
     4         .777777777777779d+0*gammaaa*t15*t18*t3*t4-2.6666666666666
     5         666d+0*gammaaa*t18*t21*t22*t23*t3)*wght+Amat2(iq,D2_RA_RA
     6         )
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(3.333333333333333d-1*t15*t50*t
     1         6+5.0d-1*t22*t48*t6-7.5d-1*gammaaa*t21*t28*t3*t31*t4*t6+6
     2         .666666666666666d-1*gammaaa*t22*t24*t28*t3*t6+t15*t18*t49
     3         -1.3333333333333333d+0*gammaaa*t18*t22*t23*t28*t3)*wght+C
     4         mat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*(t15*t6*t60-7.5d-1*gammaaa*t3*
     1         t31*t4*t59*t6+5.0d-1*t22*t58*t6)*wght+Cmat2(iq,D2_GAA_GAA
     2         )
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(5.0d-1*gammaaa*t22*t3*t4*t6*
     1         (1.8d+1*gammaaa*t12*t20*t6*t7*(-1.8251851851851847d+2*t16
     2         *t64*t8/t5**1.1333333333333333d+1+5.866666666666666d+1*rh
     3         oa*t16*t38*t8-8.296296296296295d+1*gammaaa*t64*t66-5.6888
     4         88888888889d+1*gammaaa**6*t64*t65/t5**2.2d+1+1.5644444444
     5         444444d+2*t34*t36*t64/t5**1.6666666666666663d+1+3.7333333
     6         33333333d+1*gammaaa*rhoa*t39-2.1333333333333332d+1*rhoa*t
     7         34*t36*t37)-7.466666666666665d+2*gammaaa*t13*t64*t66*t7+3
     8         .6d+1*gammaaa*t12*t6*t62*t63*t7-5.4d+1*gammaaa*t6*t62*t63
     9         *t7-5.4d+1*gammaaa*t12*t19*t33*t40*t6*t7+5.4d+1*gammaaa*t
     :         19*t33*t40*t6*t7-1.44d+2*gammaaa*rhoa*t12*t18*t20*t40*t7+
     ;         6.719999999999999d+2*gammaaa*t12*t19*t20*t35*t39*t7+3.359
     <         9999999999997d+2*gammaaa*rhoa*t13*t39*t7+1.44d+2*gammaaa*
     =         rhoa*t12*t18*t32*t33*t7-1.44d+2*gammaaa*rhoa*t18*t32*t33*
     >         t7-1.44d+2*gammaaa*t12*t18*t19*t20*t7)+8.296296296296295d
     ?         +1*gammaaa*rhoa**4.333333333333333d+0*t15*t3*t66+1.875d+0
     @         *gammaaa*t21**3*t3*t4*t6*t61+6.666666666666666d-1*gammaaa
     1         *t21*t22*t3*t43*t6-2.25d+0*gammaaa*t21*t3*t31*t4*t41*t6+2
     2         *gammaaa*t22*t24*t3*t41*t6-3*gammaaa*t24*t3*t30*t31*t6+2.
     3         962962962962963d-1*gammaaa*t15*t3*t6/rhoa**1.666666666666
     4         6669d+0+1.8666666666666665d+1*gammaaa*t21*t22*t3*t39*t42+
     5         t18*t22*(-4.0d+0*gammaaa*t23*t3*t41-1.4666666666666667d+1
     6         *gammaaa*t21*t3*t4)-8.711111111111112d+1*gammaaa*t15*t23*
     7         t3*t39+6.0d+0*gammaaa*t18*t23*t3*t30*t31+1.42222222222222
     8         22d+1*gammaaa*t15*t18*t24*t3)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(5.0d-1*t22*t6*(gammaaa*t3*t
     1         4*(1.8d+1*gammaaa*t12*t20*t6*(2.1333333333333332d+1*gamma
     2         aa**5*t35*t65/t5**2.0d+1-1.9555555555555554d+1*t35*t36*t4
     3         4/t5**1.4666666666666664d+1+2.6666666666666666d+0*t36*t44
     4         *t45-2.8444444444444444d+1*t35*t36*t37*t44+1.244444444444
     5         4443d+1*t35*t39+3.911111111111111d+1*gammaaa*t16*t35*t38-
     6         2.6666666666666666d+0*t18-5.333333333333333d+0*gammaaa*t1
     7         6*t17)*t7+3.6d+1*gammaaa*t12*t25*t32*t6*t63*t7-5.4d+1*gam
     8         maaa*t25*t32*t6*t63*t7-3.6d+1*gammaaa*t12*t19*t33*t46*t6*
     9         t7+3.6d+1*gammaaa*t19*t33*t46*t6*t7-1.8d+1*gammaaa*t12*t2
     :         5*t33*t40*t6*t7+1.8d+1*gammaaa*t25*t33*t40*t6*t7+1.8d+1*t
     ;         12*t20*t40*t6*t7-1.8d+1*t12*t32*t33*t6*t7+1.8d+1*t32*t33*
     <         t6*t7-9.6d+1*gammaaa*rhoa*t12*t18*t20*t46*t7+2.2399999999
     =         999998d+2*gammaaa*t12*t20*t25*t35*t39*t7+1.11999999999999
     >         99d+2*t13*t35*t39*t7+9.6d+1*gammaaa*rhoa*t12*t18*t19*t25*
     ?         t33*t7-9.6d+1*gammaaa*rhoa*t18*t19*t25*t33*t7-4.8d+1*gamm
     @         aaa*t12*t18*t20*t25*t7-9.6d+1*rhoa*t12*t18*t19*t20*t7-2.4
     1         d+1*t13*t18*t7+3.6d+1*gammaaa*t12*t20*t27*t3*t40*t6-3.6d+
     2         1*gammaaa*t12*t27*t3*t32*t33*t6+3.6d+1*gammaaa*t27*t3*t32
     3         *t33*t6+2.2399999999999998d+2*gammaaa*t13*t27*t3*t35*t39-
     4         1.92d+2*gammaaa*rhoa*t12*t18*t19*t20*t27*t3-4.8d+1*gammaa
     5         a*t13*t18*t27*t3)+t3*t4*t41+gammaaa*t27*t4*t41)+1.875d+0*
     6         gammaaa*t28*t3*t30*t4*t6*t61+3.333333333333333d-1*t22*(4*
     7         gammaaa*t24*t3*t47+4*t21*t24*t3+4*gammaaa*t21*t24*t27)*t6
     8         -1.5d+0*gammaaa*t21*t3*t31*t4*t47*t6+2.222222222222222d-1
     9         *gammaaa*t22*t28*t3*t43*t6+1.111111111111111d-1*t15*(-4*t
     :         3-4*gammaaa*t27)*t43*t6+2.5d-1*t31*(-3*gammaaa*t28*t3*t4*
     ;         t41-3*t3*t30*t4-3*gammaaa*t27*t30*t4)*t6-2*gammaaa*t21*t2
     <         4*t28*t3*t31*t6+t18*t22*(-2.6666666666666666d+0*gammaaa*t
     =         23*t3*t47-4.8888888888888893d+0*gammaaa*t28*t3*t4-2.66666
     >         66666666666d+0*t21*t23*t3-2.6666666666666666d+0*gammaaa*t
     ?         21*t23*t27)+t15*t39*(-1.2444444444444443d+1*t3*t42-1.2444
     @         444444444443d+1*gammaaa*t27*t42)+6.222222222222221d+0*gam
     1         maaa*t22*t28*t3*t39*t42+t15*t18*(9.777777777777779d+0*t3*
     2         t4+9.777777777777779d+0*gammaaa*t27*t4)+4.0d+0*gammaaa*t1
     3         8*t21*t23*t28*t3*t31)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(5.0d-1*t22*t6*(gammaaa*t3*
     1         t4*(1.8d+1*gammaaa*t12*t20*t6*t7*(1.3333333333333333d+1*r
     2         hoa*t36*t45*t8-8.0d+0*rhoa*t34*t65/t5**1.8d+1-5.333333333
     3         333333d+0*rhoa*t16*t17)+3.6d+1*gammaaa*t12*t19*t51*t6*t63
     4         *t7-5.4d+1*gammaaa*t19*t51*t6*t63*t7-1.8d+1*gammaaa*t12*t
     5         19*t33*t53*t6*t7+1.8d+1*gammaaa*t19*t33*t53*t6*t7-3.6d+1*
     6         gammaaa*t12*t25*t33*t46*t6*t7+3.6d+1*gammaaa*t25*t33*t46*
     7         t6*t7+3.6d+1*t12*t20*t46*t6*t7-3.6d+1*t12*t19*t25*t33*t6*
     8         t7+3.6d+1*t19*t25*t33*t6*t7-4.8d+1*gammaaa*rhoa*t12*t18*t
     9         20*t53*t7+4.8d+1*gammaaa*rhoa*t12*t18*t33*t51*t7-4.8d+1*g
     :         ammaaa*rhoa*t18*t33*t51*t7-9.6d+1*rhoa*t12*t18*t20*t25*t7
     ;         +3.6d+1*gammaaa*t12*t19*t20*t56*t6+3.6d+1*gammaaa*t12*t19
     <         *t20*t3*t55*t6+7.2d+1*gammaaa*t12*t20*t27*t3*t46*t6-7.2d+
     =         1*gammaaa*t12*t19*t25*t27*t3*t33*t6+7.2d+1*gammaaa*t19*t2
     >         5*t27*t3*t33*t6+7.2d+1*t12*t19*t20*t27*t3*t6-4.8d+1*gamma
     ?         aa*rhoa*t13*t18*t56-4.8d+1*gammaaa*rhoa*t13*t18*t3*t55-1.
     @         92d+2*gammaaa*rhoa*t12*t18*t20*t25*t27*t3-9.6d+1*rhoa*t13
     1         *t18*t27*t3)+gammaaa*t21*t4*t55+2*t3*t4*t47+2*gammaaa*t27
     2         *t4*t47+2*t21*t27*t4)+1.875d+0*gammaaa*t21*t3*t4*t59*t6*t
     3         61-gammaaa*t24*t3*t31*t59*t6+2.5d-1*t31*(-3*gammaaa*t21*t
     4         3*t4*t57-3*t28*t48-3*gammaaa*t28*t3*t4*t47-3*t21*t28*t3*t
     5         4-3*gammaaa*t21*t27*t28*t4)*t6+3.333333333333333d-1*t22*(
     6         2*gammaaa*t24*t3*t57+2*t24*t28*t3+2*gammaaa*t24*t27*t28)*
     7         t6+3.333333333333333d-1*t15*(-4*gammaaa*t24*t55-8*t24*t27
     8         )*t6-1.6666666666666666d-1*t22*t28*t50*t6+2.0d+0*gammaaa*
     9         t18*t23*t3*t31*t59+t18*t22*(-1.3333333333333333d+0*gammaa
     :         a*t23*t3*t57-1.3333333333333333d+0*t23*t28*t3-1.333333333
     ;         3333333d+0*gammaaa*t23*t27*t28)+t15*t18*(2.66666666666666
     <         66d+0*gammaaa*t23*t55+5.333333333333333d+0*t23*t27)-5.0d-
     =         1*t18*t22*t28*t49)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*(5.0d-1*t22*t6*(gammaaa*t3
     1         *t4*(3.6d+1*gammaaa*t12*t6*t63*t67*t7-5.4d+1*gammaaa*t6*t
     2         63*t67*t7+1.8d+1*gammaaa*t12*t20*t6*(3*t44*t65/t5**1.6d+1
     3         -3*gammaaa*t36*t52)*t7-5.4d+1*gammaaa*t12*t25*t33*t53*t6*
     4         t7+5.4d+1*gammaaa*t25*t33*t53*t6*t7+5.4d+1*t12*t20*t53*t6
     5         *t7-5.4d+1*t12*t33*t51*t6*t7+5.4d+1*t33*t51*t6*t7+1.8d+1*
     6         gammaaa*t13*t3*t6*t68+1.08d+2*gammaaa*t12*t20*t25*t56*t6+
     7         5.4d+1*t13*t56*t6+1.08d+2*gammaaa*t12*t20*t25*t3*t55*t6+5
     8         .4d+1*t13*t3*t55*t6+5.4d+1*gammaaa*t13*t27*t55*t6+1.08d+2
     9         *gammaaa*t12*t20*t27*t3*t53*t6-1.08d+2*gammaaa*t12*t27*t3
     :         *t33*t51*t6+1.08d+2*gammaaa*t27*t3*t33*t51*t6+2.16d+2*t12
     ;         *t20*t25*t27*t3*t6)-2*t28*t60+2*t3*t4*t57+2*gammaaa*t27*t
     <         4*t57-t29*t57+gammaaa*t28*t4*t55+2*t27*t28*t4)+t15*t6*(-g
     =         ammaaa*t4*t68-3*t4*t55)+1.875d+0*gammaaa*t28**3*t3*t4*t6*
     >         t61+2.5d-1*t31*(-3*t3*t4*t59-3*gammaaa*t27*t4*t59-3*t28*t
     ?         58)*t6-1.5d+0*gammaaa*t28*t3*t31*t4*t57*t6)*wght+Cmat3(iq
     @         ,D3_GAA_GAA_GAA)
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
            t32 = 1/t5**7.333333333333333d+0
            t33 = 1/t5**4.666666666666666d+0
            t34 = -2.6666666666666666d+0*rhoa*t31*t32*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t33
            t35 = 1/t11
            t36 = 1.8d+1*gammaaa*t12*t34*t35*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t33*t7
            t37 = 1/t14**3
            t38 = rhoa**2.3333333333333334d+0
            t39 = rhoa**3.333333333333333d-1
            t40 = 1/t25
            t41 = 1/t20**7.333333333333333d+0
            t42 = 1/t20**4.666666666666666d+0
            t43 = -2.6666666666666666d+0*gammabb*rhob*t42-2.666666666666
     1         6666d+0*rhob*t23*t40*t41
            t44 = 1/t26
            t45 = 1.8d+1*gammabb*t21*t22*t27*t43*t44-2.4d+1*gammabb*rhob
     1         *t22*t28*t42
            t46 = 1/t29**3
            t47 = rhob**2.3333333333333334d+0
            t48 = rhob**3.333333333333333d-1
            t49 = gammaaa*t31*t9+t6
            t50 = 1/t1**2
            t51 = 9.474169d-4*t2-9.474169d-4*gammaaa*t50
            t52 = 1.8d+1*gammaaa*t12*t35*t49*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t3*t51*t6
            t53 = -gammaaa*t4*t51-t3*t4
            t54 = gammabb*t24*t40+t21
            t55 = 1/t16**2
            t56 = 9.474169d-4*t17-9.474169d-4*gammabb*t55
            t57 = 1.8d+1*gammabb*t18*t21*t28*t56+1.8d+1*gammabb*t21*t22*
     1         t27*t44*t54+9.0d+0*t21*t22*t28
            t58 = -gammabb*t19*t56-t18*t19
            t59 = t36**2
            t60 = 1/t14**5
            t61 = t34**2
            t62 = 1/t11**2
            t63 = gammaaa**4
            t64 = rhoa**2
            t65 = 1/t10**3
            t66 = 1/t5**1.4666666666666667d+1
            t67 = 1/t5**9.333333333333332d+0
            t68 = 1/t5**6.666666666666666d+0
            t69 = 1.9555555555555554d+1*t31*t64*t67*t8-2.666666666666666
     1         6d+0*t31*t32*t8+1.2444444444444443d+1*gammaaa*t64*t68-7.1
     2         1111111111111d+0*t63*t64*t65*t66-2.6666666666666666d+0*ga
     3         mmaaa*t33
            t70 = 1.8d+1*gammaaa*t12*t35*t6*t69*t7+1.1199999999999999d+2
     1         *gammaaa*t13*t64*t68*t7-1.8d+1*gammaaa*t12*t6*t61*t62*t7+
     2         1.8d+1*gammaaa*t6*t61*t62*t7-9.6d+1*gammaaa*rhoa*t12*t33*
     3         t34*t35*t7-2.4d+1*gammaaa*t13*t33*t7
            t71 = rhoa**3.3333333333333337d+0
            t72 = 1/rhoa**6.666666666666666d-1
            t73 = t45**2
            t74 = 1/t29**5
            t75 = t43**2
            t76 = 1/t26**2
            t77 = gammabb**4
            t78 = rhob**2
            t79 = 1/t25**3
            t80 = 1/t20**1.4666666666666667d+1
            t81 = 1/t20**9.333333333333332d+0
            t82 = 1/t20**6.666666666666666d+0
            t83 = 1.2444444444444443d+1*gammabb*t78*t82+1.95555555555555
     1         54d+1*t23*t40*t78*t81-7.11111111111111d+0*t77*t78*t79*t80
     2         -2.6666666666666666d+0*gammabb*t42-2.6666666666666666d+0*
     3         t23*t40*t41
            t84 = 1.8d+1*gammabb*t21*t22*t27*t44*t83+1.1199999999999999d
     1         +2*gammabb*t22*t28*t78*t82-1.8d+1*gammabb*t21*t22*t27*t75
     2         *t76+1.8d+1*gammabb*t21*t22*t75*t76-9.6d+1*gammabb*rhob*t
     3         22*t27*t42*t43*t44-2.4d+1*gammabb*t22*t28*t42
            t85 = rhob**3.3333333333333337d+0
            t86 = 1/rhob**6.666666666666666d-1
            t87 = gammaaa**3
            t88 = 1/t5**1.2666666666666665d+1
            t89 = 2.6666666666666666d+0*rhoa*t65*t87*t88-2.6666666666666
     1         666d+0*rhoa*t33-5.333333333333333d+0*gammaaa*rhoa*t31*t32
            t90 = 1.8d+1*gammaaa*t12*t35*t6*t7*t89-1.8d+1*gammaaa*t12*t3
     1         4*t49*t6*t62*t7+1.8d+1*gammaaa*t34*t49*t6*t62*t7+1.8d+1*t
     2         12*t34*t35*t6*t7-4.8d+1*gammaaa*rhoa*t12*t33*t35*t49*t7-2
     3         .4d+1*rhoa*t13*t33*t7+3.6d+1*gammaaa*t12*t3*t34*t35*t51*t
     4         6-4.8d+1*gammaaa*rhoa*t13*t3*t33*t51
            t91 = gammaaa*t3*t4*t90+gammaaa*t36*t4*t51+t3*t36*t4
            t92 = 2.6666666666666666d+0*gammaaa*t38*t51+2.66666666666666
     1         66d+0*t3*t38
            t93 = -4*gammaaa*t39*t51-4*t3*t39
            t94 = gammabb**3
            t95 = 1/t20**1.2666666666666665d+1
            t96 = 2.6666666666666666d+0*rhob*t79*t94*t95-2.6666666666666
     1         666d+0*rhob*t42-5.333333333333333d+0*gammabb*rhob*t40*t41
            t97 = 1.8d+1*gammabb*t21*t22*t27*t44*t96-1.8d+1*gammabb*t21*
     1         t22*t27*t43*t54*t76+1.8d+1*gammabb*t21*t22*t43*t54*t76+3.
     2         6d+1*gammabb*t18*t21*t27*t43*t44*t56-4.8d+1*gammabb*rhob*
     3         t18*t28*t42*t56-4.8d+1*gammabb*rhob*t22*t27*t42*t44*t54+1
     4         .8d+1*t21*t22*t27*t43*t44-2.4d+1*rhob*t22*t28*t42
            t98 = gammabb*t18*t19*t97+gammabb*t19*t45*t56+t18*t19*t45
            t99 = 2.6666666666666666d+0*gammabb*t47*t56+2.66666666666666
     1         66d+0*t18*t47
            t100 = -4*gammabb*t48*t56-4*t18*t48
            t101 = t49**2
            t102 = 1/t5**1.0666666666666666d+1
            t103 = t31*t9-t102*t65*t8
            t104 = 1/t1**3
            t105 = 1.8948337999999998d-3*gammaaa*t104-1.8948337999999998
     1         d-3*t50
            t106 = t51**2
            t107 = -1.8d+1*gammaaa*t101*t12*t6*t62*t7+1.8d+1*gammaaa*t10
     1         1*t6*t62*t7+3.6d+1*t12*t35*t49*t6*t7+1.8d+1*gammaaa*t103*
     2         t12*t35*t6*t7+7.2d+1*gammaaa*t12*t3*t35*t49*t51*t6+3.6d+1
     3         *t13*t3*t51*t6+1.8d+1*gammaaa*t105*t13*t3*t6+1.8d+1*gamma
     4         aa*t106*t13*t6
            t108 = -t52*t53+gammaaa*t4*t51*t52+t3*t4*t52+gammaaa*t107*t3
     1         *t4
            t109 = t52**2
            t110 = -2*t4*t51-gammaaa*t105*t4
            t111 = t54**2
            t112 = 1/t20**1.0666666666666666d+1
            t113 = t24*t40-t112*t23*t79
            t114 = 1/t16**3
            t115 = 1.8948337999999998d-3*gammabb*t114-1.8948337999999998
     1         d-3*t55
            t116 = t56**2
            t117 = -1.8d+1*gammabb*t111*t21*t22*t27*t76+1.8d+1*gammabb*t
     1         111*t21*t22*t76+7.2d+1*gammabb*t18*t21*t27*t44*t54*t56+3.
     2         6d+1*t18*t21*t28*t56+3.6d+1*t21*t22*t27*t44*t54+1.8d+1*ga
     3         mmabb*t113*t21*t22*t27*t44+1.8d+1*gammabb*t115*t18*t21*t2
     4         8+1.8d+1*gammabb*t116*t21*t28
            t118 = -t57*t58+gammabb*t19*t56*t57+t18*t19*t57+gammabb*t117
     1         *t18*t19
            t119 = t57**2
            t120 = -2*t19*t56-gammabb*t115*t19
            t121 = 1/t14**7
            t122 = t34**3
            t123 = 1/t11**3
            t124 = rhoa**3
            t125 = 1/t10**5
            t126 = 1/t5**8.666666666666666d+0
            t127 = 1/t29**7
            t128 = t43**3
            t129 = 1/t26**3
            t130 = rhob**3
            t131 = 1/t25**5
            t132 = 1/t20**8.666666666666666d+0
            t133 = t49**3
            t134 = 5.684501399999999d-3*t104-5.684501399999999d-3*gammaa
     1         a/t1**4
            t135 = t54**3
            t136 = 5.684501399999999d-3*t114-5.684501399999999d-3*gammab
     1         b/t16**4
            fnc(iq) = 1.0d+0*(-gammaaa*t15*t3*t4*t6-gammabb*t18*t19*t21*
     1         t30)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(5.0d-1*gammaaa*t3*t36*t37*t4*t6-1.3
     1         333333333333333d+0*gammaaa*t15*t3*t39*t6+2.66666666666666
     2         66d+0*gammaaa*t15*t3*t33*t38)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-1.3333333333333333d+0*gammabb*t18*
     1         t21*t30*t48+2.6666666666666666d+0*gammabb*t18*t30*t42*t47
     2         +5.0d-1*gammabb*t18*t19*t21*t45*t46)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 1.0d+0*(t15*t53*t6+5.0d-1*gammaaa*t3*t37*t
     1         4*t52*t6)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*(t21*t30*t58+5.0d-1*gammabb*t18*t19
     1         *t21*t46*t57)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-4.444444444444444d-1*gammaaa*t
     1         15*t3*t6*t72-1.2444444444444443d+1*gammaaa*t15*t3*t68*t71
     2         +5.0d-1*gammaaa*t3*t37*t4*t6*t70-7.5d-1*gammaaa*t3*t4*t59
     3         *t6*t60+1.3333333333333333d+0*gammaaa*t3*t36*t37*t39*t6+9
     4         .777777777777779d+0*gammaaa*t15*t3*t33*t4-2.6666666666666
     5         666d+0*gammaaa*t3*t33*t36*t37*t38)*wght+Amat2(iq,D2_RA_RA
     6         )
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-4.444444444444444d-1*gammabb*t
     1         18*t21*t30*t86-1.2444444444444443d+1*gammabb*t18*t30*t82*
     2         t85+5.0d-1*gammabb*t18*t19*t21*t46*t84-7.5d-1*gammabb*t18
     3         *t19*t21*t73*t74+1.3333333333333333d+0*gammabb*t18*t21*t4
     4         5*t46*t48-2.6666666666666666d+0*gammabb*t18*t42*t45*t46*t
     5         47+9.777777777777779d+0*gammabb*t18*t19*t30*t42)*wght+Ama
     6         t2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(3.333333333333333d-1*t15*t6*t9
     1         3+t15*t33*t92+5.0d-1*t37*t6*t91-7.5d-1*gammaaa*t3*t36*t4*
     2         t52*t6*t60+6.666666666666666d-1*gammaaa*t3*t37*t39*t52*t6
     3         -1.3333333333333333d+0*gammaaa*t3*t33*t37*t38*t52)*wght+C
     4         mat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(t30*t42*t99+5.0d-1*t21*t46*t98
     1         -7.5d-1*gammabb*t18*t19*t21*t45*t57*t74+6.666666666666666
     2         d-1*gammabb*t18*t21*t46*t48*t57-1.3333333333333333d+0*gam
     3         mabb*t18*t42*t46*t47*t57+3.333333333333333d-1*t100*t21*t3
     4         0)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*(-7.5d-1*gammaaa*t109*t3*t4*t6
     1         *t60+5.0d-1*t108*t37*t6+t110*t15*t6)*wght+Cmat2(iq,D2_GAA
     2         _GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*(-7.5d-1*gammabb*t119*t18*t19*
     1         t21*t74+5.0d-1*t118*t21*t46+t120*t21*t30)*wght+Cmat2(iq,D
     2         2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(5.0d-1*gammaaa*t3*t37*t4*t6*
     1         (1.8d+1*gammaaa*t12*t35*t6*t7*(5.866666666666666d+1*rhoa*
     2         t31*t67*t8-1.8251851851851847d+2*t124*t31*t8/t5**1.133333
     3         3333333333d+1+3.733333333333333d+1*gammaaa*rhoa*t68-2.133
     4         3333333333332d+1*rhoa*t63*t65*t66+1.5644444444444444d+2*t
     5         124*t63*t65/t5**1.6666666666666663d+1-5.688888888888889d+
     6         1*gammaaa**6*t124*t125/t5**2.2d+1-8.296296296296295d+1*ga
     7         mmaaa*t124*t126)-5.4d+1*gammaaa*t12*t34*t6*t62*t69*t7+5.4
     8         d+1*gammaaa*t34*t6*t62*t69*t7-1.44d+2*gammaaa*rhoa*t12*t3
     9         3*t35*t69*t7+6.719999999999999d+2*gammaaa*t12*t34*t35*t64
     :         *t68*t7+3.3599999999999997d+2*gammaaa*rhoa*t13*t68*t7+1.4
     ;         4d+2*gammaaa*rhoa*t12*t33*t61*t62*t7-1.44d+2*gammaaa*rhoa
     <         *t33*t61*t62*t7+3.6d+1*gammaaa*t12*t122*t123*t6*t7-5.4d+1
     =         *gammaaa*t122*t123*t6*t7-1.44d+2*gammaaa*t12*t33*t34*t35*
     >         t7-7.466666666666665d+2*gammaaa*t124*t126*t13*t7)+6.66666
     ?         6666666666d-1*gammaaa*t3*t36*t37*t6*t72+1.866666666666666
     @         5d+1*gammaaa*t3*t36*t37*t68*t71+t33*t37*(-4.0d+0*gammaaa*
     1         t3*t38*t70-1.4666666666666667d+1*gammaaa*t3*t36*t4)-2.25d
     2         +0*gammaaa*t3*t36*t4*t6*t60*t70+2*gammaaa*t3*t37*t39*t6*t
     3         70-8.711111111111112d+1*gammaaa*t15*t3*t38*t68-3*gammaaa*
     4         t3*t39*t59*t6*t60+6.0d+0*gammaaa*t3*t33*t38*t59*t60+1.875
     5         d+0*gammaaa*t121*t3*t36**3*t4*t6+2.962962962962963d-1*gam
     6         maaa*t15*t3*t6/rhoa**1.6666666666666669d+0+1.422222222222
     7         2222d+1*gammaaa*t15*t3*t33*t39+8.296296296296295d+1*gamma
     8         aa*rhoa**4.333333333333333d+0*t126*t15*t3)*wght+Amat3(iq,
     9         D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(6.666666666666666d-1*gammabb
     1         *t18*t21*t45*t46*t86+1.8666666666666665d+1*gammabb*t18*t4
     2         5*t46*t82*t85+t42*t46*(-4.0d+0*gammabb*t18*t47*t84-1.4666
     3         666666666667d+1*gammabb*t18*t19*t45)-2.25d+0*gammabb*t18*
     4         t19*t21*t45*t74*t84+2*gammabb*t18*t21*t46*t48*t84+5.0d-1*
     5         gammabb*t18*t19*t21*t46*(-5.4d+1*gammabb*t21*t22*t27*t43*
     6         t76*t83+5.4d+1*gammabb*t21*t22*t43*t76*t83-1.44d+2*gammab
     7         b*rhob*t22*t27*t42*t44*t83+1.8d+1*gammabb*t21*t22*t27*t44
     8         *(3.733333333333333d+1*gammabb*rhob*t82+5.866666666666666
     9         d+1*rhob*t23*t40*t81-2.1333333333333332d+1*rhob*t77*t79*t
     :         80+1.5644444444444444d+2*t130*t77*t79/t20**1.666666666666
     ;         6663d+1-1.8251851851851847d+2*t130*t23*t40/t20**1.1333333
     <         333333333d+1-5.688888888888889d+1*gammabb**6*t130*t131/t2
     =         0**2.2d+1-8.296296296296295d+1*gammabb*t130*t132)+6.71999
     >         9999999999d+2*gammabb*t22*t27*t43*t44*t78*t82+3.359999999
     ?         9999997d+2*gammabb*rhob*t22*t28*t82+1.44d+2*gammabb*rhob*
     @         t22*t27*t42*t75*t76-1.44d+2*gammabb*rhob*t22*t42*t75*t76-
     1         1.44d+2*gammabb*t22*t27*t42*t43*t44-7.466666666666665d+2*
     2         gammabb*t130*t132*t22*t28+3.6d+1*gammabb*t128*t129*t21*t2
     3         2*t27-5.4d+1*gammabb*t128*t129*t21*t22)-8.711111111111112
     4         d+1*gammabb*t18*t30*t47*t82-3*gammabb*t18*t21*t48*t73*t74
     5         +6.0d+0*gammabb*t18*t42*t47*t73*t74+1.4222222222222222d+1
     6         *gammabb*t18*t30*t42*t48+1.875d+0*gammabb*t127*t18*t19*t2
     7         1*t45**3+2.962962962962963d-1*gammabb*t18*t21*t30/rhob**1
     8         .6666666666666669d+0+8.296296296296295d+1*gammabb*rhob**4
     9         .333333333333333d+0*t132*t18*t30)*wght+Amat3(iq,D3_RB_RB_
     :         RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(3.333333333333333d-1*t37*t6
     1         *(4*gammaaa*t3*t39*t90+4*gammaaa*t36*t39*t51+4*t3*t36*t39
     2         )+t33*t37*(-2.6666666666666666d+0*gammaaa*t3*t38*t90-4.88
     3         88888888888893d+0*gammaaa*t3*t4*t52-2.6666666666666666d+0
     4         *gammaaa*t36*t38*t51-2.6666666666666666d+0*t3*t36*t38)-1.
     5         5d+0*gammaaa*t3*t36*t4*t6*t60*t90+5.0d-1*t37*t6*(gammaaa*
     6         t3*t4*(-3.6d+1*gammaaa*t12*t34*t6*t62*t7*t89+3.6d+1*gamma
     7         aa*t34*t6*t62*t7*t89-9.6d+1*gammaaa*rhoa*t12*t33*t35*t7*t
     8         89+1.8d+1*gammaaa*t12*t35*t6*t7*(2.6666666666666666d+0*t6
     9         5*t87*t88-2.8444444444444444d+1*t64*t65*t66*t87-1.9555555
     :         555555554d+1*t64*t65*t87/t5**1.4666666666666664d+1+1.2444
     ;         444444444443d+1*t64*t68+3.911111111111111d+1*gammaaa*t31*
     <         t64*t67+2.1333333333333332d+1*gammaaa**5*t125*t64/t5**2.0
     =         d+1-2.6666666666666666d+0*t33-5.333333333333333d+0*gammaa
     >         a*t31*t32)-1.8d+1*gammaaa*t12*t49*t6*t62*t69*t7+1.8d+1*ga
     ?         mmaaa*t49*t6*t62*t69*t7+1.8d+1*t12*t35*t6*t69*t7+2.239999
     @         9999999998d+2*gammaaa*t12*t35*t49*t64*t68*t7+1.1199999999
     1         999999d+2*t13*t64*t68*t7-1.8d+1*t12*t6*t61*t62*t7+1.8d+1*
     2         t6*t61*t62*t7+9.6d+1*gammaaa*rhoa*t12*t33*t34*t49*t62*t7-
     3         9.6d+1*gammaaa*rhoa*t33*t34*t49*t62*t7+3.6d+1*gammaaa*t12
     4         *t123*t49*t6*t61*t7-5.4d+1*gammaaa*t123*t49*t6*t61*t7-4.8
     5         d+1*gammaaa*t12*t33*t35*t49*t7-9.6d+1*rhoa*t12*t33*t34*t3
     6         5*t7-2.4d+1*t13*t33*t7+3.6d+1*gammaaa*t12*t3*t35*t51*t6*t
     7         69+2.2399999999999998d+2*gammaaa*t13*t3*t51*t64*t68-3.6d+
     8         1*gammaaa*t12*t3*t51*t6*t61*t62+3.6d+1*gammaaa*t3*t51*t6*
     9         t61*t62-1.92d+2*gammaaa*rhoa*t12*t3*t33*t34*t35*t51-4.8d+
     :         1*gammaaa*t13*t3*t33*t51)+gammaaa*t4*t51*t70+t3*t4*t70)+2
     ;         .222222222222222d-1*gammaaa*t3*t37*t52*t6*t72+1.111111111
     <         111111d-1*t15*(-4*gammaaa*t51-4*t3)*t6*t72+t15*t68*(-1.24
     =         44444444444443d+1*gammaaa*t51*t71-1.2444444444444443d+1*t
     >         3*t71)+6.222222222222221d+0*gammaaa*t3*t37*t52*t68*t71+2.
     ?         5d-1*t6*t60*(-3*gammaaa*t3*t4*t52*t70-3*gammaaa*t4*t51*t5
     @         9-3*t3*t4*t59)-2*gammaaa*t3*t36*t39*t52*t6*t60+4.0d+0*gam
     1         maaa*t3*t33*t36*t38*t52*t60+1.875d+0*gammaaa*t121*t3*t4*t
     2         52*t59*t6+t15*t33*(9.777777777777779d+0*gammaaa*t4*t51+9.
     3         777777777777779d+0*t3*t4))*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*(3.333333333333333d-1*t21*t4
     1         6*(4*gammabb*t18*t48*t97+4*gammabb*t45*t48*t56+4*t18*t45*
     2         t48)+t42*t46*(-2.6666666666666666d+0*gammabb*t18*t47*t97-
     3         4.8888888888888893d+0*gammabb*t18*t19*t57-2.6666666666666
     4         666d+0*gammabb*t45*t47*t56-2.6666666666666666d+0*t18*t45*
     5         t47)-1.5d+0*gammabb*t18*t19*t21*t45*t74*t97+5.0d-1*t21*t4
     6         6*(gammabb*t18*t19*(-3.6d+1*gammabb*t21*t22*t27*t43*t76*t
     7         96+3.6d+1*gammabb*t21*t22*t43*t76*t96-9.6d+1*gammabb*rhob
     8         *t22*t27*t42*t44*t96+1.8d+1*gammabb*t21*t22*t27*t44*(2.66
     9         66666666666666d+0*t79*t94*t95-2.8444444444444444d+1*t78*t
     :         79*t80*t94-1.9555555555555554d+1*t78*t79*t94/t20**1.46666
     ;         66666666664d+1+1.2444444444444443d+1*t78*t82+3.9111111111
     <         11111d+1*gammabb*t40*t78*t81+2.1333333333333332d+1*gammab
     =         b**5*t131*t78/t20**2.0d+1-2.6666666666666666d+0*t42-5.333
     >         333333333333d+0*gammabb*t40*t41)-1.8d+1*gammabb*t21*t22*t
     ?         27*t54*t76*t83+1.8d+1*gammabb*t21*t22*t54*t76*t83+3.6d+1*
     @         gammabb*t18*t21*t27*t44*t56*t83+1.8d+1*t21*t22*t27*t44*t8
     1         3+2.2399999999999998d+2*gammabb*t18*t28*t56*t78*t82+2.239
     2         9999999999998d+2*gammabb*t22*t27*t44*t54*t78*t82+1.119999
     3         9999999999d+2*t22*t28*t78*t82-3.6d+1*gammabb*t18*t21*t27*
     4         t56*t75*t76+3.6d+1*gammabb*t18*t21*t56*t75*t76-1.8d+1*t21
     5         *t22*t27*t75*t76+1.8d+1*t21*t22*t75*t76+9.6d+1*gammabb*rh
     6         ob*t22*t27*t42*t43*t54*t76-9.6d+1*gammabb*rhob*t22*t42*t4
     7         3*t54*t76+3.6d+1*gammabb*t129*t21*t22*t27*t54*t75-5.4d+1*
     8         gammabb*t129*t21*t22*t54*t75-1.92d+2*gammabb*rhob*t18*t27
     9         *t42*t43*t44*t56-4.8d+1*gammabb*t18*t28*t42*t56-4.8d+1*ga
     :         mmabb*t22*t27*t42*t44*t54-9.6d+1*rhob*t22*t27*t42*t43*t44
     ;         -2.4d+1*t22*t28*t42)+gammabb*t19*t56*t84+t18*t19*t84)+2.2
     <         22222222222222d-1*gammabb*t18*t21*t46*t57*t86+1.111111111
     =         111111d-1*t21*t30*(-4*gammabb*t56-4*t18)*t86+t30*t82*(-1.
     >         2444444444444443d+1*gammabb*t56*t85-1.2444444444444443d+1
     ?         *t18*t85)+6.222222222222221d+0*gammabb*t18*t46*t57*t82*t8
     @         5+2.5d-1*t21*t74*(-3*gammabb*t18*t19*t57*t84-3*gammabb*t1
     1         9*t56*t73-3*t18*t19*t73)-2*gammabb*t18*t21*t45*t48*t57*t7
     2         4+4.0d+0*gammabb*t18*t42*t45*t47*t57*t74+1.875d+0*gammabb
     3         *t127*t18*t19*t21*t57*t73+t30*t42*(9.777777777777779d+0*g
     4         ammabb*t19*t56+9.777777777777779d+0*t18*t19))*wght+Cmat3(
     5         iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(-1.6666666666666666d-1*t37
     1         *t52*t6*t93-5.0d-1*t33*t37*t52*t92+2.5d-1*t6*t60*(-3*t52*
     2         t91-3*gammaaa*t3*t4*t52*t90-3*gammaaa*t36*t4*t51*t52-3*t3
     3         *t36*t4*t52-3*gammaaa*t107*t3*t36*t4)+5.0d-1*t37*t6*(2*ga
     4         mmaaa*t4*t51*t90+2*t3*t4*t90+gammaaa*t3*t4*(-3.6d+1*gamma
     5         aa*t12*t49*t6*t62*t7*t89+3.6d+1*gammaaa*t49*t6*t62*t7*t89
     6         +3.6d+1*t12*t35*t6*t7*t89+7.2d+1*gammaaa*t12*t3*t35*t51*t
     7         6*t89+1.8d+1*gammaaa*t12*t35*t6*t7*(1.3333333333333333d+1
     8         *rhoa*t65*t8*t88-8.0d+0*rhoa*t125*t63/t5**1.8d+1-5.333333
     9         333333333d+0*rhoa*t31*t32)-3.6d+1*t12*t34*t49*t6*t62*t7+3
     :         .6d+1*t34*t49*t6*t62*t7-1.8d+1*gammaaa*t103*t12*t34*t6*t6
     ;         2*t7+1.8d+1*gammaaa*t103*t34*t6*t62*t7+4.8d+1*gammaaa*rho
     <         a*t101*t12*t33*t62*t7-4.8d+1*gammaaa*rhoa*t101*t33*t62*t7
     =         +3.6d+1*gammaaa*t101*t12*t123*t34*t6*t7-5.4d+1*gammaaa*t1
     >         01*t123*t34*t6*t7-9.6d+1*rhoa*t12*t33*t35*t49*t7-4.8d+1*g
     ?         ammaaa*rhoa*t103*t12*t33*t35*t7-7.2d+1*gammaaa*t12*t3*t34
     @         *t49*t51*t6*t62+7.2d+1*gammaaa*t3*t34*t49*t51*t6*t62+7.2d
     1         +1*t12*t3*t34*t35*t51*t6+3.6d+1*gammaaa*t105*t12*t3*t34*t
     2         35*t6+3.6d+1*gammaaa*t106*t12*t34*t35*t6-1.92d+2*gammaaa*
     3         rhoa*t12*t3*t33*t35*t49*t51-9.6d+1*rhoa*t13*t3*t33*t51-4.
     4         8d+1*gammaaa*rhoa*t105*t13*t3*t33-4.8d+1*gammaaa*rhoa*t10
     5         6*t13*t33)+2*t36*t4*t51+gammaaa*t105*t36*t4)-gammaaa*t109
     6         *t3*t39*t6*t60+2.0d+0*gammaaa*t109*t3*t33*t38*t60+3.33333
     7         3333333333d-1*t37*(2*gammaaa*t39*t51*t52+2*t3*t39*t52+2*g
     8         ammaaa*t107*t3*t39)*t6+3.333333333333333d-1*t15*(-8*t39*t
     9         51-4*gammaaa*t105*t39)*t6+1.875d+0*gammaaa*t109*t121*t3*t
     :         36*t4*t6+t33*t37*(-1.3333333333333333d+0*gammaaa*t38*t51*
     ;         t52-1.3333333333333333d+0*t3*t38*t52-1.3333333333333333d+
     <         0*gammaaa*t107*t3*t38)+t15*t33*(5.333333333333333d+0*t38*
     =         t51+2.6666666666666666d+0*gammaaa*t105*t38))*wght+Cmat3(i
     >         q,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*(-5.0d-1*t42*t46*t57*t99+2.
     1         5d-1*t21*t74*(-3*t57*t98-3*gammabb*t18*t19*t57*t97-3*gamm
     2         abb*t19*t45*t56*t57-3*t18*t19*t45*t57-3*gammabb*t117*t18*
     3         t19*t45)+5.0d-1*t21*t46*(2*gammabb*t19*t56*t97+2*t18*t19*
     4         t97+gammabb*t18*t19*(-3.6d+1*gammabb*t21*t22*t27*t54*t76*
     5         t96+3.6d+1*gammabb*t21*t22*t54*t76*t96+7.2d+1*gammabb*t18
     6         *t21*t27*t44*t56*t96+3.6d+1*t21*t22*t27*t44*t96+1.8d+1*ga
     7         mmabb*t21*t22*t27*t44*(1.3333333333333333d+1*rhob*t23*t79
     8         *t95-8.0d+0*rhob*t131*t77/t20**1.8d+1-5.333333333333333d+
     9         0*rhob*t40*t41)-7.2d+1*gammabb*t18*t21*t27*t43*t54*t56*t7
     :         6+7.2d+1*gammabb*t18*t21*t43*t54*t56*t76-3.6d+1*t21*t22*t
     ;         27*t43*t54*t76+3.6d+1*t21*t22*t43*t54*t76-1.8d+1*gammabb*
     <         t113*t21*t22*t27*t43*t76+1.8d+1*gammabb*t113*t21*t22*t43*
     =         t76+4.8d+1*gammabb*rhob*t111*t22*t27*t42*t76-4.8d+1*gamma
     >         bb*rhob*t111*t22*t42*t76-1.92d+2*gammabb*rhob*t18*t27*t42
     ?         *t44*t54*t56+7.2d+1*t18*t21*t27*t43*t44*t56-9.6d+1*rhob*t
     @         18*t28*t42*t56-9.6d+1*rhob*t22*t27*t42*t44*t54+3.6d+1*gam
     1         mabb*t115*t18*t21*t27*t43*t44+3.6d+1*gammabb*t116*t21*t27
     2         *t43*t44-4.8d+1*gammabb*rhob*t113*t22*t27*t42*t44+3.6d+1*
     3         gammabb*t111*t129*t21*t22*t27*t43-5.4d+1*gammabb*t111*t12
     4         9*t21*t22*t43-4.8d+1*gammabb*rhob*t115*t18*t28*t42-4.8d+1
     5         *gammabb*rhob*t116*t28*t42)+2*t19*t45*t56+gammabb*t115*t1
     6         9*t45)-gammabb*t119*t18*t21*t48*t74+2.0d+0*gammabb*t119*t
     7         18*t42*t47*t74+3.333333333333333d-1*t21*t46*(2*gammabb*t4
     8         8*t56*t57+2*t18*t48*t57+2*gammabb*t117*t18*t48)+t42*t46*(
     9         -1.3333333333333333d+0*gammabb*t47*t56*t57-1.333333333333
     :         3333d+0*t18*t47*t57-1.3333333333333333d+0*gammabb*t117*t1
     ;         8*t47)-1.6666666666666666d-1*t100*t21*t46*t57+3.333333333
     <         333333d-1*t21*t30*(-8*t48*t56-4*gammabb*t115*t48)+t30*t42
     =         *(5.333333333333333d+0*t47*t56+2.6666666666666666d+0*gamm
     >         abb*t115*t47)+1.875d+0*gammabb*t119*t127*t18*t19*t21*t45)
     ?         *wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*(5.0d-1*t37*t6*(gammaaa*t3
     1         *t4*(1.8d+1*gammaaa*t12*t35*t6*t7*(3*t125*t87/t5**1.6d+1-
     2         3*gammaaa*t102*t65)-5.4d+1*gammaaa*t103*t12*t49*t6*t62*t7
     3         +5.4d+1*gammaaa*t103*t49*t6*t62*t7-5.4d+1*t101*t12*t6*t62
     4         *t7+5.4d+1*t101*t6*t62*t7+5.4d+1*t103*t12*t35*t6*t7+3.6d+
     5         1*gammaaa*t12*t123*t133*t6*t7-5.4d+1*gammaaa*t123*t133*t6
     6         *t7-1.08d+2*gammaaa*t101*t12*t3*t51*t6*t62+1.08d+2*gammaa
     7         a*t101*t3*t51*t6*t62+2.16d+2*t12*t3*t35*t49*t51*t6+1.08d+
     8         2*gammaaa*t103*t12*t3*t35*t51*t6+5.4d+1*gammaaa*t105*t13*
     9         t51*t6+1.08d+2*gammaaa*t105*t12*t3*t35*t49*t6+1.08d+2*gam
     :         maaa*t106*t12*t35*t49*t6+1.8d+1*gammaaa*t13*t134*t3*t6+5.
     ;         4d+1*t105*t13*t3*t6+5.4d+1*t106*t13*t6)-t107*t53+2*t4*t51
     <         *t52+gammaaa*t105*t4*t52-2*t110*t52+2*gammaaa*t107*t4*t51
     =         +2*t107*t3*t4)+2.5d-1*(-3*t108*t52-3*gammaaa*t109*t4*t51-
     >         3*t109*t3*t4)*t6*t60-1.5d+0*gammaaa*t107*t3*t4*t52*t6*t60
     ?         +1.875d+0*gammaaa*t121*t3*t4*t52**3*t6+t15*(-gammaaa*t134
     @         *t4-3*t105*t4)*t6)*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 1.0d+0*(5.0d-1*t21*t46*(gammabb*t
     1         18*t19*(1.8d+1*gammabb*t21*t22*t27*t44*(3*t131*t94/t20**1
     2         .6d+1-3*gammabb*t112*t79)-1.08d+2*gammabb*t111*t18*t21*t2
     3         7*t56*t76+1.08d+2*gammabb*t111*t18*t21*t56*t76-5.4d+1*gam
     4         mabb*t113*t21*t22*t27*t54*t76+5.4d+1*gammabb*t113*t21*t22
     5         *t54*t76-5.4d+1*t111*t21*t22*t27*t76+5.4d+1*t111*t21*t22*
     6         t76+2.16d+2*t18*t21*t27*t44*t54*t56+1.08d+2*gammabb*t113*
     7         t18*t21*t27*t44*t56+5.4d+1*gammabb*t115*t21*t28*t56+1.08d
     8         +2*gammabb*t115*t18*t21*t27*t44*t54+1.08d+2*gammabb*t116*
     9         t21*t27*t44*t54+5.4d+1*t113*t21*t22*t27*t44+1.8d+1*gammab
     :         b*t136*t18*t21*t28+5.4d+1*t115*t18*t21*t28+5.4d+1*t116*t2
     ;         1*t28+3.6d+1*gammabb*t129*t135*t21*t22*t27-5.4d+1*gammabb
     <         *t129*t135*t21*t22)-t117*t58+2*t19*t56*t57+gammabb*t115*t
     =         19*t57-2*t120*t57+2*gammabb*t117*t19*t56+2*t117*t18*t19)+
     >         2.5d-1*t21*(-3*t118*t57-3*gammabb*t119*t19*t56-3*t119*t18
     ?         *t19)*t74-1.5d+0*gammabb*t117*t18*t19*t21*t57*t74+1.875d+
     @         0*gammabb*t127*t18*t19*t21*t57**3+(-gammabb*t136*t19-3*t1
     1         15*t19)*t21*t30)*wght+Cmat3(iq,D3_GBB_GBB_GBB)
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhoa*t16*t17*t8-2.6666666666666
     1         666d+0*gammaaa*rhoa*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammaaa*t12*t19*t20*t6*t7-2.4d+1*gammaaa*rhoa*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhoa**2.3333333333333334d+0
            t24 = rhoa**3.333333333333333d-1
            t25 = gammaaa*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammaaa*t26
            t28 = 1.8d+1*gammaaa*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammaaa*t13*t27*t3*t6
            t29 = -1.0d+0*t3*t4*wght-1.0d+0*gammaaa*t27*t4*wght
            t30 = t21**2
            t31 = 1/t14**5
            t32 = rhoa**3.3333333333333337d+0
            t33 = 1/t5**6.666666666666666d+0
            t34 = 1/rhoa**6.666666666666666d-1
            t35 = t19**2
            t36 = 1/t11**2
            t37 = gammaaa**4
            t38 = rhoa**2
            t39 = 1/t10**3
            t40 = 1/t5**1.4666666666666667d+1
            t41 = 1/t5**9.333333333333332d+0
            t42 = 1.9555555555555554d+1*t16*t38*t41*t8-2.666666666666666
     1         6d+0*t16*t17*t8-7.11111111111111d+0*t37*t38*t39*t40+1.244
     2         4444444444443d+1*gammaaa*t33*t38-2.6666666666666666d+0*ga
     3         mmaaa*t18
            t43 = 1.8d+1*gammaaa*t12*t20*t42*t6*t7-1.8d+1*gammaaa*t12*t3
     1         5*t36*t6*t7+1.8d+1*gammaaa*t35*t36*t6*t7+1.11999999999999
     2         99d+2*gammaaa*t13*t33*t38*t7-9.6d+1*gammaaa*rhoa*t12*t18*
     3         t19*t20*t7-2.4d+1*gammaaa*t13*t18*t7
            t44 = 5.0d-1*gammaaa*t3*t4*t43*wght+1.3333333333333333d+0*ga
     1         mmaaa*t21*t24*t3*wght
            t45 = -1.3333333333333333d+0*t24*t3*wght-1.3333333333333333d
     1         +0*gammaaa*t24*t27*wght
            t46 = 2.6666666666666666d+0*t23*t3*wght+2.6666666666666666d+
     1         0*gammaaa*t23*t27*wght
            t47 = gammaaa**3
            t48 = 1/t5**1.2666666666666665d+1
            t49 = 2.6666666666666666d+0*rhoa*t39*t47*t48-2.6666666666666
     1         666d+0*rhoa*t18-5.333333333333333d+0*gammaaa*rhoa*t16*t17
            t50 = 1.8d+1*gammaaa*t12*t20*t49*t6*t7-1.8d+1*gammaaa*t12*t1
     1         9*t25*t36*t6*t7+1.8d+1*gammaaa*t19*t25*t36*t6*t7+1.8d+1*t
     2         12*t19*t20*t6*t7-4.8d+1*gammaaa*rhoa*t12*t18*t20*t25*t7-2
     3         .4d+1*rhoa*t13*t18*t7+3.6d+1*gammaaa*t12*t19*t20*t27*t3*t
     4         6-4.8d+1*gammaaa*rhoa*t13*t18*t27*t3
            t51 = 5.0d-1*gammaaa*t3*t4*t50*wght+5.0d-1*t21*t3*t4*wght+5.
     1         0d-1*gammaaa*t21*t27*t4*wght+6.666666666666666d-1*gammaaa
     2         *t24*t28*t3*wght
            t52 = t28**2
            t53 = 1/t1**3
            t54 = 1.8948337999999998d-3*gammaaa*t53-1.8948337999999998d-
     1         3*t26
            t55 = -1.0d+0*gammaaa*t4*t54*wght-2.0d+0*t27*t4*wght
            t56 = t25**2
            t57 = 1/t5**1.0666666666666666d+1
            t58 = t16*t9-t39*t57*t8
            t59 = t27**2
            t60 = 1.8d+1*gammaaa*t12*t20*t58*t6*t7-1.8d+1*gammaaa*t12*t3
     1         6*t56*t6*t7+1.8d+1*gammaaa*t36*t56*t6*t7+3.6d+1*t12*t20*t
     2         25*t6*t7+1.8d+1*gammaaa*t13*t59*t6+1.8d+1*gammaaa*t13*t3*
     3         t54*t6+7.2d+1*gammaaa*t12*t20*t25*t27*t3*t6+3.6d+1*t13*t2
     4         7*t3*t6
            t61 = 5.0d-1*gammaaa*t3*t4*t60*wght+5.0d-1*t28*t3*t4*wght+5.
     1         0d-1*gammaaa*t27*t28*t4*wght
            t62 = 1/t14**7
            t63 = 1/t5**8.666666666666666d+0
            t64 = t19**3
            t65 = 1/t11**3
            t66 = rhoa**3
            t67 = 1/t10**5
            t68 = 5.684501399999999d-3*t53-5.684501399999999d-3*gammaaa/
     1         t1**4
            t69 = t25**3
            fnc(iq) = fnc(iq)-1.0d+0*gammaaa*t15*t3*t4*t6*wght
            Amat(iq,D1_RA) = 5.0d-1*gammaaa*t21*t22*t3*t4*t6*wght-1.3333
     1         333333333333d+0*gammaaa*t15*t24*t3*t6*wght+2.666666666666
     2         6666d+0*gammaaa*t15*t18*t23*t3*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 5.0d-1*gammaaa*t22*t28*t3*t4*t6*wght+t15*t
     1         29*t6+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = -7.5d-1*gammaaa*t3*t30*t31*t4*t6*wght-4
     1         .444444444444444d-1*gammaaa*t15*t3*t34*t6*wght+9.77777777
     2         7777779d+0*gammaaa*t15*t18*t3*t4*wght-1.2444444444444443d
     3         +1*gammaaa*t15*t3*t32*t33*wght-2.6666666666666666d+0*gamm
     4         aaa*t18*t21*t22*t23*t3*wght+t22*t44*t6+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -7.5d-1*gammaaa*t21*t28*t3*t31*t4*t6*w
     1         ght-1.3333333333333333d+0*gammaaa*t18*t22*t23*t28*t3*wght
     2         +t22*t51*t6+t15*t45*t6+t15*t18*t46+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -7.5d-1*gammaaa*t3*t31*t4*t52*t6*wght
     1         +t22*t6*t61+t15*t55*t6-5.0d-1*t22*t28*t29*t6+Cmat2(iq,D2_
     2         GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d+1
     1         *gammaaa*t12*t20*t6*t7*(-1.8251851851851847d+2*t16*t66*t8
     2         /t5**1.1333333333333333d+1+5.866666666666666d+1*rhoa*t16*
     3         t41*t8-5.688888888888889d+1*gammaaa**6*t66*t67/t5**2.2d+1
     4         -8.296296296296295d+1*gammaaa*t63*t66+1.5644444444444444d
     5         +2*t37*t39*t66/t5**1.6666666666666663d+1-2.13333333333333
     6         32d+1*rhoa*t37*t39*t40+3.733333333333333d+1*gammaaa*rhoa*
     7         t33)-7.466666666666665d+2*gammaaa*t13*t63*t66*t7+3.6d+1*g
     8         ammaaa*t12*t6*t64*t65*t7-5.4d+1*gammaaa*t6*t64*t65*t7-5.4
     9         d+1*gammaaa*t12*t19*t36*t42*t6*t7+5.4d+1*gammaaa*t19*t36*
     :         t42*t6*t7-1.44d+2*gammaaa*rhoa*t12*t18*t20*t42*t7+6.71999
     ;         9999999999d+2*gammaaa*t12*t19*t20*t33*t38*t7+1.44d+2*gamm
     <         aaa*rhoa*t12*t18*t35*t36*t7-1.44d+2*gammaaa*rhoa*t18*t35*
     =         t36*t7+3.3599999999999997d+2*gammaaa*rhoa*t13*t33*t7-1.44
     >         d+2*gammaaa*t12*t18*t19*t20*t7)*wght+2.0d+0*gammaaa*t24*t
     ?         3*t43*wght+4.444444444444444d-1*gammaaa*t21*t3*t34*wght)+
     @         t31*t6*(-1.5d+0*gammaaa*t21*t3*t4*t43*wght-1.0d+0*gammaaa
     1         *t24*t3*t30*wght)+t18*t22*(-2.6666666666666666d+0*gammaaa
     2         *t23*t3*t43*wght-1.111111111111111d+1*gammaaa*t21*t3*t4*w
     3         ght-2.6666666666666666d+0*rhoa*t44)+8.296296296296295d+1*
     4         gammaaa*rhoa**4.333333333333333d+0*t15*t3*t63*wght+1.875d
     5         +0*gammaaa*t21**3*t3*t4*t6*t62*wght+2.222222222222222d-1*
     6         gammaaa*t21*t22*t3*t34*t6*wght+2.962962962962963d-1*gamma
     7         aa*t15*t3*t6*wght/rhoa**1.6666666666666669d+0+1.866666666
     8         6666665d+1*gammaaa*t21*t22*t3*t32*t33*wght-8.711111111111
     9         112d+1*gammaaa*t15*t23*t3*t33*wght+6.0d+0*gammaaa*t18*t23
     :         *t3*t30*t31*wght+1.4222222222222222d+1*gammaaa*t15*t18*t2
     ;         4*t3*wght-1.5d+0*t21*t31*t44*t6+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d+
     1         1*gammaaa*t12*t20*t6*(2.1333333333333332d+1*gammaaa**5*t3
     2         8*t67/t5**2.0d+1-1.9555555555555554d+1*t38*t39*t47/t5**1.
     3         4666666666666664d+1+2.6666666666666666d+0*t39*t47*t48-2.8
     4         444444444444444d+1*t38*t39*t40*t47+3.911111111111111d+1*g
     5         ammaaa*t16*t38*t41+1.2444444444444443d+1*t33*t38-2.666666
     6         6666666666d+0*t18-5.333333333333333d+0*gammaaa*t16*t17)*t
     7         7+3.6d+1*gammaaa*t12*t25*t35*t6*t65*t7-5.4d+1*gammaaa*t25
     8         *t35*t6*t65*t7-3.6d+1*gammaaa*t12*t19*t36*t49*t6*t7+3.6d+
     9         1*gammaaa*t19*t36*t49*t6*t7-1.8d+1*gammaaa*t12*t25*t36*t4
     :         2*t6*t7+1.8d+1*gammaaa*t25*t36*t42*t6*t7+1.8d+1*t12*t20*t
     ;         42*t6*t7-1.8d+1*t12*t35*t36*t6*t7+1.8d+1*t35*t36*t6*t7-9.
     <         6d+1*gammaaa*rhoa*t12*t18*t20*t49*t7+2.2399999999999998d+
     =         2*gammaaa*t12*t20*t25*t33*t38*t7+1.1199999999999999d+2*t1
     >         3*t33*t38*t7+9.6d+1*gammaaa*rhoa*t12*t18*t19*t25*t36*t7-9
     ?         .6d+1*gammaaa*rhoa*t18*t19*t25*t36*t7-4.8d+1*gammaaa*t12*
     @         t18*t20*t25*t7-9.6d+1*rhoa*t12*t18*t19*t20*t7-2.4d+1*t13*
     1         t18*t7+3.6d+1*gammaaa*t12*t20*t27*t3*t42*t6-3.6d+1*gammaa
     2         a*t12*t27*t3*t35*t36*t6+3.6d+1*gammaaa*t27*t3*t35*t36*t6+
     3         2.2399999999999998d+2*gammaaa*t13*t27*t3*t33*t38-1.92d+2*
     4         gammaaa*rhoa*t12*t18*t19*t20*t27*t3-4.8d+1*gammaaa*t13*t1
     5         8*t27*t3)*wght+1.3333333333333333d+0*gammaaa*t24*t3*t50*w
     6         ght+5.0d-1*t3*t4*t43*wght+5.0d-1*gammaaa*t27*t4*t43*wght+
     7         1.3333333333333333d+0*t21*t24*t3*wght+1.3333333333333333d
     8         +0*gammaaa*t21*t24*t27*wght)+t31*t6*(-1.5d+0*gammaaa*t21*
     9         t3*t4*t50*wght-7.5d-1*t3*t30*t4*wght-7.5d-1*gammaaa*t27*t
     :         30*t4*wght)+t18*t22*(-2.6666666666666666d+0*gammaaa*t23*t
     ;         3*t50*wght-4.8888888888888893d+0*gammaaa*t28*t3*t4*wght-2
     <         .6666666666666666d+0*t21*t23*t3*wght-2.6666666666666666d+
     =         0*gammaaa*t21*t23*t27*wght)+t15*t18*(9.777777777777779d+0
     >         *t3*t4*wght+9.777777777777779d+0*gammaaa*t27*t4*wght)+t15
     ?         *t33*(-1.2444444444444443d+1*t3*t32*wght-1.24444444444444
     @         43d+1*gammaaa*t27*t32*wght)+t15*t34*t6*(-4.44444444444444
     1         4d-1*t3*wght-4.444444444444444d-1*gammaaa*t27*wght)+1.875
     2         d+0*gammaaa*t28*t3*t30*t4*t6*t62*wght+2.222222222222222d-
     3         1*gammaaa*t22*t28*t3*t34*t6*wght+6.222222222222221d+0*gam
     4         maaa*t22*t28*t3*t32*t33*wght+4.0d+0*gammaaa*t18*t21*t23*t
     5         28*t3*t31*wght-1.5d+0*t28*t31*t44*t6+Cmat3(iq,D3_RA_RA_GA
     6         A)
            Cmat3(iq,D3_RA_GAA_GAA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(1.8d
     1         +1*gammaaa*t12*t20*t6*t7*(1.3333333333333333d+1*rhoa*t39*
     2         t48*t8-8.0d+0*rhoa*t37*t67/t5**1.8d+1-5.333333333333333d+
     3         0*rhoa*t16*t17)+3.6d+1*gammaaa*t12*t19*t56*t6*t65*t7-5.4d
     4         +1*gammaaa*t19*t56*t6*t65*t7-1.8d+1*gammaaa*t12*t19*t36*t
     5         58*t6*t7+1.8d+1*gammaaa*t19*t36*t58*t6*t7-3.6d+1*gammaaa*
     6         t12*t25*t36*t49*t6*t7+3.6d+1*gammaaa*t25*t36*t49*t6*t7+3.
     7         6d+1*t12*t20*t49*t6*t7-3.6d+1*t12*t19*t25*t36*t6*t7+3.6d+
     8         1*t19*t25*t36*t6*t7-4.8d+1*gammaaa*rhoa*t12*t18*t20*t58*t
     9         7+4.8d+1*gammaaa*rhoa*t12*t18*t36*t56*t7-4.8d+1*gammaaa*r
     :         hoa*t18*t36*t56*t7-9.6d+1*rhoa*t12*t18*t20*t25*t7+3.6d+1*
     ;         gammaaa*t12*t19*t20*t59*t6+3.6d+1*gammaaa*t12*t19*t20*t3*
     <         t54*t6+7.2d+1*gammaaa*t12*t20*t27*t3*t49*t6-7.2d+1*gammaa
     =         a*t12*t19*t25*t27*t3*t36*t6+7.2d+1*gammaaa*t19*t25*t27*t3
     >         *t36*t6+7.2d+1*t12*t19*t20*t27*t3*t6-4.8d+1*gammaaa*rhoa*
     ?         t13*t18*t59-4.8d+1*gammaaa*rhoa*t13*t18*t3*t54-1.92d+2*ga
     @         mmaaa*rhoa*t12*t18*t20*t25*t27*t3-9.6d+1*rhoa*t13*t18*t27
     1         *t3)*wght+6.666666666666666d-1*gammaaa*t24*t3*t60*wght+5.
     2         0d-1*gammaaa*t21*t4*t54*wght+1.0d+0*t3*t4*t50*wght+1.0d+0
     3         *gammaaa*t27*t4*t50*wght+1.0d+0*t21*t27*t4*wght+6.6666666
     4         66666666d-1*t24*t28*t3*wght+6.666666666666666d-1*gammaaa*
     5         t24*t27*t28*wght)+t31*t6*(-7.5d-1*gammaaa*t21*t3*t4*t60*w
     6         ght-7.5d-1*gammaaa*t28*t3*t4*t50*wght-7.5d-1*t21*t28*t3*t
     7         4*wght-7.5d-1*gammaaa*t21*t27*t28*t4*wght)+t18*t22*(-1.33
     8         33333333333333d+0*gammaaa*t23*t3*t60*wght-1.3333333333333
     9         333d+0*t23*t28*t3*wght-1.3333333333333333d+0*gammaaa*t23*
     :         t27*t28*wght)+t15*t6*(-1.3333333333333333d+0*gammaaa*t24*
     ;         t54*wght-2.6666666666666666d+0*t24*t27*wght)+t15*t18*(2.6
     <         666666666666666d+0*gammaaa*t23*t54*wght+5.333333333333333
     =         d+0*t23*t27*wght)+1.875d+0*gammaaa*t21*t3*t4*t52*t6*t62*w
     >         ght+2.0d+0*gammaaa*t18*t23*t3*t31*t52*wght-1.5d+0*t28*t31
     ?         *t51*t6-5.0d-1*t22*t28*t45*t6-5.0d-1*t18*t22*t28*t46+Cmat
     @         3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t22*t6*(5.0d-1*gammaaa*t3*t4*(3.6
     1         d+1*gammaaa*t12*t6*t65*t69*t7-5.4d+1*gammaaa*t6*t65*t69*t
     2         7+1.8d+1*gammaaa*t12*t20*t6*(3*t47*t67/t5**1.6d+1-3*gamma
     3         aa*t39*t57)*t7-5.4d+1*gammaaa*t12*t25*t36*t58*t6*t7+5.4d+
     4         1*gammaaa*t25*t36*t58*t6*t7+5.4d+1*t12*t20*t58*t6*t7-5.4d
     5         +1*t12*t36*t56*t6*t7+5.4d+1*t36*t56*t6*t7+1.8d+1*gammaaa*
     6         t13*t3*t6*t68+1.08d+2*gammaaa*t12*t20*t25*t59*t6+5.4d+1*t
     7         13*t59*t6+1.08d+2*gammaaa*t12*t20*t27*t3*t58*t6-1.08d+2*g
     8         ammaaa*t12*t27*t3*t36*t56*t6+1.08d+2*gammaaa*t27*t3*t36*t
     9         56*t6+1.08d+2*gammaaa*t12*t20*t25*t3*t54*t6+5.4d+1*t13*t3
     :         *t54*t6+5.4d+1*gammaaa*t13*t27*t54*t6+2.16d+2*t12*t20*t25
     ;         *t27*t3*t6)*wght+1.0d+0*t3*t4*t60*wght+1.0d+0*gammaaa*t27
     <         *t4*t60*wght+5.0d-1*gammaaa*t28*t4*t54*wght+1.0d+0*t27*t2
     =         8*t4*wght-t28*t55)+t15*t6*(-1.0d+0*gammaaa*t4*t68*wght-3.
     >         0d+0*t4*t54*wght)+t31*t6*(-1.5d+0*gammaaa*t28*t3*t4*t60*w
     ?         ght-7.5d-1*t3*t4*t52*wght-7.5d-1*gammaaa*t27*t4*t52*wght)
     @         +1.875d+0*gammaaa*t28**3*t3*t4*t6*t62*wght-1.5d+0*t28*t31
     1         *t6*t61-5.0d-1*t22*t29*t6*t60+7.5d-1*t29*t31*t52*t6+Cmat3
     2         (iq,D3_GAA_GAA_GAA)
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
            t17 = 1/t5**7.333333333333333d+0
            t18 = 1/t5**4.666666666666666d+0
            t19 = -2.6666666666666666d+0*rhob*t16*t17*t8-2.6666666666666
     1         666d+0*gammabb*rhob*t18
            t20 = 1/t11
            t21 = 1.8d+1*gammabb*t12*t19*t20*t6*t7-2.4d+1*gammabb*rhob*t
     1         13*t18*t7
            t22 = 1/t14**3
            t23 = rhob**2.3333333333333334d+0
            t24 = rhob**3.333333333333333d-1
            t25 = gammabb*t16*t9+t6
            t26 = 1/t1**2
            t27 = 9.474169d-4*t2-9.474169d-4*gammabb*t26
            t28 = 1.8d+1*gammabb*t12*t20*t25*t6*t7+9.0d+0*t13*t6*t7+1.8d
     1         +1*gammabb*t13*t27*t3*t6
            t29 = -1.0d+0*t3*t4*wght-1.0d+0*gammabb*t27*t4*wght
            t30 = t21**2
            t31 = 1/t14**5
            t32 = rhob**3.3333333333333337d+0
            t33 = 1/t5**6.666666666666666d+0
            t34 = 1/rhob**6.666666666666666d-1
            t35 = t19**2
            t36 = 1/t11**2
            t37 = gammabb**4
            t38 = rhob**2
            t39 = 1/t10**3
            t40 = 1/t5**1.4666666666666667d+1
            t41 = 1/t5**9.333333333333332d+0
            t42 = 1.9555555555555554d+1*t16*t38*t41*t8-2.666666666666666
     1         6d+0*t16*t17*t8-7.11111111111111d+0*t37*t38*t39*t40+1.244
     2         4444444444443d+1*gammabb*t33*t38-2.6666666666666666d+0*ga
     3         mmabb*t18
            t43 = 1.8d+1*gammabb*t12*t20*t42*t6*t7-1.8d+1*gammabb*t12*t3
     1         5*t36*t6*t7+1.8d+1*gammabb*t35*t36*t6*t7+1.11999999999999
     2         99d+2*gammabb*t13*t33*t38*t7-9.6d+1*gammabb*rhob*t12*t18*
     3         t19*t20*t7-2.4d+1*gammabb*t13*t18*t7
            t44 = 5.0d-1*gammabb*t3*t4*t43*wght+1.3333333333333333d+0*ga
     1         mmabb*t21*t24*t3*wght
            t45 = -1.3333333333333333d+0*t24*t3*wght-1.3333333333333333d
     1         +0*gammabb*t24*t27*wght
            t46 = 2.6666666666666666d+0*t23*t3*wght+2.6666666666666666d+
     1         0*gammabb*t23*t27*wght
            t47 = gammabb**3
            t48 = 1/t5**1.2666666666666665d+1
            t49 = 2.6666666666666666d+0*rhob*t39*t47*t48-2.6666666666666
     1         666d+0*rhob*t18-5.333333333333333d+0*gammabb*rhob*t16*t17
            t50 = 1.8d+1*gammabb*t12*t20*t49*t6*t7-1.8d+1*gammabb*t12*t1
     1         9*t25*t36*t6*t7+1.8d+1*gammabb*t19*t25*t36*t6*t7+1.8d+1*t
     2         12*t19*t20*t6*t7-4.8d+1*gammabb*rhob*t12*t18*t20*t25*t7-2
     3         .4d+1*rhob*t13*t18*t7+3.6d+1*gammabb*t12*t19*t20*t27*t3*t
     4         6-4.8d+1*gammabb*rhob*t13*t18*t27*t3
            t51 = 5.0d-1*gammabb*t3*t4*t50*wght+5.0d-1*t21*t3*t4*wght+5.
     1         0d-1*gammabb*t21*t27*t4*wght+6.666666666666666d-1*gammabb
     2         *t24*t28*t3*wght
            t52 = t28**2
            t53 = 1/t1**3
            t54 = 1.8948337999999998d-3*gammabb*t53-1.8948337999999998d-
     1         3*t26
            t55 = -1.0d+0*gammabb*t4*t54*wght-2.0d+0*t27*t4*wght
            t56 = t25**2
            t57 = 1/t5**1.0666666666666666d+1
            t58 = t16*t9-t39*t57*t8
            t59 = t27**2
            t60 = 1.8d+1*gammabb*t12*t20*t58*t6*t7-1.8d+1*gammabb*t12*t3
     1         6*t56*t6*t7+1.8d+1*gammabb*t36*t56*t6*t7+3.6d+1*t12*t20*t
     2         25*t6*t7+1.8d+1*gammabb*t13*t59*t6+1.8d+1*gammabb*t13*t3*
     3         t54*t6+7.2d+1*gammabb*t12*t20*t25*t27*t3*t6+3.6d+1*t13*t2
     4         7*t3*t6
            t61 = 5.0d-1*gammabb*t3*t4*t60*wght+5.0d-1*t28*t3*t4*wght+5.
     1         0d-1*gammabb*t27*t28*t4*wght
            t62 = 1/t14**7
            t63 = 1/t5**8.666666666666666d+0
            t64 = t19**3
            t65 = 1/t11**3
            t66 = rhob**3
            t67 = 1/t10**5
            t68 = 5.684501399999999d-3*t53-5.684501399999999d-3*gammabb/
     1         t1**4
            t69 = t25**3
            fnc(iq) = fnc(iq)-1.0d+0*gammabb*t15*t3*t4*t6*wght
            Amat(iq,D1_RB) = 5.0d-1*gammabb*t21*t22*t3*t4*t6*wght-1.3333
     1         333333333333d+0*gammabb*t15*t24*t3*t6*wght+2.666666666666
     2         6666d+0*gammabb*t15*t18*t23*t3*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 5.0d-1*gammabb*t22*t28*t3*t4*t6*wght+t15*t
     1         29*t6+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = -7.5d-1*gammabb*t3*t30*t31*t4*t6*wght-4
     1         .444444444444444d-1*gammabb*t15*t3*t34*t6*wght+9.77777777
     2         7777779d+0*gammabb*t15*t18*t3*t4*wght-1.2444444444444443d
     3         +1*gammabb*t15*t3*t32*t33*wght-2.6666666666666666d+0*gamm
     4         abb*t18*t21*t22*t23*t3*wght+t22*t44*t6+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -7.5d-1*gammabb*t21*t28*t3*t31*t4*t6*w
     1         ght-1.3333333333333333d+0*gammabb*t18*t22*t23*t28*t3*wght
     2         +t22*t51*t6+t15*t45*t6+t15*t18*t46+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.5d-1*gammabb*t3*t31*t4*t52*t6*wght
     1         +t22*t6*t61+t15*t55*t6-5.0d-1*t22*t28*t29*t6+Cmat2(iq,D2_
     2         GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d+1
     1         *gammabb*t12*t20*t6*t7*(-1.8251851851851847d+2*t16*t66*t8
     2         /t5**1.1333333333333333d+1+5.866666666666666d+1*rhob*t16*
     3         t41*t8-5.688888888888889d+1*gammabb**6*t66*t67/t5**2.2d+1
     4         -8.296296296296295d+1*gammabb*t63*t66+1.5644444444444444d
     5         +2*t37*t39*t66/t5**1.6666666666666663d+1-2.13333333333333
     6         32d+1*rhob*t37*t39*t40+3.733333333333333d+1*gammabb*rhob*
     7         t33)-7.466666666666665d+2*gammabb*t13*t63*t66*t7+3.6d+1*g
     8         ammabb*t12*t6*t64*t65*t7-5.4d+1*gammabb*t6*t64*t65*t7-5.4
     9         d+1*gammabb*t12*t19*t36*t42*t6*t7+5.4d+1*gammabb*t19*t36*
     :         t42*t6*t7-1.44d+2*gammabb*rhob*t12*t18*t20*t42*t7+6.71999
     ;         9999999999d+2*gammabb*t12*t19*t20*t33*t38*t7+1.44d+2*gamm
     <         abb*rhob*t12*t18*t35*t36*t7-1.44d+2*gammabb*rhob*t18*t35*
     =         t36*t7+3.3599999999999997d+2*gammabb*rhob*t13*t33*t7-1.44
     >         d+2*gammabb*t12*t18*t19*t20*t7)*wght+2.0d+0*gammabb*t24*t
     ?         3*t43*wght+4.444444444444444d-1*gammabb*t21*t3*t34*wght)+
     @         t31*t6*(-1.5d+0*gammabb*t21*t3*t4*t43*wght-1.0d+0*gammabb
     1         *t24*t3*t30*wght)+t18*t22*(-2.6666666666666666d+0*gammabb
     2         *t23*t3*t43*wght-1.111111111111111d+1*gammabb*t21*t3*t4*w
     3         ght-2.6666666666666666d+0*rhob*t44)+8.296296296296295d+1*
     4         gammabb*rhob**4.333333333333333d+0*t15*t3*t63*wght+1.875d
     5         +0*gammabb*t21**3*t3*t4*t6*t62*wght+2.222222222222222d-1*
     6         gammabb*t21*t22*t3*t34*t6*wght+2.962962962962963d-1*gamma
     7         bb*t15*t3*t6*wght/rhob**1.6666666666666669d+0+1.866666666
     8         6666665d+1*gammabb*t21*t22*t3*t32*t33*wght-8.711111111111
     9         112d+1*gammabb*t15*t23*t3*t33*wght+6.0d+0*gammabb*t18*t23
     :         *t3*t30*t31*wght+1.4222222222222222d+1*gammabb*t15*t18*t2
     ;         4*t3*wght-1.5d+0*t21*t31*t44*t6+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d+
     1         1*gammabb*t12*t20*t6*(2.1333333333333332d+1*gammabb**5*t3
     2         8*t67/t5**2.0d+1-1.9555555555555554d+1*t38*t39*t47/t5**1.
     3         4666666666666664d+1+2.6666666666666666d+0*t39*t47*t48-2.8
     4         444444444444444d+1*t38*t39*t40*t47+3.911111111111111d+1*g
     5         ammabb*t16*t38*t41+1.2444444444444443d+1*t33*t38-2.666666
     6         6666666666d+0*t18-5.333333333333333d+0*gammabb*t16*t17)*t
     7         7+3.6d+1*gammabb*t12*t25*t35*t6*t65*t7-5.4d+1*gammabb*t25
     8         *t35*t6*t65*t7-3.6d+1*gammabb*t12*t19*t36*t49*t6*t7+3.6d+
     9         1*gammabb*t19*t36*t49*t6*t7-1.8d+1*gammabb*t12*t25*t36*t4
     :         2*t6*t7+1.8d+1*gammabb*t25*t36*t42*t6*t7+1.8d+1*t12*t20*t
     ;         42*t6*t7-1.8d+1*t12*t35*t36*t6*t7+1.8d+1*t35*t36*t6*t7-9.
     <         6d+1*gammabb*rhob*t12*t18*t20*t49*t7+2.2399999999999998d+
     =         2*gammabb*t12*t20*t25*t33*t38*t7+1.1199999999999999d+2*t1
     >         3*t33*t38*t7+9.6d+1*gammabb*rhob*t12*t18*t19*t25*t36*t7-9
     ?         .6d+1*gammabb*rhob*t18*t19*t25*t36*t7-4.8d+1*gammabb*t12*
     @         t18*t20*t25*t7-9.6d+1*rhob*t12*t18*t19*t20*t7-2.4d+1*t13*
     1         t18*t7+3.6d+1*gammabb*t12*t20*t27*t3*t42*t6-3.6d+1*gammab
     2         b*t12*t27*t3*t35*t36*t6+3.6d+1*gammabb*t27*t3*t35*t36*t6+
     3         2.2399999999999998d+2*gammabb*t13*t27*t3*t33*t38-1.92d+2*
     4         gammabb*rhob*t12*t18*t19*t20*t27*t3-4.8d+1*gammabb*t13*t1
     5         8*t27*t3)*wght+1.3333333333333333d+0*gammabb*t24*t3*t50*w
     6         ght+5.0d-1*t3*t4*t43*wght+5.0d-1*gammabb*t27*t4*t43*wght+
     7         1.3333333333333333d+0*t21*t24*t3*wght+1.3333333333333333d
     8         +0*gammabb*t21*t24*t27*wght)+t31*t6*(-1.5d+0*gammabb*t21*
     9         t3*t4*t50*wght-7.5d-1*t3*t30*t4*wght-7.5d-1*gammabb*t27*t
     :         30*t4*wght)+t18*t22*(-2.6666666666666666d+0*gammabb*t23*t
     ;         3*t50*wght-4.8888888888888893d+0*gammabb*t28*t3*t4*wght-2
     <         .6666666666666666d+0*t21*t23*t3*wght-2.6666666666666666d+
     =         0*gammabb*t21*t23*t27*wght)+t15*t18*(9.777777777777779d+0
     >         *t3*t4*wght+9.777777777777779d+0*gammabb*t27*t4*wght)+t15
     ?         *t33*(-1.2444444444444443d+1*t3*t32*wght-1.24444444444444
     @         43d+1*gammabb*t27*t32*wght)+t15*t34*t6*(-4.44444444444444
     1         4d-1*t3*wght-4.444444444444444d-1*gammabb*t27*wght)+1.875
     2         d+0*gammabb*t28*t3*t30*t4*t6*t62*wght+2.222222222222222d-
     3         1*gammabb*t22*t28*t3*t34*t6*wght+6.222222222222221d+0*gam
     4         mabb*t22*t28*t3*t32*t33*wght+4.0d+0*gammabb*t18*t21*t23*t
     5         28*t3*t31*wght-1.5d+0*t28*t31*t44*t6+Cmat3(iq,D3_RB_RB_GB
     6         B)
            Cmat3(iq,D3_RB_GBB_GBB) = t22*t6*(5.0d-1*gammabb*t3*t4*(1.8d
     1         +1*gammabb*t12*t20*t6*t7*(1.3333333333333333d+1*rhob*t39*
     2         t48*t8-8.0d+0*rhob*t37*t67/t5**1.8d+1-5.333333333333333d+
     3         0*rhob*t16*t17)+3.6d+1*gammabb*t12*t19*t56*t6*t65*t7-5.4d
     4         +1*gammabb*t19*t56*t6*t65*t7-1.8d+1*gammabb*t12*t19*t36*t
     5         58*t6*t7+1.8d+1*gammabb*t19*t36*t58*t6*t7-3.6d+1*gammabb*
     6         t12*t25*t36*t49*t6*t7+3.6d+1*gammabb*t25*t36*t49*t6*t7+3.
     7         6d+1*t12*t20*t49*t6*t7-3.6d+1*t12*t19*t25*t36*t6*t7+3.6d+
     8         1*t19*t25*t36*t6*t7-4.8d+1*gammabb*rhob*t12*t18*t20*t58*t
     9         7+4.8d+1*gammabb*rhob*t12*t18*t36*t56*t7-4.8d+1*gammabb*r
     :         hob*t18*t36*t56*t7-9.6d+1*rhob*t12*t18*t20*t25*t7+3.6d+1*
     ;         gammabb*t12*t19*t20*t59*t6+3.6d+1*gammabb*t12*t19*t20*t3*
     <         t54*t6+7.2d+1*gammabb*t12*t20*t27*t3*t49*t6-7.2d+1*gammab
     =         b*t12*t19*t25*t27*t3*t36*t6+7.2d+1*gammabb*t19*t25*t27*t3
     >         *t36*t6+7.2d+1*t12*t19*t20*t27*t3*t6-4.8d+1*gammabb*rhob*
     ?         t13*t18*t59-4.8d+1*gammabb*rhob*t13*t18*t3*t54-1.92d+2*ga
     @         mmabb*rhob*t12*t18*t20*t25*t27*t3-9.6d+1*rhob*t13*t18*t27
     1         *t3)*wght+6.666666666666666d-1*gammabb*t24*t3*t60*wght+5.
     2         0d-1*gammabb*t21*t4*t54*wght+1.0d+0*t3*t4*t50*wght+1.0d+0
     3         *gammabb*t27*t4*t50*wght+1.0d+0*t21*t27*t4*wght+6.6666666
     4         66666666d-1*t24*t28*t3*wght+6.666666666666666d-1*gammabb*
     5         t24*t27*t28*wght)+t31*t6*(-7.5d-1*gammabb*t21*t3*t4*t60*w
     6         ght-7.5d-1*gammabb*t28*t3*t4*t50*wght-7.5d-1*t21*t28*t3*t
     7         4*wght-7.5d-1*gammabb*t21*t27*t28*t4*wght)+t18*t22*(-1.33
     8         33333333333333d+0*gammabb*t23*t3*t60*wght-1.3333333333333
     9         333d+0*t23*t28*t3*wght-1.3333333333333333d+0*gammabb*t23*
     :         t27*t28*wght)+t15*t6*(-1.3333333333333333d+0*gammabb*t24*
     ;         t54*wght-2.6666666666666666d+0*t24*t27*wght)+t15*t18*(2.6
     <         666666666666666d+0*gammabb*t23*t54*wght+5.333333333333333
     =         d+0*t23*t27*wght)+1.875d+0*gammabb*t21*t3*t4*t52*t6*t62*w
     >         ght+2.0d+0*gammabb*t18*t23*t3*t31*t52*wght-1.5d+0*t28*t31
     ?         *t51*t6-5.0d-1*t22*t28*t45*t6-5.0d-1*t18*t22*t28*t46+Cmat
     @         3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t22*t6*(5.0d-1*gammabb*t3*t4*(3.6
     1         d+1*gammabb*t12*t6*t65*t69*t7-5.4d+1*gammabb*t6*t65*t69*t
     2         7+1.8d+1*gammabb*t12*t20*t6*(3*t47*t67/t5**1.6d+1-3*gamma
     3         bb*t39*t57)*t7-5.4d+1*gammabb*t12*t25*t36*t58*t6*t7+5.4d+
     4         1*gammabb*t25*t36*t58*t6*t7+5.4d+1*t12*t20*t58*t6*t7-5.4d
     5         +1*t12*t36*t56*t6*t7+5.4d+1*t36*t56*t6*t7+1.8d+1*gammabb*
     6         t13*t3*t6*t68+1.08d+2*gammabb*t12*t20*t25*t59*t6+5.4d+1*t
     7         13*t59*t6+1.08d+2*gammabb*t12*t20*t27*t3*t58*t6-1.08d+2*g
     8         ammabb*t12*t27*t3*t36*t56*t6+1.08d+2*gammabb*t27*t3*t36*t
     9         56*t6+1.08d+2*gammabb*t12*t20*t25*t3*t54*t6+5.4d+1*t13*t3
     :         *t54*t6+5.4d+1*gammabb*t13*t27*t54*t6+2.16d+2*t12*t20*t25
     ;         *t27*t3*t6)*wght+1.0d+0*t3*t4*t60*wght+1.0d+0*gammabb*t27
     <         *t4*t60*wght+5.0d-1*gammabb*t28*t4*t54*wght+1.0d+0*t27*t2
     =         8*t4*wght-t28*t55)+t15*t6*(-1.0d+0*gammabb*t4*t68*wght-3.
     >         0d+0*t4*t54*wght)+t31*t6*(-1.5d+0*gammabb*t28*t3*t4*t60*w
     ?         ght-7.5d-1*t3*t4*t52*wght-7.5d-1*gammabb*t27*t4*t52*wght)
     @         +1.875d+0*gammabb*t28**3*t3*t4*t6*t62*wght-1.5d+0*t28*t31
     1         *t6*t61-5.0d-1*t22*t29*t6*t60+7.5d-1*t29*t31*t52*t6+Cmat3
     2         (iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
