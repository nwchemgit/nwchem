C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_b97.F 
C> The nwxcm_x_b97 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,
C>    \sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{
C>    \it t_1}\right)-0.9305257363490993\,\rho_\alpha^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\alpha\alpha}^2}
C>    \over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{
C>    \it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{
C>    \it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{
C>    \it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="https://doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b97(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t10*t9*wght
            Amat(iq,D1_RA) = 1.0d+0*(-1.240700981798799d+0*rhoa**3.33333
     1         3333333333d-1*t9-9.305257363490993d-1*t10*(-1.06666666666
     2         66666d-2*gammaaa*t7*t8/rhoa**3.6666666666666664d+0+4.2666
     3         66666666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t
     4         2*t5+3.413333333333333d-7*gammaaa**3*t1*t11/rhoa**9))*wgh
     5         t+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         4.0d-3*t3*t7*t8+t5*t6*(3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t
     2         7)-1.28d-7*t1*t11*t2/rhoa**8)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 4.0d-3*gammaaa*t4+1.0d+0
            t6 = 1/t5**2
            t7 = 1/rhoa**5.333333333333333d+0
            t8 = param(3)
            t9 = 1/t5
            t10 = 4.0d-3*gammaaa*t4*t8*t9+1.6d-5*t2*t3*t6*t7+t1
            t11 = rhoa**1.3333333333333333d+0
            t12 = gammabb**2
            t13 = 1/rhob**2.6666666666666666d+0
            t14 = 4.0d-3*gammabb*t13+1.0d+0
            t15 = 1/t14**2
            t16 = 1/rhob**5.333333333333333d+0
            t17 = 1/t14
            t18 = 4.0d-3*gammabb*t13*t17*t8+1.6d-5*t12*t15*t16*t2+t1
            t19 = rhob**1.3333333333333333d+0
            t20 = 1/t5**3
            t21 = 1/rhoa**6.333333333333333d+0
            t22 = 1/t14**3
            t23 = 1/rhob**6.333333333333333d+0
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t18*t19-9.3052573634
     1         90993d-1*t10*t11)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t11*(-1.06666
     1         66666666666d-2*gammaaa*t8*t9/rhoa**3.6666666666666664d+0+
     2         4.266666666666667d-5*t21*t3*t6*t8-8.533333333333333d-5*t2
     3         *t21*t3*t6+3.413333333333333d-7*gammaaa**3*t2*t20/rhoa**9
     4         )-1.240700981798799d+0*rhoa**3.333333333333333d-1*t10)*wg
     5         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t19*(4.266666
     1         666666667d-5*t12*t15*t23*t8-1.0666666666666666d-2*gammabb
     2         *t17*t8/rhob**3.6666666666666664d+0-8.533333333333333d-5*
     3         t12*t15*t2*t23+3.413333333333333d-7*gammabb**3*t2*t22/rho
     4         b**9)-1.240700981798799d+0*rhob**3.333333333333333d-1*t18
     5         )*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t11*(
     1         4.0d-3*t4*t8*t9+t6*t7*(3.2d-5*gammaaa*t2-1.6d-5*gammaaa*t
     2         8)-1.28d-7*t2*t20*t3/rhoa**8)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t19*(
     1         t15*t16*(3.2d-5*gammabb*t2-1.6d-5*gammabb*t8)+4.0d-3*t13*
     2         t17*t8-1.28d-7*t12*t2*t22/rhob**8)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*rhoa**3.3333333333333
     1         33d-1*t9*wght-9.305257363490993d-1*t10*(-1.06666666666666
     2         66d-2*gammaaa*t7*t8/rhoa**3.6666666666666664d+0+4.2666666
     3         66666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t2*t
     4         5+3.413333333333333d-7*gammaaa**3*t1*t11/rhoa**9)*wght+Am
     5         at(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         4.0d-3*t3*t7*t8+t5*t6*(3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t
     2         7)-1.28d-7*t1*t11*t2/rhoa**8)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(4)
            t2 = gammabb**2
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 4.0d-3*gammabb*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammabb*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhob**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhob**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*rhob**3.3333333333333
     1         33d-1*t9*wght-9.305257363490993d-1*t10*(-1.06666666666666
     2         66d-2*gammabb*t7*t8/rhob**3.6666666666666664d+0+4.2666666
     3         66666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t2*t
     4         5+3.413333333333333d-7*gammabb**3*t1*t11/rhob**9)*wght+Am
     5         at(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t10*(
     1         4.0d-3*t3*t7*t8+t5*t6*(3.2d-5*gammabb*t1-1.6d-5*gammabb*t
     2         7)-1.28d-7*t1*t11*t2/rhob**8)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,
C>    \sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{
C>    \it t_1}\right)-0.9305257363490993\,\rho_\alpha^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\alpha\alpha}^2}
C>    \over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{
C>    \it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{
C>    \it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{
C>    \it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="https://doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b97_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**3
            t13 = 1/t4**3
            t14 = 1/rhoa**9
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammaaa*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhoa**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/t4**4
            t23 = 1/rhoa**10
            t24 = 1/rhoa**7.333333333333333d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t10*t9*wght
            Amat(iq,D1_RA) = 1.0d+0*(-1.240700981798799d+0*t11*t9-9.3052
     1         57363490993d-1*t10*t17)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-4.13566993932933d-1*t9/rhoa**6
     1         .666666666666666d-1-9.305257363490993d-1*t10*(3.911111111
     2         111111d-2*gammaaa*t7*t8/rhoa**4.666666666666667d+0-3.84d-
     3         4*t2*t24*t5*t7+9.102222222222223d-7*t12*t13*t23*t7+5.4044
     4         44444444444d-4*t1*t2*t24*t5-4.8924444444444437d-6*t1*t12*
     5         t13*t23+1.0922666666666667d-8*gammaaa**4*t1*t22/rhoa**1.2
     6         666666666666665d+1)-2.481401963597598d+0*t11*t17)*wght+Am
     7         at2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t10*(-1.
     1         0666666666666666d-2*t16*t7*t8+1.2800000000000003d-4*gamma
     2         aa*t15*t5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.7066666
     3         666666666d-4*gammaaa*t1*t15*t5-4.0959999999999996d-9*t1*t
     4         12*t22/rhoa**1.1666666666666665d+1+1.7066666666666666d-6*
     5         t1*t13*t14*t2)-1.240700981798799d+0*t11*(-1.6d-5*gammaaa*
     6         t5*t6*t7+3.2d-5*gammaaa*t1*t5*t6+t21+t19))*wght+Cmat2(iq,
     7         D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t10*(t5*t6*(3.2d-5*t1-3.2d-5*t7)+1.5360000000000001d
     2         -9*t1*t2*t22/rhoa**1.0666666666666666d+1+t13*t18*(-8.0d-3
     3         *t20-2.56d-7*gammaaa*t1))*wght
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
            t1 = param(2)
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 4.0d-3*gammaaa*t4+1.0d+0
            t6 = 1/t5**2
            t7 = 1/rhoa**5.333333333333333d+0
            t8 = param(3)
            t9 = 1/t5
            t10 = 4.0d-3*gammaaa*t4*t8*t9+1.6d-5*t2*t3*t6*t7+t1
            t11 = rhoa**1.3333333333333333d+0
            t12 = gammabb**2
            t13 = 1/rhob**2.6666666666666666d+0
            t14 = 4.0d-3*gammabb*t13+1.0d+0
            t15 = 1/t14**2
            t16 = 1/rhob**5.333333333333333d+0
            t17 = 1/t14
            t18 = 4.0d-3*gammabb*t13*t17*t8+1.6d-5*t12*t15*t16*t2+t1
            t19 = rhob**1.3333333333333333d+0
            t20 = rhoa**3.333333333333333d-1
            t21 = gammaaa**3
            t22 = 1/t5**3
            t23 = 1/rhoa**9
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = -1.0666666666666666d-2*gammaaa*t25*t8*t9+4.26666666666
     1         6667d-5*t24*t3*t6*t8-8.533333333333333d-5*t2*t24*t3*t6+3.
     2         413333333333333d-7*t2*t21*t22*t23
            t27 = rhob**3.333333333333333d-1
            t28 = gammabb**3
            t29 = 1/t14**3
            t30 = 1/rhob**9
            t31 = 1/rhob**6.333333333333333d+0
            t32 = 1/rhob**3.6666666666666664d+0
            t33 = -1.0666666666666666d-2*gammabb*t17*t32*t8+4.2666666666
     1         66667d-5*t12*t15*t31*t8-8.533333333333333d-5*t12*t15*t2*t
     2         31+3.413333333333333d-7*t2*t28*t29*t30
            t34 = 1/rhoa**8
            t35 = -1.28d-7*t2*t22*t3*t34
            t36 = 3.2d-5*gammaaa*t2-1.6d-5*gammaaa*t8
            t37 = 4.0d-3*t4*t8*t9
            t38 = 1/rhob**8
            t39 = -1.28d-7*t12*t2*t29*t38
            t40 = 3.2d-5*gammabb*t2-1.6d-5*gammabb*t8
            t41 = 4.0d-3*t13*t17*t8
            t42 = 1/t5**4
            t43 = 1/rhoa**10
            t44 = 1/rhoa**7.333333333333333d+0
            t45 = 1/t14**4
            t46 = 1/rhob**10
            t47 = 1/rhob**7.333333333333333d+0
            t48 = 3.2d-5*t2-3.2d-5*t8
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t18*t19-9.3052573634
     1         90993d-1*t10*t11)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t11*t26-1.240
     1         700981798799d+0*t10*t20)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t19*t33-1.240
     1         700981798799d+0*t18*t27)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t11*(
     1         t36*t6*t7+t37+t35)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t19*(
     1         t41+t15*t16*t40+t39)*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t11*(3.91
     1         1111111111111d-2*gammaaa*t8*t9/rhoa**4.666666666666667d+0
     2         -3.84d-4*t3*t44*t6*t8+9.102222222222223d-7*t21*t22*t43*t8
     3         +5.404444444444444d-4*t2*t3*t44*t6-4.8924444444444437d-6*
     4         t2*t21*t22*t43+1.0922666666666667d-8*gammaaa**4*t2*t42/rh
     5         oa**1.2666666666666665d+1)-2.481401963597598d+0*t20*t26-4
     6         .13566993932933d-1*t10/rhoa**6.666666666666666d-1)*wght+A
     7         mat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363490993d-1*t19*(-3.8
     1         4d-4*t12*t15*t47*t8+9.102222222222223d-7*t28*t29*t46*t8+3
     2         .911111111111111d-2*gammabb*t17*t8/rhob**4.66666666666666
     3         7d+0+5.404444444444444d-4*t12*t15*t2*t47-4.89244444444444
     4         37d-6*t2*t28*t29*t46+1.0922666666666667d-8*gammabb**4*t2*
     5         t45/rhob**1.2666666666666665d+1)-2.481401963597598d+0*t27
     6         *t33-4.13566993932933d-1*t18/rhob**6.666666666666666d-1)*
     7         wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t11*(-1.
     1         0666666666666666d-2*t25*t8*t9+1.2800000000000003d-4*gamma
     2         aa*t24*t6*t8-3.413333333333334d-7*t22*t23*t3*t8-1.7066666
     3         666666666d-4*gammaaa*t2*t24*t6-4.0959999999999996d-9*t2*t
     4         21*t42/rhoa**1.1666666666666665d+1+1.7066666666666666d-6*
     5         t2*t22*t23*t3)-1.240700981798799d+0*t20*(-1.6d-5*gammaaa*
     6         t6*t7*t8+3.2d-5*gammaaa*t2*t6*t7+t37+t35))*wght+Cmat2(iq,
     7         D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-9.305257363490993d-1*t19*(-1.
     1         0666666666666666d-2*t17*t32*t8+1.2800000000000003d-4*gamm
     2         abb*t15*t31*t8-3.413333333333334d-7*t12*t29*t30*t8-4.0959
     3         999999999996d-9*t2*t28*t45/rhob**1.1666666666666665d+1-1.
     4         7066666666666666d-4*gammabb*t15*t2*t31+1.7066666666666666
     5         d-6*t12*t2*t29*t30)-1.240700981798799d+0*t27*(-1.6d-5*gam
     6         mabb*t15*t16*t8+t41+t39+3.2d-5*gammabb*t15*t16*t2))*wght+
     7         Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t11*(t48*t6*t7+1.5360000000000001d-9*t2*t3*t42/rhoa*
     2         *1.0666666666666666d+1+t22*t34*(-8.0d-3*t36-2.56d-7*gamma
     3         aa*t2))*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t19*(t15*t16*t48+1.5360000000000001d-9*t12*t2*t45/rh
     2         ob**1.0666666666666666d+1+t29*t38*(-8.0d-3*t40-2.56d-7*ga
     3         mmabb*t2))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**3
            t13 = 1/t4**3
            t14 = 1/rhoa**9
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammaaa*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhoa**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/t4**4
            t23 = 1/rhoa**10
            t24 = 1/rhoa**7.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*t11*t9*wght-9.3052573
     1         63490993d-1*t10*t17*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Amat2(iq,D2_RA_RA) = -4.13566993932933d-1*t9*wght/rhoa**6.66
     1         6666666666666d-1-9.305257363490993d-1*t10*(3.911111111111
     2         111d-2*gammaaa*t7*t8/rhoa**4.666666666666667d+0-3.84d-4*t
     3         2*t24*t5*t7+9.102222222222223d-7*t12*t13*t23*t7+5.4044444
     4         44444444d-4*t1*t2*t24*t5-4.8924444444444437d-6*t1*t12*t13
     5         *t23+1.0922666666666667d-8*gammaaa**4*t1*t22/rhoa**1.2666
     6         666666666665d+1)*wght-2.481401963597598d+0*t11*t17*wght+A
     7         mat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490993d-1*t10*(-1.06666666
     1         66666666d-2*t16*t7*t8+1.2800000000000003d-4*gammaaa*t15*t
     2         5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.706666666666666
     3         6d-4*gammaaa*t1*t15*t5-4.0959999999999996d-9*t1*t12*t22/r
     4         hoa**1.1666666666666665d+1+1.7066666666666666d-6*t1*t13*t
     5         14*t2)*wght-1.240700981798799d+0*t11*(-1.6d-5*gammaaa*t5*
     6         t6*t7+3.2d-5*gammaaa*t1*t5*t6+t21+t19)*wght+Cmat2(iq,D2_R
     7         A_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t10*(t5*t6*(3.2d-5*t1-3.2d-5*t7)+1.5360000000000001d
     2         -9*t1*t2*t22/rhoa**1.0666666666666666d+1+t13*t18*(-8.0d-3
     3         *t20-2.56d-7*gammaaa*t1))*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(4)
            t2 = gammabb**2
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 4.0d-3*gammabb*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammabb*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhob**1.3333333333333333d+0
            t11 = rhob**3.333333333333333d-1
            t12 = gammabb**3
            t13 = 1/t4**3
            t14 = 1/rhob**9
            t15 = 1/rhob**6.333333333333333d+0
            t16 = 1/rhob**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammabb*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhob**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammabb*t1-1.6d-5*gammabb*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/t4**4
            t23 = 1/rhob**10
            t24 = 1/rhob**7.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*t11*t9*wght-9.3052573
     1         63490993d-1*t10*t17*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Amat2(iq,D2_RB_RB) = -4.13566993932933d-1*t9*wght/rhob**6.66
     1         6666666666666d-1-9.305257363490993d-1*t10*(3.911111111111
     2         111d-2*gammabb*t7*t8/rhob**4.666666666666667d+0-3.84d-4*t
     3         2*t24*t5*t7+9.102222222222223d-7*t12*t13*t23*t7+5.4044444
     4         44444444d-4*t1*t2*t24*t5-4.8924444444444437d-6*t1*t12*t13
     5         *t23+1.0922666666666667d-8*gammabb**4*t1*t22/rhob**1.2666
     6         666666666665d+1)*wght-2.481401963597598d+0*t11*t17*wght+A
     7         mat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490993d-1*t10*(-1.06666666
     1         66666666d-2*t16*t7*t8+1.2800000000000003d-4*gammabb*t15*t
     2         5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.706666666666666
     3         6d-4*gammabb*t1*t15*t5-4.0959999999999996d-9*t1*t12*t22/r
     4         hob**1.1666666666666665d+1+1.7066666666666666d-6*t1*t13*t
     5         14*t2)*wght-1.240700981798799d+0*t11*(-1.6d-5*gammabb*t5*
     6         t6*t7+3.2d-5*gammabb*t1*t5*t6+t21+t19)*wght+Cmat2(iq,D2_R
     7         B_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t10*(t5*t6*(3.2d-5*t1-3.2d-5*t7)+1.5360000000000001d
     2         -9*t1*t2*t22/rhob**1.0666666666666666d+1+t13*t18*(-8.0d-3
     3         *t20-2.56d-7*gammabb*t1))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-0.9305257363490993\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,
C>    \sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{
C>    \it t_1}\right)-0.9305257363490993\,\rho_\alpha^{{{4}
C>    \over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}
C>    \,{\it t_2}\,\sigma_{\alpha\alpha}^2}
C>    \over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{
C>    \it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490993\,\rho_s^{{{4}\over{3}}}\,
C>    \left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{
C>    \it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{
C>    \it t_1}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="https://doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b97_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**3
            t13 = 1/t4**3
            t14 = 1/rhoa**9
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammaaa*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhoa**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/rhoa**6.666666666666666d-1
            t23 = gammaaa**4
            t24 = 1/t4**4
            t25 = 1/rhoa**1.2666666666666665d+1
            t26 = 1/rhoa**10
            t27 = 1/rhoa**7.333333333333333d+0
            t28 = 1/rhoa**4.666666666666667d+0
            t29 = 3.911111111111111d-2*gammaaa*t28*t7*t8-3.84d-4*t2*t27*
     1         t5*t7+9.102222222222223d-7*t12*t13*t26*t7+5.4044444444444
     2         44d-4*t1*t2*t27*t5-4.8924444444444437d-6*t1*t12*t13*t26+1
     3         .0922666666666667d-8*t1*t23*t24*t25
            t30 = -1.6d-5*gammaaa*t5*t6*t7+3.2d-5*gammaaa*t1*t5*t6+t21+t
     1         19
            t31 = 1/rhoa**1.1666666666666665d+1
            t32 = -1.0666666666666666d-2*t16*t7*t8+1.2800000000000003d-4
     1         *gammaaa*t15*t5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.7
     2         066666666666666d-4*gammaaa*t1*t15*t5-4.0959999999999996d-
     3         9*t1*t12*t24*t31+1.7066666666666666d-6*t1*t13*t14*t2
            t33 = 1/rhoa**1.0666666666666666d+1
            t34 = 1.5360000000000001d-9*t1*t2*t24*t33
            t35 = -8.0d-3*t20-2.56d-7*gammaaa*t1
            t36 = 3.2d-5*t1
            t37 = t36-3.2d-5*t7
            t38 = 1/t4**5
            t39 = 1/rhoa**1.3666666666666666d+1
            t40 = 1/rhoa**11
            t41 = 1/rhoa**8.333333333333334d+0
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t10*t9*wght
            Amat(iq,D1_RA) = 1.0d+0*(-1.240700981798799d+0*t11*t9-9.3052
     1         57363490993d-1*t10*t17)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-4.13566993932933d-1*t22*t9-9.3
     1         05257363490993d-1*t10*t29-2.481401963597598d+0*t11*t17)*w
     2         ght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t10*t32-
     1         1.240700981798799d+0*t11*t30)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t10*(t37*t5*t6+t13*t18*t35+t34)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(2.7571132928862196d-1*t9/rho
     1         a**1.6666666666666669d+0-9.305257363490993d-1*t10*(-1.825
     2         1851851851852d-1*gammaaa*t7*t8/rhoa**5.666666666666667d+0
     3         +3.233185185185185d-3*t2*t41*t5*t7-1.7294222222222225d-5*
     4         t12*t13*t40*t7+2.9127111111111115d-8*t23*t24*t39*t7-3.963
     5         259259259259d-3*t1*t2*t41*t5+6.045392592592592d-5*t1*t12*
     6         t13*t40-2.94912d-7*t1*t23*t24*t39+4.660337777777778d-10*g
     7         ammaaa**5*t1*t38/rhoa**1.6333333333333333d+1)-3.722102945
     8         396397d+0*t11*t29-1.240700981798799d+0*t17*t22)*wght+Amat
     9         3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t10*(
     1         3.911111111111111d-2*t28*t7*t8-9.244444444444444d-4*gamma
     2         aa*t27*t5*t7+5.802666666666667d-6*t13*t2*t26*t7-1.0922666
     3         666666667d-8*t12*t24*t25*t7+1.0808888888888889d-3*gammaaa
     4         *t1*t27*t5-1.7476266666666668d-10*t1*t23*t38/rhoa**1.5333
     5         333333333335d+1-1.9000888888888887d-5*t1*t13*t2*t26+1.023
     6         9999999999999d-7*t1*t12*t24*t25)-2.481401963597598d+0*t11
     7         *t32-4.13566993932933d-1*t22*t30)*wght+Cmat3(iq,D3_RA_RA_
     8         GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(-1.240700981798799d+0*t11*
     1         (-3.2d-5*t5*t6*t7+1.28d-7*gammaaa*t13*t18*t7+3.2d-5*t1*t5
     2         *t6+t34-5.12d-7*gammaaa*t1*t13*t18)-9.305257363490993d-1*
     3         t10*(1.706666666666667d-4*t15*t5*t7+4.096000000000001d-9*
     4         t2*t24*t31*t7-1.7066666666666672d-6*gammaaa*t13*t14*t7-1.
     5         7066666666666666d-4*t1*t15*t5+6.5536d-11*t1*t12*t38/rhoa*
     6         *1.4333333333333334d+1-3.2767999999999997d-8*t1*t2*t24*t3
     7         1+4.778666666666666d-6*gammaaa*t1*t13*t14))*wght+Cmat3(iq
     8         ,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t10*(t13*t18*(-8.0d-3*(t36-1.6d-5*t7)-8.0d-3
     2         *t37-2.56d-7*t1)-2.4576d-11*t1*t2*t38/rhoa**1.33333333333
     3         33333d+1+t24*t33*(3.0720000000000003d-9*gammaaa*t1-1.2d-2
     4         *t35))*wght
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
            t1 = param(2)
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 4.0d-3*gammaaa*t4+1.0d+0
            t6 = 1/t5**2
            t7 = 1/rhoa**5.333333333333333d+0
            t8 = param(3)
            t9 = 1/t5
            t10 = 4.0d-3*gammaaa*t4*t8*t9+1.6d-5*t2*t3*t6*t7+t1
            t11 = rhoa**1.3333333333333333d+0
            t12 = gammabb**2
            t13 = 1/rhob**2.6666666666666666d+0
            t14 = 4.0d-3*gammabb*t13+1.0d+0
            t15 = 1/t14**2
            t16 = 1/rhob**5.333333333333333d+0
            t17 = 1/t14
            t18 = 4.0d-3*gammabb*t13*t17*t8+1.6d-5*t12*t15*t16*t2+t1
            t19 = rhob**1.3333333333333333d+0
            t20 = rhoa**3.333333333333333d-1
            t21 = gammaaa**3
            t22 = 1/t5**3
            t23 = 1/rhoa**9
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = -1.0666666666666666d-2*gammaaa*t25*t8*t9+4.26666666666
     1         6667d-5*t24*t3*t6*t8-8.533333333333333d-5*t2*t24*t3*t6+3.
     2         413333333333333d-7*t2*t21*t22*t23
            t27 = rhob**3.333333333333333d-1
            t28 = gammabb**3
            t29 = 1/t14**3
            t30 = 1/rhob**9
            t31 = 1/rhob**6.333333333333333d+0
            t32 = 1/rhob**3.6666666666666664d+0
            t33 = -1.0666666666666666d-2*gammabb*t17*t32*t8+4.2666666666
     1         66667d-5*t12*t15*t31*t8-8.533333333333333d-5*t12*t15*t2*t
     2         31+3.413333333333333d-7*t2*t28*t29*t30
            t34 = 1/rhoa**8
            t35 = -1.28d-7*t2*t22*t3*t34
            t36 = 3.2d-5*gammaaa*t2-1.6d-5*gammaaa*t8
            t37 = 4.0d-3*t4*t8*t9
            t38 = 1/rhob**8
            t39 = -1.28d-7*t12*t2*t29*t38
            t40 = 3.2d-5*gammabb*t2-1.6d-5*gammabb*t8
            t41 = 4.0d-3*t13*t17*t8
            t42 = 1/rhoa**6.666666666666666d-1
            t43 = gammaaa**4
            t44 = 1/t5**4
            t45 = 1/rhoa**1.2666666666666665d+1
            t46 = 1/rhoa**10
            t47 = 1/rhoa**7.333333333333333d+0
            t48 = 1/rhoa**4.666666666666667d+0
            t49 = 3.911111111111111d-2*gammaaa*t48*t8*t9-3.84d-4*t3*t47*
     1         t6*t8+9.102222222222223d-7*t21*t22*t46*t8+5.4044444444444
     2         44d-4*t2*t3*t47*t6-4.8924444444444437d-6*t2*t21*t22*t46+1
     3         .0922666666666667d-8*t2*t43*t44*t45
            t50 = 1/rhob**6.666666666666666d-1
            t51 = gammabb**4
            t52 = 1/t14**4
            t53 = 1/rhob**1.2666666666666665d+1
            t54 = 1/rhob**10
            t55 = 1/rhob**7.333333333333333d+0
            t56 = 1/rhob**4.666666666666667d+0
            t57 = 3.911111111111111d-2*gammabb*t17*t56*t8-3.84d-4*t12*t1
     1         5*t55*t8+9.102222222222223d-7*t28*t29*t54*t8+5.4044444444
     2         44444d-4*t12*t15*t2*t55-4.8924444444444437d-6*t2*t28*t29*
     3         t54+1.0922666666666667d-8*t2*t51*t52*t53
            t58 = -1.6d-5*gammaaa*t6*t7*t8+3.2d-5*gammaaa*t2*t6*t7+t37+t
     1         35
            t59 = 1/rhoa**1.1666666666666665d+1
            t60 = -1.0666666666666666d-2*t25*t8*t9+1.2800000000000003d-4
     1         *gammaaa*t24*t6*t8-3.413333333333334d-7*t22*t23*t3*t8-1.7
     2         066666666666666d-4*gammaaa*t2*t24*t6-4.0959999999999996d-
     3         9*t2*t21*t44*t59+1.7066666666666666d-6*t2*t22*t23*t3
            t61 = -1.6d-5*gammabb*t15*t16*t8+t41+t39+3.2d-5*gammabb*t15*
     1         t16*t2
            t62 = 1/rhob**1.1666666666666665d+1
            t63 = -1.0666666666666666d-2*t17*t32*t8+1.2800000000000003d-
     1         4*gammabb*t15*t31*t8-3.413333333333334d-7*t12*t29*t30*t8-
     2         4.0959999999999996d-9*t2*t28*t52*t62-1.7066666666666666d-
     3         4*gammabb*t15*t2*t31+1.7066666666666666d-6*t12*t2*t29*t30
            t64 = 1/rhoa**1.0666666666666666d+1
            t65 = 1.5360000000000001d-9*t2*t3*t44*t64
            t66 = -8.0d-3*t36-2.56d-7*gammaaa*t2
            t67 = 3.2d-5*t2
            t68 = t67-3.2d-5*t8
            t69 = 1/rhob**1.0666666666666666d+1
            t70 = 1.5360000000000001d-9*t12*t2*t52*t69
            t71 = -8.0d-3*t40-2.56d-7*gammabb*t2
            t72 = 1/t5**5
            t73 = 1/rhoa**1.3666666666666666d+1
            t74 = 1/rhoa**11
            t75 = 1/rhoa**8.333333333333334d+0
            t76 = 1/t14**5
            t77 = 1/rhob**1.3666666666666666d+1
            t78 = 1/rhob**11
            t79 = 1/rhob**8.333333333333334d+0
            t80 = -8.0d-3*(t67-1.6d-5*t8)-8.0d-3*t68-2.56d-7*t2
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t18*t19-9.3052573634
     1         90993d-1*t10*t11)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t11*t26-1.240
     1         700981798799d+0*t10*t20)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t19*t33-1.240
     1         700981798799d+0*t18*t27)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t11*(
     1         t36*t6*t7+t37+t35)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t19*(
     1         t41+t15*t16*t40+t39)*wght
            Amat2(iq,D2_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t11*t49-4
     1         .13566993932933d-1*t10*t42-2.481401963597598d+0*t20*t26)*
     2         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 1.0d+0*(-9.305257363490993d-1*t19*t57-4
     1         .13566993932933d-1*t18*t50-2.481401963597598d+0*t27*t33)*
     2         wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t11*t60-
     1         1.240700981798799d+0*t20*t58)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*(-9.305257363490993d-1*t19*t63-
     1         1.240700981798799d+0*t27*t61)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t11*(t6*t68*t7+t22*t34*t66+t65)*wght
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t19*(t29*t38*t71+t70+t15*t16*t68)*wght
            Amat3(iq,D3_RA_RA_RA) = 1.0d+0*(-9.305257363490993d-1*t11*(-
     1         1.8251851851851852d-1*gammaaa*t8*t9/rhoa**5.6666666666666
     2         67d+0+3.233185185185185d-3*t3*t6*t75*t8-1.729422222222222
     3         5d-5*t21*t22*t74*t8+2.9127111111111115d-8*t43*t44*t73*t8-
     4         3.963259259259259d-3*t2*t3*t6*t75+6.045392592592592d-5*t2
     5         *t21*t22*t74-2.94912d-7*t2*t43*t44*t73+4.660337777777778d
     6         -10*gammaaa**5*t2*t72/rhoa**1.6333333333333333d+1)-3.7221
     7         02945396397d+0*t20*t49-1.240700981798799d+0*t26*t42+2.757
     8         1132928862196d-1*t10/rhoa**1.6666666666666669d+0)*wght+Am
     9         at3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = 1.0d+0*(-9.305257363490993d-1*t19*(3
     1         .233185185185185d-3*t12*t15*t79*t8-1.7294222222222225d-5*
     2         t28*t29*t78*t8+2.9127111111111115d-8*t51*t52*t77*t8-1.825
     3         1851851851852d-1*gammabb*t17*t8/rhob**5.666666666666667d+
     4         0-3.963259259259259d-3*t12*t15*t2*t79+6.045392592592592d-
     5         5*t2*t28*t29*t78-2.94912d-7*t2*t51*t52*t77+4.660337777777
     6         778d-10*gammabb**5*t2*t76/rhob**1.6333333333333333d+1)-3.
     7         722102945396397d+0*t27*t57-1.240700981798799d+0*t33*t50+2
     8         .7571132928862196d-1*t18/rhob**1.6666666666666669d+0)*wgh
     9         t+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 1.0d+0*(-9.305257363490993d-1*t11*(
     1         3.911111111111111d-2*t48*t8*t9-9.244444444444444d-4*gamma
     2         aa*t47*t6*t8+5.802666666666667d-6*t22*t3*t46*t8-1.0922666
     3         666666667d-8*t21*t44*t45*t8-1.7476266666666668d-10*t2*t43
     4         *t72/rhoa**1.5333333333333335d+1+1.0808888888888889d-3*ga
     5         mmaaa*t2*t47*t6-1.9000888888888887d-5*t2*t22*t3*t46+1.023
     6         9999999999999d-7*t2*t21*t44*t45)-2.481401963597598d+0*t20
     7         *t60-4.13566993932933d-1*t42*t58)*wght+Cmat3(iq,D3_RA_RA_
     8         GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = 1.0d+0*(-9.305257363490993d-1*t19*(
     1         3.911111111111111d-2*t17*t56*t8-9.244444444444444d-4*gamm
     2         abb*t15*t55*t8+5.802666666666667d-6*t12*t29*t54*t8-1.0922
     3         666666666667d-8*t28*t52*t53*t8-1.7476266666666668d-10*t2*
     4         t51*t76/rhob**1.5333333333333335d+1+1.0808888888888889d-3
     5         *gammabb*t15*t2*t55-1.9000888888888887d-5*t12*t2*t29*t54+
     6         1.0239999999999999d-7*t2*t28*t52*t53)-2.481401963597598d+
     7         0*t27*t63-4.13566993932933d-1*t50*t61)*wght+Cmat3(iq,D3_R
     8         B_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = 1.0d+0*(-1.240700981798799d+0*t20*
     1         (-3.2d-5*t6*t7*t8+1.28d-7*gammaaa*t22*t34*t8+3.2d-5*t2*t6
     2         *t7+t65-5.12d-7*gammaaa*t2*t22*t34)-9.305257363490993d-1*
     3         t11*(1.706666666666667d-4*t24*t6*t8+4.096000000000001d-9*
     4         t3*t44*t59*t8-1.7066666666666672d-6*gammaaa*t22*t23*t8+6.
     5         5536d-11*t2*t21*t72/rhoa**1.4333333333333334d+1-1.7066666
     6         666666666d-4*t2*t24*t6-3.2767999999999997d-8*t2*t3*t44*t5
     7         9+4.778666666666666d-6*gammaaa*t2*t22*t23))*wght+Cmat3(iq
     8         ,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = 1.0d+0*(-9.305257363490993d-1*t19*
     1         (4.096000000000001d-9*t12*t52*t62*t8+1.706666666666667d-4
     2         *t15*t31*t8-1.7066666666666672d-6*gammabb*t29*t30*t8+6.55
     3         36d-11*t2*t28*t76/rhob**1.4333333333333334d+1-3.276799999
     4         9999997d-8*t12*t2*t52*t62-1.7066666666666666d-4*t15*t2*t3
     5         1+4.778666666666666d-6*gammabb*t2*t29*t30)-1.240700981798
     6         799d+0*t27*(1.28d-7*gammabb*t29*t38*t8-3.2d-5*t15*t16*t8+
     7         t70-5.12d-7*gammabb*t2*t29*t38+3.2d-5*t15*t16*t2))*wght+C
     8         mat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t11*(t22*t34*t80-2.4576d-11*t2*t3*t72/rhoa**
     2         1.3333333333333333d+1+t44*t64*(3.0720000000000003d-9*gamm
     3         aaa*t2-1.2d-2*t66))*wght
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490993d-1*t19*(t29*t38*t80-2.4576d-11*t12*t2*t76/rhob*
     2         *1.3333333333333333d+1+t52*t69*(3.0720000000000003d-9*gam
     3         mabb*t2-1.2d-2*t71))*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = gammaaa**3
            t13 = 1/t4**3
            t14 = 1/rhoa**9
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammaaa*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhoa**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/rhoa**6.666666666666666d-1
            t23 = gammaaa**4
            t24 = 1/t4**4
            t25 = 1/rhoa**1.2666666666666665d+1
            t26 = 1/rhoa**10
            t27 = 1/rhoa**7.333333333333333d+0
            t28 = 1/rhoa**4.666666666666667d+0
            t29 = 3.911111111111111d-2*gammaaa*t28*t7*t8-3.84d-4*t2*t27*
     1         t5*t7+9.102222222222223d-7*t12*t13*t26*t7+5.4044444444444
     2         44d-4*t1*t2*t27*t5-4.8924444444444437d-6*t1*t12*t13*t26+1
     3         .0922666666666667d-8*t1*t23*t24*t25
            t30 = -1.6d-5*gammaaa*t5*t6*t7+3.2d-5*gammaaa*t1*t5*t6+t21+t
     1         19
            t31 = 1/rhoa**1.1666666666666665d+1
            t32 = -1.0666666666666666d-2*t16*t7*t8+1.2800000000000003d-4
     1         *gammaaa*t15*t5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.7
     2         066666666666666d-4*gammaaa*t1*t15*t5-4.0959999999999996d-
     3         9*t1*t12*t24*t31+1.7066666666666666d-6*t1*t13*t14*t2
            t33 = 1/rhoa**1.0666666666666666d+1
            t34 = 1.5360000000000001d-9*t1*t2*t24*t33
            t35 = -8.0d-3*t20-2.56d-7*gammaaa*t1
            t36 = 3.2d-5*t1
            t37 = t36-3.2d-5*t7
            t38 = 1/t4**5
            t39 = 1/rhoa**1.3666666666666666d+1
            t40 = 1/rhoa**11
            t41 = 1/rhoa**8.333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*t11*t9*wght-9.3052573
     1         63490993d-1*t10*t17*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Amat2(iq,D2_RA_RA) = -4.13566993932933d-1*t22*t9*wght-9.3052
     1         57363490993d-1*t10*t29*wght-2.481401963597598d+0*t11*t17*
     2         wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -9.305257363490993d-1*t10*t32*wght-1.2
     1         40700981798799d+0*t11*t30*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)-9.30525736349099
     1         3d-1*t10*(t37*t5*t6+t13*t18*t35+t34)*wght
            Amat3(iq,D3_RA_RA_RA) = 2.7571132928862196d-1*t9*wght/rhoa**
     1         1.6666666666666669d+0-9.305257363490993d-1*t10*(-1.825185
     2         1851851852d-1*gammaaa*t7*t8/rhoa**5.666666666666667d+0+3.
     3         233185185185185d-3*t2*t41*t5*t7-1.7294222222222225d-5*t12
     4         *t13*t40*t7+2.9127111111111115d-8*t23*t24*t39*t7-3.963259
     5         259259259d-3*t1*t2*t41*t5+6.045392592592592d-5*t1*t12*t13
     6         *t40-2.94912d-7*t1*t23*t24*t39+4.660337777777778d-10*gamm
     7         aaa**5*t1*t38/rhoa**1.6333333333333333d+1)*wght-3.7221029
     8         45396397d+0*t11*t29*wght-1.240700981798799d+0*t17*t22*wgh
     9         t+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -9.305257363490993d-1*t10*(3.911111
     1         111111111d-2*t28*t7*t8-9.244444444444444d-4*gammaaa*t27*t
     2         5*t7+5.802666666666667d-6*t13*t2*t26*t7-1.092266666666666
     3         7d-8*t12*t24*t25*t7+1.0808888888888889d-3*gammaaa*t1*t27*
     4         t5-1.7476266666666668d-10*t1*t23*t38/rhoa**1.533333333333
     5         3335d+1-1.9000888888888887d-5*t1*t13*t2*t26+1.02399999999
     6         99999d-7*t1*t12*t24*t25)*wght-2.481401963597598d+0*t11*t3
     7         2*wght-4.13566993932933d-1*t22*t30*wght+Cmat3(iq,D3_RA_RA
     8         _GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -1.240700981798799d+0*t11*(-3.2d-5
     1         *t5*t6*t7+1.28d-7*gammaaa*t13*t18*t7+3.2d-5*t1*t5*t6+t34-
     2         5.12d-7*gammaaa*t1*t13*t18)*wght-9.305257363490993d-1*t10
     3         *(1.706666666666667d-4*t15*t5*t7+4.096000000000001d-9*t2*
     4         t24*t31*t7-1.7066666666666672d-6*gammaaa*t13*t14*t7-1.706
     5         6666666666666d-4*t1*t15*t5+6.5536d-11*t1*t12*t38/rhoa**1.
     6         4333333333333334d+1-3.2767999999999997d-8*t1*t2*t24*t31+4
     7         .778666666666666d-6*gammaaa*t1*t13*t14)*wght+Cmat3(iq,D3_
     8         RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-9.305257
     1         363490993d-1*t10*(t13*t18*(-8.0d-3*(t36-1.6d-5*t7)-8.0d-3
     2         *t37-2.56d-7*t1)-2.4576d-11*t1*t2*t38/rhoa**1.33333333333
     3         33333d+1+t24*t33*(3.0720000000000003d-9*gammaaa*t1-1.2d-2
     4         *t35))*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(4)
            t2 = gammabb**2
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 4.0d-3*gammabb*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammabb*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhob**1.3333333333333333d+0
            t11 = rhob**3.333333333333333d-1
            t12 = gammabb**3
            t13 = 1/t4**3
            t14 = 1/rhob**9
            t15 = 1/rhob**6.333333333333333d+0
            t16 = 1/rhob**3.6666666666666664d+0
            t17 = -1.0666666666666666d-2*gammabb*t16*t7*t8+4.26666666666
     1         6667d-5*t15*t2*t5*t7-8.533333333333333d-5*t1*t15*t2*t5+3.
     2         413333333333333d-7*t1*t12*t13*t14
            t18 = 1/rhob**8
            t19 = -1.28d-7*t1*t13*t18*t2
            t20 = 3.2d-5*gammabb*t1-1.6d-5*gammabb*t7
            t21 = 4.0d-3*t3*t7*t8
            t22 = 1/rhob**6.666666666666666d-1
            t23 = gammabb**4
            t24 = 1/t4**4
            t25 = 1/rhob**1.2666666666666665d+1
            t26 = 1/rhob**10
            t27 = 1/rhob**7.333333333333333d+0
            t28 = 1/rhob**4.666666666666667d+0
            t29 = 3.911111111111111d-2*gammabb*t28*t7*t8-3.84d-4*t2*t27*
     1         t5*t7+9.102222222222223d-7*t12*t13*t26*t7+5.4044444444444
     2         44d-4*t1*t2*t27*t5-4.8924444444444437d-6*t1*t12*t13*t26+1
     3         .0922666666666667d-8*t1*t23*t24*t25
            t30 = -1.6d-5*gammabb*t5*t6*t7+3.2d-5*gammabb*t1*t5*t6+t21+t
     1         19
            t31 = 1/rhob**1.1666666666666665d+1
            t32 = -1.0666666666666666d-2*t16*t7*t8+1.2800000000000003d-4
     1         *gammabb*t15*t5*t7-3.413333333333334d-7*t13*t14*t2*t7-1.7
     2         066666666666666d-4*gammabb*t1*t15*t5-4.0959999999999996d-
     3         9*t1*t12*t24*t31+1.7066666666666666d-6*t1*t13*t14*t2
            t33 = 1/rhob**1.0666666666666666d+1
            t34 = 1.5360000000000001d-9*t1*t2*t24*t33
            t35 = -8.0d-3*t20-2.56d-7*gammabb*t1
            t36 = 3.2d-5*t1
            t37 = t36-3.2d-5*t7
            t38 = 1/t4**5
            t39 = 1/rhob**1.3666666666666666d+1
            t40 = 1/rhob**11
            t41 = 1/rhob**8.333333333333334d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*t11*t9*wght-9.3052573
     1         63490993d-1*t10*t17*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t10*(
     1         t20*t5*t6+t21+t19)*wght
            Amat2(iq,D2_RB_RB) = -4.13566993932933d-1*t22*t9*wght-9.3052
     1         57363490993d-1*t10*t29*wght-2.481401963597598d+0*t11*t17*
     2         wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -9.305257363490993d-1*t10*t32*wght-1.2
     1         40700981798799d+0*t11*t30*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)-9.30525736349099
     1         3d-1*t10*(t37*t5*t6+t13*t18*t35+t34)*wght
            Amat3(iq,D3_RB_RB_RB) = 2.7571132928862196d-1*t9*wght/rhob**
     1         1.6666666666666669d+0-9.305257363490993d-1*t10*(-1.825185
     2         1851851852d-1*gammabb*t7*t8/rhob**5.666666666666667d+0+3.
     3         233185185185185d-3*t2*t41*t5*t7-1.7294222222222225d-5*t12
     4         *t13*t40*t7+2.9127111111111115d-8*t23*t24*t39*t7-3.963259
     5         259259259d-3*t1*t2*t41*t5+6.045392592592592d-5*t1*t12*t13
     6         *t40-2.94912d-7*t1*t23*t24*t39+4.660337777777778d-10*gamm
     7         abb**5*t1*t38/rhob**1.6333333333333333d+1)*wght-3.7221029
     8         45396397d+0*t11*t29*wght-1.240700981798799d+0*t17*t22*wgh
     9         t+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -9.305257363490993d-1*t10*(3.911111
     1         111111111d-2*t28*t7*t8-9.244444444444444d-4*gammabb*t27*t
     2         5*t7+5.802666666666667d-6*t13*t2*t26*t7-1.092266666666666
     3         7d-8*t12*t24*t25*t7+1.0808888888888889d-3*gammabb*t1*t27*
     4         t5-1.7476266666666668d-10*t1*t23*t38/rhob**1.533333333333
     5         3335d+1-1.9000888888888887d-5*t1*t13*t2*t26+1.02399999999
     6         99999d-7*t1*t12*t24*t25)*wght-2.481401963597598d+0*t11*t3
     7         2*wght-4.13566993932933d-1*t22*t30*wght+Cmat3(iq,D3_RB_RB
     8         _GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -1.240700981798799d+0*t11*(-3.2d-5
     1         *t5*t6*t7+1.28d-7*gammabb*t13*t18*t7+3.2d-5*t1*t5*t6+t34-
     2         5.12d-7*gammabb*t1*t13*t18)*wght-9.305257363490993d-1*t10
     3         *(1.706666666666667d-4*t15*t5*t7+4.096000000000001d-9*t2*
     4         t24*t31*t7-1.7066666666666672d-6*gammabb*t13*t14*t7-1.706
     5         6666666666666d-4*t1*t15*t5+6.5536d-11*t1*t12*t38/rhob**1.
     6         4333333333333334d+1-3.2767999999999997d-8*t1*t2*t24*t31+4
     7         .778666666666666d-6*gammabb*t1*t13*t14)*wght+Cmat3(iq,D3_
     8         RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-9.305257
     1         363490993d-1*t10*(t13*t18*(-8.0d-3*(t36-1.6d-5*t7)-8.0d-3
     2         *t37-2.56d-7*t1)-2.4576d-11*t1*t2*t38/rhob**1.33333333333
     3         33333d+1+t24*t33*(3.0720000000000003d-9*gammabb*t1-1.2d-2
     4         *t35))*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
