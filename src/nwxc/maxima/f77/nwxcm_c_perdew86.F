C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_perdew86.F 
C> The nwxcm_c_perdew86 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_perdew86 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_5}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{0.01443307452126544\,{\it t_6}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_4}+0.002568}
C>    \over{0.01763993811759022\,{\it t_4}\,{\it t_6}
C>    +5.411317332115466\,{\it t_6}+0.1816419932959077\,{\it t_4}
C>    +1.0}}+0.001667\\\\ 
C>   {\it t_8} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_{10}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{12}} &=& {{0.01443307452126544\,{\it t_{11}}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_9}+0.002568}
C>    \over{0.01763993811759022\,{\it t_9}\,{\it t_{11}}
C>    +5.411317332115466\,{\it t_{11}}+0.1816419932959077\,{\it t_9}
C>    +1.0}}+0.001667\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_7}\,{\it t_6}
C>    \,{\it t_8}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{
C>    \it t_3}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}
C>    +0.3149802624737183\,\left(1.0-1.0\,{\it t_3}\,{
C>    \it t_2}\right)^{{{5}\over{3}}}}
C>    \,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{
C>    \it t_8}}}\over{{\it t_7}\,{\it t_5}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.7937005259840999\,{\it t_{12}}\,{\it t_{11}}
C>    \,\sigma_{ss}}\over{\rho_s\,e^{{{8.129082500000001 \times 10^{
C>    -4}\,\sqrt{\sigma_{ss}}}\over{{\it t_{12}}\,{
C>    \it t_{10}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew86(param,tol_rho,ipol,nq,wght, 
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
            t1 = 3.968502629920499d-1
            t2 = 1/rhoa
            t3 = gammaaa**5.0d-1
            t4 = 6.299605249474366d-1
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = 7.937005259840998d-1
            t7 = abs(rhoa)
            t8 = 1/t7**3.333333333333333d-1
            t9 = 1.443307452126544d-2*t6*t8+2.8435438314903855d-6*t4*t5+
     1         2.568d-3
            t10 = 5.411317332115466d+0*t6*t8+8.819969058795109d-3*t5*t8+
     1         1.816419932959077d-1*t4*t5+1.0d+0
            t11 = 1/t10
            t12 = t11*t9+1.6670000000000001d-3
            t13 = 1/t12
            t14 = 1/t7**1.166666666666667d+0
            t15 = exp(-7.242189178907065d-4*t13*t14*t3)
            t16 = 1/t7
            t17 = 3.149802624737183d-1
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa**3.333333333333333d-1
            t20 = t11*(-4.8110248404218126d-3*t1*t16*t19-1.8956958876602
     1         57d-6*t17*t18)-(-2.939989686265036d-3*t18*t8-1.4699948431
     2         32518d-3*t16*t2-1.8037724440384886d+0*t1*t16*t19-1.210946
     3         6219727179d-1*t17*t18)*t9/t10**2
            t21 = 1/t7**1.5000000000000002d+0
            fnc(iq) = 4.000000000000001d+0*gammaaa*t1*t12*t2*exp(-7.2421
     1         89178907064d-4*t13*t14*t3)*t8*wght+fnc(iq)
            Amat(iq,D1_RA) = 6.299605249474365d-1*t2*t6*t8*(3.1748021039
     1         363994d+0*gammaaa*t12*t15*(4.224610354362455d-4*rhoa*t13*
     2         t3/t7**3.166666666666667d+0+7.242189178907065d-4*t14*t20*
     3         t3/t12**2)*wght+3.1748021039363994d+0*gammaaa*t15*t20*wgh
     4         t)-1.0d+0*gammaaa*t12*t15*t6*t8*wght/rhoa**2-2.6456684199
     5         46998d-1*gammaaa*t12*t15*t16*wght/rhoa**1.333333333333333
     6         3d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 5.0d-1*t12*t15*t2*t6*t8*wght-1.43703234014
     1         37232d-4*t15*t2*t21*t3*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 1.0d+0*t12*t15*t2*t6*t8*wght-2.87406468028
     1         74464d-4*t15*t2*t21*t3*wght+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb+2.0d+0*gammaab+gammaaa
            t2 = rhob+rhoa
            t3 = 1/t2
            t4 = rhoa-rhob
            t5 = 1.0d+0-1.0d+0*t3*t4
            t6 = 1.0d+0*t3*t4+1.0d+0
            t7 = (3.149802624737183d-1*t6**1.6666666666666669d+0+3.14980
     1         2624737183d-1*t5**1.6666666666666669d+0)**5.0d-1
            t8 = 1/t7
            t9 = t1**5.0d-1
            t10 = 1/t2**6.666666666666666d-1
            t11 = abs(t2)
            t12 = 1/t11**3.333333333333333d-1
            t13 = 1.443307452126544d-2*t12+2.8435438314903855d-6*t10+2.5
     1         68d-3
            t14 = 1.7639938117590215d-2*t10*t12+5.411317332115466d+0*t12
     1         +1.816419932959077d-1*t10+1.0d+0
            t15 = 1/t14
            t16 = t13*t15+1.6670000000000001d-3
            t17 = 1/t16
            t18 = 1/t11**1.166666666666667d+0
            t19 = exp(-8.129082500000001d-4*t17*t18*t9)
            t20 = 1/t11
            t21 = -2.6456684199469993d-1*t1*t16*t19*t20*t8*wght/t2**1.33
     1         33333333333333d+0
            t22 = 1/t2**2
            t23 = 1.0d+0*t22*t4
            t24 = -1.0d+0*t3
            t25 = t5**6.666666666666666d-1
            t26 = -1.0d+0*t22*t4
            t27 = 1.0d+0*t3
            t28 = t6**6.666666666666666d-1
            t29 = 1/t7**3
            t30 = -7.937005259840999d-1*t1*t12*t16*t19*t22*t8*wght
            t31 = 1/t2**1.6666666666666669d+0
            t32 = 1/t2**3.333333333333333d-1
            t33 = t15*(-4.8110248404218126d-3*t20*t32-1.895695887660257d
     1         -6*t31)-t13*(-1.8037724440384886d+0*t20*t32-1.17599587450
     2         60142d-2*t12*t31-1.2109466219727179d-1*t31-5.879979372530
     3         072d-3*t20*t3)/t14**2
            t34 = t12*t3*t8*(7.937005259840999d-1*t1*t16*t19*(8.12908250
     1         0000001d-4*t18*t33*t9/t16**2+9.483929583333336d-4*t17*t2*
     2         t9/t11**3.166666666666667d+0)*wght+7.937005259840999d-1*t
     3         1*t19*t33*wght)
            t35 = 1/t11**1.5000000000000002d+0
            t36 = 7.937005259840999d-1*t12*t16*t19*t3*t8*wght-3.22602852
     1         8009071d-4*t19*t3*t35*t8*t9*wght
            fnc(iq) = 7.937005259840999d-1*t1*t12*t16*t19*t3*t8*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = -3.9685026299204995d-1*t1*t12*t16*t19*(5.24
     1         9671041228638d-1*(t27+t26)*t28+5.249671041228638d-1*(t24+
     2         t23)*t25)*t29*t3*wght+t34+t30+t21+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = -3.9685026299204995d-1*t1*t12*t16*t19*(5.24
     1         9671041228638d-1*(t26+t24)*t28+5.249671041228638d-1*t25*(
     2         t27+t23))*t29*t3*wght+t34+t30+t21+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t36+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = -6.452057056018142d-4*t19*t3*t35*t8*t9*wgh
     1         t+1.5874010519681997d+0*t12*t16*t19*t3*t8*wght+Cmat(iq,D1
     2         _GAB)
            Cmat(iq,D1_GBB) = t36+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa
            t2 = gammaaa**5.0d-1
            t3 = 1/rhoa**6.666666666666666d-1
            t4 = abs(rhoa)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/t4
            t14 = 1/rhoa**1.6666666666666669d+0
            t15 = 1/rhoa**3.333333333333333d-1
            t16 = (-4.8110248404218126d-3*t13*t15-1.895695887660257d-6*t
     1         14)*t8-(-1.1759958745060142d-2*t14*t5-1.8037724440384886d
     2         +0*t13*t15-1.2109466219727179d-1*t14-5.879979372530072d-3
     3         *t1*t13)*t6/t7**2
            fnc(iq) = 7.937005259840999d-1*gammaaa*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = t1*t5*(7.937005259840999d-1*gammaaa*t12*(8.
     1         129082500000001d-4*t11*t16*t2/t9**2+9.483929583333336d-4*
     2         rhoa*t10*t2/t4**3.166666666666667d+0)*t9*wght+7.937005259
     3         840999d-1*gammaaa*t12*t16*wght)-7.937005259840999d-1*gamm
     4         aaa*t12*t5*t9*wght/rhoa**2-2.6456684199469993d-1*gammaaa*
     5         t12*t13*t9*wght/rhoa**1.3333333333333333d+0+Amat(iq,D1_RA
     6         )
            Cmat(iq,D1_GAA) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*wght/t4**1.5000000000000002d+
     2         0+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob
            t2 = gammabb**5.0d-1
            t3 = 1/rhob**6.666666666666666d-1
            t4 = abs(rhob)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/t4
            t14 = 1/rhob**1.6666666666666669d+0
            t15 = 1/rhob**3.333333333333333d-1
            t16 = (-4.8110248404218126d-3*t13*t15-1.895695887660257d-6*t
     1         14)*t8-(-1.1759958745060142d-2*t14*t5-1.8037724440384886d
     2         +0*t13*t15-1.2109466219727179d-1*t14-5.879979372530072d-3
     3         *t1*t13)*t6/t7**2
            fnc(iq) = 7.937005259840999d-1*gammabb*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RB) = t1*t5*(7.937005259840999d-1*gammabb*t12*(8.
     1         129082500000001d-4*t11*t16*t2/t9**2+9.483929583333336d-4*
     2         rhob*t10*t2/t4**3.166666666666667d+0)*t9*wght+7.937005259
     3         840999d-1*gammabb*t12*t16*wght)-7.937005259840999d-1*gamm
     4         abb*t12*t5*t9*wght/rhob**2-2.6456684199469993d-1*gammabb*
     5         t12*t13*t9*wght/rhob**1.3333333333333333d+0+Amat(iq,D1_RB
     6         )
            Cmat(iq,D1_GBB) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*wght/t4**1.5000000000000002d+
     2         0+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_perdew86 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_5}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{0.01443307452126544\,{\it t_6}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_4}+0.002568}
C>    \over{0.01763993811759022\,{\it t_4}\,{\it t_6}
C>    +5.411317332115466\,{\it t_6}+0.1816419932959077\,{\it t_4}
C>    +1.0}}+0.001667\\\\ 
C>   {\it t_8} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_{10}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{12}} &=& {{0.01443307452126544\,{\it t_{11}}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_9}+0.002568}
C>    \over{0.01763993811759022\,{\it t_9}\,{\it t_{11}}
C>    +5.411317332115466\,{\it t_{11}}+0.1816419932959077\,{\it t_9}
C>    +1.0}}+0.001667\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_7}\,{\it t_6}
C>    \,{\it t_8}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{
C>    \it t_3}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}
C>    +0.3149802624737183\,\left(1.0-1.0\,{\it t_3}\,{
C>    \it t_2}\right)^{{{5}\over{3}}}}
C>    \,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{
C>    \it t_8}}}\over{{\it t_7}\,{\it t_5}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.7937005259840999\,{\it t_{12}}\,{\it t_{11}}
C>    \,\sigma_{ss}}\over{\rho_s\,e^{{{8.129082500000001 \times 10^{
C>    -4}\,\sqrt{\sigma_{ss}}}\over{{\it t_{12}}\,{
C>    \it t_{10}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew86_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = 3.968502629920499d-1
            t2 = 1/rhoa
            t3 = gammaaa**5.0d-1
            t4 = 6.299605249474366d-1
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = 7.937005259840998d-1
            t7 = abs(rhoa)
            t8 = 1/t7**3.333333333333333d-1
            t9 = 1.443307452126544d-2*t6*t8+2.8435438314903855d-6*t4*t5+
     1         2.568d-3
            t10 = 5.411317332115466d+0*t6*t8+8.819969058795109d-3*t5*t8+
     1         1.816419932959077d-1*t4*t5+1.0d+0
            t11 = 1/t10
            t12 = t11*t9+1.6670000000000001d-3
            t13 = 1/t12
            t14 = 1/t7**1.166666666666667d+0
            t15 = 1.9842513149602492d-1
            t16 = 1/rhoa**1.3333333333333333d+0
            t17 = exp(-7.242189178907065d-4*t13*t14*t3)
            t18 = 1/t7
            t19 = rhoa**2
            t20 = 1/t19
            t21 = 3.149802624737183d-1
            t22 = 1/rhoa**1.6666666666666669d+0
            t23 = 1/rhoa**3.333333333333333d-1
            t24 = -2.939989686265036d-3*t22*t8-1.8037724440384886d+0*t1*
     1         t18*t23-1.2109466219727179d-1*t21*t22-1.469994843132518d-
     2         3*t18*t2
            t25 = 1/t10**2
            t26 = -4.8110248404218126d-3*t1*t18*t23-1.895695887660257d-6
     1         *t21*t22
            t27 = t11*t26-t24*t25*t9
            t28 = 1/t7**3.166666666666667d+0
            t29 = 1/t12**2
            t30 = 7.242189178907065d-4*t14*t27*t29*t3+4.224610354362455d
     1         -4*rhoa*t13*t28*t3
            t31 = 1/t7**1.5000000000000002d+0
            t32 = 4.409447366578332d-1*gammaaa*t12*t17*t18*wght/rhoa**2.
     1         3333333333333334d+0
            t33 = 1/rhoa**3
            t34 = 1.5749013123685918d-1
            t35 = 1/rhoa**2.6666666666666666d+0
            t36 = -t25*(2.4499914052208632d-3*t35*t8+2.01824436995453d-1
     1         *t34*t35+1.9599931241766905d-3*t18*t20+2.4050299253846513
     2         d+0*t15*t16*t18)*t9+2*t24**2*t9/t10**3+t11*(3.15949314610
     3         04283d-6*t34*t35+6.414699787229083d-3*t15*t16*t18)-2*t24*
     4         t25*t26
            t37 = 6.299605249474365d-1*t2*t6*t8*(3.1748021039363994d+0*g
     1         ammaaa*t12*t17*(-6.688966394407221d-4*t13*t19*t3/t7**5.16
     2         6666666666667d+0+7.242189178907065d-4*t14*t29*t3*t36-8.44
     3         922070872491d-4*rhoa*t27*t28*t29*t3+2.1123051771812276d-4
     4         *t13*t28*t3-1.448437835781413d-3*t14*t27**2*t3/t12**3)*wg
     5         ht+3.1748021039363994d+0*gammaaa*t17*t36*wght+3.174802103
     6         9363994d+0*gammaaa*t12*t17*t30**2*wght+6.349604207872799d
     7         +0*gammaaa*t17*t27*t30*wght)
            t38 = 3.1498026247371824d-1*t20*t6*t8*(-6.349604207872799d+0
     1         *gammaaa*t12*t17*t30*wght-6.349604207872799d+0*gammaaa*t1
     2         7*t27*wght)
            t39 = 4.19973683298291d-1*t15*t16*t18*(-3.1748021039363994d+
     1         0*gammaaa*t12*t17*t30*wght-3.1748021039363994d+0*gammaaa*
     2         t17*t27*wght)
            t40 = 1.259921049894873d+0*t15*t16*t18*(-1.0582673679787998d
     1         +0*gammaaa*t12*t17*t30*wght-1.0582673679787998d+0*gammaaa
     2         *t17*t27*wght)
            t41 = 1/t7**2.166666666666667d+0
            t42 = 1/t3
            t43 = 6.299605249474365d-1*t2*t6*t8*(3.1748021039363994d+0*g
     1         ammaaa*t12*t17*(9.052736473633831d-5*t14*t27*t29*t42+5.28
     2         0762942953069d-5*rhoa*t13*t28*t42)*wght-2.874064680287446
     3         4d-4*t14*t17*t3*t30*wght+7.937005259840999d-1*t12*t17*t30
     4         *wght-2.8740646802874464d-4*t13*t14*t17*t27*t3*wght+7.937
     5         005259840999d-1*t17*t27*wght)-2.5d-1*t12*t17*t20*t6*t8*wg
     6         ht+6.035157649089219d-5*t1*t16*t17*t3*t41*wght+7.18516170
     7         0718616d-5*t17*t20*t3*t31*wght-3.333333333333333d-1*t12*t
     8         15*t16*t17*t18*wght
            t44 = 1/t7**2.666666666666667d+0
            t45 = 1.300907507941046d-8*t13*t17*t2*t44*wght-5.38887127553
     1         8962d-5*t17*t2*t31*t42*wght
            fnc(iq) = 4.000000000000001d+0*gammaaa*t1*t12*t2*exp(-7.2421
     1         89178907064d-4*t13*t14*t3)*t8*wght+fnc(iq)
            Amat(iq,D1_RA) = 6.299605249474365d-1*t2*t6*t8*(3.1748021039
     1         363994d+0*gammaaa*t12*t17*t30*wght+3.1748021039363994d+0*
     2         gammaaa*t17*t27*wght)-1.0d+0*gammaaa*t12*t17*t20*t6*t8*wg
     3         ht-1.3333333333333333d+0*gammaaa*t12*t15*t16*t17*t18*wght
     4         +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 5.0d-1*t12*t17*t2*t6*t8*wght-1.43703234014
     1         37232d-4*t17*t2*t3*t31*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 1.0d+0*t12*t17*t2*t6*t8*wght-2.87406468028
     1         74464d-4*t17*t2*t3*t31*wght+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 7.222222222222222d-1*gammaaa*t12*t17*t3
     1         3*t6*t8*wght+t40+t39+t38+t37+t32+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.2777777777777777d+0*gammaaa*t12*t17*t
     1         33*t6*t8*wght+t40+t39+t38+t37+t32+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t43+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 6.299605249474365d-1*t2*t6*t8*(3.17480
     1         21039363994d+0*gammaaa*t12*t17*(1.8105472947267662d-4*t14
     2         *t27*t29*t42+1.0561525885906138d-4*rhoa*t13*t28*t42)*wght
     3         -5.748129360574893d-4*t14*t17*t3*t30*wght+1.5874010519681
     4         997d+0*t12*t17*t30*wght-5.748129360574893d-4*t13*t14*t17*
     5         t27*t3*wght+1.5874010519681997d+0*t17*t27*wght)-5.0d-1*t1
     6         2*t17*t20*t6*t8*wght+1.2070315298178438d-4*t1*t16*t17*t3*
     7         t41*wght+1.4370323401437232d-4*t17*t20*t3*t31*wght-6.6666
     8         66666666666d-1*t12*t15*t16*t17*t18*wght+Cmat2(iq,D2_RA_GA
     9         B)
            Cmat2(iq,D2_RA_GBB) = t43+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t45+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = 2.601815015882092d-8*t13*t17*t2*t44*w
     1         ght-1.0777742551077923d-4*t17*t2*t31*t42*wght+Cmat2(iq,D2
     2         _GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t45+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 5.203630031764184d-8*t13*t17*t2*t44*w
     1         ght-2.1555485102155847d-4*t17*t2*t31*t42*wght+Cmat2(iq,D2
     2         _GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb+2.0d+0*gammaab+gammaaa
            t2 = rhob+rhoa
            t3 = 1/t2
            t4 = rhoa-rhob
            t5 = 1.0d+0-1.0d+0*t3*t4
            t6 = 1.0d+0*t3*t4+1.0d+0
            t7 = (3.149802624737183d-1*t6**1.6666666666666669d+0+3.14980
     1         2624737183d-1*t5**1.6666666666666669d+0)**5.0d-1
            t8 = 1/t7
            t9 = t1**5.0d-1
            t10 = 1/t2**6.666666666666666d-1
            t11 = abs(t2)
            t12 = 1/t11**3.333333333333333d-1
            t13 = 1.443307452126544d-2*t12+2.8435438314903855d-6*t10+2.5
     1         68d-3
            t14 = 1.7639938117590215d-2*t10*t12+5.411317332115466d+0*t12
     1         +1.816419932959077d-1*t10+1.0d+0
            t15 = 1/t14
            t16 = t13*t15+1.6670000000000001d-3
            t17 = 1/t16
            t18 = 1/t11**1.166666666666667d+0
            t19 = exp(-8.129082500000001d-4*t17*t18*t9)
            t20 = 1/t2**1.3333333333333333d+0
            t21 = 1/t11
            t22 = -2.6456684199469993d-1*t1*t16*t19*t20*t21*t8*wght
            t23 = t2**2
            t24 = 1/t23
            t25 = 1.0d+0*t24*t4
            t26 = -1.0d+0*t3
            t27 = t26+t25
            t28 = t5**6.666666666666666d-1
            t29 = -1.0d+0*t24*t4
            t30 = 1.0d+0*t3
            t31 = t30+t29
            t32 = t6**6.666666666666666d-1
            t33 = 5.249671041228638d-1*t31*t32+5.249671041228638d-1*t27*
     1         t28
            t34 = 1/t7**3
            t35 = -7.937005259840999d-1*t1*t12*t16*t19*t24*t8*wght
            t36 = 1/t2**1.6666666666666669d+0
            t37 = 1/t2**3.333333333333333d-1
            t38 = -1.8037724440384886d+0*t21*t37-1.1759958745060142d-2*t
     1         12*t36-1.2109466219727179d-1*t36-5.879979372530072d-3*t21
     2         *t3
            t39 = 1/t14**2
            t40 = -4.8110248404218126d-3*t21*t37-1.895695887660257d-6*t3
     1         6
            t41 = t15*t40-t13*t38*t39
            t42 = 1/t11**3.166666666666667d+0
            t43 = 1/t16**2
            t44 = 8.129082500000001d-4*t18*t41*t43*t9+9.483929583333336d
     1         -4*t17*t2*t42*t9
            t45 = 7.937005259840999d-1*t1*t16*t19*t44*wght+7.93700525984
     1         0999d-1*t1*t19*t41*wght
            t46 = t12*t3*t45*t8
            t47 = t30+t25
            t48 = t29+t26
            t49 = 5.249671041228638d-1*t32*t48+5.249671041228638d-1*t28*
     1         t47
            t50 = 1/t11**1.5000000000000002d+0
            t51 = 7.937005259840999d-1*t12*t16*t19*t3*t8*wght-3.22602852
     1         8009071d-4*t19*t3*t50*t8*t9*wght
            t52 = 8.818894733156665d-1*t1*t16*t19*t21*t8*wght/t2**2.3333
     1         333333333334d+0
            t53 = 1/t7**5
            t54 = 1/t2**3
            t55 = 1.5874010519681997d+0*t1*t12*t16*t19*t54*t8*wght
            t56 = 1/t2**2.6666666666666666d+0
            t57 = -t13*t39*(1.9599931241766905d-2*t12*t56+2.018244369954
     1         53d-1*t56+1.5679944993413522d-2*t21*t24+2.405029925384651
     2         3d+0*t20*t21)+t15*(3.1594931461004283d-6*t56+6.4146997872
     3         29083d-3*t20*t21)-2*t38*t39*t40+2*t13*t38**2/t14**3
            t58 = t12*t3*t8*(7.937005259840999d-1*t1*t16*t19*(8.12908250
     1         0000001d-4*t18*t43*t57*t9-1.8967859166666673d-3*t2*t41*t4
     2         2*t43*t9+9.483929583333336d-4*t17*t42*t9-1.62581650000000
     3         02d-3*t18*t41**2*t9/t16**3-3.003244368055557d-3*t17*t23*t
     4         9/t11**5.166666666666667d+0)*wght+7.937005259840999d-1*t1
     5         *t19*t57*wght+7.937005259840999d-1*t1*t16*t19*t44**2*wght
     6         +1.5874010519681997d+0*t1*t19*t41*t44*wght)
            t59 = t12*t24*t8*(-1.5874010519681997d+0*t1*t16*t19*t44*wght
     1         -1.5874010519681997d+0*t1*t19*t41*wght)
            t60 = 3.333333333333333d-1*t20*t21*t8*(-7.937005259840999d-1
     1         *t1*t16*t19*t44*wght-7.937005259840999d-1*t1*t19*t41*wght
     2         )
            t61 = t20*t21*t8*(-2.6456684199469993d-1*t1*t16*t19*t44*wght
     1         -2.6456684199469993d-1*t1*t19*t41*wght)
            t62 = -3.9685026299204995d-1*t1*t19*t33*t41*wght
            t63 = 1/t5**3.333333333333333d-1
            t64 = -2.0d+0*t4*t54
            t65 = 2.0d+0*t24
            t66 = 1/t6**3.333333333333333d-1
            t67 = 2.0d+0*t4*t54
            t68 = -2.0d+0*t24
            t69 = -3.9685026299204995d-1*t1*t16*t19*t33*t44*wght
            t70 = -5.0d-1*t12*t3*t34*t45*t49
            t71 = 1/t11**2.166666666666667d+0
            t72 = 1.0753428426696904d-4*t19*t20*t71*t8*t9*wght
            t73 = 3.226028528009071d-4*t19*t24*t50*t8*t9*wght
            t74 = -2.6456684199469993d-1*t16*t19*t20*t21*t8*wght
            t75 = -7.937005259840999d-1*t12*t16*t19*t24*t8*wght
            t76 = 1/t9
            t77 = t12*t3*t8*(-3.226028528009071d-4*t18*t19*t44*t9*wght-3
     1         .226028528009071d-4*t17*t18*t19*t41*t9*wght+7.93700525984
     2         0999d-1*t1*t16*t19*(4.0645412500000005d-4*t18*t41*t43*t76
     3         +4.741964791666668d-4*t17*t2*t42*t76)*wght+7.937005259840
     4         999d-1*t16*t19*t44*wght+7.937005259840999d-1*t19*t41*wght
     5         )
            t78 = 1.6130142640045356d-4*t19*t3*t33*t34*t50*t9*wght-3.968
     1         5026299204995d-1*t12*t16*t19*t3*t33*t34*wght+t77+t75+t74+
     2         t73+t72
            t79 = 2.1506856853393808d-4*t19*t20*t71*t8*t9*wght
            t80 = 6.452057056018142d-4*t19*t24*t50*t8*t9*wght
            t81 = -5.291336839893999d-1*t16*t19*t20*t21*t8*wght
            t82 = -1.5874010519681997d+0*t12*t16*t19*t24*t8*wght
            t83 = t12*t3*t8*(-6.452057056018142d-4*t18*t19*t44*t9*wght-6
     1         .452057056018142d-4*t17*t18*t19*t41*t9*wght+7.93700525984
     2         0999d-1*t1*t16*t19*(8.129082500000001d-4*t18*t41*t43*t76+
     3         9.483929583333336d-4*t17*t2*t42*t76)*wght+1.5874010519681
     4         997d+0*t16*t19*t44*wght+1.5874010519681997d+0*t19*t41*wgh
     5         t)
            t84 = 1.6130142640045356d-4*t19*t3*t34*t49*t50*t9*wght-3.968
     1         5026299204995d-1*t12*t16*t19*t3*t34*t49*wght+t77+t75+t74+
     2         t73+t72
            t85 = 1/t11**2.666666666666667d+0
            t86 = 1.311232602576965d-7*t17*t19*t3*t8*t85*wght-4.83904279
     1         20136064d-4*t19*t3*t50*t76*t8*wght
            t87 = 2.62246520515393d-7*t17*t19*t3*t8*t85*wght-9.678085584
     1         027213d-4*t19*t3*t50*t76*t8*wght
            fnc(iq) = 7.937005259840999d-1*t1*t12*t16*t19*t3*t8*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = -3.9685026299204995d-1*t1*t12*t16*t19*t3*t3
     1         3*t34*wght+t46+t35+t22+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = -3.9685026299204995d-1*t1*t12*t16*t19*t3*t3
     1         4*t49*wght+t46+t35+t22+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t51+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = -6.452057056018142d-4*t19*t3*t50*t8*t9*wgh
     1         t+1.5874010519681997d+0*t12*t16*t19*t3*t8*wght+Cmat(iq,D1
     2         _GAB)
            Cmat(iq,D1_GBB) = t51+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = t12*t3*t34*(-3.9685026299204995d-1*t1*t
     1         16*t19*(5.249671041228638d-1*t32*(t68+t67)+3.499780694152
     2         4257d-1*t31**2*t66+5.249671041228638d-1*t28*(t65+t64)+3.4
     3         997806941524257d-1*t27**2*t63)*wght+t69+t62)+5.9527539448
     4         80749d-1*t1*t12*t16*t19*t3*t33**2*t53*wght+7.937005259840
     5         999d-1*t1*t12*t16*t19*t24*t33*t34*wght+2.6456684199469993
     6         d-1*t1*t16*t19*t20*t21*t33*t34*wght+t61+t60+t59+t58+t55+t
     7         52-5.0d-1*t12*t3*t33*t34*t45+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t12*t3*t34*(-3.9685026299204995d-1*t1*t
     1         16*t19*(3.4997806941524257d-1*t31*t48*t66+3.4997806941524
     2         257d-1*t27*t47*t63+1.0499342082457278d+0*t32*t4*t54-1.049
     3         9342082457278d+0*t28*t4*t54)*wght+t69+t62)+t12*t24*t34*(3
     4         .9685026299204995d-1*t1*t16*t19*t49*wght+3.96850262992049
     5         95d-1*t1*t16*t19*t33*wght)+t20*t21*t34*(1.322834209973499
     6         7d-1*t1*t16*t19*t49*wght+1.3228342099734997d-1*t1*t16*t19
     7         *t33*wght)+5.952753944880749d-1*t1*t12*t16*t19*t3*t33*t49
     8         *t53*wght+t70+t61+t60+t59+t58+t55+t52+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = t12*t3*t34*(-3.9685026299204995d-1*t1*t
     1         16*t19*(5.249671041228638d-1*t28*(t68+t64)+5.249671041228
     2         638d-1*t32*(t67+t65)+3.4997806941524257d-1*t48**2*t66+3.4
     3         997806941524257d-1*t47**2*t63)*wght-3.9685026299204995d-1
     4         *t1*t16*t19*t44*t49*wght-3.9685026299204995d-1*t1*t19*t41
     5         *t49*wght)+5.952753944880749d-1*t1*t12*t16*t19*t3*t49**2*
     6         t53*wght+7.937005259840999d-1*t1*t12*t16*t19*t24*t34*t49*
     7         wght+2.6456684199469993d-1*t1*t16*t19*t20*t21*t34*t49*wgh
     8         t+t70+t61+t60+t59+t58+t55+t52+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t78+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 3.226028528009071d-4*t19*t3*t33*t34*t5
     1         0*t9*wght-7.937005259840999d-1*t12*t16*t19*t3*t33*t34*wgh
     2         t+t83+t82+t81+t80+t79+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t78+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t84+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 3.226028528009071d-4*t19*t3*t34*t49*t5
     1         0*t9*wght-7.937005259840999d-1*t12*t16*t19*t3*t34*t49*wgh
     2         t+t83+t82+t81+t80+t79+Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t84+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t86+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t87+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t86+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 5.24493041030786d-7*t17*t19*t3*t8*t85
     1         *wght-1.9356171168054426d-3*t19*t3*t50*t76*t8*wght+Cmat2(
     2         iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t87+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t86+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa
            t2 = gammaaa**5.0d-1
            t3 = 1/rhoa**6.666666666666666d-1
            t4 = abs(rhoa)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/rhoa**1.3333333333333333d+0
            t14 = 1/t4
            t15 = rhoa**2
            t16 = 1/t15
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/rhoa**3.333333333333333d-1
            t19 = -1.1759958745060142d-2*t17*t5-1.8037724440384886d+0*t1
     1         4*t18-1.2109466219727179d-1*t17-5.879979372530072d-3*t1*t
     2         14
            t20 = 1/t7**2
            t21 = -4.8110248404218126d-3*t14*t18-1.895695887660257d-6*t1
     1         7
            t22 = t21*t8-t19*t20*t6
            t23 = 1/t4**3.166666666666667d+0
            t24 = 1/t9**2
            t25 = 8.129082500000001d-4*t11*t2*t22*t24+9.483929583333336d
     1         -4*rhoa*t10*t2*t23
            t26 = 1/t4**1.5000000000000002d+0
            t27 = 1/rhoa**2.6666666666666666d+0
            t28 = (3.1594931461004283d-6*t27+6.414699787229083d-3*t13*t1
     1         4)*t8+2*t19**2*t6/t7**3-t20*(1.9599931241766905d-2*t27*t5
     2         +2.01824436995453d-1*t27+1.5679944993413522d-2*t14*t16+2.
     3         4050299253846513d+0*t13*t14)*t6-2*t19*t20*t21
            t29 = 1/t2
            fnc(iq) = 7.937005259840999d-1*gammaaa*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = t1*t5*(7.937005259840999d-1*gammaaa*t12*t25
     1         *t9*wght+7.937005259840999d-1*gammaaa*t12*t22*wght)-7.937
     2         005259840999d-1*gammaaa*t12*t16*t5*t9*wght-2.645668419946
     3         9993d-1*gammaaa*t12*t13*t14*t9*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*t26*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t1*t5*(7.937005259840999d-1*gammaaa*t12
     1         *(-1.6258165000000002d-3*t11*t2*t22**2/t9**3-3.0032443680
     2         55557d-3*t10*t15*t2/t4**5.166666666666667d+0+8.1290825000
     3         00001d-4*t11*t2*t24*t28-1.8967859166666673d-3*rhoa*t2*t22
     4         *t23*t24+9.483929583333336d-4*t10*t2*t23)*t9*wght+7.93700
     5         5259840999d-1*gammaaa*t12*t25**2*t9*wght+7.93700525984099
     6         9d-1*gammaaa*t12*t28*wght+1.5874010519681997d+0*gammaaa*t
     7         12*t22*t25*wght)+t13*t14*(-2.6456684199469993d-1*gammaaa*
     8         t12*t25*t9*wght-2.6456684199469993d-1*gammaaa*t12*t22*wgh
     9         t)+3.333333333333333d-1*t13*t14*(-7.937005259840999d-1*ga
     :         mmaaa*t12*t25*t9*wght-7.937005259840999d-1*gammaaa*t12*t2
     ;         2*wght)+t16*t5*(-1.5874010519681997d+0*gammaaa*t12*t25*t9
     <         *wght-1.5874010519681997d+0*gammaaa*t12*t22*wght)+1.58740
     =         10519681997d+0*gammaaa*t12*t5*t9*wght/rhoa**3+8.818894733
     >         156665d-1*gammaaa*t12*t14*t9*wght/rhoa**2.333333333333333
     ?         4d+0+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t1*t5*(7.937005259840999d-1*gammaaa*t1
     1         2*(4.0645412500000005d-4*t11*t22*t24*t29+4.74196479166666
     2         8d-4*rhoa*t10*t23*t29)*t9*wght+7.937005259840999d-1*t12*t
     3         25*t9*wght-3.226028528009071d-4*t11*t12*t2*t25*wght-3.226
     4         028528009071d-4*t10*t11*t12*t2*t22*wght+7.937005259840999
     5         d-1*t12*t22*wght)-7.937005259840999d-1*t12*t16*t5*t9*wght
     6         -2.6456684199469993d-1*t12*t13*t14*t9*wght+1.075342842669
     7         6904d-4*t12*t13*t2*wght/t4**2.166666666666667d+0+3.226028
     8         528009071d-4*t12*t16*t2*t26*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.311232602576965d-7*t1*t10*t12*wght/
     1         t4**2.666666666666667d+0-4.8390427920136064d-4*t1*t12*t26
     2         *t29*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob
            t2 = gammabb**5.0d-1
            t3 = 1/rhob**6.666666666666666d-1
            t4 = abs(rhob)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/rhob**1.3333333333333333d+0
            t14 = 1/t4
            t15 = rhob**2
            t16 = 1/t15
            t17 = 1/rhob**1.6666666666666669d+0
            t18 = 1/rhob**3.333333333333333d-1
            t19 = -1.1759958745060142d-2*t17*t5-1.8037724440384886d+0*t1
     1         4*t18-1.2109466219727179d-1*t17-5.879979372530072d-3*t1*t
     2         14
            t20 = 1/t7**2
            t21 = -4.8110248404218126d-3*t14*t18-1.895695887660257d-6*t1
     1         7
            t22 = t21*t8-t19*t20*t6
            t23 = 1/t4**3.166666666666667d+0
            t24 = 1/t9**2
            t25 = 8.129082500000001d-4*t11*t2*t22*t24+9.483929583333336d
     1         -4*rhob*t10*t2*t23
            t26 = 1/t4**1.5000000000000002d+0
            t27 = 1/rhob**2.6666666666666666d+0
            t28 = (3.1594931461004283d-6*t27+6.414699787229083d-3*t13*t1
     1         4)*t8+2*t19**2*t6/t7**3-t20*(1.9599931241766905d-2*t27*t5
     2         +2.01824436995453d-1*t27+1.5679944993413522d-2*t14*t16+2.
     3         4050299253846513d+0*t13*t14)*t6-2*t19*t20*t21
            t29 = 1/t2
            fnc(iq) = 7.937005259840999d-1*gammabb*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RB) = t1*t5*(7.937005259840999d-1*gammabb*t12*t25
     1         *t9*wght+7.937005259840999d-1*gammabb*t12*t22*wght)-7.937
     2         005259840999d-1*gammabb*t12*t16*t5*t9*wght-2.645668419946
     3         9993d-1*gammabb*t12*t13*t14*t9*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*t26*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t1*t5*(7.937005259840999d-1*gammabb*t12
     1         *(-1.6258165000000002d-3*t11*t2*t22**2/t9**3-3.0032443680
     2         55557d-3*t10*t15*t2/t4**5.166666666666667d+0+8.1290825000
     3         00001d-4*t11*t2*t24*t28-1.8967859166666673d-3*rhob*t2*t22
     4         *t23*t24+9.483929583333336d-4*t10*t2*t23)*t9*wght+7.93700
     5         5259840999d-1*gammabb*t12*t25**2*t9*wght+7.93700525984099
     6         9d-1*gammabb*t12*t28*wght+1.5874010519681997d+0*gammabb*t
     7         12*t22*t25*wght)+t13*t14*(-2.6456684199469993d-1*gammabb*
     8         t12*t25*t9*wght-2.6456684199469993d-1*gammabb*t12*t22*wgh
     9         t)+3.333333333333333d-1*t13*t14*(-7.937005259840999d-1*ga
     :         mmabb*t12*t25*t9*wght-7.937005259840999d-1*gammabb*t12*t2
     ;         2*wght)+t16*t5*(-1.5874010519681997d+0*gammabb*t12*t25*t9
     <         *wght-1.5874010519681997d+0*gammabb*t12*t22*wght)+1.58740
     =         10519681997d+0*gammabb*t12*t5*t9*wght/rhob**3+8.818894733
     >         156665d-1*gammabb*t12*t14*t9*wght/rhob**2.333333333333333
     ?         4d+0+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t1*t5*(7.937005259840999d-1*gammabb*t1
     1         2*(4.0645412500000005d-4*t11*t22*t24*t29+4.74196479166666
     2         8d-4*rhob*t10*t23*t29)*t9*wght+7.937005259840999d-1*t12*t
     3         25*t9*wght-3.226028528009071d-4*t11*t12*t2*t25*wght-3.226
     4         028528009071d-4*t10*t11*t12*t2*t22*wght+7.937005259840999
     5         d-1*t12*t22*wght)-7.937005259840999d-1*t12*t16*t5*t9*wght
     6         -2.6456684199469993d-1*t12*t13*t14*t9*wght+1.075342842669
     7         6904d-4*t12*t13*t2*wght/t4**2.166666666666667d+0+3.226028
     8         528009071d-4*t12*t16*t2*t26*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.311232602576965d-7*t1*t10*t12*wght/
     1         t4**2.666666666666667d+0-4.8390427920136064d-4*t1*t12*t26
     2         *t29*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_perdew86 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& \rho_\alpha-\rho_\beta\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_5}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{0.01443307452126544\,{\it t_6}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_4}+0.002568}
C>    \over{0.01763993811759022\,{\it t_4}\,{\it t_6}
C>    +5.411317332115466\,{\it t_6}+0.1816419932959077\,{\it t_4}
C>    +1.0}}+0.001667\\\\ 
C>   {\it t_8} &=& \sigma_{\beta\beta}+2.0\,\sigma_{\alpha\beta}
C>    +\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_9} &=& {{1}\over{\rho_s^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{10}} &=& \left| \rho_s\right| \\\\ 
C>   {\it t_{11}} &=& {{1}\over{{\it t_{10}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{12}} &=& {{0.01443307452126544\,{\it t_{11}}
C>    +2.843543831490385 \times 10^{-6}\,{\it t_9}+0.002568}
C>    \over{0.01763993811759022\,{\it t_9}\,{\it t_{11}}
C>    +5.411317332115466\,{\it t_{11}}+0.1816419932959077\,{\it t_9}
C>    +1.0}}+0.001667\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_7}\,{\it t_6}
C>    \,{\it t_8}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{
C>    \it t_3}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}
C>    +0.3149802624737183\,\left(1.0-1.0\,{\it t_3}\,{
C>    \it t_2}\right)^{{{5}\over{3}}}}
C>    \,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{
C>    \it t_8}}}\over{{\it t_7}\,{\it t_5}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{0.7937005259840999\,{\it t_{12}}\,{\it t_{11}}
C>    \,\sigma_{ss}}\over{\rho_s\,e^{{{8.129082500000001 \times 10^{
C>    -4}\,\sqrt{\sigma_{ss}}}\over{{\it t_{12}}\,{
C>    \it t_{10}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_perdew86_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = 3.968502629920499d-1
            t2 = 1/rhoa
            t3 = gammaaa**5.0d-1
            t4 = 6.299605249474366d-1
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = 7.937005259840998d-1
            t7 = abs(rhoa)
            t8 = 1/t7**3.333333333333333d-1
            t9 = 1.443307452126544d-2*t6*t8+2.8435438314903855d-6*t4*t5+
     1         2.568d-3
            t10 = 5.411317332115466d+0*t6*t8+8.819969058795109d-3*t5*t8+
     1         1.816419932959077d-1*t4*t5+1.0d+0
            t11 = 1/t10
            t12 = t11*t9+1.6670000000000001d-3
            t13 = 1/t12
            t14 = 1/t7**1.166666666666667d+0
            t15 = 1.9842513149602492d-1
            t16 = 1/rhoa**1.3333333333333333d+0
            t17 = exp(-7.242189178907065d-4*t13*t14*t3)
            t18 = 1/t7
            t19 = rhoa**2
            t20 = 1/t19
            t21 = 3.149802624737183d-1
            t22 = 1/rhoa**1.6666666666666669d+0
            t23 = 1/rhoa**3.333333333333333d-1
            t24 = -2.939989686265036d-3*t22*t8-1.8037724440384886d+0*t1*
     1         t18*t23-1.2109466219727179d-1*t21*t22-1.469994843132518d-
     2         3*t18*t2
            t25 = 1/t10**2
            t26 = -4.8110248404218126d-3*t1*t18*t23-1.895695887660257d-6
     1         *t21*t22
            t27 = t11*t26-t24*t25*t9
            t28 = 1/t7**3.166666666666667d+0
            t29 = 1/t12**2
            t30 = 7.242189178907065d-4*t14*t27*t29*t3+4.224610354362455d
     1         -4*rhoa*t13*t28*t3
            t31 = 3.1748021039363994d+0*gammaaa*t12*t17*t30*wght+3.17480
     1         21039363994d+0*gammaaa*t17*t27*wght
            t32 = 1/t7**1.5000000000000002d+0
            t33 = 9.921256574801247d-2
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = 4.444444444444445d+0*gammaaa*t12*t17*t18*t33*t34*wght
            t36 = rhoa**3
            t37 = 1/t36
            t38 = t30**2
            t39 = t24**2
            t40 = 1/t10**3
            t41 = 1.5749013123685918d-1
            t42 = 1/rhoa**2.6666666666666666d+0
            t43 = 2.4499914052208632d-3*t42*t8+2.01824436995453d-1*t41*t
     1         42+1.9599931241766905d-3*t18*t20+2.4050299253846513d+0*t1
     2         5*t16*t18
            t44 = 3.1594931461004283d-6*t41*t42+6.414699787229083d-3*t15
     1         *t16*t18
            t45 = -t25*t43*t9+2*t39*t40*t9+t11*t44-2*t24*t25*t26
            t46 = 1/t7**5.166666666666667d+0
            t47 = t27**2
            t48 = 1/t12**3
            t49 = -1.448437835781413d-3*t14*t3*t47*t48-6.688966394407221
     1         d-4*t13*t19*t3*t46+7.242189178907065d-4*t14*t29*t3*t45-8.
     2         44922070872491d-4*rhoa*t27*t28*t29*t3+2.1123051771812276d
     3         -4*t13*t28*t3
            t50 = 6.299605249474365d-1*t2*t6*t8*(3.1748021039363994d+0*g
     1         ammaaa*t12*t17*t49*wght+3.1748021039363994d+0*gammaaa*t17
     2         *t45*wght+3.1748021039363994d+0*gammaaa*t12*t17*t38*wght+
     3         6.349604207872799d+0*gammaaa*t17*t27*t30*wght)
            t51 = -6.349604207872799d+0*gammaaa*t12*t17*t30*wght-6.34960
     1         4207872799d+0*gammaaa*t17*t27*wght
            t52 = 3.1498026247371824d-1*t20*t51*t6*t8
            t53 = -3.1748021039363994d+0*gammaaa*t12*t17*t30*wght-3.1748
     1         021039363994d+0*gammaaa*t17*t27*wght
            t54 = 4.19973683298291d-1*t15*t16*t18*t53
            t55 = -1.0582673679787998d+0*gammaaa*t12*t17*t30*wght-1.0582
     1         673679787998d+0*gammaaa*t17*t27*wght
            t56 = 1.259921049894873d+0*t15*t16*t18*t55
            t57 = 1/t7**2.166666666666667d+0
            t58 = 1/t3
            t59 = 9.052736473633831d-5*t14*t27*t29*t58+5.280762942953069
     1         d-5*rhoa*t13*t28*t58
            t60 = 6.299605249474365d-1*t2*t6*t8*(3.1748021039363994d+0*g
     1         ammaaa*t12*t17*t59*wght-2.8740646802874464d-4*t14*t17*t3*
     2         t30*wght+7.937005259840999d-1*t12*t17*t30*wght-2.87406468
     3         02874464d-4*t13*t14*t17*t27*t3*wght+7.937005259840999d-1*
     4         t17*t27*wght)-2.5d-1*t12*t17*t20*t6*t8*wght+6.03515764908
     5         9219d-5*t1*t16*t17*t3*t57*wght+7.185161700718616d-5*t17*t
     6         20*t3*t32*wght-3.333333333333333d-1*t12*t15*t16*t17*t18*w
     7         ght
            t61 = 1.8105472947267662d-4*t14*t27*t29*t58+1.05615258859061
     1         38d-4*rhoa*t13*t28*t58
            t62 = 1/t7**2.666666666666667d+0
            t63 = 1.300907507941046d-8*t13*t17*t2*t62*wght-5.38887127553
     1         8962d-5*t17*t2*t32*t58*wght
            t64 = 1/rhoa**3.3333333333333337d+0
            t65 = -8.671913154270722d-1*gammaaa*t12*t17*t18*t64*wght
            t66 = 1/rhoa**4
            t67 = 3.1498026247371824d-1*t20*t6*t8*(-9.524406311809198d+0
     1         *gammaaa*t12*t17*t49*wght-9.524406311809198d+0*gammaaa*t1
     2         7*t45*wght-9.524406311809198d+0*gammaaa*t12*t17*t38*wght-
     3         1.9048812623618397d+1*gammaaa*t17*t27*t30*wght)
            t68 = 4.19973683298291d-1*t15*t16*t18*(-6.349604207872799d+0
     1         *gammaaa*t12*t17*t49*wght-6.349604207872799d+0*gammaaa*t1
     2         7*t45*wght-6.349604207872799d+0*gammaaa*t12*t17*t38*wght-
     3         1.26992084157456d+1*gammaaa*t17*t27*t30*wght)
            t69 = 1.259921049894873d+0*t15*t16*t18*(-1.0582673679787998d
     1         +0*gammaaa*t12*t17*t49*wght-1.0582673679787998d+0*gammaaa
     2         *t17*t45*wght-1.0582673679787998d+0*gammaaa*t12*t17*t38*w
     3         ght-2.1165347359575995d+0*gammaaa*t17*t27*t30*wght)
            t70 = -9.799385943626789d-1*t18*t33*t34*t53
            t71 = 1.259921049894873d+0*t18*t33*t34*(3.527557893262666d+0
     1         *gammaaa*t12*t17*t30*wght+3.527557893262666d+0*gammaaa*t1
     2         7*t27*wght)
            t72 = 7.874506561842959d-2
            t73 = 1/rhoa**3.6666666666666664d+0
            t74 = -t25*(-3.266655206961151d-3*t73*t8-5.381984986545413d-
     1         1*t72*t73-3.3483215871351796d-3*t18*t37-5.611736492564187
     2         d+0*t18*t33*t34)*t9+6*t24*t40*t43*t9-6*t24**3*t9/t10**4+t
     3         11*(-8.425315056267809d-6*t72*t73-1.496763283686786d-2*t1
     4         8*t33*t34)-3*t24*t25*t44-3*t25*t26*t43+6*t26*t39*t40
            t75 = 6.299605249474365d-1*t2*t6*t8*(3.1748021039363994d+0*g
     1         ammaaa*t12*t17*(7.242189178907065d-4*t14*t29*t3*t74+1.727
     2         982985221866d-3*t13*t3*t36/t7**7.166666666666667d+0+2.534
     3         766212617473d-3*rhoa*t28*t3*t47*t48-4.3453135073442384d-3
     4         *t14*t27*t3*t45*t48+2.006689918322166d-3*t19*t27*t29*t3*t
     5         46-1.003344959161083d-3*rhoa*t13*t3*t46-1.267383106308736
     6         5d-3*rhoa*t28*t29*t3*t45-6.336915531543683d-4*t27*t28*t29
     7         *t3+4.3453135073442384d-3*t14*t27**3*t3/t12**4)*wght+3.17
     8         48021039363994d+0*gammaaa*t17*t74*wght+9.524406311809198d
     9         +0*gammaaa*t12*t17*t30*t49*wght+9.524406311809198d+0*gamm
     :         aaa*t17*t27*t49*wght+9.524406311809198d+0*gammaaa*t17*t30
     ;         *t45*wght+9.524406311809198d+0*gammaaa*t17*t27*t38*wght+3
     <         .1748021039363994d+0*gammaaa*t12*t17*t30**3*wght)
            t76 = 6.349604207872799d+0*gammaaa*t17*t27*wght
            t77 = 6.349604207872799d+0*gammaaa*t12*t17*t30*wght
            t78 = 1.5749013123685912d-1*t37*t6*(t77+t76-2*t51)*t8
            t79 = 4.19973683298291d-1*t18*t33*t34*(t77+t76-7*t55)
            t80 = -2.011719216363073d-4*t15*t17*t3*t34*t57*wght
            t81 = -7.185161700718616d-5*t17*t3*t32*t37*wght
            t82 = 1.1111111111111112d+0*t12*t17*t18*t33*t34*wght
            t83 = 2.5d-1*t12*t17*t37*t6*t8*wght
            t84 = -1.5874010519681997d+0*t17*t27*wght
            t85 = 5.748129360574893d-4*t13*t14*t17*t27*t3*wght
            t86 = 5.748129360574893d-4*t14*t17*t3*t30*wght
            t87 = -1.5874010519681997d+0*t12*t17*t30*wght
            t88 = 3.1498026247371824d-1*t20*t6*t8*(-6.349604207872799d+0
     1         *gammaaa*t12*t17*t59*wght+t87+t86+t85+t84)
            t89 = 4.19973683298291d-1*t15*t16*t18*(-3.1748021039363994d+
     1         0*gammaaa*t12*t17*t59*wght+2.8740646802874464d-4*t14*t17*
     2         t3*t30*wght-7.937005259840999d-1*t12*t17*t30*wght+2.87406
     3         46802874464d-4*t13*t14*t17*t27*t3*wght-7.937005259840999d
     4         -1*t17*t27*wght)
            t90 = 1.259921049894873d+0*t15*t16*t18*(-1.0582673679787998d
     1         +0*gammaaa*t12*t17*t59*wght+9.580215600958156d-5*t14*t17*
     2         t3*t30*wght-2.6456684199469993d-1*t12*t17*t30*wght+9.5802
     3         15600958156d-5*t13*t14*t17*t27*t3*wght-2.6456684199469993
     4         d-1*t17*t27*wght)
            t91 = 6.299605249474365d-1*t2*t6*t8*(6.349604207872799d+0*ga
     1         mmaaa*t12*t17*t30*t59*wght+6.349604207872799d+0*gammaaa*t
     2         17*t27*t59*wght+3.1748021039363994d+0*gammaaa*t12*t17*(-1
     3         .8105472947267662d-4*t14*t47*t48*t58-8.361207993009026d-5
     4         *t13*t19*t46*t58+9.052736473633831d-5*t14*t29*t45*t58-1.0
     5         561525885906138d-4*rhoa*t27*t28*t29*t58+2.640381471476534
     6         5d-5*t13*t28*t58)*wght-2.8740646802874464d-4*t14*t17*t3*t
     7         49*wght+7.937005259840999d-1*t12*t17*t49*wght-2.874064680
     8         2874464d-4*t13*t14*t17*t3*t45*wght+7.937005259840999d-1*t
     9         17*t45*wght-2.8740646802874464d-4*t14*t17*t3*t38*wght+7.9
     :         37005259840999d-1*t12*t17*t38*wght-5.748129360574893d-4*t
     ;         13*t14*t17*t27*t3*t30*wght+1.5874010519681997d+0*t17*t27*
     <         t30*wght)
            t92 = 9.999999999999998d-1*t2*t6*t8*(2.5146490204538424d-5*t
     1         14*t17*t20*t3*wght-6.944444444444446d-2*t12*t17*t20*wght)
     2         +t91+t90+t89+t88+t83+t82+t81+t80
            t93 = -4.023438432726146d-4*t15*t17*t3*t34*t57*wght
            t94 = -1.4370323401437232d-4*t17*t3*t32*t37*wght
            t95 = 2.2222222222222224d+0*t12*t17*t18*t33*t34*wght
            t96 = 5.0d-1*t12*t17*t37*t6*t8*wght
            t97 = 3.1498026247371824d-1*t20*t6*t8*(-6.349604207872799d+0
     1         *gammaaa*t12*t17*t61*wght+1.1496258721149785d-3*t14*t17*t
     2         3*t30*wght-3.1748021039363994d+0*t12*t17*t30*wght+1.14962
     3         58721149785d-3*t13*t14*t17*t27*t3*wght-3.1748021039363994
     4         d+0*t17*t27*wght)
            t98 = 4.19973683298291d-1*t15*t16*t18*(-3.1748021039363994d+
     1         0*gammaaa*t12*t17*t61*wght+t87+t86+t85+t84)
            t99 = 1.259921049894873d+0*t15*t16*t18*(-1.0582673679787998d
     1         +0*gammaaa*t12*t17*t61*wght+1.9160431201916311d-4*t14*t17
     2         *t3*t30*wght-5.291336839893999d-1*t12*t17*t30*wght+1.9160
     3         431201916311d-4*t13*t14*t17*t27*t3*wght-5.291336839893999
     4         d-1*t17*t27*wght)
            t100 = 6.299605249474365d-1*t2*t6*t8*(6.349604207872799d+0*g
     1         ammaaa*t12*t17*t30*t61*wght+6.349604207872799d+0*gammaaa*
     2         t17*t27*t61*wght+3.1748021039363994d+0*gammaaa*t12*t17*(-
     3         3.6210945894535324d-4*t14*t47*t48*t58-1.6722415986018052d
     4         -4*t13*t19*t46*t58+1.8105472947267662d-4*t14*t29*t45*t58-
     5         2.1123051771812276d-4*rhoa*t27*t28*t29*t58+5.280762942953
     6         069d-5*t13*t28*t58)*wght-5.748129360574893d-4*t14*t17*t3*
     7         t49*wght+1.5874010519681997d+0*t12*t17*t49*wght-5.7481293
     8         60574893d-4*t13*t14*t17*t3*t45*wght+1.5874010519681997d+0
     9         *t17*t45*wght-5.748129360574893d-4*t14*t17*t3*t38*wght+1.
     :         5874010519681997d+0*t12*t17*t38*wght-1.1496258721149785d-
     ;         3*t13*t14*t17*t27*t3*t30*wght+3.1748021039363994d+0*t17*t
     <         27*t30*wght)
            t101 = 1/t7**3.3333333333333337d+0
            t102 = 1/t7**2.333333333333334d+0
            t103 = -5.748129360574893d-4*t14*t17*t3*t59*wght
            t104 = 1/t3**3
            t105 = 1.5874010519681997d+0*t12*t17*t59*wght
            t106 = 6.299605249474365d-1*t2*t6*t8*(-1.0777742551077925d-4
     1         *t14*t17*t30*t58*wght-1.0777742551077925d-4*t13*t14*t17*t
     2         27*t58*wght+2.601815015882092d-8*t102*t13*t17*t30*wght+3.
     3         1748021039363994d+0*gammaaa*t12*t17*(-1.1315920592042289d
     4         -5*t104*t14*t27*t29-6.600953678691336d-6*rhoa*t104*t13*t2
     5         8)*wght+2.601815015882092d-8*t102*t17*t27*t29*wght+t105+t
     6         103)-6.504537539705231d-9*t13*t17*t20*t62*wght+2.26318411
     7         84084575d-5*t1*t16*t17*t57*t58*wght+2.694435637769481d-5*
     8         t17*t20*t32*t58*wght-5.463469177404018d-9*t1*t101*t13*t16
     9         *t17*wght
            t107 = 6.299605249474365d-1*t2*t6*t8*(-2.8740646802874464d-4
     1         *t14*t17*t3*t61*wght+7.937005259840999d-1*t12*t17*t61*wgh
     2         t-2.155548510215585d-4*t14*t17*t30*t58*wght-2.15554851021
     3         5585d-4*t13*t14*t17*t27*t58*wght+5.203630031764184d-8*t10
     4         2*t13*t17*t30*wght+3.1748021039363994d+0*gammaaa*t12*t17*
     5         (-2.2631841184084578d-5*t104*t14*t27*t29-1.32019073573826
     6         72d-5*rhoa*t104*t13*t28)*wght+5.203630031764184d-8*t102*t
     7         17*t27*t29*wght+t105+t103)-1.300907507941046d-8*t13*t17*t
     8         20*t62*wght+4.526368236816915d-5*t1*t16*t17*t57*t58*wght+
     9         5.388871275538962d-5*t17*t20*t32*t58*wght-1.0926938354808
     :         036d-8*t1*t101*t13*t16*t17*wght
            t108 = 1/t7**3.833333333333334d+0
            t109 = 1/gammaaa
            t110 = 4.878403154778923d-9*t109*t13*t17*t2*t62*wght-1.17767
     1         72845961996d-12*t108*t17*t2*t29*t58*wght+6.73608909442370
     2         2d-6*t104*t17*t2*t32*wght
            t111 = 9.756806309557846d-9*t109*t13*t17*t2*t62*wght-2.35535
     1         4569192399d-12*t108*t17*t2*t29*t58*wght+1.347217818884740
     2         6d-5*t104*t17*t2*t32*wght
            fnc(iq) = 4.000000000000001d+0*gammaaa*t1*t12*t2*exp(-7.2421
     1         89178907064d-4*t13*t14*t3)*t8*wght+fnc(iq)
            Amat(iq,D1_RA) = -1.0d+0*gammaaa*t12*t17*t20*t6*t8*wght-1.33
     1         33333333333333d+0*gammaaa*t12*t15*t16*t17*t18*wght+6.2996
     2         05249474365d-1*t2*t31*t6*t8+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 5.0d-1*t12*t17*t2*t6*t8*wght-1.43703234014
     1         37232d-4*t17*t2*t3*t32*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = 1.0d+0*t12*t17*t2*t6*t8*wght-2.87406468028
     1         74464d-4*t17*t2*t3*t32*wght+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 7.222222222222222d-1*gammaaa*t12*t17*t3
     1         7*t6*t8*wght+t56+t54+t52+t50+t35+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 1.2777777777777777d+0*gammaaa*t12*t17*t
     1         37*t6*t8*wght+t56+t54+t52+t50+t35+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t60+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 6.299605249474365d-1*t2*t6*t8*(3.17480
     1         21039363994d+0*gammaaa*t12*t17*t61*wght-5.748129360574893
     2         d-4*t14*t17*t3*t30*wght+1.5874010519681997d+0*t12*t17*t30
     3         *wght-5.748129360574893d-4*t13*t14*t17*t27*t3*wght+1.5874
     4         010519681997d+0*t17*t27*wght)-5.0d-1*t12*t17*t20*t6*t8*wg
     5         ht+1.2070315298178438d-4*t1*t16*t17*t3*t57*wght+1.4370323
     6         401437232d-4*t17*t20*t3*t32*wght-6.666666666666666d-1*t12
     7         *t15*t16*t17*t18*wght+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t60+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t63+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = 2.601815015882092d-8*t13*t17*t2*t62*w
     1         ght-1.0777742551077923d-4*t17*t2*t32*t58*wght+Cmat2(iq,D2
     2         _GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t63+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = 5.203630031764184d-8*t13*t17*t2*t62*w
     1         ght-2.1555485102155847d-4*t17*t2*t32*t58*wght+Cmat2(iq,D2
     2         _GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 9.999999999999998d-1*t2*t6*t8*(8.333
     1         333333333336d-1*gammaaa*t12*t17*t37*wght-5.55555555555555
     2         7d-1*gammaaa*t12*t17*t20*t30*wght-5.555555555555557d-1*ga
     3         mmaaa*t17*t20*t27*wght)-1.0833333333333335d+0*gammaaa*t12
     4         *t17*t6*t66*t8*wght+5.555555555555556d-1*gammaaa*t12*t15*
     5         t17*t18*t64*wght-8.749451735381063d-2*t31*t37*t6*t8+t79+t
     6         78+t75+t71+t70+t69+t68+t67+t65+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = -1.9166666666666665d+0*gammaaa*t12*t
     1         17*t6*t66*t8*wght-1.8518518518518515d-1*gammaaa*t12*t15*t
     2         17*t18*t64*wght+8.749451735381063d-2*t31*t37*t6*t8+t79+t7
     3         8+t75+t71+t70+t69+t68+t67+t65+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t92+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = 9.999999999999998d-1*t2*t6*t8*(5.02
     1         9298040907685d-5*t14*t17*t20*t3*wght-1.3888888888888892d-
     2         1*t12*t17*t20*wght)+t99+t98+t97+t96+t95+t94+t93+t100+Cmat
     3         3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t92+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 9.999999999999998d-1*t2*t6*t8*(6.94
     1         4444444444446d-2*t12*t17*t20*wght-2.5146490204538424d-5*t
     2         14*t17*t20*t3*wght)+t91+t90+t89+t88+t83+t82+t81+t80+Cmat3
     3         (iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = 9.999999999999998d-1*t2*t6*t8*(1.38
     1         88888888888892d-1*t12*t17*t20*wght-5.029298040907685d-5*t
     2         14*t17*t20*t3*wght)+t99+t98+t97+t96+t95+t94+t93+t100+Cmat
     3         3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = t106+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = t107+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t106+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = 6.299605249474365d-1*t2*t6*t8*(-1.
     1         1496258721149785d-3*t14*t17*t3*t61*wght+3.174802103936399
     2         4d+0*t12*t17*t61*wght-4.31109702043117d-4*t14*t17*t30*t58
     3         *wght-4.31109702043117d-4*t13*t14*t17*t27*t58*wght+1.0407
     4         260063528369d-7*t102*t13*t17*t30*wght+3.1748021039363994d
     5         +0*gammaaa*t12*t17*(-4.5263682368169156d-5*t104*t14*t27*t
     6         29-2.6403814714765345d-5*rhoa*t104*t13*t28)*wght+1.040726
     7         0063528369d-7*t102*t17*t27*t29*wght)-2.601815015882092d-8
     8         *t13*t17*t20*t62*wght+9.05273647363383d-5*t1*t16*t17*t57*
     9         t58*wght+1.0777742551077923d-4*t17*t20*t32*t58*wght-2.185
     :         387670961607d-8*t1*t101*t13*t16*t17*wght+Cmat3(iq,D3_RA_G
     ;         AB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = t107+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t106+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t110+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = t111+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t110+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = 1.9513612619115692d-8*t109*t13*t1
     1         7*t2*t62*wght-4.710709138384798d-12*t108*t17*t2*t29*t58*w
     2         ght+2.694435637769481d-5*t104*t17*t2*t32*wght+Cmat3(iq,D3
     3         _GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t111+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = 3.9027225238231383d-8*t109*t13*t1
     1         7*t2*t62*wght-9.421418276769596d-12*t108*t17*t2*t29*t58*w
     2         ght+5.388871275538962d-5*t104*t17*t2*t32*wght+Cmat3(iq,D3
     3         _GAB_GAB_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb+2.0d+0*gammaab+gammaaa
            t2 = rhob+rhoa
            t3 = 1/t2
            t4 = rhoa-rhob
            t5 = 1.0d+0-1.0d+0*t3*t4
            t6 = 1.0d+0*t3*t4+1.0d+0
            t7 = (3.149802624737183d-1*t6**1.6666666666666669d+0+3.14980
     1         2624737183d-1*t5**1.6666666666666669d+0)**5.0d-1
            t8 = 1/t7
            t9 = t1**5.0d-1
            t10 = 1/t2**6.666666666666666d-1
            t11 = abs(t2)
            t12 = 1/t11**3.333333333333333d-1
            t13 = 1.443307452126544d-2*t12+2.8435438314903855d-6*t10+2.5
     1         68d-3
            t14 = 1.7639938117590215d-2*t10*t12+5.411317332115466d+0*t12
     1         +1.816419932959077d-1*t10+1.0d+0
            t15 = 1/t14
            t16 = t13*t15+1.6670000000000001d-3
            t17 = 1/t16
            t18 = 1/t11**1.166666666666667d+0
            t19 = exp(-8.129082500000001d-4*t17*t18*t9)
            t20 = 1/t2**1.3333333333333333d+0
            t21 = 1/t11
            t22 = -2.6456684199469993d-1*t1*t16*t19*t20*t21*t8*wght
            t23 = t2**2
            t24 = 1/t23
            t25 = 1.0d+0*t24*t4
            t26 = -1.0d+0*t3
            t27 = t26+t25
            t28 = t5**6.666666666666666d-1
            t29 = -1.0d+0*t24*t4
            t30 = 1.0d+0*t3
            t31 = t30+t29
            t32 = t6**6.666666666666666d-1
            t33 = 5.249671041228638d-1*t31*t32+5.249671041228638d-1*t27*
     1         t28
            t34 = 1/t7**3
            t35 = -7.937005259840999d-1*t1*t12*t16*t19*t24*t8*wght
            t36 = 1/t2**1.6666666666666669d+0
            t37 = 1/t2**3.333333333333333d-1
            t38 = -1.8037724440384886d+0*t21*t37-1.1759958745060142d-2*t
     1         12*t36-1.2109466219727179d-1*t36-5.879979372530072d-3*t21
     2         *t3
            t39 = 1/t14**2
            t40 = -4.8110248404218126d-3*t21*t37-1.895695887660257d-6*t3
     1         6
            t41 = t15*t40-t13*t38*t39
            t42 = 1/t11**3.166666666666667d+0
            t43 = 1/t16**2
            t44 = 8.129082500000001d-4*t18*t41*t43*t9+9.483929583333336d
     1         -4*t17*t2*t42*t9
            t45 = 7.937005259840999d-1*t1*t16*t19*t44*wght+7.93700525984
     1         0999d-1*t1*t19*t41*wght
            t46 = t12*t3*t45*t8
            t47 = t30+t25
            t48 = t29+t26
            t49 = 5.249671041228638d-1*t32*t48+5.249671041228638d-1*t28*
     1         t47
            t50 = 1/t11**1.5000000000000002d+0
            t51 = 7.937005259840999d-1*t12*t16*t19*t3*t8*wght-3.22602852
     1         8009071d-4*t19*t3*t50*t8*t9*wght
            t52 = 1/t2**2.3333333333333334d+0
            t53 = 8.818894733156665d-1*t1*t16*t19*t21*t52*t8*wght
            t54 = t33**2
            t55 = 1/t7**5
            t56 = t2**3
            t57 = 1/t56
            t58 = 1.5874010519681997d+0*t1*t12*t16*t19*t57*t8*wght
            t59 = t44**2
            t60 = t38**2
            t61 = 1/t14**3
            t62 = 1/t2**2.6666666666666666d+0
            t63 = 1.9599931241766905d-2*t12*t62+2.01824436995453d-1*t62+
     1         1.5679944993413522d-2*t21*t24+2.4050299253846513d+0*t20*t
     2         21
            t64 = 3.1594931461004283d-6*t62+6.414699787229083d-3*t20*t21
            t65 = t15*t64-t13*t39*t63+2*t13*t60*t61-2*t38*t39*t40
            t66 = 1/t11**5.166666666666667d+0
            t67 = t41**2
            t68 = 1/t16**3
            t69 = -1.6258165000000002d-3*t18*t67*t68*t9-3.00324436805555
     1         7d-3*t17*t23*t66*t9+8.129082500000001d-4*t18*t43*t65*t9-1
     2         .8967859166666673d-3*t2*t41*t42*t43*t9+9.483929583333336d
     3         -4*t17*t42*t9
            t70 = 7.937005259840999d-1*t1*t16*t19*t69*wght+7.93700525984
     1         0999d-1*t1*t19*t65*wght+7.937005259840999d-1*t1*t16*t19*t
     2         59*wght+1.5874010519681997d+0*t1*t19*t41*t44*wght
            t71 = t12*t3*t70*t8
            t72 = -1.5874010519681997d+0*t1*t16*t19*t44*wght-1.587401051
     1         9681997d+0*t1*t19*t41*wght
            t73 = t12*t24*t72*t8
            t74 = -7.937005259840999d-1*t1*t16*t19*t44*wght-7.9370052598
     1         40999d-1*t1*t19*t41*wght
            t75 = 3.333333333333333d-1*t20*t21*t74*t8
            t76 = -2.6456684199469993d-1*t1*t16*t19*t44*wght-2.645668419
     1         9469993d-1*t1*t19*t41*wght
            t77 = t20*t21*t76*t8
            t78 = -3.9685026299204995d-1*t1*t19*t33*t41*wght
            t79 = t27**2
            t80 = 1/t5**3.333333333333333d-1
            t81 = -2.0d+0*t4*t57
            t82 = 2.0d+0*t24
            t83 = t82+t81
            t84 = t31**2
            t85 = 1/t6**3.333333333333333d-1
            t86 = 2.0d+0*t4*t57
            t87 = -2.0d+0*t24
            t88 = t87+t86
            t89 = 5.249671041228638d-1*t32*t88+3.4997806941524257d-1*t84
     1         *t85+5.249671041228638d-1*t28*t83+3.4997806941524257d-1*t
     2         79*t80
            t90 = -3.9685026299204995d-1*t1*t16*t19*t33*t44*wght
            t91 = -3.9685026299204995d-1*t1*t16*t19*t89*wght+t90+t78
            t92 = 1.3228342099734997d-1*t1*t16*t19*t49*wght+1.3228342099
     1         734997d-1*t1*t16*t19*t33*wght
            t93 = 3.9685026299204995d-1*t1*t16*t19*t49*wght+3.9685026299
     1         204995d-1*t1*t16*t19*t33*wght
            t94 = -5.0d-1*t12*t3*t34*t45*t49
            t95 = 3.4997806941524257d-1*t31*t48*t85+3.4997806941524257d-
     1         1*t27*t47*t80+1.0499342082457278d+0*t32*t4*t57-1.04993420
     2         82457278d+0*t28*t4*t57
            t96 = -3.9685026299204995d-1*t1*t16*t19*t95*wght+t90+t78
            t97 = t49**2
            t98 = t47**2
            t99 = t87+t81
            t100 = t48**2
            t101 = t86+t82
            t102 = 5.249671041228638d-1*t28*t99+3.4997806941524257d-1*t8
     1         0*t98+3.4997806941524257d-1*t100*t85+5.249671041228638d-1
     2         *t101*t32
            t103 = -3.9685026299204995d-1*t1*t16*t19*t44*t49*wght-3.9685
     1         026299204995d-1*t1*t19*t41*t49*wght-3.9685026299204995d-1
     2         *t1*t102*t16*t19*wght
            t104 = 1/t11**2.166666666666667d+0
            t105 = 1.0753428426696904d-4*t104*t19*t20*t8*t9*wght
            t106 = 3.226028528009071d-4*t19*t24*t50*t8*t9*wght
            t107 = -2.6456684199469993d-1*t16*t19*t20*t21*t8*wght
            t108 = -7.937005259840999d-1*t12*t16*t19*t24*t8*wght
            t109 = 1/t9
            t110 = 4.0645412500000005d-4*t109*t18*t41*t43+4.741964791666
     1         668d-4*t109*t17*t2*t42
            t111 = -3.226028528009071d-4*t18*t19*t44*t9*wght-3.226028528
     1         009071d-4*t17*t18*t19*t41*t9*wght+7.937005259840999d-1*t1
     2         6*t19*t44*wght+7.937005259840999d-1*t19*t41*wght+7.937005
     3         259840999d-1*t1*t110*t16*t19*wght
            t112 = t111*t12*t3*t8
            t113 = 1.6130142640045356d-4*t19*t3*t33*t34*t50*t9*wght-3.96
     1         85026299204995d-1*t12*t16*t19*t3*t33*t34*wght+t112+t108+t
     2         107+t106+t105
            t114 = 2.1506856853393808d-4*t104*t19*t20*t8*t9*wght
            t115 = 6.452057056018142d-4*t19*t24*t50*t8*t9*wght
            t116 = -5.291336839893999d-1*t16*t19*t20*t21*t8*wght
            t117 = -1.5874010519681997d+0*t12*t16*t19*t24*t8*wght
            t118 = 8.129082500000001d-4*t109*t18*t41*t43+9.4839295833333
     1         36d-4*t109*t17*t2*t42
            t119 = -6.452057056018142d-4*t18*t19*t44*t9*wght-6.452057056
     1         018142d-4*t17*t18*t19*t41*t9*wght+1.5874010519681997d+0*t
     2         16*t19*t44*wght+1.5874010519681997d+0*t19*t41*wght+7.9370
     3         05259840999d-1*t1*t118*t16*t19*wght
            t120 = t119*t12*t3*t8
            t121 = 1.6130142640045356d-4*t19*t3*t34*t49*t50*t9*wght-3.96
     1         85026299204995d-1*t12*t16*t19*t3*t34*t49*wght+t112+t108+t
     2         107+t106+t105
            t122 = 1/t11**2.666666666666667d+0
            t123 = 1.311232602576965d-7*t122*t17*t19*t3*t8*wght-4.839042
     1         7920136064d-4*t109*t19*t3*t50*t8*wght
            t124 = 2.62246520515393d-7*t122*t17*t19*t3*t8*wght-9.6780855
     1         84027213d-4*t109*t19*t3*t50*t8*wght
            t125 = -3.468765261708288d+0*t1*t16*t19*t21*t8*wght/t2**3.33
     1         33333333333337d+0
            t126 = 1/t7**7
            t127 = 1/t2**4
            t128 = -4.762203155904599d+0*t1*t12*t127*t16*t19*t8*wght
            t129 = t12*t24*t8*(-2.3811015779522995d+0*t1*t16*t19*t69*wgh
     1         t-2.3811015779522995d+0*t1*t19*t65*wght-2.381101577952299
     2         5d+0*t1*t16*t19*t59*wght-4.762203155904599d+0*t1*t19*t41*
     3         t44*wght)
            t130 = 3.333333333333333d-1*t20*t21*t8*(-1.5874010519681997d
     1         +0*t1*t16*t19*t69*wght-1.5874010519681997d+0*t1*t19*t65*w
     2         ght-1.5874010519681997d+0*t1*t16*t19*t59*wght-3.174802103
     3         9363994d+0*t1*t19*t41*t44*wght)
            t131 = t20*t21*t8*(-2.6456684199469993d-1*t1*t16*t19*t69*wgh
     1         t-2.6456684199469993d-1*t1*t19*t65*wght-2.645668419946999
     2         3d-1*t1*t16*t19*t59*wght-5.291336839893999d-1*t1*t19*t41*
     3         t44*wght)
            t132 = -7.777777777777778d-1*t21*t52*t74*t8
            t133 = t21*t52*t8*(8.818894733156665d-1*t1*t16*t19*t44*wght+
     1         8.818894733156665d-1*t1*t19*t41*wght)
            t134 = 2.6456684199469993d-1*t1*t19*t33*t41*wght
            t135 = 2.6456684199469993d-1*t1*t16*t19*t33*t44*wght
            t136 = 3.9685026299204995d-1*t1*t19*t33*t41*wght
            t137 = 3.9685026299204995d-1*t1*t16*t19*t89*wght
            t138 = 3.9685026299204995d-1*t1*t16*t19*t33*t44*wght
            t139 = 3.333333333333333d-1*(t138+t137+t136)*t20*t21*t34
            t140 = 1.1905507889761498d+0*t1*t19*t33*t41*wght
            t141 = 1.1905507889761498d+0*t1*t16*t19*t33*t44*wght
            t142 = 5.952753944880749d-1*t1*t19*t41*t54*wght
            t143 = 5.952753944880749d-1*t1*t16*t19*t44*t54*wght
            t144 = 1/t2**3.6666666666666664d+0
            t145 = -3*t38*t39*t64+6*t13*t38*t61*t63-3*t39*t40*t63+6*t40*
     1         t60*t61-t13*t39*(-5.357314539416287d-2*t21*t57-5.61173649
     2         2564187d+0*t21*t52-5.226648331137841d-2*t12*t144-5.381984
     3         986545413d-1*t144)+t15*(-1.496763283686786d-2*t21*t52-8.4
     4         25315056267809d-6*t144)-6*t13*t38**3/t14**4
            t146 = t12*t3*t8*(7.937005259840999d-1*t1*t16*t19*(5.6903577
     1         50000001d-3*t2*t42*t67*t68*t9-4.8774495000000007d-3*t18*t
     2         41*t65*t68*t9+9.009733104166671d-3*t23*t41*t43*t66*t9-9.0
     3         09733104166671d-3*t17*t2*t66*t9-2.8451788750000007d-3*t2*
     4         t42*t43*t65*t9+1.5516762568287046d-2*t17*t56*t9/t11**7.16
     5         6666666666667d+0-2.8451788750000007d-3*t41*t42*t43*t9+8.1
     6         29082500000001d-4*t145*t18*t43*t9+4.8774495000000007d-3*t
     7         18*t41**3*t9/t16**4)*wght+2.3811015779522995d+0*t1*t16*t1
     8         9*t44*t69*wght+2.3811015779522995d+0*t1*t19*t41*t69*wght+
     9         2.3811015779522995d+0*t1*t19*t44*t65*wght+2.3811015779522
     :         995d+0*t1*t19*t41*t59*wght+7.937005259840999d-1*t1*t16*t1
     ;         9*t44**3*wght+7.937005259840999d-1*t1*t145*t19*wght)
            t147 = -3.9685026299204995d-1*t1*t16*t19*t33*t59*wght
            t148 = -3.9685026299204995d-1*t1*t19*t33*t65*wght
            t149 = 1/t5**1.3333333333333333d+0
            t150 = 6.0d+0*t127*t4
            t151 = -6.0d+0*t57
            t152 = 1/t6**1.3333333333333333d+0
            t153 = -6.0d+0*t127*t4
            t154 = 6.0d+0*t57
            t155 = -3.9685026299204995d-1*t1*t16*t19*t33*t69*wght
            t156 = -7.937005259840999d-1*t1*t19*t33*t41*t44*wght
            t157 = -t33*t70
            t158 = 1.5874010519681997d+0*t1*t19*t41*wght
            t159 = 1.5874010519681997d+0*t1*t16*t19*t44*wght
            t160 = t12*t57*(-2*t72+t159+t158)*t8
            t161 = 3.333333333333333d-1*t21*t52*(-7*t76+t159+t158)*t8
            t162 = t33*t45
            t163 = -7.937005259840999d-1*t1*t16*t19*t49*wght
            t164 = -5.0d-1*t20*t21*t34*t49*t76
            t165 = 7.937005259840999d-1*t1*t16*t19*t95*wght
            t166 = -t49*t70
            t167 = -t49*t72
            t168 = -t49*t74
            t169 = 7.5d-1*t12*t3*t45*t55*t97
            t170 = -5.0d-1*t102*t12*t3*t34*t45
            t171 = 3.9685026299204995d-1*t1*t19*t41*t49*wght
            t172 = 3.9685026299204995d-1*t1*t102*t16*t19*wght
            t173 = 3.9685026299204995d-1*t1*t16*t19*t44*t49*wght
            t174 = t45*t49
            t175 = 5.0d-1*t12*(t174+t167)*t24*t34
            t176 = 1.6666666666666666d-1*(t174+t168)*t20*t21*t34
            t177 = -3.5844761422323007d-4*t104*t19*t52*t8*t9*wght
            t178 = -6.452057056018142d-4*t19*t50*t57*t8*t9*wght
            t179 = 8.818894733156665d-1*t16*t19*t21*t52*t8*wght
            t180 = 1.5874010519681997d+0*t12*t16*t19*t57*t8*wght
            t181 = -1.5874010519681997d+0*t19*t41*wght
            t182 = 6.452057056018142d-4*t17*t18*t19*t41*t9*wght
            t183 = 6.452057056018142d-4*t18*t19*t44*t9*wght
            t184 = -1.5874010519681997d+0*t16*t19*t44*wght
            t185 = t12*t24*t8*(-1.5874010519681997d+0*t1*t110*t16*t19*wg
     1         ht+t184+t183+t182+t181)
            t186 = 3.333333333333333d-1*t20*t21*t8*(3.226028528009071d-4
     1         *t18*t19*t44*t9*wght+3.226028528009071d-4*t17*t18*t19*t41
     2         *t9*wght-7.937005259840999d-1*t16*t19*t44*wght-7.93700525
     3         9840999d-1*t19*t41*wght-7.937005259840999d-1*t1*t110*t16*
     4         t19*wght)
            t187 = t20*t21*t8*(1.0753428426696904d-4*t18*t19*t44*t9*wght
     1         +1.0753428426696904d-4*t17*t18*t19*t41*t9*wght-2.64566841
     2         99469993d-1*t16*t19*t44*wght-2.6456684199469993d-1*t19*t4
     3         1*wght-2.6456684199469993d-1*t1*t110*t16*t19*wght)
            t188 = -3.9685026299204995d-1*t19*t33*t41*wght
            t189 = 1.6130142640045356d-4*t17*t18*t19*t33*t41*t9*wght
            t190 = 1.6130142640045356d-4*t18*t19*t33*t44*t9*wght
            t191 = -3.9685026299204995d-1*t1*t110*t16*t19*t33*wght
            t192 = -3.9685026299204995d-1*t16*t19*t33*t44*wght
            t193 = t12*t3*t8*(-3.226028528009071d-4*t18*t19*t69*t9*wght-
     1         3.226028528009071d-4*t17*t18*t19*t65*t9*wght-3.2260285280
     2         09071d-4*t18*t19*t59*t9*wght-6.452057056018142d-4*t17*t18
     3         *t19*t41*t44*t9*wght+7.937005259840999d-1*t16*t19*t69*wgh
     4         t+7.937005259840999d-1*t1*t16*t19*(-8.129082500000001d-4*
     5         t109*t18*t67*t68-1.5016221840277785d-3*t109*t17*t23*t66+4
     6         .0645412500000005d-4*t109*t18*t43*t65-9.483929583333336d-
     7         4*t109*t2*t41*t42*t43+4.741964791666668d-4*t109*t17*t42)*
     8         wght+7.937005259840999d-1*t19*t65*wght+7.937005259840999d
     9         -1*t16*t19*t59*wght+1.5874010519681997d+0*t19*t41*t44*wgh
     :         t+1.5874010519681997d+0*t1*t110*t16*t19*t44*wght+1.587401
     ;         0519681997d+0*t1*t110*t19*t41*wght)
            t194 = t12*t3*t34*(1.6130142640045356d-4*t18*t19*t89*t9*wght
     1         -3.9685026299204995d-1*t16*t19*t89*wght+t192+t191+t190+t1
     2         89+t188)-2.4195213960068032d-4*t19*t3*t50*t54*t55*t9*wght
     3         -3.226028528009071d-4*t19*t24*t33*t34*t50*t9*wght-1.07534
     4         28426696904d-4*t104*t19*t20*t33*t34*t9*wght+5.95275394488
     5         0749d-1*t12*t16*t19*t3*t54*t55*wght+7.937005259840999d-1*
     6         t12*t16*t19*t24*t33*t34*wght+2.6456684199469993d-1*t16*t1
     7         9*t20*t21*t33*t34*wght-5.0d-1*t111*t12*t3*t33*t34+t193+t1
     8         87+t186+t185+t180+t179+t178+t177
            t195 = -7.168952284464601d-4*t104*t19*t52*t8*t9*wght
            t196 = -1.2904114112036286d-3*t19*t50*t57*t8*t9*wght
            t197 = 1.763778946631333d+0*t16*t19*t21*t52*t8*wght
            t198 = 3.1748021039363994d+0*t12*t16*t19*t57*t8*wght
            t199 = t12*t24*t8*(1.2904114112036286d-3*t18*t19*t44*t9*wght
     1         +1.2904114112036286d-3*t17*t18*t19*t41*t9*wght-3.17480210
     2         39363994d+0*t16*t19*t44*wght-3.1748021039363994d+0*t19*t4
     3         1*wght-1.5874010519681997d+0*t1*t118*t16*t19*wght)
            t200 = 3.333333333333333d-1*t20*t21*t8*(-7.937005259840999d-
     1         1*t1*t118*t16*t19*wght+t184+t183+t182+t181)
            t201 = t20*t21*t8*(2.1506856853393808d-4*t18*t19*t44*t9*wght
     1         +2.1506856853393808d-4*t17*t18*t19*t41*t9*wght-5.29133683
     2         9893999d-1*t16*t19*t44*wght-5.291336839893999d-1*t19*t41*
     3         wght-2.6456684199469993d-1*t1*t118*t16*t19*wght)
            t202 = -7.937005259840999d-1*t19*t33*t41*wght
            t203 = 3.226028528009071d-4*t17*t18*t19*t33*t41*t9*wght
            t204 = 3.226028528009071d-4*t18*t19*t33*t44*t9*wght
            t205 = -3.9685026299204995d-1*t1*t118*t16*t19*t33*wght
            t206 = -7.937005259840999d-1*t16*t19*t33*t44*wght
            t207 = t12*t3*t8*(-6.452057056018142d-4*t18*t19*t69*t9*wght-
     1         6.452057056018142d-4*t17*t18*t19*t65*t9*wght-6.4520570560
     2         18142d-4*t18*t19*t59*t9*wght-1.2904114112036286d-3*t17*t1
     3         8*t19*t41*t44*t9*wght+1.5874010519681997d+0*t16*t19*t69*w
     4         ght+7.937005259840999d-1*t1*t16*t19*(-1.6258165000000002d
     5         -3*t109*t18*t67*t68-3.003244368055557d-3*t109*t17*t23*t66
     6         +8.129082500000001d-4*t109*t18*t43*t65-1.8967859166666673
     7         d-3*t109*t2*t41*t42*t43+9.483929583333336d-4*t109*t17*t42
     8         )*wght+1.5874010519681997d+0*t19*t65*wght+1.5874010519681
     9         997d+0*t16*t19*t59*wght+3.1748021039363994d+0*t19*t41*t44
     :         *wght+1.5874010519681997d+0*t1*t118*t16*t19*t44*wght+1.58
     ;         74010519681997d+0*t1*t118*t19*t41*wght)
            t208 = -5.0d-1*t111*t12*t3*t34*t49
            t209 = t12*t3*t34*(1.6130142640045356d-4*t18*t19*t9*t95*wght
     1         -3.9685026299204995d-1*t16*t19*t95*wght+t192+t191+t190+t1
     2         89+t188)+t20*t21*t34*(-5.376714213348451d-5*t18*t19*t49*t
     3         9*wght-5.376714213348451d-5*t18*t19*t33*t9*wght+1.3228342
     4         099734997d-1*t16*t19*t49*wght+1.3228342099734997d-1*t16*t
     5         19*t33*wght)+t12*t24*t34*(-1.6130142640045356d-4*t18*t19*
     6         t49*t9*wght-1.6130142640045356d-4*t18*t19*t33*t9*wght+3.9
     7         685026299204995d-1*t16*t19*t49*wght+3.9685026299204995d-1
     8         *t16*t19*t33*wght)-2.4195213960068032d-4*t19*t3*t33*t49*t
     9         50*t55*t9*wght+5.952753944880749d-1*t12*t16*t19*t3*t33*t4
     :         9*t55*wght+t208+t193+t187+t186+t185+t180+t179+t178+t177
            t210 = -5.0d-1*t119*t12*t3*t34*t49
            t211 = t12*t3*t34*(1.6130142640045356d-4*t18*t19*t44*t49*t9*
     1         wght+1.6130142640045356d-4*t17*t18*t19*t41*t49*t9*wght+1.
     2         6130142640045356d-4*t102*t18*t19*t9*wght-3.96850262992049
     3         95d-1*t16*t19*t44*t49*wght-3.9685026299204995d-1*t19*t41*
     4         t49*wght-3.9685026299204995d-1*t1*t110*t16*t19*t49*wght-3
     5         .9685026299204995d-1*t102*t16*t19*wght)-2.419521396006803
     6         2d-4*t19*t3*t50*t55*t9*t97*wght+5.952753944880749d-1*t12*
     7         t16*t19*t3*t55*t97*wght-3.226028528009071d-4*t19*t24*t34*
     8         t49*t50*t9*wght-1.0753428426696904d-4*t104*t19*t20*t34*t4
     9         9*t9*wght+7.937005259840999d-1*t12*t16*t19*t24*t34*t49*wg
     :         ht+2.6456684199469993d-1*t16*t19*t20*t21*t34*t49*wght+t20
     ;         8+t193+t187+t186+t185+t180+t179+t178+t177
            t212 = 1/t11**3.3333333333333337d+0
            t213 = -4.370775341923217d-8*t17*t19*t20*t212*t8*wght
            t214 = -1.311232602576965d-7*t122*t17*t19*t24*t8*wght
            t215 = 1.6130142640045356d-4*t104*t109*t19*t20*t8*wght
            t216 = 4.8390427920136064d-4*t109*t19*t24*t50*t8*wght
            t217 = 1/t11**2.333333333333334d+0
            t218 = -6.452057056018142d-4*t110*t18*t19*t9*wght
            t219 = 1/t9**3
            t220 = 1.5874010519681997d+0*t110*t16*t19*wght
            t221 = t12*t3*t8*(1.311232602576965d-7*t17*t19*t217*t44*wght
     1         -4.8390427920136064d-4*t109*t18*t19*t44*wght+7.9370052598
     2         40999d-1*t1*t16*t19*(-2.0322706250000003d-4*t18*t219*t41*
     3         t43-2.370982395833334d-4*t17*t2*t219*t42)*wght+1.31123260
     4         2576965d-7*t19*t217*t41*t43*wght-4.8390427920136064d-4*t1
     5         09*t17*t18*t19*t41*wght+t220+t218)
            t222 = 2.4195213960068032d-4*t109*t19*t3*t33*t34*t50*wght-6.
     1         556163012884825d-8*t122*t17*t19*t3*t33*t34*wght+t221+t216
     2         +t215+t214+t213
            t223 = -8.741550683846434d-8*t17*t19*t20*t212*t8*wght
            t224 = -2.62246520515393d-7*t122*t17*t19*t24*t8*wght
            t225 = 3.226028528009071d-4*t104*t109*t19*t20*t8*wght
            t226 = 9.678085584027213d-4*t109*t19*t24*t50*t8*wght
            t227 = t12*t3*t8*(-3.226028528009071d-4*t118*t18*t19*t9*wght
     1         +2.62246520515393d-7*t17*t19*t217*t44*wght-9.678085584027
     2         213d-4*t109*t18*t19*t44*wght+7.937005259840999d-1*t1*t16*
     3         t19*(-4.0645412500000005d-4*t18*t219*t41*t43-4.7419647916
     4         66668d-4*t17*t2*t219*t42)*wght+2.62246520515393d-7*t19*t2
     5         17*t41*t43*wght-9.678085584027213d-4*t109*t17*t18*t19*t41
     6         *wght+7.937005259840999d-1*t118*t16*t19*wght+t220+t218)
            t228 = 4.8390427920136064d-4*t109*t19*t3*t33*t34*t50*wght-1.
     1         311232602576965d-7*t122*t17*t19*t3*t33*t34*wght+t227+t226
     2         +t225+t224+t223
            t229 = -1.748310136769287d-7*t17*t19*t20*t212*t8*wght
            t230 = -5.24493041030786d-7*t122*t17*t19*t24*t8*wght
            t231 = 6.452057056018142d-4*t104*t109*t19*t20*t8*wght
            t232 = 1.9356171168054426d-3*t109*t19*t24*t50*t8*wght
            t233 = t12*t3*t8*(-1.2904114112036286d-3*t118*t18*t19*t9*wgh
     1         t+5.24493041030786d-7*t17*t19*t217*t44*wght-1.93561711680
     2         54426d-3*t109*t18*t19*t44*wght+7.937005259840999d-1*t1*t1
     3         6*t19*(-8.129082500000001d-4*t18*t219*t41*t43-9.483929583
     4         333336d-4*t17*t2*t219*t42)*wght+5.24493041030786d-7*t19*t
     5         217*t41*t43*wght-1.9356171168054426d-3*t109*t17*t18*t19*t
     6         41*wght+3.1748021039363994d+0*t118*t16*t19*wght)
            t234 = 2.4195213960068032d-4*t109*t19*t3*t34*t49*t50*wght-6.
     1         556163012884825d-8*t122*t17*t19*t3*t34*t49*wght+t221+t216
     2         +t215+t214+t213
            t235 = 4.8390427920136064d-4*t109*t19*t3*t34*t49*t50*wght-1.
     1         311232602576965d-7*t122*t17*t19*t3*t34*t49*wght+t227+t226
     2         +t225+t224+t223
            t236 = 1/t11**3.833333333333334d+0
            t237 = 1/t1
            t238 = 2.4195213960068032d-4*t19*t219*t3*t50*t8*wght-5.32955
     1         9001518931d-11*t109*t19*t236*t3*t43*t8*wght+1.96684890386
     2         54477d-7*t122*t17*t19*t237*t3*t8*wght
            t239 = 4.8390427920136064d-4*t19*t219*t3*t50*t8*wght-1.06591
     1         18003037862d-10*t109*t19*t236*t3*t43*t8*wght+3.9336978077
     2         308954d-7*t122*t17*t19*t237*t3*t8*wght
            t240 = 9.678085584027213d-4*t19*t219*t3*t50*t8*wght-2.131823
     1         6006075725d-10*t109*t19*t236*t3*t43*t8*wght+7.86739561546
     2         1791d-7*t122*t17*t19*t237*t3*t8*wght
            fnc(iq) = 7.937005259840999d-1*t1*t12*t16*t19*t3*t8*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = -3.9685026299204995d-1*t1*t12*t16*t19*t3*t3
     1         3*t34*wght+t46+t35+t22+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = -3.9685026299204995d-1*t1*t12*t16*t19*t3*t3
     1         4*t49*wght+t46+t35+t22+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = t51+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = -6.452057056018142d-4*t19*t3*t50*t8*t9*wgh
     1         t+1.5874010519681997d+0*t12*t16*t19*t3*t8*wght+Cmat(iq,D1
     2         _GAB)
            Cmat(iq,D1_GBB) = t51+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = 5.952753944880749d-1*t1*t12*t16*t19*t3*
     1         t54*t55*wght+7.937005259840999d-1*t1*t12*t16*t19*t24*t33*
     2         t34*wght+2.6456684199469993d-1*t1*t16*t19*t20*t21*t33*t34
     3         *wght+t12*t3*t34*t91+t77+t75+t73+t71+t58+t53-5.0d-1*t12*t
     4         3*t33*t34*t45+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 5.952753944880749d-1*t1*t12*t16*t19*t3*
     1         t33*t49*t55*wght+t12*t3*t34*t96+t94+t12*t24*t34*t93+t20*t
     2         21*t34*t92+t77+t75+t73+t71+t58+t53+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = 5.952753944880749d-1*t1*t12*t16*t19*t3*
     1         t55*t97*wght+7.937005259840999d-1*t1*t12*t16*t19*t24*t34*
     2         t49*wght+2.6456684199469993d-1*t1*t16*t19*t20*t21*t34*t49
     3         *wght+t94+t77+t75+t73+t71+t58+t53+t103*t12*t3*t34+Amat2(i
     4         q,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t113+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = 3.226028528009071d-4*t19*t3*t33*t34*t5
     1         0*t9*wght-7.937005259840999d-1*t12*t16*t19*t3*t33*t34*wgh
     2         t+t120+t117+t116+t115+t114+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t113+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = t121+Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = 3.226028528009071d-4*t19*t3*t34*t49*t5
     1         0*t9*wght-7.937005259840999d-1*t12*t16*t19*t3*t34*t49*wgh
     2         t+t120+t117+t116+t115+t114+Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t121+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t123+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t124+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t123+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = -1.9356171168054426d-3*t109*t19*t3*t5
     1         0*t8*wght+5.24493041030786d-7*t122*t17*t19*t3*t8*wght+Cma
     2         t2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = t124+Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t123+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = t12*t3*t34*(-7.937005259840999d-1*t1
     1         *t16*t19*t44*t89*wght-7.937005259840999d-1*t1*t19*t41*t89
     2         *wght-3.9685026299204995d-1*t1*t16*t19*(1.049934208245727
     3         8d+0*t31*t85*t88+1.0499342082457278d+0*t27*t80*t83+5.2496
     4         71041228638d-1*(t154+t153)*t32-1.1665935647174752d-1*t152
     5         *t31**3+5.249671041228638d-1*(t151+t150)*t28-1.1665935647
     6         174752d-1*t149*t27**3)*wght+t157+t156+t155+t148+t147)+t12
     7         *t3*t55*(1.1905507889761498d+0*t1*t16*t19*t33*t89*wght+t1
     8         43+t142)+t12*t24*t34*(1.1905507889761498d+0*t1*t16*t19*t8
     9         9*wght+t141+t140)+t20*t21*t34*(2.6456684199469993d-1*t1*t
     :         16*t19*t89*wght+t135+t134)-2.3811015779522995d+0*t1*t12*t
     ;         16*t19*t33*t34*t57*wght-1.7858261834642247d+0*t1*t12*t16*
     <         t19*t24*t54*t55*wght-5.952753944880749d-1*t1*t16*t19*t20*
     =         t21*t54*t55*wght-1.3228342099734997d+0*t1*t16*t19*t21*t33
     >         *t34*t52*wght-1.4881884862201872d+0*t1*t12*t126*t16*t19*t
     ?         3*t33**3*wght-1.5d+0*t12*t3*t33*t55*t91-5.0d-1*t12*t3*t34
     @         *t45*t89-5.0d-1*t20*t21*t33*t34*t76+1.6666666666666666d-1
     1         *t20*t21*t34*(t162-t33*t74)+5.0d-1*t12*t24*t34*(t162-t33*
     2         t72)+7.5d-1*t12*t3*t45*t54*t55+t161+t160+t146+t139+t133+t
     3         132+t131+t130+t129+t128+t125+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t12*t3*t34*(-3.9685026299204995d-1*t
     1         1*t16*t19*t44*t95*wght-3.9685026299204995d-1*t1*t19*t41*t
     2         95*wght-3.9685026299204995d-1*t1*t16*t19*t44*t89*wght-3.9
     3         685026299204995d-1*t1*t19*t41*t89*wght-3.9685026299204995
     4         d-1*t1*t16*t19*(3.4997806941524257d-1*t48*t85*t88+1.39991
     5         22776609702d+0*t31*t4*t57*t85-1.1665935647174752d-1*t152*
     6         t48*t84+3.4997806941524257d-1*t47*t80*t83-1.3999122776609
     7         702d+0*t27*t4*t57*t80-1.1665935647174752d-1*t149*t47*t79+
     8         5.249671041228638d-1*t32*(2.0d+0*t57+t153)+5.249671041228
     9         638d-1*t28*(t150-2.0d+0*t57))*wght+t156+t155+t148+t147)+t
     :         12*t3*t55*(1.1905507889761498d+0*t1*t16*t19*t33*t95*wght+
     ;         t143+t142)+t20*t21*t34*(2.6456684199469993d-1*t1*t16*t19*
     <         t95*wght+t135+t134)+t20*t21*t55*(-1.9842513149602495d-1*t
     =         1*t16*t19*t54*wght-3.968502629920499d-1*t1*t16*t19*t33*t4
     >         9*wght)+t12*t24*t55*(-5.952753944880749d-1*t1*t16*t19*t54
     ?         *wght-1.1905507889761498d+0*t1*t16*t19*t33*t49*wght)+t21*
     @         t34*t52*(-4.4094473665783324d-1*t1*t16*t19*t49*wght-8.818
     1         894733156665d-1*t1*t16*t19*t33*wght)+t12*t34*t57*(t163-1.
     2         5874010519681997d+0*t1*t16*t19*t33*wght)-1.48818848622018
     3         72d+0*t1*t12*t126*t16*t19*t3*t49*t54*wght+5.0d-1*t12*t3*t
     4         34*(-t45*t95+t166+t157)-1.5d+0*t12*t3*t49*t55*t91+7.5d-1*
     5         t12*t3*t33*t45*t49*t55+5.0d-1*t12*(t167+t162)*t24*t34+t12
     6         *(t165+t141+t140+t137)*t24*t34+1.6666666666666666d-1*(t16
     7         8+t162)*t20*t21*t34+t164+t161+t160+t146+t139+t133+t132+t1
     8         31+t130+t129+t128+t125+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = t12*t3*t34*(-3.9685026299204995d-1*t
     1         1*t16*t19*(3.4997806941524257d-1*t27*t80*t99-1.1665935647
     2         174752d-1*t149*t27*t98+1.3999122776609702d+0*t4*t48*t57*t
     3         85+3.4997806941524257d-1*t101*t31*t85-1.3999122776609702d
     4         +0*t4*t47*t57*t80-1.0499342082457278d+0*t32*t57+1.0499342
     5         082457278d+0*t28*t57-3.149802624737183d+0*t127*t32*t4+3.1
     6         49802624737183d+0*t127*t28*t4-1.1665935647174752d-1*t100*
     7         t152*t31)*wght-7.937005259840999d-1*t1*t16*t19*t44*t95*wg
     8         ht-7.937005259840999d-1*t1*t19*t41*t95*wght+t166+t156+t15
     9         5+t148+t147)+t12*t3*t55*(5.952753944880749d-1*t1*t16*t19*
     :         t49*t95*wght+5.952753944880749d-1*t1*t16*t19*t33*t44*t49*
     ;         wght+5.952753944880749d-1*t1*t19*t33*t41*t49*wght+5.95275
     <         3944880749d-1*t1*t102*t16*t19*t33*wght)+3.333333333333333
     =         d-1*t20*t21*t34*(3.9685026299204995d-1*t1*t16*t19*t95*wgh
     >         t+t138+t136)+t20*t21*t34*(1.3228342099734997d-1*t1*t16*t1
     ?         9*t95*wght+1.3228342099734997d-1*t1*t16*t19*t44*t49*wght+
     @         1.3228342099734997d-1*t1*t19*t41*t49*wght+1.3228342099734
     1         997d-1*t1*t16*t19*t33*t44*wght+1.3228342099734997d-1*t1*t
     2         19*t33*t41*wght+1.3228342099734997d-1*t1*t102*t16*t19*wgh
     3         t)+3.333333333333333d-1*t21*t34*t52*(-3.9685026299204995d
     4         -1*t1*t16*t19*t49*wght-3.9685026299204995d-1*t1*t16*t19*t
     5         33*wght-7*t92)+t12*t24*t34*(7.937005259840999d-1*t1*t16*t
     6         19*t33*t44*wght+7.937005259840999d-1*t1*t19*t33*t41*wght+
     7         t173+t172+t171+t165)-1.4881884862201872d+0*t1*t12*t126*t1
     8         6*t19*t3*t33*t97*wght-5.952753944880749d-1*t1*t12*t16*t19
     9         *t24*t33*t49*t55*wght-1.9842513149602495d-1*t1*t16*t19*t2
     :         0*t21*t33*t49*t55*wght-4.4094473665783324d-1*t1*t16*t19*t
     ;         21*t34*t49*t52*wght-1.5d+0*t12*t3*t49*t55*t96-1.5d+0*t12*
     <         t24*t49*t55*t93+t12*t34*t57*(t163-2*t93)-1.5d+0*t20*t21*t
     =         49*t55*t92+t176+t175+t170+t169+t164+t161+t160+t146+t133+t
     >         132+t131+t130+t129+t128+t125+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = t12*t3*t34*(-3.9685026299204995d-1*t
     1         1*t16*t19*(1.0499342082457278d+0*t47*t80*t99+1.0499342082
     2         457278d+0*t101*t48*t85-1.1665935647174752d-1*t152*t48**3-
     3         1.1665935647174752d-1*t149*t47**3+5.249671041228638d-1*(t
     4         153+t151)*t32+5.249671041228638d-1*(t154+t150)*t28)*wght-
     5         3.9685026299204995d-1*t1*t16*t19*t49*t69*wght-3.968502629
     6         9204995d-1*t1*t19*t49*t65*wght-3.9685026299204995d-1*t1*t
     7         16*t19*t49*t59*wght-7.937005259840999d-1*t1*t19*t41*t44*t
     8         49*wght-7.937005259840999d-1*t1*t102*t16*t19*t44*wght-7.9
     9         37005259840999d-1*t1*t102*t19*t41*wght+t166)+t12*t3*t55*(
     :         5.952753944880749d-1*t1*t16*t19*t44*t97*wght+5.9527539448
     ;         80749d-1*t1*t19*t41*t97*wght+1.1905507889761498d+0*t1*t10
     <         2*t16*t19*t49*wght)+t12*t24*t34*(1.1905507889761498d+0*t1
     =         *t16*t19*t44*t49*wght+1.1905507889761498d+0*t1*t19*t41*t4
     >         9*wght+1.1905507889761498d+0*t1*t102*t16*t19*wght)+t20*t2
     ?         1*t34*(2.6456684199469993d-1*t1*t16*t19*t44*t49*wght+2.64
     @         56684199469993d-1*t1*t19*t41*t49*wght+2.6456684199469993d
     1         -1*t1*t102*t16*t19*wght)-1.7858261834642247d+0*t1*t12*t16
     2         *t19*t24*t55*t97*wght-5.952753944880749d-1*t1*t16*t19*t20
     3         *t21*t55*t97*wght-2.3811015779522995d+0*t1*t12*t16*t19*t3
     4         4*t49*t57*wght-1.3228342099734997d+0*t1*t16*t19*t21*t34*t
     5         49*t52*wght-1.4881884862201872d+0*t1*t12*t126*t16*t19*t3*
     6         t49**3*wght-1.5d+0*t103*t12*t3*t49*t55+3.333333333333333d
     7         -1*(t173+t172+t171)*t20*t21*t34+t176+t175+t170+t169+t164+
     8         t161+t160+t146+t133+t132+t131+t130+t129+t128+t125+Amat3(i
     9         q,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t194+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = t12*t3*t34*(3.226028528009071d-4*t1
     1         8*t19*t89*t9*wght-7.937005259840999d-1*t16*t19*t89*wght+t
     2         206+t205+t204+t203+t202)-4.8390427920136064d-4*t19*t3*t50
     3         *t54*t55*t9*wght-6.452057056018142d-4*t19*t24*t33*t34*t50
     4         *t9*wght-2.1506856853393808d-4*t104*t19*t20*t33*t34*t9*wg
     5         ht+1.1905507889761498d+0*t12*t16*t19*t3*t54*t55*wght+1.58
     6         74010519681997d+0*t12*t16*t19*t24*t33*t34*wght+5.29133683
     7         9893999d-1*t16*t19*t20*t21*t33*t34*wght-5.0d-1*t119*t12*t
     8         3*t33*t34+t207+t201+t200+t199+t198+t197+t196+t195+Cmat3(i
     9         q,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t194+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t209+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = t12*t3*t34*(3.226028528009071d-4*t1
     1         8*t19*t9*t95*wght-7.937005259840999d-1*t16*t19*t95*wght+t
     2         206+t205+t204+t203+t202)+t20*t21*t34*(-1.0753428426696904
     3         d-4*t18*t19*t49*t9*wght-1.0753428426696904d-4*t18*t19*t33
     4         *t9*wght+2.6456684199469993d-1*t16*t19*t49*wght+2.6456684
     5         199469993d-1*t16*t19*t33*wght)+t12*t24*t34*(-3.2260285280
     6         09071d-4*t18*t19*t49*t9*wght-3.226028528009071d-4*t18*t19
     7         *t33*t9*wght+7.937005259840999d-1*t16*t19*t49*wght+7.9370
     8         05259840999d-1*t16*t19*t33*wght)-4.8390427920136064d-4*t1
     9         9*t3*t33*t49*t50*t55*t9*wght+1.1905507889761498d+0*t12*t1
     :         6*t19*t3*t33*t49*t55*wght+t210+t207+t201+t200+t199+t198+t
     ;         197+t196+t195+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = t209+Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = t211+Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = t12*t3*t34*(3.226028528009071d-4*t1
     1         8*t19*t44*t49*t9*wght+3.226028528009071d-4*t17*t18*t19*t4
     2         1*t49*t9*wght+3.226028528009071d-4*t102*t18*t19*t9*wght-7
     3         .937005259840999d-1*t16*t19*t44*t49*wght-7.93700525984099
     4         9d-1*t19*t41*t49*wght-3.9685026299204995d-1*t1*t118*t16*t
     5         19*t49*wght-7.937005259840999d-1*t102*t16*t19*wght)-4.839
     6         0427920136064d-4*t19*t3*t50*t55*t9*t97*wght+1.19055078897
     7         61498d+0*t12*t16*t19*t3*t55*t97*wght-6.452057056018142d-4
     8         *t19*t24*t34*t49*t50*t9*wght-2.1506856853393808d-4*t104*t
     9         19*t20*t34*t49*t9*wght+1.5874010519681997d+0*t12*t16*t19*
     :         t24*t34*t49*wght+5.291336839893999d-1*t16*t19*t20*t21*t34
     ;         *t49*wght+t210+t207+t201+t200+t199+t198+t197+t196+t195+Cm
     <         at3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = t211+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = t222+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = t228+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t222+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = 9.678085584027213d-4*t109*t19*t3*t
     1         33*t34*t50*wght-2.62246520515393d-7*t122*t17*t19*t3*t33*t
     2         34*wght+t233+t232+t231+t230+t229+Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = t228+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t222+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = t234+Cmat3(iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = t235+Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = t234+Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = 9.678085584027213d-4*t109*t19*t3*t
     1         34*t49*t50*wght-2.62246520515393d-7*t122*t17*t19*t3*t34*t
     2         49*wght+t233+t232+t231+t230+t229+Cmat3(iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = t235+Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t234+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t238+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = t239+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t238+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = t240+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t239+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = t238+Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = 1.9356171168054426d-3*t19*t219*t3
     1         *t50*t8*wght-4.263647201215145d-10*t109*t19*t236*t3*t43*t
     2         8*wght+1.573479123092358d-6*t122*t17*t19*t237*t3*t8*wght+
     3         Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = t240+Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = t239+Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t238+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa
            t2 = gammaaa**5.0d-1
            t3 = 1/rhoa**6.666666666666666d-1
            t4 = abs(rhoa)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/rhoa**1.3333333333333333d+0
            t14 = 1/t4
            t15 = rhoa**2
            t16 = 1/t15
            t17 = 1/rhoa**1.6666666666666669d+0
            t18 = 1/rhoa**3.333333333333333d-1
            t19 = -1.1759958745060142d-2*t17*t5-1.8037724440384886d+0*t1
     1         4*t18-1.2109466219727179d-1*t17-5.879979372530072d-3*t1*t
     2         14
            t20 = 1/t7**2
            t21 = -4.8110248404218126d-3*t14*t18-1.895695887660257d-6*t1
     1         7
            t22 = t21*t8-t19*t20*t6
            t23 = 1/t4**3.166666666666667d+0
            t24 = 1/t9**2
            t25 = 8.129082500000001d-4*t11*t2*t22*t24+9.483929583333336d
     1         -4*rhoa*t10*t2*t23
            t26 = 1/t4**1.5000000000000002d+0
            t27 = 1/rhoa**2.3333333333333334d+0
            t28 = rhoa**3
            t29 = 1/t28
            t30 = t25**2
            t31 = t19**2
            t32 = 1/t7**3
            t33 = 1/rhoa**2.6666666666666666d+0
            t34 = 1.9599931241766905d-2*t33*t5+2.01824436995453d-1*t33+1
     1         .5679944993413522d-2*t14*t16+2.4050299253846513d+0*t13*t1
     2         4
            t35 = 3.1594931461004283d-6*t33+6.414699787229083d-3*t13*t14
            t36 = t35*t8-t20*t34*t6+2*t31*t32*t6-2*t19*t20*t21
            t37 = 1/t4**5.166666666666667d+0
            t38 = t22**2
            t39 = 1/t9**3
            t40 = -1.6258165000000002d-3*t11*t2*t38*t39-3.00324436805555
     1         7d-3*t10*t15*t2*t37+8.129082500000001d-4*t11*t2*t24*t36-1
     2         .8967859166666673d-3*rhoa*t2*t22*t23*t24+9.48392958333333
     3         6d-4*t10*t2*t23
            t41 = -1.5874010519681997d+0*gammaaa*t12*t25*t9*wght-1.58740
     1         10519681997d+0*gammaaa*t12*t22*wght
            t42 = -7.937005259840999d-1*gammaaa*t12*t25*t9*wght-7.937005
     1         259840999d-1*gammaaa*t12*t22*wght
            t43 = -2.6456684199469993d-1*gammaaa*t12*t25*t9*wght-2.64566
     1         84199469993d-1*gammaaa*t12*t22*wght
            t44 = 1/t4**2.166666666666667d+0
            t45 = 1/t2
            t46 = 4.0645412500000005d-4*t11*t22*t24*t45+4.74196479166666
     1         8d-4*rhoa*t10*t23*t45
            t47 = 1/t4**2.666666666666667d+0
            t48 = 1/rhoa**3.6666666666666664d+0
            t49 = (-8.425315056267809d-6*t48-1.496763283686786d-2*t14*t2
     1         7)*t8-6*t19**3*t6/t7**4-t20*(-5.226648331137841d-2*t48*t5
     2         -5.381984986545413d-1*t48-5.357314539416287d-2*t14*t29-5.
     3         611736492564187d+0*t14*t27)*t6+6*t19*t32*t34*t6-3*t19*t20
     4         *t35-3*t20*t21*t34+6*t21*t31*t32
            t50 = 1.5874010519681997d+0*gammaaa*t12*t22*wght
            t51 = 1.5874010519681997d+0*gammaaa*t12*t25*t9*wght
            t52 = 1/t4**2.333333333333334d+0
            t53 = 1/t2**3
            fnc(iq) = 7.937005259840999d-1*gammaaa*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = t1*t5*(7.937005259840999d-1*gammaaa*t12*t25
     1         *t9*wght+7.937005259840999d-1*gammaaa*t12*t22*wght)-7.937
     2         005259840999d-1*gammaaa*t12*t16*t5*t9*wght-2.645668419946
     3         9993d-1*gammaaa*t12*t13*t14*t9*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*t26*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t1*t5*(7.937005259840999d-1*gammaaa*t12
     1         *t40*t9*wght+7.937005259840999d-1*gammaaa*t12*t30*t9*wght
     2         +7.937005259840999d-1*gammaaa*t12*t36*wght+1.587401051968
     3         1997d+0*gammaaa*t12*t22*t25*wght)+1.5874010519681997d+0*g
     4         ammaaa*t12*t29*t5*t9*wght+8.818894733156665d-1*gammaaa*t1
     5         2*t14*t27*t9*wght+t16*t41*t5+t13*t14*t43+3.33333333333333
     6         3d-1*t13*t14*t42+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t1*t5*(7.937005259840999d-1*gammaaa*t1
     1         2*t46*t9*wght+7.937005259840999d-1*t12*t25*t9*wght-3.2260
     2         28528009071d-4*t11*t12*t2*t25*wght-3.226028528009071d-4*t
     3         10*t11*t12*t2*t22*wght+7.937005259840999d-1*t12*t22*wght)
     4         -7.937005259840999d-1*t12*t16*t5*t9*wght-2.64566841994699
     5         93d-1*t12*t13*t14*t9*wght+1.0753428426696904d-4*t12*t13*t
     6         2*t44*wght+3.226028528009071d-4*t12*t16*t2*t26*wght+Cmat2
     7         (iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.311232602576965d-7*t1*t10*t12*t47*w
     1         ght-4.8390427920136064d-4*t1*t12*t26*t45*wght+Cmat2(iq,D2
     2         _GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t1*t5*(7.937005259840999d-1*gammaaa*
     1         t12*(4.8774495000000007d-3*t11*t2*t22**3/t9**4+8.12908250
     2         0000001d-4*t11*t2*t24*t49+1.5516762568287046d-2*t10*t2*t2
     3         8/t4**7.166666666666667d+0+5.690357750000001d-3*rhoa*t2*t
     4         23*t38*t39-4.8774495000000007d-3*t11*t2*t22*t36*t39+9.009
     5         733104166671d-3*t15*t2*t22*t24*t37-9.009733104166671d-3*r
     6         hoa*t10*t2*t37-2.8451788750000007d-3*rhoa*t2*t23*t24*t36-
     7         2.8451788750000007d-3*t2*t22*t23*t24)*t9*wght+2.381101577
     8         9522995d+0*gammaaa*t12*t25*t40*t9*wght+7.937005259840999d
     9         -1*gammaaa*t12*t25**3*t9*wght+7.937005259840999d-1*gammaa
     :         a*t12*t49*wght+2.3811015779522995d+0*gammaaa*t12*t22*t40*
     ;         wght+2.3811015779522995d+0*gammaaa*t12*t25*t36*wght+2.381
     <         1015779522995d+0*gammaaa*t12*t22*t30*wght)+t13*t14*(-2.64
     =         56684199469993d-1*gammaaa*t12*t40*t9*wght-2.6456684199469
     >         993d-1*gammaaa*t12*t30*t9*wght-2.6456684199469993d-1*gamm
     ?         aaa*t12*t36*wght-5.291336839893999d-1*gammaaa*t12*t22*t25
     @         *wght)+3.333333333333333d-1*t13*t14*(-1.5874010519681997d
     1         +0*gammaaa*t12*t40*t9*wght-1.5874010519681997d+0*gammaaa*
     2         t12*t30*t9*wght-1.5874010519681997d+0*gammaaa*t12*t36*wgh
     3         t-3.1748021039363994d+0*gammaaa*t12*t22*t25*wght)+t16*t5*
     4         (-2.3811015779522995d+0*gammaaa*t12*t40*t9*wght-2.3811015
     5         779522995d+0*gammaaa*t12*t30*t9*wght-2.3811015779522995d+
     6         0*gammaaa*t12*t36*wght-4.762203155904599d+0*gammaaa*t12*t
     7         22*t25*wght)+t14*t27*(8.818894733156665d-1*gammaaa*t12*t2
     8         5*t9*wght+8.818894733156665d-1*gammaaa*t12*t22*wght)-4.76
     9         2203155904599d+0*gammaaa*t12*t5*t9*wght/rhoa**4-3.4687652
     :         61708288d+0*gammaaa*t12*t14*t9*wght/rhoa**3.3333333333333
     ;         337d+0+3.333333333333333d-1*t14*t27*(t51+t50-7*t43)+t29*t
     <         5*(t51+t50-2*t41)-7.777777777777778d-1*t14*t27*t42+Amat3(
     =         iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t1*t5*(1.5874010519681997d+0*gammaa
     1         a*t12*t25*t46*t9*wght+7.937005259840999d-1*gammaaa*t12*(-
     2         8.129082500000001d-4*t11*t38*t39*t45-1.5016221840277785d-
     3         3*t10*t15*t37*t45+4.0645412500000005d-4*t11*t24*t36*t45-9
     4         .483929583333336d-4*rhoa*t22*t23*t24*t45+4.74196479166666
     5         8d-4*t10*t23*t45)*t9*wght+7.937005259840999d-1*t12*t40*t9
     6         *wght+7.937005259840999d-1*t12*t30*t9*wght+1.587401051968
     7         1997d+0*gammaaa*t12*t22*t46*wght-3.226028528009071d-4*t11
     8         *t12*t2*t40*wght-3.226028528009071d-4*t10*t11*t12*t2*t36*
     9         wght+7.937005259840999d-1*t12*t36*wght-3.226028528009071d
     :         -4*t11*t12*t2*t30*wght-6.452057056018142d-4*t10*t11*t12*t
     ;         2*t22*t25*wght+1.5874010519681997d+0*t12*t22*t25*wght)+t1
     <         3*t14*(-2.6456684199469993d-1*gammaaa*t12*t46*t9*wght-2.6
     =         456684199469993d-1*t12*t25*t9*wght+1.0753428426696904d-4*
     >         t11*t12*t2*t25*wght+1.0753428426696904d-4*t10*t11*t12*t2*
     ?         t22*wght-2.6456684199469993d-1*t12*t22*wght)+3.3333333333
     @         33333d-1*t13*t14*(-7.937005259840999d-1*gammaaa*t12*t46*t
     1         9*wght-7.937005259840999d-1*t12*t25*t9*wght+3.22602852800
     2         9071d-4*t11*t12*t2*t25*wght+3.226028528009071d-4*t10*t11*
     3         t12*t2*t22*wght-7.937005259840999d-1*t12*t22*wght)+t16*t5
     4         *(-1.5874010519681997d+0*gammaaa*t12*t46*t9*wght-1.587401
     5         0519681997d+0*t12*t25*t9*wght+6.452057056018142d-4*t11*t1
     6         2*t2*t25*wght+6.452057056018142d-4*t10*t11*t12*t2*t22*wgh
     7         t-1.5874010519681997d+0*t12*t22*wght)+1.5874010519681997d
     8         +0*t12*t29*t5*t9*wght+8.818894733156665d-1*t12*t14*t27*t9
     9         *wght-3.5844761422323007d-4*t12*t2*t27*t44*wght-6.4520570
     :         56018142d-4*t12*t2*t26*t29*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t1*t5*(7.937005259840999d-1*gammaa
     1         a*t12*(-2.0322706250000003d-4*t11*t22*t24*t53-2.370982395
     2         833334d-4*rhoa*t10*t23*t53)*t9*wght+1.5874010519681997d+0
     3         *t12*t46*t9*wght+1.311232602576965d-7*t10*t12*t25*t52*wgh
     4         t+1.311232602576965d-7*t12*t22*t24*t52*wght-6.45205705601
     5         8142d-4*t11*t12*t2*t46*wght-4.8390427920136064d-4*t11*t12
     6         *t25*t45*wght-4.8390427920136064d-4*t10*t11*t12*t22*t45*w
     7         ght)-1.311232602576965d-7*t10*t12*t16*t47*wght+1.61301426
     8         40045356d-4*t12*t13*t44*t45*wght+4.8390427920136064d-4*t1
     9         2*t16*t26*t45*wght-4.370775341923217d-8*t10*t12*t13*wght/
     :         t4**3.3333333333333337d+0+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 2.4195213960068032d-4*t1*t12*t26*
     1         t53*wght+1.9668489038654477d-7*t1*t10*t12*t47*wght/gammaa
     2         a-5.329559001518931d-11*t1*t12*t24*t45*wght/t4**3.8333333
     3         33333334d+0+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob
            t2 = gammabb**5.0d-1
            t3 = 1/rhob**6.666666666666666d-1
            t4 = abs(rhob)
            t5 = 1/t4**3.333333333333333d-1
            t6 = 1.443307452126544d-2*t5+2.8435438314903855d-6*t3+2.568d
     1         -3
            t7 = 1.7639938117590215d-2*t3*t5+5.411317332115466d+0*t5+1.8
     1         16419932959077d-1*t3+1.0d+0
            t8 = 1/t7
            t9 = t6*t8+1.6670000000000001d-3
            t10 = 1/t9
            t11 = 1/t4**1.166666666666667d+0
            t12 = exp(-8.129082500000001d-4*t10*t11*t2)
            t13 = 1/rhob**1.3333333333333333d+0
            t14 = 1/t4
            t15 = rhob**2
            t16 = 1/t15
            t17 = 1/rhob**1.6666666666666669d+0
            t18 = 1/rhob**3.333333333333333d-1
            t19 = -1.1759958745060142d-2*t17*t5-1.8037724440384886d+0*t1
     1         4*t18-1.2109466219727179d-1*t17-5.879979372530072d-3*t1*t
     2         14
            t20 = 1/t7**2
            t21 = -4.8110248404218126d-3*t14*t18-1.895695887660257d-6*t1
     1         7
            t22 = t21*t8-t19*t20*t6
            t23 = 1/t4**3.166666666666667d+0
            t24 = 1/t9**2
            t25 = 8.129082500000001d-4*t11*t2*t22*t24+9.483929583333336d
     1         -4*rhob*t10*t2*t23
            t26 = 1/t4**1.5000000000000002d+0
            t27 = 1/rhob**2.3333333333333334d+0
            t28 = rhob**3
            t29 = 1/t28
            t30 = t25**2
            t31 = t19**2
            t32 = 1/t7**3
            t33 = 1/rhob**2.6666666666666666d+0
            t34 = 1.9599931241766905d-2*t33*t5+2.01824436995453d-1*t33+1
     1         .5679944993413522d-2*t14*t16+2.4050299253846513d+0*t13*t1
     2         4
            t35 = 3.1594931461004283d-6*t33+6.414699787229083d-3*t13*t14
            t36 = t35*t8-t20*t34*t6+2*t31*t32*t6-2*t19*t20*t21
            t37 = 1/t4**5.166666666666667d+0
            t38 = t22**2
            t39 = 1/t9**3
            t40 = -1.6258165000000002d-3*t11*t2*t38*t39-3.00324436805555
     1         7d-3*t10*t15*t2*t37+8.129082500000001d-4*t11*t2*t24*t36-1
     2         .8967859166666673d-3*rhob*t2*t22*t23*t24+9.48392958333333
     3         6d-4*t10*t2*t23
            t41 = -1.5874010519681997d+0*gammabb*t12*t25*t9*wght-1.58740
     1         10519681997d+0*gammabb*t12*t22*wght
            t42 = -7.937005259840999d-1*gammabb*t12*t25*t9*wght-7.937005
     1         259840999d-1*gammabb*t12*t22*wght
            t43 = -2.6456684199469993d-1*gammabb*t12*t25*t9*wght-2.64566
     1         84199469993d-1*gammabb*t12*t22*wght
            t44 = 1/t4**2.166666666666667d+0
            t45 = 1/t2
            t46 = 4.0645412500000005d-4*t11*t22*t24*t45+4.74196479166666
     1         8d-4*rhob*t10*t23*t45
            t47 = 1/t4**2.666666666666667d+0
            t48 = 1/rhob**3.6666666666666664d+0
            t49 = (-8.425315056267809d-6*t48-1.496763283686786d-2*t14*t2
     1         7)*t8-6*t19**3*t6/t7**4-t20*(-5.226648331137841d-2*t48*t5
     2         -5.381984986545413d-1*t48-5.357314539416287d-2*t14*t29-5.
     3         611736492564187d+0*t14*t27)*t6+6*t19*t32*t34*t6-3*t19*t20
     4         *t35-3*t20*t21*t34+6*t21*t31*t32
            t50 = 1.5874010519681997d+0*gammabb*t12*t22*wght
            t51 = 1.5874010519681997d+0*gammabb*t12*t25*t9*wght
            t52 = 1/t4**2.333333333333334d+0
            t53 = 1/t2**3
            fnc(iq) = 7.937005259840999d-1*gammabb*t1*t12*t5*t9*wght+fnc
     1         (iq)
            Amat(iq,D1_RB) = t1*t5*(7.937005259840999d-1*gammabb*t12*t25
     1         *t9*wght+7.937005259840999d-1*gammabb*t12*t22*wght)-7.937
     2         005259840999d-1*gammabb*t12*t16*t5*t9*wght-2.645668419946
     3         9993d-1*gammabb*t12*t13*t14*t9*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 7.937005259840999d-1*t1*t12*t5*t9*wght-3.2
     1         26028528009071d-4*t1*t12*t2*t26*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t1*t5*(7.937005259840999d-1*gammabb*t12
     1         *t40*t9*wght+7.937005259840999d-1*gammabb*t12*t30*t9*wght
     2         +7.937005259840999d-1*gammabb*t12*t36*wght+1.587401051968
     3         1997d+0*gammabb*t12*t22*t25*wght)+1.5874010519681997d+0*g
     4         ammabb*t12*t29*t5*t9*wght+8.818894733156665d-1*gammabb*t1
     5         2*t14*t27*t9*wght+t16*t41*t5+t13*t14*t43+3.33333333333333
     6         3d-1*t13*t14*t42+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t1*t5*(7.937005259840999d-1*gammabb*t1
     1         2*t46*t9*wght+7.937005259840999d-1*t12*t25*t9*wght-3.2260
     2         28528009071d-4*t11*t12*t2*t25*wght-3.226028528009071d-4*t
     3         10*t11*t12*t2*t22*wght+7.937005259840999d-1*t12*t22*wght)
     4         -7.937005259840999d-1*t12*t16*t5*t9*wght-2.64566841994699
     5         93d-1*t12*t13*t14*t9*wght+1.0753428426696904d-4*t12*t13*t
     6         2*t44*wght+3.226028528009071d-4*t12*t16*t2*t26*wght+Cmat2
     7         (iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.311232602576965d-7*t1*t10*t12*t47*w
     1         ght-4.8390427920136064d-4*t1*t12*t26*t45*wght+Cmat2(iq,D2
     2         _GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t1*t5*(7.937005259840999d-1*gammabb*
     1         t12*(4.8774495000000007d-3*t11*t2*t22**3/t9**4+8.12908250
     2         0000001d-4*t11*t2*t24*t49+1.5516762568287046d-2*t10*t2*t2
     3         8/t4**7.166666666666667d+0+5.690357750000001d-3*rhob*t2*t
     4         23*t38*t39-4.8774495000000007d-3*t11*t2*t22*t36*t39+9.009
     5         733104166671d-3*t15*t2*t22*t24*t37-9.009733104166671d-3*r
     6         hob*t10*t2*t37-2.8451788750000007d-3*rhob*t2*t23*t24*t36-
     7         2.8451788750000007d-3*t2*t22*t23*t24)*t9*wght+2.381101577
     8         9522995d+0*gammabb*t12*t25*t40*t9*wght+7.937005259840999d
     9         -1*gammabb*t12*t25**3*t9*wght+7.937005259840999d-1*gammab
     :         b*t12*t49*wght+2.3811015779522995d+0*gammabb*t12*t22*t40*
     ;         wght+2.3811015779522995d+0*gammabb*t12*t25*t36*wght+2.381
     <         1015779522995d+0*gammabb*t12*t22*t30*wght)+t13*t14*(-2.64
     =         56684199469993d-1*gammabb*t12*t40*t9*wght-2.6456684199469
     >         993d-1*gammabb*t12*t30*t9*wght-2.6456684199469993d-1*gamm
     ?         abb*t12*t36*wght-5.291336839893999d-1*gammabb*t12*t22*t25
     @         *wght)+3.333333333333333d-1*t13*t14*(-1.5874010519681997d
     1         +0*gammabb*t12*t40*t9*wght-1.5874010519681997d+0*gammabb*
     2         t12*t30*t9*wght-1.5874010519681997d+0*gammabb*t12*t36*wgh
     3         t-3.1748021039363994d+0*gammabb*t12*t22*t25*wght)+t16*t5*
     4         (-2.3811015779522995d+0*gammabb*t12*t40*t9*wght-2.3811015
     5         779522995d+0*gammabb*t12*t30*t9*wght-2.3811015779522995d+
     6         0*gammabb*t12*t36*wght-4.762203155904599d+0*gammabb*t12*t
     7         22*t25*wght)+t14*t27*(8.818894733156665d-1*gammabb*t12*t2
     8         5*t9*wght+8.818894733156665d-1*gammabb*t12*t22*wght)-4.76
     9         2203155904599d+0*gammabb*t12*t5*t9*wght/rhob**4-3.4687652
     :         61708288d+0*gammabb*t12*t14*t9*wght/rhob**3.3333333333333
     ;         337d+0+3.333333333333333d-1*t14*t27*(t51+t50-7*t43)+t29*t
     <         5*(t51+t50-2*t41)-7.777777777777778d-1*t14*t27*t42+Amat3(
     =         iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t1*t5*(1.5874010519681997d+0*gammab
     1         b*t12*t25*t46*t9*wght+7.937005259840999d-1*gammabb*t12*(-
     2         8.129082500000001d-4*t11*t38*t39*t45-1.5016221840277785d-
     3         3*t10*t15*t37*t45+4.0645412500000005d-4*t11*t24*t36*t45-9
     4         .483929583333336d-4*rhob*t22*t23*t24*t45+4.74196479166666
     5         8d-4*t10*t23*t45)*t9*wght+7.937005259840999d-1*t12*t40*t9
     6         *wght+7.937005259840999d-1*t12*t30*t9*wght+1.587401051968
     7         1997d+0*gammabb*t12*t22*t46*wght-3.226028528009071d-4*t11
     8         *t12*t2*t40*wght-3.226028528009071d-4*t10*t11*t12*t2*t36*
     9         wght+7.937005259840999d-1*t12*t36*wght-3.226028528009071d
     :         -4*t11*t12*t2*t30*wght-6.452057056018142d-4*t10*t11*t12*t
     ;         2*t22*t25*wght+1.5874010519681997d+0*t12*t22*t25*wght)+t1
     <         3*t14*(-2.6456684199469993d-1*gammabb*t12*t46*t9*wght-2.6
     =         456684199469993d-1*t12*t25*t9*wght+1.0753428426696904d-4*
     >         t11*t12*t2*t25*wght+1.0753428426696904d-4*t10*t11*t12*t2*
     ?         t22*wght-2.6456684199469993d-1*t12*t22*wght)+3.3333333333
     @         33333d-1*t13*t14*(-7.937005259840999d-1*gammabb*t12*t46*t
     1         9*wght-7.937005259840999d-1*t12*t25*t9*wght+3.22602852800
     2         9071d-4*t11*t12*t2*t25*wght+3.226028528009071d-4*t10*t11*
     3         t12*t2*t22*wght-7.937005259840999d-1*t12*t22*wght)+t16*t5
     4         *(-1.5874010519681997d+0*gammabb*t12*t46*t9*wght-1.587401
     5         0519681997d+0*t12*t25*t9*wght+6.452057056018142d-4*t11*t1
     6         2*t2*t25*wght+6.452057056018142d-4*t10*t11*t12*t2*t22*wgh
     7         t-1.5874010519681997d+0*t12*t22*wght)+1.5874010519681997d
     8         +0*t12*t29*t5*t9*wght+8.818894733156665d-1*t12*t14*t27*t9
     9         *wght-3.5844761422323007d-4*t12*t2*t27*t44*wght-6.4520570
     :         56018142d-4*t12*t2*t26*t29*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t1*t5*(7.937005259840999d-1*gammab
     1         b*t12*(-2.0322706250000003d-4*t11*t22*t24*t53-2.370982395
     2         833334d-4*rhob*t10*t23*t53)*t9*wght+1.5874010519681997d+0
     3         *t12*t46*t9*wght+1.311232602576965d-7*t10*t12*t25*t52*wgh
     4         t+1.311232602576965d-7*t12*t22*t24*t52*wght-6.45205705601
     5         8142d-4*t11*t12*t2*t46*wght-4.8390427920136064d-4*t11*t12
     6         *t25*t45*wght-4.8390427920136064d-4*t10*t11*t12*t22*t45*w
     7         ght)-1.311232602576965d-7*t10*t12*t16*t47*wght+1.61301426
     8         40045356d-4*t12*t13*t44*t45*wght+4.8390427920136064d-4*t1
     9         2*t16*t26*t45*wght-4.370775341923217d-8*t10*t12*t13*wght/
     :         t4**3.3333333333333337d+0+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 2.4195213960068032d-4*t1*t12*t26*
     1         t53*wght+1.9668489038654477d-7*t1*t10*t12*t47*wght/gammab
     2         b-5.329559001518931d-11*t1*t12*t24*t45*wght/t4**3.8333333
     3         33333334d+0+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
