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
C>   {\it t_3} &=& -\rho_\beta\\\\ 
C>   {\it t_4} &=& {\it t_3}+\rho_\alpha\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& {{0.01443307452126544\,{\it t_7}+2.843543831490385 \times 10^{-6}\,{\it t_5}+0.002568}\over{0.01763993811759022\,{\it t_5}\,{\it t_7}+5.411317332115466\,{\it t_7}+0.1816419932959077\,{\it t_5}+1.0}}+0.001667\\\\ 
C>   {\it t_9} &=& 2.0\,\sigma_{\alpha\beta}\\\\ 
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+{\it t_9}+\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{11}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \left| {\it t_{11}}\right| \\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{16}} &=& {{0.01443307452126544\,{\it t_{15}}+2.843543831490385 \times 10^{-6}\,{\it t_{13}}+0.002568}\over{0.01763993811759022\,{\it t_{13}}\,{\it t_{15}}+5.411317332115466\,{\it t_{15}}+0.1816419932959077\,{\it t_{13}}+1.0}}+0.001667\\\\ 
C>   {\it t_{17}} &=& \sigma_{ss}+\sigma_{\beta\beta}+{\it t_9}\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_8}\,{\it t_7}\,{\it t_{10}}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{\it t_4}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}+0.3149802624737183\,\left(1.0-1.0\,{\it t_4}\,{\it t_2}\right)^{{{5}\over{3}}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{10}}}}\over{{\it t_8}\,{\it t_6}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.414213562373095\,{\it t_{12}}\,{\it t_{16}}\,{\it t_{15}}\,{\it t_{17}}}\over{\left(1.0\,\left(\rho_s+{\it t_3}\right)\,{\it t_{12}}+1.0\right)^{{{5}\over{6}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{17}}}}\over{{\it t_{16}}\,{\it t_{14}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
      subroutine nwxcm_c_perdew86(param,tol_rho,ipol,nq,wght, 
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
            t1 = 3.968502629920499d-1
            t2 = 1/rhoa
            t3 = gammaaa**5.0d-1
            t4 = 6.299605249474366d-1
            t5 = 1/rhoa**6.666666666666666d-1
            t6 = 1.2599210498948732d+0
            t7 = 1/t6
            t8 = abs(rhoa)
            t9 = 1/t8**3.333333333333333d-1
            t10 = 1.443307452126544d-2*t7*t9+2.8435438314903855d-6*t4*t5
     1         +2.568d-3
            t11 = 5.411317332115466d+0*t7*t9+8.819969058795109d-3*t5*t9+
     1         1.816419932959077d-1*t4*t5+1.0d+0
            t12 = 1/t11
            t13 = t10*t12+1.6670000000000001d-3
            t14 = 1/t13
            t15 = 1/t8**1.166666666666667d+0
            t16 = exp(-7.242189178907064d-4*t14*t15*t3)
            t17 = 1/t8
            t18 = 1/rhoa**1.6666666666666669d+0
            t19 = 1/rhoa**3.333333333333333d-1
            t20 = t12*(-4.8110248404218126d-3*t17*t19*t7-9.4784794383012
     1         85d-7*t18*t6)-t10*(-5.879979372530072d-3*t18*t9-1.8037724
     2         440384886d+0*t17*t19*t7-6.054733109863589d-2*t18*t6-2.939
     3         989686265036d-3*t17*t2)/t11**2
            t21 = 4.000000000000001d+0*t1*t13*t16*t2*t9*wght-1.448437835
     1         781413d-3*t1*t16*t2*t3*wght/t8**1.5000000000000002d+0
            fnc(iq) = 4.000000000000001d+0*gammaaa*t1*t13*t16*t2*t9*wght
     1         +fnc(iq)
            Amat(iq,D1_RA) = t1*t2*t9*(4.000000000000001d+0*gammaaa*t13*
     1         t16*(8.449220708724909d-4*rhoa*t14*t3/t8**3.1666666666666
     2         67d+0+7.242189178907064d-4*t15*t20*t3/t13**2)*wght+4.0000
     3         00000000001d+0*gammaaa*t16*t20*wght)-4.000000000000001d+0
     4         *gammaaa*t1*t13*t16*t9*wght/rhoa**2-1.3333333333333336d+0
     5         *gammaaa*t1*t13*t16*t17*wght/rhoa**1.3333333333333333d+0+
     6         Amat(iq,D1_RA)
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
            t1 = rhob+rhoa
            t2 = 1/t1
            t3 = rhoa-rhob
            t4 = 1.0d+0*t2*t3+1.0d+0
            t5 = 1/t4**8.333333333333334d-1
            t6 = 1/t1**6.666666666666666d-1
            t7 = abs(t1)
            t8 = 1/t7**3.333333333333333d-1
            t9 = 1.443307452126544d-2*t8+2.8435438314903855d-6*t6+2.568d
     1         -3
            t10 = 1.7639938117590215d-2*t6*t8+5.411317332115466d+0*t8+1.
     1         816419932959077d-1*t6+1.0d+0
            t11 = 1/t10
            t12 = t11*t9+1.6670000000000001d-3
            t13 = sigmabb+2.0d+0*sigmaab+gammaaa
            t14 = 1/t12
            t15 = 1/t7**1.166666666666667d+0
            t16 = t13**5.0d-1
            t17 = exp(-8.129082500000001d-4*t14*t15*t16)
            t18 = 1/t7
            t19 = 1/t1**2
            t20 = 1/t1**1.6666666666666669d+0
            t21 = 1/t1**3.333333333333333d-1
            t22 = t11*(-4.8110248404218126d-3*t18*t21-1.895695887660257d
     1         -6*t20)-(-1.1759958745060142d-2*t20*t8-1.8037724440384886
     2         d+0*t18*t21-1.2109466219727179d-1*t20-5.879979372530072d-
     3         3*t18*t2)*t9/t10**2
            fnc(iq) = 1.414213562373095d+0*t12*t13*t17*t2*t5*t8*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = t2*t5*t8*(1.414213562373095d+0*t12*t13*t17*
     1         (9.483929583333336d-4*t1*t14*t16/t7**3.166666666666667d+0
     2         +8.129082500000001d-4*t15*t16*t22/t12**2)*wght+1.41421356
     3         2373095d+0*t13*t17*t22*wght)-1.414213562373095d+0*t12*t13
     4         *t17*t19*t5*t8*wght-1.1785113019775792d+0*t12*t13*t17*t2*
     5         (1.0d+0*t2-1.0d+0*t19*t3)*t8*wght/t4**1.8333333333333332d
     6         +0-4.714045207910317d-1*t12*t13*t17*t18*t5*wght/t1**1.333
     7         3333333333333d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.414213562373095d+0*t12*t17*t2*t5*t8*wght
     1         -5.748129360574894d-4*t16*t17*t2*t5*wght/t7**1.5000000000
     2         000002d+0+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 3.968502629920499d-1
            t2 = 1/rhob
            t3 = 6.299605249474366d-1
            t4 = 1/rhob**6.666666666666666d-1
            t5 = 1.2599210498948732d+0
            t6 = 1/t5
            t7 = abs(rhob)
            t8 = 1/t7**3.333333333333333d-1
            t9 = 1.443307452126544d-2*t6*t8+2.8435438314903855d-6*t3*t4+
     1         2.568d-3
            t10 = 5.411317332115466d+0*t6*t8+8.819969058795109d-3*t4*t8+
     1         1.816419932959077d-1*t3*t4+1.0d+0
            t11 = 1/t10
            t12 = t11*t9+1.6670000000000001d-3
            t13 = sigmabb+2.0d+0*sigmaab+gammabb
            t14 = 1/t12
            t15 = 1/t7**1.166666666666667d+0
            t16 = t13**5.0d-1
            t17 = exp(-3.621094589453532d-4*t14*t15*t16)
            t18 = 1/t7
            t19 = 1/rhob**1.6666666666666669d+0
            t20 = 1/rhob**3.333333333333333d-1
            t21 = t11*(-4.8110248404218126d-3*t18*t20*t6-9.4784794383012
     1         85d-7*t19*t5)-(-5.879979372530072d-3*t19*t8-1.80377244403
     2         84886d+0*t18*t20*t6-6.054733109863589d-2*t19*t5-2.9399896
     3         86265036d-3*t18*t2)*t9/t10**2
            fnc(iq) = 1.414213562373095d+0*t1*t12*t13*t17*t2*t8*wght+fnc
     1         (iq)
            Amat(iq,D1_RB) = t1*t2*t8*(1.414213562373095d+0*t12*t13*t17*
     1         (4.2246103543624547d-4*rhob*t14*t16/t7**3.166666666666667
     2         d+0+3.621094589453532d-4*t15*t16*t21/t12**2)*wght+1.41421
     3         3562373095d+0*t13*t17*t21*wght)-1.414213562373095d+0*t1*t
     4         12*t13*t17*t8*wght/rhob**2-4.714045207910317d-1*t1*t12*t1
     5         3*t17*t18*wght/rhob**1.3333333333333333d+0+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.414213562373095d+0*t1*t12*t17*t2*t8*wght
     1         -2.5605005395205094d-4*t1*t16*t17*t2*wght/t7**1.500000000
     2         0000002d+0+Cmat(iq,D1_GBB)
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
C>   {\it t_3} &=& -\rho_\beta\\\\ 
C>   {\it t_4} &=& {\it t_3}+\rho_\alpha\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& {{0.01443307452126544\,{\it t_7}+2.843543831490385 \times 10^{-6}\,{\it t_5}+0.002568}\over{0.01763993811759022\,{\it t_5}\,{\it t_7}+5.411317332115466\,{\it t_7}+0.1816419932959077\,{\it t_5}+1.0}}+0.001667\\\\ 
C>   {\it t_9} &=& 2.0\,\sigma_{\alpha\beta}\\\\ 
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+{\it t_9}+\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{11}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \left| {\it t_{11}}\right| \\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{16}} &=& {{0.01443307452126544\,{\it t_{15}}+2.843543831490385 \times 10^{-6}\,{\it t_{13}}+0.002568}\over{0.01763993811759022\,{\it t_{13}}\,{\it t_{15}}+5.411317332115466\,{\it t_{15}}+0.1816419932959077\,{\it t_{13}}+1.0}}+0.001667\\\\ 
C>   {\it t_{17}} &=& \sigma_{ss}+\sigma_{\beta\beta}+{\it t_9}\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_8}\,{\it t_7}\,{\it t_{10}}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{\it t_4}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}+0.3149802624737183\,\left(1.0-1.0\,{\it t_4}\,{\it t_2}\right)^{{{5}\over{3}}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{10}}}}\over{{\it t_8}\,{\it t_6}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.414213562373095\,{\it t_{12}}\,{\it t_{16}}\,{\it t_{15}}\,{\it t_{17}}}\over{\left(1.0\,\left(\rho_s+{\it t_3}\right)\,{\it t_{12}}+1.0\right)^{{{5}\over{6}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{17}}}}\over{{\it t_{16}}\,{\it t_{14}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
      subroutine nwxcm_c_perdew86_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_c_perdew86 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\beta+\rho_\alpha\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& -\rho_\beta\\\\ 
C>   {\it t_4} &=& {\it t_3}+\rho_\alpha\\\\ 
C>   {\it t_5} &=& {{1}\over{{\it t_1}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \left| {\it t_1}\right| \\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_6}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& {{0.01443307452126544\,{\it t_7}+2.843543831490385 \times 10^{-6}\,{\it t_5}+0.002568}\over{0.01763993811759022\,{\it t_5}\,{\it t_7}+5.411317332115466\,{\it t_7}+0.1816419932959077\,{\it t_5}+1.0}}+0.001667\\\\ 
C>   {\it t_9} &=& 2.0\,\sigma_{\alpha\beta}\\\\ 
C>   {\it t_{10}} &=& \sigma_{\beta\beta}+{\it t_9}+\sigma_{\alpha\alpha}\\\\ 
C>   {\it t_{11}} &=& \rho_s+\rho_\beta\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_{11}}}}\\\\ 
C>   {\it t_{13}} &=& {{1}\over{{\it t_{11}}^{{{2}\over{3}}}}}\\\\ 
C>   {\it t_{14}} &=& \left| {\it t_{11}}\right| \\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}^{{{1}\over{3}}}}}\\\\ 
C>   {\it t_{16}} &=& {{0.01443307452126544\,{\it t_{15}}+2.843543831490385 \times 10^{-6}\,{\it t_{13}}+0.002568}\over{0.01763993811759022\,{\it t_{13}}\,{\it t_{15}}+5.411317332115466\,{\it t_{15}}+0.1816419932959077\,{\it t_{13}}+1.0}}+0.001667\\\\ 
C>   {\it t_{17}} &=& \sigma_{ss}+\sigma_{\beta\beta}+{\it t_9}\\\\ 
C>   f &=& {{0.7937005259840999\,{\it t_2}\,{\it t_8}\,{\it t_7}\,{\it t_{10}}}\over{\sqrt{0.3149802624737183\,\left(1.0\,{\it t_4}\,{\it t_2}+1.0\right)^{{{5}\over{3}}}+0.3149802624737183\,\left(1.0-1.0\,{\it t_4}\,{\it t_2}\right)^{{{5}\over{3}}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{10}}}}\over{{\it t_8}\,{\it t_6}^{1.166666666666667}}}}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.414213562373095\,{\it t_{12}}\,{\it t_{16}}\,{\it t_{15}}\,{\it t_{17}}}\over{\left(1.0\,\left(\rho_s+{\it t_3}\right)\,{\it t_{12}}+1.0\right)^{{{5}\over{6}}}\,e^{{{8.129082500000001 \times 10^{-4}\,\sqrt{{\it t_{17}}}}\over{{\it t_{16}}\,{\it t_{14}}^{1.166666666666667}}}}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, Phys.Rev. B 33, 8822 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.33.8822 "> 
C> 10.1103/PhysRevB.33.8822 </a> 
C> 
      subroutine nwxcm_c_perdew86_d3(param,tol_rho,ipol,nq,wght, 
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
