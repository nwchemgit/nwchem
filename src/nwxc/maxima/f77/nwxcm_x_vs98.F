C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_vs98.F 
C> The nwxcm_x_vs98 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_vs98 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{\it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({\it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{\it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}\,\sigma_{\beta\beta}^2\,{\it t_{16}}}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0\,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}\,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{\it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{\it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({\it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{\it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,\sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
      subroutine nwxcm_x_vs98(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 2.0d+0*t3*taua-9.11559972d+0
            t5 = 1.86726d-3*t4+1.86726d-3*gammaaa*t2+1.0d+0
            t6 = 1/t5
            t7 = param(4)
            t8 = gammaaa**2
            t9 = 1/rhoa**5.333333333333333d+0
            t10 = 1/t5**2
            t11 = param(5)
            t12 = param(6)
            t13 = t4**2
            t14 = param(2)
            t15 = param(3)
            t16 = t10*t7*t8*t9+t15*t4*t6+gammaaa*t14*t2*t6+gammaaa*t10*t
     1         11*t2*t4+t10*t12*t13+param(1)
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -6.2242d-3*t2*taua-4.9793599999999993d-3*gammaaa*t17
            t19 = 1/t5**3
            t20 = rhoa**3.333333333333333d-1
            t21 = 2.0d+0*t1*t6*(2*gammaaa*t10*t7*t9-3.73452d-3*gammaaa*t
     1         11*t19*t4*t9-1.86726d-3*gammaaa*t10*t14*t9-3.73452d-3*t19
     2         *t7*t8/rhoa**8+t14*t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11
     3         *t2*t4-3.73452d-3*t12*t13*t19*t2)*wght-3.73452d-3*t10*t16
     4         *wght/t1
            t22 = 1/rhoa**4.333333333333333d+0
            fnc(iq) = 2.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = t6*(2.0d+0*t1*(-3.3333333333333337d+0*gamma
     1         aa*t10*t11*t9*taua-3.3333333333333337d+0*t15*t2*t6*taua-6
     2         .666666666666667d+0*t10*t12*t2*t4*taua-2*t18*t19*t7*t8*t9
     3         -5.333333333333333d+0*t10*t7*t8/rhoa**6.333333333333333d+
     4         0-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gammaaa*t11*
     5         t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+0*gammaa
     6         a*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13*t18*t19
     7         )*wght+2.6666666666666666d+0*t16*t20*wght)-2.0d+0*t1*t10*
     8         t16*t18*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t21+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t21+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = 2.0d+0*t1*t6*(-7.46904d-3*t19*t7*t8/rhoa**7
     1         +2.0d+0*t15*t3*t6-3.73452d-3*t10*t15*t3*t4+4.0d+0*t10*t12
     2         *t3*t4-7.46904d-3*gammaaa*t11*t19*t22*t4-7.46904d-3*t12*t
     3         13*t19*t3-3.73452d-3*gammaaa*t10*t14*t22+2.0d+0*gammaaa*t
     4         10*t11*t22)*wght-7.46904d-3*t10*t16*wght/t20+Amat(iq,D1_T
     5         A)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 2.0d+0*t3*taua-9.11559972d+0
            t5 = 1.86726d-3*t4+1.86726d-3*gammaaa*t2+1.0d+0
            t6 = 1/t5
            t7 = param(1)
            t8 = param(4)
            t9 = gammaaa**2
            t10 = 1/rhoa**5.333333333333333d+0
            t11 = 1/t5**2
            t12 = param(5)
            t13 = param(6)
            t14 = t4**2
            t15 = param(2)
            t16 = param(3)
            t17 = t10*t11*t8*t9+t7+t16*t4*t6+gammaaa*t15*t2*t6+gammaaa*t
     1         11*t12*t2*t4+t11*t13*t14
            t18 = rhob**1.3333333333333333d+0
            t19 = 1/rhob**2.6666666666666666d+0
            t20 = 1/rhob**1.6666666666666669d+0
            t21 = 2.0d+0*t20*taub-9.11559972d+0
            t22 = 1.86726d-3*t21+1.86726d-3*gammabb*t19+1.0d+0
            t23 = 1/t22
            t24 = gammabb**2
            t25 = 1/rhob**5.333333333333333d+0
            t26 = 1/t22**2
            t27 = t21**2
            t28 = t24*t25*t26*t8+t7+t13*t26*t27+gammabb*t12*t19*t21*t26+
     1         t16*t21*t23+gammabb*t15*t19*t23
            t29 = 1/rhoa**3.6666666666666664d+0
            t30 = -6.2242d-3*t2*taua-4.9793599999999993d-3*gammaaa*t29
            t31 = 1/t5**3
            t32 = rhoa**3.333333333333333d-1
            t33 = 1/rhob**3.6666666666666664d+0
            t34 = -6.2242d-3*t19*taub-4.9793599999999993d-3*gammabb*t33
            t35 = 1/t22**3
            t36 = rhob**3.333333333333333d-1
            t37 = 1/rhoa**4.333333333333333d+0
            t38 = 1/rhob**4.333333333333333d+0
            fnc(iq) = (1.0d+0*t1*t17*t6+1.0d+0*t18*t23*t28)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t6*(1.0d+0*t1*(-3.3333333333333337d+0*t16*
     1         t2*t6*taua-6.666666666666667d+0*t11*t13*t2*t4*taua-3.3333
     2         333333333337d+0*gammaaa*t10*t11*t12*taua-2*t10*t30*t31*t8
     3         *t9-5.333333333333333d+0*t11*t8*t9/rhoa**6.33333333333333
     4         3d+0-2.6666666666666666d+0*gammaaa*t15*t29*t6-2*gammaaa*t
     5         12*t2*t30*t31*t4-t11*t16*t30*t4-2.6666666666666666d+0*gam
     6         maaa*t11*t12*t29*t4-2*t13*t14*t30*t31-gammaaa*t11*t15*t2*
     7         t30)+1.3333333333333333d+0*t17*t32)-1.0d+0*t1*t11*t17*t30
     8         )*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (t23*(1.0d+0*t18*(-3.3333333333333337d+0*ga
     1         mmabb*t12*t25*t26*taub-6.666666666666667d+0*t13*t19*t21*t
     2         26*taub-3.3333333333333337d+0*t16*t19*t23*taub-2*t24*t25*
     3         t34*t35*t8-5.333333333333333d+0*t24*t26*t8/rhob**6.333333
     4         333333333d+0-2*t13*t27*t34*t35-2*gammabb*t12*t19*t21*t34*
     5         t35-t16*t21*t26*t34-gammabb*t15*t19*t26*t34-2.66666666666
     6         66666d+0*gammabb*t12*t21*t26*t33-2.6666666666666666d+0*ga
     7         mmabb*t15*t23*t33)+1.3333333333333333d+0*t28*t36)-1.0d+0*
     8         t18*t26*t28*t34)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t6*(-3.73452d-3*t31*t8*t9/rhoa*
     1         *8+2*gammaaa*t10*t11*t8+t15*t2*t6-3.73452d-3*gammaaa*t10*
     2         t12*t31*t4-1.86726d-3*t11*t16*t2*t4+t11*t12*t2*t4-3.73452
     3         d-3*t13*t14*t2*t31-1.86726d-3*gammaaa*t10*t11*t15)-1.8672
     4         6d-3*t11*t17/t1)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t18*t23*(-3.73452d-3*t24*t35*t8/rh
     1         ob**8+2*gammabb*t25*t26*t8-3.73452d-3*t13*t19*t27*t35-3.7
     2         3452d-3*gammabb*t12*t21*t25*t35-1.86726d-3*gammabb*t15*t2
     3         5*t26-1.86726d-3*t16*t19*t21*t26+t12*t19*t21*t26+t15*t19*
     4         t23)-1.86726d-3*t26*t28/t18)*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t6*(-7.46904d-3*t31*t8*t9/rhoa**
     1         7+2.0d+0*t16*t3*t6-7.46904d-3*gammaaa*t12*t31*t37*t4-3.73
     2         452d-3*t11*t16*t3*t4+4.0d+0*t11*t13*t3*t4-3.73452d-3*gamm
     3         aaa*t11*t15*t37+2.0d+0*gammaaa*t11*t12*t37-7.46904d-3*t13
     4         *t14*t3*t31)-3.73452d-3*t11*t17/t32)*wght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = (1.0d+0*t18*t23*(-7.46904d-3*t24*t35*t8/rho
     1         b**7-7.46904d-3*gammabb*t12*t21*t35*t38-3.73452d-3*gammab
     2         b*t15*t26*t38+2.0d+0*gammabb*t12*t26*t38-7.46904d-3*t13*t
     3         20*t27*t35-3.73452d-3*t16*t20*t21*t26+4.0d+0*t13*t20*t21*
     4         t26+2.0d+0*t16*t20*t23)-3.73452d-3*t26*t28/t36)*wght+Mmat
     5         (iq,D1_TB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = 1/rhoa**2.6666666666666666d+0
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 2.0d+0*t3*taua-9.11559972d+0
            t5 = 1.86726d-3*t4+1.86726d-3*gammaaa*t2+1.0d+0
            t6 = 1/t5
            t7 = param(4)
            t8 = gammaaa**2
            t9 = 1/rhoa**5.333333333333333d+0
            t10 = 1/t5**2
            t11 = param(5)
            t12 = param(6)
            t13 = t4**2
            t14 = param(2)
            t15 = param(3)
            t16 = t10*t7*t8*t9+t15*t4*t6+gammaaa*t14*t2*t6+gammaaa*t10*t
     1         11*t2*t4+t10*t12*t13+param(1)
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -6.2242d-3*t2*taua-4.9793599999999993d-3*gammaaa*t17
            t19 = 1/t5**3
            t20 = rhoa**3.333333333333333d-1
            t21 = 1/rhoa**4.333333333333333d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = t6*(1.0d+0*t1*(-3.3333333333333337d+0*gamma
     1         aa*t10*t11*t9*taua-3.3333333333333337d+0*t15*t2*t6*taua-6
     2         .666666666666667d+0*t10*t12*t2*t4*taua-2*t18*t19*t7*t8*t9
     3         -5.333333333333333d+0*t10*t7*t8/rhoa**6.333333333333333d+
     4         0-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gammaaa*t11*
     5         t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+0*gammaa
     6         a*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13*t18*t19
     7         )*wght+1.3333333333333333d+0*t16*t20*wght)-1.0d+0*t1*t10*
     8         t16*t18*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t1*t6*(2*gammaaa*t10*t7*t9-3.73452d
     1         -3*gammaaa*t11*t19*t4*t9-1.86726d-3*gammaaa*t10*t14*t9-3.
     2         73452d-3*t19*t7*t8/rhoa**8+t14*t2*t6-1.86726d-3*t10*t15*t
     3         2*t4+t10*t11*t2*t4-3.73452d-3*t12*t13*t19*t2)*wght-1.8672
     4         6d-3*t10*t16*wght/t1+Cmat(iq,D1_GAA)
            Mmat(iq,D1_TA) = 1.0d+0*t1*t6*(-7.46904d-3*t19*t7*t8/rhoa**7
     1         +2.0d+0*t15*t3*t6-3.73452d-3*t10*t15*t3*t4+4.0d+0*t10*t12
     2         *t3*t4-7.46904d-3*gammaaa*t11*t19*t21*t4-7.46904d-3*t12*t
     3         13*t19*t3-3.73452d-3*gammaaa*t10*t14*t21+2.0d+0*gammaaa*t
     4         10*t11*t21)*wght-3.73452d-3*t10*t16*wght/t20+Mmat(iq,D1_T
     5         A)
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = 1/rhob**2.6666666666666666d+0
            t3 = 1/rhob**1.6666666666666669d+0
            t4 = 2.0d+0*t3*taub-9.11559972d+0
            t5 = 1.86726d-3*t4+1.86726d-3*gammabb*t2+1.0d+0
            t6 = 1/t5
            t7 = param(4)
            t8 = gammabb**2
            t9 = 1/rhob**5.333333333333333d+0
            t10 = 1/t5**2
            t11 = param(5)
            t12 = param(6)
            t13 = t4**2
            t14 = param(2)
            t15 = param(3)
            t16 = t10*t7*t8*t9+t15*t4*t6+gammabb*t14*t2*t6+gammabb*t10*t
     1         11*t2*t4+t10*t12*t13+param(1)
            t17 = 1/rhob**3.6666666666666664d+0
            t18 = -6.2242d-3*t2*taub-4.9793599999999993d-3*gammabb*t17
            t19 = 1/t5**3
            t20 = rhob**3.333333333333333d-1
            t21 = 1/rhob**4.333333333333333d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = t6*(1.0d+0*t1*(-3.3333333333333337d+0*gamma
     1         bb*t10*t11*t9*taub-3.3333333333333337d+0*t15*t2*t6*taub-6
     2         .666666666666667d+0*t10*t12*t2*t4*taub-2*t18*t19*t7*t8*t9
     3         -5.333333333333333d+0*t10*t7*t8/rhob**6.333333333333333d+
     4         0-2.6666666666666666d+0*gammabb*t14*t17*t6-2*gammabb*t11*
     5         t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+0*gammab
     6         b*t10*t11*t17*t4-gammabb*t10*t14*t18*t2-2*t12*t13*t18*t19
     7         )*wght+1.3333333333333333d+0*t16*t20*wght)-1.0d+0*t1*t10*
     8         t16*t18*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t1*t6*(2*gammabb*t10*t7*t9-3.73452d
     1         -3*gammabb*t11*t19*t4*t9-1.86726d-3*gammabb*t10*t14*t9-3.
     2         73452d-3*t19*t7*t8/rhob**8+t14*t2*t6-1.86726d-3*t10*t15*t
     3         2*t4+t10*t11*t2*t4-3.73452d-3*t12*t13*t19*t2)*wght-1.8672
     4         6d-3*t10*t16*wght/t1+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TB) = 1.0d+0*t1*t6*(-7.46904d-3*t19*t7*t8/rhob**7
     1         +2.0d+0*t15*t3*t6-3.73452d-3*t10*t15*t3*t4+4.0d+0*t10*t12
     2         *t3*t4-7.46904d-3*gammabb*t11*t19*t21*t4-7.46904d-3*t12*t
     3         13*t19*t3-3.73452d-3*gammabb*t10*t14*t21+2.0d+0*gammabb*t
     4         10*t11*t21)*wght-3.73452d-3*t10*t16*wght/t20+Mmat(iq,D1_T
     5         B)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_vs98 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{\it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({\it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{\it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}\,\sigma_{\beta\beta}^2\,{\it t_{16}}}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0\,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}\,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{\it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{\it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({\it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{\it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,\sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
      subroutine nwxcm_x_vs98_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
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
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_vs98 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{\it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{\it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({\it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{\it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}\,\sigma_{\beta\beta}^2\,{\it t_{16}}}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0\,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}\,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{\it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{\it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{\it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({\it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{\it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,\sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,\sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
      subroutine nwxcm_x_vs98_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,NCOL_RHO)      !< [Input] The density 
      double precision rgamma(nq,NCOL_GAMMA) !< [Input] The norm of the density 
                                             !< gradients 
      double precision tau(nq,NCOL_TAU)      !< [Input] The kinetic energy density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,NCOL_AMAT)   !< [Output] The derivative wrt rho 
      double precision Cmat(nq,NCOL_CMAT)   !< [Output] The derivative wrt rgamma 
      double precision Mmat(nq,NCOL_MMAT)   !< [Output] The derivative wrt tau 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,NCOL_AMAT2)  !< [Output] The 2nd derivative wrt rho 
      double precision Cmat2(nq,NCOL_CMAT2)  !< [Output] The 2nd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat2(nq,NCOL_MMAT2)  !< [Output] The 2nd derivative wrt tau 
                                             !< and possibly rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,NCOL_AMAT3)  !< [Output] The 3rd derivative wrt rho 
      double precision Cmat3(nq,NCOL_CMAT3)  !< [Output] The 3rd derivative wrt rgamma 
                                             !< and possibly rho 
      double precision Mmat3(nq,NCOL_MMAT3)  !< [Output] The 3rd derivative wrt tau 
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
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          taua    = tau(iq,T_A) 
          taub    = tau(iq,T_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
