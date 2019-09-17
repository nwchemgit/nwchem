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
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}
C>    -9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{
C>    \it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({
C>    \it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{
C>    \it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,
C>    \sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}
C>    \,\sigma_{\beta\beta}^2\,{\it t_{16}}}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0
C>    \,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}
C>    \,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{
C>    \it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{
C>    \it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{
C>    \it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{
C>    \it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({
C>    \it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{
C>    \it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{
C>    \it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,
C>    \sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,
C>    \sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}
C>    +{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_vs98(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
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
            t21 = 1/rhoa**4.333333333333333d+0
            fnc(iq) = 2.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = (t6*(1.0d+0*t1*(-3.3333333333333337d+0*gamm
     1         aaa*t10*t11*t9*taua-3.3333333333333337d+0*t15*t2*t6*taua-
     2         6.666666666666667d+0*t10*t12*t2*t4*taua-2*t18*t19*t7*t8*t
     3         9-5.333333333333333d+0*t10*t7*t8/rhoa**6.333333333333333d
     4         +0-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gammaaa*t11
     5         *t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+0*gamma
     6         aa*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13*t18*t1
     7         9)+1.3333333333333333d+0*t16*t20)-1.0d+0*t1*t10*t16*t18)*
     8         wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t6*(2*gammaaa*t10*t7*t9-3.73452
     1         d-3*gammaaa*t11*t19*t4*t9-1.86726d-3*gammaaa*t10*t14*t9-3
     2         .73452d-3*t19*t7*t8/rhoa**8+t14*t2*t6-1.86726d-3*t10*t15*
     3         t2*t4+t10*t11*t2*t4-3.73452d-3*t12*t13*t19*t2)-1.86726d-3
     4         *t10*t16/t1)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t6*(-7.46904d-3*t19*t7*t8/rhoa**
     1         7+2.0d+0*t15*t3*t6-3.73452d-3*t10*t15*t3*t4+4.0d+0*t10*t1
     2         2*t3*t4-7.46904d-3*gammaaa*t11*t19*t21*t4-7.46904d-3*t12*
     3         t13*t19*t3-3.73452d-3*gammaaa*t10*t14*t21+2.0d+0*gammaaa*
     4         t10*t11*t21)-3.73452d-3*t10*t16/t20)*wght+Mmat(iq,D1_TA)
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
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}
C>    -9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{
C>    \it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({
C>    \it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{
C>    \it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,
C>    \sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}
C>    \,\sigma_{\beta\beta}^2\,{\it t_{16}}}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0
C>    \,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}
C>    \,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{
C>    \it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{
C>    \it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{
C>    \it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{
C>    \it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({
C>    \it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{
C>    \it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{
C>    \it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,
C>    \sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,
C>    \sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}
C>    +{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_vs98_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
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
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammaaa*t10*t11*t9*taua-3.33333
     1         33333333337d+0*t15*t2*t6*taua-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taua-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gamm
     4         aaa*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammaaa*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhoa**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22+1.0d+0*t1*t21
            t24 = 1/rhoa**8
            t25 = 2*gammaaa*t10*t7*t9-3.73452d-3*gammaaa*t11*t19*t4*t9-1
     1         .86726d-3*gammaaa*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhoa**7
            t28 = 1/rhoa**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammaaa*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammaaa*
     3         t10*t14*t28+2.0d+0*gammaaa*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhoa**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taua+1.8257653333333332d-2*g
     1         ammaaa*t33
            t35 = 1/rhoa**6
            t36 = 1/rhoa**3.3333333333333337d+0
            fnc(iq) = 2.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = (t23*t6-1.0d+0*t1*t10*t16*t18)*wght+Amat(iq
     1         ,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t25*t6-1.86726d-3*t10*t16*t26)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t29*t6-3.73452d-3*t10*t16*t30)*w
     1         ght+Mmat(iq,D1_TA)
            Amat2(iq,D2_RA_RA) = (t6*(1.0d+0*t1*(2.2222222222222224d+1*t
     1         10*t12*t9*taua**2+1.3333333333333333d+1*gammaaa*t11*t18*t
     2         19*t9*taua+8.88888888888889d+0*t15*t17*t6*taua+2.66666666
     3         66666666d+1*t12*t18*t19*t2*t4*taua+1.7777777777777778d+1*
     4         t10*t12*t17*t4*taua+2.6666666666666666d+1*gammaaa*t10*t11
     5         *t20*taua+6.666666666666667d+0*t10*t15*t18*t2*taua-2*t19*
     6         t34*t7*t8*t9+6*t31*t32*t7*t8*t9+2.1333333333333332d+1*t18
     7         *t19*t20*t7*t8+3.377777777777778d+1*t10*t7*t8/rhoa**7.333
     8         333333333333d+0+9.777777777777779d+0*gammaaa*t14*t33*t6-2
     9         *gammaaa*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.7777777777777
     :         79d+0*gammaaa*t10*t11*t33*t4+6*gammaaa*t11*t2*t31*t32*t4+
     ;         2*t15*t19*t31*t4+1.0666666666666666d+1*gammaaa*t11*t17*t1
     <         8*t19*t4-gammaaa*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t
     =         13*t31*t32+2*gammaaa*t14*t19*t2*t31+5.333333333333333d+0*
     >         gammaaa*t10*t14*t17*t18)+2.6666666666666666d+0*t21*t22+4.
     ?         444444444444444d-1*t16/rhoa**6.666666666666666d-1)+t10*(-
     @         1.0d+0*t1*t16*t34-t18*t23-1.3333333333333333d+0*t16*t18*t
     1         22-1.0d+0*t1*t18*t21)+2.0d+0*t1*t16*t19*t31)*wght+Amat2(i
     2         q,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t6*(1.0d+0*t1*(2.48968d-2*t12*t19*t4*
     1         t9*taua+6.2242d-3*t10*t15*t9*taua-3.3333333333333337d+0*t
     2         10*t11*t9*taua+1.24484d-2*gammaaa*t11*t19*t24*taua-4*gamm
     3         aaa*t18*t19*t7*t9+1.120356d-2*gammaaa*t11*t18*t32*t4*t9+3
     4         .73452d-3*gammaaa*t14*t18*t19*t9+1.120356d-2*t18*t24*t32*
     5         t7*t8+2.9876159999999996d-2*t19*t7*t8/rhoa**9-1.066666666
     6         6666666d+1*gammaaa*t10*t20*t7-2.6666666666666666d+0*t14*t
     7         17*t6+1.9917439999999997d-2*gammaaa*t11*t19*t20*t4+3.7345
     8         2d-3*t15*t18*t19*t2*t4-2*t11*t18*t19*t2*t4+4.979359999999
     9         9993d-3*t10*t15*t17*t4-2.6666666666666666d+0*t10*t11*t17*
     :         t4+1.120356d-2*t12*t13*t18*t2*t32+9.958719999999999d-3*ga
     ;         mmaaa*t10*t14*t20-t10*t14*t18*t2+9.958719999999999d-3*t12
     <         *t13*t17*t19)+1.3333333333333333d+0*t22*t25)+3.73452d-3*t
     =         16*t18*t19*t26-1.0d+0*t1*t10*t18*t25-1.86726d-3*t10*t2*t2
     >         3+4.9793599999999993d-3*t10*t16/rhoa**2.3333333333333334d
     ?         +0)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t1*t6*(2*t10*t7*t9+6.97331981
     1         5199999d-6*t15*t19*t4*t9-7.46904d-3*t11*t19*t4*t9+2.09199
     2         59445599998d-5*t12*t13*t32*t9-3.73452d-3*t10*t14*t9+2.091
     3         9959445599998d-5*t32*t7*t8/rhoa**1.0666666666666666d+1-1.
     4         4938079999999998d-2*gammaaa*t19*t24*t7+2.0919959445599998
     5         d-5*gammaaa*t11*t24*t32*t4+6.973319815199999d-6*gammaaa*t
     6         14*t19*t24)-3.73452d-3*t10*t25*t26+6.973319815199999d-6*t
     7         16*t19/rhoa**4)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (t6*(1.0d+0*t1*(4.97936d-2*t12*t19*t28*
     1         t4*taua+1.24484d-2*t10*t15*t28*taua-1.3333333333333333d+1
     2         *t10*t12*t28*taua+2.48968d-2*gammaaa*t11*t19*t27*taua+3.2
     3         36584d-2*gammaaa*t11*t19*t4*t9+1.618292d-2*gammaaa*t10*t1
     4         4*t9-8.666666666666666d+0*gammaaa*t10*t11*t9+2.240712d-2*
     5         t18*t27*t32*t7*t8+5.228328d-2*t19*t24*t7*t8-3.33333333333
     6         33337d+0*t15*t2*t6+2.240712d-2*gammaaa*t11*t18*t28*t32*t4
     7         +7.46904d-3*t15*t18*t19*t3*t4-8.0d+0*t12*t18*t19*t3*t4+6.
     8         2242d-3*t10*t15*t2*t4-6.666666666666667d+0*t10*t12*t2*t4+
     9         2.240712d-2*t12*t13*t18*t3*t32-2.0d+0*t10*t15*t18*t3+7.46
     :         904d-3*gammaaa*t14*t18*t19*t28-4.0d+0*gammaaa*t11*t18*t19
     ;         *t28+1.24484d-2*t12*t13*t19*t2)+1.3333333333333333d+0*t22
     <         *t29)+7.46904d-3*t16*t18*t19*t30-3.73452d-3*t10*t23*t3-1.
     =         0d+0*t1*t10*t18*t29+6.2242d-3*t10*t16*t26)*wght+Mmat2(iq,
     >         D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = (1.0d+0*t1*t6*(4.1839918891199995d-5*t
     1         32*t7*t8/rhoa**9.666666666666666d+0-1.4938079999999998d-2
     2         *gammaaa*t19*t27*t7+4.1839918891199995d-5*gammaaa*t11*t27
     3         *t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.4938079999
     4         999998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.1839
     5         918891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28-3.
     6         73452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.394663963039999
     7         8d-5*gammaaa*t14*t19*t27-7.46904d-3*gammaaa*t11*t19*t27)-
     8         3.73452d-3*t10*t25*t30-1.86726d-3*t10*t26*t29+1.394663963
     9         0399998d-5*t16*t19/rhoa**3)*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = (1.0d+0*t1*t6*(8.367983778239999d-5*t32
     1         *t7*t8/rhoa**8.666666666666666d+0+2.7893279260799997d-5*t
     2         15*t19*t36*t4-5.975231999999999d-2*t12*t19*t36*t4+8.36798
     3         3778239999d-5*gammaaa*t11*t32*t35*t4+8.367983778239999d-5
     4         *t12*t13*t32*t36-1.4938079999999998d-2*t10*t15*t36+8.0d+0
     5         *t10*t12*t36+2.7893279260799997d-5*gammaaa*t14*t19*t35-2.
     6         9876159999999996d-2*gammaaa*t11*t19*t35)-7.46904d-3*t10*t
     7         29*t30+2.7893279260799997d-5*t16*t19/rhoa**2)*wght+Mmat2(
     8         iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
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
            t32 = 1/rhoa**6.333333333333333d+0
            t33 = -3.3333333333333337d+0*t16*t2*t6*taua-6.66666666666666
     1         7d+0*t11*t13*t2*t4*taua-3.3333333333333337d+0*gammaaa*t10
     2         *t11*t12*taua-5.333333333333333d+0*t11*t32*t8*t9-2*t10*t3
     3         0*t31*t8*t9-2.6666666666666666d+0*gammaaa*t15*t29*t6-2*ga
     4         mmaaa*t12*t2*t30*t31*t4-t11*t16*t30*t4-2.6666666666666666
     5         d+0*gammaaa*t11*t12*t29*t4-2*t13*t14*t30*t31-gammaaa*t11*
     6         t15*t2*t30
            t34 = rhoa**3.333333333333333d-1
            t35 = 1.3333333333333333d+0*t17*t34+1.0d+0*t1*t33
            t36 = 1/rhob**3.6666666666666664d+0
            t37 = -6.2242d-3*t19*taub-4.9793599999999993d-3*gammabb*t36
            t38 = 1/t22**3
            t39 = 1/rhob**6.333333333333333d+0
            t40 = -3.3333333333333337d+0*gammabb*t12*t25*t26*taub-6.6666
     1         66666666667d+0*t13*t19*t21*t26*taub-3.3333333333333337d+0
     2         *t16*t19*t23*taub-5.333333333333333d+0*t24*t26*t39*t8-2*t
     3         24*t25*t37*t38*t8-2*t13*t27*t37*t38-2*gammabb*t12*t19*t21
     4         *t37*t38-t16*t21*t26*t37-gammabb*t15*t19*t26*t37-2.666666
     5         6666666666d+0*gammabb*t12*t21*t26*t36-2.6666666666666666d
     6         +0*gammabb*t15*t23*t36
            t41 = rhob**3.333333333333333d-1
            t42 = 1.3333333333333333d+0*t28*t41+1.0d+0*t18*t40
            t43 = 1/rhoa**8
            t44 = -3.73452d-3*t31*t43*t8*t9+2*gammaaa*t10*t11*t8+t15*t2*
     1         t6-3.73452d-3*gammaaa*t10*t12*t31*t4-1.86726d-3*t11*t16*t
     2         2*t4+t11*t12*t2*t4-3.73452d-3*t13*t14*t2*t31-1.86726d-3*g
     3         ammaaa*t10*t11*t15
            t45 = 1/t1
            t46 = 1/rhob**8
            t47 = -3.73452d-3*t24*t38*t46*t8+2*gammabb*t25*t26*t8-3.7345
     1         2d-3*t13*t19*t27*t38-3.73452d-3*gammabb*t12*t21*t25*t38-1
     2         .86726d-3*gammabb*t15*t25*t26-1.86726d-3*t16*t19*t21*t26+
     3         t12*t19*t21*t26+t15*t19*t23
            t48 = 1/t18
            t49 = 1/rhoa**7
            t50 = 1/rhoa**4.333333333333333d+0
            t51 = -7.46904d-3*t31*t49*t8*t9+2.0d+0*t16*t3*t6-7.46904d-3*
     1         gammaaa*t12*t31*t4*t50-3.73452d-3*gammaaa*t11*t15*t50+2.0
     2         d+0*gammaaa*t11*t12*t50-3.73452d-3*t11*t16*t3*t4+4.0d+0*t
     3         11*t13*t3*t4-7.46904d-3*t13*t14*t3*t31
            t52 = 1/t34
            t53 = 1/rhob**7
            t54 = 1/rhob**4.333333333333333d+0
            t55 = -7.46904d-3*t24*t38*t53*t8-7.46904d-3*gammabb*t12*t21*
     1         t38*t54-3.73452d-3*gammabb*t15*t26*t54+2.0d+0*gammabb*t12
     2         *t26*t54-7.46904d-3*t13*t20*t27*t38-3.73452d-3*t16*t20*t2
     3         1*t26+4.0d+0*t13*t20*t21*t26+2.0d+0*t16*t20*t23
            t56 = 1/t41
            t57 = t30**2
            t58 = 1/t5**4
            t59 = 1/rhoa**4.666666666666667d+0
            t60 = 1.6597866666666666d-2*t29*taua+1.8257653333333332d-2*g
     1         ammaaa*t59
            t61 = t37**2
            t62 = 1/t22**4
            t63 = 1/rhob**4.666666666666667d+0
            t64 = 1.6597866666666666d-2*t36*taub+1.8257653333333332d-2*g
     1         ammabb*t63
            t65 = 1/rhoa**6
            t66 = 1/rhoa**3.3333333333333337d+0
            t67 = 1/rhob**6
            t68 = 1/rhob**3.3333333333333337d+0
            fnc(iq) = (1.0d+0*t1*t17*t6+1.0d+0*t18*t23*t28)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t35*t6-1.0d+0*t1*t11*t17*t30)*wght+Amat(iq
     1         ,D1_RA)
            Amat(iq,D1_RB) = (t23*t42-1.0d+0*t18*t26*t28*t37)*wght+Amat(
     1         iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t44*t6-1.86726d-3*t11*t17*t45)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t18*t23*t47-1.86726d-3*t26*t28*t48
     1         )*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t51*t6-3.73452d-3*t11*t17*t52)*w
     1         ght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = (1.0d+0*t18*t23*t55-3.73452d-3*t26*t28*t56)
     1         *wght+Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = (t6*(1.0d+0*t1*(2.2222222222222224d+1*t
     1         10*t11*t13*taua**2+8.88888888888889d+0*t16*t29*t6*taua+2.
     2         6666666666666666d+1*t13*t2*t30*t31*t4*taua+1.777777777777
     3         7778d+1*t11*t13*t29*t4*taua+2.6666666666666666d+1*gammaaa
     4         *t11*t12*t32*taua+1.3333333333333333d+1*gammaaa*t10*t12*t
     5         30*t31*taua+6.666666666666667d+0*t11*t16*t2*t30*taua-2*t1
     6         0*t31*t60*t8*t9+6*t10*t57*t58*t8*t9+2.1333333333333332d+1
     7         *t30*t31*t32*t8*t9+3.377777777777778d+1*t11*t8*t9/rhoa**7
     8         .333333333333333d+0-2*gammaaa*t12*t2*t31*t4*t60-t11*t16*t
     9         4*t60-2*t13*t14*t31*t60-gammaaa*t11*t15*t2*t60+9.77777777
     :         7777779d+0*gammaaa*t15*t59*t6+9.777777777777779d+0*gammaa
     ;         a*t11*t12*t4*t59+6*gammaaa*t12*t2*t4*t57*t58+6*t13*t14*t5
     <         7*t58+2*t16*t31*t4*t57+2*gammaaa*t15*t2*t31*t57+1.0666666
     =         666666666d+1*gammaaa*t12*t29*t30*t31*t4+5.333333333333333
     >         d+0*gammaaa*t11*t15*t29*t30)+2.6666666666666666d+0*t33*t3
     ?         4+4.444444444444444d-1*t17/rhoa**6.666666666666666d-1)+t1
     @         1*(-1.0d+0*t1*t17*t60-t30*t35-1.3333333333333333d+0*t17*t
     1         30*t34-1.0d+0*t1*t30*t33)+2.0d+0*t1*t17*t31*t57)*wght+Ama
     2         t2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t23*(1.0d+0*t18*(2.2222222222222224d+1
     1         *t13*t25*t26*taub**2+2.6666666666666666d+1*gammabb*t12*t2
     2         6*t39*taub+1.3333333333333333d+1*gammabb*t12*t25*t37*t38*
     3         taub+2.6666666666666666d+1*t13*t19*t21*t37*t38*taub+6.666
     4         666666666667d+0*t16*t19*t26*t37*taub+1.7777777777777778d+
     5         1*t13*t21*t26*t36*taub+8.88888888888889d+0*t16*t23*t36*ta
     6         ub-2*t24*t25*t38*t64*t8+6*t24*t25*t61*t62*t8+2.1333333333
     7         333332d+1*t24*t37*t38*t39*t8+3.377777777777778d+1*t24*t26
     8         *t8/rhob**7.333333333333333d+0-2*t13*t27*t38*t64-2*gammab
     9         b*t12*t19*t21*t38*t64-t16*t21*t26*t64-gammabb*t15*t19*t26
     :         *t64+9.777777777777779d+0*gammabb*t12*t21*t26*t63+9.77777
     ;         7777777779d+0*gammabb*t15*t23*t63+6*t13*t27*t61*t62+6*gam
     <         mabb*t12*t19*t21*t61*t62+2*t16*t21*t38*t61+2*gammabb*t15*
     =         t19*t38*t61+1.0666666666666666d+1*gammabb*t12*t21*t36*t37
     >         *t38+5.333333333333333d+0*gammabb*t15*t26*t36*t37)+2.6666
     ?         666666666666d+0*t40*t41+4.444444444444444d-1*t28/rhob**6.
     @         666666666666666d-1)+t26*(-1.0d+0*t18*t28*t64-t37*t42-1.33
     1         33333333333333d+0*t28*t37*t41-1.0d+0*t18*t37*t40)+2.0d+0*
     2         t18*t28*t38*t61)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t6*(1.0d+0*t1*(1.24484d-2*gammaaa*t12
     1         *t31*t43*taua+2.48968d-2*t10*t13*t31*t4*taua+6.2242d-3*t1
     2         0*t11*t16*taua-3.3333333333333337d+0*t10*t11*t12*taua+1.1
     3         20356d-2*t30*t43*t58*t8*t9+2.9876159999999996d-2*t31*t8*t
     4         9/rhoa**9-1.0666666666666666d+1*gammaaa*t11*t32*t8-4*gamm
     5         aaa*t10*t30*t31*t8-2.6666666666666666d+0*t15*t29*t6+1.120
     6         356d-2*gammaaa*t10*t12*t30*t4*t58+1.120356d-2*t13*t14*t2*
     7         t30*t58+1.9917439999999997d-2*gammaaa*t12*t31*t32*t4+3.73
     8         452d-3*t16*t2*t30*t31*t4-2*t12*t2*t30*t31*t4+4.9793599999
     9         999993d-3*t11*t16*t29*t4-2.6666666666666666d+0*t11*t12*t2
     :         9*t4+9.958719999999999d-3*gammaaa*t11*t15*t32+3.73452d-3*
     ;         gammaaa*t10*t15*t30*t31+9.958719999999999d-3*t13*t14*t29*
     <         t31-t11*t15*t2*t30)+1.3333333333333333d+0*t34*t44)+3.7345
     =         2d-3*t17*t30*t31*t45-1.0d+0*t1*t11*t30*t44-1.86726d-3*t11
     >         *t2*t35+4.9793599999999993d-3*t11*t17/rhoa**2.33333333333
     ?         33334d+0)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t23*(1.0d+0*t18*(1.24484d-2*gammabb*t
     1         12*t38*t46*taub+2.48968d-2*t13*t21*t25*t38*taub+6.2242d-3
     2         *t16*t25*t26*taub-3.3333333333333337d+0*t12*t25*t26*taub+
     3         1.120356d-2*t24*t37*t46*t62*t8-1.0666666666666666d+1*gamm
     4         abb*t26*t39*t8-4*gammabb*t25*t37*t38*t8+2.987615999999999
     5         6d-2*t24*t38*t8/rhob**9+1.120356d-2*t13*t19*t27*t37*t62+1
     6         .120356d-2*gammabb*t12*t21*t25*t37*t62+1.9917439999999997
     7         d-2*gammabb*t12*t21*t38*t39+9.958719999999999d-3*gammabb*
     8         t15*t26*t39+3.73452d-3*gammabb*t15*t25*t37*t38+3.73452d-3
     9         *t16*t19*t21*t37*t38-2*t12*t19*t21*t37*t38+9.958719999999
     :         999d-3*t13*t27*t36*t38-t15*t19*t26*t37+4.9793599999999993
     ;         d-3*t16*t21*t26*t36-2.6666666666666666d+0*t12*t21*t26*t36
     <         -2.6666666666666666d+0*t15*t23*t36)+1.3333333333333333d+0
     =         *t41*t47)+3.73452d-3*t28*t37*t38*t48-1.0d+0*t18*t26*t37*t
     >         47-1.86726d-3*t19*t26*t42+4.9793599999999993d-3*t26*t28/r
     ?         hob**2.3333333333333334d+0)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t1*t6*(2.0919959445599998d-5*
     1         t58*t8*t9/rhoa**1.0666666666666666d+1-1.4938079999999998d
     2         -2*gammaaa*t31*t43*t8+2*t10*t11*t8+2.0919959445599998d-5*
     3         gammaaa*t12*t4*t43*t58+2.0919959445599998d-5*t10*t13*t14*
     4         t58+6.973319815199999d-6*gammaaa*t15*t31*t43+6.9733198151
     5         99999d-6*t10*t16*t31*t4-7.46904d-3*t10*t12*t31*t4-3.73452
     6         d-3*t10*t11*t15)-3.73452d-3*t11*t44*t45+6.973319815199999
     7         d-6*t17*t31/rhoa**4)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (1.0d+0*t18*t23*(2.0919959445599998d-
     1         5*t24*t62*t8/rhob**1.0666666666666666d+1-1.49380799999999
     2         98d-2*gammabb*t38*t46*t8+2*t25*t26*t8+2.0919959445599998d
     3         -5*gammabb*t12*t21*t46*t62+2.0919959445599998d-5*t13*t25*
     4         t27*t62+6.973319815199999d-6*gammabb*t15*t38*t46+6.973319
     5         815199999d-6*t16*t21*t25*t38-7.46904d-3*t12*t21*t25*t38-3
     6         .73452d-3*t15*t25*t26)-3.73452d-3*t26*t47*t48+6.973319815
     7         199999d-6*t28*t38/rhob**4)*wght+Cmat2(iq,D2_GBB_GBB)
            Mmat2(iq,D2_RA_TA) = (t6*(1.0d+0*t1*(4.97936d-2*t13*t31*t4*t
     1         50*taua+1.24484d-2*t11*t16*t50*taua-1.3333333333333333d+1
     2         *t11*t13*t50*taua+2.48968d-2*gammaaa*t12*t31*t49*taua+2.2
     3         40712d-2*t30*t49*t58*t8*t9+5.228328d-2*t31*t43*t8*t9-3.33
     4         33333333333337d+0*t16*t2*t6+2.240712d-2*gammaaa*t12*t30*t
     5         4*t50*t58+2.240712d-2*t13*t14*t3*t30*t58+7.46904d-3*gamma
     6         aa*t15*t30*t31*t50-4.0d+0*gammaaa*t12*t30*t31*t50+7.46904
     7         d-3*t16*t3*t30*t31*t4-8.0d+0*t13*t3*t30*t31*t4+3.236584d-
     8         2*gammaaa*t10*t12*t31*t4+6.2242d-3*t11*t16*t2*t4-6.666666
     9         666666667d+0*t11*t13*t2*t4+1.24484d-2*t13*t14*t2*t31-2.0d
     :         +0*t11*t16*t3*t30+1.618292d-2*gammaaa*t10*t11*t15-8.66666
     ;         6666666666d+0*gammaaa*t10*t11*t12)+1.3333333333333333d+0*
     <         t34*t51)+7.46904d-3*t17*t30*t31*t52-1.0d+0*t1*t11*t30*t51
     =         +6.2242d-3*t11*t17*t45-3.73452d-3*t11*t3*t35)*wght+Mmat2(
     >         iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (t23*(1.0d+0*t18*(4.97936d-2*t13*t21*t3
     1         8*t54*taub+1.24484d-2*t16*t26*t54*taub-1.3333333333333333
     2         d+1*t13*t26*t54*taub+2.48968d-2*gammabb*t12*t38*t53*taub+
     3         2.240712d-2*t24*t37*t53*t62*t8+5.228328d-2*t24*t38*t46*t8
     4         +2.240712d-2*gammabb*t12*t21*t37*t54*t62+2.240712d-2*t13*
     5         t20*t27*t37*t62+7.46904d-3*gammabb*t15*t37*t38*t54-4.0d+0
     6         *gammabb*t12*t37*t38*t54+7.46904d-3*t16*t20*t21*t37*t38-8
     7         .0d+0*t13*t20*t21*t37*t38+1.24484d-2*t13*t19*t27*t38+3.23
     8         6584d-2*gammabb*t12*t21*t25*t38-2.0d+0*t16*t20*t26*t37+1.
     9         618292d-2*gammabb*t15*t25*t26-8.666666666666666d+0*gammab
     :         b*t12*t25*t26+6.2242d-3*t16*t19*t21*t26-6.666666666666667
     ;         d+0*t13*t19*t21*t26-3.3333333333333337d+0*t16*t19*t23)+1.
     <         3333333333333333d+0*t41*t55)+7.46904d-3*t28*t37*t38*t56-1
     =         .0d+0*t18*t26*t37*t55+6.2242d-3*t26*t28*t48-3.73452d-3*t2
     >         0*t26*t42)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = (1.0d+0*t1*t6*(4.1839918891199995d-5*t
     1         58*t8*t9/rhoa**9.666666666666666d+0-1.4938079999999998d-2
     2         *gammaaa*t31*t49*t8+4.1839918891199995d-5*t13*t14*t50*t58
     3         +4.1839918891199995d-5*gammaaa*t12*t4*t49*t58+1.394663963
     4         0399998d-5*t16*t31*t4*t50-1.4938079999999998d-2*t13*t31*t
     5         4*t50-7.46904d-3*t12*t31*t4*t50-3.73452d-3*t11*t16*t50-3.
     6         73452d-3*t11*t15*t50+2.0d+0*t11*t12*t50+1.394663963039999
     7         8d-5*gammaaa*t15*t31*t49-7.46904d-3*gammaaa*t12*t31*t49)-
     8         3.73452d-3*t11*t44*t52-1.86726d-3*t11*t45*t51+1.394663963
     9         0399998d-5*t17*t31/rhoa**3)*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = (1.0d+0*t18*t23*(4.1839918891199995d-5
     1         *t24*t62*t8/rhob**9.666666666666666d+0-1.4938079999999998
     2         d-2*gammabb*t38*t53*t8+4.1839918891199995d-5*t13*t27*t54*
     3         t62+4.1839918891199995d-5*gammabb*t12*t21*t53*t62+1.39466
     4         39630399998d-5*t16*t21*t38*t54-1.4938079999999998d-2*t13*
     5         t21*t38*t54-7.46904d-3*t12*t21*t38*t54-3.73452d-3*t16*t26
     6         *t54-3.73452d-3*t15*t26*t54+2.0d+0*t12*t26*t54+1.39466396
     7         30399998d-5*gammabb*t15*t38*t53-7.46904d-3*gammabb*t12*t3
     8         8*t53)-3.73452d-3*t26*t47*t56-1.86726d-3*t26*t48*t55+1.39
     9         46639630399998d-5*t28*t38/rhob**3)*wght+Mmat2(iq,D2_GBB_T
     :         B)
            Mmat2(iq,D2_TA_TA) = (1.0d+0*t1*t6*(8.367983778239999d-5*t58
     1         *t8*t9/rhoa**8.666666666666666d+0+8.367983778239999d-5*t1
     2         3*t14*t58*t66+2.7893279260799997d-5*t16*t31*t4*t66-5.9752
     3         31999999999d-2*t13*t31*t4*t66-1.4938079999999998d-2*t11*t
     4         16*t66+8.0d+0*t11*t13*t66+8.367983778239999d-5*gammaaa*t1
     5         2*t4*t58*t65+2.7893279260799997d-5*gammaaa*t15*t31*t65-2.
     6         9876159999999996d-2*gammaaa*t12*t31*t65)-7.46904d-3*t11*t
     7         51*t52+2.7893279260799997d-5*t17*t31/rhoa**2)*wght+Mmat2(
     8         iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = (1.0d+0*t18*t23*(8.367983778239999d-5*t
     1         24*t62*t8/rhob**8.666666666666666d+0+8.367983778239999d-5
     2         *t13*t27*t62*t68+2.7893279260799997d-5*t16*t21*t38*t68-5.
     3         975231999999999d-2*t13*t21*t38*t68-1.4938079999999998d-2*
     4         t16*t26*t68+8.0d+0*t13*t26*t68+8.367983778239999d-5*gamma
     5         bb*t12*t21*t62*t67+2.7893279260799997d-5*gammabb*t15*t38*
     6         t67-2.9876159999999996d-2*gammabb*t12*t38*t67)-7.46904d-3
     7         *t26*t55*t56+2.7893279260799997d-5*t28*t38/rhob**2)*wght+
     8         Mmat2(iq,D2_TB_TB)
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
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammaaa*t10*t11*t9*taua-3.33333
     1         33333333337d+0*t15*t2*t6*taua-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taua-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gamm
     4         aaa*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammaaa*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhoa**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22*wght+1.0d+0*t1*t21*wght
            t24 = 1/rhoa**8
            t25 = 2*gammaaa*t10*t7*t9-3.73452d-3*gammaaa*t11*t19*t4*t9-1
     1         .86726d-3*gammaaa*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhoa**7
            t28 = 1/rhoa**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammaaa*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammaaa*
     3         t10*t14*t28+2.0d+0*gammaaa*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhoa**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taua+1.8257653333333332d-2*g
     1         ammaaa*t33
            t35 = 1/rhoa**6
            t36 = 1/rhoa**3.3333333333333337d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = -1.0d+0*t1*t10*t16*t18*wght+t23*t6+Amat(iq,
     1         D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t1*t25*t6*wght-1.86726d-3*t10*t16*t
     1         26*wght+Cmat(iq,D1_GAA)
            Mmat(iq,D1_TA) = 1.0d+0*t1*t29*t6*wght-3.73452d-3*t10*t16*t3
     1         0*wght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = t6*(1.0d+0*t1*(2.2222222222222224d+1*t1
     1         0*t12*t9*taua**2+1.3333333333333333d+1*gammaaa*t11*t18*t1
     2         9*t9*taua+8.88888888888889d+0*t15*t17*t6*taua+2.666666666
     3         6666666d+1*t12*t18*t19*t2*t4*taua+1.7777777777777778d+1*t
     4         10*t12*t17*t4*taua+2.6666666666666666d+1*gammaaa*t10*t11*
     5         t20*taua+6.666666666666667d+0*t10*t15*t18*t2*taua-2*t19*t
     6         34*t7*t8*t9+6*t31*t32*t7*t8*t9+2.1333333333333332d+1*t18*
     7         t19*t20*t7*t8+3.377777777777778d+1*t10*t7*t8/rhoa**7.3333
     8         33333333333d+0+9.777777777777779d+0*gammaaa*t14*t33*t6-2*
     9         gammaaa*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.77777777777777
     :         9d+0*gammaaa*t10*t11*t33*t4+6*gammaaa*t11*t2*t31*t32*t4+2
     ;         *t15*t19*t31*t4+1.0666666666666666d+1*gammaaa*t11*t17*t18
     <         *t19*t4-gammaaa*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t1
     =         3*t31*t32+2*gammaaa*t14*t19*t2*t31+5.333333333333333d+0*g
     >         ammaaa*t10*t14*t17*t18)*wght+2.6666666666666666d+0*t21*t2
     ?         2*wght+4.444444444444444d-1*t16*wght/rhoa**6.666666666666
     @         666d-1)+t10*(-1.0d+0*t1*t16*t34*wght-1.3333333333333333d+
     1         0*t16*t18*t22*wght-1.0d+0*t1*t18*t21*wght-t18*t23)+2.0d+0
     2         *t1*t16*t19*t31*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t6*(1.0d+0*t1*(2.48968d-2*t12*t19*t4*t
     1         9*taua+6.2242d-3*t10*t15*t9*taua-3.3333333333333337d+0*t1
     2         0*t11*t9*taua+1.24484d-2*gammaaa*t11*t19*t24*taua-4*gamma
     3         aa*t18*t19*t7*t9+1.120356d-2*gammaaa*t11*t18*t32*t4*t9+3.
     4         73452d-3*gammaaa*t14*t18*t19*t9+1.120356d-2*t18*t24*t32*t
     5         7*t8+2.9876159999999996d-2*t19*t7*t8/rhoa**9-1.0666666666
     6         666666d+1*gammaaa*t10*t20*t7-2.6666666666666666d+0*t14*t1
     7         7*t6+1.9917439999999997d-2*gammaaa*t11*t19*t20*t4+3.73452
     8         d-3*t15*t18*t19*t2*t4-2*t11*t18*t19*t2*t4+4.9793599999999
     9         993d-3*t10*t15*t17*t4-2.6666666666666666d+0*t10*t11*t17*t
     :         4+1.120356d-2*t12*t13*t18*t2*t32+9.958719999999999d-3*gam
     ;         maaa*t10*t14*t20-t10*t14*t18*t2+9.958719999999999d-3*t12*
     <         t13*t17*t19)*wght+1.3333333333333333d+0*t22*t25*wght)+3.7
     =         3452d-3*t16*t18*t19*t26*wght-1.0d+0*t1*t10*t18*t25*wght+4
     >         .9793599999999993d-3*t10*t16*wght/rhoa**2.333333333333333
     ?         4d+0-1.86726d-3*t10*t2*t23+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*t6*(2*t10*t7*t9+6.973319815
     1         199999d-6*t15*t19*t4*t9-7.46904d-3*t11*t19*t4*t9+2.091995
     2         9445599998d-5*t12*t13*t32*t9-3.73452d-3*t10*t14*t9+2.0919
     3         959445599998d-5*t32*t7*t8/rhoa**1.0666666666666666d+1-1.4
     4         938079999999998d-2*gammaaa*t19*t24*t7+2.0919959445599998d
     5         -5*gammaaa*t11*t24*t32*t4+6.973319815199999d-6*gammaaa*t1
     6         4*t19*t24)*wght-3.73452d-3*t10*t25*t26*wght+6.97331981519
     7         9999d-6*t16*t19*wght/rhoa**4+Cmat2(iq,D2_GAA_GAA)
            Mmat2(iq,D2_RA_TA) = t6*(1.0d+0*t1*(4.97936d-2*t12*t19*t28*t
     1         4*taua+1.24484d-2*t10*t15*t28*taua-1.3333333333333333d+1*
     2         t10*t12*t28*taua+2.48968d-2*gammaaa*t11*t19*t27*taua+3.23
     3         6584d-2*gammaaa*t11*t19*t4*t9+1.618292d-2*gammaaa*t10*t14
     4         *t9-8.666666666666666d+0*gammaaa*t10*t11*t9+2.240712d-2*t
     5         18*t27*t32*t7*t8+5.228328d-2*t19*t24*t7*t8-3.333333333333
     6         3337d+0*t15*t2*t6+2.240712d-2*gammaaa*t11*t18*t28*t32*t4+
     7         7.46904d-3*t15*t18*t19*t3*t4-8.0d+0*t12*t18*t19*t3*t4+6.2
     8         242d-3*t10*t15*t2*t4-6.666666666666667d+0*t10*t12*t2*t4+2
     9         .240712d-2*t12*t13*t18*t3*t32-2.0d+0*t10*t15*t18*t3+7.469
     :         04d-3*gammaaa*t14*t18*t19*t28-4.0d+0*gammaaa*t11*t18*t19*
     ;         t28+1.24484d-2*t12*t13*t19*t2)*wght+1.3333333333333333d+0
     <         *t22*t29*wght)+7.46904d-3*t16*t18*t19*t30*wght-1.0d+0*t1*
     =         t10*t18*t29*wght+6.2242d-3*t10*t16*t26*wght-3.73452d-3*t1
     >         0*t23*t3+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*t6*(4.1839918891199995d-5*t3
     1         2*t7*t8/rhoa**9.666666666666666d+0-1.4938079999999998d-2*
     2         gammaaa*t19*t27*t7+4.1839918891199995d-5*gammaaa*t11*t27*
     3         t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.49380799999
     4         99998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.18399
     5         18891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28-3.7
     6         3452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.3946639630399998
     7         d-5*gammaaa*t14*t19*t27-7.46904d-3*gammaaa*t11*t19*t27)*w
     8         ght-3.73452d-3*t10*t25*t30*wght-1.86726d-3*t10*t26*t29*wg
     9         ht+1.3946639630399998d-5*t16*t19*wght/rhoa**3+Mmat2(iq,D2
     :         _GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*t6*(8.367983778239999d-5*t32*
     1         t7*t8/rhoa**8.666666666666666d+0+2.7893279260799997d-5*t1
     2         5*t19*t36*t4-5.975231999999999d-2*t12*t19*t36*t4+8.367983
     3         778239999d-5*gammaaa*t11*t32*t35*t4+8.367983778239999d-5*
     4         t12*t13*t32*t36-1.4938079999999998d-2*t10*t15*t36+8.0d+0*
     5         t10*t12*t36+2.7893279260799997d-5*gammaaa*t14*t19*t35-2.9
     6         876159999999996d-2*gammaaa*t11*t19*t35)*wght-7.46904d-3*t
     7         10*t29*t30*wght+2.7893279260799997d-5*t16*t19*wght/rhoa**
     8         2+Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
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
            t20 = 1/rhob**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammabb*t10*t11*t9*taub-3.33333
     1         33333333337d+0*t15*t2*t6*taub-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taub-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammabb*t14*t17*t6-2*gamm
     4         abb*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammabb*t10*t11*t17*t4-gammabb*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhob**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22*wght+1.0d+0*t1*t21*wght
            t24 = 1/rhob**8
            t25 = 2*gammabb*t10*t7*t9-3.73452d-3*gammabb*t11*t19*t4*t9-1
     1         .86726d-3*gammabb*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhob**7
            t28 = 1/rhob**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammabb*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammabb*
     3         t10*t14*t28+2.0d+0*gammabb*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhob**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taub+1.8257653333333332d-2*g
     1         ammabb*t33
            t35 = 1/rhob**6
            t36 = 1/rhob**3.3333333333333337d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = -1.0d+0*t1*t10*t16*t18*wght+t23*t6+Amat(iq,
     1         D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t1*t25*t6*wght-1.86726d-3*t10*t16*t
     1         26*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TB) = 1.0d+0*t1*t29*t6*wght-3.73452d-3*t10*t16*t3
     1         0*wght+Mmat(iq,D1_TB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = t6*(1.0d+0*t1*(2.2222222222222224d+1*t1
     1         0*t12*t9*taub**2+1.3333333333333333d+1*gammabb*t11*t18*t1
     2         9*t9*taub+8.88888888888889d+0*t15*t17*t6*taub+2.666666666
     3         6666666d+1*t12*t18*t19*t2*t4*taub+1.7777777777777778d+1*t
     4         10*t12*t17*t4*taub+2.6666666666666666d+1*gammabb*t10*t11*
     5         t20*taub+6.666666666666667d+0*t10*t15*t18*t2*taub-2*t19*t
     6         34*t7*t8*t9+6*t31*t32*t7*t8*t9+2.1333333333333332d+1*t18*
     7         t19*t20*t7*t8+3.377777777777778d+1*t10*t7*t8/rhob**7.3333
     8         33333333333d+0+9.777777777777779d+0*gammabb*t14*t33*t6-2*
     9         gammabb*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.77777777777777
     :         9d+0*gammabb*t10*t11*t33*t4+6*gammabb*t11*t2*t31*t32*t4+2
     ;         *t15*t19*t31*t4+1.0666666666666666d+1*gammabb*t11*t17*t18
     <         *t19*t4-gammabb*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t1
     =         3*t31*t32+2*gammabb*t14*t19*t2*t31+5.333333333333333d+0*g
     >         ammabb*t10*t14*t17*t18)*wght+2.6666666666666666d+0*t21*t2
     ?         2*wght+4.444444444444444d-1*t16*wght/rhob**6.666666666666
     @         666d-1)+t10*(-1.0d+0*t1*t16*t34*wght-1.3333333333333333d+
     1         0*t16*t18*t22*wght-1.0d+0*t1*t18*t21*wght-t18*t23)+2.0d+0
     2         *t1*t16*t19*t31*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t6*(1.0d+0*t1*(2.48968d-2*t12*t19*t4*t
     1         9*taub+6.2242d-3*t10*t15*t9*taub-3.3333333333333337d+0*t1
     2         0*t11*t9*taub+1.24484d-2*gammabb*t11*t19*t24*taub-4*gamma
     3         bb*t18*t19*t7*t9+1.120356d-2*gammabb*t11*t18*t32*t4*t9+3.
     4         73452d-3*gammabb*t14*t18*t19*t9+1.120356d-2*t18*t24*t32*t
     5         7*t8+2.9876159999999996d-2*t19*t7*t8/rhob**9-1.0666666666
     6         666666d+1*gammabb*t10*t20*t7-2.6666666666666666d+0*t14*t1
     7         7*t6+1.9917439999999997d-2*gammabb*t11*t19*t20*t4+3.73452
     8         d-3*t15*t18*t19*t2*t4-2*t11*t18*t19*t2*t4+4.9793599999999
     9         993d-3*t10*t15*t17*t4-2.6666666666666666d+0*t10*t11*t17*t
     :         4+1.120356d-2*t12*t13*t18*t2*t32+9.958719999999999d-3*gam
     ;         mabb*t10*t14*t20-t10*t14*t18*t2+9.958719999999999d-3*t12*
     <         t13*t17*t19)*wght+1.3333333333333333d+0*t22*t25*wght)+3.7
     =         3452d-3*t16*t18*t19*t26*wght-1.0d+0*t1*t10*t18*t25*wght+4
     >         .9793599999999993d-3*t10*t16*wght/rhob**2.333333333333333
     ?         4d+0-1.86726d-3*t10*t2*t23+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*t6*(2*t10*t7*t9+6.973319815
     1         199999d-6*t15*t19*t4*t9-7.46904d-3*t11*t19*t4*t9+2.091995
     2         9445599998d-5*t12*t13*t32*t9-3.73452d-3*t10*t14*t9+2.0919
     3         959445599998d-5*t32*t7*t8/rhob**1.0666666666666666d+1-1.4
     4         938079999999998d-2*gammabb*t19*t24*t7+2.0919959445599998d
     5         -5*gammabb*t11*t24*t32*t4+6.973319815199999d-6*gammabb*t1
     6         4*t19*t24)*wght-3.73452d-3*t10*t25*t26*wght+6.97331981519
     7         9999d-6*t16*t19*wght/rhob**4+Cmat2(iq,D2_GBB_GBB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = t6*(1.0d+0*t1*(4.97936d-2*t12*t19*t28*t
     1         4*taub+1.24484d-2*t10*t15*t28*taub-1.3333333333333333d+1*
     2         t10*t12*t28*taub+2.48968d-2*gammabb*t11*t19*t27*taub+3.23
     3         6584d-2*gammabb*t11*t19*t4*t9+1.618292d-2*gammabb*t10*t14
     4         *t9-8.666666666666666d+0*gammabb*t10*t11*t9+2.240712d-2*t
     5         18*t27*t32*t7*t8+5.228328d-2*t19*t24*t7*t8-3.333333333333
     6         3337d+0*t15*t2*t6+2.240712d-2*gammabb*t11*t18*t28*t32*t4+
     7         7.46904d-3*t15*t18*t19*t3*t4-8.0d+0*t12*t18*t19*t3*t4+6.2
     8         242d-3*t10*t15*t2*t4-6.666666666666667d+0*t10*t12*t2*t4+2
     9         .240712d-2*t12*t13*t18*t3*t32-2.0d+0*t10*t15*t18*t3+7.469
     :         04d-3*gammabb*t14*t18*t19*t28-4.0d+0*gammabb*t11*t18*t19*
     ;         t28+1.24484d-2*t12*t13*t19*t2)*wght+1.3333333333333333d+0
     <         *t22*t29*wght)+7.46904d-3*t16*t18*t19*t30*wght-1.0d+0*t1*
     =         t10*t18*t29*wght+6.2242d-3*t10*t16*t26*wght-3.73452d-3*t1
     >         0*t23*t3+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*t6*(4.1839918891199995d-5*t3
     1         2*t7*t8/rhob**9.666666666666666d+0-1.4938079999999998d-2*
     2         gammabb*t19*t27*t7+4.1839918891199995d-5*gammabb*t11*t27*
     3         t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.49380799999
     4         99998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.18399
     5         18891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28-3.7
     6         3452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.3946639630399998
     7         d-5*gammabb*t14*t19*t27-7.46904d-3*gammabb*t11*t19*t27)*w
     8         ght-3.73452d-3*t10*t25*t30*wght-1.86726d-3*t10*t26*t29*wg
     9         ht+1.3946639630399998d-5*t16*t19*wght/rhob**3+Mmat2(iq,D2
     :         _GBB_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*t6*(8.367983778239999d-5*t32*
     1         t7*t8/rhob**8.666666666666666d+0+2.7893279260799997d-5*t1
     2         5*t19*t36*t4-5.975231999999999d-2*t12*t19*t36*t4+8.367983
     3         778239999d-5*gammabb*t11*t32*t35*t4+8.367983778239999d-5*
     4         t12*t13*t32*t36-1.4938079999999998d-2*t10*t15*t36+8.0d+0*
     5         t10*t12*t36+2.7893279260799997d-5*gammabb*t14*t19*t35-2.9
     6         876159999999996d-2*gammabb*t11*t19*t35)*wght-7.46904d-3*t
     7         10*t29*t30*wght+2.7893279260799997d-5*t16*t19*wght/rhob**
     8         2+Mmat2(iq,D2_TB_TB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_vs98 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{2.0\,\tau_\alpha}\over{\rho_\alpha^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_3} &=& 0.00186726\,{\it t_2}+0.00186726\,{\it t_1}
C>    \,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_3}}}\\\\ 
C>   {\it t_5} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_6} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_3}^2}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_9} &=& {\it param}\left(6\right)\\\\ 
C>   {\it t_{10}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{11}} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_{12}} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{13}} &=& {{2.0\,\tau_\beta}\over{\rho_\beta^{{{5}
C>    \over{3}}}}}-9.11559972\\\\ 
C>   {\it t_{14}} &=& 0.00186726\,{\it t_{13}}+0.00186726\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_{14}}^2}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_{18}} &=& {{2.0\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}
C>    -9.11559972\\\\ 
C>   {\it t_{19}} &=& 0.00186726\,{\it t_{18}}+0.00186726\,{
C>    \it t_{17}}\,\sigma_{ss}+1.0\\\\ 
C>   {\it t_{20}} &=& {{1}\over{{\it t_{19}}}}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{19}}^2}}\\\\ 
C>   f &=& 1.0\,\rho_\beta^{{{4}\over{3}}}\,{\it t_{15}}\,\left({
C>    \it t_{11}}\,{\it t_{13}}\,{\it t_{15}}+{\it t_{10}}\,{
C>    \it t_{12}}\,\sigma_{\beta\beta}\,{\it t_{15}}+{\it t_9}\,{
C>    \it t_{13}}^2\,{\it t_{16}}+{\it t_8}\,{\it t_{12}}\,
C>    \sigma_{\beta\beta}\,{\it t_{13}}\,{\it t_{16}}+{{{\it t_6}
C>    \,\sigma_{\beta\beta}^2\,{\it t_{16}}}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_5}\right)+1.0
C>    \,\rho_\alpha^{{{4}\over{3}}}\,{\it t_4}\,\left({\it t_{11}}
C>    \,{\it t_2}\,{\it t_4}+{\it t_{10}}\,{\it t_1}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{\it t_9}\,{\it t_2}^2\,{
C>    \it t_7}+{\it t_8}\,{\it t_1}\,\sigma_{\alpha\alpha}\,{
C>    \it t_2}\,{\it t_7}+{{{\it t_6}\,\sigma_{\alpha\alpha}^2\,{
C>    \it t_7}}\over{\rho_\alpha^{{{16}\over{3}}}}}+{
C>    \it t_5}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 1.0\,\rho_s^{{{4}\over{3}}}\,{\it t_{20}}\,\left({
C>    \it t_{11}}\,{\it t_{18}}\,{\it t_{20}}+{\it t_{10}}\,{
C>    \it t_{17}}\,\sigma_{ss}\,{\it t_{20}}+{\it t_9}\,{
C>    \it t_{18}}^2\,{\it t_{21}}+{\it t_8}\,{\it t_{17}}\,
C>    \sigma_{ss}\,{\it t_{18}}\,{\it t_{21}}+{{{\it t_6}\,
C>    \sigma_{ss}^2\,{\it t_{21}}}\over{\rho_s^{{{16}\over{3}}}}}
C>    +{\it t_5}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T van Voorhis, GE Scuseria, J.Chem.Phys. 109, 400 (1998)  , DOI: 
C> <a href="https://doi.org/10.1063/1.3005348 "> 
C> 10.1063/1.3005348 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_vs98_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
CDIR$ NOVECTOR
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
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammaaa*t10*t11*t9*taua-3.33333
     1         33333333337d+0*t15*t2*t6*taua-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taua-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gamm
     4         aaa*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammaaa*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhoa**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22+1.0d+0*t1*t21
            t24 = 1/rhoa**8
            t25 = 2*gammaaa*t10*t7*t9-3.73452d-3*gammaaa*t11*t19*t4*t9-1
     1         .86726d-3*gammaaa*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhoa**7
            t28 = 1/rhoa**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammaaa*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammaaa*
     3         t10*t14*t28+2.0d+0*gammaaa*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhoa**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taua+1.8257653333333332d-2*g
     1         ammaaa*t33
            t35 = 1/rhoa**7.333333333333333d+0
            t36 = taua**2
            t37 = 1.3333333333333333d+1*gammaaa*t11*t18*t19*t9*taua+8.88
     1         888888888889d+0*t15*t17*t6*taua+2.6666666666666666d+1*t12
     2         *t18*t19*t2*t4*taua+1.7777777777777778d+1*t10*t12*t17*t4*
     3         taua+2.6666666666666666d+1*gammaaa*t10*t11*t20*taua+6.666
     4         666666666667d+0*t10*t15*t18*t2*taua-2*t19*t34*t7*t8*t9+6*
     5         t31*t32*t7*t8*t9+2.2222222222222224d+1*t10*t12*t36*t9+3.3
     6         77777777777778d+1*t10*t35*t7*t8+2.1333333333333332d+1*t18
     7         *t19*t20*t7*t8+9.777777777777779d+0*gammaaa*t14*t33*t6-2*
     8         gammaaa*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.77777777777777
     9         9d+0*gammaaa*t10*t11*t33*t4+6*gammaaa*t11*t2*t31*t32*t4+2
     :         *t15*t19*t31*t4+1.0666666666666666d+1*gammaaa*t11*t17*t18
     ;         *t19*t4-gammaaa*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t1
     <         3*t31*t32+2*gammaaa*t14*t19*t2*t31+5.333333333333333d+0*g
     =         ammaaa*t10*t14*t17*t18
            t38 = 1/rhoa**6.666666666666666d-1
            t39 = 4.444444444444444d-1*t16*t38+1.0d+0*t1*t37+2.666666666
     1         6666666d+0*t21*t22
            t40 = -1.0d+0*t1*t16*t34-t18*t23-1.3333333333333333d+0*t16*t
     1         18*t22-1.0d+0*t1*t18*t21
            t41 = 1/rhoa**2.3333333333333334d+0
            t42 = 1/rhoa**9
            t43 = 2.48968d-2*t12*t19*t4*t9*taua+6.2242d-3*t10*t15*t9*tau
     1         a-3.3333333333333337d+0*t10*t11*t9*taua+1.24484d-2*gammaa
     2         a*t11*t19*t24*taua-4*gammaaa*t18*t19*t7*t9+1.120356d-2*ga
     3         mmaaa*t11*t18*t32*t4*t9+3.73452d-3*gammaaa*t14*t18*t19*t9
     4         +2.9876159999999996d-2*t19*t42*t7*t8+1.120356d-2*t18*t24*
     5         t32*t7*t8-1.0666666666666666d+1*gammaaa*t10*t20*t7-2.6666
     6         666666666666d+0*t14*t17*t6+1.9917439999999997d-2*gammaaa*
     7         t11*t19*t20*t4+3.73452d-3*t15*t18*t19*t2*t4-2*t11*t18*t19
     8         *t2*t4+4.9793599999999993d-3*t10*t15*t17*t4-2.66666666666
     9         66666d+0*t10*t11*t17*t4+1.120356d-2*t12*t13*t18*t2*t32+9.
     :         958719999999999d-3*gammaaa*t10*t14*t20-t10*t14*t18*t2+9.9
     ;         58719999999999d-3*t12*t13*t17*t19
            t44 = 1.0d+0*t1*t43+1.3333333333333333d+0*t22*t25
            t45 = 1/rhoa**1.0666666666666666d+1
            t46 = 2*t10*t7*t9+6.973319815199999d-6*t15*t19*t4*t9-7.46904
     1         d-3*t11*t19*t4*t9+2.0919959445599998d-5*t12*t13*t32*t9-3.
     2         73452d-3*t10*t14*t9+2.0919959445599998d-5*t32*t45*t7*t8-1
     3         .4938079999999998d-2*gammaaa*t19*t24*t7+2.091995944559999
     4         8d-5*gammaaa*t11*t24*t32*t4+6.973319815199999d-6*gammaaa*
     5         t14*t19*t24
            t47 = 1/rhoa**4
            t48 = 4.97936d-2*t12*t19*t28*t4*taua+1.24484d-2*t10*t15*t28*
     1         taua-1.3333333333333333d+1*t10*t12*t28*taua+2.48968d-2*ga
     2         mmaaa*t11*t19*t27*taua+3.236584d-2*gammaaa*t11*t19*t4*t9+
     3         1.618292d-2*gammaaa*t10*t14*t9-8.666666666666666d+0*gamma
     4         aa*t10*t11*t9+2.240712d-2*t18*t27*t32*t7*t8+5.228328d-2*t
     5         19*t24*t7*t8-3.3333333333333337d+0*t15*t2*t6+2.240712d-2*
     6         gammaaa*t11*t18*t28*t32*t4+7.46904d-3*t15*t18*t19*t3*t4-8
     7         .0d+0*t12*t18*t19*t3*t4+6.2242d-3*t10*t15*t2*t4-6.6666666
     8         66666667d+0*t10*t12*t2*t4+2.240712d-2*t12*t13*t18*t3*t32-
     9         2.0d+0*t10*t15*t18*t3+7.46904d-3*gammaaa*t14*t18*t19*t28-
     :         4.0d+0*gammaaa*t11*t18*t19*t28+1.24484d-2*t12*t13*t19*t2
            t49 = 1.0d+0*t1*t48+1.3333333333333333d+0*t22*t29
            t50 = 1/rhoa**9.666666666666666d+0
            t51 = 4.1839918891199995d-5*t32*t50*t7*t8-1.4938079999999998
     1         d-2*gammaaa*t19*t27*t7+4.1839918891199995d-5*gammaaa*t11*
     2         t27*t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.4938079
     3         999999998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.1
     4         839918891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28
     5         -3.73452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.394663963039
     6         9998d-5*gammaaa*t14*t19*t27-7.46904d-3*gammaaa*t11*t19*t2
     7         7
            t52 = 1/rhoa**3
            t53 = 1/rhoa**6
            t54 = 1/rhoa**3.3333333333333337d+0
            t55 = 8.367983778239999d-5*t32*t7*t8/rhoa**8.666666666666666
     1         d+0+2.7893279260799997d-5*t15*t19*t4*t54-5.97523199999999
     2         9d-2*t12*t19*t4*t54+8.367983778239999d-5*t12*t13*t32*t54-
     3         1.4938079999999998d-2*t10*t15*t54+8.0d+0*t10*t12*t54+8.36
     4         7983778239999d-5*gammaaa*t11*t32*t4*t53+2.789327926079999
     5         7d-5*gammaaa*t14*t19*t53-2.9876159999999996d-2*gammaaa*t1
     6         1*t19*t53
            t56 = 1/rhoa**2
            t57 = t18**3
            t58 = 1/t5**5
            t59 = 1/rhoa**5.666666666666667d+0
            t60 = -6.085884444444444d-2*t33*taua-8.520238222222222d-2*ga
     1         mmaaa*t59
            t61 = 1/rhoa**5
            t62 = 1/rhoa**7.666666666666667d+0
            fnc(iq) = 2.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = (t23*t6-1.0d+0*t1*t10*t16*t18)*wght+Amat(iq
     1         ,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t25*t6-1.86726d-3*t10*t16*t26)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t29*t6-3.73452d-3*t10*t16*t30)*w
     1         ght+Mmat(iq,D1_TA)
            Amat2(iq,D2_RA_RA) = (t39*t6+t10*t40+2.0d+0*t1*t16*t19*t31)*
     1         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t44*t6+4.9793599999999993d-3*t10*t16*
     1         t41+3.73452d-3*t16*t18*t19*t26-1.0d+0*t1*t10*t18*t25-1.86
     2         726d-3*t10*t2*t23)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (1.0d+0*t1*t46*t6+6.973319815199999d-
     1         6*t16*t19*t47-3.73452d-3*t10*t25*t26)*wght+Cmat2(iq,D2_GA
     2         A_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (t49*t6+7.46904d-3*t16*t18*t19*t30-3.73
     1         452d-3*t10*t23*t3-1.0d+0*t1*t10*t18*t29+6.2242d-3*t10*t16
     2         *t26)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = (1.0d+0*t1*t51*t6+1.3946639630399998d-
     1         5*t16*t19*t52-3.73452d-3*t10*t25*t30-1.86726d-3*t10*t26*t
     2         29)*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = (1.0d+0*t1*t55*t6+2.7893279260799997d-5
     1         *t16*t19*t56-7.46904d-3*t10*t29*t30)*wght+Mmat2(iq,D2_TA_
     2         TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Amat3(iq,D3_RA_RA_RA) = (t6*(1.0d+0*t1*(2.0d+1*gammaaa*t11*t
     1         19*t34*t9*taua-6.0d+1*gammaaa*t11*t31*t32*t9*taua-3.25925
     2         92592592595d+1*t15*t33*t6*taua+4.0d+1*t12*t19*t2*t34*t4*t
     3         aua-6.518518518518519d+1*t10*t12*t33*t4*taua-1.2d+2*t12*t
     4         2*t31*t32*t4*taua-1.0666666666666667d+2*t12*t17*t18*t19*t
     5         4*taua-2.0148148148148146d+2*gammaaa*t10*t11*t35*taua+1.0
     6         d+1*t10*t15*t2*t34*taua-2.0d+1*t15*t19*t2*t31*taua-1.6d+2
     7         *gammaaa*t11*t18*t19*t20*taua-2.6666666666666666d+1*t10*t
     8         15*t17*t18*taua-2*t19*t60*t7*t8*t9-24*t57*t58*t7*t8*t9+18
     9         *t18*t32*t34*t7*t8*t9-1.3333333333333333d+2*t12*t18*t19*t
     :         36*t9-2.0266666666666666d+2*t18*t19*t35*t7*t8+32*t19*t20*
     ;         t34*t7*t8-96*t20*t31*t32*t7*t8-2.477037037037037d+2*t10*t
     <         7*t8/rhoa**8.333333333333334d+0-2*gammaaa*t11*t19*t2*t4*t
     =         60-t10*t15*t4*t60-gammaaa*t10*t14*t2*t60-2*t12*t13*t19*t6
     >         0-4.5629629629629626d+1*gammaaa*t14*t59*t6-4.562962962962
     ?         9626d+1*gammaaa*t10*t11*t4*t59-24*gammaaa*t11*t2*t4*t57*t
     @         58-24*t12*t13*t57*t58-6*t15*t32*t4*t57-6*gammaaa*t14*t2*t
     1         32*t57+18*gammaaa*t11*t18*t2*t32*t34*t4+6*t15*t18*t19*t34
     2         *t4+16*gammaaa*t11*t17*t19*t34*t4-5.866666666666667d+1*ga
     3         mmaaa*t11*t18*t19*t33*t4-48*gammaaa*t11*t17*t31*t32*t4-1.
     4         777777777777778d+2*t10*t12*t20*t36+18*t12*t13*t18*t32*t34
     5         +6*gammaaa*t14*t18*t19*t2*t34+8*gammaaa*t10*t14*t17*t34-2
     6         .9333333333333333d+1*gammaaa*t10*t14*t18*t33-16*gammaaa*t
     7         14*t17*t19*t31)+1.3333333333333333d+0*t21*t38+4.0d+0*t22*
     8         t37-2.962962962962963d-1*t16*t3)+t10*(-1.0d+0*t1*t16*t60-
     9         2*t18*t39-4.444444444444444d-1*t16*t18*t38-1.0d+0*t1*t18*
     :         t37-t23*t34-2.6666666666666666d+0*t16*t22*t34-2.0d+0*t1*t
     ;         21*t34-2.6666666666666666d+0*t18*t21*t22)-6.0d+0*t1*t16*t
     <         32*t57+t19*(-2*t18*t40+4.0d+0*t1*t16*t18*t34+2.6666666666
     =         666666d+0*t16*t22*t31+2.0d+0*t1*t21*t31))*wght+Amat3(iq,D
     >         3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t6*(1.0d+0*t1*(-1.493808d-1*t12*t1
     1         8*t32*t4*t9*taua-2.48968d-2*t15*t18*t19*t9*taua+1.3333333
     2         333333333d+1*t11*t18*t19*t9*taua-1.6597866666666666d-1*ga
     3         mmaaa*t11*t19*t42*taua-1.9917439999999997d-1*t12*t19*t20*
     4         t4*taua-7.46904d-2*gammaaa*t11*t18*t24*t32*taua-4.9793599
     5         999999993d-2*t10*t15*t20*taua+2.6666666666666666d+1*t10*t
     6         11*t20*taua-4*gammaaa*t19*t34*t7*t9+12*gammaaa*t31*t32*t7
     7         *t9-4.481424d-2*gammaaa*t11*t31*t4*t58*t9+1.120356d-2*gam
     8         maaa*t11*t32*t34*t4*t9+3.73452d-3*gammaaa*t14*t19*t34*t9-
     9         1.120356d-2*gammaaa*t14*t31*t32*t9-4.481424d-2*t24*t31*t5
     :         8*t7*t8-1.7925696d-1*t18*t32*t42*t7*t8+1.120356d-2*t24*t3
     ;         2*t34*t7*t8-2.6888544000000003d-1*t19*t7*t8/rhoa**10+6.75
     <         5555555555556d+1*gammaaa*t10*t35*t7+4.2666666666666664d+1
     =         *gammaaa*t18*t19*t20*t7+9.777777777777779d+0*t14*t33*t6-4
     >         .481424d-2*t12*t13*t2*t31*t58-1.2614378666666665d-1*gamma
     ?         aa*t11*t19*t35*t4+3.73452d-3*t15*t19*t2*t34*t4-2*t11*t19*
     @         t2*t34*t4-1.8257653333333332d-2*t10*t15*t33*t4+9.77777777
     1         7777779d+0*t10*t11*t33*t4-1.120356d-2*t15*t2*t31*t32*t4+6
     2         *t11*t2*t31*t32*t4-1.1950464d-1*gammaaa*t11*t18*t20*t32*t
     3         4-1.9917439999999997d-2*t15*t17*t18*t19*t4+1.066666666666
     4         6666d+1*t11*t17*t18*t19*t4-8.298933333333335d-2*t12*t19*t
     5         24*t36-6.307189333333332d-2*gammaaa*t10*t14*t35+1.120356d
     6         -2*t12*t13*t2*t32*t34-t10*t14*t2*t34-3.6515306666666664d-
     7         2*t12*t13*t19*t33-5.975231999999999d-2*t12*t13*t17*t18*t3
     8         2+2*t14*t19*t2*t31-3.9834879999999995d-2*gammaaa*t14*t18*
     9         t19*t20+5.333333333333333d+0*t10*t14*t17*t18)+2.666666666
     :         6666666d+0*t22*t43+4.444444444444444d-1*t25*t38)+t10*(-1.
     ;         1618506666666667d-2*t16*t54-t18*t44-1.0d+0*t1*t18*t43+4.9
     <         793599999999993d-3*t21*t41-1.0d+0*t1*t25*t34-1.3333333333
     =         333333d+0*t18*t22*t25+4.9793599999999993d-3*t17*t23)-1.99
     >         17439999999997d-2*t16*t18*t19*t41-3.73452d-3*t19*t2*t40-1
     ?         .86726d-3*t10*t2*t39-1.120356d-2*t16*t26*t31*t32+2.0d+0*t
     @         1*t19*t25*t31)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t6*(1.0d+0*t1*(-6.9733198152d-5*g
     1         ammaaa*t11*t32*t45*taua-1.3946639630399998d-4*t12*t24*t32
     2         *t4*taua-2.3244399384d-5*t15*t19*t24*taua+2.48968d-2*t11*
     3         t19*t24*taua-4*t18*t19*t7*t9-8.367983778239999d-5*t12*t13
     4         *t18*t58*t9-2.0919959445599998d-5*t15*t18*t32*t4*t9+2.240
     5         712d-2*t11*t18*t32*t4*t9+7.46904d-3*t14*t18*t19*t9-8.3679
     6         83778239999d-5*t18*t45*t58*t7*t8-2.2314623408639997d-4*t3
     7         2*t7*t8/rhoa**1.1666666666666665d+1+1.1950464d-1*gammaaa*
     8         t19*t42*t7+4.481424d-2*gammaaa*t18*t24*t32*t7-1.066666666
     9         6666666d+1*t10*t20*t7-8.367983778239999d-5*gammaaa*t11*t1
     :         8*t24*t4*t58-1.6735967556479997d-4*gammaaa*t11*t32*t4*t42
     ;         -5.578655852159999d-5*gammaaa*t14*t19*t42-3.7191039014399
     <         996d-5*t15*t19*t20*t4+3.9834879999999995d-2*t11*t19*t20*t
     =         4-2.0919959445599998d-5*gammaaa*t14*t18*t24*t32-1.1157311
     >         704319998d-4*t12*t13*t20*t32+1.9917439999999997d-2*t10*t1
     ?         4*t20)+1.3333333333333333d+0*t22*t46)+6.973319815199999d-
     @         6*t19*t23*t9-3.7191039014399996d-5*t16*t19*t61-2.09199594
     1         45599998d-5*t16*t18*t32*t47-1.0d+0*t1*t10*t18*t46-3.73452
     2         d-3*t10*t2*t44+9.958719999999999d-3*t10*t25*t41+7.46904d-
     3         3*t18*t19*t25*t26)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (1.0d+0*t1*t6*(-1.562520138975642
     1         d-7*t58*t7*t8/rhoa**1.3333333333333333d+1+1.255197566736d
     2         -4*gammaaa*t32*t45*t7-2.240712d-2*t19*t24*t7-1.5625201389
     3         75642d-7*gammaaa*t11*t4*t45*t58-1.562520138975642d-7*t12*
     4         t13*t24*t58-3.9063003474391045d-8*gammaaa*t14*t32*t45-3.9
     5         063003474391045d-8*t15*t24*t32*t4+6.275987833679999d-5*t1
     6         1*t24*t32*t4+2.0919959445599995d-5*t14*t19*t24)+2.0919959
     7         445599995d-5*t19*t25*t47-5.60178d-3*t10*t26*t46-3.9063003
     8         474391045d-8*t16*t32/rhoa**6.666666666666667d+0)*wght+Cma
     9         t3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Mmat3(iq,D3_RA_RA_TA) = (t6*(1.0d+0*t1*(-2.9876159999999996d
     1         -1*t12*t19*t4*t9*taua-7.469039999999999d-2*t10*t15*t9*tau
     2         a+8.0d+1*t10*t12*t9*taua-2.987616d-1*t12*t18*t28*t32*t4*t
     3         aua-1.493808d-1*gammaaa*t11*t18*t27*t32*taua-4.97936d-2*t
     4         15*t18*t19*t28*taua+5.333333333333333d+1*t12*t18*t19*t28*
     5         taua-2.8216373333333333d-1*gammaaa*t11*t19*t24*taua-1.941
     6         9503999999999d-1*gammaaa*t11*t18*t32*t4*t9-6.473168d-2*ga
     7         mmaaa*t14*t18*t19*t9+3.466666666666666d+1*gammaaa*t11*t18
     8         *t19*t9-8.962848d-2*t27*t31*t58*t7*t8-4.1826624d-1*t19*t4
     9         2*t7*t8+2.240712d-2*t27*t32*t34*t7*t8-3.1369968d-1*t18*t2
     :         4*t32*t7*t8+8.88888888888889d+0*t15*t17*t6-8.962848d-2*ga
     ;         mmaaa*t11*t28*t31*t4*t58-8.962848d-2*t12*t13*t3*t31*t58+2
     <         .240712d-2*gammaaa*t11*t28*t32*t34*t4+7.46904d-3*t15*t19*
     =         t3*t34*t4-8.0d+0*t12*t19*t3*t34*t4-2.240712d-2*t15*t3*t31
     >         *t32*t4+2.4d+1*t12*t3*t31*t32*t4-1.7261781333333331d-1*ga
     ?         mmaaa*t11*t19*t20*t4-2.48968d-2*t15*t18*t19*t2*t4+2.66666
     @         66666666666d+1*t12*t18*t19*t2*t4-1.6597866666666666d-2*t1
     1         0*t15*t17*t4+1.7777777777777778d+1*t10*t12*t17*t4-1.65978
     2         6666666667d-1*t12*t19*t27*t36+2.240712d-2*t12*t13*t3*t32*
     3         t34-2.0d+0*t10*t15*t3*t34+7.46904d-3*gammaaa*t14*t19*t28*
     4         t34-4.0d+0*gammaaa*t11*t19*t28*t34-2.240712d-2*gammaaa*t1
     5         4*t28*t31*t32+1.2d+1*gammaaa*t11*t28*t31*t32-7.46904d-2*t
     6         12*t13*t18*t2*t32+4.0d+0*t15*t19*t3*t31-8.630890666666666
     7         d-2*gammaaa*t10*t14*t20+4.622222222222223d+1*gammaaa*t10*
     8         t11*t20+6.666666666666667d+0*t10*t15*t18*t2-3.31957333333
     9         3333d-2*t12*t13*t17*t19)+2.6666666666666666d+0*t22*t48+4.
     :         444444444444444d-1*t29*t38)+t10*(-t18*t49-1.0d+0*t1*t18*t
     ;         48-8.298933333333333d-3*t16*t41-1.0d+0*t1*t29*t34-1.33333
     <         33333333333d+0*t18*t22*t29+6.2242d-3*t21*t26+6.2242d-3*t2
     =         *t23)-7.46904d-3*t19*t3*t40-3.73452d-3*t10*t3*t39-2.24071
     >         2d-2*t16*t30*t31*t32+2.0d+0*t1*t19*t29*t31-2.48968d-2*t16
     ?         *t18*t19*t26)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_GAA_TA) = (t6*(1.0d+0*t1*(-1.3946639630399998
     1         d-4*gammaaa*t11*t32*t50*taua-2.7893279260799997d-4*t12*t2
     2         7*t32*t4*taua-4.6488798768d-5*t15*t19*t27*taua+4.97936d-2
     3         *t12*t19*t27*taua+2.48968d-2*t11*t19*t27*taua-6.043543839
     4         839999d-5*t15*t19*t4*t9+6.473168d-2*t12*t19*t4*t9+3.23658
     5         4d-2*t11*t19*t4*t9-1.8130631519519996d-4*t12*t13*t32*t9+1
     6         .618292d-2*t10*t15*t9+1.618292d-2*t10*t14*t9-8.6666666666
     7         66666d+0*t10*t11*t9-1.673596755648d-4*t18*t50*t58*t7*t8-4
     8         .044525492816d-4*t32*t45*t7*t8+4.481424d-2*gammaaa*t18*t2
     9         7*t32*t7+1.0456655999999999d-1*gammaaa*t19*t24*t7-1.67359
     :         6755648d-4*gammaaa*t11*t18*t27*t4*t58-1.673596755648d-4*t
     ;         12*t13*t18*t28*t58-4.1839918891199995d-5*t15*t18*t28*t32*
     <         t4+4.481424d-2*t12*t18*t28*t32*t4+2.240712d-2*t11*t18*t28
     =         *t32*t4-2.9287943223839996d-4*gammaaa*t11*t24*t32*t4-4.18
     >         39918891199995d-5*gammaaa*t14*t18*t27*t32+2.240712d-2*gam
     ?         maaa*t11*t18*t27*t32+7.46904d-3*t15*t18*t19*t28+7.46904d-
     @         3*t14*t18*t19*t28-4.0d+0*t11*t18*t19*t28-9.76264774128d-5
     1         *gammaaa*t14*t19*t24+5.228328d-2*gammaaa*t11*t19*t24)+1.3
     2         333333333333333d+0*t22*t51)-4.1839918891199995d-5*t16*t18
     3         *t32*t52-1.0d+0*t1*t10*t18*t51-1.86726d-3*t10*t2*t49-6.04
     4         3543839839999d-5*t16*t19*t47-3.73452d-3*t10*t3*t44+4.9793
     5         599999999993d-3*t10*t29*t41+7.46904d-3*t18*t19*t25*t30+3.
     6         73452d-3*t18*t19*t26*t29+1.3946639630399998d-5*t19*t23*t2
     7         8+6.2242d-3*t10*t25*t26)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = (1.0d+0*t1*t6*(-3.125040277951284d
     1         -7*t58*t7*t8/rhoa**1.2333333333333334d+1+1.673596755648d-
     2         4*gammaaa*t32*t50*t7-1.4938079999999998d-2*t19*t27*t7-3.1
     3         25040277951284d-7*gammaaa*t11*t4*t50*t58-3.12504027795128
     4         4d-7*t12*t13*t27*t58-7.812600694878209d-8*gammaaa*t14*t32
     5         *t50+4.1839918891199995d-5*gammaaa*t11*t32*t50-7.81260069
     6         4878209d-8*t15*t27*t32*t4+8.367983778239999d-5*t12*t27*t3
     7         2*t4+8.367983778239999d-5*t11*t27*t32*t4+1.39466396303999
     8         98d-5*t15*t19*t27+2.7893279260799997d-5*t14*t19*t27-1.493
     9         8079999999998d-2*t11*t19*t27)-7.812600694878209d-8*t16*t3
     :         2*t59+2.7893279260799997d-5*t19*t25*t52-3.73452d-3*t10*t2
     ;         6*t51+6.973319815199999d-6*t19*t29*t47-3.73452d-3*t10*t30
     <         *t46)*wght+Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_TA_TA_TA) = (1.0d+0*t1*t6*(-1.2500161111805136d-
     1         6*t58*t7*t8/rhoa**1.0333333333333333d+1-1.250016111180513
     2         6d-6*gammaaa*t11*t4*t58*t62-3.1250402779512837d-7*gammaaa
     3         *t14*t32*t62+5.020790266944d-4*gammaaa*t11*t32*t62-1.2500
     4         161111805136d-6*t12*t13*t58*t61-3.1250402779512837d-7*t15
     5         *t32*t4*t61+1.0041580533887999d-3*t12*t32*t4*t61+1.673596
     6         7556479997d-4*t15*t19*t61-1.7925696d-1*t12*t19*t61)+8.367
     7         983778239998d-5*t19*t29*t56-1.120356d-2*t10*t30*t55-3.125
     8         0402779512837d-7*t16*t17*t32)*wght+Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
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
            t32 = 1/rhoa**6.333333333333333d+0
            t33 = -3.3333333333333337d+0*t16*t2*t6*taua-6.66666666666666
     1         7d+0*t11*t13*t2*t4*taua-3.3333333333333337d+0*gammaaa*t10
     2         *t11*t12*taua-5.333333333333333d+0*t11*t32*t8*t9-2*t10*t3
     3         0*t31*t8*t9-2.6666666666666666d+0*gammaaa*t15*t29*t6-2*ga
     4         mmaaa*t12*t2*t30*t31*t4-t11*t16*t30*t4-2.6666666666666666
     5         d+0*gammaaa*t11*t12*t29*t4-2*t13*t14*t30*t31-gammaaa*t11*
     6         t15*t2*t30
            t34 = rhoa**3.333333333333333d-1
            t35 = 1.3333333333333333d+0*t17*t34+1.0d+0*t1*t33
            t36 = 1/rhob**3.6666666666666664d+0
            t37 = -6.2242d-3*t19*taub-4.9793599999999993d-3*gammabb*t36
            t38 = 1/t22**3
            t39 = 1/rhob**6.333333333333333d+0
            t40 = -3.3333333333333337d+0*gammabb*t12*t25*t26*taub-6.6666
     1         66666666667d+0*t13*t19*t21*t26*taub-3.3333333333333337d+0
     2         *t16*t19*t23*taub-5.333333333333333d+0*t24*t26*t39*t8-2*t
     3         24*t25*t37*t38*t8-2*t13*t27*t37*t38-2*gammabb*t12*t19*t21
     4         *t37*t38-t16*t21*t26*t37-gammabb*t15*t19*t26*t37-2.666666
     5         6666666666d+0*gammabb*t12*t21*t26*t36-2.6666666666666666d
     6         +0*gammabb*t15*t23*t36
            t41 = rhob**3.333333333333333d-1
            t42 = 1.3333333333333333d+0*t28*t41+1.0d+0*t18*t40
            t43 = 1/rhoa**8
            t44 = -3.73452d-3*t31*t43*t8*t9+2*gammaaa*t10*t11*t8+t15*t2*
     1         t6-3.73452d-3*gammaaa*t10*t12*t31*t4-1.86726d-3*t11*t16*t
     2         2*t4+t11*t12*t2*t4-3.73452d-3*t13*t14*t2*t31-1.86726d-3*g
     3         ammaaa*t10*t11*t15
            t45 = 1/t1
            t46 = 1/rhob**8
            t47 = -3.73452d-3*t24*t38*t46*t8+2*gammabb*t25*t26*t8-3.7345
     1         2d-3*t13*t19*t27*t38-3.73452d-3*gammabb*t12*t21*t25*t38-1
     2         .86726d-3*gammabb*t15*t25*t26-1.86726d-3*t16*t19*t21*t26+
     3         t12*t19*t21*t26+t15*t19*t23
            t48 = 1/t18
            t49 = 1/rhoa**7
            t50 = 1/rhoa**4.333333333333333d+0
            t51 = -7.46904d-3*t31*t49*t8*t9+2.0d+0*t16*t3*t6-7.46904d-3*
     1         gammaaa*t12*t31*t4*t50-3.73452d-3*gammaaa*t11*t15*t50+2.0
     2         d+0*gammaaa*t11*t12*t50-3.73452d-3*t11*t16*t3*t4+4.0d+0*t
     3         11*t13*t3*t4-7.46904d-3*t13*t14*t3*t31
            t52 = 1/t34
            t53 = 1/rhob**7
            t54 = 1/rhob**4.333333333333333d+0
            t55 = -7.46904d-3*t24*t38*t53*t8-7.46904d-3*gammabb*t12*t21*
     1         t38*t54-3.73452d-3*gammabb*t15*t26*t54+2.0d+0*gammabb*t12
     2         *t26*t54-7.46904d-3*t13*t20*t27*t38-3.73452d-3*t16*t20*t2
     3         1*t26+4.0d+0*t13*t20*t21*t26+2.0d+0*t16*t20*t23
            t56 = 1/t41
            t57 = t30**2
            t58 = 1/t5**4
            t59 = 1/rhoa**4.666666666666667d+0
            t60 = 1.6597866666666666d-2*t29*taua+1.8257653333333332d-2*g
     1         ammaaa*t59
            t61 = 1/rhoa**7.333333333333333d+0
            t62 = taua**2
            t63 = 8.88888888888889d+0*t16*t29*t6*taua+2.6666666666666666
     1         d+1*t13*t2*t30*t31*t4*taua+1.7777777777777778d+1*t11*t13*
     2         t29*t4*taua+2.6666666666666666d+1*gammaaa*t11*t12*t32*tau
     3         a+1.3333333333333333d+1*gammaaa*t10*t12*t30*t31*taua+6.66
     4         6666666666667d+0*t11*t16*t2*t30*taua+3.377777777777778d+1
     5         *t11*t61*t8*t9-2*t10*t31*t60*t8*t9+6*t10*t57*t58*t8*t9+2.
     6         1333333333333332d+1*t30*t31*t32*t8*t9+2.2222222222222224d
     7         +1*t10*t11*t13*t62-2*gammaaa*t12*t2*t31*t4*t60-t11*t16*t4
     8         *t60-2*t13*t14*t31*t60-gammaaa*t11*t15*t2*t60+9.777777777
     9         777779d+0*gammaaa*t15*t59*t6+9.777777777777779d+0*gammaaa
     :         *t11*t12*t4*t59+6*gammaaa*t12*t2*t4*t57*t58+6*t13*t14*t57
     ;         *t58+2*t16*t31*t4*t57+2*gammaaa*t15*t2*t31*t57+1.06666666
     <         66666666d+1*gammaaa*t12*t29*t30*t31*t4+5.333333333333333d
     =         +0*gammaaa*t11*t15*t29*t30
            t64 = 1/rhoa**6.666666666666666d-1
            t65 = 4.444444444444444d-1*t17*t64+1.0d+0*t1*t63+2.666666666
     1         6666666d+0*t33*t34
            t66 = -1.0d+0*t1*t17*t60-t30*t35-1.3333333333333333d+0*t17*t
     1         30*t34-1.0d+0*t1*t30*t33
            t67 = t37**2
            t68 = 1/t22**4
            t69 = 1/rhob**4.666666666666667d+0
            t70 = 1.6597866666666666d-2*t36*taub+1.8257653333333332d-2*g
     1         ammabb*t69
            t71 = 1/rhob**7.333333333333333d+0
            t72 = taub**2
            t73 = 2.6666666666666666d+1*gammabb*t12*t26*t39*taub+1.33333
     1         33333333333d+1*gammabb*t12*t25*t37*t38*taub+2.66666666666
     2         66666d+1*t13*t19*t21*t37*t38*taub+6.666666666666667d+0*t1
     3         6*t19*t26*t37*taub+1.7777777777777778d+1*t13*t21*t26*t36*
     4         taub+8.88888888888889d+0*t16*t23*t36*taub+3.3777777777777
     5         78d+1*t24*t26*t71*t8-2*t24*t25*t38*t70*t8+6*t24*t25*t67*t
     6         68*t8+2.1333333333333332d+1*t24*t37*t38*t39*t8+2.22222222
     7         22222224d+1*t13*t25*t26*t72-2*t13*t27*t38*t70-2*gammabb*t
     8         12*t19*t21*t38*t70-t16*t21*t26*t70-gammabb*t15*t19*t26*t7
     9         0+9.777777777777779d+0*gammabb*t12*t21*t26*t69+9.77777777
     :         7777779d+0*gammabb*t15*t23*t69+6*t13*t27*t67*t68+6*gammab
     ;         b*t12*t19*t21*t67*t68+2*t16*t21*t38*t67+2*gammabb*t15*t19
     <         *t38*t67+1.0666666666666666d+1*gammabb*t12*t21*t36*t37*t3
     =         8+5.333333333333333d+0*gammabb*t15*t26*t36*t37
            t74 = 1/rhob**6.666666666666666d-1
            t75 = 4.444444444444444d-1*t28*t74+1.0d+0*t18*t73+2.66666666
     1         66666666d+0*t40*t41
            t76 = -1.0d+0*t18*t28*t70-t37*t42-1.3333333333333333d+0*t28*
     1         t37*t41-1.0d+0*t18*t37*t40
            t77 = 1/rhoa**2.3333333333333334d+0
            t78 = 1/rhoa**9
            t79 = 1.24484d-2*gammaaa*t12*t31*t43*taua+2.48968d-2*t10*t13
     1         *t31*t4*taua+6.2242d-3*t10*t11*t16*taua-3.333333333333333
     2         7d+0*t10*t11*t12*taua+2.9876159999999996d-2*t31*t78*t8*t9
     3         +1.120356d-2*t30*t43*t58*t8*t9-1.0666666666666666d+1*gamm
     4         aaa*t11*t32*t8-4*gammaaa*t10*t30*t31*t8-2.666666666666666
     5         6d+0*t15*t29*t6+1.120356d-2*gammaaa*t10*t12*t30*t4*t58+1.
     6         120356d-2*t13*t14*t2*t30*t58+1.9917439999999997d-2*gammaa
     7         a*t12*t31*t32*t4+3.73452d-3*t16*t2*t30*t31*t4-2*t12*t2*t3
     8         0*t31*t4+4.9793599999999993d-3*t11*t16*t29*t4-2.666666666
     9         6666666d+0*t11*t12*t29*t4+9.958719999999999d-3*gammaaa*t1
     :         1*t15*t32+3.73452d-3*gammaaa*t10*t15*t30*t31+9.9587199999
     ;         99999d-3*t13*t14*t29*t31-t11*t15*t2*t30
            t80 = 1.0d+0*t1*t79+1.3333333333333333d+0*t34*t44
            t81 = 1/rhob**2.3333333333333334d+0
            t82 = 1/rhob**9
            t83 = 1.24484d-2*gammabb*t12*t38*t46*taub+2.48968d-2*t13*t21
     1         *t25*t38*taub+6.2242d-3*t16*t25*t26*taub-3.33333333333333
     2         37d+0*t12*t25*t26*taub+2.9876159999999996d-2*t24*t38*t8*t
     3         82+1.120356d-2*t24*t37*t46*t68*t8-1.0666666666666666d+1*g
     4         ammabb*t26*t39*t8-4*gammabb*t25*t37*t38*t8+1.120356d-2*t1
     5         3*t19*t27*t37*t68+1.120356d-2*gammabb*t12*t21*t25*t37*t68
     6         +1.9917439999999997d-2*gammabb*t12*t21*t38*t39+9.95871999
     7         9999999d-3*gammabb*t15*t26*t39+3.73452d-3*gammabb*t15*t25
     8         *t37*t38+3.73452d-3*t16*t19*t21*t37*t38-2*t12*t19*t21*t37
     9         *t38+9.958719999999999d-3*t13*t27*t36*t38-t15*t19*t26*t37
     :         +4.9793599999999993d-3*t16*t21*t26*t36-2.6666666666666666
     ;         d+0*t12*t21*t26*t36-2.6666666666666666d+0*t15*t23*t36
            t84 = 1.0d+0*t18*t83+1.3333333333333333d+0*t41*t47
            t85 = 1/rhoa**1.0666666666666666d+1
            t86 = 2.0919959445599998d-5*t58*t8*t85*t9-1.4938079999999998
     1         d-2*gammaaa*t31*t43*t8+2*t10*t11*t8+2.0919959445599998d-5
     2         *gammaaa*t12*t4*t43*t58+2.0919959445599998d-5*t10*t13*t14
     3         *t58+6.973319815199999d-6*gammaaa*t15*t31*t43+6.973319815
     4         199999d-6*t10*t16*t31*t4-7.46904d-3*t10*t12*t31*t4-3.7345
     5         2d-3*t10*t11*t15
            t87 = 1/rhoa**4
            t88 = 1/rhob**1.0666666666666666d+1
            t89 = 2.0919959445599998d-5*t24*t68*t8*t88-1.493807999999999
     1         8d-2*gammabb*t38*t46*t8+2*t25*t26*t8+2.0919959445599998d-
     2         5*gammabb*t12*t21*t46*t68+2.0919959445599998d-5*t13*t25*t
     3         27*t68+6.973319815199999d-6*gammabb*t15*t38*t46+6.9733198
     4         15199999d-6*t16*t21*t25*t38-7.46904d-3*t12*t21*t25*t38-3.
     5         73452d-3*t15*t25*t26
            t90 = 1/rhob**4
            t91 = 4.97936d-2*t13*t31*t4*t50*taua+1.24484d-2*t11*t16*t50*
     1         taua-1.3333333333333333d+1*t11*t13*t50*taua+2.48968d-2*ga
     2         mmaaa*t12*t31*t49*taua+2.240712d-2*t30*t49*t58*t8*t9+5.22
     3         8328d-2*t31*t43*t8*t9-3.3333333333333337d+0*t16*t2*t6+2.2
     4         40712d-2*gammaaa*t12*t30*t4*t50*t58+2.240712d-2*t13*t14*t
     5         3*t30*t58+7.46904d-3*gammaaa*t15*t30*t31*t50-4.0d+0*gamma
     6         aa*t12*t30*t31*t50+7.46904d-3*t16*t3*t30*t31*t4-8.0d+0*t1
     7         3*t3*t30*t31*t4+3.236584d-2*gammaaa*t10*t12*t31*t4+6.2242
     8         d-3*t11*t16*t2*t4-6.666666666666667d+0*t11*t13*t2*t4+1.24
     9         484d-2*t13*t14*t2*t31-2.0d+0*t11*t16*t3*t30+1.618292d-2*g
     :         ammaaa*t10*t11*t15-8.666666666666666d+0*gammaaa*t10*t11*t
     ;         12
            t92 = 1.0d+0*t1*t91+1.3333333333333333d+0*t34*t51
            t93 = 4.97936d-2*t13*t21*t38*t54*taub+1.24484d-2*t16*t26*t54
     1         *taub-1.3333333333333333d+1*t13*t26*t54*taub+2.48968d-2*g
     2         ammabb*t12*t38*t53*taub+2.240712d-2*t24*t37*t53*t68*t8+5.
     3         228328d-2*t24*t38*t46*t8+2.240712d-2*gammabb*t12*t21*t37*
     4         t54*t68+2.240712d-2*t13*t20*t27*t37*t68+7.46904d-3*gammab
     5         b*t15*t37*t38*t54-4.0d+0*gammabb*t12*t37*t38*t54+7.46904d
     6         -3*t16*t20*t21*t37*t38-8.0d+0*t13*t20*t21*t37*t38+1.24484
     7         d-2*t13*t19*t27*t38+3.236584d-2*gammabb*t12*t21*t25*t38-2
     8         .0d+0*t16*t20*t26*t37+1.618292d-2*gammabb*t15*t25*t26-8.6
     9         66666666666666d+0*gammabb*t12*t25*t26+6.2242d-3*t16*t19*t
     :         21*t26-6.666666666666667d+0*t13*t19*t21*t26-3.33333333333
     ;         33337d+0*t16*t19*t23
            t94 = 1.0d+0*t18*t93+1.3333333333333333d+0*t41*t55
            t95 = 1/rhoa**9.666666666666666d+0
            t96 = 4.1839918891199995d-5*t58*t8*t9*t95-1.4938079999999998
     1         d-2*gammaaa*t31*t49*t8+4.1839918891199995d-5*t13*t14*t50*
     2         t58+4.1839918891199995d-5*gammaaa*t12*t4*t49*t58+1.394663
     3         9630399998d-5*t16*t31*t4*t50-1.4938079999999998d-2*t13*t3
     4         1*t4*t50-7.46904d-3*t12*t31*t4*t50-3.73452d-3*t11*t16*t50
     5         -3.73452d-3*t11*t15*t50+2.0d+0*t11*t12*t50+1.394663963039
     6         9998d-5*gammaaa*t15*t31*t49-7.46904d-3*gammaaa*t12*t31*t4
     7         9
            t97 = 1/rhoa**3
            t98 = 1/rhob**9.666666666666666d+0
            t99 = 4.1839918891199995d-5*t24*t68*t8*t98-1.493807999999999
     1         8d-2*gammabb*t38*t53*t8+4.1839918891199995d-5*t13*t27*t54
     2         *t68+4.1839918891199995d-5*gammabb*t12*t21*t53*t68+1.3946
     3         639630399998d-5*t16*t21*t38*t54-1.4938079999999998d-2*t13
     4         *t21*t38*t54-7.46904d-3*t12*t21*t38*t54-3.73452d-3*t16*t2
     5         6*t54-3.73452d-3*t15*t26*t54+2.0d+0*t12*t26*t54+1.3946639
     6         630399998d-5*gammabb*t15*t38*t53-7.46904d-3*gammabb*t12*t
     7         38*t53
            t100 = 1/rhob**3
            t101 = 1/rhoa**8.666666666666666d+0
            t102 = 1/rhoa**6
            t103 = 1/rhoa**3.3333333333333337d+0
            t104 = 8.367983778239999d-5*t101*t58*t8*t9+8.367983778239999
     1         d-5*gammaaa*t102*t12*t4*t58+8.367983778239999d-5*t103*t13
     2         *t14*t58+2.7893279260799997d-5*t103*t16*t31*t4-5.97523199
     3         9999999d-2*t103*t13*t31*t4+2.7893279260799997d-5*gammaaa*
     4         t102*t15*t31-2.9876159999999996d-2*gammaaa*t102*t12*t31-1
     5         .4938079999999998d-2*t103*t11*t16+8.0d+0*t103*t11*t13
            t105 = 1/rhoa**2
            t106 = 1/rhob**8.666666666666666d+0
            t107 = 1/rhob**6
            t108 = 1/rhob**3.3333333333333337d+0
            t109 = 8.367983778239999d-5*t106*t24*t68*t8+8.36798377823999
     1         9d-5*t108*t13*t27*t68+8.367983778239999d-5*gammabb*t107*t
     2         12*t21*t68+2.7893279260799997d-5*t108*t16*t21*t38-5.97523
     3         1999999999d-2*t108*t13*t21*t38+2.7893279260799997d-5*gamm
     4         abb*t107*t15*t38-2.9876159999999996d-2*gammabb*t107*t12*t
     5         38-1.4938079999999998d-2*t108*t16*t26+8.0d+0*t108*t13*t26
            t110 = 1/rhob**2
            t111 = t30**3
            t112 = 1/t5**5
            t113 = 1/rhoa**5.666666666666667d+0
            t114 = -6.085884444444444d-2*t59*taua-8.520238222222222d-2*g
     1         ammaaa*t113
            t115 = t37**3
            t116 = 1/t22**5
            t117 = 1/rhob**5.666666666666667d+0
            t118 = -6.085884444444444d-2*t69*taub-8.520238222222222d-2*g
     1         ammabb*t117
            t119 = 1/rhoa**5
            t120 = 1/rhob**5
            t121 = 1/rhoa**7.666666666666667d+0
            t122 = 1/rhob**7.666666666666667d+0
            fnc(iq) = (1.0d+0*t1*t17*t6+1.0d+0*t18*t23*t28)*wght+fnc(iq)
            Amat(iq,D1_RA) = (t35*t6-1.0d+0*t1*t11*t17*t30)*wght+Amat(iq
     1         ,D1_RA)
            Amat(iq,D1_RB) = (t23*t42-1.0d+0*t18*t26*t28*t37)*wght+Amat(
     1         iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t1*t44*t6-1.86726d-3*t11*t17*t45)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t18*t23*t47-1.86726d-3*t26*t28*t48
     1         )*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TA) = (1.0d+0*t1*t51*t6-3.73452d-3*t11*t17*t52)*w
     1         ght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = (1.0d+0*t18*t23*t55-3.73452d-3*t26*t28*t56)
     1         *wght+Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = (t11*t66+t6*t65+2.0d+0*t1*t17*t31*t57)*
     1         wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t26*t76+t23*t75+2.0d+0*t18*t28*t38*t67
     1         )*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t6*t80+4.9793599999999993d-3*t11*t17*
     1         t77+3.73452d-3*t17*t30*t31*t45-1.0d+0*t1*t11*t30*t44-1.86
     2         726d-3*t11*t2*t35)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t23*t84+4.9793599999999993d-3*t26*t28
     1         *t81+3.73452d-3*t28*t37*t38*t48-1.0d+0*t18*t26*t37*t47-1.
     2         86726d-3*t19*t26*t42)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (6.973319815199999d-6*t17*t31*t87+1.0
     1         d+0*t1*t6*t86-3.73452d-3*t11*t44*t45)*wght+Cmat2(iq,D2_GA
     2         A_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (6.973319815199999d-6*t28*t38*t90+1.0
     1         d+0*t18*t23*t89-3.73452d-3*t26*t47*t48)*wght+Cmat2(iq,D2_
     2         GBB_GBB)
            Mmat2(iq,D2_RA_TA) = (t6*t92+7.46904d-3*t17*t30*t31*t52-1.0d
     1         +0*t1*t11*t30*t51+6.2242d-3*t11*t17*t45-3.73452d-3*t11*t3
     2         *t35)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (t23*t94+7.46904d-3*t28*t37*t38*t56-1.0
     1         d+0*t18*t26*t37*t55+6.2242d-3*t26*t28*t48-3.73452d-3*t20*
     2         t26*t42)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = (1.3946639630399998d-5*t17*t31*t97+1.0
     1         d+0*t1*t6*t96-3.73452d-3*t11*t44*t52-1.86726d-3*t11*t45*t
     2         51)*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = (1.0d+0*t18*t23*t99-3.73452d-3*t26*t47
     1         *t56-1.86726d-3*t26*t48*t55+1.3946639630399998d-5*t100*t2
     2         8*t38)*wght+Mmat2(iq,D2_GBB_TB)
            Mmat2(iq,D2_TA_TA) = (1.0d+0*t1*t104*t6-7.46904d-3*t11*t51*t
     1         52+2.7893279260799997d-5*t105*t17*t31)*wght+Mmat2(iq,D2_T
     2         A_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = (-7.46904d-3*t26*t55*t56+2.789327926079
     1         9997d-5*t110*t28*t38+1.0d+0*t109*t18*t23)*wght+Mmat2(iq,D
     2         2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = (t6*(1.0d+0*t1*(-2.0148148148148146d
     1         +2*gammaaa*t11*t12*t61*taua+4.0d+1*t13*t2*t31*t4*t60*taua
     2         +2.0d+1*gammaaa*t10*t12*t31*t60*taua+1.0d+1*t11*t16*t2*t6
     3         0*taua-3.2592592592592595d+1*t16*t59*t6*taua-6.5185185185
     4         18519d+1*t11*t13*t4*t59*taua-1.2d+2*t13*t2*t4*t57*t58*tau
     5         a-6.0d+1*gammaaa*t10*t12*t57*t58*taua-2.0d+1*t16*t2*t31*t
     6         57*taua-1.0666666666666667d+2*t13*t29*t30*t31*t4*taua-1.6
     7         d+2*gammaaa*t12*t30*t31*t32*taua-2.6666666666666666d+1*t1
     8         1*t16*t29*t30*taua-2.0266666666666666d+2*t30*t31*t61*t8*t
     9         9+18*t10*t30*t58*t60*t8*t9+32*t31*t32*t60*t8*t9-96*t32*t5
     :         7*t58*t8*t9-2*t10*t114*t31*t8*t9-24*t10*t111*t112*t8*t9-2
     ;         .477037037037037d+2*t11*t8*t9/rhoa**8.333333333333334d+0-
     <         1.777777777777778d+2*t11*t13*t32*t62-1.3333333333333333d+
     =         2*t10*t13*t30*t31*t62+18*gammaaa*t12*t2*t30*t4*t58*t60+18
     >         *t13*t14*t30*t58*t60+6*t16*t30*t31*t4*t60+16*gammaaa*t12*
     ?         t29*t31*t4*t60+6*gammaaa*t15*t2*t30*t31*t60+8*gammaaa*t11
     @         *t15*t29*t60-4.5629629629629626d+1*gammaaa*t113*t15*t6-5.
     1         866666666666667d+1*gammaaa*t12*t30*t31*t4*t59-2.933333333
     2         3333333d+1*gammaaa*t11*t15*t30*t59-48*gammaaa*t12*t29*t4*
     3         t57*t58-6*t111*t16*t4*t58-6*gammaaa*t111*t15*t2*t58-16*ga
     4         mmaaa*t15*t29*t31*t57-2*gammaaa*t114*t12*t2*t31*t4-24*gam
     5         maaa*t111*t112*t12*t2*t4-t11*t114*t16*t4-4.56296296296296
     6         26d+1*gammaaa*t11*t113*t12*t4-2*t114*t13*t14*t31-gammaaa*
     7         t11*t114*t15*t2-24*t111*t112*t13*t14)+1.3333333333333333d
     8         +0*t33*t64+4.0d+0*t34*t63-2.962962962962963d-1*t17*t3)+t3
     9         1*(-2*t30*t66+4.0d+0*t1*t17*t30*t60+2.6666666666666666d+0
     :         *t17*t34*t57+2.0d+0*t1*t33*t57)+t11*(-2*t30*t65-4.4444444
     ;         44444444d-1*t17*t30*t64-1.0d+0*t1*t30*t63-t35*t60-2.66666
     <         66666666666d+0*t17*t34*t60-2.0d+0*t1*t33*t60-2.6666666666
     =         666666d+0*t30*t33*t34-1.0d+0*t1*t114*t17)-6.0d+0*t1*t111*
     >         t17*t58)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t23*(1.0d+0*t18*(-2.014814814814814
     1         6d+2*gammabb*t12*t26*t71*taub+2.0d+1*gammabb*t12*t25*t38*
     2         t70*taub+4.0d+1*t13*t19*t21*t38*t70*taub+1.0d+1*t16*t19*t
     3         26*t70*taub-6.518518518518519d+1*t13*t21*t26*t69*taub-3.2
     4         592592592592595d+1*t16*t23*t69*taub-6.0d+1*gammabb*t12*t2
     5         5*t67*t68*taub-1.2d+2*t13*t19*t21*t67*t68*taub-2.0d+1*t16
     6         *t19*t38*t67*taub-1.6d+2*gammabb*t12*t37*t38*t39*taub-1.0
     7         666666666666667d+2*t13*t21*t36*t37*t38*taub-2.66666666666
     8         66666d+1*t16*t26*t36*t37*taub-2.0266666666666666d+2*t24*t
     9         37*t38*t71*t8+18*t24*t25*t37*t68*t70*t8+32*t24*t38*t39*t7
     :         0*t8-96*t24*t39*t67*t68*t8-2*t118*t24*t25*t38*t8-2.477037
     ;         037037037d+2*t24*t26*t8/rhob**8.333333333333334d+0-24*t11
     <         5*t116*t24*t25*t8-1.777777777777778d+2*t13*t26*t39*t72-1.
     =         3333333333333333d+2*t13*t25*t37*t38*t72+18*t13*t27*t37*t6
     >         8*t70+18*gammabb*t12*t19*t21*t37*t68*t70+6*t16*t21*t37*t3
     ?         8*t70+6*gammabb*t15*t19*t37*t38*t70+16*gammabb*t12*t21*t3
     @         6*t38*t70+8*gammabb*t15*t26*t36*t70-5.866666666666667d+1*
     1         gammabb*t12*t21*t37*t38*t69-2.9333333333333333d+1*gammabb
     2         *t15*t26*t37*t69-48*gammabb*t12*t21*t36*t67*t68-6*t115*t1
     3         6*t21*t68-6*gammabb*t115*t15*t19*t68-16*gammabb*t15*t36*t
     4         38*t67-2*t118*t13*t27*t38-2*gammabb*t118*t12*t19*t21*t38-
     5         24*t115*t116*t13*t27-t118*t16*t21*t26-4.5629629629629626d
     6         +1*gammabb*t117*t12*t21*t26-gammabb*t118*t15*t19*t26-4.56
     7         29629629629626d+1*gammabb*t117*t15*t23-24*gammabb*t115*t1
     8         16*t12*t19*t21)+1.3333333333333333d+0*t40*t74+4.0d+0*t41*
     9         t73-2.962962962962963d-1*t20*t28)+t38*(-2*t37*t76+4.0d+0*
     :         t18*t28*t37*t70+2.6666666666666666d+0*t28*t41*t67+2.0d+0*
     ;         t18*t40*t67)+t26*(-2*t37*t75-4.444444444444444d-1*t28*t37
     <         *t74-1.0d+0*t18*t37*t73-t42*t70-2.6666666666666666d+0*t28
     =         *t41*t70-2.0d+0*t18*t40*t70-2.6666666666666666d+0*t37*t40
     >         *t41-1.0d+0*t118*t18*t28)-6.0d+0*t115*t18*t28*t68)*wght+A
     ?         mat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t6*(1.0d+0*t1*(-1.6597866666666666
     1         d-1*gammaaa*t12*t31*t78*taua-7.46904d-2*gammaaa*t12*t30*t
     2         43*t58*taua-1.493808d-1*t10*t13*t30*t4*t58*taua-1.9917439
     3         999999997d-1*t13*t31*t32*t4*taua-4.9793599999999993d-2*t1
     4         1*t16*t32*taua+2.6666666666666666d+1*t11*t12*t32*taua-2.4
     5         8968d-2*t10*t16*t30*t31*taua+1.3333333333333333d+1*t10*t1
     6         2*t30*t31*taua-1.7925696d-1*t30*t58*t78*t8*t9+1.120356d-2
     7         *t43*t58*t60*t8*t9-4.481424d-2*t112*t43*t57*t8*t9-2.68885
     8         44000000003d-1*t31*t8*t9/rhoa**10+6.755555555555556d+1*ga
     9         mmaaa*t11*t61*t8-4*gammaaa*t10*t31*t60*t8+12*gammaaa*t10*
     :         t57*t58*t8+4.2666666666666664d+1*gammaaa*t30*t31*t32*t8-8
     ;         .298933333333335d-2*t13*t31*t43*t62-1.2614378666666665d-1
     <         *gammaaa*t12*t31*t4*t61-6.307189333333332d-2*gammaaa*t11*
     =         t15*t61+1.120356d-2*gammaaa*t10*t12*t4*t58*t60+1.120356d-
     >         2*t13*t14*t2*t58*t60+3.73452d-3*t16*t2*t31*t4*t60-2*t12*t
     ?         2*t31*t4*t60+3.73452d-3*gammaaa*t10*t15*t31*t60-t11*t15*t
     @         2*t60+9.777777777777779d+0*t15*t59*t6-1.8257653333333332d
     1         -2*t11*t16*t4*t59+9.777777777777779d+0*t11*t12*t4*t59-3.6
     2         515306666666664d-2*t13*t14*t31*t59-1.120356d-2*t16*t2*t4*
     3         t57*t58+6*t12*t2*t4*t57*t58-1.120356d-2*gammaaa*t10*t15*t
     4         57*t58-1.1950464d-1*gammaaa*t12*t30*t32*t4*t58-5.97523199
     5         9999999d-2*t13*t14*t29*t30*t58-4.481424d-2*gammaaa*t10*t1
     6         12*t12*t4*t57+2*t15*t2*t31*t57-4.481424d-2*t112*t13*t14*t
     7         2*t57-1.9917439999999997d-2*t16*t29*t30*t31*t4+1.06666666
     8         66666666d+1*t12*t29*t30*t31*t4-3.9834879999999995d-2*gamm
     9         aaa*t15*t30*t31*t32+5.333333333333333d+0*t11*t15*t29*t30)
     :         +2.6666666666666666d+0*t34*t79+4.444444444444444d-1*t44*t
     ;         64)+t11*(-t30*t80-1.0d+0*t1*t30*t79+4.9793599999999993d-3
     <         *t33*t77-1.0d+0*t1*t44*t60-1.3333333333333333d+0*t30*t34*
     =         t44+4.9793599999999993d-3*t29*t35-1.1618506666666667d-2*t
     >         103*t17)-1.9917439999999997d-2*t17*t30*t31*t77-3.73452d-3
     ?         *t2*t31*t66-1.86726d-3*t11*t2*t65-1.120356d-2*t17*t45*t57
     @         *t58+2.0d+0*t1*t31*t44*t57)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t23*(1.0d+0*t18*(-1.65978666666666
     1         66d-1*gammabb*t12*t38*t82*taub-7.46904d-2*gammabb*t12*t37
     2         *t46*t68*taub-1.493808d-1*t13*t21*t25*t37*t68*taub-1.9917
     3         439999999997d-1*t13*t21*t38*t39*taub-4.9793599999999993d-
     4         2*t16*t26*t39*taub+2.6666666666666666d+1*t12*t26*t39*taub
     5         -2.48968d-2*t16*t25*t37*t38*taub+1.3333333333333333d+1*t1
     6         2*t25*t37*t38*taub-1.7925696d-1*t24*t37*t68*t8*t82+6.7555
     7         55555555556d+1*gammabb*t26*t71*t8+1.120356d-2*t24*t46*t68
     8         *t70*t8-4*gammabb*t25*t38*t70*t8+12*gammabb*t25*t67*t68*t
     9         8-4.481424d-2*t116*t24*t46*t67*t8+4.2666666666666664d+1*g
     :         ammabb*t37*t38*t39*t8-2.6888544000000003d-1*t24*t38*t8/rh
     ;         ob**10-8.298933333333335d-2*t13*t38*t46*t72-1.26143786666
     <         66665d-1*gammabb*t12*t21*t38*t71-6.307189333333332d-2*gam
     =         mabb*t15*t26*t71+1.120356d-2*t13*t19*t27*t68*t70+1.120356
     >         d-2*gammabb*t12*t21*t25*t68*t70+3.73452d-3*gammabb*t15*t2
     ?         5*t38*t70+3.73452d-3*t16*t19*t21*t38*t70-2*t12*t19*t21*t3
     @         8*t70-t15*t19*t26*t70-3.6515306666666664d-2*t13*t27*t38*t
     1         69-1.8257653333333332d-2*t16*t21*t26*t69+9.77777777777777
     2         9d+0*t12*t21*t26*t69+9.777777777777779d+0*t15*t23*t69-1.1
     3         20356d-2*gammabb*t15*t25*t67*t68-1.120356d-2*t16*t19*t21*
     4         t67*t68+6*t12*t19*t21*t67*t68-1.1950464d-1*gammabb*t12*t2
     5         1*t37*t39*t68-5.975231999999999d-2*t13*t27*t36*t37*t68+2*
     6         t15*t19*t38*t67-4.481424d-2*t116*t13*t19*t27*t67-4.481424
     7         d-2*gammabb*t116*t12*t21*t25*t67-3.9834879999999995d-2*ga
     8         mmabb*t15*t37*t38*t39-1.9917439999999997d-2*t16*t21*t36*t
     9         37*t38+1.0666666666666666d+1*t12*t21*t36*t37*t38+5.333333
     :         333333333d+0*t15*t26*t36*t37)+2.6666666666666666d+0*t41*t
     ;         83+4.444444444444444d-1*t47*t74)+t26*(-t37*t84-1.0d+0*t18
     <         *t37*t83+4.9793599999999993d-3*t40*t81-1.0d+0*t18*t47*t70
     =         -1.3333333333333333d+0*t37*t41*t47+4.9793599999999993d-3*
     >         t36*t42-1.1618506666666667d-2*t108*t28)-1.991743999999999
     ?         7d-2*t28*t37*t38*t81-3.73452d-3*t19*t38*t76-1.86726d-3*t1
     @         9*t26*t75-1.120356d-2*t28*t48*t67*t68+2.0d+0*t18*t38*t47*
     1         t67)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t6*(1.0d+0*t1*(-6.9733198152d-5*g
     1         ammaaa*t12*t58*t85*taua-1.3946639630399998d-4*t13*t4*t43*
     2         t58*taua-2.3244399384d-5*t16*t31*t43*taua+2.48968d-2*t12*
     3         t31*t43*taua-8.367983778239999d-5*t112*t30*t8*t85*t9-2.23
     4         14623408639997d-4*t58*t8*t9/rhoa**1.1666666666666665d+1+1
     5         .1950464d-1*gammaaa*t31*t78*t8+4.481424d-2*gammaaa*t30*t4
     6         3*t58*t8-1.0666666666666666d+1*t11*t32*t8-4*t10*t30*t31*t
     7         8-1.6735967556479997d-4*gammaaa*t12*t4*t58*t78-5.57865585
     8         2159999d-5*gammaaa*t15*t31*t78-2.0919959445599998d-5*gamm
     9         aaa*t15*t30*t43*t58-2.0919959445599998d-5*t10*t16*t30*t4*
     :         t58+2.240712d-2*t10*t12*t30*t4*t58-1.1157311704319998d-4*
     ;         t13*t14*t32*t58-8.367983778239999d-5*gammaaa*t112*t12*t30
     <         *t4*t43-3.7191039014399996d-5*t16*t31*t32*t4+3.9834879999
     =         999995d-2*t12*t31*t32*t4+1.9917439999999997d-2*t11*t15*t3
     >         2+7.46904d-3*t10*t15*t30*t31-8.367983778239999d-5*t10*t11
     ?         2*t13*t14*t30)+1.3333333333333333d+0*t34*t86)-2.091995944
     @         5599998d-5*t17*t30*t58*t87-1.0d+0*t1*t11*t30*t86-3.73452d
     1         -3*t11*t2*t80+9.958719999999999d-3*t11*t44*t77+7.46904d-3
     2         *t30*t31*t44*t45+6.973319815199999d-6*t10*t31*t35-3.71910
     3         39014399996d-5*t119*t17*t31)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t23*(1.0d+0*t18*(-6.9733198152d-5
     1         *gammabb*t12*t68*t88*taub-1.3946639630399998d-4*t13*t21*t
     2         46*t68*taub-2.3244399384d-5*t16*t38*t46*taub+2.48968d-2*t
     3         12*t38*t46*taub-8.367983778239999d-5*t116*t24*t37*t8*t88+
     4         1.1950464d-1*gammabb*t38*t8*t82-1.6735967556479997d-4*gam
     5         mabb*t12*t21*t68*t82-5.578655852159999d-5*gammabb*t15*t38
     6         *t82+4.481424d-2*gammabb*t37*t46*t68*t8-2.231462340863999
     7         7d-4*t24*t68*t8/rhob**1.1666666666666665d+1-1.06666666666
     8         66666d+1*t26*t39*t8-4*t25*t37*t38*t8-2.0919959445599998d-
     9         5*gammabb*t15*t37*t46*t68-1.1157311704319998d-4*t13*t27*t
     :         39*t68-2.0919959445599998d-5*t16*t21*t25*t37*t68+2.240712
     ;         d-2*t12*t21*t25*t37*t68-8.367983778239999d-5*gammabb*t116
     <         *t12*t21*t37*t46-3.7191039014399996d-5*t16*t21*t38*t39+3.
     =         9834879999999995d-2*t12*t21*t38*t39+1.9917439999999997d-2
     >         *t15*t26*t39+7.46904d-3*t15*t25*t37*t38-8.367983778239999
     ?         d-5*t116*t13*t25*t27*t37)+1.3333333333333333d+0*t41*t89)-
     @         2.0919959445599998d-5*t28*t37*t68*t90-1.0d+0*t18*t26*t37*
     1         t89-3.73452d-3*t19*t26*t84+9.958719999999999d-3*t26*t47*t
     2         81+7.46904d-3*t37*t38*t47*t48+6.973319815199999d-6*t25*t3
     3         8*t42-3.7191039014399996d-5*t120*t28*t38)*wght+Cmat3(iq,D
     4         3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (1.0d+0*t1*t6*(-1.562520138975642
     1         d-7*t112*t8*t9/rhoa**1.3333333333333333d+1+1.255197566736
     2         d-4*gammaaa*t58*t8*t85-3.9063003474391045d-8*gammaaa*t15*
     3         t58*t85-1.562520138975642d-7*gammaaa*t112*t12*t4*t85-2.24
     4         0712d-2*t31*t43*t8-3.9063003474391045d-8*t16*t4*t43*t58+6
     5         .275987833679999d-5*t12*t4*t43*t58+2.0919959445599995d-5*
     6         t15*t31*t43-1.562520138975642d-7*t112*t13*t14*t43)+2.0919
     7         959445599995d-5*t31*t44*t87-5.60178d-3*t11*t45*t86-3.9063
     8         003474391045d-8*t17*t58/rhoa**6.666666666666667d+0)*wght+
     9         Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (2.0919959445599995d-5*t38*t47*t9
     1         0-5.60178d-3*t26*t48*t89+1.0d+0*t18*t23*(1.255197566736d-
     2         4*gammabb*t68*t8*t88-3.9063003474391045d-8*gammabb*t15*t6
     3         8*t88-1.562520138975642d-7*gammabb*t116*t12*t21*t88-2.240
     4         712d-2*t38*t46*t8-1.562520138975642d-7*t116*t24*t8/rhob**
     5         1.3333333333333333d+1-3.9063003474391045d-8*t16*t21*t46*t
     6         68+6.275987833679999d-5*t12*t21*t46*t68+2.091995944559999
     7         5d-5*t15*t38*t46-1.562520138975642d-7*t116*t13*t27*t46)-3
     8         .9063003474391045d-8*t28*t68/rhob**6.666666666666667d+0)*
     9         wght+Cmat3(iq,D3_GBB_GBB_GBB)
            Mmat3(iq,D3_RA_RA_TA) = (t6*(1.0d+0*t1*(-2.987616d-1*t13*t30
     1         *t4*t50*t58*taua-1.493808d-1*gammaaa*t12*t30*t49*t58*taua
     2         -4.97936d-2*t16*t30*t31*t50*taua+5.333333333333333d+1*t13
     3         *t30*t31*t50*taua-2.8216373333333333d-1*gammaaa*t12*t31*t
     4         43*taua-2.9876159999999996d-1*t10*t13*t31*t4*taua-7.46903
     5         9999999999d-2*t10*t11*t16*taua+8.0d+1*t10*t11*t13*taua-4.
     6         1826624d-1*t31*t78*t8*t9+2.240712d-2*t49*t58*t60*t8*t9-3.
     7         1369968d-1*t30*t43*t58*t8*t9-8.962848d-2*t112*t49*t57*t8*
     8         t9-1.659786666666667d-1*t13*t31*t49*t62+2.240712d-2*gamma
     9         aa*t12*t4*t50*t58*t60+2.240712d-2*t13*t14*t3*t58*t60+7.46
     :         904d-3*gammaaa*t15*t31*t50*t60-4.0d+0*gammaaa*t12*t31*t50
     ;         *t60+7.46904d-3*t16*t3*t31*t4*t60-8.0d+0*t13*t3*t31*t4*t6
     <         0-2.0d+0*t11*t16*t3*t60+8.88888888888889d+0*t16*t29*t6-2.
     =         240712d-2*gammaaa*t15*t50*t57*t58+1.2d+1*gammaaa*t12*t50*
     >         t57*t58-2.240712d-2*t16*t3*t4*t57*t58+2.4d+1*t13*t3*t4*t5
     ?         7*t58-1.9419503999999999d-1*gammaaa*t10*t12*t30*t4*t58-7.
     @         46904d-2*t13*t14*t2*t30*t58-8.962848d-2*gammaaa*t112*t12*
     1         t4*t50*t57+4.0d+0*t16*t3*t31*t57-8.962848d-2*t112*t13*t14
     2         *t3*t57-1.7261781333333331d-1*gammaaa*t12*t31*t32*t4-2.48
     3         968d-2*t16*t2*t30*t31*t4+2.6666666666666666d+1*t13*t2*t30
     4         *t31*t4-1.6597866666666666d-2*t11*t16*t29*t4+1.7777777777
     5         777778d+1*t11*t13*t29*t4-8.630890666666666d-2*gammaaa*t11
     6         *t15*t32+4.622222222222223d+1*gammaaa*t11*t12*t32-6.47316
     7         8d-2*gammaaa*t10*t15*t30*t31+3.466666666666666d+1*gammaaa
     8         *t10*t12*t30*t31-3.319573333333333d-2*t13*t14*t29*t31+6.6
     9         66666666666667d+0*t11*t16*t2*t30)+2.6666666666666666d+0*t
     :         34*t91+4.444444444444444d-1*t51*t64)+t11*(-t30*t92-1.0d+0
     ;         *t1*t30*t91-8.298933333333333d-3*t17*t77-1.0d+0*t1*t51*t6
     <         0-1.3333333333333333d+0*t30*t34*t51+6.2242d-3*t33*t45+6.2
     =         242d-3*t2*t35)-7.46904d-3*t3*t31*t66-3.73452d-3*t11*t3*t6
     >         5-2.240712d-2*t17*t52*t57*t58+2.0d+0*t1*t31*t51*t57-2.489
     ?         68d-2*t17*t30*t31*t45)*wght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (t23*(1.0d+0*t18*(-2.987616d-1*t13*t
     1         21*t37*t54*t68*taub-1.493808d-1*gammabb*t12*t37*t53*t68*t
     2         aub-4.97936d-2*t16*t37*t38*t54*taub+5.333333333333333d+1*
     3         t13*t37*t38*t54*taub-2.8216373333333333d-1*gammabb*t12*t3
     4         8*t46*taub-2.9876159999999996d-1*t13*t21*t25*t38*taub-7.4
     5         69039999999999d-2*t16*t25*t26*taub+8.0d+1*t13*t25*t26*tau
     6         b-4.1826624d-1*t24*t38*t8*t82+2.240712d-2*t24*t53*t68*t70
     7         *t8-3.1369968d-1*t24*t37*t46*t68*t8-8.962848d-2*t116*t24*
     8         t53*t67*t8-1.659786666666667d-1*t13*t38*t53*t72+2.240712d
     9         -2*gammabb*t12*t21*t54*t68*t70+2.240712d-2*t13*t20*t27*t6
     :         8*t70+7.46904d-3*gammabb*t15*t38*t54*t70-4.0d+0*gammabb*t
     ;         12*t38*t54*t70+7.46904d-3*t16*t20*t21*t38*t70-8.0d+0*t13*
     <         t20*t21*t38*t70-2.0d+0*t16*t20*t26*t70-2.240712d-2*gammab
     =         b*t15*t54*t67*t68+1.2d+1*gammabb*t12*t54*t67*t68-2.240712
     >         d-2*t16*t20*t21*t67*t68+2.4d+1*t13*t20*t21*t67*t68-7.4690
     ?         4d-2*t13*t19*t27*t37*t68-1.9419503999999999d-1*gammabb*t1
     @         2*t21*t25*t37*t68-8.962848d-2*gammabb*t116*t12*t21*t54*t6
     1         7+4.0d+0*t16*t20*t38*t67-8.962848d-2*t116*t13*t20*t27*t67
     2         -1.7261781333333331d-1*gammabb*t12*t21*t38*t39-8.63089066
     3         6666666d-2*gammabb*t15*t26*t39+4.622222222222223d+1*gamma
     4         bb*t12*t26*t39-6.473168d-2*gammabb*t15*t25*t37*t38+3.4666
     5         66666666666d+1*gammabb*t12*t25*t37*t38-2.48968d-2*t16*t19
     6         *t21*t37*t38+2.6666666666666666d+1*t13*t19*t21*t37*t38-3.
     7         319573333333333d-2*t13*t27*t36*t38+6.666666666666667d+0*t
     8         16*t19*t26*t37-1.6597866666666666d-2*t16*t21*t26*t36+1.77
     9         77777777777778d+1*t13*t21*t26*t36+8.88888888888889d+0*t16
     :         *t23*t36)+2.6666666666666666d+0*t41*t93+4.444444444444444
     ;         d-1*t55*t74)+t26*(-t37*t94-1.0d+0*t18*t37*t93-8.298933333
     <         333333d-3*t28*t81-1.0d+0*t18*t55*t70-1.3333333333333333d+
     =         0*t37*t41*t55+6.2242d-3*t40*t48+6.2242d-3*t19*t42)-7.4690
     >         4d-3*t20*t38*t76-3.73452d-3*t20*t26*t75-2.240712d-2*t28*t
     ?         56*t67*t68+2.0d+0*t18*t38*t55*t67-2.48968d-2*t28*t37*t38*
     @         t48)*wght+Mmat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (t6*(1.0d+0*t1*(-1.3946639630399998
     1         d-4*gammaaa*t12*t58*t95*taua-2.7893279260799997d-4*t13*t4
     2         *t49*t58*taua-4.6488798768d-5*t16*t31*t49*taua+4.97936d-2
     3         *t13*t31*t49*taua+2.48968d-2*t12*t31*t49*taua-1.673596755
     4         648d-4*t112*t30*t8*t9*t95-4.044525492816d-4*t58*t8*t85*t9
     5         +4.481424d-2*gammaaa*t30*t49*t58*t8+1.0456655999999999d-1
     6         *gammaaa*t31*t43*t8-4.1839918891199995d-5*t16*t30*t4*t50*
     7         t58+4.481424d-2*t13*t30*t4*t50*t58+2.240712d-2*t12*t30*t4
     8         *t50*t58-4.1839918891199995d-5*gammaaa*t15*t30*t49*t58+2.
     9         240712d-2*gammaaa*t12*t30*t49*t58-2.9287943223839996d-4*g
     :         ammaaa*t12*t4*t43*t58-1.8130631519519996d-4*t10*t13*t14*t
     ;         58+7.46904d-3*t16*t30*t31*t50+7.46904d-3*t15*t30*t31*t50-
     <         4.0d+0*t12*t30*t31*t50-1.673596755648d-4*t112*t13*t14*t30
     =         *t50-1.673596755648d-4*gammaaa*t112*t12*t30*t4*t49-9.7626
     >         4774128d-5*gammaaa*t15*t31*t43+5.228328d-2*gammaaa*t12*t3
     ?         1*t43-6.043543839839999d-5*t10*t16*t31*t4+6.473168d-2*t10
     @         *t13*t31*t4+3.236584d-2*t10*t12*t31*t4+1.618292d-2*t10*t1
     1         1*t16+1.618292d-2*t10*t11*t15-8.666666666666666d+0*t10*t1
     2         1*t12)+1.3333333333333333d+0*t34*t96)-4.1839918891199995d
     3         -5*t17*t30*t58*t97-1.0d+0*t1*t11*t30*t96-1.86726d-3*t11*t
     4         2*t92-6.043543839839999d-5*t17*t31*t87-3.73452d-3*t11*t3*
     5         t80+4.9793599999999993d-3*t11*t51*t77+7.46904d-3*t30*t31*
     6         t44*t52+3.73452d-3*t30*t31*t45*t51+1.3946639630399998d-5*
     7         t31*t35*t50+6.2242d-3*t11*t44*t45)*wght+Mmat3(iq,D3_RA_GA
     8         A_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_RB_GAA_TA) = Mmat3(iq,D3_RB_GAA_TA)
            Mmat3(iq,D3_RB_GAA_TB) = Mmat3(iq,D3_RB_GAA_TB)
            Mmat3(iq,D3_RB_GAB_TA) = Mmat3(iq,D3_RB_GAB_TA)
            Mmat3(iq,D3_RB_GAB_TB) = Mmat3(iq,D3_RB_GAB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = (t23*(1.0d+0*t18*(-1.39466396303999
     1         98d-4*gammabb*t12*t68*t98*taub-2.7893279260799997d-4*t13*
     2         t21*t53*t68*taub-4.6488798768d-5*t16*t38*t53*taub+4.97936
     3         d-2*t13*t38*t53*taub+2.48968d-2*t12*t38*t53*taub-1.673596
     4         755648d-4*t116*t24*t37*t8*t98-4.044525492816d-4*t24*t68*t
     5         8*t88+4.481424d-2*gammabb*t37*t53*t68*t8+1.04566559999999
     6         99d-1*gammabb*t38*t46*t8-4.1839918891199995d-5*t16*t21*t3
     7         7*t54*t68+4.481424d-2*t13*t21*t37*t54*t68+2.240712d-2*t12
     8         *t21*t37*t54*t68-4.1839918891199995d-5*gammabb*t15*t37*t5
     9         3*t68+2.240712d-2*gammabb*t12*t37*t53*t68-2.9287943223839
     :         996d-4*gammabb*t12*t21*t46*t68-1.8130631519519996d-4*t13*
     ;         t25*t27*t68+7.46904d-3*t16*t37*t38*t54+7.46904d-3*t15*t37
     <         *t38*t54-4.0d+0*t12*t37*t38*t54-1.673596755648d-4*t116*t1
     =         3*t27*t37*t54-1.673596755648d-4*gammabb*t116*t12*t21*t37*
     >         t53-9.76264774128d-5*gammabb*t15*t38*t46+5.228328d-2*gamm
     ?         abb*t12*t38*t46-6.043543839839999d-5*t16*t21*t25*t38+6.47
     @         3168d-2*t13*t21*t25*t38+3.236584d-2*t12*t21*t25*t38+1.618
     1         292d-2*t16*t25*t26+1.618292d-2*t15*t25*t26-8.666666666666
     2         666d+0*t12*t25*t26)+1.3333333333333333d+0*t41*t99)-1.0d+0
     3         *t18*t26*t37*t99-1.86726d-3*t19*t26*t94-6.043543839839999
     4         d-5*t28*t38*t90-3.73452d-3*t20*t26*t84+4.9793599999999993
     5         d-3*t26*t55*t81-4.1839918891199995d-5*t100*t28*t37*t68+7.
     6         46904d-3*t37*t38*t47*t56+3.73452d-3*t37*t38*t48*t55+1.394
     7         6639630399998d-5*t38*t42*t54+6.2242d-3*t26*t47*t48)*wght+
     8         Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = (2.7893279260799997d-5*t31*t44*t97-3
     1         .73452d-3*t11*t45*t96+1.0d+0*t1*t6*(1.673596755648d-4*gam
     2         maaa*t58*t8*t95-7.812600694878209d-8*gammaaa*t15*t58*t95+
     3         4.1839918891199995d-5*gammaaa*t12*t58*t95-3.1250402779512
     4         84d-7*gammaaa*t112*t12*t4*t95-3.125040277951284d-7*t112*t
     5         8*t9/rhoa**1.2333333333333334d+1-1.4938079999999998d-2*t3
     6         1*t49*t8-7.812600694878209d-8*t16*t4*t49*t58+8.3679837782
     7         39999d-5*t13*t4*t49*t58+8.367983778239999d-5*t12*t4*t49*t
     8         58+1.3946639630399998d-5*t16*t31*t49+2.7893279260799997d-
     9         5*t15*t31*t49-1.4938079999999998d-2*t12*t31*t49-3.1250402
     :         77951284d-7*t112*t13*t14*t49)+6.973319815199999d-6*t31*t5
     ;         1*t87-3.73452d-3*t11*t52*t86-7.812600694878209d-8*t113*t1
     <         7*t58)*wght+Mmat3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAA_GBB_TB) = (-3.73452d-3*t26*t48*t99+1.0d+0*t1
     1         8*t23*(1.673596755648d-4*gammabb*t68*t8*t98-7.81260069487
     2         8209d-8*gammabb*t15*t68*t98+4.1839918891199995d-5*gammabb
     3         *t12*t68*t98-3.125040277951284d-7*gammabb*t116*t12*t21*t9
     4         8-1.4938079999999998d-2*t38*t53*t8-3.125040277951284d-7*t
     5         116*t24*t8/rhob**1.2333333333333334d+1-7.812600694878209d
     6         -8*t16*t21*t53*t68+8.367983778239999d-5*t13*t21*t53*t68+8
     7         .367983778239999d-5*t12*t21*t53*t68+1.3946639630399998d-5
     8         *t16*t38*t53+2.7893279260799997d-5*t15*t38*t53-1.49380799
     9         99999998d-2*t12*t38*t53-3.125040277951284d-7*t116*t13*t27
     :         *t53)+6.973319815199999d-6*t38*t55*t90-3.73452d-3*t26*t56
     ;         *t89-7.812600694878209d-8*t117*t28*t68+2.7893279260799997
     <         d-5*t100*t38*t47)*wght+Mmat3(iq,D3_GAA_GBB_TB)
            Mmat3(iq,D3_GAB_GAB_TA) = (2.7893279260799997d-5*t31*t51*t97
     1         -7.46904d-3*t11*t52*t96+1.0d+0*t1*t6*(-6.250080555902569d
     2         -7*t112*t8*t9/rhoa**1.1333333333333334d+1+1.673596755648d
     3         -4*gammaaa*t101*t58*t8-1.5625201389756418d-7*t102*t16*t4*
     4         t58+3.347193511296d-4*t102*t13*t4*t58+8.367983778239999d-
     5         5*t102*t12*t4*t58-1.5625201389756418d-7*gammaaa*t101*t15*
     6         t58+1.673596755648d-4*gammaaa*t101*t12*t58-6.250080555902
     7         569d-7*gammaaa*t101*t112*t12*t4+5.578655852159999d-5*t102
     8         *t16*t31+2.7893279260799997d-5*t102*t15*t31-2.98761599999
     9         99996d-2*t102*t13*t31-2.9876159999999996d-2*t102*t12*t31-
     :         6.250080555902569d-7*t102*t112*t13*t14)-1.562520138975641
     ;         8d-7*t17*t58*t59-1.86726d-3*t104*t11*t45+2.78932792607999
     <         97d-5*t105*t31*t44)*wght+Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = (-7.46904d-3*t26*t56*t99+1.0d+0*t18
     1         *t23*(1.673596755648d-4*gammabb*t106*t68*t8-6.25008055590
     2         2569d-7*t116*t24*t8/rhob**1.1333333333333334d+1-1.5625201
     3         389756418d-7*t107*t16*t21*t68+3.347193511296d-4*t107*t13*
     4         t21*t68+8.367983778239999d-5*t107*t12*t21*t68-1.562520138
     5         9756418d-7*gammabb*t106*t15*t68+1.673596755648d-4*gammabb
     6         *t106*t12*t68+5.578655852159999d-5*t107*t16*t38+2.7893279
     7         260799997d-5*t107*t15*t38-2.9876159999999996d-2*t107*t13*
     8         t38-2.9876159999999996d-2*t107*t12*t38-6.250080555902569d
     9         -7*t107*t116*t13*t27-6.250080555902569d-7*gammabb*t106*t1
     :         16*t12*t21)-1.5625201389756418d-7*t28*t68*t69+2.789327926
     ;         0799997d-5*t100*t38*t55-1.86726d-3*t109*t26*t48+2.7893279
     <         260799997d-5*t110*t38*t47)*wght+Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_GAB_TA_TA) = (1.0d+0*t1*t6*(-1.2500161111805136d
     1         -6*t112*t8*t9/rhoa**1.0333333333333333d+1-3.1250402779512
     2         837d-7*t119*t16*t4*t58+1.0041580533887999d-3*t119*t13*t4*
     3         t58-3.1250402779512837d-7*gammaaa*t121*t15*t58+5.02079026
     4         6944d-4*gammaaa*t12*t121*t58-1.2500161111805136d-6*gammaa
     5         a*t112*t12*t121*t4+1.6735967556479997d-4*t119*t16*t31-1.7
     6         925696d-1*t119*t13*t31-1.2500161111805136d-6*t112*t119*t1
     7         3*t14)-3.1250402779512837d-7*t17*t29*t58-1.120356d-2*t104
     8         *t11*t52+8.367983778239998d-5*t105*t31*t51)*wght+Mmat3(iq
     9         ,D3_GAB_TA_TA)
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = (1.0d+0*t18*t23*(-1.250016111180513
     1         6d-6*t116*t24*t8/rhob**1.0333333333333333d+1-3.1250402779
     2         512837d-7*t120*t16*t21*t68+1.0041580533887999d-3*t120*t13
     3         *t21*t68-3.1250402779512837d-7*gammabb*t122*t15*t68+5.020
     4         790266944d-4*gammabb*t12*t122*t68+1.6735967556479997d-4*t
     5         120*t16*t38-1.7925696d-1*t120*t13*t38-1.2500161111805136d
     6         -6*t116*t120*t13*t27-1.2500161111805136d-6*gammabb*t116*t
     7         12*t122*t21)-3.1250402779512837d-7*t28*t36*t68-1.120356d-
     8         2*t109*t26*t56+8.367983778239998d-5*t110*t38*t55)*wght+Mm
     9         at3(iq,D3_GBB_TA_TA)
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
            t20 = 1/rhoa**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammaaa*t10*t11*t9*taua-3.33333
     1         33333333337d+0*t15*t2*t6*taua-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taua-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammaaa*t14*t17*t6-2*gamm
     4         aaa*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammaaa*t10*t11*t17*t4-gammaaa*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhoa**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22*wght+1.0d+0*t1*t21*wght
            t24 = 1/rhoa**8
            t25 = 2*gammaaa*t10*t7*t9-3.73452d-3*gammaaa*t11*t19*t4*t9-1
     1         .86726d-3*gammaaa*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhoa**7
            t28 = 1/rhoa**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammaaa*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammaaa*
     3         t10*t14*t28+2.0d+0*gammaaa*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhoa**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taua+1.8257653333333332d-2*g
     1         ammaaa*t33
            t35 = 1/rhoa**7.333333333333333d+0
            t36 = taua**2
            t37 = 1.3333333333333333d+1*gammaaa*t11*t18*t19*t9*taua+8.88
     1         888888888889d+0*t15*t17*t6*taua+2.6666666666666666d+1*t12
     2         *t18*t19*t2*t4*taua+1.7777777777777778d+1*t10*t12*t17*t4*
     3         taua+2.6666666666666666d+1*gammaaa*t10*t11*t20*taua+6.666
     4         666666666667d+0*t10*t15*t18*t2*taua-2*t19*t34*t7*t8*t9+6*
     5         t31*t32*t7*t8*t9+2.2222222222222224d+1*t10*t12*t36*t9+3.3
     6         77777777777778d+1*t10*t35*t7*t8+2.1333333333333332d+1*t18
     7         *t19*t20*t7*t8+9.777777777777779d+0*gammaaa*t14*t33*t6-2*
     8         gammaaa*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.77777777777777
     9         9d+0*gammaaa*t10*t11*t33*t4+6*gammaaa*t11*t2*t31*t32*t4+2
     :         *t15*t19*t31*t4+1.0666666666666666d+1*gammaaa*t11*t17*t18
     ;         *t19*t4-gammaaa*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t1
     <         3*t31*t32+2*gammaaa*t14*t19*t2*t31+5.333333333333333d+0*g
     =         ammaaa*t10*t14*t17*t18
            t38 = 1/rhoa**6.666666666666666d-1
            t39 = 4.444444444444444d-1*t16*t38*wght+1.0d+0*t1*t37*wght+2
     1         .6666666666666666d+0*t21*t22*wght
            t40 = -1.0d+0*t1*t16*t34*wght-1.3333333333333333d+0*t16*t18*
     1         t22*wght-1.0d+0*t1*t18*t21*wght-t18*t23
            t41 = 1/rhoa**2.3333333333333334d+0
            t42 = 1/rhoa**9
            t43 = 2.48968d-2*t12*t19*t4*t9*taua+6.2242d-3*t10*t15*t9*tau
     1         a-3.3333333333333337d+0*t10*t11*t9*taua+1.24484d-2*gammaa
     2         a*t11*t19*t24*taua-4*gammaaa*t18*t19*t7*t9+1.120356d-2*ga
     3         mmaaa*t11*t18*t32*t4*t9+3.73452d-3*gammaaa*t14*t18*t19*t9
     4         +2.9876159999999996d-2*t19*t42*t7*t8+1.120356d-2*t18*t24*
     5         t32*t7*t8-1.0666666666666666d+1*gammaaa*t10*t20*t7-2.6666
     6         666666666666d+0*t14*t17*t6+1.9917439999999997d-2*gammaaa*
     7         t11*t19*t20*t4+3.73452d-3*t15*t18*t19*t2*t4-2*t11*t18*t19
     8         *t2*t4+4.9793599999999993d-3*t10*t15*t17*t4-2.66666666666
     9         66666d+0*t10*t11*t17*t4+1.120356d-2*t12*t13*t18*t2*t32+9.
     :         958719999999999d-3*gammaaa*t10*t14*t20-t10*t14*t18*t2+9.9
     ;         58719999999999d-3*t12*t13*t17*t19
            t44 = 1.0d+0*t1*t43*wght+1.3333333333333333d+0*t22*t25*wght
            t45 = 1/rhoa**1.0666666666666666d+1
            t46 = 2*t10*t7*t9+6.973319815199999d-6*t15*t19*t4*t9-7.46904
     1         d-3*t11*t19*t4*t9+2.0919959445599998d-5*t12*t13*t32*t9-3.
     2         73452d-3*t10*t14*t9+2.0919959445599998d-5*t32*t45*t7*t8-1
     3         .4938079999999998d-2*gammaaa*t19*t24*t7+2.091995944559999
     4         8d-5*gammaaa*t11*t24*t32*t4+6.973319815199999d-6*gammaaa*
     5         t14*t19*t24
            t47 = 1/rhoa**4
            t48 = 4.97936d-2*t12*t19*t28*t4*taua+1.24484d-2*t10*t15*t28*
     1         taua-1.3333333333333333d+1*t10*t12*t28*taua+2.48968d-2*ga
     2         mmaaa*t11*t19*t27*taua+3.236584d-2*gammaaa*t11*t19*t4*t9+
     3         1.618292d-2*gammaaa*t10*t14*t9-8.666666666666666d+0*gamma
     4         aa*t10*t11*t9+2.240712d-2*t18*t27*t32*t7*t8+5.228328d-2*t
     5         19*t24*t7*t8-3.3333333333333337d+0*t15*t2*t6+2.240712d-2*
     6         gammaaa*t11*t18*t28*t32*t4+7.46904d-3*t15*t18*t19*t3*t4-8
     7         .0d+0*t12*t18*t19*t3*t4+6.2242d-3*t10*t15*t2*t4-6.6666666
     8         66666667d+0*t10*t12*t2*t4+2.240712d-2*t12*t13*t18*t3*t32-
     9         2.0d+0*t10*t15*t18*t3+7.46904d-3*gammaaa*t14*t18*t19*t28-
     :         4.0d+0*gammaaa*t11*t18*t19*t28+1.24484d-2*t12*t13*t19*t2
            t49 = 1.0d+0*t1*t48*wght+1.3333333333333333d+0*t22*t29*wght
            t50 = 1/rhoa**9.666666666666666d+0
            t51 = 4.1839918891199995d-5*t32*t50*t7*t8-1.4938079999999998
     1         d-2*gammaaa*t19*t27*t7+4.1839918891199995d-5*gammaaa*t11*
     2         t27*t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.4938079
     3         999999998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.1
     4         839918891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28
     5         -3.73452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.394663963039
     6         9998d-5*gammaaa*t14*t19*t27-7.46904d-3*gammaaa*t11*t19*t2
     7         7
            t52 = 1/rhoa**3
            t53 = 1/rhoa**8.666666666666666d+0
            t54 = 1/rhoa**6
            t55 = 1/rhoa**3.3333333333333337d+0
            t56 = 8.367983778239999d-5*t32*t53*t7*t8+2.7893279260799997d
     1         -5*t15*t19*t4*t55-5.975231999999999d-2*t12*t19*t4*t55+8.3
     2         67983778239999d-5*t12*t13*t32*t55-1.4938079999999998d-2*t
     3         10*t15*t55+8.0d+0*t10*t12*t55+8.367983778239999d-5*gammaa
     4         a*t11*t32*t4*t54+2.7893279260799997d-5*gammaaa*t14*t19*t5
     5         4-2.9876159999999996d-2*gammaaa*t11*t19*t54
            t57 = 1/rhoa**2
            t58 = t18**3
            t59 = 1/t5**5
            t60 = 1/rhoa**5.666666666666667d+0
            t61 = -6.085884444444444d-2*t33*taua-8.520238222222222d-2*ga
     1         mmaaa*t60
            t62 = 1/rhoa**5
            t63 = 1/rhoa**7.666666666666667d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = -1.0d+0*t1*t10*t16*t18*wght+t23*t6+Amat(iq,
     1         D1_RA)
            Cmat(iq,D1_GAA) = 1.0d+0*t1*t25*t6*wght-1.86726d-3*t10*t16*t
     1         26*wght+Cmat(iq,D1_GAA)
            Mmat(iq,D1_TA) = 1.0d+0*t1*t29*t6*wght-3.73452d-3*t10*t16*t3
     1         0*wght+Mmat(iq,D1_TA)
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = 2.0d+0*t1*t16*t19*t31*wght+t39*t6+t10*t
     1         40+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 4.9793599999999993d-3*t10*t16*t41*wght
     1         +3.73452d-3*t16*t18*t19*t26*wght-1.0d+0*t1*t10*t18*t25*wg
     2         ht+t44*t6-1.86726d-3*t10*t2*t23+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.0d+0*t1*t46*t6*wght+6.9733198151999
     1         99d-6*t16*t19*t47*wght-3.73452d-3*t10*t25*t26*wght+Cmat2(
     2         iq,D2_GAA_GAA)
            Mmat2(iq,D2_RA_TA) = 7.46904d-3*t16*t18*t19*t30*wght-1.0d+0*
     1         t1*t10*t18*t29*wght+6.2242d-3*t10*t16*t26*wght+t49*t6-3.7
     2         3452d-3*t10*t23*t3+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 1.0d+0*t1*t51*t6*wght+1.39466396303999
     1         98d-5*t16*t19*t52*wght-3.73452d-3*t10*t25*t30*wght-1.8672
     2         6d-3*t10*t26*t29*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = 1.0d+0*t1*t56*t6*wght+2.789327926079999
     1         7d-5*t16*t19*t57*wght-7.46904d-3*t10*t29*t30*wght+Mmat2(i
     2         q,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = t6*(1.0d+0*t1*(2.0d+1*gammaaa*t11*t1
     1         9*t34*t9*taua-6.0d+1*gammaaa*t11*t31*t32*t9*taua-3.259259
     2         2592592595d+1*t15*t33*t6*taua+4.0d+1*t12*t19*t2*t34*t4*ta
     3         ua-6.518518518518519d+1*t10*t12*t33*t4*taua-1.2d+2*t12*t2
     4         *t31*t32*t4*taua-1.0666666666666667d+2*t12*t17*t18*t19*t4
     5         *taua-2.0148148148148146d+2*gammaaa*t10*t11*t35*taua+1.0d
     6         +1*t10*t15*t2*t34*taua-2.0d+1*t15*t19*t2*t31*taua-1.6d+2*
     7         gammaaa*t11*t18*t19*t20*taua-2.6666666666666666d+1*t10*t1
     8         5*t17*t18*taua-2*t19*t61*t7*t8*t9-24*t58*t59*t7*t8*t9+18*
     9         t18*t32*t34*t7*t8*t9-1.3333333333333333d+2*t12*t18*t19*t3
     :         6*t9-2.0266666666666666d+2*t18*t19*t35*t7*t8+32*t19*t20*t
     ;         34*t7*t8-96*t20*t31*t32*t7*t8-2.477037037037037d+2*t10*t7
     <         *t8/rhoa**8.333333333333334d+0-2*gammaaa*t11*t19*t2*t4*t6
     =         1-t10*t15*t4*t61-gammaaa*t10*t14*t2*t61-2*t12*t13*t19*t61
     >         -4.5629629629629626d+1*gammaaa*t14*t6*t60-4.5629629629629
     ?         626d+1*gammaaa*t10*t11*t4*t60-24*gammaaa*t11*t2*t4*t58*t5
     @         9-24*t12*t13*t58*t59-6*t15*t32*t4*t58-6*gammaaa*t14*t2*t3
     1         2*t58+18*gammaaa*t11*t18*t2*t32*t34*t4+6*t15*t18*t19*t34*
     2         t4+16*gammaaa*t11*t17*t19*t34*t4-5.866666666666667d+1*gam
     3         maaa*t11*t18*t19*t33*t4-48*gammaaa*t11*t17*t31*t32*t4-1.7
     4         77777777777778d+2*t10*t12*t20*t36+18*t12*t13*t18*t32*t34+
     5         6*gammaaa*t14*t18*t19*t2*t34+8*gammaaa*t10*t14*t17*t34-2.
     6         9333333333333333d+1*gammaaa*t10*t14*t18*t33-16*gammaaa*t1
     7         4*t17*t19*t31)*wght+1.3333333333333333d+0*t21*t38*wght+4.
     8         0d+0*t22*t37*wght-2.962962962962963d-1*t16*t3*wght)+t10*(
     9         -1.0d+0*t1*t16*t61*wght-4.444444444444444d-1*t16*t18*t38*
     :         wght-1.0d+0*t1*t18*t37*wght-2.6666666666666666d+0*t16*t22
     ;         *t34*wght-2.0d+0*t1*t21*t34*wght-2.6666666666666666d+0*t1
     <         8*t21*t22*wght-2*t18*t39-t23*t34)+t19*(4.0d+0*t1*t16*t18*
     =         t34*wght+2.6666666666666666d+0*t16*t22*t31*wght+2.0d+0*t1
     >         *t21*t31*wght-2*t18*t40)-6.0d+0*t1*t16*t32*t58*wght+Amat3
     ?         (iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t6*(1.0d+0*t1*(-1.493808d-1*t12*t18
     1         *t32*t4*t9*taua-2.48968d-2*t15*t18*t19*t9*taua+1.33333333
     2         33333333d+1*t11*t18*t19*t9*taua-1.6597866666666666d-1*gam
     3         maaa*t11*t19*t42*taua-1.9917439999999997d-1*t12*t19*t20*t
     4         4*taua-7.46904d-2*gammaaa*t11*t18*t24*t32*taua-4.97935999
     5         99999993d-2*t10*t15*t20*taua+2.6666666666666666d+1*t10*t1
     6         1*t20*taua-4*gammaaa*t19*t34*t7*t9+12*gammaaa*t31*t32*t7*
     7         t9-4.481424d-2*gammaaa*t11*t31*t4*t59*t9+1.120356d-2*gamm
     8         aaa*t11*t32*t34*t4*t9+3.73452d-3*gammaaa*t14*t19*t34*t9-1
     9         .120356d-2*gammaaa*t14*t31*t32*t9-4.481424d-2*t24*t31*t59
     :         *t7*t8-1.7925696d-1*t18*t32*t42*t7*t8+1.120356d-2*t24*t32
     ;         *t34*t7*t8-2.6888544000000003d-1*t19*t7*t8/rhoa**10+6.755
     <         555555555556d+1*gammaaa*t10*t35*t7+4.2666666666666664d+1*
     =         gammaaa*t18*t19*t20*t7+9.777777777777779d+0*t14*t33*t6-4.
     >         481424d-2*t12*t13*t2*t31*t59-1.2614378666666665d-1*gammaa
     ?         a*t11*t19*t35*t4+3.73452d-3*t15*t19*t2*t34*t4-2*t11*t19*t
     @         2*t34*t4-1.8257653333333332d-2*t10*t15*t33*t4+9.777777777
     1         777779d+0*t10*t11*t33*t4-1.120356d-2*t15*t2*t31*t32*t4+6*
     2         t11*t2*t31*t32*t4-1.1950464d-1*gammaaa*t11*t18*t20*t32*t4
     3         -1.9917439999999997d-2*t15*t17*t18*t19*t4+1.0666666666666
     4         666d+1*t11*t17*t18*t19*t4-8.298933333333335d-2*t12*t19*t2
     5         4*t36-6.307189333333332d-2*gammaaa*t10*t14*t35+1.120356d-
     6         2*t12*t13*t2*t32*t34-t10*t14*t2*t34-3.6515306666666664d-2
     7         *t12*t13*t19*t33-5.975231999999999d-2*t12*t13*t17*t18*t32
     8         +2*t14*t19*t2*t31-3.9834879999999995d-2*gammaaa*t14*t18*t
     9         19*t20+5.333333333333333d+0*t10*t14*t17*t18)*wght+2.66666
     :         66666666666d+0*t22*t43*wght+4.444444444444444d-1*t25*t38*
     ;         wght)+t10*(-1.1618506666666667d-2*t16*t55*wght-1.0d+0*t1*
     <         t18*t43*wght+4.9793599999999993d-3*t21*t41*wght-1.0d+0*t1
     =         *t25*t34*wght-1.3333333333333333d+0*t18*t22*t25*wght-t18*
     >         t44+4.9793599999999993d-3*t17*t23)-1.9917439999999997d-2*
     ?         t16*t18*t19*t41*wght-1.120356d-2*t16*t26*t31*t32*wght+2.0
     @         d+0*t1*t19*t25*t31*wght-3.73452d-3*t19*t2*t40-1.86726d-3*
     1         t10*t2*t39+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t6*(1.0d+0*t1*(-6.9733198152d-5*ga
     1         mmaaa*t11*t32*t45*taua-1.3946639630399998d-4*t12*t24*t32*
     2         t4*taua-2.3244399384d-5*t15*t19*t24*taua+2.48968d-2*t11*t
     3         19*t24*taua-4*t18*t19*t7*t9-8.367983778239999d-5*t12*t13*
     4         t18*t59*t9-2.0919959445599998d-5*t15*t18*t32*t4*t9+2.2407
     5         12d-2*t11*t18*t32*t4*t9+7.46904d-3*t14*t18*t19*t9-8.36798
     6         3778239999d-5*t18*t45*t59*t7*t8-2.2314623408639997d-4*t32
     7         *t7*t8/rhoa**1.1666666666666665d+1+1.1950464d-1*gammaaa*t
     8         19*t42*t7+4.481424d-2*gammaaa*t18*t24*t32*t7-1.0666666666
     9         666666d+1*t10*t20*t7-8.367983778239999d-5*gammaaa*t11*t18
     :         *t24*t4*t59-1.6735967556479997d-4*gammaaa*t11*t32*t4*t42-
     ;         5.578655852159999d-5*gammaaa*t14*t19*t42-3.71910390143999
     <         96d-5*t15*t19*t20*t4+3.9834879999999995d-2*t11*t19*t20*t4
     =         -2.0919959445599998d-5*gammaaa*t14*t18*t24*t32-1.11573117
     >         04319998d-4*t12*t13*t20*t32+1.9917439999999997d-2*t10*t14
     ?         *t20)*wght+1.3333333333333333d+0*t22*t46*wght)-3.71910390
     @         14399996d-5*t16*t19*t62*wght-2.0919959445599998d-5*t16*t1
     1         8*t32*t47*wght-1.0d+0*t1*t10*t18*t46*wght+9.9587199999999
     2         99d-3*t10*t25*t41*wght+7.46904d-3*t18*t19*t25*t26*wght+6.
     3         973319815199999d-6*t19*t23*t9-3.73452d-3*t10*t2*t44+Cmat3
     4         (iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 1.0d+0*t1*t6*(-1.562520138975642d
     1         -7*t59*t7*t8/rhoa**1.3333333333333333d+1+1.255197566736d-
     2         4*gammaaa*t32*t45*t7-2.240712d-2*t19*t24*t7-1.56252013897
     3         5642d-7*gammaaa*t11*t4*t45*t59-1.562520138975642d-7*t12*t
     4         13*t24*t59-3.9063003474391045d-8*gammaaa*t14*t32*t45-3.90
     5         63003474391045d-8*t15*t24*t32*t4+6.275987833679999d-5*t11
     6         *t24*t32*t4+2.0919959445599995d-5*t14*t19*t24)*wght+2.091
     7         9959445599995d-5*t19*t25*t47*wght-5.60178d-3*t10*t26*t46*
     8         wght-3.9063003474391045d-8*t16*t32*wght/rhoa**6.666666666
     9         666667d+0+Cmat3(iq,D3_GAA_GAA_GAA)
            Mmat3(iq,D3_RA_RA_TA) = t6*(1.0d+0*t1*(-2.9876159999999996d-
     1         1*t12*t19*t4*t9*taua-7.469039999999999d-2*t10*t15*t9*taua
     2         +8.0d+1*t10*t12*t9*taua-2.987616d-1*t12*t18*t28*t32*t4*ta
     3         ua-1.493808d-1*gammaaa*t11*t18*t27*t32*taua-4.97936d-2*t1
     4         5*t18*t19*t28*taua+5.333333333333333d+1*t12*t18*t19*t28*t
     5         aua-2.8216373333333333d-1*gammaaa*t11*t19*t24*taua-1.9419
     6         503999999999d-1*gammaaa*t11*t18*t32*t4*t9-6.473168d-2*gam
     7         maaa*t14*t18*t19*t9+3.466666666666666d+1*gammaaa*t11*t18*
     8         t19*t9-8.962848d-2*t27*t31*t59*t7*t8-4.1826624d-1*t19*t42
     9         *t7*t8+2.240712d-2*t27*t32*t34*t7*t8-3.1369968d-1*t18*t24
     :         *t32*t7*t8+8.88888888888889d+0*t15*t17*t6-8.962848d-2*gam
     ;         maaa*t11*t28*t31*t4*t59-8.962848d-2*t12*t13*t3*t31*t59+2.
     <         240712d-2*gammaaa*t11*t28*t32*t34*t4+7.46904d-3*t15*t19*t
     =         3*t34*t4-8.0d+0*t12*t19*t3*t34*t4-2.240712d-2*t15*t3*t31*
     >         t32*t4+2.4d+1*t12*t3*t31*t32*t4-1.7261781333333331d-1*gam
     ?         maaa*t11*t19*t20*t4-2.48968d-2*t15*t18*t19*t2*t4+2.666666
     @         6666666666d+1*t12*t18*t19*t2*t4-1.6597866666666666d-2*t10
     1         *t15*t17*t4+1.7777777777777778d+1*t10*t12*t17*t4-1.659786
     2         666666667d-1*t12*t19*t27*t36+2.240712d-2*t12*t13*t3*t32*t
     3         34-2.0d+0*t10*t15*t3*t34+7.46904d-3*gammaaa*t14*t19*t28*t
     4         34-4.0d+0*gammaaa*t11*t19*t28*t34-2.240712d-2*gammaaa*t14
     5         *t28*t31*t32+1.2d+1*gammaaa*t11*t28*t31*t32-7.46904d-2*t1
     6         2*t13*t18*t2*t32+4.0d+0*t15*t19*t3*t31-8.630890666666666d
     7         -2*gammaaa*t10*t14*t20+4.622222222222223d+1*gammaaa*t10*t
     8         11*t20+6.666666666666667d+0*t10*t15*t18*t2-3.319573333333
     9         333d-2*t12*t13*t17*t19)*wght+2.6666666666666666d+0*t22*t4
     :         8*wght+4.444444444444444d-1*t29*t38*wght)+t10*(-1.0d+0*t1
     ;         *t18*t48*wght-8.298933333333333d-3*t16*t41*wght-1.0d+0*t1
     <         *t29*t34*wght-1.3333333333333333d+0*t18*t22*t29*wght+6.22
     =         42d-3*t21*t26*wght-t18*t49+6.2242d-3*t2*t23)-2.240712d-2*
     >         t16*t30*t31*t32*wght+2.0d+0*t1*t19*t29*t31*wght-2.48968d-
     ?         2*t16*t18*t19*t26*wght-7.46904d-3*t19*t3*t40-3.73452d-3*t
     @         10*t3*t39+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_GAA_TA) = t6*(1.0d+0*t1*(-1.3946639630399998d
     1         -4*gammaaa*t11*t32*t50*taua-2.7893279260799997d-4*t12*t27
     2         *t32*t4*taua-4.6488798768d-5*t15*t19*t27*taua+4.97936d-2*
     3         t12*t19*t27*taua+2.48968d-2*t11*t19*t27*taua-6.0435438398
     4         39999d-5*t15*t19*t4*t9+6.473168d-2*t12*t19*t4*t9+3.236584
     5         d-2*t11*t19*t4*t9-1.8130631519519996d-4*t12*t13*t32*t9+1.
     6         618292d-2*t10*t15*t9+1.618292d-2*t10*t14*t9-8.66666666666
     7         6666d+0*t10*t11*t9-1.673596755648d-4*t18*t50*t59*t7*t8-4.
     8         044525492816d-4*t32*t45*t7*t8+4.481424d-2*gammaaa*t18*t27
     9         *t32*t7+1.0456655999999999d-1*gammaaa*t19*t24*t7-1.673596
     :         755648d-4*gammaaa*t11*t18*t27*t4*t59-1.673596755648d-4*t1
     ;         2*t13*t18*t28*t59-4.1839918891199995d-5*t15*t18*t28*t32*t
     <         4+4.481424d-2*t12*t18*t28*t32*t4+2.240712d-2*t11*t18*t28*
     =         t32*t4-2.9287943223839996d-4*gammaaa*t11*t24*t32*t4-4.183
     >         9918891199995d-5*gammaaa*t14*t18*t27*t32+2.240712d-2*gamm
     ?         aaa*t11*t18*t27*t32+7.46904d-3*t15*t18*t19*t28+7.46904d-3
     @         *t14*t18*t19*t28-4.0d+0*t11*t18*t19*t28-9.76264774128d-5*
     1         gammaaa*t14*t19*t24+5.228328d-2*gammaaa*t11*t19*t24)*wght
     2         +1.3333333333333333d+0*t22*t51*wght)-4.1839918891199995d-
     3         5*t16*t18*t32*t52*wght-1.0d+0*t1*t10*t18*t51*wght-6.04354
     4         3839839999d-5*t16*t19*t47*wght+4.9793599999999993d-3*t10*
     5         t29*t41*wght+7.46904d-3*t18*t19*t25*t30*wght+3.73452d-3*t
     6         18*t19*t26*t29*wght+6.2242d-3*t10*t25*t26*wght-1.86726d-3
     7         *t10*t2*t49-3.73452d-3*t10*t3*t44+1.3946639630399998d-5*t
     8         19*t23*t28+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = 1.0d+0*t1*t6*(-3.125040277951284d-
     1         7*t59*t7*t8/rhoa**1.2333333333333334d+1+1.673596755648d-4
     2         *gammaaa*t32*t50*t7-1.4938079999999998d-2*t19*t27*t7-3.12
     3         5040277951284d-7*gammaaa*t11*t4*t50*t59-3.125040277951284
     4         d-7*t12*t13*t27*t59-7.812600694878209d-8*gammaaa*t14*t32*
     5         t50+4.1839918891199995d-5*gammaaa*t11*t32*t50-7.812600694
     6         878209d-8*t15*t27*t32*t4+8.367983778239999d-5*t12*t27*t32
     7         *t4+8.367983778239999d-5*t11*t27*t32*t4+1.394663963039999
     8         8d-5*t15*t19*t27+2.7893279260799997d-5*t14*t19*t27-1.4938
     9         079999999998d-2*t11*t19*t27)*wght-7.812600694878209d-8*t1
     :         6*t32*t60*wght+2.7893279260799997d-5*t19*t25*t52*wght-3.7
     ;         3452d-3*t10*t26*t51*wght+6.973319815199999d-6*t19*t29*t47
     <         *wght-3.73452d-3*t10*t30*t46*wght+Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_RA_TA_TA) = t6*(1.0d+0*t1*(-5.578655852159999d-4
     1         *t12*t32*t4*t54*taua-9.2977597536d-5*t15*t19*t54*taua+1.9
     2         91744d-1*t12*t19*t54*taua-2.7893279260799997d-4*gammaaa*t
     3         11*t32*t53*taua-3.347193511296d-4*t18*t53*t59*t7*t8-7.252
     4         252607808d-4*t32*t50*t7*t8-3.347193511296d-4*t12*t13*t18*
     5         t55*t59-3.347193511296d-4*gammaaa*t11*t18*t4*t54*t59-8.36
     6         7983778239999d-5*t15*t18*t32*t4*t55+1.7925696d-1*t12*t18*
     7         t32*t4*t55+2.9876159999999996d-2*t15*t18*t19*t55-1.6d+1*t
     8         12*t18*t19*t55-8.367983778239999d-5*gammaaa*t14*t18*t32*t
     9         54+8.962848d-2*gammaaa*t11*t18*t32*t54-5.020790266944d-4*
     :         gammaaa*t11*t27*t32*t4-9.2977597536d-5*t15*t19*t28*t4+1.9
     ;         91744d-1*t12*t19*t28*t4-2.7893279260799997d-4*t12*t13*t28
     <         *t32+4.97936d-2*t10*t15*t28-2.6666666666666666d+1*t10*t12
     =         *t28-1.673596755648d-4*gammaaa*t14*t19*t27+1.792569600000
     >         0002d-1*gammaaa*t11*t19*t27)*wght+1.3333333333333333d+0*t
     ?         22*t56*wght)-8.367983778239999d-5*t16*t18*t32*t57*wght-1.
     @         0d+0*t1*t10*t18*t56*wght-9.2977597536d-5*t16*t19*t52*wght
     1         +1.4938079999999998d-2*t18*t19*t29*t30*wght+1.24484d-2*t1
     2         0*t26*t29*wght+2.7893279260799997d-5*t19*t23*t55-7.46904d
     3         -3*t10*t3*t49+Mmat3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = 1.0d+0*t1*t6*(-6.250080555902569d-7
     1         *t59*t7*t8/rhoa**1.1333333333333334d+1+1.673596755648d-4*
     2         gammaaa*t32*t53*t7-6.250080555902569d-7*t12*t13*t54*t59-6
     3         .250080555902569d-7*gammaaa*t11*t4*t53*t59-1.562520138975
     4         6418d-7*t15*t32*t4*t54+3.347193511296d-4*t12*t32*t4*t54+8
     5         .367983778239999d-5*t11*t32*t4*t54+5.578655852159999d-5*t
     6         15*t19*t54+2.7893279260799997d-5*t14*t19*t54-2.9876159999
     7         999996d-2*t12*t19*t54-2.9876159999999996d-2*t11*t19*t54-1
     8         .5625201389756418d-7*gammaaa*t14*t32*t53+1.673596755648d-
     9         4*gammaaa*t11*t32*t53)*wght+2.7893279260799997d-5*t19*t25
     :         *t57*wght-1.86726d-3*t10*t26*t56*wght+2.7893279260799997d
     ;         -5*t19*t29*t52*wght-7.46904d-3*t10*t30*t51*wght-1.5625201
     <         389756418d-7*t16*t32*t33*wght+Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = 1.0d+0*t1*t6*(-1.2500161111805136d-6
     1         *t59*t7*t8/rhoa**1.0333333333333333d+1-1.2500161111805136
     2         d-6*gammaaa*t11*t4*t59*t63-3.1250402779512837d-7*gammaaa*
     3         t14*t32*t63+5.020790266944d-4*gammaaa*t11*t32*t63-1.25001
     4         61111805136d-6*t12*t13*t59*t62-3.1250402779512837d-7*t15*
     5         t32*t4*t62+1.0041580533887999d-3*t12*t32*t4*t62+1.6735967
     6         556479997d-4*t15*t19*t62-1.7925696d-1*t12*t19*t62)*wght+8
     7         .367983778239998d-5*t19*t29*t57*wght-1.120356d-2*t10*t30*
     8         t56*wght-3.1250402779512837d-7*t16*t17*t32*wght+Mmat3(iq,
     9         D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
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
            t20 = 1/rhob**6.333333333333333d+0
            t21 = -3.3333333333333337d+0*gammabb*t10*t11*t9*taub-3.33333
     1         33333333337d+0*t15*t2*t6*taub-6.666666666666667d+0*t10*t1
     2         2*t2*t4*taub-2*t18*t19*t7*t8*t9-5.333333333333333d+0*t10*
     3         t20*t7*t8-2.6666666666666666d+0*gammabb*t14*t17*t6-2*gamm
     4         abb*t11*t18*t19*t2*t4-t10*t15*t18*t4-2.6666666666666666d+
     5         0*gammabb*t10*t11*t17*t4-gammabb*t10*t14*t18*t2-2*t12*t13
     6         *t18*t19
            t22 = rhob**3.333333333333333d-1
            t23 = 1.3333333333333333d+0*t16*t22*wght+1.0d+0*t1*t21*wght
            t24 = 1/rhob**8
            t25 = 2*gammabb*t10*t7*t9-3.73452d-3*gammabb*t11*t19*t4*t9-1
     1         .86726d-3*gammabb*t10*t14*t9-3.73452d-3*t19*t24*t7*t8+t14
     2         *t2*t6-1.86726d-3*t10*t15*t2*t4+t10*t11*t2*t4-3.73452d-3*
     3         t12*t13*t19*t2
            t26 = 1/t1
            t27 = 1/rhob**7
            t28 = 1/rhob**4.333333333333333d+0
            t29 = -7.46904d-3*t19*t27*t7*t8+2.0d+0*t15*t3*t6-3.73452d-3*
     1         t10*t15*t3*t4+4.0d+0*t10*t12*t3*t4-7.46904d-3*gammabb*t11
     2         *t19*t28*t4-7.46904d-3*t12*t13*t19*t3-3.73452d-3*gammabb*
     3         t10*t14*t28+2.0d+0*gammabb*t10*t11*t28
            t30 = 1/t22
            t31 = t18**2
            t32 = 1/t5**4
            t33 = 1/rhob**4.666666666666667d+0
            t34 = 1.6597866666666666d-2*t17*taub+1.8257653333333332d-2*g
     1         ammabb*t33
            t35 = 1/rhob**7.333333333333333d+0
            t36 = taub**2
            t37 = 1.3333333333333333d+1*gammabb*t11*t18*t19*t9*taub+8.88
     1         888888888889d+0*t15*t17*t6*taub+2.6666666666666666d+1*t12
     2         *t18*t19*t2*t4*taub+1.7777777777777778d+1*t10*t12*t17*t4*
     3         taub+2.6666666666666666d+1*gammabb*t10*t11*t20*taub+6.666
     4         666666666667d+0*t10*t15*t18*t2*taub-2*t19*t34*t7*t8*t9+6*
     5         t31*t32*t7*t8*t9+2.2222222222222224d+1*t10*t12*t36*t9+3.3
     6         77777777777778d+1*t10*t35*t7*t8+2.1333333333333332d+1*t18
     7         *t19*t20*t7*t8+9.777777777777779d+0*gammabb*t14*t33*t6-2*
     8         gammabb*t11*t19*t2*t34*t4-t10*t15*t34*t4+9.77777777777777
     9         9d+0*gammabb*t10*t11*t33*t4+6*gammabb*t11*t2*t31*t32*t4+2
     :         *t15*t19*t31*t4+1.0666666666666666d+1*gammabb*t11*t17*t18
     ;         *t19*t4-gammabb*t10*t14*t2*t34-2*t12*t13*t19*t34+6*t12*t1
     <         3*t31*t32+2*gammabb*t14*t19*t2*t31+5.333333333333333d+0*g
     =         ammabb*t10*t14*t17*t18
            t38 = 1/rhob**6.666666666666666d-1
            t39 = 4.444444444444444d-1*t16*t38*wght+1.0d+0*t1*t37*wght+2
     1         .6666666666666666d+0*t21*t22*wght
            t40 = -1.0d+0*t1*t16*t34*wght-1.3333333333333333d+0*t16*t18*
     1         t22*wght-1.0d+0*t1*t18*t21*wght-t18*t23
            t41 = 1/rhob**2.3333333333333334d+0
            t42 = 1/rhob**9
            t43 = 2.48968d-2*t12*t19*t4*t9*taub+6.2242d-3*t10*t15*t9*tau
     1         b-3.3333333333333337d+0*t10*t11*t9*taub+1.24484d-2*gammab
     2         b*t11*t19*t24*taub-4*gammabb*t18*t19*t7*t9+1.120356d-2*ga
     3         mmabb*t11*t18*t32*t4*t9+3.73452d-3*gammabb*t14*t18*t19*t9
     4         +2.9876159999999996d-2*t19*t42*t7*t8+1.120356d-2*t18*t24*
     5         t32*t7*t8-1.0666666666666666d+1*gammabb*t10*t20*t7-2.6666
     6         666666666666d+0*t14*t17*t6+1.9917439999999997d-2*gammabb*
     7         t11*t19*t20*t4+3.73452d-3*t15*t18*t19*t2*t4-2*t11*t18*t19
     8         *t2*t4+4.9793599999999993d-3*t10*t15*t17*t4-2.66666666666
     9         66666d+0*t10*t11*t17*t4+1.120356d-2*t12*t13*t18*t2*t32+9.
     :         958719999999999d-3*gammabb*t10*t14*t20-t10*t14*t18*t2+9.9
     ;         58719999999999d-3*t12*t13*t17*t19
            t44 = 1.0d+0*t1*t43*wght+1.3333333333333333d+0*t22*t25*wght
            t45 = 1/rhob**1.0666666666666666d+1
            t46 = 2*t10*t7*t9+6.973319815199999d-6*t15*t19*t4*t9-7.46904
     1         d-3*t11*t19*t4*t9+2.0919959445599998d-5*t12*t13*t32*t9-3.
     2         73452d-3*t10*t14*t9+2.0919959445599998d-5*t32*t45*t7*t8-1
     3         .4938079999999998d-2*gammabb*t19*t24*t7+2.091995944559999
     4         8d-5*gammabb*t11*t24*t32*t4+6.973319815199999d-6*gammabb*
     5         t14*t19*t24
            t47 = 1/rhob**4
            t48 = 4.97936d-2*t12*t19*t28*t4*taub+1.24484d-2*t10*t15*t28*
     1         taub-1.3333333333333333d+1*t10*t12*t28*taub+2.48968d-2*ga
     2         mmabb*t11*t19*t27*taub+3.236584d-2*gammabb*t11*t19*t4*t9+
     3         1.618292d-2*gammabb*t10*t14*t9-8.666666666666666d+0*gamma
     4         bb*t10*t11*t9+2.240712d-2*t18*t27*t32*t7*t8+5.228328d-2*t
     5         19*t24*t7*t8-3.3333333333333337d+0*t15*t2*t6+2.240712d-2*
     6         gammabb*t11*t18*t28*t32*t4+7.46904d-3*t15*t18*t19*t3*t4-8
     7         .0d+0*t12*t18*t19*t3*t4+6.2242d-3*t10*t15*t2*t4-6.6666666
     8         66666667d+0*t10*t12*t2*t4+2.240712d-2*t12*t13*t18*t3*t32-
     9         2.0d+0*t10*t15*t18*t3+7.46904d-3*gammabb*t14*t18*t19*t28-
     :         4.0d+0*gammabb*t11*t18*t19*t28+1.24484d-2*t12*t13*t19*t2
            t49 = 1.0d+0*t1*t48*wght+1.3333333333333333d+0*t22*t29*wght
            t50 = 1/rhob**9.666666666666666d+0
            t51 = 4.1839918891199995d-5*t32*t50*t7*t8-1.4938079999999998
     1         d-2*gammabb*t19*t27*t7+4.1839918891199995d-5*gammabb*t11*
     2         t27*t32*t4+1.3946639630399998d-5*t15*t19*t28*t4-1.4938079
     3         999999998d-2*t12*t19*t28*t4-7.46904d-3*t11*t19*t28*t4+4.1
     4         839918891199995d-5*t12*t13*t28*t32-3.73452d-3*t10*t15*t28
     5         -3.73452d-3*t10*t14*t28+2.0d+0*t10*t11*t28+1.394663963039
     6         9998d-5*gammabb*t14*t19*t27-7.46904d-3*gammabb*t11*t19*t2
     7         7
            t52 = 1/rhob**3
            t53 = 1/rhob**8.666666666666666d+0
            t54 = 1/rhob**6
            t55 = 1/rhob**3.3333333333333337d+0
            t56 = 8.367983778239999d-5*t32*t53*t7*t8+2.7893279260799997d
     1         -5*t15*t19*t4*t55-5.975231999999999d-2*t12*t19*t4*t55+8.3
     2         67983778239999d-5*t12*t13*t32*t55-1.4938079999999998d-2*t
     3         10*t15*t55+8.0d+0*t10*t12*t55+8.367983778239999d-5*gammab
     4         b*t11*t32*t4*t54+2.7893279260799997d-5*gammabb*t14*t19*t5
     5         4-2.9876159999999996d-2*gammabb*t11*t19*t54
            t57 = 1/rhob**2
            t58 = t18**3
            t59 = 1/t5**5
            t60 = 1/rhob**5.666666666666667d+0
            t61 = -6.085884444444444d-2*t33*taub-8.520238222222222d-2*ga
     1         mmabb*t60
            t62 = 1/rhob**5
            t63 = 1/rhob**7.666666666666667d+0
            fnc(iq) = 1.0d+0*t1*t16*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = -1.0d+0*t1*t10*t16*t18*wght+t23*t6+Amat(iq,
     1         D1_RB)
            Cmat(iq,D1_GBB) = 1.0d+0*t1*t25*t6*wght-1.86726d-3*t10*t16*t
     1         26*wght+Cmat(iq,D1_GBB)
            Mmat(iq,D1_TB) = 1.0d+0*t1*t29*t6*wght-3.73452d-3*t10*t16*t3
     1         0*wght+Mmat(iq,D1_TB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = 2.0d+0*t1*t16*t19*t31*wght+t39*t6+t10*t
     1         40+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 4.9793599999999993d-3*t10*t16*t41*wght
     1         +3.73452d-3*t16*t18*t19*t26*wght-1.0d+0*t1*t10*t18*t25*wg
     2         ht+t44*t6-1.86726d-3*t10*t2*t23+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.0d+0*t1*t46*t6*wght+6.9733198151999
     1         99d-6*t16*t19*t47*wght-3.73452d-3*t10*t25*t26*wght+Cmat2(
     2         iq,D2_GBB_GBB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = 7.46904d-3*t16*t18*t19*t30*wght-1.0d+0*
     1         t1*t10*t18*t29*wght+6.2242d-3*t10*t16*t26*wght+t49*t6-3.7
     2         3452d-3*t10*t23*t3+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 1.0d+0*t1*t51*t6*wght+1.39466396303999
     1         98d-5*t16*t19*t52*wght-3.73452d-3*t10*t25*t30*wght-1.8672
     2         6d-3*t10*t26*t29*wght+Mmat2(iq,D2_GBB_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = 1.0d+0*t1*t56*t6*wght+2.789327926079999
     1         7d-5*t16*t19*t57*wght-7.46904d-3*t10*t29*t30*wght+Mmat2(i
     2         q,D2_TB_TB)
            Amat3(iq,D3_RB_RB_RB) = t6*(1.0d+0*t1*(2.0d+1*gammabb*t11*t1
     1         9*t34*t9*taub-6.0d+1*gammabb*t11*t31*t32*t9*taub-3.259259
     2         2592592595d+1*t15*t33*t6*taub+4.0d+1*t12*t19*t2*t34*t4*ta
     3         ub-6.518518518518519d+1*t10*t12*t33*t4*taub-1.2d+2*t12*t2
     4         *t31*t32*t4*taub-1.0666666666666667d+2*t12*t17*t18*t19*t4
     5         *taub-2.0148148148148146d+2*gammabb*t10*t11*t35*taub+1.0d
     6         +1*t10*t15*t2*t34*taub-2.0d+1*t15*t19*t2*t31*taub-1.6d+2*
     7         gammabb*t11*t18*t19*t20*taub-2.6666666666666666d+1*t10*t1
     8         5*t17*t18*taub-2*t19*t61*t7*t8*t9-24*t58*t59*t7*t8*t9+18*
     9         t18*t32*t34*t7*t8*t9-1.3333333333333333d+2*t12*t18*t19*t3
     :         6*t9-2.0266666666666666d+2*t18*t19*t35*t7*t8+32*t19*t20*t
     ;         34*t7*t8-96*t20*t31*t32*t7*t8-2.477037037037037d+2*t10*t7
     <         *t8/rhob**8.333333333333334d+0-2*gammabb*t11*t19*t2*t4*t6
     =         1-t10*t15*t4*t61-gammabb*t10*t14*t2*t61-2*t12*t13*t19*t61
     >         -4.5629629629629626d+1*gammabb*t14*t6*t60-4.5629629629629
     ?         626d+1*gammabb*t10*t11*t4*t60-24*gammabb*t11*t2*t4*t58*t5
     @         9-24*t12*t13*t58*t59-6*t15*t32*t4*t58-6*gammabb*t14*t2*t3
     1         2*t58+18*gammabb*t11*t18*t2*t32*t34*t4+6*t15*t18*t19*t34*
     2         t4+16*gammabb*t11*t17*t19*t34*t4-5.866666666666667d+1*gam
     3         mabb*t11*t18*t19*t33*t4-48*gammabb*t11*t17*t31*t32*t4-1.7
     4         77777777777778d+2*t10*t12*t20*t36+18*t12*t13*t18*t32*t34+
     5         6*gammabb*t14*t18*t19*t2*t34+8*gammabb*t10*t14*t17*t34-2.
     6         9333333333333333d+1*gammabb*t10*t14*t18*t33-16*gammabb*t1
     7         4*t17*t19*t31)*wght+1.3333333333333333d+0*t21*t38*wght+4.
     8         0d+0*t22*t37*wght-2.962962962962963d-1*t16*t3*wght)+t10*(
     9         -1.0d+0*t1*t16*t61*wght-4.444444444444444d-1*t16*t18*t38*
     :         wght-1.0d+0*t1*t18*t37*wght-2.6666666666666666d+0*t16*t22
     ;         *t34*wght-2.0d+0*t1*t21*t34*wght-2.6666666666666666d+0*t1
     <         8*t21*t22*wght-2*t18*t39-t23*t34)+t19*(4.0d+0*t1*t16*t18*
     =         t34*wght+2.6666666666666666d+0*t16*t22*t31*wght+2.0d+0*t1
     >         *t21*t31*wght-2*t18*t40)-6.0d+0*t1*t16*t32*t58*wght+Amat3
     ?         (iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t6*(1.0d+0*t1*(-1.493808d-1*t12*t18
     1         *t32*t4*t9*taub-2.48968d-2*t15*t18*t19*t9*taub+1.33333333
     2         33333333d+1*t11*t18*t19*t9*taub-1.6597866666666666d-1*gam
     3         mabb*t11*t19*t42*taub-1.9917439999999997d-1*t12*t19*t20*t
     4         4*taub-7.46904d-2*gammabb*t11*t18*t24*t32*taub-4.97935999
     5         99999993d-2*t10*t15*t20*taub+2.6666666666666666d+1*t10*t1
     6         1*t20*taub-4*gammabb*t19*t34*t7*t9+12*gammabb*t31*t32*t7*
     7         t9-4.481424d-2*gammabb*t11*t31*t4*t59*t9+1.120356d-2*gamm
     8         abb*t11*t32*t34*t4*t9+3.73452d-3*gammabb*t14*t19*t34*t9-1
     9         .120356d-2*gammabb*t14*t31*t32*t9-4.481424d-2*t24*t31*t59
     :         *t7*t8-1.7925696d-1*t18*t32*t42*t7*t8+1.120356d-2*t24*t32
     ;         *t34*t7*t8-2.6888544000000003d-1*t19*t7*t8/rhob**10+6.755
     <         555555555556d+1*gammabb*t10*t35*t7+4.2666666666666664d+1*
     =         gammabb*t18*t19*t20*t7+9.777777777777779d+0*t14*t33*t6-4.
     >         481424d-2*t12*t13*t2*t31*t59-1.2614378666666665d-1*gammab
     ?         b*t11*t19*t35*t4+3.73452d-3*t15*t19*t2*t34*t4-2*t11*t19*t
     @         2*t34*t4-1.8257653333333332d-2*t10*t15*t33*t4+9.777777777
     1         777779d+0*t10*t11*t33*t4-1.120356d-2*t15*t2*t31*t32*t4+6*
     2         t11*t2*t31*t32*t4-1.1950464d-1*gammabb*t11*t18*t20*t32*t4
     3         -1.9917439999999997d-2*t15*t17*t18*t19*t4+1.0666666666666
     4         666d+1*t11*t17*t18*t19*t4-8.298933333333335d-2*t12*t19*t2
     5         4*t36-6.307189333333332d-2*gammabb*t10*t14*t35+1.120356d-
     6         2*t12*t13*t2*t32*t34-t10*t14*t2*t34-3.6515306666666664d-2
     7         *t12*t13*t19*t33-5.975231999999999d-2*t12*t13*t17*t18*t32
     8         +2*t14*t19*t2*t31-3.9834879999999995d-2*gammabb*t14*t18*t
     9         19*t20+5.333333333333333d+0*t10*t14*t17*t18)*wght+2.66666
     :         66666666666d+0*t22*t43*wght+4.444444444444444d-1*t25*t38*
     ;         wght)+t10*(-1.1618506666666667d-2*t16*t55*wght-1.0d+0*t1*
     <         t18*t43*wght+4.9793599999999993d-3*t21*t41*wght-1.0d+0*t1
     =         *t25*t34*wght-1.3333333333333333d+0*t18*t22*t25*wght-t18*
     >         t44+4.9793599999999993d-3*t17*t23)-1.9917439999999997d-2*
     ?         t16*t18*t19*t41*wght-1.120356d-2*t16*t26*t31*t32*wght+2.0
     @         d+0*t1*t19*t25*t31*wght-3.73452d-3*t19*t2*t40-1.86726d-3*
     1         t10*t2*t39+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t6*(1.0d+0*t1*(-6.9733198152d-5*ga
     1         mmabb*t11*t32*t45*taub-1.3946639630399998d-4*t12*t24*t32*
     2         t4*taub-2.3244399384d-5*t15*t19*t24*taub+2.48968d-2*t11*t
     3         19*t24*taub-4*t18*t19*t7*t9-8.367983778239999d-5*t12*t13*
     4         t18*t59*t9-2.0919959445599998d-5*t15*t18*t32*t4*t9+2.2407
     5         12d-2*t11*t18*t32*t4*t9+7.46904d-3*t14*t18*t19*t9-8.36798
     6         3778239999d-5*t18*t45*t59*t7*t8-2.2314623408639997d-4*t32
     7         *t7*t8/rhob**1.1666666666666665d+1+1.1950464d-1*gammabb*t
     8         19*t42*t7+4.481424d-2*gammabb*t18*t24*t32*t7-1.0666666666
     9         666666d+1*t10*t20*t7-8.367983778239999d-5*gammabb*t11*t18
     :         *t24*t4*t59-1.6735967556479997d-4*gammabb*t11*t32*t4*t42-
     ;         5.578655852159999d-5*gammabb*t14*t19*t42-3.71910390143999
     <         96d-5*t15*t19*t20*t4+3.9834879999999995d-2*t11*t19*t20*t4
     =         -2.0919959445599998d-5*gammabb*t14*t18*t24*t32-1.11573117
     >         04319998d-4*t12*t13*t20*t32+1.9917439999999997d-2*t10*t14
     ?         *t20)*wght+1.3333333333333333d+0*t22*t46*wght)-3.71910390
     @         14399996d-5*t16*t19*t62*wght-2.0919959445599998d-5*t16*t1
     1         8*t32*t47*wght-1.0d+0*t1*t10*t18*t46*wght+9.9587199999999
     2         99d-3*t10*t25*t41*wght+7.46904d-3*t18*t19*t25*t26*wght+6.
     3         973319815199999d-6*t19*t23*t9-3.73452d-3*t10*t2*t44+Cmat3
     4         (iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 1.0d+0*t1*t6*(-1.562520138975642d
     1         -7*t59*t7*t8/rhob**1.3333333333333333d+1+1.255197566736d-
     2         4*gammabb*t32*t45*t7-2.240712d-2*t19*t24*t7-1.56252013897
     3         5642d-7*gammabb*t11*t4*t45*t59-1.562520138975642d-7*t12*t
     4         13*t24*t59-3.9063003474391045d-8*gammabb*t14*t32*t45-3.90
     5         63003474391045d-8*t15*t24*t32*t4+6.275987833679999d-5*t11
     6         *t24*t32*t4+2.0919959445599995d-5*t14*t19*t24)*wght+2.091
     7         9959445599995d-5*t19*t25*t47*wght-5.60178d-3*t10*t26*t46*
     8         wght-3.9063003474391045d-8*t16*t32*wght/rhob**6.666666666
     9         666667d+0+Cmat3(iq,D3_GBB_GBB_GBB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = t6*(1.0d+0*t1*(-2.9876159999999996d-
     1         1*t12*t19*t4*t9*taub-7.469039999999999d-2*t10*t15*t9*taub
     2         +8.0d+1*t10*t12*t9*taub-2.987616d-1*t12*t18*t28*t32*t4*ta
     3         ub-1.493808d-1*gammabb*t11*t18*t27*t32*taub-4.97936d-2*t1
     4         5*t18*t19*t28*taub+5.333333333333333d+1*t12*t18*t19*t28*t
     5         aub-2.8216373333333333d-1*gammabb*t11*t19*t24*taub-1.9419
     6         503999999999d-1*gammabb*t11*t18*t32*t4*t9-6.473168d-2*gam
     7         mabb*t14*t18*t19*t9+3.466666666666666d+1*gammabb*t11*t18*
     8         t19*t9-8.962848d-2*t27*t31*t59*t7*t8-4.1826624d-1*t19*t42
     9         *t7*t8+2.240712d-2*t27*t32*t34*t7*t8-3.1369968d-1*t18*t24
     :         *t32*t7*t8+8.88888888888889d+0*t15*t17*t6-8.962848d-2*gam
     ;         mabb*t11*t28*t31*t4*t59-8.962848d-2*t12*t13*t3*t31*t59+2.
     <         240712d-2*gammabb*t11*t28*t32*t34*t4+7.46904d-3*t15*t19*t
     =         3*t34*t4-8.0d+0*t12*t19*t3*t34*t4-2.240712d-2*t15*t3*t31*
     >         t32*t4+2.4d+1*t12*t3*t31*t32*t4-1.7261781333333331d-1*gam
     ?         mabb*t11*t19*t20*t4-2.48968d-2*t15*t18*t19*t2*t4+2.666666
     @         6666666666d+1*t12*t18*t19*t2*t4-1.6597866666666666d-2*t10
     1         *t15*t17*t4+1.7777777777777778d+1*t10*t12*t17*t4-1.659786
     2         666666667d-1*t12*t19*t27*t36+2.240712d-2*t12*t13*t3*t32*t
     3         34-2.0d+0*t10*t15*t3*t34+7.46904d-3*gammabb*t14*t19*t28*t
     4         34-4.0d+0*gammabb*t11*t19*t28*t34-2.240712d-2*gammabb*t14
     5         *t28*t31*t32+1.2d+1*gammabb*t11*t28*t31*t32-7.46904d-2*t1
     6         2*t13*t18*t2*t32+4.0d+0*t15*t19*t3*t31-8.630890666666666d
     7         -2*gammabb*t10*t14*t20+4.622222222222223d+1*gammabb*t10*t
     8         11*t20+6.666666666666667d+0*t10*t15*t18*t2-3.319573333333
     9         333d-2*t12*t13*t17*t19)*wght+2.6666666666666666d+0*t22*t4
     :         8*wght+4.444444444444444d-1*t29*t38*wght)+t10*(-1.0d+0*t1
     ;         *t18*t48*wght-8.298933333333333d-3*t16*t41*wght-1.0d+0*t1
     <         *t29*t34*wght-1.3333333333333333d+0*t18*t22*t29*wght+6.22
     =         42d-3*t21*t26*wght-t18*t49+6.2242d-3*t2*t23)-2.240712d-2*
     >         t16*t30*t31*t32*wght+2.0d+0*t1*t19*t29*t31*wght-2.48968d-
     ?         2*t16*t18*t19*t26*wght-7.46904d-3*t19*t3*t40-3.73452d-3*t
     @         10*t3*t39+Mmat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = t6*(1.0d+0*t1*(-1.3946639630399998d
     1         -4*gammabb*t11*t32*t50*taub-2.7893279260799997d-4*t12*t27
     2         *t32*t4*taub-4.6488798768d-5*t15*t19*t27*taub+4.97936d-2*
     3         t12*t19*t27*taub+2.48968d-2*t11*t19*t27*taub-6.0435438398
     4         39999d-5*t15*t19*t4*t9+6.473168d-2*t12*t19*t4*t9+3.236584
     5         d-2*t11*t19*t4*t9-1.8130631519519996d-4*t12*t13*t32*t9+1.
     6         618292d-2*t10*t15*t9+1.618292d-2*t10*t14*t9-8.66666666666
     7         6666d+0*t10*t11*t9-1.673596755648d-4*t18*t50*t59*t7*t8-4.
     8         044525492816d-4*t32*t45*t7*t8+4.481424d-2*gammabb*t18*t27
     9         *t32*t7+1.0456655999999999d-1*gammabb*t19*t24*t7-1.673596
     :         755648d-4*gammabb*t11*t18*t27*t4*t59-1.673596755648d-4*t1
     ;         2*t13*t18*t28*t59-4.1839918891199995d-5*t15*t18*t28*t32*t
     <         4+4.481424d-2*t12*t18*t28*t32*t4+2.240712d-2*t11*t18*t28*
     =         t32*t4-2.9287943223839996d-4*gammabb*t11*t24*t32*t4-4.183
     >         9918891199995d-5*gammabb*t14*t18*t27*t32+2.240712d-2*gamm
     ?         abb*t11*t18*t27*t32+7.46904d-3*t15*t18*t19*t28+7.46904d-3
     @         *t14*t18*t19*t28-4.0d+0*t11*t18*t19*t28-9.76264774128d-5*
     1         gammabb*t14*t19*t24+5.228328d-2*gammabb*t11*t19*t24)*wght
     2         +1.3333333333333333d+0*t22*t51*wght)-4.1839918891199995d-
     3         5*t16*t18*t32*t52*wght-1.0d+0*t1*t10*t18*t51*wght-6.04354
     4         3839839999d-5*t16*t19*t47*wght+4.9793599999999993d-3*t10*
     5         t29*t41*wght+7.46904d-3*t18*t19*t25*t30*wght+3.73452d-3*t
     6         18*t19*t26*t29*wght+6.2242d-3*t10*t25*t26*wght-1.86726d-3
     7         *t10*t2*t49-3.73452d-3*t10*t3*t44+1.3946639630399998d-5*t
     8         19*t23*t28+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = 1.0d+0*t1*t6*(-3.125040277951284d-
     1         7*t59*t7*t8/rhob**1.2333333333333334d+1+1.673596755648d-4
     2         *gammabb*t32*t50*t7-1.4938079999999998d-2*t19*t27*t7-3.12
     3         5040277951284d-7*gammabb*t11*t4*t50*t59-3.125040277951284
     4         d-7*t12*t13*t27*t59-7.812600694878209d-8*gammabb*t14*t32*
     5         t50+4.1839918891199995d-5*gammabb*t11*t32*t50-7.812600694
     6         878209d-8*t15*t27*t32*t4+8.367983778239999d-5*t12*t27*t32
     7         *t4+8.367983778239999d-5*t11*t27*t32*t4+1.394663963039999
     8         8d-5*t15*t19*t27+2.7893279260799997d-5*t14*t19*t27-1.4938
     9         079999999998d-2*t11*t19*t27)*wght-7.812600694878209d-8*t1
     :         6*t32*t60*wght+2.7893279260799997d-5*t19*t25*t52*wght-3.7
     ;         3452d-3*t10*t26*t51*wght+6.973319815199999d-6*t19*t29*t47
     <         *wght-3.73452d-3*t10*t30*t46*wght+Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = t6*(1.0d+0*t1*(-5.578655852159999d-4
     1         *t12*t32*t4*t54*taub-9.2977597536d-5*t15*t19*t54*taub+1.9
     2         91744d-1*t12*t19*t54*taub-2.7893279260799997d-4*gammabb*t
     3         11*t32*t53*taub-3.347193511296d-4*t18*t53*t59*t7*t8-7.252
     4         252607808d-4*t32*t50*t7*t8-3.347193511296d-4*t12*t13*t18*
     5         t55*t59-3.347193511296d-4*gammabb*t11*t18*t4*t54*t59-8.36
     6         7983778239999d-5*t15*t18*t32*t4*t55+1.7925696d-1*t12*t18*
     7         t32*t4*t55+2.9876159999999996d-2*t15*t18*t19*t55-1.6d+1*t
     8         12*t18*t19*t55-8.367983778239999d-5*gammabb*t14*t18*t32*t
     9         54+8.962848d-2*gammabb*t11*t18*t32*t54-5.020790266944d-4*
     :         gammabb*t11*t27*t32*t4-9.2977597536d-5*t15*t19*t28*t4+1.9
     ;         91744d-1*t12*t19*t28*t4-2.7893279260799997d-4*t12*t13*t28
     <         *t32+4.97936d-2*t10*t15*t28-2.6666666666666666d+1*t10*t12
     =         *t28-1.673596755648d-4*gammabb*t14*t19*t27+1.792569600000
     >         0002d-1*gammabb*t11*t19*t27)*wght+1.3333333333333333d+0*t
     ?         22*t56*wght)-8.367983778239999d-5*t16*t18*t32*t57*wght-1.
     @         0d+0*t1*t10*t18*t56*wght-9.2977597536d-5*t16*t19*t52*wght
     1         +1.4938079999999998d-2*t18*t19*t29*t30*wght+1.24484d-2*t1
     2         0*t26*t29*wght+2.7893279260799997d-5*t19*t23*t55-7.46904d
     3         -3*t10*t3*t49+Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
            Mmat3(iq,D3_GBB_TB_TB) = 1.0d+0*t1*t6*(-6.250080555902569d-7
     1         *t59*t7*t8/rhob**1.1333333333333334d+1+1.673596755648d-4*
     2         gammabb*t32*t53*t7-6.250080555902569d-7*t12*t13*t54*t59-6
     3         .250080555902569d-7*gammabb*t11*t4*t53*t59-1.562520138975
     4         6418d-7*t15*t32*t4*t54+3.347193511296d-4*t12*t32*t4*t54+8
     5         .367983778239999d-5*t11*t32*t4*t54+5.578655852159999d-5*t
     6         15*t19*t54+2.7893279260799997d-5*t14*t19*t54-2.9876159999
     7         999996d-2*t12*t19*t54-2.9876159999999996d-2*t11*t19*t54-1
     8         .5625201389756418d-7*gammabb*t14*t32*t53+1.673596755648d-
     9         4*gammabb*t11*t32*t53)*wght+2.7893279260799997d-5*t19*t25
     :         *t57*wght-1.86726d-3*t10*t26*t56*wght+2.7893279260799997d
     ;         -5*t19*t29*t52*wght-7.46904d-3*t10*t30*t51*wght-1.5625201
     <         389756418d-7*t16*t32*t33*wght+Mmat3(iq,D3_GBB_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = 1.0d+0*t1*t6*(-1.2500161111805136d-6
     1         *t59*t7*t8/rhob**1.0333333333333333d+1-1.2500161111805136
     2         d-6*gammabb*t11*t4*t59*t63-3.1250402779512837d-7*gammabb*
     3         t14*t32*t63+5.020790266944d-4*gammabb*t11*t32*t63-1.25001
     4         61111805136d-6*t12*t13*t59*t62-3.1250402779512837d-7*t15*
     5         t32*t4*t62+1.0041580533887999d-3*t12*t32*t4*t62+1.6735967
     6         556479997d-4*t15*t19*t62-1.7925696d-1*t12*t19*t62)*wght+8
     7         .367983778239998d-5*t19*t29*t57*wght-1.120356d-2*t10*t30*
     8         t56*wght-3.1250402779512837d-7*t16*t17*t32*wght+Mmat3(iq,
     9         D3_TB_TB_TB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
