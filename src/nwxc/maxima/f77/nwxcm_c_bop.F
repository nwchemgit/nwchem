C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_c_bop.F 
C> The nwxcm_c_bop functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_c_bop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}\,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{\it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110, 10664 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
      subroutine nwxcm_c_bop(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = gammaaa**5.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = asinh(t2*t3)
            t5 = 2.52d-2*t2*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 8.4d-3*gammaaa*t6*t7+1.8610514726982d+0
            t9 = rhoa**3.333333333333333d-1
            t10 = t1*t8*t9
            t11 = t10**2.0d+0
            t12 = 1/t11
            t13 = 5.642d-1*t1*t8*t9+2.5d-1*t11+3.183d-1
            t14 = 1/t13
            t15 = 7.607d-1*t1*t8*t9+5.764d-1
            t16 = rhoa**2
            t17 = 1/t10**3.0d+0
            t18 = rhoa**6.666666666666666d-1
            t19 = 1/t18
            t20 = 1/rhoa**3.6666666666666664d+0
            t21 = 1/(gammaaa*t7+1)**5.0d-1
            t22 = rhoa**2.3333333333333334d+0
            t23 = 1/t5**2
            t24 = -8.4d-3*gammaaa*t23*(-3.3599999999999997d-2*t2*t4/t22-
     1         3.3599999999999997d-2*gammaaa*t20*t21)*t7-2.2399999999999
     2         998d-2*gammaaa*t20*t6
            t25 = t1*t24*t9+3.333333333333333d-1*t1*t19*t8
            t26 = 1/t13**2
            t27 = 8.4d-3*t6*t7-8.4d-3*gammaaa*t23*t7*(1.26d-2*t21*t7+1.2
     1         6d-2*t3*t4/t2)
            t28 = 4.0d+0*t12*t15*t16*t26*(5.642d-1*t1*t27*t9+5.0d-1*t1**
     1         2*t18*t27*t8)*wght+8.0d+0*t1*t14*t15*t17*t22*t27*wght-3.0
     2         428d+0*t1*t12*t14*t22*t27*wght
            fnc(iq) = fnc(iq)-4.0d+0*t12*t14*t15*t16*wght
            Amat(iq,D1_RA) = t12*t14*(-4.0d+0*t16*(7.607d-1*t1*t24*t9+2.
     1         5356666666666666d-1*t1*t19*t8)*wght-8.0d+0*rhoa*t15*wght)
     2         +4.0d+0*t12*t15*t16*t26*(5.0d-1*t1*t25*t8*t9+5.642d-1*t1*
     3         t24*t9+1.8806666666666666d-1*t1*t19*t8)*wght+8.0d+0*t14*t
     4         15*t16*t17*t25*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t28+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t28+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = gammaaa**5.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = asinh(t2*t3)
            t5 = 2.52d-2*t2*t3*t4+1.0d+0
            t6 = 1/t5
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 8.4d-3*gammaaa*t6*t7+1.8610514726982d+0
            t9 = rhoa**3.333333333333333d-1
            t10 = gammabb**5.0d-1
            t11 = 1/rhob**1.3333333333333333d+0
            t12 = asinh(t10*t11)
            t13 = 2.52d-2*t10*t11*t12+1.0d+0
            t14 = 1/t13
            t15 = 1/rhob**2.6666666666666666d+0
            t16 = 8.4d-3*gammabb*t14*t15+1.8610514726982d+0
            t17 = rhob**3.333333333333333d-1
            t18 = t8*t9+t16*t17
            t19 = 1/t18
            t20 = t1*t16*t17*t19*t8*t9
            t21 = t20**2.0d+0
            t22 = 1/t21
            t23 = 1.1284000000000001d+0*t1*t16*t17*t19*t8*t9+t21+3.183d-
     1         1
            t24 = 1/t23
            t25 = 1.5214d+0*t1*t16*t17*t19*t8*t9+5.764d-1
            t26 = 1/t20**3.0d+0
            t27 = rhoa**6.666666666666666d-1
            t28 = 1/t27
            t29 = 1/rhoa**3.6666666666666664d+0
            t30 = 1/(gammaaa*t7+1)**5.0d-1
            t31 = 1/t5**2
            t32 = -8.4d-3*gammaaa*t31*(-3.3599999999999997d-2*t2*t4/rhoa
     1         **2.3333333333333334d+0-3.3599999999999997d-2*gammaaa*t29
     2         *t30)*t7-2.2399999999999998d-2*gammaaa*t29*t6
            t33 = t32*t9+3.333333333333333d-1*t28*t8
            t34 = 1/t18**2
            t35 = -t1*t16*t17*t33*t34*t8*t9+t1*t16*t17*t19*t32*t9+3.3333
     1         33333333333d-1*t1*t16*t17*t19*t28*t8
            t36 = 1/t23**2
            t37 = rhob**6.666666666666666d-1
            t38 = 1/t37
            t39 = 1/rhob**3.6666666666666664d+0
            t40 = 1/(gammabb*t15+1)**5.0d-1
            t41 = 1/t13**2
            t42 = -8.4d-3*gammabb*t15*(-3.3599999999999997d-2*gammabb*t3
     1         9*t40-3.3599999999999997d-2*t10*t12/rhob**2.3333333333333
     2         334d+0)*t41-2.2399999999999998d-2*gammabb*t14*t39
            t43 = t17*t42+3.333333333333333d-1*t16*t38
            t44 = -t1*t16*t17*t34*t43*t8*t9+t1*t17*t19*t42*t8*t9+3.33333
     1         3333333333d-1*t1*t16*t19*t38*t8*t9
            t45 = 8.4d-3*t6*t7-8.4d-3*gammaaa*t31*t7*(1.26d-2*t30*t7+1.2
     1         6d-2*t3*t4/t2)
            t46 = t1*t16*t17*t19*t45*t9-t1*t16*t17*t27*t34*t45*t8
            t47 = 8.4d-3*t14*t15-8.4d-3*gammabb*t15*(1.26d-2*t15*t40+1.2
     1         6d-2*t11*t12/t10)*t41
            t48 = t1*t17*t19*t47*t8*t9-t1*t16*t34*t37*t47*t8*t9
            fnc(iq) = fnc(iq)-1.0d+0*rhoa*rhob*t22*t24*t25*wght
            Amat(iq,D1_RA) = t22*t24*(-1.0d+0*rhoa*rhob*(-1.5214d+0*t1*t
     1         16*t17*t33*t34*t8*t9+1.5214d+0*t1*t16*t17*t19*t32*t9+5.07
     2         1333333333333d-1*t1*t16*t17*t19*t28*t8)*wght-1.0d+0*rhob*
     3         t25*wght)+1.0d+0*rhoa*rhob*t22*t25*t36*(2.0d+0*t1*t16*t17
     4         *t19*t35*t8*t9-1.1284000000000001d+0*t1*t16*t17*t33*t34*t
     5         8*t9+1.1284000000000001d+0*t1*t16*t17*t19*t32*t9+3.761333
     6         333333333d-1*t1*t16*t17*t19*t28*t8)*wght+2.0d+0*rhoa*rhob
     7         *t24*t25*t26*t35*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = t22*t24*(-1.0d+0*rhoa*rhob*(-1.5214d+0*t1*t
     1         16*t17*t34*t43*t8*t9+1.5214d+0*t1*t17*t19*t42*t8*t9+5.071
     2         333333333333d-1*t1*t16*t19*t38*t8*t9)*wght-1.0d+0*rhoa*t2
     3         5*wght)+1.0d+0*rhoa*rhob*t22*t25*t36*(2.0d+0*t1*t16*t17*t
     4         19*t44*t8*t9-1.1284000000000001d+0*t1*t16*t17*t34*t43*t8*
     5         t9+1.1284000000000001d+0*t1*t17*t19*t42*t8*t9+3.761333333
     6         333333d-1*t1*t16*t19*t38*t8*t9)*wght+2.0d+0*rhoa*rhob*t24
     7         *t25*t26*t44*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = 1.0d+0*rhoa*rhob*t22*t25*t36*(2.0d+0*t1*t1
     1         6*t17*t19*t46*t8*t9+1.1284000000000001d+0*t1*t16*t17*t19*
     2         t45*t9-1.1284000000000001d+0*t1*t16*t17*t27*t34*t45*t8)*w
     3         ght-1.0d+0*rhoa*rhob*t22*t24*(1.5214d+0*t1*t16*t17*t19*t4
     4         5*t9-1.5214d+0*t1*t16*t17*t27*t34*t45*t8)*wght+2.0d+0*rho
     5         a*rhob*t24*t25*t26*t46*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*rhoa*rhob*t22*t25*t36*(2.0d+0*t1*t1
     1         6*t17*t19*t48*t8*t9-1.1284000000000001d+0*t1*t16*t34*t37*
     2         t47*t8*t9+1.1284000000000001d+0*t1*t17*t19*t47*t8*t9)*wgh
     3         t-1.0d+0*rhoa*rhob*t22*t24*(1.5214d+0*t1*t17*t19*t47*t8*t
     4         9-1.5214d+0*t1*t16*t34*t37*t47*t8*t9)*wght+2.0d+0*rhoa*rh
     5         ob*t24*t25*t26*t48*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_c_bop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}\,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{\it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110, 10664 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
      subroutine nwxcm_c_bop_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_c_bop functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& \rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_3} &=& \rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}\,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}+1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{\it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284\,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110, 10664 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
      subroutine nwxcm_c_bop_d3(param,tol_rho,ipol,nq,wght, 
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
