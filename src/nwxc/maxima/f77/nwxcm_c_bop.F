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
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}
C>    \,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{
C>    \it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}
C>    \,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}
C>    +0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284
C>    \,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{
C>    \it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110
C>    , 10664 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_bop(param,tol_rho,ipol,nq,wght, 
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
            t17 = 1/t13**2
            t18 = rhoa**6.666666666666666d-1
            t19 = 1/t18
            t20 = 1/rhoa**3.6666666666666664d+0
            t21 = 1/(gammaaa*t7+1)**5.0d-1
            t22 = rhoa**2.3333333333333334d+0
            t23 = 1/t5**2
            t24 = -8.4d-3*gammaaa*t23*(-3.3599999999999997d-2*t2*t4/t22-
     1         3.3599999999999997d-2*gammaaa*t20*t21)*t7-2.2399999999999
     2         998d-2*gammaaa*t20*t6
            t25 = t24*t9+3.333333333333333d-1*t19*t8
            t26 = 5.0d-1*t1*t24*t9+1.6666666666666666d-1*t1*t19*t8-2.5d-
     1         1*t1*t25
            t27 = 1/t10**3.0d+0
            t28 = 8.4d-3*t6*t7-8.4d-3*gammaaa*t23*t7*(1.26d-2*t21*t7+1.2
     1         6d-2*t3*t4/t2)
            fnc(iq) = fnc(iq)-4.0d+0*t12*t14*t15*t16*wght
            Amat(iq,D1_RA) = 4.0d+0*t12*t14*(-1.0d+0*t16*(7.607d-1*t1*t2
     1         4*t9+2.5356666666666666d-1*t1*t19*t8-3.8035d-1*t1*t25)*wg
     2         ht-1.0d+0*rhoa*t15*wght)+4.0d+0*t12*t15*t16*t17*(1.0d+0*t
     3         1*t26*t8*t9+5.642d-1*t1*t24*t9+1.8806666666666666d-1*t1*t
     4         19*t8-2.821d-1*t1*t25)*wght+1.6d+1*t14*t15*t16*t26*t27*wg
     5         ht+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 4.0d+0*t12*t15*t16*t17*(2.821d-1*t1*t28*t9
     1         +2.5d-1*t1**2*t18*t28*t8)*wght+4.0d+0*t1*t14*t15*t22*t27*
     2         t28*wght-1.5214d+0*t1*t12*t14*t22*t28*wght+Cmat(iq,D1_GAA
     3         )
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
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
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}
C>    \,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{
C>    \it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}
C>    \,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}
C>    +0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284
C>    \,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{
C>    \it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110
C>    , 10664 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_bop_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = gammaaa**5.0d-1
            t3 = rhoa**1.3333333333333333d+0
            t4 = 1/t3
            t5 = asinh(t2*t4)
            t6 = 2.52d-2*t2*t4*t5+1.0d+0
            t7 = 1/t6
            t8 = rhoa**2.6666666666666666d+0
            t9 = 1/t8
            t10 = 8.4d-3*gammaaa*t7*t9+1.8610514726982d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = t1*t10*t11
            t13 = t12**2.0d+0
            t14 = 1/t13
            t15 = 2.5d-1*t13+5.642d-1*t1*t10*t11+3.183d-1
            t16 = 1/t15
            t17 = 7.607d-1*t1*t10*t11+5.764d-1
            t18 = rhoa**2
            t19 = 1/t15**2
            t20 = rhoa**6.666666666666666d-1
            t21 = 1/t20
            t22 = 1/rhoa**3.6666666666666664d+0
            t23 = (gammaaa*t9+1)**5.0d-1
            t24 = 1/t23
            t25 = rhoa**2.3333333333333334d+0
            t26 = 1/t25
            t27 = -3.3599999999999997d-2*t2*t26*t5-3.3599999999999997d-2
     1         *gammaaa*t22*t24
            t28 = 1/t6**2
            t29 = -8.4d-3*gammaaa*t27*t28*t9-2.2399999999999998d-2*gamma
     1         aa*t22*t7
            t30 = t11*t29+3.333333333333333d-1*t10*t21
            t31 = -2.5d-1*t1*t30+5.0d-1*t1*t11*t29+1.6666666666666666d-1
     1         *t1*t10*t21
            t32 = 1.0d+0*t1*t10*t11*t31-2.821d-1*t1*t30+5.642d-1*t1*t11*
     1         t29+1.8806666666666666d-1*t1*t10*t21
            t33 = 1/t12**3.0d+0
            t34 = -3.8035d-1*t1*t30+7.607d-1*t1*t11*t29+2.53566666666666
     1         66d-1*t1*t10*t21
            t35 = -1.0d+0*t18*t34*wght-1.0d+0*rhoa*t17*wght
            t36 = 1/t2
            t37 = 1.26d-2*t24*t9+1.26d-2*t36*t4*t5
            t38 = 8.4d-3*t7*t9-8.4d-3*gammaaa*t28*t37*t9
            t39 = t1**2
            t40 = 2.5d-1*t10*t20*t38*t39+2.821d-1*t1*t11*t38
            t41 = 1/t15**3
            t42 = -8.0d+0*t14*t17*t18*t32**2*t41*wght
            t43 = 1/t12**4.0d+0
            t44 = -9.6d+1*t16*t17*t18*t31**2*t43*wght
            t45 = -3.2d+1*t17*t18*t19*t31*t32*t33*wght
            t46 = -2.0d+0*rhoa*t34*wght
            t47 = 1/rhoa**1.6666666666666669d+0
            t48 = 1/rhoa
            t49 = -2.5356666666666666d-1*t1*t30*t48
            t50 = 5.071333333333333d-1*t1*t21*t29
            t51 = 1/t10
            t52 = t30**2
            t53 = 1/t11
            t54 = 3.8035d-1*t1*t51*t52*t53
            t55 = 1/rhoa**4.666666666666667d+0
            t56 = 1/t6**3
            t57 = 1/t23**3
            t58 = -8.4d-3*gammaaa*t28*(-4.4799999999999995d-2*gammaaa**2
     1         *t57/rhoa**7.333333333333333d+0+1.6799999999999998d-1*gam
     2         maaa*t24*t55+7.84d-2*t2*t5/rhoa**3.3333333333333337d+0)*t
     3         9+1.6799999999999998d-2*gammaaa*t27**2*t56*t9+8.213333333
     4         333332d-2*gammaaa*t55*t7+4.4799999999999995d-2*gammaaa*t2
     5         2*t27*t28
            t59 = t11*t58-2.222222222222222d-1*t10*t47+6.666666666666666
     1         d-1*t21*t29
            t60 = -7.607d-1*t1*t29*t30*t51
            t61 = 1.0d+0*rhoa*t17*t32*wght
            t62 = 1.0d+0*t18*t32*t34*wght
            t63 = -1.8806666666666666d-1*t1*t30*t48
            t64 = 3.761333333333333d-1*t1*t21*t29
            t65 = 3.333333333333333d-1*t1*t10*t21*t31
            t66 = 2.821d-1*t1*t51*t52*t53
            t67 = 1.0d+0*t1*t11*t29*t31
            t68 = -1.6666666666666666d-1*t1*t30*t48
            t69 = 3.333333333333333d-1*t1*t21*t29
            t70 = 2.5d-1*t1*t51*t52*t53
            t71 = -5.0d-1*t1*t29*t30*t51
            t72 = t71+t70+t69+t68-2.5d-1*t1*t59+5.0d-1*t1*t11*t58-1.1111
     1         11111111111d-1*t1*t10*t47
            t73 = -5.642d-1*t1*t29*t30*t51
            t74 = -5.0d-1*t1*t30*t31
            t75 = -t32*t35
            t76 = 2.0d+0*rhoa*t17*t31*wght
            t77 = 2.0d+0*t18*t31*t34*wght
            t78 = -2.0d+0*t31*t35
            t79 = t29**2
            t80 = 5.0d-1*t1*t11*t51*t79+t71+t70+t69+t68+5.55555555555555
     1         6d-2*t1*t10*t47
            t81 = -8.0d+0*t14*t17*t18*t32*t40*t41*wght
            t82 = -2.4d+1*t1*t16*t17*t25*t31*t38*t43*wght
            t83 = -3.8035d-1*t1*t3*t38*wght
            t84 = 1.2678333333333333d-1*t1*t21*t38
            t85 = -8.4d-3*gammaaa*t28*(1.6799999999999998d-2*gammaaa*t57
     1         /rhoa**6.333333333333333d+0-1.6799999999999998d-2*t26*t36
     2         *t5-5.04d-2*t22*t24)*t9+1.6799999999999998d-2*gammaaa*t27
     3         *t37*t56*t9-8.4d-3*t27*t28*t9-2.2399999999999998d-2*t22*t
     4         7+2.2399999999999998d-2*gammaaa*t22*t28*t37
            t86 = t11*t85+3.333333333333333d-1*t21*t38
            t87 = 8.0d+0*t19*t33*(-2.0d+0*t17*t18*t31*t40*wght-5.0d-1*t1
     1         *t17*t25*t32*t38*wght)
            t88 = 8.333333333333333d-2*t1*t21*t38
            t89 = t88-2.5d-1*t1*t86+5.0d-1*t1*t11*t85-2.5d-1*t1*t11*t29*
     1         t38*t51
            t90 = 7.607d-1*t1*t25*t31*t38*wght
            t91 = -5.0d-1*t1*t11*t35*t38
            t92 = 9.403333333333332d-2*t1*t21*t38
            t93 = 5.0d-1*t1*t11*t31*t38
            t94 = 3.8035d-1*t1*t25*t32*t38*wght
            t95 = -t35*t40
            t96 = t88+2.5d-1*t1*t11*t29*t38*t51
            t97 = t38**2
            t98 = -8.4d-3*gammaaa*t28*t9*(6.3d-3*t24*t9/gammaaa-6.3d-3*t
     1         57/rhoa**5.333333333333333d+0-6.3d-3*t4*t5/t2**3)+1.67999
     2         99999999998d-2*gammaaa*t37**2*t56*t9-1.6799999999999998d-
     3         2*t28*t37*t9
            t99 = -8.0d+0*t14*t17*t18*t40**2*t41*wght
            t100 = -8.0d+0*t1*t17*t19*t25*t33*t38*t40*wght
            t101 = -6.0d+0*t16*t17*t39*t43*t8*t97*wght
            t102 = 2.5d-1*t1*t11*t98-2.5d-1*t1*t11*t51*t97
            t103 = 7.607d-1*t1*t25*t38*t40*wght
            t104 = 3.8035d-1*t39*t8*t97*wght
            fnc(iq) = fnc(iq)-4.0d+0*t14*t16*t17*t18*wght
            Amat(iq,D1_RA) = 1.6d+1*t16*t17*t18*t31*t33*wght+4.0d+0*t14*
     1         t17*t18*t19*t32*wght+4.0d+0*t14*t16*t35+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 4.0d+0*t14*t17*t18*t19*t40*wght+4.0d+0*t1*
     1         t16*t17*t25*t33*t38*wght-1.5214d+0*t1*t14*t16*t25*t38*wgh
     2         t+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(t74+t73
     1         +1.0d+0*t1*t10*t11*t72+t67+t66+t65+t64+t63-2.821d-1*t1*t5
     2         9+5.642d-1*t1*t11*t58-1.2537777777777775d-1*t1*t10*t47)*w
     3         ght+t75+t62+t61)+8.0d+0*t16*t33*(2.0d+0*t17*t18*t72*wght+
     4         t78+t77+t76)+4.0d+0*t14*t16*(t46-1.0d+0*t18*(t60-3.8035d-
     5         1*t1*t59+7.607d-1*t1*t11*t58+t54+t50+t49-1.69044444444444
     6         44d-1*t1*t10*t47)*wght)+t45+t44+t42+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(1.0d+0*
     1         t1*t10*t11*t80+5.642d-1*t1*t11*t51*t79+t74+t73+t67+t66+t6
     2         5+t64+t63+6.268888888888888d-2*t1*t10*t47)*wght+t75+t62+t
     3         61)+8.0d+0*t16*t33*(2.0d+0*t17*t18*t80*wght+t78+t77+t76)+
     4         4.0d+0*t14*t16*(-1.0d+0*t18*(7.607d-1*t1*t11*t51*t79+t60+
     5         t54+t50+t49+8.452222222222222d-2*t1*t10*t47)*wght-1.0d+0*
     6         t17*wght+t46)+t45+t44+t42+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(t93+t9
     1         2+1.0d+0*t1*t10*t11*t89-2.821d-1*t1*t86+5.642d-1*t1*t11*t
     2         85-2.821d-1*t1*t11*t29*t38*t51)*wght+t95+t94)+8.0d+0*t16*
     3         t33*(2.0d+0*t17*t18*t89*wght+t91+t90)+4.0d+0*t14*t16*(t83
     4         -1.0d+0*t18*(-3.8035d-1*t1*t86+7.607d-1*t1*t11*t85+t84-3.
     5         8035d-1*t1*t11*t29*t38*t51)*wght)+t87+t82+t81+Cmat2(iq,D2
     6         _RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(1.0d+0
     1         *t1*t10*t11*t96+t93+t92+2.821d-1*t1*t11*t29*t38*t51)*wght
     2         +t95+t94)+8.0d+0*t16*t33*(2.0d+0*t17*t18*t96*wght+t91+t90
     3         )+4.0d+0*t14*t16*(t83-1.0d+0*t18*(t84+3.8035d-1*t1*t11*t2
     4         9*t38*t51)*wght)+t87+t82+t81+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(2.821
     1         d-1*t1*t11*t98-2.821d-1*t1*t11*t51*t97+1.25d-1*t20*t39*t9
     2         7+1.0d+0*t1*t10*t102*t11)*wght+t103)+8.0d+0*t16*t33*(2.0d
     3         +0*t102*t17*t18*wght+t104)-4.0d+0*t14*t16*t18*(3.8035d-1*
     4         t1*t11*t98-3.8035d-1*t1*t11*t51*t97)*wght+t99+t101+t100+C
     5         mat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(2.821
     1         d-1*t1*t11*t51*t97+3.75d-1*t20*t39*t97)*wght+t103)+8.0d+0
     2         *t16*t33*(5.0d-1*t1*t17*t25*t51*t97*wght+t104)-1.5214d+0*
     3         t1*t14*t16*t25*t51*t97*wght+t99+t101+t100+Cmat2(iq,D2_GAA
     4         _GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
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
            t30 = (gammaaa*t7+1)**5.0d-1
            t31 = 1/t30
            t32 = 1/rhoa**2.3333333333333334d+0
            t33 = -3.3599999999999997d-2*t2*t32*t4-3.3599999999999997d-2
     1         *gammaaa*t29*t31
            t34 = 1/t5**2
            t35 = -8.4d-3*gammaaa*t33*t34*t7-2.2399999999999998d-2*gamma
     1         aa*t29*t6
            t36 = t35*t9+3.333333333333333d-1*t28*t8
            t37 = 1/t18**2
            t38 = -t1*t16*t17*t36*t37*t8*t9+t1*t16*t17*t19*t35*t9+3.3333
     1         33333333333d-1*t1*t16*t17*t19*t28*t8
            t39 = 1/t23**2
            t40 = 2.0d+0*t1*t16*t17*t19*t38*t8*t9-1.1284000000000001d+0*
     1         t1*t16*t17*t36*t37*t8*t9+1.1284000000000001d+0*t1*t16*t17
     2         *t19*t35*t9+3.761333333333333d-1*t1*t16*t17*t19*t28*t8
            t41 = -1.5214d+0*t1*t16*t17*t36*t37*t8*t9+1.5214d+0*t1*t16*t
     1         17*t19*t35*t9+5.071333333333333d-1*t1*t16*t17*t19*t28*t8
            t42 = -1.0d+0*rhoa*rhob*t41*wght-1.0d+0*rhob*t25*wght
            t43 = rhob**6.666666666666666d-1
            t44 = 1/t43
            t45 = 1/rhob**3.6666666666666664d+0
            t46 = (gammabb*t15+1)**5.0d-1
            t47 = 1/t46
            t48 = 1/rhob**2.3333333333333334d+0
            t49 = -3.3599999999999997d-2*t10*t12*t48-3.3599999999999997d
     1         -2*gammabb*t45*t47
            t50 = 1/t13**2
            t51 = -8.4d-3*gammabb*t15*t49*t50-2.2399999999999998d-2*gamm
     1         abb*t14*t45
            t52 = t17*t51+3.333333333333333d-1*t16*t44
            t53 = -t1*t16*t17*t37*t52*t8*t9+t1*t17*t19*t51*t8*t9+3.33333
     1         3333333333d-1*t1*t16*t19*t44*t8*t9
            t54 = 2.0d+0*t1*t16*t17*t19*t53*t8*t9-1.1284000000000001d+0*
     1         t1*t16*t17*t37*t52*t8*t9+1.1284000000000001d+0*t1*t17*t19
     2         *t51*t8*t9+3.761333333333333d-1*t1*t16*t19*t44*t8*t9
            t55 = -1.5214d+0*t1*t16*t17*t37*t52*t8*t9+1.5214d+0*t1*t17*t
     1         19*t51*t8*t9+5.071333333333333d-1*t1*t16*t19*t44*t8*t9
            t56 = -1.0d+0*rhoa*rhob*t55*wght-1.0d+0*rhoa*t25*wght
            t57 = 1/t2
            t58 = 1.26d-2*t31*t7+1.26d-2*t3*t4*t57
            t59 = 8.4d-3*t6*t7-8.4d-3*gammaaa*t34*t58*t7
            t60 = 1.5214d+0*t1*t16*t17*t19*t59*t9-1.5214d+0*t1*t16*t17*t
     1         27*t37*t59*t8
            t61 = t1*t16*t17*t19*t59*t9-t1*t16*t17*t27*t37*t59*t8
            t62 = 2.0d+0*t1*t16*t17*t19*t61*t8*t9+1.1284000000000001d+0*
     1         t1*t16*t17*t19*t59*t9-1.1284000000000001d+0*t1*t16*t17*t2
     2         7*t37*t59*t8
            t63 = 1/t10
            t64 = 1.26d-2*t11*t12*t63+1.26d-2*t15*t47
            t65 = 8.4d-3*t14*t15-8.4d-3*gammabb*t15*t50*t64
            t66 = 1.5214d+0*t1*t17*t19*t65*t8*t9-1.5214d+0*t1*t16*t37*t4
     1         3*t65*t8*t9
            t67 = t1*t17*t19*t65*t8*t9-t1*t16*t37*t43*t65*t8*t9
            t68 = 2.0d+0*t1*t16*t17*t19*t67*t8*t9-1.1284000000000001d+0*
     1         t1*t16*t37*t43*t65*t8*t9+1.1284000000000001d+0*t1*t17*t19
     2         *t65*t8*t9
            t69 = 1/t23**3
            t70 = 1/t20**4.0d+0
            t71 = t36**2
            t72 = 1/t18**3
            t73 = 1/rhoa**1.6666666666666669d+0
            t74 = 1/rhoa**4.666666666666667d+0
            t75 = 1/t5**3
            t76 = 1/t30**3
            t77 = -8.4d-3*gammaaa*t34*t7*(-4.4799999999999995d-2*gammaaa
     1         **2*t76/rhoa**7.333333333333333d+0+1.6799999999999998d-1*
     2         gammaaa*t31*t74+7.84d-2*t2*t4/rhoa**3.3333333333333337d+0
     3         )+1.6799999999999998d-2*gammaaa*t33**2*t7*t75+8.213333333
     4         333332d-2*gammaaa*t6*t74+4.4799999999999995d-2*gammaaa*t2
     5         9*t33*t34
            t78 = t77*t9-2.222222222222222d-1*t73*t8+6.666666666666666d-
     1         1*t28*t35
            t79 = -t1*t16*t17*t37*t78*t8*t9+2*t1*t16*t17*t71*t72*t8*t9+t
     1         1*t16*t17*t19*t77*t9-2*t1*t16*t17*t35*t36*t37*t9-2.222222
     2         222222222d-1*t1*t16*t17*t19*t73*t8-6.666666666666666d-1*t
     3         1*t16*t17*t28*t36*t37*t8+6.666666666666666d-1*t1*t16*t17*
     4         t19*t28*t35
            t80 = 2*t1*t16*t17*t36*t52*t72*t8*t9-t1*t17*t36*t37*t51*t8*t
     1         9-3.333333333333333d-1*t1*t16*t36*t37*t44*t8*t9-t1*t16*t1
     2         7*t35*t37*t52*t9+t1*t17*t19*t35*t51*t9+3.333333333333333d
     3         -1*t1*t16*t19*t35*t44*t9-3.333333333333333d-1*t1*t16*t17*
     4         t28*t37*t52*t8+3.333333333333333d-1*t1*t17*t19*t28*t51*t8
     5         +1.111111111111111d-1*t1*t16*t19*t28*t44*t8
            t81 = 1/rhob**1.6666666666666669d+0
            t82 = t52**2
            t83 = 1/rhob**4.666666666666667d+0
            t84 = 1/t13**3
            t85 = 1/t46**3
            t86 = -8.4d-3*gammabb*t15*t50*(-4.4799999999999995d-2*gammab
     1         b**2*t85/rhob**7.333333333333333d+0+1.6799999999999998d-1
     2         *gammabb*t47*t83+7.84d-2*t10*t12/rhob**3.3333333333333337
     3         d+0)+1.6799999999999998d-2*gammabb*t15*t49**2*t84+8.21333
     4         3333333332d-2*gammabb*t14*t83+4.4799999999999995d-2*gamma
     5         bb*t45*t49*t50
            t87 = t17*t86-2.222222222222222d-1*t16*t81+6.666666666666666
     1         d-1*t44*t51
            t88 = -t1*t16*t17*t37*t8*t87*t9+t1*t17*t19*t8*t86*t9+2*t1*t1
     1         6*t17*t72*t8*t82*t9-2.222222222222222d-1*t1*t16*t19*t8*t8
     2         1*t9-2*t1*t17*t37*t51*t52*t8*t9-6.666666666666666d-1*t1*t
     3         16*t37*t44*t52*t8*t9+6.666666666666666d-1*t1*t19*t44*t51*
     4         t8*t9
            t89 = 1/t9
            t90 = -8.4d-3*gammaaa*t34*t7*(1.6799999999999998d-2*gammaaa*
     1         t76/rhoa**6.333333333333333d+0-1.6799999999999998d-2*t32*
     2         t4*t57-5.04d-2*t29*t31)+1.6799999999999998d-2*gammaaa*t33
     3         *t58*t7*t75-8.4d-3*t33*t34*t7-2.2399999999999998d-2*t29*t
     4         6+2.2399999999999998d-2*gammaaa*t29*t34*t58
            t91 = t9*t90+3.333333333333333d-1*t28*t59
            t92 = -t1*t16*t17*t37*t8*t9*t91+t1*t16*t17*t19*t9*t90-t1*t16
     1         *t17*t36*t37*t59*t9-3.333333333333333d-1*t1*t16*t17*t37*t
     2         59*t8*t89+2*t1*t16*t17*t27*t36*t59*t72*t8-t1*t16*t17*t27*
     3         t35*t37*t59+3.333333333333333d-1*t1*t16*t17*t19*t28*t59
            t93 = 2*t1*t16*t36*t43*t65*t72*t8*t9-t1*t17*t36*t37*t65*t8*t
     1         9-t1*t16*t35*t37*t43*t65*t9+t1*t17*t19*t35*t65*t9-3.33333
     2         3333333333d-1*t1*t16*t28*t37*t43*t65*t8+3.333333333333333
     3         d-1*t1*t17*t19*t28*t65*t8
            t94 = -t1*t16*t17*t37*t52*t59*t9+t1*t17*t19*t51*t59*t9+3.333
     1         333333333333d-1*t1*t16*t19*t44*t59*t9+2*t1*t16*t17*t27*t5
     2         2*t59*t72*t8-t1*t17*t27*t37*t51*t59*t8-3.333333333333333d
     3         -1*t1*t16*t27*t37*t44*t59*t8
            t95 = 1/t17
            t96 = -8.4d-3*gammabb*t15*t50*(1.6799999999999998d-2*gammabb
     1         *t85/rhob**6.333333333333333d+0-1.6799999999999998d-2*t12
     2         *t48*t63-5.04d-2*t45*t47)+1.6799999999999998d-2*gammabb*t
     3         15*t49*t64*t84+2.2399999999999998d-2*gammabb*t45*t50*t64-
     4         8.4d-3*t15*t49*t50-2.2399999999999998d-2*t14*t45
            t97 = t17*t96+3.333333333333333d-1*t44*t65
            t98 = -t1*t16*t17*t37*t8*t9*t97+t1*t17*t19*t8*t9*t96-3.33333
     1         3333333333d-1*t1*t16*t37*t65*t8*t9*t95+2*t1*t16*t43*t52*t
     2         65*t72*t8*t9-t1*t17*t37*t52*t65*t8*t9-t1*t37*t43*t51*t65*
     3         t8*t9+3.333333333333333d-1*t1*t19*t44*t65*t8*t9
            t99 = t59**2
            t100 = -8.4d-3*gammaaa*t34*t7*(-6.3d-3*t76/rhoa**5.333333333
     1         333333d+0+6.3d-3*t31*t7/gammaaa-6.3d-3*t3*t4/t2**3)+1.679
     2         9999999999998d-2*gammaaa*t58**2*t7*t75-1.6799999999999998
     3         d-2*t34*t58*t7
            t101 = 2*rhoa*t1*t16*t17*t72*t8*t99-2*t1*t16*t17*t27*t37*t99
     1         +t1*t100*t16*t17*t19*t9-t1*t100*t16*t17*t27*t37*t8
            t102 = -t1*t16*t37*t43*t59*t65*t9+t1*t17*t19*t59*t65*t9+2*t1
     1         *t16*t27*t43*t59*t65*t72*t8-t1*t17*t27*t37*t59*t65*t8
            t103 = -8.4d-3*gammabb*t15*t50*(-6.3d-3*t85/rhob**5.33333333
     1         3333333d+0+6.3d-3*t15*t47/gammabb-6.3d-3*t11*t12/t10**3)+
     2         1.6799999999999998d-2*gammabb*t15*t64**2*t84-1.6799999999
     3         999998d-2*t15*t50*t64
            t104 = t65**2
            t105 = 2*rhob*t1*t104*t16*t72*t8*t9-t1*t103*t16*t37*t43*t8*t
     1         9-2*t1*t104*t37*t43*t8*t9+t1*t103*t17*t19*t8*t9
            fnc(iq) = fnc(iq)-1.0d+0*rhoa*rhob*t22*t24*t25*wght
            Amat(iq,D1_RA) = 1.0d+0*rhoa*rhob*t22*t25*t39*t40*wght+2.0d+
     1         0*rhoa*rhob*t24*t25*t26*t38*wght+t22*t24*t42+Amat(iq,D1_R
     2         A)
            Amat(iq,D1_RB) = 1.0d+0*rhoa*rhob*t22*t25*t39*t54*wght+2.0d+
     1         0*rhoa*rhob*t24*t25*t26*t53*wght+t22*t24*t56+Amat(iq,D1_R
     2         B)
            Cmat(iq,D1_GAA) = 1.0d+0*rhoa*rhob*t22*t25*t39*t62*wght+2.0d
     1         +0*rhoa*rhob*t24*t25*t26*t61*wght-1.0d+0*rhoa*rhob*t22*t2
     2         4*t60*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*rhoa*rhob*t22*t25*t39*t68*wght+2.0d
     1         +0*rhoa*rhob*t24*t25*t26*t67*wght-1.0d+0*rhoa*rhob*t22*t2
     2         4*t66*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*t
     1         1*t16*t17*t19*t79*t8*t9-1.1284000000000001d+0*t1*t16*t17*
     2         t37*t78*t8*t9+2.2568000000000002d+0*t1*t16*t17*t71*t72*t8
     3         *t9-2.0d+0*t1*t16*t17*t36*t37*t38*t8*t9+1.128400000000000
     4         1d+0*t1*t16*t17*t19*t77*t9+2.0d+0*t1*t16*t17*t19*t35*t38*
     5         t9-2.2568000000000002d+0*t1*t16*t17*t35*t36*t37*t9-2.5075
     6         55555555555d-1*t1*t16*t17*t19*t73*t8+6.666666666666666d-1
     7         *t1*t16*t17*t19*t28*t38*t8-7.522666666666666d-1*t1*t16*t1
     8         7*t28*t36*t37*t8+7.522666666666666d-1*t1*t16*t17*t19*t28*
     9         t35)*wght+1.0d+0*rhoa*rhob*t40*t41*wght+1.0d+0*rhob*t25*t
     :         40*wght-t40*t42)+t22*t24*(-1.0d+0*rhoa*rhob*(-1.5214d+0*t
     ;         1*t16*t17*t37*t78*t8*t9+3.0428d+0*t1*t16*t17*t71*t72*t8*t
     <         9+1.5214d+0*t1*t16*t17*t19*t77*t9-3.0428d+0*t1*t16*t17*t3
     =         5*t36*t37*t9-3.380888888888889d-1*t1*t16*t17*t19*t73*t8-1
     >         .0142666666666666d+0*t1*t16*t17*t28*t36*t37*t8+1.01426666
     ?         66666666d+0*t1*t16*t17*t19*t28*t35)*wght-2.0d+0*rhob*t41*
     @         wght)+t24*t26*(2.0d+0*rhoa*rhob*t25*t79*wght+2.0d+0*rhoa*
     1         rhob*t38*t41*wght+2.0d+0*rhob*t25*t38*wght-2.0d+0*t38*t42
     2         )-6.0d+0*rhoa*rhob*t24*t25*t38**2*t70*wght-2.0d+0*rhoa*rh
     3         ob*t22*t25*t40**2*t69*wght-4.0d+0*rhoa*rhob*t25*t26*t38*t
     4         39*t40*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*t
     1         1*t16*t17*t19*t8*t80*t9+2.2568000000000002d+0*t1*t16*t17*
     2         t36*t52*t72*t8*t9-2.0d+0*t1*t16*t17*t37*t38*t52*t8*t9+2.0
     3         d+0*t1*t17*t19*t38*t51*t8*t9-1.1284000000000001d+0*t1*t17
     4         *t36*t37*t51*t8*t9+6.666666666666666d-1*t1*t16*t19*t38*t4
     5         4*t8*t9-3.761333333333333d-1*t1*t16*t36*t37*t44*t8*t9-1.1
     6         284000000000001d+0*t1*t16*t17*t35*t37*t52*t9+1.1284000000
     7         000001d+0*t1*t17*t19*t35*t51*t9+3.761333333333333d-1*t1*t
     8         16*t19*t35*t44*t9-3.761333333333333d-1*t1*t16*t17*t28*t37
     9         *t52*t8+3.761333333333333d-1*t1*t17*t19*t28*t51*t8+1.2537
     :         777777777775d-1*t1*t16*t19*t28*t44*t8)*wght+1.0d+0*rhoa*r
     ;         hob*t40*t55*wght+1.0d+0*rhoa*t25*t40*wght-t42*t54)+t22*t2
     <         4*(-1.0d+0*rhoa*rhob*(3.0428d+0*t1*t16*t17*t36*t52*t72*t8
     =         *t9-1.5214d+0*t1*t17*t36*t37*t51*t8*t9-5.071333333333333d
     >         -1*t1*t16*t36*t37*t44*t8*t9-1.5214d+0*t1*t16*t17*t35*t37*
     ?         t52*t9+1.5214d+0*t1*t17*t19*t35*t51*t9+5.071333333333333d
     @         -1*t1*t16*t19*t35*t44*t9-5.071333333333333d-1*t1*t16*t17*
     1         t28*t37*t52*t8+5.071333333333333d-1*t1*t17*t19*t28*t51*t8
     2         +1.6904444444444444d-1*t1*t16*t19*t28*t44*t8)*wght-1.0d+0
     3         *rhob*t55*wght-1.0d+0*rhoa*t41*wght-1.0d+0*t25*wght)+t24*
     4         t26*(2.0d+0*rhoa*rhob*t25*t80*wght+2.0d+0*rhoa*rhob*t38*t
     5         55*wght+2.0d+0*rhoa*t25*t38*wght-2.0d+0*t42*t53)+t26*t39*
     6         (-2.0d+0*rhoa*rhob*t25*t38*t54*wght-2.0d+0*rhoa*rhob*t25*
     7         t40*t53*wght)-6.0d+0*rhoa*rhob*t24*t25*t38*t53*t70*wght-2
     8         .0d+0*rhoa*rhob*t22*t25*t40*t54*t69*wght+Amat2(iq,D2_RA_R
     9         B)
            Amat2(iq,D2_RB_RB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*t
     1         1*t16*t17*t19*t8*t88*t9-1.1284000000000001d+0*t1*t16*t17*
     2         t37*t8*t87*t9+1.1284000000000001d+0*t1*t17*t19*t8*t86*t9+
     3         2.2568000000000002d+0*t1*t16*t17*t72*t8*t82*t9-2.50755555
     4         5555555d-1*t1*t16*t19*t8*t81*t9-2.0d+0*t1*t16*t17*t37*t52
     5         *t53*t8*t9+2.0d+0*t1*t17*t19*t51*t53*t8*t9+6.666666666666
     6         666d-1*t1*t16*t19*t44*t53*t8*t9-2.2568000000000002d+0*t1*
     7         t17*t37*t51*t52*t8*t9-7.522666666666666d-1*t1*t16*t37*t44
     8         *t52*t8*t9+7.522666666666666d-1*t1*t19*t44*t51*t8*t9)*wgh
     9         t+1.0d+0*rhoa*rhob*t54*t55*wght+1.0d+0*rhoa*t25*t54*wght-
     :         t54*t56)+t22*t24*(-1.0d+0*rhoa*rhob*(-1.5214d+0*t1*t16*t1
     ;         7*t37*t8*t87*t9+1.5214d+0*t1*t17*t19*t8*t86*t9+3.0428d+0*
     <         t1*t16*t17*t72*t8*t82*t9-3.380888888888889d-1*t1*t16*t19*
     =         t8*t81*t9-3.0428d+0*t1*t17*t37*t51*t52*t8*t9-1.0142666666
     >         666666d+0*t1*t16*t37*t44*t52*t8*t9+1.0142666666666666d+0*
     ?         t1*t19*t44*t51*t8*t9)*wght-2.0d+0*rhoa*t55*wght)+t24*t26*
     @         (2.0d+0*rhoa*rhob*t25*t88*wght+2.0d+0*rhoa*rhob*t53*t55*w
     1         ght+2.0d+0*rhoa*t25*t53*wght-2.0d+0*t53*t56)-6.0d+0*rhoa*
     2         rhob*t24*t25*t53**2*t70*wght-2.0d+0*rhoa*rhob*t22*t25*t54
     3         **2*t69*wght-4.0d+0*rhoa*rhob*t25*t26*t39*t53*t54*wght+Am
     4         at2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*
     1         t1*t16*t17*t19*t8*t9*t92-1.1284000000000001d+0*t1*t16*t17
     2         *t37*t8*t9*t91+1.1284000000000001d+0*t1*t16*t17*t19*t9*t9
     3         0+2.0d+0*t1*t16*t17*t19*t38*t59*t9-1.1284000000000001d+0*
     4         t1*t16*t17*t36*t37*t59*t9-3.761333333333333d-1*t1*t16*t17
     5         *t37*t59*t8*t89+2.2568000000000002d+0*t1*t16*t17*t27*t36*
     6         t59*t72*t8-2.0d+0*t1*t16*t17*t27*t37*t38*t59*t8-1.1284000
     7         000000001d+0*t1*t16*t17*t27*t35*t37*t59+3.761333333333333
     8         d-1*t1*t16*t17*t19*t28*t59)*wght+1.0d+0*rhoa*rhob*t40*t60
     9         *wght-t42*t62)+t24*t26*(2.0d+0*rhoa*rhob*t25*t92*wght+2.0
     :         d+0*rhoa*rhob*t38*t60*wght-2.0d+0*t42*t61)+t22*t24*(-1.0d
     ;         +0*rhoa*rhob*(-1.5214d+0*t1*t16*t17*t37*t8*t9*t91+1.5214d
     <         +0*t1*t16*t17*t19*t9*t90-1.5214d+0*t1*t16*t17*t36*t37*t59
     =         *t9-5.071333333333333d-1*t1*t16*t17*t37*t59*t8*t89+3.0428
     >         d+0*t1*t16*t17*t27*t36*t59*t72*t8-1.5214d+0*t1*t16*t17*t2
     ?         7*t35*t37*t59+5.071333333333333d-1*t1*t16*t17*t19*t28*t59
     @         )*wght-1.0d+0*rhob*t60*wght)+t26*t39*(-2.0d+0*rhoa*rhob*t
     1         25*t38*t62*wght-2.0d+0*rhoa*rhob*t25*t40*t61*wght)-6.0d+0
     2         *rhoa*rhob*t24*t25*t38*t61*t70*wght-2.0d+0*rhoa*rhob*t22*
     3         t25*t40*t62*t69*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*
     1         t1*t16*t17*t19*t8*t9*t93+2.2568000000000002d+0*t1*t16*t36
     2         *t43*t65*t72*t8*t9-2.0d+0*t1*t16*t37*t38*t43*t65*t8*t9+2.
     3         0d+0*t1*t17*t19*t38*t65*t8*t9-1.1284000000000001d+0*t1*t1
     4         7*t36*t37*t65*t8*t9-1.1284000000000001d+0*t1*t16*t35*t37*
     5         t43*t65*t9+1.1284000000000001d+0*t1*t17*t19*t35*t65*t9-3.
     6         761333333333333d-1*t1*t16*t28*t37*t43*t65*t8+3.7613333333
     7         33333d-1*t1*t17*t19*t28*t65*t8)*wght+1.0d+0*rhoa*rhob*t40
     8         *t66*wght-t42*t68)+t24*t26*(2.0d+0*rhoa*rhob*t25*t93*wght
     9         +2.0d+0*rhoa*rhob*t38*t66*wght-2.0d+0*t42*t67)+t22*t24*(-
     :         1.0d+0*rhoa*rhob*(3.0428d+0*t1*t16*t36*t43*t65*t72*t8*t9-
     ;         1.5214d+0*t1*t17*t36*t37*t65*t8*t9-1.5214d+0*t1*t16*t35*t
     <         37*t43*t65*t9+1.5214d+0*t1*t17*t19*t35*t65*t9-5.071333333
     =         333333d-1*t1*t16*t28*t37*t43*t65*t8+5.071333333333333d-1*
     >         t1*t17*t19*t28*t65*t8)*wght-1.0d+0*rhob*t66*wght)+t26*t39
     ?         *(-2.0d+0*rhoa*rhob*t25*t38*t68*wght-2.0d+0*rhoa*rhob*t25
     @         *t40*t67*wght)-6.0d+0*rhoa*rhob*t24*t25*t38*t67*t70*wght-
     1         2.0d+0*rhoa*rhob*t22*t25*t40*t68*t69*wght+Cmat2(iq,D2_RA_
     2         GBB)
            Cmat2(iq,D2_RB_GAA) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*
     1         t1*t16*t17*t19*t8*t9*t94+2.0d+0*t1*t16*t17*t19*t53*t59*t9
     2         -1.1284000000000001d+0*t1*t16*t17*t37*t52*t59*t9+1.128400
     3         0000000001d+0*t1*t17*t19*t51*t59*t9+3.761333333333333d-1*
     4         t1*t16*t19*t44*t59*t9+2.2568000000000002d+0*t1*t16*t17*t2
     5         7*t52*t59*t72*t8-2.0d+0*t1*t16*t17*t27*t37*t53*t59*t8-1.1
     6         284000000000001d+0*t1*t17*t27*t37*t51*t59*t8-3.7613333333
     7         33333d-1*t1*t16*t27*t37*t44*t59*t8)*wght+1.0d+0*rhoa*rhob
     8         *t54*t60*wght-t56*t62)+t24*t26*(2.0d+0*rhoa*rhob*t25*t94*
     9         wght+2.0d+0*rhoa*rhob*t53*t60*wght-2.0d+0*t56*t61)+t22*t2
     :         4*(-1.0d+0*rhoa*rhob*(-1.5214d+0*t1*t16*t17*t37*t52*t59*t
     ;         9+1.5214d+0*t1*t17*t19*t51*t59*t9+5.071333333333333d-1*t1
     <         *t16*t19*t44*t59*t9+3.0428d+0*t1*t16*t17*t27*t52*t59*t72*
     =         t8-1.5214d+0*t1*t17*t27*t37*t51*t59*t8-5.071333333333333d
     >         -1*t1*t16*t27*t37*t44*t59*t8)*wght-1.0d+0*rhoa*t60*wght)+
     ?         t26*t39*(-2.0d+0*rhoa*rhob*t25*t53*t62*wght-2.0d+0*rhoa*r
     @         hob*t25*t54*t61*wght)-6.0d+0*rhoa*rhob*t24*t25*t53*t61*t7
     1         0*wght-2.0d+0*rhoa*rhob*t22*t25*t54*t62*t69*wght+Cmat2(iq
     2         ,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.0d+0*
     1         t1*t16*t17*t19*t8*t9*t98-1.1284000000000001d+0*t1*t16*t17
     2         *t37*t8*t9*t97+1.1284000000000001d+0*t1*t17*t19*t8*t9*t96
     3         -3.761333333333333d-1*t1*t16*t37*t65*t8*t9*t95+2.25680000
     4         00000002d+0*t1*t16*t43*t52*t65*t72*t8*t9-2.0d+0*t1*t16*t3
     5         7*t43*t53*t65*t8*t9+2.0d+0*t1*t17*t19*t53*t65*t8*t9-1.128
     6         4000000000001d+0*t1*t17*t37*t52*t65*t8*t9-1.1284000000000
     7         001d+0*t1*t37*t43*t51*t65*t8*t9+3.761333333333333d-1*t1*t
     8         19*t44*t65*t8*t9)*wght+1.0d+0*rhoa*rhob*t54*t66*wght-t56*
     9         t68)+t24*t26*(2.0d+0*rhoa*rhob*t25*t98*wght+2.0d+0*rhoa*r
     :         hob*t53*t66*wght-2.0d+0*t56*t67)+t22*t24*(-1.0d+0*rhoa*rh
     ;         ob*(-1.5214d+0*t1*t16*t17*t37*t8*t9*t97+1.5214d+0*t1*t17*
     <         t19*t8*t9*t96-5.071333333333333d-1*t1*t16*t37*t65*t8*t9*t
     =         95+3.0428d+0*t1*t16*t43*t52*t65*t72*t8*t9-1.5214d+0*t1*t1
     >         7*t37*t52*t65*t8*t9-1.5214d+0*t1*t37*t43*t51*t65*t8*t9+5.
     ?         071333333333333d-1*t1*t19*t44*t65*t8*t9)*wght-1.0d+0*rhoa
     @         *t66*wght)+t26*t39*(-2.0d+0*rhoa*rhob*t25*t53*t68*wght-2.
     1         0d+0*rhoa*rhob*t25*t54*t67*wght)-6.0d+0*rhoa*rhob*t24*t25
     2         *t53*t67*t70*wght-2.0d+0*rhoa*rhob*t22*t25*t54*t68*t69*wg
     3         ht+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.2568
     1         000000000002d+0*rhoa*t1*t16*t17*t72*t8*t99-2.256800000000
     2         0002d+0*t1*t16*t17*t27*t37*t99+2.0d+0*t1*t101*t16*t17*t19
     3         *t8*t9+2.0d+0*t1*t16*t17*t19*t59*t61*t9+1.128400000000000
     4         1d+0*t1*t100*t16*t17*t19*t9-2.0d+0*t1*t16*t17*t27*t37*t59
     5         *t61*t8-1.1284000000000001d+0*t1*t100*t16*t17*t27*t37*t8)
     6         *wght+2.0d+0*rhoa*rhob*t60*t62*wght)+t24*t26*(4.0d+0*rhoa
     7         *rhob*t60*t61*wght+2.0d+0*rhoa*rhob*t101*t25*wght)-1.0d+0
     8         *rhoa*rhob*t22*t24*(3.0428d+0*rhoa*t1*t16*t17*t72*t8*t99-
     9         3.0428d+0*t1*t16*t17*t27*t37*t99+1.5214d+0*t1*t100*t16*t1
     :         7*t19*t9-1.5214d+0*t1*t100*t16*t17*t27*t37*t8)*wght-6.0d+
     ;         0*rhoa*rhob*t24*t25*t61**2*t70*wght-2.0d+0*rhoa*rhob*t22*
     <         t25*t62**2*t69*wght-4.0d+0*rhoa*rhob*t25*t26*t39*t61*t62*
     =         wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(-2.0d+
     1         0*t1*t16*t37*t43*t61*t65*t8*t9+2.0d+0*t1*t17*t19*t61*t65*
     2         t8*t9+2.0d+0*t1*t102*t16*t17*t19*t8*t9-1.1284000000000001
     3         d+0*t1*t16*t37*t43*t59*t65*t9+1.1284000000000001d+0*t1*t1
     4         7*t19*t59*t65*t9+2.2568000000000002d+0*t1*t16*t27*t43*t59
     5         *t65*t72*t8-1.1284000000000001d+0*t1*t17*t27*t37*t59*t65*
     6         t8)*wght+1.0d+0*rhoa*rhob*t60*t68*wght+1.0d+0*rhoa*rhob*t
     7         62*t66*wght)+t26*t39*(-2.0d+0*rhoa*rhob*t25*t61*t68*wght-
     8         2.0d+0*rhoa*rhob*t25*t62*t67*wght)+t24*t26*(2.0d+0*rhoa*r
     9         hob*t60*t67*wght+2.0d+0*rhoa*rhob*t61*t66*wght+2.0d+0*rho
     :         a*rhob*t102*t25*wght)-1.0d+0*rhoa*rhob*t22*t24*(-1.5214d+
     ;         0*t1*t16*t37*t43*t59*t65*t9+1.5214d+0*t1*t17*t19*t59*t65*
     <         t9+3.0428d+0*t1*t16*t27*t43*t59*t65*t72*t8-1.5214d+0*t1*t
     =         17*t27*t37*t59*t65*t8)*wght-6.0d+0*rhoa*rhob*t24*t25*t61*
     >         t67*t70*wght-2.0d+0*rhoa*rhob*t22*t25*t62*t68*t69*wght+Cm
     ?         at2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t22*t39*(1.0d+0*rhoa*rhob*t25*(2.2568
     1         000000000002d+0*rhob*t1*t104*t16*t72*t8*t9-2.0d+0*t1*t16*
     2         t37*t43*t65*t67*t8*t9+2.0d+0*t1*t17*t19*t65*t67*t8*t9-1.1
     3         284000000000001d+0*t1*t103*t16*t37*t43*t8*t9-2.2568000000
     4         000002d+0*t1*t104*t37*t43*t8*t9+2.0d+0*t1*t105*t16*t17*t1
     5         9*t8*t9+1.1284000000000001d+0*t1*t103*t17*t19*t8*t9)*wght
     6         +2.0d+0*rhoa*rhob*t66*t68*wght)+t24*t26*(4.0d+0*rhoa*rhob
     7         *t66*t67*wght+2.0d+0*rhoa*rhob*t105*t25*wght)-1.0d+0*rhoa
     8         *rhob*t22*t24*(3.0428d+0*rhob*t1*t104*t16*t72*t8*t9-1.521
     9         4d+0*t1*t103*t16*t37*t43*t8*t9-3.0428d+0*t1*t104*t37*t43*
     :         t8*t9+1.5214d+0*t1*t103*t17*t19*t8*t9)*wght-6.0d+0*rhoa*r
     ;         hob*t24*t25*t67**2*t70*wght-2.0d+0*rhoa*rhob*t22*t25*t68*
     <         *2*t69*wght-4.0d+0*rhoa*rhob*t25*t26*t39*t67*t68*wght+Cma
     =         t2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
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
C>   {\it t_6} &=& {{0.0084\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}\,\left(0.0252\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_8} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_9} &=& {{0.0084\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}\,\left(0.0252\,{\it t_7}
C>    \,{\rm asinh}\; \left({\it t_7}\,{\it t_8}\right)\,{\it t_8}
C>    +1.0\right)}}+1.8610514726982\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_3}\,{\it t_9}+{\it t_2}\,{
C>    \it t_6}}}\\\\ 
C>   {\it t_{11}} &=& \left({\it t_1}\,{\it t_2}\,{\it t_3}\,{
C>    \it t_6}\,{\it t_9}\,{\it t_{10}}\right)^{2.0}\\\\ 
C>   f &=& -{{1.0\,\rho_\alpha\,\rho_\beta\,\left(1.5214\,{\it t_1}
C>    \,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{\it t_{10}}
C>    +0.5764\right)}\over{{\it t_{11}}\,\left({\it t_{11}}+1.1284
C>    \,{\it t_1}\,{\it t_2}\,{\it t_3}\,{\it t_6}\,{\it t_9}\,{
C>    \it t_{10}}+0.3183\right)}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& 0.0\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] T Tsuneda, T Suzumura, K Hirao, J.Chem.Phys. 110
C>    , 10664 (1999)  , DOI: 
C> <a href="https://doi.org/10.1063/1.479012 "> 
C> 10.1063/1.479012 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_c_bop_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = gammaaa**5.0d-1
            t3 = rhoa**1.3333333333333333d+0
            t4 = 1/t3
            t5 = asinh(t2*t4)
            t6 = 2.52d-2*t2*t4*t5+1.0d+0
            t7 = 1/t6
            t8 = rhoa**2.6666666666666666d+0
            t9 = 1/t8
            t10 = 8.4d-3*gammaaa*t7*t9+1.8610514726982d+0
            t11 = rhoa**3.333333333333333d-1
            t12 = t1*t10*t11
            t13 = t12**2.0d+0
            t14 = 1/t13
            t15 = 2.5d-1*t13+5.642d-1*t1*t10*t11+3.183d-1
            t16 = 1/t15
            t17 = 7.607d-1*t1*t10*t11+5.764d-1
            t18 = rhoa**2
            t19 = 1/t15**2
            t20 = rhoa**6.666666666666666d-1
            t21 = 1/t20
            t22 = 1/rhoa**3.6666666666666664d+0
            t23 = (gammaaa*t9+1)**5.0d-1
            t24 = 1/t23
            t25 = rhoa**2.3333333333333334d+0
            t26 = 1/t25
            t27 = -3.3599999999999997d-2*t2*t26*t5-3.3599999999999997d-2
     1         *gammaaa*t22*t24
            t28 = 1/t6**2
            t29 = -8.4d-3*gammaaa*t27*t28*t9-2.2399999999999998d-2*gamma
     1         aa*t22*t7
            t30 = t11*t29+3.333333333333333d-1*t10*t21
            t31 = -2.5d-1*t1*t30+5.0d-1*t1*t11*t29+1.6666666666666666d-1
     1         *t1*t10*t21
            t32 = 1.0d+0*t1*t10*t11*t31-2.821d-1*t1*t30+5.642d-1*t1*t11*
     1         t29+1.8806666666666666d-1*t1*t10*t21
            t33 = 1/t12**3.0d+0
            t34 = -3.8035d-1*t1*t30+7.607d-1*t1*t11*t29+2.53566666666666
     1         66d-1*t1*t10*t21
            t35 = -1.0d+0*t18*t34*wght-1.0d+0*rhoa*t17*wght
            t36 = 1/t2
            t37 = 1.26d-2*t24*t9+1.26d-2*t36*t4*t5
            t38 = 8.4d-3*t7*t9-8.4d-3*gammaaa*t28*t37*t9
            t39 = t1**2
            t40 = 2.5d-1*t10*t20*t38*t39+2.821d-1*t1*t11*t38
            t41 = 1/t15**3
            t42 = t32**2
            t43 = -8.0d+0*t14*t17*t18*t41*t42*wght
            t44 = 1/t12**4.0d+0
            t45 = t31**2
            t46 = -9.6d+1*t16*t17*t18*t44*t45*wght
            t47 = -3.2d+1*t17*t18*t19*t31*t32*t33*wght
            t48 = -2.0d+0*rhoa*t34*wght
            t49 = 1/rhoa**1.6666666666666669d+0
            t50 = 1/rhoa
            t51 = -2.5356666666666666d-1*t1*t30*t50
            t52 = 5.071333333333333d-1*t1*t21*t29
            t53 = 1/t10
            t54 = t30**2
            t55 = 1/t11
            t56 = 3.8035d-1*t1*t53*t54*t55
            t57 = 1/rhoa**4.666666666666667d+0
            t58 = t27**2
            t59 = 1/t6**3
            t60 = gammaaa**2
            t61 = 1/t23**3
            t62 = 1/rhoa**7.333333333333333d+0
            t63 = 1/rhoa**3.3333333333333337d+0
            t64 = 7.84d-2*t2*t5*t63-4.4799999999999995d-2*t60*t61*t62+1.
     1         6799999999999998d-1*gammaaa*t24*t57
            t65 = -8.4d-3*gammaaa*t28*t64*t9+1.6799999999999998d-2*gamma
     1         aa*t58*t59*t9+8.213333333333332d-2*gammaaa*t57*t7+4.47999
     2         99999999995d-2*gammaaa*t22*t27*t28
            t66 = t11*t65-2.222222222222222d-1*t10*t49+6.666666666666666
     1         d-1*t21*t29
            t67 = -7.607d-1*t1*t29*t30*t53
            t68 = t67-3.8035d-1*t1*t66+7.607d-1*t1*t11*t65+t56+t52+t51-1
     1         .6904444444444444d-1*t1*t10*t49
            t69 = t48-1.0d+0*t18*t68*wght
            t70 = 1.0d+0*rhoa*t17*t32*wght
            t71 = 1.0d+0*t18*t32*t34*wght
            t72 = -1.8806666666666666d-1*t1*t30*t50
            t73 = 3.761333333333333d-1*t1*t21*t29
            t74 = 3.333333333333333d-1*t1*t10*t21*t31
            t75 = 2.821d-1*t1*t53*t54*t55
            t76 = 1.0d+0*t1*t11*t29*t31
            t77 = -1.6666666666666666d-1*t1*t30*t50
            t78 = 3.333333333333333d-1*t1*t21*t29
            t79 = 2.5d-1*t1*t53*t54*t55
            t80 = -5.0d-1*t1*t29*t30*t53
            t81 = t80+t79+t78+t77-2.5d-1*t1*t66+5.0d-1*t1*t11*t65-1.1111
     1         11111111111d-1*t1*t10*t49
            t82 = -5.642d-1*t1*t29*t30*t53
            t83 = -5.0d-1*t1*t30*t31
            t84 = t83+t82+1.0d+0*t1*t10*t11*t81+t76+t75+t74+t73+t72-2.82
     1         1d-1*t1*t66+5.642d-1*t1*t11*t65-1.2537777777777775d-1*t1*
     2         t10*t49
            t85 = -t32*t35
            t86 = 1.0d+0*t17*t18*t84*wght+t85+t71+t70
            t87 = 2.0d+0*rhoa*t17*t31*wght
            t88 = 2.0d+0*t18*t31*t34*wght
            t89 = -2.0d+0*t31*t35
            t90 = 2.0d+0*t17*t18*t81*wght+t89+t88+t87
            t91 = t29**2
            t92 = 7.607d-1*t1*t11*t53*t91+t67+t56+t52+t51+8.452222222222
     1         222d-2*t1*t10*t49
            t93 = -1.0d+0*t18*t92*wght-1.0d+0*t17*wght+t48
            t94 = 5.0d-1*t1*t11*t53*t91+t80+t79+t78+t77+5.55555555555555
     1         6d-2*t1*t10*t49
            t95 = 1.0d+0*t1*t10*t11*t94+5.642d-1*t1*t11*t53*t91+t83+t82+
     1         t76+t75+t74+t73+t72+6.268888888888888d-2*t1*t10*t49
            t96 = 1.0d+0*t17*t18*t95*wght+t85+t71+t70
            t97 = 2.0d+0*t17*t18*t94*wght+t89+t88+t87
            t98 = -8.0d+0*t14*t17*t18*t32*t40*t41*wght
            t99 = -2.4d+1*t1*t16*t17*t25*t31*t38*t44*wght
            t100 = -3.8035d-1*t1*t3*t38*wght
            t101 = 1.2678333333333333d-1*t1*t21*t38
            t102 = 1/rhoa**6.333333333333333d+0
            t103 = 1.6799999999999998d-2*gammaaa*t102*t61-1.679999999999
     1         9998d-2*t26*t36*t5-5.04d-2*t22*t24
            t104 = 1.6799999999999998d-2*gammaaa*t27*t37*t59*t9-8.4d-3*t
     1         27*t28*t9-8.4d-3*gammaaa*t103*t28*t9-2.2399999999999998d-
     2         2*t22*t7+2.2399999999999998d-2*gammaaa*t22*t28*t37
            t105 = 3.333333333333333d-1*t21*t38+t104*t11
            t106 = -3.8035d-1*t1*t11*t29*t38*t53+7.607d-1*t1*t104*t11-3.
     1         8035d-1*t1*t105+t101
            t107 = t100-1.0d+0*t106*t18*wght
            t108 = -2.0d+0*t17*t18*t31*t40*wght-5.0d-1*t1*t17*t25*t32*t3
     1         8*wght
            t109 = 8.0d+0*t108*t19*t33
            t110 = 8.333333333333333d-2*t1*t21*t38
            t111 = -2.5d-1*t1*t11*t29*t38*t53+t110+5.0d-1*t1*t104*t11-2.
     1         5d-1*t1*t105
            t112 = 7.607d-1*t1*t25*t31*t38*wght
            t113 = -5.0d-1*t1*t11*t35*t38
            t114 = 2.0d+0*t111*t17*t18*wght+t113+t112
            t115 = 9.403333333333332d-2*t1*t21*t38
            t116 = 5.0d-1*t1*t11*t31*t38
            t117 = -2.821d-1*t1*t11*t29*t38*t53+t116+t115+1.0d+0*t1*t10*
     1         t11*t111+5.642d-1*t1*t104*t11-2.821d-1*t1*t105
            t118 = 3.8035d-1*t1*t25*t32*t38*wght
            t119 = -t35*t40
            t120 = 1.0d+0*t117*t17*t18*wght+t119+t118
            t121 = 3.8035d-1*t1*t11*t29*t38*t53+t101
            t122 = t100-1.0d+0*t121*t18*wght
            t123 = 2.5d-1*t1*t11*t29*t38*t53+t110
            t124 = 2.0d+0*t123*t17*t18*wght+t113+t112
            t125 = 2.821d-1*t1*t11*t29*t38*t53+1.0d+0*t1*t10*t11*t123+t1
     1         16+t115
            t126 = 1.0d+0*t125*t17*t18*wght+t119+t118
            t127 = t38**2
            t128 = t37**2
            t129 = 1/rhoa**5.333333333333333d+0
            t130 = 1/gammaaa
            t131 = 1/t2**3
            t132 = 6.3d-3*t130*t24*t9-6.3d-3*t129*t61-6.3d-3*t131*t4*t5
            t133 = 1.6799999999999998d-2*gammaaa*t128*t59*t9-1.679999999
     1         9999998d-2*t28*t37*t9-8.4d-3*gammaaa*t132*t28*t9
            t134 = 3.8035d-1*t1*t11*t133-3.8035d-1*t1*t11*t127*t53
            t135 = t40**2
            t136 = -8.0d+0*t135*t14*t17*t18*t41*wght
            t137 = -8.0d+0*t1*t17*t19*t25*t33*t38*t40*wght
            t138 = -6.0d+0*t127*t16*t17*t39*t44*t8*wght
            t139 = 2.5d-1*t1*t11*t133-2.5d-1*t1*t11*t127*t53
            t140 = -2.821d-1*t1*t11*t127*t53+1.25d-1*t127*t20*t39+1.0d+0
     1         *t1*t10*t11*t139+2.821d-1*t1*t11*t133
            t141 = 7.607d-1*t1*t25*t38*t40*wght
            t142 = 1.0d+0*t140*t17*t18*wght+t141
            t143 = 3.8035d-1*t127*t39*t8*wght
            t144 = 2.0d+0*t139*t17*t18*wght+t143
            t145 = 2.821d-1*t1*t11*t127*t53+3.75d-1*t127*t20*t39
            t146 = 1/t15**4
            t147 = 2.4d+1*t14*t146*t17*t18*t32**3*wght
            t148 = 2.88d+2*t17*t18*t19*t32*t44*t45*wght
            t149 = 1/t12**5.0d+0
            t150 = 7.68d+2*t149*t16*t17*t18*t31**3*wght
            t151 = 9.6d+1*t17*t18*t31*t33*t41*t42*wght
            t152 = 1/t18
            t153 = 3.8034999999999997d-1*t1*t4*t53*t54
            t154 = -7.606999999999999d-1*t1*t29*t30*t50*t53
            t155 = 1/t10**2
            t156 = t30**3
            t157 = -5.705250000000001d-1*t1*t155*t156*t21
            t158 = 1.1410500000000001d+0*t1*t155*t29*t54*t55
            t159 = 1/rhoa**5.666666666666667d+0
            t160 = 1/t6**4
            t161 = 1/t23**5
            t162 = 5.04d-2*gammaaa*t27*t59*t64*t9-8.4d-3*gammaaa*t28*(5.
     1         525333333333333d-1*t60*t61/rhoa**8.333333333333334d+0-2.6
     2         133333333333336d-1*t2*t5/rhoa**4.333333333333333d+0-8.885
     3         333333333333d-1*gammaaa*t159*t24-1.7919999999999997d-1*ga
     4         mmaaa**3*t161/rhoa**11)*t9-5.04d-2*gammaaa*t160*t27**3*t9
     5         -3.8328888888888885d-1*gammaaa*t159*t7+6.719999999999998d
     6         -2*gammaaa*t22*t28*t64-1.3439999999999996d-1*gammaaa*t22*
     7         t58*t59-2.4639999999999995d-1*gammaaa*t27*t28*t57
            t163 = 3.7037037037037035d-1*t10*t9+t21*t65-6.66666666666666
     1         6d-1*t29*t49+t11*t162
            t164 = 2.0d+0*rhoa*t32*t34*wght
            t165 = 2.821d-1*t1*t4*t53*t54
            t166 = -5.642d-1*t1*t29*t30*t50*t53
            t167 = -3.333333333333333d-1*t1*t30*t31*t50
            t168 = -4.2315d-1*t1*t155*t156*t21
            t169 = 6.666666666666666d-1*t1*t21*t29*t31
            t170 = 8.463d-1*t1*t155*t29*t54*t55
            t171 = 5.0d-1*t1*t31*t53*t54*t55
            t172 = 2.5d-1*t1*t4*t53*t54
            t173 = -5.0d-1*t1*t29*t30*t50*t53
            t174 = -3.75d-1*t1*t155*t156*t21
            t175 = 7.5d-1*t1*t155*t29*t54*t55
            t176 = 1.8518518518518517d-1*t1*t10*t9+7.5d-1*t1*t30*t53*t55
     1         *t66-7.5d-1*t1*t29*t53*t66-2.5d-1*t1*t50*t66-7.5d-1*t1*t3
     2         0*t53*t65+5.0d-1*t1*t21*t65-3.333333333333333d-1*t1*t29*t
     3         49+1.6666666666666666d-1*t1*t152*t30+t175+t174+t173+t172-
     4         2.5d-1*t1*t163+5.0d-1*t1*t11*t162
            t177 = -1.0d+0*t1*t29*t30*t31*t53
            t178 = 4.0d+0*rhoa*t31*t34*wght
            t179 = -2.0d+0*rhoa*t17*t42*wght
            t180 = -2.0d+0*t18*t34*t42*wght
            t181 = -2*t32*t86
            t182 = -4.0d+0*rhoa*t17*t31*t32*wght
            t183 = -4.0d+0*t18*t31*t32*t34*wght
            t184 = -2.0d+0*t31*t86
            t185 = -t32*t90
            t186 = -6.0d+0*rhoa*t17*t45*wght
            t187 = -6.0d+0*t18*t34*t45*wght
            t188 = -3.0d+0*t31*t90
            t189 = 3.333333333333333d-1*t1*t21*t53*t91-5.0d-1*t1*t155*t3
     1         0*t91-3.7037037037037035d-2*t1*t10*t9+2.5d-1*t1*t30*t53*t
     2         55*t66-2.5d-1*t1*t29*t53*t66-8.333333333333333d-2*t1*t50*
     3         t66-2.5d-1*t1*t30*t53*t65+5.0d-1*t1*t11*t29*t53*t65+1.666
     4         6666666666666d-1*t1*t21*t65+t175+t174+t173+t172
            t190 = 2.4d+1*t14*t146*t17*t18*t40*t42*wght
            t191 = 1.92d+2*t1*t149*t16*t17*t25*t38*t45*wght
            t192 = -8.452222222222222d-2*t1*t38*t49
            t193 = -1.9017500000000004d-1*t1*t155*t38*t54*t55
            t194 = 1.6799999999999998d-2*gammaaa*t37*t59*t64*t9-8.4d-3*t
     1         28*t64*t9-8.4d-3*gammaaa*t28*(3.92d-2*t36*t5*t63-1.735999
     2         9999999998d-1*gammaaa*t61*t62+6.719999999999998d-2*t161*t
     3         60/rhoa**10+2.072d-1*t24*t57)*t9+1.6799999999999998d-2*t5
     4         8*t59*t9+3.3599999999999997d-2*gammaaa*t103*t27*t59*t9-5.
     5         04d-2*gammaaa*t160*t37*t58*t9+8.213333333333332d-2*t57*t7
     6         -8.959999999999999d-2*gammaaa*t22*t27*t37*t59-8.213333333
     7         333332d-2*gammaaa*t28*t37*t57+4.4799999999999995d-2*t22*t
     8         27*t28+4.4799999999999995d-2*gammaaa*t103*t22*t28
            t195 = -2.222222222222222d-1*t38*t49+6.666666666666666d-1*t1
     1         04*t21+t11*t194
            t196 = 8.0d+0*t33*t41*(1.0d+0*t1*t17*t25*t38*t42*wght+8.0d+0
     1         *t17*t18*t31*t32*t40*wght)
            t197 = 1.6d+1*t19*t44*(6.0d+0*t17*t18*t40*t45*wght+3.0d+0*t1
     1         *t17*t25*t31*t32*t38*wght)
            t198 = 2.0d+0*rhoa*t111*t17*wght
            t199 = 7.607d-1*t1*t3*t31*t38*wght
            t200 = 2.0d+0*t111*t18*t34*wght
            t201 = -5.555555555555556d-2*t1*t38*t49
            t202 = -1.25d-1*t1*t155*t38*t54*t55
            t203 = -2.5d-1*t1*t11*t38*t53*t65+5.0d-1*t1*t105*t30*t53*t55
     1         -1.6666666666666666d-1*t1*t21*t29*t38*t53-5.0d-1*t1*t104*
     2         t30*t53-5.0d-1*t1*t105*t29*t53-1.6666666666666666d-1*t1*t
     3         105*t50+5.0d-1*t1*t155*t29*t30*t38+3.333333333333333d-1*t
     4         1*t104*t21+t202+t201-2.5d-1*t1*t195+5.0d-1*t1*t11*t194
            t204 = 7.607d-1*t1*t25*t38*t81*wght
            t205 = -2.0d+0*t107*t31
            t206 = -5.0d-1*t1*t11*t38*t69
            t207 = 1.0d+0*rhoa*t117*t17*wght
            t208 = 3.8035d-1*t1*t3*t32*t38*wght
            t209 = 1.0d+0*t117*t18*t34*wght
            t210 = -6.268888888888888d-2*t1*t38*t49
            t211 = 3.333333333333333d-1*t1*t10*t111*t21
            t212 = 1.6666666666666666d-1*t1*t21*t31*t38
            t213 = -1.4105d-1*t1*t155*t38*t54*t55
            t214 = 1.0d+0*t1*t11*t111*t29
            t215 = 5.0d-1*t1*t11*t38*t81
            t216 = -5.0d-1*t1*t111*t30
            t217 = 3.8035d-1*t1*t25*t38*t84*wght
            t218 = -t107*t32
            t219 = -t40*t69
            t220 = -7.607d-1*t1*t25*t38*t42*wght
            t221 = -2*t40*t86
            t222 = -2.2821000000000002d+0*t1*t25*t38*t45*wght
            t223 = -7.5d-1*t1*t11*t38*t90
            t224 = -1.5214d+0*t1*t25*t31*t32*t38*wght
            t225 = -5.0d-1*t1*t11*t38*t86
            t226 = -t40*t90
            t227 = 2.5d-1*t1*t11*t38*t53*t65+1.6666666666666666d-1*t1*t2
     1         1*t29*t38*t53+t202+t201
            t228 = 2.0d+0*t121*t18*t31*wght
            t229 = -2.0d+0*t123*t35
            t230 = 1.0d+0*t121*t18*t32*wght
            t231 = 5.0d-1*t1*t11*t29*t31*t38*t53
            t232 = -t125*t35
            t233 = -2.5d-1*t1*t11*t155*t38*t91+2.5d-1*t1*t105*t30*t53*t5
     1         5-2.5d-1*t1*t104*t30*t53+5.0d-1*t1*t104*t11*t29*t53-2.5d-
     2         1*t1*t105*t29*t53-8.333333333333333d-2*t1*t105*t50+2.7777
     3         77777777778d-2*t1*t38*t49+2.5d-1*t1*t155*t29*t30*t38+1.66
     4         66666666666666d-1*t1*t104*t21+t202
            t234 = 2.4d+1*t135*t14*t146*t17*t18*t32*wght
            t235 = 4.8d+1*t127*t149*t16*t17*t31*t39*t8*wght
            t236 = -1.0d+0*rhoa*t134*wght
            t237 = -1.2678333333333333d-1*t1*t127*t21*t53
            t238 = 1.2678333333333333d-1*t1*t133*t21
            t239 = -8.4d-3*gammaaa*t28*(4.1999999999999993d-2*t102*t61+8
     1         .4d-3*t131*t26*t5-8.4d-3*t130*t22*t24-2.52d-2*gammaaa*t16
     2         1/rhoa**9)*t9+3.3599999999999997d-2*t27*t37*t59*t9+3.3599
     3         999999999997d-2*gammaaa*t103*t37*t59*t9+1.679999999999999
     4         8d-2*gammaaa*t132*t27*t59*t9-1.6799999999999998d-2*t103*t
     5         28*t9-5.04d-2*gammaaa*t128*t160*t27*t9-4.4799999999999995
     6         d-2*gammaaa*t128*t22*t59+4.4799999999999995d-2*t22*t28*t3
     7         7+2.2399999999999998d-2*gammaaa*t132*t22*t28
            t240 = t11*t239+3.333333333333333d-1*t133*t21
            t241 = 1.9017499999999998d-1*t1*t127*t155*t30
            t242 = 2.0d+0*t134*t18*t31*wght
            t243 = -8.333333333333333d-2*t1*t127*t21*t53
            t244 = 8.333333333333333d-2*t1*t133*t21
            t245 = 1.25d-1*t1*t127*t155*t30
            t246 = -5.0d-1*t1*t104*t11*t38*t53-2.5d-1*t1*t11*t133*t29*t5
     1         3+2.5d-1*t1*t11*t127*t155*t29+t245+t244+t243-2.5d-1*t1*t2
     2         40+5.0d-1*t1*t11*t239
            t247 = -2.0d+0*t139*t35
            t248 = 1.0d+0*t134*t18*t32*wght
            t249 = -9.403333333333332d-2*t1*t127*t21*t53
            t250 = 9.403333333333332d-2*t1*t133*t21
            t251 = -5.0d-1*t1*t11*t127*t31*t53
            t252 = 5.0d-1*t1*t11*t133*t31
            t253 = 1.4105d-1*t1*t127*t155*t30
            t254 = -t140*t35
            t255 = 1.6d+1*t19*t44*(1.5d+0*t1*t17*t25*t31*t38*t40*wght-7.
     1         5d-1*t1*t108*t11*t38)
            t256 = 8.0d+0*t33*t41*(1.0d+0*t1*t17*t25*t32*t38*t40*wght-2*
     1         t108*t40)
            t257 = -6.0d+0*t139*t17*t18*t31*wght
            t258 = -5.705250000000001d-1*t127*t31*t39*t8*wght
            t259 = -7.5d-1*t1*t11*t114*t38
            t260 = -2.0d+0*t139*t17*t18*t32*wght
            t261 = -2.0d+0*t140*t17*t18*t31*wght
            t262 = -7.607d-1*t1*t25*t31*t38*t40*wght
            t263 = -1.90175d-1*t127*t32*t39*t8*wght
            t264 = -t114*t40
            t265 = -5.0d-1*t1*t11*t120*t38
            t266 = -2.0d+0*t140*t17*t18*t32*wght
            t267 = -7.607d-1*t1*t25*t32*t38*t40*wght
            t268 = -2*t120*t40
            t269 = 2.5d-1*t1*t104*t11*t38*t53+8.333333333333333d-2*t1*t1
     1         27*t21*t53-1.25d-1*t1*t127*t155*t30
            t270 = -1.5d+0*t1*t123*t17*t25*t38*wght
            t271 = -2.0d+0*t123*t17*t18*t40*wght
            t272 = -5.0d-1*t1*t125*t17*t25*t38*wght
            t273 = -2.0d+0*t125*t17*t18*t40*wght
            t274 = 2.5d-1*t1*t11*t133*t29*t53-2.5d-1*t1*t11*t127*t155*t2
     1         9+t245+t244+t243
            t275 = t38**3
            t276 = -8.4d-3*gammaaa*t28*t9*(-9.450000000000001d-3*t24*t9/
     1         t60-3.15d-3*t129*t130*t61+9.450000000000001d-3*t4*t5/t2**
     2         5+9.450000000000001d-3*t161/rhoa**8)+5.04d-2*gammaaa*t132
     3         *t37*t59*t9+5.04d-2*t128*t59*t9-5.04d-2*gammaaa*t160*t37*
     4         *3*t9-2.52d-2*t132*t28*t9
            t277 = 2.4d+1*t14*t146*t17*t18*t40**3*wght
            t278 = 2.4d+1*t1*t135*t17*t25*t33*t38*t41*wght
            t279 = 1.8d+1*t127*t17*t19*t39*t40*t44*t8*wght
            t280 = t1**3
            t281 = rhoa**3
            t282 = 1.2d+1*t149*t16*t17*t275*t280*t281*wght
            t283 = -7.5d-1*t1*t11*t133*t38*t53+2.5d-1*t1*t11*t276+3.75d-
     1         1*t1*t11*t155*t275
            t284 = -7.607d-1*t1*t135*t25*t38*wght
            t285 = -2*t142*t40
            t286 = -1.4263125000000001d-1*t275*t280*t281*wght
            t287 = -7.5d-1*t1*t11*t144*t38
            t288 = -3.8035d-1*t127*t39*t40*t8*wght
            t289 = -5.0d-1*t1*t11*t142*t38
            t290 = -t144*t40
            t291 = 2.5d-1*t1*t11*t133*t38*t53-1.25d-1*t1*t11*t155*t275
            fnc(iq) = fnc(iq)-4.0d+0*t14*t16*t17*t18*wght
            Amat(iq,D1_RA) = 1.6d+1*t16*t17*t18*t31*t33*wght+4.0d+0*t14*
     1         t17*t18*t19*t32*wght+4.0d+0*t14*t16*t35+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 4.0d+0*t14*t17*t18*t19*t40*wght+4.0d+0*t1*
     1         t16*t17*t25*t33*t38*wght-1.5214d+0*t1*t14*t16*t25*t38*wgh
     2         t+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = 8.0d+0*t16*t33*t90+4.0d+0*t14*t19*t86+4
     1         .0d+0*t14*t16*t69+t47+t46+t43+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = 8.0d+0*t16*t33*t97+4.0d+0*t14*t19*t96+4
     1         .0d+0*t14*t16*t93+t47+t46+t43+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t99+t98+8.0d+0*t114*t16*t33+4.0d+0*t12
     1         0*t14*t19+4.0d+0*t107*t14*t16+t109+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t99+t98+8.0d+0*t124*t16*t33+4.0d+0*t12
     1         6*t14*t19+4.0d+0*t122*t14*t16+t109+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = -4.0d+0*t134*t14*t16*t18*wght+8.0d+0*
     1         t144*t16*t33+4.0d+0*t14*t142*t19+t138+t137+t136+Cmat2(iq,
     2         D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = 8.0d+0*t16*t33*(5.0d-1*t1*t127*t17*t2
     1         5*t53*wght+t143)+4.0d+0*t14*t19*(1.0d+0*t145*t17*t18*wght
     2         +t141)-1.5214d+0*t1*t127*t14*t16*t25*t53*wght+t138+t137+t
     3         136+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = 4.0d+0*t14*t16*(-1.0d+0*t18*(2.81740
     1         74074074074d-1*t1*t10*t9+1.1410500000000001d+0*t1*t30*t53
     2         *t55*t66-1.1410500000000001d+0*t1*t29*t53*t66-3.803499999
     3         9999997d-1*t1*t50*t66-1.1410500000000001d+0*t1*t30*t53*t6
     4         5+7.606999999999999d-1*t1*t21*t65-5.071333333333333d-1*t1
     5         *t29*t49+2.5356666666666666d-1*t1*t152*t30-3.8035d-1*t1*t
     6         163+7.607d-1*t1*t11*t162+t158+t157+t154+t153)*wght-3.0d+0
     7         *rhoa*t68*wght)+4.0d+0*t14*t19*(1.0d+0*t17*t18*(2.0896296
     8         296296293d-1*t1*t10*t9-1.0d+0*t1*t30*t81+2.0d+0*t1*t11*t2
     9         9*t81+6.666666666666666d-1*t1*t10*t21*t81+8.463d-1*t1*t30
     :         *t53*t55*t66-8.463d-1*t1*t29*t53*t66-2.821d-1*t1*t50*t66-
     ;         5.0d-1*t1*t31*t66-8.463d-1*t1*t30*t53*t65+1.0d+0*t1*t11*t
     <         31*t65+5.642d-1*t1*t21*t65-2.222222222222222d-1*t1*t10*t3
     =         1*t49-3.7613333333333326d-1*t1*t29*t49+1.8806666666666663
     >         d-1*t1*t152*t30+t177+1.0d+0*t1*t10*t11*t176+t171+t170+t16
     ?         9+t168+t167+t166+t165-2.821d-1*t1*t163+5.642d-1*t1*t11*t1
     @         62)*wght+2.0d+0*t18*t34*t84*wght+2.0d+0*rhoa*t17*t84*wght
     1         +1.0d+0*t18*t32*t68*wght-t35*t84-2*t32*t69+t164)+4.0d+0*t
     2         14*t41*(-4.0d+0*t17*t18*t32*t84*wght+t181+t180+t179)+8.0d
     3         +0*t19*t33*(-4.0d+0*t17*t18*t31*t84*wght-4.0d+0*t17*t18*t
     4         32*t81*wght+t185+t184+t183+t182)+8.0d+0*t16*t33*(4.0d+0*t
     5         18*t34*t81*wght+4.0d+0*rhoa*t17*t81*wght+2.0d+0*t18*t31*t
     6         68*wght+2.0d+0*t17*t176*t18*wght-2.0d+0*t35*t81-4.0d+0*t3
     7         1*t69+t178)+1.6d+1*t16*t44*(-1.2d+1*t17*t18*t31*t81*wght+
     8         t188+t187+t186)+t151+t150+t148+t147+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = 4.0d+0*t14*t19*(1.0d+0*t18*t34*t95*w
     1         ght+1.0d+0*rhoa*t17*t95*wght+1.0d+0*t17*t18*(-5.0d-1*t1*t
     2         30*t94+1.0d+0*t1*t11*t29*t94+3.333333333333333d-1*t1*t10*
     3         t21*t94+1.0d+0*t1*t11*t31*t53*t91+3.761333333333333d-1*t1
     4         *t21*t53*t91-5.642d-1*t1*t155*t30*t91-4.179259259259258d-
     5         2*t1*t10*t9-5.0d-1*t1*t30*t81+1.0d+0*t1*t11*t29*t81+3.333
     6         333333333333d-1*t1*t10*t21*t81+2.821d-1*t1*t30*t53*t55*t6
     7         6-2.821d-1*t1*t29*t53*t66-9.403333333333332d-2*t1*t50*t66
     8         -2.821d-1*t1*t30*t53*t65+5.642d-1*t1*t11*t29*t53*t65+1.88
     9         06666666666666d-1*t1*t21*t65+1.111111111111111d-1*t1*t10*
     :         t31*t49+1.0d+0*t1*t10*t11*t189+t177+t171+t170+t169+t168+t
     ;         167+t166+t165)*wght+1.0d+0*t18*t32*t92*wght+1.0d+0*t18*t3
     <         4*t84*wght+1.0d+0*rhoa*t17*t84*wght+1.0d+0*t17*t32*wght-t
     =         35*t95-t32*t93-t32*t69+t164)+4.0d+0*t14*t41*(-4.0d+0*t17*
     >         t18*t32*t95*wght+t181+t180+t179)+8.0d+0*t19*t33*(-4.0d+0*
     ?         t17*t18*t31*t95*wght-4.0d+0*t17*t18*t32*t94*wght+t185+t18
     @         4+t183+t182)+8.0d+0*t16*t33*(2.0d+0*t18*t34*t94*wght+2.0d
     1         +0*rhoa*t17*t94*wght+2.0d+0*t18*t31*t92*wght+2.0d+0*t18*t
     2         34*t81*wght+2.0d+0*rhoa*t17*t81*wght+2.0d+0*t17*t31*wght+
     3         2.0d+0*t17*t18*t189*wght-2.0d+0*t35*t94-2.0d+0*t31*t93-2.
     4         0d+0*t31*t69+t178)+1.6d+1*t16*t44*(-1.2d+1*t17*t18*t31*t9
     5         4*wght+t188+t187+t186)+4.0d+0*t14*t16*(-2.0d+0*rhoa*t92*w
     6         ght-1.0d+0*t18*(5.071333333333333d-1*t1*t21*t53*t91-7.607
     7         d-1*t1*t155*t30*t91-5.634814814814815d-2*t1*t10*t9+3.8035
     8         d-1*t1*t30*t53*t55*t66-3.8035d-1*t1*t29*t53*t66-1.2678333
     9         333333333d-1*t1*t50*t66-3.8035d-1*t1*t30*t53*t65+7.607d-1
     :         *t1*t11*t29*t53*t65+2.5356666666666666d-1*t1*t21*t65+t158
     ;         +t157+t154+t153)*wght-1.0d+0*rhoa*t68*wght-2.0d+0*t34*wgh
     <         t)+t151+t150+t148+t147+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(-2.
     1         821d-1*t1*t11*t38*t53*t65+5.642d-1*t1*t105*t30*t53*t55-5.
     2         0d-1*t1*t11*t29*t31*t38*t53-1.8806666666666666d-1*t1*t21*
     3         t29*t38*t53-5.642d-1*t1*t104*t30*t53-5.642d-1*t1*t105*t29
     4         *t53-1.8806666666666666d-1*t1*t105*t50+5.642d-1*t1*t155*t
     5         29*t30*t38+1.0d+0*t1*t104*t11*t31-5.0d-1*t1*t105*t31+t216
     6         +t215+t214+t213+t212+t211+t210+3.761333333333333d-1*t1*t1
     7         04*t21+1.0d+0*t1*t10*t11*t203-2.821d-1*t1*t195+5.642d-1*t
     8         1*t11*t194)*wght+1.0d+0*t106*t18*t32*wght-t117*t35+t219+t
     9         218+t217+t209+t208+t207)+4.0d+0*t14*t16*(-1.0d+0*t18*(-3.
     :         8035d-1*t1*t11*t38*t53*t65+7.607d-1*t1*t105*t30*t53*t55-2
     ;         .5356666666666666d-1*t1*t21*t29*t38*t53-7.607d-1*t1*t104*
     <         t30*t53-7.607d-1*t1*t105*t29*t53-2.5356666666666666d-1*t1
     =         *t105*t50+7.607d-1*t1*t155*t29*t30*t38+5.071333333333333d
     >         -1*t1*t104*t21-3.8035d-1*t1*t195+7.607d-1*t1*t11*t194+t19
     ?         3+t192)*wght-2.0d+0*rhoa*t106*wght)+4.0d+0*t14*t41*(-4.0d
     @         +0*t117*t17*t18*t32*wght+t221+t220)+8.0d+0*t19*t33*(-4.0d
     1         +0*t111*t17*t18*t32*wght-4.0d+0*t117*t17*t18*t31*wght+t22
     2         6+t225+t224)+1.6d+1*t16*t44*(-1.2d+1*t111*t17*t18*t31*wgh
     3         t+t223+t222)+8.0d+0*t16*t33*(2.0d+0*t106*t18*t31*wght+2.0
     4         d+0*t17*t18*t203*wght-2.0d+0*t111*t35+t206+t205+t204+t200
     5         +t199+t198)+t197+t196+t191+t190+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = 4.0d+0*t14*t16*(-1.0d+0*t18*(3.8035
     1         d-1*t1*t11*t38*t53*t65+2.5356666666666666d-1*t1*t21*t29*t
     2         38*t53+t193+t192)*wght-2.0d+0*rhoa*t121*wght)+4.0d+0*t14*
     3         t19*(1.0d+0*t17*t18*(2.821d-1*t1*t11*t38*t53*t65+1.880666
     4         6666666666d-1*t1*t21*t29*t38*t53-5.0d-1*t1*t123*t30+1.0d+
     5         0*t1*t11*t123*t29+t231+1.0d+0*t1*t10*t11*t227+t215+t213+t
     6         212+t210+3.333333333333333d-1*t1*t10*t123*t21)*wght+1.0d+
     7         0*t125*t18*t34*wght+1.0d+0*rhoa*t125*t17*wght-t122*t32+t2
     8         32+t230+t219+t217+t208)+8.0d+0*t16*t33*(2.0d+0*t123*t18*t
     9         34*wght+2.0d+0*t17*t18*t227*wght+2.0d+0*rhoa*t123*t17*wgh
     :         t-2.0d+0*t122*t31+t229+t228+t206+t204+t199)+4.0d+0*t14*t4
     ;         1*(-4.0d+0*t125*t17*t18*t32*wght+t221+t220)+8.0d+0*t19*t3
     <         3*(-4.0d+0*t123*t17*t18*t32*wght-4.0d+0*t125*t17*t18*t31*
     =         wght+t226+t225+t224)+1.6d+1*t16*t44*(-1.2d+1*t123*t17*t18
     >         *t31*wght+t223+t222)+t197+t196+t191+t190+Cmat3(iq,D3_RA_R
     ?         A_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = 4.0d+0*t14*t19*(3.8035d-1*t1*t25*t3
     1         8*t95*wght+1.0d+0*t17*t18*(5.0d-1*t1*t11*t38*t94-2.821d-1
     2         *t1*t11*t155*t38*t91+2.821d-1*t1*t105*t30*t53*t55-2.821d-
     3         1*t1*t104*t30*t53+5.642d-1*t1*t104*t11*t29*t53-2.821d-1*t
     4         1*t105*t29*t53-9.403333333333332d-2*t1*t105*t50+3.1344444
     5         44444444d-2*t1*t38*t49+2.821d-1*t1*t155*t29*t30*t38+1.0d+
     6         0*t1*t10*t11*t233+t231+t216+t214+t213+t212+t211+1.8806666
     7         666666666d-1*t1*t104*t21)*wght-t40*t93+t232+t230+t218+t20
     8         9+t208+t207)+8.0d+0*t16*t33*(7.607d-1*t1*t25*t38*t94*wght
     9         +2.0d+0*t17*t18*t233*wght-5.0d-1*t1*t11*t38*t93+t229+t228
     :         +t205+t200+t199+t198)+4.0d+0*t14*t16*(-1.0d+0*t18*(-3.803
     ;         5d-1*t1*t11*t155*t38*t91+3.8035d-1*t1*t105*t30*t53*t55-3.
     <         8035d-1*t1*t104*t30*t53+7.607d-1*t1*t104*t11*t29*t53-3.80
     =         35d-1*t1*t105*t29*t53-1.2678333333333333d-1*t1*t105*t50+4
     >         .226111111111111d-2*t1*t38*t49+3.8035d-1*t1*t155*t29*t30*
     ?         t38+2.5356666666666666d-1*t1*t104*t21+t193)*wght-3.8035d-
     @         1*t1*t11*t38*wght-1.0d+0*rhoa*t121*wght-1.0d+0*rhoa*t106*
     1         wght)+4.0d+0*t14*t41*(-2.0d+0*t125*t17*t18*t32*wght-2.0d+
     2         0*t117*t17*t18*t32*wght-2*t40*t96+t220)+8.0d+0*t19*t33*(-
     3         2.0d+0*t123*t17*t18*t32*wght-2.0d+0*t111*t17*t18*t32*wght
     4         -2.0d+0*t125*t17*t18*t31*wght-2.0d+0*t117*t17*t18*t31*wgh
     5         t-t40*t97-5.0d-1*t1*t11*t38*t96+t224)+1.6d+1*t16*t44*(-6.
     6         0d+0*t123*t17*t18*t31*wght-6.0d+0*t111*t17*t18*t31*wght-7
     7         .5d-1*t1*t11*t38*t97+t222)+t197+t196+t191+t190+Cmat3(iq,D
     8         3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(-5
     1         .642d-1*t1*t104*t11*t38*t53-2.821d-1*t1*t11*t133*t29*t53+
     2         1.0d+0*t1*t11*t111*t38+2.821d-1*t1*t11*t127*t155*t29+t253
     3         +t252+t251+t250+t249+1.0d+0*t1*t10*t11*t246-2.821d-1*t1*t
     4         240+5.642d-1*t1*t11*t239)*wght+7.607d-1*t1*t117*t25*t38*w
     5         ght-2*t107*t40+t254+t248)+4.0d+0*t14*t16*(t236-1.0d+0*t18
     6         *(-7.607d-1*t1*t104*t11*t38*t53-3.8035d-1*t1*t11*t133*t29
     7         *t53+3.8035d-1*t1*t11*t127*t155*t29+t241-3.8035d-1*t1*t24
     8         0+7.607d-1*t1*t11*t239+t238+t237)*wght)+4.0d+0*t14*t41*(-
     9         2.0d+0*t117*t17*t18*t40*wght+t268+t267+t266)+8.0d+0*t19*t
     :         33*(-2.0d+0*t111*t17*t18*t40*wght-5.0d-1*t1*t117*t17*t25*
     ;         t38*wght+t265+t264+t263+t262+t261+t260)+1.6d+1*t16*t44*(-
     <         1.5d+0*t1*t111*t17*t25*t38*wght+t259+t258+t257)+8.0d+0*t1
     =         6*t33*(1.5214d+0*t1*t111*t25*t38*wght+2.0d+0*t17*t18*t246
     >         *wght-1.0d+0*t1*t107*t11*t38+t247+t242)+t256+t255+t235+t2
     ?         34+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = 4.0d+0*t14*t16*(-1.0d+0*t18*(3.803
     1         5d-1*t1*t104*t11*t38*t53+1.2678333333333333d-1*t1*t127*t2
     2         1*t53-1.9017500000000004d-1*t1*t127*t155*t30)*wght-3.8035
     3         d-1*t1*t127*t3*t53*wght)+4.0d+0*t14*t19*(1.0d+0*t17*t18*(
     4         2.821d-1*t1*t104*t11*t38*t53+5.0d-1*t1*t11*t127*t31*t53+9
     5         .403333333333332d-2*t1*t127*t21*t53+5.0d-1*t1*t11*t123*t3
     6         8+5.0d-1*t1*t11*t111*t38-1.4105d-1*t1*t127*t155*t30+1.0d+
     7         0*t1*t10*t11*t269)*wght+3.8035d-1*t1*t127*t25*t32*t53*wgh
     8         t+3.8035d-1*t1*t125*t25*t38*wght+3.8035d-1*t1*t117*t25*t3
     9         8*wght-t122*t40-t107*t40-t145*t35)+8.0d+0*t19*t33*(-5.0d-
     :         1*t1*t127*t17*t25*t32*t53*wght-2.0d+0*t145*t17*t18*t31*wg
     ;         ht+t272+t271+t265+t264+t263+t262)+1.6d+1*t16*t44*(-1.5d+0
     <         *t1*t127*t17*t25*t31*t53*wght+t270+t259+t258)+8.0d+0*t16*
     =         t33*(7.607d-1*t1*t127*t25*t31*t53*wght+7.607d-1*t1*t123*t
     >         25*t38*wght+7.607d-1*t1*t111*t25*t38*wght+2.0d+0*t17*t18*
     ?         t269*wght-5.0d-1*t1*t11*t127*t35*t53-5.0d-1*t1*t11*t122*t
     @         38-5.0d-1*t1*t107*t11*t38)+4.0d+0*t14*t41*(-2.0d+0*t145*t
     1         17*t18*t32*wght+t273+t268+t267)+t256+t255+t235+t234+Cmat3
     2         (iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = 4.0d+0*t14*t16*(t236-1.0d+0*t18*(3
     1         .8035d-1*t1*t11*t133*t29*t53-3.8035d-1*t1*t11*t127*t155*t
     2         29+t241+t238+t237)*wght)+4.0d+0*t14*t19*(1.0d+0*t17*t18*(
     3         2.821d-1*t1*t11*t133*t29*t53+1.0d+0*t1*t11*t123*t38-2.821
     4         d-1*t1*t11*t127*t155*t29+1.0d+0*t1*t10*t11*t274+t253+t252
     5         +t251+t250+t249)*wght+7.607d-1*t1*t125*t25*t38*wght-2*t12
     6         2*t40+t254+t248)+8.0d+0*t16*t33*(1.5214d+0*t1*t123*t25*t3
     7         8*wght+2.0d+0*t17*t18*t274*wght-1.0d+0*t1*t11*t122*t38+t2
     8         47+t242)+1.6d+1*t16*(-7.5d-1*t1*t11*t124*t38+t270+t258+t2
     9         57)*t44+4.0d+0*t14*(-2*t126*t40+t273+t267+t266)*t41+8.0d+
     :         0*t19*t33*(-t124*t40-5.0d-1*t1*t11*t126*t38+t272+t271+t26
     ;         3+t262+t261+t260)+t256+t255+t235+t234+Cmat3(iq,D3_RA_GBB_
     <         GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = 4.0d+0*t14*t19*(1.0d+0*t17*t18*(-
     1         1.25d-1*t20*t275*t39*t53-8.463d-1*t1*t11*t133*t38*t53+1.2
     2         5d-1*t133*t20*t38*t39+1.0d+0*t1*t11*t139*t38+1.0d+0*t1*t1
     3         0*t11*t283+2.821d-1*t1*t11*t276+4.2315d-1*t1*t11*t155*t27
     4         5)*wght+3.0d+0*t134*t18*t40*wght+1.1410500000000001d+0*t1
     5         *t140*t25*t38*wght)+4.0d+0*t14*t41*(-4.0d+0*t140*t17*t18*
     6         t40*wght+t285+t284)+8.0d+0*t19*t33*(-4.0d+0*t139*t17*t18*
     7         t40*wght-1.0d+0*t1*t140*t17*t25*t38*wght+t290+t289+t288)+
     8         1.6d+1*t16*t44*(-3.0d+0*t1*t139*t17*t25*t38*wght+t287+t28
     9         6)+8.0d+0*t16*t33*(2.2821000000000002d+0*t1*t139*t25*t38*
     :         wght+1.5d+0*t1*t134*t25*t38*wght+2.0d+0*t17*t18*t283*wght
     ;         )-4.0d+0*t14*t16*t18*(-1.1410500000000001d+0*t1*t11*t133*
     <         t38*t53+3.8035d-1*t1*t11*t276+5.705250000000001d-1*t1*t11
     =         *t155*t275)*wght+t282+t279+t278+t277+Cmat3(iq,D3_GAA_GAA_
     >         GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = 1.6d+1*t16*t44*(-7.5d-1*t17*t275*
     1         t39*t53*t8*wght+t287+t286)+8.0d+0*t16*t33*(7.607d-1*t275*
     2         t39*t53*t8*wght+7.607d-1*t1*t139*t25*t38*wght+5.0d-1*t1*t
     3         134*t25*t38*wght+2.0d+0*t17*t18*t291*wght)+4.0d+0*t14*t19
     4         *(1.0d+0*t17*t18*(2.5d-1*t20*t275*t39*t53+2.821d-1*t1*t11
     5         *t133*t38*t53+5.0d-1*t1*t11*t139*t38+1.0d+0*t1*t10*t11*t2
     6         91-1.4105d-1*t1*t11*t155*t275)*wght+7.607d-1*t1*t127*t25*
     7         t40*t53*wght+1.0d+0*t134*t18*t40*wght+7.607d-1*t1*t145*t2
     8         5*t38*wght+3.8035d-1*t1*t140*t25*t38*wght)+8.0d+0*t19*t33
     9         *(-1.0d+0*t1*t127*t17*t25*t40*t53*wght-1.0d+0*t1*t145*t17
     :         *t25*t38*wght+t290+t289+t288)+4.0d+0*t14*t41*(-4.0d+0*t14
     ;         5*t17*t18*t40*wght+t285+t284)-4.0d+0*t14*t16*t18*(3.8035d
     <         -1*t1*t11*t133*t38*t53-1.901750000000001d-1*t1*t11*t155*t
     =         275)*wght+t282+t279+t278+t277+Cmat3(iq,D3_GAA_GAA_GBB)
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
            t1 = param(1)
            t2 = gammaaa**5.0d-1
            t3 = rhoa**1.3333333333333333d+0
            t4 = 1/t3
            t5 = asinh(t2*t4)
            t6 = 2.52d-2*t2*t4*t5+1.0d+0
            t7 = 1/t6
            t8 = 1/rhoa**2.6666666666666666d+0
            t9 = 8.4d-3*gammaaa*t7*t8+1.8610514726982d+0
            t10 = rhoa**3.333333333333333d-1
            t11 = gammabb**5.0d-1
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t12
            t14 = asinh(t11*t13)
            t15 = 2.52d-2*t11*t13*t14+1.0d+0
            t16 = 1/t15
            t17 = 1/rhob**2.6666666666666666d+0
            t18 = 8.4d-3*gammabb*t16*t17+1.8610514726982d+0
            t19 = rhob**3.333333333333333d-1
            t20 = t10*t9+t18*t19
            t21 = 1/t20
            t22 = t1*t10*t18*t19*t21*t9
            t23 = t22**2.0d+0
            t24 = 1/t23
            t25 = 1.1284000000000001d+0*t1*t10*t18*t19*t21*t9+t23+3.183d
     1         -1
            t26 = 1/t25
            t27 = 1.5214d+0*t1*t10*t18*t19*t21*t9+5.764d-1
            t28 = 1/t22**3.0d+0
            t29 = rhoa**6.666666666666666d-1
            t30 = 1/t29
            t31 = 1/rhoa**3.6666666666666664d+0
            t32 = (gammaaa*t8+1)**5.0d-1
            t33 = 1/t32
            t34 = 1/rhoa**2.3333333333333334d+0
            t35 = -3.3599999999999997d-2*t2*t34*t5-3.3599999999999997d-2
     1         *gammaaa*t31*t33
            t36 = 1/t6**2
            t37 = -8.4d-3*gammaaa*t35*t36*t8-2.2399999999999998d-2*gamma
     1         aa*t31*t7
            t38 = 3.333333333333333d-1*t30*t9+t10*t37
            t39 = 1/t20**2
            t40 = -t1*t10*t18*t19*t38*t39*t9+3.333333333333333d-1*t1*t18
     1         *t19*t21*t30*t9+t1*t10*t18*t19*t21*t37
            t41 = 1/t25**2
            t42 = 2.0d+0*t1*t10*t18*t19*t21*t40*t9-1.1284000000000001d+0
     1         *t1*t10*t18*t19*t38*t39*t9+3.761333333333333d-1*t1*t18*t1
     2         9*t21*t30*t9+1.1284000000000001d+0*t1*t10*t18*t19*t21*t37
            t43 = -1.5214d+0*t1*t10*t18*t19*t38*t39*t9+5.071333333333333
     1         d-1*t1*t18*t19*t21*t30*t9+1.5214d+0*t1*t10*t18*t19*t21*t3
     2         7
            t44 = -1.0d+0*rhoa*rhob*t43*wght-1.0d+0*rhob*t27*wght
            t45 = rhob**6.666666666666666d-1
            t46 = 1/t45
            t47 = 1/rhob**3.6666666666666664d+0
            t48 = (gammabb*t17+1)**5.0d-1
            t49 = 1/t48
            t50 = 1/rhob**2.3333333333333334d+0
            t51 = -3.3599999999999997d-2*t11*t14*t50-3.3599999999999997d
     1         -2*gammabb*t47*t49
            t52 = 1/t15**2
            t53 = -8.4d-3*gammabb*t17*t51*t52-2.2399999999999998d-2*gamm
     1         abb*t16*t47
            t54 = t19*t53+3.333333333333333d-1*t18*t46
            t55 = -t1*t10*t18*t19*t39*t54*t9+t1*t10*t19*t21*t53*t9+3.333
     1         333333333333d-1*t1*t10*t18*t21*t46*t9
            t56 = 2.0d+0*t1*t10*t18*t19*t21*t55*t9-1.1284000000000001d+0
     1         *t1*t10*t18*t19*t39*t54*t9+1.1284000000000001d+0*t1*t10*t
     2         19*t21*t53*t9+3.761333333333333d-1*t1*t10*t18*t21*t46*t9
            t57 = -1.5214d+0*t1*t10*t18*t19*t39*t54*t9+1.5214d+0*t1*t10*
     1         t19*t21*t53*t9+5.071333333333333d-1*t1*t10*t18*t21*t46*t9
            t58 = -1.0d+0*rhoa*rhob*t57*wght-1.0d+0*rhoa*t27*wght
            t59 = 1/t2
            t60 = 1.26d-2*t33*t8+1.26d-2*t4*t5*t59
            t61 = 8.4d-3*t7*t8-8.4d-3*gammaaa*t36*t60*t8
            t62 = 1.5214d+0*t1*t10*t18*t19*t21*t61-1.5214d+0*t1*t18*t19*
     1         t29*t39*t61*t9
            t63 = t1*t10*t18*t19*t21*t61-t1*t18*t19*t29*t39*t61*t9
            t64 = 2.0d+0*t1*t10*t18*t19*t21*t63*t9-1.1284000000000001d+0
     1         *t1*t18*t19*t29*t39*t61*t9+1.1284000000000001d+0*t1*t10*t
     2         18*t19*t21*t61
            t65 = 1/t11
            t66 = 1.26d-2*t13*t14*t65+1.26d-2*t17*t49
            t67 = 8.4d-3*t16*t17-8.4d-3*gammabb*t17*t52*t66
            t68 = 1.5214d+0*t1*t10*t19*t21*t67*t9-1.5214d+0*t1*t10*t18*t
     1         39*t45*t67*t9
            t69 = t1*t10*t19*t21*t67*t9-t1*t10*t18*t39*t45*t67*t9
            t70 = 2.0d+0*t1*t10*t18*t19*t21*t69*t9-1.1284000000000001d+0
     1         *t1*t10*t18*t39*t45*t67*t9+1.1284000000000001d+0*t1*t10*t
     2         19*t21*t67*t9
            t71 = 1/t25**3
            t72 = t42**2
            t73 = 1/t22**4.0d+0
            t74 = t40**2
            t75 = t38**2
            t76 = 1/t20**3
            t77 = 1/rhoa**1.6666666666666669d+0
            t78 = 1/rhoa**4.666666666666667d+0
            t79 = t35**2
            t80 = 1/t6**3
            t81 = gammaaa**2
            t82 = 1/t32**3
            t83 = 1/rhoa**7.333333333333333d+0
            t84 = 1/rhoa**3.3333333333333337d+0
            t85 = 7.84d-2*t2*t5*t84-4.4799999999999995d-2*t81*t82*t83+1.
     1         6799999999999998d-1*gammaaa*t33*t78
            t86 = -8.4d-3*gammaaa*t36*t8*t85+1.6799999999999998d-2*gamma
     1         aa*t79*t8*t80+8.213333333333332d-2*gammaaa*t7*t78+4.47999
     2         99999999995d-2*gammaaa*t31*t35*t36
            t87 = -2.222222222222222d-1*t77*t9+t10*t86+6.666666666666666
     1         d-1*t30*t37
            t88 = -1.5214d+0*t1*t10*t18*t19*t39*t87*t9-3.380888888888889
     1         d-1*t1*t18*t19*t21*t77*t9+3.0428d+0*t1*t10*t18*t19*t75*t7
     2         6*t9-1.0142666666666666d+0*t1*t18*t19*t30*t38*t39*t9+1.52
     3         14d+0*t1*t10*t18*t19*t21*t86-3.0428d+0*t1*t10*t18*t19*t37
     4         *t38*t39+1.0142666666666666d+0*t1*t18*t19*t21*t30*t37
            t89 = -1.0d+0*rhoa*rhob*t88*wght-2.0d+0*rhob*t43*wght
            t90 = -t1*t10*t18*t19*t39*t87*t9-2.222222222222222d-1*t1*t18
     1         *t19*t21*t77*t9+2*t1*t10*t18*t19*t75*t76*t9-6.66666666666
     2         6666d-1*t1*t18*t19*t30*t38*t39*t9+t1*t10*t18*t19*t21*t86-
     3         2*t1*t10*t18*t19*t37*t38*t39+6.666666666666666d-1*t1*t18*
     4         t19*t21*t30*t37
            t91 = 2.0d+0*rhoa*rhob*t27*t90*wght+2.0d+0*rhoa*rhob*t40*t43
     1         *wght+2.0d+0*rhob*t27*t40*wght-2.0d+0*t40*t44
            t92 = 2.0d+0*t1*t10*t18*t19*t21*t9*t90-1.1284000000000001d+0
     1         *t1*t10*t18*t19*t39*t87*t9-2.507555555555555d-1*t1*t18*t1
     2         9*t21*t77*t9+2.2568000000000002d+0*t1*t10*t18*t19*t75*t76
     3         *t9-2.0d+0*t1*t10*t18*t19*t38*t39*t40*t9+6.66666666666666
     4         6d-1*t1*t18*t19*t21*t30*t40*t9-7.522666666666666d-1*t1*t1
     5         8*t19*t30*t38*t39*t9+1.1284000000000001d+0*t1*t10*t18*t19
     6         *t21*t86+2.0d+0*t1*t10*t18*t19*t21*t37*t40-2.256800000000
     7         0002d+0*t1*t10*t18*t19*t37*t38*t39+7.522666666666666d-1*t
     8         1*t18*t19*t21*t30*t37
            t93 = 1.0d+0*rhoa*rhob*t27*t92*wght+1.0d+0*rhoa*rhob*t42*t43
     1         *wght+1.0d+0*rhob*t27*t42*wght-t42*t44
            t94 = 3.0428d+0*t1*t10*t18*t19*t38*t54*t76*t9-5.071333333333
     1         333d-1*t1*t18*t19*t30*t39*t54*t9-1.5214d+0*t1*t10*t19*t38
     2         *t39*t53*t9+5.071333333333333d-1*t1*t19*t21*t30*t53*t9-5.
     3         071333333333333d-1*t1*t10*t18*t38*t39*t46*t9+1.6904444444
     4         444444d-1*t1*t18*t21*t30*t46*t9-1.5214d+0*t1*t10*t18*t19*
     5         t37*t39*t54+1.5214d+0*t1*t10*t19*t21*t37*t53+5.0713333333
     6         33333d-1*t1*t10*t18*t21*t37*t46
            t95 = -1.0d+0*rhoa*rhob*t94*wght-1.0d+0*rhob*t57*wght-1.0d+0
     1         *rhoa*t43*wght-1.0d+0*t27*wght
            t96 = -2.0d+0*rhoa*rhob*t27*t40*t56*wght-2.0d+0*rhoa*rhob*t2
     1         7*t42*t55*wght
            t97 = 2*t1*t10*t18*t19*t38*t54*t76*t9-3.333333333333333d-1*t
     1         1*t18*t19*t30*t39*t54*t9-t1*t10*t19*t38*t39*t53*t9+3.3333
     2         33333333333d-1*t1*t19*t21*t30*t53*t9-3.333333333333333d-1
     3         *t1*t10*t18*t38*t39*t46*t9+1.111111111111111d-1*t1*t18*t2
     4         1*t30*t46*t9-t1*t10*t18*t19*t37*t39*t54+t1*t10*t19*t21*t3
     5         7*t53+3.333333333333333d-1*t1*t10*t18*t21*t37*t46
            t98 = 2.0d+0*rhoa*rhob*t27*t97*wght+2.0d+0*rhoa*rhob*t40*t57
     1         *wght+2.0d+0*rhoa*t27*t40*wght-2.0d+0*t44*t55
            t99 = 2.0d+0*t1*t10*t18*t19*t21*t9*t97+2.2568000000000002d+0
     1         *t1*t10*t18*t19*t38*t54*t76*t9-2.0d+0*t1*t10*t18*t19*t39*
     2         t40*t54*t9-3.761333333333333d-1*t1*t18*t19*t30*t39*t54*t9
     3         +2.0d+0*t1*t10*t19*t21*t40*t53*t9-1.1284000000000001d+0*t
     4         1*t10*t19*t38*t39*t53*t9+3.761333333333333d-1*t1*t19*t21*
     5         t30*t53*t9+6.666666666666666d-1*t1*t10*t18*t21*t40*t46*t9
     6         -3.761333333333333d-1*t1*t10*t18*t38*t39*t46*t9+1.2537777
     7         777777775d-1*t1*t18*t21*t30*t46*t9-1.1284000000000001d+0*
     8         t1*t10*t18*t19*t37*t39*t54+1.1284000000000001d+0*t1*t10*t
     9         19*t21*t37*t53+3.761333333333333d-1*t1*t10*t18*t21*t37*t4
     :         6
            t100 = 1.0d+0*rhoa*rhob*t27*t99*wght+1.0d+0*rhoa*rhob*t42*t5
     1         7*wght+1.0d+0*rhoa*t27*t42*wght-t44*t56
            t101 = t56**2
            t102 = t55**2
            t103 = 1/rhob**1.6666666666666669d+0
            t104 = t54**2
            t105 = 1/rhob**4.666666666666667d+0
            t106 = t51**2
            t107 = 1/t15**3
            t108 = gammabb**2
            t109 = 1/t48**3
            t110 = 1/rhob**7.333333333333333d+0
            t111 = 1/rhob**3.3333333333333337d+0
            t112 = 1.6799999999999998d-1*gammabb*t105*t49+7.84d-2*t11*t1
     1         11*t14-4.4799999999999995d-2*t108*t109*t110
            t113 = 4.4799999999999995d-2*gammabb*t47*t51*t52-8.4d-3*gamm
     1         abb*t112*t17*t52+1.6799999999999998d-2*gammabb*t106*t107*
     2         t17+8.213333333333332d-2*gammabb*t105*t16
            t114 = 6.666666666666666d-1*t46*t53+t113*t19-2.2222222222222
     1         22d-1*t103*t18
            t115 = 3.0428d+0*t1*t10*t104*t18*t19*t76*t9-3.0428d+0*t1*t10
     1         *t19*t39*t53*t54*t9-1.0142666666666666d+0*t1*t10*t18*t39*
     2         t46*t54*t9+1.0142666666666666d+0*t1*t10*t21*t46*t53*t9-1.
     3         5214d+0*t1*t10*t114*t18*t19*t39*t9+1.5214d+0*t1*t10*t113*
     4         t19*t21*t9-3.380888888888889d-1*t1*t10*t103*t18*t21*t9
            t116 = -2.0d+0*rhoa*t57*wght-1.0d+0*rhoa*rhob*t115*wght
            t117 = 2*t1*t10*t104*t18*t19*t76*t9-2*t1*t10*t19*t39*t53*t54
     1         *t9-6.666666666666666d-1*t1*t10*t18*t39*t46*t54*t9+6.6666
     2         66666666666d-1*t1*t10*t21*t46*t53*t9-t1*t10*t114*t18*t19*
     3         t39*t9+t1*t10*t113*t19*t21*t9-2.222222222222222d-1*t1*t10
     4         *t103*t18*t21*t9
            t118 = 2.0d+0*rhoa*rhob*t55*t57*wght+2.0d+0*rhoa*t27*t55*wgh
     1         t+2.0d+0*rhoa*rhob*t117*t27*wght-2.0d+0*t55*t58
            t119 = 2.2568000000000002d+0*t1*t10*t104*t18*t19*t76*t9-2.0d
     1         +0*t1*t10*t18*t19*t39*t54*t55*t9+2.0d+0*t1*t10*t19*t21*t5
     2         3*t55*t9+6.666666666666666d-1*t1*t10*t18*t21*t46*t55*t9-2
     3         .2568000000000002d+0*t1*t10*t19*t39*t53*t54*t9-7.52266666
     4         6666666d-1*t1*t10*t18*t39*t46*t54*t9+7.522666666666666d-1
     5         *t1*t10*t21*t46*t53*t9-1.1284000000000001d+0*t1*t10*t114*
     6         t18*t19*t39*t9+2.0d+0*t1*t10*t117*t18*t19*t21*t9+1.128400
     7         0000000001d+0*t1*t10*t113*t19*t21*t9-2.507555555555555d-1
     8         *t1*t10*t103*t18*t21*t9
            t120 = 1.0d+0*rhoa*rhob*t56*t57*wght+1.0d+0*rhoa*t27*t56*wgh
     1         t+1.0d+0*rhoa*rhob*t119*t27*wght-t56*t58
            t121 = 1/t10
            t122 = 1/rhoa**6.333333333333333d+0
            t123 = 1.6799999999999998d-2*gammaaa*t122*t82-1.679999999999
     1         9998d-2*t34*t5*t59-5.04d-2*t31*t33
            t124 = 1.6799999999999998d-2*gammaaa*t35*t60*t8*t80-8.4d-3*t
     1         35*t36*t8-8.4d-3*gammaaa*t123*t36*t8-2.2399999999999998d-
     2         2*t31*t7+2.2399999999999998d-2*gammaaa*t31*t36*t60
            t125 = 3.333333333333333d-1*t30*t61+t10*t124
            t126 = 3.0428d+0*t1*t18*t19*t29*t38*t61*t76*t9-5.07133333333
     1         3333d-1*t1*t121*t18*t19*t39*t61*t9-1.5214d+0*t1*t10*t125*
     2         t18*t19*t39*t9-1.5214d+0*t1*t10*t18*t19*t38*t39*t61-1.521
     3         4d+0*t1*t18*t19*t29*t37*t39*t61+5.071333333333333d-1*t1*t
     4         18*t19*t21*t30*t61+1.5214d+0*t1*t10*t124*t18*t19*t21
            t127 = -1.0d+0*rhob*t62*wght-1.0d+0*rhoa*rhob*t126*wght
            t128 = -2.0d+0*rhoa*rhob*t27*t40*t64*wght-2.0d+0*rhoa*rhob*t
     1         27*t42*t63*wght
            t129 = 2*t1*t18*t19*t29*t38*t61*t76*t9-3.333333333333333d-1*
     1         t1*t121*t18*t19*t39*t61*t9-t1*t10*t125*t18*t19*t39*t9-t1*
     2         t10*t18*t19*t38*t39*t61-t1*t18*t19*t29*t37*t39*t61+3.3333
     3         33333333333d-1*t1*t18*t19*t21*t30*t61+t1*t10*t124*t18*t19
     4         *t21
            t130 = 2.0d+0*rhoa*rhob*t40*t62*wght+2.0d+0*rhoa*rhob*t129*t
     1         27*wght-2.0d+0*t44*t63
            t131 = 2.2568000000000002d+0*t1*t18*t19*t29*t38*t61*t76*t9-2
     1         .0d+0*t1*t18*t19*t29*t39*t40*t61*t9-3.761333333333333d-1*
     2         t1*t121*t18*t19*t39*t61*t9-1.1284000000000001d+0*t1*t10*t
     3         125*t18*t19*t39*t9+2.0d+0*t1*t10*t129*t18*t19*t21*t9+2.0d
     4         +0*t1*t10*t18*t19*t21*t40*t61-1.1284000000000001d+0*t1*t1
     5         0*t18*t19*t38*t39*t61-1.1284000000000001d+0*t1*t18*t19*t2
     6         9*t37*t39*t61+3.761333333333333d-1*t1*t18*t19*t21*t30*t61
     7         +1.1284000000000001d+0*t1*t10*t124*t18*t19*t21
            t132 = 1.0d+0*rhoa*rhob*t42*t62*wght+1.0d+0*rhoa*rhob*t131*t
     1         27*wght-t44*t64
            t133 = 3.0428d+0*t1*t10*t18*t38*t45*t67*t76*t9-5.07133333333
     1         3333d-1*t1*t18*t30*t39*t45*t67*t9-1.5214d+0*t1*t10*t19*t3
     2         8*t39*t67*t9+5.071333333333333d-1*t1*t19*t21*t30*t67*t9-1
     3         .5214d+0*t1*t10*t18*t37*t39*t45*t67+1.5214d+0*t1*t10*t19*
     4         t21*t37*t67
            t134 = -1.0d+0*rhob*t68*wght-1.0d+0*rhoa*rhob*t133*wght
            t135 = -2.0d+0*rhoa*rhob*t27*t40*t70*wght-2.0d+0*rhoa*rhob*t
     1         27*t42*t69*wght
            t136 = 2*t1*t10*t18*t38*t45*t67*t76*t9-3.333333333333333d-1*
     1         t1*t18*t30*t39*t45*t67*t9-t1*t10*t19*t38*t39*t67*t9+3.333
     2         333333333333d-1*t1*t19*t21*t30*t67*t9-t1*t10*t18*t37*t39*
     3         t45*t67+t1*t10*t19*t21*t37*t67
            t137 = 2.2568000000000002d+0*t1*t10*t18*t38*t45*t67*t76*t9-2
     1         .0d+0*t1*t10*t18*t39*t40*t45*t67*t9-3.761333333333333d-1*
     2         t1*t18*t30*t39*t45*t67*t9+2.0d+0*t1*t10*t19*t21*t40*t67*t
     3         9-1.1284000000000001d+0*t1*t10*t19*t38*t39*t67*t9+3.76133
     4         3333333333d-1*t1*t19*t21*t30*t67*t9+2.0d+0*t1*t10*t136*t1
     5         8*t19*t21*t9-1.1284000000000001d+0*t1*t10*t18*t37*t39*t45
     6         *t67+1.1284000000000001d+0*t1*t10*t19*t21*t37*t67
            t138 = 1.0d+0*rhoa*rhob*t42*t68*wght+1.0d+0*rhoa*rhob*t137*t
     1         27*wght-t44*t70
            t139 = 2.0d+0*rhoa*rhob*t40*t68*wght+2.0d+0*rhoa*rhob*t136*t
     1         27*wght-2.0d+0*t44*t69
            t140 = 3.0428d+0*t1*t18*t19*t29*t54*t61*t76*t9-1.5214d+0*t1*
     1         t19*t29*t39*t53*t61*t9-5.071333333333333d-1*t1*t18*t29*t3
     2         9*t46*t61*t9-1.5214d+0*t1*t10*t18*t19*t39*t54*t61+1.5214d
     3         +0*t1*t10*t19*t21*t53*t61+5.071333333333333d-1*t1*t10*t18
     4         *t21*t46*t61
            t141 = -1.0d+0*rhoa*t62*wght-1.0d+0*rhoa*rhob*t140*wght
            t142 = -2.0d+0*rhoa*rhob*t27*t55*t64*wght-2.0d+0*rhoa*rhob*t
     1         27*t56*t63*wght
            t143 = 2*t1*t18*t19*t29*t54*t61*t76*t9-t1*t19*t29*t39*t53*t6
     1         1*t9-3.333333333333333d-1*t1*t18*t29*t39*t46*t61*t9-t1*t1
     2         0*t18*t19*t39*t54*t61+t1*t10*t19*t21*t53*t61+3.3333333333
     3         33333d-1*t1*t10*t18*t21*t46*t61
            t144 = 2.0d+0*rhoa*rhob*t55*t62*wght+2.0d+0*rhoa*rhob*t143*t
     1         27*wght-2.0d+0*t58*t63
            t145 = 2.2568000000000002d+0*t1*t18*t19*t29*t54*t61*t76*t9-2
     1         .0d+0*t1*t18*t19*t29*t39*t55*t61*t9-1.1284000000000001d+0
     2         *t1*t19*t29*t39*t53*t61*t9-3.761333333333333d-1*t1*t18*t2
     3         9*t39*t46*t61*t9+2.0d+0*t1*t10*t143*t18*t19*t21*t9+2.0d+0
     4         *t1*t10*t18*t19*t21*t55*t61-1.1284000000000001d+0*t1*t10*
     5         t18*t19*t39*t54*t61+1.1284000000000001d+0*t1*t10*t19*t21*
     6         t53*t61+3.761333333333333d-1*t1*t10*t18*t21*t46*t61
            t146 = 1.0d+0*rhoa*rhob*t56*t62*wght+1.0d+0*rhoa*rhob*t145*t
     1         27*wght-t58*t64
            t147 = 1/t19
            t148 = 1/rhob**6.333333333333333d+0
            t149 = -1.6799999999999998d-2*t14*t50*t65-5.04d-2*t47*t49+1.
     1         6799999999999998d-2*gammabb*t109*t148
            t150 = 2.2399999999999998d-2*gammabb*t47*t52*t66+1.679999999
     1         9999998d-2*gammabb*t107*t17*t51*t66-8.4d-3*t17*t51*t52-8.
     2         4d-3*gammabb*t149*t17*t52-2.2399999999999998d-2*t16*t47
            t151 = 3.333333333333333d-1*t46*t67+t150*t19
            t152 = 3.0428d+0*t1*t10*t18*t45*t54*t67*t76*t9-1.5214d+0*t1*
     1         t10*t19*t39*t54*t67*t9-1.5214d+0*t1*t10*t39*t45*t53*t67*t
     2         9+5.071333333333333d-1*t1*t10*t21*t46*t67*t9-5.0713333333
     3         33333d-1*t1*t10*t147*t18*t39*t67*t9-1.5214d+0*t1*t10*t151
     4         *t18*t19*t39*t9+1.5214d+0*t1*t10*t150*t19*t21*t9
            t153 = -1.0d+0*rhoa*t68*wght-1.0d+0*rhoa*rhob*t152*wght
            t154 = -2.0d+0*rhoa*rhob*t27*t55*t70*wght-2.0d+0*rhoa*rhob*t
     1         27*t56*t69*wght
            t155 = 2*t1*t10*t18*t45*t54*t67*t76*t9-t1*t10*t19*t39*t54*t6
     1         7*t9-t1*t10*t39*t45*t53*t67*t9+3.333333333333333d-1*t1*t1
     2         0*t21*t46*t67*t9-3.333333333333333d-1*t1*t10*t147*t18*t39
     3         *t67*t9-t1*t10*t151*t18*t19*t39*t9+t1*t10*t150*t19*t21*t9
            t156 = 2.2568000000000002d+0*t1*t10*t18*t45*t54*t67*t76*t9-2
     1         .0d+0*t1*t10*t18*t39*t45*t55*t67*t9+2.0d+0*t1*t10*t19*t21
     2         *t55*t67*t9-1.1284000000000001d+0*t1*t10*t19*t39*t54*t67*
     3         t9-1.1284000000000001d+0*t1*t10*t39*t45*t53*t67*t9+3.7613
     4         33333333333d-1*t1*t10*t21*t46*t67*t9-3.761333333333333d-1
     5         *t1*t10*t147*t18*t39*t67*t9-1.1284000000000001d+0*t1*t10*
     6         t151*t18*t19*t39*t9+2.0d+0*t1*t10*t155*t18*t19*t21*t9+1.1
     7         284000000000001d+0*t1*t10*t150*t19*t21*t9
            t157 = 1.0d+0*rhoa*rhob*t56*t68*wght+1.0d+0*rhoa*rhob*t156*t
     1         27*wght-t58*t70
            t158 = 2.0d+0*rhoa*rhob*t55*t68*wght+2.0d+0*rhoa*rhob*t155*t
     1         27*wght-2.0d+0*t58*t69
            t159 = t61**2
            t160 = t60**2
            t161 = 1/rhoa**5.333333333333333d+0
            t162 = 1/gammaaa
            t163 = 1/t2**3
            t164 = -6.3d-3*t161*t82+6.3d-3*t162*t33*t8-6.3d-3*t163*t4*t5
            t165 = 1.6799999999999998d-2*gammaaa*t160*t8*t80-1.679999999
     1         9999998d-2*t36*t60*t8-8.4d-3*gammaaa*t164*t36*t8
            t166 = 3.0428d+0*rhoa*t1*t159*t18*t19*t76*t9-1.5214d+0*t1*t1
     1         65*t18*t19*t29*t39*t9-3.0428d+0*t1*t159*t18*t19*t29*t39+1
     2         .5214d+0*t1*t10*t165*t18*t19*t21
            t167 = t64**2
            t168 = t63**2
            t169 = 2*rhoa*t1*t159*t18*t19*t76*t9-t1*t165*t18*t19*t29*t39
     1         *t9-2*t1*t159*t18*t19*t29*t39+t1*t10*t165*t18*t19*t21
            t170 = 4.0d+0*rhoa*rhob*t62*t63*wght+2.0d+0*rhoa*rhob*t169*t
     1         27*wght
            t171 = 2.2568000000000002d+0*rhoa*t1*t159*t18*t19*t76*t9-2.0
     1         d+0*t1*t18*t19*t29*t39*t61*t63*t9-1.1284000000000001d+0*t
     2         1*t165*t18*t19*t29*t39*t9+2.0d+0*t1*t10*t169*t18*t19*t21*
     3         t9+2.0d+0*t1*t10*t18*t19*t21*t61*t63-2.2568000000000002d+
     4         0*t1*t159*t18*t19*t29*t39+1.1284000000000001d+0*t1*t10*t1
     5         65*t18*t19*t21
            t172 = 2.0d+0*rhoa*rhob*t62*t64*wght+1.0d+0*rhoa*rhob*t171*t
     1         27*wght
            t173 = 3.0428d+0*t1*t18*t29*t45*t61*t67*t76*t9-1.5214d+0*t1*
     1         t19*t29*t39*t61*t67*t9-1.5214d+0*t1*t10*t18*t39*t45*t61*t
     2         67+1.5214d+0*t1*t10*t19*t21*t61*t67
            t174 = 2*t1*t18*t29*t45*t61*t67*t76*t9-t1*t19*t29*t39*t61*t6
     1         7*t9-t1*t10*t18*t39*t45*t61*t67+t1*t10*t19*t21*t61*t67
            t175 = 2.0d+0*rhoa*rhob*t62*t69*wght+2.0d+0*rhoa*rhob*t63*t6
     1         8*wght+2.0d+0*rhoa*rhob*t174*t27*wght
            t176 = -2.0d+0*rhoa*rhob*t27*t63*t70*wght-2.0d+0*rhoa*rhob*t
     1         27*t64*t69*wght
            t177 = 2.2568000000000002d+0*t1*t18*t29*t45*t61*t67*t76*t9-2
     1         .0d+0*t1*t10*t18*t39*t45*t63*t67*t9+2.0d+0*t1*t10*t19*t21
     2         *t63*t67*t9-1.1284000000000001d+0*t1*t19*t29*t39*t61*t67*
     3         t9+2.0d+0*t1*t10*t174*t18*t19*t21*t9-1.1284000000000001d+
     4         0*t1*t10*t18*t39*t45*t61*t67+1.1284000000000001d+0*t1*t10
     5         *t19*t21*t61*t67
            t178 = 1.0d+0*rhoa*rhob*t62*t70*wght+1.0d+0*rhoa*rhob*t64*t6
     1         8*wght+1.0d+0*rhoa*rhob*t177*t27*wght
            t179 = t70**2
            t180 = t69**2
            t181 = t66**2
            t182 = 1/rhob**5.333333333333333d+0
            t183 = 1/gammabb
            t184 = 1/t11**3
            t185 = 6.3d-3*t17*t183*t49-6.3d-3*t13*t14*t184-6.3d-3*t109*t
     1         182
            t186 = -1.6799999999999998d-2*t17*t52*t66-8.4d-3*gammabb*t17
     1         *t185*t52+1.6799999999999998d-2*gammabb*t107*t17*t181
            t187 = t67**2
            t188 = 3.0428d+0*rhob*t1*t10*t18*t187*t76*t9-3.0428d+0*t1*t1
     1         0*t187*t39*t45*t9-1.5214d+0*t1*t10*t18*t186*t39*t45*t9+1.
     2         5214d+0*t1*t10*t186*t19*t21*t9
            t189 = 2*rhob*t1*t10*t18*t187*t76*t9-2*t1*t10*t187*t39*t45*t
     1         9-t1*t10*t18*t186*t39*t45*t9+t1*t10*t186*t19*t21*t9
            t190 = 4.0d+0*rhoa*rhob*t68*t69*wght+2.0d+0*rhoa*rhob*t189*t
     1         27*wght
            t191 = 2.2568000000000002d+0*rhob*t1*t10*t18*t187*t76*t9-2.0
     1         d+0*t1*t10*t18*t39*t45*t67*t69*t9+2.0d+0*t1*t10*t19*t21*t
     2         67*t69*t9-2.2568000000000002d+0*t1*t10*t187*t39*t45*t9-1.
     3         1284000000000001d+0*t1*t10*t18*t186*t39*t45*t9+2.0d+0*t1*
     4         t10*t18*t189*t19*t21*t9+1.1284000000000001d+0*t1*t10*t186
     5         *t19*t21*t9
            t192 = 2.0d+0*rhoa*rhob*t68*t70*wght+1.0d+0*rhoa*rhob*t191*t
     1         27*wght
            t193 = 1/t25**4
            t194 = 1/t22**5.0d+0
            t195 = t38**3
            t196 = 1/t20**4
            t197 = 1/rhoa**5.666666666666667d+0
            t198 = 1/t6**4
            t199 = 1/t32**5
            t200 = 5.04d-2*gammaaa*t35*t8*t80*t85+6.719999999999998d-2*g
     1         ammaaa*t31*t36*t85-8.4d-3*gammaaa*t36*t8*(5.5253333333333
     2         33d-1*t81*t82/rhoa**8.333333333333334d+0-2.61333333333333
     3         36d-1*t2*t5/rhoa**4.333333333333333d+0-8.885333333333333d
     4         -1*gammaaa*t197*t33-1.7919999999999997d-1*gammaaa**3*t199
     5         /rhoa**11)-1.3439999999999996d-1*gammaaa*t31*t79*t80-5.04
     6         d-2*gammaaa*t198*t35**3*t8-2.4639999999999995d-1*gammaaa*
     7         t35*t36*t78-3.8328888888888885d-1*gammaaa*t197*t7
            t201 = 3.7037037037037035d-1*t8*t9+t30*t86-6.666666666666666
     1         d-1*t37*t77+t10*t200
            t202 = 6*t1*t10*t18*t19*t38*t76*t87*t9-t1*t18*t19*t30*t39*t8
     1         7*t9+3.7037037037037035d-1*t1*t18*t19*t21*t8*t9+6.6666666
     2         66666666d-1*t1*t18*t19*t38*t39*t77*t9+2*t1*t18*t19*t30*t7
     3         5*t76*t9-t1*t10*t18*t19*t201*t39*t9-6*t1*t10*t18*t19*t195
     4         *t196*t9-3*t1*t10*t18*t19*t37*t39*t87-3*t1*t10*t18*t19*t3
     5         8*t39*t86+t1*t18*t19*t21*t30*t86-6.666666666666666d-1*t1*
     6         t18*t19*t21*t37*t77+6*t1*t10*t18*t19*t37*t75*t76-2*t1*t18
     7         *t19*t30*t37*t38*t39+t1*t10*t18*t19*t200*t21
            t203 = 2*t1*t10*t18*t19*t54*t76*t87*t9-t1*t10*t19*t39*t53*t8
     1         7*t9-3.333333333333333d-1*t1*t10*t18*t39*t46*t87*t9+2.222
     2         222222222222d-1*t1*t18*t19*t39*t54*t77*t9-2.2222222222222
     3         22d-1*t1*t19*t21*t53*t77*t9-7.407407407407407d-2*t1*t18*t
     4         21*t46*t77*t9+2*t1*t10*t19*t53*t75*t76*t9+6.6666666666666
     5         66d-1*t1*t10*t18*t46*t75*t76*t9+1.3333333333333333d+0*t1*
     6         t18*t19*t30*t38*t54*t76*t9-6*t1*t10*t18*t19*t196*t54*t75*
     7         t9-6.666666666666666d-1*t1*t19*t30*t38*t39*t53*t9-2.22222
     8         2222222222d-1*t1*t18*t30*t38*t39*t46*t9-t1*t10*t18*t19*t3
     9         9*t54*t86+t1*t10*t19*t21*t53*t86+3.333333333333333d-1*t1*
     :         t10*t18*t21*t46*t86+4*t1*t10*t18*t19*t37*t38*t54*t76-6.66
     ;         6666666666666d-1*t1*t18*t19*t30*t37*t39*t54-2*t1*t10*t19*
     <         t37*t38*t39*t53+6.666666666666666d-1*t1*t19*t21*t30*t37*t
     =         53-6.666666666666666d-1*t1*t10*t18*t37*t38*t39*t46+2.2222
     >         22222222222d-1*t1*t18*t21*t30*t37*t46
            t204 = 4*t1*t10*t19*t38*t53*t54*t76*t9+1.3333333333333333d+0
     1         *t1*t10*t18*t38*t46*t54*t76*t9+2*t1*t10*t114*t18*t19*t38*
     2         t76*t9+6.666666666666666d-1*t1*t104*t18*t19*t30*t76*t9-6.
     3         666666666666666d-1*t1*t19*t30*t39*t53*t54*t9-2.2222222222
     4         22222d-1*t1*t18*t30*t39*t46*t54*t9-6.666666666666666d-1*t
     5         1*t10*t38*t39*t46*t53*t9+2.222222222222222d-1*t1*t21*t30*
     6         t46*t53*t9-t1*t10*t113*t19*t38*t39*t9+2.222222222222222d-
     7         1*t1*t10*t103*t18*t38*t39*t9-3.333333333333333d-1*t1*t114
     8         *t18*t19*t30*t39*t9-6*t1*t10*t104*t18*t19*t196*t38*t9+3.3
     9         33333333333333d-1*t1*t113*t19*t21*t30*t9-7.40740740740740
     :         7d-2*t1*t103*t18*t21*t30*t9+2*t1*t10*t104*t18*t19*t37*t76
     ;         -2*t1*t10*t19*t37*t39*t53*t54-6.666666666666666d-1*t1*t10
     <         *t18*t37*t39*t46*t54+6.666666666666666d-1*t1*t10*t21*t37*
     =         t46*t53-t1*t10*t114*t18*t19*t37*t39+t1*t10*t113*t19*t21*t
     >         37-2.222222222222222d-1*t1*t10*t103*t18*t21*t37
            t205 = t54**3
            t206 = 1/rhob**5.666666666666667d+0
            t207 = 1/t15**4
            t208 = 1/t48**5
            t209 = -2.4639999999999995d-1*gammabb*t105*t51*t52-8.4d-3*ga
     1         mmabb*t17*(-8.885333333333333d-1*gammabb*t206*t49-1.79199
     2         99999999997d-1*gammabb**3*t208/rhob**11-2.613333333333333
     3         6d-1*t11*t14/rhob**4.333333333333333d+0+5.525333333333333
     4         d-1*t108*t109/rhob**8.333333333333334d+0)*t52+6.719999999
     5         999998d-2*gammabb*t112*t47*t52-5.04d-2*gammabb*t17*t207*t
     6         51**3+5.04d-2*gammabb*t107*t112*t17*t51-1.343999999999999
     7         6d-1*gammabb*t106*t107*t47-3.8328888888888885d-1*gammabb*
     8         t16*t206
            t210 = -6.666666666666666d-1*t103*t53+t113*t46+t19*t209+3.70
     1         37037037037035d-1*t17*t18
            t211 = 6*t1*t10*t114*t18*t19*t54*t76*t9+6*t1*t10*t104*t19*t5
     1         3*t76*t9+2*t1*t10*t104*t18*t46*t76*t9-2*t1*t10*t39*t46*t5
     2         3*t54*t9-3*t1*t10*t113*t19*t39*t54*t9+6.666666666666666d-
     3         1*t1*t10*t103*t18*t39*t54*t9-3*t1*t10*t114*t19*t39*t53*t9
     4         -6.666666666666666d-1*t1*t10*t103*t21*t53*t9-t1*t10*t114*
     5         t18*t39*t46*t9+t1*t10*t113*t21*t46*t9-t1*t10*t18*t19*t210
     6         *t39*t9+t1*t10*t19*t209*t21*t9+3.7037037037037035d-1*t1*t
     7         10*t17*t18*t21*t9-6*t1*t10*t18*t19*t196*t205*t9
            t212 = 1.6799999999999998d-2*gammaaa*t60*t8*t80*t85-8.4d-3*t
     1         36*t8*t85-8.4d-3*gammaaa*t36*t8*(3.92d-2*t5*t59*t84-1.735
     2         9999999999998d-1*gammaaa*t82*t83+6.719999999999998d-2*t19
     3         9*t81/rhoa**10+2.072d-1*t33*t78)+1.6799999999999998d-2*t7
     4         9*t8*t80+3.3599999999999997d-2*gammaaa*t123*t35*t8*t80-8.
     5         959999999999999d-2*gammaaa*t31*t35*t60*t80-5.04d-2*gammaa
     6         a*t198*t60*t79*t8+8.213333333333332d-2*t7*t78-8.213333333
     7         333332d-2*gammaaa*t36*t60*t78+4.4799999999999995d-2*t31*t
     8         35*t36+4.4799999999999995d-2*gammaaa*t123*t31*t36
            t213 = -2.222222222222222d-1*t61*t77+6.666666666666666d-1*t1
     1         24*t30+t10*t212
            t214 = 2*t1*t18*t19*t29*t61*t76*t87*t9+1.3333333333333333d+0
     1         *t1*t121*t18*t19*t38*t61*t76*t9+4*t1*t10*t125*t18*t19*t38
     2         *t76*t9-6*t1*t18*t19*t196*t29*t61*t75*t9+2.22222222222222
     3         2d-1*t1*t18*t19*t39*t4*t61*t9-6.666666666666666d-1*t1*t12
     4         5*t18*t19*t30*t39*t9-t1*t10*t18*t19*t213*t39*t9-t1*t10*t1
     5         8*t19*t39*t61*t87-t1*t18*t19*t29*t39*t61*t86-2.2222222222
     6         22222d-1*t1*t18*t19*t21*t61*t77+2*t1*t10*t18*t19*t61*t75*
     7         t76+4*t1*t18*t19*t29*t37*t38*t61*t76-6.666666666666666d-1
     8         *t1*t18*t19*t30*t38*t39*t61-6.666666666666666d-1*t1*t121*
     9         t18*t19*t37*t39*t61-2*t1*t10*t124*t18*t19*t38*t39-2*t1*t1
     :         0*t125*t18*t19*t37*t39+6.666666666666666d-1*t1*t124*t18*t
     ;         19*t21*t30+t1*t10*t18*t19*t21*t212
            t215 = 2*t1*t10*t18*t45*t67*t76*t87*t9-t1*t10*t19*t39*t67*t8
     1         7*t9+2.222222222222222d-1*t1*t18*t39*t45*t67*t77*t9-2.222
     2         222222222222d-1*t1*t19*t21*t67*t77*t9+2*t1*t10*t19*t67*t7
     3         5*t76*t9+1.3333333333333333d+0*t1*t18*t30*t38*t45*t67*t76
     4         *t9-6*t1*t10*t18*t196*t45*t67*t75*t9-6.666666666666666d-1
     5         *t1*t19*t30*t38*t39*t67*t9-t1*t10*t18*t39*t45*t67*t86+t1*
     6         t10*t19*t21*t67*t86+4*t1*t10*t18*t37*t38*t45*t67*t76-6.66
     7         6666666666666d-1*t1*t18*t30*t37*t39*t45*t67-2*t1*t10*t19*
     8         t37*t38*t39*t67+6.666666666666666d-1*t1*t19*t21*t30*t37*t
     9         67
            t216 = 6.666666666666666d-1*t1*t121*t18*t19*t54*t61*t76*t9+2
     1         *t1*t19*t29*t38*t53*t61*t76*t9+6.666666666666666d-1*t1*t1
     2         8*t29*t38*t46*t61*t76*t9+2*t1*t10*t125*t18*t19*t54*t76*t9
     3         -6*t1*t18*t19*t196*t29*t38*t54*t61*t9-3.333333333333333d-
     4         1*t1*t121*t19*t39*t53*t61*t9-1.111111111111111d-1*t1*t121
     5         *t18*t39*t46*t61*t9-t1*t10*t125*t19*t39*t53*t9-3.33333333
     6         3333333d-1*t1*t10*t125*t18*t39*t46*t9+2*t1*t10*t18*t19*t3
     7         8*t54*t61*t76+2*t1*t18*t19*t29*t37*t54*t61*t76-3.33333333
     8         3333333d-1*t1*t18*t19*t30*t39*t54*t61-t1*t10*t19*t38*t39*
     9         t53*t61-t1*t19*t29*t37*t39*t53*t61+3.333333333333333d-1*t
     :         1*t19*t21*t30*t53*t61-3.333333333333333d-1*t1*t10*t18*t38
     ;         *t39*t46*t61-3.333333333333333d-1*t1*t18*t29*t37*t39*t46*
     <         t61+1.111111111111111d-1*t1*t18*t21*t30*t46*t61-t1*t10*t1
     =         24*t18*t19*t39*t54+t1*t10*t124*t19*t21*t53+3.333333333333
     >         333d-1*t1*t10*t124*t18*t21*t46
            t217 = 6.666666666666666d-1*t1*t18*t30*t45*t54*t67*t76*t9+2*
     1         t1*t10*t19*t38*t54*t67*t76*t9+2*t1*t10*t38*t45*t53*t67*t7
     2         6*t9+6.666666666666666d-1*t1*t10*t147*t18*t38*t67*t76*t9+
     3         2*t1*t10*t151*t18*t19*t38*t76*t9-6*t1*t10*t18*t196*t38*t4
     4         5*t54*t67*t9-3.333333333333333d-1*t1*t19*t30*t39*t54*t67*
     5         t9-3.333333333333333d-1*t1*t30*t39*t45*t53*t67*t9-3.33333
     6         3333333333d-1*t1*t10*t38*t39*t46*t67*t9+1.111111111111111
     7         d-1*t1*t21*t30*t46*t67*t9-1.111111111111111d-1*t1*t147*t1
     8         8*t30*t39*t67*t9-t1*t10*t150*t19*t38*t39*t9-3.33333333333
     9         3333d-1*t1*t151*t18*t19*t30*t39*t9+3.333333333333333d-1*t
     :         1*t150*t19*t21*t30*t9+2*t1*t10*t18*t37*t45*t54*t67*t76-t1
     ;         *t10*t19*t37*t39*t54*t67-t1*t10*t37*t39*t45*t53*t67+3.333
     <         333333333333d-1*t1*t10*t21*t37*t46*t67-3.333333333333333d
     =         -1*t1*t10*t147*t18*t37*t39*t67-t1*t10*t151*t18*t19*t37*t3
     >         9+t1*t10*t150*t19*t21*t37
            t218 = 4*t1*t19*t29*t53*t54*t61*t76*t9+1.3333333333333333d+0
     1         *t1*t18*t29*t46*t54*t61*t76*t9+2*t1*t114*t18*t19*t29*t61*
     2         t76*t9-6.666666666666666d-1*t1*t29*t39*t46*t53*t61*t9-t1*
     3         t113*t19*t29*t39*t61*t9+2.222222222222222d-1*t1*t103*t18*
     4         t29*t39*t61*t9-6*t1*t104*t18*t19*t196*t29*t61*t9+2*t1*t10
     5         *t104*t18*t19*t61*t76-2*t1*t10*t19*t39*t53*t54*t61-6.6666
     6         66666666666d-1*t1*t10*t18*t39*t46*t54*t61+6.6666666666666
     7         66d-1*t1*t10*t21*t46*t53*t61-t1*t10*t114*t18*t19*t39*t61+
     8         t1*t10*t113*t19*t21*t61-2.222222222222222d-1*t1*t10*t103*
     9         t18*t21*t61
            t219 = -8.213333333333332d-2*gammabb*t105*t52*t66-8.95999999
     1         9999999d-2*gammabb*t107*t47*t51*t66-5.04d-2*gammabb*t106*
     2         t17*t207*t66+1.6799999999999998d-2*gammabb*t107*t112*t17*
     3         t66-8.4d-3*gammabb*t17*t52*(3.92d-2*t111*t14*t65+2.072d-1
     4         *t105*t49+6.719999999999998d-2*t108*t208/rhob**10-1.73599
     5         99999999998d-1*gammabb*t109*t110)+4.4799999999999995d-2*t
     6         47*t51*t52+4.4799999999999995d-2*gammabb*t149*t47*t52-8.4
     7         d-3*t112*t17*t52+3.3599999999999997d-2*gammabb*t107*t149*
     8         t17*t51+1.6799999999999998d-2*t106*t107*t17+8.21333333333
     9         3332d-2*t105*t16
            t220 = -2.222222222222222d-1*t103*t67+6.666666666666666d-1*t
     1         150*t46+t19*t219
            t221 = 4*t1*t10*t45*t53*t54*t67*t76*t9+1.3333333333333333d+0
     1         *t1*t10*t147*t18*t54*t67*t76*t9+2*t1*t10*t114*t18*t45*t67
     2         *t76*t9+2*t1*t10*t104*t19*t67*t76*t9+4*t1*t10*t151*t18*t1
     3         9*t54*t76*t9-6.666666666666666d-1*t1*t10*t39*t46*t54*t67*
     4         t9-6.666666666666666d-1*t1*t10*t147*t39*t53*t67*t9-t1*t10
     5         *t113*t39*t45*t67*t9-6*t1*t10*t104*t18*t196*t45*t67*t9-t1
     6         *t10*t114*t19*t39*t67*t9+2.222222222222222d-1*t1*t10*t13*
     7         t18*t39*t67*t9-2.222222222222222d-1*t1*t10*t103*t21*t67*t
     8         9-2*t1*t10*t150*t19*t39*t54*t9-2*t1*t10*t151*t19*t39*t53*
     9         t9-6.666666666666666d-1*t1*t10*t151*t18*t39*t46*t9+6.6666
     :         66666666666d-1*t1*t10*t150*t21*t46*t9-t1*t10*t18*t19*t220
     ;         *t39*t9+t1*t10*t19*t21*t219*t9
            t222 = -8.4d-3*gammaaa*t36*t8*(4.1999999999999993d-2*t122*t8
     1         2+8.4d-3*t163*t34*t5-8.4d-3*t162*t31*t33-2.52d-2*gammaaa*
     2         t199/rhoa**9)+3.3599999999999997d-2*t35*t60*t8*t80+3.3599
     3         999999999997d-2*gammaaa*t123*t60*t8*t80+1.679999999999999
     4         8d-2*gammaaa*t164*t35*t8*t80-4.4799999999999995d-2*gammaa
     5         a*t160*t31*t80-1.6799999999999998d-2*t123*t36*t8-5.04d-2*
     6         gammaaa*t160*t198*t35*t8+4.4799999999999995d-2*t31*t36*t6
     7         0+2.2399999999999998d-2*gammaaa*t164*t31*t36
            t223 = 3.333333333333333d-1*t165*t30+t10*t222
            t224 = 4*t1*t125*t18*t19*t29*t61*t76*t9+2*t1*t165*t18*t19*t2
     1         9*t38*t76*t9+6.666666666666666d-1*t1*t159*t18*t19*t76*t9-
     2         t1*t10*t18*t19*t223*t39*t9-3.333333333333333d-1*t1*t121*t
     3         165*t18*t19*t39*t9-6*rhoa*t1*t159*t18*t19*t196*t38*t9+4*t
     4         1*t159*t18*t19*t29*t38*t76+2*rhoa*t1*t159*t18*t19*t37*t76
     5         -2*t1*t124*t18*t19*t29*t39*t61-2*t1*t10*t125*t18*t19*t39*
     6         t61-t1*t10*t165*t18*t19*t38*t39-t1*t165*t18*t19*t29*t37*t
     7         39-6.666666666666666d-1*t1*t121*t159*t18*t19*t39+3.333333
     8         333333333d-1*t1*t165*t18*t19*t21*t30+t1*t10*t18*t19*t21*t
     9         222
            t225 = 6.666666666666666d-1*t1*t121*t18*t45*t61*t67*t76*t9+2
     1         *t1*t19*t29*t38*t61*t67*t76*t9+2*t1*t10*t125*t18*t45*t67*
     2         t76*t9-6*t1*t18*t196*t29*t38*t45*t61*t67*t9-3.33333333333
     3         3333d-1*t1*t121*t19*t39*t61*t67*t9-t1*t10*t125*t19*t39*t6
     4         7*t9+2*t1*t10*t18*t38*t45*t61*t67*t76+2*t1*t18*t29*t37*t4
     5         5*t61*t67*t76-3.333333333333333d-1*t1*t18*t30*t39*t45*t61
     6         *t67-t1*t10*t19*t38*t39*t61*t67-t1*t19*t29*t37*t39*t61*t6
     7         7+3.333333333333333d-1*t1*t19*t21*t30*t61*t67-t1*t10*t124
     8         *t18*t39*t45*t67+t1*t10*t124*t19*t21*t67
            t226 = 4*t1*t10*t187*t38*t45*t76*t9+2*t1*t10*t18*t186*t38*t4
     1         5*t76*t9+6.666666666666666d-1*rhob*t1*t18*t187*t30*t76*t9
     2         -6.666666666666666d-1*t1*t187*t30*t39*t45*t9-3.3333333333
     3         33333d-1*t1*t18*t186*t30*t39*t45*t9-t1*t10*t186*t19*t38*t
     4         39*t9-6*rhob*t1*t10*t18*t187*t196*t38*t9+3.33333333333333
     5         3d-1*t1*t186*t19*t21*t30*t9+2*rhob*t1*t10*t18*t187*t37*t7
     6         6-2*t1*t10*t187*t37*t39*t45-t1*t10*t18*t186*t37*t39*t45+t
     7         1*t10*t186*t19*t21*t37
            t227 = 2*t1*t165*t18*t19*t29*t54*t76*t9+2*rhoa*t1*t159*t19*t
     1         53*t76*t9+6.666666666666666d-1*rhoa*t1*t159*t18*t46*t76*t
     2         9-6*rhoa*t1*t159*t18*t19*t196*t54*t9-t1*t165*t19*t29*t39*
     3         t53*t9-3.333333333333333d-1*t1*t165*t18*t29*t39*t46*t9+4*
     4         t1*t159*t18*t19*t29*t54*t76-t1*t10*t165*t18*t19*t39*t54-2
     5         *t1*t159*t19*t29*t39*t53+t1*t10*t165*t19*t21*t53-6.666666
     6         666666666d-1*t1*t159*t18*t29*t39*t46+3.333333333333333d-1
     7         *t1*t10*t165*t18*t21*t46
            t228 = 2*t1*t19*t29*t54*t61*t67*t76*t9+2*t1*t29*t45*t53*t61*
     1         t67*t76*t9+6.666666666666666d-1*t1*t147*t18*t29*t61*t67*t
     2         76*t9+2*t1*t151*t18*t19*t29*t61*t76*t9-6*t1*t18*t196*t29*
     3         t45*t54*t61*t67*t9-3.333333333333333d-1*t1*t29*t39*t46*t6
     4         1*t67*t9-t1*t150*t19*t29*t39*t61*t9+2*t1*t10*t18*t45*t54*
     5         t61*t67*t76-t1*t10*t19*t39*t54*t61*t67-t1*t10*t39*t45*t53
     6         *t61*t67+3.333333333333333d-1*t1*t10*t21*t46*t61*t67-3.33
     7         3333333333333d-1*t1*t10*t147*t18*t39*t61*t67-t1*t10*t151*
     8         t18*t19*t39*t61+t1*t10*t150*t19*t21*t61
            t229 = 4.4799999999999995d-2*t47*t52*t66+3.3599999999999997d
     1         -2*t107*t17*t51*t66+3.3599999999999997d-2*gammabb*t107*t1
     2         49*t17*t66-8.4d-3*gammabb*t17*(8.4d-3*t14*t184*t50-8.4d-3
     3         *t183*t47*t49-2.52d-2*gammabb*t208/rhob**9+4.199999999999
     4         9993d-2*t109*t148)*t52+2.2399999999999998d-2*gammabb*t185
     5         *t47*t52-1.6799999999999998d-2*t149*t17*t52-5.04d-2*gamma
     6         bb*t17*t181*t207*t51+1.6799999999999998d-2*gammabb*t107*t
     7         17*t185*t51-4.4799999999999995d-2*gammabb*t107*t181*t47
            t230 = 3.333333333333333d-1*t186*t46+t19*t229
            t231 = 4*t1*t10*t151*t18*t45*t67*t76*t9+4*t1*t10*t187*t45*t5
     1         4*t76*t9+2*t1*t10*t18*t186*t45*t54*t76*t9+2*rhob*t1*t10*t
     2         187*t53*t76*t9+6.666666666666666d-1*t1*t10*t18*t187*t76*t
     3         9-2*t1*t10*t150*t39*t45*t67*t9-2*t1*t10*t151*t19*t39*t67*
     4         t9-t1*t10*t186*t19*t39*t54*t9-6*rhob*t1*t10*t18*t187*t196
     5         *t54*t9-t1*t10*t186*t39*t45*t53*t9+3.333333333333333d-1*t
     6         1*t10*t186*t21*t46*t9-t1*t10*t18*t19*t230*t39*t9-6.666666
     7         666666666d-1*t1*t10*t147*t187*t39*t9-3.333333333333333d-1
     8         *t1*t10*t147*t18*t186*t39*t9+t1*t10*t19*t21*t229*t9
            t232 = t61**3
            t233 = -8.4d-3*gammaaa*t36*t8*(-3.15d-3*t161*t162*t82-9.4500
     1         00000000001d-3*t33*t8/t81+9.450000000000001d-3*t4*t5/t2**
     2         5+9.450000000000001d-3*t199/rhoa**8)+5.04d-2*gammaaa*t164
     3         *t60*t8*t80+5.04d-2*t160*t8*t80-5.04d-2*gammaaa*t198*t60*
     4         *3*t8-2.52d-2*t164*t36*t8
            t234 = 6*rhoa*t1*t165*t18*t19*t61*t76*t9-t1*t18*t19*t233*t29
     1         *t39*t9-6*t1*t18*t19*t196*t232*t3*t9+6*rhoa*t1*t18*t19*t2
     2         32*t76-6*t1*t165*t18*t19*t29*t39*t61+t1*t10*t18*t19*t21*t
     3         233
            t235 = 2*t1*t165*t18*t29*t45*t67*t76*t9+2*rhoa*t1*t159*t19*t
     1         67*t76*t9-6*rhoa*t1*t159*t18*t196*t45*t67*t9-t1*t165*t19*
     2         t29*t39*t67*t9+4*t1*t159*t18*t29*t45*t67*t76-t1*t10*t165*
     3         t18*t39*t45*t67-2*t1*t159*t19*t29*t39*t67+t1*t10*t165*t19
     4         *t21*t67
            t236 = 4*t1*t187*t29*t45*t61*t76*t9+2*t1*t18*t186*t29*t45*t6
     1         1*t76*t9-t1*t186*t19*t29*t39*t61*t9-6*rhob*t1*t18*t187*t1
     2         96*t29*t61*t9+2*rhob*t1*t10*t18*t187*t61*t76-2*t1*t10*t18
     3         7*t39*t45*t61-t1*t10*t18*t186*t39*t45*t61+t1*t10*t186*t19
     4         *t21*t61
            t237 = -5.04d-2*gammabb*t17*t207*t66**3+5.04d-2*gammabb*t107
     1         *t17*t185*t66-8.4d-3*gammabb*t17*(-9.450000000000001d-3*t
     2         17*t49/t108+9.450000000000001d-3*t208/rhob**8-3.15d-3*t10
     3         9*t182*t183+9.450000000000001d-3*t13*t14/t11**5)*t52-2.52
     4         d-2*t17*t185*t52+5.04d-2*t107*t17*t181
            t238 = t67**3
            t239 = 6*rhob*t1*t10*t18*t186*t67*t76*t9+6*rhob*t1*t10*t238*
     1         t76*t9-6*t1*t10*t186*t39*t45*t67*t9-t1*t10*t18*t237*t39*t
     2         45*t9-6*t1*t10*t12*t18*t196*t238*t9+t1*t10*t19*t21*t237*t
     3         9
            fnc(iq) = fnc(iq)-1.0d+0*rhoa*rhob*t24*t26*t27*wght
            Amat(iq,D1_RA) = 1.0d+0*rhoa*rhob*t24*t27*t41*t42*wght+2.0d+
     1         0*rhoa*rhob*t26*t27*t28*t40*wght+t24*t26*t44+Amat(iq,D1_R
     2         A)
            Amat(iq,D1_RB) = 1.0d+0*rhoa*rhob*t24*t27*t41*t56*wght+2.0d+
     1         0*rhoa*rhob*t26*t27*t28*t55*wght+t24*t26*t58+Amat(iq,D1_R
     2         B)
            Cmat(iq,D1_GAA) = 1.0d+0*rhoa*rhob*t24*t27*t41*t64*wght+2.0d
     1         +0*rhoa*rhob*t26*t27*t28*t63*wght-1.0d+0*rhoa*rhob*t24*t2
     2         6*t62*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = 1.0d+0*rhoa*rhob*t24*t27*t41*t70*wght+2.0d
     1         +0*rhoa*rhob*t26*t27*t28*t69*wght-1.0d+0*rhoa*rhob*t24*t2
     2         6*t68*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = -6.0d+0*rhoa*rhob*t26*t27*t73*t74*wght-
     1         2.0d+0*rhoa*rhob*t24*t27*t71*t72*wght-4.0d+0*rhoa*rhob*t2
     2         7*t28*t40*t41*t42*wght+t24*t41*t93+t26*t28*t91+t24*t26*t8
     3         9+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = -6.0d+0*rhoa*rhob*t26*t27*t40*t55*t73*w
     1         ght-2.0d+0*rhoa*rhob*t24*t27*t42*t56*t71*wght+t26*t28*t98
     2         +t28*t41*t96+t24*t26*t95+t100*t24*t41+Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = -6.0d+0*rhoa*rhob*t102*t26*t27*t73*wght
     1         -2.0d+0*rhoa*rhob*t101*t24*t27*t71*wght-4.0d+0*rhoa*rhob*
     2         t27*t28*t41*t55*t56*wght+t120*t24*t41+t118*t26*t28+t116*t
     3         24*t26+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = -6.0d+0*rhoa*rhob*t26*t27*t40*t63*t73*
     1         wght-2.0d+0*rhoa*rhob*t24*t27*t42*t64*t71*wght+t128*t28*t
     2         41+t132*t24*t41+t130*t26*t28+t127*t24*t26+Cmat2(iq,D2_RA_
     3         GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = -6.0d+0*rhoa*rhob*t26*t27*t40*t69*t73*
     1         wght-2.0d+0*rhoa*rhob*t24*t27*t42*t70*t71*wght+t135*t28*t
     2         41+t138*t24*t41+t139*t26*t28+t134*t24*t26+Cmat2(iq,D2_RA_
     3         GBB)
            Cmat2(iq,D2_RB_GAA) = -6.0d+0*rhoa*rhob*t26*t27*t55*t63*t73*
     1         wght-2.0d+0*rhoa*rhob*t24*t27*t56*t64*t71*wght+t142*t28*t
     2         41+t146*t24*t41+t144*t26*t28+t141*t24*t26+Cmat2(iq,D2_RB_
     3         GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = -6.0d+0*rhoa*rhob*t26*t27*t55*t69*t73*
     1         wght-2.0d+0*rhoa*rhob*t24*t27*t56*t70*t71*wght+t154*t28*t
     2         41+t157*t24*t41+t158*t26*t28+t153*t24*t26+Cmat2(iq,D2_RB_
     3         GBB)
            Cmat2(iq,D2_GAA_GAA) = -6.0d+0*rhoa*rhob*t168*t26*t27*t73*wg
     1         ht-2.0d+0*rhoa*rhob*t167*t24*t27*t71*wght-4.0d+0*rhoa*rho
     2         b*t27*t28*t41*t63*t64*wght-1.0d+0*rhoa*rhob*t166*t24*t26*
     3         wght+t172*t24*t41+t170*t26*t28+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = -6.0d+0*rhoa*rhob*t26*t27*t63*t69*t73
     1         *wght-2.0d+0*rhoa*rhob*t24*t27*t64*t70*t71*wght-1.0d+0*rh
     2         oa*rhob*t173*t24*t26*wght+t176*t28*t41+t178*t24*t41+t175*
     3         t26*t28+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -6.0d+0*rhoa*rhob*t180*t26*t27*t73*wg
     1         ht-2.0d+0*rhoa*rhob*t179*t24*t27*t71*wght-4.0d+0*rhoa*rho
     2         b*t27*t28*t41*t69*t70*wght-1.0d+0*rhoa*rhob*t188*t24*t26*
     3         wght+t192*t24*t41+t190*t26*t28+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = t24*t41*(2.0d+0*rhoa*rhob*t43*t92*wg
     1         ht+2.0d+0*rhob*t27*t92*wght+1.0d+0*rhoa*rhob*t27*(-4.0d+0
     2         *t1*t10*t18*t19*t38*t39*t9*t90+1.3333333333333333d+0*t1*t
     3         18*t19*t21*t30*t9*t90+4.0d+0*t1*t10*t18*t19*t21*t37*t90+6
     4         .770400000000001d+0*t1*t10*t18*t19*t38*t76*t87*t9-2.0d+0*
     5         t1*t10*t18*t19*t39*t40*t87*t9-1.1284000000000001d+0*t1*t1
     6         8*t19*t30*t39*t87*t9+4.1792592592592587d-1*t1*t18*t19*t21
     7         *t8*t9-4.444444444444444d-1*t1*t18*t19*t21*t40*t77*t9+7.5
     8         22666666666665d-1*t1*t18*t19*t38*t39*t77*t9+4.0d+0*t1*t10
     9         *t18*t19*t40*t75*t76*t9+2.2568000000000002d+0*t1*t18*t19*
     :         t30*t75*t76*t9-1.3333333333333333d+0*t1*t18*t19*t30*t38*t
     ;         39*t40*t9-1.1284000000000001d+0*t1*t10*t18*t19*t201*t39*t
     <         9+2.0d+0*t1*t10*t18*t19*t202*t21*t9-6.770400000000001d+0*
     =         t1*t10*t18*t19*t195*t196*t9-3.3852000000000004d+0*t1*t10*
     >         t18*t19*t37*t39*t87+2.0d+0*t1*t10*t18*t19*t21*t40*t86-3.3
     ?         852000000000004d+0*t1*t10*t18*t19*t38*t39*t86+1.128400000
     @         0000001d+0*t1*t18*t19*t21*t30*t86-7.522666666666665d-1*t1
     1         *t18*t19*t21*t37*t77+6.770400000000001d+0*t1*t10*t18*t19*
     2         t37*t75*t76-4.0d+0*t1*t10*t18*t19*t37*t38*t39*t40+1.33333
     3         33333333333d+0*t1*t18*t19*t21*t30*t37*t40-2.2568000000000
     4         002d+0*t1*t18*t19*t30*t37*t38*t39+1.1284000000000001d+0*t
     5         1*t10*t18*t19*t200*t21)*wght+1.0d+0*rhoa*rhob*t42*t88*wgh
     6         t+2.0d+0*rhob*t42*t43*wght-t44*t92-2*t42*t89)+t24*t71*(-4
     7         .0d+0*rhoa*rhob*t27*t42*t92*wght-2.0d+0*rhoa*rhob*t43*t72
     8         *wght-2.0d+0*rhob*t27*t72*wght-2*t42*t93)+t28*t41*(-4.0d+
     9         0*rhoa*rhob*t27*t40*t92*wght-4.0d+0*rhoa*rhob*t27*t42*t90
     :         *wght-4.0d+0*rhoa*rhob*t40*t42*t43*wght-4.0d+0*rhob*t27*t
     ;         40*t42*wght-2.0d+0*t40*t93-t42*t91)+t26*t28*(4.0d+0*rhoa*
     <         rhob*t43*t90*wght+4.0d+0*rhob*t27*t90*wght+2.0d+0*rhoa*rh
     =         ob*t40*t88*wght+4.0d+0*rhob*t40*t43*wght+2.0d+0*rhoa*rhob
     >         *t202*t27*wght-2.0d+0*t44*t90-4.0d+0*t40*t89)+t26*t73*(-1
     ?         .2d+1*rhoa*rhob*t27*t40*t90*wght-6.0d+0*rhoa*rhob*t43*t74
     @         *wght-6.0d+0*rhob*t27*t74*wght-3.0d+0*t40*t91)+t24*t26*(-
     1         1.0d+0*rhoa*rhob*(9.128400000000001d+0*t1*t10*t18*t19*t38
     2         *t76*t87*t9-1.5214d+0*t1*t18*t19*t30*t39*t87*t9+5.6348148
     3         14814815d-1*t1*t18*t19*t21*t8*t9+1.0142666666666666d+0*t1
     4         *t18*t19*t38*t39*t77*t9+3.0428d+0*t1*t18*t19*t30*t75*t76*
     5         t9-1.5214d+0*t1*t10*t18*t19*t201*t39*t9-9.128400000000001
     6         d+0*t1*t10*t18*t19*t195*t196*t9-4.5642000000000005d+0*t1*
     7         t10*t18*t19*t37*t39*t87-4.5642000000000005d+0*t1*t10*t18*
     8         t19*t38*t39*t86+1.5214d+0*t1*t18*t19*t21*t30*t86-1.014266
     9         6666666666d+0*t1*t18*t19*t21*t37*t77+9.128400000000001d+0
     :         *t1*t10*t18*t19*t37*t75*t76-3.0428d+0*t1*t18*t19*t30*t37*
     ;         t38*t39+1.5214d+0*t1*t10*t18*t19*t200*t21)*wght-3.0d+0*rh
     <         ob*t88*wght)+1.8d+1*rhoa*rhob*t27*t41*t42*t73*t74*wght+1.
     =         2d+1*rhoa*rhob*t27*t28*t40*t71*t72*wght+6.0d+0*rhoa*rhob*
     >         t193*t24*t27*t42**3*wght+2.4d+1*rhoa*rhob*t194*t26*t27*t4
     ?         0**3*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t24*t41*(1.0d+0*rhoa*rhob*t43*t99*wg
     1         ht+1.0d+0*rhob*t27*t99*wght+1.0d+0*rhoa*rhob*t27*(-2.0d+0
     2         *t1*t10*t18*t19*t38*t39*t9*t97+6.666666666666666d-1*t1*t1
     3         8*t19*t21*t30*t9*t97+2.0d+0*t1*t10*t18*t19*t21*t37*t97-2.
     4         0d+0*t1*t10*t18*t19*t39*t54*t9*t90+2.0d+0*t1*t10*t19*t21*
     5         t53*t9*t90+6.666666666666666d-1*t1*t10*t18*t21*t46*t9*t90
     6         +2.2568000000000002d+0*t1*t10*t18*t19*t54*t76*t87*t9-1.12
     7         84000000000001d+0*t1*t10*t19*t39*t53*t87*t9-3.76133333333
     8         3333d-1*t1*t10*t18*t39*t46*t87*t9+2.507555555555555d-1*t1
     9         *t18*t19*t39*t54*t77*t9-2.507555555555555d-1*t1*t19*t21*t
     :         53*t77*t9-8.358518518518516d-2*t1*t18*t21*t46*t77*t9+2.25
     ;         68000000000002d+0*t1*t10*t19*t53*t75*t76*t9+7.52266666666
     <         6666d-1*t1*t10*t18*t46*t75*t76*t9+4.0d+0*t1*t10*t18*t19*t
     =         38*t40*t54*t76*t9+1.5045333333333333d+0*t1*t18*t19*t30*t3
     >         8*t54*t76*t9-6.770400000000001d+0*t1*t10*t18*t19*t196*t54
     ?         *t75*t9-6.666666666666666d-1*t1*t18*t19*t30*t39*t40*t54*t
     @         9-2.0d+0*t1*t10*t19*t38*t39*t40*t53*t9+6.666666666666666d
     1         -1*t1*t19*t21*t30*t40*t53*t9-7.522666666666666d-1*t1*t19*
     2         t30*t38*t39*t53*t9-6.666666666666666d-1*t1*t10*t18*t38*t3
     3         9*t40*t46*t9+2.222222222222222d-1*t1*t18*t21*t30*t40*t46*
     4         t9-2.507555555555555d-1*t1*t18*t30*t38*t39*t46*t9+2.0d+0*
     5         t1*t10*t18*t19*t203*t21*t9-1.1284000000000001d+0*t1*t10*t
     6         18*t19*t39*t54*t86+1.1284000000000001d+0*t1*t10*t19*t21*t
     7         53*t86+3.761333333333333d-1*t1*t10*t18*t21*t46*t86+4.5136
     8         000000000004d+0*t1*t10*t18*t19*t37*t38*t54*t76-2.0d+0*t1*
     9         t10*t18*t19*t37*t39*t40*t54-7.522666666666666d-1*t1*t18*t
     :         19*t30*t37*t39*t54+2.0d+0*t1*t10*t19*t21*t37*t40*t53-2.25
     ;         68000000000002d+0*t1*t10*t19*t37*t38*t39*t53+7.5226666666
     <         66666d-1*t1*t19*t21*t30*t37*t53+6.666666666666666d-1*t1*t
     =         10*t18*t21*t37*t40*t46-7.522666666666666d-1*t1*t10*t18*t3
     >         7*t38*t39*t46+2.507555555555555d-1*t1*t18*t21*t30*t37*t46
     ?         )*wght+1.0d+0*rhoa*rhob*t42*t94*wght+1.0d+0*rhoa*rhob*t57
     @         *t92*wght+1.0d+0*rhoa*t27*t92*wght+1.0d+0*rhob*t42*t57*wg
     1         ht+1.0d+0*rhoa*t42*t43*wght+1.0d+0*t27*t42*wght-t44*t99-t
     2         42*t95-t56*t89)+t24*t71*(-4.0d+0*rhoa*rhob*t27*t42*t99*wg
     3         ht-2.0d+0*rhoa*rhob*t57*t72*wght-2.0d+0*rhoa*t27*t72*wght
     4         -2*t56*t93)+t28*t41*(-4.0d+0*rhoa*rhob*t27*t40*t99*wght-4
     5         .0d+0*rhoa*rhob*t27*t42*t97*wght-4.0d+0*rhoa*rhob*t40*t42
     6         *t57*wght-4.0d+0*rhoa*t27*t40*t42*wght-2.0d+0*t55*t93-t56
     7         *t91)+t26*t28*(2.0d+0*rhoa*rhob*t43*t97*wght+2.0d+0*rhob*
     8         t27*t97*wght+2.0d+0*rhoa*rhob*t40*t94*wght+2.0d+0*rhoa*rh
     9         ob*t57*t90*wght+2.0d+0*rhoa*t27*t90*wght+2.0d+0*rhob*t40*
     :         t57*wght+2.0d+0*rhoa*t40*t43*wght+2.0d+0*t27*t40*wght+2.0
     ;         d+0*rhoa*rhob*t203*t27*wght-2.0d+0*t44*t97-2.0d+0*t40*t95
     <         -2.0d+0*t55*t89)+t26*t73*(-1.2d+1*rhoa*rhob*t27*t40*t97*w
     =         ght-6.0d+0*rhoa*rhob*t57*t74*wght-6.0d+0*rhoa*t27*t74*wgh
     >         t-3.0d+0*t55*t91)+t24*t26*(-2.0d+0*rhob*t94*wght-1.0d+0*r
     ?         hoa*rhob*(3.0428d+0*t1*t10*t18*t19*t54*t76*t87*t9-1.5214d
     @         +0*t1*t10*t19*t39*t53*t87*t9-5.071333333333333d-1*t1*t10*
     1         t18*t39*t46*t87*t9+3.380888888888889d-1*t1*t18*t19*t39*t5
     2         4*t77*t9-3.380888888888889d-1*t1*t19*t21*t53*t77*t9-1.126
     3         962962962963d-1*t1*t18*t21*t46*t77*t9+3.0428d+0*t1*t10*t1
     4         9*t53*t75*t76*t9+1.0142666666666666d+0*t1*t10*t18*t46*t75
     5         *t76*t9+2.0285333333333333d+0*t1*t18*t19*t30*t38*t54*t76*
     6         t9-9.128400000000001d+0*t1*t10*t18*t19*t196*t54*t75*t9-1.
     7         0142666666666666d+0*t1*t19*t30*t38*t39*t53*t9-3.380888888
     8         888889d-1*t1*t18*t30*t38*t39*t46*t9-1.5214d+0*t1*t10*t18*
     9         t19*t39*t54*t86+1.5214d+0*t1*t10*t19*t21*t53*t86+5.071333
     :         333333333d-1*t1*t10*t18*t21*t46*t86+6.0856d+0*t1*t10*t18*
     ;         t19*t37*t38*t54*t76-1.0142666666666666d+0*t1*t18*t19*t30*
     <         t37*t39*t54-3.0428d+0*t1*t10*t19*t37*t38*t39*t53+1.014266
     =         6666666666d+0*t1*t19*t21*t30*t37*t53-1.0142666666666666d+
     >         0*t1*t10*t18*t37*t38*t39*t46+3.380888888888889d-1*t1*t18*
     ?         t21*t30*t37*t46)*wght-1.0d+0*rhoa*t88*wght-2.0d+0*t43*wgh
     @         t)+t41*t73*(6.0d+0*rhoa*rhob*t27*t56*t74*wght+1.2d+1*rhoa
     1         *rhob*t27*t40*t42*t55*wght)+t28*t71*(4.0d+0*rhoa*rhob*t27
     2         *t55*t72*wght+8.0d+0*rhoa*rhob*t27*t40*t42*t56*wght)+2.4d
     3         +1*rhoa*rhob*t194*t26*t27*t55*t74*wght+6.0d+0*rhoa*rhob*t
     4         193*t24*t27*t56*t72*wght+Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = t24*t41*(2.0d+0*rhoa*rhob*t57*t99*wg
     1         ht+2.0d+0*rhoa*t27*t99*wght+1.0d+0*rhoa*rhob*t27*(-4.0d+0
     2         *t1*t10*t18*t19*t39*t54*t9*t97+4.0d+0*t1*t10*t19*t21*t53*
     3         t9*t97+1.3333333333333333d+0*t1*t10*t18*t21*t46*t9*t97+4.
     4         5136000000000004d+0*t1*t10*t19*t38*t53*t54*t76*t9+1.50453
     5         33333333333d+0*t1*t10*t18*t38*t46*t54*t76*t9+4.0d+0*t1*t1
     6         0*t104*t18*t19*t40*t76*t9+2.2568000000000002d+0*t1*t10*t1
     7         14*t18*t19*t38*t76*t9+7.522666666666666d-1*t1*t104*t18*t1
     8         9*t30*t76*t9-4.0d+0*t1*t10*t19*t39*t40*t53*t54*t9-7.52266
     9         6666666666d-1*t1*t19*t30*t39*t53*t54*t9-1.333333333333333
     :         3d+0*t1*t10*t18*t39*t40*t46*t54*t9-2.507555555555555d-1*t
     ;         1*t18*t30*t39*t46*t54*t9+1.3333333333333333d+0*t1*t10*t21
     <         *t40*t46*t53*t9-7.522666666666666d-1*t1*t10*t38*t39*t46*t
     =         53*t9+2.507555555555555d-1*t1*t21*t30*t46*t53*t9-2.0d+0*t
     >         1*t10*t114*t18*t19*t39*t40*t9+2.0d+0*t1*t10*t113*t19*t21*
     ?         t40*t9-4.444444444444444d-1*t1*t10*t103*t18*t21*t40*t9-1.
     @         1284000000000001d+0*t1*t10*t113*t19*t38*t39*t9+2.50755555
     1         5555555d-1*t1*t10*t103*t18*t38*t39*t9-3.761333333333333d-
     2         1*t1*t114*t18*t19*t30*t39*t9-6.770400000000001d+0*t1*t10*
     3         t104*t18*t19*t196*t38*t9+3.761333333333333d-1*t1*t113*t19
     4         *t21*t30*t9-8.358518518518516d-2*t1*t103*t18*t21*t30*t9+2
     5         .0d+0*t1*t10*t18*t19*t204*t21*t9+2.2568000000000002d+0*t1
     6         *t10*t104*t18*t19*t37*t76-2.2568000000000002d+0*t1*t10*t1
     7         9*t37*t39*t53*t54-7.522666666666666d-1*t1*t10*t18*t37*t39
     8         *t46*t54+7.522666666666666d-1*t1*t10*t21*t37*t46*t53-1.12
     9         84000000000001d+0*t1*t10*t114*t18*t19*t37*t39+1.128400000
     :         0000001d+0*t1*t10*t113*t19*t21*t37-2.507555555555555d-1*t
     ;         1*t10*t103*t18*t21*t37)*wght+2.0d+0*rhoa*t42*t57*wght+1.0
     <         d+0*rhoa*rhob*t115*t42*wght-2*t56*t95-t119*t44)+t24*t71*(
     =         -2.0d+0*rhoa*rhob*t27*t56*t99*wght-2.0d+0*rhoa*rhob*t42*t
     >         56*t57*wght-2.0d+0*rhoa*t27*t42*t56*wght-2.0d+0*rhoa*rhob
     ?         *t119*t27*t42*wght-2*t100*t56)+t28*t41*(-2.0d+0*rhoa*rhob
     @         *t27*t55*t99*wght-2.0d+0*rhoa*rhob*t27*t56*t97*wght-2.0d+
     1         0*rhoa*rhob*t40*t56*t57*wght-2.0d+0*rhoa*rhob*t42*t55*t57
     2         *wght-2.0d+0*rhoa*t27*t40*t56*wght-2.0d+0*rhoa*t27*t42*t5
     3         5*wght-2.0d+0*rhoa*rhob*t117*t27*t42*wght-2.0d+0*rhoa*rho
     4         b*t119*t27*t40*wght-t56*t98-2.0d+0*t100*t55)+t26*t28*(4.0
     5         d+0*rhoa*rhob*t57*t97*wght+4.0d+0*rhoa*t27*t97*wght+4.0d+
     6         0*rhoa*t40*t57*wght+2.0d+0*rhoa*rhob*t115*t40*wght+2.0d+0
     7         *rhoa*rhob*t204*t27*wght-4.0d+0*t55*t95-2.0d+0*t117*t44)+
     8         t26*t73*(-6.0d+0*rhoa*rhob*t27*t55*t97*wght-6.0d+0*rhoa*r
     9         hob*t40*t55*t57*wght-6.0d+0*rhoa*t27*t40*t55*wght-6.0d+0*
     :         rhoa*rhob*t117*t27*t40*wght-3.0d+0*t55*t98)+t24*t26*(-2.0
     ;         d+0*rhoa*t94*wght-1.0d+0*rhoa*rhob*(6.0856d+0*t1*t10*t19*
     <         t38*t53*t54*t76*t9+2.0285333333333333d+0*t1*t10*t18*t38*t
     =         46*t54*t76*t9+3.0428d+0*t1*t10*t114*t18*t19*t38*t76*t9+1.
     >         0142666666666666d+0*t1*t104*t18*t19*t30*t76*t9-1.01426666
     ?         66666666d+0*t1*t19*t30*t39*t53*t54*t9-3.380888888888889d-
     @         1*t1*t18*t30*t39*t46*t54*t9-1.0142666666666666d+0*t1*t10*
     1         t38*t39*t46*t53*t9+3.380888888888889d-1*t1*t21*t30*t46*t5
     2         3*t9-1.5214d+0*t1*t10*t113*t19*t38*t39*t9+3.3808888888888
     3         89d-1*t1*t10*t103*t18*t38*t39*t9-5.071333333333333d-1*t1*
     4         t114*t18*t19*t30*t39*t9-9.128400000000001d+0*t1*t10*t104*
     5         t18*t19*t196*t38*t9+5.071333333333333d-1*t1*t113*t19*t21*
     6         t30*t9-1.126962962962963d-1*t1*t103*t18*t21*t30*t9+3.0428
     7         d+0*t1*t10*t104*t18*t19*t37*t76-3.0428d+0*t1*t10*t19*t37*
     8         t39*t53*t54-1.0142666666666666d+0*t1*t10*t18*t37*t39*t46*
     9         t54+1.0142666666666666d+0*t1*t10*t21*t37*t46*t53-1.5214d+
     :         0*t1*t10*t114*t18*t19*t37*t39+1.5214d+0*t1*t10*t113*t19*t
     ;         21*t37-3.380888888888889d-1*t1*t10*t103*t18*t21*t37)*wght
     <         -2.0d+0*t57*wght-1.0d+0*rhob*t115*wght)+t28*t71*(4.0d+0*r
     =         hoa*rhob*t27*t42*t55*t56*wght-2*t56*t96)+t41*t73*(6.0d+0*
     >         rhoa*rhob*t27*t40*t55*t56*wght-3.0d+0*t55*t96)+6.0d+0*rho
     ?         a*rhob*t101*t193*t24*t27*t42*wght+2.4d+1*rhoa*rhob*t102*t
     @         194*t26*t27*t40*wght+Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0d+
     1         0*t1*t10*t104*t18*t19*t55*t76*t9+6.770400000000001d+0*t1*
     2         t10*t114*t18*t19*t54*t76*t9+6.770400000000001d+0*t1*t10*t
     3         104*t19*t53*t76*t9+2.2568000000000002d+0*t1*t10*t104*t18*
     4         t46*t76*t9-4.0d+0*t1*t10*t19*t39*t53*t54*t55*t9-1.3333333
     5         333333333d+0*t1*t10*t18*t39*t46*t54*t55*t9+1.333333333333
     6         3333d+0*t1*t10*t21*t46*t53*t55*t9-2.0d+0*t1*t10*t114*t18*
     7         t19*t39*t55*t9+2.0d+0*t1*t10*t113*t19*t21*t55*t9-4.444444
     8         444444444d-1*t1*t10*t103*t18*t21*t55*t9-2.256800000000000
     9         2d+0*t1*t10*t39*t46*t53*t54*t9-4.0d+0*t1*t10*t117*t18*t19
     :         *t39*t54*t9-3.3852000000000004d+0*t1*t10*t113*t19*t39*t54
     ;         *t9+7.522666666666665d-1*t1*t10*t103*t18*t39*t54*t9-3.385
     <         2000000000004d+0*t1*t10*t114*t19*t39*t53*t9+4.0d+0*t1*t10
     =         *t117*t19*t21*t53*t9-7.522666666666665d-1*t1*t10*t103*t21
     >         *t53*t9-1.1284000000000001d+0*t1*t10*t114*t18*t39*t46*t9+
     ?         1.3333333333333333d+0*t1*t10*t117*t18*t21*t46*t9+1.128400
     @         0000000001d+0*t1*t10*t113*t21*t46*t9-1.1284000000000001d+
     1         0*t1*t10*t18*t19*t210*t39*t9+2.0d+0*t1*t10*t18*t19*t21*t2
     2         11*t9+1.1284000000000001d+0*t1*t10*t19*t209*t21*t9+4.1792
     3         592592592587d-1*t1*t10*t17*t18*t21*t9-6.770400000000001d+
     4         0*t1*t10*t18*t19*t196*t205*t9)*wght+2.0d+0*rhoa*t56*t57*w
     5         ght+2.0d+0*rhoa*rhob*t119*t57*wght+1.0d+0*rhoa*rhob*t115*
     6         t56*wght+2.0d+0*rhoa*t119*t27*wght-t119*t58-2*t116*t56)+t
     7         24*t26*(-1.0d+0*rhoa*rhob*(9.128400000000001d+0*t1*t10*t1
     8         14*t18*t19*t54*t76*t9+9.128400000000001d+0*t1*t10*t104*t1
     9         9*t53*t76*t9+3.0428d+0*t1*t10*t104*t18*t46*t76*t9-3.0428d
     :         +0*t1*t10*t39*t46*t53*t54*t9-4.5642000000000005d+0*t1*t10
     ;         *t113*t19*t39*t54*t9+1.0142666666666666d+0*t1*t10*t103*t1
     <         8*t39*t54*t9-4.5642000000000005d+0*t1*t10*t114*t19*t39*t5
     =         3*t9-1.0142666666666666d+0*t1*t10*t103*t21*t53*t9-1.5214d
     >         +0*t1*t10*t114*t18*t39*t46*t9+1.5214d+0*t1*t10*t113*t21*t
     ?         46*t9-1.5214d+0*t1*t10*t18*t19*t210*t39*t9+1.5214d+0*t1*t
     @         10*t19*t209*t21*t9+5.634814814814815d-1*t1*t10*t17*t18*t2
     1         1*t9-9.128400000000001d+0*t1*t10*t18*t19*t196*t205*t9)*wg
     2         ht-3.0d+0*rhoa*t115*wght)+t28*t41*(-4.0d+0*rhoa*rhob*t55*
     3         t56*t57*wght-4.0d+0*rhoa*t27*t55*t56*wght-4.0d+0*rhoa*rho
     4         b*t117*t27*t56*wght-4.0d+0*rhoa*rhob*t119*t27*t55*wght-t1
     5         18*t56-2.0d+0*t120*t55)+t26*t28*(4.0d+0*rhoa*t55*t57*wght
     6         +4.0d+0*rhoa*rhob*t117*t57*wght+2.0d+0*rhoa*rhob*t115*t55
     7         *wght+2.0d+0*rhoa*rhob*t211*t27*wght+4.0d+0*rhoa*t117*t27
     8         *wght-2.0d+0*t117*t58-4.0d+0*t116*t55)+t26*t73*(-6.0d+0*r
     9         hoa*rhob*t102*t57*wght-1.2d+1*rhoa*rhob*t117*t27*t55*wght
     :         -6.0d+0*rhoa*t102*t27*wght-3.0d+0*t118*t55)+t24*t71*(-2.0
     ;         d+0*rhoa*rhob*t101*t57*wght-4.0d+0*rhoa*rhob*t119*t27*t56
     <         *wght-2.0d+0*rhoa*t101*t27*wght-2*t120*t56)+1.8d+1*rhoa*r
     =         hob*t102*t27*t41*t56*t73*wght+1.2d+1*rhoa*rhob*t101*t27*t
     >         28*t55*t71*wght+6.0d+0*rhoa*rhob*t193*t24*t27*t56**3*wght
     ?         +2.4d+1*rhoa*rhob*t194*t26*t27*t55**3*wght+Amat3(iq,D3_RB
     @         _RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t24*t41*(1.0d+0*rhoa*rhob*t62*t92*w
     1         ght+1.0d+0*rhoa*rhob*t27*(-2.0d+0*t1*t18*t19*t29*t39*t61*
     2         t9*t90+2.0d+0*t1*t10*t18*t19*t21*t61*t90+2.25680000000000
     3         02d+0*t1*t18*t19*t29*t61*t76*t87*t9+4.0d+0*t1*t18*t19*t29
     4         *t38*t40*t61*t76*t9+1.5045333333333333d+0*t1*t121*t18*t19
     5         *t38*t61*t76*t9+4.5136000000000004d+0*t1*t10*t125*t18*t19
     6         *t38*t76*t9-6.770400000000001d+0*t1*t18*t19*t196*t29*t61*
     7         t75*t9-6.666666666666666d-1*t1*t121*t18*t19*t39*t40*t61*t
     8         9+2.507555555555555d-1*t1*t18*t19*t39*t4*t61*t9-2.0d+0*t1
     9         *t10*t125*t18*t19*t39*t40*t9-2.0d+0*t1*t10*t129*t18*t19*t
     :         38*t39*t9-7.522666666666666d-1*t1*t125*t18*t19*t30*t39*t9
     ;         -1.1284000000000001d+0*t1*t10*t18*t19*t213*t39*t9+6.66666
     <         6666666666d-1*t1*t129*t18*t19*t21*t30*t9+2.0d+0*t1*t10*t1
     =         8*t19*t21*t214*t9-1.1284000000000001d+0*t1*t10*t18*t19*t3
     >         9*t61*t87-1.1284000000000001d+0*t1*t18*t19*t29*t39*t61*t8
     ?         6-2.507555555555555d-1*t1*t18*t19*t21*t61*t77+2.256800000
     @         0000002d+0*t1*t10*t18*t19*t61*t75*t76+4.5136000000000004d
     1         +0*t1*t18*t19*t29*t37*t38*t61*t76-2.0d+0*t1*t10*t18*t19*t
     2         38*t39*t40*t61-2.0d+0*t1*t18*t19*t29*t37*t39*t40*t61+6.66
     3         6666666666666d-1*t1*t18*t19*t21*t30*t40*t61-7.52266666666
     4         6666d-1*t1*t18*t19*t30*t38*t39*t61-7.522666666666666d-1*t
     5         1*t121*t18*t19*t37*t39*t61+2.0d+0*t1*t10*t124*t18*t19*t21
     6         *t40-2.2568000000000002d+0*t1*t10*t124*t18*t19*t38*t39-2.
     7         2568000000000002d+0*t1*t10*t125*t18*t19*t37*t39+2.0d+0*t1
     8         *t10*t129*t18*t19*t21*t37+7.522666666666666d-1*t1*t124*t1
     9         8*t19*t21*t30+1.1284000000000001d+0*t1*t10*t18*t19*t21*t2
     :         12)*wght+1.0d+0*rhob*t42*t62*wght+1.0d+0*rhoa*rhob*t131*t
     ;         43*wght+1.0d+0*rhoa*rhob*t126*t42*wght+1.0d+0*rhob*t131*t
     <         27*wght-t64*t89-t131*t44-t127*t42)+t26*t28*(2.0d+0*rhoa*r
     =         hob*t62*t90*wght+2.0d+0*rhob*t40*t62*wght+2.0d+0*rhoa*rho
     >         b*t129*t43*wght+2.0d+0*rhoa*rhob*t126*t40*wght+2.0d+0*rho
     ?         a*rhob*t214*t27*wght+2.0d+0*rhob*t129*t27*wght-2.0d+0*t63
     @         *t89-2.0d+0*t129*t44-2.0d+0*t127*t40)+t24*t26*(-1.0d+0*rh
     1         oa*rhob*(3.0428d+0*t1*t18*t19*t29*t61*t76*t87*t9+2.028533
     2         3333333333d+0*t1*t121*t18*t19*t38*t61*t76*t9+6.0856d+0*t1
     3         *t10*t125*t18*t19*t38*t76*t9-9.128400000000001d+0*t1*t18*
     4         t19*t196*t29*t61*t75*t9+3.380888888888889d-1*t1*t18*t19*t
     5         39*t4*t61*t9-1.0142666666666666d+0*t1*t125*t18*t19*t30*t3
     6         9*t9-1.5214d+0*t1*t10*t18*t19*t213*t39*t9-1.5214d+0*t1*t1
     7         0*t18*t19*t39*t61*t87-1.5214d+0*t1*t18*t19*t29*t39*t61*t8
     8         6-3.380888888888889d-1*t1*t18*t19*t21*t61*t77+3.0428d+0*t
     9         1*t10*t18*t19*t61*t75*t76+6.0856d+0*t1*t18*t19*t29*t37*t3
     :         8*t61*t76-1.0142666666666666d+0*t1*t18*t19*t30*t38*t39*t6
     ;         1-1.0142666666666666d+0*t1*t121*t18*t19*t37*t39*t61-3.042
     <         8d+0*t1*t10*t124*t18*t19*t38*t39-3.0428d+0*t1*t10*t125*t1
     =         8*t19*t37*t39+1.0142666666666666d+0*t1*t124*t18*t19*t21*t
     >         30+1.5214d+0*t1*t10*t18*t19*t21*t212)*wght-2.0d+0*rhob*t1
     ?         26*wght)+t41*t73*(6.0d+0*rhoa*rhob*t27*t64*t74*wght+1.2d+
     @         1*rhoa*rhob*t27*t40*t42*t63*wght)+t26*t73*(-6.0d+0*rhoa*r
     1         hob*t62*t74*wght-1.2d+1*rhoa*rhob*t129*t27*t40*wght-3.0d+
     2         0*t63*t91)+t28*t71*(4.0d+0*rhoa*rhob*t27*t63*t72*wght+8.0
     3         d+0*rhoa*rhob*t27*t40*t42*t64*wght)+t24*t71*(-2.0d+0*rhoa
     4         *rhob*t62*t72*wght-4.0d+0*rhoa*rhob*t131*t27*t42*wght-2*t
     5         64*t93)+t28*t41*(-4.0d+0*rhoa*rhob*t40*t42*t62*wght-4.0d+
     6         0*rhoa*rhob*t129*t27*t42*wght-4.0d+0*rhoa*rhob*t131*t27*t
     7         40*wght-2.0d+0*t63*t93-t64*t91)+2.4d+1*rhoa*rhob*t194*t26
     8         *t27*t63*t74*wght+6.0d+0*rhoa*rhob*t193*t24*t27*t64*t72*w
     9         ght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t68*t92*w
     1         ght+1.0d+0*rhoa*rhob*t27*(-2.0d+0*t1*t10*t18*t39*t45*t67*
     2         t9*t90+2.0d+0*t1*t10*t19*t21*t67*t9*t90+2.256800000000000
     3         2d+0*t1*t10*t18*t45*t67*t76*t87*t9-1.1284000000000001d+0*
     4         t1*t10*t19*t39*t67*t87*t9+2.507555555555555d-1*t1*t18*t39
     5         *t45*t67*t77*t9-2.507555555555555d-1*t1*t19*t21*t67*t77*t
     6         9+2.2568000000000002d+0*t1*t10*t19*t67*t75*t76*t9+4.0d+0*
     7         t1*t10*t18*t38*t40*t45*t67*t76*t9+1.5045333333333333d+0*t
     8         1*t18*t30*t38*t45*t67*t76*t9-6.770400000000001d+0*t1*t10*
     9         t18*t196*t45*t67*t75*t9-6.666666666666666d-1*t1*t18*t30*t
     :         39*t40*t45*t67*t9-2.0d+0*t1*t10*t19*t38*t39*t40*t67*t9+6.
     ;         666666666666666d-1*t1*t19*t21*t30*t40*t67*t9-7.5226666666
     <         66666d-1*t1*t19*t30*t38*t39*t67*t9-2.0d+0*t1*t10*t136*t18
     =         *t19*t38*t39*t9+6.666666666666666d-1*t1*t136*t18*t19*t21*
     >         t30*t9+2.0d+0*t1*t10*t18*t19*t21*t215*t9-1.12840000000000
     ?         01d+0*t1*t10*t18*t39*t45*t67*t86+1.1284000000000001d+0*t1
     @         *t10*t19*t21*t67*t86+4.5136000000000004d+0*t1*t10*t18*t37
     1         *t38*t45*t67*t76-2.0d+0*t1*t10*t18*t37*t39*t40*t45*t67-7.
     2         522666666666666d-1*t1*t18*t30*t37*t39*t45*t67+2.0d+0*t1*t
     3         10*t19*t21*t37*t40*t67-2.2568000000000002d+0*t1*t10*t19*t
     4         37*t38*t39*t67+7.522666666666666d-1*t1*t19*t21*t30*t37*t6
     5         7+2.0d+0*t1*t10*t136*t18*t19*t21*t37)*wght+1.0d+0*rhob*t4
     6         2*t68*wght+1.0d+0*rhoa*rhob*t137*t43*wght+1.0d+0*rhoa*rho
     7         b*t133*t42*wght+1.0d+0*rhob*t137*t27*wght-t70*t89-t137*t4
     8         4-t134*t42)+t26*t28*(2.0d+0*rhoa*rhob*t68*t90*wght+2.0d+0
     9         *rhob*t40*t68*wght+2.0d+0*rhoa*rhob*t136*t43*wght+2.0d+0*
     :         rhoa*rhob*t133*t40*wght+2.0d+0*rhoa*rhob*t215*t27*wght+2.
     ;         0d+0*rhob*t136*t27*wght-2.0d+0*t69*t89-2.0d+0*t136*t44-2.
     <         0d+0*t134*t40)+t24*t26*(-1.0d+0*rhoa*rhob*(3.0428d+0*t1*t
     =         10*t18*t45*t67*t76*t87*t9-1.5214d+0*t1*t10*t19*t39*t67*t8
     >         7*t9+3.380888888888889d-1*t1*t18*t39*t45*t67*t77*t9-3.380
     ?         888888888889d-1*t1*t19*t21*t67*t77*t9+3.0428d+0*t1*t10*t1
     @         9*t67*t75*t76*t9+2.0285333333333333d+0*t1*t18*t30*t38*t45
     1         *t67*t76*t9-9.128400000000001d+0*t1*t10*t18*t196*t45*t67*
     2         t75*t9-1.0142666666666666d+0*t1*t19*t30*t38*t39*t67*t9-1.
     3         5214d+0*t1*t10*t18*t39*t45*t67*t86+1.5214d+0*t1*t10*t19*t
     4         21*t67*t86+6.0856d+0*t1*t10*t18*t37*t38*t45*t67*t76-1.014
     5         2666666666666d+0*t1*t18*t30*t37*t39*t45*t67-3.0428d+0*t1*
     6         t10*t19*t37*t38*t39*t67+1.0142666666666666d+0*t1*t19*t21*
     7         t30*t37*t67)*wght-2.0d+0*rhob*t133*wght)+t41*t73*(6.0d+0*
     8         rhoa*rhob*t27*t70*t74*wght+1.2d+1*rhoa*rhob*t27*t40*t42*t
     9         69*wght)+t26*t73*(-6.0d+0*rhoa*rhob*t68*t74*wght-1.2d+1*r
     :         hoa*rhob*t136*t27*t40*wght-3.0d+0*t69*t91)+t28*t71*(4.0d+
     ;         0*rhoa*rhob*t27*t69*t72*wght+8.0d+0*rhoa*rhob*t27*t40*t42
     <         *t70*wght)+t24*t71*(-2.0d+0*rhoa*rhob*t68*t72*wght-4.0d+0
     =         *rhoa*rhob*t137*t27*t42*wght-2*t70*t93)+t28*t41*(-4.0d+0*
     >         rhoa*rhob*t40*t42*t68*wght-4.0d+0*rhoa*rhob*t136*t27*t42*
     ?         wght-4.0d+0*rhoa*rhob*t137*t27*t40*wght-2.0d+0*t69*t93-t7
     @         0*t91)+2.4d+1*rhoa*rhob*t194*t26*t27*t69*t74*wght+6.0d+0*
     1         rhoa*rhob*t193*t24*t27*t70*t72*wght+Cmat3(iq,D3_RA_RA_GBB
     2         )
            Cmat3(iq,D3_RA_RB_GAA) = t24*t41*(1.0d+0*rhoa*rhob*t62*t99*w
     1         ght+1.0d+0*rhoa*rhob*t27*(-2.0d+0*t1*t18*t19*t29*t39*t61*
     2         t9*t97+2.0d+0*t1*t10*t18*t19*t21*t61*t97+4.0d+0*t1*t18*t1
     3         9*t29*t40*t54*t61*t76*t9+7.522666666666666d-1*t1*t121*t18
     4         *t19*t54*t61*t76*t9+2.2568000000000002d+0*t1*t19*t29*t38*
     5         t53*t61*t76*t9+7.522666666666666d-1*t1*t18*t29*t38*t46*t6
     6         1*t76*t9+2.2568000000000002d+0*t1*t10*t125*t18*t19*t54*t7
     7         6*t9-6.770400000000001d+0*t1*t18*t19*t196*t29*t38*t54*t61
     8         *t9-2.0d+0*t1*t19*t29*t39*t40*t53*t61*t9-3.76133333333333
     9         3d-1*t1*t121*t19*t39*t53*t61*t9-6.666666666666666d-1*t1*t
     :         18*t29*t39*t40*t46*t61*t9-1.2537777777777775d-1*t1*t121*t
     ;         18*t39*t46*t61*t9-2.0d+0*t1*t10*t129*t18*t19*t39*t54*t9-1
     <         .1284000000000001d+0*t1*t10*t125*t19*t39*t53*t9+2.0d+0*t1
     =         *t10*t129*t19*t21*t53*t9-3.761333333333333d-1*t1*t10*t125
     >         *t18*t39*t46*t9+6.666666666666666d-1*t1*t10*t129*t18*t21*
     ?         t46*t9+2.0d+0*t1*t10*t18*t19*t21*t216*t9+2.25680000000000
     @         02d+0*t1*t10*t18*t19*t38*t54*t61*t76+2.2568000000000002d+
     1         0*t1*t18*t19*t29*t37*t54*t61*t76-2.0d+0*t1*t10*t18*t19*t3
     2         9*t40*t54*t61-3.761333333333333d-1*t1*t18*t19*t30*t39*t54
     3         *t61+2.0d+0*t1*t10*t19*t21*t40*t53*t61-1.1284000000000001
     4         d+0*t1*t10*t19*t38*t39*t53*t61-1.1284000000000001d+0*t1*t
     5         19*t29*t37*t39*t53*t61+3.761333333333333d-1*t1*t19*t21*t3
     6         0*t53*t61+6.666666666666666d-1*t1*t10*t18*t21*t40*t46*t61
     7         -3.761333333333333d-1*t1*t10*t18*t38*t39*t46*t61-3.761333
     8         333333333d-1*t1*t18*t29*t37*t39*t46*t61+1.253777777777777
     9         5d-1*t1*t18*t21*t30*t46*t61-1.1284000000000001d+0*t1*t10*
     :         t124*t18*t19*t39*t54+1.1284000000000001d+0*t1*t10*t124*t1
     ;         9*t21*t53+3.761333333333333d-1*t1*t10*t124*t18*t21*t46)*w
     <         ght+1.0d+0*rhoa*t42*t62*wght+1.0d+0*rhoa*rhob*t131*t57*wg
     =         ht+1.0d+0*rhoa*rhob*t140*t42*wght+1.0d+0*rhoa*t131*t27*wg
     >         ht-t64*t95-t127*t56-t145*t44)+t26*t28*(2.0d+0*rhoa*rhob*t
     ?         62*t97*wght+2.0d+0*rhoa*t40*t62*wght+2.0d+0*rhoa*rhob*t12
     @         9*t57*wght+2.0d+0*rhoa*rhob*t140*t40*wght+2.0d+0*rhoa*rho
     1         b*t216*t27*wght+2.0d+0*rhoa*t129*t27*wght-2.0d+0*t63*t95-
     2         2.0d+0*t127*t55-2.0d+0*t143*t44)+t24*t26*(-1.0d+0*rhoa*rh
     3         ob*(1.0142666666666666d+0*t1*t121*t18*t19*t54*t61*t76*t9+
     4         3.0428d+0*t1*t19*t29*t38*t53*t61*t76*t9+1.014266666666666
     5         6d+0*t1*t18*t29*t38*t46*t61*t76*t9+3.0428d+0*t1*t10*t125*
     6         t18*t19*t54*t76*t9-9.128400000000001d+0*t1*t18*t19*t196*t
     7         29*t38*t54*t61*t9-5.071333333333333d-1*t1*t121*t19*t39*t5
     8         3*t61*t9-1.6904444444444444d-1*t1*t121*t18*t39*t46*t61*t9
     9         -1.5214d+0*t1*t10*t125*t19*t39*t53*t9-5.071333333333333d-
     :         1*t1*t10*t125*t18*t39*t46*t9+3.0428d+0*t1*t10*t18*t19*t38
     ;         *t54*t61*t76+3.0428d+0*t1*t18*t19*t29*t37*t54*t61*t76-5.0
     <         71333333333333d-1*t1*t18*t19*t30*t39*t54*t61-1.5214d+0*t1
     =         *t10*t19*t38*t39*t53*t61-1.5214d+0*t1*t19*t29*t37*t39*t53
     >         *t61+5.071333333333333d-1*t1*t19*t21*t30*t53*t61-5.071333
     ?         333333333d-1*t1*t10*t18*t38*t39*t46*t61-5.071333333333333
     @         d-1*t1*t18*t29*t37*t39*t46*t61+1.6904444444444444d-1*t1*t
     1         18*t21*t30*t46*t61-1.5214d+0*t1*t10*t124*t18*t19*t39*t54+
     2         1.5214d+0*t1*t10*t124*t19*t21*t53+5.071333333333333d-1*t1
     3         *t10*t124*t18*t21*t46)*wght-1.0d+0*t62*wght-1.0d+0*rhob*t
     4         140*wght-1.0d+0*rhoa*t126*wght)+t41*t73*(6.0d+0*rhoa*rhob
     5         *t27*t40*t55*t64*wght-3.0d+0*t63*t96)+t28*t71*(4.0d+0*rho
     6         a*rhob*t27*t42*t56*t63*wght-2*t64*t96)+t24*t71*(-2.0d+0*r
     7         hoa*rhob*t42*t56*t62*wght-2.0d+0*rhoa*rhob*t131*t27*t56*w
     8         ght-2.0d+0*rhoa*rhob*t145*t27*t42*wght-2*t100*t64)+t28*t4
     9         1*(-2.0d+0*rhoa*rhob*t40*t56*t62*wght-2.0d+0*rhoa*rhob*t4
     :         2*t55*t62*wght-2.0d+0*rhoa*rhob*t129*t27*t56*wght-2.0d+0*
     ;         rhoa*rhob*t131*t27*t55*wght-2.0d+0*rhoa*rhob*t143*t27*t42
     <         *wght-2.0d+0*rhoa*rhob*t145*t27*t40*wght-t64*t98-2.0d+0*t
     =         100*t63)+t26*t73*(-6.0d+0*rhoa*rhob*t40*t55*t62*wght-6.0d
     >         +0*rhoa*rhob*t129*t27*t55*wght-6.0d+0*rhoa*rhob*t143*t27*
     ?         t40*wght-3.0d+0*t63*t98)+6.0d+0*rhoa*rhob*t193*t24*t27*t4
     @         2*t56*t64*wght+2.4d+1*rhoa*rhob*t194*t26*t27*t40*t55*t63*
     1         wght+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t68*t99*w
     1         ght+1.0d+0*rhoa*rhob*t27*(-2.0d+0*t1*t10*t18*t39*t45*t67*
     2         t9*t97+2.0d+0*t1*t10*t19*t21*t67*t9*t97+4.0d+0*t1*t10*t18
     3         *t40*t45*t54*t67*t76*t9+7.522666666666666d-1*t1*t18*t30*t
     4         45*t54*t67*t76*t9+2.2568000000000002d+0*t1*t10*t19*t38*t5
     5         4*t67*t76*t9+2.2568000000000002d+0*t1*t10*t38*t45*t53*t67
     6         *t76*t9+7.522666666666666d-1*t1*t10*t147*t18*t38*t67*t76*
     7         t9+2.2568000000000002d+0*t1*t10*t151*t18*t19*t38*t76*t9-6
     8         .770400000000001d+0*t1*t10*t18*t196*t38*t45*t54*t67*t9-2.
     9         0d+0*t1*t10*t19*t39*t40*t54*t67*t9-3.761333333333333d-1*t
     :         1*t19*t30*t39*t54*t67*t9-2.0d+0*t1*t10*t39*t40*t45*t53*t6
     ;         7*t9-3.761333333333333d-1*t1*t30*t39*t45*t53*t67*t9+6.666
     <         666666666666d-1*t1*t10*t21*t40*t46*t67*t9-3.7613333333333
     =         33d-1*t1*t10*t38*t39*t46*t67*t9+1.2537777777777775d-1*t1*
     >         t21*t30*t46*t67*t9-6.666666666666666d-1*t1*t10*t147*t18*t
     ?         39*t40*t67*t9-1.2537777777777775d-1*t1*t147*t18*t30*t39*t
     @         67*t9-2.0d+0*t1*t10*t136*t18*t19*t39*t54*t9+2.0d+0*t1*t10
     1         *t136*t19*t21*t53*t9+6.666666666666666d-1*t1*t10*t136*t18
     2         *t21*t46*t9-2.0d+0*t1*t10*t151*t18*t19*t39*t40*t9+2.0d+0*
     3         t1*t10*t150*t19*t21*t40*t9-1.1284000000000001d+0*t1*t10*t
     4         150*t19*t38*t39*t9-3.761333333333333d-1*t1*t151*t18*t19*t
     5         30*t39*t9+3.761333333333333d-1*t1*t150*t19*t21*t30*t9+2.0
     6         d+0*t1*t10*t18*t19*t21*t217*t9+2.2568000000000002d+0*t1*t
     7         10*t18*t37*t45*t54*t67*t76-1.1284000000000001d+0*t1*t10*t
     8         19*t37*t39*t54*t67-1.1284000000000001d+0*t1*t10*t37*t39*t
     9         45*t53*t67+3.761333333333333d-1*t1*t10*t21*t37*t46*t67-3.
     :         761333333333333d-1*t1*t10*t147*t18*t37*t39*t67-1.12840000
     ;         00000001d+0*t1*t10*t151*t18*t19*t37*t39+1.128400000000000
     <         1d+0*t1*t10*t150*t19*t21*t37)*wght+1.0d+0*rhoa*t42*t68*wg
     =         ht+1.0d+0*rhoa*rhob*t137*t57*wght+1.0d+0*rhoa*rhob*t152*t
     >         42*wght+1.0d+0*rhoa*t137*t27*wght-t70*t95-t134*t56-t156*t
     ?         44)+t26*t28*(2.0d+0*rhoa*rhob*t68*t97*wght+2.0d+0*rhoa*t4
     @         0*t68*wght+2.0d+0*rhoa*rhob*t136*t57*wght+2.0d+0*rhoa*rho
     1         b*t152*t40*wght+2.0d+0*rhoa*rhob*t217*t27*wght+2.0d+0*rho
     2         a*t136*t27*wght-2.0d+0*t69*t95-2.0d+0*t134*t55-2.0d+0*t15
     3         5*t44)+t24*t26*(-1.0d+0*rhoa*rhob*(1.0142666666666666d+0*
     4         t1*t18*t30*t45*t54*t67*t76*t9+3.0428d+0*t1*t10*t19*t38*t5
     5         4*t67*t76*t9+3.0428d+0*t1*t10*t38*t45*t53*t67*t76*t9+1.01
     6         42666666666666d+0*t1*t10*t147*t18*t38*t67*t76*t9+3.0428d+
     7         0*t1*t10*t151*t18*t19*t38*t76*t9-9.128400000000001d+0*t1*
     8         t10*t18*t196*t38*t45*t54*t67*t9-5.071333333333333d-1*t1*t
     9         19*t30*t39*t54*t67*t9-5.071333333333333d-1*t1*t30*t39*t45
     :         *t53*t67*t9-5.071333333333333d-1*t1*t10*t38*t39*t46*t67*t
     ;         9+1.6904444444444444d-1*t1*t21*t30*t46*t67*t9-1.690444444
     <         4444444d-1*t1*t147*t18*t30*t39*t67*t9-1.5214d+0*t1*t10*t1
     =         50*t19*t38*t39*t9-5.071333333333333d-1*t1*t151*t18*t19*t3
     >         0*t39*t9+5.071333333333333d-1*t1*t150*t19*t21*t30*t9+3.04
     ?         28d+0*t1*t10*t18*t37*t45*t54*t67*t76-1.5214d+0*t1*t10*t19
     @         *t37*t39*t54*t67-1.5214d+0*t1*t10*t37*t39*t45*t53*t67+5.0
     1         71333333333333d-1*t1*t10*t21*t37*t46*t67-5.07133333333333
     2         3d-1*t1*t10*t147*t18*t37*t39*t67-1.5214d+0*t1*t10*t151*t1
     3         8*t19*t37*t39+1.5214d+0*t1*t10*t150*t19*t21*t37)*wght-1.0
     4         d+0*t68*wght-1.0d+0*rhob*t152*wght-1.0d+0*rhoa*t133*wght)
     5         +t41*t73*(6.0d+0*rhoa*rhob*t27*t40*t55*t70*wght-3.0d+0*t6
     6         9*t96)+t28*t71*(4.0d+0*rhoa*rhob*t27*t42*t56*t69*wght-2*t
     7         70*t96)+t24*t71*(-2.0d+0*rhoa*rhob*t42*t56*t68*wght-2.0d+
     8         0*rhoa*rhob*t137*t27*t56*wght-2.0d+0*rhoa*rhob*t156*t27*t
     9         42*wght-2*t100*t70)+t28*t41*(-2.0d+0*rhoa*rhob*t40*t56*t6
     :         8*wght-2.0d+0*rhoa*rhob*t42*t55*t68*wght-2.0d+0*rhoa*rhob
     ;         *t136*t27*t56*wght-2.0d+0*rhoa*rhob*t137*t27*t55*wght-2.0
     <         d+0*rhoa*rhob*t155*t27*t42*wght-2.0d+0*rhoa*rhob*t156*t27
     =         *t40*wght-t70*t98-2.0d+0*t100*t69)+t26*t73*(-6.0d+0*rhoa*
     >         rhob*t40*t55*t68*wght-6.0d+0*rhoa*rhob*t136*t27*t55*wght-
     ?         6.0d+0*rhoa*rhob*t155*t27*t40*wght-3.0d+0*t69*t98)+6.0d+0
     @         *rhoa*rhob*t193*t24*t27*t42*t56*t70*wght+2.4d+1*rhoa*rhob
     1         *t194*t26*t27*t40*t55*t69*wght+Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0d
     1         +0*t1*t18*t19*t29*t54*t55*t61*t76*t9+4.5136000000000004d+
     2         0*t1*t19*t29*t53*t54*t61*t76*t9+1.5045333333333333d+0*t1*
     3         t18*t29*t46*t54*t61*t76*t9+2.2568000000000002d+0*t1*t114*
     4         t18*t19*t29*t61*t76*t9-2.0d+0*t1*t19*t29*t39*t53*t55*t61*
     5         t9-6.666666666666666d-1*t1*t18*t29*t39*t46*t55*t61*t9-7.5
     6         22666666666666d-1*t1*t29*t39*t46*t53*t61*t9-2.0d+0*t1*t11
     7         7*t18*t19*t29*t39*t61*t9-1.1284000000000001d+0*t1*t113*t1
     8         9*t29*t39*t61*t9+2.507555555555555d-1*t1*t103*t18*t29*t39
     9         *t61*t9-6.770400000000001d+0*t1*t104*t18*t19*t196*t29*t61
     :         *t9-2.0d+0*t1*t10*t143*t18*t19*t39*t54*t9+2.0d+0*t1*t10*t
     ;         143*t19*t21*t53*t9+6.666666666666666d-1*t1*t10*t143*t18*t
     <         21*t46*t9+2.0d+0*t1*t10*t18*t19*t21*t218*t9+2.25680000000
     =         00002d+0*t1*t10*t104*t18*t19*t61*t76-2.0d+0*t1*t10*t18*t1
     >         9*t39*t54*t55*t61+2.0d+0*t1*t10*t19*t21*t53*t55*t61+6.666
     ?         666666666666d-1*t1*t10*t18*t21*t46*t55*t61-2.256800000000
     @         0002d+0*t1*t10*t19*t39*t53*t54*t61-7.522666666666666d-1*t
     1         1*t10*t18*t39*t46*t54*t61+7.522666666666666d-1*t1*t10*t21
     2         *t46*t53*t61-1.1284000000000001d+0*t1*t10*t114*t18*t19*t3
     3         9*t61+2.0d+0*t1*t10*t117*t18*t19*t21*t61+1.12840000000000
     4         01d+0*t1*t10*t113*t19*t21*t61-2.507555555555555d-1*t1*t10
     5         *t103*t18*t21*t61)*wght+1.0d+0*rhoa*t56*t62*wght+1.0d+0*r
     6         hoa*rhob*t119*t62*wght+1.0d+0*rhoa*rhob*t145*t57*wght+1.0
     7         d+0*rhoa*rhob*t140*t56*wght+1.0d+0*rhoa*t145*t27*wght-t11
     8         6*t64-t145*t58-t141*t56)+t24*t26*(-1.0d+0*rhoa*rhob*(6.08
     9         56d+0*t1*t19*t29*t53*t54*t61*t76*t9+2.0285333333333333d+0
     :         *t1*t18*t29*t46*t54*t61*t76*t9+3.0428d+0*t1*t114*t18*t19*
     ;         t29*t61*t76*t9-1.0142666666666666d+0*t1*t29*t39*t46*t53*t
     <         61*t9-1.5214d+0*t1*t113*t19*t29*t39*t61*t9+3.380888888888
     =         889d-1*t1*t103*t18*t29*t39*t61*t9-9.128400000000001d+0*t1
     >         *t104*t18*t19*t196*t29*t61*t9+3.0428d+0*t1*t10*t104*t18*t
     ?         19*t61*t76-3.0428d+0*t1*t10*t19*t39*t53*t54*t61-1.0142666
     @         666666666d+0*t1*t10*t18*t39*t46*t54*t61+1.014266666666666
     1         6d+0*t1*t10*t21*t46*t53*t61-1.5214d+0*t1*t10*t114*t18*t19
     2         *t39*t61+1.5214d+0*t1*t10*t113*t19*t21*t61-3.380888888888
     3         889d-1*t1*t10*t103*t18*t21*t61)*wght-2.0d+0*rhoa*t140*wgh
     4         t)+t28*t71*(8.0d+0*rhoa*rhob*t27*t55*t56*t64*wght+4.0d+0*
     5         rhoa*rhob*t101*t27*t63*wght)+t41*t73*(6.0d+0*rhoa*rhob*t1
     6         02*t27*t64*wght+1.2d+1*rhoa*rhob*t27*t55*t56*t63*wght)+t2
     7         8*t41*(-4.0d+0*rhoa*rhob*t55*t56*t62*wght-4.0d+0*rhoa*rho
     8         b*t143*t27*t56*wght-4.0d+0*rhoa*rhob*t145*t27*t55*wght-t1
     9         18*t64-2.0d+0*t120*t63)+t26*t28*(2.0d+0*rhoa*t55*t62*wght
     :         +2.0d+0*rhoa*rhob*t117*t62*wght+2.0d+0*rhoa*rhob*t143*t57
     ;         *wght+2.0d+0*rhoa*rhob*t140*t55*wght+2.0d+0*rhoa*rhob*t21
     <         8*t27*wght+2.0d+0*rhoa*t143*t27*wght-2.0d+0*t116*t63-2.0d
     =         +0*t143*t58-2.0d+0*t141*t55)+t26*t73*(-6.0d+0*rhoa*rhob*t
     >         102*t62*wght-1.2d+1*rhoa*rhob*t143*t27*t55*wght-3.0d+0*t1
     ?         18*t63)+t24*t71*(-2.0d+0*rhoa*rhob*t101*t62*wght-4.0d+0*r
     @         hoa*rhob*t145*t27*t56*wght-2*t120*t64)+6.0d+0*rhoa*rhob*t
     1         101*t193*t24*t27*t64*wght+2.4d+1*rhoa*rhob*t102*t194*t26*
     2         t27*t63*wght+Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0d
     1         +0*t1*t10*t18*t45*t54*t55*t67*t76*t9+4.5136000000000004d+
     2         0*t1*t10*t45*t53*t54*t67*t76*t9+1.5045333333333333d+0*t1*
     3         t10*t147*t18*t54*t67*t76*t9+2.2568000000000002d+0*t1*t10*
     4         t114*t18*t45*t67*t76*t9+2.2568000000000002d+0*t1*t10*t104
     5         *t19*t67*t76*t9+4.5136000000000004d+0*t1*t10*t151*t18*t19
     6         *t54*t76*t9-2.0d+0*t1*t10*t19*t39*t54*t55*t67*t9-2.0d+0*t
     7         1*t10*t39*t45*t53*t55*t67*t9+6.666666666666666d-1*t1*t10*
     8         t21*t46*t55*t67*t9-6.666666666666666d-1*t1*t10*t147*t18*t
     9         39*t55*t67*t9-7.522666666666666d-1*t1*t10*t39*t46*t54*t67
     :         *t9-7.522666666666666d-1*t1*t10*t147*t39*t53*t67*t9-2.0d+
     ;         0*t1*t10*t117*t18*t39*t45*t67*t9-1.1284000000000001d+0*t1
     <         *t10*t113*t39*t45*t67*t9-6.770400000000001d+0*t1*t10*t104
     =         *t18*t196*t45*t67*t9-1.1284000000000001d+0*t1*t10*t114*t1
     >         9*t39*t67*t9+2.507555555555555d-1*t1*t10*t13*t18*t39*t67*
     ?         t9+2.0d+0*t1*t10*t117*t19*t21*t67*t9-2.507555555555555d-1
     @         *t1*t10*t103*t21*t67*t9-2.0d+0*t1*t10*t151*t18*t19*t39*t5
     1         5*t9+2.0d+0*t1*t10*t150*t19*t21*t55*t9-2.0d+0*t1*t10*t155
     2         *t18*t19*t39*t54*t9-2.2568000000000002d+0*t1*t10*t150*t19
     3         *t39*t54*t9-2.2568000000000002d+0*t1*t10*t151*t19*t39*t53
     4         *t9+2.0d+0*t1*t10*t155*t19*t21*t53*t9-7.522666666666666d-
     5         1*t1*t10*t151*t18*t39*t46*t9+6.666666666666666d-1*t1*t10*
     6         t155*t18*t21*t46*t9+7.522666666666666d-1*t1*t10*t150*t21*
     7         t46*t9-1.1284000000000001d+0*t1*t10*t18*t19*t220*t39*t9+2
     8         .0d+0*t1*t10*t18*t19*t21*t221*t9+1.1284000000000001d+0*t1
     9         *t10*t19*t21*t219*t9)*wght+1.0d+0*rhoa*t56*t68*wght+1.0d+
     :         0*rhoa*rhob*t119*t68*wght+1.0d+0*rhoa*rhob*t156*t57*wght+
     ;         1.0d+0*rhoa*rhob*t152*t56*wght+1.0d+0*rhoa*t156*t27*wght-
     <         t116*t70-t156*t58-t153*t56)+t24*t26*(-1.0d+0*rhoa*rhob*(6
     =         .0856d+0*t1*t10*t45*t53*t54*t67*t76*t9+2.0285333333333333
     >         d+0*t1*t10*t147*t18*t54*t67*t76*t9+3.0428d+0*t1*t10*t114*
     ?         t18*t45*t67*t76*t9+3.0428d+0*t1*t10*t104*t19*t67*t76*t9+6
     @         .0856d+0*t1*t10*t151*t18*t19*t54*t76*t9-1.014266666666666
     1         6d+0*t1*t10*t39*t46*t54*t67*t9-1.0142666666666666d+0*t1*t
     2         10*t147*t39*t53*t67*t9-1.5214d+0*t1*t10*t113*t39*t45*t67*
     3         t9-9.128400000000001d+0*t1*t10*t104*t18*t196*t45*t67*t9-1
     4         .5214d+0*t1*t10*t114*t19*t39*t67*t9+3.380888888888889d-1*
     5         t1*t10*t13*t18*t39*t67*t9-3.380888888888889d-1*t1*t10*t10
     6         3*t21*t67*t9-3.0428d+0*t1*t10*t150*t19*t39*t54*t9-3.0428d
     7         +0*t1*t10*t151*t19*t39*t53*t9-1.0142666666666666d+0*t1*t1
     8         0*t151*t18*t39*t46*t9+1.0142666666666666d+0*t1*t10*t150*t
     9         21*t46*t9-1.5214d+0*t1*t10*t18*t19*t220*t39*t9+1.5214d+0*
     :         t1*t10*t19*t21*t219*t9)*wght-2.0d+0*rhoa*t152*wght)+t28*t
     ;         71*(8.0d+0*rhoa*rhob*t27*t55*t56*t70*wght+4.0d+0*rhoa*rho
     <         b*t101*t27*t69*wght)+t41*t73*(6.0d+0*rhoa*rhob*t102*t27*t
     =         70*wght+1.2d+1*rhoa*rhob*t27*t55*t56*t69*wght)+t28*t41*(-
     >         4.0d+0*rhoa*rhob*t55*t56*t68*wght-4.0d+0*rhoa*rhob*t155*t
     ?         27*t56*wght-4.0d+0*rhoa*rhob*t156*t27*t55*wght-t118*t70-2
     @         .0d+0*t120*t69)+t26*t28*(2.0d+0*rhoa*t55*t68*wght+2.0d+0*
     1         rhoa*rhob*t117*t68*wght+2.0d+0*rhoa*rhob*t155*t57*wght+2.
     2         0d+0*rhoa*rhob*t152*t55*wght+2.0d+0*rhoa*rhob*t221*t27*wg
     3         ht+2.0d+0*rhoa*t155*t27*wght-2.0d+0*t116*t69-2.0d+0*t155*
     4         t58-2.0d+0*t153*t55)+t26*t73*(-6.0d+0*rhoa*rhob*t102*t68*
     5         wght-1.2d+1*rhoa*rhob*t155*t27*t55*wght-3.0d+0*t118*t69)+
     6         t24*t71*(-2.0d+0*rhoa*rhob*t101*t68*wght-4.0d+0*rhoa*rhob
     7         *t156*t27*t56*wght-2*t120*t70)+6.0d+0*rhoa*rhob*t101*t193
     8         *t24*t27*t70*wght+2.4d+1*rhoa*rhob*t102*t194*t26*t27*t69*
     9         wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = t24*t26*(-1.0d+0*rhoa*rhob*(6.0856
     1         d+0*t1*t125*t18*t19*t29*t61*t76*t9+3.0428d+0*t1*t165*t18*
     2         t19*t29*t38*t76*t9+1.0142666666666666d+0*t1*t159*t18*t19*
     3         t76*t9-1.5214d+0*t1*t10*t18*t19*t223*t39*t9-5.07133333333
     4         3333d-1*t1*t121*t165*t18*t19*t39*t9-9.128400000000001d+0*
     5         rhoa*t1*t159*t18*t19*t196*t38*t9+6.0856d+0*t1*t159*t18*t1
     6         9*t29*t38*t76+3.0428d+0*rhoa*t1*t159*t18*t19*t37*t76-3.04
     7         28d+0*t1*t124*t18*t19*t29*t39*t61-3.0428d+0*t1*t10*t125*t
     8         18*t19*t39*t61-1.5214d+0*t1*t10*t165*t18*t19*t38*t39-1.52
     9         14d+0*t1*t165*t18*t19*t29*t37*t39-1.0142666666666666d+0*t
     :         1*t121*t159*t18*t19*t39+5.071333333333333d-1*t1*t165*t18*
     ;         t19*t21*t30+1.5214d+0*t1*t10*t18*t19*t21*t222)*wght-1.0d+
     <         0*rhob*t166*wght)+t24*t41*(1.0d+0*rhoa*rhob*t27*(4.513600
     =         0000000004d+0*t1*t125*t18*t19*t29*t61*t76*t9+4.0d+0*rhoa*
     >         t1*t159*t18*t19*t40*t76*t9+2.2568000000000002d+0*t1*t165*
     ?         t18*t19*t29*t38*t76*t9+7.522666666666666d-1*t1*t159*t18*t
     @         19*t76*t9-4.0d+0*t1*t129*t18*t19*t29*t39*t61*t9-2.0d+0*t1
     1         *t165*t18*t19*t29*t39*t40*t9-1.1284000000000001d+0*t1*t10
     2         *t18*t19*t223*t39*t9-3.761333333333333d-1*t1*t121*t165*t1
     3         8*t19*t39*t9-6.770400000000001d+0*rhoa*t1*t159*t18*t19*t1
     4         96*t38*t9+2.0d+0*t1*t10*t18*t19*t21*t224*t9+4.51360000000
     5         00004d+0*t1*t159*t18*t19*t29*t38*t76+2.2568000000000002d+
     6         0*rhoa*t1*t159*t18*t19*t37*t76-2.2568000000000002d+0*t1*t
     7         124*t18*t19*t29*t39*t61-2.2568000000000002d+0*t1*t10*t125
     8         *t18*t19*t39*t61+4.0d+0*t1*t10*t129*t18*t19*t21*t61-4.0d+
     9         0*t1*t159*t18*t19*t29*t39*t40+2.0d+0*t1*t10*t165*t18*t19*
     :         t21*t40-1.1284000000000001d+0*t1*t10*t165*t18*t19*t38*t39
     ;         -1.1284000000000001d+0*t1*t165*t18*t19*t29*t37*t39-7.5226
     <         66666666666d-1*t1*t121*t159*t18*t19*t39+3.761333333333333
     =         d-1*t1*t165*t18*t19*t21*t30+1.1284000000000001d+0*t1*t10*
     >         t18*t19*t21*t222)*wght+2.0d+0*rhoa*rhob*t131*t62*wght+1.0
     ?         d+0*rhoa*rhob*t166*t42*wght-2*t127*t64-t171*t44)+t28*t71*
     @         (4.0d+0*rhoa*rhob*t27*t42*t63*t64*wght-2*t128*t64)+t41*t7
     1         3*(6.0d+0*rhoa*rhob*t27*t40*t63*t64*wght-3.0d+0*t128*t63)
     2         +t24*t71*(-2.0d+0*rhoa*rhob*t42*t62*t64*wght-2.0d+0*rhoa*
     3         rhob*t131*t27*t64*wght-2.0d+0*rhoa*rhob*t171*t27*t42*wght
     4         -2*t132*t64)+t28*t41*(-2.0d+0*rhoa*rhob*t40*t62*t64*wght-
     5         2.0d+0*rhoa*rhob*t129*t27*t64*wght-2.0d+0*rhoa*rhob*t42*t
     6         62*t63*wght-2.0d+0*rhoa*rhob*t131*t27*t63*wght-2.0d+0*rho
     7         a*rhob*t169*t27*t42*wght-2.0d+0*rhoa*rhob*t171*t27*t40*wg
     8         ht-t130*t64-2.0d+0*t132*t63)+t26*t73*(-6.0d+0*rhoa*rhob*t
     9         40*t62*t63*wght-6.0d+0*rhoa*rhob*t129*t27*t63*wght-6.0d+0
     :         *rhoa*rhob*t169*t27*t40*wght-3.0d+0*t130*t63)+t26*t28*(4.
     ;         0d+0*rhoa*rhob*t129*t62*wght+2.0d+0*rhoa*rhob*t166*t40*wg
     <         ht+2.0d+0*rhoa*rhob*t224*t27*wght-4.0d+0*t127*t63-2.0d+0*
     =         t169*t44)+6.0d+0*rhoa*rhob*t167*t193*t24*t27*t42*wght+2.4
     >         d+1*rhoa*rhob*t168*t194*t26*t27*t40*wght+Cmat3(iq,D3_RA_G
     ?         AA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0
     1         d+0*t1*t18*t29*t40*t45*t61*t67*t76*t9+7.522666666666666d-
     2         1*t1*t121*t18*t45*t61*t67*t76*t9+2.2568000000000002d+0*t1
     3         *t19*t29*t38*t61*t67*t76*t9+2.2568000000000002d+0*t1*t10*
     4         t125*t18*t45*t67*t76*t9-6.770400000000001d+0*t1*t18*t196*
     5         t29*t38*t45*t61*t67*t9-2.0d+0*t1*t19*t29*t39*t40*t61*t67*
     6         t9-3.761333333333333d-1*t1*t121*t19*t39*t61*t67*t9-2.0d+0
     7         *t1*t10*t129*t18*t39*t45*t67*t9-1.1284000000000001d+0*t1*
     8         t10*t125*t19*t39*t67*t9+2.0d+0*t1*t10*t129*t19*t21*t67*t9
     9         -2.0d+0*t1*t136*t18*t19*t29*t39*t61*t9+2.0d+0*t1*t10*t18*
     :         t19*t21*t225*t9+2.2568000000000002d+0*t1*t10*t18*t38*t45*
     ;         t61*t67*t76+2.2568000000000002d+0*t1*t18*t29*t37*t45*t61*
     <         t67*t76-2.0d+0*t1*t10*t18*t39*t40*t45*t61*t67-3.761333333
     =         333333d-1*t1*t18*t30*t39*t45*t61*t67+2.0d+0*t1*t10*t19*t2
     >         1*t40*t61*t67-1.1284000000000001d+0*t1*t10*t19*t38*t39*t6
     ?         1*t67-1.1284000000000001d+0*t1*t19*t29*t37*t39*t61*t67+3.
     @         761333333333333d-1*t1*t19*t21*t30*t61*t67-1.1284000000000
     1         001d+0*t1*t10*t124*t18*t39*t45*t67+1.1284000000000001d+0*
     2         t1*t10*t124*t19*t21*t67+2.0d+0*t1*t10*t136*t18*t19*t21*t6
     3         1)*wght+1.0d+0*rhoa*rhob*t131*t68*wght+1.0d+0*rhoa*rhob*t
     4         137*t62*wght+1.0d+0*rhoa*rhob*t173*t42*wght-t127*t70-t134
     5         *t64-t177*t44)+t24*t26*(-1.0d+0*rhoa*rhob*(1.014266666666
     6         6666d+0*t1*t121*t18*t45*t61*t67*t76*t9+3.0428d+0*t1*t19*t
     7         29*t38*t61*t67*t76*t9+3.0428d+0*t1*t10*t125*t18*t45*t67*t
     8         76*t9-9.128400000000001d+0*t1*t18*t196*t29*t38*t45*t61*t6
     9         7*t9-5.071333333333333d-1*t1*t121*t19*t39*t61*t67*t9-1.52
     :         14d+0*t1*t10*t125*t19*t39*t67*t9+3.0428d+0*t1*t10*t18*t38
     ;         *t45*t61*t67*t76+3.0428d+0*t1*t18*t29*t37*t45*t61*t67*t76
     <         -5.071333333333333d-1*t1*t18*t30*t39*t45*t61*t67-1.5214d+
     =         0*t1*t10*t19*t38*t39*t61*t67-1.5214d+0*t1*t19*t29*t37*t39
     >         *t61*t67+5.071333333333333d-1*t1*t19*t21*t30*t61*t67-1.52
     ?         14d+0*t1*t10*t124*t18*t39*t45*t67+1.5214d+0*t1*t10*t124*t
     @         19*t21*t67)*wght-1.0d+0*rhob*t173*wght)+t41*t73*(6.0d+0*r
     1         hoa*rhob*t27*t40*t63*t70*wght-3.0d+0*t128*t69)+t28*t71*(4
     2         .0d+0*rhoa*rhob*t27*t42*t64*t69*wght-2*t128*t70)+t24*t71*
     3         (-2.0d+0*rhoa*rhob*t42*t64*t68*wght-2.0d+0*rhoa*rhob*t137
     4         *t27*t64*wght-2.0d+0*rhoa*rhob*t177*t27*t42*wght-2*t132*t
     5         70)+t28*t41*(-2.0d+0*rhoa*rhob*t40*t64*t68*wght-2.0d+0*rh
     6         oa*rhob*t42*t63*t68*wght-2.0d+0*rhoa*rhob*t136*t27*t64*wg
     7         ht-2.0d+0*rhoa*rhob*t137*t27*t63*wght-2.0d+0*rhoa*rhob*t1
     8         74*t27*t42*wght-2.0d+0*rhoa*rhob*t177*t27*t40*wght-t130*t
     9         70-2.0d+0*t132*t69)+t26*t73*(-6.0d+0*rhoa*rhob*t40*t63*t6
     :         8*wght-6.0d+0*rhoa*rhob*t136*t27*t63*wght-6.0d+0*rhoa*rho
     ;         b*t174*t27*t40*wght-3.0d+0*t130*t69)+t26*t28*(2.0d+0*rhoa
     <         *rhob*t129*t68*wght+2.0d+0*rhoa*rhob*t136*t62*wght+2.0d+0
     =         *rhoa*rhob*t173*t40*wght+2.0d+0*rhoa*rhob*t225*t27*wght-2
     >         .0d+0*t127*t69-2.0d+0*t134*t63-2.0d+0*t174*t44)+6.0d+0*rh
     ?         oa*rhob*t193*t24*t27*t42*t64*t70*wght+2.4d+1*rhoa*rhob*t1
     @         94*t26*t27*t40*t63*t69*wght+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t24*t26*(-1.0d+0*rhoa*rhob*(6.0856
     1         d+0*t1*t10*t187*t38*t45*t76*t9+3.0428d+0*t1*t10*t18*t186*
     2         t38*t45*t76*t9+1.0142666666666666d+0*rhob*t1*t18*t187*t30
     3         *t76*t9-1.0142666666666666d+0*t1*t187*t30*t39*t45*t9-5.07
     4         1333333333333d-1*t1*t18*t186*t30*t39*t45*t9-1.5214d+0*t1*
     5         t10*t186*t19*t38*t39*t9-9.128400000000001d+0*rhob*t1*t10*
     6         t18*t187*t196*t38*t9+5.071333333333333d-1*t1*t186*t19*t21
     7         *t30*t9+3.0428d+0*rhob*t1*t10*t18*t187*t37*t76-3.0428d+0*
     8         t1*t10*t187*t37*t39*t45-1.5214d+0*t1*t10*t18*t186*t37*t39
     9         *t45+1.5214d+0*t1*t10*t186*t19*t21*t37)*wght-1.0d+0*rhob*
     :         t188*wght)+t24*t41*(1.0d+0*rhoa*rhob*t27*(4.5136000000000
     ;         004d+0*t1*t10*t187*t38*t45*t76*t9+2.2568000000000002d+0*t
     <         1*t10*t18*t186*t38*t45*t76*t9+4.0d+0*rhob*t1*t10*t18*t187
     =         *t40*t76*t9+7.522666666666666d-1*rhob*t1*t18*t187*t30*t76
     >         *t9-4.0d+0*t1*t10*t136*t18*t39*t45*t67*t9+4.0d+0*t1*t10*t
     ?         136*t19*t21*t67*t9-4.0d+0*t1*t10*t187*t39*t40*t45*t9-2.0d
     @         +0*t1*t10*t18*t186*t39*t40*t45*t9-7.522666666666666d-1*t1
     1         *t187*t30*t39*t45*t9-3.761333333333333d-1*t1*t18*t186*t30
     2         *t39*t45*t9+2.0d+0*t1*t10*t186*t19*t21*t40*t9-1.128400000
     3         0000001d+0*t1*t10*t186*t19*t38*t39*t9-6.770400000000001d+
     4         0*rhob*t1*t10*t18*t187*t196*t38*t9+3.761333333333333d-1*t
     5         1*t186*t19*t21*t30*t9+2.0d+0*t1*t10*t18*t19*t21*t226*t9+2
     6         .2568000000000002d+0*rhob*t1*t10*t18*t187*t37*t76-2.25680
     7         00000000002d+0*t1*t10*t187*t37*t39*t45-1.1284000000000001
     8         d+0*t1*t10*t18*t186*t37*t39*t45+1.1284000000000001d+0*t1*
     9         t10*t186*t19*t21*t37)*wght+2.0d+0*rhoa*rhob*t137*t68*wght
     :         +1.0d+0*rhoa*rhob*t188*t42*wght-2*t134*t70-t191*t44)+t28*
     ;         t71*(4.0d+0*rhoa*rhob*t27*t42*t69*t70*wght-2*t135*t70)+t4
     <         1*t73*(6.0d+0*rhoa*rhob*t27*t40*t69*t70*wght-3.0d+0*t135*
     =         t69)+t24*t71*(-2.0d+0*rhoa*rhob*t42*t68*t70*wght-2.0d+0*r
     >         hoa*rhob*t137*t27*t70*wght-2.0d+0*rhoa*rhob*t191*t27*t42*
     ?         wght-2*t138*t70)+t28*t41*(-2.0d+0*rhoa*rhob*t40*t68*t70*w
     @         ght-2.0d+0*rhoa*rhob*t136*t27*t70*wght-2.0d+0*rhoa*rhob*t
     1         42*t68*t69*wght-2.0d+0*rhoa*rhob*t137*t27*t69*wght-2.0d+0
     2         *rhoa*rhob*t189*t27*t42*wght-2.0d+0*rhoa*rhob*t191*t27*t4
     3         0*wght-t139*t70-2.0d+0*t138*t69)+t26*t73*(-6.0d+0*rhoa*rh
     4         ob*t40*t68*t69*wght-6.0d+0*rhoa*rhob*t136*t27*t69*wght-6.
     5         0d+0*rhoa*rhob*t189*t27*t40*wght-3.0d+0*t139*t69)+t26*t28
     6         *(4.0d+0*rhoa*rhob*t136*t68*wght+2.0d+0*rhoa*rhob*t188*t4
     7         0*wght+2.0d+0*rhoa*rhob*t226*t27*wght-4.0d+0*t134*t69-2.0
     8         d+0*t189*t44)+6.0d+0*rhoa*rhob*t179*t193*t24*t27*t42*wght
     9         +2.4d+1*rhoa*rhob*t180*t194*t26*t27*t40*wght+Cmat3(iq,D3_
     :         RA_GBB_GBB)
            Cmat3(iq,D3_RB_GAA_GAA) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0
     1         d+0*rhoa*t1*t159*t18*t19*t55*t76*t9+2.2568000000000002d+0
     2         *t1*t165*t18*t19*t29*t54*t76*t9+2.2568000000000002d+0*rho
     3         a*t1*t159*t19*t53*t76*t9+7.522666666666666d-1*rhoa*t1*t15
     4         9*t18*t46*t76*t9-4.0d+0*t1*t143*t18*t19*t29*t39*t61*t9-2.
     5         0d+0*t1*t165*t18*t19*t29*t39*t55*t9-6.770400000000001d+0*
     6         rhoa*t1*t159*t18*t19*t196*t54*t9-1.1284000000000001d+0*t1
     7         *t165*t19*t29*t39*t53*t9-3.761333333333333d-1*t1*t165*t18
     8         *t29*t39*t46*t9+2.0d+0*t1*t10*t18*t19*t21*t227*t9+4.51360
     9         00000000004d+0*t1*t159*t18*t19*t29*t54*t76+4.0d+0*t1*t10*
     :         t143*t18*t19*t21*t61-4.0d+0*t1*t159*t18*t19*t29*t39*t55+2
     ;         .0d+0*t1*t10*t165*t18*t19*t21*t55-1.1284000000000001d+0*t
     <         1*t10*t165*t18*t19*t39*t54-2.2568000000000002d+0*t1*t159*
     =         t19*t29*t39*t53+1.1284000000000001d+0*t1*t10*t165*t19*t21
     >         *t53-7.522666666666666d-1*t1*t159*t18*t29*t39*t46+3.76133
     ?         3333333333d-1*t1*t10*t165*t18*t21*t46)*wght+2.0d+0*rhoa*r
     @         hob*t145*t62*wght+1.0d+0*rhoa*rhob*t166*t56*wght-2*t141*t
     1         64-t171*t58)+t24*t26*(-1.0d+0*rhoa*rhob*(3.0428d+0*t1*t16
     2         5*t18*t19*t29*t54*t76*t9+3.0428d+0*rhoa*t1*t159*t19*t53*t
     3         76*t9+1.0142666666666666d+0*rhoa*t1*t159*t18*t46*t76*t9-9
     4         .128400000000001d+0*rhoa*t1*t159*t18*t19*t196*t54*t9-1.52
     5         14d+0*t1*t165*t19*t29*t39*t53*t9-5.071333333333333d-1*t1*
     6         t165*t18*t29*t39*t46*t9+6.0856d+0*t1*t159*t18*t19*t29*t54
     7         *t76-1.5214d+0*t1*t10*t165*t18*t19*t39*t54-3.0428d+0*t1*t
     8         159*t19*t29*t39*t53+1.5214d+0*t1*t10*t165*t19*t21*t53-1.0
     9         142666666666666d+0*t1*t159*t18*t29*t39*t46+5.071333333333
     :         333d-1*t1*t10*t165*t18*t21*t46)*wght-1.0d+0*rhoa*t166*wgh
     ;         t)+t28*t71*(4.0d+0*rhoa*rhob*t27*t56*t63*t64*wght-2*t142*
     <         t64)+t41*t73*(6.0d+0*rhoa*rhob*t27*t55*t63*t64*wght-3.0d+
     =         0*t142*t63)+t24*t71*(-2.0d+0*rhoa*rhob*t56*t62*t64*wght-2
     >         .0d+0*rhoa*rhob*t145*t27*t64*wght-2.0d+0*rhoa*rhob*t171*t
     ?         27*t56*wght-2*t146*t64)+t28*t41*(-2.0d+0*rhoa*rhob*t55*t6
     @         2*t64*wght-2.0d+0*rhoa*rhob*t143*t27*t64*wght-2.0d+0*rhoa
     1         *rhob*t56*t62*t63*wght-2.0d+0*rhoa*rhob*t145*t27*t63*wght
     2         -2.0d+0*rhoa*rhob*t169*t27*t56*wght-2.0d+0*rhoa*rhob*t171
     3         *t27*t55*wght-t144*t64-2.0d+0*t146*t63)+t26*t73*(-6.0d+0*
     4         rhoa*rhob*t55*t62*t63*wght-6.0d+0*rhoa*rhob*t143*t27*t63*
     5         wght-6.0d+0*rhoa*rhob*t169*t27*t55*wght-3.0d+0*t144*t63)+
     6         t26*t28*(4.0d+0*rhoa*rhob*t143*t62*wght+2.0d+0*rhoa*rhob*
     7         t166*t55*wght+2.0d+0*rhoa*rhob*t227*t27*wght-4.0d+0*t141*
     8         t63-2.0d+0*t169*t58)+6.0d+0*rhoa*rhob*t167*t193*t24*t27*t
     9         56*wght+2.4d+1*rhoa*rhob*t168*t194*t26*t27*t55*wght+Cmat3
     :         (iq,D3_RB_GAA_GAA)
            Cmat3(iq,D3_RB_GAA_GAB) = Cmat3(iq,D3_RB_GAA_GAB)
            Cmat3(iq,D3_RB_GAA_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.0
     1         d+0*t1*t18*t29*t45*t55*t61*t67*t76*t9+2.2568000000000002d
     2         +0*t1*t19*t29*t54*t61*t67*t76*t9+2.2568000000000002d+0*t1
     3         *t29*t45*t53*t61*t67*t76*t9+7.522666666666666d-1*t1*t147*
     4         t18*t29*t61*t67*t76*t9+2.2568000000000002d+0*t1*t151*t18*
     5         t19*t29*t61*t76*t9-2.0d+0*t1*t19*t29*t39*t55*t61*t67*t9-6
     6         .770400000000001d+0*t1*t18*t196*t29*t45*t54*t61*t67*t9-3.
     7         761333333333333d-1*t1*t29*t39*t46*t61*t67*t9-2.0d+0*t1*t1
     8         0*t143*t18*t39*t45*t67*t9+2.0d+0*t1*t10*t143*t19*t21*t67*
     9         t9-2.0d+0*t1*t155*t18*t19*t29*t39*t61*t9-1.12840000000000
     :         01d+0*t1*t150*t19*t29*t39*t61*t9+2.0d+0*t1*t10*t18*t19*t2
     ;         1*t228*t9+2.2568000000000002d+0*t1*t10*t18*t45*t54*t61*t6
     <         7*t76-2.0d+0*t1*t10*t18*t39*t45*t55*t61*t67+2.0d+0*t1*t10
     =         *t19*t21*t55*t61*t67-1.1284000000000001d+0*t1*t10*t19*t39
     >         *t54*t61*t67-1.1284000000000001d+0*t1*t10*t39*t45*t53*t61
     ?         *t67+3.761333333333333d-1*t1*t10*t21*t46*t61*t67-3.761333
     @         333333333d-1*t1*t10*t147*t18*t39*t61*t67-1.12840000000000
     1         01d+0*t1*t10*t151*t18*t19*t39*t61+2.0d+0*t1*t10*t155*t18*
     2         t19*t21*t61+1.1284000000000001d+0*t1*t10*t150*t19*t21*t61
     3         )*wght+1.0d+0*rhoa*rhob*t145*t68*wght+1.0d+0*rhoa*rhob*t1
     4         56*t62*wght+1.0d+0*rhoa*rhob*t173*t56*wght-t141*t70-t153*
     5         t64-t177*t58)+t24*t26*(-1.0d+0*rhoa*rhob*(3.0428d+0*t1*t1
     6         9*t29*t54*t61*t67*t76*t9+3.0428d+0*t1*t29*t45*t53*t61*t67
     7         *t76*t9+1.0142666666666666d+0*t1*t147*t18*t29*t61*t67*t76
     8         *t9+3.0428d+0*t1*t151*t18*t19*t29*t61*t76*t9-9.1284000000
     9         00001d+0*t1*t18*t196*t29*t45*t54*t61*t67*t9-5.07133333333
     :         3333d-1*t1*t29*t39*t46*t61*t67*t9-1.5214d+0*t1*t150*t19*t
     ;         29*t39*t61*t9+3.0428d+0*t1*t10*t18*t45*t54*t61*t67*t76-1.
     <         5214d+0*t1*t10*t19*t39*t54*t61*t67-1.5214d+0*t1*t10*t39*t
     =         45*t53*t61*t67+5.071333333333333d-1*t1*t10*t21*t46*t61*t6
     >         7-5.071333333333333d-1*t1*t10*t147*t18*t39*t61*t67-1.5214
     ?         d+0*t1*t10*t151*t18*t19*t39*t61+1.5214d+0*t1*t10*t150*t19
     @         *t21*t61)*wght-1.0d+0*rhoa*t173*wght)+t41*t73*(6.0d+0*rho
     1         a*rhob*t27*t55*t63*t70*wght-3.0d+0*t142*t69)+t28*t71*(4.0
     2         d+0*rhoa*rhob*t27*t56*t64*t69*wght-2*t142*t70)+t24*t71*(-
     3         2.0d+0*rhoa*rhob*t56*t64*t68*wght-2.0d+0*rhoa*rhob*t156*t
     4         27*t64*wght-2.0d+0*rhoa*rhob*t177*t27*t56*wght-2*t146*t70
     5         )+t28*t41*(-2.0d+0*rhoa*rhob*t55*t64*t68*wght-2.0d+0*rhoa
     6         *rhob*t56*t63*t68*wght-2.0d+0*rhoa*rhob*t155*t27*t64*wght
     7         -2.0d+0*rhoa*rhob*t156*t27*t63*wght-2.0d+0*rhoa*rhob*t174
     8         *t27*t56*wght-2.0d+0*rhoa*rhob*t177*t27*t55*wght-t144*t70
     9         -2.0d+0*t146*t69)+t26*t73*(-6.0d+0*rhoa*rhob*t55*t63*t68*
     :         wght-6.0d+0*rhoa*rhob*t155*t27*t63*wght-6.0d+0*rhoa*rhob*
     ;         t174*t27*t55*wght-3.0d+0*t144*t69)+t26*t28*(2.0d+0*rhoa*r
     <         hob*t143*t68*wght+2.0d+0*rhoa*rhob*t155*t62*wght+2.0d+0*r
     =         hoa*rhob*t173*t55*wght+2.0d+0*rhoa*rhob*t228*t27*wght-2.0
     >         d+0*t141*t69-2.0d+0*t153*t63-2.0d+0*t174*t58)+6.0d+0*rhoa
     ?         *rhob*t193*t24*t27*t56*t64*t70*wght+2.4d+1*rhoa*rhob*t194
     @         *t26*t27*t55*t63*t69*wght+Cmat3(iq,D3_RB_GAA_GBB)
            Cmat3(iq,D3_RB_GAB_GAB) = Cmat3(iq,D3_RB_GAB_GAB)
            Cmat3(iq,D3_RB_GAB_GBB) = Cmat3(iq,D3_RB_GAB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t24*t26*(-1.0d+0*rhoa*rhob*(6.0856
     1         d+0*t1*t10*t151*t18*t45*t67*t76*t9+6.0856d+0*t1*t10*t187*
     2         t45*t54*t76*t9+3.0428d+0*t1*t10*t18*t186*t45*t54*t76*t9+3
     3         .0428d+0*rhob*t1*t10*t187*t53*t76*t9+1.0142666666666666d+
     4         0*t1*t10*t18*t187*t76*t9-3.0428d+0*t1*t10*t150*t39*t45*t6
     5         7*t9-3.0428d+0*t1*t10*t151*t19*t39*t67*t9-1.5214d+0*t1*t1
     6         0*t186*t19*t39*t54*t9-9.128400000000001d+0*rhob*t1*t10*t1
     7         8*t187*t196*t54*t9-1.5214d+0*t1*t10*t186*t39*t45*t53*t9+5
     8         .071333333333333d-1*t1*t10*t186*t21*t46*t9-1.5214d+0*t1*t
     9         10*t18*t19*t230*t39*t9-1.0142666666666666d+0*t1*t10*t147*
     :         t187*t39*t9-5.071333333333333d-1*t1*t10*t147*t18*t186*t39
     ;         *t9+1.5214d+0*t1*t10*t19*t21*t229*t9)*wght-1.0d+0*rhoa*t1
     <         88*wght)+t24*t41*(1.0d+0*rhoa*rhob*t27*(4.513600000000000
     =         4d+0*t1*t10*t151*t18*t45*t67*t76*t9+4.0d+0*rhob*t1*t10*t1
     >         8*t187*t55*t76*t9+4.5136000000000004d+0*t1*t10*t187*t45*t
     ?         54*t76*t9+2.2568000000000002d+0*t1*t10*t18*t186*t45*t54*t
     @         76*t9+2.2568000000000002d+0*rhob*t1*t10*t187*t53*t76*t9+7
     1         .522666666666666d-1*t1*t10*t18*t187*t76*t9-4.0d+0*t1*t10*
     2         t155*t18*t39*t45*t67*t9-2.2568000000000002d+0*t1*t10*t150
     3         *t39*t45*t67*t9-2.2568000000000002d+0*t1*t10*t151*t19*t39
     4         *t67*t9+4.0d+0*t1*t10*t155*t19*t21*t67*t9-4.0d+0*t1*t10*t
     5         187*t39*t45*t55*t9-2.0d+0*t1*t10*t18*t186*t39*t45*t55*t9+
     6         2.0d+0*t1*t10*t186*t19*t21*t55*t9-1.1284000000000001d+0*t
     7         1*t10*t186*t19*t39*t54*t9-6.770400000000001d+0*rhob*t1*t1
     8         0*t18*t187*t196*t54*t9-1.1284000000000001d+0*t1*t10*t186*
     9         t39*t45*t53*t9+3.761333333333333d-1*t1*t10*t186*t21*t46*t
     :         9-1.1284000000000001d+0*t1*t10*t18*t19*t230*t39*t9-7.5226
     ;         66666666666d-1*t1*t10*t147*t187*t39*t9-3.761333333333333d
     <         -1*t1*t10*t147*t18*t186*t39*t9+2.0d+0*t1*t10*t18*t19*t21*
     =         t231*t9+1.1284000000000001d+0*t1*t10*t19*t21*t229*t9)*wgh
     >         t+2.0d+0*rhoa*rhob*t156*t68*wght+1.0d+0*rhoa*rhob*t188*t5
     ?         6*wght-2*t153*t70-t191*t58)+t28*t71*(4.0d+0*rhoa*rhob*t27
     @         *t56*t69*t70*wght-2*t154*t70)+t41*t73*(6.0d+0*rhoa*rhob*t
     1         27*t55*t69*t70*wght-3.0d+0*t154*t69)+t24*t71*(-2.0d+0*rho
     2         a*rhob*t56*t68*t70*wght-2.0d+0*rhoa*rhob*t156*t27*t70*wgh
     3         t-2.0d+0*rhoa*rhob*t191*t27*t56*wght-2*t157*t70)+t28*t41*
     4         (-2.0d+0*rhoa*rhob*t55*t68*t70*wght-2.0d+0*rhoa*rhob*t155
     5         *t27*t70*wght-2.0d+0*rhoa*rhob*t56*t68*t69*wght-2.0d+0*rh
     6         oa*rhob*t156*t27*t69*wght-2.0d+0*rhoa*rhob*t189*t27*t56*w
     7         ght-2.0d+0*rhoa*rhob*t191*t27*t55*wght-t158*t70-2.0d+0*t1
     8         57*t69)+t26*t73*(-6.0d+0*rhoa*rhob*t55*t68*t69*wght-6.0d+
     9         0*rhoa*rhob*t155*t27*t69*wght-6.0d+0*rhoa*rhob*t189*t27*t
     :         55*wght-3.0d+0*t158*t69)+t26*t28*(4.0d+0*rhoa*rhob*t155*t
     ;         68*wght+2.0d+0*rhoa*rhob*t188*t55*wght+2.0d+0*rhoa*rhob*t
     <         231*t27*wght-4.0d+0*t153*t69-2.0d+0*t189*t58)+6.0d+0*rhoa
     =         *rhob*t179*t193*t24*t27*t56*wght+2.4d+1*rhoa*rhob*t180*t1
     >         94*t26*t27*t55*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.
     1         0d+0*rhoa*t1*t159*t18*t19*t63*t76*t9+6.770400000000001d+0
     2         *rhoa*t1*t165*t18*t19*t61*t76*t9-2.0d+0*t1*t165*t18*t19*t
     3         29*t39*t63*t9-4.0d+0*t1*t169*t18*t19*t29*t39*t61*t9-1.128
     4         4000000000001d+0*t1*t18*t19*t233*t29*t39*t9-6.77040000000
     5         0001d+0*t1*t18*t19*t196*t232*t3*t9+2.0d+0*t1*t10*t18*t19*
     6         t21*t234*t9+6.770400000000001d+0*rhoa*t1*t18*t19*t232*t76
     7         -4.0d+0*t1*t159*t18*t19*t29*t39*t63+2.0d+0*t1*t10*t165*t1
     8         8*t19*t21*t63-6.770400000000001d+0*t1*t165*t18*t19*t29*t3
     9         9*t61+4.0d+0*t1*t10*t169*t18*t19*t21*t61+1.12840000000000
     :         01d+0*t1*t10*t18*t19*t21*t233)*wght+3.0d+0*rhoa*rhob*t166
     ;         *t64*wght+3.0d+0*rhoa*rhob*t171*t62*wght)+t28*t41*(-4.0d+
     <         0*rhoa*rhob*t62*t63*t64*wght-4.0d+0*rhoa*rhob*t169*t27*t6
     =         4*wght-4.0d+0*rhoa*rhob*t171*t27*t63*wght-t170*t64-2.0d+0
     >         *t172*t63)+t24*t71*(-4.0d+0*rhoa*rhob*t171*t27*t64*wght-2
     ?         .0d+0*rhoa*rhob*t167*t62*wght-2*t172*t64)+t26*t73*(-1.2d+
     @         1*rhoa*rhob*t169*t27*t63*wght-6.0d+0*rhoa*rhob*t168*t62*w
     1         ght-3.0d+0*t170*t63)+t26*t28*(6.0d+0*rhoa*rhob*t166*t63*w
     2         ght+6.0d+0*rhoa*rhob*t169*t62*wght+2.0d+0*rhoa*rhob*t234*
     3         t27*wght)-1.0d+0*rhoa*rhob*t24*t26*(9.128400000000001d+0*
     4         rhoa*t1*t165*t18*t19*t61*t76*t9-1.5214d+0*t1*t18*t19*t233
     5         *t29*t39*t9-9.128400000000001d+0*t1*t18*t19*t196*t232*t3*
     6         t9+9.128400000000001d+0*rhoa*t1*t18*t19*t232*t76-9.128400
     7         000000001d+0*t1*t165*t18*t19*t29*t39*t61+1.5214d+0*t1*t10
     8         *t18*t19*t21*t233)*wght+1.8d+1*rhoa*rhob*t168*t27*t41*t64
     9         *t73*wght+1.2d+1*rhoa*rhob*t167*t27*t28*t63*t71*wght+6.0d
     :         +0*rhoa*rhob*t193*t24*t27*t64**3*wght+2.4d+1*rhoa*rhob*t1
     ;         94*t26*t27*t63**3*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.
     1         0d+0*t1*t18*t29*t45*t61*t63*t67*t76*t9+2.2568000000000002
     2         d+0*t1*t165*t18*t29*t45*t67*t76*t9+2.2568000000000002d+0*
     3         rhoa*t1*t159*t19*t67*t76*t9-2.0d+0*t1*t19*t29*t39*t61*t63
     4         *t67*t9-2.0d+0*t1*t10*t169*t18*t39*t45*t67*t9-6.770400000
     5         000001d+0*rhoa*t1*t159*t18*t196*t45*t67*t9-1.128400000000
     6         0001d+0*t1*t165*t19*t29*t39*t67*t9+2.0d+0*t1*t10*t169*t19
     7         *t21*t67*t9-2.0d+0*t1*t174*t18*t19*t29*t39*t61*t9+2.0d+0*
     8         t1*t10*t18*t19*t21*t235*t9+4.5136000000000004d+0*t1*t159*
     9         t18*t29*t45*t67*t76-2.0d+0*t1*t10*t18*t39*t45*t61*t63*t67
     :         +2.0d+0*t1*t10*t19*t21*t61*t63*t67-1.1284000000000001d+0*
     ;         t1*t10*t165*t18*t39*t45*t67-2.2568000000000002d+0*t1*t159
     <         *t19*t29*t39*t67+1.1284000000000001d+0*t1*t10*t165*t19*t2
     =         1*t67+2.0d+0*t1*t10*t174*t18*t19*t21*t61)*wght+1.0d+0*rho
     >         a*rhob*t166*t70*wght+1.0d+0*rhoa*rhob*t171*t68*wght+2.0d+
     ?         0*rhoa*rhob*t173*t64*wght+2.0d+0*rhoa*rhob*t177*t62*wght)
     @         +t28*t71*(8.0d+0*rhoa*rhob*t27*t63*t64*t70*wght+4.0d+0*rh
     1         oa*rhob*t167*t27*t69*wght)+t41*t73*(6.0d+0*rhoa*rhob*t168
     2         *t27*t70*wght+1.2d+1*rhoa*rhob*t27*t63*t64*t69*wght)+t26*
     3         t28*(2.0d+0*rhoa*rhob*t166*t69*wght+2.0d+0*rhoa*rhob*t169
     4         *t68*wght+4.0d+0*rhoa*rhob*t173*t63*wght+4.0d+0*rhoa*rhob
     5         *t174*t62*wght+2.0d+0*rhoa*rhob*t235*t27*wght)+t28*t41*(-
     6         4.0d+0*rhoa*rhob*t63*t64*t68*wght-4.0d+0*rhoa*rhob*t174*t
     7         27*t64*wght-4.0d+0*rhoa*rhob*t177*t27*t63*wght-t170*t70-2
     8         .0d+0*t172*t69)+t26*t73*(-6.0d+0*rhoa*rhob*t168*t68*wght-
     9         1.2d+1*rhoa*rhob*t174*t27*t63*wght-3.0d+0*t170*t69)+t24*t
     :         71*(-2.0d+0*rhoa*rhob*t167*t68*wght-4.0d+0*rhoa*rhob*t177
     ;         *t27*t64*wght-2*t172*t70)-1.0d+0*rhoa*rhob*t24*t26*(3.042
     <         8d+0*t1*t165*t18*t29*t45*t67*t76*t9+3.0428d+0*rhoa*t1*t15
     =         9*t19*t67*t76*t9-9.128400000000001d+0*rhoa*t1*t159*t18*t1
     >         96*t45*t67*t9-1.5214d+0*t1*t165*t19*t29*t39*t67*t9+6.0856
     ?         d+0*t1*t159*t18*t29*t45*t67*t76-1.5214d+0*t1*t10*t165*t18
     @         *t39*t45*t67-3.0428d+0*t1*t159*t19*t29*t39*t67+1.5214d+0*
     1         t1*t10*t165*t19*t21*t67)*wght+6.0d+0*rhoa*rhob*t167*t193*
     2         t24*t27*t70*wght+2.4d+1*rhoa*rhob*t168*t194*t26*t27*t69*w
     3         ght+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.
     1         0d+0*rhob*t1*t10*t18*t187*t63*t76*t9+4.5136000000000004d+
     2         0*t1*t187*t29*t45*t61*t76*t9+2.2568000000000002d+0*t1*t18
     3         *t186*t29*t45*t61*t76*t9-4.0d+0*t1*t10*t174*t18*t39*t45*t
     4         67*t9+4.0d+0*t1*t10*t174*t19*t21*t67*t9-4.0d+0*t1*t10*t18
     5         7*t39*t45*t63*t9-2.0d+0*t1*t10*t18*t186*t39*t45*t63*t9+2.
     6         0d+0*t1*t10*t186*t19*t21*t63*t9-1.1284000000000001d+0*t1*
     7         t186*t19*t29*t39*t61*t9-6.770400000000001d+0*rhob*t1*t18*
     8         t187*t196*t29*t61*t9+2.0d+0*t1*t10*t18*t19*t21*t236*t9+2.
     9         2568000000000002d+0*rhob*t1*t10*t18*t187*t61*t76-2.256800
     :         0000000002d+0*t1*t10*t187*t39*t45*t61-1.1284000000000001d
     ;         +0*t1*t10*t18*t186*t39*t45*t61+1.1284000000000001d+0*t1*t
     <         10*t186*t19*t21*t61)*wght+2.0d+0*rhoa*rhob*t173*t70*wght+
     =         2.0d+0*rhoa*rhob*t177*t68*wght+1.0d+0*rhoa*rhob*t188*t64*
     >         wght+1.0d+0*rhoa*rhob*t191*t62*wght)+t28*t71*(4.0d+0*rhoa
     ?         *rhob*t27*t64*t69*t70*wght-2*t176*t70)+t41*t73*(6.0d+0*rh
     @         oa*rhob*t27*t63*t69*t70*wght-3.0d+0*t176*t69)+t24*t71*(-2
     1         .0d+0*rhoa*rhob*t64*t68*t70*wght-2.0d+0*rhoa*rhob*t177*t2
     2         7*t70*wght-2.0d+0*rhoa*rhob*t191*t27*t64*wght-2*t178*t70)
     3         +t28*t41*(-2.0d+0*rhoa*rhob*t63*t68*t70*wght-2.0d+0*rhoa*
     4         rhob*t174*t27*t70*wght-2.0d+0*rhoa*rhob*t64*t68*t69*wght-
     5         2.0d+0*rhoa*rhob*t177*t27*t69*wght-2.0d+0*rhoa*rhob*t189*
     6         t27*t64*wght-2.0d+0*rhoa*rhob*t191*t27*t63*wght-t175*t70-
     7         2.0d+0*t178*t69)+t26*t73*(-6.0d+0*rhoa*rhob*t63*t68*t69*w
     8         ght-6.0d+0*rhoa*rhob*t174*t27*t69*wght-6.0d+0*rhoa*rhob*t
     9         189*t27*t63*wght-3.0d+0*t175*t69)+t26*t28*(4.0d+0*rhoa*rh
     :         ob*t173*t69*wght+4.0d+0*rhoa*rhob*t174*t68*wght+2.0d+0*rh
     ;         oa*rhob*t188*t63*wght+2.0d+0*rhoa*rhob*t189*t62*wght+2.0d
     <         +0*rhoa*rhob*t236*t27*wght)-1.0d+0*rhoa*rhob*t24*t26*(6.0
     =         856d+0*t1*t187*t29*t45*t61*t76*t9+3.0428d+0*t1*t18*t186*t
     >         29*t45*t61*t76*t9-1.5214d+0*t1*t186*t19*t29*t39*t61*t9-9.
     ?         128400000000001d+0*rhob*t1*t18*t187*t196*t29*t61*t9+3.042
     @         8d+0*rhob*t1*t10*t18*t187*t61*t76-3.0428d+0*t1*t10*t187*t
     1         39*t45*t61-1.5214d+0*t1*t10*t18*t186*t39*t45*t61+1.5214d+
     2         0*t1*t10*t186*t19*t21*t61)*wght+6.0d+0*rhoa*rhob*t179*t19
     3         3*t24*t27*t64*wght+2.4d+1*rhoa*rhob*t180*t194*t26*t27*t63
     4         *wght+Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t24*t41*(1.0d+0*rhoa*rhob*t27*(4.
     1         0d+0*rhob*t1*t10*t18*t187*t69*t76*t9+6.770400000000001d+0
     2         *rhob*t1*t10*t18*t186*t67*t76*t9+6.770400000000001d+0*rho
     3         b*t1*t10*t238*t76*t9-4.0d+0*t1*t10*t187*t39*t45*t69*t9-2.
     4         0d+0*t1*t10*t18*t186*t39*t45*t69*t9+2.0d+0*t1*t10*t186*t1
     5         9*t21*t69*t9-4.0d+0*t1*t10*t18*t189*t39*t45*t67*t9-6.7704
     6         00000000001d+0*t1*t10*t186*t39*t45*t67*t9+4.0d+0*t1*t10*t
     7         189*t19*t21*t67*t9-1.1284000000000001d+0*t1*t10*t18*t237*
     8         t39*t45*t9+2.0d+0*t1*t10*t18*t19*t21*t239*t9-6.7704000000
     9         00001d+0*t1*t10*t12*t18*t196*t238*t9+1.1284000000000001d+
     :         0*t1*t10*t19*t21*t237*t9)*wght+3.0d+0*rhoa*rhob*t188*t70*
     ;         wght+3.0d+0*rhoa*rhob*t191*t68*wght)+t28*t41*(-4.0d+0*rho
     <         a*rhob*t68*t69*t70*wght-4.0d+0*rhoa*rhob*t189*t27*t70*wgh
     =         t-4.0d+0*rhoa*rhob*t191*t27*t69*wght-t190*t70-2.0d+0*t192
     >         *t69)+t24*t71*(-4.0d+0*rhoa*rhob*t191*t27*t70*wght-2.0d+0
     ?         *rhoa*rhob*t179*t68*wght-2*t192*t70)+t26*t73*(-1.2d+1*rho
     @         a*rhob*t189*t27*t69*wght-6.0d+0*rhoa*rhob*t180*t68*wght-3
     1         .0d+0*t190*t69)+t26*t28*(6.0d+0*rhoa*rhob*t188*t69*wght+6
     2         .0d+0*rhoa*rhob*t189*t68*wght+2.0d+0*rhoa*rhob*t239*t27*w
     3         ght)-1.0d+0*rhoa*rhob*t24*t26*(9.128400000000001d+0*rhob*
     4         t1*t10*t18*t186*t67*t76*t9+9.128400000000001d+0*rhob*t1*t
     5         10*t238*t76*t9-9.128400000000001d+0*t1*t10*t186*t39*t45*t
     6         67*t9-1.5214d+0*t1*t10*t18*t237*t39*t45*t9-9.128400000000
     7         001d+0*t1*t10*t12*t18*t196*t238*t9+1.5214d+0*t1*t10*t19*t
     8         21*t237*t9)*wght+1.8d+1*rhoa*rhob*t180*t27*t41*t70*t73*wg
     9         ht+1.2d+1*rhoa*rhob*t179*t27*t28*t69*t71*wght+6.0d+0*rhoa
     :         *rhob*t193*t24*t27*t70**3*wght+2.4d+1*rhoa*rhob*t194*t26*
     ;         t27*t69**3*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)
            Amat(iq,D1_RB) = Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
