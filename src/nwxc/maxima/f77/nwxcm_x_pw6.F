C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_pw6.F 
C> The nwxcm_x_pw6 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_pw6 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}^{{\it t_6}}}}\\\\ 
C>   {\it t_8} &=& {{{\it t_6}}\over{2}}\\\\ 
C>   {\it t_9} &=& \sigma_{\alpha\alpha}^{{\it t_8}}\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_3}-0.001890381166699926\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{13}} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^{4.0}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{14}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{18}} &=& \sigma_{\beta\beta}^{{\it t_8}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{13}}^{8.0}}}\\\\ 
C>   {\it t_{20}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{21}} &=& {\it t_{20}}^{4.0}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{{\it t_{21}}}}\\\\ 
C>   {\it t_{23}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{24}} &=& {{1}\over{{\it t_{21}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{25}} &=& \sigma_{ss}^{{\it t_8}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{20}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{
C>    \it t_{18}}-{\it t_3}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}
C>    \,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{
C>    \it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{
C>    \it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,
C>    \sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}
C>    \,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,
C>    \sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}
C>    -{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}
C>    \over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{
C>    \it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; 
C>    \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="https://doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pw6(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = 1.0d+0*rhoa**3.333333333333333d-1
            t5 = t4**4.0d+0
            t6 = 1/t5**t1
            t7 = param(1)
            t8 = gammaaa**5.0d-1
            t9 = 1/t5
            t10 = asinh(t8*t9)
            t11 = 6.0d+0*t10*t7*t8*t9+1.0746613026776465d-6*t3*t6+1.0d+0
            t12 = 1/t11
            t13 = 1/t4**8.0d+0
            t14 = t7-1.890381166699926d-3
            t15 = param(2)
            t16 = exp(-gammaaa*t13*t15)
            t17 = -gammaaa*t13*t7+1.0d-6*t3*t6+gammaaa*t13*t14*t16
            t18 = 1/t11**2
            t19 = 1/rhoa
            t20 = 1/(gammaaa*t13+1)**5.0d-1
            t21 = 1/t4**9.0d+0
            t22 = 1/rhoa**6.666666666666666d-1
            t23 = gammaaa**(t2-1)
            fnc(iq) = 2.0d+0*t12*t17*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.0d+0*t17*t18*t5*(-8.0d+0*t10*t22*t7*t8/
     1         t4**5.0d+0-8.0d+0*gammaaa*t20*t21*t22*t7-1.43288173690352
     2         88d-6*t1*t19*t3*t6)+1.0d+0*t12*t5*(2.6666666666666666d+0*
     3         gammaaa*t21*t22*t7-1.3333333333333333d-6*t1*t19*t3*t6+2.6
     4         666666666666666d+0*gammaaa**2*t14*t15*t16*t22/t4**1.7d+1-
     5         2.6666666666666666d+0*gammaaa*t14*t16*t21*t22)+1.33333333
     6         33333333d+0*t12*t17*t22*t4**3.0d+0)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t12*t5*(-t13*t7+5.0d-7*t1*t23*t6-g
     1         ammaaa*t14*t15*t16/t4**1.6d+1+t13*t14*t16)-1.0d+0*t17*t18
     2         *t5*(3.0d+0*t10*t7*t9/t8+3.0d+0*t13*t20*t7+5.373306513388
     3         233d-7*t1*t23*t6))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = gammabb**t2
            t20 = rhob**3.333333333333333d-1
            t21 = 1.0d+0*t20
            t22 = t21**4.0d+0
            t23 = 1/t22**t1
            t24 = gammabb**5.0d-1
            t25 = 1/t22
            t26 = asinh(t24*t25)
            t27 = 6.0d+0*t24*t25*t26*t8+1.0746613026776465d-6*t19*t23+1.
     1         0d+0
            t28 = 1/t27
            t29 = 1/t21**8.0d+0
            t30 = exp(-gammabb*t16*t29)
            t31 = -gammabb*t29*t8+gammabb*t15*t29*t30+1.0d-6*t19*t23
            t32 = 1/t12**2
            t33 = 1/(gammaaa*t14+1)**5.0d-1
            t34 = 1/t5**9.0d+0
            t35 = 1/rhoa**6.666666666666666d-1
            t36 = 1.0d+0/t4
            t37 = 1/t27**2
            t38 = 1/(gammabb*t29+1)**5.0d-1
            t39 = 1/t21**9.0d+0
            t40 = 1/rhob**6.666666666666666d-1
            t41 = 1.0d+0/t20
            t42 = t2-1
            t43 = gammaaa**t42
            t44 = gammabb**t42
            fnc(iq) = (1.0d+0*t13*t18*t6+1.0d+0*t22*t28*t31)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = (-1.0d+0*t18*t32*t6*(-8.0d+0*t11*t35*t8*t9/
     1         t5**5.0d+0-8.0d+0*gammaaa*t33*t34*t35*t8-1.43288173690352
     2         88d-6*t1*t3*t35*t36*t7)+1.0d+0*t13*t6*(2.6666666666666666
     3         d+0*gammaaa*t34*t35*t8-1.3333333333333333d-6*t1*t3*t35*t3
     4         6*t7+2.6666666666666666d+0*gammaaa**2*t15*t16*t17*t35/t5*
     5         *1.7d+1-2.6666666666666666d+0*gammaaa*t15*t17*t34*t35)+1.
     6         3333333333333333d+0*t13*t18*t35*t5**3.0d+0)*wght+Amat(iq,
     7         D1_RA)
            Amat(iq,D1_RB) = (-1.0d+0*t22*t31*t37*(-8.0d+0*gammabb*t38*t
     1         39*t40*t8-8.0d+0*t24*t26*t40*t8/t21**5.0d+0-1.43288173690
     2         35288d-6*t1*t19*t23*t40*t41)+1.0d+0*t22*t28*(2.6666666666
     3         666666d+0*gammabb*t39*t40*t8-1.3333333333333333d-6*t1*t19
     4         *t23*t40*t41-2.6666666666666666d+0*gammabb*t15*t30*t39*t4
     5         0+2.6666666666666666d+0*gammabb**2*t15*t16*t30*t40/t21**1
     6         .7d+1)+1.3333333333333333d+0*t21**3.0d+0*t28*t31*t40)*wgh
     7         t+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t6*(-t14*t8+5.0d-7*t1*t43*t7-g
     1         ammaaa*t15*t16*t17/t5**1.6d+1+t14*t15*t17)-1.0d+0*t18*t32
     2         *t6*(3.0d+0*t10*t11*t8/t9+3.0d+0*t14*t33*t8+5.37330651338
     3         8233d-7*t1*t43*t7))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t22*t28*(-t29*t8+5.0d-7*t1*t23*t44
     1         +t15*t29*t30-gammabb*t15*t16*t30/t21**1.6d+1)-1.0d+0*t22*
     2         t31*t37*(3.0d+0*t29*t38*t8+3.0d+0*t25*t26*t8/t24+5.373306
     3         513388233d-7*t1*t23*t44))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = 1/t12**2
            t20 = 1/(gammaaa*t14+1)**5.0d-1
            t21 = 1/t5**9.0d+0
            t22 = 1/rhoa**6.666666666666666d-1
            t23 = 1.0d+0/t4
            t24 = gammaaa**(t2-1)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = -1.0d+0*t18*t19*t6*(-8.0d+0*t11*t22*t8*t9/t
     1         5**5.0d+0-8.0d+0*gammaaa*t20*t21*t22*t8-1.432881736903528
     2         8d-6*t1*t22*t23*t3*t7)*wght+1.0d+0*t13*t6*(2.666666666666
     3         6666d+0*gammaaa*t21*t22*t8-1.3333333333333333d-6*t1*t22*t
     4         23*t3*t7+2.6666666666666666d+0*gammaaa**2*t15*t16*t17*t22
     5         /t5**1.7d+1-2.6666666666666666d+0*gammaaa*t15*t17*t21*t22
     6         )*wght+1.3333333333333333d+0*t13*t18*t22*t5**3.0d+0*wght+
     7         Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -1.0d+0*t18*t19*t6*(3.0d+0*t10*t11*t8/t9+3
     1         .0d+0*t14*t20*t8+5.373306513388233d-7*t1*t24*t7)*wght+1.0
     2         d+0*t13*t6*(-t14*t8+5.0d-7*t1*t24*t7-gammaaa*t15*t16*t17/
     3         t5**1.6d+1+t14*t15*t17)*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammabb**t2
            t4 = rhob**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammabb**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammabb*t14*t16)
            t18 = -gammabb*t14*t8+1.0d-6*t3*t7+gammabb*t14*t15*t17
            t19 = 1/t12**2
            t20 = 1/(gammabb*t14+1)**5.0d-1
            t21 = 1/t5**9.0d+0
            t22 = 1/rhob**6.666666666666666d-1
            t23 = 1.0d+0/t4
            t24 = gammabb**(t2-1)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = -1.0d+0*t18*t19*t6*(-8.0d+0*t11*t22*t8*t9/t
     1         5**5.0d+0-8.0d+0*gammabb*t20*t21*t22*t8-1.432881736903528
     2         8d-6*t1*t22*t23*t3*t7)*wght+1.0d+0*t13*t6*(2.666666666666
     3         6666d+0*gammabb*t21*t22*t8-1.3333333333333333d-6*t1*t22*t
     4         23*t3*t7+2.6666666666666666d+0*gammabb**2*t15*t16*t17*t22
     5         /t5**1.7d+1-2.6666666666666666d+0*gammabb*t15*t17*t21*t22
     6         )*wght+1.3333333333333333d+0*t13*t18*t22*t5**3.0d+0*wght+
     7         Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -1.0d+0*t18*t19*t6*(3.0d+0*t10*t11*t8/t9+3
     1         .0d+0*t14*t20*t8+5.373306513388233d-7*t1*t24*t7)*wght+1.0
     2         d+0*t13*t6*(-t14*t8+5.0d-7*t1*t24*t7-gammabb*t15*t16*t17/
     3         t5**1.6d+1+t14*t15*t17)*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pw6 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}^{{\it t_6}}}}\\\\ 
C>   {\it t_8} &=& {{{\it t_6}}\over{2}}\\\\ 
C>   {\it t_9} &=& \sigma_{\alpha\alpha}^{{\it t_8}}\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_3}-0.001890381166699926\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{13}} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^{4.0}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{14}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{18}} &=& \sigma_{\beta\beta}^{{\it t_8}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{13}}^{8.0}}}\\\\ 
C>   {\it t_{20}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{21}} &=& {\it t_{20}}^{4.0}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{{\it t_{21}}}}\\\\ 
C>   {\it t_{23}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{24}} &=& {{1}\over{{\it t_{21}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{25}} &=& \sigma_{ss}^{{\it t_8}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{20}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{
C>    \it t_{18}}-{\it t_3}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}
C>    \,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{
C>    \it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{
C>    \it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,
C>    \sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}
C>    \,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,
C>    \sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}
C>    -{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}
C>    \over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{
C>    \it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; 
C>    \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="https://doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pw6_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = 1.0d+0*rhoa**3.333333333333333d-1
            t5 = t4**4.0d+0
            t6 = 1/t5**t1
            t7 = param(1)
            t8 = gammaaa**5.0d-1
            t9 = 1/t5
            t10 = asinh(t8*t9)
            t11 = 6.0d+0*t10*t7*t8*t9+1.0746613026776465d-6*t3*t6+1.0d+0
            t12 = 1/t11
            t13 = 1/t4**8.0d+0
            t14 = t7-1.890381166699926d-3
            t15 = param(2)
            t16 = exp(-gammaaa*t13*t15)
            t17 = -gammaaa*t13*t7+1.0d-6*t3*t6+gammaaa*t13*t14*t16
            t18 = 1/t11**2
            t19 = 1/rhoa
            t20 = (gammaaa*t13+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t4**9.0d+0
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = 1/t4**5.0d+0
            t25 = -8.0d+0*t10*t23*t24*t7*t8-8.0d+0*gammaaa*t21*t22*t23*t
     1         7-1.4328817369035288d-6*t1*t19*t3*t6
            t26 = t4**3.0d+0
            t27 = gammaaa**2
            t28 = 1/t4**1.7d+1
            t29 = 2.6666666666666666d+0*gammaaa*t22*t23*t7-1.33333333333
     1         33333d-6*t1*t19*t3*t6+2.6666666666666666d+0*t14*t15*t16*t
     2         23*t27*t28-2.6666666666666666d+0*gammaaa*t14*t16*t22*t23
            t30 = t2-1
            t31 = gammaaa**t30
            t32 = 1/t4**1.6d+1
            t33 = -t13*t7+5.0d-7*t1*t31*t6-gammaaa*t14*t15*t16*t32+t13*t
     1         14*t16
            t34 = 1/t8
            t35 = 3.0d+0*t10*t34*t7*t9+3.0d+0*t13*t21*t7+5.3733065133882
     1         33d-7*t1*t31*t6
            t36 = 1/t11**3
            t37 = 1/rhoa**1.6666666666666669d+0
            t38 = 1/rhoa**1.3333333333333333d+0
            t39 = 1/rhoa**2
            t40 = t1**2
            t41 = 1/t4**1.0d+1
            t42 = t15**2
            t43 = 1/t4**1.8d+1
            t44 = 1/t20**3
            t45 = gammaaa**(t2-2)
            fnc(iq) = 2.0d+0*t12*t17*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t12*t29*t5-1.0d+0*t17*t18*t25*t5+1.
     1         3333333333333333d+0*t12*t17*t23*t26)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t12*t33*t5-1.0d+0*t17*t18*t35*t5)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t18*(-1.0d+0*t17*t5*(1.333333333333333
     1         3d+1*t10*t38*t7*t8/t4**6.0d+0+5.333333333333333d+0*t10*t2
     2         4*t37*t7*t8-1.0666666666666666d+1*t27*t38*t43*t44*t7+3.46
     3         6666666666666d+1*gammaaa*t21*t38*t41*t7+5.333333333333333
     4         d+0*gammaaa*t21*t22*t37*t7+1.9105089825380384d-6*t3*t39*t
     5         40*t6+1.4328817369035288d-6*t1*t3*t39*t6)-2.0d+0*t25*t29*
     6         t5)+1.0d+0*t12*t5*(-8.0d+0*gammaaa*t38*t41*t7-1.777777777
     7         7777776d+0*gammaaa*t22*t37*t7+1.7777777777777776d-6*t3*t3
     8         9*t40*t6+1.3333333333333333d-6*t1*t3*t39*t6-2.22222222222
     9         2222d+1*t14*t15*t16*t27*t38*t43+7.11111111111111d+0*gamma
     :         aa**3*t14*t16*t38*t42/t4**2.6d+1+8.0d+0*gammaaa*t14*t16*t
     ;         38*t41-1.7777777777777776d+0*t14*t15*t16*t27*t28*t37+1.77
     <         77777777777776d+0*gammaaa*t14*t16*t22*t37)+2.0d+0*t17*t25
     =         **2*t36*t5+1.3333333333333333d+0*t12*t17*t38*t4**2.0d+0-8
     >         .888888888888888d-1*t12*t17*t26*t37+2.6666666666666666d+0
     ?         *t12*t23*t26*t29-2.6666666666666666d+0*t17*t18*t23*t25*t2
     @         6)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t18*(-1.0d+0*t17*t5*(4.0d+0*gammaaa*t
     1         23*t28*t44*t7-4.0d+0*t10*t23*t24*t34*t7-1.2d+1*t21*t22*t2
     2         3*t7-7.164408684517644d-7*t19*t31*t40*t6)-1.0d+0*t29*t35*
     3         t5-1.0d+0*t25*t33*t5)+1.0d+0*t12*t5*(2.6666666666666666d+
     4         0*t22*t23*t7-6.666666666666666d-7*t19*t31*t40*t6-2.666666
     5         6666666666d+0*t14*t16*t23*t27*t42/t4**2.5d+1+8.0d+0*gamma
     6         aa*t14*t15*t16*t23*t28-2.6666666666666666d+0*t14*t16*t22*
     7         t23)+2.0d+0*t17*t25*t35*t36*t5-1.3333333333333333d+0*t17*
     8         t18*t23*t26*t35+1.3333333333333333d+0*t12*t23*t26*t33)*wg
     9         ht+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t18*(-1.0d+0*t17*t5*(-1.5d+0*t10*t7*
     1         t9/t8**3-1.5d+0*t32*t44*t7+1.5d+0*t13*t21*t7/gammaaa+5.37
     2         3306513388233d-7*t1*t30*t45*t6)-2.0d+0*t33*t35*t5)+1.0d+0
     3         *t12*t5*(5.0d-7*t1*t30*t45*t6+gammaaa*t14*t16*t42/t4**2.4
     4         d+1-2*t14*t15*t16*t32)+2.0d+0*t17*t35**2*t36*t5)*wght+Cma
     5         t2(iq,D2_GAA_GAA)
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
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = gammabb**t2
            t20 = rhob**3.333333333333333d-1
            t21 = 1.0d+0*t20
            t22 = t21**4.0d+0
            t23 = 1/t22**t1
            t24 = gammabb**5.0d-1
            t25 = 1/t22
            t26 = asinh(t24*t25)
            t27 = 6.0d+0*t24*t25*t26*t8+1.0746613026776465d-6*t19*t23+1.
     1         0d+0
            t28 = 1/t27
            t29 = 1/t21**8.0d+0
            t30 = exp(-gammabb*t16*t29)
            t31 = -gammabb*t29*t8+gammabb*t15*t29*t30+1.0d-6*t19*t23
            t32 = 1/t12**2
            t33 = (gammaaa*t14+1)**5.0d-1
            t34 = 1/t33
            t35 = 1/t5**9.0d+0
            t36 = 1/rhoa**6.666666666666666d-1
            t37 = -8.0d+0*gammaaa*t34*t35*t36*t8
            t38 = 1/t5**5.0d+0
            t39 = -8.0d+0*t11*t36*t38*t8*t9
            t40 = 1.0d+0/t4
            t41 = -1.4328817369035288d-6*t1*t3*t36*t40*t7+t39+t37
            t42 = t5**3.0d+0
            t43 = 2.6666666666666666d+0*gammaaa*t35*t36*t8
            t44 = gammaaa**2
            t45 = 1/t5**1.7d+1
            t46 = 2.6666666666666666d+0*t15*t16*t17*t36*t44*t45
            t47 = -2.6666666666666666d+0*gammaaa*t15*t17*t35*t36
            t48 = -1.3333333333333333d-6*t1*t3*t36*t40*t7+t47+t46+t43
            t49 = 1/t27**2
            t50 = (gammabb*t29+1)**5.0d-1
            t51 = 1/t50
            t52 = 1/t21**9.0d+0
            t53 = 1/rhob**6.666666666666666d-1
            t54 = -8.0d+0*gammabb*t51*t52*t53*t8
            t55 = 1/t21**5.0d+0
            t56 = -8.0d+0*t24*t26*t53*t55*t8
            t57 = 1.0d+0/t20
            t58 = -1.4328817369035288d-6*t1*t19*t23*t53*t57+t56+t54
            t59 = t21**3.0d+0
            t60 = 2.6666666666666666d+0*gammabb*t52*t53*t8
            t61 = gammabb**2
            t62 = 1/t21**1.7d+1
            t63 = 2.6666666666666666d+0*t15*t16*t30*t53*t61*t62
            t64 = -2.6666666666666666d+0*gammabb*t15*t30*t52*t53
            t65 = t64+t63+t60-1.3333333333333333d-6*t1*t19*t23*t53*t57
            t66 = t2-1
            t67 = gammaaa**t66
            t68 = 1/t5**1.6d+1
            t69 = -t14*t8+5.0d-7*t1*t67*t7-gammaaa*t15*t16*t17*t68+t14*t
     1         15*t17
            t70 = 1/t9
            t71 = 3.0d+0*t10*t11*t70*t8+3.0d+0*t14*t34*t8+5.373306513388
     1         233d-7*t1*t67*t7
            t72 = gammabb**t66
            t73 = 1/t21**1.6d+1
            t74 = -t29*t8-gammabb*t15*t16*t30*t73+5.0d-7*t1*t23*t72+t15*
     1         t29*t30
            t75 = 1/t24
            t76 = 3.0d+0*t25*t26*t75*t8+3.0d+0*t29*t51*t8+5.373306513388
     1         233d-7*t1*t23*t72
            t77 = 1/t12**3
            t78 = 1/rhoa
            t79 = -1.4328817369035288d-6*t1*t3*t7*t78+t39+t37
            t80 = 1/rhoa**1.6666666666666669d+0
            t81 = 1/rhoa**1.3333333333333333d+0
            t82 = 1/rhoa**2
            t83 = t1**2
            t84 = 1/t5**1.0d+1
            t85 = t16**2
            t86 = 1/t5**1.8d+1
            t87 = -1.3333333333333333d-6*t1*t3*t7*t78+t47+t46+t43
            t88 = 1/t33**3
            t89 = 1/t27**3
            t90 = 1/rhob
            t91 = -1.4328817369035288d-6*t1*t19*t23*t90+t56+t54
            t92 = 1/rhob**1.6666666666666669d+0
            t93 = 1/rhob**1.3333333333333333d+0
            t94 = 1/rhob**2
            t95 = 1/t21**1.0d+1
            t96 = 1/t21**1.8d+1
            t97 = -1.3333333333333333d-6*t1*t19*t23*t90+t64+t63+t60
            t98 = 1/t50**3
            t99 = t2-2
            t100 = gammaaa**t99
            t101 = gammabb**t99
            fnc(iq) = (1.0d+0*t13*t18*t6+1.0d+0*t22*t28*t31)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = (1.0d+0*t13*t48*t6-1.0d+0*t18*t32*t41*t6+1.
     1         3333333333333333d+0*t13*t18*t36*t42)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t22*t28*t65+1.3333333333333333d+0*t
     1         28*t31*t53*t59-1.0d+0*t22*t31*t49*t58)*wght+Amat(iq,D1_RB
     2         )
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t6*t69-1.0d+0*t18*t32*t6*t71)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t22*t28*t74-1.0d+0*t22*t31*t49*t76
     1         )*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t32*(-1.0d+0*t18*t6*(1.333333333333333
     1         3d+1*t11*t8*t81*t9/t5**6.0d+0+5.333333333333333d+0*t11*t3
     2         8*t8*t80*t9-1.0666666666666666d+1*t44*t8*t81*t86*t88+3.46
     3         6666666666666d+1*gammaaa*t34*t8*t81*t84+1.910508982538038
     4         4d-6*t3*t40*t7*t80*t83+1.4328817369035288d-6*t1*t3*t7*t82
     5         +5.333333333333333d+0*gammaaa*t34*t35*t8*t80)-1.0d+0*t41*
     6         t6*t87-1.0d+0*t48*t6*t79)+t13*t36*(1.3333333333333333d+0*
     7         t42*t87+1.3333333333333333d+0*t42*t48)+1.0d+0*t13*t6*(-2.
     8         222222222222222d+1*t15*t16*t17*t44*t81*t86+7.111111111111
     9         11d+0*gammaaa**3*t15*t17*t81*t85/t5**2.6d+1-8.0d+0*gammaa
     :         a*t8*t81*t84+8.0d+0*gammaaa*t15*t17*t81*t84+1.77777777777
     ;         77776d-6*t3*t40*t7*t80*t83+1.3333333333333333d-6*t1*t3*t7
     <         *t82-1.7777777777777776d+0*gammaaa*t35*t8*t80-1.777777777
     =         7777776d+0*t15*t16*t17*t44*t45*t80+1.7777777777777776d+0*
     >         gammaaa*t15*t17*t35*t80)+1.3333333333333333d+0*t13*t18*t5
     ?         **2.0d+0*t81-8.888888888888888d-1*t13*t18*t42*t80+t32*t36
     @         *(-1.3333333333333333d+0*t18*t42*t79-1.3333333333333333d+
     1         0*t18*t41*t42)+2.0d+0*t18*t41*t6*t77*t79)*wght+Amat2(iq,D
     2         2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t49*(-1.0d+0*t22*t31*(-1.0666666666666
     1         666d+1*t61*t8*t93*t96*t98+3.466666666666666d+1*gammabb*t5
     2         1*t8*t93*t95+1.4328817369035288d-6*t1*t19*t23*t94+1.33333
     3         33333333333d+1*t24*t26*t8*t93/t21**6.0d+0+1.9105089825380
     4         384d-6*t19*t23*t57*t83*t92+5.333333333333333d+0*t24*t26*t
     5         55*t8*t92+5.333333333333333d+0*gammabb*t51*t52*t8*t92)-1.
     6         0d+0*t22*t58*t97-1.0d+0*t22*t65*t91)+t28*t53*(1.333333333
     7         3333333d+0*t59*t97+1.3333333333333333d+0*t59*t65)+1.0d+0*
     8         t22*t28*(-2.222222222222222d+1*t15*t16*t30*t61*t93*t96-8.
     9         0d+0*gammabb*t8*t93*t95+8.0d+0*gammabb*t15*t30*t93*t95+1.
     :         3333333333333333d-6*t1*t19*t23*t94+7.11111111111111d+0*ga
     ;         mmabb**3*t15*t30*t85*t93/t21**2.6d+1+1.7777777777777776d-
     <         6*t19*t23*t57*t83*t92-1.7777777777777776d+0*gammabb*t52*t
     =         8*t92-1.7777777777777776d+0*t15*t16*t30*t61*t62*t92+1.777
     >         7777777777776d+0*gammabb*t15*t30*t52*t92)+1.3333333333333
     ?         333d+0*t21**2.0d+0*t28*t31*t93-8.888888888888888d-1*t28*t
     @         31*t59*t92+t49*t53*(-1.3333333333333333d+0*t31*t59*t91-1.
     1         3333333333333333d+0*t31*t58*t59)+2.0d+0*t22*t31*t58*t89*t
     2         91)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t32*(-1.0d+0*t18*t6*(4.0d+0*gammaaa*t
     1         36*t45*t8*t88-7.164408684517644d-7*t67*t7*t78*t83-4.0d+0*
     2         t11*t36*t38*t70*t8-1.2d+1*t34*t35*t36*t8)-1.0d+0*t6*t71*t
     3         87-1.0d+0*t6*t69*t79)+1.0d+0*t13*t6*(-2.6666666666666666d
     4         +0*t15*t17*t36*t44*t85/t5**2.5d+1-6.666666666666666d-7*t6
     5         7*t7*t78*t83+2.6666666666666666d+0*t35*t36*t8+8.0d+0*gamm
     6         aaa*t15*t16*t17*t36*t45-2.6666666666666666d+0*t15*t17*t35
     7         *t36)+2.0d+0*t18*t6*t71*t77*t79-1.3333333333333333d+0*t18
     8         *t32*t36*t42*t71+1.3333333333333333d+0*t13*t36*t42*t69)*w
     9         ght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t49*(-1.0d+0*t22*t31*(4.0d+0*gammabb*
     1         t53*t62*t8*t98-7.164408684517644d-7*t23*t72*t83*t90-4.0d+
     2         0*t26*t53*t55*t75*t8-1.2d+1*t51*t52*t53*t8)-1.0d+0*t22*t7
     3         6*t97-1.0d+0*t22*t74*t91)+2.0d+0*t22*t31*t76*t89*t91+1.0d
     4         +0*t22*t28*(-6.666666666666666d-7*t23*t72*t83*t90-2.66666
     5         66666666666d+0*t15*t30*t53*t61*t85/t21**2.5d+1+2.66666666
     6         66666666d+0*t52*t53*t8+8.0d+0*gammabb*t15*t16*t30*t53*t62
     7         -2.6666666666666666d+0*t15*t30*t52*t53)-1.333333333333333
     8         3d+0*t31*t49*t53*t59*t76+1.3333333333333333d+0*t28*t53*t5
     9         9*t74)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t32*(-1.0d+0*t18*t6*(-1.5d+0*t10*t11
     1         *t8/t9**3-1.5d+0*t68*t8*t88+1.5d+0*t14*t34*t8/gammaaa+5.3
     2         73306513388233d-7*t1*t100*t66*t7)-2.0d+0*t6*t69*t71)+1.0d
     3         +0*t13*t6*(gammaaa*t15*t17*t85/t5**2.4d+1+5.0d-7*t1*t100*
     4         t66*t7-2*t15*t16*t17*t68)+2.0d+0*t18*t6*t71**2*t77)*wght+
     5         Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t49*(-1.0d+0*t22*t31*(-1.5d+0*t73*t8
     1         *t98+1.5d+0*t29*t51*t8/gammabb-1.5d+0*t25*t26*t8/t24**3+5
     2         .373306513388233d-7*t1*t101*t23*t66)-2.0d+0*t22*t74*t76)+
     3         2.0d+0*t22*t31*t76**2*t89+1.0d+0*t22*t28*(gammabb*t15*t30
     4         *t85/t21**2.4d+1-2*t15*t16*t30*t73+5.0d-7*t1*t101*t23*t66
     5         ))*wght+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = 1/t12**2
            t20 = (gammaaa*t14+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t5**9.0d+0
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = -8.0d+0*gammaaa*t21*t22*t23*t8
            t25 = 1/t5**5.0d+0
            t26 = -8.0d+0*t11*t23*t25*t8*t9
            t27 = 1.0d+0/t4
            t28 = -1.4328817369035288d-6*t1*t23*t27*t3*t7+t26+t24
            t29 = t5**3.0d+0
            t30 = 2.6666666666666666d+0*gammaaa*t22*t23*t8
            t31 = gammaaa**2
            t32 = 1/t5**1.7d+1
            t33 = 2.6666666666666666d+0*t15*t16*t17*t23*t31*t32
            t34 = -2.6666666666666666d+0*gammaaa*t15*t17*t22*t23
            t35 = -1.3333333333333333d-6*t1*t23*t27*t3*t7+t34+t33+t30
            t36 = t2-1
            t37 = gammaaa**t36
            t38 = 1/t5**1.6d+1
            t39 = -t14*t8+5.0d-7*t1*t37*t7-gammaaa*t15*t16*t17*t38+t14*t
     1         15*t17
            t40 = 1/t9
            t41 = 3.0d+0*t10*t11*t40*t8+3.0d+0*t14*t21*t8+5.373306513388
     1         233d-7*t1*t37*t7
            t42 = 1/t12**3
            t43 = 1/rhoa
            t44 = -1.4328817369035288d-6*t1*t3*t43*t7+t26+t24
            t45 = 1/rhoa**1.6666666666666669d+0
            t46 = 1/rhoa**1.3333333333333333d+0
            t47 = 1/rhoa**2
            t48 = t1**2
            t49 = 1/t5**1.0d+1
            t50 = t16**2
            t51 = 1/t5**1.8d+1
            t52 = -1.3333333333333333d-6*t1*t3*t43*t7+t34+t33+t30
            t53 = 1/t20**3
            t54 = gammaaa**(t2-2)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t13*t35*t6*wght-1.0d+0*t18*t19*t28*t
     1         6*wght+1.3333333333333333d+0*t13*t18*t23*t29*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = -1.0d+0*t18*t19*t41*t6*wght+1.0d+0*t13*t39
     1         *t6*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t19*(-1.0d+0*t18*t6*(1.3333333333333333
     1         d+1*t11*t46*t8*t9/t5**6.0d+0+5.333333333333333d+0*t11*t25
     2         *t45*t8*t9-1.0666666666666666d+1*t31*t46*t51*t53*t8+3.466
     3         666666666666d+1*gammaaa*t21*t46*t49*t8+5.333333333333333d
     4         +0*gammaaa*t21*t22*t45*t8+1.9105089825380384d-6*t27*t3*t4
     5         5*t48*t7+1.4328817369035288d-6*t1*t3*t47*t7)*wght-1.0d+0*
     6         t28*t52*t6*wght-1.0d+0*t35*t44*t6*wght)+t13*t23*(1.333333
     7         3333333333d+0*t29*t52*wght+1.3333333333333333d+0*t29*t35*
     8         wght)+t19*t23*(-1.3333333333333333d+0*t18*t29*t44*wght-1.
     9         3333333333333333d+0*t18*t28*t29*wght)+1.0d+0*t13*t6*(-8.0
     :         d+0*gammaaa*t46*t49*t8-1.7777777777777776d+0*gammaaa*t22*
     ;         t45*t8+1.7777777777777776d-6*t27*t3*t45*t48*t7+1.33333333
     <         33333333d-6*t1*t3*t47*t7-2.222222222222222d+1*t15*t16*t17
     =         *t31*t46*t51+7.11111111111111d+0*gammaaa**3*t15*t17*t46*t
     >         50/t5**2.6d+1+8.0d+0*gammaaa*t15*t17*t46*t49-1.7777777777
     ?         777776d+0*t15*t16*t17*t31*t32*t45+1.7777777777777776d+0*g
     @         ammaaa*t15*t17*t22*t45)*wght+2.0d+0*t18*t28*t42*t44*t6*wg
     1         ht+1.3333333333333333d+0*t13*t18*t46*t5**2.0d+0*wght-8.88
     2         8888888888888d-1*t13*t18*t29*t45*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t19*(-1.0d+0*t18*t6*(4.0d+0*gammaaa*t2
     1         3*t32*t53*t8-4.0d+0*t11*t23*t25*t40*t8-1.2d+1*t21*t22*t23
     2         *t8-7.164408684517644d-7*t37*t43*t48*t7)*wght-1.0d+0*t41*
     3         t52*t6*wght-1.0d+0*t39*t44*t6*wght)+1.0d+0*t13*t6*(2.6666
     4         666666666666d+0*t22*t23*t8-6.666666666666666d-7*t37*t43*t
     5         48*t7-2.6666666666666666d+0*t15*t17*t23*t31*t50/t5**2.5d+
     6         1+8.0d+0*gammaaa*t15*t16*t17*t23*t32-2.6666666666666666d+
     7         0*t15*t17*t22*t23)*wght+2.0d+0*t18*t41*t42*t44*t6*wght-1.
     8         3333333333333333d+0*t18*t19*t23*t29*t41*wght+1.3333333333
     9         333333d+0*t13*t23*t29*t39*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = t19*(-1.0d+0*t18*t6*(-1.5d+0*t10*t11*
     1         t8/t9**3-1.5d+0*t38*t53*t8+1.5d+0*t14*t21*t8/gammaaa+5.37
     2         3306513388233d-7*t1*t36*t54*t7)*wght-2.0d+0*t39*t41*t6*wg
     3         ht)+1.0d+0*t13*t6*(5.0d-7*t1*t36*t54*t7+gammaaa*t15*t17*t
     4         50/t5**2.4d+1-2*t15*t16*t17*t38)*wght+2.0d+0*t18*t41**2*t
     5         42*t6*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammabb**t2
            t4 = rhob**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammabb**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammabb*t14*t16)
            t18 = -gammabb*t14*t8+1.0d-6*t3*t7+gammabb*t14*t15*t17
            t19 = 1/t12**2
            t20 = (gammabb*t14+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t5**9.0d+0
            t23 = 1/rhob**6.666666666666666d-1
            t24 = -8.0d+0*gammabb*t21*t22*t23*t8
            t25 = 1/t5**5.0d+0
            t26 = -8.0d+0*t11*t23*t25*t8*t9
            t27 = 1.0d+0/t4
            t28 = -1.4328817369035288d-6*t1*t23*t27*t3*t7+t26+t24
            t29 = t5**3.0d+0
            t30 = 2.6666666666666666d+0*gammabb*t22*t23*t8
            t31 = gammabb**2
            t32 = 1/t5**1.7d+1
            t33 = 2.6666666666666666d+0*t15*t16*t17*t23*t31*t32
            t34 = -2.6666666666666666d+0*gammabb*t15*t17*t22*t23
            t35 = -1.3333333333333333d-6*t1*t23*t27*t3*t7+t34+t33+t30
            t36 = t2-1
            t37 = gammabb**t36
            t38 = 1/t5**1.6d+1
            t39 = -t14*t8+5.0d-7*t1*t37*t7-gammabb*t15*t16*t17*t38+t14*t
     1         15*t17
            t40 = 1/t9
            t41 = 3.0d+0*t10*t11*t40*t8+3.0d+0*t14*t21*t8+5.373306513388
     1         233d-7*t1*t37*t7
            t42 = 1/t12**3
            t43 = 1/rhob
            t44 = -1.4328817369035288d-6*t1*t3*t43*t7+t26+t24
            t45 = 1/rhob**1.6666666666666669d+0
            t46 = 1/rhob**1.3333333333333333d+0
            t47 = 1/rhob**2
            t48 = t1**2
            t49 = 1/t5**1.0d+1
            t50 = t16**2
            t51 = 1/t5**1.8d+1
            t52 = -1.3333333333333333d-6*t1*t3*t43*t7+t34+t33+t30
            t53 = 1/t20**3
            t54 = gammabb**(t2-2)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t13*t35*t6*wght-1.0d+0*t18*t19*t28*t
     1         6*wght+1.3333333333333333d+0*t13*t18*t23*t29*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GBB) = -1.0d+0*t18*t19*t41*t6*wght+1.0d+0*t13*t39
     1         *t6*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t19*(-1.0d+0*t18*t6*(1.3333333333333333
     1         d+1*t11*t46*t8*t9/t5**6.0d+0+5.333333333333333d+0*t11*t25
     2         *t45*t8*t9-1.0666666666666666d+1*t31*t46*t51*t53*t8+3.466
     3         666666666666d+1*gammabb*t21*t46*t49*t8+5.333333333333333d
     4         +0*gammabb*t21*t22*t45*t8+1.9105089825380384d-6*t27*t3*t4
     5         5*t48*t7+1.4328817369035288d-6*t1*t3*t47*t7)*wght-1.0d+0*
     6         t28*t52*t6*wght-1.0d+0*t35*t44*t6*wght)+t13*t23*(1.333333
     7         3333333333d+0*t29*t52*wght+1.3333333333333333d+0*t29*t35*
     8         wght)+t19*t23*(-1.3333333333333333d+0*t18*t29*t44*wght-1.
     9         3333333333333333d+0*t18*t28*t29*wght)+1.0d+0*t13*t6*(-8.0
     :         d+0*gammabb*t46*t49*t8-1.7777777777777776d+0*gammabb*t22*
     ;         t45*t8+1.7777777777777776d-6*t27*t3*t45*t48*t7+1.33333333
     <         33333333d-6*t1*t3*t47*t7-2.222222222222222d+1*t15*t16*t17
     =         *t31*t46*t51+7.11111111111111d+0*gammabb**3*t15*t17*t46*t
     >         50/t5**2.6d+1+8.0d+0*gammabb*t15*t17*t46*t49-1.7777777777
     ?         777776d+0*t15*t16*t17*t31*t32*t45+1.7777777777777776d+0*g
     @         ammabb*t15*t17*t22*t45)*wght+2.0d+0*t18*t28*t42*t44*t6*wg
     1         ht+1.3333333333333333d+0*t13*t18*t46*t5**2.0d+0*wght-8.88
     2         8888888888888d-1*t13*t18*t29*t45*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t19*(-1.0d+0*t18*t6*(4.0d+0*gammabb*t2
     1         3*t32*t53*t8-4.0d+0*t11*t23*t25*t40*t8-1.2d+1*t21*t22*t23
     2         *t8-7.164408684517644d-7*t37*t43*t48*t7)*wght-1.0d+0*t41*
     3         t52*t6*wght-1.0d+0*t39*t44*t6*wght)+1.0d+0*t13*t6*(2.6666
     4         666666666666d+0*t22*t23*t8-6.666666666666666d-7*t37*t43*t
     5         48*t7-2.6666666666666666d+0*t15*t17*t23*t31*t50/t5**2.5d+
     6         1+8.0d+0*gammabb*t15*t16*t17*t23*t32-2.6666666666666666d+
     7         0*t15*t17*t22*t23)*wght+2.0d+0*t18*t41*t42*t44*t6*wght-1.
     8         3333333333333333d+0*t18*t19*t23*t29*t41*wght+1.3333333333
     9         333333d+0*t13*t23*t29*t39*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = t19*(-1.0d+0*t18*t6*(-1.5d+0*t10*t11*
     1         t8/t9**3-1.5d+0*t38*t53*t8+1.5d+0*t14*t21*t8/gammabb+5.37
     2         3306513388233d-7*t1*t36*t54*t7)*wght-2.0d+0*t39*t41*t6*wg
     3         ht)+1.0d+0*t13*t6*(5.0d-7*t1*t36*t54*t7+gammabb*t15*t17*t
     4         50/t5**2.4d+1-2*t15*t16*t17*t38)*wght+2.0d+0*t18*t41**2*t
     5         42*t6*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pw6 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_4} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_5} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_6} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_7} &=& {{1}\over{{\it t_2}^{{\it t_6}}}}\\\\ 
C>   {\it t_8} &=& {{{\it t_6}}\over{2}}\\\\ 
C>   {\it t_9} &=& \sigma_{\alpha\alpha}^{{\it t_8}}\\\\ 
C>   {\it t_{10}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{11}} &=& {\it t_3}-0.001890381166699926\\\\ 
C>   {\it t_{12}} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_{13}} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_{14}} &=& {\it t_{13}}^{4.0}\\\\ 
C>   {\it t_{15}} &=& {{1}\over{{\it t_{14}}}}\\\\ 
C>   {\it t_{16}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{17}} &=& {{1}\over{{\it t_{14}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{18}} &=& \sigma_{\beta\beta}^{{\it t_8}}\\\\ 
C>   {\it t_{19}} &=& {{1}\over{{\it t_{13}}^{8.0}}}\\\\ 
C>   {\it t_{20}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{21}} &=& {\it t_{20}}^{4.0}\\\\ 
C>   {\it t_{22}} &=& {{1}\over{{\it t_{21}}}}\\\\ 
C>   {\it t_{23}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{24}} &=& {{1}\over{{\it t_{21}}^{{\it t_6}}}}\\\\ 
C>   {\it t_{25}} &=& \sigma_{ss}^{{\it t_8}}\\\\ 
C>   {\it t_{26}} &=& {{1}\over{{\it t_{20}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}
C>    \,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{
C>    \it t_{18}}-{\it t_3}\,{\it t_{19}}\,
C>    \sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}
C>    \,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{
C>    \it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{
C>    \it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{
C>    \it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,
C>    \sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{
C>    -6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}
C>    \,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}
C>    +1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}
C>    \,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,
C>    \sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}
C>    -{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}
C>    \over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{
C>    \it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; 
C>    \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="https://doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pw6_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = 1.0d+0*rhoa**3.333333333333333d-1
            t5 = t4**4.0d+0
            t6 = 1/t5**t1
            t7 = param(1)
            t8 = gammaaa**5.0d-1
            t9 = 1/t5
            t10 = asinh(t8*t9)
            t11 = 6.0d+0*t10*t7*t8*t9+1.0746613026776465d-6*t3*t6+1.0d+0
            t12 = 1/t11
            t13 = 1/t4**8.0d+0
            t14 = t7-1.890381166699926d-3
            t15 = param(2)
            t16 = exp(-gammaaa*t13*t15)
            t17 = -gammaaa*t13*t7+1.0d-6*t3*t6+gammaaa*t13*t14*t16
            t18 = 1/t11**2
            t19 = 1/rhoa
            t20 = (gammaaa*t13+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t4**9.0d+0
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = 1/t4**5.0d+0
            t25 = -8.0d+0*t10*t23*t24*t7*t8-8.0d+0*gammaaa*t21*t22*t23*t
     1         7-1.4328817369035288d-6*t1*t19*t3*t6
            t26 = t4**3.0d+0
            t27 = gammaaa**2
            t28 = 1/t4**1.7d+1
            t29 = 2.6666666666666666d+0*gammaaa*t22*t23*t7-1.33333333333
     1         33333d-6*t1*t19*t3*t6+2.6666666666666666d+0*t14*t15*t16*t
     2         23*t27*t28-2.6666666666666666d+0*gammaaa*t14*t16*t22*t23
            t30 = t2-1
            t31 = gammaaa**t30
            t32 = 1/t4**1.6d+1
            t33 = -t13*t7+5.0d-7*t1*t31*t6-gammaaa*t14*t15*t16*t32+t13*t
     1         14*t16
            t34 = 1/t8
            t35 = 3.0d+0*t10*t34*t7*t9+3.0d+0*t13*t21*t7+5.3733065133882
     1         33d-7*t1*t31*t6
            t36 = 1/t11**3
            t37 = t25**2
            t38 = 1/rhoa**1.6666666666666669d+0
            t39 = t4**2.0d+0
            t40 = 1/rhoa**1.3333333333333333d+0
            t41 = 1/rhoa**2
            t42 = t1**2
            t43 = 1/t4**1.0d+1
            t44 = t15**2
            t45 = gammaaa**3
            t46 = 1/t4**2.6d+1
            t47 = 1/t4**1.8d+1
            t48 = -8.0d+0*gammaaa*t40*t43*t7-1.7777777777777776d+0*gamma
     1         aa*t22*t38*t7+1.7777777777777776d-6*t3*t41*t42*t6+1.33333
     2         33333333333d-6*t1*t3*t41*t6-2.222222222222222d+1*t14*t15*
     3         t16*t27*t40*t47+7.11111111111111d+0*t14*t16*t40*t44*t45*t
     4         46+8.0d+0*gammaaa*t14*t16*t40*t43-1.7777777777777776d+0*t
     5         14*t15*t16*t27*t28*t38+1.7777777777777776d+0*gammaaa*t14*
     6         t16*t22*t38
            t49 = 1/t20**3
            t50 = 1/t4**6.0d+0
            t51 = 1.3333333333333333d+1*t10*t40*t50*t7*t8+5.333333333333
     1         333d+0*t10*t24*t38*t7*t8-1.0666666666666666d+1*t27*t40*t4
     2         7*t49*t7+3.466666666666666d+1*gammaaa*t21*t40*t43*t7+5.33
     3         3333333333333d+0*gammaaa*t21*t22*t38*t7+1.910508982538038
     4         4d-6*t3*t41*t42*t6+1.4328817369035288d-6*t1*t3*t41*t6
            t52 = -1.0d+0*t17*t5*t51-2.0d+0*t25*t29*t5
            t53 = 1/t4**2.5d+1
            t54 = 2.6666666666666666d+0*t22*t23*t7-6.666666666666666d-7*
     1         t19*t31*t42*t6-2.6666666666666666d+0*t14*t16*t23*t27*t44*
     2         t53+8.0d+0*gammaaa*t14*t15*t16*t23*t28-2.6666666666666666
     3         d+0*t14*t16*t22*t23
            t55 = 4.0d+0*gammaaa*t23*t28*t49*t7-4.0d+0*t10*t23*t24*t34*t
     1         7-1.2d+1*t21*t22*t23*t7-7.164408684517644d-7*t19*t31*t42*
     2         t6
            t56 = -1.0d+0*t17*t5*t55-1.0d+0*t29*t35*t5-1.0d+0*t25*t33*t5
            t57 = t2-2
            t58 = gammaaa**t57
            t59 = 1/t4**2.4d+1
            t60 = 5.0d-7*t1*t30*t58*t6+gammaaa*t14*t16*t44*t59-2*t14*t15
     1         *t16*t32
            t61 = t35**2
            t62 = 1/gammaaa
            t63 = 1/t8**3
            t64 = -1.5d+0*t10*t63*t7*t9+1.5d+0*t13*t21*t62*t7-1.5d+0*t32
     1         *t49*t7+5.373306513388233d-7*t1*t30*t58*t6
            t65 = -1.0d+0*t17*t5*t64-2.0d+0*t33*t35*t5
            t66 = 1/t11**4
            t67 = 1/rhoa**2.6666666666666666d+0
            t68 = 1/rhoa**2.3333333333333334d+0
            t69 = 1/rhoa**3
            t70 = t1**3
            t71 = 1/t4**1.1d+1
            t72 = t15**3
            t73 = 1/t4**2.7d+1
            t74 = 1/t4**1.9d+1
            t75 = 1/t20**5
            t76 = gammaaa**(t2-3)
            fnc(iq) = 2.0d+0*t12*t17*t5*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t12*t29*t5-1.0d+0*t17*t18*t25*t5+1.
     1         3333333333333333d+0*t12*t17*t23*t26)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t12*t33*t5-1.0d+0*t17*t18*t35*t5)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (t18*t52+1.0d+0*t12*t48*t5+2.0d+0*t17*t
     1         36*t37*t5+1.3333333333333333d+0*t12*t17*t39*t40-8.8888888
     2         88888888d-1*t12*t17*t26*t38+2.6666666666666666d+0*t12*t23
     3         *t26*t29-2.6666666666666666d+0*t17*t18*t23*t25*t26)*wght+
     4         Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t18*t56+1.0d+0*t12*t5*t54+2.0d+0*t17*
     1         t25*t35*t36*t5-1.3333333333333333d+0*t17*t18*t23*t26*t35+
     2         1.3333333333333333d+0*t12*t23*t26*t33)*wght+Cmat2(iq,D2_R
     3         A_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t18*t65+2.0d+0*t17*t36*t5*t61+1.0d+0
     1         *t12*t5*t60)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (t18*(-1.0d+0*t17*t5*(-2.66666666666
     1         66666d+1*t10*t50*t68*t7*t8-8.88888888888889d+0*t10*t24*t6
     2         7*t7*t8-2.6666666666666666d+1*t10*t41*t7*t8/t4**7.0d+0-4.
     3         2666666666666664d+1*t41*t45*t7*t73*t75+1.1022222222222221
     4         d+2*t27*t41*t49*t7*t74-1.333333333333333d+2*gammaaa*t21*t
     5         41*t7*t71-2.5473453100507176d-6*t3*t6*t69*t70+2.133333333
     6         3333332d+1*t27*t47*t49*t68*t7-6.933333333333332d+1*gammaa
     7         a*t21*t43*t68*t7-8.88888888888889d+0*gammaaa*t21*t22*t67*
     8         t7-5.731526947614115d-6*t3*t42*t6*t69-2.8657634738070575d
     9         -6*t1*t3*t6*t69)-3.0d+0*t29*t5*t51-1.3333333333333333d+0*
     :         t17*t23*t26*t51-3.0d+0*t25*t48*t5-2.6666666666666666d+0*t
     ;         23*t25*t26*t29)+1.0d+0*t12*t5*(1.5466666666666665d+2*t14*
     <         t15*t16*t27*t41*t74-1.2088888888888887d+2*t14*t16*t41*t44
     =         *t45*t73+1.8962962962962962d+1*gammaaa**4*t14*t16*t41*t72
     >         /t4**3.5d+1+2.6666666666666666d+1*gammaaa*t41*t7*t71-2.66
     ?         66666666666666d+1*gammaaa*t14*t16*t41*t71-2.3703703703703
     @         7d-6*t3*t6*t69*t70+1.6d+1*gammaaa*t43*t68*t7+2.9629629629
     1         62963d+0*gammaaa*t22*t67*t7-5.333333333333333d-6*t3*t42*t
     2         6*t69-2.6666666666666666d-6*t1*t3*t6*t69+4.44444444444444
     3         4d+1*t14*t15*t16*t27*t47*t68-1.4222222222222222d+1*t14*t1
     4         6*t44*t45*t46*t68-1.6d+1*gammaaa*t14*t16*t43*t68+2.962962
     5         962962963d+0*t14*t15*t16*t27*t28*t67-2.962962962962963d+0
     6         *gammaaa*t14*t16*t22*t67)-2.6666666666666666d+0*t12*t17*t
     7         39*t68+1.4814814814814814d+0*t12*t17*t26*t67-6.0d+0*t17*t
     8         25**3*t5*t66+t36*(-2*t25*t52+4.0d+0*t17*t25*t5*t51+2.0d+0
     9         *t29*t37*t5)+t18*t23*(-2.6666666666666666d+0*t17*t26*t51-
     :         2.6666666666666666d+0*t17*t23*t25*t39-5.333333333333333d+
     ;         0*t25*t26*t29)+t12*t23*(4.0d+0*t26*t48+2.6666666666666666
     <         d+0*t23*t29*t39)+1.3333333333333333d+0*t12*t29*t39*t40-1.
     =         3333333333333333d+0*t17*t18*t25*t39*t40-2.666666666666666
     >         6d+0*t12*t26*t29*t38+2.6666666666666666d+0*t17*t18*t25*t2
     ?         6*t38+8.888888888888888d-1*t12*t17*t38+8.0d+0*t17*t23*t26
     @         *t36*t37)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t18*(-1.0d+0*t17*t5*(1.6d+1*t27*t4
     1         0*t46*t7*t75+9.552544912690191d-7*t31*t41*t6*t70+6.666666
     2         666666666d+0*t10*t34*t40*t50*t7-3.8666666666666666d+1*gam
     3         maaa*t40*t47*t49*t7-2.6666666666666666d+0*gammaaa*t28*t38
     4         *t49*t7+4.133333333333333d+1*t21*t40*t43*t7+2.66666666666
     5         66666d+0*t10*t24*t34*t38*t7+8.0d+0*t21*t22*t38*t7+7.16440
     6         8684517644d-7*t31*t41*t42*t6)-2.0d+0*t29*t5*t55-2.0d+0*t2
     7         5*t5*t54-1.0d+0*t33*t5*t51-1.0d+0*t35*t48*t5)+1.0d+0*t12*
     8         t5*(-7.11111111111111d+0*t14*t16*t40*t45*t72/t4**3.4d+1+8
     9         .888888888888887d-7*t31*t41*t6*t70-8.0d+0*t40*t43*t7-1.77
     :         77777777777776d+0*t22*t38*t7+6.666666666666666d-7*t31*t41
     ;         *t42*t6+1.7777777777777776d+0*t14*t16*t27*t38*t44*t53-5.2
     <         44444444444444d+1*gammaaa*t14*t15*t16*t40*t47+4.355555555
     =         555556d+1*t14*t16*t27*t40*t44*t46+8.0d+0*t14*t16*t40*t43-
     >         5.333333333333333d+0*gammaaa*t14*t15*t16*t28*t38+1.777777
     ?         7777777776d+0*t14*t16*t22*t38)-6.0d+0*t17*t35*t37*t5*t66+
     @         t36*(4.0d+0*t17*t25*t5*t55-2*t35*t52+2.0d+0*t33*t37*t5)+t
     1         18*t23*(-2.6666666666666666d+0*t17*t26*t55-2.666666666666
     2         6666d+0*t26*t29*t35-2.6666666666666666d+0*t25*t26*t33)+2.
     3         6666666666666666d+0*t12*t23*t26*t54-1.3333333333333333d+0
     4         *t17*t18*t35*t39*t40+1.3333333333333333d+0*t12*t33*t39*t4
     5         0+8.888888888888888d-1*t17*t18*t26*t35*t38-8.888888888888
     6         888d-1*t12*t26*t33*t38+5.333333333333333d+0*t17*t23*t25*t
     7         26*t35*t36)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t18*(-1.0d+0*t17*t5*(-6.0d+0*gamm
     1         aaa*t23*t53*t7*t75+2.0d+0*t10*t23*t24*t63*t7-2.0d+0*t21*t
     2         22*t23*t62*t7+1.0d+1*t23*t28*t49*t7-7.164408684517644d-7*
     3         t19*t30*t42*t58*t6)-1.0d+0*t29*t5*t64-1.0d+0*t25*t5*t60-2
     4         .0d+0*t33*t5*t55-2.0d+0*t35*t5*t54)+1.0d+0*t12*t5*(2.6666
     5         666666666666d+0*t14*t16*t23*t27*t72/t4**3.3d+1-6.66666666
     6         6666666d-7*t19*t30*t42*t58*t6-1.3333333333333333d+1*gamma
     7         aa*t14*t16*t23*t44*t53+1.0666666666666666d+1*t14*t15*t16*
     8         t23*t28)-6.0d+0*t17*t25*t5*t61*t66+t36*(2.0d+0*t17*t25*t5
     9         *t64-2*t35*t56+2.0d+0*t17*t35*t5*t55+2.0d+0*t25*t33*t35*t
     :         5)+t18*t23*(-1.3333333333333333d+0*t17*t26*t64-2.66666666
     ;         66666666d+0*t26*t33*t35)+2.6666666666666666d+0*t17*t23*t2
     <         6*t36*t61+1.3333333333333333d+0*t12*t23*t26*t60)*wght+Cma
     =         t3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t18*(-1.0d+0*t17*t5*(2.25d+0*t10
     1         *t7*t9/t8**5+5.373306513388233d-7*t1*t30*t57*t6*t76+2.25d
     2         +0*t59*t7*t75-7.5d-1*t32*t49*t62*t7-2.25d+0*t13*t21*t7/t2
     3         7)-3.0d+0*t33*t5*t64-3.0d+0*t35*t5*t60)+1.0d+0*t12*t5*(5.
     4         0d-7*t1*t30*t57*t6*t76-gammaaa*t14*t16*t72/t4**3.2d+1+3*t
     5         14*t16*t44*t59)-6.0d+0*t17*t35**3*t5*t66+t36*(-2*t35*t65+
     6         4.0d+0*t17*t35*t5*t64+2.0d+0*t33*t5*t61))*wght+Cmat3(iq,D
     7         3_GAA_GAA_GAA)
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
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = gammabb**t2
            t20 = rhob**3.333333333333333d-1
            t21 = 1.0d+0*t20
            t22 = t21**4.0d+0
            t23 = 1/t22**t1
            t24 = gammabb**5.0d-1
            t25 = 1/t22
            t26 = asinh(t24*t25)
            t27 = 6.0d+0*t24*t25*t26*t8+1.0746613026776465d-6*t19*t23+1.
     1         0d+0
            t28 = 1/t27
            t29 = 1/t21**8.0d+0
            t30 = exp(-gammabb*t16*t29)
            t31 = -gammabb*t29*t8+gammabb*t15*t29*t30+1.0d-6*t19*t23
            t32 = 1/t12**2
            t33 = (gammaaa*t14+1)**5.0d-1
            t34 = 1/t33
            t35 = 1/t5**9.0d+0
            t36 = 1/rhoa**6.666666666666666d-1
            t37 = -8.0d+0*gammaaa*t34*t35*t36*t8
            t38 = 1/t5**5.0d+0
            t39 = -8.0d+0*t11*t36*t38*t8*t9
            t40 = 1.0d+0/t4
            t41 = -1.4328817369035288d-6*t1*t3*t36*t40*t7+t39+t37
            t42 = t5**3.0d+0
            t43 = 2.6666666666666666d+0*gammaaa*t35*t36*t8
            t44 = gammaaa**2
            t45 = 1/t5**1.7d+1
            t46 = 2.6666666666666666d+0*t15*t16*t17*t36*t44*t45
            t47 = -2.6666666666666666d+0*gammaaa*t15*t17*t35*t36
            t48 = -1.3333333333333333d-6*t1*t3*t36*t40*t7+t47+t46+t43
            t49 = 1/t27**2
            t50 = (gammabb*t29+1)**5.0d-1
            t51 = 1/t50
            t52 = 1/t21**9.0d+0
            t53 = 1/rhob**6.666666666666666d-1
            t54 = -8.0d+0*gammabb*t51*t52*t53*t8
            t55 = 1/t21**5.0d+0
            t56 = -8.0d+0*t24*t26*t53*t55*t8
            t57 = 1.0d+0/t20
            t58 = -1.4328817369035288d-6*t1*t19*t23*t53*t57+t56+t54
            t59 = t21**3.0d+0
            t60 = 2.6666666666666666d+0*gammabb*t52*t53*t8
            t61 = gammabb**2
            t62 = 1/t21**1.7d+1
            t63 = 2.6666666666666666d+0*t15*t16*t30*t53*t61*t62
            t64 = -2.6666666666666666d+0*gammabb*t15*t30*t52*t53
            t65 = t64+t63+t60-1.3333333333333333d-6*t1*t19*t23*t53*t57
            t66 = t2-1
            t67 = gammaaa**t66
            t68 = 1/t5**1.6d+1
            t69 = -t14*t8+5.0d-7*t1*t67*t7-gammaaa*t15*t16*t17*t68+t14*t
     1         15*t17
            t70 = 1/t9
            t71 = 3.0d+0*t10*t11*t70*t8+3.0d+0*t14*t34*t8+5.373306513388
     1         233d-7*t1*t67*t7
            t72 = gammabb**t66
            t73 = 1/t21**1.6d+1
            t74 = -t29*t8-gammabb*t15*t16*t30*t73+5.0d-7*t1*t23*t72+t15*
     1         t29*t30
            t75 = 1/t24
            t76 = 3.0d+0*t25*t26*t75*t8+3.0d+0*t29*t51*t8+5.373306513388
     1         233d-7*t1*t23*t72
            t77 = 1/t12**3
            t78 = 1/rhoa
            t79 = -1.4328817369035288d-6*t1*t3*t7*t78+t39+t37
            t80 = 1/rhoa**1.6666666666666669d+0
            t81 = t5**2.0d+0
            t82 = 1/rhoa**1.3333333333333333d+0
            t83 = 1/rhoa**2
            t84 = 1.3333333333333333d-6*t1*t3*t7*t83
            t85 = -1.7777777777777776d+0*gammaaa*t35*t8*t80
            t86 = t1**2
            t87 = 1/t5**1.0d+1
            t88 = -8.0d+0*gammaaa*t8*t82*t87
            t89 = -1.7777777777777776d+0*t15*t16*t17*t44*t45*t80
            t90 = 1.7777777777777776d+0*gammaaa*t15*t17*t35*t80
            t91 = t16**2
            t92 = gammaaa**3
            t93 = 1/t5**2.6d+1
            t94 = 7.11111111111111d+0*t15*t17*t82*t91*t92*t93
            t95 = 1/t5**1.8d+1
            t96 = -2.222222222222222d+1*t15*t16*t17*t44*t82*t95
            t97 = 8.0d+0*gammaaa*t15*t17*t82*t87
            t98 = t97+t96+t94+t90+t89+t88+1.7777777777777776d-6*t3*t40*t
     1         7*t80*t86+t85+t84
            t99 = -1.3333333333333333d-6*t1*t3*t7*t78+t47+t46+t43
            t100 = 1.4328817369035288d-6*t1*t3*t7*t83
            t101 = 5.333333333333333d+0*gammaaa*t34*t35*t8*t80
            t102 = 5.333333333333333d+0*t11*t38*t8*t80*t9
            t103 = 1/t33**3
            t104 = -1.0666666666666666d+1*t103*t44*t8*t82*t95
            t105 = 3.466666666666666d+1*gammaaa*t34*t8*t82*t87
            t106 = 1/t5**6.0d+0
            t107 = 1.3333333333333333d+1*t106*t11*t8*t82*t9
            t108 = 1.9105089825380384d-6*t3*t40*t7*t80*t86+t107+t105+t10
     1         4+t102+t101+t100
            t109 = 1/t27**3
            t110 = 1/rhob
            t111 = t56+t54-1.4328817369035288d-6*t1*t110*t19*t23
            t112 = 1/rhob**1.6666666666666669d+0
            t113 = t21**2.0d+0
            t114 = 1/rhob**1.3333333333333333d+0
            t115 = 1/rhob**2
            t116 = 1.3333333333333333d-6*t1*t115*t19*t23
            t117 = -1.7777777777777776d+0*gammabb*t112*t52*t8
            t118 = 1/t21**1.0d+1
            t119 = -8.0d+0*gammabb*t114*t118*t8
            t120 = -1.7777777777777776d+0*t112*t15*t16*t30*t61*t62
            t121 = 1.7777777777777776d+0*gammabb*t112*t15*t30*t52
            t122 = gammabb**3
            t123 = 1/t21**2.6d+1
            t124 = 7.11111111111111d+0*t114*t122*t123*t15*t30*t91
            t125 = 1/t21**1.8d+1
            t126 = -2.222222222222222d+1*t114*t125*t15*t16*t30*t61
            t127 = 8.0d+0*gammabb*t114*t118*t15*t30
            t128 = 1.7777777777777776d-6*t112*t19*t23*t57*t86+t127+t126+
     1         t124+t121+t120+t119+t117+t116
            t129 = t64+t63+t60-1.3333333333333333d-6*t1*t110*t19*t23
            t130 = 1.4328817369035288d-6*t1*t115*t19*t23
            t131 = 5.333333333333333d+0*gammabb*t112*t51*t52*t8
            t132 = 5.333333333333333d+0*t112*t24*t26*t55*t8
            t133 = 1/t50**3
            t134 = -1.0666666666666666d+1*t114*t125*t133*t61*t8
            t135 = 3.466666666666666d+1*gammabb*t114*t118*t51*t8
            t136 = 1/t21**6.0d+0
            t137 = 1.3333333333333333d+1*t114*t136*t24*t26*t8
            t138 = 1.9105089825380384d-6*t112*t19*t23*t57*t86+t137+t135+
     1         t134+t132+t131+t130
            t139 = 1/t5**2.5d+1
            t140 = -2.6666666666666666d+0*t139*t15*t17*t36*t44*t91-6.666
     1         666666666666d-7*t67*t7*t78*t86+2.6666666666666666d+0*t35*
     2         t36*t8+8.0d+0*gammaaa*t15*t16*t17*t36*t45-2.6666666666666
     3         666d+0*t15*t17*t35*t36
            t141 = -7.164408684517644d-7*t67*t7*t78*t86-4.0d+0*t11*t36*t
     1         38*t70*t8+4.0d+0*gammaaa*t103*t36*t45*t8-1.2d+1*t34*t35*t
     2         36*t8
            t142 = -1.0d+0*t6*t71*t99-1.0d+0*t6*t69*t79-1.0d+0*t141*t18*
     1         t6
            t143 = 1/t21**2.5d+1
            t144 = -2.6666666666666666d+0*t143*t15*t30*t53*t61*t91-6.666
     1         666666666666d-7*t110*t23*t72*t86+2.6666666666666666d+0*t5
     2         2*t53*t8+8.0d+0*gammabb*t15*t16*t30*t53*t62-2.66666666666
     3         66666d+0*t15*t30*t52*t53
            t145 = -7.164408684517644d-7*t110*t23*t72*t86-4.0d+0*t26*t53
     1         *t55*t75*t8+4.0d+0*gammabb*t133*t53*t62*t8-1.2d+1*t51*t52
     2         *t53*t8
            t146 = -1.0d+0*t129*t22*t76-1.0d+0*t111*t22*t74-1.0d+0*t145*
     1         t22*t31
            t147 = t2-2
            t148 = gammaaa**t147
            t149 = 1/t5**2.4d+1
            t150 = gammaaa*t149*t15*t17*t91+5.0d-7*t1*t148*t66*t7-2*t15*
     1         t16*t17*t68
            t151 = t71**2
            t152 = 1/gammaaa
            t153 = 1/t9**3
            t154 = -1.5d+0*t103*t68*t8+1.5d+0*t14*t152*t34*t8-1.5d+0*t10
     1         *t11*t153*t8+5.373306513388233d-7*t1*t148*t66*t7
            t155 = -2.0d+0*t6*t69*t71-1.0d+0*t154*t18*t6
            t156 = gammabb**t147
            t157 = 1/t21**2.4d+1
            t158 = gammabb*t15*t157*t30*t91-2*t15*t16*t30*t73+5.0d-7*t1*
     1         t156*t23*t66
            t159 = t76**2
            t160 = 1/gammabb
            t161 = 1/t24**3
            t162 = -1.5d+0*t133*t73*t8+1.5d+0*t160*t29*t51*t8-1.5d+0*t16
     1         1*t25*t26*t8+5.373306513388233d-7*t1*t156*t23*t66
            t163 = -2.0d+0*t22*t74*t76-1.0d+0*t162*t22*t31
            t164 = 1/t12**4
            t165 = t79**2
            t166 = 1/rhoa**2.6666666666666666d+0
            t167 = 1/rhoa**2.3333333333333334d+0
            t168 = 1/rhoa**3
            t169 = t1**3
            t170 = 1/t5**1.1d+1
            t171 = t16**3
            t172 = 1/t5**2.7d+1
            t173 = 1/t5**1.9d+1
            t174 = t97+t96+t94+t90+t89+t88+1.7777777777777776d-6*t3*t7*t
     1         83*t86+t85+t84
            t175 = 1/t33**5
            t176 = 1.9105089825380384d-6*t3*t7*t83*t86+t107+t105+t104+t1
     1         02+t101+t100
            t177 = -2.0d+0*t6*t79*t99-1.0d+0*t176*t18*t6
            t178 = 1/t27**4
            t179 = t111**2
            t180 = 1/rhob**2.6666666666666666d+0
            t181 = 1/rhob**2.3333333333333334d+0
            t182 = 1/rhob**3
            t183 = 1/t21**1.1d+1
            t184 = 1/t21**2.7d+1
            t185 = 1/t21**1.9d+1
            t186 = 1.7777777777777776d-6*t115*t19*t23*t86+t127+t126+t124
     1         +t121+t120+t119+t117+t116
            t187 = 1/t50**5
            t188 = 1.9105089825380384d-6*t115*t19*t23*t86+t137+t135+t134
     1         +t132+t131+t130
            t189 = -1.0d+0*t188*t22*t31-2.0d+0*t111*t129*t22
            t190 = t2-3
            t191 = gammaaa**t190
            t192 = gammabb**t190
            fnc(iq) = (1.0d+0*t13*t18*t6+1.0d+0*t22*t28*t31)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = (1.0d+0*t13*t48*t6-1.0d+0*t18*t32*t41*t6+1.
     1         3333333333333333d+0*t13*t18*t36*t42)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t22*t28*t65+1.3333333333333333d+0*t
     1         28*t31*t53*t59-1.0d+0*t22*t31*t49*t58)*wght+Amat(iq,D1_RB
     2         )
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t6*t69-1.0d+0*t18*t32*t6*t71)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t22*t28*t74-1.0d+0*t22*t31*t49*t76
     1         )*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t32*(-1.0d+0*t41*t6*t99-1.0d+0*t48*t6*
     1         t79-1.0d+0*t108*t18*t6)+t13*t36*(1.3333333333333333d+0*t4
     2         2*t99+1.3333333333333333d+0*t42*t48)+1.0d+0*t13*t6*t98+1.
     3         3333333333333333d+0*t13*t18*t81*t82-8.888888888888888d-1*
     4         t13*t18*t42*t80+t32*t36*(-1.3333333333333333d+0*t18*t42*t
     5         79-1.3333333333333333d+0*t18*t41*t42)+2.0d+0*t18*t41*t6*t
     6         77*t79)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t28*t53*(1.3333333333333333d+0*t59*t65
     1         +1.3333333333333333d+0*t129*t59)+t49*(-1.0d+0*t111*t22*t6
     2         5-1.0d+0*t129*t22*t58-1.0d+0*t138*t22*t31)+t49*t53*(-1.33
     3         33333333333333d+0*t31*t58*t59-1.3333333333333333d+0*t111*
     4         t31*t59)-8.888888888888888d-1*t112*t28*t31*t59+2.0d+0*t10
     5         9*t111*t22*t31*t58+1.3333333333333333d+0*t113*t114*t28*t3
     6         1+1.0d+0*t128*t22*t28)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0d+0*t18*t6*t71*t77*t79-1.333333333
     1         3333333d+0*t18*t32*t36*t42*t71+1.3333333333333333d+0*t13*
     2         t36*t42*t69+1.0d+0*t13*t140*t6+t142*t32)*wght+Cmat2(iq,D2
     3         _RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.3333333333333333d+0*t31*t49*t53*t5
     1         9*t76+2.0d+0*t109*t111*t22*t31*t76+1.3333333333333333d+0*
     2         t28*t53*t59*t74+t146*t49+1.0d+0*t144*t22*t28)*wght+Cmat2(
     3         iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2.0d+0*t151*t18*t6*t77+1.0d+0*t13*t1
     1         50*t6+t155*t32)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t163*t49+2.0d+0*t109*t159*t22*t31+1.
     1         0d+0*t158*t22*t28)*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (t13*t36*(2.6666666666666666d+0*t36*
     1         t81*t99+2.6666666666666666d+0*t42*t98+1.3333333333333333d
     2         +0*t174*t42)+t32*(-2.6666666666666666d+0*t36*t42*t79*t99-
     3         2.0d+0*t108*t6*t99-2.0d+0*t6*t79*t98-1.0d+0*t18*t6*(2.133
     4         3333333333332d+1*t103*t167*t44*t8*t95-4.2666666666666664d
     5         +1*t172*t175*t8*t83*t92-2.6666666666666666d+1*t11*t8*t83*
     6         t9/t5**7.0d+0-8.88888888888889d+0*t11*t166*t38*t8*t9-2.66
     7         66666666666666d+1*t106*t11*t167*t8*t9-6.933333333333332d+
     8         1*gammaaa*t167*t34*t8*t87-1.9105089825380384d-6*t166*t3*t
     9         40*t7*t86-3.8210179650760767d-6*t168*t3*t7*t86+1.10222222
     :         22222221d+2*t103*t173*t44*t8*t83-1.333333333333333d+2*gam
     ;         maaa*t170*t34*t8*t83-8.88888888888889d+0*gammaaa*t166*t34
     <         *t35*t8-2.5473453100507176d-6*t166*t169*t3*t40*t7-2.86576
     =         34738070575d-6*t1*t168*t3*t7)-1.0d+0*t176*t48*t6-1.0d+0*t
     >         174*t41*t6-1.3333333333333333d+0*t176*t18*t36*t42)+t32*t3
     ?         6*(-2.6666666666666666d+0*t41*t42*t99-2.6666666666666666d
     @         +0*t18*t36*t79*t81-2.6666666666666666d+0*t42*t48*t79-2.66
     1         66666666666666d+0*t108*t18*t42)-1.7777777777777776d+0*t13
     2         *t42*t80*t99+1.0d+0*t13*t6*(4.444444444444444d+1*t15*t16*
     3         t167*t17*t44*t95-1.4222222222222222d+1*t15*t167*t17*t91*t
     4         92*t93-1.2088888888888887d+2*t15*t17*t172*t83*t91*t92+1.6
     5         d+1*gammaaa*t167*t8*t87-1.6d+1*gammaaa*t15*t167*t17*t87-1
     6         .7777777777777776d-6*t166*t3*t40*t7*t86-3.555555555555555
     7         d-6*t168*t3*t7*t86+2.6666666666666666d+1*gammaaa*t170*t8*
     8         t83+1.8962962962962962d+1*gammaaa**4*t15*t17*t171*t83/t5*
     9         *3.5d+1+1.5466666666666665d+2*t15*t16*t17*t173*t44*t83-2.
     :         6666666666666666d+1*gammaaa*t15*t17*t170*t83+2.9629629629
     ;         62963d+0*gammaaa*t166*t35*t8-2.37037037037037d-6*t166*t16
     <         9*t3*t40*t7-2.6666666666666666d-6*t1*t168*t3*t7+2.9629629
     =         62962963d+0*t15*t16*t166*t17*t44*t45-2.962962962962963d+0
     >         *gammaaa*t15*t166*t17*t35)+1.3333333333333333d+0*t13*t48*
     ?         t81*t82-1.3333333333333333d+0*t18*t32*t41*t81*t82-2.66666
     @         66666666666d+0*t13*t167*t18*t81+1.7777777777777776d+0*t18
     1         *t32*t42*t79*t80-8.888888888888888d-1*t13*t42*t48*t80+8.8
     2         88888888888888d-1*t18*t32*t41*t42*t80+8.888888888888888d-
     3         1*t13*t18*t80+t77*(4.0d+0*t108*t18*t6*t79+2.0d+0*t165*t48
     4         *t6-2*t177*t41)+t36*t77*(5.333333333333333d+0*t18*t41*t42
     5         *t79+2.6666666666666666d+0*t165*t18*t42)-6.0d+0*t164*t165
     6         *t18*t41*t6+1.4814814814814814d+0*t13*t166*t18*t42)*wght+
     7         Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (1.0d+0*t22*t28*(-1.2088888888888887
     1         d+2*t115*t122*t15*t184*t30*t91-1.4222222222222222d+1*t122
     2         *t123*t15*t181*t30*t91-1.7777777777777776d-6*t180*t19*t23
     3         *t57*t86-3.555555555555555d-6*t182*t19*t23*t86+2.96296296
     4         2962963d+0*gammabb*t180*t52*t8+2.6666666666666666d+1*gamm
     5         abb*t115*t183*t8+1.6d+1*gammabb*t118*t181*t8+2.9629629629
     6         62963d+0*t15*t16*t180*t30*t61*t62+1.5466666666666665d+2*t
     7         115*t15*t16*t185*t30*t61+4.444444444444444d+1*t125*t15*t1
     8         6*t181*t30*t61-2.37037037037037d-6*t169*t180*t19*t23*t57-
     9         2.962962962962963d+0*gammabb*t15*t180*t30*t52+1.896296296
     :         2962962d+1*gammabb**4*t115*t15*t171*t30/t21**3.5d+1-2.666
     ;         6666666666666d+1*gammabb*t115*t15*t183*t30-1.6d+1*gammabb
     <         *t118*t15*t181*t30-2.6666666666666666d-6*t1*t182*t19*t23)
     =         +t49*(-1.0d+0*t22*t31*(-1.9105089825380384d-6*t180*t19*t2
     >         3*t57*t86-3.8210179650760767d-6*t182*t19*t23*t86+1.102222
     ?         2222222221d+2*t115*t133*t185*t61*t8+2.1333333333333332d+1
     @         *t125*t133*t181*t61*t8-8.88888888888889d+0*t180*t24*t26*t
     1         55*t8-8.88888888888889d+0*gammabb*t180*t51*t52*t8-1.33333
     2         3333333333d+2*gammabb*t115*t183*t51*t8-6.933333333333332d
     3         +1*gammabb*t118*t181*t51*t8-2.6666666666666666d+1*t115*t2
     4         4*t26*t8/t21**7.0d+0-2.6666666666666666d+1*t136*t181*t24*
     5         t26*t8-4.2666666666666664d+1*t115*t122*t184*t187*t8-2.547
     6         3453100507176d-6*t169*t180*t19*t23*t57-2.8657634738070575
     7         d-6*t1*t182*t19*t23)-1.0d+0*t188*t22*t65-1.33333333333333
     8         33d+0*t188*t31*t53*t59-2.6666666666666666d+0*t111*t129*t5
     9         3*t59-1.0d+0*t186*t22*t58-2.0d+0*t129*t138*t22-2.0d+0*t11
     :         1*t128*t22)+t49*t53*(-2.6666666666666666d+0*t111*t59*t65-
     ;         2.6666666666666666d+0*t129*t58*t59-2.6666666666666666d+0*
     <         t138*t31*t59-2.6666666666666666d+0*t111*t113*t31*t53)+t10
     =         9*(2.0d+0*t179*t22*t65-2*t189*t58+4.0d+0*t111*t138*t22*t3
     >         1)-8.888888888888888d-1*t112*t28*t59*t65+1.33333333333333
     ?         33d+0*t113*t114*t28*t65+t109*t53*(5.333333333333333d+0*t1
     @         11*t31*t58*t59+2.6666666666666666d+0*t179*t31*t59)+t28*t5
     1         3*(1.3333333333333333d+0*t186*t59+2.6666666666666666d+0*t
     2         128*t59+2.6666666666666666d+0*t113*t129*t53)+8.8888888888
     3         88888d-1*t112*t31*t49*t58*t59+1.7777777777777776d+0*t111*
     4         t112*t31*t49*t59+1.4814814814814814d+0*t180*t28*t31*t59-1
     5         .7777777777777776d+0*t112*t129*t28*t59-1.3333333333333333
     6         d+0*t113*t114*t31*t49*t58-6.0d+0*t178*t179*t22*t31*t58-2.
     7         6666666666666666d+0*t113*t181*t28*t31+8.888888888888888d-
     8         1*t112*t28*t31)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t32*t36*(-2.6666666666666666d+0*t4
     1         2*t71*t99-2.6666666666666666d+0*t42*t69*t79-2.66666666666
     2         66666d+0*t141*t18*t42)+t32*(-2.0d+0*t141*t6*t99-1.0d+0*t1
     3         8*t6*(-3.8666666666666666d+1*gammaaa*t103*t8*t82*t95+1.6d
     4         +1*t175*t44*t8*t82*t93+4.133333333333333d+1*t34*t8*t82*t8
     5         7+7.164408684517644d-7*t67*t7*t83*t86+9.552544912690191d-
     6         7*t169*t67*t7*t83+6.666666666666666d+0*t106*t11*t70*t8*t8
     7         2+2.6666666666666666d+0*t11*t38*t70*t8*t80-2.666666666666
     8         6666d+0*gammaaa*t103*t45*t8*t80+8.0d+0*t34*t35*t8*t80)-2.
     9         0d+0*t140*t6*t79-1.0d+0*t174*t6*t71-1.0d+0*t176*t6*t69)+1
     :         .0d+0*t13*t6*(-5.244444444444444d+1*gammaaa*t15*t16*t17*t
     ;         82*t95+4.355555555555556d+1*t15*t17*t44*t82*t91*t93-7.111
     <         11111111111d+0*t15*t17*t171*t82*t92/t5**3.4d+1+1.77777777
     =         77777776d+0*t139*t15*t17*t44*t80*t91-8.0d+0*t8*t82*t87+8.
     >         0d+0*t15*t17*t82*t87+6.666666666666666d-7*t67*t7*t83*t86+
     ?         8.888888888888887d-7*t169*t67*t7*t83-1.7777777777777776d+
     @         0*t35*t8*t80-5.333333333333333d+0*gammaaa*t15*t16*t17*t45
     1         *t80+1.7777777777777776d+0*t15*t17*t35*t80)-1.33333333333
     2         33333d+0*t18*t32*t71*t81*t82+1.3333333333333333d+0*t13*t6
     3         9*t81*t82+8.888888888888888d-1*t18*t32*t42*t71*t80-8.8888
     4         88888888888d-1*t13*t42*t69*t80+t77*(4.0d+0*t141*t18*t6*t7
     5         9-2*t177*t71+2.0d+0*t165*t6*t69)+5.333333333333333d+0*t18
     6         *t36*t42*t71*t77*t79-6.0d+0*t164*t165*t18*t6*t71+2.666666
     7         6666666666d+0*t13*t140*t36*t42)*wght+Cmat3(iq,D3_RA_RA_GA
     8         A)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (1.0d+0*t22*t28*(1.7777777777777776
     1         d+0*t112*t143*t15*t30*t61*t91+4.355555555555556d+1*t114*t
     2         123*t15*t30*t61*t91+6.666666666666666d-7*t115*t23*t72*t86
     3         -1.7777777777777776d+0*t112*t52*t8-8.0d+0*t114*t118*t8+8.
     4         888888888888887d-7*t115*t169*t23*t72-5.333333333333333d+0
     5         *gammabb*t112*t15*t16*t30*t62+1.7777777777777776d+0*t112*
     6         t15*t30*t52-7.11111111111111d+0*t114*t122*t15*t171*t30/t2
     7         1**3.4d+1-5.244444444444444d+1*gammabb*t114*t125*t15*t16*
     8         t30+8.0d+0*t114*t118*t15*t30)+t49*(-1.0d+0*t22*t31*(7.164
     9         408684517644d-7*t115*t23*t72*t86+2.6666666666666666d+0*t1
     :         12*t26*t55*t75*t8+6.666666666666666d+0*t114*t136*t26*t75*
     ;         t8-2.6666666666666666d+0*gammabb*t112*t133*t62*t8+1.6d+1*
     <         t114*t123*t187*t61*t8+8.0d+0*t112*t51*t52*t8+4.1333333333
     =         33333d+1*t114*t118*t51*t8-3.8666666666666666d+1*gammabb*t
     >         114*t125*t133*t8+9.552544912690191d-7*t115*t169*t23*t72)-
     ?         1.0d+0*t186*t22*t76-1.0d+0*t188*t22*t74-2.0d+0*t129*t145*
     @         t22-2.0d+0*t111*t144*t22)+t49*t53*(-2.6666666666666666d+0
     1         *t129*t59*t76-2.6666666666666666d+0*t111*t59*t74-2.666666
     2         6666666666d+0*t145*t31*t59)+t109*(-2*t189*t76+2.0d+0*t179
     3         *t22*t74+4.0d+0*t111*t145*t22*t31)+5.333333333333333d+0*t
     4         109*t111*t31*t53*t59*t76+8.888888888888888d-1*t112*t31*t4
     5         9*t59*t76-1.3333333333333333d+0*t113*t114*t31*t49*t76-6.0
     6         d+0*t178*t179*t22*t31*t76-8.888888888888888d-1*t112*t28*t
     7         59*t74+1.3333333333333333d+0*t113*t114*t28*t74+2.66666666
     8         66666666d+0*t144*t28*t53*t59)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t32*(-1.0d+0*t154*t6*t99-1.0d+0*t
     1         18*t6*(-7.164408684517644d-7*t148*t66*t7*t78*t86+1.0d+1*t
     2         103*t36*t45*t8+2.0d+0*t11*t153*t36*t38*t8-2.0d+0*t152*t34
     3         *t35*t36*t8-6.0d+0*gammaaa*t139*t175*t36*t8)-1.0d+0*t150*
     4         t6*t79-2.0d+0*t140*t6*t71-2.0d+0*t141*t6*t69)+1.0d+0*t13*
     5         t6*(-1.3333333333333333d+1*gammaaa*t139*t15*t17*t36*t91-6
     6         .666666666666666d-7*t148*t66*t7*t78*t86+2.666666666666666
     7         6d+0*t15*t17*t171*t36*t44/t5**3.3d+1+1.0666666666666666d+
     8         1*t15*t16*t17*t36*t45)+t77*(2.0d+0*t6*t69*t71*t79+2.0d+0*
     9         t154*t18*t6*t79+2.0d+0*t141*t18*t6*t71-2*t142*t71)-6.0d+0
     :         *t151*t164*t18*t6*t79+2.6666666666666666d+0*t151*t18*t36*
     ;         t42*t77+t32*t36*(-2.6666666666666666d+0*t42*t69*t71-1.333
     <         3333333333333d+0*t154*t18*t42)+1.3333333333333333d+0*t13*
     =         t150*t36*t42)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (1.0d+0*t22*t28*(-1.33333333333333
     1         33d+1*gammabb*t143*t15*t30*t53*t91-6.666666666666666d-7*t
     2         110*t156*t23*t66*t86+1.0666666666666666d+1*t15*t16*t30*t5
     3         3*t62+2.6666666666666666d+0*t15*t171*t30*t53*t61/t21**3.3
     4         d+1)+t49*(-1.0d+0*t22*t31*(-7.164408684517644d-7*t110*t15
     5         6*t23*t66*t86+1.0d+1*t133*t53*t62*t8+2.0d+0*t161*t26*t53*
     6         t55*t8-2.0d+0*t160*t51*t52*t53*t8-6.0d+0*gammabb*t143*t18
     7         7*t53*t8)-2.0d+0*t144*t22*t76-2.0d+0*t145*t22*t74-1.0d+0*
     8         t129*t162*t22-1.0d+0*t111*t158*t22)+t49*t53*(-2.666666666
     9         6666666d+0*t59*t74*t76-1.3333333333333333d+0*t162*t31*t59
     :         )+t109*(2.0d+0*t111*t22*t74*t76+2.0d+0*t145*t22*t31*t76-2
     ;         *t146*t76+2.0d+0*t111*t162*t22*t31)+2.6666666666666666d+0
     <         *t109*t159*t31*t53*t59+1.3333333333333333d+0*t158*t28*t53
     =         *t59-6.0d+0*t111*t159*t178*t22*t31)*wght+Cmat3(iq,D3_RB_G
     >         BB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (1.0d+0*t13*t6*(3*t149*t15*t17*t9
     1         1+5.0d-7*t1*t147*t191*t66*t7-gammaaa*t15*t17*t171/t5**3.2
     2         d+1)+t32*(-1.0d+0*t18*t6*(2.25d+0*t10*t11*t8/t9**5-7.5d-1
     3         *t103*t152*t68*t8-2.25d+0*t14*t34*t8/t44+2.25d+0*t149*t17
     4         5*t8+5.373306513388233d-7*t1*t147*t191*t66*t7)-3.0d+0*t15
     5         0*t6*t71-3.0d+0*t154*t6*t69)+(4.0d+0*t154*t18*t6*t71-2*t1
     6         55*t71+2.0d+0*t151*t6*t69)*t77-6.0d+0*t164*t18*t6*t71**3)
     7         *wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (1.0d+0*t22*t28*(3*t15*t157*t30*t
     1         91+5.0d-7*t1*t147*t192*t23*t66-gammabb*t15*t171*t30/t21**
     2         3.2d+1)+t49*(-1.0d+0*t22*t31*(-7.5d-1*t133*t160*t73*t8-2.
     3         25d+0*t29*t51*t8/t61+2.25d+0*t25*t26*t8/t24**5+2.25d+0*t1
     4         57*t187*t8+5.373306513388233d-7*t1*t147*t192*t23*t66)-3.0
     5         d+0*t158*t22*t76-3.0d+0*t162*t22*t74)-6.0d+0*t178*t22*t31
     6         *t76**3+t109*(4.0d+0*t162*t22*t31*t76-2*t163*t76+2.0d+0*t
     7         159*t22*t74))*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammaaa**t2
            t4 = rhoa**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammaaa**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammaaa*t14*t16)
            t18 = -gammaaa*t14*t8+1.0d-6*t3*t7+gammaaa*t14*t15*t17
            t19 = 1/t12**2
            t20 = (gammaaa*t14+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t5**9.0d+0
            t23 = 1/rhoa**6.666666666666666d-1
            t24 = -8.0d+0*gammaaa*t21*t22*t23*t8
            t25 = 1/t5**5.0d+0
            t26 = -8.0d+0*t11*t23*t25*t8*t9
            t27 = 1.0d+0/t4
            t28 = -1.4328817369035288d-6*t1*t23*t27*t3*t7+t26+t24
            t29 = t5**3.0d+0
            t30 = 2.6666666666666666d+0*gammaaa*t22*t23*t8
            t31 = gammaaa**2
            t32 = 1/t5**1.7d+1
            t33 = 2.6666666666666666d+0*t15*t16*t17*t23*t31*t32
            t34 = -2.6666666666666666d+0*gammaaa*t15*t17*t22*t23
            t35 = -1.3333333333333333d-6*t1*t23*t27*t3*t7+t34+t33+t30
            t36 = t2-1
            t37 = gammaaa**t36
            t38 = 1/t5**1.6d+1
            t39 = -t14*t8+5.0d-7*t1*t37*t7-gammaaa*t15*t16*t17*t38+t14*t
     1         15*t17
            t40 = 1/t9
            t41 = 3.0d+0*t10*t11*t40*t8+3.0d+0*t14*t21*t8+5.373306513388
     1         233d-7*t1*t37*t7
            t42 = 1/t12**3
            t43 = 1/rhoa
            t44 = -1.4328817369035288d-6*t1*t3*t43*t7+t26+t24
            t45 = 1/rhoa**1.6666666666666669d+0
            t46 = t5**2.0d+0
            t47 = 1/rhoa**1.3333333333333333d+0
            t48 = 1/rhoa**2
            t49 = 1.3333333333333333d-6*t1*t3*t48*t7
            t50 = -1.7777777777777776d+0*gammaaa*t22*t45*t8
            t51 = t1**2
            t52 = 1/t5**1.0d+1
            t53 = -8.0d+0*gammaaa*t47*t52*t8
            t54 = -1.7777777777777776d+0*t15*t16*t17*t31*t32*t45
            t55 = 1.7777777777777776d+0*gammaaa*t15*t17*t22*t45
            t56 = t16**2
            t57 = gammaaa**3
            t58 = 1/t5**2.6d+1
            t59 = 7.11111111111111d+0*t15*t17*t47*t56*t57*t58
            t60 = 1/t5**1.8d+1
            t61 = -2.222222222222222d+1*t15*t16*t17*t31*t47*t60
            t62 = 8.0d+0*gammaaa*t15*t17*t47*t52
            t63 = 1.7777777777777776d-6*t27*t3*t45*t51*t7+t62+t61+t59+t5
     1         5+t54+t53+t50+t49
            t64 = -1.3333333333333333d-6*t1*t3*t43*t7+t34+t33+t30
            t65 = 1.4328817369035288d-6*t1*t3*t48*t7
            t66 = 5.333333333333333d+0*gammaaa*t21*t22*t45*t8
            t67 = 5.333333333333333d+0*t11*t25*t45*t8*t9
            t68 = 1/t20**3
            t69 = -1.0666666666666666d+1*t31*t47*t60*t68*t8
            t70 = 3.466666666666666d+1*gammaaa*t21*t47*t52*t8
            t71 = 1/t5**6.0d+0
            t72 = 1.3333333333333333d+1*t11*t47*t71*t8*t9
            t73 = t72+t70+1.9105089825380384d-6*t27*t3*t45*t51*t7+t69+t6
     1         7+t66+t65
            t74 = 1/t5**2.5d+1
            t75 = 2.6666666666666666d+0*t22*t23*t8-2.6666666666666666d+0
     1         *t15*t17*t23*t31*t56*t74-6.666666666666666d-7*t37*t43*t51
     2         *t7+8.0d+0*gammaaa*t15*t16*t17*t23*t32-2.6666666666666666
     3         d+0*t15*t17*t22*t23
            t76 = 4.0d+0*gammaaa*t23*t32*t68*t8-4.0d+0*t11*t23*t25*t40*t
     1         8-1.2d+1*t21*t22*t23*t8-7.164408684517644d-7*t37*t43*t51*
     2         t7
            t77 = -1.0d+0*t18*t6*t76*wght-1.0d+0*t41*t6*t64*wght-1.0d+0*
     1         t39*t44*t6*wght
            t78 = t2-2
            t79 = gammaaa**t78
            t80 = 1/t5**2.4d+1
            t81 = gammaaa*t15*t17*t56*t80+5.0d-7*t1*t36*t7*t79-2*t15*t16
     1         *t17*t38
            t82 = t41**2
            t83 = 1/gammaaa
            t84 = 1/t9**3
            t85 = -1.5d+0*t10*t11*t8*t84+1.5d+0*t14*t21*t8*t83-1.5d+0*t3
     1         8*t68*t8+5.373306513388233d-7*t1*t36*t7*t79
            t86 = -1.0d+0*t18*t6*t85*wght-2.0d+0*t39*t41*t6*wght
            t87 = 1/t12**4
            t88 = t44**2
            t89 = 1/rhoa**2.6666666666666666d+0
            t90 = 1/rhoa**2.3333333333333334d+0
            t91 = 1/rhoa**3
            t92 = t1**3
            t93 = 1/t5**1.1d+1
            t94 = t16**3
            t95 = 1/t5**2.7d+1
            t96 = 1/t5**1.9d+1
            t97 = 1.7777777777777776d-6*t3*t48*t51*t7+t62+t61+t59+t55+t5
     1         4+t53+t50+t49
            t98 = 1/t20**5
            t99 = t72+t70+1.9105089825380384d-6*t3*t48*t51*t7+t69+t67+t6
     1         6+t65
            t100 = -1.0d+0*t18*t6*t99*wght-2.0d+0*t44*t6*t64*wght
            t101 = gammaaa**(t2-3)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*t13*t35*t6*wght-1.0d+0*t18*t19*t28*t
     1         6*wght+1.3333333333333333d+0*t13*t18*t23*t29*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = -1.0d+0*t18*t19*t41*t6*wght+1.0d+0*t13*t39
     1         *t6*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = t19*(-1.0d+0*t18*t6*t73*wght-1.0d+0*t28
     1         *t6*t64*wght-1.0d+0*t35*t44*t6*wght)+t13*t23*(1.333333333
     2         3333333d+0*t29*t64*wght+1.3333333333333333d+0*t29*t35*wgh
     3         t)+t19*t23*(-1.3333333333333333d+0*t18*t29*t44*wght-1.333
     4         3333333333333d+0*t18*t28*t29*wght)+1.0d+0*t13*t6*t63*wght
     5         +2.0d+0*t18*t28*t42*t44*t6*wght+1.3333333333333333d+0*t13
     6         *t18*t46*t47*wght-8.888888888888888d-1*t13*t18*t29*t45*wg
     7         ht+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 1.0d+0*t13*t6*t75*wght+2.0d+0*t18*t41*
     1         t42*t44*t6*wght-1.3333333333333333d+0*t18*t19*t23*t29*t41
     2         *wght+1.3333333333333333d+0*t13*t23*t29*t39*wght+t19*t77+
     3         Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 2.0d+0*t18*t42*t6*t82*wght+1.0d+0*t13
     1         *t6*t81*wght+t19*t86+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t19*(-1.0d+0*t35*t6*t99*wght-1.33333
     1         33333333333d+0*t18*t23*t29*t99*wght-1.0d+0*t18*t6*(-4.266
     2         6666666666664d+1*t48*t57*t8*t95*t98+1.1022222222222221d+2
     3         *t31*t48*t68*t8*t96-1.333333333333333d+2*gammaaa*t21*t48*
     4         t8*t93-2.5473453100507176d-6*t27*t3*t7*t89*t92-3.82101796
     5         50760767d-6*t3*t51*t7*t91-2.8657634738070575d-6*t1*t3*t7*
     6         t91-2.6666666666666666d+1*t11*t71*t8*t9*t90+2.13333333333
     7         33332d+1*t31*t60*t68*t8*t90-6.933333333333332d+1*gammaaa*
     8         t21*t52*t8*t90-8.88888888888889d+0*t11*t25*t8*t89*t9-2.66
     9         66666666666666d+1*t11*t48*t8*t9/t5**7.0d+0-8.888888888888
     :         89d+0*gammaaa*t21*t22*t8*t89-1.9105089825380384d-6*t27*t3
     ;         *t51*t7*t89)*wght-1.0d+0*t28*t6*t97*wght-2.0d+0*t6*t64*t7
     <         3*wght-2.6666666666666666d+0*t23*t29*t44*t64*wght-2.0d+0*
     =         t44*t6*t63*wght)+t13*t23*(1.3333333333333333d+0*t29*t97*w
     >         ght+2.6666666666666666d+0*t23*t46*t64*wght+2.666666666666
     ?         6666d+0*t29*t63*wght)+t42*(2.0d+0*t35*t6*t88*wght+4.0d+0*
     @         t18*t44*t6*t73*wght-2*t100*t28)+t23*t42*(2.66666666666666
     1         66d+0*t18*t29*t88*wght+5.333333333333333d+0*t18*t28*t29*t
     2         44*wght)+t19*t23*(-2.6666666666666666d+0*t18*t29*t73*wght
     3         -2.6666666666666666d+0*t28*t29*t64*wght-2.666666666666666
     4         6d+0*t18*t23*t44*t46*wght-2.6666666666666666d+0*t29*t35*t
     5         44*wght)+1.0d+0*t13*t6*(1.5466666666666665d+2*t15*t16*t17
     6         *t31*t48*t96-1.2088888888888887d+2*t15*t17*t48*t56*t57*t9
     7         5+1.8962962962962962d+1*gammaaa**4*t15*t17*t48*t94/t5**3.
     8         5d+1+2.6666666666666666d+1*gammaaa*t48*t8*t93-2.666666666
     9         6666666d+1*gammaaa*t15*t17*t48*t93-2.37037037037037d-6*t2
     :         7*t3*t7*t89*t92-3.555555555555555d-6*t3*t51*t7*t91-2.6666
     ;         666666666666d-6*t1*t3*t7*t91+1.6d+1*gammaaa*t52*t8*t90+4.
     <         444444444444444d+1*t15*t16*t17*t31*t60*t90-1.422222222222
     =         2222d+1*t15*t17*t56*t57*t58*t90-1.6d+1*gammaaa*t15*t17*t5
     >         2*t90+2.962962962962963d+0*gammaaa*t22*t8*t89-1.777777777
     ?         7777776d-6*t27*t3*t51*t7*t89+2.962962962962963d+0*t15*t16
     @         *t17*t31*t32*t89-2.962962962962963d+0*gammaaa*t15*t17*t22
     1         *t89)*wght-2.6666666666666666d+0*t13*t18*t46*t90*wght+1.4
     2         814814814814814d+0*t13*t18*t29*t89*wght-6.0d+0*t18*t28*t6
     3         *t87*t88*wght-1.7777777777777776d+0*t13*t29*t45*t64*wght+
     4         1.3333333333333333d+0*t13*t35*t46*t47*wght-1.333333333333
     5         3333d+0*t18*t19*t28*t46*t47*wght+1.7777777777777776d+0*t1
     6         8*t19*t29*t44*t45*wght-8.888888888888888d-1*t13*t29*t35*t
     7         45*wght+8.888888888888888d-1*t18*t19*t28*t29*t45*wght+8.8
     8         88888888888888d-1*t13*t18*t45*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t19*(-1.0d+0*t39*t6*t99*wght-1.0d+0
     1         *t18*t6*(1.6d+1*t31*t47*t58*t8*t98+9.552544912690191d-7*t
     2         37*t48*t7*t92+6.666666666666666d+0*t11*t40*t47*t71*t8-3.8
     3         666666666666666d+1*gammaaa*t47*t60*t68*t8-2.6666666666666
     4         666d+0*gammaaa*t32*t45*t68*t8+4.133333333333333d+1*t21*t4
     5         7*t52*t8+2.6666666666666666d+0*t11*t25*t40*t45*t8+8.0d+0*
     6         t21*t22*t45*t8+7.164408684517644d-7*t37*t48*t51*t7)*wght-
     7         1.0d+0*t41*t6*t97*wght-2.0d+0*t6*t64*t76*wght-2.0d+0*t44*
     8         t6*t75*wght)+t42*(2.0d+0*t39*t6*t88*wght+4.0d+0*t18*t44*t
     9         6*t76*wght-2*t100*t41)+t19*t23*(-2.6666666666666666d+0*t1
     :         8*t29*t76*wght-2.6666666666666666d+0*t29*t41*t64*wght-2.6
     ;         666666666666666d+0*t29*t39*t44*wght)+1.0d+0*t13*t6*(-7.11
     <         111111111111d+0*t15*t17*t47*t57*t94/t5**3.4d+1+8.88888888
     =         8888887d-7*t37*t48*t7*t92-8.0d+0*t47*t52*t8-1.77777777777
     >         77776d+0*t22*t45*t8+1.7777777777777776d+0*t15*t17*t31*t45
     ?         *t56*t74+6.666666666666666d-7*t37*t48*t51*t7-5.2444444444
     @         44444d+1*gammaaa*t15*t16*t17*t47*t60+4.355555555555556d+1
     1         *t15*t17*t31*t47*t56*t58+8.0d+0*t15*t17*t47*t52-5.3333333
     2         33333333d+0*gammaaa*t15*t16*t17*t32*t45+1.777777777777777
     3         6d+0*t15*t17*t22*t45)*wght-6.0d+0*t18*t41*t6*t87*t88*wght
     4         +2.6666666666666666d+0*t13*t23*t29*t75*wght-1.33333333333
     5         33333d+0*t18*t19*t41*t46*t47*wght+1.3333333333333333d+0*t
     6         13*t39*t46*t47*wght+8.888888888888888d-1*t18*t19*t29*t41*
     7         t45*wght-8.888888888888888d-1*t13*t29*t39*t45*wght+5.3333
     8         33333333333d+0*t18*t23*t29*t41*t42*t44*wght+Cmat3(iq,D3_R
     9         A_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t19*(-1.0d+0*t18*t6*(-6.0d+0*gamma
     1         aa*t23*t74*t8*t98+2.0d+0*t11*t23*t25*t8*t84-2.0d+0*t21*t2
     2         2*t23*t8*t83+1.0d+1*t23*t32*t68*t8-7.164408684517644d-7*t
     3         36*t43*t51*t7*t79)*wght-1.0d+0*t6*t64*t85*wght-1.0d+0*t44
     4         *t6*t81*wght-2.0d+0*t39*t6*t76*wght-2.0d+0*t41*t6*t75*wgh
     5         t)+t42*(2.0d+0*t18*t44*t6*t85*wght+2.0d+0*t18*t41*t6*t76*
     6         wght+2.0d+0*t39*t41*t44*t6*wght-2*t41*t77)+t19*t23*(-1.33
     7         33333333333333d+0*t18*t29*t85*wght-2.6666666666666666d+0*
     8         t29*t39*t41*wght)+1.0d+0*t13*t6*(2.6666666666666666d+0*t1
     9         5*t17*t23*t31*t94/t5**3.3d+1-6.666666666666666d-7*t36*t43
     :         *t51*t7*t79-1.3333333333333333d+1*gammaaa*t15*t17*t23*t56
     ;         *t74+1.0666666666666666d+1*t15*t16*t17*t23*t32)*wght-6.0d
     <         +0*t18*t44*t6*t82*t87*wght+2.6666666666666666d+0*t18*t23*
     =         t29*t42*t82*wght+1.3333333333333333d+0*t13*t23*t29*t81*wg
     >         ht+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t19*(-1.0d+0*t18*t6*(2.25d+0*t8*t
     1         80*t98+2.25d+0*t10*t11*t8/t9**5-7.5d-1*t38*t68*t8*t83-2.2
     2         5d+0*t14*t21*t8/t31+5.373306513388233d-7*t1*t101*t36*t7*t
     3         78)*wght-3.0d+0*t39*t6*t85*wght-3.0d+0*t41*t6*t81*wght)+t
     4         42*(4.0d+0*t18*t41*t6*t85*wght+2.0d+0*t39*t6*t82*wght-2*t
     5         41*t86)+1.0d+0*t13*t6*(-gammaaa*t15*t17*t94/t5**3.2d+1+3*
     6         t15*t17*t56*t80+5.0d-7*t1*t101*t36*t7*t78)*wght-6.0d+0*t1
     7         8*t41**3*t6*t87*wght+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(3)
            t2 = 5.0d-1*t1
            t3 = gammabb**t2
            t4 = rhob**3.333333333333333d-1
            t5 = 1.0d+0*t4
            t6 = t5**4.0d+0
            t7 = 1/t6**t1
            t8 = param(1)
            t9 = gammabb**5.0d-1
            t10 = 1/t6
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t3*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t5**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = param(2)
            t17 = exp(-gammabb*t14*t16)
            t18 = -gammabb*t14*t8+1.0d-6*t3*t7+gammabb*t14*t15*t17
            t19 = 1/t12**2
            t20 = (gammabb*t14+1)**5.0d-1
            t21 = 1/t20
            t22 = 1/t5**9.0d+0
            t23 = 1/rhob**6.666666666666666d-1
            t24 = -8.0d+0*gammabb*t21*t22*t23*t8
            t25 = 1/t5**5.0d+0
            t26 = -8.0d+0*t11*t23*t25*t8*t9
            t27 = 1.0d+0/t4
            t28 = -1.4328817369035288d-6*t1*t23*t27*t3*t7+t26+t24
            t29 = t5**3.0d+0
            t30 = 2.6666666666666666d+0*gammabb*t22*t23*t8
            t31 = gammabb**2
            t32 = 1/t5**1.7d+1
            t33 = 2.6666666666666666d+0*t15*t16*t17*t23*t31*t32
            t34 = -2.6666666666666666d+0*gammabb*t15*t17*t22*t23
            t35 = -1.3333333333333333d-6*t1*t23*t27*t3*t7+t34+t33+t30
            t36 = t2-1
            t37 = gammabb**t36
            t38 = 1/t5**1.6d+1
            t39 = -t14*t8+5.0d-7*t1*t37*t7-gammabb*t15*t16*t17*t38+t14*t
     1         15*t17
            t40 = 1/t9
            t41 = 3.0d+0*t10*t11*t40*t8+3.0d+0*t14*t21*t8+5.373306513388
     1         233d-7*t1*t37*t7
            t42 = 1/t12**3
            t43 = 1/rhob
            t44 = -1.4328817369035288d-6*t1*t3*t43*t7+t26+t24
            t45 = 1/rhob**1.6666666666666669d+0
            t46 = t5**2.0d+0
            t47 = 1/rhob**1.3333333333333333d+0
            t48 = 1/rhob**2
            t49 = 1.3333333333333333d-6*t1*t3*t48*t7
            t50 = -1.7777777777777776d+0*gammabb*t22*t45*t8
            t51 = t1**2
            t52 = 1/t5**1.0d+1
            t53 = -8.0d+0*gammabb*t47*t52*t8
            t54 = -1.7777777777777776d+0*t15*t16*t17*t31*t32*t45
            t55 = 1.7777777777777776d+0*gammabb*t15*t17*t22*t45
            t56 = t16**2
            t57 = gammabb**3
            t58 = 1/t5**2.6d+1
            t59 = 7.11111111111111d+0*t15*t17*t47*t56*t57*t58
            t60 = 1/t5**1.8d+1
            t61 = -2.222222222222222d+1*t15*t16*t17*t31*t47*t60
            t62 = 8.0d+0*gammabb*t15*t17*t47*t52
            t63 = 1.7777777777777776d-6*t27*t3*t45*t51*t7+t62+t61+t59+t5
     1         5+t54+t53+t50+t49
            t64 = -1.3333333333333333d-6*t1*t3*t43*t7+t34+t33+t30
            t65 = 1.4328817369035288d-6*t1*t3*t48*t7
            t66 = 5.333333333333333d+0*gammabb*t21*t22*t45*t8
            t67 = 5.333333333333333d+0*t11*t25*t45*t8*t9
            t68 = 1/t20**3
            t69 = -1.0666666666666666d+1*t31*t47*t60*t68*t8
            t70 = 3.466666666666666d+1*gammabb*t21*t47*t52*t8
            t71 = 1/t5**6.0d+0
            t72 = 1.3333333333333333d+1*t11*t47*t71*t8*t9
            t73 = t72+t70+1.9105089825380384d-6*t27*t3*t45*t51*t7+t69+t6
     1         7+t66+t65
            t74 = 1/t5**2.5d+1
            t75 = 2.6666666666666666d+0*t22*t23*t8-2.6666666666666666d+0
     1         *t15*t17*t23*t31*t56*t74-6.666666666666666d-7*t37*t43*t51
     2         *t7+8.0d+0*gammabb*t15*t16*t17*t23*t32-2.6666666666666666
     3         d+0*t15*t17*t22*t23
            t76 = 4.0d+0*gammabb*t23*t32*t68*t8-4.0d+0*t11*t23*t25*t40*t
     1         8-1.2d+1*t21*t22*t23*t8-7.164408684517644d-7*t37*t43*t51*
     2         t7
            t77 = -1.0d+0*t18*t6*t76*wght-1.0d+0*t41*t6*t64*wght-1.0d+0*
     1         t39*t44*t6*wght
            t78 = t2-2
            t79 = gammabb**t78
            t80 = 1/t5**2.4d+1
            t81 = gammabb*t15*t17*t56*t80+5.0d-7*t1*t36*t7*t79-2*t15*t16
     1         *t17*t38
            t82 = t41**2
            t83 = 1/gammabb
            t84 = 1/t9**3
            t85 = -1.5d+0*t10*t11*t8*t84+1.5d+0*t14*t21*t8*t83-1.5d+0*t3
     1         8*t68*t8+5.373306513388233d-7*t1*t36*t7*t79
            t86 = -1.0d+0*t18*t6*t85*wght-2.0d+0*t39*t41*t6*wght
            t87 = 1/t12**4
            t88 = t44**2
            t89 = 1/rhob**2.6666666666666666d+0
            t90 = 1/rhob**2.3333333333333334d+0
            t91 = 1/rhob**3
            t92 = t1**3
            t93 = 1/t5**1.1d+1
            t94 = t16**3
            t95 = 1/t5**2.7d+1
            t96 = 1/t5**1.9d+1
            t97 = 1.7777777777777776d-6*t3*t48*t51*t7+t62+t61+t59+t55+t5
     1         4+t53+t50+t49
            t98 = 1/t20**5
            t99 = t72+t70+1.9105089825380384d-6*t3*t48*t51*t7+t69+t67+t6
     1         6+t65
            t100 = -1.0d+0*t18*t6*t99*wght-2.0d+0*t44*t6*t64*wght
            t101 = gammabb**(t2-3)
            fnc(iq) = 1.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RB) = 1.0d+0*t13*t35*t6*wght-1.0d+0*t18*t19*t28*t
     1         6*wght+1.3333333333333333d+0*t13*t18*t23*t29*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GBB) = -1.0d+0*t18*t19*t41*t6*wght+1.0d+0*t13*t39
     1         *t6*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = t19*(-1.0d+0*t18*t6*t73*wght-1.0d+0*t28
     1         *t6*t64*wght-1.0d+0*t35*t44*t6*wght)+t13*t23*(1.333333333
     2         3333333d+0*t29*t64*wght+1.3333333333333333d+0*t29*t35*wgh
     3         t)+t19*t23*(-1.3333333333333333d+0*t18*t29*t44*wght-1.333
     4         3333333333333d+0*t18*t28*t29*wght)+1.0d+0*t13*t6*t63*wght
     5         +2.0d+0*t18*t28*t42*t44*t6*wght+1.3333333333333333d+0*t13
     6         *t18*t46*t47*wght-8.888888888888888d-1*t13*t18*t29*t45*wg
     7         ht+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 1.0d+0*t13*t6*t75*wght+2.0d+0*t18*t41*
     1         t42*t44*t6*wght-1.3333333333333333d+0*t18*t19*t23*t29*t41
     2         *wght+1.3333333333333333d+0*t13*t23*t29*t39*wght+t19*t77+
     3         Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 2.0d+0*t18*t42*t6*t82*wght+1.0d+0*t13
     1         *t6*t81*wght+t19*t86+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t19*(-1.0d+0*t35*t6*t99*wght-1.33333
     1         33333333333d+0*t18*t23*t29*t99*wght-1.0d+0*t18*t6*(-4.266
     2         6666666666664d+1*t48*t57*t8*t95*t98+1.1022222222222221d+2
     3         *t31*t48*t68*t8*t96-1.333333333333333d+2*gammabb*t21*t48*
     4         t8*t93-2.5473453100507176d-6*t27*t3*t7*t89*t92-3.82101796
     5         50760767d-6*t3*t51*t7*t91-2.8657634738070575d-6*t1*t3*t7*
     6         t91-2.6666666666666666d+1*t11*t71*t8*t9*t90+2.13333333333
     7         33332d+1*t31*t60*t68*t8*t90-6.933333333333332d+1*gammabb*
     8         t21*t52*t8*t90-8.88888888888889d+0*t11*t25*t8*t89*t9-2.66
     9         66666666666666d+1*t11*t48*t8*t9/t5**7.0d+0-8.888888888888
     :         89d+0*gammabb*t21*t22*t8*t89-1.9105089825380384d-6*t27*t3
     ;         *t51*t7*t89)*wght-1.0d+0*t28*t6*t97*wght-2.0d+0*t6*t64*t7
     <         3*wght-2.6666666666666666d+0*t23*t29*t44*t64*wght-2.0d+0*
     =         t44*t6*t63*wght)+t13*t23*(1.3333333333333333d+0*t29*t97*w
     >         ght+2.6666666666666666d+0*t23*t46*t64*wght+2.666666666666
     ?         6666d+0*t29*t63*wght)+t42*(2.0d+0*t35*t6*t88*wght+4.0d+0*
     @         t18*t44*t6*t73*wght-2*t100*t28)+t23*t42*(2.66666666666666
     1         66d+0*t18*t29*t88*wght+5.333333333333333d+0*t18*t28*t29*t
     2         44*wght)+t19*t23*(-2.6666666666666666d+0*t18*t29*t73*wght
     3         -2.6666666666666666d+0*t28*t29*t64*wght-2.666666666666666
     4         6d+0*t18*t23*t44*t46*wght-2.6666666666666666d+0*t29*t35*t
     5         44*wght)+1.0d+0*t13*t6*(1.5466666666666665d+2*t15*t16*t17
     6         *t31*t48*t96-1.2088888888888887d+2*t15*t17*t48*t56*t57*t9
     7         5+1.8962962962962962d+1*gammabb**4*t15*t17*t48*t94/t5**3.
     8         5d+1+2.6666666666666666d+1*gammabb*t48*t8*t93-2.666666666
     9         6666666d+1*gammabb*t15*t17*t48*t93-2.37037037037037d-6*t2
     :         7*t3*t7*t89*t92-3.555555555555555d-6*t3*t51*t7*t91-2.6666
     ;         666666666666d-6*t1*t3*t7*t91+1.6d+1*gammabb*t52*t8*t90+4.
     <         444444444444444d+1*t15*t16*t17*t31*t60*t90-1.422222222222
     =         2222d+1*t15*t17*t56*t57*t58*t90-1.6d+1*gammabb*t15*t17*t5
     >         2*t90+2.962962962962963d+0*gammabb*t22*t8*t89-1.777777777
     ?         7777776d-6*t27*t3*t51*t7*t89+2.962962962962963d+0*t15*t16
     @         *t17*t31*t32*t89-2.962962962962963d+0*gammabb*t15*t17*t22
     1         *t89)*wght-2.6666666666666666d+0*t13*t18*t46*t90*wght+1.4
     2         814814814814814d+0*t13*t18*t29*t89*wght-6.0d+0*t18*t28*t6
     3         *t87*t88*wght-1.7777777777777776d+0*t13*t29*t45*t64*wght+
     4         1.3333333333333333d+0*t13*t35*t46*t47*wght-1.333333333333
     5         3333d+0*t18*t19*t28*t46*t47*wght+1.7777777777777776d+0*t1
     6         8*t19*t29*t44*t45*wght-8.888888888888888d-1*t13*t29*t35*t
     7         45*wght+8.888888888888888d-1*t18*t19*t28*t29*t45*wght+8.8
     8         88888888888888d-1*t13*t18*t45*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t19*(-1.0d+0*t39*t6*t99*wght-1.0d+0
     1         *t18*t6*(1.6d+1*t31*t47*t58*t8*t98+9.552544912690191d-7*t
     2         37*t48*t7*t92+6.666666666666666d+0*t11*t40*t47*t71*t8-3.8
     3         666666666666666d+1*gammabb*t47*t60*t68*t8-2.6666666666666
     4         666d+0*gammabb*t32*t45*t68*t8+4.133333333333333d+1*t21*t4
     5         7*t52*t8+2.6666666666666666d+0*t11*t25*t40*t45*t8+8.0d+0*
     6         t21*t22*t45*t8+7.164408684517644d-7*t37*t48*t51*t7)*wght-
     7         1.0d+0*t41*t6*t97*wght-2.0d+0*t6*t64*t76*wght-2.0d+0*t44*
     8         t6*t75*wght)+t42*(2.0d+0*t39*t6*t88*wght+4.0d+0*t18*t44*t
     9         6*t76*wght-2*t100*t41)+t19*t23*(-2.6666666666666666d+0*t1
     :         8*t29*t76*wght-2.6666666666666666d+0*t29*t41*t64*wght-2.6
     ;         666666666666666d+0*t29*t39*t44*wght)+1.0d+0*t13*t6*(-7.11
     <         111111111111d+0*t15*t17*t47*t57*t94/t5**3.4d+1+8.88888888
     =         8888887d-7*t37*t48*t7*t92-8.0d+0*t47*t52*t8-1.77777777777
     >         77776d+0*t22*t45*t8+1.7777777777777776d+0*t15*t17*t31*t45
     ?         *t56*t74+6.666666666666666d-7*t37*t48*t51*t7-5.2444444444
     @         44444d+1*gammabb*t15*t16*t17*t47*t60+4.355555555555556d+1
     1         *t15*t17*t31*t47*t56*t58+8.0d+0*t15*t17*t47*t52-5.3333333
     2         33333333d+0*gammabb*t15*t16*t17*t32*t45+1.777777777777777
     3         6d+0*t15*t17*t22*t45)*wght-6.0d+0*t18*t41*t6*t87*t88*wght
     4         +2.6666666666666666d+0*t13*t23*t29*t75*wght-1.33333333333
     5         33333d+0*t18*t19*t41*t46*t47*wght+1.3333333333333333d+0*t
     6         13*t39*t46*t47*wght+8.888888888888888d-1*t18*t19*t29*t41*
     7         t45*wght-8.888888888888888d-1*t13*t29*t39*t45*wght+5.3333
     8         33333333333d+0*t18*t23*t29*t41*t42*t44*wght+Cmat3(iq,D3_R
     9         B_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t19*(-1.0d+0*t18*t6*(-6.0d+0*gamma
     1         bb*t23*t74*t8*t98+2.0d+0*t11*t23*t25*t8*t84-2.0d+0*t21*t2
     2         2*t23*t8*t83+1.0d+1*t23*t32*t68*t8-7.164408684517644d-7*t
     3         36*t43*t51*t7*t79)*wght-1.0d+0*t6*t64*t85*wght-1.0d+0*t44
     4         *t6*t81*wght-2.0d+0*t39*t6*t76*wght-2.0d+0*t41*t6*t75*wgh
     5         t)+t42*(2.0d+0*t18*t44*t6*t85*wght+2.0d+0*t18*t41*t6*t76*
     6         wght+2.0d+0*t39*t41*t44*t6*wght-2*t41*t77)+t19*t23*(-1.33
     7         33333333333333d+0*t18*t29*t85*wght-2.6666666666666666d+0*
     8         t29*t39*t41*wght)+1.0d+0*t13*t6*(2.6666666666666666d+0*t1
     9         5*t17*t23*t31*t94/t5**3.3d+1-6.666666666666666d-7*t36*t43
     :         *t51*t7*t79-1.3333333333333333d+1*gammabb*t15*t17*t23*t56
     ;         *t74+1.0666666666666666d+1*t15*t16*t17*t23*t32)*wght-6.0d
     <         +0*t18*t44*t6*t82*t87*wght+2.6666666666666666d+0*t18*t23*
     =         t29*t42*t82*wght+1.3333333333333333d+0*t13*t23*t29*t81*wg
     >         ht+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t19*(-1.0d+0*t18*t6*(2.25d+0*t8*t
     1         80*t98+2.25d+0*t10*t11*t8/t9**5-7.5d-1*t38*t68*t8*t83-2.2
     2         5d+0*t14*t21*t8/t31+5.373306513388233d-7*t1*t101*t36*t7*t
     3         78)*wght-3.0d+0*t39*t6*t85*wght-3.0d+0*t41*t6*t81*wght)+t
     4         42*(4.0d+0*t18*t41*t6*t85*wght+2.0d+0*t39*t6*t82*wght-2*t
     5         41*t86)+1.0d+0*t13*t6*(-gammabb*t15*t17*t94/t5**3.2d+1+3*
     6         t15*t17*t56*t80+5.0d-7*t1*t101*t36*t7*t78)*wght-6.0d+0*t1
     7         8*t41**3*t6*t87*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
