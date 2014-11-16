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
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}\,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}-{\it t_3}\,{\it t_{19}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}\,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{\it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}\,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}-{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw6(param,tol_rho,ipol,nq,wght, 
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
            t25 = 2.0d+0*t13*t6*(-t14*t8+5.0d-7*t1*t24*t7-gammaaa*t15*t1
     1         6*t17/t5**1.6d+1+t14*t15*t17)*wght-2.0d+0*t18*t19*t6*(3.0
     2         d+0*t10*t11*t8/t9+3.0d+0*t14*t20*t8+5.373306513388233d-7*
     3         t1*t24*t7)*wght
            fnc(iq) = 2.0d+0*t13*t18*t6*wght+fnc(iq)
            Amat(iq,D1_RA) = -2.0d+0*t18*t19*t6*(-8.0d+0*t11*t22*t8*t9/t
     1         5**5.0d+0-8.0d+0*gammaaa*t20*t21*t22*t8-1.432881736903528
     2         8d-6*t1*t22*t23*t3*t7)*wght+2.0d+0*t13*t6*(2.666666666666
     3         6666d+0*gammaaa*t21*t22*t8-1.3333333333333333d-6*t1*t22*t
     4         23*t3*t7+2.6666666666666666d+0*gammaaa**2*t15*t16*t17*t22
     5         /t5**1.7d+1-2.6666666666666666d+0*gammaaa*t15*t17*t21*t22
     6         )*wght+2.6666666666666666d+0*t13*t18*t22*t5**3.0d+0*wght+
     7         Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t25+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t25+Cmat(iq,D1_GAB)
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
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}\,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}-{\it t_3}\,{\it t_{19}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}\,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{\it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}\,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}-{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw6_d2(param,tol_rho,ipol,nq,wght, 
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
C>   f &=& {{1.0\,{\it t_{14}}\,\left({{{\it t_{11}}\,{\it t_{19}}\,\sigma_{\beta\beta}}\over{e^{{\it t_{12}}\,{\it t_{19}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}-{\it t_3}\,{\it t_{19}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{17}}\,{\it t_{18}}+6.0\,{\it t_3}\,{\it t_{15}}\,{\rm asinh}\; \left({\it t_{15}}\,{\it t_{16}}\right)\,{\it t_{16}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{11}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}\over{e^{{\it t_{12}}\,{\it t_{10}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_7}\,{\it t_9}-{\it t_3}\,{\it t_{10}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_7}\,{\it t_9}+6.0\,{\it t_3}\,{\it t_4}\,{\rm asinh}\; \left({\it t_4}\,{\it t_5}\right)\,{\it t_5}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{21}}\,\left({{{\it t_{11}}\,{\it t_{26}}\,\sigma_{ss}}\over{e^{{\it t_{12}}\,{\it t_{26}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}-{\it t_3}\,{\it t_{26}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{24}}\,{\it t_{25}}+6.0\,{\it t_3}\,{\it t_{22}}\,{\rm asinh}\; \left({\it t_{22}}\,{\it t_{23}}\right)\,{\it t_{23}}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] Y Zhao, DG Truhlar, J.Phys.Chem. A 109, 5656 (2005)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw6_d3(param,tol_rho,ipol,nq,wght, 
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
