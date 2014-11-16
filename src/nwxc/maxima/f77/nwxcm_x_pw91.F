C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_pw91.F 
C> The nwxcm_x_pw91 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_pw91 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}^{4.0}\\\\ 
C>   {\it t_5} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_7} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_2}^{{\it t_8}}}}\\\\ 
C>   {\it t_{10}} &=& {{{\it t_8}}\over{2}}\\\\ 
C>   {\it t_{11}} &=& \sigma_{\alpha\alpha}^{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{13}} &=& {\it t_5}-.001890381166699926\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_{15}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_4}^{{\it t_8}}}}\\\\ 
C>   {\it t_{17}} &=& \sigma_{\beta\beta}^{{\it t_{10}}}\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_3}^{8.0}}}\\\\ 
C>   {\it t_{19}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{20}} &=& {\it t_{19}}^{4.0}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{20}}}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{20}}^{{\it t_8}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^{{\it t_{10}}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{{\it t_{19}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, C. Fiolhais, Phys.Rev. B 46, 6671 (1992)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw91(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**3.333333333333333d-1
            t2 = 1.0d+0*t1
            t3 = t2**4.0d+0
            t4 = param(1)
            t5 = 5.0d-1*t4
            t6 = gammaaa**t5
            t7 = 1/t3**t4
            t8 = param(2)
            t9 = gammaaa**5.0d-1
            t10 = 1/t3
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t6*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t2**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = exp(-1.6455d+0*gammaaa*t14)
            t17 = -gammaaa*t14*t8+1.0d-6*t6*t7+gammaaa*t14*t15*t16
            t18 = t2**3.0d+0
            t19 = 1/rhoa**6.666666666666666d-1
            t20 = 1/t12**2
            t21 = 1/(gammaaa*t14+1)**5.0d-1
            t22 = 1/t2**9.0d+0
            t23 = 1.0d+0/t1
            t24 = gammaaa**(t5-1)
            t25 = (2.0d+0*t13*t3*(-t14*t8+5.0d-7*t24*t4*t7-1.6455d+0*gam
     1         maaa*t15*t16/t2**1.6d+1+t14*t15*t16)-2.0d+0*t17*t20*t3*(3
     2         .0d+0*t10*t11*t8/t9+3.0d+0*t14*t21*t8+5.373306513388233d-
     3         7*t24*t4*t7))*wght
            fnc(iq) = (2.0d+0*t13*t17*t3-1.8610514726982d+0*t3)*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = (-2.0d+0*t17*t20*t3*(-8.0d+0*t11*t19*t8*t9/
     1         t2**5.0d+0-8.0d+0*gammaaa*t19*t21*t22*t8-1.43288173690352
     2         88d-6*t19*t23*t4*t6*t7)+2.0d+0*t13*t3*(2.6666666666666666
     3         d+0*gammaaa*t19*t22*t8-1.3333333333333333d-6*t19*t23*t4*t
     4         6*t7-2.6666666666666666d+0*gammaaa*t15*t16*t19*t22+4.3879
     5         999999999997d+0*gammaaa**2*t15*t16*t19/t2**1.7d+1)+2.6666
     6         666666666666d+0*t13*t17*t18*t19-2.4814019635975998d+0*t18
     7         *t19)*wght+Amat(iq,D1_RA)
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
            t1 = rhoa**3.333333333333333d-1
            t2 = 1.0d+0*t1
            t3 = t2**4.0d+0
            t4 = param(1)
            t5 = 5.0d-1*t4
            t6 = gammaaa**t5
            t7 = 1/t3**t4
            t8 = param(2)
            t9 = gammaaa**5.0d-1
            t10 = 1/t3
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t6*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t2**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = exp(-1.6455d+0*gammaaa*t14)
            t17 = -gammaaa*t14*t8+1.0d-6*t6*t7+gammaaa*t14*t15*t16
            t18 = rhob**3.333333333333333d-1
            t19 = 1.0d+0*t18
            t20 = t19**4.0d+0
            t21 = gammabb**t5
            t22 = 1/t20**t4
            t23 = gammabb**5.0d-1
            t24 = 1/t20
            t25 = asinh(t23*t24)
            t26 = 6.0d+0*t23*t24*t25*t8+1.0746613026776465d-6*t21*t22+1.
     1         0d+0
            t27 = 1/t26
            t28 = 1/t19**8.0d+0
            t29 = exp(-1.6455d+0*gammabb*t28)
            t30 = -gammabb*t28*t8+gammabb*t15*t28*t29+1.0d-6*t21*t22
            t31 = t2**3.0d+0
            t32 = 1/rhoa**6.666666666666666d-1
            t33 = 1/t12**2
            t34 = 1/(gammaaa*t14+1)**5.0d-1
            t35 = 1/t2**9.0d+0
            t36 = 1.0d+0/t1
            t37 = t19**3.0d+0
            t38 = 1/rhob**6.666666666666666d-1
            t39 = 1/t26**2
            t40 = 1/(gammabb*t28+1)**5.0d-1
            t41 = 1/t19**9.0d+0
            t42 = 1.0d+0/t18
            t43 = t5-1
            t44 = gammaaa**t43
            t45 = gammabb**t43
            fnc(iq) = (1.0d+0*t20*t27*t30+1.0d+0*t13*t17*t3-9.3052573634
     1         91d-1*t3-9.305257363491d-1*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.0d+0*t17*t3*t33*(-8.0d+0*t11*t32*t8*t9/
     1         t2**5.0d+0-8.0d+0*gammaaa*t32*t34*t35*t8-1.43288173690352
     2         88d-6*t32*t36*t4*t6*t7)+1.0d+0*t13*t3*(2.6666666666666666
     3         d+0*gammaaa*t32*t35*t8-1.3333333333333333d-6*t32*t36*t4*t
     4         6*t7-2.6666666666666666d+0*gammaaa*t15*t16*t32*t35+4.3879
     5         999999999997d+0*gammaaa**2*t15*t16*t32/t2**1.7d+1)+1.3333
     6         333333333333d+0*t13*t17*t31*t32-1.2407009817987999d+0*t31
     7         *t32)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-1.0d+0*t20*t30*t39*(-8.0d+0*gammabb*t38*t
     1         40*t41*t8-8.0d+0*t23*t25*t38*t8/t19**5.0d+0-1.43288173690
     2         35288d-6*t21*t22*t38*t4*t42)+1.0d+0*t20*t27*(2.6666666666
     3         666666d+0*gammabb*t38*t41*t8-1.3333333333333333d-6*t21*t2
     4         2*t38*t4*t42-2.6666666666666666d+0*gammabb*t15*t29*t38*t4
     5         1+4.3879999999999997d+0*gammabb**2*t15*t29*t38/t19**1.7d+
     6         1)+1.3333333333333333d+0*t27*t30*t37*t38-1.24070098179879
     7         99d+0*t37*t38)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*(-t14*t8+5.0d-7*t4*t44*t7-1
     1         .6455d+0*gammaaa*t15*t16/t2**1.6d+1+t14*t15*t16)-1.0d+0*t
     2         17*t3*t33*(3.0d+0*t10*t11*t8/t9+3.0d+0*t14*t34*t8+5.37330
     3         6513388233d-7*t4*t44*t7))*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t20*t27*(-t28*t8+5.0d-7*t22*t4*t45
     1         +t15*t28*t29-1.6455d+0*gammabb*t15*t29/t19**1.6d+1)-1.0d+
     2         0*t20*t30*t39*(3.0d+0*t28*t40*t8+3.0d+0*t24*t25*t8/t23+5.
     3         373306513388233d-7*t22*t4*t45))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**3.333333333333333d-1
            t2 = 1.0d+0*t1
            t3 = t2**4.0d+0
            t4 = param(1)
            t5 = 5.0d-1*t4
            t6 = gammaaa**t5
            t7 = 1/t3**t4
            t8 = param(2)
            t9 = gammaaa**5.0d-1
            t10 = 1/t3
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t6*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t2**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = exp(-1.6455d+0*gammaaa*t14)
            t17 = -gammaaa*t14*t8+1.0d-6*t6*t7+gammaaa*t14*t15*t16
            t18 = t2**3.0d+0
            t19 = 1/rhoa**6.666666666666666d-1
            t20 = 1/t12**2
            t21 = 1/(gammaaa*t14+1)**5.0d-1
            t22 = 1/t2**9.0d+0
            t23 = 1.0d+0/t1
            t24 = gammaaa**(t5-1)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RA) = (-1.0d+0*t17*t20*t3*(-8.0d+0*t11*t19*t8*t9/
     1         t2**5.0d+0-8.0d+0*gammaaa*t19*t21*t22*t8-1.43288173690352
     2         88d-6*t19*t23*t4*t6*t7)+1.0d+0*t13*t3*(2.6666666666666666
     3         d+0*gammaaa*t19*t22*t8-1.3333333333333333d-6*t19*t23*t4*t
     4         6*t7-2.6666666666666666d+0*gammaaa*t15*t16*t19*t22+4.3879
     5         999999999997d+0*gammaaa**2*t15*t16*t19/t2**1.7d+1)+1.3333
     6         333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0*t18
     7         *t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*(-t14*t8+5.0d-7*t24*t4*t7-1
     1         .6455d+0*gammaaa*t15*t16/t2**1.6d+1+t14*t15*t16)-1.0d+0*t
     2         17*t20*t3*(3.0d+0*t10*t11*t8/t9+3.0d+0*t14*t21*t8+5.37330
     3         6513388233d-7*t24*t4*t7))*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**3.333333333333333d-1
            t2 = 1.0d+0*t1
            t3 = t2**4.0d+0
            t4 = param(1)
            t5 = 5.0d-1*t4
            t6 = gammabb**t5
            t7 = 1/t3**t4
            t8 = param(2)
            t9 = gammabb**5.0d-1
            t10 = 1/t3
            t11 = asinh(t10*t9)
            t12 = 6.0d+0*t10*t11*t8*t9+1.0746613026776465d-6*t6*t7+1.0d+
     1         0
            t13 = 1/t12
            t14 = 1/t2**8.0d+0
            t15 = t8-1.890381166699926d-3
            t16 = exp(-1.6455d+0*gammabb*t14)
            t17 = -gammabb*t14*t8+1.0d-6*t6*t7+gammabb*t14*t15*t16
            t18 = t2**3.0d+0
            t19 = 1/rhob**6.666666666666666d-1
            t20 = 1/t12**2
            t21 = 1/(gammabb*t14+1)**5.0d-1
            t22 = 1/t2**9.0d+0
            t23 = 1.0d+0/t1
            t24 = gammabb**(t5-1)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RB) = (-1.0d+0*t17*t20*t3*(-8.0d+0*t11*t19*t8*t9/
     1         t2**5.0d+0-8.0d+0*gammabb*t19*t21*t22*t8-1.43288173690352
     2         88d-6*t19*t23*t4*t6*t7)+1.0d+0*t13*t3*(2.6666666666666666
     3         d+0*gammabb*t19*t22*t8-1.3333333333333333d-6*t19*t23*t4*t
     4         6*t7-2.6666666666666666d+0*gammabb*t15*t16*t19*t22+4.3879
     5         999999999997d+0*gammabb**2*t15*t16*t19/t2**1.7d+1)+1.3333
     6         333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0*t18
     7         *t19)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = (1.0d+0*t13*t3*(-t14*t8+5.0d-7*t24*t4*t7-1
     1         .6455d+0*gammabb*t15*t16/t2**1.6d+1+t14*t15*t16)-1.0d+0*t
     2         17*t20*t3*(3.0d+0*t10*t11*t8/t9+3.0d+0*t14*t21*t8+5.37330
     3         6513388233d-7*t24*t4*t7))*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pw91 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}^{4.0}\\\\ 
C>   {\it t_5} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_7} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_2}^{{\it t_8}}}}\\\\ 
C>   {\it t_{10}} &=& {{{\it t_8}}\over{2}}\\\\ 
C>   {\it t_{11}} &=& \sigma_{\alpha\alpha}^{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{13}} &=& {\it t_5}-.001890381166699926\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_{15}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_4}^{{\it t_8}}}}\\\\ 
C>   {\it t_{17}} &=& \sigma_{\beta\beta}^{{\it t_{10}}}\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_3}^{8.0}}}\\\\ 
C>   {\it t_{19}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{20}} &=& {\it t_{19}}^{4.0}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{20}}}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{20}}^{{\it t_8}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^{{\it t_{10}}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{{\it t_{19}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, C. Fiolhais, Phys.Rev. B 46, 6671 (1992)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw91_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_pw91 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& 1.0\,\rho_\alpha^{{{1}\over{3}}}\\\\ 
C>   {\it t_2} &=& {\it t_1}^{4.0}\\\\ 
C>   {\it t_3} &=& 1.0\,\rho_\beta^{{{1}\over{3}}}\\\\ 
C>   {\it t_4} &=& {\it t_3}^{4.0}\\\\ 
C>   {\it t_5} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{{\it t_2}}}\\\\ 
C>   {\it t_7} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_8} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_9} &=& {{1}\over{{\it t_2}^{{\it t_8}}}}\\\\ 
C>   {\it t_{10}} &=& {{{\it t_8}}\over{2}}\\\\ 
C>   {\it t_{11}} &=& \sigma_{\alpha\alpha}^{{\it t_{10}}}\\\\ 
C>   {\it t_{12}} &=& {{1}\over{{\it t_1}^{8.0}}}\\\\ 
C>   {\it t_{13}} &=& {\it t_5}-.001890381166699926\\\\ 
C>   {\it t_{14}} &=& {{1}\over{{\it t_4}}}\\\\ 
C>   {\it t_{15}} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_{16}} &=& {{1}\over{{\it t_4}^{{\it t_8}}}}\\\\ 
C>   {\it t_{17}} &=& \sigma_{\beta\beta}^{{\it t_{10}}}\\\\ 
C>   {\it t_{18}} &=& {{1}\over{{\it t_3}^{8.0}}}\\\\ 
C>   {\it t_{19}} &=& 1.0\,\rho_s^{{{1}\over{3}}}\\\\ 
C>   {\it t_{20}} &=& {\it t_{19}}^{4.0}\\\\ 
C>   {\it t_{21}} &=& {{1}\over{{\it t_{20}}}}\\\\ 
C>   {\it t_{22}} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   {\it t_{23}} &=& {{1}\over{{\it t_{20}}^{{\it t_8}}}}\\\\ 
C>   {\it t_{24}} &=& \sigma_{ss}^{{\it t_{10}}}\\\\ 
C>   {\it t_{25}} &=& {{1}\over{{\it t_{19}}^{8.0}}}\\\\ 
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.0746613026776465 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, C. Fiolhais, Phys.Rev. B 46, 6671 (1992)  , DOI: 
C> <a href="http://dx.doi.org/10.1021/jp050536c "> 
C> 10.1021/jp050536c </a> 
C> 
      subroutine nwxcm_x_pw91_d3(param,tol_rho,ipol,nq,wght, 
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
