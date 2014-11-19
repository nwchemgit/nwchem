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
C>   {\it t_{13}} &=& {\it t_5}-0.001890381166699926\\\\ 
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
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
C>   {\it t_{13}} &=& {\it t_5}-0.001890381166699926\\\\ 
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
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
            t21 = (gammaaa*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammaaa*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammaaa*t19*t23*t8
            t30 = gammaaa**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammaaa*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammaaa**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammaaa*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = (2.0d+0*t13*t3*t38-2.0d+0*t17*t20*t3*t40)*wght
            t42 = 1/rhoa**1.6666666666666669d+0
            t43 = t2**2.0d+0
            t44 = 1/rhoa**1.3333333333333333d+0
            t45 = 1/t12**3
            t46 = 1/rhoa
            t47 = -1.4328817369035288d-6*t4*t46*t6*t7+t26+t24
            t48 = 1/rhoa**2
            t49 = t4**2
            t50 = 1/t2**1.0d+1
            t51 = 1/t2**1.8d+1
            t52 = -1.3333333333333333d-6*t4*t46*t6*t7+t33+t32+t29
            t53 = 1/t21**3
            t54 = (t20*(-2.0d+0*t17*t3*(1.3333333333333333d+1*t11*t44*t8
     1         *t9/t2**6.0d+0+5.333333333333333d+0*t11*t25*t42*t8*t9-1.0
     2         666666666666666d+1*t30*t44*t51*t53*t8+3.466666666666666d+
     3         1*gammaaa*t22*t44*t50*t8+5.333333333333333d+0*gammaaa*t22
     4         *t23*t42*t8+1.9105089825380384d-6*t27*t42*t49*t6*t7+1.432
     5         8817369035288d-6*t4*t48*t6*t7)-2.0d+0*t28*t3*t52-2.0d+0*t
     6         3*t34*t47)+2.0d+0*t13*t3*(-8.0d+0*gammaaa*t44*t50*t8-1.77
     7         77777777777776d+0*gammaaa*t23*t42*t8+1.7777777777777776d-
     8         6*t27*t42*t49*t6*t7+1.3333333333333333d-6*t4*t48*t6*t7-3.
     9         6566666666666664d+1*t15*t16*t30*t44*t51+8.0d+0*gammaaa*t1
     :         5*t16*t44*t50+1.9254543999999998d+1*gammaaa**3*t15*t16*t4
     ;         4/t2**2.6d+1-2.925333333333333d+0*t15*t16*t30*t31*t42+1.7
     <         777777777777776d+0*gammaaa*t15*t16*t23*t42)+t13*t19*(2.66
     =         66666666666666d+0*t18*t52+2.6666666666666666d+0*t18*t34)+
     >         t19*t20*(-2.6666666666666666d+0*t17*t18*t47-2.66666666666
     ?         66666d+0*t17*t18*t28)+4.0d+0*t17*t28*t3*t45*t47+2.6666666
     @         666666666d+0*t13*t17*t43*t44-2.4814019635975998d+0*t43*t4
     1         4-1.7777777777777776d+0*t13*t17*t18*t42+1.654267975731733
     2         d+0*t18*t42)*wght
            t55 = (t20*(-2.0d+0*t17*t3*(4.0d+0*gammaaa*t19*t31*t53*t8-4.
     1         0d+0*t11*t19*t25*t39*t8-1.2d+1*t19*t22*t23*t8-7.164408684
     2         517644d-7*t36*t46*t49*t7)-2.0d+0*t3*t40*t52-2.0d+0*t3*t38
     3         *t47)+2.0d+0*t13*t3*(2.6666666666666666d+0*t19*t23*t8-6.6
     4         66666666666666d-7*t36*t46*t49*t7+1.3164d+1*gammaaa*t15*t1
     5         6*t19*t31-7.220454d+0*t15*t16*t19*t30/t2**2.5d+1-2.666666
     6         6666666666d+0*t15*t16*t19*t23)+4.0d+0*t17*t3*t40*t45*t47-
     7         2.6666666666666666d+0*t17*t18*t19*t20*t40+2.6666666666666
     8         666d+0*t13*t18*t19*t38)*wght
            t56 = gammaaa**(t5-2)
            t57 = (t20*(-2.0d+0*t17*t3*(-1.5d+0*t10*t11*t8/t9**3-1.5d+0*
     1         t37*t53*t8+1.5d+0*t14*t22*t8/gammaaa+5.373306513388233d-7
     2         *t35*t4*t56*t7)-4.0d+0*t3*t38*t40)+2.0d+0*t13*t3*(5.0d-7*
     3         t35*t4*t56*t7-3.291d+0*t15*t16*t37+2.70767025d+0*gammaaa*
     4         t15*t16/t2**2.4d+1)+4.0d+0*t17*t3*t40**2*t45)*wght
            fnc(iq) = (2.0d+0*t13*t17*t3-1.8610514726982d+0*t3)*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = (2.0d+0*t13*t3*t34-2.0d+0*t17*t20*t28*t3+2.
     1         6666666666666666d+0*t13*t17*t18*t19-2.4814019635975998d+0
     2         *t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t41+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t41+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = t54+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t54+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t55+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = t55+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t55+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t57+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t57+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t57+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = t57+Cmat2(iq,D2_GAB_GAB)
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
            t34 = (gammaaa*t14+1)**5.0d-1
            t35 = 1/t34
            t36 = 1/t2**9.0d+0
            t37 = -8.0d+0*gammaaa*t32*t35*t36*t8
            t38 = 1/t2**5.0d+0
            t39 = -8.0d+0*t11*t32*t38*t8*t9
            t40 = 1.0d+0/t1
            t41 = -1.4328817369035288d-6*t32*t4*t40*t6*t7+t39+t37
            t42 = 2.6666666666666666d+0*gammaaa*t32*t36*t8
            t43 = gammaaa**2
            t44 = 1/t2**1.7d+1
            t45 = 4.3879999999999997d+0*t15*t16*t32*t43*t44
            t46 = -2.6666666666666666d+0*gammaaa*t15*t16*t32*t36
            t47 = -1.3333333333333333d-6*t32*t4*t40*t6*t7+t46+t45+t42
            t48 = t19**3.0d+0
            t49 = 1/rhob**6.666666666666666d-1
            t50 = 1/t26**2
            t51 = (gammabb*t28+1)**5.0d-1
            t52 = 1/t51
            t53 = 1/t19**9.0d+0
            t54 = -8.0d+0*gammabb*t49*t52*t53*t8
            t55 = 1/t19**5.0d+0
            t56 = -8.0d+0*t23*t25*t49*t55*t8
            t57 = 1.0d+0/t18
            t58 = -1.4328817369035288d-6*t21*t22*t4*t49*t57+t56+t54
            t59 = 2.6666666666666666d+0*gammabb*t49*t53*t8
            t60 = gammabb**2
            t61 = 1/t19**1.7d+1
            t62 = 4.3879999999999997d+0*t15*t29*t49*t60*t61
            t63 = -2.6666666666666666d+0*gammabb*t15*t29*t49*t53
            t64 = t63+t62+t59-1.3333333333333333d-6*t21*t22*t4*t49*t57
            t65 = t5-1
            t66 = gammaaa**t65
            t67 = 1/t2**1.6d+1
            t68 = -t14*t8+5.0d-7*t4*t66*t7-1.6455d+0*gammaaa*t15*t16*t67
     1         +t14*t15*t16
            t69 = 1/t9
            t70 = 3.0d+0*t10*t11*t69*t8+3.0d+0*t14*t35*t8+5.373306513388
     1         233d-7*t4*t66*t7
            t71 = gammabb**t65
            t72 = 1/t19**1.6d+1
            t73 = -t28*t8-1.6455d+0*gammabb*t15*t29*t72+5.0d-7*t22*t4*t7
     1         1+t15*t28*t29
            t74 = 1/t23
            t75 = 3.0d+0*t24*t25*t74*t8+3.0d+0*t28*t52*t8+5.373306513388
     1         233d-7*t22*t4*t71
            t76 = 1/rhoa**1.6666666666666669d+0
            t77 = t2**2.0d+0
            t78 = 1/rhoa**1.3333333333333333d+0
            t79 = 1/t12**3
            t80 = 1/rhoa
            t81 = -1.4328817369035288d-6*t4*t6*t7*t80+t39+t37
            t82 = 1/rhoa**2
            t83 = t4**2
            t84 = 1/t2**1.0d+1
            t85 = 1/t2**1.8d+1
            t86 = -1.3333333333333333d-6*t4*t6*t7*t80+t46+t45+t42
            t87 = 1/t34**3
            t88 = 1/rhob**1.6666666666666669d+0
            t89 = t19**2.0d+0
            t90 = 1/rhob**1.3333333333333333d+0
            t91 = 1/t26**3
            t92 = 1/rhob
            t93 = -1.4328817369035288d-6*t21*t22*t4*t92+t56+t54
            t94 = 1/rhob**2
            t95 = 1/t19**1.0d+1
            t96 = 1/t19**1.8d+1
            t97 = -1.3333333333333333d-6*t21*t22*t4*t92+t63+t62+t59
            t98 = 1/t51**3
            t99 = t5-2
            t100 = gammaaa**t99
            t101 = gammabb**t99
            fnc(iq) = (1.0d+0*t20*t27*t30+1.0d+0*t13*t17*t3-9.3052573634
     1         91d-1*t3-9.305257363491d-1*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t13*t3*t47-1.0d+0*t17*t3*t33*t41+1.
     1         3333333333333333d+0*t13*t17*t31*t32-1.2407009817987999d+0
     2         *t31*t32)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t20*t27*t64-1.0d+0*t20*t30*t50*t58+
     1         1.3333333333333333d+0*t27*t30*t48*t49-1.2407009817987999d
     2         +0*t48*t49)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*t68-1.0d+0*t17*t3*t33*t70)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t20*t27*t73-1.0d+0*t20*t30*t50*t75
     1         )*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t33*(-1.0d+0*t17*t3*(1.333333333333333
     1         3d+1*t11*t78*t8*t9/t2**6.0d+0+5.333333333333333d+0*t11*t3
     2         8*t76*t8*t9-1.0666666666666666d+1*t43*t78*t8*t85*t87+3.46
     3         6666666666666d+1*gammaaa*t35*t78*t8*t84+1.910508982538038
     4         4d-6*t40*t6*t7*t76*t83+1.4328817369035288d-6*t4*t6*t7*t82
     5         +5.333333333333333d+0*gammaaa*t35*t36*t76*t8)-1.0d+0*t3*t
     6         41*t86-1.0d+0*t3*t47*t81)+t13*t32*(1.3333333333333333d+0*
     7         t31*t86+1.3333333333333333d+0*t31*t47)+1.0d+0*t13*t3*(-3.
     8         6566666666666664d+1*t15*t16*t43*t78*t85-8.0d+0*gammaaa*t7
     9         8*t8*t84+8.0d+0*gammaaa*t15*t16*t78*t84+1.777777777777777
     :         6d-6*t40*t6*t7*t76*t83+1.3333333333333333d-6*t4*t6*t7*t82
     ;         -1.7777777777777776d+0*gammaaa*t36*t76*t8+1.9254543999999
     <         998d+1*gammaaa**3*t15*t16*t78/t2**2.6d+1-2.92533333333333
     =         3d+0*t15*t16*t43*t44*t76+1.7777777777777776d+0*gammaaa*t1
     >         5*t16*t36*t76)+t32*t33*(-1.3333333333333333d+0*t17*t31*t8
     ?         1-1.3333333333333333d+0*t17*t31*t41)+2.0d+0*t17*t3*t41*t7
     @         9*t81+1.3333333333333333d+0*t13*t17*t77*t78-1.24070098179
     1         87999d+0*t77*t78-8.888888888888888d-1*t13*t17*t31*t76+8.2
     2         71339878658666d-1*t31*t76)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t50*(-1.0d+0*t20*t30*(-1.0666666666666
     1         666d+1*t60*t8*t90*t96*t98+3.466666666666666d+1*gammabb*t5
     2         2*t8*t90*t95+1.4328817369035288d-6*t21*t22*t4*t94+1.33333
     3         33333333333d+1*t23*t25*t8*t90/t19**6.0d+0+1.9105089825380
     4         384d-6*t21*t22*t57*t83*t88+5.333333333333333d+0*t23*t25*t
     5         55*t8*t88+5.333333333333333d+0*gammabb*t52*t53*t8*t88)-1.
     6         0d+0*t20*t58*t97-1.0d+0*t20*t64*t93)+t27*t49*(1.333333333
     7         3333333d+0*t48*t97+1.3333333333333333d+0*t48*t64)+1.0d+0*
     8         t20*t27*(-3.6566666666666664d+1*t15*t29*t60*t90*t96-8.0d+
     9         0*gammabb*t8*t90*t95+8.0d+0*gammabb*t15*t29*t90*t95+1.333
     :         3333333333333d-6*t21*t22*t4*t94+1.9254543999999998d+1*gam
     ;         mabb**3*t15*t29*t90/t19**2.6d+1+1.7777777777777776d-6*t21
     <         *t22*t57*t83*t88-1.7777777777777776d+0*gammabb*t53*t8*t88
     =         -2.925333333333333d+0*t15*t29*t60*t61*t88+1.7777777777777
     >         776d+0*gammabb*t15*t29*t53*t88)+t49*t50*(-1.3333333333333
     ?         333d+0*t30*t48*t93-1.3333333333333333d+0*t30*t48*t58)+2.0
     @         d+0*t20*t30*t58*t91*t93+1.3333333333333333d+0*t27*t30*t89
     1         *t90-1.2407009817987999d+0*t89*t90-8.888888888888888d-1*t
     2         27*t30*t48*t88+8.271339878658666d-1*t48*t88)*wght+Amat2(i
     3         q,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t33*(-1.0d+0*t17*t3*(4.0d+0*gammaaa*t
     1         32*t44*t8*t87-7.164408684517644d-7*t66*t7*t80*t83-4.0d+0*
     2         t11*t32*t38*t69*t8-1.2d+1*t32*t35*t36*t8)-1.0d+0*t3*t70*t
     3         86-1.0d+0*t3*t68*t81)+1.0d+0*t13*t3*(-6.666666666666666d-
     4         7*t66*t7*t80*t83+2.6666666666666666d+0*t32*t36*t8+1.3164d
     5         +1*gammaaa*t15*t16*t32*t44-7.220454d+0*t15*t16*t32*t43/t2
     6         **2.5d+1-2.6666666666666666d+0*t15*t16*t32*t36)+2.0d+0*t1
     7         7*t3*t70*t79*t81-1.3333333333333333d+0*t17*t31*t32*t33*t7
     8         0+1.3333333333333333d+0*t13*t31*t32*t68)*wght+Cmat2(iq,D2
     9         _RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t50*(-1.0d+0*t20*t30*(4.0d+0*gammabb*
     1         t49*t61*t8*t98-7.164408684517644d-7*t22*t71*t83*t92-4.0d+
     2         0*t25*t49*t55*t74*t8-1.2d+1*t49*t52*t53*t8)-1.0d+0*t20*t7
     3         5*t97-1.0d+0*t20*t73*t93)+2.0d+0*t20*t30*t75*t91*t93+1.0d
     4         +0*t20*t27*(-6.666666666666666d-7*t22*t71*t83*t92+2.66666
     5         66666666666d+0*t49*t53*t8+1.3164d+1*gammabb*t15*t29*t49*t
     6         61-7.220454d+0*t15*t29*t49*t60/t19**2.5d+1-2.666666666666
     7         6666d+0*t15*t29*t49*t53)-1.3333333333333333d+0*t30*t48*t4
     8         9*t50*t75+1.3333333333333333d+0*t27*t48*t49*t73)*wght+Cma
     9         t2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (t33*(-1.0d+0*t17*t3*(-1.5d+0*t10*t11
     1         *t8/t9**3-1.5d+0*t67*t8*t87+1.5d+0*t14*t35*t8/gammaaa+5.3
     2         73306513388233d-7*t100*t4*t65*t7)-2.0d+0*t3*t68*t70)+2.0d
     3         +0*t17*t3*t70**2*t79+1.0d+0*t13*t3*(5.0d-7*t100*t4*t65*t7
     4         -3.291d+0*t15*t16*t67+2.70767025d+0*gammaaa*t15*t16/t2**2
     5         .4d+1))*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t50*(-1.0d+0*t20*t30*(-1.5d+0*t72*t8
     1         *t98+1.5d+0*t28*t52*t8/gammabb-1.5d+0*t24*t25*t8/t23**3+5
     2         .373306513388233d-7*t101*t22*t4*t65)-2.0d+0*t20*t73*t75)+
     3         2.0d+0*t20*t30*t75**2*t91+1.0d+0*t20*t27*(-3.291d+0*t15*t
     4         29*t72+5.0d-7*t101*t22*t4*t65+2.70767025d+0*gammabb*t15*t
     5         29/t19**2.4d+1))*wght+Cmat2(iq,D2_GBB_GBB)
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
            t21 = (gammaaa*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammaaa*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammaaa*t19*t23*t8
            t30 = gammaaa**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammaaa*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammaaa**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammaaa*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = 1/rhoa**1.6666666666666669d+0
            t42 = t2**2.0d+0
            t43 = 1/rhoa**1.3333333333333333d+0
            t44 = 1/t12**3
            t45 = 1/rhoa
            t46 = -1.4328817369035288d-6*t4*t45*t6*t7+t26+t24
            t47 = 1/rhoa**2
            t48 = t4**2
            t49 = 1/t2**1.0d+1
            t50 = 1/t2**1.8d+1
            t51 = -1.3333333333333333d-6*t4*t45*t6*t7+t33+t32+t29
            t52 = 1/t21**3
            t53 = gammaaa**(t5-2)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RA) = (1.0d+0*t13*t3*t34-1.0d+0*t17*t20*t28*t3+1.
     1         3333333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0
     2         *t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*t38-1.0d+0*t17*t20*t3*t40)*
     1         wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = (t20*(-1.0d+0*t17*t3*(1.333333333333333
     1         3d+1*t11*t43*t8*t9/t2**6.0d+0+5.333333333333333d+0*t11*t2
     2         5*t41*t8*t9-1.0666666666666666d+1*t30*t43*t50*t52*t8+3.46
     3         6666666666666d+1*gammaaa*t22*t43*t49*t8+5.333333333333333
     4         d+0*gammaaa*t22*t23*t41*t8+1.9105089825380384d-6*t27*t41*
     5         t48*t6*t7+1.4328817369035288d-6*t4*t47*t6*t7)-1.0d+0*t28*
     6         t3*t51-1.0d+0*t3*t34*t46)+1.0d+0*t13*t3*(-8.0d+0*gammaaa*
     7         t43*t49*t8-1.7777777777777776d+0*gammaaa*t23*t41*t8+1.777
     8         7777777777776d-6*t27*t41*t48*t6*t7+1.3333333333333333d-6*
     9         t4*t47*t6*t7-3.6566666666666664d+1*t15*t16*t30*t43*t50+8.
     :         0d+0*gammaaa*t15*t16*t43*t49+1.9254543999999998d+1*gammaa
     ;         a**3*t15*t16*t43/t2**2.6d+1-2.925333333333333d+0*t15*t16*
     <         t30*t31*t41+1.7777777777777776d+0*gammaaa*t15*t16*t23*t41
     =         )+t13*t19*(1.3333333333333333d+0*t18*t51+1.33333333333333
     >         33d+0*t18*t34)+t19*t20*(-1.3333333333333333d+0*t17*t18*t4
     ?         6-1.3333333333333333d+0*t17*t18*t28)+2.0d+0*t17*t28*t3*t4
     @         4*t46+1.3333333333333333d+0*t13*t17*t42*t43-1.24070098179
     1         87999d+0*t42*t43-8.888888888888888d-1*t13*t17*t18*t41+8.2
     2         71339878658666d-1*t18*t41)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (t20*(-1.0d+0*t17*t3*(4.0d+0*gammaaa*t
     1         19*t31*t52*t8-4.0d+0*t11*t19*t25*t39*t8-1.2d+1*t19*t22*t2
     2         3*t8-7.164408684517644d-7*t36*t45*t48*t7)-1.0d+0*t3*t40*t
     3         51-1.0d+0*t3*t38*t46)+1.0d+0*t13*t3*(2.6666666666666666d+
     4         0*t19*t23*t8-6.666666666666666d-7*t36*t45*t48*t7+1.3164d+
     5         1*gammaaa*t15*t16*t19*t31-7.220454d+0*t15*t16*t19*t30/t2*
     6         *2.5d+1-2.6666666666666666d+0*t15*t16*t19*t23)+2.0d+0*t17
     7         *t3*t40*t44*t46-1.3333333333333333d+0*t17*t18*t19*t20*t40
     8         +1.3333333333333333d+0*t13*t18*t19*t38)*wght+Cmat2(iq,D2_
     9         RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = (t20*(-1.0d+0*t17*t3*(-1.5d+0*t10*t11
     1         *t8/t9**3-1.5d+0*t37*t52*t8+1.5d+0*t14*t22*t8/gammaaa+5.3
     2         73306513388233d-7*t35*t4*t53*t7)-2.0d+0*t3*t38*t40)+1.0d+
     3         0*t13*t3*(5.0d-7*t35*t4*t53*t7-3.291d+0*t15*t16*t37+2.707
     4         67025d+0*gammaaa*t15*t16/t2**2.4d+1)+2.0d+0*t17*t3*t40**2
     5         *t44)*wght+Cmat2(iq,D2_GAA_GAA)
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
            t21 = (gammabb*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammabb*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammabb*t19*t23*t8
            t30 = gammabb**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammabb*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammabb**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammabb*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = 1/rhob**1.6666666666666669d+0
            t42 = t2**2.0d+0
            t43 = 1/rhob**1.3333333333333333d+0
            t44 = 1/t12**3
            t45 = 1/rhob
            t46 = -1.4328817369035288d-6*t4*t45*t6*t7+t26+t24
            t47 = 1/rhob**2
            t48 = t4**2
            t49 = 1/t2**1.0d+1
            t50 = 1/t2**1.8d+1
            t51 = -1.3333333333333333d-6*t4*t45*t6*t7+t33+t32+t29
            t52 = 1/t21**3
            t53 = gammabb**(t5-2)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RB) = (1.0d+0*t13*t3*t34-1.0d+0*t17*t20*t28*t3+1.
     1         3333333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0
     2         *t18*t19)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = (1.0d+0*t13*t3*t38-1.0d+0*t17*t20*t3*t40)*
     1         wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = (t20*(-1.0d+0*t17*t3*(1.333333333333333
     1         3d+1*t11*t43*t8*t9/t2**6.0d+0+5.333333333333333d+0*t11*t2
     2         5*t41*t8*t9-1.0666666666666666d+1*t30*t43*t50*t52*t8+3.46
     3         6666666666666d+1*gammabb*t22*t43*t49*t8+5.333333333333333
     4         d+0*gammabb*t22*t23*t41*t8+1.9105089825380384d-6*t27*t41*
     5         t48*t6*t7+1.4328817369035288d-6*t4*t47*t6*t7)-1.0d+0*t28*
     6         t3*t51-1.0d+0*t3*t34*t46)+1.0d+0*t13*t3*(-8.0d+0*gammabb*
     7         t43*t49*t8-1.7777777777777776d+0*gammabb*t23*t41*t8+1.777
     8         7777777777776d-6*t27*t41*t48*t6*t7+1.3333333333333333d-6*
     9         t4*t47*t6*t7-3.6566666666666664d+1*t15*t16*t30*t43*t50+8.
     :         0d+0*gammabb*t15*t16*t43*t49+1.9254543999999998d+1*gammab
     ;         b**3*t15*t16*t43/t2**2.6d+1-2.925333333333333d+0*t15*t16*
     <         t30*t31*t41+1.7777777777777776d+0*gammabb*t15*t16*t23*t41
     =         )+t13*t19*(1.3333333333333333d+0*t18*t51+1.33333333333333
     >         33d+0*t18*t34)+t19*t20*(-1.3333333333333333d+0*t17*t18*t4
     ?         6-1.3333333333333333d+0*t17*t18*t28)+2.0d+0*t17*t28*t3*t4
     @         4*t46+1.3333333333333333d+0*t13*t17*t42*t43-1.24070098179
     1         87999d+0*t42*t43-8.888888888888888d-1*t13*t17*t18*t41+8.2
     2         71339878658666d-1*t18*t41)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (t20*(-1.0d+0*t17*t3*(4.0d+0*gammabb*t
     1         19*t31*t52*t8-4.0d+0*t11*t19*t25*t39*t8-1.2d+1*t19*t22*t2
     2         3*t8-7.164408684517644d-7*t36*t45*t48*t7)-1.0d+0*t3*t40*t
     3         51-1.0d+0*t3*t38*t46)+1.0d+0*t13*t3*(2.6666666666666666d+
     4         0*t19*t23*t8-6.666666666666666d-7*t36*t45*t48*t7+1.3164d+
     5         1*gammabb*t15*t16*t19*t31-7.220454d+0*t15*t16*t19*t30/t2*
     6         *2.5d+1-2.6666666666666666d+0*t15*t16*t19*t23)+2.0d+0*t17
     7         *t3*t40*t44*t46-1.3333333333333333d+0*t17*t18*t19*t20*t40
     8         +1.3333333333333333d+0*t13*t18*t19*t38)*wght+Cmat2(iq,D2_
     9         RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t20*(-1.0d+0*t17*t3*(-1.5d+0*t10*t11
     1         *t8/t9**3-1.5d+0*t37*t52*t8+1.5d+0*t14*t22*t8/gammabb+5.3
     2         73306513388233d-7*t35*t4*t53*t7)-2.0d+0*t3*t38*t40)+1.0d+
     3         0*t13*t3*(5.0d-7*t35*t4*t53*t7-3.291d+0*t15*t16*t37+2.707
     4         67025d+0*gammabb*t15*t16/t2**2.4d+1)+2.0d+0*t17*t3*t40**2
     5         *t44)*wght+Cmat2(iq,D2_GBB_GBB)
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
C>   {\it t_{13}} &=& {\it t_5}-0.001890381166699926\\\\ 
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
C>   f &=& {{1.0\,{\it t_4}\,\left({{{\it t_{13}}\,{\it t_{18}}\,\sigma_{\beta\beta}}\over{e^{1.6455\,{\it t_{18}}\,\sigma_{\beta\beta}}}}+1.0 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}-{\it t_5}\,{\it t_{18}}\,\sigma_{\beta\beta}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{16}}\,{\it t_{17}}+6.0\,{\it t_5}\,{\it t_{14}}\,{\rm asinh}\; \left({\it t_{14}}\,{\it t_{15}}\right)\,{\it t_{15}}+1.0}}+{{1.0\,{\it t_2}\,\left({{{\it t_{13}}\,{\it t_{12}}\,\sigma_{\alpha\alpha}}\over{e^{1.6455\,{\it t_{12}}\,\sigma_{\alpha\alpha}}}}+1.0 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}-{\it t_5}\,{\it t_{12}}\,\sigma_{\alpha\alpha}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_9}\,{\it t_{11}}+6.0\,{\it t_5}\,{\it t_6}\,{\rm asinh}\; \left({\it t_6}\,{\it t_7}\right)\,{\it t_7}+1.0}}-0.9305257363491\,{\it t_4}-0.9305257363491\,{\it t_2}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& {{1.0\,{\it t_{20}}\,\left({{{\it t_{13}}\,{\it t_{25}}\,\sigma_{ss}}\over{e^{1.6455\,{\it t_{25}}\,\sigma_{ss}}}}+1.0 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}-{\it t_5}\,{\it t_{25}}\,\sigma_{ss}\right)}\over{1.074661302677647 \times 10^{-6}\,{\it t_{23}}\,{\it t_{24}}+6.0\,{\it t_5}\,{\it t_{21}}\,{\rm asinh}\; \left({\it t_{21}}\,{\it t_{22}}\right)\,{\it t_{22}}+1.0}}-0.9305257363491\,{\it t_{20}}\\\\ 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
            t21 = (gammaaa*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammaaa*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammaaa*t19*t23*t8
            t30 = gammaaa**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammaaa*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammaaa**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammaaa*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = (2.0d+0*t13*t3*t38-2.0d+0*t17*t20*t3*t40)*wght
            t42 = 1/rhoa**1.6666666666666669d+0
            t43 = t2**2.0d+0
            t44 = 1/rhoa**1.3333333333333333d+0
            t45 = 1/t12**3
            t46 = 1/rhoa
            t47 = -1.4328817369035288d-6*t4*t46*t6*t7+t26+t24
            t48 = 1/rhoa**2
            t49 = 1.3333333333333333d-6*t4*t48*t6*t7
            t50 = -1.7777777777777776d+0*gammaaa*t23*t42*t8
            t51 = t4**2
            t52 = 1/t2**1.0d+1
            t53 = -8.0d+0*gammaaa*t44*t52*t8
            t54 = -2.925333333333333d+0*t15*t16*t30*t31*t42
            t55 = 1.7777777777777776d+0*gammaaa*t15*t16*t23*t42
            t56 = gammaaa**3
            t57 = 1/t2**2.6d+1
            t58 = 1.9254543999999998d+1*t15*t16*t44*t56*t57
            t59 = 1/t2**1.8d+1
            t60 = -3.6566666666666664d+1*t15*t16*t30*t44*t59
            t61 = 8.0d+0*gammaaa*t15*t16*t44*t52
            t62 = 1.7777777777777776d-6*t27*t42*t51*t6*t7+t61+t60+t58+t5
     1         5+t54+t53+t50+t49
            t63 = -1.3333333333333333d-6*t4*t46*t6*t7+t33+t32+t29
            t64 = 1.4328817369035288d-6*t4*t48*t6*t7
            t65 = 5.333333333333333d+0*gammaaa*t22*t23*t42*t8
            t66 = 5.333333333333333d+0*t11*t25*t42*t8*t9
            t67 = 1/t21**3
            t68 = -1.0666666666666666d+1*t30*t44*t59*t67*t8
            t69 = 3.466666666666666d+1*gammaaa*t22*t44*t52*t8
            t70 = 1/t2**6.0d+0
            t71 = 1.3333333333333333d+1*t11*t44*t70*t8*t9
            t72 = t71+1.9105089825380384d-6*t27*t42*t51*t6*t7+t69+t68+t6
     1         6+t65+t64
            t73 = (t20*(-2.0d+0*t17*t3*t72-2.0d+0*t28*t3*t63-2.0d+0*t3*t
     1         34*t47)+t13*t19*(2.6666666666666666d+0*t18*t63+2.66666666
     2         66666666d+0*t18*t34)+2.0d+0*t13*t3*t62+t19*t20*(-2.666666
     3         6666666666d+0*t17*t18*t47-2.6666666666666666d+0*t17*t18*t
     4         28)+4.0d+0*t17*t28*t3*t45*t47+2.6666666666666666d+0*t13*t
     5         17*t43*t44-2.4814019635975998d+0*t43*t44-1.77777777777777
     6         76d+0*t13*t17*t18*t42+1.654267975731733d+0*t18*t42)*wght
            t74 = 1/t2**2.5d+1
            t75 = 2.6666666666666666d+0*t19*t23*t8-7.220454d+0*t15*t16*t
     1         19*t30*t74-6.666666666666666d-7*t36*t46*t51*t7+1.3164d+1*
     2         gammaaa*t15*t16*t19*t31-2.6666666666666666d+0*t15*t16*t19
     3         *t23
            t76 = 4.0d+0*gammaaa*t19*t31*t67*t8-4.0d+0*t11*t19*t25*t39*t
     1         8-1.2d+1*t19*t22*t23*t8-7.164408684517644d-7*t36*t46*t51*
     2         t7
            t77 = -2.0d+0*t17*t3*t76-2.0d+0*t3*t40*t63-2.0d+0*t3*t38*t47
            t78 = (t20*t77+2.0d+0*t13*t3*t75+4.0d+0*t17*t3*t40*t45*t47-2
     1         .6666666666666666d+0*t17*t18*t19*t20*t40+2.66666666666666
     2         66d+0*t13*t18*t19*t38)*wght
            t79 = t5-2
            t80 = gammaaa**t79
            t81 = 1/t2**2.4d+1
            t82 = 2.70767025d+0*gammaaa*t15*t16*t81+5.0d-7*t35*t4*t7*t80
     1         -3.291d+0*t15*t16*t37
            t83 = t40**2
            t84 = 1/gammaaa
            t85 = 1/t9**3
            t86 = -1.5d+0*t10*t11*t8*t85+1.5d+0*t14*t22*t8*t84+5.3733065
     1         13388233d-7*t35*t4*t7*t80-1.5d+0*t37*t67*t8
            t87 = -2.0d+0*t17*t3*t86-4.0d+0*t3*t38*t40
            t88 = (t20*t87+4.0d+0*t17*t3*t45*t83+2.0d+0*t13*t3*t82)*wght
            t89 = 1/rhoa**2.6666666666666666d+0
            t90 = 1/rhoa**2.3333333333333334d+0
            t91 = 1/t12**4
            t92 = t47**2
            t93 = 1/rhoa**3
            t94 = t4**3
            t95 = 1/t2**1.1d+1
            t96 = 1/t2**2.7d+1
            t97 = 1/t2**1.9d+1
            t98 = 1.7777777777777776d-6*t48*t51*t6*t7+t61+t60+t58+t55+t5
     1         4+t53+t50+t49
            t99 = 1/t21**5
            t100 = t71+1.9105089825380384d-6*t48*t51*t6*t7+t69+t68+t66+t
     1         65+t64
            t101 = -4.0d+0*t3*t47*t63-2.0d+0*t100*t17*t3
            t102 = (t20*(-2.0d+0*t17*t3*(-4.2666666666666664d+1*t48*t56*
     1         t8*t96*t99+1.1022222222222221d+2*t30*t48*t67*t8*t97-1.333
     2         333333333333d+2*gammaaa*t22*t48*t8*t95-2.5473453100507176
     3         d-6*t27*t6*t7*t89*t94-3.8210179650760767d-6*t51*t6*t7*t93
     4         -2.8657634738070575d-6*t4*t6*t7*t93-2.6666666666666666d+1
     5         *t11*t70*t8*t9*t90+2.1333333333333332d+1*t30*t59*t67*t8*t
     6         90-6.933333333333332d+1*gammaaa*t22*t52*t8*t90-8.88888888
     7         888889d+0*t11*t25*t8*t89*t9-2.6666666666666666d+1*t11*t48
     8         *t8*t9/t2**7.0d+0-8.88888888888889d+0*gammaaa*t22*t23*t8*
     9         t89-1.9105089825380384d-6*t27*t51*t6*t7*t89)-2.0d+0*t28*t
     :         3*t98-4.0d+0*t3*t63*t72-5.333333333333333d+0*t18*t19*t47*
     ;         t63-4.0d+0*t3*t47*t62-2.0d+0*t100*t3*t34-2.66666666666666
     <         66d+0*t100*t17*t18*t19)+t13*t19*(2.6666666666666666d+0*t1
     =         8*t98+5.333333333333333d+0*t19*t43*t63+5.333333333333333d
     >         +0*t18*t62)+2.0d+0*t13*t3*(2.5450399999999995d+2*t15*t16*
     ?         t30*t48*t97-3.2732724799999996d+2*t15*t16*t48*t56*t96+2.6
     @         666666666666666d+1*gammaaa*t48*t8*t95-2.6666666666666666d
     1         +1*gammaaa*t15*t16*t48*t95-2.37037037037037d-6*t27*t6*t7*
     2         t89*t94-3.555555555555555d-6*t51*t6*t7*t93-2.666666666666
     3         6666d-6*t4*t6*t7*t93+1.6d+1*gammaaa*t52*t8*t90+7.31333333
     4         3333333d+1*t15*t16*t30*t59*t90-3.8509087999999997d+1*t15*
     5         t16*t56*t57*t90-1.6d+1*gammaaa*t15*t16*t52*t90+2.96296296
     6         2962963d+0*gammaaa*t23*t8*t89-1.7777777777777776d-6*t27*t
     7         51*t6*t7*t89+4.875555555555556d+0*t15*t16*t30*t31*t89-2.9
     8         62962962962963d+0*gammaaa*t15*t16*t23*t89+8.4488939071999
     9         99d+1*gammaaa**4*t15*t16*t48/t2**3.5d+1)+t45*(4.0d+0*t3*t
     :         34*t92+8.0d+0*t17*t3*t47*t72-2*t101*t28)+t19*t45*(5.33333
     ;         3333333333d+0*t17*t18*t92+1.0666666666666666d+1*t17*t18*t
     <         28*t47)-1.2d+1*t17*t28*t3*t91*t92-5.333333333333333d+0*t1
     =         3*t17*t43*t90+4.9628039271951996d+0*t43*t90+2.96296296296
     >         2963d+0*t13*t17*t18*t89-2.757113292886222d+0*t18*t89+t19*
     ?         t20*(-5.333333333333333d+0*t17*t18*t72-5.333333333333333d
     @         +0*t18*t28*t63-5.333333333333333d+0*t17*t19*t43*t47-5.333
     1         333333333333d+0*t18*t34*t47)-3.555555555555555d+0*t13*t18
     2         *t42*t63+3.555555555555555d+0*t17*t18*t20*t42*t47+2.66666
     3         66666666666d+0*t13*t34*t43*t44-2.6666666666666666d+0*t17*
     4         t20*t28*t43*t44-1.7777777777777776d+0*t13*t18*t34*t42+1.7
     5         777777777777776d+0*t17*t18*t20*t28*t42+1.7777777777777776
     6         d+0*t13*t17*t42-1.654267975731733d+0*t42)*wght
            t103 = (t20*(-2.0d+0*t17*t3*(1.6d+1*t30*t44*t57*t8*t99+9.552
     1         544912690191d-7*t36*t48*t7*t94+6.666666666666666d+0*t11*t
     2         39*t44*t70*t8-3.8666666666666666d+1*gammaaa*t44*t59*t67*t
     3         8-2.6666666666666666d+0*gammaaa*t31*t42*t67*t8+4.13333333
     4         3333333d+1*t22*t44*t52*t8+2.6666666666666666d+0*t11*t25*t
     5         39*t42*t8+8.0d+0*t22*t23*t42*t8+7.164408684517644d-7*t36*
     6         t48*t51*t7)-2.0d+0*t3*t40*t98-4.0d+0*t3*t63*t76-4.0d+0*t3
     7         *t47*t75-2.0d+0*t100*t3*t38)+2.0d+0*t13*t3*(8.88888888888
     8         8887d-7*t36*t48*t7*t94-8.0d+0*t44*t52*t8-1.77777777777777
     9         76d+0*t23*t42*t8+4.813636d+0*t15*t16*t30*t42*t74+6.666666
     :         666666666d-7*t36*t48*t51*t7-8.629733333333333d+1*gammaaa*
     ;         t15*t16*t44*t59+1.17934082d+2*t15*t16*t30*t44*t57-3.16833
     <         52152d+1*t15*t16*t44*t56/t2**3.4d+1+8.0d+0*t15*t16*t44*t5
     =         2-8.775999999999999d+0*gammaaa*t15*t16*t31*t42+1.77777777
     >         77777776d+0*t15*t16*t23*t42)+t45*(4.0d+0*t3*t38*t92+8.0d+
     ?         0*t17*t3*t47*t76-2*t101*t40)-1.2d+1*t17*t3*t40*t91*t92+t1
     @         9*t20*(-5.333333333333333d+0*t17*t18*t76-5.33333333333333
     1         3d+0*t18*t40*t63-5.333333333333333d+0*t18*t38*t47)+5.3333
     2         33333333333d+0*t13*t18*t19*t75+1.0666666666666666d+1*t17*
     3         t18*t19*t40*t45*t47-2.6666666666666666d+0*t17*t20*t40*t43
     4         *t44+2.6666666666666666d+0*t13*t38*t43*t44+1.777777777777
     5         7776d+0*t17*t18*t20*t40*t42-1.7777777777777776d+0*t13*t18
     6         *t38*t42)*wght
            t104 = (t20*(-2.0d+0*t17*t3*(-6.0d+0*gammaaa*t19*t74*t8*t99+
     1         2.0d+0*t11*t19*t25*t8*t85-2.0d+0*t19*t22*t23*t8*t84-7.164
     2         408684517644d-7*t35*t46*t51*t7*t80+1.0d+1*t19*t31*t67*t8)
     3         -2.0d+0*t3*t63*t86-2.0d+0*t3*t47*t82-4.0d+0*t3*t38*t76-4.
     4         0d+0*t3*t40*t75)-1.2d+1*t17*t3*t47*t83*t91+t45*(4.0d+0*t1
     5         7*t3*t47*t86-2*t40*t77+4.0d+0*t17*t3*t40*t76+4.0d+0*t3*t3
     6         8*t40*t47)+t19*t20*(-2.6666666666666666d+0*t17*t18*t86-5.
     7         333333333333333d+0*t18*t38*t40)+5.333333333333333d+0*t17*
     8         t18*t19*t45*t83+2.6666666666666666d+0*t13*t18*t19*t82+2.0
     9         d+0*t13*t3*(-6.666666666666666d-7*t35*t46*t51*t7*t80-3.61
     :         0227d+1*gammaaa*t15*t16*t19*t74+1.7552d+1*t15*t16*t19*t31
     ;         +1.1881257057d+1*t15*t16*t19*t30/t2**3.3d+1))*wght
            t105 = gammaaa**(t5-3)
            t106 = (t20*(-2.0d+0*t17*t3*(2.25d+0*t8*t81*t99+2.25d+0*t10*
     1         t11*t8/t9**5-7.5d-1*t37*t67*t8*t84-2.25d+0*t14*t22*t8/t30
     2         +5.373306513388233d-7*t105*t35*t4*t7*t79)-6.0d+0*t3*t38*t
     3         86-6.0d+0*t3*t40*t82)-1.2d+1*t17*t3*t40**3*t91+t45*(-2*t4
     4         0*t87+8.0d+0*t17*t3*t40*t86+4.0d+0*t3*t38*t83)+2.0d+0*t13
     5         *t3*(8.12301075d+0*t15*t16*t81+5.0d-7*t105*t35*t4*t7*t79-
     6         4.455471396375d+0*gammaaa*t15*t16/t2**3.2d+1))*wght
            fnc(iq) = (2.0d+0*t13*t17*t3-1.8610514726982d+0*t3)*wght+fnc
     1         (iq)
            Amat(iq,D1_RA) = (2.0d+0*t13*t3*t34-2.0d+0*t17*t20*t28*t3+2.
     1         6666666666666666d+0*t13*t17*t18*t19-2.4814019635975998d+0
     2         *t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t41+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t41+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = t73+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t73+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t78+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = t78+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t78+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t88+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t88+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t88+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = t88+Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = t102+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t102+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t103+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = t103+Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t103+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t103+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = t103+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = t104+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = t104+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t104+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = t104+Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = t104+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t104+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t106+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = t106+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t106+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = t106+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t106+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = t106+Cmat3(iq,D3_GAB_GAB_GAB)
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
            t34 = (gammaaa*t14+1)**5.0d-1
            t35 = 1/t34
            t36 = 1/t2**9.0d+0
            t37 = -8.0d+0*gammaaa*t32*t35*t36*t8
            t38 = 1/t2**5.0d+0
            t39 = -8.0d+0*t11*t32*t38*t8*t9
            t40 = 1.0d+0/t1
            t41 = -1.4328817369035288d-6*t32*t4*t40*t6*t7+t39+t37
            t42 = 2.6666666666666666d+0*gammaaa*t32*t36*t8
            t43 = gammaaa**2
            t44 = 1/t2**1.7d+1
            t45 = 4.3879999999999997d+0*t15*t16*t32*t43*t44
            t46 = -2.6666666666666666d+0*gammaaa*t15*t16*t32*t36
            t47 = -1.3333333333333333d-6*t32*t4*t40*t6*t7+t46+t45+t42
            t48 = t19**3.0d+0
            t49 = 1/rhob**6.666666666666666d-1
            t50 = 1/t26**2
            t51 = (gammabb*t28+1)**5.0d-1
            t52 = 1/t51
            t53 = 1/t19**9.0d+0
            t54 = -8.0d+0*gammabb*t49*t52*t53*t8
            t55 = 1/t19**5.0d+0
            t56 = -8.0d+0*t23*t25*t49*t55*t8
            t57 = 1.0d+0/t18
            t58 = -1.4328817369035288d-6*t21*t22*t4*t49*t57+t56+t54
            t59 = 2.6666666666666666d+0*gammabb*t49*t53*t8
            t60 = gammabb**2
            t61 = 1/t19**1.7d+1
            t62 = 4.3879999999999997d+0*t15*t29*t49*t60*t61
            t63 = -2.6666666666666666d+0*gammabb*t15*t29*t49*t53
            t64 = t63+t62+t59-1.3333333333333333d-6*t21*t22*t4*t49*t57
            t65 = t5-1
            t66 = gammaaa**t65
            t67 = 1/t2**1.6d+1
            t68 = -t14*t8+5.0d-7*t4*t66*t7-1.6455d+0*gammaaa*t15*t16*t67
     1         +t14*t15*t16
            t69 = 1/t9
            t70 = 3.0d+0*t10*t11*t69*t8+3.0d+0*t14*t35*t8+5.373306513388
     1         233d-7*t4*t66*t7
            t71 = gammabb**t65
            t72 = 1/t19**1.6d+1
            t73 = -t28*t8-1.6455d+0*gammabb*t15*t29*t72+5.0d-7*t22*t4*t7
     1         1+t15*t28*t29
            t74 = 1/t23
            t75 = 3.0d+0*t24*t25*t74*t8+3.0d+0*t28*t52*t8+5.373306513388
     1         233d-7*t22*t4*t71
            t76 = 1/rhoa**1.6666666666666669d+0
            t77 = t2**2.0d+0
            t78 = 1/rhoa**1.3333333333333333d+0
            t79 = 1/t12**3
            t80 = 1/rhoa
            t81 = -1.4328817369035288d-6*t4*t6*t7*t80+t39+t37
            t82 = 1/rhoa**2
            t83 = 1.3333333333333333d-6*t4*t6*t7*t82
            t84 = -1.7777777777777776d+0*gammaaa*t36*t76*t8
            t85 = t4**2
            t86 = 1/t2**1.0d+1
            t87 = -8.0d+0*gammaaa*t78*t8*t86
            t88 = -2.925333333333333d+0*t15*t16*t43*t44*t76
            t89 = 1.7777777777777776d+0*gammaaa*t15*t16*t36*t76
            t90 = gammaaa**3
            t91 = 1/t2**2.6d+1
            t92 = 1.9254543999999998d+1*t15*t16*t78*t90*t91
            t93 = 1/t2**1.8d+1
            t94 = -3.6566666666666664d+1*t15*t16*t43*t78*t93
            t95 = 8.0d+0*gammaaa*t15*t16*t78*t86
            t96 = t95+t94+t92+t89+t88+t87+1.7777777777777776d-6*t40*t6*t
     1         7*t76*t85+t84+t83
            t97 = -1.3333333333333333d-6*t4*t6*t7*t80+t46+t45+t42
            t98 = 1.4328817369035288d-6*t4*t6*t7*t82
            t99 = 5.333333333333333d+0*gammaaa*t35*t36*t76*t8
            t100 = 5.333333333333333d+0*t11*t38*t76*t8*t9
            t101 = 1/t34**3
            t102 = -1.0666666666666666d+1*t101*t43*t78*t8*t93
            t103 = 3.466666666666666d+1*gammaaa*t35*t78*t8*t86
            t104 = 1/t2**6.0d+0
            t105 = 1.3333333333333333d+1*t104*t11*t78*t8*t9
            t106 = t99+t98+1.9105089825380384d-6*t40*t6*t7*t76*t85+t105+
     1         t103+t102+t100
            t107 = 1/rhob**1.6666666666666669d+0
            t108 = t19**2.0d+0
            t109 = 1/rhob**1.3333333333333333d+0
            t110 = 1/t26**3
            t111 = 1/rhob
            t112 = t56+t54-1.4328817369035288d-6*t111*t21*t22*t4
            t113 = 1/rhob**2
            t114 = 1.3333333333333333d-6*t113*t21*t22*t4
            t115 = -1.7777777777777776d+0*gammabb*t107*t53*t8
            t116 = 1/t19**1.0d+1
            t117 = -8.0d+0*gammabb*t109*t116*t8
            t118 = -2.925333333333333d+0*t107*t15*t29*t60*t61
            t119 = 1.7777777777777776d+0*gammabb*t107*t15*t29*t53
            t120 = gammabb**3
            t121 = 1/t19**2.6d+1
            t122 = 1.9254543999999998d+1*t109*t120*t121*t15*t29
            t123 = 1/t19**1.8d+1
            t124 = -3.6566666666666664d+1*t109*t123*t15*t29*t60
            t125 = 8.0d+0*gammabb*t109*t116*t15*t29
            t126 = 1.7777777777777776d-6*t107*t21*t22*t57*t85+t125+t124+
     1         t122+t119+t118+t117+t115+t114
            t127 = t63+t62+t59-1.3333333333333333d-6*t111*t21*t22*t4
            t128 = 1.4328817369035288d-6*t113*t21*t22*t4
            t129 = 5.333333333333333d+0*gammabb*t107*t52*t53*t8
            t130 = 5.333333333333333d+0*t107*t23*t25*t55*t8
            t131 = 1/t51**3
            t132 = -1.0666666666666666d+1*t109*t123*t131*t60*t8
            t133 = 3.466666666666666d+1*gammabb*t109*t116*t52*t8
            t134 = 1/t19**6.0d+0
            t135 = 1.3333333333333333d+1*t109*t134*t23*t25*t8
            t136 = 1.9105089825380384d-6*t107*t21*t22*t57*t85+t135+t133+
     1         t132+t130+t129+t128
            t137 = 1/t2**2.5d+1
            t138 = -6.666666666666666d-7*t66*t7*t80*t85+2.66666666666666
     1         66d+0*t32*t36*t8+1.3164d+1*gammaaa*t15*t16*t32*t44-7.2204
     2         54d+0*t137*t15*t16*t32*t43-2.6666666666666666d+0*t15*t16*
     3         t32*t36
            t139 = -7.164408684517644d-7*t66*t7*t80*t85-4.0d+0*t11*t32*t
     1         38*t69*t8+4.0d+0*gammaaa*t101*t32*t44*t8-1.2d+1*t32*t35*t
     2         36*t8
            t140 = -1.0d+0*t3*t70*t97-1.0d+0*t3*t68*t81-1.0d+0*t139*t17*
     1         t3
            t141 = 1/t19**2.5d+1
            t142 = -6.666666666666666d-7*t111*t22*t71*t85+2.666666666666
     1         6666d+0*t49*t53*t8+1.3164d+1*gammabb*t15*t29*t49*t61-7.22
     2         0454d+0*t141*t15*t29*t49*t60-2.6666666666666666d+0*t15*t2
     3         9*t49*t53
            t143 = -7.164408684517644d-7*t111*t22*t71*t85-4.0d+0*t25*t49
     1         *t55*t74*t8+4.0d+0*gammabb*t131*t49*t61*t8-1.2d+1*t49*t52
     2         *t53*t8
            t144 = -1.0d+0*t127*t20*t75-1.0d+0*t112*t20*t73-1.0d+0*t143*
     1         t20*t30
            t145 = t5-2
            t146 = gammaaa**t145
            t147 = 1/t2**2.4d+1
            t148 = 5.0d-7*t146*t4*t65*t7-3.291d+0*t15*t16*t67+2.70767025
     1         d+0*gammaaa*t147*t15*t16
            t149 = t70**2
            t150 = 1/gammaaa
            t151 = 1/t9**3
            t152 = -1.5d+0*t101*t67*t8+1.5d+0*t14*t150*t35*t8-1.5d+0*t10
     1         *t11*t151*t8+5.373306513388233d-7*t146*t4*t65*t7
            t153 = -2.0d+0*t3*t68*t70-1.0d+0*t152*t17*t3
            t154 = gammabb**t145
            t155 = 1/t19**2.4d+1
            t156 = -3.291d+0*t15*t29*t72+5.0d-7*t154*t22*t4*t65+2.707670
     1         25d+0*gammabb*t15*t155*t29
            t157 = t75**2
            t158 = 1/gammabb
            t159 = 1/t23**3
            t160 = -1.5d+0*t131*t72*t8+1.5d+0*t158*t28*t52*t8-1.5d+0*t15
     1         9*t24*t25*t8+5.373306513388233d-7*t154*t22*t4*t65
            t161 = -2.0d+0*t20*t73*t75-1.0d+0*t160*t20*t30
            t162 = 1/rhoa**2.6666666666666666d+0
            t163 = 1/rhoa**2.3333333333333334d+0
            t164 = 1/t12**4
            t165 = t81**2
            t166 = 1/rhoa**3
            t167 = t4**3
            t168 = 1/t2**1.1d+1
            t169 = 1/t2**2.7d+1
            t170 = 1/t2**1.9d+1
            t171 = t95+t94+t92+t89+t88+t87+1.7777777777777776d-6*t6*t7*t
     1         82*t85+t84+t83
            t172 = 1/t34**5
            t173 = t99+t98+1.9105089825380384d-6*t6*t7*t82*t85+t105+t103
     1         +t102+t100
            t174 = -2.0d+0*t3*t81*t97-1.0d+0*t17*t173*t3
            t175 = 1/rhob**2.6666666666666666d+0
            t176 = 1/rhob**2.3333333333333334d+0
            t177 = 1/t26**4
            t178 = t112**2
            t179 = 1/rhob**3
            t180 = 1/t19**1.1d+1
            t181 = 1/t19**2.7d+1
            t182 = 1/t19**1.9d+1
            t183 = 1.7777777777777776d-6*t113*t21*t22*t85+t125+t124+t122
     1         +t119+t118+t117+t115+t114
            t184 = 1/t51**5
            t185 = 1.9105089825380384d-6*t113*t21*t22*t85+t135+t133+t132
     1         +t130+t129+t128
            t186 = -1.0d+0*t185*t20*t30-2.0d+0*t112*t127*t20
            t187 = t5-3
            t188 = gammaaa**t187
            t189 = gammabb**t187
            fnc(iq) = (1.0d+0*t20*t27*t30+1.0d+0*t13*t17*t3-9.3052573634
     1         91d-1*t3-9.305257363491d-1*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (1.0d+0*t13*t3*t47-1.0d+0*t17*t3*t33*t41+1.
     1         3333333333333333d+0*t13*t17*t31*t32-1.2407009817987999d+0
     2         *t31*t32)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (1.0d+0*t20*t27*t64-1.0d+0*t20*t30*t50*t58+
     1         1.3333333333333333d+0*t27*t30*t48*t49-1.2407009817987999d
     2         +0*t48*t49)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*t68-1.0d+0*t17*t3*t33*t70)*
     1         wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (1.0d+0*t20*t27*t73-1.0d+0*t20*t30*t50*t75
     1         )*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (t33*(-1.0d+0*t3*t41*t97-1.0d+0*t3*t47*
     1         t81-1.0d+0*t106*t17*t3)+t13*t32*(1.3333333333333333d+0*t3
     2         1*t97+1.3333333333333333d+0*t31*t47)+1.0d+0*t13*t3*t96+t3
     3         2*t33*(-1.3333333333333333d+0*t17*t31*t81-1.3333333333333
     4         333d+0*t17*t31*t41)+2.0d+0*t17*t3*t41*t79*t81+1.333333333
     5         3333333d+0*t13*t17*t77*t78-1.2407009817987999d+0*t77*t78-
     6         8.888888888888888d-1*t13*t17*t31*t76+8.271339878658666d-1
     7         *t31*t76)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t27*t49*(1.3333333333333333d+0*t48*t64
     1         +1.3333333333333333d+0*t127*t48)+t50*(-1.0d+0*t112*t20*t6
     2         4-1.0d+0*t127*t20*t58-1.0d+0*t136*t20*t30)+t49*t50*(-1.33
     3         33333333333333d+0*t30*t48*t58-1.3333333333333333d+0*t112*
     4         t30*t48)+2.0d+0*t110*t112*t20*t30*t58-8.888888888888888d-
     5         1*t107*t27*t30*t48+8.271339878658666d-1*t107*t48+1.333333
     6         3333333333d+0*t108*t109*t27*t30+1.0d+0*t126*t20*t27-1.240
     7         7009817987999d+0*t108*t109)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0d+0*t17*t3*t70*t79*t81-1.333333333
     1         3333333d+0*t17*t31*t32*t33*t70+1.3333333333333333d+0*t13*
     2         t31*t32*t68+t140*t33+1.0d+0*t13*t138*t3)*wght+Cmat2(iq,D2
     3         _RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (-1.3333333333333333d+0*t30*t48*t49*t5
     1         0*t75+2.0d+0*t110*t112*t20*t30*t75+1.3333333333333333d+0*
     2         t27*t48*t49*t73+t144*t50+1.0d+0*t142*t20*t27)*wght+Cmat2(
     3         iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (2.0d+0*t149*t17*t3*t79+t153*t33+1.0d
     1         +0*t13*t148*t3)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t161*t50+2.0d+0*t110*t157*t20*t30+1.
     1         0d+0*t156*t20*t27)*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (t33*(-2.6666666666666666d+0*t31*t32
     1         *t81*t97-2.0d+0*t106*t3*t97-2.0d+0*t3*t81*t96-1.0d+0*t17*
     2         t3*(2.1333333333333332d+1*t101*t163*t43*t8*t93-4.26666666
     3         66666664d+1*t169*t172*t8*t82*t90-2.6666666666666666d+1*t1
     4         1*t8*t82*t9/t2**7.0d+0-8.88888888888889d+0*t11*t162*t38*t
     5         8*t9-2.6666666666666666d+1*t104*t11*t163*t8*t9-6.93333333
     6         3333332d+1*gammaaa*t163*t35*t8*t86-1.9105089825380384d-6*
     7         t162*t40*t6*t7*t85-3.8210179650760767d-6*t166*t6*t7*t85+1
     8         .1022222222222221d+2*t101*t170*t43*t8*t82-1.3333333333333
     9         33d+2*gammaaa*t168*t35*t8*t82-8.88888888888889d+0*gammaaa
     :         *t162*t35*t36*t8-2.5473453100507176d-6*t162*t167*t40*t6*t
     ;         7-2.8657634738070575d-6*t166*t4*t6*t7)-1.0d+0*t173*t3*t47
     <         -1.0d+0*t171*t3*t41-1.3333333333333333d+0*t17*t173*t31*t3
     =         2)+t13*t32*(2.6666666666666666d+0*t32*t77*t97+2.666666666
     >         6666666d+0*t31*t96+1.3333333333333333d+0*t171*t31)+t32*t3
     ?         3*(-2.6666666666666666d+0*t31*t41*t97-2.6666666666666666d
     @         +0*t17*t32*t77*t81-2.6666666666666666d+0*t31*t47*t81-2.66
     1         66666666666666d+0*t106*t17*t31)-1.7777777777777776d+0*t13
     2         *t31*t76*t97+1.0d+0*t13*t3*(7.313333333333333d+1*t15*t16*
     3         t163*t43*t93-3.8509087999999997d+1*t15*t16*t163*t90*t91-3
     4         .2732724799999996d+2*t15*t16*t169*t82*t90+1.6d+1*gammaaa*
     5         t163*t8*t86-1.6d+1*gammaaa*t15*t16*t163*t86-1.77777777777
     6         77776d-6*t162*t40*t6*t7*t85-3.555555555555555d-6*t166*t6*
     7         t7*t85+2.6666666666666666d+1*gammaaa*t168*t8*t82+2.545039
     8         9999999995d+2*t15*t16*t170*t43*t82+8.448893907199999d+1*g
     9         ammaaa**4*t15*t16*t82/t2**3.5d+1-2.6666666666666666d+1*ga
     :         mmaaa*t15*t16*t168*t82+2.962962962962963d+0*gammaaa*t162*
     ;         t36*t8-2.37037037037037d-6*t162*t167*t40*t6*t7-2.66666666
     <         66666666d-6*t166*t4*t6*t7+4.875555555555556d+0*t15*t16*t1
     =         62*t43*t44-2.962962962962963d+0*gammaaa*t15*t16*t162*t36)
     >         +t32*t79*(5.333333333333333d+0*t17*t31*t41*t81+2.66666666
     ?         66666666d+0*t165*t17*t31)+t79*(4.0d+0*t106*t17*t3*t81+2.0
     @         d+0*t165*t3*t47-2*t174*t41)+1.7777777777777776d+0*t17*t31
     1         *t33*t76*t81+1.3333333333333333d+0*t13*t47*t77*t78-1.3333
     2         333333333333d+0*t17*t33*t41*t77*t78-2.6666666666666666d+0
     3         *t13*t163*t17*t77+2.4814019635975998d+0*t163*t77-8.888888
     4         888888888d-1*t13*t31*t47*t76+8.888888888888888d-1*t17*t31
     5         *t33*t41*t76+8.888888888888888d-1*t13*t17*t76-8.271339878
     6         658666d-1*t76-6.0d+0*t164*t165*t17*t3*t41+1.4814814814814
     7         814d+0*t13*t162*t17*t31-1.378556646443111d+0*t162*t31)*wg
     8         ht+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t50*(-1.0d+0*t20*t30*(-1.9105089825
     1         380384d-6*t175*t21*t22*t57*t85-3.8210179650760767d-6*t179
     2         *t21*t22*t85+1.1022222222222221d+2*t113*t131*t182*t60*t8+
     3         2.1333333333333332d+1*t123*t131*t176*t60*t8-8.88888888888
     4         889d+0*t175*t23*t25*t55*t8-8.88888888888889d+0*gammabb*t1
     5         75*t52*t53*t8-1.333333333333333d+2*gammabb*t113*t180*t52*
     6         t8-6.933333333333332d+1*gammabb*t116*t176*t52*t8-2.666666
     7         6666666666d+1*t113*t23*t25*t8/t19**7.0d+0-2.6666666666666
     8         666d+1*t134*t176*t23*t25*t8-4.2666666666666664d+1*t113*t1
     9         20*t181*t184*t8-2.5473453100507176d-6*t167*t175*t21*t22*t
     :         57-2.8657634738070575d-6*t179*t21*t22*t4)-1.0d+0*t185*t20
     ;         *t64-1.0d+0*t183*t20*t58-1.3333333333333333d+0*t185*t30*t
     <         48*t49-2.6666666666666666d+0*t112*t127*t48*t49-2.0d+0*t12
     =         7*t136*t20-2.0d+0*t112*t126*t20)+1.0d+0*t20*t27*(-1.77777
     >         77777777776d-6*t175*t21*t22*t57*t85-3.555555555555555d-6*
     ?         t179*t21*t22*t85+2.962962962962963d+0*gammabb*t175*t53*t8
     @         +2.6666666666666666d+1*gammabb*t113*t180*t8+1.6d+1*gammab
     1         b*t116*t176*t8+4.875555555555556d+0*t15*t175*t29*t60*t61+
     2         2.5450399999999995d+2*t113*t15*t182*t29*t60+7.31333333333
     3         3333d+1*t123*t15*t176*t29*t60-2.37037037037037d-6*t167*t1
     4         75*t21*t22*t57-2.962962962962963d+0*gammabb*t15*t175*t29*
     5         t53-2.6666666666666666d-6*t179*t21*t22*t4+8.4488939071999
     6         99d+1*gammabb**4*t113*t15*t29/t19**3.5d+1-3.2732724799999
     7         996d+2*t113*t120*t15*t181*t29-2.6666666666666666d+1*gamma
     8         bb*t113*t15*t180*t29-3.8509087999999997d+1*t120*t121*t15*
     9         t176*t29-1.6d+1*gammabb*t116*t15*t176*t29)+t49*t50*(-2.66
     :         66666666666666d+0*t112*t48*t64-2.6666666666666666d+0*t127
     ;         *t48*t58-2.6666666666666666d+0*t108*t112*t30*t49-2.666666
     <         6666666666d+0*t136*t30*t48)+t110*(2.0d+0*t178*t20*t64-2*t
     =         186*t58+4.0d+0*t112*t136*t20*t30)-8.888888888888888d-1*t1
     >         07*t27*t48*t64+1.3333333333333333d+0*t108*t109*t27*t64+t1
     ?         10*t49*(5.333333333333333d+0*t112*t30*t48*t58+2.666666666
     @         6666666d+0*t178*t30*t48)+8.888888888888888d-1*t107*t30*t4
     1         8*t50*t58-1.3333333333333333d+0*t108*t109*t30*t50*t58-6.0
     2         d+0*t177*t178*t20*t30*t58+1.7777777777777776d+0*t107*t112
     3         *t30*t48*t50+t27*t49*(2.6666666666666666d+0*t108*t127*t49
     4         +1.3333333333333333d+0*t183*t48+2.6666666666666666d+0*t12
     5         6*t48)+1.4814814814814814d+0*t175*t27*t30*t48-1.777777777
     6         7777776d+0*t107*t127*t27*t48-1.378556646443111d+0*t175*t4
     7         8-2.6666666666666666d+0*t108*t176*t27*t30+8.8888888888888
     8         88d-1*t107*t27*t30+2.4814019635975998d+0*t108*t176-8.2713
     9         39878658666d-1*t107)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t32*t33*(-2.6666666666666666d+0*t3
     1         1*t70*t97-2.6666666666666666d+0*t31*t68*t81-2.66666666666
     2         66666d+0*t139*t17*t31)+t33*(-2.0d+0*t139*t3*t97-1.0d+0*t1
     3         7*t3*(-3.8666666666666666d+1*gammaaa*t101*t78*t8*t93+1.6d
     4         +1*t172*t43*t78*t8*t91+4.133333333333333d+1*t35*t78*t8*t8
     5         6+7.164408684517644d-7*t66*t7*t82*t85+9.552544912690191d-
     6         7*t167*t66*t7*t82+6.666666666666666d+0*t104*t11*t69*t78*t
     7         8+2.6666666666666666d+0*t11*t38*t69*t76*t8-2.666666666666
     8         6666d+0*gammaaa*t101*t44*t76*t8+8.0d+0*t35*t36*t76*t8)-2.
     9         0d+0*t138*t3*t81-1.0d+0*t171*t3*t70-1.0d+0*t173*t3*t68)+1
     :         .0d+0*t13*t3*(-8.629733333333333d+1*gammaaa*t15*t16*t78*t
     ;         93+1.17934082d+2*t15*t16*t43*t78*t91-3.1683352152d+1*t15*
     <         t16*t78*t90/t2**3.4d+1-8.0d+0*t78*t8*t86+8.0d+0*t15*t16*t
     =         78*t86+6.666666666666666d-7*t66*t7*t82*t85+8.888888888888
     >         887d-7*t167*t66*t7*t82-1.7777777777777776d+0*t36*t76*t8-8
     ?         .775999999999999d+0*gammaaa*t15*t16*t44*t76+4.813636d+0*t
     @         137*t15*t16*t43*t76+1.7777777777777776d+0*t15*t16*t36*t76
     1         )+t79*(4.0d+0*t139*t17*t3*t81-2*t174*t70+2.0d+0*t165*t3*t
     2         68)+5.333333333333333d+0*t17*t31*t32*t70*t79*t81-1.333333
     3         3333333333d+0*t17*t33*t70*t77*t78+1.3333333333333333d+0*t
     4         13*t68*t77*t78+8.888888888888888d-1*t17*t31*t33*t70*t76-8
     5         .888888888888888d-1*t13*t31*t68*t76-6.0d+0*t164*t165*t17*
     6         t3*t70+2.6666666666666666d+0*t13*t138*t31*t32)*wght+Cmat3
     7         (iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t50*(-1.0d+0*t20*t30*(7.1644086845
     1         17644d-7*t113*t22*t71*t85+2.6666666666666666d+0*t107*t25*
     2         t55*t74*t8+6.666666666666666d+0*t109*t134*t25*t74*t8-2.66
     3         66666666666666d+0*gammabb*t107*t131*t61*t8+1.6d+1*t109*t1
     4         21*t184*t60*t8+8.0d+0*t107*t52*t53*t8+4.133333333333333d+
     5         1*t109*t116*t52*t8-3.8666666666666666d+1*gammabb*t109*t12
     6         3*t131*t8+9.552544912690191d-7*t113*t167*t22*t71)-1.0d+0*
     7         t183*t20*t75-1.0d+0*t185*t20*t73-2.0d+0*t127*t143*t20-2.0
     8         d+0*t112*t142*t20)+1.0d+0*t20*t27*(6.666666666666666d-7*t
     9         113*t22*t71*t85-1.7777777777777776d+0*t107*t53*t8-8.0d+0*
     :         t109*t116*t8+8.888888888888887d-7*t113*t167*t22*t71-8.775
     ;         999999999999d+0*gammabb*t107*t15*t29*t61+4.813636d+0*t107
     <         *t141*t15*t29*t60+1.17934082d+2*t109*t121*t15*t29*t60+1.7
     =         777777777777776d+0*t107*t15*t29*t53-3.1683352152d+1*t109*
     >         t120*t15*t29/t19**3.4d+1-8.629733333333333d+1*gammabb*t10
     ?         9*t123*t15*t29+8.0d+0*t109*t116*t15*t29)+t49*t50*(-2.6666
     @         666666666666d+0*t127*t48*t75-2.6666666666666666d+0*t112*t
     1         48*t73-2.6666666666666666d+0*t143*t30*t48)+t110*(-2*t186*
     2         t75+2.0d+0*t178*t20*t73+4.0d+0*t112*t143*t20*t30)+8.88888
     3         8888888888d-1*t107*t30*t48*t50*t75-1.3333333333333333d+0*
     4         t108*t109*t30*t50*t75+5.333333333333333d+0*t110*t112*t30*
     5         t48*t49*t75-6.0d+0*t177*t178*t20*t30*t75-8.88888888888888
     6         8d-1*t107*t27*t48*t73+1.3333333333333333d+0*t108*t109*t27
     7         *t73+2.6666666666666666d+0*t142*t27*t48*t49)*wght+Cmat3(i
     8         q,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t33*(-1.0d+0*t152*t3*t97-1.0d+0*t
     1         17*t3*(-7.164408684517644d-7*t146*t65*t7*t80*t85+1.0d+1*t
     2         101*t32*t44*t8+2.0d+0*t11*t151*t32*t38*t8-2.0d+0*t150*t32
     3         *t35*t36*t8-6.0d+0*gammaaa*t137*t172*t32*t8)-1.0d+0*t148*
     4         t3*t81-2.0d+0*t138*t3*t70-2.0d+0*t139*t3*t68)+1.0d+0*t13*
     5         t3*(-6.666666666666666d-7*t146*t65*t7*t80*t85+1.7552d+1*t
     6         15*t16*t32*t44+1.1881257057d+1*t15*t16*t32*t43/t2**3.3d+1
     7         -3.610227d+1*gammaaa*t137*t15*t16*t32)+t79*(2.0d+0*t3*t68
     8         *t70*t81+2.0d+0*t152*t17*t3*t81+2.0d+0*t139*t17*t3*t70-2*
     9         t140*t70)-6.0d+0*t149*t164*t17*t3*t81+2.6666666666666666d
     :         +0*t149*t17*t31*t32*t79+t32*t33*(-2.6666666666666666d+0*t
     ;         31*t68*t70-1.3333333333333333d+0*t152*t17*t31)+1.33333333
     <         33333333d+0*t13*t148*t31*t32)*wght+Cmat3(iq,D3_RA_GAA_GAA
     =         )
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t50*(-1.0d+0*t20*t30*(-7.16440868
     1         4517644d-7*t111*t154*t22*t65*t85+1.0d+1*t131*t49*t61*t8+2
     2         .0d+0*t159*t25*t49*t55*t8-2.0d+0*t158*t49*t52*t53*t8-6.0d
     3         +0*gammabb*t141*t184*t49*t8)-2.0d+0*t142*t20*t75-2.0d+0*t
     4         143*t20*t73-1.0d+0*t127*t160*t20-1.0d+0*t112*t156*t20)+1.
     5         0d+0*t20*t27*(-6.666666666666666d-7*t111*t154*t22*t65*t85
     6         +1.7552d+1*t15*t29*t49*t61+1.1881257057d+1*t15*t29*t49*t6
     7         0/t19**3.3d+1-3.610227d+1*gammabb*t141*t15*t29*t49)+t49*t
     8         50*(-2.6666666666666666d+0*t48*t73*t75-1.3333333333333333
     9         d+0*t160*t30*t48)+t110*(2.0d+0*t112*t20*t73*t75+2.0d+0*t1
     :         43*t20*t30*t75-2*t144*t75+2.0d+0*t112*t160*t20*t30)+2.666
     ;         6666666666666d+0*t110*t157*t30*t48*t49+1.3333333333333333
     <         d+0*t156*t27*t48*t49-6.0d+0*t112*t157*t177*t20*t30)*wght+
     =         Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t33*(-1.0d+0*t17*t3*(2.25d+0*t10
     1         *t11*t8/t9**5-7.5d-1*t101*t150*t67*t8-2.25d+0*t14*t35*t8/
     2         t43+2.25d+0*t147*t172*t8+5.373306513388233d-7*t145*t188*t
     3         4*t65*t7)-3.0d+0*t148*t3*t70-3.0d+0*t152*t3*t68)+(4.0d+0*
     4         t152*t17*t3*t70-2*t153*t70+2.0d+0*t149*t3*t68)*t79-6.0d+0
     5         *t164*t17*t3*t70**3+1.0d+0*t13*t3*(5.0d-7*t145*t188*t4*t6
     6         5*t7-4.455471396375d+0*gammaaa*t15*t16/t2**3.2d+1+8.12301
     7         075d+0*t147*t15*t16))*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (t50*(-1.0d+0*t20*t30*(-7.5d-1*t1
     1         31*t158*t72*t8-2.25d+0*t28*t52*t8/t60+2.25d+0*t24*t25*t8/
     2         t23**5+2.25d+0*t155*t184*t8+5.373306513388233d-7*t145*t18
     3         9*t22*t4*t65)-3.0d+0*t156*t20*t75-3.0d+0*t160*t20*t73)-6.
     4         0d+0*t177*t20*t30*t75**3+t110*(4.0d+0*t160*t20*t30*t75-2*
     5         t161*t75+2.0d+0*t157*t20*t73)+1.0d+0*t20*t27*(5.0d-7*t145
     6         *t189*t22*t4*t65-4.455471396375d+0*gammabb*t15*t29/t19**3
     7         .2d+1+8.12301075d+0*t15*t155*t29))*wght+Cmat3(iq,D3_GBB_G
     8         BB_GBB)
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
            t21 = (gammaaa*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammaaa*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammaaa*t19*t23*t8
            t30 = gammaaa**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammaaa*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammaaa**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammaaa*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = 1/rhoa**1.6666666666666669d+0
            t42 = t2**2.0d+0
            t43 = 1/rhoa**1.3333333333333333d+0
            t44 = 1/t12**3
            t45 = 1/rhoa
            t46 = -1.4328817369035288d-6*t4*t45*t6*t7+t26+t24
            t47 = 1/rhoa**2
            t48 = 1.3333333333333333d-6*t4*t47*t6*t7
            t49 = -1.7777777777777776d+0*gammaaa*t23*t41*t8
            t50 = t4**2
            t51 = 1/t2**1.0d+1
            t52 = -8.0d+0*gammaaa*t43*t51*t8
            t53 = -2.925333333333333d+0*t15*t16*t30*t31*t41
            t54 = 1.7777777777777776d+0*gammaaa*t15*t16*t23*t41
            t55 = gammaaa**3
            t56 = 1/t2**2.6d+1
            t57 = 1.9254543999999998d+1*t15*t16*t43*t55*t56
            t58 = 1/t2**1.8d+1
            t59 = -3.6566666666666664d+1*t15*t16*t30*t43*t58
            t60 = 8.0d+0*gammaaa*t15*t16*t43*t51
            t61 = 1.7777777777777776d-6*t27*t41*t50*t6*t7+t60+t59+t57+t5
     1         4+t53+t52+t49+t48
            t62 = -1.3333333333333333d-6*t4*t45*t6*t7+t33+t32+t29
            t63 = 1.4328817369035288d-6*t4*t47*t6*t7
            t64 = 5.333333333333333d+0*gammaaa*t22*t23*t41*t8
            t65 = 5.333333333333333d+0*t11*t25*t41*t8*t9
            t66 = 1/t21**3
            t67 = -1.0666666666666666d+1*t30*t43*t58*t66*t8
            t68 = 3.466666666666666d+1*gammaaa*t22*t43*t51*t8
            t69 = 1/t2**6.0d+0
            t70 = 1.3333333333333333d+1*t11*t43*t69*t8*t9
            t71 = t70+1.9105089825380384d-6*t27*t41*t50*t6*t7+t68+t67+t6
     1         5+t64+t63
            t72 = 1/t2**2.5d+1
            t73 = 2.6666666666666666d+0*t19*t23*t8-7.220454d+0*t15*t16*t
     1         19*t30*t72-6.666666666666666d-7*t36*t45*t50*t7+1.3164d+1*
     2         gammaaa*t15*t16*t19*t31-2.6666666666666666d+0*t15*t16*t19
     3         *t23
            t74 = 4.0d+0*gammaaa*t19*t31*t66*t8-4.0d+0*t11*t19*t25*t39*t
     1         8-1.2d+1*t19*t22*t23*t8-7.164408684517644d-7*t36*t45*t50*
     2         t7
            t75 = -1.0d+0*t17*t3*t74-1.0d+0*t3*t40*t62-1.0d+0*t3*t38*t46
            t76 = t5-2
            t77 = gammaaa**t76
            t78 = 1/t2**2.4d+1
            t79 = 2.70767025d+0*gammaaa*t15*t16*t78+5.0d-7*t35*t4*t7*t77
     1         -3.291d+0*t15*t16*t37
            t80 = t40**2
            t81 = 1/gammaaa
            t82 = 1/t9**3
            t83 = -1.5d+0*t10*t11*t8*t82+1.5d+0*t14*t22*t8*t81-1.5d+0*t3
     1         7*t66*t8+5.373306513388233d-7*t35*t4*t7*t77
            t84 = -1.0d+0*t17*t3*t83-2.0d+0*t3*t38*t40
            t85 = 1/rhoa**2.6666666666666666d+0
            t86 = 1/rhoa**2.3333333333333334d+0
            t87 = 1/t12**4
            t88 = t46**2
            t89 = 1/rhoa**3
            t90 = t4**3
            t91 = 1/t2**1.1d+1
            t92 = 1/t2**2.7d+1
            t93 = 1/t2**1.9d+1
            t94 = 1.7777777777777776d-6*t47*t50*t6*t7+t60+t59+t57+t54+t5
     1         3+t52+t49+t48
            t95 = 1/t21**5
            t96 = t70+1.9105089825380384d-6*t47*t50*t6*t7+t68+t67+t65+t6
     1         4+t63
            t97 = -1.0d+0*t17*t3*t96-2.0d+0*t3*t46*t62
            t98 = gammaaa**(t5-3)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RA) = (1.0d+0*t13*t3*t34-1.0d+0*t17*t20*t28*t3+1.
     1         3333333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0
     2         *t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (1.0d+0*t13*t3*t38-1.0d+0*t17*t20*t3*t40)*
     1         wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = (t20*(-1.0d+0*t17*t3*t71-1.0d+0*t28*t3*
     1         t62-1.0d+0*t3*t34*t46)+t13*t19*(1.3333333333333333d+0*t18
     2         *t62+1.3333333333333333d+0*t18*t34)+1.0d+0*t13*t3*t61+t19
     3         *t20*(-1.3333333333333333d+0*t17*t18*t46-1.33333333333333
     4         33d+0*t17*t18*t28)+2.0d+0*t17*t28*t3*t44*t46+1.3333333333
     5         333333d+0*t13*t17*t42*t43-1.2407009817987999d+0*t42*t43-8
     6         .888888888888888d-1*t13*t17*t18*t41+8.271339878658666d-1*
     7         t18*t41)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (t20*t75+1.0d+0*t13*t3*t73+2.0d+0*t17*
     1         t3*t40*t44*t46-1.3333333333333333d+0*t17*t18*t19*t20*t40+
     2         1.3333333333333333d+0*t13*t18*t19*t38)*wght+Cmat2(iq,D2_R
     3         A_GAA)
            Cmat2(iq,D2_GAA_GAA) = (t20*t84+2.0d+0*t17*t3*t44*t80+1.0d+0
     1         *t13*t3*t79)*wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = (t44*(-2*t28*t97+2.0d+0*t3*t34*t88+4
     1         .0d+0*t17*t3*t46*t71)+t20*(-1.0d+0*t3*t34*t96-1.333333333
     2         3333333d+0*t17*t18*t19*t96-1.0d+0*t17*t3*(-4.266666666666
     3         6664d+1*t47*t55*t8*t92*t95+1.1022222222222221d+2*t30*t47*
     4         t66*t8*t93-1.333333333333333d+2*gammaaa*t22*t47*t8*t91-2.
     5         5473453100507176d-6*t27*t6*t7*t85*t90-2.6666666666666666d
     6         +1*t11*t69*t8*t86*t9-8.88888888888889d+0*t11*t25*t8*t85*t
     7         9-2.6666666666666666d+1*t11*t47*t8*t9/t2**7.0d+0-3.821017
     8         9650760767d-6*t50*t6*t7*t89-2.8657634738070575d-6*t4*t6*t
     9         7*t89+2.1333333333333332d+1*t30*t58*t66*t8*t86-6.93333333
     :         3333332d+1*gammaaa*t22*t51*t8*t86-8.88888888888889d+0*gam
     ;         maaa*t22*t23*t8*t85-1.9105089825380384d-6*t27*t50*t6*t7*t
     <         85)-1.0d+0*t28*t3*t94-2.0d+0*t3*t62*t71-2.666666666666666
     =         6d+0*t18*t19*t46*t62-2.0d+0*t3*t46*t61)+t13*t19*(1.333333
     >         3333333333d+0*t18*t94+2.6666666666666666d+0*t19*t42*t62+2
     ?         .6666666666666666d+0*t18*t61)+1.0d+0*t13*t3*(2.5450399999
     @         999995d+2*t15*t16*t30*t47*t93-3.2732724799999996d+2*t15*t
     1         16*t47*t55*t92+2.6666666666666666d+1*gammaaa*t47*t8*t91-2
     2         .6666666666666666d+1*gammaaa*t15*t16*t47*t91-2.3703703703
     3         7037d-6*t27*t6*t7*t85*t90-3.555555555555555d-6*t50*t6*t7*
     4         t89-2.6666666666666666d-6*t4*t6*t7*t89+1.6d+1*gammaaa*t51
     5         *t8*t86+7.313333333333333d+1*t15*t16*t30*t58*t86-3.850908
     6         7999999997d+1*t15*t16*t55*t56*t86-1.6d+1*gammaaa*t15*t16*
     7         t51*t86+2.962962962962963d+0*gammaaa*t23*t8*t85-1.7777777
     8         777777776d-6*t27*t50*t6*t7*t85+4.875555555555556d+0*t15*t
     9         16*t30*t31*t85-2.962962962962963d+0*gammaaa*t15*t16*t23*t
     :         85+8.448893907199999d+1*gammaaa**4*t15*t16*t47/t2**3.5d+1
     ;         )+t19*t44*(2.6666666666666666d+0*t17*t18*t88+5.3333333333
     <         33333d+0*t17*t18*t28*t46)-6.0d+0*t17*t28*t3*t87*t88-2.666
     =         6666666666666d+0*t13*t17*t42*t86+2.4814019635975998d+0*t4
     >         2*t86+1.4814814814814814d+0*t13*t17*t18*t85-1.37855664644
     ?         3111d+0*t18*t85+t19*t20*(-2.6666666666666666d+0*t17*t18*t
     @         71-2.6666666666666666d+0*t18*t28*t62-2.6666666666666666d+
     1         0*t17*t19*t42*t46-2.6666666666666666d+0*t18*t34*t46)-1.77
     2         77777777777776d+0*t13*t18*t41*t62+1.7777777777777776d+0*t
     3         17*t18*t20*t41*t46+1.3333333333333333d+0*t13*t34*t42*t43-
     4         1.3333333333333333d+0*t17*t20*t28*t42*t43-8.8888888888888
     5         88d-1*t13*t18*t34*t41+8.888888888888888d-1*t17*t18*t20*t2
     6         8*t41+8.888888888888888d-1*t13*t17*t41-8.271339878658666d
     7         -1*t41)*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = (t44*(-2*t40*t97+2.0d+0*t3*t38*t88+
     1         4.0d+0*t17*t3*t46*t74)+t20*(-1.0d+0*t3*t38*t96-1.0d+0*t17
     2         *t3*(1.6d+1*t30*t43*t56*t8*t95+9.552544912690191d-7*t36*t
     3         47*t7*t90+6.666666666666666d+0*t11*t39*t43*t69*t8-3.86666
     4         66666666666d+1*gammaaa*t43*t58*t66*t8-2.6666666666666666d
     5         +0*gammaaa*t31*t41*t66*t8+4.133333333333333d+1*t22*t43*t5
     6         1*t8+2.6666666666666666d+0*t11*t25*t39*t41*t8+8.0d+0*t22*
     7         t23*t41*t8+7.164408684517644d-7*t36*t47*t50*t7)-1.0d+0*t3
     8         *t40*t94-2.0d+0*t3*t62*t74-2.0d+0*t3*t46*t73)+1.0d+0*t13*
     9         t3*(8.888888888888887d-7*t36*t47*t7*t90-8.0d+0*t43*t51*t8
     :         -1.7777777777777776d+0*t23*t41*t8+4.813636d+0*t15*t16*t30
     ;         *t41*t72+6.666666666666666d-7*t36*t47*t50*t7-8.6297333333
     <         33333d+1*gammaaa*t15*t16*t43*t58+1.17934082d+2*t15*t16*t3
     =         0*t43*t56-3.1683352152d+1*t15*t16*t43*t55/t2**3.4d+1+8.0d
     >         +0*t15*t16*t43*t51-8.775999999999999d+0*gammaaa*t15*t16*t
     ?         31*t41+1.7777777777777776d+0*t15*t16*t23*t41)-6.0d+0*t17*
     @         t3*t40*t87*t88+t19*t20*(-2.6666666666666666d+0*t17*t18*t7
     1         4-2.6666666666666666d+0*t18*t40*t62-2.6666666666666666d+0
     2         *t18*t38*t46)+2.6666666666666666d+0*t13*t18*t19*t73+5.333
     3         333333333333d+0*t17*t18*t19*t40*t44*t46-1.333333333333333
     4         3d+0*t17*t20*t40*t42*t43+1.3333333333333333d+0*t13*t38*t4
     5         2*t43+8.888888888888888d-1*t17*t18*t20*t40*t41-8.88888888
     6         8888888d-1*t13*t18*t38*t41)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = (t20*(-1.0d+0*t17*t3*(-6.0d+0*gamm
     1         aaa*t19*t72*t8*t95+2.0d+0*t11*t19*t25*t8*t82-2.0d+0*t19*t
     2         22*t23*t8*t81+1.0d+1*t19*t31*t66*t8-7.164408684517644d-7*
     3         t35*t45*t50*t7*t77)-1.0d+0*t3*t62*t83-1.0d+0*t3*t46*t79-2
     4         .0d+0*t3*t38*t74-2.0d+0*t3*t40*t73)-6.0d+0*t17*t3*t46*t80
     5         *t87+t44*(2.0d+0*t17*t3*t46*t83-2*t40*t75+2.0d+0*t17*t3*t
     6         40*t74+2.0d+0*t3*t38*t40*t46)+t19*t20*(-1.333333333333333
     7         3d+0*t17*t18*t83-2.6666666666666666d+0*t18*t38*t40)+2.666
     8         6666666666666d+0*t17*t18*t19*t44*t80+1.3333333333333333d+
     9         0*t13*t18*t19*t79+1.0d+0*t13*t3*(-6.666666666666666d-7*t3
     :         5*t45*t50*t7*t77-3.610227d+1*gammaaa*t15*t16*t19*t72+1.75
     ;         52d+1*t15*t16*t19*t31+1.1881257057d+1*t15*t16*t19*t30/t2*
     <         *3.3d+1))*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = (t20*(-1.0d+0*t17*t3*(5.373306513
     1         388233d-7*t35*t4*t7*t76*t98+2.25d+0*t78*t8*t95+2.25d+0*t1
     2         0*t11*t8/t9**5-7.5d-1*t37*t66*t8*t81-2.25d+0*t14*t22*t8/t
     3         30)-3.0d+0*t3*t38*t83-3.0d+0*t3*t40*t79)+1.0d+0*t13*t3*(5
     4         .0d-7*t35*t4*t7*t76*t98+8.12301075d+0*t15*t16*t78-4.45547
     5         1396375d+0*gammaaa*t15*t16/t2**3.2d+1)-6.0d+0*t17*t3*t40*
     6         *3*t87+t44*(-2*t40*t84+4.0d+0*t17*t3*t40*t83+2.0d+0*t3*t3
     7         8*t80))*wght+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t21 = (gammabb*t14+1)**5.0d-1
            t22 = 1/t21
            t23 = 1/t2**9.0d+0
            t24 = -8.0d+0*gammabb*t19*t22*t23*t8
            t25 = 1/t2**5.0d+0
            t26 = -8.0d+0*t11*t19*t25*t8*t9
            t27 = 1.0d+0/t1
            t28 = -1.4328817369035288d-6*t19*t27*t4*t6*t7+t26+t24
            t29 = 2.6666666666666666d+0*gammabb*t19*t23*t8
            t30 = gammabb**2
            t31 = 1/t2**1.7d+1
            t32 = 4.3879999999999997d+0*t15*t16*t19*t30*t31
            t33 = -2.6666666666666666d+0*gammabb*t15*t16*t19*t23
            t34 = -1.3333333333333333d-6*t19*t27*t4*t6*t7+t33+t32+t29
            t35 = t5-1
            t36 = gammabb**t35
            t37 = 1/t2**1.6d+1
            t38 = -t14*t8+5.0d-7*t36*t4*t7-1.6455d+0*gammabb*t15*t16*t37
     1         +t14*t15*t16
            t39 = 1/t9
            t40 = 3.0d+0*t10*t11*t39*t8+3.0d+0*t14*t22*t8+5.373306513388
     1         233d-7*t36*t4*t7
            t41 = 1/rhob**1.6666666666666669d+0
            t42 = t2**2.0d+0
            t43 = 1/rhob**1.3333333333333333d+0
            t44 = 1/t12**3
            t45 = 1/rhob
            t46 = -1.4328817369035288d-6*t4*t45*t6*t7+t26+t24
            t47 = 1/rhob**2
            t48 = 1.3333333333333333d-6*t4*t47*t6*t7
            t49 = -1.7777777777777776d+0*gammabb*t23*t41*t8
            t50 = t4**2
            t51 = 1/t2**1.0d+1
            t52 = -8.0d+0*gammabb*t43*t51*t8
            t53 = -2.925333333333333d+0*t15*t16*t30*t31*t41
            t54 = 1.7777777777777776d+0*gammabb*t15*t16*t23*t41
            t55 = gammabb**3
            t56 = 1/t2**2.6d+1
            t57 = 1.9254543999999998d+1*t15*t16*t43*t55*t56
            t58 = 1/t2**1.8d+1
            t59 = -3.6566666666666664d+1*t15*t16*t30*t43*t58
            t60 = 8.0d+0*gammabb*t15*t16*t43*t51
            t61 = 1.7777777777777776d-6*t27*t41*t50*t6*t7+t60+t59+t57+t5
     1         4+t53+t52+t49+t48
            t62 = -1.3333333333333333d-6*t4*t45*t6*t7+t33+t32+t29
            t63 = 1.4328817369035288d-6*t4*t47*t6*t7
            t64 = 5.333333333333333d+0*gammabb*t22*t23*t41*t8
            t65 = 5.333333333333333d+0*t11*t25*t41*t8*t9
            t66 = 1/t21**3
            t67 = -1.0666666666666666d+1*t30*t43*t58*t66*t8
            t68 = 3.466666666666666d+1*gammabb*t22*t43*t51*t8
            t69 = 1/t2**6.0d+0
            t70 = 1.3333333333333333d+1*t11*t43*t69*t8*t9
            t71 = t70+1.9105089825380384d-6*t27*t41*t50*t6*t7+t68+t67+t6
     1         5+t64+t63
            t72 = 1/t2**2.5d+1
            t73 = 2.6666666666666666d+0*t19*t23*t8-7.220454d+0*t15*t16*t
     1         19*t30*t72-6.666666666666666d-7*t36*t45*t50*t7+1.3164d+1*
     2         gammabb*t15*t16*t19*t31-2.6666666666666666d+0*t15*t16*t19
     3         *t23
            t74 = 4.0d+0*gammabb*t19*t31*t66*t8-4.0d+0*t11*t19*t25*t39*t
     1         8-1.2d+1*t19*t22*t23*t8-7.164408684517644d-7*t36*t45*t50*
     2         t7
            t75 = -1.0d+0*t17*t3*t74-1.0d+0*t3*t40*t62-1.0d+0*t3*t38*t46
            t76 = t5-2
            t77 = gammabb**t76
            t78 = 1/t2**2.4d+1
            t79 = 2.70767025d+0*gammabb*t15*t16*t78+5.0d-7*t35*t4*t7*t77
     1         -3.291d+0*t15*t16*t37
            t80 = t40**2
            t81 = 1/gammabb
            t82 = 1/t9**3
            t83 = -1.5d+0*t10*t11*t8*t82+1.5d+0*t14*t22*t8*t81-1.5d+0*t3
     1         7*t66*t8+5.373306513388233d-7*t35*t4*t7*t77
            t84 = -1.0d+0*t17*t3*t83-2.0d+0*t3*t38*t40
            t85 = 1/rhob**2.6666666666666666d+0
            t86 = 1/rhob**2.3333333333333334d+0
            t87 = 1/t12**4
            t88 = t46**2
            t89 = 1/rhob**3
            t90 = t4**3
            t91 = 1/t2**1.1d+1
            t92 = 1/t2**2.7d+1
            t93 = 1/t2**1.9d+1
            t94 = 1.7777777777777776d-6*t47*t50*t6*t7+t60+t59+t57+t54+t5
     1         3+t52+t49+t48
            t95 = 1/t21**5
            t96 = t70+1.9105089825380384d-6*t47*t50*t6*t7+t68+t67+t65+t6
     1         4+t63
            t97 = -1.0d+0*t17*t3*t96-2.0d+0*t3*t46*t62
            t98 = gammabb**(t5-3)
            fnc(iq) = (1.0d+0*t13*t17*t3-9.305257363491d-1*t3)*wght+fnc(
     1         iq)
            Amat(iq,D1_RB) = (1.0d+0*t13*t3*t34-1.0d+0*t17*t20*t28*t3+1.
     1         3333333333333333d+0*t13*t17*t18*t19-1.2407009817987999d+0
     2         *t18*t19)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = (1.0d+0*t13*t3*t38-1.0d+0*t17*t20*t3*t40)*
     1         wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = (t20*(-1.0d+0*t17*t3*t71-1.0d+0*t28*t3*
     1         t62-1.0d+0*t3*t34*t46)+t13*t19*(1.3333333333333333d+0*t18
     2         *t62+1.3333333333333333d+0*t18*t34)+1.0d+0*t13*t3*t61+t19
     3         *t20*(-1.3333333333333333d+0*t17*t18*t46-1.33333333333333
     4         33d+0*t17*t18*t28)+2.0d+0*t17*t28*t3*t44*t46+1.3333333333
     5         333333d+0*t13*t17*t42*t43-1.2407009817987999d+0*t42*t43-8
     6         .888888888888888d-1*t13*t17*t18*t41+8.271339878658666d-1*
     7         t18*t41)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (t20*t75+1.0d+0*t13*t3*t73+2.0d+0*t17*
     1         t3*t40*t44*t46-1.3333333333333333d+0*t17*t18*t19*t20*t40+
     2         1.3333333333333333d+0*t13*t18*t19*t38)*wght+Cmat2(iq,D2_R
     3         B_GBB)
            Cmat2(iq,D2_GBB_GBB) = (t20*t84+2.0d+0*t17*t3*t44*t80+1.0d+0
     1         *t13*t3*t79)*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = (t44*(-2*t28*t97+2.0d+0*t3*t34*t88+4
     1         .0d+0*t17*t3*t46*t71)+t20*(-1.0d+0*t3*t34*t96-1.333333333
     2         3333333d+0*t17*t18*t19*t96-1.0d+0*t17*t3*(-4.266666666666
     3         6664d+1*t47*t55*t8*t92*t95+1.1022222222222221d+2*t30*t47*
     4         t66*t8*t93-1.333333333333333d+2*gammabb*t22*t47*t8*t91-2.
     5         5473453100507176d-6*t27*t6*t7*t85*t90-2.6666666666666666d
     6         +1*t11*t69*t8*t86*t9-8.88888888888889d+0*t11*t25*t8*t85*t
     7         9-2.6666666666666666d+1*t11*t47*t8*t9/t2**7.0d+0-3.821017
     8         9650760767d-6*t50*t6*t7*t89-2.8657634738070575d-6*t4*t6*t
     9         7*t89+2.1333333333333332d+1*t30*t58*t66*t8*t86-6.93333333
     :         3333332d+1*gammabb*t22*t51*t8*t86-8.88888888888889d+0*gam
     ;         mabb*t22*t23*t8*t85-1.9105089825380384d-6*t27*t50*t6*t7*t
     <         85)-1.0d+0*t28*t3*t94-2.0d+0*t3*t62*t71-2.666666666666666
     =         6d+0*t18*t19*t46*t62-2.0d+0*t3*t46*t61)+t13*t19*(1.333333
     >         3333333333d+0*t18*t94+2.6666666666666666d+0*t19*t42*t62+2
     ?         .6666666666666666d+0*t18*t61)+1.0d+0*t13*t3*(2.5450399999
     @         999995d+2*t15*t16*t30*t47*t93-3.2732724799999996d+2*t15*t
     1         16*t47*t55*t92+2.6666666666666666d+1*gammabb*t47*t8*t91-2
     2         .6666666666666666d+1*gammabb*t15*t16*t47*t91-2.3703703703
     3         7037d-6*t27*t6*t7*t85*t90-3.555555555555555d-6*t50*t6*t7*
     4         t89-2.6666666666666666d-6*t4*t6*t7*t89+1.6d+1*gammabb*t51
     5         *t8*t86+7.313333333333333d+1*t15*t16*t30*t58*t86-3.850908
     6         7999999997d+1*t15*t16*t55*t56*t86-1.6d+1*gammabb*t15*t16*
     7         t51*t86+2.962962962962963d+0*gammabb*t23*t8*t85-1.7777777
     8         777777776d-6*t27*t50*t6*t7*t85+4.875555555555556d+0*t15*t
     9         16*t30*t31*t85-2.962962962962963d+0*gammabb*t15*t16*t23*t
     :         85+8.448893907199999d+1*gammabb**4*t15*t16*t47/t2**3.5d+1
     ;         )+t19*t44*(2.6666666666666666d+0*t17*t18*t88+5.3333333333
     <         33333d+0*t17*t18*t28*t46)-6.0d+0*t17*t28*t3*t87*t88-2.666
     =         6666666666666d+0*t13*t17*t42*t86+2.4814019635975998d+0*t4
     >         2*t86+1.4814814814814814d+0*t13*t17*t18*t85-1.37855664644
     ?         3111d+0*t18*t85+t19*t20*(-2.6666666666666666d+0*t17*t18*t
     @         71-2.6666666666666666d+0*t18*t28*t62-2.6666666666666666d+
     1         0*t17*t19*t42*t46-2.6666666666666666d+0*t18*t34*t46)-1.77
     2         77777777777776d+0*t13*t18*t41*t62+1.7777777777777776d+0*t
     3         17*t18*t20*t41*t46+1.3333333333333333d+0*t13*t34*t42*t43-
     4         1.3333333333333333d+0*t17*t20*t28*t42*t43-8.8888888888888
     5         88d-1*t13*t18*t34*t41+8.888888888888888d-1*t17*t18*t20*t2
     6         8*t41+8.888888888888888d-1*t13*t17*t41-8.271339878658666d
     7         -1*t41)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = (t44*(-2*t40*t97+2.0d+0*t3*t38*t88+
     1         4.0d+0*t17*t3*t46*t74)+t20*(-1.0d+0*t3*t38*t96-1.0d+0*t17
     2         *t3*(1.6d+1*t30*t43*t56*t8*t95+9.552544912690191d-7*t36*t
     3         47*t7*t90+6.666666666666666d+0*t11*t39*t43*t69*t8-3.86666
     4         66666666666d+1*gammabb*t43*t58*t66*t8-2.6666666666666666d
     5         +0*gammabb*t31*t41*t66*t8+4.133333333333333d+1*t22*t43*t5
     6         1*t8+2.6666666666666666d+0*t11*t25*t39*t41*t8+8.0d+0*t22*
     7         t23*t41*t8+7.164408684517644d-7*t36*t47*t50*t7)-1.0d+0*t3
     8         *t40*t94-2.0d+0*t3*t62*t74-2.0d+0*t3*t46*t73)+1.0d+0*t13*
     9         t3*(8.888888888888887d-7*t36*t47*t7*t90-8.0d+0*t43*t51*t8
     :         -1.7777777777777776d+0*t23*t41*t8+4.813636d+0*t15*t16*t30
     ;         *t41*t72+6.666666666666666d-7*t36*t47*t50*t7-8.6297333333
     <         33333d+1*gammabb*t15*t16*t43*t58+1.17934082d+2*t15*t16*t3
     =         0*t43*t56-3.1683352152d+1*t15*t16*t43*t55/t2**3.4d+1+8.0d
     >         +0*t15*t16*t43*t51-8.775999999999999d+0*gammabb*t15*t16*t
     ?         31*t41+1.7777777777777776d+0*t15*t16*t23*t41)-6.0d+0*t17*
     @         t3*t40*t87*t88+t19*t20*(-2.6666666666666666d+0*t17*t18*t7
     1         4-2.6666666666666666d+0*t18*t40*t62-2.6666666666666666d+0
     2         *t18*t38*t46)+2.6666666666666666d+0*t13*t18*t19*t73+5.333
     3         333333333333d+0*t17*t18*t19*t40*t44*t46-1.333333333333333
     4         3d+0*t17*t20*t40*t42*t43+1.3333333333333333d+0*t13*t38*t4
     5         2*t43+8.888888888888888d-1*t17*t18*t20*t40*t41-8.88888888
     6         8888888d-1*t13*t18*t38*t41)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = (t20*(-1.0d+0*t17*t3*(-6.0d+0*gamm
     1         abb*t19*t72*t8*t95+2.0d+0*t11*t19*t25*t8*t82-2.0d+0*t19*t
     2         22*t23*t8*t81+1.0d+1*t19*t31*t66*t8-7.164408684517644d-7*
     3         t35*t45*t50*t7*t77)-1.0d+0*t3*t62*t83-1.0d+0*t3*t46*t79-2
     4         .0d+0*t3*t38*t74-2.0d+0*t3*t40*t73)-6.0d+0*t17*t3*t46*t80
     5         *t87+t44*(2.0d+0*t17*t3*t46*t83-2*t40*t75+2.0d+0*t17*t3*t
     6         40*t74+2.0d+0*t3*t38*t40*t46)+t19*t20*(-1.333333333333333
     7         3d+0*t17*t18*t83-2.6666666666666666d+0*t18*t38*t40)+2.666
     8         6666666666666d+0*t17*t18*t19*t44*t80+1.3333333333333333d+
     9         0*t13*t18*t19*t79+1.0d+0*t13*t3*(-6.666666666666666d-7*t3
     :         5*t45*t50*t7*t77-3.610227d+1*gammabb*t15*t16*t19*t72+1.75
     ;         52d+1*t15*t16*t19*t31+1.1881257057d+1*t15*t16*t19*t30/t2*
     <         *3.3d+1))*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (t20*(-1.0d+0*t17*t3*(5.373306513
     1         388233d-7*t35*t4*t7*t76*t98+2.25d+0*t78*t8*t95+2.25d+0*t1
     2         0*t11*t8/t9**5-7.5d-1*t37*t66*t8*t81-2.25d+0*t14*t22*t8/t
     3         30)-3.0d+0*t3*t38*t83-3.0d+0*t3*t40*t79)+1.0d+0*t13*t3*(5
     4         .0d-7*t35*t4*t7*t76*t98+8.12301075d+0*t15*t16*t78-4.45547
     5         1396375d+0*gammabb*t15*t16/t2**3.2d+1)-6.0d+0*t17*t3*t40*
     6         *3*t87+t44*(-2*t40*t84+4.0d+0*t17*t3*t40*t83+2.0d+0*t3*t3
     7         8*t80))*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
