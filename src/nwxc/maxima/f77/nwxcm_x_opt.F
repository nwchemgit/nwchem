C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_opt.F 
C> The nwxcm_x_opt functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006
C>    \,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}
C>    -3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}
C>    \,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{
C>    \it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}
C>    \,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}
C>    \,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="https://doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_opt(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = (gammaaa*t1*t3)**2.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = 1/t2**2
            t7 = 1/t5
            fnc(iq) = fnc(iq)-7.199999999999999d-5*t4*t5*wght
            Amat(iq,D1_RA) = (-7.199999999999999d-5*gammaaa*t3*(1.6d-2*g
     1         ammaaa**2*t6/rhoa**6.333333333333333d+0-2.666666666666666
     2         6d+0*gammaaa*t3/rhoa**3.6666666666666664d+0)*t7-4.7999999
     3         999999987d-5*rhoa**3.333333333333333d-1*t4)*wght+Amat(iq,
     4         D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*gamma
     1         aa*t3*(t1*t3-6.0d-3*gammaaa*t6/rhoa**5.333333333333333d+0
     2         )*t7*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 6.0d-3*gammabb*t7+1.0d+0
            t9 = 1/t8
            t10 = gammabb*t7*t9
            t11 = t10**2.0d+0
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t2**2
            t14 = 1/t8**2
            fnc(iq) = (-3.5999999999999993d-5*t5*t6-3.5999999999999993d-
     1         5*t11*t12)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-7.199999999999999d-5*gammaaa*t3*(1.6d-2*g
     1         ammaaa**2*t13/rhoa**6.333333333333333d+0-2.66666666666666
     2         66d+0*gammaaa*t3/rhoa**3.6666666666666664d+0)/t6-4.799999
     3         9999999987d-5*rhoa**3.333333333333333d-1*t5)*wght+Amat(iq
     4         ,D1_RA)
            Amat(iq,D1_RB) = (-7.199999999999999d-5*gammabb*t9*(1.6d-2*g
     1         ammabb**2*t14/rhob**6.333333333333333d+0-2.66666666666666
     2         66d+0*gammabb*t9/rhob**3.6666666666666664d+0)/t12-4.79999
     3         99999999987d-5*rhob**3.333333333333333d-1*t11)*wght+Amat(
     4         iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*(t1*t
     1         3-6.0d-3*gammaaa*t13/rhoa**5.333333333333333d+0)*t4*t6*wg
     2         ht
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t10*t
     1         12*(t7*t9-6.0d-3*gammabb*t14/rhob**5.333333333333333d+0)*
     2         wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t2**2
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RA) = -7.199999999999999d-5*gammaaa*t3*(1.6d-2*ga
     1         mmaaa**2*t7/rhoa**6.333333333333333d+0-2.6666666666666666
     2         d+0*gammaaa*t3/rhoa**3.6666666666666664d+0)*wght/t6-4.799
     3         9999999999987d-5*rhoa**3.333333333333333d-1*t5*wght+Amat(
     4         iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t4*t6
     1         *(t1*t3-6.0d-3*gammaaa*t7/rhoa**5.333333333333333d+0)*wgh
     2         t
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 6.0d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = gammabb*t1*t3
            t5 = t4**2.0d+0
            t6 = rhob**1.3333333333333333d+0
            t7 = 1/t2**2
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RB) = -7.199999999999999d-5*gammabb*t3*(1.6d-2*ga
     1         mmabb**2*t7/rhob**6.333333333333333d+0-2.6666666666666666
     2         d+0*gammabb*t3/rhob**3.6666666666666664d+0)*wght/t6-4.799
     3         9999999999987d-5*rhob**3.333333333333333d-1*t5*wght+Amat(
     4         iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t4*t6
     1         *(t1*t3-6.0d-3*gammabb*t7/rhob**5.333333333333333d+0)*wgh
     2         t
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006
C>    \,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}
C>    -3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}
C>    \,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{
C>    \it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}
C>    \,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}
C>    \,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="https://doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_opt_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = (gammaaa*t1*t3)**2.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = gammaaa**2
            t7 = 1/t2**2
            t8 = 1/rhoa**6.333333333333333d+0
            t9 = 1/rhoa**3.6666666666666664d+0
            t10 = 1.6d-2*t6*t7*t8-2.6666666666666666d+0*gammaaa*t3*t9
            t11 = 1/t5
            t12 = 1/rhoa**5.333333333333333d+0
            t13 = t1*t3-6.0d-3*gammaaa*t12*t7
            t14 = 1/t2**3
            fnc(iq) = fnc(iq)-7.199999999999999d-5*t4*t5*wght
            Amat(iq,D1_RA) = (-4.7999999999999987d-5*rhoa**3.33333333333
     1         3333d-1*t4-7.199999999999999d-5*gammaaa*t10*t11*t3)*wght+
     2         Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*gamma
     1         aa*t11*t13*t3*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-7.199999999999999d-5*gammaaa*t11*t3*(
     1         -1.44d-1*t6*t7/rhoa**7.333333333333333d+0+9.7777777777777
     2         79d+0*gammaaa*t3/rhoa**4.666666666666667d+0+5.12d-4*gamma
     3         aa**3*t14/rhoa**10)-1.152d-6*t10*t6*t7/rhoa**5-1.59999999
     4         99999995d-5*t4/rhoa**6.666666666666666d-1)*wght+Amat2(iq,
     5         D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t11*t3*(-7.199999999999999d-5*gammaaa
     1         *(-2.6666666666666666d+0*t3*t9+4.8d-2*gammaaa*t7*t8-1.92d
     2         -4*t14*t6/rhoa**9)-7.199999999999999d-5*t10)+4.3199999999
     3         999994d-7*gammaaa*t10*t7/rhoa**4-9.599999999999997d-5*gam
     4         maaa*t13*t3/rhoa**2.3333333333333334d+0)*wght+Cmat2(iq,D2
     5         _RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*gammaaa*t11*t3*
     1         (7.2d-5*gammaaa*t14/rhoa**8-1.2d-2*t12*t7)*wght-7.1999999
     2         99999999d-5*t13**2*t5*wght+Cmat2(iq,D2_GAA_GAA)
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 6.0d-3*gammabb*t7+1.0d+0
            t9 = 1/t8
            t10 = gammabb*t7*t9
            t11 = t10**2.0d+0
            t12 = rhob**1.3333333333333333d+0
            t13 = gammaaa**2
            t14 = 1/t2**2
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = 1.6d-2*t13*t14*t15-2.6666666666666666d+0*gammaaa*t16*t
     1         3
            t18 = 1/t6
            t19 = gammabb**2
            t20 = 1/t8**2
            t21 = 1/rhob**6.333333333333333d+0
            t22 = 1/rhob**3.6666666666666664d+0
            t23 = 1.6d-2*t19*t20*t21-2.6666666666666666d+0*gammabb*t22*t
     1         9
            t24 = 1/t12
            t25 = 1/rhoa**5.333333333333333d+0
            t26 = t1*t3-6.0d-3*gammaaa*t14*t25
            t27 = 1/rhob**5.333333333333333d+0
            t28 = t7*t9-6.0d-3*gammabb*t20*t27
            t29 = 1/t2**3
            t30 = 1/t8**3
            fnc(iq) = (-3.5999999999999993d-5*t5*t6-3.5999999999999993d-
     1         5*t11*t12)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-4.7999999999999987d-5*rhoa**3.33333333333
     1         3333d-1*t5-7.199999999999999d-5*gammaaa*t17*t18*t3)*wght+
     2         Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-7.199999999999999d-5*gammabb*t23*t24*t9-4
     1         .7999999999999987d-5*rhob**3.333333333333333d-1*t11)*wght
     2         +Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t26*t
     1         4*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t10*t
     1         12*t28*wght
            Amat2(iq,D2_RA_RA) = (-1.5999999999999995d-5*t5/rhoa**6.6666
     1         66666666666d-1-7.199999999999999d-5*gammaaa*t18*t3*(9.777
     2         777777777779d+0*gammaaa*t3/rhoa**4.666666666666667d+0+5.1
     3         2d-4*gammaaa**3*t29/rhoa**10-1.44d-1*t13*t14/rhoa**7.3333
     4         33333333333d+0)-1.152d-6*t13*t14*t17/rhoa**5)*wght+Amat2(
     5         iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-7.199999999999999d-5*gammabb*t24*t9*(
     1         9.777777777777779d+0*gammabb*t9/rhob**4.666666666666667d+
     2         0+5.12d-4*gammabb**3*t30/rhob**10-1.44d-1*t19*t20/rhob**7
     3         .333333333333333d+0)-1.152d-6*t19*t20*t23/rhob**5-1.59999
     4         99999999995d-5*t11/rhob**6.666666666666666d-1)*wght+Amat2
     5         (iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t18*t3*(-7.199999999999999d-5*gammaaa
     1         *(-2.6666666666666666d+0*t16*t3-1.92d-4*t13*t29/rhoa**9+4
     2         .8d-2*gammaaa*t14*t15)-7.199999999999999d-5*t17)-9.599999
     3         999999997d-5*gammaaa*t26*t3/rhoa**2.3333333333333334d+0+4
     4         .3199999999999994d-7*gammaaa*t14*t17/rhoa**4)*wght+Cmat2(
     5         iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t24*t9*(-7.199999999999999d-5*gammabb
     1         *(-2.6666666666666666d+0*t22*t9-1.92d-4*t19*t30/rhob**9+4
     2         .8d-2*gammabb*t20*t21)-7.199999999999999d-5*t23)-9.599999
     3         999999997d-5*gammabb*t28*t9/rhob**2.3333333333333334d+0+4
     4         .3199999999999994d-7*gammabb*t20*t23/rhob**4)*wght+Cmat2(
     5         iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*t26**2*t6*wght-
     1         7.199999999999999d-5*gammaaa*t18*(7.2d-5*gammaaa*t29/rhoa
     2         **8-1.2d-2*t14*t25)*t3*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.199999999999999d-5*gammabb*t24*(7.
     1         2d-5*gammabb*t30/rhob**8-1.2d-2*t20*t27)*t9*wght-7.199999
     2         999999999d-5*t12*t28**2*wght+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = gammaaa**2
            t8 = 1/t2**2
            t9 = 1/rhoa**6.333333333333333d+0
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = 1.6d-2*t7*t8*t9-2.6666666666666666d+0*gammaaa*t10*t3
            t12 = 1/t6
            t13 = 1/rhoa**5.333333333333333d+0
            t14 = t1*t3-6.0d-3*gammaaa*t13*t8
            t15 = 1/t2**3
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RA) = -4.7999999999999987d-5*rhoa**3.333333333333
     1         333d-1*t5*wght-7.199999999999999d-5*gammaaa*t11*t12*t3*wg
     2         ht+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t14*t
     1         4*t6*wght
            Amat2(iq,D2_RA_RA) = -7.199999999999999d-5*gammaaa*t12*t3*(-
     1         1.44d-1*t7*t8/rhoa**7.333333333333333d+0+9.77777777777777
     2         9d+0*gammaaa*t3/rhoa**4.666666666666667d+0+5.12d-4*gammaa
     3         a**3*t15/rhoa**10)*wght-1.152d-6*t11*t7*t8*wght/rhoa**5-1
     4         .5999999999999995d-5*t5*wght/rhoa**6.666666666666666d-1+A
     5         mat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t12*t3*(-7.199999999999999d-5*gammaaa*
     1         (4.8d-2*gammaaa*t8*t9-1.92d-4*t15*t7/rhoa**9-2.6666666666
     2         666666d+0*t10*t3)*wght-7.199999999999999d-5*t11*wght)+4.3
     3         199999999999994d-7*gammaaa*t11*t8*wght/rhoa**4-9.59999999
     4         9999997d-5*gammaaa*t14*t3*wght/rhoa**2.3333333333333334d+
     5         0+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*gammaaa*t12*t3*
     1         (7.2d-5*gammaaa*t15/rhoa**8-1.2d-2*t13*t8)*wght-7.1999999
     2         99999999d-5*t14**2*t6*wght+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 6.0d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = gammabb*t1*t3
            t5 = t4**2.0d+0
            t6 = rhob**1.3333333333333333d+0
            t7 = gammabb**2
            t8 = 1/t2**2
            t9 = 1/rhob**6.333333333333333d+0
            t10 = 1/rhob**3.6666666666666664d+0
            t11 = 1.6d-2*t7*t8*t9-2.6666666666666666d+0*gammabb*t10*t3
            t12 = 1/t6
            t13 = 1/rhob**5.333333333333333d+0
            t14 = t1*t3-6.0d-3*gammabb*t13*t8
            t15 = 1/t2**3
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RB) = -4.7999999999999987d-5*rhob**3.333333333333
     1         333d-1*t5*wght-7.199999999999999d-5*gammabb*t11*t12*t3*wg
     2         ht+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t14*t
     1         4*t6*wght
            Amat2(iq,D2_RB_RB) = -7.199999999999999d-5*gammabb*t12*t3*(-
     1         1.44d-1*t7*t8/rhob**7.333333333333333d+0+9.77777777777777
     2         9d+0*gammabb*t3/rhob**4.666666666666667d+0+5.12d-4*gammab
     3         b**3*t15/rhob**10)*wght-1.152d-6*t11*t7*t8*wght/rhob**5-1
     4         .5999999999999995d-5*t5*wght/rhob**6.666666666666666d-1+A
     5         mat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t12*t3*(-7.199999999999999d-5*gammabb*
     1         (4.8d-2*gammabb*t8*t9-1.92d-4*t15*t7/rhob**9-2.6666666666
     2         666666d+0*t10*t3)*wght-7.199999999999999d-5*t11*wght)+4.3
     3         199999999999994d-7*gammabb*t11*t8*wght/rhob**4-9.59999999
     4         9999997d-5*gammabb*t14*t3*wght/rhob**2.3333333333333334d+
     5         0+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.199999999999999d-5*gammabb*t12*t3*
     1         (7.2d-5*gammabb*t15/rhob**8-1.2d-2*t13*t8)*wght-7.1999999
     2         99999999d-5*t14**2*t6*wght+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006
C>    \,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}
C>    -3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}
C>    \,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{
C>    \it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}
C>    \,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}
C>    \,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="https://doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_opt_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = (gammaaa*t1*t3)**2.0d+0
            t5 = rhoa**1.3333333333333333d+0
            t6 = gammaaa**2
            t7 = 1/t2**2
            t8 = 1/rhoa**6.333333333333333d+0
            t9 = 1/rhoa**3.6666666666666664d+0
            t10 = 1.6d-2*t6*t7*t8-2.6666666666666666d+0*gammaaa*t3*t9
            t11 = 1/t5
            t12 = 1/rhoa**5.333333333333333d+0
            t13 = t1*t3-6.0d-3*gammaaa*t12*t7
            t14 = 1/rhoa**5
            t15 = gammaaa**3
            t16 = 1/t2**3
            t17 = 1/rhoa**10
            t18 = 1/rhoa**7.333333333333333d+0
            t19 = 1/rhoa**4.666666666666667d+0
            t20 = -1.44d-1*t18*t6*t7+9.777777777777779d+0*gammaaa*t19*t3
     1         +5.12d-4*t15*t16*t17
            t21 = 1/rhoa**4
            t22 = 1/rhoa**2.3333333333333334d+0
            t23 = 1/rhoa**9
            t24 = -2.6666666666666666d+0*t3*t9+4.8d-2*gammaaa*t7*t8-1.92
     1         d-4*t16*t23*t6
            t25 = -7.199999999999999d-5*gammaaa*t24-7.199999999999999d-5
     1         *t10
            t26 = 1/rhoa**8
            t27 = 7.2d-5*gammaaa*t16*t26-1.2d-2*t12*t7
            t28 = 1/rhoa**3.3333333333333337d+0
            t29 = 1/t2**4
            fnc(iq) = fnc(iq)-7.199999999999999d-5*t4*t5*wght
            Amat(iq,D1_RA) = (-4.7999999999999987d-5*rhoa**3.33333333333
     1         3333d-1*t4-7.199999999999999d-5*gammaaa*t10*t11*t3)*wght+
     2         Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*gamma
     1         aa*t11*t13*t3*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-1.152d-6*t10*t14*t6*t7-1.599999999999
     1         9995d-5*t4/rhoa**6.666666666666666d-1-7.199999999999999d-
     2         5*gammaaa*t11*t20*t3)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (4.3199999999999994d-7*gammaaa*t10*t21
     1         *t7+t11*t25*t3-9.599999999999997d-5*gammaaa*t13*t22*t3)*w
     2         ght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*t13**2*t5*wght-
     1         7.199999999999999d-5*gammaaa*t11*t27*t3*wght+Cmat2(iq,D2_
     2         GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-7.199999999999999d-5*gammaaa*t11*t
     1         3*(1.2124444444444442d+0*t6*t7/rhoa**8.333333333333334d+0
     2         -4.5629629629629626d+1*gammaaa*t3/rhoa**5.666666666666667
     3         d+0+2.4576d-5*gammaaa**4*t29/rhoa**1.3666666666666666d+1-
     4         9.728000000000001d-3*t15*t16/rhoa**11)-2.304d-6*t14*t20*t
     5         6*t7+5.76d-6*t10*t6*t7/rhoa**6+1.0666666666666665d-5*t4/r
     6         hoa**1.6666666666666669d+0-3.199999999999999d-5*gammaaa*t
     7         10*t28*t3+9.599999999999997d-5*gammaaa*t20*t22*t3-3.68639
     8         99999999997d-8*t10*t15*t16/rhoa**8.666666666666666d+0)*wg
     9         ht+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t11*t3*(-7.199999999999999d-5*gamm
     1         aaa*(-3.466666666666667d-1*gammaaa*t18*t7+3.264d-3*t16*t1
     2         7*t6+9.777777777777779d+0*t19*t3-9.216000000000001d-6*t15
     3         *t29/rhoa**1.2666666666666665d+1)-7.199999999999999d-5*t2
     4         0)+t14*(-1.152d-6*t24*t6-2.304d-6*gammaaa*t10)*t7+4.31999
     5         99999999994d-7*gammaaa*t20*t21*t7+1.3824d-8*t10*t16*t6/rh
     6         oa**7.666666666666667d+0-3.199999999999999d-5*gammaaa*t13
     7         *t28*t3)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t11*t3*(-7.199999999999999d-5*gam
     1         maaa*(6.4d-2*t7*t8+3.456d-6*t29*t6/rhoa**1.16666666666666
     2         65d+1-9.600000000000001d-4*gammaaa*t16*t23)-1.44d-4*t24)+
     3         t21*(-6.0d-3*t25+4.3199999999999994d-7*gammaaa*t24+4.3199
     4         999999999994d-7*t10)*t7+5.759999999999998d-7*gammaaa*t13*
     5         t14*t7+t22*(-9.599999999999997d-5*gammaaa*t27-9.599999999
     6         999997d-5*t13)*t3-5.184d-9*gammaaa*t10*t16/rhoa**6.666666
     7         666666667d+0)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t11*t3*(-7.199999999999999d-5*gam
     1         maaa*(2.16d-4*t16*t26-1.2960000000000002d-6*gammaaa*t29/r
     2         hoa**1.0666666666666666d+1)*wght-7.199999999999999d-5*t27
     3         *wght)+4.3199999999999994d-7*gammaaa*t21*t27*t7*wght-1.44
     4         d-4*t13*t27*t5*wght+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 6.0d-3*gammabb*t7+1.0d+0
            t9 = 1/t8
            t10 = gammabb*t7*t9
            t11 = t10**2.0d+0
            t12 = rhob**1.3333333333333333d+0
            t13 = gammaaa**2
            t14 = 1/t2**2
            t15 = 1/rhoa**6.333333333333333d+0
            t16 = 1/rhoa**3.6666666666666664d+0
            t17 = 1.6d-2*t13*t14*t15-2.6666666666666666d+0*gammaaa*t16*t
     1         3
            t18 = 1/t6
            t19 = gammabb**2
            t20 = 1/t8**2
            t21 = 1/rhob**6.333333333333333d+0
            t22 = 1/rhob**3.6666666666666664d+0
            t23 = 1.6d-2*t19*t20*t21-2.6666666666666666d+0*gammabb*t22*t
     1         9
            t24 = 1/t12
            t25 = 1/rhoa**5.333333333333333d+0
            t26 = t1*t3-6.0d-3*gammaaa*t14*t25
            t27 = 1/rhob**5.333333333333333d+0
            t28 = t7*t9-6.0d-3*gammabb*t20*t27
            t29 = 1/rhoa**5
            t30 = gammaaa**3
            t31 = 1/t2**3
            t32 = 1/rhoa**10
            t33 = 1/rhoa**7.333333333333333d+0
            t34 = 1/rhoa**4.666666666666667d+0
            t35 = 9.777777777777779d+0*gammaaa*t3*t34-1.44d-1*t13*t14*t3
     1         3+5.12d-4*t30*t31*t32
            t36 = 1/rhob**5
            t37 = gammabb**3
            t38 = 1/t8**3
            t39 = 1/rhob**10
            t40 = 1/rhob**7.333333333333333d+0
            t41 = 1/rhob**4.666666666666667d+0
            t42 = 9.777777777777779d+0*gammabb*t41*t9-1.44d-1*t19*t20*t4
     1         0+5.12d-4*t37*t38*t39
            t43 = 1/rhoa**4
            t44 = 1/rhoa**2.3333333333333334d+0
            t45 = 1/rhoa**9
            t46 = -1.92d-4*t13*t31*t45-2.6666666666666666d+0*t16*t3+4.8d
     1         -2*gammaaa*t14*t15
            t47 = -7.199999999999999d-5*gammaaa*t46-7.199999999999999d-5
     1         *t17
            t48 = 1/rhob**4
            t49 = 1/rhob**2.3333333333333334d+0
            t50 = 1/rhob**9
            t51 = -2.6666666666666666d+0*t22*t9-1.92d-4*t19*t38*t50+4.8d
     1         -2*gammabb*t20*t21
            t52 = -7.199999999999999d-5*gammabb*t51-7.199999999999999d-5
     1         *t23
            t53 = 1/rhoa**8
            t54 = 7.2d-5*gammaaa*t31*t53-1.2d-2*t14*t25
            t55 = 1/rhob**8
            t56 = 7.2d-5*gammabb*t38*t55-1.2d-2*t20*t27
            t57 = 1/rhoa**3.3333333333333337d+0
            t58 = 1/t2**4
            t59 = 1/rhob**3.3333333333333337d+0
            t60 = 1/t8**4
            fnc(iq) = (-3.5999999999999993d-5*t5*t6-3.5999999999999993d-
     1         5*t11*t12)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-4.7999999999999987d-5*rhoa**3.33333333333
     1         3333d-1*t5-7.199999999999999d-5*gammaaa*t17*t18*t3)*wght+
     2         Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-7.199999999999999d-5*gammabb*t23*t24*t9-4
     1         .7999999999999987d-5*rhob**3.333333333333333d-1*t11)*wght
     2         +Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t26*t
     1         4*t6*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t10*t
     1         12*t28*wght
            Amat2(iq,D2_RA_RA) = (-1.5999999999999995d-5*t5/rhoa**6.6666
     1         66666666666d-1-7.199999999999999d-5*gammaaa*t18*t3*t35-1.
     2         152d-6*t13*t14*t17*t29)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-7.199999999999999d-5*gammabb*t24*t42*
     1         t9-1.152d-6*t19*t20*t23*t36-1.5999999999999995d-5*t11/rho
     2         b**6.666666666666666d-1)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t18*t3*t47-9.599999999999997d-5*gamma
     1         aa*t26*t3*t44+4.3199999999999994d-7*gammaaa*t14*t17*t43)*
     2         wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t24*t52*t9-9.599999999999997d-5*gamma
     1         bb*t28*t49*t9+4.3199999999999994d-7*gammabb*t20*t23*t48)*
     2         wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*t26**2*t6*wght-
     1         7.199999999999999d-5*gammaaa*t18*t3*t54*wght+Cmat2(iq,D2_
     2         GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.199999999999999d-5*gammabb*t24*t56
     1         *t9*wght-7.199999999999999d-5*t12*t28**2*wght+Cmat2(iq,D2
     2         _GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (-7.199999999999999d-5*gammaaa*t18*t
     1         3*(2.4576d-5*gammaaa**4*t58/rhoa**1.3666666666666666d+1-9
     2         .728000000000001d-3*t30*t31/rhoa**11-4.5629629629629626d+
     3         1*gammaaa*t3/rhoa**5.666666666666667d+0+1.212444444444444
     4         2d+0*t13*t14/rhoa**8.333333333333334d+0)-3.19999999999999
     5         9d-5*gammaaa*t17*t3*t57+1.0666666666666665d-5*t5/rhoa**1.
     6         6666666666666669d+0+9.599999999999997d-5*gammaaa*t3*t35*t
     7         44-2.304d-6*t13*t14*t29*t35-3.6863999999999997d-8*t17*t30
     8         *t31/rhoa**8.666666666666666d+0+5.76d-6*t13*t14*t17/rhoa*
     9         *6)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (-7.199999999999999d-5*gammabb*t24*t
     1         9*(-4.5629629629629626d+1*gammabb*t9/rhob**5.666666666666
     2         667d+0+2.4576d-5*gammabb**4*t60/rhob**1.3666666666666666d
     3         +1-9.728000000000001d-3*t37*t38/rhob**11+1.21244444444444
     4         42d+0*t19*t20/rhob**8.333333333333334d+0)-3.1999999999999
     5         99d-5*gammabb*t23*t59*t9+9.599999999999997d-5*gammabb*t42
     6         *t49*t9-2.304d-6*t19*t20*t36*t42-3.6863999999999997d-8*t2
     7         3*t37*t38/rhob**8.666666666666666d+0+5.76d-6*t19*t20*t23/
     8         rhob**6+1.0666666666666665d-5*t11/rhob**1.666666666666666
     9         9d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t18*t3*(-7.199999999999999d-5*gamm
     1         aaa*(-9.216000000000001d-6*t30*t58/rhoa**1.26666666666666
     2         65d+1+9.777777777777779d+0*t3*t34-3.466666666666667d-1*ga
     3         mmaaa*t14*t33+3.264d-3*t13*t31*t32)-7.199999999999999d-5*
     4         t35)-3.199999999999999d-5*gammaaa*t26*t3*t57+t14*t29*(-1.
     5         152d-6*t13*t46-2.304d-6*gammaaa*t17)+4.3199999999999994d-
     6         7*gammaaa*t14*t35*t43+1.3824d-8*t13*t17*t31/rhoa**7.66666
     7         6666666667d+0)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t24*t9*(-7.199999999999999d-5*gamm
     1         abb*(9.777777777777779d+0*t41*t9-9.216000000000001d-6*t37
     2         *t60/rhob**1.2666666666666665d+1-3.466666666666667d-1*gam
     3         mabb*t20*t40+3.264d-3*t19*t38*t39)-7.199999999999999d-5*t
     4         42)-3.199999999999999d-5*gammabb*t28*t59*t9+t20*t36*(-1.1
     5         52d-6*t19*t51-2.304d-6*gammabb*t23)+4.3199999999999994d-7
     6         *gammabb*t20*t42*t48+1.3824d-8*t19*t23*t38/rhob**7.666666
     7         666666667d+0)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t18*t3*(-7.199999999999999d-5*gam
     1         maaa*(3.456d-6*t13*t58/rhoa**1.1666666666666665d+1-9.6000
     2         00000000001d-4*gammaaa*t31*t45+6.4d-2*t14*t15)-1.44d-4*t4
     3         6)+t3*t44*(-9.599999999999997d-5*gammaaa*t54-9.5999999999
     4         99997d-5*t26)+t14*t43*(-6.0d-3*t47+4.3199999999999994d-7*
     5         gammaaa*t46+4.3199999999999994d-7*t17)-5.184d-9*gammaaa*t
     6         17*t31/rhoa**6.666666666666667d+0+5.759999999999998d-7*ga
     7         mmaaa*t14*t26*t29)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t24*(-7.199999999999999d-5*gammab
     1         b*(3.456d-6*t19*t60/rhob**1.1666666666666665d+1-9.6000000
     2         00000001d-4*gammabb*t38*t50+6.4d-2*t20*t21)-1.44d-4*t51)*
     3         t9+t49*(-9.599999999999997d-5*gammabb*t56-9.5999999999999
     4         97d-5*t28)*t9+t20*t48*(-6.0d-3*t52+4.3199999999999994d-7*
     5         gammabb*t51+4.3199999999999994d-7*t23)-5.184d-9*gammabb*t
     6         23*t38/rhob**6.666666666666667d+0+5.759999999999998d-7*ga
     7         mmabb*t20*t28*t36)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t18*t3*(-7.199999999999999d-5*gam
     1         maaa*(2.16d-4*t31*t53-1.2960000000000002d-6*gammaaa*t58/r
     2         hoa**1.0666666666666666d+1)*wght-7.199999999999999d-5*t54
     3         *wght)-1.44d-4*t26*t54*t6*wght+4.3199999999999994d-7*gamm
     4         aaa*t14*t43*t54*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t24*t9*(-7.199999999999999d-5*gam
     1         mabb*(2.16d-4*t38*t55-1.2960000000000002d-6*gammabb*t60/r
     2         hob**1.0666666666666666d+1)*wght-7.199999999999999d-5*t56
     3         *wght)+4.3199999999999994d-7*gammabb*t20*t48*t56*wght-1.4
     4         4d-4*t12*t28*t56*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = gammaaa**2
            t8 = 1/t2**2
            t9 = 1/rhoa**6.333333333333333d+0
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = 1.6d-2*t7*t8*t9-2.6666666666666666d+0*gammaaa*t10*t3
            t12 = 1/t6
            t13 = 1/rhoa**5.333333333333333d+0
            t14 = t1*t3-6.0d-3*gammaaa*t13*t8
            t15 = 1/rhoa**5
            t16 = gammaaa**3
            t17 = 1/t2**3
            t18 = 1/rhoa**10
            t19 = 1/rhoa**7.333333333333333d+0
            t20 = 1/rhoa**4.666666666666667d+0
            t21 = -1.44d-1*t19*t7*t8+9.777777777777779d+0*gammaaa*t20*t3
     1         +5.12d-4*t16*t17*t18
            t22 = 1/rhoa**4
            t23 = 1/rhoa**2.3333333333333334d+0
            t24 = 1/rhoa**9
            t25 = 4.8d-2*gammaaa*t8*t9-1.92d-4*t17*t24*t7-2.666666666666
     1         6666d+0*t10*t3
            t26 = -7.199999999999999d-5*gammaaa*t25*wght-7.1999999999999
     1         99d-5*t11*wght
            t27 = 1/rhoa**8
            t28 = 7.2d-5*gammaaa*t17*t27-1.2d-2*t13*t8
            t29 = 1/rhoa**3.3333333333333337d+0
            t30 = 1/t2**4
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RA) = -4.7999999999999987d-5*rhoa**3.333333333333
     1         333d-1*t5*wght-7.199999999999999d-5*gammaaa*t11*t12*t3*wg
     2         ht+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t14*t
     1         4*t6*wght
            Amat2(iq,D2_RA_RA) = -1.152d-6*t11*t15*t7*t8*wght-1.59999999
     1         99999995d-5*t5*wght/rhoa**6.666666666666666d-1-7.19999999
     2         9999999d-5*gammaaa*t12*t21*t3*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 4.3199999999999994d-7*gammaaa*t11*t22*
     1         t8*wght-9.599999999999997d-5*gammaaa*t14*t23*t3*wght+t12*
     2         t26*t3+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -7.199999999999999d-5*t14**2*t6*wght-
     1         7.199999999999999d-5*gammaaa*t12*t28*t3*wght+Cmat2(iq,D2_
     2         GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = -7.199999999999999d-5*gammaaa*t12*t3
     1         *(1.2124444444444442d+0*t7*t8/rhoa**8.333333333333334d+0+
     2         2.4576d-5*gammaaa**4*t30/rhoa**1.3666666666666666d+1-4.56
     3         29629629629626d+1*gammaaa*t3/rhoa**5.666666666666667d+0-9
     4         .728000000000001d-3*t16*t17/rhoa**11)*wght-2.304d-6*t15*t
     5         21*t7*t8*wght+5.76d-6*t11*t7*t8*wght/rhoa**6+1.0666666666
     6         666665d-5*t5*wght/rhoa**1.6666666666666669d+0-3.199999999
     7         999999d-5*gammaaa*t11*t29*t3*wght+9.599999999999997d-5*ga
     8         mmaaa*t21*t23*t3*wght-3.6863999999999997d-8*t11*t16*t17*w
     9         ght/rhoa**8.666666666666666d+0+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t12*t3*(-7.199999999999999d-5*gamma
     1         aa*(-3.466666666666667d-1*gammaaa*t19*t8+3.264d-3*t17*t18
     2         *t7-9.216000000000001d-6*t16*t30/rhoa**1.2666666666666665
     3         d+1+9.777777777777779d+0*t20*t3)*wght-7.199999999999999d-
     4         5*t21*wght)+t15*t8*(-1.152d-6*t25*t7*wght-2.304d-6*gammaa
     5         a*t11*wght)+4.3199999999999994d-7*gammaaa*t21*t22*t8*wght
     6         +1.3824d-8*t11*t17*t7*wght/rhoa**7.666666666666667d+0-3.1
     7         99999999999999d-5*gammaaa*t14*t29*t3*wght+Cmat3(iq,D3_RA_
     8         RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t12*t3*(-7.199999999999999d-5*gamm
     1         aaa*(6.4d-2*t8*t9+3.456d-6*t30*t7/rhoa**1.166666666666666
     2         5d+1-9.600000000000001d-4*gammaaa*t17*t24)*wght-1.44d-4*t
     3         25*wght)+t23*t3*(-9.599999999999997d-5*gammaaa*t28*wght-9
     4         .599999999999997d-5*t14*wght)+t22*t8*(4.3199999999999994d
     5         -7*gammaaa*t25*wght+4.3199999999999994d-7*t11*wght-6.0d-3
     6         *t26)+5.759999999999998d-7*gammaaa*t14*t15*t8*wght-5.184d
     7         -9*gammaaa*t11*t17*wght/rhoa**6.666666666666667d+0+Cmat3(
     8         iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t12*t3*(-7.199999999999999d-5*gam
     1         maaa*(2.16d-4*t17*t27-1.2960000000000002d-6*gammaaa*t30/r
     2         hoa**1.0666666666666666d+1)*wght-7.199999999999999d-5*t28
     3         *wght)+4.3199999999999994d-7*gammaaa*t22*t28*t8*wght-1.44
     4         d-4*t14*t28*t6*wght+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 6.0d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = gammabb*t1*t3
            t5 = t4**2.0d+0
            t6 = rhob**1.3333333333333333d+0
            t7 = gammabb**2
            t8 = 1/t2**2
            t9 = 1/rhob**6.333333333333333d+0
            t10 = 1/rhob**3.6666666666666664d+0
            t11 = 1.6d-2*t7*t8*t9-2.6666666666666666d+0*gammabb*t10*t3
            t12 = 1/t6
            t13 = 1/rhob**5.333333333333333d+0
            t14 = t1*t3-6.0d-3*gammabb*t13*t8
            t15 = 1/rhob**5
            t16 = gammabb**3
            t17 = 1/t2**3
            t18 = 1/rhob**10
            t19 = 1/rhob**7.333333333333333d+0
            t20 = 1/rhob**4.666666666666667d+0
            t21 = -1.44d-1*t19*t7*t8+9.777777777777779d+0*gammabb*t20*t3
     1         +5.12d-4*t16*t17*t18
            t22 = 1/rhob**4
            t23 = 1/rhob**2.3333333333333334d+0
            t24 = 1/rhob**9
            t25 = 4.8d-2*gammabb*t8*t9-1.92d-4*t17*t24*t7-2.666666666666
     1         6666d+0*t10*t3
            t26 = -7.199999999999999d-5*gammabb*t25*wght-7.1999999999999
     1         99d-5*t11*wght
            t27 = 1/rhob**8
            t28 = 7.2d-5*gammabb*t17*t27-1.2d-2*t13*t8
            t29 = 1/rhob**3.3333333333333337d+0
            t30 = 1/t2**4
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RB) = -4.7999999999999987d-5*rhob**3.333333333333
     1         333d-1*t5*wght-7.199999999999999d-5*gammabb*t11*t12*t3*wg
     2         ht+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t14*t
     1         4*t6*wght
            Amat2(iq,D2_RB_RB) = -1.152d-6*t11*t15*t7*t8*wght-1.59999999
     1         99999995d-5*t5*wght/rhob**6.666666666666666d-1-7.19999999
     2         9999999d-5*gammabb*t12*t21*t3*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 4.3199999999999994d-7*gammabb*t11*t22*
     1         t8*wght-9.599999999999997d-5*gammabb*t14*t23*t3*wght+t12*
     2         t26*t3+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.199999999999999d-5*t14**2*t6*wght-
     1         7.199999999999999d-5*gammabb*t12*t28*t3*wght+Cmat2(iq,D2_
     2         GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = -7.199999999999999d-5*gammabb*t12*t3
     1         *(1.2124444444444442d+0*t7*t8/rhob**8.333333333333334d+0+
     2         2.4576d-5*gammabb**4*t30/rhob**1.3666666666666666d+1-4.56
     3         29629629629626d+1*gammabb*t3/rhob**5.666666666666667d+0-9
     4         .728000000000001d-3*t16*t17/rhob**11)*wght-2.304d-6*t15*t
     5         21*t7*t8*wght+5.76d-6*t11*t7*t8*wght/rhob**6+1.0666666666
     6         666665d-5*t5*wght/rhob**1.6666666666666669d+0-3.199999999
     7         999999d-5*gammabb*t11*t29*t3*wght+9.599999999999997d-5*ga
     8         mmabb*t21*t23*t3*wght-3.6863999999999997d-8*t11*t16*t17*w
     9         ght/rhob**8.666666666666666d+0+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t12*t3*(-7.199999999999999d-5*gamma
     1         bb*(-3.466666666666667d-1*gammabb*t19*t8+3.264d-3*t17*t18
     2         *t7-9.216000000000001d-6*t16*t30/rhob**1.2666666666666665
     3         d+1+9.777777777777779d+0*t20*t3)*wght-7.199999999999999d-
     4         5*t21*wght)+t15*t8*(-1.152d-6*t25*t7*wght-2.304d-6*gammab
     5         b*t11*wght)+4.3199999999999994d-7*gammabb*t21*t22*t8*wght
     6         +1.3824d-8*t11*t17*t7*wght/rhob**7.666666666666667d+0-3.1
     7         99999999999999d-5*gammabb*t14*t29*t3*wght+Cmat3(iq,D3_RB_
     8         RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t12*t3*(-7.199999999999999d-5*gamm
     1         abb*(6.4d-2*t8*t9+3.456d-6*t30*t7/rhob**1.166666666666666
     2         5d+1-9.600000000000001d-4*gammabb*t17*t24)*wght-1.44d-4*t
     3         25*wght)+t23*t3*(-9.599999999999997d-5*gammabb*t28*wght-9
     4         .599999999999997d-5*t14*wght)+t22*t8*(4.3199999999999994d
     5         -7*gammabb*t25*wght+4.3199999999999994d-7*t11*wght-6.0d-3
     6         *t26)+5.759999999999998d-7*gammabb*t14*t15*t8*wght-5.184d
     7         -9*gammabb*t11*t17*wght/rhob**6.666666666666667d+0+Cmat3(
     8         iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t12*t3*(-7.199999999999999d-5*gam
     1         mabb*(2.16d-4*t17*t27-1.2960000000000002d-6*gammabb*t30/r
     2         hob**1.0666666666666666d+1)*wght-7.199999999999999d-5*t28
     3         *wght)+4.3199999999999994d-7*gammabb*t22*t28*t8*wght-1.44
     4         d-4*t14*t28*t6*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
