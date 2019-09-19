C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_s12.F 
C> The nwxcm_x_s12 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_s12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)
C>    +1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}
C>    \over{{{{\it t_5}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{
C>    \it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)
C>    \,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}
C>    +1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{
C>    \it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_s12(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/t7**2
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t13*t8*t9
     1         /rhoa**2.3333333333333334d+0-1.2407009817988002d+0*rhoa**
     2         3.333333333333333d-1*t1*t11*t8-9.305257363491002d-1*t1*t1
     3         1*t12*t14*(-2.6666666666666666d+0*gammaaa*t5/rhoa**3.6666
     4         666666666664d+0-5.333333333333333d+0*t2*t3/rhoa**6.333333
     5         333333333d+0))*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t13*t8*t9/t12-9.
     1         305257363491002d-1*t1*t11*t12*t14*(t5*t6+2*gammaaa*t2*t4)
     2         )*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = gammabb**2
            t14 = 1/rhob**5.333333333333333d+0
            t15 = 1/rhob**2.6666666666666666d+0
            t16 = gammabb*t15*t5+t13*t14*t2+1.0d+0
            t17 = 1.0d+0-1.0d+0/t16
            t18 = gammabb*t15*t9+1.0d+0
            t19 = 1.0d+0-1.0d+0/t18
            t20 = rhob**1.3333333333333333d+0
            t21 = 1/t10**2
            t22 = 1/t7**2
            t23 = 1/t18**2
            t24 = 1/t16**2
            fnc(iq) = (-9.305257363491002d-1*t1*t11*t12*t8-9.30525736349
     1         1002d-1*t1*t17*t19*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t21*t8*t9
     1         /rhoa**2.3333333333333334d+0-1.2407009817988002d+0*rhoa**
     2         3.333333333333333d-1*t1*t11*t8-9.305257363491002d-1*t1*t1
     3         1*t12*t22*(-2.6666666666666666d+0*gammaaa*t5/rhoa**3.6666
     4         666666666664d+0-5.333333333333333d+0*t2*t3/rhoa**6.333333
     5         333333333d+0))*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (2.4814019635976003d+0*gammabb*t1*t17*t23*t
     1         9/rhob**2.3333333333333334d+0-9.305257363491002d-1*t1*t19
     2         *t20*t24*(-2.6666666666666666d+0*gammabb*t5/rhob**3.66666
     3         66666666664d+0-5.333333333333333d+0*t13*t2/rhob**6.333333
     4         333333333d+0)-1.2407009817988002d+0*rhob**3.3333333333333
     5         33d-1*t1*t17*t19)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t21*t8*t9/t12-9.
     1         305257363491002d-1*t1*t11*t12*t22*(t5*t6+2*gammaaa*t2*t4)
     2         )*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (-9.305257363491002d-1*t1*t17*t23*t9/t20-9
     1         .305257363491002d-1*t1*t19*t20*t24*(t15*t5+2*gammabb*t14*
     2         t2))*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/t7**2
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = 2.4814019635976003d+0*gammaaa*t1*t13*t8*t9*
     1         wght/rhoa**2.3333333333333334d+0-1.2407009817988002d+0*rh
     2         oa**3.333333333333333d-1*t1*t11*t8*wght-9.305257363491002
     3         d-1*t1*t11*t12*t14*(-2.6666666666666666d+0*gammaaa*t5/rho
     4         a**3.6666666666666664d+0-5.333333333333333d+0*t2*t3/rhoa*
     5         *6.333333333333333d+0)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -9.305257363491002d-1*t1*t13*t8*t9*wght/t1
     1         2-9.305257363491002d-1*t1*t11*t12*t14*(t5*t6+2*gammaaa*t2
     2         *t4)*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammabb**2
            t4 = 1/rhob**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhob**2.6666666666666666d+0
            t7 = gammabb*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammabb*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/t7**2
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RB) = 2.4814019635976003d+0*gammabb*t1*t13*t8*t9*
     1         wght/rhob**2.3333333333333334d+0-1.2407009817988002d+0*rh
     2         ob**3.333333333333333d-1*t1*t11*t8*wght-9.305257363491002
     3         d-1*t1*t11*t12*t14*(-2.6666666666666666d+0*gammabb*t5/rho
     4         b**3.6666666666666664d+0-5.333333333333333d+0*t2*t3/rhob*
     5         *6.333333333333333d+0)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -9.305257363491002d-1*t1*t13*t8*t9*wght/t1
     1         2-9.305257363491002d-1*t1*t11*t12*t14*(t5*t6+2*gammabb*t2
     2         *t4)*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_s12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)
C>    +1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}
C>    \over{{{{\it t_5}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{
C>    \it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)
C>    \,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}
C>    +1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{
C>    \it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_s12_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhoa**2.3333333333333334d+0
            t15 = rhoa**3.333333333333333d-1
            t16 = 1/rhoa**6.333333333333333d+0
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammaaa*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammaaa*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/t7**3
            t25 = 1/rhoa**4
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t13*t14*t
     1         8*t9-1.2407009817988002d+0*t1*t11*t15*t8-9.30525736349100
     2         2d-1*t1*t11*t12*t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t13*t20*t8*t9-9.
     1         305257363491002d-1*t1*t11*t12*t19*t21)*wght+Cmat(iq,D1_GA
     2         A)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-2.481401963597601d+0*gammaaa*t1*t13*t
     1         8*t9/rhoa**3.3333333333333337d+0+4.9628039271952007d+0*ga
     2         mmaaa*t1*t13*t14*t18*t19*t9+1.3234143805853868d+1*t1*t22*
     3         t23*t3*t8/rhoa**6-4.135669939329334d-1*t1*t11*t8/rhoa**6.
     4         666666666666666d-1+t19*(-9.305257363491002d-1*t1*t11*t12*
     5         (9.777777777777779d+0*gammaaa*t5/rhoa**4.666666666666667d
     6         +0+3.377777777777778d+1*t2*t3/rhoa**7.333333333333333d+0)
     7         -2.4814019635976003d+0*t1*t11*t15*t18)+1.8610514726982003
     8         d+0*t1*t11*t12*t18**2*t24)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (1.2407009817988002d+0*t1*t13*t14*t8*t
     1         9+2.4814019635976003d+0*gammaaa*t1*t13*t14*t19*t21*t9-9.3
     2         05257363491002d-1*t1*t13*t18*t19*t20*t9-4.962803927195200
     3         7d+0*gammaaa*t1*t22*t23*t8/rhoa**5+t19*(-9.30525736349100
     4         2d-1*t1*t11*t12*(-2.6666666666666666d+0*t17*t5-1.06666666
     5         66666666d+1*gammaaa*t16*t2)-1.2407009817988002d+0*t1*t11*
     6         t15*t21)+1.8610514726982003d+0*t1*t11*t12*t18*t21*t24)*wg
     7         ht+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (-1.8610514726982003d+0*t1*t13*t19*t2
     1         0*t21*t9+1.8610514726982003d+0*t1*t22*t23*t25*t8-1.861051
     2         4726982003d+0*t1*t11*t19*t2*t25+1.8610514726982003d+0*t1*
     3         t11*t12*t21**2*t24)*wght+Cmat2(iq,D2_GAA_GAA)
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
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = gammabb**2
            t14 = 1/rhob**5.333333333333333d+0
            t15 = 1/rhob**2.6666666666666666d+0
            t16 = gammabb*t15*t5+t13*t14*t2+1.0d+0
            t17 = 1.0d+0-1.0d+0/t16
            t18 = gammabb*t15*t9+1.0d+0
            t19 = 1.0d+0-1.0d+0/t18
            t20 = rhob**1.3333333333333333d+0
            t21 = 1/t10**2
            t22 = 1/rhoa**2.3333333333333334d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = -2.6666666666666666d+0*gammaaa*t25*t5-5.33333333333333
     1         3d+0*t2*t24*t3
            t27 = 1/t7**2
            t28 = 1/t18**2
            t29 = 1/rhob**2.3333333333333334d+0
            t30 = rhob**3.333333333333333d-1
            t31 = 1/rhob**6.333333333333333d+0
            t32 = 1/rhob**3.6666666666666664d+0
            t33 = -2.6666666666666666d+0*gammabb*t32*t5-5.33333333333333
     1         3d+0*t13*t2*t31
            t34 = 1/t16**2
            t35 = 1/t12
            t36 = t5*t6+2*gammaaa*t2*t4
            t37 = 1/t20
            t38 = t15*t5+2*gammabb*t14*t2
            t39 = t9**2
            t40 = 1/t10**3
            t41 = 1/t7**3
            t42 = 1/t18**3
            t43 = 1/t16**3
            t44 = 1/rhoa**4
            t45 = 1/rhob**4
            fnc(iq) = (-9.305257363491002d-1*t1*t11*t12*t8-9.30525736349
     1         1002d-1*t1*t17*t19*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t21*t22*t
     1         8*t9-1.2407009817988002d+0*t1*t11*t23*t8-9.30525736349100
     2         2d-1*t1*t11*t12*t26*t27)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (2.4814019635976003d+0*gammabb*t1*t17*t28*t
     1         29*t9-9.305257363491002d-1*t1*t19*t20*t33*t34-1.240700981
     2         7988002d+0*t1*t17*t19*t30)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t21*t35*t8*t9-9.
     1         305257363491002d-1*t1*t11*t12*t27*t36)*wght+Cmat(iq,D1_GA
     2         A)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (-9.305257363491002d-1*t1*t17*t28*t37*t9-9
     1         .305257363491002d-1*t1*t19*t20*t34*t38)*wght+Cmat(iq,D1_G
     2         BB)
            Amat2(iq,D2_RA_RA) = (-2.481401963597601d+0*gammaaa*t1*t21*t
     1         8*t9/rhoa**3.3333333333333337d+0+4.9628039271952007d+0*ga
     2         mmaaa*t1*t21*t22*t26*t27*t9+1.3234143805853868d+1*t1*t3*t
     3         39*t40*t8/rhoa**6-4.135669939329334d-1*t1*t11*t8/rhoa**6.
     4         666666666666666d-1+t27*(-9.305257363491002d-1*t1*t11*t12*
     5         (9.777777777777779d+0*gammaaa*t5/rhoa**4.666666666666667d
     6         +0+3.377777777777778d+1*t2*t3/rhoa**7.333333333333333d+0)
     7         -2.4814019635976003d+0*t1*t11*t23*t26)+1.8610514726982003
     8         d+0*t1*t11*t12*t26**2*t41)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (4.9628039271952007d+0*gammabb*t1*t28*t
     1         29*t33*t34*t9-2.481401963597601d+0*gammabb*t1*t17*t28*t9/
     2         rhob**3.3333333333333337d+0+t34*(-9.305257363491002d-1*t1
     3         *t19*t20*(9.777777777777779d+0*gammabb*t5/rhob**4.6666666
     4         66666667d+0+3.377777777777778d+1*t13*t2/rhob**7.333333333
     5         333333d+0)-2.4814019635976003d+0*t1*t19*t30*t33)+1.861051
     6         4726982003d+0*t1*t19*t20*t33**2*t43+1.3234143805853868d+1
     7         *t1*t13*t17*t39*t42/rhob**6-4.135669939329334d-1*t1*t17*t
     8         19/rhob**6.666666666666666d-1)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (1.2407009817988002d+0*t1*t21*t22*t8*t
     1         9+2.4814019635976003d+0*gammaaa*t1*t21*t22*t27*t36*t9-9.3
     2         05257363491002d-1*t1*t21*t26*t27*t35*t9-4.962803927195200
     3         7d+0*gammaaa*t1*t39*t40*t8/rhoa**5+t27*(-9.30525736349100
     4         2d-1*t1*t11*t12*(-2.6666666666666666d+0*t25*t5-1.06666666
     5         66666666d+1*gammaaa*t2*t24)-1.2407009817988002d+0*t1*t11*
     6         t23*t36)+1.8610514726982003d+0*t1*t11*t12*t26*t36*t41)*wg
     7         ht+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (2.4814019635976003d+0*gammabb*t1*t28*
     1         t29*t34*t38*t9-9.305257363491002d-1*t1*t28*t33*t34*t37*t9
     2         +1.2407009817988002d+0*t1*t17*t28*t29*t9+t34*(-9.30525736
     3         3491002d-1*t1*t19*t20*(-2.6666666666666666d+0*t32*t5-1.06
     4         66666666666666d+1*gammabb*t2*t31)-1.2407009817988002d+0*t
     5         1*t19*t30*t38)+1.8610514726982003d+0*t1*t19*t20*t33*t38*t
     6         43-4.9628039271952007d+0*gammabb*t1*t17*t39*t42/rhob**5)*
     7         wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (-1.8610514726982003d+0*t1*t21*t27*t3
     1         5*t36*t9+1.8610514726982003d+0*t1*t39*t40*t44*t8-1.861051
     2         4726982003d+0*t1*t11*t2*t27*t44+1.8610514726982003d+0*t1*
     3         t11*t12*t36**2*t41)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (-1.8610514726982003d+0*t1*t28*t34*t3
     1         7*t38*t9+1.8610514726982003d+0*t1*t17*t39*t42*t45-1.86105
     2         14726982003d+0*t1*t19*t2*t34*t45+1.8610514726982003d+0*t1
     3         *t19*t20*t38**2*t43)*wght+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhoa**2.3333333333333334d+0
            t15 = rhoa**3.333333333333333d-1
            t16 = 1/rhoa**6.333333333333333d+0
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammaaa*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammaaa*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/t7**3
            t25 = 1/rhoa**4
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = 2.4814019635976003d+0*gammaaa*t1*t13*t14*t8
     1         *t9*wght-1.2407009817988002d+0*t1*t11*t15*t8*wght-9.30525
     2         7363491002d-1*t1*t11*t12*t18*t19*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -9.305257363491002d-1*t1*t13*t20*t8*t9*wgh
     1         t-9.305257363491002d-1*t1*t11*t12*t19*t21*wght+Cmat(iq,D1
     2         _GAA)
            Amat2(iq,D2_RA_RA) = t19*(-9.305257363491002d-1*t1*t11*t12*(
     1         9.777777777777779d+0*gammaaa*t5/rhoa**4.666666666666667d+
     2         0+3.377777777777778d+1*t2*t3/rhoa**7.333333333333333d+0)*
     3         wght-2.4814019635976003d+0*t1*t11*t15*t18*wght)-2.4814019
     4         63597601d+0*gammaaa*t1*t13*t8*t9*wght/rhoa**3.33333333333
     5         33337d+0+4.9628039271952007d+0*gammaaa*t1*t13*t14*t18*t19
     6         *t9*wght+1.3234143805853868d+1*t1*t22*t23*t3*t8*wght/rhoa
     7         **6-4.135669939329334d-1*t1*t11*t8*wght/rhoa**6.666666666
     8         666666d-1+1.8610514726982003d+0*t1*t11*t12*t18**2*t24*wgh
     9         t+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = t19*(-9.305257363491002d-1*t1*t11*t12*
     1         (-2.6666666666666666d+0*t17*t5-1.0666666666666666d+1*gamm
     2         aaa*t16*t2)*wght-1.2407009817988002d+0*t1*t11*t15*t21*wgh
     3         t)+1.2407009817988002d+0*t1*t13*t14*t8*t9*wght+2.48140196
     4         35976003d+0*gammaaa*t1*t13*t14*t19*t21*t9*wght-9.30525736
     5         3491002d-1*t1*t13*t18*t19*t20*t9*wght-4.9628039271952007d
     6         +0*gammaaa*t1*t22*t23*t8*wght/rhoa**5+1.8610514726982003d
     7         +0*t1*t11*t12*t18*t21*t24*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -1.8610514726982003d+0*t1*t13*t19*t20
     1         *t21*t9*wght+1.8610514726982003d+0*t1*t22*t23*t25*t8*wght
     2         -1.8610514726982003d+0*t1*t11*t19*t2*t25*wght+1.861051472
     3         6982003d+0*t1*t11*t12*t21**2*t24*wght+Cmat2(iq,D2_GAA_GAA
     4         )
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammabb**2
            t4 = 1/rhob**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhob**2.6666666666666666d+0
            t7 = gammabb*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammabb*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhob**2.3333333333333334d+0
            t15 = rhob**3.333333333333333d-1
            t16 = 1/rhob**6.333333333333333d+0
            t17 = 1/rhob**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammabb*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammabb*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/t7**3
            t25 = 1/rhob**4
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RB) = 2.4814019635976003d+0*gammabb*t1*t13*t14*t8
     1         *t9*wght-1.2407009817988002d+0*t1*t11*t15*t8*wght-9.30525
     2         7363491002d-1*t1*t11*t12*t18*t19*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -9.305257363491002d-1*t1*t13*t20*t8*t9*wgh
     1         t-9.305257363491002d-1*t1*t11*t12*t19*t21*wght+Cmat(iq,D1
     2         _GBB)
            Amat2(iq,D2_RB_RB) = t19*(-9.305257363491002d-1*t1*t11*t12*(
     1         9.777777777777779d+0*gammabb*t5/rhob**4.666666666666667d+
     2         0+3.377777777777778d+1*t2*t3/rhob**7.333333333333333d+0)*
     3         wght-2.4814019635976003d+0*t1*t11*t15*t18*wght)-2.4814019
     4         63597601d+0*gammabb*t1*t13*t8*t9*wght/rhob**3.33333333333
     5         33337d+0+4.9628039271952007d+0*gammabb*t1*t13*t14*t18*t19
     6         *t9*wght+1.3234143805853868d+1*t1*t22*t23*t3*t8*wght/rhob
     7         **6-4.135669939329334d-1*t1*t11*t8*wght/rhob**6.666666666
     8         666666d-1+1.8610514726982003d+0*t1*t11*t12*t18**2*t24*wgh
     9         t+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = t19*(-9.305257363491002d-1*t1*t11*t12*
     1         (-2.6666666666666666d+0*t17*t5-1.0666666666666666d+1*gamm
     2         abb*t16*t2)*wght-1.2407009817988002d+0*t1*t11*t15*t21*wgh
     3         t)+1.2407009817988002d+0*t1*t13*t14*t8*t9*wght+2.48140196
     4         35976003d+0*gammabb*t1*t13*t14*t19*t21*t9*wght-9.30525736
     5         3491002d-1*t1*t13*t18*t19*t20*t9*wght-4.9628039271952007d
     6         +0*gammabb*t1*t22*t23*t8*wght/rhob**5+1.8610514726982003d
     7         +0*t1*t11*t12*t18*t21*t24*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -1.8610514726982003d+0*t1*t13*t19*t20
     1         *t21*t9*wght+1.8610514726982003d+0*t1*t22*t23*t25*t8*wght
     2         -1.8610514726982003d+0*t1*t11*t19*t2*t25*wght+1.861051472
     3         6982003d+0*t1*t11*t12*t21**2*t24*wght+Cmat2(iq,D2_GBB_GBB
     4         )
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_s12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)
C>    +1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}
C>    \over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}
C>    \over{{{{\it t_5}\,\sigma_{\beta\beta}^2}
C>    \over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}
C>    \,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{
C>    \it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)
C>    \,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}
C>    +{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}
C>    +1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{
C>    \it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="https://doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_s12_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhoa**2.3333333333333334d+0
            t15 = rhoa**3.333333333333333d-1
            t16 = 1/rhoa**6.333333333333333d+0
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammaaa*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammaaa*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/rhoa**6
            t25 = 1/rhoa**3.3333333333333337d+0
            t26 = 1/rhoa**6.666666666666666d-1
            t27 = t18**2
            t28 = 1/t7**3
            t29 = 1/rhoa**7.333333333333333d+0
            t30 = 1/rhoa**4.666666666666667d+0
            t31 = 9.777777777777779d+0*gammaaa*t30*t5+3.377777777777778d
     1         +1*t2*t29*t3
            t32 = -9.305257363491002d-1*t1*t11*t12*t31-2.481401963597600
     1         3d+0*t1*t11*t15*t18
            t33 = 1/rhoa**5
            t34 = -2.6666666666666666d+0*t17*t5-1.0666666666666666d+1*ga
     1         mmaaa*t16*t2
            t35 = -9.305257363491002d-1*t1*t11*t12*t34-1.240700981798800
     1         2d+0*t1*t11*t15*t21
            t36 = 1/rhoa**4
            t37 = t21**2
            t38 = t9**3
            t39 = 1/t10**4
            t40 = 1/t7**4
            t41 = 1/rhoa**7.666666666666667d+0
            t42 = 1/rhoa**6.666666666666667d+0
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t13*t14*t
     1         8*t9-1.2407009817988002d+0*t1*t11*t15*t8-9.30525736349100
     2         2d-1*t1*t11*t12*t18*t19)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t13*t20*t8*t9-9.
     1         305257363491002d-1*t1*t11*t12*t19*t21)*wght+Cmat(iq,D1_GA
     2         A)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-2.481401963597601d+0*gammaaa*t1*t13*t
     1         25*t8*t9+4.9628039271952007d+0*gammaaa*t1*t13*t14*t18*t19
     2         *t9+1.3234143805853868d+1*t1*t22*t23*t24*t3*t8-4.13566993
     3         9329334d-1*t1*t11*t26*t8+t19*t32+1.8610514726982003d+0*t1
     4         *t11*t12*t27*t28)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (1.2407009817988002d+0*t1*t13*t14*t8*t
     1         9+2.4814019635976003d+0*gammaaa*t1*t13*t14*t19*t21*t9-9.3
     2         05257363491002d-1*t1*t13*t18*t19*t20*t9-4.962803927195200
     3         7d+0*gammaaa*t1*t22*t23*t33*t8+t19*t35+1.8610514726982003
     4         d+0*t1*t11*t12*t18*t21*t28)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (-1.8610514726982003d+0*t1*t13*t19*t2
     1         0*t21*t9+1.8610514726982003d+0*t1*t22*t23*t36*t8+1.861051
     2         4726982003d+0*t1*t11*t12*t28*t37-1.8610514726982003d+0*t1
     3         *t11*t19*t2*t36)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (t19*(2.4814019635976003d+0*gammaaa*
     1         t1*t13*t14*t31*t9+6.617071902926934d+0*gammaaa*t1*t13*t18
     2         *t25*t9-9.305257363491002d-1*t1*t11*t12*(-4.5629629629629
     3         626d+1*gammaaa*t5/rhoa**5.666666666666667d+0-2.4770370370
     4         37037d+2*t2*t3/rhoa**8.333333333333334d+0)-3.722102945396
     5         4005d+0*t1*t11*t15*t31-8.271339878658668d-1*t1*t11*t18*t2
     6         6)+9.374185195813161d+0*gammaaa*t1*t13*t8*t9/rhoa**4.3333
     7         33333333333d+0+4.9628039271952007d+0*gammaaa*t1*t13*t14*t
     8         19*t31*t9-1.4888411781585603d+1*gammaaa*t1*t13*t14*t27*t2
     9         8*t9-1.4061277793719737d+1*gammaaa*t1*t13*t18*t19*t25*t9+
     :         1.0587315044683095d+2*gammaaa**3*t1*t38*t39*t8/rhoa**9.66
     ;         6666666666666d+0-9.263900664097708d+1*t1*t22*t23*t3*t8/rh
     <         oa**7+2.7571132928862224d-1*t1*t11*t8/rhoa**1.66666666666
     =         66669d+0-5.583154418094601d+0*t1*t11*t12*t18**3*t40+t28*(
     >         -2*t18*t32+3.7221029453964005d+0*t1*t11*t12*t18*t31+2.481
     ?         4019635976003d+0*t1*t11*t15*t27)+3.9702431417561607d+1*t1
     @         *t18*t19*t22*t23*t24*t3-4.135669939329334d-1*t1*t11*t18*t
     1         19*t26)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t13*t14*t19*(4.9628039271952007d+0
     1         *gammaaa*t1*t34*t9+4.9628039271952007d+0*t1*t18*t9)+t19*(
     2         -9.305257363491002d-1*t1*t13*t20*t31*t9-2.481401963597600
     3         3d+0*t1*t13*t14*t18*t9-9.305257363491002d-1*t1*t11*t12*(9
     4         .777777777777779d+0*t30*t5+6.755555555555556d+1*gammaaa*t
     5         2*t29)-2.4814019635976003d+0*t1*t11*t15*t34)-2.8949689575
     6         305343d+0*t1*t13*t25*t8*t9+1.8610514726982003d+0*t1*t13*t
     7         20*t27*t28*t9-9.925607854390401d+0*gammaaa*t1*t13*t14*t18
     8         *t21*t28*t9-2.481401963597601d+0*gammaaa*t1*t13*t19*t21*t
     9         25*t9-3.9702431417561607d+1*t1*t3*t38*t39*t8/rhoa**8.6666
     :         66666666666d+0+3.143109153890294d+1*gammaaa*t1*t22*t23*t2
     ;         4*t8-5.583154418094601d+0*t1*t11*t12*t21*t27*t40+t28*(3.7
     <         221029453964005d+0*t1*t11*t12*t18*t34-2*t21*t32)-9.925607
     =         854390401d+0*gammaaa*t1*t18*t19*t22*t23*t33+1.32341438058
     >         53868d+1*t1*t19*t21*t22*t23*t24*t3-4.135669939329334d-1*t
     ?         1*t11*t19*t21*t26)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t19*(-9.305257363491002d-1*t1*t13
     1         *t20*t34*t9-1.2407009817988002d+0*t1*t13*t14*t21*t9+7.444
     2         205890792801d+0*t1*t11*t2*t33)+4.9628039271952007d+0*gamm
     3         aaa*t1*t13*t19*t2*t41*t9-4.9628039271952007d+0*gammaaa*t1
     4         *t13*t14*t28*t37*t9-9.305257363491002d-1*t1*t13*t19*t20*t
     5         34*t9+3.7221029453964005d+0*t1*t13*t18*t20*t21*t28*t9+3.7
     6         221029453964005d+0*t1*t13*t14*t19*t21*t9+1.48884117815856
     7         03d+1*gammaaa*t1*t38*t39*t41*t8-7.444205890792801d+0*t1*t
     8         22*t23*t33*t8-5.583154418094601d+0*t1*t11*t12*t18*t37*t40
     9         +3.7221029453964005d+0*t1*t11*t18*t2*t28*t36+1.8610514726
     :         982003d+0*t1*t18*t19*t22*t23*t36+t28*(1.8610514726982003d
     ;         +0*t1*t11*t12*t21*t34-2*t21*t35)-9.925607854390401d+0*gam
     <         maaa*t1*t19*t21*t22*t23*t33)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (-5.583154418094601d+0*t1*t13*t19
     1         *t2*t42*t9+5.583154418094601d+0*t1*t13*t20*t28*t37*t9-5.5
     2         83154418094601d+0*t1*t38*t39*t42*t8-5.583154418094601d+0*
     3         t1*t11*t12*t21**3*t40+1.1166308836189202d+1*t1*t11*t2*t21
     4         *t28*t36+5.583154418094601d+0*t1*t19*t21*t22*t23*t36)*wgh
     5         t+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = gammabb**2
            t14 = 1/rhob**5.333333333333333d+0
            t15 = 1/rhob**2.6666666666666666d+0
            t16 = gammabb*t15*t5+t13*t14*t2+1.0d+0
            t17 = 1.0d+0-1.0d+0/t16
            t18 = gammabb*t15*t9+1.0d+0
            t19 = 1.0d+0-1.0d+0/t18
            t20 = rhob**1.3333333333333333d+0
            t21 = 1/t10**2
            t22 = 1/rhoa**2.3333333333333334d+0
            t23 = rhoa**3.333333333333333d-1
            t24 = 1/rhoa**6.333333333333333d+0
            t25 = 1/rhoa**3.6666666666666664d+0
            t26 = -2.6666666666666666d+0*gammaaa*t25*t5-5.33333333333333
     1         3d+0*t2*t24*t3
            t27 = 1/t7**2
            t28 = 1/t18**2
            t29 = 1/rhob**2.3333333333333334d+0
            t30 = rhob**3.333333333333333d-1
            t31 = 1/rhob**6.333333333333333d+0
            t32 = 1/rhob**3.6666666666666664d+0
            t33 = -2.6666666666666666d+0*gammabb*t32*t5-5.33333333333333
     1         3d+0*t13*t2*t31
            t34 = 1/t16**2
            t35 = 1/t12
            t36 = t5*t6+2*gammaaa*t2*t4
            t37 = 1/t20
            t38 = t15*t5+2*gammabb*t14*t2
            t39 = t9**2
            t40 = 1/t10**3
            t41 = 1/rhoa**6
            t42 = 1/rhoa**3.3333333333333337d+0
            t43 = 1/rhoa**6.666666666666666d-1
            t44 = t26**2
            t45 = 1/t7**3
            t46 = 1/rhoa**7.333333333333333d+0
            t47 = 1/rhoa**4.666666666666667d+0
            t48 = 9.777777777777779d+0*gammaaa*t47*t5+3.377777777777778d
     1         +1*t2*t3*t46
            t49 = -9.305257363491002d-1*t1*t11*t12*t48-2.481401963597600
     1         3d+0*t1*t11*t23*t26
            t50 = 1/t18**3
            t51 = 1/rhob**6
            t52 = 1/rhob**3.3333333333333337d+0
            t53 = 1/rhob**6.666666666666666d-1
            t54 = t33**2
            t55 = 1/t16**3
            t56 = 1/rhob**7.333333333333333d+0
            t57 = 1/rhob**4.666666666666667d+0
            t58 = 9.777777777777779d+0*gammabb*t5*t57+3.377777777777778d
     1         +1*t13*t2*t56
            t59 = -9.305257363491002d-1*t1*t19*t20*t58-2.481401963597600
     1         3d+0*t1*t19*t30*t33
            t60 = 1/rhoa**5
            t61 = -2.6666666666666666d+0*t25*t5-1.0666666666666666d+1*ga
     1         mmaaa*t2*t24
            t62 = -9.305257363491002d-1*t1*t11*t12*t61-1.240700981798800
     1         2d+0*t1*t11*t23*t36
            t63 = 1/rhob**5
            t64 = -2.6666666666666666d+0*t32*t5-1.0666666666666666d+1*ga
     1         mmabb*t2*t31
            t65 = -9.305257363491002d-1*t1*t19*t20*t64-1.240700981798800
     1         2d+0*t1*t19*t30*t38
            t66 = 1/rhoa**4
            t67 = t36**2
            t68 = 1/rhob**4
            t69 = t38**2
            t70 = t9**3
            t71 = 1/t10**4
            t72 = 1/t7**4
            t73 = 1/t18**4
            t74 = 1/t16**4
            t75 = 1/rhoa**7.666666666666667d+0
            t76 = 1/rhob**7.666666666666667d+0
            t77 = 1/rhoa**6.666666666666667d+0
            t78 = 1/rhob**6.666666666666667d+0
            fnc(iq) = (-9.305257363491002d-1*t1*t11*t12*t8-9.30525736349
     1         1002d-1*t1*t17*t19*t20)*wght+fnc(iq)
            Amat(iq,D1_RA) = (2.4814019635976003d+0*gammaaa*t1*t21*t22*t
     1         8*t9-1.2407009817988002d+0*t1*t11*t23*t8-9.30525736349100
     2         2d-1*t1*t11*t12*t26*t27)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (2.4814019635976003d+0*gammabb*t1*t17*t28*t
     1         29*t9-9.305257363491002d-1*t1*t19*t20*t33*t34-1.240700981
     2         7988002d+0*t1*t17*t19*t30)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (-9.305257363491002d-1*t1*t21*t35*t8*t9-9.
     1         305257363491002d-1*t1*t11*t12*t27*t36)*wght+Cmat(iq,D1_GA
     2         A)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (-9.305257363491002d-1*t1*t17*t28*t37*t9-9
     1         .305257363491002d-1*t1*t19*t20*t34*t38)*wght+Cmat(iq,D1_G
     2         BB)
            Amat2(iq,D2_RA_RA) = (-2.481401963597601d+0*gammaaa*t1*t21*t
     1         42*t8*t9+4.9628039271952007d+0*gammaaa*t1*t21*t22*t26*t27
     2         *t9-4.135669939329334d-1*t1*t11*t43*t8+1.3234143805853868
     3         d+1*t1*t3*t39*t40*t41*t8+t27*t49+1.8610514726982003d+0*t1
     4         *t11*t12*t44*t45)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-2.481401963597601d+0*gammabb*t1*t17*t
     1         28*t52*t9+4.9628039271952007d+0*gammabb*t1*t28*t29*t33*t3
     2         4*t9+t34*t59+1.8610514726982003d+0*t1*t19*t20*t54*t55-4.1
     3         35669939329334d-1*t1*t17*t19*t53+1.3234143805853868d+1*t1
     4         *t13*t17*t39*t50*t51)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (1.2407009817988002d+0*t1*t21*t22*t8*t
     1         9+2.4814019635976003d+0*gammaaa*t1*t21*t22*t27*t36*t9-9.3
     2         05257363491002d-1*t1*t21*t26*t27*t35*t9-4.962803927195200
     3         7d+0*gammaaa*t1*t39*t40*t60*t8+t27*t62+1.8610514726982003
     4         d+0*t1*t11*t12*t26*t36*t45)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (2.4814019635976003d+0*gammabb*t1*t28*
     1         t29*t34*t38*t9-9.305257363491002d-1*t1*t28*t33*t34*t37*t9
     2         +1.2407009817988002d+0*t1*t17*t28*t29*t9+t34*t65-4.962803
     3         9271952007d+0*gammabb*t1*t17*t39*t50*t63+1.86105147269820
     4         03d+0*t1*t19*t20*t33*t38*t55)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (-1.8610514726982003d+0*t1*t21*t27*t3
     1         5*t36*t9+1.8610514726982003d+0*t1*t39*t40*t66*t8+1.861051
     2         4726982003d+0*t1*t11*t12*t45*t67-1.8610514726982003d+0*t1
     3         *t11*t2*t27*t66)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (-1.8610514726982003d+0*t1*t28*t34*t3
     1         7*t38*t9+1.8610514726982003d+0*t1*t19*t20*t55*t69+1.86105
     2         14726982003d+0*t1*t17*t39*t50*t68-1.8610514726982003d+0*t
     3         1*t19*t2*t34*t68)*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (t27*(2.4814019635976003d+0*gammaaa*
     1         t1*t21*t22*t48*t9+6.617071902926934d+0*gammaaa*t1*t21*t26
     2         *t42*t9-9.305257363491002d-1*t1*t11*t12*(-4.5629629629629
     3         626d+1*gammaaa*t5/rhoa**5.666666666666667d+0-2.4770370370
     4         37037d+2*t2*t3/rhoa**8.333333333333334d+0)-3.722102945396
     5         4005d+0*t1*t11*t23*t48-8.271339878658668d-1*t1*t11*t26*t4
     6         3)+9.374185195813161d+0*gammaaa*t1*t21*t8*t9/rhoa**4.3333
     7         33333333333d+0+4.9628039271952007d+0*gammaaa*t1*t21*t22*t
     8         27*t48*t9-1.4888411781585603d+1*gammaaa*t1*t21*t22*t44*t4
     9         5*t9-1.4061277793719737d+1*gammaaa*t1*t21*t26*t27*t42*t9+
     :         1.0587315044683095d+2*gammaaa**3*t1*t70*t71*t8/rhoa**9.66
     ;         6666666666666d+0-9.263900664097708d+1*t1*t3*t39*t40*t8/rh
     <         oa**7+2.7571132928862224d-1*t1*t11*t8/rhoa**1.66666666666
     =         66669d+0-5.583154418094601d+0*t1*t11*t12*t26**3*t72+t45*(
     >         -2*t26*t49+3.7221029453964005d+0*t1*t11*t12*t26*t48+2.481
     ?         4019635976003d+0*t1*t11*t23*t44)-4.135669939329334d-1*t1*
     @         t11*t26*t27*t43+3.9702431417561607d+1*t1*t26*t27*t3*t39*t
     1         40*t41)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t34*(2.4814019635976003d+0*gammabb*
     1         t1*t28*t29*t58*t9+6.617071902926934d+0*gammabb*t1*t28*t33
     2         *t52*t9-3.7221029453964005d+0*t1*t19*t30*t58-8.2713398786
     3         58668d-1*t1*t19*t33*t53-9.305257363491002d-1*t1*t19*t20*(
     4         -4.5629629629629626d+1*gammabb*t5/rhob**5.666666666666667
     5         d+0-2.477037037037037d+2*t13*t2/rhob**8.333333333333334d+
     6         0))+4.9628039271952007d+0*gammabb*t1*t28*t29*t34*t58*t9-1
     7         .4888411781585603d+1*gammabb*t1*t28*t29*t54*t55*t9-1.4061
     8         277793719737d+1*gammabb*t1*t28*t33*t34*t52*t9+9.374185195
     9         813161d+0*gammabb*t1*t17*t28*t9/rhob**4.333333333333333d+
     :         0-5.583154418094601d+0*t1*t19*t20*t33**3*t74+1.0587315044
     ;         683095d+2*gammabb**3*t1*t17*t70*t73/rhob**9.6666666666666
     <         66d+0+t55*(-2*t33*t59+3.7221029453964005d+0*t1*t19*t20*t3
     =         3*t58+2.4814019635976003d+0*t1*t19*t30*t54)-4.13566993932
     >         9334d-1*t1*t19*t33*t34*t53+3.9702431417561607d+1*t1*t13*t
     ?         33*t34*t39*t50*t51-9.263900664097708d+1*t1*t13*t17*t39*t5
     @         0/rhob**7+2.7571132928862224d-1*t1*t17*t19/rhob**1.666666
     1         6666666669d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t21*t22*t27*(4.9628039271952007d+0
     1         *gammaaa*t1*t61*t9+4.9628039271952007d+0*t1*t26*t9)+t27*(
     2         -9.305257363491002d-1*t1*t21*t35*t48*t9-2.481401963597600
     3         3d+0*t1*t21*t22*t26*t9-2.4814019635976003d+0*t1*t11*t23*t
     4         61-9.305257363491002d-1*t1*t11*t12*(9.777777777777779d+0*
     5         t47*t5+6.755555555555556d+1*gammaaa*t2*t46))-2.8949689575
     6         305343d+0*t1*t21*t42*t8*t9+1.8610514726982003d+0*t1*t21*t
     7         35*t44*t45*t9-9.925607854390401d+0*gammaaa*t1*t21*t22*t26
     8         *t36*t45*t9-2.481401963597601d+0*gammaaa*t1*t21*t27*t36*t
     9         42*t9-3.9702431417561607d+1*t1*t3*t70*t71*t8/rhoa**8.6666
     :         66666666666d+0+3.143109153890294d+1*gammaaa*t1*t39*t40*t4
     ;         1*t8-5.583154418094601d+0*t1*t11*t12*t36*t44*t72+t45*(3.7
     <         221029453964005d+0*t1*t11*t12*t26*t61-2*t36*t49)-9.925607
     =         854390401d+0*gammaaa*t1*t26*t27*t39*t40*t60-4.13566993932
     >         9334d-1*t1*t11*t27*t36*t43+1.3234143805853868d+1*t1*t27*t
     ?         3*t36*t39*t40*t41)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (t28*t29*t34*(4.9628039271952007d+0
     1         *gammabb*t1*t64*t9+4.9628039271952007d+0*t1*t33*t9)+t34*(
     2         -9.305257363491002d-1*t1*t28*t37*t58*t9-2.481401963597600
     3         3d+0*t1*t28*t29*t33*t9-2.4814019635976003d+0*t1*t19*t30*t
     4         64-9.305257363491002d-1*t1*t19*t20*(9.777777777777779d+0*
     5         t5*t57+6.755555555555556d+1*gammabb*t2*t56))+1.8610514726
     6         982003d+0*t1*t28*t37*t54*t55*t9-9.925607854390401d+0*gamm
     7         abb*t1*t28*t29*t33*t38*t55*t9-2.481401963597601d+0*gammab
     8         b*t1*t28*t34*t38*t52*t9-2.8949689575305343d+0*t1*t17*t28*
     9         t52*t9-5.583154418094601d+0*t1*t19*t20*t38*t54*t74-3.9702
     :         431417561607d+1*t1*t13*t17*t70*t73/rhob**8.66666666666666
     ;         6d+0+t55*(3.7221029453964005d+0*t1*t19*t20*t33*t64-2*t38*
     <         t59)-9.925607854390401d+0*gammabb*t1*t33*t34*t39*t50*t63-
     =         4.135669939329334d-1*t1*t19*t34*t38*t53+1.323414380585386
     >         8d+1*t1*t13*t34*t38*t39*t50*t51+3.143109153890294d+1*gamm
     ?         abb*t1*t17*t39*t50*t51)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (t27*(-9.305257363491002d-1*t1*t21
     1         *t35*t61*t9-1.2407009817988002d+0*t1*t21*t22*t36*t9+7.444
     2         205890792801d+0*t1*t11*t2*t60)+4.9628039271952007d+0*gamm
     3         aaa*t1*t2*t21*t27*t75*t9-4.9628039271952007d+0*gammaaa*t1
     4         *t21*t22*t45*t67*t9-9.305257363491002d-1*t1*t21*t27*t35*t
     5         61*t9+3.7221029453964005d+0*t1*t21*t26*t35*t36*t45*t9+3.7
     6         221029453964005d+0*t1*t21*t22*t27*t36*t9+1.48884117815856
     7         03d+1*gammaaa*t1*t70*t71*t75*t8-7.444205890792801d+0*t1*t
     8         39*t40*t60*t8-5.583154418094601d+0*t1*t11*t12*t26*t67*t72
     9         +3.7221029453964005d+0*t1*t11*t2*t26*t45*t66+1.8610514726
     :         982003d+0*t1*t26*t27*t39*t40*t66+t45*(1.8610514726982003d
     ;         +0*t1*t11*t12*t36*t61-2*t36*t62)-9.925607854390401d+0*gam
     <         maaa*t1*t27*t36*t39*t40*t60)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (t34*(-9.305257363491002d-1*t1*t28
     1         *t37*t64*t9-1.2407009817988002d+0*t1*t28*t29*t38*t9+7.444
     2         205890792801d+0*t1*t19*t2*t63)+4.9628039271952007d+0*gamm
     3         abb*t1*t2*t28*t34*t76*t9-4.9628039271952007d+0*gammabb*t1
     4         *t28*t29*t55*t69*t9-9.305257363491002d-1*t1*t28*t34*t37*t
     5         64*t9+3.7221029453964005d+0*t1*t28*t33*t37*t38*t55*t9+3.7
     6         221029453964005d+0*t1*t28*t29*t34*t38*t9+1.48884117815856
     7         03d+1*gammabb*t1*t17*t70*t73*t76-5.583154418094601d+0*t1*
     8         t19*t20*t33*t69*t74+3.7221029453964005d+0*t1*t19*t2*t33*t
     9         55*t68+1.8610514726982003d+0*t1*t33*t34*t39*t50*t68+t55*(
     :         1.8610514726982003d+0*t1*t19*t20*t38*t64-2*t38*t65)-9.925
     ;         607854390401d+0*gammabb*t1*t34*t38*t39*t50*t63-7.44420589
     <         0792801d+0*t1*t17*t39*t50*t63)*wght+Cmat3(iq,D3_RB_GBB_GB
     =         B)
            Cmat3(iq,D3_GAA_GAA_GAA) = (-5.583154418094601d+0*t1*t2*t21*
     1         t27*t77*t9+5.583154418094601d+0*t1*t21*t35*t45*t67*t9-5.5
     2         83154418094601d+0*t1*t70*t71*t77*t8-5.583154418094601d+0*
     3         t1*t11*t12*t36**3*t72+1.1166308836189202d+1*t1*t11*t2*t36
     4         *t45*t66+5.583154418094601d+0*t1*t27*t36*t39*t40*t66)*wgh
     5         t+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (-5.583154418094601d+0*t1*t2*t28*
     1         t34*t78*t9+5.583154418094601d+0*t1*t28*t37*t55*t69*t9-5.5
     2         83154418094601d+0*t1*t17*t70*t73*t78-5.583154418094601d+0
     3         *t1*t19*t20*t38**3*t74+1.1166308836189202d+1*t1*t19*t2*t3
     4         8*t55*t68+5.583154418094601d+0*t1*t34*t38*t39*t50*t68)*wg
     5         ht+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhoa**2.6666666666666666d+0
            t7 = gammaaa*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammaaa*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhoa**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhoa**2.3333333333333334d+0
            t15 = rhoa**3.333333333333333d-1
            t16 = 1/rhoa**6.333333333333333d+0
            t17 = 1/rhoa**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammaaa*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammaaa*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/rhoa**6
            t25 = 1/rhoa**3.3333333333333337d+0
            t26 = 1/rhoa**6.666666666666666d-1
            t27 = t18**2
            t28 = 1/t7**3
            t29 = 1/rhoa**7.333333333333333d+0
            t30 = 1/rhoa**4.666666666666667d+0
            t31 = 9.777777777777779d+0*gammaaa*t30*t5+3.377777777777778d
     1         +1*t2*t29*t3
            t32 = -9.305257363491002d-1*t1*t11*t12*t31*wght-2.4814019635
     1         976003d+0*t1*t11*t15*t18*wght
            t33 = 1/rhoa**5
            t34 = -2.6666666666666666d+0*t17*t5-1.0666666666666666d+1*ga
     1         mmaaa*t16*t2
            t35 = -9.305257363491002d-1*t1*t11*t12*t34*wght-1.2407009817
     1         988002d+0*t1*t11*t15*t21*wght
            t36 = 1/rhoa**4
            t37 = t21**2
            t38 = t9**3
            t39 = 1/t10**4
            t40 = 1/t7**4
            t41 = 1/rhoa**7.666666666666667d+0
            t42 = 1/rhoa**6.666666666666667d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = 2.4814019635976003d+0*gammaaa*t1*t13*t14*t8
     1         *t9*wght-1.2407009817988002d+0*t1*t11*t15*t8*wght-9.30525
     2         7363491002d-1*t1*t11*t12*t18*t19*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = -9.305257363491002d-1*t1*t13*t20*t8*t9*wgh
     1         t-9.305257363491002d-1*t1*t11*t12*t19*t21*wght+Cmat(iq,D1
     2         _GAA)
            Amat2(iq,D2_RA_RA) = -2.481401963597601d+0*gammaaa*t1*t13*t2
     1         5*t8*t9*wght+4.9628039271952007d+0*gammaaa*t1*t13*t14*t18
     2         *t19*t9*wght+1.3234143805853868d+1*t1*t22*t23*t24*t3*t8*w
     3         ght-4.135669939329334d-1*t1*t11*t26*t8*wght+1.86105147269
     4         82003d+0*t1*t11*t12*t27*t28*wght+t19*t32+Amat2(iq,D2_RA_R
     5         A)
            Cmat2(iq,D2_RA_GAA) = 1.2407009817988002d+0*t1*t13*t14*t8*t9
     1         *wght+2.4814019635976003d+0*gammaaa*t1*t13*t14*t19*t21*t9
     2         *wght-9.305257363491002d-1*t1*t13*t18*t19*t20*t9*wght-4.9
     3         628039271952007d+0*gammaaa*t1*t22*t23*t33*t8*wght+1.86105
     4         14726982003d+0*t1*t11*t12*t18*t21*t28*wght+t19*t35+Cmat2(
     5         iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -1.8610514726982003d+0*t1*t13*t19*t20
     1         *t21*t9*wght+1.8610514726982003d+0*t1*t22*t23*t36*t8*wght
     2         +1.8610514726982003d+0*t1*t11*t12*t28*t37*wght-1.86105147
     3         26982003d+0*t1*t11*t19*t2*t36*wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = t19*(2.4814019635976003d+0*gammaaa*t
     1         1*t13*t14*t31*t9*wght+6.617071902926934d+0*gammaaa*t1*t13
     2         *t18*t25*t9*wght-9.305257363491002d-1*t1*t11*t12*(-4.5629
     3         629629629626d+1*gammaaa*t5/rhoa**5.666666666666667d+0-2.4
     4         77037037037037d+2*t2*t3/rhoa**8.333333333333334d+0)*wght-
     5         3.7221029453964005d+0*t1*t11*t15*t31*wght-8.2713398786586
     6         68d-1*t1*t11*t18*t26*wght)+t28*(3.7221029453964005d+0*t1*
     7         t11*t12*t18*t31*wght+2.4814019635976003d+0*t1*t11*t15*t27
     8         *wght-2*t18*t32)+9.374185195813161d+0*gammaaa*t1*t13*t8*t
     9         9*wght/rhoa**4.333333333333333d+0+4.9628039271952007d+0*g
     :         ammaaa*t1*t13*t14*t19*t31*t9*wght-1.4888411781585603d+1*g
     ;         ammaaa*t1*t13*t14*t27*t28*t9*wght-1.4061277793719737d+1*g
     <         ammaaa*t1*t13*t18*t19*t25*t9*wght+1.0587315044683095d+2*g
     =         ammaaa**3*t1*t38*t39*t8*wght/rhoa**9.666666666666666d+0-9
     >         .263900664097708d+1*t1*t22*t23*t3*t8*wght/rhoa**7+2.75711
     ?         32928862224d-1*t1*t11*t8*wght/rhoa**1.6666666666666669d+0
     @         -5.583154418094601d+0*t1*t11*t12*t18**3*t40*wght+3.970243
     1         1417561607d+1*t1*t18*t19*t22*t23*t24*t3*wght-4.1356699393
     2         29334d-1*t1*t11*t18*t19*t26*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = t13*t14*t19*(4.9628039271952007d+0*
     1         gammaaa*t1*t34*t9*wght+4.9628039271952007d+0*t1*t18*t9*wg
     2         ht)+t19*(-9.305257363491002d-1*t1*t13*t20*t31*t9*wght-2.4
     3         814019635976003d+0*t1*t13*t14*t18*t9*wght-9.3052573634910
     4         02d-1*t1*t11*t12*(9.777777777777779d+0*t30*t5+6.755555555
     5         555556d+1*gammaaa*t2*t29)*wght-2.4814019635976003d+0*t1*t
     6         11*t15*t34*wght)+t28*(3.7221029453964005d+0*t1*t11*t12*t1
     7         8*t34*wght-2*t21*t32)-2.8949689575305343d+0*t1*t13*t25*t8
     8         *t9*wght+1.8610514726982003d+0*t1*t13*t20*t27*t28*t9*wght
     9         -9.925607854390401d+0*gammaaa*t1*t13*t14*t18*t21*t28*t9*w
     :         ght-2.481401963597601d+0*gammaaa*t1*t13*t19*t21*t25*t9*wg
     ;         ht-3.9702431417561607d+1*t1*t3*t38*t39*t8*wght/rhoa**8.66
     <         6666666666666d+0+3.143109153890294d+1*gammaaa*t1*t22*t23*
     =         t24*t8*wght-5.583154418094601d+0*t1*t11*t12*t21*t27*t40*w
     >         ght-9.925607854390401d+0*gammaaa*t1*t18*t19*t22*t23*t33*w
     ?         ght+1.3234143805853868d+1*t1*t19*t21*t22*t23*t24*t3*wght-
     @         4.135669939329334d-1*t1*t11*t19*t21*t26*wght+Cmat3(iq,D3_
     1         RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = t19*(-9.305257363491002d-1*t1*t13*
     1         t20*t34*t9*wght-1.2407009817988002d+0*t1*t13*t14*t21*t9*w
     2         ght+7.444205890792801d+0*t1*t11*t2*t33*wght)+t28*(1.86105
     3         14726982003d+0*t1*t11*t12*t21*t34*wght-2*t21*t35)+4.96280
     4         39271952007d+0*gammaaa*t1*t13*t19*t2*t41*t9*wght-4.962803
     5         9271952007d+0*gammaaa*t1*t13*t14*t28*t37*t9*wght-9.305257
     6         363491002d-1*t1*t13*t19*t20*t34*t9*wght+3.722102945396400
     7         5d+0*t1*t13*t18*t20*t21*t28*t9*wght+3.7221029453964005d+0
     8         *t1*t13*t14*t19*t21*t9*wght+1.4888411781585603d+1*gammaaa
     9         *t1*t38*t39*t41*t8*wght-7.444205890792801d+0*t1*t22*t23*t
     :         33*t8*wght-5.583154418094601d+0*t1*t11*t12*t18*t37*t40*wg
     ;         ht+3.7221029453964005d+0*t1*t11*t18*t2*t28*t36*wght+1.861
     <         0514726982003d+0*t1*t18*t19*t22*t23*t36*wght-9.9256078543
     =         90401d+0*gammaaa*t1*t19*t21*t22*t23*t33*wght+Cmat3(iq,D3_
     >         RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = -5.583154418094601d+0*t1*t13*t19*
     1         t2*t42*t9*wght+5.583154418094601d+0*t1*t13*t20*t28*t37*t9
     2         *wght-5.583154418094601d+0*t1*t38*t39*t42*t8*wght-5.58315
     3         4418094601d+0*t1*t11*t12*t21**3*t40*wght+1.11663088361892
     4         02d+1*t1*t11*t2*t21*t28*t36*wght+5.583154418094601d+0*t1*
     5         t19*t21*t22*t23*t36*wght+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)-param(1)+1.0d+0
            t2 = param(4)
            t3 = gammabb**2
            t4 = 1/rhob**5.333333333333333d+0
            t5 = param(3)
            t6 = 1/rhob**2.6666666666666666d+0
            t7 = gammabb*t5*t6+t2*t3*t4+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = param(5)
            t10 = gammabb*t6*t9+1.0d+0
            t11 = 1.0d+0-1.0d+0/t10
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t10**2
            t14 = 1/rhob**2.3333333333333334d+0
            t15 = rhob**3.333333333333333d-1
            t16 = 1/rhob**6.333333333333333d+0
            t17 = 1/rhob**3.6666666666666664d+0
            t18 = -2.6666666666666666d+0*gammabb*t17*t5-5.33333333333333
     1         3d+0*t16*t2*t3
            t19 = 1/t7**2
            t20 = 1/t12
            t21 = t5*t6+2*gammabb*t2*t4
            t22 = t9**2
            t23 = 1/t10**3
            t24 = 1/rhob**6
            t25 = 1/rhob**3.3333333333333337d+0
            t26 = 1/rhob**6.666666666666666d-1
            t27 = t18**2
            t28 = 1/t7**3
            t29 = 1/rhob**7.333333333333333d+0
            t30 = 1/rhob**4.666666666666667d+0
            t31 = 9.777777777777779d+0*gammabb*t30*t5+3.377777777777778d
     1         +1*t2*t29*t3
            t32 = -9.305257363491002d-1*t1*t11*t12*t31*wght-2.4814019635
     1         976003d+0*t1*t11*t15*t18*wght
            t33 = 1/rhob**5
            t34 = -2.6666666666666666d+0*t17*t5-1.0666666666666666d+1*ga
     1         mmabb*t16*t2
            t35 = -9.305257363491002d-1*t1*t11*t12*t34*wght-1.2407009817
     1         988002d+0*t1*t11*t15*t21*wght
            t36 = 1/rhob**4
            t37 = t21**2
            t38 = t9**3
            t39 = 1/t10**4
            t40 = 1/t7**4
            t41 = 1/rhob**7.666666666666667d+0
            t42 = 1/rhob**6.666666666666667d+0
            fnc(iq) = fnc(iq)-9.305257363491002d-1*t1*t11*t12*t8*wght
            Amat(iq,D1_RB) = 2.4814019635976003d+0*gammabb*t1*t13*t14*t8
     1         *t9*wght-1.2407009817988002d+0*t1*t11*t15*t8*wght-9.30525
     2         7363491002d-1*t1*t11*t12*t18*t19*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = -9.305257363491002d-1*t1*t13*t20*t8*t9*wgh
     1         t-9.305257363491002d-1*t1*t11*t12*t19*t21*wght+Cmat(iq,D1
     2         _GBB)
            Amat2(iq,D2_RB_RB) = -2.481401963597601d+0*gammabb*t1*t13*t2
     1         5*t8*t9*wght+4.9628039271952007d+0*gammabb*t1*t13*t14*t18
     2         *t19*t9*wght+1.3234143805853868d+1*t1*t22*t23*t24*t3*t8*w
     3         ght-4.135669939329334d-1*t1*t11*t26*t8*wght+1.86105147269
     4         82003d+0*t1*t11*t12*t27*t28*wght+t19*t32+Amat2(iq,D2_RB_R
     5         B)
            Cmat2(iq,D2_RB_GBB) = 1.2407009817988002d+0*t1*t13*t14*t8*t9
     1         *wght+2.4814019635976003d+0*gammabb*t1*t13*t14*t19*t21*t9
     2         *wght-9.305257363491002d-1*t1*t13*t18*t19*t20*t9*wght-4.9
     3         628039271952007d+0*gammabb*t1*t22*t23*t33*t8*wght+1.86105
     4         14726982003d+0*t1*t11*t12*t18*t21*t28*wght+t19*t35+Cmat2(
     5         iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -1.8610514726982003d+0*t1*t13*t19*t20
     1         *t21*t9*wght+1.8610514726982003d+0*t1*t22*t23*t36*t8*wght
     2         +1.8610514726982003d+0*t1*t11*t12*t28*t37*wght-1.86105147
     3         26982003d+0*t1*t11*t19*t2*t36*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = t19*(2.4814019635976003d+0*gammabb*t
     1         1*t13*t14*t31*t9*wght+6.617071902926934d+0*gammabb*t1*t13
     2         *t18*t25*t9*wght-9.305257363491002d-1*t1*t11*t12*(-4.5629
     3         629629629626d+1*gammabb*t5/rhob**5.666666666666667d+0-2.4
     4         77037037037037d+2*t2*t3/rhob**8.333333333333334d+0)*wght-
     5         3.7221029453964005d+0*t1*t11*t15*t31*wght-8.2713398786586
     6         68d-1*t1*t11*t18*t26*wght)+t28*(3.7221029453964005d+0*t1*
     7         t11*t12*t18*t31*wght+2.4814019635976003d+0*t1*t11*t15*t27
     8         *wght-2*t18*t32)+9.374185195813161d+0*gammabb*t1*t13*t8*t
     9         9*wght/rhob**4.333333333333333d+0+4.9628039271952007d+0*g
     :         ammabb*t1*t13*t14*t19*t31*t9*wght-1.4888411781585603d+1*g
     ;         ammabb*t1*t13*t14*t27*t28*t9*wght-1.4061277793719737d+1*g
     <         ammabb*t1*t13*t18*t19*t25*t9*wght+1.0587315044683095d+2*g
     =         ammabb**3*t1*t38*t39*t8*wght/rhob**9.666666666666666d+0-9
     >         .263900664097708d+1*t1*t22*t23*t3*t8*wght/rhob**7+2.75711
     ?         32928862224d-1*t1*t11*t8*wght/rhob**1.6666666666666669d+0
     @         -5.583154418094601d+0*t1*t11*t12*t18**3*t40*wght+3.970243
     1         1417561607d+1*t1*t18*t19*t22*t23*t24*t3*wght-4.1356699393
     2         29334d-1*t1*t11*t18*t19*t26*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = t13*t14*t19*(4.9628039271952007d+0*
     1         gammabb*t1*t34*t9*wght+4.9628039271952007d+0*t1*t18*t9*wg
     2         ht)+t19*(-9.305257363491002d-1*t1*t13*t20*t31*t9*wght-2.4
     3         814019635976003d+0*t1*t13*t14*t18*t9*wght-9.3052573634910
     4         02d-1*t1*t11*t12*(9.777777777777779d+0*t30*t5+6.755555555
     5         555556d+1*gammabb*t2*t29)*wght-2.4814019635976003d+0*t1*t
     6         11*t15*t34*wght)+t28*(3.7221029453964005d+0*t1*t11*t12*t1
     7         8*t34*wght-2*t21*t32)-2.8949689575305343d+0*t1*t13*t25*t8
     8         *t9*wght+1.8610514726982003d+0*t1*t13*t20*t27*t28*t9*wght
     9         -9.925607854390401d+0*gammabb*t1*t13*t14*t18*t21*t28*t9*w
     :         ght-2.481401963597601d+0*gammabb*t1*t13*t19*t21*t25*t9*wg
     ;         ht-3.9702431417561607d+1*t1*t3*t38*t39*t8*wght/rhob**8.66
     <         6666666666666d+0+3.143109153890294d+1*gammabb*t1*t22*t23*
     =         t24*t8*wght-5.583154418094601d+0*t1*t11*t12*t21*t27*t40*w
     >         ght-9.925607854390401d+0*gammabb*t1*t18*t19*t22*t23*t33*w
     ?         ght+1.3234143805853868d+1*t1*t19*t21*t22*t23*t24*t3*wght-
     @         4.135669939329334d-1*t1*t11*t19*t21*t26*wght+Cmat3(iq,D3_
     1         RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = t19*(-9.305257363491002d-1*t1*t13*
     1         t20*t34*t9*wght-1.2407009817988002d+0*t1*t13*t14*t21*t9*w
     2         ght+7.444205890792801d+0*t1*t11*t2*t33*wght)+t28*(1.86105
     3         14726982003d+0*t1*t11*t12*t21*t34*wght-2*t21*t35)+4.96280
     4         39271952007d+0*gammabb*t1*t13*t19*t2*t41*t9*wght-4.962803
     5         9271952007d+0*gammabb*t1*t13*t14*t28*t37*t9*wght-9.305257
     6         363491002d-1*t1*t13*t19*t20*t34*t9*wght+3.722102945396400
     7         5d+0*t1*t13*t18*t20*t21*t28*t9*wght+3.7221029453964005d+0
     8         *t1*t13*t14*t19*t21*t9*wght+1.4888411781585603d+1*gammabb
     9         *t1*t38*t39*t41*t8*wght-7.444205890792801d+0*t1*t22*t23*t
     :         33*t8*wght-5.583154418094601d+0*t1*t11*t12*t18*t37*t40*wg
     ;         ht+3.7221029453964005d+0*t1*t11*t18*t2*t28*t36*wght+1.861
     <         0514726982003d+0*t1*t18*t19*t22*t23*t36*wght-9.9256078543
     =         90401d+0*gammabb*t1*t19*t21*t22*t23*t33*wght+Cmat3(iq,D3_
     >         RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = -5.583154418094601d+0*t1*t13*t19*
     1         t2*t42*t9*wght+5.583154418094601d+0*t1*t13*t20*t28*t37*t9
     2         *wght-5.583154418094601d+0*t1*t38*t39*t42*t8*wght-5.58315
     3         4418094601d+0*t1*t11*t12*t21**3*t40*wght+1.11663088361892
     4         02d+1*t1*t11*t2*t21*t28*t36*wght+5.583154418094601d+0*t1*
     5         t19*t21*t22*t23*t36*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
