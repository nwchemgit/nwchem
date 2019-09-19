C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_pkzb99.F 
C> The nwxcm_x_pkzb99 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_pkzb99 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{0.09873184707612337\,\tau_\alpha}
C>    \over{\rho_\alpha^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{0.09873184707612337\,\tau_\beta}
C>    \over{\rho_\beta^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{0.09873184707612337\,\tau_s}
C>    \over{\rho_s^{{{5}\over{3}}}}}-0.001371275653835046\,{\it t_8}
C>    \,\sigma_{ss}-0.45\\\\ 
C>   f &=& -0.9305257363490997\,{\it t_2}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_6}^2-0.002966018451258026\,{\it t_5}\,
C>    \sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_1}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_4}^2-0.002966018451258026\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}
C>    \over{3}}}}}+0.002031519487163032\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\right)+1.0}}\right)-0.9305257363490997
C>    \,{\it t_2}-0.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,{\it t_7}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_9}^2-0.002966018451258026\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pkzb99(param,tol_rho,ipol,nq,wght, 
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
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = 3.656735076893457d-3*gammaaa*t10-1.6455307846020562d-1
     1         *t4*taua
            t12 = 1/t7**2
            fnc(iq) = (-1.8610514726981994d+0*t1*t8-1.8610514726981994d+
     1         0*t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t12*(1.44197530864197
     2         54d-1*t11*t6+7.909382536688069d-3*gammaaa*t10*t6-2.966018
     3         451258026d-3*gammaaa*t11*t4-1.9056529108825535d-4*t2/rhoa
     4         **6.333333333333333d+0-5.417385299101418d-3*gammaaa*t10))
     5         *wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         2*(t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+7.55
     2         2921304913142d-5*gammaaa*t3)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t12*
     1         (1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaaa/
     2         rhoa**4.333333333333333d+0)*wght
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
            t2 = rhob**1.3333333333333333d+0
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 1/rhoa**1.6666666666666669d+0
            t7 = 9.873184707612337d-2*t6*taua-1.3712756538350462d-3*gamm
     1         aaa*t5-4.5d-1
            t8 = 1.2437810945273631d+0*(7.209876543209877d-2*t7**2-2.966
     1         018451258026d-3*gammaaa*t5*t7+2.031519487163032d-3*gammaa
     2         a*t5+3.573099207904788d-5*t3*t4)+1.0d+0
            t9 = 8.04d-1-8.04d-1/t8
            t10 = gammabb**2
            t11 = 1/rhob**5.333333333333333d+0
            t12 = 1/rhob**2.6666666666666666d+0
            t13 = 1/rhob**1.6666666666666669d+0
            t14 = 9.873184707612337d-2*t13*taub-1.3712756538350462d-3*ga
     1         mmabb*t12-4.5d-1
            t15 = 1.2437810945273631d+0*(7.209876543209877d-2*t14**2-2.9
     1         66018451258026d-3*gammabb*t12*t14+2.031519487163032d-3*ga
     2         mmabb*t12+3.573099207904788d-5*t10*t11)+1.0d+0
            t16 = 8.04d-1-8.04d-1/t15
            t17 = rhoa**3.333333333333333d-1
            t18 = 1/rhoa**3.6666666666666664d+0
            t19 = 3.656735076893457d-3*gammaaa*t18-1.6455307846020562d-1
     1         *t5*taua
            t20 = 1/t8**2
            t21 = rhob**3.333333333333333d-1
            t22 = 1/rhob**3.6666666666666664d+0
            t23 = 3.656735076893457d-3*gammabb*t22-1.6455307846020562d-1
     1         *t12*taub
            t24 = 1/t15**2
            fnc(iq) = (-9.305257363490997d-1*t1*t9-9.305257363490997d-1*
     1         t16*t2-9.305257363490997d-1*t2-9.305257363490997d-1*t1)*w
     2         ght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t17*t9-9.3052573634
     1         90997d-1*t1*t20*(1.4419753086419754d-1*t19*t7+7.909382536
     2         688069d-3*gammaaa*t18*t7-2.966018451258026d-3*gammaaa*t19
     3         *t5-1.9056529108825535d-4*t3/rhoa**6.333333333333333d+0-5
     4         .417385299101418d-3*gammaaa*t18)-1.2407009817987995d+0*t1
     5         7)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t2*(1.44197530864197
     1         54d-1*t14*t23-2.966018451258026d-3*gammabb*t12*t23+7.9093
     2         82536688069d-3*gammabb*t14*t22-5.417385299101418d-3*gamma
     3         bb*t22-1.9056529108825535d-4*t10/rhob**6.333333333333333d
     4         +0)*t24-1.2407009817987995d+0*t16*t21-1.2407009817987995d
     5         +0*t21)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t2
     1         0*(t5*(2.031519487163032d-3-3.163753014675228d-3*t7)+7.55
     2         2921304913142d-5*gammaaa*t4)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*(t12*
     1         (2.031519487163032d-3-3.163753014675228d-3*t14)+7.5529213
     2         04913142d-5*gammabb*t11)*t2*t24*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t20*
     1         (1.423688856603853d-2*t6*t7-2.928404801545677d-4*gammaaa/
     2         rhoa**4.333333333333333d+0)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*(1.4236
     1         88856603853d-2*t13*t14-2.928404801545677d-4*gammabb/rhob*
     2         *4.333333333333333d+0)*t2*t24*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**3.6666666666666664d+0
            t11 = 3.656735076893457d-3*gammaaa*t10-1.6455307846020562d-1
     1         *t4*taua
            t12 = 1/t7**2
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t12*(1.44197530864197
     2         54d-1*t11*t6+7.909382536688069d-3*gammaaa*t10*t6-2.966018
     3         451258026d-3*gammaaa*t11*t4-1.9056529108825535d-4*t2/rhoa
     4         **6.333333333333333d+0-5.417385299101418d-3*gammaaa*t10))
     5         *wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         2*(t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+7.55
     2         2921304913142d-5*gammaaa*t3)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t12*
     1         (1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaaa/
     2         rhoa**4.333333333333333d+0)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = gammabb**2
            t3 = 1/rhob**5.333333333333333d+0
            t4 = 1/rhob**2.6666666666666666d+0
            t5 = 1/rhob**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taub-1.3712756538350462d-3*gamm
     1         abb*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammabb*t4*t6+2.031519487163032d-3*gammab
     2         b*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhob**3.333333333333333d-1
            t10 = 1/rhob**3.6666666666666664d+0
            t11 = 3.656735076893457d-3*gammabb*t10-1.6455307846020562d-1
     1         *t4*taub
            t12 = 1/t7**2
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RB) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t12*(1.44197530864197
     2         54d-1*t11*t6+7.909382536688069d-3*gammabb*t10*t6-2.966018
     3         451258026d-3*gammabb*t11*t4-1.9056529108825535d-4*t2/rhob
     4         **6.333333333333333d+0-5.417385299101418d-3*gammabb*t10))
     5         *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t1*t1
     1         2*(t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+7.55
     2         2921304913142d-5*gammabb*t3)*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t1*t12*
     1         (1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammabb/
     2         rhob**4.333333333333333d+0)*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pkzb99 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{0.09873184707612337\,\tau_\alpha}
C>    \over{\rho_\alpha^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{0.09873184707612337\,\tau_\beta}
C>    \over{\rho_\beta^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{0.09873184707612337\,\tau_s}
C>    \over{\rho_s^{{{5}\over{3}}}}}-0.001371275653835046\,{\it t_8}
C>    \,\sigma_{ss}-0.45\\\\ 
C>   f &=& -0.9305257363490997\,{\it t_2}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_6}^2-0.002966018451258026\,{\it t_5}\,
C>    \sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_1}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_4}^2-0.002966018451258026\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}
C>    \over{3}}}}}+0.002031519487163032\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\right)+1.0}}\right)-0.9305257363490997
C>    \,{\it t_2}-0.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,{\it t_7}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_9}^2-0.002966018451258026\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pkzb99_d2(param,tol_rho,ipol,nq,wght, 
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
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**6.333333333333333d+0
            t11 = 1/rhoa**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammaaa*t11-1.6455307846020562d-1
     1         *t4*taua
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         aaa*t11*t6-2.966018451258026d-3*gammaaa*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammaaa*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammaaa*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t18 = 1/t7**3
            t19 = 1/rhoa**4.666666666666667d+0
            t20 = 4.38808209227215d-1*t11*taua-1.3408028615276008d-2*gam
     1         maaa*t19
            t21 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            fnc(iq) = (-1.8610514726981994d+0*t1*t8-1.8610514726981994d+
     1         0*t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Amat2(iq,D2_RA_RA) = (t14*(-2.481401963597599d+0*t13*t9-9.30
     1         5257363490997d-1*t1*(1.4419753086419754d-1*t20*t6-2.90010
     2         69301189586d-2*gammaaa*t19*t6-2.966018451258026d-3*gammaa
     3         a*t20*t4+1.2069135102256172d-3*t2/rhoa**7.333333333333333
     4         d+0+1.98637460967052d-2*gammaaa*t19+1.4419753086419754d-1
     5         *t12**2+1.5818765073376137d-2*gammaaa*t11*t12))+(-4.13566
     6         99393293317d-1*t8-4.1356699393293317d-1)/rhoa**6.66666666
     7         6666666d-1+2.3147406376843277d+0*t1*t13**2*t18)*wght+Amat
     8         2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t14*(-1.2407009817987995d+0*t21*t9-9.
     1         305257363490997d-1*t1*(8.436674705800608d-3*t11*t6-3.1637
     2         53014675228d-3*t12*t4-5.417385299101418d-3*t11-4.02822469
     3         59536755d-4*gammaaa*t10))+2.3147406376843277d+0*t1*t13*t1
     4         8*t21)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 2.3147406376843277d+0*t1*t16*t18*t21*
     1         wght-7.431884849105774d-5*t14*wght/rhoa**4+Cmat2(iq,D2_GA
     2         A_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (t14*(-1.2407009817987995d+0*t17*t9-9.3
     1         05257363490997d-1*t1*(-2.372814761006422d-2*t4*t6+1.42368
     2         8856603853d-2*t12*t5+1.2689754140031267d-3*gammaaa*t3))+2
     3         .3147406376843277d+0*t1*t13*t17*t18)*wght+Mmat2(iq,D2_RA_
     4         TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 2.3147406376843277d+0*t1*t16*t17*t18*w
     1         ght+2.906619769905633d-4*t14*wght/rhoa**3+Mmat2(iq,D2_GAA
     2         _TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = 2.3147406376843277d+0*t1*t17**2*t18*wgh
     1         t-1.3079788964575353d-3*t14*wght/rhoa**2+Mmat2(iq,D2_TA_T
     2         A)
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
            t2 = rhob**1.3333333333333333d+0
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 1/rhoa**1.6666666666666669d+0
            t7 = 9.873184707612337d-2*t6*taua-1.3712756538350462d-3*gamm
     1         aaa*t5-4.5d-1
            t8 = 1.2437810945273631d+0*(7.209876543209877d-2*t7**2-2.966
     1         018451258026d-3*gammaaa*t5*t7+2.031519487163032d-3*gammaa
     2         a*t5+3.573099207904788d-5*t3*t4)+1.0d+0
            t9 = 8.04d-1-8.04d-1/t8
            t10 = gammabb**2
            t11 = 1/rhob**5.333333333333333d+0
            t12 = 1/rhob**2.6666666666666666d+0
            t13 = 1/rhob**1.6666666666666669d+0
            t14 = 9.873184707612337d-2*t13*taub-1.3712756538350462d-3*ga
     1         mmabb*t12-4.5d-1
            t15 = 1.2437810945273631d+0*(7.209876543209877d-2*t14**2-2.9
     1         66018451258026d-3*gammabb*t12*t14+2.031519487163032d-3*ga
     2         mmabb*t12+3.573099207904788d-5*t10*t11)+1.0d+0
            t16 = 8.04d-1-8.04d-1/t15
            t17 = rhoa**3.333333333333333d-1
            t18 = 1/rhoa**6.333333333333333d+0
            t19 = 1/rhoa**3.6666666666666664d+0
            t20 = 3.656735076893457d-3*gammaaa*t19-1.6455307846020562d-1
     1         *t5*taua
            t21 = 1.4419753086419754d-1*t20*t7+7.909382536688069d-3*gamm
     1         aaa*t19*t7-2.966018451258026d-3*gammaaa*t20*t5-1.90565291
     2         08825535d-4*t18*t3-5.417385299101418d-3*gammaaa*t19
            t22 = 1/t8**2
            t23 = rhob**3.333333333333333d-1
            t24 = 1/rhob**6.333333333333333d+0
            t25 = 1/rhob**3.6666666666666664d+0
            t26 = 3.656735076893457d-3*gammabb*t25-1.6455307846020562d-1
     1         *t12*taub
            t27 = 1.4419753086419754d-1*t14*t26-2.966018451258026d-3*gam
     1         mabb*t12*t26+7.909382536688069d-3*gammabb*t14*t25-5.41738
     2         5299101418d-3*gammabb*t25-1.9056529108825535d-4*t10*t24
            t28 = 1/t15**2
            t29 = 7.552921304913142d-5*gammaaa*t4
            t30 = t5*(2.031519487163032d-3-3.163753014675228d-3*t7)+t29
            t31 = 7.552921304913142d-5*gammabb*t11
            t32 = t31+t12*(2.031519487163032d-3-3.163753014675228d-3*t14
     1         )
            t33 = 1.423688856603853d-2*t6*t7-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t34 = 1.423688856603853d-2*t13*t14-2.928404801545677d-4*gamm
     1         abb/rhob**4.333333333333333d+0
            t35 = 1/t8**3
            t36 = 1/rhoa**4.666666666666667d+0
            t37 = 4.38808209227215d-1*t19*taua-1.3408028615276008d-2*gam
     1         maaa*t36
            t38 = 1/t15**3
            t39 = 1/rhob**4.666666666666667d+0
            t40 = 4.38808209227215d-1*t25*taub-1.3408028615276008d-2*gam
     1         mabb*t39
            t41 = -3.163753014675228d-3*t5*t7+2.031519487163032d-3*t5+t2
     1         9
            t42 = t31-3.163753014675228d-3*t12*t14+2.031519487163032d-3*
     1         t12
            fnc(iq) = (-9.305257363490997d-1*t1*t9-9.305257363490997d-1*
     1         t16*t2-9.305257363490997d-1*t2-9.305257363490997d-1*t1)*w
     2         ght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t17*t9-9.3052573634
     1         90997d-1*t1*t21*t22-1.2407009817987995d+0*t17)*wght+Amat(
     2         iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t2*t27*t28-1.2407009
     1         817987995d+0*t16*t23-1.2407009817987995d+0*t23)*wght+Amat
     2         (iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t2
     1         2*t30*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t2*t2
     1         8*t32*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t22*
     1         t33*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t2*t28*
     1         t34*wght
            Amat2(iq,D2_RA_RA) = ((-4.1356699393293317d-1*t9-4.135669939
     1         3293317d-1)/rhoa**6.666666666666666d-1+t22*(-9.3052573634
     2         90997d-1*t1*(1.4419753086419754d-1*t37*t7-2.9001069301189
     3         586d-2*gammaaa*t36*t7-2.966018451258026d-3*gammaaa*t37*t5
     4         +1.98637460967052d-2*gammaaa*t36+1.2069135102256172d-3*t3
     5         /rhoa**7.333333333333333d+0+1.4419753086419754d-1*t20**2+
     6         1.5818765073376137d-2*gammaaa*t19*t20)-2.481401963597599d
     7         +0*t17*t21)+2.3147406376843277d+0*t1*t21**2*t35)*wght+Ama
     8         t2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t28*(-9.305257363490997d-1*t2*(1.44197
     1         53086419754d-1*t14*t40-2.966018451258026d-3*gammabb*t12*t
     2         40-2.9001069301189586d-2*gammabb*t14*t39+1.98637460967052
     3         d-2*gammabb*t39+1.4419753086419754d-1*t26**2+1.5818765073
     4         376137d-2*gammabb*t25*t26+1.2069135102256172d-3*t10/rhob*
     5         *7.333333333333333d+0)-2.481401963597599d+0*t23*t27)+2.31
     6         47406376843277d+0*t2*t27**2*t38+(-4.1356699393293317d-1*t
     7         16-4.1356699393293317d-1)/rhob**6.666666666666666d-1)*wgh
     8         t+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t22*(-9.305257363490997d-1*t1*(8.4366
     1         74705800608d-3*t19*t7-3.163753014675228d-3*t20*t5-5.41738
     2         5299101418d-3*t19-4.0282246959536755d-4*gammaaa*t18)-1.24
     3         07009817987995d+0*t17*t41)+2.3147406376843277d+0*t1*t21*t
     4         35*t41)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t28*(-1.2407009817987995d+0*t23*t42-9
     1         .305257363490997d-1*t2*(-3.163753014675228d-3*t12*t26+8.4
     2         36674705800608d-3*t14*t25-5.417385299101418d-3*t25-4.0282
     3         246959536755d-4*gammabb*t24))+2.3147406376843277d+0*t2*t2
     4         7*t38*t42)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 2.3147406376843277d+0*t1*t30*t35*t41*
     1         wght-7.431884849105774d-5*t22*wght/rhoa**4+Cmat2(iq,D2_GA
     2         A_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 2.3147406376843277d+0*t2*t32*t38*t42*
     1         wght-7.431884849105774d-5*t28*wght/rhob**4+Cmat2(iq,D2_GB
     2         B_GBB)
            Mmat2(iq,D2_RA_TA) = (t22*(-9.305257363490997d-1*t1*(-2.3728
     1         14761006422d-2*t5*t7+1.423688856603853d-2*t20*t6+1.268975
     2         4140031267d-3*gammaaa*t4)-1.2407009817987995d+0*t17*t33)+
     3         2.3147406376843277d+0*t1*t21*t33*t35)*wght+Mmat2(iq,D2_RA
     4         _TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (2.3147406376843277d+0*t2*t27*t34*t38+t
     1         28*(-1.2407009817987995d+0*t23*t34-9.305257363490997d-1*t
     2         2*(1.423688856603853d-2*t13*t26-2.372814761006422d-2*t12*
     3         t14+1.2689754140031267d-3*gammabb*t11)))*wght+Mmat2(iq,D2
     4         _RB_TB)
            Mmat2(iq,D2_GAA_TA) = 2.3147406376843277d+0*t1*t30*t33*t35*w
     1         ght+2.906619769905633d-4*t22*wght/rhoa**3+Mmat2(iq,D2_GAA
     2         _TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 2.3147406376843277d+0*t2*t32*t34*t38*w
     1         ght+2.906619769905633d-4*t28*wght/rhob**3+Mmat2(iq,D2_GBB
     2         _TB)
            Mmat2(iq,D2_TA_TA) = 2.3147406376843277d+0*t1*t33**2*t35*wgh
     1         t-1.3079788964575353d-3*t22*wght/rhoa**2+Mmat2(iq,D2_TA_T
     2         A)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = 2.3147406376843277d+0*t2*t34**2*t38*wgh
     1         t-1.3079788964575353d-3*t28*wght/rhob**2+Mmat2(iq,D2_TB_T
     2         B)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**6.333333333333333d+0
            t11 = 1/rhoa**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammaaa*t11-1.6455307846020562d-1
     1         *t4*taua
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         aaa*t11*t6-2.966018451258026d-3*gammaaa*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammaaa*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammaaa*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t18 = 1/t7**3
            t19 = 1/rhoa**4.666666666666667d+0
            t20 = 4.38808209227215d-1*t11*taua-1.3408028615276008d-2*gam
     1         maaa*t19
            t21 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = (t14*(-2.481401963597599d+0*t13*t9-9.30
     1         5257363490997d-1*t1*(1.4419753086419754d-1*t20*t6-2.90010
     2         69301189586d-2*gammaaa*t19*t6-2.966018451258026d-3*gammaa
     3         a*t20*t4+1.2069135102256172d-3*t2/rhoa**7.333333333333333
     4         d+0+1.98637460967052d-2*gammaaa*t19+1.4419753086419754d-1
     5         *t12**2+1.5818765073376137d-2*gammaaa*t11*t12))+(-4.13566
     6         99393293317d-1*t8-4.1356699393293317d-1)/rhoa**6.66666666
     7         6666666d-1+2.3147406376843277d+0*t1*t13**2*t18)*wght+Amat
     8         2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (t14*(-1.2407009817987995d+0*t21*t9-9.
     1         305257363490997d-1*t1*(8.436674705800608d-3*t11*t6-3.1637
     2         53014675228d-3*t12*t4-5.417385299101418d-3*t11-4.02822469
     3         59536755d-4*gammaaa*t10))+2.3147406376843277d+0*t1*t13*t1
     4         8*t21)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 2.3147406376843277d+0*t1*t16*t18*t21*
     1         wght-7.431884849105774d-5*t14*wght/rhoa**4+Cmat2(iq,D2_GA
     2         A_GAA)
            Mmat2(iq,D2_RA_TA) = (t14*(-1.2407009817987995d+0*t17*t9-9.3
     1         05257363490997d-1*t1*(-2.372814761006422d-2*t4*t6+1.42368
     2         8856603853d-2*t12*t5+1.2689754140031267d-3*gammaaa*t3))+2
     3         .3147406376843277d+0*t1*t13*t17*t18)*wght+Mmat2(iq,D2_RA_
     4         TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 2.3147406376843277d+0*t1*t16*t17*t18*w
     1         ght+2.906619769905633d-4*t14*wght/rhoa**3+Mmat2(iq,D2_GAA
     2         _TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = 2.3147406376843277d+0*t1*t17**2*t18*wgh
     1         t-1.3079788964575353d-3*t14*wght/rhoa**2+Mmat2(iq,D2_TA_T
     2         A)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = gammabb**2
            t3 = 1/rhob**5.333333333333333d+0
            t4 = 1/rhob**2.6666666666666666d+0
            t5 = 1/rhob**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taub-1.3712756538350462d-3*gamm
     1         abb*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammabb*t4*t6+2.031519487163032d-3*gammab
     2         b*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhob**3.333333333333333d-1
            t10 = 1/rhob**6.333333333333333d+0
            t11 = 1/rhob**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammabb*t11-1.6455307846020562d-1
     1         *t4*taub
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         abb*t11*t6-2.966018451258026d-3*gammabb*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammabb*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammabb*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammab
     1         b/rhob**4.333333333333333d+0
            t18 = 1/t7**3
            t19 = 1/rhob**4.666666666666667d+0
            t20 = 4.38808209227215d-1*t11*taub-1.3408028615276008d-2*gam
     1         mabb*t19
            t21 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RB) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = (t14*(-2.481401963597599d+0*t13*t9-9.30
     1         5257363490997d-1*t1*(1.4419753086419754d-1*t20*t6-2.90010
     2         69301189586d-2*gammabb*t19*t6-2.966018451258026d-3*gammab
     3         b*t20*t4+1.2069135102256172d-3*t2/rhob**7.333333333333333
     4         d+0+1.98637460967052d-2*gammabb*t19+1.4419753086419754d-1
     5         *t12**2+1.5818765073376137d-2*gammabb*t11*t12))+(-4.13566
     6         99393293317d-1*t8-4.1356699393293317d-1)/rhob**6.66666666
     7         6666666d-1+2.3147406376843277d+0*t1*t13**2*t18)*wght+Amat
     8         2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (t14*(-1.2407009817987995d+0*t21*t9-9.
     1         305257363490997d-1*t1*(8.436674705800608d-3*t11*t6-3.1637
     2         53014675228d-3*t12*t4-5.417385299101418d-3*t11-4.02822469
     3         59536755d-4*gammabb*t10))+2.3147406376843277d+0*t1*t13*t1
     4         8*t21)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 2.3147406376843277d+0*t1*t16*t18*t21*
     1         wght-7.431884849105774d-5*t14*wght/rhob**4+Cmat2(iq,D2_GB
     2         B_GBB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (t14*(-1.2407009817987995d+0*t17*t9-9.3
     1         05257363490997d-1*t1*(-2.372814761006422d-2*t4*t6+1.42368
     2         8856603853d-2*t12*t5+1.2689754140031267d-3*gammabb*t3))+2
     3         .3147406376843277d+0*t1*t13*t17*t18)*wght+Mmat2(iq,D2_RB_
     4         TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 2.3147406376843277d+0*t1*t16*t17*t18*w
     1         ght+2.906619769905633d-4*t14*wght/rhob**3+Mmat2(iq,D2_GBB
     2         _TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = 2.3147406376843277d+0*t1*t17**2*t18*wgh
     1         t-1.3079788964575353d-3*t14*wght/rhob**2+Mmat2(iq,D2_TB_T
     2         B)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pkzb99 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& \rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   {\it t_2} &=& \rho_\beta^{{{4}\over{3}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {{0.09873184707612337\,\tau_\alpha}
C>    \over{\rho_\alpha^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{0.09873184707612337\,\tau_\beta}
C>    \over{\rho_\beta^{{{5}\over{3}}}}}-0.001371275653835046\,{
C>    \it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{0.09873184707612337\,\tau_s}
C>    \over{\rho_s^{{{5}\over{3}}}}}-0.001371275653835046\,{\it t_8}
C>    \,\sigma_{ss}-0.45\\\\ 
C>   f &=& -0.9305257363490997\,{\it t_2}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_6}^2-0.002966018451258026\,{\it t_5}\,
C>    \sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_1}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_4}^2-0.002966018451258026\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{
C>    -5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}
C>    \over{3}}}}}+0.002031519487163032\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}\right)+1.0}}\right)-0.9305257363490997
C>    \,{\it t_2}-0.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363490997\,{\it t_7}\,\left(0.804
C>    -{{0.804}\over{1.243781094527363\,\left(0.07209876543209877
C>    \,{\it t_9}^2-0.002966018451258026\,{\it t_8}\,\sigma_{ss}
C>    \,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,
C>    \sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}
C>    +0.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)
C>    +1.0}}\right)-0.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82
C>    , 2544 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pkzb99_d3(param,tol_rho,ipol,nq,wght, 
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
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**6.333333333333333d+0
            t11 = 1/rhoa**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammaaa*t11-1.6455307846020562d-1
     1         *t4*taua
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         aaa*t11*t6-2.966018451258026d-3*gammaaa*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammaaa*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammaaa*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t18 = t13**2
            t19 = 1/t7**3
            t20 = 1/rhoa**7.333333333333333d+0
            t21 = 1/rhoa**4.666666666666667d+0
            t22 = 4.38808209227215d-1*t11*taua-1.3408028615276008d-2*gam
     1         maaa*t21
            t23 = 1.4419753086419754d-1*t22*t6-2.9001069301189586d-2*gam
     1         maaa*t21*t6-2.966018451258026d-3*gammaaa*t22*t4+1.9863746
     2         0967052d-2*gammaaa*t21+1.2069135102256172d-3*t2*t20+1.441
     3         9753086419754d-1*t12**2+1.5818765073376137d-2*gammaaa*t11
     4         *t12
            t24 = -2.481401963597599d+0*t13*t9-9.305257363490997d-1*t1*t
     1         23
            t25 = 1/rhoa**6.666666666666666d-1
            t26 = -4.1356699393293317d-1*t8-4.1356699393293317d-1
            t27 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            t28 = 8.436674705800608d-3*t11*t6-3.163753014675228d-3*t12*t
     1         4-5.417385299101418d-3*t11-4.0282246959536755d-4*gammaaa*
     2         t10
            t29 = -1.2407009817987995d+0*t27*t9-9.305257363490997d-1*t1*
     1         t28
            t30 = 1/rhoa**4
            t31 = -2.372814761006422d-2*t4*t6+1.423688856603853d-2*t12*t
     1         5+1.2689754140031267d-3*gammaaa*t3
            t32 = 1/rhoa**3
            t33 = 1/rhoa**2
            t34 = 1/t7**4
            t35 = 1/rhoa**5.666666666666667d+0
            t36 = 6.257080020462137d-2*gammaaa*t35-1.6089634338331216d+0
     1         *t21*taua
            t37 = t27**2
            fnc(iq) = (-1.8610514726981994d+0*t1*t8-1.8610514726981994d+
     1         0*t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Amat2(iq,D2_RA_RA) = (t25*t26+t14*t24+2.3147406376843277d+0*
     1         t1*t18*t19)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (t14*t29+2.3147406376843277d+0*t1*t13*
     1         t19*t27)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.431884849105774d-5*t14*t30*wght+2.
     1         3147406376843277d+0*t1*t16*t19*t27*wght+Cmat2(iq,D2_GAA_G
     2         AA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Mmat2(iq,D2_RA_TA) = (t14*(-1.2407009817987995d+0*t17*t9-9.3
     1         05257363490997d-1*t1*t31)+2.3147406376843277d+0*t1*t13*t1
     2         7*t19)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 2.906619769905633d-4*t14*t32*wght+2.31
     1         47406376843277d+0*t1*t16*t17*t19*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_TA_TA) = -1.3079788964575353d-3*t14*t33*wght+2.3
     1         147406376843277d+0*t1*t17**2*t19*wght+Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Amat3(iq,D3_RA_RA_RA) = (t14*(-3.7221029453963983d+0*t23*t9-
     1         9.305257363490997d-1*t1*(1.4419753086419754d-1*t36*t6+1.3
     2         53383234055514d-1*gammaaa*t35*t6-2.966018451258026d-3*gam
     3         maaa*t36*t4-9.269748178462428d-2*gammaaa*t35+4.3259259259
     4         259264d-1*t12*t22+2.3728147610064204d-2*gammaaa*t11*t22-8
     5         .700320790356875d-2*gammaaa*t12*t21-8.85069907498786d-3*t
     6         2/rhoa**8.333333333333334d+0)-8.271339878658663d-1*t13*t2
     7         5)+t19*(3.08632085024577d+0*t18*t9-2.4875621890547261d+0*
     8         t13*t24+4.6294812753686554d+0*t1*t13*t23)-6.6666666666666
     9         66d-1*t26*t5-8.637091931657939d+0*t1*t13**3*t34-4.1356699
     :         393293317d-1*t13*t14*t25)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t14*(-2.481401963597599d+0*t28*t9-
     1         9.305257363490997d-1*t1*(-3.093447392126889d-2*t21*t6-3.1
     2         63753014675228d-3*t22*t4+1.98637460967052d-2*t21+2.551208
     3         9741039945d-3*gammaaa*t20+1.6873349411601213d-2*t11*t12))
     4         -8.637091931657939d+0*t1*t18*t27*t34+t19*(4.6294812753686
     5         554d+0*t1*t13*t28-2.4875621890547261d+0*t24*t27)-4.135669
     6         9393293317d-1*t14*t25*t27)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-8.637091931657939d+0*t1*t13*t34*
     1         t37+1.8487275744044215d-4*t13*t19*t30+t19*(2.314740637684
     2         3277d+0*t1*t27*t28-2.4875621890547261d+0*t27*t29)+2.97275
     3         393964231d-4*t14/rhoa**5)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t19*t30*(3.697455148808843d-4*t27
     1         *wght+1.8487275744044215d-4*t16*wght)-8.637091931657939d+
     2         0*t1*t16*t34*t37*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Mmat3(iq,D3_RA_RA_TA) = (t14*(-2.481401963597599d+0*t31*t9-9
     1         .305257363490997d-1*t1*(6.327506029350458d-2*t11*t6+1.423
     2         688856603853d-2*t22*t5-4.745629522012844d-2*t12*t4-6.7678
     3         68874683343d-3*gammaaa*t10))-8.637091931657939d+0*t1*t17*
     4         t18*t34+t19*(4.6294812753686554d+0*t1*t13*t31-2.487562189
     5         0547261d+0*t17*t24)-4.1356699393293317d-1*t14*t17*t25)*wg
     6         ht+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_GAA_TA) = (-8.637091931657939d+0*t1*t13*t17*t
     1         27*t34-7.2303974375762d-4*t13*t19*t32+t19*(2.314740637684
     2         3277d+0*t1*t27*t31-2.4875621890547261d+0*t17*t29)-8.71985
     3         9309716897d-4*t14*t30)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_RA_GAB_TA) = Mmat3(iq,D3_RA_GAB_TA)
            Mmat3(iq,D3_RA_GAB_TB) = Mmat3(iq,D3_RA_GAB_TB)
            Mmat3(iq,D3_RA_GBB_TA) = Mmat3(iq,D3_RA_GBB_TA)
            Mmat3(iq,D3_RA_GBB_TB) = Mmat3(iq,D3_RA_GBB_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = t19*t32*(-7.2303974375762d-4*t27*w
     1         ght-7.2303974375762d-4*t16*wght)-8.637091931657939d+0*t1*
     2         t16*t17*t27*t34*wght+1.8487275744044215d-4*t17*t19*t30*wg
     3         ht+Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_GAA_GAB_TA) = Mmat3(iq,D3_GAA_GAB_TA)
            Mmat3(iq,D3_GAA_GAB_TB) = Mmat3(iq,D3_GAA_GAB_TB)
            Mmat3(iq,D3_GAA_GBB_TA) = Mmat3(iq,D3_GAA_GBB_TA)
            Mmat3(iq,D3_GAB_GAB_TA) = Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_TA_TA_TA) = -8.637091931657939d+0*t1*t17**3*t34*
     1         wght+9.761036540727874d-3*t17*t19*t33*wght+Mmat3(iq,D3_TA
     2         _TA_TA)
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
            t2 = rhob**1.3333333333333333d+0
            t3 = gammaaa**2
            t4 = 1/rhoa**5.333333333333333d+0
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 1/rhoa**1.6666666666666669d+0
            t7 = 9.873184707612337d-2*t6*taua-1.3712756538350462d-3*gamm
     1         aaa*t5-4.5d-1
            t8 = 1.2437810945273631d+0*(7.209876543209877d-2*t7**2-2.966
     1         018451258026d-3*gammaaa*t5*t7+2.031519487163032d-3*gammaa
     2         a*t5+3.573099207904788d-5*t3*t4)+1.0d+0
            t9 = 8.04d-1-8.04d-1/t8
            t10 = gammabb**2
            t11 = 1/rhob**5.333333333333333d+0
            t12 = 1/rhob**2.6666666666666666d+0
            t13 = 1/rhob**1.6666666666666669d+0
            t14 = 9.873184707612337d-2*t13*taub-1.3712756538350462d-3*ga
     1         mmabb*t12-4.5d-1
            t15 = 1.2437810945273631d+0*(7.209876543209877d-2*t14**2-2.9
     1         66018451258026d-3*gammabb*t12*t14+2.031519487163032d-3*ga
     2         mmabb*t12+3.573099207904788d-5*t10*t11)+1.0d+0
            t16 = 8.04d-1-8.04d-1/t15
            t17 = rhoa**3.333333333333333d-1
            t18 = 1/rhoa**6.333333333333333d+0
            t19 = 1/rhoa**3.6666666666666664d+0
            t20 = 3.656735076893457d-3*gammaaa*t19-1.6455307846020562d-1
     1         *t5*taua
            t21 = 1.4419753086419754d-1*t20*t7+7.909382536688069d-3*gamm
     1         aaa*t19*t7-2.966018451258026d-3*gammaaa*t20*t5-1.90565291
     2         08825535d-4*t18*t3-5.417385299101418d-3*gammaaa*t19
            t22 = 1/t8**2
            t23 = rhob**3.333333333333333d-1
            t24 = 1/rhob**6.333333333333333d+0
            t25 = 1/rhob**3.6666666666666664d+0
            t26 = 3.656735076893457d-3*gammabb*t25-1.6455307846020562d-1
     1         *t12*taub
            t27 = 1.4419753086419754d-1*t14*t26-2.966018451258026d-3*gam
     1         mabb*t12*t26+7.909382536688069d-3*gammabb*t14*t25-5.41738
     2         5299101418d-3*gammabb*t25-1.9056529108825535d-4*t10*t24
            t28 = 1/t15**2
            t29 = 7.552921304913142d-5*gammaaa*t4
            t30 = t5*(2.031519487163032d-3-3.163753014675228d-3*t7)+t29
            t31 = 7.552921304913142d-5*gammabb*t11
            t32 = t31+t12*(2.031519487163032d-3-3.163753014675228d-3*t14
     1         )
            t33 = 1.423688856603853d-2*t6*t7-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t34 = 1.423688856603853d-2*t13*t14-2.928404801545677d-4*gamm
     1         abb/rhob**4.333333333333333d+0
            t35 = t21**2
            t36 = 1/t8**3
            t37 = 1/rhoa**7.333333333333333d+0
            t38 = 1/rhoa**4.666666666666667d+0
            t39 = 4.38808209227215d-1*t19*taua-1.3408028615276008d-2*gam
     1         maaa*t38
            t40 = 1.4419753086419754d-1*t39*t7-2.9001069301189586d-2*gam
     1         maaa*t38*t7-2.966018451258026d-3*gammaaa*t39*t5+1.9863746
     2         0967052d-2*gammaaa*t38+1.2069135102256172d-3*t3*t37+1.441
     3         9753086419754d-1*t20**2+1.5818765073376137d-2*gammaaa*t19
     4         *t20
            t41 = -9.305257363490997d-1*t1*t40-2.481401963597599d+0*t17*
     1         t21
            t42 = 1/rhoa**6.666666666666666d-1
            t43 = -4.1356699393293317d-1*t9-4.1356699393293317d-1
            t44 = t27**2
            t45 = 1/t15**3
            t46 = 1/rhob**7.333333333333333d+0
            t47 = 1/rhob**4.666666666666667d+0
            t48 = 4.38808209227215d-1*t25*taub-1.3408028615276008d-2*gam
     1         mabb*t47
            t49 = 1.4419753086419754d-1*t14*t48-2.966018451258026d-3*gam
     1         mabb*t12*t48-2.9001069301189586d-2*gammabb*t14*t47+1.9863
     2         7460967052d-2*gammabb*t47+1.2069135102256172d-3*t10*t46+1
     3         .4419753086419754d-1*t26**2+1.5818765073376137d-2*gammabb
     4         *t25*t26
            t50 = -9.305257363490997d-1*t2*t49-2.481401963597599d+0*t23*
     1         t27
            t51 = 1/rhob**6.666666666666666d-1
            t52 = -4.1356699393293317d-1*t16-4.1356699393293317d-1
            t53 = -3.163753014675228d-3*t5*t7+2.031519487163032d-3*t5+t2
     1         9
            t54 = 8.436674705800608d-3*t19*t7-3.163753014675228d-3*t20*t
     1         5-5.417385299101418d-3*t19-4.0282246959536755d-4*gammaaa*
     2         t18
            t55 = -9.305257363490997d-1*t1*t54-1.2407009817987995d+0*t17
     1         *t53
            t56 = t31-3.163753014675228d-3*t12*t14+2.031519487163032d-3*
     1         t12
            t57 = -3.163753014675228d-3*t12*t26+8.436674705800608d-3*t14
     1         *t25-5.417385299101418d-3*t25-4.0282246959536755d-4*gamma
     2         bb*t24
            t58 = -9.305257363490997d-1*t2*t57-1.2407009817987995d+0*t23
     1         *t56
            t59 = 1/rhoa**4
            t60 = 1/rhob**4
            t61 = -2.372814761006422d-2*t5*t7+1.423688856603853d-2*t20*t
     1         6+1.2689754140031267d-3*gammaaa*t4
            t62 = 1.423688856603853d-2*t13*t26-2.372814761006422d-2*t12*
     1         t14+1.2689754140031267d-3*gammabb*t11
            t63 = 1/rhoa**3
            t64 = 1/rhob**3
            t65 = t33**2
            t66 = 1/rhoa**2
            t67 = t34**2
            t68 = 1/rhob**2
            t69 = 1/t8**4
            t70 = 1/rhoa**5.666666666666667d+0
            t71 = 6.257080020462137d-2*gammaaa*t70-1.6089634338331216d+0
     1         *t38*taua
            t72 = 1/t15**4
            t73 = 1/rhob**5.666666666666667d+0
            t74 = 6.257080020462137d-2*gammabb*t73-1.6089634338331216d+0
     1         *t47*taub
            t75 = t53**2
            t76 = t56**2
            fnc(iq) = (-9.305257363490997d-1*t1*t9-9.305257363490997d-1*
     1         t16*t2-9.305257363490997d-1*t2-9.305257363490997d-1*t1)*w
     2         ght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t17*t9-9.3052573634
     1         90997d-1*t1*t21*t22-1.2407009817987995d+0*t17)*wght+Amat(
     2         iq,D1_RA)
            Amat(iq,D1_RB) = (-9.305257363490997d-1*t2*t27*t28-1.2407009
     1         817987995d+0*t16*t23-1.2407009817987995d+0*t23)*wght+Amat
     2         (iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t2
     1         2*t30*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t2*t2
     1         8*t32*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t22*
     1         t33*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t2*t28*
     1         t34*wght
            Amat2(iq,D2_RA_RA) = (t42*t43+t22*t41+2.3147406376843277d+0*
     1         t1*t35*t36)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (t51*t52+t28*t50+2.3147406376843277d+0*
     1         t2*t44*t45)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (t22*t55+2.3147406376843277d+0*t1*t21*
     1         t36*t53)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (t28*t58+2.3147406376843277d+0*t2*t27*
     1         t45*t56)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = -7.431884849105774d-5*t22*t59*wght+2.
     1         3147406376843277d+0*t1*t30*t36*t53*wght+Cmat2(iq,D2_GAA_G
     2         AA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.431884849105774d-5*t28*t60*wght+2.
     1         3147406376843277d+0*t2*t32*t45*t56*wght+Cmat2(iq,D2_GBB_G
     2         BB)
            Mmat2(iq,D2_RA_TA) = (t22*(-9.305257363490997d-1*t1*t61-1.24
     1         07009817987995d+0*t17*t33)+2.3147406376843277d+0*t1*t21*t
     2         33*t36)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (t28*(-9.305257363490997d-1*t2*t62-1.24
     1         07009817987995d+0*t23*t34)+2.3147406376843277d+0*t2*t27*t
     2         34*t45)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GAA_TA) = 2.906619769905633d-4*t22*t63*wght+2.31
     1         47406376843277d+0*t1*t30*t33*t36*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_GAB_TA) = Mmat2(iq,D2_GAB_TA)
            Mmat2(iq,D2_GAB_TB) = Mmat2(iq,D2_GAB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 2.906619769905633d-4*t28*t64*wght+2.31
     1         47406376843277d+0*t2*t32*t34*t45*wght+Mmat2(iq,D2_GBB_TB)
            Mmat2(iq,D2_TA_TA) = -1.3079788964575353d-3*t22*t66*wght+2.3
     1         147406376843277d+0*t1*t36*t65*wght+Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = -1.3079788964575353d-3*t28*t68*wght+2.3
     1         147406376843277d+0*t2*t45*t67*wght+Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = (t22*(-9.305257363490997d-1*t1*(1.44
     1         19753086419754d-1*t7*t71-2.966018451258026d-3*gammaaa*t5*
     2         t71+1.353383234055514d-1*gammaaa*t7*t70-9.269748178462428
     3         d-2*gammaaa*t70+4.3259259259259264d-1*t20*t39+2.372814761
     4         0064204d-2*gammaaa*t19*t39-8.700320790356875d-2*gammaaa*t
     5         20*t38-8.85069907498786d-3*t3/rhoa**8.333333333333334d+0)
     6         -8.271339878658663d-1*t21*t42-3.7221029453963983d+0*t17*t
     7         40)-8.637091931657939d+0*t1*t21**3*t69-6.666666666666666d
     8         -1*t43*t6-4.1356699393293317d-1*t21*t22*t42+t36*(-2.48756
     9         21890547261d+0*t21*t41+4.6294812753686554d+0*t1*t21*t40+3
     :         .08632085024577d+0*t17*t35))*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (t28*(-9.305257363490997d-1*t2*(1.44
     1         19753086419754d-1*t14*t74-2.966018451258026d-3*gammabb*t1
     2         2*t74+1.353383234055514d-1*gammabb*t14*t73-9.269748178462
     3         428d-2*gammabb*t73+4.3259259259259264d-1*t26*t48+2.372814
     4         7610064204d-2*gammabb*t25*t48-8.700320790356875d-2*gammab
     5         b*t26*t47-8.85069907498786d-3*t10/rhob**8.333333333333334
     6         d+0)-8.271339878658663d-1*t27*t51-3.7221029453963983d+0*t
     7         23*t49)-8.637091931657939d+0*t2*t27**3*t72-6.666666666666
     8         666d-1*t13*t52-4.1356699393293317d-1*t27*t28*t51+t45*(-2.
     9         4875621890547261d+0*t27*t50+4.6294812753686554d+0*t2*t27*
     :         t49+3.08632085024577d+0*t23*t44))*wght+Amat3(iq,D3_RB_RB_
     ;         RB)
            Cmat3(iq,D3_RA_RA_GAA) = (t22*(-9.305257363490997d-1*t1*(-3.
     1         093447392126889d-2*t38*t7-3.163753014675228d-3*t39*t5+1.9
     2         8637460967052d-2*t38+2.5512089741039945d-3*gammaaa*t37+1.
     3         6873349411601213d-2*t19*t20)-2.481401963597599d+0*t17*t54
     4         )-8.637091931657939d+0*t1*t35*t53*t69+t36*(4.629481275368
     5         6554d+0*t1*t21*t54-2.4875621890547261d+0*t41*t53)-4.13566
     6         99393293317d-1*t22*t42*t53)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-8.637091931657939d+0*t2*t44*t56*t
     1         72+t45*(4.6294812753686554d+0*t2*t27*t57-2.48756218905472
     2         61d+0*t50*t56)+t28*(-2.481401963597599d+0*t23*t57-9.30525
     3         7363490997d-1*t2*(-3.163753014675228d-3*t12*t48-3.0934473
     4         92126889d-2*t14*t47+1.98637460967052d-2*t47+2.55120897410
     5         39945d-3*gammabb*t46+1.6873349411601213d-2*t25*t26))-4.13
     6         56699393293317d-1*t28*t51*t56)*wght+Cmat3(iq,D3_RB_RB_GBB
     7         )
            Cmat3(iq,D3_RA_GAA_GAA) = (-8.637091931657939d+0*t1*t21*t69*
     1         t75+1.8487275744044215d-4*t21*t36*t59+t36*(2.314740637684
     2         3277d+0*t1*t53*t54-2.4875621890547261d+0*t53*t55)+2.97275
     3         393964231d-4*t22/rhoa**5)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-8.637091931657939d+0*t2*t27*t72*
     1         t76+1.8487275744044215d-4*t27*t45*t60+t45*(2.314740637684
     2         3277d+0*t2*t56*t57-2.4875621890547261d+0*t56*t58)+2.97275
     3         393964231d-4*t28/rhob**5)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t36*t59*(3.697455148808843d-4*t53
     1         *wght+1.8487275744044215d-4*t30*wght)-8.637091931657939d+
     2         0*t1*t30*t69*t75*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t45*t60*(3.697455148808843d-4*t56
     1         *wght+1.8487275744044215d-4*t32*wght)-8.637091931657939d+
     2         0*t2*t32*t72*t76*wght+Cmat3(iq,D3_GBB_GBB_GBB)
            Mmat3(iq,D3_RA_RA_TA) = (t22*(-9.305257363490997d-1*t1*(6.32
     1         7506029350458d-2*t19*t7+1.423688856603853d-2*t39*t6-4.745
     2         629522012844d-2*t20*t5-6.767868874683343d-3*gammaaa*t18)-
     3         2.481401963597599d+0*t17*t61)-8.637091931657939d+0*t1*t33
     4         *t35*t69+t36*(4.6294812753686554d+0*t1*t21*t61-2.48756218
     5         90547261d+0*t33*t41)-4.1356699393293317d-1*t22*t33*t42)*w
     6         ght+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_RB_TA) = Mmat3(iq,D3_RA_RB_TA)
            Mmat3(iq,D3_RA_RB_TB) = Mmat3(iq,D3_RA_RB_TB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (-8.637091931657939d+0*t2*t34*t44*t7
     1         2+t45*(4.6294812753686554d+0*t2*t27*t62-2.487562189054726
     2         1d+0*t34*t50)+t28*(-2.481401963597599d+0*t23*t62-9.305257
     3         363490997d-1*t2*(1.423688856603853d-2*t13*t48-4.745629522
     4         012844d-2*t12*t26+6.327506029350458d-2*t14*t25-6.76786887
     5         4683343d-3*gammabb*t24))-4.1356699393293317d-1*t28*t34*t5
     6         1)*wght+Mmat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (-8.637091931657939d+0*t1*t21*t33*t
     1         53*t69-7.2303974375762d-4*t21*t36*t63+t36*(2.314740637684
     2         3277d+0*t1*t53*t61-2.4875621890547261d+0*t33*t55)-8.71985
     3         9309716897d-4*t22*t59)*wght+Mmat3(iq,D3_RA_GAA_TA)
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
            Mmat3(iq,D3_RB_GBB_TB) = (-8.637091931657939d+0*t2*t27*t34*t
     1         56*t72-7.2303974375762d-4*t27*t45*t64+t45*(2.314740637684
     2         3277d+0*t2*t56*t62-2.4875621890547261d+0*t34*t58)-8.71985
     3         9309716897d-4*t28*t60)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_RA_TA_TA) = t36*t63*(-7.2303974375762d-4*t53*wgh
     1         t-7.2303974375762d-4*t30*wght)-8.637091931657939d+0*t1*t3
     2         0*t33*t53*t69*wght+1.8487275744044215d-4*t33*t36*t59*wght
     3         +Mmat3(iq,D3_RA_TA_TA)
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
            Mmat3(iq,D3_GAA_GBB_TB) = t45*t64*(-7.2303974375762d-4*t56*w
     1         ght-7.2303974375762d-4*t32*wght)-8.637091931657939d+0*t2*
     2         t32*t34*t56*t72*wght+1.8487275744044215d-4*t34*t45*t60*wg
     3         ht+Mmat3(iq,D3_GAA_GBB_TB)
            Mmat3(iq,D3_GAB_GAB_TA) = -8.637091931657939d+0*t1*t30*t65*t
     1         69*wght+3.2536788469092914d-3*t30*t36*t66*wght-1.44607948
     2         751524d-3*t33*t36*t63*wght+Mmat3(iq,D3_GAB_GAB_TA)
            Mmat3(iq,D3_GAB_GAB_TB) = Mmat3(iq,D3_GAB_GAB_TB)
            Mmat3(iq,D3_GAB_GBB_TA) = Mmat3(iq,D3_GAB_GBB_TA)
            Mmat3(iq,D3_GAB_GBB_TB) = Mmat3(iq,D3_GAB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = -8.637091931657939d+0*t2*t32*t67*t7
     1         2*wght+3.2536788469092914d-3*t32*t45*t68*wght-1.446079487
     2         51524d-3*t34*t45*t64*wght+Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_GAB_TA_TA) = -8.637091931657939d+0*t1*t33**3*t69
     1         *wght+9.761036540727874d-3*t33*t36*t66*wght+Mmat3(iq,D3_G
     2         AB_TA_TA)
            Mmat3(iq,D3_GAB_TA_TB) = Mmat3(iq,D3_GAB_TA_TB)
            Mmat3(iq,D3_GAB_TB_TB) = Mmat3(iq,D3_GAB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = -8.637091931657939d+0*t2*t34**3*t72
     1         *wght+9.761036540727874d-3*t34*t45*t68*wght+Mmat3(iq,D3_G
     2         BB_TA_TA)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = gammaaa**2
            t3 = 1/rhoa**5.333333333333333d+0
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taua-1.3712756538350462d-3*gamm
     1         aaa*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammaaa*t4*t6+2.031519487163032d-3*gammaa
     2         a*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhoa**3.333333333333333d-1
            t10 = 1/rhoa**6.333333333333333d+0
            t11 = 1/rhoa**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammaaa*t11-1.6455307846020562d-1
     1         *t4*taua
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         aaa*t11*t6-2.966018451258026d-3*gammaaa*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammaaa*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammaaa*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaa
     1         a/rhoa**4.333333333333333d+0
            t18 = t13**2
            t19 = 1/t7**3
            t20 = 1/rhoa**7.333333333333333d+0
            t21 = 1/rhoa**4.666666666666667d+0
            t22 = 4.38808209227215d-1*t11*taua-1.3408028615276008d-2*gam
     1         maaa*t21
            t23 = 1.4419753086419754d-1*t22*t6-2.9001069301189586d-2*gam
     1         maaa*t21*t6-2.966018451258026d-3*gammaaa*t22*t4+1.9863746
     2         0967052d-2*gammaaa*t21+1.2069135102256172d-3*t2*t20+1.441
     3         9753086419754d-1*t12**2+1.5818765073376137d-2*gammaaa*t11
     4         *t12
            t24 = -2.481401963597599d+0*t13*t9-9.305257363490997d-1*t1*t
     1         23
            t25 = 1/rhoa**6.666666666666666d-1
            t26 = -4.1356699393293317d-1*t8-4.1356699393293317d-1
            t27 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            t28 = 8.436674705800608d-3*t11*t6-3.163753014675228d-3*t12*t
     1         4-5.417385299101418d-3*t11-4.0282246959536755d-4*gammaaa*
     2         t10
            t29 = -1.2407009817987995d+0*t27*t9-9.305257363490997d-1*t1*
     1         t28
            t30 = 1/rhoa**4
            t31 = -2.372814761006422d-2*t4*t6+1.423688856603853d-2*t12*t
     1         5+1.2689754140031267d-3*gammaaa*t3
            t32 = -1.2407009817987995d+0*t17*t9-9.305257363490997d-1*t1*
     1         t31
            t33 = 1/rhoa**3
            t34 = t17**2
            t35 = 1/rhoa**2
            t36 = 1/t7**4
            t37 = 1/rhoa**5.666666666666667d+0
            t38 = 6.257080020462137d-2*gammaaa*t37-1.6089634338331216d+0
     1         *t21*taua
            t39 = t27**2
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)
            Amat2(iq,D2_RA_RA) = (t25*t26+t14*t24+2.3147406376843277d+0*
     1         t1*t18*t19)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = (t14*t29+2.3147406376843277d+0*t1*t13*
     1         t19*t27)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -7.431884849105774d-5*t14*t30*wght+2.
     1         3147406376843277d+0*t1*t16*t19*t27*wght+Cmat2(iq,D2_GAA_G
     2         AA)
            Mmat2(iq,D2_RA_TA) = (t14*t32+2.3147406376843277d+0*t1*t13*t
     1         17*t19)*wght+Mmat2(iq,D2_RA_TA)
            Mmat2(iq,D2_RA_TB) = Mmat2(iq,D2_RA_TB)
            Mmat2(iq,D2_GAA_TA) = 2.906619769905633d-4*t14*t33*wght+2.31
     1         47406376843277d+0*t1*t16*t17*t19*wght+Mmat2(iq,D2_GAA_TA)
            Mmat2(iq,D2_GAA_TB) = Mmat2(iq,D2_GAA_TB)
            Mmat2(iq,D2_TA_TA) = -1.3079788964575353d-3*t14*t35*wght+2.3
     1         147406376843277d+0*t1*t19*t34*wght+Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RA_RA_RA) = (t14*(-3.7221029453963983d+0*t23*t9-
     1         9.305257363490997d-1*t1*(1.4419753086419754d-1*t38*t6+1.3
     2         53383234055514d-1*gammaaa*t37*t6-2.966018451258026d-3*gam
     3         maaa*t38*t4-9.269748178462428d-2*gammaaa*t37+4.3259259259
     4         259264d-1*t12*t22+2.3728147610064204d-2*gammaaa*t11*t22-8
     5         .700320790356875d-2*gammaaa*t12*t21-8.85069907498786d-3*t
     6         2/rhoa**8.333333333333334d+0)-8.271339878658663d-1*t13*t2
     7         5)+t19*(3.08632085024577d+0*t18*t9-2.4875621890547261d+0*
     8         t13*t24+4.6294812753686554d+0*t1*t13*t23)-6.6666666666666
     9         66d-1*t26*t5-8.637091931657939d+0*t1*t13**3*t36-4.1356699
     :         393293317d-1*t13*t14*t25)*wght+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = (t14*(-2.481401963597599d+0*t28*t9-
     1         9.305257363490997d-1*t1*(-3.093447392126889d-2*t21*t6-3.1
     2         63753014675228d-3*t22*t4+1.98637460967052d-2*t21+2.551208
     3         9741039945d-3*gammaaa*t20+1.6873349411601213d-2*t11*t12))
     4         -8.637091931657939d+0*t1*t18*t27*t36+t19*(4.6294812753686
     5         554d+0*t1*t13*t28-2.4875621890547261d+0*t24*t27)-4.135669
     6         9393293317d-1*t14*t25*t27)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = (-8.637091931657939d+0*t1*t13*t36*
     1         t39+1.8487275744044215d-4*t13*t19*t30+t19*(2.314740637684
     2         3277d+0*t1*t27*t28-2.4875621890547261d+0*t27*t29)+2.97275
     3         393964231d-4*t14/rhoa**5)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = t19*t30*(3.697455148808843d-4*t27
     1         *wght+1.8487275744044215d-4*t16*wght)-8.637091931657939d+
     2         0*t1*t16*t36*t39*wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Mmat3(iq,D3_RA_RA_TA) = (t14*(-2.481401963597599d+0*t31*t9-9
     1         .305257363490997d-1*t1*(6.327506029350458d-2*t11*t6+1.423
     2         688856603853d-2*t22*t5-4.745629522012844d-2*t12*t4-6.7678
     3         68874683343d-3*gammaaa*t10))-8.637091931657939d+0*t1*t17*
     4         t18*t36+t19*(4.6294812753686554d+0*t1*t13*t31-2.487562189
     5         0547261d+0*t17*t24)-4.1356699393293317d-1*t14*t17*t25)*wg
     6         ht+Mmat3(iq,D3_RA_RA_TA)
            Mmat3(iq,D3_RA_RA_TB) = Mmat3(iq,D3_RA_RA_TB)
            Mmat3(iq,D3_RA_GAA_TA) = (-8.637091931657939d+0*t1*t13*t17*t
     1         27*t36-7.2303974375762d-4*t13*t19*t33+t19*(2.314740637684
     2         3277d+0*t1*t27*t31-2.4875621890547261d+0*t17*t29)-8.71985
     3         9309716897d-4*t14*t30)*wght+Mmat3(iq,D3_RA_GAA_TA)
            Mmat3(iq,D3_RA_GAA_TB) = Mmat3(iq,D3_RA_GAA_TB)
            Mmat3(iq,D3_GAA_GAA_TA) = t19*t33*(-7.2303974375762d-4*t27*w
     1         ght-7.2303974375762d-4*t16*wght)-8.637091931657939d+0*t1*
     2         t16*t17*t27*t36*wght+1.8487275744044215d-4*t17*t19*t30*wg
     3         ht+Mmat3(iq,D3_GAA_GAA_TA)
            Mmat3(iq,D3_GAA_GAA_TB) = Mmat3(iq,D3_GAA_GAA_TB)
            Mmat3(iq,D3_RA_TA_TA) = (-8.637091931657939d+0*t1*t13*t34*t3
     1         6+3.2536788469092914d-3*t13*t19*t35+2.6159577929150707d-3
     2         *t14*t33+t19*(2.3147406376843277d+0*t1*t17*t31-2.48756218
     3         90547261d+0*t17*t32))*wght+Mmat3(iq,D3_RA_TA_TA)
            Mmat3(iq,D3_RA_TA_TB) = Mmat3(iq,D3_RA_TA_TB)
            Mmat3(iq,D3_RA_TB_TB) = Mmat3(iq,D3_RA_TB_TB)
            Mmat3(iq,D3_GAA_TA_TA) = -8.637091931657939d+0*t1*t16*t34*t3
     1         6*wght+3.2536788469092914d-3*t16*t19*t35*wght-1.446079487
     2         51524d-3*t17*t19*t33*wght+Mmat3(iq,D3_GAA_TA_TA)
            Mmat3(iq,D3_GAA_TA_TB) = Mmat3(iq,D3_GAA_TA_TB)
            Mmat3(iq,D3_GAA_TB_TB) = Mmat3(iq,D3_GAA_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = -8.637091931657939d+0*t1*t17**3*t36*
     1         wght+9.761036540727874d-3*t17*t19*t35*wght+Mmat3(iq,D3_TA
     2         _TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = Mmat3(iq,D3_TB_TB_TB)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = gammabb**2
            t3 = 1/rhob**5.333333333333333d+0
            t4 = 1/rhob**2.6666666666666666d+0
            t5 = 1/rhob**1.6666666666666669d+0
            t6 = 9.873184707612337d-2*t5*taub-1.3712756538350462d-3*gamm
     1         abb*t4-4.5d-1
            t7 = 1.2437810945273631d+0*(7.209876543209877d-2*t6**2-2.966
     1         018451258026d-3*gammabb*t4*t6+2.031519487163032d-3*gammab
     2         b*t4+3.573099207904788d-5*t2*t3)+1.0d+0
            t8 = 8.04d-1-8.04d-1/t7
            t9 = rhob**3.333333333333333d-1
            t10 = 1/rhob**6.333333333333333d+0
            t11 = 1/rhob**3.6666666666666664d+0
            t12 = 3.656735076893457d-3*gammabb*t11-1.6455307846020562d-1
     1         *t4*taub
            t13 = 1.4419753086419754d-1*t12*t6+7.909382536688069d-3*gamm
     1         abb*t11*t6-2.966018451258026d-3*gammabb*t12*t4-1.90565291
     2         08825535d-4*t10*t2-5.417385299101418d-3*gammabb*t11
            t14 = 1/t7**2
            t15 = 7.552921304913142d-5*gammabb*t3
            t16 = t4*(2.031519487163032d-3-3.163753014675228d-3*t6)+t15
            t17 = 1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammab
     1         b/rhob**4.333333333333333d+0
            t18 = t13**2
            t19 = 1/t7**3
            t20 = 1/rhob**7.333333333333333d+0
            t21 = 1/rhob**4.666666666666667d+0
            t22 = 4.38808209227215d-1*t11*taub-1.3408028615276008d-2*gam
     1         mabb*t21
            t23 = 1.4419753086419754d-1*t22*t6-2.9001069301189586d-2*gam
     1         mabb*t21*t6-2.966018451258026d-3*gammabb*t22*t4+1.9863746
     2         0967052d-2*gammabb*t21+1.2069135102256172d-3*t2*t20+1.441
     3         9753086419754d-1*t12**2+1.5818765073376137d-2*gammabb*t11
     4         *t12
            t24 = -2.481401963597599d+0*t13*t9-9.305257363490997d-1*t1*t
     1         23
            t25 = 1/rhob**6.666666666666666d-1
            t26 = -4.1356699393293317d-1*t8-4.1356699393293317d-1
            t27 = -3.163753014675228d-3*t4*t6+2.031519487163032d-3*t4+t1
     1         5
            t28 = 8.436674705800608d-3*t11*t6-3.163753014675228d-3*t12*t
     1         4-5.417385299101418d-3*t11-4.0282246959536755d-4*gammabb*
     2         t10
            t29 = -1.2407009817987995d+0*t27*t9-9.305257363490997d-1*t1*
     1         t28
            t30 = 1/rhob**4
            t31 = -2.372814761006422d-2*t4*t6+1.423688856603853d-2*t12*t
     1         5+1.2689754140031267d-3*gammabb*t3
            t32 = -1.2407009817987995d+0*t17*t9-9.305257363490997d-1*t1*
     1         t31
            t33 = 1/rhob**3
            t34 = t17**2
            t35 = 1/rhob**2
            t36 = 1/t7**4
            t37 = 1/rhob**5.666666666666667d+0
            t38 = 6.257080020462137d-2*gammabb*t37-1.6089634338331216d+0
     1         *t21*taub
            t39 = t27**2
            fnc(iq) = (-9.305257363490997d-1*t1*t8-9.305257363490997d-1*
     1         t1)*wght+fnc(iq)
            Amat(iq,D1_RB) = (-1.2407009817987995d+0*t8*t9-1.24070098179
     1         87995d+0*t9-9.305257363490997d-1*t1*t13*t14)*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490997d-1*t1*t1
     1         4*t16*wght
            Mmat(iq,D1_TB) = Mmat(iq,D1_TB)-9.305257363490997d-1*t1*t14*
     1         t17*wght
            Mmat(iq,D1_TA) = Mmat(iq,D1_TA)
            Amat2(iq,D2_RB_RB) = (t25*t26+t14*t24+2.3147406376843277d+0*
     1         t1*t18*t19)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = (t14*t29+2.3147406376843277d+0*t1*t13*
     1         t19*t27)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -7.431884849105774d-5*t14*t30*wght+2.
     1         3147406376843277d+0*t1*t16*t19*t27*wght+Cmat2(iq,D2_GBB_G
     2         BB)
            Mmat2(iq,D2_RB_TA) = Mmat2(iq,D2_RB_TA)
            Mmat2(iq,D2_RB_TB) = (t14*t32+2.3147406376843277d+0*t1*t13*t
     1         17*t19)*wght+Mmat2(iq,D2_RB_TB)
            Mmat2(iq,D2_GBB_TA) = Mmat2(iq,D2_GBB_TA)
            Mmat2(iq,D2_GBB_TB) = 2.906619769905633d-4*t14*t33*wght+2.31
     1         47406376843277d+0*t1*t16*t17*t19*wght+Mmat2(iq,D2_GBB_TB)
            Mmat2(iq,D2_TA_TA) = Mmat2(iq,D2_TA_TA)
            Mmat2(iq,D2_TA_TB) = Mmat2(iq,D2_TA_TB)
            Mmat2(iq,D2_TB_TB) = -1.3079788964575353d-3*t14*t35*wght+2.3
     1         147406376843277d+0*t1*t19*t34*wght+Mmat2(iq,D2_TB_TB)
            Amat3(iq,D3_RB_RB_RB) = (t14*(-3.7221029453963983d+0*t23*t9-
     1         9.305257363490997d-1*t1*(1.4419753086419754d-1*t38*t6+1.3
     2         53383234055514d-1*gammabb*t37*t6-2.966018451258026d-3*gam
     3         mabb*t38*t4-9.269748178462428d-2*gammabb*t37+4.3259259259
     4         259264d-1*t12*t22+2.3728147610064204d-2*gammabb*t11*t22-8
     5         .700320790356875d-2*gammabb*t12*t21-8.85069907498786d-3*t
     6         2/rhob**8.333333333333334d+0)-8.271339878658663d-1*t13*t2
     7         5)+t19*(3.08632085024577d+0*t18*t9-2.4875621890547261d+0*
     8         t13*t24+4.6294812753686554d+0*t1*t13*t23)-6.6666666666666
     9         66d-1*t26*t5-8.637091931657939d+0*t1*t13**3*t36-4.1356699
     :         393293317d-1*t13*t14*t25)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = (t14*(-2.481401963597599d+0*t28*t9-
     1         9.305257363490997d-1*t1*(-3.093447392126889d-2*t21*t6-3.1
     2         63753014675228d-3*t22*t4+1.98637460967052d-2*t21+2.551208
     3         9741039945d-3*gammabb*t20+1.6873349411601213d-2*t11*t12))
     4         -8.637091931657939d+0*t1*t18*t27*t36+t19*(4.6294812753686
     5         554d+0*t1*t13*t28-2.4875621890547261d+0*t24*t27)-4.135669
     6         9393293317d-1*t14*t25*t27)*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = (-8.637091931657939d+0*t1*t13*t36*
     1         t39+1.8487275744044215d-4*t13*t19*t30+t19*(2.314740637684
     2         3277d+0*t1*t27*t28-2.4875621890547261d+0*t27*t29)+2.97275
     3         393964231d-4*t14/rhob**5)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = t19*t30*(3.697455148808843d-4*t27
     1         *wght+1.8487275744044215d-4*t16*wght)-8.637091931657939d+
     2         0*t1*t16*t36*t39*wght+Cmat3(iq,D3_GBB_GBB_GBB)
            Mmat3(iq,D3_RB_RB_TA) = Mmat3(iq,D3_RB_RB_TA)
            Mmat3(iq,D3_RB_RB_TB) = (t14*(-2.481401963597599d+0*t31*t9-9
     1         .305257363490997d-1*t1*(6.327506029350458d-2*t11*t6+1.423
     2         688856603853d-2*t22*t5-4.745629522012844d-2*t12*t4-6.7678
     3         68874683343d-3*gammabb*t10))-8.637091931657939d+0*t1*t17*
     4         t18*t36+t19*(4.6294812753686554d+0*t1*t13*t31-2.487562189
     5         0547261d+0*t17*t24)-4.1356699393293317d-1*t14*t17*t25)*wg
     6         ht+Mmat3(iq,D3_RB_RB_TB)
            Mmat3(iq,D3_RB_GBB_TA) = Mmat3(iq,D3_RB_GBB_TA)
            Mmat3(iq,D3_RB_GBB_TB) = (-8.637091931657939d+0*t1*t13*t17*t
     1         27*t36-7.2303974375762d-4*t13*t19*t33+t19*(2.314740637684
     2         3277d+0*t1*t27*t31-2.4875621890547261d+0*t17*t29)-8.71985
     3         9309716897d-4*t14*t30)*wght+Mmat3(iq,D3_RB_GBB_TB)
            Mmat3(iq,D3_GBB_GBB_TA) = Mmat3(iq,D3_GBB_GBB_TA)
            Mmat3(iq,D3_GBB_GBB_TB) = t19*t33*(-7.2303974375762d-4*t27*w
     1         ght-7.2303974375762d-4*t16*wght)-8.637091931657939d+0*t1*
     2         t16*t17*t27*t36*wght+1.8487275744044215d-4*t17*t19*t30*wg
     3         ht+Mmat3(iq,D3_GBB_GBB_TB)
            Mmat3(iq,D3_RB_TA_TA) = Mmat3(iq,D3_RB_TA_TA)
            Mmat3(iq,D3_RB_TA_TB) = Mmat3(iq,D3_RB_TA_TB)
            Mmat3(iq,D3_RB_TB_TB) = (-8.637091931657939d+0*t1*t13*t34*t3
     1         6+3.2536788469092914d-3*t13*t19*t35+2.6159577929150707d-3
     2         *t14*t33+t19*(2.3147406376843277d+0*t1*t17*t31-2.48756218
     3         90547261d+0*t17*t32))*wght+Mmat3(iq,D3_RB_TB_TB)
            Mmat3(iq,D3_GBB_TA_TA) = Mmat3(iq,D3_GBB_TA_TA)
            Mmat3(iq,D3_GBB_TA_TB) = Mmat3(iq,D3_GBB_TA_TB)
            Mmat3(iq,D3_GBB_TB_TB) = -8.637091931657939d+0*t1*t16*t34*t3
     1         6*wght+3.2536788469092914d-3*t16*t19*t35*wght-1.446079487
     2         51524d-3*t17*t19*t33*wght+Mmat3(iq,D3_GBB_TB_TB)
            Mmat3(iq,D3_TA_TA_TA) = Mmat3(iq,D3_TA_TA_TA)
            Mmat3(iq,D3_TA_TA_TB) = Mmat3(iq,D3_TA_TA_TB)
            Mmat3(iq,D3_TA_TB_TB) = Mmat3(iq,D3_TA_TB_TB)
            Mmat3(iq,D3_TB_TB_TB) = -8.637091931657939d+0*t1*t17**3*t36*
     1         wght+9.761036540727874d-3*t17*t19*t35*wght+Mmat3(iq,D3_TB
     2         _TB_TB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
