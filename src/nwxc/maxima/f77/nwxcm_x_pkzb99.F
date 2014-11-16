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
C>   {\it t_4} &=& {{.09873184707612337\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{.09873184707612337\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{.09873184707612337\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_8}\,\sigma_{ss}-0.45\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_6}^2-.002966018451258026\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_4}^2-.002966018451258026\,{\it t_3}\,\sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_3}\,\sigma_{\alpha\alpha}\right)+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_7}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_9}^2-.002966018451258026\,{\it t_8}\,\sigma_{ss}\,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)+1.0}}\right)-.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82, 2544 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
      subroutine nwxcm_x_pkzb99(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Cmat,Mmat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
            t13 = -1.8610514726981994d+0*t1*t12*(t4*(2.031519487163032d-
     1         3-3.163753014675228d-3*t6)+7.552921304913142d-5*gammaaa*t
     2         3)*wght
            fnc(iq) = (-1.8610514726981994d+0*t1*t8-1.8610514726981994d+
     1         0*t1)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-2.481401963597599d+0*t8*t9-2.481401963597
     1         599d+0*t9-1.8610514726981994d+0*t1*t12*(1.441975308641975
     2         4d-1*t11*t6+7.909382536688069d-3*gammaaa*t10*t6-2.9660184
     3         51258026d-3*gammaaa*t11*t4-1.9056529108825535d-4*t2/rhoa*
     4         *6.333333333333333d+0-5.417385299101418d-3*gammaaa*t10))*
     5         wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t13+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t13+Cmat(iq,D1_GAB)
            Amat(iq,D1_TA) = Amat(iq,D1_TA)-1.8610514726981994d+0*t1*t12
     1         *(1.423688856603853d-2*t5*t6-2.928404801545677d-4*gammaaa
     2         /rhoa**4.333333333333333d+0)*wght
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
C>   {\it t_4} &=& {{.09873184707612337\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{.09873184707612337\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{.09873184707612337\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_8}\,\sigma_{ss}-0.45\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_6}^2-.002966018451258026\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_4}^2-.002966018451258026\,{\it t_3}\,\sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_3}\,\sigma_{\alpha\alpha}\right)+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_7}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_9}^2-.002966018451258026\,{\it t_8}\,\sigma_{ss}\,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)+1.0}}\right)-.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82, 2544 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
      subroutine nwxcm_x_pkzb99_d2(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
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
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
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
C>   {\it t_4} &=& {{.09873184707612337\,\tau_\alpha}\over{\rho_\alpha^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_3}\,\sigma_{\alpha\alpha}-0.45\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& {{.09873184707612337\,\tau_\beta}\over{\rho_\beta^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_5}\,\sigma_{\beta\beta}-0.45\\\\ 
C>   {\it t_7} &=& \rho_s^{{{4}\over{3}}}\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& {{.09873184707612337\,\tau_s}\over{\rho_s^{{{5}\over{3}}}}}-.001371275653835046\,{\it t_8}\,\sigma_{ss}-0.45\\\\ 
C>   f &=& -.9305257363490997\,{\it t_2}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_6}^2-.002966018451258026\,{\it t_5}\,\sigma_{\beta\beta}\,{\it t_6}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_5}\,\sigma_{\beta\beta}\right)+1.0}}\right)-.9305257363490997\,{\it t_1}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_4}^2-.002966018451258026\,{\it t_3}\,\sigma_{\alpha\alpha}\,{\it t_4}+{{3.573099207904788 \times 10^{-5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_3}\,\sigma_{\alpha\alpha}\right)+1.0}}\right)-.9305257363490997\,{\it t_2}-.9305257363490997\,{\it t_1}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490997\,{\it t_7}\,\left(0.804-{{0.804}\over{1.243781094527363\,\left(.07209876543209877\,{\it t_9}^2-.002966018451258026\,{\it t_8}\,\sigma_{ss}\,{\it t_9}+{{3.573099207904788 \times 10^{-5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+.002031519487163032\,{\it t_8}\,\sigma_{ss}\right)+1.0}}\right)-.9305257363490997\,{\it t_7}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, S Kurth, A Zupan, P Blaha, Phys.Rev.Lett. 82, 2544 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.82.2544 "> 
C> 10.1103/PhysRevLett.82.2544 </a> 
C> 
      subroutine nwxcm_x_pkzb99_d3(param,tol_rho,ipol,nq,wght, 
     +rho,rgamma,tau,fnc,Amat,Amat2,Amat3, 
     +Cmat,Cmat2,Cmat3,Mmat,Mmat2,Mmat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
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
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          taua    = 0.5d0*tau(iq,T_T) 
          if (rhoa.gt.tol_rho) then 
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
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
