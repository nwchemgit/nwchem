C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_b97.F 
C> The nwxcm_x_b97 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,\sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
      subroutine nwxcm_x_b97(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhoa**6.333333333333333d+0
            t13 = -1.8610514726981986d+0*t10*(4.0d-3*t3*t7*t8+t5*t6*(3.2
     1         d-5*gammaaa*t1-1.6d-5*gammaaa*t7)-1.28d-7*t1*t11*t2/rhoa*
     2         *8)*wght
            fnc(iq) = fnc(iq)-1.8610514726981986d+0*t10*t9*wght
            Amat(iq,D1_RA) = -2.481401963597598d+0*rhoa**3.3333333333333
     1         33d-1*t9*wght-1.8610514726981986d+0*t10*(-1.0666666666666
     2         666d-2*gammaaa*t7*t8/rhoa**3.6666666666666664d+0+4.266666
     3         666666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t2*
     4         t5+3.413333333333333d-7*gammaaa**3*t1*t11/rhoa**9)*wght+A
     5         mat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t13+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t13+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(2)
            t2 = param(4)
            t3 = gammaaa**2
            t4 = 1/rhoa**2.6666666666666666d+0
            t5 = 4.0d-3*gammaaa*t4+1.0d+0
            t6 = 1/t5**2
            t7 = 1/rhoa**5.333333333333333d+0
            t8 = param(3)
            t9 = 1/t5
            t10 = 4.0d-3*gammaaa*t4*t8*t9+1.6d-5*t2*t3*t6*t7+t1
            t11 = rhoa**1.3333333333333333d+0
            t12 = gammabb**2
            t13 = 1/rhob**2.6666666666666666d+0
            t14 = 4.0d-3*gammabb*t13+1.0d+0
            t15 = 1/t14**2
            t16 = 1/rhob**5.333333333333333d+0
            t17 = 1/t14
            t18 = 4.0d-3*gammabb*t13*t17*t8+1.6d-5*t12*t15*t16*t2+t1
            t19 = rhob**1.3333333333333333d+0
            t20 = 1/t5**3
            t21 = 1/rhoa**6.333333333333333d+0
            t22 = 1/t14**3
            t23 = 1/rhob**6.333333333333333d+0
            fnc(iq) = 1.0d+0*(-9.305257363490993d-1*t18*t19-9.3052573634
     1         90993d-1*t10*t11)*wght+fnc(iq)
            Amat(iq,D1_RA) = 1.0d+0*(-9.305257363490993d-1*t11*(-1.06666
     1         66666666666d-2*gammaaa*t8*t9/rhoa**3.6666666666666664d+0+
     2         4.266666666666667d-5*t21*t3*t6*t8-8.533333333333333d-5*t2
     3         *t21*t3*t6+3.413333333333333d-7*gammaaa**3*t2*t20/rhoa**9
     4         )-1.240700981798799d+0*rhoa**3.333333333333333d-1*t10)*wg
     5         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 1.0d+0*(-9.305257363490993d-1*t19*(4.266666
     1         666666667d-5*t12*t15*t23*t8-1.0666666666666666d-2*gammabb
     2         *t17*t8/rhob**3.6666666666666664d+0-8.533333333333333d-5*
     3         t12*t15*t2*t23+3.413333333333333d-7*gammabb**3*t2*t22/rho
     4         b**9)-1.240700981798799d+0*rhob**3.333333333333333d-1*t18
     5         )*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t11*(
     1         4.0d-3*t4*t8*t9+t6*t7*(3.2d-5*gammaaa*t2-1.6d-5*gammaaa*t
     2         8)-1.28d-7*t2*t20*t3/rhoa**8)*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t19*(
     1         t15*t16*(3.2d-5*gammabb*t2-1.6d-5*gammabb*t8)+4.0d-3*t13*
     2         t17*t8-1.28d-7*t12*t2*t22/rhob**8)*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(4)
            t2 = gammaaa**2
            t3 = 1/rhoa**2.6666666666666666d+0
            t4 = 4.0d-3*gammaaa*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhoa**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammaaa*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhoa**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhoa**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RA) = -1.240700981798799d+0*rhoa**3.3333333333333
     1         33d-1*t9*wght-9.305257363490993d-1*t10*(-1.06666666666666
     2         66d-2*gammaaa*t7*t8/rhoa**3.6666666666666664d+0+4.2666666
     3         66666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t2*t
     4         5+3.413333333333333d-7*gammaaa**3*t1*t11/rhoa**9)*wght+Am
     5         at(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-9.305257363490993d-1*t10*(
     1         4.0d-3*t3*t7*t8+t5*t6*(3.2d-5*gammaaa*t1-1.6d-5*gammaaa*t
     2         7)-1.28d-7*t1*t11*t2/rhoa**8)*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(4)
            t2 = gammabb**2
            t3 = 1/rhob**2.6666666666666666d+0
            t4 = 4.0d-3*gammabb*t3+1.0d+0
            t5 = 1/t4**2
            t6 = 1/rhob**5.333333333333333d+0
            t7 = param(3)
            t8 = 1/t4
            t9 = 4.0d-3*gammabb*t3*t7*t8+1.6d-5*t1*t2*t5*t6+param(2)
            t10 = rhob**1.3333333333333333d+0
            t11 = 1/t4**3
            t12 = 1/rhob**6.333333333333333d+0
            fnc(iq) = fnc(iq)-9.305257363490993d-1*t10*t9*wght
            Amat(iq,D1_RB) = -1.240700981798799d+0*rhob**3.3333333333333
     1         33d-1*t9*wght-9.305257363490993d-1*t10*(-1.06666666666666
     2         66d-2*gammabb*t7*t8/rhob**3.6666666666666664d+0+4.2666666
     3         66666667d-5*t12*t2*t5*t7-8.533333333333333d-5*t1*t12*t2*t
     4         5+3.413333333333333d-7*gammabb**3*t1*t11/rhob**9)*wght+Am
     5         at(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-9.305257363490993d-1*t10*(
     1         4.0d-3*t3*t7*t8+t5*t6*(3.2d-5*gammabb*t1-1.6d-5*gammabb*t
     2         7)-1.28d-7*t1*t11*t2/rhob**8)*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,\sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
      subroutine nwxcm_x_b97_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_b97 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)\\\\ 
C>   {\it t_2} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& 0.004\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0\\\\ 
C>   {\it t_5} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& 0.004\,{\it t_6}\,\sigma_{\beta\beta}+1.0\\\\ 
C>   {\it t_8} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_9} &=& 0.004\,{\it t_8}\,\sigma_{ss}+1.0\\\\ 
C>   f &=& 1.0\,\left(-.9305257363490993\,\rho_\beta^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_6}\,\sigma_{\beta\beta}}\over{{\it t_7}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}\,{\it t_7}^2}}+{\it t_1}\right)-.9305257363490993\,\rho_\alpha^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{{\it t_4}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}\,{\it t_4}^2}}+{\it t_1}\right)\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -.9305257363490993\,\rho_s^{{{4}\over{3}}}\,\left({{0.004\,{\it t_5}\,{\it t_8}\,\sigma_{ss}}\over{{\it t_9}}}+{{1.6 \times 10^{-5}\,{\it t_2}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}\,{\it t_9}^2}}+{\it t_1}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 107, 8554 (1997)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.475007 "> 
C> 10.1063/1.475007 </a> 
C> 
      subroutine nwxcm_x_b97_d3(param,tol_rho,ipol,nq,wght, 
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
