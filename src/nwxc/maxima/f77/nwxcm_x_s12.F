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
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)+1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_s12(param,tol_rho,ipol,nq,wght, 
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
            t15 = -1.8610514726982003d+0*t1*t13*t8*t9*wght/t12-1.8610514
     1         726982003d+0*t1*t11*t12*t14*(t5*t6+2*gammaaa*t2*t4)*wght
            fnc(iq) = fnc(iq)-1.8610514726982003d+0*t1*t11*t12*t8*wght
            Amat(iq,D1_RA) = 4.9628039271952007d+0*gammaaa*t1*t13*t8*t9*
     1         wght/rhoa**2.3333333333333334d+0-2.4814019635976003d+0*rh
     2         oa**3.333333333333333d-1*t1*t11*t8*wght-1.861051472698200
     3         3d+0*t1*t11*t12*t14*(-2.6666666666666666d+0*gammaaa*t5/rh
     4         oa**3.6666666666666664d+0-5.333333333333333d+0*t2*t3/rhoa
     5         **6.333333333333333d+0)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t15+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t15+Cmat(iq,D1_GAB)
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
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)+1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_s12_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_s12 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(2\right)-{\it param}\left(1\right)+1.0\\\\ 
C>   {\it t_2} &=& {\it param}\left(5\right)\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& {\it param}\left(3\right)\\\\ 
C>   {\it t_5} &=& {\it param}\left(4\right)\\\\ 
C>   {\it t_6} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_7} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -0.9305257363491002\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\beta\beta}^2}\over{\rho_\beta^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_6}\,\sigma_{\beta\beta}+1.0}}\right)-0.9305257363491002\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{\alpha\alpha}^2}\over{\rho_\alpha^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_3}\,\sigma_{\alpha\alpha}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491002\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{\it t_2}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\,\left(1.0-{{1.0}\over{{{{\it t_5}\,\sigma_{ss}^2}\over{\rho_s^{{{16}\over{3}}}}}+{\it t_4}\,{\it t_7}\,\sigma_{ss}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] M Swart, Chem.Phys.Lett. 580, 166 (2013)  , DOI: 
C> <a href="http://dx.doi.org/10.1016/j.cplett.2013.06.045 "> 
C> 10.1016/j.cplett.2013.06.045 </a> 
C> 
      subroutine nwxcm_x_s12_d3(param,tol_rho,ipol,nq,wght, 
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
