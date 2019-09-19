C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_pbe.F 
C> The nwxcm_x_pbe functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\right)-0.9305257363491\,{\it t_1}
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pbe(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020564d-2*gammaaa*t2/(rhoa**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhoa**1.3333333333333333d+0
            t6 = 1/t3**2
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t4*t5*wght
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t2*t6/rhoa**
     1         2.3333333333333334d+0-1.2407009817987999d+0*rhoa**3.33333
     2         3333333333d-1*t1*t4)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t2*t6
     1         *wght/t5
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1.6455307846020564d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = rhob**1.3333333333333333d+0
            t10 = 1/t4**2
            t11 = 1/t7**2
            fnc(iq) = (-9.305257363491d-1*t1*t8*t9-9.305257363491d-1*t1*
     1         t5*t6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t10*t3/rhoa*
     1         *2.3333333333333334d+0-1.2407009817987999d+0*rhoa**3.3333
     2         33333333333d-1*t1*t5)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.0832233200718415d-2*gammabb*t11*t3/rhob*
     1         *2.3333333333333334d+0-1.2407009817987999d+0*rhob**3.3333
     2         33333333333d-1*t1*t8)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t10*t
     1         3*wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t11*t
     1         3*wght/t9
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020564d-2*gammaaa*t2/(rhoa**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhoa**1.3333333333333333d+0
            t6 = 1/t3**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t4*t5*wght
            Amat(iq,D1_RA) = 4.0832233200718415d-2*gammaaa*t2*t6*wght/rh
     1         oa**2.3333333333333334d+0-1.2407009817987999d+0*rhoa**3.3
     2         33333333333333d-1*t1*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t2*t6
     1         *wght/t5
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020564d-2*gammabb*t2/(rhob**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhob**1.3333333333333333d+0
            t6 = 1/t3**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t4*t5*wght
            Amat(iq,D1_RB) = 4.0832233200718415d-2*gammabb*t2*t6*wght/rh
     1         ob**2.3333333333333334d+0-1.2407009817987999d+0*rhob**3.3
     2         33333333333333d-1*t1*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t2*t6
     1         *wght/t5
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\right)-0.9305257363491\,{\it t_1}
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pbe_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhoa**2.3333333333333334d+0
            t9 = t3**2
            t10 = 1/t4**3
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t5*t6*wght
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t3*t7*t8-1.2
     1         407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5)*wght
     2         +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (3.5835038259110563d-3*gammaaa**2*t10*t
     1         2*t9/rhoa**6-4.0832233200718415d-2*gammaaa*t3*t7/rhoa**3.
     2         3333333333333337d+0-4.135669939329333d-1*t1*t5/rhoa**6.66
     3         6666666666666d-1)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0416116600359205d-2*t3*t7*t8-1.3438
     1         139347166464d-3*gammaaa*t10*t2*t9/rhoa**5)*wght+Cmat2(iq,
     2         D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t10*t2*t9*wght/r
     1         hoa**4+Cmat2(iq,D2_GAA_GAA)
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
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1.6455307846020564d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = rhob**1.3333333333333333d+0
            t10 = 1/t4**2
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/t7**2
            t13 = 1/rhob**2.3333333333333334d+0
            t14 = t3**2
            t15 = 1/t4**3
            t16 = 1/t7**3
            fnc(iq) = (-9.305257363491d-1*t1*t8*t9-9.305257363491d-1*t1*
     1         t5*t6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t10*t11*t3-1
     1         .2407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5)*wg
     2         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.0832233200718415d-2*gammabb*t12*t13*t3-1
     1         .2407009817987999d+0*rhob**3.333333333333333d-1*t1*t8)*wg
     2         ht+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t10*t
     1         3*wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t12*t
     1         3*wght/t9
            Amat2(iq,D2_RA_RA) = (-4.135669939329333d-1*t1*t5/rhoa**6.66
     1         6666666666666d-1-4.0832233200718415d-2*gammaaa*t10*t3/rho
     2         a**3.3333333333333337d+0+3.5835038259110563d-3*gammaaa**2
     3         *t14*t15*t2/rhoa**6)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-4.135669939329333d-1*t1*t8/rhob**6.66
     1         6666666666666d-1-4.0832233200718415d-2*gammabb*t12*t3/rho
     2         b**3.3333333333333337d+0+3.5835038259110563d-3*gammabb**2
     3         *t14*t16*t2/rhob**6)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0416116600359205d-2*t10*t11*t3-1.34
     1         38139347166464d-3*gammaaa*t14*t15*t2/rhoa**5)*wght+Cmat2(
     2         iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (2.0416116600359205d-2*t12*t13*t3-1.34
     1         38139347166464d-3*gammabb*t14*t16*t2/rhob**5)*wght+Cmat2(
     2         iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t14*t15*t2*wght/
     1         rhoa**4+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 5.039302255187424d-4*t14*t16*t2*wght/
     1         rhob**4+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhoa**2.3333333333333334d+0
            t9 = t3**2
            t10 = 1/t4**3
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t5*t6*wght
            Amat(iq,D1_RA) = 4.0832233200718415d-2*gammaaa*t3*t7*t8*wght
     1         -1.2407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Amat2(iq,D2_RA_RA) = 3.5835038259110563d-3*gammaaa**2*t10*t2
     1         *t9*wght/rhoa**6-4.0832233200718415d-2*gammaaa*t3*t7*wght
     2         /rhoa**3.3333333333333337d+0-4.135669939329333d-1*t1*t5*w
     3         ght/rhoa**6.666666666666666d-1+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -1.3438139347166464d-3*gammaaa*t10*t2*
     1         t9*wght/rhoa**5+2.0416116600359205d-2*t3*t7*t8*wght+Cmat2
     2         (iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t10*t2*t9*wght/r
     1         hoa**4+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhob**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhob**2.3333333333333334d+0
            t9 = t3**2
            t10 = 1/t4**3
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t5*t6*wght
            Amat(iq,D1_RB) = 4.0832233200718415d-2*gammabb*t3*t7*t8*wght
     1         -1.2407009817987999d+0*rhob**3.333333333333333d-1*t1*t5*w
     2         ght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Amat2(iq,D2_RB_RB) = 3.5835038259110563d-3*gammabb**2*t10*t2
     1         *t9*wght/rhob**6-4.0832233200718415d-2*gammabb*t3*t7*wght
     2         /rhob**3.3333333333333337d+0-4.135669939329333d-1*t1*t5*w
     3         ght/rhob**6.666666666666666d-1+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -1.3438139347166464d-3*gammabb*t10*t2*
     1         t9*wght/rhob**5+2.0416116600359205d-2*t3*t7*t8*wght+Cmat2
     2         (iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 5.039302255187424d-4*t10*t2*t9*wght/r
     1         hob**4+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0}}\right)-0.9305257363491\,{\it t_1}
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}
C>    \over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}
C>    \,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{
C>    \it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0}}\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77
C>    , 3865 (1996)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_pbe_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhoa**2.3333333333333334d+0
            t9 = t3**2
            t10 = gammaaa**2
            t11 = 1/t4**3
            t12 = 1/rhoa**6
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**5
            t15 = 1/t1**2
            t16 = t3**3
            t17 = 1/t4**4
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t5*t6*wght
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t3*t7*t8-1.2
     1         407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5)*wght
     2         +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (3.5835038259110563d-3*t10*t11*t12*t2*t
     1         9-4.0832233200718415d-2*gammaaa*t13*t3*t7-4.1356699393293
     2         33d-1*t1*t5/rhoa**6.666666666666666d-1)*wght+Amat2(iq,D2_
     3         RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0416116600359205d-2*t3*t7*t8-1.3438
     1         139347166464d-3*gammaaa*t11*t14*t2*t9)*wght+Cmat2(iq,D2_R
     2         A_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t11*t2*t9*wght/r
     1         hoa**4+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-2.50845267813774d-2*t10*t11*t2*t9/
     1         rhoa**7+1.5425510320271402d-1*gammaaa*t3*t7/rhoa**4.33333
     2         3333333333d+0+2.757113292886222d-1*t1*t5/rhoa**1.66666666
     3         66666669d+0+4.717412689820714d-4*gammaaa**3*t15*t16*t17/r
     4         hoa**9.666666666666666d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (8.51082158653876d-3*gammaaa*t11*t1
     1         2*t2*t9-4.7637605400838157d-2*t13*t3*t7-1.769029758682767
     2         7d-4*t10*t15*t16*t17/rhoa**8.666666666666666d+0)*wght+Cma
     3         t3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (6.633861595060381d-5*gammaaa*t15*
     1         t16*t17/rhoa**7.666666666666667d+0-2.0157209020749695d-3*
     2         t11*t14*t2*t9)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-2.487698
     1         0981476426d-5*t15*t16*t17*wght/rhoa**6.666666666666667d+0
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
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1.6455307846020564d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = rhob**1.3333333333333333d+0
            t10 = 1/t4**2
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/t7**2
            t13 = 1/rhob**2.3333333333333334d+0
            t14 = t3**2
            t15 = gammaaa**2
            t16 = 1/t4**3
            t17 = 1/rhoa**6
            t18 = 1/rhoa**3.3333333333333337d+0
            t19 = gammabb**2
            t20 = 1/t7**3
            t21 = 1/rhob**6
            t22 = 1/rhob**3.3333333333333337d+0
            t23 = 1/rhoa**5
            t24 = 1/rhob**5
            t25 = 1/t1**2
            t26 = t3**3
            t27 = 1/t4**4
            t28 = 1/t7**4
            fnc(iq) = (-9.305257363491d-1*t1*t8*t9-9.305257363491d-1*t1*
     1         t5*t6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t10*t11*t3-1
     1         .2407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5)*wg
     2         ht+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.0832233200718415d-2*gammabb*t12*t13*t3-1
     1         .2407009817987999d+0*rhob**3.333333333333333d-1*t1*t8)*wg
     2         ht+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t10*t
     1         3*wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t12*t
     1         3*wght/t9
            Amat2(iq,D2_RA_RA) = (-4.135669939329333d-1*t1*t5/rhoa**6.66
     1         6666666666666d-1-4.0832233200718415d-2*gammaaa*t10*t18*t3
     2         +3.5835038259110563d-3*t14*t15*t16*t17*t2)*wght+Amat2(iq,
     3         D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-4.135669939329333d-1*t1*t8/rhob**6.66
     1         6666666666666d-1-4.0832233200718415d-2*gammabb*t12*t22*t3
     2         +3.5835038259110563d-3*t14*t19*t2*t20*t21)*wght+Amat2(iq,
     3         D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (2.0416116600359205d-2*t10*t11*t3-1.34
     1         38139347166464d-3*gammaaa*t14*t16*t2*t23)*wght+Cmat2(iq,D
     2         2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (2.0416116600359205d-2*t12*t13*t3-1.34
     1         38139347166464d-3*gammabb*t14*t2*t20*t24)*wght+Cmat2(iq,D
     2         2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t14*t16*t2*wght/
     1         rhoa**4+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 5.039302255187424d-4*t14*t2*t20*wght/
     1         rhob**4+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (2.757113292886222d-1*t1*t5/rhoa**1.
     1         6666666666666669d+0+1.5425510320271402d-1*gammaaa*t10*t3/
     2         rhoa**4.333333333333333d+0+4.717412689820714d-4*gammaaa**
     3         3*t25*t26*t27/rhoa**9.666666666666666d+0-2.50845267813774
     4         d-2*t14*t15*t16*t2/rhoa**7)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (2.757113292886222d-1*t1*t8/rhob**1.
     1         6666666666666669d+0+1.5425510320271402d-1*gammabb*t12*t3/
     2         rhob**4.333333333333333d+0+4.717412689820714d-4*gammabb**
     3         3*t25*t26*t28/rhob**9.666666666666666d+0-2.50845267813774
     4         d-2*t14*t19*t2*t20/rhob**7)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-4.7637605400838157d-2*t10*t18*t3-
     1         1.7690297586827677d-4*t15*t25*t26*t27/rhoa**8.66666666666
     2         6666d+0+8.51082158653876d-3*gammaaa*t14*t16*t17*t2)*wght+
     3         Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-4.7637605400838157d-2*t12*t22*t3-
     1         1.7690297586827677d-4*t19*t25*t26*t28/rhob**8.66666666666
     2         6666d+0+8.51082158653876d-3*gammabb*t14*t2*t20*t21)*wght+
     3         Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (6.633861595060381d-5*gammaaa*t25*
     1         t26*t27/rhoa**7.666666666666667d+0-2.0157209020749695d-3*
     2         t14*t16*t2*t23)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (6.633861595060381d-5*gammabb*t25*
     1         t26*t28/rhob**7.666666666666667d+0-2.0157209020749695d-3*
     2         t14*t2*t20*t24)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-2.487698
     1         0981476426d-5*t25*t26*t27*wght/rhoa**6.666666666666667d+0
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-2.487698
     1         0981476426d-5*t25*t26*t28*wght/rhob**6.666666666666667d+0
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhoa**2.3333333333333334d+0
            t9 = t3**2
            t10 = gammaaa**2
            t11 = 1/t4**3
            t12 = 1/rhoa**6
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**5
            t15 = 1/t1**2
            t16 = t3**3
            t17 = 1/t4**4
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t5*t6*wght
            Amat(iq,D1_RA) = 4.0832233200718415d-2*gammaaa*t3*t7*t8*wght
     1         -1.2407009817987999d+0*rhoa**3.333333333333333d-1*t1*t5*w
     2         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Amat2(iq,D2_RA_RA) = 3.5835038259110563d-3*t10*t11*t12*t2*t9
     1         *wght-4.0832233200718415d-2*gammaaa*t13*t3*t7*wght-4.1356
     2         69939329333d-1*t1*t5*wght/rhoa**6.666666666666666d-1+Amat
     3         2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -1.3438139347166464d-3*gammaaa*t11*t14
     1         *t2*t9*wght+2.0416116600359205d-2*t3*t7*t8*wght+Cmat2(iq,
     2         D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 5.039302255187424d-4*t11*t2*t9*wght/r
     1         hoa**4+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = -2.50845267813774d-2*t10*t11*t2*t9*w
     1         ght/rhoa**7+1.5425510320271402d-1*gammaaa*t3*t7*wght/rhoa
     2         **4.333333333333333d+0+2.757113292886222d-1*t1*t5*wght/rh
     3         oa**1.6666666666666669d+0+4.717412689820714d-4*gammaaa**3
     4         *t15*t16*t17*wght/rhoa**9.666666666666666d+0+Amat3(iq,D3_
     5         RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = 8.51082158653876d-3*gammaaa*t11*t12
     1         *t2*t9*wght-4.7637605400838157d-2*t13*t3*t7*wght-1.769029
     2         7586827677d-4*t10*t15*t16*t17*wght/rhoa**8.66666666666666
     3         6d+0+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -2.0157209020749695d-3*t11*t14*t2*
     1         t9*wght+6.633861595060381d-5*gammaaa*t15*t16*t17*wght/rho
     2         a**7.666666666666667d+0+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-2.487698
     1         0981476426d-5*t15*t16*t17*wght/rhoa**6.666666666666667d+0
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020564d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhob**1.3333333333333333d+0
            t7 = 1/t4**2
            t8 = 1/rhob**2.3333333333333334d+0
            t9 = t3**2
            t10 = gammabb**2
            t11 = 1/t4**3
            t12 = 1/rhob**6
            t13 = 1/rhob**3.3333333333333337d+0
            t14 = 1/rhob**5
            t15 = 1/t1**2
            t16 = t3**3
            t17 = 1/t4**4
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t5*t6*wght
            Amat(iq,D1_RB) = 4.0832233200718415d-2*gammabb*t3*t7*t8*wght
     1         -1.2407009817987999d+0*rhob**3.333333333333333d-1*t1*t5*w
     2         ght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.531208745026941d-2*t3*t7
     1         *wght/t6
            Amat2(iq,D2_RB_RB) = 3.5835038259110563d-3*t10*t11*t12*t2*t9
     1         *wght-4.0832233200718415d-2*gammabb*t13*t3*t7*wght-4.1356
     2         69939329333d-1*t1*t5*wght/rhob**6.666666666666666d-1+Amat
     3         2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -1.3438139347166464d-3*gammabb*t11*t14
     1         *t2*t9*wght+2.0416116600359205d-2*t3*t7*t8*wght+Cmat2(iq,
     2         D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 5.039302255187424d-4*t11*t2*t9*wght/r
     1         hob**4+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = -2.50845267813774d-2*t10*t11*t2*t9*w
     1         ght/rhob**7+1.5425510320271402d-1*gammabb*t3*t7*wght/rhob
     2         **4.333333333333333d+0+2.757113292886222d-1*t1*t5*wght/rh
     3         ob**1.6666666666666669d+0+4.717412689820714d-4*gammabb**3
     4         *t15*t16*t17*wght/rhob**9.666666666666666d+0+Amat3(iq,D3_
     5         RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = 8.51082158653876d-3*gammabb*t11*t12
     1         *t2*t9*wght-4.7637605400838157d-2*t13*t3*t7*wght-1.769029
     2         7586827677d-4*t10*t15*t16*t17*wght/rhob**8.66666666666666
     3         6d+0+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -2.0157209020749695d-3*t11*t14*t2*
     1         t9*wght+6.633861595060381d-5*gammabb*t15*t16*t17*wght/rho
     2         b**7.666666666666667d+0+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-2.487698
     1         0981476426d-5*t15*t16*t17*wght/rhob**6.666666666666667d+0
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
