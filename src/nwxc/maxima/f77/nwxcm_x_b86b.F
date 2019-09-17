C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_b86b.F 
C> The nwxcm_x_b86b functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_b86b functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}
C>    \,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="https://doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b86b(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 1/t1**1.8d+0
            fnc(iq) = fnc(iq)-7.5d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t2/rhoa**2.3
     1         333333333333334d+0-5.599999999999999d-5*gammaaa**2*t4/rho
     2         a**5)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t4/rhoa**4-
     1         3.75d-3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t5 = 1/t4**8.0d-1
            t6 = 1/rhob**1.3333333333333333d+0
            t7 = 1/t1**1.8d+0
            t8 = 1/t4**1.8d+0
            fnc(iq) = (-3.75d-3*gammabb*t5*t6-3.75d-3*gammaaa*t2*t3)*wgh
     1         t+fnc(iq)
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t2/rhoa**2.3
     1         333333333333334d+0-5.599999999999999d-5*gammaaa**2*t7/rho
     2         a**5)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.9999999999999994d-3*gammabb*t5/rhob**2.3
     1         333333333333334d+0-5.599999999999999d-5*gammabb**2*t8/rho
     2         b**5)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t7/rhoa**4-
     1         3.75d-3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (2.1000000000000002d-5*gammabb*t8/rhob**4-
     1         3.75d-3*t5*t6)*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 1/t1**1.8d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = -5.599999999999999d-5*gammaaa**2*t4*wght/rh
     1         oa**5+4.9999999999999994d-3*gammaaa*t2*wght/rhoa**2.33333
     2         33333333334d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 2.1000000000000002d-5*gammaaa*t4*wght/rhoa
     1         **4-3.75d-3*t2*t3*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhob**1.3333333333333333d+0
            t4 = 1/t1**1.8d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammabb*t2*t3*wght
            Amat(iq,D1_RB) = -5.599999999999999d-5*gammabb**2*t4*wght/rh
     1         ob**5+4.9999999999999994d-3*gammabb*t2*wght/rhob**2.33333
     2         33333333334d+0+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 2.1000000000000002d-5*gammabb*t4*wght/rhob
     1         **4-3.75d-3*t2*t3*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b86b functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}
C>    \,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="https://doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b86b_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = gammaaa**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhoa**5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhoa**4
            t9 = 1/t1**2.7999999999999997d+0
            fnc(iq) = fnc(iq)-7.5d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t2*t7-5.5999
     1         99999999999d-5*t4*t5*t6)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t5*t8-3.75d
     1         -3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-1.8816d-6*gammaaa**3*t9/rhoa**8.66666
     1         6666666666d+0+3.5466666666666663d-4*t4*t5/rhoa**6-1.16666
     2         66666666665d-2*gammaaa*t2/rhoa**3.3333333333333337d+0)*wg
     3         ht+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (7.056d-7*t4*t9/rhoa**7.66666666666666
     1         7d+0+4.9999999999999994d-3*t2*t7-1.3999999999999999d-4*ga
     2         mmaaa*t5*t6)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (4.2000000000000004d-5*t5*t8-2.646d-7
     1         *gammaaa*t9/rhoa**6.666666666666667d+0)*wght+Cmat2(iq,D2_
     2         GAA_GAA)
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t5 = 1/t4**8.0d-1
            t6 = 1/rhob**1.3333333333333333d+0
            t7 = gammaaa**2
            t8 = 1/t1**1.8d+0
            t9 = 1/rhoa**5
            t10 = 1/rhoa**2.3333333333333334d+0
            t11 = gammabb**2
            t12 = 1/t4**1.8d+0
            t13 = 1/rhob**5
            t14 = 1/rhob**2.3333333333333334d+0
            t15 = 1/rhoa**4
            t16 = 1/rhob**4
            t17 = 1/t1**2.7999999999999997d+0
            t18 = 1/t4**2.7999999999999997d+0
            fnc(iq) = (-3.75d-3*gammabb*t5*t6-3.75d-3*gammaaa*t2*t3)*wgh
     1         t+fnc(iq)
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t10*t2-5.599
     1         999999999999d-5*t7*t8*t9)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.9999999999999994d-3*gammabb*t14*t5-5.599
     1         999999999999d-5*t11*t12*t13)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t15*t8-3.75
     1         d-3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (2.1000000000000002d-5*gammabb*t12*t16-3.7
     1         5d-3*t5*t6)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (3.5466666666666663d-4*t7*t8/rhoa**6-1.
     1         1666666666666665d-2*gammaaa*t2/rhoa**3.3333333333333337d+
     2         0-1.8816d-6*gammaaa**3*t17/rhoa**8.666666666666666d+0)*wg
     3         ht+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-1.1666666666666665d-2*gammabb*t5/rhob
     1         **3.3333333333333337d+0-1.8816d-6*gammabb**3*t18/rhob**8.
     2         666666666666666d+0+3.5466666666666663d-4*t11*t12/rhob**6)
     3         *wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.3999999999999999d-4*gammaaa*t8*t9+
     1         7.056d-7*t17*t7/rhoa**7.666666666666667d+0+4.999999999999
     2         9994d-3*t10*t2)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.9999999999999994d-3*t14*t5+7.056d-7
     1         *t11*t18/rhob**7.666666666666667d+0-1.3999999999999999d-4
     2         *gammabb*t12*t13)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (4.2000000000000004d-5*t15*t8-2.646d-
     1         7*gammaaa*t17/rhoa**6.666666666666667d+0)*wght+Cmat2(iq,D
     2         2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (4.2000000000000004d-5*t12*t16-2.646d
     1         -7*gammabb*t18/rhob**6.666666666666667d+0)*wght+Cmat2(iq,
     2         D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = gammaaa**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhoa**5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhoa**4
            t9 = 1/t1**2.7999999999999997d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = 4.9999999999999994d-3*gammaaa*t2*t7*wght-5.
     1         599999999999999d-5*t4*t5*t6*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 2.1000000000000002d-5*gammaaa*t5*t8*wght-3
     1         .75d-3*t2*t3*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = -1.8816d-6*gammaaa**3*t9*wght/rhoa**8.6
     1         66666666666666d+0+3.5466666666666663d-4*t4*t5*wght/rhoa**
     2         6-1.1666666666666665d-2*gammaaa*t2*wght/rhoa**3.333333333
     3         3333337d+0+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 7.056d-7*t4*t9*wght/rhoa**7.6666666666
     1         66667d+0+4.9999999999999994d-3*t2*t7*wght-1.3999999999999
     2         999d-4*gammaaa*t5*t6*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = -2.646d-7*gammaaa*t9*wght/rhoa**6.666
     1         666666666667d+0+4.2000000000000004d-5*t5*t8*wght+Cmat2(iq
     2         ,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhob**1.3333333333333333d+0
            t4 = gammabb**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhob**5
            t7 = 1/rhob**2.3333333333333334d+0
            t8 = 1/rhob**4
            t9 = 1/t1**2.7999999999999997d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammabb*t2*t3*wght
            Amat(iq,D1_RB) = 4.9999999999999994d-3*gammabb*t2*t7*wght-5.
     1         599999999999999d-5*t4*t5*t6*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 2.1000000000000002d-5*gammabb*t5*t8*wght-3
     1         .75d-3*t2*t3*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = -1.8816d-6*gammabb**3*t9*wght/rhob**8.6
     1         66666666666666d+0+3.5466666666666663d-4*t4*t5*wght/rhob**
     2         6-1.1666666666666665d-2*gammabb*t2*wght/rhob**3.333333333
     3         3333337d+0+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 7.056d-7*t4*t9*wght/rhob**7.6666666666
     1         66667d+0+4.9999999999999994d-3*t2*t7*wght-1.3999999999999
     2         999d-4*gammabb*t5*t6*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = -2.646d-7*gammabb*t9*wght/rhob**6.666
     1         666666666667d+0+4.2000000000000004d-5*t5*t8*wght+Cmat2(iq
     2         ,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b86b functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,
C>    \sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}
C>    \,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}
C>    +1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="https://doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_b86b_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = gammaaa**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhoa**5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhoa**4
            t9 = gammaaa**3
            t10 = 1/t1**2.7999999999999997d+0
            t11 = 1/rhoa**8.666666666666666d+0
            t12 = 1/rhoa**6
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**7.666666666666667d+0
            t15 = 1/rhoa**6.666666666666667d+0
            t16 = 1/t1**3.8d+0
            fnc(iq) = fnc(iq)-7.5d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t2*t7-5.5999
     1         99999999999d-5*t4*t5*t6)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t5*t8-3.75d
     1         -3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-1.8816d-6*t10*t11*t9+3.54666666666666
     1         63d-4*t12*t4*t5-1.1666666666666665d-2*gammaaa*t13*t2)*wgh
     2         t+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (4.9999999999999994d-3*t2*t7-1.3999999
     1         999999999d-4*gammaaa*t5*t6+7.056d-7*t10*t14*t4)*wght+Cmat
     2         2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = (4.2000000000000004d-5*t5*t8-2.646d-7
     1         *gammaaa*t10*t15)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (2.8224d-5*t10*t9/rhoa**9.6666666666
     1         66666d+0-2.302222222222222d-3*t4*t5/rhoa**7+3.88888888888
     2         8889d-2*gammaaa*t2/rhoa**4.333333333333333d+0-9.834495999
     3         999997d-8*gammaaa**4*t16/rhoa**1.2333333333333334d+1)*wgh
     4         t+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (3.687936d-8*t16*t9/rhoa**1.1333333
     1         333333334d+1+7.746666666666666d-4*gammaaa*t12*t5-1.011359
     2         9999999999d-5*t10*t11*t4-1.1666666666666665d-2*t13*t2)*wg
     3         ht+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-1.6799999999999998d-4*t5*t6-1.38
     1         2976d-8*t16*t4/rhoa**1.0333333333333333d+1+3.175199999999
     2         9997d-6*gammaaa*t10*t14)*wght+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (5.186160000000001d-9*gammaaa*t16
     1         /rhoa**9.333333333333333d+0-7.938000000000002d-7*t10*t15)
     2         *wght+Cmat3(iq,D3_GAA_GAA_GAA)
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t5 = 1/t4**8.0d-1
            t6 = 1/rhob**1.3333333333333333d+0
            t7 = gammaaa**2
            t8 = 1/t1**1.8d+0
            t9 = 1/rhoa**5
            t10 = 1/rhoa**2.3333333333333334d+0
            t11 = gammabb**2
            t12 = 1/t4**1.8d+0
            t13 = 1/rhob**5
            t14 = 1/rhob**2.3333333333333334d+0
            t15 = 1/rhoa**4
            t16 = 1/rhob**4
            t17 = gammaaa**3
            t18 = 1/t1**2.7999999999999997d+0
            t19 = 1/rhoa**8.666666666666666d+0
            t20 = 1/rhoa**6
            t21 = 1/rhoa**3.3333333333333337d+0
            t22 = gammabb**3
            t23 = 1/t4**2.7999999999999997d+0
            t24 = 1/rhob**8.666666666666666d+0
            t25 = 1/rhob**6
            t26 = 1/rhob**3.3333333333333337d+0
            t27 = 1/rhoa**7.666666666666667d+0
            t28 = 1/rhob**7.666666666666667d+0
            t29 = 1/rhoa**6.666666666666667d+0
            t30 = 1/rhob**6.666666666666667d+0
            t31 = 1/t1**3.8d+0
            t32 = 1/t4**3.8d+0
            fnc(iq) = (-3.75d-3*gammabb*t5*t6-3.75d-3*gammaaa*t2*t3)*wgh
     1         t+fnc(iq)
            Amat(iq,D1_RA) = (4.9999999999999994d-3*gammaaa*t10*t2-5.599
     1         999999999999d-5*t7*t8*t9)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.9999999999999994d-3*gammabb*t14*t5-5.599
     1         999999999999d-5*t11*t12*t13)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (2.1000000000000002d-5*gammaaa*t15*t8-3.75
     1         d-3*t2*t3)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (2.1000000000000002d-5*gammabb*t12*t16-3.7
     1         5d-3*t5*t6)*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RA_RA) = (3.5466666666666663d-4*t20*t7*t8-1.1666
     1         666666666665d-2*gammaaa*t2*t21-1.8816d-6*t17*t18*t19)*wgh
     2         t+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-1.1666666666666665d-2*gammabb*t26*t5+
     1         3.5466666666666663d-4*t11*t12*t25-1.8816d-6*t22*t23*t24)*
     2         wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (-1.3999999999999999d-4*gammaaa*t8*t9+
     1         7.056d-7*t18*t27*t7+4.9999999999999994d-3*t10*t2)*wght+Cm
     2         at2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.9999999999999994d-3*t14*t5+7.056d-7
     1         *t11*t23*t28-1.3999999999999999d-4*gammabb*t12*t13)*wght+
     2         Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = (4.2000000000000004d-5*t15*t8-2.646d-
     1         7*gammaaa*t18*t29)*wght+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = (4.2000000000000004d-5*t12*t16-2.646d
     1         -7*gammabb*t23*t30)*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (-2.302222222222222d-3*t7*t8/rhoa**7
     1         -9.834495999999997d-8*gammaaa**4*t31/rhoa**1.233333333333
     2         3334d+1+3.888888888888889d-2*gammaaa*t2/rhoa**4.333333333
     3         333333d+0+2.8224d-5*t17*t18/rhoa**9.666666666666666d+0)*w
     4         ght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (3.888888888888889d-2*gammabb*t5/rho
     1         b**4.333333333333333d+0-9.834495999999997d-8*gammabb**4*t
     2         32/rhob**1.2333333333333334d+1+2.8224d-5*t22*t23/rhob**9.
     3         666666666666666d+0-2.302222222222222d-3*t11*t12/rhob**7)*
     4         wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (7.746666666666666d-4*gammaaa*t20*t
     1         8-1.0113599999999999d-5*t18*t19*t7+3.687936d-8*t17*t31/rh
     2         oa**1.1333333333333334d+1-1.1666666666666665d-2*t2*t21)*w
     3         ght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-1.1666666666666665d-2*t26*t5+3.68
     1         7936d-8*t22*t32/rhob**1.1333333333333334d+1+7.74666666666
     2         6666d-4*gammabb*t12*t25-1.0113599999999999d-5*t11*t23*t24
     3         )*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RA_GAA_GAA) = (-1.6799999999999998d-4*t8*t9-1.38
     1         2976d-8*t31*t7/rhoa**1.0333333333333333d+1+3.175199999999
     2         9997d-6*gammaaa*t18*t27)*wght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (-1.382976d-8*t11*t32/rhob**1.0333
     1         333333333333d+1+3.1751999999999997d-6*gammabb*t23*t28-1.6
     2         799999999999998d-4*t12*t13)*wght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = (5.186160000000001d-9*gammaaa*t31
     1         /rhoa**9.333333333333333d+0-7.938000000000002d-7*t18*t29)
     2         *wght+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = (5.186160000000001d-9*gammabb*t32
     1         /rhob**9.333333333333333d+0-7.938000000000002d-7*t23*t30)
     2         *wght+Cmat3(iq,D3_GBB_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = gammaaa**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhoa**5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhoa**4
            t9 = gammaaa**3
            t10 = 1/t1**2.7999999999999997d+0
            t11 = 1/rhoa**8.666666666666666d+0
            t12 = 1/rhoa**6
            t13 = 1/rhoa**3.3333333333333337d+0
            t14 = 1/rhoa**7.666666666666667d+0
            t15 = 1/rhoa**6.666666666666667d+0
            t16 = 1/t1**3.8d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = 4.9999999999999994d-3*gammaaa*t2*t7*wght-5.
     1         599999999999999d-5*t4*t5*t6*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 2.1000000000000002d-5*gammaaa*t5*t8*wght-3
     1         .75d-3*t2*t3*wght+Cmat(iq,D1_GAA)
            Amat2(iq,D2_RA_RA) = -1.8816d-6*t10*t11*t9*wght+3.5466666666
     1         666663d-4*t12*t4*t5*wght-1.1666666666666665d-2*gammaaa*t1
     2         3*t2*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 4.9999999999999994d-3*t2*t7*wght-1.399
     1         9999999999999d-4*gammaaa*t5*t6*wght+7.056d-7*t10*t14*t4*w
     2         ght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 4.2000000000000004d-5*t5*t8*wght-2.64
     1         6d-7*gammaaa*t10*t15*wght+Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = 2.8224d-5*t10*t9*wght/rhoa**9.666666
     1         666666666d+0-2.302222222222222d-3*t4*t5*wght/rhoa**7+3.88
     2         8888888888889d-2*gammaaa*t2*wght/rhoa**4.333333333333333d
     3         +0-9.834495999999997d-8*gammaaa**4*t16*wght/rhoa**1.23333
     4         33333333334d+1+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = 3.687936d-8*t16*t9*wght/rhoa**1.133
     1         3333333333334d+1+7.746666666666666d-4*gammaaa*t12*t5*wght
     2         -1.0113599999999999d-5*t10*t11*t4*wght-1.1666666666666665
     3         d-2*t13*t2*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -1.6799999999999998d-4*t5*t6*wght-
     1         1.382976d-8*t16*t4*wght/rhoa**1.0333333333333333d+1+3.175
     2         1999999999997d-6*gammaaa*t10*t14*wght+Cmat3(iq,D3_RA_GAA_
     3         GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = 5.186160000000001d-9*gammaaa*t16*
     1         wght/rhoa**9.333333333333333d+0-7.938000000000002d-7*t10*
     2         t15*wght+Cmat3(iq,D3_GAA_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 7.000000000000001d-3*gammabb/rhob**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhob**1.3333333333333333d+0
            t4 = gammabb**2
            t5 = 1/t1**1.8d+0
            t6 = 1/rhob**5
            t7 = 1/rhob**2.3333333333333334d+0
            t8 = 1/rhob**4
            t9 = gammabb**3
            t10 = 1/t1**2.7999999999999997d+0
            t11 = 1/rhob**8.666666666666666d+0
            t12 = 1/rhob**6
            t13 = 1/rhob**3.3333333333333337d+0
            t14 = 1/rhob**7.666666666666667d+0
            t15 = 1/rhob**6.666666666666667d+0
            t16 = 1/t1**3.8d+0
            fnc(iq) = fnc(iq)-3.75d-3*gammabb*t2*t3*wght
            Amat(iq,D1_RB) = 4.9999999999999994d-3*gammabb*t2*t7*wght-5.
     1         599999999999999d-5*t4*t5*t6*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 2.1000000000000002d-5*gammabb*t5*t8*wght-3
     1         .75d-3*t2*t3*wght+Cmat(iq,D1_GBB)
            Amat2(iq,D2_RB_RB) = -1.8816d-6*t10*t11*t9*wght+3.5466666666
     1         666663d-4*t12*t4*t5*wght-1.1666666666666665d-2*gammabb*t1
     2         3*t2*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 4.9999999999999994d-3*t2*t7*wght-1.399
     1         9999999999999d-4*gammabb*t5*t6*wght+7.056d-7*t10*t14*t4*w
     2         ght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 4.2000000000000004d-5*t5*t8*wght-2.64
     1         6d-7*gammabb*t10*t15*wght+Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = 2.8224d-5*t10*t9*wght/rhob**9.666666
     1         666666666d+0-2.302222222222222d-3*t4*t5*wght/rhob**7+3.88
     2         8888888888889d-2*gammabb*t2*wght/rhob**4.333333333333333d
     3         +0-9.834495999999997d-8*gammabb**4*t16*wght/rhob**1.23333
     4         33333333334d+1+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = 3.687936d-8*t16*t9*wght/rhob**1.133
     1         3333333333334d+1+7.746666666666666d-4*gammabb*t12*t5*wght
     2         -1.0113599999999999d-5*t10*t11*t4*wght-1.1666666666666665
     3         d-2*t13*t2*wght+Cmat3(iq,D3_RB_RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -1.6799999999999998d-4*t5*t6*wght-
     1         1.382976d-8*t16*t4*wght/rhob**1.0333333333333333d+1+3.175
     2         1999999999997d-6*gammabb*t10*t14*wght+Cmat3(iq,D3_RB_GBB_
     3         GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = 5.186160000000001d-9*gammabb*t16*
     1         wght/rhob**9.333333333333333d+0-7.938000000000002d-7*t10*
     2         t15*wght+Cmat3(iq,D3_GBB_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
