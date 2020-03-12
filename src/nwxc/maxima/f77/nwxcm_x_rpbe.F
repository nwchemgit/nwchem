C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_rpbe.F 
C> The nwxcm_x_rpbe functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_rpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0
C>    -e^ {- {{0.00449276922095889\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{ss}}\over{\rho_s^{{{8}
C>    \over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59
C>    , 7413 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_rpbe(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2/rhoa**2.33
     1         33333333333334d+0-9.975235893662352d-1*rhoa**3.3333333333
     2         33333d-1*t3)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = rhob**1.3333333333333333d+0
            t5 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t6 = 1.0d+0-t5
            fnc(iq) = (-7.481426920246764d-1*t4*t6-7.481426920246764d-1*
     1         t1*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2/rhoa**2.33
     1         33333333333334d+0-9.975235893662352d-1*rhoa**3.3333333333
     2         33333d-1*t3)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (8.963286558970113d-3*gammabb*t5/rhob**2.33
     1         33333333333334d+0-9.975235893662352d-1*rhob**3.3333333333
     2         33333d-1*t6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t5*wg
     1         ht/t4
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = -9.975235893662352d-1*rhoa**3.3333333333333
     1         33d-1*t3*wght+8.963286558970113d-3*gammaaa*t2*wght/rhoa**
     2         2.3333333333333334d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = -9.975235893662352d-1*rhob**3.3333333333333
     1         33d-1*t3*wght+8.963286558970113d-3*gammabb*t2*wght/rhob**
     2         2.3333333333333334d+0+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t2*wg
     1         ht/t1
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_rpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0
C>    -e^ {- {{0.00449276922095889\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{ss}}\over{\rho_s^{{{8}
C>    \over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59
C>    , 7413 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_rpbe_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2*t4-9.97523
     1         5893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-3.3250786312207836d-1*t3/rhoa**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammaaa*t2/rhoa**3.33
     2         33333333333337d+0+1.0738660792206783d-4*gammaaa**2*t2/rho
     3         a**6)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (4.4816432794850564d-3*t2*t4-4.0269977
     1         97077544d-5*gammaaa*t2/rhoa**5)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = rhob**1.3333333333333333d+0
            t5 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t6 = 1.0d+0-t5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = (-7.481426920246764d-1*t4*t6-7.481426920246764d-1*
     1         t1*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2*t7-9.97523
     1         5893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq
     2         ,D1_RA)
            Amat(iq,D1_RB) = (8.963286558970113d-3*gammabb*t5*t8-9.97523
     1         5893662352d-1*rhob**3.333333333333333d-1*t6)*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t5*wg
     1         ht/t4
            Amat2(iq,D2_RA_RA) = (-3.3250786312207836d-1*t3/rhoa**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammaaa*t2/rhoa**3.33
     2         33333333333337d+0+1.0738660792206783d-4*gammaaa**2*t2/rho
     3         a**6)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-3.3250786312207836d-1*t6/rhob**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammabb*t5/rhob**3.33
     2         33333333333337d+0+1.0738660792206783d-4*gammabb**2*t5/rho
     3         b**6)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (4.4816432794850564d-3*t2*t7-4.0269977
     1         97077544d-5*gammaaa*t2/rhoa**5)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.4816432794850564d-3*t5*t8-4.0269977
     1         97077544d-5*gammabb*t5/rhob**5)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.510124173904079d-5*t5*wght/rhob**4+
     1         Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = 8.963286558970113d-3*gammaaa*t2*t4*wght-9.9
     1         75235893662352d-1*rhoa**3.333333333333333d-1*t3*wght+Amat
     2         (iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Amat2(iq,D2_RA_RA) = -3.3250786312207836d-1*t3*wght/rhoa**6.
     1         666666666666666d-1-8.963286558970114d-3*gammaaa*t2*wght/r
     2         hoa**3.3333333333333337d+0+1.0738660792206783d-4*gammaaa*
     3         *2*t2*wght/rhoa**6+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 4.4816432794850564d-3*t2*t4*wght-4.026
     1         997797077544d-5*gammaaa*t2*wght/rhoa**5+Cmat2(iq,D2_RA_GA
     2         A)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = 8.963286558970113d-3*gammabb*t2*t4*wght-9.9
     1         75235893662352d-1*rhob**3.333333333333333d-1*t3*wght+Amat
     2         (iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Amat2(iq,D2_RB_RB) = -3.3250786312207836d-1*t3*wght/rhob**6.
     1         666666666666666d-1-8.963286558970114d-3*gammabb*t2*wght/r
     2         hob**3.3333333333333337d+0+1.0738660792206783d-4*gammabb*
     3         *2*t2*wght/rhob**6+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 4.4816432794850564d-3*t2*t4*wght-4.026
     1         997797077544d-5*gammabb*t2*wght/rhob**5+Cmat2(iq,D2_RB_GB
     2         B)
            Cmat2(iq,D2_GBB_GBB) = 1.510124173904079d-5*t2*wght/rhob**4+
     1         Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_rpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0
C>    -e^ {- {{0.00449276922095889\,\sigma_{\beta\beta}}
C>    \over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764
C>    \,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{\alpha\alpha}}
C>    \over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {
C>    - {{0.00449276922095889\,\sigma_{ss}}\over{\rho_s^{{{8}
C>    \over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59
C>    , 7413 (1999)  , DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_rpbe_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = gammaaa**2
            t6 = 1/rhoa**6
            t7 = 1/rhoa**3.3333333333333337d+0
            t8 = 1/rhoa**5
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2*t4-9.97523
     1         5893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq
     2         ,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-8.963286558970114d-3*gammaaa*t2*t7+1.
     1         0738660792206783d-4*t2*t5*t6-3.3250786312207836d-1*t3/rho
     2         a**6.666666666666666d-1)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = (4.4816432794850564d-3*t2*t4-4.0269977
     1         97077544d-5*gammaaa*t2*t8)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (-7.517062554544749d-4*t2*t5/rhoa**7
     1         +2.2167190874805223d-1*t3/rhoa**1.6666666666666669d+0+3.3
     2         86130477833154d-2*gammaaa*t2/rhoa**4.333333333333333d+0+1
     3         .286568658174524d-6*gammaaa**3*t2/rhoa**9.666666666666666
     4         d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-1.0457167652131799d-2*t2*t7+2.550
     1         431938149111d-4*gammaaa*t2*t6-4.8246324681544644d-7*t2*t5
     2         /rhoa**8.666666666666666d+0)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = (1.8092371755579242d-7*gammaaa*t2/
     1         rhoa**7.666666666666667d+0-6.040496695616315d-5*t2*t8)*wg
     2         ht+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-6.784639
     1         408342217d-8*t2*wght/rhoa**6.666666666666667d+0
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = rhob**1.3333333333333333d+0
            t5 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t6 = 1.0d+0-t5
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhob**2.3333333333333334d+0
            t9 = gammaaa**2
            t10 = 1/rhoa**6
            t11 = 1/rhoa**3.3333333333333337d+0
            t12 = gammabb**2
            t13 = 1/rhob**6
            t14 = 1/rhob**3.3333333333333337d+0
            t15 = 1/rhoa**5
            t16 = 1/rhob**5
            fnc(iq) = (-7.481426920246764d-1*t4*t6-7.481426920246764d-1*
     1         t1*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (8.963286558970113d-3*gammaaa*t2*t7-9.97523
     1         5893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq
     2         ,D1_RA)
            Amat(iq,D1_RB) = (8.963286558970113d-3*gammabb*t5*t8-9.97523
     1         5893662352d-1*rhob**3.333333333333333d-1*t6)*wght+Amat(iq
     2         ,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t5*wg
     1         ht/t4
            Amat2(iq,D2_RA_RA) = (1.0738660792206783d-4*t10*t2*t9-3.3250
     1         786312207836d-1*t3/rhoa**6.666666666666666d-1-8.963286558
     2         970114d-3*gammaaa*t11*t2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-3.3250786312207836d-1*t6/rhob**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammabb*t14*t5+1.0738
     2         660792206783d-4*t12*t13*t5)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (4.4816432794850564d-3*t2*t7-4.0269977
     1         97077544d-5*gammaaa*t15*t2)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.4816432794850564d-3*t5*t8-4.0269977
     1         97077544d-5*gammabb*t16*t5)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.510124173904079d-5*t5*wght/rhob**4+
     1         Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (-7.517062554544749d-4*t2*t9/rhoa**7
     1         +2.2167190874805223d-1*t3/rhoa**1.6666666666666669d+0+3.3
     2         86130477833154d-2*gammaaa*t2/rhoa**4.333333333333333d+0+1
     3         .286568658174524d-6*gammaaa**3*t2/rhoa**9.666666666666666
     4         d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (2.2167190874805223d-1*t6/rhob**1.66
     1         66666666666669d+0-7.517062554544749d-4*t12*t5/rhob**7+3.3
     2         86130477833154d-2*gammabb*t5/rhob**4.333333333333333d+0+1
     3         .286568658174524d-6*gammabb**3*t5/rhob**9.666666666666666
     4         d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-4.8246324681544644d-7*t2*t9/rhoa*
     1         *8.666666666666666d+0-1.0457167652131799d-2*t11*t2+2.5504
     2         31938149111d-4*gammaaa*t10*t2)*wght+Cmat3(iq,D3_RA_RA_GAA
     3         )
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-1.0457167652131799d-2*t14*t5+2.55
     1         0431938149111d-4*gammabb*t13*t5-4.8246324681544644d-7*t12
     2         *t5/rhob**8.666666666666666d+0)*wght+Cmat3(iq,D3_RB_RB_GB
     3         B)
            Cmat3(iq,D3_RA_GAA_GAA) = (1.8092371755579242d-7*gammaaa*t2/
     1         rhoa**7.666666666666667d+0-6.040496695616315d-5*t15*t2)*w
     2         ght+Cmat3(iq,D3_RA_GAA_GAA)
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
            Cmat3(iq,D3_RB_GBB_GBB) = (1.8092371755579242d-7*gammabb*t5/
     1         rhob**7.666666666666667d+0-6.040496695616315d-5*t16*t5)*w
     2         ght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-6.784639
     1         408342217d-8*t2*wght/rhoa**6.666666666666667d+0
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-6.784639
     1         408342217d-8*t5*wght/rhob**6.666666666666667d+0
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammaaa/rhoa**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = gammaaa**2
            t6 = 1/rhoa**6
            t7 = 1/rhoa**3.3333333333333337d+0
            t8 = 1/rhoa**5
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = 8.963286558970113d-3*gammaaa*t2*t4*wght-9.9
     1         75235893662352d-1*rhoa**3.333333333333333d-1*t3*wght+Amat
     2         (iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Amat2(iq,D2_RA_RA) = -8.963286558970114d-3*gammaaa*t2*t7*wgh
     1         t+1.0738660792206783d-4*t2*t5*t6*wght-3.3250786312207836d
     2         -1*t3*wght/rhoa**6.666666666666666d-1+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -4.026997797077544d-5*gammaaa*t2*t8*wg
     1         ht+4.4816432794850564d-3*t2*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.510124173904079d-5*t2*wght/rhoa**4+
     1         Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = -7.517062554544749d-4*t2*t5*wght/rho
     1         a**7+2.2167190874805223d-1*t3*wght/rhoa**1.66666666666666
     2         69d+0+3.386130477833154d-2*gammaaa*t2*wght/rhoa**4.333333
     3         333333333d+0+1.286568658174524d-6*gammaaa**3*t2*wght/rhoa
     4         **9.666666666666666d+0+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -1.0457167652131799d-2*t2*t7*wght+2
     1         .550431938149111d-4*gammaaa*t2*t6*wght-4.8246324681544644
     2         d-7*t2*t5*wght/rhoa**8.666666666666666d+0+Cmat3(iq,D3_RA_
     3         RA_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -6.040496695616315d-5*t2*t8*wght+1
     1         .8092371755579242d-7*gammaaa*t2*wght/rhoa**7.666666666666
     2         667d+0+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-6.784639
     1         408342217d-8*t2*wght/rhoa**6.666666666666667d+0
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.49276922095889d-3*gammabb/rhob**2.66666666666666
     1         66d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhob**2.3333333333333334d+0
            t5 = gammabb**2
            t6 = 1/rhob**6
            t7 = 1/rhob**3.3333333333333337d+0
            t8 = 1/rhob**5
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = 8.963286558970113d-3*gammabb*t2*t4*wght-9.9
     1         75235893662352d-1*rhob**3.333333333333333d-1*t3*wght+Amat
     2         (iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.361232459613792d-3*t2*wg
     1         ht/t1
            Amat2(iq,D2_RB_RB) = -8.963286558970114d-3*gammabb*t2*t7*wgh
     1         t+1.0738660792206783d-4*t2*t5*t6*wght-3.3250786312207836d
     2         -1*t3*wght/rhob**6.666666666666666d-1+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -4.026997797077544d-5*gammabb*t2*t8*wg
     1         ht+4.4816432794850564d-3*t2*t4*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.510124173904079d-5*t2*wght/rhob**4+
     1         Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = -7.517062554544749d-4*t2*t5*wght/rho
     1         b**7+2.2167190874805223d-1*t3*wght/rhob**1.66666666666666
     2         69d+0+3.386130477833154d-2*gammabb*t2*wght/rhob**4.333333
     3         333333333d+0+1.286568658174524d-6*gammabb**3*t2*wght/rhob
     4         **9.666666666666666d+0+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -1.0457167652131799d-2*t2*t7*wght+2
     1         .550431938149111d-4*gammabb*t2*t6*wght-4.8246324681544644
     2         d-7*t2*t5*wght/rhob**8.666666666666666d+0+Cmat3(iq,D3_RB_
     3         RB_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -6.040496695616315d-5*t2*t8*wght+1
     1         .8092371755579242d-7*gammabb*t2*wght/rhob**7.666666666666
     2         667d+0+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-6.784639
     1         408342217d-8*t2*wght/rhob**6.666666666666667d+0
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
