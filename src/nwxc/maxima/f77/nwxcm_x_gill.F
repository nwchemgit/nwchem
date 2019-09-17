C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_gill.F 
C> The nwxcm_x_gill functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_gill functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}
C>    \over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}
C>    -{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}
C>    \over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}
C>    -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}
C>    \over{\rho_s^{{{2}\over{3}}}}}-0.9305257363491002
C>    \,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="https://doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_gill(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            fnc(iq) = (-1.45985401459854d-2*t1*t2-1.8610514726982003d+0*
     1         rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1/rhoa**1.666666666
     1         6666669d+0-1.2407009817988002d+0*rhoa**3.333333333333333d
     2         -1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*wg
     1         ht/gammaaa**2.5d-1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = gammabb**7.5d-1
            t4 = 1/rhob**6.666666666666666d-1
            fnc(iq) = (-7.2992700729927d-3*t3*t4-7.2992700729927d-3*t1*t
     1         2-9.305257363491002d-1*rhob**1.3333333333333333d+0-9.3052
     2         57363491002d-1*rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1/rhoa**1.666666666
     1         6666669d+0-1.2407009817988002d+0*rhoa**3.333333333333333d
     2         -1)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t3/rhob**1.666666666
     1         6666669d+0-1.2407009817988002d+0*rhob**3.333333333333333d
     2         -1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*wg
     1         ht/gammaaa**2.5d-1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t4*wg
     1         ht/gammabb**2.5d-1
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         oa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1/rhoa**1.666666666
     1         6666669d+0-1.2407009817988002d+0*rhoa**3.333333333333333d
     2         -1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*wg
     1         ht/gammaaa**2.5d-1
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**7.5d-1
            t2 = 1/rhob**6.666666666666666d-1
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         ob**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t1/rhob**1.666666666
     1         6666669d+0-1.2407009817988002d+0*rhob**3.333333333333333d
     2         -1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t2*wg
     1         ht/gammabb**2.5d-1
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_gill functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}
C>    \over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}
C>    -{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}
C>    \over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}
C>    -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}
C>    \over{\rho_s^{{{2}\over{3}}}}}-0.9305257363491002
C>    \,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="https://doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_gill_d2(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 1/gammaaa**2.5d-1
            fnc(iq) = (-1.45985401459854d-2*t1*t2-1.8610514726982003d+0*
     1         rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t4
     1         *wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-4.135669939329334d-1*t2-8.11030008110
     1         3d-3*t1/rhoa**2.6666666666666666d+0)*wght+Amat2(iq,D2_RA_
     2         RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t2*wght/gammaaa*
     1         *1.25d+0+Cmat2(iq,D2_GAA_GAA)
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = gammabb**7.5d-1
            t4 = 1/rhob**6.666666666666666d-1
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 1/rhob**1.6666666666666669d+0
            t7 = 1/gammaaa**2.5d-1
            t8 = 1/gammabb**2.5d-1
            fnc(iq) = (-7.2992700729927d-3*t3*t4-7.2992700729927d-3*t1*t
     1         2-9.305257363491002d-1*rhob**1.3333333333333333d+0-9.3052
     2         57363491002d-1*rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t5-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t3*t6-1.240700981798
     1         8002d+0*rhob**3.333333333333333d-1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t7
     1         *wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t4*t8
     1         *wght
            Amat2(iq,D2_RA_RA) = (-4.135669939329334d-1*t2-8.11030008110
     1         3d-3*t1/rhoa**2.6666666666666666d+0)*wght+Amat2(iq,D2_RA_
     2         RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-4.135669939329334d-1*t4-8.11030008110
     1         3d-3*t3/rhob**2.6666666666666666d+0)*wght+Amat2(iq,D2_RB_
     2         RB)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t5*t7*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 3.64963503649635d-3*t6*t8*wght+Cmat2(i
     1         q,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t2*wght/gammaaa*
     1         *1.25d+0+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.368613138686131d-3*t4*wght/gammabb*
     1         *1.25d+0+Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 1/gammaaa**2.5d-1
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         oa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t4
     1         *wght
            Amat2(iq,D2_RA_RA) = (-4.135669939329334d-1*t2-8.11030008110
     1         3d-3*t1/rhoa**2.6666666666666666d+0)*wght+Amat2(iq,D2_RA_
     2         RA)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t2*wght/gammaaa*
     1         *1.25d+0+Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**7.5d-1
            t2 = 1/rhob**6.666666666666666d-1
            t3 = 1/rhob**1.6666666666666669d+0
            t4 = 1/gammabb**2.5d-1
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         ob**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhob**3.333333333333333d-1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t2*t4
     1         *wght
            Amat2(iq,D2_RB_RB) = (-4.135669939329334d-1*t2-8.11030008110
     1         3d-3*t1/rhob**2.6666666666666666d+0)*wght+Amat2(iq,D2_RB_
     2         RB)
            Cmat2(iq,D2_RB_GBB) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.368613138686131d-3*t2*wght/gammabb*
     1         *1.25d+0+Cmat2(iq,D2_GBB_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_gill functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}
C>    \over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}
C>    -{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}
C>    \over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}
C>    -0.9305257363491002\,\rho_\beta^{{{4}\over{3}}}
C>    -0.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}
C>    \over{\rho_s^{{{2}\over{3}}}}}-0.9305257363491002
C>    \,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> Code generated with Maxima 5.34.0 [2] 
C> driven by autoxc [3]. 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="https://doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
C> [2] Maxima, a computer algebra system, 
C> <a href="http://maxima.sourceforge.net/"> 
C> http://maxima.sourceforge.net/</a> 
C> 
C> [3] autoxc, revision 27097 2015-05-08 
C> 
      subroutine nwxcm_x_gill_d3(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 1/gammaaa**2.5d-1
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 1/gammaaa**1.25d+0
            fnc(iq) = (-1.45985401459854d-2*t1*t2-1.8610514726982003d+0*
     1         rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t4
     1         *wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = (-8.110300081103d-3*t1*t5-4.13566993932
     1         9334d-1*t2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t2*t6*wght+Cmat2
     1         (iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = (2.7571132928862224d-1*t3+2.16274668
     1         82941332d-2*t1/rhoa**3.6666666666666664d+0)*wght+Amat3(iq
     2         ,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)-6.0827250608
     1         2725d-3*t4*t5*wght
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)-9.12408759
     1         1240875d-4*t3*t6*wght
            Cmat3(iq,D3_RA_GAA_GAB) = Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-1.710766
     1         4233576642d-3*t2*wght/gammaaa**2.25d+0
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = gammabb**7.5d-1
            t4 = 1/rhob**6.666666666666666d-1
            t5 = 1/rhoa**1.6666666666666669d+0
            t6 = 1/rhob**1.6666666666666669d+0
            t7 = 1/gammaaa**2.5d-1
            t8 = 1/gammabb**2.5d-1
            t9 = 1/rhoa**2.6666666666666666d+0
            t10 = 1/rhob**2.6666666666666666d+0
            t11 = 1/gammaaa**1.25d+0
            t12 = 1/gammabb**1.25d+0
            fnc(iq) = (-7.2992700729927d-3*t3*t4-7.2992700729927d-3*t1*t
     1         2-9.305257363491002d-1*rhob**1.3333333333333333d+0-9.3052
     2         57363491002d-1*rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t5-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t3*t6-1.240700981798
     1         8002d+0*rhob**3.333333333333333d-1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t7
     1         *wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t4*t8
     1         *wght
            Amat2(iq,D2_RA_RA) = (-8.110300081103d-3*t1*t9-4.13566993932
     1         9334d-1*t2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-4.135669939329334d-1*t4-8.11030008110
     1         3d-3*t10*t3)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t5*t7*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = 3.64963503649635d-3*t6*t8*wght+Cmat2(i
     1         q,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t11*t2*wght+Cmat
     1         2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.368613138686131d-3*t12*t4*wght+Cmat
     1         2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (2.7571132928862224d-1*t5+2.16274668
     1         82941332d-2*t1/rhoa**3.6666666666666664d+0)*wght+Amat3(iq
     2         ,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (2.7571132928862224d-1*t6+2.16274668
     1         82941332d-2*t3/rhob**3.6666666666666664d+0)*wght+Amat3(iq
     2         ,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)-6.0827250608
     1         2725d-3*t7*t9*wght
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)-6.0827250608
     1         2725d-3*t10*t8*wght
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)-9.12408759
     1         1240875d-4*t11*t5*wght
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
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)-9.12408759
     1         1240875d-4*t12*t6*wght
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-1.710766
     1         4233576642d-3*t2*wght/gammaaa**2.25d+0
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-1.710766
     1         4233576642d-3*t4*wght/gammabb**2.25d+0
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = 1/rhoa**1.6666666666666669d+0
            t4 = 1/gammaaa**2.5d-1
            t5 = 1/rhoa**2.6666666666666666d+0
            t6 = 1/gammaaa**1.25d+0
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         oa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhoa**3.333333333333333d-1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-5.474452554744524d-3*t2*t4
     1         *wght
            Amat2(iq,D2_RA_RA) = (-8.110300081103d-3*t1*t5-4.13566993932
     1         9334d-1*t2)*wght+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.368613138686131d-3*t2*t6*wght+Cmat2
     1         (iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = (2.7571132928862224d-1*t3+2.16274668
     1         82941332d-2*t1/rhoa**3.6666666666666664d+0)*wght+Amat3(iq
     2         ,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = Cmat3(iq,D3_RA_RA_GAA)-6.0827250608
     1         2725d-3*t4*t5*wght
            Cmat3(iq,D3_RA_GAA_GAA) = Cmat3(iq,D3_RA_GAA_GAA)-9.12408759
     1         1240875d-4*t3*t6*wght
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-1.710766
     1         4233576642d-3*t2*wght/gammaaa**2.25d+0
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**7.5d-1
            t2 = 1/rhob**6.666666666666666d-1
            t3 = 1/rhob**1.6666666666666669d+0
            t4 = 1/gammabb**2.5d-1
            t5 = 1/rhob**2.6666666666666666d+0
            t6 = 1/gammabb**1.25d+0
            fnc(iq) = (-7.2992700729927d-3*t1*t2-9.305257363491002d-1*rh
     1         ob**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RB) = (4.8661800486617995d-3*t1*t3-1.240700981798
     1         8002d+0*rhob**3.333333333333333d-1)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-5.474452554744524d-3*t2*t4
     1         *wght
            Amat2(iq,D2_RB_RB) = (-8.110300081103d-3*t1*t5-4.13566993932
     1         9334d-1*t2)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 3.64963503649635d-3*t3*t4*wght+Cmat2(i
     1         q,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.368613138686131d-3*t2*t6*wght+Cmat2
     1         (iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = (2.7571132928862224d-1*t3+2.16274668
     1         82941332d-2*t1/rhob**3.6666666666666664d+0)*wght+Amat3(iq
     2         ,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = Cmat3(iq,D3_RB_RB_GBB)-6.0827250608
     1         2725d-3*t4*t5*wght
            Cmat3(iq,D3_RB_GBB_GBB) = Cmat3(iq,D3_RB_GBB_GBB)-9.12408759
     1         1240875d-4*t3*t6*wght
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-1.710766
     1         4233576642d-3*t2*wght/gammabb**2.25d+0
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
