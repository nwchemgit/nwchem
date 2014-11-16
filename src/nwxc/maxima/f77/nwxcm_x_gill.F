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
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}\over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}-{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}\over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}-.9305257363491002\,\rho_\beta^{{{4}\over{3}}}-.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}\over{\rho_s^{{{2}\over{3}}}}}-.9305257363491002\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
      subroutine nwxcm_x_gill(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**7.5d-1
            t2 = 1/rhoa**6.666666666666666d-1
            t3 = -1.094890510948905d-2*t2*wght/gammaaa**2.5d-1
            fnc(iq) = (-1.45985401459854d-2*t1*t2-1.8610514726982003d+0*
     1         rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = (9.732360097323599d-3*t1/rhoa**1.6666666666
     1         666669d+0-2.4814019635976003d+0*rhoa**3.333333333333333d-
     2         1)*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t3+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t3+Cmat(iq,D1_GAB)
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
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}\over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}-{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}\over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}-.9305257363491002\,\rho_\beta^{{{4}\over{3}}}-.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}\over{\rho_s^{{{2}\over{3}}}}}-.9305257363491002\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
      subroutine nwxcm_x_gill_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_gill functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.0072992700729927\,\sigma_{\beta\beta}^{{{3}\over{4}}}}\over{\rho_\beta^{{{2}\over{3}}}}}-{{0.0072992700729927\,\sigma_{\alpha\alpha}^{{{3}\over{4}}}}\over{\rho_\alpha^{{{2}\over{3}}}}}-.9305257363491002\,\rho_\beta^{{{4}\over{3}}}-.9305257363491002\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0072992700729927\,\sigma_{ss}^{{{3}\over{4}}}}\over{\rho_s^{{{2}\over{3}}}}}-.9305257363491002\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] PMW Gill, Mol.Phys. 89, 433 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/002689796173813 "> 
C> 10.1080/002689796173813 </a> 
C> 
      subroutine nwxcm_x_gill_d3(param,tol_rho,ipol,nq,wght, 
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
