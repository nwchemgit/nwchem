C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_xc_kt1.F 
C> The nwxcm_xc_kt1 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
      subroutine nwxcm_xc_kt1(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = -1.2d-2*t2*wght
            fnc(iq) = fnc(iq)-1.2d-2*gammaaa*t2*wght
            Amat(iq,D1_RA) = 1.6d-2*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
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
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            t3 = rhob**1.3333333333333333d+0+1.0d-1
            t4 = 1/t3
            fnc(iq) = (-6.0d-3*gammabb*t4-6.0d-3*gammaaa*t2)*wght+fnc(iq
     1         )
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = 8.0d-3*gammabb*rhob**3.333333333333333d-1*w
     1         ght/t3**2+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t4*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            fnc(iq) = fnc(iq)-6.0d-3*gammaaa*t2*wght
            Amat(iq,D1_RA) = 8.0d-3*gammaaa*rhoa**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-6.0d-3*t2*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0+1.0d-1
            t2 = 1/t1
            fnc(iq) = fnc(iq)-6.0d-3*gammabb*t2*wght
            Amat(iq,D1_RB) = 8.0d-3*gammabb*rhob**3.333333333333333d-1*w
     1         ght/t1**2+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-6.0d-3*t2*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
      subroutine nwxcm_xc_kt1_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_xc_kt1 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.006\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}+0.1}}-{{0.006\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}+0.1}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.006\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}+0.1}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] TW Keal, DJ Tozer, J.Chem.Phys. 119, 3015 (2003)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.1590634 "> 
C> 10.1063/1.1590634 </a> 
C> 
      subroutine nwxcm_xc_kt1_d3(param,tol_rho,ipol,nq,wght, 
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
