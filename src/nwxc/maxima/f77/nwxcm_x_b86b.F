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
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
      subroutine nwxcm_x_b86b(param,tol_rho,ipol,nq,wght, 
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
            t1 = 7.000000000000001d-3*gammaaa/rhoa**2.6666666666666666d+
     1         0+1.0d+0
            t2 = 1/t1**8.0d-1
            t3 = 1/rhoa**1.3333333333333333d+0
            t4 = 1/t1**1.8d+0
            t5 = 4.2000000000000004d-5*gammaaa*t4*wght/rhoa**4-7.5d-3*t2
     1         *t3*wght
            fnc(iq) = fnc(iq)-7.5d-3*gammaaa*t2*t3*wght
            Amat(iq,D1_RA) = -1.12d-4*gammaaa**2*t4*wght/rhoa**5+9.99999
     1         9999999999d-3*gammaaa*t2*wght/rhoa**2.3333333333333334d+0
     2         +Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t5+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t5+Cmat(iq,D1_GAB)
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
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
      subroutine nwxcm_x_b86b_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_b86b functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -{{0.00375\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}-{{0.00375\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.00375\,\sigma_{ss}}\over{\rho_s^{{{4}\over{3}}}\,\left({{0.007\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0\right)^{0.8}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, J.Chem.Phys. 85, 7184 (1986)  , DOI: 
C> <a href="http://dx.doi.org/10.1063/1.451353 "> 
C> 10.1063/1.451353 </a> 
C> 
      subroutine nwxcm_x_b86b_d3(param,tol_rho,ipol,nq,wght, 
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
