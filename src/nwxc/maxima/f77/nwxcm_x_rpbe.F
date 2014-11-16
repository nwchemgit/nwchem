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
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59, 7413 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
      subroutine nwxcm_x_rpbe(param,tol_rho,ipol,nq,wght, 
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
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = -6.722464919227583d-3*t2*wght/t1
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = -1.9950471787324703d+0*rhoa**3.333333333333
     1         333d-1*t3*wght+1.792657311794022d-2*gammaaa*t2*wght/rhoa*
     2         *2.3333333333333334d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t4+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t4+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = rhob**1.3333333333333333d+0
            t5 = exp(-4.492769220958889d-3*gammabb/rhob**2.6666666666666
     1         666d+0)
            t6 = 1.0d+0-t5
            fnc(iq) = (-7.481426920246764d-1*t4*t6-7.481426920246764d-1*
     1         t1*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (8.96328655897011d-3*gammaaa*t2/rhoa**2.333
     1         3333333333334d+0-9.975235893662352d-1*rhoa**3.33333333333
     2         3333d-1*t3)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (8.96328655897011d-3*gammabb*t5/rhob**2.333
     1         3333333333334d+0-9.975235893662352d-1*rhob**3.33333333333
     2         3333d-1*t6)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t5*w
     1         ght/t4
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = -9.975235893662352d-1*rhoa**3.3333333333333
     1         33d-1*t3*wght+8.96328655897011d-3*gammaaa*t2*wght/rhoa**2
     2         .3333333333333334d+0+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammabb/rhob**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = -9.975235893662352d-1*rhob**3.3333333333333
     1         33d-1*t3*wght+8.96328655897011d-3*gammabb*t2*wght/rhob**2
     2         .3333333333333334d+0+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t2*w
     1         ght/t1
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_rpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59, 7413 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
      subroutine nwxcm_x_rpbe_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_rpbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.7481426920246764\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}} }\right)-0.7481426920246764\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}} }\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.7481426920246764\,\rho_s^{{{4}\over{3}}}\,\left(1.0-e^ {- {{0.004492769220958889\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}} }\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] B Hammer, LB Hanssen, JK Norskov, Phys.Rev.B 59, 7413 (1999)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevB.59.7413 "> 
C> 10.1103/PhysRevB.59.7413 </a> 
C> 
      subroutine nwxcm_x_rpbe_d3(param,tol_rho,ipol,nq,wght, 
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
