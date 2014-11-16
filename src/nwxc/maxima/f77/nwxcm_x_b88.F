C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_b88.F 
C> The nwxcm_x_b88 functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252\,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)\,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}\,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}\,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
      subroutine nwxcm_x_b88(param,tol_rho,ipol,nq,wght, 
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
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 1/(gammaaa*t7+1)**5.0d-1
            t9 = 1/t4**2
            t10 = 8.4d-3*gammaaa*t2*(1.26d-2*t7*t8+1.26d-2*t2*t3/t1)*t9*
     1         wght-8.4d-3*t2*t5*wght
            fnc(iq) = fnc(iq)-8.4d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = 8.4d-3*gammaaa*t2*(-3.3599999999999997d-2*g
     1         ammaaa*t8/rhoa**3.6666666666666664d+0-3.3599999999999997d
     2         -2*t1*t3*t6)*t9*wght+1.1199999999999999d-2*gammaaa*t5*t6*
     3         wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t10+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t10+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = gammabb**5.0d-1
            t7 = 1/rhob**1.3333333333333333d+0
            t8 = asinh(t6*t7)
            t9 = 2.52d-2*t6*t7*t8+1.0d+0
            t10 = 1/t9
            t11 = 1/rhoa**2.3333333333333334d+0
            t12 = 1/rhoa**2.6666666666666666d+0
            t13 = 1/(gammaaa*t12+1)**5.0d-1
            t14 = 1/t4**2
            t15 = 1/rhob**2.3333333333333334d+0
            t16 = 1/rhob**2.6666666666666666d+0
            t17 = 1/(gammabb*t16+1)**5.0d-1
            t18 = 1/t9**2
            fnc(iq) = (-4.2d-3*gammabb*t10*t7-4.2d-3*gammaaa*t2*t5)*wght
     1         +fnc(iq)
            Amat(iq,D1_RA) = (5.599999999999999d-3*gammaaa*t11*t5+4.2d-3
     1         *gammaaa*t14*t2*(-3.3599999999999997d-2*t1*t11*t3-3.35999
     2         99999999997d-2*gammaaa*t13/rhoa**3.6666666666666664d+0))*
     3         wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.2d-3*gammabb*t18*t7*(-3.3599999999999997
     1         d-2*t15*t6*t8-3.3599999999999997d-2*gammabb*t17/rhob**3.6
     2         666666666666664d+0)+5.599999999999999d-3*gammabb*t10*t15)
     3         *wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = (4.2d-3*gammaaa*t14*t2*(1.26d-2*t2*t3/t1+1
     1         .26d-2*t12*t13)-4.2d-3*t2*t5)*wght+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = (4.2d-3*gammabb*t18*t7*(1.26d-2*t7*t8/t6+1
     1         .26d-2*t16*t17)-4.2d-3*t10*t7)*wght+Cmat(iq,D1_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = gammaaa**5.0d-1
            t2 = 1/rhoa**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhoa**2.3333333333333334d+0
            t7 = 1/rhoa**2.6666666666666666d+0
            t8 = 1/(gammaaa*t7+1)**5.0d-1
            t9 = 1/t4**2
            fnc(iq) = fnc(iq)-4.2d-3*gammaaa*t2*t5*wght
            Amat(iq,D1_RA) = 4.2d-3*gammaaa*t2*(-3.3599999999999997d-2*g
     1         ammaaa*t8/rhoa**3.6666666666666664d+0-3.3599999999999997d
     2         -2*t1*t3*t6)*t9*wght+5.599999999999999d-3*gammaaa*t5*t6*w
     3         ght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = 4.2d-3*gammaaa*t2*(1.26d-2*t7*t8+1.26d-2*t
     1         2*t3/t1)*t9*wght-4.2d-3*t2*t5*wght+Cmat(iq,D1_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = gammabb**5.0d-1
            t2 = 1/rhob**1.3333333333333333d+0
            t3 = asinh(t1*t2)
            t4 = 2.52d-2*t1*t2*t3+1.0d+0
            t5 = 1/t4
            t6 = 1/rhob**2.3333333333333334d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 1/(gammabb*t7+1)**5.0d-1
            t9 = 1/t4**2
            fnc(iq) = fnc(iq)-4.2d-3*gammabb*t2*t5*wght
            Amat(iq,D1_RB) = 4.2d-3*gammabb*t2*(-3.3599999999999997d-2*g
     1         ammabb*t8/rhob**3.6666666666666664d+0-3.3599999999999997d
     2         -2*t1*t3*t6)*t9*wght+5.599999999999999d-3*gammabb*t5*t6*w
     3         ght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = 4.2d-3*gammabb*t2*(1.26d-2*t7*t8+1.26d-2*t
     1         2*t3/t1)*t9*wght-4.2d-3*t2*t5*wght+Cmat(iq,D1_GBB)
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252\,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)\,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}\,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}\,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
      subroutine nwxcm_x_b88_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_b88 functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& \sqrt{\sigma_{\alpha\alpha}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_\beta^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_4} &=& \sqrt{\sigma_{\beta\beta}}\\\\ 
C>   {\it t_5} &=& {{1}\over{\rho_s^{{{4}\over{3}}}}}\\\\ 
C>   {\it t_6} &=& \sqrt{\sigma_{ss}}\\\\ 
C>   f &=& -{{0.0042\,{\it t_3}\,\sigma_{\beta\beta}}\over{0.0252\,{\it t_3}\,{\rm asinh}\; \left({\it t_3}\,{\it t_4}\right)\,{\it t_4}+1.0}}-{{0.0042\,{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.0252\,{\it t_1}\,{\rm asinh}\; \left({\it t_1}\,{\it t_2}\right)\,{\it t_2}+1.0}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -{{0.0042\,{\it t_5}\,\sigma_{ss}}\over{0.0252\,{\it t_5}\,{\rm asinh}\; \left({\it t_5}\,{\it t_6}\right)\,{\it t_6}+1.0}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] AD Becke, Phys.Rev.A 38, 3098 (1988)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevA.38.3098 "> 
C> 10.1103/PhysRevA.38.3098 </a> 
C> 
      subroutine nwxcm_x_b88_d3(param,tol_rho,ipol,nq,wght, 
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
