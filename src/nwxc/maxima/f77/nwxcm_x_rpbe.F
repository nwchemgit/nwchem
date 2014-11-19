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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = -6.722464919227583d-3*t2*wght/t1
            t6 = -6.650157262441567d-1*t3*wght/rhoa**6.666666666666666d-
     1         1-1.7926573117940228d-2*gammaaa*t2*wght/rhoa**3.333333333
     2         3333337d+0+2.1477321584413558d-4*gammaaa**2*t2*wght/rhoa*
     3         *6
            t7 = 8.96328655897011d-3*t2*t4*wght-8.053995594155084d-5*gam
     1         maaa*t2*wght/rhoa**5
            t8 = 3.0202483478081565d-5*t2*wght/rhoa**4
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = 1.792657311794022d-2*gammaaa*t2*t4*wght-1.9
     1         950471787324703d+0*rhoa**3.333333333333333d-1*t3*wght+Ama
     2         t(iq,D1_RA)
            Cmat(iq,D1_GAA) = t5+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t5+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = t6+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t6+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t7+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = t7+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t7+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t8+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t8+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t8+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = t8+Cmat2(iq,D2_GAB_GAB)
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
            t7 = 1/rhoa**2.3333333333333334d+0
            t8 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = (-7.481426920246764d-1*t4*t6-7.481426920246764d-1*
     1         t1*t3)*wght+fnc(iq)
            Amat(iq,D1_RA) = (8.96328655897011d-3*gammaaa*t2*t7-9.975235
     1         893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq,
     2         D1_RA)
            Amat(iq,D1_RB) = (8.96328655897011d-3*gammabb*t5*t8-9.975235
     1         893662352d-1*rhob**3.333333333333333d-1*t6)*wght+Amat(iq,
     2         D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t5*w
     1         ght/t4
            Amat2(iq,D2_RA_RA) = (-3.3250786312207836d-1*t3/rhoa**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammaaa*t2/rhoa**3.33
     2         33333333333337d+0+1.0738660792206779d-4*gammaaa**2*t2/rho
     3         a**6)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-3.3250786312207836d-1*t6/rhob**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammabb*t5/rhob**3.33
     2         33333333333337d+0+1.0738660792206779d-4*gammabb**2*t5/rho
     3         b**6)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (4.481643279485055d-3*t2*t7-4.02699779
     1         7077542d-5*gammaaa*t2/rhoa**5)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.481643279485055d-3*t5*t8-4.02699779
     1         7077542d-5*gammabb*t5/rhob**5)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.5101241739040783d-5*t2*wght/rhoa**4
     1         +Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5101241739040783d-5*t5*wght/rhob**4
     1         +Cmat2(iq,D2_GBB_GBB)
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = 8.96328655897011d-3*gammaaa*t2*t4*wght-9.97
     1         5235893662352d-1*rhoa**3.333333333333333d-1*t3*wght+Amat(
     2         iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Amat2(iq,D2_RA_RA) = -3.3250786312207836d-1*t3*wght/rhoa**6.
     1         666666666666666d-1-8.963286558970114d-3*gammaaa*t2*wght/r
     2         hoa**3.3333333333333337d+0+1.0738660792206779d-4*gammaaa*
     3         *2*t2*wght/rhoa**6+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = 4.481643279485055d-3*t2*t4*wght-4.0269
     1         97797077542d-5*gammaaa*t2*wght/rhoa**5+Cmat2(iq,D2_RA_GAA
     2         )
            Cmat2(iq,D2_GAA_GAA) = 1.5101241739040783d-5*t2*wght/rhoa**4
     1         +Cmat2(iq,D2_GAA_GAA)
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammabb/rhob**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhob**2.3333333333333334d+0
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = 8.96328655897011d-3*gammabb*t2*t4*wght-9.97
     1         5235893662352d-1*rhob**3.333333333333333d-1*t3*wght+Amat(
     2         iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Amat2(iq,D2_RB_RB) = -3.3250786312207836d-1*t3*wght/rhob**6.
     1         666666666666666d-1-8.963286558970114d-3*gammabb*t2*wght/r
     2         hob**3.3333333333333337d+0+1.0738660792206779d-4*gammabb*
     3         *2*t2*wght/rhob**6+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = 4.481643279485055d-3*t2*t4*wght-4.0269
     1         97797077542d-5*gammabb*t2*wght/rhob**5+Cmat2(iq,D2_RB_GBB
     2         )
            Cmat2(iq,D2_GBB_GBB) = 1.5101241739040783d-5*t2*wght/rhob**4
     1         +Cmat2(iq,D2_GBB_GBB)
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
      double precision nwxcm_heaviside 
      external         nwxcm_heaviside 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          gammaaa = 0.25d0*rgamma(iq,G_TT) 
          if (rhoa.gt.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = -6.722464919227583d-3*t2*wght/t1
            t6 = gammaaa**2
            t7 = 1/rhoa**6
            t8 = 1/rhoa**3.3333333333333337d+0
            t9 = -1.7926573117940228d-2*gammaaa*t2*t8*wght+2.14773215844
     1         13558d-4*t2*t6*t7*wght-6.650157262441567d-1*t3*wght/rhoa*
     2         *6.666666666666666d-1
            t10 = 1/rhoa**5
            t11 = 8.96328655897011d-3*t2*t4*wght-8.053995594155084d-5*ga
     1         mmaaa*t10*t2*wght
            t12 = 3.0202483478081565d-5*t2*wght/rhoa**4
            t13 = -1.503412510908949d-3*t2*t6*wght/rhoa**7+4.43343817496
     1         10446d-1*t3*wght/rhoa**1.6666666666666669d+0+6.7722609556
     2         66308d-2*gammaaa*t2*wght/rhoa**4.333333333333333d+0+2.573
     3         137316349046d-6*gammaaa**3*t2*wght/rhoa**9.66666666666666
     4         6d+0
            t14 = -2.0914335304263598d-2*t2*t8*wght+5.10086387629822d-4*
     1         gammaaa*t2*t7*wght-9.649264936308923d-7*t2*t6*wght/rhoa**
     2         8.666666666666666d+0
            t15 = 3.618474351115847d-7*gammaaa*t2*wght/rhoa**7.666666666
     1         666667d+0-1.2080993391232628d-4*t10*t2*wght
            t16 = -1.3569278816684424d-7*t2*wght/rhoa**6.666666666666667
     1         d+0
            fnc(iq) = fnc(iq)-1.4962853840493529d+0*t1*t3*wght
            Amat(iq,D1_RA) = 1.792657311794022d-2*gammaaa*t2*t4*wght-1.9
     1         950471787324703d+0*rhoa**3.333333333333333d-1*t3*wght+Ama
     2         t(iq,D1_RA)
            Cmat(iq,D1_GAA) = t5+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t5+Cmat(iq,D1_GAB)
            Amat2(iq,D2_RA_RA) = t9+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = t9+Amat2(iq,D2_RA_RB)
            Cmat2(iq,D2_RA_GAA) = t11+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = t11+Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = t11+Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_GAA_GAA) = t12+Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = t12+Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = t12+Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = t12+Cmat2(iq,D2_GAB_GAB)
            Amat3(iq,D3_RA_RA_RA) = t13+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = t13+Amat3(iq,D3_RA_RA_RB)
            Cmat3(iq,D3_RA_RA_GAA) = t14+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = t14+Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = t14+Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = t14+Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = t14+Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_GAA_GAA) = t15+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_RA_GAA_GAB) = t15+Cmat3(iq,D3_RA_GAA_GAB)
            Cmat3(iq,D3_RA_GAA_GBB) = t15+Cmat3(iq,D3_RA_GAA_GBB)
            Cmat3(iq,D3_RA_GAB_GAB) = t15+Cmat3(iq,D3_RA_GAB_GAB)
            Cmat3(iq,D3_RA_GAB_GBB) = t15+Cmat3(iq,D3_RA_GAB_GBB)
            Cmat3(iq,D3_RA_GBB_GBB) = t15+Cmat3(iq,D3_RA_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = t16+Cmat3(iq,D3_GAA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAB) = t16+Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = t16+Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = t16+Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = t16+Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = t16+Cmat3(iq,D3_GAB_GAB_GAB)
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
            Amat(iq,D1_RA) = (8.96328655897011d-3*gammaaa*t2*t7-9.975235
     1         893662352d-1*rhoa**3.333333333333333d-1*t3)*wght+Amat(iq,
     2         D1_RA)
            Amat(iq,D1_RB) = (8.96328655897011d-3*gammabb*t5*t8-9.975235
     1         893662352d-1*rhob**3.333333333333333d-1*t6)*wght+Amat(iq,
     2         D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t5*w
     1         ght/t4
            Amat2(iq,D2_RA_RA) = (1.0738660792206779d-4*t10*t2*t9-3.3250
     1         786312207836d-1*t3/rhoa**6.666666666666666d-1-8.963286558
     2         970114d-3*gammaaa*t11*t2)*wght+Amat2(iq,D2_RA_RA)
            Amat2(iq,D2_RA_RB) = Amat2(iq,D2_RA_RB)
            Amat2(iq,D2_RB_RB) = (-3.3250786312207836d-1*t6/rhob**6.6666
     1         66666666666d-1-8.963286558970114d-3*gammabb*t14*t5+1.0738
     2         660792206779d-4*t12*t13*t5)*wght+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RA_GAA) = (4.481643279485055d-3*t2*t7-4.02699779
     1         7077542d-5*gammaaa*t15*t2)*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_RA_GAB) = Cmat2(iq,D2_RA_GAB)
            Cmat2(iq,D2_RA_GBB) = Cmat2(iq,D2_RA_GBB)
            Cmat2(iq,D2_RB_GAA) = Cmat2(iq,D2_RB_GAA)
            Cmat2(iq,D2_RB_GAB) = Cmat2(iq,D2_RB_GAB)
            Cmat2(iq,D2_RB_GBB) = (4.481643279485055d-3*t5*t8-4.02699779
     1         7077542d-5*gammabb*t16*t5)*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GAA_GAA) = 1.5101241739040783d-5*t2*wght/rhoa**4
     1         +Cmat2(iq,D2_GAA_GAA)
            Cmat2(iq,D2_GAA_GAB) = Cmat2(iq,D2_GAA_GAB)
            Cmat2(iq,D2_GAA_GBB) = Cmat2(iq,D2_GAA_GBB)
            Cmat2(iq,D2_GAB_GAB) = Cmat2(iq,D2_GAB_GAB)
            Cmat2(iq,D2_GAB_GBB) = Cmat2(iq,D2_GAB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5101241739040783d-5*t5*wght/rhob**4
     1         +Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RA_RA_RA) = (-7.517062554544747d-4*t2*t9/rhoa**7
     1         +2.2167190874805223d-1*t3/rhoa**1.6666666666666669d+0+3.3
     2         86130477833154d-2*gammaaa*t2/rhoa**4.333333333333333d+0+1
     3         .286568658174523d-6*gammaaa**3*t2/rhoa**9.666666666666666
     4         d+0)*wght+Amat3(iq,D3_RA_RA_RA)
            Amat3(iq,D3_RA_RA_RB) = Amat3(iq,D3_RA_RA_RB)
            Amat3(iq,D3_RA_RB_RB) = Amat3(iq,D3_RA_RB_RB)
            Amat3(iq,D3_RB_RB_RB) = (2.2167190874805223d-1*t6/rhob**1.66
     1         66666666666669d+0-7.517062554544747d-4*t12*t5/rhob**7+3.3
     2         86130477833154d-2*gammabb*t5/rhob**4.333333333333333d+0+1
     3         .286568658174523d-6*gammabb**3*t5/rhob**9.666666666666666
     4         d+0)*wght+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RA_RA_GAA) = (-4.8246324681544617d-7*t2*t9/rhoa*
     1         *8.666666666666666d+0-1.0457167652131799d-2*t11*t2+2.5504
     2         3193814911d-4*gammaaa*t10*t2)*wght+Cmat3(iq,D3_RA_RA_GAA)
            Cmat3(iq,D3_RA_RA_GAB) = Cmat3(iq,D3_RA_RA_GAB)
            Cmat3(iq,D3_RA_RA_GBB) = Cmat3(iq,D3_RA_RA_GBB)
            Cmat3(iq,D3_RA_RB_GAA) = Cmat3(iq,D3_RA_RB_GAA)
            Cmat3(iq,D3_RA_RB_GAB) = Cmat3(iq,D3_RA_RB_GAB)
            Cmat3(iq,D3_RA_RB_GBB) = Cmat3(iq,D3_RA_RB_GBB)
            Cmat3(iq,D3_RB_RB_GAA) = Cmat3(iq,D3_RB_RB_GAA)
            Cmat3(iq,D3_RB_RB_GAB) = Cmat3(iq,D3_RB_RB_GAB)
            Cmat3(iq,D3_RB_RB_GBB) = (-1.0457167652131799d-2*t14*t5+2.55
     1         043193814911d-4*gammabb*t13*t5-4.8246324681544617d-7*t12*
     2         t5/rhob**8.666666666666666d+0)*wght+Cmat3(iq,D3_RB_RB_GBB
     3         )
            Cmat3(iq,D3_RA_GAA_GAA) = (1.8092371755579234d-7*gammaaa*t2/
     1         rhoa**7.666666666666667d+0-6.040496695616313d-5*t15*t2)*w
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
            Cmat3(iq,D3_RB_GBB_GBB) = (1.8092371755579234d-7*gammabb*t5/
     1         rhob**7.666666666666667d+0-6.040496695616313d-5*t16*t5)*w
     2         ght+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-6.784639
     1         408342211d-8*t2*wght/rhoa**6.666666666666667d+0
            Cmat3(iq,D3_GAA_GAA_GAB) = Cmat3(iq,D3_GAA_GAA_GAB)
            Cmat3(iq,D3_GAA_GAA_GBB) = Cmat3(iq,D3_GAA_GAA_GBB)
            Cmat3(iq,D3_GAA_GAB_GAB) = Cmat3(iq,D3_GAA_GAB_GAB)
            Cmat3(iq,D3_GAA_GAB_GBB) = Cmat3(iq,D3_GAA_GAB_GBB)
            Cmat3(iq,D3_GAA_GBB_GBB) = Cmat3(iq,D3_GAA_GBB_GBB)
            Cmat3(iq,D3_GAB_GAB_GAB) = Cmat3(iq,D3_GAB_GAB_GAB)
            Cmat3(iq,D3_GAB_GAB_GBB) = Cmat3(iq,D3_GAB_GAB_GBB)
            Cmat3(iq,D3_GAB_GBB_GBB) = Cmat3(iq,D3_GAB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-6.784639
     1         408342211d-8*t5*wght/rhob**6.666666666666667d+0
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = rhoa**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammaaa/rhoa**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhoa**2.3333333333333334d+0
            t5 = gammaaa**2
            t6 = 1/rhoa**6
            t7 = 1/rhoa**3.3333333333333337d+0
            t8 = 1/rhoa**5
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RA) = 8.96328655897011d-3*gammaaa*t2*t4*wght-9.97
     1         5235893662352d-1*rhoa**3.333333333333333d-1*t3*wght+Amat(
     2         iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Amat2(iq,D2_RA_RA) = -8.963286558970114d-3*gammaaa*t2*t7*wgh
     1         t+1.0738660792206779d-4*t2*t5*t6*wght-3.3250786312207836d
     2         -1*t3*wght/rhoa**6.666666666666666d-1+Amat2(iq,D2_RA_RA)
            Cmat2(iq,D2_RA_GAA) = -4.026997797077542d-5*gammaaa*t2*t8*wg
     1         ht+4.481643279485055d-3*t2*t4*wght+Cmat2(iq,D2_RA_GAA)
            Cmat2(iq,D2_GAA_GAA) = 1.5101241739040783d-5*t2*wght/rhoa**4
     1         +Cmat2(iq,D2_GAA_GAA)
            Amat3(iq,D3_RA_RA_RA) = -7.517062554544747d-4*t2*t5*wght/rho
     1         a**7+2.2167190874805223d-1*t3*wght/rhoa**1.66666666666666
     2         69d+0+3.386130477833154d-2*gammaaa*t2*wght/rhoa**4.333333
     3         333333333d+0+1.286568658174523d-6*gammaaa**3*t2*wght/rhoa
     4         **9.666666666666666d+0+Amat3(iq,D3_RA_RA_RA)
            Cmat3(iq,D3_RA_RA_GAA) = -1.0457167652131799d-2*t2*t7*wght+2
     1         .55043193814911d-4*gammaaa*t2*t6*wght-4.8246324681544617d
     2         -7*t2*t5*wght/rhoa**8.666666666666666d+0+Cmat3(iq,D3_RA_R
     3         A_GAA)
            Cmat3(iq,D3_RA_GAA_GAA) = -6.040496695616313d-5*t2*t8*wght+1
     1         .8092371755579234d-7*gammaaa*t2*wght/rhoa**7.666666666666
     2         667d+0+Cmat3(iq,D3_RA_GAA_GAA)
            Cmat3(iq,D3_GAA_GAA_GAA) = Cmat3(iq,D3_GAA_GAA_GAA)-6.784639
     1         408342211d-8*t2*wght/rhoa**6.666666666666667d+0
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = rhob**1.3333333333333333d+0
            t2 = exp(-4.492769220958889d-3*gammabb/rhob**2.6666666666666
     1         666d+0)
            t3 = 1.0d+0-t2
            t4 = 1/rhob**2.3333333333333334d+0
            t5 = gammabb**2
            t6 = 1/rhob**6
            t7 = 1/rhob**3.3333333333333337d+0
            t8 = 1/rhob**5
            fnc(iq) = fnc(iq)-7.481426920246764d-1*t1*t3*wght
            Amat(iq,D1_RB) = 8.96328655897011d-3*gammabb*t2*t4*wght-9.97
     1         5235893662352d-1*rhob**3.333333333333333d-1*t3*wght+Amat(
     2         iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-3.3612324596137916d-3*t2*w
     1         ght/t1
            Amat2(iq,D2_RB_RB) = -8.963286558970114d-3*gammabb*t2*t7*wgh
     1         t+1.0738660792206779d-4*t2*t5*t6*wght-3.3250786312207836d
     2         -1*t3*wght/rhob**6.666666666666666d-1+Amat2(iq,D2_RB_RB)
            Cmat2(iq,D2_RB_GBB) = -4.026997797077542d-5*gammabb*t2*t8*wg
     1         ht+4.481643279485055d-3*t2*t4*wght+Cmat2(iq,D2_RB_GBB)
            Cmat2(iq,D2_GBB_GBB) = 1.5101241739040783d-5*t2*wght/rhob**4
     1         +Cmat2(iq,D2_GBB_GBB)
            Amat3(iq,D3_RB_RB_RB) = -7.517062554544747d-4*t2*t5*wght/rho
     1         b**7+2.2167190874805223d-1*t3*wght/rhob**1.66666666666666
     2         69d+0+3.386130477833154d-2*gammabb*t2*wght/rhob**4.333333
     3         333333333d+0+1.286568658174523d-6*gammabb**3*t2*wght/rhob
     4         **9.666666666666666d+0+Amat3(iq,D3_RB_RB_RB)
            Cmat3(iq,D3_RB_RB_GBB) = -1.0457167652131799d-2*t2*t7*wght+2
     1         .55043193814911d-4*gammabb*t2*t6*wght-4.8246324681544617d
     2         -7*t2*t5*wght/rhob**8.666666666666666d+0+Cmat3(iq,D3_RB_R
     3         B_GBB)
            Cmat3(iq,D3_RB_GBB_GBB) = -6.040496695616313d-5*t2*t8*wght+1
     1         .8092371755579234d-7*gammabb*t2*wght/rhob**7.666666666666
     2         667d+0+Cmat3(iq,D3_RB_GBB_GBB)
            Cmat3(iq,D3_GBB_GBB_GBB) = Cmat3(iq,D3_GBB_GBB_GBB)-6.784639
     1         408342211d-8*t2*wght/rhob**6.666666666666667d+0
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
