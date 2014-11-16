C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_opt.F 
C> The nwxcm_x_opt functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}\over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}-3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}\,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}\,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}\,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
      subroutine nwxcm_x_opt(param,tol_rho,ipol,nq,wght, 
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
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t2**2
            t8 = -1.44d-4*t4*t6*(t1*t3-6.0d-3*gammaaa*t7/rhoa**5.3333333
     1         33333333d+0)*wght
            fnc(iq) = fnc(iq)-7.199999999999999d-5*t5*t6*wght
            Amat(iq,D1_RA) = -1.44d-4*gammaaa*t3*(1.6d-2*gammaaa**2*t7/r
     1         hoa**6.333333333333333d+0-2.6666666666666666d+0*gammaaa*t
     2         3/rhoa**3.6666666666666664d+0)*wght/t6-9.599999999999997d
     3         -5*rhoa**3.333333333333333d-1*t5*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t8+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t8+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/rhob**2.6666666666666666d+0
            t8 = 6.0d-3*gammabb*t7+1.0d+0
            t9 = 1/t8
            t10 = gammabb*t7*t9
            t11 = t10**2.0d+0
            t12 = rhob**1.3333333333333333d+0
            t13 = 1/t2**2
            t14 = 1/t8**2
            fnc(iq) = (-3.5999999999999993d-5*t5*t6-3.5999999999999993d-
     1         5*t11*t12)*wght+fnc(iq)
            Amat(iq,D1_RA) = (-7.199999999999999d-5*gammaaa*t3*(1.6d-2*g
     1         ammaaa**2*t13/rhoa**6.333333333333333d+0-2.66666666666666
     2         66d+0*gammaaa*t3/rhoa**3.6666666666666664d+0)/t6-4.799999
     3         9999999987d-5*rhoa**3.333333333333333d-1*t5)*wght+Amat(iq
     4         ,D1_RA)
            Amat(iq,D1_RB) = (-7.199999999999999d-5*gammabb*t9*(1.6d-2*g
     1         ammabb**2*t14/rhob**6.333333333333333d+0-2.66666666666666
     2         66d+0*gammabb*t9/rhob**3.6666666666666664d+0)/t12-4.79999
     3         99999999987d-5*rhob**3.333333333333333d-1*t11)*wght+Amat(
     4         iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*(t1*t
     1         3-6.0d-3*gammaaa*t13/rhoa**5.333333333333333d+0)*t4*t6*wg
     2         ht
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t10*t
     1         12*(t7*t9-6.0d-3*gammabb*t14/rhob**5.333333333333333d+0)*
     2         wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = 1/rhoa**2.6666666666666666d+0
            t2 = 6.0d-3*gammaaa*t1+1.0d+0
            t3 = 1/t2
            t4 = gammaaa*t1*t3
            t5 = t4**2.0d+0
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1/t2**2
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RA) = -7.199999999999999d-5*gammaaa*t3*(1.6d-2*ga
     1         mmaaa**2*t7/rhoa**6.333333333333333d+0-2.6666666666666666
     2         d+0*gammaaa*t3/rhoa**3.6666666666666664d+0)*wght/t6-4.799
     3         9999999999987d-5*rhoa**3.333333333333333d-1*t5*wght+Amat(
     4         iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-7.199999999999999d-5*t4*t6
     1         *(t1*t3-6.0d-3*gammaaa*t7/rhoa**5.333333333333333d+0)*wgh
     2         t
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = 1/rhob**2.6666666666666666d+0
            t2 = 6.0d-3*gammabb*t1+1.0d+0
            t3 = 1/t2
            t4 = gammabb*t1*t3
            t5 = t4**2.0d+0
            t6 = rhob**1.3333333333333333d+0
            t7 = 1/t2**2
            fnc(iq) = fnc(iq)-3.5999999999999993d-5*t5*t6*wght
            Amat(iq,D1_RB) = -7.199999999999999d-5*gammabb*t3*(1.6d-2*ga
     1         mmabb**2*t7/rhob**6.333333333333333d+0-2.6666666666666666
     2         d+0*gammabb*t3/rhob**3.6666666666666664d+0)*wght/t6-4.799
     3         9999999999987d-5*rhob**3.333333333333333d-1*t5*wght+Amat(
     4         iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-7.199999999999999d-5*t4*t6
     1         *(t1*t3-6.0d-3*gammabb*t7/rhob**5.333333333333333d+0)*wgh
     2         t
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}\over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}-3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}\,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}\,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}\,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
      subroutine nwxcm_x_opt_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_opt functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {{1}\over{\rho_\alpha^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_2} &=& {{1}\over{\rho_\beta^{{{8}\over{3}}}}}\\\\ 
C>   {\it t_3} &=& {{1}\over{\rho_s^{{{8}\over{3}}}}}\\\\ 
C>   f &=& -3.599999999999999 \times 10^{-5}\,\rho_\beta^{{{4}\over{3}}}\,\left({{{\it t_2}\,\sigma_{\beta\beta}}\over{0.006\,{\it t_2}\,\sigma_{\beta\beta}+1.0}}\right)^{2.0}-3.599999999999999 \times 10^{-5}\,\rho_\alpha^{{{4}\over{3}}}\,\left({{{\it t_1}\,\sigma_{\alpha\alpha}}\over{0.006\,{\it t_1}\,\sigma_{\alpha\alpha}+1.0}}\right)^{2.0}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -3.599999999999999 \times 10^{-5}\,\rho_s^{{{4}\over{3}}}\,\left({{{\it t_3}\,\sigma_{ss}}\over{0.006\,{\it t_3}\,\sigma_{ss}+1.0}}\right)^{2.0}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] NC Handy, AJ Cohen, Mol.Phys. 99, 403 (2001)  , DOI: 
C> <a href="http://dx.doi.org/10.1080/00268970010018431 "> 
C> 10.1080/00268970010018431 </a> 
C> 
      subroutine nwxcm_x_opt_d3(param,tol_rho,ipol,nq,wght, 
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
