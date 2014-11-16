C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_pbe.F 
C> The nwxcm_x_pbe functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\right)-0.9305257363491\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77, 3865 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
      subroutine nwxcm_x_pbe(param,tol_rho,ipol,nq,wght, 
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
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020562d-2*gammaaa*t2/(rhoa**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhoa**1.3333333333333333d+0
            t6 = 1/t3**2
            t7 = -3.062417490053881d-2*t2*t6*wght/t5
            fnc(iq) = fnc(iq)-1.8610514726982d+0*t1*t4*t5*wght
            Amat(iq,D1_RA) = 8.166446640143683d-2*gammaaa*t2*t6*wght/rho
     1         a**2.3333333333333334d+0-2.4814019635975998d+0*rhoa**3.33
     2         3333333333333d-1*t1*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = t7+Cmat(iq,D1_GAA)
            Cmat(iq,D1_GAB) = t7+Cmat(iq,D1_GAB)
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          gammaaa = rgamma(iq,G_AA) 
          gammaab = rgamma(iq,G_AB) 
          gammabb = rgamma(iq,G_BB) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = 1/t1
            t3 = param(2)
            t4 = 1.6455307846020562d-2*gammaaa*t2*t3/rhoa**2.66666666666
     1         66666d+0+1.0d+0
            t5 = 1.0d+0-1.0d+0/t4
            t6 = rhoa**1.3333333333333333d+0
            t7 = 1.6455307846020562d-2*gammabb*t2*t3/rhob**2.66666666666
     1         66666d+0+1.0d+0
            t8 = 1.0d+0-1.0d+0/t7
            t9 = rhob**1.3333333333333333d+0
            t10 = 1/t4**2
            t11 = 1/t7**2
            fnc(iq) = (-9.305257363491d-1*t1*t8*t9-9.305257363491d-1*t1*
     1         t5*t6)*wght+fnc(iq)
            Amat(iq,D1_RA) = (4.0832233200718415d-2*gammaaa*t10*t3/rhoa*
     1         *2.3333333333333334d+0-1.2407009817987999d+0*rhoa**3.3333
     2         33333333333d-1*t1*t5)*wght+Amat(iq,D1_RA)
            Amat(iq,D1_RB) = (4.0832233200718415d-2*gammabb*t11*t3/rhob*
     1         *2.3333333333333334d+0-1.2407009817987999d+0*rhob**3.3333
     2         33333333333d-1*t1*t8)*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.5312087450269404d-2*t10*
     1         t3*wght/t6
            Cmat(iq,D1_GAB) = Cmat(iq,D1_GAB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.5312087450269404d-2*t11*
     1         t3*wght/t9
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020562d-2*gammaaa*t2/(rhoa**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhoa**1.3333333333333333d+0
            t6 = 1/t3**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t4*t5*wght
            Amat(iq,D1_RA) = 4.0832233200718415d-2*gammaaa*t2*t6*wght/rh
     1         oa**2.3333333333333334d+0-1.2407009817987999d+0*rhoa**3.3
     2         33333333333333d-1*t1*t4*wght+Amat(iq,D1_RA)
            Cmat(iq,D1_GAA) = Cmat(iq,D1_GAA)-1.5312087450269404d-2*t2*t
     1         6*wght/t5
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            t1 = param(1)
            t2 = param(2)
            t3 = 1.6455307846020562d-2*gammabb*t2/(rhob**2.6666666666666
     1         666d+0*t1)+1.0d+0
            t4 = 1.0d+0-1.0d+0/t3
            t5 = rhob**1.3333333333333333d+0
            t6 = 1/t3**2
            fnc(iq) = fnc(iq)-9.305257363491d-1*t1*t4*t5*wght
            Amat(iq,D1_RB) = 4.0832233200718415d-2*gammabb*t2*t6*wght/rh
     1         ob**2.3333333333333334d+0-1.2407009817987999d+0*rhob**3.3
     2         33333333333333d-1*t1*t4*wght+Amat(iq,D1_RB)
            Cmat(iq,D1_GBB) = Cmat(iq,D1_GBB)-1.5312087450269404d-2*t2*t
     1         6*wght/t5
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\right)-0.9305257363491\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77, 3865 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
      subroutine nwxcm_x_pbe_d2(param,tol_rho,ipol,nq,wght, 
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
C> \brief Evaluate the nwxcm_x_pbe functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   {\it t_1} &=& {\it param}\left(1\right)\\\\ 
C>   {\it t_2} &=& {{1}\over{{\it t_1}}}\\\\ 
C>   {\it t_3} &=& {\it param}\left(2\right)\\\\ 
C>   f &=& -0.9305257363491\,{\it t_1}\,\rho_\beta^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\beta\beta}}\over{\rho_\beta^{{{8}\over{3}}}}}+1.0}}\right)-0.9305257363491\,{\it t_1}\,\rho_\alpha^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{\alpha\alpha}}\over{\rho_\alpha^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,{\it t_1}\,\rho_s^{{{4}\over{3}}}\,\left(1.0-{{1.0}\over{{{0.01645530784602056\,{\it t_2}\,{\it t_3}\,\sigma_{ss}}\over{\rho_s^{{{8}\over{3}}}}}+1.0}}\right)\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JP Perdew, K Burke, M Ernzerhof, Phys.Rev.Lett. 77, 3865 (1996)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRevLett.77.3865 "> 
C> 10.1103/PhysRevLett.77.3865 </a> 
C> 
      subroutine nwxcm_x_pbe_d3(param,tol_rho,ipol,nq,wght, 
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
