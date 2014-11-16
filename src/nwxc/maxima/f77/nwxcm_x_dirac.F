C> \ingroup nwxc 
C> @{ 
C> 
C> \file nwxcm_x_dirac.F 
C> The nwxcm_x_dirac functional 
C> 
C> @} 
C> 
C> \ingroup nwxc_priv 
C> @{ 
C> 
C> \brief Evaluate the nwxcm_x_dirac functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JC Slater, Phys.Rev. 81, 385 (1951)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRev.81.385 "> 
C> 10.1103/PhysRev.81.385 </a> 
C> 
      subroutine nwxcm_x_dirac(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
            fnc(iq) = fnc(iq)-1.8610514726982d+0*rhoa**1.333333333333333
     1         3d+0*wght
            Amat(iq,D1_RA) = Amat(iq,D1_RA)-2.4814019635975998d+0*rhoa**
     1         3.333333333333333d-1*wght
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = (-9.305257363491d-1*rhob**1.3333333333333333d+0-9.
     1         305257363491d-1*rhoa**1.3333333333333333d+0)*wght+fnc(iq)
            Amat(iq,D1_RA) = Amat(iq,D1_RA)-1.2407009817987999d+0*rhoa**
     1         3.333333333333333d-1*wght
            Amat(iq,D1_RB) = Amat(iq,D1_RB)-1.2407009817987999d+0*rhob**
     1         3.333333333333333d-1*wght
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
            fnc(iq) = fnc(iq)-9.305257363491d-1*rhoa**1.3333333333333333
     1         d+0*wght
            Amat(iq,D1_RA) = Amat(iq,D1_RA)-1.2407009817987999d+0*rhoa**
     1         3.333333333333333d-1*wght
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
            fnc(iq) = fnc(iq)-9.305257363491d-1*rhob**1.3333333333333333
     1         d+0*wght
            Amat(iq,D1_RB) = Amat(iq,D1_RB)-1.2407009817987999d+0*rhob**
     1         3.333333333333333d-1*wght
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_dirac functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JC Slater, Phys.Rev. 81, 385 (1951)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRev.81.385 "> 
C> 10.1103/PhysRev.81.385 </a> 
C> 
      subroutine nwxcm_x_dirac_d2(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> 
C> \brief Evaluate the nwxcm_x_dirac functional [1] 
C> 
C> \f{eqnarray*}{ 
C>   f &=& -0.9305257363491\,\rho_\beta^{{{4}\over{3}}}-0.9305257363491\,\rho_\alpha^{{{4}\over{3}}}\\\\ 
C>   g &=& 0\\\\ 
C>   G &=& -0.9305257363491\,\rho_s^{{{4}\over{3}}}\\\\ 
C> \f} 
C> 
C> ### References ### 
C> 
C> [1] JC Slater, Phys.Rev. 81, 385 (1951)  , DOI: 
C> <a href="http://dx.doi.org/10.1103/PhysRev.81.385 "> 
C> 10.1103/PhysRev.81.385 </a> 
C> 
      subroutine nwxcm_x_dirac_d3(param,tol_rho,ipol,nq,wght, 
     +rho,fnc,Amat,Amat2,Amat3) 
c $Id: $ 
      implicit double precision (a-h,o-z), integer (i-n) 
#include "nwxc_param.fh" 
      double precision param(*)     !< [Input] Parameters of functional 
      double precision tol_rho      !< [Input] The lower limit on the density 
      integer ipol                  !< [Input] The number of spin channels 
      integer nq                    !< [Input] The number of points 
      double precision wght         !< [Input] The weight of the functional 
      double precision rho(nq,*)    !< [Input] The density 
      double precision fnc(nq)      !< [Output] The value of the functional 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho 
c 
c     Sampling Matrices for the XC Kernel 
c 
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho 
c 
      integer iq 
      double precision tmp 
      double precision rhoa,rhob 
      double precision gammaaa,gammaab,gammabb 
      double precision taua,taub 
      do iq = 1, nq 
        if (ipol.eq.1) then 
          rhoa    = 0.5d0*rho(iq,R_T) 
          if (rhoa.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho 
        else  ! ipol.eq.1 
          rhoa    = rho(iq,R_A) 
          rhob    = rho(iq,R_B) 
          if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then 
          elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then 
          elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then 
          endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho 
        endif ! ipol.eq.1 
      enddo ! iq 
      end 
C> @} 
