      double precision function f_xs_new(n, c_n, rho, gradrho,
     A     mylda, mu)
      implicit none
      integer n
      double precision c_n,rho,gradrho
      double precision mu
      external mylda
      double precision e,eprime
#include "xc_hcth_extra.fh"
      call mylda(mu, rho,e,eprime)
      f_xs_new=e*c_n*
     *     ((gamma*gradrho*gradrho)/dg(gamma,rho,gradrho))**n
      return
      end
      double precision function df_xs_by_drhos_new(n, c_n, rho, gradrho,
     A     mylda, mu)
      implicit none
      integer n
      double precision mu
      double precision c_n,rho,gradrho
      external mylda
      double precision e,eprime
#include "xc_hcth_extra.fh"      
c
      call mylda(mu, rho,e,eprime)
      if(n.eq.0) then
         df_xs_by_drhos_new=eprime*c_n
      else
         df_xs_by_drhos_new=(eprime -
     -        e*n* 
     *        (8d0/3d0)*rho**(5d0/3d0)/
     /        dg(gamma,rho,gradrho))*
     /        c_n*
     *        ((gamma*gradrho*gradrho)/dg(gamma,rho,gradrho))**n
      endif
      return
      end
      double precision function df_xs_by_dzs_new(n, c_n, rho, gradrho,
     A     mylda, mu)
      implicit none
      integer n
      double precision c_n,rho,gradrho
      double precision mu
      double precision e,eprime
#include "xc_hcth_extra.fh"      
      call mylda(mu, rho,e,eprime)
      if(n.eq.0) then
         df_xs_by_dzs_new=0d0
      else
         df_xs_by_dzs_new=e*c_n*n*
     *        (gamma*gradrho*gradrho/dg(gamma,rho,gradrho))**(n-1)*
     *        2d0*gamma*gradrho*rho**(8d0/3d0)/
     /        (dg(gamma,rho,gradrho)**2)
      endif
      return
      end
c ldas      
      subroutine xc_kernel_lda(in1, rho,e,eprime)
      implicit none
      double precision in1      ! ignore [in]
      double precision rho      ! [in]
      double precision e        ! [out]
      double precision eprime   ! [out]
c
      double precision pi
      double precision gradrho
#include "xc_hcth_extra.fh"
      pi=acos(-1d0)
      e=-3d0/2d0*(3d0/(4d0*pi))**(1d0/3d0)*rho**(4d0/3d0)
      eprime=(4d0/3d0)*e/rho
      return
      end
      subroutine xc_kernel_lrlda(mu,rho,e,eprime)
      implicit none
      double precision rho      ! [in]
      double precision e        ! [out]
      double precision eprime   ! [out]
      double precision mu       ! omega [in]
c
      call LRCLSDA(mu, rho, e, eprime)
      return
      end
