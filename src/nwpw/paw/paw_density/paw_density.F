!$Id: paw_density.F,v 1.1 2003-02-16 21:40:00 marat Exp $
!
!$Log: not supported by cvs2svn $
**********************************
*
*
*
*
*
**********************************
      subroutine init_paw_density()
      implicit none

#include "paw_density_data.fh"
#include "paw_basis.fh"
#include "mafdecls.fh"      
#include "paw_ma.fh"
#include "paw_geom.fh"

      logical ok
      integer in,ia
      integer nion
      integer paw_rho_size

      nion = ion_nion()
*     *** lmax for multipole expansion of the density***      
      paw_rho_lmax = 0
      
*     *** size for multipole density arrays ***
      ok = my_alloc(mt_int,nion,"i_paw_rho",
     >                                i_paw_rho)

      paw_rho_size = 0     
      do in=1,nion
        ia = ion_katm(in)
        int_mb(i_paw_rho(1) + in - 1) = paw_rho_size
        write(*,*) "paw_density:",paw_rho_size
        paw_rho_size = paw_rho_size + paw_basis_i_cut(ia)*
     >                                (paw_rho_lmax+1)**2
      end do

*     *** allocate multipole density arrays ****
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ae",
     >                               paw_rho_u_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ae",
     >                               paw_rho_d_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ps",
     >                               paw_rho_u_ps) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ps",
     >                               paw_rho_d_ps) 




      end 
**********************************
*
*
*
*
*
**********************************
      subroutine paw_density_solve()
      implicit none

      end
