!$Id: paw_xc.F,v 1.2 2003-02-17 05:43:37 marat Exp $
!
!$Log: not supported by cvs2svn $
!Revision 1.1  2003/02/17 05:28:19  marat
!initial revision
!MV
!

**********************************
*
*
*
*
*
**********************************
      subroutine init_paw_xc(lmax)
      implicit none

      integer lmax

#include "paw_xc_data.fh"
#include "paw_basis.fh"
#include "mafdecls.fh"      
#include "paw_ma.fh"
#include "paw_geom.fh"

      logical ok
      integer in,ia
      integer nion
      integer paw_vxc_size

      nion = ion_nion()

*     *** lmax for multipole expansion of the density***      
      paw_xc_lmax = lmax
      
*     *** size for multipole density arrays ***
      ok = my_alloc(mt_int,nion,"i_paw_vxc",
     >                                i_paw_vxc)

      paw_vxc_size = 0     
      do in=1,nion
        ia = ion_katm(in)
        int_mb(i_paw_vxc(1) + in - 1) = paw_vxc_size
        paw_vxc_size = paw_vxc_size + paw_basis_i_cut(ia)*
     >                                (paw_xc_lmax+1)**2
      end do

*     *** allocate multipole density arrays ****
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_u_ae",
     >                               paw_vxc_u_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_d_ae",
     >                               paw_vxc_d_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_u_ps",
     >                               paw_vxc_u_ps) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_d_ps",
     >                               paw_vxc_d_ps) 


      write(*,*) "paw_xc initialized"


      end 

**********************************
*
*
*
*
*
**********************************
      subroutine paw_xc_end()
      implicit none

#include "paw_xc_data.fh"
#include "paw_ma.fh"

      logical ok
      
      ok = my_dealloc(paw_vxc_d_ps) 
      ok = ok .and. my_dealloc(paw_vxc_u_ps) 
      ok = ok .and. my_dealloc(paw_vxc_d_ae) 
      ok = ok .and. my_dealloc(paw_vxc_u_ae) 
      ok = ok .and. my_dealloc(i_paw_vxc)

      if (.not.ok)
     > call errquit("paw_vxc_end: error freeing heap",0)


      end 
