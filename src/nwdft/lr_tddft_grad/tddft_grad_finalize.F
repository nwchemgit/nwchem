c
c     Finalize TDDFT gradient calculation
c
      subroutine tddft_grad_finalize(nsroots,ihdl_sroots,
     +                           ntroots,ihdl_troots,
     +                           nroots, ihdl_roots,
     +                           ihdl_geom, ihdl_bfao)
c
      implicit none
c
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
c
      integer nsroots     ! the number of singlet roots
      integer ihdl_sroots ! the memory handle for the singlet roots
      integer ntroots     ! the number of triplet roots
      integer ihdl_troots ! the memory handle for the triplet roots
      integer nroots      ! the number of roots
      integer ihdl_roots  ! the memory handle for the roots
c
      integer ihdl_geom   ! the geometry handle
      integer ihdl_bfao   ! the basis set handle
c
      character*32 pname
      pname = "tddft_grad_finalize: "
c
      if (nsroots.ne.0) then
        if (.not.ma_pop_stack(ihdl_sroots)) 
     +    call errquit(pname//'failed to deallocate sroots')
      endif
      if (ntroots.ne.0) then
        if (.not.ma_pop_stack(ihdl_troots)) 
     +    call errquit(pname//'failed to deallocate troots')
      endif
      if (nroots.ne.0) then
        if (.not.ma_pop_stack(ihdl_roots)) 
     +    call errquit(pname//'failed to deallocate roots')
      endif
c
      if (.not.geom_destroy(ihdl_geom)) 
     +    call errquit(pname//'failed to destroy geometry')
      if (.not.bas_destroy(ihdl_bfao)) 
     +    call errquit(pname//'failed to destroy basis set')
c
      end
c $Id$
