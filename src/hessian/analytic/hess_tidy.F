C
C $Id: hess_tidy.F,v 1.4 2000-06-01 23:02:10 windus Exp $
C
      subroutine hess_tidy()
c
c  This routine cleans up any of the memory that was used and any
c  other details that need to be taken care of.
c
      implicit none
c
#include "dra.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
c
      integer j
c
c  get rid of geometry and basis set handles and associated memory
c
      if (.not. geom_destroy(geom))
     *  call errquit('hess_tidy: geom_destroy problem', 555)
      if (.not. bas_destroy(basis))
     *  call errquit('hess_tidy: bas_destroy problem', 555)
c
c  get rid of local memory for the nuclear hessian Exy
c
      if (.not.ma_pop_stack(l_act))
     *   call errquit('hess_tidy: cannot deallocate oactive',555)
      if (.not.ma_pop_stack(l_hess_xc))
     *   call errquit('hess_tidy: cannot deallocate hess_xc',555)
      if (.not.ma_pop_stack(l_exy))
     *   call errquit('hess_tidy: cannot deallocate Exy',555)
c
c  Release the density and global fock derivative matrices (Qs or RHS).
c
      if (scftype .eq. 'UHF') then
        if (.not. ga_destroy(g_dens(1))) call errquit
     *     ('hess_tidy: failed to free g_dens(1)', 555)
        if (.not. ga_destroy(g_dens(2))) call errquit
     *     ('hess_tidy: failed to free g_dens(2)', 555)
      else
        if (.not. ga_destroy(g_dens(1))) call errquit
     *     ('hess_tidy: failed to free g_dens(1)', 555)
      endif

      do j = 1, ndens*n3xyz
        if (.not.ga_destroy(g_rhs_xc(j))) 
     *    call errquit('hess_tidy: problem destroying g_rhs_xc',1)
      enddo
c
c  Get rid of temporary DRAs and terminate DRA interface
c
      if (dra_delete(sx_dra_handle).ne.0) call errquit
     *   ('hess_tidy: unable to delete the sx DRA',555)
      if (dra_terminate().ne.0) call errquit
     *   ('hess_tidy: unable to terminate DRA interface',555)
c
      return
      end
