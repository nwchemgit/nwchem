C
C $Id$
C
      subroutine hess_tidy()
c
c  This routine cleans up any of the memory that was used and any
c  other details that need to be taken care of.
c
      implicit none
#include "errquit.fh"
c
#include "dra.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
c
c  get rid of geometry and basis set handles and associated memory
c
      if (.not. geom_destroy(geom))
     *  call errquit('hess_tidy: geom_destroy problem', 555, GEOM_ERR)
      if (.not. bas_destroy(basis))
     *  call errquit('hess_tidy: bas_destroy problem', 555, BASIS_ERR)
c
c  get rid of local memory for the nuclear hessian Exy
c
      if (.not.ma_pop_stack(l_act))
     *   call errquit('hess_tidy: cannot deallocate oactive',555,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_dder))
     *   call errquit('hess_tidy: cannot deallocate hess_dder',555,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_hess_xc))
     *   call errquit('hess_tidy: cannot deallocate hess_xc',555,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_exy))
     *   call errquit('hess_tidy: cannot deallocate Exy',555,
     &       MA_ERR)
c
c  Get rid of temporary DRAs and terminate DRA interface
c
c      if (dra_delete(fx_dra_handle).ne.0) call errquit
c     *   ('hess_tidy: unable to delete the fx DRA',555, UNKNOWN_ERR)
c      if (dra_delete(sx_dra_handle).ne.0) call errquit
c     *   ('hess_tidy: unable to delete the sx DRA',555, UNKNOWN_ERR)
c      if (dra_terminate().ne.0) call errquit
c     *   ('hess_tidy: unable to terminate DRA interface',555,
c     &       UNKNOWN_ERR)
c
c
      return
      end
c
      subroutine hess_energytidy(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c  This routine sets the rtdb up the way it was before entering the
c  Hessian code.
c
      integer rtdb
c
      if (.not.geom_rtdb_store(rtdb,geomold,'geometry'))
     $  call errquit('hess_energytidy: geom_rtdb_store failed',555,
     &       RTDB_ERR)
      if (.not.geom_destroy(geomold))
     $  call errquit('hess_energytidy: geom_destroy failed',555,
     &       GEOM_ERR)
      if (theory.eq.'scf') then
         if (.not. rtdb_put(rtdb, 'scf:thresh', MT_DBL, 1, tthresh))
     $     call errquit('hess_energytidy: rtdb_put failed',tthresh,
     &       RTDB_ERR)
      else
         if (.not. rtdb_put(rtdb, 'dft:g_conv', MT_DBL, 1, tthresh))
     $     call errquit('hess_energytidy: rtdb_put failed',tthresh,
     &       RTDB_ERR)
      endif
c
      return
      end

