      subroutine hess_wrt
*     $Id: hess_wrt.F,v 1.1 2000-05-31 21:49:34 edo Exp $
      implicit none
#include "mafdecls.fh"
#include "hess_info.fh"
c
      character*255 filehess,fileddipole
      integer nat
      integer h_ddipole, k_ddipole
      logical dipole_okay
c
      if (.not. MA_Push_Get(MT_DBL,3*n3xyz,'derivative dipole moment',
     &     h_ddipole, k_ddipole)) call errquit(
     &     'hess_wrt: ma_push_get failed',15)
      dipole_okay=.false.
      call util_file_name('dipole',  .false., .false.,fileddipole)
      call util_file_name('hess',  .false., .false.,filehess)
      call stpr_wrt_fd_from_sq(dbl_mb(k_exy),n3xyz,filehess)
      if (dipole_okay) then
         call stpr_wrt_fd_dipole(dbl_mb(k_ddipole),nat,fileddipole)
      else
         call util_file_unlink(fileddipole)
      endif
      if(.not.ma_pop_stack(h_ddipole)) call errquit
     &     ('hess_wrt: ma_pop_stack(h_ddipole) failed?',43)
      return
      end
