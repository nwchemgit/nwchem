      subroutine hess_wrt
*     $Id: hess_wrt.F,v 1.2 2002-02-06 01:11:16 windus Exp $
      implicit none
#include "mafdecls.fh"
#include "hess_info.fh"
c
      character*255 filehess,fileddipole
      logical dipole_okay
c
      dipole_okay=.true.
      call util_file_name('fd_ddipole', .false., .false.,fileddipole)
      call util_file_name('hess', .false., .false.,filehess)
      call stpr_wrt_fd_from_sq(dbl_mb(k_exy),n3xyz,filehess)
      if (dipole_okay) then
         call stpr_wrt_fd_dipole(dbl_mb(k_dder),ncent,fileddipole)
      else
         call util_file_unlink(fileddipole)
      endif
      return
      end
