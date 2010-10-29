      integer function tce_diis_switch(iter)
c
c $Id$
c
c     Returns DIIS order
c
      implicit none
#include "tce_diis.fh"
c      integer iter
c
      if (diis .le. 0) then
        tce_diis_switch = 0
        return
      endif
      if (mod(iter,2) .ne. 0) then
        tce_diis_switch = 0
        return
      endif
      if (iter .lt. 20) then
        tce_diis_switch = min(diis,iter)
      else
        tce_diis_switch = diis
      endif
      return
      end
