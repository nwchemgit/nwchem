      integer function tce_diis_switch(iter)
c
c $Id: tce_diis_switch.F,v 1.1 2008-09-30 18:35:43 jhammond Exp $
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
