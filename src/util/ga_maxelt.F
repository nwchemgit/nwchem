       subroutine ga_maxelt(g_a, value)
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
c
c     Determine the absolute value of the element with largest 
c     absolute magnitude
c
       integer ilo, ihi, jlo, jhi, idim, jdim, ld
       integer j, ind, ind_max
       integer idamax
       external idamax
c
c     Each process first determines max of data local to it
c
       call ga_sync()
       value = 0.d0
       call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
       if (ihi.gt.0 .and. ihi.ge.ilo) then
          idim = ihi - ilo + 1
          jdim = jhi - jlo + 1
          call ga_access(g_a, ilo, ihi, jlo, jhi, ind, ld)
          do j = 0, jdim-1
             ind_max = idamax(idim, dbl_mb(ind), 1) + ind - 1
             value = max(value, abs(dbl_mb(ind_max)))
             ind = ind + ld
          enddo
          call ga_release(g_a, ilo, ihi, jlo, jhi)
       endif
c
c     Finally a global max via message passing
c
       call ga_sync()
       call ga_dgop(Msg_MxElt, value, 1, 'absmax')
       call ga_sync()
c
       end
