       subroutine ga_screen(g_a, value)
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
c
c    
c     set all elements whose absolute value is less than
c     value to a hard zero
c
       integer ilo, ihi, jlo, jhi, idim, jdim, ld
       integer i, j, ij, ind
       double precision zero
       parameter (zero = 0.0d0)
c
       call ga_sync()
       call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
       if (ihi.gt.0 .and. jhi.gt.0) then
          idim = ihi - ilo + 1
          jdim = jhi - jlo + 1
          call ga_access(g_a, ilo, ihi, jlo, jhi, ind, ld)
          do j = 0, jdim-1
             do i = 0, idim-1
                ij = i + ind
                if (abs(dbl_mb(ij)) .lt. value) dbl_mb(ij) = zero
             enddo
             ind = ind + ld
          enddo
          call ga_release_update(g_a, ilo, ihi, jlo, jhi)
       endif
c
       call ga_sync()
c
       end
