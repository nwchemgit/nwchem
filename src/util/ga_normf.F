       subroutine GA_NormF(g_a, value)
C$Id: ga_normf.F,v 1.1 1998-02-02 15:57:40 bernhold Exp $
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
c
C     Frobenius norm of a global array
c
       integer ilo, ihi, jlo, jhi, idim, jdim, ld
       integer ind
       Double Precision Work
       Double Precision dlange
       External dlange
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
          Value = dlange('F', IDim, JDim, Dbl_MB(Ind), LD, Work)
          call ga_release(g_a, ilo, ihi, jlo, jhi)
       endif
c
c     Combine partial results
C     NOTE:  dlange is careful about overflows, but what's below is not!
c
       Value = Value * Value
       call ga_sync()
       call ga_dgop(Msg_NormF, value, 1, '+')
       call ga_sync()
       Value = Sqrt( Value )
c
       Return
       end
