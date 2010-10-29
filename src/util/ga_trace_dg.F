      double precision function ga_trace_diag(g_a)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
c
c     Take the trace of diagonal elements of the
c     real global square array  in a 'scalable' fashion
c
c     Everyone ends up with the trace.
c
      integer ma_type, dim1, dim2, n
      integer i, ilo, ihi, jlo, jhi
      double precision xx,tr
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a, GA_ERR)
      n = dim1
c
c
c     Extract and sum the diags local to each process
c
      tr = 0.d0
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
         do i = ilo,ihi
            call ga_get(g_a, i, i, i, i, xx, 1)
            tr = tr + xx
         enddo
      endif
c
      call ga_sync()
      call ga_dgop(Msg_trace_diag, tr, 1, '+')
      call ga_sync()
      ga_trace_diag = tr
c
      end
