      subroutine ga_get_diagonal(g_a, diags)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      double precision diags(*)
c
c     Extract out the diagonal elements of the
c     real global array  in a 'scalable' fashion
c
c     Everyone ends up with the diagonal.
c
      integer ma_type, dim1, dim2, n
      integer i, ilo, ihi, jlo, jhi
c
      call ga_sync
c
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_get_diagonal: array is not square/real', g_a, GA_ERR)
      n = dim1
c
      call dfill(n, 0.0d0, diags, 1)
c
c     Extract out the diags local to each process
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
         do i = ilo,ihi
            call ga_get(g_a, i, i, i, i, diags(i), 1)
         enddo
      endif
c
      call ga_dgop(333, diags, n, '+')
c
      end
