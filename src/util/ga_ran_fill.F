      subroutine ga_ran_fill(g_a, ilo, ihi, jlo, jhi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
       integer g_a
       integer ilo, ihi, jlo, jhi
       integer i, j, bufsiz
       parameter(bufsiz = 1024)
       double precision buf(bufsiz)
       integer ibot, itop, icount, me, nproc
       double precision drand
       external drand
c
c     Fill patch of global array with random numbers
c
       call ga_sync()
c
       icount = 0
       me = ga_nodeid()
       nproc = ga_nnodes()
c
       do j = jlo, jhi
          do ibot = ilo, ihi, bufsiz
             if (mod(icount,nproc) .eq.  me) then
                itop = min(ihi, ibot+bufsiz-1)
                do i = 1, itop-ibot+1
                   buf(i) = drand(0)
                enddo
                call ga_put(g_a, ibot, itop, j, j, buf, bufsiz)
             endif
             icount = icount + 1
          enddo
       enddo
c
       call ga_sync()
c
       end
