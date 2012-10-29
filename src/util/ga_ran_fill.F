      subroutine ga_ran_fill(g_a, ilo, ihi, jlo, jhi)
C$Id$
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
       integer g_a
       integer ilo, ihi, jlo, jhi
       integer itypa
       integer i, j, bufsiz
       parameter(bufsiz = 1024)
       double precision buf(bufsiz)
       double complex zbuf(bufsiz)
       integer ibot, itop, icount, me, nproc
       double precision util_random
       external util_random
c
c     Fill patch of global array with random numbers
c
       call ga_sync()
       call ga_inquire(g_a,itypa,i,j)
c
       icount = 0
       me = ga_nodeid()
       nproc = ga_nnodes()
c
       do j = jlo, jhi
          do ibot = ilo, ihi, bufsiz
             if (mod(icount,nproc) .eq.  me) then
                itop = min(ihi, ibot+bufsiz-1)
                if (itypa.eq.MT_DBL) then
                  do i = 1, itop-ibot+1
                     buf(i) = util_random(0)
                  enddo
                  call ga_put(g_a, ibot, itop, j, j, buf, bufsiz)
                else if (itypa.eq.MT_DCPL) then
                  do i = 1, itop-ibot+1
                     zbuf(i) = cmplx(util_random(0),util_random(0))
                  enddo
                  call ga_put(g_a, ibot, itop, j, j, zbuf, bufsiz)
                else
                  call errquit('ga_ran_fill: illegal type',itypa,UERR)
                endif
             endif
             icount = icount + 1
          enddo
       enddo
c
       call ga_sync()
c
       end
