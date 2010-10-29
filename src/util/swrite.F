      subroutine swrite(itape,a,n)
#include "errquit.fh"
*
* $Id$
*
      double precision a(n)
      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call sswrit(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('swrite: left .ne. 0',left, DISK_ERR)
c
      end
      subroutine sswrit(itape,a,n)
      double precision a(n)
c
      write(itape,err=20) a
      return
 20   call errquit('sswrit: error writing to unit ', itape, DISK_ERR)
c
      end
