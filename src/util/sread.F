      subroutine sread(itape,a,n)
#include "errquit.fh"
*
* $Id: sread.F,v 1.1 2003-10-17 22:00:10 carlfahl Exp $
*
      double precision a(n)
      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call ssread(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('sread: left .ne. 0',left, INPUT_ERR)
c
      end
      subroutine ssread(itape,a,n)
      double precision a(n)
c
      read(itape,end=10,err=20) a
      return
 10   call errquit('ssread: end of file reading from unit ', itape,
     &       DISK_ERR)
 20   call errquit('ssread: error reading from unit ', itape,
     &       DISK_ERR)
c
      end
