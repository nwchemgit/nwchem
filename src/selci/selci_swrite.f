      subroutine selci_swrite(itape,a,n)
*
* $Id: selci_swrite.f,v 1.2 2003-10-27 23:21:39 marat Exp $
*
      real*8 a(n)
      parameter (lenbuf = 512)
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call selci_sswrit(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('swrite: left .ne. 0',left,0)
c
      end
      subroutine selci_sswrit(itape,a,n)
      real*8 a(n)
c
      write(itape) a
c
      end
