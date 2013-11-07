C> \ingroup selci
C> @{
C>
      subroutine selci_sread(itape,a,n)
*
* $Id$
*
      real*8 a(n)
      parameter (lenbuf = 512)
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call sselci_sread(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('swrite: left .ne. 0',left,0)
c
      end
      subroutine sselci_sread(itape,a,n)
      real*8 a(n)
c
      read(itape) a
c
      end
C>
C> @}
