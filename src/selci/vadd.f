      subroutine selci_vadd(n,a,ia,b,ib,c,ic)
*
* $Id: vadd.f,v 1.2 1997-10-31 23:42:32 d3e129 Exp $
*
      real*8 a(ia,*),b(ib,*),c(ic,*)
c
c     c(*) = b(*) + a(*)
C
      do 10 m = 1,n
         c(1,m) = b(1,m) + a(1,m)
 10   continue
c
      end
