      integer function selci_isum(n,m,im)
*
* $Id: isum.f,v 1.2 1997-10-31 23:42:07 d3e129 Exp $
*
      dimension m(im,*)
c
c     return sum of integer array
c
      is = 0
      do 10 i = 1,n
         is = is + m(1,i)
 10   continue
c
      selci_isum = is
c
      end
