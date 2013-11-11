C> \ingroup selci
C> @{
      integer function selci_isum(n,m,im)
*
* $Id$
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
C> @}
