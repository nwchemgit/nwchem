      subroutine selci_icopy(n,x,ix,y,iy)
      implicit integer (a-z)
*
* $Id: selci_icopy.f,v 1.1 2003-04-07 21:58:55 windus Exp $
*
      dimension x(*),y(*)
c
      ixx = 1
      iyy = 1
      do 10 i = 1,n
         y(iyy) = x(ixx)
         ixx = ixx + ix
         iyy = iyy + iy
 10   continue
c
      return
      end
