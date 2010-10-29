*
* $Id$
*
      

      double precision function simp(n,y,h)
      implicit none
      integer n
      double precision y(n)
      double precision h,s
      integer ne,no
      double precision dsum
      external         dsum

      ne=n/2
      no=ne+1
      S=2.0d0*dsum(no,y(1),2) + 4.0d0*dsum(ne,y(2),2)-y(1)-y(n)
      simp=s*h/3.0d0
      return
      end

