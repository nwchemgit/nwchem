      integer function idmin(n,a,ia)
*
* $Id: idmin.f,v 1.2 1997-10-31 23:42:04 d3e129 Exp $
*
      real *8 a(ia,*)
c
c     return index of minimum value in array a
c
      val = a(1,1)
      ind = 1
c
      do 10 i = 1,n
         if (a(1,i).lt.val) then
            val = a(1,i)
            ind = i
         endif
10    continue
c
      idmin = ind
      end
