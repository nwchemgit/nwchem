      integer function selci_idmin(n,a,ia)
*
* $Id: idmin.f,v 1.3 2009-02-10 03:27:12 jhammond Exp $
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
      selci_idmin = ind
      end
