      integer function ishft(i,j)
      implicit none
c
c     routine to return i shifted by j bits
c
      if(j.eq.0) ishft=i
      if(j.gt.0) ishft=lshift(i,j)
      if(j.lt.0) ishift=rshift(i,-j)
      return
      end
