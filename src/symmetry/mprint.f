      subroutine mprint(rmat,idim,jdim)
      implicit real*8 (a-h,o-z)
      dimension rmat(idim,jdim)
      do 100 i=1,idim
         write(*,10) (rmat(i,j), j=1,jdim)
100   continue
10    format(19x,4(f10.6))
      return
      end
