      subroutine mprint(rmat,idim,jdim)
C$Id: mprint.f,v 1.2 1995-02-02 23:23:11 d3g681 Exp $
      implicit real*8 (a-h,o-z)
      dimension rmat(idim,jdim)
      do 100 i=1,idim
         write(*,10) (rmat(i,j), j=1,jdim)
100   continue
10    format(19x,4(f10.6))
      return
      end
