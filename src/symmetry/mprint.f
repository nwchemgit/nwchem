      subroutine mprint(rmat,idim,jdim)
C$Id: mprint.f,v 1.3 1995-12-15 12:13:03 d3g681 Exp $
      implicit none
      integer idim, jdim
      double precision rmat(idim,jdim)
      integer i, j
c
      do 100 i=1,idim
         write(*,10) (rmat(i,j), j=1,jdim)
100   continue
10    format(19x,4(f10.6))
c
      end
