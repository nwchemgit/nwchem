      subroutine pre_wrtnoe(lfndrs,fildrs,mdrs,ndrs,idrs,rdrs)
c
      implicit none
c
      integer lfndrs
      character*(*) fildrs
      integer mdrs,ndrs
      integer idrs(2,mdrs)
      real*8 rdrs(5,mdrs)
c
      integer i,j
c
c     write solute NOE restraints
c
      if(ndrs.gt.0) then
      open(unit=lfndrs,file=fildrs(1:index(fildrs,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(lfndrs,1000) ndrs
 1000 format(i5)
      do 1 j=1,ndrs
      write(lfndrs,1001)
     + (idrs(i,j),i=1,2),(rdrs(i,j),i=1,5)
 1001 format(2i10,3f12.6,/,20x,2e12.5)
    1 continue
      close(unit=lfndrs)
      endif
c
      return
c
 9999 continue
      call errquit('Could not open fildrs',0)
      return
      end
