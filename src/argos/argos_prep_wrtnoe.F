      subroutine argos_prep_wrtnoe(lfntop,filtop,mdrs,ndrs,idrs,rdrs)
c
      implicit none
c
      integer lfntop
      character*(*) filtop
      integer mdrs,ndrs
      integer idrs(2,mdrs)
      real*8 rdrs(5,mdrs)
c
      integer i,j
      character*3 string
c
c     write solute NOE restraints
c
      if(ndrs.gt.0) then
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='unknown',err=9999)
    9 continue
      read(lfntop,100,end=10,err=9999) string
  100 format(a3)
      if(string.ne.'noe') goto 9
      goto 11
   10 continue
      write(lfntop,101)
  101 format('noe')
   11 continue
      write(lfntop,1000) ndrs
 1000 format(i5)
      do 1 j=1,ndrs
      write(lfntop,1001)
     + (idrs(i,j),i=1,2),(rdrs(i,j),i=1,5)
 1001 format(2i10,3f12.6,/,20x,2e12.5)
    1 continue
      close(unit=lfntop)
      endif
c
      return
c
 9999 continue
      call md_abort('Could not open filtop to write noe',0)
      return
      end
c $Id$
