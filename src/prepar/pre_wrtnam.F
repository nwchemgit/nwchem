      logical function pre_wrtnam(lfntmp,filtmp,catm,matm,natm)
c
      implicit none
c
      integer lfntmp,matm,natm
      character*255 filtmp
      character*6 catm(3,matm)
c
      integer i
c
      open(unit=lfntmp,file=filtmp(1:index(filtmp,' ')-1),
     + form='formatted',status='new',err=9999)
c
      do 1 i=1,natm
      write(unit=lfntmp,1000) catm(1,i),catm(2,i)
 1000 format(2a6)
    1 continue
c
      close(unit=lfntmp)
c
      pre_wrtnam=.true.
      return
c
 9999 continue
      pre_wrtnam=.false.
      return
      end
