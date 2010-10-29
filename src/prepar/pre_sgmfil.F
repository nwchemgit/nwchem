      logical function pre_sgmfil(sgmnam,lfnsgm,filsgm)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfnsgm
      character*10 sgmnam
      character*255 filsgm
c
      integer lens,lend,lenf,i
c
      lens=index(sgmnam,' ')-1
      if(lens.le.0) lens=10
c
      do 1 i=mdirpar,1,-1
      if(nfilpar(i).gt.0) then
      filsgm=
     + dirpar(i)(1:index(dirpar(i),' ')-1)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c     
      if(util_print('where',print_debug)) then
      write(*,1001) filsgm(1:lenf)
 1001 format(' Checking existence of file ',a)
      endif
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=1)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
      endif
    1 continue
c
      pre_sgmfil=.false.
      return
      end

