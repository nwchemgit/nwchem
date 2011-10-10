      subroutine argos_space_gthdl(ibbl,lhandl,lself,
     + iwfr,iwto,isfr,isto,
     + lbbl,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2)
      logical lself
c
      integer ndx,i,ibox
      character*255 string
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
      lself=ibox.eq.lbbl(ibbl,2)
c
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1001) ibox,me
 1001 format('SP0001: Could not find box',i5,' on node',i5)
      call md_abort(string,0)
      endif
c
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
c
      return
      end
c $Id$
