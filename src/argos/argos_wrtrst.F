      subroutine argos_wrtrst(lfn,fil,lveloc)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 timer_wall
      external timer_wall
c
      integer lfn
      character*255 fil
      logical lveloc
      integer i,left
      character*255 filn
c
      if(me.eq.0) then
      if(keepr.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown',err=9999)
      else
      if(ntype.eq.3.and.npg.le.1) then
      print*,npg,ntype,irun
      print*,fil(1:index(fil,' ')-1)
      print*,fil(1:index(fil,'.rst')-1)
      write(filn,'(a,i5.5,a)') fil(1:index(fil,'.rst')-1),
     + irun,'.rst '
      else
      write(filn,'(a,a,i5.5,a)') fil(1:index(fil,'.rst')-1),'-',
     + keepr,'.rst '
      endif
      open(unit=lfn,file=filn(1:index(filn,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      endif
c
      call argos_space_wrtrst(lfn,fil,lveloc,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_gw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),dbl_mb(i_xscr),projct,int_mb(i_lseq))
c
      call argos_wtrest(lfn)
c
      call argos_prop_wtrest(lfn)
c
      call argos_space_wtrest(lfn)
c
      if(me.eq.0) then
      close(unit=lfn)
      endif
c
      if(keepr.gt.0) keepr=keepr+1
c
      call timer_stop(205)
      tneed=timer_wall(205)
      call timer_reset(205)
      call timer_start(205)
      left=util_time_remaining(irtdb)
      tleft=dble(left)
      i=0
      if(left.lt.0) then
      i=1
      elseif(tleft.gt.two*tneed) then
      i=1
      endif
      call ga_brdcst(mrg_d48,i,ma_sizeof(mt_int,1,mt_byte),0)
      lstop=i.eq.0
c
      return
c
 9999 continue
      call md_abort('Unable to open restart for writing',me)
      return
      end
c $Id$
