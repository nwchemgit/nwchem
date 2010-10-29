      subroutine argos_diana_gethdr
c
c $Id$
c
c     read the trajectory file header
c     -------------------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      read(card(8:17),'(2i5)') ifrst,ilast
      filtrj=card(18:80)
c
      fmttrj='trj'
      if(index(filtrj,'.sco').gt.0) fmttrj='sco'
      if(index(filtrj,'.coo').gt.0) fmttrj='coo'
c
      call argos_diana_rdhdr(byte_mb(i_snam))
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      call argos_diana_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call argos_diana_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
      lsel=.false.
c
      if(me.eq.0) then
      write(*,'(/,a,a,/)') ' Trajectory file header from ',
     + filtrj(1:index(filtrj,' ')-1)
      endif
c
      return
      end
