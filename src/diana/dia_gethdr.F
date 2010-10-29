      subroutine dia_gethdr
c
c $Id$
c
c     read the trajectory file header
c     -------------------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      read(card(8:22),'(3i5)') ifrst,ilast,nframi
      filtrj=card(23:80)
c
      fmttrj='trj'
      if(index(filtrj,'.sco').gt.0) fmttrj='sco'
      if(index(filtrj,'.coo').gt.0) fmttrj='coo'
c
      call dia_rdhdr(byte_mb(i_snam))
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      call dia_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call dia_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
      lsel=.false.
c
      if(me.eq.0) then
      write(*,'(/,a,a,/)') ' Trajectory file header from ',
     + filtrj(1:index(filtrj,' ')-1)
      endif
c
      return
      end
