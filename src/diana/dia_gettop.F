      subroutine dia_gettop
c
c $Id: dia_gettop.F,v 1.1 2007-01-27 02:26:03 d3j191 Exp $
c
c     read the topology file
c     ----------------------
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
      filtop=card(8:80)
      call dia_siztop()
      if(.not.ma_push_get(mt_int,7*msgm,'iram',l_iram,i_iram))
     + call md_abort('Could not allocate iram',0)
      if(.not.ma_push_get(mt_int,msa,'imol',l_imol,i_imol))
     + call md_abort('Could not allocate imol',0)
      if(.not.ma_push_get(mt_int,2*msb,'sbnd',l_sbnd,i_sbnd))
     + call md_abort('Could not allocate sbnd',0)
      call dia_rdtop(byte_mb(i_snam),dbl_mb(i_qdat),dbl_mb(i_qwdat),
     + int_mb(i_iram),int_mb(i_imol),int_mb(i_sbnd))
c
      if(me.eq.0) then
      if(ltop) then
      write(*,'(a,a)') ' Topology read from ',
     + filtop(1:index(filtop,' ')-1)
      else
      write(*,'(a,a)') ' Topology could not be read from ',
     + filtop(1:index(filtop,' ')-1)
      endif
      endif
c
      return
      end
