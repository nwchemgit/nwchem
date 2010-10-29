      subroutine argos_diana_gettop
c
c $Id$
c
c     read the topology file
c     ----------------------
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
      integer l_tmp,i_tmp
c
      filtop=card(8:80)
      call argos_diana_siztop()
      if(.not.ma_push_get(mt_int,7*msgm,'iram',l_iram,i_iram))
     + call md_abort('Could not allocate iram',0)
      if(.not.ma_push_get(mt_int,msa,'imol',l_imol,i_imol))
     + call md_abort('Could not allocate imol',0)
      if(.not.ma_push_get(mt_int,msa,'iequi',l_iequi,i_iequi))
     + call md_abort('Could not allocate iequi',0)
      if(.not.ma_push_get(mt_int,2*msb,'sbnd',l_sbnd,i_sbnd))
     + call md_abort('Could not allocate sbnd',0)
      call argos_diana_rdtop(byte_mb(i_snam),dbl_mb(i_qdat),
     + dbl_mb(i_qwdat),
     + int_mb(i_iram),int_mb(i_imol),int_mb(i_sbnd))
      if(.not.ma_push_get(mt_int,2*nsa,'tmp',l_tmp,i_tmp))
     + call md_abort('Failed to allocate temp array',0)
      call argos_diana_equivalent(byte_mb(i_snam),int_mb(i_imol),
     + int_mb(i_iequi),int_mb(i_tmp))
      if(.not.ma_pop_stack(l_tmp))
     + call md_abort('Failed to deallocate temp array',0)
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
