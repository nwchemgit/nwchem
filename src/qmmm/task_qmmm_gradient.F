c
c $Id: task_qmmm_gradient.F,v 1.1 2004-10-09 01:09:01 marat Exp $
c
c
      function task_qmmm_gradient(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
c
      integer rtdb
      logical  task_qmmm_gradient
c
      logical  qmmm_gradient
      external qmmm_gradient
c     local variables
      character*32 interface
      double precision energy

      if (.not. rtdb_cget(rtdb, 'qmmm:interface', 1, interface))
     $     call errquit('qmmm_interface: failed ',0,0)
c
      
      task_qmmm_gradient = .true.
      if (interface.eq.'mm') then
         call qmmm_get_qm_coord_qm(rtdb)
         call qmmm_push_qm_coord_mm(rtdb)
         call md_sp()
         if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

         if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy/cau2kj))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      else if (interface.eq.'qm') then
         task_qmmm_gradient = qmmm_gradient(rtdb)
      endif
c
      end



