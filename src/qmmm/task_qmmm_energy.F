c
c $Id: task_qmmm_energy.F,v 1.3 2005-03-04 01:16:25 marat Exp $
c
c
      function task_qmmm_energy(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      logical task_qmmm_energy
c     external functions
      logical qmmm_energy
      external qmmm_energy 

c     local variables
      logical ignore
      character*32 interface
      double precision energy
c
      interface = qmmm_get_interface()
c
      task_qmmm_energy = .true.
c      
      if (interface.eq.'mm') then
         call qmmm_get_qm_coord_qm(rtdb)
         call qmmm_push_qm_coord_mm(rtdb)
         call md_sp()
         if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

         if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy/cau2kj))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      else if (interface.eq.'qm') then
         task_qmmm_energy= qmmm_energy(rtdb)
      endif
c
      end



