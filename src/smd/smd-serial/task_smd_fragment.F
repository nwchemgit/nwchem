      function task_smd_fragment(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_fragment
      logical ogroups,result,ofragment
c
      logical task_smd_energy_groups
      external task_smd_energy_groups
      logical task_smd_energy_sp
      external task_smd_energy_sp
      character*30 operation,pname
c
      double precision energy
c
      energy = 0.0d0
      pname = "task_smd_fragment"
c
      call smd_groups_inquire(rtdb,ogroups)
c
      call smd_init("fragment",rtdb)
c
      if (.not. rtdb_cget(rtdb,'task:operation', 1,operation))
     $      call errquit(pname,0, RTDB_ERR)

      write(*,*) "found operation",operation
       write(*,*) "Goodbye from task_smd_fragment"
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      task_smd_fragment = .true.

      return
      END

