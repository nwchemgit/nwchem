c
c $Id: task_qmmm_energy.F,v 1.10 2005-08-25 23:52:22 marat Exp $
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
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      logical ignore
      character*32 interface
      logical save_rst
      logical do_print
      character*30 operation
      double precision energy

      if(qmmm_print_debug()) 
     $  write(*,*) "in task_qmmm_energy"
c
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
c
      if (.not. rtdb_cput(rtdb, 'qmmm:operation', 1, "energy"))
     $     call errquit('qmmm: failed setting operation', 0, RTDB_ERR)
c
      save_rst = .false.
      if (operation.eq.'optimize' .or.
     >    operation.eq.'saddle') save_rst = .true.
c
      do_print = .false.
      if (operation.eq.'energy') do_print = .true.
c
      interface = qmmm_get_interface()
c
      task_qmmm_energy = .true.
c      
      if (interface.eq.'mm') then
         call qmmm_geom_push_active(rtdb)
         call md_sp()
         if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)
c
         if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy/cau2kj))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c
         call qmmm_print_energy2(rtdb)
c
      else if (interface.eq.'qm') then
         call qmmm_geom_push_active(rtdb)
         call md_sp_qmmm()
         task_qmmm_energy= qmmm_energy_gradient(rtdb,.false.)
         call qmmm_energy_rtdb_push(rtdb)
         call qmmm_print_energy(rtdb)
         if(save_rst)
     >    call mm_write_restart()
         if(do_print)
     >     call mm_print_info()
      endif
c
      if (.not. rtdb_delete(rtdb, 'qmmm:operation'))
     $     call errquit('qmmm: failed deleting operation', 0, RTDB_ERR)
c
      end
c
      function qmmm_energy(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "global.fh"
c
      integer rtdb
      logical qmmm_energy
c     local variables
      character*30 operation
      logical status
      logical do_print
      character*255 theory
      double precision energy
      double precision qm_energy
      double precision mm_energy
      double precision eatoms
c     external functions
      logical ignore_mm
      logical ignore_qm
      logical  task_energy_doit
      external task_energy_doit
      character*30 pname

      pname ="qmmm energy"

      if(qmmm_print_debug())
     >       write(*,*) "in "//pname

c      call util_print_push() 
c      call util_print_rtdb_load(rtdb, 'qmmm')
      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

      do_print = .false.
      if (operation.eq.'energy') do_print = .true.

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)
      call qmmm_geom_print_current(rtdb)
c
      call qmmm_geom_push_active(rtdb)
      call mm_task_sp()
c
c     create full geometry
c     -------------------
      call qmmm_geom_create_full(rtdb)
      call qmmm_bq_coord_update()

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,qm_energy)
c
      call qmmm_geom_restore(rtdb)

c     combine QM and MM energies
c     -------------------------
      mm_energy=0.0d0
      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)
      mm_energy = mm_energy/cau2kj
      
      if (.not. rtdb_get(rtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm)) then
         ignore_mm = .false.
      end if

      if (.not. rtdb_get(rtdb,'qmmm:ignore_qm',mt_log,1,ignore_qm)) then
         ignore_qm = .false.
      end if

      if(ignore_mm) mm_energy = 0.0d0
      if(ignore_qm) qm_energy = 0.0d0

      energy = qm_energy + mm_energy-eatoms

      if (.not. rtdb_put(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)


      call qmmm_print_energy(rtdb)

      if(do_print)
     >  call mm_print_info()

      qmmm_energy=status

c      call util_print_pop()

      if(qmmm_print_debug())
     >       write(*,*) "out "//pname

      end

      subroutine qmmm_energy_rtdb_push(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      double precision eatoms
      double precision energy
      double precision mm_energy
      double precision qm_energy

      if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy',0,RTDB_ERR)

      if (.not. rtdb_get(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed getting  md energy',0,RTDB_ERR)

      mm_energy = mm_energy/cau2kj
      energy = qm_energy + mm_energy-eatoms

      if (.not. rtdb_put(rtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      end



