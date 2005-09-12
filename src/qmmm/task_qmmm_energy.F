c
c $Id: task_qmmm_energy.F,v 1.13 2005-09-12 06:00:42 marat Exp $
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
      logical  task_qmmm_energy_sp
      external task_qmmm_energy_sp
c
      logical  task_qmmm_energy_spi
      external task_qmmm_energy_spi
c     local variables
      logical traj
      character*30 operation

      if(qmmm_print_debug()) 
     $  write(*,*) "in task_qmmm_energy"
c
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
c
      if (.not. rtdb_get(rtdb, 'qmmm:trajectory', mt_log,1,traj))
     $     traj = .false.
c
      if(operation.ne.'energy') traj = .false.
c
      if(traj) then
        task_qmmm_energy = task_qmmm_energy_spi(rtdb)
      else
        task_qmmm_energy = task_qmmm_energy_sp(rtdb)
      end if
c
      end
c
      function task_qmmm_energy_spi(rtdb)
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
      logical task_qmmm_energy_spi
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      integer i,iframe
      logical ignore
      character*32 interface
      character*255 filename
      logical save_rst
      logical do_print
      character*30 operation
      double precision energy
      integer fn
      integer trn(3)
      character*30 pname
      character*80 geom_file
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "task_qmmm_energy_spi"
c
      task_qmmm_energy_spi = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
      if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
c
      if (.not.rtdb_cget(rtdb,"qmmm:trajectory_file",1,filename)) 
     >      call errquit(pname//'failed to get trajectory filename',0,
     >       RTDB_ERR)
c
      open(unit=fn,file=filename,
     + form='formatted',status='old',err=9999)
      rewind(fn)
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)

      if (rtdb_cget(rtdb,"qmmm:geom_file",1,geom_file)) 
     >    geom_file = " "

      if(.not.mm_read_frame(fn,trn(1)))
     >    call errquit(pname//'failed to get skip frames',
     >                 0,0)


      do i=trn(1),trn(2),trn(3)
        if(geom_file.ne." ") 
     >     call mm_set_solute_coord_file(geom_file)
        call md_sp()
        task_qmmm_energy_spi = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)

        if(.not.mm_read_frame(fn,trn(3))) go to 23

      end do
      close(fn)
23    continue
      return

 9999 continue
      call errquit('Failed to open file '//filename,0,0)

      end
c
      function task_qmmm_energy_sp(rtdb)
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
      logical task_qmmm_energy_sp
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
      task_qmmm_energy_sp = .true.
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
         task_qmmm_energy_sp= qmmm_energy_gradient(rtdb,.false.)
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



