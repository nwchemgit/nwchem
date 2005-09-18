c
c $Id: task_qmmm_energy.F,v 1.15 2005-09-18 23:06:38 marat Exp $
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
      character*255 buf
      logical save_rst
      logical do_print
      character*30 operation
      double precision energy
      integer fn
      integer trn(3)
      character*30 pname
      character*80 geom_file
      character*3  ftype
      logical master
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
      master = qmmm_master()
c
      if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
c
      if (.not.rtdb_cget(rtdb,"qmmm:trajectory_file",1,filename)) 
     >      call errquit(pname//'failed to get trajectory filename',0,
     >       RTDB_ERR)
c
      if(index(filename,".xyzi").ne.0) then
        ftype = "xyz"
      else if(index(filename,".trj").ne.0) then
        ftype = "trj"
      else if(index(filename,".tri").ne.0) then
        ftype = "trj"
      else
        call errquit(pname//" unknown trajectory format "//
     >               filename,0,0)
      end if
c
c     convert xyz trajectory into mm trajectory
c     -----------------------------------------
      if(ftype.eq."xyz") then
        i = index(filename,".xyzi")
        buf = filename
        filename = buf(1:i)//"tri"
        call mm_create_trj_from_xyz(buf,filename)
      end if
c
      call util_file_name_resolve(filename,.false.)
      if(master) then
        open(unit=fn,file=filename,
     +   form='formatted',status='old',err=9999)
      end if
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)

      if (.not.rtdb_cget(rtdb,"qmmm:geom_file",1,geom_file)) 
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
 23   continue
c
      if(master) 
     >  close(fn)
      return

 9999 continue
      call errquit('Failed to open file '//filename,0,0)

      end
cc
c      function task_qmmm_energy_path(rtdb)
c      implicit none
cc
c#include "rtdb.fh"
c#include "util.fh"
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "qmmm_params.fh"
c#include "qmmm.fh"
cc
c      integer rtdb
c      logical task_qmmm_energy_path
cc     external functions
c      logical qmmm_energy_gradient
c      external qmmm_energy_gradient 
c
cc     local variables
c      integer i,iframe
c      logical ignore
c      character*32 interface
c      character*255 filename
c      character*255 buf
c      logical save_rst
c      logical do_print
c      character*30 operation
c      double precision energy
c      integer fn
c      integer trn(3)
c      character*30 pname
c      character*80 geom_file
c      character*3  ftype
c      logical master
c      character*255 prefix
c      integer nbeads
cc
c      logical  mm_read_frame
c      external mm_read_frame
c      logical  mm_skip_frame
c      external mm_skip_frame
cc
c      pname = "task_qmmm_energy_path"
cc
c      task_qmmm_energy_path = .true.
cc
c      if(qmmm_print_debug()) 
c     $  write(*,*) "in ",pname
cc
c      master = qmmm_master()
cc
cc     file prefix from mm module
cc     --------------------------
cc      call mm_system_prefix(prefix)
cc
cc      if(.not.qmmm_get_io_unit(fn)) 
cc     >   call errquit("cannot get file number",0,0)
cc
cc     get pathway dimension (aka number of beads) and common
cc     filename prefix
c      if (.not.rtdb_get(rtdb,"qmmm:pathway_size",mt_int,1,nbeads)) 
c     >      call errquit(pname//'failed to get trajectory filename',0,
c     >       RTDB_ERR)
c      if (.not.rtdb_cget(rtdb,"qmmm:pathway_prefix",1,prefix)) 
c     >      call errquit(pname//'failed to get trajectory filename',0,
c     >       RTDB_ERR)
cc
c      do i=1,nbeads
c        if(geom_file.ne." ") 
c     >     call mm_set_solute_coord_file(geom_file)
c        call md_sp()
c        task_qmmm_energy_path = qmmm_energy_gradient(rtdb,.false.)
c        call qmmm_energy_rtdb_push(rtdb)
c        call qmmm_print_energy(rtdb)
c      end do
cc
c      if(index(filename,".xyzi").ne.0) then
c        ftype = "xyz"
c      else if(index(filename,".trj").ne.0) then
c        ftype = "trj"
c      else if(index(filename,".tri").ne.0) then
c        ftype = "trj"
c      else
c        call errquit(pname//" unknown trajectory format "//
c     >               filename,0,0)
c      end if
cc
cc     convert xyz trajectory into mm trajectory
cc     -----------------------------------------
c      if(ftype.eq."xyz") then
c        i = index(filename,".xyzi")
c        buf = filename
c        filename = buf(1:i)//"tri"
c        call mm_create_trj_from_xyz(buf,filename)
c      end if
cc
c      call util_file_name_resolve(filename,.false.)
c      if(master) then
c        open(unit=fn,file=filename,
c     +   form='formatted',status='old',err=9999)
c      end if
cc
c      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
c     >      call errquit(pname//'failed to get trn',0,
c     >       RTDB_ERR)
c
c      if (.not.rtdb_cget(rtdb,"qmmm:geom_file",1,geom_file)) 
c     >    geom_file = " "
c
c      if(.not.mm_read_frame(fn,trn(1)))
c     >    call errquit(pname//'failed to get skip frames',
c     >                 0,0)
c
c
c      do i=trn(1),trn(2),trn(3)
c        if(geom_file.ne." ") 
c     >     call mm_set_solute_coord_file(geom_file)
c        call md_sp()
c        task_qmmm_energy_path = qmmm_energy_gradient(rtdb,.false.)
c        call qmmm_energy_rtdb_push(rtdb)
c        call qmmm_print_energy(rtdb)
c
c        if(.not.mm_read_frame(fn,trn(3))) go to 23
c
c      end do
c 23   continue
cc
c      if(master) 
c     >  close(fn)
c      return
c
c 9999 continue
c      call errquit('Failed to open file '//filename,0,0)
c
c      end
cc
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



