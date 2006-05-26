c
c $Id: task_qmmm_abc.F,v 1.1 2006-05-26 18:55:58 marat Exp $
c
      function task_qmmm_abc(rtdb)
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
      logical task_qmmm_abc
      logical ostatus
      logical ignore
c     external functions
      logical  task_qmmm_pes 
      external task_qmmm_pes 
      logical  qmmm_fep_deriv
      external qmmm_fep_deriv
      logical  qmmm_dynamics
      external qmmm_dynamics
      character*30 pname
      character*255 prefix
      character*255 filetrj
      integer in,ncycles
      integer nf,nfs
      logical master
      character*255 rst_in
      character*255 filerest
      character*255 filetop(2)
      character*255 atheory(2)
      character*255 filegeom(2)
      character*255 fileesp(2)
      character*255 esp_in(3)
      character*255 geom_in(3)
      character*255 oldtop
      character*255 logfile
      double precision lambda(3)
      logical oesp,ogeom
      integer mequi 
      logical rtdb_mode
      integer fn
      integer pes_npoints
      integer nfr

      logical qmmm_end
      external qmmm_end

      logical qmmm_init
      external qmmm_init

      logical qmmm_abc_generate_frames
      external qmmm_abc_generate_frames

      pname = "task_qmmm_abc"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
c     get reference and target topology files
c     ---------------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:abc_theory', 2, atheory))
     $    call errquit(pname//" no theory found",0,0)
      if (.not. rtdb_cget(rtdb, 'qmmm:abc_top', 2, filetop))
     $    call errquit(pname//" no topology file found",0,0)
      call util_file_name_resolve(filetop(1),.false.)
      call util_file_name_resolve(filetop(2),.false.)
c
c     get default topology filename
c     -----------------------------
      call mm_top_file(oldtop)
c
c     prefix for numbered restart file
c     --------------------------------
      call mm_restart_filename(filerest)
      rst_in = filerest(1:index(filerest,'.rst')-1)//
     +         "-"//"00"
c
      if (.not. rtdb_cput(rtdb,'task:operation', 1,'energy'))
     $      call errquit(pname,0, RTDB_ERR)
c
c     open file to dump energy information
c     ------------------------------------
      rtdb_mode = rtdb_parallel(.false.)
      if(qmmm_master()) then
      if(.not.qmmm_get_io_unit(fn))
     $  call errquit(pname//" no free unit number",0,0)
        call util_file_name('abc', .false., .false.,logfile)
        open(unit=fn,status="unknown",form="formatted",file=logfile)
        if (.not. rtdb_put(rtdb,'qmmm:pes_logfile',mt_int,1,fn))
     $    call errquit(pname,0,0)
      end if
      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)
c
      ostatus = qmmm_abc_generate_frames(nfr,filetop(1),atheory(1),rtdb)
      write(*,*) "nfr=",nfr
      if(master) then 
          write(*,*) "copying filename", filetop(2)
          call util_file_copy(filetop(2),oldtop)
      end if
      call ga_sync()
      ignore = qmmm_init(rtdb)
      call qmmm_bq_data_load()

c     do the perturbation resampling
      if (.not. rtdb_cput(rtdb,'qmmm:pes_rst_input',
     +                    1,rst_in))
     +      call errquit(pname,0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'qmmm:pes_npoints',mt_int,
     +                    1,nfr))
     +     call errquit(pname//'Failed to set pes points',
     +     1, RTDB_ERR)
      ostatus = task_qmmm_pes(rtdb)
      task_qmmm_abc = .true.

      return
      end
c
      function qmmm_abc_generate_frames(nfr,topfile,atheory,rtdb)
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
      integer nfr
      character*(*) topfile
      character*(*) atheory
      integer rtdb
      logical qmmm_abc_generate_frames
      logical ostatus
      logical ignore
c     external functions
      logical  qmmm_dynamics
      external qmmm_dynamics
      character*30 pname
      character*255 oldtheory
      character*255 oldoperation
      integer in,ncycles
      integer nf,nfs
      logical master
      character*255 oldtop
      character*255 savetop

      logical qmmm_init
      external qmmm_init

      logical qmmm_end
      external qmmm_end


      pname = "qmmm_abc_generate_frames"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname

c
c     destroy qmmm structures (if any)
c     --------------------------------
      ignore = qmmm_end(rtdb)
      call ga_sync()
c
c     deactivate qmmm for now
c     -----------------------
      if(.not. rtdb_cget(rtdb,'task:theory',1,oldtheory))
     $      call errquit(pname,0, RTDB_ERR)
      if(.not. rtdb_cput(rtdb,'task:theory',1,"md"))
     $      call errquit(pname,0, RTDB_ERR)
c
c     enforce dynamics run
c     --------------------
      if (.not. rtdb_cget(rtdb,'task:operation', 1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
      if (.not. rtdb_cput(rtdb,'task:operation', 1,'dynamics'))
     $      call errquit(pname,0, RTDB_ERR)
c
c     process rtdb through md module
c     ------------------------------
      call mm_set_irtdb(rtdb)
      call md_rdinp()
c
c     save old topology file
c     ---------------------
      call mm_top_file(oldtop)
      savetop = "tmp.abc.top"
      if(master) then 
        call util_file_copy(oldtop,savetop)
      end if
c
c     copy new topology file from perm
c     --------------------------------
      call util_file_name_resolve(topfile,.false.)
      if(master) then 
        write(*,*) "copying filename", topfile
        call util_file_copy(topfile,oldtop)
      end if
c
c     enforce numbered restart files
c     ------------------------------
      call mm_set_keep_restart()
c
c     commence dynamics
c     -----------------
      call md_start()
      call md_md()
c
c     how many frames were generated
c     ------------------------------
      call mm_get_keep_restart_index(nfr)
c
c     finish
c     ------
      call md_finish()
      call ga_sync()
c
c     restore old theory and operation
c     --------------------------------
      if(.not. rtdb_cput(rtdb,'task:theory',1,oldtheory))
     $      call errquit(pname,0, RTDB_ERR)
      if(.not. rtdb_cput(rtdb,'task:operation',1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
c
c     restore old topology file
c     ---------------------
      if(master) then 
        call util_file_copy(savetop,oldtop)
        call util_file_unlink(savetop)
      end if
c
      return

      end
cc
c      function task_qmmm_abc(rtdb)
c      implicit none
cc
c#include "rtdb.fh"
c#include "util.fh"
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "qmmm_params.fh"
c#include "qmmm.fh"
c#include "global.fh"
cc
c      integer rtdb
c      logical task_qmmm_abc
c      logical ostatus
c      logical ignore
cc     external functions
c      logical  task_qmmm_pes 
c      external task_qmmm_pes 
c      logical  qmmm_fep_deriv
c      external qmmm_fep_deriv
c      logical  qmmm_dynamics
c      external qmmm_dynamics
c      character*30 pname
c      character*255 prefix
c      character*255 filetrj
c      integer in,ncycles
c      integer nf,nfs
c      logical master
c      character*255 rst_in
c      character*255 filerest
c      character*255 filetop(2)
c      character*255 filegeom(2)
c      character*255 fileesp(2)
c      character*255 esp_in(3)
c      character*255 geom_in(3)
c      character*255 oldtop
c      character*255 logfile
c      double precision lambda(3)
c      logical oesp,ogeom
c      integer mequi 
c      logical rtdb_mode
c      integer fn
c      integer pes_npoints
c
c      logical qmmm_end
c      external qmmm_end
c
c      logical qmmm_init
c      external qmmm_init
c
c      pname = "task_qmmm_abc"
c      master = ga_nodeid().eq.0
c      if(qmmm_print_debug()) 
c     $  write(*,*) "in "//pname
cc
cc     get reference and target topology files
cc     ---------------------------------------
c      if (.not. rtdb_cget(rtdb, 'qmmm:abc_theory', 2, atheory))
c     $    call errquit(pname//" no theory found",0,0)
c      if (.not. rtdb_cget(rtdb, 'qmmm:abc_top', 2, filetop))
c     $    call errquit(pname//" no topology file found",0,0)
c      call util_file_name_resolve(filetop(1),.false.)
c      call util_file_name_resolve(filetop(2),.false.)
cc
cc     open file to dump energy information
cc     ------------------------------------
c      rtdb_mode = rtdb_parallel(.false.)
c      if(qmmm_master()) then
c      if(.not.qmmm_get_io_unit(fn))
c     $  call errquit(pname//" no free unit number",0,0)
c        call util_file_name('abc', .false., .false.,logfile)
c        open(unit=fn,status="unknown",form="formatted",file=logfile)
c        if (.not. rtdb_put(rtdb,'qmmm:pes_log',mt_int,1,fn))
c     $    call errquit(pname,0,0)
c      end if
c      call ga_sync()
c      rtdb_mode = rtdb_parallel(rtdb_mode)
cc
c      call mm_system_prefix(prefix)
c
c      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
c     +  ncycles = 1
cc
cc     get default topology filename
cc     -----------------------------
c      call mm_top_file(oldtop)
cc
cc     prefix for numbered restart file
cc     --------------------------------
c      call mm_restart_filename(filerest)
c      rst_in = filerest(1:index(filerest,'.rst')-1)//
c     +         "-"//"00"
cc
cc     commence the actual run
cc     -----------------------
c      do in=1,ncycles
c        ignore =  qmmm_end(rtdb)
c        if (.not. rtdb_cput(rtdb,'task:operation', 1,"dynamics"))
c     $      call errquit(pname,0, RTDB_ERR)
c        if(master) then 
c          write(*,*) "copying filename", filetop(1)
c          call util_file_copy(filetop(1),oldtop)
c        end if
c        call ga_sync()
c        ignore =  qmmm_init(rtdb)
cc       enforce numbered restart files
c        call mm_set_keep_restart()
cc       dynamics task to generate restart files
c        ostatus = qmmm_dynamics(rtdb)
c        call mm_get_keep_restart_index(pes_npoints)
c        write(*,*) "pes_npoints",pes_npoints
c        ignore =  qmmm_end(rtdb)
cc       start processing restart file sequence
c        if(master) then 
c          write(*,*) "copying filename", filetop(2)
c          call util_file_copy(filetop(2),oldtop)
c        end if
c        call ga_sync()
c
c        if (.not. rtdb_cput(rtdb,'task:operation', 1,"energy"))
c     $      call errquit(pname,0, RTDB_ERR)
c        ignore = qmmm_init(rtdb)
c        call qmmm_bq_data_load()
c
cc       do the perturbation resampling
c        if (.not. rtdb_cput(rtdb,'qmmm:pes_rst_input',
c     +                    1,rst_in))
c     +      call errquit(pname,0, RTDB_ERR)
c        if (.not. rtdb_put(rtdb,'qmmm:pes_npoints',mt_int,
c     +                    1,pes_npoints))
c     +     call errquit(pname//'Failed to set pes points',
c     +     1, RTDB_ERR)
c        ostatus = task_qmmm_pes(rtdb)
cc       request to reuse collected abc data for ncycle.gt.1
c        if (.not.rtdb_put(rtdb,"qmmm:extend",mt_log,1,.true.)) 
c     >      call errquit(pname,0, RTDB_ERR)
c      end do
cc
cc     reset any region constraints
c      call qmmm_cons_reset()
c      task_qmmm_abc = ostatus
c      return
c
c      end
c

