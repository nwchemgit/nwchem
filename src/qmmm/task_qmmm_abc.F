c
c $Id: task_qmmm_abc.F,v 1.3 2006-05-28 01:18:41 marat Exp $
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
#include "msgids.fh"
c
      integer rtdb
      logical task_qmmm_abc
c     local variables
      logical ostatus
      character*30 pname
      logical master
      character*255 rst_in
      character*255 rst_out
      character*255 filerest
      character*255 filetop(2)
      character*255 atheory(2)
      character*255 oldtop
      character*255 logfile(2)
      logical otop
      integer fn(2)
      integer nfr
      integer in
      integer ncycles

c     external functions
      logical qmmm_end
      external qmmm_end
      logical qmmm_init
      external qmmm_init
      logical qmmm_abc_generate_frames
      external qmmm_abc_generate_frames
      logical qmmm_abc_resample
      external qmmm_abc_resample

      pname = "task_qmmm_abc"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
c     get reference and target theories
c     ---------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:abc_theory', 2, atheory)) then
        atheory(1) = "espmm"
        atheory(2) = "qmmm"
      end if
c
c     get reference and target topology files
c     ---------------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:abc_top', 2, filetop)) then
        filetop(1)=" "
        filetop(2)=" "     
        otop = .false.
      else
        call util_file_name_resolve(filetop(1),.false.)
        call util_file_name_resolve(filetop(2),.false.)
      end if
c
        if (.not.rtdb_get(rtdb,"qmmm:abc_ncycles",mt_int,1,ncycles)) 
     >      ncycles = 1
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
c     open files to dump energy information
c     -------------------------------------
      if(qmmm_master()) then
        if(.not.qmmm_get_io_unit(fn(1)))
     $    call errquit(pname//" no free unit number",0,0)
        call util_file_name('low.abc', .false., .false.,logfile(1))
        open(unit=fn(1),status="unknown",
     >          form="formatted",file=logfile(1))
        if(.not.qmmm_get_io_unit(fn(2)))
     $    call errquit(pname//" no free unit number",0,0)
        call util_file_name('high.abc', .false., .false.,logfile(2))
        open(unit=fn(2),status="unknown",
     >          form="formatted",file=logfile(2))
      end if
      call ga_brdcst(msg_qmmm_misc,
     >     fn,
     >     2*ma_sizeof(mt_int,1,mt_byte),
     >     0)
      write(*,*) "fn=",fn(1),fn(2)
      call ga_sync()
c
      if (.not. rtdb_put(rtdb,'qmmm:pes_rst_perm',mt_log,
     +                    1,.false.))
     +    call errquit(pname,0,0)

      do in=1,ncycles
c
c       generate frames for resampling
c       ------------------------------
        ostatus = qmmm_abc_generate_frames(nfr,filetop(1),
     >                                     atheory(1),rtdb)
c
c       resample with low level theory
c       ------------------------------
        if (.not. rtdb_put(rtdb,'qmmm:pes_logfile',mt_int,1,fn(1)))
     $      call errquit(pname,0,0)
        ostatus = qmmm_abc_resample(nfr,filetop(1),
     >                              rst_in,atheory(1),rtdb)
c
c       resample with low level theory
c       ------------------------------
        if (.not. rtdb_put(rtdb,'qmmm:pes_logfile',mt_int,1,fn(2)))
     $      call errquit(pname,0,0)
        ostatus = qmmm_abc_resample(nfr,filetop(2),
     >                              rst_in,atheory(2),rtdb)
      
        write(rst_out,21) rst_in(1:inp_strlen(rst_in)),nfr 
 21     format(a,i3.3,'.rst')
        if(master) then
           call util_file_copy(rst_out(1:inp_strlen(rst_out)),
     >          filerest(1:inp_strlen(filerest)))

        end if
      end do

      task_qmmm_abc = .true.
      
      close(fn(1))
      close(fn(2))
      return
      end
c
      function qmmm_abc_resample(nfr,topfile,rstfile,atheory,rtdb)
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
      character*(*) rstfile
      character*(*) atheory
      integer rtdb
      logical qmmm_abc_resample
c
      logical ostatus
      logical ignore
      logical master
      logical otop
      logical oesp
      character*255 oldtop
      character*255 oldoperation
      character*30 pname
      character*30 tmp_name
c     
      logical qmmm_end
      external qmmm_end
      logical qmmm_init
      external qmmm_init
c
      logical task_qmmm_pes
      external task_qmmm_pes
c
      pname = "qmmm_abc_resample"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
      call mm_top_file(oldtop)
      otop = topfile.ne." "
c
      if (.not. rtdb_get(rtdb,'qmmm:esp_density',
     $                   mt_log,1,oesp))
     $     oesp=.false.
c
      if(atheory.eq."mm") then
         if(.not.otop) 
     >     call errquit(pname//" need mm topology file",0,0)
      else if (atheory.eq."espmm") then
         if (.not. rtdb_put(rtdb,'qmmm:esp_density',
     $                      mt_log,1,.true.))
     $         call errquit(pname,0, RTDB_ERR)
      end if

      if(otop) then
        ignore = qmmm_end(rtdb)
        if(master) then 
            write(*,*) "copying filename", topfile
            call util_file_copy(topfile,oldtop)
        end if
        call ga_sync()
      end if
c
      if (.not. rtdb_cget(rtdb,'task:operation', 1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
      if (.not. rtdb_cput(rtdb,'task:operation', 1,'energy'))
     $      call errquit(pname,0, RTDB_ERR)
c
      ignore = qmmm_init(rtdb)
      call qmmm_bq_data_load()
c
      if (.not. rtdb_cput(rtdb,'qmmm:pes_rst_input',
     +                    1,rstfile))
     +      call errquit(pname,0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'qmmm:pes_npoints',mt_int,
     +                    1,nfr))
     +     call errquit(pname//'Failed to set pes points',
     +     1, RTDB_ERR)
c
      ostatus = task_qmmm_pes(rtdb)
c
c     restore old settings
c     --------------------
      if(.not. rtdb_cput(rtdb,'task:operation',1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'qmmm:esp_density',
     $                   mt_log,1,oesp))
     $      call errquit(pname,0, RTDB_ERR)
c
      qmmm_abc_resample = ostatus


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
c
      logical ostatus
      logical ignore
      logical master
      character*30 pname
      character*30 tmp_name
c     external functions
      logical  qmmm_abc_generate_frames_mm
      external qmmm_abc_generate_frames_mm
      logical qmmm_abc_generate_frames_esp
      external qmmm_abc_generate_frames_esp
c     
      pname = "qmmm_abc_generate_frames"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
      if(atheory.eq."mm") then
        ostatus = qmmm_abc_generate_frames_mm(nfr,topfile,rtdb)
      else if (atheory.eq."espmm") then
        ostatus = qmmm_abc_generate_frames_esp(nfr,rtdb)
      else
        tmp_name = atheory
        call errquit(pname//"unknown theory "//tmp_name,0,0)  
      end if
      return

      end
c
      function qmmm_abc_generate_frames_mm(nfr,topfile,rtdb)
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
      integer rtdb
      logical qmmm_abc_generate_frames_mm
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
      logical otop

      logical qmmm_init
      external qmmm_init

      logical qmmm_end
      external qmmm_end


      pname = "qmmm_abc_generate_frames"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c 
      otop = topfile.ne." "
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
c     tinker with topology file
c     -------------------------
      if(otop) then
        call mm_top_file(oldtop)
        savetop = "tmp.abc.top"
        if(master) then 
          call util_file_copy(oldtop,savetop)
        end if
        call util_file_name_resolve(topfile,.false.)
        if(master) then 
          write(*,*) "copying filename", topfile
          call util_file_copy(topfile,oldtop)
        end if
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
      if(otop) then
        if(master) then 
          call util_file_copy(savetop,oldtop)
          call util_file_unlink(savetop)
        end if
        call ga_sync()
      end if
c
      qmmm_abc_generate_frames_mm = .true.
c
      return

      end

      function qmmm_abc_generate_frames_esp(nfr,rtdb)
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
      integer rtdb
      logical qmmm_abc_generate_frames_esp
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
      logical oesp

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
c      ignore = qmmm_end(rtdb)
c      call ga_sync()
c
c     enforce esp density and save old setting
c     ----------------------------------------
      if (.not. rtdb_get(rtdb,'qmmm:esp_density',
     $                   mt_log,1,oesp))
     $     oesp=.false.
c
      if (.not. rtdb_put(rtdb,'qmmm:esp_density',
     $                   mt_log,1,.true.))
     $      call errquit(pname,0, RTDB_ERR)
c
c     enforce dynamics run
c     --------------------
      if (.not. rtdb_cget(rtdb,'task:operation', 1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
      if (.not. rtdb_cput(rtdb,'task:operation', 1,'dynamics'))
     $      call errquit(pname,0, RTDB_ERR)
c
c     enforce numbered restart files
c     ------------------------------
      call mm_set_keep_restart()
c
c     commence dynamics
c     -----------------
      ostatus = qmmm_dynamics(rtdb)
c
c     how many frames were generated
c     ------------------------------
      call mm_get_keep_restart_index(nfr)
c
c     restore old settings
c     --------------------
      if(.not. rtdb_cput(rtdb,'task:operation',1,oldoperation))
     $      call errquit(pname,0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'qmmm:esp_density',
     $                   mt_log,1,oesp))
     $      call errquit(pname,0, RTDB_ERR)
c
      qmmm_abc_generate_frames_esp = ostatus

      return

      end

