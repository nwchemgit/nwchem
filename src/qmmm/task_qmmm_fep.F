c
c $Id: task_qmmm_fep.F,v 1.1 2006-02-14 19:02:03 marat Exp $
c
c
      function task_qmmm_fep(rtdb)
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
      logical task_qmmm_fep
      logical ostatus
c     external functions
      logical  qmmm_fep
      external qmmm_fep
      logical  qmmm_dynamics
      external qmmm_dynamics
      character*30 pname
      character*255 prefix
      character*255 filetrj
      integer in,ncycles
      integer nf,nfs
      logical master

      pname = "task_qmmm_fep"
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
      if (.not. rtdb_cput(rtdb, 'task:operation', 1, "dynamics"))
     $    call errquit(pname,0, RTDB_ERR)
c
      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1
c
      master = ga_nodeid().eq.0
c
c     trajectory filename
c     -------------------
      call mm_trajectory_filename(filetrj)
c      filetrj=prefix(1:inp_strlen(prefix))//".trj"
      
      call mm_get_nfc(nf,nfs)
      do in=1,ncycles
c        
c
c       First we do dynamics to generate trajectory
c       -------------------------------------------
        ostatus = qmmm_dynamics(rtdb)
c
        if (.not.rtdb_cput(rtdb,"qmmm:trajectory_file",1,filetrj)) 
     >        call errquit(pname//'failed to put trajectory filename',0,
     >         RTDB_ERR)
c
        if(master)
     >  call util_file_copy(filetrj,"tmp.trj")
        call ga_sync()

        call mm_set_nfc(0,0)
        ostatus = qmmm_fep(rtdb)     
        call mm_set_nfc(nf,nfs)
        if (.not.rtdb_put(rtdb,"qmmm:extend",mt_log,1,.true.)) 
     >      call errquit(pname,0, RTDB_ERR)
        call mm_end()
        call mm_init(rtdb)
c        if(master) 
c     >   call mm_rewind_trajectory()
      end do
c
      return

      end
c
      function qmmm_fep(rtdb)
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
      logical qmmm_fep
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision stime
      integer i
      logical ignore
      character*255 filetrj
      character*255 filename
      character*255 buf
      integer fn_trj,fn
      integer trn(3)
      integer in,nf,nf0
      character*30 pname
      character*255 geom_file
      character*255 esp_file
      character*255 buffer
      character*3  ftype
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      integer i_t,h_t
      integer offset 
      double precision de
      double precision tsum, fsum
      double precision taver,faver
      double precision lambda
      logical oextend
c
      logical ti_geom, ti_esp
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "qmmm_fep"
c
      qmmm_fep = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
      master = qmmm_master()
c
c     fep filename
c     ---------------
      call mm_system_prefix(buffer)
      filename=buffer(1:inp_strlen(buffer))//".thm"
c
      if (.not.rtdb_cget(rtdb,"qmmm:trajectory_file",1,filetrj)) 
     >      call errquit(pname//'failed to get trajectory filetrj',0,
     >       RTDB_ERR)
c
      if(index(filetrj,".xyzi").ne.0) then
        ftype = "xyz"
      else if(index(filetrj,".trj").ne.0) then
        ftype = "trj"
      else if(index(filetrj,".tri").ne.0) then
        ftype = "trj"
      else
        call errquit(pname//" unknown trajectory format "//
     >               filetrj,0,0)
      end if
c
c     convert xyz trajectory into mm trajectory
c     -----------------------------------------
      if(ftype.eq."xyz") then
        i = index(filetrj,".xyzi")
        buf = filetrj
        filetrj = buf(1:i)//"tri"
        call mm_create_trj_from_xyz(buf,filetrj)
      end if
c
      if(master) then
        if(.not.qmmm_get_io_unit(fn_trj)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn_trj,file=filetrj,
     +   form='formatted',status='old',err=998)
        if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn,file=filename,
     +   form='formatted',status='unknown',err=999)
      end if
c
      if (.not.rtdb_get(rtdb,"qmmm:lambda",mt_dbl,1,lambda)) 
     >      lambda = 1.0d0
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)
c
      if (.not.rtdb_get(rtdb,"qmmm:extend",mt_log,1,oextend)) 
     >      oextend = .false.
c
      ti_geom = .true.
      if (.not.rtdb_cget(rtdb,"qmmm:ti_geom_file",1,geom_file)) 
     >       ti_geom=.false.

      ti_esp  = .true.
      if (.not.rtdb_cget(rtdb,"qmmm:ti_esp_file",1,esp_file)) 
     >    ti_esp = .false.

      if( (.not.ti_geom).and.(.not.ti_esp))  
     > call errquit(pname//'neither esp or geom file was specified',0,0)

c
c     total number of frames to process
c     ---------------------------------
      nf = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nf
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti e1',
     &      h_e1, i_e1) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_e1),1)
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti e2',
     &      h_e2, i_e2) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_e2),1)
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti t',
     &      h_t, i_t) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_t),1)
c
        if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $     call errquit('qmmm ti: failed ', 0, RTDB_ERR)
c
      nf0  = 0
      tsum = 0.0d0
      fsum = 0.0d0
c
     
      if(master) then
       if(oextend) then
       do i=1,10000
         read(fn,*,END=10)  nf0,taver,faver,de
       end do
10     continue
       tsum = taver*nf0
       fsum = faver*nf0
       end if
      end if

      call ga_sync()

      offset = trn(1)
      in = 0
      do i=trn(1),trn(2),trn(3)
        if(.not.mm_read_frame(fn_trj,offset))
     >      call errquit(pname//'failed to get skip frames',
     >                   0,0)
        in = in + 1
        call mm_get_temp(dbl_mb(i_t+in-1))
        call mm_get_stime(stime)
        call md_sp()
        qmmm_fep = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        if(ti_geom)  
     >    call mm_mix_solute_coord_file(geom_file,lambda)
        if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_load_mix(rtdb,lambda)
        end if

        call md_sp()
        qmmm_fep = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        call ga_sync()
        if(master) write(*,'(A,I6,4F12.6)') 
     >      "ti:",
     >      in,
     >      stime,
     >      dbl_mb(i_t+in-1),
     >      dbl_mb(i_e1+in-1),
     >      dbl_mb(i_e2+in-1)

        if(master) call util_flush(6)
        ignore = rtdb_delete(rtdb,"qmmm:espfilename")
        call qmmm_esp_reset(rtdb)
        offset = trn(3)
        if(master) then
          tsum = tsum + dbl_mb(i_t+in-1)
          fsum = fsum + exp
     >                  (
     >         (
     >         dbl_mb(i_e1+in-1)-
     >         dbl_mb(i_e2+in-1)
     >         )
     >         /(kb_au*dbl_mb(i_t+in-1))
     >                  )
          taver = tsum/(in+nf0)
          faver = fsum/(in+nf0)
          de = -kb_au*taver*log(faver)
          write(*,*) "current free energy difference",de*627.51
        end if
        call ga_sync()

      end do
 23   continue
c
      if(master) then
        write(fn,*) nf0+nf,taver,faver,de
        write(*,*) "debug",nf0+nf,tsum,fsum,de
      end if
c
      if(.not.ma_free_heap(h_t))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

      if(.not.ma_free_heap(h_e2))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

      if(.not.ma_free_heap(h_e1))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

c
      if(master) then
        close(fn_trj)
        close(fn)
      end if

      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
      if(master) call util_flush(6)
      return
 998  continue
      call errquit('Failed to open trajectory file '//filetrj,0,0)

 999  continue
      call errquit('Failed to open fep file '//filename,0,0)

      end

      function qmmm_fep_save(rtdb)
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
      logical qmmm_fep_save
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision stime
      integer i
      logical ignore
      character*255 filetrj
      character*255 filename
      character*255 buf
      integer fn_trj,fn
      integer trn(3)
      integer in,nf,nf0
      character*30 pname
      character*255 geom_file
      character*255 esp_file
      character*255 buffer
      character*3  ftype
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      integer i_t,h_t
      integer offset 
      double precision de
      double precision tsum, fsum
      logical oextend
c
      logical ti_geom, ti_esp
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "qmmm_fep_save"
c
      qmmm_fep_save = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
      master = qmmm_master()
c
c     fep filename
c     ---------------
      call mm_system_prefix(buffer)
      filename=buffer(1:inp_strlen(buffer))//".thm"
c
      if (.not.rtdb_cget(rtdb,"qmmm:trajectory_file",1,filetrj)) 
     >      call errquit(pname//'failed to get trajectory filetrj',0,
     >       RTDB_ERR)
c
      if(index(filetrj,".xyzi").ne.0) then
        ftype = "xyz"
      else if(index(filetrj,".trj").ne.0) then
        ftype = "trj"
      else if(index(filetrj,".tri").ne.0) then
        ftype = "trj"
      else
        call errquit(pname//" unknown trajectory format "//
     >               filetrj,0,0)
      end if
c
c     convert xyz trajectory into mm trajectory
c     -----------------------------------------
      if(ftype.eq."xyz") then
        i = index(filetrj,".xyzi")
        buf = filetrj
        filetrj = buf(1:i)//"tri"
        call mm_create_trj_from_xyz(buf,filetrj)
      end if
c
      if(master) then
        if(.not.qmmm_get_io_unit(fn_trj)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn_trj,file=filetrj,
     +   form='formatted',status='old',err=998)
        if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn,file=filename,
     +   form='formatted',status='unknown',err=999)
      end if
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)
c
      if (.not.rtdb_get(rtdb,"qmmm:extend",mt_log,1,oextend)) 
     >      oextend = .false.
c
      ti_geom = .true.
      if (.not.rtdb_cget(rtdb,"qmmm:ti_geom_file",1,geom_file)) 
     >       ti_geom=.false.

      ti_esp  = .true.
      if (.not.rtdb_cget(rtdb,"qmmm:ti_esp_file",1,esp_file)) 
     >    ti_esp = .false.

      if( (.not.ti_geom).and.(.not.ti_esp))  
     > call errquit(pname//'neither esp or geom file was specified',0,0)

c
c     total number of frames to process
c     ---------------------------------
      nf = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nf
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti e1',
     &      h_e1, i_e1) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_e1),1)
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti e2',
     &      h_e2, i_e2) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_e2),1)
c
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti t',
     &      h_t, i_t) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_t),1)
c
        if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $     call errquit('qmmm ti: failed ', 0, RTDB_ERR)
c
      nf0  = 0
      tsum = 0.0d0
      fsum = 0.0d0
c
     
      if(master) then
       if(oextend) then
       do i=1,10000
         read(fn,*,END=10)  nf0,tsum,fsum,de
       end do
10     continue
       end if
      end if

      call ga_sync()

      offset = trn(1)
      in = 0
      do i=trn(1),trn(2),trn(3)
        if(.not.mm_read_frame(fn_trj,offset))
     >      call errquit(pname//'failed to get skip frames',
     >                   0,0)
        in = in + 1
        call mm_get_temp(dbl_mb(i_t+in-1))
        call mm_get_stime(stime)
        call md_sp()
        qmmm_fep_save = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        if(ti_geom)  
     >    call mm_set_solute_coord_file(geom_file)
        if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
        end if

        call md_sp()
        qmmm_fep_save = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        call ga_sync()
        if(master) write(*,'(A,I6,4F12.6)') 
     >      "ti:",
     >      in,
     >      stime,
     >      dbl_mb(i_t+in-1),
     >      dbl_mb(i_e1+in-1),
     >      dbl_mb(i_e2+in-1)

        if(master) call util_flush(6)
        ignore = rtdb_delete(rtdb,"qmmm:espfilename")
        call qmmm_esp_reset(rtdb)
        offset = trn(3)
        if(master) then
          tsum = tsum + dbl_mb(i_t+in-1)
          fsum = fsum + exp
     >                  (
     >         (
     >         dbl_mb(i_e1+in-1)-
     >         dbl_mb(i_e2+in-1)
     >         )
     >         /(kb_au*dbl_mb(i_t+in-1))
     >                  )
          de = -kb_au*(tsum/(in+nf0))*log(fsum/(in+nf0))
          write(*,*) "current free energy difference",de
        end if
        call ga_sync()

      end do
 23   continue
c
      if(master) then
        write(fn,*) nf0+nf,tsum,fsum,de
        write(*,*) "debug",nf0+nf,tsum,fsum,de
      end if
c
      if(.not.ma_free_heap(h_t))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

      if(.not.ma_free_heap(h_e2))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

      if(.not.ma_free_heap(h_e1))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

c
      if(master) then
        close(fn_trj)
        close(fn)
      end if

      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
      if(master) call util_flush(6)
      return
 998  continue
      call errquit('Failed to open trajectory file '//filetrj,0,0)

 999  continue
      call errquit('Failed to open fep file '//filename,0,0)

      end
c

