c
c $Id$
c
c
C> \ingroup task
C> @{
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
      logical  task_qmmm_energy_pes
      external task_qmmm_energy_pes
c
      logical  task_qmmm_energy_spi
      logical  task_qmmm_energy_ti
      external task_qmmm_energy_ti
      external task_qmmm_energy_spi
c     local variables
      logical traj
      logical ti
      logical opes
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
      if (.not. rtdb_get(rtdb, 'qmmm:pes', mt_log,1,opes))
     $     opes = .false.
c
      if (.not. rtdb_get(rtdb, 'qmmm:ti', mt_log,1,ti))
     $     ti = .false.
c
      if(operation.ne.'energy') traj = .false.
c
      if(opes) then
        task_qmmm_energy = task_qmmm_energy_pes(rtdb)
        goto 1
      end if
c
      if(ti) then
        task_qmmm_energy = task_qmmm_energy_ti(rtdb)
        goto 1
      end if
c
      if(traj) then
        task_qmmm_energy = task_qmmm_energy_spi(rtdb)
        goto 1
      end if
c
      task_qmmm_energy = task_qmmm_energy_sp(rtdb)
c
 1    continue
      return

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
#include "global.fh"

c
      integer rtdb
      logical task_qmmm_energy_spi
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c     local variables
      integer i,iframe
      logical ignore
      character*32 interface
      character*255 filename
      character*255 prefix
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
      character*30 region(3)
      integer nregion
      character*84 tag
      integer ifr,nfr 
      integer offset
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
c     query rtdb for parameters
c     -------------------------
      if (.not.rtdb_cget(rtdb,"qmmm:trajectory_file",1,filename)) 
     >      call errquit(pname//'failed to get trajectory filename',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)
      if (.not.rtdb_cget(rtdb,"qmmm:geom_file",1,geom_file)) 
     >    geom_file = " "
c
c     figure out trajectory type
c     --------------------------
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
      if(qmmm_print_debug()) 
     $  write(*,*) "starting file preparation ",pname
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
c     open trajectory file
c     --------------------
      if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
      call util_file_name_resolve(filename,.false.)
      if(master) then
        open(unit=fn,file=filename,
     +   form='formatted',status='old',err=9999)
      end if
c
c     construct restart file name prefix
c     ----------------------------------
      i = index(filename,".tri")
      prefix = filename(1:i)
c
      if(qmmm_print_debug()) 
     $  write(*,*) "finished file preparation ",pname
c
c     process region definitions
c     --------------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
      if(nregion.gt.0) then
        if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >        call errquit(pname//tag,0,RTDB_ERR)
        call qmmm_cons_free_exclusive(region(1))
      end if
c
c     total number of frames to process
c     ---------------------------------
      nfr = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      offset = trn(1)
      ifr = 0
      do i=trn(1),trn(2),trn(3)
        ignore = rtdb_delete(rtdb,"dft:converged")
        ignore = rtdb_delete(rtdb,"scf:converged")
        ifr = ifr + 1
c       read frame
        if(.not.mm_read_frame(fn,offset))
     >      call errquit(pname//'failed to get read frame',
     >                   0,0)
c       modify coords if requested
        if(geom_file.ne." ") 
     >     call mm_set_solute_coord_file(geom_file)
c       calculate energy
        call md_sp()
        task_qmmm_energy_spi = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
c       save numbered restart file
        buf = ' '
        write(buf,14) prefix(1:inp_strlen(prefix)),ifr 
        call mm_write_restart_named(buf,.false.)
        offset = trn(3)
      end do
 23   continue
c
      if(master) 
     >  close(fn)
c
      if(qmmm_print_debug()) 
     $  write(*,*) "out ",pname
c
      return

 9999 continue
      call errquit('Failed to open file '//filename,0,0)
14    format(a,i3.3,'.rst')

      end
c
      function task_qmmm_energy_ti0(rtdb)
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
      logical task_qmmm_energy_ti0
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision temp
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
      integer nfr,nfra
      character*30 pname
      character*255 geom_file
      character*255 esp_file
      character*3  ftype
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      character*30 region(3)
      integer nregion
      character*84 tag
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "task_qmmm_energy_ti0"
c
      task_qmmm_energy_ti0 = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     no region so exiting out
c     -----------------------
      if(nregion.eq.0) return

      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
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

      if (.not.rtdb_cget(rtdb,"qmmm:ti_geom_file",1,geom_file)) 
     >    geom_file = " "

      if (.not.rtdb_cget(rtdb,"qmmm:ti_esp_file",1,esp_file)) 
     >    geom_file = " "

c
c     total number of frames to process
c     ---------------------------------
      nfr = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nfr
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e1',
     &      h_e1, i_e1) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e1),1)
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e2',
     &      h_e2, i_e2) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e2),1)
c
      nfra = 0

      if(.not.mm_read_frame(fn,trn(1)))
     >    call errquit(pname//'failed to get skip frames',
     >                 0,0)


      do i=trn(1),trn(2),trn(3)
        nfra = nfra + 1
        call mm_get_temp(temp)
        call md_sp()
        task_qmmm_energy_ti0 = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+i-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        if(master) write(*,*) "ti:e1",dbl_mb(i_e1+i-1),temp
        if(.not.mm_read_frame(fn,trn(3))) go to 23

      end do
 23   continue

      rewind(fn)

      if(.not.mm_read_frame(fn,trn(1)))
     >    call errquit(pname//'failed to get skip frames',
     >                 0,0)


      do i=trn(1),trn(2),trn(3)
        call mm_get_temp(temp)
        if(geom_file.ne." ") 
     >     call mm_set_solute_coord_file(geom_file)
        call md_sp()
        task_qmmm_energy_ti0 = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+i-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        if(master) write(*,*) "ti:e2",dbl_mb(i_e2+i-1),temp

        if(.not.mm_read_frame(fn,trn(3))) go to 24

      end do
 24   continue


c
      if(master) 
     >  close(fn)
      return

 9999 continue
      call errquit('Failed to open file '//filename,0,0)

      end
c
      function task_qmmm_energy_ti(rtdb)
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
      logical task_qmmm_energy_ti
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision temp
      double precision stime
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
      integer nfr,nfra
      character*30 pname
      character*255 geom_file
      character*255 esp_file
      character*3  ftype
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      integer i_t,h_t
      double precision taver
      character*30 region(3)
      integer nregion
      character*84 tag
      integer offset 
      double precision de
      double precision a
      double precision tsum, fsum
c
      logical ti_geom, ti_esp
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "task_qmmm_energy_ti"
c
      task_qmmm_energy_ti = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     no region so exiting out
c     -----------------------
      if(nregion.eq.0) return

      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
c
      master = qmmm_master()
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
        if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn,file=filename,
     +   form='formatted',status='old',err=9999)
      end if
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)

      
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
      nfr = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nfr
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e1',
     &      h_e1, i_e1) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e1),1)
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e2',
     &      h_e2, i_e2) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e2),1)
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti t',
     &      h_t, i_t) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_t),1)
c
        if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $     call errquit('qmmm ti: failed ', 0, RTDB_ERR)
c
      nfra = 0
      tsum = 0.0d0
      fsum = 0.0d0
      offset = trn(1)
      do i=trn(1),trn(2),trn(3)
        if(.not.mm_read_frame(fn,offset))
     >      call errquit(pname//'failed to get skip frames',
     >                   0,0)
        nfra = nfra + 1
c
c       get current simulation temperature and time
c       -------------------------------------------
        call mm_get_temp(dbl_mb(i_t+nfra-1))
        call mm_get_stime(stime)
c
c       ENERGY AT THE REFERENCE CONFIGURATION
c       -------------------------------------
        call md_sp()
        task_qmmm_energy_ti = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+nfra-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

c
c       ENERGY AT THE TARGET CONFIGURATION
c       ----------------------------------
        if(ti_geom)  
     >    call mm_set_solute_coord_file(geom_file)
        if(ti_esp) then
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file))
     >        call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if

        call md_sp()
        task_qmmm_energy_ti = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+nfra-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        call ga_sync()
c
        ignore = rtdb_delete(rtdb,"qmmm:espfilename")
        call qmmm_esp_reset(rtdb)
        offset = trn(3)
c
        if(master) write(*,'(A,I4,4F12.6)') 
     >      "ti:",
     >      nfra,
     >      stime,
     >      dbl_mb(i_t+nfra-1),
     >      dbl_mb(i_e1+nfra-1),
     >      dbl_mb(i_e2+nfra-1)

        if(master) call util_flush(6)
        tsum = tsum + dbl_mb(i_t+nfra-1)
c        fsum = fsum + exp
c     >                (
c     >       627.51*(
c     >       dbl_mb(i_e1+nfra-1)-
c     >       dbl_mb(i_e2+nfra-1)
c     >       )
c     >       /(kb*dbl_mb(i_t+nfra-1))
c     >                )
c        de = -kb*(tsum/nfra)*log(fsum/nfra)
        fsum = fsum + exp
     >                (
     >       (
     >       dbl_mb(i_e1+nfra-1)-
     >       dbl_mb(i_e2+nfra-1)
     >       )
     >       /(kb_au*dbl_mb(i_t+nfra-1))
     >                )
        de = -kb_au*(tsum/nfra)*log(fsum/nfra)
        if(master) then
          write(*,*) "current free energy difference",de*627.51,
     >               nfra,tsum,fsum
        end if

      end do
c
      taver = 0.0d0
      do i=1,nfra
        taver = taver + dbl_mb(i_t+i-1)
      end do
      taver = taver/nfra
      de = 0.0d0
      do i=1,nfra
         
        a = 627.51*(
     >       dbl_mb(i_e1+i-1)-
     >       dbl_mb(i_e2+i-1)
     >       )
c     >       /(kb*dbl_mb(i_t+i-1))
     >       /(kb*taver)

        de = de + exp(a)
      end do
      de = de/nfra
      de = -kb*taver*log(de)

      if(master) then
        write(*,*) "total free energy difference",de
      end if
c
      if(.not.ma_free_heap(h_t))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nfr, MA_ERR)

      if(.not.ma_free_heap(h_e2))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nfr, MA_ERR)

      if(.not.ma_free_heap(h_e1))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nfr, MA_ERR)

c
      if(master) 
     >  close(fn)

      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
      if(master) call util_flush(6)
      return

 9999 continue
      call errquit('Failed to open file '//filename,0,0)

      end
c
      function task_qmmm_energy_ti1(rtdb)
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
      logical task_qmmm_energy_ti1
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision temp
      double precision stime
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
      integer nfr,nfra
      character*30 pname
      character*255 geom_file
      character*255 esp_file
      character*3  ftype
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      character*30 region(3)
      integer nregion
      character*84 tag
      integer offset 
      double precision de
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "task_qmmm_energy_ti1"
c
      task_qmmm_energy_ti1 = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     no region so exiting out
c     -----------------------
      if(nregion.eq.0) return

      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
c
      master = qmmm_master()
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
        if(.not.qmmm_get_io_unit(fn)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn,file=filename,
     +   form='formatted',status='old',err=9999)
      end if
c
      if (.not.rtdb_get(rtdb,"qmmm:trajectory_n",mt_int,3,trn)) 
     >      call errquit(pname//'failed to get trn',0,
     >       RTDB_ERR)

      if (.not.rtdb_cget(rtdb,"qmmm:ti_geom_file",1,geom_file)) 
     >       call errquit(pname//'failed to get ti_geom_file ',
     >                 0,0)

      if (.not.rtdb_cget(rtdb,"qmmm:ti_esp_file",1,esp_file)) 
     >    call errquit(pname//'failed to get qmmm:ti_esp_file',
     >                 0,0)

c
c     total number of frames to process
c     ---------------------------------
      nfr = MAX((trn(2)-trn(1)+trn(3))/trn(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nfr
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e1',
     &      h_e1, i_e1) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e1),1)
c
      if(.not.ma_alloc_get(MT_DBL, nfr, 'qmmm ti e2',
     &      h_e2, i_e2) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nfr, MA_ERR)
      call dfill(nfr,0,dbl_mb(i_e2),1)
c
      nfra = 0


      offset = trn(1)
      do i=trn(1),trn(2),trn(3)
        if(.not.mm_read_frame(fn,offset))
     >      call errquit(pname//'failed to get skip frames',
     >                   0,0)
        nfra = nfra + 1
        if(master) then
        write(*,*) "frame ",nfra
        call util_flush(6)
        end if
       
        call mm_get_temp(temp)
        call mm_get_stime(stime)
        call md_sp()
        task_qmmm_energy_ti1 = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+nfra-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        call mm_set_solute_coord_file(geom_file)
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)

        call md_sp()
        task_qmmm_energy_ti1 = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+nfra-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        call ga_sync()
        if(master) write(*,'(A,I4,4F12.6)') 
     >      "ti:",
     >      nfra,
     >      stime,
     >      temp,
     >      dbl_mb(i_e1+nfra-1),
     >      dbl_mb(i_e2+nfra-1)

        if(master) call util_flush(6)
        ignore = rtdb_delete(rtdb,"qmmm:espfilename")
        offset = trn(3)

      end do
 23   continue
c
      de = 0.0d0
      do i=1,nfra
         
        de = de + dbl_mb(i_e2+i-1)-
     >       dbl_mb(i_e1+i-1)

      end do
      de = de/nfra

      if(master) then
        write(*,*) "free energy difference",de
      end if
c
      if(.not.ma_free_heap(h_e2))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nfr, MA_ERR)

      if(.not.ma_free_heap(h_e1))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nfr, MA_ERR)

c
      if(master) 
     >  close(fn)

      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
      if(master) call util_flush(6)
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
#include "global.fh"
#include "stdio.fh"

c
      integer rtdb
      logical task_qmmm_energy_sp
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      integer i
      integer un
      logical ignore
      logical oregion
      logical oprint
      character*32 interface
      logical save_rst
      logical do_print
      character*30 operation
      double precision energy
      character*30 region(3)
      integer nregion
      character*84 tag
      character*30 pname
      character*255 prefix
      character*255 xyzfile

      pname = "task_qmmm_energy_sp"

c
      if(qmmm_print_debug()) 
     $  write(*,*) "in task_qmmm_energy_sp"
c
      if (.not. rtdb_get(rtdb, 'qmmm:active_region', mt_log,1,oregion))
     $     oregion = .false.
c
      if (.not. rtdb_get(rtdb, 'qmmm:region_print', mt_log,1,oprint))
     $     oprint = .false.
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
c     process region definitions
c     --------------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >   nregion = 0

      if(nregion.gt.0) then
        if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >        call errquit(pname//tag,0,RTDB_ERR)
        call qmmm_cons_free_exclusive(region(1))
        if(oprint) then
          call mm_restart_filename(prefix)
          i = index(prefix,".rst")-1
          xyzfile = prefix(1:i)//".xyzi"
          call util_file_name_resolve(xyzfile,.false.)
          if(ga_nodeid().eq.0) then
            if(.not.qmmm_get_io_unit(un))
     >       call errquit('no free unit', 0, RTDB_ERR) 
            open(unit=un,file=xyzfile,
     >         form='formatted',status='unknown')
          end if
          call qmmm_cons_print_region(un,"xyzi")
          if(ga_nodeid().eq.0) then
            close(unit=un)
          end if
          xyzfile = prefix(1:i)//".pdb"
          call util_file_name_resolve(xyzfile,.false.)
          if(ga_nodeid().eq.0) then
            if(.not.qmmm_get_io_unit(un))
     >       call errquit('no free unit', 0, RTDB_ERR) 
            open(unit=un,file=xyzfile,
     >         form='formatted',status='unknown')
          end if
          call qmmm_cons_print_region(un,"pdbi")
          if(ga_nodeid().eq.0) then
            close(unit=un)
          end if
        end if

        if(.not.oregion)   
     >  call qmmm_cons_reset()
      end if

      call ga_sync()
   
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
1     continue
      if (.not. rtdb_delete(rtdb, 'qmmm:operation'))
     $     call errquit('qmmm: failed deleting operation', 0, RTDB_ERR)
c
      end
C> @}
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

C> \ingroup task
C> @{
      function task_qmmm_energy_pes0(irtdb)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "qmmm.fh"
#include "inp.fh"
#include "msgids.fh"
#include "tcgmsg.fh"
      integer irtdb
      logical task_qmmm_energy_pes0
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical   task_qmmm_energy_sp
      external  task_qmmm_energy_sp
      logical ignore
      character*255 filename, dir,xyz,prefix
      character*255 filename0
      character*255 movecs
      character*255 movecs0
      double precision energy
      double precision rb(2)
      integer np
      integer ip
      integer in
      double precision dr
      double precision r
      integer i1,i2,i3
      logical ofile
      logical backward
      logical loadrst
      logical loadmovecs
      integer ncycles
      integer istart
      integer iend
      logical pes_spring
      character*30 region(3)
      integer nregion
      character*84 tag
 
      pname = "task_qmmm_energy_pes0"

c
c     get lower/upper  limits for pes scan
c     ------------------------------------
      pes_spring = .true.
      if (.not. rtdb_get(irtdb,'qmmm:pes_bounds',mt_dbl,
     +                    2,rb))
     + pes_spring = .false.

c
c     get number of points for pes scan
c     ------------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_ncycles',mt_int,
     +                    1,ncycles))
     + ncycles = 1
c
c     get number of points for pes scan
c     ------------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_npoints',mt_int,
     +                    1,np))
     + call errquit(pname//'Failed to get pes points',
     + 1, RTDB_ERR)
c
c     get starting point
c     ------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_istart',mt_int,
     +                    1,istart))
     +    istart = 1
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_iend',mt_int,
     +                    1,iend))
     +    iend = np
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_backward',mt_log,
     +                    1,backward))
     +  backward = .false.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_loadrst',mt_log,
     +                    1,loadrst))
     +  loadrst = .false.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_loadmovecs',mt_log,
     +                    1,loadmovecs))
     +  loadmovecs = .false.
c
c     directory name to dump restart files
c     ------------------------------------
      dir      = ' '
      call util_directory_name(dir, .false., 0)
c
c     file prefix from mm module
c     --------------------------
      call mm_system_prefix(prefix)
c
c     construct movecs file name  
      ignore = rtdb_delete(irtdb,'dft:input vectors')
      ignore = rtdb_delete(irtdb,'scf:input vectors')
      movecs0 = ' '
      call util_file_name("tmp.movecs",.false.,.false.,movecs0)        
      if (.not. rtdb_cput(irtdb,'scf:input vectors',
     +                    1,movecs0))
     + call errquit(pname//'Failed to set scf:input vectors',
     + 1, RTDB_ERR)
      if (.not. rtdb_cput(irtdb,'dft:input vectors',
     +                    1,movecs0))
     + call errquit(pname//'Failed to set dft:input vectors',
     + 1, RTDB_ERR)

c 
c     set increments and initial value for spring
c     -------------------------------------------
      if(pes_spring) then
        dr = (rb(2)-rb(1))/(np-1)
        r  = rb(1)
      end if

      ncycles = 1
      do in=1,ncycles
      if(backward) then
        i1=iend
        i2=istart
        i3=-1
      else
        i1=istart
        i2=iend
        i3=1
      end if
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(irtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     no region so exiting out
c     -----------------------
      if(nregion.eq.0) return

      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(irtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
c
c     actually do the scan
c     --------------------
      do ip=i1,i2,i3
        if(pes_spring) then
          r  = rb(1)+dr*(ip-1) 
          if(ga_nodeid().eq.0) then
            write(*,15) ip,r
          end if
          if (.not. rtdb_put(irtdb,'qmmm:r_spring',mt_dbl,
     +                       1,r))
     +    call errquit(pname//'Failed to set r_spring',
     +    1, RTDB_ERR)
        end if

        filename0 = ' '
        write(filename0,14) 
     $     prefix(1:inp_strlen(prefix)),ip
        filename = ' '
        write(filename,13) dir(1:inp_strlen(dir)), 
     $     prefix(1:inp_strlen(prefix)),ip

c       start load restart file if it exists and requested     
        if(loadrst) then
          if(ga_nodeid().eq.0) then 
            inquire(file=filename,exist=ofile)
            if(ofile) then 
             write(*,*) "pes: found restart file",
     >                 filename(1:inp_strlen(filename))
            call util_file_copy(filename(1:inp_strlen(filename)),
     >                          filename0(1:inp_strlen(filename0)))
            end if
          end if
          call ga_brdcst(msg_qmmm_misc, ofile, 
     $     MA_sizeof(MT_INT,1,MT_BYTE),0)
          
          call ga_sync()
          if(ofile) 
     >      call mm_reload_rst(irtdb,filename0)
        end if
c       end load restart 
c       start load movecs if it exists and requested     
        movecs = ' '
        write(movecs,17) dir(1:inp_strlen(dir)), 
     $   prefix(1:inp_strlen(prefix)),ip
        if(loadmovecs) then
          if(ga_nodeid().eq.0) then 
            inquire(file=movecs,exist=ofile)
            if(ofile) then 
             write(*,*) "pes: found movecs file",
     >                 movecs(1:inp_strlen(movecs))
            call util_file_copy(movecs(1:inp_strlen(movecs)),
     >                          movecs0(1:inp_strlen(movecs0)))
            end if
          end if
          call ga_sync()
        end if
c       end load restart 
        call qmmm_esp_reset(irtdb)
        status = task_qmmm_energy_sp(irtdb)
        if(ga_nodeid().eq.0) then 
           write(*,*) "pes: saving restart file",
     >                 filename
        end if
        call mm_write_restart_named(filename,.true.)
        if(ga_nodeid().eq.0) then 
           write(*,*) "pes: saving movecs file",
     >                 movecs0(1:inp_strlen(movecs0)),
     >                 movecs(1:inp_strlen(movecs))
          call util_file_copy(movecs0(1:inp_strlen(movecs0)),
     >                        movecs(1:inp_strlen(movecs)))
        end if
        call ga_sync()

      end do
      backward = .not.backward
      loadrst = .true.
      loadmovecs = .true.
      end do
      task_qmmm_energy_pes0 = status

 13      format(a,'/',a,'_pes',i3.3,'.rst')
 14      format(a,'_pes',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)
 17      format(a,'/',a,'_pes',i3.3,'.movecs')

      return
      end

      function task_qmmm_energy_pes(irtdb)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "qmmm.fh"
#include "inp.fh"
#include "msgids.fh"
      integer irtdb
      logical task_qmmm_energy_pes
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical   task_qmmm_energy_sp
      external  task_qmmm_energy_sp
      logical ignore
      character*255 filename, dir,xyz,prefix
      character*255 filename0
      character*255 movecs
      character*255 movecs0
      double precision energy
      double precision rb(2)
      integer np
      integer ip
      integer in
      double precision dr
      double precision r
      integer i1,i2,i3
      logical ofile
      logical backward
      logical loadrst
      logical loadmovecs
      integer ncycles
      integer istart
      integer iend
      logical pes_spring
      character*30 region(3)
      integer nregion
      character*84 tag
      character*3 suffix
 
      pname = "task_qmmm_energy_pes"
      suffix = "neb"
c
c     get lower/upper  limits for pes scan
c     ------------------------------------
      pes_spring = .true.
      if (.not. rtdb_get(irtdb,'qmmm:pes_bounds',mt_dbl,
     +                    2,rb))
     + pes_spring = .false.

c
c     get number of points for pes scan
c     ------------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_ncycles',mt_int,
     +                    1,ncycles))
     + ncycles = 1
c
c     get number of points for pes scan
c     ------------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_npoints',mt_int,
     +                    1,np))
     + call errquit(pname//'Failed to get pes points',
     + 1, RTDB_ERR)
c
c     get starting point
c     ------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_istart',mt_int,
     +                    1,istart))
     +    istart = 1
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_iend',mt_int,
     +                    1,iend))
     +    iend = np
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_backward',mt_log,
     +                    1,backward))
     +  backward = .false.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_loadrst',mt_log,
     +                    1,loadrst))
     +  loadrst = .false.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_loadmovecs',mt_log,
     +                    1,loadmovecs))
     +  loadmovecs = .false.
c
c     directory name to dump restart files
c     ------------------------------------
      dir      = ' '
      call util_directory_name(dir, .false., 0)
c
c     file prefix from mm module
c     --------------------------
      call mm_system_prefix(prefix)
c
c     construct movecs file name  
      ignore = rtdb_delete(irtdb,'dft:input vectors')
      ignore = rtdb_delete(irtdb,'scf:input vectors')
      movecs0 = ' '
      call util_file_name("tmp.movecs",.false.,.false.,movecs0)        
      if (.not. rtdb_cput(irtdb,'scf:input vectors',
     +                    1,movecs0))
     + call errquit(pname//'Failed to set scf:input vectors',
     + 1, RTDB_ERR)
      if (.not. rtdb_cput(irtdb,'dft:input vectors',
     +                    1,movecs0))
     + call errquit(pname//'Failed to set dft:input vectors',
     + 1, RTDB_ERR)

c 
c     set increments and initial value for spring
c     -------------------------------------------
      if(pes_spring) then
        dr = (rb(2)-rb(1))/(np-1)
        r  = rb(1)
      end if

      ncycles = 1
      do in=1,ncycles
      if(backward) then
        i1=iend
        i2=istart
        i3=-1
      else
        i1=istart
        i2=iend
        i3=1
      end if
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(irtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     no region so exiting out
c     -----------------------
      if(nregion.eq.0) return

      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(irtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
c
c     actually do the scan
c     --------------------
      do ip=i1,i2,i3
        if(pes_spring) then
          r  = rb(1)+dr*(ip-1) 
          if(ga_nodeid().eq.0) then
            write(*,15) ip,r
          end if
          if (.not. rtdb_put(irtdb,'qmmm:r_spring',mt_dbl,
     +                       1,r))
     +    call errquit(pname//'Failed to set r_spring',
     +    1, RTDB_ERR)
        end if

        filename0 = ' '
        write(filename0,14) 
     $     prefix(1:inp_strlen(prefix)),ip
        filename = ' '
        write(filename,13) dir(1:inp_strlen(dir)), 
     $     prefix(1:inp_strlen(prefix)),ip

c       start load restart file if it exists and requested     
        if(loadrst) then
          if(ga_nodeid().eq.0) then 
            inquire(file=filename,exist=ofile)
            if(ofile) then 
             write(*,*) "pes: found restart file",
     >                 filename(1:inp_strlen(filename))
            call util_file_copy(filename(1:inp_strlen(filename)),
     >                          filename0(1:inp_strlen(filename0)))
            end if
          end if
          call ga_brdcst(msg_qmmm_misc, ofile,
     $     MA_sizeof(MT_INT,1,MT_BYTE),0)
          
          call ga_sync()
          if(ofile) 
     >      call mm_reload_rst(irtdb,filename0)
        end if
c       end load restart 
c       start load movecs if it exists and requested     
        movecs = ' '
        write(movecs,17) dir(1:inp_strlen(dir)), 
     $   prefix(1:inp_strlen(prefix)),ip
        if(loadmovecs) then
          if(ga_nodeid().eq.0) then 
            inquire(file=movecs,exist=ofile)
            if(ofile) then 
             write(*,*) "pes: found movecs file",
     >                 movecs(1:inp_strlen(movecs))
            call util_file_copy(movecs(1:inp_strlen(movecs)),
     >                          movecs0(1:inp_strlen(movecs0)))
            end if
          end if
          call ga_sync()
        end if
c       end load restart 
        call qmmm_esp_reset(irtdb)
        status = task_qmmm_energy_sp(irtdb)
        if(ga_nodeid().eq.0) then 
           write(*,*) "pes: saving restart file",
     >                 filename
        end if
        call mm_write_restart_named(filename,.true.)
        if(ga_nodeid().eq.0) then 
           write(*,*) "pes: saving movecs file",
     >                 movecs0(1:inp_strlen(movecs0)),
     >                 movecs(1:inp_strlen(movecs))
          call util_file_copy(movecs0(1:inp_strlen(movecs0)),
     >                        movecs(1:inp_strlen(movecs)))
        end if
        call ga_sync()

      end do
      backward = .not.backward
      loadrst = .true.
      loadmovecs = .true.
      end do
      task_qmmm_energy_pes = status

 13      format(a,'/',a,'_neb',i3.3,'.rst')
 14      format(a,'_neb',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)
 17      format(a,'/',a,'_pes',i3.3,'.movecs')

      return
      end
C> @}


