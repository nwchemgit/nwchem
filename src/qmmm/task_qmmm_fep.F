c
c $Id$
c
c
C> \ingroup task
C> @{
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
#include "stdio.fh"
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
      character*255 filegeom(2)
      character*255 fileesp(2)
      character*255 esp_in(2)
      character*255 geom_in(2)
      double precision lambda(2)
      logical oesp,ogeom,oderiv,task_qmmm_fep_deriv
      logical oesp_density
      integer mequi 
      integer nsamples
      integer trn(3)
      character*30 region(3)
      integer nregion
      character*84 tag

      pname = "task_qmmm_fep"

      if (.not. rtdb_cput(rtdb, 'task:operation', 1, "dynamics"))
     $    call errquit(pname,0, RTDB_ERR)
c
c     check for correct region
c     ------------------------
      tag ="qmmm:region"
      if (.not.rtdb_cget(rtdb,tag,nregion,region)) then
        nregion = 1
        region(1)="mm"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
        if (.not.rtdb_cput(rtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
      else if(region(1).ne."mm_solute".or.
     >        region(1).ne."solvent".or.
     >        region(1).ne."mm") then
            call errquit("only classical regions supported",0,RTDB_ERR)
      end if
c
      oesp_density = .false.
      if (.not. rtdb_cget(rtdb, 'qmmm:density_treatment', 
     $                   1, tag)) then
        
        oesp_density = .true.
        if (.not. rtdb_put(rtdb, 'qmmm:esp_density', mt_log, 
     $                   1, oesp_density))   
     $       call errquit(pname,0,RTDB_ERR)
      tag = "espfit"
      if (.not. rtdb_cput(rtdb, 'qmmm:density_treatment', 
     $                   1, tag))   
     $       call errquit(pname,0,RTDB_ERR)
      end if
      
      if(tag.ne."espfit") then
        if(ga_nodeid().eq.0) then
        call banner(luout,"density espfit is strongly reccomended ",
     $              "*",.true.,.true.,.true.)
        end if
      end if
      if (rtdb_get(rtdb,'qmmm:nsamples',mt_int,1,nsamples)) then
        call mm_set_high_precision()
        call mm_set_nfc(1,1)
        call mm_set_mdacq(nsamples)
        trn(1) = 1
        trn(2) = nsamples
        trn(3) = 1
        if (.not.rtdb_put(rtdb,"qmmm:trajectory_n",mt_int,3,trn))
     >      call errquit(pname//'failed to put trn',0,
     >       RTDB_ERR)
        call mm_init(rtdb)
      else
        call errquit(pname//'define nsamples now',0,RTDB_ERR)
      end if

      if (.not. rtdb_get(rtdb,'qmmm:fep_deriv',mt_log,1,oderiv))
     +  oderiv = .false.

      if(oderiv) then
c change task_qmmm_fep_deriv from a subroutine call to a function call
        ostatus=task_qmmm_fep_deriv(rtdb)
c        call task_qmmm_fep_deriv(rtdb)
c        ostatus = .true.
        goto  40
      end if
      
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
      oesp  = .true.
      ogeom = .true.
      esp_in(1) =" "
      esp_in(2) =" "
      geom_in(1)=" "
      geom_in(2)=" "
c
      if (.not. rtdb_cput(rtdb, 'task:operation', 1, "dynamics"))
     $    call errquit(pname,0, RTDB_ERR)
c
c     reference and target esp files
c     ------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:fep_esp', 2, fileesp))
     $    oesp = .false.
c
c     reference and target geometries
c     -------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:fep_geom', 2, filegeom))
     $    ogeom=.false.
c
      if((.not.oesp).and.(.not.ogeom))
     > call errquit(pname//"need esp or geometry files ",0, RTDB_ERR)
c
      call mm_system_prefix(prefix)
c
c     coupling parameter for reference and target states
c     --------------------------------------------------
      if (.not.rtdb_get(rtdb,'qmmm:fep_lambda',mt_dbl,2,lambda)) then
         lambda(1) = 0.0d0
         lambda(2) = 1.0d0
      end if
c
c     construct reference and target states using lambda
c     --------------------------------------------------
      if(oesp) then
       esp_in(1)=prefix(1:inp_strlen(prefix))//"-a.esp"
       esp_in(2)=prefix(1:inp_strlen(prefix))//"-b.esp"
       if(master) then
         call qmmm_interp_esp_file(fileesp(1),fileesp(2),
     >                              esp_in(1),lambda(1))
         call qmmm_interp_esp_file(fileesp(1),fileesp(2),
     >                              esp_in(2),lambda(2))
       end if
      end if

      if(ogeom) then
       geom_in(1)=prefix(1:inp_strlen(prefix))//"-a.xyzi"
       geom_in(2)=prefix(1:inp_strlen(prefix))//"-b.xyzi"
       if(master) then
         call qmmm_interp_xyzi_file(filegeom(1),filegeom(2),
     >                              geom_in(1),lambda(1))

         call qmmm_interp_xyzi_file(filegeom(1),filegeom(2),
     >                              geom_in(2),lambda(2))
       end if
      end if
c
      call ga_sync()
c
      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1
c
c
c     trajectory filename
c     -------------------
      call mm_trajectory_filename(filetrj)
c
c     frequency of solvent(nf) and solute(nfs) coord dumps
c     ----------------------------------------------------
      call mm_get_nfc(nf,nfs)
c
c     number of equilibration steps
c     -----------------------------
      call mm_get_mequi(mequi)
c
c     commence the actual run
c     -----------------------
      do in=1,ncycles
c       dynamics task to generate trajectory
c       push reference geometry charges if any
        if(ogeom)  
     >    call mm_set_solute_coord_file(geom_in(1))
        if(oesp) then
          if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $      call errquit('qmmm ti: failed ', 0, RTDB_ERR)
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_in(1)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if


c
c        if (.not.rtdb_get(irtdb,"qmmm:update_geom",
c     >                   mt_log,1,.true.)) 
c     >   call errquit(pname,0, RTDB_ERR)
c        call qmmm_geom_create_full(irtdb)
c        if (.not.rtdb_get(irtdb,"qmmm:update_geom",
c     >                   mt_log,1,.false.)) 
c     >   call errquit(pname,0, RTDB_ERR)
c 

c       note that region will be set in dynamics
        ostatus = qmmm_dynamics(rtdb)
c
        call mm_end()
        call mm_set_nfc(0,0)
        if (.not. rtdb_cput(rtdb,'task:operation', 1,"energy"))
     $      call errquit(pname,0, RTDB_ERR)
        call mm_init(rtdb)
        call qmmm_bq_data_load()
c       disable coord dumps
        call mm_set_nfc(0,0)
c       do the perturbation resampling
        ostatus = qmmm_fep(rtdb,geom_in(2),esp_in(1),
     >                     esp_in(2),filetrj)     
c       restore coord dumps
        call mm_set_nfc(nf,nfs)
c       request to reuse collected fep data for ncycle.gt.1
        if (.not.rtdb_put(rtdb,"qmmm:extend",mt_log,1,.true.)) 
     >      call errquit(pname,0, RTDB_ERR)
c       reinitialize md module
        if(in.ne.ncycles) then
          call mm_end()
        if (.not. rtdb_cput(rtdb,'task:operation', 1,"dynamics"))
     $      call errquit(pname,0, RTDB_ERR)
          call mm_init(rtdb)
        end if
c       disable equlibration for future cycles
        call mm_set_mequi(0)
      end do
c
c     reset any region constraints
      call qmmm_cons_reset()
40    task_qmmm_fep = ostatus
      return

      end

      function task_qmmm_fep_deriv(rtdb)
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
      logical task_qmmm_fep_deriv
      logical ostatus
c     external functions
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
      character*255 filegeom(2)
      character*255 fileesp(2)
      character*255 esp_in(3)
      character*255 geom_in(3)
      double precision lambda(3)
      logical oesp,ogeom
      integer mequi 

      pname = "task_qmmm_fep_deriv"
      master = ga_nodeid().eq.0
      if(qmmm_print_debug()) 
     $  write(*,*) "in "//pname
c
      oesp  = .true.
      ogeom = .true.
      esp_in(1) =" "
      esp_in(2) =" "
      esp_in(3) =" "
      geom_in(1)=" "
      geom_in(2)=" "
      geom_in(3)=" "
c
      if (.not. rtdb_cput(rtdb, 'task:operation', 1, "dynamics"))
     $    call errquit(pname,0, RTDB_ERR)
c
c     reference and target esp files
c     ------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:fep_esp', 2, fileesp))
     $    oesp = .false.
c
c     reference and target geometries
c     -------------------------------
      if (.not. rtdb_cget(rtdb, 'qmmm:fep_geom', 2, filegeom))
     $    ogeom=.false.
c
      if((.not.oesp).and.(.not.ogeom))
     > call errquit(pname//"need esp or geometry files ",0, RTDB_ERR)
c
      call mm_system_prefix(prefix)
c
c     coupling parameter for reference and target states
c     --------------------------------------------------
      if (.not.rtdb_get(rtdb,'qmmm:fep_lambda',mt_dbl,2,lambda)) then
         lambda(1) = 0.0d0
         lambda(2) = 0.5d0
      end if
      lambda(3) = lambda(1)-(lambda(2)-lambda(1))
c
c     construct reference and target states using lambda
c     --------------------------------------------------
      if(oesp) then
       esp_in(1)=prefix(1:inp_strlen(prefix))//"-a.esp"
       esp_in(2)=prefix(1:inp_strlen(prefix))//"-b.esp"
       esp_in(3)=prefix(1:inp_strlen(prefix))//"-c.esp"
       if(master) then
         call qmmm_interp_esp_file(fileesp(1),fileesp(2),
     >                              esp_in(1),lambda(1))
         call qmmm_interp_esp_file(fileesp(1),fileesp(2),
     >                              esp_in(2),lambda(2))
         call qmmm_interp_esp_file(fileesp(1),fileesp(2),
     >                              esp_in(3),lambda(3))
       end if
      end if

      if(ogeom) then
       geom_in(1)=prefix(1:inp_strlen(prefix))//"-a.xyzi"
       geom_in(2)=prefix(1:inp_strlen(prefix))//"-b.xyzi"
       geom_in(3)=prefix(1:inp_strlen(prefix))//"-c.xyzi"
       if(master) then
         call qmmm_interp_xyzi_file(filegeom(1),filegeom(2),
     >                              geom_in(1),lambda(1))
         call qmmm_interp_xyzi_file(filegeom(1),filegeom(2),
     >                              geom_in(2),lambda(2))
         call qmmm_interp_xyzi_file(filegeom(1),filegeom(2),
     >                              geom_in(3),lambda(3))
       end if
      end if
c
      call ga_sync()
c
      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1
c
c
c     trajectory filename
c     -------------------
      call mm_trajectory_filename(filetrj)
c
c     frequency of solvent(nf) and solute(nfs) coord dumps
c     ----------------------------------------------------
      call mm_get_nfc(nf,nfs)
c
c     number of equilibration steps
c     -----------------------------
      call mm_get_mequi(mequi)
c
c     commence the actual run
c     -----------------------
      do in=1,ncycles
c       dynamics task to generate trajectory
c       push reference geometry if any
        if(ogeom)  
     >    call mm_set_solute_coord_file(geom_in(1))
        if(oesp) then
          if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $      call errquit('qmmm ti: failed ', 0, RTDB_ERR)
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_in(1)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if


c
c        if (.not.rtdb_get(irtdb,"qmmm:update_geom",
c     >                   mt_log,1,.true.)) 
c     >   call errquit(pname,0, RTDB_ERR)
c        call qmmm_geom_create_full(irtdb)
c        if (.not.rtdb_get(irtdb,"qmmm:update_geom",
c     >                   mt_log,1,.false.)) 
c     >   call errquit(pname,0, RTDB_ERR)
c 

c       note that region will be set in dynamics
        if(ga_nodeid().eq.0) write(*,*) "GENERATING MD TRAJECTORY"
        ostatus = qmmm_dynamics(rtdb)
c
        call mm_end()
        call mm_set_nfc(0,0)
        if (.not. rtdb_cput(rtdb,'task:operation', 1,"energy"))
     $      call errquit(pname,0, RTDB_ERR)
        call mm_init(rtdb)
        call qmmm_bq_data_load()
c       disable coord dumps
        call mm_set_nfc(0,0)
c       do the perturbation resampling
        if(ga_nodeid().eq.0) then
          write(*,*) "STARTING RESAMPLING OF MD TRAJECTORY"
          write(*,*) "trajectory file", filetrj
        end if
        ostatus = qmmm_fep_deriv(rtdb,geom_in,esp_in,
     >                     filetrj)     
c       restore coord dumps
        call mm_set_nfc(nf,nfs)
c       request to reuse collected fep data for ncycle.gt.1
        if (.not.rtdb_put(rtdb,"qmmm:extend",mt_log,1,.true.)) 
     >      call errquit(pname,0, RTDB_ERR)
c       reinitialize md module
        if(in.ne.ncycles) then
          call mm_end()
        if (.not. rtdb_cput(rtdb,'task:operation', 1,"dynamics"))
     $      call errquit(pname,0, RTDB_ERR)
          call mm_init(rtdb)
        end if
c       disable equlibration for future cycles
        call mm_set_mequi(0)
      end do
c
c     reset any region constraints
      call qmmm_cons_reset()
      task_qmmm_fep_deriv = ostatus
      return

      end
C> @}
cc
      function qmmm_fep_deriv(rtdb,geom_file,esp_file,filetrj)
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
      character*(*) geom_file(3)
      character*(*) esp_file(3)
      character*(*) filetrj
      logical qmmm_fep_deriv
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision stime
      integer i,i0
      character*255 filename
      integer fn_trj,fn
      integer trn(3)
      integer in,nf,nf0
      character*30 pname
      character*255 buffer
      logical master
      integer i_e1,h_e1
      integer i_e2,h_e2
      integer i_e3,h_e3
      integer i_t,h_t
      integer offset 
      double precision de1,de2
      double precision tsum, fsum1,fsum2
      double precision taver,faver1,faver2
      logical oextend
c
      logical ti_geom, ti_esp
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "qmmm_fep_deriv"
c
      qmmm_fep_deriv = .true.
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
      ti_geom = geom_file(1).ne." "

      ti_esp  = esp_file(1).ne." "

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
      if(.not.ma_alloc_get(MT_DBL, nf, 'qmmm ti e3',
     &      h_e3, i_e3) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      nf, MA_ERR)
      call dfill(nf,0,dbl_mb(i_e3),1)
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
      fsum1 = 0.0d0
      fsum2 = 0.0d0
c
     
      if(master) then
       if(oextend) then
       do i=1,10000
         read(fn,*,END=10)  nf0,taver,faver1,de1,faver2,de2
       end do
10     continue
c this was done to get rid off past EOF problem on chinook
       i0=i-1
       write(*,*) "rereading", i0
       call util_flush(6)
       rewind(fn)
       do i=1,i0
         read(fn,*)  nf0,taver,faver1,de1,faver2,de2
         write(*,*) i,nf0,taver,faver1,de1,faver2,de2
       end do
c
       tsum = taver*nf0
       fsum1 = faver1*nf0
       fsum2 = faver2*nf0
       end if
      end if

      call ga_sync()
      if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(1)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
      end if


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
        qmmm_fep_deriv = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e1+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        if(ti_geom)  
     >    call mm_set_solute_coord_file(geom_file(2))
        if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(2)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
        end if

        call md_sp()
        qmmm_fep_deriv = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e2+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        if(ti_geom)  
     >    call mm_set_solute_coord_file(geom_file(3))
        if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(3)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
        end if

        call md_sp()
        qmmm_fep_deriv = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,dbl_mb(i_e3+in-1)))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)



        call ga_sync()
        if(master) write(*,'(A,I6,5F12.6)') 
     >      "ti:",
     >      in,
     >      stime,
     >      dbl_mb(i_t+in-1),
     >      dbl_mb(i_e1+in-1),
     >      dbl_mb(i_e2+in-1),
     >      dbl_mb(i_e3+in-1)

        if(master) call util_flush(6)
        if(ti_esp) then
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(1)))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if

        offset = trn(3)
        if(master) then
          tsum = tsum + dbl_mb(i_t+in-1)
          fsum1 = fsum1 + exp
     >                  (
     >         (
     >         dbl_mb(i_e1+in-1)-
     >         dbl_mb(i_e2+in-1)
     >         )
     >         /(kb_au*dbl_mb(i_t+in-1))
     >                  )
           fsum2 = fsum2 + exp
     >                  (
     >         (
     >         dbl_mb(i_e1+in-1)-
     >         dbl_mb(i_e3+in-1)
     >         )
     >         /(kb_au*dbl_mb(i_t+in-1))
     >                  )
          taver = tsum/(in+nf0)
          faver1 = fsum1/(in+nf0)
          faver2 = fsum2/(in+nf0)
          de1 = -kb_au*taver*log(faver1)
          de2 = -kb_au*taver*log(faver2)
          write(*,66) de1*627.51,de2*627.51
 66       format('current free energy difference ',2F12.6)
        end if
        call ga_sync()

      end do
 23   continue
c
      if(master) then
        write(fn,*) nf0+nf,taver,faver1,de1,faver2,de2
        write(*,*) "Total free energy difference (kcal/mol):",
     +               (de1-de2)*627.51
      end if
c
      if(.not.ma_free_heap(h_t))
     &      call errquit(
     &      pname//'unable to deallocate heap space',
     &      nf, MA_ERR)

      if(.not.ma_free_heap(h_e3))
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
      call errquit('Failed to open trajectory file ',0,0)

 999  continue
      call errquit('Failed to open fep file '//filename,0,0)

      end

      function qmmm_fep_gen(rtdb,
     >                      ngeom,geom_file,
     >                      nesp,esp_file,
     >                      ntr,filetrj,
     >                      ns,tav,de)
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
      integer ngeom
      character*(*) geom_file(ngeom)
      integer nesp
      character*(*) esp_file(nesp)
      integer ntr(3)
      character*(*) filetrj
      integer ns
      double precision tav
      double precision de(*)
      logical qmmm_fep_gen
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      integer n,np
      double precision stime
      integer i
      integer fn_trj
      integer in,nf,nf0
      character*30 pname
      logical master
      integer offset 
      double precision eref
      double precision taver
      double precision tsum
      double precision e(100)
      double precision t
      double precision fsum(100)
c
      logical ti_geom, ti_esp, oextend
c
      logical  mm_read_frame
      external mm_read_frame
      logical  mm_skip_frame
      external mm_skip_frame
c
      pname = "qmmm_fep_gen"
c
      qmmm_fep_gen = .true.
c
      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
c
      master = qmmm_master()
c
      if(master) then
        if(.not.qmmm_get_io_unit(fn_trj)) 
     >   call errquit("cannot get file number",0,0)
        open(unit=fn_trj,file=filetrj,
     +   form='formatted',status='old',err=998)
      end if
c
      ti_geom = ngeom.gt.1
      ti_esp  = nesp.gt.1
      if( (.not.ti_geom).and.(.not.ti_esp))  
     > call errquit(pname//'neither esp or geom file was specified',0,0)
      np = max(ngeom,nesp)-1
c
      oextend = ns.gt.0
      if(oextend) then
        tsum = taver*ns
        do n=1,np
          fsum(n) = exp(-de(n)/(kb_au*tav))*ns
        end do
      else
        tsum = 0.0d0
        do n=1,np
          fsum(n) = 0.0d0
        end do
      end if
c
c     total number of frames to process
c     ---------------------------------
      nf = MAX((ntr(2)-ntr(1)+ntr(3))/ntr(3),0)
      if(master) 
     >   write(*,*) "ti: total number of frames",nf
c
      offset = ntr(1)
      in = 0
      do i=ntr(1),ntr(2),ntr(3)
        in = in + 1
        if(.not.mm_read_frame(fn_trj,offset))
     >      call errquit(pname//'failed to get skip frames',
     >                   0,0)
        call mm_get_temp(t)
        call mm_get_stime(stime)
c
        if(ti_esp) then
          if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $       call errquit('qmmm ti: failed ', 0, RTDB_ERR)
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(n)))
     >        call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if
c
        call md_sp()
        qmmm_fep_gen = qmmm_energy_gradient(rtdb,.false.)
        call qmmm_energy_rtdb_push(rtdb)
        call qmmm_print_energy(rtdb)
        if (.not. rtdb_get(rtdb,'qmmm:energy',
     >      mt_dbl,1,eref))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

        do n=1,np
          if(ti_geom)  
     >      call mm_set_solute_coord_file(geom_file(n))
          if(ti_esp) then
            if(.not.rtdb_put(rtdb,'qmmm:readesp',mt_log,1,.true.))
     $         call errquit('qmmm ti: failed ', 0, RTDB_ERR)
            if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_file(n)))
     >          call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
            call qmmm_esp_reset(rtdb)
          end if
          call md_sp()
          qmmm_fep_gen = qmmm_energy_gradient(rtdb,.false.)
          call qmmm_energy_rtdb_push(rtdb)
          call qmmm_print_energy(rtdb)
          if (.not. rtdb_get(rtdb,'qmmm:energy',
     >        mt_dbl,1,e(n)))
     $       call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        end do
        offset = ntr(3)

        if(master) write(*,*) 
     >      "ti:",
     >      in,
     >      stime,
     >      t,
     >      (e(n),n=1,np)

         tsum = tsum + t
         tav = tsum/(in+ns)
         do n=1,np
           fsum(n) = fsum(n) + 
     +             exp((eref-e(n))/(t*kb_au))
           de(n) = -kb_au*tav*log(fsum(n)/(in+nf0))
         end do
         write(*,*) 'current free energy difference ',
     +        (de(n)*627.51,n=1,np)

      end do
 23   continue

      if(master) then
        close(fn_trj)
      end if

      if(qmmm_print_debug()) 
     $  write(*,*) "in ",pname
      if(master) call util_flush(6)
      return
 998  continue
      call errquit('Failed to open trajectory file ',0,0)
      end

      function qmmm_fep(rtdb,geom_file,esp_ref,esp_pert,filetrj)
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
      character*(*) geom_file
      character*(*) esp_ref
      character*(*) esp_pert
      character*(*) filetrj
      logical qmmm_fep
c     external functions
      logical qmmm_energy_gradient
      external qmmm_energy_gradient 

c     local variables
      double precision stime
      integer i,i0
      character*255 filename
      integer fn_trj,fn
      integer trn(3)
      integer in,nf,nf0
      character*30 pname
      character*255 buffer
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
      ti_geom = geom_file.ne." "

      ti_esp  = esp_pert.ne." "

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
c this was done to get rid off past EOF problem on chinook
       i0=i-1
       write(*,*) "rereading", i0
       call util_flush(6)
       rewind(fn)
       do i=1,i0
         read(fn,*)  nf0,taver,faver,de
         write(*,*)  i,nf0,taver,faver,de
       end do
c
       tsum = taver*nf0
       fsum = faver*nf0
       end if
      end if

      call ga_sync()
      if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_ref))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
      end if


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
     >    call mm_set_solute_coord_file(geom_file)
        if(ti_esp) then
        if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_pert))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
        call qmmm_esp_reset(rtdb)
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
        if(ti_esp) then
          if(.not.rtdb_cput(rtdb,"qmmm:espfilename",1,esp_ref))
     >      call errquit('qmmm: failed set espfile', 0, RTDB_ERR)
          call qmmm_esp_reset(rtdb)
        end if

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
        write(*,*) "Total free energy difference (kcal/mol):",
     +               de*627.51
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
      call errquit('Failed to open trajectory file ',0,0)

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

