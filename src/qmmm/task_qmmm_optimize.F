c**************************************************
c     Main driver for all qmmm optimization tasks *
c**************************************************
C> \ingroup task
C> @{
      function task_qmmm_optimize(rtdb)
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
      integer rtdb
      logical task_qmmm_optimize
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
      character*30 operation
c     
      logical   qmmm_mixed_optimization
      external  qmmm_mixed_optimization
c
      logical   qmmm_pes
      external  qmmm_pes
c
      logical divcon
      logical opes
c
      pname = "task_qmmm_optimize"
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if(qmmm_print_debug())
     >  write(*,*) "in " //pname, ga_nodeid()

      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0,
     &       RTDB_ERR)

c
c     figure out the optimization task
c     --------------------------------
      if (.not. rtdb_get(rtdb, 'qmmm:pes', mt_log, 
     $                   1, opes))
     $     opes = .false.

      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '

      if (inp_compare(.false., 'saddle', operation)) then
        if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,1))
     $        call errquit('task_saddle: rtdb corrupt?', 0, RTDB_ERR)
      end if

      if(qmmm_print_debug())
     >  write(*,*) "opes",opes


      if (.not. rtdb_get(rtdb, 'driver:divcon', mt_log, 
     $                   1, divcon))
     $     divcon = .false.


      status = .true.
c
c     pes task
c     --------
      if(opes) then
          status=qmmm_pes(rtdb)
        go to 110
      end if
c
c     default multiregion optimization
c     --------------------------------
      status = qmmm_mixed_optimization(rtdb)
c
110   continue
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_optimize: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_optimize: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_optimize: failed to set status',0,
     &       RTDB_ERR)
c     
      call mm_write_restart()
      call mm_write_restart_default()
      task_qmmm_optimize = status
c     
      end
C> @}
c
      function qmmm_mixed_optimization(rtdb)
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
      integer rtdb
      logical qmmm_mixed_optimization
      double precision cpu, wall
c     
      logical  status
c     
      integer   i,k
      character*84 tag
      character*32 pname
      character*48 string
c     
      logical   qmmm_optimizer
      external  qmmm_optimizer
c
      logical   qmmm_energy_gradient
      external  qmmm_energy_gradient
c
      logical oesp
c
      logical oconv
      character*30 region(3)
      character*30 method(3)
      character*30 xyz(3)
      integer maxiter(3)
      integer nmethod
      integer nxyz
      integer nmaxiter
      integer nregion
      integer ncycles
      character*60 xyzfile,xyztmp
      double precision enew,eold,etol
c
      character*255 rfile
      character*255 stub
c
      pname = "qmmm_mixed_optimization"
c
      qmmm_mixed_optimization = .true.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if(qmmm_print_debug())
     >  write(*,*) "in " //pname, ga_nodeid()
c
c     region definitions
c     ------------------
      tag ="qmmm:convergence"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag)),
     >                 mt_dbl,1,etol)) 
     >      etol = 1.0d-4
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
c     method definitions if any
c     -------------------------
      tag ="qmmm:method"
      do i=1,nregion
        call qmmm_optimization_method(method(i),
     >                                region(i))
      end do
      if (rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nmethod)) then
        if(nmethod.gt.3) 
     >       call errquit(pname//"too many methods",0,0)
        if (.not.rtdb_cget(rtdb,tag,nmethod,method)) 
     >       call errquit(pname//tag,0,RTDB_ERR)
      end if
c
c     xyz files
c     ---------
      tag ="qmmm:xyz"
      do i=1,nregion
        xyz(i) = " "
      end do
      if (rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nxyz)) then
        if(nxyz.gt.3) 
     >       call errquit(pname//"too many xyz",0,0)
        if (.not.rtdb_cget(rtdb,tag,nxyz,xyz)) 
     >       call errquit(pname//tag,0,RTDB_ERR)
      end if
c
c
c     max iterations if any
c     ---------------------
      tag ="qmmm:maxiter"
      do i=1,nregion
        call qmmm_optimization_maxiter(maxiter(i),
     >                                region(i))
      end do
      if (rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nmaxiter)) then
        if(nmaxiter.gt.3) 
     >       call errquit(pname//"too many maxiters",0,0)
        if (.not.rtdb_get(rtdb,tag,mt_int,nmaxiter,maxiter)) 
     >       call errquit(pname//tag,0,RTDB_ERR)
      end if
c
      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1

      if (.not. rtdb_get(rtdb, 'qmmm:esp_density', mt_log, 
     $                   1, oesp))
     $     oesp = .false.


      status = .true.
c
c     print informational header
c     --------------------------
      if(qmmm_master()) then
        write(luout,10) 
        write(luout,12)
        write(luout,12) " QM/MM Multiple Region Optimization"

        do k=1,nregion
            write(luout,13) k,region(k),method(k),maxiter(k)
        end do 
        write(luout,12)
10      format("@",48("-"))
12      format("@",a)
13      format("@ region",i1,":  ",a10,"with ",a10,"maxiter =",i3)
        write(luout,12)
        write(luout,10)
      end if

c
c     start the optimization cycle
c     ----------------------------
      eold = 0.0d0
      do i=1,ncycles
       if(ga_nodeid().eq.0) then
         write(luout,14) i
         write(luout,12)
14     format("@ ncycle =",I3)
       end if
       if(oesp) call qmmm_esp_reset(rtdb)
       do k=1,nregion
c
         if(maxiter(k).eq.0) goto 111
         xyzfile = " "
         if(xyz(k).ne." ") then
           xyztmp  = xyz(k)
           write(xyzfile,'(A,"-",I3.3)') 
     >       xyztmp(1:inp_strlen(xyztmp)),i 
         end if

         status = qmmm_optimizer(rtdb,region(k),method(k),
     >                         xyzfile,maxiter(k))
111    continue
       end do
c
c     evaluate energy each cycle
c     updating pair lists and periodic
c     conditions at the same time
c     ---------------------------
      call qmmm_cons_reset()
      call  qmmm_bq_data_reload()
      call qmmm_set_interface("qm")
      call md_sp()
      status = qmmm_energy_gradient(rtdb,.false.) 
      call qmmm_energy_rtdb_push(rtdb)
      call qmmm_print_energy3(rtdb)
      call qmmm_print_spring_energy(rtdb)
      if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,enew))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
c
c      save restart file
c      -----------------
       stub = " "
       rfile = " "
       write(stub,'(i3.3,".rst")') i
       call util_file_name(stub(1:inp_strlen(stub)),.false.,
     >                     .false.,rfile)
       call mm_write_restart_named(rfile,.false.)
c
c      test for convergence
c      --------------------
       oconv = i.gt.1.and.
     >         abs(enew-eold).lt.etol
       if(oconv) goto 109
       eold = enew
      end do
109   continue
      if(ga_nodeid().eq.0) then
         write(luout,10) 
         write(luout,12) 
         write(luout,15) 
         write(luout,12) 
15     format("@  OPTIMIZATION FINISHED ")
      end if

      status = .true.         
c     
110   continue
c--- commenting this out because energy
c--- is already calculated during each cycle
cc
cc     final energy evaluation
cc     ---------------------------
c      call qmmm_cons_reset()
c      call qmmm_bq_data_update_active()
c      call qmmm_set_interface("qm")
c      call md_sp()
c      status = qmmm_energy_gradient(rtdb,.false.) 
c      call qmmm_energy_rtdb_push(rtdb)
c      call qmmm_print_energy3(rtdb)
cc
c----------------------------------------------
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c     
      call mm_write_restart()
      qmmm_mixed_optimization = status
c     
      end
c
      function driver1(rtdb)
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
      integer rtdb
      logical driver1
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical   driver
      external  driver
      logical   qmmm_energy_gradient
      external  qmmm_energy_gradient
      integer nr,resid(2),i,k,resall(3)
      logical ignore
      character*255 filename, dir,xyz
      double precision energy

      xyz = "test"
      nr = 1
      resall(1)=1
      resall(2)=2
      resall(3)=3

      pname = "driver1"

      do i=1,10
      do k=1,58,2
      resid(1) = k
      resid(2) = k+1
      call qmmm_geom_end()
      ignore = rtdb_delete(rtdb,'geometry')
      call mm_free_solute_res(2,resid)
      call qmmm_geom_init(rtdb,.true.)
      status = driver(rtdb)
      call mm_write_restart()
         dir      = ' '
         filename = ' ' 
         call util_directory_name(dir, .false., 0)
         write(filename,13) dir(1:inp_strlen(dir)), 
     $        xyz(1:inp_strlen(xyz)),i,k
 13      format(a,'/',a,'-',i3.3,'-',i3.3,'.pdb')
         open(88,file=filename,form='formatted',status='unknown',
     $        err=133)

         call mm_print_system_file(88)
         write(*,*) "filename",filename(1:inp_strlen(filename))
         close(88,status='keep',err=133)
c
c      call mm_free_solute_res(3,resall)
c      call mm_task_sp()
c      status = qmmm_energy_gradient(rtdb,.false.) 
c      call qmmm_energy_rtdb_push(rtdb)
c      if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,energy))
c     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
c      write(*,*) "qmmm_energy",i,k,energy
      end do
      end do
      call mm_print_info()

      driver1 = status

      return

 133  call errquit(pname//'error opening/closing xyz file',0, GEOM_ERR)
     
      end
c
      subroutine qmmm_optimization_method(method,region)
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
#include "mm_utils.fh"

      character*(*) region,method
      character*(40) buffer
c
      if (region.eq."qm") then
        method = "bfgs"
      else if (region.eq."qmlink") then
        method = "bfgs"
      else if (region.eq."mm_solute") then
        method = "lbfgs"
      else if (region.eq."mm") then
        method = "sd"
      else if (region.eq."solute") then
        method = "sd"
      else if (region.eq."solvent") then
        method = "sd"
      else if (region.eq."all") then
        method = "sd"
      else
        buffer = region
        call errquit("Unknown region "//buffer,0,0)
      end if
c
c      if (region.eq."all") then
c        method = "sd"
c      else if (region.eq."solvent") then
c        method = "sd"
c      else if (region.eq."solvent") then
c        method = "sd"
c      else
c        method = "lbfgs"
c      end if
c
      return

      end
c
      subroutine qmmm_optimization_maxiter(maxiter,region)
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
#include "mm_utils.fh"

      integer       maxiter
      character*(*) region
c
      if(region.eq."qm") then
        maxiter = 20
      else if (region.eq."qmlink") then
        maxiter = 20
      else if (region.eq."mm") then
        maxiter = 100
      else if (region.eq."solvent") then
        maxiter = 100
      else
        maxiter = 50
      end if

      return

      end
c
      function qmmm_optimizer(rtdb,region,method,xyz,maxiter)
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
#include "mm_utils.fh"

      integer rtdb
      logical qmmm_optimizer
      character*(*) region,method,xyz
      integer maxiter
c
      logical  status
c     
      character*32 pname
      character*48 string
c     
      logical   driver
      external  driver
      logical   qmmm_lbfgsb_driver
      external  qmmm_lbfgsb_driver
      logical   qmmm_energy_gradient
      external  qmmm_energy_gradient
      logical task_qmmm_gradient
      external task_qmmm_gradient
      integer   i,k
      logical ignore
      double precision energy
      double precision grms_quantum
      double precision grms_classical
      double precision grms_solvent
      integer ncycles
      logical oesp,ofrozen
      integer nact

      pname = "qmmm_optimizer"

      qmmm_optimizer = .true.
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region)
c
      if (.not. rtdb_get(rtdb, 'qmmm:esp_density', mt_log, 
     $                   1, oesp))
     $     oesp = .false.
c
      if (.not. rtdb_get(rtdb, 'qmmm:frozen_density', mt_log, 
     $                   1, ofrozen))
     $     ofrozen = .false.
c
c      call qmmm_bq_data_update_active()
c
      call qmmm_cons_get_nacts_qm(nact)
c
      oesp = oesp.and.(nact.eq.0)
      ofrozen = ofrozen.and.(nact.eq.0)
c
      if(qmmm_master()) then
        if(oesp) then
         write(luout,14) region,
     >                    method
        else if (ofrozen) then
         write(luout,15) region,
     >                    method
        else
         write(luout,16) region,
     >                    method
        end if

      end if

14       format("@",/,"@ Optimizing " ,a10," region with ",a6,
     >          /,"@ and esp fitted charges for the qm region",
     >           /,"@")
15       format("@",/,"@ Optimizing " ,a10," region with ",a6,
     >          /,"@ and static density for the qm region",
     >          /,"@")
16       format("@",/,"@ Optimizing " ,a10," region with ",a6,
     >          /,"@")
c
c     perform an optimization as requested
c     ------------------------------------
      if(method.eq."bfgs") then
         call qmmm_geom_init_aux(rtdb,"geometry")
         if (.not. rtdb_cput(rtdb,'driver:xyz',1,xyz))
     $        call errquit(pname,0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter))
     $        call errquit(pname,0, RTDB_ERR)
         call qmmm_set_interface("qm")
         status = driver(rtdb)
         ignore = rtdb_delete(rtdb,'geometry')
         call qmmm_geom_aux_end(rtdb)
       else if(method.eq."lbfgs") then
         if (.not. rtdb_cput(rtdb,'driver:xyz',1,xyz))
     $        call errquit(pname,0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter))
     $        call errquit(pname,0, RTDB_ERR)
         call qmmm_set_interface("qm")
         status =  qmmm_lbfgsb_driver(rtdb)
       else if(method.eq."sd") then
         call qmmm_set_interface("mm")
         call mm_set_sd_maxiter(maxiter)
         call md_em()
       else
         call errquit(pname//'unknown method',0, RTDB_ERR)
       end if


      return

      end
c
      function qmmm_pes(irtdb)
c      function qmmm_pes(irtdb,optimizer)
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
      logical optimizer
      external optimizer
      logical qmmm_pes
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical   qmmm_mixed_optimization
      external  qmmm_mixed_optimization
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

      pname = "qmmm_pes"

c
c     get lower/upper  limits for pes scan
c     ------------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_bounds',mt_dbl,
     +                    2,rb))
     + call errquit(pname//'Failed to get pes bounds',
     + 1, RTDB_ERR)

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
c      write(*,*) "pes, mm_system_prefix", prefix(1:inp_strlen(prefix))
c
c     construct movecs file name  
      movecs0 = ' '
      call util_file_name("movecs",.false.,.false.,movecs0)        
      ignore = rtdb_delete(irtdb,'dft:input vectors')
      ignore = rtdb_delete(irtdb,'scf:input vectors')
c 
c     set increments and initial value
c     --------------------------------
      dr = (rb(2)-rb(1))/(np-1)
      r  = rb(1)

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
c     actually do the scan
c     --------------------
      do ip=i1,i2,i3
        r  = rb(1)+dr*(ip-1) 
        if(ga_nodeid().eq.0) then
          write(*,15) ip,r
        end if
        if (.not. rtdb_put(irtdb,'qmmm:r_spring',mt_dbl,
     +                     1,r))
     +  call errquit(pname//'Failed to set r_spring',
     +  1, RTDB_ERR)

c       construct restart file name  
        filename0 = ' '
        call mm_restart_filename(filename0)
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
        status = qmmm_mixed_optimization(irtdb)
        if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $       call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        if(ga_nodeid().eq.0) 
     >     write(*,16) energy,r
        call mm_write_restart_named(filename,.false.)
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
      end do
      qmmm_pes = status

 13      format(a,'/',a,'_pes',i3.3,'.rst')
 14      format(a,'_pes',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)
 17      format(a,'/',a,'_pes',i3.3,'.movecs')

      return
      end


