c**************************************************
c     Main driver for all qmmm optimization tasks *
c**************************************************
      function task_qmmm_optimize(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      task_qmmm_optimize = status
c     
      end
c
      function qmmm_mixed_optimization(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      character*30 region(3)
      character*30 method(3)
      integer maxiter(3)
      integer nmethod
      integer nmaxiter
      integer nregion
      integer ncycles
c
      pname = "qmmm_mixed_optimization"
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if(qmmm_print_debug())
     >  write(*,*) "in " //pname, ga_nodeid()
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
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
c
         status = qmmm_optimizer(rtdb,region(k),method(k),
     >                         maxiter(k))
111    continue
       end do
      end do

      if(ga_nodeid().eq.0) then
         write(luout,10) 
         write(luout,12) 
         write(luout,15) 
         write(luout,12) 
15     format("@  OPTIMIZATION FINISHED ")
      end if

      call qmmm_cons_reset()
      call qmmm_bq_data_update_active()
      call qmmm_set_interface("qm")
      call mm_task_sp()
      status = qmmm_energy_gradient(rtdb,.false.) 
      call qmmm_energy_rtdb_push(rtdb)
      call qmmm_print_energy3(rtdb)

      status = .true.         
c     
110   continue
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c     
      call mm_write_restart()
      qmmm_mixed_optimization = status
c     
      end
c
      function driver1(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      write(*,*) "debug 3"
      write(*,*) "debug 4"
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

      return

 133  call errquit(pname//'error opening/closing xyz file',0, GEOM_ERR)
     
      end
c
      subroutine qmmm_optimization_method(method,region)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
c
      if(region.eq."qm") then
        method = "bfgs"
      else if (region.eq."qmlink") then
        method = "bfgs"
      else if (region.eq."mm") then
        method = "sd"
      else if (region.eq."solvent") then
        method = "sd"
      else
        method = "lbfgs"
      end if

      return

      end
c
      subroutine qmmm_optimization_maxiter(maxiter,region)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      function qmmm_optimizer(rtdb,region,method,maxiter)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      character*(*) region,method
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
      call qmmm_cons_set(rtdb,"fix","solute")
      call qmmm_cons_set(rtdb,"fix","solvent")
      call qmmm_cons_set(rtdb,"free",region)
c
      if (.not. rtdb_get(rtdb, 'qmmm:esp_density', mt_log, 
     $                   1, oesp))
     $     oesp = .false.
c
      if (.not. rtdb_get(rtdb, 'qmmm:frozen_density', mt_log, 
     $                   1, ofrozen))
     $     ofrozen = .false.
c
      write(*,*) "oesp,ofrozen",oesp,ofrozen
      call qmmm_bq_data_update_active()
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
         call qmmm_geom_init_aux(rtdb)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter))
     $        call errquit(pname,0, RTDB_ERR)
         call qmmm_set_interface("qm")
         status = driver(rtdb)
         ignore = rtdb_delete(rtdb,'geometry')
         call qmmm_geom_aux_end(rtdb)
       else if(method.eq."lbfgs") then
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
      function qmmm_pes(irtdb,optimizer)
c     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
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
      double precision energy
      double precision rb(2)
      integer np
      integer ip
      double precision dr
      double precision r
      logical ofile

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
      if (.not. rtdb_get(irtdb,'qmmm:pes_npoints',mt_int,
     +                    1,np))
     + call errquit(pname//'Failed to get pes points',
     + 1, RTDB_ERR)

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
c     set increments and initial value
c     --------------------------------
      dr = (rb(2)-rb(1))/(np-1)
      r  = rb(1)
c
c     actually do the scan
c     --------------------
      do ip=1,np
        if (.not. rtdb_put(irtdb,'qmmm:r_spring',mt_dbl,
     +                     1,r))
     +  call errquit(pname//'Failed to set r_spring',
     +  1, RTDB_ERR)

c       construct restart file name  
        filename0 = ' '
        call mm_restart_filename(filename0)
c        write(filename0,14)  prefix(1:inp_strlen(prefix)),
c     $                       ip
        filename = ' '
        write(filename,13) dir(1:inp_strlen(dir)), 
     $     prefix(1:inp_strlen(prefix)),ip

c       load restart file if it exists      
        if(ga_nodeid().eq.0) then 
          inquire(file=filename,exist=ofile)
          write(*,*) "pes found file",
     >               filename(1:inp_strlen(filename)),
     >               ofile
          if(ofile) 
     >     call util_file_copy(filename(1:inp_strlen(filename)),
     >                        filename0(1:inp_strlen(filename0)))
          write(*,*) "pes point ",ip,r,dr,
     >             filename0(1:inp_strlen(filename0))
        end if
        call ga_brdcst(msg_qmmm_misc, ofile, mitob(1), 0)
        call ga_sync()
        if(ofile) then
          call mm_reload_rst(irtdb,filename0)
        end if
        status = qmmm_mixed_optimization(irtdb)
        if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $       call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        write(*,*) "pes energy",ip,r,energy
        call mm_write_restart_named(filename)
         r = r + dr
      end do

 13      format(a,'/',a,'_pes',i3.3,'.rst')
 14      format(a,'_pes',i3.3,'.rst')

      return
      end

c      function qmmm_pes(irtdb,optimizer)
cc     $Id: task_qmmm_optimize.F,v 1.36 2005-09-02 20:28:27 marat Exp $
c      implicit none
c#include "errquit.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "geom.fh"
c#include "global.fh"
c#include "stdio.fh"
c#include "util.fh"
c#include "qmmm.fh"
c#include "inp.fh"
c#include "msgids.fh"
c#include "tcgmsg.fh"
c      integer irtdb
c      logical optimizer
c      external optimizer
c      logical qmmm_pes
c      double precision cpu, wall
cc     
c      logical  status
cc     
c      character*32 pname
cc     
c      logical   driver
c      external  driver
c      logical ignore
c      character*255 filename, dir,xyz,prefix
c      character*255 filename0
c      double precision energy
c      double precision rb(2)
c      integer np
c      integer ip
c      double precision dr
c      double precision r
c      logical ofile
c
c      pname = "qmmm_pes"
c
cc
cc     get lower/upper  limits for pes scan
cc     ------------------------------------
c      if (.not. rtdb_get(irtdb,'qmmm:pes_bounds',mt_dbl,
c     +                    2,rb))
c     + call errquit(pname//'Failed to get pes bounds',
c     + 1, RTDB_ERR)
c
cc
cc     get number of points for pes scan
cc     ------------------------------------
c      if (.not. rtdb_get(irtdb,'qmmm:pes_npoints',mt_int,
c     +                    1,np))
c     + call errquit(pname//'Failed to get pes points',
c     + 1, RTDB_ERR)
c
cc
cc     directory name to dump restart files
cc     ------------------------------------
c      dir      = ' '
c      call util_directory_name(dir, .false., 0)
cc
cc     file prefix from mm module
cc     --------------------------
c      call mm_system_prefix(prefix)
cc      write(*,*) "pes, mm_system_prefix", prefix(1:inp_strlen(prefix))
cc 
cc     set increments and initial value
cc     --------------------------------
c      dr = (rb(2)-rb(1))/(np-1)
c      r  = rb(1)
cc
cc     actually do the scan
cc     --------------------
c      do ip=1,np
c        if (.not. rtdb_put(irtdb,'qmmm:r_spring',mt_dbl,
c     +                     1,r))
c     +  call errquit(pname//'Failed to set r_spring',
c     +  1, RTDB_ERR)
c
cc       construct restart file name  
c        filename0 = ' '
c        call mm_restart_filename(filename0)
cc        write(filename0,14)  prefix(1:inp_strlen(prefix)),
cc     $                       ip
c        filename = ' '
c        write(filename,13) dir(1:inp_strlen(dir)), 
c     $     prefix(1:inp_strlen(prefix)),ip
c
cc       load restart file if it exists      
c        if(ga_nodeid().eq.0) then 
c          inquire(file=filename,exist=ofile)
c          write(*,*) "pes found file",
c     >               filename(1:inp_strlen(filename)),
c     >               ofile
c          if(ofile) 
c     >     call util_file_copy(filename(1:inp_strlen(filename)),
c     >                        filename0(1:inp_strlen(filename0)))
c          write(*,*) "pes point ",ip,r,dr,
c     >             filename0(1:inp_strlen(filename0))
c        end if
c        call ga_brdcst(msg_qmmm_misc, ofile, mitob(1), 0)
c        call ga_sync()
c        if(ofile) then
c          call mm_reload_rst(irtdb,filename0)
c        end if
cc       recreate active geometry
c        call qmmm_geom_create_active(irtdb)
c        status = optimizer(irtdb)
c        if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
c     $       call errquit('qmmm: failed get energy', 0, RTDB_ERR)
c        write(*,*) "pes energy",ip,r,energy
c        call mm_write_restart_named(filename)
c         r = r + dr
c      end do
c
c 13      format(a,'/',a,'_pes',i3.3,'.rst')
c 14      format(a,'_pes',i3.3,'.rst')
c
c      return
c      end
