      function task_qmmm_optimize(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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
      logical ouselbfgs
      character*32 theory
      character*32 optimization
      character*84 tag
      character*32 pname
      character*48 string
c     
      logical   driver
      external  driver
c
      logical   driver1
      external  driver1
c
      logical   qmmm_lbfgsb_driver
      external  qmmm_lbfgsb_driver
c
      logical   qmmm_mixed_optimization
      external  qmmm_mixed_optimization
c
      logical   qmmm_optimizer
      external  qmmm_optimizer
c
      logical   qmmm_tandem_optimization
      external  qmmm_tandem_optimization
      logical   qmmm_pes
      external  qmmm_pes
      logical   lbfgsb_driver
      external  lbfgsb_driver
      integer resid(1),i,k
      logical ignore
      logical omixed
      logical ofrozen
      logical divcon
      logical opes
c
      character*30 region(3)
      character*30 method(3)
      integer maxiter(3)
      integer nmethod
      integer nmaxiter
      integer nregion
      integer ncycles
c
      pname = "task_qmmm_optimize"
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
        method(i) = " "
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
        maxiter(i) = -1
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
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0,
     &       RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'qmmm:pes', mt_log, 
     $                   1, opes))
     $     opes = .false.

      if (.not. rtdb_get(rtdb, 'driver:divcon', mt_log, 
     $                   1, divcon))
     $     divcon = .false.

      if (.not. rtdb_get(rtdb, 'qmmm:frozen_density', mt_log, 
     $                   1, ofrozen))
     $     ofrozen = .false.


      if (.not. rtdb_get(rtdb, 'qmmm:mixed', mt_log, 
     $                   1, omixed))
     $     omixed = .false.


      status = .true.
c
c     pes task
c     --------
      if(opes) then
        if(ouselbfgs) then
          status=qmmm_pes(rtdb,qmmm_lbfgsb_driver)
        else
          status=qmmm_pes(rtdb,driver)
        end if
        go to 110
      end if
      if(omixed) then
        if(qmmm_master()) write(*,*) "mixed optimization"
        status=qmmm_mixed_optimization(rtdb)
        go to 110
      end if
c
c     print informational header
c     --------------------------
      if(qmmm_master()) then
        string = "QM/MM Optimization"
        call util_print_centered(luout,
     >      string(1:inp_strlen(string)), 
     >      0,.true.)

        string = ' '
        do k=1,nregion
          if(maxiter(k).ne.0)
     >     write(luout,13) k,
     >         region(k),method(k),maxiter(k)
        end do 
13      format("region",i1,":  ",a10,"with ",a10,"maxiter =",i3)
        write(luout,'(//)')
      end if

      do i=1,ncycles
       do k=1,nregion
c
       if(maxiter(k).eq.0) goto 111
c
       status = qmmm_optimizer(rtdb,region(k),method(k),
     >                         maxiter(k),ofrozen)
111    continue
       end do
       end do

c
c
c      if (optimization.eq.'all') then
c         if(qmmm_print_debug())
c     >      write(*,*) "going into md_em"
c         call md_em()
c      else if (optimization.eq.'tandem') then
c         status =  qmmm_mixed_optimization(rtdb)
c      else if (optimization.eq.'mm') then
c         if(qmmm_print_debug())
c     >      write(*,*) "going into md_em"
c         call md_em()
c      else if (optimization.eq.'qm') then
c         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int, 1, 0))
c     $        call errquit('task_optimize: rtdb?',0, RTDB_ERR)
c         if(qmmm_print_debug())
c     >      write(*,*) "going into driver"
c          if(ouselbfgs) then
c            status =  qmmm_lbfgsb_driver(rtdb)
c          else if(divcon) then
c            status = driver1(rtdb)
c          else
c            status = driver(rtdb)
c          end if
c         call mm_print_info()
c
c        forcing driver to never fail
c        ---------------------------
         status = .true.         
c      endif
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
c      call util_print_pop()
c     
      end
c
      function driver1(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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

      function qmmm_tandem_optimization(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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
      logical qmmm_tandem_optimization
c
      logical  status
c     
      character*32 pname
c     
      logical   driver
      external  driver
      logical   qmmm_energy_gradient
      external  qmmm_energy_gradient
      integer   i
      logical ignore
      double precision energy
      double precision grms_quantum
      integer ncycles

      pname = "qmmm_tandem_optimization"

      qmmm_tandem_optimization = .true.

      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1

      write(*,*) "tandem optimization",ncycles

      do i=1,ncycles
       call qmmm_cons_reset()
       call qmmm_cons_set(rtdb,"fix","mm_solute")
       ignore = rtdb_delete(rtdb,'geometry')
       call qmmm_geom_init_aux(rtdb)
       call qmmm_set_interface("qm")
       status = driver(rtdb)
       call mm_get_rms_solute_force_gen(mm_quant,grms_quantum)

       call qmmm_cons_reset()
       call qmmm_cons_set(rtdb,"fix","qmlink")
       call qmmm_set_interface("mm")
       if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $                   mt_log,1,.true.))
     $     call errquit(pname,0, RTDB_ERR)
       call md_em()
       if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $                   mt_log,1,.false.))
     $     call errquit(pname,0, RTDB_ERR)
       call mm_free_solute(mm_quant)

       status = qmmm_energy_gradient(rtdb,.false.) 
       call qmmm_energy_rtdb_push(rtdb)
       if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
       write(*,*) "qmmm_energy",i,energy,grms_quantum

      end do
      return

      end

      function qmmm_mixed_optimization(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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
      logical qmmm_mixed_optimization
c
      logical  status
c     
      character*32 pname
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
      integer maxiter(3)
      character*30 region(3)
      character*30 interface(3)
      character*30 method(3)
      character*48 string

      pname = "qmmm_mixed_optimization"

      qmmm_mixed_optimization = .true.

      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1

      if (.not. rtdb_get(rtdb,'qmmm:maxiter',mt_int,3,maxiter)) then
         maxiter(1) = 20
         maxiter(2) = 100
         maxiter(3) = 100
      end if

      region(1) = "qmlink"
      region(2) = "mm_solute"
      region(3) = "solvent"

      interface(1) = "qm"
      interface(2) = "qm"
      interface(3) = "mm"

      method(1) = "bfgs"
      method(2) = "lbfgs"
      method(3) = "sd"
c
c     print informational header
c     --------------------------
      if(qmmm_master()) then
        string = "Mixed Mode QM/MM Optimization"
        call util_print_centered(luout,
     >      string(1:inp_strlen(string)), 
     >      0,.true.)

        string = ' '
        do k=1,3
          if(maxiter(k).ne.0)
     >     write(luout,13) k,
     >         region(k),method(k),maxiter(k)
        end do 
13      format("region",i1,":  ",a10,"with ",a10,"maxiter =",i3)
        write(luout,'(//)')
      end if

      do i=1,ncycles
       do k=1,3
c
       if(maxiter(k).eq.0) goto 111
c
       if(qmmm_master()) then
         write(string,14) region(k),
     >                    method(k)
14       format("Optimizing " ,a10," region with ",a10)
         call util_print_centered(luout,
     >        string(1:inp_strlen(string)),
     >        30,.true.)

       end if

       call qmmm_cons_reset()
       call qmmm_cons_set(rtdb,"fix","solute")
       call qmmm_cons_set(rtdb,"fix","solvent")
c
       call qmmm_cons_set(rtdb,"free",region(k))
c
       if(method(k).eq."bfgs") then
         call qmmm_geom_init_aux(rtdb)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter(k)))
     $        call errquit(pname,0, RTDB_ERR)
         call qmmm_set_interface("qm")
         if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $       mt_log,1,.false.))
     $      call errquit(pname,0, RTDB_ERR)
         status = driver(rtdb)
         ignore = rtdb_delete(rtdb,'geometry')
         call qmmm_geom_aux_end(rtdb)
       else if(method(k).eq."lbfgs") then
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter(k)))
     $        call errquit(pname,0, RTDB_ERR)
         call qmmm_set_interface("qm")
         if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $       mt_log,1,.true.))
     $      call errquit(pname,0, RTDB_ERR)
         status =  qmmm_lbfgsb_driver(rtdb)
       else if(method(k).eq."sd") then
         call qmmm_set_interface("mm")
         if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $       mt_log,1,.true.))
     $      call errquit(pname,0, RTDB_ERR)
         if(.not.rtdb_put(rtdb,'md:msdit',mt_int,1,maxiter(k)))
     +      call errquit(pname,0,RTDB_ERR)
         call md_em()
       else
         call errquit(pname//'unknown method',0, RTDB_ERR)
       end if
111    continue
       end do
       call qmmm_set_interface("qm")
       call qmmm_cons_reset()
       if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $       mt_log,1,.false.))
     $      call errquit(pname,0, RTDB_ERR)
       status = task_qmmm_gradient(rtdb) 
       if (.not. rtdb_get(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)
       call mm_get_rms_solute_force_gen(mm_quant_link,grms_quantum)
       call mm_get_rms_solute_force_gen(mm_classical,grms_classical)
       call mm_get_rms_solvent_force_gen(mm_active_w,grms_solvent)
c
       if(qmmm_master()) 
     $  write(*,*) "@@",i,energy,grms_quantum,grms_classical,
     $                          grms_solvent


       end do

      return

      end

      function qmmm_optimizer(rtdb,region,method,maxiter,rhofix)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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
      logical rhofix 
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

      pname = "qmmm_optimizer"

      qmmm_optimizer = .true.

      if (.not. rtdb_put(rtdb,'qmmm:frozen_density',
     $       mt_log,1,rhofix))
     $      call errquit(pname,0, RTDB_ERR)
c
      if(method.eq." ") then
        if(region.eq."qm") then
          method = "bfgs"
        else if (region.eq."qmlink") then
          method = "bfgs"
        else if (region.eq."solvent") then
          method = "sd"
        else
          method = "lbfgs"
        end if
      endif

c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set(rtdb,"fix","solute")
      call qmmm_cons_set(rtdb,"fix","solvent")
      call qmmm_cons_set(rtdb,"free",region)
c
       if(qmmm_master()) then
         write(string,14) region,
     >                    method
14       format("Optimizing " ,a10," region with ",a10)
         call util_print_centered(luout,
     >        string(1:inp_strlen(string)),
     >        30,.true.)

       end if

c
c     perform an optimization as requested
c     ------------------------------------
      if(method.eq."bfgs") then
         call qmmm_geom_init_aux(rtdb)
         if(maxiter.ge.0) then
           if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter))
     $        call errquit(pname,0, RTDB_ERR)
         end if
         call qmmm_set_interface("qm")
         status = driver(rtdb)
         ignore = rtdb_delete(rtdb,'geometry')
         call qmmm_geom_aux_end(rtdb)
       else if(method.eq."lbfgs") then
         if(maxiter.ge.0) then
           if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,maxiter))
     $        call errquit(pname,0, RTDB_ERR)
         end if
         call qmmm_set_interface("qm")
         status =  qmmm_lbfgsb_driver(rtdb)
       else if(method.eq."sd") then
         call qmmm_set_interface("mm")
         if(maxiter.ge.0) then
           if(.not.rtdb_put(rtdb,'md:msdit',mt_int,1,maxiter))
     +      call errquit(pname,0,RTDB_ERR)
         end if
         call md_em()
       else
         call errquit(pname//'unknown method',0, RTDB_ERR)
       end if


      return

      end

      function qmmm_pes(irtdb,optimizer)
c     $Id: task_qmmm_optimize.F,v 1.28 2005-08-25 23:54:40 marat Exp $
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
      logical   driver
      external  driver
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
c       recreate active geometry
        call qmmm_geom_create_active(irtdb)
        status = optimizer(irtdb)
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
