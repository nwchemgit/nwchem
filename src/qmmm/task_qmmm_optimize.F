      logical function task_qmmm_optimize(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.18 2005-08-12 01:12:11 marat Exp $
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "qmmm.fh"
      integer rtdb
      double precision cpu, wall
c     
      logical  status
c     
      logical ouselbfgs
      character*32 theory
      character*32 optimization
      character*32 pname
c     
      logical   driver
      external  driver
c
      logical   driver1
      external  driver1
c
      logical   qmmm_lbfgsb_driver
      external  qmmm_lbfgsb_driver

      logical   qmmm_tandem_optimization
      external  qmmm_tandem_optimization
      logical   qmmm_pes
      external  qmmm_pes
      logical   lbfgsb_driver
      external  lbfgsb_driver
      integer nr,resid(1),i,k
      logical ignore
      logical divcon
      logical opes

      nr = 1


c      call util_print_push() 
c      call util_print_rtdb_load(rtdb, 'qmmm')
c
      pname = "task_qmmm_optimize"
c
      if(qmmm_print_debug())
     >  write(*,*) "in " //pname, ga_nodeid()
c
      if (.not. rtdb_cget(rtdb, 'qmmm:optimization', 1, optimization))
     $     call errquit('qmmm_optimize: failed ',0,0)
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0,
     &       RTDB_ERR)
     
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +     call errquit('task_optimize: failed rtdb_cget task:theory',0,
     &       RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'driver:lbfgs', mt_log, 
     $                   1, ouselbfgs))
     $     ouselbfgs = .false.

      if (.not. rtdb_get(rtdb, 'qmmm:pes', mt_log, 
     $                   1, opes))
     $     opes = .false.

      if (.not. rtdb_get(rtdb, 'driver:divcon', mt_log, 
     $                   1, divcon))
     $     divcon = .false.


      status = .true.

      if(opes) then
        if(ouselbfgs) then
          status=qmmm_pes(rtdb,qmmm_lbfgsb_driver)
        else
          status=qmmm_pes(rtdb,driver)
        end if
        go to 110
      end if

      if (optimization.eq.'all') then
         if(qmmm_print_debug())
     >      write(*,*) "going into md_em"
         call md_em()
      else if (optimization.eq.'tandem') then
         status =  qmmm_tandem_optimization(rtdb)
      else if (optimization.eq.'mm') then
         if(qmmm_print_debug())
     >      write(*,*) "going into md_em"
         call md_em()
      else if (optimization.eq.'qm') then
         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int, 1, 0))
     $        call errquit('task_optimize: rtdb?',0, RTDB_ERR)
         if(qmmm_print_debug())
     >      write(*,*) "going into driver"
          if(ouselbfgs) then
            status =  qmmm_lbfgsb_driver(rtdb)
          else if(divcon) then
            status = driver1(rtdb)
          else
            status = driver(rtdb)
          end if
c         call mm_print_info()
c
c        forcing driver to never fail
c        ---------------------------
         status = .true.         
      endif
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
c     $Id: task_qmmm_optimize.F,v 1.18 2005-08-12 01:12:11 marat Exp $
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
      call qmmm_geom_data_release()
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
c     $Id: task_qmmm_optimize.F,v 1.18 2005-08-12 01:12:11 marat Exp $
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
      integer ncycles

      pname = "qmmm_tandem_optimization"

      qmmm_tandem_optimization = .true.

      if (.not. rtdb_get(rtdb,'qmmm:ncycles',mt_int,1,ncycles))
     +  ncycles = 1

      write(*,*) "tandem optimization",ncycles

      do i=1,ncycles
       call mm_fix_solute(mm_classical) 
       ignore = rtdb_delete(rtdb,'geometry')
       call qmmm_geom_create_active(rtdb)
       call qmmm_set_interface("qm")
       status = driver(rtdb)
       call mm_free_solute(mm_classical)

       call mm_fix_solute(mm_quant)
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
       write(*,*) "qmmm_energy",i,energy

      end do
      return

      end

      function qmmm_pes(irtdb,optimizer)
c     $Id: task_qmmm_optimize.F,v 1.18 2005-08-12 01:12:11 marat Exp $
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
