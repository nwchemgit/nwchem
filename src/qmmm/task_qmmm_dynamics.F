       function task_qmmm_dynamics(rtdb)
*
* $Id: task_qmmm_dynamics.F,v 1.6 2005-09-22 20:54:07 marat Exp $
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
      logical  task_qmmm_dynamics
c
      logical status
      double precision cpu,wall
c
      character*84 tag
      character*30 region(3)
      character*30 pname
      integer nregion
c
      pname = "task_qmmm_dynamics"

c
c
c      call util_print_push() 
c      call util_print_rtdb_load(rtdb, 'qmmm')
c
      call ecce_print_module_entry('task dynamics')
c
      task_qmmm_dynamics = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()

c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:dynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dynamics: failed to invalidate status',0,
     &       RTDB_ERR)
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
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region(1))
c
      call qmmm_bq_data_update_active()
c
      call qmmm_set_interface("mm")
c
c     Actually do the deed
c
      call md_md()
      status = .true.
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_dynamics: failed to set status',0, RTDB_ERR)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dynamics: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dynamics: failed storing walltime',0,
     &       RTDB_ERR)
c
      endif
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task dynamics', 'ok')
      else
         call ecce_print_module_exit('task dynamics', 'failed')
      endif
c         
      task_qmmm_dynamics = status
c
c      call util_print_pop()
c
      end

       function qmmm_dynamics(rtdb,region)
*
* $Id: task_qmmm_dynamics.F,v 1.6 2005-09-22 20:54:07 marat Exp $
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
      character*(*) region
      logical  qmmm_dynamics
c
      character*30 pname
c
      pname = "qmmm_dynamics"

c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set("fix","solute")
      call qmmm_cons_set("fix","solvent")
      call qmmm_cons_set("free",region)
c
      call qmmm_bq_data_update_active()
c
      call qmmm_set_interface("mm")
c
c     Actually do the deed
c
      call md_md()
c
c     restore interface
      call qmmm_set_interface("qm")
c
      qmmm_dynamics = .true.
c
      end

      function qmmm_pes_dynamics(irtdb)
c     $Id: task_qmmm_dynamics.F,v 1.6 2005-09-22 20:54:07 marat Exp $
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
      logical qmmm_pes_dynamics
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical ignore
      character*255 filename, dir,xyz,prefix
      character*255 filename0
      double precision energy
      double precision rb(2)
      integer np
      integer ip
      double precision dr
      double precision r
      integer i1,i2,i3
      logical ofile
      logical backward
      logical loadrst
      logical qmmm_dynamics
      external qmmm_dynamics

      pname = "qmmm_pes_dynamics"

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
      if (.not. rtdb_get(irtdb,'qmmm:pes_backward',mt_log,
     +                    1,backward))
     +  backward = .false.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_loadrst',mt_log,
     +                    1,loadrst))
     +  loadrst = .false.
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

      if(backward) then
        i1=np
        i2=1
        i3=-1
      else
        i1=1
        i2=np
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
          call ga_brdcst(msg_qmmm_misc, ofile, mitob(1), 0)
          call ga_sync()
          if(ofile) 
     >      call mm_reload_rst(irtdb,filename0)
        end if
c       end load restart 
        status = qmmm_dynamics(irtdb,"mm")
        if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $       call errquit('qmmm: failed get energy', 0, RTDB_ERR)
        if(ga_nodeid().eq.0) 
     >     write(*,16) energy,r
        call mm_write_restart_named(filename)
      end do

 13      format(a,'/',a,'_pes',i3.3,'.rst')
 14      format(a,'_pes',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)

      return
      end

