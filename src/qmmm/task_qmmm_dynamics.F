C> \ingroup task
C> @{
       function task_qmmm_dynamics(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb !< [Input] The RTDB handle
      logical  task_qmmm_dynamics
c
      logical status
      double precision cpu,wall
c
      character*84 tag
      character*30 pname
c
      logical opes
c
      logical qmmm_dynamics
      external qmmm_dynamics
      logical  qmmm_pes_dynamics
      external qmmm_pes_dynamics
c
      pname = "task_qmmm_dynamics"
c
      if (.not. rtdb_get(rtdb, 'qmmm:pes', mt_log, 
     $                   1, opes))
     $     opes = .false.
c
      task_qmmm_dynamics = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     set output parameters as if failed
c     ----------------------------------
      status = rtdb_delete(rtdb, 'task:dynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dynamics: failed to invalidate status',0,
     &       RTDB_ERR)
c
      if(opes) then
          status = qmmm_pes_dynamics(rtdb)
      else 
          status = qmmm_dynamics(rtdb)
      end if
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_dynamics: failed to set status',0, RTDB_ERR)
c
      if (status) then
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dynamics: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dynamics: failed storing walltime',0,
     &       RTDB_ERR)
      endif
c         
      task_qmmm_dynamics = status
      end
C> @}

       function qmmm_dynamics(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "qmmm.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
      logical  qmmm_dynamics
c
      character*30 pname
c
      character*32 interface
      character*84 tag
      character*30 region(3)
      integer nregion
      logical oesp
c
      pname = "qmmm_dynamics"
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
      if (.not. rtdb_get(rtdb, 'qmmm:esp_density', mt_log, 
     $                   1, oesp))
     $     oesp = .false.
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_free_exclusive(region(1))
c
c      call qmmm_bq_data_update_active()
c
      interface = qmmm_get_interface()
      call qmmm_set_interface("mm")
c
c     Actually do the deed
c
      if(oesp) call qmmm_esp_reset(rtdb)
      call md_md()
      call mm_write_restart_default_dynamics()
c
c     restore interface
      call qmmm_set_interface(interface)
c
      qmmm_dynamics = .true.
c
      end

      function qmmm_pes_dynamics(irtdb)
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
      logical qmmm_pes_dynamics
      double precision cpu, wall
c     
      logical  status
c     
      character*32 pname
c     
      logical   qmmm_dynamics
      external  qmmm_dynamics
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
c        filename0 = ' '
c        call mm_restart_filename(filename0)
        filename = ' '
        write(filename0,14) 
     $     prefix(1:inp_strlen(prefix)),ip
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
        status = qmmm_dynamics(irtdb)
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
      end do
      qmmm_pes_dynamics = status

 13      format(a,'/',a,'_pes',i3.3,'.rst')
 14      format(a,'_pes',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)
 17      format(a,'/',a,'_pes',i3.3,'.movecs')

      return
      end
