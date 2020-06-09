c
c $Id$
c
C> \ingroup task
C> @{
      function task_qmmm_pes(irtdb)
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
      logical task_qmmm_pes
c     
      character*30 operation
      character*30 theory
      logical  status
c     
      character*32 pname
c     
      logical   task_qmmm_energy_sp
      external  task_qmmm_energy_sp

      logical   qmmm_dynamics
      external  qmmm_dynamics

      logical   task_qmmm_gradient
      external  task_qmmm_gradient

      logical   task_qmmm_property
      external  task_qmmm_property

      logical   task_qmmm_esp
      external  task_qmmm_esp

      logical   qmmm_mixed_optimization
      external  qmmm_mixed_optimization

      logical ignore
      character*255 filename, prefix
      character*255 rst_in, movecs_in, espfile
      character*255 rst_perm,rst_tmp,movecs_perm
      character*255 rst_old
      character*255 geom_in,geom_tmp
      character*255 logfile
      double precision rb(2)
      integer ii
      integer np
      integer ip
      integer in
      double precision dr
      double precision r
      integer i1,i2,i3
      integer ilog
      logical ofile,ofile1
      logical backward
      logical loadrst
      logical loadmovecs
      integer ncycles
      integer istart
      integer iend
      logical pes_spring
      character*84 inmovecs_tag
      character*84 outmovecs_tag
      character*84 espmovecs_tag
      double precision energy
      double precision bq_el_energy
      double precision bq_nuc_energy
      double precision mm_energy
      double precision qm_energy
      double precision temp
      logical lvel
      logical operm
      integer tce_ref
 
      pname = "task_qmmm_pes"
c
c     which operation or theory?
c     -------------------------
      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = ' '
      if(.not.rtdb_cget(irtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)
      
      if(theory.eq.'esp') operation=theory
c
c     default file prefix from mm module
c     ----------------------------------
      call mm_system_prefix(prefix)
      prefix = prefix(1:inp_strlen(prefix))//"_pes"
c
      call mm_restart_filename(rst_old)
c
c     get number of cycles 
c     --------------------
      if (.not. rtdb_get(irtdb,'qmmm:pes_ncycles',mt_int,
     +                    1,ncycles))
     + ncycles = 1
c
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
      if (.not. rtdb_cget(irtdb,'qmmm:pes_geom_input',
     +                    1,geom_in))
     +  geom_in = " "
c
      if (.not. rtdb_cget(irtdb,'qmmm:pes_rst_input',
     +                    1,rst_in))
     +  rst_in = " "
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_rst_perm',mt_log,
     +                    1,operm))
     +  operm = .true.
c
      if (.not. rtdb_get(irtdb,'qmmm:pes_vel',mt_log,
     +                    1,lvel))
     +  lvel = .false.
c
      if (.not. rtdb_cget(irtdb,'qmmm:pes_vectors_input',
     +                    1,movecs_in))
     +  movecs_in = " "
c
c     was logfile unit number provided
c     (created for abc calculations)
c     -------------------------------- 
      if (.not. rtdb_get(irtdb,'qmmm:pes_logfile',mt_int,1,ilog))
     +  ilog = -1
c      
      if(rst_in .eq. " ") then
       loadrst = .false.
       rst_in = prefix
      else
       loadrst = .true.
      end if

      if(movecs_in .eq. " ") then
       loadmovecs = .false.
       movecs_in = prefix 
      else
       loadmovecs = .true.
      end if
c
      if(theory.eq."dft") then
        inmovecs_tag = 'dft:input vectors'
        outmovecs_tag = 'dft:output vectors'
      else if(theory.eq."tce") then
        if (.not.rtdb_get(irtdb,'tce:reference',mt_int,1,tce_ref))
     +    call errquit('failed to get  tce:reference',0,RTDB_ERR)

        if(tce_ref.eq.0) then
        inmovecs_tag = 'dft:input vectors'
        outmovecs_tag = 'dft:output vectors'
        else
        inmovecs_tag = "scf:input vectors"
        outmovecs_tag = "scf:output vectors"
        end if
      else if(theory.eq."mp2") then
        inmovecs_tag = "scf:input vectors"
        outmovecs_tag = "scf:output vectors"
       else if(theory.eq."scf") then
        inmovecs_tag = "scf:input vectors"
        outmovecs_tag = "scf:output vectors"
      else if(theory.eq."esp") then
        inmovecs_tag = "esp:input vectors"
        outmovecs_tag = "esp:input vectors"
      else
        call errquit(pname//'theory not supported'//theory,0,0)
      end if
c
      espmovecs_tag = "esp:input vectors"
c
c     this will force property module
c     not to recalculate movecs
c     -------------------------------
      if(operation.eq."property") then
        inmovecs_tag = "prop:vectors"
        outmovecs_tag = "prop:vectors"
      end if

c 
c     set increments and initial value for spring
c     -------------------------------------------
      if(pes_spring) then
        dr = (rb(2)-rb(1))/(np-1)
        r  = rb(1)
      end if

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


      ignore = rtdb_delete(irtdb,outmovecs_tag)
      ignore = rtdb_delete(irtdb,inmovecs_tag)
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
        else
          if(ga_nodeid().eq.0) then
            write(*,23) ip
          end if
        end if

c
c       initiating restart files
c       -----------------------
        rst_tmp  = " "
        rst_perm = " "
        write(rst_tmp,21) rst_in(1:inp_strlen(rst_in)), ip
 21     format(a,i3.3,'.rst')
        rst_perm = rst_tmp
        ii = index(rst_perm,".rst")-1
        espfile = rst_perm(1:ii)//".esp"
        if (.not.rtdb_cput(irtdb,'qmmm:esp:save:file',1,espfile)) 
     >   call errquit(pname,0,RTDB_ERR)
        if(operm)  
     >    call util_file_name_resolve(rst_perm, .false.)
        if(loadrst) then
          if(ga_nodeid().eq.0) then 
            inquire(file=rst_perm,exist=ofile)
            if(ofile) then 
             write(*,*) "pes: found restart file ",
     >                 rst_perm(1:inp_strlen(rst_perm))
             if(operm)  
     >          call util_file_copy(rst_perm(1:inp_strlen(rst_perm)),
     >                          rst_tmp(1:inp_strlen(rst_tmp)))
            else
             write(*,*) "pes: cannot find restart file "
            end if
          end if
          call ga_brdcst(msg_qmmm_misc, ofile,
     $         MA_sizeof(MT_INT,1,MT_BYTE), 0)
          call ga_sync()
          if(ofile) then
            call mm_reload_rst(irtdb,rst_tmp)
          else
            call errquit(pname//'cannot find rst file'//rst_perm,
     >                   0,0) 
          end if
        end if
        call mm_get_temp(temp)
c
c       initiating movecs files
c       ----------------------
        movecs_perm = ' '
        write(movecs_perm,22) movecs_in(1:inp_strlen(movecs_in)), ip
 22     format(a,i3.3,'.movecs')
        if(loadmovecs) then
          ignore = rtdb_delete(irtdb,inmovecs_tag)
          if (.not. rtdb_cput(irtdb,inmovecs_tag,
     +                      1,movecs_perm))
     +    call errquit(pname//'Failed to set movecs',0,0)
        end if
        ignore = rtdb_delete(irtdb,outmovecs_tag)
        if (.not. rtdb_cput(irtdb,espmovecs_tag,
     +                    1,movecs_perm))
     +  call errquit(pname//'Failed to set movecs',0,0)
        if (.not. rtdb_cput(irtdb,outmovecs_tag,
     +                    1,movecs_perm))
     +  call errquit(pname//'Failed to set movecs',0,0)
c
c       loading geometries if any
c       -------------------------
        if(geom_in.ne." ") then
          write(geom_tmp,34) geom_in(1:inp_strlen(geom_in)), ip
 34       format(a,i3.3,'.xyzi')
c          call util_file_name_resolve(geom_tmp, .false.)
          write(*,*) geom_tmp
          call mm_set_solute_coord_file(geom_tmp)
        end if

        if(operation.eq."energy") then
          status = task_qmmm_energy_sp(irtdb) 
        else if (operation.eq."gradient") then
          status = task_qmmm_gradient(irtdb)
        else if (operation.eq."dynamics") then
          status = qmmm_dynamics(irtdb)
          lvel =.true.
        else if (operation.eq."optimize") then
          status = qmmm_mixed_optimization(irtdb)
        else if (operation.eq."property") then
          status = task_qmmm_property(irtdb)
        else if (operation.eq."esp") then
          status = task_qmmm_esp(irtdb)
        end if 
 
        if(ga_nodeid().eq.0) then 
           write(*,*) "pes: saving restart file",
     >                 rst_perm
        end if
        call mm_write_restart_named(rst_perm,lvel)
        if(operm) then
        if(ga_nodeid().eq.0) then 
          if(rst_tmp.ne.rst_old) then
            write(*,*) "pes: removing file",rst_tmp
            call util_file_unlink(rst_tmp)
          end if
        end if
        end if

        if(operation.eq."energy".or.
     >     operation.eq."gradient".or.
     >     operation.eq."optimize") then
           if (.not. rtdb_get(irtdb,'qmmm:qm_energy',
     $                        mt_dbl,1,qm_energy))
     $          call errquit(pname//" qm energy", 0, RTDB_ERR)
           if (.not. rtdb_get(irtdb,'qmmm:mm_energy',
     $                        mt_dbl,1,mm_energy))
     $          call errquit(pname//" mm energy", 0, RTDB_ERR)
           if (.not. rtdb_get(irtdb,'qmmm:energy',
     $                        mt_dbl,1,energy))
     $         call errquit('qmmm: failed get energy', 0, RTDB_ERR)
           if (.not. rtdb_get(irtdb,'dft:bq_energy',
     $                        mt_dbl,1,bq_el_energy))
     $         bq_el_energy = 0.0d0
           if (.not. rtdb_get(irtdb,'qmmm:bq_nuc_energy',mt_dbl,1,
     >                   bq_nuc_energy))
     $     bq_nuc_energy = 0.0d0

          if(ga_nodeid().eq.0) then
            if(ilog.ne.-1)
     >       write(ilog,*) temp,energy,mm_energy,qm_energy,
     >                     bq_el_energy,bq_nuc_energy
            if(pes_spring) then
              write(*,24) ip ,energy, r
            else
              write(*,25) ip ,energy
            end if
          end if
        end if

        if(.not.loadmovecs) then
          ignore = rtdb_delete(irtdb,inmovecs_tag)
          if (.not. rtdb_cput(irtdb,inmovecs_tag,
     +                      1,movecs_perm))
     +    call errquit(pname//'Failed to set movecs',0,0)
        end if
 
        call ga_sync()

      end do
      backward = .not.backward
      loadrst = .true.
      loadmovecs = .true.
      end do
c
c     restore old restart filename
c     in case of more tasks
c     -----------------------------
      call mm_set_restart_file(irtdb,rst_old)
      ignore = rtdb_delete(irtdb,outmovecs_tag)
      ignore = rtdb_delete(irtdb,inmovecs_tag)
      

      task_qmmm_pes = status

 13      format(a,'/',a,'_neb',i3.3,'.rst')
 14      format(a,'_neb',i3.3,'.rst')
 15      format("@pes: point ", I2,3X,"spring length = ",F12.6)
 16      format("@pes: energy =  ", F18.12,3X,"length = ",F12.6)
 17      format(a,'/',a,'_pes',i3.3,'.movecs')
 23      format("@pes: processing point ", I2)
 24      format("@pes: point = ",I2,
     >          " energy =  ", F18.12,3X,
     >          " spring length = ",F12.6)
 25      format("@pes: point = ",I2,
     >          " energy =  ", F18.12,3X)

      return
      end
C> @}


