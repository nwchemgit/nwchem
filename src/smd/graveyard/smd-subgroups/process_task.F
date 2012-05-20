c     Wrapper method to invoke calculation on an SMD monomer calculation.
c     @author Sriram Krishnamoorthy
c
c     This method is invoked from the dynamic process group client
C     -server implementation.  


c     Wrapper method invoked from the C client-server scheduler
c     @param tskid IN Id of task to execute (here monomer id)
c     @param p_grp IN Handle to GA group on which to execute
c

      block data process_data
      implicit none
#include "process_task.fh"
      data populate_first_time /.true./
      data populate_done /.false./
      end

      subroutine process_task(tskid, p_grp)
      implicit none
      integer tskid, p_grp 

#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"

      integer rtdb_grp !rtdb of p_grp (to be created)
      double precision etime1,etime2
c
      external duplicate_rtdb !duplicate an rtdb file
c
      if(ga_pgroup_get_default().ne.p_grp) 
     $  call errquit('Default group not properly set',p_grp,0)
c

      call cpu_time(etime1)
      call duplicate_rtdb(rtdb_grp, p_grp)
      call cpu_time(etime2)
c      write(*,*) "time monomer rtdb",etime2-etime1

      call smd_rtdb_set_handle(rtdb_grp)
c
      call ga_pgroup_sync(ga_pgroup_get_default())
      call smd_monomer(tskid)
      call ga_pgroup_sync(ga_pgroup_get_default())
c
c      call smd_charge_ga_print()
      if(.not.rtdb_close(rtdb_grp, 'delete'))
     + call errquit('Failed to close group rtdb',0, GEOM_ERR)
c
      end !process_task

c     Create a temporary file name
c     fname_out OUT File name
c     world_ldr_ir IN proc id of leader of this grp in world
c     grp_id IN id of this group
      subroutine create_file_name(base, fname_out, world_ldr_id, 
     *     grp_id)
      implicit none
      character*(*) base, fname_out
      integer world_ldr_id, grp_id

#include "global.fh"
      character*100 name
      character*100 pname

      pname = "create_file_name"

      if(grp_id.ge.9999) 
     &   call errquit('too big a group id', grp_id, 0)
      
      if(world_ldr_id.ge.9999) 
     &   call errquit('too big a world ldr id', world_ldr_id, 0)


      name = " "
      write(name,'(A7,A1,I4.4,A1,I4.4)') base, '.', world_ldr_id, '.',
     *    grp_id

      call util_file_name(name,.true.,.false.,fname_out)


      end
      subroutine create_file_name_perm(base, fname_out, world_ldr_id, 
     *     grp_id)
      implicit none
      character*(*) base, fname_out
      integer world_ldr_id, grp_id

#include "global.fh"
      character*100 name
      character*100 pname

      pname = "create_file_name"

      if(grp_id.ge.9999) 
     &   call errquit('too big a group id', grp_id, 0)
      
      if(world_ldr_id.ge.9999) 
     &   call errquit('too big a world ldr id', world_ldr_id, 0)


      name = " "
      write(name,'(A7,A1,I4.4,A1,I4.4)') base, '.', world_ldr_id, '.',
     *    grp_id

      call util_file_name(name,.false.,.false.,fname_out)


      end

c     Replaces util_sgmovecs to make the file names group aware
      subroutine create_movecs(rtdb)
      implicit none
      integer rtdb

#include "util.fh"      
#include "rtdb.fh"
#include "global.fh"

      character*100 fname, base,fname_in
      character*255 theory, tag
      integer p_grp
      
      base = "test.movecs"
      p_grp = ga_pgroup_get_default()

      call create_file_name(base,fname, 
     *     ga_pgroup_absolute_id(p_grp, 0),
     *     p_grp)
c
c Set output movecs
c
c      if (.not.(rtdb_cget(rtdb,'task:theory',1,theory)))
c     *     call errquit('rtdb_get task:theory failed',0,0)
      if (.not.rtdb_cget(rtdb,"smd:theory",1,
     >                   theory)) then
       write(0,*) "PLEASE SET smd:theory"
       call errquit('failed to get theory',0,0)
      end if
      if (theory.eq.'dft') then
       tag='dft:output vectors'
      else
       tag='scf:output vectors'
      endif
c
      if (.not. rtdb_cput(rtdb,tag,1,fname)) then
       write (6,*) 'error in rtdb_cput'
      endif
c
      if (theory.eq.'dft') then
       tag='dft:input vectors'
      else
       tag='scf:input vectors'
      endif
      if (.not. rtdb_cget(rtdb,tag,1,fname_in)) then
       write (6,*) 'error in rtdb_cput'
      endif
      if(ga_nodeid().eq.0) then
        write(*,*) "movecs copy",fname_in,fname
        call util_file_copy(fname_in,fname)
      end if
      if (.not. rtdb_cput(rtdb,tag,1,fname)) then
       write (6,*) 'error in rtdb_cput'
      endif
ccc
      end !create_movecs


c     @param rtdb_in IN
c     @param rtdb_out OUT
c     @param p_out IN
      subroutine duplicate_rtdb(rtdb_out, p_out)
      implicit none
      integer rtdb_out, p_out

#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util_sgroup.fh"
      
      character*100 fname_in, fname_out, name, base
      integer myproc, i, world_ldr_id
      logical status, mode

      base = "dyn_db"
      world_ldr_id=0
      myproc = ga_nodeid()

      call ga_pgroup_sync(p_out)
C
C     Node 0 of the group now copies the original RTDB to a group specific one
C
      if (myproc .eq. 0) then
         world_ldr_id = ga_pgroup_nodeid(ga_pgroup_get_world())
         call create_file_name(base,fname_out, world_ldr_id, p_out)

         call smd_rtdb_get_name(fname_in)
c
c         write (*,*) 'rtdb filename_in', fname_in
c         write (*,*) 'rtdb filename_out', fname_out
c         call util_flush(6)
         inquire(file=fname_in,exist=status)
         if(.not.status) then
c     most likely cause of error: permdir not available on this node
            write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(0,*) ' ! please change permanent_dir (2) to a   !'
            write(0,*) ' ! filesystem available to all the    !'
            write(0,*) ' ! processors (e.g. NFS, GPFS, Lustre)!'
            write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            call errquit(' util_sgrtdb: failed to open rtdb',0,0)
         endif
c
         call util_file_copy(fname_in, fname_out)
        
        inquire(file=fname_out, exist=status)
        if (.not. status) then
          write (0,*) 'copy failed'
          call errquit('duplicate_file: problem with system call', 0, 0)
        else if (sgprint) then
c          write(6,*) ' rtdb copied ', ga_pgroup_get_default()
c          call flush(6)
        endif
      endif


      world_ldr_id = ga_pgroup_absolute_id(p_out, 0);
c      call ga_pgroup_brdcst(p_out, MT_INT, world_ldr_id, 1, 0);
c      call ga_igop(world_ldr_id,1,'+')
      call create_file_name(base,fname_out, world_ldr_id, p_out)

      if (.not. rtdb_open(fname_out,'old',rtdb_out)) call errquit
     *     (' rtdb_open old failed ',0,myproc)
c
c     Create movecs for each subgroup
c
      call create_movecs(rtdb_out)
c
      end !duplicate_rtdb



c     task_list IN Task list to insert tasks into
      subroutine populate_tasks(task_list)
      integer task_list

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "process_task.fh"

      integer nr, na, i, ntsks, ga_charges
      logical changed_flag
      integer lo(1), hi(1), ld(1)
      double precision diff
      integer istart,iend,nrpoc

c      istart=1
c      iend=4
      call smd_rtdb_get_istart(istart)
      call smd_rtdb_get_iend(iend)
      call smd_rtdb_get_nproc(nproc)
      if(populate_done)
     &  call errquit('trying to populate_tasks'// 
     &    ' after done',0,0)

      call smd_system_ntot(na)
      call smd_system_nres(nr)


      ntsks = 0
      if(populate_first_time) then

         call smd_charge_get_ptr(i_q)

         if(.not.ma_alloc_get(MT_DBL,
     +        na,
     +        "tmp q",
     +        tmp_h_q,
     +        tmp_i_q))
     +        call errquit('memory allocation',
     +        na,
     +        MA_ERR)

         old_i_q = i_q
         new_i_q = tmp_i_q

c         do i = 1, nr
        do i = istart,iend
            call sched_grp_insert_task(task_list, i, nproc)
            ntsks = ntsks + 1
         enddo
         
         populate_first_time = .false.
      else
c
c         call smd_charge_ga_handle(ga_charges)
c
c         lo(1) = 1
c         hi(1) = na
c         ld(1) = 1
c         call nga_get(ga_charges,lo, hi, dbl_mb(new_i_q), ld)
c
c         changed_flag = .false.
c         do i = 1, na
c            diff = dabs(dbl_mb(new_i_q+i-1)-dbl_mb(old_i_q+i-1))
c            if(diff.gt.0.00001) then
c              write(*,*) "charge difference",
c     >                   diff
c              changed_flag = .true.
c              goto 100
c            end if
c         end do
c 100     continue
c
c         if(changed_flag) then
cc            do i = 1, nr
c           do i = istart,iend
c               call sched_grp_insert_task(task_list, i, 1)
c               ntsks = ntsks + 1
c            enddo
c         endif
c
c         !Flip new and old pointers
c         i = new_i_q
c         new_i_q = old_i_q
c         old_i_q = i
      end if

      if(ntsks.eq.0) then
         if(.not.ma_free_heap(tmp_h_q))
     $        call errquit('memory deallocation',0,tmp_h_q)
         populate_done=.true.
      end if

      end !populate_tasks

c $Id$
