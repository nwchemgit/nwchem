      logical function task_save_state(rtdb,suffix)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb              ! [input]
      character*(*) suffix      ! [input]
c     
c     Each module saves any files/databse entries neccessary
c     to restart the calculation at its current point by appending the
c     given suffix to any names.
c     
c     The exact (and perhaps only) application of this routine is in
c     computation of derivatives by finite difference.  The energy/gradient
c     is computed at a reference geometry (or zero field) and then
c     the wavefunction is saved by calling this routine.  Subsequent
c     calculations at displaced geometries (or non-zero fields) call
c     task_restore_state() in order to use the wavefunction at the
c     reference geometry as a starting guess for the calculation
c     at the displaced geometry.  Thus, there is no need to save basis
c     or geometry (or field) information.  E.g., in the SCF only the
c     MO vector file is saved.
c     
c     Input argument ... the suffix
c     
c     RTDB arguments ... the theory name
c     
c     Output ... function value T/F on success/failure
c     
      logical status
      character*20 theory
      character*(nw_max_path_len) vectors, name
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_save_state: theory not specified',0)
c     
      if ( (theory .eq. 'scf') .or.
     $     (theory .eq. 'dft') .or.
     $     (theory .eq. 'mcscf') ) then
c
         if (ga_nodeid() .eq. 0) then
            status = rtdb_parallel(.false.)
c     
c     Get the name of the input MO vector file
c     
            name = ' '
            write(name,'(a,'':input vectors'')') 
     $           theory(1:inp_strlen(theory))
            if (.not. rtdb_cget(rtdb,name,1,vectors)) call errquit
     $           ('task_save_state: no DFT/SCF/MCSCF vectors',0)
c     
c     Save the filename in the database as "theory:vectors <suffix>"
c     
            name = ' '
            write(name,'(a,'':vectors '',a)') 
     $           theory(1:inp_strlen(theory)),
     $           suffix(1:inp_strlen(suffix))
            if (.not. rtdb_cput(rtdb,name,1,vectors)) call errquit
     $           ('task_save_state: rtdb corrupt',0)
c     
c     Copy the file to <filename>.<suffix>
c     
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(vectors, name)
c
            if (util_print('task_save_state',print_default)) then
               write(6,1) theory(1:inp_strlen(theory)),
     $              vectors(1:inp_strlen(vectors)),
     $              name(1:inp_strlen(name))
 1             format(/' Saving state for ',a/
     $              ' Original vectors = ',a/
     $              ' Backup vectors   = ',a//)
               call util_flush(6)
            endif
            status = rtdb_parallel(.true.)
         endif
      else if (ga_nodeid() .eq. 0) then
         write(6,*) ' task_save_state: ', 
     $        theory(1:inp_strlen(theory)),
     $        suffix(1:inp_strlen(suffix))
         call errquit('task_save_state: unknown theory',0)
      endif
c
      call ga_sync()
c     
      task_save_state = .true.
c     
      end
      logical function task_restore_state(rtdb,suffix)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb              ! [input]
      character*(*) suffix      ! [input]
c     
c     See comments in task_save_state
c     
      logical status
      character*20 theory
      character*(nw_max_path_len) vectors, name
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_restore_state: theory not specified',0)
c     
      if ( (theory .eq. 'scf') .or.
     $     (theory .eq. 'dft') .or.
     $     (theory .eq. 'mcscf') ) then
         if (ga_nodeid() .eq. 0) then
            status = rtdb_parallel(.false.)
c     
c     Get original name from the database entry "theory:vectors <suffix>"
c     
            name = ' '
            write(name,'(a,'':vectors '',a)') 
     $           theory(1:inp_strlen(theory)),
     $           suffix(1:inp_strlen(suffix))
            if (.not. rtdb_cget(rtdb,name,1,vectors)) call errquit
     $           ('task_restore_state: no saved DFT/SCF/MCSCF MOs',0)
c     
c     Now make the name of the backup file and copy back to the original
c     
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(name, vectors)
c
            if (util_print('task_state',print_high)) then
               write(6,1) theory(1:inp_strlen(theory)),
     $              vectors(1:inp_strlen(vectors)),
     $              name(1:inp_strlen(name))
 1             format(/' Restoring state for ',a/
     $              ' Original vectors = ',a/
     $              ' Backup vectors   = ',a//)
               call util_flush(6)
            endif
c
c     Also set the entry in the database for the input vectors
c     ... should we also do the output vectors?  Probably should
c     have saved this.
c
            name = ' '
            write(name,'(a,'':input vectors'')') 
     $           theory(1:inp_strlen(theory))
            if (.not. rtdb_cput(rtdb,name,1,vectors)) call errquit
     $           ('task_restore_state: rtdb corrupt',0)
            status = rtdb_parallel(.true.)
         endif
      else if (ga_nodeid() .eq. 0) then
         write(6,*) ' task_restore_state: ', 
     $        theory(1:inp_strlen(theory)),
     $        suffix(1:inp_strlen(suffix))
         call errquit('task_restore_state: unknown theory',0)
      endif
c     
      call ga_sync()
c     
      task_restore_state = .true.
c     
      end
      logical function task_delete_state(rtdb,suffix)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb              ! [input]
      character*(*) suffix      ! [input]
c     
c     See comments in task_save_state
c     
      logical status
      character*20 theory
      character*(nw_max_path_len) vectors, name
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_delete_state: theory not specified',0)
c     
      if ( (theory .eq. 'scf') .or.
     $     (theory .eq. 'dft') .or.
     $     (theory .eq. 'mcscf') ) then
c
         if (ga_nodeid() .eq. 0) then
            status = rtdb_parallel(.false.)
c     
c     Get original name from the database entry "theory:vectors <suffix>"
c     
            name = ' '
            write(name,'(a,'':vectors '',a)') 
     $           theory(1:inp_strlen(theory)),
     $           suffix(1:inp_strlen(suffix))
            if (.not. rtdb_cget(rtdb,name,1,vectors)) call errquit
     $           ('task_delete_state: no saved DFT/SCF/MCSCF vectors',0)
c     
c     Now delete the backup file
c     
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_unlink(name)
c     
            if (util_print('task_state',print_default)) then
               write(6,1) theory(1:inp_strlen(theory)),
     $              vectors(1:inp_strlen(vectors)),
     $              name(1:inp_strlen(name))
 1             format(/' Deleting saved state for ',a/
     $              ' Original vectors = ',a/
     $              ' Backup vectors   = ',a//)
               call util_flush(6)
            endif
            status = rtdb_parallel(.true.)
         endif
      else if (ga_nodeid() .eq. 0) then
         write(6,*) ' task_delete_state: ', 
     $        theory(1:inp_strlen(theory)),
     $        suffix(1:inp_strlen(suffix))
         call errquit('task_delete_state: unknown theory',0)
      endif
c     
      call ga_sync()
c     
      task_delete_state = .true.
c     
      end
