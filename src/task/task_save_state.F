      logical function task_save_state(rtdb,suffix)
*
* $Id: task_save_state.F,v 1.7 2000-01-05 00:25:53 windus Exp $
*
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
      character*32 theory
      character*(nw_max_path_len) vectors, name
      logical oexist, oprint
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_save_state: theory not specified',0)
      oprint = util_print('task_save_state',print_default)
c
      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         if (oprint) then
            write(6,1) theory(1:inp_strlen(theory)), suffix
 1          format(/' Saving state for ',a, ' with suffix ', a)
            call util_flush(6)
         endif
c     
c     Get the name of the input MO vector file ... look for 
c     "theory:input vectors" otherwise use the default
c     
         name = ' '
         write(name,'(a,'':input vectors'')') 
     $        theory(1:inp_strlen(theory))
         if (rtdb_cget(rtdb,name,1,vectors)) then
c     
c     Default is overriden ... save the filename in the database as 
c     "theory:vectors <suffix>"
c     
            name = ' '
            write(name,'(a,'':vectors '',a)') 
     $           theory(1:inp_strlen(theory)),
     $           suffix(1:inp_strlen(suffix))
            if (.not. rtdb_cput(rtdb,name,1,vectors)) call errquit
     $           ('task_save_state: rtdb corrupt',0)
         else
            call util_file_name('movecs',.false.,.false.,vectors)
         endif
c     
c     Copy the file to <filename>.<suffix>
c     
         name = ' '
         write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $        suffix(1:inp_strlen(suffix))
         call util_file_copy(vectors, name)
c     
         if (oprint) then
            write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
            call util_flush(6)
         endif
c
         if (theory.eq.'mcscf' .or. theory.eq.'selci') then
c     
c     Also need to save the CI vectors for MCSCF.  Don't know if
c     the selci used SCF or MCSCF.
c     
            call util_file_name('civec',.false.,.false.,vectors)
            inquire(file=vectors,exist=oexist)
            if (oexist) then
               name = ' '
               write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $              suffix(1:inp_strlen(suffix))
               call util_file_copy(vectors, name)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            endif
         endif
c
         if (theory.eq.'ccsd' .or. theory.eq.'ccsd(t)'
     $        .or. theory.eq.'ccsd+t(ccsd)') then
c     
c     Also need to save the CCSD amplitudes
c     
            call util_file_name('t2',.false.,.false.,vectors)
            inquire(file=vectors,exist=oexist)
            if (oexist) then
               name = ' '
               write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $              suffix(1:inp_strlen(suffix))
               call util_file_copy(vectors, name)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            endif
         endif
c
         if (theory .eq. 'selci') then
c     
c     Also need to save the CI vectors, and ciconf file (moints?)
c     
            call util_file_name('civecs',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(vectors, name)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
c
            call util_file_name('ciconf',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(vectors, name)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
         endif
c
         if (oprint) then
            write(6,*)
            call util_flush(6)
         endif
         status = rtdb_parallel(.true.)
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
      character*32 theory
      character*(nw_max_path_len) vectors, name
      logical oexist, oprint
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_restore_state: theory not specified',0)
      oprint = util_print('task_restore_state',print_high)
c     
      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         if (oprint) then
            write(6,1) theory(1:inp_strlen(theory)), suffix
 1          format(/' Restoring state for ',a, ' with suffix ', a)
            call util_flush(6)
         endif
c     
c     Get original MOvecs name from the database entry 
c     "theory:vectors <suffix>" or use the default if this is not present.
c     
         name = ' '
         write(name,'(a,'':vectors '',a)') 
     $        theory(1:inp_strlen(theory)),
     $        suffix(1:inp_strlen(suffix))
         if (rtdb_cget(rtdb,name,1,vectors)) then
c     
c     Set the entry in the database for the input vectors since it was
c     set when we saved state and is therefore probably not the default
c     
            name = ' '
            write(name,'(a,'':input vectors'')') 
     $           theory(1:inp_strlen(theory))
            if (.not. rtdb_cput(rtdb,name,1,vectors)) call errquit
     $           ('task_restore_state: rtdb corrupt',0)
         else
            call util_file_name('movecs',.false.,.false.,vectors)
         endif
c     
c     Now make the name of the backup file and copy back to the original
c     
         name = ' '
         write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $        suffix(1:inp_strlen(suffix))
         call util_file_copy(name, vectors)
c     
         if (oprint) then
            write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
            call util_flush(6)
         endif
c
         if (theory.eq.'mcscf' .or. theory.eq.'selci') then
c     
c     Also need to restore the CI vectors for MCSCF.  Don't know if
c     the selci used SCF or MCSCF so just copy the file if it exists.
c     
            call util_file_name('civec',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            inquire(file=name,exist=oexist)
            if (oexist) then
               call util_file_copy(name,vectors)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            endif
         endif
c
         if (theory.eq.'ccsd' .or. theory.eq.'ccsd(t)'
     $        .or. theory.eq.'ccsd+t(ccsd)') then
c     
c     Also need to restore the CCSD amplitudes
c     
            call util_file_name('t2',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            inquire(file=name,exist=oexist)
            if (oexist) then
               call util_file_copy(name,vectors)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            endif
         endif
c
         if (theory .eq. 'selci') then
c     
c     Also need to restore the CI vectors, and ciconf file (moints?)
c     
            call util_file_name('civecs',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(name, vectors)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
c
            call util_file_name('ciconf',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_copy(name, vectors)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
         endif
c
         if (oprint) then
            write(6,*)
            call util_flush(6)
         endif
         status = rtdb_parallel(.true.)
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
      character*32 theory
      character*(nw_max_path_len) vectors, name
      logical oexist, oprint
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_delete_state: theory not specified',0)
      oprint = util_print('task_delete_state',print_default)
c     
      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         if (oprint) then
            write(6,1) theory(1:inp_strlen(theory)), suffix
 1          format(/' Deleting state for ',a, ' with suffix ', a)
            call util_flush(6)
         endif
c     
c     Get original name from the database entry "theory:vectors <suffix>"
c     or use the default if this is not present
c     
         name = ' '
         write(name,'(a,'':vectors '',a)') 
     $        theory(1:inp_strlen(theory)),
     $        suffix(1:inp_strlen(suffix))
         if (rtdb_cget(rtdb,name,1,vectors)) then
c     Delete the database entry
            if (.not. rtdb_delete(rtdb, name)) call errquit
     $           ('task_delete_state: rtdb corrupt',0)
         else
            call util_file_name('movecs',.false.,.false.,vectors)
         endif
c     
c     Now delete the backup file
c     
         name = ' '
         write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $        suffix(1:inp_strlen(suffix))
         call util_file_unlink(name)
c     
         if (oprint) then
            write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
            call util_flush(6)
         endif
c     
         if (theory.eq.'mcscf' .or. theory.eq.'selci') then
c     
c     Delete the MCSCF CI vectors
c     
            call util_file_name('civec',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            inquire(file=name,exist=oexist)
            if (oexist) then
               call util_file_unlink(name)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            endif
         endif
c
         if (theory .eq. 'selci') then
c     
c     Delete the selected ci files ciconf and civecs
c     
            call util_file_name('civecs',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_unlink(name)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
c
            call util_file_name('ciconf',.false.,.false.,vectors)
            name = ' '
            write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_unlink(name)
            if (oprint) then
               write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
               call util_flush(6)
            endif
         endif
c     
         if (oprint) then
            write(6,*)
            call util_flush(6)
         endif
         status = rtdb_parallel(.true.)
      endif
c     
      call ga_sync()
c     
      task_delete_state = .true.
c     
      end
