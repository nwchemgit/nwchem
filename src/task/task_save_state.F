C> \ingroup task
C> @{
C>
C> \brief A generic interface to save restart information
C>
C> Each module saves any files/database entries neccessary
C> to restart the calculation at its current point by appending the
C> given suffix to any names.
C>    
C> The exact (and perhaps only) application of this routine is in
C> computation of derivatives by finite difference.  The energy/gradient
C> is computed at a reference geometry (or zero field) and then
C> the wavefunction is saved by calling this routine.  Subsequent
C> calculations at displaced geometries (or non-zero fields) call
C> task_restore_state() in order to use the wavefunction at the
C> reference geometry as a starting guess for the calculation
C> at the displaced geometry.  Thus, there is no need to save basis
C> or geometry (or field) information.  E.g., in the SCF only the
C> MO vector file is saved.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function task_save_state(rtdb,suffix)
*     
*     $Id$
*     
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb              !< [Input] The RTDB handle
      character*(*) suffix      !< [Input] The file name suffix
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
      character*(nw_max_path_len) vectors, name, done(3)
      character*8 methods(3)
      logical oexist, oprint, odone
      integer nmethods, method, i
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_save_state: theory not specified',0,
     &       RTDB_ERR)
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
c     Save the MO vector file ... rather than trying to figure out
c     which theory generated the MO vectors just save all MO vectors
c     known by theories that generate MOs (currently SCF, DFT and MCSCF)
c     or the default name.  Also, don't complain if they are not there.
c     
         if ( inp_compare(.false.,'scf',theory) .or.
     $        inp_compare(.false.,'dft',theory) .or.
     $        inp_compare(.false.,'mcscf',theory) ) then
            nmethods = 1
            methods(1) = theory
         else
            nmethods = 3
            methods(1) = 'scf'
            methods(2) = 'dft'
            methods(3) = 'mcscf'
         end if
         do method = 1, nmethods
            write(name,'(a,'':input vectors'')') 
     $           methods(method)(1:inp_strlen(methods(method)))
            done(method) = ' ' 
            if (rtdb_cget(rtdb,name,1,vectors)) then
               call util_file_name_resolve(vectors, .false.)
               done(method) = vectors
               odone = .false.
               do i = 1, method-1
                  odone = odone .or. inp_compare(.true.,done(i),vectors)
               end do
               if (.not. odone) then
                  inquire(file=vectors,exist=oexist)
                  if (oexist) then
                     name = ' '
                     write(name,'(a,''.'',a)') 
     $                    vectors(1:inp_strlen(vectors)),
     $                    suffix(1:inp_strlen(suffix))
                     call util_file_name_resolve(name, .false.)
                     if (oprint) then
                        write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                        call util_flush(6)
                     endif
                     call util_file_copy(vectors, name)
                  end if
               end if
            end if
         end do
c     
c     Also do the default name
c     
         call util_file_name('movecs',.false.,.false.,vectors)
         call util_file_name_resolve(vectors, .false.)
         odone = .false.
         do i = 1, nmethods
            odone = odone .or. inp_compare(.true.,done(i),vectors)
         end do
         if (.not. odone) then
            inquire(file=vectors,exist=oexist)
            if (oexist) then
               name = ' '
               write(name,'(a,''.'',a)') 
     $              vectors(1:inp_strlen(vectors)),
     $              suffix(1:inp_strlen(suffix))
               call util_file_name_resolve(name, .false.)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
               call util_file_copy(vectors, name)
            end if
         end if
c     
c     Selci and MCSCF have CI vectors
c     
         if (theory.eq.'mcscf' .or. theory.eq.'selci') then
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
         if (theory(1:4).eq.'ccsd') then
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
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      character*(*) suffix      ! [input]
c     
c     See comments in task_save_state
c     
      logical ignore, status
      character*32 theory
      character*(nw_max_path_len) vectors, name, done(3)
      character*8 methods(3)
      logical oexist, oprint, odone
      integer nmethods, method, i
c     
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_restore_state: theory not specified',0,
     &       RTDB_ERR)
      oprint = util_print('task_restore_state',print_high)
c
      status = .true.           ! Set to false if there was nothing to restore
c     
      if (ga_nodeid() .eq. 0) then
         ignore = rtdb_parallel(.false.)
         if (oprint) then
            write(6,1) theory(1:inp_strlen(theory)), suffix
 1          format(/' Restoring state for ',a, ' with suffix ', a)
            call util_flush(6)
         endif
c     
c     Restore the MO vector file ... rather than trying to figure out
c     which theory generated the MO vectors just save all MO vectors
c     known by theories that generate MOs (currently SCF, DFT and MCSCF)
c     or the default name.  Also, don't complain if they are not there.
c     
         if ( inp_compare(.false.,'scf',theory) .or.
     $        inp_compare(.false.,'dft',theory) .or.
     $        inp_compare(.false.,'mcscf',theory) ) then
            nmethods = 1
            methods(1) = theory
         else
            nmethods = 3
            methods(1) = 'scf'
            methods(2) = 'dft'
            methods(3) = 'mcscf'
         end if
         do method = 1, nmethods
            write(name,'(a,'':input vectors'')') 
     $           methods(method)(1:inp_strlen(methods(method)))
            done(method) = ' ' 
            if (rtdb_cget(rtdb,name,1,vectors)) then
               call util_file_name_resolve(vectors, .false.)
               done(method) = vectors
               odone = .false.
               do i = 1, method-1
                  odone = odone .or. inp_compare(.true.,done(i),vectors)
               end do
               if (.not. odone) then
                  name = ' '
                  write(name,'(a,''.'',a)') 
     $                 vectors(1:inp_strlen(vectors)),
     $                 suffix(1:inp_strlen(suffix))
                  call util_file_name_resolve(name, .false.)
                  inquire(file=name,exist=oexist)
                  if (oexist) then
                     if (oprint) then
                        write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                        call util_flush(6)
                     endif
                     call util_file_copy(name, vectors)
                  end if
               end if
            end if
         end do
c     
c     Also do the default name
c     
         call util_file_name('movecs',.false.,.false.,vectors)
         call util_file_name_resolve(vectors, .false.)
         odone = .false.
         do i = 1, nmethods
            odone = odone .or. inp_compare(.true.,done(i),vectors)
         end do
         if (.not. odone) then
            name = ' '
            write(name,'(a,''.'',a)') 
     $           vectors(1:inp_strlen(vectors)),
     $           suffix(1:inp_strlen(suffix))
            call util_file_name_resolve(name, .false.)
            inquire(file=name,exist=oexist)
            if (oexist) then
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
               call util_file_copy(name, vectors)
            end if
         end if
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
            inquire(file=name,exist=oexist)
            if (oexist) then
               call util_file_copy(name, vectors)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
               call util_file_name('ciconf',.false.,.false.,vectors)
               name = ' '
               write(name,'(a,''.'',a)') vectors(1:inp_strlen(vectors)),
     $              suffix(1:inp_strlen(suffix))
               call util_file_copy(name, vectors)
               if (oprint) then
                  write(6,'(8x,a)') vectors(1:inp_strlen(vectors))
                  call util_flush(6)
               endif
            else
               if (oprint) write(6,*) name(1:inp_strlen(name)),
     $              ' does not exist'
               status = .false.
               goto 555
            endif
         endif
c
 555     if (oprint) then
            write(6,*)
            call util_flush(6)
         endif
         ignore = rtdb_parallel(.true.)
      endif
c     
      call ga_sync()
      call ga_brdcst(88, status, MA_sizeof(MT_INT,1,MT_BYTE), 0)
c     
      task_restore_state = status
c     
      end
      logical function task_delete_state(rtdb,suffix)
      implicit none
#include "errquit.fh"
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
     $     call errquit('task_delete_state: theory not specified',0,
     &       RTDB_ERR)
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
     $           ('task_delete_state: rtdb corrupt',0, RTDB_ERR)
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
C> @}
