      subroutine selci(rtdb)
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cselcifiles.fh"
      integer rtdb              ! [input]
c
c     Driver for the selected CI code
c
c     1) The configuration generate can be run at either input
c        time or now from default (SCF/MCSCF with FZC) generated input.
c
c     2) cimode ... run a one-shot CI+DIAG on the current
c        ciconf file.  This is the default mode with one root.
c
c     3) select ... run a sequence of selection calculations
c        given the list of tolerances (truncate list as each
c        one is done to facilitate restart)
c
c     4) if an moints file does not exist or if the moints
c        directive is given with no filename then do the
c        transformation
c
      logical status, oprint, oexist, oforce
      integer ma_type
      character*26 date
      integer l_mem, k_mem, lword, nelem, i
      character*20 mode
      character*255 ctmp
      integer max_thresh, nthresh
      parameter (max_thresh = 50)
      double precision stols(max_thresh)
      integer nroot, iwpt
      double precision diag_tol
c
      call ga_sync()
c
      status = rtdb_parallel(.true.)
      if (.not. context_push('selci'))
     $     call errquit('selci: context_push failed', 0)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'selci')
      oprint = util_print('information', print_low)
c
c     What is the name of the moints file and what is its format?
c     Should we regenerate it?  If it does not exist make it.
c
      if (.not. rtdb_cget(rtdb, 'selci:moints', 1, moints)) then
         if (.not. rtdb_cget(rtdb, 'fourindex:moints file', 1, moints))
     $        call util_file_prefix('moints', moints)
      endif
      if (.not. rtdb_get(rtdb, 'selci:moints:force', mt_log, 1, oforce))
     $     oforce = .false.
      inquire(file=moints, exist=oexist)
      if (oforce .or. (.not. oexist)) call moints_full(rtdb)
      call ga_sync()
c
c     What is the name of the ciconf file?  Should we delete the
c     existing one?  (not checking this now). 
c     If it does not exisit then generate a default.
c
c     What to do about frozen core? Must be consistent with 4-index
c     which is why we must do it first and have it record info somewhere
c     (which it does not do right now).
c
      if (.not. rtdb_cget(rtdb, 'selci:ciconf', 1, ciconf))
     $     call util_file_prefix('ciconf', ciconf)
      inquire(file=ciconf, exist=oexist)
      if (.not. oexist) call errquit('selci: no default conf',0)
c
c     Allocate all available MA space
c
      lword = ma_inquire_avail(mt_dbl)
      if (.not. ma_push_get(mt_dbl, lword, 'selci', l_mem, k_mem))
     $     call errquit('selci: could not allocate all memory',lword)
c
c     What about #roots, davids conv., ??????
c     Passing down file names ?????
c     Passing down selection thresholds ????
c     Passing back error conditions ????
c     tcgmsg vs GA numbering of processes.
c     use of Ga global operations
c
c     Get remaining parameters and filenames ...
c     What should we do?  Default is ci+davids from current ciconf.
c     Other alternative is a sequence of selected CI calculations
c     (ci+davids+select) with given sequence of thresholds.
c
      if (.not. rtdb_cget(rtdb, 'selci:mode', 1, mode))
     $     mode = 'ci+davids'
      if (.not. rtdb_get(rtdb, 'selci:selection criterion',
     $     mt_int, 1, iwpt)) iwpt = 0
      if (.not. rtdb_get(rtdb, 'selci:nroot',
     $     mt_int, 1, nroot)) nroot = 1
      if (.not. rtdb_get(rtdb, 'selci:diag tol',
     $     mt_dbl, 1, diag_tol)) diag_tol = 1d-6
c
c     These files only opened by process 0
c
      if (.not. rtdb_cget(rtdb, 'selci:civecs', 1, civecs))
     $     call util_file_prefix('civecs',civecs)
      if (.not. rtdb_cget(rtdb, 'selci:wmatrx', 1, wmatrx))
     $     call util_file_prefix('wmatrx',wmatrx)
      if (.not. rtdb_cget(rtdb, 'selci:eigins', 1, eigins))
     $     call util_file_prefix('eigins',eigins)
      if (.not. rtdb_cget(rtdb, 'selci:rtname', 1, names(1)))
     $     call util_file_prefix('rtname',names(1))
      if (.not. rtdb_cget(rtdb, 'selci:roname', 1, names(2)))
     $     call util_file_prefix('roname',names(2))
      if (.not. rtdb_cget(rtdb, 'selci:hamdg', 1, hamdg))
     $     call util_file_prefix('hamdg',hamdg)
c
c     This file is opened by all processes
c
      if (.not. rtdb_cget(rtdb, 'selci:hamil', 1, ctmp))
     $     call util_file_prefix('hamil',ctmp)
      call util_pname(ctmp, hamil)
c      
      if (mode .eq. 'ci+davids') then
         call selci_ci(dbl_mb(k_mem), lword)
         call selci_davids(dbl_mb(k_mem), lword, nroot, diag_tol)
      else if (mode .eq. 'select') then
         if (.not. rtdb_get_info(rtdb, 'selci:selection thresholds', 
     $        ma_type, nthresh, date)) call errquit
     $        ('selci: no selection thresholds in database', 0)
         if (nelem .gt. max_thresh) call errquit
     $        ('selci: too many selection thresholds', max_thresh)
         if (.not. rtdb_get(rtdb, 'selci:selection thresholds', 
     $        mt_dbl, max_thresh, stols)) call errquit
     $        ('selci: no selection thresholds in database?', 0)
         if (.not. rtdb_get(rtdb, 'selci:selection criterion',
     $        mt_int, 1, iwpt)) iwpt = 0
         
         do i = 1, nthresh
            call selci_ci(dbl_mb(k_mem), lword)
            call selci_davids(dbl_mb(k_mem),lword,nroot,diag_tol)
            call selci_select(dbl_mb(k_mem),lword,nroot,iwpt,stols(i))
c     
            if (i .ne. nthresh) then
               if (.not. rtdb_put(rtdb, 'selci:selection thresholds', 
     $              mt_dbl, nthresh-i, stols(i+1))) call errquit
     $              ('selci: failed writing thresholds to database?', 0)
            else
               if (.not. rtdb_delete(rtdb,'selci:selection thresholds'))
     $              call errquit('selci: failed deleting thresholds',0)
            endif
         enddo
      else
         call errquit('selci: unknown mode', 0)
      endif
c
c     What files need deleting?
c
c     That's it!
c
      if (.not. ma_pop_stack(l_mem)) 
     $     call errquit('selci: ma corrupted', 0)
c
      call ga_sync()
c
      end
