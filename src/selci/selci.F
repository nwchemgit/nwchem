C> \file selci.F
C> The main driver for the selected CI module
C>
C> \defgroup selci Selected Configuration Interaction
C>
C> \brief The Selected Configuration Interaction Module
C>
C> The selected CI code dates back to Robert Harrison then based at
C> Argonne [1]. 
C>
C> [1] R.J. Harrison, "Approximating full configuration interaction
C>     with selected configuration interaction and perturbation theory",
C>     J. Chem. Phys. <b>94</b> (1991), 5021-5031,
C>     <a href="https://doi.org/10.1063/1.460537">
C>     10.1063/1.460537</a>.

C> \ingroup selci
C> @{
C>
C> \brief The main Selected CI driver
C>
C> This routine drives the selected CI calculations.
C>
C> \return .True. if the calculation was successfull, and
C> .False. otherwise.
C>
      logical function selci(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "cselcifiles.fh"
      integer rtdb  !< [Input] The RTDB handle
c
c     Driver for the selected CI code
c
c     0) Run SCF program unless MCSCF keyword is present in which
c        case run that.
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
      logical status, oprint, oexist, oforce, oupdate
      integer ma_type
      character*26 date
      integer l_mem, k_mem, lword, i
      character*20 mode
      character*255 ctmp
      integer max_thresh, nthresh, max_roots
      parameter (max_thresh = 50, max_roots = 50)
      double precision stols(max_thresh)
      integer nroot, iwpt
      double precision diag_tol, roots(max_roots), rootspt(max_roots)
      double precision findtol
      character*10 reference
      logical scf, mcscf, moints_full
      external scf, mcscf, moints_full
      logical cckohn_ci
      if (.not. rtdb_get(rtdb,'selci:cckohn_ci',mt_log,1,cckohn_ci))
     $     cckohn_ci=.false.
c
      call ga_sync()
c
      status = rtdb_parallel(.true.)
c
c     First generate the orbitals using SCF/MCSCF
c
      if (.not. rtdb_cget(rtdb, 'selci:reference', 1, reference))
     $     reference = 'scf'
      if (reference .eq. 'scf') then
         selci = scf(rtdb)
      else if (reference .eq. 'kohn') then
         selci = .true. 
      else
         selci = mcscf(rtdb)
      endif
c
      if (.not. selci) return
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'selci')
      oprint = util_print('information', print_low)
c
c     What is the name of the moints file and what is its format?
c     Should we regenerate it?  If it does not exist make it.
c
      if (.not. rtdb_cget(rtdb, 'selci:moints', 1, moints)) then
         if (.not. rtdb_cget(rtdb, 'fourindex:moints file', 1, moints))
     $        call util_file_name('moints', .true.,.false.,moints)
      endif
      if (.not. rtdb_get(rtdb, 'selci:moints:force', mt_log, 1, oforce))
     $     oforce = .false.
      if (.not. rtdb_get(rtdb, 'fourindex:tol', mt_dbl, 1, findtol))
     $     findtol=1d-12
      findtol=max(findtol,1d-12)
      if (.not. rtdb_put(rtdb, 'fourindex:tol', mt_dbl, 1, findtol))
     $     call errquit('selci: rtdb_put failed',0,RTDB_ERR)
      if (ga_nodeid() .eq. 0) inquire(file=moints, exist=oexist)
      call ga_brdcst(303, oexist, mitob(1), 0)
      if (oforce .or. (.not. oexist)) then
         if (.not. moints_full(rtdb)) then
            selci = .false.
            return
         endif
      endif
      call ga_sync()
c
c     Allocate all available MA space
c
      lword = ma_inquire_avail(mt_dbl)
      if (.not. ma_push_get(mt_dbl, lword, 'selci', l_mem, k_mem))
     $     call errquit('selci: could not allocate all memory',lword,
     &       MA_ERR)
c
c     tcgmsg vs GA numbering of processes ?
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
      if (.not. rtdb_get(rtdb, 'selci:update', mt_log, 1, oupdate))
     $     oupdate = .true.
c
c     These files only opened by process 0
c
      if (.not. rtdb_cget(rtdb, 'selci:ciconf', 1, ciconf))
     $     call util_file_name('ciconf', .false.,.false., ciconf)
      if (.not. rtdb_cget(rtdb, 'selci:civecs', 1, civecs))
     $     call util_file_name('civecs', .false.,.false., civecs)
      if (.not. rtdb_cget(rtdb, 'selci:wmatrx', 1, wmatrx))
     $     call util_file_name('wmatrx', .true.,.false., wmatrx)
      if (.not. rtdb_cget(rtdb, 'selci:eigins', 1, eigins))
     $     call util_file_name('eigins', .true.,.false., eigins)
      if (.not. rtdb_cget(rtdb, 'selci:rtname', 1, names(1)))
     $     call util_file_name('rtname', .false.,.false., names(1))
      if (.not. rtdb_cget(rtdb, 'selci:roname', 1, names(2)))
     $     call util_file_name('roname', .false.,.false., names(2))
      if (.not. rtdb_cget(rtdb, 'selci:hamdg', 1, hamdg))
     $     call util_file_name('hamdg', .false.,.false., hamdg)
c
c     This file is opened by all processes
c
      if (.not. rtdb_cget(rtdb, 'selci:hamil', 1, ctmp)) then
         call util_file_name('hamil', .false.,.true., hamil)
      else
         call util_pname(ctmp, hamil)
      endif
C      write(*,*)"hamiltonian matrices changed to permanant"
c
c     Should we delete the existing ciconf file?  (not checking this now). 
c     If it does not exisit then generate a default.
c
      if (ga_nodeid() .eq. 0) inquire(file=ciconf, exist=oexist)
      call ga_brdcst(304, oexist, mitob(1), 0)
      if (.not. oexist) call selci_conf(rtdb, dbl_mb(k_mem), lword)
c      
      if (mode .eq. 'ci+davids') then
         call selci_ci(dbl_mb(k_mem), lword)
         call selci_davids(dbl_mb(k_mem), lword, nroot, diag_tol, roots)
         if (.not. rtdb_put(rtdb, 'selci:roots', mt_dbl, nroot, roots))
     $        call errquit('selci: failed to store CI roots',nroot,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'selci:energy',mt_dbl,1,roots(nroot)))
     $        call errquit('selci: failed to store CI energy',nroot,
     &       RTDB_ERR)
      else if (mode .eq. 'select') then
         if (.not. rtdb_get_info(rtdb, 'selci:selection thresholds', 
     $        ma_type, nthresh, date)) then
            nthresh = 1
            stols(1) = 1.0d99   ! i.e., don't select any new configuration
         else
            if (nthresh .gt. max_thresh) call errquit
     $           ('selci: too many selection thresholds', max_thresh,
     &       INPUT_ERR)
            if (.not. rtdb_get(rtdb, 'selci:selection thresholds', 
     $           mt_dbl, max_thresh, stols)) call errquit
     $           ('selci: no selection thresholds in database?', 0,
     &       RTDB_ERR)
         endif
         if (.not. rtdb_get(rtdb, 'selci:selection criterion',
     $        mt_int, 1, iwpt)) iwpt = 0
c         
         do i = 1, nthresh
            call selci_ci(dbl_mb(k_mem), lword)
            call selci_davids(dbl_mb(k_mem),lword,nroot,diag_tol,roots)
            call selci_select(dbl_mb(k_mem),lword,nroot,iwpt,stols(i),
     $           oupdate,rootspt)
c     
            if (i .ne. nthresh) then
               if (.not. rtdb_put(rtdb, 'selci:selection thresholds', 
     $              mt_dbl, nthresh-i, stols(i+1))) call errquit
     $              ('selci: failed writing thresholds to database?', 0,
     &       RTDB_ERR)
            else
               status = rtdb_delete(rtdb,'selci:selection thresholds')
            endif
         enddo
         if (.not. rtdb_put(rtdb, 'selci:roots', mt_dbl, nroot, roots))
     $        call errquit('selci: failed to store CI roots',nroot,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'selci:roots+pt',
     $        mt_dbl, nroot, rootspt))
     $        call errquit('selci: failed to store CI+PT roots',nroot,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'selci:energy',mt_dbl,1,
     $        rootspt(nroot)))
     $        call errquit('selci: failed to store CI+PT energy',nroot,
     &       RTDB_ERR)
      else
         call errquit('selci: unknown mode', 0, INPUT_ERR)
      endif
c
c     What files need deleting?
c
      if(.not.cckohn_ci)
     &     call util_file_unlink(hamil)
      if (ga_nodeid() .eq. 0) then
         if(.not.cckohn_ci)call util_file_unlink(hamdg)
         call util_file_unlink(wmatrx)
      endif
c
c     That's it!
c
      if (.not. ma_pop_stack(l_mem)) 
     $     call errquit('selci: ma corrupted', 0, MA_ERR)
      call util_print_pop
c
      call ga_sync()
c
c     If anything went wrong a hard fail would have resulted
c
      selci = .true.
c
      end
C>
C> @}
