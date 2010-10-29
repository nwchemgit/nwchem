      logical function cckohn(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
c#include "cgridfile.fh" 
      integer rtdb              ! [input]
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
      character*10 reference
      logical scf, mcscf, moints_full
      external scf, mcscf, moints_full
c
c     declarations for free functions
c
      integer lmax,nper,nint,np,lmaxnp,npd,lbig
      parameter(npd=4015,lbig=6)
      double precision xstart,xmax,alpha
      integer l_r,k_r,l_wwl,k_wwl,l_derl,k_derl,l_csc,k_csc,l_csl,k_csl
c
c     declaration for cckohn
c     ltop:    highest angular momentum for channel functions 
c     mtop:    highest angular momentum projection for channel functions
c     lmtop:   largest number of l, m pairs for channel functions, determined 
c              by the values of ltop and mtop
c     nchnl:   maximum number of channels allowed. 
c     nstate:  nstate = nchnl**2
c     nsymst:  nsymst = nchnl*(nchnl+1)/2
c     nbfmax:  maximum number of basis functions
c     nbfcmx:  maximum number of basis functions allowed for each chammel.
c     maxene:  maximum number of incident energies
c
      integer ltop, mtop, lmtop, nchnl, nstate, nsymst, nbfmax, nbfcmx,
     &     maxene
      parameter (ltop=6, mtop=4, lmtop=43, nchnl=20, nstate=nchnl**2, 
     &     nsymst=nchnl*(nchnl+1)/2, nbfmax=130, nbfcmx=130, maxene=20)
      double precision echan(nchnl),energy(maxene),kchan(nchnl)
c     
c     declarations for free-free and bound-free matrix elements
c     
      double complex 
     1     ovbf(lmtop,nbfcmx,nchnl),
     2     hpvb(lmtop,nbfcmx,nstate),
     3     hpvhp(lmtop,lmtop,nsymst),
     4     hpvhm(lmtop,lmtop,nstate)
c
c     declarations for overlap and bound-bound integrals
c
      double precision obbtrn(nbfmax,nbfmax,nstate)
      double precision hbbtrnd(nbfmax,nbfmax,nstate)
c
      integer ilm,ibf
      integer g_movecs, g_scat, g_borb
c     
      call ga_sync()
c
      status = rtdb_parallel(.true.)
c
c     First generate the orbitals using SCF/MCSCF
c     
c     if (.not. rtdb_cget(rtdb, 'selci:reference', 1, reference))
c     $     reference = 'scf'
c     if (reference .eq. 'scf') then
      cckohn = scf(rtdb)
c     else
c     cckohn = mcscf(rtdb)
c     endif
c     
c     generate the free functions
c     
      if(.not.rtdb_get(rtdb,'cckohn:lmax',mt_int,1,lmax))
     &     call errquit("can not get lmax from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:nper',mt_int,1,nper))
     &     call errquit("can not get nper from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:xstart',mt_dbl,1,xstart))
     &     call errquit("can not get xstart from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:xmax',mt_dbl,1,xmax))
     &     call errquit("can not get xmax from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:nint',mt_int,1,nint))
     &     call errquit("can not get nint from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:alpha',mt_dbl,1,alpha))
     &     call errquit("can not get alpha from rtdb", 100, RTDB_ERR)
c      
      write(6,*)' Greens Function Code '
      write(6,101)lmax,nper,xstart,xmax,nint,alpha
 101  format(/,' L-Max  = ',i5,/,
     2     ' Nper   = ',i5,/,
     3     ' XStart = ',f20.12,/,
     4     ' XMax   = ',f20.12,/,
     5     ' Nint   = ',i5,/,
     6     ' Alpha  = ',f20.12)
c
      np=nper*xmax+1
      lmaxnp=npd*(1+lbig)
c     
c     real 
c
      if(.not.MA_push_get(mt_dbl,npd,'r', l_r, k_r))
     $     call errquit('cckohn: could not allocate r',npd, MA_ERR)
      call dfill(npd,0d0,dbl_mb(k_r),1)
c
c     complex
c      
      if(.not.MA_push_get(mt_dcpl,lmaxnp,'wwl', l_wwl, k_wwl))
     $     call errquit('cckohn: could not allocate wwl',lmaxnp, MA_ERR)
      call dfill(2*lmaxnp,0d0,dcpl_mb(k_wwl),1)
c
c     real 
c      
      if(.not.MA_push_get(mt_dbl,lmaxnp,'derl', l_derl, k_derl))
     $     call errquit('cckohn: could not allocate derl',lmaxnp,
     &       MA_ERR)
      call dfill(lmaxnp,0d0,dbl_mb(k_derl),1)
c
c     complex
c     
      if(.not.MA_push_get(mt_dcpl,lmaxnp,'csc', l_csc, k_csc))
     $     call errquit('cckohn: could not allocate csc',lmaxnp,
     &       MA_ERR)
      call dfill(2*lmaxnp,0d0,dcpl_mb(k_csc),1)
c     
c     real 
c     
      if(.not.MA_push_get(mt_dbl,lmaxnp,'csl', l_csl, k_csl))
     $     call errquit('cckohn: could not allocate csl',lmaxnp,
     &       MA_ERR)
      call dfill(lmaxnp,0d0,dbl_mb(k_csl),1)

      call ubesplg(rtdb,npd,lmax,nper,nint,alpha,
     &     dbl_mb(k_r),
     &     dcpl_mb(k_wwl),dbl_mb(k_derl),
     &     dcpl_mb(k_csc),dbl_mb(k_csl))
c
c     Build the grid. grid_cckohn will build the grid and do all 
c     deeds that are evaluated on the grid. First zero the matrix 
c     elements
c
      call dfill(2*lmtop*nbfcmx*nchnl,0d0,ovbf,1)
      call dfill(2*lmtop*nbfcmx*nstate,0d0,hpvb,1)
      call dfill(2*lmtop*lmtop*nsymst,0d0,hpvhp,1)
      call dfill(2*lmtop*lmtop*nstate,0d0,hpvhm,1)
c
      call grid_cckohn(rtdb,
     &     dbl_mb(k_r),
     &     dcpl_mb(k_wwl),dbl_mb(k_derl),
     &     dcpl_mb(k_csc),dbl_mb(k_csl),
     &     ovbf,hpvb,hpvhp,hpvhm)
c
      write(*,*)"ovbf"
      write(*,"(6f10.5)")((ovbf(ilm,ibf,1),ilm=1,3),ibf=1,70)
c
      write(*,*)"hpvb"
      write(*,"(6f10.5)")((hpvb(ilm,ibf,1),ilm=1,3),ibf=1,70)
c
      write(*,*)"hpvhp"
      write(*,"(6f10.5)")((hpvhp(ilm,ibf,1),ilm=1,3),ibf=1,3)
c
      write(*,*)"hpvhm"
      write(*,"(6f10.5)")((hpvhm(ilm,ibf,1),ilm=1,3),ibf=1,3)
c
c
c     define the scattering orbitals and evaluate bound-bound integrals
c
      call rohf_scat(rtdb, g_movecs, g_scat, g_borb,obbtrn,hbbtrnd)
c
c     solve the scattering problem
c
      call cknew(rtdb, g_borb, ovbf,hpvb,hpvhp,hpvhm,obbtrn,hbbtrnd)
c
c
c
c      if (.not. selci) return
c
c      call util_print_push
c      call util_print_rtdb_load(rtdb, 'selci')
c      oprint = util_print('information', print_low)
c
c     What is the name of the moints file and what is its format?
c     Should we regenerate it?  If it does not exist make it.
c
c      if (.not. rtdb_cget(rtdb, 'selci:moints', 1, moints)) then
c         if (.not. rtdb_cget(rtdb, 'fourindex:moints file', 1, moints))
c     $        call util_file_name('moints', .true.,.false.,moints)
c      endif
c      if (.not. rtdb_get(rtdb, 'selci:moints:force', mt_log, 1, oforce))
c     $     oforce = .false.
c      if (ga_nodeid() .eq. 0) inquire(file=moints, exist=oexist)
c      call ga_brdcst(303, oexist, mitob(1), 0)
c      if (oforce .or. (.not. oexist)) then
c         if (.not. moints_full(rtdb)) then
c            selci = .false.
c            return
c         endif
c      endif
c      call ga_sync()
c
c     Allocate all available MA space
c
c      lword = ma_inquire_avail(mt_dbl)
c      if (.not. ma_push_get(mt_dbl, lword, 'selci', l_mem, k_mem))
c     $     call errquit('selci: could not allocate all memory',lword)
c
c     tcgmsg vs GA numbering of processes ?
c
c     Get remaining parameters and filenames ...
c     What should we do?  Default is ci+davids from current ciconf.
c     Other alternative is a sequence of selected CI calculations
c     (ci+davids+select) with given sequence of thresholds.
c
c      if (.not. rtdb_cget(rtdb, 'selci:mode', 1, mode))
c     $     mode = 'ci+davids'
c      if (.not. rtdb_get(rtdb, 'selci:selection criterion',
c     $     mt_int, 1, iwpt)) iwpt = 0
c      if (.not. rtdb_get(rtdb, 'selci:nroot',
c     $     mt_int, 1, nroot)) nroot = 1
c      if (.not. rtdb_get(rtdb, 'selci:diag tol',
c     $     mt_dbl, 1, diag_tol)) diag_tol = 1d-6
c      if (.not. rtdb_get(rtdb, 'selci:update', mt_log, 1, oupdate))
c     $     oupdate = .true.
c
c     These files only opened by process 0
c
c      if (.not. rtdb_cget(rtdb, 'selci:ciconf', 1, ciconf))
c     $     call util_file_name('ciconf', .false.,.false., ciconf)
c      if (.not. rtdb_cget(rtdb, 'selci:civecs', 1, civecs))
c     $     call util_file_name('civecs', .false.,.false., civecs)
c      if (.not. rtdb_cget(rtdb, 'selci:wmatrx', 1, wmatrx))
c     $     call util_file_name('wmatrx', .true.,.false., wmatrx)
c      if (.not. rtdb_cget(rtdb, 'selci:eigins', 1, eigins))
c     $     call util_file_name('eigins', .true.,.false., eigins)
c      if (.not. rtdb_cget(rtdb, 'selci:rtname', 1, names(1)))
c     $     call util_file_name('rtname', .false.,.false., names(1))
c      if (.not. rtdb_cget(rtdb, 'selci:roname', 1, names(2)))
c     $     call util_file_name('roname', .false.,.false., names(2))
c      if (.not. rtdb_cget(rtdb, 'selci:hamdg', 1, hamdg))
c     $     call util_file_name('hamdg', .true.,.false., hamdg)
c
c     This file is opened by all processes
c
c      if (.not. rtdb_cget(rtdb, 'selci:hamil', 1, ctmp)) then
c         call util_file_name('hamil', .true.,.true., hamil)
c      else
c         call util_pname(ctmp, hamil)
c      endif
c
c     Should we delete the existing ciconf file?  (not checking this now). 
c     If it does not exisit then generate a default.
c
c      if (ga_nodeid() .eq. 0) inquire(file=ciconf, exist=oexist)
c      call ga_brdcst(304, oexist, mitob(1), 0)
c      if (.not. oexist) call selci_conf(rtdb, dbl_mb(k_mem), lword)
c      
c      if (mode .eq. 'ci+davids') then
c         call selci_ci(dbl_mb(k_mem), lword)
c         call selci_davids(dbl_mb(k_mem), lword, nroot, diag_tol, roots)
c         if (.not. rtdb_put(rtdb, 'selci:roots', mt_dbl, nroot, roots))
c     $        call errquit('selci: failed to store CI roots',nroot)
c         if (.not. rtdb_put(rtdb,'selci:energy',mt_dbl,1,roots(nroot)))
c     $        call errquit('selci: failed to store CI energy',nroot)
c      else if (mode .eq. 'select') then
c         if (.not. rtdb_get_info(rtdb, 'selci:selection thresholds', 
c     $        ma_type, nthresh, date)) then
c            nthresh = 1
c            stols(1) = 1.0d99   ! i.e., don't select any new configuration
c         else
c            if (nthresh .gt. max_thresh) call errquit
c     $           ('selci: too many selection thresholds', max_thresh)
c            if (.not. rtdb_get(rtdb, 'selci:selection thresholds', 
c     $           mt_dbl, max_thresh, stols)) call errquit
c     $           ('selci: no selection thresholds in database?', 0)
c         endif
c         if (.not. rtdb_get(rtdb, 'selci:selection criterion',
c     $        mt_int, 1, iwpt)) iwpt = 0
c         
c         do i = 1, nthresh
c            call selci_ci(dbl_mb(k_mem), lword)
c            call selci_davids(dbl_mb(k_mem),lword,nroot,diag_tol,roots)
c            call selci_select(dbl_mb(k_mem),lword,nroot,iwpt,stols(i),
c     $           oupdate,rootspt)
c     
c            if (i .ne. nthresh) then
c               if (.not. rtdb_put(rtdb, 'selci:selection thresholds', 
c     $              mt_dbl, nthresh-i, stols(i+1))) call errquit
c     $              ('selci: failed writing thresholds to database?', 0)
c            else
c               status = rtdb_delete(rtdb,'selci:selection thresholds')
c            endif
c         enddo
c         if (.not. rtdb_put(rtdb, 'selci:roots', mt_dbl, nroot, roots))
c     $        call errquit('selci: failed to store CI roots',nroot)
c         if (.not. rtdb_put(rtdb, 'selci:roots+pt',
c     $        mt_dbl, nroot, rootspt))
c     $        call errquit('selci: failed to store CI+PT roots',nroot)
c         if (.not. rtdb_put(rtdb, 'selci:energy',mt_dbl,1,
c     $        rootspt(nroot)))
c c    $        call errquit('selci: failed to store CI+PT energy',nroot)
c      else
c         call errquit('selci: unknown mode', 0)
c      endif
c
c     What files need deleting?
c
c      call util_file_unlink(hamil)
c      if (ga_nodeid() .eq. 0) then
c         call util_file_unlink(hamdg)
c         call util_file_unlink(wmatrx)
c      endif
c
c     That's it!
c
c      if (.not. ma_pop_stack(l_mem)) 
c     $     call errquit('selci: ma corrupted', 0)
c      call util_print_pop
c
c      call ga_sync()
c
c     If anything went wrong a hard fail would have resulted
c
      cckohn = .true.
c
      end
