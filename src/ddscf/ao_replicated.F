      subroutine ao_replicated(geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, vg_dens, vg_fock, asym)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "cscfps.fh"
#include "cfock.fh"
#include "util.fh"
#include "stdio.fh"
c     
c     Replicated-data AO 2e-Fock construction routine using the
c     blocked integral interface.
c     
c     F(i) = jfac(i)*J[D(i)] + kfac(i)*K[D(i)]
c     
c     arguments
c     
      integer geom, basis       ! [input] parameter handles
      integer nfock             ! [input] number of Fock matrices
      double precision jfac(nfock) ! [input] Coulomb prefactor
      double precision kfac(nfock) ! [input] exchange prefactor
      double precision tol2e    ! [input] integral selection threshold
      logical oskel             ! [input] toggle skeleton Fock matrix
      integer vg_dens(nfock)    ! [input] array of handles to densities
      integer vg_fock(nfock)    ! [output] array of handles to Fock matrices
      logical asym              ! [input] flag to anisymmetrize Fock matrices
c     
c     local variables
c     
      integer natom, nbf, nsh
      integer i
      integer l_fock, k_fock, l_dens, k_dens, l_rdens, k_rdens
      integer l_q4, k_q4, k_ijkl, l_ijkl
      integer tdim
      parameter (tdim = 100000)  ! Max no. of distinct pair types
      integer l_plist, k_plist, l_ps, k_ps, l_pinfo, k_pinfo
      integer words,memlist(10),lmemlist
      logical status
      logical int2e_test_mem, int2e_buf_write, int2e_file_rewind
      external int2e_test_mem, int2e_buf_write, int2e_file_rewind
      double precision tol2e_local
c
      if (.not. geom_ncent(geom, natom)) call errquit
     $     ('ao_replicated: geom_ncent failed', 0, GEOM_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ao_replicated: bas_numbf failed', 0, BASIS_ERR)
      if ( .not. bas_numcont(basis, nsh) ) call errquit(
     $     'ao_replicated: problem with call to bas_numcont', basis,
     &       BASIS_ERR)
c     
c     If we can hold most of the integrals in core then we should
c     do this ... i.e., the default should be semi-direct. However,
c     if we came in thinking that we were not going to store the
c     integrals then tol2e might be set too high.  In which case
c     reset it (locally)
c     
      tol2e_local = tol2e
      if (.not. (oreadfile .or. owritefile)) then
         owritefile = int2e_test_mem(geom, basis,
     $     min(tol2e*0.1d0,1d-9))
         if (owritefile) tol2e_local = min(tol2e,1d-8)
      end if
c
*      call ga_print(vg_dens)
c
c     If reading/writing to/from cache/file then rewind
c     
      if (oreadfile .or. owritefile) then
         if (.not. int2e_file_rewind()) call errquit
     $        ('ao_replicated: failed rewinding integral file',0,
     &       UNKNOWN_ERR)
      endif
c     
c     Halve exchange factor to conform to internal definition
c     
      call dscal(nfock, 0.5d0, kfac, 1)
c     
c     allocate necessary local temporary arrays on the stack
c     
      lmemlist = 0
      status = .true.
      status = status .and. ma_push_get(MT_DBL, nbf*nbf*nfock, 
     $     'ao_rep: density', l_dens, k_dens)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_dens
      endif
      status = status .and. ma_push_get(MT_DBL, nbf*nbf*nfock, 
     $     'ao_rep: fock',    l_fock, k_fock)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_fock
      endif
      status = status .and. ma_push_get(MT_DBL, nsh*nsh, 
     $     'ao_rep: fock',    l_rdens, k_rdens)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_rdens
      endif
      status = status .and. ma_push_get(MT_DBL, maxquartet,
     $     'ao_rep: q4',    l_q4, k_q4)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_q4
      endif
      status = status .and. ma_push_get(MT_INT, 4*maxquartet,
     $     'ao_rep: ijkl',    l_ijkl, k_ijkl)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_ijkl
      endif
      status = status .and. ma_push_get(MT_INT, nbf*(nbf+1),
     $     'ao_rep: plist',    l_plist, k_plist)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_plist
      endif
      status = status .and. ma_push_get(MT_DBL, nbf*(nbf+1)/2,
     $     'ao_rep: plist',    l_ps, k_ps)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_ps
      endif
      status = status .and. ma_push_get(MT_INT, 2*tdim,
     $     'ao_rep: pinfo',    l_pinfo, k_pinfo)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_pinfo
      endif
      if (.not. status) then
        do i = lmemlist, 1, -1
          if (.not.ma_pop_stack(memlist(i))) then
          endif
        enddo
        if (ga_nodeid().eq.0) then
          write(luout,*)"lmemlist   = ",lmemlist
          write(luout,*)"nbf        = ",nbf
          write(luout,*)"nsh        = ",nsh
          write(luout,*)"nfock      = ",nfock
          write(luout,*)"maxquartet = ",maxquartet
          write(luout,*)"tdim       = ",tdim
          words = 2*nbf*nbf*nfock+nsh*nsh+maxquartet+nbf*(nbf+1)/2+
     $            ma_sizeof(MT_INT,4*maxquartet+nbf*(nbf+1)+2*tdim,
     $                      MT_DBL)+8*ma_sizeof_overhead(MT_DBL)
          write(luout,*)"needing ",words," words or ",
     $                  ma_sizeof(MT_DBL,words,MT_BYTE)," bytes of ",
     $                  "memory"
          write(luout,*)"available are ",ma_inquire_avail(MT_DBL),
     $                  "words or ",ma_inquire_avail(MT_BYTE)," bytes"
        endif
        call errquit
     $     ('ao_replicated: insufficient memory', nfock*nbf*nbf*2,
     &       MA_ERR)
      endif
c     
c     Process 0 gets the density matrices and then broadcasts them
c     (next version should use triangles)
c
c     Load the densities into the space for the Fock matrix and then
c     reorder from (nbf,nbf,nfock) into (nfock,nbf,nbf)
c     
      call ga_sync()
      if (ga_nodeid() .eq. 0) then
         do i = 1, nfock
            call ga_get(vg_dens(i), 1, nbf, 1, nbf, 
     $           dbl_mb(k_fock+(i-1)*nbf*nbf), nbf)
         end do
         call scf_dens_reorder(nfock, nbf,
     $        dbl_mb(k_fock), dbl_mb(k_dens))
      end if
      call ga_sync()
      call ga_brdcst(373, dbl_mb(k_dens), mdtob(nbf*nbf*nfock), 0)
      call ga_sync()
      call dfill(nbf*nbf*nfock, 0.0d0, dbl_mb(k_fock), 1)
c     
      if (util_print('semi_direct',print_debug) .and.
     $     ga_nodeid().eq.0) then
         write(6,*) ' oread, owrite, incmax, incmin, incomplete ',
     $        oreadfile, owritefile, incmax, incmin, incomplete
         call util_flush(6)
      endif
c
      if (oreadfile) then
         if (util_print('semi_direct',print_debug) .and.
     $        ga_nodeid().eq.0) then
            write(6,*) ' Building from file ', tol2e_local
            call util_flush(6)
         endif
         call fock_2e_rep_from_file(geom, basis, nfock, nbf,
     $        jfac, kfac, tol2e_local, oskel,
     $        dbl_mb(k_dens), dbl_mb(k_fock))
      endif
c
c     If (fully_direct .OR. writing_to_disk .OR. semi_direct)
c     need to do a (partial) direct fock build
c
      if ((.not. (oreadfile .or. owritefile)) .or. 
     $    (owritefile) .or. (incmax.gt.-1)) then
         iscreen(6) = iscreen(6) + 1 ! Counts direct fock builds
         if (util_print('semi_direct',print_debug) .and.
     $        ga_nodeid().eq.0) then
            write(6,*) ' Building direct ', tol2e_local
            call util_flush(6)
         endif
c     
         if (util_print('debuginteg',print_never)) then
            open(66,file='/scratch/ints',form='formatted',
     $           status='unknown')
         endif
         call fock_2e_rep( geom, basis, nfock, nbf, nsh, natom,
     $        jfac, kfac, tol2e_local, oskel,
     $        dbl_mb(k_dens), dbl_mb(k_fock), dbl_mb(k_rdens),
     $        int_mb(k_plist), dbl_mb(k_ps), maxquartet,
     $        dbl_mb(k_q4), int_mb(k_ijkl),
     $        tdim, int_mb(k_pinfo))
         if (util_print('debuginteg',print_never)) then
            close(66)
            call errquit('done',0, UNKNOWN_ERR)
         endif
c     
         if (util_print('screening statistics', print_debug)) then
            call schwarz_print(natom, nsh)
         end if
c     
      end if
c     
      call ga_sync()
      call ga_dgop(374, dbl_mb(k_fock), nbf*nbf*nfock, '+')
      call ga_sync()
c     
c     Process 0 PUTs stuff into the global fock matrix
c
c     Reorder from (nfock,nbf,nbf) into (nbf,nbf,nfock)
c     
      if (ga_nodeid() .eq. 0) then
         call scf_fock_reorder(nfock, nbf,
     $        dbl_mb(k_fock), dbl_mb(k_dens))
         do i = 1, nfock
            call ga_put(vg_fock(i), 1, nbf, 1, nbf, 
     $           dbl_mb(k_dens+(i-1)*nbf*nbf), nbf)
         end do
      end if
      status = ma_pop_stack(l_pinfo)
      status = ma_pop_stack(l_ps) .and. status
      status = ma_pop_stack(l_plist) .and. status
      status = ma_pop_stack(l_ijkl)  .and. status
      status = ma_pop_stack(l_q4)    .and. status
      status = ma_pop_stack(l_rdens) .and. status
      status = ma_pop_stack(l_fock)  .and. status
      status = ma_pop_stack(l_dens)  .and. status
      if (.not. status) call errquit('fock_2e: ma_pop?', 0, MA_ERR)
c     
      call ga_sync()
c     
      do i=1,nfock
         call ga_dscal(vg_fock(i), 4.0d0) ! Undo scaling
         if (.not. oskel .and. .not. asym) 
     &        call ga_symmetrize(vg_fock(i))
      end do
c     
c     Reset exchange factors
c     
      call dscal(nfock, 2.d0, kfac, 1)
c     
c     Disable writing integrals to cache/file and enable reading for
c     next call if we were just writing.  Also flush trailing block.
c     
      if (owritefile) then
         oreadfile = .true.
         owritefile = .false.
         if (.not. int2e_buf_write(.true.)) 
     $        call errquit('ao_replicated: failed writing last block',0,
     &       UNKNOWN_ERR)
         if (incmax.eq.-1) call int2e_file_balance()
      end if
c     
      end
      subroutine fock_2e_rep(geom, basis, nfock, nbf, nsh, natom,
     $     jfac, kfac, tol2e, oskel, dens, fock, rdens,
     $     plist, ps, qdim, q4, ijkl, tdim, pinfo)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "cfock.fh"
#include "sym.fh"
#include "cscfps.fh"
#include "util.fh"
#include "mafdecls.fh"
c     
      integer geom, basis
      integer nfock
      integer nbf
      integer nsh
      integer natom
      double precision jfac(nfock), kfac(nfock)
      double precision tol2e
      logical oskel
      double precision dens(nfock, nbf*nbf)
      double precision fock(nfock, nbf*nbf)
      double precision rdens(nsh,nsh)
      integer plist(2,*)        ! Pairs sorted by type
      double precision ps(*)    ! Schwarz screening for pair
      integer qdim
      double precision q4(qdim) ! Q4 factors for each quartet
      integer ijkl(qdim, 4)     ! Indices for quartet
      integer tdim              ! Max no. of pair types (guess)
      integer pinfo(2,tdim)     ! Start and end of each pair type
      integer ntp               ! No. of tasks per processor
c     
      integer natom_unique
      integer nshpair
      integer ntype             ! No. of pair types
      integer max_p             ! Max pairs per block
      integer chunk             ! Chunking factor for tasks
c     
      integer ijklt, ijt, klt, ij, kl, ish, jsh, ksh, lsh
      integer  iat, jat, kat, lat, imat
      double precision q2, qq4, qjunk, sijkl, denmax, dijkl
      double precision dentol, integ_acc
      character*80 buff
c
      double precision jfacmax, kfacmax
      double precision nqraw, nqsch, nqsym
c
      integer nq                ! No. of quartets generated
c
      integer ntreq
      double precision neri_direct, neri_cached
      double precision nqdone, neri_done
      common/fred/nqdone, neri_done, neri_direct, neri_cached, ntreq
c
      integer cnmax
      parameter (cnmax=nw_max_shells) ! Max no. of shells 
      integer cntoce(cnmax)           ! Map from contractions to centers
      integer cntobfr(2,cnmax)        ! Map from contractions to bf range
c
      integer nops              ! No. of group operators = max value of q4
      double precision dnops, dentol_pairs
c
      integer next, current, nproc, incomplete_top
c
      logical oprint, osavewritefile, odotask
      logical odonea
c     
      logical int2e_set_bf_range
      external int2e_set_bf_range
      integer nxtask
      external nxtask
c     
      next = 1          ! take care of compiler warnings
      incomplete_top = 0
c
      if (oscfps) call pstat_on(ps_fock_2e)
c
      oprint = util_print('ao_replicated',print_never)
c
      nops = sym_number_ops(geom) + 1
      dnops = nops
c
      call fock_init_rep_cmul(nbf) ! lookup table for f build
c
      osavewritefile = owritefile
      if (owritefile) then
         if (.not. int2e_set_bf_range(
     $        1,nbf,1,nbf,1,nbf,1,nbf)) call errquit
     $        ('fock_2e_rep: set range failed',0, UNKNOWN_ERR)
         incomplete = -1
      end if
      if (.not. geom_ncent_unique(geom, natom_unique)) call errquit
     $     ('fock_2e_rep: bad geom handle?',0, GEOM_ERR)
      if (nsh .gt. cnmax)   call errquit('fock_2e_rep: cnmax',nsh,
     &       UNKNOWN_ERR)
c     
c     Build maps (for speed) and construct density matrix over atoms
c     
      do ish = 1, nsh
         if (.not. bas_cn2ce(basis, ish, cntoce(ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis,ish,cntobfr(1,ish),cntobfr(2,ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do
c
      call dfill(nsh*nsh,0.0d0,rdens,1)
      jfacmax = 0.0d0
      kfacmax = 0.0d0
      do imat = 1, nfock
         jfacmax = max(jfacmax,abs(jfac(imat)))
         kfacmax = max(kfacmax,abs(kfac(imat)))
         call util_mat_reduce(nbf, nsh, cntobfr, dens(imat,1), 
     $        nfock, rdens, 'rms')
      end do
      denmax = 0.0d0
      do ish = 1, nsh
         do jsh = 1, ish
            denmax = max(denmax,rdens(jsh,ish))
         end do
      end do
      if (oprint) then
         write(6,*) ' JFACMAX ', jfacmax
         write(6,*) ' KFACMAX ', kfacmax
         write(6,*) '  DENMAX ', denmax
***         write(6,*) '   RDENS'
***         call doutput(rdens,1,nsh, 1, nsh, nsh, nsh,1)
      endif

*      if (owritefile) denmax = max(denmax,10.0d0) ! If do this not consistent
c
      dentol = min(dentolmax,tol2e/denmax,tol2e/denmax**2) ! Threshold to screen integs only
c
c     If this is the first call for a potential semi-direct calculation 
c     then we must save the tol2e since the definition of the tasks 
c     depends on the tolerance and the tasks must be identical.
c
      if (owritefile) semi_direct_tol2e = dentol
c     
c     There is a BIG speed gain from aggressive reduction of accuracy
c     in the actual computation of the integrals. However, the screening
c     on primitive integrals is done on radial prefactors not using
c     the Schwarz inequality, hence the extra factor of 0.01d0.
c     (note that we also have at least another 1/10 from tol2e->dentol)
c
c     Michel uses a factor of 1d-6 in HONDO (compared to our 10^-3) ->
c     more experiments are needed.
c
c     C60 has trouble with max value for integ_acc less than 1d-9
c     (using texas integrals, do not know about nwchem)
c
c     More testing shows that NWchem integrals can be wrong with
c     an accuracy coarser than 1d-15 which is now the standard accuracy.
c
c     Thus now use the default integral accuracy unless set by the
c     user.  Sigh.
c     
      if (intacc .eq. 0.0d0) then
cedo         integ_acc = 1d-15 ! Not used
         integ_acc = 1d-11 !   This seems fine with txs EA 1/2003
         continue
*         integ_acc = min(1d-15,max(0.01d0*dentol,1d-30)) ! Variable
*         call int_acc_set(integ_acc)
      else
         integ_acc = intacc     ! User controlled
         call int_acc_set(integ_acc)
      end if
      if (util_print('ao_rep_info', print_debug) .and.
     $       ga_nodeid().eq.0) then
         write(6,111) tol2e, denmax, dentol, integ_acc, odensityscreen,
     $        oreadfile, owritefile
 111     format(1x,1p,' tol2e=',d8.1,' denmax=',d8.1,
     $        ' dentol=',d8.1,' integ_acc=',d8.1,
     $        ' odenscr=',l1,' oread=',l1,' owrite=',l1)
      endif
c     
c     assume a dense integral list but with symmetry and figure out
c     no. of pairs per pair block dividing quartets between processors
c     and allowing for the max no. of quartets to be processed per block
c     
      ntp = 10
      nshpair = nsh*(nsh+1)/2
      if (ga_nnodes() .gt. 1) then
         max_p = int(dble(nshpair)*dsqrt(dnops))
         max_p = max_p/sqrt(dble(ntp)*ga_nnodes())
         max_p = min(nint(sqrt(dble(nops*qdim))),max_p)
         max_p = max(10,max_p)
      else
         max_p = int(sqrt(dble(nops*qdim)))
      end if
c     
c     generate info about the pairs
c
c     If not first call of semi-direct then MUST reuse the 
c     old tol2e so that end up with consistent task definitions
c     
      dentol_pairs = dentol/dnops
      if (oreadfile .and. incmax.ne.-1) 
     $     dentol_pairs = semi_direct_tol2e/dnops
c
      odonea = .false.          ! To avoid looping in block adjust logic
 10   call fock_pairs(basis, max_p, .true., 1, natom, 1, natom, 
     $     dentol_pairs, plist, ps, pinfo, tdim, ntype)
      if (ga_nodeid() .eq. 0) call util_flush(6)
c
c     Figure out no. of tasks (chunk) to do on each request to the
c     shared counter.  Memory permitting we want this to be close to 1.
c
c     Also, try to make sure we have enough tasks for load-balance.
c
      chunk = max(1, (ntype*ntype/2)/(ntp*ga_nnodes()))
      if (max_p .lt. nint(sqrt(dble(nops*qdim))) .and.
     $     chunk .gt. 4) then
         max_p = min(nint(sqrt(dble(nops*qdim))),2*max_p)
         if (oprint .and. ga_nodeid().eq.0) write(6,*) ' A ', max_p
         odonea = .true.
         goto 10
      else if ((ntype*ntype/2)/ga_nnodes().lt.ntp  .and.
     $        max_p.gt.10 .and. (.not. odonea)) then
         max_p = max(10,max_p/2)
         if (oprint .and. ga_nodeid().eq.0) write(6,*) ' B ', max_p
         goto 10
      end if
      if ((oprint  .or. util_print('semi_direct',print_debug))
     $     .and. ga_nodeid().eq.0) then
         write(6,*) ' max_p ', max_p, ' ntype ', ntype,
     $     ' chunk ',chunk
         call util_flush(6)
      end if
c     
c     loop over pair blocks and compute integrals ... want most expensive
c     integrals first, hence outer loop over ijklt ... assume that
c     fock_pairs orders most expensive pairs last.
c     
c     We could constrain the pair loops to run over a wierd unique
c     list of quartets, but we need the conventional 'canonical'
c     unique list to use petite-list symmetry.
c     
      nproc = ga_nnodes()
c
      if (oreadfile .and. incomplete.eq.-1) 
     $     goto 12321           ! Semidirect ... I have none to compute
c
      if (oreadfile) then
c     
c     Semi-direct.  This process filled up the file at task incomplete
c     but task chunking means it did other tasks that
c     would not have been assigned to another process.  This process
c     must do tasks in the range incomplete...incomplete_top
c
         incomplete_top = (incomplete/chunk + 1)*chunk - 1
      else
         next = nxtask(nproc, chunk)
      endif
      current = 0
      nq = 0
c
      nqraw = 0
      nqsch = 0
      nqsym = 0
c
      nqdone = 0
      neri_done = 0
      ntreq  = 1
      neri_direct = 0
      neri_cached = 0
c
      do ijklt = 2*ntype, 2, -1
         do ijt = min(ntype,ijklt-1), max(1,ijklt-ntype), -1
            klt = ijklt - ijt
            if (oreadfile .and. current.le.incmax) then
               odotask = (current.ge.incomplete) .and. 
     $              (current.le.incomplete_top)
            else  
               odotask = next.eq.current
            endif
            if (odotask) then
c
               if (owritefile) call int2e_file_record_position()
               do ij = pinfo(1,ijt), pinfo(2,ijt)
                  ish = plist(1,ij)
                  jsh = plist(2,ij)
                  iat = cntoce(ish)
                  jat = cntoce(jsh)
                  q2  = 1.0d0
                  if (ish .eq. jsh) q2 = 0.5d0
                  do kl = pinfo(1,klt),pinfo(2,klt)
                     qq4 = q2
                     ksh = plist(1,kl)
                     lsh = plist(2,kl)
                     kat = cntoce(ksh)
                     lat = cntoce(lsh)
c     enforce canonical indices
                     if (ish.lt.ksh) then
                        goto 200 ! Next kl
                     else if (ish.eq.ksh) then
                        if (jsh.lt.lsh) then
                           goto 200 ! Next kl
                        else if (jsh.eq.lsh) then
                           qq4 = qq4 * 0.5d0 ! i==k && j==l diagonal
                        end if
                     end if
c     kl diagonal
                     if (ksh.eq.lsh) qq4 = qq4*0.5d0
c
                     nqraw = nqraw + 1.0d0
c     schwarz screening
                     sijkl = ps(ij)*ps(kl)
                     if (sijkl*dnops*qq4 .lt. dentol) goto 200 ! Next kl
                     nqsch = nqsch + 1.0d0
c     symmetry test
                     if (oskel) then
                        if (.not. sym_atom_quartet(geom,iat,jat,kat,lat,
     $                       qjunk)) goto 200
                        qq4 = qq4*qjunk
                     end if
                     nqsym = nqsym + 1.0d0
c     density screening
                     if (odensityscreen .and. (.not. owritefile)) then
                        dijkl = jfacmax*max(rdens(ish,jsh),
     $                       rdens(ksh,lsh)) + 
     $                       kfacmax*max(rdens(ish,ksh),rdens(ish,lsh), 
     $                       rdens(jsh,ksh),rdens(jsh,lsh))
                        dijkl = max(dijkl,dijkl*dijkl)
                        if (qq4*sijkl*dijkl .lt. tol2e) goto 200 ! Next kl
                     end if
c     save it for later
                     nq = nq + 1
                     ijkl(nq,1) = ish
                     ijkl(nq,2) = jsh
                     ijkl(nq,3) = ksh
                     ijkl(nq,4) = lsh
                     q4(nq)     = qq4
                     if (nq .eq. qdim) then
                        nqdone = nqdone + nq
                        call fock_rep_txs(basis, nfock, nbf,
     $                       jfac, kfac, dentol, dens, fock, 
     $                       nq, ijkl(1,1), ijkl(1,2), ijkl(1,3),
     $                       ijkl(1,4), q4, current)
                        nq = 0
                     end if
 200              end do        ! Next KL pair
 100           end do           ! Next IJ pair
c     
c     Empty out interactions for this type
c     
               if (nq .gt. 0) then
                  nqdone = nqdone + nq
                  call fock_rep_txs(basis, nfock, nbf,
     $                 jfac, kfac, dentol, dens, fock, 
     $                 nq, ijkl(1,1), ijkl(1,2), ijkl(1,3),
     $                 ijkl(1,4), q4, current)
                  nq = 0
               end if
c
c     Get next task ... messy because of semidirect
c
               if (oreadfile) then ! Semi-direct
                  if (current .ge. incomplete_top) then
                     if (incmin.eq.-1) then
                        goto 12321 ! Did my only task
                     endif
                     next = nxtask(nproc,chunk)+incmax+1
                     ntreq = ntreq + 1
                  endif
               else
                  next = nxtask(nproc,chunk)
                  ntreq = ntreq + 1
               endif
            end if
            current = current + 1
         end do                 ! Next IJ type
      end do                    ! Next IJKL type
c
12321 owritefile = osavewritefile
      if (owritefile) then
         if (incomplete.eq.-1) then
            call ga_sync        ! Since did not in fock_rep_txs
            incmax = -1
            incmin = -1
         else
            incmax = (incomplete/chunk + 1)*chunk - 1 ! Round to task top
            incmin = incmax
         endif
         call ga_igop(123, incmax, 1, 'max')
         call ga_igop(321, incmin, 1, 'min')
c
      endif
c     
      if (oscfps) call pstat_off(ps_fock_2e)
      next = nxtask(-nproc,1)
c
      if (oprint) then
         call begin_seq_output
         buff = ' '
         write(buff,99) ga_nodeid(), ntreq, nqdone, neri_done
 99      format(i4,': tasks=',i4,1p,': nqdone=',d9.1,': neri=',d9.1)
         call write_seq(6,buff)
         call end_seq_output
      end if
c     
      call ga_dgop(11,nqdone,4,'+')
      neri_direct = neri_done - neri_cached
      if (owritefile .and. util_print('ao_replicated',print_default)
     $     .and. ga_nodeid().eq.0) then
         neri_direct = 100.0*neri_direct/neri_done
         neri_cached = 100.0*neri_cached/neri_done
         write(6,109) nqdone, neri_done, neri_direct, neri_cached
 109     format(/1p,' #quartets =',d10.3, ' #integrals =',d10.3,
     $        0p,' #direct =',f5.1,'% #cached =',f5.1,'%'/)
         call util_flush(6)
      endif
c
      if (util_print('ao_rep_info', print_debug)) then
         write(6,77) nqraw, nqsch, nqsym, nqdone
 77      format(1p,' raw=',d8.1,' schwarz=',d8.1,' sym=',d8.1,
     $        ' dens=',d8.1)
      endif
c
c     done
c     
      call int_acc_std()
c     
      end
      subroutine fock_rep_txs(basis, nfock, nbf,
     $     jfac, kfac, tol2e, dens, fock, 
     $     nq, iq, jq, kq, lq, q4, current)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cfock.fh"
#include "util.fh"
      integer basis
      integer nfock, nbf
      double precision jfac(nfock), kfac(nfock)
      double precision tol2e
      double precision dens(nfock, nbf*nbf)
      double precision fock(nfock, nbf*nbf)
      integer nq
      integer iq(nq), jq(nq), kq(nq), lq(nq)
      double precision q4(nq)
      integer current           ! Current task number for semi-direct
c     
c     Given a list of integral quartets, notionally all of the
c     same type, compute the integrals and add them into the
c     fock matrix, optionally caching them on disk/memory.
c     
      integer lenscr, leneri, lenlab
      integer l_eri, k_eri, l_lab, k_lab, l_scr, k_scr
      integer k_i, k_j, k_k, k_l 
*      integer kkk
c     
      integer neri, ntries
      double precision block_eff
      logical more
c
      integer maxg
c
      integer ntreq
      double precision neri_direct, neri_cached
      double precision nqdone, neri_done
      common/fred/nqdone, neri_done, neri_direct, neri_cached,ntreq
c     
      logical int2e_file_write, intb_2e4c, intb_init4c
      external int2e_file_write, intb_2e4c, intb_init4c
c     
      if (nq .le. 0) return
c
c     Figure out dimensions for scratch and eri buffer space.
c     The eri buffer space is defaulted in scf_get_info to some
c     suitable non-zero number.  The scratch space defaults to zero
c     in which case int_mem_2e4c is used, otherwise the user input
c     is taken.
c
      call intb_mem_2e4c(maxg, lenscr) ! blocking algorithm
      if (maxscr .ne. 0) then
         lenscr = maxscr
      end if
      leneri = max(maxeri,maxg)
      lenlab = 4*leneri
c
      if (util_print('fock_rep_txs', print_never)) then
         write(6,11) nq, lenscr, leneri
 11      format(' nq=',i6,', lenscr=',i8,', leneri=',i8)
         call util_flush(6)
      end if
c
c     Allow two goes at reducing memory before fail completely
c
c     Ooops ... cannot reduced the texas scratch memory.  Only
c     try reducing the ERI and LAB memory
c
      ntries = 1                ! No. of tries at allocating memory
 12   if (ntries.eq.1) then
         continue
      else if (ntries.le.3) then
****         lenscr = lenscr/2
         leneri = leneri/2
         lenlab = 4*leneri
      else
         call errquit('fock_rep_txs: insufficient memory ',
     $        lenscr+leneri+lenlab, MEM_ERR)
      endif
      ntries = ntries + 1
c 
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri)) then
         goto 12
      endif
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr)) then
         if (.not. ma_pop_stack(l_eri))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         goto 12
      endif
      if (.not. ma_push_get(mt_int, lenlab, 'lab', l_lab, k_lab)) then
         if (.not. ma_pop_stack(l_scr))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         if (.not. ma_pop_stack(l_eri))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         goto 12
      endif
c     
      k_i = k_lab
      k_j = k_i   + leneri
      k_k = k_j   + leneri
      k_l = k_k   + leneri
c     
      if (.not. intb_init4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., lenscr, dbl_mb(k_scr), leneri, block_eff))
     $     call errquit('intb_init?',nq, UNKNOWN_ERR)
c     
 50   continue
c
      more = intb_2e4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., tol2e, .false.,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), leneri, neri, lenscr, dbl_mb(k_scr))
*
*
*      do kkk = 0, neri-1
*         if (int_mb(k_i+kkk).le.0 .or. int_mb(k_i+kkk).gt.nbf)
*     $        call errquit('bad i label ', kkk)
*         if (int_mb(k_j+kkk).le.0 .or. int_mb(k_j+kkk).gt.nbf)
*     $        call errquit('bad j label ', kkk)
*         if (int_mb(k_k+kkk).le.0 .or. int_mb(k_k+kkk).gt.nbf)
*     $        call errquit('bad k label ', kkk)
*         if (int_mb(k_l+kkk).le.0 .or. int_mb(k_l+kkk).gt.nbf)
*     $        call errquit('bad l label ', kkk)
*      enddo
*
*
c     
      if (neri .gt. 0) then
c
         neri_done = neri_done + neri
c     
c     Cache integrals if desired
c     
         if (owritefile) then
            if (.not. int2e_file_write(neri,
     $           int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $           dbl_mb(k_eri))) then
               call int2e_file_reposition_truncate()
               owritefile = .false.
               incomplete = current
               call ga_sync()
            else
               neri_cached = neri_done
            endif
         end if
c     
         if (util_print('debuginteg',print_never)) then
            call print_integ_list(neri,
     $           int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $           dbl_mb(k_eri))
         endif
c     
         call fock_2e_rep_label(nfock, nbf, jfac, kfac, tol2e, neri,
     $        int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $        dbl_mb(k_eri), dens, fock)
      endif
c     
      if (more) then
         if (util_print('fock_rep_txs', print_never)) then
            write(6,*) ' Texas has split the request '
         end if
         goto 50
      end if
c     
 100  if (.not. ma_pop_stack(l_lab))
     $     call errquit('fock_rep_txs: cannot free lab',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('fock_rep_txs: cannot free scr',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_eri))
     $     call errquit('fock_rep_txs: cannot free eri',lenlab, MA_ERR)
c     
      end
      subroutine fock_pairs(basis, max_p, oij, 
     $     iatlo, iathi, jatlo, jathi, 
     $     tol2e, ijlist, sij, ijinfo, max_ij, nt_ij)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "util.fh"
#include "global.fh"
      integer basis
      integer max_p, max_ij
      logical oij
      integer iatlo, iathi, jatlo, jathi
      integer ijlist(2,*), ijinfo(2,max_ij), nt_ij
      double precision sij(*)
c     
      integer junk, iat, jat, iatsh, jatsh, shlo, shhi
      integer nij, ij, ijt, i, j, jtop, spch
      integer jlo, jhi, jnprim, jngen, jtype
      integer ilo, ihi, inprim, ingen, itype
      double precision smax, test, tol2e
      integer n_in_b, old_nt
      integer iuat, juat        ! Unique atom numbers
c
      double precision price
      double precision int_2e4c_price
      external int_2e4c_price
      double precision thetype, it, jt, maxtype, prev
      integer iprice
#include "itri.fh"
c
      iprice(i,j) = nint(100.0d0*int_2e4c_price(basis,i,j,basis,i,j))
      thetype(iuat,iatsh,itype,inprim,ingen) = 
     $     iuat-1 + 32.0d0*(iatsh-1 + 32.0d0*(inprim-1 + 32.0d0*(
     $     ingen-1 + 32.0d0*(itype+1))))
c
c     Max value of the above function will be about 10*32**4 = 1.05d+07
c
      maxtype = 1.05d7
c
c     We want to block pairs of shells so that within each block
c     all pairs are of the same type (contraction level, angular
c     momentum, general contraction, and for Texas 95 the exponents
c     and coeffs need to be the same also).  Blocks need to be
c     less than a certain size and also if blocks are too small 
c     join them together since this further reduces precomputation cost.
c
c     Furthermore, for good load balance, we want the most expensive
c     shell-quartets computed first.  These will probably be very
c     high-angular momentum shells, or deep general contractions.
c     But, for best gain from semidirect we want the most expensive
c     single integrals computed first, but these will probably be
c     just the heavily contracted S and P shells.
c
c     Resolve this by constructing the price per integral between 0 and 1.
c     This typically has a range between 1e-4 and 1.0d0.
c     Scaling by 100 to make an integer, 0 - 100.  Use this
c     as the leading value in the sort criterion.  The most expensive 
c     single integrals will have a non-zero cost and will be done first,
c     but most others will have a zero integer cost and will be sorted
c     by other criteria. And we will make these in the order of angular
c     momentum.
c
c     Fun and games here to handle centers with no functions
c
      do iat = iatlo, iathi
         if (.not. bas_ce2cnr(basis, iat, ilo, junk)) call errquit
     $        ('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (ilo .gt. 0) goto 10
      enddo
 10   do iat = iathi,iatlo,-1
         if (.not. bas_ce2cnr(basis, iat, junk, ihi)) call errquit
     $        ('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (ihi .gt. 0) goto 20
      enddo
 20   do jat = jatlo, jathi
         if (.not. bas_ce2cnr(basis, jat, jlo, junk)) call errquit
     $        ('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (jlo .gt. 0) goto 30
      enddo
 30   do jat = jathi,jatlo,-1
         if (.not. bas_ce2cnr(basis, jat, junk, jhi)) call errquit
     $        ('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (jhi .gt. 0) goto 40
      enddo
 40   continue
c     
      smax = schwarz_max()
c     
      nt_ij = 0
      nij = 0
      do i = ilo, ihi
         if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (.not. bas_cn2uce(basis,i,iuat))
     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (.not. bas_cn2ce(basis,i,iat))
     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
         if (.not. bas_ce2cnr(basis,iat,shlo,shhi))
     $        call errquit('fock_pairs: basis corrupt ',0, BASIS_ERR)
         iatsh = i-shlo+1
c     reorder L,s,p (-1,0,1) into order by complexity s,p,L (-1,0,1)
         if (itype.eq.-1) then
            itype = 1
         else if (itype.eq.0) then
            itype = -1
         else if (itype.eq.1) then
            itype = 0
         end if
         it = thetype(iuat,iatsh,itype,inprim,ingen)
         jtop = jhi
         if (oij) jtop = i
         do j = jlo, jtop
            test = schwarz_shell(i,j)
            if (test*smax .gt. tol2e) then
               if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               if (.not. bas_cn2uce(basis,j,juat))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               if (.not. bas_cn2ce(basis,j,jat))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               if (.not. bas_ce2cnr(basis,jat,shlo,shhi))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               jatsh = j-shlo+1
c     reorder L,s,p (-1,0,1) into order by complexity s,p,L
               if (jtype.eq.-1) then
                  jtype = 1
               else if (jtype.eq.0) then
                  jtype = -1
               else if (jtype.eq.1) then
                  jtype = 0
               end if
               jt = thetype(juat,jatsh,jtype,jnprim,jngen)
               nij           = nij + 1
               ijlist(1,nij) = i
               ijlist(2,nij) = j
*               write(6,910) nij, it, jt, iprice(i,j)
* 910           format(1x,i8,1p,2d9.2,0p,i8)
               sij(nij) =  min(it,jt) +
     $              maxtype*(max(it,jt) + maxtype*dble(iprice(i,j)))
            end if
         end do
      end do
      if (nij .eq. 0) return
c     
      call fock_pair_sort(nij,ijlist,sij)
c     
c     Split pair blocks on type boundaries with regard to size
c     (first do type only then combine)
c     
      nt_ij = 0
      prev = -1.0d0
      n_in_b = 0
      do ij = 1, nij
         if (sij(ij).ne.prev .or. n_in_b.ge.max_p) then
            if (nt_ij .ne. 0) ijinfo(2,nt_ij) = ij-1
            nt_ij = nt_ij + 1
            if (nt_ij .gt. max_ij) call errquit
     $           ('fock_pairs: dimension failure', nt_ij, UNKNOWN_ERR)
            ijinfo(1,nt_ij) = ij
            prev = sij(ij)
            n_in_b = 0
         end if
         n_in_b = n_in_b + 1
      end do
      ijinfo(2,nt_ij) = nij
c
      if (util_print('fock pairs', print_debug)
     $     .and. ga_nodeid().eq.0) then
         write(6,*) ' Before combining, #pair blocks = ', nt_ij
         do ijt = 1, nt_ij
            write(6,*) ' ij type ', ijt, ijinfo(1,ijt),ijinfo(2,ijt)
         end do
      end if
c     
      old_nt = nt_ij
      nt_ij = 1
      n_in_b = (ijinfo(2,1)-ijinfo(1,1)+1)
      do ijt = 2, old_nt
         n_in_b = n_in_b + (ijinfo(2,ijt)-ijinfo(1,ijt)+1)
         if (n_in_b .le. max_p) then
            ijinfo(2,nt_ij) = ijinfo(2,ijt)
         else
            nt_ij = nt_ij + 1
            if (nt_ij .gt. max_ij) call errquit
     $           ('fock_pairs: dimension failure(2)', nt_ij,
     &       UNKNOWN_ERR)
            ijinfo(1,nt_ij) = ijinfo(1,ijt)
            ijinfo(2,nt_ij) = ijinfo(2,ijt)
            n_in_b = (ijinfo(2,ijt)-ijinfo(1,ijt)+1)
         end if
      end do
      ijinfo(2,nt_ij) = nij
c     
c     Gather schwarz info again
c     
      do ij = 1, nij
         i = ijlist(1,ij)
         j = ijlist(2,ij)
         sij(ij) = schwarz_shell(i,j)
      end do
c     
      if (util_print('fock pairs', print_debug) .and.
     $     ga_nodeid().eq.0) then
         write(6,*) ' After combining, #pair blocks = ', nt_ij
         do ijt = 1, nt_ij
            write(6,*) ' ij type ', ijt, ijinfo(1,ijt),ijinfo(2,ijt)
            write(6,*) '    ij    i  j   it jt   ip jp   ig jg  iu  ju'
            write(6,*) '    --    -  -   -- --   -- --   -- --  --  --'
            do ij = ijinfo(1,ijt),ijinfo(2,ijt)
               i = ijlist(1,ij)
               j = ijlist(2,ij)
               price = int_2e4c_price(basis, i, j, basis, i, j)
              if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
     $              call errquit('fock_pairs: basis bad ',0,
     &       BASIS_ERR)
               if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
     $              call errquit('fock_pairs: basis bad ',0,
     &       BASIS_ERR)
               if (.not. bas_cn2uce(basis,i,iuat))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               if (.not. bas_cn2uce(basis,j,juat))
     $              call errquit('fock_pairs: basis corrupt ',0,
     &       BASIS_ERR)
               write(6,11) ij,i,j,itype,jtype,inprim,jnprim,
     $              ingen,jngen, iuat, juat , price, iprice(i,j)
 11            format(1x,i6,2x,2i3,2x,2i3,2x,2i3,2x,2i3,2x,2i3,2x,
     $              1p,d9.2,2x,i4)
            end do
         end do
      end if
c     
      end
      subroutine fock_pair_sort(n,ijlist,ijtype)
      implicit integer (a-z)
      dimension ijlist(2,n)
      double precision ijtype(n)
      double precision rra
c     
      if (n .eq. 1) return
c     
      l=n/2+1
      ir=n
 10   continue
      if(l.gt.1)then
         l=l-1
         rra=ijtype(l)
         rrb1=ijlist(1,l)
         rrb2=ijlist(2,l)
      else
         rra=ijtype(ir)
         rrb1=ijlist(1,ir)
         rrb2=ijlist(2,ir)
         ijtype(ir)=ijtype(1)
         ijlist(1,ir)=ijlist(1,1)
         ijlist(2,ir)=ijlist(2,1)
         ir=ir-1
         if(ir.eq.1)then
            ijtype(1)=rra
            ijlist(1,1)=rrb1
            ijlist(2,1)=rrb2
            return
         end if
      end if
      i=l
      j=l+l
 20   if(j.le.ir)then
         if(j.lt.ir)then
            if(ijtype(j).lt.ijtype(j+1))j=j+1
         end if
         if(rra.lt.ijtype(j))then
            ijtype(i)=ijtype(j)
            ijlist(1,i)=ijlist(1,j)
            ijlist(2,i)=ijlist(2,j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         go to 20
      end if
      ijtype(i)=rra
      ijlist(1,i)=rrb1
      ijlist(2,i)=rrb2
      go to 10
      end
      subroutine util_mat_reduce(n, nr, map, a, lda, r, op)
      implicit none
#include "errquit.fh"
c     
      integer n                 ! Original size [input]
      integer nr                ! Reduced size  [input]
      integer map(2,nr)         ! map(1,*)=lo, map(2,*)=hi [input]
      integer lda               ! leading dimension of A (1 if square)
      double precision a(lda, n,n)   ! Original matrix [input]
      double precision r(nr,nr) ! Reduced matrix [output]
      character*(*) op          ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, i, j
      double precision sum
c     
*     write(6,*) ' util_mat_reduce: input matrix '
*     call output(a, 1, n, 1, n, n, n, 1)
c     
      do jr = 1, nr
         do ir = 1, nr
            sum = 0.0d0
            if (op .eq. 'abssum') then
               do j = map(1,jr), map(2,jr)
                  do i = map(1,ir), map(2,ir)
                     sum = sum + abs(a(1,i,j))
                  end do
               end do
            else if (op .eq. 'absmax') then
               do j = map(1,jr), map(2,jr)
                  do i = map(1,ir), map(2,ir)
                     sum = max(sum, abs(a(1,i,j)))
                  end do
               end do
            else if (op .eq. 'rms') then
               do j = map(1,jr), map(2,jr)
                  do i = map(1,ir), map(2,ir)
                     sum = sum + a(1,i,j)*a(1,i,j)
                  end do
               enddo
               sum = sqrt(sum)
            else
               call errquit('util_mat_reduce: unknown op', 0,
     &       UNKNOWN_ERR)
            end if
            r(ir,jr) = max(r(ir,jr),sum)
         end do
      end do
c     
*     write(6,*) ' util_mat_reduce: reduced matrix '
*     call output(r, 1, nr, 1, nr, nr, nr, 1)
      end
      subroutine scf_dens_reorder(nfock, nbf, old, new)
      implicit none
      integer nfock, nbf
      double precision old(nbf*nbf,nfock), new(nfock,nbf*nbf)
c
      integer i, j
c
      do i = 1, nfock
         do j = 1, nbf*nbf
            new(i,j) = old(j,i)
         end do
      end do
c
      end
      subroutine scf_fock_reorder(nfock, nbf, old, new)
      implicit none
      integer nfock, nbf
      double precision old(nfock,nbf*nbf), new(nbf*nbf,nfock)
c
      integer i, j
c
      do i = 1, nfock
         do j = 1, nbf*nbf
            new(j,i) = old(i,j)
         end do
      end do
c
      end
      double precision function int_2e4c_price(
     $     bra, ish, jsh, ket, ksh, lsh)
      implicit none
c
      integer bra, ish, jsh, ket, ksh, lsh
c
      double precision price
c
      logical cando_txs
      external cando_txs
c
      if (cando_txs(bra, ish, jsh) .and. 
     $    cando_txs(ket, ksh, lsh)) then
         call get_int_price(ish,jsh,ksh,lsh,price)
         int_2e4c_price = price
      else
         int_2e4c_price = 0.0d0
      endif
c
      end

