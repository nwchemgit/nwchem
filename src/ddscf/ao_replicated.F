      subroutine ao_replicated( geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, vg_dens, vg_fock )
c     
c     $Id: ao_replicated.F,v 1.4 1996-07-03 17:11:35 d3g681 Exp $
c     
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "cscfps.fh"
#include "cfock.fh"
#include "util.fh"
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
c     
c     local variables
c     
      integer natom, nbf, nsh
      integer i
      integer l_fock, k_fock, l_dens, k_dens, l_rdens, k_rdens
      integer l_q4, k_q4, k_ijkl, l_ijkl
      integer l_plist, k_plist, l_ps, k_ps
      logical status
      logical int2e_test_mem
      external int2e_test_mem
      double precision tol2e_local
c
      if (.not. geom_ncent(geom, natom)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ao_replicated: bas_numbf failed', 0)
c     
c     If we can hold all of the integrals in core then we should
c     do this ... i.e., the default should be semi-direct, though
c     right now we must store all of the integrals.  However, 
c     if we came in thinking that we were not going to store the
c     integrals then tol2e might be set too high.  In which case
c     reset it (locally)
c     
      tol2e_local = tol2e
      if (.not. (oreadfile .or. owritefile)) then
         owritefile = int2e_test_mem(geom, basis,
     $     min(tol2e*0.1d0,1d-8))
         if (owritefile) tol2e_local = min(tol2e,1d-8)
      end if
c     
c     If reading/writing to/from cache/file then rewind
c     
      if (oreadfile .or. owritefile) call int2e_file_rewind
c     
c     Halve exchange factor to conform to internal definition
c     
      call dscal(nfock, 0.5d0, kfac, 1)
c     
c     allocate necessary local temporary arrays on the stack
c     
      status = .true.
      status = status .and. ma_push_get(MT_DBL, nbf*nbf*nfock, 
     $     'ao_rep: density', l_dens, k_dens)
      status = status .and. ma_push_get(MT_DBL, nbf*nbf*nfock, 
     $     'ao_rep: fock',    l_fock, k_fock)
      status = status .and. ma_push_get(MT_DBL, natom*natom, 
     $     'ao_rep: fock',    l_rdens, k_rdens)
      status = status .and. ma_push_get(MT_DBL, maxquartet,
     $     'ao_rep: q4',    l_q4, k_q4)
      status = status .and. ma_push_get(MT_INT, 4*maxquartet,
     $     'ao_rep: ijkl',    l_ijkl, k_ijkl)
      status = status .and. ma_push_get(MT_INT, nbf*(nbf+1),
     $     'ao_rep: plist',    l_plist, k_plist)
      status = status .and. ma_push_get(MT_DBL, nbf*(nbf+1)/2,
     $     'ao_rep: plist',    l_ps, k_ps)
      if (.not. status) call errquit
     $     ('ao_replicated: insufficient memory', nfock*nbf*nbf*2)
c     
c     Process 0 gets the density matrices and then broadcasts them
c     (next version should use triangles)
c
c     Load the densities into the space for the Fock matrix and then
c     reorder from (nbf,nbf,nfock) into (nfock,nbf,nbf)
c     
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
      if (oreadfile) then
         call fock_2e_rep_from_file(geom, basis, nfock, nbf,
     $        jfac, kfac, tol2e_local, oskel,
     $        dbl_mb(k_dens), dbl_mb(k_fock))
      else
c     
         iscreen(6) = iscreen(6) + 1 ! Counts direct fock builds
c     
*         open(66,file='/scratch1/ints',form='formatted',
*     $        status='unknown')
         call fock_2e_rep( geom, basis, nfock, nbf, natom,
     $        jfac, kfac, tol2e_local, oskel,
     $        dbl_mb(k_dens), dbl_mb(k_fock), dbl_mb(k_rdens),
     $        int_mb(k_plist), dbl_mb(k_ps), maxquartet,
     $        dbl_mb(k_q4), int_mb(k_ijkl))
*         close(66)
*         call errquit('done',0)
c     
         if (util_print('screening statistics', print_debug)) then
            if ( .not. bas_numcont(basis, nsh) ) call errquit(
     $           'ao_fock_2e: problem with call to bas_numcont', basis)
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
      status = ma_pop_stack(l_ps)
      status = ma_pop_stack(l_plist) .and. status
      status = ma_pop_stack(l_ijkl)  .and. status
      status = ma_pop_stack(l_q4)    .and. status
      status = ma_pop_stack(l_rdens) .and. status
      status = ma_pop_stack(l_fock)  .and. status
      status = ma_pop_stack(l_dens)  .and. status
      if (.not. status) call errquit('fock_2e: ma_pop?', 0)
c     
      call ga_sync()
c     
      do i=1,nfock
         call ga_dscal(vg_fock(i), 4.0d0) ! Undo scaling
         if (.not. oskel) call ga_symmetrize(vg_fock(i))
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
         call int2e_buf_write
c     insert munge of incomplete blocks here.
      end if
c     
      end
      subroutine fock_2e_rep(geom, basis, nfock, nbf, natom,
     $     jfac, kfac, tol2e, oskel, dens, fock, rdens,
     $     plist, ps, qdim, q4, ijkl)
      implicit none
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
      integer natom
      double precision jfac(nfock), kfac(nfock)
      double precision tol2e
      logical oskel
      double precision dens(nfock, nbf*nbf)
      double precision fock(nfock, nbf*nbf)
      double precision rdens(natom,natom)
      integer plist(2,*)        ! Pairs sorted by type
      double precision ps(*)    ! Schwarz screening for pair
      integer qdim
      double precision q4(qdim) ! Q4 factors for each quartet
      integer ijkl(qdim, 4)     ! Indices for quartet
c     
      integer natom_unique
      integer nshpair, nsh
      integer tdim
      parameter (tdim=1024)
      integer pinfo(2,tdim)     ! Start and end of each pair type
      integer ntype             ! No. of pair types
      integer max_p             ! Max pairs per block
c     
      integer nq                ! No. of quartets generated
c     
      integer ijklt, ijt, klt, ij, kl, ish, jsh, ksh, lsh
      integer  iat, jat, kat, lat, imat
      double precision q2, qq4, qjunk, sijkl, denmax, dijkl
      double precision dentol, integ_acc
      character*80 buff
c
      integer l_ijtype, k_ijtype
c
      integer ntreq
      double precision nqdone, nidone
      common/fred/nqdone, nidone, ntreq
c
      integer cnmax, cemax
      parameter (cnmax=1024, cemax=1024) ! Max no. of shells and centers
      integer cntoce(cnmax)     ! Map from contractions to centers
      integer cetobfr(2,cemax)  ! Map from centers to bf range
c
      integer next, current, nproc
c
      logical oprint
c     
      logical int2e_set_bf_range
      external int2e_set_bf_range
      integer nxtask
      external nxtask
c     
      if (oscfps) call pstat_on(ps_fock_2e)
c
      oprint = util_print('ao_replicated',print_never)
c
      call fock_init_rep_cmul(nbf) ! lookup table for f build
      if (owritefile) then
         if (.not. int2e_set_bf_range(
     $        1,nbf,1,nbf,1,nbf,1,nbf)) call errquit
     $        ('fock_2e_rep: set range failed',0)
      end if
      if (.not. geom_ncent_unique(geom, natom_unique)) call errquit
     $     ('fock_2e_rep: bad geom handle?',0)
      if ( .not. bas_numcont(basis, nsh) ) call errquit(
     $     'ao_fock_2e: problem with call to bas_numcont', basis)
      if (natom .gt. cemax) call errquit('fock_2e_rep: cemax',natom)
      if (nsh .gt. cnmax)   call errquit('fock_2e_rep: cnmax',nsh)
c     
c     Build maps (for speed) and construct density matrix over atoms
c     
      do ish = 1, nsh
         if (.not. bas_cn2ce(basis, ish, cntoce(ish)))
     $        call errquit('ao_rep: bad basis', 0)
      end do
      do iat = 1, natom
         if (.not. bas_ce2bfr(basis,iat,cetobfr(1,iat),cetobfr(2,iat)))
     $        call errquit('ao_rep: bad basis', 0)
      end do
      call dfill(natom*natom,0.0d0,rdens,1)
      do imat = 1, nfock
         call util_mat_reduce(nbf, natom, cetobfr, dens(1,imat), 
     $        rdens, 'abssum')
      end do
*      write(6,*) ' RDENS'
*      call output(rdens,1,natom,1,natom,natom,natom,1)
      denmax = 0.0d0
      do iat = 1, natom
         do jat = 1, iat
            denmax = max(denmax,rdens(jat,iat))
         end do
      end do
      if (owritefile) denmax = max(denmax,10.0d0)
c
      dentol = min(1d-3,tol2e/denmax) ! Threshold to screen integs only
c     
c     There is a BIG speed gain from agressive reduction of accuracy
c     in the actual computation of the integrals. However, the screening
c     on primitive integrals is done on radial prefactors not using
c     the Schwarz inequality, hence the extra factor of 0.1d0
c     (also tried 0.01 on C12H14 but was more expensive with no gain)
c     
      if (intacc .eq. 0.0d0) then
         integ_acc = min(1d-6,max(0.1d0*dentol,1d-15)) ! Variable
      else
         integ_acc = intacc     ! User controlled
      end if
      call int_acc_set(integ_acc)
*      write(6,111) dentol, integ_acc
* 111  format(1x,1p,' FRED  dentol=',d8.1,'  integ_acc=',d8.1)
c     
c     assume a dense integral list but with symmetry and figure out
c     no. of pairs per pair block dividing quartets between processors
c     and allowing for the max no. of quartets to be processed per block
c     
      nshpair = nsh*(nsh+1)/2
      max_p   = sqrt(
     $     dfloat(nshpair)**2*dfloat(natom_unique)) / 
     $     dfloat(natom)
      if (ga_nnodes() .gt. 1) then
         max_p = max_p/sqrt(10.0d0*ga_nnodes()) ! 10 tasks per processor
      end if
      max_p   = min(nint(sqrt(dfloat(qdim))),max_p)
      max_p   = max(10,max_p)
c     
c     generate info about the pairs
c     
      if (.not. ma_push_get(mt_int, nsh*(nsh+1)/2, 'fock2e_rep',
     $     l_ijtype, k_ijtype)) call errquit
     $     ('fock_2e_rep: ma failed', nsh*(nsh+1)/2)
 10   call fock_pairs(basis, max_p, .true., 1, natom, 1, natom, 
     $     dentol, plist, ps, pinfo, tdim, ntype, int_mb(k_ijtype))
c
c     Figure out no. of tasks to do on each request to the
c     shared counter.  Memory permitting we want this to be 
c     close to 1.
c
c     I worry about an infinite loop here ... is there always a 
c     solution to these constraints?  Should we limit the no. of passes?
c
      task_bf = max(ntype*ntype/(10*ga_nnodes()),1)
      if (max_p .lt. nint(sqrt(dfloat(qdim))) .and.
     $     task_bf .ge. 4) then
         max_p = min(nint(sqrt(dfloat(qdim))),2*max_p)
         goto 10
      end if
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*) ' max_p ', max_p, ' ntype ', ntype,
     $     ' task_bf ',task_bf
         call util_flush(6)
      end if
      if (.not. ma_pop_stack(l_ijtype))
     $     call errquit('fock_2e_rep: ma?',0)
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
      next = nxtask(nproc, task_bf)
      current = 0
      nq = 0
c
      nqdone = 0
      nidone = 0
      ntreq  = 1
c
      do ijklt = 2*ntype, 2, -1
         do ijt = min(ntype,ijklt-1), max(1,ijklt-ntype), -1
            if (next .eq. current) then
               klt = ijklt - ijt
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
c     schwarz screening
                     sijkl = ps(ij)*ps(kl)
                     if (sijkl .lt. dentol) goto 200 ! Next kl
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
c     symmetry test
                     if (oskel) then
                        if (.not. sym_atom_quartet(geom,iat,jat,kat,lat,
     $                       qjunk)) goto 200
                        qq4 = qq4*qjunk
                     end if
c     density screening
                     if (.not. owritefile) then
                        dijkl = max(
     $                       rdens(iat,jat),rdens(kat,lat),
     $                       rdens(iat,kat),rdens(iat,lat), 
     $                       rdens(jat,kat),rdens(jat,lat))
                        if (sijkl*dijkl .lt. tol2e) goto 200 ! Next kl
                     end if
c     kl diagonal
                     if (ksh.eq.lsh) qq4 = qq4*0.5d0
c     save it for later
                     if (nq .lt. qdim) then
                        nq = nq + 1
                        ijkl(nq,1) = ish
                        ijkl(nq,2) = jsh
                        ijkl(nq,3) = ksh
                        ijkl(nq,4) = lsh
                        q4(nq)     = qq4
                     else
                        call errquit('ao_replicated: qdim/max_p?',nq)
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
     $                 ijkl(1,4), q4)
                  nq = 0
               end if
               next = nxtask(nproc,task_bf)
               ntreq = ntreq + 1
            end if
            current = current + 1
         end do                 ! Next IJ type
      end do                    ! Next IJKL type
c     
      if (oscfps) call pstat_off(ps_fock_2e)
      next = nxtask(-nproc,1)
c
      if (oprint) then
         call begin_seq_output
         buff = ' '
         write(buff,99) ga_nodeid(), ntreq, nqdone, nidone
 99      format(i4,': tasks=',i4,1p,': nqdone=',d9.1,': neri=',d9.1)
         call write_seq(6,buff)
         call end_seq_output
      end if
c     
c     done
c     
      call int_acc_std()
c     
      end
      subroutine fock_rep_txs(basis, nfock, nbf,
     $     jfac, kfac, tol2e, dens, fock, 
     $     nq, iq, jq, kq, lq, q4)
      implicit none
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
c     
c     Given a list of integral quartets, notionally all of the
c     same type, compute the integrals and add them into the
c     fock matrix, optionally caching them on disk/memory.
c     
      integer lenscr, leneri, lenlab
      integer l_eri, k_eri, l_lab, k_lab, l_scr, k_scr
      integer k_i, k_j, k_k, k_l
c     
      integer neri
      double precision block_eff
      logical more
c
      integer maxg
c
      integer ntreq
      double precision nqdone, nidone
      common/fred/nqdone, nidone, ntreq
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
      if (maxscr .eq. 0) then
         call int_mem_2e4c(maxg, lenscr)
      else
         lenscr = maxscr
      end if
      leneri = maxeri
      lenlab = 4*maxeri
c
      if (util_print('fock_rep_txs', print_never)) then
         write(6,11) nq, lenscr, leneri
 11      format(' nq=',i6,', lenscr=',i8,', leneri=',i8)
         call util_flush(6)
      end if
c 
c     Could make this more intelligent so that if run out
c     of space request a smaller amount
c
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('fock_rep_txs: no room for eri',leneri)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('fock_rep_txs: no room for scr',lenscr)
      if (.not. ma_push_get(mt_int, lenlab, 'lab', l_lab, k_lab))
     $     call errquit('fock_rep_txs: no room for lab',lenlab)
c     
      k_i = k_lab
      k_j = k_i   + leneri
      k_k = k_j   + leneri
      k_l = k_k   + leneri
c     
      if (.not. intb_init4c(basis, iq, jq, basis, kq, lq,
     $     nq, lenscr, dbl_mb(k_scr), leneri, block_eff))
     $     call errquit('intb_init?',nq)
c     
 50   more = intb_2e4c(basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., tol2e, .false.,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), leneri, neri, lenscr, dbl_mb(k_scr))
c     
      if (neri .eq. 0) goto 100 ! Free memory before going home
c
      nidone = nidone + neri
c     
c     Cache integrals if desired
c     
      if (owritefile) then
         if (.not. int2e_file_write(neri,
     $        int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $        dbl_mb(k_eri))) call errquit('Arrgggghh',0)
      end if
c     
*      call print_integ_list(neri,
*     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
*     $     dbl_mb(k_eri))
c     
      call fock_2e_rep_label(nfock, nbf, jfac, kfac, tol2e, neri,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), dens, fock)
c
      if (more) then
         if (util_print('fock_rep_txs', print_never)) then
            write(6,*) ' Texas has split the request '
         end if
         goto 50
      end if
c     
 100  if (.not. ma_pop_stack(l_lab))
     $     call errquit('fock_rep_txs: cannot free lab',lenlab)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('fock_rep_txs: cannot free scr',lenlab)
      if (.not. ma_pop_stack(l_eri))
     $     call errquit('fock_rep_txs: cannot free eri',lenlab)
c     
      end
      subroutine fock_pairs(basis, max_p, oij, 
     $     iatlo, iathi, jatlo, jathi, 
     $     tol2e, ijlist, sij, ijinfo, max_ij, nt_ij, ijtype)
      implicit none
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
      integer ijtype(*)
c     
      integer junk
      integer nij, ij, ijt, i, j, jtop, prev, spch
      integer jlo, jhi, jnprim, jngen, jt, jtype
      integer ilo, ihi, inprim, ingen, it, itype
      double precision smax, test, tol2e
      integer n_in_b, old_nt
      integer iuat, juat        ! Unique atom numbers
c
      integer thetype
#include "itri.fh"
c     Below relies on #unique atoms<10, angmax<10, ngen<10.
      thetype(iuat,itype,inprim,ingen) = 
     $     iuat + 10*(itype+1 + 10*(ingen + 10*inprim))
c     
      if (.not. bas_ce2cnr(basis, iatlo, ilo, junk)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, iathi, junk, ihi)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, jatlo, jlo, junk)) call errquit
     $     ('fock_pairs: basis screwed ',0)
      if (.not. bas_ce2cnr(basis, jathi, junk, jhi)) call errquit
     $     ('fock_pairs: basis screwed ',0)
c     
      smax = schwarz_max()
c     
      nt_ij = 0
      nij = 0
      do i = ilo, ihi
         if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
     $        call errquit('fock_pairs: basis screwed ',0)
         if (.not. bas_cn2uce(basis,i,iuat))
     $        call errquit('fock_pairs: basis screwed ',0)
c     reorder L,s,p (-1,0,1) into order by complexity s,p,L
         if (itype.eq.-1) then
            itype = 1
         else if (itype.eq.0) then
            itype = -1
         else if (itype.eq.1) then
            itype = 0
         end if
         it = thetype(iuat,itype,inprim,ingen)
         jtop = jhi
         if (oij) jtop = i
         do j = jlo, jtop
            test = schwarz_shell(i,j)
            if (test*smax .gt. tol2e) then
               if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
     $              call errquit('fock_pairs: basis screwed ',0)
               if (.not. bas_cn2uce(basis,j,juat))
     $              call errquit('fock_pairs: basis screwed ',0)
c     reorder L,s,p (-1,0,1) into order by complexity s,p,L
               if (jtype.eq.-1) then
                  jtype = 1
               else if (jtype.eq.0) then
                  jtype = -1
               else if (jtype.eq.1) then
                  jtype = 0
               end if
               jt = thetype(juat,jtype,jnprim,jngen)
               nij           = nij + 1
               ijlist(1,nij) = i
               ijlist(2,nij) = j
               ijtype(nij)   = itri(it,jt)
            end if
         end do
      end do
      if (nij .eq. 0) return
c     
      call fock_pair_sort(nij,ijlist,ijtype)
c     
c     Split pair blocks on type boundaries with regard to size
c     (first do type only then combine)
c     
      nt_ij = 0
      prev = -99999999
      n_in_b = 0
      do ij = 1, nij
         if (ijtype(ij).ne.prev .or. n_in_b.ge.max_p) then
            if (nt_ij .ne. 0) ijinfo(2,nt_ij) = ij-1
            nt_ij = nt_ij + 1
            if (nt_ij .gt. max_ij) call errquit
     $           ('fock_pairs: dimension failure', nt_ij)
            ijinfo(1,nt_ij) = ij
            prev = ijtype(ij)
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
     $           ('fock_pairs: dimension failure(2)', nt_ij)
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
               if (.not. bas_continfo(basis,i,itype,inprim,ingen,spch))
     $              call errquit('fock_pairs: basis bad ',0)
               if (.not. bas_continfo(basis,j,jtype,jnprim,jngen,spch))
     $              call errquit('fock_pairs: basis bad ',0)
               if (.not. bas_cn2uce(basis,i,iuat))
     $              call errquit('fock_pairs: basis screwed ',0)
               if (.not. bas_cn2uce(basis,j,juat))
     $              call errquit('fock_pairs: basis screwed ',0)
               write(6,11) ij,i,j,itype,jtype,inprim,jnprim,
     $              ingen,jngen, iuat, juat
 11            format(1x,i6,2x,2i3,2x,2i3,2x,2i3,2x,2i3,2x,2i3)
            end do
         end do
      end if
c     
      end
      subroutine fock_pair_sort(n,ijlist,ijtype)
      implicit integer (a-z)
      dimension ijlist(2,n), ijtype(n)
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
      subroutine util_mat_reduce(n, nr, map, a, r, op)
      implicit none
c     
      integer n                 ! Original size [input]
      integer nr                ! Reduced size  [input]
      integer map(2,nr)         ! map(1,*)=lo, map(2,*)=hi [input]
      double precision a(n,n)   ! Original matrix [input]
      double precision r(nr,nr) ! Reduced matrix [output]
      character*(*) op          ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax' (extend as necessary)
c     
      integer ir, jr, i, j
      double precision sum
c     
*     write(6,*) ' util_mat_reduce: input matrix '
*     call output(a, 1, n, 1, n, n, n, 1)
c     
      do jr = 1, nr
         do ir = 1, nr
            sum = abs(r(ir,jr))
            if (op .eq. 'abssum') then
               do j = map(1,jr), map(2,jr)
                  do i = map(1,ir), map(2,ir)
                     sum = sum + abs(a(i,j))
                  end do
               end do
            else if (op .eq. 'absmax') then
               do j = map(1,jr), map(2,jr)
                  do i = map(1,ir), map(2,ir)
                     sum = max(sum , abs(a(i,j)))
                  end do
               end do
            else
               call errquit('util_mat_reduce: unknown op', 0)
            end if
            r(ir,jr) = sum
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
