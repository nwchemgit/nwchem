      logical function uccsdt(rtdb)
      implicit none
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
      integer rtdb
c
      character*80 title
      integer alo, ahi, blo, bhi, clo, chi, ctop
      integer ablocksize, bblocksize, cblocksize
      integer spina, spinb, spinc
      double precision dummy(5000)
      double precision uccsdt_triples_sub_blocked
      double precision energy(2,2), total
      logical oprint_super
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (ga_nodeid().eq.0 .and. 
     $     util_print('information', print_low)) then
         call util_print_centered(6,'NWChem UCCSD(T) Module',40,.true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
      end if
c
      call uccsdt_init(rtdb)
      oprint_super = util_print('superblock', print_debug)
c
c     For debug purposes generate the full list of Dirac integrals
c     in a common block
c
      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'most', l_most, k_most))
     $     call errquit('memory for most?',nbf*nbf)
      if (.not. ma_push_get(mt_dbl, nbf**4, 'moint', l_moint, k_moint))
     $     call errquit('memory for moint?',nbf**4)
      call util_transpose(dbl_mb(k_movecs(1)),nbf,
     $     dbl_mb(k_most),nbf,nbf,nbf)
      call jan_full_transform_noinit(
     $     rtdb, basis, 
     $     nbf, nbf, nbf, nbf,
     $     nbf, nbf, nbf, nbf,
     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
     $     dbl_mb(k_moint), 'Dirac')
c
c     Crude loop over super blocks for debugging
c
      ablocksize = 2
      bblocksize = 2
      cblocksize = 3
      total = 0.0d0
      do spina = 1, 2
        spinb = spina
        do spinc = 1, 2
          energy(spina,spinc) = 0.0d0
c         
          do alo = nc(spina)+no(spina)+1, nmo, ablocksize
            ahi = min(alo+ablocksize-1,nmo)
            do blo = nc(spinb)+no(spinb)+1, alo, bblocksize
              bhi = min(blo+bblocksize-1,nmo)
              ctop = nmo
              if (spina.eq.spinc) ctop = blo
              do clo = nc(spinc)+no(spinc)+1, ctop, cblocksize
                chi = min(clo+cblocksize-1,nmo)
c               
                call uccsdt_triples_set_super_block
     $               (alo, ahi, blo, bhi, clo, chi, spina, spinb, spinc,
     $               oprint_super)
                call uccsdt_triples_moints(rtdb)
                call uccsdt_triples_dummy_amplitudes()
                call dfill(5000, 0.0d0, dummy, 1)
                energy(spina,spinc) = energy(spina,spinc) + 
     $               uccsdt_triples_sub_blocked(dummy)
c               
                call uccsdt_triples_unset_super_block()
c               
                write(6,*) ' SUPERBLOCK: ', spina, spinc, alo, blo, clo,
     $               energy(spina,spinc)
c               
              end do
            end do
          end do
          write(6,*) ' ENERGY ', spina, spinc, energy(spina,spinc)
          total = total + energy(spina, spinc)
        end do
      end do
c
      write(6,77) total
 77   format(//' The total (T) energy is ', f20.12,//)
c
      if (.not. ma_pop_stack(l_moint)) call errquit('moint?',0)
      if (.not. ma_pop_stack(l_most)) call errquit('most?',0)
c
      call uccsdt_tidy(rtdb)
c
      uccsdt = .false.
c
      end
      subroutine uccsdt_init(rtdb)
      implicit none
#include "cuccsdtP.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
      integer rtdb
c
c     Load the common blocks /cuccsdt/ and /cbassym/
c     Initialize the integrals and schwarz
c
      character*255 movecs, title_vecs, basis_vecs
      character*8 scftypevecs, scftype, group
      integer nbf_vecs, nsets, nmos(2)
      integer num(0:7), first(0:7), last(0:7)
      integer nalpha, nbeta, nocc, nmixed
      integer l_occ_a, k_occ_a, k_occ_b, l_irs_a, k_irs_a,
     $     l_irs_b, k_irs_b, k_irs(2), l_work, k_work
      integer i, ispin, lo, hi, ir
      character*4 uline(8)
c
      logical int_normalize
      external int_normalize
      logical movecs_read,movecs_read_header, status, oprint_sym_info
      external movecs_read,movecs_read_header
      integer sym_number_ops
      external sym_number_ops
      data uline /8*'----'/
c
      call ifill(8*2,0,nc_sym,1)
      call ifill(8*2,0,no_sym,1)
      call ifill(8*2,0,nv_sym,1)
      call ifill(8*2, 0,c_sym(1,0,1),2)
      call ifill(8*2,-1,c_sym(2,0,1),2)
      call ifill(8*2, 0,o_sym(1,0,1),2)
      call ifill(8*2,-1,o_sym(2,0,1),2)
      call ifill(8*2, 0,v_sym(1,0,1),2)
      call ifill(8*2,-1,v_sym(2,0,1),2)
c     
c     load the geometry/basis set, initialize the integrals
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('uccsdt: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('uccsdt: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('uccsdt: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $     call errquit('uccsdt: no ao basis set', 0)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('uccsdt: bas_numbf?',0)
      if (.not. bas_numcont(basis, nsh))
     $     call errquit('uccsdt: bas_numcont?',0)
c
c     Print basis if requested before normalization
c
      nir = sym_number_ops(geom) + 1 ! For Abelian groups
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('uccsdt: geom_print', 0)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis summary', print_default)) then
            if (.not.bas_summary_print(basis))
     &           call errquit('uccsdt: basis summary print failed',911)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('uccsdt: bas_print', 0)
            write(6,*)
         endif
         if (nir.gt.1 .and. util_print('basis sym', print_default)) then
c     note that bf_per_ir is recomputed below
            call sym_bas_irreps(basis, .true., bf_per_ir)
         endif
         if (util_print('basis labels',print_high) .and.
     $        ga_nodeid().eq.0) then
            call bas_print_labels(basis)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(basis)) 
     $           call errquit('uccsdt:gmp', 0)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
      if (.not.int_normalize(rtdb, basis))
     $           call errquit('uccsdt:int_normalize failed', 0)
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get symmetry info about the AO basis set
c
      oprint_sym_info = util_print('sym_info', print_debug)
      call sym_abelian_bas_info(basis,
     $     sh_op_map, bf_op_map, bf_phase,
     $     sh_n_uq_op, bf_n_uq_op,
     $     sh_uq_op, bf_uq_op,
     $     sh_n_uq, bf_n_uq,
     $     sh_uq, bf_uq,
     $     char_tab, 
     $     bf_per_ir, bf_per_ir_cum,
     $     bf_so_ir,
     $     zir,
     $     so_uq_to_lab, so_lab1, so_lab2, oprint_sym_info)
c
c     Read the header on the MO vecs file
c
      if (.not. rtdb_cget(rtdb,'scf:input vectors',1,movecs)) 
     $     call errquit('uccsdt: no SCF MO vectors',0)
      if(.not.movecs_read_header(movecs,title_vecs,basis_vecs,
     $     scftypevecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('uccsdt: problem with movecs_read_header',1)
      if(nbf.ne.nbf_vecs)call errquit('mp2: MO vectors mismatch?',0)
      nmo = nmos(1)
c
c     Figure out frozen core options
c
      if (rtdb_get(rtdb,'uccsdt:frozen core:freeze by atoms',mt_log, 1, 
     $     status)) then
         if (.not. geom_num_core(geom,nc(1))) 
     $        call errquit('semi_direct_uccsd: geom_num_core?',0)
      else if (rtdb_get(rtdb, 'uccsdt:frozen core', mt_int, 1,
     $        nc(1))) then
      else
         nc(1) = 0
      endif
      nc(2) = nc(1)
c
c     From the no. of electrons etc, figure out the no. of occupied
c     and virtuals for each spin.
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('mp2:scftype absent',1)
      if(scftype.eq.'UHF')then
         if(.not.rtdb_get(rtdb,'scf:nalpha',mt_int,1,nalpha))
     $        call errquit('no alphas',1)
         if(.not.rtdb_get(rtdb,'scf:nbeta',mt_int,1,nbeta))
     $        call errquit('no betas',1)
      else if(scftype.eq.'RHF')then
         if(.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nocc))
     $        call errquit('no occ',1)
         nalpha=nocc
         nbeta=nocc
      else if(scftype.eq.'ROHF')then
         call errquit('uccsdt: ROHF not yet supported',0)
      else
         call errquit('uccsdt: unknown SCFTYPE ', 0)
      end if
      no(1) = nalpha - nc(1)
      no(2) = nbeta  - nc(2)
      nv(1) = nmo - nc(1) - no(1)
      nv(2) = nmo - nc(2) - no(2)
      nact  = nmo - nc(1)
c
c     Get the MO vectors, swap and then symmetry block them
c     
      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha vec',nbf,0,g_movecs(1)))
     $     call errquit('mp2_energy: ga_create failed',0)
      if(.not.ga_create(mt_dbl,nbf,nmo,'beta vec',nbf,0, g_movecs(2)))
     $     call errquit('mp2_energy: ga_create failed',0)
c
      if (.not. ma_push_get(mt_dbl,nmo*nbf,'alpha mo',
     $     l_movecs(1), k_movecs(1))) 
     $     call errquit('uccsdt:cannot allocate MA for vectors',nbf*nmo)
      if (.not. ma_push_get(mt_dbl,nmo*nbf,'beta mo',
     $     l_movecs(2), k_movecs(2))) 
     $     call errquit('uccsdt:cannot allocate MA for vectors',nbf*nmo)
c      
c     Temporary allocations for reading MOs and symmetry reordering
c
      if(.not.ma_push_get(mt_dbl,nbf*2,'occ alpha',
     $     l_occ_a,k_occ_a))      
     $     call errquit('uccsdt: ma problem with alpha occ',1)
      k_occ_b = k_occ_a + nbf   ! So that are contiguous for swap. Ugh.
      if(.not.ma_push_get(mt_int,nmo,'irs_a',l_irs_a,k_irs_a))
     $     call errquit('uccsdt: memory?',nmo)
      if(.not.ma_push_get(mt_int,nmo,'irs_b',l_irs_b,k_irs_b))
     $     call errquit('uccsdt: memory?',nmo)
c     
      k_irs(1) = k_irs_a ! for symmetry adaption indexed by spin
      k_irs(2) = k_irs_b
c
      if(.not.ma_push_get(mt_dbl,nmo,'work',l_work,k_work))
     $     call errquit('uccsdt: work?',nbf)
c     
      if(.not.movecs_read(movecs,1,dbl_mb(k_occ_a),
     $     eps(1,1),g_movecs(1)))
     $     call errquit('uccsdt:could not read alpha vecs',1)
      if(scftype.eq.'UHF')then
         if(.not.movecs_read(movecs,2,dbl_mb(k_occ_b),
     $        eps(1,2), g_movecs(2)))
     $        call errquit('could not read beta vecs',1)
      endif
c
      call movecs_swap(rtdb, 'uccsdt', scftype, g_movecs,
     $     dbl_mb(k_occ_a), eps(1,1))
c
      if (scftype .ne. 'UHF') then
         call ga_copy(g_movecs(1),g_movecs(2)) ! RHF just duplicate vectors
         call dcopy(nmo, eps(1,1), 1, eps(1,2), 1)
      end if
c     
      call sym_movecs_adapt(basis,1d-8,g_movecs(1),int_mb(k_irs_a),
     $     nmixed)
      call sym_movecs_adapt(basis,1d-8,g_movecs(2),int_mb(k_irs_b),
     $     nmixed)      
c
c     Change irreps to start counting at 0
c     
      do i=1,nmo
         int_mb(k_irs_a+i-1)=int_mb(k_irs_a+i-1)-1
         int_mb(k_irs_b+i-1)=int_mb(k_irs_b+i-1)-1
      end do
c
      do ispin = 1, 2
c
c     Core
c
         lo = 1
         hi = lo + nc(ispin) - 1
         if (nc(ispin) .gt. 0) then
            call moints_vecs_sym_sort(g_movecs(ispin),nbf,lo,hi,
     $           int_mb(k_irs(ispin)),map(1,ispin),
     $           dbl_mb(k_work), num, first, last)
            do ir = 0,7
               nc_sym(ir,ispin)  = num(ir)
               c_sym(1,ir,ispin) = first(ir)
               c_sym(2,ir,ispin) = last(ir)
            end do
         endif
c     
c     Occupied
c
         lo = hi + 1
         hi = lo + no(ispin) - 1
         if (no(ispin) .gt. 0) then
            call moints_vecs_sym_sort(g_movecs(ispin),nbf,lo,hi,
     $           int_mb(k_irs(ispin)),map(1,ispin),
     $           dbl_mb(k_work), num, first, last)
            do ir = 0,7
               no_sym(ir,ispin)  = num(ir)
               o_sym(1,ir,ispin) = first(ir)
               o_sym(2,ir,ispin) = last(ir)
            end do
         endif
c
c     Virtual
c
         lo = hi + 1
         hi = lo + nv(ispin) - 1
         if (hi .ne. nmo) call errquit('nmo is confused',0)
         if (nv(ispin) .gt. 0) then
            call moints_vecs_sym_sort(g_movecs(ispin),nbf,lo,hi,
     $           int_mb(k_irs(ispin)),map(1,ispin),
     $           dbl_mb(k_work), num, first, last)
            do ir = 0,7
               nv_sym(ir,ispin)  = num(ir)
               v_sym(1,ir,ispin) = first(ir)
               v_sym(2,ir,ispin) = last(ir)
            end do
         end if
c
c     Also reorder the eigenvalues
c
         call dcopy(nbf,eps(1,ispin),1,dbl_mb(k_work),1)
         do i=1,nmo
            eps(map(i,ispin),ispin) = dbl_mb(k_work+i-1)
         end do
      end do
c
c     Get rid of all the temporary workspace
c
      if (.not. ma_chop_stack(l_occ_a)) call errquit
     $     ('uccsdt_init: ma corrupted?',0)
c
c     Get local copied of the MO coefficients
c
      call ga_get(g_movecs(1), 1, nbf, 1, nmo, dbl_mb(k_movecs(1)), nbf)
      call ga_get(g_movecs(2), 1, nbf, 1, nmo, dbl_mb(k_movecs(2)), nbf)
c
c     Let the user know what's going on
c
      call sym_group_name(geom, group)
      if (ga_nodeid().eq.0 .and. 
     $     util_print('information', print_low)) then
         write(6,123) scftype, movecs(1:inp_strlen(movecs)), 
     $        nbf, nmo, nc(1), no(1), no(2), nv(1), nv(2), 
     $        group(1:inp_strlen(group)), nir, (zir(ir),ir=0,nir-1)
 123     format( 
     $        '  scftype            = ', a/
     $        '  MO vector file     = ', a/
     $        '  basis functions    = ', i5/
     $        '  molecular orbitals = ', i5/
     $        '  frozen core        = ', i5/
     $        '  alpha occupied     = ', i5/
     $        '  beta  occupied     = ', i5/
     $        '  alpha virtual      = ', i5/
     $        '  beta  virtual      = ', i5/
     $        '  point group        = ', a/
     $        '  irreducible reps   = ', i5/
     $        '  names of irreps    =  ', 8(a4,1x))
         write(6,1233) (uline(ir),ir=1,nir)
 1233    format('                       ',8(a4,1x:))
         write(6,124) '  no. of alpha core  = ',
     $        (nc_sym(ir,1),ir=0,nir-1)
         write(6,124) '  first  alpha core  = ',
     $        (c_sym(1,ir,1),ir=0,nir-1)
         write(6,124) '  last   alpha core  = ',
     $        (c_sym(2,ir,1),ir=0,nir-1)
c     
         write(6,124) '  no. of beta  core  = ',
     $        (nc_sym(ir,2),ir=0,nir-1)
         write(6,124) '  first  beta  core  = ',
     $        (c_sym(1,ir,2),ir=0,nir-1)
         write(6,124) '  last   beta  core  = ',
     $        (c_sym(2,ir,2),ir=0,nir-1)
c     
         write(6,124) '  no. of alpha occ   = ',
     $        (no_sym(ir,1),ir=0,nir-1)
         write(6,124) '  first  alpha occ   = ',
     $        (o_sym(1,ir,1),ir=0,nir-1)
         write(6,124) '  last   alpha occ   = ',
     $        (o_sym(2,ir,1),ir=0,nir-1)
c     
         write(6,124) '  no. of  beta occ   = ',
     $        (no_sym(ir,2),ir=0,nir-1)
         write(6,124) '  first   beta occ   = ',
     $        (o_sym(1,ir,2),ir=0,nir-1)
         write(6,124) '  last    beta occ   = ',
     $        (o_sym(2,ir,2),ir=0,nir-1)
c     
         write(6,124) '  no. of alpha virt  = ',
     $        (nv_sym(ir,1),ir=0,nir-1)
         write(6,124) '  first  alpha virt  = ',
     $        (v_sym(1,ir,1),ir=0,nir-1)
         write(6,124) '  last   alpha virt  = ',
     $        (v_sym(2,ir,1),ir=0,nir-1)
c     
         write(6,124) '  no. of  beta virt  = ',
     $        (nv_sym(ir,2),ir=0,nir-1)
         write(6,124) '  first   beta virt  = ',
     $        (v_sym(1,ir,2),ir=0,nir-1)
         write(6,124) '  last    beta virt  = ',
     $        (v_sym(2,ir,2),ir=0,nir-1)
 124     format(a,8(i4,1x))
      end if
      if (ga_nodeid().eq.0) then
         if (util_print('eigenvalues', print_high)) then
            write(6,*)
            write(6,*) '  Alpha eigenvalues (symmetry blocked)'
            write(6,*)
            call output(eps(1,1), 1, nmo, 1, 1, nmo, 1, 1)
            write(6,*)
            write(6,*) '  Beta  eigenvalues (symmetry blocked)'
            write(6,*)
            call output(eps(1,2), 1, nmo, 1, 1, nmo, 1, 1)
         end if
         if (util_print('eigenvectors', print_never)) then
            write(6,*)
            write(6,*) '  Alpha eigenvectors (symmetry blocked)'
            write(6,*)
            call output(dbl_mb(k_movecs(1)),1,nbf,1,nmo,nbf,nmo,1)
            write(6,*)
            write(6,*) '  Beta  eigenvectors (symmetry blocked)'
            write(6,*)
            call output(dbl_mb(k_movecs(2)),1,nbf,1,nmo,nbf,nmo,1)
         end if
      end if
c
      end
      subroutine uccsdt_tidy(rtdb)
      implicit none
#include "cuccsdtP.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb
c
      call schwarz_tidy()
      call int_terminate
      if (.not. ga_destroy(g_movecs(1))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',1)
      if (.not. ga_destroy(g_movecs(2))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',2)
c
*      call ma_summarize_allocated_blocks()
c
      if (.not. ma_pop_stack(l_movecs(2)))
     $     call errquit('uccsdt_tidy: ma_pop?',2)
      if (.not. ma_pop_stack(l_movecs(1)))
     $     call errquit('uccsdt_tidy: ma_pop?',1)
c
      if (.not. bas_destroy(basis))
     $     call errquit('uccsdt_tidy: basis?', 0)
      if (.not. geom_destroy(geom))
     $     call errquit('uccsdt_tidy: geom?', 0)
c
      call ma_summarize_allocated_blocks()
      call ga_summarize(0)
c
      end
      double precision function uccsdt_triples_sub_blocked(t1a)
      implicit none
#include "cuccsdtP.fh"
#include "util.fh"
      double precision t1a(*)
c
      integer ablk, alo, ahi, syma
      integer bblk, blo, bhi, symb
      integer cblk, clo, chi, symc
      integer iblk, ilo, ihi, symi
      integer jblk, jlo, jhi, symj
      integer kblk, klo, khi, symk
      integer mlo, mhi, msym
      integer elo, ehi, esym
      integer spine, spinm, test_sym
c
      double precision energy
c
      double precision int_abij(5000),
     &                 int_abjk(5000),
     &                 int_abik(5000),
     &                 int_bcij(5000),
     &                 int_acij(5000),
     &                 int_bcjk(5000),
     &                 int_acjk(5000),
     &                 int_bcik(5000),
     &                 int_acik(5000)
c
      double precision int_buf(5000), t_buf(5000), w(100000)
      integer btop, ctop, ktop, vlo, vhi, olo, ohi
      integer lentbuf, lenibuf
      logical oprint
c
      integer spin, spinck
      integer sym3, sym4, i, j, k, l, int_sym, t1_sym
      sym3(i,j,k) = ieor(ieor(i,j),k)
      sym4(i,j,k,l) = ieor(ieor(ieor(i,j),k),l)
c
      spin   = asuper_spin
      spinck = csuper_spin
c
      oprint = util_print('triples_sub_blocked',print_debug)
c
      energy = 0.0d0
      lentbuf = 5000
      lenibuf = 5000
c
      if (oprint) write(6,*)'TSB:spin, spinck: ', spin, spinck
      do ablk = 1, nablock
         alo = ablock(1,ablk)
         ahi = ablock(2,ablk)
         syma = ablock(3,ablk)
         if (oprint) write(6,*)'TSB:ablk, alo, ahi, syma ',
     &            ablk, alo, ahi, syma
         btop = nbblock
         if (asuper(1).eq.bsuper(1))btop = ablk
         do bblk = 1, btop
            blo = bblock(1,bblk)
            bhi = bblock(2,bblk)
            symb = bblock(3,bblk)
            if (oprint) write(6,*)'TSB:bblk, blo, bhi, symb ',
     &           bblk, blo, bhi, symb
            ctop = ncblock
            if (spin.eq.spinck .and. 
     $           bsuper(1).eq.csuper(1))ctop = bblk
            do cblk = 1, ctop
               clo = cblock(1,cblk)
               chi = cblock(2,cblk)
               symc = cblock(3,cblk)
               if (oprint) write(6,*)'TSB:cblk, clo, chi, symc ',
     &            cblk, clo, chi, symc
               do iblk = 1, noblock(spin)
                  ilo = oblock(1,iblk,spin)
                  ihi = oblock(2,iblk,spin)
                  symi = oblock(3,iblk,spin)
                  if (oprint) write(6,*)'TSB:iblk, ilo, ihi, symi ',
     &            iblk, ilo, ihi, symi
                  do jblk = 1, iblk
                     jlo = oblock(1,jblk,spin)
                     jhi = oblock(2,jblk,spin)
                     symj = oblock(3,jblk,spin)
                     if (oprint) write(6,*)'TSB:jblk, jlo, jhi, symj ',
     &                    jblk, jlo, jhi, symj
                     ktop = noblock(spinck)
                     if (spin .eq. spinck) ktop = jblk
                     do kblk = 1, ktop
                        klo = oblock(1,kblk,spin)
                        khi = oblock(2,kblk,spin)
                        symk = oblock(3,kblk,spin)
                        if (oprint) write(6,*)
     $                       'TSB:kblk, klo, khi, symk ',
     &                       kblk, klo, khi, symk
c
                        call dfill (100000, 0.0d0, w, 1)
c
                        spine = spin
                        esym = sym3(symc,symj,symi)
                        test_sym  = sym3(symk,syma,symb)
                        if ( spin.eq.spinck .and.
     $                       esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(14, 1, 
     &                         cblk, jblk, iblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(1, 1, 
     &                         kblk, ablk, bblk, int_buf, lenibuf)
                          call w_t2_ecji_int_ekab
     &                         (w, t_buf, int_buf,
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        endif
c
                        spine = spin
                        esym = sym3(symc,symj,symk)
                        test_sym  = sym3(symi,syma,symb)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(14, 1,
     &                         cblk, jblk, kblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(1, 1,
     &                         iblk, ablk, bblk, int_buf, lenibuf)
                          call w_t2_ecjk_int_eiab
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spin
                        esym = sym3(symc,symi,symk)
                        test_sym  = sym3(symj,syma,symb)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(14, 1,
     &                         cblk, iblk, kblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(1, 1,
     &                         jblk, ablk, bblk, int_buf, lenibuf)
                          call w_t2_ecik_int_ejab
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spin
                        esym = sym3(syma,symj,symi)
                        test_sym  = sym3(symk,symb,symc)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(15, 1,
     &                         ablk, jblk, iblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(3, 1,
     &                         kblk, bblk, cblk, int_buf, lenibuf)
                          call w_t2_eaji_int_ekbc
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spin
                        esym = sym3(symb,symj,symi)
                        test_sym  = sym3(symk,syma,symc)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(16, 1,
     &                         bblk, jblk, iblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(2, 1,
     &                         kblk, ablk, cblk, int_buf, lenibuf)
                          call w_t2_ebji_int_ekac
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spinck
                        esym = sym3(syma,symk,symj)
                        test_sym  = sym3(symi,symc,symb)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(17, 1,
     &                         ablk, kblk, jblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(5, 1,
     &                         iblk, cblk, bblk, int_buf, lenibuf)
                          call w_t2_eakj_int_eicb
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spinck
                        esym = sym3(symb,symk,symj)
                        test_sym  = sym3(symi,symc,syma)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(18, 1,
     &                         bblk, kblk, jblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(4, 1,
     &                         iblk, cblk, ablk, int_buf, lenibuf)
                          call w_t2_ebkj_int_eica
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spinck
                        esym = sym3(syma,symk,symi)
                        test_sym  = sym3(symj,symc,symb)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(17, 1,
     &                         ablk, kblk, iblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(5, 1,
     &                         jblk, cblk, bblk, int_buf, lenibuf)
                          call w_t2_eaki_int_ejcb
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spine = spinck
                        esym = sym3(symb,symk,symi)
                        test_sym  = sym3(symj,symc,syma)
                        if ( esym.eq.test_sym .and. 
     $                       nv_sym(esym,spine).gt.0 ) then
                          elo = v_sym(1,esym,spine)
                          ehi = v_sym(2,esym,spine)
                          call uccsdt_triples_get_list(18, 1,
     &                         bblk, kblk, iblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(4, 1,
     &                         jblk, cblk, ablk, int_buf, lenibuf)
                          call w_t2_ebki_int_ejca
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi, elo, ehi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi,
     &                         spin, spinck)
                        end if
c
                        spinm = spin
                        msym = sym3(symk,syma,symb)
                        test_sym  = sym3(symc,symj,symi)
                        if ( spin.eq.spinck .and.
     $                       msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(19, 1,
     &                         kblk, ablk, bblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(6, 1,
     &                         cblk, jblk, iblk, int_buf, lenibuf)
                          call w_t2_mkab_int_mcji
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spin
                        msym = sym3(symi,syma,symb)
                        test_sym  = sym3(symc,symj,symk)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(19, 1,
     &                         iblk, ablk, bblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(6, 1,
     &                         cblk, jblk, kblk, int_buf, lenibuf)
                          call w_t2_miab_int_mcjk
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spin
                        msym = sym3(symj,syma,symb)
                        test_sym  = sym3(symc,symi,symk)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(19, 1,
     &                         jblk, ablk, bblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(6, 1,
     &                         cblk, iblk, kblk, int_buf, lenibuf)
                          call w_t2_mjab_int_mcik
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spin
                        msym = sym3(symk,symb,symc)
                        test_sym  = sym3(syma,symj,symi)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(21, 1,
     &                         kblk, bblk, cblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(7, 1,
     &                         ablk, jblk, iblk, int_buf, lenibuf)
                          call w_t2_mkbc_int_maji
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spin
                        msym = sym3(symk,syma,symc)
                        test_sym  = sym3(symb,symj,symi)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(20, 1,
     &                         kblk, ablk, cblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(8, 1,
     &                         bblk, jblk, iblk, int_buf, lenibuf)
                          call w_t2_mkac_int_mbji
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spinck
                        msym = sym3(symi,symc,symb)
                        test_sym  = sym3(syma,symk,symj)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(23, 1,
     &                         iblk, cblk, bblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(9, 1,
     &                         ablk, kblk, jblk, int_buf, lenibuf)
                          call w_t2_micb_int_makj
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spinck
                        msym = sym3(symi,symc,syma)
                        test_sym  = sym3(symb,symk,symj)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(22, 1,
     &                         iblk, cblk, ablk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(10, 1,
     &                         bblk, kblk, jblk, int_buf, lenibuf)
                          call w_t2_mica_int_mbkj
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spinck
                        msym = sym3(symj,symc,symb)
                        test_sym  = sym3(syma,symk,symi)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(23, 1,
     &                         jblk, cblk, bblk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(9, 1,
     &                         ablk, kblk, iblk, int_buf, lenibuf)
                          call w_t2_mjcb_int_maki
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        spinm = spinck
                        msym = sym3(symj,symc,syma)
                        test_sym  = sym3(symb,symk,symi)
                        if ( msym.eq.test_sym .and. 
     $                       no_sym(msym,spinm).gt.0 ) then
                          mlo = o_sym(1,msym,spinm)
                          mhi = o_sym(2,msym,spinm)
                          call uccsdt_triples_get_list(22, 1,
     &                         jblk, cblk, ablk, t_buf, lentbuf)
                          call uccsdt_triples_get_list(10, 1,
     &                         bblk, kblk, iblk, int_buf, lenibuf)
                          call w_t2_mjca_int_mbki
     &                         (w, t_buf, int_buf, 
     &                         alo, ahi, blo, bhi, 
     &                         clo, chi,
     &                         ilo, ihi, jlo, jhi, 
     &                         klo, khi, mlo, mhi,
     &                         spin, spinck)
                        end if
c
                        int_sym = sym4(syma, symb, symi, symj)
                        t1_sym = ieor(symc, symk)
                        if ( int_sym.eq.0 .and. t1_sym.eq.0 ) then
                          call uccsdt_triples_get_list(11, 
     &                         ablk, bblk, iblk, jblk, int_abij,
     $                         lenibuf)
                        else
                          call dfill(lenibuf, 0.0d0, int_abij, 1)
                        end if
                        if (spin .eq. spinck) then
                          int_sym = sym4(syma, symb, symj, symk)
                          t1_sym = ieor(symc, symi)
                          if (int_sym.eq.0 .and. t1_sym.eq.0) then
                            call uccsdt_triples_get_list(11,
     &                           ablk, bblk, jblk, kblk, int_abjk,
     $                           lenibuf)
                          else
                            call dfill(lenibuf, 0d0, int_abjk, 1)
                          endif
                          int_sym = sym4(syma, symb, symi, symk)
                          t1_sym = ieor(symc, symj)
                          if (int_sym.eq.0 .and. t1_sym.eq.0) then
                            call uccsdt_triples_get_list(11,
     &                           ablk, bblk, iblk, kblk, int_abik,
     $                           lenibuf)
                          else
                            call dfill(lenibuf, 0d0, int_abik, 1)
                          endif
                          int_sym = sym4(symb, symc, symi, symj)
                          t1_sym = ieor(syma, symk)
                          if (int_sym.eq.0 .and. t1_sym.eq.0) then
                            call uccsdt_triples_get_list(13,
     &                           bblk, cblk, iblk, jblk, int_bcij,
     $                           lenibuf)
                          else
                            call dfill(lenibuf, 0d0, int_bcij, 1)
                          endif
                          int_sym = sym4(syma, symc, symi, symj)
                          t1_sym = ieor(symb, symk)
                          if (int_sym.eq.0 .and. t1_sym.eq.0) then
                            call uccsdt_triples_get_list(12,
     &                           ablk, cblk, iblk, jblk, int_acij,
     $                           lenibuf)
                          else
                            call dfill(lenibuf, 0d0, int_acij, 1)
                          endif
                        end if
                        int_sym = sym4(symb, symc, symj, symk)
                        t1_sym = ieor(syma, symi)
                        if (int_sym.eq.0 .and. t1_sym.eq.0) then
                          call uccsdt_triples_get_list(13,
     &                         bblk, cblk, jblk, kblk, int_bcjk, 
     $                         lenibuf)
                        else
                          call dfill(lenibuf, 0d0, int_bcjk, 1)
                        end if
c
                        int_sym = sym4(syma, symc, symj, symk)
                        t1_sym = ieor(symb, symi)
                        if (int_sym.eq.0 .and. t1_sym.eq.0) then
                          call uccsdt_triples_get_list(12,
     &                         ablk, cblk, jblk, kblk, int_acjk,
     $                         lenibuf)
                        else
                          call dfill(lenibuf, 0d0, int_acjk, 1)
                        endif
                        int_sym = sym4(symb, symc, symi, symk)
                        t1_sym = ieor(syma, symj)
                        if (int_sym.eq.0 .and. t1_sym.eq.0) then
                          call uccsdt_triples_get_list(13,
     &                         bblk, cblk, iblk, kblk, int_bcik,
     $                         lenibuf)
                        else
                          call dfill(lenibuf, 0d0, int_bcik, 1)
                        endif
                        int_sym = sym4(syma, symc, symi, symk)
                        t1_sym = ieor(symb, symj)
                        if (int_sym.eq.0 .and. t1_sym.eq.0) then
                          call uccsdt_triples_get_list(12,
     &                         ablk, cblk, iblk, kblk, int_acik,
     $                         lenibuf)
                        else
                          call dfill(lenibuf, 0d0, int_acik, 1)
                        endif
c
                        vlo = nc(spin) + no(spin) + 1
                        vhi = nmo
                        olo = nc(spin) + 1
                        ohi = nc(spin) + no(spin)
                        call uccsdt_w_t1a_int_energy
     &                    (w, t1a, eps, energy,
     &                     int_abij,  int_abjk, int_abik,  int_bcij,
     &                     int_acij,  int_bcjk, int_acjk,  int_bcik,
     &                     int_acik,
     &                     alo, ahi, blo, bhi, 
     &                     clo, chi,
     &                     ilo, ihi, jlo, jhi, 
     &                     klo, khi,
     &                     vlo, vhi, olo, ohi, spin, spinck, nw_max_nbf)
c
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      uccsdt_triples_sub_blocked = energy
c
      end
      subroutine uccsdt_w_t1a_int_energy
     &  (w, t1a, eps, energy,
     &   int_abij, int_abjk, int_abik, int_bcij,
     &   int_acij, int_bcjk, int_acjk, int_bcik,
     &   int_acik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   vlo, vhi, olo, ohi, spin, spinck, nmo)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        vlo, vhi, olo, ohi, spin, spinck, nmo
c
      double precision t1a(vlo:vhi,olo:ohi)
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
      double precision eps(nmo,2), energy, v, d
      integer i, j, k, a, b, c, ctop, ktop
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
           ctop = chi
           if (spin.eq.spinck) ctop = min(chi,b-1)
            do c = clo, ctop
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                    ktop = khi
                    if (spin.eq.spinck) ktop = min(khi,j-1)
                     do k = klo, ktop
c
                        v = w(k,j,i,c,b,a)
     &                        + t1a(c,k)*
     &                          int_abij(a,b,i,j)
     &                        + t1a(a,i)*
     &                          int_bcjk(b,c,j,k)
     &                        - t1a(b,i)*
     &                          int_acjk(a,c,j,k)
     &                        - t1a(a,j)*
     &                          int_bcik(b,c,i,k)
     &                        + t1a(b,j)*
     &                          int_acik(a,c,i,k)
c
                        if (spin .eq. spinck) then
                          v = v +
     &                        + t1a(c,i)*
     &                          int_abjk(a,b,j,k)
     &                        - t1a(c,j)*
     &                          int_abik(a,b,i,k)
     &                        + t1a(a,k)*
     &                          int_bcij(b,c,i,j)
     &                        - t1a(b,k)*
     &                          int_acij(a,c,i,j)
                        end if
c
                        d = eps(a,spin)+
     &                      eps(b,spin)+
     &                      eps(c,spinck)-
     &                      eps(i,spin)-
     &                      eps(j,spin)-
     &                      eps(k,spin)
c
                        energy = energy - 
     &                           v*w(k,j,i,c,b,a)/d
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      end
