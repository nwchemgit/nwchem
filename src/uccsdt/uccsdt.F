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
      integer alo, ahi, blo, bhi, clo, chi
      integer ablocksize, bblocksize, cblocksize
      integer spina, spinb, spinc
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
      call jan_full_transform(
     $     rtdb, basis, 
     $     nbf, nbf, nbf, nbf,
     $     nbf, nbf, nbf, nbf,
     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
     $     dbl_mb(k_moint), 'Dirac')
c
c     Crude loop over super blocks for debugging
c
      spina = 1
      spinb = 1
      spinc = 1
      ablocksize = 3
      bblocksize = 5
      cblocksize = 7
      do alo = nc(spina)+no(spina), nmo, ablocksize
         ahi = min(alo+ablocksize-1,nmo)
         do blo = nc(spinb)+no(spinb), nmo, bblocksize
            bhi = min(blo+bblocksize-1,nmo)
            do clo = nc(spinc)+no(spinc), nmo, cblocksize
               chi = min(clo+cblocksize-1,nmo)
c
               call uccsdt_triples_set_super_block
     $              (alo, ahi, blo, bhi, clo, chi, spina, spinb, spinc)
               call uccsdt_triples_moints(rtdb)
               call uccsdt_triples_dummy_amplitudes()
               call uccsdt_triples_unset_super_block()
c
            end do
         end do
      end do
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
     $     l_irs_b, k_irs_b, l_work, k_work
      integer i, ispin, lo, hi, ir
      character*4 uline(8)
c
      logical int_normalize
      external int_normalize
      logical movecs_read,movecs_read_header, status
      external movecs_read,movecs_read_header
      integer sym_number_ops
      external sym_number_ops
      data uline /8*'----'/
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
     $     so_uq_to_lab, so_lab1, so_lab2)
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
     $           int_mb(k_irs_a),map(1,ispin),
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
     $           int_mb(k_irs_a),map(1,ispin),
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
     $           int_mb(k_irs_a),map(1,ispin),
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
