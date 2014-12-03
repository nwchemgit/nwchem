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
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'uccsdt')
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
c     Initialize performance statistics
c
      call uccsdt_pstat_init(rtdb)
c
      call uccsdt_triples(rtdb)
c     
      call uccsdt_tidy(rtdb)
c
      call uccsdt_pstat_print
c
      call util_print_pop
      uccsdt = .true.
c
      end
      subroutine uccsdt_init(rtdb)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
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
c$$$c
c$$$c     For debug
c$$$c
c$$$      integer
c$$$     $     ishuqlo, ishuqhi, ioff(nw_max_nbf_small), 
c$$$     $     irange(2,0:7), 
c$$$     $     jshuqlo, jshuqhi, joff(nw_max_nbf_small,0:7),
c$$$     $     jrange(2,0:7),
c$$$     $     kshuqlo, kshuqhi, koff(nw_max_nbf_small,0:7), 
c$$$     $     krange(2,0:7),
c$$$     $     lshuqlo, lshuqhi, loff(nw_max_nbf_small,0:7), 
c$$$     $     lrange(2,0:7),
c$$$     $     lenij(0:7), lenkl(0:7), leng, symkloff(0:7)
c$$$      double precision g(100000)
c$$$c     
c$$$c     end debug
c$$$c

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
     $     call errquit('uccsdt: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('uccsdt: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('uccsdt: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $     call errquit('uccsdt: no ao basis set', 0, RTDB_ERR)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('uccsdt: bas_numbf?',0, BASIS_ERR)
      if (.not. bas_numcont(basis, nsh))
     $     call errquit('uccsdt: bas_numcont?',0, BASIS_ERR)
c
      if (nbf.gt.nw_max_nbf_small .or. nsh.gt.nw_max_shells_small)
     $     call errquit('uccsdt: nbf or nsh too large - modify'//
     $     ' *_small in nwc_const.fh',0, BASIS_ERR)
c
c     Other CCSD(T) options
c
      if (.not. rtdb_get(rtdb, 'uccsdt:ocheck', mt_log, 1, ocheck))
     $     ocheck = .false.
c
c     Print basis if requested before normalization
c
      nir = sym_number_ops(geom) + 1 ! For Abelian groups
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('uccsdt: geom_print', 0, GEOM_ERR)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis summary', print_default)) then
            if (.not.bas_summary_print(basis))
     $           call errquit('uccsdt: basis summary print failed',911,
     &       BASIS_ERR)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('uccsdt: bas_print', 0, BASIS_ERR)
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
     $           call errquit('uccsdt:gmp', 0, UNKNOWN_ERR)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
      if (.not.int_normalize(rtdb, basis))
     $           call errquit('uccsdt:int_normalize failed', 0, INT_ERR)
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
     $     sh_nat, bf_nat, sh_uq_bf,
     $     char_tab, 
     $     bf_per_ir, bf_per_ir_cum,
     $     bf_so_ir,
     $     zir,
     $     so_uq_to_lab, so_lab1, so_lab2, oprint_sym_info)
c
c     Read the header on the MO vecs file
c
      if (.not. rtdb_cget(rtdb,'scf:input vectors',1,movecs)) 
     $     call errquit('uccsdt: no SCF MO vectors',0, RTDB_ERR)
      if(.not.movecs_read_header(movecs,title_vecs,basis_vecs,
     $     scftypevecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('uccsdt: problem with movecs_read_header',1,
     &       DISK_ERR)
      if(nbf.ne.nbf_vecs)call errquit('mp2: MO vectors mismatch?',0,
     &       INPUT_ERR)
      nmo = nmos(1)
c
c     Figure out frozen core options
c
      if(.not. geom_freeze(rtdb, geom, 'uccsdt', nc(1))) then
         call errquit('semi_direct_uccsd: unable get get frozen core',0,
     &       GEOM_ERR)
      endif
 
      nc(2) = nc(1)
c
c     From the no. of electrons etc, figure out the no. of occupied
c     and virtuals for each spin.
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('mp2:scftype absent',1, RTDB_ERR)
      if(scftype.eq.'UHF')then
         if(.not.rtdb_get(rtdb,'scf:nalpha',mt_int,1,nalpha))
     $        call errquit('no alphas',1, RTDB_ERR)
         if(.not.rtdb_get(rtdb,'scf:nbeta',mt_int,1,nbeta))
     $        call errquit('no betas',1, RTDB_ERR)
      else if(scftype.eq.'RHF')then
         if(.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nocc))
     $        call errquit('no occ',1, RTDB_ERR)
         nalpha=nocc
         nbeta=nocc
      else if(scftype.eq.'ROHF')then
         call errquit('uccsdt: ROHF not yet supported',0, CAPMIS_ERR)
      else
         call errquit('uccsdt: unknown SCFTYPE ', 0, INPUT_ERR)
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
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nmo,'beta vec',nbf,0, g_movecs(2)))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
c
      if (.not. ma_push_get(mt_dbl,nmo*nbf,'alpha mo',
     $     l_movecs(1), k_movecs(1))) 
     $     call errquit('uccsdt:cannot allocate MA for vectors',nbf*nmo,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,nmo*nbf,'beta mo',
     $     l_movecs(2), k_movecs(2))) 
     $     call errquit('uccsdt:cannot allocate MA for vectors',nbf*nmo,
     &       MA_ERR)
c      
c     Temporary allocations for reading MOs and symmetry reordering
c
      if(.not.ma_push_get(mt_dbl,nbf*2,'occ alpha',
     $     l_occ_a,k_occ_a))      
     $     call errquit('uccsdt: ma problem with alpha occ',1, MA_ERR)
      k_occ_b = k_occ_a + nbf   ! So that are contiguous for swap. Ugh.
      if(.not.ma_push_get(mt_int,nbf,'irs_a',l_irs_a,k_irs_a))
     $     call errquit('uccsdt: memory?',nmo, MA_ERR)
      if(.not.ma_push_get(mt_int,nbf,'irs_b',l_irs_b,k_irs_b))
     $     call errquit('uccsdt: memory?',nmo, MA_ERR)
c     
      k_irs(1) = k_irs_a ! for symmetry adaption indexed by spin
      k_irs(2) = k_irs_b
c
      if(.not.ma_push_get(mt_dbl,nbf,'work',l_work,k_work))
     $     call errquit('uccsdt: work?',nbf, MA_ERR)
c     
      if(.not.movecs_read(movecs,1,dbl_mb(k_occ_a),
     $     eps(1,1),g_movecs(1)))
     $     call errquit('uccsdt:could not read alpha vecs',1, DISK_ERR)
      if(scftype.eq.'UHF')then
         if(.not.movecs_read(movecs,2,dbl_mb(k_occ_b),
     $        eps(1,2), g_movecs(2)))
     $        call errquit('could not read beta vecs',1, DISK_ERR)
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
         if (hi .ne. nmo) call errquit('nmo is confused',0, UNKNOWN_ERR)
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
     $     ('uccsdt_init: ma corrupted?',0, MA_ERR)
c
c     Get local copies of the MO coefficients
c
      call ga_get(g_movecs(1), 1, nbf, 1, nmo, dbl_mb(k_movecs(1)), nbf)
      call ga_get(g_movecs(2), 1, nbf, 1, nmo, dbl_mb(k_movecs(2)), nbf)
c
c     Allocate the hole/particle matrices
c
      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha hole',nbf,0,g_hole(1)))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nmo,'beta hole',nbf,0, g_hole(2)))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha part',nbf,0,g_part(1)))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nmo,'beta part',nbf,0, g_part(2)))
     $     call errquit('mp2_energy: ga_create failed',0, GA_ERR)

c
c     Let the user know what is going on
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
c     Initialize the amplitude access info and file handling
c
      if (.not. uccsdt_ampfile_init(rtdb)) call errquit
     $     ('uccsdt: failed initializing amplitude info',0, CALC_ERR)
c
c$$$c$$$      call soints_sovecs(nbf,nmo,dbl_mb(k_movecs(1)))
c$$$c$$$      stop 
c$$$c
c$$$c     Test the SO integral code
c$$$c
c$$$      do lshuqlo = 1, sh_n_uq, 4
c$$$         lshuqhi = min(sh_n_uq,lshuqlo+4-1)
c$$$         do kshuqlo = 1, lshuqlo, 4
c$$$            kshuqhi = min(sh_n_uq,kshuqlo+4-1)
c$$$            do jshuqlo = 1, sh_n_uq, 4
c$$$               jshuqhi = min(sh_n_uq,jshuqlo+4-1)
c$$$               do ishuqlo = 1, jshuqlo, 4
c$$$                  ishuqhi = min(sh_n_uq,ishuqlo+4-1)
c$$$
c$$$                  write(6,*) ishuqlo, jshuqlo, kshuqlo, lshuqlo
c$$$
c$$$                  call soints_block_info(basis, nbf,
c$$$     $                 ishuqlo, ishuqhi, ioff, irange, 
c$$$     $                 jshuqlo, jshuqhi, joff, jrange,
c$$$     $                 kshuqlo, kshuqhi, koff, krange,
c$$$     $                 lshuqlo, lshuqhi, loff, lrange,
c$$$     $                 lenij, lenkl, symkloff, leng)
c$$$                  if (leng .gt. 100000) call errquit('leng',0)
c$$$c
c$$$                  call soints_block(basis, nbf,
c$$$     $                 ishuqlo, ishuqhi, ioff, irange, 
c$$$     $                 jshuqlo, jshuqhi, joff, jrange,
c$$$     $                 kshuqlo, kshuqhi, koff, krange,
c$$$     $                 lshuqlo, lshuqhi, loff, lrange,
c$$$     $                 1d-12, g)
c$$$
c$$$c
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$      write(6,*) ' SOINTS seem OK'
c
      end
      subroutine uccsdt_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "amplitudes.fh"
      integer rtdb
c
      call schwarz_tidy()
      call int_terminate
c
      if (.not. ga_destroy(g_movecs(1))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',1, GA_ERR)
      if (.not. ga_destroy(g_movecs(2))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',2, GA_ERR)
      if (.not. ga_destroy(g_hole(1))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',1, GA_ERR)
      if (.not. ga_destroy(g_hole(2))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',2, GA_ERR)
      if (.not. ga_destroy(g_part(1))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',1, GA_ERR)
      if (.not. ga_destroy(g_part(2))) 
     $     call errquit('uccsdt_tidy: ga_destroy?',2, GA_ERR)
c
*      call ma_summarize_allocated_blocks()
c
      if (.not. ma_pop_stack(l_movecs(2)))
     $     call errquit('uccsdt_tidy: ma_pop?',2, MA_ERR)
      if (.not. ma_pop_stack(l_movecs(1)))
     $     call errquit('uccsdt_tidy: ma_pop?',1, MA_ERR)
c
      if (.not. bas_destroy(basis))
     $     call errquit('uccsdt_tidy: basis?', 0, BASIS_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('uccsdt_tidy: geom?', 0, GEOM_ERR)
c
      if (.not. uccsdt_ampfile_tidy()) 
     $     call errquit('uccsdt_tidy: ampfile?',0, UNKNOWN_ERR)
c
      call ma_summarize_allocated_blocks()
      call ga_summarize(0)
c
      end
      double precision function uccsdt_triples_sub_blocked()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
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
c      double precision int_abij(5000),
c     $     int_abjk(5000),
c     $     int_abik(5000),
c     $     int_bcij(5000),
c     $     int_acij(5000),
c     $     int_bcjk(5000),
c     $     int_acjk(5000),
c     $     int_bcik(5000),
c     $     int_acik(5000)
c     
      integer btop, ctop, ktop, vlo, vhi, olo, ohi,
     $        vlo_ck, vhi_ck, olo_ck, ohi_ck
      integer z_block_max
      integer lentbuf, l_tbuf, k_tbuf
      integer lenibuf, l_ibuf, k_ibuf
      integer lenwbuf, l_wbuf, k_wbuf
      integer lenrbuf, l_rbuf, k_rbuf
      integer lenxbuf, l_xbuf, k_xbuf
      integer len9ibuf, l_ibuf_abij, k_ibuf_abij,
     $                  l_ibuf_abjk, k_ibuf_abjk,
     $                  l_ibuf_abik, k_ibuf_abik,
     $                  l_ibuf_bcij, k_ibuf_bcij,
     $                  l_ibuf_acij, k_ibuf_acij,
     $                  l_ibuf_bcjk, k_ibuf_bcjk,
     $                  l_ibuf_acjk, k_ibuf_acjk,
     $                  l_ibuf_bcik, k_ibuf_bcik,
     $                  l_ibuf_acik, k_ibuf_acik
c  
      logical oprint
      integer nxtask, next, count, me, nproc
      external nxtask
c     
      integer spin, spinck, k_t1(2), k_t, k_t_ck
      integer sym3, sym4, i, j, k, l, int_sym, t1_sym
      sym3(i,j,k) = ieor(ieor(i,j),k)
      sym4(i,j,k,l) = ieor(ieor(ieor(i,j),k),l)
      spin   = asuper_spin
      spinck = csuper_spin
c     
      oprint = util_print('triples_sub_blocked',print_debug) 
     $     .and. ga_nodeid().eq.0
c     
      energy = 0.0d0
c
      lentbuf = max(e_block_max*vblock_max*oblock_max**2,
     &              m_block_max*oblock_max*vblock_max**2)
      if (.not.ma_push_get(mt_dbl,lentbuf,'tbuf',l_tbuf,k_tbuf))
     $   call errquit('uccsdt: tbuf?',lentbuf, MA_ERR)
c
      lenibuf = max(e_block_max*oblock_max*vblock_max**2,
     &              m_block_max*vblock_max*oblock_max**2)
      if (.not.ma_push_get(mt_dbl,lenibuf,'ibuf',l_ibuf,k_ibuf))
     $   call errquit('uccsdt: ibuf?',lenibuf, MA_ERR)
c
      lenwbuf = (vblock_max*oblock_max)**3
      if (.not.ma_push_get(mt_dbl,lenwbuf,'wbuf',l_wbuf,k_wbuf))
     $   call errquit('uccsdt: wbuf?',lenwbuf, MA_ERR)
c
      len9ibuf = (vblock_max*oblock_max)**2
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_abij',l_ibuf_abij,
     $   k_ibuf_abij))call errquit('uccsdt: ibuf_abij?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_abjk',l_ibuf_abjk,
     $   k_ibuf_abjk))call errquit('uccsdt: ibuf_abjk?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_abik',l_ibuf_abik,
     $   k_ibuf_abik))call errquit('uccsdt: ibuf_abik?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_bcij',l_ibuf_bcij,
     $   k_ibuf_bcij))call errquit('uccsdt: ibuf_bcij?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_acij',l_ibuf_acij,
     $   k_ibuf_acij))call errquit('uccsdt: ibuf_acij?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_bcjk',l_ibuf_bcjk,
     $   k_ibuf_bcjk))call errquit('uccsdt: ibuf_bcjk?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_acjk',l_ibuf_acjk,
     $   k_ibuf_acjk))call errquit('uccsdt: ibuf_acjk?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_bcik',l_ibuf_bcik,
     $   k_ibuf_bcik))call errquit('uccsdt: ibuf_bcik?',len9ibuf,
     &       MA_ERR)
      if (.not.ma_push_get(mt_dbl,len9ibuf,'ibuf_acik',l_ibuf_acik,
     $   k_ibuf_acik))call errquit('uccsdt: ibuf_acik?',len9ibuf,
     &       MA_ERR)
c
      z_block_max = max(e_block_max,m_block_max)
      lenrbuf = z_block_max*vblock_max**2
      if (.not.ma_push_get(mt_dbl,lenrbuf,'rbuf',l_rbuf,k_rbuf))
     $   call errquit('uccsdt: rbuf?',lenrbuf, MA_ERR)
      lenxbuf = oblock_max*vblock_max**3
      if (.not.ma_push_get(mt_dbl,lenxbuf,'xbuf',l_xbuf,k_xbuf))
     $   call errquit('uccsdt: xbuf?',lenxbuf, MA_ERR)
c
**      write(6,*)' e_block_max, m_block_max, oblock_max, vblock_max: ',
**     &            e_block_max, m_block_max, oblock_max, vblock_max

**      write(6,*)' lentbuf, lenibuf, lenwbuf, len9ibuf: ',
**     &            lentbuf, lenibuf, lenwbuf, len9ibuf
c
c note t1 is still dummy memory
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      next = nxtask(nproc, 1)
      count = 0
c     
      if (oprint) write(6,*)'TSB:spin, spinck: ', spin, spinck
      do ablk = nablock, 1, -1
         alo = ablock(1,ablk)
         ahi = ablock(2,ablk)
         syma = ablock(3,ablk)
         if (oprint) write(6,*)'TSB:ablk, alo, ahi, syma ',
     $        ablk, alo, ahi, syma
c     
c     A & B will be blocked in the same fashion if they are
c     the same super-blocks ... so can just restrict the
c     block indices
c     
         btop = nbblock
         if (asuper(1).eq.bsuper(1))btop = ablk
         do bblk = btop, 1, -1
            blo = bblock(1,bblk)
            bhi = bblock(2,bblk)
            symb = bblock(3,bblk)
            if (oprint) write(6,*)'TSB:bblk, blo, bhi, symb ',
     $           bblk, blo, bhi, symb
c     
c     The C blocking may not be the same as that for A & B, so
c     must look at the actual orbital indices to apply the
c     restriction for pure spin.
c     
            ctop = ncblock
            if (spin.eq.spinck) then
               do cblk = 1, ncblock
                  if (cblock(1,cblk) .ge. bhi) then 
                     ctop = cblk - 1
                     goto 10
                  end if
               end do
 10            continue
            endif
            do cblk = ctop, 1, -1
               clo = cblock(1,cblk)
               chi = cblock(2,cblk)
               symc = cblock(3,cblk)
               if (oprint) write(6,*)'TSB:cblk, clo, chi, symc ',
     $              cblk, clo, chi, symc
               do iblk = noblock(spin), 1, -1
                  ilo = oblock(1,iblk,spin)
                  ihi = oblock(2,iblk,spin)
                  symi = oblock(3,iblk,spin)
                  if (oprint) write(6,*)'TSB:iblk, ilo, ihi, symi ',
     $                 iblk, ilo, ihi, symi
                  do jblk = iblk, 1, -1
                     jlo = oblock(1,jblk,spin)
                     jhi = oblock(2,jblk,spin)
                     symj = oblock(3,jblk,spin)
                     if (oprint) write(6,*)'TSB:jblk, jlo, jhi, symj ',
     $                    jblk, jlo, jhi, symj
                     ktop = noblock(spinck)
                     if (spin .eq. spinck) ktop = jblk
                     do kblk = ktop, 1, -1
                        klo = oblock(1,kblk,spinck)
                        khi = oblock(2,kblk,spinck)
                        symk = oblock(3,kblk,spinck)
c     
                        test_sym = ieor(sym3(symi,symj,symk),
     $                       sym3(syma,symb,symc))
                        if (test_sym .eq. 0) then
                           if (next .eq. count) then
**                              write(6,*) ga_nodeid(), ' doing ',
**     $                             iblk, jblk, kblk, ablk, bblk, cblk
                              call util_flush(6)
                              if (oprint) write(6,*)
     $                             'TSB:kblk, klo, khi, symk ',
     $                             kblk, klo, khi, symk
c     
                              call dfill (lenwbuf, 0.0d0, 
     $                                    dbl_mb(k_wbuf),1)
c     
c     Now that we are testing on sym(ijkabc) above we do not
c     need to test on esym.eq.test_sym below
c
                              spine = spin
                              esym = sym3(symc,symj,symi)
c                              test_sym  = sym3(symk,syma,symb)
                              if ( spin.eq.spinck .and.
c     $                             esym.eq.test_sym .and. 
     $                             nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(14, 1, 
     $                                cblk, jblk, iblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(1, 1, 
     $                                kblk, ablk, bblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcji_zkab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              endif
c     
                              spine = spin
                              esym = sym3(symc,symj,symk)
c                              test_sym  = sym3(symi,syma,symb)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(14, 1,
     $                                cblk, jblk, kblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(1, 1,
     $                                iblk, ablk, bblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcjk_ziab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spine = spin
                              esym = sym3(symc,symi,symk)
c                              test_sym  = sym3(symj,syma,symb)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(14, 1,
     $                                cblk, iblk, kblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(1, 1,
     $                                jblk, ablk, bblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcik_zjab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spine = spin
                              esym = sym3(syma,symj,symi)
c                              test_sym  = sym3(symk,symb,symc)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(15, 1,
     $                                ablk, jblk, iblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(3, 1,
     $                                kblk, bblk, cblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zaji_zkbc
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spine = spin
                              esym = sym3(symb,symj,symi)
c                              test_sym  = sym3(symk,syma,symc)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(16, 1,
     $                                bblk, jblk, iblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(2, 1,
     $                                kblk, ablk, cblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbji_zkac
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spine = spinck
                              esym = sym3(syma,symk,symj)
c                              test_sym  = sym3(symi,symc,symb)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(17, 1,
     $                                ablk, kblk, jblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(5, 1,
     $                                iblk, cblk, bblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zakj_zicb
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spine = spinck
                              esym = sym3(symb,symk,symj)
c                              test_sym  = sym3(symi,symc,syma)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(18, 1,
     $                                bblk, kblk, jblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(4, 1,
     $                                iblk, cblk, ablk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbkj_zica
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spine = spinck
                              esym = sym3(syma,symk,symi)
c                              test_sym  = sym3(symj,symc,symb)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(17, 1,
     $                                ablk, kblk, iblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(5, 1,
     $                                jblk, cblk, bblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zaki_zjcb
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spine = spinck
                              esym = sym3(symb,symk,symi)
c                              test_sym  = sym3(symj,symc,syma)
c                              if ( esym.eq.test_sym .and. 
c     $                             nv_sym(esym,spine).gt.0 ) then
                              if ( nv_sym(esym,spine).gt.0 ) then
                                 elo = v_sym(1,esym,spine)
                                 ehi = v_sym(2,esym,spine)
                                 call uccsdt_triples_get_list(18, 1,
     $                                bblk, kblk, iblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(4, 1,
     $                                jblk, cblk, ablk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbki_zjca
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_tbuf), dbl_mb(k_ibuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                elo, ehi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spin
                              msym = sym3(symk,syma,symb)
c                              test_sym  = sym3(symc,symj,symi)
                              if ( spin.eq.spinck .and.
c     $                             msym.eq.test_sym .and. 
     $                             no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(19, 1,
     $                                kblk, ablk, bblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(6, 1,
     $                                cblk, jblk, iblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcji_zkab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spin
                              msym = sym3(symi,syma,symb)
c                              test_sym  = sym3(symc,symj,symk)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(19, 1,
     $                                iblk, ablk, bblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(6, 1,
     $                                cblk, jblk, kblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcjk_ziab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spinm = spin
                              msym = sym3(symj,syma,symb)
c                              test_sym  = sym3(symc,symi,symk)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(19, 1,
     $                                jblk, ablk, bblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(6, 1,
     $                                cblk, iblk, kblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zcik_zjab
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spin
                              msym = sym3(symk,symb,symc)
c                              test_sym  = sym3(syma,symj,symi)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(21, 1,
     $                                kblk, bblk, cblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(7, 1,
     $                                ablk, jblk, iblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zaji_zkbc
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spin
                              msym = sym3(symk,syma,symc)
c                              test_sym  = sym3(symb,symj,symi)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(20, 1,
     $                                kblk, ablk, cblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(8, 1,
     $                                bblk, jblk, iblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbji_zkac
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spinm = spinck
                              msym = sym3(symi,symc,symb)
c                              test_sym  = sym3(syma,symk,symj)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(23, 1,
     $                                iblk, cblk, bblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(9, 1,
     $                                ablk, kblk, jblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zakj_zicb
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              spinm = spinck
                              msym = sym3(symi,symc,syma)
c                              test_sym  = sym3(symb,symk,symj)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(22, 1,
     $                                iblk, cblk, ablk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(10, 1,
     $                                bblk, kblk, jblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbkj_zica
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spinck
                              msym = sym3(symj,symc,symb)
c                              test_sym  = sym3(syma,symk,symi)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(23, 1,
     $                                jblk, cblk, bblk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(9, 1,
     $                                ablk, kblk, iblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zaki_zjcb
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), -1.0d0)
                              end if
c     
                              spinm = spinck
                              msym = sym3(symj,symc,syma)
c                              test_sym  = sym3(symb,symk,symi)
c                              if ( msym.eq.test_sym .and. 
c     $                             no_sym(msym,spinm).gt.0 ) then
                              if ( no_sym(msym,spinm).gt.0 ) then
                                 mlo = o_sym(1,msym,spinm)
                                 mhi = o_sym(2,msym,spinm)
                                 call uccsdt_triples_get_list(22, 1,
     $                                jblk, cblk, ablk, dbl_mb(k_tbuf), 
     $                                lentbuf)
                                 call uccsdt_triples_get_list(10, 1,
     $                                bblk, kblk, iblk, dbl_mb(k_ibuf), 
     $                                lenibuf)
                                 call w_zbki_zjca
     $                                (dbl_mb(k_wbuf), 
     $                                dbl_mb(k_ibuf), dbl_mb(k_tbuf), 
     $                                alo, ahi, blo, bhi, clo, chi,
     $                                ilo, ihi, jlo, jhi, klo, khi, 
     $                                mlo, mhi, 
     $                                spin, spinck, dbl_mb(k_rbuf), 
     $                                dbl_mb(k_xbuf), 1.0d0)
                              end if
c     
                              int_sym = sym4(syma, symb, symi, symj)
                              t1_sym = ieor(symc, symk)
                              if ( int_sym.eq.0 .and. t1_sym.eq.0) then
                                 call uccsdt_triples_get_list(11, 
     $                                ablk, bblk, iblk, jblk, 
     $                                dbl_mb(k_ibuf_abij), len9ibuf)
                              else
                                 call dfill(len9ibuf, 0.0d0,
     $                                   dbl_mb(k_ibuf_abij), 1)
                              end if
                              if (spin .eq. spinck) then
                                 int_sym = sym4(syma, symb, symj, symk)
                                 t1_sym = ieor(symc, symi)
                                 if (int_sym.eq.0.and.t1_sym.eq.0) then
                                    call uccsdt_triples_get_list(11,
     $                                   ablk, bblk, jblk, kblk,
     $                                   dbl_mb(k_ibuf_abjk), len9ibuf)
                                 else
                                    call dfill(len9ibuf, 0d0,
     $                                   dbl_mb(k_ibuf_abjk), 1)
                                 endif
                                 int_sym = sym4(syma, symb, symi, symk)
                                 t1_sym = ieor(symc, symj)
                                 if (int_sym.eq.0.and.t1_sym.eq.0) then
                                    call uccsdt_triples_get_list(11,
     $                                   ablk, bblk, iblk, kblk,
     $                                   dbl_mb(k_ibuf_abik), len9ibuf)
                                 else
                                    call dfill(len9ibuf, 0d0,
     $                                   dbl_mb(k_ibuf_abik), 1)
                                 endif
                                 int_sym = sym4(symb, symc, symi, symj)
                                 t1_sym = ieor(syma, symk)
                                 if (int_sym.eq.0.and.t1_sym.eq.0) then
                                    call uccsdt_triples_get_list(13,
     $                                   bblk, cblk, iblk, jblk,
     $                                   dbl_mb(k_ibuf_bcij), len9ibuf)
                                 else
                                    call dfill(len9ibuf, 0d0,
     $                                   dbl_mb(k_ibuf_bcij), 1)
                                 endif
                                 int_sym = sym4(syma, symc, symi, symj)
                                 t1_sym = ieor(symb, symk)
                                 if (int_sym.eq.0.and.t1_sym.eq.0) then
                                    call uccsdt_triples_get_list(12,
     $                                   ablk, cblk, iblk, jblk,
     $                                   dbl_mb(k_ibuf_acij), len9ibuf)
                                 else
                                    call dfill(len9ibuf, 0d0,
     $                                   dbl_mb(k_ibuf_acij), 1)
                                 endif
                              end if
                              int_sym = sym4(symb, symc, symj, symk)
                              t1_sym = ieor(syma, symi)
                              if (int_sym.eq.0 .and. t1_sym.eq.0) then
                                 call uccsdt_triples_get_list(13,
     $                                bblk, cblk, jblk, kblk,
     $                                dbl_mb(k_ibuf_bcjk), len9ibuf)
                              else
                                 call dfill(len9ibuf, 0d0, 
     $                                dbl_mb(k_ibuf_bcjk), 1)
                              end if
c     
                              int_sym = sym4(syma, symc, symj, symk)
                              t1_sym = ieor(symb, symi)
                              if (int_sym.eq.0 .and. t1_sym.eq.0) then
                                 call uccsdt_triples_get_list(12,
     $                                ablk, cblk, jblk, kblk, 
     $                                dbl_mb(k_ibuf_acjk), len9ibuf)
                              else
                                 call dfill(len9ibuf, 0d0, 
     $                                dbl_mb(k_ibuf_acjk), 1)
                              endif
                              int_sym = sym4(symb, symc, symi, symk)
                              t1_sym = ieor(syma, symj)
                              if (int_sym.eq.0 .and. t1_sym.eq.0) then
                                 call uccsdt_triples_get_list(13,
     $                                bblk, cblk, iblk, kblk, 
     $                                dbl_mb(k_ibuf_bcik), len9ibuf)
                              else
                                 call dfill(len9ibuf, 0d0, 
     $                                dbl_mb(k_ibuf_bcik), 1)
                              endif
                              int_sym = sym4(syma, symc, symi, symk)
                              t1_sym = ieor(symb, symj)
                              if (int_sym.eq.0 .and. t1_sym.eq.0) then
                                 call uccsdt_triples_get_list(12,
     $                                ablk, cblk, iblk, kblk, 
     $                                dbl_mb(k_ibuf_acik), len9ibuf)
                              else
                                 call dfill(len9ibuf, 0d0, 
     $                                dbl_mb(k_ibuf_acik), 1)
                              endif
c
                              k_t1(1) = k_t1a
                              k_t1(2) = k_t1b
c                              write(6,*)' k_t1: ',k_t1
c
                              vlo = nc(spin) + no(spin) + 1
                              vhi = nmo
                              olo = nc(spin) + 1
                              ohi = nc(spin) + no(spin)
                              k_t = k_t1(spin)
c                              write(6,*)' for spin: ',spin
c                              write(6,*)' vlo, vhi, olo, ohi, k_t: ',
c     &                                    vlo, vhi, olo, ohi, k_t
c
                              vlo_ck = nc(spinck) + no(spinck) + 1
                              vhi_ck = nmo
                              olo_ck = nc(spinck) + 1
                              ohi_ck = nc(spinck) + no(spinck)
                              k_t_ck = k_t1(spinck)
c                              write(6,*)' for spinck: ',spinck
c              write(6,*)' vlo_ck, vhi_ck, olo_ck, ohi_ck, k_t_ck: ',
c     &                    vlo_ck, vhi_ck, olo_ck, ohi_ck, k_t_ck
c
                              call w_energy
     $                             (dbl_mb(k_wbuf), dbl_mb(k_t),
     $                             dbl_mb(k_t_ck), eps, energy,
     $                             dbl_mb(k_ibuf_abij),
     $                             dbl_mb(k_ibuf_abjk),
     $                             dbl_mb(k_ibuf_abik),
     $                             dbl_mb(k_ibuf_bcij),
     $                             dbl_mb(k_ibuf_acij),
     $                             dbl_mb(k_ibuf_bcjk),
     $                             dbl_mb(k_ibuf_acjk),
     $                             dbl_mb(k_ibuf_bcik),
     $                             dbl_mb(k_ibuf_acik),
     $                             alo, ahi, blo, bhi, 
     $                             clo, chi,
     $                             ilo, ihi, jlo, jhi, 
     $                             klo, khi,
     $                             vlo, vhi, olo, ohi, 
     $                             vlo_ck, vhi_ck, olo_ck, ohi_ck, 
     $                             spin, spinck, 
     $                             nw_max_nbf_small)
c     
                              next = nxtask(nproc, 1)
                           endif
                           count = count + 1
                        end if
c     
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      if (.not. ma_chop_stack(l_tbuf)) call errquit('ma tbuf chop?',0,
     &       MA_ERR)
c
      next = nxtask(-nproc, 1)
      call ga_sync()
c     
      call ga_dgop(1, energy, 1, '+')
      call ga_sync()
c     
      uccsdt_triples_sub_blocked = energy
      end
      subroutine uccsdt_triples_super_blocking(spina,spinc,
     $     max_nsuper_block, nasuper, ncsuper, a_super_block,
     $     c_super_block)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
      integer spina, spinc      ! [input]
      integer max_nsuper_block  ! [input]
      integer nasuper, ncsuper  ! [output] 
      integer a_super_block(2,max_nsuper_block) ! [output]
      integer c_super_block(2,max_nsuper_block) ! [output]
c
c     From the available memory, calculation parameters, 
c     machine characteristics, determine the various
c     blocksizes.
c
c     The triples algorithm proceeds in the following 
c     steps.
c
c     loop over spin cases
c     .   loop over super-blocks of a, b, c
c     .       allocate space for lists
c     .       compute and store integrals, load amplitudes
c     .       loop over sub-blocks of a, b, c, i, k, 
c     .           make W, V, form energy
c     .       free lists
c
c     The super-blocking is so that we can fit the necessary
c     integral+amplitude lists & workspace in memory.
c     
c     The routine uccsdt_triples_sub_blocking() given the
c     super-block ranges determines the sub-blocking and
c     allocates space for the lists.  It needs oblock_max
c     and vblock_max to control the sub-blocking.  These are
c     chosen small enough to provide sufficient parallelism
c     but large enough to be efficient according to the performance
c     model, but not so large as to blow out memory.
c
c     Loop through plausible choices for the superblocking
c     .  Loop thru the superblocks
c     .     Compute the maximum amount of memory usage
c     .     Count the no. of times the integrals are evaluated
c
c     Choose the superblocking with lowest no. of integral
c     evaluations that fits in memory.
c
      logical oprint, oprint_more
      integer maxsh, wblockdim, ublockdim
      integer a1, a2, c1, c2, numasuper, numcsuper, 
     $     a_block_max, c_block_max
      integer numasuper_best, numcsuper_best, npass_best
      double precision avail_ma, avail_ga, ga_max, ma_max,
     $     ma_moints, ga_moints, ma_list, ga_list, ma_wv,
     $     ga_moints_max, ma_list_max, 
     $     ga_list_max, ma_wv_max, ma_test
      integer leneri, lenscr, maxq, npass, ablk, bblk, cblk
      integer adim, bdim, cdim, ctop, ir
      logical ok
c
      oprint = util_print('information', print_low) .and. 
     $     (ga_nodeid().eq.0)
      oprint_more = util_print('super-blocking', print_debug) .and. 
     $     (ga_nodeid().eq.0)
      oprint = oprint .or. oprint_more
c
c     Determine amount of memory available
c
      call util_ma_ga_avail(avail_ma, avail_ga)
      call ga_dgop(11,avail_ma,1,'min')
      call ga_dgop(12,avail_ga,1,'+') ! Sum available GA across all nodes
      if (oprint) then
         write(6,2) spina, spinc, ga_uses_ma(), ga_memory_limited(), 
     $        avail_ga, avail_ma
 2       format(/' Super-blocking information for spin-case ', 2i2,/,
     $        '  GA uses MA       ', l5,4x,4x,
     $        '  GA limited       ', l5,4x,4x,1p,/,
     $        '  GA available     ', d9.2,4x,
     $        '  MA available     ', d9.2)
         call util_flush(6)
      end if
c
      avail_ma = avail_ma * 0.95d0
      avail_ga = avail_ga * 0.95d0
c
c     MAXSH is the largest shell dimension
c
      if (.not. bas_nbf_cn_max(basis, maxsh)) call errquit
     $     ('uccsdt: bas_nbf_cn_max?', 0, BASIS_ERR)
c
c     WBLOCKDIM is used to limit the amount of memory needed
c     to store AO integrals ... the bigger it is the more
c     efficiently the integrals will run but the minimum is
c     MAXSH.  Give it room for at least 10 of these. The final
c     value is determined by the 4-index code which will use more 
c     if possible.
c
      wblockdim = min(10*maxsh, nbf)
c
c     E_BLOCK_MAX is just the largest no. of virtuals in any irrep
c
            e_block_max = 0
            do ir = 0, nir-1
               e_block_max = max(e_block_max,nv_sym(ir,spina),
     $              nv_sym(ir,spinc))
            end do
c
c     M_BLOCK_MAX is just the largest no. of occupied in any irrep
c
            m_block_max = 0
            do ir = 0, nir-1
               m_block_max = max(m_block_max,no_sym(ir,spina),
     $              no_sym(ir,spinc))
            end do
c
c     local memory required by first part of 4-index does
c     not depend upon the blocking and the subsequent 4-index
c     steps are always smaller.
c
      maxq = wblockdim*wblockdim
      call intb_mem_2e4c(leneri, lenscr)
      leneri = max(leneri,100*(6**4))
      ma_moints = 
     $     maxsh*maxsh*nbf*(no(spina)+no(spinc)) + 
     $     maxsh*maxsh*wblockdim*wblockdim + 
     $     maxq + leneri + lenscr + 
     $     ma_sizeof(mt_int, 4*maxq+4*leneri, mt_dbl)
c
c     a & b must be identically blocked
c     
      a1 = nc(spina)+no(spina)+1
      a2 = nmo
      c1 = nc(spinc)+no(spinc)+1
      c2 = nmo
c
      numasuper_best = -1
      numcsuper_best = -1
      npass_best = 99999999
      do numasuper = 1, max_nsuper_block
         a_block_max = (nv(spina)-1)/numasuper + 1
         call uccsdt_block_range(a1, a2, a_block_max, nasuper,
     $        a_super_block, 2, max_nsuper_block)
         do numcsuper = 1, max_nsuper_block
            c_block_max = (nv(spinc)-1)/numcsuper + 1
            call uccsdt_block_range(c1, c2, c_block_max, ncsuper,
     $           c_super_block, 2, max_nsuper_block)
c
c     UBLOCKDIM is used to limit the amount of GA used by
c     the 4-index code.  Smaller means less memory, but 
c     also means less parallelism and more communication.
c     The minimum is MAXSH.  As a compromise set this to 
c     the same as the A virtual blocksize.  The actual size
c     will be determined by the 4-index code and more will 
c     be used if possible.
c
            ublockdim = max(maxsh, a_block_max)
c
            npass = 0
            ga_max = 0d0
            ma_max = 0d0
            do ablk = 1, nasuper
               adim = a_super_block(2,ablk)-a_super_block(1,ablk)+1
               do bblk = 1, ablk
                  bdim = a_super_block(2,bblk)-a_super_block(1,bblk)+1
c
c     Since C is blocked differently from A & B we must compare
c     the actual orbital indices to apply the restriction
c
                  ctop = ncsuper
                  if (spina .eq. spinc) then 
                     do cblk = 1, ncsuper
                        if ( c_super_block(1,cblk).ge.
     $                       a_super_block(2,bblk) ) then
                           ctop = cblk - 1
                           goto 10
                        end if
                     end do
 10                  continue
                  end if
                  do cblk = 1, ctop
                     cdim = c_super_block(2,cblk)-
     $                    c_super_block(1,cblk)+1
                     npass = npass + 1
c
c     Compute the memory used for list pointers & values
c
                     call uccsdt_triples_sub_blocking(
     $                    a_super_block(1,ablk), a_super_block(2,ablk),
     $                    a_super_block(1,bblk), a_super_block(2,bblk),
     $                    c_super_block(1,cblk), c_super_block(2,cblk),
     $                    spina, spinc,
     $                    .false., .true., ma_list, ga_list)
c
c     MO integral transform needs some scratch GAs
c
                     ga_moints = ublockdim*nbf*
     $                    ((adim+bdim)*no(spina) + cdim*no(spinc))
c
c     WV just uses some MA
c
                     ma_wv = 9*(oblock_max*vblock_max)**2 + 
     $                    (oblock_max*vblock_max)**3 + 
     $                    e_block_max*oblock_max*vblock_max**2 + 
     $                    e_block_max*vblock_max*oblock_max**2
c
                     if (ga_list+ga_moints .gt. ga_max) then
                        ga_max = ga_list+ga_moints
                        ga_list_max = ga_list
                        ga_moints_max = ga_moints
                     end if
                     ma_test = max(ma_max, ma_list +
     $                    max(ma_wv, ma_moints))
                     if (ma_test .gt. ma_max) then
                        ma_max = ma_test
                        ma_list_max = ma_list
                        ma_wv_max = ma_wv
                     end if
c
                  end do
               end do
            end do
c
c     Does it fit in memory?
c     If it fits & uses the fewest integral passes then record it
c
            call util_ma_ga_test(ma_max, ga_max/ga_nnodes(), ok)
****c
****c temp force multiple super blocks
****c
****            if (numasuper.eq.1.or.numcsuper.eq.1)ok=.false.            
c
            if (ok .and. npass.lt.npass_best) then
               npass_best = npass
               numasuper_best = numasuper
               numcsuper_best = numcsuper
c     
               if (oprint_more) then
                  write(6,1) nasuper, ncsuper, 
     $                 oblock_max, vblock_max, npass,
     $                 ga_max, ga_list_max, ga_moints_max, 
     $                 ma_max, ma_list_max, ma_moints, ma_wv_max
 1                format(/
     $                 '  AC super blocks   ', 2i5,/,
     $                 '  OV max sub-blocks ', 2i5,/,
     $                 '  integral passes  ', i5,/,
     $                 1p,
     $                 '  GA max used      ', d9.2,/,
     $                 '      lists        ', d9.2,/,
     $                 '      transform    ', d9.2,/,
     $                 '  MA max used      ', d9.2,/,
     $                 '      lists        ', d9.2,/,
     $                 '      transform    ', d9.2,/,
     $                 '      w & v        ', d9.2,//)
               end if
            end if
         end do
      end do
c
      if (numasuper_best .eq. -1) 
     $     call errquit('uccsdt: failed to find a good blocking', 0,
     &       UNKNOWN_ERR)
c
c     Restore the best choice
c
      numasuper = numasuper_best
      numcsuper = numcsuper_best
      a_block_max = (nv(spina)-1)/numasuper + 1
      call uccsdt_block_range(a1, a2, a_block_max, nasuper,
     $     a_super_block, 2, max_nsuper_block)
      c_block_max = (nv(spinc)-1)/numcsuper + 1
      call uccsdt_block_range(c1, c2, c_block_max, ncsuper,
     $     c_super_block, 2, max_nsuper_block)
c
      if (oprint) then
         write(6,3) nasuper, ncsuper, npass_best
 3       format(
     $        '  No. of AC blocks ', 2i5,3x,
     $        '  Integral passes  ', i5/)
         call util_flush(6)
      end if
      if (oprint_more) then
         write(6,*) ' A super-blocking'
         write(6,4) (a_super_block(1,ablk),a_super_block(2,ablk),
     $        ablk=1,nasuper)
         write(6,*) ' C super-blocking'
         write(6,4) (c_super_block(1,cblk),c_super_block(2,cblk),
     $        cblk=1,ncsuper)
 4       format(3x,5(2i5,5x))
         call util_flush(6)
      end if
c     
      end
      subroutine uccsdt_triples_sub_blocking
     $     (alo, ahi, blo, bhi, clo, chi, spina, spinc,
     $     oprint, odryrun, ma_used, ga_used)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer alo, ahi, blo, bhi, clo, chi ! [input]
      integer spina, spinc
      logical oprint
      logical odryrun
      double precision ma_used, ga_used
c
c     The outer-3-most loops of the triples are over blocks of virtual 
c     indices.  The blocking is determined outside of this routine
c     so as to minimize recomputation of the integrals within memory 
c     constraints.
c
c     This routine, given the specified super blocks and spins of 
c     abc, sets up all of the blocking info in the common block
c
c     1. The sub-blocking info for both virtual and occupied indices
c     2. The addressing info for indexing into the 5 global
c     .  arrays by sub-block indices
c     3. Allocation of the global arrays themselves
c
c     if (odryrun) then
c     .  do not allocate memory but return the maximum amount
c     .  of MA and GA allocated in this call tree (in doubles)
c     endif
c
      integer ir, a1, a2, b1, b2, c1, c2, spin, list, i
      integer nblk              ! No. of blocks in a symmetry
      integer spini, spinj, spink, spinb
      integer asub, bsub, csub, osub, a, b, c, o
c
      ma_used = 0.0d0
      ga_used = 0.0d0
c
      asuper(1) = alo
      asuper(2) = ahi
      bsuper(1) = blo
      bsuper(2) = bhi
      csuper(1) = clo
      csuper(2) = chi
      asuper_spin = spina
      csuper_spin = spinc
      spinb     = spina
c
c     Determine the optimal sub-blocking for the occupied
c     and virtual indices ... these should be determined by the
c     machine characteristics and the performance model, 
c     adequate parallelism, ...
c
c     Must also ensure that there are not so many sub-blocks as
c     to overflow the dimensions of the adressing arrays.
c
c     For debug set them to something sensible and then
c     recompute their actual values.
c
c      oblock_max = 20
c
c tmp set oblock for testing multiple blocks.
c
      oblock_max = 4
c
c tmp set vblock for testing multiple blocks.
c
      vblock_max = 10
c
c     Sub-block the 'a' super-block ... sub-blocks will be restricted
c     to orbitals of the same symmetry whereas super-blocks can contain
c     multiple symmetries and split symmetry blocks.
c
      nablock = 0
      do ir = 0, nir-1
         a1 = max(alo,v_sym(1,ir,spina))
         a2 = min(ahi,v_sym(2,ir,spina))
         call uccsdt_block_range(a1, a2, vblock_max, nblk, 
     $        ablock(1,nablock+1), 3, nw_max_nbf_small-nablock)
         do i = nablock+1,nablock+nblk
            ablock(3,i) = ir
         end do
         nablock = nablock + nblk
      end do
      do asub = 1, nablock
         do a = ablock(1,asub),ablock(2,asub)
            ablock_inv(a) = asub
         end do
      end do
c
c     Ditto for b
c
      nbblock = 0
      do ir = 0, nir-1
         b1 = max(blo,v_sym(1,ir,spinb))
         b2 = min(bhi,v_sym(2,ir,spinb))
         call uccsdt_block_range(b1, b2, vblock_max, nblk, 
     $        bblock(1,nbblock+1), 3, nw_max_nbf_small-nbblock)
         do i = nbblock+1,nbblock+nblk
            bblock(3,i) = ir
         end do
         nbblock = nbblock + nblk
      end do
      do bsub = 1, nbblock
         do b = bblock(1,bsub),bblock(2,bsub)
            bblock_inv(b) = bsub
         end do
      end do
c
c     Ditto for c
c
      ncblock = 0
      do ir = 0, nir-1
         c1 = max(clo,v_sym(1,ir,spinc))
         c2 = min(chi,v_sym(2,ir,spinc))
         call uccsdt_block_range(c1, c2, vblock_max, nblk, 
     $        cblock(1,ncblock+1), 3, nw_max_nbf_small-ncblock)
         do i = ncblock+1,ncblock+nblk
            cblock(3,i) = ir
         end do
         ncblock = ncblock + nblk
      end do
      do csub = 1, ncblock
         do c = cblock(1,csub),cblock(2,csub)
            cblock_inv(c) = csub
         end do
      end do
c
c     The blocking of the occupied never changes but just re-do it 
c     every time for convenience
c
      do spin = 1, 2
         noblock(spin) = 0
         do ir = 0, nir-1
            call uccsdt_block_range(
     $           o_sym(1,ir,spin), o_sym(2,ir,spin), 
     $           oblock_max, nblk, 
     $           oblock(1,noblock(spin)+1,spin), 3, 
     $           nw_max_nbf_small-noblock(spin))
            do i = noblock(spin)+1,noblock(spin)+nblk
               oblock(3,i,spin) = ir
            end do
            noblock(spin) = noblock(spin) + nblk
         end do
         do osub = 1, noblock(spin)
            do o = oblock(1,osub,spin),oblock(2,osub,spin)
               oblock_inv(o,spin) = osub
            end do
         end do
      end do
c
c     Recompute oblock_max and vblock_max
c
      vblock_max = 0
      oblock_max = 0
      do i = 1, nablock
         vblock_max = max(vblock_max, ablock(2,i)-ablock(1,i)+1)
      end do
      do i = 1, nbblock
         vblock_max = max(vblock_max, bblock(2,i)-bblock(1,i)+1)
      end do
      do i = 1, ncblock
         vblock_max = max(vblock_max, cblock(2,i)-cblock(1,i)+1)
      end do
      do spin = 1, 2
         do i = 1, noblock(spin)
            oblock_max = max(oblock_max,
     $           oblock(2,i,spin)-oblock(1,i,spin)+1)
         end do
      end do
c
      if (oprint) then
         write(6,*) ' A sub-blocking ', alo, ahi, spina, nablock
         write(6,1) (i,ablock(1,i),ablock(2,i),ablock(3,i),i=1,nablock)
         write(6,*) ' B sub-blocking ', blo, bhi, spinb, nbblock
         write(6,1) (i,bblock(1,i),bblock(2,i),bblock(3,i),i=1,nbblock)
         write(6,*) ' C sub-blocking ', clo, chi, spinc, ncblock
         write(6,1) (i,cblock(1,i),cblock(2,i),cblock(3,i),i=1,ncblock)
         write(6,*) ' Alpha occupied blocking'
         write(6,1) (i,oblock(1,i,1),oblock(2,i,1),oblock(3,i,1),
     $        i=1,noblock(1))
         write(6,*) ' Beta  occupied blocking'
         write(6,1) (i,oblock(1,i,2),oblock(2,i,2),oblock(3,i,2),
     $        i=1,noblock(2))
 1       format(4x,i5,4x,3i5)
      endif
c
c     Now set up the addressing arrays for the by-sub-block
c     addressing of the 1-d global arrays with contiguous sub-blocks
c
c     Mixed spin case has
c     .   spin(a)=spin(b)=spin(i)=spin(j)
c     .   spin(c)=spin(k)
c
c     Pure spin case has all spins the same but the same lists
c     are created as for mixed spin, except that all integrals
c     are antisymmetrized.
c
c     3-x integrals
c     1.  <ei||ab>  spin(e)=spin(a)  pure
c     2.  <ek|ac>   spin(e)=spin(a)  mixed 
c     3.  <ek|bc>   spin(e)=spin(a)  mixed
c     4.  <ei|ca>   spin(e)=spin(c)  mixed
c     5.  <ei|cb>   spin(e)=spin(c)  mixed
c
c     1-x integrals
c     6.  <mc|ik>  spin(m)=spin(a)  mixed
c     7.  <ma||ij> spin(m)=spin(a)  pure
c     8.  <mb||ij> spin(m)=spin(a)  pure
c     9.  <ma|ki>  spin(m)=spin(c)  mixed
c     10. <mb|ki>  spin(m)=spin(c)  mixed
c
c     2-x integrals (for V)
c     11. <ab||ij> pure
c     12. <ac|ik>  mixed
c     13. <bc|ik>  mixed
c
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c     15.  t(e,a,i,j) spin(e)=spin(a) pure
c     16.  t(e,b,i,j) spin(e)=spin(a) pure
c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
c
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
c     19.  t(m,i,a,b) spin(m)=spin(a) pure
c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
c
      spini = spina
      spinj = spina
      spink = spinc
c
c     Setup the names of the lists for easier printing & debugging
c
c     3-x integrals
      listname(1) = '<ei||ab>'
      listname(2) = '<ek|ac>'
      listname(3) = '<ek|bc>'
      listname(4) = '<ei|ca>'
      listname(5) = '<ei|cb>'
c     1-x integrals
      listname(6) = '<mc|ik>'
      listname(7) = '<ma||ij>'
      listname(8) = '<mb||ij>'
      listname(9) = '<ma|ki>'
      listname(10) = '<mb|ki>'
c     2-x integrals (for V)
      listname(11) = '<ab||ij>'
      listname(12) = '<ac|ik>'
      listname(13) = '<bc|ik>'
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
      listname(14) = 't(e,c,i,k)'
      listname(15) = 't(e,a,i,j)'
      listname(16) = 't(e,b,i,j)'
      listname(17) = 't(e,a,k,i)'
      listname(18) = 't(e,b,k,i)'
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
      listname(19) = 't(m,i,a,b)'
      listname(20) = 't(m,k,a,c)'
      listname(21) = 't(m,k,b,c)'
      listname(22) = 't(m,i,c,a)'
      listname(23) = 't(m,i,c,b)'
c
c     Compute the size of each dimension of the pointer arrays
c     and then computer their size so they can be MA'd
c     
c     3-x integrals
      listinfo(6,1) = 1                         ! 1.  <ei||ab> 
      listinfo(7,1) = noblock(spini)
      listinfo(8,1) = nablock
      listinfo(9,1) = nbblock                   
      listinfo(6,2) = 1                         ! 2.  <ek|ac>  
      listinfo(7,2) = noblock(spink)
      listinfo(8,2) = nablock
      listinfo(9,2) = ncblock                   
      listinfo(6,3) = 1                         ! 3.  <ek|bc>  
      listinfo(7,3) = noblock(spink)
      listinfo(8,3) = nbblock
      listinfo(9,3) = ncblock                   
      listinfo(6,4) = 1                         ! 4.  <ei|ca>  
      listinfo(7,4) = noblock(spini)
      listinfo(8,4) = ncblock
      listinfo(9,4) = nablock                   
      listinfo(6,5) = 1                         ! 5.  <ei|cb>  
      listinfo(7,5) = noblock(spini)
      listinfo(8,5) = ncblock
      listinfo(9,5) = nbblock                   
c     1-x integrals
      listinfo(6,6) = 1                         ! 6.  <mc|ik>
      listinfo(7,6) = ncblock
      listinfo(8,6) = noblock(spini)
      listinfo(9,6) = noblock(spink)            
      listinfo(6,7) = 1                         ! 7.  <ma||ij>
      listinfo(7,7) = nablock
      listinfo(8,7) = noblock(spini)
      listinfo(9,7) = noblock(spinj)            
      listinfo(6,8) = 1                         ! 8.  <mb||ij>
      listinfo(7,8) = nbblock
      listinfo(8,8) = noblock(spini)
      listinfo(9,8) = noblock(spinj)            
      listinfo(6,9) = 1                         ! 9.  <ma|ki>
      listinfo(7,9) = nablock
      listinfo(8,9) = noblock(spink)
      listinfo(9,9) = noblock(spini)            
      listinfo(6,10)= 1                         ! 10. <mb|ki> 
      listinfo(7,10)= nbblock
      listinfo(8,10)= noblock(spink)
      listinfo(9,10)= noblock(spini)            
c     2-x integrals (for V)
      listinfo(6,11)= nablock                   ! 11. <ab||ij>
      listinfo(7,11)= nbblock
      listinfo(8,11)= noblock(spini)
      listinfo(9,11)= noblock(spinj)            
      listinfo(6,12)= nablock                   ! 12. <ac|ik>
      listinfo(7,12)= ncblock
      listinfo(8,12)= noblock(spini)
      listinfo(9,12)= noblock(spink)            
      listinfo(6,13)= nbblock                   ! 13. <bc|ik>
      listinfo(7,13)= ncblock
      listinfo(8,13)= noblock(spini)
      listinfo(9,13)= noblock(spink)            
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
      listinfo(6,14)= 1                         ! 14.  t(e,c,i,k)
      listinfo(7,14)= ncblock
      listinfo(8,14)= noblock(spini)
      listinfo(9,14)= noblock(spink)            
      listinfo(6,15)= 1                         ! 15.  t(e,a,i,j)
      listinfo(7,15)= nablock
      listinfo(8,15)= noblock(spini)
      listinfo(9,15)= noblock(spinj)            
      listinfo(6,16)= 1                         ! 16.  t(e,b,i,j)
      listinfo(7,16)= nbblock
      listinfo(8,16)= noblock(spini)
      listinfo(9,16)= noblock(spinj)            
      listinfo(6,17)= 1                         ! 17.  t(e,a,k,i)
      listinfo(7,17)= nablock
      listinfo(8,17)= noblock(spink)
      listinfo(9,17)= noblock(spini)            
      listinfo(6,18)= 1                         ! 18.  t(e,b,k,i)
      listinfo(7,18)= nbblock
      listinfo(8,18)= noblock(spink)
      listinfo(9,18)= noblock(spini)            
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
      listinfo(6,19)= 1                         ! 19.  t(m,i,a,b)
      listinfo(7,19)= noblock(spini)
      listinfo(8,19)= nablock
      listinfo(9,19)= nbblock                   
      listinfo(6,20)= 1                         ! 20.  t(m,k,a,c)
      listinfo(7,20)= noblock(spink)
      listinfo(8,20)= nablock
      listinfo(9,20)= ncblock                   
      listinfo(6,21)= 1                         ! 21.  t(m,k,b,c)
      listinfo(7,21)= noblock(spink)
      listinfo(8,21)= nbblock
      listinfo(9,21)= ncblock                   
      listinfo(6,22)= 1                         ! 22.  t(m,i,c,a)
      listinfo(7,22)= noblock(spini)
      listinfo(8,22)= ncblock
      listinfo(9,22)= nablock                   
      listinfo(6,23)= 1                         ! 23.  t(m,i,c,b)
      listinfo(7,23)= noblock(spini)
      listinfo(8,23)= ncblock
      listinfo(9,23)= nbblock                   
c
c     Allocate and then create the pointer/reclen arrays
c
      do list = 1, 23
         listinfo(1,list) = listinfo(6,list)*listinfo(7,list)*
     $        listinfo(8,list)*listinfo(9,list)
         ma_used = ma_used + 
     $        ma_sizeof(mt_int, listinfo(1,list), mt_dbl)
         if (.not. odryrun) then
            if (.not. ma_push_get(mt_int,2*listinfo(1,list),
     $           listname(list), listinfo(3,list), listinfo(2,list))) 
     $           call errquit
     $           ('uccsdt_triples_set_super_block: ma of list?', list,
     &       MA_ERR)
         end if
      end do
c
c     1.  <ei||ab>  spin(e)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spini), oblock(1,1,spini), nablock, ablock, 
     $     nbblock, bblock, int_mb(listinfo(2,1)), listinfo(4,1),
     $     odryrun)
      ga_used = ga_used + listinfo(4,1)
c
c     2.  <ek|ac>   spin(e)=spin(a)  mixed 
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nablock, ablock, 
     $     ncblock, cblock, int_mb(listinfo(2,2)), listinfo(4,2),
     $     odryrun)
      ga_used = ga_used + listinfo(4,2)
c
c     3.  <ek|bc>   spin(e)=spin(a)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nbblock, bblock, 
     $     ncblock, cblock, int_mb(listinfo(2,3)), listinfo(4,3),
     $     odryrun)
      ga_used = ga_used + listinfo(4,3)
c      
c     4.  <ei|ca>   spin(e)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock, 
     $     nablock, ablock, int_mb(listinfo(2,4)), listinfo(4,4),
     $     odryrun)
      ga_used = ga_used + listinfo(4,4)
c
c     5.  <ei|cb>   spin(e)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock, 
     $     nbblock, bblock, int_mb(listinfo(2,5)), listinfo(4,5),
     $     odryrun)
      ga_used = ga_used + listinfo(4,5)
c
c     1-x integrals
c
c     6.  <mc|ik>  spin(m)=spin(a)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     ncblock, cblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,6)), 
     $     listinfo(4,6), odryrun)
      ga_used = ga_used + listinfo(4,6)
c
c     7.  <ma||ij> spin(m)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     nablock, ablock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,7)), 
     $     listinfo(4,7), odryrun)
      ga_used = ga_used + listinfo(4,7)
c
c     8.  <mb||ij> spin(m)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     nbblock, bblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,8)), 
     $     listinfo(4,8), odryrun)
      ga_used = ga_used + listinfo(4,8)
c
c     9.  <ma|ki>  spin(m)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     nablock, ablock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,9)), 
     $     listinfo(4,9), odryrun)
      ga_used = ga_used + listinfo(4,9)
c
c     10. <mb|ki>  spin(m)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     nbblock, bblock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,10)), 
     $     listinfo(4,10), odryrun)
      ga_used = ga_used + listinfo(4,10)
c
c     2-x integrals (for V)
c
c     11. <ab||ij> pure
c
      call uccsdt_2x_pointers(nablock, ablock, nbblock, bblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,11)), 
     $     listinfo(4,11), odryrun)
      ga_used = ga_used + listinfo(4,11)
c
c     12. <ac|ik>  mixed
c
      call uccsdt_2x_pointers(nablock, ablock, ncblock, cblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,12)), 
     $     listinfo(4,12), odryrun)
      ga_used = ga_used + listinfo(4,12)
c
c     13. <bc|ik>  mixed
c
      call uccsdt_2x_pointers(nbblock, bblock, ncblock, cblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,13)), 
     $     listinfo(4,13), odryrun)
      ga_used = ga_used + listinfo(4,13)
c
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
c
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     ncblock, cblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,14)), 
     $     listinfo(4,14), odryrun)
      ga_used = ga_used + listinfo(4,14)
c
c     15.  t(e,a,i,j) spin(e)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     nablock, ablock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,15)), 
     $     listinfo(4,15), odryrun)
      ga_used = ga_used + listinfo(4,15)
c      
c     16.  t(e,b,i,j) spin(e)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     nbblock, bblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,16)), 
     $     listinfo(4,16), odryrun)
      ga_used = ga_used + listinfo(4,16)
c
c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     nablock, ablock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,17)), 
     $     listinfo(4,17), odryrun)
      ga_used = ga_used + listinfo(4,17)
c
c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     nbblock, bblock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,18)), 
     $     listinfo(4,18), odryrun)
      ga_used = ga_used + listinfo(4,18)
c
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
c
c     19.  t(m,i,a,b) spin(m)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spini), oblock(1,1,spini), nablock, ablock,
     $     nbblock, bblock, int_mb(listinfo(2,19)), listinfo(4,19),
     $     odryrun)
      ga_used = ga_used + listinfo(4,19)
c
c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nablock, ablock,
     $     ncblock, cblock, int_mb(listinfo(2,20)), listinfo(4,20),
     $     odryrun)
      ga_used = ga_used + listinfo(4,20)
c
c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nbblock, bblock,
     $     ncblock, cblock, int_mb(listinfo(2,21)), listinfo(4,21),
     $     odryrun)
      ga_used = ga_used + listinfo(4,21)
c
c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock,
     $     nablock, ablock, int_mb(listinfo(2,22)), listinfo(4,22),
     $     odryrun)
      ga_used = ga_used + listinfo(4,22)
c
c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock,
     $     nbblock, bblock, int_mb(listinfo(2,23)), listinfo(4,23),
     $     odryrun)
      ga_used = ga_used + listinfo(4,23)
c
c     Print a summary of the lists
c
      if (oprint) then
         write(6,2)
 2       format('  Summary of block-sorted list of integrals',
     $        ' and amplitudes',//,
     $        '  list    name    #blocks  #elements   block',
     $        ' dimensions',/
     $        ' ----- ---------- ------- ----------',
     $        ' --------------------')
         do list = 1, 23
            write(6,3) list, listname(list), listinfo(1,list), 
     $           listinfo(4,list), 
     $           listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $           listinfo(9,list)
 3          format(1x,i5,1x,a10,i8,i11,1x,4i5)
         end do
      endif
c     
c     Finally, we can allocate the arrays for the lists themselves
c
      if (.not. odryrun) then
         do list = 1, 23
            if (listinfo(4,list) .gt. 0) then
               if (.not. ga_create(mt_dbl, listinfo(4,list), 1, 
     $              listname(list), 0, 0, listinfo(5,list))) 
     $              call errquit
     $              ('uccsdt_triples_set_super_block: ga_create', list,
     &       GA_ERR)
               if (list .lt. 14) then
                 call ga_zero(listinfo(5,list))
               else
                 call ga_fill(listinfo(5,list), 1d300)
               endif
            else
               listinfo(5,list) = 99999999 ! invalid handle for empty array
            endif
         end do
      end if
c     
      end
      subroutine uccsdt_1x_3x_pointers(
     $     ne_sym,
     $     niblock, iblock,
     $     nablock, ablock, 
     $     nbblock, bblock,
     $     pointer, length, 
     $     odryrun)
      implicit none
      integer ne_sym(0:7)       ! [input] No. of e orbitals by symmetry
      integer niblock, iblock(3,*) ! [input] Blocking for i
      integer nablock, ablock(3,*) ! [input] Blocking for a
      integer nbblock, bblock(3,*) ! [input] Blocking for b
      integer pointer(2,niblock, nablock, nbblock) ! [output]
      integer length            ! [output] Total length of the array
c
c     Generate the pointer to the start of a block of <ei|ab>
c     for all e (of appropriate symmetry) and i, a, b blocked
c     as specified
c
c     pointer(iblock,ablock,bblock) = start of block of <ei|ab>
c
c     Eventually the lists must be stored with randomization to avoid
c     hot spts.  Also, eventually the lists must store double precision
c     pointers to use the long GA interface for > 16 GB arrays.
c
c     If (odryrun) then
c     .   do not actually make the pointer array ... just compute
c     .   the size of the GA
c     endif
c
      integer ptr, reclen
      integer a, b, i
      integer syma, symb, symi, syme
      integer dima, dimb, dimi, dime
      logical odryrun
c
      ptr = 1
      do b = 1, nbblock
         dimb = bblock(2,b) - bblock(1,b) + 1
         symb = bblock(3,b)
         do a = 1, nablock
            dima = ablock(2,a) - ablock(1,a) + 1
            syma = ablock(3,a)
            do i = 1, niblock
               dimi = iblock(2,i) - iblock(1,i) + 1
               symi = iblock(3,i)
               syme = ieor(symi,ieor(syma,symb))
               dime = ne_sym(syme)
c
               reclen = dime*dimi*dima*dimb
               if (.not. odryrun) then
                  pointer(1,i,a,b) = ptr
                  pointer(2,i,a,b) = reclen
               endif
               ptr = ptr + reclen
            end do
         end do
      end do
      length = ptr - 1
c
      end
      subroutine uccsdt_2x_pointers(
     $     nablock, ablock,
     $     nbblock, bblock,
     $     niblock, iblock, 
     $     njblock, jblock,
     $     pointer, length, odryrun)
      implicit none
      integer nablock, ablock(3,*) ! [input] Blocking for a
      integer nbblock, bblock(3,*) ! [input] Blocking for b
      integer niblock, iblock(3,*) ! [input] Blocking for i
      integer njblock, jblock(3,*) ! [input] Blocking for i
      integer pointer(2,nablock, nbblock, niblock, njblock) ! [output]
      integer length            ! [output] Total length of the array
      logical odryrun
c
c     Generate the pointer to the start of a block of <ab|ij>
c     for i, j, a, b blocked as specified
c
c     pointer(ablock,bblock,iblock,jblock) = start of block of <ab|ij>
c
c     Eventually the lists must be stored with randomization to avoid
c     hot spts.  Also, eventually the lists must store double precision
c     pointers to use the long GA interface for > 16 GB arrays.
c
c     THESE POINTER LISTS GET TOO BIG ... NEED TO REDUCE TO JUST
C     LAST 3 BLOCK DIMENSIONS AND USE LOGIC TO FIGURE OUT FIRST.
c
c     If (odryrun) then
c     .   do not actually make the pointer array ... just compute
c     .   the size of the GA
c     endif
c
      integer ptr, reclen
      integer a, b, i, j
      integer syma, symb, symi, symj
      integer dima, dimb, dimi, dimj
c
      ptr = 1
      do j = 1, njblock
         dimj = jblock(2,j) - jblock(1,j) + 1
         symj = jblock(3,j)
         do i = 1, niblock
            dimi = iblock(2,i) - iblock(1,i) + 1
            symi = iblock(3,i)
            do b = 1, nbblock
               dimb = bblock(2,b) - bblock(1,b) + 1
               symb = bblock(3,b)
               do a = 1, nablock
                  dima = ablock(2,a) - ablock(1,a) + 1
                  syma = ablock(3,a)
c
                  if (ieor(ieor(ieor(symi,symj),syma),symb) .eq. 0) then
                     reclen = dimi*dimj*dima*dimb
                     if (.not. odryrun) then
                        pointer(1,a,b,i,j) = ptr
                        pointer(2,a,b,i,j) = reclen
                     endif
                     ptr = ptr + reclen
                  else if (.not. odryrun) then
                     pointer(1,a,b,i,j) = -999999999 ! Zero by symmetry
                     pointer(2,a,b,i,j) = -999999999 ! Zero by symmetry
                  end if
               end do
            end do
         end do
      end do
      length = ptr - 1
c
      end
      subroutine uccsdt_triples_sub_blocking_end()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
c     
c     Clean up any MA/GA allocated by set_super_block()
c
      integer list
c
      do list = 23, 1, -1
         if (.not. ma_pop_stack(listinfo(3,list))) call errquit
     $        ('uccsdt_triples_unset_super_block: ma overwrite?',list,
     &       MA_ERR)
         if (listinfo(5,list) .ne. 99999999) then
            if (.not. ga_destroy(listinfo(5,list))) call errquit
     $           ('uccsdt_triples_unset_super_block: ga destroy?', list,
     &       GA_ERR)
         end if
      end do
c
      end
      subroutine uccsdt_triples(rtdb)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cuccsdtP.fh"
#include "cuccsdtps.fh"
#include "amplitudes.fh"
#include "mafdecls.fh"
#include "dra.fh"
      integer rtdb
c     
      integer alo, ahi, blo, bhi, clo, chi, ctop
      integer asuperblk, bsuperblk, csuperblk
      integer spina, spinc
      double precision dummy(5000)
      double precision uccsdt_triples_sub_blocked
      double precision energy(2,2), total
      logical oprint_sub
      double precision ma_used, ga_used
c
      integer max_nsuper_block
      parameter (max_nsuper_block = 10)
      integer nasuper, ncsuper  ! No. of super blocks for current spin case
      integer a_super_block(2,max_nsuper_block) ! Super block info
      integer c_super_block(2,max_nsuper_block)
      integer d_amp
      integer l_ttemp, k_ttemp
      integer noa, nob, nva, nvb
      integer oloa, ohia, olob, ohib, vloa, vhia, vlob, vhib
      integer k_t1ar, l_t1ar, k_t1br, l_t1br
      integer k_t2aar, l_t2aar, k_t2bbr, l_t2bbr
      integer k_t2abr, l_t2abr
      integer lent1a, lent1b, lent2aa, lent2bb, lent2ab
      logical status
c     
      oprint_sub = util_print('sub-blocking', print_debug) .and.
     $     ga_nodeid().eq.0
c
      if (ouccsdtps) call pstat_on(ps_triples)
c     
c     For debug purposes generate the full list of Dirac integrals
c     in a common block
c     
c$$$      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'most', l_most, k_most))
c$$$     $     call errquit('memory for most?',nbf*nbf)
c$$$      if (.not. ma_push_get(mt_dbl, nbf**4, 'moint', l_moint, k_moint))
c$$$     $     call errquit('memory for moint?',nbf**4)
c$$$      call util_transpose(dbl_mb(k_movecs(1)),nbf,
c$$$     $     dbl_mb(k_most),nbf,nbf,nbf)
c$$$      call jan_full_transform_noinit(
c$$$     $     rtdb, basis, 
c$$$     $     nbf, nbf, nbf, nbf,
c$$$     $     nbf, nbf, nbf, nbf,
c$$$     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
c$$$     $     dbl_mb(k_moint), 'Dirac')
C****c     
C****c     For debug purposes generate the full list of t1 and t2 amplitudes
C****c
C****      noa = nc(1) + no(1)
C****      nob = nc(2) + no(2)
C****      nva = nv(1)
C****      nvb = nv(2)
C****      write(6,*)' nva, nvb, noa, nob, nbf: ',
C****     &            nva, nvb, noa, nob, nbf
C****c
C****      if (.not. ma_push_get(mt_dbl, (noa*nva), 't1a', l_t1a, k_t1a))
C****     $     call errquit('ma t1a', (noa*nva))
C****      if (.not. ma_push_get(mt_dbl, (nob*nvb), 't1b', l_t1b, k_t1b))
C****     $     call errquit('ma t1b', (nob*nvb))
C****      if (.not. ma_push_get(mt_dbl, nbf**4, 't2aa', l_t2aa, k_t2aa))
C****     $     call errquit('ma t2_aa', nbf**4)
C****      if (.not. ma_push_get(mt_dbl, nbf**4, 't2bb', l_t2bb, k_t2bb))
C****     $     call errquit('ma t2_bb', nbf**4)
C****      if (.not. ma_push_get(mt_dbl, nbf**4, 't2ab', l_t2ab, k_t2ab))
C****     $     call errquit('ma t2_ab', nbf**4)
C****c
C****      call read_dummy_t(nbf,
C****     $     noa, nob, nva, nvb, 
C****     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
C****     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
C****c      call jan_debug_print('Ta', dbl_mb(k_t1a), nva, noa, 1, 1)
C****c      call jan_debug_print('Tb', dbl_mb(k_t1b), nvb, nob, 1, 1)
C****c      call jan_debug_print('k_t2aa', dbl_mb(k_t2aa), nbf, nbf, nbf, nbf)
C****c      call jan_debug_print('k_t2bb', dbl_mb(k_t2bb), nbf, nbf, nbf, nbf)
C****c      call jan_debug_print('k_t2ab', dbl_mb(k_t2ab), nbf, nbf, nbf, nbf)
C****c
C****c     re-map these t1 and t2 amplitudes into symmetry blocks
C****c
C****      if (.not. ma_push_get(mt_dbl, nbf**4, 'ttemp', l_ttemp, k_ttemp))
C****     $     call errquit('ma ttemp', nbf**4)
C****      call sym_map_amps_t1(
C****     $     dbl_mb(k_t1a), dbl_mb(k_ttemp), nva, noa, 1)
C****      write(6,*)' re-mapped t1a '
C****      call sym_map_amps_t1(
C****     $     dbl_mb(k_t1b), dbl_mb(k_ttemp), nvb, nob, 2)
C****      write(6,*)' re-mapped t1b '
C****      call sym_map_amps_t2(
C****     $     dbl_mb(k_t2aa), dbl_mb(k_ttemp), 1, 1)
C****      write(6,*)' re-mapped t2aa '
C****      call sym_map_amps_t2(
C****     $     dbl_mb(k_t2bb), dbl_mb(k_ttemp), 2, 2)
C****      write(6,*)' re-mapped t2bb '
C****      call sym_map_amps_t2(
C****     $     dbl_mb(k_t2ab), dbl_mb(k_ttemp), 1, 2)
C****      write(6,*)' re-mapped t2ab '
C****      if (.not. ma_pop_stack(l_ttemp)) call errquit('ttemp?',0)
C****c
C****c      call jan_debug_print('Ta', dbl_mb(k_t1a), nva, noa, 1, 1)
C****c      call jan_debug_print('Tb', dbl_mb(k_t1b), nvb, nob, 1, 1)
C****c      call jan_debug_print('k_t2aa', dbl_mb(k_t2aa), nbf, nbf, nbf, nbf)
C****c      call jan_debug_print('k_t2bb', dbl_mb(k_t2bb), nbf, nbf, nbf, nbf)
C****c      call jan_debug_print('k_t2ab', dbl_mb(k_t2ab), nbf, nbf, nbf, nbf)
C****c
C****c     check reading and writing t1 and t2 amps.
C****c
C****c     debug note: write real t1 and t2 to dra in uccsdt_ampfile_test_write_t1
C****c     and uccsdt_ampfile_test_write_t2 read back in uccsdt_ampfile_test_read_t1
C****c     and uccsdt_ampfile_test_read_t2
C****c
C****      oloa = nc(1) + 1
C****      ohia = nc(1) + no(1)
C****      vloa = nc(1) + no(1) + 1
C****      vhia = nmo
C****      olob = nc(2) + 1
C****      ohib = nc(2) + no(2)
C****      vlob = nc(2) + no(2) + 1
C****      vhib = nmo
C****c
C****      lent1a = (vhia-vloa+1)*(ohia-oloa+1)
C****      lent1b = (vhib-vlob+1)*(ohib-olob+1)
C****      lent2aa = lent1a**2
C****      lent2bb = lent1b**2
C****      lent2ab = lent1a*lent1b
C****c
C****      if (.not. ma_push_get(mt_dbl, lent1a, 't1ar', l_t1ar, k_t1ar))
C****     $     call errquit('ma t1ar', lent1a)
C****      if (.not. ma_push_get(mt_dbl, lent1b, 't1br', l_t1br, k_t1br))
C****     $     call errquit('ma t1br', lent1b)
C****      if (.not. ma_push_get(mt_dbl, lent2aa, 't2aar', l_t2aar, k_t2aar))
C****     $     call errquit('ma t2_aar', lent2aa)
C****      if (.not. ma_push_get(mt_dbl, lent2bb, 't2bbr', l_t2bbr, k_t2bbr))
C****     $     call errquit('ma t2_bbr', lent2bb)
C****      if (.not. ma_push_get(mt_dbl, lent2ab, 't2abr', l_t2abr, k_t2abr))
C****     $     call errquit('ma t2_abr', lent2ab)
C****c
C****c     copy t1 and t2 into more compact local arrays (for now)
C****c
C****      call amps_copy(dbl_mb(k_t1a), dbl_mb(k_t1b),
C****     &               dbl_mb(k_t2aa), dbl_mb(k_t2bb), dbl_mb(k_t2ab),
C****     &               dbl_mb(k_t1ar), dbl_mb(k_t1br),
C****     &               dbl_mb(k_t2aar), dbl_mb(k_t2bbr), dbl_mb(k_t2abr),
C****     &               nva, nvb, noa, nob, nbf,
C****     &               oloa, ohia, vloa, vhia, 
C****     &               olob, ohib, vlob, vhib)
C****c
C****      if (.not. uccsdt_ampfile_open('fred','/tmp/fred','new',d_amp))
C****     $     call errquit('open failed',0)
C****      write(6,*)' after opening fred: d_amp = ',d_amp
C****c
C****      call uccsdt_ampfile_test_write_t1(dbl_mb(k_t1ar), 
C****     &                                  dbl_mb(k_t1br),
C****     &                                  oloa, ohia, vloa, vhia, 
C****     &                                  olob, ohib, vlob, vhib,
C****     &                                  d_amp)
C****c
C****      call uccsdt_ampfile_test_write_t2(dbl_mb(k_t2aar), 
C****     &                                  dbl_mb(k_t2bbr), 
C****     &                                  dbl_mb(k_t2abr),
C****     &                                  oloa, ohia, vloa, vhia, 
C****     &                                  olob, ohib, vlob, vhib,
C****     &                                  d_amp)
C****c
C****c     zero out test t2 arrays and try to read in correctly
C****c
C****      call dfill(lent1a, 0.0d0, dbl_mb(k_t1ar), 1)
C****      call dfill(lent1b, 0.0d0, dbl_mb(k_t1br), 1)
C****      call dfill(lent2aa, 0.0d0, dbl_mb(k_t2aar), 1)
C****      call dfill(lent2bb, 0.0d0, dbl_mb(k_t2bbr), 1)
C****      call dfill(lent2ab, 0.0d0, dbl_mb(k_t2abr), 1)
C****c
C****      call uccsdt_ampfile_test_read_t1(dbl_mb(k_t1ar), 
C****     &                                 dbl_mb(k_t1br),
C****     &                                 oloa, ohia, vloa, vhia, 
C****     &                                 olob, ohib, vlob, vhib,
C****     &                                 d_amp)
C****c
C****      call uccsdt_ampfile_test_read_t2(dbl_mb(k_t2aar), 
C****     &                                 dbl_mb(k_t2bbr), 
C****     &                                 dbl_mb(k_t2abr),
C****     &                                 oloa, ohia, vloa, vhia, 
C****     &                                 olob, ohib, vlob, vhib,
C****     &                                 d_amp)
C****c
C****c     diff new t1 and t2 with prior values
C****c
C****      call amps_diff(dbl_mb(k_t1a), dbl_mb(k_t1b),
C****     &               dbl_mb(k_t2aa), dbl_mb(k_t2bb), dbl_mb(k_t2ab),
C****     &               dbl_mb(k_t1ar), dbl_mb(k_t1br),
C****     &               dbl_mb(k_t2aar), dbl_mb(k_t2bbr), dbl_mb(k_t2abr),
C****     &               nva, nvb, noa, nob, nbf,
C****     &               oloa, ohia, vloa, vhia, 
C****     &               olob, ohib, vlob, vhib)
c
      if (.not. uccsdt_ampfile_open('fred','/tmp/fred','old',d_amp))
     $     call errquit('open failed',0, DISK_ERR)
c   
c     Get t1 (both spins) temporarily into local arrays
c
      lent1a = (nmo - nc(1) - no(1)) * no(1)
      if (.not. ma_push_get(mt_dbl, lent1a, 't1a', l_t1a, k_t1a))
     $     call errquit('ma t1a', lent1a, MA_ERR)
      lent1b = (nmo - nc(2) - no(2)) * no(2)
      if (.not. ma_push_get(mt_dbl, lent1b, 't1b', l_t1b, k_t1b))
     $     call errquit('ma t1b', lent1b, MA_ERR)
      call uccsdt_get_unpacked_t1(d_amp, 1, dbl_mb(k_t1a))
      call uccsdt_get_unpacked_t1(d_amp, 2, dbl_mb(k_t1b))
c
      call uccsdt_triples_moints_init()
      call uccsdt_flops_init()
c
c     Loop over the four spin cases
c     
      total = 0.0d0
      do spina = 1, 2
         do spinc = 1, 2
            energy(spina,spinc) = 0.0d0
c     
            if (ouccsdtps) call pstat_on(ps_triples_super_blocking)
c     
            call uccsdt_triples_super_blocking(spina, spinc,
     $           max_nsuper_block, nasuper, ncsuper, a_super_block,
     $           c_super_block)
c     
            if (ouccsdtps) call pstat_off(ps_triples_super_blocking)
c     
            do asuperblk = 1, nasuper
               alo = a_super_block(1,asuperblk)
               ahi = a_super_block(2,asuperblk)
               do bsuperblk = 1, asuperblk
                  blo = a_super_block(1,bsuperblk)
                  bhi = a_super_block(2,bsuperblk)
                  ctop = ncsuper
                  if (spina .eq. spinc) then 
                     do csuperblk = 1, ncsuper
                        clo = c_super_block(1,csuperblk)
                        if (clo .ge. bhi) then
                           ctop = csuperblk - 1
                           goto 10
                        end if
                     end do
 10                  continue
                  end if
                  do csuperblk = 1, ctop
                     clo = c_super_block(1,csuperblk)
                     chi = c_super_block(2,csuperblk)
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_triples_sub_blocking)
c     
                     call uccsdt_triples_sub_blocking
     $                    (alo, ahi, blo, bhi, clo, chi, spina, spinc,
     $                    oprint_sub, .false., ma_used, ga_used)
c     
                     if (ouccsdtps) 
     $                    call pstat_off(ps_triples_sub_blocking)
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_triples_moints)
c     
                     call uccsdt_triples_moints()
c
                     if (ouccsdtps) 
     $                    call pstat_off(ps_triples_moints)
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_triples_amplitudes)
c     
                     call uccsdt_triples_amplitudes(d_amp, spina, spinc)
c
                     if (ouccsdtps) 
     $                    call pstat_off(ps_triples_amplitudes)
c     
                     call dfill(5000, 0.0d0, dummy, 1)
c
                     if (ouccsdtps) 
     $                    call pstat_on(ps_triples_sub_blocked)
c     
                     energy(spina,spinc) = energy(spina,spinc) + 
     $                    uccsdt_triples_sub_blocked()
c     
                     if (ouccsdtps) 
     $                    call pstat_off(ps_triples_sub_blocked)
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_triples_sub_blocking_end)
c     
                     call uccsdt_triples_sub_blocking_end()
c     
                     if (ouccsdtps) 
     $                    call pstat_off(ps_triples_sub_blocking_end)
c     
                  end do
               end do
            end do
            if (ga_nodeid() .eq. 0) then
               write(6,88) spina, spinc, energy(spina,spinc)
 88            format(' spina =',i2,'  spinc =',i2,'  energy =',f20.12/)
               call util_flush(6)
            end if
            total = total + energy(spina, spinc)
         end do
      end do
c
      if (ga_nodeid() .eq. 0) then
         write(6,77) total
 77      format(//' The total (T) energy is ', f20.12,//)
      end if
c
c     Temp write some energy to the RTDB
c
      if ( .not. rtdb_put( rtdb, 'uccsdt:energy', mt_dbl, 1, total) ) 
     $     call errquit('uccsdt: can''t write to rtdb', 14, RTDB_ERR)
c
c     keep dra around for now (testing purposes only)
c
      if (.not. uccsdt_ampfile_close(d_amp,.true.)) call errquit
     $     ('close save failed',0, DISK_ERR)
c
      call uccsdt_triples_moints_tidy()
      call uccsdt_flops_print()
c     
C****      if (.not. ma_pop_stack(l_t2abr)) call errquit('t2abr?',0)
C****      if (.not. ma_pop_stack(l_t2bbr)) call errquit('t2bbr?',0)
C****      if (.not. ma_pop_stack(l_t2aar)) call errquit('t2aar?',0)
C****      if (.not. ma_pop_stack(l_t1br)) call errquit('t1br?',0)
C****      if (.not. ma_pop_stack(l_t1ar)) call errquit('t1ar?',0)
C****      if (.not. ma_pop_stack(l_t2ab)) call errquit('t2ab?',0)
C****      if (.not. ma_pop_stack(l_t2bb)) call errquit('t2bb?',0)
C****      if (.not. ma_pop_stack(l_t2aa)) call errquit('t2aa?',0)
      if (.not. ma_pop_stack(l_t1b)) call errquit('t1b?',0, MA_ERR)
      if (.not. ma_pop_stack(l_t1a)) call errquit('t1a?',0, MA_ERR)
*      if (.not. ma_pop_stack(l_moint)) call errquit('moint?',0)
*      if (.not. ma_pop_stack(l_most)) call errquit('most?',0)
c     
      if (ouccsdtps) call pstat_off(ps_triples)
c     
      end
      subroutine sym_map_amps_t1(
     $     t1, t_buf, nv1, no1, spin)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
      integer nv1, no1
      integer v, o
      integer spin
      integer v_old_abs, v_new_rel
      double precision t1(nv1, no1)
      double precision t_buf(nv1, no1)
c
c     copy t1 to t_buf, zero t1 and move t_buf back 
c     to proper symmetry locations in t1
c
      call dfill(nv1*no1, 0.0d0, t_buf,1)
      call dcopy(nv1*no1, t1, 1, t_buf, 1)
      call dfill(nv1*no1, 0.0d0, t1, 1)
      do v = 1, nv1
         do o = 1, no1
c            write(6,*)' v, o, map(v,spin), map(o,spin): ',
c     &                  v, o, map(v,spin), map(o,spin)
            v_old_abs = v + no1
            v_new_rel = map(v_old_abs,spin) - no1
c            write(6,*)' v_old_abs, o, v_new_rel, map(o,spin): ',
c     &                  v_old_abs, o, v_new_rel, map(o,spin)
            t1(v_new_rel,map(o,spin)) = t_buf(v,o)
         end do
      end do
c      if (.not. ma_verify_allocator_stuff())
c     $     call errquit('x',1)
      return
      end
      subroutine sym_map_amps_t2(
     $     t2, t_buf, spin1, spin2)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
      double precision t2(nbf, nbf, nbf, nbf)
      double precision t_buf(nbf, nbf, nbf, nbf)
      integer v1, v2, o1, o2
      integer vv2, oo1, oo2
      integer spin1, spin2
      call dfill(nbf**4, 0.0d0, t_buf, 1)
      call dcopy(nbf**4, t2, 1, t_buf, 1)
      call dfill(nbf**4, 0.0d0, t2, 1)
      do o2 = 1, nbf
         oo2 = map(o2,spin2)
         do o1 = 1, nbf
            oo1 = map(o1,spin1)
            do v2 = 1, nbf
               vv2 = map(v2,spin2)
               do v1 = 1, nbf
c
c            write(6,*)' v1, v2, map(v1,spin1), map(v2,spin2): ',
c     &                  v1, v2, map(v1,spin1), map(v2,spin2)
c            write(6,*)' o1, o2, map(o1,spin1), map(o2,spin2): ',
c     &                  o1, o2, map(o1,spin1), map(o2,spin2)
                  t2(map(v1,spin1),vv2,oo1,oo2) = t_buf(v1,v2,o1,o2)
c
               end do
            end do
         end do
      end do
c      if (.not. ma_verify_allocator_stuff())
c     $     call errquit('x',2)
      return
      end
      subroutine uccsdt_flops_init()
      implicit none
#include "cflopsP.fh"
c
      call dfill(19, 0d0, flops, 1)
      call dfill(19, 0d0, times, 1)
      call dfill(19, 0d0, calls, 1)
c
      end
      subroutine uccsdt_flops_print()
      implicit none
#include "cflopsP.fh"
#include "util.fh"
c
      double precision overhead, start, mflops, used
      integer ntime, i, j
c
c     Estimate the timing overhead
c
      start = util_cpusec()
      ntime = 0
      do i = 1, 100
         do j = 1, 1000
            call uccsdt_flops_enter()
            call uccsdt_flops_exit(19,0d0)
         end do
         ntime = ntime + 100
         overhead = util_cpusec() - start
         if (overhead .gt. 1.0d0) goto 10
      enddo
 10   overhead = overhead / dble(ntime)
      write(6,1) overhead
 1    format(//,'Triples kernel performance information ',//,
     $     ' Timing overhead = ', 1p,d9.1,'s per call',/)
c
      write(6,4) 
 4    format(' Routine   Calls     FLOPs     Time     MFLOP/s',/,
     $       ' -------  -------  --------  --------  --------')
      do i = 1, 18
         used = times(i) - calls(i)*overhead
         if (used .gt. 0d0) then
            mflops = 1d-6*flops(i)/used
            write(6,2) i, calls(i), flops(i), used, mflops
 2          format(1x,i5,2x,f9.0, 1x, 1p, d9.1, 0p, 1x, f9.2, 1x, f9.2)
         end if
      end do
c
      end
      subroutine uccsdt_flops_enter()
      implicit none
#include "cflopsP.fh"
#include "util.fh"
c
      enter = util_cpusec()
c
      end
      subroutine uccsdt_flops_exit(i, flopcount)
      implicit none
#include "cflopsP.fh"
#include "util.fh"
      integer i
      double precision flopcount
c
      times(i) = times(i) + util_cpusec() - enter
      flops(i) = flops(i) + flopcount
      calls(i) = calls(i) + 1
c
      end
      subroutine uccsdt_lambda(d_amp)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer d_amp             ! [input] Disk file with amplitudes
c
c     Make the particle-hole transformed MO coefficients
c     putting the results into the pre-allocated global arrays
c     g_hole(1:2), g_part(1:2)
c
      integer l_t1, k_t1, l_p, k_p, l_h, k_h, l_work, k_work, spin
c
      if (.not. ma_push_get(mt_dbl, nmo*nmo, 'work',
     $     l_work, k_work)) call errquit('uccsdt_lambda: ma?',1, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'work',
     $     l_h, k_h)) call errquit('uccsdt_lambda: ma?',2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'work',
     $     l_p, k_p)) call errquit('uccsdt_lambda: ma?',3, MA_ERR)
c
      do spin = 1, 2
c     
c     Get the symmetry unpacked T1 into a local array, transform
c     and store the results.
c
         if (.not. ma_push_get(mt_dbl, no(spin)*nv(spin), 'tmp',
     $        l_t1, k_t1)) call errquit('uccsdt_lambda: MA?',0, MA_ERR)
         call uccsdt_get_unpacked_t1(d_amp, spin, dbl_mb(k_t1))
         if (ga_nodeid() .eq. 0) then
            call uccsdt_do_lambda(nbf, nmo, nc(spin), no(spin),
     $           nv(spin), dbl_mb(k_t1), dbl_mb(k_movecs(spin)),
     $           dbl_mb(k_p), dbl_mb(k_h), dbl_mb(k_work))
            call ga_put(g_hole(spin), 1, nbf, 1, nmo, dbl_mb(k_h), nbf)
            call ga_put(g_part(spin), 1, nbf, 1, nmo, dbl_mb(k_p), nbf)
         end if
         if (.not. ma_pop_stack(l_t1)) call errquit
     $        ('uccsdt_lambda: ma?', 55, MA_ERR)
         call ga_sync()
      end do
      if (.not. ma_pop_stack(l_p)) call errquit
     $     ('uccsdt_lambda: ma?', 56, MA_ERR)
      if (.not. ma_pop_stack(l_h)) call errquit
     $     ('uccsdt_lambda: ma?', 57, MA_ERR)
      if (.not. ma_pop_stack(l_work)) call errquit
     $     ('uccsdt_lambda: ma?', 58, MA_ERR)
c
      do spin = 1, 2
         write(6,*) ' HOLE MOs ', spin
         call ga_print(g_hole(spin))
         write(6,*) ' PARTICLE MOs ', spin
         call ga_print(g_part(spin))
      end do
c
      end
      subroutine uccsdt_do_lambda(nbf, nmo, nc, no, nv, 
     $        t1, c, p, h, work)
      implicit none
c
      integer nbf, nmo, nc, no, nv
      double precision t1(nv, no)
      double precision c(nbf, nmo)
      double precision p(nbf, nmo)
      double precision h(nbf, nmo)
      double precision work(nmo, nmo)
c
      integer i, a
c
c     Particle transformation
c
      call dfill(nmo*nmo, 0.0d0, work, 1)
      call dfill(nmo, 1.0d0, work, nmo+1)
      do i = 1, no
         do a = 1, nv
            work(i+nc,a+no+nc) = work(i+nc,a+no+nc) - t1(a,i)
         end do
      end do
      call dgemm('n','n',nbf,nmo,nmo,
     $     1.0d0, c, nbf, work, nmo, 
     $     0.0d0, p, nmo)
c
c     Hole transformation
c
      call dfill(nmo*nmo, 0.0d0, work, 1)
      call dfill(nmo, 1.0d0, work, nmo+1)
      do i = 1, no
         do a = 1, nv
            work(a+nc+no,i+nc) = work(a+nc+no,i+nc) + t1(a,i)
         end do
      end do
      call dgemm('n','n',nbf,nmo,nmo,
     $     1.0d0, c, nmo, work, nmo, 
     $     0.0d0, h, nmo)
c
      end
      subroutine uccsdt_get_unpacked_t1(d_amp, spin, t1)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer d_amp             ! [input]
      integer spin              ! [input]
      double precision t1(*)    ! [output]
c
c     Convenience routine to read the symmetry packed single
c     excitation amplitudes from disk and unpack them into
c     a local array of dimension (nv(spin),no(spin)).
c
      integer g_t1, sym, l_tmp, k_tmp, ai
c
c     Read the packed amplitudes from disk into a GA and 
c     fetch into a local array
c     
      if (.not. ma_push_get(mt_dbl, t1_len(spin), 'tmp',
     $     l_tmp, k_tmp)) call errquit('uccsdt_get_t1: MA?',0, MA_ERR)
      if (.not. uccsdt_ampfile_read_t1(d_amp, spin, g_t1, .true.))
     $     call errquit('uccsdt_get_t1: reading t1 failed', d_amp,
     &       DISK_ERR)
      if (ga_nodeid() .eq. 0) call ga_get
     $     (g_t1, 1, t1_len(spin), 1, 1, dbl_mb(k_tmp), 1)
      call ga_brdcst(1, dbl_mb(k_tmp), 8*t1_len(spin), 0)
      if (.not. ga_destroy(g_t1)) call errquit
     $     ('uccsdt_get_t1: GA?',0, GA_ERR)
c
c     Transpose into the symmetry unpacked output
c
      call dfill(no(spin)*nv(spin), 0d0, t1, 1)
      do sym = 0, 7
         if (no_sym(sym,spin)*nv_sym(sym,spin) .gt. 0) then
            ai = 1 + v_sym(1,sym,spin)-v_sym(1,0,spin) + 
     $           (o_sym(1,sym,spin)-o_sym(1,0,spin))*nv(spin)
c            write(6,*) ' T1 block ', sym
c            call output(dbl_mb(k_tmp), 1, no_sym(sym,spin), 
c     $           1, nv_sym(sym,spin), no_sym(sym,spin), 
c     $           nv_sym(sym,spin), 1)
            call util_transpose(
     $           dbl_mb(k_tmp), no_sym(sym,spin), 
     $           t1(ai), nv(spin),
     $           no_sym(sym,spin), nv_sym(sym,spin))
            k_tmp = k_tmp + no_sym(sym,spin)*nv_sym(sym,spin)
         end if
      end do
c
c      write(6,*) ' Full T1'
c      call output(t1, 1, nv(spin), 1, no(spin), nv(spin), no(spin), 1)
c
      if (.not. ma_pop_stack(l_tmp)) call errquit
     $     ('uccsdt_get_t1: MA corrupt?',0, MA_ERR)
c
      end

c $Id$
