      logical function uccsdtest(rtdb)
*     
*     $Id: uccsdtest.F,v 1.8 2000-04-26 00:34:13 d3h449 Exp $
*
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
c
      integer basis, geom, nbf
      character*255 movecs      ! Name of movector file
      character*80 title, name_of_basis, scftype
      integer nbf_file, nsets, nmo_file(2)
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      integer nmo, nalpha, nbeta, g_tmp, l_occ, k_occ,
     $     l_aeval, k_aeval, l_beval, k_beval, l_amos, k_amos,
     $     l_bmos, k_bmos
      integer noa, nob, nva, nvb
      integer l_aa, l_ab, l_ba, l_bb, k_iaa, k_iab, k_iba, k_ibb
      integer l_t2aa, l_t2ab, l_t2bb, l_t1a, l_t1b, l_t2ba
      integer k_t2aa, k_t2ab, k_t2bb, k_t1a, k_t1b, k_t2ba
      integer l_r2aa, l_r2ab, l_r2bb, l_r1a, l_r1b
      integer k_r2aa, k_r2ab, k_r2bb, k_r1a, k_r1b
      integer l_fa, k_fa, l_fb, k_fb
      integer l_alp, l_alh, l_blp, l_blh, l_ha, l_hb, l_work
      integer k_alp, k_alh, k_blp, k_blh, k_ha, k_hb, k_work
      integer l_z1, l_z2, l_z3, l_z4, l_z5, l_z6, l_z7, l_z8, 
     $     l_x1, l_x2, l_x3
      integer k_z1, k_z2, k_z3, k_z4, k_z5, k_z6, k_z7, k_z8, 
     $     k_x1, k_x2, k_x3
      integer iter
c
      double precision energyaaa, energybbb, energyaab, energybba
      double precision uccsdtest_triples_pure, uccsdtest_triples_mixed
      double precision uccsdtest_triples_pure_blocked
      logical int_normalize
      external int_normalize
c     
c     load the geometry/basis set and get info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0)
      if (.not.int_normalize(rtdb, basis))
     $           call errquit('scf:int_normalize failed', 0)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c     Read the MO vectors and evals from a UHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('rjhtest: failed to read movecs header',0)
      if (nsets.ne.2) call errquit('UHF ONLY!',0)
      nmo = nmo_file(1)
      if (nmo .ne. nbf) call errquit('NO LINEAR DEPENDNCE PERMITTED',0)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      if (nmo_file(1).ne.nmo_file(2)) call errquit('Uh?',0)
      if (.not. rtdb_get(rtdb, 'scf:nalpha', mt_int, 1, nalpha))
     $     call errquit('nalpha?',0)
      if (.not. rtdb_get(rtdb, 'scf:nbeta', mt_int, 1, nbeta))
     $     call errquit('nbeta?',0)

      noa = nalpha
      nob = nbeta
      nva = nmo - nalpha
      nvb = nmo - nbeta

      write(6,*) ' No. of alpha electrons:    ', nalpha
      write(6,*) ' No. of beta  electrons:    ', nbeta
      write(6,*) ' No. of molecular orbitals: ', nmo
      write(6,*) ' No. of basis functions:    ', nbf
c
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 1, 1, g_tmp))
     &     call errquit('scf_v_g: tmp', 0)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf)
      if (.not. ma_push_get(mt_dbl, nbf,'aeval',l_aeval, k_aeval))
     $     call errquit('ma eval', nbf)
      if (.not. ma_push_get(mt_dbl, nbf,'beval',l_beval, k_beval))
     $     call errquit('ma eval', nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'amos', l_amos, k_amos))
     $     call errquit('ma amos', nbf*nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'bmos', l_bmos, k_bmos))
     $     call errquit('ma bmos', nbf*nbf)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_aeval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_amos), nbf)
      if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), dbl_mb(k_beval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_bmos), nbf)
c
      write(6,*) ' Alpha eigenvalues '
      call output(dbl_mb(k_aeval),1,nmo,1,1,nmo,1,1)
      write(6,*) ' Beta  eigenvalues '
      call output(dbl_mb(k_beval),1,nmo,1,1,nmo,1,1)
c$$$      write(6,*) ' Alpha MOs'
c$$$      call output(dbl_mb(k_amos),1,nbf,1,nmo,nbf,nmo,1)
c$$$      write(6,*) ' Beta  MOs'
c$$$      call output(dbl_mb(k_bmos),1,nbf,1,nmo,nbf,nmo,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0)
c
      call util_inplace_transpose(nmo, dbl_mb(k_amos))
      call util_inplace_transpose(nmo, dbl_mb(k_bmos))
c
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_aa, k_iaa))
     $     call errquit('ma aa moints', nbf**4)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_ab, k_iab))
     $     call errquit('ma ab moints', nbf**4)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_ba, k_iba))
     $     call errquit('ma ba moints', nbf**4)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_bb, k_ibb))
     $     call errquit('ma bb moints', nbf**4)
c
c     Canonical UHF MP2 as a test and to get starting amplitudes.
c
      if (.not. ma_push_get(mt_dbl, (noa*nva), 't1', l_t1a, k_t1a))
     $     call errquit('ma t1a', (noa*nva))
      if (.not. ma_push_get(mt_dbl, (nob*nvb), 't1', l_t1b, k_t1b))
     $     call errquit('ma t1b', (nob*nvb))
      if (.not. ma_push_get(mt_dbl, (noa*nva)**2, 't2', l_t2aa, k_t2aa))
     $     call errquit('ma t2_aa', (noa*nva)**2)
      if (.not. ma_push_get(mt_dbl, (nob*nvb)**2, 't2', l_t2bb, k_t2bb))
     $     call errquit('ma t2_bb', (nob*nvb)**2)
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     't2', l_t2ab, k_t2ab))
     $     call errquit('ma t2_ab', (noa*nva)*(nob*nvb))
c
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_amos),dbl_mb(k_amos),
     $        dbl_mb(k_amos),dbl_mb(k_amos), 
     $        dbl_mb(k_iaa), 'RightAsymDirac')
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_bmos),dbl_mb(k_bmos),
     $        dbl_mb(k_bmos),dbl_mb(k_bmos), 
     $        dbl_mb(k_ibb), 'RightAsymDirac')
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_amos),dbl_mb(k_bmos),
     $        dbl_mb(k_amos),dbl_mb(k_bmos), 
     $        dbl_mb(k_iab), 'Dirac')
      call uccsdtest_mp2(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_aeval), dbl_mb(k_beval),
     $     dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
c
      call dfill((noa*nva)**2,0.0d0,dbl_mb(k_t2aa),1)
      call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_t2bb),1)
      call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_t2ab),1)
      call dfill((noa*nva), 0.0d0, dbl_mb(k_t1a), 1)
      call dfill((nob*nvb), 0.0d0, dbl_mb(k_t1b), 1)
c
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fa', l_fa, k_fa))
     $     call errquit('ma fa', nmo**2)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fb', l_fb, k_fb))
     $     call errquit('ma fb', nmo**2)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fa', l_ha, k_ha))
     $     call errquit('ma ha', nmo**2)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fb', l_hb, k_hb))
     $     call errquit('ma hb', nmo**2)
c$$$      call dfill(nmo*nmo, 0.0d0, dbl_mb(k_fa), 1)
c$$$      call dfill(nmo*nmo, 0.0d0, dbl_mb(k_fb), 1)
c$$$      call dcopy(nmo, dbl_mb(k_aeval), 1, dbl_mb(k_fa), nmo+1)
c$$$      call dcopy(nmo, dbl_mb(k_beval), 1, dbl_mb(k_fb), nmo+1)
c
c     Space for the residual vectors
c
      if (.not. ma_push_get(mt_dbl, (noa*nva), 'r1', l_r1a, k_r1a))
     $     call errquit('ma r1a', (noa*nva))
      if (.not. ma_push_get(mt_dbl, (nob*nvb), 'r1', l_r1b, k_r1b))
     $     call errquit('ma r1b', (nob*nvb))
      if (.not. ma_push_get(mt_dbl, (noa*nva)**2, 'r2', l_r2aa, k_r2aa))
     $     call errquit('ma r2_aa', (noa*nva)**2)
      if (.not. ma_push_get(mt_dbl, (nob*nvb)**2, 'r2', l_r2bb, k_r2bb))
     $     call errquit('ma r2_bb', (nob*nvb)**2)
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     'r2', l_r2ab, k_r2ab))
     $     call errquit('ma r2_ab', (noa*nva)*(nob*nvb))
c
c     Temporary workspace for transposed t2ab and
c     particle & hole transformation matrices and workspace
c
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     't2', l_t2ba, k_t2ba))
     $     call errquit('ma t2_ba', (noa*nva)*(nob*nvb))
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'a lamda p', l_alp, k_alp))
     $     call errquit('ma alp', nbf*nmo)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'a lamda h', l_alh, k_alh))
     $     call errquit('ma alh', nbf*nmo)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'b lamda p', l_blp, k_blp))
     $     call errquit('ma alp', nmo*nmo)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'b lamda h', l_blh, k_blh))
     $     call errquit('ma alh', nbf*nmo)
      if (.not. ma_push_get(mt_dbl, nmo*nmo, 'work', l_work, k_work))
     $     call errquit('ma work', nmo*nmo)
c
c     For the n^6 algorithm we need some more intermediates
c
      if (.not. ma_push_get(mt_dbl, noa*nva*noa*nva, 'z1', l_z1, k_z1))
     $     call errquit(' ma z1 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*nva*noa*nva, 'z2', l_z2, k_z2))
     $     call errquit(' ma z2 ', 0)
      if (.not. ma_push_get(mt_dbl, nob*nvb*noa*nva, 'z3', l_z3, k_z3))
     $     call errquit(' ma z3 ', 0)
      if (.not. ma_push_get(mt_dbl, nob*nvb*nob*nvb, 'z4', l_z4, k_z4))
     $     call errquit(' ma z4 ', 0)
      if (.not. ma_push_get(mt_dbl, nob*nvb*nob*nvb, 'z5', l_z5, k_z5))
     $     call errquit(' ma z5 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*nva*nob*nvb, 'z6', l_z6, k_z6))
     $     call errquit(' ma z6 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*nva*nob*nvb, 'z7', l_z7, k_z7))
     $     call errquit(' ma z7 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*nvb*noa*nvb, 'z8', l_z8, k_z8))
     $     call errquit(' ma z8 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*noa*noa*noa, 'x1', l_x1, k_x1))
     $     call errquit(' ma x1 ', 0)
      if (.not. ma_push_get(mt_dbl, nob*nob*nob*nob, 'x2', l_x2, k_x2))
     $     call errquit(' ma x2 ', 0)
      if (.not. ma_push_get(mt_dbl, noa*nob*noa*nob, 'x3', l_x3, k_x3))
     $     call errquit(' ma x3 ', 0)
c
c     Iterate
c
      call dfill((noa*nva)**2,0.0d0,dbl_mb(k_t2aa),1)
      call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_t2bb),1)
      call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_t2ab),1)
      call dfill((noa*nva), 0.0d0, dbl_mb(k_t1a), 1)
      call dfill((nob*nvb), 0.0d0, dbl_mb(k_t1b), 1)
c
      do iter = 1, 3
c
c     Transform the MO coefficients with the particle and hole matrices
c
         call uccsdtest_lambda(nbf, nmo, noa, nva, 
     $        dbl_mb(k_t1a), dbl_mb(k_amos), 
     $        dbl_mb(k_alp), dbl_mb(k_alh), dbl_mb(k_work))
         call uccsdtest_lambda(nbf, nmo, nob, nvb, 
     $        dbl_mb(k_t1b), dbl_mb(k_bmos), 
     $        dbl_mb(k_blp), dbl_mb(k_blh), dbl_mb(k_work))
c
c     Transform the integrals
c
c     Generate the 3 sets of integrals in Dirac order
c
c     For simplicity also generate transpose of mixed spin integrals
c
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_alp),dbl_mb(k_alp),
     $        dbl_mb(k_alh),dbl_mb(k_alh), 
     $        dbl_mb(k_iaa), 'RightAsymDirac')
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_blp),dbl_mb(k_blp),
     $        dbl_mb(k_blh),dbl_mb(k_blh), 
     $        dbl_mb(k_ibb), 'RightAsymDirac')
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_alp),dbl_mb(k_blp),
     $        dbl_mb(k_alh),dbl_mb(k_blh), 
     $        dbl_mb(k_iab), 'Dirac')
         call jan_full_transform(
     $        rtdb, basis, 
     $        nmo, nmo, nmo, nmo,
     $        nmo, nmo, nmo, nmo,
     $        dbl_mb(k_blp),dbl_mb(k_alp),
     $        dbl_mb(k_blh),dbl_mb(k_alh), 
     $        dbl_mb(k_iba), 'Dirac')
c
         call jan_h(rtdb, basis, nmo, nmo, nmo, nmo,
     $        dbl_mb(k_alp), dbl_mb(k_alh), dbl_mb(k_ha))
         call jan_h(rtdb, basis, nmo, nmo, nmo, nmo,
     $        dbl_mb(k_blp), dbl_mb(k_blh), dbl_mb(k_hb))
c
c     Assemble the transformed Fock matrices
c
         call uccsdtest_fock(nmo, noa, nob,
     $        dbl_mb(k_ha), dbl_mb(k_iaa), dbl_mb(k_iab), 
     $        dbl_mb(k_fa))
         call uccsdtest_fock(nmo, nob, noa,
     $        dbl_mb(k_hb), dbl_mb(k_ibb), dbl_mb(k_iba), 
     $        dbl_mb(k_fb))
c
         call dfill((noa*nva)**2,0.0d0,dbl_mb(k_r2aa),1)
         call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_r2bb),1)
         call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_r2ab),1)
         call dfill((noa*nva), 0.0d0, dbl_mb(k_r1a), 1)
         call dfill((nob*nvb), 0.0d0, dbl_mb(k_r1b), 1)
c
c     Alpha pure spin
c
         call uccsdtest_product_pure(
     $        nmo,
     $        noa, nob, nva, nvb, 
     $        dbl_mb(k_fa), dbl_mb(k_fb),
     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
     $        dbl_mb(k_r2aa))
c
c     Beta pure spin ... use the transposed mixed spin amplitudes
c     and integrals and then call the pure spin routine with spins flipped

         call uccsdtest_transpose_tab(noa, nob, nva, nvb, 
     $        dbl_mb(k_t2ab), dbl_mb(k_t2ba))
         call uccsdtest_product_pure(
     $        nmo,
     $        nob, noa, nvb, nva, 
     $        dbl_mb(k_fb), dbl_mb(k_fa),
     $        dbl_mb(k_ibb), dbl_mb(k_iaa), dbl_mb(k_iba),
     $        dbl_mb(k_t1b), dbl_mb(k_t1a),
     $        dbl_mb(k_t2bb),dbl_mb(k_t2aa),dbl_mb(k_t2ba),
     $        dbl_mb(k_r1b), dbl_mb(k_r1a),
     $        dbl_mb(k_r2bb))
c
c     Mixed
c
         call uccsdtest_product_mixed(
     $        nmo,
     $        noa, nob, nva, nvb, 
     $        dbl_mb(k_fa), dbl_mb(k_fb),
     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
     $        dbl_mb(k_r2ab))
c
c$$$         call jan_debug_print('Ta', dbl_mb(k_t1a), nva, noa, 1, 1)
c$$$         call jan_debug_print('Tb', dbl_mb(k_t1b), nvb, nob, 1, 1)
c$$$         call jan_debug_print('Taa',dbl_mb(k_t2aa), nva, nva, noa, noa)
c$$$         call jan_debug_print('Tbb',dbl_mb(k_t2bb), nvb, nvb, nob, nob)
c$$$         call jan_debug_print('Tab',dbl_mb(k_t2ab), nva, nvb, noa, nob)
c$$$         call jan_debug_print('Ra', dbl_mb(k_r1a), nva, noa, 1, 1)
c$$$         call jan_debug_print('Rb', dbl_mb(k_r1b), nvb, nob, 1, 1)
c$$$         call jan_debug_print('Raa',dbl_mb(k_r2aa), nva, nva, noa, noa)
c$$$         call jan_debug_print('Rbb',dbl_mb(k_r2bb), nvb, nvb, nob, nob)
c$$$         call jan_debug_print('Rab',dbl_mb(k_r2ab), nva, nvb, noa, nob)
c
c     Now run the n^6 algorithm
c
         call uccsdtest_intermediates(
     $        nmo,
     $        noa, nob, nva, nvb, 
     $        dbl_mb(k_fa), dbl_mb(k_fb),
     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $        dbl_mb(k_z1), dbl_mb(k_z2), dbl_mb(k_z3), dbl_mb(k_z4), 
     $        dbl_mb(k_z5), dbl_mb(k_z6), dbl_mb(k_z7), dbl_mb(k_z8), 
     $        dbl_mb(k_x1), dbl_mb(k_x2), dbl_mb(k_x3))
         call dfill((noa*nva)**2,0.0d0,dbl_mb(k_r2aa),1)
         call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_r2bb),1)
         call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_r2ab),1)
         call dfill((noa*nva), 0.0d0, dbl_mb(k_r1a), 1)
         call dfill((nob*nvb), 0.0d0, dbl_mb(k_r1b), 1)
         call uccsdtest_new_product(
     $        nmo,
     $        noa, nob, nva, nvb, 
     $        dbl_mb(k_fa), dbl_mb(k_fb),
     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $        dbl_mb(k_z1), dbl_mb(k_z2), dbl_mb(k_z3), dbl_mb(k_z4), 
     $        dbl_mb(k_z5), dbl_mb(k_z6), dbl_mb(k_z7), dbl_mb(k_z8), 
     $        dbl_mb(k_x1), dbl_mb(k_x2), dbl_mb(k_x3),
     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
     $        dbl_mb(k_r2aa),dbl_mb(k_r2bb),dbl_mb(k_r2ab))
c$$$         call jan_debug_print('Na', dbl_mb(k_r1a), nva, noa, 1, 1)
c$$$         call jan_debug_print('Nb', dbl_mb(k_r1b), nvb, nob, 1, 1)
c$$$         call jan_debug_print('Naa',dbl_mb(k_r2aa), nva, nva, noa, noa)
c$$$         call jan_debug_print('Nbb',dbl_mb(k_r2bb), nvb, nvb, nob, nob)
c$$$         call jan_debug_print('Nab',dbl_mb(k_r2ab), nva, nvb, noa, nob)
c
c     Update
c
         call uccsdtest_update(
     $        nmo,
     $        noa, nob, nva, nvb, 
     $        dbl_mb(k_fa), dbl_mb(k_fb),
     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
     $        dbl_mb(k_r2aa),dbl_mb(k_r2bb),dbl_mb(k_r2ab),
     $        dbl_mb(k_aeval), dbl_mb(k_beval))
c
c$$$         call jan_debug_print('Taa',dbl_mb(k_t2aa), nva, nva, noa, noa)
c$$$         call jan_debug_print('Tbb',dbl_mb(k_t2bb), nvb, nvb, nob, nob)
c$$$         call jan_debug_print('Tab',dbl_mb(k_t2ab), nva, nvb, noa, nob)
c
      end do
c
c     Triples
c
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_amos),dbl_mb(k_amos),
     $     dbl_mb(k_amos),dbl_mb(k_amos), 
     $     dbl_mb(k_iaa), 'RightAsymDirac')
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_bmos),dbl_mb(k_bmos),
     $     dbl_mb(k_bmos),dbl_mb(k_bmos), 
     $     dbl_mb(k_ibb), 'RightAsymDirac')
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_amos),dbl_mb(k_bmos),
     $     dbl_mb(k_amos),dbl_mb(k_bmos), 
     $     dbl_mb(k_iab), 'Dirac')
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_bmos),dbl_mb(k_amos),
     $     dbl_mb(k_bmos),dbl_mb(k_amos), 
     $     dbl_mb(k_iba), 'Dirac')
c
      energyaaa = uccsdtest_triples_pure(
     $     nmo,
     $     noa, nva, 
     $     dbl_mb(k_iaa), 
     $     dbl_mb(k_t1a), 
     $     dbl_mb(k_t2aa),
     $     dbl_mb(k_aeval))
      energybbb = uccsdtest_triples_pure(
     $     nmo,
     $     nob, nvb, 
     $     dbl_mb(k_ibb), 
     $     dbl_mb(k_t1b), 
     $     dbl_mb(k_t2bb),
     $     dbl_mb(k_beval))
      energyaab = uccsdtest_triples_mixed(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $     dbl_mb(k_aeval), dbl_mb(k_beval))
      call uccsdtest_transpose_tab(noa, nob, nva, nvb, 
     $     dbl_mb(k_t2ab), dbl_mb(k_t2ba))
      energybba = uccsdtest_triples_mixed(
     $     nmo,
     $     nob, noa, nvb, nva, 
     $     dbl_mb(k_ibb), dbl_mb(k_iaa), dbl_mb(k_iba),
     $     dbl_mb(k_t1b), dbl_mb(k_t1a),
     $     dbl_mb(k_t2bb),dbl_mb(k_t2aa),dbl_mb(k_t2ba),
     $     dbl_mb(k_beval), dbl_mb(k_aeval))
c
      write(6,*) energyaaa, energybbb
      write(6,*) energyaab, energybba
      write(6,*) energyaaa+energybbb+energyaab+energybba
c
      call uccsdt_initialize(nmo, noa, nob, nva, nvb, 
     $     dbl_mb(k_beval), dbl_mb(k_aeval),
     $     k_ibb, k_iaa, k_iba,
     $     k_t1b, k_t1a,
     $     k_t2bb,k_t2aa,k_t2ab)
c
      energyaaa = uccsdtest_triples_pure_blocked(
     $     dbl_mb(k_iaa), 
     $     dbl_mb(k_t1a), 
     $     dbl_mb(k_t2aa),
     $     dbl_mb(k_aeval), 1, nva, noa)
c
      write(6,*) energyaaa
c
      energybbb = uccsdtest_triples_pure_blocked(
     $     dbl_mb(k_ibb), 
     $     dbl_mb(k_t1b), 
     $     dbl_mb(k_t2bb),
     $     dbl_mb(k_beval), 2, nvb, nob)
c
      write(6,*) energybbb
c
      if (.not. ma_chop_stack(l_occ)) call errquit('ma?',0)
c
      uccsdtest = .true.
c
      end
      subroutine uccsdtest_product_pure(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     fa, fb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     r1a, r1b,
     $     r2aa)
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision fa(nmo,nmo), fb(nmo,nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
      double precision r1a(nva, noa)
      double precision r1b(nvb, nob)
      double precision r2aa(nva, nva, noa, noa)
c
      integer i, j, m, n
      integer a, b, e, f
c
c     The formulae in comments denote beta spin with an asterisk
c
c     R1 
c
      do m = 1, noa
         do e = 1, nva
c
            r1a(e,m) = fa(e+noa,m)
c
            do j = 1, noa
               do b = 1, nva
                  r1a(e,m) = r1a(e,m) + fa(j,b+noa)*t2aa(e,b,m,j)
               end do
            end do
            do j = 1, nob
               do b = 1, nvb
                  r1a(e,m) = r1a(e,m) + fb(j,b+nob)*t2ab(e,b,m,j)
               end do
            end do
c
            do j = 1, noa
               do i = 1, noa
                  do b = 1, nva
                     r1a(e,m) = r1a(e,m) 
     $                    - 0.5d0*iaa(i,j,m,b+noa)*t2aa(e,b,i,j)
                  end do
               end do
            end do
c
            do j = 1, nob
               do i = 1, noa
                  do b = 1, nvb
                     r1a(e,m) = r1a(e,m)
     $                    - iab(i,j,m,b+nob)*t2ab(e,b,i,j)
                  end do
               end do
            end do
c
            do i = 1, noa
               do b = 1, nva
                  do a = 1, nva
                     r1a(e,m) = r1a(e,m)
     $                    + 0.5d0*iaa(i,e+noa,a+noa,b+noa)*t2aa(a,b,i,m)
                  end do
               end do
            end do
c
            do i = 1, nob
               do b = 1, nva
                  do a = 1, nvb
                     r1a(e,m) = r1a(e,m)
     $                    + iab(e+noa,i,b+noa,a+nob)*t2ab(b,a,m,i)
                  end do
               end do
            end do
c
         end do
      end do
c
c     R2 
c
      do n = 1, noa
         do m = 1, noa
            do f = 1, nva
               do e = 1, nva
c
c     <ef||mn>
c
                  r2aa(e,f,m,n) = iaa(e+noa,f+noa,m,n)
c
c     Ffa*t(e,a,m,n)-Fea*t(f,a,m,n)
c
                  do a = 1, nva
                     r2aa(e,f,m,n) = r2aa(e,f,m,n)
     $                    + fa(f+noa,a+noa)*t2aa(e,a,m,n) 
     $                    - fa(e+noa,a+noa)*t2aa(f,a,m,n) 
                  end do
c
c     -Fjn*t(e,f,m,j)+Fjm*t(e,f,n,j)
c
                  do j = 1, noa
                     r2aa(e,f,m,n) = r2aa(e,f,m,n)
     $                    - fa(j,n)*t2aa(e,f,m,j) 
     $                    + fa(j,m)*t2aa(e,f,n,j) 
                  end do
c
c     (1-(ef))(1-(mn)) <if||an>*t(e,a,m,i)
c
                  do i = 1, noa
                     do a = 1, nva
                        r2aa(e,f,m,n) = r2aa(e,f,m,n)
     $                       + iaa(i,f+noa,a+noa,n)*t2aa(e,a,m,i)
     $                       - iaa(i,e+noa,a+noa,n)*t2aa(f,a,m,i)
     $                       - iaa(i,f+noa,a+noa,m)*t2aa(e,a,n,i)
     $                       + iaa(i,e+noa,a+noa,m)*t2aa(f,a,n,i)
                     end do
                  end do
c
c     (1-(ef))(1-(mn)) <fi*|na*>*t(e,a*,m,i*)
c
                  do i = 1, nob
                     do a = 1, nvb
                        r2aa(e,f,m,n) = r2aa(e,f,m,n)
     $                       + iab(f+noa,i,n,a+nob)*t2ab(e,a,m,i)
     $                       - iab(e+noa,i,n,a+nob)*t2ab(f,a,m,i)
     $                       - iab(f+noa,i,m,a+nob)*t2ab(e,a,n,i)
     $                       + iab(e+noa,i,m,a+nob)*t2ab(f,a,n,i)
                     end do
                  end do
c
c     (1/2)<ij||mn>*t(e,f,i,j)
c
                  do j = 1, noa
                     do i = 1, noa
                        r2aa(e,f,m,n) = r2aa(e,f,m,n) 
     $                       + 0.5d0*iaa(i,j,m,n)*t2aa(e,f,i,j)
                     end do
                  end do
c
c     (1/2)<ef||ab>*t(a,b,m,n)
c
                  do b = 1, nva
                     do a = 1, nva
                        r2aa(e,f,m,n) = r2aa(e,f,m,n)
     $                       + 0.5d0*iaa(e+noa,f+noa,a+noa,b+noa)*
     $                       t2aa(a,b,m,n)
                     end do
                  end do
c     
c     (1/2)*<ij||ab>[+(1-(ef))(1-(mn))t(f,b,j,n)t(e,a,i,m)
c     .              -(1-(ef))t(f,b,m,n)*t(e,a,i,j)
c     .              -(1-(mn))t(e,f,m,i)*t(a,b,n,j)
c     .              +(1/2)t(e,f,i,j)*t(a,b,m,n)]
c
                  do j = 1, noa
                     do i = 1, noa
                        do b = 1, nva
                           do a = 1, nva
                              r2aa(e,f,m,n) = r2aa(e,f,m,n) + 
     $                             0.5d0*iaa(i,j,a+noa,b+noa)*(
                              
     $                             + t2aa(f,b,j,n)*t2aa(e,a,i,m)
     $                             - t2aa(e,b,j,n)*t2aa(f,a,i,m)
     $                             - t2aa(f,b,j,m)*t2aa(e,a,i,n)
     $                             + t2aa(e,b,j,m)*t2aa(f,a,i,n)
                              
     $                             - t2aa(f,b,m,n)*t2aa(e,a,i,j)
     $                             + t2aa(e,b,m,n)*t2aa(f,a,i,j)
                              
     $                             - t2aa(e,f,m,i)*t2aa(a,b,n,j)
     $                             + t2aa(e,f,n,i)*t2aa(a,b,m,j)
                              
     $                             + 0.5d0*t2aa(e,f,i,j)*t2aa(a,b,m,n))
                           end do
                        end do
                     end do
                  end do
c
c     (1/2)<i*j*||a*b*>(1-(ef))(1-(mn))t(f,b*,n,j*)t(e,a*,m,i*)
c
                  do j = 1, nob
                     do i = 1, nob
                        do b = 1, nvb
                           do a = 1, nvb
                              r2aa(e,f,m,n) = r2aa(e,f,m,n) + 
     $                             0.5d0*ibb(i,j,a+nob,b+nob)*(0.0d0
     $                             + t2ab(f,b,n,j)*t2ab(e,a,m,i)
     $                             - t2ab(e,b,n,j)*t2ab(f,a,m,i)
     $                             - t2ab(f,b,m,j)*t2ab(e,a,n,i)
     $                             + t2ab(e,b,m,j)*t2ab(f,a,n,i))
                           end do
                        end do
                     end do
                  end do
c
c     <ij*||ab*>[(1-(ef))t(f,a,m,n)t(e,b*,i,j*)
c     .         -(1-(mn))t(e,f,m,i)*t(a,b*,n,j*)
c     .         -(1-(ef))(1-(mn))t(f,b*,n,j*)t(e,a,i,m)]
c
                  do j = 1, nob
                     do i = 1, noa
                        do b = 1, nvb
                           do a = 1, nva
                              r2aa(e,f,m,n) = r2aa(e,f,m,n) + 
     $                             iab(i,j,a+noa,b+nob)*(0.0d0

     $                             + t2aa(f,a,m,n)*t2ab(e,b,i,j)
     $                             - t2aa(e,a,m,n)*t2ab(f,b,i,j)

     $                             - t2aa(e,f,m,i)*t2ab(a,b,n,j)
     $                             + t2aa(e,f,n,i)*t2ab(a,b,m,j)

     $                             - t2ab(f,b,n,j)*t2aa(e,a,i,m)
     $                             + t2ab(e,b,n,j)*t2aa(f,a,i,m)
     $                             + t2ab(f,b,m,j)*t2aa(e,a,i,n)
     $                             - t2ab(e,b,m,j)*t2aa(f,a,i,n))
                           end do
                        end do
                     end do
                  end do
c
               end do
            end do
         end do
      end do
c                  
      end
      subroutine uccsdtest_product_mixed(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     fa, fb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     r1a, r1b,
     $     r2ab)
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision fa(nmo,nmo), fb(nmo,nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
      double precision r1a(nva, noa)
      double precision r1b(nvb, nob)
      double precision r2ab(nva, nvb, noa, nob)
c
      integer i, j, m, n
      integer a, b, e, f
c
c     The formulae in comments denote beta spin with an asterisk
c
      do n = 1, nob
         do m = 1, noa
            do f = 1, nvb
               do e = 1, nva
c
c     <ef||mn>
c
                  r2ab(e,f,m,n) = iab(e+noa,f+nob,m,n)
c
c     .  F(f*,a*)t(e,a*,m,n*) + F(e,a)t(a,f*,m,n*) 
c     . -F(i*,n*)t(e,f*,m,i*) - F(i,m)t(e,f*,i,n*)
c
                  do a = 1, nvb
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                    + fb(f+nob,a+nob)*t2ab(e,a,m,n) 
                  end do
                  do a = 1, nva
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                    + fa(e+noa,a+noa)*t2ab(a,f,m,n) 
                  end do
                  do i = 1, nob
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                 - fb(i,n)*t2ab(e,f,m,i) 
                  end do
                  do i = 1, noa
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                 - fa(i,m)*t2ab(e,f,i,n)
                  end do
c
c     <i*f*||a*n*>t(e,a*,m,i*)
c
                  do i = 1, nob
                     do a = 1, nvb
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       + ibb(i,f+nob,a+nob,n)*t2ab(e,a,m,i)
                     end do
                  end do
c
c     <if*|an*>t(e,a,m,i)
c
                  do i = 1, noa
                     do a = 1, nva
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       + iab(i,f+nob,a+noa,n)*t2aa(e,a,m,i)
                     end do
                  end do
c
c     -<ei*|an*>t(a,f*,m,i*)
c
                  do i = 1, nob
                     do a = 1, nva
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       - iab(e+noa,i,a+noa,n)*t2ab(a,f,m,i)
                     end do
                  end do
c
c     -<if*|ma*>t(e,a*,i,n*)
c
                  do i = 1, noa
                     do a = 1, nvb
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       - iab(i,f+nob,m,a+nob)*t2ab(e,a,i,n)
                     end do
                  end do
c
c     +<ie||am>t(a,f,i,n)
c
                  do i = 1, noa
                     do a = 1, nva
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       + iaa(i,e+noa,a+noa,m)*t2ab(a,f,i,n)
                     end do
                  end do
c
c     +<ei*|ma*>t(f*,a*,n*,i*)
c
                  do i = 1, nob
                     do a = 1, nvb
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       + iab(e+noa,i,m,a+nob)*t2bb(f,a,n,i)
                     end do
                  end do
c
c     <ij*|mn*>t(e,f*,i,j*)
c
                  do i = 1, noa
                     do j = 1, nob
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       + iab(i,j,m,n)*t2ab(e,f,i,j)
                     end do
                  end do
c
c     <ef*|ab*>t(a,b*,m,n*)
c
                  do b = 1, nvb
                     do a = 1, nva
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                       +iab(e+noa,f+nob,a+noa,b+nob)*t2ab(a,b,m,n)
                     end do
                  end do
c
c     <ij||ab>[-t(b,f*,j,n*)t(e,a,i,m) 
c     .   +(1/2)t(b,f*,m,n*)t(e,a,i,j)
c     .   -(1/2)t(e,f*,i,n*)t(a,b,m,j)
c
                  do j = 1, noa
                     do i = 1, noa
                        do b = 1, nva
                           do a = 1, nva
                              r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                             + iaa(i,j,a+noa,b+noa)*(0.0d0
                              
     $                             - t2ab(b,f,j,n)*t2aa(e,a,i,m)
     $                             + 0.5d0*t2ab(b,f,m,n)*t2aa(e,a,i,j)
     $                             - 0.5d0*t2ab(e,f,i,n)*t2aa(a,b,m,j))
                           end do
                        end do
                     end do
                  end do
c
c     <i*j*||a*b*>[-t(e,a*,m,i*)t(f*,b*,j*,n*)
c     .      +(1/2)t(e,b*,m,n*)t(f*,a*,i*,j*)
c     .      -(1/2)t(e,f*,m,i*)t(a*,b*,n*,j*)]
c
                  do b = 1, nvb
                     do a = 1, nvb
                        do j = 1, nob
                           do i = 1, nob
                              r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                             + ibb(i,j,a+nob,b+nob)*(0.0d0
     $                             - t2ab(e,a,m,i)*t2bb(f,b,j,n)
     $                             + t2ab(e,b,m,n)*t2bb(f,a,i,j)*0.5d0
     $                             - t2ab(e,f,m,i)*t2bb(a,b,n,j)*0.5d0)
                           end do
                        end do
                     end do
                  end do
c
c     <ij*|ab*>[-t(a,f*,m,n*)t(e,b*,i,j*)
c     .         -t(e,f*,m,j*)t(a,b*,i,n*)
c     .         -t(e,b*,m,n*)t(a,f*,i,j*)
c     .         -t(e,f*,i,n*)t(a,b*,m,j*)
c     .         +t(e,f*,i,j*)t(a,b*,m,n*)
c     .         +t(f*,b*,j*,n*)t(e,a,i,m)
c     .         +t(a,f*,i,n*)t(e,b*,m,j*)
c     .         +t(a,f*,m,j*)t(e,b*,i,n*)
c
                  do b = 1, nvb
                     do a = 1, nva
                        do j = 1, nob
                           do i = 1, noa
                              r2ab(e,f,m,n) = r2ab(e,f,m,n) 
     $                             + iab(i,j,a+noa,b+nob)*(0.0d0
     $                             - t2ab(a,f,m,n)*t2ab(e,b,i,j)
     $                             - t2ab(e,f,m,j)*t2ab(a,b,i,n)
     $                             - t2ab(e,b,m,n)*t2ab(a,f,i,j)
     $                             - t2ab(e,f,i,n)*t2ab(a,b,m,j)
     $                             + t2ab(e,f,i,j)*t2ab(a,b,m,n)
     $                             + t2bb(f,b,j,n)*t2aa(e,a,i,m)
     $                             + t2ab(a,f,i,n)*t2ab(e,b,m,j)
     $                             + t2ab(a,f,m,j)*t2ab(e,b,i,n))
                           end do
                        end do
                     end do
                  end do
c
               end do
            end do
         end do
      end do
c
      end
      subroutine uccsdtest_mp2(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     ea, eb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab)
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision ea(nmo), eb(nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
c
c     t1 = 0
c     t2 = <ef||mn>/(em+en-ee-ef)
c     
      integer m, n, e, f
      double precision eaa, ebb, eab, energy
c
      call dfill(noa*nva, 0.0d0, t1a, 1)
      call dfill(nob*nvb, 0.0d0, t1b, 1)
c
c     Pure alpha
c
      eaa = 0.0d0
      do n = 1, noa
         do m = 1, noa
            do f = 1, nva
               do e = 1, nva
                  t2aa(e,f,m,n) = iaa(noa+e,noa+f,m,n) / 
     $                 (ea(m)+ea(n)-ea(noa+e)-ea(noa+f))
                  eaa = eaa + t2aa(e,f,m,n)*iaa(m,n,noa+e,noa+f)
               end do
            end do
         end do
      end do
c
c     Pure beta
c
      ebb = 0.0d0
      do n = 1, nob
         do m = 1, nob
            do f = 1, nvb
               do e = 1, nvb
                  t2bb(e,f,m,n) = ibb(nob+e,nob+f,m,n) / 
     $                 (eb(m)+eb(n)-eb(nob+e)-eb(nob+f))
                  ebb = ebb + t2bb(e,f,m,n)*ibb(m,n,nob+e,nob+f)
               end do
            end do
         end do
      end do
c
c     Mixed
c
      eab = 0.0d0
      do n = 1, nob
         do m = 1, noa
            do f = 1, nvb
               do e = 1, nva
                  t2ab(e,f,m,n) = iab(noa+e,nob+f,m,n) / 
     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
                  eab = eab + t2ab(e,f,m,n)*iab(m,n,noa+e,nob+f)
               end do
            end do
         end do
      end do
c
      energy = eaa/4.0d0 + ebb/4.0d0 + eab
c
      write(6,1) energy
 1    format(//,'  The UMP2 energy is ', f20.8)
      write(6,*) eaa/4.0d0, ebb/4.0d0, eab
c
      end
      subroutine uccsdtest_transpose_tab(noa, nob, nva, nvb, t2ab, t2ba)
      implicit none
      integer noa, nob, nva, nvb
      double precision t2ab(nva, nvb, noa, nob)
      double precision t2ba(nvb, nva, nob, noa)
c
      integer i, j, a, b
c
      do j = 1, nob
         do i = 1, noa
            do b = 1, nvb
               do a = 1, nva
                  t2ba(b,a,j,i) = t2ab(a,b,i,j)
               end do
            end do
         end do
      end do
c
      end
      subroutine uccsdtest_update(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     fa, fb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     r1a, r1b,
     $     r2aa,r2bb,r2ab, ea, eb)
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision fa(nmo,nmo), fb(nmo,nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
      double precision r1a(nva, noa)
      double precision r1b(nvb, nob)
      double precision r2aa(nva, nva, noa, noa)
      double precision r2bb(nvb, nvb, nob, nob)
      double precision r2ab(nva, nvb, noa, nob)
      double precision ea(nmo), eb(nmo)
c
c     t1 = t1 + R(e,m)/(em-ee)
c     t2 = t2 + R(e,f,m,n)/(em+en-ee-ef)
c     
      integer m, n, e, f
      double precision eaa, ebb, eab, energy
      double precision delta, r1norm, r2norm 
c
      r1norm = 0.0d0
      r2norm = 0.0d0
c
*      call dfill(noa*nva, 0.0d0, r1a, 1)
*      call dfill(nob*nvb, 0.0d0, r1b, 1)
c
c     Singles
c
      do m = 1, noa
         do e = 1, nva
            delta = r1a(e,m) / (ea(m) - ea(e+noa))
            r1norm = r1norm + delta**2
            t1a(e,m) = t1a(e,m) + delta
         end do
      end do
      do m = 1, nob
         do e = 1, nvb
            delta = r1b(e,m) / (eb(m) - eb(e+nob))
            r1norm = r1norm + delta**2
            t1b(e,m) = t1b(e,m) + delta
         end do
      end do
c
c     Pure alpha
c
      eaa = 0.0d0
      do n = 1, noa
         do m = 1, noa
            do f = 1, nva
               do e = 1, nva
                  delta = r2aa(e,f,m,n) / 
     $                 (ea(m)+ea(n)-ea(noa+e)-ea(noa+f))
                  r2norm = r2norm + delta*delta
                  t2aa(e,f,m,n) = t2aa(e,f,m,n) + delta
                  eaa = eaa + (
     $                 t2aa(e,f,m,n) + 
     $                 t1a(e,m)*t1a(f,n) - t1a(e,n)*t1a(f,m))*
     $                 iaa(m,n,noa+e,noa+f)
               end do
            end do
         end do
      end do
c
c     Pure beta
c
      ebb = 0.0d0
      do n = 1, nob
         do m = 1, nob
            do f = 1, nvb
               do e = 1, nvb
                  delta = r2bb(e,f,m,n) /
     $                 (eb(m)+eb(n)-eb(nob+e)-eb(nob+f))
                  r2norm = r2norm + delta*delta
                  t2bb(e,f,m,n) = t2bb(e,f,m,n) + delta
                  ebb = ebb + (
     $                 t2bb(e,f,m,n) + 
     $                 t1b(e,m)*t1b(f,n) - t1b(e,n)*t1b(f,m))*
     $                 ibb(m,n,nob+e,nob+f)
               end do
            end do
         end do
      end do
c
c     Mixed
c
      eab = 0.0d0
      do n = 1, nob
         do m = 1, noa
            do f = 1, nvb
               do e = 1, nva
                  delta = r2ab(e,f,m,n) /
     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
*                  write(6,*) ' denom ', 
*     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
                  r2norm = r2norm + delta*delta
                  t2ab(e,f,m,n) = t2ab(e,f,m,n) + delta
                  eab = eab + (
     $                 t2ab(e,f,m,n) + 
     $                 t1a(e,m)*t1b(f,n))*
     $                 iab(m,n,noa+e,nob+f)
               end do
            end do
         end do
      end do
c
      energy = eaa/4.0d0 + ebb/4.0d0 + eab
c
      write(6,1) energy, sqrt(r2norm), sqrt(r1norm)
 1    format('  Energy = ', f20.8,'  residual norm = ',1p,2d9.2)
c
      end
      subroutine uccsdtest_fock(nmo, noa, nob, ha, iaa, iab, fa)
      implicit none
      integer nmo, noa, nob
      double precision ha(nmo,nmo)
      double precision iaa(nmo,nmo,nmo,nmo), iab(nmo,nmo,nmo,nmo)
      double precision fa(nmo,nmo)
c
c     fij = hij + <ik||jk> = hij + <ik||jk>(same spin) + <ik|jk>(diff spin)
c
      integer i, j, k
c
      call dfill(nmo*nmo, 0.0d0, fa, 1)
c
      do j = 1, nmo
         do i = 1, nmo
            fa(i,j) = ha(i,j)
            do k = 1, noa
               fa(i,j) = fa(i,j) + iaa(i,k,j,k)
            end do
            do k = 1, nob
               fa(i,j) = fa(i,j) + iab(i,k,j,k)
            end do
         end do
      end do
c
c$$$      write(6,*) ' FOCK'
c$$$      call output(fa, 1, nmo, 1, nmo, nmo, nmo, 1)
c
      end
      subroutine uccsdtest_intermediates(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     fa, fb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     z1, z2, z3, z4, z5, z6, z7, z8, x1, x2, x3)
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision fa(nmo,nmo), fb(nmo,nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
c
      double precision z1(noa, nva, noa, nva)
      double precision z2(noa, nva, noa, nva)
      double precision z3(nob, nvb, noa, nva)
      double precision z4(nob, nvb, nob, nvb)
      double precision z5(nob, nvb, nob, nvb)
      double precision z6(noa, nva, nob, nvb)
      double precision z7(noa, nva, nob, nvb)
      double precision z8(noa, nvb, noa, nvb)
      double precision x1(noa, noa, noa, noa)
      double precision x2(nob, nob, nob, nob)
      double precision x3(noa, nob, noa, nob)
c
      integer i, j, a, b, f, m, n
c
      call dfill(noa*nva*noa*nva, 0.0d0, z1, 1)
      call dfill(noa*nva*noa*nva, 0.0d0, z2, 1)
      call dfill(nob*nvb*noa*nva, 0.0d0, z3, 1)
      call dfill(nob*nvb*nob*nvb, 0.0d0, z4, 1)
      call dfill(nob*nvb*nob*nvb, 0.0d0, z5, 1)
      call dfill(noa*nva*nob*nvb, 0.0d0, z6, 1)
      call dfill(noa*nva*nob*nvb, 0.0d0, z7, 1)
      call dfill(noa*nvb*noa*nvb, 0.0d0, z8, 1)
      call dfill(noa*noa*noa*noa, 0.0d0, x1, 1)
      call dfill(nob*nob*nob*nob, 0.0d0, x2, 1)
      call dfill(noa*nob*noa*nob, 0.0d0, x3, 1)
c
      do i = 1, noa
         do a = 1, nva
            do n = 1, noa
               do f = 1, nva
                  do j = 1, noa
                     do b = 1, nva
                        z1(i,a,n,f) = z1(i,a,n,f) + 
     $                       iaa(i,j,a+noa,b+noa)*t2aa(f,b,n,j)
                     end do
                  end do
                  do j = 1, nob
                     do b = 1, nvb
                        z2(i,a,n,f) = z2(i,a,n,f) + 
     $                       iab(i,j,a+noa,b+nob)*t2ab(f,b,n,j)
                     enddo
                  enddo
               end do
            end do
         end do
      end do
c
      do i = 1, nob
         do a = 1, nvb
            do n = 1, noa
               do f = 1, nva
                  do b = 1, nvb
                     do j = 1, nob
                        z3(i,a,n,f) = z3(i,a,n,f) + 
     $                       ibb(i,j,a+nob,b+nob)*t2ab(f,b,n,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, nob
         do a = 1, nvb
            do n = 1, nob
               do f = 1, nvb
                  do j = 1, nob
                     do b = 1, nvb
                        z4(i,a,n,f) = z4(i,a,n,f) + 
     $                       ibb(i,j,a+nob,b+nob)*t2bb(f,b,n,j)
                     end do
                  end do
                  do j = 1, noa
                     do b = 1, nva
                        z5(i,a,n,f) = z5(i,a,n,f) + 
     $                       iab(j,i,b+noa,a+nob)*t2ab(b,f,j,n)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, noa
         do a = 1, nva
            do n = 1, nob
               do f = 1, nvb
                  do j = 1, noa
                     do b = 1, nva
                        z6(i,a,n,f) = z6(i,a,n,f) + 
     $                       iaa(i,j,a+noa,b+noa)*t2ab(b,f,j,n)
                     end do
                  end do
                  do j = 1, nob
                     do b = 1, nvb
                        z7(i,a,n,f) = z7(i,a,n,f) + 
     $                       iab(i,j,a+noa,b+nob)*t2bb(f,b,n,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, noa
         do a = 1, nvb
            do m = 1, noa
               do f = 1, nvb
                  do j = 1, nob
                     do b = 1, nva
                        z8(i,a,m,f) = z8(i,a,m,f) + 
     $                       iab(i,j,b+noa,a+nob)*t2ab(b,f,m,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, noa
         do j = 1, noa
            do m = 1, noa
               do n = 1, noa
                  do a = 1, nva
                     do b = 1, nva
                        x1(i,j,m,n) = x1(i,j,m,n) + 
     $                       iaa(i,j,a+noa,b+noa)*t2aa(a,b,m,n)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, nob
         do j = 1, nob
            do m = 1, nob
               do n = 1, nob
                  do a = 1, nvb
                     do b = 1, nvb
                        x2(i,j,m,n) = x2(i,j,m,n) + 
     $                       ibb(i,j,a+nob,b+nob)*t2bb(a,b,m,n)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, noa
         do j = 1, nob
            do m = 1, noa
               do n = 1, nob
                  do a = 1, nva
                     do b = 1, nvb
                        x3(i,j,m,n) = x3(i,j,m,n) + 
     $                       iab(i,j,a+noa,b+nob)*t2ab(a,b,m,n)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      end
      subroutine uccsdtest_new_product(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     fa, fb,
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     z1, z2, z3, z4, z5, z6, z7, z8, x1, x2, x3,
     $     r1a, r1b,
     $     r2aa,r2bb,r2ab)
      implicit none
#include "mafdecls.fh"
      integer nmo, noa, nob, nva, nvb
      double precision fa(nmo,nmo), fb(nmo,nmo)
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
c
      double precision z1(noa, nva, noa, nva)
      double precision z2(noa, nva, noa, nva)
      double precision z3(nob, nvb, noa, nva)
      double precision z4(nob, nvb, nob, nvb)
      double precision z5(nob, nvb, nob, nvb)
      double precision z6(noa, nva, nob, nvb)
      double precision z7(noa, nva, nob, nvb)
      double precision z8(noa, nvb, noa, nvb)
      double precision x1(noa, noa, noa, noa)
      double precision x2(nob, nob, nob, nob)
      double precision x3(noa, nob, noa, nob)
c
      double precision r1a(nva, noa)
      double precision r1b(nvb, nob)
      double precision r2aa(nva, nva, noa, noa)
      double precision r2bb(nvb, nvb, nob, nob)
      double precision r2ab(nva, nvb, noa, nob)
c
      integer i, j, a, b, e, f, m, n
      double precision tmp
c
c$$$      call jan_debug_print('z1', z1, noa, nva, noa, nva)
c$$$      call jan_debug_print('z2', z2, noa, nva, noa, nva)
c$$$      call jan_debug_print('z3', z3, nob, nvb, noa, nva)
c$$$      call jan_debug_print('z4', z4, nob, nvb, nob, nvb)
c$$$      call jan_debug_print('z5', z5, nob, nvb, nob, nvb)
c$$$      call jan_debug_print('z6', z6, noa, nva, nob, nvb)
c$$$      call jan_debug_print('z7', z7, noa, nva, nob, nvb)
c$$$      call jan_debug_print('z8', z8, noa, nvb, noa, nvb)
c$$$      call jan_debug_print('x1', x1, noa, noa, noa, noa)
c$$$      call jan_debug_print('x2', x2, nob, nob, nob, nob)
c$$$      call jan_debug_print('x3', x3, noa, nob, noa, nob)
c
c     R1 alpha
c
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('x',1)
      do m = 1, noa
         do e = 1, nva
c
            r1a(e,m) = fa(e+noa,m)
c
            do j = 1, noa
               do b = 1, nva
                  r1a(e,m) = r1a(e,m) + fa(j,b+noa)*t2aa(e,b,m,j)
               end do
            end do
            do j = 1, nob
               do b = 1, nvb
                  r1a(e,m) = r1a(e,m) + fb(j,b+nob)*t2ab(e,b,m,j)
               end do
            end do
c
            do j = 1, noa
               do i = 1, noa
                  do b = 1, nva
                     r1a(e,m) = r1a(e,m) 
     $                    - 0.5d0*iaa(i,j,m,b+noa)*t2aa(e,b,i,j)
                  end do
               end do
            end do
c
            do j = 1, nob
               do i = 1, noa
                  do b = 1, nvb
                     r1a(e,m) = r1a(e,m)
     $                    - iab(i,j,m,b+nob)*t2ab(e,b,i,j)
                  end do
               end do
            end do
c
            do i = 1, noa
               do b = 1, nva
                  do a = 1, nva
                     r1a(e,m) = r1a(e,m)
     $                    + 0.5d0*iaa(i,e+noa,a+noa,b+noa)*t2aa(a,b,i,m)
                  end do
               end do
            end do
c
            do i = 1, nob
               do b = 1, nva
                  do a = 1, nvb
                     r1a(e,m) = r1a(e,m)
     $                    + iab(e+noa,i,b+noa,a+nob)*t2ab(b,a,m,i)
                  end do
               end do
            end do
c
         end do
      end do
c
c     R1 beta
c
      do m = 1, nob
         do e = 1, nvb
c
            r1b(e,m) = fb(e+nob,m)
c
            do j = 1, nob
               do b = 1, nvb
                  r1b(e,m) = r1b(e,m) + fb(j,b+nob)*t2bb(e,b,m,j)
               end do
            end do
            do j = 1, noa
               do b = 1, nva
                  r1b(e,m) = r1b(e,m) + fa(j,b+noa)*t2ab(b,e,j,m) !(e,b,m,j)
               end do
            end do
c
            do j = 1, nob
               do i = 1, nob
                  do b = 1, nvb
                     r1b(e,m) = r1b(e,m) 
     $                    - 0.5d0*ibb(i,j,m,b+nob)*t2bb(e,b,i,j)
                  end do
               end do
            end do
c
            do j = 1, noa
               do i = 1, nob
                  do b = 1, nva
                     r1b(e,m) = r1b(e,m)
     $                    - iab(j,i,b+noa,m)*t2ab(b,e,j,i)
                  end do
               end do
            end do
c
            do i = 1, nob
               do b = 1, nvb
                  do a = 1, nvb
                     r1b(e,m) = r1b(e,m)
     $                    + 0.5d0*ibb(i,e+nob,a+nob,b+nob)*t2bb(a,b,i,m)
                  end do
               end do
            end do
c
            do i = 1, noa
               do b = 1, nvb
                  do a = 1, nva
                     r1b(e,m) = r1b(e,m)
     $                    + iab(i,e+nob,a+noa,b+nob)*t2ab(a,b,i,m)
                  end do
               end do
            end do
c
         end do
      end do
c
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('x',2)
c
c     Pure spin alpha
c
      do m = 1, noa
         do n = 1, noa
            do e = 1, nva
               do f = 1, nva
c
                  do i = 1, noa ! A
                     do j = 1, noa
                        r2aa(e,f,m,n) = r2aa(e,f,m,n) + 
     $                       0.5d0*t2aa(e,f,i,j)*
     $                       (iaa(i,j,m,n) + 0.5d0*x1(i,j,m,n))
                     end do
                  end do
c
                  do a = 1, nva ! B
                     do b = 1, nva
                        r2aa(e,f,m,n) = r2aa(e,f,m,n) + 0.5d0*
     $                       iaa(e+noa,f+noa,a+noa,b+noa)*t2aa(a,b,m,n)
                     end do
                  end do
c
                  do i = 1, noa ! C1
                     do a = 1, nva
                        tmp = t2aa(e,a,m,i) * (iaa(i,f+noa,a+noa,n) + 
     $                       0.5d0*z1(i,a,n,f) + z2(i,a,n,f))
                        r2aa(e,f,m,n) = r2aa(e,f,m,n) + tmp
                        r2aa(f,e,m,n) = r2aa(f,e,m,n) - tmp
                        r2aa(e,f,n,m) = r2aa(e,f,n,m) - tmp
                        r2aa(f,e,n,m) = r2aa(f,e,n,m) + tmp
                     end do
                  end do
c
                  do i = 1, nob ! C2
                     do a = 1, nvb
                        tmp = t2ab(e,a,m,i) * (iab(f+noa,i,n,a+nob) + 
     $                       0.5d0*z3(i,a,n,f))
                        r2aa(e,f,m,n) = r2aa(e,f,m,n) + tmp
                        r2aa(f,e,m,n) = r2aa(f,e,m,n) - tmp
                        r2aa(e,f,n,m) = r2aa(e,f,n,m) - tmp
                        r2aa(f,e,n,m) = r2aa(f,e,n,m) + tmp
                     end do
                  end do
c
                  do a = 1, nva ! E1
                     tmp = 0.0d0
                     do i = 1, noa
                        tmp = tmp - 0.5d0*z1(i,a,i,f) - z2(i,a,i,f)
                     enddo
                     tmp = t2aa(e,a,m,n) * (fa(f+noa,a+noa) + tmp)
                     r2aa(e,f,m,n) = r2aa(e,f,m,n) + tmp
                     r2aa(f,e,m,n) = r2aa(f,e,m,n) - tmp
                  end do
c
                  do i = 1, noa ! E2
                     tmp = 0.0d0
                     do a = 1, nva
                        tmp = tmp + 0.5d0*z1(i,a,n,a) + z2(i,a,n,a)
                     end do
                     tmp = t2aa(e,f,m,i) * (fa(i,n) + tmp)
                     r2aa(e,f,m,n) = r2aa(e,f,m,n) - tmp
                     r2aa(e,f,n,m) = r2aa(e,f,n,m) + tmp
                  enddo
c
                  r2aa(e,f,m,n) = r2aa(e,f,m,n) + iaa(e+noa,f+noa,m,n) ! F
c
               end do
            end do
         end do
      end do
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('x',4)
c
c     Pure spin beta
c
c     This is just the alpha code with spins flipped, x1 -> x2, and
c     z1/2/3 -> z4/5/6.  For mixed spin quantities have to transpose
c     the indices so that end up with beta spin in the correct position.
c
      do m = 1, nob
         do n = 1, nob
            do e = 1, nvb
               do f = 1, nvb
c
                  do i = 1, nob ! A
                     do j = 1, nob
                        r2bb(e,f,m,n) = r2bb(e,f,m,n) + 
     $                       0.5d0*t2bb(e,f,i,j)*
     $                       (ibb(i,j,m,n) + 0.5d0*x2(i,j,m,n))
                     end do
                  end do
c
                  do a = 1, nvb ! B
                     do b = 1, nvb
                        r2bb(e,f,m,n) = r2bb(e,f,m,n) + 0.5d0*
     $                       ibb(e+nob,f+nob,a+nob,b+nob)*t2bb(a,b,m,n)
                     end do
                  end do
c
                  do i = 1, nob ! C1
                     do a = 1, nvb
                        tmp = t2bb(e,a,m,i) * (ibb(i,f+nob,a+nob,n) + 
     $                       0.5d0*z4(i,a,n,f) + z5(i,a,n,f))
                        r2bb(e,f,m,n) = r2bb(e,f,m,n) + tmp
                        r2bb(f,e,m,n) = r2bb(f,e,m,n) - tmp
                        r2bb(e,f,n,m) = r2bb(e,f,n,m) - tmp
                        r2bb(f,e,n,m) = r2bb(f,e,n,m) + tmp
                     end do
                  end do
c
                  do i = 1, noa ! C2
                     do a = 1, nva
                        tmp = t2ab(a,e,i,m) * (iab(i,f+nob,a+noa,n) + 
     $                       0.5d0*z6(i,a,n,f))
                        r2bb(e,f,m,n) = r2bb(e,f,m,n) + tmp
                        r2bb(f,e,m,n) = r2bb(f,e,m,n) - tmp
                        r2bb(e,f,n,m) = r2bb(e,f,n,m) - tmp
                        r2bb(f,e,n,m) = r2bb(f,e,n,m) + tmp
                     end do
                  end do
c
                  do a = 1, nvb ! E1
                     tmp = 0.0d0
                     do i = 1, nob
                        tmp = tmp - 0.5d0*z4(i,a,i,f) - z5(i,a,i,f)
                     enddo
                     tmp = t2bb(e,a,m,n) * (fb(f+nob,a+nob) + tmp)
                     r2bb(e,f,m,n) = r2bb(e,f,m,n) + tmp
                     r2bb(f,e,m,n) = r2bb(f,e,m,n) - tmp
                  end do
c
                  do i = 1, nob ! E2
                     tmp = 0.0d0
                     do a = 1, nvb
                        tmp = tmp + 0.5d0*z4(i,a,n,a) + z5(i,a,n,a)
                     end do
                     tmp = t2bb(e,f,m,i) * (fb(i,n) + tmp)
                     r2bb(e,f,m,n) = r2bb(e,f,m,n) - tmp
                     r2bb(e,f,n,m) = r2bb(e,f,n,m) + tmp
                  enddo
c
                  r2bb(e,f,m,n) = r2bb(e,f,m,n) + ibb(e+nob,f+nob,m,n) ! F
c
               end do
            end do
         end do
      end do
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('x',5)
c
c     Mixed spin
c
      do m = 1, noa
         do n = 1, nob
            do e = 1, nva
               do f = 1, nvb
c
                  do i = 1, noa ! A
                     do j = 1, nob
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) + 
     $                       t2ab(e,f,i,j)*
     $                       (iab(i,j,m,n) + x3(i,j,m,n))
                     end do
                  end do
c
                  do a = 1, nva ! B
                     do b = 1, nvb
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) + 
     $                       iab(e+noa,f+nob,a+noa,b+nob)*t2ab(a,b,m,n)
                     end do
                  end do
c
                  do a = 1, nva ! C1
                     do i = 1, noa
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) - t2aa(e,a,i,m)*
     $                       (iab(i,f+nob,a+noa,n) + 
     $                       z6(i,a,n,f) + z7(i,a,n,f))
                     end do
                  end do
c
                  do a = 1, nvb ! C2
                     do i = 1, nob
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) + t2ab(e,a,m,i)*
     $                       (ibb(i,f+nob,a+nob,n) + z4(i,a,n,f))
                     end do
                  end do
c
                  do a = 1, nvb ! C3
                     do i = 1, noa
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) - t2ab(e,a,i,n)*
     $                       (iab(i,f+nob,m,a+nob) - z8(i,a,m,f))
                     end do
                  end do
c
                  do a = 1, nva ! C4
                     do i = 1, noa
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) + t2ab(a,f,i,n)*
     $                       (iaa(i,e+noa,a+noa,m) + z2(i,a,m,e))
                     end do
                  end do
c
                  do a = 1, nvb ! C5
                     do i = 1, nob
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) + 
     $                       iab(e+noa,i,m,a+nob)*t2bb (f,a,n,i)
                     end do
                  end do
c
                  do a = 1, nva ! C5
                     do i = 1, nob
                        r2ab(e,f,m,n) = r2ab(e,f,m,n) -
     $                       iab(e+noa,i,a+noa,n)*t2ab(a,f,m,i)
                     end do
                  end do
c
                  do a = 1, nva ! E1
                     tmp = 0.0d0
                     do i = 1, noa
                        tmp = tmp - z2(i,a,i,e) - 0.5d0*z1(i,a,i,e)
                     end do
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) + t2ab(a,f,m,n)*
     $                       (fa(e+noa,a+noa) + tmp)
                  end do
c
                  do a = 1, nvb ! E2
                     tmp = 0.0d0
                     do i = 1, noa
                        tmp = tmp - z8(i,a,i,f)
                     end do
                     do i = 1, nob
                        tmp = tmp - 0.5d0*z4(i,a,i,f)
                     end do
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) + t2ab(e,a,m,n)*
     $                    (fb(f+nob,a+nob) + tmp)
                  end do
c
                  do i = 1, nob ! E3
                     tmp = 0.0d0
                     do a = 1, nvb
                        tmp = tmp + z5(i,a,n,a) + 0.5d0*z4(i,a,n,a)
                     end do
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) - t2ab(e,f,m,i)*
     $                    (fb(i,n) + tmp)
                  end do
c
                  do i = 1, noa ! E4
                     tmp = 0.0d0
                     do a = 1, nva
                        tmp = tmp + z2(i,a,m,a) + 0.5d0*z1(i,a,m,a)
                     end do
                     r2ab(e,f,m,n) = r2ab(e,f,m,n) - t2ab(e,f,i,n)*
     $                    (fa(i,m) + tmp)
                  end do
c
                  r2ab(e,f,m,n) = r2ab(e,f,m,n) + iab(e+noa,f+nob,m,n) ! F
c
               end do
            end do
         end do
      end do
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('x',6)
c
      end
      double precision function uccsdtest_triples_pure(
     $     nmo,
     $     noa, nva, 
     $     iaa, 
     $     t1a,
     $     t2aa,
     $     evals)
      implicit none
      integer nmo, noa, nva
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t2aa(nva, nva, noa, noa)
      double precision evals(nmo)
c
      integer i, j, k, m, a, b, c, e
      double precision w, v, d
      double precision energy
c
      double precision ge, gm, gs
      ge(i,j,k,a,b,c,e) = t2aa(c,e,i,j)*iaa(a+noa,b+noa,e+noa,k)
      gm(i,j,k,a,b,c,m) =-t2aa(a,b,m,k)*iaa(c+noa,m,i,j)
      gs(i,j,k,a,b,c)   = t1a(c,k)*iaa(i,j,a+noa,b+noa)
c
      call output(evals,1,nmo,1,1,nmo,1,1)
c
      energy = 0.0d0
      do a = 1, nva
         do b = 1, a-1
            do c = 1, b-1
               do i = 1, noa
                  do j = 1, i-1
                     do k = 1, j-1
                        w = 0.0d0
                        do e = 1, nva
                           w = w 
     $                          + ge(i,j,k,a,b,c,e)
     $                          - ge(k,j,i,a,b,c,e)
     $                          - ge(i,k,j,a,b,c,e)
     $                          - ge(i,j,k,c,b,a,e)
     $                          - ge(i,j,k,a,c,b,e)
     $                          + ge(k,j,i,c,b,a,e)
     $                          + ge(k,j,i,a,c,b,e)
     $                          + ge(i,k,j,c,b,a,e)
     $                          + ge(i,k,j,a,c,b,e)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e): ',
c     &                              a+noa,b+noa,c+noa,i,j,k,w
                        do m = 1, noa
                           w = w 
     $                          + gm(i,j,k,a,b,c,m)
     $                          - gm(k,j,i,a,b,c,m)
     $                          - gm(i,k,j,a,b,c,m)
     $                          - gm(i,j,k,c,b,a,m)
     $                          - gm(i,j,k,a,c,b,m)
     $                          + gm(k,j,i,c,b,a,m)
     $                          + gm(k,j,i,a,c,b,m)
     $                          + gm(i,k,j,c,b,a,m)
     $                          + gm(i,k,j,a,c,b,m)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e+m): ',
c     &                              a+noa,b+noa,c+noa,i,j,k,w
                        v = w 
     $                       + gs(i,j,k,a,b,c)
     $                       - gs(k,j,i,a,b,c)
     $                       - gs(i,k,j,a,b,c)
     $                       - gs(i,j,k,c,b,a)
     $                       - gs(i,j,k,a,c,b)
     $                       + gs(k,j,i,c,b,a)
     $                       + gs(k,j,i,a,c,b)
     $                       + gs(i,k,j,c,b,a)
     $                       + gs(i,k,j,a,c,b)
                        d = evals(a+noa)+evals(b+noa)+evals(c+noa)
     $                       -evals(i)-evals(j)-evals(k)
c                        write(6,*)' a,b,c,i,j,k,w,v,d,energy: ',
c     &                        a+noa,b+noa,c+noa,i,j,k,w,v,d,energy
                        energy = energy - v*w/d
c
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      write(6,*) ' E ', energy
c
      uccsdtest_triples_pure = energy
c
      end
      double precision function uccsdtest_triples_pure_blocked(
     $     iaa,
     $     t1a,
     $     t2aa,
     $     evals, spin, nv_dim, no_dim)
      implicit none
#include "cuccsdtP.fh"
      integer spin
      integer nv_dim, no_dim
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision t1a(nv_dim, no_dim)
      double precision t2aa(nv_dim, nv_dim, no_dim, no_dim)
      double precision evals(nmo)
c
      integer a_blk, a_blk_lo, a_blk_hi, a_blk_sym
      integer b_blk, b_blk_lo, b_blk_hi, b_blk_sym
      integer c_blk, c_blk_lo, c_blk_hi, c_blk_sym
      integer i_blk, i_blk_lo, i_blk_hi, i_blk_sym
      integer j_blk, j_blk_lo, j_blk_hi, j_blk_sym
      integer k_blk, k_blk_lo, k_blk_hi, k_blk_sym
      integer m_blk, m_blk_lo, m_blk_hi, m_blk_sym
      integer e_blk, e_blk_lo, e_blk_hi, e_blk_sym
      double precision energy
c
c     temp buffers for t2 and ints (max_vblock = 3, max_oblock = 2)
c
      double precision t2aa_ecji(8,3,2,2), iaa_ekab(8,2,3,3),
     &                 t2aa_ecjk(8,3,2,2), iaa_eiab(8,2,3,3),
     &                 t2aa_ecki(8,3,2,2), iaa_ejab(8,2,3,3),
     &                 t2aa_eaji(8,3,2,2), iaa_ekcb(8,2,3,3),
     &                 t2aa_ebji(8,3,2,2), iaa_ekac(8,2,3,3),
     &                 t2aa_eajk(8,3,2,2), iaa_eicb(8,2,3,3),
     &                 t2aa_ebjk(8,3,2,2), iaa_eiac(8,2,3,3),
     &                 t2aa_eaki(8,3,2,2), iaa_ejcb(8,2,3,3),
     &                 t2aa_ebki(8,3,2,2), iaa_ejac(8,2,3,3)
c
      double precision t2aa_mkab(3,2,3,3), iaa_mcji(3,3,2,2),
     &                 t2aa_miab(3,2,3,3), iaa_mcjk(3,3,2,2),
     &                 t2aa_mjab(3,2,3,3), iaa_mcki(3,3,2,2),
     &                 t2aa_mkcb(3,2,3,3), iaa_maji(3,3,2,2),
     &                 t2aa_mkac(3,2,3,3), iaa_mbji(3,3,2,2),
     &                 t2aa_micb(3,2,3,3), iaa_majk(3,3,2,2),
     &                 t2aa_miac(3,2,3,3), iaa_mbjk(3,3,2,2),
     &                 t2aa_mjcb(3,2,3,3), iaa_maki(3,3,2,2),
     &                 t2aa_mjac(3,2,3,3), iaa_mbki(3,3,2,2)
c
      double precision iaa_ijab(2,2,3,3),
     &                 iaa_kjab(2,2,3,3),
     &                 iaa_ikab(2,2,3,3),
     &                 iaa_ijcb(2,2,3,3),
     &                 iaa_ijac(2,2,3,3),
     &                 iaa_kjcb(2,2,3,3),
     &                 iaa_kjac(2,2,3,3),
     &                 iaa_ikcb(2,2,3,3),
     &                 iaa_ikac(2,2,3,3)
c
      double precision int_buf(8*3*3*2)
c
c      double precision ge, gm, gs
c      ge(i,j,k,a,b,c,e) = t2aa(c,e,i,j)*
c     &                    iaa(a+no(spin),b+no(spin),e+no(spin),k)
c      gm(i,j,k,a,b,c,m) =-t2aa(a,b,m,k)*iaa(c+no(spin),m,i,j)
c      gs(i,j,k,a,b,c)   = t1a(c,k)*iaa(i,j,a+no(spin),b+no(spin))
c
c      ge(i,j,k,a,b,c,e) = t2aa(c-no(spin),e-no(spin),i,j)*
c     &                    iaa(a,b,e,k)
c      gm(i,j,k,a,b,c,m) =-t2aa(a-no(spin),b-no(spin),m,k)*
c     &                    iaa(c,m,i,j)
c      gs(i,j,k,a,b,c)   = t1a(c-no(spin),k)*iaa(i,j,a,b)
c
      call output(evals,1,nmo,1,1,nmo,1,1)
c
      energy = 0.0d0
c
      write(6,*)' spin, nvblock(spin): ',
     &            spin, nvblock(spin)
      do a_blk = 1, nvblock(spin)
         a_blk_lo = vblock(1,a_blk,spin)
         a_blk_hi = vblock(2,a_blk,spin)
         a_blk_sym = vblock_sym(a_blk,spin)
c
      write(6,*)' a_blk, a_blk_lo, a_blk_hi, a_blk_sym ',
     &            a_blk, a_blk_lo, a_blk_hi, a_blk_sym
         do b_blk = 1, a_blk
            b_blk_lo = vblock(1,b_blk,spin)
            b_blk_hi = vblock(2,b_blk,spin)
            b_blk_sym = vblock_sym(b_blk,spin)
c
      write(6,*)' b_blk, b_blk_lo, b_blk_hi, b_blk_sym ',
     &            b_blk, b_blk_lo, b_blk_hi, b_blk_sym
            do c_blk = 1, b_blk
               c_blk_lo = vblock(1,c_blk,spin)
               c_blk_hi = vblock(2,c_blk,spin)
               c_blk_sym = vblock_sym(c_blk,spin)
c
      write(6,*)' c_blk, c_blk_lo, c_blk_hi, c_blk_sym ',
     &            c_blk, c_blk_lo, c_blk_hi, c_blk_sym
               do i_blk = 1, noblock(spin)
                  i_blk_lo = oblock(1,i_blk,spin)
                  i_blk_hi = oblock(2,i_blk,spin)
                  i_blk_sym = oblock_sym(i_blk,spin)
c
      write(6,*)' i_blk, i_blk_lo, i_blk_hi, i_blk_sym ',
     &            i_blk, i_blk_lo, i_blk_hi, i_blk_sym
                  do j_blk = 1, i_blk
                     j_blk_lo = oblock(1,j_blk,spin)
                     j_blk_hi = oblock(2,j_blk,spin)
                     j_blk_sym = oblock_sym(j_blk,spin)
c
      write(6,*)' j_blk, j_blk_lo, j_blk_hi, j_blk_sym ',
     &            j_blk, j_blk_lo, j_blk_hi, j_blk_sym
                     do k_blk = 1, j_blk
                        k_blk_lo = oblock(1,k_blk,spin)
                        k_blk_hi = oblock(2,k_blk,spin)
                        k_blk_sym = oblock_sym(k_blk,spin)
c
      write(6,*)' k_blk, k_blk_lo, k_blk_hi, k_blk_sym ',
     &            k_blk, k_blk_lo, k_blk_hi, k_blk_sym
c
c     now get T's and Ints to build given permutation using the APIs
c     temporarily build these in scratch buffers to compare to 
c     what is coming out of the APIs
c
c     note e runs over all virtual and m over all occupied
c     later worry about symmetry blocking
c
                        e_blk_lo = no(spin)+1
                        e_blk_hi = no(spin)+nv(spin)
                        m_blk_lo = 1
                        m_blk_hi = no(spin)
c
c     zero all the buffers
c
                        call dfill (3*8*2*2, 0.0d0, t2aa_ecji, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ekab, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_ecjk, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_eiab, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_ecki, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ejab, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_eaji, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ekcb, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_ebji, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ekac, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_eajk, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_eicb, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_ebjk, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_eiac, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_eaki, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ejcb, 1)
                        call dfill (3*8*2*2, 0.0d0, t2aa_ebki, 1)
                        call dfill (3*3*8*2, 0.0d0, iaa_ejac, 1)
c
                        call dfill (3*3*3*2, 0.0d0, t2aa_mkab, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mcji, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_miab, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mcjk, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_mjab, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mcki, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_mkcb, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_maji, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_mkac, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mbji, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_micb, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_majk, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_miac, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mbjk, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_mjcb, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_maki, 1)
                        call dfill (3*3*3*2, 0.0d0, t2aa_mjac, 1)
                        call dfill (3*3*2*2, 0.0d0, iaa_mbki, 1)
c
                        call dfill (2*2*3*3, 0.0d0, iaa_ijab, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_kjab, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_ikab, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_ijcb, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_ijac, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_kjcb, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_kjac, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_ikcb, 1)
                        call dfill (2*2*3*3, 0.0d0, iaa_ikac, 1)
c
                        write(6,*)' after butt-load of dfills '
c
                        call triples_temp_fill_buffers(iaa, t2aa,
     &                     t2aa_ecji, iaa_ekab, t2aa_ecjk, iaa_eiab,
     &                     t2aa_ecki, iaa_ejab, t2aa_eaji, iaa_ekcb,
     &                     t2aa_ebji, iaa_ekac, t2aa_eajk, iaa_eicb,
     &                     t2aa_ebjk, iaa_eiac, t2aa_eaki, iaa_ejcb,
     &                     t2aa_ebki, iaa_ejac, t2aa_mkab, iaa_mcji,
     &                     t2aa_miab, iaa_mcjk, t2aa_mjab, iaa_mcki,
     &                     t2aa_mkcb, iaa_maji, t2aa_mkac, iaa_mbji,
     &                     t2aa_micb, iaa_majk, t2aa_miac, iaa_mbjk,
     &                     t2aa_mjcb, iaa_maki, t2aa_mjac, iaa_mbki,
     &                     iaa_ijab,  iaa_kjab, iaa_ikab,  iaa_ijcb,
     &                     iaa_ijac,  iaa_kjcb, iaa_kjac,  iaa_ikcb,
     &                     iaa_ikac,
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     nv_dim, no_dim, spin)
c
c     do api block reads to fill integral buffers: iaa_ekab, iaa_eiab, 
c                                                  iaa_ejab, iaa_ekcb
c                                                  iaa_ekac, iaa_eicb
c                                                  iaa_eiac, iaa_ejcb
c                                                  iaa_ejac
c
c     I(nv_sym(esym),dim(kblock),dim(ablock),dim(bblock)). 
c
                        call uccsdt_get_3x(k_blk, a_blk, b_blk, 
     &                     spin, spin, spin, iaa_ekab)
                        call uccsdt_get_3x(i_blk, a_blk, b_blk, 
     &                     spin, spin, spin, iaa_eiab)
                        call uccsdt_get_3x(j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, iaa_ejab)
                        call uccsdt_get_3x(k_blk, c_blk, b_blk, 
     &                     spin, spin, spin, iaa_ekcb)
                        call uccsdt_get_3x(k_blk, a_blk, c_blk, 
     &                     spin, spin, spin, iaa_ekac)
                        call uccsdt_get_3x(i_blk, c_blk, b_blk, 
     &                     spin, spin, spin, iaa_eicb)
                        call uccsdt_get_3x(i_blk, a_blk, c_blk, 
     &                     spin, spin, spin, iaa_eiac)
                        call uccsdt_get_3x(j_blk, c_blk, b_blk, 
     &                     spin, spin, spin, iaa_ejcb)
                        call uccsdt_get_3x(j_blk, a_blk, c_blk, 
     &                     spin, spin, spin, iaa_ejac)
                        call uccsdt_get_2x(i_blk, j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_ijab)
                        call uccsdt_get_2x(k_blk, j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_kjab)
                        call uccsdt_get_2x(i_blk, k_blk, a_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_ikab)
                        call uccsdt_get_2x(i_blk, j_blk, c_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_ijcb)
                        call uccsdt_get_2x(i_blk, j_blk, a_blk, c_blk, 
     &                     spin, spin, spin, spin, iaa_ijac)
                        call uccsdt_get_2x(k_blk, j_blk, c_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_kjcb)
                        call uccsdt_get_2x(k_blk, j_blk, a_blk, c_blk, 
     &                     spin, spin, spin, spin, iaa_kjac)
                        call uccsdt_get_2x(i_blk, k_blk, c_blk, b_blk, 
     &                     spin, spin, spin, spin, iaa_ikcb)
                        call uccsdt_get_2x(i_blk, k_blk, a_blk, c_blk, 
     &                     spin, spin, spin, spin, iaa_ikac)
                        call uccsdt_get_1x(c_blk, j_blk, i_blk, 
     &                     spin, spin, spin, iaa_mcji)
                        call uccsdt_get_1x(c_blk, j_blk, k_blk, 
     &                     spin, spin, spin, iaa_mcjk)
                        call uccsdt_get_1x(c_blk, k_blk, i_blk, 
     &                     spin, spin, spin, iaa_mcki)
                        call uccsdt_get_1x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, iaa_maji)
                        call uccsdt_get_1x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, iaa_mbji)
                        call uccsdt_get_1x(a_blk, j_blk, k_blk, 
     &                     spin, spin, spin, iaa_majk)
                        call uccsdt_get_1x(b_blk, j_blk, k_blk, 
     &                     spin, spin, spin, iaa_mbjk)
                        call uccsdt_get_1x(a_blk, k_blk, i_blk, 
     &                     spin, spin, spin, iaa_maki)
                        call uccsdt_get_1x(b_blk, k_blk, i_blk, 
     &                     spin, spin, spin, iaa_mbki)
                        call uccsdt_get_t3x(c_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t2aa_ecji)
                        call uccsdt_get_t3x(c_blk, j_blk, k_blk, 
     &                     spin, spin, spin, t2aa_ecjk)
                        call uccsdt_get_t3x(c_blk, k_blk, i_blk, 
     &                     spin, spin, spin, t2aa_ecki)
                        call uccsdt_get_t3x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t2aa_eaji)
                        call uccsdt_get_t3x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t2aa_ebji)
                        call uccsdt_get_t3x(a_blk, j_blk, k_blk, 
     &                     spin, spin, spin, t2aa_eajk)
                        call uccsdt_get_t3x(b_blk, j_blk, k_blk, 
     &                     spin, spin, spin, t2aa_ebjk)
                        call uccsdt_get_t3x(a_blk, k_blk, i_blk, 
     &                     spin, spin, spin, t2aa_eaki)
                        call uccsdt_get_t3x(b_blk, k_blk, i_blk, 
     &                     spin, spin, spin, t2aa_ebki)
                        call uccsdt_get_t1x(k_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call writeblock ('int_buf', int_buf, 
     &                     m_blk_lo, m_blk_hi, k_blk_lo, k_blk_hi, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi)
c                        call uccsdt_get_t1x(k_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_mkab)
c                        call uccsdt_get_t1x(i_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_miab)
c                        call uccsdt_get_t1x(j_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_mjab)
c                        call uccsdt_get_t1x(k_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_mkcb)
c                        call uccsdt_get_t1x(k_blk, a_blk, c_blk, 
c     &                     spin, spin, spin, t2aa_mkac)
c                        call uccsdt_get_t1x(i_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_micb)
c                        call uccsdt_get_t1x(i_blk, a_blk, c_blk, 
c     &                     spin, spin, spin, t2aa_miac)
c                        call uccsdt_get_t1x(j_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, t2aa_mjcb)
c                        call uccsdt_get_t1x(j_blk, a_blk, c_blk, 
c     &                     spin, spin, spin, t2aa_mjac)
c
c                        call writeblock ('int_buf', int_buf, 
c     &                     e_blk_lo, e_blk_hi, k_blk_lo, k_blk_hi, 
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi)
c                        call errquit(' temp stop ', 911)
c
                        call uccsdtest_triples_pure_blocked_energy
     &                    (iaa, t2aa, t1a, evals, energy,
     &                     t2aa_ecji, iaa_ekab, t2aa_ecjk, iaa_eiab,
     &                     t2aa_ecki, iaa_ejab, t2aa_eaji, iaa_ekcb,
     &                     t2aa_ebji, iaa_ekac, t2aa_eajk, iaa_eicb,
     &                     t2aa_ebjk, iaa_eiac, t2aa_eaki, iaa_ejcb,
     &                     t2aa_ebki, iaa_ejac, t2aa_mkab, iaa_mcji,
     &                     t2aa_miab, iaa_mcjk, t2aa_mjab, iaa_mcki,
     &                     t2aa_mkcb, iaa_maji, t2aa_mkac, iaa_mbji,
     &                     t2aa_micb, iaa_majk, t2aa_miac, iaa_mbjk,
     &                     t2aa_mjcb, iaa_maki, t2aa_mjac, iaa_mbki,
     &                     iaa_ijab,  iaa_kjab, iaa_ikab,  iaa_ijcb,
     &                     iaa_ijac,  iaa_kjcb, iaa_kjac,  iaa_ikcb,
     &                     iaa_ikac,
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     nv_dim, no_dim, spin)
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      write(6,*) ' E ', energy
c
      uccsdtest_triples_pure_blocked = energy
c
      end
      double precision function uccsdtest_triples_mixed(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     evalsa, evalsb)
c
      implicit none
      integer nmo, noa, nob, nva, nvb
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob) 
      double precision evalsa(nmo), evalsb(nmo)
c
      integer i, j, k, m, a, b, c, e
      double precision w, v, d, energy
c
      call output(evalsa,1,nmo,1,1,nmo,1,1)
      call output(evalsb,1,nmo,1,1,nmo,1,1)
c
      energy = 0.0d0
      do a = 1, nva
         do b = 1, a-1
            do c = 1, nvb
               do i = 1, noa
                  do j = 1, i-1
                     do k = 1, nob
                        w = 0.0d0
                        do e = 1, nva
                           w = w 
     $                          - t2ab(e,c,j,k)*iaa(a+noa,b+noa,e+noa,i)
     $                          + t2ab(e,c,i,k)*iaa(a+noa,b+noa,e+noa,j)
     $                          + t2aa(a,e,i,j)*iab(b+noa,c+nob,e+noa,k)
     $                          - t2aa(b,e,i,j)*iab(a+noa,c+nob,e+noa,k)
                        end do
                        do e = 1, nvb
                           w = w 
     $                          - t2ab(a,e,j,k)*iab(b+noa,c+nob,i,e+nob)
     $                          + t2ab(b,e,j,k)*iab(a+noa,c+nob,i,e+nob)
     $                          + t2ab(a,e,i,k)*iab(b+noa,c+nob,j,e+nob)
     $                          - t2ab(b,e,i,k)*iab(a+noa,c+nob,j,e+nob)
                        end do
                        do m = 1, noa
                           w = w
     $                          + t2aa(a,b,m,i)*iab(m,c+nob,j,k)
     $                          - t2aa(a,b,m,j)*iab(m,c+nob,i,k)
     $                          - t2ab(b,c,m,k)*iaa(a+noa,m,i,j)
     $                          + t2ab(a,c,m,k)*iaa(b+noa,m,i,j)
                        end do
                        do m = 1, nob
                           w = w 
     $                          + t2ab(b,c,i,m)*iab(a+noa,m,j,k)
     $                          - t2ab(a,c,i,m)*iab(b+noa,m,j,k)
     $                          - t2ab(b,c,j,m)*iab(a+noa,m,i,k)
     $                          + t2ab(a,c,j,m)*iab(b+noa,m,i,k)
                        end do
                        v = w +
     $                       + t1b(c,k)*iaa(i,j,a+noa,b+noa)
     $                       + t1a(a,i)*iab(j,k,b+noa,c+nob)
     $                       - t1a(b,i)*iab(j,k,a+noa,c+nob)
     $                       - t1a(a,j)*iab(i,k,b+noa,c+nob)
     $                       + t1a(b,j)*iab(i,k,a+noa,c+nob)
c
                        d = evalsa(a+noa)+evalsa(b+noa)+evalsb(c+noa)
     $                       -evalsa(i)-evalsa(j)-evalsb(k)
c
                        energy = energy - v*w/d
c
c$$$                        if (abs(v*w/d) .gt. 1e-10) then
c$$$                           write(6,1) i,j,k,a,b,c,v,w,d,energy
c$$$ 1                         format(6i5,1p,4d16.8)
c$$$                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      write(6,*) ' energy mixed ', energy
c
      uccsdtest_triples_mixed = energy
c
      end
      subroutine uccsdt_initialize(nummo, noa, nob, nva, nvb, 
     $     epsa, epsb, kk_iaa, kk_ibb, kk_iab, kk_t1a, kk_t1b,
     $     kk_t2aa,kk_t2bb,kk_t2ab)
      implicit none
#include "cuccsdtP.fh"
      integer nummo, noa, nob, nva, nvb
      double precision epsa(*), epsb(*)
      integer kk_iaa, kk_ibb, kk_iab, kk_t1a, kk_t1b,
     $     kk_t2aa,kk_t2bb,kk_t2ab
c
      integer i, ir, spin, block
c
      nmo = nummo
      nact = nmo
      nir = 1
      nc(1) = 0
      nc(2) = 0
      no(1) = noa
      no(2) = nob
      nv(1) = nva
      nv(2) = nvb
c
      call ifill(8*2, 0, nc_sym, 1)
      call ifill(8*2, 0, no_sym, 1)
      call ifill(8*2, 0, nv_sym, 1)
      no_sym(0,1) = noa
      no_sym(0,2) = nob
      nv_sym(0,1) = nva
      nv_sym(0,2) = nvb
c
      do spin = 1, 2
         call ifill(8*2, 0, c_sym(1,0,spin), 2)
         call ifill(8*2,-1, c_sym(2,0,spin), 2)
         call ifill(8*2, 0, o_sym(1,0,spin), 2)
         call ifill(8*2,-1, o_sym(2,0,spin), 2)
         call ifill(8*2, 0, v_sym(1,0,spin), 2)
         call ifill(8*2,-1, v_sym(2,0,spin), 2)
      end do
      o_sym(1,0,1) = 1
      o_sym(2,0,1) = noa
      o_sym(1,0,2) = 1
      o_sym(2,0,2) = nob
      v_sym(1,0,1) = noa+1
      v_sym(2,0,1) = nmo
      v_sym(1,0,2) = nob+1
      v_sym(2,0,2) = nmo
c
      call ifill(nmo, -999999999, map, 1)
c
      oblock_max = 2
      vblock_max = 3
      do spin = 1, 2
         nvblock(spin) = 0
         noblock(spin) = 0
         do ir = 0, nir-1
            do i = v_sym(1,ir,spin),v_sym(2,ir,spin),vblock_max
               nvblock(spin) = nvblock(spin) + 1
               vblock_sym(nvblock(spin),spin) = 1
               vblock(1,nvblock(spin),spin) = i
               vblock(2,nvblock(spin),spin) = 
     $              min(v_sym(2,ir,spin),i+vblock_max-1)
            end do
            do i = o_sym(1,ir,spin),o_sym(2,ir,spin),oblock_max
               noblock(spin) = noblock(spin) + 1
               oblock_sym(noblock(spin),spin) = 0
               oblock(1,noblock(spin),spin) = i
               oblock(2,noblock(spin),spin) = 
     $              min(o_sym(2,ir,spin),i+oblock_max-1)
            end do
         end do
      end do
c
c     This just temporary stuff for development
c
      call dcopy(nmo, epsa, 1, eps(1,1), 1)
      call dcopy(nmo, epsb, 1, eps(1,2), 1)
      k_moints_aa = kk_iaa
      k_moints_ab = kk_iab
      k_moints_bb = kk_ibb
      k_taa = kk_t2aa
      k_tab = kk_t2ab
      k_tbb = kk_t2bb
      k_ta  = kk_t1a
      k_tb  = kk_t1b
c
      write(6,*) ' nact ', nact
      write(6,*) ' nmo  ', nmo
      write(6,*) ' nir  ', nir
      write(6,*) ' nc   ', nc
      write(6,*) ' no   ', no
      write(6,*) ' nv   ', nv
      write(6,*) ' nc_sym ', ((nc_sym(ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' no_sym ', ((no_sym(ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' nv_sym ', ((nv_sym(ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' c_sym(1) ', ((c_sym(1,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' c_sym(2) ', ((c_sym(2,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' o_sym(1) ', ((o_sym(1,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' o_sym(2) ', ((o_sym(2,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' v_sym(1) ', ((v_sym(1,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' v_sym(2) ', ((v_sym(2,ir,spin),ir=0,nir-1),spin=1,2)
      write(6,*) ' oblock_max ', oblock_max
      write(6,*) ' vblock_max ', vblock_max
      write(6,*) ' noblock ', noblock
      write(6,*) ' nvblock ', nvblock
      write(6,*) ' oblock(1) ', 
     $     ((oblock(1,block,spin),block=1,noblock(spin)),spin=1,2)
      write(6,*) ' oblock(2) ', 
     $     ((oblock(2,block,spin),block=1,noblock(spin)),spin=1,2)
      write(6,*) ' vblock(1) ', 
     $     ((vblock(1,block,spin),block=1,nvblock(spin)),spin=1,2)
      write(6,*) ' vblock(2) ', 
     $     ((vblock(2,block,spin),block=1,nvblock(spin)),spin=1,2)
c
      write(6,*) ' eps ', ((eps(i,spin),i=1,nmo),spin=1,2)
c
      end
      subroutine uccsdt_get_3x(kblock, ablock, bblock, 
     $     spink, spina, spinb, ints)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
      integer kblock, ablock, bblock, spink, spina, spinb
      double precision ints(*)
c
      logical oflip_ab, oflip_ek
      integer base,spine,syme
      integer alo, ahi, blo, bhi, klo, khi, elo, ehi
c
c     Returns <ek||ab> for orbitals in the specified triples blocks and for given
c     spins of each block.  The integrals will be in the array
c     I(nv_sym(esym),dim(kblock),dim(ablock),dim(bblock)).  
c
c     Mixed spin integrals will get one of <ek|ab> or -<ek|ba>.  
c
      oflip_ab = .false.
      oflip_ek = .false.
      if (spink*spina*spinb .eq. 1) then ! <ek||ab> pure alpha
         base = k_moints_aa
         spine = 1
      else if (spink*spina*spinb .eq. 8) then ! <ek||ab> pure beta
         base = k_moints_bb
         spine = 2
      else
         base = k_moints_ab
         oflip_ek = spink.eq.1
         oflip_ab = spinb.eq.1
         spine = 1
         if (spink.eq.1) spine = 2
      endif
c
      syme = ieor(oblock_sym(kblock,spink),
     $     ieor(vblock_sym(ablock,spina),vblock_sym(bblock,spinb)))
c
      alo = vblock(1,ablock,spina)
      ahi = vblock(2,ablock,spina)
      blo = vblock(1,bblock,spinb)
      bhi = vblock(2,bblock,spinb)
      klo = oblock(1,kblock,spink)
      khi = oblock(2,kblock,spink)
      elo = v_sym(1,syme,spine)
      ehi = v_sym(2,syme,spine)
c
      write(6,*)' alo, ahi, blo, bhi, klo, khi, elo, ehi: ',
     &            alo, ahi, blo, bhi, klo, khi, elo, ehi
c
      call uccsdt_copy_int_block
     $     (ints,dbl_mb(base),nmo,nmo,nmo,
     $     elo,ehi,klo,khi,alo,ahi,blo,bhi,
     $     oflip_ek,oflip_ab,.false.)
c
      end
      subroutine uccsdt_get_2x(iblock, jblock, ablock, bblock, 
     $     spini, spinj, spina, spinb, ints)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
      integer iblock, jblock, ablock, bblock, spini, spinj, spina, spinb
      double precision ints(*)
c
      logical oflip_ab, oflip_ij
      integer base
      integer alo, ahi, blo, bhi, ilo, ihi, jlo, jhi
c
c     Returns <ij||ab> for orbitals in the specified triples blocks and for given
c     spins of each block.  The integrals will be in the array
c     I(dim(iblock),dim(jblock),dim(ablock),dim(bblock)).  
c
c     Mixed spin integrals will get one of <ij|ab> or -<ij|ba>.  
c
      oflip_ab = .false.
      oflip_ij = .false.
      if (spini*spinj*spina*spinb .eq. 1) then ! <ij||ab> pure alpha
         base = k_moints_aa
      else if (spini*spinj*spina*spinb .eq. 16) then ! <ij||ab> pure beta
         base = k_moints_bb
      else
         base = k_moints_ab
         oflip_ij = spinj.eq.1
         oflip_ab = spinb.eq.1
      endif
c
      alo = vblock(1,ablock,spina)
      ahi = vblock(2,ablock,spina)
      blo = vblock(1,bblock,spinb)
      bhi = vblock(2,bblock,spinb)
      ilo = oblock(1,iblock,spini)
      ihi = oblock(2,iblock,spini)
      jlo = oblock(1,jblock,spinj)
      jhi = oblock(2,jblock,spinj)
c
      call uccsdt_copy_int_block
     $     (ints,dbl_mb(base),nmo,nmo,nmo,
     $     ilo,ihi,jlo,jhi,alo,ahi,blo,bhi,
     $     oflip_ij,oflip_ab,.false.)
c
      end
      subroutine uccsdt_get_1x(cblock, iblock, jblock, 
     $     spinc, spini, spinj, ints)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
      integer cblock, iblock, jblock, spinc, spini, spinj
      double precision ints(*)
c
      logical oflip_ij, oflip_mc
      integer base,spinm,symm
      integer ilo, ihi, jlo, jhi, clo, chi, mlo, mhi
c
c     Returns <mc||ij> for orbitals in the specified triples blocks and for given
c     spins of each block.  The integrals will be in the array
c     I(no_sym(msym,spin1),dim(cblock),dim(iblock),dim(jblock)).
c
c     Mixed spin integrals will get one of <mc|ij> or -<mc|ji>.  
c
      oflip_ij = .false.
      oflip_mc = .false.
      if (spinc*spini*spinj .eq. 1) then ! <mc||ij> pure alpha
         base = k_moints_aa
         spinm = 1
      else if (spinc*spini*spinj .eq. 8) then ! <mc||ij> pure beta
         base = k_moints_bb
         spinm = 2
      else
         base = k_moints_ab
         oflip_mc = spinc.eq.1
         oflip_ij = spinj.eq.1
         spinm = 1
         if (spinc.eq.1) spinm = 2
      endif
c
      symm = ieor(vblock_sym(cblock,spinc),
     $     ieor(oblock_sym(iblock,spini),oblock_sym(jblock,spinj)))
c
      ilo = oblock(1,iblock,spini)
      ihi = oblock(2,iblock,spini)
      jlo = oblock(1,jblock,spinj)
      jhi = oblock(2,jblock,spinj)
      clo = vblock(1,cblock,spinc)
      chi = vblock(2,cblock,spinc)
      mlo = o_sym(1,symm,spinm)
      mhi = o_sym(2,symm,spinm)
c
      call uccsdt_copy_int_block
     $     (ints,dbl_mb(base),nmo,nmo,nmo,
     $     mlo,mhi,clo,chi,ilo,ihi,jlo,jhi,
     $     oflip_mc,oflip_ij,.false.)
c
      end
      subroutine uccsdt_copy_int_block
     $     (ints,g,dim1,dim2,dim3,
     $     elo,ehi,klo,khi,alo,ahi,blo,bhi,
     $     oflip_ek,oflip_ab,oflip_ek_ab)
      implicit none
      integer dim1, dim2, dim3
      double precision g(dim1,dim2,dim3,*)
      integer elo,ehi,klo,khi,alo,ahi,blo,bhi
      logical oflip_ab,oflip_ek,oflip_ek_ab
      double precision ints(*)
c
      integer e,k,a,b,ind
c
      ind = 0
      do b = blo, bhi
         do a = alo, ahi
            do k = klo, khi
               if (oflip_ek_ab) then
                  if (oflip_ek) then
                     if (oflip_ab) then
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = g(b,a,k,e)
                        end do
                     else
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = -g(a,b,k,e)
                        end do
                     end if
                  else
                     if (oflip_ab) then
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = -g(b,a,e,k)
                        end do
                     else
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = g(a,b,e,k)
                        end do
                     end if
                  end if
               else
                  if (oflip_ek) then
                     if (oflip_ab) then
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = g(k,e,b,a)
                        end do
                     else
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = -g(k,e,a,b)
                        end do
                     end if
                  else
                     if (oflip_ab) then
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = -g(e,k,b,a)
                        end do
                     else
                        do e = elo, ehi
                           ind = ind + 1
                           ints(ind) = g(e,k,a,b)
                        end do
                     end if
                  end if
               end if
            end do
         end do
      end do
c
      end
      subroutine uccsdt_get_t3x(cblock, iblock, jblock, 
     $     spinc, spini, spinj, T)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
      integer cblock, iblock, jblock, spinc, spini, spinj
      double precision T(*)
c
c     Returns T(e,c,i,j) with orbitals i & e corresponding to electon 1.
c
c     T(nv_sym(esym,spin1),dim(cblock),dim(iblock),dim(jblock))
c
      logical oflip_ij, oflip_ec
      integer base,spine,syme
      integer ilo, ihi, jlo, jhi, clo, chi, elo, ehi
      integer dim1, dim2, dim3
c
c     Returns <mc||ij> for orbitals in the specified triples blocks and for given
c     spins of each block.  The integrals will be in the array
c     I(no_sym(msym,spin1),dim(cblock),dim(iblock),dim(jblock)).
c
c     Mixed spin integrals will get one of <mc|ij> or -<mc|ji>.  
c
      oflip_ij = .false.
      oflip_ec = .false.
      if (spinc*spini*spinj .eq. 1) then ! <ec||ij> pure alpha
         base = k_taa
         spine = 1
         dim1 = nv(1)
         dim2 = nv(1)
         dim3 = no(1)
      else if (spinc*spini*spinj .eq. 8) then ! <ec||ij> pure beta
         base = k_tbb
         spine = 2
         dim1 = nv(2)
         dim2 = nv(2)
         dim3 = no(2)
      else
         base = k_tab
         oflip_ec = spinc.eq.1
         oflip_ij = spinj.eq.1
         spine = 1
         if (spinc.eq.1) spine = 2
         dim1 = nv(1)
         dim2 = nv(2)
         dim3 = no(1)
      endif
c
      syme = ieor(vblock_sym(cblock,spinc),
     $     ieor(oblock_sym(iblock,spini),oblock_sym(jblock,spinj)))
c
      ilo = oblock(1,iblock,spini)
      ihi = oblock(2,iblock,spini)
      jlo = oblock(1,jblock,spinj)
      jhi = oblock(2,jblock,spinj)
      clo = vblock(1,cblock,spinc) - no(spinc)
      chi = vblock(2,cblock,spinc) - no(spinc)
      elo = v_sym(1,syme,spine)    - no(spine)
      ehi = v_sym(2,syme,spine)    - no(spine)
c
      call uccsdt_copy_int_block
     $     (T,dbl_mb(base),dim1,dim2,dim3,
     $     elo,ehi,clo,chi,ilo,ihi,jlo,jhi,
     $     oflip_ec,oflip_ij,.false.)
c
      end
      subroutine uccsdt_get_t1x(kblock, ablock, bblock, 
     $     spink, spina, spinb, T)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
      integer kblock, ablock, bblock, spink, spina, spinb
      double precision T(*)
c
      logical oflip_ab, oflip_mk
      integer base,spinm,symm,dim1,dim2,dim3
      integer alo, ahi, blo, bhi, klo, khi, mlo, mhi
c
c     Returns <mk||ab> for orbitals in the specified triples blocks and for given
c     spins of each block.  The integrals will be in the array
c     I(nv_sym(esym),dim(kblock),dim(ablock),dim(bblock)).  
c
c     Mixed spin integrals will get one of <mk|ab> or -<mk|ba>.  
c
      oflip_ab = .false.
      oflip_mk = .false.
      if (spink*spina*spinb .eq. 1) then ! <mk||ab> pure alpha
         base = k_taa
         spinm = 1
         dim1 = nv(1)
         dim2 = nv(1)
         dim3 = no(1)
      else if (spink*spina*spinb .eq. 8) then ! <mk||ab> pure beta
         base = k_tbb
         spinm = 2
         dim1 = nv(2)
         dim2 = nv(2)
         dim3 = no(2)
      else
         base = k_tab
         oflip_mk = spink.eq.1
         oflip_ab = spinb.eq.1
         spinm = 1
         if (spink.eq.1) spinm = 2
         dim1 = nv(1)
         dim2 = nv(2)
         dim3 = no(1)
      endif
c
      symm = ieor(oblock_sym(kblock,spink),
     $     ieor(vblock_sym(ablock,spina),vblock_sym(bblock,spinb)))
c
      alo = vblock(1,ablock,spina)
      ahi = vblock(2,ablock,spina)
      blo = vblock(1,bblock,spinb)
      bhi = vblock(2,bblock,spinb)
      klo = oblock(1,kblock,spink)
      khi = oblock(2,kblock,spink)
      mlo = o_sym(1,symm,spinm)
      mhi = o_sym(2,symm,spinm)
c
      call uccsdt_copy_int_block
     $     (T,dbl_mb(base),dim1,dim2,dim3,
     $     mlo,mhi,klo,khi,alo,ahi,blo,bhi,
     $     oflip_mk,oflip_ab,.true.)
c
      end
      subroutine util_iswap2(i,j)
      integer i, j, tmp
c
c     Swap the arguments
c
      tmp = i
      i = j
      j = tmp
c     
      end
      subroutine writeblock (buf_title, ints, 
     &                       wlo, whi, xlo, xhi, ylo, yhi, zlo, zhi)
      integer wlo, whi, xlo, xhi, ylo, yhi, zlo, zhi, w, x, y, z
      double precision ints(wlo:whi, xlo:xhi, ylo:yhi, zlo:zhi)
      character*(*) buf_title
      write(6,*)buf_title
      do w = wlo, whi
         do x = xlo, xhi
            do y = ylo, yhi
               do z = zlo, zhi
                  if (abs(ints(w,x,y,z)).gt.1d-10)
     &               write(6,*)' w, x, y, z, ints(w,x,y,z): ',
     &                           w, x, y, z, ints(w,x,y,z)
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine triples_temp_fill_buffers(iaa, t2aa,
     &   t2aa_ecji, iaa_ekab, t2aa_ecjk, iaa_eiab,
     &   t2aa_ecki, iaa_ejab, t2aa_eaji, iaa_ekcb,
     &   t2aa_ebji, iaa_ekac, t2aa_eajk, iaa_eicb,
     &   t2aa_ebjk, iaa_eiac, t2aa_eaki, iaa_ejcb,
     &   t2aa_ebki, iaa_ejac, t2aa_mkab, iaa_mcji,
     &   t2aa_miab, iaa_mcjk, t2aa_mjab, iaa_mcki,
     &   t2aa_mkcb, iaa_maji, t2aa_mkac, iaa_mbji,
     &   t2aa_micb, iaa_majk, t2aa_miac, iaa_mbjk,
     &   t2aa_mjcb, iaa_maki, t2aa_mjac, iaa_mbki,
     &   iaa_ijab,  iaa_kjab, iaa_ikab,  iaa_ijcb,
     &   iaa_ijac,  iaa_kjcb, iaa_kjac,  iaa_ikcb,
     &   iaa_ikac,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   nv_dim, no_dim, spin)
c
      implicit none
#include "cuccsdtP.fh"
      integer          alo, ahi,
     &                 blo, bhi,
     &                 clo, chi,
     &                 ilo, ihi,
     &                 jlo, jhi,
     &                 klo, khi,
     &                 elo, ehi,
     &                 mlo, mhi,
     &                 nv_dim, no_dim, spin
c
      double precision iaa(nmo, nmo, nmo, nmo)
c
      double precision t2aa(nv_dim, nv_dim, no_dim, no_dim)
c
      double precision t2aa_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision iaa_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      double precision t2aa_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision iaa_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      double precision t2aa_ecki(elo:ehi,clo:chi,klo:khi,ilo:ihi)
      double precision iaa_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      double precision t2aa_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision iaa_ekcb(elo:ehi,klo:khi,clo:chi,blo:bhi)
      double precision t2aa_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision iaa_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      double precision t2aa_eajk(elo:ehi,alo:ahi,jlo:jhi,klo:khi)
      double precision iaa_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision t2aa_ebjk(elo:ehi,blo:bhi,jlo:jhi,klo:khi)
      double precision iaa_eiac(elo:ehi,ilo:ihi,alo:ahi,clo:chi)
      double precision t2aa_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision iaa_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      double precision t2aa_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision iaa_ejac(elo:ehi,jlo:jhi,alo:ahi,clo:chi)
      double precision t2aa_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision iaa_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision t2aa_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision iaa_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      double precision t2aa_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_mcki(mlo:mhi,clo:chi,klo:khi,ilo:ihi)
      double precision t2aa_mkcb(mlo:mhi,klo:khi,clo:chi,blo:bhi)
      double precision iaa_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision t2aa_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision iaa_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision t2aa_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision iaa_majk(mlo:mhi,alo:ahi,jlo:jhi,klo:khi)
      double precision t2aa_miac(mlo:mhi,ilo:ihi,alo:ahi,clo:chi)
      double precision iaa_mbjk(mlo:mhi,blo:bhi,jlo:jhi,klo:khi)
      double precision t2aa_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      double precision t2aa_mjac(mlo:mhi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      double precision iaa_ijab(ilo:ihi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_kjab(klo:khi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_ikab(ilo:ihi,klo:khi,alo:ahi,blo:bhi)
      double precision iaa_ijcb(ilo:ihi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_ijac(ilo:ihi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_kjcb(klo:khi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_kjac(klo:khi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_ikcb(ilo:ihi,klo:khi,clo:chi,blo:bhi)
      double precision iaa_ikac(ilo:ihi,klo:khi,alo:ahi,clo:chi)
      integer i, j, k, m, a, b, c, e
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            do c = clo, min(chi,b-1)
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     do k = klo, min(khi,j-1)
                        do e = no(spin)+1, no(spin)+nv(spin)
c                           write(6,*)' i, j, k, m, a, b, c, e ',
c     &                                 i, j, k, m, a, b, c, e
c
                           t2aa_ecji(e,c,j,i) = 
     &                     t2aa(c-no(spin),e-no(spin),i,j)
c
c                  if (abs(iaa(a,b,e,k)).gt.1d-10)
c     &               write(6,*)' a, b, e, k, iaa(a,b,e,k): ',
c     &                           a, b, e, k, iaa(a,b,e,k)
c
                           
                           iaa_ekab(e,k,a,b) = iaa(a,b,e,k)
c
c       if (abs(iaa_abek(a,b,e,k)).gt.1d-10)
c     &    write(6,*)' a, b, e, k, iaa_abek(a,b,e,k): ',
c     &                a, b, e, k, iaa_abek(a,b,e,k)
c
                           t2aa_ecjk(e,c,j,k) = 
     &                     t2aa(c-no(spin),e-no(spin),k,j)
c
                           iaa_eiab(e,i,a,b) = iaa(a,b,e,i)
c
                           t2aa_ecki(e,c,k,i) = 
     &                     t2aa(c-no(spin),e-no(spin),i,k)
c
                           iaa_ejab(e,j,a,b) = iaa(a,b,e,j)
c
                           t2aa_eaji(e,a,j,i) = 
     &                     t2aa(a-no(spin),e-no(spin),i,j)
c
                           iaa_ekcb(e,k,c,b) = iaa(c,b,e,k)
c
                           t2aa_ebji(e,b,j,i) = 
     &                     t2aa(b-no(spin),e-no(spin),i,j)
c
                           iaa_ekac(e,k,a,c) = iaa(a,c,e,k)
c
                           t2aa_eajk(e,a,j,k) = 
     &                     t2aa(a-no(spin),e-no(spin),k,j)
c
                           iaa_eicb(e,i,c,b) = iaa(c,b,e,i)
c
                           t2aa_ebjk(e,b,j,k) = 
     &                     t2aa(b-no(spin),e-no(spin),k,j)
c
                           iaa_eiac(e,i,a,c) = iaa(a,c,e,i)
c
                           t2aa_eaki(e,a,k,i) = 
     &                     t2aa(a-no(spin),e-no(spin),i,k)
c
                           iaa_ejcb(e,j,c,b) = iaa(c,b,e,j)
c
                           t2aa_ebki(e,b,k,i) = 
     &                     t2aa(b-no(spin),e-no(spin),i,k)
c
                           iaa_ejac(e,j,a,c) = iaa(a,c,e,j)
c
                        end do
c
                        do m = 1, no(spin)
c
                           t2aa_mkab(m,k,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,k)
c
                           iaa_mcji(m,c,j,i) = iaa(c,m,i,j)
c
                           t2aa_miab(m,i,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,i)
c
                           iaa_mcjk(m,c,j,k) = iaa(c,m,k,j)
c
                           t2aa_mjab(m,j,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,j)
c
                           iaa_mcki(m,c,k,i) = iaa(c,m,i,k)
c
                           t2aa_mkcb(m,k,c,b) = 
     &                     t2aa(c-no(spin),b-no(spin),m,k)
c
                           iaa_maji(m,a,j,i) = iaa(a,m,i,j)
c
                           t2aa_mkac(m,k,a,c) = 
     &                     t2aa(a-no(spin),c-no(spin),m,k)
c
                           iaa_mbji(m,b,j,i) = iaa(b,m,i,j)
c
                           t2aa_micb(m,i,c,b) = 
     &                     t2aa(c-no(spin),b-no(spin),m,i)
c
                           iaa_majk(m,a,j,k) = iaa(a,m,k,j)
c
                           t2aa_miac(m,i,a,c) = 
     &                     t2aa(a-no(spin),c-no(spin),m,i)
c
                           iaa_mbjk(m,b,j,k) = iaa(b,m,k,j)
c
                           t2aa_mjcb(m,j,c,b) = 
     &                     t2aa(c-no(spin),b-no(spin),m,j)
c
                           iaa_maki(m,a,k,i) = iaa(a,m,i,k)
c
                           t2aa_mjac(m,j,a,c) = 
     &                     t2aa(a-no(spin),c-no(spin),m,j)
c
                           iaa_mbki(m,b,k,i) = iaa(b,m,i,k)
c
                        end do
c
c                       assume t1 in memory set up integral blocks only
c
                        iaa_ijab(i,j,a,b) = iaa(i,j,a,b)
c
                        iaa_kjab(k,j,a,b) = iaa(k,j,a,b)
c
                        iaa_ikab(i,k,a,b) = iaa(i,k,a,b)
c
                        iaa_ijcb(i,j,c,b) = iaa(i,j,c,b)
c
                        iaa_ijac(i,j,a,c) = iaa(i,j,a,c)
c
                        iaa_kjcb(k,j,c,b) = iaa(k,j,c,b)
c
                        iaa_kjac(k,j,a,c) = iaa(k,j,a,c)
c
                        iaa_ikcb(i,k,c,b) = iaa(i,k,c,b)
c
                        iaa_ikac(i,k,a,c) = iaa(i,k,a,c)
c
                     end do
                  end do
               end do
            end do
         end do
      end do
c
c end block loops used to fill buffers
c
c temp write buffers
c
c      call writeblock('t2aa_ceij', t2aa_ceij, 
c     &                clo, chi, elo, ehi, ilo, ihi, jlo, jhi)
c      call writeblock('iaa_ekab', iaa_ekab, 
c     &                elo, ehi, klo, khi, alo, ahi, blo, bhi)
c      call writeblock('t2aa_cekj', t2aa_cekj, 
c     &                clo, chi, elo, ehi, klo, khi, jlo, jhi)
c      call writeblock('iaa_eiab', iaa_eiab, 
c     &                elo, ehi, ilo, ihi, alo, ahi, blo, bhi)
c      call writeblock('t2aa_ceik', t2aa_ceik, 
c     &                clo, chi, elo, ehi, ilo, ihi, klo, khi)
c      call writeblock('iaa_ejab', iaa_ejab, 
c     &                elo, ehi, jlo, jhi, alo, ahi, blo, bhi)
c      call writeblock('t2aa_aeij', t2aa_aeij, 
c     &                alo, ahi, elo, ehi, ilo, ihi, jlo, jhi)
c      call writeblock('iaa_ekcb', iaa_ekcb, 
c     &                elo, ehi, klo, khi, clo, chi, blo, bhi)
c      call writeblock('t2aa_beij', t2aa_beij, 
c     &                blo, bhi, elo, ehi, ilo, ihi, jlo, jhi)
c      call writeblock('iaa_ekac', iaa_ekac, 
c     &                elo, ehi, klo, khi, alo, ahi, clo, chi)
c      call writeblock('t2aa_aekj', t2aa_aekj, 
c     &                alo, ahi, elo, ehi, klo, khi, jlo, jhi)
c      call writeblock('iaa_eicb', iaa_eicb, 
c     &                elo, ehi, ilo, ihi, clo, chi, blo, bhi)
c      call writeblock('t2aa_bekj', t2aa_bekj, 
c     &                blo, bhi, elo, ehi, klo, khi, jlo, jhi)
c      call writeblock('iaa_eiac', iaa_eiac, 
c     &                elo, ehi, ilo, ihi, alo, ahi, clo, chi)
c      call writeblock('t2aa_aeik', t2aa_aeik, 
c     &                alo, ahi, elo, ehi, ilo, ihi, klo, khi)
c      call writeblock('iaa_ejcb', iaa_ejcb, 
c     &                elo, ehi, jlo, jhi, clo, chi, blo, bhi)
c      call writeblock('t2aa_beik', t2aa_beik, 
c     &                blo, bhi, elo, ehi, ilo, ihi, klo, khi)
c      call writeblock('iaa_ejac', iaa_ejac, 
c     &                elo, ehi, jlo, jhi, alo, ahi, clo, chi)
      call writeblock('t2aa_mkab', t2aa_mkab, 
     &                 mlo, mhi, klo, khi, alo, ahi, blo, bhi)
c      call writeblock('iaa_cmij', iaa_cmij, 
c     &                clo, chi, mlo, mhi, ilo, ihi, jlo, jhi)
c      call writeblock('t2aa_abmi', t2aa_abmi, 
c     &                alo, ahi, blo, bhi, mlo, mhi, ilo, ihi)
c      call writeblock('iaa_cmkj', iaa_cmkj, 
c     &                clo, chi, mlo, mhi, klo, khi, jlo, jhi)
c      call writeblock('t2aa_abmj', t2aa_abmj, 
c     &                alo, ahi, blo, bhi, mlo, mhi, jlo, jhi)
c      call writeblock('iaa_cmik', iaa_cmik, 
c     &                clo, chi, mlo, mhi, ilo, ihi, klo, khi)
c      call writeblock('t2aa_cbmk', t2aa_cbmk, 
c     &                clo, chi, blo, bhi, mlo, mhi, klo, khi)
c      call writeblock('iaa_amij', iaa_amij, 
c     &                alo, ahi, mlo, mhi, ilo, ihi, jlo, jhi)
c      call writeblock('t2aa_acmk', t2aa_acmk, 
c     &                alo, ahi, clo, chi, mlo, mhi, klo, khi)
c      call writeblock('iaa_bmij', iaa_bmij, 
c     &                blo, bhi, mlo, mhi, ilo, ihi, jlo, jhi)
c      call writeblock('t2aa_cbmi', t2aa_cbmi, 
c     &                clo, chi, blo, bhi, mlo, mhi, ilo, ihi)
c      call writeblock('iaa_amkj', iaa_amkj, 
c     &                alo, ahi, mlo, mhi, klo, khi, jlo, jhi)
c      call writeblock('t2aa_acmi', t2aa_acmi, 
c     &                alo, ahi, clo, chi, mlo, mhi, ilo, ihi)
c      call writeblock('iaa_bmkj', iaa_bmkj, 
c     &                blo, bhi, mlo, mhi, klo, khi, jlo, jhi)
c      call writeblock('t2aa_cbmj', t2aa_cbmj, 
c     &                clo, chi, blo, bhi, mlo, mhi, jlo, jhi)
c      call writeblock('iaa_amik', iaa_amik, 
c     &                alo, ahi, mlo, mhi, ilo, ihi, klo, khi)
c      call writeblock('t2aa_acmj', t2aa_acmj, 
c     &                alo, ahi, clo, chi, mlo, mhi, jlo, jhi)
c      call writeblock('iaa_bmik', iaa_bmik, 
c     &                blo, bhi, mlo, mhi, ilo, ihi, klo, khi)
c      call writeblock('iaa_ijab', iaa_ijab, 
c     &                ilo, ihi, jlo, jhi, alo, ahi, blo, bhi)
c      call writeblock('iaa_kjab', iaa_kjab, 
c     &                klo, khi, jlo, jhi, alo, ahi, blo, bhi)
c      call writeblock('iaa_ikab', iaa_ikab, 
c     &                ilo, ihi, klo, khi, alo, ahi, blo, bhi)
c      call writeblock('iaa_ijcb', iaa_ijcb, 
c     &                ilo, ihi, jlo, jhi, clo, chi, blo, bhi)
c      call writeblock('iaa_ijac', iaa_ijac, 
c     &                ilo, ihi, jlo, jhi, alo, ahi, clo, chi)
c      call writeblock('iaa_kjcb', iaa_kjcb, 
c     &                klo, khi, jlo, jhi, clo, chi, blo, bhi)
c      call writeblock('iaa_kjac', iaa_kjac, 
c     &                klo, khi, jlo, jhi, alo, ahi, clo, chi)
c      call writeblock('iaa_ikcb', iaa_ikcb, 
c     &                ilo, ihi, klo, khi, clo, chi, blo, bhi)
c      call writeblock('iaa_ikac', iaa_ikac, 
c     &                ilo, ihi, klo, khi, alo, ahi, clo, chi)
c
      return
      end
      subroutine uccsdtest_triples_pure_blocked_energy
     &  (iaa, t2aa, t1a, evals, energy,
     &   t2aa_ecji, iaa_ekab, t2aa_ecjk, iaa_eiab,
     &   t2aa_ecki, iaa_ejab, t2aa_eaji, iaa_ekcb,
     &   t2aa_ebji, iaa_ekac, t2aa_eajk, iaa_eicb,
     &   t2aa_ebjk, iaa_eiac, t2aa_eaki, iaa_ejcb,
     &   t2aa_ebki, iaa_ejac, t2aa_mkab, iaa_mcji,
     &   t2aa_miab, iaa_mcjk, t2aa_mjab, iaa_mcki,
     &   t2aa_mkcb, iaa_maji, t2aa_mkac, iaa_mbji,
     &   t2aa_micb, iaa_majk, t2aa_miac, iaa_mbjk,
     &   t2aa_mjcb, iaa_maki, t2aa_mjac, iaa_mbki,
     &   iaa_ijab,  iaa_kjab, iaa_ikab,  iaa_ijcb,
     &   iaa_ijac,  iaa_kjcb, iaa_kjac,  iaa_ikcb,
     &   iaa_ikac,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   nv_dim, no_dim, spin)
c
      implicit none
#include "cuccsdtP.fh"
      integer          alo, ahi,
     &                 blo, bhi,
     &                 clo, chi,
     &                 ilo, ihi,
     &                 jlo, jhi,
     &                 klo, khi,
     &                 elo, ehi,
     &                 mlo, mhi,
     &                 nv_dim, no_dim, spin
c
      double precision iaa(nmo, nmo, nmo, nmo)
c
      double precision t1a(nv_dim, no_dim)
      double precision t2aa(nv_dim, nv_dim, no_dim, no_dim)
      double precision evals(nmo)
c
      double precision t2aa_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision iaa_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      double precision t2aa_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision iaa_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      double precision t2aa_ecki(elo:ehi,clo:chi,klo:khi,ilo:ihi)
      double precision iaa_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      double precision t2aa_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision iaa_ekcb(elo:ehi,klo:khi,clo:chi,blo:bhi)
      double precision t2aa_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision iaa_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      double precision t2aa_eajk(elo:ehi,alo:ahi,jlo:jhi,klo:khi)
      double precision iaa_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision t2aa_ebjk(elo:ehi,blo:bhi,jlo:jhi,klo:khi)
      double precision iaa_eiac(elo:ehi,ilo:ihi,alo:ahi,clo:chi)
      double precision t2aa_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision iaa_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      double precision t2aa_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision iaa_ejac(elo:ehi,jlo:jhi,alo:ahi,clo:chi)
      double precision t2aa_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision iaa_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision t2aa_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision iaa_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      double precision t2aa_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_mcki(mlo:mhi,clo:chi,klo:khi,ilo:ihi)
      double precision t2aa_mkcb(mlo:mhi,klo:khi,clo:chi,blo:bhi)
      double precision iaa_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision t2aa_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision iaa_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision t2aa_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision iaa_majk(mlo:mhi,alo:ahi,jlo:jhi,klo:khi)
      double precision t2aa_miac(mlo:mhi,ilo:ihi,alo:ahi,clo:chi)
      double precision iaa_mbjk(mlo:mhi,blo:bhi,jlo:jhi,klo:khi)
      double precision t2aa_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      double precision t2aa_mjac(mlo:mhi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      double precision iaa_ijab(ilo:ihi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_kjab(klo:khi,jlo:jhi,alo:ahi,blo:bhi)
      double precision iaa_ikab(ilo:ihi,klo:khi,alo:ahi,blo:bhi)
      double precision iaa_ijcb(ilo:ihi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_ijac(ilo:ihi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_kjcb(klo:khi,jlo:jhi,clo:chi,blo:bhi)
      double precision iaa_kjac(klo:khi,jlo:jhi,alo:ahi,clo:chi)
      double precision iaa_ikcb(ilo:ihi,klo:khi,clo:chi,blo:bhi)
      double precision iaa_ikac(ilo:ihi,klo:khi,alo:ahi,clo:chi)
      integer i, j, k, m, a, b, c, e
      double precision w, v, d
      double precision energy
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            do c = clo, min(chi,b-1)
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     do k = klo, min(khi,j-1)
                        w = 0.0d0
                        do e = no(spin)+1, no(spin)+nv(spin)
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w = w 
     &                           + t2aa_ecji(e,c,j,i)*
     &                             iaa_ekab(e,k,a,b)
     &                           - t2aa_ecjk(e,c,j,k)*
     &                             iaa_eiab(e,i,a,b)
     &                           - t2aa_ecki(e,c,k,i)*
     &                             iaa_ejab(e,j,a,b)
     &                           - t2aa_eaji(e,a,j,i)*
     &                             iaa_ekcb(e,k,c,b)
     &                           - t2aa_ebji(e,b,j,i)*
     &                             iaa_ekac(e,k,a,c)
     &                           + t2aa_eajk(e,a,j,k)*
     &                             iaa_eicb(e,i,c,b)
     &                           + t2aa_ebjk(e,b,j,k)*
     &                             iaa_eiac(e,i,a,c)
     &                           + t2aa_eaki(e,a,k,i)*
     &                             iaa_ejcb(e,j,c,b)
     &                           + t2aa_ebki(e,b,k,i)*
     &                             iaa_ejac(e,j,a,c)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e): ',
c     &                              a,b,c,i,j,k,w
c
c  later not 1 but number of orbitals in symmetry block                       
c
                        do m = 1, no(spin)
c
                           w = w 
     &                           - t2aa_mkab(m,k,a,b)*
     &                             iaa_mcji(m,c,j,i)
     &                           + t2aa_miab(m,i,a,b)*
     &                             iaa_mcjk(m,c,j,k)
     &                           + t2aa_mjab(m,j,a,b)*
     &                             iaa_mcki(m,c,k,i)
     &                           + t2aa_mkcb(m,k,c,b)*
     &                             iaa_maji(m,a,j,i)
     &                           + t2aa_mkac(m,k,a,c)*
     &                             iaa_mbji(m,b,j,i)
     &                           - t2aa_micb(m,i,c,b)*
     &                             iaa_majk(m,a,j,k)
     &                           - t2aa_miac(m,i,a,c)*
     &                             iaa_mbjk(m,b,j,k)
     &                           - t2aa_mjcb(m,j,c,b)*
     &                             iaa_maki(m,a,k,i)
     &                           - t2aa_mjac(m,j,a,c)*
     &                             iaa_mbki(m,b,k,i)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e+m): ',
c     &                              a,b,c,i,j,k,w
                        v = w 
     &                        + t1a(c-no(spin),k)*
     &                          iaa_ijab(i,j,a,b)
     &                        - t1a(c-no(spin),i)*
     &                          iaa_kjab(k,j,a,b)
     &                        - t1a(c-no(spin),j)*
     &                          iaa_ikab(i,k,a,b)
     &                        - t1a(a-no(spin),k)*
     &                          iaa_ijcb(i,j,c,b)
     &                        - t1a(b-no(spin),k)*
     &                          iaa_ijac(i,j,a,c)
     &                        + t1a(a-no(spin),i)*
     &                          iaa_kjcb(k,j,c,b)
     &                        + t1a(b-no(spin),i)*
     &                          iaa_kjac(k,j,a,c)
     &                        + t1a(a-no(spin),j)*
     &                          iaa_ikcb(i,k,c,b)
     &                        + t1a(b-no(spin),j)*
     &                          iaa_ikac(i,k,a,c)
c
                        d = evals(a)+
     &                      evals(b)+
     &                      evals(c)-
     &                      evals(i)-
     &                      evals(j)-
     &                      evals(k)
c
                        energy = energy - v*w/d
c                        write(6,*)' a,b,c,i,j,k,w,v,d,energy: ',
c     &                              a,b,c,i,j,k,w,v,d,energy
c
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      write(6,*) ' E ', energy
c
c end block loops
      return
      end
