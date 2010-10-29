      logical function uccsdtest(rtdb)
*     
*     $Id$
*
      implicit none
#include "errquit.fh"
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
      integer l_x, l_df, l_delta, k_x, k_df, k_delta, nvar
      integer iter, i
      logical converged, t2_rd_test, t2_mp2
      double precision xxx, shift
c
      double precision energyaaa, energybbb, energyaab, energybba
      double precision uccsdtest_triples_pure, uccsdtest_triples_mixed
      double precision uccsdtest_triples_pure_blocked
      double precision uccsdtest_triples_mixed_blocked
      logical int_normalize
      external int_normalize
c
c     load the geometry/basis set and get info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0, RTDB_ERR)
      if (.not.int_normalize(rtdb, basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
c
c     Read the MO vectors and evals from a UHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('rjhtest: failed to read movecs header',0, DISK_ERR)
      if (nsets.ne.2) call errquit('UHF ONLY!',0, INPUT_ERR)
      nmo = nmo_file(1)
      if (nmo .ne. nbf) call errquit('NO LINEAR DEPENDNCE PERMITTED',0,
     &       INPUT_ERR)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      if (nmo_file(1).ne.nmo_file(2)) call errquit('Uh?',0, BASIS_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nalpha', mt_int, 1, nalpha))
     $     call errquit('nalpha?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nbeta', mt_int, 1, nbeta))
     $     call errquit('nbeta?',0, RTDB_ERR)

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
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'aeval',l_aeval, k_aeval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'beval',l_beval, k_beval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'amos', l_amos, k_amos))
     $     call errquit('ma amos', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'bmos', l_bmos, k_bmos))
     $     call errquit('ma bmos', nbf*nbf, MA_ERR)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_aeval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_amos), nbf)
      if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), dbl_mb(k_beval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_bmos), nbf)
c
c     Change the phase and order of some of the beta MOs so that they
c     are different from alpha even for closed shell
c
c      do i = 1, nmo, 2          ! Change phase of even MOs
c         call dscal(nbf, -1d0, dbl_mb(k_bmos+(i-1)*nbf), 1)
c      end do
c      do i = 2, nbeta, 2        ! Swap alternate beta occupied orbitals
c         call dcopy(nbf, dbl_mb(k_bmos+(i-1)*nbf), 1, dbl_mb(k_occ), 1) 
c         call dcopy(nbf, dbl_mb(k_bmos+(i-2)*nbf), 1, 
c     $        dbl_mb(k_bmos+(i-1)*nbf), 1)
c         call dcopy(nbf, dbl_mb(k_occ), 1, dbl_mb(k_bmos+(i-2)*nbf), 1)
c         xxx = dbl_mb(k_beval+i-1)
c         dbl_mb(k_beval+i-1) = dbl_mb(k_beval+i-2)
c         dbl_mb(k_beval+i-2) = xxx
c      end do
c      do i = nbeta+2, nmo, 2        ! Swap alternate beta virtual orbitals
c         call dcopy(nbf, dbl_mb(k_bmos+(i-1)*nbf), 1, dbl_mb(k_occ), 1) 
c         call dcopy(nbf, dbl_mb(k_bmos+(i-2)*nbf), 1, 
c     $        dbl_mb(k_bmos+(i-1)*nbf), 1)
c         call dcopy(nbf, dbl_mb(k_occ), 1, dbl_mb(k_bmos+(i-2)*nbf), 1)
c         xxx = dbl_mb(k_beval+i-1)
c         dbl_mb(k_beval+i-1) = dbl_mb(k_beval+i-2)
c         dbl_mb(k_beval+i-2) = xxx
c      end do
c
      write(6,*) ' Alpha eigenvalues '
      call output(dbl_mb(k_aeval),1,nmo,1,1,nmo,1,1)
      write(6,*) ' Beta  eigenvalues '
      call output(dbl_mb(k_beval),1,nmo,1,1,nmo,1,1)
c
      if (.not. rtdb_get(rtdb, 'uccsdtest:shift', mt_dbl, 1, shift))then
c
c     Shift is added to the denominators in the preconditioner
c     ... default is to make the smallest denominator 0.3
c
         shift = 2d0 * min(
     $        dbl_mb(k_aeval+nalpha)-dbl_mb(k_aeval+nalpha-1),
     $        dbl_mb(k_beval+nbeta)-dbl_mb(k_beval+nbeta-1))
         shift = 0.3d0 - shift
      end if
      write(6,*) ' Preconditioner denominator shift ', shift
c$$$      write(6,*) ' Alpha MOs'
c$$$      call output(dbl_mb(k_amos),1,nbf,1,nmo,nbf,nmo,1)
c$$$      write(6,*) ' Beta  MOs'
c$$$      call output(dbl_mb(k_bmos),1,nbf,1,nmo,nbf,nmo,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
c
      call util_inplace_transpose(nmo, dbl_mb(k_amos))
      call util_inplace_transpose(nmo, dbl_mb(k_bmos))
c
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_aa, k_iaa))
     $     call errquit('ma aa moints', nbf**4, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_ab, k_iab))
     $     call errquit('ma ab moints', nbf**4, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_ba, k_iba))
     $     call errquit('ma ba moints', nbf**4, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf**4,'moints', l_bb, k_ibb))
     $     call errquit('ma bb moints', nbf**4, MA_ERR)
c
c     Canonical UHF MP2 as a test and to get starting amplitudes.
c
      if (.not. ma_push_get(mt_dbl, (noa*nva), 't1', l_t1a, k_t1a))
     $     call errquit('ma t1a', (noa*nva), MA_ERR)
      if (.not. ma_push_get(mt_dbl, (nob*nvb), 't1', l_t1b, k_t1b))
     $     call errquit('ma t1b', (nob*nvb), MA_ERR)
      if (.not. ma_push_get(mt_dbl, (noa*nva)**2, 't2', l_t2aa, k_t2aa))
     $     call errquit('ma t2_aa', (noa*nva)**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, (nob*nvb)**2, 't2', l_t2bb, k_t2bb))
     $     call errquit('ma t2_bb', (nob*nvb)**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     't2', l_t2ab, k_t2ab))
     $     call errquit('ma t2_ab', (noa*nva)*(nob*nvb), MA_ERR)
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
c$$$      call dfill((noa*nva)**2,0.0d0,dbl_mb(k_t2aa),1)
c$$$      call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_t2bb),1)
c$$$      call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_t2ab),1)
c$$$      call dfill((noa*nva), 0.0d0, dbl_mb(k_t1a), 1)
c$$$      call dfill((nob*nvb), 0.0d0, dbl_mb(k_t1b), 1)
c
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fa', l_fa, k_fa))
     $     call errquit('ma fa', nmo**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fb', l_fb, k_fb))
     $     call errquit('ma fb', nmo**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fa', l_ha, k_ha))
     $     call errquit('ma ha', nmo**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo**2, 'fb', l_hb, k_hb))
     $     call errquit('ma hb', nmo**2, MA_ERR)
c$$$      call dfill(nmo*nmo, 0.0d0, dbl_mb(k_fa), 1)
c$$$      call dfill(nmo*nmo, 0.0d0, dbl_mb(k_fb), 1)
c$$$      call dcopy(nmo, dbl_mb(k_aeval), 1, dbl_mb(k_fa), nmo+1)
c$$$      call dcopy(nmo, dbl_mb(k_beval), 1, dbl_mb(k_fb), nmo+1)
c
c     Space for the residual vectors
c
      if (.not. ma_push_get(mt_dbl, (noa*nva), 'r1', l_r1a, k_r1a))
     $     call errquit('ma r1a', (noa*nva), MA_ERR)
      if (.not. ma_push_get(mt_dbl, (nob*nvb), 'r1', l_r1b, k_r1b))
     $     call errquit('ma r1b', (nob*nvb), MA_ERR)
      if (.not. ma_push_get(mt_dbl, (noa*nva)**2, 'r2', l_r2aa, k_r2aa))
     $     call errquit('ma r2_aa', (noa*nva)**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, (nob*nvb)**2, 'r2', l_r2bb, k_r2bb))
     $     call errquit('ma r2_bb', (nob*nvb)**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     'r2', l_r2ab, k_r2ab))
     $     call errquit('ma r2_ab', (noa*nva)*(nob*nvb), MA_ERR)
c
c     Temporary workspace for transposed t2ab and
c     particle & hole transformation matrices and workspace
c
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     't2', l_t2ba, k_t2ba))
     $     call errquit('ma t2_ba', (noa*nva)*(nob*nvb), MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'a lamda p', l_alp, k_alp))
     $     call errquit('ma alp', nbf*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'a lamda h', l_alh, k_alh))
     $     call errquit('ma alh', nbf*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'b lamda p', l_blp, k_blp))
     $     call errquit('ma alp', nmo*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'b lamda h', l_blh, k_blh))
     $     call errquit('ma alh', nbf*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo*nmo, 'work', l_work, k_work))
     $     call errquit('ma work', nmo*nmo, MA_ERR)
c
c     For the n^6 algorithm we need some more intermediates
c
      if (.not. ma_push_get(mt_dbl, noa*nva*noa*nva, 'z1', l_z1, k_z1))
     $     call errquit(' ma z1 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*nva*noa*nva, 'z2', l_z2, k_z2))
     $     call errquit(' ma z2 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nob*nvb*noa*nva, 'z3', l_z3, k_z3))
     $     call errquit(' ma z3 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nob*nvb*nob*nvb, 'z4', l_z4, k_z4))
     $     call errquit(' ma z4 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nob*nvb*nob*nvb, 'z5', l_z5, k_z5))
     $     call errquit(' ma z5 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*nva*nob*nvb, 'z6', l_z6, k_z6))
     $     call errquit(' ma z6 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*nva*nob*nvb, 'z7', l_z7, k_z7))
     $     call errquit(' ma z7 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*nvb*noa*nvb, 'z8', l_z8, k_z8))
     $     call errquit(' ma z8 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*noa*noa*noa, 'x1', l_x1, k_x1))
     $     call errquit(' ma x1 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nob*nob*nob*nob, 'x2', l_x2, k_x2))
     $     call errquit(' ma x2 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, noa*nob*noa*nob, 'x3', l_x3, k_x3))
     $     call errquit(' ma x3 ', 0, MA_ERR)
c
c     Space for the new solver
c
      nvar = noa*noa*nva*nva + nob*nob*nvb*nvb + noa*nob*nva*nvb + 
     $     noa*nva + nob*nvb
      if (.not. ma_push_get(mt_dbl, nvar*10, 'x', l_x, k_x))
     $     call errquit(' ma x4 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar*10, 'x', l_df, k_df))
     $     call errquit(' ma x5 ', 0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar, 'x', l_delta, k_delta))
     $     call errquit(' ma x6 ', 0, MA_ERR)
c
c     if t1 and t2 on disk jump and read
c
      if (.not. rtdb_get(rtdb, 'uccsdt:t2_rd_test', mt_log, 
     &   1, t2_rd_test)) t2_rd_test = .false.
      if(t2_rd_test)goto 3768
c
c     Iterate
c
c$$$      call dfill((noa*nva)**2,0.0d0,dbl_mb(k_t2aa),1)
c$$$      call dfill((nob*nvb)**2,0.0d0,dbl_mb(k_t2bb),1)
c$$$      call dfill((noa*nva)*(nob*nvb),0.0d0,dbl_mb(k_t2ab),1)
c$$$      call dfill((noa*nva), 0.0d0, dbl_mb(k_t1a), 1)
c$$$      call dfill((nob*nvb), 0.0d0, dbl_mb(k_t1b), 1)
c
      do iter = 1, 10           ! NOTE DIMENSIONS OF X/DF
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
c$$$         call uccsdtest_product_pure(
c$$$     $        nmo,
c$$$     $        noa, nob, nva, nvb, 
c$$$     $        dbl_mb(k_fa), dbl_mb(k_fb),
c$$$     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
c$$$     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
c$$$     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
c$$$     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
c$$$     $        dbl_mb(k_r2aa))
c$$$c
c$$$c     Beta pure spin ... use the transposed mixed spin amplitudes
c$$$c     and integrals and then call the pure spin routine with spins flipped
c$$$
c$$$         call uccsdtest_transpose_tab(noa, nob, nva, nvb, 
c$$$     $        dbl_mb(k_t2ab), dbl_mb(k_t2ba))
c$$$         call uccsdtest_product_pure(
c$$$     $        nmo,
c$$$     $        nob, noa, nvb, nva, 
c$$$     $        dbl_mb(k_fb), dbl_mb(k_fa),
c$$$     $        dbl_mb(k_ibb), dbl_mb(k_iaa), dbl_mb(k_iba),
c$$$     $        dbl_mb(k_t1b), dbl_mb(k_t1a),
c$$$     $        dbl_mb(k_t2bb),dbl_mb(k_t2aa),dbl_mb(k_t2ba),
c$$$     $        dbl_mb(k_r1b), dbl_mb(k_r1a),
c$$$     $        dbl_mb(k_r2bb))
c$$$c
c$$$c     Mixed
c$$$c
c$$$         call uccsdtest_product_mixed(
c$$$     $        nmo,
c$$$     $        noa, nob, nva, nvb, 
c$$$     $        dbl_mb(k_fa), dbl_mb(k_fb),
c$$$     $        dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
c$$$     $        dbl_mb(k_t1a), dbl_mb(k_t1b),
c$$$     $        dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
c$$$     $        dbl_mb(k_r1a), dbl_mb(k_r1b),
c$$$     $        dbl_mb(k_r2ab))
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
     $        dbl_mb(k_aeval), dbl_mb(k_beval), converged,
     $        dbl_mb(k_x), dbl_mb(k_df), dbl_mb(k_delta),
     $        nvar, iter, shift)
c
c$$$         call jan_debug_print('T1a',dbl_mb(k_t1a),  nva, noa, 1, 1)
c$$$         call jan_debug_print('T1b',dbl_mb(k_t1b),  nvb, nob, 1, 1)
c$$$         call jan_debug_print('Taa',dbl_mb(k_t2aa), nva, nva, noa, noa)
c$$$         call jan_debug_print('Tbb',dbl_mb(k_t2bb), nvb, nvb, nob, nob)
c$$$         call jan_debug_print('Tab',dbl_mb(k_t2ab), nva, nvb, noa, nob)
c
         if (converged) goto 7764
      end do
 7764 continue
c     
      write(6,*) ' CONVERGED amplitudes'
      call jan_debug_print('T1a',dbl_mb(k_t1a),  nva, noa, 1, 1)
      call jan_debug_print('T1b',dbl_mb(k_t1b),  nvb, nob, 1, 1)
      call jan_debug_print('Taa',dbl_mb(k_t2aa), nva, nva, noa, noa)
      call jan_debug_print('Tbb',dbl_mb(k_t2bb), nvb, nvb, nob, nob)
      call jan_debug_print('Tab',dbl_mb(k_t2ab), nva, nvb, noa, nob)
c
c     write converged t1 and t2 to disk
c
      call write_t(
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
c
 3768 continue
c
      call read_t(
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
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
c  test amplitudes with mp2 t2's.
c
      if (.not. rtdb_get(rtdb, 'uccsdt:t2_mp2', mt_log, 
     &   1, t2_mp2)) t2_mp2 = .false.
      if (t2_mp2)then
         write(6,*)' ************************************** '
         write(6,*)' testing triples with mp2 t2 amplitudes '
         write(6,*)' ************************************** '
         call uccsdtest_mp2(
     $     nmo,
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_aeval), dbl_mb(k_beval),
     $     dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
c     
         write(6,*) ' MP2 amplitudes'
         call jan_debug_print('T1a',dbl_mb(k_t1a),  nva, noa, 1, 1)
         call jan_debug_print('T1b',dbl_mb(k_t1b),  nvb, nob, 1, 1)
         call jan_debug_print('Taa',dbl_mb(k_t2aa), nva, nva, noa, noa)
         call jan_debug_print('Tbb',dbl_mb(k_t2bb), nvb, nvb, nob, nob)
         call jan_debug_print('Tab',dbl_mb(k_t2ab), nva, nvb, noa, nob)
c
c        write mp2 t1 and t2 to disk
c
         call write_t(
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab))
c
         write(6,*)' t1 and t2 amplitudes reset '
c
      endif
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
      write(6,*) ' ORIGINAL TEST CODE (T) '
      write(6,*) energyaaa, energybbb
      write(6,*) energyaab, energybba
      write(6,*) energyaaa+energybbb+energyaab+energybba
c
      call uccsdt_initialize(nmo, noa, nob, nva, nvb, 
     $     dbl_mb(k_aeval), dbl_mb(k_beval),
     $     k_iaa, k_ibb, k_iab,
     $     k_t1a, k_t1b,
     $     k_t2aa,k_t2bb,k_t2ab)
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
      energyaab = uccsdtest_triples_mixed_blocked(
     $     noa, nob, nva, nvb, 
     $     dbl_mb(k_iaa), dbl_mb(k_ibb), dbl_mb(k_iab),
     $     dbl_mb(k_t1a), dbl_mb(k_t1b),
     $     dbl_mb(k_t2aa),dbl_mb(k_t2bb),dbl_mb(k_t2ab),
     $     dbl_mb(k_aeval), dbl_mb(k_beval), 1)
c
      write(6,*) energyaab
c
      if (.not. ma_chop_stack(l_occ)) call errquit('ma?',0, MA_ERR)
c
      uccsdtest = .true.
c
      if (.not. bas_destroy(basis)) call errquit('basis',0, BASIS_ERR)
      if (.not. geom_destroy(geom)) call errquit('geom',0, GEOM_ERR)
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
     $     r2aa,r2bb,r2ab, ea, eb, converged,
     $     x, df, delta, nvar, iter, shift)
      implicit none
#include "errquit.fh"
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
      integer nvar, iter           ! #variables, iteration number
      double precision x(nvar,*), df(nvar,*), delta(nvar)
      double precision shift
      logical converged
c
c     t1 = t1 + R(e,m)/(em-ee)
c     t2 = t2 + R(e,f,m,n)/(em+en-ee-ef)
c     
      integer m, n, e, f, ind
      double precision eaa, ebb, eab, energy
      double precision delt, r1norm, r2norm 
c
      r1norm = 0.0d0
      r2norm = 0.0d0
c
*      call dfill(noa*nva, 0.0d0, r1a, 1)
*      call dfill(nob*nvb, 0.0d0, r1b, 1)
c
c     Singles
c
      ind = 1
      do m = 1, noa
         do e = 1, nva
            x(ind,iter) = t1a(e,m)
            df(ind,iter) = -r1a(e,m) / (ea(m) - ea(e+noa) - shift*0.5d0)
            ind = ind + 1
         end do
      end do
      do m = 1, nob
         do e = 1, nvb
            x(ind,iter) = t1b(e,m)
            df(ind,iter) = -r1b(e,m) / (eb(m) - eb(e+nob) - shift*0.5d0)
            ind = ind + 1
         end do
      end do
c
c     Pure alpha
c
      do n = 1, noa
         do m = 1, noa
            do f = 1, nva
               do e = 1, nva
                  x(ind,iter) = t2aa(e,f,m,n)
                  df(ind,iter) = -r2aa(e,f,m,n) / 
     $                 (ea(m)+ea(n)-ea(noa+e)-ea(noa+f) - shift)
                  ind = ind + 1
               end do
            end do
         end do
      end do
c
c     Pure beta
c
      do n = 1, nob
         do m = 1, nob
            do f = 1, nvb
               do e = 1, nvb
                  x(ind,iter) = t2bb(e,f,m,n)
                  df(ind,iter) = -r2bb(e,f,m,n) / 
     $                 (eb(m)+eb(n)-eb(nob+e)-eb(nob+f) - shift)
                  ind = ind + 1
               end do
            end do
         end do
      end do
c
c     Mixed
c
      do n = 1, nob
         do m = 1, noa
            do f = 1, nvb
               do e = 1, nva
                  x(ind,iter) = t2ab(e,f,m,n)
                  df(ind,iter) = -r2ab(e,f,m,n) /
     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f) - shift)
                  ind = ind + 1
               end do
            end do
         end do
      end do
      if (ind .ne. nvar+1) call errquit('nvar ????? ', ind, UNKNOWN_ERR)
c
      call uccsdt_solver(nvar,iter,x,df,delta)
c
      ind = 1
      do m = 1, noa
         do e = 1, nva
            delt = delta(ind)
            ind = ind + 1
            r1norm = r1norm + delt**2
            t1a(e,m) = t1a(e,m) + delt
         end do
      end do
      do m = 1, nob
         do e = 1, nvb
            delt = delta(ind)
            ind = ind + 1
            r1norm = r1norm + delt**2
            t1b(e,m) = t1b(e,m) + delt
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
                  delt = delta(ind)
                  ind = ind + 1
                  r2norm = r2norm + delt*delt
                  t2aa(e,f,m,n) = t2aa(e,f,m,n) + delt
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
                  delt = delta(ind)
                  ind = ind + 1
                  r2norm = r2norm + delt*delt
                  t2bb(e,f,m,n) = t2bb(e,f,m,n) + delt
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
                  delt = delta(ind)
                  ind = ind + 1
*                  write(6,*) ' denom ', 
*     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
                  r2norm = r2norm + delt*delt
                  t2ab(e,f,m,n) = t2ab(e,f,m,n) + delt
                  eab = eab + (
     $                 t2ab(e,f,m,n) + 
     $                 t1a(e,m)*t1b(f,n))*
     $                 iab(m,n,noa+e,nob+f)
               end do
            end do
         end do
      end do
      if (ind .ne. nvar+1) call errquit('nvar 2????? ', ind,
     &       UNKNOWN_ERR)
c
c$$$      do m = 1, noa
c$$$         do e = 1, nva
c$$$            delt = r1a(e,m) / (ea(m) - ea(e+noa))
c$$$            r1norm = r1norm + delt**2
c$$$            t1a(e,m) = t1a(e,m) + delt
c$$$         end do
c$$$      end do
c$$$      do m = 1, nob
c$$$         do e = 1, nvb
c$$$            delt = r1b(e,m) / (eb(m) - eb(e+nob))
c$$$            r1norm = r1norm + delt**2
c$$$            t1b(e,m) = t1b(e,m) + delt
c$$$         end do
c$$$      end do
c$$$c
c$$$c     Pure alpha
c$$$c
c$$$      eaa = 0.0d0
c$$$      do n = 1, noa
c$$$         do m = 1, noa
c$$$            do f = 1, nva
c$$$               do e = 1, nva
c$$$                  delt = r2aa(e,f,m,n) / 
c$$$     $                 (ea(m)+ea(n)-ea(noa+e)-ea(noa+f))
c$$$                  r2norm = r2norm + delt*delt
c$$$                  t2aa(e,f,m,n) = t2aa(e,f,m,n) + delt
c$$$                  eaa = eaa + (
c$$$     $                 t2aa(e,f,m,n) + 
c$$$     $                 t1a(e,m)*t1a(f,n) - t1a(e,n)*t1a(f,m))*
c$$$     $                 iaa(m,n,noa+e,noa+f)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$c     Pure beta
c$$$c
c$$$      ebb = 0.0d0
c$$$      do n = 1, nob
c$$$         do m = 1, nob
c$$$            do f = 1, nvb
c$$$               do e = 1, nvb
c$$$                  delt = r2bb(e,f,m,n) /
c$$$     $                 (eb(m)+eb(n)-eb(nob+e)-eb(nob+f))
c$$$                  r2norm = r2norm + delt*delt
c$$$                  t2bb(e,f,m,n) = t2bb(e,f,m,n) + delt
c$$$                  ebb = ebb + (
c$$$     $                 t2bb(e,f,m,n) + 
c$$$     $                 t1b(e,m)*t1b(f,n) - t1b(e,n)*t1b(f,m))*
c$$$     $                 ibb(m,n,nob+e,nob+f)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$c     Mixed
c$$$c
c$$$      eab = 0.0d0
c$$$      do n = 1, nob
c$$$         do m = 1, noa
c$$$            do f = 1, nvb
c$$$               do e = 1, nva
c$$$                  delt = r2ab(e,f,m,n) /
c$$$     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
c$$$*                  write(6,*) ' denom ', 
c$$$*     $                 (ea(m)+eb(n)-ea(noa+e)-eb(nob+f))
c$$$                  r2norm = r2norm + delt*delt
c$$$                  t2ab(e,f,m,n) = t2ab(e,f,m,n) + delt
c$$$                  eab = eab + (
c$$$     $                 t2ab(e,f,m,n) + 
c$$$     $                 t1a(e,m)*t1b(f,n))*
c$$$     $                 iab(m,n,noa+e,nob+f)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c
      energy = eaa/4.0d0 + ebb/4.0d0 + eab
c     
      r1norm = sqrt(r1norm)
      r2norm = sqrt(r2norm)
      write(6,1) iter, energy, r1norm, r2norm
 1    format(i4,' Energy = ', f20.8,'  update norms = ',1p,2d9.2)
c
      converged = max(r1norm,r2norm) .lt. 1e-8
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
#include "errquit.fh"
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
     $     call errquit('x',1, MA_ERR)
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
     $     call errquit('x',2, MA_ERR)
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
     $     call errquit('x',4, MA_ERR)
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
     $     call errquit('x',5, MA_ERR)
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
     $                       (ibb(i,f+nob,a+nob,n) + z4(i,a,n,f)
     $                       +z5(i,a,n,f))
c    $                       (ibb(i,f+nob,a+nob,n) + z4(i,a,n,f))
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
     $                       iaa(i,e+noa,a+noa,m)
c    $                       (iaa(i,e+noa,a+noa,m) + z2(i,a,m,e))
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
     $     call errquit('x',6, MA_ERR)
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
#include "errquit.fh"
#include "cuccsdtP.fh"
      integer spin, spinck
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
      double precision energy, energy2
c
c     temp buffers for t2 and ints 
c
c      double precision t2_ecji(5000), int_ekab(5000),
c     &                 t2_ecjk(5000), int_eiab(5000),
c     &                 t2_ecik(5000), int_ejab(5000),
c     &                 t2_eaji(5000), int_ekbc(5000),
c     &                 t2_ebji(5000), int_ekac(5000),
c     &                 t2_eakj(5000), int_eicb(5000),
c     &                 t2_ebkj(5000), int_eica(5000),
c     &                 t2_eaki(5000), int_ejcb(5000),
c     &                 t2_ebki(5000), int_ejca(5000)
cc
c      double precision t2_mkab(5000), int_mcji(5000),
c     &                 t2_miab(5000), int_mcjk(5000),
c     &                 t2_mjab(5000), int_mcik(5000),
c     &                 t2_mkbc(5000), int_maji(5000),
c     &                 t2_mkac(5000), int_mbji(5000),
c     &                 t2_micb(5000), int_makj(5000),
c     &                 t2_mica(5000), int_mbkj(5000),
c     &                 t2_mjcb(5000), int_maki(5000),
c     &                 t2_mjca(5000), int_mbki(5000)
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
c      double precision int_buf(5000), t_buf(5000), w(100000), v(100000)
      double precision int_buf(5000), t_buf(5000), w(100000)
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
c      gs(i,j,k,a,b,c)   = t1a(c-no(spin),k)*iaa(a,b,i,j)
c
      call output(evals,1,nmo,1,1,nmo,1,1)
c
c
c     note c&k have opposite spins of a&i
c
      if (spin.eq.1)then
         spinck=2
      elseif (spin.eq.2)then
         spinck=1
      else
         call errquit(' unrecognized spin ',spin, UNKNOWN_ERR)
      endif
      energy = 0.0d0
      energy2 = 0.0d0
c
c     debug: temp zero t1a
c
c      call dfill (nv_dim*no_dim, 0.0d0, t1a, 1)
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
c                        call dfill (5000, 0.0d0, t2_ecji, 1)
c                        call dfill (5000, 0.0d0, int_ekab, 1)
c                        call dfill (5000, 0.0d0, t2_ecjk, 1)
c                        call dfill (5000, 0.0d0, int_eiab, 1)
c                        call dfill (5000, 0.0d0, t2_ecik, 1)
c                        call dfill (5000, 0.0d0, int_ejab, 1)
c                        call dfill (5000, 0.0d0, t2_eaji, 1)
c                        call dfill (5000, 0.0d0, int_ekbc, 1)
c                        call dfill (5000, 0.0d0, t2_ebji, 1)
c                        call dfill (5000, 0.0d0, int_ekac, 1)
c                        call dfill (5000, 0.0d0, t2_eakj, 1)
c                        call dfill (5000, 0.0d0, int_eicb, 1)
c                        call dfill (5000, 0.0d0, t2_ebkj, 1)
c                        call dfill (5000, 0.0d0, int_eica, 1)
c                        call dfill (5000, 0.0d0, t2_eaki, 1)
c                        call dfill (5000, 0.0d0, int_ejcb, 1)
c                        call dfill (5000, 0.0d0, t2_ebki, 1)
c                        call dfill (5000, 0.0d0, int_ejca, 1)
cc
c                        call dfill (5000, 0.0d0, t2_mkab, 1)
c                        call dfill (5000, 0.0d0, int_mcji, 1)
c                        call dfill (5000, 0.0d0, t2_miab, 1)
c                        call dfill (5000, 0.0d0, int_mcjk, 1)
c                        call dfill (5000, 0.0d0, t2_mjab, 1)
c                        call dfill (5000, 0.0d0, int_mcik, 1)
c                        call dfill (5000, 0.0d0, t2_mkbc, 1)
c                        call dfill (5000, 0.0d0, int_maji, 1)
c                        call dfill (5000, 0.0d0, t2_mkac, 1)
c                        call dfill (5000, 0.0d0, int_mbji, 1)
c                        call dfill (5000, 0.0d0, t2_micb, 1)
c                        call dfill (5000, 0.0d0, int_makj, 1)
c                        call dfill (5000, 0.0d0, t2_mica, 1)
c                        call dfill (5000, 0.0d0, int_mbkj, 1)
c                        call dfill (5000, 0.0d0, t2_mjcb, 1)
c                        call dfill (5000, 0.0d0, int_maki, 1)
c                        call dfill (5000, 0.0d0, t2_mjca, 1)
c                        call dfill (5000, 0.0d0, int_mbki, 1)
c
                        call dfill (5000, 0.0d0, int_abij, 1)
                        call dfill (5000, 0.0d0, int_abjk, 1)
                        call dfill (5000, 0.0d0, int_abik, 1)
                        call dfill (5000, 0.0d0, int_bcij, 1)
                        call dfill (5000, 0.0d0, int_acij, 1)
                        call dfill (5000, 0.0d0, int_bcjk, 1)
                        call dfill (5000, 0.0d0, int_acjk, 1)
                        call dfill (5000, 0.0d0, int_bcik, 1)
                        call dfill (5000, 0.0d0, int_acik, 1)
c
                        call dfill (100000, 0.0d0, w, 1)
c                        call dfill (100000, 0.0d0, v, 1)
c
c                        write(6,*)' after butt-load of dfills '
c
c                        call triples_temp_fill_buffers(iaa, t2aa,
c     &                     t2_ecji, int_ekab, t2_ecjk, int_eiab,
c     &                     t2_ecik, int_ejab, t2_eaji, int_ekbc,
c     &                     t2_ebji, int_ekac, t2_eakj, int_eicb,
c     &                     t2_ebkj, int_eica, t2_eaki, int_ejcb,
c     &                     t2_ebki, int_ejca, t2_mkab, int_mcji,
c     &                     t2_miab, int_mcjk, t2_mjab, int_mcik,
c     &                     t2_mkbc, int_maji, t2_mkac, int_mbji,
c     &                     t2_micb, int_makj, t2_mica, int_mbkj,
c     &                     t2_mjcb, int_maki, t2_mjca, int_mbki,
c     &                     int_abij,  int_abjk, int_abik,  int_bcij,
c     &                     int_acij,  int_bcjk, int_acjk,  int_bcik,
c     &                     int_acik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
c     &                     nv_dim, no_dim, spin)
c
c                        call uccsdt_get_3x(k_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, int_ekab)
c                        call uccsdt_get_3x(i_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, int_eiab)
c                        call uccsdt_get_3x(j_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, int_ejab)
c                        call uccsdt_get_3x(k_blk, b_blk, c_blk, 
c     &                     spin, spin, spin, int_ekbc)
c                        call uccsdt_get_3x(k_blk, a_blk, c_blk, 
c     &                     spin, spin, spin, int_ekac)
c                        call uccsdt_get_3x(i_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, int_eicb)
c                        call uccsdt_get_3x(i_blk, c_blk, a_blk, 
c     &                     spin, spin, spin, int_eica)
c                        call uccsdt_get_3x(j_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, int_ejcb)
c                        call uccsdt_get_3x(j_blk, c_blk, a_blk,
c     &                     spin, spin, spin, int_ejca)
c                        call uccsdt_get_2x(a_blk, b_blk, i_blk, j_blk,
c     &                     spin, spin, spin, spin, int_abij)
c                        call uccsdt_get_2x(a_blk, b_blk, j_blk, k_blk,
c     &                     spin, spin, spin, spin, int_abjk)
c                        call uccsdt_get_2x(a_blk, b_blk, i_blk, k_blk,
c     &                     spin, spin, spin, spin, int_abik)
c                        call uccsdt_get_2x(b_blk, c_blk, i_blk, j_blk,
c     &                     spin, spin, spin, spin, int_bcij)
c                        call uccsdt_get_2x(a_blk, c_blk, i_blk, j_blk,
c     &                     spin, spin, spin, spin, int_acij)
c                        call uccsdt_get_2x(b_blk, c_blk, j_blk, k_blk,
c     &                     spin, spin, spin, spin, int_bcjk)
c                        call uccsdt_get_2x(a_blk, c_blk, j_blk, k_blk,
c     &                     spin, spin, spin, spin, int_acjk)
c                        call uccsdt_get_2x(b_blk, c_blk, i_blk, k_blk,
c     &                     spin, spin, spin, spin, int_bcik)
c                        call uccsdt_get_2x(a_blk, c_blk, i_blk, k_blk, 
c     &                     spin, spin, spin, spin, int_acik)
c                        call uccsdt_get_1x(c_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, int_mcji)
c                        call uccsdt_get_1x(c_blk, j_blk, k_blk, 
c     &                     spin, spin, spin, int_mcjk)
c                        call uccsdt_get_1x(c_blk, i_blk, k_blk,
c     &                     spin, spin, spin, int_mcik)
c                        call uccsdt_get_1x(a_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, int_maji)
c                        call uccsdt_get_1x(b_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, int_mbji)
c                        call uccsdt_get_1x(a_blk, k_blk, j_blk,
c     &                     spin, spin, spin, int_makj)
c                        call uccsdt_get_1x(b_blk, k_blk, j_blk,
c     &                     spin, spin, spin, int_mbkj)
c                        call uccsdt_get_1x(a_blk, k_blk, i_blk, 
c     &                     spin, spin, spin, int_maki)
c                        call uccsdt_get_1x(b_blk, k_blk, i_blk, 
c     &                     spin, spin, spin, int_mbki)
c                        call uccsdt_get_t3x(c_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, t2_ecji)
c                        call uccsdt_get_t3x(c_blk, j_blk, k_blk, 
c     &                     spin, spin, spin, t2_ecjk)
c                        call uccsdt_get_t3x(c_blk, i_blk, k_blk, 
c     &                     spin, spin, spin, t2_ecik)
c                        call uccsdt_get_t3x(a_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, t2_eaji)
c                        call uccsdt_get_t3x(b_blk, j_blk, i_blk, 
c     &                     spin, spin, spin, t2_ebji)
c                        call uccsdt_get_t3x(a_blk, k_blk, j_blk, 
c     &                     spin, spin, spin, t2_eakj)
c                        call uccsdt_get_t3x(b_blk, k_blk, j_blk, 
c     &                     spin, spin, spin, t2_ebkj)
c                        call uccsdt_get_t3x(a_blk, k_blk, i_blk, 
c     &                     spin, spin, spin, t2_eaki)
c                        call uccsdt_get_t3x(b_blk, k_blk, i_blk, 
c     &                     spin, spin, spin, t2_ebki)
c                        call uccsdt_get_t1x(k_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2_mkab)
c                        call uccsdt_get_t1x(i_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2_miab)
c                        call uccsdt_get_t1x(j_blk, a_blk, b_blk, 
c     &                     spin, spin, spin, t2_mjab)
c                        call uccsdt_get_t1x(k_blk, b_blk, c_blk,
c     &                     spin, spin, spin, t2_mkbc)
c                        call uccsdt_get_t1x(k_blk, a_blk, c_blk, 
c     &                     spin, spin, spin, t2_mkac)
c                        call uccsdt_get_t1x(i_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, t2_micb)
c                        call uccsdt_get_t1x(i_blk, c_blk, a_blk,
c     &                     spin, spin, spin, t2_mica)
c                        call uccsdt_get_t1x(j_blk, c_blk, b_blk, 
c     &                     spin, spin, spin, t2_mjcb)
c                        call uccsdt_get_t1x(j_blk, c_blk, a_blk,
c     &                     spin, spin, spin, t2_mjca)
c
c                        call triples_temp_dump_buffers(
c     &                     t2_ecji, int_ekab, t2_ecjk, int_eiab,
c     &                     t2_ecik, int_ejab, t2_eaji, int_ekbc,
c     &                     t2_ebji, int_ekac, t2_eakj, int_eicb,
c     &                     t2_ebkj, int_eica, t2_eaki, int_ejcb,
c     &                     t2_ebki, int_ejca, t2_mkab, int_mcji,
c     &                     t2_miab, int_mcjk, t2_mjab, int_mcik,
c     &                     t2_mkbc, int_maji, t2_mkac, int_mbji,
c     &                     t2_micb, int_makj, t2_mica, int_mbkj,
c     &                     t2_mjcb, int_maki, t2_mjca, int_mbki,
c     &                     int_abij,  int_abjk, int_abik,  int_bcij,
c     &                     int_acij,  int_bcjk, int_acjk,  int_bcik,
c     &                     int_acik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi)
c
c                        call uccsdtest_triples_pure_blocked_energy_45
c     &                    (t1a, evals, energy,
c     &                     t2_ecji, int_ekab, t2_ecjk, int_eiab,
c     &                     t2_ecik, int_ejab, t2_eaji, int_ekbc,
c     &                     t2_ebji, int_ekac, t2_eakj, int_eicb,
c     &                     t2_ebkj, int_eica, t2_eaki, int_ejcb,
c     &                     t2_ebki, int_ejca, t2_mkab, int_mcji,
c     &                     t2_miab, int_mcjk, t2_mjab, int_mcik,
c     &                     t2_mkbc, int_maji, t2_mkac, int_mbji,
c     &                     t2_micb, int_makj, t2_mica, int_mbkj,
c     &                     t2_mjcb, int_maki, t2_mjca, int_mbki,
c     &                     int_abij,  int_abjk, int_abik,  int_bcij,
c     &                     int_acij,  int_bcjk, int_acjk,  int_bcik,
c     &                     int_acik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
c     &                     nv_dim, no_dim, spin)
c
                        call uccsdt_get_t3x(c_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(k_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
C****      call writeblock('t2_ecji', t_buf, 
C****     &   e_blk_lo, e_blk_hi, c_blk_lo, c_blk_hi, 
C****     &   j_blk_lo, j_blk_hi, i_blk_lo, i_blk_hi)
C****      call writeblock('int_ekab', int_buf, 
C****     &   e_blk_lo, e_blk_hi, k_blk_lo, k_blk_hi, 
C****     &   a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi)
                        call w_t2_ecji_int_ekab
     &                    (w, t_buf, int_buf,
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(c_blk, j_blk, k_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(i_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_ecjk_int_eiab
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(c_blk, i_blk, k_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_ecik_int_ejab
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(k_blk, b_blk, c_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_eaji_int_ekbc
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(k_blk, a_blk, c_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_ebji_int_ekac
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(a_blk, k_blk, j_blk,
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(i_blk, c_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_eakj_int_eicb
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(b_blk, k_blk, j_blk,
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(i_blk, c_blk, a_blk,
     &                     spin, spin, spin, int_buf)
                        call w_t2_ebkj_int_eica
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(a_blk, k_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(j_blk, c_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_eaki_int_ejcb
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(b_blk, k_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(j_blk, c_blk, a_blk,
     &                     spin, spin, spin, int_buf)
                        call w_t2_ebki_int_ejca
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(k_blk, a_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(c_blk, j_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mkab_int_mcji
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(i_blk, a_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(c_blk, j_blk, k_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_miab_int_mcjk
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(c_blk, i_blk, k_blk,
     &                     spin, spin, spin, int_buf)
                        call w_t2_mjab_int_mcik
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(k_blk, b_blk, c_blk,
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mkbc_int_maji
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(k_blk, a_blk, c_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mkac_int_mbji
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(i_blk, c_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(a_blk, k_blk, j_blk,
     &                     spin, spin, spin, int_buf)
                        call w_t2_micb_int_makj
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(i_blk, c_blk, a_blk,
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(b_blk, k_blk, j_blk,
     &                     spin, spin, spin, int_buf)
                        call w_t2_mica_int_mbkj
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(j_blk, c_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(a_blk, k_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mjcb_int_maki
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(j_blk, c_blk, a_blk,
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(b_blk, k_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mjca_int_mbki
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
c                        call dcopy(100000, w, 1, v, 1)
c
c                        call w_t1a_int_abij
c     &                    (v, w, t1a, int_abij,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_abjk
c     &                    (v, w, t1a, int_abjk,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_abik
c     &                    (v, w, t1a, int_abik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_bcij
c     &                    (v, w, t1a, int_bcij,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_acij
c     &                    (v, w, t1a, int_acij,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_bcjk
c     &                    (v, w, t1a, int_bcjk,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_acjk
c     &                    (v, w, t1a, int_acjk,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_bcik
c     &                    (v, w, t1a, int_bcik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c                        call w_t1a_int_acik
c     &                    (v, w, t1a, int_acik,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi,
c     &                     nv_dim, no_dim, spin, spinck)
c
                        call uccsdt_get_2x(a_blk, b_blk, i_blk, j_blk,
     &                     spin, spin, spin, spin, int_abij)
                        call uccsdt_get_2x(a_blk, b_blk, j_blk, k_blk,
     &                     spin, spin, spin, spin, int_abjk)
                        call uccsdt_get_2x(a_blk, b_blk, i_blk, k_blk,
     &                     spin, spin, spin, spin, int_abik)
                        call uccsdt_get_2x(b_blk, c_blk, i_blk, j_blk,
     &                     spin, spin, spin, spin, int_bcij)
                        call uccsdt_get_2x(a_blk, c_blk, i_blk, j_blk,
     &                     spin, spin, spin, spin, int_acij)
                        call uccsdt_get_2x(b_blk, c_blk, j_blk, k_blk,
     &                     spin, spin, spin, spin, int_bcjk)
                        call uccsdt_get_2x(a_blk, c_blk, j_blk, k_blk,
     &                     spin, spin, spin, spin, int_acjk)
                        call uccsdt_get_2x(b_blk, c_blk, i_blk, k_blk, 
     &                     spin, spin, spin, spin, int_bcik)
                        call uccsdt_get_2x(a_blk, c_blk, i_blk, k_blk, 
     &                     spin, spin, spin, spin, int_acik)
c
                        call w_t1a_int_energy
     &                    (w, t1a, evals, energy,
     &                     int_abij,  int_abjk, int_abik,  int_bcij,
     &                     int_acij,  int_bcjk, int_acjk,  int_bcik,
     &                     int_acik,
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     nv_dim, no_dim, spin)
c
c                        call w_v_d_energy
c     &                    (v, w, evals, energy2,
c     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
c     &                     c_blk_lo, c_blk_hi,
c     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
c     &                     k_blk_lo, k_blk_hi)
                     end do
                  end do
               end do
            end do
         end do
      end do
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
                        v = w 
     $                       + t1b(c,k)*iaa(i,j,a+noa,b+noa)
     $                       + t1a(a,i)*iab(j,k,b+noa,c+nob)
     $                       - t1a(b,i)*iab(j,k,a+noa,c+nob)
     $                       - t1a(a,j)*iab(i,k,b+noa,c+nob)
     $                       + t1a(b,j)*iab(i,k,a+noa,c+nob)
c
                        d = evalsa(a+noa)+evalsa(b+noa)+evalsb(c+nob)
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
      double precision function uccsdtest_triples_mixed_blocked(
     $     noa, nob, nva, nvb, 
     $     iaa, ibb, iab,
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab,
     $     evalsa, evalsb, spin)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
      integer noa, nob, nva, nvb
      double precision iaa(nmo, nmo, nmo, nmo)
      double precision ibb(nmo, nmo, nmo, nmo)
      double precision iab(nmo, nmo, nmo, nmo)
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob) 
      double precision evalsa(nmo), evalsb(nmo)
      integer spin, spinck
c
      integer a_blk, a_blk_lo, a_blk_hi, a_blk_sym
      integer b_blk, b_blk_lo, b_blk_hi, b_blk_sym
      integer c_blk, c_blk_lo, c_blk_hi, c_blk_sym
      integer i_blk, i_blk_lo, i_blk_hi, i_blk_sym
      integer j_blk, j_blk_lo, j_blk_hi, j_blk_sym
      integer k_blk, k_blk_lo, k_blk_hi, k_blk_sym
      integer m_blk, m_blk_lo, m_blk_hi, m_blk_sym
      integer e_blk, e_blk_lo, e_blk_hi, e_blk_sym
      double precision energy, energy2
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
      double precision tmpr(20*20*20*20)
c
      call output(evalsa,1,nmo,1,1,nmo,1,1)
      call output(evalsb,1,nmo,1,1,nmo,1,1)
c
      energy = 0.0d0
c
c     note c&k have opposite spins of a&i
c
      if (spin.eq.1)then
         spinck=2
      elseif (spin.eq.2)then
         spinck=1
      else
         call errquit(' unrecognized spin ',spin, UNKNOWN_ERR)
      endif
c
      write(6,*)' spin, spinck, nvblock(spin), nvblock(spinck): ',
     &            spin, spinck, nvblock(spin), nvblock(spinck)
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
            do c_blk = 1, nvblock(spinck)
               c_blk_lo = vblock(1,c_blk,spinck)
               c_blk_hi = vblock(2,c_blk,spinck)
               c_blk_sym = vblock_sym(c_blk,spinck)
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
                     do k_blk = 1, noblock(spinck)
                        k_blk_lo = oblock(1,k_blk,spinck)
                        k_blk_hi = oblock(2,k_blk,spinck)
                        k_blk_sym = oblock_sym(k_blk,spinck)
c
      write(6,*)' k_blk, k_blk_lo, k_blk_hi, k_blk_sym ',
     &            k_blk, k_blk_lo, k_blk_hi, k_blk_sym
c
c                       need two loops over e and m (spin and spinck).
c
c                       e has spin
c
                        e_blk_lo = no(spin)+1
                        e_blk_hi = no(spin)+nv(spin)
c
c                       m has spin 
c
                        m_blk_lo = 1
                        m_blk_hi = no(spin)
c
                        call dfill (100000, 0.0d0, w, 1)
c
                        call uccsdt_get_t3x(c_blk, j_blk, k_blk, 
     &                     spinck, spin, spinck, t_buf)
                        call uccsdt_get_3x(i_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_ecjk_int_eiab
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(c_blk, i_blk, k_blk, 
     &                     spinck, spin, spinck, t_buf)
                        call uccsdt_get_3x(j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_ecik_int_ejab
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(k_blk, b_blk, c_blk, 
     &                     spinck, spin, spinck, int_buf)
                        call w_t2_eaji_int_ekbc
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_3x(k_blk, a_blk, c_blk, 
     &                     spinck, spin, spinck, int_buf)
                        call w_t2_ebji_int_ekac
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
c
c                       need two loops over e and m (spin and spinck).
c
c                       e has spinck
c
                        e_blk_lo = no(spinck)+1
                        e_blk_hi = no(spinck)+nv(spinck)
c
                        call uccsdt_get_t3x(a_blk, k_blk, j_blk,
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_3x(i_blk, c_blk, b_blk, 
     &                     spin, spinck, spin, int_buf)
                        call w_t2_eakj_int_eicb
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck, tmpr)
c
                        call uccsdt_get_t3x(b_blk, k_blk, j_blk,
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_3x(i_blk, c_blk, a_blk,
     &                     spin, spinck, spin, int_buf)
                        call w_t2_ebkj_int_eica
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(a_blk, k_blk, i_blk, 
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_3x(j_blk, c_blk, b_blk, 
     &                     spin, spinck, spin, int_buf)
                        call w_t2_eaki_int_ejcb
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t3x(b_blk, k_blk, i_blk, 
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_3x(j_blk, c_blk, a_blk,
     &                     spin, spinck, spin, int_buf)
                        call w_t2_ebki_int_ejca
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi, e_blk_lo, e_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     spin, spinck)
c
c                       m has spin 
c
                        m_blk_lo = 1
                        m_blk_hi = no(spin)
c
                        call uccsdt_get_t1x(i_blk, a_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(c_blk, j_blk, k_blk, 
     &                     spinck, spin, spinck, int_buf)
                        call w_t2_miab_int_mcjk
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(j_blk, a_blk, b_blk, 
     &                     spin, spin, spin, t_buf)
                        call uccsdt_get_1x(c_blk, i_blk, k_blk,
     &                     spinck, spin, spinck, int_buf)
                        call w_t2_mjab_int_mcik
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(k_blk, b_blk, c_blk,
     &                     spinck, spin, spinck, t_buf)
                        call uccsdt_get_1x(a_blk, j_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mkbc_int_maji
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(k_blk, a_blk, c_blk, 
     &                     spinck, spin, spinck, t_buf)
                        call uccsdt_get_1x(b_blk, j_blk, i_blk, 
     &                     spin, spin, spin, int_buf)
                        call w_t2_mkac_int_mbji
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
c                       m has spinck
c
                        m_blk_lo = 1
                        m_blk_hi = no(spinck)
c
                        call uccsdt_get_t1x(i_blk, c_blk, b_blk, 
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_1x(a_blk, k_blk, j_blk,
     &                     spin, spinck, spin, int_buf)
                        call w_t2_micb_int_makj
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call uccsdt_get_t1x(i_blk, c_blk, a_blk,
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_1x(b_blk, k_blk, j_blk,
     &                     spin, spinck, spin, int_buf)
                        call w_t2_mica_int_mbkj
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin)
c
                        call uccsdt_get_t1x(j_blk, c_blk, b_blk, 
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_1x(a_blk, k_blk, i_blk, 
     &                     spin, spinck, spin, int_buf)
                        call w_t2_mjcb_int_maki
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin)
c
                        call uccsdt_get_t1x(j_blk, c_blk, a_blk,
     &                     spin, spinck, spin, t_buf)
                        call uccsdt_get_1x(b_blk, k_blk, i_blk, 
     &                     spin, spinck, spin, int_buf)
                        call w_t2_mjca_int_mbki
     &                    (w, t_buf, int_buf, 
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi, m_blk_lo, m_blk_hi,
     &                     spin, spinck)
c
                        call dfill (5000, 0.0d0, int_abij, 1)
                        call dfill (5000, 0.0d0, int_bcjk, 1)
                        call dfill (5000, 0.0d0, int_acjk, 1)
                        call dfill (5000, 0.0d0, int_bcik, 1)
                        call dfill (5000, 0.0d0, int_acik, 1)
c
                        call uccsdt_get_2x(a_blk, b_blk, i_blk, j_blk,
     &                     spin, spin, spin, spin, int_abij)
c
                        call uccsdt_get_2x(b_blk, c_blk, j_blk, k_blk,
     &                     spin, spinck, spin, spinck, int_bcjk)
c
                        call uccsdt_get_2x(a_blk, c_blk, j_blk, k_blk,
     &                     spin, spinck, spin, spinck, int_acjk)
c
                        call uccsdt_get_2x(b_blk, c_blk, i_blk, k_blk, 
     &                     spin, spinck, spin, spinck, int_bcik)
c
                        call uccsdt_get_2x(a_blk, c_blk, i_blk, k_blk, 
     &                     spin, spinck, spin, spinck, int_acik)
c
                        call w_mixed_t1a_int_energy
     &                    (w, t1a, t1b, evalsa, evalsb, energy,
     &                     int_abij,  
     &                     int_bcjk, int_acjk,  int_bcik,
     &                     int_acik,
     &                     a_blk_lo, a_blk_hi, b_blk_lo, b_blk_hi, 
     &                     c_blk_lo, c_blk_hi,
     &                     i_blk_lo, i_blk_hi, j_blk_lo, j_blk_hi, 
     &                     k_blk_lo, k_blk_hi,
     &                     noa, nob, nva, nvb, spin)
c
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      uccsdtest_triples_mixed_blocked = energy
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
               vblock_sym(nvblock(spin),spin) = 0
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
c      write(6,*)' alo, ahi, blo, bhi, klo, khi, elo, ehi: ',
c     &            alo, ahi, blo, bhi, klo, khi, elo, ehi
c
      call uccsdt_copy_int_block
     $     (ints,dbl_mb(base),nmo,nmo,nmo,
     $     elo,ehi,klo,khi,alo,ahi,blo,bhi,
     $     oflip_ek,oflip_ab,.false.)
c
      end
      subroutine uccsdt_get_2x(ablock, bblock, iblock, jblock, 
     $     spina, spinb, spini, spinj, ints)
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
     $     alo,ahi,blo,bhi,ilo,ihi,jlo,jhi,
     $     oflip_ab,oflip_ij,.false.)
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
c                           write(6,*)' e, k, a, b, g: ',
c     &                                 e, k, a, b, g(e,k,a,b)
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
c      write(6,*)' elo, ehi, clo, chi, ilo, ihi, jlo, jhi: ',
c     &            elo, ehi, clo, chi, ilo, ihi, jlo, jhi
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
      alo = vblock(1,ablock,spina) - no(spina)
      ahi = vblock(2,ablock,spina) - no(spina)
      blo = vblock(1,bblock,spinb) - no(spinb)
      bhi = vblock(2,bblock,spinb) - no(spinb)
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
     &   t2_ecji, int_ekab, t2_ecjk, int_eiab,
     &   t2_ecik, int_ejab, t2_eaji, int_ekbc,
     &   t2_ebji, int_ekac, t2_eakj, int_eicb,
     &   t2_ebkj, int_eica, t2_eaki, int_ejcb,
     &   t2_ebki, int_ejca, t2_mkab, int_mcji,
     &   t2_miab, int_mcjk, t2_mjab, int_mcik,
     &   t2_mkbc, int_maji, t2_mkac, int_mbji,
     &   t2_micb, int_makj, t2_mica, int_mbkj,
     &   t2_mjcb, int_maki, t2_mjca, int_mbki,
     &   int_abij,  int_abjk, int_abik,  int_bcij,
     &   int_acij,  int_bcjk, int_acjk,  int_bcik,
     &   int_acik,
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
      double precision t2_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision int_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      double precision t2_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision int_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      double precision t2_ecik(elo:ehi,clo:chi,ilo:ihi,klo:khi)
      double precision int_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      double precision t2_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision int_ekbc(elo:ehi,klo:khi,blo:bhi,clo:chi)
      double precision t2_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision int_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      double precision t2_eakj(elo:ehi,alo:ahi,klo:khi,jlo:jhi)
      double precision int_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision t2_ebkj(elo:ehi,blo:bhi,klo:khi,jlo:jhi)
      double precision int_eica(elo:ehi,ilo:ihi,clo:chi,alo:ahi)
      double precision t2_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision int_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      double precision t2_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_ejca(elo:ehi,jlo:jhi,clo:chi,alo:ahi)
      double precision t2_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision int_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision t2_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision int_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      double precision t2_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision int_mcik(mlo:mhi,clo:chi,ilo:ihi,klo:khi)
      double precision t2_mkbc(mlo:mhi,klo:khi,blo:bhi,clo:chi)
      double precision int_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision t2_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision int_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision t2_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision int_makj(mlo:mhi,alo:ahi,klo:khi,jlo:jhi)
      double precision t2_mica(mlo:mhi,ilo:ihi,clo:chi,alo:ahi)
      double precision int_mbkj(mlo:mhi,blo:bhi,klo:khi,jlo:jhi)
      double precision t2_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision int_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      double precision t2_mjca(mlo:mhi,jlo:jhi,clo:chi,alo:ahi)
      double precision int_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
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
                           t2_ecji(e,c,j,i) = 
     &                     t2aa(c-no(spin),e-no(spin),i,j)
c
c                  if (abs(iaa(a,b,e,k)).gt.1d-10)
c     &               write(6,*)' a, b, e, k, iaa(a,b,e,k): ',
c     &                           a, b, e, k, iaa(a,b,e,k)
c
                           
                           int_ekab(e,k,a,b) = iaa(a,b,e,k)
c
c       if (abs(int_abek(a,b,e,k)).gt.1d-10)
c     &    write(6,*)' a, b, e, k, int_abek(a,b,e,k): ',
c     &                a, b, e, k, int_abek(a,b,e,k)
c
                           t2_ecjk(e,c,j,k) = 
     &                     t2aa(c-no(spin),e-no(spin),k,j)
c
                           int_eiab(e,i,a,b) = iaa(a,b,e,i)
c
                           t2_ecik(e,c,i,k) = 
     &                     t2aa(c-no(spin),e-no(spin),i,k)
c
                           int_ejab(e,j,a,b) = iaa(a,b,e,j)
c
                           t2_eaji(e,a,j,i) = 
     &                     t2aa(a-no(spin),e-no(spin),i,j)
c
                           int_ekbc(e,k,b,c) = iaa(b,c,e,k)
c
                           t2_ebji(e,b,j,i) = 
     &                     t2aa(b-no(spin),e-no(spin),i,j)
c
                           int_ekac(e,k,a,c) = iaa(a,c,e,k)
c
                           t2_eakj(e,a,k,j) = 
     &                     t2aa(a-no(spin),e-no(spin),j,k)
c
                           int_eicb(e,i,c,b) = iaa(c,b,e,i)
c
                           t2_ebkj(e,b,k,j) = 
     &                     t2aa(b-no(spin),e-no(spin),j,k)
c
                           int_eica(e,i,c,a) = iaa(c,a,e,i)
c
                           t2_eaki(e,a,k,i) = 
     &                     t2aa(a-no(spin),e-no(spin),i,k)
c
                           int_ejcb(e,j,c,b) = iaa(c,b,e,j)
c
                           t2_ebki(e,b,k,i) = 
     &                     t2aa(b-no(spin),e-no(spin),i,k)
c
                           int_ejca(e,j,c,a) = iaa(c,a,e,j)
c
                        end do
c
                        do m = 1, no(spin)
c
                           t2_mkab(m,k,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,k)
c
                           int_mcji(m,c,j,i) = iaa(c,m,i,j)
c
                           t2_miab(m,i,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,i)
c
                           int_mcjk(m,c,j,k) = iaa(c,m,k,j)
c
                           t2_mjab(m,j,a,b) = 
     &                     t2aa(a-no(spin),b-no(spin),m,j)
c
                           int_mcik(m,c,i,k) = iaa(m,c,i,k)
c
                           t2_mkbc(m,k,b,c) = 
     &                     t2aa(b-no(spin),c-no(spin),m,k)
c
                           int_maji(m,a,j,i) = iaa(a,m,i,j)
c
                           t2_mkac(m,k,a,c) = 
     &                     t2aa(a-no(spin),c-no(spin),m,k)
c
                           int_mbji(m,b,j,i) = iaa(b,m,i,j)
c
                           t2_micb(m,i,c,b) = 
     &                     t2aa(c-no(spin),b-no(spin),m,i)
c
                           int_makj(m,a,k,j) = iaa(a,m,j,k)
c
                           t2_mica(m,i,c,a) = 
     &                     t2aa(a-no(spin),c-no(spin),i,m)
c
                           int_mbkj(m,b,k,j) = iaa(b,m,j,k)
c
                           t2_mjcb(m,j,c,b) = 
     &                     t2aa(c-no(spin),b-no(spin),m,j)
c
                           int_maki(m,a,k,i) = iaa(a,m,i,k)
c
                           t2_mjca(m,j,c,a) = 
     &                     t2aa(a-no(spin),c-no(spin),j,m)
c
                           int_mbki(m,b,k,i) = iaa(b,m,i,k)
c
                        end do
c
c                       assume t1 in memory set up integral blocks only
c
                        int_abij(a,b,i,j) = iaa(a,b,i,j)
c
                        int_abjk(a,b,j,k) = iaa(a,b,j,k)
c
                        int_abik(a,b,i,k) = iaa(a,b,i,k)
c
                        int_bcij(b,c,i,j) = iaa(b,c,i,j)
c
                        int_acij(a,c,i,j) = iaa(a,c,i,j)
c
                        int_bcjk(b,c,j,k) = iaa(b,c,j,k)
c
                        int_acjk(a,c,j,k) = iaa(a,c,j,k)
c
                        int_bcik(b,c,i,k) = iaa(b,c,i,k)
c
                        int_acik(a,c,i,k) = iaa(a,c,i,k)
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
      return
      end
      subroutine triples_temp_dump_buffers(
     &   t2_ecji, int_ekab, t2_ecjk, int_eiab,
     &   t2_ecik, int_ejab, t2_eaji, int_ekbc,
     &   t2_ebji, int_ekac, t2_eakj, int_eicb,
     &   t2_ebkj, int_eica, t2_eaki, int_ejcb,
     &   t2_ebki, int_ejca, t2_mkab, int_mcji,
     &   t2_miab, int_mcjk, t2_mjab, int_mcik,
     &   t2_mkbc, int_maji, t2_mkac, int_mbji,
     &   t2_micb, int_makj, t2_mica, int_mbkj,
     &   t2_mjcb, int_maki, t2_mjca, int_mbki,
     &   int_abij,  int_abjk, int_abik,  int_bcij,
     &   int_acij,  int_bcjk, int_acjk,  int_bcik,
     &   int_acik,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi)
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
     &                 mlo, mhi
c
      double precision t2_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision int_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      double precision t2_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision int_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      double precision t2_ecik(elo:ehi,clo:chi,ilo:ihi,klo:khi)
      double precision int_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      double precision t2_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision int_ekbc(elo:ehi,klo:khi,blo:bhi,clo:chi)
      double precision t2_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision int_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      double precision t2_eakj(elo:ehi,alo:ahi,klo:khi,jlo:jhi)
      double precision int_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision t2_ebkj(elo:ehi,blo:bhi,klo:khi,jlo:jhi)
      double precision int_eica(elo:ehi,ilo:ihi,clo:chi,alo:ahi)
      double precision t2_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision int_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      double precision t2_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_ejca(elo:ehi,jlo:jhi,clo:chi,alo:ahi)
      double precision t2_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision int_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision t2_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision int_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      double precision t2_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision int_mcik(mlo:mhi,clo:chi,ilo:ihi,klo:khi)
      double precision t2_mkbc(mlo:mhi,klo:khi,blo:bhi,clo:chi)
      double precision int_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision t2_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision int_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision t2_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision int_makj(mlo:mhi,alo:ahi,klo:khi,jlo:jhi)
      double precision t2_mica(mlo:mhi,ilo:ihi,clo:chi,alo:ahi)
      double precision int_mbkj(mlo:mhi,blo:bhi,klo:khi,jlo:jhi)
      double precision t2_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision int_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      double precision t2_mjca(mlo:mhi,jlo:jhi,clo:chi,alo:ahi)
      double precision int_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
c
c temp write buffers
c
      call writeblock('t2_ecji', t2_ecji, 
     &                elo, ehi, clo, chi, jlo, jhi, ilo, ihi)
      call writeblock('int_ekab', int_ekab, 
     &                elo, ehi, klo, khi, alo, ahi, blo, bhi)
      call writeblock('t2_ecjk', t2_ecjk, 
     &                elo, ehi, clo, chi, jlo, jhi, klo, khi)
      call writeblock('int_eiab', int_eiab, 
     &                elo, ehi, ilo, ihi, alo, ahi, blo, bhi)
      call writeblock('t2_ecik', t2_ecik, 
     &                elo, ehi, clo, chi, ilo, ihi, klo, khi)
      call writeblock('int_ejab', int_ejab, 
     &                elo, ehi, jlo, jhi, alo, ahi, blo, bhi)
      call writeblock('t2_eaji', t2_eaji, 
     &                elo, ehi, alo, ahi, jlo, jhi, ilo, ihi)
      call writeblock('int_ekbc', int_ekbc, 
     &                elo, ehi, klo, khi, blo, bhi, clo, chi)
      call writeblock('t2_ebji', t2_ebji, 
     &                elo, ehi, blo, bhi, jlo, jhi, ilo, ihi)
      call writeblock('int_ekac', int_ekac, 
     &                elo, ehi, klo, khi, alo, ahi, clo, chi)
      call writeblock('t2_eakj', t2_eakj, 
     &                elo, ehi, alo, ahi,  klo, khi, jlo, jhi)
      call writeblock('int_eicb', int_eicb, 
     &                elo, ehi, ilo, ihi, clo, chi, blo, bhi)
      call writeblock('t2_ebkj', t2_ebkj, 
     &                elo, ehi, blo, bhi, klo, khi, jlo, jhi)
      call writeblock('int_eica', int_eica, 
     &                elo, ehi, ilo, ihi, clo, chi, alo, ahi)
      call writeblock('t2_eaki', t2_eaki, 
     &                elo, ehi, alo, ahi, klo, khi, ilo, ihi)
      call writeblock('int_ejcb', int_ejcb, 
     &                elo, ehi, jlo, jhi, clo, chi, blo, bhi)
      call writeblock('t2_ebki', t2_ebki, 
     &                elo, ehi, blo, bhi, klo, khi, ilo, ihi)
      call writeblock('int_ejca', int_ejca, 
     &                elo, ehi, jlo, jhi, clo, chi, alo, ahi)
      call writeblock('t2_mkab', t2_mkab, 
     &                mlo, mhi, klo, khi, alo, ahi, blo, bhi)
      call writeblock('int_mcji', int_mcji, 
     &                mlo, mhi, clo, chi, jlo, jhi, ilo, ihi)
      call writeblock('t2_miab', t2_miab, 
     &                mlo, mhi, ilo, ihi, alo, ahi, blo, bhi)
      call writeblock('int_mcjk', int_mcjk, 
     &                mlo, mhi, clo, chi, jlo, jhi, klo, khi)
      call writeblock('t2_mjab', t2_mjab, 
     &                mlo, mhi, jlo, jhi, alo, ahi, blo, bhi)
      call writeblock('int_mcik', int_mcik, 
     &                mlo, mhi, clo, chi, ilo, ihi, klo, khi)
      call writeblock('t2_mkbc', t2_mkbc, 
     &                mlo, mhi, klo, khi, blo, bhi, clo, chi)
      call writeblock('int_maji', int_maji, 
     &                mlo, mhi, alo, ahi, jlo, jhi, ilo, ihi)
      call writeblock('t2_mkac', t2_mkac, 
     &                mlo, mhi, klo, khi, alo, ahi, clo, chi)
      call writeblock('int_mbji', int_mbji, 
     &                mlo, mhi, blo, bhi, jlo, jhi, ilo, ihi)
      call writeblock('t2_micb', t2_micb, 
     &                mlo, mhi, ilo, ihi, clo, chi, blo, bhi)
      call writeblock('int_makj', int_makj, 
     &                mlo, mhi, alo, ahi, klo, khi, jlo, jhi)
      call writeblock('t2_mica', t2_mica, 
     &                mlo, mhi, ilo, ihi, clo, chi, alo, ahi)
      call writeblock('int_mbkj', int_mbkj, 
     &                mlo, mhi, blo, bhi, klo, khi, jlo, jhi)
      call writeblock('t2_mjcb', t2_mjcb, 
     &                mlo, mhi, jlo, jhi, clo, chi, blo, bhi)
      call writeblock('int_maki', int_maki, 
     &                mlo, mhi, alo, ahi, klo, khi, ilo, ihi)
      call writeblock('t2_mjca', t2_mjca, 
     &                mlo, mhi, jlo, jhi, clo, chi, alo, ahi)
      call writeblock('int_mbki', int_mbki, 
     &                mlo, mhi, blo, bhi, klo, khi, ilo, ihi)
      call writeblock('int_abij', int_abij, 
     &                alo, ahi, blo, bhi, ilo, ihi, jlo, jhi)
      call writeblock('int_abjk', int_abjk, 
     &                alo, ahi, blo, bhi, jlo, jhi, klo, khi)
      call writeblock('int_abik', int_abik, 
     &                alo, ahi, blo, bhi, ilo, ihi, klo, khi)
      call writeblock('int_bcij', int_bcij, 
     &                blo, bhi, clo, chi, ilo, ihi, jlo, jhi)
      call writeblock('int_acij', int_acij, 
     &                alo, ahi, clo, chi, ilo, ihi, jlo, jhi)
      call writeblock('int_bcjk', int_bcjk, 
     &                blo, bhi, clo, chi, jlo, jhi, klo, khi)
      call writeblock('int_acjk', int_acjk, 
     &                alo, ahi, clo, chi, jlo, jhi, klo, khi)
      call writeblock('int_bcik', int_bcik, 
     &                blo, bhi, clo, chi, ilo, ihi, klo, khi)
      call writeblock('int_acik', int_acik, 
     &                alo, ahi, clo, chi, ilo, ihi, klo, khi)
c
      return
      end
      subroutine uccsdtest_triples_pure_blocked_energy_45
     &  (t1a, evals, energy,
     &   t2_ecji, int_ekab, t2_ecjk, int_eiab,
     &   t2_ecik, int_ejab, t2_eaji, int_ekbc,
     &   t2_ebji, int_ekac, t2_eakj, int_eicb,
     &   t2_ebkj, int_eica, t2_eaki, int_ejcb,
     &   t2_ebki, int_ejca, t2_mkab, int_mcji,
     &   t2_miab, int_mcjk, t2_mjab, int_mcik,
     &   t2_mkbc, int_maji, t2_mkac, int_mbji,
     &   t2_micb, int_makj, t2_mica, int_mbkj,
     &   t2_mjcb, int_maki, t2_mjca, int_mbki,
     &   int_abij,  int_abjk, int_abik,  int_bcij,
     &   int_acij,  int_bcjk, int_acjk,  int_bcik,
     &   int_acik,
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
      double precision t1a(nv_dim, no_dim)
      double precision evals(nmo)
c
      double precision t2_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision int_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      double precision t2_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision int_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      double precision t2_ecik(elo:ehi,clo:chi,ilo:ihi,klo:khi)
      double precision int_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      double precision t2_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision int_ekbc(elo:ehi,klo:khi,blo:bhi,clo:chi)
      double precision t2_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision int_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      double precision t2_eakj(elo:ehi,alo:ahi,klo:khi,jlo:jhi)
      double precision int_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision t2_ebkj(elo:ehi,blo:bhi,klo:khi,jlo:jhi)
      double precision int_eica(elo:ehi,ilo:ihi,clo:chi,alo:ahi)
      double precision t2_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision int_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      double precision t2_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_ejca(elo:ehi,jlo:jhi,clo:chi,alo:ahi)
      double precision t2_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision int_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision t2_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision int_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      double precision t2_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision int_mcik(mlo:mhi,clo:chi,ilo:ihi,klo:khi)
      double precision t2_mkbc(mlo:mhi,klo:khi,blo:bhi,clo:chi)
      double precision int_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision t2_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision int_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision t2_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision int_makj(mlo:mhi,alo:ahi,klo:khi,jlo:jhi)
      double precision t2_mica(mlo:mhi,ilo:ihi,clo:chi,alo:ahi)
      double precision int_mbkj(mlo:mhi,blo:bhi,klo:khi,jlo:jhi)
      double precision t2_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision int_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      double precision t2_mjca(mlo:mhi,jlo:jhi,clo:chi,alo:ahi)
      double precision int_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
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
                        do e = elo, ehi
c                        do e = no(spin)+1, no(spin)+nv(spin)
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w = w 
     &                           + t2_ecji(e,c,j,i)*
     &                             int_ekab(e,k,a,b)
     &                           - t2_ecjk(e,c,j,k)*
     &                             int_eiab(e,i,a,b)
     &                           + t2_ecik(e,c,i,k)*
     &                             int_ejab(e,j,a,b)
     &                           + t2_eaji(e,a,j,i)*
     &                             int_ekbc(e,k,b,c)
     &                           - t2_ebji(e,b,j,i)*
     &                             int_ekac(e,k,a,c)
     &                           - t2_eakj(e,a,k,j)*
     &                             int_eicb(e,i,c,b)
     &                           + t2_ebkj(e,b,k,j)*
     &                             int_eica(e,i,c,a)
     &                           + t2_eaki(e,a,k,i)*
     &                             int_ejcb(e,j,c,b)
     &                           - t2_ebki(e,b,k,i)*
     &                             int_ejca(e,j,c,a)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e): ',
c     &                              a,b,c,i,j,k,w
c
c  later not 1 but number of orbitals in symmetry block                       
c
                        do m = mlo, mhi
c                        do m = 1, no(spin)
c
                           w = w 
     &                           - t2_mkab(m,k,a,b)*
     &                             int_mcji(m,c,j,i)
     &                           + t2_miab(m,i,a,b)*
     &                             int_mcjk(m,c,j,k)
     &                           - t2_mjab(m,j,a,b)*
     &                             int_mcik(m,c,i,k)
     &                           - t2_mkbc(m,k,b,c)*
     &                             int_maji(m,a,j,i)
     &                           + t2_mkac(m,k,a,c)*
     &                             int_mbji(m,b,j,i)
     &                           + t2_micb(m,i,c,b)*
     &                             int_makj(m,a,k,j)
     &                           - t2_mica(m,i,c,a)*
     &                             int_mbkj(m,b,k,j)
     &                           - t2_mjcb(m,j,c,b)*
     &                             int_maki(m,a,k,i)
     &                           + t2_mjca(m,j,c,a)*
     &                             int_mbki(m,b,k,i)
                        end do
c                        write(6,*)' a,b,c,i,j,k,w(e+m): ',
c     &                              a,b,c,i,j,k,w
                        v = w 
     &                        + t1a(c-no(spin),k)*
     &                          int_abij(a,b,i,j)
     &                        + t1a(c-no(spin),i)*
     &                          int_abjk(a,b,j,k)
     &                        - t1a(c-no(spin),j)*
     &                          int_abik(a,b,i,k)
     &                        + t1a(a-no(spin),k)*
     &                          int_bcij(b,c,i,j)
     &                        - t1a(b-no(spin),k)*
     &                          int_acij(a,c,i,j)
     &                        + t1a(a-no(spin),i)*
     &                          int_bcjk(b,c,j,k)
     &                        - t1a(b-no(spin),i)*
     &                          int_acjk(a,c,j,k)
     &                        - t1a(a-no(spin),j)*
     &                          int_bcik(b,c,i,k)
     &                        + t1a(b-no(spin),j)*
     &                          int_acik(a,c,i,k)
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
      write(6,*) ' E(45) block', energy
c
      return
      end
      subroutine w_t2_ecji_int_ekab
     &  (w, t2_ecji, int_ekab, 
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ecji(elo:ehi,clo:chi,jlo:jhi,ilo:ihi)
      double precision int_ekab(elo:ehi,klo:khi,alo:ahi,blo:bhi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_ecji(e,c,j,i)*
     &                             int_ekab(e,k,a,b)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(1, flopcount)
c
      end
      subroutine w_t2_ecjk_int_eiab
     &  (w, t2_ecjk, int_eiab,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ecjk(elo:ehi,clo:chi,jlo:jhi,klo:khi)
      double precision int_eiab(elo:ehi,ilo:ihi,alo:ahi,blo:bhi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_ecjk(e,c,j,k)*
     &                             int_eiab(e,i,a,b)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(2,flopcount)
c
      end
      subroutine w_t2_ecik_int_ejab
     &  (w, t2_ecik, int_ejab,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ecik(elo:ehi,clo:chi,ilo:ihi,klo:khi)
      double precision int_ejab(elo:ehi,jlo:jhi,alo:ahi,blo:bhi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_ecik(e,c,i,k)*
     &                             int_ejab(e,j,a,b)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(3, flopcount)
c
      end
      subroutine w_t2_eaji_int_ekbc
     &  (w, t2_eaji, int_ekbc,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_eaji(elo:ehi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision int_ekbc(elo:ehi,klo:khi,blo:bhi,clo:chi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_eaji(e,a,j,i)*
     &                             int_ekbc(e,k,b,c)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(4, flopcount)
c
      end
      subroutine w_t2_ebji_int_ekac
     &  (w, t2_ebji, int_ekac,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ebji(elo:ehi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision int_ekac(elo:ehi,klo:khi,alo:ahi,clo:chi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_ebji(e,b,j,i)*
     &                             int_ekac(e,k,a,c)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(5, flopcount)
c      
      end
      subroutine w_t2_eakj_int_eicb
     &  (w, t2_eakj, int_eicb,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck, r)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_eakj(elo:ehi,alo:ahi,klo:khi,jlo:jhi)
      double precision int_eicb(elo:ehi,ilo:ihi,clo:chi,blo:bhi)
      double precision r(ilo:ihi,alo:ahi,clo:chi,klo:khi)
      integer i, j, k, a, b, c, e, cend, kend, abegin, ibegin, cend4
      integer ihi4
      double precision sum, sum0, sum1, sum2, sum3
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
c     Near original code with re-ordered loops
c
c     jbkcaie = 124 Mflop
c     bacjkie = slower
c
      do j = jlo, jhi
         kend = khi
         if (spin.eq.spinck) kend = min(khi,j-1)
         ibegin = max(j+1,ilo)
         do b = blo, bhi
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            abegin = max(b+1,alo)
            do k = klo, kend
               do c = clo, cend
                  do a = abegin, ahi
                     do i = ibegin, ihi
                        sum = 0d0
                        do e = elo, ehi
                           sum = sum + t2_eakj(e,a,k,j)* 
     $                          int_eicb(e,i,c,b)
                        end do
                        w(k,j,i,c,b,a) = w(k,j,i,c,b,a) - sum
                     end do
                  end do
               end do
            end do
            flopcount = flopcount + 2d0*dble(kend-klo+1)*
     $           dble(ahi-abegin+1)*dble(cend-clo+1)*
     $           dble(ihi-ibegin+1)*dble(ehi-elo+1)
         end do
      end do
c
c     Near original code with re-ordered loops and buffer for sum
c
c$$$      do j = jlo, jhi
c$$$         kend = khi
c$$$         if (spin.eq.spinck) kend = min(khi,j-1)
c$$$         ibegin = max(j+1,ilo)
c$$$         do b = blo, bhi
c$$$            cend = chi
c$$$            if (spin.eq.spinck)cend = min(chi,b-1)
c$$$            abegin = max(b+1,alo)
c$$$            do k = klo, kend
c$$$               do c = clo, cend
c$$$                  do a = abegin, ahi
c$$$                     do i = ibegin, ihi
c$$$                        sum = 0d0
c$$$                        do e = elo, ehi
c$$$                           sum = sum + t2_eakj(e,a,k,j)* 
c$$$     $                          int_eicb(e,i,c,b)
c$$$                        end do
c$$$                        r(i,a,c,k) = sum
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$            do a = abegin, ahi
c$$$               do c = clo, cend
c$$$                  do k = klo, kend
c$$$                     do i = ibegin, ihi
c$$$                        w(k,j,i,c,b,a) = w(k,j,i,c,b,a) - r(i,a,c,k)
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$            flopcount = flopcount + 2d0*dble(kend-klo+1)*
c$$$     $           dble(ahi-abegin+1)*dble(cend-clo+1)*
c$$$     $           dble(ihi-ibegin+1)*dble(ehi-elo+1)
c$$$         end do
c$$$      end do
c
c     4-way unrolled code using i and buffer for sum
c
c$$$      do j = jlo, jhi
c$$$         kend = khi
c$$$         if (spin.eq.spinck) kend = min(khi,j-1)
c$$$         ibegin = max(j+1,ilo)
c$$$         ihi4 = ibegin + ((ihi-ibegin+1)/4)*4 - 1
c$$$         do b = blo, bhi
c$$$            cend = chi
c$$$            if (spin.eq.spinck)cend = min(chi,b-1)
c$$$            abegin = max(b+1,alo)
c$$$            do k = klo, kend
c$$$               do c = clo, cend
c$$$                  do a = abegin, ahi
c$$$                     do i = ibegin, ihi4, 4
c$$$                        sum0 = 0d0
c$$$                        sum1 = 0d0
c$$$                        sum2 = 0d0
c$$$                        sum3 = 0d0
c$$$                        do e = elo, ehi
c$$$                           sum0 = sum0 + t2_eakj(e,a,k,j)*
c$$$     $                          int_eicb(e,i  ,c,b)
c$$$                           sum1 = sum1 + t2_eakj(e,a,k,j)*
c$$$     $                          int_eicb(e,i+1,c,b)
c$$$                           sum2 = sum2 + t2_eakj(e,a,k,j)*
c$$$     $                          int_eicb(e,i+2,c,b)
c$$$                           sum3 = sum3 + t2_eakj(e,a,k,j)*
c$$$     $                          int_eicb(e,i+3,c,b)
c$$$                        end do
c$$$                        r(i  ,a,c,k) = sum0
c$$$                        r(i+1,a,c,k) = sum1
c$$$                        r(i+2,a,c,k) = sum2
c$$$                        r(i+3,a,c,k) = sum3
c$$$c$$$                        w(k,j,i  ,c,b,a) = w(k,j,i  ,c,b,a) - sum0
c$$$c$$$                        w(k,j,i+1,c,b,a) = w(k,j,i+1,c,b,a) - sum1
c$$$c$$$                        w(k,j,i+2,c,b,a) = w(k,j,i+2,c,b,a) - sum2
c$$$c$$$                        w(k,j,i+3,c,b,a) = w(k,j,i+3,c,b,a) - sum3
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$            if (ihi4 .lt. ihi) then
c$$$               do i = ihi4+1, ihi
c$$$                  do k = klo, kend
c$$$                     do c = clo, cend
c$$$                        do a = abegin, ahi
c$$$                           sum = 0d0
c$$$                           do e = elo, ehi
c$$$                              sum = sum + t2_eakj(e,a,k,j)* 
c$$$     $                             int_eicb(e,i,c,b)
c$$$                           end do
c$$$                           r(i,a,c,k) = sum
c$$$c$$$                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a) - sum
c$$$                        end do
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end if
c$$$            do a = abegin, ahi
c$$$               do c = clo, cend
c$$$                  do k = klo, kend
c$$$                     do i = ibegin, ihi
c$$$                        w(k,j,i,c,b,a) = w(k,j,i,c,b,a) - r(i,a,c,k)
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$            flopcount = flopcount + 2d0*dble(kend-klo+1)*
c$$$     $           dble(ahi-abegin+1)*dble(cend-clo+1)*
c$$$     $           dble(ihi-ibegin+1)*dble(ehi-elo+1)
c$$$         end do
c$$$      end do
c
c     4-way unrolled code using C
c
c$$$      do j = jlo, jhi
c$$$         kend = khi
c$$$         if (spin.eq.spinck) kend = min(khi,j-1)
c$$$         ibegin = max(j+1,ilo)
c$$$         do b = blo, bhi
c$$$            cend = chi
c$$$            if (spin.eq.spinck)cend = min(chi,b-1)
c$$$            abegin = max(b+1,alo)
c$$$            cend4 = clo + ((cend-clo+1)/4)*4 - 1
c$$$            do k = klo, kend
c$$$               do c = clo, cend4, 4
c$$$                  do a = abegin, ahi
c$$$                     do i = ibegin, ihi
c$$$                        sum0 = 0d0
c$$$                        sum1 = 0d0
c$$$                        sum2 = 0d0
c$$$                        sum3 = 0d0
c$$$                        do e = elo, ehi
c$$$                           sum0 = sum0 + t2_eakj(e,a,k,j)*
c$$$     $                          int_eicb(e,i,c  ,b)
c$$$                           sum1 = sum1 + t2_eakj(e,a,k,j)* 
c$$$     $                          int_eicb(e,i,c+1,b)
c$$$                           sum2 = sum2 + t2_eakj(e,a,k,j)* 
c$$$     $                          int_eicb(e,i,c+2,b)
c$$$                           sum3 = sum3 + t2_eakj(e,a,k,j)* 
c$$$     $                          int_eicb(e,i,c+3,b)
c$$$                        end do
c$$$                        w(k,j,i,c  ,b,a) = w(k,j,i,c  ,b,a) - sum0
c$$$                        w(k,j,i,c+1,b,a) = w(k,j,i,c+1,b,a) - sum1
c$$$                        w(k,j,i,c+2,b,a) = w(k,j,i,c+2,b,a) - sum2
c$$$                        w(k,j,i,c+3,b,a) = w(k,j,i,c+3,b,a) - sum3
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$            if (cend .gt. cend4) then
c$$$               do k = klo, kend
c$$$                  do c = cend4+1, cend
c$$$                     do a = abegin, ahi
c$$$                        do i = ibegin, ihi
c$$$                           sum = 0d0
c$$$                           do e = elo, ehi
c$$$                              sum = sum + t2_eakj(e,a,k,j)* 
c$$$     $                             int_eicb(e,i,c,b)
c$$$                           end do
c$$$                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a) - sum
c$$$                        end do
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end if
c$$$            flopcount = flopcount + 2d0*dble(kend-klo+1)*
c$$$     $           dble(ahi-abegin+1)*dble(cend-clo+1)*
c$$$     $           dble(ihi-ibegin+1)*dble(ehi-elo+1)
c$$$         end do
c$$$      end do
c
c     Original code
c
c$$$      do a = alo, ahi
c$$$         do b = blo, min(bhi,a-1)
c$$$            cend = chi
c$$$            if (spin.eq.spinck)cend = min(chi,b-1)
c$$$            do c = clo, cend
c$$$               do i = ilo, ihi
c$$$                  do j = jlo, min(jhi,i-1)
c$$$                     kend = khi
c$$$                     if (spin.eq.spinck)kend = min(khi,j-1)
c$$$                     do k = klo, kend
c$$$                        do e = elo, ehi
c$$$c
c$$$c                        write(6,*)' a,b,c,i,j,k,e: ',
c$$$c     &                              a,b,c,i,j,k,e
c$$$c
c$$$                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
c$$$     &                           - t2_eakj(e,a,k,j)*
c$$$     &                             int_eicb(e,i,c,b)
c$$$                        end do
c$$$                     end do
c$$$                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c
      call uccsdt_flops_exit(6, flopcount)
c
      end
      subroutine w_t2_ebkj_int_eica
     &  (w, t2_ebkj, int_eica,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ebkj(elo:ehi,blo:bhi,klo:khi,jlo:jhi)
      double precision int_eica(elo:ehi,ilo:ihi,clo:chi,alo:ahi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_ebkj(e,b,k,j)*
     &                             int_eica(e,i,c,a)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(7, flopcount)
c
      end
      subroutine w_t2_eaki_int_ejcb
     &  (w, t2_eaki, int_ejcb,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_eaki(elo:ehi,alo:ahi,klo:khi,ilo:ihi)
      double precision int_ejcb(elo:ehi,jlo:jhi,clo:chi,blo:bhi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_eaki(e,a,k,i)*
     &                             int_ejcb(e,j,c,b)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(8, flopcount)
c
      end
      subroutine w_t2_ebki_int_ejca
     &  (w, t2_ebki, int_ejca,
     &   alo, ahi, blo, bhi, clo, chi, elo, ehi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        elo, ehi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_ebki(elo:ehi,blo:bhi,klo:khi,ilo:ihi)
      double precision int_ejca(elo:ehi,jlo:jhi,clo:chi,alo:ahi)
      integer i, j, k, a, b, c, e, cend, kend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do e = elo, ehi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_ebki(e,b,k,i)*
     &                             int_ejca(e,j,c,a)
                        end do
                     end do
                     flopcount = flopcount + 2*(ehi-elo+1)*(kend-klo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(9, flopcount)
c
      end
      subroutine w_t2_mkab_int_mcji
     &  (w, t2_mkab, int_mcji,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mkab(mlo:mhi,klo:khi,alo:ahi,blo:bhi)
      double precision int_mcji(mlo:mhi,clo:chi,jlo:jhi,ilo:ihi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_mkab(m,k,a,b)*
     &                             int_mcji(m,c,j,i)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(10,flopcnt)
c
      end
      subroutine w_t2_miab_int_mcjk
     &  (w, t2_miab, int_mcjk,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_miab(mlo:mhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision int_mcjk(mlo:mhi,clo:chi,jlo:jhi,klo:khi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_miab(m,i,a,b)*
     &                             int_mcjk(m,c,j,k)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(11,flopcnt)
c
      end
      subroutine w_t2_mjab_int_mcik
     &  (w, t2_mjab, int_mcik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mjab(mlo:mhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision int_mcik(mlo:mhi,clo:chi,ilo:ihi,klo:khi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_mjab(m,j,a,b)*
     &                             int_mcik(m,c,i,k)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(12,flopcnt)
c
      end
      subroutine w_t2_mkbc_int_maji
     &  (w, t2_mkbc, int_maji,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mkbc(mlo:mhi,klo:khi,blo:bhi,clo:chi)
      double precision int_maji(mlo:mhi,alo:ahi,jlo:jhi,ilo:ihi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_mkbc(m,k,b,c)*
     &                             int_maji(m,a,j,i)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(13,flopcnt)
c
      end
      subroutine w_t2_mkac_int_mbji
     &  (w, t2_mkac, int_mbji,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mkac(mlo:mhi,klo:khi,alo:ahi,clo:chi)
      double precision int_mbji(mlo:mhi,blo:bhi,jlo:jhi,ilo:ihi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_mkac(m,k,a,c)*
     &                             int_mbji(m,b,j,i)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(14,flopcnt)
c
      end
      subroutine w_t2_micb_int_makj
     &  (w, t2_micb, int_makj,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_micb(mlo:mhi,ilo:ihi,clo:chi,blo:bhi)
      double precision int_makj(mlo:mhi,alo:ahi,klo:khi,jlo:jhi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_micb(m,i,c,b)*
     &                             int_makj(m,a,k,j)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(15,flopcnt)
c
      end
      subroutine w_t2_mica_int_mbkj
     &  (w, t2_mica, int_mbkj,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mica(mlo:mhi,ilo:ihi,clo:chi,alo:ahi)
      double precision int_mbkj(mlo:mhi,blo:bhi,klo:khi,jlo:jhi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_mica(m,i,c,a)*
     &                             int_mbkj(m,b,k,j)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(16,flopcnt)
c
      end
      subroutine w_t2_mjcb_int_maki
     &  (w, t2_mjcb, int_maki,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mjcb(mlo:mhi,jlo:jhi,clo:chi,blo:bhi)
      double precision int_maki(mlo:mhi,alo:ahi,klo:khi,ilo:ihi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           - t2_mjcb(m,j,c,b)*
     &                             int_maki(m,a,k,i)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(17,flopcnt)
c
      end
      subroutine w_t2_mjca_int_mbki
     &  (w, t2_mjca, int_mbki,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi, mlo, mhi,
     &   spin, spinck)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        mlo, mhi, spin, spinck
c
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision t2_mjca(mlo:mhi,jlo:jhi,clo:chi,alo:ahi)
      double precision int_mbki(mlo:mhi,blo:bhi,klo:khi,ilo:ihi)
      integer i, j, k, m, a, b, c, cend, kend
      double precision flopcnt
c
      call uccsdt_flops_enter()
      flopcnt = 0d0
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
                        do m = mlo, mhi
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           w(k,j,i,c,b,a) = w(k,j,i,c,b,a)
     &                           + t2_mjca(m,j,c,a)*
     &                             int_mbki(m,b,k,i)
                        end do
                     end do
                     flopcnt = flopcnt + (kend-klo+1)*(mhi-mlo+1)
                  end do
               end do
            end do
         end do
      end do
c
      call uccsdt_flops_exit(18,flopcnt)
c
      end
      subroutine w_t1a_int_abij
     &  (v, t1a, int_abij,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           + t1a(c-no(spin),k)*
     &                             int_abij(a,b,i,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_abjk
     &  (v, t1a, int_abjk,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           + t1a(c-no(spin),i)*
     &                             int_abjk(a,b,j,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_abik
     &  (v, t1a, int_abik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           - t1a(c-no(spin),j)*
     &                             int_abik(a,b,i,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_bcij
     &  (v, t1a, int_bcij,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           + t1a(a-no(spin),k)*
     &                             int_bcij(b,c,i,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_acij
     &  (v, t1a, int_acij,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           - t1a(b-no(spin),k)*
     &                             int_acij(a,c,i,j)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_bcjk
     &  (v, t1a, int_bcjk,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           + t1a(a-no(spin),i)*
     &                             int_bcjk(b,c,j,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_acjk
     &  (v, t1a, int_acjk,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           - t1a(b-no(spin),i)*
     &                             int_acjk(a,c,j,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_bcik
     &  (v, t1a, int_bcik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           - t1a(a-no(spin),j)*
     &                             int_bcik(b,c,i,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_t1a_int_acik
     &  (v, t1a, int_acik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin, spinck)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin, spinck
c
      double precision t1a(nv_dim, no_dim)
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
      integer i, j, k, a, b, c, cend, kend
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            cend = chi
            if (spin.eq.spinck)cend = min(chi,b-1)
            do c = clo, cend
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     kend = khi
                     if (spin.eq.spinck)kend = min(khi,j-1)
                     do k = klo, kend
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                           v(k,j,i,c,b,a) = v(k,j,i,c,b,a)
     &                           + t1a(b-no(spin),j)*
     &                             int_acik(a,c,i,k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
      subroutine w_mixed_t1a_int_energy
     &  (w, t1a, t1b, evalsa, evalsb, energy,
     &   int_abij, 
     &   int_bcjk, int_acjk, int_bcik,
     &   int_acik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   noa, nob, nva, nvb, spin)
c
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nva, nvb, noa, nob, spin
c
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision evalsa(nmo), evalsb(nmo)
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
      double precision energy, v, d
      integer i, j, k, a, b, c, spinck
c
c     note c&k have opposite spins of a&i
c
      if (spin.eq.1)then
         spinck=2
      elseif (spin.eq.2)then
         spinck=1
      else
         call errquit(' unrecognized spin ',spin, UNKNOWN_ERR)
      endif
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            do c = clo, chi
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     do k = klo, khi
c
c                        write(6,*)' a,b,c,i,j,k: ',
c     &                              a,b,c,i,j,k
c
c                        v = w 
c     $                       + t1b(c,k)*iaa(i,j,a+noa,b+noa)
c     $                       + t1a(a,i)*iab(j,k,b+noa,c+nob)
c     $                       - t1a(b,i)*iab(j,k,a+noa,c+nob)
c     $                       - t1a(a,j)*iab(i,k,b+noa,c+nob)
c     $                       + t1a(b,j)*iab(i,k,a+noa,c+nob)
c
c                        d = evalsa(a+noa)+evalsa(b+noa)+evalsb(c+noa)
c     $                       -evalsa(i)-evalsa(j)-evalsb(k)
c
c                        energy = energy - v*w/d
cc
                        v = w(k,j,i,c,b,a)
     &                        + t1b(c-no(spinck),k)*
     &                          int_abij(a,b,i,j)
     &                        + t1a(a-no(spin),i)*
     &                          int_bcjk(b,c,j,k)
     &                        - t1a(b-no(spin),i)*
     &                          int_acjk(a,c,j,k)
     &                        - t1a(a-no(spin),j)*
     &                          int_bcik(b,c,i,k)
     &                        + t1a(b-no(spin),j)*
     &                          int_acik(a,c,i,k)
c
c   needs to be modified for generic case
c
                        d = evalsa(a)+
     &                      evalsa(b)+
     &                      evalsb(c)-
     &                      evalsa(i)-
     &                      evalsa(j)-
     &                      evalsb(k)
c
                        energy = energy - 
     &                           v*w(k,j,i,c,b,a)/d
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
      write(6,*) ' E(2) block', energy
c
      return
      end
      subroutine w_t1a_int_energy
     &  (w, t1a, evals, energy,
     &   int_abij, int_abjk, int_abik, int_bcij,
     &   int_acij, int_bcjk, int_acjk, int_bcik,
     &   int_acik,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   nv_dim, no_dim, spin)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        nv_dim, no_dim, spin
c
      double precision t1a(nv_dim, no_dim)
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
      double precision evals(nmo), energy, v, d
      integer i, j, k, a, b, c
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            do c = clo, min(chi,b-1)
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     do k = klo, min(khi,j-1)
c
c                        write(6,*)' a,b,c,i,j,k: ',
c     &                              a,b,c,i,j,k
c
                        v = w(k,j,i,c,b,a)
     &                        + t1a(c-no(spin),k)*
     &                          int_abij(a,b,i,j)
     &                        + t1a(c-no(spin),i)*
     &                          int_abjk(a,b,j,k)
     &                        - t1a(c-no(spin),j)*
     &                          int_abik(a,b,i,k)
     &                        + t1a(a-no(spin),k)*
     &                          int_bcij(b,c,i,j)
     &                        - t1a(b-no(spin),k)*
     &                          int_acij(a,c,i,j)
     &                        + t1a(a-no(spin),i)*
     &                          int_bcjk(b,c,j,k)
     &                        - t1a(b-no(spin),i)*
     &                          int_acjk(a,c,j,k)
     &                        - t1a(a-no(spin),j)*
     &                          int_bcik(b,c,i,k)
     &                        + t1a(b-no(spin),j)*
     &                          int_acik(a,c,i,k)
c
                        d = evals(a)+
     &                      evals(b)+
     &                      evals(c)-
     &                      evals(i)-
     &                      evals(j)-
     &                      evals(k)
c
                        energy = energy - 
     &                           v*w(k,j,i,c,b,a)/d
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
      write(6,*) ' E(2) block', energy
c
      return
      end
      subroutine w_v_d_energy
     &  (v, w, evals, energy,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi)
c
      implicit none
#include "cuccsdtP.fh"
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi
c
      double precision v(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision w(klo:khi,jlo:jhi,ilo:ihi,
     &                   clo:chi,blo:bhi,alo:ahi)
      double precision evals(nmo), energy, d
      integer i, j, k, a, b, c
c
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
            do c = clo, min(chi,b-1)
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                     do k = klo, min(khi,j-1)
c
c                        write(6,*)' a,b,c,i,j,k,e: ',
c     &                              a,b,c,i,j,k,e
c
                        d = evals(a)+
     &                      evals(b)+
     &                      evals(c)-
     &                      evals(i)-
     &                      evals(j)-
     &                      evals(k)
c
                        energy = energy - 
     &                           v(k,j,i,c,b,a)*w(k,j,i,c,b,a)/d
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
      write(6,*) ' E(2) block', energy
c
      return
      end
      subroutine uccsdt_solver(n,m,x,df,delta)
      implicit none
#include "errquit.fh"
      integer n, m
      double precision x(n,m), df(n,m), delta(n)
c
      integer maxm
      parameter (maxm=100)
      double precision a(maxm,maxm), b(maxm), c(maxm)
      integer ipiv(maxm)
      integer i, k, info
      double precision ddot
      external ddot
c
      if (m .lt. 1) call errquit('m?', m, UNKNOWN_ERR)
      if (m .gt. maxm) call errquit('maxm', m, UNKNOWN_ERR)
c
      call dcopy(n, df(1,m), 1, delta, 1)
      call dscal(n, -1d0, delta, 1)
c
      if (m .eq. 1) return      ! use the jacobi update
c
      do i = 1, m-1
         do k = 1, m-1
            a(i,k) = 
     $           + ddot(n,x(1,i),1,df(1,k),1)
     $           - ddot(n,x(1,m),1,df(1,k),1)
     $           - ddot(n,x(1,i),1,df(1,m),1)
     $           + ddot(n,x(1,m),1,df(1,m),1)
         end do
         b(i) = 
     $        - ddot(n,x(1,i),1,df(1,m),1)
     $        + ddot(n,x(1,m),1,df(1,m),1)
      end do
c
c$$$      write(6,*) ' REDUCED A'
c$$$      call doutput(a, 1, m-1, 1, m-1, maxm, maxm, 1)
c$$$      write(6,*) ' REDUCED B'
c$$$      call doutput(b, 1, m-1, 1, 1, m-1, 1, 1)
c
      call dgesv(m-1, 1, a, maxm, ipiv, b, maxm, info)
      if (info .ne. 0) call errquit(' info ?  ', 0, UNKNOWN_ERR)
      call dcopy(m-1, b, 1, c, 1)
c$$$      write(6,*) ' REDUCED C'
c$$$      call doutput(b, 1, m-1, 1, 1, m-1, 1, 1)
c
      do k = 1, m-1
         call daxpy(n, c(k), x(1,k), 1, delta, 1)
         call daxpy(n,-c(k), x(1,m), 1, delta, 1)
         call daxpy(n,-c(k),df(1,k), 1, delta, 1)
         call daxpy(n, c(k),df(1,m), 1, delta, 1)
      end do
c
      end
      subroutine write_t(
     $     noa, nob, nva, nvb, 
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab)
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "inp.fh"
      integer noa, nob, nva, nvb
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
      character*256 actualname  ! Name of file
      integer fd                ! CHEMIO fd
      double precision offset
      integer ierr
      character*80 errmsg
      actualname = 't_file'
      if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     $   call errquit('write_t: eaf_open failed', 0, DISK_ERR)
      offset = 0
      ierr = eaf_write(fd, offset, t1a, 8*nva*noa) ! write t1a
      offset = offset + 8*nva*noa
      ierr = eaf_write(fd, offset, t1b, 8*nvb*nob) ! write t1b
      offset = offset + 8*nvb*nob
      ierr = eaf_write(fd, offset, t2aa, 8*(nva*noa)**2) ! write t2aa
c      call jan_debug_print('t2aa', t2aa, nva, nva, noa, noa)
      offset = offset + 8*(nva*noa)**2
      ierr = eaf_write(fd, offset, t2bb, 8*(nvb*nob)**2) ! write t2bb
c      call jan_debug_print('t2bb', t2bb, nvb, nvb, nob, nob)
      offset = offset + 8*(nvb*nob)**2
      ierr = eaf_write(fd, offset, t2ab, 8*nva*noa*nvb*nob) ! write t2ab
c      call jan_debug_print('t2ab', t2ab, nva, nvb, noa, nob)
      if (ierr .ne. 0) then
         call eaf_errmsg(ierr, errmsg)
         write(6,*) ' IO offset ', offset
         write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
         call errquit(' write_t write: write failed',0, DISK_ERR)
      end if
      return
      end
      subroutine read_t(
     $     noa, nob, nva, nvb, 
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab)
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "inp.fh"
      integer noa, nob, nva, nvb
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nva, nva, noa, noa)
      double precision t2bb(nvb, nvb, nob, nob)
      double precision t2ab(nva, nvb, noa, nob)
      character*256 actualname  ! Name of file
      integer fd                ! CHEMIO fd
      double precision offset
      integer ierr
      character*80 errmsg
      actualname = 't_file'
      if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     $   call errquit('read_t: eaf_open failed', 0, DISK_ERR)
      offset = 0
      ierr = eaf_read(fd, offset, t1a, 8*nva*noa) ! read t1a
      offset = offset + 8*nva*noa
      ierr = eaf_read(fd, offset, t1b, 8*nvb*nob) ! read t1b
      offset = offset + 8*nvb*nob
      ierr = eaf_read(fd, offset, t2aa, 8*(nva*noa)**2) ! read t2aa
      offset = offset + 8*(nva*noa)**2
      ierr = eaf_read(fd, offset, t2bb, 8*(nvb*nob)**2) ! read t2bb
      offset = offset + 8*(nvb*nob)**2
      ierr = eaf_read(fd, offset, t2ab, 8*nva*noa*nvb*nob) ! read t2ab
      if (ierr .ne. 0) then
         call eaf_errmsg(ierr, errmsg)
         write(6,*) ' IO offset ', offset
         write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
         call errquit(' read_t read: read failed',0, DISK_ERR)
      end if
      return
      end
      subroutine read_dummy_t(nbf,
     $     noa, nob, nva, nvb, 
     $     t1a, t1b,
     $     t2aa,t2bb,t2ab)
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "inp.fh"
#include "mafdecls.fh"
      integer noa, nob, nva, nvb, nbf
      double precision t1a(nva, noa)
      double precision t1b(nvb, nob)
      double precision t2aa(nbf, nbf, nbf, nbf)
      double precision t2bb(nbf, nbf, nbf, nbf)
      double precision t2ab(nbf, nbf, nbf, nbf)
      character*256 actualname  ! Name of file
      integer fd                ! CHEMIO fd
      double precision offset
      integer ierr
      integer k_temp, l_temp
      character*80 errmsg
      actualname = 't_file'
      if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     $   call errquit('read_t: eaf_open failed', 0, DISK_ERR)
c
c      if (.not. ma_push_get(mt_dbl, (noa*nva), 
c     $    't1_temp', l_temp, k_temp))
c     $     call errquit('ma t1a', (noa*nva))
      offset = 0
      ierr = eaf_read(fd, offset, t1a, 8*nva*noa) ! read t1a
c      ierr = eaf_read(fd, offset, dbl_mb(k_temp), 8*nva*noa) ! read t1a
c      call repack_t1 (dbl_mb(k_temp), t1a, nva, noa, nbf)
c      if (.not. ma_pop_stack(l_temp)) call errquit('t1_temp?',0)
c
c      if (.not. ma_push_get(mt_dbl, (nob*nvb), 
c     $    't1_temp', l_temp, k_temp))
c     $     call errquit('ma t1b', (nob*nvb))
      offset = offset + 8*nva*noa
      ierr = eaf_read(fd, offset, t1b, 8*nvb*nob) ! read t1b
c      ierr = eaf_read(fd, offset, dbl_mb(k_temp), 8*nvb*nob) ! read t1b
c      call repack_t1 (dbl_mb(k_temp), t1b, nvb, nob, nbf)
c      if (.not. ma_pop_stack(l_temp)) call errquit('t1_temp?',0)
c
      if (.not. ma_push_get(mt_dbl, (noa*nva)**2, 
     $    't2_temp', l_temp, k_temp))
     $     call errquit('ma t2_aa', (noa*nva)**2, MA_ERR)
      offset = offset + 8*nvb*nob
      ierr = eaf_read(fd, offset, dbl_mb(k_temp), 8*(nva*noa)**2) ! read t2aa
      call repack_t2 (dbl_mb(k_temp), t2aa, nva, nva, noa, noa, nbf)
      if (.not. ma_pop_stack(l_temp)) call errquit('t1_temp?',0,
     &       MA_ERR)
c
      if (.not. ma_push_get(mt_dbl, (nob*nvb)**2, 
     $     't2_temp', l_temp, k_temp))
     $     call errquit('ma t2_bb', (nob*nvb)**2, MA_ERR)
      offset = offset + 8*(nva*noa)**2
      ierr = eaf_read(fd, offset, dbl_mb(k_temp), 8*(nvb*nob)**2) ! read t2bb
      call repack_t2 (dbl_mb(k_temp), t2bb, nvb, nvb, nob, nob, nbf)
      if (.not. ma_pop_stack(l_temp)) call errquit('t1_temp?',0, MA_ERR)
c
      if (.not. ma_push_get(mt_dbl, (noa*nva)*(nob*nvb), 
     $     't2_temp', l_temp, k_temp))
     $     call errquit('ma t2_ab', (noa*nva)*(nob*nvb), MA_ERR)
      offset = offset + 8*(nvb*nob)**2
      ierr = eaf_read(fd, offset, dbl_mb(k_temp), 8*nva*noa*nvb*nob) ! read t2ab
      call repack_t2 (dbl_mb(k_temp), t2ab, nva, nvb, noa, nob, nbf)
      if (.not. ma_pop_stack(l_temp)) call errquit('t1_temp?',0, MA_ERR)
c
      if (ierr .ne. 0) then
         call eaf_errmsg(ierr, errmsg)
         write(6,*) ' IO offset ', offset
         write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
         call errquit(' read_t read: read failed',0, DISK_ERR)
      end if
      return
      end
      subroutine repack_t1(
     $     rd_buf, sv_buf, nv1, no1, nbf)
      implicit none
      integer nv1, no1, nbf
      integer nva, noa
      double precision rd_buf(nv1, no1)
      double precision sv_buf(nbf, nbf)
      call dfill(nbf**2,0.0d0,sv_buf,1)
      do nva = 1, nv1
         do noa = 1, no1
            sv_buf(no1+nva,noa) = rd_buf(nva,noa)
            sv_buf(noa,no1+nva) = rd_buf(nva,noa)
         end do
      end do
      return
      end
      subroutine repack_t2(
     $     rd_buf, sv_buf, nv1, nv2, no1, no2, nbf)
      implicit none
      integer nv1, nv2, no1, no2, nbf
      integer nva, noa, nvb, nob
      double precision rd_buf(nv1, nv2, no1, no2)
      double precision sv_buf(nbf, nbf, nbf, nbf)
      call dfill(nbf**4,0.0d0,sv_buf,1)
c      call jan_debug_print('rd_buf', rd_buf, nv1, nv2, no1, no2)
      do nva = 1, nv1
         do nvb = 1, nv2
            do noa = 1, no1
               do nob = 1, no2
c
                  sv_buf(no1+nva,no2+nvb,noa,nob) = 
     &                rd_buf(nva,nvb,noa,nob)
c
                  sv_buf(noa,nob,no1+nva,no2+nvb) = 
     &                rd_buf(nva,nvb,noa,nob)
c
               end do
            end do
         end do
      end do
c      call jan_debug_print('sv_buf', sv_buf, nbf, nbf, nbf, nbf)
      return
      end
