      subroutine rohf_scat(rtdb, g_vecs, g_scatmo, g_borb,
     &     obbtrn,hbbtrnd)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "cscfps.fh"
c     
c     Build scattering orbitals. The scattering orbitals are obtained
c     from diagonaling the space spanned by the diffuse basis set with
c     the scf core Hamiltonian.
c     
      integer rtdb
      integer g_vecs            ! [input] MO vectors
      integer g_all             ! [output] scattering orbitals
      integer g_tran,g_tmp,g_scat,g_over, g_tmpscat,g_scatmo
c
c     The following variables are defined in different header files
c
c     integer geom, basis       ! 
c     integer nclosed, nopen    ! no. of closed and open shells. see cscf.fh
c     integer nmo               ! no. of MO vectors. see cscf.fh
c     double precision tol2e    ! Integral tolerance see cscf.fh
c
      double precision eone, etwo ! [output] Energy terms
      integer g_fcv             ! [output] Closed-virtual Fock matrix
      integer g_fpv             ! [output] Open-virtual Fock matrix
      integer g_fcp             ! [output] Closed-open Fock matrix
c
c     logical oskel             ! use skeleton symmetry see cscf.fh
c     
c     Local variables
c     
      double precision eone_closed, eone_open, etwo_closed, etwo_open
      double precision errmax
      integer nset
      parameter(nset=3)
      integer g_cdens, g_odens, g_cfock, g_hcore, g_ocoul, g_oexch 
      integer nvir, nfock, nscat
c     integer nbf               ! see cscf.fh
      double precision two, one, zero, mone
      parameter (two=2.0d0, one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision sym_tol  ! Threshold for detecting symmetry breaking
      double precision jfactor(nset),kfactor(nset)
      integer iv_dens(nset),iv_fock(nset)
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oforcesym
      logical movecs_read, movecs_write, movecs_read_header
      external movecs_read, movecs_write, movecs_read_header
      logical int_normalize
      external int_normalize
c
c     Variables for vector projection
c
      integer natoms, iat, ihi, ilo, jlo, jhi, i, j
      double precision temp
c     
c     variables for reading movecs
c
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs, kohn_movecs
      integer nbf_vecs, nsets, mbf, noc, nbf_all
c     integer nmo   ! see cscf.fh 
      logical status, all_bas_han
c     integer k_occ, l_occ
      integer k_evals, l_evals, l_ran, k_ran, l_vec, k_vec
c
      data oforcesym /.false./
c
      integer g_cos1e
      dimension g_cos1e(3)
      logical odbug
      logical osome
c
      integer itype_wfn
c
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_prep2 , k_txs_prep2 
c
      integer nchnl,nbfmax,nstate,istat,ibf,nstat
      parameter (nchnl=20,nbfmax=130,nstate=nchnl**2)
      integer nsc_orb(nchnl), nsch(nbfmax,nchnl)
      double precision obbtrn(nbfmax,nbfmax,nstate)
      double precision hbbtrnd(nbfmax,nbfmax,nstate)
      integer g_tr, g_tmpfock, g_bvb, g_borb
c
      integer nswap, swap(20)
c
      data jfactor / 1.0d0, 1.0d0, 0.0d0/
      data kfactor /-0.5d0, 0.0d0,-0.5d0/
c     
c     Initialise
c     
      status = rtdb_parallel(.true.)
      call scf_init(rtdb)
      call rohf_init(rtdb)
      if (.not.int_normalize(rtdb,basis))
     $     call errquit('rohf_scat:int_normalize failed', 0, INT_ERR)
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_scat: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
c     
      if (.not. bas_create(all_bas_han, 'all'))
     &     call errquit('rohf_scat: bas_create?', 86, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, all_bas_han, 'all'))
     &     call errquit('rohf_scat: no cckohn basis ', 86, BASIS_ERR)
      if (.not. bas_numbf(all_bas_han, nbf_all))
     &   call errquit('rohf_scat: basis set error:', 86, BASIS_ERR)
      nscat = nbf_all - nbf
      if(.not.bas_print(all_bas_han))
     &     call errquit('cant print basis in rohf_sct', 1, BASIS_ERR)
      if (.not. int_normalize(rtdb,all_bas_han))
     &     call errquit('rohf_scat: int_normalize failed', 86, INT_ERR)
      call int_init(rtdb,1,all_bas_han)
      if(.not.bas_print(all_bas_han))
     &     call errquit('cant print basis in rohf_sct', 1, BASIS_ERR)
      call schwarz_init(geom, all_bas_han)
c     
c     Read the converged scf movecs
c     
      g_vecs = ga_create_atom_blocked(geom,basis,'Molecular orbital')
      if (.not. movecs_read_header(movecs_out, title_vecs,
     $     basis_vecs, scftype_vecs, nbf_vecs, nsets, nmo, 2))
     $     call errquit('rohf_scat: cannot read MO vectors header',0,
     &       UNKNOWN_ERR)
      status = MA_push_get(MT_DBL, nbf_all, 'Evals', l_evals, k_evals)
      status = ma_push_get(MT_DBL, nbf_all, 'MO occ', l_occ, k_occ)
      call dfill(nbf_all    , 0.0d0, dbl_mb(k_occ), 1)
      if (.not. movecs_read(movecs_out, 1, dbl_mb(k_occ),
     $     dbl_mb(k_evals), g_vecs)) call
     $     errquit('rohf_scat: cannot read MO vectors',0,
     &       UNKNOWN_ERR)
c
      call ga_orthog_mos(basis, g_vecs)
c
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('rohf_scat: bas_numbf failed',0, BASIS_ERR)
c     
c     Project the MOs to the whole orbital space
c
      if(.not.ga_create(MT_DBL, nbf_all, nbf_all, 'Total MO',
     &     nbf_all, 0, g_all))
     &     call errquit('rohf_scat: cant create g_all', nbf_all, GA_ERR)
      if(.not.ga_create(MT_DBL,nbf_all,nscat,'rohf:tran',
     &     nbf_all,0,g_tran))
     &     call errquit('can not create transformation matrix',0,
     &       GA_ERR)
      if(.not.ga_create(MT_DBL,nbf_all,nscat,'rohf:tmp',
     &     nbf_all,0,g_tmp))
     &     call errquit('can not create tmp matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nscat,nscat,'rohf:tmp scat',
     &     nscat,0,g_tmpscat))
     &     call errquit('can not create scat tmp matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nscat,nscat,'rohf:scat',
     &     nscat,0,g_scat))
     &     call errquit('can not create scat matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nbf_all,nscat,'cckohn:scatmo',
     &     nbf_all,0,g_scatmo))
     &     call errquit('can not create scatmo matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nscat,nscat,'rohf:over',
     &     nscat,0,g_over))
     &     call errquit('can not create overlap matrix',0, GA_ERR)
      call ga_zero(g_all)
      call ga_zero(g_tran)
      call ga_zero(g_tmp)
      call ga_zero(g_tmpscat)
      call ga_zero(g_scat)
      call ga_zero(g_scatmo)
      call ga_zero(g_over)
c
c     The following definition is for the transformation matrix, which 
c     includes all the occupied as well as the target and scattering orbitals.
c
      nsc_orb(1)=9
      if(.not.ga_create(MT_DBL,nbf_all,nclosed+nsc_orb(1),
     &     'rohf:g_borb',nbf_all,0,g_borb))
     &     call errquit('can not create bound orbital matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nbf_all,nsc_orb(1),
     &     'rohf:g_tmpfock',nbf_all,0,g_tmpfock))
     &     call errquit('can not create tmp fock matrix',0, GA_ERR)
      if(.not.ga_create(MT_DBL,nbf_all,nsc_orb(1),
     &     'rohf:g_tr',nbf_all,0,g_tr))
     &     call errquit('can not create transformation matrix',0,
     &       GA_ERR)
      if(.not.ga_create(MT_DBL,nsc_orb(1),nsc_orb(1),
     &     'rohf:g_tr',nsc_orb(1),0,g_bvb))
     &     call errquit('can not create transformation matrix',0,
     &       GA_ERR)
      call ga_zero(g_borb)
      call ga_zero(g_tmpfock)
      call ga_zero(g_tr)
      call ga_zero(g_bvb)
      nsch(1,1)=8
      nsch(2,1)=16
      nsch(3,1)=26
      nsch(4,1)=31
      nsch(5,1)=50
      nsch(6,1)=52
      nsch(7,1)=57
      nsch(8,1)=63
      nsch(9,1)=70
c
c     Projected the scf orbitals to the full space. The scattering 
c     orbitals are initilized with random numbers.
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('rohf_scat: geom_ncent ?', 0, GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl, nbf_all, 'diag elements',
     $     l_ran, k_ran)) call errquit
     $     ('rohf_scat: ma diag elements', nbf_all, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'vecs',
     $     l_vec, k_vec)) call errquit
     $     ('rohf_scat: vecs', nbf, MA_ERR)
      do i=1+ga_nodeid(),nbf,ga_nnodes()
         call dfill(nbf_all, 0.0d0, dbl_mb(k_ran), 1)
         call dfill(nbf    , 0.0d0, dbl_mb(k_vec), 1)
         call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_vec), 1)
         do iat = 1, natoms
            if (.not. bas_ce2bfr(basis, iat, ilo, ihi)) 
     &           call errquit('rohf_scat: bas_ce2bfr basis?', 0,
     &       BASIS_ERR)
            if (.not. bas_ce2bfr(all_bas_han, iat, jlo, jhi)) 
     &           call errquit('rohf_scat: bas_ce2bfr all_bas_han?', 0,
     &       BASIS_ERR)
            if(ihi-ilo+1.gt.0)then
               do j = 0, ihi-ilo
                  dbl_mb(k_ran+jlo+j-1)=dbl_mb(k_vec+ilo+j-1)
               enddo
            endif
         enddo
         call ga_put(g_all, 1, nbf_all, i, i, dbl_mb(k_ran), 1)
      enddo
c     
c      call ga_copy_patch('n',g_vecs,1,nbf,1,nbf,g_all,1,nbf,1,nbf)
c      if (.not. ma_push_get(mt_dbl, nbf_all, 'diag elements',
c     $     l_ran, k_ran)) call errquit
c     $     ('rohf_scat: ma diag elements', nbf_all)
c
      temp = util_random(ga_nodeid()*314159+1)
      do i = nbf+1+ga_nodeid(), nbf_all, ga_nnodes()
         do j = 1, nbf_all
            dbl_mb(k_ran+j-1) = (util_random(0) - 0.5d0)
         enddo
         call ga_put(g_all, 1, nbf_all, i, i, dbl_mb(k_ran), 1)
      end do
      if (.not. ma_pop_stack(l_vec)) call errquit
     $     ('rohf_scat: ma pop l_vec?', 0, MA_ERR)
      if (.not. ma_pop_stack(l_ran)) call errquit
     $     ('rohf_scat: ma pop l_ran?', 0, MA_ERR)
c
      call ga_orthog_mos(all_bas_han, g_all)
c
      call ga_copy_patch('n',g_all,1,nbf_all,nbf+1,nbf_all,
     &     g_tran,1,nbf_all,1,nscat)
c
      nvir = nmo - nclosed - nopen
c     
c     hcore = 1-e integrals
c     cdens = closed shell AO density
c     cfock = closed shell AO fock
c     odens = open shell density
c     ocoul = open shell coulomb
c     oexch = open shell exchange
c     
      g_cdens = 
     &     ga_create_atom_blocked(geom, all_bas_han, 'rohf: cdens')
      g_cfock = 
     &     ga_create_atom_blocked(geom, all_bas_han, 'rohf: cfock')
      g_hcore = 
     &     ga_create_atom_blocked(geom, all_bas_han, 'rohf: hcore')
      if (nopen .gt. 0) then
         g_odens = 
     &        ga_create_atom_blocked(geom, all_bas_han, 'rohf: odens')
         g_ocoul = 
     &        ga_create_atom_blocked(geom, all_bas_han, 'rohf: oexch')
         g_oexch = 
     &        ga_create_atom_blocked(geom, all_bas_han, 'rohf: ocoul')
      endif
c     
c     Compute closed and open shell density matrices
c     
      call rohf_dens(g_all, nbf_all, nclosed, nopen, g_cdens, g_odens)
      if (util_print('rohf_fock',print_debug)) then
         call ga_print(g_vecs)
         call ga_print(g_cdens)
         if (nopen .gt. 0) call ga_print(g_odens)
      endif
c
c     Debug ... verify that the density is symmetric
c
      if (oforcesym .and. oskel) then
         call ga_copy(g_cdens, g_hcore)
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_cdens)
         if (oscfps) call pstat_off(ps_sym_sym)
         call ga_dadd(one, g_hcore, mone, g_cdens, g_hcore)
         call ga_maxelt(g_hcore, errmax)
         sym_tol = tol2e*10.0d0
         if (errmax .gt. sym_tol) then
            write(6,*) ' Symmetry error ', errmax
            call  errquit
     $           ('rohf_fock:closed density is not totally symmetric',0,
     &       UNKNOWN_ERR)
         endif
         if (nopen .gt. 0) then
            call ga_copy(g_odens, g_hcore)
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .true., g_odens)
            if (oscfps) call pstat_off(ps_sym_sym)
            call ga_dadd(one, g_hcore, mone, g_odens, g_hcore)
            call ga_maxelt(g_hcore, errmax)
            if (errmax .gt. sym_tol) then
               write(6,*) ' Symmetry error ', errmax
               call errquit
     $              ('rohf_fock:open density not totally symmetric',0,
     &       UNKNOWN_ERR)
            endif
         endif
      endif
c
c     Form the overlap matrix g_over
c
      call ga_zero(g_over)
      if (.not. ma_push_get(mt_dbl, nscat, 'diag elements',
     $     l_ran, k_ran)) call errquit
     $     ('rohf_scat: ma diag elements', nbf_all, MA_ERR)
      do i = 1+ga_nodeid(), nscat, ga_nnodes()
         call dfill(nscat, 0.0d0, dbl_mb(k_ran), 1)
         dbl_mb(k_ran+i-1)=1.0
         call ga_put(g_over, 1, nscat, i, i, dbl_mb(k_ran), 1)
      end do
      if (.not. ma_pop_stack(l_ran)) call errquit
     $     ('rohf_scat: ma pop?', 0, MA_ERR)
c
c     One-electron integrals and energies
c     
      call ga_zero(g_hcore)
c      call int_1e_ga(all_bas_han,all_bas_han, 
c     &     g_hcore, 'kinetic', oskel)
      call int_1e_ga(all_bas_han,all_bas_han, 
     &     g_hcore, 'kinetic', .false.)
c      call int_1e_ga(all_bas_han,all_bas_han,
c     &     g_hcore,'potential', oskel)
      call int_1e_ga(all_bas_han,all_bas_han,
     &     g_hcore,'potential', .false.)
      if (util_print('rohf_scat',print_debug)) call ga_print(g_hcore)
      call ga_sync()
c
      eone_closed = ga_ddot(g_hcore, g_cdens)
      eone_open = 0.0d0
      if (nopen .gt. 0) eone_open = ga_ddot(g_hcore, g_odens)
      eone = eone_closed + eone_open
c
c     Form the (one electron?) core Hamiltonian in the scattering space
c
      call ga_dgemm('n','n',nbf_all,nscat,nbf_all,
     &     1.0d0, g_hcore, g_tran, 0.0d0, g_tmp)
      call ga_dgemm('t','n',nscat,nscat,nbf_all,
     &     1.0d0, g_tran, g_tmp, 1.0d0, g_scat)
c
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call ga_pdsygv(g_scat, g_over, g_tmpscat, dbl_mb(k_evals+nbf))
#else
      call ga_diag(g_scat, g_over, g_tmpscat, dbl_mb(k_evals+nbf))
#endif
#else 
      call ga_diag_seq(g_scat, g_over, g_tmpscat, dbl_mb(k_evals+nbf))
#endif
c
c     Transform the MOs to atomic basis
c
      call ga_dgemm('n','n',nbf_all,nscat,nscat,
     &     1.0d0, g_tran, g_tmpscat, 0.0d0, g_scatmo)
      call ga_copy_patch('n',g_scatmo,1,nbf_all,1,nscat,
     &     g_all,1,nbf_all,nbf+1,nbf_all)
      call ga_orthog_mos(all_bas_han, g_all)
      if(ga_nodeid().eq.0)then
         write(*,*)"eigenvalues of the scattering space"
         write(*,"(6f10.5)")(dbl_mb(k_evals+nbf+i-1),
     &        i=1,nscat)
         call ga_print(g_scatmo)
         call ga_print(g_all)
      endif
      if (.not. ma_push_get(mt_int, nbf_all, 'irreps',
     $     l_irs, k_irs)) call errquit('crap',0, MA_ERR)
      if(oskel)         
     &     call scf_movecs_sym_adapt(all_bas_han, g_all, .true.,
     &     nbf_all, '- final', .true., int_mb(k_irs))
      if (.not. ma_pop_stack(l_irs)) call errquit
     $     ('rohf_scat: ma pop?', 0, MA_ERR)
      write(*,*)"final movecs vectors used"
      call ga_print(g_all)
c
c     Form the transformation and bound orbital matrix
c
      call ga_copy_patch('n',g_all,1,nbf_all,1,nclosed,
     &     g_borb,1,nbf_all,1,nclosed)
      do i = 1, nsc_orb(1)
         call ga_copy_patch('n',g_all,1,nbf_all,nsch(i,1),nsch(i,1),
     &        g_tr,1,nbf_all,i,i)
      enddo
      call ga_copy_patch('n',g_tr,1,nbf_all,1,nsc_orb(1),
     &     g_borb,1,nbf_all,nclosed+1,nclosed+nsc_orb(1))
      call ga_print(g_tr)
      call ga_print(g_borb)
      write(*,*)"swap the mo's and write to file"
      nswap = 16
      swap(1)=9
      swap(2)=16
      swap(3)=10
      swap(4)=26
      swap(5)=11
      swap(6)=31
      swap(7)=12
      swap(8)=50
      swap(9)=13
      swap(10)=52
      swap(11)=14
      swap(12)=57
      swap(13)=15
      swap(14)=63
      swap(15)=16
      swap(16)=70
      if (.not. rtdb_put(rtdb, "kohn:swap", mt_int, nswap,swap))
     $     call errquit('rohf_scat: put of swap?', nswap, RTDB_ERR)
      if(ga_nodeid().eq.0)then
         write(*,*)"eigenvalues before swap"
         write(*,"(6f10.5)")(dbl_mb(k_evals+i-1),
     &        i=1,nbf+nscat)
         write(*,*)"occupation numbers before swap"
         write(*,"(6f10.5)")(dbl_mb(k_occ+i-1),
     &        i=1,nbf+nscat)
         call ga_print(g_all)
      endif
      call movecs_swap(rtdb,"kohn","kohn",g_all,
     &     dbl_mb(k_occ),dbl_mb(k_evals))
      write(*,*)"after swap"
      call ga_print(g_all)
      if(ga_nodeid().eq.0)then
         write(*,*)"eigenvalues after swap"
         write(*,"(6f10.5)")(dbl_mb(k_evals+i-1),
     &        i=1,nbf+nscat)
         write(*,*)"occupation numbers after swap"
         write(*,"(6f10.5)")(dbl_mb(k_occ+i-1),
     &        i=1,nbf+nscat)
      endif
      call util_file_name('kohn_movecs', .true.,.false.,kohn_movecs)
      if (.not. movecs_write(rtdb, all_bas_han, kohn_movecs, 
     $     'kohn', 'kohn scf orbitals', 
     $     nbf_all, 1, nbf_all,
     $     dbl_mb(k_occ), nbf_all, dbl_mb(k_evals),
     $     nbf_all, g_all)) call errquit
     $     ('rohf_scat: kohn_movecs write failed', 0, UNKNOWN_ERR)
c     
c     Generate two-electron closed- and open-shell Fock components
c     
c     Fock              1.0*F[ D ]      (g_cfock)
c     
c     Open coulomb      1.0*J[ D ]      (g_ocoul) 
c     
c     Open exchange    -0.5*K[ D ]      (g_oexch) 
c     
      iv_dens(1) =  g_cdens
      iv_fock(1) =  g_cfock
      call ga_zero(g_cfock)
      nfock = 1
      if (nopen .gt. 0) then
         iv_dens(2) =  g_odens
         iv_dens(3) =  g_odens
         iv_fock(2) =  g_ocoul
         iv_fock(3) =  g_oexch
         call ga_zero(g_oexch)
         call ga_zero(g_ocoul)
         nfock = 3
      endif
c
C     don't do it if RI HESSIAN is active
c
      call do_riscf (.false.)
c
C     Dyalls method for RHF
c
      if ( nopen .eq. 0 ) call riscf_dovecs(.true., g_all)
c
      kfactor(1)=0.0
c
      call fock_2e(geom, all_bas_han, nfock, jfactor, kfactor,
     $     tol2e, .false., iv_dens, iv_fock, .false.)
      call do_riscf (.true.)
      if ( nopen .eq. 0 ) call riscf_dovecs(.false., g_all)
c     
c     Energy evaluation
c     
      write(*,*)"cdens and cfock"
c
c      call ga_print(g_cdens)
c      call ga_print(g_cfock)
c
      etwo_closed = 0.5d0*ga_ddot(g_cdens,g_cfock)
      etwo_open = 0.0d0
      if (nopen .gt. 0) etwo_open = 
     $     ga_ddot(g_cdens,g_ocoul) + 
     $     ga_ddot(g_cdens,g_oexch) +
     $     0.5d0*ga_ddot(g_odens,g_ocoul) +
     $     ga_ddot(g_odens,g_oexch)
      etwo = etwo_closed + etwo_open
c     
      call ga_dadd(One, g_hcore, One, g_cfock, g_cfock)
      call ga_print(g_cfock)
      call ga_dgemm('n','n',nbf_all,nsc_orb(1),nbf_all,
     &     1.0d0, g_cfock, g_tr, 0.0d0, g_tmpfock)
      call ga_dgemm('t','n',nsc_orb(1),nsc_orb(1),nbf_all,
     &     1.0d0, g_tr, g_tmpfock, 0.0d0, g_bvb)
      call ga_print(g_bvb)
      write(*,*)"The bound bound matrix dimentions are hardwired"
      nstat = 1
      do istat = 1, nstat
         do ibf = 1, nsc_orb(1)
            call ga_get(g_bvb,1,nsc_orb(1),ibf,ibf,
     &           hbbtrnd(1,ibf,istat),1)
            write(6,100)ibf,(hbbtrnd(i,ibf,istat),i=1,nsc_orb(1))
 100        format(i5,6f12.5/(5x,6f12.5))
         enddo
      enddo
c
      call ga_zero(g_hcore)
      call int_1e_ga(all_bas_han,all_bas_han, 
     &     g_hcore, 'overlap', .false.)
      call ga_print(g_hcore)
      call ga_dgemm('n','n',nbf_all,nsc_orb(1),nbf_all,
     &     1.0d0, g_hcore, g_tr, 0.0d0, g_tmpfock)
      call ga_dgemm('t','n',nsc_orb(1),nsc_orb(1),nbf_all,
     &     1.0d0, g_tr, g_tmpfock, 0.0d0, g_bvb)
      call ga_print(g_bvb)
c
      do istat = 1, nstat
         do ibf = 1, nsc_orb(1)
            call ga_get(g_bvb,1,nsc_orb(1),ibf,ibf,
     &           obbtrn(1,ibf,istat),1)
            write(6,100)ibf,(obbtrn(i,ibf,istat),i=1,nsc_orb(1))
         enddo
      enddo
c     
c     Assemble Fock components
c     F(closed-open)    = F(closed) + Fcoul(open)                   (fcp)
c     F(closed-virtual) = F(closed) + Fcoul(open) + Fexch(open)     (fcv)
c     F(open-virtual)   = F(closed) + Fcoul(open) + 2*Fexch(open)   (fpv)
c     
      if (nopen .gt. 0) then
c     
c     Strive for g_cfock=FCV(AO), g_ocoul=FCP(AO), g_oexch=FPV(AO)
c     
         call ga_dadd(one, g_ocoul, one, g_cfock, g_ocoul)
         call ga_dadd(one, g_ocoul, one, g_oexch, g_cfock)
         call ga_dadd(one, g_cfock, one, g_oexch, g_oexch)
      endif
c     
c     Symmetrize the AO components
c     
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., g_cfock)
         if (oscfps) call pstat_off(ps_sym_sym)
         if (nopen .gt. 0) then
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .false., g_ocoul)
            if (oscfps) call pstat_off(ps_sym_sym)
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .false., g_oexch)
            if (oscfps) call pstat_off(ps_sym_sym)
         endif
      endif
c
      if (util_print('rohf_fock',print_debug)) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' Symmetrized combined AO fock matrices '
            call util_flush(6)
         endif
         call ga_print(g_cfock)
         if (nopen .gt. 0) then
            call ga_print(g_ocoul)
            call ga_print(g_oexch)
         endif
      endif
c     
c     Transform Fock matrices to MO representation using hcore as
c     scratch space
c     
      return
      call two_index_transf(g_cfock, g_vecs, g_vecs, g_hcore, g_fcv)
      if (nopen .gt. 0) then
         call two_index_transf(g_ocoul, g_vecs, g_vecs, g_hcore, g_fcp)
         call two_index_transf(g_oexch, g_vecs, g_vecs, g_hcore, g_fpv)
      endif
c     
      if (util_print('rohf_fock',print_debug)) then
         call ga_print(g_fcv)
         if (nopen .gt. 0) then
            call ga_print(g_fcp)
            call ga_print(g_fpv)
         endif
         if (ga_nodeid().eq.0) then
            write(6,77) eone_open, eone_closed, etwo_open, etwo_closed
 77         format(1p,'  e1o=',d16.8,', e1c=',d16.8,', e2o=',d16.8,
     $           ', e2c=',d16.8)
            call util_flush(6)
         endif
      endif
c     
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_cdens)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_cfock)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_ocoul)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
         if (.not. ga_destroy(g_oexch)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
         if (.not. ga_destroy(g_odens)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
      endif
c     
      end











       



c $Id$
