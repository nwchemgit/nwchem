      subroutine rohf_fock(rtdb, geom, basis, nclosed, nopen, nmo,
     $     tol2e, g_vecs,
     $     eone, etwo, g_fcv, g_fpv,
     $     g_fcp, oskel, iter)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "cscfps.fh"
#include "cosmo.fh"
#include "case.fh"
#include "zora.fh"
#include "crohf.fh"
#include "frozemb.fh"
c     
c     Build fock matrices and energies for ROHF wavefunction
c     
c     Arguments
c     
      integer rtdb
      integer geom, basis       ! [input] Handles to environment
      integer nclosed, nopen    ! [input] no. of closed and open shells
      integer nmo               ! [input] no. of MO vectors
      integer g_vecs            ! [input] MO vectors
      double precision tol2e    ! [input] Integral tolerance
      double precision eone, etwo ! [output] Energy terms
      integer g_fcv             ! [output] Closed-virtual Fock matrix
      integer g_fpv             ! [output] Open-virtual Fock matrix
      integer g_fcp             ! [output] Closed-open Fock matrix
      logical oskel             ! [input] use skeleton symmetry
      integer iter              ! [input] iteration number
c     
c     Local variables
c     
      double precision eone_closed, eone_open, etwo_closed, etwo_open
      double precision errmax
      integer nset
      parameter(nset=3)
      integer g_cdens, g_odens, g_cfock, g_hcore, g_ocoul, g_oexch
      integer g_dens
      integer nbf, nvir, nfock
      double precision two, one, zero, mone
      parameter (two=2.0d0, one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision sym_tol  ! Threshold for detecting symmetry breaking
      double precision jfactor(nset),kfactor(nset)
      integer iv_dens(nset),iv_fock(nset)
      integer g_tmp(2)
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_cos1e
      dimension g_cos1e(3)
      double precision ecosmo
      logical osome
c
      logical odebug
      integer itype_wfn
      character*255 cosmo_file
c     JEM DIM/QM
      integer g_vdim
      logical ldimqm
      logical noskew_uhf
c
c     Initialise
c     
      odebug = util_print('rohf_fock',print_debug)
c
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('rohf_solve: bas_numbf failed',0, BASIS_ERR)
      call rohf_jkfac(jfactor,kfactor)
      nvir = nmo - nclosed - nopen
c     
c     hcore = 1-e integrals
c     cdens = closed shell AO density
c     cfock = closed shell AO fock
c     odens = open shell density
c     ocoul = open shell coulomb
c     oexch = open shell exchange
c     
      g_cdens = ga_create_atom_blocked(geom, basis, 'rohf: cdens')
      g_cfock = ga_create_atom_blocked(geom, basis, 'rohf: cfock')
      g_hcore = ga_create_atom_blocked(geom, basis, 'rohf: hcore')
c
      if (nopen .gt. 0) then
         g_odens = ga_create_atom_blocked(geom, basis, 'rohf: odens')
         g_ocoul = ga_create_atom_blocked(geom, basis, 'rohf: oexch')
         g_oexch = ga_create_atom_blocked(geom, basis, 'rohf: ocoul')
      endif
c     
c     Compute closed and open shell density matrices
c     
      call rohf_dens(g_vecs, nbf, nclosed, nopen, g_cdens, g_odens)
      if (util_print('rohf_fock',print_debug)) then
         call ga_print(g_vecs)
         call ga_print(g_cdens)
         if (nopen .gt. 0) call ga_print(g_odens)
      endif
c
c     Symmetrize the density matrix.
c     If the HOMO and LUMO orbitals are of the same irrep and degenerate
c     then they can rotate in arbitrary ways among eachother. As a
c     result it cannot be guaranteed that a symmetric density matrix
c     can be calculated even if the orbitals are symmetry adapted.
c     Hence density matrix symmetrization is essential. This phenomenon
c     was observed when running the c2h4 test case using D2d symmetry.
c
      if (oskel) then
         call ga_copy(g_cdens, g_hcore)
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_cdens)
         if (oscfps) call pstat_off(ps_sym_sym)
         call ga_dadd(one, g_hcore, mone, g_cdens, g_hcore)
         call ga_maxelt(g_hcore, errmax)
         sym_tol = tol2e*10.0d0
         if (errmax .gt. sym_tol) then
           if (ga_nodeid().eq.0) then
             write(luout,*)'Warning: symmetry breaking in RHF: ',
     &                     errmax
             call util_flush(luout)
           endif
         endif
         if (nopen .gt. 0) then
            call ga_copy(g_odens, g_hcore)
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .true., g_odens)
            if (oscfps) call pstat_off(ps_sym_sym)
            call ga_dadd(one, g_hcore, mone, g_odens, g_hcore)
            call ga_maxelt(g_hcore, errmax)
            if (errmax .gt. sym_tol) then
              if (ga_nodeid().eq.0) then
                write(luout,*)'Warning: symmetry breaking in ROHF: ',
     &                        errmax
                call util_flush(luout)
              endif
            endif
         endif
      endif
c
c     One-electron integrals and energies
c     
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel) ! kinetic
      if (do_zora .and. .not. (do_NonRel)) then
c     is rohf() called from cphf_solve2 (say from hnd_hyperfine_zora)?
c     hack: use 'cphf:skew' to check this condition
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1,
     $        noskew_uhf)) then
            call errquit('ZORA not available yet for SCF ',0,CAPMIS_ERR)
         endif
        call ga_dadd(1.d0,g_hcore,1.d0,g_zora_Kinetic(1),g_hcore) ! zora kinetic
      endif
      call int_1e_ga(basis, basis, g_hcore, 'potential', oskel) ! potential
      if (util_print('rohf_fock',print_debug)) call ga_print(g_hcore)
c
c     cosmo charges, potential and energy contribution
c
      osome=.false.
      osome=osome.and.ga_nodeid().eq.0
      if(cosmo_on.and.cosmo_phase.eq.2) then
       cosmo_file = "cosmo.xyz"
       g_dens = ga_create_atom_blocked(geom,basis,'density matrix')
       call ga_zero(g_dens)
       call ga_copy(g_cdens,g_dens)
       if (nopen.gt.0) call ga_add(1.d0,g_cdens,one,g_odens,g_dens)
        call cosmo_charges_from_dmat(rtdb, basis, geom, ecosmo, odebug,
     &                                  1,           ! 1 = closed  shell
     &                                  g_dens,      ! input density
     &                                  cosmo_file)  ! cosmo charges file name
       call int_1e_ga(basis,basis,g_hcore,'cos_chg_pot',.false.)
       if (.not.ga_destroy(g_dens))
     &  call errquit('dft_roks_fock: ga_destroy failed g_dens',0,GA_ERR)
      endif  ! cosmo check
c
c     add in frozen embedding
c
      if (frozemb) then
        call ga_dadd(1.d0,g_hcore,1.d0,g_frozemb,g_hcore)
      end if
c
c     DIM/QM JEM
c     Add DIM potential
      if (.not.rtdb_get(rtdb, 'dimqm:ldimqm', mt_log, 1, ldimqm))
     $  ldimqm = .false.
      if (ldimqm) then
        call dimqm_rohf_wrap(rtdb, geom, basis, nbf, g_cdens, g_vdim)
        call ga_dadd(one, g_hcore, one, g_vdim, g_hcore)
      end if
c
      call ga_sync()
c
      eone_closed = ga_ddot(g_hcore, g_cdens)
      eone_open = 0.0d0
      if (nopen .gt. 0) eone_open = ga_ddot(g_hcore, g_odens)
      eone = eone_closed + eone_open
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
      call do_riscf (.false.)
C     Dyalls method for RHF
      if ( nopen .eq. 0 ) call riscf_dovecs(.true., g_vecs)
c
      if (.not.cam_exch) then
        call fock_2e(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, iv_fock, .false.)
      else 
c       for attenuated calculations
c
c       get some work space
        g_tmp(1) = ga_create_atom_blocked(geom, basis,'rohf: tmp1')
c
c       calculate the CAM exchange
c
        call ga_zero(g_tmp(1))
        call case_setflags(.true.)
        jfactor(1)=0d0
        jfactor(2)=0d0
        jfactor(3)=0d0
        kfactor(1)=-.5d0
        kfactor(2)=0d0
        kfactor(3)=-.5d0
        call fock_2e_cam(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, g_tmp(1), .false., .false.)
        call ga_dadd(1d0,iv_fock,1d0,g_tmp(1),iv_fock)
c
c       calculate the full Coulomb
c
        call ga_zero(g_tmp(1))
        call case_setflags(.false.)
        jfactor(1)=1d0
        jfactor(2)=1d0
        jfactor(3)=0d0
        kfactor(1)=0d0
        kfactor(2)=0d0
        kfactor(3)=0d0
        call fock_2e_cam(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, g_tmp(1), .false., .true.)
        call ga_dadd(1d0,iv_fock,1d0,g_tmp(1),iv_fock)
c
c       destroy work space
        if (.not. ga_destroy(g_tmp(1))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
      end if
      call do_riscf (.true.)
      if ( nopen .eq. 0 ) call riscf_dovecs(.false., g_vecs)
c     
c     Energy evaluation
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
      call rohf_dumpfock(nopen,g_oexch,g_cfock)
      if (.not. ga_destroy(g_cfock)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
c     DIM/QM JEM
      if (ldimqm) then
        if(.not. ga_destroy(g_vdim))
     $    call errquit('rohf_fock: destroy vdim failed',0,GA_ERR)
      end if
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_ocoul)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
         if (.not. ga_destroy(g_oexch)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
 1       if (.not. ga_destroy(g_odens)) call errquit
     $        ('rohf_fock: ga_destroy failed',0, GA_ERR)
      endif
c     
      end
c
      subroutine rohf_get_fock(g_fock)
      implicit none
#include "cscf.fh"
#include "crohf.fh"
      integer g_fock
c
c     Form the ROHF 'Fock' matrix from FCV, FCP, FPV
c
      call ga_copy(crohf_g_fcv, g_fock)
      if (nopen .gt. 0) then
         call ga_copy_patch('n',
     $        crohf_g_fcp, nclosed+1, nclosed+nopen, 1, nclosed,
     $        g_fock,      nclosed+1, nclosed+nopen, 1, nclosed)
         call ga_copy_patch('n',
     $        crohf_g_fcp, 1, nclosed, nclosed+1, nclosed+nopen,  
     $        g_fock,      1, nclosed, nclosed+1, nclosed+nopen)
         call ga_copy_patch('n',
     $        crohf_g_fpv,nclosed+1,nclosed+nopen,nclosed+nopen+1,nmo,
     $        g_fock,     nclosed+1,nclosed+nopen,nclosed+nopen+1,nmo)
         call ga_copy_patch('n',
     $        crohf_g_fpv,nclosed+nopen+1,nmo,nclosed+1,nclosed+nopen, 
     $        g_fock,     nclosed+nopen+1,nmo,nclosed+1,nclosed+nopen)
      endif
c
      end
