      subroutine dft_uks_energy( rtdb, g_vecs, eone, etwo, exc, enrep,
     $                           energy, edisp, g_grad, nexc )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "util.fh"
#include "cscfps.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "case.fh"
c     
c     $Id$
c
      integer g_vecs(2)
      double precision energy
      integer g_grad
      integer rtdb
c     
      double precision eone, etwo, enrep, Exc(2)
      integer nExc
      integer gtype, grow, gcol
      integer d(4), f(8), nfock
      integer g_a_dens, g_a_coul, g_a_exch, g_a_xc
      integer g_b_dens, g_b_coul, g_b_exch, g_b_xc
      integer g_hcore
      integer g_tmp(6),ifock
      double precision jfac(4), kfac(4), one, zero, mone
      parameter (one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision e_a_coul, e_a_exch, e_b_coul, e_b_exch,
     &     e_a_xc, e_b_xc
      double precision errmaxa, errmaxb
      double precision Edisp, dum, Ebq
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical odebug
      logical disp
      logical cphf_uhf
      logical xc_gotxc
      external xc_gotxc
      integer dims(3),chunk(3)
c
      double precision xc_hfexch
      external xc_hfexch
c
      data f/8*0/
c     
c     Check
c     
      odebug = util_print('uks_debug', print_debug)
      call uhf_jkfac(jfac,kfac)
      if (.not.cuhf_init_flag)
     $     call errquit('dft_uks_energy: UKS internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.cuhf_vlen).or.(gcol.ne.1))
     $     call errquit('dft_uks_energy: invalid vector length',0,
     $                  GA_ERR)
      cphf_uhf = .false.
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_uhf', mt_log, 1, cphf_uhf)) then
         cphf_uhf = .false.
      endif
c     
c     Arrays for AO density, coulomb and exchange matrices
c
      g_a_coul = ga_create_atom_blocked(geom, basis, 'uks:a coul')
      g_b_coul = ga_create_atom_blocked(geom, basis, 'uks:b coul')
      g_a_exch = ga_create_atom_blocked(geom, basis, 'uks:a exch')
      g_b_exch = ga_create_atom_blocked(geom, basis, 'uks:b exch')
      if(cphf_uhf .or. xc_gotxc())then
         g_a_xc   = ga_create_atom_blocked(geom, basis, 'uks:a xc')
         g_b_xc   = ga_create_atom_blocked(geom, basis, 'uks:b xc')
      endif
      g_a_dens = ga_create_atom_blocked(geom, basis, 'uks:a dens')
      g_b_dens = ga_create_atom_blocked(geom, basis, 'uks:b dens')
      call ga_zero(g_a_dens)
      call ga_zero(g_b_dens)
c
c     Make the densites and build the fock matrices
c
      call ga_dgemm('n', 't', nbf, nbf, nalpha, one, g_vecs(1),
     $     g_vecs(1), zero, g_a_dens)
      if (nbeta .gt. 0) then
         call ga_dgemm('n', 't', nbf, nbf, nbeta, one, g_vecs(2),
     $        g_vecs(2), zero, g_b_dens)
      else
         call ga_zero(g_b_dens)
      endif
c
c     Since UHF can break spatial symmetry by localizing the orbitals
c     the densities may not be totally symmetric, but since the
c     Hamiltonian is symmetric contraction with the integrals projects
c     out the totally symmetric component ... hence we can symmetrize
c     the densities and exploit symmetry. Compute the max change in any 
c     element due to symmetrizing and print a warning if it is big.
c
c     !! If this is the case then where does the 'force' for symmetry
c     breaking come from?  Must be missing something?
c     !! Yes, the symmetrization does not take to irrep of the wave-
c     function into account. Instead it generates some totally symmetric
c     density, which is wrong if the wavefunction has a different
c     irrep.
c
      call ga_copy(g_a_dens,g_a_coul)
      call ga_copy(g_b_dens,g_b_coul)
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_a_dens)
         if (oscfps) call pstat_off(ps_sym_sym)
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_b_dens)
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
      call ga_dadd(one, g_a_dens, mone, g_a_coul, g_a_coul)
      call ga_dadd(one, g_b_dens, mone, g_b_coul, g_b_coul)
      call ga_maxelt(g_a_coul, errmaxa)
      call ga_maxelt(g_b_coul, errmaxb)
      if (max(errmaxa,errmaxb).gt.1d-4) then
         if (ga_nodeid().eq.0) then
            write(6,77) errmaxa,errmaxb
 77         format(' Warning: spatial symmetry breaking in UKS: ',
     $           1p,2d9.2)
            call util_flush(6)
         endif
      endif
c
      if (odebug) then
         call ga_print(g_vecs(1))
         call ga_print(g_vecs(2))
         call ga_print(g_a_dens)
         call ga_print(g_b_dens)
      endif
 
      call ga_zero(g_a_coul)
      call ga_zero(g_b_coul)
      call ga_zero(g_a_exch)
      call ga_zero(g_b_exch)
      if(cphf_uhf .or. xc_gotxc())then
         call ga_zero(g_a_xc)
         call ga_zero(g_b_xc)
      endif
      d(1) = g_a_dens
      d(2) = g_a_dens
      d(3) = g_b_dens
      d(4) = g_b_dens
      f(1) = g_a_coul
      f(2) = g_a_exch
      f(3) = g_b_coul
      f(4) = g_b_exch
      if(cphf_uhf .or. xc_gotxc())then
         f(5) = g_a_xc
         f(6) = g_b_xc
      endif
c
c     two extra ga's are passed to fock_2e to get the xc matrix
c
      if(cphf_uhf)then
         nfock = 6
      else
         nfock = 4
      endif
      call do_riscf (.false.)
      if (.not.cam_exch) then
        kfac(1)=0d0
        kfac(2)=xc_hfexch()
        kfac(3)=0d0
        kfac(4)=xc_hfexch()
        call fock_2e(geom, basis, nfock, jfac, kfac, tol2e,
     $     oskel, d, f, .false.)
      else ! for attenuated calculations
c
c       calculate the CAM exchange
c
        do ifock = 1,nfock
           g_tmp(ifock) = ga_create_atom_blocked(geom, basis, 'tmp')
        end do
        call ga_zero(g_tmp)
c
        call case_setflags(.true.)
        jfac(1)=0d0
        jfac(2)=0d0
        jfac(3)=0d0
        jfac(4)=0d0
        kfac(1)=0d0
        kfac(2)=xc_hfexch()
        kfac(3)=0d0
        kfac(4)=xc_hfexch()
        call fock_2e_cam(geom, basis, nfock, jfac, kfac, tol2e, 
     &     oskel, d, g_tmp, .false., .false.)
        call ga_dadd(1d0,f,1d0,g_tmp,f)
c
c       calculate the full Coulomb
c
        call case_setflags(.false.)
        jfac(1)=1d0
        jfac(2)=0d0
        jfac(3)=1d0
        jfac(4)=0d0
        kfac(1)=0d0
        kfac(2)=0d0
        kfac(3)=0d0
        kfac(4)=0d0
        call fock_2e_cam(geom, basis, nfock, jfac, kfac, tol2e, 
     &     oskel, d, g_tmp, .false., .true.)
        call ga_dadd(1d0,f,1d0,g_tmp,f)
c
c       destroy work space
        if (.not. ga_destroy(g_tmp)) 
     &   call errquit('uhf: ga corrupt?',0, GA_ERR)
      end if
      call do_riscf (.true.)
c
c     do DFT stuff
c
      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      if (xc_gotxc()) then
        call fock_xc(geom, nbf, basis, 4*2, d, f(5), Exc, nExc)
      endif
c
c     enddo DFT stuff
c
      e_a_coul = 0.5d0*
     $     (ga_ddot(g_a_dens,g_a_coul) + ga_ddot(g_a_dens,g_b_coul))
      e_b_coul = 0.5d0*
     $     (ga_ddot(g_b_dens,g_a_coul) + ga_ddot(g_b_dens,g_b_coul))
      e_a_exch = 0.5d0*ga_ddot(g_a_dens,g_a_exch)
      e_b_exch = 0.5d0*ga_ddot(g_b_dens,g_b_exch)
      if (xc_gotxc()) then
        etwo = e_a_coul + e_b_coul + Edisp
        exc(1) = exc(1) - e_a_exch - e_b_exch
      else
        etwo = e_a_coul + e_b_coul - e_a_exch - e_b_exch + Edisp
      endif
c take this out? it seems that the xc_gotxc part does essentially the
c same thing apart from the fact that the energy is added on properly.
      if(cphf_uhf)then
         e_a_xc = ga_ddot(g_a_dens,g_a_xc)
         e_b_xc = ga_ddot(g_b_dens,g_b_xc)
         etwo = etwo + e_a_xc + e_b_xc        
      endif
c
      if (odebug .and. ga_nodeid().eq.0) then
         write(6,*) ' coulomb energies', e_a_coul, e_b_coul
         write(6,*) ' exchang energies', e_a_exch, e_b_exch
         call util_flush(6)
      endif
      if (odebug) then
         call ga_print(g_a_coul)
         call ga_print(g_a_exch)
      endif
c
c     Form energies and AO fock matrices
c
c     Fa (in g_a_coul) = h + J(a) + J(b) - K(a)
c     Fb (in g_b_coul) = h + J(a) + J(b) - K(b)
c
c     E = ((Da + Db)*h + Da*Fa + Db*Fb) / 2
c     Eone = h * (Da + Db)
c
c     2e denotes 2-electron components only
c
      call ga_dadd(one, g_a_coul, one, g_b_coul, g_a_coul)
      call ga_copy(g_a_coul, g_b_coul)
      call ga_dadd(one, g_a_coul, mone, g_a_exch, g_a_coul)
      call ga_dadd(one, g_b_coul, mone, g_b_exch, g_b_coul)
      if(cphf_uhf.or.xc_gotxc())then
         call ga_dadd(one, g_a_coul, one, g_a_xc, g_a_coul)
         call ga_dadd(one, g_b_coul, one, g_b_xc, g_b_coul)
      endif
c
c     reuse g_a_exch to hold the 1-e integrals
c
      g_hcore = g_a_exch
      if (rtdb_get(rtdb,'dft:bq_energy',mt_dbl,1,Ebq)) then
        call ga_zero(g_hcore)
        call int_1e_ga(basis, basis, g_hcore, 'bq_pot', oskel)
        Ebq = ga_ddot(g_a_dens,g_hcore)+ga_ddot(g_b_dens,g_hcore)
        if (.not. rtdb_put(rtdb,'dft:bq_energy',mt_dbl,1,Ebq))
     $     call errquit('dft_roks_fock: writing bq_energy failed',
     $     0, RTDB_ERR)
      endif
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_hcore, 'potential', oskel)
      eone = 
     $     (ga_ddot(g_a_dens,g_hcore) + ga_ddot(g_b_dens,g_hcore))
      call ga_dadd(one, g_hcore, one, g_a_coul, g_a_coul)
      call ga_dadd(one, g_hcore, one, g_b_coul, g_b_coul)
c
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., g_a_coul)
         if (oscfps) call pstat_off(ps_sym_sym)
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., g_b_coul)
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
      if (odebug) then
         call ga_print(g_a_coul)
         call ga_print(g_b_coul)
      endif
c
c     Transform the Fock matrices to the MO basis using g_a_dens
c     for scratch
c
      call two_index_transf(g_a_coul, g_vecs(1), g_vecs(1),
     $     g_a_dens, cuhf_g_falpha)
      call two_index_transf(g_b_coul, g_vecs(2), g_vecs(2),
     $     g_a_dens, cuhf_g_fbeta)
c
      if (odebug) then
         call ga_print(cuhf_g_falpha)
         call ga_print(cuhf_g_fbeta)
      endif
c
c     Free up dead global arrays
c
      if (.not. ga_destroy(g_a_dens)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_dens)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      if(cphf_uhf .or. xc_gotxc())then
         if (.not. ga_destroy(g_a_xc)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_b_xc)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      endif
      if (.not. ga_destroy(g_a_exch)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_exch)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_a_coul)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_coul)) call errquit('uks_e: destroy',0,
     &       GA_ERR)
c
c     extract the gradient
c
      call uhf_get_grad(g_grad)
c
      if (odebug) call ga_print(g_grad)
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('dft_uks_energy: no repulsion energy?', 0,
     $                  GEOM_ERR)
      energy = eone + etwo + enrep
      if (xc_gotxc()) then
        energy = energy + Exc(1) + Exc(2)
      endif
c
      if (odebug .and. ga_nodeid().eq.0) then
         write(6,*) ' eone, etwo, enrep, energy ',
     $        eone, etwo, enrep, energy
      endif
c
      end
