      subroutine mrcc_uhf_energy(g_vecs, eone, etwo, enrep, energy,
     + logzorafromdft)
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
#include "stdio.fh"
#include "zora.fh"


c     
cjb Modyfied code from ddscf/uhf.F for MR CC purposes 
c
      integer g_vecs(2)
      double precision energy
      integer g_grad
ckbn zora
      logical logzorafromdft
c     
      double precision eone, etwo, enrep
      integer gtype, grow, gcol
      integer d(4), f(6), nfock
      integer g_a_dens, g_a_coul, g_a_exch, g_a_xc
      integer g_b_dens, g_b_coul, g_b_exch, g_b_xc
      integer g_hcore
      double precision jfac(4), kfac(4), one, zero, mone
      parameter (one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision e_a_coul, e_a_exch, e_b_coul, e_b_exch,
     &     e_a_xc, e_b_xc
      double precision errmaxa, errmaxb
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical odebug
      logical cphf_uhf
      logical nodezero
c      double precision cpu,wall

      nodezero = (ga_nodeid().eq.0)

cjb Some parameters for one step energy only

      odebug = .false.
      cphf_uhf = .false.
      energy = 0.0d0
      eone = 0.0d0
      etwo = 0.0d0


cjb Weight of energy contributions
ckbn     cpu = - util_cpusec()
ckbn      wall = - util_wallsec()

      call uhf_jkfac(jfac,kfac)
c     
c     Arrays for AO density, coulomb and exchange matrices
c
      g_a_coul = ga_create_atom_blocked(geom, basis, 'uhf:a coul')
      g_b_coul = ga_create_atom_blocked(geom, basis, 'uhf:b coul')
      g_a_exch = ga_create_atom_blocked(geom, basis, 'uhf:a exch')
      g_b_exch = ga_create_atom_blocked(geom, basis, 'uhf:b exch')
      if(cphf_uhf)then
         g_a_xc   = ga_create_atom_blocked(geom, basis, 'uhf:a xc')
         g_b_xc   = ga_create_atom_blocked(geom, basis, 'uhf:b xc')
      endif
      g_a_dens = ga_create_atom_blocked(geom, basis, 'uhf:a dens')
      g_b_dens = ga_create_atom_blocked(geom, basis, 'uhf:b dens')
      call ga_zero(g_a_dens)
      call ga_zero(g_b_dens)

c
c     Make the densites and build the fock matrices
c
ckbn      cpu = cpu + util_cpusec()
ckbn     wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf1", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
      call ga_dgemm('n', 't', nbf, nbf, nalpha, one, g_vecs(1),
     $     g_vecs(1), zero, g_a_dens)
      if (nbeta .gt. 0) then
         call ga_dgemm('n', 't', nbf, nbf, nbeta, one, g_vecs(2),
     $        g_vecs(2), zero, g_b_dens)
      else
         call ga_zero(g_b_dens)
      endif

ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf2", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
c kk ---
ckbn       if(nodezero) then
ckbn         write(6,*)'tce_mrcc_uhf_step1'
ckbn         call util_flush(6)
ckbn       end if
c ------
      call ga_copy(g_a_dens,g_a_coul)
c kk ---
ckbn       if(nodezero) then
ckbn         write(6,*)'tce_mrcc_uhf_step2'
ckbn         call util_flush(6)
ckbn       end if
c ------
      call ga_copy(g_b_dens,g_b_coul)
c kk ---
ckbn       if(nodezero) then
ckbn         write(6,*)'tce_mrcc_uhf_step3'
ckbn         call util_flush(6)
ckbn       end if
c ------
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf3", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
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
 77         format(' Warning: spatial symmetry breaking in UHF: ',
     $           1p,2d9.2)
            call util_flush(6)
         endif
      endif

c
c      if (odebug) then
c         call ga_print(g_vecs(1))
c         call ga_print(g_vecs(2))
c         call ga_print(g_a_dens)
c         call ga_print(g_b_dens)
c      endif

      call ga_zero(g_a_coul)
      call ga_zero(g_b_coul)
      call ga_zero(g_a_exch)
      call ga_zero(g_b_exch)
      if(cphf_uhf)then
         call ga_zero(g_a_xc)
         call ga_zero(g_b_xc)
      endif
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf4", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
      d(1) = g_a_dens
      d(2) = g_a_dens
      d(3) = g_b_dens
      d(4) = g_b_dens
      f(1) = g_a_coul
      f(2) = g_a_exch
      f(3) = g_b_coul
      f(4) = g_b_exch
      if(cphf_uhf)then
         f(5) = g_a_xc
         f(6) = g_b_xc
      endif
      if(cphf_uhf)then
         nfock = 6
      else
         nfock = 4
      endif
      call do_riscf (.false.)
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf4.5", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
      call fock_2e(geom, basis, nfock, jfac, kfac, tol2e,
     $     oskel, d, f, .false.)
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf4.7", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
      call do_riscf (.true.)
c
      e_a_coul = 0.5d0*
     $     (ga_ddot(g_a_dens,g_a_coul) + ga_ddot(g_a_dens,g_b_coul))
      e_b_coul = 0.5d0*
     $     (ga_ddot(g_b_dens,g_a_coul) + ga_ddot(g_b_dens,g_b_coul))
      e_a_exch = 0.5d0*ga_ddot(g_a_dens,g_a_exch)
      e_b_exch = 0.5d0*ga_ddot(g_b_dens,g_b_exch)
      etwo = e_a_coul + e_b_coul - e_a_exch - e_b_exch
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
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf5", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
c      if (odebug) then
c         call ga_print(g_a_coul)
c         call ga_print(g_a_exch)
c      endif
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
c kk ---
ckbn       if(nodezero) then
ckbn         write(6,*)'tce_mrcc_uhf_step4'
ckbn         call util_flush(6)
ckbn       end if
c ------
      call ga_copy(g_a_coul, g_b_coul)
c kk ---
ckbn       if(nodezero) then
ckbn         write(6,*)'tce_mrcc_uhf_step5'
ckbn         call util_flush(6)
ckbn       end if
c ------
      call ga_dadd(one, g_a_coul, mone, g_a_exch, g_a_coul)
      call ga_dadd(one, g_b_coul, mone, g_b_exch, g_b_coul)
      if(cphf_uhf)then
         call ga_dadd(one, g_a_coul, one, g_a_xc, g_a_coul)
         call ga_dadd(one, g_b_coul, one, g_b_xc, g_b_coul)
      endif
c
c     reuse g_a_exch to hold the 1-e integrals
c
      g_hcore = g_a_exch
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel)
ckbn zora This is a hack to add zora in mrcc it assumes g_hcore is not deleted after dft calculation for core level paper
      if(logzorafromdft) then
       call ga_dadd(1.d0,g_hcore,1.d0,g_zora_Kinetic(1),g_hcore) ! zora kinetic
      endif
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
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf6", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif
ckbn      cpu = - util_cpusec()
ckbn      wall = - util_wallsec()
c
cjb      if (odebug) then
cjb         call ga_print(g_a_coul)
cjb         call ga_print(g_b_coul)
cjb      endif
c
c     Transform the Fock matrices to the MO basis using g_a_dens for scratch
c
c
c     Free up dead global arrays
c
      if (.not. ga_destroy(g_a_dens)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_dens)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      if(cphf_uhf)then
         if (.not. ga_destroy(g_a_xc)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_b_xc)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      endif
      if (.not. ga_destroy(g_a_exch)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_exch)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_a_coul)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_b_coul)) call errquit('uhf_e: destroy',0,
     &       GA_ERR)

c
c     extract the gradient
c
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('uhf_energy: no repulsion energy?', 0, GEOM_ERR)
      energy = eone + etwo + enrep
c
      if (ga_nodeid().eq.0) then
c         write(6,*) ' eone, etwo, enrep, energy ',
c     $        eone, etwo, enrep, energy
       write(6,'(A,4F20.12)')"eone,etwo,enrep,energy",
     +                       eone,etwo,enrep,energy
      endif
c
ckbn      cpu = cpu + util_cpusec()
ckbn      wall = wall + util_wallsec()
ckbn      if(nodezero) then
ckbn        write(LuOut,'(A,2f15.5)')"tce_mrcc_uhf7", cpu, wall
ckbn        call util_flush(LuOut)
ckbn      endif

      end
c $Id$
