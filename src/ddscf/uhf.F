      logical function uhf(rtdb, energy)
C$Id: uhf.F,v 1.13 1995-11-13 06:41:21 d3g681 Exp $
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "geom.fh"
#include "bas.fh"
c
c     Compute UHF wavefunction.
c
c     Initial vectors are assumed to be loaded in common/cscf/
c
      integer rtdb
      double precision energy   ! [output] Total SCF energy
c
c     local variables
c
      logical oprint, oprint_sym
      logical converged
      integer nprint
      double precision eone     ! One electron energy
      double precision etwo     ! Two electron energy
      double precision enrep    ! Effective nuclear repulsion energy
      double precision uhf_time
      double precision sz, s2
C
      Integer NAtoms, NShells
c
      logical uhf_nr_solve
      external ga_create_atom_blocked, uhf_nr_solve
c
      oprint = util_print('information',print_low)
c
c     Initialize the UHF data structures and allocate memory
c
      call uhf_init(rtdb)
c
c     Start timer
c
      call ga_sync()
      uhf_time = -util_cpusec()
c
c     Newton-Raphson based SCF
c
      converged = uhf_nr_solve(rtdb, energy,  eone, etwo, enrep)
c
      call ga_sync()
      uhf_time = uhf_time + util_cpusec()
c
c     Print out results
c
      call uhf_spin(sz,s2)
c
      if (ga_nodeid().eq.0 .and. (oprint .or. .not.converged)) then
         if (.not. converged) then
            write(6,*)
            call util_print_centered(6,
     $           'Calculation failed to converge', 20, .true.)
            write(6,*)
         end if
         write(6,2) scftype, energy, eone, etwo, enrep, 
     $        sz, s2, s2-sz*(sz+1), uhf_time
 2       format(//
     $        '       Final ',a4,' results '/
     $        '       ------------------ '//
     $        '         Total SCF energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '      Two electron energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12//
     $        '                       Sz =', f12.4/
     $        '                      S^2 =', f12.4/
     $        '            Contamination =', f12.4//
     $        '        Time for solution =', f9.1,'s'//)
         call util_flush(6)
      end if
c
      if (ga_nodeid() .eq. 0) then
C
         if (util_print('schwarz',print_high)) then
            if ( .NOT. geom_ncent(geom, natoms) ) Call ErrQuit(
     $         'rohf: problem with call to geom_ncent', geom)
            if ( .NOT. bas_numcont(basis, nshells) ) Call ErrQuit(
     $         'rohf: problem with call to bas_numcont', basis)
            call schwarz_print(natoms, nshells)
         Endif
C
         if (util_print('final evals', print_default)) then
            call util_print_centered(6,'Final alpha eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval),
     $           1, min(nalpha+10,nbf),
     $           1, 1, nbf, 1, 1)
            write(6,*)
            call util_print_centered(6,'Final beta eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval+nbf),
     $           1, min(nalpha+10,nbf),
     $           1, 1, nbf, 1, 1)
            call util_flush(6)
         end if
         if (util_print('final vectors', print_debug)) then
            write(6,*)
            call util_print_centered(6,'Final MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         end if
      end if
      if (oadapt) then
         nprint = min(nalpha+10,nbf)
         if (util_print('all vector symmetries', print_high))
     $        nprint = nbf
         oprint_sym = util_print('final vector symmetries',
     $        print_default)
         call scf_movecs_sym_adapt(basis, g_movecs, oprint_sym,
     $        nprint, '- alpha')
         call scf_movecs_sym_adapt(basis, g_movecs(2), oprint_sym,
     $        nprint, '- beta')
      endif
      if (util_print('final vectors', print_debug)) then
         call ga_print(g_movecs)
         call ga_print(g_movecs(2))
      end if
c
c     Clean up common blocks and allocated data and return
c
      call uhf_tidy(rtdb)
c
      uhf = converged
c
      end
      subroutine uhf_init(rtdb)
      implicit none
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb
c
      cuhf_vlen = nalpha*(nbf-nalpha) + nbeta*(nbf-nbeta)
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_init: Falpha',
     $     1, 1, cuhf_g_falpha)) call errquit('uhf_init: Fa', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_init: Fbeta',
     $     1, 1, cuhf_g_fbeta)) call errquit('uhf_init: Fb', 0)
c
      cuhf_init_flag = .true.
c
      end
      subroutine uhf_tidy(rtdb)
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      if (.not. ga_destroy(cuhf_g_falpha)) call errquit
     $     ('cuhf_tidy: failed to destroy fock matrix',0)
      if (.not. ga_destroy(cuhf_g_fbeta)) call errquit
     $     ('cuhf_tidy: failed to destroy fock matrix',0)
      cuhf_init_flag = .false.
c
      end
      logical function uhf_nr_solve(rtdb, energy, eone, etwo, enrep)
C     $Id: uhf.F,v 1.13 1995-11-13 06:41:21 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "cuhf.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "util.fh"
#include "cscf.fh"
c     
c     Solve the UHF equations using a hybrid NR/PCG method
c     
c     Arguments
c     
      integer rtdb
      double precision energy   ! Return UHF energy
      double precision eone, etwo, enrep ! Return UHF energy contribs
c     
c     Local GA handles 
c     
      integer g_grad            ! gradient vector  ... cuhf_vlen
      integer g_search          ! search direction ... cuhf_vlen
      integer g_work            ! scratch vector   ... cuhf_vlen
c     
c     Local variables
c     
      integer ls_max            ! Max. no. of points for line search
      parameter (ls_max=10)
      integer iwork(10)         ! PCG routine internal info
      double precision dwork(3*ls_max+5) ! PCG routine internal info
      character*16 mode         ! For communication with PCG routine 
c
      double precision step     ! Step to take in search direction
c
      double precision gnorm    ! measures for convergence
      logical converged
      logical oprint_parm, oprint_conv, oprint_eval, oprint_vecs
      logical oprint_ls
c     
c     Functions
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
C      external uhf_hessv, uhf_hessv_precond
c     
      if (ga_nodeid().eq.0) call util_flush(6)
c     
c     Allocate global arrays.  Gradient, search direction
c     and work space for the PCG routine.
c     
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: grad',
     $     1, 1, g_grad)) call errquit('uhf_nr_solve: ga_cre grad', 0)
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: serch',
     $     1, 1, g_search)) call errquit('uhf_nr_solve: ga_cre grad',0)
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: work',
     $     1, 1, g_work)) call errquit('uhf_nr_solve: ga_cre work', 0)
c     
      oprint_parm = util_print('parameters', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_ls   = util_print('line search', print_high)
      oprint_vecs = util_print('intermediate vectors', print_debug)
      oprint_eval = util_print('intermediate evals', print_debug)
c     
      if (ga_nodeid().eq.0.and. oprint_parm) then
         write(6,1) gnorm_tol, maxiter, tol2e
         if (ouser_changed_conv) write(6,11) shifts, nr_gswitch
         write(6,111)
 1       format(//,
     $        1x,'----------------------------------------------',/
     $        1x,'    Not yet quadratically convergent UHF',//,
     $        1x,'Convergence threshold     :',9x,1p,e10.3,0p,/,
     $        1x,'Maximum no. of iterations :',9x,i4,/,
     $        1x,'Integral*density screening:',9x,1p,e10.3,0p)
 11      format(/,
     $        1x,'PCG initial level shift   :',9x,f10.3,/,
     $        1x,'PCG change shift at maxg  :',9x,f10.3,/,
     $        1x,'PCG final level shift     :',9x,f10.3,/,
     $        1x,'NR  initial level shift   :',9x,f10.3,/,
     $        1x,'NR  change shift at maxg  :',9x,f10.3,/,
     $        1x,'NR  final level shift     :',9x,f10.3,/,
     $        1x,'NR  enabled at maxg       :',9x,f10.3)
 111     format(
     $        1x,'----------------------------------------------',/)
         call util_flush(6)
      end if
c     
c     PCG iterative loop
c     
      mode = 'startup'
      converged = .false.
 10   if (.not. converged) then
         call ga_pcg_minim(cuhf_vlen, iter, energy, g_grad, g_work,
     $        g_search, step, ls_tol, ls_max, eprec, oconjugacy, 
     $        oprint_conv, oprint_ls, iwork, dwork, mode)
c     
         if (mode .eq. 'energy+gradient') then
c     
c     Compute the energy and gradient at step*search direction
c     
            call uhf_step_energy(rtdb, step, g_search,
     $           eone, etwo, enrep, energy, g_grad)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            goto 10
c     
         else if (mode .eq. 'precondition') then
c     
c     Precondition the gradient direction with approx./exact hessian
c     
            call uhf_search_precond(rtdb, g_grad, g_work)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            goto 10
c     
         else if (mode .eq. 'accept step') then
c     
c     Apply the rotation and check for overall convergence.
c     The current fock matrices and energy correspond to the
c     updated MO coefficients.
c     
            call uhf_rotate_movecs(step, g_search)
            converged = gnorm .lt. gnorm_tol
            if ((iter.le.maxiter).and.(.not. converged)) then
c     
c     Not finished ... dump mo-vectors to disk for restart and
c     go back to top of loop
c     
               if (oprint_vecs) then
                  if (ga_nodeid() .eq. 0) then
                     write(6,*)
                     call util_print_centered(6,
     $                    'Intermediate MO vectors',40,.true.)
                     write(6,*)
                     call util_flush(6)
                  end if
                  call ga_print(g_movecs)
                  call ga_print(g_movecs(2))
               end if
               call scf_movecs_write
               goto 10
            end if
c            
         end if
      end if
c     
c     End SCF minimisation
c     
      if (ga_nodeid().eq.0) call util_flush(6)
c     
c     *ALWAYS* return canonical MOs ... do not do aufbau
c     
      call uhf_canon(.false., oprint_eval)
      call movecs_fix_phase(g_movecs)
      call movecs_fix_phase(g_movecs(2))
      call scf_movecs_write
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('uhf_nr_solve: ga_destroy grad', 0)
      if (.not. ga_destroy(g_search)) call errquit
     $     ('uhf_nr_solve: ga_destroy search', 0)
      if (.not. ga_destroy(g_work)) call errquit
     $     ('uhf_nr_solve: ga_destroy work', 0)
c     
      if (ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c     
      uhf_nr_solve = converged
c     
      end
      subroutine uhf_step_energy(rtdb, step, g_search,
     $     eone, etwo, enrep, energy, g_grad)
      implicit none
#include "cuhf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
c
      integer rtdb
      double precision step
      integer g_search
      double precision eone, etwo, enrep, energy
      integer g_grad
c
c     Compute the energy and gradient at the given step
c
      integer g_tmp(2)
      if (oscfps) call pstat_on(ps_search)
c
c     Apply rotation
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_s_e: alpha MOs',
     $     1, 1, g_tmp(1))) call errquit('uhf_s_e: MOs', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_s_e: beta MOs',
     $     1, 1, g_tmp(2))) call errquit('uhf_s_e: MOs', 0)
      call uhf_rotate(step, g_search, g_movecs, g_tmp)
c
c     Compute energy and gradient
c     
      call uhf_energy(g_tmp, eone, etwo, enrep, energy, g_grad)
c     
c     Tidy up
c
      if (.not. ga_destroy(g_tmp(1))) call errquit
     $     ('uhf_nr_solve: tmp?', 0)
      if (.not. ga_destroy(g_tmp(2))) call errquit
     $     ('uhf_nr_solve: tmp?', 0)
c
      if (oscfps) call pstat_off(ps_search)
c
      end
      subroutine uhf_rotate_movecs(step, g_search)
      implicit none
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
      double precision step
      integer g_search
c
      integer g_tmp(2)
c
c     Apply rotation
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_s_e: alpha MOs',
     $     1, 1, g_tmp(1))) call errquit('uhf_s_e: MOs', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_s_e: beta MOs',
     $     1, 1, g_tmp(2))) call errquit('uhf_s_e: MOs', 0)
c
      call uhf_rotate(step, g_search, g_movecs, g_tmp)
      call ga_copy(g_tmp(1), g_movecs(1))
      call ga_copy(g_tmp(2), g_movecs(2))
c
      if (.not. ga_destroy(g_tmp(1))) call errquit('urm: destroy?',0)
      if (.not. ga_destroy(g_tmp(2))) call errquit('urm: destroy?',0)
c
      end
      subroutine uhf_rotate(step, g_k, g_mos_in, g_mos_out)
      implicit none
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
      double precision step
      integer g_k, g_mos_in(2), g_mos_out(2)
c
c     Generate MOs * Exp(-step*K) in temporary arrays
c     
      double precision scale
c
      if (step .ne. 0.0d0) then
         call ga_dscal(g_k, step)
         call uhf_k2cf(basis,
     $        nbf, nalpha, nbeta,
     $        g_k, g_mos_in, g_mos_out)
         scale = 1.0d0 / step
         call ga_dscal(g_k, scale)
      else
         call ga_copy(g_mos_in(1), g_mos_out(1))
         call ga_copy(g_mos_in(2), g_mos_out(2))
      end if
c
      end
      subroutine uhf_k2cf(basis, nbf, nalpha, nbeta, g_kvec,
     $     g_mocf_in, g_mocf_out)
C$Id: uhf.F,v 1.13 1995-11-13 06:41:21 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
c     given the step in g_kvec generate the exponential rotation
c     and apply it to the input vectors
c
c     currently the step is interpreted as being the alpha rotations
c     followed by the beta rotations ... this will change when the
c     singlet/triplet rotations are adopted.
c     
      integer nbf, basis, nalpha, nbeta
      integer g_kvec
      integer g_mocf_out(2), g_mocf_in(2)
c     
      integer nocc, nvir, nocc1, ioff
      integer g_tmp
      integer i, l_tmp, k_tmp, iset
      double precision mone, one, zero
      parameter(mone=-1.0d0, one=1.0d0, zero=0.0d0)
c     
      call ga_sync()
c     
      if (.not. ma_push_get(mt_dbl,nbf,'temp k',l_tmp,k_tmp))
     $     call errquit('uhf_k2cf: ma failed on tmp',  nbf)
c
      do iset = 1, 2            ! 1= alpha, 2= beta
c
         if (iset .eq. 1) then
            nocc = nalpha
            ioff = 0
         else
            nocc = nbeta
            ioff = nalpha*(nbf - nalpha)
         end if
         nvir  = nbf - nocc
         nocc1 = nocc + 1
c     
c     form k-matrix from vector
c
         call ga_zero(g_mocf_out(iset))
         do i=ga_nodeid()+1,nocc,ga_nnodes()
            call ga_get(g_kvec, ioff+1+(i-1)*nvir, ioff+i*nvir, 1, 1,
     $           dbl_mb(k_tmp), 1)
            call ga_put(g_mocf_out(iset), nocc1, nbf, i, i,
     $           dbl_mb(k_tmp), 1)
            call dscal(nvir, mone, dbl_mb(k_tmp), 1)
            call ga_put(g_mocf_out(iset), i, i, nocc1, nbf,
     $           dbl_mb(k_tmp), 1)
         end do
c
c     Make near zeoes exactly zero so that sparsity tests in the
c     matrix multiply work well
c
         call ga_screen(g_mocf_out(iset), 1d-12)
c
c     exponentiate the matrix
c
         call matrix_exp(g_mocf_out(iset))
c
c     apply the rotation
c
         if (.not. ga_create(mt_dbl, nbf, nbf, 'uhf_k2cf:tmp',
     $        1, 1, g_tmp)) call errquit('uhf_k2cf:ga_create?', nbf)
c
         call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_mocf_in(iset),
     $        g_mocf_out(iset),
     $        zero, g_tmp)
         call ga_copy(g_tmp, g_mocf_out(iset))
c
         if (.not. ga_destroy(g_tmp))
     $        call errquit('uhf_k2cf: ga_destroy?', 0)
c     
*         call ga_orthog_mos(basis, g_mocf_out(iset))
      end do
c
      if (.not. ma_pop_stack(l_tmp))
     $     call errquit('uhf_k2cf: pop failed', 0)
c     
      end
      subroutine uhf_energy( g_vecs, eone, etwo, enrep, energy,
     $     g_grad )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "util.fh"
c     
c     $Id: uhf.F,v 1.13 1995-11-13 06:41:21 d3g681 Exp $
c
      integer g_vecs(2)
      double precision energy
      integer g_grad
c     
      double precision eone, etwo, enrep
      integer gtype, grow, gcol
      integer d(4), f(4)
      integer g_a_dens, g_a_coul, g_a_exch
      integer g_b_dens, g_b_coul, g_b_exch
      integer g_hcore
      double precision jfac(4), kfac(4), one, zero, mone
      parameter (one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision e_a_coul, e_a_exch, e_b_coul, e_b_exch
      double precision errmaxa, errmaxb
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical odebug
      data jfac / 1.0d0, 0.0d0, 1.0d0, 0.0d0/
      data kfac / 0.0d0, 1.0d0, 0.0d0, 1.0d0/
*      data jfac / 2.0d0, 0.0d0, 2.0d0, 0.0d0/
*      data kfac / 0.0d0, 2.0d0, 0.0d0, 2.0d0/
c     
c     Check
c     
      odebug = util_print('uhf_debug', print_debug)
      if (.not.cuhf_init_flag)
     $     call errquit('uhf_energy: UHF internal block invalid',0)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.cuhf_vlen).or.(gcol.ne.1))
     $     call errquit('uhf_energy: invalid vector length',0)
c     
c     Arrays for AO density, coulomb and exchange matrices
c
      g_a_coul = ga_create_atom_blocked(geom, basis, 'uhf:a coul')
      g_b_coul = ga_create_atom_blocked(geom, basis, 'uhf:b coul')
      g_a_exch = ga_create_atom_blocked(geom, basis, 'uhf:a exch')
      g_b_exch = ga_create_atom_blocked(geom, basis, 'uhf:b exch')
      g_a_dens = ga_create_atom_blocked(geom, basis, 'uhf:a dens')
      g_b_dens = ga_create_atom_blocked(geom, basis, 'uhf:b dens')
      call ga_zero(g_a_dens)
      call ga_zero(g_b_dens)
c
c     Make the densites and build the fock matrices
c
      call ga_dgemm('n', 't', nbf, nbf, nalpha, one, g_vecs(1),
     $     g_vecs(1), zero, g_a_dens)
      call ga_dgemm('n', 't', nbf, nbf, nbeta, one, g_vecs(2),
     $     g_vecs(2), zero, g_b_dens)
c
c     Since UHF can break spatial symmetry by localizing the orbitals
c     the densities may not be totally symmetric, but since the Hamiltonian
c     is symmetric contraction with the integrals projects out the totally 
c     symmetric component ... hence we can symmetrize the densities and
c     exploit symmetry.  Compute the max change in any element due to
c     symmetrizing and print a warning if it is big.
c
c     If this is the case then where does the 'force' for symmetry breaking
c     come from?  Must be missing something?
c
      call ga_copy(g_a_dens,g_a_coul)
      call ga_copy(g_b_dens,g_b_coul)
      call sym_symmetrize(geom, basis, .true., g_a_dens)
      call sym_symmetrize(geom, basis, .true., g_b_dens)
      call ga_dadd(one, g_a_dens, mone, g_a_coul, g_a_coul)
      call ga_dadd(one, g_b_dens, mone, g_b_coul, g_b_coul)
      call ga_maxelt(g_a_coul, errmaxa)
      call ga_maxelt(g_b_coul, errmaxb)
      if (max(errmaxa,errmaxb).gt.1d-10) then
         if (ga_nodeid().eq.0) then
            write(6,77) errmaxa,errmaxb
 77         format(' Warning: spatial symmetry breaking in UHF: ',
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
c
      call ga_zero(g_a_coul)
      call ga_zero(g_b_coul)
      call ga_zero(g_a_exch)
      call ga_zero(g_b_exch)
      d(1) = g_a_dens
      d(2) = g_a_dens
      d(3) = g_b_dens
      d(4) = g_b_dens
      f(1) = g_a_coul
      f(2) = g_a_exch
      f(3) = g_b_coul
      f(4) = g_b_exch
c
      call do_ri (.false.)
      call fock_2e(geom, basis, 4, jfac, kfac, tol2e,
     $     oskel, d, f)
      call do_ri (.true.)
c
      e_a_coul = 0.5d0*
     $     (ga_ddot(g_a_dens,g_a_coul) + ga_ddot(g_a_dens,g_b_coul))
      e_b_coul = 0.5d0*
     $     (ga_ddot(g_b_dens,g_a_coul) + ga_ddot(g_b_dens,g_b_coul))
      e_a_exch = 0.5d0*ga_ddot(g_a_dens,g_a_exch)
      e_b_exch = 0.5d0*ga_ddot(g_b_dens,g_b_exch)
      etwo = e_a_coul + e_b_coul - e_a_exch - e_b_exch
c
      if (odebug .and. ga_nodeid().eq.0) then
         write(6,*) ' coulomb energies', e_a_coul, e_b_coul
         write(6,*) ' exchang energies', e_a_exch, e_b_exch
         call util_flush(6)
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
c
c     reuse g_a_exch to hold the 1-e integrals
c
      g_hcore = g_a_exch
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_hcore, 'potential', oskel)
      eone = 
     $     (ga_ddot(g_a_dens,g_hcore) + ga_ddot(g_b_dens,g_hcore))
      call ga_dadd(one, g_hcore, one, g_a_coul, g_a_coul)
      call ga_dadd(one, g_hcore, one, g_b_coul, g_b_coul)
c
      call sym_symmetrize(geom, basis, .false., g_a_coul)
      call sym_symmetrize(geom, basis, .false., g_b_coul)
c
      if (odebug) then
         call ga_print(g_a_coul)
         call ga_print(g_b_coul)
      endif
c
c     Transform the Fock matrices to the MO basis using g_a_dens for scratch
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
      if (.not. ga_destroy(g_a_dens)) call errquit('uhf_e: destroy',0)
      if (.not. ga_destroy(g_b_dens)) call errquit('uhf_e: destroy',0)
      if (.not. ga_destroy(g_a_exch)) call errquit('uhf_e: destroy',0)
      if (.not. ga_destroy(g_b_exch)) call errquit('uhf_e: destroy',0)
      if (.not. ga_destroy(g_a_coul)) call errquit('uhf_e: destroy',0)
      if (.not. ga_destroy(g_b_coul)) call errquit('uhf_e: destroy',0)
c
c     extract the gradient
c
      call uhf_get_grad(g_grad)
c
      if (odebug) call ga_print(g_grad)
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('uhf_energy: no repulsion energy?', 0)
      energy = eone + etwo + enrep
c
      if (odebug .and. ga_nodeid().eq.0) then
         write(6,*) ' eone, etwo, enrep, energy ',
     $        eone, etwo, enrep, energy
      endif
c
      end
      subroutine uhf_get_grad(g_grad)
      implicit none
#include "cscf.fh"
#include "cuhf.fh"
      integer g_grad
c
      integer ioff
c
      call ga_copy_patch('n',
     $     cuhf_g_falpha, nalpha+1, nbf, 1, nalpha,
     $     g_grad, 1, (nbf-nalpha)*nalpha, 1, 1)
c
      ioff = (nbf-nalpha)*nalpha
      call ga_copy_patch('n',
     $     cuhf_g_fbeta, nbeta+1, nbf, 1, nbeta,
     $     g_grad, ioff+1, ioff+(nbf-nbeta)*nbeta, 1, 1)
c
      call ga_dscal(g_grad, 2.0d0)
c
      end
      subroutine uhf_search_precond(rtdb, g_grad, g_work)
      implicit none
#include "cuhf.fh"      
#include "cscfps.fh"
#include "global.fh"
#include "cscf.fh"
#include "util.fh"
      integer rtdb
      integer g_grad
      integer g_work
c
c     Precondition the gradient with an approximation to the Hessian
c     return the result in g_work.
c
c     Right now only divide by the level-shifted diagonals which 
c     is in uhf the preconditioner to the preconditioner!
c
      logical oprint_conv, oprint_eval, oprint_debug
      double precision gnorm
      double precision max_grad
      double precision min_shift
c
      oprint_conv = util_print('convergence', print_default)
      oprint_eval = util_print('intermediate evals', print_debug)
      oprint_debug= util_print('uhf_debug', print_debug)
c
      call ga_maxelt(g_grad, max_grad)
      gnorm = sqrt(ga_ddot(g_grad, g_grad))
      eprec = max(gnorm*0.01d0, tol2e*100.0d0) ! Approx precision available
c
c     Determine shift necessary to make the preconditioner
c     positive definite
c
      call uhf_make_shift(min_shift)
c
c     Figure out which preconditioner to use, convergence
c     thresholds, shifts, ...
c
      call scf_precond_select(gnorm, max_grad, min_shift,
     $     oprint_conv)
c
c     Do a diagonalization if it will make the Hessian more
c     diagonally dominant
c
      if (ododiag) then
         call uhf_canon(.false., oprint_eval)
         call uhf_get_grad(g_grad)
         call uhf_make_shift(min_shift) ! Recompute
      end if
c
c     After doing the diag make sure that the shift is still OK
c
      if (lshift .lt. min_shift) then
         lshift = min_shift + 2.0d0
         if (ga_nodeid().eq.0 .and. oprint_conv) then
            write(6,3131) lshift
 3131       format('  Setting level-shift to ', f6.2,
     $           ' to force positive preconditioner')
            call util_flush(6)
         end if
      end if
c
      call ga_copy(g_grad,g_work)
      call uhf_precond(g_work, 0.0d0)
c
      end
      subroutine uhf_precond(g_x, solveshift)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
#include "cuhf.fh"
      integer g_x
      double precision solveshift
c      
      integer nvir, nocc(2), iset
      integer i, j, ioff, ibase
      integer l_diag, k_diag, l_x, k_x
      integer g_f(2)
      double precision diag, denominator, shift
      diag(i) = dbl_mb(k_diag + i - 1)
c     
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta
c
      if (.not.ma_push_get(MT_DBL,nbf,'uhf: tmp',l_diag,k_diag))
     $     call errquit('uhf_precond: cannot allocate',0)
      if (.not.ma_push_get(MT_DBL,nbf,'uhf: tmp',l_x,k_x))
     $     call errquit('uhf_hdiag_scale: cannot allocate',0)
c
      shift = lshift - solveshift
c
      ibase = 1
      do iset = 1, 2
         nvir = nbf - nocc(iset)
         call ga_get_diagonal(g_f(iset), dbl_mb(k_diag))
c
         do i = ga_nodeid()+1, nocc(iset), ga_nnodes()
            ioff = (i-1)*nvir + ibase
            call ga_get(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
            do j=1,nvir
               denominator = 2.d0*(diag(j+nocc(iset)) - diag(i)) + shift
               if (denominator .lt. 0.1d0) denominator = 0.1d0
               dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
            end do
            call ga_put(g_x, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_x),1)
         end do
         ibase = ibase + nocc(1)*(nbf-nocc(1))
      end do
c
      if (.not. ma_pop_stack(l_x)) call errquit('uhf:pop x',0)
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf:pop ',0)
c
      end
      subroutine uhf_make_shift(shift)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "cuhf.fh"
#include "cscf.fh"
c
c     Return minimum shift necessary for diagonal approximation
c     to 1-e Hessian to be positive definite
c
      double precision shift
c
      integer l_diag, k_diag
      double precision shift_a, shift_b
c
      if (.not. ma_push_get(mt_dbl, nbf, 'uhf_m_s: diag', 
     $     l_diag, k_diag)) call errquit
     $     ('uhf_make_shift: insufficient memory', nbf)
c
      call ga_sync()
c
      shift_a = 0.0d0
      shift_b = 0.0d0
c
      call rohf_make_shift_2(cuhf_g_falpha, dbl_mb(k_diag),
     $     1, nalpha, nalpha+1, nbf,
     $     2.0d0, shift_a)
      call rohf_make_shift_2(cuhf_g_fbeta, dbl_mb(k_diag),
     $     1, nbeta, nbeta+1, nbf,
     $     2.0d0, shift_b)
c
      shift = max(shift_a, shift_b)
c
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf_n_s: e',0)
c
      end
      subroutine uhf_canon(oaufbau, oprint)
C$Id: uhf.F,v 1.13 1995-11-13 06:41:21 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "cscfps.fh"
#include "cscf.fh"
#include "cuhf.fh"
      logical oaufbau
      logical oprint
c
      integer nocc(2), ioff, iset
      integer g_u, g_fock, g_f(2), g_tmp
      double precision one, zero
      data one, zero/1.d0, 0.d0/
c
c     This routine assumes that uhf_energy/uhf_fock have been called
c     so that the contents of /cuhf/ are current.
c
c     Diagonalize the UHF 'Fock' matrices
c     
c     If (oaufbau) 
c        diagonalize the whole thing and allow mixing of closed-open-virt
c     else
c        diagonalize separately the occ-occ, and virt-virt parts
c
c     Transform Fock matrices and MO coefficients into the new canonical basis
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_canon: u',
     $     1, 1, g_u)) call errquit('uhf_canon: ga failed for u', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_canon: fock',
     $     1, 1, g_fock)) call errquit
     $     ('uhf_canon: ga failed for fock', 0)
c
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta
c
      ioff = 0
      do iset = 1, 2
c
         if (oscfps) call pstat_on(ps_diag)
c
         call ga_zero(g_u)
         call ga_copy(g_f(iset), g_fock)
c     
         if (oaufbau) then
#if defined(PARALLEL_DIAG)
            call ga_diag_std(g_fock, g_u, dbl_mb(k_eval+ioff))
#else 
            call ga_diag_std_seq(g_fock, g_u, dbl_mb(k_eval+ioff))
#endif
         else
c     
c     occ-occ piece
c
            call rohf_canon_subspace(g_fock, g_u,
     $           dbl_mb(k_eval+ioff),
     $           1, nocc(iset))
c
c     virt-virt piece
c     
            call rohf_canon_subspace(g_fock, g_u,
     $           dbl_mb(k_eval+ioff+nocc(iset)),
     $           nocc(iset)+1, nbf)

         end if
c
         ioff = ioff + nbf
c
         if (oscfps) call pstat_off(ps_diag)
c
         call movecs_fix_phase(g_u)
c
c     Apply rotation to orbitals and fock matrix
c
         if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_canon: tmp',
     $        1, 1, g_tmp)) call errquit
     $        ('uhf_canon: ga failed for tmp', 0)
c
         call ga_copy(g_movecs(iset), g_tmp)
         call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_tmp, g_u,
     $        zero, g_movecs(iset))
c
         call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_f(iset), g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nbf, nbf, nbf, one, g_u, g_tmp, 
     $        zero, g_f(iset))
c
         if (.not. ga_destroy(g_tmp))
     $        call errquit('uhf_canon: destroy', 0)
      end do
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6, 'Alpha-spin eigenvalues',
     $        20, .true.)
         call output(dbl_mb(k_eval), 1, min(nalpha+5,nbf), 
     $        1, 1, nbf, 1, 1)
         write(6,*)
         call util_print_centered(6, 'Beta-spin eigenvalues',
     $        20, .true.)
         call output(dbl_mb(k_eval+nbf), 1, min(nalpha+5,nbf), 
     $        1, 1, nbf, 1, 1)
         call util_flush(6)
      end if
c
      if (.not. ga_destroy(g_u))
     $     call errquit('uhf_canon: destroy', 0)
      if (.not. ga_destroy(g_fock))
     $     call errquit('uhf_canon: destroy', 0)
c
      end
      subroutine uhf_spin(sz, s2)
      implicit none
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      double precision sz, s2
c
c     Return the expection values of Sz and S^2
c
      integer g_over, g_work, l_x, k_x, i, j
c
      sz = 0.5d0 * (nalpha - nbeta)
c     
c     Form Sab
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_spin: overlap',
     $     1, 1, g_over)) call errquit('uhf_spin: overlap', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_spin: work',
     $     1, 1, g_work)) call errquit('uhf_spin: work', 0)
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      call two_index_transf(g_over, g_movecs(1), g_movecs(2),
     $     g_work, g_over)
      if (.not. ga_destroy(g_work)) call errquit('uhf_spin:ga?',0)
c
c     Sum overlap of occupied orbitals
c
      if (.not. ma_push_get(mt_dbl,nbf,'uhf_spin:x',l_x,k_x))
     $     call errquit('uhf_spin: ma failed on x',  nalpha)
c
      s2 = 0
      do i = 1+ga_nodeid(), nbeta, ga_nnodes()
         call ga_get(g_over, 1, nalpha, i, i, dbl_mb(k_x), 1)
         do j = 1, nalpha
            s2 = s2 + dbl_mb(k_x+j-1)**2
         enddo
      enddo
      call ga_dgop(msg_uhf, s2, 1, '+')
      s2 = sz*(sz+1) + nbeta - s2
c
      if (.not. ma_pop_stack(l_x)) call errquit('uhf_spin:ma?',0)
      if (.not. ga_destroy(g_over)) call errquit('uhf_spin:ga?',0)
c
      end
