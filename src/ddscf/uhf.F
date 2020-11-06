      logical function uhf(rtdb, energy)
C$Id$
      implicit none
#include "errquit.fh"
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
c     Compute UHF wavefunction .
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
      double precision ecosmo   ! COSMO solvation energy
      double precision uhf_time
      double precision sz, s2
c
      integer ilo, ihi          ! For printing movecs analysis
      double precision eval_pr_tol_lo, eval_pr_tol_hi
      parameter (eval_pr_tol_lo = -1.5d0, eval_pr_tol_hi=0.5)
C
      Integer NAtoms, NShells
c
      logical uhf_nr_solve
      external uhf_nr_solve
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
      ecosmo = 0.0d0
      if (maxiter .le. 0) then
         converged = .false.
         energy = 0.0d0
         eone = 0.0d0
         etwo = 0.0d0
      else
         converged = uhf_nr_solve(rtdb, energy, eone, etwo, enrep,
     &                            ecosmo)
      endif
c
      call ga_sync()
      uhf_time = uhf_time + util_cpusec()
c
c     Print out results
c
      call uhf_spin(sz,s2)
      if (.not. rtdb_put(rtdb,'scf:sz', mt_dbl, 1, sz)) call errquit
     $     ('uhf: storing sz failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'scf:s2', mt_dbl, 1, s2)) call errquit
     $     ('uhf: storing s2 failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'uhf:coulomb', mt_dbl, 1, etwo)) call 
     $     errquit('uhf: writing etwo failed', 0, RTDB_ERR)
c
      if (ga_nodeid().eq.0 .and. (oprint .or. .not.converged)) then
         if (.not. converged) then
            write(6,*)
            call util_print_centered(6,
     $           'Calculation failed to converge', 20, .true.)
            write(6,*)
         end if
         write(6,2) scftype, energy, eone, etwo, enrep  
         if (abs(ecosmo).gt.0.0d0) then
           write(6,3) ecosmo
         endif
         write(6,4) sz, sz*(sz+1), s2, uhf_time
 2       format(//
     $        '       Final ',a4,' results '/
     $        '       ------------------ '//
     $        '         Total SCF energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '      Two electron energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12) 
 3       format(
     $        '             COSMO energy =', f20.12) 
 4       format(/
     $        '                       Sz =', f12.4/
     $        '                 Sz(Sz+1) =', f12.4/
     $        '                      S^2 =', f12.4//
     $        '        Time for solution =', f9.1,'s'//)
         call util_flush(6)
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print1('one-electron energy', mt_dbl, eone, 1)
         call ecce_print1('two-electron energy', mt_dbl, etwo, 1)
         call ecce_print1('nuclear repulsion energy', mt_dbl, enrep, 1)
         call ecce_print1('Sz', mt_dbl, sz, 1)
         call ecce_print1('S^2', mt_dbl, s2, 1)
      end if
C
      if (util_print('schwarz',print_high)) then
         if ( .NOT. geom_ncent(geom, natoms) ) Call ErrQuit(
     $        'uhf: problem with call to geom_ncent', geom, 0)
         if ( .NOT. bas_numcont(basis, nshells) ) Call ErrQuit(
     $        'uhf: problem with call to bas_numcont', basis,0)
         call schwarz_print(natoms, nshells)
      Endif
C
      nprint = min(nalpha+10,nmo)
      if (util_print('all vector symmetries', print_high))
     $     nprint = nmo
      if (oadapt) then
         oprint_sym = util_print('final vector symmetries',
     $        print_default)
         call scf_movecs_sym_adapt(basis, g_movecs, oprint_sym,
     $        nprint, '- alpha', .true., int_mb(k_irs))
         call scf_movecs_sym_adapt(basis, g_movecs(2), oprint_sym,
     $        nprint, '- beta', .true., int_mb(k_irs+nmo))
      endif
      if (ga_nodeid() .eq. 0) then
         if (util_print('final evals', print_default)) then
            call util_print_centered(6,'Final alpha eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval), 1, nprint, 1, 1, nmo, 1, 1)
            write(6,*)
            call util_print_centered(6,'Final beta eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval+nbf), 1, nprint, 1, 1, nmo, 1, 1)
            call util_flush(6)
         end if
      endif
      if (util_print('final vectors analysis', print_default)) then
         do ilo = 1,max(1,min(nalpha,nbeta)-10)
            if (dbl_mb(k_eval+ilo-1) .ge. eval_pr_tol_lo) goto 961
         enddo
 961     do ihi = min(max(nalpha,nbeta)+10,nmo),nmo
            if (dbl_mb(k_eval+ihi-1) .ge. eval_pr_tol_hi) goto 9611
         enddo
         ihi = max(ihi-1,1)
 9611    continue
         call movecs_print_anal(basis, ilo, ihi, 0.15d0, g_movecs, 
     $        'UHF Final Alpha Molecular Orbital Analysis', 
     $        .true., dbl_mb(k_eval), oadapt, int_mb(k_irs),
     $        .true., dbl_mb(k_occ))
         call movecs_print_anal(basis, ilo, ihi, 0.15d0, g_movecs(2), 
     $        'UHF Final Beta Molecular Orbital Analysis', 
     $        .true., dbl_mb(k_eval+nbf), oadapt, int_mb(k_irs+nmo),
     $        .true., dbl_mb(k_occ+nbf))
      endif
      if (ga_nodeid() .eq. 0) then
         if (util_print('final vectors', print_debug)) then
            write(6,*)
            call util_print_centered(6,'Final MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         end if
      end if
      if (util_print('final vectors', print_debug)) then
         call ga_print(g_movecs)
         call ga_print(g_movecs(2))
      end if
c
      call uhf_analyze(rtdb)
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
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb
c
      cuhf_vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
c
c     Special case for no degrees of freedom ... avoid breaking GA
c
      if (cuhf_vlen .eq. 0) cuhf_vlen = 1
c
      if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_init: Falpha',
     $     32, 32, cuhf_g_falpha)) call errquit('uhf_init: Fa', 0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_init: Fbeta',
     $     32, 32, cuhf_g_fbeta)) call errquit('uhf_init: Fb', 0,
     &       GA_ERR)
c
      cuhf_init_flag = .true.
      noskew_uhf = .true.
c
      end
      subroutine uhf_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb
c
      if (.not. ga_destroy(cuhf_g_falpha)) call errquit
     $     ('uhf_tidy: failed to destroy fock matrix',0, GA_ERR)
      if (.not. ga_destroy(cuhf_g_fbeta)) call errquit
     $     ('uhf_tidy: failed to destroy fock matrix',0, GA_ERR)
      cuhf_init_flag = .false.
c
      end
      logical function uhf_nr_solve(rtdb, energy, eone, etwo, enrep,
     $                              ecosmo)
C     $Id$
      implicit none
#include "errquit.fh"
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
      double precision ecosmo !< [Output] The COSMO solvation energy
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
      parameter (ls_max=20)
      integer iwork(10)         ! PCG routine internal info
      double precision dwork(3*ls_max+5) ! PCG routine internal info
      character*16 mode         ! For communication with PCG routine 
c
      double precision tlastwrite,ttest  ! time last wrote the MO vectors
      double precision step     ! Step to take in search direction
      double precision time_iter_start
      integer time_prev_iter
      logical noscf ! one-shot energy option
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
*ga:1:0
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: grad',
     $     0, 1, g_grad)) call errquit('uhf_nr_solve: ga_cre grad', 0,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: serch',
     $     0, 0, g_search)) call errquit('uhf_nr_solve: ga_cre grad',0,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, cuhf_vlen, 1, 'uhf_nr_solve: work',
     $     0, 0, g_work)) call errquit('uhf_nr_solve: ga_cre work', 0,
     &       GA_ERR)
c     
      oprint_parm = util_print('parameters', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_ls   = util_print('line search', print_high)
      oprint_vecs = util_print('intermediate vectors', print_debug)
      oprint_eval = util_print('intermediate evals', print_debug)
c
      tlastwrite = util_wallsec()
c     
      if (ga_nodeid().eq.0.and. oprint_parm) then
         write(6,1) gnorm_tol, maxiter, tol2e
         if (ouser_changed_conv) write(6,11) shifts, nr_gswitch
         write(6,111)
 1       format(//,
     $        1x,'----------------------------------------------',/
     $        1x,'        Quadratically convergent UHF',//,
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
      ododiag = .true.
      odisable_nr = .false.
c     
c     PCG iterative loop
c     
      time_prev_iter = 0
      time_iter_start = util_wallsec()
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
     $           eone, etwo, enrep, ecosmo, energy, g_grad)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            call ecce_print1('total energy', mt_dbl, energy, 1)
            call ecce_print1('orbital gradient norm', mt_dbl, gnorm, 1)
c
c           == is this a oneshot "noscf" type calculation ==
            noscf = .false.
            if (.not.rtdb_get(rtdb, 'scf:noscf', mt_log, 1, noscf))
     &                      noscf=.false.
            if (noscf) converged = .true.
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
            time_prev_iter = util_wallsec() - time_iter_start
            time_iter_start = util_wallsec()
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
               ttest = util_wallsec()-tlastwrite
               call ga_dgop(1324,ttest,1,'max')
               if (ttest .gt. 300) then
*     if (olock) call scf_movecs_lock
                 call scf_movecs_write(rtdb)
                 tlastwrite = util_cpusec()
               endif
               if (util_test_time_remaining(rtdb,time_prev_iter*3))
     $              goto 10
            end if
c            
         end if
      end if
c     
c     End SCF minimisation
c     
      if (ga_nodeid().eq.0) call util_flush(6)
c
      if (.not.noscf) then  ! bypass for noscf
c     
c     *ALWAYS* return canonical MOs ... do not do aufbau
c     
      call uhf_canon(.false., oprint_eval)
      if (oadapt) then
        call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     $        0, 'after canon', .true., int_mb(k_irs))
        call scf_movecs_sym_adapt(basis, g_movecs(2), .false.,
     $        0, 'after canon', .true., int_mb(k_irs+nmo))
      endif
      call movecs_fix_phase(g_movecs)
      call movecs_fix_phase(g_movecs(2))
      if (olock) call scf_movecs_lock
      call scf_movecs_write(rtdb)
c
      end if  ! bypass for noscf
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('uhf_nr_solve: ga_destroy grad', 0, GA_ERR)
      if (.not. ga_destroy(g_search)) call errquit
     $     ('uhf_nr_solve: ga_destroy search', 0, GA_ERR)
      if (.not. ga_destroy(g_work)) call errquit
     $     ('uhf_nr_solve: ga_destroy work', 0, GA_ERR)
c     
      if (ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c     
      uhf_nr_solve = converged
c     
      end
      subroutine uhf_step_energy(rtdb, step, g_search,
     $     eone, etwo, enrep, ecosmo, energy, g_grad)
      implicit none
#include "errquit.fh"
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
      double precision ecosmo !< [Output] The COSMO solvation energy
      integer g_grad
c
c     Compute the energy and gradient at the given step
c
      integer g_tmp(2)
      if (oscfps) call pstat_on(ps_search)
c
c     Apply rotation
c
      if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_s_e: alpha MOs',
     $     32, 32, g_tmp(1))) call errquit('uhf_s_e: MOs', 0, GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_s_e: beta MOs',
     $     32, 32, g_tmp(2))) call errquit('uhf_s_e: MOs', 0, GA_ERR)
c
      call uhf_rotate(step, g_search, g_movecs, g_tmp)
c
c     Compute energy and gradient
c     
      call uhf_energy(rtdb, g_tmp, eone, etwo, enrep, ecosmo, energy,
     $                g_grad)
c     
c     Tidy up
c
      if (.not. ga_destroy(g_tmp(1))) call errquit
     $     ('uhf_nr_solve: tmp?', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp(2))) call errquit
     $     ('uhf_nr_solve: tmp?', 0, GA_ERR)
c
      if (oscfps) call pstat_off(ps_search)
c
      end
      subroutine uhf_rotate_movecs(step, g_search)
      implicit none
#include "errquit.fh"
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
      if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_s_e: alpha MOs',
     $     32, 32, g_tmp(1))) call errquit('uhf_s_e: MOs', 0, GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_s_e: beta MOs',
     $     32, 32, g_tmp(2))) call errquit('uhf_s_e: MOs', 0, GA_ERR)
c
      call uhf_rotate(step, g_search, g_movecs, g_tmp)
      call ga_copy(g_tmp(1), g_movecs(1))
      call ga_copy(g_tmp(2), g_movecs(2))
c
      if (.not. ga_destroy(g_tmp(1))) call errquit('urm: destroy?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp(2))) call errquit('urm: destroy?',0,
     &       GA_ERR)
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
     $        nbf, nmo, nalpha, nbeta,
     $        g_k, g_mos_in, g_mos_out)
         scale = 1.0d0 / step
         call ga_dscal(g_k, scale)
      else
         call ga_copy(g_mos_in(1), g_mos_out(1))
         call ga_copy(g_mos_in(2), g_mos_out(2))
      end if
c
      end
      subroutine uhf_k2cf(basis, nbf, nmo, nalpha, nbeta, g_kvec,
     $     g_mocf_in, g_mocf_out)
C$Id$
      implicit none
#include "errquit.fh"
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
      integer nbf, nmo, basis, nalpha, nbeta
      integer g_kvec
      integer g_mocf_out(2), g_mocf_in(2)
c     
      integer nocc, nvir, nocc1, ioff
      integer g_kmat
      integer i, l_tmp, k_tmp, iset
      double precision mone, one, zero, maxk
      parameter(mone=-1.0d0, one=1.0d0, zero=0.0d0)
c     
      call ga_sync()
c     
      if (.not. ma_push_get(mt_dbl,nmo,'temp k',l_tmp,k_tmp))
     $     call errquit('uhf_k2cf: ma failed on tmp',  nmo, MA_ERR)
      if (.not. ga_create(mt_dbl, nmo, nmo, 'uhf_k2cf:kmat',
     $     32, 32, g_kmat)) call errquit('uhf_k2cf:ga_create?', nmo,
     &       MA_ERR)
c
      do iset = 1, 2            ! 1= alpha, 2= beta
c
         if (iset .eq. 1) then
            nocc = nalpha
            ioff = 0
         else
            nocc = nbeta
            ioff = nalpha*(nmo - nalpha)
         end if
         nvir  = nmo - nocc
         nocc1 = nocc + 1
c     
c     form k-matrix from vector
c
         call ga_zero(g_kmat)
         if (nvir*nocc .gt.0) then
            do i=ga_nodeid()+1,nocc,ga_nnodes()
               call ga_get(g_kvec, ioff+1+(i-1)*nvir, ioff+i*nvir, 1, 1,
     $              dbl_mb(k_tmp), nvir)
               call ga_put(g_kmat, nocc1, nmo, i, i,
     $              dbl_mb(k_tmp), nmo-nocc1+1)
               call dscal(nvir, mone, dbl_mb(k_tmp), 1)
               call ga_put(g_kmat, i, i, nocc1, nmo,
     $              dbl_mb(k_tmp), 1)
            end do
         endif
c
c     Make near zeoes exactly zero so that sparsity tests in the
c     matrix multiply work well.  Can't screen with too large
c     a value since this breaks symmetry in non-abelian groups.
c     Must also allow value to tend to zero to allow for tight
c     convergence.
c
         call ga_maxelt(g_kmat, maxk)
         call ga_screen(g_kmat, min(maxk*maxk,maxk*1d-2,1d-12))
c
c     Exponentiate the matrix and apply the rotation
c
         call matrix_exp(g_kmat)
         call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_mocf_in(iset),
     $        g_kmat, zero, g_mocf_out(iset))
c
*     call ga_orthog_mos(basis, g_mocf_out(iset)) ! No longer necessary
      end do
c
      if (.not. ga_destroy(g_kmat))
     $     call errquit('uhf_k2cf: ga_destroy?', 0, GA_ERR)
      if (.not. ma_pop_stack(l_tmp))
     $     call errquit('uhf_k2cf: pop failed', 0, MA_ERR)
c     
      end
      subroutine uhf_energy(rtdb, g_vecs, eone, etwo, enrep, ecosmo, 
     $                      energy, g_grad )
c
      implicit none
c
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
#include "zora.fh"
#include "frozemb.fh"
#include "cosmo.fh"
c     
c     $Id$
c
      integer rtdb
      integer g_vecs(2)
      double precision energy
      integer g_grad
c     
      double precision eone, etwo, enrep
      integer gtype, grow, gcol
      integer d(4), f(6), nfock
      integer g_a_dens, g_a_coul, g_a_exch, g_a_xc
      integer g_b_dens, g_b_coul, g_b_exch, g_b_xc
cc AJL/Begin/SPIN ECPs
c      integer g_hcore
      integer g_a_hcore, g_b_hcore
c AJL/End
      integer g_tmp(6),ifock
      double precision jfac(4), kfac(4), one, zero, mone
      parameter (one=1.0d0, zero=0.0d0, mone=-1.0d0)
      double precision e_a_coul, e_a_exch, e_b_coul, e_b_exch,
     &     e_a_xc, e_b_xc
      double precision errmaxa, errmaxb
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical odebug
      logical cphf_uhf
      double precision ecosmo !< [Output] The COSMO solvation energy
      integer g_dens(2),i
      character*255 cosmo_file
c     DIM/QM JEM
      logical ldimqm
      integer g_vdim
      double precision edimqm
c
cc AJL/Begin/SPIN ECPs
      integer ecp_handle
      integer ecp_channels
      logical bas_get_ecp_handle
      external bas_get_ecp_handle
      logical ecp_get_high_chan
      external ecp_get_high_chan
cc AJL/End
c     
c     Check
c     
      odebug = util_print('uhf_debug', print_debug)
      call uhf_jkfac(jfac,kfac)
      if (.not.cuhf_init_flag)
     $     call errquit('uhf_energy: UHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_grad, gtype, grow, gcol)
      if ((grow.ne.cuhf_vlen).or.(gcol.ne.1))
     $     call errquit('uhf_energy: invalid vector length',0, GA_ERR)
      cphf_uhf = .false.
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_uhf', mt_log, 1, cphf_uhf)) then
         cphf_uhf = .false.
      endif
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
      if (nalpha .gt. 0) then
         call ga_dgemm('n', 't', nbf, nbf, nalpha, one, g_vecs(1),
     $        g_vecs(1), zero, g_a_dens)
      else
         call ga_zero(g_a_dens)
      endif
      if (nbeta .gt. 0) then
         call ga_dgemm('n', 't', nbf, nbf, nbeta, one, g_vecs(2),
     $        g_vecs(2), zero, g_b_dens)
      else
         call ga_zero(g_b_dens)
      endif
c
c     Since UHF can break spatial symmetry by localizing the orbitals
c     the densities may not be totally symmetric, but since the Hamiltonian
c     is symmetric contraction with the integrals projects out the totally 
c     symmetric component ... hence we can symmetrize the densities and
c     exploit symmetry.  Compute the max change in any element due to
c     symmetrizing and print a warning if it is big.
c
c     !! If this is the case then where does the 'force' for symmetry breaking
c     come from?  Must be missing something?
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
 
      call ga_zero(g_a_coul)
      call ga_zero(g_b_coul)
      call ga_zero(g_a_exch)
      call ga_zero(g_b_exch)
      if(cphf_uhf)then
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
      if(cphf_uhf)then
         f(5) = g_a_xc
         f(6) = g_b_xc
      endif
c
czz   two extra ga's are passed to fock_2e to get the xc matrix
c
      if(cphf_uhf)then
         nfock = 6
      else
         nfock = 4
      endif
      call do_riscf (.false.)
      if (.not.cam_exch) then
        call fock_2e(geom, basis, nfock, jfac, kfac, tol2e,
     $     oskel, d, f, .false.)
      else ! for attenuated calculations
c
c       calculate the CAM exchange
c
        do ifock = 1,nfock
           g_tmp(ifock) = ga_create_atom_blocked(geom, basis, 'tmp')
           call ga_zero(g_tmp(ifock))
        end do
c
        call case_setflags(.true.)  ! set LC flag
        jfac(1)=0d0
        jfac(2)=0d0
        jfac(3)=0d0
        jfac(4)=0d0
        kfac(1)=0d0
        kfac(2)=1d0
        kfac(3)=0d0
        kfac(4)=1d0
        call fock_2e_cam(geom, basis, nfock, jfac, kfac, tol2e, 
     &     oskel, d, g_tmp, .false., .false.)
        call ga_dadd(1d0,f(2),1d0,g_tmp(2),f(2))  ! LC exchange part
        call ga_dadd(1d0,f(4),1d0,g_tmp(4),f(4))  ! LC exchange part
c
c       calculate the full Coulomb
c
        do ifock = 1,nfock
           call ga_zero(g_tmp(ifock))
        end do
        call case_setflags(.false.) ! turn off LC flag for full Coulomb
        jfac(1)=1d0
        jfac(2)=0d0
        jfac(3)=1d0
        jfac(4)=0d0
        kfac(1)=0d0
        kfac(2)=0d0
        kfac(3)=0d0
        kfac(4)=0d0
        call fock_2e_cam(geom, basis, nfock, jfac, kfac, tol2e, 
     &     oskel, d, g_tmp, .false., .true.)   ! last argument toggles xc
        call ga_dadd(1d0,f(1),1d0,g_tmp(1),f(1))   ! full Coulomb part
        call ga_dadd(1d0,f(3),1d0,g_tmp(3),f(3))   ! full Coulomb part
        call ga_dadd(1d0,f(2),1d0,g_tmp(2),f(2))   ! DFT xc part
        call ga_dadd(1d0,f(4),1d0,g_tmp(4),f(4))   ! DFT xc part
c
c       destroy work space
        if (.not. ga_destroy(g_tmp)) 
     &   call errquit('uhf: ga corrupt?',0, GA_ERR)
      end if  ! cam_exch
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
      if(cphf_uhf)then
         call ga_dadd(one, g_a_coul, one, g_a_xc, g_a_coul)
         call ga_dadd(one, g_b_coul, one, g_b_xc, g_b_coul)
      endif
c
c     reuse g_a_exch to hold the 1-e integrals
c
cc AJL/Begin/SPIN ECPs
c      g_hcore = g_a_exch
c      call ga_zero(g_hcore)

      ecp_channels = 1 
      if (bas_get_ecp_handle(basis,ecp_handle)) then 
        if (ecp_get_high_chan(ecp_handle,ecp_channels)) 
     &    continue
      end if

      g_a_hcore = g_a_exch
      g_b_hcore = g_b_exch
      call ga_zero(g_a_hcore)
      call ga_zero(g_b_hcore)
c
      call int_1e_ga(basis, basis, g_a_hcore, 'kinetic', oskel)  ! kinetic
      if (do_zora .and. .not.(do_NonRel)) then
c     is uhf() called from cphf_solve2 (say from hnd_hyperfine_zora)?
c     hack: use 'cphf:skew' to check this condition
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1,
     $        noskew_uhf)) then
            call errquit('ZORA not available yet for SCF ',0,CAPMIS_ERR)
         endif
        call ga_dadd(1.d0,g_a_hcore,1.d0,g_zora_Kinetic(1),g_a_hcore) ! zora kinetic
      endif

cc If spin polarised ECP, split g_hcore
      if (ecp_channels.gt.1)
     &  call ga_copy(g_a_hcore,g_b_hcore)

cc Put normal potential in to alpha hcore
c      call int_1e_ga(basis, basis, g_hcore, 'potential', oskel) !potential
      call int_1e_ga(basis, basis, g_a_hcore, 'potential', oskel) !potential
      if (ecp_channels.gt.1) then ! Put beta potential in to beta hcore
        call int_1e_ga(basis, basis, g_b_hcore, 'potential_beta', oskel)!potential_beta
      else ! it is spin independent, so copy alpha to beta
        call ga_copy(g_a_hcore,g_b_hcore)
      endif
cc AJL/End
c
c     cosmo charges, potential and energy contribution
c
      if (cosmo_on) then
         call errquit('uhf_energy not compatible with COSMO',0,
     C     CAPMIS_ERR)
      endif
      if(cosmo_on.and.cosmo_phase.eq.2) then
        cosmo_file = "cosmo.xyz"
        do i = 1, 2
         g_dens(i) = ga_create_atom_blocked(geom,basis,'density matrix')
         call ga_zero(g_dens(i))
        end do
        call ga_copy(g_a_dens,g_dens(1))  ! alpha = 1
        call ga_copy(g_b_dens,g_dens(2))  ! beta = 2
        call cosmo_charges_from_dmat(rtdb, basis, geom, ecosmo, odebug, 
     &                                  2,           ! 2 = open shell 
     &                                  g_dens,      ! input density
     &                                  cosmo_file)  ! cosmo charges file name
cc AJL/Begin/SPIN ECPs
        call int_1e_ga(basis,basis,g_a_hcore,'cos_chg_pot',.false.)
        call int_1e_ga(basis,basis,g_b_hcore,'cos_chg_pot',.false.)
cc AJL/End
        do i = 1, 2
         if (.not.ga_destroy(g_dens(i)))
     &    call errquit('uhf_energy: ga_destroy failed g_dens',0,GA_ERR)
        enddo
      endif  ! cosmo check
c
c     DIM/QM JEM
c     Calculate DIM potential
      if (.not. rtdb_get(rtdb, 'dimqm:lrsp', mt_log, 1, ldimqm))
     $  ldimqm = .false.
      if (ldimqm) then
c       We need to combine the alpha and beta densities to calculate the total DIM dipoles
        g_dens(1) = ga_create_atom_blocked(geom,
     $                                     basis, 'dim density matrix')
        call ga_copy(g_a_dens, g_dens(1))
        call ga_dadd(one, g_dens(1), one, g_b_dens, g_dens(1))
        call dimqm_rohf_wrap(rtdb, geom, basis, nbf, g_dens(1), g_vdim)
        if(.not.ga_destroy(g_dens(1)))
     $    call errquit('uhf_energy: ga_destroy failed DIM den',0,GA_ERR)
      end if
c
c     add in frozen embedding
      if (frozemb) then
cc AJL/Begin/Spin ECPs
        call ga_dadd(1.d0,g_a_hcore,1.d0,g_frozemb,g_a_hcore)
        call ga_dadd(1.d0,g_b_hcore,1.d0,g_frozemb,g_b_hcore)
cc AJL/End
      end if
c
      call ga_sync()
c
      eone = 
     $     (ga_ddot(g_a_dens,g_a_hcore) + ga_ddot(g_b_dens,g_b_hcore))
      call ga_dadd(one, g_a_hcore, one, g_a_coul, g_a_coul)
      call ga_dadd(one, g_b_hcore, one, g_b_coul, g_b_coul)
c
c     DIM/QM JEM
c     Add DIM potential
      if (ldimqm) then
        call ga_dadd(one, g_vdim, one, g_a_coul, g_a_coul)
        call ga_dadd(one, g_vdim, one, g_b_coul, g_b_coul)
      end if
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
c      DIM/QM JEM
      if (ldimqm) then
        if(.not.ga_destroy(g_vdim))
     $    call errquit('uhf_energy: destroy vdim failed',0,GA_ERR)
      end if
c
c     extract the gradient
c
      call uhf_get_grad(g_grad)
c
      if (odebug) call ga_print(g_grad)
c
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $     call errquit('uhf_energy: no repulsion energy?', 0, GEOM_ERR)
      energy = eone + etwo + enrep
c
c     cosmo energy term 
c
      if(cosmo_on.and.cosmo_phase.eq.2) then
         if (.not. rtdb_get(rtdb,'cosmo:energy',mt_dbl,1,ecosmo))
     $     call errquit('rohf_energy: rtdb get failed for ecosmo',911,
     &       RTDB_ERR)
         energy=energy+ecosmo
      endif  ! cosmo check
c
c     DIM/QM JEM
c     DIM energy term
      if (ldimqm) then
        if (.not.rtdb_get(rtdb, 'dimqm:edimqm', mt_dbl, 1, edimqm))
     $    call errquit('uhf_energy get edimqm failed', 1, RTDB_ERR)
        energy = energy + edimqm
      end if
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
#include "mafdecls.fh"
      integer g_grad
c
      integer ioff
      integer nvir, ivoff
      double precision maxelt
c
      call ga_copy_patch('n',
     $     cuhf_g_falpha, nalpha+1, nmo, 1, nalpha,
     $     g_grad, 1, (nmo-nalpha)*nalpha, 1, 1)
c
      ioff = (nmo-nalpha)*nalpha
      call ga_copy_patch('n',
     $     cuhf_g_fbeta, nbeta+1, nmo, 1, nbeta,
     $     g_grad, ioff+1, ioff+(nmo-nbeta)*nbeta, 1, 1)
c
      call ga_dscal(g_grad, 2.0d0)
c
      if(oadapt) then
        nvir = nmo - nalpha
        call scf_sym_screen(nalpha, 0, nvir, int_mb(k_irs),
     $       g_grad, 0, .true., maxelt)
        nvir = nmo - nbeta
        ivoff = nalpha*(nmo-nalpha)
        call scf_sym_screen(nbeta, 0, nvir, int_mb(k_irs+nmo),
     $       g_grad, ivoff, .true., maxelt)
      endif
c
      end
      subroutine uhf_search_precond(rtdb, g_grad, g_work)
      implicit none
#include "cuhf.fh"      
#include "cscfps.fh"
#include "global.fh"
#include "cscf.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      integer g_grad
      integer g_work
c
c     Precondition the gradient with an approximation to the Hessian
c     return the result in g_work.
c
c     Old version only divided by the level-shifted diagonals which 
c     is in rohf the preconditioner to the preconditioner!
c
      logical oprint_conv, oprint_eval, oprint_debug, oprint_solve
      double precision gnorm
      double precision max_grad
      double precision min_shift
      integer max_precond_iter  ! Max. no. of iters for preconditioner
      integer max_precond_sub   ! Max. dim. of precond. iterative space
      parameter (max_precond_iter=20)
c
      double precision precond_acc ! Accuracy attained by solver
      integer precond_iter       ! No. of iterations used by solver
c
      integer nvir, ivoff
      double precision maxelt
c
      external uhf_precond, uhf_hessv
      logical ga_iter_lsolve
      external ga_iter_lsolve
c
      oprint_conv = util_print('convergence', print_default)
      oprint_eval = util_print('intermediate evals', print_debug)
      oprint_debug= util_print('uhf_debug', print_debug)
      oprint_solve= util_print('solve', print_high)
c
      if (.not.rtdb_get(rtdb, 'scf:maxsub', MT_INT, 1, 
     $     max_precond_sub)) max_precond_sub = 10
c
      call ga_maxelt(g_grad, max_grad)
      gnorm = sqrt(ga_ddot(g_grad, g_grad))
      eprec = max(1d-7,gnorm*0.01d0, tol2e*100.0d0) ! Approx precision available
c
c     Determine shift necessary to make the preconditioner
c     positive definite
c
 30   continue
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
         if (oadapt) then
           call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     $           0, 'after canon', .true., int_mb(k_irs))
           call scf_movecs_sym_adapt(basis, g_movecs(2), .false.,
     $           0, 'after canon', .true., int_mb(k_irs+nmo))
         endif
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
c     Next to lines precondition with just the eigenvalue differences
c
*      call ga_copy(g_grad,g_work)
*      call uhf_precond(g_work, 0.0d0)
c
      call ga_zero(g_work)      ! ESSENTIAL
c     
c     Attempt to solve the damned equations
c
*      write(6,*) ' pflg ', pflg
c
 20   if (.not. ga_iter_lsolve(cuhf_vlen, max_precond_iter,
     $     max_precond_sub, precond_tol, 
     $     uhf_precond, uhf_hessv, .true., oprint_solve, 
     $     g_grad, g_work, precond_acc, precond_iter)) then
c     
c     Iterative solution did not converge.  If we got at least
c     1 sig. fig. then just continue.  Otherwise for the 1-e
c     approximation increase the level-shift and restart with
c     the current vector as guess.  If the 2-e approximation did
c     not converge disable it forever and resort to the 1-e method.
c     
         if (precond_acc .gt. 0.1d0) then
            if (pflg .eq. 2) then
               odisable_nr = .true.
               maxiter = maxiter + 10
               if (ga_nodeid() .eq. 0 .and. oprint_conv) then
                  write(6,22) maxiter
 22               format(/' Disabled NR: increased maxiter to ',i3/)
                  call util_flush(6)
               endif
               goto 30
            else
               if (lshift .eq. 0.0d0) then
                  lshift = 2.0d0
                  call ga_zero(g_work)
               else
                  lshift = lshift*2.0d0
               endif
            endif
            if (ga_nodeid() .eq. 0 .and. oprint_conv) then
               write(6,2) lshift
 2             format(' Increased level shift to ', f8.2)
               call util_flush(6)
            endif
            goto 20
         endif
      endif
c
      if(oadapt) then
        nvir = nmo - nalpha
        call scf_sym_screen(nalpha, 0, nvir, int_mb(k_irs),
     $       g_work, 0, .true., maxelt)
        nvir = nmo - nbeta
        ivoff = nalpha*(nmo-nalpha)
        call scf_sym_screen(nbeta, 0, nvir, int_mb(k_irs+nmo),
     $       g_work, ivoff, .true., maxelt)
      endif
c
      end
      subroutine uhf_precond(g_x, solveshift)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
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
      integer g_f(2), ivec, nvec, gtype, vlen
      double precision diag, denominator, shift, dnrm
      integer ilo(2), ihi(2)
c
      logical oprint, olprint
c
      diag(i) = dbl_mb(k_diag + i - 1)
c     
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta
c
      oprint = util_print('precond',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: tmp',l_diag,k_diag))
     $     call errquit('uhf_precond: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: tmp',l_x,k_x))
     $     call errquit('uhf_hdiag_scale: cannot allocate',0, MA_ERR)
c
      call ga_sync()
      call ga_inquire(g_x, gtype, vlen, nvec)
c
      shift = lshift - solveshift
c
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"uhf: in g_x = ",i4,f24.8)')
     +      ivec,dnrm
          endif
        enddo
      endif
c
      ibase = 1
      do iset = 1, 2
         nvir = nmo - nocc(iset)
         if (nvir .gt. 0) then
            call ga_get_diagonal(g_f(iset), dbl_mb(k_diag))
            do ivec = 1, nvec
               do i = ga_nodeid()+1, nocc(iset), ga_nnodes()
                  ioff = (i-1)*nvir + ibase
                  call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),nvir)
                  do j=1,nvir
                     denominator = 2.d0 * (diag(j+nocc(iset)) - diag(i))
     $                    + shift
                     if (denominator .lt. 0.1d0) denominator = 0.1d0
                     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
                  end do
                  call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),nvir)
               end do
            enddo
         endif
         ibase = ibase + nocc(1)*(nmo-nocc(1))
      end do
c
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"uhf: out g_x = ",i4,f24.8)')
     +      ivec,dnrm
          endif
        enddo
      endif
c
      if (.not. ma_pop_stack(l_x)) call errquit('uhf:pop x',0, MA_ERR)
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf:pop ',0, MA_ERR)
      call ga_sync()
c
      end
      subroutine uhf_make_shift(shift)
      implicit none
#include "errquit.fh"
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
      if (.not. ma_push_get(mt_dbl, nmo, 'uhf_m_s: diag', 
     $     l_diag, k_diag)) call errquit
     $     ('uhf_make_shift: insufficient memory', nmo, MA_ERR)
c
      call ga_sync()
c
      shift_a = 0.0d0
      shift_b = 0.0d0
c
      call rohf_make_shift_2(cuhf_g_falpha, dbl_mb(k_diag),
     $     1, nalpha, nalpha+1, nmo, oadapt, int_mb(k_irs),
     $     2.0d0, shift_a)
      call rohf_make_shift_2(cuhf_g_fbeta, dbl_mb(k_diag),
     $     1, nbeta, nbeta+1, nmo, oadapt, int_mb(k_irs+nmo),
     $     2.0d0, shift_b)
c
      shift = max(shift_a, shift_b)
c
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf_n_s: e',0,
     &       MA_ERR)
c
      end
      subroutine uhf_canon(oaufbau, oprint)
C$Id$
      implicit none
#include "errquit.fh"
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
c        diagonalize the whole thing and allow mixing of occ-virt
c     else
c        diagonalize separately the occ-occ, and virt-virt parts
c
c     Transform Fock matrices and MO coefficients into the new canonical basis
c
      if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_canon: u',
     $     32, 32, g_u)) call errquit('uhf_canon: ga failed for u', 0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_canon: fock',
     $     32, 32, g_fock)) call errquit
     $     ('uhf_canon: ga failed for fock', 0, GA_ERR)
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
#ifdef SCALAPACK
            call ga_pdsyev(g_fock, g_u, dbl_mb(k_eval+ioff), 0)
#else
            call ga_diag_std(g_fock, g_u, dbl_mb(k_eval+ioff))
#endif
#else 
            call ga_diag_std_seq(g_fock, g_u, dbl_mb(k_eval+ioff))
#endif
         else
c     
c     occ-occ piece
c
            if (nocc(iset) .gt. 0) call rohf_canon_subspace
     $           (g_fock, g_u, dbl_mb(k_eval+ioff),
     $           1, nocc(iset), int_mb(k_irs+(iset-1)*nmo))
c
c     virt-virt piece
c     
            if (nmo-nocc(iset) .gt. 0)
     $           call rohf_canon_subspace(g_fock, g_u,
     $           dbl_mb(k_eval+ioff+nocc(iset)),
     $           nocc(iset)+1, nmo, int_mb(k_irs+(iset-1)*nmo))

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
         if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_canon: tmp',
     $        32, 32, g_tmp)) call errquit
     $        ('uhf_canon: ga failed for tmp', 0, GA_ERR)
         call ga_copy(g_movecs(iset), g_tmp)
         call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_tmp, g_u,
     $        zero, g_movecs(iset))
c
         if (nbf .ne. nmo) then
            if (.not. ga_destroy(g_tmp))
     $           call errquit('uhf_canon: destroy', 0,
     &       GA_ERR)
            if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_canon: tmp',
     $           32, 32, g_tmp)) call errquit
     $           ('uhf_canon: ga failed for tmp', 0, GA_ERR)
         endif
c
         call ga_dgemm('n', 'n', nmo, nmo, nmo, one, g_f(iset), g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nmo, nmo, nmo, one, g_u, g_tmp, 
     $        zero, g_f(iset))
c
         if (.not. ga_destroy(g_tmp))
     $        call errquit('uhf_canon: destroy', 0, GA_ERR)
      end do
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6, 'Alpha-spin eigenvalues',
     $        20, .true.)
         call output(dbl_mb(k_eval), 1, min(nalpha+5,nmo), 
     $        1, 1, nmo, 1, 1)
         write(6,*)
         call util_print_centered(6, 'Beta-spin eigenvalues',
     $        20, .true.)
         call output(dbl_mb(k_eval+nbf), 1, min(nalpha+5,nmo), 
     $        1, 1, nmo, 1, 1)
         call util_flush(6)
      end if
c
      if (.not. ga_destroy(g_u))
     $     call errquit('uhf_canon: destroy', 0, GA_ERR)
      if (.not. ga_destroy(g_fock))
     $     call errquit('uhf_canon: destroy', 0, GA_ERR)
c
      end
      subroutine uhf_spin(sz, s2)
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      double precision sz, s2
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c     Return the expection values of Sz and S^2
c
      integer g_smo, g_over, g_work, l_x, k_x, i, j
c
      sz = 0.5d0 * (nalpha - nbeta)
c     
c     Form Sab
c
      g_over = ga_create_atom_blocked(geom, basis, 'uhf_spin: over')
*      if (.not. ga_create(MT_DBL, nbf, nbf, 'uhf_spin: overlap',
*     $     32, 32, g_over)) call errquit('uhf_spin: overlap', 0)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'uhf_spin: work',
     $     32, 32, g_work)) call errquit('uhf_spin: work', 0, GA_ERR)
      if (.not. ga_create(MT_DBL, nmo, nmo, 'uhf_spin: overlap',
     $     32, 32, g_smo)) call errquit('uhf_spin: overlap', 0, GA_ERR)
c      
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      call two_index_transf(g_over, g_movecs(1), g_movecs(2),
     $     g_work, g_smo)
      if (.not. ga_destroy(g_work)) call errquit('uhf_spin:ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_over)) call errquit('uhf_spin:ga?',0,
     &       GA_ERR)
c
c     Sum overlap of occupied orbitals
c
      if (.not. ma_push_get(mt_dbl,nalpha,'uhf_spin:x',l_x,k_x))
     $     call errquit('uhf_spin: ma failed on x',  nalpha, MA_ERR)
c
      s2 = 0
      do i = 1+ga_nodeid(), nbeta, ga_nnodes()
         call ga_get(g_smo, 1, nalpha, i, i, dbl_mb(k_x), nalpha)
         do j = 1, nalpha
            s2 = s2 + dbl_mb(k_x+j-1)**2
         enddo
      enddo
      call ga_dgop(msg_uhf, s2, 1, '+')
      s2 = sz*(sz+1) + nbeta - s2
c
      if (.not. ma_pop_stack(l_x)) call errquit('uhf_spin:ma?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_smo)) call errquit('uhf_spin:ga?',0,
     &       GA_ERR)
c
      end
      subroutine uhf_analyze(rtdb)
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
      integer rtdb
c
c     Analyze the UHF wavefunction.  Optionally print the mulliken
c     analysis and multipoles, always store the dipole in the RTDB
c     
      integer g_dens, g_adens, g_bdens, g_over
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oprintmulliken, oprintmultipole,oprintinertia
c
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments) 
      double precision alphamoments (len_moments)
      double precision betamoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      data center/3*0.0d0/
c
      oprintmulliken  = util_print('mulliken', print_default)
      oprintmultipole = util_print('multipole', print_default)
      oprintinertia = util_print('inertia', print_default).and.
     . ga_nodeid().eq.0
c
c     moments of inertia
c
      if(oprintinertia) call geom_momint(geom)

c     
c     First analyze the total density
c     
      g_dens = ga_create_atom_blocked(geom, basis, 'uhf_mull:dens')
      g_adens = ga_create_atom_blocked(geom, basis, 'uhf_mull:dens')
      g_bdens = ga_create_atom_blocked(geom, basis, 'uhf_mull:dens')
      g_over  = ga_create_atom_blocked(geom, basis, 'uhf_mull:over')
      call ga_zero(g_dens)
      call ga_zero(g_adens)
      call ga_zero(g_bdens)
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
c     
      if (nalpha .gt. 0) 
     $    call ga_dgemm('n', 't', nbf, nbf, nalpha, 1.0d0, g_movecs(1),
     $         g_movecs(1), 0.0d0, g_adens)
      if (nbeta .gt. 0) 
     $    call ga_dgemm('n', 't', nbf, nbf, nbeta, 1.0d0, g_movecs(2),
     $         g_movecs(2), 0.0d0, g_bdens)
      call ga_dadd(1.0d0, g_adens, 1.0d0, g_bdens, g_dens)
c
      if (oprintmulliken) then
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,
     $           'Mulliken analysis of the total density', 20,.true.)
         endif
         call mull_pop(geom, basis, g_dens, g_over, 'total')
c
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,
     $           'Mulliken analysis of the alpha density', 20,.true.)
         endif
         call mull_pop(geom, basis, g_adens, g_over, 'alpha')
c     
         if (nbeta .gt. 0) then
            if (ga_nodeid() .eq. 0) then
               write(6,*)
               call util_print_centered(6,
     $              'Mulliken analysis of the beta density', 
     $              20,.true.)
            endif
            call mull_pop(geom, basis, g_bdens, g_over, 'beta')
c
            if (ga_nodeid() .eq. 0) then
               write(6,*)
               call util_print_centered(6,
     $              'Mulliken analysis of the spin density', 
     $              20,.true.)
            endif
            call ga_dadd(1.0d0, g_adens, -1.0d0, g_bdens, g_adens)
            call mull_pop(geom, basis, g_adens, g_over, 'spin')
            call ga_dadd(1.0d0, g_adens,  1.0d0, g_bdens, g_adens)
         endif
      endif
c
      if (.not. bas_cando_mpoles(basis)) goto 100
c     
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, alphamoments, 1)
      call dfill(len_moments, 0.0d0, betamoments, 1)
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, alphamoments, 1)
      call dfill(len_moments, 0.0d0, betamoments, 1)
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      if (nalpha .gt. 0)
     $     call multipole_density(basis, center, lmax, g_adens, 
     $     alphamoments, len_moments)
      if (nbeta .gt. 0) 
     $     call multipole_density(basis, center, lmax, g_bdens, 
     $     betamoments, len_moments)
      do k = 1, len_moments
         totalmoments(k) = alphamoments(k) + betamoments(k) + 
     $        nuclmoments(k)
      enddo
c
      if (oprintmultipole .and. ga_nodeid().eq.0) then
         write(6,*)
         call util_print_centered(6,
     $     'Multipole analysis of the density wrt the origin',
     $      30, .true.)
         write(6,*)
         write(6,*) '    L   x y z        total         alpha',
     $                 '         beta         nuclear'
         write(6,*) '    -   - - -        -----         -----',
     $                 '         ----         -------'
         ind = 0
         do ltotal = 0, lmax
            do k = ltotal, 0, -1
               do l = ltotal-k, 0, -1
                  m = ltotal - k - l
                  ind = ind + 1
                  write(6,12) ltotal, k, l, m, totalmoments(ind),
     $                 alphamoments(ind), betamoments(ind),
     $                 nuclmoments(ind)
 12               format(4x,i2,2x,3i2,4f14.6)
               enddo
            enddo
            write(6,*)
         enddo
      endif
c
      call ecce_print1('total dipole', mt_dbl, totalmoments(2), 3)
      call ecce_print1('alpha electronic dipole', mt_dbl, 
     $     alphamoments(2), 3)
      call ecce_print1('beta electronic dipole', mt_dbl, 
     $     betamoments(2), 3)
      call ecce_print1('nuclear dipole', mt_dbl, nuclmoments(2), 3)
      call ecce_print1('total quadrupole', mt_dbl, totalmoments(5), 6)
      call ecce_print1('alpha electronic quadrupole',mt_dbl,
     $     alphamoments(5),6)
      call ecce_print1('beta electronic quadrupole',mt_dbl,
     $     betamoments(5),6)
      call ecce_print1('nuclear quadrupole', mt_dbl, nuclmoments(5), 6)
c
c     Store the dipole moment
c
      if (.not. rtdb_put(rtdb, 'scf:dipole', mt_dbl, 
     $     3, totalmoments(2))) call errquit('multipole: rtdb?',0,
     &       RTDB_ERR)
c     
 100  if (.not. ga_destroy(g_over)) call errquit('uhf_mull: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_adens)) call errquit('uhf_mull: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_bdens)) call errquit('uhf_mull: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dens)) call errquit('uhf_mull: ga?',0,
     &       GA_ERR)
c     
      end




