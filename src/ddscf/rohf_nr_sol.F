      logical function rohf_nr_solve(rtdb, energy, eone, etwo, enrep)
C     $Id: rohf_nr_sol.F,v 1.2 1996-03-13 22:22:59 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "crohf.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "util.fh"
#include "cscf.fh"
c     
c     Solve the ROHF equations using a hybrid NR/PCG method
c     
c     Arguments
c     
      integer rtdb
      double precision energy   ! Return ROHF energy
      double precision eone, etwo, enrep ! Return ROHF energy contribs
c     
c     Local GA handles 
c     
      integer g_grad            ! gradient vector  ... crohf_vlen
      integer g_search          ! search direction ... crohf_vlen
      integer g_work            ! scratch vector   ... crohf_vlen
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
      logical ga_iter_lsolve
      external ga_create_atom_blocked, ga_iter_lsolve
      external rohf_hessv, rohf_hessv_precond
c     
      if (ga_nodeid().eq.0) call util_flush(6)
c     
c     Allocate global arrays.  Gradient, search direction
c     and work space for the PCG routine.
c     
      if (.not. ga_create(MT_DBL, crohf_vlen, 1, 'rohf_nr_solve: grad',
     $     1, 1, g_grad)) call errquit('rohf_nr_solve: ga_cre grad', 0)
      if (.not. ga_create(MT_DBL, crohf_vlen, 1, 'rohf_nr_solve: serch',
     $     1, 1, g_search)) call errquit('rohf_nr_solve: ga_cre grad',0)
      if (.not. ga_create(MT_DBL, crohf_vlen, 1, 'rohf_nr_solve: work',
     $     1, 1, g_work)) call errquit('rohf_nr_solve: ga_cre work', 0)
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
     $        1x,'        Quadratically convergent ROHF',//,
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
      endif
c     
c     PCG iterative loop
c     
      mode = 'startup'
      converged = .false.
 10   if (.not. converged) then
         call ga_pcg_minim(crohf_vlen, iter, energy, g_grad, g_work,
     $        g_search, step, ls_tol, ls_max, eprec, oconjugacy, 
     $        oprint_conv, oprint_ls, iwork, dwork, mode)
c     
         if (mode .eq. 'energy+gradient') then
c     
c     Compute the energy and gradient at step*search direction
c     
            call rohf_step_energy(rtdb, step, g_search,
     $           eone, etwo, enrep, energy, g_grad)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            goto 10
c     
         else if (mode .eq. 'precondition') then
c     
c     Precondition the gradient direction with approx./exact hessian
c     
            call rohf_search_precond(rtdb, g_grad, g_work)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            goto 10
c     
         else if (mode .eq. 'accept step') then
c     
c     Apply the rotation and check for overall convergence.
c     The current fock matrices and energy correspond to the
c     updated MO coefficients.
c     
            call rohf_rotate_movecs(step, g_search)
            converged = gnorm .lt. gnorm_tol
            if ((iter.le.maxiter).and.(.not. converged)) then
               if (mod(iter,5).eq.0) 
     $              call ga_orthog_mos(basis, g_movecs)
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
                  endif
                  call ga_print(g_movecs)
               endif
               call scf_movecs_write
               goto 10
            endif
c            
         endif
      endif
c     
c     End SCF minimisation
c     
      if (ga_nodeid().eq.0) call util_flush(6)
c     
c     *ALWAYS* return canonical MOs ... don't do aufbau
c     
      call rohf_canon(.false., oprint_eval)
      call movecs_fix_phase(g_movecs)
      if (olock) call scf_movecs_lock
      call scf_movecs_write
      if (olagr) call rohf_lagr()
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('rohf_nr_solve: ga_destroy grad', 0)
      if (.not. ga_destroy(g_search)) call errquit
     $     ('rohf_nr_solve: ga_destroy search', 0)
      if (.not. ga_destroy(g_work)) call errquit
     $     ('rohf_nr_solve: ga_destroy work', 0)
c     
      if (ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c     
      rohf_nr_solve = converged
c     
      end
