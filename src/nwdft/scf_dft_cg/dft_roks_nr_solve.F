      logical function dft_roks_nr_solve(rtdb, energy, eone, etwo, exc,
     +                                   enrep, edisp, nexc)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "crohf.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "util.fh"
#include "cscf.fh"
#include "cfock.fh"
#include "stdio.fh"
#include "cosmo.fh"
c     
c     Solve the ROKS equations using a hybrid NR/PCG method
c     
c     Arguments
c     
      integer rtdb
      double precision energy   ! Return ROKS energy
      double precision eone, etwo, enrep ! Return ROKS energy contribs
      double precision exc(2) ! Return DFT energy contribs
      double precision edisp ! Dispersion correction
      integer nexc
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
      parameter (ls_max=20)
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
      logical ovariable
c
      double precision tol2e_save ! Used with variable integral
                                  ! tolerance
*     double precision djunk
      double precision eprec_tol
      double precision tlastwrite,ttest ! time last wrote the MO vectors
      double precision time_iter_start
      double precision time1_rot, time2_rot
      integer time_prev_iter
c     
c     Functions
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      external rohf_hessv, rohf_hessv_precond
c
c     cosmo code
c
      logical odebug
      logical osome
c
      double precision egas
      double precision esol
c
      if (ga_nodeid().eq.0) call util_flush(LuOut)
c
c     Determine if using variable precision.  Turn it off if
c     there is the slightest hint we might have problems
c
c     TURN OFF VARIABLE PRECISION FOCK BUILD FOR NOW
c
      ovariable = .false. 
*      if (.not. rtdb_get(rtdb, 'scf:variable', mt_log, 1, ovariable))
*     $     ovariable = (nmo .eq. nbf)  ! Need high acc. if linear dependence
*      if (rtdb_get(rtdb,'tolguess',mt_dbl, 1, djunk) .or.
*     $     rtdb_get(rtdb,'fock:dentolmax',mt_dbl,1,djunk))
*     $     ovariable = .false.
*      if (((nmo-nclosed).eq.0) .or. ((nmo-nopen).eq.0))
*     $     ovariable = .false.
c
c     
c     Allocate global arrays.  Gradient, search direction
c     and work space for the PCG routine.
c     
*ga:1:0
      if (.not. ga_create(MT_DBL, crohf_vlen, 1, 
     $          'dft_roks_nr_solve: grad', 0, 1, g_grad))
     $    call errquit('dft_roks_nr_solve: ga_cre grad', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, crohf_vlen, 1, 
     $          'dft_roks_nr_solve: serch', 0, 0, g_search))
     $    call errquit('dft_roks_nr_solve: ga_cre grad',0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, crohf_vlen, 1,
     $          'dft_roks_nr_solve: work', 0, 0, g_work))
     $    call errquit('dft_roks_nr_solve: ga_cre work', 0, GA_ERR)
c
      tlastwrite = util_wallsec()
c
      oprint_parm = util_print('parameters', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_ls   = util_print('line search', print_high)
      oprint_vecs = util_print('intermediate vectors', print_debug)
      oprint_eval = util_print('intermediate evals', print_debug)
c
      if (ga_nodeid().eq.0.and. oprint_parm) then
         write(LuOut,1) gnorm_tol, maxiter, tol2e
         if (ouser_changed_conv) write(LuOut,11) shifts,nr_gswitch
         write(LuOut,111)
 1       format(//,
     $        1x,'----------------------------------------------',/
     $        1x,'        Quadratically convergent ROKS',//,
     $        1x,'Convergence threshold     :',9x,1p,e10.3,0p,/,
     $        1x,'Maximum no. of iterations :',9x,i4,/,
     $        1x,'Final Fock-matrix accuracy:',9x,1p,e10.3,0p)
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
         call util_flush(LuOut)
      endif
c
      call ga_sync()
c
c     ----- -nr- options ... -----
c
      ododiag = .true.
      odisable_nr = .false.
c
c     Only use low accuracy at outset for startup jobs
c
      tol2e_save = tol2e
c
      if (ovariable) then
         if (.not. (nelec.eq.1 .or. owritefile .or. 
     $        movecs_guess.eq.'restart')) 
     $        tol2e = max(1d-7,tol2e) ! Initial low accuracy
      endif
      gnorm = 0.0d0             ! So is not undefined
      eprec_tol = max(1d-7,100.0d0*tol2e) !To avoid noise in line search
      eprec = eprec_tol
c     
c     PCG iterative loop
c     
      time_iter_start = util_wallsec()
c
      mode = 'startup'
      converged = .false.
c
   9  continue
      odebug=.false.
      odebug=odebug.and.ga_nodeid().eq.0
      if(cosmo_on.and.cosmo_phase.eq.1) then
         osome=ga_nodeid().eq.0
         if(osome) then
            write(LuOut,909)
         endif
         if(odebug) then
            write(LuOut,*) 'mode = ',mode
         endif
      elseif(cosmo_on.and.cosmo_phase.eq.2) then
         iter     = 0
         energy   = 0.0d0
         eone     = 0.0d0
         etwo     = 0.0d0
         mode     = 'startup'
         converged=.false.
         osome=ga_nodeid().eq.0
         if(osome) then
            write(LuOut,910)
         endif
         if(odebug) then
            write(LuOut,*) 'mode = ',mode
         endif
      endif
      if (ga_nodeid().eq.0) call util_flush(6)
 909    format(2x,'   COSMO gas phase')
 910    format(2x,'   COSMO solvation phase')
c
 10   if (.not. converged) then
c
c        If running disk resident, and especially if semi-direct
c        then variable accuracy hinders things
c
         if (oreadfile) tol2e = tol2e_save
c
         if(odebug) then
            if(ga_nodeid().eq.0) then
               write(LuOut,*) 'before calling ga_pcg_minim ...'
               write(LuOut,*) 'mode = ',mode
               write(LuOut,*) 'iter, energy = ',iter,energy
            endif
         endif
         if (ga_nodeid().eq.0) call util_flush(6)
c
         call ga_pcg_minim(crohf_vlen, iter, energy, g_grad, g_work,
     $        g_search, step, ls_tol, ls_max, eprec, oconjugacy, 
     $        oprint_conv, oprint_ls, iwork, dwork, mode)
c
         if(odebug) then
            if(ga_nodeid().eq.0) then
               write(LuOut,*) 'after  calling ga_pcg_minim ...'
               write(LuOut,*) 'mode = ',mode
               write(LuOut,*) 'iter, energy = ',iter,energy
            endif
         endif
         if (ga_nodeid().eq.0) call util_flush(6)
c     
         if (mode .eq. 'energy+gradient') then
c     
c     Compute the energy and gradient at step*search direction
c     
*            write(LuOut,*) ' TOL2E ', tol2e
*            call util_flush(LuOut)
c
            if(odebug) then
               if(ga_nodeid().eq.0) then
                  write(LuOut,*) 'before calling roks_step_energy ...'
                  write(LuOut,*) 'mode = ',mode
                  write(LuOut,*) 'iter, energy = ',iter,energy
                  write(LuOut,*) 'iter, eone   = ',iter,eone  
                  write(LuOut,*) 'iter, etwo   = ',iter,etwo
                  write(LuOut,*) 'iter, enrep  = ',iter,enrep 
               endif
            endif
            if (ga_nodeid().eq.0) call util_flush(6)

            call dft_roks_step_energy(rtdb, step, g_search,
     $           eone, etwo, exc, enrep, energy, edisp, g_grad, nexc)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
c
            if(odebug) then
               if(ga_nodeid().eq.0) then
                  write(LuOut,*) 'after  calling roks_step_energy ...'
                  write(LuOut,*) 'mode = ',mode
                  write(LuOut,*) 'iter, energy = ',iter,energy
                  write(LuOut,*) 'iter, eone   = ',iter,eone  
                  write(LuOut,*) 'iter, etwo   = ',iter,etwo
                  write(LuOut,*) 'iter, enrep  = ',iter,enrep 
               endif
            endif
            if (ga_nodeid().eq.0) call util_flush(6)
c
            call ecce_print1('total energy', mt_dbl, energy, 1)
            call ecce_print1('orbital gradient norm', mt_dbl, gnorm, 1)
            eprec = max(1d-7,gnorm*0.01d0, eprec_tol) ! Avoid noise in
                                                      ! LS
            goto 10
c     
         else if (mode .eq. 'precondition') then
c     
c           Precondition the gradient direction with approx./exact
c           hessian
c     
            call dft_roks_search_precond(rtdb, g_grad, g_work)
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            goto 10
c     
         else if (mode .eq. 'accept step') then
c     
c           Apply the rotation and check for overall convergence.
c           The current fock matrices and energy correspond to the
c           updated MO coefficients.
c     
            if (util_print('dft timings', print_high))
     &         time1_rot=util_cpusec()   ! start rotation time
            call rohf_rotate_movecs(step, g_search)
            if (util_print('dft timings', print_high))
     &         time2_rot=util_cpusec()   ! end rotation time
c
c           print rotation time
c
            if(util_print('dft timings', print_high)) then
              if (ga_nodeid().eq.0) then
                write(*,"(4x,'Rotation Time:',F13.1,'s')")
     &                      time2_rot-time1_rot
              endif
            end if

            converged = (gnorm.lt.gnorm_tol)
            time_prev_iter = util_wallsec() - time_iter_start
            time_iter_start = util_wallsec()
c
            if (tol2e .ne. tol2e_save) then
c
c              Switch of variable screening will happen in the middle of
c              a line search so must ensure that we have sufficient
c              precision in the energy so that the line search does not
c              fail. Also must allow for quadratic convergence.
c     
               if (converged .or. pflg .eq. 2) then
                  tol2e = tol2e_save
               else
                  tol2e = tol2e*0.1d0
                  if (tol2e .lt. 2.0d0*tol2e_save) tol2e = tol2e_save
               endif
               converged = .false. ! Only converge with full precision
               if (util_print('parameters', print_high) .and.
     $              ga_nodeid() .eq. 0) then
                  write(LuOut,2310) tol2e
 2310             format(/' Increasing Fock-build accuracy to ',1p,d9.2)
                  call util_flush(LuOut)
               endif
            else
c
c              Reset eprec AFTER we have got a new start at the revised
c              tolerance
c
               eprec_tol = max(1d-7,100.0d0*tol2e) ! To avoid noise in
                                                   ! line search
            endif
            if ((iter.le.maxiter).and.(.not. converged)) then
               if (mod(iter,5).eq.0) 
     $              call ga_orthog_mos(basis, g_movecs)
c     
c              Not finished ... dump mo-vectors to disk for restart
c              and go back to top of loop
c     
               if (oprint_vecs) then
                  if (ga_nodeid() .eq. 0) then
                     write(LuOut,*)
                     call util_print_centered(LuOut,
     $                    'Intermediate MO vectors',40,.true.)
                     write(LuOut,*)
                     call util_flush(LuOut)
                  endif
                  call ga_print(g_movecs)
               endif
               ttest = util_wallsec()-tlastwrite
               call ga_dgop(1324, ttest, 1, 'max')  ! So all make same
                                                    ! choice
               if (ttest .gt. 300) then
                 call scf_movecs_write(rtdb)
                 tlastwrite = util_wallsec()
               endif
               if (util_test_time_remaining(rtdb,time_prev_iter*3))
     $              goto 10
            endif
c            
         endif
      endif
c     
c     End SCF minimisation
c     
      if (ga_nodeid().eq.0) call util_flush(LuOut)
c
c     ----- if sol_phase=.T. go back and do solution phase -----
c
      odebug=.false.
      odebug=odebug.and.ga_nodeid().eq.0
      if(cosmo_on.and.cosmo_phase.eq.1) then
         converged =.false.
         egas      =energy
         if (.not. rtdb_put(rtdb, 'cosmo:egas', MT_DBL, 1, egas))
     &     call errquit('dft_roks_nr_solve: failed to store egas', 0,
     &       RTDB_ERR)
c
         cosmo_phase = 2
         if(odebug) then
            write(LuOut,*) 'gas_phase done, do sol-phase now ...'
         endif
         call ga_sync()
         go to 9
      elseif(cosmo_on.and.cosmo_phase.eq.2) then
         esol      =energy
         if (.not. rtdb_put(rtdb, 'cosmo:esol', MT_DBL, 1, esol))
     &     call errquit('dft_roks_nr_solve: failed to store esol', 0,
     &       RTDB_ERR)
c
      endif
      if (ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c     
c     *ALWAYS* return canonical MOs ... don't do aufbau
c     
      call rohf_canon(.false., oprint_eval) ! Messes up symmetry of orbitals
      if (oadapt) call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     $        0, 'after canon', .true., int_mb(k_irs))
c
      call movecs_fix_phase(g_movecs)
      if (olock) call scf_movecs_lock
      call scf_movecs_write(rtdb)
      if (converged .and. nopen.gt.0 .and. olagr) call rohf_lagr()
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('dft_roks_nr_solve: ga_destroy grad', 0, GA_ERR)
      if (.not. ga_destroy(g_search)) call errquit
     $     ('dft_roks_nr_solve: ga_destroy search', 0, GA_ERR)
      if (.not. ga_destroy(g_work)) call errquit
     $     ('dft_roks_nr_solve: ga_destroy work', 0, GA_ERR)
c     
      if (ga_nodeid().eq.0) call util_flush(LuOut)
      call ga_sync()
c     
      dft_roks_nr_solve = converged
c     
      end
