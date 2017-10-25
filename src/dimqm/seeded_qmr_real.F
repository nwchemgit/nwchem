      subroutine seeded_qmr_real(rtdb, nlen, g_work, niter, tol, ierr,
     $                    unknowns, matvec, xyz, ldebug, dir)
c
c     Fancy QMR algorithm.  Will document later.
c      use constants
      implicit None
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "dimqm_constants.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
c
c     Input Variables
      integer rtdb
      integer nlen
      integer g_work(6)
      integer niter
      double precision tol
      integer ierr
      integer g_unknowns
      external matvec
      double precision xyz(3,*)
      double precision unknowns(nlen)
      logical ldebug
      character dir
c
c     Local variables
      integer maxiter
      double precision d_n, e_n
      double precision sinn, cosn
      double precision rhs_n, rhs_n_plus1
      double precision l_n, l_n_plus1
      double precision res_0, res_n, res_n_minus1, res_n_plus1
      double precision omegamax, omega
      double precision res_nrm, res_nrm_upper_bound, nrm_check
      double precision nrm_tol, min_tol, max_tol
      double precision scv, scpn, scs
      double precision u_n_minus1, dtmp
      integer i, j, k
      integer id, ld, k_seed, l_seed
      double precision time
      double precision g_n
      logical stat
c
c     BLAS/LAPACK routines
      external          dlamch
      double precision  dlamch
      double precision dnrm2
      external dnrm2
      external drotg
      call ga_sync()
c
c     Get node ID
      id = ga_nodeid()
      if(id.eq.0.and.ldebug) then
        write(LuOut,*) "Start non-seed QMR routine"
        call util_flush(LuOut)
      end if
c     Check tolerances of the machine
c     Not sure if this is needed.
      nrm_tol = dlamch('E') * TEN
      min_tol = SQRT(SQRT(dlamch('S')))
      max_tol = ONE / min_tol
      if (tol <= ZERO) tol = SQRT(dlamch('E'))
      if(ldebug) then
        if(id .eq. 0) write(LuOut,*) "Tol: ", tol
      end if
c
c     Zero work arrays
      do i = 1, 6
        call ga_zero(g_work(i))
      end do
c
c     Allocate temp storage
      if(id.eq.0) then
        if(.not.ma_push_get(mt_dbl, nLen, 'seed:temp', l_seed, k_seed))
     $    call errquit('seeded temp malloc failed', 9, MA_ERR)
      end if
c
c     Initalize by copying input field to work 2 and 3
c     We do this since the inital guess to the solution is always zero
      call ga_put(g_work(2), 1, nLen, 1, 1, unknowns, 1)
      call ga_copy(g_work(2), g_work(3))
c   Save the initial residual vector in work(6)
      call ga_copy(g_work(3), g_work(6))
      call ga_sync()
c
c     Initial residual
c      res_0 = dnrm2(nlen,unknowns,1)
      res_0 = SQRT(ga_ddot(g_work(3), g_work(3)))
      if(id.eq.0.and.ldebug) write(luout,*) "res_0:", res_0
c
c     If already converged, exit here
      if((tol >= ONE) .or. (res_0 <= tol)) then
        niter = 0
        tol = res_0
        if(id.eq.0) then
          stat = ma_pop_stack(l_seed)
        end if
        return
      end if
c     Pull in how many seed vectors we have
      if(.not.rtdb_get(rtdb,'seed:k'//dir, mt_int, 1, k))
     $      k = 0
c     $      call errquit('get seed:k failed', 1, RTDB_ERR)
c
c     Now project the seed space onto the iterate and residual
      rhs_n     = res_0
      if(id.eq.0.and.ldebug) write(luout,*) "Projecting", k, "steps"
      do j = 1, k
c       Pull in scalars
        if(id.eq.0) then
          stat = rtdb_parallel(.false.)
          if(.not.rtdb_get(rtdb,'seed:sin'//CHAR(j)//dir,
     $                     mt_dbl, 1, sinn))
     $      call errquit('get seed sinn failed', j, RTDB_ERR)
          if(.not.rtdb_get(rtdb,'seed:cos'//CHAR(j)//dir,
     $                     mt_dbl, 1, cosn))
     $      call errquit('get seed cosn failed', j, RTDB_ERR)
          if(.not.rtdb_get(rtdb,'seed:scp'//CHAR(j)//dir,
     $                     mt_dbl, 1, scpn))
     $      call errquit('get seed scpn failed', j, RTDB_ERR)
          if(.not.rtdb_get(rtdb,'seed:scs'//CHAR(j)//dir,
     $                     mt_dbl, 1, scs))
     $      call errquit('get seed scsn failed', j, RTDB_ERR)
c
c       Pull in arrays
c        call ga_access(g_work(4), 1, nlen, 1, 1, k_seed, ld)
        if(.not.rtdb_get(rtdb,'seed:v'//CHAR(j)//dir,
     $                   mt_dbl, nlen, dbl_mb(k_seed)))
     $      call errquit('get seed work(3) failed', j, RTDB_ERR)
        call ga_put(g_work(4), 1, nLen, 1, 1, dbl_mb(k_seed), 1)
c        call ga_release_update(g_work(4), 1, nlen, 1, 1)
c        call ga_access(g_work(5), 1, nlen, 1, 1, k_seed, ld)
        if(.not.rtdb_get(rtdb,'seed:s'//CHAR(j)//dir,
     $                   mt_dbl, nlen, dbl_mb(k_seed)))
     $      call errquit('get seed s failed', j, RTDB_ERR)
        call ga_put(g_work(5), 1, nLen, 1, 1, dbl_mb(k_seed), 1)
c        call ga_release_update(g_work(5), 1, nlen, 1, 1)
          stat = rtdb_parallel(.true.)
        end if
        call ga_brdcst(msg_dimqm, sinn, mdtob(1), 0)
        call ga_brdcst(msg_dimqm, cosn, mdtob(1), 0)
        call ga_brdcst(msg_dimqm, scpn, mdtob(1), 0)
        call ga_brdcst(msg_dimqm, scs,  mdtob(1), 0)
        call ga_sync()
        g_n = ga_ddot(g_work(6), g_work(4)) * scpn
        rhs_n_plus1 = -sinn * rhs_n + cosn * g_n
        rhs_n       =  cosn * rhs_n + sinn * g_n
        dtmp = scs * rhs_n
        call ga_add(ONE, g_work(1), dtmp, g_work(5), g_work(1))
        rhs_n = rhs_n_plus1
      end do
      call matvec(rtdb, nlen, g_work(1), g_work(3), xyz, ldebug)
      call ga_add(ONE, g_work(2), -ONE, g_work(3), g_work(3))
      dtmp = SQRT(ga_ddot(g_work(3), g_work(3)))
c
c     Initialize the variables
      maxiter   = niter
      scv       = dtmp
      e_n       = ONE
      cosn      = ONE
      res_nrm   = dtmp/res_0
      scpn      = ONE
      scs       = ZERO
      sinn      = ZERO
      l_n_plus1 = ZERO
      omega     = ONE
      rhs_n     = omega * dtmp
      omegamax  = ONE / omega
c
c     Re-zero work arrays used in seeding
      call ga_zero(g_work(4))
      call ga_zero(g_work(5))
c
c     Begin the algorithm
      do niter = 1, maxiter
c        k = k + 1
        time = util_timer()
c
c       Check if E_n is nonsingular
        if(e_n == ZERO) then
          ierr = 4
          return
        end if
c
c       Compute scale factor for the vector w_{n}
c       Check for invariant subspaces, and scale the vectors if needed.
        ierr = 0
        if (scpn * scv < nrm_tol) ierr = 5
        if (ierr /= 0) return ! A-invarient subspace
c
        d_n = ga_ddot(g_work(3), g_work(3)) / (scv**2)
        if((scv >= max_tol) .or. (scv <= min_tol)) then
          dtmp = ONE / scv
          call ga_scale(g_work(3), dtmp)
          scv = ONE
        end if
        scv = ONE / scv
c
c       Build the vector p_n
        u_n_minus1 = d_n * l_n_plus1 / e_n
        dtmp       = u_n_minus1 * scpn / scv
        call ga_add(ONE, g_work(3), -dtmp, g_work(4), g_work(4))
        scpn = scv
c
c       Check if D_n is nonsingular
        if(d_n == ZERO) then
          ierr = 4
          return
        end if
c
c       Multiply current residual by the matrix
        call matvec(rtdb, nlen, g_work(4), g_work(6), xyz, ldebug)
c
c       Compute p_n^T A p_n
        e_n = ga_ddot(g_work(4), g_work(6)) * scpn**2
c
c       Build the vector v_{n+1}
        l_n = e_n / d_n
        call ga_add(ONE, g_work(6), -l_n, g_work(3), g_work(3))
c
c       Compute the scale factor for v_{n+1}
        scv = SQRT(ga_ddot(g_work(3), g_work(3)))
        l_n_plus1 = scpn * scv
c       
c       The QMR code starts here.
c       Multiply the new column by the previous omeags
c       Get the next scaling factor omega(i) and update omegamax
        res_n       = omega * l_n
        omega       = ONE
        res_n_plus1 = omega * l_n_plus1
        omegamax    = MAX(omegamax, ONE/omega)
c
c       Apply the previous rotation
        res_n_minus1 = sinn * res_n
        res_n        = cosn * res_n 
c
c       Compute and apply the rotation for the last element
        call drotg(res_n, res_n_plus1, cosn, sinn)
c
c       Apply the new rotation to the right-hand side vector
        rhs_n_plus1 = -sinn * rhs_n
        rhs_n       =  cosn * rhs_n
c
c       Compute the next search direction s_i
        dtmp = res_n_minus1 * scs / scpn
c       g_work(:,5) = g_work(:,4) + -dtmp * g_work(:,5)
        call ga_add(ONE, g_work(4), -dtmp, g_work(5), g_work(5))
c
c       Compute the new QMR iterate, then scale the search direction
        scs  = scpn / res_n
        dtmp = scs * rhs_n
c       g_work(:,1) = g_work(:,1) + dtmp * g_work(:,5)
        call ga_add(ONE, g_work(1), dtmp, g_work(5), g_work(1))
        if((ABS(scs) >= max_tol) .or. (ABS(scs) <= min_tol)) then
c         g_work(:,5) = scs * g_work(:,5)
          call ga_scale(g_work(5), scs)
          scs = ONE
        end if
c
c       Compute the residual norm upper bound
c       If the scaled upper bound is within one order of magnitude of
c       the targer convergence norm, compute the true residual norm.
        rhs_n = rhs_n_plus1
        res_nrm_upper_bound =
     $    SQRT(REAL(niter+1))*omegamax*ABS(rhs_n_plus1)/res_0
        nrm_check = res_nrm_upper_bound
        if((res_nrm_upper_bound/tol <= TEN).or.(niter >= maxiter)) then
c         Multiply the current residual by the matrix
           call matvec(rtdb, nlen, g_work(1), g_work(6), xyz, ldebug)
c          g_work(:,6) = g_work(:,2) - g_work(:,6)
           call ga_add(ONE, g_work(2), -ONE, g_work(6), g_work(6))
           call ga_sync()
           res_nrm = SQRT(ga_ddot(g_work(6), g_work(6))) / res_0
           call ga_sync()
           nrm_check = res_nrm
           if(id.eq.0.and.ldebug) write(LuOut,*) "Res_nrm: ", res_nrm
        else
          if(id.eq.0.and.ldebug) then
            write(LuOut,*) "Res_nrm_upper_bound:", res_nrm_upper_bound
          end if
        end if
        time = util_timer() - time
        if(id.eq.0.and.ldebug) then
          write(LuOut,*) "Iteration", niter
          write(LuOut,*) "Time (s):", time
          call util_flush(LuOut)
        end if
c
c       Check for convergece or termination.  Stop if:
c         1. algorithm converged;
c         2. there is an error condition;
c         3. the residual norm upper bound is smaller than the computed
c            residual norm by a factor of at least 100;
c         4. algorithm exceeded the iterations limit
        if(res_nrm <= tol) then
          ierr = 0
          exit
        else if(ierr /= 0) then
          exit
        else if(res_nrm_upper_bound < nrm_check / HUNDRED) then
          ierr = 3
          exit
        end if
c        call ga_get(g_work(1), 1, nLen, 1, 1, unknowns, 1)
c        write(luout,*) "END QMR DIPOLES"
c        write(luout,*) unknowns
        call ga_sync()
      end do
c
c     Put proper values into output variables
      if (niter > maxiter) ierr = 3
      tol = res_nrm
      call ga_get(g_work(1), 1, nLen, 1, 1, unknowns, 1)
      if(id.eq.0) then
        stat = ma_pop_stack(l_seed)
      end if
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "End QMR Routine"
      end subroutine seeded_qmr_real
