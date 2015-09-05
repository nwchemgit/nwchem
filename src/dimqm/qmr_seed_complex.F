      subroutine qmr_seed_complex(rtdb, nlen, g_work, niter, tol, ierr,
     $                            unknowns, matvec, xyz, ldebug, dir)
c
c     Fancy QMR algorithm.  Will document later.
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
c
c     Input Variables
      integer rtdb
      integer nlen
      integer g_work(7)
      integer niter
      double precision tol
      integer ierr
      integer g_unknowns
      external matvec
      double precision xyz(3,*)
      double complex unknowns(nlen)
      logical ldebug
      character dir
c
c     Local variables
      integer maxiter
      complex*16 d_n, e_n, l_n
      complex*16 sinn, scs
      complex*16 rhs_n, rhs_n_plus1
      complex*16  res_n, res_n_minus1, res_n_plus1
      complex*16 u_n_minus1, dtmp
      double precision l_n_plus1
      double precision omegamax, omega, cosn, res_0
      double precision res_nrm, res_nrm_upper_bound, nrm_check
      double precision nrm_tol, min_tol, max_tol
      double precision scv, scpn, tmp
      integer i
      integer id, ld, k_seed, l_seed
      double precision time
      logical stat
c
c     BLAS/LAPACK routines
      external          dlamch
      double precision  dlamch
      double precision dnrm2
      external dnrm2
      external zrotg
      call ga_sync()
c
c     Check for valid input
      ierr = 0
      if(nlen < 1 .or. niter < 1) then
        ierr = 2
        return
      end if
c
c     Get node ID
      id = ga_nodeid()
      if(id.eq.0.and.ldebug) then
        write(LuOut,*) "Start seed QMR complex routine"
        call util_flush(LuOut)
      end if
c     Check tolerances of the machine
c     Not sure if this is needed.
      nrm_tol = dlamch('E') * TEN
      min_tol = SQRT(SQRT(dlamch('S')))
c      write(luout,*) "mintol:", min_tol
c      min_tol = 1.0d-30
      max_tol = ONE / min_tol
      if (tol <= ZERO) tol = SQRT(dlamch('E'))
      if(ldebug) then
        if(id .eq. 0) write(LuOut,*) "Tol: ", tol
      end if
c
c     Zero work arrays
      do i = 1, 7
        call ga_zero(g_work(i))
      end do
      if(id.eq.0) then
        if(.not.ma_push_get(mt_dcpl, nlen, 'qmr:seedtmp', l_seed,
     $                      k_seed))
     $    call errquit('qmr complex malloc k_seed failed', 9, MA_ERR)
      end if
c
c     Initalize by placing input field into work 2 and 3
      call ga_put(g_work(2), 1, nLen, 1, 1, unknowns, 1)
      call ga_copy(g_work(2), g_work(3))
      call ga_sync()
c
c     Initial residual
c      res_0 = dnrm2(nlen,unknowns,1)
      call ga_conjg(g_work(3), g_work(7))
      res_0 = SQRT(REAL(ga_zdot(g_work(7), g_work(3)), KIND=8))
c
c     If already converged, exit here
      if((tol >= ONE) .or. (res_0 <= tol )) then
        niter = 0
        tol = res_0
        if(id.eq.0) then
          stat = ma_pop_stack(l_seed)
        end if
        return
      end if
c
c     Initialize the variables
      maxiter   = niter
      scv       = res_0
      e_n       = ONE_C
      cosn      = ONE
      res_nrm   = ONE
      scpn      = ONE
      scs       = ZERO_C
      sinn      = ZERO_C
      l_n_plus1 = ZERO
      omega     = ONE
      rhs_n     = omega * res_0
      omegamax  = ONE / omega
c
c     Begin the algorithm
      do niter = 1, maxiter
        time = util_timer()
c
c       Check if E_n is nonsingular
        if(ABS(e_n) == ZERO) then
          ierr = 3
          return
        end if
c
c       Compute scale factor for the vector w_{n}
c       Check for invariant subspaces, and scale the vectors if needed.
        ierr = 0
        if (scpn * scv < nrm_tol) then
          ierr = 5 ! A-invarient subspace
          return
        end if
c
        d_n = ga_zdot(g_work(3), g_work(3))
        d_n = d_n / (scv**2)
        if((scv >= max_tol) .or. (scv <= min_tol)) then
          dtmp = ONE_C / scv
          call ga_scale(g_work(3), dtmp)
          scv = ONE
        end if
        scv = ONE / scv
c
c       Build the vector p_n
        u_n_minus1 = d_n * l_n_plus1 / e_n
        dtmp       = u_n_minus1 * scpn / scv
        call ga_add(ONE_C, g_work(3), -dtmp, g_work(4), g_work(4))
        scpn = scv
c
c       Check if D_n is nonsingular
        if(ABS(d_n) == ZERO) then
          ierr = 4
          return
        end if
c
c       Multiply current residual by the matrix
        call matvec(rtdb, nlen, g_work(4), g_work(6), xyz, ldebug)
c
c       Compute p_n^T A p_n
        e_n = ga_zdot(g_work(4), g_work(6))
        e_n = e_n * scpn**2
c
c       Build the vector v_{n+1}
        l_n = e_n / d_n
        call ga_add(ONE_C, g_work(6), -l_n, g_work(3), g_work(3))
        call ga_sync()
c
c       Compute the scale factor for v_{n+1}
        call ga_conjg(g_work(3), g_work(7))
        scv = SQRT(REAL(ga_zdot(g_work(7), g_work(3)), KIND=8))
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
        call zrotg(res_n, res_n_plus1, cosn, sinn)
c
c       Save rotations and vector for seeding
        if(id.eq.0) then
          stat = rtdb_parallel(.false.)
          if(.not.rtdb_put(rtdb,'seed:sin'//CHAR(niter)//dir//'c',
     $                   mt_dcpl, 1, sinn))
     $      call errquit('put sinn failed', niter, RTDB_ERR)
          if(.not.rtdb_put(rtdb,'seed:cos'//CHAR(niter)//dir//'c',
     $                   mt_dbl, 1, cosn))
     $      call errquit('put cosn failed', niter, RTDB_ERR)
          if(.not.rtdb_put(rtdb,'seed:scp'//CHAR(niter)//dir//'c',
     $                   mt_dbl, 1, scpn))
     $      call errquit('put scpn failed', niter, RTDB_ERR)
         call ga_get(g_work(3), 1, nlen, 1, 1, dcpl_mb(k_seed), 1)
c        call ga_access(g_work(3), 1, nlen, 1, 1, k_seed, ld)
          if(.not.rtdb_put(rtdb,'seed:v'//CHAR(niter)//dir//'c',
     $                   mt_dcpl, nlen, dcpl_mb(k_seed)))
     $      call errquit('put v failed', niter, RTDB_ERR)
c        call ga_release(g_work(3), 1, nlen, 1, 1)
          stat = rtdb_parallel(.true.)
        end if
        call ga_sync()
c
c       Apply the new rotation to the right-hand side vector
        rhs_n_plus1 = -CONJG(sinn) * rhs_n
        rhs_n       =  cosn * rhs_n
c
c       Compute the next search direction s_i
        dtmp = res_n_minus1 * scs / scpn
c       g_work(:,5) = g_work(:,4) + -dtmp * g_work(:,5)
        call ga_add(ONE_C, g_work(4), -dtmp, g_work(5), g_work(5))
c
c       Compute the new QMR iterate, then scale the search direction
        scs  = scpn / res_n
        dtmp = scs * rhs_n
c
c       Save search direction for seeding
        if(id.eq.0) then
           stat = rtdb_parallel(.false.)
           call ga_get(g_work(5), 1, nlen, 1, 1, dcpl_mb(k_seed), 1)
c          call ga_access(g_work(5), 1, nlen, 1, 1, k_seed, ld)
          if(.not.rtdb_put(rtdb,'seed:s'//CHAR(niter)//dir//'c',
     $                   mt_dcpl, nlen, dcpl_mb(k_seed)))
     $      call errquit('put s failed', niter, RTDB_ERR)
c          call ga_release(g_work(5), 1, nlen, 1, 1)
          if(.not.rtdb_put(rtdb,'seed:scs'//CHAR(niter)//dir//'c',
     $                   mt_dcpl, 1, scs))
     $      call errquit('put scs failed', niter, RTDB_ERR)
          stat = rtdb_parallel(.true.)
        end if
        call ga_sync()
c       g_work(:,1) = g_work(:,1) + dtmp * g_work(:,5)
        call ga_add(ONE_C, g_work(1), dtmp, g_work(5), g_work(1))
        if((ABS(scs) >= max_tol) .or. (ABS(scs) <= min_tol)) then
c         g_work(:,5) = scs * g_work(:,5)
          call ga_scale(g_work(5), scs)
          scs = ONE_C
        end if
c
c       Compute the residual norm upper bound
c       If the scaled upper bound is within one order of magnitude of
c       the targer convergence norm, compute the true residual norm.
        rhs_n = rhs_n_plus1
        res_nrm_upper_bound =
     $    SQRT(REAL(niter+1, KIND=8))*omegamax*ABS(rhs_n_plus1)/res_0
        nrm_check = res_nrm_upper_bound
        if((res_nrm_upper_bound/tol <= TEN).or.(niter >= maxiter)) then
c         Multiply the current residual by the matrix
           call matvec(rtdb, nlen, g_work(1), g_work(6), xyz, ldebug)
           call ga_add(ONE_C, g_work(2), -ONE, g_work(6), g_work(6))
           call ga_sync()
           call ga_conjg(g_work(6), g_work(7))
           res_nrm = SQRT(REAL(ga_zdot(g_work(7), g_work(6)),KIND=8))
           res_nrm = res_nrm / res_0
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
          write(luout,*) ""
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
        call ga_sync()
      end do
c
c     Put proper values into output variables
      if (niter > maxiter) ierr = 3
      tol = res_nrm
c      call ga_copy(g_work(1), g_unknowns)
      call ga_get(g_work(1), 1, nLen, 1, 1, unknowns, 1)
c
c     Save num,ber of iterations for projection
      if(.not.rtdb_put(rtdb,'seed:k'//dir//'c', mt_int, 1, niter))
     $  call errquit('put seed:k failed', 1, RTDB_ERR)
      if(id.eq.0) then
        stat = ma_pop_stack(l_seed)
      end if
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "End seed QMR Complex Routine"
      end subroutine qmr_seed_complex
