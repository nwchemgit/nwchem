C
C     prop_magnus_ord2_step.F
C
C     Steps dens mat forward dt in time (in MO basis) using second order
C     (single term) Magnus propagator.
C
C     P'(t+dt) = e^W P'(t) (e^W)^+  (even if non-Hermitian)
C              = e^W P'(t) e^-W  (assuming W is Hermitian, which is not true if we have a complex absorbing potential)
C
C     where W = -i F'(t+dt/2).
C
C     To exponentiate, we can diagonalize, use a power series, or use
C     the Baker-Campbell-Hausdorff expansion to evaluate the whole
C     expression via commutators (note only for Hermitian W).  P', F'
C     and W are all ns_mo x ns_mo.  Note that dt can be different from
C     the simulation time step!
C
      subroutine prop_magnus_ord2_step (params, dt,
     $     g_zfock_mo_mid, g_zdens_mo)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
      
      
C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      double precision, intent(in)   :: dt              ! time step dt
      integer, intent(in)            :: g_zfock_mo_mid  ! F'(t+dt/2), in MO basis

      
C     == In/out ==
      integer, intent(in) :: g_zdens_mo           ! in = P'(t), out = P'(t+dt)


C     == Parameters ==
      character(*), parameter   :: pname = "prop_magnus_ord2_step: "
C      integer, parameter        :: diag_max = 64 ! max sys size to diag

      
C     == Variables ==
      integer g_zscr1, g_zscr2  ! scratch
      integer g_w               ! W = Omega_1 = -i F'(t+dt/2)
      integer g_expw            ! e^W
      logical ok
      integer me
      double complex zdt
      double precision spec_range, max_elem
      logical imag_time


C     (junk)
      integer ivals, lvals
      integer g_vecs
      double complex zval

      
      me = ga_nodeid ()


      if (.not.rtdb_get(params%rtdb, "rt_tddft:imag_time",
     $     mt_log, 1, imag_time))
     $     imag_time = .false.

      if (imag_time)
     $     call errquit ("Imaginary time propagation disabled",0,0)

C
C     Allocation.
C
      ok = .true.
      ok = ok .and. ga_duplicate(g_zdens_mo, g_w, "W_1")
      ok = ok .and. ga_duplicate(g_zdens_mo, g_zscr1, "zscr1")
      ok = ok .and. ga_duplicate(g_zdens_mo, g_zscr2, "zscr2")
      if (.not. ok)
     $     call errquit (pname//"failed to allocate GAs", 0, GA_ERR)

         
C
C      W = -i F'(t+dt/2) dt.
C     

C     Real time propagation:
      zdt = dcmplx (dt, 0d0)

C     Imaginary time propagation:
      if (imag_time) then
         zdt = zni * zdt
         call rt_tddft_print_notice ("Doing imaginary time propagation")
      endif

      call ga_copy (g_zfock_mo_mid, g_w)
      call ga_scale (g_w, zni*zdt)


      
C     Correct way:
c$$$      zdt = dcmplx (dt, 0d0)
c$$$      call ga_copy (g_zfock_mo_mid, g_w)
c$$$      call ga_scale (g_w, zni*zdt)
C     IMAG TIME PROP
c$$$      zdt = dcmplx (-1d0*dt, 0d0) !XXX NEGATIVE?
c$$$      call ga_copy (g_zfock_mo_mid, g_w)
c$$$      call ga_scale (g_w, zdt)
c$$$      if (me.eq.0) then
c$$$         write (luout, *) pname//"IMAG TIME!"
c$$$      endif


C     (BCH exp disabled since it assumes e^W A e^-W, which might not be
C     true) if e^W is not Hermitian.
      if (params%exp_method.eq.3)
     $     call errquit (pname//"BCH exp disabled",0,0)

      
      g_expw = g_zscr2
      
      if (params%exp_method.eq.1) then
         call exp_pseries (params, g_w, g_expw)
      elseif (params%exp_method.eq.2) then
         call exp_diag (params, g_w, g_expw)
      else
         call errquit (pname//"bad exp_method", 0, 0)
      endif


      
C
C     Compute eigenvalues of e^W (for diagnostics only)
C     
#if 0
      ok = .true.
      ok = ok .and. ga_duplicate(g_zdens_mo, g_vecs, "vecs")
      ok = ok .and. ma_push_get(mt_dcpl, params%ns_mo,
     $     "evals", lvals, ivals)
      if (.not. ok) call errquit (pname//"alloc failed", 0, 0)
      
      call zmat_diag (g_expw, g_vecs, dcpl_mb(ivals), .false., .true.)
C      call zmat_diag_check (g_expw, g_vecs,
C     $     dcpl_mb(ivals), params%tol_zero)
C     call zmat_diag_print_eigeninfo (g_expw, g_vecs, dcpl_mb(ivals))
C     call zmat_diag_print_eigenvals (g_expw, g_vecs, dcpl_mb(ivals))
      
      ok = .true.
      ok = ok .and. ga_destroy (g_vecs)
      ok = ok .and. ma_pop_stack (lvals)
      if (.not. ok) call errquit (pname//"free failed", 0, 0)
#endif



C
C     Since F is Hermitian, W = -iF is skew Hermitian, and thus e^W must
C     be unitary.  Note, this will throw an error if doing non-Hermitian
C     propagation (e.g., using a CAP).
C
      if (params%checklvl .ge. 2) then
         if (.not. mat_is_hermitian (g_zfock_mo_mid,
     $        params%tol_zero))
     $        call errquit (pname//"fock mid not hermitian", 0, 0)
         
         if (.not. mat_is_unitary (g_expw, params%tol_zero))
     $        call errquit (pname//"e^W not unitary", 0, 0)
      endif


C        
C     P'(t+dt) = e^W P'(t) (e^W)^+.  Recycle g_w as scratch space to
C     hold (e^W)^+.
C
      call ga_zero (g_w)
      call mat_hermconj (g_expw, g_w)
      
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens_mo, g_w, z0, g_zscr1)
      
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_expw, g_zscr1, z0, g_zdens_mo)


C     
C     Clean up.
C     
      ok = .true.
      ok = ok .and. ga_destroy(g_w)
      ok = ok .and. ga_destroy(g_zscr1)
      ok = ok .and. ga_destroy(g_zscr2)
      
      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs", 0, GA_ERR)


      end subroutine
c $Id$
