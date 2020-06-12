C
C     Time-dependent absorbing "boundary" by selectively draining
C     virtual MOs via phenomenological complex absorbing potential
C     (CAP).
C
C     see: K. Lopata and N. Govind, J. Chem. Theory Comput. 9, 4939 (2013)
C
C     This routine creates the TD complex absorbing potential
C     in the orthogonal basis.
C
C     The functional form is a clamped exponential in the eigenvalues.
C
      subroutine rt_tddft_mocap (params, g_zfock_mo, g_zcap)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zfock_mo


C     == Outputs ==
      integer, intent(in)           :: g_zcap  !complex MO basis absorbing potential (nmo x nmo)

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_mocap: "


C     == Variables ==
      integer :: me
      double precision :: elapsed
      double precision :: val, val2
      double complex :: zval
      integer :: g_zcap_diag
      integer :: g_zevecs
      integer :: lvals, ivals
      double precision :: evalue, evalue_max, ecut, capmax
      double precision :: sigw, ee, sigarg
      integer :: icapval, lcapval !only for printing
      integer :: i

      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()


C
C     Allocation
C      
      if (.not. ga_duplicate (g_zcap, g_zcap_diag, "cap diag"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_duplicate (g_zcap, g_zevecs, "zevecs"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not.ma_push_get(mt_dcpl, params%ns_mo,'vals',lvals,ivals))
     $     call errquit(pname//'cannot allocate vals',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_mo,
     $     'capval', lcapval, icapval))
     $     call errquit(pname//'cannot allocate capval',0, MA_ERR)


C
C     Generate CAP in diagonal representation of instantaneous eigenstates
C     by diagonalizing TD complex density matrix.
C
C
C          [ 0   0 ]
C     G0 = [       ]
C          [ 0 -ig ]
C
      call ga_zero (g_zcap_diag)

      call ga_zero (g_zevecs)
      call zmat_diag (g_zfock_mo, g_zevecs,
     $     dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy

      
C     
C     Construct phenomenological absorbing potential, which is exponential
C     in the eigenvalues.
C     
C     This assumes "exp" is positive in input deck, such that overall
C     potential is negative imaginary.
C
C     Done on all nodes (not costly).
C
      if (params%mocap_print) then
         if (me.eq.0) then
            write(luout,*) ""
            write(luout, *) "i           Re[e_i]               ",
     $           "Im[e_i]               Re[CAP_i]             Im[CAP_i]"
            write(luout, *) "-----------------------------------------",
     $           "----------------------------------------------------"
         endif
      endif

      do i = 1, params%ns_mo
         ee = dble (dcpl_mb (ivals + i - 1))
         
         if (ee .le. params%mocap_emin) then
            val = 0d0
         else
            val = 1d0 - 
     $           exp (params%mocap_expconst * (ee - params%mocap_emin))
         endif
         
         val = val * params%mocap_prefac 

C     clamp large values
         if (val .lt. -1d0*params%mocap_maxval)
     $        val = -1d0*params%mocap_maxval
         
         dbl_mb (icapval + i - 1) = val
         zval = dcmplx (0d0, val)
         
         call ga_put (g_zcap_diag, i, i, i, i, zval, 1)

C     print to screen
         if (params%mocap_print) then
            if (me.eq.0) then
               write(6, "(1x,i6,4es22.12e3)") i, ee, 0d0, !xxx hardcoded 0 imag part eval
     $              dble(zval), aimag(zval)
            endif
         endif
         
      enddo

      if (params%mocap_print) then
         if (me.eq.0) then
            write(luout, *) ""
            write(luout, *) ""
         endif
      endif
      
      call ga_sync ()


C
C     Project onto time-dependent adiabats of the TD Fock matrix
C
C
      if (params%mocap_checks) then
         call rt_tddft_project_check (params, g_zcap_diag,
     $        g_zfock_mo, .true.) 
         if (me.eq.0) write (luout,*) pname//"projection check passed"
      endif
      
      call rt_tddft_project (params, "R",
     $     g_zcap_diag, g_zfock_mo, g_zcap, .true.)



C
C      Clean up
C      
      if (.not. ga_destroy (g_zcap_diag))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ma_chop_stack (lvals))
     $   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (params%prof)
     $     call prof_end (elapsed, "Molecular orbital CAP construction")

      end subroutine

