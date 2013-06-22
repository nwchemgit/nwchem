C
C     exp_diag.F
C
C     Uses diagonalization to exponentiate a complex operator.
C
C     e^A = X e^a X^-1 ,
C
C     where X is the eigenvector matrix and a is the eigenvalues.  Works
C     for any square complex matrix--XXX not extensively checked.
C
      subroutine exp_diag (params, g_za, g_zexpa)
      implicit none


#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
      

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_za 

      
C     == Outputs ==
      integer, intent(in) :: g_zexpa


C     == Parameters ==
      character(*), parameter   :: pname = "exp_diag: "
      
      
C     == Variables ==
      integer levals, ievals
      integer g_evecs
      integer g_evalmat
      integer g_zscr1
      integer me
      integer j
      logical ok
      integer dtype, n1, n2
      double complex zval, exp_eval
      double precision rval, ival
      integer stat
      double precision elapsed
      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid ()
      
C
C     Check the GA.
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//" not a valid GA")
      
      call ga_inquire (g_za, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")

      
C
C     Allocate.
C
      if (.not.ma_push_get(mt_dcpl, n1, "evals", levals, ievals))
     $     call errquit (pname//"failed to push evals", 0, MA_ERR)

      ok = .true.
      ok = ok .and. ga_duplicate(g_za, g_evecs, "evecs")
      ok = ok .and. ga_duplicate(g_za, g_evalmat, "evalmat")
      ok = ok .and. ga_duplicate(g_za, g_zscr1, "zscr1")
      if (.not.ok) call errquit (pname//"alloc failed", 0, GA_ERR)

      
C
C     Diagonalize.  We can't use ga_diag_std() as it is only for real
C     data.
C
      call zmat_diag (g_za, g_evecs, dcpl_mb(ievals), .false., .true.)
      
C      call zmat_diag_check (g_za, g_evecs,
C     $     dcpl_mb(ievals), params%tol_zero)
C      call zmat_diag_print_eigeninfo (g_za, g_evecs, dcpl_mb(ievals))



C
C     Construct matrix with e^evals as diagonals.  GA diagonal routines
C     might have bugs, so I'll do it all manually.  Shouldn't be an
C     issue as we only diagonalize for very small systems anyways.
C
      call ga_zero (g_evalmat)
      do j = 1, n1
         zval = dcpl_mb(ievals + j - 1)
         exp_eval = exp (zval)
         call ga_put (g_evalmat, j,j,j,j, exp_eval, 1)
      enddo

      call ga_sync ()           !XXX needed?


C
C     Find inverse of eigenvector matrix: X^-1; store as g_zscr1
C
      call zmat_inv (g_evecs, g_zscr1)
      call zmat_inv_check (g_evecs, g_zscr1)


C
C     Compute e^a X^-1, store in g_zscr1
C
      call ga_zgemm ("N", "N", n1, n1, n1,
     $     z1, g_evalmat, g_zscr1, z0, g_zscr1)
      

C     
C     Finally, e^A = X e^a X^-1.
C     
      call ga_zgemm ("N", "N", n1, n1, n1,
     $     z1, g_evecs, g_zscr1, z0, g_zexpa)

      
C
C     Clean up.
C      
      if (.not.ma_pop_stack(levals))
     $     call errquit (pname//"failed to pop evals", 0, MA_ERR)

      ok = .true.
      ok = ok .and. ga_destroy(g_evecs)
      ok = ok .and. ga_destroy(g_evalmat)
      ok = ok .and. ga_destroy(g_zscr1)
      if (.not.ok) call errquit (pname//"destroy failed", 0, GA_ERR)


      if (params%prof) then 
         call prof_end (elapsed, "Diagonalization")
      endif

      end subroutine

c $Id$
