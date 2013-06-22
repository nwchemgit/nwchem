C
C     zmat_spectral_range
C
C     Determines the spectral range of a complex matrix W:
C
C     p(W) = max |eigenvalue|
C
C     This is especially useful for determining how much to scale W by
C     when doing power series exponentiation since ||W|| = p(W) (I
C     forget the name of this type of norm).
C
C
      subroutine zmat_spectral_range (g_w, max_abs_eval)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"
      

C     == Inputs ==
      integer, intent(in) :: g_w


C     == Outputs ==
      double precision, intent(out) :: max_abs_eval


C     == Parameters ==
      character, parameter :: pname = "zmat_spectral_range: "


C     == Variables ==
      integer g_evecs
      integer lvals, ivals
      integer dtype, n1, n2
      integer i
      double precision abs_eval


C
C     Check the GA.
C
      call ga_check_handle (g_w,
     $     "first argument of "//pname//" not a valid GA")
      
      call ga_inquire (g_w, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")

      
C
C     Allocation
C      
      if (.not. ga_duplicate (g_w, g_evecs, "evecs"))
     $     call errquit ("failed to create evecs", 0, 0)
      

      if (.not.ma_push_get (mt_dcpl, n1, 'vals', lvals, ivals))
     $     call errquit (pname//'cannot allocate vals',0, MA_ERR)



C
C     Diagonalize
C
      call zmat_diag (g_w, g_evecs, dcpl_mb(ivals), .false., .true.)
C      call zmat_diag_check (g_w, g_evecs, dcpl_mb(ivals), 1d-8)


C
C     Determine largest eigenvalue (abs value)
C     
      max_abs_eval = -99d0
      
      do i = 1, n1
         abs_eval = abs (dcpl_mb(ivals + i - 1))
         if (abs_eval .gt. max_abs_eval) max_abs_eval = abs_eval
      enddo


C
C     Make sure we found a legit max eval
C
      if (max_abs_eval < 0d0)
     $     call errquit (pname//"failed to find max eval", 0, 0)


C
C     Clean up
C
      if (.not. ga_destroy (g_evecs))
     $     call errquit ("failed to destroy evecs", 0, 0)

      if (.not. ma_pop_stack (lvals))
     $     call errquit(pname//'failed to pop vals',0, MA_ERR)

      end subroutine
c $Id$
