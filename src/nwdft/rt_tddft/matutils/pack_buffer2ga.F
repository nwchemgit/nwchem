C
C     Pack real-valued input buffer into a GA.  Determines size of
C     buffer from the size of the GA (must be square matrix).
C
      subroutine pack_buffer2ga_dbl (buff, g_a)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      double precision, intent(in) :: buff(*)

C     == Outputs ==
      integer, intent(in) :: g_a
      
C     == Parameters ==
      character(*), parameter :: pname = "pack_buffer2ga_dbl: "

      
C     == Variables ==
      logical ok
      integer dtype, dim1, dim2
      integer n
      integer i, j
      integer itmp, ltmp

C
C     Check the output GA (input 2).
C
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"not a valid GA")
      call ga_inquire (g_a, dtype, dim1, dim2)

      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      n = dim1
      
      if (dtype .ne. mt_dbl)
     $     call errquit (pname//"expecting real GA", 0, 0)


C
C     Allocation
C
      ok = .true.
      ok = ok.and.ma_push_get(mt_dbl, n, "tmp buffer", ltmp,itmp)
      if (.not.ok) call errquit (pname//"alloc failed", 0, MA_ERR)


C
C     Pack into GA. Note pointer vs fortran indexing.
C
      do i = 1, n
         do j = 1, n
            dbl_mb(itmp + j - 1) = buff(n*(i-1) + j)
         enddo
         call ga_put (g_a, 1, n, i, i, dbl_mb(itmp), n)
      enddo


C
C     Clean up
C
      if (.not.ma_pop_stack(ltmp))
     $     call errquit (pname//"failed to pop stack", 0, MA_ERR)
      
      end subroutine 




C====================================================================
C
C     Pack complex-valued input buffer into a GA.  Determines size of
C     buffer from the size of the GA (must be square matrix).
C
      subroutine pack_buffer2ga_dcpl (buff, g_a)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      double complex, intent(in) :: buff(*)

      
C     == Outputs ==
      integer, intent(in) :: g_a

      
C     == Parameters ==
      character(*), parameter :: pname = "pack_buffer2ga_dcpl: "

      
C     == Variables ==
      logical ok
      integer dtype, dim1, dim2
      integer n
      integer i, j
      integer itmp, ltmp

C
C     Check the output GA (input 2).
C
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"not a valid GA")
      call ga_inquire (g_a, dtype, dim1, dim2)

      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      n = dim1
      
      if (dtype .ne. mt_dcpl)
     $     call errquit (pname//"expecting complex GA", 0, 0)


C
C     Allocation
C
      ok = .true.
      ok = ok.and.ma_push_get(mt_dcpl, n, "tmp buffer", ltmp,itmp)
      if (.not.ok) call errquit (pname//"alloc failed", 0, MA_ERR)


C
C     Pack into GA. Note pointer vs fortran indexing.
C
      do i = 1, n
         do j = 1, n
            dcpl_mb(itmp + j - 1) = buff(n*(i-1) + j)
         enddo
         call ga_put (g_a, 1, n, i, i, dcpl_mb(itmp), n)
      enddo


C
C     Clean up
C
      if (.not.ma_pop_stack(ltmp))
     $     call errquit (pname//"failed to pop stack", 0, MA_ERR)
      
      end subroutine 


c $Id$
