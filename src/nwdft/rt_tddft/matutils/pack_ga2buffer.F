C
C     Pack real-valued input GA into buffer.  Determines size of buffer
C     from the size of the GA (must be square matrix).
C
      subroutine pack_ga2buffer_dbl (g_a, buff)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_a

C     == Outputs ==
      double precision, intent(out) :: buff(*)
      
C     == Parameters ==
      character(*), parameter :: pname = "pack_ga2buffer_dbl: "

      
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
     $     "first argument of "//pname//"not a valid GA")
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
C     Pack into buffer. Note pointer vs fortran indexing.
C
      do i = 1, n
         call ga_get (g_a, 1, n, i, i, dbl_mb(itmp), n)
         do j = 1, n
            buff(n*(i-1) + j) = dbl_mb(itmp + j - 1)
         enddo
      enddo


C
C     Clean up
C
      if (.not.ma_pop_stack(ltmp))
     $     call errquit (pname//"failed to pop stack", 0, MA_ERR)
      
      end subroutine 



C====================================================================
C
C     Pack complex-valued input GA into buffer.  Determines size of buffer
C     from the size of the GA (must be square matrix).
C
      subroutine pack_ga2buffer_dcpl (g_a, buff)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_a

C     == Outputs ==
      double complex, intent(out) :: buff(*)
      
C     == Parameters ==
      character(*), parameter :: pname = "pack_ga2buffer_dcpl: "

      
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
     $     "first argument of "//pname//"not a valid GA")
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
C     Pack into buffer. Note pointer vs fortran indexing.
C
      do i = 1, n
         call ga_get (g_a, 1, n, i, i, dcpl_mb(itmp), n)
         do j = 1, n
            buff(n*(i-1) + j) = dcpl_mb(itmp + j - 1)
         enddo
      enddo


C
C     Clean up
C
      if (.not.ma_pop_stack(ltmp))
     $     call errquit (pname//"failed to pop stack", 0, MA_ERR)
      
      end subroutine 




C
C     example driver
C
c$$$      if (.not. ga_create(mt_dbl,3,3,"test" , 0, 0, g_test))
c$$$     $     call errquit ("failed to create test", 0, 0)
c$$$
c$$$      call ga_zero (g_test)
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            junk(i,j) = dble(i) + 10.0*dble(j)
c$$$         enddo
c$$$      enddo
c$$$      call pack_buffer2ga_dbl (junk, g_test)
c$$$      call ga_print (g_test)
c$$$
c$$$
c$$$      call ga_zero (g_zmat)
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            zjunk(i,j) = dcmplx (dble(i), 10.0*dble(j))
c$$$         enddo
c$$$      enddo
c$$$      call pack_buffer2ga_dcpl (zjunk, g_zmat)
c$$$      call ga_print (g_zmat)
c$$$      
c$$$
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            junk(i,j) = -99d0
c$$$         enddo
c$$$      enddo
c$$$      call pack_ga2buffer_dbl (g_test, junk)
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            write (6, *) i, j, junk(i,j)
c$$$         enddo
c$$$      enddo
c$$$
c$$$
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            zjunk(i,j) = dcmplx (-99d0, -88d0)
c$$$         enddo
c$$$      enddo
c$$$      call pack_ga2buffer_dcpl (g_zmat, zjunk)
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            write (6, *) i, j, zjunk(i,j)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      call errquit ("XXX", 0, 0)
c$$$

c $Id$
