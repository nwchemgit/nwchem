C
C     zmat_max_abs_elem
C
C     Determines the maximum absolute element of a complex matrix.
C
      subroutine zmat_max_abs_elem (g_w, max_abs_elem)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"
      

C     == Inputs ==
      integer, intent(in) :: g_w


C     == Outputs ==
      double precision, intent(out) :: max_abs_elem


C     == Parameters ==
      character, parameter :: pname = "zmat_spectral_elem: "


C     == Variables ==
      integer dtype, n1, n2
      integer i, j
      integer ibuff, lbuff
      double precision abs_elem


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
      if (.not.ma_push_get (mt_dcpl, n1*n1, 'vals', lbuff, ibuff))
     $     call errquit (pname//'cannot allocate buff',0, MA_ERR)

      
C
C     Determine largest element (abs value) by unpacking GA into buffer
C     and looping over it.
C     
      call pack_ga2buffer_dcpl (g_w, dcpl_mb(ibuff))


      max_abs_elem = -99d0
      do i = 1, n1
         do j = 1, n1
            abs_elem = abs (dcpl_mb(ibuff + (n1-1)*i + j - 1))
            if (abs_elem > max_abs_elem) max_abs_elem = abs_elem
         enddo
      enddo


      
C
C     Make sure we found a legit max val
C
      if (max_abs_elem < 0d0)
     $     call errquit (pname//"failed to find max abs elem", 0, 0)


C
C     Clean up
C
      if (.not. ma_pop_stack (lbuff))
     $     call errquit(pname//'failed to pop buff',0, MA_ERR)

      end subroutine
c $Id$
