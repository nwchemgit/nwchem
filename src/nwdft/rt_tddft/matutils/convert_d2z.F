C====================================================================
C
C     convert_d2z.F
C     
C     Packs real part g_a (scaled by aa) and imaginary part g_b (scaled
C     by bb) into complex g_z.  All GA's must be square and the same
C     size.
C      
C     Z = a*A + i*b*B .
C
      subroutine convert_d2z (aa, g_a, bb, g_b, g_z)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a     ! real part of GA
      integer, intent(in)          :: g_b     ! imag part of GA
      double precision, intent(in) :: aa, bb  ! prefactors for real and im parts

      
C     == Outputs ==
      integer, intent(in) :: g_z               ! complex valued GA


C     == Parameters ==
      character(*), parameter :: pname   = "convert_d2z: "

      
C     == Variables ==
      integer g_ascal, g_bscal
      integer dim1, dim2, dtype, n

      
C     == Get dims of GAs and check that they are correct types ==
C
C     Check the real part (input 1).
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"not a valid GA")
      call ga_inquire (g_a, dtype, dim1, dim2)
      if (dtype .ne. mt_dbl) call errquit (pname//
     $     "expecting real-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2")

      
C      
C     The size of all matricies is n x n.
C
      n = dim1

      
C     
C     Check the imaginary part (input 1), and make sure its the
C     same size as real part.
C
      call ga_check_handle (g_b,
     $     "second argument of "//pname//"not a valid GA")
      call ga_inquire (g_b, dtype, dim1, dim2)
      if (dtype .ne. mt_dbl) call errquit (pname//
     $     "expecting real-valued GA as fourth argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2 for arg 2")
      if (dim1 .ne. n)
     $     call errquit (pname//"all GAs must have same size", 0, 0)

      
C
C     Check the complex GA (output).
C
      call ga_check_handle (g_z,
     $     "third argument of "//pname//"not a valid GA")
      call ga_inquire (g_z, dtype, dim1, dim2)

      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as fifth argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2 for arg 3")
      if (dim1 .ne. n)
     $     call errquit (pname//"all GAs must have same size", 0, 0)


C
C     Allocate some work space for scaling the matricies.
C
      if (.not.ga_duplicate(g_a, g_ascal, "A scaled"))
     $     call errquit (pname//"alloc ascal failed", 0, GA_ERR)

      if (.not.ga_duplicate(g_b, g_bscal, "B scaled"))
     $     call errquit (pname//"alloc bscal failed", 0, GA_ERR)


C
C     Scale copies of input A and B by aa and bb.
C
      call ga_copy (g_a, g_ascal)
      call ga_scale (g_ascal, aa)

      call ga_copy (g_b, g_bscal)
      call ga_scale (g_bscal, bb)

      
C
C     All GAs are the right data type, shape, and size.  Call the
C     parallel conversion routine.
C
      
      call convert_d2z_parallel (n, g_ascal, g_bscal, g_z)
C      call convert_d2z_serial (n, g_ascal, g_bscal, g_z)

      
C
C     Clean up.
C
      if (.not.ga_destroy(g_ascal))
     $     call errquit (pname//"destroy ascal failed", 0, GA_ERR)
      if (.not.ga_destroy(g_bscal))
     $     call errquit (pname//"destroy bscal failed", 0, GA_ERR)

      end subroutine




C====================================================================
C
C     Driver routine for packing g_z <- (g_a, g_b), note g_a and g_b
C     were previously scaled.  Does all work on processor 0.
C
      subroutine convert_d2z_serial (n, g_a, g_b, g_z)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_a               ! real part of GA, already scaled
      integer, intent(in) :: g_b               ! imag part of GA, already scaled
      integer, intent(in) :: n                 ! size of g_a, g_b, g_z

      
C     == Outputs ==
      integer, intent(in) :: g_z


C     == Parameters ==
      character(*), parameter :: pname = "convert_d2z_serial: "


C     == Variables ==
      integer me
      integer i, j
      double precision val_a, val_b
      double complex val_z
      
      me = ga_nodeid ()

      if (me.eq.0) then
         do i = 1, n
            do j = 1, n
               call ga_get (g_a, i, i, j, j, val_a, 1)
               call ga_get (g_b, i, i, j, j, val_b, 1)
               
               val_z = dcmplx (val_a, val_b)
               call ga_put (g_z, i, i, j, j, val_z, 1)
            enddo
         enddo
      endif
      call ga_sync ()
      
      end subroutine




C====================================================================
C
C     Driver routine for packing g_z <- (g_a, g_b), note g_a and g_b
C     were previously scaled. 
C
      subroutine convert_d2z_parallel (n, g_a, g_b, g_z)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"


C     == Inputs ==
      integer, intent(in) :: g_a               ! real part of GA, already scaled
      integer, intent(in) :: g_b               ! imag part of GA, already scaled
      integer, intent(in) :: n                 ! size of g_a, g_b, g_z

      
C     == Outputs ==
      integer, intent(in) :: g_z


C     == Parameters ==
      character(*), parameter :: pname = "convert_d2z_parallel: "

 
C     == Variables ==
      integer m
      integer istart, iend
      integer bfr_size
      integer la, ia, lb, ib    !real buffer handles and indicies
      integer lz, iz            !complex buffer handle and index
      logical ok
      integer j
      double precision aval, bval

 
C     not needed I think
C      call ga_sync ()

      call calc_array_partitioning (n, m, istart, iend, bfr_size)
      
      if (bfr_size .gt. 0) then
         ok = .true.
         ok = ok.and.ma_push_get(mt_dbl,  bfr_size, pname//"a", la, ia)
         ok = ok.and.ma_push_get(mt_dbl,  bfr_size, pname//"b", lb, ib)
         ok = ok.and.ma_push_get(mt_dcpl, bfr_size, pname//"z", lz, iz)
         if (.not.ok) call errquit(pname//"failed to allocate",0,MA_ERR)
         
         call ga_get (g_a, 1, n, istart, iend, dbl_mb(ia), n)
         call ga_get (g_b, 1, n, istart, iend, dbl_mb(ib), n)
         
         do j = 1, bfr_size
            aval = dbl_mb(ia+j-1)
            bval = dbl_mb(ib+j-1)
            dcpl_mb(iz + j - 1) = dcmplx (aval, bval)
         enddo
         
         call ga_put (g_z, 1, n, istart, iend, dcpl_mb(iz), n)

         if (.not. ma_chop_stack (la))
     $        call errquit (pname//"failed to chop stack", 0, MA_ERR)
      endif

C     not needed I think
C      call ga_sync ()

      end subroutine
