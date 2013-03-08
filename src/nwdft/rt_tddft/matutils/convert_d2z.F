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
      
      call convert_d2z_parallel (n, g_ascal, g_bscal, g_z)

      
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
      integer np                !number of processors
      integer me                !current processor
      integer m0                !this proc gets a m x n chunk to work with
      integer m                 !proc 0 gets a chunk this big to work with
      integer istart, iend
      integer bfr_size
      integer la, ia, lb, ib    !real buffer handles and indicies
      integer lz, iz            !complex buffer handle and index
      integer jslice, j
      logical ok

      call ga_sync ()

      me = ga_nodeid ()
      np = ga_nnodes ()

C
C     Determine how much work this current processor does, which is a
C     slice of size m x n.  Since the total matrix size is n x n, each
C     processor will get an n/np x n slice.  If n/np does not divide
C     evenly, the extra will be given to processor 0.
C
      m0 = n/np + mod (n, np)
      
      if (me.eq.0) then
         m = m0
      else
         m = n/np
      endif


C      write (luout, "(1x,a,i0,a,i0,a,i0,a)")
C     $     "Processor ", me, " gets ", m, " x ", n, "."


C
C     Find the start and end index for the column index.  That is, this
C     processor is working with an m x n slice, so we need to figure out
C     what part of the full matrix we are working with.  We use fortran
C     notation (starting at 1).  Proc 0 goes from 1 to m0, proc 1 goes
C     from m0 + 1 to m0 + 1 + m, etc.
C
      if (me.eq.0) then
         istart = 1
         iend = m0
      else
         istart = m0 + 1 + (me-1)*m
         iend = istart + m - 1
      endif

      
C      write (luout, "(1x,a,i0,a,i0,a,i0,a, i0, a, i0, a)")
C     $     "Processor ", me, " starts at ", istart,
C     $     " and ends at ", iend, "; length ", m, " x ", n, "."

      
C
C     Now allocate the two real local buffers (a,b) and the one complex
C     local buffer (z).  Each buffer is size n.  
C
      bfr_size = n

      ok = .true.
      ok = ok .and. ma_push_get(mt_dbl,  bfr_size, pname//"a", la, ia)
      ok = ok .and. ma_push_get(mt_dbl,  bfr_size, pname//"b", lb, ib)
      ok = ok .and. ma_push_get(mt_dcpl, bfr_size, pname//"z", lz, iz)
      if (.not.ok) call errquit (pname//"failed to allocate", 0, MA_ERR)


      do jslice = istart, iend
         call ga_get (g_a, 1, n, jslice, jslice, dbl_mb(ia), 1)
         call ga_get (g_b, 1, n, jslice, jslice, dbl_mb(ib), 1)


C
C     For the current slice on this proc, convert and load a row of the
C     matrix.  Huub's idea to split loops seems to yield about 25%
C     speedup.
C
         do j = 1, n
            dcpl_mb(iz+j-1) = dcmplx (dbl_mb(ia+j-1), 0d0)
         enddo

         do j = 1, n
            dcpl_mb(iz+j-1) = dcpl_mb(iz+j-1) +
     $           dcmplx (0d0, dbl_mb(ib+j-1))
         enddo

         call ga_put (g_z, 1, n, jslice, jslice, dcpl_mb(iz), 1)
         
      enddo                     ! end loop over slices on this proc
      

C
C     Clean up
C
      if (.not. ma_chop_stack (la))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)

      call ga_sync ()
      
      end subroutine
