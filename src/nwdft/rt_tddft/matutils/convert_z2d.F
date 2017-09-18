C
C     convert_z2d.F
C
C     Unpacks complex g_z into real part g_a and imaginary part g_b.
C     All GA's must be square and the same size.
C      
C     A = Re[Z], B = Im[Z] .
C
      subroutine convert_z2d (g_z, g_a, g_b)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_z               ! complex valued GA


C     == Outputs ==
      integer, intent(in) :: g_a              ! real part of GA
      integer, intent(in) :: g_b              ! imag part of GA

      
C     == Parameters ==
      character(*), parameter :: pname = "convert_z2d: "


C     == Variables ==
C      integer dim1, dim2, dtype, n !XXXASB
      integer dim1, dim2, dtype, dim1_chk, dim2_chk !XXXASB
      integer np 


C     == Get dims of GAs and check that they are correct types ==
C
C     Check the real part (input 1).
C
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"not a valid GA")
      call ga_inquire (g_z, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as first argument", 0, 0)
C      if (dim1 .ne. dim2)
C     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      dim1_chk = dim1 !XXXASB
      dim2_chk = dim2 !XXXASB

C      
C     The size of all matricies is n x n. !n1 x n2 XXXASB
C
C      n = dim1 !XXXASB

      
C     
C     Check the imaginary part (input 1), and make sure its the
C     same size as real part.
C
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"not a valid GA")
      call ga_inquire (g_a, dtype, dim1, dim2)
      if (dtype .ne. mt_dbl) call errquit (pname//
     $     "expecting real-valued GA as second argument", 0, 0)
C      if (dim1 .ne. dim2) !XXXASB
C     $     call errquit (pname//"dim1 must equal dim2 for arg 2", 0, 0) !XXXASB
C      if (dim1 .ne. n) !XXXASB
C     $     call errquit (pname//"all GAs must have same size", 0, 0) !XXXASB
      if (dim1 .ne. dim1_chk) !XXXASB
     $     call errquit (pname//"all GAs must have same size x", 0, 0) !XXXASB
      if (dim2 .ne. dim2_chk) !XXXASB
     $     call errquit (pname//"all GAs must have same size y", 0, 0) !XXXASB

C
C     Check the complex GA (output).
C
      call ga_check_handle (g_b,
     $     "third argument of "//pname//"not a valid GA")
      call ga_inquire (g_b, dtype, dim1, dim2)

      if (dtype .ne. mt_dbl) call errquit (pname//
     $     "expecting complex-valued GA as third argument", 0, 0)
C      if (dim1 .ne. dim2) !XXXASB
C     $     call errquit (pname//"dim1 must equal dim2 for arg 3", 0, 0) !XXXASB
C      if (dim1 .ne. n) !XXXASB
C     $     call errquit (pname//"all GAs must have same size", 0, 0) !XXXASB
      if (dim1 .ne. dim1_chk)
     $     call errquit (pname//"all GAs must have same size z", 0, 0)
      if (dim2 .ne. dim2_chk)
     $     call errquit (pname//"all GAs must have same size l", 0, 0)

      
C     Call the conversion routine.
C      call convert_z2d_parallel (n, g_z, g_a, g_b) !XXXASB
      call convert_z2d_parallel (dim1, dim2, g_z, g_a, g_b) !XXXASB
C      call convert_z2d_serial (n, g_z, g_a, g_b)

      end subroutine


      
C====================================================================
C
C     Driver routine for unpacking g_z into real and imag parts.  Serial
C     version, does all work on proc 0.
C
      subroutine convert_z2d_serial (n, g_z, g_a, g_b)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_z
      integer, intent(in) :: n                 ! size of g_a, g_b, g_z

      
C     == Outputs ==
      integer, intent(in) :: g_a               ! real part of GA
      integer, intent(in) :: g_b               ! imag part of GA


C     == Parameters ==
      character(*), parameter :: pname = "convert_z2d_serial: "

C     == Variables ==
      integer me
      integer i, j
      double precision val_a, val_b
      double complex val_z
      
      me = ga_nodeid ()

      if (me.eq.0) then
         do i = 1, n
            do j = 1, n
               call ga_get (g_z, i, i, j, j, val_z, 1)
               
               val_a = dble (val_z)
               val_b = aimag (val_z)
               call ga_put (g_a, i, i, j, j, val_a, 1)
               call ga_put (g_b, i, i, j, j, val_b, 1)

            enddo
         enddo
      endif
      call ga_sync ()

      
      end subroutine


C#if 0
CC====================================================================
CC
CC     Driver routine for unpacking g_z into real and imag parts.
CC
C      subroutine convert_z2d_parallel (n, g_z, g_a, g_b)
C      implicit none
C
C#include "errquit.fh"
C#include "mafdecls.fh"
C#include "global.fh"
C#include "stdio.fh"
C
CC     == Inputs ==
C      integer, intent(in) :: g_z
C      integer, intent(in) :: n                 ! size of g_a, g_b, g_z
C
C      
CC     == Outputs ==
C      integer, intent(in) :: g_a               ! real part of GA
C      integer, intent(in) :: g_b               ! imag part of GA
C
C
CC     == Parameters ==
C      character(*), parameter :: pname = "convert_z2d_parallel: "
C
C 
CC     == Variables ==
C      integer m                 !proc 0 gets a chunk this big to work with
C      integer istart, iend
C      integer bfr_size
C      integer la, ia, lb, ib    !real buffer handles and indicies
C      integer lz, iz            !complex buffer handle and index
C      integer j
C      logical ok
C      double complex zval
C
C      
CC     not needed I think
CC      call ga_sync ()
C
C      call calc_array_partitioning (n, m, istart, iend, bfr_size)
C
C      if (bfr_size .gt. 0) then
C         ok = .true.
C         ok = ok.and.ma_push_get(mt_dbl,  bfr_size, pname//"a", la, ia)
C         ok = ok.and.ma_push_get(mt_dbl,  bfr_size, pname//"b", lb, ib)
C         ok = ok.and.ma_push_get(mt_dcpl, bfr_size, pname//"z", lz, iz)
C         if (.not.ok) call errquit(pname//"failed to allocate",0,MA_ERR)
C
C         call ga_get (g_z, 1, n, istart, iend, dcpl_mb(iz), n)
C      
C         do j = 1, bfr_size
C            zval = dcpl_mb(iz + j - 1)
C            dbl_mb(ia + j - 1) = dble (zval)
C            dbl_mb(ib + j - 1) = aimag (zval)
C         enddo
C         
C         call ga_put (g_a, 1, n, istart, iend, dbl_mb(ia), n)
C         call ga_put (g_b, 1, n, istart, iend, dbl_mb(ib), n)
C         
C         if (.not. ma_chop_stack (la))
C     $        call errquit (pname//"failed to chop stack", 0, MA_ERR)
C      endif
C
CC     not needed I think
CC      call ga_sync ()
C      
C      end subroutine
C#endif
C

      

C====================================================================
C
C     Driver routine for unpacking g_z into real and imag parts.
C
C      subroutine convert_z2d_parallel (n, g_z, g_a, g_b) !XXXASB
      subroutine convert_z2d_parallel (n1, n2, g_z, g_a, g_b) !XXXASB
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_z
C      integer, intent(in) :: n                 ! size of g_a, g_b, g_z
      integer, intent(in) :: n1, n2            ! size of g_a, g_b, g_z: n1 x n2

      
C     == Outputs ==
      integer, intent(in) :: g_a               ! real part of GA
      integer, intent(in) :: g_b               ! imag part of GA


C     == Parameters ==
      character(*), parameter :: pname = "convert_z2d_parallel: "

 
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
C     slice of size m x n.  Since the total matrix size is n1 x n2, where
C     n2 may be smaller than n1 (linear dependencies), each
C     processor will get an n2/np x n slice.  If n2/np does not divide
C     evenly, the extra will be given to processor 0.
C
C      m0 = n/np + mod (n, np) !XXXASB
      m0 = n2/np + mod (n2, np) !XXXASB
      
      if (me.eq.0) then
         m = m0
      else
C         m = n/np !XXXASB
         m = n2/np !XXXASB
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
C      bfr_size = n !XXXASB
      bfr_size = n1 !XXXASB

      ok = .true.
      ok = ok .and. ma_push_get(mt_dbl,  bfr_size, pname//"a", la, ia)
      ok = ok .and. ma_push_get(mt_dbl,  bfr_size, pname//"b", lb, ib)
      ok = ok .and. ma_push_get(mt_dcpl, bfr_size, pname//"z", lz, iz)
      if (.not.ok) call errquit (pname//"failed to allocate", 0, MA_ERR)

      
C
C     See notes in convert_d2z ().
C
      do jslice = istart, iend
C         call ga_get (g_z, 1, n, jslice, jslice, dcpl_mb(iz), n) !XXXASB
         call ga_get (g_z, 1, n1, jslice, jslice, dcpl_mb(iz), 1) !XXXASB

         do j = 1, n1 !XXXASB n to n1
            dbl_mb(ia+j-1) = dble (dcpl_mb(iz+j-1))
         enddo
         call ga_put (g_a, 1, n1, jslice, jslice, dbl_mb(ia), 1) !XXXASB

         do j = 1, n1 !XXXASB n to n1
            dbl_mb(ib+j-1) = aimag (dcpl_mb(iz+j-1))
         enddo
         call ga_put (g_b, 1, n1, jslice, jslice, dbl_mb(ib), 1) !XXXASB
      enddo

      
C
C     Clean up
C
      if (.not. ma_chop_stack (la))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)

      call ga_sync ()
      
      end subroutine
c $Id$
