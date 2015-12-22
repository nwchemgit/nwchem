C
C     dmat_io.F
C
C     Input/output routines for dumping and reading double precision
C     matrices to/from file.
C
      
C====================================================================
C
C     Dump matrix to file.
C
      logical function dmat_io_dump (g_a, fname)
      
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in)      :: g_a   ! matrix to output
      character(*), intent(in) :: fname ! file name to output to


C     == Parameters ==
      character(*), parameter :: pname = "dmat_io_dump: "
      integer, parameter      :: unitno = 1893747    !XXX TERRIBLE: HARDCODED UNIT #

      
C     == Variables ==
      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios

      dmat_io_dump = .true.


C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (dtype .ne. mt_dbl) call errquit (
     $     pname//" only valid for double precision matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")
      


C
C     Do all work on processor 0.
C     
      if (ga_nodeid() .eq. 0) then

C
C     Allocate memory for the mat slice.
C     
         if (.not. ma_push_get (mt_dbl, n1, "mat slice",
     $        l_matsl, k_matsl))
     $        call errquit (pname//"failed to alloc mat slice", 0, 0)


C
C     Open output file.
C         
         open (unitno, status='replace', form='unformatted',
     $        file=fname, iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname)
            dmat_io_dump = .false.
            return
         endif


C
C     Loop over slices of mat and print to file.
C         
         do i = 1, n1
            call ga_get (g_a, 1, n1, i, i, dbl_mb(k_matsl), n1)
            call swrite (unitno, dbl_mb(k_matsl), n1)
         enddo

         
C
C     Close output file and deallocate.
C
         close (unitno, iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to close: "//trim(fname)
            dmat_io_dump = .false.
            return
         endif

         if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)
      endif
      call ga_sync ()

      end function


C====================================================================
C
C     Read matrix from file.
C
      logical function dmat_io_read (g_a, fname)
      
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in)      :: g_a   ! matrix to read in
      character(*), intent(in) :: fname ! file name to read from


C     == Parameters ==
      character(*), parameter :: pname = "dmat_io_read: "
      integer, parameter      :: unitno = 1893747    !XXX TERRIBLE: HARDCODED UNIT #

      
C     == Variables ==
      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios

      dmat_io_read = .true.


C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (dtype .ne. mt_dbl) call errquit (
     $     pname//" only valid for double precision  matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")
      


C
C     Do all work on processor 0.
C     
      if (ga_nodeid() .eq. 0) then

C
C     Allocate memory for the mat slice.
C     
         if (.not. ma_push_get (mt_dbl, n1, "mat slice",
     $        l_matsl, k_matsl))
     $        call errquit (pname//"failed to alloc mat slice", 0, 0)


C
C     Open output file.
C         
         open (unitno, status='old', form='unformatted',
     $        file=fname, iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname)
            dmat_io_read = .false.
            return
         endif


C
C     Loop over slices of mat and store in GA.
C         
         do i = 1, n1
            call sread (unitno, dbl_mb(k_matsl), n1)
            call ga_put (g_a, 1, n1, i, i, dbl_mb(k_matsl), n1)
         enddo

         
C
C     Close output file and deallocate.
C
         close (unitno, iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to close: "//trim(fname)
            dmat_io_read = .false.
            return
         endif

         if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)
      endif
      call ga_sync ()

      end function
c $Id$
