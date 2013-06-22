C
C     dmat_to_file.F
C
C     Input/output routines for dumping and reading double precision
C     matrices to/from file.
C
C     Write matrix to file.
C
      logical function dmat_to_file(g_a, fname)
c      
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in)      :: g_a   ! matrix to output
      character(*), intent(in) :: fname        ! file name to output to
c
C     == Parameters ==
      character(*), parameter :: pname = "dmat_to_file: "
      integer, parameter      :: unitno = 1893747    !XXX TERRIBLE: HARDCODED UNIT #
c
C     == Variables ==
      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios
c
      dmat_to_file = .true.
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
         open (unitno, status='unknown', form='unformatted',
     $        file=fname, iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname)
            dmat_to_file = .false.
            return
         endif
c
C     Loop over slices of mat and print to file.
C         
         do i = 1, n1
            call ga_get (g_a, 1, n1, i, i, dbl_mb(k_matsl), 1)
            call swrite (unitno, dbl_mb(k_matsl), n1)
         enddo
C
C     Close output file and deallocate.
C
      close (unitno, iostat=ios)
c
      if (ios .ne. 0) then
            write (luout, *) pname//"failed to close: "//trim(fname)
            dmat_to_file = .false.
            return
      endif
c
      if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)
c
      endif  ! ga_nodeid()
c
      call ga_sync ()
c
      end function
c
C====================================================================
C
C     Read matrix from file.
C
      logical function dmat_from_file(g_a, fname)
      
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in)      :: g_a(*)   ! matrix to read in
      character(*), intent(in) :: fname ! file name to read from
C
C     == Parameters ==
      character(*), parameter :: pname = "dmat_from_file: "
      integer, parameter      :: unitno = 1893747    !XXX TERRIBLE: HARDCODED UNIT #
c
C     == Variables ==
      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios

      dmat_from_file = .true.
C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"not a valid GA")
C      
      call ga_inquire (g_a, dtype, n1, n2)
C      
      if (dtype .ne. mt_dbl) call errquit (
     $     pname//" only valid for double precision  matricies", 0, 0)
C      
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
c
      if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname)
            dmat_from_file = .false.
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
            dmat_from_file = .false.
            return
      endif
C
      if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)
c
      endif ! ga_nodeid()
c
      call ga_sync ()
C
      end function
C
C====================================================================
c
c     Map supplied index to a character string
c
      character*7 function id_to_char(idx)
      integer idx
c
      integer itmp,j0,j1,j2,j3,j4,j5
      character*7 name
c
      itmp = idx
c
      j5 = itmp/100000
      itmp = itmp - j5*100000
      j4 = itmp/10000
      itmp = itmp - j4*10000
      j3 = itmp/1000
      itmp = itmp - j3*1000
      j2 = itmp/100
      itmp = itmp - j2*100
      j1 = itmp/10
      itmp = itmp - j1*10
      j0 = itmp/1
C
      name(1:1) = '_'
      name(2:2) = CHAR(j5+ICHAR('0'))
      name(3:3) = CHAR(j4+ICHAR('0'))
      name(4:4) = CHAR(j3+ICHAR('0'))
      name(5:5) = CHAR(j2+ICHAR('0'))
      name(6:6) = CHAR(j1+ICHAR('0'))
      name(7:7) = CHAR(j0+ICHAR('0'))
c
      id_to_char = name
c
      return
      end
c
C====================================================================
c
c     Create a unique file name using processor indices
c
      subroutine createuniquefilename(filename,idx)

      implicit none

      character*(*) filename
      integer idx
      character*7 idxstring

      character*7 id_to_char
      external id_to_char

      idxstring= id_to_char(idx)
      filename=filename(1:index(filename,' ')-1)//idxstring

      return
      end
c $Id$
