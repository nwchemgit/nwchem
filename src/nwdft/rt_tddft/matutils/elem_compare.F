C
C     elem_compare.F
C
C     Compare two real-valued GA's element by element, and print differing quantities.
C
C
      subroutine elem_compare (g_a, g_b, tol)

      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Inputs ==
      integer, intent(in)                   :: g_a, g_b
      double precision, intent(in)          :: tol


C     == Parameters ==
      character(len=20), parameter :: pname = "elem_compare: "
      

C     == Variables ==
      integer me
      integer dtypea, n1a, n2a
      integer dtypeb, n1b, n2b
      integer i, j
      double complex va, vb, diff
      character*100 namea, nameb


      me = ga_nodeid()
      

C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"() is not a valid GA")

      call ga_check_handle (g_b,
     $     "second argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtypea, n1a, n2a)
      call ga_inquire (g_b, dtypeb, n1b, n2b)

      call ga_inquire_name (g_a, namea)
      call ga_inquire_name (g_b, nameb)

      if (dtypea .ne. mt_dcpl) call errquit (
     $     pname//"only valid for complex matricies", 0, 0)

      if (dtypeb .ne. mt_dcpl) call errquit (
     $     pname//"only valid for complex matricies", 0, 0)
      
      if (n1a .ne. n1b) call errquit (
     $     pname//"mats must be same size",0,0)

      if (n2a .ne. n2b) call errquit (
     $     pname//"mats must be same size",0,0)

      
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "   === Element-wise matrix comparison ===   "
         write (luout, *) ""
         write (luout, *) "Matrix A: ",trim(namea)
         write (luout, *) "Matrix B: ",trim(nameb)
         write (luout, *) ""
         write (luout, *)
     $        "    i      j     ",
     $        "     Re[A]       Re[B]            Im[A]       Im[B]     "
         write (luout, *)
     $        "------------------------------------"//
     $        "------------------------------------"         
      endif


      do i = 1, n1a 
         do j = 1, n2a
            call ga_get (g_a, i, i, j, j, va, 1)
            call ga_get (g_b, i, i, j, j, vb, 1)
            diff = va - vb
            if (abs(diff) > tol) then
               if (me.eq.0) then
                  write (luout, "(i6,x,i6,5x,2e12.3,5x,2e12.3)")
     $                 i, j,
     $                 dble(va), dble(vb), 
     $                 aimag(va), aimag(vb)
C     $                 i, j, va, vb, diff
               endif
            endif
         enddo
      enddo

      end subroutine
c $Id$
