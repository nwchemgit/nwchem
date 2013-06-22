C
C     rt_tddft_debug.F
C
C     Debugging routines.
C
      subroutine halt_fl (fname, lnum)
      implicit none

#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
      
C     == Inputs ==
      character(len=*), intent(in) :: fname
      integer, intent(in)          :: lnum

C     == Variables ==
      integer me
      character*10 alnum
      character*100 str
      integer leng
      integer i

      me = ga_nodeid ()

      write (alnum, "(i0)") lnum
      
      str = "Halted at " // trim (fname) // ":" // trim (alnum)
      leng = len_trim (str)


      if (me.eq.0) then
         write (luout, *) "########## " // trim (str) // " ##########"
         call util_flush (luout)
      endif

      call ga_sync ()
      call errquit ("(halted)", 0, 0)
      end subroutine


C====================================================================
      subroutine errquit_fl (msg, i1, i2, fname, lnum)
      implicit none

#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
      
C     == Inputs ==
      character(len=*), intent(in) :: msg
      integer, intent(in)          :: i1, i2
      character(len=*), intent(in) :: fname
      integer, intent(in)          :: lnum

      
C     == Variables ==
      integer me

      
      me = ga_nodeid ()
      
      if (me.eq.0) then
         write (luout, *)
     $        "errquit () called near line ", lnum, " of ", fname, "."
         call util_flush (luout)
      endif
      
      call errquit (msg, i1, i2)

      end subroutine



      subroutine prof_start (time)
      implicit none

#include "util.fh"

      double precision time
      
      time = util_cpusec ()
      end subroutine


      subroutine prof_end (time, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      double precision, intent(in) :: time
      character(len=*), intent(in) ::  msg


C     == Variables ==
      double precision elapsed
      integer me
      integer sec, msec

      me = ga_nodeid ()
      elapsed = util_cpusec () - time

      if (elapsed .lt. 10d0) then
         msec = ceiling (elapsed*1000)
         
         if (me.eq.0) then
            write (luout, "(1x, a, i0, a)")
C     $           out_tag//"[prof] "//trim (msg)//" took ", msec, " ms."
     $           "[prof] "//trim (msg)//" took ", msec, " ms."
            call util_flush (luout)
         endif
         
      else
         sec = ceiling (elapsed)
         
         if (me.eq.0) then
            write (luout, "(1x, a, i0, a)")
C     $           out_tag//"[prof] "//trim (msg)//" took ", sec, " s."
     $           "[prof] "//trim (msg)//" took ", sec, " s."
            call util_flush (luout)
         endif

      endif

      end subroutine




C
C     Graveyard
C
#if 0
      if (me.eq.0) then

         write (luout, *) ""
         write (luout, "(a)", advance="no") " "
         do i = 1, leng+6
            write (luout, "(a)", advance="no") "#"
         enddo
         write (luout, *) ""

         call util_flush (luout)

         write (luout, "(1x,a)") "#  " // trim (str) // "  #"

         call util_flush (luout)

         write (luout, "(a)", advance="no") " "
         do i = 1, leng+6
            write (luout, "(a)", advance="no") "#"
         enddo
         write (luout, *) ""
         write (luout, *) ""

         call util_flush (luout)
      endif

      call ga_sync ()
      call errquit ("(halted)", 0, 0)
#endif      
      
      
#if 0
      if (me.eq.0) then
         write (luout, *) "----------------------------------------"
         write (luout, "(1x,a,i0,a)")
     $        "Halted at "//fname//":", lnum, ""
         write (luout, *) "----------------------------------------"
         call util_flush (luout)
      endif
#endif

c $Id$
