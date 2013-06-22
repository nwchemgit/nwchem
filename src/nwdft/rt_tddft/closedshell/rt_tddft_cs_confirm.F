C
C     Checks that parameters are consistent with closed shell
C     calculation.  Designed to be called as "rt_tddft_cs_confirm
C     (params)", which uses preprocessor magic to supply filename and
C     linenumber where it was called from (see rt_tddft.fh).
C
C     Note strange routine name--short so that preprocessor doesnt go
C     beyond silly fortran line width limit.
C
      subroutine rtcscfl (params, fname, lnum)

      implicit none

#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      character(len=*), intent(in)  :: fname
      integer, intent(in)           :: lnum


C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_cs_confirm: "


C     == Variables ==
      integer me
      character(len=256)   :: outstr
      logical ok

      
      me = ga_nodeid ()

      write (outstr, "(a,i0)")
     $     pname//"checks failed at "//trim(fname)//":", lnum
      

C
C     Actual checks.
C
      ok = .true.

      if (params%theory .ne. "CS") then
         if (me.eq.0) then
            write (luout, "(a)") pname//
     $           "theory should be CS, found: "//params%theory
         endif
         ok = .false.
      endif

      if (params%ipol.ne.1) then
         if (me.eq.0) then
            write (luout, "(a,i0)") pname//
     $           "ipol should equal 1, found: ", params%ipol
         endif
         ok = .false.
      endif

      if (params%nbf_ao .ne. params%ns_ao) then
         if (me.eq.0) then
            write (luout, "(a,i0,a,i0)") pname//
     $           "ns_ao should equal nbf_ao, found: ",
     $           params%ns_ao, ", ", params%nbf_ao
         endif
         ok = .false.
      endif

      if (params%nbf_mo .ne. params%ns_mo) then
         if (me.eq.0) then
            write (luout, "(a,i0,a,i0)") pname//
     $           "ns_mo should equal nbf_mo, found: ",
     $           params%ns_mo, ", ", params%nbf_mo
         endif
         ok = .false.
      endif

      
C
C     If all checks didnt pass, exit with error message.
C
      if (.not. ok) then
         if (me.eq.0) then
            write (luout, "(a)") trim(outstr)
            write (luout, *) ""
            call util_flush (luout)
         endif
         call errquit (outstr, 0, 0)
      endif

      end subroutine 
c $Id$
