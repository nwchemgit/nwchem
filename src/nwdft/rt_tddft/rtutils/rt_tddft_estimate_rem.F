C
C     rt_tddft_estimate_rem.F
C
C     Prints out roughly how much time remaining in simulation.  This is
C     based on how long it has taken since the start of the propagation
C     (i.e, the number of time steps completed and the cpu time taken).
C     
      subroutine rt_tddft_estimate_rem (params, cputime_start, it_start,
     $     it, tt)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: cputime_start
      integer, intent(in)           :: it_start       !starting index (>1 if did a restart)
      integer, intent(in)           :: it             !current index
      double precision, intent(in)  :: tt             !current time


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_estimate_rem: "


C     == Variables ==
      double precision cputime_curr
      double precision cputime_taken
      double precision time_per_step, time_remaining
      integer stepsdid, stepsleft
      character(100) outstr
      integer nsec_tot
      integer nsec, nmin, nhour, nday

      
      cputime_curr = util_cpusec ()
      cputime_taken = cputime_curr - cputime_start

      stepsdid = it - it_start
      stepsleft = params%nt - it

      if (stepsdid .le. 0)
     $     call errquit (pname//"bad num stepsdid", 0, 0)

      if (stepsleft .lt. 0)
     $     call errquit (pname//"bad num stepsleft", 0, 0)

      time_per_step = cputime_taken / stepsdid
      time_remaining = time_per_step * stepsleft

      
      nsec_tot = nint (time_remaining)


C
C     Generate string to print (human readable time format)
C
      if (nsec_tot .le. 60) then !seconds
         write (outstr, "(a,i0,a)") "Estimated time remaining: ",
     $        nsec_tot, " second(s)"
         
      elseif (nsec_tot.le.3600) then !minutes
         nmin = nsec_tot / 60
         nsec = mod (nsec_tot, 60)
         
         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nmin, " minute(s) ", nsec, " second(s)"

      elseif (nsec_tot.le.86400) then !hours
         nhour = nsec_tot / 3600
         nsec = mod (nsec_tot, 3600)
         nmin = nsec / 60

         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nhour, " hour(s) ", nmin, " minute(s)"
         
      elseif (nsec_tot.le.604800) then !days
         nday = nsec_tot / 86400
         nsec = mod (nsec_tot, 86400)
         nhour = nsec / 3600

         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nday, " day(s) ", nhour, " hour(s)"

      else !weeks
         nday = nsec_tot / 86400         

         write (outstr, "(a,i0,a)") "Estimated time remaining: ",
     $        nday, " day(s)"

      endif
         
      call rt_tddft_print_msg (params, tt, outstr)
      
      end subroutine



#if 0
C
C     rt_tddft_estimate_rem.F
C
C     Prints out roughly how much time remaining in simulation.  This is
C     based on how long it has taken since the start of the propagation
C     (i.e, the number of time steps completed and the cpu time taken).
C     
      subroutine rt_tddft_estimate_rem (params, cputime_start, it_start,
     $     it, tt)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: cputime_start
      integer, intent(in)           :: itstart        !starting index (>1 if did a restart)
      integer, intent(in)           :: it             !current index
      double precision, intent(in)  :: tt             !current time


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_estimate_rem: "


C     == Variables ==
      double precision cputime_curr
      double precision cputime_taken
      double precision time_per_step, time_remaining
      integer stepsdid, stepsleft
      character(100) outstr
      integer nsec, nmin, nhour, nday

      cputime_curr = util_cpusec ()
      cputime_taken = cputime_curr - cputime_start

      stepsdid = it - it_start
      stepsleft = params%nt - it

      if (stepsdid .le. 0)
     $     call errquit (pname//"bad num stepsdid", 0, 0)

      if (stepsleft .lt. 0)
     $     call errquit (pname//"bad num stepsleft", 0, 0)

      time_per_step = cputime_taken / stepsdid
      time_remaining = time_per_step * stepsleft

      
C
C     Generate string to print (human readable time format)
C
      if (time_remaining .le. 60d0) then !seconds
         nsec = nint (time_remaining)
         write (outstr, "(a,i0,a)") "Estimated time remaining: ",
     $        nsec, " second(s)"
         
      elseif (time_remaining.le.3600d0) then !minutes
         nmin = nint (time_remaining / 60d0)
         nsec = nint (time_remaining - nmin*60)

         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nmin, " minute(s) ", nsec, " second(s)"
         
      elseif (time_remaining.le.86400d0) then !hours
         nhour = nint (time_remaining / 3600d0)
         nmin = nint ( (time_remaining - nhour*3600)/60d0 )

         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nhour, " hour(s) ", nmin, " minute(s)"

      elseif (time_remaining.le.604800d0) then !days
         nday = nint (time_remaining / 86400d0)
         nhour = nint ( (time_remaining - nday*86400d0) / 3600d0 )

         write (outstr, "(a,i0,a,i0,a)") "Estimated time remaining: ",
     $        nday, " day(s) ", nhour, " hour(s)"

      else !weeks
         
         write (outstr, "(a,a)") "Estimated time remaining: ",
     $        "greater than one week"
         
      endif

      call rt_tddft_print_msg (params, tt, outstr)
      
      end subroutine
#endif
