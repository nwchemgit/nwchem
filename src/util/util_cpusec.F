      double precision function util_cpusec()
      implicit none
! we implement cpu time as wall time because cpu time
! is often unhelpful, especially when threads are used.
! there were issues with the implementation of cpu time
! as well, that motivated deleting the code.
      double precision util_wallsec
      external util_wallsec
      util_cpusec=util_wallsec()
      end
