      double precision function util_cpusec()
      implicit none
! we implement cpu time as wall time because cpu time
! is often unhelpful, especially when threads are used.
! there were issues with the implementation of cpu time
! as well, that motivated deleting the code.
      return util_wallsec()
      end
