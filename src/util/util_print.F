      logical function util_print(name, level)
      implicit none
#include "printlevels.fh"
      character*(*) name
      integer level
c
      util_print = level .le. print_default
c
      end
