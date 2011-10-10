
      subroutine print_elapsed_time(autime)
      implicit none
      real*8 autime

#include "stdio.fh"

      real*8 sectime

      sectime = autime*2.41889d-17

      if (sectime.lt.1.0d-12) then
         write(luout,1800) (sectime/1.0d-15)," fs"
      else if (sectime.lt.1.0d-9) then
         write(luout,1800) (sectime/1.0d-12)," ps"
      else
         write(luout,1800) (sectime/1.0d-9 )," ns"
      end if

      return
 1800 format(//' Elapsed time of simulation was',F8.3,A)
      end

c $Id$
