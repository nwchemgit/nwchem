*
* $Id: nwpw_timing.F,v 1.1 2001-12-17 23:31:14 bylaska Exp $
*

* nwpw_timing.F
* Author - Eric Bylaska
*
*   These routines are to be used to time the nwpw module
*

      subroutine nwpw_timing_init()
      implicit none

#include "nwpw_timing_common.fh"

	  call dcopy(50,0.0d0,0,times,1)
      return 
      end

      subroutine nwpw_timing_start(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"

	  call current_second(nwpw_tim1(counter))
      return
      end

      subroutine nwpw_timing_end(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"

	  call current_second(nwpw_tim2(counter))
      times(counter) = times(counter) 
     >               + (nwpw_tim2(counter)-nwpw_tim1(counter))
      return
      end

      real*8 function nwpw_timing(counter)
      implicit none
      integer counter
#include "nwpw_timing_common.fh"

      nwpw_timing = times(counter)
      return
      end

