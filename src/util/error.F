      block data error_data
C$Id: error.F,v 1.6 1995-02-02 18:09:21 d3g681 Exp $
      implicit none
#include "errorP.h"
c
      data fatal_level / 0 /
      data warn_level / 0 /
c
      end
      subroutine error(level, message, status)
      implicit none
      integer level             ! [input]
      character*(*) message     ! [input]
      integer status            ! [input]
c
c     Central error handling facility
c
c     Fatal_level and warn_level are internally maintained
c
c     0 <= fatal_level <= warn_level <= 2
c
c     if (level <= fatal_level) then
c        print out level, message and status as an error message
c        terminate execution abnormally with given status
c     else if (level <= warn_level || level > 2) then
c        print out level, message and status as a warning message
c        continue execution
c     else
c        quietly continue execution
c     endif
c
      if (level .le. fatal_level) then
         write(6, 1) message, icode
 1       format(//' fatal error(',i2,'): ',a,'(',i,')'/)
#ifdef TCGMSG
         call parerr(status)
#endif
#if defined(PARAGON)
	Call Kill(0, 9)
#elif defined(DELTA)
         call killproc(-1,0)
#elif defined(IPSC)
         call killcube(-1,-1)
#endif
#if defined(SUN) || defined(SGI)
         call abort
#else
         stop status
#endif
      else if (level .le. warn_level .or. level .gt. 2) then
         write(6, 1) level, message, icode
 1       format(' warning(',i2,'): ',a,'(',i,')')
      endif
c
      end
