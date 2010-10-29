      double precision function util_wallsec()
      implicit none
c $Id$
*
* routine to give wall clock seconds since start of execution
*  delta times are handled by all applications
*
*
      double precision tcgtime
      external tcgtime
      util_wallsec = tcgtime()
*
      end


