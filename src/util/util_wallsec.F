      double precision function util_wallsec()
      implicit none
c $Id: util_wallsec.F,v 1.7 2004-02-10 17:34:18 edo Exp $
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


