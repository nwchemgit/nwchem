      double precision function util_wallsec()
      implicit none
c $Id: util_wallsec.F,v 1.5 2001-10-31 21:10:43 windus Exp $
*
* routine to give wall clock seconds since start of execution
*  delta times are handled by all applications
*
* tcgtime is more accurate on all (?) platforms
*
#define USETCGTIME
#if defined(USETCGTIME)
      double precision tcgtime
      external tcgtime
      util_wallsec = tcgtime()
#elif defined(IBM) || defined(SP1)
      integer time
      external time
      util_wallsec = dble(time(0))
#elif defined(KSR)
      real all_seconds
      external all_seconds
      util_wallsec = all_seconds()
#elif defined(SUN) || defined(SGI) || defined(SGITFP) || defined(SOLARIS) || defined(FUJITSU_VPP) || defined(FUJITSU_SOLARIS)
      integer time
      external time
      util_wallsec = dble(time())
#elif defined(__PARAGON__) || defined(DELTA)
      double precision dclock
      external dclock
      util_wallsec = dclock()
#else
      util_wallsec = 0.0d00
#endif
*
      end


