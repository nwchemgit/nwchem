      double precision function util_wallsec()
      implicit none
c $Id: util_wallsec.F,v 1.1 1995-07-10 17:02:22 d3e129 Exp $
*
* routine to give wall clock seconds since start of execution
*  delta times are handled by all applications
*
#if defined(IBM) || defined(SP1)
#define __UTIL_WALLSEC__
      integer time
      external time
      util_wallsec = dble(time(0))
#endif
#if defined(KSR)
#define __UTIL_WALLSEC__
*--KSR-start
      real all_seconds
      external all_seconds
c
      util_wallsec = all_seconds()
*--KSR-end
#endif
#if defined(SUN) || defined(SGI) || defined(SGITFP)
#define __UTIL_WALLSEC__
      integer time
      external time
      util_wallsec = dble(time())
#endif
#if defined(__PARAGON__)
#define __UTIL_WALLSEC__
      double precision dclock
      external dclock
      util_wallsec = dclock()
#endif
*
#if !defined(__UTIL_WALLSEC__)
      util_wallsec = 0.0d00
#endif
*
      end


