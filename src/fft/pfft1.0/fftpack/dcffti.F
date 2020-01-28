#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine dcffti (n,wsave)
*
* $Id$
*
      double precision wsave(1)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcfti1 (n,wsave(iw1),wsave(iw2))
c
      return
      end
