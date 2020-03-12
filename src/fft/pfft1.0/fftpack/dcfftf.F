#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine dcfftf (n,c,wsave)
*
* $Id$
*
      double precision c(1), wsave(1)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
      end
