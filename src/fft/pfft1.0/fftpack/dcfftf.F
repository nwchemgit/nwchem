      subroutine dcfftf (n,c,wsave)
*
* $Id: dcfftf.F,v 1.2 1997-11-04 10:07:17 d3e129 Exp $
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
