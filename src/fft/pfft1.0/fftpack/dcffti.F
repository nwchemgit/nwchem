      subroutine dcffti (n,wsave)
*
* $Id: dcffti.F,v 1.2 1997-11-04 10:07:17 d3e129 Exp $
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
