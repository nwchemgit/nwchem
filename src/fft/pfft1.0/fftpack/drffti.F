      subroutine drffti (n,wsave)
*
* $Id: drffti.F,v 1.2 1997-11-04 10:07:29 d3e129 Exp $
*
      double precision wsave(1)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
