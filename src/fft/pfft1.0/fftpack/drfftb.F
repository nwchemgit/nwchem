      subroutine drfftb (n,r,wsave)
      double precision r(1), wsave(1)
c
      if (n .eq. 1) return
c
      call drftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
