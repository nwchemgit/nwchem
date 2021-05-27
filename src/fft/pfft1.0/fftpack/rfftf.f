      subroutine rfftf (n,r,wsave)
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
