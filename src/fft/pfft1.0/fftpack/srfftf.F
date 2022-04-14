#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine srfftf (n,r,wsave)
      IMPLICIT REAL(A-H), INTEGER(I-N), REAL(O-Z)
      dimension       r(*)       ,wsave(*)

      if (n .eq. 1) return
      call srftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
