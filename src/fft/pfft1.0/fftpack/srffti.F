#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine srffti (n,wsave)
      IMPLICIT REAL(A-H), INTEGER(I-N), REAL(O-Z)
      dimension  wsave(*)
      if (n .eq. 1) return
      call srfti1 (n,wsave(n+1),wsave(2*n+1))

      return
      end
