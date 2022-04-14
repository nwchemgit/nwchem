#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine scfftf (n,c,wsave)
      IMPLICIT REAL(A-H), INTEGER(I-N), REAL(O-Z)
      dimension  c(*),wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call scftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
