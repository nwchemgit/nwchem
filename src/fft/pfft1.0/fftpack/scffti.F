#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine scffti (n,wsave)
      dimension       wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call scfti1 (n,wsave(iw1),wsave(iw2))
      return
      end
