      SUBROUTINE stpr_wrtsq (IUNIT,A,LEN,IERR)
c $Id: stpr_wrtsq.F,v 1.1 1994-10-26 23:26:14 d3e129 Exp $
      INTEGER LEN, IUNIT, IERR
      REAL*8 A(LEN)
#if defined (SGITFP)
      WRITE (IUNIT, ERR=10) (A(IERR),IERR=1,LEN)
#else
      WRITE (IUNIT, ERR=10) A
#endif
      IERR = 0
      RETURN
   10 CONTINUE
      IERR = 1
      RETURN
      END
