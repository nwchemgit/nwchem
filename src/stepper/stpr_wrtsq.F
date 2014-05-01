C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_wrtsq (IUNIT,A,LEN,IERR)
c $Id$
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
C> @}
