C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_readsq (IUNIT,A,LEN,IERR)
c $Id$
      INTEGER LEN, IUNIT, IERR
      REAL*8 A(LEN)
      READ (IUNIT, END=10, ERR=10) A
      IERR = 0
      RETURN
   10 CONTINUE
      IERR = 1
      RETURN
      END
C> @}
