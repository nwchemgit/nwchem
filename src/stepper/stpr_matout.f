      SUBROUTINE stpr_matout(A,NR,NC,M)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
C
C     Rectangular array output.
C     input:
C       A    Array to be written to for006.
C       NR   The row order of A.
C       NC   The column order of A.
C       M    The row dimension of A in the calling routine.
C
      DIMENSION A(M,NC)
      CHARACTER*15 LINE(5)
      DATA LINE /5*' _____________ '/
      MAXCOL=0
      NSETS=(NC-1)/5+1
      DO 100 NS=1,NSETS
        MINCOL=MAXCOL
        IF (NS.EQ.NSETS) THEN
          NUMCOL=NC-MINCOL
        ELSE
          NUMCOL=5
        END IF
        MAXCOL=MINCOL+NUMCOL
        MINCOL=MINCOL+1
        WRITE (6,1000) (I,I=MINCOL,MAXCOL)
        WRITE (6,1010) (LINE(I),I=1,NUMCOL)
        DO 90 I=1,NR
          WRITE (6,1020) I,(A(I,J),J=MINCOL,MAXCOL)
   90   CONTINUE
  100 CONTINUE
      WRITE (6,1030)
 1000 FORMAT (/,5X,5(6X,I3,6X))
 1010 FORMAT (5X,5A15)
 1020 FORMAT (1X,I3,1X,5(E14.7,1X))
 1030 FORMAT (/)
      RETURN
      END
