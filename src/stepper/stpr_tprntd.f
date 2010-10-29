      SUBROUTINE stpr_tprntd(A,NR)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      CHARACTER*15 LINE(5)
C
C     This routine writes a triangular packed real array, A,
C     of order NR to for006.
C
      DIMENSION A(NR*(NR+1)/2)
      DATA LINE /5*' _____________ '/
      ISYM2(I,J)=I*(I-1)/2+J
      MAXCOL=0
      NSETS=(NR-1)/5+1
      DO 100 NS=1,NSETS
        MINCOL=MAXCOL
        IF (NS.EQ.NSETS) THEN
          NUMCOL=NR-MINCOL
        ELSE
          NUMCOL=5
        END IF
        MAXCOL=MINCOL+NUMCOL
        MINCOL=MINCOL+1
        WRITE(6,1000)(I,I=MINCOL,MAXCOL)
        WRITE(6,1010)(LINE(I),I=1,NUMCOL)
        DO 90 I=MINCOL,NR
          MXCOL = MIN0(MAXCOL,I)
          WRITE (6,1020) I,(A(ISYM2(I,J)),J=MINCOL,MXCOL)
   90   CONTINUE
  100 CONTINUE
      WRITE (6,1030)
 1000 FORMAT (/,5X,5(6X,I3,6X))
 1010 FORMAT (5X,5A15)
 1020 FORMAT (1X,I3,1X,5(D14.7,1X))
 1030 FORMAT (/)
      RETURN
      END
