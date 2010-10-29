      SUBROUTINE stpr_uprcase (ANSWER)
c $Id$
C
C     Converts ANSWER to upper case letters.
C
      CHARACTER*(*) ANSWER
      IA=ICHAR('a')
      IZ=ICHAR('z')
      IDIF=ICHAR('A')-IA
      DO 10 I=1,LEN(ANSWER)
        J=ICHAR(ANSWER(I:I))
        IF ((J.GE.IA).AND.(J.LE.IZ)) ANSWER(I:I)=CHAR(J+IDIF)
   10 CONTINUE
      RETURN
      END
