C $Id$
C**********************************************************************
C CFLOAT
C**********************************************************************
C
      double precision FUNCTION drdy_cfloat(STRING)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "errquit.fh"
C
#include "drdyP.fh"
C
      CHARACTER*(*) STRING
      CHARACTER*80  NUMBER
      CHARACTER CH
      LOGICAL LEXP,LDEC
C
*      write(fu6,*)' drdy_cfloat:string: <',string,'>'
      LEXP = .FALSE.
      LDEC = .FALSE.
      LENGTH = LEN(STRING)
      IF (LENGTH .EQ. 0) THEN
         DRDY_CFLOAT = 0.0D0
         RETURN
      ENDIF
C     WRITE(fu6,*) LENGTH,STRING
C
C     Find the first nonblank character
C
      I = 1
10    IF (STRING(I:I) .EQ. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 10
      ENDIF
C
C     If it is a blank string set function to zero
C
      IF (I .GT. LENGTH) THEN
         DRDY_CFLOAT = 0.0D0
         RETURN
      ENDIF
      IBEG = I
C
C     Find the first blank character after the number
C
      I = IBEG+1
20    IF (STRING(I:I) .NE. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 20
      ENDIF
      IEND = I-1
C
C     Stripe the blanks before and after the number
C
      NUMBER = STRING(IBEG:IEND)
      LENGTH = IEND - IBEG + 1
C   
C     Make sure there is no blank left
C
      IF (INDEX(NUMBER,' ') .LE. LENGTH) THEN
         WRITE(fu6,1000) STRING
         call errquit ('drdy_cfloat:1: fatal error',911, UNKNOWN_ERR)
      ENDIF
C
C     Find the decimal point
C
      IDEC = INDEX(NUMBER,'.')
      IF (IDEC .NE. 0) LDEC = .TRUE.
C
C     Find the exponential symbol
C
      IUE = INDEX(NUMBER,'E')
      ILE = INDEX(NUMBER,'e')
      IUD = INDEX(NUMBER,'D')
      ILD = INDEX(NUMBER,'d')
      ISUM = IUE + ILE + IUD + ILD
      IEXP = MAX0(IUE,ILE,IUD,ILD)
      IF (ISUM .GT. IEXP) THEN
         WRITE(fu6,1000) STRING
         call errquit('drdy_cfloat:2: fatal error',911, UNKNOWN_ERR)
      ENDIF
      IF (IEXP .NE. 0) THEN
         LEXP = .TRUE.
      ELSE
         IEXP = LENGTH + 1
      ENDIF
C
      IF (.NOT. LDEC) IDEC = IEXP
C
C     Get the number before decimal
C
      IBEG = 2
      IF (NUMBER(1:1) .EQ. '+') THEN
         SIGN = 1.0D0
      ELSEIF(NUMBER(1:1) .EQ. '-') THEN
         SIGN = -1.0D0
      ELSE
         SIGN = 1.0D0
         IBEG = 1
      ENDIF
      IF (IBEG .EQ. IEXP) THEN
         F1 = 1.0D0
      ELSE
         F1 = 0.0D0
      ENDIF
      DO 50 I = IBEG,IDEC-1
         CH = NUMBER(I:I)
         IF (CH .GE. '0' .AND. CH .LE. '9') THEN
            N = ICHAR(CH) - ICHAR('0')
            F1 = F1 * 10.0D0 + DBLE(N)
         ELSE
            WRITE(fu6,1000) STRING
            call errquit('drdy_cfloat:3: fatal error',911, UNKNOWN_ERR)
         ENDIF
50    CONTINUE
C
C     Get the number after decimal 
C
      F2 = 0.0D0
      IF (LDEC) THEN
         J = 0
         DO 60 I = IDEC+1,IEXP-1
            CH = NUMBER(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F2 = F2 * 10.0D0 + DBLE(N)
               J = J + 1
            ELSE
               WRITE(fu6,1000) STRING
               call errquit('drdy_cfloat:4: fatal error',911,
     &       UNKNOWN_ERR)
            ENDIF
60       CONTINUE
         F2 = F2 / 10.0D0 ** DBLE(J)
      ENDIF
C
C    Get the exponent
C
      ESIGN = 1.0D0
      F3 = 0.0D0
      IF (LEXP) THEN 
         IBEG = IEXP + 2
         IF (NUMBER(IEXP+1:IEXP+1) .EQ. '+') THEN
            ESIGN = 1.0D0
         ELSEIF(NUMBER(IEXP+1:IEXP+1) .EQ. '-') THEN
            ESIGN = -1.0D0
         ELSE
            ESIGN = 1.0D0
            IBEG = IEXP + 1
         ENDIF
         DO 70 I = IBEG,LENGTH
            CH = NUMBER(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F3 = F3 * 10.0D0 + DBLE(N)
            ELSE
               WRITE(fu6,1000) STRING
               call errquit('drdy_cfloat:5: fatal error',911,
     &       UNKNOWN_ERR)
            ENDIF
70       CONTINUE
      ENDIF 
C
      DRDY_CFLOAT = (SIGN * (F1 + F2)) * 10.0D0 ** (ESIGN*F3)
C
      RETURN
C
1000  FORMAT(/1X,'Illegal number: ',A80)
C
      END
