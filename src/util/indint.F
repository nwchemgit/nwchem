       SUBROUTINE INDINT (STRING, START, END)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Find the position of the first valid integer in STRING
C
C Arguments:    STRING   character string (input only)
C               START    index of beginning of integer (output only)
C               END      index of last digit of integer (output only)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Header: /tmp/mss/nwchem/src/util/indint.F,v 1.1 1997-03-24 02:00:36 gg502 Exp $
C
C $Log: not supported by cvs2svn $
C Revision 1.1  1996/08/27 21:56:14  d3e129
C removed requirments of LLE, LGE, LLT, LGT for SGITFP.
C rak
C
C Revision 1.3  1995/12/16 21:06:34  gg502
C The great implicit none-ification :-)
C
c Revision 1.2  1995/02/02  23:21:10  d3g681
c RJH: A CVS ID for every file and automated generation of a version output
c
c Revision 1.1  1994/06/14  21:54:18  gg502
c First cut at RIMP2.
c
c Revision 1.1  91/08/26  23:09:37  bernhold
c Initial revision
c 
c Revision 1.1  90/05/14  17:13:09  bernhold
c Initial revision
c 
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C$Id$
       IMPLICIT NONE
       CHARACTER*(*) STRING
       CHARACTER*1 C
       INTEGER START, END
       INTEGER I
       LOGICAL INONE, SIGN
C
       START = 0
       END = 0
       INONE = .FALSE.
       SIGN = .FALSE.
C
       DO 100 I = 1, LEN(STRING)
          C = STRING(I:I)
          IF (INONE .AND. (LLT(C,'0') .OR. LGT(C,'9'))) THEN
C               End of an integer...
             INONE = .FALSE.
C               Check to make sure it wasn't just a + or a -
             IF (SIGN .AND. (START .EQ. I-1)) THEN
                SIGN = .FALSE.
                START = 0
             ELSE
                END = I - 1
                RETURN
             ENDIF
          ELSEIF (.NOT.INONE .AND. 
     1           (C.EQ.'-' .OR. C.EQ.'+')) THEN
C               This is a sign character, which may begin an integer
             SIGN = .TRUE.
             INONE = .TRUE.
             START = I
          ELSEIF (.NOT.INONE .AND. 
     1           (LGE(C,'0') .AND. LLE(C,'9'))) THEN
C               This is a digit, which does begin an integer
             INONE = .TRUE.
             START = I
          ENDIF
 100   CONTINUE
C
C      If we as still INONE here, the integer is at the end of STRING
C
       IF (INONE) END = LEN(STRING)
       RETURN
       END
