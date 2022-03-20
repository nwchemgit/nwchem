       INTEGER FUNCTION ATOI (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Convert STRING to an integer value
C
C Arguments:    STRING   character string (input only)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Header: /tmp/mss/nwchem/src/util/atoi.F,v 1.1 1997-03-24 02:00:35 gg502 Exp $
C
C $Log: not supported by cvs2svn $
C Revision 1.1  1996/08/27 21:56:09  d3e129
C removed requirments of LLE, LGE, LLT, LGT for SGITFP.
C rak
C
C Revision 1.3  1995/12/16 21:06:29  gg502
C The great implicit none-ification :-)
C
c Revision 1.2  1995/02/02  23:21:06  d3g681
c RJH: A CVS ID for every file and automated generation of a version output
c
c Revision 1.1  1994/06/14  21:54:07  gg502
c First cut at RIMP2.
c
c Revision 1.1  91/08/26  23:09:29  bernhold
c Initial revision
c 
C     Revision 1.1  88/12/07  13:38:52  bernhold
C     Initial revision
C     
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C$Id$
       IMPLICIT NONE
       CHARACTER*(*) STRING
       CHARACTER*1 C
       INTEGER I
       LOGICAL NEG
C
       ATOI = 0
       NEG = .FALSE.
       I = 1
       C = STRING(I:I)
C
C      Pass over any leading spaces
C
 100   IF (C .EQ. ' ') THEN
          I = I + 1
          C = STRING(I:I)
          GOTO 100
       ENDIF
C
C      See if first character makes the number negative
C      Accept + as valid character before the digits start
C
       IF (C .EQ. '-') THEN
          NEG = .TRUE.
          I = I + 1
          C = STRING(I:I)
       ELSEIF (C .EQ. '+') THEN
          I = I + 1
          C = STRING(I:I)
       ENDIF
C
C      Continue as long as its a digit...
C
 200   IF (LGE(C, '0') .AND. LLE(C,'9')) THEN
C            Shift number over & add new digit
          ATOI = 10*ATOI + ICHAR(C)-48
          I = I + 1
          C = STRING(I:I)
          GOTO 200
       ENDIF
C
C      Negate the result if necessary
C
       IF (NEG) ATOI = -ATOI
       RETURN
       END
