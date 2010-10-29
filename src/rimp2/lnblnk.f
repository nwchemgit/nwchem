       INTEGER FUNCTION LNBLNK (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:   Returns the position of the last non-blank character
C
C Arguments: STRING   character string (input only)
C
C Remarks:   All FORTRAN 77 character variables are blank padded on the
C            right.  The intrinsic function LEN returns the dimension
C            of the character object, not the length of the contents.
C            
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Header: /tmp/mss/nwchem/src/rimp2/lnblnk.f,v 1.4 1995-12-16 21:06:36 gg502 Exp $
C
C    Revision 0.0  87/07/24  bernholdt (VAX)
C $Log: not supported by cvs2svn $
c Revision 1.3  1995/02/02  23:21:12  d3g681
c RJH: A CVS ID for every file and automated generation of a version output
c
c Revision 1.2  1994/09/01  21:07:41  d3e129
c removed call to char from parameter statement.  not allowed by
c fortran standard.
c RAK
c
c Revision 1.1  1994/06/14  21:54:19  gg502
c First cut at RIMP2.
c
c Revision 1.1  91/08/26  23:11:19  bernhold
c Initial revision
c 
C    Revision 1.1  88/01/11  22:08:15  bernhold
C    Initial revision
C    
C
C System:     Standard FORTRAN 77
C
C Copyright 1987 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C$Id$
       IMPLICIT NONE
       INTEGER I
       CHARACTER*(*) STRING
       CHARACTER*1 BLANK, NULL
       PARAMETER (BLANK = ' ')
C
C      Start at the end and work backwards
C
       NULL=CHAR(0)
       DO 100 I = LEN(STRING), 1, -1
C         Look for first non-whitespace character
          IF (STRING(I:I) .NE. BLANK .AND. String(I:I) .ne. NULL) THEN
             LNBLNK = I
             RETURN
          ENDIF
  100  CONTINUE
C
C      If we get this far, the string is empty
       LNBLNK = 0
       RETURN
       END
