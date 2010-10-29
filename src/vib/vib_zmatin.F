      SUBROUTINE ZMATIN(IZMAT,NELS)
C*MODULE VIB  *DECK ZMATIN
* $Id$
C
      IMPLICIT NONE
#include "nwc_const.fh"
*: input elsewhere now :      LOCICAL NVAROK
C
      INTEGER MXATM
      PARAMETER (MXATM=nw_max_atom)
      double precision zan,c
      COMMON /cvib_INFOA /  ZAN(MXATM),C(3,MXATM)
c
      integer nrow, ndeg
      COMMON /cvib_VBDCMP/ NROW,NDEG
c
      LOGICAL LINEAR
      integer nzmat, nzvar, nvar
      COMMON /cvib_ZMAT  / NZMAT,NZVAR,NVAR,LINEAR
c
      integer nzmat2, nzvar2, nvar2, nzmtrd
      COMMON /cvib_ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD
c
      integer nat, nat3, nhess, nhesst
      COMMON /cvib_HESS/ NAT,NAT3,NHESS,NHESST   ! Hessian information
C
      INTEGER NELS
      INTEGER IZMAT(NELS)
C
*: input elsewhere now :      INTEGER NNAM
*: input elsewhere now :      PARAMETER (NNAM=1)
*: input elsewhere now :      DIMENSION QNAM(NNAM),KQNAM(NNAM)
C
*: input elsewhere now :      DATA QNAM /8HIZMAT   /
C
C     ----- THIS ROUTINE READS THE CODED Z-MATRIX -----
C     IN COMMON /cvib_ZMAT/, NZVAR IS THE NUMBER OF INPUT
C     COORDINATES, NVAR IS THE TOTAL NUMBER OF COORDINATES,
C     AND NZMAT IS THE TOTAL NUMBER OF CODE INTEGERS INPUT.
C
C----------Read NVAR, number of internal coordinates from
C          standard input.
C
*rak flags from input
*: input elsewhere now :      LINEAR = .TRUE.
*: input elsewhere now :      CALL FIND('LINEAR',IERR)
*: input elsewhere now :      IF(IERR .NE. 0)LINEAR = .FALSE.
*: input elsewhere now :      CALL FIND('NVAR  ',IERR)
*: input elsewhere now :      IF(IERR .NE. 0)THEN
*: input elsewhere now :        WRITE(6,*)'Did not find flag NVAR.  Call Bomb.'
*: input elsewhere now :        CALL BOMB(0)
*: input elsewhere now :      ENDIF
*: input elsewhere now :      READ(5,*)NZVAR
*: input elsewhere now :C
*: input elsewhere now :      KQNAM(1) = 10*7*NZVAR + 1
*: input elsewhere now :      DO 100 I=1,NELS
*: input elsewhere now :         IZMAT(I) = 0
*: input elsewhere now :  100 CONTINUE
*: input elsewhere now :C
*: input elsewhere now :C----------READ IN DESCRIPTION OF INTERNAL COORDINATES.
*: input elsewhere now :C
*: input elsewhere now :       CALL FIND('INTCRD',IERR)
*: input elsewhere now :       IF(IERR.NE.0)THEN
*: input elsewhere now :         WRITE(6,*)'Could not find flag INTCRD.  Call Bomb.'
*: input elsewhere now :         CALL BOMB(0)
*: input elsewhere now :       ENDIF
*: input elsewhere now :       IPT = 1
*: input elsewhere now :       DO 7117 IVAR = 1,NZVAR
*: input elsewhere now :         READ(5,*)IZMAT(IPT)
*: input elsewhere now :         ITYPE = IZMAT(IPT)
*: input elsewhere now :         IF(ITYPE .EQ. 1)THEN
*: input elsewhere now :          IPT2 = IPT + 2
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 2)THEN
*: input elsewhere now :          IPT2 = IPT + 3
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 3)THEN
*: input elsewhere now :          IPT2 = IPT + 4
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 4)THEN
*: input elsewhere now :          IPT2 = IPT + 4
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 5)THEN
*: input elsewhere now :          IPT2 = IPT + 3
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 6)THEN
*: input elsewhere now :          IPT2 = IPT + 5
*: input elsewhere now :         ELSEIF(ITYPE .EQ. 7)THEN
*: input elsewhere now :          IPT2 = IPT + 6
*: input elsewhere now :         ELSE
*: input elsewhere now :          WRITE(6,9001)ITYPE
*: input elsewhere now :          CALL BOMB(0)
*: input elsewhere now :         ENDIF
*: input elsewhere now :         READ(5,*)(IZMAT(I),I=IPT+1,IPT2)
*: input elsewhere now :         IPT = IPT2 + 1
*: input elsewhere now : 7117  CONTINUE
*: input elsewhere now :C
*: input elsewhere now :C         INPUT OF INTERNALS REQUIRED IF COORD=UNIQUE, HINT, OR CART,
*: input elsewhere now :C         BUT FOR COORD=ZMT, ZMTMPC WE CAN USE $DATA INFO INSTEAD
*: input elsewhere now :C
*: input elsewhere now :C
*: input elsewhere now :C     ----- SCAN AND ANALYZE IZMAT -----
*: input elsewhere now :C     THE VALID TYPES ARE
*: input elsewhere now :C        1= STRETCH I-J
*: input elsewhere now :C        2= BEND I-J-K
*: input elsewhere now :C        3= TORSION I-J-K-L
*: input elsewhere now :C           (BETWEEN TWO PLANES HAVING TWO ATOMS IN COMMON)
*: input elsewhere now :C        4= OUT OF PLANE BEND I-J-K-L
*: input elsewhere now :C        5= LINEAR BEND I-J-K
*: input elsewhere now :C        6= 5 ATOM DIHEDRAL ANGLE I-J-K-L-M
*: input elsewhere now :C           (BETWEEN TWO PLANES HAVING ONE ATOM IN COMMON)
*: input elsewhere now :C        7= 6 ATOM GHOST TORSION ANGLE I1-I2-J-K-L1-L2
*: input elsewhere now :C
*: input elsewhere now :  200 CONTINUE
*: input elsewhere now :c      WRITE(6,9048)
*: input elsewhere now :      NVAR=0
*: input elsewhere now :      IZ = 1
*: input elsewhere now :      DO 230 I=1,NZVAR
*: input elsewhere now :         ITYPE=IZMAT(IZ)
*: input elsewhere now :         NCHK=0
*: input elsewhere now :         IF(ITYPE.EQ.1) NCHK=2
*: input elsewhere now :         IF(ITYPE.EQ.2) NCHK=3
*: input elsewhere now :         IF(ITYPE.EQ.3) NCHK=4
*: input elsewhere now :         IF(ITYPE.EQ.4) NCHK=4
*: input elsewhere now :         IF(ITYPE.EQ.5) NCHK=3
*: input elsewhere now :         IF(ITYPE.EQ.6) NCHK=5
*: input elsewhere now :         IF(ITYPE.EQ.7) NCHK=6
*: input elsewhere now :         IF(NCHK.EQ.0) THEN
*: input elsewhere now :            WRITE(6,9068) I,IZ,IZMAT(IZ)
*: input elsewhere now :            CALL BOMB(0)
*: input elsewhere now :         END IF
*: input elsewhere now :C
*: input elsewhere now :c         WRITE(6,9058) I,ITYPE,(IZMAT(IZ+J),J=1,NCHK)
*: input elsewhere now :         DO 210 J=1,NCHK
*: input elsewhere now :            IAT=IZMAT(IZ+J)
*: input elsewhere now :            IF(IAT.LE.0  .OR.  IAT.GT.NAT) THEN
*: input elsewhere now :               WRITE(6,9088) I,IZ,IAT,J,NAT
*: input elsewhere now :               CALL BOMB(0)
*: input elsewhere now :            END IF
*: input elsewhere now :  210    CONTINUE
*: input elsewhere now :C
*: input elsewhere now :         IZ=IZ+NCHK+1
*: input elsewhere now :         NVAR=NVAR+1
*: input elsewhere now :         IF(ITYPE.EQ.5) NVAR=NVAR+1
*: input elsewhere now :  230 CONTINUE
*: input elsewhere now :      NZMAT = IZ - 1
*: input elsewhere now :C
*: input elsewhere now :C --- IF NOT DOING VIBRATIONAL ANALYSIS, MUST SPECIFY EXACTLY 3N-NDEG
*: input elsewhere now :C --- INTERNAL COORDINATES.  OTHERWISE, MORE THAN 3N-NDEG IS OK.
*: input elsewhere now :C
*: input elsewhere now :      NDEG = 6
*: input elsewhere now :      IF (LINEAR) NDEG = 5
*: input elsewhere now :      NROW = NVAR + NDEG
*: input elsewhere now :      NVAROK = NVAR .EQ. (3*NAT-NDEG)
*: input elsewhere now :      IF (.NOT.NVAROK) THEN
*: input elsewhere now :         WRITE(6,9148) NDEG,NAT,NVAR
*: input elsewhere now :         CALL BOMB(0)
*: input elsewhere now :      END IF
*: input elsewhere now :C
*: input elsewhere now :C     ----- SAVE IZMAT ON TAPE 10 ----
*: input elsewhere now :C
*: input elsewhere now :c       CALL ADDARRAY(10,IZMAT,INTOWP(NELS),393)
*: input elsewhere now :C
*: input elsewhere now :      NVAR2  = NVAR
*: input elsewhere now :      NZMAT2 = NZMAT
*: input elsewhere now :      RETURN
*: input elsewhere now :C
*: input elsewhere now : 9000 FORMAT(1X,'COULD NOT LOCATE U:VIB')
*: input elsewhere now : 9001 FORMAT(1X,'UNKNOWN INTERNAL COORDINATE TYPE IN ',
*: input elsewhere now :     *          'INPUT DECK. ITYPE = ',I5)
*: input elsewhere now : 9008 FORMAT(1X,'ERROR IN $ZMAT INPUT - ERROR')
*: input elsewhere now : 9048 FORMAT(/3X,'--- ENCODED Z MATRIX ---'/
*: input elsewhere now :     *        1X,'COORD',2X,'TYPE',3X,'I',3X,'J',3X,'K',3X,'L',3X,'M',
*: input elsewhere now :     *        3X,'N')
*: input elsewhere now : 9058 FORMAT(1X,I3,3X,I4,1X,6I4)
*: input elsewhere now : 9068 FORMAT(1X,'INTERNAL COORD.',I3,' AT SUBSCRIPT',I4,
*: input elsewhere now :     *     ' IS OF ILLEGAL TYPE',I5)
*: input elsewhere now : 9088 FORMAT(1X,'INTERNAL COORD.',I3,' AT SUBSCRIPT',I4,
*: input elsewhere now :     *    ' IS FOLLOWED BY ILLEGAL ATOM',I4,' AT OFFSET',I4/
*: input elsewhere now :     *    ' HIGHEST LEGAL ATOM NUMBER IS',I4)
*: input elsewhere now : 9148 FORMAT(1X,'YOU MUST INPUT 3N-',I1,' INTERNAL COORDS WHERE N ='
*: input elsewhere now :     *      ,I3/' YOU HAVE PROVIDED ',I6)
      END
