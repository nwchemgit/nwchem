      SUBROUTINE ZMATIN(IZMAT,NELS)
C*MODULE VIB  *DECK ZMATIN
* $Id: vib_zmatin.F,v 1.1 1995-07-29 23:19:51 d3e129 Exp $
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR,NVAROK
C
      PARAMETER (MXATM=50)
      PARAMETER (NNAM=1)
C
      COMMON /INFOA /  ZAN(MXATM),C(3,MXATM)
      COMMON /VBDCMP/ NROW,NDEG
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,LINEAR
      COMMON /ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD
      COMMON /HESS/ NAT,NAT3,NHESS,NHESST   ! Hessian information
C
      DIMENSION IZMAT(NELS)
      DIMENSION QNAM(NNAM),KQNAM(NNAM),SZMAT(1)
C
      DATA SZMAT/8HZMAT    /
      DATA QNAM /8HIZMAT   /
C
C     ----- THIS ROUTINE READS THE CODED Z-MATRIX -----
C     IN COMMON /ZMAT/, NZVAR IS THE NUMBER OF INPUT
C     COORDINATES, NVAR IS THE TOTAL NUMBER OF COORDINATES,
C     AND NZMAT IS THE TOTAL NUMBER OF CODE INTEGERS INPUT.
C
C----------Read NVAR, number of internal coordinates from
C          standard input.
C
*rak flags from input
*rak:      LINEAR = .TRUE.
*rak:      CALL FIND('LINEAR',IERR)
*rak:      IF(IERR .NE. 0)LINEAR = .FALSE.
*rak:      CALL FIND('NVAR  ',IERR)
*rak:      IF(IERR .NE. 0)THEN
*rak:        WRITE(6,*)'Did not find flag NVAR.  Call Bomb.'
*rak:        CALL BOMB(0)
*rak:      ENDIF
*rak:      READ(5,*)NZVAR
*rak:C
*rak:      KQNAM(1) = 10*7*NZVAR + 1
*rak:      DO 100 I=1,NELS
*rak:         IZMAT(I) = 0
*rak:  100 CONTINUE
*rak:C
*rak:C----------READ IN DESCRIPTION OF INTERNAL COORDINATES.
*rak:C
*rak:       CALL FIND('INTCRD',IERR)
*rak:       IF(IERR.NE.0)THEN
*rak:         WRITE(6,*)'Could not find flag INTCRD.  Call Bomb.'
*rak:         CALL BOMB(0)
*rak:       ENDIF
*rak:       IPT = 1
*rak:       DO 7117 IVAR = 1,NZVAR
*rak:         READ(5,*)IZMAT(IPT)
*rak:         ITYPE = IZMAT(IPT)
*rak:         IF(ITYPE .EQ. 1)THEN
*rak:          IPT2 = IPT + 2
*rak:         ELSEIF(ITYPE .EQ. 2)THEN
*rak:          IPT2 = IPT + 3
*rak:         ELSEIF(ITYPE .EQ. 3)THEN
*rak:          IPT2 = IPT + 4
*rak:         ELSEIF(ITYPE .EQ. 4)THEN
*rak:          IPT2 = IPT + 4
*rak:         ELSEIF(ITYPE .EQ. 5)THEN
*rak:          IPT2 = IPT + 3
*rak:         ELSEIF(ITYPE .EQ. 6)THEN
*rak:          IPT2 = IPT + 5
*rak:         ELSEIF(ITYPE .EQ. 7)THEN
*rak:          IPT2 = IPT + 6
*rak:         ELSE
*rak:          WRITE(6,9001)ITYPE
*rak:          CALL BOMB(0)
*rak:         ENDIF
*rak:         READ(5,*)(IZMAT(I),I=IPT+1,IPT2)
*rak:         IPT = IPT2 + 1
*rak: 7117  CONTINUE
*rak:C
*rak:C         INPUT OF INTERNALS REQUIRED IF COORD=UNIQUE, HINT, OR CART,
*rak:C         BUT FOR COORD=ZMT, ZMTMPC WE CAN USE $DATA INFO INSTEAD
*rak:C
*rak:C
*rak:C     ----- SCAN AND ANALYZE IZMAT -----
*rak:C     THE VALID TYPES ARE
*rak:C        1= STRETCH I-J
*rak:C        2= BEND I-J-K
*rak:C        3= TORSION I-J-K-L
*rak:C           (BETWEEN TWO PLANES HAVING TWO ATOMS IN COMMON)
*rak:C        4= OUT OF PLANE BEND I-J-K-L
*rak:C        5= LINEAR BEND I-J-K
*rak:C        6= 5 ATOM DIHEDRAL ANGLE I-J-K-L-M
*rak:C           (BETWEEN TWO PLANES HAVING ONE ATOM IN COMMON)
*rak:C        7= 6 ATOM GHOST TORSION ANGLE I1-I2-J-K-L1-L2
*rak:C
*rak:  200 CONTINUE
*rak:c      WRITE(6,9048)
*rak:      NVAR=0
*rak:      IZ = 1
*rak:      DO 230 I=1,NZVAR
*rak:         ITYPE=IZMAT(IZ)
*rak:         NCHK=0
*rak:         IF(ITYPE.EQ.1) NCHK=2
*rak:         IF(ITYPE.EQ.2) NCHK=3
*rak:         IF(ITYPE.EQ.3) NCHK=4
*rak:         IF(ITYPE.EQ.4) NCHK=4
*rak:         IF(ITYPE.EQ.5) NCHK=3
*rak:         IF(ITYPE.EQ.6) NCHK=5
*rak:         IF(ITYPE.EQ.7) NCHK=6
*rak:         IF(NCHK.EQ.0) THEN
*rak:            WRITE(6,9068) I,IZ,IZMAT(IZ)
*rak:            CALL BOMB(0)
*rak:         END IF
*rak:C
*rak:c         WRITE(6,9058) I,ITYPE,(IZMAT(IZ+J),J=1,NCHK)
*rak:         DO 210 J=1,NCHK
*rak:            IAT=IZMAT(IZ+J)
*rak:            IF(IAT.LE.0  .OR.  IAT.GT.NAT) THEN
*rak:               WRITE(6,9088) I,IZ,IAT,J,NAT
*rak:               CALL BOMB(0)
*rak:            END IF
*rak:  210    CONTINUE
*rak:C
*rak:         IZ=IZ+NCHK+1
*rak:         NVAR=NVAR+1
*rak:         IF(ITYPE.EQ.5) NVAR=NVAR+1
*rak:  230 CONTINUE
*rak:      NZMAT = IZ - 1
*rak:C
*rak:C --- IF NOT DOING VIBRATIONAL ANALYSIS, MUST SPECIFY EXACTLY 3N-NDEG
*rak:C --- INTERNAL COORDINATES.  OTHERWISE, MORE THAN 3N-NDEG IS OK.
*rak:C
*rak:      NDEG = 6
*rak:      IF (LINEAR) NDEG = 5
*rak:      NROW = NVAR + NDEG
*rak:      NVAROK = NVAR .EQ. (3*NAT-NDEG)
*rak:      IF (.NOT.NVAROK) THEN
*rak:         WRITE(6,9148) NDEG,NAT,NVAR
*rak:         CALL BOMB(0)
*rak:      END IF
*rak:C
*rak:C     ----- SAVE IZMAT ON TAPE 10 ----
*rak:C
*rak:c       CALL ADDARRAY(10,IZMAT,INTOWP(NELS),393)
*rak:C
*rak:      NVAR2  = NVAR
*rak:      NZMAT2 = NZMAT
*rak:      RETURN
*rak:C
*rak: 9000 FORMAT(1X,'COULD NOT LOCATE U:VIB')
*rak: 9001 FORMAT(1X,'UNKNOWN INTERNAL COORDINATE TYPE IN ',
*rak:     *          'INPUT DECK. ITYPE = ',I5)
*rak: 9008 FORMAT(1X,'ERROR IN $ZMAT INPUT - ERROR')
*rak: 9048 FORMAT(/3X,'--- ENCODED Z MATRIX ---'/
*rak:     *        1X,'COORD',2X,'TYPE',3X,'I',3X,'J',3X,'K',3X,'L',3X,'M',
*rak:     *        3X,'N')
*rak: 9058 FORMAT(1X,I3,3X,I4,1X,6I4)
*rak: 9068 FORMAT(1X,'INTERNAL COORD.',I3,' AT SUBSCRIPT',I4,
*rak:     *     ' IS OF ILLEGAL TYPE',I5)
*rak: 9088 FORMAT(1X,'INTERNAL COORD.',I3,' AT SUBSCRIPT',I4,
*rak:     *    ' IS FOLLOWED BY ILLEGAL ATOM',I4,' AT OFFSET',I4/
*rak:     *    ' HIGHEST LEGAL ATOM NUMBER IS',I4)
*rak: 9148 FORMAT(1X,'YOU MUST INPUT 3N-',I1,' INTERNAL COORDS WHERE N ='
*rak:     *      ,I3/' YOU HAVE PROVIDED ',I6)
      END
