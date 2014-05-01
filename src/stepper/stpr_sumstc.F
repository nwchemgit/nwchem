C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_sumstc(STEP,COORD,ATMASS,CMASS,TENIN,CNVGRD)
c     $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL CNVGRD
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
#include "util.fh"
#include "chow.fh"
      DIMENSION STEP(NAT3),COORD(3,NATOM),ATMASS(NATOM),CMASS(3),
     &     TENIN(3,3)
      IF(.NOT.CNVGRD)THEN
C     
C     Add STEP to COORD.
C     
         IJ = 0
         DO 20 I = 1, NATOM
            DO 10 J = 1,3
               IJ = IJ + 1
               COORD(J,I) = COORD(J,I) + STEP(IJ)
 10         CONTINUE
 20      CONTINUE
      ENDIF
C     
C     Write new coordinates.
C     
*rak:      if (util_print('new coordinates',print_low)) then
*rak:         IF(CNVGRD)THEN
*rak:            WRITE(6,1003)
*rak:         ELSE
*rak:            WRITE(6,1000)
*rak:         ENDIF
*rak:         WRITE(6,1001)
*rak:         DO 30 I = 1,NATOM
*rak:            WRITE(6,1002)I,(COORD(J,I),J=1,3)
*rak: 30      CONTINUE
*rak:         write(6,*)
*rak:      endif
*rak:      if (util_print('distances',print_default)) then
*rak:         CALL DISTAN(NATOM,COORD)
*rak:         write(6,*)
*rak:      endif
*rak:      if (util_print('angles',print_default)) then
*rak:         CALL ANGLE(NATOM,COORD)
*rak:         write(6,*)
*rak:      endif
C     
C     Calculate the vector of center of mass and the inertia
C     tensor for the new geometry.
C     
      IF (MWC) THEN
         CALL stpr_mgeom(COORD,ATMASS,CMASS,TENIN)
C     
C     Shift the mass center to the origin of the coordinate system.
C     
         DO 50 I=1,3
            IF (ABS(CMASS(I)).LT.1.0D-06) GOTO 50
            DO 60 J=1,NATOM
               COORD(I,J)=COORD(I,J)-CMASS(I)
 60         CONTINUE
 50      CONTINUE
         WRITE(6,1004)
         WRITE(6,1001)
         DO 70 I = 1,NATOM
            WRITE(6,1002)I,(COORD(J,I),J=1,3)
 70      CONTINUE
         write(6,*)
      ENDIF
      RETURN
 1000 FORMAT('  ',19('-'),' New coordinates ',19('-'))
 1001 FORMAT('  ATOM:',12x,'X',14x,'Y',14x,'Z')
 1002 FORMAT('  ',I3,6X,3(1X,F14.10))
 1003 FORMAT('  ',16('-'),' Converged coordinates ',16('-'))
 1004 FORMAT('  ',19('-'),' Shifted coordinates ',19('-'))
      END
C> @}
