C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_fndmde(STEP,SELECT,VC,NINTER,NAT3,NMODE)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL FIND
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
#include "cstprfiles.fh"
      DIMENSION STEP(NAT3),VC(NAT3,NINTER),SELECT(NAT3)
C
C     If first step determine intial tracking vector.
C
      IF(IWCTR.EQ.1)THEN
C
C       Find norm of SELECT.  If norm is zero then use the user
C       specified vector NMODE.  If norm is nonzero then normalize
C       and find the hessian eigenvector (having non-zero
C       corresponding gradient) with greatest overlap.
C
        SLNGTH = 0.0D0
        DO 10 I = 1,NAT3
          SLNGTH = SLNGTH + SELECT(I)**2
   10   CONTINUE
        SLNGTH = SQRT(SLNGTH)
        IF(SLNGTH.LT.1.D-12)THEN
          FIND = .FALSE.
        ELSE
          FIND = .TRUE.
          DO 20 I = 1,NAT3
            SELECT(I) = SELECT(I)/SLNGTH
   20     CONTINUE
          WRITE(6,*)'  Normalized SELECT: '
c          CALL stpr_matout(SELECT,3,NATOM,3)
          WRITE(6,1001)
          DO 30 I = 1,NATOM
            JSTART = 3*(I-1) + 1
            JEND = JSTART + 2
            WRITE(6,1002)I,(SELECT(J),J=JSTART,JEND)
   30     CONTINUE
        ENDIF
        IF(FIND)THEN
C
C         Use the VC eigenvectors (which are in the cartesian (mw)
C         basis) and dot the SELECT
C         vector into these to determine maximum overlap ... this
C         will define the vector being tracked.
C
          WRITE(6,*)'  Tracking; overlaps of Hessian',
     &          '  eigenvectors and "SELECT" vector. '
          OVRLPM = 0.0D0
          DO 60 J = 1,NINTER
            OVRLP = 0.0D0
            DO 50 K = 1,NAT3
              OVRLP = OVRLP + VC(K,J)*SELECT(K)
   50       CONTINUE
            OVRLP = DABS(OVRLP)
            IF(OVRLP.GT.OVRLPM)THEN
              OVRLPM = OVRLP
              MAXJ = J
            ENDIF
            WRITE(6,1000)J,OVRLP
   60     CONTINUE
          WRITE(6,1010)MAXJ
          NMODE = MAXJ
        ENDIF
      ELSE
C
C       Read previous tracked mode and find current hessian
C       eigenvector with maximum overlap.
C
        WRITE(6,*)'  Tracking; overlaps of Hessian',
     &            '  eigenvectors and previous mode. '
#if defined(OPENCLOSE)
        OPEN (UNIT=ITAP44, FORM='UNFORMATTED', FILE=FILE44,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        REWIND (UNIT=ITAP44)
#else
        REWIND ITAP44
#endif
        READ(ITAP44)(STEP(I),I=1,NAT3)
        OVRLPM = 0.0D0
        DO 80 I = 1,NINTER
          OVRLP = 0.0D0
          DO 70 K = 1,NAT3
            OVRLP = OVRLP + STEP(K)*VC(K,I)
   70     CONTINUE
          OVRLP = DABS(OVRLP)
          IF(OVRLP.GT.OVRLPM)THEN
            OVRLPM = OVRLP
            MAXI = I
          ENDIF
          WRITE(6,1000)I,OVRLP
   80   CONTINUE
        WRITE(6,1010)MAXI
        NMODE = MAXI
#if defined(OPENCLOSE)
        close (unit=itap44,status='keep')
#endif
      ENDIF
C
C     Write this chosen vector to disk.
C
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP44, FORM='UNFORMATTED', FILE=FILE44,
     &    ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP44)
      WRITE(ITAP44)(VC(I,NMODE),I=1,NAT3)
      close (unit=itap44,status='keep')
#else
      REWIND ITAP44
      WRITE(ITAP44)(VC(I,NMODE),I=1,NAT3)
#endif
      WRITE(6,*)'  Eigenvector : ',NMODE,' written to disk. '
      WRITE(6,1001)
      DO 90 I = 1,NATOM
        JSTART = 3*(I-1) + 1
        JEND = JSTART + 2
        WRITE(6,1002)I,(VC(J,NMODE),J=JSTART,JEND)
   90 CONTINUE
c      CALL stpr_matout(VC(1,NMODE),3,NATOM,3)
      RETURN
 1001 FORMAT('  ATOM:',12x,'X',14x,'Y',14x,'Z')
 1002 FORMAT('  ',I3,6X,3(1X,F14.10))
 1000 FORMAT('  eigenvector: ',I2,' has overlap: ',F12.7)
 1010 FORMAT('  eigenvector: ',I2,' defines NMODE ')
      END
C> @}
