      SUBROUTINE stpr_mgeom(COORD,ATMASS,CMASS,TENIN)
c $Id: stpr_mgeom.f,v 1.1 1994-06-24 20:42:31 d3e129 Exp $
C
C     This routine calculates vector of the center of mass
C     and tensor of inertia.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / MASS / TOTM, NUMAS
      DIMENSION COORD(3,NATOM),ATMASS(NATOM),CMASS(3),TENIN(3,3)
      DIMENSION SHIT(3)
      DO 10 I=1,3
        CMASS(I)=0.D0
        DO 20 J=1,NATOM
           CMASS(I)=CMASS(I)+ATMASS(J)*COORD(I,J)
   20   CONTINUE
        CMASS(I) = CMASS(I)/TOTM
   10 CONTINUE
      WRITE(6,*)' STEPPER| Vector of the center of mass:'
      WRITE(6,*) (CMASS(I),I=1,3)
      DO 30 I=1,3
        DO 40 J=1,2
          TENIN(I,J)=0.D0
          DO 50 K=1,NATOM
            TENIN(I,J)=TENIN(I,J)-ATMASS(K)*COORD(I,K)*COORD(J,K)
   50     CONTINUE
          TENIN(J,I)=TENIN(I,J)
   40   CONTINUE
   30 CONTINUE
      DO 60 I=1,3
        SHIT(I)=0.D0
        DO 70 J=1,NATOM
          SHIT(I)=SHIT(I)+ATMASS(J)*COORD(I,J)**2
   70   CONTINUE
   60 CONTINUE
      TENIN(1,1)=SHIT(2)+SHIT(3)
      TENIN(2,2)=SHIT(1)+SHIT(3)
      TENIN(3,3)=SHIT(2)+SHIT(1)
      WRITE(6,*)' STEPPER| Moment of inertia'
      DO 80 I=1,3
        WRITE(6,*)(TENIN(I,J),J=1,3)
   80 CONTINUE
      RETURN
      END
