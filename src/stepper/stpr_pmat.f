      SUBROUTINE stpr_pmat(VC,P,NEXTER)
c $Id$
C
C     This routine constructs P matrix.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION VC(NAT3,NAT3),P(NAT3TR)
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
      DO 10 I=1,NAT3
        DO 20 J=1,I
          IDUM = ISYM2(I,J)
          P(IDUM)=0.D0
          DO 30 K=1,NEXTER
            P(IDUM)=P(IDUM)+VC(I,K)*VC(J,K)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
