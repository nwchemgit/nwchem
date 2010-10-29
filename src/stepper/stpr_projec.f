      SUBROUTINE stpr_projec(HESS,P,HMP)
c $Id$
C
C     This routine constructs HMP = (1-P)H(1-P)
C     where H is the mass weighted hessian.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION HESS(NAT3TR),P(NAT3TR),HMP(NAT3TR)
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
      DO 10 I=1,NAT3
        DO 20 J=1,I
          HMP(ISYM2(I,J))=HESS(ISYM2(I,J))
          DO 30 K=1,NAT3
            HMP(ISYM2(I,J))=HMP(ISYM2(I,J))-HESS(ISYM2(I,K))*
     &      P(ISYM2(K,J))-P(ISYM2(I,K))*HESS(ISYM2(K,J))
            DO 40 L=1,NAT3
              HMP(ISYM2(I,J))=HMP(ISYM2(I,J))+P(ISYM2(I,K))*
     &        HESS(ISYM2(K,L))*P(ISYM2(L,J))
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
