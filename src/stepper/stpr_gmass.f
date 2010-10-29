      SUBROUTINE stpr_gmass(GRAD,ATMASS)
c $Id$
C
C     This routine mass weights the gradient vector.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION GRAD(NAT3)
      DIMENSION ATMASS(NATOM)
      IATOM(I)  = (I+2)/3
C
      DO 00100 II = 1,NAT3
         IATII = IATOM(II)
         FACT = SQRT(ATMASS(IATII))
         GRAD(II) = GRAD(II)/FACT
00100 CONTINUE
      WRITE(6,9000)
      CALL stpr_matout(GRAD,NAT3,1,NAT3)
c      WRITE(6,*)(GRAD(I),I=1,NAT3)
 9000 FORMAT(/,10X,52('-'),
     &  /,10X,'Mass-weighted gradient (atomic units)',
     &  /,10X,52('-'),//)
      RETURN
      END
