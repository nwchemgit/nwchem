      SUBROUTINE stpr_hmass(HESS,ATMASS)
c $Id$
C
C     This routine mass weights and scales the Hessian matrix.
C     The scaling is done to avoid numerical problems in the
C     diagonalization routine.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION HESS(NAT3TR)
      DIMENSION ATMASS(NATOM)
C
C     Set up function for locating i,j elements packed canonically
C     as i.
C
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
      IATOM(I)  = (I+2)/3
C
      DO 00100 II = 1,NAT3
         JJEND = II
         IATII = IATOM(II)
         DO 00100 JJ = 1,JJEND
            IDUM = ISYM2(II,JJ)
            IATJJ = IATOM(JJ)
            FACT = SQRT(ATMASS(IATII)*ATMASS(IATJJ))
            HESS(IDUM) = HESS(IDUM)/FACT
00100 CONTINUE
      WRITE(6,9000)
      CALL stpr_prntpd(HESS,NAT3TR,NAT3,6)

 9000 FORMAT(/,10X,52('-'),
     &  /,10X,'Mass-weighted nuclear hessian (atomic units)',
     &  /,10X,52('-'),//)
      RETURN
      END
