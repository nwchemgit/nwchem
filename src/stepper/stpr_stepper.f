      SUBROUTINE stpr_stepper(A,ILEFT,G,C,ETOT,NATD,CONVGE,CONVGG,
     &    converged)
c $Id: stpr_stepper.f,v 1.1 1994-06-24 20:42:46 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      logical converged
      DIMENSION G(NATD*3),A(ILEFT)
      DIMENSION C(3,NATD)
C
C     Calculate number of 64 bit words that core will need for
C     the calculation.
C
      CALL stpr_cneed (INEED, NATD)
C
C     Convergence tolerances.
C
      if(convgg.eq.0.0d0)CONVGG=1.D-5
      if(convge.eq.0.0d0)CONVGE=1.D-8
      IF (INEED.GT.ILEFT) THEN
        WRITE(6,*)'  STEPPER| NOT ENOUGH MEMORY NEED ',INEED-ILEFT,
     *         ' MORE WORDS!'
        stop
      ELSE
        CALL stpr_stepcor(A,INEED,G,C,ETOT,NATD,CONVGE,CONVGG,converged)
      ENDIF
      RETURN
      END
