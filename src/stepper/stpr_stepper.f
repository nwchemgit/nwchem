      SUBROUTINE stpr_stepper(A,ILEFT,G,C,ETOT,NATD,CONVGE,CONVGG,
     &    converged, rtdb)
c $Id: stpr_stepper.f,v 1.2 1994-07-28 15:52:03 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      integer rtdb
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
      IF (INEED.GT.ILEFT) THEN
        WRITE(6,*)'  STEPPER| NOT ENOUGH MEMORY NEED ',INEED-ILEFT,
     *         ' MORE WORDS!'
        stop
      ELSE
        CALL stpr_stepcor(A,INEED,G,C,ETOT,NATD,CONVGE,CONVGG,converged,
     &      rtdb)
      ENDIF
      RETURN
      END
