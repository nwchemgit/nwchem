      SUBROUTINE stpr_clnstp(STEP,THRES)
c $Id: stpr_clnstp.F,v 1.1 1994-06-24 20:42:16 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION STEP(NAT3)
      DO 10 I=1,NAT3
        STEPX = STEP(I)/THRES
        ISTEPX = NINT(STEPX)
#if defined(KSR) || defined(CRAY)
        STEPX = real(ISTEPX)*THRES
#else
        STEPX = DBLE(ISTEPX)*THRES
#endif
        STEP(I) = STEPX
   10 CONTINUE
      RETURN
      END
