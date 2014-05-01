C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_clnstp(STEP,THRES)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION STEP(NAT3)
      DO 10 I=1,NAT3
        STEPX = STEP(I)/THRES
        ISTEPX = NINT(STEPX)
        STEPX = dble(ISTEPX)*THRES
        STEP(I) = STEPX
   10 CONTINUE
      RETURN
      END
C> @}
