      SUBROUTINE drdy_DSDT (N, G0, XLAM, TS, F)
      implicit real*8 (a-h,o-z)
      DIMENSION G0(N), XLAM(N)
C Evaluate ds/dt for Page-McIver gradient following algorithm
      F = 0.0d00
      DO 10 I = 1,N
         T = G0(I)*EXP(-XLAM(I)*TS)
         F = F + T*T
   10 CONTINUE
      F = SQRT(F)
      RETURN
      END
