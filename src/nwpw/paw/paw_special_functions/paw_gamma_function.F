!  *************************************************************************
!  
!        This is a slightly modified version of Log(Gamma) function
!        program from Num. Rec.
!       Serves as backup if intrinic Gamma function
!       is not available
!  ************************************************************************
 
      DOUBLE PRECISION FUNCTION PAW_LN_GAMMA (XX)
      implicit none
      DOUBLE PRECISION XX
      INTEGER J
      DOUBLE PRECISION SER
      DOUBLE PRECISION STP
      DOUBLE PRECISION TMP
      DOUBLE PRECISION X
      DOUBLE PRECISION Y, COF(6)
      SAVE STP, COF
      DATA COF, STP/ 76.18009172947146D0, -86.50532032941677D0, 
     1   24.01409824083091D0, -1.231739572450155D0, 
     2   0.1208650973866179D-2, -.5395239384953D-5, 2.5066282746310005D0
     3   / 
      X = XX
      Y = X
      TMP = X + 5.5D0
      TMP = (X + 0.5D0)*DLOG(TMP) - TMP
      SER = 1.000000000190015D0
      DO J = 1, 6
         Y = Y + 1.0D0
         SER = SER + COF(J)/Y
      END DO
 
      PAW_LN_GAMMA = TMP + DLOG(STP*SER/X)
 
      RETURN 
      END 
 
!  *************************************************
!  
!     Name    :
!  
!  
!     Purpose :
!  
!  
!     Created :
!  
!  *************************************************
      DOUBLE PRECISION FUNCTION PAW_GAMMA(X)
      implicit none
      DOUBLE PRECISION X
 
      DOUBLE PRECISION XX
      DOUBLE PRECISION PAW_LN_GAMMA
 
      XX = X
      PAW_GAMMA = DEXP(PAW_LN_GAMMA(XX))
 
      return
      END 
 
 
 
!  *************************************************
!  
!     Name    :
!  
!  
!     Purpose :
!  
!  
!     Created :
!  
!  *************************************************
      DOUBLE PRECISION FUNCTION PAW_GAMMP (A, X)
      implicit none
      DOUBLE PRECISION A, X
 
      DOUBLE PRECISION GAMMCF, GAMSER, GLN
 
      IF (X .LT. A+1.0D0) THEN
 
         CALL PAW_GSER(GAMSER, A, X, GLN)
         PAW_GAMMP = GAMSER
 
      ELSE
 
         CALL PAW_GCF(GAMMCF, A, X, GLN)
         PAW_GAMMP = 1.0D0 - GAMMCF
 
      ENDIF
 
      return
      END 
 
!  *************************************************
!  
!     Name    :
!  
!  
!     Purpose :
!  
!  
!     Created :
!  
!  *************************************************
      SUBROUTINE PAW_GCF(GAMMCF, A, X, GLN)
      implicit none
      INTEGER ITMAX
      DOUBLE PRECISION A, GAMMCF, GLN, X, EPS, FPMIN
      PARAMETER (ITMAX = 100, EPS = 3.D-16, FPMIN = 1.D-30)
      DOUBLE PRECISION AN, B, C, D, DEL, H
      INTEGER I

      DOUBLE PRECISION PAW_LN_GAMMA
      EXTERNAL         PAW_LN_GAMMA
 
      GLN = PAW_LN_GAMMA(A)
      B = X + 1.0D0 - A
      C = 1.0D0/1.D-30
      D = 1.0D0/B
      H = D
      DO I = 1, 100
         AN = -I*(I - A)
         B = B + 2.0D0
         D = AN*D + B
         IF (DABS(D) .LT. 1.D-30) D = 1.D-30
         C = B + AN/C
         IF (DABS(C) .LT. 1.D-30) C = 1.D-30
         D = 1.0D0/D
         DEL = D*C
         H = H*DEL
         IF (DABS(DEL - 1.0D0) .LT. 3.D-16) GO TO 1
      END DO
      PAUSE 'a too large, ITMAX too small in gcf'
    1 CONTINUE
      GAMMCF = DEXP((-X) + A*DLOG(X) - GLN)*H
 
      return
      END 
 
!  *************************************************
!  
!     Name    :
!  
!  
!     Purpose :
!  
!  
!     Created :
!  
!  *************************************************
      SUBROUTINE PAW_GSER(GAMSER, A, X, GLN)
      implicit none
      DOUBLE PRECISION A, X
      DOUBLE PRECISION GAMSER, GLN
 
!    *** local variables ***
      INTEGER ITMAX
      PARAMETER (ITMAX = 100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 3.0D-16)
      INTEGER N
      DOUBLE PRECISION AP, DEL, SUM

      DOUBLE PRECISION PAW_LN_GAMMA
      EXTERNAL         PAW_LN_GAMMA
 
 
      GLN = PAW_LN_GAMMA(A)
 
      IF (X .LE. 0.0D0) THEN
      IF(X.lt.0.0d0) CALL errquit("x < 0 in PAW_GSER",0,1)
         GAMSER = 0.0D0
         RETURN 
      ENDIF
 
      AP = A
      SUM = 1.0D0/A
      DEL = SUM
      DO N = 1, 100
         AP = AP + 1.0D0
         DEL = DEL*X/AP
         SUM = SUM + DEL
         IF (DABS(DEL) .LT. DABS(SUM)*3.0D-16) GO TO 1
 
      END DO
 
      CALL errquit
     >     ("a too large,ITMAX too small in PAW_GSER",0,1)
 
    1 CONTINUE
      GAMSER = SUM*DEXP((-X) + A*DLOG(X) - GLN)
 
      return 
      END 
 
 
 
c $Id$
