      subroutine vib_wmatrix(A,MR,MC,NR,NC)
* $Id: vib_wmatrix.F,v 1.1 1995-07-29 23:19:47 d3e129 Exp $
C*****************************************************************
C                   SUBROUTINE WMATRIX
C*****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(MR,MC)
      lfnpr=6
      NCL=1
      NCU=10
      NLOOPS=NC/10+1
      DO 20 L=1,NLOOPS
      IF (NCU.GT.NC)NCU=NC
      WRITE(LFNPR,1100)(J,J=NCL,NCU)
      DO 10 I=1,NR
   10 WRITE(LFNPR,1200)I,(A(I,J),J=NCL,NCU)
      IF(NCU.GE.NC)RETURN
      NCL=NCU+1
   20 NCU=NCU+10
      RETURN
 1100 FORMAT(/11X,10(I3,9X))
 1200 FORMAT(1X,I3,10F12.5)
      END
