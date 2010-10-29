      subroutine vib_wmatrix(A,MR,MC,NR,NC)
* $Id$
C*****************************************************************
C                   SUBROUTINE WMATRIX
C*****************************************************************
      IMPLICIT none
c
#include "stdio.fh"
c
      integer mr, mc, nr, nc
      double precision A(MR,MC)
c
      integer i, j, l
      integer ncl, ncu, nloops
c
      NCL=1
      NCU=10
      NLOOPS=NC/10+1
      DO 20 L=1,NLOOPS
      IF (NCU.GT.NC)NCU=NC
      WRITE(LUOUT,1100)(J,J=NCL,NCU)
      DO 10 I=1,NR
   10 WRITE(LUOUT,1200)I,(A(I,J),J=NCL,NCU)
      IF(NCU.GE.NC)RETURN
      NCL=NCU+1
   20 NCU=NCU+10
      RETURN
 1100 FORMAT(/11X,10(I3,9X))
 1200 FORMAT(1X,I3,10F12.5)
      END
