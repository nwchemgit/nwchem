      SUBROUTINE vib_PRDATA(ARRAY,VECT,NR,NC,NRP,NCP,PVEC,NPRI)
* $Id: vib_prdata.F,v 1.2 1995-10-21 09:16:16 d3e129 Exp $
C*MODULE VIB  *DECK PRDATA
c
      IMPLICIT NONE ! DOUBLE PRECISION(A-H,O-Z)
#include "stdio.fh"
c
      integer nr, nc, nrp, ncp, npri
      double precision ARRAY(NR,NC), VECT(NR)
      CHARACTER*9 FREQ
      LOGICAL PVEC
c
      integer i, j, l
      integer ncl, ncu, nloops
c
      FREQ = 'Frequency'
*      ITAP96=96
C
C
      NCL=1
*      NCU=10
*      NLOOPS=NCP/10+1
      NCU=6
      NLOOPS=NCP/6+1
      DO 20 L=1,NLOOPS
        IF (NCU.GT.NCP)NCU=NCP
C
           WRITE(luout,9110)(J,J=NCL,NCU)
*           IF (NPRI.EQ.0) WRITE(ITAP96,9110)(J,J=NCL,NCU)
           IF(PVEC) THEN
             WRITE(luout,9210)FREQ,(VECT(J),J=NCL,NCU)
             IF (NPRI.EQ.0) THEN
*               WRITE(ITAP96,9210)FREQ,(VECT(J),J=NCL,NCU)
*               WRITE(ITAP96,9300)
             ENDIF
           ENDIF
           WRITE(luout,9300)
           DO 12 I=1,NRP
              WRITE(luout,9402)I,(ARRAY(I,J),J=NCL,NCU)
              IF(NPRI.EQ.0) THEN
*                WRITE (ITAP96,9402) I,(ARRAY(I,J),J=NCL,NCU)
              ENDIF
   12      CONTINUE
        IF(NCU.GE.NCP)RETURN
        NCL=NCU+1
*        NCU=NCU+10
        NCU=NCU+6
   20 CONTINUE
C
      RETURN
C
* 9110 FORMAT(/18X,10(I3,9X))
* 9210 FORMAT(/,1X,A9,2X,10F12.2)
 9110 FORMAT(/18X,6(I3,9X))
 9210 FORMAT(/,1X,A9,2X,6F12.2)
 9300 FORMAT(' ')
* 9402 FORMAT(9X,I3,10F12.5)
 9402 FORMAT(9X,I3,6F12.5)
C
      END
