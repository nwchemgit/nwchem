      SUBROUTINE vib_PRDATA(ARRAY,VECT,NR,NC,NRP,NCP,PVEC,NPRI
     &    ,first_pass)
* $Id: vib_prdata.F,v 1.4 1997-04-22 20:01:43 d3e129 Exp $
C*MODULE VIB  *DECK PRDATA
c
      IMPLICIT NONE ! DOUBLE PRECISION(A-H,O-Z)
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer nr, nc, nrp, ncp, npri
      double precision ARRAY(NR,NC), VECT(NR)
      CHARACTER*11 FREQ
      LOGICAL PVEC
      logical first_pass
      logical projected
c
      integer i, j, l
      integer ncl, ncu, nloops
c
      FREQ = 'Frequency'
      projected = .not.first_pass
      if (projected) FREQ = 'P.Frequency'
C
C
      NCL=1
      NCU=6
      NLOOPS=NCP/6+1
      DO 20 L=1,NLOOPS
        IF (NCU.GT.NCP)NCU=NCP
C
        WRITE(luout,9110)(J,J=NCL,NCU)        
        write(luout,9300)
        IF(PVEC) THEN
          WRITE(luout,9210)FREQ,(VECT(J),J=NCL,NCU)
        ENDIF
        write(luout,9300)
        DO 12 I=1,NRP
          WRITE(luout,9402)I,(ARRAY(I,J),J=NCL,NCU)
 12     CONTINUE
        IF(NCU.GE.NCP) goto 00001
        NCL=NCU+1
        NCU=NCU+6
   20 CONTINUE
C
00001 continue
*
* print ecce information
*      
      if (first_pass) then
        call ecce_print1('frequencies',mt_dbl,vect,nr)
        call ecce_print2('normal modes',mt_dbl,
     &      array,nr,nr,nc)
      else
        call ecce_print1('projected frequencies',mt_dbl,vect,nr)
        call ecce_print2('projected normal modes',mt_dbl,
     &      array,nr,nr,nc)
      endif
      RETURN
C
 9110 FORMAT(/18X,6(I3,9X))
 9111 format(1x,'Projected')
 9210 FORMAT(1X,A11,6F12.2)
 9300 FORMAT(' ')
 9402 FORMAT(9X,I3,6F12.5)
C
      END
