      SUBROUTINE vib_PRDATA(ARRAY,VECT,NR,NC,NRP,NCP,PVEC,NPRI
     &    ,first_pass)
* $Id: vib_prdata.F,v 1.6 2000-03-23 22:51:50 edo Exp $
C*MODULE VIB  *DECK PRDATA
c
      IMPLICIT NONE ! DOUBLE PRECISION(A-H,O-Z)
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
c
      integer nr, nc, nrp, ncp, npri
      double precision ARRAY(NR,NC), VECT(NR)
      CHARACTER*11 FREQ
      character*255 filename
      LOGICAL PVEC
      logical first_pass
      logical projected
c
      integer i, j, l
      integer ncl, ncu, nloops
      integer fd
      double precision offset
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
c
c print frequencies and normal modes to file in case it is needed for vscf
c
      if ((.not.first_pass).and.(ga_nodeid().eq.0)) then
        call util_file_name('nmode',.false.,.false.,filename)
        if (eaf_open(filename,eaf_rw,fd).ne.0)
     &    call errquit('vib_prdata:eaf_open failed',555)
        offset = 0
        if (eaf_write(fd,offset,vect,nr).ne.0)
     &    call errquit('vib_prdata:error writing frequencies',555)
        offset = offset + nr
        if (eaf_write(fd,offset,array,nr*nc).ne.0)
     &    call errquit('vib_prdata:error writing modes',555)
        if (eaf_close(fd).ne.0)
     &    call errquit('vib_prdata:file close problem',555)
      endif
c
      RETURN
C
 9110 FORMAT(/18X,6(I3,9X))
 9111 format(1x,'Projected')
 9210 FORMAT(1X,A11,6F12.2)
 9300 FORMAT(' ')
 9402 FORMAT(9X,I3,6F12.5)
C
      END
