      SUBROUTINE vib_PRDATA(rtdb,ARRAY,VECT,NR,NC,NRP,NCP,PVEC,NPRI
     &    ,first_pass)
* $Id$
C*MODULE VIB  *DECK PRDATA
c
      IMPLICIT NONE
#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
c
      integer rtdb, nr, nc, nrp, ncp, npri
      double precision ARRAY(NR,NC), VECT(NC)
      CHARACTER*11 FREQ
      character*255 filename
      LOGICAL PVEC
      logical first_pass
      logical projected
c
      integer i, j, l
      integer ncl, ncu, nloops
      integer fd
      double precision offset, nrd(2)
c
      FREQ = 'Frequency'
      projected = .not.first_pass
      if (projected) FREQ = 'P.Frequency'
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
        call ecce_print1('frequencies',mt_dbl,vect,nc)
        if(.not. rtdb_put(rtdb,'vib:frequencies',MT_DBL,nc,
     &      vect))
     &      call errquit('vib_prdata:rtdb_put of frequencies failed',
     &      555, RTDB_ERR)
        call ecce_print2('normal modes',mt_dbl,
     &      array,nr,nr,nc)
      else
        call ecce_print1('projected frequencies',mt_dbl,vect,nc)
        if(.not. rtdb_put(rtdb,'vib:projected frequencies',
     &      MT_DBL,nc,vect))
     &      call errquit(
     &      'vib_prdata:rtdb_put of projected frequencies failed',
     &      555, RTDB_ERR)
        call ecce_print2('projected normal modes',mt_dbl,
     &      array,nr,nr,nc)
      endif
c
c print frequencies and normal modes to file in case it is needed for vscf
c
      if ((.not.first_pass).and.(ga_nodeid().eq.0)) then
        call util_file_name('nmode',.false.,.false.,filename)
        if (eaf_open(filename,eaf_rw,fd).ne.0)
     &    call errquit('vib_prdata:eaf_open failed',555, DISK_ERR)
c
c     Format is (with all data as reals)
c
c     nr = no. of rows or variables
c     nc = no. of columns or modes
c     freq(nc) = frequencies
c     evecs(nr,nc) = eigenvectors
c     
c     To avoid having to worry about integer sizes between languages
c     just write the integer dimensions as doubles
c
        nrd(1) = nr
        nrd(2) = nc
c
        offset = 0
        if (eaf_write(fd,offset,nrd,16).ne.0)
     &    call errquit('vib_prdata:error writing nmodes',555, DISK_ERR)
        offset = 16
        if (eaf_write(fd,offset,vect,nr*8).ne.0)
     &    call errquit('vib_prdata:error writing frequencies',555,
     &       DISK_ERR)
        offset = offset + nr*8
        if (eaf_write(fd,offset,array,nr*nc*8).ne.0)
     &    call errquit('vib_prdata:error writing modes',555, DISK_ERR)
        if (eaf_close(fd).ne.0)
     &    call errquit('vib_prdata:file close problem',555, DISK_ERR)
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
