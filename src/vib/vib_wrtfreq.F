      SUBROUTINE vib_WRTFREQ(rtdb,EIGVAL,NAT3,ZEROPE,NPRI)
* $Id: vib_wrtfreq.F,v 1.6 2001-07-27 21:34:01 windus Exp $
C
C This routine scales the eigenvalues properly and then converts them
C      to wavenumbers
C
      IMPLICIT none
c
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer nat3
      LOGICAL ZEROPE
      integer npri
      double precision EIGVAL(NAT3)    ! eigenvalues
c
      integer maxfrq
      PARAMETER (MAXFRQ = 2000)
      double precision d1, d2
      PARAMETER(D1=1.0D00,D2=2.0D00)
      double precision ams, wave
      integer numans
      COMMON /cvib_SETCON/ AMS(36),WAVE,NUMANS      ! setup parameters
      double precision PRFREQ(MAXFRQ)    ! freqencies to be written to tape10
c
      double precision scale, dd1, zero, sumfreq, rate, xdum, freq, zpe
      integer i
c
*rak:      ITAP96 = 96
c: error check :
      if (nat3.gt.maxfrq) then
c        dimension problems
        write(luout,*)' maxfrq = ',maxfrq
        write(luout,*)' nat3   = ',nat3
        write(luout,*)' increase maxfrq to ',nat3
        call errquit('vib_wrtfreq:',nat3)
        stop ' error'
      endif
C
C Scale eigenvalues down by 1000
C
      SCALE = 1.0D-03   ! Scaling factor after diagonaization
      call dscal(nat3,scale,eigval,1) ! Scale eigenvalues
C
      DD1 = D1          ! used with dsign function
C      WRITE(6,10000)    ! write header
10000 FORMAT(//,1X,31('-'),' Vibrations ',31('-'),/,
     +5X,'#',6X,'Eigenvalue',14X,'Frequency (cm**-1)',/,1X,74('-'))
10001 FORMAT(1X,I5,D20.10,F20.3)
      ZERO=0.0D00
      SUMFREQ=zero
*rak: where does this come from?
      RATE=2.860D-03  
      DO 00100 I=1,NAT3 ! loop over coordinates
         XDUM = EIGVAL(I) ! eigenvalue
         FREQ = SIGN(DD1,XDUM)*SQRT(ABS(XDUM))*WAVE ! form frequency
         PRFREQ(I) = FREQ
C         WRITE(6,10001,ERR=00099)I,XDUM,FREQ ! write it out
         EIGVAL(I) = FREQ
         IF(ZEROPE) THEN
            IF (FREQ.GT.ZERO)THEN
              SUMFREQ=SUMFREQ+FREQ
            ENDIF
         END IF
         GO TO 00100
C00099    WRITE(6,*)EIGVAL(I),FREQ  ! in case of error with format
00100 CONTINUE
      IF(ZEROPE)THEN
        ZPE=0.5D00*SUMFREQ*RATE
        WRITE (6,10002)ZPE
        call ecce_print1('zero point energy',mt_dbl,zpe,1)
        if(.not. rtdb_put(rtdb,'vib:zpe',MT_DBL,1,zpe))
     &      call errquit('vib_wrtfreq:rtdb_put of ZPE failed',
     &      555)
      ENDIF
10002 FORMAT(1X,//,'From the projected analysis ',/,
     &    'The Zero-Point Energy (Kcal/mol) = ',
     &    f20.8)
*rak:      IF (NPRI.EQ.0) WRITE (ITAP96,10002)ZPE
*rak:      GO TO 00096
*rak:      ELSE
*rak:      END IF
C
00096 RETURN                    ! go home
      END
