      SUBROUTINE vib_WRTFREQ(EIGVAL,NAT3,ZEROPE,NPRI)
* $Id: vib_wrtfreq.F,v 1.1 1995-07-29 23:19:49 d3e129 Exp $
C
C This routine scales the eigenvalues properly and then converts them
C      to wavenumbers
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ZEROPE
      PARAMETER (MAXFRQ = 200)
      PARAMETER(D1=1.0D00,D2=2.0D00)
      COMMON /SETCON/ AMS(36),WAVE,NUMANS      ! setup parameters
      DIMENSION EIGVAL(NAT3)    ! eigenvalues
      DIMENSION PRFREQ(MAXFRQ)    ! freqencies to be written to tape10
*rak:      ITAP96 = 96
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
*rak: where does this come from??
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
      ENDIF
10002 FORMAT(1X,//,'THE ZERO-POINT ENERGY (KCAL/MOL) = ', D15.8)
*rak:      IF (NPRI.EQ.0) WRITE (ITAP96,10002)ZPE
*rak:      GO TO 00096
*rak:      ELSE
*rak:      END IF
C
00096 RETURN                    ! go home
      END
