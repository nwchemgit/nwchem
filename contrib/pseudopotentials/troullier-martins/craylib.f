CFrom lyang@tazdevil.llnl.gov Mon Sep 23 17:44:01 1991
CReturn-Path: <lyang@tazdevil.llnl.gov>
CReceived: from tazdevil.llnl.gov by  icose  (NeXT-1.0 (From Sendmail 5.52)/NeXT-2.0)
C	id AA06758; Mon, 23 Sep 91 17:43:58 CDT
CReceived: by tazdevil.llnl.gov (4.1/1.15)
C	id AA06980; Mon, 23 Sep 91 15:42:04 PDT
CDate: Mon, 23 Sep 91 15:42:04 PDT
CFrom: lyang@tazdevil.llnl.gov (Lin H. Yang)
CMessage-Id: <9109232242.AA06980@tazdevil.llnl.gov>
CTo: asmith@icose.msd.anl.gov
CSubject: library routines
CStatus: R
C
C
CI am including some library routines you will need for non-Cray machines.
C
CMakesure you also need to modify your Makefile to include the library routines.
C
C
C ********* SSCAL
C
      SUBROUTINE SSCAL(N,SCALE,ARRAY,INC)
      INTEGER*4 I,N,INC
      REAL*8 SCALE,ARRAY(*)
      DO 101 I = 1, N, INC
      ARRAY(I) = ARRAY(I) * SCALE
101   CONTINUE
      RETURN
      END
C
C ********* CSSCAL
C
      SUBROUTINE CSSCAL(N,SCALE,CARRAY,INC)
      INTEGER*4 I,N,INC
      REAL*8 SCALE
      DOUBLE COMPLEX CARRAY(*)
      DO 101 I = 1, N, INC
      CARRAY(I) = CARRAY(I) * SCALE
101   CONTINUE
      RETURN
      END
C
C *********** SSUM
C
      FUNCTION SSUM(N,ARRAY,INC)
      INTEGER*4 I,N,INC
      REAL*8 SSUM,ARRAY(*)
      SSUM = 0.0E0
      DO 101 I = 1, N, INC
      SSUM = SSUM + ARRAY(I)
101   CONTINUE
      RETURN
      END
C
C *********** CSUM
C
      FUNCTION CSUM(N,CARRAY,INC)
      INTEGER*4 I,N,INC
      DOUBLE COMPLEX CSUM,CARRAY(*)
      CSUM = CMPLX(0.0E0, 0.0E0)
      DO 101 I = 1, N, INC
      CSUM = CSUM + CARRAY(I)
101   CONTINUE
      RETURN
      END
C
C ********** SDOT
C
      FUNCTION SDOT(N,ARRAY,INC1,BRRAY,INC2)
      INTEGER*4 I,N,INC1,INC2,I1,I2
      REAL*8 SDOT,ARRAY(*),BRRAY(*)
      SDOT = 0.0E0
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      SDOT = SDOT + ARRAY(I1) * BRRAY(I2)
101   CONTINUE
      RETURN
      END
C
C ********** CDOTC
C
      FUNCTION CDOTC(N,CARRAY,INC1,CBRRAY,INC2)
      INTEGER*4 I,N,INC1,INC2,I1,I2
      DOUBLE COMPLEX CDOTC,CARRAY(*),CBRRAY(*)
      CDOTC = CMPLX(0.0E0, 0.0E0)
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      CDOTC = CDOTC + CONJG(CARRAY(I1)) * CBRRAY(I2)
101   CONTINUE
      RETURN
      END
C
C ********** SCNRM2
C
      FUNCTION SCNRM2(N,CARRAY,INC)
      INTEGER*4 I,N,INC
      REAL*8 SCNRM2
      DOUBLE COMPLEX CARRAY(*)
      SCNRM2 = 0.0E0
      DO 101 I = 1, N, INC
      SCNRM2 = SCNRM2 + ABS(CARRAY(I))**2
101   CONTINUE
      SCNRM2 = DSQRT(SCNRM2)
      RETURN
      END
C
C ********** CAXPY
C
      SUBROUTINE CAXPY(N,CNST,CARRAY,INC1,CBRRAY,INC2)
      INTEGER*4 I,N,INC1,INC2,I1,I2
      DOUBLE COMPLEX CNST,CARRAY(*),CBRRAY(*)
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      CBRRAY(I2) = CBRRAY(I2) + CNST * CARRAY(I1)
101   CONTINUE
      RETURN
      END
C
C ********** SCOPY
C
      SUBROUTINE SCOPY(N,ARRAY,INC1,BRRAY,INC2)
      INTEGER*4 I,N,INC1,INC2,I1,I2
      REAL*8 ARRAY(*),BRRAY(*)
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      BRRAY(I2) = ARRAY(I1)
101   CONTINUE
      RETURN
      END
C
C ********** CCOPY
C
      SUBROUTINE CCOPY(N,CARRAY,INC1,CBRRAY,INC2)
      INTEGER*4 I,N,INC1,INC2,I1,I2
      DOUBLE COMPLEX CARRAY(*),CBRRAY(*)
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      CBRRAY(I2) = CARRAY(I1)
101   CONTINUE
      RETURN
      END
C
C ********** ICOPY
C
      SUBROUTINE ICOPY(N,ARRAY,INC1,BRRAY,INC2)
      INTEGER I,N,INC1,INC2,I1,I2
      INTEGER ARRAY(*),BRRAY(*)
      DO 101 I = 1, N
      I1 = 1 + (I-1)*INC1
      I2 = 1 + (I-1)*INC2
      BRRAY(I2) = ARRAY(I1)
101   CONTINUE
      RETURN
      END
C
      SUBROUTINE GATHER(NSIZE,A_AR,B_AR,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A_AR(*),B_AR(*),INDEX(*)
C
      DO 101 IA = 1, NSIZE
      A_AR(IA) = B_AR(INDEX(IA))
101   CONTINUE
      RETURN
      END
C
C ********* SECOND
C
      SUBROUTINE SECOND(ELAPSE)
C
C --- Returns the CPU time (seconds) since beginning of program.
C     Resolution of 0.01 second.
C
      REAL*8 ELAPSE
      REAL TARRAY(2),TIME_T
C
      TIME_T = ETIME(TARRAY)
      TIME_T = TARRAY(1)
      ELAPSE = TIME_T
C
      RETURN
      END
C
C ********* CLOCK
C
      SUBROUTINE CLOCK(TCLOCK)
C
C --- Returns the CLOCK time 
C
      CHARACTER*24 TJOB
      CHARACTER*1 TCLOCK(10),TJOB_T(24)
      EQUIVALENCE (TJOB,TJOB_T)
      INTEGER I, II
C
      CALL FDATE(TJOB)
      DO 101 I = 11, 20
      II = I - 10
      TCLOCK(II) = TJOB_T(I)
101   CONTINUE
C
      RETURN
      END
C
C ********* DATE
C
      SUBROUTINE DATE(TDATE)
C
C --- Returns the CLOCK time 
C
      CHARACTER*24 TJOB
      CHARACTER*1 TDATE(10),TJOB_T(24)
      EQUIVALENCE (TJOB,TJOB_T)
C
      CALL FDATE(TJOB)
      DO 101 I = 5, 10
      II = I - 4
      TDATE(II) = TJOB_T(I)
101   CONTINUE
      TDATE(7) = ' '
      TDATE(8) = ' '
      TDATE(9) = TJOB_T(23)
      TDATE(10) = TJOB_T(24)
C
      RETURN
      END
C
C ********* AIMAG
C
      DOUBLE PRECISION FUNCTION AIMAG(CX)
      DOUBLE COMPLEX CX
      AIMAG = IMAG(CX)
      RETURN
      END


