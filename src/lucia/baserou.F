      SUBROUTINE FRMDSCO(ARRAY,NDIM,MBLOCK,IFILE,IMZERO)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'rou_stat.inc'
      DIMENSION ARRAY(*)
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(IMZERO,1,MBLOCK,IFILE)
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
          CALL SETVEC(ARRAY,ZERO,NDIM)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      IREST=NDIM
      IBASE=0
  100 CONTINUE
       IF(IREST.GT.NBLOCK) THEN
        READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
        IBASE=IBASE+NBLOCK
        IREST=IREST-NBLOCK
        XOP_FRMDSCO = XOP_FRMDSCO + NBLOCK
       ELSE
        READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
        XOP_FRMDSCO = XOP_FRMDSCO + IREST
        IREST=0
       END IF
      IF( IREST .GT. 0 ) GOTO 100
      END IF
 1001 CONTINUE
*
      RETURN
      END
      SUBROUTINE SKPRCD2(NDIM,MBLOCK,IFILE)
C
C     Skip record in file IFILE           
C
*. Version allowing zero and packed blocks 
*
* Dos not work with FASTIO - I expect 
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
 
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
*. Loop over packed records of dimension LPBLK
*. The next LPBLK elements 
  999   CONTINUE
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) 
            READ(IFILE) 
          END IF
          READ(IFILE) ISTOP
        IF(ISTOP.EQ.0) GOTO 999
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) 
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) 
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
C       IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
C        CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
C       END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
      SUBROUTINE FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,
     &                   NO_ZEROING)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks 
*
* If NO_ZEROING = 1, the elements of zero blocks
*    are not set to zero, the routine just returns with 
*    IMZERO = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'rou_stat.inc'
      DIMENSION ARRAY(*)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
*
      NTEST = 0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from FRMDSC2'
        WRITE(6,*) ' IFILE, NDIM, MBLOCK = ', IFILE,NDIM,MBLOCK
      END IF
 
      IMZERO = 0
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          IF(NO_ZEROING.EQ.0) THEN
            ZERO = 0.0D0
            CALL SETVEC(ARRAY,ZERO,NDIM)
          END IF
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements 
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
            XOP_FRMDSC2 = XOP_FRMDSC2 + LBATCH
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT 
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE    
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              STOP ' problem in FRMDSC '
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
          XOP_FRMDSC2 = XOP_FRMDSC2 + NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          XOP_FRMDSC2 = XOP_FRMDSC2 + IREST
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
        IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
         CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
      SUBROUTINE FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks 
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'rou_stat.inc'
      DIMENSION ARRAY(*)
*
      PARAMETER (NTEST=00)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
 
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
C?        write(6,*) ' frmdsc, length of zero block',NDIM
          CALL SETVEC(ARRAY,ZERO,NDIM)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements 
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
            XOP_FRMDSC = XOP_FRMDSC + LBATCH
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT 
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE    
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              STOP ' problem in FRMDSC '
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
          XOP_FRMDSC = XOP_FRMDSC + NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          XOP_FRMDSC = XOP_FRMDSC + IREST
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
        IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
         CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
      SUBROUTINE FRMDSCE
     &     (ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,IERR)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
C     version with error code
C
*. Version allowing zero and packed blocks 
*
      IMPLICIT REAL*8(A-H,O-Z) 
      INCLUDE 'rou_stat.inc'
      DIMENSION ARRAY(*)
*
      PARAMETER (NTEST=00)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
 
C
      IERR = 0
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
        CALL IFRMDSE(ISCR,2,2,IFILE,IERR)
        IF (IERR.NE.0) RETURN
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
          CALL SETVEC(ARRAY,ZERO,NDIM)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements 
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE,END=201,ERR=202) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
            XOP_FRMDSCE = XOP_FRMDSCE + LBATCH
          END IF
          READ(IFILE,END=201,ERR=202) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT 
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE    
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              STOP ' problem in FRMDSC '
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE,END=201,ERR=202) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
          XOP_FRMDSCE = XOP_FRMDSCE + NBLOCK
         ELSE
          READ(IFILE,END=201,ERR=202) (ARRAY(IBASE+I),I=1,IREST)
          XOP_FRMDSCE = XOP_FRMDSCE + IREST
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
 201  IERR = 1  ! end of file
      RETURN
 202  IERR = 2
      RETURN
      END
      SUBROUTINE TODSC(A,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'rou_stat.inc'
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
*
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. No packing 
        ISCR(2) = 0
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
       XOP_TODSC = XOP_TODSC + NTRANS
      IF(NBACK.NE.0) GOTO 100
      END IF
C
      IF( ICRAY.EQ.0.AND.MBLOCK.LT.0.AND.NDIM.GT.0) THEN
       CALL SQFILE(IFIL,1,A,2*NDIM)
      END IF
*
 1001 CONTINUE
C
C?    write(6,*) ' leaving TODSC '
      RETURN
      END
      SUBROUTINE TODSCP(A,NDIM,MBLOCK,IFIL)
*
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
*
* Packed version : Store only nonzero elements
*. Small elements should be zeroed outside
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'rou_stat.inc'
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
* 
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
*
*
C?    write(6,*) ' entering TODSCP, file = ', IFIL
C?    CALL UTIL_FLUSH(6)
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. Packing 
        ISCR(2) = 1
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
*. Loop over packed records of dimension LPBLK
      IELMNT = 0
 1000 CONTINUE
*. The next LPBLK elements 
      LBATCH = 0
*. Obtain next batch of elemnts
  999 CONTINUE
       IF(NDIM.GE.1) THEN
       IELMNT = IELMNT+1
       IF(A(IELMNT).NE.0.0D0) THEN
         LBATCH=LBATCH+1
         IPAK(LBATCH) = IELMNT
         XPAK(LBATCH) = A(IELMNT)
       END IF
       END IF
       IF(LBATCH.EQ.LPBLK.OR.IELMNT.EQ.NDIM) goto 998
       GOTO 999
*. Send to DISC
 998   CONTINUE   
       WRITE(IFIL) LBATCH
       IF(LBATCH.GT.0) THEN 
         WRITE(IFIL) (IPAK(I),I=1, LBATCH)
         WRITE(IFIL) (XPAK(I),I=1, LBATCH)
         XOP_TODSCP = XOP_TODSCP + LBATCH
       END IF
       IF(IELMNT.EQ.NDIM) THEN
         WRITE(IFIL) -1
       ELSE
         WRITE(IFIL) 0
         GOTO 1000
       END IF
*. End of loop over records of truncated elements
      END IF
 1001 CONTINUE
*
C?    CALL UTIL_FLUSH(6)
      RETURN
      END
      SUBROUTINE ADDDIA(A,FACTOR,NDIM,IPACK)
*
* add factor to diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      DO 100 I = 1,NDIM
        IF(IPACK .EQ. 0 ) THEN
          II = (I-1)*NDIM + I
        ELSE
          II = I*(I+1)/2
        END IF
        A(II) = A(II) + FACTOR
  100 CONTINUE
*
      RETURN
      END
        SUBROUTINE BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
C
C       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      DOUBLE PRECISION E,F
C*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(NSIZE,1),EL(NSIZE,1)
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO 30 I=1,N
        DO 10 J=1,N
 10     EL(I,J)=0.0D0
 30     EL(I,I)=1.0D0
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO 50 J=1,N1
        DO 45 I=M,N
        IF(A(I,J).EQ.0.0D0)GO TO 45
        D=DSQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
        C=A(J,J)/D
        S=A(I,J)/D
 38     DO 39 K=J,N
        D=C*A(J,K)+S*A(I,K)
        A(I,K)=C*A(I,K)-S*A(J,K)
        A(J,K)=D
 39     CONTINUE
        DO 40 K=1,N
        D=C*EL(J,K)+S*EL(I,K)
        EL(I,K)=C*EL(I,K)-S*EL(J,K)
        EL(J,K)=D
 40     CONTINUE
 45     CONTINUE
 50     M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
 51     DETERP=A(1,1)
        DO 52 I=2,N
 52     DETERP=DETERP*A(I,I)
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
 520    F=A(1,1)
        E=A(1,1)
        DO 58 I=2,N
        IF(DABS(F).LT.DABS(A(I,I)))F=A(I,I)
        IF(DABS(E).GT.DABS(A(I,I)))E=A(I,I)
 58     CONTINUE
        EPSILP=EPSIL
        IF(EPSILP.LE.0)EPSILP=1.0E-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO 100 J1=1,N
C       CALL SLITE(2)
        I=J
        ISL2=1
        DO 90 I1=1,J
C       CALL SLITET(2,K000FX)
        IF(ISL2.EQ.0)K000FX=2
        IF(ISL2.EQ.1)K000FX=1
        IF(ISL2.EQ.1)ISL2=0
        GO TO (70,75),K000FX
 70     A(I,J)=1.0D0/A(I,I)
        GO TO 90
 75     KS=I+1
        D=0.0D0
        DO 80 K=KS,J
 80     D=D+A(I,K)*A(K,J)
        A(I,J)=-D/A(I,I)
 90     I=I-1
 100    J=J-1
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX
 
C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
 105    M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
 107    D=D+A(I,K)*EL(K,J)
        EL(I,J)=D
 118    CONTINUE
 120    M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
 123    DO 124 I=1,N
        DO 124 J=1,N
 124    A(I,J)=EL(I,J)
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
      INTEGER FUNCTION CANIND(I,J)
C
      IF(I.GT.J) THEN
       CANIND=I*(I-1)/2 + J
      ELSE
       CANIND=J*(J-1)/2 + I
      END IF
      RETURN
      END
      SUBROUTINE CHLFC1(AL,NDIM)
C
C FACTORIZE A SYMMETRIX MATRIX IN AL TO GIVE
C CHOLESKY FACTOR , ALSO IN AL .
C
C INPUT MATRIX AND FACTORIZED MATRIX ARE ASSUMED GIVEN IN
C LOWER TRIANGULAR FORM WITH INDEXING (I,J) = I*(I-1)-2 + J
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(*)
      REAL * 8   INPROD
C
      DO 100 J = 1, NDIM
        JJ = J*(J-1)/2
        AL(JJ+J) = SQRT( AL(JJ+J) - INPROD(AL(JJ+1),AL(JJ+1),J-1) )
        ALJJI = 1.0D0/AL(JJ+J)
        DO 80 I = J+1, NDIM
           II = I*(I-1)/2
           AL(II+J) = (AL(II+J) -
     &                 INPROD( AL(II+1), AL(JJ+1), J-1 ) ) * ALJJI
  80   CONTINUE
 100  CONTINUE
C
      NTEST = 00
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' CHOLESKY FACTORIZATION '
        CALL PRSYM(AL,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE CHLFCB(AL,NDIM,IB,INDEF)
C
C FACTORIZE A SYMMETRIC POSITIVE DEFINITE BAND  MATRIX,AL,TO GIVE
C CHOLESKY FACTOR , ALSO IN AL .
C
C BANDWIDTH IS IB SO 2*IB + 1 ELEMENTS IN EACH ROW ARE NONVANISHING
C ( IN COMPLETE MATRIX )
C
C
C THE MATRIX IS PACKED IN THE FOLLOWING FORM
C         FIRST INDEX J : NONVANISHING COLUMN ELEMENTS FOR ROW NUMBER
C                         CORRESPONDING TO SECOND INDEX
C                         FIRST ELEMENT IS FIRST NONVANISHING ELEMENT
C         SECOND INDEX I : ROW NUMBER
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(IB+1,NDIM)
      REAL * 8   INPROD
C
      NTEST = 00
      INDEF = 0
      LROW = IB + 1
      DO 100 J = 1, NDIM
        KTERMJ = MIN(IB,J-1)
C       WRITE(6,*) ' KTERMJ  ',KTERMJ
        JEFF = MIN(LROW,J)
C       AL(JEFF,J) =
C    &  SQRT(AL(JEFF,J)-INPROD(AL(1,J),AL(1,J),KTERMJ) )
        XXX        =
     &  AL(JEFF,J)-INPROD(AL(1,J),AL(1,J),KTERMJ)
        IF(XXX.LE.0.0D0 ) THEN
          WRITE(6,*) ' NEGATIVE DIAGONAL ELEMENT IN CHLFCB,J = ',J
          WRITE(6,*) ' VALUE ', XXX
          INDEF = 1
          RETURN
        ELSE
          AL(JEFF,J) = SQRT(XXX)
        END IF
C
        ALJJI = 1.0D0/AL(JEFF,J)
C       WRITE(6,*) ' ALJJI ',ALJJI
        IMIN = J+1
        IMAX = MIN(NDIM,J+IB)
C       WRITE(6,*) ' IMIN IMAX       ',IMIN,IMAX
        DO 80 I = IMIN,IMAX
           IABSTR = MAX(1,I-IB)
           JABSTR = MAX(1,J-IB)
           KSTRJ = IABSTR-JABSTR + 1
           KMAX = MIN(IB + 1 - KSTRJ,J-KSTRJ)
           JEFFI = J + 1 - IABSTR
C          WRITE(6,*) ' I IABSTR JABSTR KSTRJ '
C          WRITE(6,*)   I,IABSTR,JABSTR,KSTRJ
C          WRITE(6,*) ' KMAX ,JEFFI ', KMAX,JEFFI
           AL(JEFFI,I) = (AL(JEFFI,I) -
     &                  INPROD( AL(1,I),AL(KSTRJ,J),KMAX ) )*ALJJI
  80   CONTINUE
C       WRITE(6,*) ' CHOLESKY FACTORIZATION AFTER J ',J
C       CALL WRTMAT(AL,IB+1,NDIM,IB+1,NDIM)
 100  CONTINUE
C
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' CHOLESKY FACTORIZATION '
        CALL WRTMAT(AL,IB+1,NDIM,IB+1,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE CHLFCE(AL,NDIM,IB,IALOFF,INDEF)
C
C FACTORIZE A SYMMETRIC POSITIVE DEFINITE ENVELOPE MATRIX,AL,TO GIVE
C CHOLESKY FACTOR , ALSO IN AL .
C
C Matrix AL is stored rowwise in vector AL.
C
C ILOFF(I) Adress in L of first element of row I
C IB(I)    Column number of first row of I
C
C on output L will be stored in the same format
 
C L : matrix stored rowwise in one dimensional array .
C of first nonvaninhing element in row I
C
C Bordering method is used
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(*),IB(*),IALOFF(*)
      REAL * 8  INPROD
C
      NTEST = 0
      INDEF = 0
      DO 100 I = 1, NDIM
        IOFF = IALOFF(I)
        ICSTRT = IB(I)
        NJ = I - ICSTRT
        DO 50 J = ICSTRT, I-1
            JOFF = IALOFF(J)
            JCSTRT = IB(J)
            KMIN =  MAX(ICSTRT,JCSTRT)
            IADDI = KMIN-ICSTRT
            IADDJ = KMIN-JCSTRT
            NK =  J - KMIN
            IJEFF = IOFF + J - ICSTRT
            JJEFF = JOFF + J - JCSTRT
            AL(IJEFF) =
     &      (AL(IJEFF)-INPROD(AL(JOFF+IADDJ),AL(IOFF+IADDJ),NK)) /
     &      AL(JJEFF)
   50   CONTINUE
*
        XXX =  AL(IOFF+I-ICSTRT)-INPROD(AL(IOFF),AL(IOFF),NJ)
        IF(XXX.LE.0.0D0 ) THEN
          WRITE(6,*) ' NEGATIVE DIAGONAL ELEMENT IN CHLFCB,I = ',I
          WRITE(6,*) ' VALUE ', XXX
          INDEF = 1
          RETURN
        ELSE
          AL(IOFF+NJ) = SQRT(XXX)
        END IF
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CLSKHB(AL,X,B,NDIM,IB,ITASK,INDEF)
C
C MASTER ROUTINE FOR SOLVING LINEAR EQUATIONS THROUGH
C CHOLESKY DECOMPOSITION OF POSITIVE DEFINITE BANDED MATRIX  A
C
C THE ACTUAL TASK IS DEFINED THROUGH ITASK
C
C  ITASK = 1 : FACTORIZE MATRIX AND RETURN
C        = 2 : FACTORIZATION HAVE BEEN PERFORMED ( INPUT IN AL )
C              SOLVE LINEAR EQS. MATRIX * X = B
C        = 3 : FACTORIZE AND SOLVE LINEAR EQUATIONS A X = B
C.. INPUT
C
C         AL : ITASK = 1,3 : INPUT MATRIX ( FORMAT : SEE BELOW )
C              OVERWRITTEN !
C              ITASK = 2:  L DECOMPOSITOTATION ASSUMED IN AL )
C              NOT OVERWRITTEN
C
C         X  : VECTOR FOR SOLUTION TO LINEAR EQUATIONS
C         B  : RHS VECTOR FOR LINEAR EQUATIONS( OVERWRITTEN )
C         ( FOR ITASK = 1 X AND B CAN BE DUMMY VARIABLES )
C         NDIM : ORDER OF MATRIX OF MATRICES AND VECTORS
C         IB :  HALF BANDWIDTH, I.E. 2*IB + 1 ELEMENTS IN EACH
C              ROW ARE ASSUMED NONVANISHING
C         ITASK : DEFINING TASK OF ROUTINE AS ABOVE
C
C OUTPUT :
C        ITASK = 1, 3 : AL IS L DECOMPOSITITION , I.E,
C        L IS A LOWER TRIANGULAR POSITIVE MATRIX AND
C        A = L * L ( TRANSPOSED )
C
C        ITASK = 2,3 : X IS SOLUTION TO LINEAR SET OF EQUATIONS
C        INDEF ( FOR ITASK = 1, 3 ) :
C            0 : MATRIX DECOMPOSED IS NOT INDEFINITE
C        .NE.0 : ABNORMAL TERMINATION DUE TO INDEFINITE MATRIX
C
C NOTE ON STRUCTURE OF MATRIX
C
C THE MATRIX IS ASSUMED PACKED SO ONLY LOWER HALF ELEMENTS IN
C THE BAND IS STORED . A IS STORED AS A TWO DIMENSIONAL ARRAY
CWITH
C SECOND INDEX : ROW NUMBER
C FIRST INDEX :  NONVANISHING ELEMENTS FOR THIS ROW, STARTS WITH
C                FIRST NONVANISHING ELEMENT, AND ENDS WITH
C                DIAGONAL ELEMENT.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(IB+1,NDIM),X(*),B(*)
C
      IF ( ITASK .EQ. 1  .OR. ITASK .EQ. 3 ) THEN
C
C...     CHOLESKY FACTORIZATION
C
         CALL LUCIAQENTER('CHOLF')
         CALL CHLFCB(AL,NDIM,IB,INDEF)
         CALL LUCIAQEXIT ('CHOLF')
      END IF
C
      IF( ITASK .EQ. 2 .OR. ITASK .EQ. 3 ) THEN
C         L * L (T) X = B
C         IS SOLVED IN TWO STEPS
C         1 : L Y = B TO GET Y
C         2 : L(T) X = Y TO GET X
C
         CALL LUCIAQENTER('CHOLS')
         CALL LXEBB(AL,X,B,NDIM,IB)
         CALL COPVEC(X,B,NDIM)
         CALL LTXEBB(AL,X,B,NDIM,IB)
         CALL LUCIAQEXIT('CHOLS')
      END IF
C
      RETURN
      END
      SUBROUTINE CLSKHE(AL,X,B,NDIM,IB,IALOFF,ITASK,INDEF)
C
C Master routine for envelope Cholesky routines .
C Factorize and/or solve set of linear equations for a
C positive definete matrix A.
C The envelope of A is given through IB :
C IB(I) is column number for first nonvanishing element of
C row I
C
C  ITASK = 1 : FACTORIZE MATRIX AND RETURN
C        = 2 : FACTORIZATION HAVE BEEN PERFORMED ( INPUT IN AL )
C              SOLVE LINEAR EQS. MATRIX * X = B
C        = 3 : FACTORIZE AND SOLVE LINEAR EQUATIONS A X = B
C.. INPUT
C
C         AL : ITASK = 1,3 : INPUT MATRIX ( FORMAT : SEE BELOW )
C              OVERWRITTEN !
C              ITASK = 2:  L DECOMPOSITOTATION ASSUMED IN AL )
C              NOT OVERWRITTEN
C
C         X  : VECTOR FOR SOLUTION TO LINEAR EQUATIONS
C         B  : RHS VECTOR FOR LINEAR EQUATIONS( OVERWRITTEN )
C         ( FOR ITASK = 1 X AND B CAN BE DUMMY VARIABLES )
C         NDIM : ORDER OF MATRIX OF MATRICES AND VECTORS
C         IB(I) is column number for first nonvanishing element of
C          row I
C         IALOFF : scratch array .
C         ITASK : DEFINING TASK OF ROUTINE AS ABOVE
C
C OUTPUT :
C        ITASK = 1, 3 : AL IS L DECOMPOSITITION , I.E,
C        L IS A LOWER TRIANGULAR POSITIVE MATRIX AND
C        A = L * L ( TRANSPOSED )
C
C        ITASK = 2,3 : X IS SOLUTION TO LINEAR SET OF EQUATIONS
C        INDEF ( FOR ITASK = 1, 3 ) :
C            0 : MATRIX DECOMPOSED IS NOT INDEFINITE
C        .NE.0 : ABNORMAL TERMINATION DUE TO INDEFINITE MATRIX
C
C NOTE ON STRUCTURE OF MATRIX
C
C THE MATRIX IS ASSUMED PACKED SO ONLY LOWER ELEMENTS of
C THE envelope are stored . The matrix is stored as consecutive rows
C in a one dimensional vector AL
C
C in order to ease indexing an offset vector IALOFF is constructed
C so IALOFF(I) is first adress in AL of first element in row I
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(*),X(*),B(*)
      DIMENSION IB(*),IALOFF(*)
*
      NTEST = 00 
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from CLSKHE:'
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' NDIM = ', NDIM
        WRITE(6,*) ' ITASK = ', ITASK
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Envelope array (IB) '
        CALL IWRTMA(IB,1,NDIM,1,NDIM)
      END IF
C
C
C. Pointer array IALOFF
      IALOFF(1) = 1
      DO 10 I = 1, NDIM - 1
        IALOFF(I+1) = IALOFF(I) + ( I + 1 - IB(I) )
   10 CONTINUE
C
       IF (NTEST .GE. 100) THEN
         WRITE(6,*) ' IALOFF array '
         CALL IWRTMA(IALOFF,1,NDIM,1,NDIM)
       END IF
C
      IF ( ITASK .EQ. 1  .OR. ITASK .EQ. 3 ) THEN
C
C...     CHOLESKY FACTORIZATION
C
         CALL LUCIAQENTER('CHOLF')
         CALL CHLFCE(AL,NDIM,IB,IALOFF,INDEF)
         CALL LUCIAQEXIT ('CHOLF')
*
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Cholesky factorized matrix '
           CALL PRSYM(AL,NDIM)
         END IF
      END IF
C
      IF( ITASK .EQ. 2 .OR. ITASK .EQ. 3 ) THEN
C         L * L (T) X = B
C         IS SOLVED IN TWO STEPS
C         1 : L Y = B TO GET Y
C         2 : L(T) X = Y TO GET X
C
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' The right hand side vector '
           CALL WRTMAT(B,1,NDIM,1,NDIM)
         END IF
         CALL LUCIAQENTER('CHOLS')
         CALL LXEBE(AL,X,B,NDIM,IB,IALOFF)
         CALL COPVEC(X,B,NDIM)
         CALL LTXEBE(AL,X,B,NDIM,IB,IALOFF)
         CALL LUCIAQEXIT('CHOLS')
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' The solution vector '
           CALL WRTMAT(X,1,NDIM,1,NDIM)
         END IF
      END IF
C
      RETURN
      END
      SUBROUTINE CMP2VC(VEC1,VEC2,NDIM,THRES)
C
C COMPARE TWO DOUBLE PRECISION  VECTORS VEC1,AND VEC2
C
C ONLY ELEMENTS THAT DIFFERS BY MORE THAN THRE ARE PRINTED
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
C
      XMXDIF = 0.0D0
      IMXPLC = 0
      WRITE(6,*) ' COMPARISON OF TWO VECTORS '
      WRITE(6,*) '      VECTOR1      VECTOR2        DIFFERENCE '
      DO 100 I = 1, NDIM
        DIF = VEC1(I) - VEC2 ( I )
        IF( ABS(DIF ) .GE. XMXDIF ) THEN
          XMXDIF = ABS(DIF)
          IMXPLC = I
        END IF
        IF( ABS ( DIF ) .GT. THRES ) THEN
          WRITE(6,'(2X,I5,3E15.8)') I,VEC1(I),VEC2(I),DIF
        END IF
  100 CONTINUE
C
      IF( XMXDIF .EQ. 0.0D0 ) THEN
        WRITE(6,*) ' THE TWO VECTORS ARE IDENTICAL '
      ELSE
        WRITE(6,*) ' SIZE AND LAST PLACE OF LARGEST DEVIATION ',
     &  XMXDIF,IMXPLC
      END IF
C
      RETURN
      END
      SUBROUTINE COPDSC(ARRAY,NDIM,NBLOCK,IFROM,ITO)
C
C     COPY DOUBLE PRECISION  ARRAY FROM DISC FILE IFROM TO DISCFILE ITO
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARRAY(1)
C
      IREST=NDIM
      IF ( IREST .EQ. 0 ) GOTO 101
  100 CONTINUE
C     DO 100 WHILE(IREST.GT.0)
C      WHILE (IREST.GT.0)
       IF(IREST.GT.NBLOCK) THEN
        READ(IFROM) (ARRAY(I),I=1,NBLOCK)
        WRITE(ITO) (ARRAY(I),I=1,NBLOCK)
C        IBASE=IBASE+NBLOCK
        IREST=IREST-NBLOCK
       ELSE
        READ(IFROM) (ARRAY(I),I=1,IREST)
        WRITE(ITO) (ARRAY(I),I=1,IREST)
        IREST=0
       END IF
       IF( IREST .GT. 0) GOTO 100
  101 CONTINUE
C      END WHILE
C 100 END  DO
C
      RETURN
      END
      SUBROUTINE COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
*
* Packed version 
C
C
C LBLK DEFINES STRUCTURE OF FILE
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF
 
C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
      SUBROUTINE COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILE
*
* Structure of output file is inherited by output file,
* if input file is packed, so is output file
*
*
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF
 
C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &         NO_ZEROING)
          IF(IAMPACK.NE.0) THEN
C?          WRITE(6,*) ' COPVCD, IAMPACK,FILE = ', IAMPACK,LUIN
          END IF
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
            ELSE  
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            CALL ZERORC(LBL,LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
      SUBROUTINE COPVEC(FROM,TO,NDIM)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/COPVECST/XNCALL_COPVEC, XNMOVE_COPVEC
      INCLUDE 'rou_stat.inc'
C     COMMON/ROU_STAT/NCALL_SCALVE,NCALL_SETVEC,NCALL_COPVEC,
C    &                NCALL_MATCG,NCALL_MATCAS,NCALL_ADD_SKAIIB,
C    &                NCALL_GET_CKAJJB,
C    &                XOP_SCALVE,XOP_SETVEC,XOP_COPVEC,
C    &                XOP_MATCG,XOP_MATCAS,XOP_ADD_SKAIIB,
C    &                XOP_GET_CKAJJB

      DIMENSION FROM(1),TO(1)
C
      XNCALL_COPVEC = XNCALL_COPVEC + 1
      NCALL_COPVEC = NCALL_COPVEC + 1
      XOP_COPVEC = XOP_COPVEC + NDIM
*
      XNMOVE_COPVEC = XNMOVE_COPVEC + NDIM
      DO 100 I=1,NDIM
       TO(I)=FROM(I)
  100 CONTINUE
C
      RETURN
      END
*
      SUBROUTINE DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
C
C VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECOUT(1),VECIN(1),DIAG(1)
C
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIAVC2: '
        WRITE(6,*) ' NDIM = ', NDIM
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'DIAG and VECIN: '
        CALL WRTMAT(DIAG,1,NDIM,1,NDIM)
        CALL WRTMAT(VECIN,1,NDIM,1,NDIM)
      END IF
*
      DO 100 I=1,NDIM
      DIVIDE=DIAG(I)+SHIFT
      THRES=1.0D-10
      IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
      IF(VECIN(I).EQ.0.0D0) THEN
        VECOUT(I) = 0.0D0
      ELSE
        VECOUT(I)=VECIN(I)/DIVIDE
      END IF
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
*
* VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
*
* VDSV = SUM(I) VECIN(I) ** 2 /( DIAG(I) + SHIFT )
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECOUT(1),VECIN(1),DIAG(1)
*
      THRES=1.0D-10
      VDSV = 0.0D0
      DO 100 I=1,NDIM
*
        DIVIDE=DIAG(I)+SHIFT
        IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
*
        VDSV = VDSV + VECIN(I) ** 2 /DIVIDE
        VECOUT(I)=VECIN(I)/DIVIDE
*
  100 CONTINUE
*
      NTEST =00
c     IF(NTEST.GE.100) THEN
      WRITE(6,*) 'DIAVC3 : VECIN, DIAG,VECOUT '
      DO I = 1, NDIM
        WRITE(6,'(3E15.8)') VECIN(I),DIAG(I),VECOUT(I)
      END DO
c     END IF
      RETURN
      END
*
      SUBROUTINE DIAVC3G(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
*
* VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
*
* VDSV = SUM(I) VECIN(I) ** 2 /( DIAG(I) + SHIFT )

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer VECOUT,VECIN,DIAG
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      THRES=1.0D-10
      VDSV = 0.0D0
      call ga_distribution(DIAG, ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(DIAG,ilo,ihi,1,1,idiag,idum)
      call ga_access(VECIN,ilo,ihi,1,1,ivecin,idum)
      if (VECOUT.eq.DIAG) ivecout = idiag  
      if (VECOUT.eq.VECIN) ivecout = ivecin
c
      DO 100 I=0, ihi-ilo  
        DIVIDE=dbl_mb(idiag+i)+SHIFT
        IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
*
        VDSV = VDSV + dbl_mb(ivecin+i) ** 2 /DIVIDE
        dbl_mb(ivecout+i) = dbl_mb(ivecin+i) / DIVIDE
*
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
      WRITE(6,*) 'DIAVC3 : VECIN, DIAG'
      DO I = 0, ihi-ilo
        WRITE(6,'(I5,3E15.8)') i,dbl_mb(ivecout+i),dbl_mb(idiag+i),
     &                          dbl_mb(ivecin+i)
      END DO
      END IF
      call ga_release(DIAG,ilo,ihi,1,1)
      call ga_release(VECIN,ilo,ihi,1,1)
      endif
      call ga_sync()
      call ga_dgop(1,VDSV,1,'+')
      RETURN
      END
*
      SUBROUTINE DMTVCD_OLD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
C
C  IF( INV .NE. 0 ) THEN
C    V3(I) = (V1(I)+FAC)-1 * V2(I)
C    LU3      LU1            LU2
C  IF( INV .EQ. 0 ) THEN
C    V3(I) = (V1(I)+FAC) * V2(I)
C    LU3         LU1        LU2
C WHERE V1 AND V2 ARE VECTORS ON FILES LU1 AND LU2,
C AND LU3 IS WRITTEN ON FILE LU3
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
C
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          CALL IFRMDS(LBL2,1,-1,LU2)
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IBLK = IBLK + 1
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVCD_OLD : '
     &                     , LBL1,LBL2
          WRITE(6,*) 'CURRENT SEGMENT WAS ',IBLK
          STOP ' DIFFERENT BLOCKSIZES IN DMTVCD_OLD '
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,LBL1,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
             CALL DIAVC2(VEC2,VEC2,VEC1,FAC,LBL1)
            ELSE
             CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
             CALL VECSUM(VEC2,VEC1,VEC2,1.0D0,FAC,LBL1)
           END IF
C          CALL TODSC(VEC2,LBL1,KBLK,LU3)
         END IF
         CALL TODSC(VEC2,LBL1,KBLK,LU3)
      END IF
C
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
C mod version where lu1=lu2 is allowed
C
C  IF( INV .NE. 0 ) THEN
C    V3(I) = (V1(I)+FAC)-1 * V2(I)
C    LU3      LU1            LU2
C  IF( INV .EQ. 0 ) THEN
C    V3(I) = (V1(I)+FAC) * V2(I)
C    LU3         LU1        LU2
C WHERE V1 AND V2 ARE VECTORS ON FILES LU1 AND LU2,
C AND LU3 IS WRITTEN ON FILE LU3
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
C
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          IF (LU2.NE.LU1) REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          IF (LU2.NE.LU1) CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          IF (LU1.NE.LU2) THEN
            READ(LU2) LBL2
          ELSE
            LBL2 = LBL1
          END IF
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          IF (LU1.NE.LU2) THEN
            CALL IFRMDS(LBL2,1,-1,LU2)
          ELSE
            LBL2 = LBL1
          END IF
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IBLK = IBLK + 1
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I5)') ' DIFFERENT BLOCKSIZES IN DMTVCD : '
     &                     , LBL1,LBL2
          WRITE(6,'(A,2I3,A,I3,A)') 
     &              ' UNITS: ',LU1, LU2,'(IN) - ',LU3,' (OUT)'
          WRITE(6,*) 'CURRENT SEGMENT WAS ',IBLK
          CALL UNIT_INFO(LU1)
          CALL UNIT_INFO(LU2)
          CALL UNIT_INFO(LU3)
          STOP ' DIFFERENT BLOCKSIZES IN DMTVCD '
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,LBL1,KBLK,LU1,IMZERO,IAMPACK)
          IF (LU2.NE.LU1)
     &      CALL FRMDSC(VEC2,LBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF (LU2.NE.LU1.AND.LBL1.GT.0) THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC2,VEC2,VEC1,FAC,LBL1)
            ELSE
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
              CALL VECSUM(VEC2,VEC1,VEC2,1.0D0,FAC,LBL1)
            END IF
          ELSE IF (LBL1.GT.0) THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC2,VEC1,VEC1,FAC,LBL1)
            ELSE
              CALL VVTOV(VEC1,VEC1,VEC2,LBL1)
              CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,FAC,LBL1)
            END IF
          END IF
          CALL TODSC(VEC2,LBL1,KBLK,LU3)
        END IF
C
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE DMTVCD2(VEC1,VEC2,LU1,LU2,LU3,FAC,DMP,IREW,INV,LBLK)
C mod version where lu1=lu2 is allowed
C
C  IF( INV .NE. 0 ) THEN
C    V3(I) = FAC1 * (V1(I)+DMP)-1 * V2(I)
C    LU3      LU1            LU2
C  IF( INV .EQ. 0 ) THEN
C    V3(I) = FAC1 * (V1(I)+DMP) * V2(I)
C    LU3         LU1        LU2
C WHERE V1 AND V2 ARE VECTORS ON FILES LU1 AND LU2,
C AND LU3 IS WRITTEN ON FILE LU3
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
C
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          IF (LU2.NE.LU1) REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          IF (LU2.NE.LU1) CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          IF (LU1.NE.LU2) THEN
            READ(LU2) LBL2
          ELSE
            LBL2 = LBL1
          END IF
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          IF (LU1.NE.LU2) THEN
            CALL IFRMDS(LBL2,1,-1,LU2)
          ELSE
            LBL2 = LBL1
          END IF
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IBLK = IBLK + 1
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I5)') ' DIFFERENT BLOCKSIZES IN DMTVCD2 : '
     &                     , LBL1,LBL2
          WRITE(6,'(A,2I3,A,I3,A)') 
     &              ' UNITS: ',LU1, LU2,'(IN) - ',LU3,' (OUT)'
          WRITE(6,*) 'CURRENT SEGMENT WAS ',IBLK
          CALL UNIT_INFO(LU1)
          CALL UNIT_INFO(LU2)
          CALL UNIT_INFO(LU3)
          STOP ' DIFFERENT BLOCKSIZES IN DMTVCD2 '
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,LBL1,KBLK,LU1,IMZERO,IAMPACK)
          IF (LU2.NE.LU1)
     &      CALL FRMDSC(VEC2,LBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF (LU2.NE.LU1.AND.LBL1.GT.0) THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC2,VEC2,VEC1,DMP,LBL1)
            ELSE
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
              CALL VECSUM(VEC2,VEC1,VEC2,1.0D0,DMP,LBL1)
            END IF
          ELSE IF (LBL1.GT.0) THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC2,VEC1,VEC1,DMP,LBL1)
            ELSE
              CALL VVTOV(VEC1,VEC1,VEC2,LBL1)
              CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,DMP,LBL1)
            END IF
          END IF
          IF (FAC.NE.1d0) CALL SCALVE(VEC2,FAC,LBL1)
          CALL TODSC(VEC2,LBL1,KBLK,LU3)
        END IF
C
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE EIGENL(A,R,N,MV,MFKR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer*4 ii4,jj4,kk4,ll4
      DIMENSION A(1),R(1)
      DATA TESTIT/1.D-20/
      DATA TESTX/1.D-26/
      DATA TESTY/1.D-18/
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C           MATRIX
C
C        USAGE
C           CALL EIGEN(A,R,N,MV,MFKR)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN ASSCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C       SEQUENCE)
C           MFKR=0 NO SORT
C               =1 SORT
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN ?MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS?, EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
C    1 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
    5 RANGE=1.0D-12
      ii4=mv-1
      IF(ii4) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0D+00
      jj4=i-j
      IF(jj4) 20,15,20
   15 R(IJ)=1.0D+00
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0D+00
      DO 35 I=1,N
      DO 35 J=I,N
      ii4=i-j
      IF(ii4) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(anorm) 165,165,40
   40 ANORM=1.414D+00*DSQRT(ANORM)
      ANRMX=ANORM*RANGE/DBLE(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/DBLE(N)
      IF(THR.LT.TESTY)THR=0.D0
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
      IF(DABS(A(LM)).LT.TESTY)A(LM)=0.D0
      IF(DABS(A(LM)).EQ.0.D0.AND.THR.EQ.0.D0)GO TO 130
   62 IF( DABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5D+00*(A(LL)-A(MM))
      AJUK=(A(LM)*A(LM)+X*X)
      AJUK=DSQRT(AJUK)
      IF(DABS(AJUK).LT.TESTIT)WRITE(6,3000)TESTIT,AJUK,A(LM)
 3000 FORMAT(1H0,'***DENOMINATOR LT ',D12.6,'. VALUE=',D14.8,
     ['. NUMERATOR=',D14.8)
      Y=0.D0
      IF(DABS(AJUK).LT.TESTIT)GO TO 67
      Y=-A(LM)/AJUK
   67 CONTINUE
   68 CONTINUE
C  68 Y=-A(LM)/ DSQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 AJUK=(1.D0-Y*Y)
      IF(AJUK.LT.0.D0)WRITE(6,3001) AJUK
 3001 FORMAT(1H0,'***DSQRT OF ',D14.8)
      IF(AJUK.LT.0.D0)AJUK=0.D0
      AJUK=DSQRT(AJUK)
      AJUK=2.D0*(1.D0+AJUK)
      AJUK=DSQRT(AJUK)
      SINX=Y/AJUK
   76 CONTINUE
C     SINX=Y/ DSQRT(2.0D+00*(1.0D+00+( DSQRT(1.0D+00-Y*Y))))
      SINX2=SINX*SINX
C  78 COSX= DSQRT(1.0D+00-SINX2)
   78 CONTINUE
      AJUK=1.D0-SINX2
      IF(AJUK.LT.TESTX)AJUK=0.D0
      COSX=DSQRT(AJUK)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      ii4=i-l
      IF(ii4) 80,115,80
   80 jj4=i-m
      IF(jj4) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 kk4=I-L
      IF(kk4) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 ll4=MV-1
      IF(ll4) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0D+00*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 ii4=M-N
      IF(ii4) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 jj4=L-(N-1)
      IF(jj4) 145,150,145
  145 L=L+1
      GO TO 55
  150 ll4=IND-1
      IF(ll4) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      IF(MFKR.EQ.0)GO TO 186
  166 CONTINUE
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(MM)-A(LL)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      ii4=mv-1
      IF(ii4) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
186   CONTINUE
      RETURN
      END
 
      REAL*8   FUNCTION FINDMN(VECTOR,NDIM)
C
C FIND SMALLEST ELEMENT OF DOUBLE PRECISION  VECTOR VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
C
      XMIN=VECTOR(1)
      DO 100 I=2,NDIM
       IF(VECTOR(I).LT.XMIN) XMIN=VECTOR(I)
  100 CONTINUE
      FINDMN=XMIN
C
      RETURN
      END
      SUBROUTINE FNDMN2(VEC,NDIM,NVAL,NELMNT,IPLACE,VECORD,NELPVL,
     &                  IPRT)
C
C FIND NVAL LOWEST ELEMENTS IN VEC .
C IF THE SAME VALUE OCCURS SEVERAL TIMES IT IS INCLUDED SEVERAL TIMES
C THE NUMBER OF OCCURENCIES OF THE NVAL LOWEST VALUES ARE RETURNED
C AS NELMNT , AND THEIR VALUES ARE RETURNED IN  VECORD ,AND THEIR
C ORIGINAL PLACE IN IPLACE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(1   ),IPLACE(*),VECORD(*)
      DIMENSION NELPVL(1   )
C
      THRES = 1.0D-8
      CALL ISETVC(IPLACE,0,NVAL)
      CALL ISETVC(NELPVL,0,NVAL)
C
      NELMNT = NVAL
      IELMNT = 0
      IVAL = 0
C. LARGEST ELEMENT TO START WITH
      XMAX = FNDMNX(VEC,NDIM,2)
C
C FIND NEXT LOWEST ELEMENT
1000  CONTINUE
C?    WRITE(6,*) ' START OF LOOP 1000 '
C?    WRITE(6,*) ' IVAL IELMNT ',IVAL,IELMNT
 
C
      XMIN = XMAX
      DO 100 I = 1, NDIM
 
        IF(VEC(I) .LE. XMIN ) THEN
C..  C HECK TO ENSURE THAT I HAS NOT BEEN USED YET
         INEW = 1
         DO 90 JELMNT = 1, IELMNT
           IF(I .EQ. IPLACE(JELMNT)) INEW = 0
   90    CONTINUE
C
         IF( INEW .EQ. 1 ) THEN
           XMIN = VEC(I)
           IMIN = I
         END IF
       END IF
C
C      WRITE(6,*) ' END OF 100 I XMIN IMIN INEW '
C      WRITE(6,*) I,XMIN,IMIN,INEW
C
  100 CONTINUE
C?    WRITE(6,*) ' XMIN AND IMIN ', XMIN,IMIN
C
C
 
      IF(IELMNT .GT. 0 ) THEN
C NEW VALUE ?
COLD    IF(XMIN . EQ. VECORD(IELMNT) ) THEN
        IF( ABS( XMIN-VECORD(IELMNT) ) .LT. THRES ) THEN
          NELMNT = NELMNT + 1
          IPLACE(NELMNT) = 0
          IELMNT = IELMNT + 1
          VECORD(IELMNT) = XMIN
          IPLACE(IELMNT) = IMIN
        ELSE
          IVAL = IVAL + 1
          IF( IVAL .LE. NVAL ) THEN
            IELMNT = IELMNT + 1
            VECORD(IELMNT) = XMIN
            IPLACE(IELMNT) = IMIN
          END IF
        END IF
      ELSE
        IVAL = 1
        IELMNT = 1
        VECORD(1) = XMIN
        IPLACE(1) = IMIN
      END IF
C
      NELPVL(IVAL) = NELPVL(IVAL) + 1
      NELMNT = MIN(NELMNT,NDIM)
      IF( IVAL .LE. NVAL .AND.IELMNT .LT. NDIM) GOTO 1000
C
C
      IF( IPRT  .NE. 0 ) THEN
        WRITE(6,*) ' From FNDMN2 : '
        WRITE(6,*) '   Lowest values '
        CALL WRTMAT(VECORD,1,NELMNT,1,NELMNT)
C       WRITE(6,*) '   places of lowest elements '
C       CALL IWRTMA(IPLACE,1,NELMNT,1,NELMNT)
        WRITE(6,*) '   Number of elements per value '
        CALL IWRTMA(NELPVL,1,NVAL,1,NVAL)
      END IF
C
      RETURN
      END
      REAL*8 FUNCTION FNDMNX(VECTOR,NDIM,MINMAX)
C
C     FIND SMALLEST(MINMAX=1) OR LARGEST(MINMAX=2)
C     ABSOLUTE VALUE OF ELEMENTS IN VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
C
      IF(MINMAX.EQ.1) THEN
       RESULT=ABS(VECTOR(1))
       DO I=2,NDIM
        RESULT=MIN(RESULT,ABS(VECTOR(I)))
       END DO
      END IF
C
      IF(MINMAX.EQ.2) THEN
       RESULT=ABS(VECTOR(1))
       DO I=2,NDIM
        RESULT=MAX(RESULT,ABS(VECTOR(I)))
       END DO
      END IF
C
      IF(MINMAX.EQ.-1) THEN
       RESULT=VECTOR(1)
       DO I=2,NDIM
        RESULT=MIN(RESULT,VECTOR(I))
       END DO
      END IF
C
      IF(MINMAX.EQ.-2) THEN
       RESULT=VECTOR(1)
       DO I=2,NDIM
        RESULT=MAX(RESULT,VECTOR(I))
       END DO
      END IF

      FNDMNX=RESULT
      RETURN
      END
      SUBROUTINE SGATVEC(VECO,VECI,INDEX,NDIM)
C
C GATHER VECTOR with sign encoded:
C VECO(I) = SIGN(INDEX(I))VECI(ABS(INDEX(I))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECI(1),VECO(1   ),INDEX(1   )
*
C
      DO I = 1, NDIM
        IF(INDEX(I).GT.0) THEN
          VECO(I) = VECI(INDEX(I))
        ELSE
          VECO(I) = -VECI(-INDEX(I))
        END IF
      END DO
C
      RETURN
      END
      SUBROUTINE GATVEC(VECO,VECI,INDEX,NDIM)
C
C GATHER VECTOR :
C VECO(I) = VECI(INDEX(I))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECI(1),VECO(1   ),INDEX(1   )
C
      DO 100 I = 1, NDIM
  100 VECO(I) = VECI(INDEX(I))
C
      RETURN
      END
      SUBROUTINE SCAVEC(VECO,VECI,INDEX,NDIM)
C
C SCATTER VECTOR
C VECO(INDEX(I)) = VECI(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECI(1   ),VECO(1),INDEX(1   )
C
      DO 100 I = 1, NDIM
  100 VECO(INDEX(I)) = VECI(I)
C
      RETURN
      END
      SUBROUTINE GPRCTV(DIAG,VECIN,VECUT,NVAR,NPRDIM,IPNTR,
     &                  PEIGVL,PEIGVC,SHIFT,WORK,XH0PSX )
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift ) X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of dimension NPRDIM.
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*
*  NPRDIM : Size of block matrix
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*  PEIGVL   : Eigenvalues of subblock mateix
*  PEIGVC   : Eigenvectors of subblock matrix
*
* Jeppe Olsen , Sept. 1989
*
* Input
*=======
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPRDIM,PEIGVL,PEIGVC : See above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*NPRDIM
*
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* Output
*========
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
 
 
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),PEIGVL(*),PEIGVC(*)
      DIMENSION WORK(*)
*
      IF(NPRDIM.NE.0) THEN
        CALL GATVEC(WORK(1),VECIN,IPNTR,NPRDIM)
* X(T)(DIAG+SHIFT)X in subspace, for later subtraction
        CALL GATVEC(WORK(1+NPRDIM),DIAG,IPNTR,NPRDIM)
        CALL DIAVC3(WORK(1+NPRDIM),WORK(1),
     &       WORK(1+NPRDIM),SHIFT,NPRDIM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPRDIM .NE. 0 ) THEN
         CALL SCAVEC(VECUT,WORK(1),IPNTR,NPRDIM)
         CALL SBINTV(NPRDIM,PEIGVC,PEIGVL,SHIFT,
     &              IPNTR,VECUT,VECUT,WORK(1),WORK(1+NPRDIM),X3)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3
 
 
*
      NTEST = 0
      IF(NTEST.GT. 0 ) THEN
        WRITE(6,*) ' Output vector from GPRCTV '
        WRITE(6,*) ' ========================= '
        CALL WRTMAT(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
      SUBROUTINE H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS
 
*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL LUCIAQENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*
 
      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          CALL TRIPAK(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMAT(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM,ISING)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================
 
         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      CALL LUCIAQEXIT('H0LNS')
      RETURN
      END
      SUBROUTINE H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , May 1990
 
*
*
* =====
* Input
* =====
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* ======
* Output
* ======
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      REAL * 8  INPROD
*
CNW   DIMENSION DIAG(*),VECIN(*),VECUT(*)
      integer diag, vecin, vecut
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      CALL LUCIAQENTER('H0M1T')
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*    
      KLGA = KLFREE
      KLFREE = KLGA + NPQDM
*
      KLSCR = KLFREE
*
      DO I = 0, NPQDM-1
         WORK(KLGA) = 1
      ENDDO
      IF(NPQDM.NE.0) THEN
CNW     CALL GATVEC(dbl_mb(KLV1),VECIN,IPNTR,NPQDM)
        call ga_gather(VECIN,WORK(KLV1),IPNTR,WORK(KLGA),NPQDM)
* X(T)(DIAG+SHIFT)-1 X in subspace, for later subtraction
CNW     CALL GATVEC(dbl_mb(KLV2),DIAG,IPNTR,NPQDM)
        call ga_gather(DIAG,WORK(KLV2),IPNTR,WORK(KLGA),NPQDM)
        CALL DIAVC3(WORK(KLV2),WORK(KLV1),
     &       WORK(KLV2),SHIFT,NPQDM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3G(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPQDM .NE. 0 ) THEN
C                H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
C    &           X,RHS,S,SCR)
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
         X3 = INPROD(WORK(KLV1),WORK(KLV2),NPQDM)
CNW      CALL SCAVEC(VECUT,dbl_mb(KLV2),IPNTR,NPQDM)
         call ga_scatter(VECUT,WORK(KLV2),IPNTR,WORK(KLGA),NPQDM)
         call ga_sync()
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3
 
 
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
CNW     CALL WRTMAT(VECUT,1,NVAR,1,NVAR)
        call ga_print(VECUT)
      END IF
      CALL LUCIAQEXIT('H0M1T')
*
      RETURN
      END
*
      SUBROUTINE H0TV(VECIN,VECUT,DIAG,NVAR,NPQDM,IPNTR,H0,
     &                WORK,NP1,NP2,NQ)
*
* Calculate H0 times vector , where H0 is the diagonal
* approximation plus a P1P2Q preconditioner in a subspace
*
*
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
* Diagonal Times vector
      CALL VVTOV(VECIN,DIAG,VECUT,NVAR)
*
      IF(NPQDM.NE.0) THEN
*.Extract elements belonging to subspace
        CALL GATVEC(WORK(KLV1),VECIN,IPNTR,NPQDM)
        KLPHP = 1
        KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
        KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
C             HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT,WORK)
         CALL HPQTV(NP1,NP2,NQ,H0(KLPHP),H0(KLPHQ),H0(KLQHQ),
     &               WORK(KLV1),WORK(KLV2) )
         CALL SCAVEC(VECUT,WORK(KLV2),IPNTR,NPQDM)
      END IF
*
      NTEST = 0
      IF(NTEST.GT. 0 ) THEN
        WRITE(6,*) ' Output vector from H0TV '
        WRITE(6,*) ' ========================= '
        CALL WRTMAT(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
      INTEGER FUNCTION IBION(M,N)
C
C BIONOMIAL COEFFICIENT (M / N ) = IFAC(M)/(IFAC(M-N)*IFAC(N))
C
*
      INCLUDE 'implicit.inc'
*
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*
* Good old route based on integers
*
      IB = 1
      IF(M-N.GE.N) THEN
         DO K = (M-N+1), M
           IB = IB * K
         END DO    
         IB = IB/IFAC(N)
      ELSE
         DO K = N+1,M
           IB = IB * K
         END DO     
         IB = IB/IFAC(M-N)
      END IF
      IBION = IB
*
      ELSE IF (IWAY.EQ.2) THEN
*
* Use reals
*
        XIB = 1.0D0
        IF(M-N.GE.N) THEN
          DO K = (M-N+1), M
            XK = K
            XIB = XIB * XK
          END DO     
          FACN = IFAC(N)
          XIB = XIB/FACN
        ELSE
          DO K = N+1,M
            XK = K
            XIB = XIB * XK
          END DO    
          FACMN = IFAC(M-N)
          XIB = XIB/FACMN
        END IF
        IBION = NINT(XIB)
      END IF
*
      RETURN
      END
      SUBROUTINE ICOPVE(IFROM,ITO,NDIM)
C
C COPY INTEGER ARRAY
C
      DIMENSION IFROM(1   ),ITO(1   )
C
      DO 100 I = 1,NDIM
        ITO(I) = IFROM(I)
  100 CONTINUE
C
      RETURN
      END
      FUNCTION IFAC(N)
C
C N !
C
      IF( N .LT. 0 ) THEN
       IFAC = 0
       WRITE(6,*) ' WARNING FACULTY OF NEGATIVE NUMBER SET TO ZERO '
      ELSE
C
       IFACN = 1
       DO 100 K = 2,N
        IFACN = IFACN * K
  100  CONTINUE
       IFAC = IFACN
      END IF
C
      RETURN
      END
      SUBROUTINE IFRMDS(IARRAY,NDIM,MBLOCK,IFILE)
C
C     TRANSFER INTEGER ARRAY FROM DISC FILE IFILE
C
C NBLOCK .LT. 0 INDICATES USE OF FASTIO
C
C If nblock .eq. 0 NBLOCK = NDIM
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IARRAY(1)
C
      ICRAY = 1
      NBLOCK = MBLOCK
 
      IF( ICRAY.EQ.1.OR.NBLOCK .GE. 0 ) THEN
C       DO NOT USE FASTIO
        IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
          IF(IREST.GT.NBLOCK) THEN
            READ(IFILE) (IARRAY(IBASE+I),I=1,NBLOCK)
            IBASE=IBASE+NBLOCK
            IREST=IREST-NBLOCK
          ELSE
            READ(IFILE) (IARRAY(IBASE+I),I=1,IREST)
            IREST=0
          END IF
        IF( IREST .GT. 0 ) GOTO 100
      ELSE
C       USE FAST IO
        CALL SQFILE(IFILE,2,IARRAY,NDIM)
      END IF
      RETURN
      END
      SUBROUTINE IFRMDSE(IARRAY,NDIM,MBLOCK,IFILE,IERR)
C
C     TRANSFER INTEGER ARRAY FROM DISC FILE IFILE
C
C     version with error-code
C
C NBLOCK .LT. 0 INDICATES USE OF FASTIO
C
C If nblock .eq. 0 NBLOCK = NDIM
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IARRAY(1)
C
      ICRAY = 1
      NBLOCK = MBLOCK
      IERR = 0 ! begin optimistic
 
      IF( ICRAY.EQ.1.OR.NBLOCK .GE. 0 ) THEN
C       DO NOT USE FASTIO
        IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
          IF(IREST.GT.NBLOCK) THEN
            READ(IFILE,END=201,ERR=202) (IARRAY(IBASE+I),I=1,NBLOCK)
            IBASE=IBASE+NBLOCK
            IREST=IREST-NBLOCK
          ELSE
            READ(IFILE,END=201,ERR=202) (IARRAY(IBASE+I),I=1,IREST)
            IREST=0
          END IF
        IF( IREST .GT. 0 ) GOTO 100
      ELSE
C       USE FAST IO
        CALL SQFILE(IFILE,2,IARRAY,NDIM)
      END IF
      RETURN
 201  IERR = 1 ! end of file reached
      RETURN
 202  IERR = 2 ! I/O-error
      RETURN
      END
      SUBROUTINE IIM1SU(IMAX)
C
C      CREATE ARRAY IIM1AR(I)=I*(I-1)/2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       COMMON/IIM1CM/IIM1AR(5050  )
C
      IIM1AR(1)=0
      IMAXM1=IMAX-1
      DO 100 I=1,IMAXM1
       IIM1AR(I+1)=IIM1AR(I)+I
  100 CONTINUE
C
      RETURN
      END
 
      SUBROUTINE INPACK(A,SCR,NDIM,MATDIM)
C
C PACK LOWER HALF OF MATRIX TO ARRAY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),SCR(NDIM,NDIM)
C
      DO 100 I=1,NDIM
      DO 100 J=1,NDIM
       SCR(J,I)=A(J,I)
  100 CONTINUE
C
      IROW=0
      ICOL=1
C
      DO 200 I=1,NDIM
      DO 200 J=1,I
       IROW=IROW+1
       IF(IROW.GT.MATDIM) THEN
        ICOL=ICOL+1
        IROW=1
       END IF
       A(IROW,ICOL)=SCR(I,J)
  200 CONTINUE
      RETURN
      END
 
      REAL*8 FUNCTION INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
C
C DISC VERSION OF INPROD
C
C LBLK DEFINES STRUCTURE OF FILE
C
*. Last revision, Sept 2003 : FRMDSC => FRMDSC2 to simplify handling 
*                             of vectors containing many zeo blocks
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
      LOGICAL DIFVEC
C
      X = 0.0D0
      IF( LU1 .NE. LU2 ) THEN
        DIFVEC = .TRUE.
      ELSE
        DIFVEC =  .FALSE.
      END IF
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          IF(DIFVEC) REWIND LU2
         ELSE
          CALL REWINE( LU1,LBLK)
          IF( DIFVEC ) CALL REWINE( LU2,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS OF VECTORS
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          IF( DIFVEC) READ(LU2) NBL2
        ELSE IF ( LBLK .LT. 0 ) THEN
          CALL IFRMDS(NBL1,1,-1,LU1)
          IF( DIFVEC)CALL IFRMDS(NBL2,1,-1,LU2)
        END IF
C
        NO_ZEROING = 1
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
C     FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,
C    &                   NO_ZEROING)
          IF( DIFVEC) THEN
            CALL FRMDSC2(VEC2,NBL1,KBLK,LU2,IMZERO2,IAMPACK,
     &                   NO_ZEROING)
            IF(NBL1 .GT. 0 .AND. IMZERO1.EQ.0.AND.IMZERO2.EQ.0)
     &      X = X + INPROD(VEC1,VEC2,NBL1)
          ELSE
          IF(NBL1 .GT. 0 .AND. IMZERO1.EQ.0 )
     &    X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      INPRDD = X
C
      RETURN
      END
      REAL*8 FUNCTION INPRDe(VEC1,VEC2,LU1,LU2,IREW)
C
C DISC VERSION OF INPROD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8   INPROD
      DIMENSION VEC1(*),VEC2(*)
      LOGICAL DIFVEC
C
      X = 0.0D0
      IF( LU1 .NE. LU2 ) THEN
        DIFVEC = .TRUE.
      ELSE
        DIFVEC =  .FALSE.
      END IF
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINO( LU1)
        IF( DIFVEC ) CALL REWINO( LU2)
      END IF
C
C LOOP OVER BLOCKS OF VECTORS
C
 1000 CONTINUE
C
      READ(LU1) NBL1
      IF( DIFVEC) READ(LU2) NBL2
      IF(NBL1 .GE. 0 ) THEN
        CALL FRMDSC(VEC1,NBL1,-1  ,LU1,IMZERO,IAMPACK)
        IF( DIFVEC) THEN
          CALL FRMDSC(VEC2,NBL1,-1  ,LU2,IMZERO,IAMPACK)
          IF(NBL1 .GT. 0 )
     &    X = X + INPROD(VEC1,VEC2,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &    X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IF(NBL1 .GE. 0 ) GOTO 1000
C
      INPRDD = X
      INPRDe = X
C
      RETURN
      END
      REAL*8 FUNCTION INPROD(A,B,NDIM)
C      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(2),B(2)
C
      INPROD=0.0D0
      DO 100 I=1,NDIM
       INPROD=INPROD+A(I)*B(I)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM,ISING)
C FIND INVERSE OF MATRIX A
C INPUT :
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
*
* ISING = 0 => No convergence problems
*       = 1  => Convergence problems
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      ITEST=0
      IF(NDIM.EQ.0) THEN
       RETURN
      ELSE IF(NDIM.EQ.1)THEN
        IF(A(1,1) .NE. 0.0D0 ) THEN
           A(1,1) = 1.0D0/A(1,1)
        ELSE
           ITEST = 1
        END IF
      ELSE
        DETERM=0.0D0
        EPSIL=0.0D0
        CALL BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
      END IF
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
*
      IF(ITEST.NE.0) THEN 
        ISING = 1
      ELSE 
        ISING = 0
      END IF
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE ISETVC(IVEC,IVALUE,NDIM)
C
      DIMENSION IVEC(NDIM)
C
      DO 100 I = 1, NDIM
        IVEC(I) = IVALUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE ISTVC2(IVEC,IBASE,IFACT,NDIM)
C
C IVEC(I) = IBASE + IFACT * I
C
      DIMENSION IVEC(1   )
C
      DO 100 I = 1,NDIM
        IVEC(I) = IBASE + IFACT*I
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE ITODS(IA,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY INTEGER IA(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IA(1)
      INTEGER START,STOP
*
      ICRAY = 1
      NBLOCK = MBLOCK
      IF( NBLOCK .GE.0.OR.ICRAY.EQ.1 ) THEN
C
      IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (IA(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
      END IF
C
      IF(ICRAY.EQ.0.AND. NBLOCK .LT. 0 ) THEN
       CALL SQFILE(IFIL,1,IA,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE IWRTMA10(IMAT,NROW,NCOL,MAXROW,MAXCOL)
* I10 format
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NROW
        WRITE(6,1110) (IMAT(I,J),J= 1,NCOL)
  100 CONTINUE
 1110 FORMAT(/,1X,8I10,/,(1X,8I10))
C
      RETURN
      END
      SUBROUTINE IWRTMA3(IMAT,NROW,NCOL,MAXROW,MAXCOL)
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NROW
        WRITE(6,1110) I,(IMAT(I,J),J= 1,NCOL)
  100 CONTINUE
 1110 FORMAT(/"<",I3,">",1X,20(1X,I3),/,(6X,20(1X,I3)))
C
      RETURN
      END
      SUBROUTINE IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NROW
        WRITE(6,1110) I,(IMAT(I,J),J= 1,NCOL)
  100 CONTINUE
 1110 FORMAT(/"<",I3,">",1X,10I8,/,(6X,10I8))
C
      RETURN
      END
      SUBROUTINE IWRTMA_T(IMAT,NROW,NCOL,MAXROW,MAXCOL)
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NCOL
        WRITE(6,1110) I,(IMAT(J,I),J= 1,NROW)
  100 CONTINUE
 1110 FORMAT(/"<",I3,">",1X,10I8,/,(1X,10I8))
C
      RETURN
      END
      SUBROUTINE LRMTVC(NRANK,NVAR,A,AVEC,VECIN,VECOUT,IZERO)
*
C calculate the product of a low rank matrix and a vector
C the low rank matrix is defined as
C     sum(i,j) avec(j)*a(j,i)*avec(i)t
C        (avec: COLUMN vectors)
C ( IF IZERO .NE. 0 VECOUT IS ZEROED FIRST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8   INPROD
      DIMENSION A(NRANK,NRANK), AVEC(NVAR,NRANK)
      DIMENSION VECIN(1   ),VECOUT(1   )
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from LRMTVC '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Input vector '
        CALL WRTMAT(VECIN,1,NVAR,1,NVAR)
        WRITE(6,*)
        WRITE(6,*) ' Input vectors defining subspace'
        CALL WRTMAT(AVEC,NRANK,NVAR,NRANK,NVAR)
        WRITE(6,*)
        WRITE(6,*) ' Subspace matrix '
        CALL WRTMAT(A,NRANK,NRANK,NRANK,NRANK)
      END IF
*
      IF(IZERO.NE.0) CALL SETVEC(VECOUT(1),0.0D0,NVAR)
      DO 200 I = 1,NRANK
        AVECTV = INPROD(VECIN,AVEC(1,I),NVAR)
        DO 180 J = 1,NRANK
          FACTOR = A(J,I)*AVECTV
          CALL VECSUM ( VECOUT,VECOUT,AVEC(1,J) , 1.0D0,FACTOR,NVAR)
  180   CONTINUE
  200 CONTINUE
*
      IF (NTEST.NE.0) THEN
       WRITE(6,*) ' MATRIX TIMES VECTOR FOR LRMTVC'
       CALL WRTMAT(VECOUT,1,NVAR,1,NVAR)
      END IF
 
      RETURN
      END
      SUBROUTINE LTXEBB(L,X,B,NDIM,IB)
C
C SOLVE L(TRANSPOSED ) X = B
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  L(IB+1,NDIM),X(*),B(*)
C
      CALL COPVEC(B(1),X(1),NDIM)
C
      DO 100 I = NDIM,1,-1
C
        IEFF = MIN(I,IB+1)
        RLII = L(IEFF,I)
        X(I) = X(I) / RLII
        XIM = -X(I)
C
        JMIN = MAX(1,I-IB)
        JMAX = I - 1
        NJ = JMAX - JMIN + 1
C
        CALL VECSUM(X(JMIN),X(JMIN),L(1,I),1.0D0,XIM,NJ)
C
  100 CONTINUE
C
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' X AND B FROM LTBEBB '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
        CALL WRTMAT(B,1,NDIM,1,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE LTXEBE(L,X,B,NDIM,IB,ILOFF)
C
C SOLVE L(TRANSPOSED ) X = B
C
C where L is a lower trinagular matrix stored in envelope fashion
C
C ILOFF(I) Adress in L of first element of row I
C IB(I)    Column number of first row of I
C L : matrix stores rowwise in one dimensional array .
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  L(*),X(*),B(*)
      DIMENSION IB(*),ILOFF(*)
C
C
      DO 100 I = NDIM,1,-1
C
        JMIN = IB(I)
        NJ =   I - JMIN
        IOFF = ILOFF(I)
        RLII = L(IOFF+NJ)
        X(I) = B(I) / RLII
        XIM = -X(I)
C
        CALL VECSUM(B(JMIN),B(JMIN),L(IOFF),1.0D0,XIM,NJ)
C
  100 CONTINUE
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' X AND B FROM LTBEBB '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
        CALL WRTMAT(B,1,NDIM,1,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE LXEBB(L,X,B,NDIM,IB)
C
C SOLVE L X = B
C
C WHERE L IS A LOWER TRIANGULAR MATRIX WITH BAND WIDTH IB,
C AND STORED AS DESCRIBED IN CHLFCB.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  L(IB+1,NDIM),X(*),B(*)
      REAL * 8   INPROD
C
C X AND B CAN BE THE SAME VECTOR
C
      DO 100 I = 1, NDIM
C
        JTERM  = MIN(IB,I-1)
        JSTRT = MAX ( 1, I - IB )
C?      WRITE(6,*) ' I JTERM JSTRT  ',I,JTERM,JSTRT
C
        X(I) =
     &  (B(I)-INPROD(L(1,I),X(JSTRT),JTERM) ) /L(JTERM+1,I)
C
  100 CONTINUE
C
      NTEST = 00
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' X AND B VECTOR '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
        CALL WRTMAT(B,1,NDIM,1,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE LXEBE(L,X,B,NDIM,IB,ILOFF)
C
C SOLVE L X = B
C
C where L is a lower trinagular matrix stored in envelope fashion
C
C ILOFF(I) Adress in L of first element of row I
C IB(I)    Column number of first row of I
C L : matrix stored rowwise in one dimensional array .
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  L(*),X(*),B(*)
      DIMENSION IB(*),ILOFF(*)
      REAL * 8   INPROD
C
C X AND B CAN BE THE SAME VECTOR
C
C x(i) = (b(i)-sum(j) l(i,j)*x(j)) / l(j,j)
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' B VECTOR on input to LXEBE '
        CALL WRTMAT(B,1,NDIM,1,NDIM)
        write(6,*) ' ib and iloff '
        call iwrtma(ib,1,ndim,1,ndim)
        call iwrtma(iloff,1,ndim,1,ndim)
      END IF
C
      DO 100 I = 1, NDIM
        JTERM  = I - IB(I)
        JSTRT =  IB(I)
        IOFF = ILOFF(I)
        X(I) =
     &  (B(I)-INPROD(L(IOFF),X(JSTRT),JTERM) ) /L(IOFF+JTERM)
  100 CONTINUE
C
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' X AND B VECTOR on exit from LXEBE '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
        CALL WRTMAT(B,1,NDIM,1,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE MATDIF(A,B,NMXDIM,MATDIM)
C
C     A=A-B
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMXDIM,2),B(NMXDIM,2)
C
      DO 100 J=1,MATDIM
      DO 100 I=1,MATDIM
       A(I,J)=A(I,J)-B(I,J)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE MATML2(A,B,C,SCR,MATDIM,NDIM,ITRNSP)
C
C             C=A*B.C AND A CAN OCCUPY SAME SPACE
C             LENGTH OF SCR AT LEAST NDIM
C             IF ITRANSP.NE.0 MATRIX A IS TRANSPOSED
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),SCR(2)
      DIMENSION B(MATDIM,MATDIM),C(MATDIM,MATDIM)
C
      IF(ITRNSP.NE.0) CALL TRNSPO(A,MATDIM,NDIM)
C
      DO 300 I=1,NDIM
C
       DO 250 K =1,NDIM
  250  SCR(K)=A(I,K)
C
       DO 200 J=1,NDIM
        X=0.0D0
        DO 100 K=1,NDIM
  100   X=X+SCR(K)*B(K,J)
        A(I,J)=X
  200  CONTINUE
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE MATML3(A,B,C,MATDIM,NDIM,ITRANS)
C
C ANOTHER ROUTINE FOR MATRIX MULT :
C     ITRANS = 0 : C = A*B
C     ITRANS = 1 : C = A(T) * B
C     ITRANS = 2 : C = A * B(T)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8   INPROD
      DIMENSION A(1            ),B(1            ),
     +          C(MATDIM*MATDIM)
 
      CALL SETVEC(C,0.0D0,MATDIM**2)
      IF(ITRANS.EQ.0) THEN
       DO 100 K = 1,NDIM
       DO 100 J = 1,NDIM
C        BKJ = B(K,J)
         BKJ = B( (J-1)*MATDIM + K )
         CALL VECSUM(C((J-1)*MATDIM+1),C((J-1)*MATDIM+1)
     +              ,A((K-1)*MATDIM+1),1.0D0,BKJ,NDIM)
  100  CONTINUE
      END IF
C
      IF(ITRANS.EQ.1) THEN
        DO 200 I = 1,NDIM
        DO 200 J = 1,NDIM
          C((J-1)*MATDIM + I ) =
     &    INPROD(A((I-1)*MATDIM+1),B((J-1)*MATDIM+1),NDIM)
  200   CONTINUE
      END IF
C
      IF(ITRANS.EQ.2) THEN
        DO 300 J = 1,NDIM
        DO 300 K = 1,NDIM
          BJK = B( (K-1)*MATDIM + J)
C          BJK = B(J,K)
          CALL VECSUM(C((J-1)*MATDIM+1),C((J-1)*MATDIM+1)
     +               ,A((K-1)*MATDIM+1),1.0D0,BJK,NDIM)
  300   CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C = A * B             FOR ITRNSP = 0
C
C     C = A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C = A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN, LAST REVISION JULY 24 1987
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
C
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' A AND B MATRIX FROM MATML4 '
        WRITE(6,*)
        CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
        CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
        WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL '
        WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      END IF
*
      IF(ITRNSP.LT.0.OR.ITRNSP.GT.2) THEN
        WRITE(6,*) ' Illegal value of ITRNSP in MATML4 ', ITRNSP
        STOP ' Illegal value of ITRNSP in MATML4 '
      END IF
C
      CALL SETVEC(C,0.0D0,NCROW*NCCOL)
C
      IF( ITRNSP .NE. 0 ) GOTO 001
        DO 50 J = 1,NCCOL
          DO 40 K = 1,NBROW
            BKJ = B(K,J)
            DO 30 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BKJ
  30        CONTINUE
  40      CONTINUE
  50    CONTINUE
C
C
  001 CONTINUE
C
      IF ( ITRNSP .NE. 1 ) GOTO 101
C... C = A(T) * B
         DO 150 J = 1, NCCOL
           DO 140 K = 1, NBROW
             BKJ = B(K,J)
             DO 130 I = 1, NCROW
               C(I,J) = C(I,J) + A(K,I)*BKJ
  130        CONTINUE
  140      CONTINUE
  150    CONTINUE
C
  101 CONTINUE
C
      IF ( ITRNSP .NE. 2 ) GOTO 201
C... C = A*B(T)
        DO 250 J = 1,NCCOL
          DO 240 K = 1,NBCOL
            BJK = B(J,K)
            DO 230 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BJK
 230        CONTINUE
 240      CONTINUE
 250    CONTINUE
C
C
  201 CONTINUE
C
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' C MATRIX FROM MATML4 '
        WRITE(6,*)
        CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
      RETURN
      END
       SUBROUTINE MATMUL(A,B,AB,NMXDIM,MATDIM,ITRANS)
C MULTIPLY MATRICES A AND B AND STORE IN AB
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMXDIM,2),B(NMXDIM,2),AB(NMXDIM,2)
C
      TEST=1.0D-15
      IF(ITRANS.NE.0) THEN
       STOP  1
      END IF
      DO 300 I=1,MATDIM
      DO 300 J=1,MATDIM
       AB(J,I)=0.0D0
  300 CONTINUE
C
      DO 200 K=1,MATDIM
      DO 200 J=1,MATDIM
C
       BKJ=B(K,J)
       IF(ABS(BKJ).GT.TEST) THEN
        DO 100 I=1,MATDIM
         AB(I,J)=AB(I,J)+A(I,K)*BKJ
  100  CONTINUE
       END IF
  200 CONTINUE
      RETURN
      END
      SUBROUTINE MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION   MATRIX(MATDIM,MATDIM),VECIN(2),VECOUT(2)
C
C     VECOUT=MATRIX*VECIN FOR ITRNSP=0
C     VECOUT=MATRIX(TRANSPOSED)*VECIN FOR ITRNSP .NE. 0
C
      DO 10 I=1,NDIM
   10 VECOUT(I)=0.0D0
      IF(ITRNSP.EQ.0) THEN
C
       DO 100 J=1,NDIM
        VECINJ=VECIN(J)
        DO 90 I=1,NDIM
         VECOUT(I)=VECOUT(I)+MATRIX(I,J)*VECINJ
   90   CONTINUE
  100  CONTINUE
      END IF
C
      IF(ITRNSP.NE.0) THEN
       DO 200 I=1,NDIM
        X=0.0D0
        DO 190 J=1,NDIM
         X=X+MATRIX(J,I)*VECIN(J)
  190   CONTINUE
        VECOUT(I)=X
  200  CONTINUE
      END IF
      RETURN
      END
       SUBROUTINE MGS(NDIM,NVECIN,IVCFIL,NVECUT
     +    ,X,A1,A2,B1,B2,MAXVEC)
C
C SUBROUTINE FOR MODIFIED GRAM SCHMIDT ORTHONORMALIZATION.CARE
C HAS BEEN TAKEN IN ORDER TO ASSURE STABLE NUMERICAL PERFORMANCE
C JO 10 MARCH '86
C
C THE NVECIN INPUT VECTORS RESIDE ON  DISCFILE IVCFIL WITH A SPACING OF
C THE NVECUT ORTHOGONALIZED VECTORS IS DESCRIBED BY MATRIX X: X(OLDVEC,N
C THE UNIT BASIS  IS ASSUMED ORTHOGONAL.
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8   INPROD
C
      DIMENSION A1(1  )   ,A2(1     ),B1(1   ),B2(1   )
      DIMENSION X(MAXVEC,MAXVEC)
C
      NTEST=1
      IF(NVECIN.GT.MAXVEC) THEN
       WRITE(6,1011)NVECIN,MAXVEC
 1011  FORMAT(1H0,' ACTUAL SUBSPACE DIMENSION',I3,
     +        ' GREATER THAN ALLOWED MAXIMUM',I3,'!!!!!')
       STOP
      END IF
C
C
       XMAX = 1.0D+06
       IEFF = 0
       DO 10 I = 1,NVECIN
        DO 9 J = 1,NVECIN
    9  X(I,J)=0.0D0
   10  X(I,I)=1.0D0
C
C LOOP OVER NEW VECTORS
       DO 600 I=1,NVECIN
C**     UNNORMALIZED VECTOR I
        ISTOP = I
        IADD  = 1
        IMULT = 2
        CALL SETVEC (B1,0.0D0,NDIM)
        CALL VECSMF(B1,X(1,I),B2,ISTOP,IMULT,IADD,IVCFIL,NDIM)
C**     NORMALIZE
        BNORM = INPROD(B1,B1,NDIM)
        SCALE = 1.0D0/SQRT(BNORM)
        CALL SCALVE(B1,SCALE,NDIM)
        CALL SCALVE(X(1,I),SCALE,NVECIN)
C        WRITE(6,*) ' I X(*,I) ',I
C        CALL WRITVE( X(1,I),NVECIN)
        XLARGE = FNDMNX(X(1,I),NVECIN,2)
        IF( ABS(XLARGE) .LE. XMAX ) THEN
C**      NEW VECTOR IS OKAY SO
         IEFF = IEFF + 1
         IF( IEFF .NE. I ) CALL COPVEC(X(1,I),X(1,IEFF),NVECIN)
C**      ORTHOGONALIZE REMAINING VECTORS TO THIS VECTOR
         IPL1 = I + 1
C*       OVERLAP BETWEEN NEW VECTOR AND ORIGINAL VECTORS
         CALL REWINO( IVCFIL)
         DO 500 J=1,NVECIN
          IF( J.NE.1) READ(IVCFIL)
          CALL FRMDSC(B2,NDIM,-1  ,IVCFIL,IMZERO,IAMPACK)
          A1(J)=INPROD(B1,B2,NDIM)
          CALL MATVCB(X,A1,A2,MAXVEC,NVECIN,1)
  500    CONTINUE
C*       ORTHOGONALIZE
         DO 450 K=IPL1,NVECIN
          FAC1 = 1.0D0/(1.0D0 - A2(K)**2)
          FAC2=  -A2(K)*FAC1
          CALL VECSUM(X(1,K),X(1,K),X(1,IEFF),FAC1,FAC2,NVECIN)
  450   CONTINUE
       END IF
  600 CONTINUE
C
C
      NVECUT = IEFF
      IF( NVECUT.NE.NVECIN)
     + WRITE(6,1010 ) NVECUT
 1010 FORMAT(1H0,' number of vectors reduced to..',I4)
      RETURN
      END
C the Structure of files on the following can have one of
C three structures.The type of structure is defined by
C a parameter LBLK
C
C LBLK .GT. 0 :
C==============
C Each record is a single block of length LBLK,
C file has structure
C      Record 1
C      Record 2
C        etc
C so no information about block size and end of record is
C given
C
C LBLK .EQ. 0 .
C==============
C Each record can consist of several blocks, information about
C length of block and end of record explicitly written on file
C file has structure
C Loop over records
C  Loop over blocks of record
C    LLBLK : .GE. 0 : length of next block
C            .LT. 0 : End of record
C    block of size LLBLK
C  End of loop over blocks
C End of loop over records
C
C LBLK .LT. 0
C=============
C As LBLK .EQ. 0 , but use FASTIO routines to write/read files
      SUBROUTINE MICDV4O(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
*
* Jeppe Olsen Winter of 1991
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
       REAL*8 INPRDD
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
 
       IOLSTM = 0
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
*
       CALL MEMMAN(KAPROJ,MAXVEC*(MAXVEC+1)/2,'ADDL  ',2,'KAPROJ')
CNW    KAPROJ = 1
CNW    KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINE(LU1,LBLK)
       CALL REWINE(LU2,LBLK)
       DO 10 IVEC = 1,NINVEC
         CALL REWINO(LU3)
         CALL REWINO(LU4)
C             COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
         CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU3,LU4,0,0)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINE(LU4,LBLK)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINE(LU2,LBLK)
         DO 8 JVEC = 1, IVEC
           CALL REWINE(LU3,LBLK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU3,0,LBLK)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,dbl_mb(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(dbl_mb(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        IF(IPRT  .GE. 10 ) THEN
         WRITE(6,*) ' Info from iteration .... ', ITER
        END IF
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         EIGAPR = EIG(ITER-1,IROOT)
*
         CALL REWINE(LU1,LBLK)
         CALL REWINE(LU2,LBLK)
         EIGAPR = EIG(ITER-1,IROOT)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINE( LU3, LBLK )
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(6,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
C                 WRTVCD(SEGMNT,LU,IREW,LBLK)
         END IF
*  Strange place to put convergence but ....
C                      INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
C          CALL REWINO( LUDIA)
C          CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
C          CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
C    &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ)
*. Inverted diagonal times (HX-EX) on LU3
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU4,LU3,-EIGAPR,1,1,LBLK)
 
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* Current eigen-vector on LU4
             CALL REWINE(LU1,LBLK)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 CALL REWINE( LU4, LBLK )
                 CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
               ELSE
                 CALL REWINE(LU5,LBLK)
                 CALL REWINE(LU4,LBLK)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU4,LU5,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(6,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
C             CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
C    &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ)
             CALL DMTVCD(VEC1,VEC2,LUDIA,LU5,LU4,-EIGAPR,1,1,LBLK)
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPRDD(VEC1,VEC2,LU5,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
C                          VCSMDN(VEC1,VEC2,FAC1,FAC2,LU1,LU2,IREW,LBLK)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU4,LU5,1,LBLK))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                DELTA = INPRDD(LU4,LU3,VEC1,VEC2,1,LBLK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
C                    VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU5,LU4,1,LBLK)
                CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
           CALL REWINE( LU1 ,LBLK)
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             WORK(IVEC) = INPRDD(VEC1,VEC2,LU1,LU3,0,LBLK)
   80      CONTINUE
*
           CALL REWINE(LU1,LBLK)
           DO 82 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
             CALL VECSMD(VEC1,VEC2,-WORK(IVEC),1.0D0,LU1,LU3,
     &                   LU4,0,LBLK)
             CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
   82      CONTINUE
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
        CALL MV7(VEC1,VEC2,LU3,LU4,0,0)
        CALL REWINE(LU4,LBLK)
        CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO 140 JVEC = 1, NVEC+IVEC
          CALL REWINE(LU4,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU4,0,LBLK)
  140   CONTINUE
  150 CONTINUE
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,dbl_mb(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGENL(dbl_mb(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = dbl_mb(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
*. overwrite with pert solution
        C1NRM = SQRT(C0VAR**2 + C1VAR**2)
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        dbl_mb(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
      IPULAY = 1
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3*NROOT .AND.NVEC.GE.2*NROOT) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
*. Current trial vectors
        CALL REWINE( LU5,LBLK)
        DO 421 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  421   CONTINUE
*. Previous trial vectors orthogonalized
C                ORTVCD(LUIN,LUVEC,LUOUT,LUSCR,VEC1,VEC2,NVEC,LBLK,
C    &                  SCR,INORMA)
        CALL REWINE(LU1,LBLK)
        DO 430 IROOT = 1, NROOT
          CALL ORTVCD(LU1,LU5,LU3,LU4,VEC1,VEC2,NROOT+IROOT-1,LBLK,
     &                AVEC((IROOT-1)*NVEC+1),1)
  430   CONTINUE
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
*. Transfer C vectors to LU1
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 441 IVEC = 1,2*NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  441   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 450 IROOT = 1, 2*NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  450   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 460 IVEC = 1,2*NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  460   CONTINUE
        NVEC = 2*NROOT
*
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
       CALL REWINO( LU1 )
       DO 2010 IVEC = 1,NVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL REWINO( LU2)
         DO 2008 JVEC = 1, IVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
 
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINO(LU1)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE MINDAV(VEC1,VEC2,LU1,LU2,RNRM,EIG,EIGAPR,MAXIT,NVAR,
     &                  LU3,LUDIA)
C
C MINIMAL DAVIDSON ALGORITHM WITH ONLY TWO VECTOR SEGMEMNTS IN CORE .
C
C INPUT :
C
C        VEC1,VEC2 : TWO VECTORS,EACH MUST HOLD LATGEST BLOCK OF
C        VECTOR
C
C       LU1,LU2,LU3 : TWO SCRATCH FILES
C       LUDIA : FILE CONTAINING CI DIAGONAL
C ON INPUT VEC1/LU1 IS SUPPOSED TO HOLD INITIAL GUESS TO EIGENVECTOR
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(1    ),EIG(1    )
       LOGICAL CONVER
C
       NTEST = 1
       TEST = 1.0D-5
       CONVER = .FALSE.
       DO 1234 MACRO = 1,1
C
C... INITAL ITERATION
       ITER = 1
*. Does not work ... 
       CALL MV7(VEC1,VEC2,LU1,LU2,0,0)
       EIGAPR = INPROD(VEC1,VEC2,NVAR)
       EIG(ITER) = EIGAPR
C
       CALL REWINO( LU1       )
       CALL REWINO( LU2)
       CALL TODSC(VEC1,NVAR,-1  ,LU1)
       CALL TODSC(VEC2,NVAR,-1  ,LU2)
C
C
C** LOOP OVER ITERATIONS
C
 1000 CONTINUE
      IF(NTEST .GE. 10 ) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
       WRITE(6,*) ' EIGEN VECTOR APPROXIMATION '
       CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
       WRITE(6,*)
       WRITE(6,*) 'MATRIX TIMES EIGEN VECTOR APPROXIMATION '
       CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
      END IF
 
        ITER = ITER + 1
C
C *** 1        : NEW DIRECTION TO BE INCLUDED
C
C.  1.1 : R = H*X - EIGAPR*X IN VEC2
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-EIGAPR,NVAR)
      IF(NTEST .GE. 10 ) THEN
       WRITE(6,*) ' HX - EX '
       CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
      END IF
C
C... STRANGE PLACE FOR CONVERGENCE TEST , BUT.
        RNORM = INPROD(VEC2,VEC2,NVAR)
        RNORM = SQRT(RNORM)
        RNRM(ITER-1) = RNORM
        IF(RNORM.LT. TEST ) THEN
          CONVER = .TRUE.
          GOTO 1001
        END IF
C.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
       CALL REWINO( LUDIA)
       CALL FRMDSC(VEC1,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
       CALL DIAVC2(VEC2,VEC2,VEC1,-EIGAPR,NVAR)
       IF ( NTEST .GE. 10 ) THEN
         WRITE(6,*) ' (D-E)-1 *( HX - EX ) '
         CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
       END IF
 
C.  1.3 : ORTHOGONALIZE R TO CURRENT EIGEN VECTOR APROXIMATION (VEC2)
       CALL REWINO( LU1   )
        CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        OVRLP = INPROD(VEC2,VEC1,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVRLP,NVAR)
       IF ( NTEST .GE. 10 ) THEN
         WRITE(6,*) 'ORTHOGONALIZED (D-E)-1 *( HX - EX ) '
         CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
       END IF
C
C
C.  1.4 : NORMALIZE NEW DIRECTION TO 1
       CALL SCALE2(VEC2,NVAR,FACTOR)
C.   NEW DIRECTION IS NOW IN VEC2, SAVE IN LU3
       CALL REWINO( LU3    )
       CALL TODSC(VEC2,NVAR,-1  ,LU3)
 
C
C.. 2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
C
C. 2.1: MULTIPLY NEW DIRECTION WITH MATRIX
       CALL MV7(VEC2,VEC1,LU1,LU2,0,0)
C. 2.2: 2 BY 2 PROJECTED MATRIX
       E00 = EIGAPR
       E11 = INPROD(VEC2,VEC1,NVAR)
C PREVIOUS X VECTOR IN VEC2
       CALL REWINO( LU1     )
       CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
       E01 = INPROD(VEC1,VEC2,NVAR)
C?     WRITE(6,*) ' E00,E01,E11',E00,E01,E11
C
C
C                           E00  E01
C   LOWEST EIGENVALUE OF               IS
C                           E01  E11
C
C
C  -B/2  - SQRT(B**2-4C)/2  WHERE
C
C  B = -(E00+E11), C = E00*E11-E01*E01
       B = -E00-E11
       C = E00*E11-E01**2
C      WRITE(6,*) 'B C ',B,C
C
       EIGAPR = -B/2.0D0 - SQRT(B*B - 4*C )/2.0D0
C
C NEW EIGENVECTOR IN TWO VECTOR BASE
      FAC = SQRT(1.0D0+(E00-EIGAPR)**2/E01**2 )
      FAC = DSQRT(1.0D0+((E00-EIGAPR)/E01)**2)
      X1 = 1.0D0/FAC
      X2 = -(E00-EIGAPR)/E01/FAC
      IF(ABS(E01) .LE. 1.0D-5 ) THEN
C FIRST ORDER CORRECTION
        DELTA = E01/(E00-E11)
        FAC = 1.0D0/DSQRT(1.0D0+DELTA**2)
        X1 = FAC
        X2 = DELTA * FAC
        EIGAPR =(E00*X1**2 + E11*X2**2 + 2*E01*X1*X2)/(X1**2+X2**2)
      END IF
       WRITE(6,*) ' EIGAPR',EIGAPR
C?     WRITE(6,*) ' E00 - EIGAPR ',E00-EIGAPR
       EIG(ITER) = EIGAPR
C?    WRITE(6,*) ' X1, X2', X1,X2
C?    EIGAP2 =(E00*X1**2 + E11*X2**2 + 2*E01*X1*X2)/(X1**2+X2**2)
C?    WRITE(6,*) ' ANOTHER ENERGY EVALUTION GIVES ',EIGAP2
C
C OFF DIAGONAL ELEMENT IN NEW BASIS
C ?   DELTA = (E00-E11)*X1*X2 + E01*(X2**2-X1**2)
C ?   WRITE(6,*) ' NEW OFF DIAGONAL MATRIX ELELMENT ',DELTA
C
C** 3 : PREPARE FOR NEXT ITERATION
C
C H TIMES CURRENT CI VECTOR
      CALL REWINO( LU2       )
      CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,X2,X1,NVAR)
      CALL REWINO( LU2        )
      CALL TODSC(VEC1,NVAR,-1  ,LU2)
C CURRENT CI VECTOR TO DISC
      CALL REWINO( LU1         )
      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      CALL REWINO( LU3          )
      CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,X1,X2,NVAR)
      CALL REWINO( LU1           )
      CALL TODSC(VEC1,NVAR,-1  ,LU1)
      CALL REWINO( LU2 )
      CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
C
C
      IF( ITER .LT. MAXIT ) GOTO 1000
 1001 CONTINUE
C
C
C
      IF( .NOT. CONVER ) THEN
C       CONVERGENCE WAS NOT OBTAINED
        WRITE(6,1170) MAXIT
 1170   FORMAT('0  CONVERGENCE WAS NOT OBTAINED IN ',I3,'ITERATIONS')
      ELSE
C       CONVERGENCE WAS OBTAINED
        ITER = ITER - 1
        WRITE(6,1180) ITER
 1180   FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
       END IF
C
      WRITE(6,1190) EIGAPR
 1190 FORMAT(' The final approximation to eigenvalue ',F18.10)
C     WRITE(6,1200)
C1200 FORMAT(1H0,'THE FINAL APPROXIMATION TO EIGENVECTOR')
C     WRITE(6,1030) (VEC1(I),I=1,NVAR)
      WRITE(6,1300)
 1300 FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
      WRITE(6,1310)
 1310 FORMAT(1H0,' Iteration point      Eigenvalue         Residual ')
      DO 1330 I=1,ITER
 1330 WRITE(6,1340) I,EIG(I),RNRM(I)
 1340 FORMAT(1H ,6X,I4,8X,F18.13,2X,E12.5)
 1234 CONTINUE
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE MINDV4(MV7,
     &                  VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO,CONVER,RNRM_CNV,IROOT_SEL)
*
* Davidson algorithm , requires two vectors in core
* Multi root version
*
* Allows updating of preconditioning matrix so this is
* the current eigenvector approximation
* is an eigenvector for the preconditioner
*
* Jeppe Olsen Sept 89
*             Jan  92 : MV7 entry
*             Feb. 13: IROOT_SEL added 
*
* Input :
* =======
*        MV7 : Name of routine performing matrix*vector calculation
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
* On input LU1 is supposed to hold initial guess to eigenvectors
*
* IOLSEN : Use inverse iteration modified Davidson
* IPICO  : Use perturbation estimate of new vector instead of
*          variational method
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
       integer VEC1,VEC2
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
       DIMENSION RNRM_CNV(*)
       character*100 myname,filename
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(1000)
*
       EXTERNAL MV7
*
       TEST = 1.0D-6
       IPRT= 1
       IF(IPRT.GE.1.and.ga_nodeid().eq.0) THEN
         WRITE(6,*) ' MINDV4 in action '
         WRITE(6,*) ' Convergence threshold for residual = ', TEST
       END IF
*
       IOLSTM = IOLSEN
       if (ga_nodeid().eq.0) then
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
       endif
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MINDV4 wounded , MAXVEC .LT. 2*NROOT '
         STOP ' Enforced stop in MINDV4'
       END IF
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*
       CALL MEMMAN(KAPROJ,MAXVEC*(MAXVEC+1)/2,'ADDL  ',2,'KAPROJ')
CNW    KAPROJ = 1
CNW    KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       CONVER = .FALSE.
       DO 1234 MACRO = 1,1
*
       CALL LUCIAQENTER('MINDV')
*.   INITAL ITERATION
       ITER = 1
CNW    CALL REWINO( LU1 )
CNW    CALL REWINO( LU2 )
       DO 10 IVEC = 1,NINVEC
C?       WRITE(6,*) ' Before FRMDSC, NVAR = ', NVAR
CNW      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &       IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
         if(dra_wait(itask).ne.0) call errquit('dra read err',itask,911)
*
C?       CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
*
         call ga_zero(VEC2) 
         CALL MV7(VEC1,VEC2,0,0,0,0)
CNW      CALL TODSC(VEC2,NVAR,-1  ,LU2)
!         write(*,*) ga_nodeid(), NVAR, IVEC
         if (dra_write_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &       IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
         if(dra_wait(itask).ne.0) call errquit('dra writ err',itask,911)
         
*        PROJECTED MATRIX
CNW      CALL REWINO( LU1)
         DO 8 JVEC = 1, IVEC
CNW        CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &        JVEC,JVEC,itask).ne.0) call errquit('dra error',itask,911)
           if (dra_wait(itask).ne.0) call errquit('read',itask,911)
           IJ = IVEC*(IVEC-1)/2 + JVEC
CNW        APROJ(IJ) = INPROD(VEC1,VEC2)
           APROJ(IJ) = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                               VEC2,'N',1,NVAR,1,1)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.10 .and. ga_nodeid().eq.0) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*  DIAGONALIZE INITIAL PROJECTED MATRIX
       CALL COPVEC(APROJ,dbl_mb(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(dbl_mb(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF( IPRT  .GE. 3 .and. ga_nodeid().eq.0) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F22.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NROOT)
       END IF
*. No root selection here
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
** LOOP OVER ITERATIONS
*
 1000 CONTINUE
      IF(IPRT  .GE. 5 .and. ga_nodeid().eq.0) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
      END IF
 
 
        ITER = ITER + 1
*
** 1          NEW DIRECTION TO BE INCLUDED
*
*   1.1 : R = H*X - EIGAPR*X
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
CNW      CALL SETVEC(VEC1,0.0D0,NVAR)
         call ga_zero(VEC1)
*   
CNW      CALL REWINO( LU2)
         DO 60 IVEC = 1, NVEC
CNW        CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPAC)
           if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
           if (dra_wait(itask).ne.0) call errquit('dra read err',
     &         itask,911)
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           call ga_add(FACTOR,VEC2,1.0D0,VEC1,VEC1)
CNW        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   60    CONTINUE
         EIGAPR = EIG(ITER-1,IROOT)
CNW      CALL REWINO( LU1)
         DO 50 IVEC = 1, NVEC
CNW        CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU1,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
           if (dra_wait(itask).ne.0) call errquit('dra read err',
     &         itask,911)
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           call ga_add(FACTOR,VEC2,1.0D0,VEC1,VEC1)
CNW        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             call ga_print(VEC1)
CNW          CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
CNW      RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNORM = SQRT( ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                               VEC1,'N',1,NVAR,1,1))
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
CNW        CALL REWINO( LUDIA)
CNW        CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           if (dra_read_section(.false.,VEC2,1,NVAR,1,1,
     &         LUDIA,1,NVAR,1,1,itask).ne.0)
     &         call errquit('dra error',itask,911)
           if (dra_wait(itask).ne.0) call errquit('dra read err',
     &         itask,911)
           CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             call ga_print(VEC1)
CNW          CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
CNW           CALL REWINO(LU3)
CNW           CALL TODSC(VEC1,NVAR,-1,LU3)
              if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &            1,1,itask).ne.0) call errquit('dra error',itask,911)
              if (dra_wait(itask).ne.0) call errquit('dra read err',
     &            itask,911)
* Current eigen vector
CNW           CALL REWINO( LU1)
CNW           CALL SETVEC(VEC1,0.0D0,NVAR)
              call ga_zero(VEC1)
              DO 59 IVEC = 1, NVEC
CNW             CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
               if (dra_read_section(.false.,VEC2,1,NVAR,1,1,
     &            LU1,1,NVAR,IVEC,IVEC,itask).ne.0)
     &            call errquit('dra error',itask,911)
               if (dra_wait(itask).ne.0) call errquit('dra read err',
     &             itask,911)
                FACTOR = AVEC((IROOT-1)*NVEC+ IVEC)
CNW             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
                call ga_add(FACTOR,VEC2,1.0D0,VEC1,VEC1)
   59         CONTINUE
              IF ( IPRT  .GE. 600 ) THEN
                WRITE(6,*) ' And X  '
CNW             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
              END IF
CNW           CALL TODSC(VEC1,NVAR,-1,LU3)
              if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &            2,2,itask).ne.0) call errquit('dra error',itask,911)
              if (dra_wait(itask).ne.0) call errquit('dra read err',
     &            itask,911)
* (H0 - E )-1  * X
CNW           CALL REWINO( LUDIA)
CNW           CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
               if (dra_read_section(.false.,VEC2,1,NVAR,1,1,
     &            LUDIA,1,NVAR,1,1,itask).ne.0)
     &            call errquit('dra error',itask,911)
              if (dra_wait(itask).ne.0) call errquit('dra read err',
     &            itask,911)
              CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
     &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                   IPRT)
CNW           CALL TODSC(VEC2,NVAR,-1,LU3)
              if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &            3,3,itask).ne.0) call errquit('dra error',itask,911)
              if (dra_wait(itask).ne.0) call errquit('dra read err',
     &            itask,911)
* Gamma = X(T) * (H0 - E) ** -1 * X
CNW           GAMMA = INPROD(VEC2,VEC1,NVAR)
              GAMMA = ga_ddot_patch(VEC2,'N',1,NVAR,1,1,
     &                              VEC1,'N',1,NVAR,1,1)

CBERT: H0M1TV and subsequent DDOT can be combined as we don't need
CVEC3(2) anymore after we're done

* is X an eigen vector for (H0 - 1 ) - 1
              call ga_add(GAMMA,VEC1,-1.0D0,VEC2,VEC2)
CNW           CALL VECSUM(VEC2,VEC1,VEC2,GAMMA,-1.0D0,NVAR)
              VNORM = SQRT(MAX(0.0D0,ga_ddot_patch(VEC2,'N',1,NVAR,1,1,
     &                                            VEC2,'N',1,NVAR,1,1)))
CNW           VNORM = SQRT(MAX(0.0D0,INPROD(VEC2,VEC2,NVAR)))

              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.10) WRITE(6,*) ' Olsen Correction active '
CNW             CALL REWINO(LU3)
CNW             CALL FRMDSC(VEC2,NVAR,-1,LU3,IMZERO,IAMPACK)
                if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU3,1,NVAR,
     &              1,1,itask).ne.0) call errquit('dra error',itask,911)
                if (dra_wait(itask).ne.0) call errquit('dra read err',
     &              itask,911)
CNW             DELTA = INPROD(VEC1,VEC2,NVAR)
                DELTA = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                                VEC2,'N',1,NVAR,1,1)
CNW             CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
CNW             CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &              3,3,itask).ne.0) call errquit('dra error',itask,911)
                if (dra_wait(itask).ne.0) call errquit('dra read err',
     &              itask,911)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.10.and.ga_nodeid().eq.0) 
     &            WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.10.and.ga_nodeid().eq.0) 
     &            WRITE(6,*)   DELTA,GAMMA,FACTOR
CNW             CALL VECSUM(VEC1,VEC1,VEC2,FACTOR,1.0D0,NVAR)
                call ga_add(FACTOR,VEC2,1.0D0,VEC1,VEC1)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) '  Modified new trial vector '
CNW               CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
                  call ga_print(VEC1)
                END IF
              ELSE
                IF(IPRT.GT.0.and.ga_nodeid().eq.0) WRITE(6,*) 
     &          ' Inverse correction switched of'
CNW             CALL REWINO(LU3)
CNW             CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &              1,1,itask).ne.0) call errquit('dra error',itask,911)
                if (dra_wait(itask).ne.0) call errquit('dra read err',
     &              itask,911)
              END IF
            END IF
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
           XNRMI = ga_ddot_patch(VEC1,'N',1,NVAR,1,1, 
     &                           VEC1,'N',1,NVAR,1,1)
CNW        XNRMI =    INPROD(VEC1,VEC1,NVAR)
CNW        CALL REWINO( LU1 )
 
           DO 80 IVEC = 1,NVEC+IADD-1
CNW          CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
             if (dra_read_section(.false.,VEC2,1,NVAR,1,1,
     &          LU1,1,NVAR,IVEC,IVEC,itask).ne.0)
     &          call errquit('dra error',itask,911)
             if (dra_wait(itask).ne.0) call errquit('dra read err',
     &           itask,911)
CNW          OVLAP = INPROD(VEC1,VEC2,NVAR)
             OVLAP = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                             VEC2,'N',1,NVAR,1,1)
CNW          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
             call ga_add(-OVLAP,VEC2,1.0D0,VEC1,VEC1)
   80      CONTINUE
*. 1.4 Normalize vector and check for linear dependency
CNW        SCALE = INPROD(VEC1,VEC1,NVAR)
           SCALE = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                           VEC1,'N',1,NVAR,1,1)
           IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
             IADD = IADD - 1
             IF ( IPRT  .GE. 10 .and. ga_nodeid().eq.0) THEN
               WRITE(6,*) '  Trial vector linear dependent so OUT !!! '
             END IF
           ELSE
             C1NRM = SQRT(SCALE)
             FACTOR = 1.0D0/SQRT(SCALE)
CNW          CALL SCALVE(VEC1,FACTOR,NVAR)
             call ga_scale(VEC1,FACTOR)
*
CNW          CALL TODSC(VEC1,NVAR,-1  ,LU1)
             if (dra_write_section(.false.,VEC1,1,NVAR,1,1,
     &          LU1,1,NVAR,NVEC+IADD,NVEC+IADD,itask).ne.0) 
     &          call errquit('dra error',itask,911)
             if (dra_wait(itask).ne.0) call errquit('dra read err',
     &           itask,911)
             IF ( IPRT  .GE.600 ) THEN
               WRITE(6,*) 'ORTHONORMALIZED (D-E)-1 *( HX - EX ) '
CNW            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
               call ga_print(VEC1)
             END IF
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
*  2.1: MULTIPLY NEW DIRECTION WITH MATRIX
CNW    CALL REWINO( LU1)
CNW    CALL REWINO( LU2)
CNW    DO 110 IVEC = 1, NVEC
CNW      CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
CNW      CALL FRMDSC(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
CNW  110  CONTINUE
*
      DO 150 IVEC = 1, IADD
CNW     CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        if (dra_read_section(.false.,VEC1,1,NVAR,1,1,
     &     LU1,1,NVAR,NVEC+IVEC,NVEC+IVEC,itask).ne.0) 
     &     call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
        CALL MV7(VEC1,VEC2,0,0,0,0)
CNW     CALL TODSC(VEC2,NVAR,-1  ,LU2)
        if (dra_write_section(.false.,VEC2,1,NVAR,1,1,
     &     LU2,1,NVAR,NVEC+IVEC,NVEC+IVEC,itask).ne.0)
     &     call errquit('dra error',itask,911)
        if(dra_wait(itask).ne.0) call errquit('dra read err',itask,911)
*   AUGMENT PROJECTED MATRIX
CNW     CALL REWINO( LU1)
        DO 140 JVEC = 1, NVEC+IVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
CNW       CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &        JVEC,JVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
          APROJ(IJ) = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                              VEC2,'N',1,NVAR,1,1)
  140   CONTINUE
  150 CONTINUE
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,dbl_mb(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGENL(dbl_mb(KAPROJ),AVEC,NVEC,0,1)
*. Select if required the roots to be followed
      IF(IROOT_SEL.NE.0) THEN
        ISEL_MET = IROOT_SEL 
        if (ga_nodeid().eq.0) WRITE(6,*) ' I will do root selection '
*
        IF(IPRT .GE. 30 .and. ga_nodeid().eq.0) THEN
          WRITE(6,*) ' Info before selection: '
          WRITE(6,*) ' Projected matrix and eigen vectors '
          CALL PRSYM(APROJ,NVEC)
          CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
        END IF
C       SEL_ROOT(SUBSPCVC,SUBSPCMT,ISEL_MET,NVEC,NROOT,LUC,VEC1)
        CALL SEL_ROOT(AVEC,dbl_mb(KAPROJ),ISEL_MET,NVEC,NROOT,LU1,VEC1)
      END IF

        
     
      IF(IPICO.NE.0) THEN
        E0VAR = dbl_mb(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        dbl_mb(KAPROJ) = E0PERT
        if (ga_nodeid().eq.0) then
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
        endif
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
      IF(IPRT .GE. 3.and.ga_nodeid().eq.0 ) THEN
        WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(6,'(5F22.13)')
     &  ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
      END IF
*
      IF( IPRT  .GE. 5 .and. ga_nodeid().eq.0) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 1
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3 .AND.NVEC.GE.2.
     &   .AND. .NOT.CONVER) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
CNW     CALL REWINO( LU3)
CNW     CALL REWINO( LU1)
*. Current trial vector
CNW     CALL SETVEC(VEC1,0.0D0,NVAR)
        call ga_zero(VEC1)
        DO 2200 IVEC = 1, NVEC
CNW       CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU1,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
          FACTOR =  AVEC(IVEC)
CNW       CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
          call ga_add(FACTOR,VEC2,1.0d0,VEC1,VEC1)
 2200   CONTINUE
CNW     SCALE = INPROD(VEC1,VEC1,NVAR)
        SCALE = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                        VEC1,'N',1,NVAR,1,1)
        SCALE  = 1.0D0/SQRT(SCALE)
CNW     CALL SCALVE(VEC1,SCALE,NVAR)
        call ga_scale(VEC1,SCALE)
CNW     CALL TODSC(VEC1,NVAR,-1  ,LU3)
        if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &      1,1,itask).ne.0) call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
* Previous trial vector orthonormalized
CNW     CALL REWINO(LU1)
CNW     CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
        if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU1,1,NVAR,
     &      1,1,itask).ne.0) call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
CNW     OVLAP = INPROD(VEC1,VEC2,NVAR)
        OVLAP = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,
     &                        VEC2,'N',1,NVAR,1,1)
CNW     CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        call ga_add(-OVLAP,VEC1,1.0d0,VEC2,VEC2)
CNW     SCALE2 = INPROD(VEC2,VEC2,NVAR)
        SCALE2 = ga_ddot_patch(VEC2,'N',1,NVAR,1,1,
     &                         VEC2,'N',1,NVAR,1,1)
        SCALE2 = 1.0D0/SQRT(SCALE2)
CNW     CALL SCALVE(VEC2,SCALE2,NVAR)
        call ga_scale(VEC2,SCALE2)
CNW     CALL TODSC(VEC2,NVAR,-1,LU3)
        if (dra_write_section(.false.,VEC2,1,NVAR,1,1,LU3,1,NVAR,
     &      2,2,itask).ne.0) call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
*
CNW     CALL REWINO( LU1)
CNW     CALL REWINO( LU3)
        DO 2411 IVEC = 1,2
          if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
          if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
CNW       CALL TODSC (VEC1,NVAR,-1,  LU1)
 2411   CONTINUE
*. Corresponding sigma vectors
CNW     CALL REWINO ( LU3)
CNW     CALL REWINO( LU2)
CNW     CALL SETVEC(VEC1,0.0D0,NVAR)
        call ga_zero(VEC1)
        DO 2250 IVEC = 1, NVEC
          if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
CNW       CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
          call ga_add(FACTOR,VEC2,1.0d0,VEC1,VEC1)
 2250   CONTINUE
*
CNW     CALL SCALVE(VEC1,SCALE,NVAR)
        call ga_scale(VEC1,SCALE)
CNW     CALL TODSC(VEC1,NVAR,-1,  LU3)  
         if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &       1,1,itask).ne.0) call errquit('dra error',itask,911)
         if (dra_wait(itask).ne.0) call errquit('dra read err',
     &       itask,911)
* Sigma vector corresponding to second vector on LU1
CNW     CALL REWINO(LU2)
CNW     CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
        if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &      1,1,itask).ne.0) call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
CNW     CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        call ga_add(-OVLAP,VEC1,1.0d0,VEC2,VEC2)
CNW     CALL SCALVE(VEC2,SCALE2,NVAR)
        call ga_scale(VEC2,SCALE2)
CNW     CALL TODSC(VEC2,NVAR,-1,LU3)
        if (dra_write_section(.false.,VEC2,1,NVAR,1,1,LU3,1,NVAR,
     &      2,2,itask).ne.0) call errquit('dra error',itask,911)
        if (dra_wait(itask).ne.0) call errquit('dra read err',
     &      itask,911)
*
CNW     CALL REWINO( LU2)
CNW     CALL REWINO( LU3)
        DO 2400 IVEC = 1,2
          if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU3,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
          if (dra_write_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
CNW       CALL TODSC (VEC2,NVAR,-1  ,LU2)
 2400   CONTINUE
        NVEC = 2
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
CNW    CALL REWINO( LU1 )
       DO 2010 IVEC = 1,NVEC
         if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &       IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
         if (dra_wait(itask).ne.0) call errquit('dra read err',
     &       itask,911)
CNW      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
CNW      CALL REWINO( LU2)
         DO 2008 JVEC = 1, IVEC
CNW        CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,IVEC,
     &         IVEC,itask).ne.0) call errquit('dra error',itask,911)
           if (dra_wait(itask).ne.0) call errquit('dra read err',
     &         itask,911)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = ga_ddot(VEC1,VEC2)
CNW        APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
CNW     CALL REWINO( LU3)
        DO 320 IROOT = 1, NROOT
CNW       CALL REWINO( LU1)
CNW       CALL SETVEC(VEC1,0.0D0,NVAR)
          call ga_zero(VEC1)
          DO 200 IVEC = 1, NVEC
            if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU1,1,NVAR,
     &      IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
            if (dra_wait(itask).ne.0) call errquit('dra read err',
     &          itask,911)
CNW         CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
CNW         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
            call ga_add(FACTOR,VEC2,1.0d0,VEC1,VEC1)
  200     CONTINUE
*
CNW       SCALE = INPROD(VEC1,VEC1,NVAR)
          SCALE = ga_ddot_patch(VEC1,'N',1,NVAR,1,1,VEC1,'N',1,NVAR,1,1)
          SCALE  = 1.0D0/SQRT(SCALE)
CNW       CALL SCALVE(VEC1,SCALE,NVAR)
          call ga_scale(VEC1,SCALE)
CNW       CALL TODSC(VEC1,NVAR,-1  ,LU3)
          if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &    IROOT,IROOT,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
  320   CONTINUE
CNW     CALL REWINO( LU1)
CNW     CALL REWINO( LU3)
        DO 411 IVEC = 1,NROOT
          if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
          if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
CNW       CALL TODSC (VEC1,NVAR,-1,  LU1)
  411   CONTINUE
* CORRESPONDING SIGMA VECTOR
CNW     CALL REWINO ( LU3)
        DO 329 IROOT = 1, NROOT
CNW       CALL REWINO( LU2)
CNW       CALL SETVEC(VEC1,0.0D0,NVAR)
          call ga_zero(VEC1)
          DO 250 IVEC = 1, NVEC
            if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &       IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
            if (dra_wait(itask).ne.0) call errquit('dra read err',
     &          itask,911)
CNW         CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
CNW         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
            call ga_add(FACTOR,VEC2,1.0d0,VEC1,VEC1)
  250     CONTINUE
*
CNW       CALL SCALVE(VEC1,SCALE,NVAR)
          call ga_scale(VEC1,SCALE)
CNW       CALL TODSC(VEC1,NVAR,-1,  LU3)
          if (dra_write_section(.false.,VEC1,1,NVAR,1,1,LU3,1,NVAR,
     &    IROOT,IROOT,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
  329   CONTINUE
* PLACE C IN LU1 AND HC IN LU2
CNW     CALL REWINO( LU2)
CNW     CALL REWINO( LU3)
        DO 400 IVEC = 1,NROOT
          if (dra_read_section(.false.,VEC2,1,NVAR,1,1,LU3,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
          if (dra_write_section(.false.,VEC2,1,NVAR,1,1,LU2,1,NVAR,
     &        IVEC,IVEC,itask).ne.0) call errquit('dra error',itask,911)
          if (dra_wait(itask).ne.0) call errquit('dra read err',
     &        itask,911)
CNW       CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
CNW       CALL TODSC (VEC2,NVAR,-1  ,LU2)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
C
      END IF
C
C     IF( ITER .LT. MAXIT .AND. .NOT. CONVER) GOTO 1000
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*. Place first eigenvector in vec1
CNW   CALL REWINO(LU1)
CNW   CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      if (dra_read_section(.false.,VEC1,1,NVAR,1,1,LU1,1,NVAR,1,1,
     &    itask).ne.0) call errquit('dra error',itask,911)
      if (dra_wait(itask).ne.0) call errquit('dra read err',
     &    itask,911)

* ( End of loop over iterations )
*
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 .and. ga_nodeid().eq.0)
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 .and. ga_nodeid().eq.0)
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*. Final eigenvalues
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           RNRM_CNV(IROOT) = RNRM(ITER,IROOT)
 1601   CONTINUE
*
      IF ( IPRT .GT. 1 ) THEN
        DO 1600 IROOT = 1, NROOT
          if (ga_nodeid().eq.0) then
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          endif
          IF(IPRT.GE.1000) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
cVOG        CALL REWINO( LU1)
cVOG        CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
          END IF
          if (ga_nodeid().eq.0) then
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
          endif
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          if(ga_nodeid().eq.0) WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
 1234 CONTINUE
C
      CALL LUCIAQEXIT('MINDV')
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE MINGCG(MV8,LU1,LU2,LU3,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPROD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CONVER = .FALSE.
      ITER = 1
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X
      CALL REWINO(LU1)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINO(LU2)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-1.0D0,NVAR)
*
      RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
      ERROR(1) = RNORM
      CALL REWINO(LU2)
      CALL TODSC(VEC1,NVAR,-1  ,LU2)
*. Preconditioner H times initial vector , H * R
*.H * R
      CALL REWINO(LUDIA)
      CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
      CALL DIAVC2(VEC2,VEC1,VEC2,W,NVAR)
      IF(LUPROJ.NE.0) THEN
        CALL REWINO(LUPROJ)
        CALL FRMDSC(VEC1,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL REWINO(LU2)
        CALL FRMDSC(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPROD(VEC1,VEC2,NVAR)
*. P = RHO * H*R
      RHO = 1.0D0
      CALL SCALVE(VEC2,RHO,NVAR)
      CALL REWINO(LU3)
      CALL TODSC(VEC2,NVAR,-1  ,LU3)
      CALL COPVEC(VEC2,VEC1,NVAR)
*.S = AP
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINO (LU3)
      CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*.    P is assumed in VEC1 and S = A*P in VEC2
 
        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPROD(VEC1,VEC2,NVAR)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S
        CALL REWINO(LU2)
        CALL FRMDSC(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-A,NVAR)
        CALL REWINO(LU2)
        CALL TODSC(VEC1,NVAR,-1  ,LU2)
*.    new residual has been obtained , check for convergence
        RNORM = INPROD(VEC1,VEC1,NVAR)
        ERROR(ITER+1) = SQRT(RNORM)
*.    X = X + A * P
        CALL REWINO(LU1)
        CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL REWINO(LU3)
        CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC2,VEC1,1.0D0,A,NVAR)
        CALL REWINO(LU1)
        CALL TODSC(VEC1,NVAR,-1  ,LU1)
*
        IF( SQRT(RNORM) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.       H * R
           CALL REWINO(LU2)
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           CALL REWINO(LUDIA)
           CALL FRMDSC(VEC1,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL DIAVC2(VEC1,VEC2,VEC1 ,W,NVAR)
           IF(LUPROJ.NE.0) THEN
             CALL REWINO(LUPROJ)
             CALL FRMDSC(VEC2,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
             CALL REWINO(LU2)
             CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
           END IF
           GAMMA = INPROD(VEC1,VEC2,NVAR)
           B = GAMMA/C
*.       P = RHO*(H*R + B*P)
           CALL REWINO(LU3)
           CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,B,NVAR)
*.       Define next RHO
           RHO = 1.0D0
           CALL SCALVE(VEC1,RHO,NVAR)
           CALL REWINO(LU3)
           CALL TODSC(VEC1,NVAR,-1  ,LU3)
*.       S = MATRIX * P
           CALL MV8(VEC1,VEC2,0,0)
           CALL REWINO(LU3)
           CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
*.End of prepations for next iteration
        END IF
*
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',F13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT +1 )
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',F13.8)
      END IF
      END IF
C
      IF(NTEST.GT. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL REWINO(LU1)
       CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
       CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call wrtmat(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
C
      RETURN
      END 
      SUBROUTINE MINPRD(VU,A,VI,IP,NPROD,NROW)
*
* VU(I) = SUM(J) A(J,IP(I))*VI(J)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VU(*),A(NROW,*),VI(*),IP(*)
*. Loop structure for IBM 3090
      CALL SETVEC(VU,0.0D0,NPROD)
      DO 50 J = 1, NROW
      DO 100 I = 1, NPROD
          VU(I) = VU(I) + A(J,IP(I))*VI(J)
  100 CONTINUE
   50 CONTINUE
*
      RETURN
      END
      SUBROUTINE MSAXPY(AX,A,X,TEST,NDIM,NVEC,INDEX,NVCEFF)
*
* AX(I) = SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
*
* New version with seperate treatment of small loop lengths
* IBM 3090 VERSION
*
* Jeppe Olsen , Spring of 1990
*
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AX(*),X(NDIM,*)
      DIMENSION A(*) ,INDEX(*)
*
      IF(NDIM.EQ.1) THEN
*. Loop length 1
        X1 = 0.0D0
        DO 11 L = 1, NVCEFF
          X1 = X1 + A(L)*X(1,INDEX(L))
   11   CONTINUE
        AX(1) = X1
        RETURN
      ELSE IF(NDIM.EQ.2) THEN
*. Loop length 2
        X1 = 0.0D0
        X2 = 0.0D0
        DO 12 L = 1, NVCEFF
          X1 = X1 + A(L)*X(1,INDEX(L))
          X2 = X2 + A(L)*X(2,INDEX(L))
   12   CONTINUE
        AX(1) = X1
        AX(2) = X2
        RETURN
      ELSE IF(NDIM.EQ.3) THEN
*. Loop length 3
        X1 = 0.0D0
        X2 = 0.0D0
        X3 = 0.0D0
        DO 13 L = 1, NVCEFF
          X1 = X1 + A(L)*X(1,INDEX(L))
          X2 = X2 + A(L)*X(2,INDEX(L))
          X3 = X3 + A(L)*X(3,INDEX(L))
   13   CONTINUE
        AX(1) = X1
        AX(2) = X2
        AX(3) = X3
        RETURN
      ELSE IF(NDIM.EQ.4) THEN
*. Loop length 4
        X1 = 0.0D0
        X2 = 0.0D0
        X3 = 0.0D0
        X4 = 0.0D0
        DO 14 L = 1, NVCEFF
          X1 = X1 + A(L)*X(1,INDEX(L))
          X2 = X2 + A(L)*X(2,INDEX(L))
          X3 = X3 + A(L)*X(3,INDEX(L))
          X4 = X4 + A(L)*X(4,INDEX(L))
   14   CONTINUE
        AX(1) = X1
        AX(2) = X2
        AX(3) = X3
        AX(4) = X4
        RETURN
      ELSE IF( NDIM .GE.5) THEN
*. Loop length atleast 5
        DO 100 I = 1, NDIM
          T = 0.0D0
          DO 80 L = 1,NVCEFF
            T = T + A(L)*X(I,INDEX(L))
   80     CONTINUE
          AX(I) = T
  100   CONTINUE
        RETURN
      END IF
*
      END
 
 
 
      SUBROUTINE MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
C
C ADD VECTORS ON FILE LUIN TIMES FACTOR AND STORE ON LUOUT
C
C LUOUT AND LUSCR ARE INITIALLY REWINDED
C
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1)
C
      IF( MOD(NVEC,2) .EQ. 0 ) THEN
        LLUOUT = LUSCR
        LLUSCR = LUOUT
      ELSE
        LLUOUT = LUOUT
        LLUSCR = LUSCR
      END IF
C
      IF(IREW .NE. 0 ) CALL REWINE(LUIN,LBLK)
C
      DO 100 IVEC = 1, NVEC
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
        IF( IVEC .EQ. 1 ) THEN
          CALL SCLVCD(LUIN,LLUOUT,FAC(IVEC),VEC1,0,LBLK)
        ELSE
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),1.0D0,LUIN,LLUSCR,LLUOUT,
     &                0,LBLK)
        END IF
C
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE NEWDIR(D,X,G,DIAG,E,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &                  INVCOR,WORK)
*
* Calculate
*
*  D = (H0-E)** (-1) * G    (INVCOR = 0 )
*
*  D = (H0-E)** (-1) * G - ALPHA * (H0 - E)**(-1) * X
*
*       ALPHA = X(T)(H0-E)**(-1)*D / X(T)(H0-E)**(-1)*X  (INVCOR .NE.0)
*
* The latter correction corresponds to inverse iteration
* correction to Davidson
*
* Where H0 consists of a diagonal Diag
* and a block matrix of dimension NPRDIM.
*
* The block matrix is defined by
* ==============================
*
*  NPRDIM : Size of block matrix
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*  PEIGVL   : Eigenvalues of subblock mateix
*  PEIGVC   : Eigenvectors of subblock matrix
*
* Input
*=======
* X : for eigenvalue problem X is current eigenvector
*     (for INVCOR = 0 X can be a dummy variable )
* G : for eigenvalue problem G = (H - E ) * X
* Diag : Diaginal of matrix
* E : Energy for shift
* NVAR : Dimension of full matrix
* NPRDIM,IPNTR,PEIGVL,PEIGVC : See above
* INVCOR : use(.NE.0) , do not use (.eq.0) inverse correction
* Modification
* Work : Scratch space , at least ??
*
* Output
* ======
* D as given above, code has been constructed so D
* can occupy the same place as either X,G,DIAG
*
* Scratch space
*===============
* Should at least be of length ???
*
* Externals   GPRCTV,INPROD
*===========
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 INPROD
*
      DIMENSION D(*)
      DIMENSION X(*),G(*),DIAG(*)
      DIMENSION IPNTR(*),PEIGVL(*),PEIGVC(*)
      DIMENSION WORK(*)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Information from NEWDIR '
        WRITE(6,*) ' ========================'
      END IF
      IF( INVCOR .EQ. 0 ) THEN
* (H0 - E ) **(-1) * G , store in D
C       SUBROUTINE GPRCTV(DIAG,VECIN,VECUT,NVAR,NPRDIM,IPNTR,
C    &                    PEIGVL,PEIGVC,SHIFT,WORK)
        CALL GPRCTV(DIAG,G,D,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &              -E,WORK,XDUMMY)
      ELSE
* (H0 - E ) **(-1) * G , store in G
        CALL GPRCTV(DIAG,G,G,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &              -E,WORK,XDUMMY)
* X(T) (H0 - E) ** (-1) X
        XH0MEG = INPROD(X,G,NVAR)
C?      write(6,*) ' XH0MEG ', XH0MEG
* (H0 - E ) **(-1) * X , store in X
        CALL GPRCTV(DIAG,X,X,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &              -E,WORK,XH0MEX)
C?      write(6,*) ' XH0MEX ', XH0MEX
*
        FACTOR = -XH0MEG/XH0MEX
        CALL VECSUM(D,G,X,1.0D0,FACTOR,NVAR)
C?      write(6,*) ' New direction '
C?      call wrtmat(D,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
      SUBROUTINE ONEMAT(A,B,NBAS,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      ONEMAT PACK THE UPPER HALF OF A TWO DIM MATRIX
C      A INTO A ONE DIM MATRIX B
      DIMENSION A(NBAS,1),B(1)
      DO 100 I=1,N
      DO 200 J=1,I
      IJ=I*(I-1)/2 + J
200   B(IJ)=A(J,I)
100   CONTINUE
      RETURN
      END
      SUBROUTINE ORTVCD(LUIN,LUVEC,LUOUT,LUSCR,VEC1,VEC2,NVEC,LBLK,
     &                  SCR,INORMA)
*
* Orthonormalize vector on file LUIN to NVEC vectors on file LUVEC
* and save result on file LUOUT.
* If INORMA .ne. 0 the vector is normalized
* The transformation vector is returned in SCR
*
*. All files are rewinded
 
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      REAL*8 INPRDD
*.Scratch
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION SCR(1)
*
      IF(INORMA.NE.0) THEN
        IF( MOD(NVEC,2) .EQ. 0 ) THEN
          LLUOUT = LUSCR
          LLUSCR = LUOUT
        ELSE
          LLUOUT = LUOUT
          LLUSCR = LUSCR
        END IF
      ELSE IF( INORMA.EQ.0) THEN
        IF( MOD(NVEC,2) .EQ. 1 ) THEN
          LLUOUT = LUSCR
          LLUSCR = LUOUT
        ELSE
          LLUOUT = LUOUT
          LLUSCR = LUSCR
        END IF
      END IF
*.Pass 1 : Obtain overlap vector
      CALL REWINE(LUVEC,LBLK)
      DO 200 IVEC = 1, NVEC
        CALL REWINE(LUIN,LBLK)
        SCR(IVEC) = INPRDD(VEC1,VEC2,LUVEC,LUIN,0,LBLK)
  200 CONTINUE
* Pass 2 : Orthogonalize
      CALL COPVCD(LUIN,LLUOUT,VEC1,1,LBLK)
      LBUF = LLUOUT
      LLUOUT = LLUSCR
      LLUSCR = LBUF
      CALL REWINE(LUVEC,LBLK)
      DO 100 IVEC = 1, NVEC
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
        CALL VECSMD(VEC1,VEC2,SCR(IVEC),1.0D0,LUVEC,LLUSCR,LLUOUT,
     &                0,LBLK)
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
  100 CONTINUE
*
      IF(INORMA.NE.0) THEN
        XNORM = INPRDD(VEC1,VEC1,LLUSCR,LLUSCR,1,LBLK)
        FACTOR = 1.0D0/SQRT(XNORM)
C            SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
        CALL SCLVCD(LLUSCR,LLUOUT,FACTOR,VEC1,1,LBLK)
        CALL SCALVE(SCR,FACTOR,NVEC)
      END IF
*
      RETURN
      END
      SUBROUTINE OUTPAK(MATRIX,NROW,NCTL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...........VERSION = 09/16/73/04
C.......................................................................
C
C OUTPAK PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        MATRIX(*)...........PACKED MATRIX
C        NROW................NUMBER OF ROWS TO BE OUTPUT
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C   2 FOR DOUBLE SPACE,
C   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C
C AND SO ON.
C
C OUTPAK IS SET UP TO HANDLE 8 COLUMNS/PAGE WITH A 8F15.7 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE, FLORIDA, AND DIVISION OF THEORETICAL
C          CHEMISTRY, DEPARTMENT OF CHEMISTRY, AARHUS UNIVERSITY,
C          AARHUS, DENMARK
C
C.......................................................................
      INTEGER BEGIN,ASA,BLANK,CTL
      LOGICAL HEADER
      DOUBLE PRECISION  MATRIX
      DIMENSION MATRIX(1),ASA(3)
      DATA KCOL/8/, COLUMN/8HCOLUMN  /, ASA/4H    ,4H0   ,4H-   /,
     X     BLANK/4H    /, ZERO/0.D+00/
      CTL = BLANK
      IF ((NCTL.LE.3).AND.(NCTL.GT.0)) CTL = ASA(NCTL)
C.......................................................................
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
C.......................................................................
      LAST = MIN(NROW,KCOL)
C.......................................................................
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINT_D.
C
C.......................................................................
      BEGIN = 1
  100 NCOL = 1
      NCOL = 1
      HEADER = .TRUE.
           DO 500 K = BEGIN,NROW
           KTOTAL = (K*(K-1))/2 + BEGIN - 1
                DO 200 I = 1,NCOL
                IF (MATRIX(KTOTAL+I) .NE. ZERO) GO TO 300
  200           CONTINUE
           GO TO 400
  300      IF (HEADER) WRITE (6,10000) (COLUMN,I, I = BEGIN,LAST)
           HEADER = .FALSE.
           WRITE (6,20000) CTL,K,(MATRIX(KTOTAL+I), I = 1,NCOL)
  400      IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
  500      CONTINUE
      LAST = MIN(LAST+KCOL,NROW)
      BEGIN = BEGIN + NCOL
      IF (BEGIN .LE. NROW) GO TO 100
      RETURN
10000 FORMAT (1H0,8X,8(5X,A6,I4))
20000 FORMAT (A1,4H ROW,I4,8F15.7)
      END
      SUBROUTINE LUCIAOUTPUT (MATRIX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     X          COLDIM,
     X          NCTL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...........VERSION = 09/16/73/03
C.......................................................................
C
C OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        MATRIX(*,*).........MATRIX TO BE OUTPUT
C        ROWLOW..................ROW NUMBER AT WHICH OUTPUT IS STARTED
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
C        ROWDIM..............ROW DIMENSION OF MATRIX(*,*)
C        COLDIM..............COLUMN DIMENSION OF MATRIX(*,*)
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE
C   2 FOR DOUBLE SPACE
C   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*6*4.  THE
C PROGRAM IS SET UP TO HANDLE 8 COLUMNS/PAGE WITH A 8F15.7 FORMAT FOR
C THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE, FLORIDA, AND DIVISION OF THEORETICAL
C          CHEMISTRY, DEPARTMENT OF CHEMISTRY, AARHUS UNIVERSITY,
C          AARHUS, DENMARK
C
C.......................................................................
      DOUBLE PRECISION  MATRIX,COLUMN
      INTEGER ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,ASA,BLANK,
     X          CTL
      LOGICAL HEADER
      DIMENSION MATRIX(ROWDIM,COLDIM),ASA(3)
      DATA KCOL/8/, COLUMN/8HCOLUMN  /, ASA/4H    ,4H0   ,4H-   /,
     X     BLANK/4H    /, ZERO/0.D+00/
      CTL = BLANK
      IF ((NCTL.LE.3).AND.(NCTL.GT.0)) CTL = ASA(NCTL)
      IF (ROWHI .LT. ROWLOW) GO TO 500
      IF (COLHI .LT. COLLOW) GO TO 500
      LAST = MIN(COLHI,COLLOW+KCOL-1)
           DO 400 BEGIN = COLLOW,COLHI,KCOL
           HEADER = .TRUE.
            DO 300 K = ROWLOW,ROWHI
            DO 100 I = BEGIN,LAST
             IF (MATRIX(K,I) .NE. ZERO) GO TO 200
  100       CONTINUE
            GO TO 300
  200           IF (HEADER) WRITE(6,10000) (COLUMN,I, I = BEGIN,LAST)
                HEADER = .FALSE.
                WRITE(6,20000) CTL,K,(MATRIX(K,I), I = BEGIN,LAST)
  300           CONTINUE
  400      LAST = MIN(LAST+KCOL,COLHI)
  500 RETURN
10000 FORMAT (1H0,8X,8(5X,A6,I4))
20000 FORMAT (A1,4H ROW,I4,8F15.7)
      END
      SUBROUTINE PACKDI(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C       PACKDI COPY THE DIAGONAL ELEMENTS OF A INTO B
      DIMENSION A(1),B(1)
      DO 100 I=1,N
      II=I*(I+1)/2
100   B(I)=A(II)
      RETURN
      END
      SUBROUTINE PACKMT(A,B,NBAS,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      PACMMAT PACK A TWO DIM MATRIX THAT IS STORED
C      AS CULOMN VECTORS IN A ONE DIM ARRAY INTO
C      A TWO DIM MATRIX.
      DIMENSION B(NBAS,1),A(1)
      IQ=-N
      DO 100 I=1,N
      IQ=IQ+N
      DO 200 J=1,N
      IJ=IQ+J
200   B(J,I)=A(IJ)
100   CONTINUE
      RETURN
      END
 
      SUBROUTINE POSIFL(NREC,IFIL)
C
C POSITION FILE IFIL AT BEGINNING OF RECORD NREC
C
      CALL REWINO( IFIL      )
      ISKIP=NREC-1
      IF(ISKIP.NE.0) THEN
       DO 100 I=1,ISKIP
        READ(IFIL)
  100  CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE PRSYM_F7(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(6,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,10(1X,F7.3),/,(1H ,5X,10(1X,F7.3)))
      END
      SUBROUTINE PRSYM_EP(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
*
* Extended precision, E22.15
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(6,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,3(1X,E22.15),/,(1H ,5X,3(1X,E22.15)))
      END
      SUBROUTINE PRSYM(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(6,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,5(1X,E24.16),/,(1H ,5X,5(1X,E24.16)))
      END
      SUBROUTINE PRSYM_GEN(A,MATDIM,IROW_OR_COL)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
*
* IROW_OR_COL = 1 => Stored rowwise
* IROW_OR_COL = 2 => Stored columnwise
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
*
      IF(IROW_OR_COL.EQ.1) THEN
        JSTART=1
        JSTOP=0
        DO 100 I=1,MATDIM
          JSTART=JSTART+I-1
          JSTOP=JSTOP +I
          WRITE(6,1010) I,(A(J),J=JSTART,JSTOP)
  100   CONTINUE
      ELSE
        DO I = 1, MATDIM
          WRITE(6,1010)  I, (A((J-1)*MATDIM-J*(J-1)/2+I),J=1,I)
        END DO
      END IF
*
 1010 FORMAT(1H0,2X,I3,5(1X,E13.7),/,(1H ,5X,5(1X,E13.7)))
      RETURN
      END
      SUBROUTINE REWINE( LU ,LBLK )
*
* LBLK .LT. 0 :  REWIND SEQ FILE LU WITH FASTIO ROUTINES
* LBLK .GE. 0 :  rewinf seq file LU with normal REWIND
      ICRAY = 1
      IF ( ICRAY.EQ.0.AND.LBLK .LT. 0 ) THEN
        IDUM = 1
        CALL SQFILE(LU,5,IDUM,IDUM)
      ELSE
        REWIND LU
      END IF
*
      RETURN
      END
      SUBROUTINE REWINO( LU )
C
C REWIND SEQ FILE LU WITH FASTIO ROUTINES
C
C?    WRITE(6,*) ' TO REWIND FILE ',LU
      IDUM = 1
C     CALL SQFILE(LU,5,IDUM,IDUM)
      REWIND (LU)
C?    WRITE(6,*) ' FILE REWOUND '
C
      RETURN
      END
      SUBROUTINE SBINTV(NSBDIM,EIGVC,EIGVL,SHIFT,INDEX,VECI,VECO,X1,X2,
     &                  XHPSX)
*
* INVERTED SHIFTED SUBSPACE MATRIX TIMES VECTOR
*
* Last revision, oct 1989
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EIGVC(NSBDIM,NSBDIM),EIGVL(NSBDIM),INDEX(NSBDIM)
      DIMENSION X1(NSBDIM),X2(NSBDIM)
      DIMENSION VECI(1),VECO(1)
*
      CALL GATVEC(X1,VECI,INDEX,NSBDIM)
      CALL MATVCB(EIGVC,X1,X2,NSBDIM,NSBDIM,1)
      CALL DIAVC3(X1,X2,EIGVL,SHIFT,NSBDIM,XHPSX)
      CALL MATVCB(EIGVC,X1,X2,NSBDIM,NSBDIM,0)
      CALL SCAVEC(VECO,X2,INDEX,NSBDIM)
C
      NTEST = 0
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' OUTPUT FROM SBINTV, VECTOR IN GATHERED FORM '
        CALL WRTMAT(X1,1,NSBDIM,1,NSBDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE SCALE2(VECTOR,NDIM,SCALE)
C
C SCALE VECTOR TO HAVE NORM 1.VECTORS WITH ELEMENTS THAT CANNOT
C BE SQARED WITHOUT OVERFLOW  CAN BE HANDLED.SCALE FACTOR
C IS RETURNED IN SCALE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
C
C FIRST FIND GREATEST ELEMENT
      GREAT=0.0D0
      DO 100 I=1,NDIM
  100 IF(ABS(VECTOR(I)).GE.GREAT) GREAT=ABS(VECTOR(I))
C
C SCALE DOWN
      FACTOR=1.0D0/GREAT
      DO 200 I=1,NDIM
  200 VECTOR(I)=VECTOR(I)*FACTOR
C
C NORM OF SCALED VECTOR
      FACTOR=0.0D0
      DO 300 I=1,NDIM
       FACTOR=FACTOR+ VECTOR(I)**2
  300 CONTINUE
C
C THEN NORMALIZE
      FACTOR=DSQRT(FACTOR)
      DO 400 I=1,NDIM
  400 VECTOR(I)=VECTOR(I)/FACTOR
C
      SCALE=1.0D0/(FACTOR*GREAT)
C
      RETURN
      END
 
      SUBROUTINE SCALVE(VECTOR,FACTOR,NDIM)
C
C CALCULATE SCALAR(FACTOR) TIMES VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
      INCLUDE 'rou_stat.inc'
C     COMMON/ROU_STAT/NCALL_SCALVE,NCALL_SETVEC,NCALL_COPVEC,
C    &                NCALL_MATCG,NCALL_MATCAS,NCALL_ADD_SKAIIB,
C    &                NCALL_GET_CKAJJB,
C    &                XOP_SCALVE,XOP_SETVEC,XOP_COPVEC,
C    &                XOP_MATCG,XOP_MATCAS,XOP_ADD_SKAIIB,
C    &                XOP_GET_CKAJJB
C
      NCALL_SCALVE = NCALL_SCALVE + 1
      XOP_SCALVE = XOP_SCALVE + NDIM
*
      DO 100 I=1,NDIM
       VECTOR(I)=VECTOR(I)*FACTOR
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE SSCAVEC(VECO,VECI,INDEX,NDIM)
C
C SCATTER VECTOR with sign encoded
C VECO(ABS(INDEX(I)) = Sign(INDEX(I)*VECI(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECI(1   ),VECO(1),INDEX(1   )
C
      DO I = 1, NDIM
        IF(INDEX(I).GT.0) THEN
         VECO(INDEX(I)) = VECI(I)
        ELSE
         VECO(-INDEX(I)) = -VECI(I)
        END IF
      END DO
C
      RETURN
      END
      SUBROUTINE SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LUIN
          REWIND LUOUT
        ELSE
          CALL REWINE( LUIN ,LBLK)
          CALL REWINE( LUOUT,LBLK)
        END IF
      END IF
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF (LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
C
        IF ( LBL .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(LBL .GT. 0 )
     &    CALL SCALVE(SEGMNT,SCALE,LBL)
          CALL TODSC(SEGMNT,LBL,KBLK,LUOUT)
        END IF
C
      IF( LBL .GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
* Set diagonal elements of matrix MATRIX to VALUE
*
* IPACK = 0 => full quadratic matrix
* IPACK = 1 => lower triangular matrix, row packed
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MATRIX(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I=1,NDIM
100     MATRIX((I-1)*NDIM+I) = VALUE
      ELSE IF (IPACK .EQ. 1 ) THEN
        DO 200 I = 1, NDIM
 200    MATRIX(I*(I+1)/2) = VALUE
      ELSE
        WRITE(6,*) ' IPACK called with IPACK = ', IPACK
        STOP ' SETDIA ,IPACK out of range '
      END IF
*
      RETURN
      END
      SUBROUTINE SETVEC(VECTOR,VALUE,NDIM)
C
C VECTOR (*) = VALUE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VECTOR(2)
      INCLUDE 'rou_stat.inc'
C     COMMON/ROU_STAT/NCALL_SCALVE,NCALL_SETVEC,NCALL_COPVEC,
C    &                NCALL_MATCG,NCALL_MATCAS,NCALL_ADD_SKAIIB,
C    &                NCALL_GET_CKAJJB,
C    &                XOP_SCALVE,XOP_SETVEC,XOP_COPVEC,
C    &                XOP_MATCG,XOP_MATCAS,XOP_ADD_SKAIIB,
C    &                XOP_GET_CKAJJB
C
C
      NCALL_SETVEC = NCALL_SETVEC + 1
      XOP_SETVEC = XOP_SETVEC + NDIM
      DO 10 I=1,NDIM
   10 VECTOR(I) = VALUE
C
      RETURN
      END
      SUBROUTINE SETVECI(IVECTOR,IVALUE,NDIM)
C
C VECTOR (*) = VALUE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IVECTOR(2)
      INCLUDE 'rou_stat.inc'
C     COMMON/ROU_STAT/NCALL_SCALVE,NCALL_SETVEC,NCALL_COPVEC,
C    &                NCALL_MATCG,NCALL_MATCAS,NCALL_ADD_SKAIIB,
C    &                NCALL_GET_CKAJJB,
C    &                XOP_SCALVE,XOP_SETVEC,XOP_COPVEC,
C    &                XOP_MATCG,XOP_MATCAS,XOP_ADD_SKAIIB,
C    &                XOP_GET_CKAJJB
C
C
      NCALL_SETVEC = NCALL_SETVEC + 1
      XOP_SETVEC = XOP_SETVEC + NDIM
      DO 10 I=1,NDIM
   10 IVECTOR(I) = IVALUE
C
      RETURN
      END
 
      SUBROUTINE SKPRC3(IREC,IFILE)
C
C SKIP IREC RECORDS OF FILE IFILE
C
      DO 100 I=1,IREC
       READ(IFILE)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
C
C SKIP OVER NVEC VECTORS ON FILE LUIN
C
C LBLK DEFINES STRUCTURE OF FILE
C (see note on structure of files )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      NTEST = 00
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' SKPVCD: LU,NVEC,IREW,LBLK',LU,NVEC,IREW,LBLK
      IF( IREW .NE. 0 ) THEN
        CALL REWINE(LU ,LBLK)
      END IF
      DO 1001 IVEC = 1, NVEC
      IF(NTEST.GE.100) WRITE(6,*) ' Start IVEC = ', IVEC
C
C LOOP OVER BLOCKS OF GIVEN VECTOR
C
 1000   CONTINUE
C
          IF( LBLK .GT. 0 ) THEN
            LBL = LBLK
          ELSE IF (LBLK .EQ. 0 ) THEN
            READ(LU) LBL
          ELSE IF (LBLK .LT. 0 ) THEN
            CALL IFRMDS(LBL,1,-1,LU)
          END IF
C?        WRITE(6,*) ' LBL = ', LBL
C
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBLK
            ELSE
              KBLK = -1
            END IF
C?          WRITE(6,*) 'Before FRMDSC '
            CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,IAMPACK)
C?          WRITE(6,*) ' After Frmdsc '
          END IF
        IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C?      IF(NTEST.GE.100) WRITE(6,*) ' Stop  IVEC = ', IVEC
 1001 CONTINUE
C
      RETURN
      END
      SUBROUTINE SLRMTV(NMAT,NVAR,A,AVEC,NRANK,VECIN,VECOUT,IZERO,
     &                  DISCH,LUHFIL)
C CALCULATE PRODUCT OF MATRIX WITH VECTOR
C MATRIX IS DEFINED AS A SUM OF NMAT NRANK-MATRICES
C
C IF DISCH THEN VECTORS ARE ASSUMED STORED ON FILE LUHFIL. LENGTH
C OF AVEC MUST THEN AT LEAST BE NRANK*NVAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),AVEC(NVAR,*),VECIN(1   ),VECOUT(1   )
      LOGICAL DISCH
C
      IF ( DISCH ) REWIND LUHFIL
C
      DO 500 I = 1,NMAT
        IF( DISCH) THEN
          DO 400 IVEC = 1,NRANK
C           CALL SQFILE(LUHFIL,2,AVEC(1,IVEC),2*NVAR)
            READ(LUHFIL) (AVEC(II,IVEC),II=1,NVAR)
  400     CONTINUE
          IAVEC = 1
        ELSE
          IAVEC = (I-1)*NRANK + 1
        END IF
        IA = (I-1)*NRANK**2 + 1
        IF ( I.GT.1)  IZERO = 0
        CALL LRMTVC(NRANK,NVAR,A(IA),AVEC(1,IAVEC),VECIN,VECOUT,IZERO)
  500 CONTINUE
C
      NTEST = 0
      IF (NTEST.NE.0) THEN
       WRITE(6,*) ' MATRIX TIMES VECTOR FOR SLRMTVC'
       CALL RECPRT(VECOUT,1,NVAR)
C      CALL WRITVE(VECOUT,NVAR)
      END IF
C
      RETURN
      END
 
      SUBROUTINE SWAPVE(VEC1,VEC2,NDIM)
C
C      SWAP ELEMENTS OF VECTORS VEC1 AND VEC2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(1   ) ,VEC2(1   )
      DO 100 I=1,NDIM
       BUF=VEC1(I)
       VEC1(I)=VEC2(I)
       VEC2(I)=BUF
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE SYMTVC(A,VECIN,VECOUT,NDIM)
C
C INPUT :
C        A : LOWER HALF OF SYMMETRIC MATRIX A
C            A(I,J) = A((I(I-1)/2+J) (I.GE.J)
C        VECIN : A VECTOR
C        NDIM  : DIMENSION OF A
C OUTPUT :
C        VECOUT: A*VECIN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8   INPROD
      DIMENSION A(2),VECIN(2),VECOUT(2)
C
C** 1 : LOWER HALF TIMES VECTOR
C
      DO 100 I = 1,NDIM
  100 VECOUT(I) = INPROD(A((I-1)*I/2+1),VECIN(1),I)
C
C** 2 : UPPER HALF TIMES VECTOR
      DO 200 J = 1,NDIM
        JBASE = J*(J-1)/2
        VECINJ = VECIN(J)
        DO 190 I = 1,(J-1)
          VECOUT(I) = VECOUT(I)+ A(JBASE + I)*VECINJ
  190   CONTINUE
  200 CONTINUE
C
      NTEST = 0
      IF ( NTEST.GT.0) THEN
       WRITE(6,*) ' MATRIX TIMES VECTOR FROM SYMTVC '
       CALL WRTMAT(VECOUT,NDIM,1,NDIM,1)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
*   THIS IS A VERY STUPIDLY CODED PROGRAM FOR TRANSFORMING            *
*   A GENERALIZED EIGENVALUE PROBLEM INTO A NORMAL EIGENVALUE PROBLEM *
*                                                                     *
*   AUTHOR: M. MASAMURA                                               *
*           J.COMP.CHEM 9 (1988) 257.                                 *
*                                                                     *
*   THE ALGORITHM MIGHT BE USEFUL, BUT THE IMPLEMENTION IS FAR FROM   *
*   PREFECT.                                                          *
*                                                                     *
*   THIS CODE HAS BEEN ALMOST DIRECTLY COPIED FROM THE JOURNAL ABOVE  *
*   BY DAGE SUNDHOLM  (29.4.1988)                                     *
*                                                                     *
***********************************************************************
 
      SUBROUTINE TRANSH(N,H,S,P,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
 
C Symmetric matrices are assumed
C Transform H to H' obtain the transformation matrix P
C (HC=ESC) => (H'C'=EC') and (C=PC')
 
C N     : Dimension of the problem
C H     : Hamilton matrix, in H out H' (full matrix)
C S     : Overlap matrix,  in S out I  (full matrix)
C P     : Transformation matrix in trash out P (full matrix)
 
      DIMENSION S(N,N),H(N,N),P(N,N),WORK(N)
 
C Neglect matrix elements less than DEPS
 
      DEPS=0.5D-14
      ONE=1.0D0
 
C Set P to unit matrix
 
      CALL SETVEC(P,0.0D0, N ** 2 )
      CALL SETDIA(P,1.0D0,N,0)
*
C First part of the transformation of the H and P matrices
 
      DO 20 K=1,N-1
        DO 20 J=N,K+1,-1
 
          D=S(K,J)/S(K,K)
          IF(ABS(D).GT.DEPS) THEN
 
            DO 30 I=K+1,J
30          S(I,J)=S(I,J)-D*S(K,I)
            DO 31 I=K+1,J
31          H(I,J)=H(I,J)-D*H(K,I)
 
            DO 40 I=1,K
40          H(I,J)=H(I,J)-D*H(I,K)
 
            DO 50 I=J,N
50          H(J,I)=H(J,I)-D*H(K,I)
 
            DO 60 I=1,K
60          P(I,J)=P(I,J)-D*P(I,K)
 
          END IF
20    CONTINUE
 
C Second part of the transformation obtaining the final H and P matrices
C but just the upper triangle.
 
      DO 70 I=1,N
        E=SQRT(S(I,I))
 
        DO 80 J=1,N
80      H(I,J)=H(I,J)/E
 
        DO 90 J=1,I
90      H(J,I)=H(J,I)/E
 
        DO 100 J=1,I
        P(J,I)=P(J,I)/E
100     CONTINUE
70    CONTINUE
 
C To be sure, copy the upper triangle to the lower triangle
C set the S matrix to be unit matrix
C (Just in case)
 
      DO 200 I=1,N-1
        DO 200 J=I+1,N
200     H(J,I)=H(I,J)
*
      CALL SETVEC(P,0.0D0, N ** 2 )
      CALL SETDIA(P,1.0D0,N,0)
*
      RETURN
      END
      SUBROUTINE TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY =-1 : FULL TO PACKED + SYMMETRIZING
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO I = 1,NDIM
          DO J = 1, I
            APAK(IJ+J) = AUTPAK(J,I)
          END DO
          IJ = IJ + I
        END DO
      ELSE IF( IWAY .EQ. -1 ) THEN
        IJ = 0
        DO I = 1,NDIM
          DO J = 1, I
           APAK(IJ+J) = 0.5*(AUTPAK(J,I)+AUTPAK(I,J))
          END DO
          IJ = IJ + I
        END DO
      ELSE IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO I = 1,NDIM
          DO J = 1, I
           AUTPAK(I,J) = APAK(IJ+J)
           AUTPAK(J,I) = APAK(IJ+J)
         END DO
          IJ = IJ + I
        END DO
      ELSE
        STOP 'WHICH WAY? UNKNOWN IWAY IN TRIPAK!'
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPAK '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE UPTRIPAK(ATRI,AFUL,IWAY,NDIM,NDIMFUL)
c
c     switch between full matrix and upper triangular matrix:
c       iway = -1   pack and symmetrize
c       iway =  1   pack
c       iway =  2   unpack
c     the algorithm allows for in-place (un)packing, i.e. ATRI and
c     AFUL may have the same start address
c
      INCLUDE "implicit.inc"
      DIMENSION ATRI(*), AFUL(NDIMFUL,*)

      IF (IWAY.EQ.-1) THEN
        DO JJ = 1, NDIM
          IDXTRI = (JJ-1)*JJ/2
          DO II = 1, JJ
            ATRI(IDXTRI+II) = 0.5D0*(AFUL(II,JJ)+AFUL(JJ,II))
          END DO
        END DO
      ELSE IF (IWAY.EQ.1) THEN
        DO JJ = 1, NDIM
          IDXTRI = (JJ-1)*JJ/2
          DO II = 1, JJ
            ATRI(IDXTRI+II) = AFUL(II,JJ)
          END DO
        END DO
      ELSE IF (IWAY.EQ.2) THEN
        DO JJ = NDIM, 1, -1
          IDXTRI = (JJ-1)*JJ/2
          DO II = JJ, 1, -1
            AFUL(II,JJ) = ATRI(IDXTRI+II)
          END DO
        END DO
        DO JJ = 1, NDIM
          IDXTRI = (JJ-1)*JJ/2
          DO II = 1, JJ
            AFUL(JJ,II) = AFUL(II,JJ)
          END DO
        END DO
      ELSE
        WRITE(6,*) 'ILLEGAL VALUE FOR IWAY (',IWAY,')'
        STOP 'UPTRIPAK'
      END IF

      RETURN
      END      

        SUBROUTINE TRNMA2(A,X,SCRA,NDIM,MATDIM,itrans)
C
C TRANSFORM MATRIX A : if( itrans .eq.1 ) X(TRANS)*A*X
c                      if( itrans .eq.2 ) x * a * x(trans)
C A IS OVERWRITTEN BY TRANSFORMED A
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(MATDIM,1),X(MATDIM,1),
     +          SCRA(MATDIM,1)
C
C
      NTEST=1
C
      IF(NTEST.GE.3) THEN
       WRITE(16,1020)
 1020  FORMAT(1H0,'*** OUTPUT FROM TRANMAT')
       WRITE(16,1030)
 1030  FORMAT(1H0,'A- AND X-MATRIX')
C      CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
C      CALL WRTMAT(X,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      if ( itrans.eq.2) then
c sloopy transpose of x
      DO 2000 I = 1, NDIM
      DO 2000 J = 1, I
       BUF = X(I,J)
       X(I,J) = X(J,I)
       X(J,I) = BUF
 2000 CONTINUE
      END IF
C A*X
      DO 1000 I=1,NDIM
       DO 900 J=1,NDIM
       AX=0.0D0
        DO 800 K=1,NDIM
         AX=AX+A(I,K)*X(K,J)
  800   CONTINUE
        SCRA(I,J)=AX
  900  CONTINUE
 1000 CONTINUE
C
      IF(NTEST.GE.2) THEN
       WRITE(16,1040)
 1040  FORMAT(1H0,' AX MATRIX')
C      CALL WRTMAT(SCRA,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
C X(TRANS)*(A*X)
      DO 600 I=1,NDIM
       DO 500 J=1,NDIM
        XAX=0.0D0
        DO 400 K=1,NDIM
         XAX=XAX+X(K,I)*SCRA(K,J)
  400   CONTINUE
        A(I,J)=XAX
  500  CONTINUE
  600 CONTINUE
C
      if ( itrans.eq.2) then
c sloopy transpose of x
      DO 2100 I = 1, NDIM
      DO 2100 J = 1, I
       BUF = X(I,J)
       X(I,J) = X(J,I)
       X(J,I) = BUF
 2100 CONTINUE
      END IF
C
      IF(NTEST.GE.1) THEN
       WRITE(16,1010)
 1010  FORMAT(1H0,' TRANSFORMED MATRIX')
C      CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE TRNMAT(A,X,SCRA,NDIM,MATDIM)
C
C TRANSFORM MATRIX A : X(TRANS)*A*X
C A IS OVERWRITTREN BY TRANSFORMED A
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,1),X(MATDIM,1),
     +          SCRA(MATDIM,1)
C
C
      NTEST=0
C
      IF(NTEST.GE.3) THEN
       WRITE(6,1020)
 1020  FORMAT(1H0,'*** OUTPUT FROM TRANMAT')
       WRITE(6,1030)
 1030  FORMAT(1H0,'A- AND X-MATRIX')
       CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
       CALL WRTMAT(X,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
C A*X
      DO 1000 I=1,NDIM
       DO 900 J=1,NDIM
       AX=0.0D0
        DO 800 K=1,NDIM
         AX=AX+A(I,K)*X(K,J)
  800   CONTINUE
        SCRA(I,J)=AX
  900  CONTINUE
 1000 CONTINUE
C
      IF(NTEST.GE.2) THEN
       WRITE(6,1040)
 1040  FORMAT(1H0,' AX MATRIX')
       CALL WRTMAT(SCRA,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
C X(TRANS)*(A*X)
      DO 600 I=1,NDIM
       DO 500 J=1,NDIM
        XAX=0.0D0
        DO 400 K=1,NDIM
         XAX=XAX+X(K,I)*SCRA(K,J)
  400   CONTINUE
        A(I,J)=XAX
  500  CONTINUE
  600 CONTINUE
C
      IF(NTEST.GE.2) THEN
       WRITE(6,1010)
 1010  FORMAT(1H0,' TRANSFORMED MATRIX')
       CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
      SUBROUTINE TRNSPO(A,MATDIM,NDIM)
C
C       TRANSPOSE MATRIX A
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM)
      DO 100 I=1,NDIM
      DO 100 J=1,I-1
       BUF=A(I,J)
       A(I,J)=A(J,I)
       A(J,I)=BUF
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRPMAT(XIN,NROW,NCOL,XOUT)
C
C XOUT(I,J) = XIN(J,I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
C
      DO 200 IROW =1, NROW
        DO 100 ICOL = 1, NCOL
          XOUT(ICOL,IROW) = XIN(IROW,ICOL)
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
       SUBROUTINE TYMPAK(AIN,AOUT,NVAR)
C
C PACK SYMMETRIC MATRIX AIN TO LOWER TRIANGULAR FORM
C FOR REASON OF ADRESSING THE UPPER HALF OF AIN IS USED TO COPY FROM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AIN(NVAR,NVAR),AOUT(NVAR)
 
      IJ = 0
      DO 100 I = 1, NVAR
      DO 100 J = 1, I
        IJ = IJ + 1
        AOUT(IJ) = AIN(J,I)
  100 CONTINUE
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
       WRITE(6,*) ' MATRIX IN EXPANDED AND PACKED FORMAT '
       CALL WRTMAT(AIN,NVAR,NVAR,NVAR,NVAR)
       CALL PRSYM(AOUT,NVAR)
      END IF
C
      RETURN
      END
      SUBROUTINE UTPAK(A,SCR,NDIM,MATDIM,NNDIM)
C
C OUTPACK PACKED MATRIX A
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1    ),SCR(NDIM,NDIM)
C
      IJ=0
      DO 100 I=1,NDIM
      DO 100 J=1,I
       IJ=IJ+1
       SCR(I,J)=A(IJ)
  100 CONTINUE
C
      IJ=0
      DO 150 I=1,MATDIM
      DO 150 J=1,MATDIM
       IJ=IJ+1
       A(IJ)=0.0
  150 CONTINUE
C
      DO 200 I=1,NDIM
      DO 200 J=1,I
       A((J-1)*MATDIM+I)= SCR(I,J)
       A((I-1)*MATDIM+J)= SCR(I,J)
  200 CONTINUE
C
      RETURN
      END
      FUNCTION  VCSMDN(VEC1,VEC2,FAC1,FAC2,LU1,LU2,IREW,LBLK)
*
* Norm of sum of two vectors residing on disc
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPROD
*
      XNORM = 0.0D0
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
      END IF
*
* LOOP OVER BLOCKS OF VECTOR
*
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VCSMDN ',
     &  NBL1,NBL2
        STOP ' INCOMPATIBLE BLOCKSIZES IN VECSMF '
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        CALL FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
        CALL FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
        IF( NBL1 .GT. 0 ) THEN
          CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          XNORM = XNORM + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
*
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      VCSMDN = XNORM
      RETURN
      END
      SUBROUTINE VECSMDP(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
*
* Packed version, May 1996
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU3,LBLK)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
          WRITE(LU3) NBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
          CALL ITODS ( NBL1,1,-1,LU3)
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMDP',
     &  NBL1,NBL2
        STOP ' INCOMPATIBLE BLOCKSIZES IN VECSMF '
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        NO_ZEROING = 1
        CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
        CALL FRMDSC2(VEC2,NBL1,KBLK,LU2,IMZERO2,IAMPACK,NO_ZEROING)
        IF( NBL1 .GT. 0 ) THEN
          IF(IMZERO1.EQ.1.AND.IMZERO2.EQ.1) THEN
*. Simple zero record
            CALL ZERORC(NBL1,LU3,IAMPACK)
          ELSE
*. Nonvanishing record
            ZERO = 0.0D0
            IF(IMZERO1.EQ.1) THEN
              CALL VECSUM(VEC1,VEC1,VEC2,ZERO,FAC2,NBL1)
            ELSE IF(IMZERO2.EQ.1) THEN
              CALL VECSUM(VEC1,VEC1,VEC2,FAC1,ZERO,NBL1)
            ELSE
              CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
            END IF
            CALL TODSCP(VEC1,NBL1,KBLK,LU3)
          END IF
        ELSE IF (NBL1.EQ.0) THEN
          CALL TODSCP(VEC1,NBL1,KBLK,LU3)
        END IF
      END IF
C
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
C
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)        
        CALL REWINE( LU3,LBLK)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
          WRITE(LU3) NBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
          CALL ITODS ( NBL1,1,-1,LU3)
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I10)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &         NBL1,NBL2
          WRITE(6,'(A,2I3,A,I3,A)') 
     &              ' UNITS: ',LU1, LU2,'(IN) - ',LU3,' (OUT)'
          CALL UNIT_INFO(LU1)
          CALL UNIT_INFO(LU2)
          CALL UNIT_INFO(LU3)
          WRITE(6,*) 'CURRENT SEGMENT WAS ',IBLK
          STOP ' INCOMPATIBLE BLOCKSIZES IN VECSMD '
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF( NBL1 .GT. 0 )
     &         CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
        
          IF(IAMPACK.EQ.0) THEN
            CALL TODSC(VEC1,NBL1,KBLK,LU3)
          ELSE
            CALL TODSCP(VEC1,NBL1,KBLK,LU3)
          END IF
        END IF
C
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE VECSMe(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINO( LU1)
        CALL REWINO( LU2)
        CALL REWINO( LU3)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
       READ(LU1) NBL1
       READ(LU2) NBL2
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSME ',
     &  NBL1,NBL2
        STOP ' INCOMPATIBLE BLOCKSIZES IN VECSMF '
      END IF
C
      WRITE(LU3) NBL1
      IF(NBL1 .GE. 0 ) THEN
        CALL FRMDSC(VEC1,NBL1,-1  ,LU1,IMZERO,IAMPACK)
        CALL FRMDSC(VEC2,NBL1,-1  ,LU2,IMZERO,IAMPACK)
        IF( NBL1 .GT. 0 )
     &  CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
        CALL TODSC(VEC1,NBL1,-1  ,LU3)
      END IF
C
      IF(NBL1 .GE. 0 ) GOTO 1000
C
      RETURN
      END
      SUBROUTINE VECSMF(Q,V,SCRA,NVEC,IMULT,IADD,IVCFIL,NDIM)
C
C CALCULATE SUM OF VECTORS  RESIDING ON DISC.
C
C       Q(J)=SUM(IVEC) V(IVEC)*VECTOR(IVEC)(J)
C        IVEC=(I-1)*IMULT+IADD,I=1,NVEC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Q(1),V(1),SCRA(1)
C
      DO 1000 I=1,NVEC
       IF(I.EQ.1) THEN
        CALL POSIFL(IADD,IVCFIL)
       ELSE
        IF(IMULT.NE.1) CALL SKPRC3((IMULT-1),IVCFIL)
       END IF
       CALL FRMDSC(SCRA,NDIM,-1  ,IVCFIL,IMZERO,IAMPACK)
       CALL VECSUM(Q,Q,SCRA,1.0D0,V(I),NDIM)
 1000 CONTINUE
C
      RETURN
      END
      SUBROUTINE VECSUM(C,A,B,FACA,FACB,NDIM)
C
C     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1   ),B(1   ),C(1   )
*
      IF(FACA.NE.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 100 I=1,NDIM
          S=FACA*A(I)+FACB*B(I)
          C(I)=S
  100   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 200 I=1,NDIM
          S=FACB*B(I)
          C(I)=S
  200   CONTINUE
*
      ELSE IF(FACA.NE.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 300 I=1,NDIM
          S=FACA*A(I)
          C(I)=S
  300   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 400 I=1,NDIM
          C(I)=0.0D0
  400   CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE VTVTOV(AB,A,B,NDIM)
C AB(*) = A(*) * B(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1   ),B(1   ),AB(1   )
      DO 100 I = 1,NDIM
        AB(I) = A(I)*B(I)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE VVTOV(VECIN1,VECIN2,VECUT,NDIM)
C
C VECUT(I) = VECIN1(I) * VECIN2(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECIN1( 1  ),VECIN2(1   ),VECUT(1   )
C
      DO 100 I = 1, NDIM
        VECUT(I) = VECIN1(I) * VECIN2(I)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE WRITVE(VEC,NDIM)
      DOUBLE PRECISION  VEC
      DIMENSION VEC(1   )
C
      WRITE(6,1010) (VEC(I),I=1,NDIM)
 1010 FORMAT(1H0,2X,4(2X,E15.8),/,(1H ,2X,4(2X,E15.8)))
      RETURN
      END
      SUBROUTINE WRTDIA(A,NDIM,IFORM)
C
C PRINT DIAGONAL OF MATRIX A
C
C IFORM = 1 : MATRIX IS SQUARE PACKED
C IFORM = 2 : MATRIX IS LOWER TRIANGULAR PACKED
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
C
      IF( IFORM .EQ.1 ) THEN
        WRITE(6,'(4(2X,E14.8))')
     &  (A((I-1)*NDIM+I),I=1,NDIM)
      ELSEIF (IFORM .EQ. 2 ) THEN
        WRITE(6,'(4(2X,E14.8))')
     &  (A((I+1)*I/2),I=1,NDIM)
       END IF
C
      RETURN
      END
      SUBROUTINE WRTMAT_EP(A,NROW,NCOL,NMROW,NMCOL)
*
* Print matrix, extended precision (E25.15)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,2(1X,E25.15),/,(1H ,5X,2(1X,E25.15)))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE WRTMAT_F7(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,10(1X,F7.3),/,(1H ,5X,10(1X,F7.3)))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE WRTMAT(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,4(1X,E24.16),/,(1H ,5X,4(1X,E24.16)))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE WRTMAT2(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      ICOLMX=4
      ICOLL=0
      ICOLH=0
      DO WHILE (ICOLH.NE.NCOL)
        ICOLL = ICOLH+1
        ICOLH = MIN(ICOLL-1+ICOLMX,NCOL)
        WRITE(6,1000) (J,J=ICOLL,ICOLH)
        DO I=1,NROW
          WRITE(6,1010) I,(A(I,J),J=ICOLL,ICOLH)
        END DO
      END DO

      RETURN
 1000 FORMAT(1H0,3X,2X,4(1X,6X,I6,6X))
 1010 FORMAT(1H0,I3,2X,4(1X,E18.10))
      END
      SUBROUTINE WRTIMAT(IA,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(IA(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,6(1X,I10),/,(1H ,5X,6(1X,I10)))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE WRTVCD_EP(SEGMNT,LU,IREW,LBLK)
C
C PRINT VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,KBLK,LU,IMZERO,IAMPACK)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMAT_EP(SEGMNT,1,LBL,1,LBL)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE WRTVCD(SEGMNT,LU,IREW,LBLK)
C
C PRINT VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,KBLK,LU,IMZERO,IAMPACK)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      SUBROUTINE WRTVSD(SEGMNT,LU,IREW,LBLK)
C
C PRINT VECTOR STRUCTURE ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      REAL(8), EXTERNAL :: INPROD
C
      WRITE(6,*)
      WRITE(6,*) 'Structure of vector on unit ',lu
      CALL UNIT_INFO(LU)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,KBLK,LU,IMZERO,IAMPACK)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             WRITE(6,'(2(A,I3),A,E20.7)') ' zero_flag: ',IMZERO,
     &                     '  pack_flag: ',IAMPACK,
     &                     '  norm: ',SQRT(INPROD(SEGMNT,SEGMNT,LBL))
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
      FUNCTION XFAC(N)
*
* N !  as double precision real
*
      IMPLICIT REAL*8(A-H,O-Z)
      IF( N .LT. 0 ) THEN
       IFAC = 0
       WRITE(6,*) ' WARNING FACULTY OF NEGATIVE NUMBER SET TO ZERO '
      ELSE
C
       XFACN = 1.0D0
       DO 100 K = 2,N
        XFACN = XFACN * DBLE(K)
  100  CONTINUE
       XFAC = XFACN
      END IF
C
      RETURN
      END
c----------------------------------------------------------------------c
      SUBROUTINE SYMMAT(AMAT,NDIM,MAXDIM)
c----------------------------------------------------------------------c
c     symmetrize NDIMxNDIM block in matrix AMAT
c----------------------------------------------------------------------c
      include "implicit.inc"
      DIMENSION AMAT(MAXDIM,*)
      
      DO I = 1, NDIM
        DO J = 1, I-1
          ELM = 0.5*(AMAT(J,I)+AMAT(I,J))
          AMAT(J,I)=ELM
          AMAT(I,J)=ELM
        END DO
      END DO

      RETURN
      END
c----------------------------------------------------------------------c
      SUBROUTINE TEST_SYMMAT(AMAT,NDIM,MAXDIM)
c----------------------------------------------------------------------c
c     test NDIMxNDIM block in matrix AMAT on symmetry
c----------------------------------------------------------------------c
      include "implicit.inc"
      DIMENSION AMAT(MAXDIM,*)
      
      DO I = 1, NDIM
        DO J = 1, I-1
          THR = EPSILON(AMAT(J,I))
          IF (ABS(AMAT(J,I)-AMAT(I,J)).GT.THR) THEN
            WRITE(6,'(X,A,2I6,A,E12.6)')
     &           'Symmetry violation in pair ',I,J,' by ',
     &           ABS(AMAT(J,I)-AMAT(I,J))
          END IF
        END DO
      END DO

      RETURN
      END
c----------------------------------------------------------------------c
      SUBROUTINE LIST_SL(IMODE,VEC,NDIM,VECLIST,IVECLIST,NLIST)
c----------------------------------------------------------------------c
c get the NLIST smallest/largest (IMODE=1/2) vectors from VEC(NDIM) 
c and put them sorted into vector VECLIST(NLIST) (indices on IVECLIST)
c----------------------------------------------------------------------c
      INCLUDE "implicit.inc"
      DIMENSION VEC(NDIM), VECLIST(NLIST), IVECLIST(NLIST)

      ILIST=0
      XEXTR=0D0
      IF(IMODE.EQ.1) XEXTR=HUGE(XEXTR)

* Initialization cycles
      DO IDX = 1, NLIST
        XEL = VEC(IDX)
        VECLIST(IDX)=XEL
        IVECLIST(IDX)=IDX
        IF ((IMODE.EQ.1.AND.XEL.GT.XEXTR).OR.
     &      (IMODE.EQ.2.AND.XEL.LT.XEXTR)) THEN
          XEXTR = XEL
          IMAX = IDX
        END IF
      END DO
* Search for further small elements
      DO IDX = NLIST+1, NDIM
        XEL = VEC(IDX)
        IF ((IMODE.EQ.1.AND.XEL.LT.XEXTR).OR.
     &      (IMODE.EQ.2.AND.XEL.GT.XEXTR)) THEN
          VECLIST(IMAX) = XEL
          IVECLIST(IMAX) = IDX
          XEXTR = 0D0
          IF(IMODE.EQ.1) XEXTR=HUGE(XEXTR)
          DO JDX = 1, NLIST
            XEL = VECLIST(JDX)
            IF ((IMODE.EQ.1.AND.XEL.GT.XEXTR).OR.
     &          (IMODE.EQ.2.AND.XEL.LT.XEXTR)) THEN
              XEXTR = XEL
              IMAX = JDX
            END IF
          END DO
        END IF
      END DO
c sort the final list
      
      DO
        ISWAP = 0
        DO IDX = 2, NLIST
          IF ((IMODE.EQ.1.AND.VECLIST(IDX-1).GT.VECLIST(IDX)).OR.
     &        (IMODE.EQ.2.AND.VECLIST(IDX-1).LT.VECLIST(IDX)) ) THEN
            XHLP=VECLIST(IDX)
            VECLIST(IDX)=VECLIST(IDX-1)
            VECLIST(IDX-1)=XHLP
            IHLP=IVECLIST(IDX)
            IVECLIST(IDX)=IVECLIST(IDX-1)
            IVECLIST(IDX-1)=IHLP
            ISWAP = 1
          END IF          

        END DO

        IF (ISWAP.EQ.0) EXIT

      END DO

      RETURN

      END
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      SUBROUTINE LIST_ASL(IMODE,VEC,NDIM,VECLIST,IVECLIST,NLIST)
c----------------------------------------------------------------------c
c get the NLIST smallest/largest (IMODE=1/2) elements (abs. value) 
c from VEC(NDIM) and put them sorted into vector VECLIST(NLIST) 
c (indices on IVECLIST)
c----------------------------------------------------------------------c
      INCLUDE "implicit.inc"
      DIMENSION VEC(NDIM), VECLIST(NLIST), IVECLIST(NLIST)

      ILIST=0
      XEXTR=0D0
      IF(IMODE.EQ.2) XEXTR=HUGE(XEXTR)

* Initialization cycles
      DO IDX = 1, NLIST
        XEL = VEC(IDX)
        AXEL = ABS(XEL)
        VECLIST(IDX)=XEL
        IVECLIST(IDX)=IDX
        IF ((IMODE.EQ.1.AND.AXEL.GT.XEXTR).OR.
     &      (IMODE.EQ.2.AND.AXEL.LT.XEXTR)) THEN
          XEXTR = AXEL
          IMAX = IDX
        END IF
      END DO
* Search for further small elements
      DO IDX = NLIST+1, NDIM
        XEL = VEC(IDX)
        AXEL = ABS(XEL)
        IF ((IMODE.EQ.1.AND.AXEL.LT.XEXTR).OR.
     &      (IMODE.EQ.2.AND.AXEL.GT.XEXTR)) THEN
          VECLIST(IMAX) = XEL
          IVECLIST(IMAX) = IDX
          XEXTR = 0D0
          IF(IMODE.EQ.2) XEXTR=HUGE(XEXTR)
          DO JDX = 1, NLIST
            XEL = VECLIST(JDX)
            AXEL = ABS(XEL)
            IF ((IMODE.EQ.1.AND.AXEL.GT.XEXTR).OR.
     &          (IMODE.EQ.2.AND.AXEL.LT.XEXTR)) THEN
              XEXTR = AXEL
              IMAX = JDX
            END IF
          END DO
        END IF
      END DO
c sort the final list
      
      DO
        ISWAP = 0
        DO IDX = 2, NLIST
          IF ((IMODE.EQ.1.AND.
     &               ABS(VECLIST(IDX-1)).GT.ABS(VECLIST(IDX))).OR.
     &        (IMODE.EQ.2.AND.
     &               ABS(VECLIST(IDX-1)).LT.ABS(VECLIST(IDX))) ) THEN
            XHLP=VECLIST(IDX)
            VECLIST(IDX)=VECLIST(IDX-1)
            VECLIST(IDX-1)=XHLP
            IHLP=IVECLIST(IDX)
            IVECLIST(IDX)=IVECLIST(IDX-1)
            IVECLIST(IDX-1)=IHLP
            ISWAP = 1
          END IF          

        END DO

        IF (ISWAP.EQ.0) EXIT

      END DO

      RETURN

      END
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      REAL*8 FUNCTION FDMNXD(LUVE,MINMAX,SEGMNT,IREW,LBLK)
C
C FIND ELEMENT WITH SMALLEST (MINMAX==1) OR LARGEST (MINMAX==2) ABSOLUTE
C VALUE OF ELEMENTS OF VECTOR ON FILE LUVE 
C OR THE SMALLEST (MINMAX=-1) OR THE LARGEST (MINMAX=-2) ELEMENT
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      LOGICAL FIRST
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LUVE
        ELSE
          CALL REWINE( LUVE ,LBLK)
        END IF
      END IF
C
      IF (MINMAX.LT.-2.OR.MINMAX.GT.2.OR.MINMAX.EQ.0) THEN
        WRITE(6,*) 'Illegal parameter MINMAX in FDMNXD!'
        STOP       'Illegal parameter MINMAX in FDMNXD!'
      END IF
      FIRST=.TRUE.
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF (LBLK .EQ. 0 ) THEN
          READ(LUVE) LBL
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUVE)
        END IF
C
        IF ( LBL .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUVE,IMZERO,IAMPACK)
          IF(LBL .GT. 0 ) THEN
            IF (FIRST) THEN
              XMNX = ABS(SEGMNT(1))
              FIRST = .FALSE.
            END IF
            XMNXBLK = FNDMNX(SEGMNT,LBL,MINMAX)
            IF (ABS(MINMAX).EQ.1) XMNX = MIN(XMNX,XMNXBLK) 
            IF (ABS(MINMAX).EQ.2) XMNX = MAX(XMNX,XMNXBLK) 
          END IF
        END IF
C
      IF( LBL .GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      FDMNXD = XMNX
C
      RETURN
      END
c----------------------------------------------------------------------c
      SUBROUTINE CMP2VCD(VEC1,VEC2,LU1,LU2,THRSH,IREW,LBLK)
C
C DISC VERSION OF CMP2VC :
C
C      COMPARE BLOCKED VECTORS ON FILES LU1 AND LU2
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
      IBLK = 0
      NBL1 = 0
C
C loop over blocks
      DO 
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
        END IF
        IBLK = IBLK+1
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN CMP2VCD',
     &         NBL1,NBL2
          STOP ' INCOMPATIBLE BLOCKSIZES IN CMP2VCD '
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
          CALL FRMDSC2(VEC2,NBL1,KBLK,LU2,IMZERO2,IAMPACK,NO_ZEROING)
          IF( NBL1 .GT. 0 ) THEN
            WRITE(6,*) 'Current segment: ',IBLK,NBL1
            IF(IMZERO1.EQ.1.AND.IMZERO2.EQ.1) THEN
              WRITE(6,*) 'Segment is zero on both files'
            ELSE
*. Nonvanishing record
              ZERO = 0.0D0
              IF(IMZERO1.EQ.1) THEN
                CALL SETVEC(VEC1,ZERO,NBL1)
              ELSE IF(IMZERO2.EQ.1) THEN
                CALL SETVEC(VEC2,ZERO,NBL1)
              END IF
              CALL CMP2VC(VEC1,VEC2,NBL1,THRSH)
            END IF
          END IF
        END IF
C
        IF (.NOT.(NBL1.GE. 0 .AND. LBLK .LE. 0)) EXIT
C
      END DO
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine prtrlt(v,m)
      implicit real*8(a-h,o-z)
c
c     ----- print out the lower triangle of a symmetric matrix (stored
c           in packed canonical form (actually an upper triangle) !) -----
c
      dimension v(m*(m+1)/2)

      max=5
      imax = 0
      do while(imax.lt.m)
        imin = imax+1
        imax = min(imax+max,m)
        write(*,'(/,5x,10(6x,i4,5x)/)') (i,i = imin,imax)
        do i=1,m
          ii = i*(i-1)/2
          mm = imin + ii
          kk = min(i,imax) + ii
          if(mm.le.kk) then
            write(*,'(i4,1x,10e15.7)') i,(v(j),j=mm,kk)
          end if
        end do
      end do
      write(*,*)
      return
      end
      SUBROUTINE CMP2VSC(VEC1,VEC2,LU1,LU2,IREW,LBLK)
C
C      COMPARE STRUCTURE OF BLOCKED VECTORS ON FILES LU1 AND LU2
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
      IBLK = 0
      NBL1 = 0
C
C loop over blocks
      DO 
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
        END IF
        IBLK = IBLK+1
C
        IF (NBL1.EQ.-1.AND.NBL2.NE.-1.OR.
     &      NBL2.EQ.-1.AND.NBL1.NE.-1) THEN
          WRITE(6,*) 'Premature end of one vector: ',NBL1, NBL2
          RETURN
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
          WRITE(6,*) 'Current segment: ',IBLK,NBL1,NBL2
          IF( NBL1 .NE. NBL2 ) THEN
            WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES !',
     &         NBL1,NBL2
          END IF
          
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK1,NO_ZEROING)
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL2
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC2(VEC2,NBL2,KBLK,LU2,IMZERO2,IAMPACK2,NO_ZEROING)
          IF( NBL1 .GT. 0 ) THEN
            WRITE(6,*) 'Current segment is non-empty on: ',LU1,IBLK,NBL1
            IF(IMZERO1.EQ.1) THEN
              WRITE(6,*) 'Segment is zero on ', LU1
            END IF
            IF(IMPACK1.EQ.1) THEN
              WRITE(6,*) 'Segment is packed on ', LU1
            END IF
          END IF
          IF( NBL2 .GT. 0 ) THEN
            WRITE(6,*) 'Current segment is non-empty on: ',LU2,IBLK,NBL2
            IF(IMZERO2.EQ.1) THEN
              WRITE(6,*) 'Segment is zero on ', LU2
            END IF
            IF(IMPACK2.EQ.1) THEN
              WRITE(6,*) 'Segment is packed on ', LU2
            END IF
          END IF
        END IF
C
        IF (.NOT.(NBL1.GE. 0 .AND. LBLK .LE. 0)) EXIT
C
      END DO
C
      RETURN
      END

      REAL*8 FUNCTION INPRDD3(VEC1,VEC2,LU1,LU2,LU3,
     &                        SHIFT,XPOT,IREW,LBLK)
C
C CALC   X = sum_i f_i (m_i+shift)**xpot g_i
C
C LBLK DEFINES STRUCTURE OF FILE
C
*. Last revision, Sept 2003 : FRMDSC => FRMDSC2 to simplify handling 
*                             of vectors containing many zeo blocks
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
      LOGICAL DIFVEC
C
      X = 0.0D0
      IF( LU1 .NE. LU2 ) THEN
        DIFVEC = .TRUE.
      ELSE
        DIFVEC =  .FALSE.
      END IF
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          IF(DIFVEC) REWIND LU2
          REWIND LU3
         ELSE
          CALL REWINE( LU1,LBLK)
          IF( DIFVEC ) CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS OF VECTORS
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
          NBL3 = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          IF( DIFVEC) READ(LU2) NBL2
          READ(LU3) NBL3
        ELSE IF ( LBLK .LT. 0 ) THEN
          CALL IFRMDS(NBL1,1,-1,LU1)
          IF( DIFVEC)CALL IFRMDS(NBL2,1,-1,LU2)
          CALL IFRMDS(NBL3,1,-1,LU3)
        END IF
C
        NO_ZEROING = 1
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
C     FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,
C    &                   NO_ZEROING)
          CALL FRMDSC2(VEC2,NBL3,KBLK,LU3,IMZERO3,IAMPACK,NO_ZEROING)
          XPOTABS = ABS(XPOT)
          DO II = 1, NBL3
            VEC2(II) = (VEC2(II)+SHIFT)**XPOTABS
          END DO
          IF (XPOT.LT.0d0) THEN
            CALL DIAVC2(VEC2,VEC1,VEC2,0d0,NBL1)
          ELSE
            CALL VVTOV(VEC1,VEC2,VEC2,NBL1)
          END IF
          IF( DIFVEC) THEN
            CALL FRMDSC2(VEC1,NBL1,KBLK,LU2,IMZERO2,IAMPACK,
     &                   NO_ZEROING)
            IF(NBL1 .GT. 0 .AND. IMZERO1.EQ.0.AND.IMZERO2.EQ.0)
     &      X = X + INPROD(VEC1,VEC2,NBL1)
          ELSE
          IF(NBL1 .GT. 0 .AND. IMZERO1.EQ.0 )
     &    X = X + INPROD(VEC1,VEC2,NBL1)
        END IF
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      INPRDD3 = X
C
      RETURN
      END

      subroutine iptchma(imat,nlin,ncol,ilin,icol,ival)

      implicit none

      integer, intent(in) ::
     &     nlin,ncol,ilin,icol,ival
      integer, intent(inout) ::
     &     imat(nlin,ncol)

      imat(ilin,icol) = ival

      return
      end 
*----------------------------------------------------------------------*
*     follows: exp and log of matrices
*----------------------------------------------------------------------*
      subroutine expgmat(ndim,expx,xmat,xscr1,xscr2,thrsh)
*----------------------------------------------------------------------*
*     calculate exp(X), returned on expx, of (ndim,ndim)-matrix X,
*     input on xmat, by Taylor expansion (threshold thrsh)
*     xscr is a scratch matrix of the same dimensions as xmat, expx
*
*     any quadratic matrix may be supplied
*
*     andreas, aug 2004
*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter :: ntest = 100, maxn = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     thrsh
      real(8), intent(inout) ::
     &     expx(ndim,ndim), xmat(ndim,ndim),
     &     xscr1(ndim,ndim),  xscr2(ndim,ndim)

      logical ::
     &     conv
      integer ::
     &     n, ndim2, ii
      real(8) ::
     &     xnrm, fac

      real(8), external ::
     &     inprod

      expx(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)
      xscr2(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)

      do ii = 1, ndim
        expx(ii,ii) = expx(ii,ii) + 1d0
      end do

      ndim2 = ndim*ndim
      n = 1
      conv = .false.

      do while (.not.conv)
        n = n+1
        if (n.gt.maxn) exit

        fac = 1d0/dble(n)

        ! Xscr = 1/N Xscr * X
        call matml7(xscr1,xscr2,xmat,
     &              ndim,ndim,
     &              ndim,ndim,
     &              ndim,ndim,0d0,fac,0)

        xnrm = sqrt(inprod(xscr1,xscr1,ndim2))
        if (xnrm.lt.thrsh) conv = .true.

        if (ntest.ge.10)
     &       write(6,*) ' N = ',n,'  |1/N! X^N| = ',xnrm

        expx(1:ndim,1:ndim) = expx(1:ndim,1:ndim) + xscr1(1:ndim,1:ndim)
c        call vecsum(expx,expx,xscr1,1d0,1d0,ndim2)

        xscr2(1:ndim,1:ndim) = xscr1(1:ndim,1:ndim)
c        call copvec(xscr1,xscr2,ndim2)

      end do

      if (.not.conv) then
        write(6,*) ' Taylor expansion of exp(X) did not converge!'
        stop 'expgmat'
      end if

      return
      end
*----------------------------------------------------------------------*
      subroutine logumat(ndim,xlogx,xmat,xscr1,xscr2,xscr3)
*----------------------------------------------------------------------*
*     calculate the logarithm of a unitary matrix
*
*     the algorithm will use the eispack-routine rg() to calculate the
*     eigenvalues of the matrix which are decompose in modulus and angle.
*     the modulus should be one always, else the routine exits.
*
*     andreas, aug 2004
*     
*----------------------------------------------------------------------*

      implicit none

      integer, parameter :: ntest = 00, maxn = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     xlogx(ndim,ndim), xmat(ndim,ndim),
     &     xscr1(ndim,ndim), xscr2(ndim,ndim), xscr3(ndim,ndim)

      integer ::
     &     ii, ierr
      real(8) ::
     &     ang, xmod     ,ang1,ang2

* O(N) scratch
      real(8) ::
     &     eigr(ndim), eigi(ndim), scr(ndim)
      integer ::
     &     iscr(ndim)


      if (ntest.gt.0) then
        write(6,*) ' ==================== '
        write(6,*) '  LOGUMAT at work !!  '
        write(6,*) ' ==================== '

      end if

      if (ntest.ge.100) then
        write(6,*) ' xmat on entry:'
        call wrtmat2(xmat,ndim,ndim,ndim,ndim)
      end if

      xscr2(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)

      ! get eigenvalues and -vectors ...
      call rg(ndim,ndim,xscr2,
     &        eigr,eigi,1,xscr1,
     &        iscr,scr,ierr)
      ! and normalize vectors (not done by rg)
      call nrmvec(ndim,xscr1,eigi)

      if (ierr.ne.0) then
        write(6,*) 'error code from rg: ',ierr
        stop 'logumat (1)'
      end if

      if (ntest.ge.10) write(6,*) ' eigenvalues of matrix:'

*----------------------------------------------------------------------*
*     the eigenvalues are v = exp(a+ib) so the logarithm log(v) yields
*     a and b. as the matrix is unitary, a is always 0 and we are left
*     with b, which is the angle in the complex plane.
*     the angles will be collected in eigr(), later referred to as 
*     matrix D
*----------------------------------------------------------------------*
      ierr = 0
      do ii = 1, ndim
        xmod = eigr(ii)*eigr(ii) + eigi(ii)*eigi(ii)
        if (abs(xmod-1d0).gt.100d0*epsilon(1d0)) ierr = ierr+1
        ang1 = atan2(eigi(ii),eigr(ii))
c        ang2 = acos(eigr(ii))*sign(1d0,eigi(ii))
        if (ntest.ge.10)
     &       write(6,'(i4,2(2x,e20.10),3(2x,f15.10))')
     &       ii,eigr(ii),eigi(ii),xmod,ang1,ang2
        eigr(ii) = ang1
      end do

      if (ierr.gt.0) then
        write(6,*) 'error: detected eigenvalues with |v| != 1'
        stop 'logumat (2)'
      end if

      ! sort components of transformation matrix into 
      ! real and imaginary part: U = A + iB
      !  A on xscr1
      !  B on xscr2

      ii = 0
      do while(ii.lt.ndim)
        ii = ii+1
        if (eigi(ii).eq.0d0) then ! real eigenvalue
          xscr2(1:ndim,ii) = 0d0
        else ! complex pair
          xscr2(1:ndim,ii) = xscr1(1:ndim,ii+1)  ! imag. part
          xscr2(1:ndim,ii+1) = -xscr2(1:ndim,ii) ! and cmplx. conj.
          xscr1(1:ndim,ii+1) = xscr1(1:ndim,ii)
          ii = ii+1                              ! add. increment
        end if
      end do

      if (ntest.ge.100) then
        write(6,*) ' eigenvectors (Re):'
        call wrtmat2(xscr1,ndim,ndim,ndim,ndim)
        write(6,*) ' eigenvectors (Im):'
        call wrtmat2(xscr1,ndim,ndim,ndim,ndim)
      end if

*----------------------------------------------------------------------*
*
*     now we have to calculate U iD U^+
*     
*     the real part is
*                A iD (iB)^+ + (iB) iD A^+ = A D B^T - B D A^T
*
*     the imaginary part is  
*                A iD A^+ + iB iD (iB)^+ = i (A D A^T + B D B^T)
*     
*     as iD has either zero or pairwise conjugate entries, the imaginary
*     part vanishes (note that we started from a real unitary matrix)
*
*----------------------------------------------------------------------*
      
      ! A on xscr1
      ! B on xscr2

      ! A D --> xscr3
      do ii = 1, ndim
        xscr3(1:ndim,ii) = xscr1(1:ndim,ii)*eigr(ii)
      end do

      ! AD B^T --> xlogx
      call matml7(xlogx,xscr3,xscr2,
     &            ndim,ndim,
     &            ndim,ndim,
     &            ndim,ndim,
     &            0d0,1d0, 2 )

      ! B D --> xscr3
      do ii = 1, ndim
        xscr3(1:ndim,ii) = xscr2(1:ndim,ii)*eigr(ii)
      end do

      !-BD A^T --> xlogx
      call matml7(xlogx,xscr3,xscr1,
     &            ndim,ndim,
     &            ndim,ndim,
     &            ndim,ndim,
     &            1d0,-1d0, 2 )

      if (ntest.ge.100) then
        write(6,*) ' result on xlogx:'
        call wrtmat2(xlogx,ndim,ndim,ndim,ndim)
      end if

      return

      end

*----------------------------------------------------------------------*
      subroutine nrmvec(ndim,eigvec,eigvi)
*----------------------------------------------------------------------*
*     normalize the eigenvectors in array eigvec(ndim,ndim)
*     imaginary pairs are handled as described in rg(), eispack.f
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     eigvi(ndim)
      real(8), intent(inout) ::
     &     eigvec(ndim,ndim)

      integer ::
     &     ivec
      real(8) ::
     &     xnrm

      real(8), external ::
     &     inprod

      ivec = 0
      do while (ivec.lt.ndim)
        ivec = ivec+1
        
        xnrm = inprod(eigvec(1,ivec),eigvec(1,ivec),ndim)

        if (eigvi(ivec).ne.0d0) then
          if (ivec+1.gt.ndim) then
            write(6,*) 'inconsistency in eigenvalue structure'
            stop 'nrmvec'
          end if

          xnrm = xnrm + inprod(eigvec(1,ivec+1),eigvec(1,ivec+1),ndim)
        end if

        xnrm = sqrt(xnrm)

        call scalve(eigvec(1,ivec),1d0/xnrm,ndim)

        if (eigvi(ivec).ne.0d0) then
          call scalve(eigvec(1,ivec+1),1d0/xnrm,ndim)
          ivec = ivec+1
        end if

      end do

      return

      end

*----------------------------------------------------------------------*

      integer function ifndmax(ivec,idxoff,lvec,inc)

      implicit none

      integer, intent(in) ::
     &     ivec(*), lvec, inc, idxoff

      integer ::
     &     i, imx, idx

      imx = -huge(imx)
      idx = idxoff
      do i = 1, lvec
        imx = max(imx,ivec(idx))
        idx = idx + inc
      end do
      
      ifndmax = imx

      return
      end

*----------------------------------------------------------------------*

      integer function ifndmin(ivec,idxoff,lvec,inc)

      implicit none

      integer, intent(in) ::
     &     ivec(*), lvec, inc, idxoff

      integer ::
     &     i, imn, idx

      imn = huge(imn)
      idx = idxoff
      do i = 1, lvec
        imn = min(imn,ivec(idx))
        idx = idx + inc
      end do
      
      ifndmin = imn

      return
      end

*----------------------------------------------------------------------*

      subroutine sweepvec(vec,ndim)

* purpose: replace numerical zeroes by real zeroes 
*          (convenient for debugging)

      implicit none

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     vec(ndim)
      
      integer ::
     &     i
      real(8) ::
     &     thr

      thr = 100d0*epsilon(1d0)
      do i = 1, ndim
        if (abs(vec(i)).lt.thr) vec(i) = 0d0
      end do

      return
      end
      SUBROUTINE WRT_2VEC(VEC1,VEC2,NDIM)
*
* Write two vectors 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(NDIM),VEC2(NDIM)
C
      DO I=1,NDIM
        WRITE(6,1010) I,VEC1(I),VEC2(I)
      END DO
 1010 FORMAT(1H0,I5,2X,E18.10,2X,E18.10)
      RETURN
      END
      FUNCTION IELSUM_IND(IACT,NACT,IVEC)
*
* IELSUM_IND = sum_I=1 ^NACT IVEC(IACT(I))
*
      INCLUDE 'implicit.inc'
      INTEGER IVEC(*), IACT(NACT)
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from IELSUM_IND '
        WRITE(6,*) ' NACT = ', NACT
        WRITE(6,*) ' IACT: '
        CALL IWRTMA(IACT,1,NACT,1,NACT)
      END IF
*
      ISUM = 0
      DO I = 1, NACT
       ISUM = ISUM + IVEC(IACT(I))
      END DO
*
      IELSUM_IND = ISUM
*
      RETURN
      END 
      SUBROUTINE MULT_MAT_SPMAT_MAT(AOUT,AIN,X,NAOUT_R,NAOUT_C,NX_C,
     &           IAINPAK )
*
* AOUT = X*AIN, where AIN is sparse
* IF IAINPAK = 1, then AIN is symmetric and delivered in standard 
*                 lower packed form(i*(i-1)/2+j )
*
*. Jeppe Olsen, June 2012 - for transforming Hamilton matrices ...
*
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION X(NAOUT_R,NX_C),AIN(*)
C  AIN(NX_C,NAOUT_C)
*. Output
      DIMENSION AOUT(NAOUT_R,NAOUT_C)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from MULT_MAT_SPMAT_MAT '
       WRITE(6,*) ' ============================='
       WRITE(6,*)
       WRITE(6,'(A,3(2X,I6))') ' NAOUT_R,NAOUT_C,NX_C = ',
     &                           NAOUT_R,NAOUT_C,NX_C 
       WRITE(6,'(A,I3)') ' IAINPAK = ', IAINPAK
      END IF
       
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' X and AIN matrices (input) '
        CALL WRTMAT(X, NAOUT_R, NX_C,NAOUT_R, NX_C)
        WRITE(6,*)
        IF(IAINPAK.EQ.0) THEN
         CALL WRTMAT(AIN, NX_C, NAOUT_C,NX_C, NAOUT_C)
        ELSE
         CALL PRSYM(AIN,NX_C)
        END IF
      END IF
*
      ZERO = 0.0D0
      CALL SETVEC(AOUT,ZERO,NAOUT_R*NAOUT_C)
*
* AOUT(I,J) = Sum(k) X(I,K) AIN(K,J)
*
      DO K = 1, NX_C
       DO J = 1,NAOUT_C
        IF(IAINPAK.EQ.0) THEN
C                 AIN(NX_C,NAOUT_C)
          AINKJ = AIN((J-1)*NX_C + K)
        ELSE
          KJ = MAX(K,J)*(MAX(K,J)-1)/2 + MIN(K,J)
          AINKJ = AIN(KJ)
        END IF
        
        IF(AINKJ.NE.0.0D0) THEN
         DO I = 1, NAOUT_R
          AOUT(I,J) = AOUT(I,J) + AINKJ*X(I,K)
         END DO
        END IF
       END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' The AOUT matrix '
       CALL WRTMAT(AOUT,NAOUT_R, NAOUT_C,NAOUT_R, NAOUT_C)
      END IF
*
      RETURN
      END
      FUNCTION IS_I1_EQ_I2(I1,I2,NDIM)
* Two integer arrays I1 and I2 are given. Are the identical
*
      INCLUDE 'implicit.inc'
      INTEGER I1(NDIM), I2(NDIM)
*
      NTEST = 000
*
      IDENT = 1
      DO I = 1, NDIM
        IF(I1(I).NE.I2(I)) IDENT = 0
      END DO
*
      IS_I1_EQ_I2 = IDENT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  'Output from IS_I1_EQ_I2 '
        IF(IDENT.EQ.1) THEN
          WRITE(6,*) ' The two integer arrays are identical '
        ELSE
          WRITE(6,*) ' The two integer arrays differs '
        END IF
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The two integer arrays '
        CALL IWRTMA3(I1,1,NDIM,1,NDIM)
        WRITE(6,*)
        CALL IWRTMA3(I2,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_XVAL_WITH_THRES(A,THRES, XVAL, NDIM,IVAL)
*
* Find first element in A with ABS(A-XVAL).LE.THRES
*
*. Jeppe Olsen, Feb 13, 2013
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(NDIM)
*
      NTEST = 100
*
      IVAL = 0
      DO I = 1, NDIM
        IF(ABS(A(I)-XVAL).LE.THRES) THEN
          IVAL = I
          GOTO 1001
        END IF
      END DO
 1001 CONTINUE
*
      IF(IVAL.EQ.0) THEN
         WRITE(6,*) ' FIND_XVAL_WITH_THRES in trouble '
         WRITE(6,*) ' Requested value and tolerance ', XVAL, THRES
         WRITE(6,*) ' No such value obtained '
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from FIND_XVAL_WITH_THRES '
        WRITE(6,*) ' target value: ', XVAL
        WRITE(6,*) ' Obtained address ', IVAL
      END IF
*
      RETURN
      END
      
 
      
c $Id$
