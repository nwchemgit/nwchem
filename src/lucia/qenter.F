      SUBROUTINE LUCIAQENTER(NAMEIN)
      IMPLICIT REAL*8 (A-H,O-Z)
*
* ROUTINE FOR TRACING PROGRAM FLOW AND ROUTINE TIMING.
*
      PARAMETER (MXROUT=200, MXLEV=20)
      DIMENSION NO(MXLEV)
      DIMENSION TCP(MXROUT),TIO(MXROUT),NCALL(MXROUT)
      CHARACTER NAME(MXROUT)*6,NAMEIN*(*)
      COMMON /LQDATA/ NO, TCP, TCP0, TIO, TIO0, NCALL, NAME, LEVEL,NROUT
*
      CALL TIMING(TCP1,DUM,TIO1,DUM)
*
* UPDATE CALLER
*
      IPRINT = 0
      IF(LEVEL.GT.0) THEN
         N=NO(LEVEL)
         TCP(N)=TCP(N)+TCP1-TCP0
         TIO(N)=TIO(N)+TIO1-TIO0
         IF(IPRINT.GT.0) THEN
            WRITE(*,*)
            WRITE(*,*) '****** ',NAME(N),' CALLING ',NAMEIN,' ******'
            WRITE(*,*)
         END IF
      ELSE IF(IPRINT.GT.0) THEN
         WRITE(*,*)
         WRITE(*,*) '****** MAIN CALLING ',NAMEIN,' ******'
         WRITE(*,*)
      END IF
*
* UPDATE CALLED
*
      LEVEL=LEVEL+1
      IF(LEVEL.GT.MXLEV) THEN
         WRITE(*,*) 'OVERFLOW IN SUBROUTINE STACK',
     *              ' IN ROUTINE QENTRY'
         STOP
      END IF
      N=0
      DO 100 I=1,NROUT
         IF(NAMEIN.EQ.NAME(I)) N=I
100   CONTINUE
      IF(N.EQ.0) THEN
         NROUT=NROUT+1
         IF(NROUT.GT.MXROUT) THEN
            WRITE(*,*) 'OVERFLOW IN SUBROUTINE LIST',
     *                 ' IN ROUTINE QENTRY'
            STOP
         END IF
         N=NROUT
         NAME(N)=NAMEIN
         TCP(N)=0.0D0
         TIO(N)=0.0D0
         NCALL(N)=0
      END IF
      NCALL(N)=NCALL(N)+1
      NO(LEVEL)=N
      CALL TIMING(TCP0,DUM,TIO0,DUM)
      RETURN
      END
      SUBROUTINE LUCIAQEXIT(NAMEIN)
      IMPLICIT REAL*8 (A-H,O-Z)
*
      PARAMETER (MXROUT=200, MXLEV=20)
      DIMENSION NO(MXLEV)
      DIMENSION TCP(MXROUT),TIO(MXROUT),NCALL(MXROUT)
      CHARACTER NAME(MXROUT)*6,NAMEIN*(*)
      COMMON /LQDATA/ NO, TCP, TCP0, TIO, TIO0, NCALL, NAME, LEVEL,NROUT
*
      IPRINT = 0
*
      CALL TIMING(TCP1,DUM,TIO1,DUM)
      IF(LEVEL.EQ.0) THEN
         WRITE(*,*) 'QEXIT CALLED WITH EMPTY SUBROUTINE STACK'
         STOP
      END IF
*
* UPDATE EXITED ROUTINE
*
      N=NO(LEVEL)
      IF(NAME(N).NE.NAMEIN) THEN
         WRITE(*,*) 'ATTEMPT TO EXIT ',NAMEIN,' WHEN IN ',NAME(N)
         STOP
      END IF
      TCP(N)=TCP(N)+TCP1-TCP0
      TIO(N)=TIO(N)+TIO1-TIO0
      IF(IPRINT.GT.0) THEN
         IF(LEVEL.GT.1) THEN
            M=NO(LEVEL-1)
            WRITE(*,*)
            WRITE(*,*) '****** ',NAMEIN,' RETURNING TO ',NAME(M),
     *                 ' ******'
            WRITE(*,*)
         ELSE
            WRITE(*,*)
            WRITE(*,*) '****** ',NAME(N),' RETURNING TO MAIN ******'
            WRITE(*,*)
         END IF
      END IF
*
*
*
      LEVEL=LEVEL-1
      CALL TIMING(TCP0,DUM,TIO0,DUM)
      RETURN
      END
      SUBROUTINE LUCIAQSTAT
      IMPLICIT REAL*8 (A-H,O-Z)
*
      PARAMETER (MXROUT=200, MXLEV=20)
      DIMENSION NO(MXLEV)
      DIMENSION TCP(MXROUT),TIO(MXROUT),NCALL(MXROUT)
      CHARACTER NAME(MXROUT)*6
      COMMON /LQDATA/ NO, TCP, TCP0, TIO, TIO0, NCALL, NAME, LEVEL,NROUT
*
      CALL TIMING(TCP2,DUM,TIO2,DUM)
      TOTTIME = 0.0D0
      WRITE(*,*)
C     WRITE(*,'(5X,A)') 'RUNTIME STATISTICS'
C     WRITE(*,'(5X,A)') '------------------'
C     WRITE(*,*)
      WRITE(*,'(5X,A)') 'ROUTINE     CALLS      CPU     I/O'
      DO 300 I=1,NROUT
         WRITE(*,'(6X,A6,I9,F10.2,F8.2)')
     *      NAME(I),NCALL(I),TCP(I),TIO(I)
         TOTTIME = TOTTIME + TCP(I)
300   CONTINUE
      WRITE(*,*)
      WRITE(6,'(A,E8.3)') ' Total time in above routines: ', TOTTIME
      CALL TIMING(TCP3,DUM,TIO3,DUM)
      TCP0=TCP0+TCP3-TCP2
      TIO0=TIO0+TIO3-TIO2
      RETURN
      END
      BLOCK DATA LQBLOCK
      IMPLICIT REAL*8 (A-H,O-Z)
*
      PARAMETER (MXROUT=200, MXLEV=20)
      DIMENSION NO(MXLEV)
      DIMENSION TCP(MXROUT),TIO(MXROUT),NCALL(MXROUT)
      CHARACTER NAME(MXROUT)*6
      COMMON /LQDATA/ NO, TCP, TCP0, TIO, TIO0, NCALL, NAME, LEVEL,NROUT
*
      DATA LEVEL,NROUT/0,0/
*
      END
      SUBROUTINE SETTIM
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE CPU0,CPU1,WALL0,WALL1
      DATA CPU0,CPU1,WALL0,WALL1/4*0.0D0/
      CPU  = 0.01D0*mClock()
      CPU  = 0.01D0*mClock()
      WALL = CPU
      CPU0 = CPU
      CPU1 = CPU
      WALL0= WALL
      WALL1= WALL
      RETURN
      END
      SUBROUTINE TIMING(CPUA,CPUE,TIOA,TIOE)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE CPU0,CPU1,WALL0,WALL1
      DATA CPU0,CPU1,WALL0,WALL1/4*0.0D0/
c at least for linux:
      CALL ATIM(CPU,WALL)
c
c      CPU  = 0.01D0*LTIMEX()
      WALL = CPU
      CPUA = (CPU-CPU0)
      CPUE = (CPU-CPU1)
      TIOA = (WALL-WALL0-CPUA)
      TIOE = (WALL-WALL1-CPUE)
      CPU1 = CPU
      WALL1= WALL
      RETURN
      END
      FUNCTION LTIMEX()
      IMPLICIT REAL*8(A-H,O-Z)
*
* CPU SEC/100
*
c      LTIMEX = 0.0D0
*. RISC 6K
C     LTIMEX = mCLOCK(*)
*. Cray
C     LTIMEX = 100.0D0*SECOND()
*.Linux environment
C     LTIMEX = 100.0D0*TIMELINUX()
*. Regatta
      CALL CPU_TIME(USTIME)
      LTIMEX = 100.0D0*USTIME
      RETURN
      END
c
      FUNCTION WIMEX()
      IMPLICIT REAL*8(A-H,O-Z)
cc     CALL WALL_TIME(USTIME)
      WIMEX = 100.0D0*USTIME
      RETURN
      END
c
c      SUBROUTINE ATIM(CPU,WALL)
c      IMPLICIT REAL*8(A-H,O-Z)
c      CALL CPU_TIME(CPU)
c      WALL=CPU
cc      CALL WALL_TIME(WALL)
c      RETURN
c      END
c
      SUBROUTINE QTRACE()
      IMPLICIT REAL*8 (A-H,O-Z)
*
      PARAMETER (MXROUT=200, MXLEV=20)
      DIMENSION NO(MXLEV)
      DIMENSION TCP(MXROUT),TIO(MXROUT),NCALL(MXROUT)
      CHARACTER NAME(MXROUT)*6
      COMMON /LQDATA/ NO, TCP, TCP0, TIO, TIO0, NCALL, NAME, LEVEL,NROUT

      WRITE(*,'(2(/,x,a),/,2(/,x,a))')
     &     'INTERNAL STACK TRACE:',
     &     '=====================',
     &     ' level  module name',
     &     '--------------------'
      DO ILEV = 1, LEVEL
        WRITE(*,'(3x,i3,4x,a)') ILEV, NAME(NO(LEVEL))
      END DO
      WRITE(*,'(x,a,/)')
     &     '--------------------'

      END

c $Id$
