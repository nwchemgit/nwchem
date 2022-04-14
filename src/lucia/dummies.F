      SUBROUTINE RECPRT
      RETURN
      END
      SUBROUTINE XUFLOW
      RETURN
      END
      SUBROUTINE TRACEBACK
      RETURN
      END
      SUBROUTINE  CICNCH
      WRITE(6,*) ' Entering dummy CICNCH is fatal for me'
      STOP ' Dummmy CICNCH'
      END
      SUBROUTINE  SYMCM2
      WRITE(6,*) ' Entering dummy SYMCM2 is fatal for me'
      STOP ' Dummmy SYMCM2'
      END
C     SUBROUTINE  DGEMUL
C     WRITE(6,*) ' Entering dummy DGEMUL is fatal for me'
C     STOP ' Dummmy DGEMUL'
C     END
C     SUBROUTINE  QSTAT 
C     RETURN             
C     END
C     SUBROUTINE  LUCIAQEXIT 
C     RETURN             
C     END
C     SUBROUTINE  LUCIAQENTER
C     RETURN             
C     END
      SUBROUTINE  SQFILE
      WRITE(6,*) ' Entering DUMMY SQFILE is FATAL for me '
      STOP'SQFILE'
      RETURN             
      END
      SUBROUTINE DANAME
C?    WRITE(6,*) ' Entering DUMMY DANAME, nonfatal '
      RETURN
      END
      SUBROUTINE DAFILE
      WRITE(6,*) ' Entering DUMMY DAFILE is FATAL to me'
      STOP'DUMMY DAFILE'
      END
      SUBROUTINE DGEMUL
      WRITE(6,*) ' Entering DUMMY DGEMUL is FATAL to me'
      STOP'DUMMY DGEMUL'     
      END
c      SUBROUTINE DGEMM 
c      WRITE(6,*) ' Entering DUMMY DGEMM  is FATAL to me'
c      STOP'DUMMY DGEMM '     
c      END
      SUBROUTINE MCLOCK
      RETURN
      END
      SUBROUTINE GFLUSH(LU)
*
*. Outer routine for flushing UNITS !
*
*. RISC/6000 
C       CALL UTIL_FLUSH_(LU)
*. CRAY and gfortran
        CALL UTIL_FLUSH(LU)
*
      RETURN
      END
      SUBROUTINE GSLIST
      WRITE(6,*) ' Entering DUMMY GSLIST is FATAL to me'
      STOP'DUMMY GSLIST'     
      END
     

      
c $Id$
