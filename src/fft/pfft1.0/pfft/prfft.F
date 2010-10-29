*
* $Id$
*
*======================================================================
*
* DISCLAIMER
*
* This material was prepared as an account of work sponsored by an
* agency of the United States Government.  Neither the United States
* Government nor the United States Department of Energy, nor Battelle,
* nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
* COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
* SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
* INFRINGE PRIVATELY OWNED RIGHTS.
*
* ACKNOWLEDGMENT
*
* This software and its documentation were produced with Government
* support under Contract Number DE-AC06-76RLO-1830 awarded by the United
* States Department of Energy.  The Government retains a paid-up
* non-exclusive, irrevocable worldwide license to reproduce, prepare
* derivative works, perform publicly and display publicly by or for the
* Government, including the right to distribute to other Government
* contractors.
*
*======================================================================
*
*  -- PFFT routine (version 1.0) --
*     Pacific Northwest Laboratory
*     April 5, 1995
*
*======================================================================
      SUBROUTINE PRFFT( IDIR, ISIGN, N, X, WORK )

      INTEGER           IDIR, ISIGN, N
      DOUBLE PRECISION  X(*), WORK(*)
*-----------------------------------------------------------------------
*
*     Do N-point, 1-dimensional, real-to-complex or complex-to-real fft.
*
*-----------------------------------------------------------------------

#ifdef FFTPACK
#include "const.h"

      INTEGER         I

      EXTERNAL          DRFFTB, DRFFTF
#endif

* EXECUTABLE CODE
************************************************************************


#ifdef TESTTRN
      RETURN
#endif

#ifdef FFTPACK

      IF( IDIR .EQ. 1 ) THEN

*        Real-to-complex fft

         CALL DRFFTF( N, X, WORK )

         DO 10 I = N, 2, -1
           X(I+1) = X(I)
   10    CONTINUE

         X(2) = ZERO

         IF( MOD( N, 2 ) .EQ. 0 )
     $      X(N+2) = ZERO

         IF( ISIGN .EQ. 1 ) THEN
           DO 20 I = 4, N+1, 2
             X(I) = -X(I)
   20      CONTINUE
         ENDIF

      ELSE

*        complex-to-real fft

         IF( ISIGN .EQ. -1 ) THEN
           DO 30 I = 4, N+1, 2
             X(I) = -X(I)
   30      CONTINUE
         ENDIF

         DO 40 I = 2, N
           X(I) = X(I+1)
   40    CONTINUE

         X(N+1) = ZERO

         IF( MOD( N, 2 ) .EQ. 0 )
     $      X(N+2) = ZERO

         CALL DRFFTB( N, X, WORK )

      ENDIF

#else

      WRITE(*,1000)
      STOP

 1000 FORMAT(  ' Subroutine PRFFT compiled without any',
     $       / ' single-processor fft library defined.'  )

#endif


      RETURN
      END
