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
      SUBROUTINE MCFFTR( FFT, FFTI, IDIR, ISIGN, N1, N2, NSETS, X,
     $                   WORK )

      INTEGER           IDIR, ISIGN, N1, N2, NSETS
      DOUBLE PRECISION  X(*), WORK(*)

      EXTERNAL          FFT, FFTI
*-----------------------------------------------------------------------
*
*     Do 1-dimensional fft on multiple complex vectors which are stored
*     by rows, second dimension, in a 3-dimensional array.
*
*-----------------------------------------------------------------------
      INTEGER           IDWORK, I1, I2, I, JIN, JOUT, KIN, KOUT
      LOGICAL           NOTFST

* EXECUTABLE CODE
************************************************************************

      IDWORK = 2 * N2 + 1

      CALL FFTI( IDIR, ISIGN, N2, WORK( IDWORK ) )

      NOTFST = .FALSE.

      JOUT = 1                 ! satisfy compiler warnings

      DO 500 I2 = 1, NSETS

        JIN  = 1 + (I2 - 1) * 2 * N1 * N2

        DO 400 I1 = 1, N1

*         Load next complex vector, 2*N2 elements, on which to do 1D FFT

          IF( NOTFST ) THEN

            KIN  = JIN
            KOUT = JOUT
            DO 50 I = 1, 2 * N2, 2

              X( KOUT )   = WORK( I )
              X( KOUT+1 ) = WORK( I+1 )

              WORK( I )   = X( KIN )
              WORK( I+1 ) = X( KIN+1 )

              KIN  = KIN + N1 * 2
              KOUT = KOUT + N1 * 2

   50       CONTINUE

          ELSE

            NOTFST = .TRUE.

            KIN = JIN
            DO 60 I = 1, 2 * N2, 2

              WORK( I )   = X( KIN )
              WORK( I+1 ) = X( KIN + 1)

              KIN = KIN + N1 * 2

   60       CONTINUE

          ENDIF

          JOUT = JIN
          JIN  = JIN + 2

          CALL FFT( IDIR, ISIGN, N2, WORK, WORK( IDWORK ) )

  400   CONTINUE
  500 CONTINUE

      KOUT = JOUT
      DO 600 I = 1, 2 * N2, 2

        X( KOUT )   = WORK( I )
        X( KOUT+1 ) = WORK( I+1 )

        KOUT = KOUT + N1 * 2
  600 CONTINUE

      RETURN
      END
