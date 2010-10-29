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
      SUBROUTINE MFFTC1( FFT, FFTI, IDIR, ISIGN, NELEVC, N, NVEC, X,
     $                   WORK )

      INTEGER           IDIR, ISIGN, NELEVC, N, NVEC
      DOUBLE PRECISION  X(*), WORK(*)

      EXTERNAL          FFT, FFTI
*-----------------------------------------------------------------------
*
*     Do 1-dimensional, N-point, fft on NVEC vectors
*     which are stored consecutively in X, each vector is assumed
*     to occupy NELEVC DOUBLE PRECISION elements of X.
*
*     Generally, 
*       NELEVC = N,         FFT does 1d real-to-real fft
*                2*N,       FFT does 1d complex-to-complex fft
*                2*(N/2+1), FFT does 1d real-to-complex or
*                                    1d complex-to-real fft and all
*                               complex data is partially packed
*-----------------------------------------------------------------------

      INTEGER           IDWORK, KIN, IVEC, I

* EXECUTABLE CODE
************************************************************************

      IDWORK = NELEVC + 1

      CALL FFTI( IDIR, ISIGN, N, WORK( IDWORK ) )

      KIN = -NELEVC

      DO 500 IVEC = 1, NVEC

*       Load next vector, NELEVC elements, on which to do 1D FFT

        KIN  = KIN + NELEVC

        DO 50 I = 1, NELEVC
          WORK( I )   = X( KIN + I )
   50   CONTINUE

        CALL FFT( IDIR, ISIGN, N, WORK, WORK( IDWORK ) )

        DO 60 I = 1, NELEVC
          X( KIN + I ) = WORK( I )
   60   CONTINUE

  500 CONTINUE

      RETURN
      END
