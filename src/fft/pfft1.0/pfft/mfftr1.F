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
      SUBROUTINE MFFTR1( FFT, FFTI, ICMPLX, IDIR, ISIGN, N1, N2, NSETS,
     $                   X, WORK )

      INTEGER           ICMPLX, IDIR, ISIGN, N1, N2, NSETS
      DOUBLE PRECISION  X(*), WORK(*)

      EXTERNAL          FFT, FFTI
*-----------------------------------------------------------------------
*
*     Multiple, FFTs on Rows, 1 dimensional
*
*     Do 1-dimensional real-to-real (sin or cos ) or complex-to-complex
*     fft on multiple vectors which are stored by
*     rows, second dimension, in a 3-dimensional array.
*
*     This routine will generally NOT WORK right if routine FFT does a
*     real-to-complex or complex-to-real fft.
*
*-----------------------------------------------------------------------

* EXECUTABLE CODE
************************************************************************

      IF( ICMPLX .EQ. 1 ) THEN

*       Real data

        CALL MRFFTR( FFT, FFTI, IDIR, ISIGN, N1, N2, NSETS, X, WORK )

      ELSE
      
*       Complex data

        CALL MCFFTR( FFT, FFTI, IDIR, ISIGN, N1, N2, NSETS, X, WORK )

      ENDIF

      RETURN
      END
