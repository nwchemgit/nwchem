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
      SUBROUTINE PFFT1( IDIR, ISIGN, N, X, WORK )

      INTEGER           IDIR, ISIGN, N
      DOUBLE PRECISION  X(*), WORK(*)
*-----------------------------------------------------------------------
*
*     Do N-point, 1-dimensional, complex-to-complex fft.
*
*-----------------------------------------------------------------------

#ifdef FFTPACK
      EXTERNAL          DCFFTB, DCFFTF
#endif

* EXECUTABLE CODE
************************************************************************


#ifdef TESTTRN
      RETURN
#endif

#ifdef FFTPACK
      IF( ISIGN .EQ. 1 ) THEN
        CALL DCFFTB( N, X, WORK )
      ELSE
        CALL DCFFTF( N, X, WORK )
      ENDIF
#else

      WRITE(*,1000)
      STOP

 1000 FORMAT(  ' Subroutine PFFT1 compiled without any',
     $       / ' single-processor fft library defined.'  )

#endif


      RETURN
      END
