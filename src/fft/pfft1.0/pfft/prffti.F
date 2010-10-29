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
      SUBROUTINE PRFFTI( IDIR, ISIGN, N, WORK )

      INTEGER           IDIR, ISIGN, N
      DOUBLE PRECISION  WORK(*)
*-----------------------------------------------------------------------
*
*     Initialize N-point 1-dimensional real fft WORK array.
*
*-----------------------------------------------------------------------

#ifdef FFTPACK
      EXTERNAL          DRFFTI
#endif

* EXECUTABLE CODE
************************************************************************


#ifdef FFTPACK
      CALL DRFFTI( N, WORK )
#else

      WRITE(*,1000)
      STOP

 1000 FORMAT(  ' Subroutine PRFFTI compiled without any',
     $       / ' single-processor fft library defined.'  )
#endif

      RETURN
      END
