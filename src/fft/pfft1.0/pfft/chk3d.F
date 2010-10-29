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
      SUBROUTINE CHK3D( CALLER, IDIR, ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $                   LWORK, WORK, LIWORK, IWORK, INFO )

      INTEGER           IDIR, ISIGN, N1, N2, N3, LWORK, LIWORK, INFO
      INTEGER           MAP1(N3), MAP2(*), IWORK(LIWORK)
      DOUBLE PRECISION  X(*), WORK(LWORK)
	character*6	     CALLER
c	external caller

*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*
*  This routine checks the input data to the parallel 3D fft routines
*  in PFFT.  CALLER is the name of the PFFT routine which called this
*  routine.  CALLER must be ALL UPPER CASE or ALL LOWER CASE letters.
*
*  Only to be called internally from PFFT routines.
*
*-----------------------------------------------------------------------
#include "datasize.h"

      INTEGER           NAPROC, ME, LINFO, IPROC, NPROC, IIWORK, NDIFF,
     $                  M2
      LOGICAL           RCFFT

      INTEGER           LSTCNT
      EXTERNAL          LSTCNT

      INTEGER           MXNPRC, MXMYND, MXLBUF, MAPCHK, IGDIFF
      EXTERNAL          MXNPRC, MXMYND, MXLBUF, MAPCHK, IGDIFF
      character *26     msg  

*-----------------------------------------------------------------------

      NAPROC = MXNPRC()
      ME     = MXMYND()

      INFO = 0

      RCFFT =  CALLER .EQ. 'PRFFT3'  .OR.  CALLER .EQ. 'prfft3'
c
*     Set M2 = number of elements in MAP2.
      
      IF( N1 .EQ. 1  .AND.  RCFFT ) THEN
         M2 = N2 / 2 + 1
      ELSE
         M2 = N2
      ENDIF
      
      IF( N2 .LT. 1 ) THEN
         INFO = -3
         
      ELSE IF( N3 .LT. 1 ) THEN
         INFO = -4
         
      ELSE IF( MAPCHK( MAP1, N3 ) .NE. 0 ) THEN
         INFO = -6
        
      ELSE IF( MAPCHK( MAP2, M2 ) .NE. 0 ) THEN
        INFO = -7

      ELSE IF( LWORK .LE. MXLBUF() / NBYTED + 1 ) THEN
*       Not enough workspace for GEXIT.
        INFO = -8

      ELSE IF( LIWORK .LE. 2 * NAPROC ) THEN
*       Not enough workspace to call PRLST2.
        INFO = -10

      ELSE IF( LSTCNT( ME, MAP1, N3 )  .EQ. 0  .AND.
     $         LSTCNT( ME, MAP2, M2 )  .EQ. 0        ) THEN
        INFO = -50

      ENDIF

      IF( INFO .NE. 0  .AND. RCFFT ) 
     $  INFO = INFO - 1 

      IF( INFO .NE. 0 ) THEN
        WRITE(*,*) ' ME = ', ME,
     $             ' Illegal input to ', CALLER, '. INFO= ', INFO
        CALL XSTOP( CALLER )
      ENDIF

*-----------------------------------------------------------------------

*     Determine set of unique processor ids in MAP1 and/or MAP2,
*     NPROC, and store them in a processor list, IWORK(1:NPROC).

      IPROC = 1

      CALL PRLST2( N3, MAP1, M2, MAP2, NPROC, IWORK( IPROC ),
     $             IWORK( NAPROC+1 ) )

      IIWORK = NPROC + 1

*-----------------------------------------------------------------------

*     Finish checking input data

      IF( ISIGN .NE. 1  .AND.  ISIGN .NE. -1 ) THEN
        INFO = -1

      ELSE IF( N1 .LT. 1 ) THEN
        INFO = -2

      ELSE IF( N2 .LT. 2 ) THEN
        INFO = -3

      ELSE IF( N3 .LT. 2 ) THEN
        INFO = -4

      ELSE IF( INFO .EQ. 0  .AND.
     $         LWORK .LT. NPROC + MAX( 10, M2, N3 ) ) THEN
*       Not enough room to do global compare of input data
        INFO = -8

      ENDIF

      IF( RCFFT ) THEN

        IF( INFO .NE. 0 ) THEN
          INFO = INFO - 1 
        ELSE IF( IDIR .NE. 1  .AND.  IDIR .NE. -1 ) THEN
          INFO = -1
        ENDIF

      ENDIF

      IF( INFO .NE. 0 )
     $  WRITE(*,*) ' ME = ', ME,
     $             ' Illegal input to ', CALLER, '. INFO= ', INFO

      LINFO = INFO
      CALL GEXIT( LINFO, NPROC, IWORK, WORK )

      IF( LINFO .NE. 0 ) THEN
        IF( INFO .EQ. 0 ) INFO = -51
        RETURN
      ENDIF

*-----------------------------------------------------------------------

*     Make sure IDIR, ISIGN, N1, N2, N3, MAP1, and MAP2 are the same on
*     all processors.

      IWORK( IIWORK   ) = ISIGN
      IWORK( IIWORK+1 ) = N1
      IWORK( IIWORK+2 ) = N2
      IWORK( IIWORK+3 ) = N3
      IWORK( IIWORK+4 ) = IDIR
c
c
      ndiff = 0
c      
      msg = 'IDIR, ISIGN, N1, N2, or N3'
      NDIFF = IGDIFF( 5,
     $     IWORK( IIWORK ), NPROC,
     $     IWORK( IPROC ), CALLER, 
     $     msg,
     $     IWORK( IIWORK + 5 )  )
c     
c     

      NDIFF = MAX( NDIFF, 
     $     IGDIFF( n3, MAP1, NPROC,
     $     IWORK( IPROC ), CALLER, 'MAP1',
     $     IWORK( IIWORK ) ) )
c     
c     
      NDIFF = MAX( NDIFF, 
     $     IGDIFF( m2, MAP2, NPROC,
     $     IWORK( IPROC ), CALLER, 'MAP2',
     $     IWORK( IIWORK ) ) )
c     
      CALL GEXIT( NDIFF, NPROC, IWORK, WORK )
      
      IF( NDIFF .NE. 0 ) THEN
         INFO = -52
        RETURN
      ENDIF

*-----------------------------------------------------------------------

      RETURN
      END
