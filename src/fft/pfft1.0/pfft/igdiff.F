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
      INTEGER FUNCTION IGDIFF( N, DATA, NPROC, PRCLST, ROUTINE,
     $                         MESSAG, WORK )

      INTEGER            N, NPROC
      INTEGER            PRCLST( NPROC )
      CHARACTER*26      MESSAG
      character*(*)	routine
      integer        DATA( N )	, WORK( N )
      
*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*
*   This routine compares the data in DATA on the calling processor
*   with the data in DATA on its neighboring processor in
*   PRCLST(1:NPROC).
*
*   A call to this routine is normally followed by a call to GEXIT.
*   GEXIT forces all processors in PRCLST to stop execution if
*   the flag, usually the value returned by IGDIFF, given to GEXIT
*   is not zero on all processor in PRCLST.
*
*   WARNING:  PRCLST must be the same on all processors calling this
*             routine or deadlock is likely to occur.  
*
*             Any processor which calls this routine but is not in
*             PRCLST prints an error message and STOPs execution without
*             telling any of the other processors what it is doing.
*
* ARGUMENTS
* =========
*    
*   N ........ (Input) INTEGER
*              Number of bytes in character array DATA, see also DATA.
*
*   DATA ..... (Input) CHARACTER*1 array, dim( N )
*              Array of data which should be the same on
*              all processors in PRCLST.  Generally DATA will
*              be an INTEGER or DOUBLE PRECISION array in the
*              calling program and N is equal to the number of BYTES
*              in this data type times the number of elements of
*              the arrray which are to be compared.
*
*   NPROC .... (Input) INTEGER
*              Number of processor ids in PRCLST
*
*   PRCLST ... (Input) INTEGER array, dim ( NPROC )
*              List of processor ids on which DATA is to be
*              compared.  PRCLST is assumed to have no
*              duplications.  Also, PRCLST must be identical
*              on all processors calling this routine.
*
*   ROUTINE .. (Input) CHARACTER*(*)
*              Name of routine to print with MESSAG.  See MESSAG below.
*
*   MESSAG ... (Input) CHARACTER*(*)
*              A message to be printed by any processors who
*              determine that their DATA differs from their
*              neighbors DATA.
*
*   WORK ..... (Workspace) CHARACTER*1 array, dim( N )
*              Workspace 
*
*   IGDIFF ... (Output) INTEGER
*              = 0 if DATA is the same on this processor and its
*                     neighbor in PRCLST
*              = > if DATA on this processor and on this processor's 
*-----------------------------------------------------------------------

      INTEGER            ME, MYINDX, LASTPR, NEXTPR, IVAL, NDIFF, I

      INTEGER            LSTIND
      EXTERNAL           LSTIND

      INTEGER            MXMYND, MXREAD, MXWRIT
      EXTERNAL           MXMYND, MXREAD, MXWRIT

      EXTERNAL           MXPEND
c
#include "datasize.h"
* EXECUTABLE CODE
*-----------------------------------------------------------------------
      
      IGDIFF = 0
      ME = MXMYND()
      if ( nproc .eq. 0 ) then
         WRITE(*,*) ' ME = ', ME, ' Error in IGDIFF nproc = ', nproc
      endif
	
      MYINDX = LSTIND( ME, PRCLST, NPROC )
      
      IF( MYINDX .EQ. 0 ) THEN
      WRITE(*,*) ' ME = ', ME, ' Error in IGDIFF.'
      WRITE(*,*) ' ME = ', ME, ' I called IGDIFF but am not in PRCLST'
      WRITE(*,*) ' ME = ', ME, ' Exiting via IGDIFF. nproc = ', nproc
      CALL MXPEND()
      STOP
      ENDIF	
      
      IF( NPROC .EQ. 1 )
     $     RETURN
      
      MYINDX = MYINDX - 1
      
      LASTPR = MOD( MYINDX + NPROC - 1, NPROC ) + 1
      LASTPR = PRCLST( LASTPR )
      
      NEXTPR = MOD( MYINDX + 1, NPROC ) + 1
      NEXTPR = PRCLST( NEXTPR )
      
      nn = n*nbytei
      IVAL = MXWRIT( DATA, nn, NEXTPR, 1 )
      IVAL = MXREAD( WORK, nn, LASTPR, 1 )

c      IF ( MOD( MYINDX, 2 ) .EQ. 0 ) THEN
c         IVAL = MXWRIT( DATA, nn, NEXTPR, 1 )
c         IVAL = MXREAD(	WORK, nn, LASTPR, 2 )
c      ELSE	         
c         IVAL = MXREAD( WORK, nn, LASTPR, 1 )
c         IVAL = MXWRIT( DATA, nn, NEXTPR, 2 )
c      ENDIF
      
      NDIFF = 0
      DO 10 I = 1, N
         IF( WORK(I) .NE. DATA(I) )
     $        NDIFF = NDIFF + 1
 10   CONTINUE
      
      IF( NDIFF .NE. 0 ) then
         write(*,*) me,  ' error  in igdiff '
c     $  WRITE(*,1000) ME, ROUTINE, MESSAG, ME, LASTPR
      endif
      IGDIFF = NDIFF
      RETURN
 1000 FORMAT(
     $     // 2X, 'ME = ', I4, 2X, 'In ', A, 1X, A
     $     / 2X, 'ME = ', I4, 2X, 'is not the same as on processor ', I4
     $     //)
      END
      
