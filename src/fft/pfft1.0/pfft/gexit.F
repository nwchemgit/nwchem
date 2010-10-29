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
      SUBROUTINE GEXIT( IFLAG, NPROC, PRCLST, WORK )

      INTEGER            IFLAG, NPROC
      INTEGER            PRCLST( NPROC )
      INTEGER            WORK(*)

*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*
*   This routine checks IFLAG in all processors listed in
*   PRCLST(1:NPROC).  If IFLAG is non-zero on any processor in
*   PRCLST(1:NPROC), then all processors in PRCLST(1:NPROC)
*   stop execution and print a note that they are exiting via GEXIT.
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
*   IFLAG .... (Input/Output) INTEGER
*              On Entry:
*               Flag such that if IFLAG is non-zero on any processor
*               in PRCLST(1:NPROC), then all processors in
*               PRCLST(1:NPROC) are to stop execution.
*
*              On Exit:
*               = 0  if input IFLAG = 0 on processors PRCLST(1:NPROC)
*               = -1 if input IFLAG <> 0 on an processor in
*                     PRCLST(1:NPROC).  In this case this routine
*                     calls XSTOP to stop program execution.
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
*   WORK ..... (Workspace) CHARACTER*1 array, dim( MXLBUF() )
*              Workspace  for use by MXCOMBV1.
*
*              Should be DOUBLE PRECISION aligned for best performance
*              in MXCOMBV1.
*
*-----------------------------------------------------------------------

#include           "datasize.h"

      INTEGER            ME, JFLAG(1), IVAL

      INTEGER            LSTCNT
      EXTERNAL           LSTCNT

      INTEGER            SUMIV, MXMYND, MXCOMBV1
      EXTERNAL           SUMIV, MXMYND, MXCOMBV1
  
      EXTERNAL           MXPEND

* EXECUTABLE CODE
*-----------------------------------------------------------------------

      ME = MXMYND()

      IF( LSTCNT( ME, PRCLST, NPROC ) .EQ. 0 ) THEN
    
        WRITE(*,*) ' ME = ', ME, ' Error in GEXIT.'
        WRITE(*,*) ' ME = ', ME, ' I called GEXIT but am not in PRCLST'
        WRITE(*,*) ' ME = ', ME, ' Exiting via GEXIT'

        CALL MXPEND()
        STOP

      ENDIF

      JFLAG(1) = ABS( IFLAG )
      
      IVAL = MXCOMBV1( JFLAG, SUMIV, NBYTEI, 1,
     $     NPROC, PRCLST, 2, WORK)
      
      IF( JFLAG(1) .NE. 0 ) THEN

        CALL XSTOP( 'GEXIT' )
        IFLAG = -1

      ENDIF

      RETURN
      END
