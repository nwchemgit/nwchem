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
      SUBROUTINE WRKSPC( JOB, N, M, NELEBK, MAP1, MAP2, LIWORK, LWORK,
     $                   INFO )

      INTEGER           JOB, N, M, NELEBK, LIWORK, LWORK, INFO
      INTEGER           MAP1(M), MAP2(N)

*-----------------------------------------------------------------------

* Compute workspace requirements of various fft routines:

*   JOB = 1, Workspace for BLKTRN


*-----------------------------------------------------------------------

      INTEGER           NAPROC, ME, N1MAX, N2MAX, I

      INTEGER           LSTCNT
      EXTERNAL          LSTCNT

      INTEGER           MXNPRC, MXMYND
      EXTERNAL          MXNPRC, MXMYND

*-----------------------------------------------------------------------

      NAPROC = MXNPRC()
      ME     = MXMYND()

      INFO = 0

      IF( JOB .LT. 1  .OR.  JOB .GT. 1 ) THEN
        INFO = -2

      ELSE IF( N .LT. 1 ) THEN
        INFO = -2

      ELSE IF( M .LT. 1 ) THEN
        INFO = -3

      ELSE IF( NELEBK .LT. 1 ) THEN
        INFO = -4

      ELSE IF( MAPCHK( MAP1, M ) .NE. 0 ) THEN
        INFO = -5
        
      ELSE IF( MAPCHK( MAP2, N ) .NE. 0 ) THEN
        INFO = -6

      ENDIF

      IF( INFO .NE. 0 ) THEN
        CALL XERPRT( 'MEMUSE', 'Illegal input argument.', INFO )
        RETURN
      ENDIF

      N1MAX = 0
      DO 10 I = 0, NAPROC-1
        N1MAX = MAX( N1MAX, LSTCNT( I, MAP1, M ) )
   10 CONTINUE

      N2MAX = 0
      DO 20 I = 0, NAPROC-1
        N2MAX = MAX( N2MAX, LSTCNT( I, MAP2, N ) )
   20 CONTINUE


*     BLKTRN
*     ======

      LIWORK = MAX( N1MAX, N2MAX )

      LWORK = NELEBK * N1MAX * LSTCNT( ME, MAP2, N ) + 10 +
     $        NELEBK * N2MAX * LSTCNT( ME, MAP1, M ) + 10

      IF( JOB .EQ. 1 )
     $  RETURN


      RETURN
      END
