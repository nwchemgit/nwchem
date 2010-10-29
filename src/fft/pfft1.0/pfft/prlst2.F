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
      SUBROUTINE PRLST2( N1, MAP1, N2, MAP2, NPROC, PRCLST, IWORK )
*
      INTEGER           N1, N2, NPROC
      INTEGER           MAP1(N1), MAP2(N2), PRCLST(*), IWORK(0:*)
*
*  Return a list of the unique processor ids appearing in MAP1 and/or
*  MAP2.
*
*  N1 ....... (Input) INTEGER
*             Number of processor ids in MAP1
*
*  MAP1 ..... (Input) INTEGER array, dim( N2 )
*             A list of N1 processor ids, possibly with duplication.
*             MUST have 0 <= MAP1(I) < MXNPRC(), I = 1 to N1
*
*  N2 ....... (Input) INTEGER
*             Number of processor ids in MAP2
*
*  MAP2 ..... (Input) INTEGER array, dim( N2 )
*             A list of N2 processor ids, possibly with duplication.
*             MUST have 0 <= MAP2(I) < MXNPRC(), I = 1 to N2
*
*  NPROC .... (Output) INTEGER
*             The number of unique processor ids appearing in
*             MAP1 and/or MAP2.
*
*  PRCLST ... (Output) INTEGER array, dim( NPROC )
*             The list of unique processor ids appearing in
*             MAP1 and/or MAP2.  This list contains no duplicates.
*
*  IWORK .... (Workspace) INTEGER, dim( MXNPRC() )
*             A Work array of length at least MXNPRC().
*-----------------------------------------------------------------------
      INTEGER          NAPROC, I

      INTEGER          MXNPRC
      EXTERNAL         MXNPRC
*
*-----------------------------
*  
      NAPROC = MXNPRC()
*
      NPROC = 0
*
      IF( NAPROC .LE. 0 )
     $  RETURN
*
      DO 10 I = 0, NAPROC - 1
        IWORK( I ) = -1
   10 CONTINUE
*
      DO 20 I = 1, N1
        IF( IWORK( MAP1( I ) ) .EQ.  -1 ) THEN
          IWORK( MAP1( I ) ) = 0
          NPROC              = NPROC + 1
          PRCLST( NPROC )    = MAP1( I )
        ENDIF
   20 CONTINUE
*
      DO 30 I = 1, N2
        IF( IWORK( MAP2( I ) ) .EQ.  -1 ) THEN
          IWORK( MAP2( I ) ) = 0
          NPROC              = NPROC + 1
          PRCLST( NPROC )    = MAP2( I )
        ENDIF
   30 CONTINUE
*
      RETURN
      END
