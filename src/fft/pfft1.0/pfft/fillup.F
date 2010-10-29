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
      SUBROUTINE FILLUP( NVECS, N, A, NELEBK, IP, MAP, BUF, NELE,
     $                   IWORK )

      INTEGER           NVECS, N, NELEBK, IP, NELE

      INTEGER           MAP(N), IWORK(*)
      DOUBLE PRECISION  A(*), BUF(*)

*-----------------------------------------------------------------------

*  Let A be an N-by-* block matrix distributed on processors
*  by columns.  This routine loads the part of A that I own
*  which needs to be sent to processor IP when transposing A
*  (and having A transpose stored by columns).

*  NVECS .... (Input) INTEGER
*             Number of columns of A that I own.

*  N ........ (Input) INTEGER
*             Number of rows in A.

*  A ........ (Input) DOUBLE PRECISION array, dim(  )
*             An N-by-* block matrix in which each block contains
*             NELEBK elements.  All entries in
*             each block are assumed to be stored consecutively.
*             The NVECS columns of A that I own should be stored
*             in packed format. Thus, A(1+(J-1)*NELEBK:J*NELEBK)
*             contains the Jth column of A that I own, J = 1 to NVECS.

*  NELEBK ... (Input) INTEGER
*             Number of (DOUBLE PREICISION) elements in each block
*             of A.

*  IP ....... (Input) INTEGER
*             The id of the processor to which I need to
*             send his part of A transpose that I currently own.

*  MAP ...... (Input) INTEGER array, dim( N )
*             List whose Ith entry is the id of the processor
*             which owns the Ith column of (A Transpose),
*             which is the Ith row of A.

*  BUF ...... (Output) DOUBLE PRECISION array, dim(    )
*             Buffer containing the blocks of (A transpose) which
*             belong to processor IP, but which I currently
*             own and must send to IP.  The blocks are stored in
*             packed format.  The packing is by columns of A
*             (NOT A transpose).

*  NELE ..... (Output) INTEGER
*             Number of elements of BUF filled with data to
*             send.  Thus, must send processor IP a message of length
*             NELE * (# of bytes in a DOUBLE PRECISION scaler).

*  IWORK .... (Workspace) INTEGER array, dim( LSTCNT( IP, MAP, N ) ).

*-----------------------------------------------------------------------

      INTEGER           NIP, KELE, K, I, J

      EXTERNAL          MAPFIL

*-----------------------------------------------------------------------

      CALL MAPFIL( IP, N, MAP, NIP, IWORK  )

      IF( NELEBK .EQ. 1 ) THEN

        KELE = 1
        DO 20 K = 1, NVECS

          J = N * ( K - 1 )
          DO 10 I = 1, NIP

            BUF( KELE ) = A( J + IWORK( I ) )
            KELE        = KELE + 1

   10     CONTINUE

   20   CONTINUE

      ELSE

        KELE = 1
        DO 120 K = 1, NVECS

          J = NELEBK * N * ( K - 1 ) + 1
          DO 110 I = 1, NIP

            CALL dcopy( NELEBK, A( J + NELEBK * (IWORK( I )-1) ),
     $                  1, BUF( KELE ), 1 )

            KELE = KELE + NELEBK

  110     CONTINUE

  120   CONTINUE

      ENDIF
  
      NELE = KELE - 1

      RETURN
      END
