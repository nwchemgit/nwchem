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
      SUBROUTINE DRAIN( NVECS, M, ATRANS, NELEBK, IP, MAP, BUF )

      INTEGER           NVECS, M, NELEBK, IP

      INTEGER           MAP(M)
      DOUBLE PRECISION  ATRANS(*), BUF(*)

*-----------------------------------------------------------------------


*  Let A be an *-by-M block matrix distributed on processors
*  by columns.  This routine unloads the part of A transpose, ATRANS,
*  that I own from BUF, which was received from processor IP when
*  transposing A (and having A transpose stored by columns).

*  NVECS .... (Input) INTEGER
*             Number of columns of A transpose that I own.

*  M ........ (Input) INTEGER
*             Number of columns in A.

*  ATRANS ... (Output) DOUBLE PRECISION array, dim(  )
*             A *-by-M block matrix in which each block contains
*             NELEBK elements.  All entries in
*             each block are assumed to be stored consecutively.
*             The NVECS columns of ATRANS that I own will be stored
*             in packed format. Thus, ATRANS(1+(J-1)*NELEBK:J*NELEBK)
*             contains the Jth column of ATRANS that I own,
*             J = 1 to NVECS.

*  NELEBK ... (Input) INTEGER
*             Number of (DOUBLE PREICISION) elements in each block
*             of A.

*  IP ....... (Input) INTEGER
*             The id of the processor from which I received that part
*             of A transpose that is currently in BUF.

*  MAP ...... (Input) INTEGER array, dim( M )
*             List whose Ith entry is the id of the processor
*             which owns the Ith column of A ,
*             which is the Ith row of A transpose.

*  BUF ...... (Input) DOUBLE PRECISION array, dim(    )
*             Buffer containing the blocks of (A transpose) which
*             belong to me and were received from processor IP.
*             The blocks should be stored in packed format.  The
*             packing is by columns of A (NOT A transpose).

*-----------------------------------------------------------------------


      INTEGER           KELE, K, I, J

*-----------------------------------------------------------------------

      IF( NELEBK .EQ. 1 ) THEN

        KELE = 1
        DO 20 K = 1, M

          IF( MAP( K ) .EQ. IP ) THEN

            J = K
            DO 10 I = 1, NVECS

              ATRANS( J ) =  BUF( KELE )

              KELE = KELE + 1
              J    = J + M

   10       CONTINUE

          ENDIF

   20   CONTINUE

      ELSE

        KELE = 1
        DO 120 K = 1, M

          IF( MAP( K ) .EQ. IP ) THEN

            J = 1 + (K - 1) * NELEBK
            DO 110 I = 1, NVECS

              DO 105 LL = J, J+NELEBK-1
                ATRANS(LL) = BUF( KELE )
                KELE = KELE + 1
  105         CONTINUE

              J    = J + M * NELEBK

  110       CONTINUE

          ENDIF

  120   CONTINUE

      ENDIF

      RETURN
      END
