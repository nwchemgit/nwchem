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
      SUBROUTINE BLKTRN( N, M, NELEBK, A, MAP1, MAP2,
     $                   ITYPE, ATRANS, WORK, IWORK, INFO )

      INTEGER           N, M, NELEBK, ITYPE, INFO
      INTEGER           MAP1(M), MAP2(N), IWORK(*)
      DOUBLE PRECISION  A(*), ATRANS(*), WORK(*)

*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*
*  Let Z be an N-by-M block matrix, i.e., Z_{ij} is a vector of length
*  NELEBK, distributed to processors by block columns.  This routine
*  does a block transpose of Z, i.e. shifts Z so that it is distributed
*  by block rows.  Thus, (Z transpose)_{ij} = Z_{ji} where
*  (Z transpose) is distributed to processors by block columns (the
*  block rows of Z).  The data within the blocks is not modified by the
*  transpose, thus (Z transpose)_{ij} = Z_{ji} is actually a vector
*  equallity for vectors of length NELEBK.
*
*  On entry A contains the part of Z owned by this processor, and on
*  exit ATRANS contains the part of (Z transpose) owned by this
*  processor.
*
*  The inverse of this routine can is gotten by interchanging:
*  N and M, A and ATRANS, and MAP1 and MAP2.  I.e.,
*
*  CALL BLKTRN( N, M, NELEBK, A, MAP1, MAP2,
*               ITYPE, ATRANS, WORK, IWORK, INFO )
*
*  followed by
*
*  CALL BLKTRN( M, N, NELEBK, ATRANS, MAP2, MAP1,
*               ITYPE, A, WORK, IWORK, INFO )
*
*  block transposes and then untransposes A.
*
*
*  !!! This routine does no input error checking.  !!!
*
* ARGUMENTS
* =========
*  Define:
*
*       ME .................... this processor's id,
*       NAPROC ................ number of allocated processors,
*                               MXNPRC(),
*       LSTCNT( I, MAP, K ) ... be the number of entries in
*                               the array MAP(1:K) equal to I
*       N1MAX ................. The maximum over I=0 to NAPROC-1
*                               of LSTCNT( I, MAP1, M ) ).
*       N2MAX ................. The maximum over I=0 to NAPROC-1
*                               of LSTCNT( I, MAP2, N ) ).
*
*  Then
*
*  N ........ (Input) INTEGER
*             Number of block rows in Z.
*
*  M ........ (Input) INTEGER
*             Number of block columns in Z.
*
*  NELEBK ... (Input) INTEGER
*             Number of (DOUBLE PRECISION) elements in each block
*             of Z (the same for all blocks of Z).
*
*  A ........ (Input) DOUBLE PRECISION array,
*                     dim( NELEBK * N * LSTCNT( ME, MAP1, M ) )
*
*             Let Z be an N-by-M block matrix in which each block
*             contains NELEBK DOUBLE PRECISION elements.  Furthermore,
*             let Z be distributed to processors by block columns
*             according to MAP1 (i.e., processor MAP1(I) owns the
*             Ith block column of Z).  Then A contains the part of Z
*             owned by this processor stored in packed form.  In
*             particular, suppose this processors owns the Ith block
*             column of Z and that this is the IVECth block column
*             of Z owned by this processor, then
*
*                A(I0+(I-1)*NELEBK+1:I0+I*NELEBK) = Z_{I,J},
*
*                where I0 = (IVEC-1)*NELEBK*N.
*
*             This amounts to assuming that that calling routine
*             declared A to be DIMENSION A(NELEBK, N, *)
*             and defined 
*
*                A(1:NELEBK, 1:N, IVEC ) = Z_{I,J}
*
*
*             If, in the calling program A is:
*
*             DOUBLE PRECISION A(N1,N2), then in this routine set
*                                        NELEBK = 1, N = N1, M = N2
*
*             DOUBLE PRECISION A(N1,N2,N3), then in this routine set
*                                           NELEBK = N1, N = N2, M = N3
*
*             DOUBLE COMPLEX A(N1,N2), then in this routine set
*                                      NELEBK = 2, N = N1, M = N2
*
*             DOUBLE COMPLEX A(N1,N2,N3), then in this routine set
*                                      NELEBK = 2*N1, N = N2, M = N3
*
*             If the A coming into this routine is 2 dimensional array,
*             then it represents a block matrix where the block size is
*             one (two when a is DOUBLE COMPLEX).
*
*             If the A coming into this routine is 3 dimensional array,
*             then it repsents a block matrix where the block size is
*             equal to one (two if A is DOUBLE COMPLEX) times the 
*             length of the array in the first dimension.
*
*  MAP1 ..... (Input) INTEGER array, dim( M )
*             List whose Ith entry is the id of the processor
*             which owns the Ith block column of A.
*
*  MAP2 ..... (Input) INTEGER array, dim( N )
*             List whose Ith entry is the id of the processor
*             which is to own the Ith block column of A transpose
*             (i.e., the Ith block row of A) after A is block transposed.
*
*  ITYPE .... (Input) INTEGER
*             Message type to use for interprocessor communitation.
*
*  ATRANS ... (Output) DOUBLE PRECISION array,
*                      dim( NELEBK * M * LSTCNT( ME, MAP2, N ) )
*
*             Let Z be as in the description of A above.  Then
*             this routine computes the block transpose of Z
*             and distributes (Z transpose) to processors to
*             by block columns (the block rows of Z)
*             according to MAP2 (i.e., processor MAP2(I) owns the
*             Ith block column of (Z transpose).
*
*             ATRANS contains the part of (Z transpose)
*             owned by this processor stored in packed form.  In
*             particular, suppose this processors owns the Ith block
*             column of (Z transpose) and that this is the IVECth
*             block column of (Z transpose) owned by this processor,
*             then
*
*             ATRANS(I0+(I-1)*NELEBK+1:I0+I*NELEBK) = (Z trans.)_{I,J}
*                                                   = Z_{J,I}
*
*                where I0 = (IVEC-1)*NELEBK*M
*
*             Z_{I,J} is a DOUBLE PRECISION vector of length NELEBK.
*
*             This amounts to assuming that that calling routine
*             declared ATRANS to be DIMENSION A(NELEBK, M, *)
*             in which case
*
*                ATRANS(1:NELEBK, 1:M, IVEC ) = (Z transpose)_{I,J}
*
*
*             If, in the calling program the input A is:
*
*             DOUBLE PRECISION A(N1,N2), then the returned ATRANS is
*                consistent with the calling program declaring
*                DOUBLE PRECISION ATRANS(N2,N1)
*
*             DOUBLE PRECISION A(N1,N2,N3), then the returned ATRANS is
*                consistent with the calling program declaring
*                DOUBLE PRECISION ATRANS(N1,N3,N2)
*
*             DOUBLE COMPLEX A(N1,N2), then the returned ATRANS is
*                consistent with the calling program declaring
*                DOUBLE COMPLEX ATRANS(N2,N1)
*
*             DOUBLE COMPLEX A(N1,N2,N3), then the returned ATRANS is
*                consistent with the calling program declaring
*                DOUBLE COMPLEX ATRANS(N1,N3,N2)
*
*  WORK ..... (Workspace) DOUBLE PRECISION array,
*                         dim( 10 +
*                              NELEBK * N1MAX * LSTCNT( ME, MAP2, N ) +
*                              NELEBK * N2MAX * LSTCNT( ME, MAP1, M )  )
*
*  IWORK .... (Workspace) INTEGER array, dim( N2MAX )
*
*  INFO ..... (Output) INTEGER
*             = 0 then no errors detected, only value
*                      currently returned.
*
*-----------------------------------------------------------------------

      INTEGER           NAPROC, ME, IOFFST, I

      INTEGER           LSTCNT
      EXTERNAL          LSTCNT

      INTEGER           MXNPRC, MXMYND
      EXTERNAL          MXNPRC, MXMYND

      EXTERNAL          BLKTR2
*-----------------------------------------------------------------------

      NAPROC = MXNPRC()
      ME     = MXMYND()

*     Compute maximum number of elements to be received in a message
*     and set IOFFST to 1 + this value.  IOFFST is then the offset
*     in WORK to use for the write buffer.

      IOFFST = 0
      DO 10 I = 0, NAPROC-1
        IOFFST = MAX( IOFFST, LSTCNT( I, MAP1, M ) )
   10 CONTINUE

      IOFFST = 1 + NELEBK * IOFFST * LSTCNT( ME, MAP2, N ) + 10

      CALL BLKTR2( N, M, NELEBK, A, MAP1, MAP2, ITYPE,
     $             WORK(1), WORK(IOFFST), ATRANS, IWORK, INFO )

      RETURN
      END
      SUBROUTINE BLKTR2( N, M, NELEBK, A, MAP1, MAP2, ITYPE,
     $                   BUFIN, BUFOUT, ATRANS, IWORK, INFO )

      INTEGER           N, M, NELEBK, ITYPE, INFO
      INTEGER           MAP1(M), MAP2(N), IWORK(*)
      DOUBLE PRECISION  A(*), ATRANS(*), BUFIN(*), BUFOUT(*)
*-----------------------------------------------------------------------

*  Let A be an N-by-M block matrix distributed on processors
*  by columns according to MAP1.  This routine transposes A
*  so that (A transpose) is distributed by columns according
*  to MAP2.

c
c     Assume processes are directly connected ... sequence of
c     exchanges is designed for maximum overlap of messages so
c     will complete in time O(nproc*lenbuf) if the hardware permits
c     (and all is done async and messages long enuf to overcome
c     startup overhead ... get real!).
c
c     This may need modification to reflect specific hardware toplogy.
c
c     e.g. for 6 processors have 5 steps each with 3 exchanges in parallel
c
c     0-1     0-2     0-3     0-4     0-5
c     5-2     1-3     2-4     3-5     4-1
c     4-3     5-4     1-5     2-1     3-2
c
*-----------------------------------------------------------------------
#include          "datasize.h"

      INTEGER          NAPROC, ME, NGROUP, IG, IP, NVECS,
     $                 NVECST, NELE, IVAL

      INTEGER          LSTCNT
      EXTERNAL         LSTCNT

      INTEGER          MXNPRC, MXMYND, MXREAD, MXWRIT
      EXTERNAL         MXNPRC, MXMYND, MXREAD, MXWRIT

      EXTERNAL         FILLUP, DRAIN, PAIRUP

*-----------------------------------------------------------------------

      INFO = 0

      NAPROC = MXNPRC()
      ME     = MXMYND()

      NGROUP = NAPROC + MOD( NAPROC, 2 )

      NVECS  = LSTCNT( ME, MAP1, M )
      NVECST = LSTCNT( ME, MAP2, N )

*     Copy any data I have in A to the part of ATRANS that I
*     own.

      CALL FILLUP( NVECS, N, A, NELEBK, ME, MAP2, BUFOUT, NELE, IWORK )

      IF( NELE .GT. 0 )
     $  CALL DRAIN( NVECST, M, ATRANS, NELEBK, ME, MAP1, BUFOUT )

      DO 10 IG = 1, NGROUP - 1

        CALL PAIRUP( NGROUP, ME, IG, IP )

        IF( IP .LT. NAPROC ) THEN

          IF( IP .GT. ME ) THEN

            CALL FILLUP( NVECS, N, A, NELEBK, IP, MAP2, BUFOUT, NELE,
     $                   IWORK )

            IF( NELE .GT. 0 )
     $        IVAL = MXWRIT( BUFOUT, NBYTED * NELE, IP, ITYPE )

            NELE = NELEBK * NVECST * LSTCNT( IP, MAP1, M )

            IF( NELE .GT. 0 ) THEN
              IVAL = MXREAD( BUFIN,  NBYTED * NELE, IP, ITYPE )
              CALL DRAIN( NVECST, M, ATRANS, NELEBK, IP, MAP1, BUFIN )
            ENDIF

           ELSE

             NELE = NELEBK * NVECST * LSTCNT( IP, MAP1, M )

             IF( NELE .GT. 0 ) THEN
               IVAL = MXREAD( BUFIN,  NBYTED * NELE, IP, ITYPE )
               CALL DRAIN( NVECST, M, ATRANS, NELEBK, IP, MAP1, BUFIN )
             ENDIF

             CALL FILLUP( NVECS, N, A, NELEBK, IP, MAP2, BUFOUT, NELE,
     $                    IWORK )

             IF( NELE .GT. 0 )
     $         IVAL = MXWRIT( BUFOUT, NBYTED * NELE, IP, ITYPE )

           ENDIF

         ENDIF
 10   CONTINUE

      RETURN
      END
