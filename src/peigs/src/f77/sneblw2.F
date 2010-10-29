*
* $Id$
*
#include "blas_lapackf.h"
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
*  -- PEIGS  routine (version 2.1) --
*     Pacific Northwest Laboratory
*     July 28, 1995
*
*======================================================================
      INTEGER FUNCTION SNEBLW2(N, V, D, E, WORK, INFO )
*     
*     
*     Same functionality as neblw1 from the Eispack library
*     but using LAPACK code.  
*     
*     WARNING:  in this code E(1) is the first subdiagonal
*     element of T and E(N) is junk.  While in
*     neblw1 E(1) is junk and E(2) is the first
*     subdiagonal element of T.
*     
*     D.M. Elwood 09/26/94
*     
*     .. Scalar Arguments ..
      INTEGER            N, INFO
      REAL   V
*     ..
*     .. Array Arguments ..
      REAL   D( * ), E( * ), WORK( * )
*     ..
*     
*     Purpose
*     =======
*     
*     NEBLW2 computes the number of eigenvalues of a symmetric
*     tridiagonal matrix, T, which are less than or equal to V.  Uses
*     code from LAPACK's DSTEBZ and DLAEBZ.
*     
*     Arguments
*     =========
*     
*     N       (input) INTEGER
*     The dimension of the tridiagonal matrix T.  N >= 0.
*
*     V       (input) REAL
*     Compute the number of eigenvalues less than or equal to V.
*     
*     D       (input) REAL array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.  To
*          avoid overflow, the matrix must be scaled so that its largest
*          entry is no greater than overflow**(1/2) * underflow**(1/4)
*          in absolute value, and for greatest accuracy, it should not
*          be much smaller than that.
*
*  E       (input) REAL array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*          To avoid overflow, the matrix must be scaled so that its
*          largest entry is no greater than overflow**(1/2) *
*          underflow**(1/4) in absolute value, and for greatest
*          accuracy, it should not be much smaller than that.
*
*  WORK    (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*  =====================================================================
*
*     .. Parameters ..
      REAL   ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, MOUT
      REAL   PIVMIN, SAFEMN, TMP1, ULP
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
#include "blas_lapack.data"
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Check for Errors
*
      IF( N.LT.0 )
     $   INFO = -1
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NEWBLW2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      SNEBLW2 = 0
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
c      SAFEMN = SLAMCH( 'S' )
c      ULP = SLAMCH( 'P' )
c
      SAFEMN = SLAMCHS
      ULP = SLAMCHP

*
*     Special Case when N=1
*
      IF( N.EQ.1 ) THEN
         IF( D( 1 ) .LE. V ) THEN
            SNEBLW2 = 1
         ELSE
            SNEBLW2 = 0
         END IF
         RETURN
      END IF
*
*     Compute Splitting Points
*     Set WORK = E * E
*
      WORK( N ) = ZERO
      PIVMIN = ONE
*
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         IF( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   10 CONTINUE
      PIVMIN = PIVMIN*SAFEMN
*
*     Compute the number of eigenvalues in the initial intervals.
*
      TMP1 = D( 1 ) - V
      IF( ABS( TMP1 ).LT.PIVMIN )
     $   TMP1 = -PIVMIN
      MOUT = 0
      IF( TMP1.LE.ZERO )
     $   MOUT = 1
*
      DO 20 J = 2, N
         TMP1 = D( J ) - WORK( J-1 ) / TMP1 - V
         IF( ABS( TMP1 ).LT.PIVMIN )
     $      TMP1 = -PIVMIN
         IF( TMP1.LE.ZERO )
     $      MOUT = MOUT + 1
   20 CONTINUE
      SNEBLW2 = MOUT
*
*     End of NEWBLW2
*
      RETURN
      END
