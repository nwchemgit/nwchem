      COMPLEX FUNCTION CLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX            X, Y
*     ..
*
c
* $Id: cladiv.f,v 1.1 1997-03-17 21:44:03 d3e129 Exp $
c
*  Purpose
*  =======
*
*  CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX
*  Y       (input) COMPLEX
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      CALL SLADIV( REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR,
     $             ZI )
      CLADIV = CMPLX( ZR, ZI )
*
      RETURN
*
*     End of CLADIV
*
      END
