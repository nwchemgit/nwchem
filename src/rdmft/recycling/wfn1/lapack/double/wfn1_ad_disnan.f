      LOGICAL FUNCTION WFN1_AD_DISNAN( DIN )
      USE WFN1_AD1
      IMPLICIT NONE
#include "lapack/double/intf_wfn1_ad_dlaisnan.fh"
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*     .. Scalar Arguments ..
      TYPE(WFN1_AD_DBLE) :: DIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN     (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
*  ..
*  .. Executable Statements ..
      WFN1_AD_DISNAN = WFN1_AD_DLAISNAN(DIN,DIN)
      RETURN
      END
