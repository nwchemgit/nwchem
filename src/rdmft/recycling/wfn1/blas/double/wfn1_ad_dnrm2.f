*> \brief \b DNRM2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    DNRM2 := sqrt( x'*x )
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to DLASSQ.
*>     Sven Hammarling, Nag Ltd.
*> \endverbatim
*>
*  =====================================================================
      FUNCTION WFN1_AD_DNRM2(N,X,INCX)
      USE WFN1_AD1
      IMPLICIT NONE
      TYPE(WFN1_AD_DBLE) :: WFN1_AD_DNRM2
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      TYPE(WFN1_AD_DBLE) :: X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      TYPE(WFN1_AD_DBLE) :: ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
c     INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*
*         The following loop is equivalent to this call to the LAPACK
*         auxiliary routine:
*         CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
*         The following loop also causes severe problems when using 
*         automatic differentiation to compute the gradient of the norm.
*         In particular if the vector starts with a few small numbers, e.g.
*           ( 1.0d-50 )
*           ( 2.0d-50 )
*           ( 1.0d+00 )
*         then differentiating the (SCALE/ABSXI) factor generates massive
*         that introduce major absolute errors in the gradients.
*
*         The scale factor approach is needed, however, to avoid/minimize
*         underflow and overflow issues associated with squaring the vector
*         elements. The general idea is to scale the largest vector element
*         to 1, then compute the norm and scale the result back. The original
*         implementation aimed to do this in a single pass through the data
*         to minimize memory access. To avoid the numerical problems in the
*         gradient evaluation we have to establish the maximum vector
*         element in a first pass, and then compute the norm in a second
*         pass.
*
          do ix = 1, 1+(n-1)*incx, incx
            scale = max(scale,abs(x(ix)))
          enddo
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
*                 IF (SCALE.LT.ABSXI) THEN
*                     SSQ = ONE + SSQ* (SCALE/ABSXI)**2.0d0
*                     SCALE = ABSXI
*                 ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2.0d0
*                 END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      WFN1_AD_DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
c $Id: wfn1_ad_dnrm2.f,v 1.3 2015/03/23 03:17:37 D3Y133 Exp D3Y133 $
