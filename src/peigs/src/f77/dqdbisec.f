*
* $Id$
*
      SUBROUTINE DQDBISEC( RANGE, N, B1, BN, VL, VU, IL, IU,
     $               L, D, LD, LLD,
     $     ABSTOL, SHIFT, M, W, WLEFT, WRIGHT,
     $               LPLUS,
     $               DPLUS, UMINUS, DMINUS, WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 0.0) -- in progress --
*     January 1996
*
*     .. Scalar Arguments ..
      CHARACTER          RANGE
      INTEGER            N, B1, BN, IL, IU, INFO, M
      DOUBLE PRECISION   VL, VU, ABSTOL, SHIFT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ),
     $                   W( * ), WLEFT( * ), WRIGHT( * ),
     $                   DMINUS( * ), LPLUS( * ), DPLUS( * ),
     $                   UMINUS( * ), WORK( * )
      INTEGER            IWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DQDBISEC computes the eigenvalues of L * D * L^T by bisection.
*  The inner loops are QD-recurrences that guarantee that the exact
*  eigenvalues of a relatively perturbed matrix are found.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix L * D * L^T.
*
*  B1      (input) INTEGER
*          Starting index of the submatrix (of L * D * L^T).
*
*  BN      (input) INTEGER
*          Last index of the submatrix (of L * D * L^T).
*
*  RANGE   (input) CHARACTER
*          = 'A': ("All")   all eigenvalues will be found.
*          = 'V': ("Value") all eigenvalues in the half-open interval
*                           (VL, VU] will be found.
*          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*                           entire matrix) will be found.
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues.  Eigenvalues less than or equal
*          to VL, or greater than VU, will not be returned.  VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of the bidiagonal matrix
*          L, in elements 1 to N-1.  L(N) need not be set.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D.
*
*  LD      (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*D(i).
*
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*L(i)*D(i).
*
*  ABSTOL  (input) DOUBLE PRECISION
*          Input tolerance.
*
*  M       (output) INTEGER
*          Number of eigenvalues found.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WLEFT   (output) DOUBLE PRECISION array, dimension (N)
*  WRIGHT  (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalue W(I) is contained in
*          [ WLEFT(I), WRIGHT(I) ]
*
*  LPLUS   (workspace) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of L+.
*
*  DPLUS   (workspace) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D+.
*
*  UMINUS  (workspace) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) superdiagonal elements of U-.
*
*  DMINUS  (workspace) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D-.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, ONE = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, KF, KFNEW, KL, KLNEW, firstind,
     $                   NGL, NGU, NLEFT, NMID, NRIGHT, save
      DOUBLE PRECISION   BNRM, DELTA, EPS, GL, GU, LEFT, MID,
     $                   RIGHT, TMP1, TMP2, TMP3, UNFL, STEP,
     $                   TMPD, TMPL, TMPT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT, ABS
      integer            doprnt1,doprnt2
      common             doprnt1,doprnt2

*     ..
*     .. Executable Statements ..
*

      M = 0
      DELTA = ZERO
      EPS = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'U' )
      IF( RANGE.EQ.'I' .AND. IL.EQ.B1 .AND. IU.EQ.BN )
     $     RANGE = 'A'
      
      IF( RANGE.EQ.'A') THEN
         IF( B1.EQ.1 ) THEN
            TMP3 = D( B1 )
         ELSE
            TMP3 = LLD( B1-1 ) + D( B1 )
         END IF
         GU = TMP3
         GL = TMP3
         TMP1 = ZERO
*     
         DO I = B1, BN-1
            TMP2 = ABS( LD( I ) )
            GU = MAX( GU, TMP3+TMP1+TMP2 )
            GL = MIN( GL, TMP3-TMP1-TMP2 )
            TMP3 = LLD( I ) + D( I+1 )
            TMP1 = TMP2
         END DO
*     
         GU = MAX( GU, TMP3+TMP1 )
         GL = MIN( GL, TMP3-TMP1 )
         BNRM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - 2.0d+0*BNRM*EPS*N
         GU = GU + 2.0d+0*BNRM*EPS*N
         NGL = 0
         NGU = N
      ELSE IF( RANGE.EQ.'V' ) THEN
         STEP = 10.0d+0
 112     CONTINUE
         GL = VL
         CALL DDSTQDS( GL, DELTA, N, B1, BN, L, D, LD,
     $        LPLUS, DPLUS, WORK )
         NGL = 0
         DO I = B1, BN
            IF( DPLUS( I ).LT.ZERO ) NGL = NGL + 1
         END DO
         
*     CALL DDQDS( GL, DELTA, N, B1, BN, L, D, LLD,
*     $                  UMINUS, DMINUS, WORK )
*     NGL = 0
*     DO I = B1, BN
*     IF( DMINUS( I ).LT.ZERO ) NGL = NGL + 1
*     END DO
         IF( NGL.GT.IL-1 ) THEN
            VL = VL - STEP * EPS * ABS( SHIFT )
            STEP = STEP * 10.0d+0
            GO TO 112
         ELSE IF( NGL.LT.IL-1 ) THEN
            print *,"NGL too small! NGL,IL = ",NGL,IL
            stop
         END IF
         
         STEP = 10.0d+0
 113     CONTINUE
         GU = VU
         CALL DDSTQDS( GU, DELTA, N, B1, BN, L, D, LD,
     $        LPLUS, DPLUS, WORK )
         NGU = 0
         DO I = B1, BN
            IF( DPLUS( I ).LT.ZERO ) NGU = NGU + 1
         END DO
         
*     CALL DDQDS( GU, DELTA, N, B1, BN, L, D, LLD,
*     $                  UMINUS, DMINUS, WORK )
*     NGU = 0
*     DO I = B1, BN
*     IF( DMINUS( I ).LT.ZERO ) NGU = NGU + 1
*     END DO
         IF( NGU.LT.IU ) THEN
            VU = VU + STEP * EPS * ABS( SHIFT )
            STEP = STEP * 10.0d+0
            GO TO 113
         ELSE IF( NGU.GT.IU ) THEN
            print *,"NGU too big! NGU,IU = ",NGU,IU
            stop
         END IF
         
         IF( NGL.GE.NGU ) THEN
            PRINT *,"NGL = ",NGL," should be smaller than NGU = ",NGU
            stop
         END IF
      END IF
      
      WORK( 1 ) = GL
      WORK( 2 ) = GU
      IWORK( 1 ) = NGL
      IWORK( 2 ) = NGU
      KF = 1
      KL = 1
110   CONTINUE
      IF( KL.GE.KF ) THEN
          KFNEW = KF
          KLNEW = KL
          DO I = KF, KL
             LEFT = WORK( 2*I-1 )
             RIGHT = WORK( 2*I )
             NLEFT = IWORK( 2*I-1 )
             NRIGHT = IWORK( 2*I )
             IF( RIGHT-LEFT.LE.MAX( MAX( 4.0d+0*EPS*MAX( ABS( LEFT ),
     $                    ABS( RIGHT ) ), UNFL ), ABSTOL ) ) THEN
                 IF( I.GT.KFNEW ) THEN
                     WORK( 2*I-1 ) = WORK( 2*KFNEW-1 )
                     WORK( 2*I ) = WORK( 2*KFNEW )
                     IWORK( 2*I-1 ) = IWORK( 2*KFNEW-1 )
                     IWORK( 2*I ) = IWORK( 2*KFNEW )
                     WORK( 2*KFNEW-1 ) = LEFT
                     WORK( 2*KFNEW ) = RIGHT
                     IWORK( 2*KFNEW-1 ) = NLEFT
                     IWORK( 2*KFNEW ) = NRIGHT
                 END IF
                 KFNEW = KFNEW + 1
             ELSE
                 MID = ( LEFT + RIGHT ) / 2.0d+0
                 NMID = 0

*                CALL DDSTQDS( MID, DELTA, N, B1, BN, L, D, LD,
*    $                         LPLUS, DPLUS, WORK( 2*N+1 ) )

                 IF( B1.EQ.1 ) THEN
                    TMPT = -MID
                 ELSE
                    TMPT = LLD( B1-1 ) - MID
                 END IF
                 DO J = B1, BN-1
                    TMPD = D( J ) + TMPT
                    IF( ( TMPD.LT.ZERO ).OR.( (TMPD.EQ.ZERO)
     $                 .AND. (ONE/TMPD.LT.ZERO) ) ) NMID = NMID + 1
                    TMPL = LD( J ) / TMPD
                    TMPT = TMPT * TMPL * L( J ) - MID
                 END DO
                 TMPD = D( BN ) + TMPT
                 IF( ( TMPD.LT.ZERO ).OR.( (TMPD.EQ.ZERO)
     $              .AND. (ONE/TMPD.LT.ZERO)) ) NMID = NMID + 1

                 IF(.NOT.( ( TMPD.GT.ZERO ) .OR.
     $                     ( TMPD.LT.ONE ) ) ) THEN
                    NMID = 0
                    IF( B1.EQ.1 ) THEN
                       TMPT = -MID
                    ELSE
                       TMPT = LLD( B1-1 ) - MID
                    END IF
                    DO J = B1, BN-1
                       TMPD = D( J ) + TMPT
                       IF( ( TMPD.LT.ZERO ).OR.( (TMPD.EQ.ZERO)
     $                    .AND. (ONE/TMPD.LT.ZERO)) ) NMID = NMID + 1
                       TMPL = LD( J ) / TMPD
*
*                      Need to check the next few lines
*
                       IF(TMPL.EQ.0 .AND. ONE/TMPT.EQ.ZERO ) THEN
                          TMPT = LLD( J ) - MID
                       ELSE IF( TMPT.EQ.ZERO .AND. ONE/TMPL.EQ.0 ) THEN
                          TMPT = -ONE / TMPD
                       ELSE
                          TMPT = TMPT * TMPL * L( J ) - MID
                       END IF
                    END DO
                    TMPD = D( BN ) + TMPT
                    IF( ( TMPD.LT.ZERO ).OR.( (TMPD.EQ.ZERO)
     $                 .AND. (ONE/TMPD.LT.ZERO)) ) NMID = NMID + 1
                 END IF

                 IF(.NOT.( ( TMPD.GT.ZERO ) .OR.
     $                     ( TMPD.LT.ONE ) ) ) THEN
                    print *,"DQDBISEC : NaN detected! Aborting!"
                    stop
                 END IF

                 if( nmid.lt.nleft .or. nmid.gt.nright) then
                    print *,"Count is nonmonotonic. Check!"
                    print *,"left,right,nleft,nright=",left,right,
     $                       nleft,nright
                    print *,"mid,nmid = ",mid,nmid
                    print *,"dplus = [ "
                    do j = b1,bn
                      write(*,101) dplus(j)
                    end do
                    print *,"];"
                    print *,"lplus = [ "
                    do j = b1,bn-1
                      write(*,101) lplus(j)
                    end do
                    print *,"];"
                    print *,"t  = [ "
                    do j = b1,bn
                      write(*,101) work(2*n+j)
                    end do
                    print *,"];"
*                   stop
                 end if
                 NMID = MAX( MIN( NRIGHT, NMID ), NLEFT )
                 IF( NMID.GT.NLEFT .AND. NMID.LT.NRIGHT ) THEN
                     WORK( 2*I ) = MID
                     IWORK( 2*I ) = NMID
                     KLNEW = KLNEW + 1
                     WORK( 2*KLNEW-1 ) = MID
                     IWORK( 2*KLNEW-1 ) = NMID
                     WORK( 2*KLNEW ) = RIGHT
                     IWORK( 2*KLNEW ) = NRIGHT
                 ELSE IF( NMID.GT.NLEFT ) THEN
                     WORK( 2*I ) = MID
                     IWORK( 2*I ) = NMID
                 ELSE IF( NMID.LT.NRIGHT ) THEN
                     WORK( 2*I-1 ) = MID
                     IWORK( 2*I-1 ) = NMID
                 END IF
             END IF
          END DO
          KF = KFNEW
          KL = KLNEW
          GO TO 110
      END IF
*
*     Sort the eigenvalues
*
      DO I = 1,KL
         DO J = IWORK( 2*I-1 )+1, IWORK( 2*I )
            W( J ) = ( WORK( 2*I-1 )+WORK( 2*I ) ) / 2.0d+0
            WLEFT( J ) = WORK( 2*I-1 )
            WRIGHT( J ) = WORK( 2*I )
            M = M + 1
         END DO
      END DO


      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*
*     End of DQDBISEC
*
      END

