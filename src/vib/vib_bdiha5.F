      SUBROUTINE vib_BDIHA5(NOINT,NDIM,IZMAT,CARTC,DIHANG,BMAT)
* $Id$
C*MODULE VIB  *DECK BDIHA5
C* 21 MAR 1983 - STE * 9 MAR 1983 - STE
C********************************************************************
C*       ROUTINE BDIHA5(IW,NOINT,NDIM,IZMAT,CARTC,DIHANG,BMAT)
C*
C*    AUTHOR: S. T. ELBERT (AMES LABORATORY - USDOE)  DATE: 25 FEB 1983
C*
C*    PURPOSE -
C*       EVALUATE B MATRIX AND DIHEDRAL ANGLE BETWEEN 5 ATOMS FORMING
C*                *            ***      *             *
C*       TWO PLANES HAVING ONE ATOM IN COMMON
C*
C*    ON ENTERING -
C*       IW     - INTEGER (LOGICAL UNIT NO.)
C*                ERROR MESSAGE FILE
C*       NOINT  - INTEGER
C*                NO. OF INTERNAL COORDINATE BEING EVALUATED
C*       NDIM   - INTEGER
C*                NO. OF COORDINATES = 3 * NATOMS
C*       IZMAT  - INTEGER (5)
C*                LIST OF ATOMS INVOLVED
C*       CARTC  - W.P. REAL (3,NATOMS)
C*                CARTESIAN COORDINATES FOR ATOMS
C*
C*    ON EXITING -
C*       DIHANG - W.P. REAL
C*                DIHEDRAL ANGLE IN DEGREES
C*       BMAT   - W.P. REAL (NDIM,NDIM)
C*                B MATRIX (TRANSFORMS INTERNAL COORD TO CART. COORD
C*                          TO FIRST ORDER)
C*
C********************************************************************
C
C                    DECLARATIONS
C
      IMPLICIT NONE 
C
      integer noint, ndim
      integer izmat(5)
      double precision cartc(3,*)
      double precision bmat(ndim,ndim)
      double precision dihang
c
c:old      DIMENSION IZMAT(5), CARTC(3,1), BMAT(NDIM,NDIM)
c:old     *         ,A(3),B(3),C(3),D(3),E1(3),E2(3),E3(3)
C
      double precision one, zero
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
c
      double precision A(3),B(3),C(3),D(3)
      double precision E1(3), E2(3), E3(3)
      double precision adote2, cdote1, e1de2, sindi
      double precision axbmag, cxdmag, axbmi, cxdmi
      double precision f1, f2, f4, f5
      double precision b1, b2, b3, b4, b5
      integer iatom, jatom, katom, latom, matom, ixyz
      integer i, j, k, l, m, n
      
C
C---------------------------------------------------------------------
C
C           THE FIRST THREE ATOMS ARE ON THE FIRST PLANE
C           THE  LAST THREE ATOMS ARE ON THE SECOND PLANE
C
      IATOM = IZMAT(1)
      JATOM = IZMAT(2)
      KATOM = IZMAT(3)
      LATOM = IZMAT(4)
      MATOM = IZMAT(5)
C
C           DEFINE VECTORS FROM PIVOT (THIRD) ATOM THAT DEFINE PLANES
C                    A = I - K   AND  B = J - K   DEFINE FIRST PLANE
C                    C = L - K   AND  D = M - K   DEFINE SECOND PLANE
C
      DO 110 IXYZ=1,3
         A(IXYZ) = CARTC(IXYZ,IATOM) - CARTC(IXYZ,KATOM)
         B(IXYZ) = CARTC(IXYZ,JATOM) - CARTC(IXYZ,KATOM)
         C(IXYZ) = CARTC(IXYZ,LATOM) - CARTC(IXYZ,KATOM)
         D(IXYZ) = CARTC(IXYZ,MATOM) - CARTC(IXYZ,KATOM)
  110 CONTINUE
C
C              E1 = A X B  IS PERPENDICULAR TO FIRST PLANE
C
      E1(1) = A(2)*B(3) - A(3)*B(2)
      E1(2) = A(3)*B(1) - A(1)*B(3)
      E1(3) = A(1)*B(2) - A(2)*B(1)
      AXBMAG = SQRT(E1(1)**2 + E1(2)**2 + E1(3)**2)
      IF(ONE+AXBMAG.EQ.ONE) GO TO 810
      AXBMI = ONE / AXBMAG
C
C              E2 = C X D  IS PERPENDICULAR TO FIRST PLANE
C
      E2(1) = C(2)*D(3) - C(3)*D(2)
      E2(2) = C(3)*D(1) - C(1)*D(3)
      E2(3) = C(1)*D(2) - C(2)*D(1)
      CXDMAG = SQRT(E2(1)**2 + E2(2)**2 + E2(3)**2)
      IF(ONE+CXDMAG.EQ.ONE) GO TO 820
      CXDMI = ONE / CXDMAG
C
C              NORMALIZE E1 AND E2
C
      DO 120 IXYZ=1,3
         E1(IXYZ) = AXBMI * E1(IXYZ)
         E2(IXYZ) = CXDMI * E2(IXYZ)
  120 CONTINUE
C
C        THE DIRECTION OF E2 IS DEFINED SUCH THAT WHEN THE PLANES ARE
C        ROTATED ALONG THE INTERSECTION E3 (= E1 X E2) SO THEY COINCIDE
C        (DIHEDRAL ANGLE = 0), THEN ATOMS 1 AND 4 LIE ON THE SAME
C        SIDE OF E3.  IF ATOMS 1 AND/OR 4 LIE ON E3, THEN
C        ATOMS 2 AND/OR 5 ARE USED.
C
      ADOTE2 = A(1)*E2(1) + A(2)*E2(2) + A(3)*E2(3)
      IF(ONE+ADOTE2.EQ.ONE) ADOTE2 = B(1)*E2(1)+B(2)*E2(2)+B(3)*E2(3)
      CDOTE1 = C(1)*E1(1) + C(2)*E1(2) + C(3)*E1(3)
      IF(ONE+CDOTE1.EQ.ONE) CDOTE1 = D(1)*E1(1)+D(2)*E1(2)+D(3)*E1(3)
      IF(ADOTE2*CDOTE1 .LT. ZERO) GO TO 140
C
C        CHANGE SENSE OF E2 SO A (OR B) WILL BE ON SAME SIDE
C        OF E1 X E2 AS C (OR D)
C
         ADOTE2=-ADOTE2
         DO 130 IXYZ=1,3
            E2(IXYZ) = -E2(IXYZ)
  130    CONTINUE
  140 CONTINUE
C
C           EVALUATE DIHEDRAL ANGLE    0 < ABS(DIHANG) < 180
C           ANGLE IS POSITIVE FOR SECOND PLANE "ABOVE" ATOM 1
C           ANGLE IS NEGATIVE FOR SECOND PLANE "BELOW" ATOM 1
C
      E1DE2 = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
      IF(ABS(E1DE2).GT.ONE) GO TO 840
      IF(ABS(E1DE2).EQ.ONE) GO TO 830
      DIHANG = ACOS(E1DE2)
      IF(ADOTE2 .GT. ZERO) DIHANG = -DIHANG
      SINDI = ONE / SIN(DIHANG)
C
C         E3 = E1 X E2  IS ALONG THE INTERSECTION LINE OF THE TWO PLANES
C
      E3(1) = SINDI * (E1(2)*E2(3) - E1(3)*E2(2))
      E3(2) = SINDI * (E1(3)*E2(1) - E1(1)*E2(3))
      E3(3) = SINDI * (E1(1)*E2(2) - E1(2)*E2(1))
C
      F1 =  AXBMI * (B(1)*E3(1) + B(2)*E3(2) + B(3)*E3(3))
      F2 = -AXBMI * (A(1)*E3(1) + A(2)*E3(2) + A(3)*E3(3))
      F4 = -CXDMI * (D(1)*E3(1) + D(2)*E3(2) + D(3)*E3(3))
      F5 =  CXDMI * (C(1)*E3(1) + C(2)*E3(2) + C(3)*E3(3))
C
      I = 3*IATOM - 3
      J = 3*JATOM - 3
      K = 3*KATOM - 3
      L = 3*LATOM - 3
      M = 3*MATOM - 3
C
C           CALCULATE B MATRIX ELEMENTS
C
      DO 160 N=1,3
         B1 = E1(N) * F1
         B2 = E1(N) * F2
         B4 = E2(N) * F4
         B5 = E2(N) * F5
         B3 = -(B1+B2+B4+B5)
         BMAT(NOINT,I+N) = B1
         BMAT(NOINT,J+N) = B2
         BMAT(NOINT,K+N) = B3
         BMAT(NOINT,L+N) = B4
         BMAT(NOINT,M+N) = B5
  160 CONTINUE
      RETURN
C
C               ERROR MESSAGE SECTION
C
  810 CONTINUE
      WRITE(6,910)
      GO TO 880
C
  820 CONTINUE
      WRITE(6,920)
      GO TO 880
C
  830 CONTINUE
      WRITE(6,930)
      GO TO 880
C
  840 CONTINUE
      WRITE(6,940)
C
  880 CONTINUE
      WRITE(6,901) NOINT,IATOM,(CARTC(IXYZ,IATOM),IXYZ=1,3)
     *                   ,JATOM,(CARTC(IXYZ,JATOM),IXYZ=1,3)
     *                   ,KATOM,(CARTC(IXYZ,KATOM),IXYZ=1,3)
     *                   ,LATOM,(CARTC(IXYZ,LATOM),IXYZ=1,3)
     *                   ,MATOM,(CARTC(IXYZ,MATOM),IXYZ=1,3)
      WRITE(6,900)
C
  900 FORMAT('0THESE MESSAGES BROUGHT TO YOU BY ROUTINE BDIHA5 *****')
  901 FORMAT('0INTERNAL COORDINATE',I4,' (DIHEDRAL ANGLE BETWEEN TWO'
     *      ,' PLANES WITH ONE ATOM IN COMMON) USES THESE ATOMS'/
     *       (I10,3F20.6))
  910 FORMAT('0***** ERROR - FIRST THREE ATOMS ARE LINEAR *****')
  920 FORMAT('0***** ERROR -- LAST THREE ATOMS ARE LINEAR *****')
  930 FORMAT('0***** ERROR - ALL FIVE ATOMS LIE IN THE SAME PLANE ****')
  940 FORMAT('0***** ERROR - COSIN OF DIHEDRAL ANGLE IS GREATER THAN 1')
      END
