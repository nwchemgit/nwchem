
      SUBROUTINE DC_Initialize
C******************************************************************

C This program initializes all variables
C for the psum functions.

      IMPLICIT NONE

#include "dc_common.fh"

      COMMON/DCI/kb,k,e,esqrd,PI,Q(4),
     +QSQRD(4,4),SGMSQ(4,4),EPS(4,4)
      
      DOUBLE PRECISION kb,k,e,esqrd,PI
      DOUBLE PRECISION Q,SIGMA(4),EPSILON(4)
      DOUBLE PRECISION QSQRD,SGMSQ,EPS
C      DOUBLE PRECISION R ! in include file
C      INTEGER NMOLS
C      INTEGER MXATOM ! in include file
C      INTEGER MXMOLS

C      NMOLS = 2 ! get NMOLS from input file

      kb = 1.987905D-03
      k = 1.293601055D+40
      e = 1.60219D-19
      PI = 3.141592653589793D+00
      esqrd = e*e  
     
      Q(1) = 0.0000
      SIGMA(1) = 3.2340
      EPSILON(1) = 0.1825
      
      Q(2) = 0.5190
      SIGMA(2) = 0.0000
      EPSILON(2) = 0.0000

      Q(3) = 0.5190
      SIGMA(3) = 0.0000
      EPSILON(3) = 0.0000

      Q(4) = -1.0380
      SIGMA(4) = 0.0000
      EPSILON(4) = 0.0000

      QSQRD(1,1) = Q(1)*Q(1)
      QSQRD(2,2) = Q(2)*Q(2)
      QSQRD(3,3) = Q(3)*Q(3)
      QSQRD(4,4) = Q(4)*Q(4)
      QSQRD(1,2) = Q(1)*Q(2)
      QSQRD(2,1) = QSQRD(1,2)
      QSQRD(1,3) =  Q(1)*Q(3)
      QSQRD(3,1) = QSQRD(1,3)
      QSQRD(1,4) = Q(1)*Q(4)
      QSQRD(4,1) = QSQRD(1,4)
      QSQRD(2,3) = Q(2)*Q(3)
      QSQRD(3,2) = QSQRD(2,3)
      QSQRD(2,4) = Q(2)*Q(4)
      QSQRD(4,2) = QSQRD(2,4)
      QSQRD(3,4) = Q(3)*Q(4)
      QSQRD(4,3) = QSQRD(3,4)

c NOTE:  Only the Oxygen's parameters matter
c Still not sure where the 1/2's get carried
c to in all the cubing etc in
c the Vljones subrountine

      SGMSQ(1,1) = SIGMA(1)*SIGMA(1)
      SGMSQ(2,2) = SIGMA(2)*SIGMA(2)
      SGMSQ(3,3) = SIGMA(3)*SIGMA(3)
      SGMSQ(4,4) = SIGMA(4)*SIGMA(4)
      SGMSQ(1,2) = 0.5D+00*( SIGMA(1) + SIGMA(2) )
     +*0.5D+00*( SIGMA(1) + SIGMA(2) )
      SGMSQ(2,1) = SGMSQ(1,2)
      SGMSQ(1,3) = 0.5D+00*( SIGMA(1) + SIGMA(3) )
     +*0.5D+00*( SIGMA(1) + SIGMA(3) )
      SGMSQ(3,1) = SGMSQ(1,3)
      SGMSQ(1,4) = 0.5D+00*( SIGMA(1) + SIGMA(4) )
     +*0.5D+00*( SIGMA(1) + SIGMA(4) )
      SGMSQ(4,1) = SGMSQ(1,4)
      SGMSQ(2,3) = 0.5D+00*( SIGMA(2) + SIGMA(3) )
     +*0.5D+00*( SIGMA(2) + SIGMA(3) )
      SGMSQ(3,2) = SGMSQ(2,3)
      SGMSQ(2,4) = 0.5D+00*( SIGMA(2) + SIGMA(4) )
     +*0.5D+00*( SIGMA(2) + SIGMA(4) )
      SGMSQ(4,2) = SGMSQ(2,4)
      SGMSQ(3,4) = 0.5D+00*( SIGMA(3) + SIGMA(4) )
     +*0.5D+00*( SIGMA(3) + SIGMA(4) )
      SGMSQ(4,3) = SGMSQ(3,4)

      EPS(1,1) = EPSILON(1)
      EPS(2,2) = EPSILON(2)
      EPS(3,3) = EPSILON(3)
      EPS(4,4) = EPSILON(4)
      EPS(1,2) = DSQRT(EPSILON(1)*EPSILON(2))
      EPS(2,1) = EPS(1,2)
      EPS(1,3) = DSQRT(EPSILON(1)*EPSILON(3))
      EPS(3,1) = EPS(1,3)
      EPS(1,4) = DSQRT(EPSILON(1)*EPSILON(4))
      EPS(4,1) = EPS(1,4)
      EPS(2,3) = DSQRT(EPSILON(2)*EPSILON(3))
      EPS(3,2) = EPS(2,3)
      EPS(2,4) = DSQRT(EPSILON(2)*EPSILON(4))
      EPS(4,2) = EPS(2,4)
      EPS(3,4) = DSQRT(EPSILON(3)*EPSILON(4))
      EPS(4,3) = EPS(3,4)

C The initial cluster coordinates for a fixed center
C of mass at the origin are:

C      R(1,1,1) =  0.39367650361336998
C      R(1,1,2) =  -1.7437946269568274
C      R(1,1,3) =  -0.76229191212927105
C      R(1,2,1) =  -0.42722715712577775
C      R(1,2,2) =  -1.2791388125263257
C      R(1,2,3) =  -0.92489827978131922
C      R(1,3,1) =  1.0754639527170682
C      R(1,3,2) = -1.0958839290750606
C      R(1,3,3) = -0.94007345986422264
C      R(1,4,1) = 0.36815081106402403
C      R(1,4,2) = -1.539655724666992
C      R(1,4,3) = -0.82474787735099175
C      R(2,1,1) = -0.38589075355937486
C      R(2,1,2) =  1.6193569594453325
C      R(2,1,3) =  0.80220298850272087
C      R(2,2,1) =  -0.52480281587511424
C      R(2,2,2) =  1.9381013476193931
C      R(2,2,3) =  -8.9612088097092979E-2
C      R(2,3,1) =  -0.2469786912436355
C      R(2,3,2) =  2.4115051693716021
C      R(2,3,3) =  1.3212726652839746
C      R(2,4,1) =  -0.38589075355937486
C      R(2,4,2) =  1.8231887241339813
C      R(2,4,3) =  0.73380992084296937

C INITIAL UU =  -0.14686911198510338  kcal/mole

      END
c $Id$
