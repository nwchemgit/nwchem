C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     GEAXPY -- Matrix Y <-- a*X + Y
C
C REVISION
C     $Id: geaxpy.f,v 2.2 1997-04-17 05:57:55 d3e129 Exp $
C
C SYNOPSIS
      Subroutine GEAXPY(M, N, Alpha, X, LDX, Y, LDY)
      Implicit NONE
      Integer M, N, LDX, LDY
      Double precision Alpha
      Double precision X(LDX, N), Y(LDY, N)
C
C ARGUMENTS
C     M       Row dimension of X, Y [IN]
C     N       Column dimension of X, Y [IN]
C     Alpha   Scale factor for X [IN]
C     X       MxN matrix [IN]
C     LDX     Leading dimension of X [IN]
C     Y       Result matrix [INOUT]
C     LDY     Leading dimension of Y [IN]
C
C DESCRIPTION
C     A simple convenience routine to simplify the addition of two
C     matrices which may have different leading dimensions (which makes
C     use of DAXPY unsafe).
C
C     Implemented as DAXPY calls on each column of X/Y.  Easy.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer J
C
      Do J = 1, N
         Call daxpy(M, Alpha, X(1, J), 1, Y(1, J), 1)
      EndDo
C
      Return
      End
