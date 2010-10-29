C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     util_detriangularize -- convert from a triangular (I>=J) index to
C     (I, J)
C
C SYNOPSIS
      Subroutine Util_Detriangularize(IJ, I, J)
*
* $Id$
*
      Implicit NONE
#include "errquit.fh"
      Integer I, J, IJ
C
C DESCRIPTION
C     Converts a triangular index ( I*(I-1)/2 + J, I >= J) to I and J.
C     IJ must be positive.
C
C BUGS/LIMITATIONS
C     This conversion is a bit tricky because of the difference between
C     real and floating point math.  Brute-force testing gives a measure
C     of confidence that things are okay, and sanity checks make sure
C     we stay in the appropriate realm.
C
C     Tested with IEEE floating point arithmetic to the point where
C     2*IJ overflows a 32-bit integer (I ~ 46340).
C
C ROUTINES REQUIRED
      Double precision DLAMch
      External DLAMch
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Be careful that we're operating within bounds.  If 2*IJ overflows
C     an integer, it will probably appear negative, hence the check.
C
      If ( IJ .lt. 0 .OR. 2*IJ .lt. 0) Call ErrQuit(
     $   'Util_Detriangularize: Bad input/possible overflow', IJ,
     &       INPUT_ERR)
C
C     The DLAMch business results in the floating point number
C     closest to 1 o this machine.  It is necessary because the
C     sqrt often produces things with exactly .5 at the end, and
C     these need to be rounded down instead of up.

      I = NInt( Sqrt( Dble( 2 * IJ) ) * ( 1.0d0 - DLAMch('P') ) )
      J = IJ - I * (I-1) / 2
C
C     Be careful some more.  If you get the first message, something
C     is definitely wrong.  If you get the second, your input is
C     beyond what has actually been tested and you should check
C     carefully before using it.
C
      If ( I .lt. 1 .OR. J .lt. 1 ) Call ErrQuit(
     $   'Util_Detriangularize: Unreasonable result for', IJ,
     &       CALC_ERR)
      If ( I .gt. 46340 ) Call ErrQuit(
     $   'Util_Detriangularize: this is un-tested territory', IJ,
     &       UNKNOWN_ERR)
C
      Return
      End
