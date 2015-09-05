      Subroutine spline_real (n, x, y, r)
      implicit none
#include "dimqm_constants.fh"
!  ============================================================================
!  Given arrays x(1:n) and y(1:n) containing a tabulated function,
!  i.e. y_i = f(x_i), with x_1 < x_2 ... < x_n, and given values yp1 and ypn
!  for the first derivative of the interpolation function at point 1 and n,
!  respectively, this routine returns an array r of length N that contains the
!  second derivatives of the interpolating function at the tabulated points x.
!  The interpolation is "natural", meaning that the second derivative on the
!  endpoints is zero.
!  Adopted from Numerical Recipes by Seth M. Morton, 2011
!  ============================================================================

       integer, Intent(in)  :: n
       double precision,   Intent(In)  :: x(n) ! Domain of data
       double precision,   Intent(In)  :: y(n) ! Range of rata (i.e. f(x))
       double precision,   Intent(Out) :: r(n) ! The second derivatives at data points
   
       double precision :: a(n)
       double precision :: b(n)
       double precision :: c(n)
   
       integer :: ier

!  Set up the tridiagonal equations
      c(1:n-1) = x(2:n) - x(1:n-1)
      r(1:n-1) = SIX * ( ( y(2:n) - y(1:n-1) ) / c(1:n-1) )
      r(2:n-1) = r(2:n-1) - r(1:n-1)
      a(2:n-1) = c(1:n-2)
      b(2:n-1) = TWO * ( c(2:n-1) + a(2:n-1) )
       b(1) = ONE
      b(n) = ONE

!  ---------------------------------------------------
!  The lower and upper boundaries are set to "natural"
!  ---------------------------------------------------

      r(1) = ZERO
      c(1) = ZERO
      r(n) = ZERO
      c(n) = ZERO

!  ---------------------------------------------------------------
!  Solve the tridiagonal system of equations with a LAPACK routine
!  ---------------------------------------------------------------

      call dgtsv (n, 1, a(2:n), b(1:n), c(1:n-1), r(1:n), n, ier)
      if (ier /= 0) then
        write(*,*) 'Error in spline'
        stop
      end if

      End Subroutine spline_real

      Subroutine spline (n, x, yr, yi, rr, ri)
      implicit none

!  =========================================================
!  An extention of spline_real to spline complex functions
!  Two real arrays representing real and imaginary are given
!  =========================================================

      Integer,            Intent(in)  :: n
      double precision,   Intent(In)  :: x(n)  ! Domain of data
      double precision,   Intent(In)  :: yr(n) ! Range of rata (i.e. f(x)) (real)
      double precision,   Intent(In)  :: yi(n) ! Range of rata (i.e. f(x)) (imag)
      double precision,   Intent(Out) :: rr(n) ! The second derivatives at data points
      double precision,   Intent(Out) :: ri(n) ! The second derivatives at data points

!  ------------------------------------
!  Spline real and imaginary separately
!  ------------------------------------

      call spline_real (n, x, yr, rr)
      call spline_real (n, x, yi, ri)

      End Subroutine spline

      Function interpolate_real (n, xa, ya, y2a, x) Result(y)
      implicit none
#include "dimqm_constants.fh"

!  =============================================================================
!  Given the arrays xa(1:n) and ya(1:n) of length n, which tabulates a function
!  (with the xa_i 's in order) and given the array y2a(1:n) which is output from
!  spline (above), and given a value of x, this routine returns a cubic-spline
!  interpolation value y.
!  iter was added to prevent a new search for each iteration in a sequencial
!  series of calls.  Use iter=1 for each call if they are unrelated from each
!  other.
!  Adopted from Numerical Recipes by Seth M. Morton, 2011
!  =============================================================================

      integer, Intent(In)  :: n
      double precision,   Intent(In)  :: xa(n)
      double precision,   Intent(In)  :: ya(n)
      double precision,   Intent(In)  :: y2a(n)
      double precision,   Intent(In)  :: x
      double precision                :: y

      integer :: klo
      integer :: khi
      double precision   :: a
      double precision   :: b
      double precision   :: h
      integer locate

!  -----------------------------
!  Find low and high index range
!  -----------------------------

      klo = MAX(MIN(locate(n, xa, x), n-1), 1)
      khi = klo + 1
   
      if (xa(klo) == xa(khi)) then
       write(*,*) 'bad xa input in interpolate'
       stop
      end if

!  -------------------
!  Evaluate the spline
!  -------------------

      h = xa(khi) - xa(klo)
      a = ( xa(khi) - x ) / h
      b = ( x - xa(klo) ) / h
      y = a * ya(klo) + b * ya(khi)
     $      + ( ( a**3 - a ) * y2a(klo) + ( b**3 - b ) * y2a(khi) )
     $      * ( h**2 ) / SIX

      End Function interpolate_real

      Function interpolate_complex1 (n, xa, yar, yai, y2ar, y2ai, x) 
     $         Result(y)
      implicit none
#include "dimqm_constants.fh"

!  ===============================================================
!  Extention of interpolate_real to interpolate a complex function
!  The original y-values are given in two separate real arrays
!  representing real and imaginary.
!  ===============================================================

      integer,  Intent(In)  :: n
      double precision,    Intent(In)  :: xa(n)
      double precision,    Intent(In)  :: yar(n)
      double precision,    Intent(In)  :: yai(n)
      double precision,    Intent(In)  :: y2ar(n)
      double precision,    Intent(In)  :: y2ai(n)
      double precision,    Intent(In)  :: x
      double complex              :: y

      integer :: klo
      integer :: khi
      double precision   :: a
      double precision   :: b
      double precision   :: h
      double precision   :: t
      integer locate

!  -----------------------------
!  Find low and high index range
!  -----------------------------`

      klo = MAX(MIN(locate(n, xa, x), n-1), 1)
      khi = klo + 1
   
      if (xa(klo) == xa(khi)) then
        write(*,*) 'bad xa input in interpolate'
        stop
      end if

!  -------------------
!  Evaluate the spline
!  -------------------

      h = xa(khi) - xa(klo)
      a = ( xa(khi) - x ) / h
      b = ( x - xa(klo) ) / h
      y = a * yar(klo) + b * yar(khi)
     $      + ( ( a**3 - a ) * y2ar(klo) + ( b**3 - b ) * y2ar(khi) )
     $      * ( h**2 ) / SIX
       t = a * yai(klo) + b * yai(khi)
     $       + ( ( a**3 - a ) * y2ai(klo) + ( b**3 - b ) * y2ai(khi) )
     $       * ( h**2 ) / SIX
       y = CMPLX(REAL(y), t)

      End Function interpolate_complex1

      Function interpolate_complex2 (n, xa, ya, y2a, x) Result(y)
      implicit none
#include "dimqm_constants.fh"

!  ===============================================================
!  Extention of interpolate_real to interpolate a complex function
!  The original y-values are given as complex.
!  ===============================================================

      integer,  Intent(In)  :: n
      double precision,    Intent(In)  :: xa(n)
      double complex, Intent(In)  :: ya(n)
      double complex, Intent(In)  :: y2a(n)
      double precision,    Intent(In)  :: x
      double complex              :: y

      integer :: klo
      integer :: khi
      double precision   :: a
      double precision   :: b
      double precision   :: h
      integer locate

!  -----------------------------
!  Find low and high index range
!  -----------------------------

      klo = MAX(MIN(locate(n, xa, x), n-1), 1)
      khi = klo + 1
   
      if (xa(klo) == xa(khi)) then
        write(*,*) 'bad xa input in interpolate'
        stop
      end if

!  -------------------
!  Evaluate the spline
!  -------------------

       h = xa(khi) - xa(klo)
       a = ( xa(khi) - x ) / h
       b = ( x - xa(klo) ) / h
       y = a * ya(klo) + b * ya(khi)
     $       + ( ( a**3 - a ) * y2a(klo) + ( b**3 - b ) * y2a(khi) )
     $       * ( h**2 ) / SIX_C

      End Function interpolate_complex2

      Function locate (n, xx, x) Result(i)
      implicit none
#include "dimqm_constants.fh"      

!  ==========================================================================
!  Given an array xx(1:n), and given a value x, returns a value j such that
!  x is between xx(j) and xx(j+1). xx must be monotonic, either increasing or
!  decreasing. j =0 or j =n is returned to indicate that x is out of range. 
!  Adopted from Numerical Recipes by Seth M. Morton, 2011
!  ==========================================================================

      integer, Intent(In)  :: n
      double precision,   Intent(In)  :: xx(n)
      double precision,   Intent(In)  :: x 
      integer              :: i

      integer :: jl
      integer :: jm
      integer :: ju 
      Logical       :: ascnd

!  ------------------------
!  Quick return if possible
!  ------------------------

      if (x == xx(1)) then
        i = 1
        return
      else if (x == xx(n)) then
        i = n - 1
        return
      end if

!  ---------------------------------
!  Initiallize limits and order flag
!  ---------------------------------

!  True if ascending order of table, false otherwise. 
      ascnd = xx(n) >= xx(1)
      jl = 0
      ju = n + 1

!  -----------------------------------------
!  Find the index for the value immeadiately
!  below or equal to the requested value
!  -----------------------------------------

!  Repeat until this condition is satisfied. 
      do while (ju - jl > 1)
!     Compute a midpoint, and replace either the lower limit
!     or the upper limit, as appropriate. 
        jm = ( ju + jl ) / 2
        if (ascnd .eqv. ( x >= xx(jm) )) then 
          jl = jm
        else 
          ju = jm
        end if 
      end do 

!  Set the output
      i = jl 

      End Function locate 
