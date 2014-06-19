!> \ingroup nwad
!> @{
!>
!> A module implementing Automatic Differentiation [1] capabilities in 
!> Fortran 90. The main aim is to enable calculating derivatives of existing
!> capabilities with minimal code changes. In particular executable statements
!> should not have to be changed.
!>
!> A target application is the calculation of derivatives of Density
!> Functionals, where this capability at least provides baselines for functional
!> implementations generated with symbolic mathematics packages.
!>
!> The package is designed to generate compact data structures to keep the
!> overhead of copying data low. In addition applications are targetted where
!> all partial derivatives of an expression are needed. The method chosen 
!> requires multiple executions of the function to obtain all its derivatives.
!> Every function execution computes only one specific partial derivative.
!>
!> ### References ###
!>
!> [1] See e.g. <a href="http://www.autodiff.org/">www.autodiff.org</a>
!>
!> $Id: $
!>
!> Huub van Dam, 2014
!>
module nwad
  !>
  !> The data type to hold a double precision value and the derivatives of this
  !> expression with respect to the sum of the active variables. I.e. the
  !> member \f$\mathrm{d}n\f$ holds
  !> \f{eqnarray}{
  !>   \mathrm{d}n &=& \frac{\mathrm{d}^n f}{\mathrm{d}(\sum_i x_i)^n}
  !> \f}
  !> where \f$x_i\f$ are all active variables.
  !>
  type :: nwad_dble
    double precision :: d0
    double precision :: d1
    double precision :: d2
    double precision :: d3
  end type nwad_dble
  interface assignment (=)
    module procedure nwad_dble_assign
  end interface
  interface operator (+)
    module procedure nwad_dble_add
    module procedure nwad_dble_addx
    module procedure nwad_dble_addy
  end interface
  interface operator (-)
    module procedure nwad_dble_sub
    module procedure nwad_dble_subx
    module procedure nwad_dble_suby
  end interface
  interface operator (*)
    module procedure nwad_dble_mult
    module procedure nwad_dble_multx
    module procedure nwad_dble_multy
  end interface
  interface operator (/)
    module procedure nwad_dble_div
    module procedure nwad_dble_divx
    module procedure nwad_dble_divy
  end interface
  interface operator (**)
    module procedure nwad_dble_pow
    module procedure nwad_dble_powx
    module procedure nwad_dble_powy
  end interface
  interface abs
    module procedure nwad_dble_abs
  end interface
  interface exp
    module procedure nwad_dble_exp
  end interface
  interface sqrt
    module procedure nwad_dble_sqrt
  end interface
  interface log
    module procedure nwad_dble_log
  end interface
! interface log10
!   module procedure nwad_dble_log10
! end interface
  interface sin
    module procedure nwad_dble_sin
  end interface
  interface cos
    module procedure nwad_dble_cos
  end interface
  interface tan
    module procedure nwad_dble_tan
  end interface
  interface asin
    module procedure nwad_dble_asin
  end interface
  interface acos
    module procedure nwad_dble_acos
  end interface
  interface atan
    module procedure nwad_dble_atan
  end interface
  interface sinh
    module procedure nwad_dble_sinh
  end interface
  interface cosh
    module procedure nwad_dble_cosh
  end interface
  interface tanh
    module procedure nwad_dble_tanh
  end interface
  interface asinh
    module procedure nwad_dble_asinh
  end interface
  interface erf
    module procedure nwad_dble_erf
  end interface
  interface erfc
    module procedure nwad_dble_erfc
  end interface
  interface active
    module procedure nwad_dble_active
  end interface
  interface active_neg
    module procedure nwad_dble_active_neg
  end interface
  interface inactive
    module procedure nwad_dble_inactive
  end interface
contains
  !>
  !> \brief Assign a value to an inactive variable
  !>
  !> Assign a floating point value to an nwad variable. This operation 
  !> generates an inactive variable. I.e. a variable for which no derivatives
  !> are calculated. In practice it means that the components are initialized as
  !> \f{eqnarray*}{
  !>   d0 &=& x \\\\
  !>   d1 &=& 0 \\\\
  !>   d2 &=& 0 \\\\
  !>   d3 &=& 0 \\\\
  !> \f}
  !>
  subroutine nwad_dble_assign(s,x)
    double precision, intent(in) :: x
    type(nwad_dble), intent(out) :: s
    s%d0 = x
    s%d1 = 0
    s%d2 = 0
    s%d3 = 0
  end subroutine nwad_dble_assign
  !>
  !> \brief Evaluate the addition operator and its derivatives
  !>
  !> The implementation of the addition operator. The chain rule is used to
  !> evaluate the derivatives. I.e. we calculate \f$s(r) = x(r) + y(r)\f$
  !> and differentiate this expression as
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     + \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     + \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}  
  !>     + \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>       \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}  
  !>     + \frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3} \\\\
  !> \f}
  !>
  function nwad_dble_add(x,y) result (s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    type(nwad_dble)             :: s
    s%d0 = x%d0 + y%d0
    s%d1 = x%d1 + y%d1
    s%d2 = x%d2 + y%d2
    s%d3 = x%d3 + y%d3
  end function nwad_dble_add
  !>
  !> \brief Evaluate the addition where \f$y\f$ is inactive
  !>
  !> This routine does the same as nwad_dble_add but \f$y\f$ is inactive
  !> and therefore all its derivatives are 0.
  !>
  function nwad_dble_addx(x,y) result (s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x%d0 + y
    s%d1 = x%d1
    s%d2 = x%d2
    s%d3 = x%d3
  end function nwad_dble_addx
  !>
  !> \brief Evaluate the addition where \f$x\f$ is inactive
  !>
  !> This routine does the same as nwad_dble_add but \f$x\f$ is inactive
  !> and therefore all its derivatives are 0.
  !>
  function nwad_dble_addy(x,y) result (s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x + y%d0
    s%d1 =     y%d1
    s%d2 =     y%d2
    s%d3 =     y%d3
  end function nwad_dble_addy
  !>
  !> \brief Evaluation of the subtraction operator
  !>
  !> The implementation of the subtraction operator. The chain rule is used to
  !> evaluate the derivatives. Obviously this is very similar to the 
  !> addition operator apart from the minus sign.
  !>
  function nwad_dble_sub(x,y) result (s)
    type(nwad_dble), intent(in) :: x, y
    type(nwad_dble)             :: s
    s%d0 = x%d0 - y%d0
    s%d1 = x%d1 - y%d1
    s%d2 = x%d2 - y%d2
    s%d3 = x%d3 - y%d3
  end function nwad_dble_sub
  !>
  !> \brief Evaluation of the subtraction operator where \f$y\f$ is inactive
  !>
  !> This function is similar to nwad_dble_sub but because \f$y\f$ is in
  !> inactive all derivatives of \f$y\f$ vanish.
  !>
  function nwad_dble_subx(x,y) result (s)
    type(nwad_dble), intent(in)  :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x%d0 - y
    s%d1 = x%d1
    s%d2 = x%d2
    s%d3 = x%d3
  end function nwad_dble_subx
  !>
  !> \brief Evaluation of the subtraction operator where \f$x\f$ is inactive
  !>
  !> This function is similar to nwad_dble_sub but because \f$x\f$ is in
  !> inactive all derivatives of \f$x\f$ vanish.
  !>
  function nwad_dble_suby(x,y) result (s)
    double precision, intent(in) :: x
    type(nwad_dble), intent(in)  :: y
    type(nwad_dble)              :: s
    s%d0 = x - y%d0
    s%d1 =   - y%d1
    s%d2 =   - y%d2
    s%d3 =   - y%d3
  end function nwad_dble_suby
  !>
  !> \brief Evaluate the multiplicition operator and its derivatives
  !>
  !> The implementation of the multiplication operator. The chain rule is used
  !> to evaluate the derivatives. I.e. the derivatives of 
  !> \f$s(r) = x(r)*y(r)\f$ are evaluated as
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} + 
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}  
  !>     * \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} +
  !>      2\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} +
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>       \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}  
  !>     * \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} +
  !>      3\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} +
  !>      3\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} +
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3} \\\\
  !> \f}
  !>
  function nwad_dble_mult(x,y) result (s)
    type(nwad_dble), intent(in) :: x, y
    type(nwad_dble)             :: s
    s%d0 =    x%d0 * y%d0
    s%d1 =    x%d1 * y%d0 +     x%d0 * y%d1
    s%d2 =    x%d2 * y%d0 + 2 * x%d1 * y%d1 + &
              x%d0 * y%d2
    s%d3 =    x%d3 * y%d0 + 3 * x%d2 * y%d1 + &
           3* x%d1 * y%d2 +     x%d0 * y%d3
  end function nwad_dble_mult
  !>
  !> \brief Evaluate the multiplication operator where \f$y\f$ is inactive
  !>
  !> This is similar to the regular multiplication operator except that all
  !> derivative of \f$y\f$ are zero.
  !>
  function nwad_dble_multx(x,y) result (s)
    type(nwad_dble), intent(in)  :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x%d0 * y
    s%d1 = x%d1 * y
    s%d2 = x%d2 * y
    s%d3 = x%d3 * y
  end function nwad_dble_multx
  !>
  !> \brief Evaluate the multiplication operator where \f$x\f$ is inactive
  !>
  !> This is similar to the regular multiplication operator except that all
  !> derivative of \f$x\f$ are zero.
  !>
  function nwad_dble_multy(x,y) result (s)
    double precision, intent(in) :: x
    type(nwad_dble), intent(in)  :: y
    type(nwad_dble)              :: s
    s%d0 = x * y%d0
    s%d1 = x * y%d1
    s%d2 = x * y%d2
    s%d3 = x * y%d3
  end function nwad_dble_multy
  !>
  !> \brief Evaluate the division operator
  !>
  !> The implementation of the division operator. The chain rule is used to
  !> evaluate the derivatives. I.e. if \f$s(r) = x(r)/y(r)\f$ then we have
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>       \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     / \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     / \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} 
  !>     - \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0} 
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \frac{1}{y(r)^2} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}  
  !>     / \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} 
  !>     -2\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \frac{1}{y(r)^2}
  !>     - \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} \frac{1}{y(r)^2} 
  !>     +2\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^2
  !>       \frac{1}{y(r)^3}  \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>       \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}  
  !>     / \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} 
  !>     -3\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \frac{1}{y(r)^2}
  !>     +6\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^2
  !>       \frac{1}{y(r)^3}
  !>     -6\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^3
  !>       \frac{1}{y(r)^4}
  !>     -3\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}  
  !>     * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>       \frac{1}{y(r)^2}
  !>     +6\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>     * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>       \frac{1}{y(r)^3}
  !>     - \frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}  
  !>     * \frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3}
  !>       \frac{1}{y(r)^2} \\\\
  !> \f}
  !>
  function nwad_dble_div(x,y) result (s)
    type(nwad_dble), intent(in) :: x, y
    type(nwad_dble)             :: s
    s%d0 = x%d0 / y%d0
    s%d1 = x%d1 / y%d0               -     x%d0 * y%d1 / (y%d0 ** 2)
    s%d2 = x%d2 / y%d0               - 2 * x%d1 * y%d1 / (y%d0 ** 2) &
         - x%d0 * y%d2 / (y%d0 ** 2) + 2 * x%d0 * y%d1 ** 2 / (y%d0 ** 3)
    s%d3 =     x%d3               /  y%d0       &
         - 3 * x%d2 * y%d1        / (y%d0 ** 2) &
         + 6 * x%d1 * y%d1 ** 2   / (y%d0 ** 3) &
         - 6 * x%d0 * y%d1 ** 3   / (y%d0 ** 4) &
         - 3 * x%d1 * y%d2        / (y%d0 ** 3) &
         + 6 * x%d0 * y%d1 * y%d2 / (y%d0 ** 3) &
         -     x%d0 * y%d3        / (y%d0 ** 2)
  end function nwad_dble_div
  !>
  !> \brief Evaluate the division operator where \f$y\f$ is inactive
  !>
  !> This function is particularly simple. It is essentially the same as
  !> multiplying with \f$1/y\f$ while \f$y\f$ is inactive.
  !>
  function nwad_dble_divx(x,y) result (s)
    type(nwad_dble), intent(in)  :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x%d0 / y
    s%d1 = x%d1 / y
    s%d2 = x%d2 / y
    s%d3 = x%d3 / y
  end function nwad_dble_divx
  !>
  !> \brief Evaluate the division operator where \f$x\f$ is inactive
  !>
  !> Similar to nwad_dble_div but now \f$x\f$ is inactive. I.e. we consider
  !> \f$s(r) = x/y(r)\f$:
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>       x / \frac{\mathrm{d}^0y(r)}{\mathrm{d}r^0} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       x * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \frac{1}{y(r)^2} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       2 * x * \left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^2
  !>               \frac{1}{y(r)^3}
  !>     -     x * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>               \frac{1}{y(r)^2} \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>     -     x * \frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3}
  !>               \frac{1}{y(r)^2} 
  !>     + 6 * x * \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>             * \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>               \frac{1}{y(r)^3} 
  !>     - 6 * x * \left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^3
  !>               \frac{1}{y(r)^4} 
  !> \f}
  !>
  function nwad_dble_divy(x,y) result (s)
    double precision, intent(in) :: x
    type(nwad_dble), intent(in)  :: y
    type(nwad_dble)              :: s
    s%d0 =   x / y%d0
    s%d1 = - x * y%d1 / (y%d0 ** 2)
    s%d2 =   2 * x * y%d1 ** 2 / (y%d0 ** 3) - x * y%d2 / (y%d0 ** 2)
    s%d3 = -     x * y%d3        / (y%d0 ** 2) &
           + 6 * x * y%d1 * y%d2 / (y%d0 ** 3) &
           - 6 * x * y%d1 ** 3   / (y%d0 ** 4)
  end function nwad_dble_divy
  !>
  !> \brief Evaluate the exponentiation operator
  !>
  !> The implementation of the exponentiation operator. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = x(r)^{y(r)}\f$:
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=&  x(r)^{y(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=&  x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>     +\frac{y(r)}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right) \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=&
  !>      x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>        +\frac{2}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                       \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        +\frac{y(r)}{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>        -\frac{y(r)}{x(r)^2}\left(
  !>             \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2
  !>      \right) \\\\
  !>  &+& x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        +\frac{y(r)}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>      \right)^2 \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=&
  !>      x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3}
  !>        +\frac{3}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                       \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>        +\frac{3}{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>                       \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        +\frac{y(r)}{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}
  !>      \right. \\\\
  !>   && \left.
  !>     -\frac{3}{x(r)^2}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2
  !>                               \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        -\frac{3y(r)}{x(r)^2}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                             \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>        +\frac{2y(r)}{x(r)^3}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                             \right)^3
  !>      \right) \\\\
  !>  &+& 3x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        + \frac{y(r)}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>      \right) \\\\
  !>   && \left(\log(x(r))\frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2}
  !>        + \frac{2}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                        \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        + \frac{y(r)}{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>        - \frac{y(r)}{x(r)^2}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                             \right)^2
  !>      \right) \\\\
  !>  &+& x(r)^{y(r)}
  !>      \left(\log(x(r))\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>        + \frac{y(r)}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>      \right)^3
  !> \f}
  !>
  function nwad_dble_pow(x,y) result (s)
    type(nwad_dble), intent(in) :: x, y
    type(nwad_dble)             :: s
    s%d0 = x%d0**y%d0
    s%d1 = x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)
    s%d2 = x%d0**y%d0*(log(x%d0)*y%d2 + (2.0d0/x%d0)*x%d1*y%d1 &
                      +(y%d0/x%d0)*x%d2 - (y%d0/x%d0**2)*x%d1**2) &
         + x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)**2
    s%d3 = x%d0**y%d0* &
           (log(x%d0)*y%d3 + (3.0d0/x%d0)*x%d1*y%d2 &
            +(3.0d0/x%d0)*x%d2*y%d1 + (y%d0/x%d0)*x%d3 &
            -(3.0d0/x%d0**2)*x%d1**2*y%d1 - (3*y%d0/x%d0**2)*x%d1*x%d2 &
            +(2*y%d0/x%d0**3)*x%d1**3) &
         + 3*x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1) * &
           (log(x%d0)*y%d2 + (2.0d0/x%d0)*x%d1*y%d1 &
            +(y%d0/x%d0)*x%d2 - (y%d0/x%d0**2)*x%d1**2) &
         + x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)**3
  end function nwad_dble_pow
  !>
  !> \brief Evaluate the exponentiation operator where \f$y\f$ is inactive
  !> 
  !> We consider \f$s(r) = x(r)^y\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>       \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      y \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-1}
  !>              \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      y (y-1) \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-2}
  !>              \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2
  !>    + y       \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-1}
  !>                    \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>    \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      y(y-1)(y-2) \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-3}
  !>                  \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>    +3y(y-1)      \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-2}
  !>                        \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                        \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>    + y           \left(\frac{\mathrm{d}^0x(r)}{\mathrm{d}r^0}\right)^{y-1}
  !>                        \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}
  !> \f}
  !>
  function nwad_dble_powx(x,y) result (s)
    type(nwad_dble), intent(in)  :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s%d0 = x%d0**y
    s%d1 = x%d1*y*x%d0**(y - 1.0d0)
    s%d2 = x%d1**2 * y*(y-1.0d0)*x%d0**(y - 2.0d0) + &
           x%d2    * y*x%d0**(y - 1.0d0)
    s%d3 = x%d1**3    * y*(y-1.0d0)*(y-2.0d0)*x%d0**(y - 3.0d0) + &
           x%d1*x%d2  * 3*y*(y-1.0d0)*x%d0**(y - 2.0d0) + &
           x%d3       * y*x%d0**(y - 1.0d0)
  end function nwad_dble_powx
  !>
  !> \brief Evaluate the exponentiation operator where \f$x\f$ is inactive
  !>
  !> We consider \f$s(r) = x^{y(r)}\f$ then the derivatives are:
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& x^{y(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=&
  !>      x^{y(r)}\log(x) \frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=&
  !>      x^{y(r)}\left(\log(x) \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} +
  !>              \log^2(x)\left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^2
  !>              \right) \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=&
  !>      x^{y(r)}\left(
  !>      \log(x)\frac{\mathrm{d}^3y(r)}{\mathrm{d}r^3} +
  !>      3\log^2(x)\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}
  !>                \frac{\mathrm{d}^2y(r)}{\mathrm{d}r^2} + 
  !>      \log^3(x)\left(\frac{\mathrm{d}^1y(r)}{\mathrm{d}r^1}\right)^3
  !>      \right) \\\\
  !> \f}
  !>
  function nwad_dble_powy(x,y) result (s)
    double precision, intent(in) :: x
    type(nwad_dble), intent(in)  :: y
    type(nwad_dble)              :: s
    s%d0 = x**y%d0
    s%d1 = x**y%d0*log(x)*y%d1
    s%d2 = x**y%d0*(log(x)*y%d2 + (log(x)*y%d1)**2)
    s%d3 = x**y%d0*(log(x)*y%d3 + 3*log(x)**2*y%d1*y%d2+(log(x)*y%d1)**3)
  end function nwad_dble_powy
  !>
  !> \brief Evaluate the \f$|\;\;|\f$ function
  !>
  !> The implementation of the \f$|\;\;|\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = |x(r)|\f$.
  !> Note that this function is continuous but not continously differentiable.
  !> The resolution is that at \f$x(r) = 0\f$ we assume all derivatives to be
  !> zero as well. The derivatives are:
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& |x(r)| \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       \frac{x(r)}{|x(r)|}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       \frac{x(r)}{|x(r)|}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>       \frac{x(r)}{|x(r)|}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} \\\\
  !> \f}
  !>
  function nwad_dble_abs(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = abs(x%d0)
    if (x%d0.lt.0.0d0) then
      s%d1 = -x%d1
      s%d2 = -x%d2
      s%d3 = -x%d3
    else if (x%d0.gt.0.0d0) then
      s%d1 =  x%d1
      s%d2 =  x%d2
      s%d3 =  x%d3
    else
      s%d1 =  0.0d0
      s%d2 =  0.0d0
      s%d3 =  0.0d0
    endif
  end function nwad_dble_abs
  !>
  !> \brief Evaluate the \f$\sqrt{\;\;\;}\f$ function
  !>
  !> The implementation of the \f$\sqrt{\;\;\;}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider 
  !> \f$s(r) = \sqrt{x(r)}\f$:
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \sqrt{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>       \frac{1}{2\sqrt{x(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>       \frac{1}{2\sqrt{x(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     - \frac{1}{4x(r)\sqrt{x(r)}}
  !>       \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>       \frac{1}{2\sqrt{x(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}
  !>     - 3\frac{1}{4x(r)\sqrt{x(r)}}
  !>       \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>       \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     + \frac{3}{8x^2(r)\sqrt{x(r)}}
  !>       \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_sqrt(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = sqrt(x%d0)
    s%d1 = 0.5d0/sqrt(x%d0)*x%d1
    s%d2 = 0.5d0/sqrt(x%d0)*x%d2-0.25d0/(x%d0*sqrt(x%d0))*x%d1**2
    s%d3 = 0.5d0/sqrt(x%d0)*x%d3-0.75d0/(x%d0*sqrt(x%d0))*x%d1*x%d2 &
         + 3.0d0/(8.0d0*x%d0**2*sqrt(x%d0))*x%d1**3
  end function nwad_dble_sqrt
  !>
  !> \brief Evaluate the \f$\exp\f$ function
  !>
  !> The implementation of the \f$\exp\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \exp(x(r))\f$.
  !>
  function nwad_dble_exp(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = x%d0*exp(x%d0)
    s%d1 = x%d1*exp(x%d0)
    s%d2 = x%d2*exp(x%d0) + x%d1**2*exp(x%d0)
    s%d3 = x%d3*exp(x%d0) + 3*x%d1*x%d2*exp(x%d0) + x%d1**3*exp(x%d0)
  end function nwad_dble_exp
  !>
  !> \brief Evaluate the \f$\log\f$ function
  !>
  !> The implementation of the \f$\log\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \log(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \log{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     -\left(\frac{1}{x(r)}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3}
  !>     -\frac{3}{x^2(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                      \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +2\left(\frac{1}{x(r)}
  !>             \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_log(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = log(x%d0)
    s%d1 = x%d1/x%d0
    s%d2 = x%d2/x%d0 - (x%d1/x%d0)**2
    s%d3 = x%d3/x%d0 - 3*(x%d1/x%d0)*(x%d2/x%d0) + 2*(x%d1/x%d0)**3
  end function nwad_dble_log
  !>
  !> \brief Evaluate the \f$\sin\f$ function
  !>
  !> The implementation of the \f$\sin\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \sin(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \sin{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \cos{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \cos{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\sin{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \cos{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -3\sin{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                 \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     -\cos{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_sin(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = sin(x%d0)
    s%d1 = x%d1*cos(x%d0)
    s%d2 = x%d2*cos(x%d0) - x%d1**2*sin(x%d0)
    s%d3 = x%d3*cos(x%d0) - 3*x%d1*x%d2*sin(x%d0) - x%d1**3*cos(x%d0)
  end function nwad_dble_sin
  !>
  !> \brief Evaluate the \f$\cos\f$ function
  !>
  !> The implementation of the \f$\cos\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \cos(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \cos{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>     -\sin{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>     -\sin{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\cos{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>     -\sin{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -3\cos{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                 \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\sin{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_cos(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = cos(x%d0)
    s%d1 = -x%d1*sin(x%d0)
    s%d2 = -x%d2*sin(x%d0) - x%d1**2*cos(x%d0)
    s%d3 = -x%d3*sin(x%d0) - 3*x%d1*x%d2*cos(x%d0) + x%d1**3*sin(x%d0)
  end function nwad_dble_cos
  !>
  !> \brief Evaluate the \f$\tan\f$ function
  !>
  !> The implementation of the \f$\tan\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \tan(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \tan{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{\cos^2{x(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{\cos^2{x(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     +\frac{2\tan{x(r)}}{\cos^2{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{\cos^2{x(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     +\frac{6\tan{x(r)}}{\cos^2{x(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                                      \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\frac{4\tan^2{x(r)}}{\cos^2{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     +\frac{2}{\cos^4{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_tan(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = tan(x%d0)
    s%d1 = x%d1/cos(x%d0)**2
    s%d2 = x%d2/cos(x%d0)**2 + x%d1**2*(2*tan(x%d0)/cos(x%d0)**2)
    s%d3 = x%d3/cos(x%d0)**2 + 6*x%d1*x%d2*tan(x%d0)/cos(x%d0)**2 &
         + 4*x%d1**3*tan(x%d0)**2/cos(x%d0)**2 + 2*x%d1**3/cos(x%d0)**4
  end function nwad_dble_tan
  !>
  !> \brief Evaluate the \f$\sinh\f$ function
  !>
  !> The implementation of the \f$\sinh\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \sinh(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \sinh{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \cosh{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \cosh{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     +\sinh{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \cosh{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     +3\sinh{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                  \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\cosh{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_sinh(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = sinh(x%d0)
    s%d1 = x%d1*cosh(x%d0)
    s%d2 = x%d2*cosh(x%d0) + x%d1**2*sinh(x%d0)
    s%d3 = x%d3*cosh(x%d0) + 3*x%d1*x%d2*sinh(x%d0) + x%d1**3*cosh(x%d0)
  end function nwad_dble_sinh
  !>
  !> \brief Evaluate the \f$\cosh\f$ function
  !>
  !> The implementation of the \f$\cosh\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \cosh(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \cosh{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \sinh{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \sinh{x(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     +\cosh{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \sinh{x(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     +3\cosh{x(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                  \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\sinh{x(r)}\left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_cosh(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = cosh(x%d0)
    s%d1 = x%d1*sinh(x%d0)
    s%d2 = x%d2*sinh(x%d0) + x%d1**2*cosh(x%d0)
    s%d3 = x%d3*sinh(x%d0) + 3*x%d1*x%d2*cosh(x%d0) + x%d1**3*sinh(x%d0)
  end function nwad_dble_cosh
  !>
  !> \brief Evaluate the \f$\tanh\f$ function
  !>
  !> The implementation of the \f$\tanh\f$ function. The chain rule is used
  !> to evaluate the derivatives. I.e. we consider \f$s(r) = \tanh(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \tanh{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{\cosh^2{x(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{\cosh^2{x(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\frac{2\tanh{x(r)}}{\cosh^2{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{\cosh^2{x(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -\frac{6\tanh{x(r)}}{\cosh^2{x(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>                                      \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\frac{4\tanh^2{x(r)}}{\cosh^2{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     -\frac{2}{\cosh^4{x(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  function nwad_dble_tanh(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = tanh(x%d0)
    s%d1 = x%d1/cosh(x%d0)**2
    s%d2 = x%d2/cosh(x%d0)**2 - x%d1**2*(2*tanh(x%d0)/cosh(x%d0)**2)
    s%d3 = x%d3/cosh(x%d0)**2 - 6*x%d1*x%d2*tanh(x%d0)/cosh(x%d0)**2 &
         + 4*x%d1**3*tanh(x%d0)**2/cosh(x%d0)**2 - 2*x%d1**3/cosh(x%d0)**4
  end function nwad_dble_tanh
  !>
  !> \brief Evaluate the \f$\mathrm{asin}\f$ function
  !>
  !> The implementation of the \f$\mathrm{asin}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{asin}(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \asin{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     +\frac{x(r)}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     +\frac{3x(r)}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>  && +\frac{1}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     +\frac{3x^2(r)}{\left(1-x^2(r)\right)^2\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_asin(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1, t12
    t1 = 1.0d0 - x%d0*x%d0
    t12 = sqrt(t1)
    s%d0 = asin(x%d0)
    s%d1 = x%d1/t12
    s%d2 = x%d2/t12 + x%d1**2*x%d0/(t1*t12)
    s%d3 = x%d3/t12 + 3*x%d1*x%d2*x%d0/(t1*t12) + x%d1**3/(t1*t12) &
         + 3*x%d1**3*x%d0**2/(t1*t1*t12)
  end function nwad_dble_asin
  !>
  !> \brief Evaluate the \f$\mathrm{acos}\f$ function
  !>
  !> The implementation of the \f$\mathrm{acos}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{acos}(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \acos{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>     -\frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>     -\frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\frac{x(r)}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>     -\frac{1}{\sqrt{1-x^2(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -\frac{3x(r)}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>  && -\frac{1}{\left(1-x^2(r)\right)\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     -\frac{3x^2(r)}{\left(1-x^2(r)\right)^2\sqrt{1-x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_acos(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1, t12
    t1 = 1.0d0 - x%d0*x%d0
    t12 = sqrt(t1)
    s%d0 = acos(x%d0)
    s%d1 = -x%d1/t12
    s%d2 = -x%d2/t12 - x%d1**2*x%d0/(t1*t12)
    s%d3 = -x%d3/t12 - 3*x%d1*x%d2*x%d0/(t1*t12) - x%d1**3/(t1*t12) &
         - 3*x%d1**3*x%d0**2/(t1*t1*t12)
  end function nwad_dble_acos
  !>
  !> \brief Evaluate the \f$\mathrm{atan}\f$ function
  !>
  !> The implementation of the \f$\mathrm{atan}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{atan}(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \atan{x(r)} \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{1+x^2(r)}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{1+x^2(r)}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\frac{2x(r)}{\left(1+x^2(r)\right)^2}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{1+x^2(r)}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -\frac{6x(r)}{\left(1+x^2(r)\right)^2}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>  && -\frac{2}{\left(1+x^2(r)\right)^2}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     +\frac{8x^2(r)}{\left(1+x^2(r)\right)^3}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_atan(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1
    t1 = 1.0d0 + x%d0*x%d0
    s%d0 = atan(x%d0)
    s%d1 = x%d1/t1
    s%d2 = x%d2/t1 - x%d1**2*x%d0/(t1*t1)
    s%d3 = x%d3/t1 - 6*x%d1*x%d2*x%d0/(t1*t1) - 2*x%d1**3/(t1*t1) &
         + 8*x%d1**3*x%d0**2/(t1*t1*t1)
  end function nwad_dble_atan
  !>
  !> \brief Evaluate the \f$\mathrm{asinh}\f$ function
  !>
  !> The implementation of the \f$\mathrm{asinh}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{asinh}(x(r))\f$,
  !> where \f$\mathrm{asinh}(a) = \log(a+\sqrt{1+a^2})\f$
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>      \log(x(r)+\sqrt{1+x^2(r)}) \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{1}{\sqrt{1+x^2(r)}}\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{1}{\sqrt{1+x^2(r)}}\frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\frac{x(r)}{\left(1+x^2(r)\right)\sqrt{1+x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{1}{\sqrt{1+x^2(r)}}\frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -\frac{3x(r)}{\left(1+x^2(r)\right)\sqrt{1+x^2(r)}}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} \\\\
  !>  && -\frac{1}{\left(1+x^2(r)\right)\sqrt{1+x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     +\frac{3x^2(r)}{\left(1+x^2(r)\right)^2\sqrt{1+x^2(r)}}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_asinh(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1, t12
    t1 = 1.0d0 + x%d0*x%d0
    t12 = sqrt(t1)
    s%d0 = log(x%d0+t12)
    s%d1 = x%d1/t12
    s%d2 = x%d2/t12 - x%d1**2*x%d0/(t1*t12)
    s%d3 = x%d3/t12 - 3*x%d1*x%d2*x%d0/(t1*t12) - x%d1**3/(t1*t12) &
         + 3*x%d1**3*x%d0**2/(t1*t1*t12)
  end function nwad_dble_asinh
  !>
  !> \brief Evaluate the \f$\mathrm{erf}\f$ function
  !>
  !> The implementation of the \f$\mathrm{erf}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{erf}(x(r))\f$, where
  !> \f$\mathrm{erf}(a)=\frac{2}{\sqrt{\pi}}\int_0^a e^{-t^2}\mathrm{d}t\f$
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>      \mathrm{erf}(x(r)) \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>      \frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>      \frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     -\frac{4x(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>      \frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     -\frac{12x(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     +\frac{8x^2(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     -\frac{4}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_erf(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1
    t1 = exp(-(x%d0**2))/sqrt(acos(-1.0d0))
    s%d0 = erf(x%d0)
    s%d1 = 2*t1*x%d1
    s%d2 = 2*t1*x%d2 - 4*x%d0*t1*x%d1**2
    s%d3 = 2*t1*x%d3 - 12*x%d0*t1*x%d1*x%d2 + 8*x%d0**2*t1*x%d1**3 &
         - 4*t1*x%d1**3
  end function nwad_dble_erf
  !>
  !> \brief Evaluate the \f$\mathrm{erfc}\f$ function
  !>
  !> The implementation of the \f$\mathrm{erfc}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{erfc}(x(r))\f$, where
  !> \f$\mathrm{erfc}(a)=1-\mathrm{erf}(a)\f$
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& 
  !>      1-\mathrm{erf}(x(r)) \\\\
  !>   \frac{\mathrm{d}^1s(r)}{\mathrm{d}r^1} &=& 
  !>     -\frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1} \\\\
  !>   \frac{\mathrm{d}^2s(r)}{\mathrm{d}r^2} &=& 
  !>     -\frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2} 
  !>     +\frac{4x(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^2 \\\\
  !>   \frac{\mathrm{d}^3s(r)}{\mathrm{d}r^3} &=& 
  !>     -\frac{2}{\sqrt{\pi}}e^{-x^2(r)}
  !>      \frac{\mathrm{d}^3x(r)}{\mathrm{d}r^3} 
  !>     +\frac{12x(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}
  !>            \frac{\mathrm{d}^2x(r)}{\mathrm{d}r^2}
  !>     -\frac{8x^2(r)}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3
  !>     +\frac{4}{\sqrt{\pi}}e^{-x^2(r)}
  !>            \left(\frac{\mathrm{d}^1x(r)}{\mathrm{d}r^1}\right)^3 \\\\
  !> \f}
  !>
  function nwad_dble_erfc(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    double precision             :: t1
    t1 = exp(-(x%d0**2))/sqrt(acos(-1.0d0))
    s%d0 = 1.0d0-erf(x%d0)
    s%d1 = -2*t1*x%d1
    s%d2 = -2*t1*x%d2 + 4*x%d0*t1*x%d1**2
    s%d3 = -2*t1*x%d3 + 12*x%d0*t1*x%d1*x%d2 - 8*x%d0**2*t1*x%d1**3 &
         + 4*t1*x%d1**3
  end function nwad_dble_erfc
  !>
  !> \brief Initialize an inactive variable
  !>
  !> Initialize an inactive variable. Inactive variables are those with respect
  !> to which no derivatives are calculated in the current evaluation of the
  !> code. In practice it means that the components are initialized as
  !> \f{eqnarray*}{
  !>   d0 &=& x \\\\
  !>   d1 &=& 0 \\\\
  !>   d2 &=& 0 \\\\
  !>   d3 &=& 0 \\\\
  !> \f}
  !>
  function nwad_dble_inactive(x) result (s)
    double precision, intent(in) :: x
    type(nwad_dble)              :: s
    s%d0 = x
    s%d1 = 0
    s%d2 = 0
    s%d3 = 0
  end function nwad_dble_inactive
  !>
  !> \brief Initialize an active variable
  !>
  !> Initialize an active variable. Active variables are those with respect
  !> to which the derivatives are calculated in the current evaluation of the
  !> code. In practice it means that the components are initialized as
  !> \f{eqnarray*}{
  !>   d0 &=& \frac{\mathrm{d}^0 x}{\mathrm{d}x^0} = x \\\\
  !>   d1 &=& \frac{\mathrm{d}^1 x}{\mathrm{d}x^1} = 1 \\\\
  !>   d2 &=& \frac{\mathrm{d}^2 x}{\mathrm{d}x^2} = 0 \\\\
  !>   d3 &=& \frac{\mathrm{d}^3 x}{\mathrm{d}x^3} = 0 \\\\
  !> \f}
  !>
  function nwad_dble_active(x) result (s)
    double precision, intent(in) :: x
    type(nwad_dble)              :: s
    s%d0 = x
    s%d1 = 1
    s%d2 = 0
    s%d3 = 0
  end function nwad_dble_active
  !>
  !> \brief Initialize an negatively active variable
  !>
  !> Initialize an negatively active variable. Active variables are those with
  !> respect to which the derivatives are calculated in the current evaluation
  !> of the code. In practice it means that the components are initialized as
  !> \f{eqnarray*}{
  !>   d0 &=& \frac{\mathrm{d}^0 x}{\mathrm{d}(-x)^0} =  x \\\\
  !>   d1 &=& \frac{\mathrm{d}^1 x}{\mathrm{d}(-x)^1} = -1 \\\\
  !>   d2 &=& \frac{\mathrm{d}^2 x}{\mathrm{d}(-x)^2} =  0 \\\\
  !>   d3 &=& \frac{\mathrm{d}^3 x}{\mathrm{d}(-x)^3} =  0 \\\\
  !> \f}
  !>
  function nwad_dble_active_neg(x) result (s)
    double precision, intent(in) :: x
    type(nwad_dble)              :: s
    s%d0 =  x
    s%d1 = -1
    s%d2 =  0
    s%d3 =  0
  end function nwad_dble_active_neg
end module nwad
!> @}
