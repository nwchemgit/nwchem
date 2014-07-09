!> @defgroup nwad NW Automatic Differentiation
!>   @ingroup nwxc_priv
!>
!> \ingroup nwad
!> @{
!>
!> \brief NW Automatic Differentation
!>
!> # NW Automatic Differentiation #
!>
!> A module implementing Automatic Differentiation [1,2] capabilities in 
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
!> Every function execution computes only one specific partial derivative [3].
!>
!> ## The chain rule and Taylor series ##
!>
!> The chain rule is fundamental to automatic differentiation. The way 
!> derivatives can be described using the chain rule separates into two
!> classes. In one class the chain rule is directly applied to expressions
!> operating on functions. In the other class the functions are first written
!> as Taylor series, and the chain rule is applied to expressions on Taylor
!> series. In this section the two approaches are compared.
!>
!> In general we have to consider two kinds of expressions. One kind involves
!> unary operators and functions of one argument. Examples are: the negation
!> operator \f$-\f$, the factorial operator \f$!\f$, the trigoniometric
!> functions \f$\sin\f$, \f$\cos\f$, \f$\tan\f$, etc. The other kind of
!> expressions involves binary operators, i.e. the addition \f$+\f$,
!> subtraction \f$-\f$, division \f$/\f$, multiplication \f$*\f$, and 
!> exponentiation \f$\fbox{}^\fbox{}\f$. These kinds of expressions can be
!> summarized in the equations
!> \f{eqnarray*}{
!>   z &=& F(x)
!> \f}
!> for the first kind and
!> \f{eqnarray*}{
!>   z &=& F(x,y)
!> \f}
!> for the second.
!>
!> To derive expressions for the derivatives of \f$z\f$ in terms of derivatives
!> of \f$x\f$ and \f$y\f$ using the chain rule first consider these quantities
!> formally as functions \f$z(t), x(t)\f$ and \f$y(t)\f$. Applying the chain
!> rule up to 3rd order one obtains
!> \f{eqnarray*}{
!>   \frac{\partial^0 z(t)}{\partial t^0}
!>   &=& \frac{\partial^0 F(x(t))}{\partial t^0} \\\\
!>   \frac{\partial^1 z(t)}{\partial t^1}
!>   &=& \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^1 x(t)}{\partial t^1} \\\\
!>   \frac{\partial^2 z(t)}{\partial t^2}
!>   &=& \frac{\partial^2 F(x(t))}{\partial x^2}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^2
!>    +  \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^2 x(t)}{\partial t^2} \\\\
!>   \frac{\partial^3 z(t)}{\partial t^3}
!>   &=& \frac{\partial^3 F(x(t))}{\partial x^3}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^3
!>    + 3\frac{\partial^2 F(x(t))}{\partial x^2}
!>       \frac{\partial^2 x(t)}{\partial t^2}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>    +  \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^3 x(t)}{\partial t^3} \\\\
!> \f}
!> and
!> \f{eqnarray*}{
!>   \frac{\partial^0 z(t)}{\partial t^0}
!>   &=& \frac{\partial^0 F(x(t),y(t))}{\partial t^0} \\\\
!>   \frac{\partial^1 z(t)}{\partial t^1}
!>   &=& \frac{\partial^1 F(x(t),y(t))}{\partial x^1}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>    +  \frac{\partial^1 F(x(t),y(t))}{\partial y^1}
!>       \frac{\partial^1 y(t)}{\partial t^1} \\\\
!>   \frac{\partial^2 z(t)}{\partial t^2}
!>   &=& \frac{\partial^2 F(x(t),y(t))}{\partial x^2}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^2
!>    + 2\frac{\partial^2 F(x(t),y(t))}{\partial x^1\partial y^1}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>       \frac{\partial^1 y(t)}{\partial t^1} 
!>    +  \frac{\partial^2 F(x(t),y(t))}{\partial y^2}
!>       \left(\frac{\partial^1 y(t)}{\partial t^1}\right)^2 \\\\
!>   &+& \frac{\partial^1 F(x(t),y(t))}{\partial x^1}
!>       \frac{\partial^2 x(t)}{\partial t^2}
!>    +  \frac{\partial^1 F(x(t),y(t))}{\partial y^1}
!>       \frac{\partial^2 y(t)}{\partial t^2} \\\\
!>   \frac{\partial^3 z(t)}{\partial t^3}
!>   &=& \frac{\partial^3 F(x(t),y(t))}{\partial x^3}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^3
!>    + 3\frac{\partial^3 F(x(t),y(t))}{\partial x^2\partial y}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^2
!>       \frac{\partial^1 y(t)}{\partial t^1}
!>    + 3\frac{\partial^3 F(x(t),y(t))}{\partial x\partial y^2}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>       \left(\frac{\partial^1 y(t)}{\partial t^1}\right)^2
!>    +  \frac{\partial^3 F(x(t),y(t))}{\partial y^3}
!>       \left(\frac{\partial^1 y(t)}{\partial t^1}\right)^3 \\\\
!>   &+&3\frac{\partial^2 F(x(t),y(t))}{\partial x^2}
!>       \frac{\partial^2 x(t)}{\partial t^2}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>    + 3\frac{\partial^2 F(x(t),y(t))}{\partial x^1\partial y^1}
!>       \frac{\partial^2 x(t)}{\partial t^2}
!>       \frac{\partial^1 y(t)}{\partial t^1}
!>    + 3\frac{\partial^2 F(x(t),y(t))}{\partial x^1\partial y^1}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>       \frac{\partial^2 y(t)}{\partial t^2}
!>    + 3\frac{\partial^2 F(x(t),y(t))}{\partial y^2}
!>       \frac{\partial^2 y(t)}{\partial t^2}
!>       \frac{\partial^1 y(t)}{\partial t^1}  \\\\
!>   &+& \frac{\partial^1 F(x(t),y(t))}{\partial x^1}
!>       \frac{\partial^3 x(t)}{\partial t^3}
!>    +  \frac{\partial^1 F(x(t),y(t))}{\partial y^1}
!>       \frac{\partial^3 y(t)}{\partial t^3}
!> \f}
!> Alternatively one can consider the functions \f$z(t), x(t)\f$ and \f$y(t)\f$
!> in terms of their Taylor expansions up to degree \f$d\f$ (where \f$d=3\f$
!> here) [4] (see page 147):
!> \f{eqnarray*}{
!>    z(t) &=& \sum_{j=0}^{d} z_j t^j \\\\
!>    x(t) &=& \sum_{j=0}^{d} x_j t^j \\\\
!>    y(t) &=& \sum_{j=0}^{d} y_j t^j \\\\
!> \f}
!> where \f$z_j, x_j, y_j \in \mathbb{R}^n\f$ and \f$t \in \mathbb{R}\f$.
!> Thus \f$z, x\f$ and \f$y\f$ are vector polynomials in the scalar variable
!> \f$t\f$. The Taylor coefficients are given by
!> \f{eqnarray*}{
!>    x_j &=& \left.\frac{1}{j!}\frac{\partial^j}{\partial t^j}x(t)\right|_{t=0}
!> \f}
!> Now the problem \f$z(t) = F(x(t))\f$ can be approached by substituting the
!> Taylor expansions, differentiating the equation with respect to \f$t\f$ and
!> evaluating the resulting expression at \f$t=0\f$ to get the Taylor
!> coefficients of \f$z\f$ in terms of those of \f$x\f$. From this we have
!> \f{eqnarray*}{
!>   z_0 &=& F(x_0) \\\\
!>   z_1 &=& \frac{\partial F(x_0)}{\partial x}x_1  \\\\
!>  2z_2 &=& \frac{\partial^2 F(x_0)}{\partial x^2}x_1^2 
!>        + 2\frac{\partial F(x_0)}{\partial x}x_2  \\\\
!>  6z_3 &=& \frac{\partial^3 F(x_0)}{\partial x^3}x_1^3
!>        + 6\frac{\partial^2 F(x_0)}{\partial x^2}x_2 x_1
!>        + 6\frac{\partial   F(x_0)}{\partial x}x_3
!> \f}
!> Similarly the problem \f$z(t) = F(x(t),y(t))\f$ can approached to yield
!> \f{eqnarray*}{
!>   z_0 &=& F(x_0,y_0) \\\\
!>   z_1 &=& \frac{\partial F(x_0,y_0)}{\partial x}x_1
!>        +  \frac{\partial F(x_0,y_0)}{\partial y}y_1 \\\\
!>  2z_2 &=& \frac{\partial^2 F(x_0,y_0)}{\partial x^2}x_1^2
!>        + 2\frac{\partial^2 F(x_0,y_0)}{\partial x \partial y}x_1 y_1
!>        +  \frac{\partial^2 F(x_0,y_0)}{\partial y^2}y_1^2 \\\\
!>       &+&2\frac{\partial F(x_0,y_0)}{\partial x}x_2
!>        + 2\frac{\partial F(x_0,y_0)}{\partial y}y_2 \\\\
!>  6z_3 &=& \frac{\partial^3 F(x_0,y_0)}{\partial x^3}x_1^3
!>        + 3\frac{\partial^3 F(x_0,y_0)}{\partial x^2 \partial y}x_1^2 y_1
!>        + 3\frac{\partial^3 F(x_0,y_0)}{\partial x \partial y^2}x_1 y_1^2
!>        +  \frac{\partial^3 F(x_0,y_0)}{\partial y^3}y_1^3 \\\\
!>       &+&6\frac{\partial^2 F(x_0,y_0)}{\partial x^2}x_2 x_1 
!>        + 6\frac{\partial^2 F(x_0,y_0)}{\partial x \partial y}x_2 y_1
!>        + 6\frac{\partial^2 F(x_0,y_0)}{\partial x \partial y}x_1 y_2
!>        + 6\frac{\partial^2 F(x_0,y_0)}{\partial y^2}y_2 y_1  \\\\
!>       &+&6\frac{\partial   F(x_0,y_0)}{\partial x}x_3
!>        + 6\frac{\partial   F(x_0,y_0)}{\partial y}y_3
!> \f}
!> Substituting the expressions for the Taylor coefficients we get for 
!> \f$z(t) = F(x(t))\f$ case (the \f$z(t) = F(x(t),y(t))\f$ case is left as 
!> an excercise for the reader)
!> \f{eqnarray*}{
!>   \frac{\partial^0 z(t)}{\partial t^0}
!>   &=& \frac{\partial^0 F(x(t))}{\partial t^0} \\\\
!>   \frac{\partial^1 z(t)}{\partial t^1}
!>   &=& \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^1 x(t)}{\partial t^1} \\\\
!>   \frac{\partial^2 z(t)}{\partial t^2}
!>   &=& \frac{\partial^2 F(x(t))}{\partial x^2}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^2
!>    +  \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^2 x(t)}{\partial t^2} \\\\
!>   \frac{\partial^3 z(t)}{\partial t^3}
!>   &=& \frac{\partial^3 F(x(t))}{\partial x^3}
!>       \left(\frac{\partial^1 x(t)}{\partial t^1}\right)^3
!>    + 3\frac{\partial^2 F(x(t))}{\partial x^2}
!>       \frac{\partial^2 x(t)}{\partial t^2}
!>       \frac{\partial^1 x(t)}{\partial t^1}
!>    +  \frac{\partial^1 F(x(t))}{\partial x^1}
!>       \frac{\partial^3 x(t)}{\partial t^3} \\\\
!> \f}
!> Hence we find that applying the chain rule directly or expressing the
!> problem in terms of Taylor series leads to exactly the same results overall.
!> 
!> The advantage of Taylor series is that the formal properties of polynomials
!> can be used as has been done by Griewank et al. [3].
!>
!> ### References ###
!>
!> [1] R. E. Wengert (1964) "A simple automatic derivative evaluation program",
!>     Communications of the ACM, <b>7</b>, pp. 463-464, DOI:
!>     <a href="http://dx.doi.org/10.1145/355586.364791">
!>     10.1145/355586.364791</a>.
!>
!> [2] See e.g. <a href="http://www.autodiff.org/">www.autodiff.org</a>
!>
!> [3] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
!>     tensors by forward propagation of univariate Taylor series",
!>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
!>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
!>     10.1090/S0025-5718-00-01120-0</a>.
!>
!> [4] A. Griewank, D. Juedes, J. Utke (1996) "Algorithm 755: ADOL-C: A package
!>     for the automatic differentiation of algorithms written in C/C++",
!>     ACM Transactions on Mathematical Software, <b>22</b>, pp. 131-167, DOI:
!>     <a href="http://dx.doi.org/10.1145/229473.229474">
!>     10.1145/229473.229474</a>.
!>
!> [5] I. Charpentier, J. Utke, "Rapsodia: User Manual", Argonne National
!>     Laboratory, <a href="http://www.mcs.anl.gov/Rapsodia/userManual.pdf">
!>     http://wwww.mcs.anl.gov/Rapsodia/userManual.pdf</a> (referenced
!>     July 3, 2014).
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
  interface max
    module procedure nwad_dble_max
  end interface
  interface min
    module procedure nwad_dble_min
  end interface
  interface operator (+)
    module procedure nwad_dble_add
    module procedure nwad_dble_addx
    module procedure nwad_dble_addy
  end interface
  interface operator (-)
    module procedure nwad_dble_minus
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
  interface operator (.eq.)
    module procedure nwad_dble_equal
    module procedure nwad_dble_equalx
    module procedure nwad_dble_equaly
  end interface operator (==)
  interface operator (.ne.)
    module procedure nwad_dble_notequal
    module procedure nwad_dble_notequalx
    module procedure nwad_dble_notequaly
  end interface operator (/=)
  interface operator (.lt.)
    module procedure nwad_dble_lessthan
    module procedure nwad_dble_lessthanx
    module procedure nwad_dble_lessthany
  end interface operator (<)
  interface operator (.le.)
    module procedure nwad_dble_lessequal
    module procedure nwad_dble_lessequalx
    module procedure nwad_dble_lessequaly
  end interface operator (<=)
  interface operator (.gt.)
    module procedure nwad_dble_greaterthan
    module procedure nwad_dble_greaterthanx
    module procedure nwad_dble_greaterthany
  end interface operator (>)
  interface operator (.ge.)
    module procedure nwad_dble_greaterequal
    module procedure nwad_dble_greaterequalx
    module procedure nwad_dble_greaterequaly
  end interface operator (>=)
  interface sign
    module procedure nwad_dble_sign
    module procedure nwad_dble_signx
    module procedure nwad_dble_signy
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
    module procedure nwad_dble_active_n
  end interface
  interface active_neg
    module procedure nwad_dble_active_neg
  end interface
  interface inactive
    module procedure nwad_dble_inactive
  end interface
  interface value
    module procedure nwad_dble_value
  end interface
  interface inter_d1_dx
    module procedure nwad_dble_inter_d1_dx
  end interface
  interface inter_d2_dx
    module procedure nwad_dble_inter_d2_dx
  end interface
  interface inter_d2_dx2
    module procedure nwad_dble_inter_d2_dx2
  end interface
  interface inter_d2_dxy
    module procedure nwad_dble_inter_d2_dxy
  end interface
  interface inter_d3_dx
    module procedure nwad_dble_inter_d3_dx
  end interface
  interface inter_d3_dx2
    module procedure nwad_dble_inter_d3_dx2
  end interface
  interface inter_d3_dxy
    module procedure nwad_dble_inter_d3_dxy
  end interface
  interface inter_d3_dx3
    module procedure nwad_dble_inter_d3_dx3
  end interface
  interface inter_d3_dx2y
    module procedure nwad_dble_inter_d3_dx2y
  end interface
  interface inter_d3_dxyz
    module procedure nwad_dble_inter_d3_dxyz
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
  !> Find the maximum value of the arguments
  !>
  !> This routine finds the maximum value of all the presented arguments.
  !> In Fortran the MAX function must have at least 2 arguments but can
  !> have any arbitrary number of arguments. Also the arguments all have to be
  !> of the same type. Here this capability is implemented using optional
  !> arguments, allowing for a maximum of 5 arguments. Whereas Fortran 90 
  !> allows optional arguments to be specified by name, this implementation
  !> will produce incorrect results if an optional argument in the middle
  !> is left out!
  !>
  function nwad_dble_max(a,b,c,d,e) result (s)
    type(nwad_dble), intent(in)           :: a
    type(nwad_dble), intent(in)           :: b
    type(nwad_dble), intent(in), optional :: c
    type(nwad_dble), intent(in), optional :: d
    type(nwad_dble), intent(in), optional :: e
    type(nwad_dble)                       :: s
    type(nwad_dble)                       :: t1
    type(nwad_dble)                       :: t2
    if (.not.present(c)) then
      if (a%d0 .ge. b%d0 ) then
        s = a
      else
        s = b
      endif
    else
      if (a%d0 .ge. b%d0) then
        t1 = a
      else
        t1 = b
      endif
      if (.not.present(d)) then
        if (t1%d0 .ge. c%d0) then
          s = t1
        else
          s = c
        endif
      else
        if (t1%d0 .ge. c%d0) then
          t2 = t1
        else
          t2 = c
        endif
        if (.not.present(e)) then
          if (t2%d0 .ge. d%d0) then
            s = t2
          else
            s = d
          endif
        else
          if (t2%d0 .ge. d%d0) then
            t1 = t2
          else
            t1 = d
          endif
          if (t1%d0 .ge. e%d0) then
            s = t1
          else
            s = e
          endif
        endif
      endif
    endif
  end function nwad_dble_max
  !>
  !> Find the minimum value of the arguments
  !>
  !> This routine finds the minimum value of all the presented arguments.
  !> In Fortran the MIN function must have at least 2 arguments but can
  !> have any arbitrary number of arguments. Also the arguments all have to be
  !> of the same type. Here this capability is implemented using optional
  !> arguments, allowing for a maximum of 5 arguments. Whereas Fortran 90 
  !> allows optional arguments to be specified by name, this implementation
  !> will produce incorrect results if an optional argument in the middle
  !> is left out!
  !>
  function nwad_dble_min(a,b,c,d,e) result (s)
    type(nwad_dble), intent(in)           :: a
    type(nwad_dble), intent(in)           :: b
    type(nwad_dble), intent(in), optional :: c
    type(nwad_dble), intent(in), optional :: d
    type(nwad_dble), intent(in), optional :: e
    type(nwad_dble)                       :: s
    type(nwad_dble)                       :: t1
    type(nwad_dble)                       :: t2
    if (.not.present(c)) then
      if (a%d0 .lt. b%d0 ) then
        s = a
      else
        s = b
      endif
    else
      if (a%d0 .lt. b%d0) then
        t1 = a
      else
        t1 = b
      endif
      if (.not.present(d)) then
        if (t1%d0 .lt. c%d0) then
          s = t1
        else
          s = c
        endif
      else
        if (t1%d0 .lt. c%d0) then
          t2 = t1
        else
          t2 = c
        endif
        if (.not.present(e)) then
          if (t2%d0 .lt. d%d0) then
            s = t2
          else
            s = d
          endif
        else
          if (t2%d0 .lt. d%d0) then
            t1 = t2
          else
            t1 = d
          endif
          if (t1%d0 .lt. e%d0) then
            s = t1
          else
            s = e
          endif
        endif
      endif
    endif
  end function nwad_dble_min
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
         - 3 * x%d1 * y%d2        / (y%d0 ** 2) &
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
  !> \brief Return whether \f$x\f$ equals \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
  !> context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_equal(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = x%d0 .eq. y%d0
  end function nwad_dble_equal
  !>
  !> \brief Return whether \f$x\f$ equals \f$y\f$ where the latter is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
  !> context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_equalx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = x%d0 .eq. y
  end function nwad_dble_equalx
  !>
  !> \brief Return whether \f$x\f$ equals \f$y\f$ where the former is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
  !> context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_equaly(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = x .eq. y%d0
  end function nwad_dble_equaly
  !>
  !> \brief Return whether \f$x\f$ does not equal \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_notequal(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = .not. (x .eq. y)
  end function nwad_dble_notequal
  !>
  !> \brief Return whether \f$x\f$ does not equal \f$y\f$ where the latter is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_notequalx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = .not. (x .eq. y)
  end function nwad_dble_notequalx
  !>
  !> \brief Return whether \f$x\f$ does not equal \f$y\f$ where the former is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_notequaly(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = .not. (x .eq. y)
  end function nwad_dble_notequaly
  !>
  !> \brief Return whether \f$x\f$ is less than \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_lessthan(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = x%d0 .lt. y%d0
  end function nwad_dble_lessthan
  !>
  !> \brief Return whether \f$x\f$ is less than \f$y\f$ where the latter is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_lessthanx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = x%d0 .lt. y
  end function nwad_dble_lessthanx
  !>
  !> \brief Return whether \f$x\f$ is less than \f$y\f$ where the former is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_lessthany(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = x .lt. y%d0
  end function nwad_dble_lessthany
  !>
  !> \brief Return whether \f$x\f$ is less than or equal to \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than or equal
  !> to\f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_lessequal(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = x%d0 .le. y%d0
  end function nwad_dble_lessequal
  !>
  !> \brief Return whether \f$x\f$ is less than or equal to \f$y\f$ where the
  !> latter is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than or equal
  !> to \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_lessequalx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = x%d0 .le. y
  end function nwad_dble_lessequalx
  !>
  !> \brief Return whether \f$x\f$ is less than or equal to \f$y\f$ where the
  !> former is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is less than or equal
  !> to \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_lessequaly(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = x .le. y%d0
  end function nwad_dble_lessequaly
  !>
  !> \brief Return whether \f$x\f$ is greater than \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_greaterthan(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = x%d0 .lt. y%d0
  end function nwad_dble_greaterthan
  !>
  !> \brief Return whether \f$x\f$ is greater than \f$y\f$ where the latter is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_greaterthanx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = x%d0 .lt. y
  end function nwad_dble_greaterthanx
  !>
  !> \brief Return whether \f$x\f$ is greater than \f$y\f$ where the former is
  !> inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_greaterthany(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = x .lt. y%d0
  end function nwad_dble_greaterthany
  !>
  !> \brief Return whether \f$x\f$ is greater than or equal to \f$y\f$
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than or equal
  !> to\f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives.
  !>
  function nwad_dble_greaterequal(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    logical                     :: s
    s = x%d0 .ge. y%d0
  end function nwad_dble_greaterequal
  !>
  !> \brief Return whether \f$x\f$ is greater than or equal to \f$y\f$ where the
  !> latter is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than or equal
  !> to \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$y\f$
  !> is an inactive variable.
  !>
  function nwad_dble_greaterequalx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    logical                      :: s
    s = x%d0 .ge. y
  end function nwad_dble_greaterequalx
  !>
  !> \brief Return whether \f$x\f$ is greater than or equal to \f$y\f$ where the
  !> former is inactive
  !>
  !> Return a logical value reflecting whether \f$x\f$ is greater than or equal
  !> to \f$y\f$.
  !> In context of derivative calculations the comparison only applies to the
  !> the values and not to any of the derivatives. In this function \f$x\f$
  !> is an inactive variable.
  !>
  function nwad_dble_greaterequaly(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    logical                      :: s
    s = x .ge. y%d0
  end function nwad_dble_greaterequaly
  !>
  !> \brief Evaluate the sign function 
  !> 
  !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
  !> the sign of \f$y\f$. This routine implements this function for the case
  !> where both \f$x\f$ and \f$y\f$ are active variables.
  !>
  function nwad_dble_sign(x,y) result(s)
    type(nwad_dble), intent(in) :: x
    type(nwad_dble), intent(in) :: y
    type(nwad_dble)             :: s
    s = abs(x) * sign(1.0d0,y%d0)
  end function nwad_dble_sign
  !>
  !> \brief Evaluate the sign function where \f$y\f$ is inactive
  !> 
  !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
  !> the sign of \f$y\f$. This routine implements this function for the case
  !> where \f$x\f$ is an active and \f$y\f$ is an inactive variable.
  !>
  function nwad_dble_signx(x,y) result(s)
    type(nwad_dble),  intent(in) :: x
    double precision, intent(in) :: y
    type(nwad_dble)              :: s
    s = abs(x) * sign(1.0d0,y)
  end function nwad_dble_signx
  !>
  !> \brief Evaluate the sign function where \f$x\f$ is inactive
  !> 
  !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
  !> the sign of \f$y\f$. This routine implements this function for the case
  !> where \f$x\f$ is an inactive and \f$y\f$ is an active variable.
  !>
  function nwad_dble_signy(x,y) result(s)
    double precision, intent(in) :: x
    type(nwad_dble),  intent(in) :: y
    double precision             :: s
    s = abs(x) * sign(1.0d0,y%d0)
  end function nwad_dble_signy
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
  !> \brief Evaluate the unary negation operator \f$-\f$
  !>
  !> The unary negation operator simply replaces the value and all the
  !> derivatives with the same with the opposite sign.
  !>
  function nwad_dble_minus(x) result (s)
    type(nwad_dble), intent(in)  :: x
    type(nwad_dble)              :: s
    s%d0 = -x%d0
    s%d1 = -x%d1
    s%d2 = -x%d2
    s%d3 = -x%d3
  end function nwad_dble_minus
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
    s%d0 = exp(x%d0)
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
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \mathrm{asin}(x(r)) \\\\
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
    s%d0 = asin(x%d0)
    s%d1 = x%d1/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
    s%d2 = x%d2/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) &
         + x%d0*(x%d1**2.0d0)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
    s%d3 = x%d3/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) &
         + 3.0d0*x%d0*x%d1*x%d2/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) &
         + (x%d1**3.0d0)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) &
         + 3.0d0*(x%d0**2.0d0)*(x%d1**3.0d0)/((1.0d0-x%d0*x%d0)**(5.0d0/2.0d0))
  end function nwad_dble_asin
  !>
  !> \brief Evaluate the \f$\mathrm{acos}\f$ function
  !>
  !> The implementation of the \f$\mathrm{acos}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{acos}(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \mathrm{acos}(x(r)) \\\\
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
    s%d0 = acos(x%d0)
    s%d1 = -x%d1/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
    s%d2 = -x%d2/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) &
         - x%d0*(x%d1**2.0d0)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
    s%d3 = -x%d3/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) &
         - 3.0d0*x%d0*x%d1*x%d2/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) &
         - (x%d1**3.0d0)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) &
         - 3.0d0*(x%d0**2.0d0)*(x%d1**3.0d0)/((1.0d0-x%d0*x%d0)**(5.0d0/2.0d0))
  end function nwad_dble_acos
  !>
  !> \brief Evaluate the \f$\mathrm{atan}\f$ function
  !>
  !> The implementation of the \f$\mathrm{atan}\f$ function. The chain rule is
  !> used to evaluate the derivatives. I.e. we consider
  !> \f$s(r) = \mathrm{atan}(x(r))\f$.
  !> \f{eqnarray*}{
  !>   \frac{\mathrm{d}^0s(r)}{\mathrm{d}r^0} &=& \mathrm{atan}(x(r)) \\\\
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
    s%d0 = atan(x%d0)
    s%d1 = x%d1/(1.0d0+x%d0*x%d0)
    s%d2 = x%d2/(1.0d0+x%d0*x%d0) &
         - 2.0d0*x%d0*(x%d1**2.0d0)/((1.0d0+x%d0*x%d0)**2.0d0)
    s%d3 = x%d3/(1.0d0+x%d0*x%d0) &
         - 6.0d0*x%d0*x%d1*x%d2/((1.0d0+x%d0*x%d0)**2.0d0) &
         - 2.0d0*(x%d1**3.0d0)/((1.0d0+x%d0*x%d0)**2.0d0) &
         + 8.0d0*(x%d0**2.0d0)*(x%d1**3.0d0)/((1.0d0+x%d0*x%d0)**3.0d0)
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
  !> \brief Initialize an active variable of a given order
  !>
  !> Initialize an active variable. Active variables are those with respect
  !> to which the derivatives are calculated in the current evaluation of the
  !> code. The code calculates the derivatives with respect to some 
  !> multi-index \f$\mathbf{i}\f$. The norm of \f$\mathbf{i}\f$ is given by
  !> \f$|\mathbf{i}|=\sum_{j=1}^n\mathbf{i}_j\f$. All multi indeces are allowed
  !> as long as \f$|\mathbf{i}| \le d\f$ where \f$d\f$ is the order of
  !> differentiation. In particular a component of \f$\mathbf{i}\f$ may occur
  !> more than once. This routine specifies the integer number of times this
  !> particular variables is differentiated with respect to.
  !>
  !> In practice it means that the components are initialized as
  !> \f{eqnarray*}{
  !>   d0 &=& \frac{\mathrm{d}^0 x}{\mathrm{d}x^0} = x \\\\
  !>   d1 &=& \frac{\mathrm{d}^1 x}{\mathrm{d}x^1} = n \\\\
  !>   d2 &=& \frac{\mathrm{d}^2 x}{\mathrm{d}x^2} = 0 \\\\
  !>   d3 &=& \frac{\mathrm{d}^3 x}{\mathrm{d}x^3} = 0 \\\\
  !> \f}
  !> See [3] for details.
  !>
  function nwad_dble_active_n(x,n) result (s)
    double precision, intent(in) :: x
    integer, intent(in)          :: n
    type(nwad_dble)              :: s
    s%d0 = x
    s%d1 = n
    s%d2 = 0
    s%d3 = 0
  end function nwad_dble_active_n
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
  !>
  !> \brief Return the value 
  !>
  !> Simply return the value, i.e. 0th order derivative, of the quantity
  !>
  function nwad_dble_value(x) result (s)
    type(nwad_dble), intent(in) :: x
    double precision            :: s
    s = x%d0
  end function nwad_dble_value
  !>
  !> \brief Interpolate the 1st order partial derivative from an evaluation
  !> up to 1st order
  !>
  !> Interpolate the 1st order partial derivative from a variable involved
  !> in a functional evaluation with up to 1st order partial derivatives.
  !> The corresponding interpolation expression is particularly simply but
  !> this function is provided for completeness.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d1_dx(dx) result (s)
    type(nwad_dble), intent(in) :: dx 
    double precision            :: s
    s = dx%d1
  end function nwad_dble_inter_d1_dx
  !>
  !> \brief Interpolate the 1st order partial derivative from an evaluation
  !> up to 2nd order
  !>
  !> Interpolate the 1st order partial derivative from a variable involved
  !> in a functional evaluation with up to 2nd order partial derivatives.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d2_dx(dx) result (s)
    type(nwad_dble), intent(in) :: dx 
    double precision            :: s
    double precision :: f12
    f12 = 1.0d0/2.0d0
    s = f12*dx%d1
  end function nwad_dble_inter_d2_dx
  !>
  !> \brief Interpolate the diagonal 2nd order partial derivative from an
  !> evaluation up to 2nd order
  !>
  !> Interpolate the 2nd order partial derivative from a variable involved
  !> in a functional evaluation with up to 2nd order partial derivatives.
  !> This routine takes care of the diagonal term, i.e. 
  !> \f$\frac{\partial^2 f}{\partial x^2}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d2_dx2(dx2) result (s)
    type(nwad_dble), intent(in) :: dx2
    double precision            :: s
    double precision :: f14
    f14 = 1.0d0/4.0d0
    s = f14*dx2%d2
  end function nwad_dble_inter_d2_dx2
  !>
  !> \brief Interpolate the off-diagonal 2nd order partial derivative from an
  !> evaluation up to 2nd order
  !>
  !> Interpolate the 2nd order partial derivative from a variable involved
  !> in a functional evaluation with up to 2nd order partial derivatives.
  !> This routine takes care of the diagonal term, i.e. 
  !> \f$\frac{\partial^2 f}{\partial x\partial y}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d2_dxy(dx2,dxy,dy2) result (s)
    type(nwad_dble), intent(in) :: dx2
    type(nwad_dble), intent(in) :: dxy
    type(nwad_dble), intent(in) :: dy2
    double precision            :: s
    double precision :: f12
    double precision :: f18
    f12 = 1.0d0/2.0d0
    f18 = 1.0d0/8.0d0
    s = f12*dxy%d2-f18*dx2%d2-f18*dy2%d2
  end function nwad_dble_inter_d2_dxy
  !>
  !> \brief Interpolate the 1st order partial derivative from an evaluation
  !> up to 3rd order
  !>
  !> Interpolate the 1st order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dx(dx) result (s)
    type(nwad_dble), intent(in) :: dx 
    double precision            :: s
    double precision :: f13
    f13 = 1.0d0/3.0d0
    s = f13*dx%d1
  end function nwad_dble_inter_d3_dx
  !>
  !> \brief Interpolate the diagonal 2nd order partial derivative from an
  !> evaluation up to 3rd order
  !>
  !> Interpolate the 2nd order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> This routine takes care of the diagonal term, i.e. 
  !> \f$\frac{\partial^2 f}{\partial x^2}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dx2(dx2) result (s)
    type(nwad_dble), intent(in) :: dx2
    double precision            :: s
    double precision :: f19
    f19 = 1.0d0/9.0d0
    s = f19*dx2%d2
  end function nwad_dble_inter_d3_dx2
  !>
  !> \brief Interpolate the off-diagonal 2nd order partial derivative from an
  !> evaluation up to 3rd order
  !>
  !> Interpolate the 2nd order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> This routine takes care of the diagonal term, i.e. 
  !> \f$\frac{\partial^2 f}{\partial x\partial y}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dxy(dx3,dx2y,dxy2,dy3) result (s)
    type(nwad_dble), intent(in) :: dx3
    type(nwad_dble), intent(in) :: dx2y
    type(nwad_dble), intent(in) :: dxy2
    type(nwad_dble), intent(in) :: dy3
    double precision            :: s
    double precision :: f18
    double precision :: f572
    f18  = 1.0d0/8.0d0
    f572 = 5.0d0/72.0d0
    s = f18*(dx2y%d2+dxy2%d2)-f572*(dx3%d2+dy3%d2)
  end function nwad_dble_inter_d3_dxy
  !>
  !> \brief Interpolate the diagonal 3rd order partial derivative from an
  !> evaluation up to 3rd order
  !>
  !> Interpolate the 3rd order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> This routine takes care of the diagonal term, i.e. 
  !> \f$\frac{\partial^3 f}{\partial x^3}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dx3(dx3) result (s)
    type(nwad_dble), intent(in) :: dx3
    double precision            :: s
    double precision :: f127
    f127 = 1.0d0/27.0d0
    s = f127*dx3%d3
  end function nwad_dble_inter_d3_dx3
  !>
  !> \brief Interpolate the semi off-diagonal 3rd order partial derivative from
  !> an evaluation up to 3rd order
  !>
  !> Interpolate the 3rd order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> This routine takes care of the semi off-diagonal term, i.e. 
  !> \f$\frac{\partial^3 f}{\partial x^2\partial y}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dx2y(dx3,dx2y,dxy2,dy3) result (s)
    type(nwad_dble), intent(in) :: dx3
    type(nwad_dble), intent(in) :: dx2y
    type(nwad_dble), intent(in) :: dxy2
    type(nwad_dble), intent(in) :: dy3
    double precision            :: s
    double precision :: f118
    double precision :: f181
    double precision :: f19
    double precision :: f5162
    f118  = 1.0d0/18.0d0
    f181  = 1.0d0/81.0d0
    f19   = 1.0d0/9.0d0
    f5162 = 5.0d0/162.0d0
    s = f19*dx2y%d3+f181*dy3%d3-f5162*dx3%d3-f118*dxy2%d3
  end function nwad_dble_inter_d3_dx2y
  !>
  !> \brief Interpolate the off-diagonal 3rd order partial derivative from
  !> an evaluation up to 3rd order
  !>
  !> Interpolate the 3rd order partial derivative from a variable involved
  !> in a functional evaluation with up to 3rd order partial derivatives.
  !> This routine takes care of the off-diagonal term, i.e. 
  !> \f$\frac{\partial^3 f}{\partial x\partial y\partial z}\f$.
  !> See also Eq.(17) [1].
  !>
  !> ### References ###
  !>
  !> [1] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
  !>     tensors by forward propagation of univariate Taylor series",
  !>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
  !>     <a href="http://dx.doi.org/10.1090/S0025-5718-00-01120-0">
  !>     10.1090/S0025-5718-00-01120-0</a>.
  !>
  function nwad_dble_inter_d3_dxyz(dx3,dy3,dz3,dx2y,dx2z,dxy2,dy2z,dxz2,dyz2, &
                                   dxyz) result (s)
    type(nwad_dble), intent(in) :: dx3, dy3, dz3
    type(nwad_dble), intent(in) :: dx2y, dx2z, dxy2, dy2z, dxz2, dyz2
    type(nwad_dble), intent(in) :: dxyz
    double precision            :: s
    double precision :: f136
    double precision :: f16
    double precision :: f181
    f136  = 1.0d0/36.0d0
    f16   = 1.0d0/6.0d0
    f181  = 1.0d0/81.0d0
    s = f16*dxyz%d3+f181*(dx3%d3+dy3%d3+dz3%d3) &
      - f136*(dx2y%d3+dx2z%d3+dxy2%d3+dy2z%d3+dxz2%d3+dyz2%d3)
  end function nwad_dble_inter_d3_dxyz
end module nwad
!> @}
