#if !defined(NWAD_DEG)
!>
!> \file wfn1_ad.F
!> wfn1 automatic differentiation for matrix function derivatives
!>
!> \ingroup wfn1_nxt2
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
!> The package was designed to generate compact data structures to keep the
!> overhead of copying data low. In addition applications are targetted where
!> all partial derivatives of an expression are needed. The method chosen 
!> required multiple executions of the function to obtain all its derivatives.
!> Every function execution computed only one specific partial derivative [3].
!> The downside of this approach was that it introduced significant overheads
!> that increased with increasing orders of derivatives.
!>
!> The performance problem was addressed by switching to a regular multi-variate
!> differentiation approach. While the code for each of the elementary steps
!> becomes more complicated duplication of work is avoided. In addition 
!> "sparcity" can be exploited, i.e. if we have an expression like
!> \f{eqnarray*}{
!>   f(x,y) &=& y*sin(x)
!> \f}
!> then the derivative of the \f$sin\f$ function wrt. \f$y\f$ is 0 and we should
!> not need to calculate that. Keeping track of these dependencies introduces
!> some overhead but we expect to save more by avoiding unnecessary work.
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
!>     <a href="https://doi.org/10.1145/355586.364791">
!>     10.1145/355586.364791</a>.
!>
!> [2] See e.g. <a href="http://www.autodiff.org/">www.autodiff.org</a>
!>
!> [3] A. Griewank, J. Utke, A. Walther (2000) "Evaluating higher derivative
!>     tensors by forward propagation of univariate Taylor series",
!>     Mathematics of Computation, <b>69</b>, pp. 1117-1130, DOI:
!>     <a href="https://doi.org/10.1090/S0025-5718-00-01120-0">
!>     10.1090/S0025-5718-00-01120-0</a>.
!>
!> [4] A. Griewank, D. Juedes, J. Utke (1996) "Algorithm 755: ADOL-C: A package
!>     for the automatic differentiation of algorithms written in C/C++",
!>     ACM Transactions on Mathematical Software, <b>22</b>, pp. 131-167, DOI:
!>     <a href="https://doi.org/10.1145/229473.229474">
!>     10.1145/229473.229474</a>.
!>
!> [5] I. Charpentier, J. Utke, "Rapsodia: User Manual", Argonne National
!>     Laboratory, <a href="http://www.mcs.anl.gov/Rapsodia/userManual.pdf">
!>     http://wwww.mcs.anl.gov/Rapsodia/userManual.pdf</a> (referenced
!>     July 3, 2014).
!>
!> $Id: wfn1_ad.F,v 1.5 2015/11/21 22:16:22 hvandam Exp $
!>
!> Huub van Dam, 2014
!>
#endif
#if !defined(NWAD_DEG)
#define NWAD_DEG 0
#endif
#if NWAD_DEG <= 0
      module wfn1_ad0
#elif NWAD_DEG <= 1
      module wfn1_ad1
#elif NWAD_DEG <= 2
      module wfn1_ad2
#else
      module wfn1_ad3
#endif
      !> The maximum number of variables which determines the sizes of the
      !> data arrays. 
      integer mxvar
      parameter (mxvar = 15*15) ! E.g. Be in 6-31G* (3S,2P,1D)
      !>
      !> The data type to hold a double precision value and the derivatives of this
      !> expression with respect to the sum of the active variables. I.e. the
      !> member \f$\mathrm{d}n\f$ holds
      !> \f{eqnarray}{
      !>   \mathrm{d}n &=& \frac{\mathrm{d}^n f}{\mathrm{d}(\sum_i x_i)^n}
      !> \f}
      !> where \f$x_i\f$ are all active variables.
      !>
      !> In addition we store how many variables this quantity depends on
      !> as well as which variables they are. The list of variables is 
      !> ordered in increasing order.
      !>
      type :: wfn1_ad_dble
        integer(kind=4)  :: nvar
        integer(kind=4)  :: ivar(mxvar)
        double precision :: d0
#if NWAD_DEG >= 1
        double precision :: d1(mxvar)
#endif
#if NWAD_DEG >= 2
        double precision :: d2(mxvar*(mxvar+1)/2)
#endif
#if NWAD_DEG >= 3
        double precision :: d3(mxvar*(mxvar+1)*(mxvar+2)/6)
#endif
      end type wfn1_ad_dble
      interface assignment (=)
        module procedure wfn1_ad_dble_assignx
        module procedure wfn1_ad_dble_assign
      end interface
      interface max
        module procedure wfn1_ad_dble_max
      end interface
      interface min
        module procedure wfn1_ad_dble_min
      end interface
      interface operator (+)
        module procedure wfn1_ad_dble_add
        module procedure wfn1_ad_dble_addx
        module procedure wfn1_ad_dble_addy
      end interface
      interface operator (-)
        module procedure wfn1_ad_dble_minus
        module procedure wfn1_ad_dble_sub
        module procedure wfn1_ad_dble_subx
        module procedure wfn1_ad_dble_suby
      end interface
      interface operator (*)
        module procedure wfn1_ad_dble_mult
        module procedure wfn1_ad_dble_multx
        module procedure wfn1_ad_dble_multy
      end interface
      interface operator (/)
        module procedure wfn1_ad_dble_div
        module procedure wfn1_ad_dble_divx
        module procedure wfn1_ad_dble_divy
      end interface
      interface operator (**)
        module procedure wfn1_ad_dble_pow
        module procedure wfn1_ad_dble_powx
        module procedure wfn1_ad_dble_powy
      end interface
      interface operator (.eq.)
        module procedure wfn1_ad_dble_equal
        module procedure wfn1_ad_dble_equalx
        module procedure wfn1_ad_dble_equaly
      end interface operator (.eq.)
      interface operator (.ne.)
        module procedure wfn1_ad_dble_notequal
        module procedure wfn1_ad_dble_notequalx
        module procedure wfn1_ad_dble_notequaly
      end interface operator (.ne.)
      interface operator (.lt.)
        module procedure wfn1_ad_dble_lessthan
        module procedure wfn1_ad_dble_lessthanx
        module procedure wfn1_ad_dble_lessthany
      end interface operator (.lt.)
      interface operator (.le.)
        module procedure wfn1_ad_dble_lessequal
        module procedure wfn1_ad_dble_lessequalx
        module procedure wfn1_ad_dble_lessequaly
      end interface operator (.le.)
      interface operator (.gt.)
        module procedure wfn1_ad_dble_greaterthan
        module procedure wfn1_ad_dble_greaterthanx
        module procedure wfn1_ad_dble_greaterthany
      end interface operator (.gt.)
      interface operator (.ge.)
        module procedure wfn1_ad_dble_greaterequal
        module procedure wfn1_ad_dble_greaterequalx
        module procedure wfn1_ad_dble_greaterequaly
      end interface operator (.ge.)
      interface sign
        module procedure wfn1_ad_dble_sign
        module procedure wfn1_ad_dble_signx
        module procedure wfn1_ad_dble_signy
      end interface
      interface abs
        module procedure wfn1_ad_dble_abs
      end interface
      interface exp
        module procedure wfn1_ad_dble_exp
      end interface
      interface sqrt
        module procedure wfn1_ad_dble_sqrt
      end interface
      interface log
        module procedure wfn1_ad_dble_log
      end interface
! interface log10
!   module procedure wfn1_ad_dble_log10
! end interface
      interface sin
        module procedure wfn1_ad_dble_sin
      end interface
      interface cos
        module procedure wfn1_ad_dble_cos
      end interface
      interface tan
        module procedure wfn1_ad_dble_tan
      end interface
      interface asin
        module procedure wfn1_ad_dble_asin
      end interface
      interface acos
        module procedure wfn1_ad_dble_acos
      end interface
      interface atan
        module procedure wfn1_ad_dble_atan
      end interface
      interface sinh
        module procedure wfn1_ad_dble_sinh
      end interface
      interface cosh
        module procedure wfn1_ad_dble_cosh
      end interface
      interface tanh
        module procedure wfn1_ad_dble_tanh
      end interface
      interface asinh
        module procedure wfn1_ad_dble_asinh
      end interface
      interface erf
        module procedure wfn1_ad_dble_erf
      end interface
      interface erfc
        module procedure wfn1_ad_dble_erfc
      end interface
      interface set_dij
        module procedure wfn1_ad_dble_dij
      end interface
      interface get_nvar
        module procedure wfn1_ad_dble_get_nvar
      end interface
      interface value
        module procedure wfn1_ad_dble_get_val
      end interface
      interface get_val
        module procedure wfn1_ad_dble_get_val
      end interface
#if NWAD_DEG >= 1
      interface get_d1
        module procedure wfn1_ad_dble_get_d1
      end interface
#endif
#if NWAD_DEG >= 2
      interface get_d2
        module procedure wfn1_ad_dble_get_d2
      end interface
#endif
#if NWAD_DEG >= 3
      interface get_d3
        module procedure wfn1_ad_dble_get_d3
      end interface
#endif
      interface print
        module procedure wfn1_ad_dble_print
      end interface
      contains
      !>
      !> \brief Merge two sets of variables into one
      !>
      !> When using multi-variate differentiation and exploiting the sparcity
      !> of the expressions we have to assume that the value of each variable
      !> depends only on a subset of all variables. The result of binary
      !> operators then depends on the union of the variables that the inputs
      !> depend on. Hence we need to merge the dependencies of the input
      !> variables to construct the set of dependencies for the result.
      !> This subroutine executes this merger and also establishes for each
      !> rank of input dependency where this contributes to in the result.
      !>
      subroutine wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
        type(wfn1_ad_dble), intent(in)  :: x !< [Input] Variable X
        type(wfn1_ad_dble), intent(in)  :: y !< [Input] Variable Y
        type(wfn1_ad_dble), intent(out) :: s !< [Output] Variable S, the result
        integer,         intent(out) :: ixto(mxvar) !< [Output] Where X contributes
        integer,         intent(out) :: iyto(mxvar) !< [Output] Where Y contributes
        integer                      :: ix ! local counter
        integer                      :: iy ! local counter
        ! The result of the addition depends on the union of the sets of 
        ! variables that X and Y depend on. Here we construct this union as
        ! an ordered list.
c       call wfn1_ad_test(x)
c       call wfn1_ad_test(y)
        ix = 1
        iy = 1
        s%nvar = 0
        do while (ix.le.x%nvar.and.iy.le.y%nvar)
          if (x%ivar(ix).eq.y%ivar(iy)) then
            s%nvar         = s%nvar + 1
            s%ivar(s%nvar) = x%ivar(ix)
            ixto(ix)       = s%nvar
            iyto(iy)       = s%nvar
            ix             = ix + 1
            iy             = iy + 1
          else if (x%ivar(ix).lt.y%ivar(iy)) then
            s%nvar         = s%nvar + 1
            s%ivar(s%nvar) = x%ivar(ix)
            ixto(ix)       = s%nvar
            ix             = ix + 1
          else ! (x%ivar(ix).gt.y%ivar(iy))
            s%nvar         = s%nvar + 1
            s%ivar(s%nvar) = y%ivar(iy)
            iyto(iy)       = s%nvar
            iy             = iy + 1
          endif
        enddo
        do while (ix.le.x%nvar)
          s%nvar         = s%nvar + 1
          s%ivar(s%nvar) = x%ivar(ix)
          ixto(ix)       = s%nvar
          ix             = ix + 1
        enddo
        do while (iy.le.y%nvar)
          s%nvar         = s%nvar + 1
          s%ivar(s%nvar) = y%ivar(iy)
          iyto(iy)       = s%nvar
          iy             = iy + 1
        enddo
      end subroutine wfn1_ad_merge_dependencies
      !>
      !> \brief Copy the dependencies
      !>
      !> Just a couple of lines of code to copy the dependencies between
      !> wfn1_ad variables.
      !>
      subroutine wfn1_ad_copy_dependencies(s,x)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble), intent(out) :: s
        integer                      :: i
c       call wfn1_ad_test(x)
        s%nvar = x%nvar
        do i = 1, x%nvar
          s%ivar(i) = x%ivar(i)
        enddo
      end subroutine wfn1_ad_copy_dependencies
      !>
      !> \brief Assign a value to an inactive variable
      !>
      !> Assign a floating point value to an wfn1_ad variable. This operation 
      !> generates an inactive variable. I.e. a variable for which no derivatives
      !> are calculated. In practice it means that the components are initialized as
      !> \f{eqnarray*}{
      !>   d0 &=& x \\\\
      !>   d1 &=& 0 \\\\
      !>   d2 &=& 0 \\\\
      !>   d3 &=& 0 \\\\
      !> \f}
      !> As this variable does not depend on any other variables \f$nvar = 0\f$.
      !>
      subroutine wfn1_ad_dble_assignx(s,x)
        double precision, intent(in) :: x
        type(wfn1_ad_dble), intent(out) :: s
        s%nvar = 0
        s%d0   = x
#if NWAD_DEG >= 1
        s%d1   = 0
#endif
#if NWAD_DEG >= 2
        s%d2   = 0
#endif
#if NWAD_DEG >= 3
        s%d3   = 0
#endif
      end subroutine wfn1_ad_dble_assignx
      !>
      !> \brief Assign a value to an active variable
      !>
      !> Assign an wfn1_ad value to an wfn1_ad variable. This operation 
      !> essentially copies the value of one variable to another minimizing 
      !> the amount of data moved.
      !>
      subroutine wfn1_ad_dble_assign(s,x)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble), intent(out) :: s
        integer i
        s%d0   = x%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do i = 1, x%nvar
          s%d1(i)   = x%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, x%nvar*(x%nvar+1)/2
          s%d2(i) = x%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(i) = x%d3(i)
        enddo
#endif
      end subroutine wfn1_ad_dble_assign
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
      function wfn1_ad_dble_max(a,b,c,d,e) result (s)
        type(wfn1_ad_dble), intent(in)           :: a
        type(wfn1_ad_dble), intent(in)           :: b
        type(wfn1_ad_dble), intent(in), optional :: c
        type(wfn1_ad_dble), intent(in), optional :: d
        type(wfn1_ad_dble), intent(in), optional :: e
        type(wfn1_ad_dble)                       :: s
        type(wfn1_ad_dble)                       :: t1
        type(wfn1_ad_dble)                       :: t2
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
      end function wfn1_ad_dble_max
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
      function wfn1_ad_dble_min(a,b,c,d,e) result (s)
        type(wfn1_ad_dble), intent(in)           :: a
        type(wfn1_ad_dble), intent(in)           :: b
        type(wfn1_ad_dble), intent(in), optional :: c
        type(wfn1_ad_dble), intent(in), optional :: d
        type(wfn1_ad_dble), intent(in), optional :: e
        type(wfn1_ad_dble)                       :: s
        type(wfn1_ad_dble)                       :: t1
        type(wfn1_ad_dble)                       :: t2
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
      end function wfn1_ad_dble_min
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
      function wfn1_ad_dble_add(x,y) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        integer                     :: ix, ix2, ix3
        integer                     :: iy, iy2, iy3
        integer                     :: ixto(mxvar) ! contribution from X goes
                                                   ! S element IXTO
        integer                     :: iyto(mxvar) ! contribution from Y goes
                                                   ! S element IYTO
        s%d0 = x%d0 + y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
        do ix = 1, s%nvar
          s%d1(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          s%d1(ixto(ix)) = s%d1(ixto(ix)) + x%d1(ix)
        enddo
        do iy = 1, y%nvar
          s%d1(iyto(iy)) = s%d1(iyto(iy)) + y%d1(iy)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, s%nvar*(s%nvar+1)/2
          s%d2(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2)) 
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      + x%d2(ix*(ix-1)/2+ix2)
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2)) 
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      + y%d2(iy*(iy-1)/2+iy2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, s%nvar*(s%nvar+1)*(s%nvar+2)/6
          s%d3(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3(ixto(ix)*(ixto(ix)-1)*(ixto(ix)+1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3)) 
     +        = s%d3(ixto(ix)*(ixto(ix)-1)*(ixto(ix)+1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        + x%d3(ix*(ix-1)*(ix+1)/6+ix2*(ix2-1)/2+ix3)
            enddo
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3(iyto(iy)*(iyto(iy)-1)*(iyto(iy)+1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3)) 
     +        = s%d3(iyto(iy)*(iyto(iy)-1)*(iyto(iy)+1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        + y%d3(iy*(iy-1)*(iy+1)/6+iy2*(iy2-1)/2+iy3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_add
      !>
      !> \brief Evaluate the addition where \f$y\f$ is inactive
      !>
      !> This routine does the same as wfn1_ad_dble_add but \f$y\f$ is inactive
      !> and therefore all its derivatives are 0.
      !>
      function wfn1_ad_dble_addx(x,y) result (s)
        type(wfn1_ad_dble),  intent(in) :: x
        double precision, intent(in) :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: i
        s%d0 = x%d0 + y
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do i = 1, x%nvar
          s%d1(i) = x%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, x%nvar*(x%nvar+1)/2
          s%d2(i) = x%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(i) = x%d3(i)
        enddo
#endif
      end function wfn1_ad_dble_addx
      !>
      !> \brief Evaluate the addition where \f$x\f$ is inactive
      !>
      !> This routine does the same as wfn1_ad_dble_add but \f$x\f$ is inactive
      !> and therefore all its derivatives are 0.
      !>
      function wfn1_ad_dble_addy(x,y) result (s)
        double precision, intent(in) :: x
        type(wfn1_ad_dble),  intent(in) :: y
        type(wfn1_ad_dble)              :: s
        s%d0 = x + y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,y)
        do i = 1, y%nvar
          s%d1(i) = y%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, y%nvar*(y%nvar+1)/2
          s%d2(i) = y%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, y%nvar*(y%nvar+1)*(y%nvar+2)/6
          s%d3(i) = y%d3(i)
        enddo
#endif
      end function wfn1_ad_dble_addy
      !>
      !> \brief Evaluation of the subtraction operator
      !>
      !> The implementation of the subtraction operator. The chain rule is used to
      !> evaluate the derivatives. Obviously this is very similar to the 
      !> addition operator apart from the minus sign.
      !>
      function wfn1_ad_dble_sub(x,y) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        integer                     :: ix, ix2, ix3
        integer                     :: iy, iy2, iy3
        integer                     :: ixto(mxvar) ! contribution from X goes
                                                   ! S element IXTO
        integer                     :: iyto(mxvar) ! contribution from Y goes
                                                   ! S element IYTO
        s%d0 = x%d0 - y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
        do ix = 1, s%nvar
          s%d1(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          s%d1(ixto(ix)) = s%d1(ixto(ix)) + x%d1(ix)
        enddo
        do iy = 1, y%nvar
          s%d1(iyto(iy)) = s%d1(iyto(iy)) - y%d1(iy)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, s%nvar*(s%nvar+1)/2
          s%d2(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2)) 
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      + x%d2(ix*(ix-1)/2+ix2)
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2)) 
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      - y%d2(iy*(iy-1)/2+iy2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, s%nvar*(s%nvar+1)*(s%nvar+2)/6
          s%d3(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3(ixto(ix)*(ixto(ix)-1)*(ixto(ix)+1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3)) 
     +        = s%d3(ixto(ix)*(ixto(ix)-1)*(ixto(ix)+1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        + x%d3(ix*(ix-1)*(ix+1)/6+ix2*(ix2-1)/2+ix3)
            enddo
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3(iyto(iy)*(iyto(iy)-1)*(iyto(iy)+1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3)) 
     +        = s%d3(iyto(iy)*(iyto(iy)-1)*(iyto(iy)+1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        - y%d3(iy*(iy-1)*(iy+1)/6+iy2*(iy2-1)/2+iy3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_sub
      !>
      !> \brief Evaluation of the subtraction operator where \f$y\f$ is inactive
      !>
      !> This function is similar to wfn1_ad_dble_sub but because \f$y\f$ is in
      !> inactive all derivatives of \f$y\f$ vanish.
      !>
      function wfn1_ad_dble_subx(x,y) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        double precision, intent(in) :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: i
        s%d0 = x%d0 - y
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do i = 1, x%nvar
          s%d1(i) = x%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, x%nvar*(x%nvar+1)/2
          s%d2(i) = x%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(i) = x%d3(i)
        enddo
#endif
      end function wfn1_ad_dble_subx
      !>
      !> \brief Evaluation of the subtraction operator where \f$x\f$ is inactive
      !>
      !> This function is similar to wfn1_ad_dble_sub but because \f$x\f$ is in
      !> inactive all derivatives of \f$x\f$ vanish.
      !>
      function wfn1_ad_dble_suby(x,y) result (s)
        double precision, intent(in) :: x
        type(wfn1_ad_dble), intent(in)  :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: i
        s%d0 = x - y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,y)
        do i = 1, y%nvar
          s%d1(i) = - y%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, y%nvar*(y%nvar+1)/2
          s%d2(i) = - y%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, y%nvar*(y%nvar+1)*(y%nvar+2)/6
          s%d3(i) = - y%d3(i)
        enddo
#endif
      end function wfn1_ad_dble_suby
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
      function wfn1_ad_dble_mult(x,y) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        integer                     :: ix, ix2, ix3
        integer                     :: iy, iy2, iy3
        integer                     :: ixto(mxvar) ! contribution from X goes
                                                   ! S element IXTO
        integer                     :: iyto(mxvar) ! contribution from Y goes
                                                   ! S element IYTO
        s%d0 =    x%d0 * y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
        do ix = 1, s%nvar
          s%d1(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          s%d1(ixto(ix)) = s%d1(ixto(ix)) + x%d1(ix) * y%d0
        enddo
        do iy = 1, y%nvar
          s%d1(iyto(iy)) = s%d1(iyto(iy)) + x%d0 * y%d1(iy)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, s%nvar*(s%nvar+1)/2
          s%d2(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      + x%d2(ix*(ix-1)/2+ix2) * y%d0
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while(iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      + x%d1(ix) * y%d1(iy)
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while(ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      + x%d1(ix) * y%d1(iy)
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      + x%d0 * y%d2(iy*(iy-1)/2+iy2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, s%nvar*(s%nvar+1)*(s%nvar+2)/6
          s%d3(ix) = 0.0d0
        enddo
c
c       Dxxx
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        + x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3) * y%d0
            enddo
          enddo
        enddo
c
c       Dxxy
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            iy = 1
            do while (iyto(iy).le.ixto(ix2).and.iy.le.y%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        + x%d2(ix*(ix-1)/2+ix2) * y%d1(iy)
              iy = iy + 1
            enddo
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            ix2 = 1
            do while (ixto(ix2).le.iyto(iy).and.ix2.le.x%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        + x%d2(ix*(ix-1)/2+ix2) * y%d1(iy)
              ix2 = ix2 + 1
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            do ix2 = 1, ix
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        + x%d2(ix*(ix-1)/2+ix2) * y%d1(iy)
            enddo
            ix = ix + 1
          enddo
        enddo
c
c       Dxyy
c
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            do iy2 = 1, iy
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        + x%d1(ix) * y%d2(iy*(iy-1)/2+iy2)
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            iy2 = 1
            do while (iyto(iy2).le.ixto(ix).and.iy2.le.y%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        + x%d1(ix) * y%d2(iy*(iy-1)/2+iy2)
              iy2 = iy2 + 1
            enddo
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            ix = 1
            do while (ixto(ix).le.iyto(iy2).and.ix.le.x%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        + x%d1(ix) * y%d2(iy*(iy-1)/2+iy2)
              ix = ix + 1
            enddo
          enddo
        enddo
c
c       Dyyy
c
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        + x%d0 * y%d3((iy+1)*iy*(iy-1)/6+iy2*(iy2-1)/2+iy3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_mult
      !>
      !> \brief Evaluate the multiplication operator where \f$y\f$ is inactive
      !>
      !> This is similar to the regular multiplication operator except that all
      !> derivative of \f$y\f$ are zero.
      !>
      function wfn1_ad_dble_multx(x,y) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        double precision, intent(in) :: y
        type(wfn1_ad_dble)              :: s
        s%d0 = x%d0 * y
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do i = 1, x%nvar
          s%d1(i) = x%d1(i) * y
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, x%nvar*(x%nvar+1)/2
          s%d2(i) = x%d2(i) * y
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(i) = x%d3(i) * y
        enddo
#endif
      end function wfn1_ad_dble_multx
      !>
      !> \brief Evaluate the multiplication operator where \f$x\f$ is inactive
      !>
      !> This is similar to the regular multiplication operator except that all
      !> derivative of \f$x\f$ are zero.
      !>
      function wfn1_ad_dble_multy(x,y) result (s)
        double precision, intent(in) :: x
        type(wfn1_ad_dble), intent(in)  :: y
        type(wfn1_ad_dble)              :: s
        s%d0 = x * y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,y)
        do i = 1, y%nvar
          s%d1(i) = x * y%d1(i)
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, y%nvar*(y%nvar+1)/2
          s%d2(i) = x * y%d2(i)
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, y%nvar*(y%nvar+1)*(y%nvar+2)/6
          s%d3(i) = x * y%d3(i)
        enddo
#endif
      end function wfn1_ad_dble_multy
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
      function wfn1_ad_dble_div(x,y) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        integer                     :: ix, ix2, ix3
        integer                     :: iy, iy2, iy3
        integer                     :: ixto(mxvar) ! contribution from X goes
                                                   ! S element IXTO
        integer                     :: iyto(mxvar) ! contribution from Y goes
                                                   ! S element IYTO
        s%d0 = x%d0/y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
!       s%d1 = x%d1/y%d0               -     x%d0 * y%d1 / (y%d0 ** 2)
        do ix = 1, s%nvar
          s%d1(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          s%d1(ixto(ix)) = s%d1(ixto(ix)) + x%d1(ix) / y%d0
        enddo
        do iy = 1, y%nvar
          s%d1(iyto(iy)) = s%d1(iyto(iy)) - x%d0 * y%d1(iy) / (y%d0**2)
        enddo
#endif
#if NWAD_DEG >= 2
!       s%d2 = x%d2/y%d0             - 2*x%d1*y%d1/(y%d0 ** 2)
!    &       - x%d0*y%d2/(y%d0 ** 2) + 2*x%d0*y%d1**2/(y%d0 ** 3)
        do ix = 1, s%nvar*(s%nvar+1)/2
          s%d2(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      + x%d2(ix*(ix-1)/2+ix2) / y%d0
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while(iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      - x%d1(ix) * y%d1(iy) / (y%d0 ** 2)
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while(ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      - x%d1(ix) * y%d1(iy) / (y%d0 ** 2)
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      - x%d0 * y%d2(iy*(iy-1)/2+iy2) / (y%d0 ** 2)
     +      + 2 * x%d0 * y%d1(iy) * y%d1(iy2) / (y%d0 ** 3)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
!       s%d3 =     x%d3               /  y%d0        
!    &       - 3 * x%d2 * y%d1        / (y%d0 ** 2)  
!    &       + 6 * x%d1 * y%d1 ** 2   / (y%d0 ** 3)  
!    &       - 6 * x%d0 * y%d1 ** 3   / (y%d0 ** 4)  
!    &       - 3 * x%d1 * y%d2        / (y%d0 ** 2)  
!    &       + 6 * x%d0 * y%d1 * y%d2 / (y%d0 ** 3)  
!    &       -     x%d0 * y%d3        / (y%d0 ** 2)
        do ix = 1, s%nvar*(s%nvar+1)*(s%nvar+2)/6
          s%d3(ix) = 0.0d0
        enddo
c
c       Dxxx
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        + x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3) / y%d0
            enddo
          enddo
        enddo
c
c       Dxxy
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            iy = 1
            do while (iyto(iy).le.ixto(ix2).and.iy.le.y%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        - x%d2(ix*(ix-1)/2+ix2) * y%d1(iy) / (y%d0 ** 2)
              iy = iy + 1
            enddo
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            ix2 = 1
            do while (ixto(ix2).le.iyto(iy).and.ix2.le.x%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        - x%d2(ix*(ix-1)/2+ix2) * y%d1(iy) / (y%d0 ** 2)
              ix2 = ix2 + 1
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            do ix2 = 1, ix
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        - x%d2(ix*(ix-1)/2+ix2) * y%d1(iy) / (y%d0 ** 2)
            enddo
            ix = ix + 1
          enddo
        enddo
c
c       Dxyy
c
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            do iy2 = 1, iy
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        + 2 * x%d1(ix) * y%d1(iy) * y%d1(iy2) / (y%d0 ** 3)
     +        -     x%d1(ix) * y%d2(iy*(iy-1)/2+iy2) / (y%d0 ** 2)
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            iy2 = 1
            do while (iyto(iy2).le.ixto(ix).and.iy2.le.y%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        + 2 * x%d1(ix) * y%d1(iy) * y%d1(iy2) / (y%d0 ** 3)
     +        -     x%d1(ix) * y%d2(iy*(iy-1)/2+iy2) / (y%d0 ** 2)
              iy2 = iy2 + 1
            enddo
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            ix = 1
            do while (ixto(ix).le.iyto(iy2).and.ix.le.x%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        + 2 * x%d1(ix) * y%d1(iy) * y%d1(iy2) / (y%d0 ** 3)
     +        -     x%d1(ix) * y%d2(iy*(iy-1)/2+iy2) / (y%d0 ** 2)
              ix = ix + 1
            enddo
          enddo
        enddo
c
c       Dyyy
c
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        - 6 * x%d0 * y%d1(iy) * y%d1(iy2) * y%d1(iy3) 
     +          / (y%d0 ** 4)
     +        + 6 * x%d0 * y%d1(iy) * y%d2(iy2*(iy2-1)/2+iy3)
     +          / (y%d0 ** 3)
     +        - x%d0 * y%d3((iy+1)*iy*(iy-1)/6+iy2*(iy2-1)/2+iy3) 
     +          / (y%d0 ** 2)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_div
      !>
      !> \brief Evaluate the division operator where \f$y\f$ is inactive
      !>
      !> This function is particularly simple. It is essentially the same as
      !> multiplying with \f$1/y\f$ while \f$y\f$ is inactive.
      !>
      function wfn1_ad_dble_divx(x,y) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        double precision, intent(in) :: y
        type(wfn1_ad_dble)              :: s
        s%d0 = x%d0 / y
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do i = 1, x%nvar
          s%d1(i) = x%d1(i) / y
        enddo
#endif
#if NWAD_DEG >= 2
        do i = 1, x%nvar*(x%nvar+1)/2
          s%d2(i) = x%d2(i) / y
        enddo
#endif
#if NWAD_DEG >= 3
        do i = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(i) = x%d3(i) / y
        enddo
#endif
      end function wfn1_ad_dble_divx
      !>
      !> \brief Evaluate the division operator where \f$x\f$ is inactive
      !>
      !> Similar to wfn1_ad_dble_div but now \f$x\f$ is inactive. I.e. we consider
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
      function wfn1_ad_dble_divy(x,y) result (s)
        double precision, intent(in) :: x
        type(wfn1_ad_dble), intent(in)  :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: iy, iy2, iy3
        s%d0 =   x / y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,y)
        do iy = 1, y%nvar
          s%d1(iy) = - x * y%d1(iy) / (y%d0 ** 2)
        enddo
#endif
#if NWAD_DEG >= 2
        do iy = 1, y%nvar
          do iy2= 1, iy
            s%d2(iy*(iy-1)/2+iy2) = 2*x*y%d1(iy)*y%d1(iy2)/(y%d0**3)
     +                            - x*y%d2(iy*(iy-1)/2+iy2)/(y%d0**2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do iy = 1, y%nvar
          do iy2= 1, iy
            do iy3= 1, iy2
              s%d3((iy-1)*iy*(iy+1)/6+iy2*(iy2-1)/2+iy3)
     +        = -x*y%d3((iy-1)*iy*(iy+1)/6+iy2*(iy2-1)/2+iy3)/(y%d0**2)
c
     &        + 2*x*y%d1(iy)*y%d2(iy2*(iy2-1)/2+iy3)/(y%d0**3)  
     &        + 2*x*y%d1(iy2)*y%d2(iy*(iy-1)/2+iy3)/(y%d0**3)  
     &        + 2*x*y%d1(iy3)*y%d2(iy*(iy-1)/2+iy2)/(y%d0**3)  
c
     &        - 6*x*y%d1(iy)*y%d1(iy2)*y%d1(iy3)/(y%d0**4)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_divy
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
      function wfn1_ad_dble_pow(x,y) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        integer                     :: ix, ix2, ix3
        integer                     :: iy, iy2, iy3
        integer                     :: ixto(mxvar) ! contribution from X goes
                                                   ! S element IXTO
        integer                     :: iyto(mxvar) ! contribution from Y goes
                                                   ! S element IYTO
        s%d0 = x%d0**y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_merge_dependencies(s,x,y,ixto,iyto)
!       s%d1 = x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)
        do ix = 1, s%nvar
          s%d1(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          s%d1(ixto(ix)) = s%d1(ixto(ix))
     +                   + (x%d0**y%d0)*(y%d0/x%d0)*x%d1(ix)
        enddo
        do iy = 1, y%nvar
          s%d1(iyto(iy)) = s%d1(iyto(iy))
     +                   + (x%d0**y%d0)*log(x%d0)*y%d1(iy)
        enddo
#endif
#if NWAD_DEG >= 2
!       s%d2 = x%d0**y%d0*(log(x%d0)*y%d2 + (2.0d0/x%d0)*x%d1*y%d1 
!    &                    +(y%d0/x%d0)*x%d2 - (y%d0/x%d0**2)*x%d1**2) 
!    &       + x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)**2
        do ix = 1, s%nvar*(s%nvar+1)/2
          s%d2(ix) = 0.0d0
        enddo
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +      + (x%d0**y%d0)*(y%d0/x%d0)*x%d2(ix*(ix-1)/2+ix2)
     +      - (x%d0**y%d0)*(y%d0/(x%d0**2))*x%d1(ix)*x%d1(ix2)
     +      + (x%d0**y%d0)*((y%d0/x%d0)**2)*x%d1(ix)*x%d1(ix2)
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while(iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      = s%d2(ixto(ix)*(ixto(ix)-1)/2+iyto(iy))
     +      + (x%d0**y%d0)/x%d0*x%d1(ix)*y%d1(iy)
     +      + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*x%d1(ix)*y%d1(iy)
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while(ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+ixto(ix))
     +      + (x%d0**y%d0)/x%d0*x%d1(ix)*y%d1(iy)
     +      + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*x%d1(ix)*y%d1(iy)
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      = s%d2(iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +      + (x%d0**y%d0)*log(x%d0)*y%d2(iy*(iy-1)/2+iy2)
     +      + (x%d0**y%d0)*(log(x%d0)**2)*y%d1(iy)*y%d1(iy2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
!       s%d3 = x%d0**y%d0* 
!    &         (log(x%d0)*y%d3 + (3.0d0/x%d0)*x%d1*y%d2 
!    &          +(3.0d0/x%d0)*x%d2*y%d1 + (y%d0/x%d0)*x%d3 
!    &          -(3.0d0/x%d0**2)*x%d1**2*y%d1 
!    &          -(3*y%d0/x%d0**2)*x%d1*x%d2 
!    &          +(2*y%d0/x%d0**3)*x%d1**3) 
!    &       + 3*x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1) * 
!    &         (log(x%d0)*y%d2 + (2.0d0/x%d0)*x%d1*y%d1 
!    &          +(y%d0/x%d0)*x%d2 - (y%d0/x%d0**2)*x%d1**2) 
!    &       + x%d0**y%d0*(log(x%d0)*y%d1 + (y%d0/x%d0)*x%d1)**3
        do ix = 1, s%nvar*(s%nvar+1)*(s%nvar+2)/6
          s%d3(ix) = 0.0d0
        enddo
c
c       Dxxx
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+ixto(ix3))
     +        + (x%d0**y%d0)*(y%d0/x%d0)*
     +          x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
c
c             d1*d2 needs permutation of arguments here
c
     +        - (x%d0**y%d0)*(y%d0/(x%d0**2))*x%d1(ix)*
     +          x%d2(ix2*(ix2-1)/2+ix3)
     +        - (x%d0**y%d0)*(y%d0/(x%d0**2))*x%d1(ix2)*
     +          x%d2(ix*(ix-1)/2+ix3)
     +        - (x%d0**y%d0)*(y%d0/(x%d0**2))*x%d1(ix3)*
     +          x%d2(ix*(ix-1)/2+ix2)
c
c             d1*d2 needs permutation of arguments here
c
     +        + (x%d0**y%d0)*((y%d0/x%d0)**2)*x%d1(ix)*
     +          x%d2(ix2*(ix2-1)/2+ix3)
     +        + (x%d0**y%d0)*((y%d0/x%d0)**2)*x%d1(ix2)*
     +          x%d2(ix*(ix-1)/2+ix3)
     +        + (x%d0**y%d0)*((y%d0/x%d0)**2)*x%d1(ix3)*
     +          x%d2(ix*(ix-1)/2+ix2)
c
     +        + 2*(x%d0**y%d0)*(y%d0/(x%d0**3))*x%d1(ix)*x%d1(ix2)*
     +          x%d1(ix3)
     +        - 3*(x%d0**y%d0)*(y%d0/x%d0)*(y%d0/x%d0**2)*x%d1(ix)*
     +          x%d1(ix2)*x%d1(ix3)
     +        + (x%d0**y%d0)*((y%d0/x%d0)**3)*x%d1(ix)*
     +          x%d1(ix2)*x%d1(ix3)
            enddo
          enddo
        enddo
c
c       Dxxy
c
        do ix = 1, x%nvar
          do ix2 = 1, ix
            iy = 1
            do while (iyto(iy).le.ixto(ix2).and.iy.le.y%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               ixto(ix2)*(ixto(ix2)-1)/2+iyto(iy))
     +        + (x%d0**y%d0)/x%d0*x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*
     +          x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        - (x%d0**y%d0)/(x%d0**2)*x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        - (x%d0**y%d0)*log(x%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*((y%d0/x%d0)**2)*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + 2*(x%d0**y%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
              iy = iy + 1
            enddo
          enddo
        enddo
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            ix2 = 1
            do while (ixto(ix2).le.iyto(iy).and.ix2.le.x%nvar)
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+ixto(ix2))
     +        + (x%d0**y%d0)/x%d0*x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*
     +          x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        - (x%d0**y%d0)/(x%d0**2)*x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        - (x%d0**y%d0)*log(x%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*((y%d0/x%d0)**2)*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + 2*(x%d0**y%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
              ix2 = ix2 + 1
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            do ix2 = 1, ix
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+ixto(ix2))
     +        + (x%d0**y%d0)/x%d0*x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*
     +          x%d2(ix*(ix-1)/2+ix2)*y%d1(iy)
     +        - (x%d0**y%d0)/(x%d0**2)*x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        - (x%d0**y%d0)*log(x%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + (x%d0**y%d0)*log(x%d0)*((y%d0/x%d0)**2)*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
     +        + 2*(x%d0**y%d0)*(y%d0/(x%d0**2))*
     +          x%d1(ix)*x%d1(ix2)*y%d1(iy)
            enddo
            ix = ix + 1
          enddo
        enddo
c
c       Dxyy
c
        do ix = 1, x%nvar
          iy = 1
          do while (iyto(iy).le.ixto(ix).and.iy.le.y%nvar)
            do iy2 = 1, iy
              s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +             iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        = s%d3((ixto(ix)+1)*ixto(ix)*(ixto(ix)-1)/6+
     +               iyto(iy)*(iyto(iy)-1)/2+iyto(iy2))
     +        + (x%d0**y%d0)/x%d0*x%d1(ix)*y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*(log(x%d0)**2)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
     +        + 2*(x%d0**y%d0)*(log(x%d0)/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
            enddo
            iy = iy + 1
          enddo
        enddo
        do iy = 1, y%nvar
          ix = 1
          do while (ixto(ix).le.iyto(iy).and.ix.le.x%nvar)
            iy2 = 1
            do while (iyto(iy2).le.ixto(ix).and.iy2.le.y%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               ixto(ix)*(ixto(ix)-1)/2+iyto(iy2))
     +        + (x%d0**y%d0)/x%d0*x%d1(ix)*y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*(log(x%d0)**2)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
     +        + 2*(x%d0**y%d0)*(log(x%d0)/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
              iy2 = iy2 + 1
            enddo
            ix = ix + 1
          enddo
        enddo
        do iy = 1, y%nvar
          do iy2 = 1, iy
            ix = 1
            do while (ixto(ix).le.iyto(iy2).and.ix.le.x%nvar)
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+ixto(ix))
     +        + (x%d0**y%d0)/x%d0*x%d1(ix)*y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*log(x%d0)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d2(iy*(iy-1)/2+iy2)
     +        + (x%d0**y%d0)*(log(x%d0)**2)*(y%d0/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
     +        + 2*(x%d0**y%d0)*(log(x%d0)/x%d0)*x%d1(ix)*
     +          y%d1(iy)*y%d1(iy2)
              ix = ix + 1
            enddo
          enddo
        enddo
c
c       Dyyy
c
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +             iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        = s%d3((iyto(iy)+1)*iyto(iy)*(iyto(iy)-1)/6+
     +               iyto(iy2)*(iyto(iy2)-1)/2+iyto(iy3))
     +        + (x%d0**y%d0)*log(x%d0)*
     +          y%d3((iy+1)*iy*(iy-1)/6+iy2*(iy2-1)/2+iy3)
c
c             d1*d2 needs permutation of arguments here
c
     +        + (x%d0**y%d0)*(log(x%d0)**2)*y%d1(iy)*
     +          y%d2(iy2*(iy2-1)/2+iy3)
     +        + (x%d0**y%d0)*(log(x%d0)**2)*y%d1(iy2)*
     +          y%d2(iy*(iy-1)/2+iy3)
     +        + (x%d0**y%d0)*(log(x%d0)**2)*y%d1(iy3)*
     +          y%d2(iy*(iy-1)/2+iy2)
c
     +        + (x%d0**y%d0)*(log(x%d0)**3)*y%d1(iy)*
     +          y%d1(iy2)*y%d1(iy3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_pow
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
      function wfn1_ad_dble_powx(x,y) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        double precision, intent(in) :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = x%d0**y
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)*y*x%d0**(y - 1.0d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2) 
     &      = x%d1(ix)*x%d1(ix2)*y*(y-1.0d0)*x%d0**(y-2.0d0)
     &      + x%d2(ix*(ix-1)/2+ix2)*y*x%d0**(y-1.0d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix-1)*ix*(ix+1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d1(ix)*x%d1(ix2)*x%d1(ix3)*y*(y-1.0d0)*(y-2.0d0)
     &          *x%d0**(y-3.0d0)
c
c             d1*d2 needs permutation of arguments here
c
     &        + x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*y*(y-1.0d0)
     &          *x%d0**(y - 2.0d0)
     &        + x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*y*(y-1.0d0)
     &          *x%d0**(y - 2.0d0)
     &        + x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*y*(y-1.0d0)
     &          *x%d0**(y - 2.0d0)
c
     &        + x%d3((ix-1)*ix*(ix+1)/6+ix2*(ix2-1)/2+ix3)*y
     &          *x%d0**(y-1.0d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_powx
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
      function wfn1_ad_dble_powy(x,y) result (s)
        double precision, intent(in) :: x
        type(wfn1_ad_dble), intent(in)  :: y
        type(wfn1_ad_dble)              :: s
        integer                      :: iy, iy2, iy3
        s%d0 = x**y%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,y)
        do iy = 1, y%nvar
          s%d1(iy) = x**y%d0*log(x)*y%d1(iy)
        enddo
#endif
#if NWAD_DEG >= 2
        do iy = 1, y%nvar
          do iy2 = 1, iy
            s%d2(iy*(iy-1)/2+iy2)
     &      = x**y%d0*(log(x)*y%d2(iy*(iy-1)/2+iy2) 
     &      + log(x)**2*y%d1(iy)*y%d1(iy2))
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do iy = 1, y%nvar
          do iy2 = 1, iy
            do iy3 = 1, iy2
              s%d3((iy-1)*iy*(iy+1)/6+iy2*(iy2-1)/2+iy3)
     &        = x**y%d0*(
     &          log(x)*y%d3((iy-1)*iy*(iy+1)/6+iy2*(iy2-1)/2+iy3)
c
c             d1*d2 needs permutation of arguments here
c
     &        + log(x)**2*y%d1(iy)*y%d2(iy2*(iy2-1)/2+iy3)
     &        + log(x)**2*y%d1(iy2)*y%d2(iy*(iy-1)/2+iy3)
     &        + log(x)**2*y%d1(iy3)*y%d2(iy*(iy-1)/2+iy2)
c
     &        + log(x)**3*y%d1(iy)*y%d1(iy2)*y%d1(iy3))
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_powy
      !>
      !> \brief Return whether \f$x\f$ equals \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
      !> context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_equal(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                     :: s
        s = x%d0 .eq. y%d0
      end function wfn1_ad_dble_equal
      !>
      !> \brief Return whether \f$x\f$ equals \f$y\f$ where the latter is inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
      !> context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$y\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_equalx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = x%d0 .eq. y
      end function wfn1_ad_dble_equalx
      !>
      !> \brief Return whether \f$x\f$ equals \f$y\f$ where the former is inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ equals \f$y\f$. In 
      !> context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$x\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_equaly(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x .eq. y%d0
      end function wfn1_ad_dble_equaly
      !>
      !> \brief Return whether \f$x\f$ does not equal \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_notequal(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = .not. (x .eq. y)
      end function wfn1_ad_dble_notequal
      !>
      !> \brief Return whether \f$x\f$ does not equal \f$y\f$ where the latter is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$y\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_notequalx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = .not. (x .eq. y)
      end function wfn1_ad_dble_notequalx
      !>
      !> \brief Return whether \f$x\f$ does not equal \f$y\f$ where the former is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ does not equal \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$x\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_notequaly(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = .not. (x .eq. y)
      end function wfn1_ad_dble_notequaly
      !>
      !> \brief Return whether \f$x\f$ is less than \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_lessthan(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x%d0 .lt. y%d0
      end function wfn1_ad_dble_lessthan
      !>
      !> \brief Return whether \f$x\f$ is less than \f$y\f$ where the latter is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$y\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_lessthanx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = x%d0 .lt. y
      end function wfn1_ad_dble_lessthanx
      !>
      !> \brief Return whether \f$x\f$ is less than \f$y\f$ where the former is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ is less than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$x\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_lessthany(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x .lt. y%d0
      end function wfn1_ad_dble_lessthany
      !>
      !> \brief Return whether \f$x\f$ is less than or equal to \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ is less than or equal
      !> to\f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_lessequal(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x%d0 .le. y%d0
      end function wfn1_ad_dble_lessequal
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
      function wfn1_ad_dble_lessequalx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = x%d0 .le. y
      end function wfn1_ad_dble_lessequalx
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
      function wfn1_ad_dble_lessequaly(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x .le. y%d0
      end function wfn1_ad_dble_lessequaly
      !>
      !> \brief Return whether \f$x\f$ is greater than \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_greaterthan(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x%d0 .gt. y%d0
      end function wfn1_ad_dble_greaterthan
      !>
      !> \brief Return whether \f$x\f$ is greater than \f$y\f$ where the latter is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$y\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_greaterthanx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = x%d0 .gt. y
      end function wfn1_ad_dble_greaterthanx
      !>
      !> \brief Return whether \f$x\f$ is greater than \f$y\f$ where the former is
      !> inactive
      !>
      !> Return a logical value reflecting whether \f$x\f$ is greater than \f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives. In this function \f$x\f$
      !> is an inactive variable.
      !>
      function wfn1_ad_dble_greaterthany(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x .gt. y%d0
      end function wfn1_ad_dble_greaterthany
      !>
      !> \brief Return whether \f$x\f$ is greater than or equal to \f$y\f$
      !>
      !> Return a logical value reflecting whether \f$x\f$ is greater than or equal
      !> to\f$y\f$.
      !> In context of derivative calculations the comparison only applies to the
      !> the values and not to any of the derivatives.
      !>
      function wfn1_ad_dble_greaterequal(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x%d0 .ge. y%d0
      end function wfn1_ad_dble_greaterequal
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
      function wfn1_ad_dble_greaterequalx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        logical                        :: s
        s = x%d0 .ge. y
      end function wfn1_ad_dble_greaterequalx
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
      function wfn1_ad_dble_greaterequaly(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        logical                        :: s
        s = x .ge. y%d0
      end function wfn1_ad_dble_greaterequaly
      !>
      !> \brief Evaluate the sign function 
      !> 
      !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
      !> the sign of \f$y\f$. This routine implements this function for the case
      !> where both \f$x\f$ and \f$y\f$ are active variables.
      !>
      function wfn1_ad_dble_sign(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        type(wfn1_ad_dble)             :: s
        s = abs(x) * sign(1.0d0,y%d0)
      end function wfn1_ad_dble_sign
      !>
      !> \brief Evaluate the sign function where \f$y\f$ is inactive
      !> 
      !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
      !> the sign of \f$y\f$. This routine implements this function for the case
      !> where \f$x\f$ is an active and \f$y\f$ is an inactive variable.
      !>
      function wfn1_ad_dble_signx(x,y) result(s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision,   intent(in) :: y
        type(wfn1_ad_dble)             :: s
        s = abs(x) * sign(1.0d0,y)
      end function wfn1_ad_dble_signx
      !>
      !> \brief Evaluate the sign function where \f$x\f$ is inactive
      !> 
      !> The function \f$\mathrm{sign}(x,y)\f$ returns the value of \f$x\f$ with
      !> the sign of \f$y\f$. This routine implements this function for the case
      !> where \f$x\f$ is an inactive and \f$y\f$ is an active variable.
      !>
      function wfn1_ad_dble_signy(x,y) result(s)
        double precision,   intent(in) :: x
        type(wfn1_ad_dble), intent(in) :: y
        double precision               :: s
        s = abs(x) * sign(1.0d0,y%d0)
      end function wfn1_ad_dble_signy
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
      function wfn1_ad_dble_abs(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix
        s%d0 = abs(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        if (x%d0.lt.0.0d0) then
          do ix = 1, x%nvar
            s%d1(ix) = -x%d1(ix)
          enddo
#if NWAD_DEG >= 2
          do ix = 1, x%nvar*(x%nvar+1)/2
            s%d2(ix) = -x%d2(ix)
          enddo
#endif
#if NWAD_DEG >= 3
          do ix = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
            s%d3(ix) = -x%d3(ix)
          enddo
#endif
        else if (x%d0.gt.0.0d0) then
          do ix = 1, x%nvar
            s%d1(ix) = x%d1(ix)
          enddo
#if NWAD_DEG >= 2
          do ix = 1, x%nvar*(x%nvar+1)/2
            s%d2(ix) = x%d2(ix)
          enddo
#endif
#if NWAD_DEG >= 3
          do ix = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
            s%d3(ix) = x%d3(ix)
          enddo
#endif
        else
          do ix = 1, x%nvar
            s%d1(ix) = 0.0d0
          enddo
#if NWAD_DEG >= 2
          do ix = 1, x%nvar*(x%nvar+1)/2
            s%d2(ix) = 0.0d0
          enddo
#endif
#if NWAD_DEG >= 3
          do ix = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
            s%d3(ix) = 0.0d0
          enddo
#endif
        endif
#endif
      end function wfn1_ad_dble_abs
      !>
      !> \brief Evaluate the unary negation operator \f$-\f$
      !>
      !> The unary negation operator simply replaces the value and all the
      !> derivatives with the same with the opposite sign.
      !>
      function wfn1_ad_dble_minus(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix
        s%d0 = -x%d0
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = -x%d1(ix)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar*(x%nvar+1)/2
          s%d2(ix) = -x%d2(ix)
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar*(x%nvar+1)*(x%nvar+2)/6
          s%d3(ix) = -x%d3(ix)
        enddo
#endif
      end function wfn1_ad_dble_minus
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
      function wfn1_ad_dble_sqrt(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = sqrt(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = 0.5d0/sqrt(x%d0)*x%d1(ix)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = 0.5d0/sqrt(x%d0)*x%d2(ix*(ix-1)/2+ix2)
     &      - 0.25d0/(x%d0*sqrt(x%d0))*x%d1(ix)*x%d1(ix2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = 0.5d0/sqrt(x%d0)
     &          *x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
c
c             d1*d2 needs permutation of arguments here
c
     &        - 0.25d0/(x%d0*sqrt(x%d0))
     &          *x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &        - 0.25d0/(x%d0*sqrt(x%d0))
     &          *x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &        - 0.25d0/(x%d0*sqrt(x%d0))
     &          *x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
c
     &        + 3.0d0/(8.0d0*x%d0**2*sqrt(x%d0))
     &          *x%d1(ix)*x%d1(ix2)*x%d1(ix3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_sqrt
      !>
      !> \brief Evaluate the \f$\exp\f$ function
      !>
      !> The implementation of the \f$\exp\f$ function. The chain rule is used
      !> to evaluate the derivatives. I.e. we consider \f$s(r) = \exp(x(r))\f$.
      !>
      function wfn1_ad_dble_exp(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = exp(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)*exp(x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)*exp(x%d0)
     &      + x%d1(ix)*x%d1(ix2)*exp(x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)*exp(x%d0)
c
c             d1*d2 requires permutation of arguments here
c
     &        + x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*exp(x%d0)
     &        + x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*exp(x%d0)
     &        + x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*exp(x%d0)
c
     &        + x%d1(ix)*x%d1(ix2)*x%d1(ix3)*exp(x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_exp
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
      function wfn1_ad_dble_log(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = log(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/x%d0
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/x%d0
     &      - (x%d1(ix)/x%d0)*(x%d1(ix2)/x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)/x%d0
c
c             d1*d2 requires permutation of arguments here
c
     &        - (x%d1(ix)/x%d0)*(x%d2(ix2*(ix2-1)/2+ix3)/x%d0)
     &        - (x%d1(ix2)/x%d0)*(x%d2(ix*(ix-1)/2+ix3)/x%d0)
     &        - (x%d1(ix3)/x%d0)*(x%d2(ix*(ix-1)/2+ix2)/x%d0)
c
     &        + 2*(x%d1(ix)/x%d0)*(x%d1(ix2)/x%d0)*(x%d1(ix3)/x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_log
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
      function wfn1_ad_dble_sin(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = sin(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)*cos(x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)*cos(x%d0)
     &      - x%d1(ix)*x%d1(ix2)*sin(x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)*cos(x%d0)
c
c             d1*d2 requires permutation of arguments here
c
     &        - x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*sin(x%d0)
     &        - x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*sin(x%d0)
     &        - x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*sin(x%d0)
c
     &        - x%d1(ix)*x%d1(ix2)*x%d1(ix3)*cos(x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_sin
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
      function wfn1_ad_dble_cos(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = cos(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = -x%d1(ix)*sin(x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = -x%d2(ix*(ix-1)/2+ix2)*sin(x%d0)
     &      -  x%d1(ix)*x%d1(ix2)*cos(x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = -x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)*sin(x%d0)
c
c             d1*d2 requires permutation of arguments here
c
     &        - x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*cos(x%d0)
     &        - x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*cos(x%d0)
     &        - x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*cos(x%d0)
c
     &        + x%d1(ix)*x%d1(ix2)*x%d1(ix3)*sin(x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_cos
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
      function wfn1_ad_dble_tan(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = tan(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/cos(x%d0)**2
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/cos(x%d0)**2
     &      + x%d1(ix)*x%d1(ix2)*(2*tan(x%d0)/cos(x%d0)**2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)/cos(x%d0)**2
c
c             d1*d2 requires permutations of all arguments
c
     &        + 2*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &          *tan(x%d0)/cos(x%d0)**2 
     &        + 2*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &          *tan(x%d0)/cos(x%d0)**2 
     &        + 2*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
     &          *tan(x%d0)/cos(x%d0)**2 
c
     &        + 4*x%d1(ix)*x%d1(ix2)*x%d1(ix3)*tan(x%d0)**2/cos(x%d0)**2
     &        + 2*x%d1(ix)*x%d1(ix2)*x%d1(ix3)/cos(x%d0)**4
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_tan
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
      function wfn1_ad_dble_sinh(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = sinh(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)*cosh(x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)*cosh(x%d0)
     &      + x%d1(ix)*x%d1(ix2)*sinh(x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)*cosh(x%d0)
c
c             d1*d2 requires argument permutations
c
     &        + x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*sinh(x%d0)
     &        + x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*sinh(x%d0)
     &        + x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*sinh(x%d0)
c
     &        + x%d1(ix)*x%d1(ix2)*x%d1(ix3)*cosh(x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_sinh
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
      function wfn1_ad_dble_cosh(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = cosh(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)*sinh(x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)*sinh(x%d0)
     &      + x%d1(ix)*x%d1(ix2)*cosh(x%d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
         do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)*sinh(x%d0)
c
c             d1*d2 requires argument permutations
c
     &        + x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*cosh(x%d0)
     &        + x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*cosh(x%d0)
     &        + x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*cosh(x%d0)
c
     &        + x%d1(ix)*x%d1(ix2)*x%d1(ix3)*sinh(x%d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_cosh
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
      function wfn1_ad_dble_tanh(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = tanh(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/cosh(x%d0)**2
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/cosh(x%d0)**2
     &      - x%d1(ix)*x%d1(ix2)*(2*tanh(x%d0)/cosh(x%d0)**2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)/cosh(x%d0)**2
c
c             d1*d2 requires argument permutations
c
     &        - 2*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &          *tanh(x%d0)/cosh(x%d0)**2
     &        - 2*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &          *tanh(x%d0)/cosh(x%d0)**2
     &        - 2*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
     &          *tanh(x%d0)/cosh(x%d0)**2
c
     &        + 4*x%d1(ix)*x%d1(ix2)*x%d1(ix3)
     &          *tanh(x%d0)**2/cosh(x%d0)**2
     &        - 2*x%d1(ix)*x%d1(ix2)*x%d1(ix3)/cosh(x%d0)**4
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_tanh
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
      function wfn1_ad_dble_asin(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = asin(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) 
     &      + x%d0*x%d1(ix)*x%d1(ix2)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &          /((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0)) 
c
c             d1*d2 requires argument permutations
c
     &        + x%d0*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) 
     &        + x%d0*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) 
     &        + x%d0*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) 
c
     &        + (x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0)) 
     &        + 3.0d0*(x%d0**2.0d0)*(x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0-x%d0*x%d0)**(5.0d0/2.0d0))
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_asin
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
      function wfn1_ad_dble_acos(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = acos(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = -x%d1(ix)/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = -x%d2(ix*(ix-1)/2+ix2)/((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
     &      - x%d0*x%d1(ix)*x%d1(ix2)/((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = -x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &           /((1.0d0-x%d0*x%d0)**(1.0d0/2.0d0))
c
c             d1*d2 requires argument permutations
c
     &        - x%d0*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
     &        - x%d0*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
     &        - x%d0*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
c
     &        - (x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0-x%d0*x%d0)**(3.0d0/2.0d0))
     &        - 3.0d0*(x%d0**2.0d0)*(x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0-x%d0*x%d0)**(5.0d0/2.0d0))
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_acos
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
      function wfn1_ad_dble_atan(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        integer                      :: ix, ix2, ix3
        s%d0 = atan(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/(1.0d0+x%d0*x%d0)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/(1.0d0+x%d0*x%d0) 
     &      - 2.0d0*x%d0*(x%d1(ix)*x%d1(ix2))/((1.0d0+x%d0*x%d0)**2.0d0)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &          /(1.0d0+x%d0*x%d0) 
c
c             d1*d2 requires argument permutations
c
     &        - 2.0d0*x%d0*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &          /((1.0d0+x%d0*x%d0)**2.0d0) 
     &        - 2.0d0*x%d0*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &          /((1.0d0+x%d0*x%d0)**2.0d0) 
     &        - 2.0d0*x%d0*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
     &          /((1.0d0+x%d0*x%d0)**2.0d0) 
c
     &        - 2.0d0*(x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0+x%d0*x%d0)**2.0d0) 
     &        + 8.0d0*(x%d0**2.0d0)*(x%d1(ix)*x%d1(ix2)*x%d1(ix3))
     &          /((1.0d0+x%d0*x%d0)**3.0d0)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_atan
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
      function wfn1_ad_dble_asinh(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        double precision             :: t1, t12
        integer                      :: ix, ix2, ix3
        t1 = 1.0d0 + x%d0*x%d0
        t12 = sqrt(t1)
        s%d0 = log(x%d0+t12)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = x%d1(ix)/t12
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = x%d2(ix*(ix-1)/2+ix2)/t12
     &      - x%d1(ix)*x%d1(ix2)*x%d0/(t1*t12)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)/t12
c
c             d1*d2 requires argument permutations
c
     &        - x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)*x%d0/(t1*t12)
     &        - x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)*x%d0/(t1*t12)
     &        - x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)*x%d0/(t1*t12)
c
     &        - x%d1(ix)*x%d1(ix2)*x%d1(ix3)/(t1*t12) 
     &        + 3*x%d1(ix)*x%d1(ix2)*x%d1(ix3)*x%d0**2/(t1*t1*t12)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_asinh
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
      function wfn1_ad_dble_erf(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        double precision             :: t1
        integer                      :: ix, ix2, ix3
        t1 = exp(-(x%d0**2))/sqrt(acos(-1.0d0))
        s%d0 = erf(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = 2*t1*x%d1(ix)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = 2*t1*x%d2(ix*(ix-1)/2+ix2)
     &      - 4*x%d0*t1*x%d1(ix)*x%d1(ix2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = 2*t1*x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
c
c             d1*d2 requires argument permutations
c
     &        - 4*x%d0*t1*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &        - 4*x%d0*t1*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &        - 4*x%d0*t1*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
c
     &        + 8*x%d0**2*t1*x%d1(ix)*x%d1(ix2)*x%d1(ix3)
     &        - 4*t1*x%d1(ix)*x%d1(ix2)*x%d1(ix3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_erf
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
      function wfn1_ad_dble_erfc(x) result (s)
        type(wfn1_ad_dble), intent(in)  :: x
        type(wfn1_ad_dble)              :: s
        double precision             :: t1
        integer                      :: ix, ix2, ix3
        t1 = exp(-(x%d0**2))/sqrt(acos(-1.0d0))
        s%d0 = 1.0d0-erf(x%d0)
#if NWAD_DEG >= 1
        call wfn1_ad_copy_dependencies(s,x)
        do ix = 1, x%nvar
          s%d1(ix) = -2*t1*x%d1(ix)
        enddo
#endif
#if NWAD_DEG >= 2
        do ix = 1, x%nvar
          do ix2 = 1, ix
            s%d2(ix*(ix-1)/2+ix2)
     &      = -2*t1*x%d2(ix*(ix-1)/2+ix2)
     &      + 4*x%d0*t1*x%d1(ix)*x%d1(ix2)
          enddo
        enddo
#endif
#if NWAD_DEG >= 3
        do ix = 1, x%nvar
          do ix2 = 1, ix
            do ix3 = 1, ix2
              s%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
     &        = -2*t1*x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
c
c             d1*d2 requires argument permutations
c
     &        + 4*x%d0*t1*x%d1(ix)*x%d2(ix2*(ix2-1)/2+ix3)
     &        + 4*x%d0*t1*x%d1(ix2)*x%d2(ix*(ix-1)/2+ix3)
     &        + 4*x%d0*t1*x%d1(ix3)*x%d2(ix*(ix-1)/2+ix2)
c
     &        - 8*x%d0**2*t1*x%d1(ix)*x%d1(ix2)*x%d1(ix3)
     &        + 4*t1*x%d1(ix)*x%d1(ix2)*x%d1(ix3)
            enddo
          enddo
        enddo
#endif
      end function wfn1_ad_dble_erfc
      !>
      !> \brief Return the number of variables
      !>
      !> Look up and return the number of variables that the
      !> active variable depends on.
      !>
      function wfn1_ad_dble_get_nvar(x) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        integer                     :: s
        s = x%nvar
      end function wfn1_ad_dble_get_nvar
      !>
      !> \brief Initialize an active variable to \f$\rho_\alpha\f$
      !>
      !> Initialize an active variable and set it to \f$\rho_\alpha\f$.
      !> Active variables are those with respect
      !> to which the derivatives are calculated in the current evaluation of the
      !> code. In practice it means that the components are initialized as
      !> \f{eqnarray*}{
      !>   d0 &=& \frac{\mathrm{d}^0 x}{\mathrm{d}x^0} = x \\\\
      !>   d1 &=& \frac{\mathrm{d}^1 x}{\mathrm{d}x^1} = 1 \\\\
      !>   d2 &=& \frac{\mathrm{d}^2 x}{\mathrm{d}x^2} = 0 \\\\
      !>   d3 &=& \frac{\mathrm{d}^3 x}{\mathrm{d}x^3} = 0 \\\\
      !> \f}
      !> where there variable set only contains the current variable.
      !>
      function wfn1_ad_dble_dij(x,ndim,i,j) result (s)
#include "errquit.fh"
        double precision, intent(in) :: x
        type(wfn1_ad_dble)           :: s
        if (i*i>mxvar) then
          call errquit("wfn1_ad_dble_dij: i too large: ",i,UERR)
        endif
        if (j*j>mxvar) then
          call errquit("wfn1_ad_dble_dij: j too large: ",j,UERR)
        endif
        if (ndim*ndim>mxvar) then
          call errquit("wfn1_ad_dble_dij: ndim too large: ",ndim,UERR)
        endif
        s%nvar = 1
        s%ivar(1) = i+ndim*(j-1)
        s%d0 = x
#if NWAD_DEG >= 1
        s%d1(1) = 1
#endif
#if NWAD_DEG >= 2
        s%d2(1) = 0
#endif
#if NWAD_DEG >= 3
        s%d3(1) = 0
#endif
      end function wfn1_ad_dble_dij
      !>
      !> \brief Return the value 
      !>
      !> Simply return the value, i.e. 0th order derivative, of the quantity
      !>
      function wfn1_ad_dble_get_val(x) result (s)
        type(wfn1_ad_dble), intent(in) :: x
        double precision            :: s
        s = x%d0
      end function wfn1_ad_dble_get_val
      !>
      !> \brief Return the 1st derivative
      !>
      !> Look up and return the value of a specific 1st order derivative
      !> component. The argument 'ix' must be in the range 1 to x%nvar.
      !> On return 's' will hold the value of the derivative whereas
      !> 'iv' contains the number of the variable, i.e. x%ivar(ix).
      !>
      subroutine wfn1_ad_dble_get_d1(x,ix,s,iv)
        type(wfn1_ad_dble),  intent(in)  :: x
        integer,          intent(in)  :: ix
        integer,          intent(out) :: iv
        double precision, intent(out) :: s
#ifndef NODEBUG
        if (ix.lt.1.or.ix.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d1: ix out of range",ix,UERR)
        endif
#endif
#if NWAD_DEG >= 1
        s  = x%d1(ix)
        iv = x%ivar(ix)
#else
        s  = 0.0d0
        iv = -1
        call errquit("wfn1_ad_dble_get_d1: 1st order derivatives not "//
     +               "supported",NWAD_DEG,UERR)
#endif
      end subroutine wfn1_ad_dble_get_d1
      !>
      !> \brief Return the 2nd derivative
      !>
      !> Look up and return the value of a specific 2nd order derivative
      !> component. The arguments 'ix' and 'ix2' must be in the range
      !> 1 to x%nvar. Also 'ix2' must be less than or equal to 'ix'.
      !> On return 's' will hold the value of the derivative, whereas
      !> 'iv' and 'iv2' hold the numbers of the corresponding variables,
      !> i.e. 'iv' contains 'x%ivar(ix)' and 'iv2' contains 'x%ivar(ix2)'.
      !>
      subroutine wfn1_ad_dble_get_d2(x,ix,ix2,s,iv,iv2)
        type(wfn1_ad_dble),  intent(in)  :: x
        integer,          intent(in)  :: ix
        integer,          intent(in)  :: ix2
        double precision, intent(out) :: s
        integer,          intent(out) :: iv
        integer,          intent(out) :: iv2
#ifndef NODEBUG
        if (ix.lt.1.or.ix.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d2: ix out of range",ix,UERR)
        endif
        if (ix2.lt.1.or.ix2.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d2: ix2 out of range",ix2,UERR)
        endif
        if (ix2.gt.ix) then
          call errquit("wfn1_ad_dble_get_d2: ix2 greater than ix",
     $                 ix2,UERR)
        endif
#endif
#if NWAD_DEG >= 2
        s   = x%d2(ix*(ix-1)/2+ix2)
        iv  = x%ivar(ix)
        iv2 = x%ivar(ix2)
#else
        s   = 0.0d0
        iv  = -1
        iv2 = -1
        call errquit("wfn1_ad_dble_get_d2: 2nd order derivatives not "//
     +               "supported",NWAD_DEG,UERR)
#endif
      end subroutine wfn1_ad_dble_get_d2
      !>
      !> \brief Return the 3rd derivative
      !>
      !> Look up and return the value of a specific 3rd order derivative
      !> component. The arguments 'ix', 'ix2', and 'ix3' must be in the range
      !> 1 to x%nvar. Also 'ix2' must be less than or equal to 'ix', and
      !> 'ix3' must be less than or equal to 'ix2'.
      !> On return 's' will hold the value of the derivative, whereas
      !> 'iv', 'iv2' and 'iv3' hold the numbers of the corresponding variables,
      !> i.e. 'iv' contains 'x%ivar(ix)', 'iv2' contains 'x%ivar(ix2)', and
      !> 'iv3' contains 'x%ivar(ix3)'.
      !>
      subroutine wfn1_ad_dble_get_d3(x,ix,ix2,ix3,s,iv,iv2,iv3)
        type(wfn1_ad_dble),  intent(in)  :: x
        integer,          intent(in)  :: ix
        integer,          intent(in)  :: ix2
        integer,          intent(in)  :: ix3
        double precision, intent(out) :: s
        integer,          intent(out) :: iv
        integer,          intent(out) :: iv2
        integer,          intent(out) :: iv3
#ifndef NODEBUG
        if (ix.lt.1.or.ix.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d3: ix out of range",ix,UERR)
        endif
        if (ix2.lt.1.or.ix2.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d3: ix2 out of range",ix2,UERR)
        endif
        if (ix3.lt.1.or.ix3.gt.x%nvar) then
          call errquit("wfn1_ad_dble_get_d3: ix3 out of range",ix3,UERR)
        endif
        if (ix2.gt.ix) then
          call errquit("wfn1_ad_dble_get_d3: ix2 greater than ix",
     +                 ix2,UERR)
        endif
        if (ix3.gt.ix2) then
          call errquit("wfn1_ad_dble_get_d3: ix3 greater than ix2",
     &                 ix3,UERR)
        endif
#endif
#if NWAD_DEG >= 3
        s   = x%d3((ix+1)*ix*(ix-1)/6+ix2*(ix2-1)/2+ix3)
        iv  = x%ivar(ix)
        iv2 = x%ivar(ix2)
        iv3 = x%ivar(ix3)
#else
        s   = 0.0d0
        iv  = -1
        iv2 = -1
        iv3 = -1
        call errquit("wfn1_ad_dble_get_d3: 3rd order derivatives not "//
     +               "supported",NWAD_DEG,UERR)
#endif
      end subroutine wfn1_ad_dble_get_d3
      !>
      !> \brief Print the value of a NWAD_DBLE variable
      !>
      !> Generate formated output on the specified unit number, detailing
      !> the contents of a NWAD_DBLE variable.
      !>
      subroutine wfn1_ad_dble_print(iunit,x)
        integer,         intent(in) :: iunit
        type(wfn1_ad_dble), intent(in) :: x
        integer                     :: i
        write(iunit,'("vv = ",7i3)')(x%ivar(i),i=1,x%nvar)
        write(iunit,'("d0 = ",e14.6)')x%d0
#if NWAD_DEG >= 1
        write(iunit,'("d1 = ",7e14.6)')(x%d1(i),i=1,x%nvar)
#endif
#if NWAD_DEG >= 2
        write(iunit,'("d2 = ",28e14.6)')
     &       (x%d2(i),i=1,x%nvar*(x%nvar+1)/2)
#endif
#if NWAD_DEG >= 3
        write(iunit,'("d3 = ",84e14.6)')
     &       (x%d3(i),i=1,x%nvar*(x%nvar+1)*(x%nvar+2)/6)
#endif
        write(iunit,*)
      end subroutine wfn1_ad_dble_print
      subroutine wfn1_ad_test(x)
         type(wfn1_ad_dble), intent(in) :: x
         integer i
         do i = 1, x%nvar
#if NWAD_DEG >= 1
           if (abs(x%d1(i)).gt.1.0d6) then
             call wfn1_ad_stop()
           endif
#endif
         enddo
      end subroutine wfn1_ad_test
      subroutine wfn1_ad_stop()
        write(0,*)
      end subroutine wfn1_ad_stop
#if NWAD_DEG <= 0
      end module wfn1_ad0
#elif NWAD_DEG <= 1
      end module wfn1_ad1
#elif NWAD_DEG <= 2
      end module wfn1_ad2
#else
      end module wfn1_ad3
#endif
#if NWAD_DEG < 1
#undef NWAD_DEG
#define NWAD_DEG 1
#include "wfn1_ad.F"
#endif
#if NWAD_DEG < 2
#undef NWAD_DEG
#define NWAD_DEG 2
#include "wfn1_ad.F"
#endif
#if NWAD_DEG < 3
#undef NWAD_DEG
#define NWAD_DEG 3
#include "wfn1_ad.F"
#endif
!> @}
