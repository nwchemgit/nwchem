C> \ingroup wfn1
C> @{
C>
C> \brief Fit a parabola given three points
C>
C> Given \f$ f(x_1) \f$, \f$ f(x_2) \f$, and \f$ f(x_3) \f$
C> find the parabola \f$ f(x)=ax^2+bx+c \f$ matching this
C> data, and predict its minimum. We assume that
C> \f$ x_1 \lt x_2 \lt x_3 \f$.
C>
      subroutine wfn1_f1f2f3(x1,x2,x3,f1,f2,f3,t,a,b,c,xm,fxm)
      implicit none
c
      double precision, intent(in)  :: x1  !< \f$x_1\f$
      double precision, intent(in)  :: x2  !< \f$x_2\f$
      double precision, intent(in)  :: x3  !< \f$x_3\f$
      double precision, intent(in)  :: f1  !< \f$f(x_1)\f$
      double precision, intent(in)  :: f2  !< \f$f(x_2)\f$
      double precision, intent(in)  :: f3  !< \f$f(x_3)\f$
      double precision, intent(in)  :: t   !< The trust region ensuring
                                           !< that \f$-t \le xm \le t\f$
      double precision, intent(out) :: a   !< The coefficient of \f$x^2\f$
      double precision, intent(out) :: b   !< The coefficient of \f$x\f$
      double precision, intent(out) :: c   !< The coefficient of \f$x^0\f$
      double precision, intent(out) :: xm  !< The value of \f$x\f$ that
                                           !< minimizes \f$f(x)\f$
      double precision, intent(out) :: fxm !< The value of \f$f(xm)\f$
c
      double precision xn, xu
      double precision f
      f(xm) = a*xm*xm+b*xm+c
c
      c = f2
      xn = x1 - x2
      xu = x3 - x2
      if (abs(xu).lt.1.0d-8.or.abs(xn).lt.1.0d-8) then
        xm = -1.0d0
        return
      endif
      b = (xu*xu*f1-xn*xn*f3-c*(xu*xu-xn*xn))/(xn*xu*xu-xu*xn*xn)
      a = (xu*xu*f1+xn*xn*f3-b*(xu*xn*xn+xu*xu*xn)-c*(xu*xu+xn*xn))/
     +    (xu*xu*xn*xn*2.0d0)
      b = b - 2.0d0*x2*a
      c = c - x2*x2*a - b*x2
c
      if (a.gt.0.0d0) then
        xm = -b/(2.0d0*a)
      else
c
c       the function has no minimum 
c
        if (2.0d0*a*x2+b.lt.0.0d0) then
          xm = t
        else
          xm = -t
        endif
      endif
      fxm = f(xm)
c
      end

C> @}

