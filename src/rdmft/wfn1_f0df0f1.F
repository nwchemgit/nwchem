C> \ingroup wfn1
C> @{
C>
C> \brief Fit a parabola given two points and a gradient
C>
C> Given \f$(0)\f$, 
C> \f$\left.\frac{\partial f(x)}{\partial x}\right|_{x=0}\f$,
C> and \f$f(x_1)\f$ find the parabola \f$f(x)=ax^2+bx+c\f$ matching this
C> data, and predict its minimum.
C>
      subroutine wfn1_f0df0f1(x1,f0,df0,f1,t,a,b,c,xm,fxm)
      implicit none
c
      double precision, intent(in)  :: x1  !< \f$x_1\f$
      double precision, intent(in)  :: f0  !< \f$f(0)\f$
      double precision, intent(in)  :: df0 !< \f$\left.\frac{df(x)}{dx}\right|_{x=0}\f$
      double precision, intent(in)  :: f1  !< \f$f(x_1)\f$
      double precision, intent(in)  :: t   !< The trust region ensuring
                                           !< that \f$-t\le xm \le t\f$
      double precision, intent(out) :: a   !< The coefficient of \f$x^2\f$
      double precision, intent(out) :: b   !< The coefficient of \f$x\f$
      double precision, intent(out) :: c   !< The coefficient of \f$x^0\f$
      double precision, intent(out) :: xm  !< The value of \f$x\f$ that
                                           !< minimizes \f$f(x)\f$
      double precision, intent(out) :: fxm !< The value of \f$f(xm)\f$
c
      double precision f
      f(xm) = a*xm*xm+b*xm+c
c
      c = f0
      b = df0
      a = (f1-f0-x1*df0)/(x1*x1)
c
      if (a.gt.0.0d0) then
        xm = -b/(2.0d0*a)
      else
c
c       the function has no minimum 
c
        if (df0.lt.0.0d0) then
          xm = t
        else
          xm = -t
        endif
      endif
      fxm = f(xm)
c
      end

C> @}

