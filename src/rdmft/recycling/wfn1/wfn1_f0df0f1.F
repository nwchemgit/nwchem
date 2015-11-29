C> \ingroup wfn1
C> @{
C>
C> \brief Fit a parabola given two points and a gradient
C>
C> Given \f$(0)\f$, 
C> \f$\left.\frac{\partial f(x)}{\partial x}\right|_{x=0}\f$,
C> and \f$f(x)\f$ find the parabola \f$f(x)=ax^2+bx+c\f$ matching this
C> data, and predict its minimum.
C>
      subroutine wfn1_f0df0f1(x,f0,df0,f1,t,a,b,c,xm,fxm)
      implicit none
c
      double precision x   !< [Input] \f$x\f$
      double precision f0  !< [Input] \f$f(0)\f$
      double precision df0 !< [Input] \f$\left.\frac{df(x)}{dx}\right|_{x=0}\f$
      double precision f1  !< [Input] \f$f(x)\f$
      double precision t   !< [Input] The trust region ensuring that
                           !< \f$ -t \le xm \le t\f$
      double precision a   !< [Output] The coefficient of \f$x^2\f$
      double precision b   !< [Output] The coefficient of \f$x\f$
      double precision c   !< [Output] The coefficient of \f$x^0\f$
      double precision xm  !< [Output] The value of \f$x\f$ that
                           !< minimizes \f$f(x)\f$
      double precision fxm !< [Output] The value of \f$f(xm)\f$
c
      double precision f
      f(xm) = a*xm*xm+b*xm+c
c
      c = f0
      b = df0
      a = (f1-f0-x*df0)/(x*x)
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

