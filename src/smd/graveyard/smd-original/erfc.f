      real*8 FUNCTION erfxc(x)

      implicit none

      real*8 a1,a2,a3,a4,a5,p

      parameter ( a1 = 0.254829592, a2 = -0.284496736 )
      parameter ( a3 = 1.421413741, a4 = -1.453152027 )
      parameter ( a5 = 1.061405429, p  =  0.327591100 )

      real*8 t,x,xsq,tp

      t=1.0/(1.0+p*x)
      xsq=x*x

      tp=t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))

      erfxc=tp*exp(-xsq)

      return

      END
c $Id$
