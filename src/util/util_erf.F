#if  defined(HPUX)
      double precision function derf(x)
      implicit none
      double precision x,util_erf,derf
      external util_erf
      derf=util_erf(x)
      return
      end

      double precision function derfc(x)
      implicit none
      double precision x,util_erfc,derfc
      external util_erfc
      derfc=util_erfc(x)
      return
      end
c
      double precision function erfc(x)
      implicit none
      double precision x,util_erfc
      external util_erfc
      erfc=util_erfc(x)
      return
      end
c
#endif                                                                          

      real*8 function util_erf(x) 
      implicit none
      real*8 x
*     **** local variables ****
      real*8 f
*     ***** external functions ****
      real*8  gammp
      external gammp
      IF(x.lt.0.0d0)then
        f = -gammp(0.5d0,x**2)
      ELSE
        f = +gammp(0.5d0,x**2)
      ENDIF
      util_erf = f
      return
      end 

      real*8 function util_erfc(x) 
      implicit none
      real*8 x
*     **** local variables ****
      real*8 f
*     ***** external functions ****
      real*8  gammp
      external gammp
      IF(x.lt.0.0d0)then
        f = -gammp(0.5d0,x**2)
      ELSE
        f = +gammp(0.5d0,x**2)
      ENDIF
      util_erfc = 1.0d0-f
      return
      end 

c*************************************************************************
c !
c !      This is a slightly modified version of Log(Gamma) function
c !      program from Num. Rec.
c !************************************************************************

      real*8 function ln_gamma(xx)
      implicit none
      real*8  xx

*     *** local variables ****
      integer j
      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     >24.01409824083091d0,-1.231739572450155d0,0.1208650973866179d-2,
     >-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      DO 11 j=1,6
        y=y+1.0d0
        ser=ser+cof(j)/y
11    CONTINUE

      ln_gamma=tmp+dlog(stp*ser/x)

      return
      END

      real*8 function util_gamma(x)
      implicit none
      real*8 x

      real*8 xx
      real*8   ln_gamma
      external ln_gamma

      XX = X
      util_gamma = dexp(ln_gamma(xx))

      return
      END

      real*8 function util_gammp(a,x)
      implicit none
      real*8 a,x

#include "errquit.fh"

*     **** external functions ****
      real*8   gammcf,gamser,gln

      IF(x.LT.0.0d0 .OR. a.LE.0.0d0) THEN
        call errquit('bad arguments in util_gammp',0, INPUT_ERR)
      END IF

      IF (x .lt. (a+1.0d0)) THEN
        call gser(gamser,a,x,gln)
        util_gammp=gamser
      ELSE
        CALL gcf(gammcf,a,x,gln)
        util_gammp=1.0d0-gammcf
      ENDIF
      return
      end


     
      real*8 function gammp(a,x)
      implicit none
      real*8 a,x

#include "errquit.fh"

*     **** external functions ****
      real*8   gammcf,gamser,gln

      IF(x.LT.0.0d0 .OR. a.LE.0.0d0) THEN
        call errquit('bad arguments in gammp',0, INPUT_ERR)
      END IF

      IF (x .lt. (a+1.0d0)) THEN
        call gser(gamser,a,x,gln)
        gammp=gamser
      ELSE
        CALL gcf(gammcf,a,x,gln)
        gammp=1.0d0-gammcf
      ENDIF
      return
      end

      SUBROUTINE gcf(gammcf,a,x,gln)
      implicit none
#include "errquit.fh"
      integer ITMAX
      real*8  a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=1000,EPS=3.0d-12,FPMIN=1.0d-30)
      real*8       an,b,c,d,del,h
      double precision undovl
      parameter(undovl=-20d0*2.3025d0)
      integer i
      real*8   ln_gamma
      external ln_gamma

      gln=ln_gamma(a)
      b=x + 1.0d0 - a
      c=1.0d0/FPMIN
      d=1.0d0/b
      h=d
      DO 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.0d0
        d=an*d+b
        IF(DABS(d).lt.FPMIN) d=FPMIN
        c=b+an/c
        IF(DABS(c).lt.FPMIN) c=FPMIN
        d=1.0d0/d
        del=d*c
        h=h*del
        IF(DABS(del-1.0d0).lt.EPS) GOTO 1
11    CONTINUE
      call errquit('a too large, ITMAX too small in gcf',0, INPUT_ERR)
1     if((-x+a*LOG(x)-gln).gt.undovl) then
         gammcf=DEXP(-x+a*LOG(x)-gln)*h
      else
         gammcf=0d0
      endif
      return
      END 



      subroutine gser(gamser,a,x,gln)
      implicit none
#include "errquit.fh"
      real*8    a,x
      real*8   gamser,gln
    
*     **** parameters ****
      integer  ITMAX 
      parameter (ITMAX=1000)
      real*8  EPS
      parameter (EPS=3.0d-12)

*     **** local variables ****
      integer n
      real*8  ap,del,sum

*     **** external functions ***
      real*8   ln_gamma
      external ln_gamma


      gln=ln_gamma(a)

      IF(x.le.0.0d0) THEN
        IF(x.lt.0.0d0) call errquit('x < 0 in gser',0, UNKNOWN_ERR)
          gamser=0.0d0
          RETURN
      ENDIF

      ap=a
      sum=1.0d0/a
      del=sum
      DO n=1,ITMAX
        ap=ap+1.0d0
        del=del*x/ap
        sum=sum+del
        IF(DABS(del).lt.DABS(sum)*EPS) GOTO 1
      END DO

      CALL errquit('a too large, ITMAX too small in gser',0, INPUT_ERR)

1     gamser=sum*DEXP(-x+a*LOG(x)-gln)

      return
      END 
c $Id$
