      Subroutine igamma(Fj,N,L)
c $Id: igamma.F,v 1.4 1996-07-11 22:21:29 d3e129 Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Parameter (EXPLIM=100.D0)
      Parameter (TMAX=30.D0,ACCY=1.D-14)

      Parameter (PI=3.141592653589793D0)
      Parameter (PI4=0.25D0*PI)

      Dimension Fj(N,0:L)

c Compute the incomplete gamma function.
c
c               /1   2j
c     Fj(T)  =  |   x  exp( - T x**2 ) dx
c               /0
c
c******************************************************************************

      do 100 m = 1,N

      Tm = Fj(m,0)

      if( Tm.GT.TMAX )then

       T2inv = 0.5D0/Tm
       expT = exp(-min(EXPLIM,Tm))

c ... upward recursion.

       Fj(m,0) = sqrt(PI4/Tm)
       do 10 j = 1,L
        Fj(m,j) = ((2*j-1)*Fj(m,j-1) - expT)*T2inv
  10   continue

      else

       T2 = 2.D0*Tm
       expT = exp(-min(EXPLIM,Tm))

c ... series expansion.

       TERM  = 1 / dble(2*L+1)

       sum  = TERM
       do 20 i = 1,1000
        TERM = TERM*T2/dble(2*(L+i)+1)
        sum = sum + TERM
        if( TERM.lt.ACCY ) go to 25
  20   continue

       write(*,*) 'IGAMMA:  Unable to achieve required accuracy.'
       write(*,*) '         TERM = ',TERM,' ACCY = ',ACCY
       stop

  25   continue

c ... downward recursion.

       Fj(m,L) = expT*sum
       do 30 j = L-1,0,-1
        Fj(m,j) = (T2*Fj(m,j+1) + expT) / dble(2*j+1)
  30   continue

      end if

  100 continue

      end
