      Subroutine igamma(Fj,N,L)
c $Id: igamma.F,v 1.8 2001-08-19 21:26:17 bjohnson Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
#ifdef BETTIS
#include "screened_op.fh"
#endif
      Parameter (EXPLIM=100.D0)
      Parameter (TMAX=30.D0,ACCY=1.D-14)
*
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
#ifdef BETTIS
      if (screened_operator_enabled()) then
         call errquit('Screened operator on in igamma',0)
      endif
#endif
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

*acc_debug:      end
*acc_debug:      Subroutine igamma_acc(Fj,N,L,accy,met)
*acc_debug:
*acc_debug:      Implicit real*8 (a-h,o-z)
*acc_debug:      Implicit integer (i-n)
*acc_debug:
*acc_debug:      Parameter (EXPLIM=100.D0)
*acc_debug:      Parameter (TMAX=30.D0)
*acc_debug:*
*acc_debug:      Parameter (PI=3.141592653589793D0)
*acc_debug:      Parameter (PI4=0.25D0*PI)
*acc_debug:
*acc_debug:      Dimension Fj(N,0:L)
*acc_debug:      Logical met
*acc_debug:
*acc_debug:c Compute the incomplete gamma function.
*acc_debug:c
*acc_debug:c               /1   2j
*acc_debug:c     Fj(T)  =  |   x  exp( - T x**2 ) dx
*acc_debug:c               /0
*acc_debug:c
*acc_debug:c******************************************************************************
*acc_debug:
*acc_debug:      met = .true.
*acc_debug:
*acc_debug:      do 100 m = 1,N
*acc_debug:
*acc_debug:      Tm = Fj(m,0)
*acc_debug:
*acc_debug:      if( Tm.GT.TMAX )then
*acc_debug:
*acc_debug:       T2inv = 0.5D0/Tm
*acc_debug:       expT = exp(-min(EXPLIM,Tm))
*acc_debug:
*acc_debug:c ... upward recursion.
*acc_debug:
*acc_debug:       Fj(m,0) = sqrt(PI4/Tm)
*acc_debug:       do 10 j = 1,L
*acc_debug:        Fj(m,j) = ((2*j-1)*Fj(m,j-1) - expT)*T2inv
*acc_debug:  10   continue
*acc_debug:
*acc_debug:      else
*acc_debug:
*acc_debug:       T2 = 2.D0*Tm
*acc_debug:       expT = exp(-min(EXPLIM,Tm))
*acc_debug:
*acc_debug:c ... series expansion.
*acc_debug:
*acc_debug:       TERM  = 1 / dble(2*L+1)
*acc_debug:
*acc_debug:       sum  = TERM
*acc_debug:       do 20 i = 1,1000
*acc_debug:        TERM = TERM*T2/dble(2*(L+i)+1)
*acc_debug:        sum = sum + TERM
*acc_debug:        if( TERM.lt.ACCY ) go to 25
*acc_debug:  20   continue
*acc_debug:
*acc_debug:       write(*,*) 'IGAMMA_ACC:  Unable to achieve required accuracy.'
*acc_debug:       write(*,*) '         TERM = ',TERM,' ACCY = ',ACCY
*acc_debug:       met = .false.
*acc_debug:       return
*acc_debug:
*acc_debug:  25   continue
*acc_debug:
*acc_debug:c ... downward recursion.
*acc_debug:
*acc_debug:       Fj(m,L) = expT*sum
*acc_debug:       do 30 j = L-1,0,-1
*acc_debug:        Fj(m,j) = (T2*Fj(m,j+1) + expT) / dble(2*j+1)
*acc_debug:  30   continue
*acc_debug:
*acc_debug:      end if
*acc_debug:
*acc_debug:  100 continue
*acc_debug:
      end
