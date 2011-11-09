      Subroutine igamma(Fj,N,L)
c $Id$

      Implicit none
      integer n
      integer l
      Double precision Fj(N,0:L)

c Compute the incomplete gamma function.
c
c               /1   2j
c     Fj(T)  =  |   x  exp( - T x**2 ) dx
c               /0
c
c******************************************************************************
      double precision x,faci
      integer lmax
      parameter (lmax=40)
      double precision ftex(0:lmax)
      integer m
      parameter(faci=1d0/1.128379167096d0)
c
      if(l.eq.0) then
        do  m = 1,N
           call fm(fj(m,0),0,ftex)
           fj(m,0)=ftex(0)*faci
        enddo
      elseif(l.eq.1) then
        do  m = 1,N
           call fm(fj(m,0),1,ftex)
           fj(m,0)=ftex(0)*faci
           fj(m+N,0)=ftex(1)*faci
        enddo
      elseif(l.eq.2) then
        do  m = 1,N
           call fm(fj(m,0),2,ftex)
           fj(m,0)=ftex(0)*faci
           fj(m+N,0)=ftex(1)*faci
           fj(m+N+N,0)=ftex(2)*faci
        enddo
      elseif(l.eq.3) then
        do  m = 1,N
           call fm(fj(m,0),3,ftex)
           fj(m,0)=ftex(0)*faci
           fj(m+N,0)=ftex(1)*faci
           fj(m+N+N,0)=ftex(2)*faci
           fj(m+N+N+N,0)=ftex(3)*faci
        enddo
      else
        do  m = 1,N
           call fm(fj(m,0),l,ftex)
           call dcopy(l+1,ftex,1,fj(m,0),N)
        enddo
        call dscal(N*(L+1),faci,fj(1,0),1)
      endif
      return
      end
      Subroutine igammao(Fj,N,L)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
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
      subroutine igamma_init
c
c     call need when txs is not active
c
      call fprep
      return
      end
