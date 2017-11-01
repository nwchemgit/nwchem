c     $Id$                       
*     ********************************
*     *                              *
*     *      linesearch_init	     *
*     *                              *
*     ********************************
      subroutine linesearch_init()
      implicit none

*     **** linesearch_counter common block ****
      integer counter,maxiter
      common / linesearch_counter / counter,maxiter

!$OMP MASTER
      counter = 0
      maxiter = 8
!$OMP END MASTER

      return
      end

*     ********************************
*     *                              *
*     *   linesearch_maxiter_set     *
*     *                              *
*     ********************************
      subroutine linesearch_maxiter_set(maxiter0)
      implicit none
      integer maxiter0

*     **** linesearch_counter common block ****
      integer counter,maxiter
      common / linesearch_counter / counter,maxiter

!$OMP MASTER
      maxiter = maxiter0
!$OMP END MASTER
      return
      end


*     ********************************
*     *                              *
*     *      linesearch_maxiter      *
*     *                              *
*     ********************************
      integer function linesearch_maxiter()
      implicit none

*     **** linesearch_counter common block ****
      integer counter,maxiter
      common / linesearch_counter / counter,maxiter

      linesearch_maxiter = maxiter
      return
      end


*     ********************************
*     *                              *
*     *      linesearch_count        *
*     *                              *
*     ********************************
      integer function linesearch_count()
      implicit none

*     **** linesearch_counter common block ****
      integer counter,maxiter
      common / linesearch_counter / counter,maxiter
 
      linesearch_count = counter
      return
      end

*     ********************************
*     *                              *
*     *         linesearch	         *
*     *                              *
*     ********************************

      real*8 function linesearch(t0,f0,df0,deltat,FUNC1,DFUNC1,
     >                           tolerance,tmin,deltaE,
     >                           stoptype)
      implicit none
      real*8 t0,f0,df0
      real*8 deltat
      real*8   FUNC1
      external FUNC1
      real*8   DFUNC1
      external DFUNC1
      real*8   tolerance
      real*8   tmin
      real*8   deltaE
      integer  stoptype


*     **** local variables ****
      !integer maxiter,iteration
      !parameter (maxiter=8)
      integer iteration
      logical secant,notfinished
      integer indx(3)
      real*8  t(3)
      real*8  f(3)
      real*8  df(3)
      real*8  t_last, f_last
      real*8  t_first,f_first,df_first
      real*8 up,down,fmin,dfmin,deltaf

*     **** linesearch_counter common block ****
      integer counter,maxiter
      common / linesearch_counter / counter,maxiter

    
!$OMP MASTER 
      counter = counter + 1
!$OMP END MASTER 

      secant = .true.

      t(1) = t0
c     f(1)  =  FUNC1(t(1))
c     df(1) = DFUNC1(t(1))
      f(1)  = f0
      df(1) = df0
      t_last = t(1)
      f_last = f(1)
      t_first = t(1)
      f_first = f(1)
      df_first = df(1)

      f(2)  =  FUNC1(t(1)+deltat)
      df(2) = DFUNC1(t(1)+deltat)
      
      iteration = 1
*     **** make sure that f2 < f1 ******
      do while ((f(2).gt.f(1)).and.(iteration.le.maxiter))
        deltat = 0.5d0*deltat
        f(2)  =  FUNC1(t(1)+deltat)
        df(2) = DFUNC1(t(1)+deltat)
        iteration = iteration + 1
      end do
      t(2) = t(1) + deltat
      t_last = t(2)
      f_last = f(2)

*     **** use secant method to generate f(3) *****
      deltat = -df(1)*(t(2)-t(1))/(df(2)-df(1))
      t(3)   = t(1) + deltat
      f(3)   =  FUNC1(t(3))
      df(3)  = DFUNC1(t(3))
      iteration = iteration + 1
      t_last = t(3)
      f_last = f(3)

*     **** sort the function values ****
      call Order_Values(3,f,indx)

      deltaf = f(indx(2)) - f(indx(1))

      if (stoptype.eq.1) then
        notfinished = (dabs(deltaf).gt.tolerance)
     >                .and.(iteration.le.maxiter)
      else
        notfinished = (dabs(df(indx(1))/df_first).gt.tolerance)
     >                .and.(iteration.le.maxiter)
      end if

      do while (notfinished) 
      

*       **** use secant interpolation to generate tmin ***
        if (secant) then

          deltat = -df(indx(1))
     >           *(t(indx(2))-t(indx(1)))
     >           /(df(indx(2))-df(indx(1)))
          tmin = t(indx(1)) + deltat
          fmin  =  FUNC1(tmin)
          dfmin = DFUNC1(tmin)
          iteration = iteration + 1
          t_last = tmin
          f_last = fmin

*         **** finish using secant method ****
          if (fmin.ge.f(indx(1))) then
            secant = .false.
            if (fmin.lt.f(indx(3))) then
              t(indx(3))  = tmin
              f(indx(3))  = fmin
              df(indx(3)) = dfmin
              call Order_Values(3,f,indx)
            end if
          end if


        end if 

*       **** use quadradic interpolation to generate tmin ***
        if (.not.secant) then
          up  = (t(2)*t(2) - t(3)*t(3))*f(1)
     >        + (t(3)*t(3) - t(1)*t(1))*f(2)
     >        + (t(1)*t(1) - t(2)*t(2))*f(3)
          down = (t(2) - t(3))*f(1)
     >         + (t(3) - t(1))*f(2)
     >         + (t(1) - t(2))*f(3)

*         **** check added by E.Apra ****
          if(abs(down).gt.tolerance**2) then
             tmin = 0.5d0*up/down
             fmin  =  FUNC1(tmin)
             dfmin = DFUNC1(tmin)
             iteration = iteration + 1
             t_last = tmin
             f_last = fmin

*         **** parabola fit failed - exit loop ****
          else
             tmin=t(indx(3))
             fmin=f(indx(3))+tolerance
             iteration = maxiter+1
          endif

        end if


*       **** tolerance check and replacement ****
        if (fmin.lt.f(indx(3))) then
           t(indx(3))  = tmin
           f(indx(3))  = fmin
           df(indx(3)) = dfmin
           call Order_Values(3,f,indx)
           deltaf = f(indx(2)) - f(indx(1))
        else
           deltaf=0.0d0
        end if

        if (stoptype.eq.1) then
          notfinished = (dabs(deltaf).gt.tolerance)
     >                .and.(iteration.le.maxiter)
        else
          notfinished = (dabs(df(indx(1))/df_first).gt.tolerance)
     >                .and.(iteration.le.maxiter)
        end if
           
      end do

*     **** make sure that tmin is last functions call *****
      tmin = t(indx(1))
      fmin = f(indx(1))
      if (tmin.ne.t_last) fmin = FUNC1(tmin)

      deltaE = (fmin-f_first)
  
      linesearch = fmin
      return
      end
