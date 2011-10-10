!      program main
!      integer i,ndata
!      double precision x,erf
!      ndata=20
!      x=0.0d0
!      do i=1,ndata
!       call merf(x,erf)
!       write(*,1) i,x,erf
!1      format('(x,erf)(',i3,')=(',f15.8,',',f15.8,')') 
!       x=x+0.10d0
!      enddo
!      end
       subroutine merf (y, erf)

!  ======================================================================
!  purpose:  evaluation of the error function.
!            two related entries are provided: the complementary
!            error function and the normal distribution function.
!            (entries merfc and mdnor)
!
!  input  :  y      - the argument of the function
!  output :  erf    - the function value
!                     merf : ranging from -1 to +1 for increasing y
!                     merfc  : ranging from +2 to 0
!                     mdnor  : max. for y=0, falling off to zero for
!                              increasing abs(y)
!
!  *=====================================================================
!  ======================================================================
      double precision c,eight,half,one,r0p47,
     &                 r5p5,small,sqr1d2,tenth,thirdm,
     &                 two,zero,a,p,q,x,y,erf,sn,xx
      integer kret,md,isw   
!
      parameter (c = 1.1283791670955d0, eight = 8.0d0, 
     &           half = 0.5d0, one = 1.0d0, r0p47 = 0.47d0,     
     &           r5p5 = 5.5d0, small = 1.0d-15, 
     &           sqr1d2 = 0.70710678118655d0, 
     &           tenth = 0.1d0, thirdm =  
     &          -one/3.0d0, ttt = -one/1320.0d0, two = 2.0d0, 
     &           zero = 0.0d0)
!
      dimension a(10), p(8), q(9)
!
      data a /one, thirdm, tenth, -0.02380923809238d0, 
     &        0.4629629629630d-2, 
     &       ttt, 1.0683760683761d-4,   
     &       -1.3227513227513d-5, 1.4589169000934d-6, 
     &       -1.45038552223150d-7/
!
      data p /883.4789426085d0, 1549.6793124037d0, 
     &        1347.1941340976d0, 723.04000277753d0,              
     &        255.50049469496d0, 59.240010112914d0, 
     &        8.3765310814197d0, 0.56418955944261d0/
!
      data q /883.4789426085d0, 2546.5785458098d0, 
     &        3337.2213699893d0, 2606.7120152651d0,              
     &        1333.5699756800d0, 460.28512369160d0, 
     &        105.50025439769d0, 14.847012237523d0, one/

             kret = 2
          !
             x  = y
             md = 0
          !
             if (x<zero) then
                isw = 1
                x   = - x
             else
                isw = 2
             end if
          !
             if (kret/=2) goto 120
             if (x<=r5p5) goto 40
          !
             erf = one
             if (isw==1) erf = - one
             goto 150
          
          !  -----------------------------------------------------------
          !  abs(x) less than .47 compute erf by taylor series expansion
          !  -----------------------------------------------------------
          
   40        if (x<=r0p47) goto 90
             kret = 1
          
          !  -----------------------------------------------------
          !  abs(x) between .47 and 10. compute complemented error
          !  function by a rational function of x
          !  -----------------------------------------------------
          
   50        sn = p(8)
            do i = 7, 1, -1
               sn = sn*x + p(i)
            end do
          !
             sd = q(9)
             do i = 8, 1, -1
                sd = sd*x + q(i)
             end do
          !
             erf = (sn/sd)*exp(-x*x)
             if (kret/=1) goto 80
          
          !  -----------------------------------
          !  compute complemented error function
          !  -----------------------------------
          
             erf = one - erf
             if (isw/=2) erf = - erf
             goto 150
          !
   80         if (isw/=2) erf = two - erf
             if (md/=0) erf = half*erf
             goto 150
          !
   90        xx = x*x
          !
             erf = a(10)
             do i = 9, 1, -1
                erf = erf*xx + a(i)
             end do
          !
             erf = c*erf*x
          !
             if (kret==2) then
                if (isw/=2) erf = - erf
                goto 150
             else
                erf = one - erf
                goto 80
             end if
          !
  120        if (x>eight) then
                erf = zero
                goto 80
          !
             else if (x>r0p47) then
                kret = 2
                goto 50
          !
             else if (x>small) then
                goto 90
          !
             else
                erf = one
                goto 80
             end if
          !
  150     return
          end
c $Id$
