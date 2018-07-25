c
c     $Id$
c

*     *********************************************
*     *                                           *
*     *              nwpw_spline                  *
*     *                                           *
*     *********************************************
      subroutine nwpw_spline(x,y,n,yp1,ypn,y2,utmp)
      implicit none
      real*8 x(*),y(*)
      integer n
      real*8 yp1,ypn
      real*8 y2(*),utmp(*)

*     **** local variables ****
      integer i,k
      real*8 sig,qn,un,p

      if (yp1.gt.0.99d30) then
        y2(1)   = 0.0d0
        utmp(1) = 0.0d0
      else
        y2(1)   = -0.5d0
        utmp(1) = 3.0d0 / (x(2)-x(1))
     >          * ( (y(2)-y(1)) / (x(2)-x(1)) - yp1 )
      end if
      do i=2,n-1
        sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.0d0
        y2(i) = (sig-1.0d0) / p
        utmp(i) = ( 
     >           6.0d0 * 
     >           ( 
     >             (y(i+1)-y(i  )) / (x(i+1)-x(i  ))
     >           - (y(i  )-y(i-1)) / (x(i  )-x(i-1))
     >           )
     >           / (x(i+1)-x(i-1)) - sig*utmp(i-1)
     >         )
     >         / p
      end do
      if (ypn.gt.0.99d30) then
        qn = 0.0d0
        un = 0.0d0
      else
        qn = 0.5d0
        un = 3.0d0 / (x(n)-x(n-1))
     >        * ( ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)) )
      end if
      y2(n) = (un-qn*utmp(n-1)) / (qn*y2(n-1)+1.0d0)
      do k=n-1,1,-1
         y2(k) = y2(k)*y2(k+1) + utmp(k)
      end do
      return
      end

*     *********************************************
*     *                                           *
*     *              nwpw_splint                  *
*     *                                           *
*     *********************************************
      real*8 function nwpw_splint(xa,ya,y2a,n,nx,x)
      implicit none
      real*8 xa(*),ya(*),y2a(*)
      integer n,nx
      real*8 x

*     **** local variables ****
      integer khi,klo
      real*8 a,b,h,y

      khi = nx+1
      klo = nx

      do while ((xa(klo).gt.x) .or. (xa(khi).lt.x))
         if (xa(klo).gt.x) then
           klo = klo - 1
           khi = khi - 1
         end if
         if (xa(khi).lt.x) then
           klo = klo + 1
           khi = khi + 1
         end if
      end do

      h = xa(khi)-xa(klo)
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) 
     >  + b*ya(khi) 
     >  + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h**2 / 6.0d0
      nwpw_splint = y
      return
      end


*     *********************************************
*     *                                           *
*     *              nwpw_splint_nx               *
*     *                                           *
*     *********************************************
      integer function nwpw_splint_nx(xa,nx,x)
      implicit none
      real*8 xa(*)
      integer nx
      real*8 x
*     **** local variables ****

      integer khi,klo
      real*8 a,b,h,y

      khi = nx+1
      klo = nx
      do while ((xa(klo).gt.x) .or. (xa(khi).lt.x))
         if (xa(klo).gt.x) then
           klo = klo - 1
           khi = khi - 1
         end if
         if (xa(khi).lt.x) then
           klo = klo + 1
           khi = khi + 1
         end if
      end do

      nwpw_splint_nx = klo
      return
      end 

*     *********************************************
*     *                                           *
*     *              nwpw_splint0                 *
*     *                                           *
*     *********************************************
      real*8 function nwpw_splint0(xa,ya,y2a,n,nx,x)
      implicit none
      real*8 xa(*),ya(*),y2a(*)
      integer n,nx
      real*8 x

*     **** local variables ****
      integer khi,klo
      real*8 a,b,h,y

      khi = nx+1
      klo = nx
      h = xa(khi)-xa(klo)
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)
     >  + b*ya(khi)
     >  + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h**2 / 6.0d0
      nwpw_splint0 = y
      return
      end



*     *********************************************
*     *                                           *
*     *              nwpw_dsplint                 *
*     *                                           *
*     *********************************************
      real*8 function nwpw_dsplint(xa,ya,y2a,n,nx,x)
      implicit none
      real*8 xa(*),ya(*),y2a(*)
      integer n,nx
      real*8 x

*     **** local variables ****
      integer khi,klo
      real*8 a,b,h,da,db,dy

      khi = nx+1
      klo = nx

      do while ((xa(klo).gt.x) .or. (xa(khi).lt.x))
         if (xa(klo).gt.x) then
           klo = klo - 1
           khi = khi - 1
         end if
         if (xa(khi).lt.x) then
           klo = klo + 1
           khi = khi + 1
         end if
      end do

      h = xa(khi)-xa(klo)
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      da = -1.0d0/h
      db =  1.0d0/h
c      y = a*ya(klo)
c     >  + b*ya(khi)
c     >  + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h**2 / 6.0d0

      dy = da*ya(klo)
     >   + db*ya(khi)
     >   + ( (da*3*a**2-da)*y2a(klo) 
     >    +  (db*3*b**2-db)*y2a(khi)) * h**2 / 6.0d0

      nwpw_dsplint = dy
      return
      end



*     *********************************************
*     *                                           *
*     *              nwpw_dsplint0                *
*     *                                           *
*     *********************************************
      real*8 function nwpw_dsplint0(xa,ya,y2a,n,nx,x)
      implicit none
      real*8 xa(*),ya(*),y2a(*)
      integer n,nx
      real*8 x

*     **** local variables ****
      integer khi,klo
      real*8 a,b,h,da,db,dy

      khi = nx+1
      klo = nx

      h = xa(khi)-xa(klo)
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      da = -1.0d0/h
      db =  1.0d0/h
c      y = a*ya(klo)
c     >  + b*ya(khi)
c     >  + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h**2 / 6.0d0

      dy = da*ya(klo)
     >   + db*ya(khi)
     >   + ( (da*3*a**2-da)*y2a(klo)
     >    +  (db*3*b**2-db)*y2a(khi)) * h**2 / 6.0d0

      nwpw_dsplint0 = dy
      return
      end

