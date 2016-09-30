!234567890
      subroutine spline(f,x,n,t,d2f)
!     constructs a spline assuming periodic boundary conditions
!     input: coordinate x, function value f, number of points n, period t
!     output: second derivatives d2f
      implicit none
      integer n,nm
      double precision x(n),f(n),d2f(n),t
      integer ii,jj,iperm(n),nperm
      double precision cmat(n,n)

      do ii=1,n
      do jj=1,n
        cmat(ii,jj)=0.d0
      enddo
      enddo
      d2f(1)=(f(2)-f(1))/(x(2)-x(1)) &
&           -(f(1)-f(n))/(x(1)-x(n)+t)
      cmat(1,n)=(x(1)-x(n)+t)/6.d0
      cmat(1,1)=(x(2)-x(n)+t)/3.d0
      cmat(1,2)=(x(2)-x(1))/6.d0
      do jj=2,n-1
        d2f(jj)=(f(jj+1)-f(jj))/(x(jj+1)-x(jj)) &
&              -(f(jj)-f(jj-1))/(x(jj)-x(jj-1))
        cmat(jj,jj-1)=(x(jj)-x(jj-1))/6.d0
        cmat(jj,jj)=(x(jj+1)-x(jj-1))/3.d0
        cmat(jj,jj+1)=(x(jj+1)-x(jj))/6.d0
      enddo
      d2f(n)=(f(1)-f(n))/(x(1)-x(n)+t) &
&           -(f(n)-f(n-1))/(x(n)-x(n-1))
      cmat(n,n-1)=(x(n)-x(n-1))/6.d0
      cmat(n,n)=(x(1)-x(n-1)+t)/3.d0
      cmat(n,1)=(x(1)-x(n)+t)/6.d0

      call loupd(cmat,n,n,iperm,nperm)

      call baksublu(cmat,d2f,n,n,iperm)

      return
      end

!*************************************************************************

      subroutine dspline(f,x,d2f,n,df)
!     constructs table of first derivatives of spline function
!     through values f on grid n with second derivatives d2f
      implicit none
      integer n
      double precision f(n),x(n),d2f(n),df(n)
      integer ii,jj

      do ii=1,n-1
        df(ii)=(f(ii+1)-f(ii))/(x(ii+1)-x(ii)) &
&             -(x(ii+1)-x(ii))*d2f(ii)/3.d0 &
&             -(x(ii+1)-x(ii))*d2f(ii+1)/6.d0
      enddo
      df(n)=(f(n)-f(n-1))/(x(n)-x(n-1)) &
&          +(x(n)-x(n-1))*d2f(n-1)/6.d0 &
&          +(x(n)-x(n-1))*d2f(n)/3.d0

      return
      end


!*************************************************************************

      subroutine dsplv(f,x,t,d2f,n,xv,dfv)
!     first derivative of spline function
!     through values f on grid n with second derivatives d2f
      implicit none
      integer n
      double precision f(n),x(n),t,d2f(n),xv,dfv
      integer ii,jj
      integer ihi,ilo
      double precision aa,bb,dx

      if (xv.lt.x(1)) then
        ihi=1
        ilo=n
        dx=x(ihi)-x(ilo)+t
      else
        ihi=n
        ilo=1
 100    if (ihi-ilo.gt.1) then
          ii=(ihi+ilo)/2
          if (x(ii).gt.xv) then
            ihi=ii
          else
            ilo=ii
          endif
          goto 100
        endif
        dx=x(ihi)-x(ilo)
      endif
      aa=(x(ihi)-xv)/dx
      bb=1-aa
      dfv=(f(ihi)-f(ilo))/dx &
&        -dx*d2f(ilo)*(3*aa**2-1.d0)/6.d0 &
&        +dx*d2f(ihi)*(3*bb**2-1.d0)/6.d0

      return
      end


!*************************************************************************

      subroutine evalspline(f,x,t,d2f,n,xv,fv)
!     evaluate spline through periodic function f on points x 
!     with second derivatives d2f at value xv and period t, returning fv
      implicit none
      integer n
      double precision f(n),x(n),d2f(n),t
      double precision xv,fv
      integer ihi,ilo,ii
      double precision aa,bb,dx,f1,f2

      if (xv.lt.x(1)) then
        ihi=1
        ilo=n
        dx=x(ihi)-x(ilo)+t
      else
        ihi=n
        ilo=1
 100    if (ihi-ilo.gt.1) then
          ii=(ihi+ilo)/2
          if (x(ii).gt.xv) then
            ihi=ii
          else
            ilo=ii
          endif
          goto 100
        endif
        dx=x(ihi)-x(ilo)
      endif
      aa=(x(ihi)-xv)/dx
      bb=1-aa
      fv=aa*f(ilo)+bb*f(ihi) &
&        +((aa**3-aa)*d2f(ilo)+(bb**3-bb)*d2f(ihi))*dx**2/6.d0

      return
      end
