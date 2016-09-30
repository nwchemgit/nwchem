subroutine trispline(f,x,y,z,n1,n2,n3,t1,t2,t3,fx,fy,fxy)
!subroutine trispline(f,x,y,z,n1,n2,n3,t1,t2,t3,fx,fy,fz,fxy,fxz,fyz,fxyz)
!     constructs interpolation coefficients for a periodic 
!     function f evaluated on a 
!     three dimensional grid (x,y,z) with size n1 x n2 x n3.
!     t1,t2,t3 are the periods.
implicit none
integer n1,n2,n3
double precision f(n1,n2,n3),x(n1),y(n2),z(n3),t1,t2,t3
double precision fx(n1,n2,n3),fy(n1,n2,n3),fxy(n1,n2,n3)   !, &
!&                fz(n1,n2,n3),fxz(n1,n2,n3),fyz(n1,n2,n3), &
!&                fxyz(n1,n2,n3)
double precision, allocatable :: coefs(:)
integer ii,jj,kk

allocate (coefs(n1))
do jj=1,n2
do kk=1,n3
  call spline(f(:,jj,kk),x,n1,t1,coefs)
  call dspline(f(:,jj,kk),x,coefs,n1,fx(:,jj,kk))
enddo
enddo

deallocate (coefs)
allocate (coefs(n2))
do ii=1,n1
do kk=1,n3
  call spline(f(ii,:,kk),y,n2,t2,coefs)
  call dspline(f(ii,:,kk),y,coefs,n2,fy(ii,:,kk))
  call spline(fx(ii,:,kk),y,n2,t2,coefs)
  call dspline(fx(ii,:,kk),y,coefs,n2,fxy(ii,:,kk))
enddo
enddo
      
!deallocate (coefs)
!allocate (coefs(n3))
!do ii=1,n1
!do jj=1,n2
!  call spline(f(ii,jj,:),z,n3,t3,coefs)
!  call dspline(f(ii,jj,:),z,coefs,n3,fz(ii,jj,:))
!  call spline(fx(ii,jj,:),z,n3,t3,coefs)
!  call dspline(fx(ii,jj,:),z,coefs,n3,fxz(ii,jj,:))
!  call spline(fy(ii,jj,:),z,n3,t3,coefs)
!  call dspline(fy(ii,jj,:),z,coefs,n3,fyz(ii,jj,:))
!  call spline(fxy(ii,jj,:),z,n3,t3,coefs)
!  call dspline(fxy(ii,jj,:),z,coefs,n3,fxyz(ii,jj,:))
!enddo
!enddo
      
return
end subroutine trispline

!***************************************************************************

subroutine xyzspline(f,x,y,z,n1,n2,n3,t1,t2,t3,fx,fy,fz,fxy,fxz,fyz,fxyz)
!     constructs interpolation coefficients for a periodic 
!     function f evaluated on a 
!     three dimensional grid (x,y,z) with size n1 x n2 x n3.
!     t1,t2,t3 are the periods.
implicit none
integer n1,n2,n3
double precision f(n1,n2,n3),x(n1),y(n2),z(n3),t1,t2,t3
double precision fx(n1,n2,n3),fy(n1,n2,n3),fxy(n1,n2,n3), &
&                fz(n1,n2,n3),fxz(n1,n2,n3),fyz(n1,n2,n3), &
&                fxyz(n1,n2,n3)
double precision, allocatable :: coefs(:)
integer ii,jj,kk

allocate (coefs(n1))
do jj=1,n2
do kk=1,n3
  call spline(f(:,jj,kk),x,n1,t1,coefs)
  call dspline(f(:,jj,kk),x,coefs,n1,fx(:,jj,kk))
enddo
enddo

deallocate (coefs)
allocate (coefs(n2))
do ii=1,n1
do kk=1,n3
  call spline(f(ii,:,kk),y,n2,t2,coefs)
  call dspline(f(ii,:,kk),y,coefs,n2,fy(ii,:,kk))
  call spline(fx(ii,:,kk),y,n2,t2,coefs)
  call dspline(fx(ii,:,kk),y,coefs,n2,fxy(ii,:,kk))
enddo
enddo
      
deallocate (coefs)
allocate (coefs(n3))
do ii=1,n1
do jj=1,n2
  call spline(f(ii,jj,:),z,n3,t3,coefs)
  call dspline(f(ii,jj,:),z,coefs,n3,fz(ii,jj,:))
  call spline(fx(ii,jj,:),z,n3,t3,coefs)
  call dspline(fx(ii,jj,:),z,coefs,n3,fxz(ii,jj,:))
  call spline(fy(ii,jj,:),z,n3,t3,coefs)
  call dspline(fy(ii,jj,:),z,coefs,n3,fyz(ii,jj,:))
  call spline(fxy(ii,jj,:),z,n3,t3,coefs)
  call dspline(fxy(ii,jj,:),z,coefs,n3,fxyz(ii,jj,:))
enddo
enddo
      
return
end subroutine xyzspline

!***************************************************************************

subroutine spl23(f,fx,fy,fxy,x,y,z,n1,n2,n3,t1,t2,t3,xv,yv,nz,zv,fv)
implicit none
integer :: n1,n2,n3,nz
double precision :: f(n1,n2,n3),x(n1),y(n2),z(n3),t1,t2,t3,xv,yv,zv(nz)
double precision :: fv(nz)
double precision :: fx(n1,n2,n3),fy(n1,n2,n3),fxy(n1,n2,n3),zvals(n3)
double precision, allocatable :: coefs(:)
double precision :: veca(4),vecb(4),dmat(4,4)
integer :: ii,jj,kk
integer :: ihi,ilo,jhi,jlo
double precision :: dx,dxg,dy,dyg

if (xv.lt.x(1)) then
  ihi=1
  ilo=n1
  dxg=x(ihi)-x(ilo)+t1
  dx=(xv-x(ilo)+t1)/dxg
else
  ihi=n1
  ilo=1
  do
    if (ihi-ilo.le.1) exit
    ii=(ihi+ilo)/2
    if (x(ii).gt.xv) then
      ihi=ii
    else
      ilo=ii
    endif
  enddo
  dxg=x(ihi)-x(ilo)
  dx=(xv-x(ilo))/dxg
endif

if (yv.lt.y(1)) then
  jhi=1
  jlo=n2
  dyg=y(jhi)-y(jlo)+t2
  dy=(yv-y(jlo)+t2)/dyg
else
  jhi=n2
  jlo=1
  do
    if (jhi-jlo.le.1) exit
    jj=(jhi+jlo)/2
    if (y(jj).gt.yv) then
      jhi=jj
    else
      jlo=jj
    endif
  enddo
  dyg=y(jhi)-y(jlo)
  dy=(yv-y(jlo))/dyg
endif

!if (zv.lt.z(1)) then
!  khi=1
!  klo=n3
!  dzg=z(khi)-z(klo)+t3
!  dz=(zv-z(klo)+t3)/dzg
!else
!  khi=n3
!  klo=1
!  do
!    if (khi-klo.le.1) exit
!    kk=(khi+klo)/2
!    if (z(kk).gt.zv) then
!      khi=kk
!    else
!      klo=kk
!    endif
!  enddo
!  dzg=z(khi)-z(klo)
!  dz=(zv-z(klo))/dzg
!endif

do kk=1,n3
!  allocate (coefs(n1))
!  do jj=1,n2
!    call spline(f(:,jj,kk),x,n1,t1,coefs)
!    call dspline(f(:,jj,kk),x,coefs,n1,fx(:,jj,kk))
!  enddo
!  deallocate (coefs)
!  
!  allocate (coefs(n2))
!  do ii=1,n1
!    call spline(f(ii,:,kk),y,n2,t2,coefs)
!    call dspline(f(ii,:,kk),y,coefs,n2,fy(ii,:,kk))
!    call spline(fx(ii,:,kk),y,n2,t2,coefs)
!    call dspline(fx(ii,:,kk),y,coefs,n2,fxy(ii,:,kk))
!  enddo
!  deallocate (coefs)
!
  veca(1)=(1.d0-dx)**2*(1.d0+2.d0*dx)
  veca(2)=dx**2*(3.d0-2.d0*dx)
  veca(3)=(1.d0-dx)**2*dx*dxg
  veca(4)=(dx-1.d0)*dx**2*dxg

  vecb(1)=(1.d0-dy)**2*(1.d0+2.d0*dy)
  vecb(2)=dy**2*(3.d0-2.d0*dy)
  vecb(3)=(1.d0-dy)**2*dy*dyg
  vecb(4)=(dy-1.d0)*dy**2*dyg

  dmat(1,1)=f(ilo,jlo,kk)
  dmat(1,2)=f(ilo,jhi,kk)
  dmat(1,3)=fy(ilo,jlo,kk)
  dmat(1,4)=fy(ilo,jhi,kk)
  dmat(2,1)=f(ihi,jlo,kk)
  dmat(2,2)=f(ihi,jhi,kk)
  dmat(2,3)=fy(ihi,jlo,kk)
  dmat(2,4)=fy(ihi,jhi,kk)
  dmat(3,1)=fx(ilo,jlo,kk)
  dmat(3,2)=fx(ilo,jhi,kk)
  dmat(3,3)=fxy(ilo,jlo,kk)
  dmat(3,4)=fxy(ilo,jhi,kk)
  dmat(4,1)=fx(ihi,jlo,kk)
  dmat(4,2)=fx(ihi,jhi,kk)
  dmat(4,3)=fxy(ihi,jlo,kk)
  dmat(4,4)=fxy(ihi,jhi,kk)

  zvals(kk)=0.d0
  do ii=1,4
  do jj=1,4
    zvals(kk)=zvals(kk)+dmat(jj,ii)*vecb(ii)*veca(jj)
  enddo
  enddo
enddo

allocate (coefs(n3))
call spline(zvals,z,n3,t3,coefs)
do kk=1,nz
  call evalspline(zvals,z,t3,coefs,n3,zv(kk),fv(kk))
enddo
deallocate (coefs)

return
end subroutine spl23
     
