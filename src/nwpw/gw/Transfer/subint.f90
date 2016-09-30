subroutine fvq2(iqq,shiftk,shiftkq,bmet,ngkpt,ipw1,ipw2,ipwlf,npwlf,vq2)
implicit none
integer :: iqq(3),ngkpt(3),npwlf,ipw1,ipw2
integer :: ipwlf(3,npwlf)
double precision :: shiftk(3),shiftkq(3),bmet(3,3)
double precision :: vq2(8)
double precision :: qq(3),qp(3),qs(3),qq2,qp2
integer :: iv,ii,jj,ikk(3),iqqp(3),iks(3),ikslf1(3),ikslf2(3)

do iv=1,8
  ikk(1)=mod(iv-1,2)
  ikk(2)=mod((iv-1)/2,2)
  ikk(3)=mod((iv-1)/4,2)
  iqqp=iqq+ikk
  iks=iqqp-ngkpt/2
  ikslf1=iks+ipwlf(:,ipw1)
  ikslf2=iks+ipwlf(:,ipw2)
  if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
    call fqq(shiftk,shiftkq,ikslf1,ngkpt,qq,qs)
    call fqq(shiftk,shiftkq,ikslf2,ngkpt,qp,qs)
  else
    call fqq(shiftk,shiftk,ikslf1,ngkpt,qq,qs)
    call fqq(shiftk,shiftk,ikslf2,ngkpt,qp,qs)
  endif
  qq2=0.d0
  qp2=0.d0
  do ii=1,3
  do jj=1,3
    qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
    qp2=qp2+qp(ii)*bmet(ii,jj)*qp(jj)
  enddo
  enddo
  vq2(iv)=sqrt(qq2*qp2)
enddo

return
end subroutine fvq2

!***************************************************************************

subroutine floss(dielf,nwpt,nqpt,iwpt,ngkpt,iqndx,iqq,lossv)
! find loss function at cube vertices from tabulated dielectric function
implicit none
integer :: iqq(3),iqqp(3),ngkpt(3),nqpt,nwpt,iwpt
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: dielf(nwpt,nqpt)
double precision :: lossv(8)
integer :: ii

do ii=1,3
  iqqp(ii)=iqq(ii)+1
  if (iqqp(ii).gt.ngkpt(ii)) iqqp(ii)=iqqp(ii)-ngkpt(ii)
enddo

lossv(1)=dimag(1/dielf(iwpt,iqndx( iqq(1), iqq(2), iqq(3))))
lossv(2)=dimag(1/dielf(iwpt,iqndx(iqqp(1), iqq(2), iqq(3))))
lossv(3)=dimag(1/dielf(iwpt,iqndx( iqq(1),iqqp(2), iqq(3))))
lossv(4)=dimag(1/dielf(iwpt,iqndx(iqqp(1),iqqp(2), iqq(3))))
lossv(5)=dimag(1/dielf(iwpt,iqndx( iqq(1), iqq(2),iqqp(3))))
lossv(6)=dimag(1/dielf(iwpt,iqndx(iqqp(1), iqq(2),iqqp(3))))
lossv(7)=dimag(1/dielf(iwpt,iqndx( iqq(1),iqqp(2),iqqp(3))))
lossv(8)=dimag(1/dielf(iwpt,iqndx(iqqp(1),iqqp(2),iqqp(3))))

return
end subroutine floss

!***************************************************************************

subroutine subgridq2(nsub,sidelength,origin,bmet,q2grid,qgrid)
! evaluate q^2 at the gridpoints of a subdivided cube
implicit none
integer :: nsub
double precision :: sidelength(3),origin(3),q2grid(nsub+1,nsub+1,nsub+1)
double precision :: qgrid(3,nsub+1,nsub+1,nsub+1)
double precision :: bmet(3,3)
double precision :: qq(3)
integer :: ii,jj,ix,iy,iz

do ix=1,nsub+1
do iy=1,nsub+1
do iz=1,nsub+1
  qq=(/dble(ix-1),dble(iy-1),dble(iz-1)/)*sidelength/dble(nsub)+origin
  q2grid(ix,iy,iz)=0.d0
  do ii=1,3
  do jj=1,3
    q2grid(ix,iy,iz)=q2grid(ix,iy,iz)+qq(ii)*bmet(ii,jj)*qq(jj)
  enddo
  enddo
  q2grid(ix,iy,iz)=max(q2grid(ix,iy,iz),1.d-16)
  qgrid(:,ix,iy,iz)=qq
enddo
enddo
enddo

return
end subroutine subgridq2

!**************************************************************************

subroutine subcube(valvtx,nsub,valgrid)
! interpolate function known at vertices of cube onto finer grid.
! Divide cube into tetrahedrons and approximate function by its
! value at tetrahedron "origin" plus gradient.
implicit none
integer :: nsub
double precision :: valvtx(8),valgrid(nsub+1,nsub+1,nsub+1)
double precision :: contragrad(3,3,6),relvtx(3,3,6),grad(3,6),coord(3),point(3)
double precision :: rr(3,8),qq(3),qq2,volume
integer :: ivndx(4,6)
integer :: ii,jj,ix,iy,iz,itet,inside
data rr(1:3,1) /0.0d0,0.0d0,0.0d0/
data rr(1:3,2) /1.0d0,0.0d0,0.0d0/
data rr(1:3,3) /0.0d0,1.0d0,0.0d0/
data rr(1:3,4) /1.0d0,1.0d0,0.0d0/
data rr(1:3,5) /0.0d0,0.0d0,1.0d0/
data rr(1:3,6) /1.0d0,0.0d0,1.0d0/
data rr(1:3,7) /0.0d0,1.0d0,1.0d0/
data rr(1:3,8) /1.0d0,1.0d0,1.0d0/
data ivndx(1:4,1) /1,2,3,5/
data ivndx(1:4,2) /3,5,6,7/
data ivndx(1:4,3) /2,3,5,6/
data ivndx(1:4,4) /2,3,4,6/
data ivndx(1:4,5) /3,4,6,7/
data ivndx(1:4,6) /4,6,7,8/

do itet=1,6
! vertex coordinates relative to tetrahedron "origin"
  relvtx(:,1,itet)=rr(:,ivndx(2,itet))-rr(:,ivndx(1,itet))
  relvtx(:,2,itet)=rr(:,ivndx(3,itet))-rr(:,ivndx(1,itet))
  relvtx(:,3,itet)=rr(:,ivndx(4,itet))-rr(:,ivndx(1,itet))
! find contragradient vector
  call cross(relvtx(:,2,itet),relvtx(:,3,itet),contragrad(:,1,itet))
  call cross(relvtx(:,3,itet),relvtx(:,1,itet),contragrad(:,2,itet))
  call cross(relvtx(:,1,itet),relvtx(:,2,itet),contragrad(:,3,itet))
  volume=dot_product(contragrad(:,1,itet),relvtx(:,1,itet))
  contragrad(:,:,itet)=contragrad(:,:,itet)/volume
! find gradient
  do ii=1,3
    grad(ii,itet)=0.d0
    do jj=1,3
      grad(ii,itet)=grad(ii,itet)+contragrad(ii,jj,itet)*(valvtx(ivndx(jj+1,itet))-valvtx(ivndx(1,itet)))
    enddo
  enddo
!  write(42,'(3es10.3,5x,es20.13)') grad(:,itet),sqrt(dot_product(grad(:,itet),grad(:,itet)))
enddo

do ix=1,nsub+1
do iy=1,nsub+1
do iz=1,nsub+1
  coord=(/dble(ix-1),dble(iy-1),dble(iz-1)/)/dble(nsub)
  do itet=1,6
    point=coord-rr(:,ivndx(1,itet))
    call tbound(relvtx(:,:,itet),contragrad(:,:,itet),point,inside)
    if (inside.eq.1) then
      valgrid(ix,iy,iz)=valvtx(ivndx(1,itet))
      do ii=1,3
        valgrid(ix,iy,iz)=valgrid(ix,iy,iz)+(coord(ii)-rr(ii,ivndx(1,itet)))*grad(ii,itet)
      enddo
      exit
    endif
  enddo
enddo
enddo
enddo

return
end subroutine subcube

!**************************************************************************

subroutine subcubec(valvtx,nsub,valgrid)
! interpolate complex function known at vertices of cube onto finer grid.
! Divide cube into tetrahedrons and approximate function by its
! value at tetrahedron "origin" plus gradient.
implicit none
integer :: nsub
double complex :: valvtx(8),valgrid(nsub+1,nsub+1,nsub+1)
double precision :: contragrad(3,3,6),relvtx(3,3,6),coord(3),point(3)
double complex :: grad(3,6)
double precision :: rr(3,8),qq(3),qq2,volume
integer :: ivndx(4,6)
integer :: ii,jj,ix,iy,iz,itet,inside
data rr(1:3,1) /0.0d0,0.0d0,0.0d0/
data rr(1:3,2) /1.0d0,0.0d0,0.0d0/
data rr(1:3,3) /0.0d0,1.0d0,0.0d0/
data rr(1:3,4) /1.0d0,1.0d0,0.0d0/
data rr(1:3,5) /0.0d0,0.0d0,1.0d0/
data rr(1:3,6) /1.0d0,0.0d0,1.0d0/
data rr(1:3,7) /0.0d0,1.0d0,1.0d0/
data rr(1:3,8) /1.0d0,1.0d0,1.0d0/
data ivndx(1:4,1) /1,2,3,5/
data ivndx(1:4,2) /3,5,6,7/
data ivndx(1:4,3) /2,3,5,6/
data ivndx(1:4,4) /2,3,4,6/
data ivndx(1:4,5) /3,4,6,7/
data ivndx(1:4,6) /4,6,7,8/

do itet=1,6
! vertex coordinates relative to tetrahedron "origin"
  relvtx(:,1,itet)=rr(:,ivndx(2,itet))-rr(:,ivndx(1,itet))
  relvtx(:,2,itet)=rr(:,ivndx(3,itet))-rr(:,ivndx(1,itet))
  relvtx(:,3,itet)=rr(:,ivndx(4,itet))-rr(:,ivndx(1,itet))
! find contragradient vector
  call cross(relvtx(:,2,itet),relvtx(:,3,itet),contragrad(:,1,itet))
  call cross(relvtx(:,3,itet),relvtx(:,1,itet),contragrad(:,2,itet))
  call cross(relvtx(:,1,itet),relvtx(:,2,itet),contragrad(:,3,itet))
  volume=dot_product(contragrad(:,1,itet),relvtx(:,1,itet))
  contragrad(:,:,itet)=contragrad(:,:,itet)/volume
! find gradient
  do ii=1,3
    grad(ii,itet)=(0.d0,0.d0)
    do jj=1,3
      grad(ii,itet)=grad(ii,itet)+contragrad(ii,jj,itet)*(valvtx(ivndx(jj+1,itet))-valvtx(ivndx(1,itet)))
    enddo
  enddo
!  write(42,'(3es10.3,5x,es20.13)') grad(:,itet),sqrt(dot_product(grad(:,itet),grad(:,itet)))
enddo

do ix=1,nsub+1
do iy=1,nsub+1
do iz=1,nsub+1
  coord=(/dble(ix-1),dble(iy-1),dble(iz-1)/)/dble(nsub)
  do itet=1,6
    point=coord-rr(:,ivndx(1,itet))
    call tbound(relvtx(:,:,itet),contragrad(:,:,itet),point,inside)
    if (inside.eq.1) then
      valgrid(ix,iy,iz)=valvtx(ivndx(1,itet))
      do ii=1,3
        valgrid(ix,iy,iz)=valgrid(ix,iy,iz)+(coord(ii)-rr(ii,ivndx(1,itet)))*grad(ii,itet)
      enddo
      exit
    endif
  enddo
enddo
enddo
enddo

return
end subroutine subcubec

!**************************************************************************

subroutine tbound(vertex,contragrad,point,inside)
! checks if point is inside tetrahedron defined by the origin and 
! the three points given by vertex array.
! contragrad(:,1)=vertex(:,2) cross vertex(:,3)/volume
! contragrad(:,2)=vertex(:,3) cross vertex(:,1)/volume
! contragrad(:,3)=vertex(:,1) cross vertex(:,2)/volume
! volume = vertex(:,1) dot contragrad(:,1)
implicit none
double precision :: vertex(3,3),contragrad(3,3),point(3)
integer :: inside
double precision :: pointprime(3)
integer :: ii,jj

! perform transformation to oblique coordinates with vectors of array vertex
! as basis vectors
pointprime=(/0,0,0/)
do ii=1,3
do jj=1,3
  pointprime(ii)=pointprime(ii)+contragrad(jj,ii)*point(jj)
enddo
enddo

! check if transformed coordinates of point are inside unit tetrahedron
inside=1
! x+y+z=1 defines a plane through vertices of unit tetrahedron
if (pointprime(1)+pointprime(2)+pointprime(3).gt.1.d0) inside=0 
if (pointprime(1).lt.0.d0) inside=0
if (pointprime(2).lt.0.d0) inside=0
if (pointprime(3).lt.0.d0) inside=0

return
end subroutine tbound

!**************************************************************************

subroutine subint(mass,ikk,igb,xk,iqq,shiftk,shiftkq,lossv,ww,iwpt,bmet, &
& iqndx,ngkpt,nwpt,nqpt,ek,ikndxq,vol,ipw1,ipwlf,npwlf,pi,rlr,abr,tseincr)
implicit none
integer :: ikk(3),igb(3),iqq(3),ngkpt(3),nwpt,nqpt,iwpt,npwlf,ipw1
integer :: ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ipwlf(3,npwlf)
double complex :: tseincr
double precision :: bmet(3,3),rlr,abr,pi,ek,shiftk(3),shiftkq(3),ww
double precision :: mass,vol,xk(3)
integer :: iks(3),ikslf(3),jka(3),jkk(3),igg(3)
double precision :: xck(3),xckq(3),qq(3),qs(3),kmq(3)
integer :: nsub,nlevel
parameter (nsub=2,nlevel=20)
double complex :: lossv(8)
double precision :: sidelength(3),origin(3)
double precision :: q2grid(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: qgrid(3,nsub+1,nsub+1,nsub+1,nlevel)
double precision :: kq2grid(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: kqgrid(3,nsub+1,nsub+1,nsub+1,nlevel)
double precision :: egrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: lossgrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: xgrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: valgrid(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: evert(8)
double complex :: xvert(8),xincr
integer :: ii,jj,ix,iy,iz,ilevel,jx,jy,jz,lmnt(nlevel)
double complex :: xsum,xsumold,xint
double precision :: diff

do ii=1,3
  iks(ii)=iqq(ii)-ngkpt(ii)/2
enddo
if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
  call fqq(shiftk,shiftkq,iks,ngkpt,qq,qs)
else
  call fqq(shiftk,shiftk,iks,ngkpt,qq,qs)
endif
qq=qq+dble(ipwlf(:,ipw1))
kmq=xk-qq

ilevel=1
xsumold=tseincr
xint=(0.d0,0.d0)

300 continue
sidelength=1.d0/dble(ngkpt*nsub**(ilevel-1))
call subgridq2(nsub,sidelength,qq,bmet,q2grid(:,:,:,ilevel),qgrid(:,:,:,:,ilevel))
call subcubec(lossv,nsub,lossgrid(:,:,:,ilevel))
xgrid(:,:,:,ilevel)=lossgrid(:,:,:,ilevel)*4.d0*pi/q2grid(:,:,:,ilevel)
call subgridq2(nsub,-sidelength,kmq,bmet,kq2grid(:,:,:,ilevel),kqgrid(:,:,:,:,ilevel))
egrid(:,:,:,ilevel)=ek-kq2grid(:,:,:,ilevel)/(2.d0*mass)
!write(42,'(4f10.3)') ek,kq2grid(1,1,1,ilevel)/2.d0,egrid(1,1,1,ilevel),ww
!write(42,'("qq = ",3es12.3)') qq
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*) 'loss'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(1:2)),dimag(lossv(1:2))
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(3:4)),dimag(lossv(3:4))
!write(42,*)
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(5:6)),dimag(lossv(5:6))
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(7:8)),dimag(lossv(7:8))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,1,ilevel)),dimag(lossgrid(1:3,1,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,1,ilevel)),dimag(lossgrid(1:3,2,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,1,ilevel)),dimag(lossgrid(1:3,3,1,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,2,ilevel)),dimag(lossgrid(1:3,1,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,2,ilevel)),dimag(lossgrid(1:3,2,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,2,ilevel)),dimag(lossgrid(1:3,3,2,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,3,ilevel)),dimag(lossgrid(1:3,1,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,3,ilevel)),dimag(lossgrid(1:3,2,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,3,ilevel)),dimag(lossgrid(1:3,3,3,ilevel))
!write(42,*)'------------------------------------------------------'
!write(42,*) 'energy'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2(3x,3es12.3))') egrid(1:3,1,1,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,2,1,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,3,1,ilevel)-ww
!write(42,*)
!write(42,'(2(3x,3es12.3))') egrid(1:3,1,2,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,2,2,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,3,2,ilevel)-ww
!write(42,*)
!write(42,'(2(3x,3es12.3))') egrid(1:3,1,3,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,2,3,ilevel)-ww
!write(42,'(2(3x,3es12.3))') egrid(1:3,3,3,ilevel)-ww
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*)'q'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,1,1,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,2,1,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,3,1,ilevel),jj=1,3)
!write(42,*)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,1,2,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,2,2,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,3,2,ilevel),jj=1,3)
!write(42,*)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,1,3,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,2,3,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (qgrid(1:3,jj,3,3,ilevel),jj=1,3)
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,1,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,1,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,1,ilevel)
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,2,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,2,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,2,ilevel)
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,3,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,3,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,3,ilevel)
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*)'k-q'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,1,1,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,2,1,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,3,1,ilevel),jj=1,3)
!write(42,*)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,1,2,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,2,2,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,3,2,ilevel),jj=1,3)
!write(42,*)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,1,3,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,2,3,ilevel),jj=1,3)
!write(42,'(3(2x,3f8.3))') (kqgrid(1:3,jj,3,3,ilevel),jj=1,3)
!write(42,*)
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,1,1,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,2,1,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,3,1,ilevel)/2
!write(42,*)
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,1,2,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,2,2,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,3,2,ilevel)/2
!write(42,*)
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,1,3,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,2,3,ilevel)/2
!write(42,'(2(3x,3es12.3))') kq2grid(1:3,3,3,ilevel)/2
!write(42,*)
!write(42,*) ek-ww
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*)'4 pi loss / q^2'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,1,ilevel)),dimag(xgrid(1:3,1,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,1,ilevel)),dimag(xgrid(1:3,2,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,1,ilevel)),dimag(xgrid(1:3,3,1,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,2,ilevel)),dimag(xgrid(1:3,1,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,2,ilevel)),dimag(xgrid(1:3,2,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,2,ilevel)),dimag(xgrid(1:3,3,2,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,3,ilevel)),dimag(xgrid(1:3,1,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,3,ilevel)),dimag(xgrid(1:3,2,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,3,ilevel)),dimag(xgrid(1:3,3,3,ilevel))
!write(42,*)

xsum=0.d0
do ix=1,nsub
do iy=1,nsub
do iz=1,nsub
  evert(1)=egrid(ix  ,iy  ,iz  ,ilevel)
  evert(2)=egrid(ix+1,iy  ,iz  ,ilevel)
  evert(3)=egrid(ix  ,iy+1,iz  ,ilevel)
  evert(4)=egrid(ix+1,iy+1,iz  ,ilevel)
  evert(5)=egrid(ix  ,iy  ,iz+1,ilevel)
  evert(6)=egrid(ix+1,iy  ,iz+1,ilevel)
  evert(7)=egrid(ix  ,iy+1,iz+1,ilevel)
  evert(8)=egrid(ix+1,iy+1,iz+1,ilevel)
  xvert(1)=xgrid(ix  ,iy  ,iz  ,ilevel)
  xvert(2)=xgrid(ix+1,iy  ,iz  ,ilevel)
  xvert(3)=xgrid(ix  ,iy+1,iz  ,ilevel)
  xvert(4)=xgrid(ix+1,iy+1,iz  ,ilevel)
  xvert(5)=xgrid(ix  ,iy  ,iz+1,ilevel)
  xvert(6)=xgrid(ix+1,iy  ,iz+1,ilevel)
  xvert(7)=xgrid(ix  ,iy+1,iz+1,ilevel)
  xvert(8)=xgrid(ix+1,iy+1,iz+1,ilevel)
  call cubeint(evert,ww,xvert,xincr)
  valgrid(ix,iy,iz,ilevel)=xincr
  xsum=xsum+xincr/dble(nsub)**3
!write(42,*)
!write(42,'(3i3)') ix,iy,iz
!write(42,*)
!write(42,'(2(3x,2es12.3))') evert (1:2)-ww,dble(xvert(1:2))
!write(42,'(2(3x,2es12.3))') evert (3:4)-ww,dble(xvert(3:4))
!write(42,*)
!write(42,'(2(3x,2es12.3))') evert (5:6)-ww,dble(xvert(5:6))
!write(42,'(2(3x,2es12.3))') evert (7:8)-ww,dble(xvert(7:8))
!write(42,*)
!write(42,'(es12.3)') dble(xincr)
!write(42,*)
!!  write(6,'(2es12.3)') xsum,xincr/dble(nsub)**3
enddo
enddo
enddo
!write(42,*)'------------------------------------------------------'
!write(42,*)'integrated values'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,1,1,ilevel)),dimag(valgrid(1:2,1,1,ilevel))
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,2,1,ilevel)),dimag(valgrid(1:2,2,1,ilevel))
!write(42,*)
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,1,2,ilevel)),dimag(valgrid(1:2,1,2,ilevel))
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,2,2,ilevel)),dimag(valgrid(1:2,2,2,ilevel))
!write(42,*)
diff=abs(xsum-xsumold)
!write(42,*) '       ****** New Value ******    ******* Old Value *******   |difference|'
!write(42,'(2i3,3(2es12.3,3x))') ilevel-1,lmnt(ilevel-1),xsum,xsumold,diff
if (diff.le.rlr*abs(xsum).or.diff.le.abr) then
  ilevel=ilevel-1
!  xint=xint+xsum/(nsub**(3*ilevel))
  xint=xint+xsum/(dble(nsub**ilevel)**3)
!write(42,'(2i3,5x,3i3,a,2es12.3)') ilevel,lmnt(ilevel),jx,jy,jz,'!',xsum,xint
  if (ilevel.le.0) then
    tseincr=xint
    return
  else
    goto 400
  endif
endif
lmnt(ilevel)=1
100 continue
  jx=mod(lmnt(ilevel)-1,nsub)+1
  jy=mod((lmnt(ilevel)-1)/nsub,nsub)+1
  jz=mod((lmnt(ilevel)-1)/(nsub**2),nsub)+1
  xsumold=valgrid(jx,jy,jz,ilevel)
  lossv(1)=lossgrid(jx  ,jy  ,jz  ,ilevel)
  lossv(2)=lossgrid(jx+1,jy  ,jz  ,ilevel)
  lossv(3)=lossgrid(jx  ,jy+1,jz  ,ilevel)
  lossv(4)=lossgrid(jx+1,jy+1,jz  ,ilevel)
  lossv(5)=lossgrid(jx  ,jy  ,jz+1,ilevel)
  lossv(6)=lossgrid(jx+1,jy  ,jz+1,ilevel)
  lossv(7)=lossgrid(jx  ,jy+1,jz+1,ilevel)
  lossv(8)=lossgrid(jx+1,jy+1,jz+1,ilevel)
  kmq=kqgrid(:,jx,jy,jz,ilevel)
  qq=qgrid(:,jx,jy,jz,ilevel)
!write(42,'(10(2i3,3x))') (ii,lmnt(ii),ii=ilevel,1,-1)
!write(42,'(2i3,5x,3i3,5x,es12.3)') ilevel,lmnt(ilevel),jx,jy,jz,xsumold
  ilevel=ilevel+1
  if (ilevel.gt.nlevel) then
    write(6,*) 'subint: too many subdivisions'
    stop
  endif
  goto 300
  400 continue
  lmnt(ilevel)=lmnt(ilevel)+1
  if (lmnt(ilevel).gt.nsub**3) then
    ilevel=ilevel-1
    if (ilevel.le.0) then
      tseincr=xint
      return
    endif
    goto 400
  endif
goto 100

return
end subroutine subint

!**************************************************************************

subroutine subint2(iqq,shiftk,shiftkq,lossv,enval,ww,iwpt,bmet,iqndx,ngkpt, &
& nwpt,nqpt,ikndxq,vol,ipw1,ipw2,ipwlf,npwlf,pi,rlr,abr,tseincr)
implicit none
integer :: iqq(3),ngkpt(3),nwpt,nqpt,iwpt,npwlf,ipw1,ipw2
integer :: ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ipwlf(3,npwlf)
double complex :: tseincr
double precision :: bmet(3,3),rlr,abr,pi,ek,shiftk(3),shiftkq(3),ww
double precision :: vol,xk(3)
integer :: iks(3),ikslf(3),jka(3),jkk(3),igg(3)
double precision :: xck(3),xckq(3),qq(3),qs(3),kmq(3)
integer :: nsub,nlevel
parameter (nsub=2,nlevel=20)
double precision :: enval(8)
double complex :: lossv(8)
double precision :: sidelength(3),origin(3)
double precision :: q2grid(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: q2gridA(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: q2gridB(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: qgrid(3,nsub+1,nsub+1,nsub+1,nlevel)
double precision :: engrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: lossgrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: xgrid(nsub+1,nsub+1,nsub+1,nlevel)
double complex :: valgrid(nsub+1,nsub+1,nsub+1,nlevel)
double precision :: evert(8)
double complex :: xvert(8),xincr
integer :: ii,jj,ix,iy,iz,ilevel,jx,jy,jz,lmnt(nlevel)
double complex :: xsum,xsumold,xint
double precision :: diff

do ii=1,3
  iks(ii)=iqq(ii)-ngkpt(ii)/2
enddo
!write(42,*) "iqq = ",iqq
!write(42,*) "iks = ",iks
!write(42,*) "shiftk = ",shiftk
!write(42,*) "shiftkq = ",shiftkq
!write(42,*) "ipw1 = ",ipwlf(:,ipw1)
!write(42,*) "ipw2 = ",ipwlf(:,ipw2)
if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
  call fqq(shiftk,shiftkq,iks,ngkpt,qq,qs)
else
  call fqq(shiftk,shiftk,iks,ngkpt,qq,qs)
endif
kmq=xk-qq

ilevel=1
xsumold=tseincr
xint=0.d0

300 continue
sidelength=1.d0/dble(ngkpt*nsub**(ilevel-1))
call subgridq2(nsub,sidelength,qq+dble(ipwlf(:,ipw1)),bmet,q2gridA(:,:,:,ilevel),qgrid(:,:,:,:,ilevel))
call subgridq2(nsub,sidelength,qq+dble(ipwlf(:,ipw2)),bmet,q2gridB(:,:,:,ilevel),qgrid(:,:,:,:,ilevel))
q2grid(:,:,:,ilevel)=sqrt(q2gridA(:,:,:,ilevel)*q2gridB(:,:,:,ilevel))
call subcube(enval,nsub,engrid(:,:,:,ilevel))
call subcubec(lossv,nsub,lossgrid(:,:,:,ilevel))
xgrid(:,:,:,ilevel)=lossgrid(:,:,:,ilevel)*4.d0*pi/q2grid(:,:,:,ilevel)
!!write(42,'(4f10.3)') ek,kq2grid(1,1,1,ilevel)/2.d0,egrid(1,1,1,ilevel),ww
!write(42,'(2i3,4x,"qq = ",3es12.3)') ilevel-1,lmnt(ilevel-1),qq
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*) 'loss'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(1:2)),dimag(lossv(1:2))
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(3:4)),dimag(lossv(3:4))
!write(42,*)
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(5:6)),dimag(lossv(5:6))
!write(42,'(2i3,2(3x,2es12.3))') ilevel-1,lmnt(ilevel-1),dble(lossv(7:8)),dimag(lossv(7:8))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,1,ilevel)),dimag(lossgrid(1:3,1,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,1,ilevel)),dimag(lossgrid(1:3,2,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,1,ilevel)),dimag(lossgrid(1:3,3,1,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,2,ilevel)),dimag(lossgrid(1:3,1,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,2,ilevel)),dimag(lossgrid(1:3,2,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,2,ilevel)),dimag(lossgrid(1:3,3,2,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,1,3,ilevel)),dimag(lossgrid(1:3,1,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,2,3,ilevel)),dimag(lossgrid(1:3,2,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(lossgrid(1:3,3,3,ilevel)),dimag(lossgrid(1:3,3,3,ilevel))
!write(42,*)'------------------------------------------------------'
!write(42,*) 'energy'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2i3,3x,2es12.3)') ilevel-1,lmnt(ilevel-1),enval(1:2)-ww
!write(42,'(2i3,3x,2es12.3)') ilevel-1,lmnt(ilevel-1),enval(3:4)-ww
!write(42,*)
!write(42,'(2i3,3x,2es12.3)') ilevel-1,lmnt(ilevel-1),enval(5:6)-ww
!write(42,'(2i3,3x,2es12.3)') ilevel-1,lmnt(ilevel-1),enval(7:8)-ww
!write(42,*)
!write(42,'(2(3x,3es12.3))') engrid(1:3,1,1,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,2,1,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,3,1,ilevel)-ww
!write(42,*)
!write(42,'(2(3x,3es12.3))') engrid(1:3,1,2,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,2,2,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,3,2,ilevel)-ww
!write(42,*)
!write(42,'(2(3x,3es12.3))') engrid(1:3,1,3,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,2,3,ilevel)-ww
!write(42,'(2(3x,3es12.3))') engrid(1:3,3,3,ilevel)-ww
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*)'q^2'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,1,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,1,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,1,ilevel)
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,2,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,2,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,2,ilevel)
!write(42,*)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,1,3,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,2,3,ilevel)
!write(42,'(2(3x,3es12.3))') q2grid(1:3,3,3,ilevel)
!write(42,*)
!write(42,*)'------------------------------------------------------'
!write(42,*)'4 pi loss / q^2'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,1,ilevel)),dimag(xgrid(1:3,1,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,1,ilevel)),dimag(xgrid(1:3,2,1,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,1,ilevel)),dimag(xgrid(1:3,3,1,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,2,ilevel)),dimag(xgrid(1:3,1,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,2,ilevel)),dimag(xgrid(1:3,2,2,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,2,ilevel)),dimag(xgrid(1:3,3,2,ilevel))
!write(42,*)
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,1,3,ilevel)),dimag(xgrid(1:3,1,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,2,3,ilevel)),dimag(xgrid(1:3,2,3,ilevel))
!write(42,'(2(3x,3es12.3))') dble(xgrid(1:3,3,3,ilevel)),dimag(xgrid(1:3,3,3,ilevel))
!write(42,*)

xsum=0.d0
do ix=1,nsub
do iy=1,nsub
do iz=1,nsub
!write(42,*)
!write(42,'(3i3)') ix,iy,iz
!write(42,*)
  evert(1)=engrid(ix  ,iy  ,iz  ,ilevel)
  evert(2)=engrid(ix+1,iy  ,iz  ,ilevel)
  evert(3)=engrid(ix  ,iy+1,iz  ,ilevel)
  evert(4)=engrid(ix+1,iy+1,iz  ,ilevel)
  evert(5)=engrid(ix  ,iy  ,iz+1,ilevel)
  evert(6)=engrid(ix+1,iy  ,iz+1,ilevel)
  evert(7)=engrid(ix  ,iy+1,iz+1,ilevel)
  evert(8)=engrid(ix+1,iy+1,iz+1,ilevel)
  xvert(1)=xgrid(ix  ,iy  ,iz  ,ilevel)
  xvert(2)=xgrid(ix+1,iy  ,iz  ,ilevel)
  xvert(3)=xgrid(ix  ,iy+1,iz  ,ilevel)
  xvert(4)=xgrid(ix+1,iy+1,iz  ,ilevel)
  xvert(5)=xgrid(ix  ,iy  ,iz+1,ilevel)
  xvert(6)=xgrid(ix+1,iy  ,iz+1,ilevel)
  xvert(7)=xgrid(ix  ,iy+1,iz+1,ilevel)
  xvert(8)=xgrid(ix+1,iy+1,iz+1,ilevel)
  call cubeint(evert,ww,xvert,xincr)
  valgrid(ix,iy,iz,ilevel)=xincr
  xsum=xsum+xincr/dble(nsub)**3
!write(42,'(2(3x,2es12.3))') dble(xvert(1:2)),dimag(xvert(1:2))
!write(42,'(2(3x,2es12.3))') dble(xvert(3:4)),dimag(xvert(3:4))
!write(42,*)
!write(42,'(2(3x,2es12.3))') dble(xvert(5:6)),dimag(xvert(5:6))
!write(42,'(2(3x,2es12.3))') dble(xvert(7:8)),dimag(xvert(7:8))
!write(42,*)
!write(42,*) '------------------------------------------------------------'
!write(42,*)
!write(42,'(2(3x,2es12.3))') evert (1:2)-ww
!write(42,'(2(3x,2es12.3))') evert (3:4)-ww
!write(42,*)
!write(42,'(2(3x,2es12.3))') evert (5:6)-ww
!write(42,'(2(3x,2es12.3))') evert (7:8)-ww
!write(42,'(es12.3)') dble(xincr)
!write(42,*)
!!write(42,'(82es12.3)') evert-ww
!!  write(6,'(2es12.3)') xsum,xincr/dble(nsub)**3
enddo
enddo
enddo
!!write(6,'(es10.3)') xsumold,xsum
!write(42,*)'------------------------------------------------------'
!write(42,*)'integrated values'
!write(42,*)'------------------------------------------------------'
!write(42,*)
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,1,1,ilevel)),dimag(valgrid(1:2,1,1,ilevel))
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,2,1,ilevel)),dimag(valgrid(1:2,2,1,ilevel))
!write(42,*)
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,1,2,ilevel)),dimag(valgrid(1:2,1,2,ilevel))
!write(42,'(10x,2(2es12.3,3x))') dble(valgrid(1:2,2,2,ilevel)),dimag(valgrid(1:2,2,2,ilevel))
!write(42,*)
diff=abs(xsum-xsumold)
!write(42,*) '       ****** New Value ******    ******* Old Value *******   |difference|'
!write(42,'(2i3,3(2es12.3,3x))') ilevel-1,lmnt(ilevel-1),xsum,xsumold,diff
if (diff.le.rlr*abs(xsum).or.diff.le.abr) then
  ilevel=ilevel-1
  xint=xint+xsum/(dble(nsub**ilevel)**3)
!write(42,'(2i3,5x,3i3,a,2(2es12.3,4x))') ilevel,lmnt(ilevel),jx,jy,jz,'!',xsum,xint
!!write(42,'(2i15)') nsub,nsub**(3*ilevel)
  if (ilevel.le.0) then
    tseincr=xint
    return
  else
    goto 400
  endif
endif
lmnt(ilevel)=1
100 continue
  jx=mod(lmnt(ilevel)-1,nsub)+1
  jy=mod((lmnt(ilevel)-1)/nsub,nsub)+1
  jz=mod((lmnt(ilevel)-1)/(nsub**2),nsub)+1
  xsumold=valgrid(jx,jy,jz,ilevel)
  lossv(1)=lossgrid(jx  ,jy  ,jz  ,ilevel)
  lossv(2)=lossgrid(jx+1,jy  ,jz  ,ilevel)
  lossv(3)=lossgrid(jx  ,jy+1,jz  ,ilevel)
  lossv(4)=lossgrid(jx+1,jy+1,jz  ,ilevel)
  lossv(5)=lossgrid(jx  ,jy  ,jz+1,ilevel)
  lossv(6)=lossgrid(jx+1,jy  ,jz+1,ilevel)
  lossv(7)=lossgrid(jx  ,jy+1,jz+1,ilevel)
  lossv(8)=lossgrid(jx+1,jy+1,jz+1,ilevel)
  enval(1)=engrid(jx  ,jy  ,jz  ,ilevel)
  enval(2)=engrid(jx+1,jy  ,jz  ,ilevel)
  enval(3)=engrid(jx  ,jy+1,jz  ,ilevel)
  enval(4)=engrid(jx+1,jy+1,jz  ,ilevel)
  enval(5)=engrid(jx  ,jy  ,jz+1,ilevel)
  enval(6)=engrid(jx+1,jy  ,jz+1,ilevel)
  enval(7)=engrid(jx  ,jy+1,jz+1,ilevel)
  enval(8)=engrid(jx+1,jy+1,jz+1,ilevel)
  qq=qgrid(:,jx,jy,jz,ilevel)-dble(ipwlf(:,ipw2))
!write(42,'(10(2i3,3x))') (ii,lmnt(ii),ii=ilevel,1,-1)
!write(42,'(2i3,5x,3i3,5x,es12.3)') ilevel,lmnt(ilevel),jx,jy,jz,xsumold
  ilevel=ilevel+1
  if (ilevel.gt.nlevel) then
    write(6,*) 'subint2: too many subdivisions'
    stop
  endif
  goto 300
  400 continue
  lmnt(ilevel)=lmnt(ilevel)+1
  if (lmnt(ilevel).gt.nsub**3) then
    ilevel=ilevel-1
    if (ilevel.le.0) then
      tseincr=xint
      return
    endif
    goto 400
  endif
goto 100

return
end subroutine subint2
