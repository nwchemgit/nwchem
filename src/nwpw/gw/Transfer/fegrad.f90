subroutine fegrad(enrgy,bantot,rprimd,ngkpt,kpt,nkpt,shiftk,ikndx,indxkbnd,ncband,egrad)
implicit none
integer :: bantot,ngkpt(3),nkpt,ncband
integer :: indxkbnd(nkpt),ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: twopi
double precision :: enrgy(bantot),rprimd(3,3),shiftk(3),kpt(3,nkpt)
double precision :: egrad(3,bantot)
double precision :: engrid(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fx(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fy(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fz(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fxy(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fxz(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fyz(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: fxyz(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: gridx(ngkpt(1)),gridy(ngkpt(2)),gridz(ngkpt(3))
double precision :: xx,yy,zz,grk(3)
integer :: ix,iy,iz,ii,jj,ie,ikk(3),iband,ikpt

! Unit cartesian vectors x(ii,1), x(ii,2), x(ii,3)
! Reciprocal lattice vector jj G(ii,jj)
! Direct lattice vector jj R(ii,jj)
! Wave vector expressed in either cartesian coordiantes or reciprocal lattice coordinates
! k(ii) = x(ii,1) p(1) + x(ii,2) p(2) + x(ii,3) p(3)
!       = G(ii,1) q(1) + G(ii,2) q(2) + G(ii,3) q(3)
! sum_ii R(ii,jj)*G(ii,kk) = 2 pi delta(jj,kk)
! implies q(jj) = sum_ii R(ii,jj)*k(ii)/(2 pi)
! which implies (partial/partial p(ii)) q(jj) = sum_kk x(kk,ii)*R(kk,jj)/(2 pi)
!                                             = R(ii,jj)/(2 pi)
! Energy known on q grid E(q(1),q(2),q(3))
! Want cartesian gradient of E
! grad E(ii) = (partial/partial p(ii)) E(q(1),q(2),q(3))
! = sum_jj [(partial/partial p(ii)) q(jj)] (partial/partial q(jj)) E(q(1),q(2),q(3))
! = sum_jj R(ii,jj) (partial/partial q(jj)) E(q(1),q(2),q(3)) / (2 pi)

twopi=2.d0*acos(-1.d0)
do ix=1,ngkpt(1)
  xx=dble(ix)/dble(ngkpt(1))+shiftk(1)/dble(ngkpt(1))
  xx=mod(xx,1.d0)
  if (xx.lt.0.d0) xx=xx+1.d0
  xx=1.d0-mod(1.d0-xx,1.d0)
  xx=xx-0.5d0
  gridx(ix)=xx
enddo
do iy=1,ngkpt(2)
  yy=dble(iy)/dble(ngkpt(2))+shiftk(2)/dble(ngkpt(2))
  yy=mod(yy,1.d0)
  if (yy.lt.0.d0) yy=yy+1.d0
  yy=1.d0-mod(1.d0-yy,1.d0)
  yy=yy-0.5d0
  gridy(iy)=yy
enddo
do iz=1,ngkpt(3)
  zz=dble(iz)/dble(ngkpt(3))+shiftk(3)/dble(ngkpt(3))
  zz=mod(zz,1.d0)
  if (zz.lt.0.d0) zz=zz+1.d0
  zz=1.d0-mod(1.d0-zz,1.d0)
  zz=zz-0.5d0
  gridz(iz)=zz
enddo
do iband=1,ncband
  do ix=1,ngkpt(1)
  do iy=1,ngkpt(2)
  do iz=1,ngkpt(3)
    ikpt=ikndx(ix,iy,iz)
    engrid(ix,iy,iz)=enrgy(indxkbnd(ikpt)+iband)
  enddo
  enddo
  enddo
  call xyzspline(engrid,gridx,gridy,gridz,ngkpt(1),ngkpt(2),ngkpt(3),1.d0,1.d0,1.d0,fx,fy,fz,fxy,fxz,fyz,fxyz)
  do ikpt=1,nkpt
    do ii=1,3
      ikk(ii)=nint((kpt(ii,ikpt)+0.5d0)*dble(ngkpt(ii))-shiftk(ii))
    enddo
    grk=(/fx(ikk(1),ikk(2),ikk(3)),fy(ikk(1),ikk(2),ikk(3)),fz(ikk(1),ikk(2),ikk(3))/)
    ie=indxkbnd(ikpt)+iband
    do ii=1,3
      egrad(ii,ie)=0.d0
      do jj=1,3
        egrad(ii,ie)=egrad(ii,ie)+rprimd(jj,ii)*grk(jj)/twopi
      enddo
    enddo
  enddo
enddo

return
end
