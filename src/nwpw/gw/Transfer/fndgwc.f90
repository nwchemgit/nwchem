subroutine fndgwc(bindir,nbcore,ncband,eigen,enrgy,bantot, &
&   kpt,nkpt,indxkbnd)
implicit none
character(80) :: bindir
integer :: nbcore,ncband,ngwb,bantot
integer :: nkpt
double precision :: eigen(bantot),kpt(3,nkpt),enrgy(bantot)
integer :: indxkbnd(nkpt)
double precision :: gwc(bantot)
double precision, allocatable :: xkgrid(:),gwa(:,:,:,:), &
& fx(:,:,:),fy(:,:,:),fxy(:,:,:)
integer :: ix,iy,iz,ix1,iy1,iz1,ilo,ihi,jlo,jhi,klo,khi,ib,npts, &
&          ii,ieval,ik,nbrm
integer :: nkc1,ist1,ind1
double precision :: xk,yk,zk,xk1,yk1,zk1,gw,t,gw0,dfv
double precision :: aa(0:3),dx
double precision :: ryd,hart
parameter (ryd=13.6056981,hart=2.d0*ryd)

open (unit=10,file=trim(bindir)//'gw.dat',status='old')
read(10,*) nkc1,ist1,ind1,ngwb
allocate(gwa(ncband-nbcore,ist1:ind1,ist1:ind1,ist1:ind1),xkgrid(ist1:ind1))
do ix=ist1,ind1
do iy=ist1,ind1
do iz=ist1,ind1
  nbrm=min(ngwb,ncband-nbcore)
  read(10,*) ix1,iy1,iz1,gwa(1:nbrm,ix,iy,iz)
  if (nbrm.lt.ncband-nbcore) then
    do ib=nbrm+1,ncband-nbcore
      gwa(ib,ix,iy,iz)=gwa(nbrm,ix,iy,iz)
    enddo
  endif
enddo
enddo
enddo

do ix=ist1,ind1
  xkgrid(ix)=dble(ix)/dble(nkc1)
enddo
t=1.d0

allocate(fx(ist1:ind1,ist1:ind1,ist1:ind1),fy(ist1:ind1,ist1:ind1,ist1:ind1), &
& fxy(ist1:ind1,ist1:ind1,ist1:ind1))

do ib=1,ncband-nbcore
  call trispline(gwa(ib,:,:,:),xkgrid,xkgrid,xkgrid,nkc1,nkc1,nkc1, &
&                t,t,t,fx,fy,fxy)

!  write(6,'(i2)') ib
!  write(72,'(i2)') ib
!  write(72,*)
enddo

do ib=1,ncband-nbcore
  call trispline(gwa(ib,:,:,:),xkgrid,xkgrid,xkgrid,nkc1,nkc1,nkc1, &
&                t,t,t,fx,fy,fxy)

!  write(6,'(i2)') ib
!  write(72,'(i2)') ib
  do ik=1,nkpt
    ii=indxkbnd(ik)+ib+nbcore
    xk=kpt(1,ik)
    yk=kpt(2,ik)
    zk=kpt(3,ik)
    call spl23(gwa(ib,:,:,:),fx,fy,fxy,xkgrid,xkgrid,xkgrid,nkc1,nkc1,nkc1,1.d0,1.d0,1.d0,xk,yk,1,zk,gw)
    gwc(ii)=gw/hart
    enrgy(ii)=eigen(ii)+gw/hart
!    write(72,'(2x,i4,3f12.6)') ik,gwc(ik,ib)*hart,enrgy(ik,ib)*hart,en(ik,ib)*hart
  enddo
!  write(72,*)
enddo

return
end subroutine fndgwc
