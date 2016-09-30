program testint
implicit none
double precision :: enval(8),ww
double complex :: fv(8)
double complex :: xint(9)
double precision :: vtx(3,8),grad(3),agrad,f0,xgrad(3,6),fact
integer :: ii
data vtx(1:3,1) /0.0d0,0.0d0,0.0d0/
data vtx(1:3,2) /1.0d0,0.0d0,0.0d0/
data vtx(1:3,3) /0.0d0,1.0d0,0.0d0/
data vtx(1:3,4) /1.0d0,1.0d0,0.0d0/
data vtx(1:3,5) /0.0d0,0.0d0,1.0d0/
data vtx(1:3,6) /1.0d0,0.0d0,1.0d0/
data vtx(1:3,7) /0.0d0,1.0d0,1.0d0/
data vtx(1:3,8) /1.0d0,1.0d0,1.0d0/

fv=(/((1.d0,0.d0),ii=1,8)/)
fact=1.d0

!enval=(/2.710E-07,2.661E-07,5.861E-09,9.025E-10,2.661E-07,2.611E-07,9.025E-10,-1.846E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(1))
!write(6,*) xint(1)
!
!enval=(/2.661E-07,2.611E-07,9.025E-10,-1.846E-07,2.611E-07,7.565E-08,-4.056E-09,-1.895E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(2))
!write(6,*) xint(2)
!
!enval=(/ 5.861E-09,  9.025E-10, -2.593E-07, -2.643E-07,  9.025E-10, -1.846E-07, -2.643E-07, -4.497E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(3))
!write(6,*) xint(3)
!
!enval=(/9.025E-10, -1.846E-07, -2.643E-07, -4.497E-07, -4.056E-09, -1.895E-07, -2.692E-07, -4.547E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(4))
!write(6,*) xint(4)
!
!enval=(/2.661E-07,  2.611E-07,  9.025E-10, -4.056E-09,  2.611E-07,  7.565E-08, -1.846E-07, -1.895E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(5))
!write(6,*) xint(5)
!
!enval=(/2.611E-07,  7.565E-08, -1.846E-07, -1.895E-07,  7.565E-08, -1.098E-07, -1.895E-07, -3.750E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(6))
!write(6,*) xint(6)
!
!enval=(/9.025E-10, -4.056E-09, -2.643E-07, -2.692E-07, -1.846E-07, -1.895E-07, -4.497E-07, -4.547E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(7))
!write(6,*) xint(7)
!
!enval=(/-1.846E-07, -1.895E-07, -4.497E-07, -4.547E-07, -1.895E-07, -3.750E-07, -4.547E-07, -6.401E-07/)*fact
!call cubeint(enval,0.d0,fv,xint(8))
!write(6,*) xint(8)

enval=(/2.710E-07,2.611E-07,-2.593E-07,-2.692E-07,2.611E-07,-1.098E-07,-2.692E-07,-6.401E-07/)*fact
call getgrad(enval,xgrad)
do ii=1,6
  write(6,'(3es12.3)') xgrad(:,ii)
enddo
call cubeint(enval,0.d0,fv,xint(9))
write(6,*) xint(9)
!write(6,*) sum(xint(1:8))/xint(3)

!!f0=-1.0000001d0
f0=enval(1)
!!grad=(/1.d0,1.d0,1.d0/)
grad=xgrad(:,1)
agrad=sqrt(dot_product(grad,grad))
enval=(/(dot_product(grad,vtx(:,ii))+f0,ii=1,8)/)
!
!!enval=(/(0.5d0*(-1)**((ii-1)/1),ii=1,8)/)
!!enval=(/(0.9d0,ii=1,4),(-0.1d0,ii=1,4)/)
write(6,*)
write(6,*) agrad
write(6,*)
write(6,'(2(3x,2es12.3))') enval(1:2)
write(6,'(2(3x,2es12.3))') enval(3:4)
write(6,*)
write(6,'(2(3x,2es12.3))') enval(5:6)
write(6,'(2(3x,2es12.3))') enval(7:8)
write(6,*)
!call cubeint(enval,0.d0,fv,xint(1))
!write(6,*) xint(1)
!write(6,*) xint(1)*agrad
!
stop
end program testint

!------------------------------------------------------------------------

subroutine getgrad(valvtx,grad)
implicit none
double precision :: valvtx(8)
double precision :: contragrad(3,3,6),relvtx(3,3,6),grad(3,6),coord(3),point(3)
double precision :: rr(3,8),volume
integer :: ivndx(4,6)
integer :: ii,jj,ix,iy,iz,itet
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
data ivndx(1:4,5) /3,4,6,8/
data ivndx(1:4,6) /3,6,7,8/

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

return
end subroutine getgrad
