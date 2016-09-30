!*************************************************************************

subroutine mkwf(rr,iband,ikpt,ncg,nkpt,npwt,cg,kpt,indxkcg,indxkpw,npwarr,kg,wf)
implicit none
integer :: ncg,iband,ikpt,nkpt,npwt
double complex :: cg(ncg)
double precision :: kpt(3,nkpt),rr(3),xk(3),pi
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: wf
double complex :: phase
!double precision :: xnorm
integer :: ii,ig,icg
parameter (pi=3.1415926535897932384626433832795d0)

wf=(0.d0,0.d0)
!xnorm=0.d0
do ig=1,npwarr(ikpt)
  do ii=1,3
    xk(ii)=kpt(ii,ikpt)+kg(ii,indxkpw(ikpt)+ig)
  enddo
  phase=2.d0*pi*(rr(1)*xk(1)+rr(2)*xk(2)+rr(3)*xk(3))*(0.d0,1.d0)
  icg=indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ig
  wf=wf+cg(icg)*exp(phase)
!  xnorm=xnorm+dble(cg(icg)*conjg(cg(icg)))
!  write(6,*) ig,wf
!  write(6,'(3i3,3(3x,2es10.3))') kg(1:3,indxkpw(ikpt)+ig), &
!&  cg(icg),wf,xnorm
enddo

return
end subroutine mkwf

!*************************************************************************


