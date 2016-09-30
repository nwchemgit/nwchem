subroutine invertdielf(pola,polb,ipwndx,npwndx,ntpwndx,ipwc,npwc,nwpt,nqpt,lossfn)
implicit none
integer :: npwndx,ntpwndx,npwc,nwpt,nqpt
integer :: ipwndx(2,ntpwndx),ipwc(3,npwc)
double complex :: dielf(npwndx),pola(npwndx),polb(npwndx),lossfn(ntpwndx)
integer :: iqpt,iwpt,ipw1,ipw2,iipw
double complex :: A(npwc,npwc),B(npwc,npwc),AA(npwc,npwc)

  do ipw1=1,npwc
  do ipw2=1,npwc
    A(ipw1,ipw2)=0.d0
  enddo
  enddo
  do iipw=1,npwndx
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    if (ipw1.ne.ipw2) then
      A(ipw1,ipw2)=-pola(iipw)-polb(iipw)
      A(ipw2,ipw1)=-conjg(pola(iipw))+conjg(polb(iipw))
    else
      A(ipw1,ipw2)=1.d0-dble(pola(iipw))-dimag(polb(iipw))*(0.d0,1.d0)
    endif
    dielf(iipw)=A(ipw1,ipw2)
  enddo
  do ipw1=1,npwc
  do ipw2=1,npwc
    AA(ipw1,ipw2)=A(ipw1,ipw2)
  enddo
  enddo
  call inversec(A,B,npwc,npwc)
  do iipw=1,ntpwndx
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    lossfn(iipw)=B(ipw1,ipw2)
  enddo

return
end subroutine invertdielf

!*****************************************************************************

subroutine invertdielfT(pola,polb,ipwndx,npwndx,ntpwndx,ipwc,npwc,nwpt,nqpt,lossfn)
use geometry
implicit none
integer :: npwndx,ntpwndx,npwc,nwpt,nqpt
integer :: npwct
integer :: ipwndx(2,ntpwndx),ipwc(3,npwc)
double complex :: pola(9,npwndx),polb(9,npwndx),lossfn(9,ntpwndx)
integer :: iqpt,iqt,iwpt,ipw1,ipw2,ipw3,iipw,ii,jj,kk,ll,ix,jx
double precision :: blatinv(3,3)
double complex :: dummy(3,3), dummy2(3,3), det
double complex :: A(3*npwc,3*npwc),B(3*npwc,3*npwc),AA(npwc,npwc),BB(npwc,npwc)
double complex :: rllossfn(3,3,npwc,npwc),rlDielf(3,3,npwc,npwc),rlAlpha(3,3,npwc,npwc) ! expressed in reciprocal lattice
double complex :: cartlossfn(3,3,npwc,npwc),cartDielf(3,3,npwc,npwc),cartAlpha(3,3,npwc,npwc) ! expressed in cartesian coordinates
double precision :: prefact
double complex :: testlossfn(3,3,npwc,npwc),testD(3,3,npwc,npwc),testdummy(3,3,npwc,npwc),testB(3,3,npwc,npwc)

! Polarizability in reciprocal lattice
  do ipw1=1,npwc
  do ipw2=1,npwc
    do ii=1,3
    do jj=1,3
      rlAlpha(ii,jj,ipw1,ipw2)=(0.d0,0.d0)
    enddo
    enddo
  enddo
  enddo
  do iipw=1,npwndx
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    do ii=1,3
    do jj=1,3
      iqt=ii+3*(jj-1)
      if (ipw1.ne.ipw1) then
        rlAlpha(ii,jj,ipw1,ipw2)=-pola(iqt,iipw)-polb(iqt,iipw)
        rlAlpha(ii,jj,ipw2,ipw1)=-conjg(pola(iqt,iipw))+conjg(polb(iqt,iipw))
      else
        rlAlpha(ii,jj,ipw1,ipw2)=-dble(pola(iqt,iipw))-dimag(polb(iqt,iipw))*(0.d0,1.d0)
      endif
    enddo
    enddo
  enddo
! Polarizability and dielectric function in cartesian coordiantes
  do ipw1=1,npwc
  do ipw2=1,npwc
    do ii=1,3
    do jj=1,3
      cartAlpha(ii,jj,ipw1,ipw2)=(0.d0,0.d0)
      do kk=1,3
      do ll=1,3
        cartAlpha(ii,jj,ipw1,ipw2)=cartAlpha(ii,jj,ipw1,ipw2)+blat(kk,ii)*rlAlpha(kk,ll,ipw1,ipw2)*blat(ll,jj)
      enddo
      enddo
      cartDielf(ii,jj,ipw1,ipw2)=cartAlpha(ii,jj,ipw1,ipw2)
    enddo
    enddo
    do ii=1,3
      if (ipw1.eq.ipw2) cartDielf(ii,ii,ipw1,ipw2)=cartDielf(ii,ii,ipw1,ipw2)+1.d0
    enddo
  enddo
  enddo
! Invert, lossfn in cartesian coordinates
  do ii=1,3
  do jj=1,3
    do ipw1=1,npwc
    do ipw2=1,npwc
      ix=(ipw1-1)*3+ii
      jx=(ipw2-1)*3+jj
      A(ix,jx)=cartDielf(ii,jj,ipw1,ipw2)
!      AA(ipw1,ipw2)=cartDielf(ii,jj,ipw1,ipw2)
    enddo
    enddo
!    call inversec(AA,BB,npwc,npwc)
  enddo
  enddo
  call inversec(A,B,3*npwc,3*npwc)
  do ii=1,3
  do jj=1,3
    do ipw1=1,npwc
    do ipw2=1,npwc
      ix=(ipw1-1)*3+ii
      jx=(ipw2-1)*3+jj
      cartlossfn(ii,jj,ipw1,ipw2)=B(ix,jx)
!      cartlossfn(ii,jj,ipw1,ipw2)=BB(ipw1,ipw2)
    enddo
    enddo
  enddo
  enddo
! lossfn in reciprocal lattice
  prefact = 1.d0/(2.d0*acos(-1.d0))**2
  do ipw1=1,npwc
  do ipw2=1,npwc
    do ii=1,3
    do jj=1,3
      rllossfn(ii,jj,ipw1,ipw2)=(0.d0,0.d0)
      do kk=1,3
      do ll=1,3
        rllossfn(ii,jj,ipw1,ipw2)=rllossfn(ii,jj,ipw1,ipw2)+rprimd(ii,kk)*cartlossfn(kk,ll,ipw1,ipw2)*rprimd(jj,ll)
      enddo
      enddo
      rllossfn(ii,jj,ipw1,ipw2)=rllossfn(ii,jj,ipw1,ipw2)*prefact
    enddo
    enddo
  enddo
  enddo
  do iipw=1,npwndx
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    do ii=1,3
    do jj=1,3
      iqt=ii+3*(jj-1)
      lossfn(iqt,iipw)=rllossfn(ii,jj,ipw1,ipw2)
    enddo
    enddo
  enddo

return
end subroutine invertdielfT

!*****************************************************************************

subroutine locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,iipw)
implicit none
integer :: nsym,npwc,npwx,ngkpt(3),pwcmin(3),pwcmax(3)
integer :: iqq(3),ipw1,ipw2,isym,ipw3,ipw4,iipw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: pwsymndx(npwc,2*nsym)
integer :: invpw2ndx(npwx,npwx)

isym=iqsymndx(iqq(1),iqq(2),iqq(3))
ipw3=pwsymndx(ipw1,isym)
ipw4=pwsymndx(ipw2,isym)
if (ipw3.ne.0.and.ipw4.ne.0) then
  iipw=invpw2ndx(ipw3,ipw4)
else
  iipw=0
endif

return
end subroutine locateelement

!*****************************************************************************

subroutine mkpwsymndx(ipwc,npwc,syminv,nsym,invpwndx,pwcmin,pwcmax,pwsymndx)
implicit none
integer :: nsym,npwc,pwcmin(3),pwcmax(3)
integer :: syminv(3,3,nsym),ipwc(3,npwc)
integer :: invpwndx(pwcmin(1):pwcmax(1),pwcmin(2):pwcmax(2),pwcmin(3):pwcmax(3))
integer :: pwsymndx(npwc,2*nsym)
integer :: igg(3),iggr(3)
integer :: ii,jj,isym,ipw,isign,iisym

do ipw=1,npwc
  igg=ipwc(:,ipw)
  do isym=1,2*nsym
    iisym=mod(isym-1,nsym)+1
    isign=-(2*((isym-1)/nsym)-1)
    iggr=(/0,0,0/)
    do ii=1,3
    do jj=1,3
      iggr(ii)=iggr(ii)+syminv(jj,ii,iisym)*isign*igg(jj)
    enddo
    enddo
    if (iggr(1).ge.pwcmin(1).and.iggr(1).le.pwcmax(1).and. &
&       iggr(2).ge.pwcmin(2).and.iggr(2).le.pwcmax(2).and. &
&       iggr(3).ge.pwcmin(3).and.iggr(3).le.pwcmax(3)) then
      pwsymndx(ipw,isym)=invpwndx(iggr(1),iggr(2),iggr(3))
    else
      pwsymndx(ipw,isym)=0
    endif
  enddo
enddo

return
end subroutine mkpwsymndx

!*****************************************************************************

subroutine mkinvpwndx(ipwc,npwc,invpwndx,pwcmax,pwcmin)
implicit none
integer :: npwc,pwcmax(3),pwcmin(3)
integer :: ipwc(3,npwc)
integer :: invpwndx(pwcmin(1):pwcmax(1),pwcmin(2):pwcmax(2),pwcmin(3):pwcmax(3))
integer :: ix,iy,iz,ii

do ix=pwcmin(1),pwcmax(1)
do iy=pwcmin(2),pwcmax(2)
do iz=pwcmin(3),pwcmax(3)
  invpwndx(ix,iy,iz)=0
enddo
enddo
enddo

do ii=1,npwc
  invpwndx(ipwc(1,ii),ipwc(2,ii),ipwc(3,ii))=ii
enddo

return
end subroutine mkinvpwndx

!*****************************************************************************

subroutine mkinvpw2ndx(ipwndx,napwndx,npwx,invpw2ndx)
implicit none
integer :: napwndx,npwx
integer :: ipwndx(2,napwndx)
integer :: invpw2ndx(npwx,npwx)
integer :: ipw1,ipw2,ii

do ipw1=1,npwx
do ipw2=1,npwx
  invpw2ndx(ipw1,ipw2)=0
enddo
enddo

do ii=1,napwndx
  invpw2ndx(ipwndx(1,ii),ipwndx(2,ii))=ii
enddo

return
end subroutine mkinvpw2ndx
