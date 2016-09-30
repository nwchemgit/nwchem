subroutine mkrho(igglf,iks,vol,pi, &
& nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,ipaw, &
& pwmatel1,tpwmatel1,pwmatel2,tpwmatel2, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,isymndx,isymndxq,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,pola,polb)
implicit none
double precision :: qq(3),qp(3),dq(3),vq,vol,pi,wmax,dw,xred(3,natom)
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npwtq,npw,npwq,iks(3)
integer :: bantot,bantotq,ncg,ncgq,nkpt,nkptq,ipw1,ipw2,ipw,nlmn,natom
integer :: nsym,nsymq,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: kg(3,npwt),kgq(3,npwtq),igglf1(3),igglf2(3)
double complex :: projwf(natom,nlmn,nkpt,ncband)
double precision :: enrgy(bantot),enrgyq(bantotq)
double complex :: cg(ncg),cgq(ncgq)
integer :: indxkpw(nkpt),indxkpwq(nkptq),indxkbnd(nkpt),indxkbndq(nkptq)
integer :: indxkcg(nkpt),indxkcgq(nkptq),npwarr(nkpt),npwarrq(nkptq)
double precision :: kpt(3,nkpt),kptq(3,nkptq),shiftk(3),shiftkq(3)
integer :: nsymk(nkpt),nsymkq(nkptq),symk(nkpt,nsym*2),symkq(nkptq,nsymq*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsymq)
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: lvtransq(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt),ihlfq(nkptq)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: pola(nwpt),polb(nwpt),dielf(npwt),tpolincr,polci(npwt),dpol
! pola - hermitian part of vq * polarization matrix element for igglf1,igglf2
! polb - antihermitian part of vq * polarization matrix element 
double complex :: cmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,cmat2v(8)
double complex :: amatel,amatellf
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double complex :: pwmatel1(nlmn,nlmn),pwmatel2(nlmn,nlmn),tpwmatel1(nlmn,nlmn),tpwmatel2(nlmn,nlmn)

  igg0=(/0,0,0/)
  rho=(0.d0,0.d0)

  do iocc=nbcore+1,nbocc 
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=1,1
!    do iz=1,1
      ikk=(/ix,iy,iz/)
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftkq(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndxq(jkk(1),jkk(2),jkk(3))
      isymq=isymndxq(jkk(1),jkk(2),jkk(3))
!      write(6,'(4(3i3,3x))') ikk,jkk,iks,igg
!      write(6,'(3(3f10.6,3x))') xck,xckq
!      write(6,'(3(3f10.6,3x))') qq,qq-igglf1
!      write(6,'(3f8.4,4x,3f8.4)') kpt(1:3,ikpt),kptq(1:3,ikptq)
!      write(6,'(3i5)') ikpt,ikptq
!      write(6,'(3i5)') isym,isymq,nsym
!      write(6,'(2i4,4x,2i4)') mod(isym-1,nsym)+1,-(2*((isym-1)/nsym)-1),mod(isymq-1,nsym)+1,-(2*((isymq-1)/nsym)-1)
      call mkmatelX(iocc,iocc,ikpt,ikpt,igglf1,igg, &
&       ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&       isym,isymq,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&       lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&       cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&       cmatel)
      if (ipaw.ne.0) then
!        call mkmatelP(pi,xred,natom,iocc,iunocc,ikpt,ikptq,qq,igg0,ngkpt,pwmatel1,tpwmatel1,projwf,nlmn,nkpt,nkptq,ncband,amatel)
      else
        amatel=(0.d0,0.d0)
      endif
      rho=rho+cmatel+amatel
    enddo
    enddo
    enddo
  enddo
  rho=rho/(vol*ngkpt(1)*ngkpt(2)*ngkpt(3))

return
end subroutine mkrho
