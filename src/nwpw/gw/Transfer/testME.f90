subroutine testME(igglf1,iks,ikk,iocc,iunocc, &
& pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
& nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
& indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
& ihlf,lvtrans,ntypepaw, &
& pwjmatel,tpwjmatel,pwmatel1,tpwmatel1, &
& igndx,ikndx,isymndx,isymg)
! test symmetries of matrix elements
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),dq(3),pi,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: igglf1(3)
double precision :: shiftk(3)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
double complex :: cmatel
! cmatel - Complex density MATrix ELement
double complex :: amatel
! amatel - contribution to matrix elements from PAW
double complex :: jmatel(3)
! jmatel - current matrix element
double complex :: vmatel(3)
! vmatel - contribution to current matrix elements from PAW
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
integer :: iwh,iwl
double precision :: testksym(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwjmatel(ntypepaw,nlmn,nlmn),tpwjmatel(ntypepaw,nlmn,nlmn)
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double precision :: xck(3),xqq(3),xcksym(3),xqqsym(3)
integer :: ikk0(3),ikksym0(3),ikc(3),ikksym(3),ikssym(3),iggsym(3)
double complex :: jsym(3)

write(22,*) "Begin symmetry tests"
write(22,*) "input"
write(22,*) "bands: ",iocc,iunocc
write(22,'(a,3i3)') "q: ",iks
write(22,'(a,3i3)') "k: ",ikk
write(22,'(a,3i3)') "G: ",igglf1
!write(22,'(3i3,3x,3i3,3x,3i3)') ikk,iks,igglf1
ikk0=ikk
do ii=1,3
  ikk(ii)=mod(ikk(ii)-1,ngkpt(ii))+1
  if (ikk(ii).le.0) ikk(ii)=ikk(ii)+ngkpt(ii)
  xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
  xck(ii)=mod(xck(ii),1.d0)
  if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
  xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
  xck(ii)=xck(ii)-0.5d0
  xqq(ii)=dble(iks(ii))/dble(ngkpt(ii))
  igglf1(ii)=igglf1(ii)+(ikk0(ii)-ikk(ii))/ngkpt(ii)
enddo
write(22,*)
write(22,*) "fixed input"
write(22,'(a,3i3,3x,3f10.6)') "q: ",iks,xqq
write(22,'(a,3i3,3x,3f10.6)') "k: ",ikk,xck
write(22,'(a,3i3)') "G: ",igglf1
do isym=1,9,8
  xcksym=(/0.d0,0.d0,0.d0/)
  xqqsym=(/0.d0,0.d0,0.d0/)
  iggsym=(/0,0,0/)
  do ii=1,3
  do jj=1,3
    xcksym(jj)=xcksym(jj)+xck(ii)*symrel(ii,jj,isym)
    xqqsym(jj)=xqqsym(jj)+xqq(ii)*symrel(ii,jj,isym)
    iggsym(jj)=iggsym(jj)+igglf1(ii)*symrel(ii,jj,isym)
  enddo
  enddo
  do ii=1,3
    ikksym(ii)=nint((xcksym(ii)+0.5d0)*ngkpt(ii))
    ikksym0(ii)=ikksym(ii)
    ikksym(ii)=mod(ikksym(ii)-1,ngkpt(ii))+1
    if (ikksym(ii).le.0) ikksym(ii)=ikksym(ii)+ngkpt(ii)
    ikssym(ii)=nint(xqqsym(ii)*ngkpt(ii))
    iggsym(ii)=iggsym(ii)+(ikksym0(ii)-ikksym(ii))/ngkpt(ii)
  enddo
write(22,*)
write(22,*) "Symmetry: ",isym
do ii=1,3
write(22,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
enddo
write(22,'(a,3i3,3x,3f10.6)') "q: ",ikssym,xqqsym
write(22,'(a,3i3,3x,3f10.6)') "k: ",ikksym,xcksym
write(22,'(a,3i3)') "G: ",iggsym
  call testMErho(iggsym,ikssym,ikksym,iocc,iunocc, &
&   pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
&   nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
&   indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
&   ihlf,lvtrans,ntypepaw,pwmatel1,tpwmatel1, &
&   igndx,ikndx,isymndx,isymg,cmatel)
  call testMEj(iggsym,ikssym,ikksym,iocc,iunocc, &
&   pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
&   nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
&   indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
&   ihlf,lvtrans,ntypepaw, &
&   pwjmatel,tpwjmatel,pwmatel1,tpwmatel1, &
&   igndx,ikndx,isymndx,isymg,jmatel)
  jsym=(/0.d0,0.d0,0.d0/)
  do ii=1,3
  do jj=1,3
    jsym(jj)=jsym(jj)+jmatel(ii)*syminv(ii,jj,isym)
  enddo
  enddo
!  write(6,'(i2,1x,4("(",f8.5,",",f8.5,")"))') isym,cmatel,jsym
  write(22,'(i2,1x,4("(",f8.5,",",f8.5,")"))') isym,cmatel,jsym
!  write(6,'(i2,1x,4("(",f8.5,",",f8.5,")"))') isym,cmatel,jmatel
!  write(6,'(i2,1x,4("(",es8.1,",",es8.1,")"))') isym,cmatel,jmatel
!  write(22,'(3i3,1x,4("(",es12.5,",",es12.5,")"))') ikk,cmatel
  do ii=1,3
    jkk(ii)=mod(ikksym(ii)-ikssym(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
    if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
  enddo
  ikptq=ikndx(jkk(1),jkk(2),jkk(3))
  do ig=1,npwarr(ikpt)
    write(60+isym,'(i4,2x,3i3,2x,"(",f8.5,",",f8.5,")")') ig, &
&             kg(1:3,indxkpw(ikptq)+ig), &
&             cg(indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+ig)
  enddo
enddo

return
end subroutine testME

!*************************************************************************

subroutine testMErho(igglf1,iks,ikk,iocc,iunocc, &
& pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
& nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
& indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
& ihlf,lvtrans,ntypepaw,pwmatel1,tpwmatel1, &
& igndx,ikndx,isymndx,isymg,cmatel)
! test to see if density matrix elements are unaffected by symmetry of lattice
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),dq(3),pi,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: igglf1(3)
double precision :: shiftk(3)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
double complex :: cmatel
! cmatel - Complex density MATrix ELement
double complex :: amatel
! amatel - contribution to matrix elements from PAW
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: testksym(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)

igg0=(/0,0,0/)

!write(6,'(a)') 'occupied band, unoccupied band'
!    write(6,*) iocc,iunocc
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)  ! integer index of k-q
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
      write(22,'(4(3i3,3x))') ikk,jkk,iks,igg
      write(22,'(5(3i3,3x))') ikk,iks,jka,jkk,igg
      write(22,*) ikpt,isym,ikptq,isymq
      write(22,'(3f10.6,4x,3f10.6)') kpt(:,ikpt),kpt(:,ikptq)
!      write(6,'(3(3f10.6,3x))') xck,xckq
!      write(6,'(3(3f10.6,3x))') qq,qq-igglf1
!write(44,*) "band (hole) = ",iunocc
!write(44,*) "band (conduction electron) = ",iocc
!write(44,'(a,3f8.4)') "k   = ",xck
!write(44,'(a,3f8.4)') "k-q = ",xckq
!write(44,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(44,'(a,3f8.4)') "k-q in irreducable Brillouin zone = ",kpt(1:3,ikptq)
!write(44,'(a,3i3,3x)') "integer index of k   = ",ikk
!write(44,'(a,3i3,3x)') "integer index of k-q = ",jkk
!write(44,'(a,3i3,3x)') "integer index of q   = ",iks
!write(44,*) "symmetry operation on k"
!write(44,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!testksym=(/0.d0,0.d0,0.d0/)
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + syminv(jj,ii,mod(isym-1,nsym)+1)*(xck(jj)-igsymk(jj))
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!testksym=igsymk
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + symrel(jj,ii,mod(isym-1,nsym)+1)*kpt(jj,ikpt)
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!write(44,*)
!write(44,*) "symmetry operation on k-q"
!write(44,*) "symmetry index ",isymq,mod(isymq-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!enddo
!testksym=(/0.d0,0.d0,0.d0/)
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + syminv(jj,ii,mod(isymq-1,nsym)+1)*(xckq(jj)-igsymq(jj))
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!testksym=igsymq
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + symrel(jj,ii,mod(isymq-1,nsym)+1)*kpt(jj,ikptq)
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!write(44,*)
!      write(6,'(a,3i5)') "ikpt  = ",ikpt
!      write(6,'(a,3i5)') "ikptq = ",ikptq
!      write(6,'(a,3i3,3x)') "integer index of k   = ",ikk
!      write(6,'(a,3i3,3x)') "integer index of k-q = ",jkk
!      write(6,*) "k-points used"
!      write(6,'(a,3f8.4)') "k   = ",xck
!      write(6,'(a,3f8.4)') "k-q = ",xckq
!      write(6,*) "k-points symmetry transform to irreducable Brillouin zone"
!      write(6,'(a,3f8.4)') "k   = ",kpt(1:3,ikpt)
!      write(6,'(a,3f8.4)') "k-q = ",kpt(1:3,ikptq)
!!      write(6,'(3f8.4,4x,3f8.4)') kpt(1:3,ikpt),kpt(1:3,ikptq)
!!      write(6,'(3i5)') isym,isymq,nsym
!      write(6,'(a,3i5)') "isym  = ",isym
!      write(6,'(a,3i5)') "isymq = ",isymq
!      write(6,'(2i4,4x,2i4)') mod(isym-1,nsym)+1,-(2*((isym-1)/nsym)-1),mod(isymq-1,nsym)+1,-(2*((isymq-1)/nsym)-1)
!      write(6,*)
!      write(6,*) "symmetry: k"
!      testksym=(/0.d0,0.d0,0.d0/)
!      do ii=1,3
!        write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!        do jj=1,3
!          testksym(ii) = testksym(ii) + symrel(ii,jj,mod(isym-1,nsym)+1)*xck(jj)
!        enddo
!      enddo
!      write(6,'(3f8.4)') testksym
!      write(6,*)
!      write(6,*) "symmetry: k-q"
!      testksym=(/0.d0,0.d0,0.d0/)
!      do ii=1,3
!        write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!        do jj=1,3
!          testksym(ii) = testksym(ii) + symrel(ii,jj,mod(isymq-1,nsym)+1)*xckq(jj)
!        enddo
!      enddo
!      write(6,'(3f8.4)') testksym
!      do ii=1,npwarr(ikptq)
!        write(91,'(i5,i10,2f12.6)') ii,indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+ii,cg(indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+ii)
!      enddo
      call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&       ncg,nkpt,npwt,igmx,igmn,igndx, &
&       isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&       lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&       cg,indxkcg,indxkpw,npwarr,kg, &
&       cmatel)
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      if (ipaw.ne.0) then
        call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
& pwmatel1,tpwmatel1,nbcore,ncband,amatel)
      else
        amatel=(0.d0,0.d0)
      endif
!write(6,*) jkk
!write(6,*) ihlf(ikpt),ihlf(ikptq)
!write(6,*) cmatel
!write(6,*) amatel
!write(6,*) cmatel+amatel
!write(6,*) (cmatel+amatel)*conjg(cmatel+amatel)
!write(6,*) 
!write(6,*) cmat2(ikk(1),ikk(2),ikk(3)),omega(ikk(1),ikk(2),ikk(3))
!write(46,'(2i3,2x,3i3,2x,f7.4,1x,"(",e10.3,",",e10.3,")",2x,2i3)') &
!&  iocc,iunocc,ikk(1:3),omega(ikk(1),ikk(2),ikk(3)), &
!&  cmat2(ikk(1),ikk(2),ikk(3)) !, isym, isymq
!write(6,*) dble(cmat2(ikk(1),ikk(2),ikk(3)))
!write(6,*) 
!      stvec(ikk(3))=enrgy(indxkbnd(ikpt)+iunocc)
!      stvec(ikk(3))=omega(ikk(1),ikk(2),ikk(3))-dw*58
!      stvec(ikk(3))=dble(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec2(ikk(3))=dimag(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec(ikk(3))=abs(cmatel)
!      stvec(ikk(3))=dble(cmatel)
!      stvec2(ikk(3))=dimag(cmatel)
!      stvec(ikk(3))=dble(cmatellf)
!      stvec2(ikk(3))=dimag(cmatellf)
      cmatel = cmatel+amatel

return
end subroutine testMErho

!*************************************************************************

subroutine testMEj(igglf1,iks,ikk,iocc,iunocc, &
& pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
& nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
& indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
& ihlf,lvtrans,ntypepaw, &
& pwjmatel,tpwjmatel,pwmatel1,tpwmatel1, &
& igndx,ikndx,isymndx,isymg,jmatel)
! test to see if current matrix elements transform as a vector under lattice symmetries
use geometry
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),qq2,qp2,dq(3),pi,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
! iks: integer index of q-vector in Brillouin zone
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: igglf1(3)
double precision :: shiftk(3)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
double complex :: jmatel(3)
! jmatel - current matrix element
double complex :: vmatel(3)
! vmatel - contribution to current matrix elements from PAW
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,ll,mm,nn,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq,iiq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwjmatel(ntypepaw,nlmn,nlmn),tpwjmatel(ntypepaw,nlmn,nlmn)
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double precision :: matdummy1(3,3),matdummy2(3,3),matdummy3(3,3),dummy
double complex :: cmatdummy1(3,3),cmatdummy2(3,3),cmatdummy3(3,3)

igg0=(/0,0,0/)

      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)  ! integer index of k-q
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
!      write(6,'(4(3i3,3x))') ikk,jkk,iks,igg
!      write(6,'(3(3f10.6,3x))') xck,xckq
!      write(6,'(3(3f10.6,3x))') qq,qq-igglf1
!write(6,'(a,3f8.4)') "k   = ",xck
!write(6,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(6,*) "symmetry operation on k"
!write(6,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!write(44,*) "band = ",iunocc
!write(44,'(a,3f8.4)') "k   = ",xck
!write(44,'(a,3f8.4)') "k-q = ",xckq
!write(44,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(44,'(a,3f8.4)') "k-q in irreducable Brillouin zone = ",kpt(1:3,ikptq)
!write(44,'(a,3i3,3x)') "integer index of k   = ",ikk
!write(44,'(a,3i3,3x)') "integer index of k-q = ",jkk
!write(44,'(a,3i3,3x)') "integer index of q   = ",iks
!write(44,*) "symmetry operation on k"
!write(44,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!write(44,*)
!write(44,*) "symmetry operation on k-q"
!write(44,*) "symmetry index ",isymq,mod(isymq-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!enddo
!write(44,*)
!write(44,'(a,3i3)') "k   folding vector: ",igsymk
!write(44,'(a,3i3)') "k-q folding vector: ",igsymq
!write(44,'(a,3i3)') "igg               : ",igg
!write(44,'(a,3i3)') "igglf             : ",igglf1
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&       ncg,nkpt,npwt,igmx,igmn,igndx, &
&       isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&       lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&       cg,indxkcg,indxkpw,npwarr,kg, &
&       jmatel)
!write(6,'(3("(",f8.4,",",f8.4,")"))') jmatel
      if (ipaw.ne.0) then
        call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&                       pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
        jmatel(:) = jmatel(:) + vmatel(:)
      endif
!write(6,'(3("(",f8.4,",",f8.4,")"))') vmatel
!write(6,'(3("(",f8.4,",",f8.4,")"))') jmatel

!stop

return
end subroutine testMEj

