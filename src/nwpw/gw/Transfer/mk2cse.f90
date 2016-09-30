subroutine mk2cse(iband,ikk,lossfn, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& pwmatel,tpwmatel, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nqpt,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,bmet,blat,ipaw, &
& ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
& npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,iqndx,isymndx,isymndxq,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,zz,cse,xse)
implicit none
integer :: iband,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,npwtq,bantot,bantotq,ncg,ncgq,nkpt,nkptq,nqpt,nsym,nsymq,nsppol
integer :: npw,npwc,npwx,npwq,ipw1,npwndx,ntpwndx,napwndx,ipaw
double precision :: vol,pi,wmax,xred(3,natom)
double complex :: lossfn(nwpt,nqpt,ntpwndx)
double complex :: projwf(natom,nlmn,nkpt,ncband)
integer :: kg(3,npwt),kgq(3,npwtq)
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
double precision :: bmet(3,3),blat(3,3)
integer :: ipwx(3,npwx),ipwndx(2,napwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex, allocatable :: ssi(:)
double complex :: ssc(-nwpt:nwpt),cse,dcse,cse1(-nwpt:nwpt),zz,ctest,cse2(-nwpt:nwpt)
double precision :: xse,ssx,ssx1,ssx2,xse2
integer :: ii,jj,kk,ix,iy,iz,iskip,isign
integer :: iqq(3),iqqp(3),jka(3),jkb(3),jkk(3),iqv(3),ictr(3)
integer :: ikpt,ikptq,iks(3),ikslf1(3),ikslf2(3)
integer :: igg(3),igh(3),igg0(3),isym,isymq,ibp,iqpt,iqsym,iw,iww,ie,je,ies
double precision :: xck(3),xckq(3),qadj(6)
double complex :: cmatel,cmatel2,amatel,amatel2,vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3),nwpt),vfv(8,nwpt),vfx(8)
double precision :: vq2(8),vqvtx0(4),vqvtx(4),qkcvt(3,3)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,www,enval(8),dw,eshift,evtx0(4),evtx(4)
double precision :: rrpyr(3,4),vpyr0(4,nwpt),kvtx(3,4),xk(3,4),rg(3,3),tvol,vpyr(4),xpyr0(4),xpyr(4)
double precision :: de21,de31,de32,de41,de42,de43,thresh
double precision :: fbx(4),fb(4),cmx(3,4),cm(3,4),xkt(3,3)
double precision :: aa0(nwpt),av(3,nwpt),xme,xv(3)
double precision :: avec(3),bgrad(3),xmult
double precision :: abr,rlr
integer :: iwh,iwl,ibmin,ibmax,iwhi,iwlo,iehi,ielo
integer :: indxe(4),iwwhi,iwwlo
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),qq(3),qp(3),qq2,qp2,qs(3),ek,ekmq
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),temparray(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iipw,jjpw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwx,npwx)
integer :: pwsymndx(npwc,2*nsym)
integer :: ipw2,jw
double precision :: eps1,eps2,eps
double precision :: eval,brd,esprd
double precision :: gfo,gamma(3),pln(3),dist(3)
double precision :: sint1,sint1a,sint1b,svec(3),sveca(3),svecb(3)
double complex :: wint,wint0,w2int,w2int0,gterm,vcentr
double precision :: wcentr(-nwpt:nwpt),wgrid(nwpt)
double complex :: pwmatel(nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3))
logical :: lqsing,lqcentr
integer :: idum
logical :: lx,lc
double precision :: xdum,vdum(3)
double precision :: linterp
double precision :: rr(3,8)
integer :: ivndx(4,6),itet,iv
character*4 :: label
external linterp
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

abr=1.d-6
rlr=1.d-6
sei=0.d0
xse=(0.d0,0.d0)
ctest=(0.d0,0.d0)
cse=(0.d0,0.d0)
zz=(0.d0,0.d0)
igg0=(/0,0,0/)
do ie=-nwpt,nwpt
  cse1(ie)=(0.d0,0.d0)
enddo
xse=0.d0
dw=wmax/dble(nwpt)
do iw=1,nwpt
  wgrid(iw)=iw*dw
enddo
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
isym=isymndx(ikk(1),ikk(2),ikk(3))
ek=enrgy(indxkbnd(ikpt)+iband)
do ii=1,3
do jj=1,3
  qkcvt(ii,jj)=blat(ii,jj)/dble(ngkpt(ii))
enddo
enddo

do ix=1,ngkpt(1)
do iy=1,ngkpt(2)
do iz=1,ngkpt(3)
  temparray(ix,iy,iz)=(0.d0,0.d0)
enddo
enddo
enddo

gamma=(blat(1,:)+blat(2,:)+blat(3,:))/2.d0
do ii=1,3
  jj=mod(ii,3)+1
  kk=mod(ii+1,3)+1
  pln(1)=blat(jj,2)*blat(kk,3)-blat(jj,3)*blat(kk,2)
  pln(2)=-blat(jj,1)*blat(kk,3)+blat(jj,3)*blat(kk,1)
  pln(3)=blat(jj,1)*blat(kk,2)-blat(jj,2)*blat(kk,1)
  dist(ii)=dot_product(gamma,pln)/sqrt(dot_product(pln,pln))
enddo
gfo=minval(dist)

!ipaw=0
!write(6,'(a)') '      finding contibutions from:'
!write(6,'(a)') 'band  plane waves    correlation                exchange'
do ibp=nbcore+1,ncband
!do ibp=1,1
!do ibp=16,16
  seibc=sei
  cse2=cse1
  xse2=xse
  if (ibp.gt.nbocc) then
    isign=-1
  else
    isign=1
  endif
  do iipw=1,napwndx
!  do iipw=1,1
!  do iipw=15,15
!    write(6,'(a,i4)') 'iipw ',iipw
    seipw=sei
    do ie=-nwpt,nwpt
      ssc(ie)=(0.d0,0.d0)
    enddo
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
!    if (ipw1.ne.ipw2) cycle
!    write(6,'(38x,i3,14x,2i3)') ibp,ipw1,ipw2
    lx=ipw1.eq.ipw2.and.ibp.le.nbocc
    lc=iipw.le.ntpwndx
    if ((.not.lx).and.(.not.lc)) cycle

    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=2,2
!    do iy=2,2
!    do iz=4,4
      iqq=(/ix,iy,iz/)
!write(6,*) '>>>>> iqq = ',iqq
!write(6,'(a,3f10.5)') 'kpt = ',kpt(:,ikpt)
      iks=iqq-ngkpt/2
      jka=ikk-iks
      jkk=mod(jka,ngkpt(ii))
      do ii=1,3
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
      enddo
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
!      if (.true.) then
        ikptq=ikndxq(jkk(1),jkk(2),jkk(3))
        jkb=nint(kptq(:,ikptq)*dble(ngkpt))+ngkpt/2
        qq=kpt(:,ikpt)-kptq(:,ikptq)
!write(6,'(a,3f10.5)') 'kptq = ',kptq(:,ikptq)
      else
        ikptq=ikndx(jkk(1),jkk(2),jkk(3))
        jkb=nint(kpt(:,ikptq)*dble(ngkpt))+ngkpt/2
        qq=kpt(:,ikpt)-kpt(:,ikptq)
!write(6,'(a,3f10.5)') 'kptq = ',kptq(:,ikptq)
      endif
      igg=nint(dble(jkb-jka)/dble(ngkpt))
      qp=qq+igg+ipwx(:,ipw2)
      qq=qq+igg+ipwx(:,ipw1)
      qq2=0.d0
      qp2=0.d0
      do ii=1,3
      do jj=1,3
        qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
        qp2=qp2+qp(ii)*bmet(ii,jj)*qp(jj)
      enddo
      enddo
      vq(ix,iy,iz)=4.d0*pi/sqrt(qq2*qp2)
!write(6,'(a,3f10.5)') 'qq = ',qq
!write(6,'(a,3f10.5)') 'qp = ',qp
!write(6,'(a,3f10.5)') 'qq2 = ',qq2
!write(6,'(a,3f10.5)') 'qp2 = ',qp2
!write(6,*) 'vq = ',vq(ix,iy,iz)
!write(6,*) 'jka = ',jka
!write(6,*) 'jkk = ',jkk
!write(6,*) 'jkb = ',jkb
!write(6,*) 'igg = ',igg
!write(6,*) 'ikpt = ',ikpt
!write(6,*) 'ikptq = ',ikptq
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
        call mkmatelX(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
&         ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&         isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&         cmatel)
        if (ipaw.ne.0) then
!          call mkmatelP(pi,xred,natom,ibp,iband,ikpt,ikptq,qq,igg0,ngkpt, &
!&               pwmatel(:,:,ipw1,ix,iy,iz),tpwmatel(:,:,ipw1,ix,iy,iz), &
!&               projwf,nlmn,nkpt,nkptq,ncband,amatel)
        else
          amatel=(0.d0,0.d0)
        endif
        if (ipw1.eq.ipw2) then
          cmatel2=cmatel
          amatel2=amatel
        else
          call mkmatelX(ibp,iband,ikpt,ikptq,ipwx(:,ipw2),igh, &
&         ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&         isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&         cmatel2)
          if (ipaw.ne.0) then
!            call mkmatelP(pi,xred,natom,ibp,iband,ikpt,ikptq,qp,igg0,ngkpt, &
!&                 pwmatel(:,:,ipw2,ix,iy,iz),tpwmatel(:,:,ipw2,ix,iy,iz), &
!&                 projwf,nlmn,nkpt,nkptq,ncband,amatel2)
          else
            amatel2=(0.d0,0.d0)
          endif
        endif
        omega(ix,iy,iz)=enrgyq(indxkbndq(ikptq)+ibp)-ek+dw/2
      else
        call mkmatelX(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
&         ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&         isym,isymq,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&         cmatel)
        if (ipaw.ne.0) then
!          call mkmatelP(pi,xred,natom,ibp,iband,ikpt,ikptq,qq,igg0,ngkpt, &
!&               pwmatel(:,:,ipw1,ix,iy,iz),tpwmatel(:,:,ipw1,ix,iy,iz), &
!&               projwf,nlmn,nkpt,nkpt,ncband,amatel)
        else
          amatel=(0.d0,0.d0)
        endif
        if (ipw1.eq.ipw2) then
          cmatel2=cmatel
          amatel2=amatel
        else
          call mkmatelX(ibp,iband,ikpt,ikptq,ipwx(:,ipw2),igh, &
&         ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&         isym,isymq,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&         cmatel2)
          if (ipaw.ne.0) then
!            call mkmatelP(pi,xred,natom,ibp,iband,ikpt,ikptq,qp,igg0,ngkpt, &
!&                 pwmatel(:,:,ipw2,ix,iy,iz),tpwmatel(:,:,ipw2,ix,iy,iz), &
!&                 projwf,nlmn,nkpt,nkpt,ncband,amatel2)
          else
            amatel2=(0.d0,0.d0)
          endif
        endif
        omega(ix,iy,iz)=enrgy(indxkbnd(ikptq)+ibp)-ek+dw/2
      endif
      vqmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)*vq(ix,iy,iz)
!write(6,*) 'cmatel  = ',cmatel
!write(6,*) 'amatel  = ',amatel
!write(6,*) 'cmatel+amatel  = ',cmatel+amatel
!write(6,*) '|cmatel+amatel|^2  = ',dble((cmatel+amatel)*conjg(cmatel+amatel))
!write(6,*) 'cmatel2 = ',cmatel2
!write(6,*) 'vqmat2 = ',vqmat2(iqq(1),iqq(2),iqq(3))
!write(56,*) '>>>>>',iqq,vqmat2(iqq(1),iqq(2),iqq(3))
!      stvec(iqq(3))=vq(ix,iy,iz)
!      stvec(iqq(3))=sqrt(qq2)
!      stvec(iqq(3))=dble(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec2(iqq(3))=dimag(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec(iqq(3))=dble(cmatel)
!      stvec2(iqq(3))=dimag(cmatel)
!      stvec(iqq(3))=dble(cmatel2)
!      stvec2(iqq(3))=dimag(cmatel2)
!      stvec(iqq(3))=dble(cmatel*conjg(cmatel))
!      stvec(iqq(3))=dble(amatel*conjg(amatel))
!      stvec(iqq(3))=dble((cmatel+amatel)*conjg(cmatel+amatel))
!      temparray(ix,iy,iz)=temparray(ix,iy,iz)+stvec(iqq(3))
    enddo
!    write(76,'(10f8.2)') stvec(:)
!    write(76,'(10es8.1)') stvec(:)
!    write(76,'(10es8.1)') stvec2(:)
!    write(76,'(10f8.5)') 10.d0**(log10(stvec(:))-int(log10(stvec(:))))
!    write(76,*) 
    enddo
!    write(76,*) iqq(1)+1,'--------------------------------------------------',ibp
    enddo
!stop
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
      iqq=(/ix,iy,iz/)
      call fhilo(omega,iqq,ngkpt,whi,wlo)
      if (ix.eq.1.and.iy.eq.1.and.iz.eq.1) then
        wwhi=whi
        wwlo=wlo
      else
        wwhi=max(wwhi,whi)
        wwlo=min(wwlo,wlo)
      endif
    enddo
    enddo
    enddo

    iwlo=nint(wwlo*dble(nwpt)/wmax)
    iwhi=nint(wwhi*dble(nwpt)/wmax)
    if (ibp.gt.nbocc) then
      ielo=iwlo+1
      iehi=iwhi+nwpt
    else
      ielo=iwlo-nwpt
      iehi=iwhi-1
    endif
    allocate (ssi(ielo:iehi))
    do ie=ielo,iehi
      ssi(ie)=0.d0
    enddo
    ssx=0.d0

    if (lc) then
      do iw=1,nwpt
        do ix=1,ngkpt(1)
        do iy=1,ngkpt(2)
        do iz=1,ngkpt(3)
          iqq=(/ix,iy,iz/)
          iqpt=iqndx(iqq(1),iqq(2),iqq(3))
          call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
          if (jjpw.ne.0) then
            vfactor(iqq(1),iqq(2),iqq(3),iw)=vqmat2(iqq(1),iqq(2),iqq(3))*dimag(lossfn(iw,iqpt,jjpw))
          else
            vfactor(iqq(1),iqq(2),iqq(3),iw)=(0.d0,0.d0)
          endif
        enddo
        enddo
        enddo
      enddo
    endif

    if (ipw1.eq.1.or.ipw2.eq.1) then
      ix=ngkpt(1)/2
      iy=ngkpt(2)/2
      iz=ngkpt(3)/2
      qadj(1)=abs(vqmat2(ix,iy,iz)/vqmat2(ix+1,iy,iz))
      qadj(2)=abs(vqmat2(ix,iy,iz)/vqmat2(ix-1,iy,iz))
      qadj(3)=abs(vqmat2(ix,iy,iz)/vqmat2(ix,iy+1,iz))
      qadj(4)=abs(vqmat2(ix,iy,iz)/vqmat2(ix,iy-1,iz))
      qadj(5)=abs(vqmat2(ix,iy,iz)/vqmat2(ix,iy,iz+1))
      qadj(6)=abs(vqmat2(ix,iy,iz)/vqmat2(ix,iy,iz-1))
      if (qadj(1).gt.1.d3.and.qadj(2).gt.1.d3.and.qadj(3).gt.1.d3.and.  &
&         qadj(4).gt.1.d3.and.qadj(5).gt.1.d3.and.qadj(6).gt.1.d3) then
        lqsing=.true.
      else
        lqsing=.false.
      endif
    else
      lqsing=.false.
    endif

    ictr=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=6,6
!    do iy=5,5
!    do iz=5,5
!    do ix=1,1
!    do iy=1,1
!    do iz=1,1
      iqq=(/ix,iy,iz/)
      ssx1=ssx
      call fval(omega,iqq,iqqp,ngkpt,enval)
      if (lqsing) call fval(vq,iqq,iqqp,ngkpt,vq2)
      if (lc) then
        do iw=1,nwpt
          call fpol(vfactor(:,:,:,iw),ngkpt,iqq,iqqp,vfv(:,iw))
        enddo
      endif
      if (lx) call fpol(vqmat2(:,:,:),ngkpt,iqq,iqqp,vfx)
      do itet=1,6
!      do itet=1,1
        ssx2=ssx
        lqcentr=.false.
        do iv=1,4
          evtx0(iv)=enval(ivndx(iv,itet))
          rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
          if (lc) then
            do iw=1,nwpt
              vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
            enddo
          endif
          if (lx) xpyr0(iv)=vfx(ivndx(iv,itet))
          if (lqsing.and..not.lqcentr) then
            iqv=iqq+nint(rrpyr(:,iv))
            if (iqv(1).eq.ictr(1).and.iqv(2).eq.ictr(2).and.iqv(3).eq.ictr(3)) then
              lqcentr=.true.
            endif
          endif
        enddo
        call indxhpsort(4,4,evtx0,indxe)
        evtx=evtx0(indxe)    ! order vertices by energy
!write(6,'(4es15.5)') evtx
!        kvtx(1:3,1:4)=rrpyr(1:3,indxe)
        do ii=1,4
          jj=indxe(ii)
          do kk=1,3
            kvtx(kk,ii)=dot_product(iqq+rrpyr(:,jj)-ictr,qkcvt(:,kk))
          enddo
        enddo
        xk(1:3,1)=(/0.d0,0.d0,0.d0/)
        do ii=2,4
          xk(1:3,ii)=kvtx(1:3,ii)-kvtx(1:3,1)
        enddo
        call cross(xk(:,3),xk(:,4),vdum)
        tvol=dot_product(vdum,xk(:,2))
        do jj=1,3  ! make contragradient
          call cross(xk(1:3,mod(jj,3)+2),xk(1:3,mod(jj+1,3)+2),rg(1:3,jj))
        enddo
        rg=rg/tvol
        tvol=abs(tvol)
        iwwlo=nint(evtx(1)/dw)
        iwwhi=nint(evtx(4)/dw)
        de21=evtx(2)-evtx(1)
        de31=evtx(3)-evtx(1)
        de32=evtx(3)-evtx(2)
        de41=evtx(4)-evtx(1)
        de42=evtx(4)-evtx(2)
        de43=evtx(4)-evtx(3)
        thresh=1.d-5*de41
        if (lqcentr) then            ! integral around coulomb singularity
!        if (.true.) then
          do iv=1,4
            vqvtx0(iv)=vq2(ivndx(iv,itet))
          enddo
          bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
          do jj=1,3
            bgrad=bgrad+(evtx(jj+1)-evtx(1))*rg(:,jj)
          enddo
          xmult=4.d0*pi/sqrt(dot_product(bgrad,bgrad))
          if (lc) then
            do iw=1,nwpt
              vpyr(:)=vpyr0(:,iw)/vqvtx0
              aa0(iw)=vpyr(indxe(1))
              av(1:3,iw)=(/0.d0,0.d0,0.d0/)
              do jj=1,3
                av(1:3,iw)=av(1:3,iw) &
&                         +(vpyr(indxe(jj+1))-vpyr(indxe(1)))*rg(1:3,jj)
              enddo
            enddo
          endif
          if (lx) then
            xpyr=xpyr0/vqvtx0
            xme=xpyr(indxe(1))
            xv=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              xv=xv+(xpyr(indxe(jj+1))-xpyr(indxe(1)))*rg(:,jj)
            enddo
          endif
          do iww=iwwlo,iwwhi
            www=dw*iww
            if (www.lt.evtx(1)) then
              cycle
            elseif (www.lt.evtx(2)) then
              call singint(1,www,xk,kvtx(:,1),evtx,sint1,svec)
            elseif (www.lt.evtx(3)) then
              if (de21.lt.thresh.and.de43.lt.thresh) then
                call singint2(www,xk,kvtx(:,1),evtx,sint1b,svecb)
                sint1a=0.d0
                sveca=(/0.d0,0.d0,0.d0/)
              elseif (de21.ge.de43) then
                call singint(2,www,xk,kvtx(:,1),evtx,sint1a,sveca)
                call singint(1,www,xk,kvtx(:,1),evtx,sint1b,svecb)
              else
                call singint(3,www,xk,kvtx(:,1),evtx,sint1a,sveca)
                call singint(4,www,xk,kvtx(:,1),evtx,sint1b,svecb)
              endif
              sint1=sint1b-sint1a
              svec=svecb-sveca
            elseif (www.lt.evtx(4)) then
              call singint(4,www,xk,kvtx(:,1),evtx,sint1,svec)
            else
              cycle
            endif
            sint1=sint1*xmult
            svec=svec*xmult
            if (lc) then
              do iw=1,nwpt
                ie=iww-isign*iw
                ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
              enddo
            endif
            if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))*dw
          enddo
        else                         ! non-singular, tetrahedron method
          fbx(1)=tvol/(2.d0*de21*de31*de41)
          if (de21.ge.de43) then
            fbx(2)=tvol/(2.d0*de21*de32*de42)
          else
            fbx(3)=tvol/(2.d0*de31*de32*de43)
          endif
          fbx(4)=tvol/(2.d0*de41*de42*de43)
          do ii=1,4
            cmx(1:3,ii)=(/0.d0,0.d0,0.d0/)
            do jj=1,4
              if (jj.ne.ii) cmx(1:3,ii)=cmx(1:3,ii)+(xk(1:3,jj)-xk(1:3,ii)) &
&                                                  /(3*(evtx(jj)-evtx(ii)))
            enddo
          enddo
          if (lc) then
            do iw=1,nwpt
              aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
              av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
              do jj=1,3
                av(1:3,iw)=av(1:3,iw) &
&                         +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
              enddo
            enddo
          endif
          if (lx) then
            xme=xpyr0(indxe(1))
            xv=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              xv=xv+(xpyr0(indxe(jj+1))-xpyr0(indxe(1)))*rg(:,jj)
            enddo
          endif
          do iww=iwwlo,iwwhi
            www=dw*iww
            if (www.lt.evtx(1)) then
              cycle
            elseif (www.lt.evtx(2)) then
label='1  F'
              sint1=fbx(1)*(www-evtx(1))**2
              svec=(xk(1:3,1)+(www-evtx(1))*cmx(1:3,1))*sint1
            elseif (www.lt.evtx(3)) then
              if (de21.lt.thresh.and.de43.lt.thresh) then
                xkt(:,1)=xk(:,3)*(www-evtx(1))/de31
                xkt(:,2)=xk(:,4)*(www-evtx(1))/de41
                xkt(:,3)=(xk(:,4)-xk(:,2))*(www-evtx(2))/de42+xk(:,2)
                sint1=tvol*(2*(www-evtx(1))/(de31*de41) &
&                     -(www-evtx(1))**2*(de31+de41)/(de31*de41)**2)/2
                svec=((xkt(:,1)+xkt(:,3))/2.d0)*sint1
              elseif (de21.ge.de43) then
                fb(1)=fbx(1)*(www-evtx(1))**2
                fb(2)=fbx(2)*(www-evtx(2))**2
                sint1=fb(1)-fb(2)
                cm(1:3,1)=xk(1:3,1)+(www-evtx(1))*cmx(1:3,1)
                cm(1:3,2)=xk(1:3,2)+(www-evtx(2))*cmx(1:3,2)
                svec=cm(:,1)*fb(1)-cm(:,2)*fb(2)
              else
                fb(3)=fbx(3)*(www-evtx(3))**2
                fb(4)=fbx(4)*(www-evtx(4))**2
                sint1=fb(4)-fb(3)
                cm(1:3,3)=xk(1:3,3)+(www-evtx(3))*cmx(1:3,3)
                cm(1:3,4)=xk(1:3,4)+(www-evtx(4))*cmx(1:3,4)
                svec=cm(:,4)*fb(4)-cm(:,3)*fb(3)
              endif
            elseif (www.lt.evtx(4)) then
              sint1=fbx(4)*(www-evtx(4))**2
              svec=(xk(1:3,4)+(www-evtx(4))*cmx(1:3,4))*sint1
            else
              cycle
            endif
            if (lc) then
              do iw=1,nwpt
                ie=iww-isign*iw
                ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
              enddo
            endif
            if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))*dw
          enddo
        endif
!if (lqcentr) then
!write(47,'(3i3,3x,i3,f16.8,3x,a)') ix,iy,iz,itet,ssx-ssx2,'T'
!else
!write(47,'(3i3,3x,i3,f16.8,3x,a)') ix,iy,iz,itet,ssx-ssx2,'F'
!endif
      enddo
!      temparray(ix,iy,iz)=(ssx-ssx1)/((2*pi)**3)
    enddo
    enddo
    enddo
!!!!!!!!!!    ssi=ssi/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    if (lc) ssi=-isign*ssi/((2*pi)**3*pi)
    if (lx) ssx=-ssx/((2*pi)**3)

    if (lc) then
      do ie=-nwpt,nwpt
        eps1=dble(ie)*dw-dw/2
!        ssc(ie)=(0.d0,0.d0)
        do je=ielo,iehi
          if (ie.eq.je) then
            ssc(ie)=ssc(ie)+(0.d0,1.d0)*pi*ssi(je)
          else
            eps2=dble(je)*dw-dw/2
            ssc(ie)=ssc(ie)+isign*ssi(je)*dw/(eps2-eps1)
          endif
!write(6,'(2i4,2f10.3,es12.3,2(2x,2es11.3))') ie,je,eps1*27.2114,eps2*27.2114,1/(eps2-eps1),ssc(ie)*27.2114,dble(ssi(je)*27.2114)
        enddo
      enddo
    endif
    deallocate (ssi)

    if (lc) cse1=cse1+ssc
    if (lx) xse=xse+ssx
!    if (lx) then
!      write(6,'(i4,4x,i5,3x,3i3,5x,"(",f10.6,",",f10.6,")",5x,f10.6)')  &
!&          ibp,ipw1,ipwx(:,ipw1),((ssc(1)+ssc(0))/2)*27.2114,dble(ssx)*27.2114
!    else
!      write(6,'(i4,4x,i5,3x,3i3,5x,"(",f10.6,",",f10.6,")")') ibp,ipw1,ipwx(:,ipw1), &
!&          ((ssc(1)+ssc(0))/2)*27.2114
!    endif

  enddo
!  if (lx) then
!    write(6,'(i4,5x,"(",f10.6,",",f10.6,")",5x,f10.6)') ibp, &
!& (cse1(1)+cse1(0)-cse2(1)-cse2(0))*27.2114/2.d0, &
!& dble(xse-xse2)*27.2114
!  else
!    write(6,'(i4,5x,"(",f10.6,",",f10.6,")",f10.6)') ibp, &
!& (cse1(1)+cse1(0)-cse2(1)-cse2(0))*27.2114/2.d0
!  endif
enddo

do ie=-nwpt,nwpt
  eps1=dble(ie)*dw-dw/2
!  write(12,'(i4,f10.3,2(3x,2es12.3))') ie,eps1*27.2114,cse1(ie)*27.2114
enddo
cse=(cse1(0)+cse1(1))/2
dcse=(cse1(1)-cse1(0))/dw
zz=1.d0/(1.d0-dcse)
!write(6,*) cse*27.2114
!write(6,*) xse*27.2114
!write(6,*) zz
!write(6,*) dble(zz)*dimag(cse)*27.2114

!temparray=log10(temparray)
!temparray=10**(temparray-int(temparray))
!do ix=1,ngkpt(1)
!  write(76,*) ix,'--------------------------------------------------'
!  do iy=1,ngkpt(2)
!!    write(76,'(10es8.1)') dble(temparray(ix,iy,:))
!    write(76,'(10f8.5)') dble(temparray(ix,iy,:))
!  enddo
!enddo

return
end subroutine mk2cse

!***************************************************************************
!
! Over a surface S defined by the triangle with vertices xkvtx+xkp(1:3,1)
!                                                        xkvtx+xkp(1:3,2)
!                                                        xkvtx+xkp(1:3,3)
! where xkp defined to give S of constant energy 
! (see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972))
! and vector k measured from tip of xkvtx
! find      sint1 = integral dS 1/q^2
! and       svec  = integral dS k/q^2
! where q=k+xkvtx
! On S, define coordinates X and Y such that k=xkp(1:3,3)+Y*xq1+X*xq2
! normal=cross product(xq1,xq2)
! dot_product(k+xkvtx,k+xkvtx) = aa*Y^2+Y*(bb+cc*X)+dd*X+ee*X^2+FF
! integral dS -> |normal|*integral^1_0 dX integral^{1-X}_0 dY
! integral^{1-X}_0 dY 1/(dot_product(k+xkvtx,k+xkvtx)) = xfn1
! integral^{1-X}_0 dY X/(dot_product(k+xkvtx,k+xkvtx)) = X*xfn1 = xfn3
! integral^{1-X}_0 dY Y/(dot_product(k+xkvtx,k+xkvtx)) = (xfn2-bb*xfn1-cc*xfn3)/(2*aa)
! numerically integrate xfn1,xfn2,xfn3
! sint1 = integral^1_0 dX xfn1
! sint2 = integral^1_0 dX xfn2
! sint3 = integral^1_0 dX xfn3
! svec = xkp(:,3)*sint1+xq1*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2*sint3

subroutine singint(ive,ww,xk,xkvtx,evtx,sint1,svec)
implicit none
integer :: ive
double precision :: xk(3,4),xkvtx(3),evtx(4),ww
double precision :: xkp(3,3),xq1(3),xq2(3),xkv0(3),xcvt
double precision :: aa,bb,cc,dd,ee,ff,qq
double precision :: xmin,xmax,abr,rlr,xsing(20),error
double precision :: xfn1,xfn2,xfn3,grater
double precision :: sint1,sint2,sint3,svec(3),root1,root2
double precision :: xnorm(3),area
double precision :: xdum
integer :: ii,jj,ll
integer :: nsing,numcal,maxns,nroots
integer :: iter  !######
common /fn/ aa,bb,cc,dd,ee,ff,iter !######
external xfn1,xfn2,xfn3,grater

!write(6,*) ive
!write(6,'(4es15.5)') evtx
!write(6,'(3f10.6)') xk(:,1)
!write(6,*)
!write(6,*) 'boundary vectors of tetrahedron'
!write(6,'(3f10.6)') xk(:,2)
!write(6,'(3f10.6)') xk(:,3)
!write(6,'(3f10.6)') xk(:,4)
  if (ive.eq.1) then
    xkp(:,1)=(ww-evtx(1))*xk(:,2)/(evtx(2)-evtx(1))
    xkp(:,2)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
    xkp(:,3)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  elseif (ive.eq.2) then
    xkp(:,1)=(ww-evtx(1))*xk(:,2)/(evtx(2)-evtx(1))
    xkp(:,2)=(ww-evtx(2))*(xk(:,3)-xk(:,2))/(evtx(3)-evtx(2)) + xk(:,2)
    xkp(:,3)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
  elseif (ive.eq.3) then
    xkp(:,1)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
    xkp(:,2)=(ww-evtx(2))*(xk(:,3)-xk(:,2))/(evtx(3)-evtx(2)) + xk(:,2)
    xkp(:,3)=(ww-evtx(3))*(xk(:,4)-xk(:,3))/(evtx(4)-evtx(3)) + xk(:,3)
  else
    xkp(:,1)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
    xkp(:,2)=(ww-evtx(3))*(xk(:,4)-xk(:,3))/(evtx(4)-evtx(3)) + xk(:,3)
    xkp(:,3)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  endif
!write(6,*)
!write(6,*) 'vertices of S'
!write(6,'(3f10.6)') xkp(:,1)
!write(6,'(3f10.6)') xkp(:,2)
!write(6,'(3f10.6)') xkp(:,3)
  xq1=xkp(:,1)-xkp(:,3)
  xq2=xkp(:,2)-xkp(:,3)
!write(6,*)
!write(6,*) 'boundary vectors xq'
!write(6,'(3f10.6)') xq1
!write(6,'(3f10.6)') xq2
  xkv0=xkp(:,3)+xkvtx
!write(6,'(3f10.6)') xkvtx
!write(6,'(3f10.6)') xkv0
  call cross(xq1,xq2,xnorm)
  area=sqrt(dot_product(xnorm,xnorm))
!write(6,*)
!write(6,'(a,es15.5)') 'area  ',area
  aa=dot_product(xq1,xq1)
  bb=2*dot_product(xq1,xkv0)
  cc=2*dot_product(xq1,xq2)
  dd=2*dot_product(xq2,xkv0)
  ee=dot_product(xq2,xq2)
  ff=dot_product(xkv0,xkv0)
!write(6,'("coef",6f10.6)') aa,bb,cc,dd,ee,ff
  nsing=0
  call rquadroots(ee,dd,ff,nroots,root1,root2)
!write(6,'(i4,2es15.5)') nroots,root1,root2
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
  call rquadroots(1.d0+cc/2+ee,bb/2+dd,ff,nroots,root1,root2)
!write(6,'(i4,2es15.5)') nroots,root1,root2
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
!write(6,'(i4,4es15.5)') nsing,xsing(1:nsing)
  call hpsort(nsing,nsing,xsing(1:nsing))
!write(6,'("sings",i4,4es15.5)') nsing,xsing(1:nsing)
  xmin=0.d0
  xmax=1.d0
  abr=1.d-12
  rlr=1.d-6
  sint1=grater(xfn1,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  sint2=grater(xfn2,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  sint3=grater(xfn3,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
!write(6,*)
!write(6,'("sints",6f10.6)') sint1,sint2,sint3
  svec=xkp(:,3)*sint1+xq1*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2*sint3
!write(6,'("svec",6f10.6)') svec
!tint=(aa0*sint1+dot_product(av,svec))/bmag

end subroutine singint

!***************************************************************************

double precision function xfn1(xx)
implicit none
double precision :: xx,xx2
double precision :: aa,bb,cc,dd,ee,ff
double precision :: xt1,xt2,xl1,qq,dx
double precision :: sqq,xterm,df,d2f
integer :: iter  
common /fn/ aa,bb,cc,dd,ee,ff,iter 
!integer :: ict  
!data ict /0/  
!save ict  

  xx2=xx*xx
  xt1=bb+cc*xx
  xt2=dd*xx+ee*xx2+ff
  dx=2.d0*aa*(1.d0-xx)
  xl1=dx+xt1
  qq=4.d0*aa*xt2-xt1**2
  if (qq.lt.0.d0) then
    sqq=sqrt(-qq)
!    xfn1=(log(abs((sqq+xl1)/(sqq-xl1)))-log(abs((sqq+xt1)/(sqq-xt1))))/sqq
    xfn1=(log(abs(((sqq+xl1)*(sqq-xt1))/((sqq-xl1)*(sqq+xt1)))))/sqq
!write(iter,'(" B",es15.5,4x,2es15.5)') xx,xfn1
  else
!    df=2.d0*dx/(qq+xt1**2)
!    d2f=df**2*xt1/2.d0
!    if (d2f/df.lt.1.d-6) then
!      xfn1=df-d2f
!    else
      sqq=sqrt(qq)
      if (abs(dx/(qq+xt1**2)).gt.abs(1.d8*(dx**2*xt1/(qq+xt1**2)**2))) then
! atan(z+d)->atan(z)+d/(1+z^2)-d^2z/(1+z^2)^2+..., d->0 (taylor expansion)
! for small dx, xfn1->dx/(qq+xt1^2)
        xfn1=2*(dx/(qq+xt1**2)-dx**2*xt1/(qq+xt1**2)**2)
      elseif (abs(1.d0/xl1-1.d0/xt1).gt.abs(1.d8*qq*(1.d0/xl1**3-1.d0/xt1**3))) then
! atan(z)->pi/2-1/z+1/(3z^3)-1/(5z^5)+..., z->infty
! => at small sqq or large xt1, xl1, xfn1->-1/xl1+1/xt1
        xfn1=2*(-1.d0/xl1+1.d0/xt1+(qq/3.d0)*(1.d0/xl1**3-1.d0/xt1**3))
      else
        xfn1=2*(atan(xl1/sqq)-atan(xt1/sqq))/sqq
      endif
!write(iter,'(" A",es15.5,3x,4es15.5)') xx,xfn1,xt1,bb,cc
!    endif
  endif
!write(6,'(es15.5,4es15.5)') xx,qq,xl1,xt1,xfn1
!write(6,'(es15.5,4es15.5)') xx-1.d0,2/(qq+xt1**2),df,df-d2f,xfn1
!write(6,'(es15.5,4es15.5)') xx-1.d0,d2f/df,df,df-d2f,xfn1
!write(47,'(es15.7,4es15.7)') xx-1.d0,2/(qq+xt1**2),df,df-d2f,xfn1
!write(6,'(es15.7,4es15.7)') xx,2/(qq+xt1**2),qq,xt1**2,4.d0*aa*xt2
!write(6,'(es15.7,4es15.7)') xx,qq,xl1,xt1,xfn1
!write(47,'(es15.7,4es15.7)') xx,qq,xl1,xt1,xfn1
!write(6,'(es15.5,4es15.5)') xx,dx,qq,xt1**2,xfn1
!write(iter,'(es15.5,4x,2es15.5)') xx,xfn1
!ict=ict+1  
!if (ict.gt.3000) stop  

end function xfn1

!***************************************************************************

double precision function xfn2(xx)
implicit none
double precision :: xx,xx2,omx
double precision :: aa,bb,cc,dd,ee,ff
double precision :: xt1,xt2,xl1,qq
double precision :: xterm,sqq
integer :: iter  !######
common /fn/ aa,bb,cc,dd,ee,ff,iter !######

  xx2=xx*xx
  omx=1.d0-xx
  xterm=(aa*omx*omx+omx*(bb+cc*xx)+dd*xx+ee*xx2+ff)/(dd*xx+ee*xx2+ff)
  xfn2=log(abs(xterm))
!  xt1=bb+cc*xx
!  xt2=dd*xx+ee*xx2+ff
!  xl1=2.d0*aa*omx+xt1
!  qq=4.d0*aa*xt2-xt1**2
!  if (qq.lt.0.d0) then
!    sqq=sqrt(-qq)
!    xfn2=log(abs(xterm))+xt1*(log(abs((sqq+xl1)/(sqq-xl1)))-log(abs((sqq+xt1)/(sqq-xt1))))/sqq
!  else
!    sqq=sqrt(qq)
!    xfn2=log(abs(xterm))+2*xt1*(atan((xl1)/sqq)-atan(xt1/sqq))/sqq
!  endif
!!write(47,'(4es15.5)') xx,xfn2,(aa+cc+ee)*xx2+(bb+dd)*xx+ff,dd*xx+ee*xx2+ff

end function xfn2

!***************************************************************************

double precision function xfn3(xx)
implicit none
double precision :: xfn1,xx
external xfn1

  xfn3=xx*xfn1(xx)

end function xfn3

!***************************************************************************
!
! Over a surface S defined by the parallelogram with vertices 
!           xkvtx+xkp(1:3,1)
!           xkvtx+xkp(1:3,2)
!           xkvtx+xkp(1:3,3)
!           xkvtx+xkp(1:3,1)+(xkp(1:3,2)-xkp(1:3,3))
! where xkp defined to give S of constant energy 
! (see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972))
! and vector k measured from tip of xkvtx
! find      sint1 = integral dS 1/q^2
! and       svec  = integral dS k/q^2
! where q=k+xkvtx
! On S, define coordinates X and Y such that k=xkp(1:3,3)+Y*xq1+X*xq2
! normal=cross product(xq1,xq2)
! dot_product(k+xkvtx,k+xkvtx) = aa*Y^2+Y*(bb+cc*X)+dd*X+ee*X^2+FF
! integral dS -> |normal|*integral^1_0 dX integral^1_0 dY
! integral^1_0 dY 1/(dot_product(k+xkvtx,k+xkvtx)) = xfn4
! integral^1_0 dY X/(dot_product(k+xkvtx,k+xkvtx)) = X*xfn4 = xfn6
! integral^1_0 dY Y/(dot_product(k+xkvtx,k+xkvtx)) = (xfn5-bb*xfn4-cc*xfn6)/(2*aa)
! numerically integrate xfn4,xfn5,xfn6
! sint1 = integral^1_0 dX xfn4
! sint2 = integral^1_0 dX xfn5
! sint3 = integral^1_0 dX xfn6
! svec = xkp(:,3)*sint1+xq1*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2*sint3

subroutine singint2(ww,xk,xkvtx,evtx,sint1,svec)
implicit none
double precision :: xk(3,4),xkvtx(3),evtx(4),ww
double precision :: xkp(3,3),xq1(3),xq2(3),xkv0(3),xcvt
double precision :: aa,bb,cc,dd,ee,ff,qq
double precision :: xmin,xmax,abr,rlr,xsing(20),error
double precision :: xfn4,xfn5,xfn6,grater
double precision :: sint1,sint2,sint3,svec(3),root1,root2
double precision :: xnorm(3),area
double precision :: xdum
integer :: ii,jj,ll
integer :: nsing,numcal,maxns,nroots
integer :: iter  !######
common /fn/ aa,bb,cc,dd,ee,ff,iter !######
external xfn4,xfn5,xfn6,grater

!write(6,*) ive
!write(6,'(4es15.5)') evtx
!write(6,'(3f10.6)') xk(:,1)
!write(6,*)
!write(6,*) 'boundary vectors of tetrahedron'
!write(6,'(3f10.6)') xk(:,2)
!write(6,'(3f10.6)') xk(:,3)
!write(6,'(3f10.6)') xk(:,4)
  xkp(:,1)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
  xkp(:,2)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  xkp(:,3)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
!write(6,*)
!write(6,*) 'vertices of S'
!write(6,'(3f10.6)') xkp(:,1)
!write(6,'(3f10.6)') xkp(:,2)
!write(6,'(3f10.6)') xkp(:,3)
  xq1=xkp(:,1)-xkp(:,3)
  xq2=xkp(:,2)-xkp(:,3)
!write(6,*)
!write(6,*) 'boundary vectors xq'
!write(6,'(3f10.6)') xq1
!write(6,'(3f10.6)') xq2
  xkv0=xkp(:,3)+xkvtx
!write(6,'(3f10.6)') xkvtx
!write(6,'(3f10.6)') xkv0
  call cross(xq1,xq2,xnorm)
  area=sqrt(dot_product(xnorm,xnorm))
!write(6,*)
!write(6,'(a,es15.5)') 'area  ',area
  aa=dot_product(xq1,xq1)
  bb=2*dot_product(xq1,xkv0)
  cc=2*dot_product(xq1,xq2)
  dd=2*dot_product(xq2,xkv0)
  ee=dot_product(xq2,xq2)
  ff=dot_product(xkv0,xkv0)
!write(6,'("coef",6f10.6)') aa,bb,cc,dd,ee,ff
  nsing=0
  call rquadroots(ee,dd,ff,nroots,root1,root2)
!write(6,'(i4,2es15.5)') nroots,root1,root2
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
  call rquadroots(1.d0+cc/2+ee,bb/2+dd,ff,nroots,root1,root2)
!write(6,'(i4,2es15.5)') nroots,root1,root2
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
!write(6,'(i4,4es15.5)') nsing,xsing(1:nsing)
  call hpsort(nsing,nsing,xsing(1:nsing))
!write(6,'("sings",i4,4es15.5)') nsing,xsing(1:nsing)
  xmin=0.d0
  xmax=1.d0
  abr=1.d-12
  rlr=1.d-6
  sint1=grater(xfn4,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  sint2=grater(xfn5,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  sint3=grater(xfn6,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
!write(6,*)
!write(6,'("sints",6f10.6)') sint1,sint2,sint3
  svec=xkp(:,3)*sint1+xq1*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2*sint3
!write(6,'("svec",6f10.6)') svec
!tint=(aa0*sint1+dot_product(av,svec))/bmag

end subroutine singint2

!***************************************************************************

double precision function xfn4(xx)
implicit none
double precision :: xx,xx2
double precision :: aa,bb,cc,dd,ee,ff
double precision :: xt1,xt2,xl1,qq
double precision :: sqq,xterm
integer :: iter  !######
common /fn/ aa,bb,cc,dd,ee,ff,iter !######

  xx2=xx*xx
  xt1=bb+cc*xx
  xt2=dd*xx+ee*xx2+ff
  xl1=2.d0*aa+xt1
  qq=4.d0*aa*xt2-xt1**2
  if (qq.lt.0.d0) then
    sqq=sqrt(-qq)
    xfn4=(log(abs((sqq+xl1)/(sqq-xl1)))-log(abs((sqq+xt1)/(sqq-xt1))))/sqq
  else
    sqq=sqrt(qq)
    xfn4=2*(atan((xl1)/sqq)-atan(xt1/sqq))/sqq
  endif
!write(47,'(es15.5,4x,2es15.5)') xx,xfn4

end function xfn4

!***************************************************************************

double precision function xfn5(xx)
implicit none
double precision :: xx,xx2
double precision :: aa,bb,cc,dd,ee,ff
double precision :: xt1,xt2,xl1,qq
double precision :: xterm,sqq
integer :: iter  !######
common /fn/ aa,bb,cc,dd,ee,ff,iter !######

  xx2=xx*xx
  xterm=(aa+bb+(cc+dd)*xx+ee*xx2+ff)/(dd*xx+ee*xx2+ff)
  xfn5=log(abs(xterm))
!  xt1=bb+cc*xx
!  xt2=dd*xx+ee*xx2+ff
!  xl1=2.d0*aa+xt1
!  qq=4.d0*aa*xt2-xt1**2
!  if (qq.lt.0.d0) then
!    sqq=sqrt(-qq)
!    xfn5=log(abs(xterm))+xt1*(log(abs((sqq+xl1)/(sqq-xl1)))-log(abs((sqq+xt1)/(sqq-xt1))))/sqq
!  else
!    sqq=sqrt(qq)
!    xfn5=log(abs(xterm))+2*xt1*(atan((xl1)/sqq)-atan(xt1/sqq))/sqq
!  endif

end function xfn5

!***************************************************************************

double precision function xfn6(xx)
implicit none
double precision :: xfn4,xx
external xfn4

  xfn6=xx*xfn4(xx)

end function xfn6

