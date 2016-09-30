subroutine mk3cse(iband,ikk,lossfn, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& test_bands_se, &
& pwmatel,tpwmatel, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nqpt,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,bmet,blat,ipaw,itetrahedron, &
& ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
& npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,iqndx,isymndx,isymndxq,isymg,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,zz,cse,xse)
implicit none
integer :: iband,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,npwtq,bantot,bantotq,ncg,ncgq,nkpt,nkptq,nqpt,nsym,nsymq,nsppol
integer :: npw,npwc,npwx,npwq,ipw1,npwndx,ntpwndx,napwndx,ipaw,itetrahedron
double precision :: vol,pi,wmax,xred(3,natom)
double complex :: lossfn(nwpt,nqpt+9,ntpwndx)
integer :: test_bands_se(2)
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
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex, allocatable :: ssi(:)
double complex :: ssc(-nwpt:nwpt),cse,dcse,cse0,dcse0,cse1(-nwpt:nwpt),zz,ctest,cse2(-nwpt:nwpt)
double precision :: xse,xse0,ssx,ssx1,ssx2,xse2,ssx3
integer :: ii,jj,kk,ix,iy,iz,iskip,ocsign
integer :: iqq(3),iqqp(3),jka(3),jkb(3),jkk(3),iqv(3),ictr(3)
integer :: ikpt,ikptq,iks(3),ikslf1(3),ikslf2(3)
integer :: igg(3),igh(3),igg0(3),isym,isymq,ibp,iqpt,iqsym,iw,iww,ie,je,ies,iqctr
double precision :: xck(3),xckq(3),qadj(6)
double complex :: cmatel,cmatel2,amatel,amatel2
double complex :: vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: xmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3),nwpt),vfv(8,nwpt),vfx(8)
double precision :: vq2(8),vqvtx0(4),vqvtx(4),qkcvt(3,3)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,www,wwwlo,wwwhi,enval(8),dw,eshift,evtx0(4),evtx(4)
double precision :: rrpyr(3,4),kvtx(3,4),xk(3,4),rg(3,3),tvol,xpyr0(4),xpyr(4)
double complex :: vpyr0(4,nwpt),vpyr(4)
double precision :: de21,de31,de32,de41,de42,de43,thresh
double precision :: fbx(4),fb(4),cmx(3,4),cm(3,4),xkt(3,3)
double complex :: aa0(nwpt),av(3,nwpt)
double precision :: xme,xv(3)
double precision :: avec(3),bgrad(3),xmult
double precision :: abr,rlr
integer :: iwh,iwl,ibmin,ibmax,iwhi,iwlo,iehi,ielo
integer :: indxe(4),iwwlim(4)
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),qq(3),qp(3),qq2,qp2,qs(3),ek,ekmq
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),temparray(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iipw,jjpw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwx,npwx)
integer :: pwsymndx(npwc,2*nsym)
integer :: ipw2,jw,iqcentr,jqcentr
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
integer :: nsing
double precision :: wsing(20)
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
double precision :: aa,bb,cc,dd,ee,ff
integer :: iter 
common /fn/ aa,bb,cc,dd,ee,ff,iter 
data iter /30/  

abr=1.d-6
rlr=1.d-6
sei=0.d0
ctest=(0.d0,0.d0)
cse=(0.d0,0.d0)
dcse=(0.d0,0.d0)
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
igsymk=isymg(:,ikpt,isym)
ek=enrgy(indxkbnd(ikpt)+iband)
do ii=1,3
do jj=1,3
  qkcvt(ii,jj)=blat(ii,jj)/dble(ngkpt(ii))
enddo
enddo
ictr=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
iqctr=iqndx(ictr(1),ictr(2),ictr(3))

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
do ibp=test_bands_se(1),test_bands_se(2)
!do ibp=2,2
!do ibp=iband,iband
  write(6,'(3(a,i4))') 'quasiparticle band ',iband,', polarized band ',ibp,' out of ',test_bands_se(2)
  seibc=sei
  cse2=cse1
  xse2=xse
  if (ibp.gt.nbocc) then
    ocsign=-1
  else
    ocsign=1
  endif
  do iipw=1,napwndx
!  do iipw=napwndx,napwndx
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
!    write(6,'(1x,i3,4x,2i4)') ibp,ipw1,ipw2 
    lx=ipw1.eq.ipw2.and.ibp.le.nbocc  ! do exchange
    lc=iipw.le.ntpwndx  ! do correlation
    if ((.not.lx).and.(.not.lc)) cycle

! Step 1: Compute matrix elements
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=10,10
!    do iz=10,10
      iqq=(/ix,iy,iz/)
!write(6,*) '>>>>> iqq = ',iqq
!write(6,'(a,3f10.5)') 'kpt = ',kpt(:,ikpt)
      iks=iqq-ngkpt/2
      jka=ikk-iks
      jkk=mod(jka,ngkpt)
      do ii=1,3
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
      enddo
      ikptq=ikndxq(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
!      if (.true.) then
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
      xmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)
      vqmat2(ix,iy,iz)=xmat2(ix,iy,iz)*vq(ix,iy,iz)
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

    if (itetrahedron.eq.0) then
! sum over points - faster, less accurate
! doesn't work yet, needs debugging
      brd=dw ! broadening for sum
      wint0=(0.d0,0.d0)
      w2int0=(0.d0,0.d0)
      cse0=(0.d0,0.d0)
      dcse0=(0.d0,0.d0)
      xse0=0.d0

! determine if singular integration is required
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
        if (qadj(1).gt.1.d1.and.qadj(2).gt.1.d1.and.qadj(3).gt.1.d1.and.  &
&           qadj(4).gt.1.d1.and.qadj(5).gt.1.d1.and.qadj(6).gt.1.d1) then
          lqsing=.true.
          iqq=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
          iqpt=iqndx(iqq(1),iqq(2),iqq(3))
          call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
          if (jjpw.ne.0) then
            do iw=1,nwpt
              ww=iw*dw
              eval=omega(iqq(1),iqq(2),iqq(3))-ocsign*ww-dw/2
              wint0=wint0+dw*dimag(lossfn(iw,iqpt,jjpw))/cmplx(eval,brd)
              w2int0=w2int0-dw*dimag(lossfn(iw,iqpt,jjpw))/cmplx(eval,brd)**2
            enddo
          endif
        else
          lqsing=.false.
        endif
      else
        lqsing=.false.
      endif
      
      if (lc) then
        do ix=1,ngkpt(1)
        do iy=1,ngkpt(2)
        do iz=1,ngkpt(3)
          iqq=(/ix,iy,iz/)
          iqpt=iqndx(iqq(1),iqq(2),iqq(3))
          call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
          wint=(0.d0,0.d0)
          w2int=(0.d0,0.d0)
          do iw=1,nwpt
            ww=iw*dw
            eval=omega(ix,iy,iz)-ocsign*ww-dw/2
            wint=wint+dw*dimag(lossfn(iw,iqpt,jjpw))/cmplx(eval,brd)
            w2int=w2int-dw*dimag(lossfn(iw,iqpt,jjpw))/cmplx(eval,brd)**2
          enddo
          if (lqsing) then
            qq2=4.d0*pi/vq(ix,iy,iz)
            if (sqrt(qq2).lt.gfo) then
              if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) cycle
              gterm=xmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)*(1+cos(pi*sqrt(qq2)/gfo))/(2.d0)
              if (lx) xse0=xse0+vqmat2(ix,iy,iz)-gterm*vq(ix,iy,iz)
              cse0=cse0+wint*vqmat2(ix,iy,iz)-wint0*gterm*vq(ix,iy,iz)
              dcse0=dcse0+w2int*vqmat2(ix,iy,iz)-w2int0*gterm*vq(ix,iy,iz)
            else
              if (lx) xse0=xse0+vqmat2(ix,iy,iz)
              cse0=cse0+wint*vqmat2(ix,iy,iz)
              dcse0=dcse0+w2int*vqmat2(ix,iy,iz)
            endif
          else
            if (lx) xse0=xse0+vqmat2(ix,iy,iz)
            cse0=cse0+wint*vqmat2(ix,iy,iz)
            dcse0=dcse0+w2int*vqmat2(ix,iy,iz)
          endif
        enddo
        enddo
        enddo
        cse0=cse0/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
        dcse0=dcse0/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
        if (lx) xse0=xse0/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      elseif (lx) then
        do ix=1,ngkpt(1)
        do iy=1,ngkpt(2)
        do iz=1,ngkpt(3)
          if (lqsing) then
            qq2=4.d0*pi/vq(ix,iy,iz)
            if (sqrt(qq2).lt.gfo) then
              if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) cycle
              gterm=xmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)*(1+cos(pi*sqrt(qq2)/gfo))/(2.d0)
              xse0=xse0+vqmat2(ix,iy,iz)-gterm*vq(ix,iy,iz)
            else
              xse0=xse0+vqmat2(ix,iy,iz)
            endif
          else
            xse0=xse0+vqmat2(ix,iy,iz)
          endif
        enddo
        enddo
        enddo
        xse0=xse0/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      endif

      if (lqsing) then
        if (lx) xse0=xse0+xmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)*8.d0*pi*gfo/((2.d0*pi)**3)
        if (lc) cse0=cse0+wint0*xmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)*8.d0*pi*gfo/((2.d0*pi)**3)
        if (lc) dcse0=dcse0+w2int0*xmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)*8.d0*pi*gfo/((2.d0*pi)**3)
      endif

      if (lc) cse=cse-cse0
      if (lc) dcse=dcse-dcse0
      if (lx) xse=xse-xse0*pi

    else 
! tetrahedron integration

! Don't waste time integrating over energies with no contribution
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
              vfactor(iqq(1),iqq(2),iqq(3),iw)=vqmat2(iqq(1),iqq(2),iqq(3))*max(dimag(lossfn(iw,iqpt,jjpw)),0.d0)
            else
              vfactor(iqq(1),iqq(2),iqq(3),iw)=(0.d0,0.d0)
            endif
          enddo
          enddo
          enddo
        enddo
      endif

! determine if singular integration is required
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
        if (qadj(1).gt.1.d1.and.qadj(2).gt.1.d1.and.qadj(3).gt.1.d1.and.  &
&           qadj(4).gt.1.d1.and.qadj(5).gt.1.d1.and.qadj(6).gt.1.d1) then
          lqsing=.true.
        else
          lqsing=.false.
        endif
      else
        lqsing=.false.
      endif
  
! Do the integration
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=1,1
!      do iy=1,1
!      do iz=9,9
!      do ix=1,1
!      do iy=1,1
!      do iz=1,1
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
!        do itet=2,2
          ssx2=ssx
          lqcentr=.false.
!          lqcentr=.true.
          do iv=1,4
            evtx0(iv)=enval(ivndx(iv,itet))
            rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
            if (lc) then
              do iw=1,nwpt
                vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
              enddo
            endif
            if (lx) xpyr0(iv)=vfx(ivndx(iv,itet))
            if (lqsing) vqvtx0(iv)=vq2(ivndx(iv,itet))
            if (lqsing.and..not.lqcentr) then
              iqv=iqq+nint(rrpyr(:,iv))
              if (iqv(1).eq.ictr(1).and.iqv(2).eq.ictr(2).and.iqv(3).eq.ictr(3)) then
                lqcentr=.true.
                iqcentr=iv
              endif
            endif
          enddo
          call indxhpsort(4,4,evtx0,indxe)
          evtx=evtx0(indxe)    ! order vertices by energy
          jqcentr=indxe(iqcentr)
!          kvtx(1:3,1:4)=rrpyr(1:3,indxe)
          do ii=1,4
            jj=indxe(ii)
!            kvtx(:,ii)=rrpyr(:,jj)+iqq-ictr
            do kk=1,3
              kvtx(kk,ii)=dot_product(iqq+rrpyr(:,jj)-ictr,qkcvt(:,kk))
!              kvtx(kk,ii)=dot_product(iqq+rrpyr(:,jj)-ictr+ipwx(:,ipw1)*ngkpt,qkcvt(:,kk))
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
          iwwlim=nint(evtx/dw)
          de21=evtx(2)-evtx(1)
          de31=evtx(3)-evtx(1)
          de32=evtx(3)-evtx(2)
          de41=evtx(4)-evtx(1)
          de42=evtx(4)-evtx(2)
          de43=evtx(4)-evtx(3)
          thresh=1.d-5*de41
          if (lqcentr) then
            bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
            do jj=1,3
              bgrad=bgrad+(evtx(jj+1)-evtx(1))*rg(:,jj)
            enddo
            xmult=4.d0*pi/sqrt(dot_product(bgrad,bgrad))
            if (lc) then
              do iw=1,nwpt
                aa0(iw)=xmat2(ictr(1),ictr(2),ictr(3))*max(dimag(lossfn(iw,iqctr,1)),0.d0)
                vpyr0(:,iw)=vpyr0(:,iw)-aa0(iw)*vqvtx0(:)
                vpyr0(iqcentr,iw)=0.d0
              enddo
            endif
            if (lx) then
              xme=dble(xmat2(ictr(1),ictr(2),ictr(3)))
              xpyr0=xpyr0-xme*vqvtx0
              xpyr0(iqcentr)=0.d0
            endif
            do iww=iwwlim(1),iwwlim(2) ! energies between evtx(1) & evtx(2)
              wwwlo=max(dw*(iww-0.5d0),evtx(1))
              wwwhi=min(dw*(iww+0.5d0),evtx(2))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              call wsinggrater(wwwlo,wwwhi,1,0,xk,kvtx(:,1),evtx,abr,rlr,sint1)
              sint1=sint1*xmult
              if (lc) then
                do iw=1,nwpt
                  ssi(iww-ocsign*iw)=ssi(iww-ocsign*iw)+aa0(iw)*sint1
                enddo
              endif
              if (lx) ssx=ssx+xme*sint1
            enddo
            do iww=iwwlim(2),iwwlim(3) ! energies between evtx(2) & evtx(3)
              wwwlo=max(dw*(iww-0.5d0),evtx(2))
              wwwhi=min(dw*(iww+0.5d0),evtx(3))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              if (de21.lt.thresh.and.de43.lt.thresh) then
                call  w2singgrater(wwwlo,wwwhi,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
                sint1a=0.d0
              elseif (de21.ge.de43) then
                call wsinggrater(wwwlo,wwwhi,2,0,xk,kvtx(:,1),evtx,abr,rlr,sint1a)
                call wsinggrater(wwwlo,wwwhi,1,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
              else
                call wsinggrater(wwwlo,wwwhi,3,0,xk,kvtx(:,1),evtx,abr,rlr,sint1a)
                call wsinggrater(wwwlo,wwwhi,4,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
              endif
              sint1=(sint1b-sint1a)*xmult
              if (lc) then
                do iw=1,nwpt
                  ssi(iww-ocsign*iw)=ssi(iww-ocsign*iw)+aa0(iw)*sint1
                enddo
              endif
              if (lx) ssx=ssx+xme*sint1
            enddo
            do iww=iwwlim(3),iwwlim(4) ! energies between evtx(3) & evtx(4)
              wwwlo=max(dw*(iww-0.5d0),evtx(3))
              wwwhi=min(dw*(iww+0.5d0),evtx(4))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              call wsinggrater(wwwlo,wwwhi,4,0,xk,kvtx(:,1),evtx,abr,rlr,sint1)
              sint1=sint1*xmult
              if (lc) then
                do iw=1,nwpt
                  ssi(iww-ocsign*iw)=ssi(iww-ocsign*iw)+aa0(iw)*sint1
                enddo
              endif
              if (lx) ssx=ssx+xme*sint1
            enddo
          endif
!          if (lqcentr) then            ! integral around coulomb singularity
!          if (.true.) then
          if (.false.) then
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
&                           +(vpyr(indxe(jj+1))-vpyr(indxe(1)))*rg(1:3,jj)
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
            do iww=iwwlim(1),iwwlim(2) ! energies between evtx(1) & evtx(2)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(1))
              wwwhi=min(dw*(iww+0.5d0),evtx(2))
              call wsinggrater(wwwlo,wwwhi,1,0,xk,kvtx(:,1),evtx,abr,rlr,sint1)
write(6,*) iww
              call wsinggrater(wwwlo,wwwhi,1,1,xk,kvtx(:,1),evtx,abr,rlr,svec(1))
              call wsinggrater(wwwlo,wwwhi,1,2,xk,kvtx(:,1),evtx,abr,rlr,svec(2))
              call wsinggrater(wwwlo,wwwhi,1,3,xk,kvtx(:,1),evtx,abr,rlr,svec(3))
              sint1=sint1*xmult
              svec=svec*xmult
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') '1  T',itet,iww,ssx-ssx3
            enddo
            do iww=iwwlim(2),iwwlim(3) ! energies between evtx(2) & evtx(3)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(2))
              wwwhi=min(dw*(iww+0.5d0),evtx(3))
              if (de21.lt.thresh.and.de43.lt.thresh) then
label='2a T'
                call  w2singgrater(wwwlo,wwwhi,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
                call  w2singgrater(wwwlo,wwwhi,1,xk,kvtx(:,1),evtx,abr,rlr,svecb(1))
                call  w2singgrater(wwwlo,wwwhi,2,xk,kvtx(:,1),evtx,abr,rlr,svecb(2))
                call  w2singgrater(wwwlo,wwwhi,3,xk,kvtx(:,1),evtx,abr,rlr,svecb(3))
                sint1a=0.d0
                sveca=(/0.d0,0.d0,0.d0/)
              elseif (de21.ge.de43) then
label='2b T'
                call wsinggrater(wwwlo,wwwhi,2,0,xk,kvtx(:,1),evtx,abr,rlr,sint1a)
                call wsinggrater(wwwlo,wwwhi,2,1,xk,kvtx(:,1),evtx,abr,rlr,sveca(1))
                call wsinggrater(wwwlo,wwwhi,2,2,xk,kvtx(:,1),evtx,abr,rlr,sveca(2))
                call wsinggrater(wwwlo,wwwhi,2,3,xk,kvtx(:,1),evtx,abr,rlr,sveca(3))
                call wsinggrater(wwwlo,wwwhi,1,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
                call wsinggrater(wwwlo,wwwhi,1,1,xk,kvtx(:,1),evtx,abr,rlr,svecb(1))
                call wsinggrater(wwwlo,wwwhi,1,2,xk,kvtx(:,1),evtx,abr,rlr,svecb(2))
                call wsinggrater(wwwlo,wwwhi,1,3,xk,kvtx(:,1),evtx,abr,rlr,svecb(3))
              else
label='2c T'
                call wsinggrater(wwwlo,wwwhi,3,0,xk,kvtx(:,1),evtx,abr,rlr,sint1a)
                call wsinggrater(wwwlo,wwwhi,3,1,xk,kvtx(:,1),evtx,abr,rlr,sveca(1))
                call wsinggrater(wwwlo,wwwhi,3,2,xk,kvtx(:,1),evtx,abr,rlr,sveca(2))
                call wsinggrater(wwwlo,wwwhi,3,3,xk,kvtx(:,1),evtx,abr,rlr,sveca(3))
                call wsinggrater(wwwlo,wwwhi,4,0,xk,kvtx(:,1),evtx,abr,rlr,sint1b)
                call wsinggrater(wwwlo,wwwhi,4,1,xk,kvtx(:,1),evtx,abr,rlr,svecb(1))
                call wsinggrater(wwwlo,wwwhi,4,2,xk,kvtx(:,1),evtx,abr,rlr,svecb(2))
                call wsinggrater(wwwlo,wwwhi,4,3,xk,kvtx(:,1),evtx,abr,rlr,svecb(3))
              endif
              sint1=(sint1b-sint1a)*xmult
              svec=(svecb-sveca)*xmult
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') label,itet,iww,ssx-ssx3
            enddo
            do iww=iwwlim(3),iwwlim(4) ! energies between evtx(3) & evtx(4)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(3))
              wwwhi=min(dw*(iww+0.5d0),evtx(4))
              call wsinggrater(wwwlo,wwwhi,4,0,xk,kvtx(:,1),evtx,abr,rlr,sint1)
              call wsinggrater(wwwlo,wwwhi,4,1,xk,kvtx(:,1),evtx,abr,rlr,svec(1))
              call wsinggrater(wwwlo,wwwhi,4,2,xk,kvtx(:,1),evtx,abr,rlr,svec(2))
              call wsinggrater(wwwlo,wwwhi,4,3,xk,kvtx(:,1),evtx,abr,rlr,svec(3))
              sint1=sint1*xmult
              svec=svec*xmult
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') '3  T',itet,iww,ssx-ssx3
            enddo
          else                         ! non-singular, tetrahedron method
            fbx(1)=tvol/(6.d0*de21*de31*de41)
            if (de21.ge.de43) then
              fbx(2)=tvol/(6.d0*de21*de32*de42)
            else
              fbx(3)=tvol/(6.d0*de31*de32*de43)
            endif
            fbx(4)=tvol/(6.d0*de41*de42*de43)
            do ii=1,4
              cmx(1:3,ii)=(/0.d0,0.d0,0.d0/)
              do jj=1,4
                if (jj.ne.ii) cmx(1:3,ii)=cmx(1:3,ii)+(xk(1:3,jj)-xk(1:3,ii)) &
&                                                    /(4*(evtx(jj)-evtx(ii)))
              enddo
            enddo
            if (lc) then
              do iw=1,nwpt
                aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
                av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
                do jj=1,3
                  av(1:3,iw)=av(1:3,iw) &
&                           +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
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
            do iww=iwwlim(1),iwwlim(2) ! energies between evtx(1) & evtx(2)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(1))
              wwwhi=min(dw*(iww+0.5d0),evtx(2))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              sint1=fbx(1)*((wwwhi-evtx(1))**3-(wwwlo-evtx(1))**3)
              svec=xk(:,1)*sint1 &
&                 +fbx(1)*cmx(:,1)*((wwwhi-evtx(1))**4-(wwwlo-evtx(1))**4)
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') '1  F',itet,iww,ssx-ssx3
            enddo
            do iww=iwwlim(2),iwwlim(3) ! energies between evtx(2) & evtx(3)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(2))
              wwwhi=min(dw*(iww+0.5d0),evtx(3))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              if (de21.lt.thresh.and.de43.lt.thresh) then
label='2a F'
                sint1=(tvol/(6.d0*de31**2)) &
&                    *(3.d0*((wwwhi-evtx(1))**2-(wwwlo-evtx(1))**2) &
&                     -2.d0*((wwwhi-evtx(1))**3-(wwwlo-evtx(1))**3)/de31)
                svec=xk(:,2)*sint1 &
&                   -(tvol/(24.d0*de31**4)) &
&                    *((wwwhi-evtx(1))**4-(wwwlo-evtx(1))**4) &
&                    *(xk(:,3)+xk(:,4)-2.d0*xk(:,2))
              elseif (de21.ge.de43) then
label='2b F'
                fb(1)=fbx(1)*((wwwhi-evtx(1))**3-(wwwlo-evtx(1))**3)
                fb(2)=fbx(2)*((wwwhi-evtx(2))**3-(wwwlo-evtx(2))**3)
                sint1=fb(1)-fb(2)
                svec=xk(:,1)*fb(1)-xk(:,2)*fb(2) &
&                   +fbx(1)*cmx(:,1)*((wwwhi-evtx(1))**4-(wwwlo-evtx(1))**4) &
&                   -fbx(2)*cmx(:,2)*((wwwhi-evtx(2))**4-(wwwlo-evtx(2))**4)
              else
label='2c F'
                fb(3)=fbx(3)*((wwwhi-evtx(3))**3-(wwwlo-evtx(3))**3)
                fb(4)=fbx(4)*((wwwhi-evtx(4))**3-(wwwlo-evtx(4))**3)
                sint1=fb(4)-fb(3)
                svec=xk(:,4)*fb(4)-xk(:,3)*fb(3) &
&                   +fbx(4)*cmx(:,4)*((wwwhi-evtx(4))**4-(wwwlo-evtx(4))**4) &
&                   -fbx(3)*cmx(:,3)*((wwwhi-evtx(3))**4-(wwwlo-evtx(3))**4)
              endif
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') label,itet,iww,ssx-ssx3
            enddo
            do iww=iwwlim(3),iwwlim(4) ! energies between evtx(3) & evtx(4)
!ssx3=ssx
              wwwlo=max(dw*(iww-0.5d0),evtx(3))
              wwwhi=min(dw*(iww+0.5d0),evtx(4))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              sint1=fbx(4)*((wwwhi-evtx(4))**3-(wwwlo-evtx(4))**3)
              svec=xk(:,4)*sint1 &
&                 +fbx(4)*cmx(:,4)*((wwwhi-evtx(4))**4-(wwwlo-evtx(4))**4)
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))
                enddo
              endif
              if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))
!write(6,'(a,2i3,f16.8)') '3  F',itet,iww,ssx-ssx3
            enddo
          endif
!if (lqcentr) then
!write(48,'(3i3,3x,i3,f16.8,3x,a)') ix,iy,iz,itet,ssx-ssx2,'T'
!else
!write(48,'(3i3,3x,i3,f16.8,3x,a)') ix,iy,iz,itet,ssx-ssx2,'F'
!endif
        enddo
!      temparray(ix,iy,iz)=(ssx-ssx1)/((2*pi)**3)
      enddo
      enddo
      enddo
!    if (lc) ssi=-ssi/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
!    if (lx) ssx=-ssx/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
      if (lc) ssi=-ocsign*ssi/((2*pi)**3*pi)
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
              ssc(ie)=ssc(ie)+ocsign*ssi(je)*dw/(eps2-eps1)
            endif
!if (ie.eq.0) write(16,'(2i4,2f10.3,es12.3,2(2x,2es11.3))') ie,je,eps1*27.2114,eps2*27.2114,1/(eps2-eps1),ssc(ie)*27.2114,dble(ssi(je)*27.2114)
          enddo
        enddo
      endif
      deallocate (ssi)

      if (lc) cse1=cse1+ssc
      if (lx) xse=xse+ssx
!      if (lx.and.lc) then
!        write(6,'(i4,4x,i5,2i3,3x,3i3,5x,"(",f10.6,",",f10.6,")",5x,f10.6)')  &
!&            ibp,iipw,ipw1,ipw2,ipwx(:,ipw1),((ssc(1)+ssc(0))/2)*27.2114,dble(ssx)*27.2114
!      elseif (lc) then
!        write(6,'(i4,4x,i5,2i3,3x,3i3,5x,"(",f10.6,",",f10.6,")")') ibp,iipw,ipw1,ipw2,ipwx(:,ipw1), &
!&            ((ssc(1)+ssc(0))/2)*27.2114
!      else
!        write(6,'(i4,4x,i5,2i3,3x,3i3,28x,5x,f10.6)')  &
!&            ibp,iipw,ipw1,ipw2,ipwx(:,ipw1),dble(ssx)*27.2114
!      endif
    endif

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
if (itetrahedron.eq.1) then
  cse=(cse1(0)+cse1(1))/2
  dcse=(cse1(1)-cse1(0))/dw
endif
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
end subroutine mk3cse

!****************************************************************************
! evaluate Integral_wmin^wmax dww sing1int

subroutine wsinggrater(wmin,wmax,ive,isn,xk,xkvtx,evtx,abr,rlr,sint)
implicit none
integer :: numcal,maxns,mx,nstack,ii,jj,kk,icount,ive,isn
parameter (mx=1500)
double precision :: xk(3,4),xkvtx(3),evtx(4),sint
double precision :: value,valu,fval(3,mx)
double precision :: wmin,wmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: wleft(mx),dw(3),wt(3),ww
double precision :: wt9(9)
logical :: atsing
save dw,wt,wt9
data dw /0.1127016653792583,0.5,0.8872983346207417/
data wt /0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9 /0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater) and the region
! is removed from consideration (and a new region is the right-most)
  nstack=1
  maxns=nstack
  error=0.d0
  sint=0.d0
! The array xleft stores the boundary points of the regions.
  wleft(1)=wmin
  wleft(2)=wmax
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=wleft(ii+1)-wleft(ii)
    do jj=1,3
      ww=wleft(ii)+del*dw(jj)
      call sing1int(ive,isn,ww,xk,xkvtx,evtx,fval(jj,ii))
    enddo
  enddo
  numcal=nstack*3
  do
    if(nstack+3.ge.mx) then
      write(6,*) 'wsinggrater: too many regions'
      stop
    endif
! Divide the rightmost region into three subregions.
    del=wleft(nstack+1)-wleft(nstack)
    wleft(nstack+3)=wleft(nstack+1)
    wleft(nstack+1)=wleft(nstack)+del*dw(1)*2.d0
    wleft(nstack+2)=wleft(nstack+3)-del*dw(1)*2.d0
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=0.d0
    valu=0.d0
    do jj=nstack,nstack+2
      del1=wleft(jj+1)-wleft(jj)
      ww=wleft(jj)+dw(1)*del1
      call sing1int(ive,isn,ww,xk,xkvtx,evtx,fval(1,jj))
      ww=wleft(jj)+dw(3)*del1
      call sing1int(ive,isn,ww,xk,xkvtx,evtx,fval(3,jj))
      numcal=numcal+2
!if (numcal.gt.1000) stop
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac=del/(wmax-wmin)
    atsing=.false.
    if (frac.le.1.0d-8) atsing=.true.
    if (dif.le.abr*frac.or.dif.le.rlr*abs(value).or. &
&        (atsing.and.(frac.le.1.0d-15.or.dif.le.abr*0.1))) then
      sint=sint+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo
end subroutine wsinggrater

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

subroutine sing1int(ive,isn,ww,xk,xkvtx,evtx,sint)
implicit none
integer :: ive,isn
double precision :: xk(3,4),xkvtx(3),evtx(4),ww
double precision :: xkp(3,3),xq1(3),xq2(3),xkv0(3),xcvt
double precision :: aa,bb,cc,dd,ee,ff,qq,delta
double precision :: xmin,xmax,abr,rlr,xsing(20),error
double precision :: xfn1,xfn2,xfn3,grater
double precision :: sint,root1,root2,sint1,sint2,sint3
double precision :: xnorm(3),area
double precision :: xdum,xmx,rlxk
integer :: ii,jj,kk,ll
integer :: nsing,numcal,maxns,nroots
integer :: iter 
common /fn/ aa,bb,cc,dd,ee,ff,iter 
external xfn1,xfn2,xfn3,grater

! broadening: 1/|q|^2 -> 1/(|q|^2+delta^2)
  delta=min(dot_product(xk(:,2),xk(:,2)),dot_product(xk(:,3),xk(:,3)),dot_product(xk(:,4),xk(:,4)))*1.d-4
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
  xq1=xkp(:,1)-xkp(:,3)
  xq2=xkp(:,2)-xkp(:,3)
  xkv0=xkp(:,3)+xkvtx
  call cross(xq1,xq2,xnorm)
  area=sqrt(dot_product(xnorm,xnorm))
  aa=dot_product(xq1,xq1)
  bb=2*dot_product(xq1,xkv0)
  cc=2*dot_product(xq1,xq2)
  dd=2*dot_product(xq2,xkv0)
  ee=dot_product(xq2,xq2)
  ff=dot_product(xkv0,xkv0)+delta**2
!  ff=dot_product(xkv0,xkv0)
  nsing=0
  call rquadroots(ee,dd,ff,nroots,root1,root2)
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
  call rquadroots(1.d0+cc/2+ee,bb/2+dd,ff,nroots,root1,root2)
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
  call rquadroots(4.d0*aa*ee-cc**2,4.d0*aa*dd-2.d0*bb*cc,4.d0*aa*ff-bb**2,nroots,root1,root2)
  if (nroots.gt.0.and.root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
  endif
  if (nroots.gt.1.and.root2.gt.0.d0.and.root2.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root2
  endif
  call hpsort(nsing,nsing,xsing(1:nsing))
  xmin=0.d0
  xmax=1.d0
  abr=1.d-12
  rlr=1.d-6
  sint1=grater(xfn1,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  if (isn.ge.1.and.isn.le.3) then
    sint2=grater(xfn2,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
    sint3=grater(xfn3,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
    sint=xkp(isn,ll)*sint1+xq1(isn)*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2(isn)*sint3
  else
    sint=sint1
  endif
!  svec=xkp(:,3)*sint1+xq1*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2*sint3
!  iter=iter+1 
end subroutine sing1int


!****************************************************************************
! evaluate Integral_wmin^wmax dww sing2int

subroutine w2singgrater(wmin,wmax,isn,xk,xkvtx,evtx,abr,rlr,sint)
implicit none
integer :: numcal,maxns,mx,nstack,ii,jj,kk,icount,isn
parameter (mx=1500)
double precision :: xk(3,4),xkvtx(3),evtx(4),sint
double precision :: value,valu,fval(3,mx)
double precision :: wmin,wmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: wleft(mx),dw(3),wt(3),ww
double precision :: wt9(9)
logical :: atsing
save dw,wt,wt9
data dw /0.1127016653792583,0.5,0.8872983346207417/
data wt /0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9 /0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater) and the region
! is removed from consideration (and a new region is the right-most)
  nstack=1
  maxns=nstack
  error=0.d0
  sint=0.d0
! The array xleft stores the boundary points of the regions.
  wleft(1)=wmin
  wleft(2)=wmax
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=wleft(ii+1)-wleft(ii)
    do jj=1,3
      ww=wleft(ii)+del*dw(jj)
      call sing2int(isn,ww,xk,xkvtx,evtx,fval(jj,ii))
    enddo
  enddo
  numcal=nstack*3
  do
    if(nstack+3.ge.mx) then
      write(6,*) 'wsinggrater: too many regions'
      stop
    endif
! Divide the rightmost region into three subregions.
    del=wleft(nstack+1)-wleft(nstack)
    wleft(nstack+3)=wleft(nstack+1)
    wleft(nstack+1)=wleft(nstack)+del*dw(1)*2.d0
    wleft(nstack+2)=wleft(nstack+3)-del*dw(1)*2.d0
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=0.d0
    valu=0.d0
    do jj=nstack,nstack+2
      del1=wleft(jj+1)-wleft(jj)
      ww=wleft(jj)+dw(1)*del1
      call sing2int(isn,ww,xk,xkvtx,evtx,fval(1,jj))
      ww=wleft(jj)+dw(3)*del1
      call sing2int(isn,ww,xk,xkvtx,evtx,fval(3,jj))
      numcal=numcal+2
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac=del/(wmax-wmin)
    atsing=.false.
    if (frac.le.1.0d-8) atsing=.true.
    if (dif.le.abr*frac.or.dif.le.rlr*abs(value).or. &
&        (atsing.and.(frac.le.1.0d-15.or.dif.le.abr*0.1))) then
      sint=sint+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo
end subroutine w2singgrater

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

subroutine sing2int(isn,xk,xkvtx,evtx,sint)
implicit none
integer :: isn
double precision :: xk(3,4),xkvtx(3),evtx(4),ww
double precision :: xkp(3,3),xq1(3),xq2(3),xkv0(3),xcvt
double precision :: aa,bb,cc,dd,ee,ff,qq,delta
double precision :: xmin,xmax,abr,rlr,xsing(20),error
double precision :: xfn4,xfn5,xfn6,grater
double precision :: sint1,sint2,sint3,svec(3),root1,root2,sint
double precision :: xnorm(3),area
double precision :: xdum
integer :: ii,jj,ll
integer :: nsing,numcal,maxns,nroots
integer :: iter 
common /fn/ aa,bb,cc,dd,ee,ff,iter 
external xfn4,xfn5,xfn6,grater

! broadening: 1/|q|^2 -> 1/(|q|^2+delta^2)
  delta=min(dot_product(xk(:,2),xk(:,2)),dot_product(xk(:,3),xk(:,3)),dot_product(xk(:,4),xk(:,4)))*1.d-4
  xkp(:,1)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
  xkp(:,2)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  xkp(:,3)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
  xq1=xkp(:,1)-xkp(:,3)
  xq2=xkp(:,2)-xkp(:,3)
  xkv0=xkp(:,3)+xkvtx
  call cross(xq1,xq2,xnorm)
  area=sqrt(dot_product(xnorm,xnorm))
  aa=dot_product(xq1,xq1)
  bb=2*dot_product(xq1,xkv0)
  cc=2*dot_product(xq1,xq2)
  dd=2*dot_product(xq2,xkv0)
  ee=dot_product(xq2,xq2)
  ff=dot_product(xkv0,xkv0)+delta**2
  nsing=0
  call rquadroots(ee,dd,ff,nroots,root1,root2)
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
  call hpsort(nsing,nsing,xsing(1:nsing))
  xmin=0.d0
  xmax=1.d0
  abr=1.d-12
  rlr=1.d-6
  sint1=grater(xfn4,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
  if (isn.ge.1.and.isn.le.3) then
    sint2=grater(xfn5,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
    sint3=grater(xfn6,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
    sint=xkp(isn,3)*sint1+xq1(isn)*(sint2-bb*sint1-cc*sint3)/(2*aa)+xq2(isn)*sint3
  else
    sint=sint1
  endif

end subroutine sing2int
