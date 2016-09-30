subroutine mkcse(iband,jband,ikk,W, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& pwmatel,tpwmatel, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nqpt,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,bmet,blat,ipaw, &
& ipwlf,npwlf,ipwndx,npwndx,ntpwndx, &
& npwup,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,iqndx,isymndx,isymndxq,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,cse,xse)
implicit none
integer :: iband,jband,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,npwtq,bantot,bantotq,ncg,ncgq,nkpt,nkptq,nqpt,nsym,nsymq,nsppol
integer :: npw,npwq,ipw1,npwlf,npwndx,ntpwndx,ipaw
double precision :: vol,pi,wmax,xred(3,natom)
double complex :: W(nwpt,nqpt,ntpwndx)
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
integer :: ipwlf(3,npwlf),ipwndx(2,ntpwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex :: ssi(nwpt),ssc(nwpt),sec(nwpt),xse1,xse2,cse(nbcore+1:ncband),cse1(nbcore+1:ncband),zz(nbcore+1:ncband),zz1(nbcore+1:ncband),ctest
double precision :: xse
integer :: ii,jj,kk,ix,iy,iz,iskip,isign,iloss
integer :: iqq(3),iqqp(3),jka(3),jkk(3),ikks(3)
integer :: ikpt,ikptq,iks(3),ikslf1(3),ikslf2(3)
integer :: igg(3),igh(3),igg0(3),icenter,isym,isymq,ibp,iqpt,iqsym,iw,ie,je,ies
double precision :: xck(3),xckq(3)
double complex :: cmatel,cmatel2,amatel,amatel2,vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3)),vfv(8),lossv(8)
double precision :: vq2(8)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,enval(8),dw,eshift
double precision :: abr,rlr
double complex :: tseincr,tseincrx,dse,dcse
integer :: iwh,iwl,ibmin,ibmax,ioff,ioff2
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),qq(3),qp(3),qq2,qp2,qs(3),xk(3),xkmq(3),ek,ekmq
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),temparray(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: npwup,iipw,jjpw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
integer :: ipw2,jw
double precision :: eps1(nbcore+1:ncband),eps2,eps
double precision :: eval,brd,esprd
double precision :: gfo,gamma(3),pln(3),dist(3)
double complex :: wint(nbcore+1:ncband),wint0(nbcore+1:ncband)
double complex :: w2int(nbcore+1:ncband),w2int0(nbcore+1:ncband)
double complex :: cmgamma2,gterm
double complex :: pwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqsing
integer :: idum
logical :: lx,lc
double precision :: xdum
double precision :: linterp
external linterp

abr=1.d-6
rlr=1.d-6
sei=0.d0
xse=(0.d0,0.d0)
ctest=(0.d0,0.d0)
cse=(0.d0,0.d0)
zz=(0.d0,0.d0)
igg0=(/0,0,0/)
do ie=1,nwpt
  sec(ie)=(0.d0,0.d0)
enddo
dw=wmax/dble(nwpt)
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
!isym=isymndx(ikk(1),ikk(2),ikk(3))
isym=1 ! symtest
ek=enrgy(indxkbnd(ikpt)+iband)

do ix=1,ngkpt(1)
do iy=1,ngkpt(2)
do iz=1,ngkpt(3)
  temparray(ix,iy,iz)=(0.d0,0.d0)
enddo
enddo
enddo
ikks=(/0,0,0/)   ! symtest
do ii=1,3   ! symtest
do jj=1,3   ! symtest
  ikks(jj)=ikks(jj)+(ikk(ii)-ngkpt(ii)/2)*symrel(ii,jj,isym)   ! symtest
enddo   ! symtest
enddo   ! symtest
ikks=ikks+ngkpt/2   ! symtest
!write(6,*) 'ikk = ',ikk
!write(6,*) 'ikks = ',ikks

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

!write(6,'(a)') '      finding contibutions from band:      plane waves:'
do ibp=nbcore+1,ncband
!do ibp=1,1
  seibc=sei
  if (ibp.gt.nbocc) then
    isign=-1
  else
    isign=1
  endif
  if (ibp.ge.ibmin.and.ibp.le.ibmax) then
    iloss=1
  else
    iloss=0
  endif
!  do iipw=1,ntpwndx
  do iipw=1,1
!  do iipw=2,48,46
!  do iipw=2,2
!  do iipw=55,55
!    write(6,'(a,i4)') 'iipw ',iipw
!    write(6,'(38x,i3,14x,2i3)') ibp,ipw1,ipw2
    seipw=sei
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    if (ipw1.ne.ipw2) cycle
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do iy=2,2
!    do ix=3,3
!    do iz=4,4
      iqq=(/ix,iy,iz/)
      do ii=1,3
        iks(ii)=iqq(ii)-ngkpt(ii)/2
        ikslf1(ii)=iks(ii)+ipwlf(ii,ipw1)*ngkpt(ii)
        ikslf2(ii)=iks(ii)+ipwlf(ii,ipw2)*ngkpt(ii)
      enddo
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
      vq(ix,iy,iz)=4.d0*pi/sqrt(qq2*qp2)
!write(6,*) '>>>>> iqq = ',iqq
!write(6,*) 'ikslf1 = ',ikslf1
!write(6,*) 'ikslf2 = ',ikslf2
!write(6,'(a,3f10.6)') ' qq = ',qq
!write(6,'(a,3f10.6)') ' qp = ',qp
!write(6,*) 'qq2 = ',qq2
!write(6,*) 'qp2 = ',qp2
!write(6,*) 'vq = ',vq(ix,iy,iz)
      do ii=1,3
!        jka(ii)=ikk(ii)-iks(ii)
        jka(ii)=ikks(ii)-iks(ii)  ! symtest
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
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))
        igh(ii)=nint(xckq(ii)-xck(ii)+qp(ii))
      enddo
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
        ikptq=ikndxq(jkk(1),jkk(2),jkk(3))
        isymq=isymndxq(jkk(1),jkk(2),jkk(3))
      else
        ikptq=ikndx(jkk(1),jkk(2),jkk(3))
        isymq=isymndx(jkk(1),jkk(2),jkk(3))
      endif
!write(6,*) 'jka = ',jka
!write(6,*) 'jkk = ',jkk
!write(6,'(a,3f10.6)') ' xck = ',xck
!write(6,'(a,3f10.6)') ' xckq = ',xckq
!write(6,*) 'igg = ',igg
!write(6,*) 'igh = ',igh
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
        call mkmatelS(ibp,iband,ikpt,ikptq,ipwlf(:,ipw1),igg, &
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
        if (ipw1.eq.ipw2.and.iband.eq.jband) then
          cmatel2=cmatel
          amatel2=amatel
        else
          call mkmatelS(ibp,jband,ikpt,ikptq,ipwlf(:,ipw2),igh, &
&         ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&         isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&         cmatel2)
          if (ipaw.ne.0) then
!            call mkmatelP(pi,xred,natom,ibp,jband,ikpt,ikptq,qp,igg0,ngkpt, &
!&                 pwmatel(:,:,ipw2,ix,iy,iz),tpwmatel(:,:,ipw2,ix,iy,iz), &
!&                 projwf,nlmn,nkpt,nkptq,ncband,amatel2)
          else
            amatel2=(0.d0,0.d0)
          endif
        endif
        if (ibp.gt.nbocc) then
          omega(ix,iy,iz)=-enrgyq(indxkbndq(ikptq)+ibp)
        else
          omega(ix,iy,iz)=enrgyq(indxkbndq(ikptq)+ibp)
        endif
      else
        call mkmatelS(ibp,iband,ikpt,ikptq,ipwlf(:,ipw1),igg, &
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
        if (ipw1.eq.ipw2.and.iband.eq.jband) then
          cmatel2=cmatel
          amatel2=amatel
        else
          call mkmatelS(ibp,jband,ikpt,ikptq,ipwlf(:,ipw2),igh, &
&         ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&         isym,isymq,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&         cmatel2)
          if (ipaw.ne.0) then
!            call mkmatelP(pi,xred,natom,ibp,jband,ikpt,ikptq,qp,igg0,ngkpt, &
!&                 pwmatel(:,:,ipw2,ix,iy,iz),tpwmatel(:,:,ipw2,ix,iy,iz), &
!&                 projwf,nlmn,nkpt,nkpt,ncband,amatel2)
          else
            amatel2=(0.d0,0.d0)
          endif
        endif
        if (ibp.gt.nbocc) then
          omega(ix,iy,iz)=-enrgy(indxkbnd(ikptq)+ibp)
        else
          omega(ix,iy,iz)=enrgy(indxkbnd(ikptq)+ibp)
        endif
      endif
      vqmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)*vq(ix,iy,iz)
      if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) then
        cmgamma2=(cmatel+amatel)*conjg(cmatel2+amatel2)
      endif
!write(6,*) 'cmatel  = ',cmatel
!write(6,*) 'cmatel2 = ',cmatel2
!write(6,*) '|cmatel+amatel|^2 = ',(cmatel+amatel)*conjg(cmatel2+amatel2)
!write(6,*) 'vqmat2 = ',vqmat2(iqq(1),iqq(2),iqq(3))
!write(56,*) '>>>>>',iqq,vqmat2(iqq(1),iqq(2),iqq(3))
!      stvec(iqq(3))=vq(ix,iy,iz)
!      stvec(iqq(3))=sqrt(qq2)
      stvec(iqq(3))=dble(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec2(iqq(3))=dimag(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec(iqq(3))=dble(cmatel)
!      stvec2(iqq(3))=dimag(cmatel)
!      stvec(iqq(3))=dble(cmatel2)
!      stvec2(iqq(3))=dimag(cmatel2)
!      stvec(iqq(3))=dble(cmatel*conjg(cmatel))
      temparray(ix,iy,iz)=stvec(iqq(3))
    enddo
!    write(76,'(10f8.2)') stvec(:)
!    write(76,'(10es8.1)') stvec(:)
!    write(76,'(10es8.1)') stvec2(:)
!    write(76,*) 
    enddo
!    write(76,*) iqq(1)+1,'--------------------------------------------------',ibp
    enddo
!stop
! Cut material for testing energy dependence of self energy at bottom of file

    lx=ipw1.eq.ipw2.and.ibp.le.nbocc
    if (lx) xse1=(0.d0,0.d0)
    do ie=nbcore+1,ncband
      cse1(ie)=(0.d0,0.d0)
      eps=enrgy(indxkbnd(ikpt)+ie)
      if (ibp.gt.nbocc) then
        eps1(ie)=eps
      else
        eps1(ie)=-eps
      endif
    enddo
    brd=dw
!write(6,*) brd*27.2114
    if (abs(vqmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)/vqmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2-1)).gt.1.d3 &
&   .and.(ipw1.eq.1.or.ipw2.eq.1)) then
      iqsing=1
      iqq=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
      iqpt=iqndx(iqq(1),iqq(2),iqq(3))
!      call fespread(iqq,ngkpt,omega,esprd)
!      brd=max(esprd,dw)
      call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
      do ie=nbcore+1,ncband
        wint0(ie)=(0.d0,0.d0)
        if (jjpw.ne.0) then
          do iw=1,nwpt
            ww=iw*dw
            eval=omega(iqq(1),iqq(2),iqq(3))+eps1(ie)-ww
            wint0(ie)=wint0(ie)+dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)
            w2int0(ie)=w2int0(ie)-dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)**2
          enddo
        endif
      enddo
    else
      iqsing=0
    endif
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
      iqq=(/ix,iy,iz/)
      xse2=xse1
      iqpt=iqndx(ix,iy,iz)
!      call fespread(iqq,ngkpt,omega,esprd)
!      brd=max(esprd,dw)
!write(6,'(3i3,es20.10)') ix,iy,iz,esprd*27.2114
      call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
      do ie=nbcore+1,ncband
        wint(ie)=(0.d0,0.d0)
        w2int(ie)=(0.d0,0.d0)
        if (jjpw.ne.0) then
          do iw=1,nwpt
            ww=iw*dw
            eval=omega(ix,iy,iz)+eps1(ie)-ww
            wint(ie)=wint(ie)+dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)
            w2int(ie)=w2int(ie)-dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)**2
          enddo
        endif
      enddo
      if (iqsing.eq.0) then
        if (lx) xse1=xse1+vqmat2(ix,iy,iz)
        cse1=cse1+wint*vqmat2(ix,iy,iz)
        zz1=zz1+w2int*vqmat2(ix,iy,iz)
      else
        qq2=4.d0*pi/vq(ix,iy,iz)
        if (sqrt(qq2).lt.gfo) then
          if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) cycle
          gterm=cmgamma2*(1+cos(pi*sqrt(qq2)/gfo))/(2.d0)
          if (lx) xse1=xse1+vqmat2(ix,iy,iz)-gterm*vq(ix,iy,iz)
          cse1=cse1+wint*vqmat2(ix,iy,iz)-wint0*gterm*vq(ix,iy,iz)
          zz1=zz1+w2int*vqmat2(ix,iy,iz)-w2int0*gterm*vq(ix,iy,iz)
!write(6,'(3i3,2f12.6)') ix,iy,iz,(wint*vqmat2(ix,iy,iz)-wint0*gterm*vq(ix,iy,iz))/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
!write(21,'(3i3,2f14.2,4x,2f14.2)') ix,iy,iz,wint*vqmat2(ix,iy,iz),wint0*gterm*vq(ix,iy,iz)
!write(22,'(3i3,2f14.2,4x,2f14.2)') ix,iy,iz,vqmat2(ix,iy,iz),gterm*vq(ix,iy,iz)
!write(23,'(3i3,2f14.2,4x,2f14.2)') ix,iy,iz,wint,wint0
        else
          if (lx) xse1=xse1+vqmat2(ix,iy,iz)
          cse1=cse1+wint*vqmat2(ix,iy,iz)
          zz1=zz1+w2int*vqmat2(ix,iy,iz)
!write(6,'(3i3,2f12.6)') ix,iy,iz,wint*vqmat2(ix,iy,iz)/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
        endif
      endif
!      temparray(ix,iy,iz)=(xse1-xse2)/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
    enddo
    enddo
    enddo
    if (lx) xse1=xse1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
    cse1=cse1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    zz1=zz1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    if (iqsing.eq.1) then
! 1/pi Integral (d^3q/(2 pi)^3) wint0 gterm vq
      if (lx) xse1=xse1+cmgamma2*8.d0*pi*pi*gfo/((2.d0*pi)**3)
      cse1=cse1+wint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
      zz1=zz1+w2int0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
!write(6,*) wint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
!write(6,*) cmgamma2
!write(6,*) wint0
    endif
    xse1=-xse1
    cse1=-cse1*isign
    zz1=-zz1*isign
    if (lx) xse=xse+dble(xse1)
!    if (lx) ctest=ctest+xse1
    cse=cse+cse1
    zz=zz+zz1
!if (lx) then 
!  write(6,'(i4,4x,2i3,5x,"(",f10.6,",",f10.6,")",5x,f10.6)') ibp,ipw1,ipw2,cse1*27.2114,dble(xse1)*27.2114
!else
!  write(6,'(i4,4x,2i3,5x,"(",f10.6,",",f10.6,")")') ibp,ipw1,ipw2,cse1*27.2114
!endif
!write(6,'("(",f10.6,",",f10.6,")")') cse1*27.2114
!if (lx) write(6,'(f10.6)') dble(xse1)*27.2114
!write(6,*) ioff2-ioff
!write(6,*) cse
!write(6,*) (ssi(ioff2-ioff)+ssi(ioff2-ioff-1))/2.d0
!write(6,*) zz
!write(6,*) (ssi(ioff2-ioff)-ssi(ioff2-ioff-1))/dw

  enddo
!  write(68,'(2i3,2f10.3)') iband,ibp,(sei-seibc)*27.2114,sei*27.2114
enddo

do ie=1,nwpt
  eps1=ek+dble(ie-nwpt/2)*dw+dw/2.d0
!  write(69,'(i4,f10.3,2(3x,2es12.3))') ie,eps1*27.2114,sec(ie)*27.2114
enddo
!cse=(0.d0,1.d0)*sei

zz=1.d0/(1.d0-zz)
!cse=zz*cse
!write(6,*) cse*27.2114
!write(6,*) xse*27.2114
!write(6,*) zz

do ix=1,ngkpt(1)
  write(76,*) ix,'--------------------------------------------------'
  do iy=1,ngkpt(2)
    write(76,'(10es8.1)') dble(temparray(ix,iy,:))
!    write(76,'(10f8.5)') dble(temparray(ix,iy,:))
  enddo
enddo

return
end subroutine mkcse

!**************************************************************************

subroutine fek(xk,bmet,ek)
! Finds the kinetic energy of a plane wave
implicit none
double precision :: xk(3),bmet(3,3)
double precision :: ek
integer ii,jj

ek=0.d0
do ii=1,3
do jj=1,3
  ek=ek+xk(ii)*bmet(ii,jj)*xk(jj)
enddo
enddo
ek=ek/2.d0

return
end subroutine fek

!**************************************************************************

subroutine fcenter(iqq,ngkpt,icenter)
implicit none
integer :: iqq(3),ngkpt(3),icenter
integer :: ii

icenter=1
do ii=1,3
  if (iqq(ii).ne.ngkpt(ii)/2.and.iqq(ii)+1.ne.ngkpt(ii)/2) icenter=0
enddo

return
end subroutine fcenter

!**************************************************************************

subroutine fespread(iqq,ngkpt,omega,esprd)
implicit none
integer :: iqq(3),ngkpt(3)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3)),esprd
integer :: ix,iy,iz,iqp(3)
double precision :: e1,e0

e0=omega(iqq(1),iqq(2),iqq(3))
!write(6,'(3i3,3x,3i3)') iqq
esprd=0.d0
do ix=-1,1
do iy=-1,1
do iz=-1,1
  iqp=mod(iqq+(/ix,iy,iz/)-1+ngkpt,ngkpt)+1
  e1=omega(iqp(1),iqp(2),iqp(3))
  esprd=max(esprd,abs(e1-e0))
!write(6,'(3i3,3x,3i3,3x,2f10.3)') ix,iy,iz,iqp,e1*27.2114,abs(e1-e0)*27.2114
enddo
enddo
enddo

return
end subroutine fespread

!**************************************************************************
! cut material from mkcse for testing energy dependence of self energy

!    do ix=1,ngkpt(1)
!    do iy=1,ngkpt(2)
!    do iz=1,ngkpt(3)
!      iqq=(/ix,iy,iz/)
!      call fhilo(omega,iqq,ngkpt,whi,wlo)
!      if (ix.eq.1.and.iy.eq.1.and.iz.eq.1) then
!        wwhi=whi
!        wwlo=wlo
!      else
!        wwhi=max(wwhi,whi)
!        wwlo=min(wwlo,wlo)
!      endif
!    enddo
!    enddo
!    enddo
!    if (ibp.gt.nbocc) then
!      ioff=nint(-wwlo*dble(nwpt)/wmax)
!    else
!      ioff=nint(wwhi*dble(nwpt)/wmax)-nwpt
!    endif
!    ioff2=nint(ek*dble(nwpt)/wmax)
!    do ie=1,nwpt
!      eps=(ie+ioff-ioff2)*dw+ek+dw/2.d0
!      ssi(ie)=(0.d0,0.d0)
!!IF (.TRUE.) THEN
!IF (.FALSE.) THEN
!      if (ibp.gt.nbocc) then
!        iwh=nint((eps+wwhi)*dble(nwpt)/wmax)
!        iwl=nint((eps+wwlo)*dble(nwpt)/wmax)
!      else
!        iwh=nint((wwhi-eps)*dble(nwpt)/wmax)
!        iwl=nint((wwlo-eps)*dble(nwpt)/wmax)
!      endif
!      do iw=max(iwl,1),min(iwh,nwpt)
!        if (ibp.gt.nbocc) then
!          ww=iw*dw-eps
!        else
!          ww=iw*dw+eps
!        endif
!!        write(6,*) ie,iw,ww*27.2114
!        do ix=1,ngkpt(1)
!        do iy=1,ngkpt(2)
!        do iz=1,ngkpt(3)
!          iqq=(/ix,iy,iz/)
!          iqpt=iqndx(iqq(1),iqq(2),iqq(3))
!          call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
!          if (jjpw.ne.0) then
!            vfactor(iqq(1),iqq(2),iqq(3))=vqmat2(iqq(1),iqq(2),iqq(3))*dimag(W(iw,iqpt,jjpw))
!          else
!            vfactor(iqq(1),iqq(2),iqq(3))=(0.d0,0.d0)
!          endif
!!          stvec(iqq(3))=dble(vfactor(iqq(1),iqq(2),iqq(3)))
!!          stvec2(iqq(3))=dimag(vfactor(iqq(1),iqq(2),iqq(3)))
!!          stvec(iqq(3))=dimag(W(iw,iqpt,jjpw))
!!          stvec(iqq(3))=dble(iqpt)
!!          stvec2(iqq(3))=dble(jjpw)
!        enddo
!!        write(77,'(10es8.1)') stvec(:)
!!        write(77,'(10es8.1)') stvec2(:)
!!        write(77,'(10i8)') nint(stvec(:))
!!        write(77,'(10i8)') nint(stvec2(:))
!!        write(77,*) 
!        enddo
!!        write(77,*) iqq(1)+1,'--------------------------------------------------',ibp
!        enddo
!        seiold=sei
!        do ix=1,ngkpt(1)
!        do iy=1,ngkpt(2)
!        do iz=1,ngkpt(3)
!          iqq=(/ix,iy,iz/)
!          call fval(omega,iqq,iqqp,ngkpt,enval)
!          call fpol(vfactor,ngkpt,iqq,iqqp,vfv)
!!          write(42,*) '>>>>',iqq,ww
!!          write(42,*)
!!          write(42,'(2es12.3,2(4x,2es12.3))') enval(1:2),dble(vfv(1:2)),dimag(vfv(1:2))
!!          write(42,'(2es12.3,2(4x,2es12.3))') enval(3:4),dble(vfv(3:4)),dimag(vfv(3:4))
!!          write(42,*)
!!          write(42,'(2es12.3,2(4x,2es12.3))') enval(5:6),dble(vfv(5:6)),dimag(vfv(5:6))
!!          write(42,'(2es12.3,2(4x,2es12.3))') enval(7:8),dble(vfv(7:8)),dimag(vfv(7:8))
!!          write(42,*)
!          call cubeint(enval,ww,vfv,tseincr)
!!          write(42,*) tseincr
!          call fcenter(iqq,ngkpt,icenter)
!          if (abs(tseincr).ne.0.d0.and.icenter.ne.0 &
!&   .and.(ipw1.eq.1.or.ipw2.eq.1)) then
!            rlr=1.d-2
!            abr=1.d-2
!!            write(42,*) '>>>>',iqq,ww
!!            write(6,*) tseincr
!            call fvq2(iqq,shiftk,shiftkq,bmet,ngkpt,ipw1,ipw2,ipwlf,npwlf,vq2)
!!            write(42,*)
!!            write(42,'(2es12.3,5x,2es12.3)') vq2(1:2),dble(vfv(1:2))
!!            write(42,'(2es12.3,5x,2es12.3)') vq2(3:4),dble(vfv(3:4))
!!            write(42,*)
!!            write(42,'(2es12.3,5x,2es12.3)') vq2(5:6),dble(vfv(5:6))
!!            write(42,'(2es12.3,5x,2es12.3)') vq2(7:8),dble(vfv(7:8))
!!            write(42,*)
!            if (ipw1.ge.ipw2) then
!              lossv=vfv*vq2/(4.d0*pi)
!              call subint2(iqq,shiftk,shiftkq,lossv,enval,ww,iw,bmet,iqndx,ngkpt,nwpt,nqpt,ikndxq,vol,ipw1,ipw2,ipwlf,npwlf,pi,rlr,abr,tseincr)
!            else
!              lossv=-conjg(vfv)*vq2/(4.d0*pi)
!              tseincrx=-conjg(tseincr)
!              call subint2(iqq,shiftk,shiftkq,lossv,enval,ww,iw,bmet,iqndx,ngkpt,nwpt,nqpt,ikndxq,vol,ipw1,ipw2,ipwlf,npwlf,pi,rlr,abr,tseincrx)
!              tseincr=-conjg(tseincrx)
!            endif
!!            write(6,*) tseincr
!          endif
!!          write(42,*) tseincr
!          dse=tseincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
!          ssi(ie)=ssi(ie)+dse*dw*isign/pi
!!          write(66,'(3i4,5x,es10.3)') iqq,dse
!!          write(66,'(3i4,5x,f20.6)') iqq,dble(dse)*27.2114
!!          write(66,'(3i4,a)') iqq,'-----------------------------------------------'
!!          write(66,'(es10.3,5x,2es10.3)') ww,dse
!!          write(66,*)
!!          write(66,'(2es10.3,5x,2es10.3)') dble(vfv(1)),dble(vfv(2)),enval(1:2)-ww
!!          write(66,'(2es10.3,5x,2es10.3)') dble(vfv(3)),dble(vfv(4)),enval(3:4)-ww
!!          write(66,*)
!!          write(66,'(2es10.3,5x,2es10.3)') dble(vfv(5)),dble(vfv(6)),enval(5:6)-ww
!!          write(66,'(2es10.3,5x,2es10.3)') dble(vfv(7)),dble(vfv(8)),enval(7:8)-ww
!!          write(42,*) ssi(iw)
!!          stvec(iqq(3))=dble(tseincr)
!!          stvec2(iqq(3))=dimag(tseincr)
!!          stvec(iqq(3))=dble(icenter)
!        enddo
!!        write(78,'(10es8.1)') stvec(:)
!!        write(78,'(10es8.1)') stvec2(:)
!!        write(78,*) 
!!        write(78,'(10i4)') nint(stvec(:))
!        enddo
!!        write(78,*) iqq(1)+1,'--------------------------------------------------',ibp
!        enddo
!      enddo
!ELSE
!      brd=dw
!      if (ibp.gt.nbocc) then
!        eps1=eps
!      else
!        eps1=-eps
!      endif
!      if (abs(vqmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)/vqmat2(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2-1)).gt.1.d3 &
!&     .and.(ipw1.eq.1.or.ipw2.eq.1)) then
!        iqsing=1
!        iqpt=iqndx(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)
!        call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
!        wint0=(0.d0,0.d0)
!        if (jjpw.ne.0) then
!          do iw=1,nwpt
!            ww=iw*dw
!            eval=omega(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)+eps1-ww
!            wint0=wint0+dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)
!          enddo
!        endif
!!        ssi(ie)=ssi(ie)+wint0*cmgamma2*gfo/pi
!!        ssi(ie)=ssi(ie)+wint0*cmgamma2*4.d0*pi*gfo**3*(1.d0/6.d0-1.d0/(pi**2))
!      else
!        iqsing=0
!      endif
!      do ix=1,ngkpt(1)
!      do iy=1,ngkpt(2)
!      do iz=1,ngkpt(3)
!        iqq=(/ix,iy,iz/)
!        iqpt=iqndx(ix,iy,iz)
!        wint=(0.d0,0.d0)
!        call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
!        if (jjpw.eq.0) cycle
!        do iw=1,nwpt
!          ww=iw*dw
!          eval=omega(ix,iy,iz)+eps1-ww
!          wint=wint+dw*dimag(W(iw,iqpt,jjpw))/(eval,brd)
!        enddo
!        if (iqsing.eq.0) then
!          ssi(ie)=ssi(ie)+wint*vqmat2(ix,iy,iz)
!        else
!          qq2=4.d0*pi/vq(ix,iy,iz)
!          if (sqrt(qq2).lt.gfo) then
!            gterm=cmgamma2*(1+cos(pi*sqrt(qq2)/gfo))/(2.d0)
!            ssi(ie)=ssi(ie)+(wint*vqmat2(ix,iy,iz)/vq(ix,iy,iz)-wint0*gterm)*vq(ix,iy,iz)
!          else
!            gterm=0.d0
!            ssi(ie)=ssi(ie)+wint*vqmat2(ix,iy,iz)
!          endif
!        endif
!      enddo
!      enddo
!      enddo
!      ssi(ie)=ssi(ie)/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
!      if (iqsing.eq.1) then
!! 1/pi Integral (d^3q/(2 pi)^3) wint0 gterm vq
!        ssi(ie)=ssi(ie)+wint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
!      endif
!ENDIF
!    enddo
!
!    do ie=1,nwpt
!!      eps1=ek+dble(ie-nwpt/2)*dw+dw/2.d0
!!      ssc(ie)=(0.d0,0.d0)
!!      do je=1,nwpt
!!        if (ie.eq.(je+ioff-ioff2+nwpt/2)) then
!!          ssc(ie)=ssc(ie)-(0.d0,pi)*ssi(je)
!!        else
!!          eps2=dble(je+ioff-ioff2)*dw+ek+dw/2.d0
!!          ssc(ie)=ssc(ie)+ssi(je)*dw/(eps2-eps1)
!!        endif
!!write(6,'(2i4,2f10.3,es12.3,2(2x,2es11.3))') ie,je,eps1*27.2114,eps2*27.2114,1/(eps2-eps1),ssc(ie)*27.2114,dble(ssi(je)*27.2114)
!!      enddo
!!      write(69,'(i4,f10.3,2(3x,2es12.3))') ie,eps1*27.2114,ssc(ie)*27.2114,ssi(ie)*27.2114
!      sec(ie)=sec(ie)+ssi(ie)
!    enddo
