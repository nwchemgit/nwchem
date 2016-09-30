subroutine mkppse(iband,jband,ikk,omegap2,xkf,lossfn,fsum, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& pwmatel,tpwmatel, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nqpt,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,bmet,blat,ipaw, &
& ipwlf,npwlf,ipwndx,npwndx,ntpwndx, &
& npwup,invpw2ndx,pwsymndx,iqsymndx,invpwndx,pwupmin,pwupmax, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,iqndx,isymndx,isymndxq,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,zz,zzpp,zzcorr,cse,csepp,csecorr,xse)
implicit none
integer :: iband,jband,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,npwtq,bantot,bantotq,ncg,ncgq,nkpt,nkptq,nqpt,nsym,nsymq,nsppol
integer :: npw,npwq,ipw1,npwlf,npwndx,ntpwndx,ipwa,ipaw
double precision :: vol,pi,wmax,xred(3,natom)
double precision :: omegap2,xkf,ppfct
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
integer :: ipwlf(3,npwlf),ipwndx(2,ntpwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex :: xse1,xse2
double complex :: cse(nbcore+1:ncband),cse1(nbcore+1:ncband),cse2(nbcore+1:ncband)
double complex :: csepp(nbcore+1:ncband),cse1pp(nbcore+1:ncband)
double complex :: csecorr(nbcore+1:ncband),cse1corr(nbcore+1:ncband)
double complex :: zz(nbcore+1:ncband),zz1(nbcore+1:ncband),zz2(nbcore+1:ncband)
double complex :: zzpp(nbcore+1:ncband),zz1pp(nbcore+1:ncband)
double complex :: zzcorr(nbcore+1:ncband),zz1corr(nbcore+1:ncband)
double precision :: xse
integer :: ii,jj,kk,ix,iy,iz,iqq(3),iqqp(3),jka(3),jkb(3),jkk(3),iskip,isign,iloss,insd
integer :: ikpt,ikptq,iks(3),ikslf1(3),ikslf2(3)
integer :: pwupmin(3),pwupmax(3)
integer :: igg(3),igpw(3),igh(3),igg0(3),icenter,isym,isymq,ibp,iqpt,iqsym,iw,ie,je,ies
double precision :: xck(3),xckq(3),qadj(6)
double complex :: cmatel,cmatel2,amatel,amatel2,vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3)),vfv(8),lossv(8)
double precision :: vq2(8)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,enval(8),dw,eshift,wwq
double precision :: abr,rlr
double complex :: tseincr,tseincrx,dse,dcse
integer :: iwh,iwl,ibmin,ibmax
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),qq(3),qp(3),qq2,qp2,qs(3),xk(3),xkmq(3),ek,ekmq,ekpw,xpw(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3))
integer :: npwup,iipw,jjpw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup),invpwndx(pwupmin(1):pwupmax(1),pwupmin(2):pwupmax(2),pwupmin(3):pwupmax(3))
integer :: pwsymndx(npwup,2*nsym)
integer :: ipw2,jw
double precision :: fsum(npwup,nqpt)
double precision :: eps1(nbcore+1:ncband),eps2,eps
double precision :: eval,brd,esprd
double precision :: gfo,gamma(3),pln(3),dist(3)
double complex :: wint(nbcore+1:ncband),wint0(nbcore+1:ncband)
double complex :: w2int(nbcore+1:ncband),w2int0(nbcore+1:ncband)
double complex :: wppint(nbcore+1:ncband),wppint0(nbcore+1:ncband)
double complex :: wpp2int(nbcore+1:ncband),wpp2int0(nbcore+1:ncband)
double complex :: cmgamma2,gterm
double complex :: pwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqsing
integer :: idum
logical :: lx,lc,lcorr,lqsing
double precision :: xdum
double precision :: linterp
external linterp

abr=1.d-6
rlr=1.d-6
sei=0.d0
xse=0.d0
cse=(0.d0,0.d0)
csecorr=(0.d0,0.d0)
zz=(0.d0,0.d0)
igg0=(/0,0,0/)
dw=wmax/dble(nwpt)
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
isym=isymndx(ikk(1),ikk(2),ikk(3))
ek=enrgy(indxkbnd(ikpt)+iband)

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
  xse2=xse
  cse2=cse
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
  do ipw1=1,npwarr(ikpt)
!  do ipw1=1,1
    ipw2=ipw1
    igpw=kg(:,indxkpw(ikpt)+ipw1)
    if (igpw(1).ge.pwupmin(1).and.igpw(1).le.pwupmax(1).and.  &
&       igpw(2).ge.pwupmin(2).and.igpw(2).le.pwupmax(2).and.  &
&       igpw(3).ge.pwupmin(3).and.igpw(3).le.pwupmax(3)) then
      ipwa=invpwndx(igpw(1),igpw(2),igpw(3))
    else
      ipwa=0
    endif
    xpw=dble(igpw)+dble(ikk)/dble(ngkpt)
    ekpw=0.d0
    do ii=1,3
    do jj=1,3
      ekpw=ekpw+xpw(ii)*bmet(ii,jj)*xpw(jj)
    enddo
    enddo
    if (ekpw.gt.8.d0) cycle
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
      qq=qq+igg+igpw
      qq2=0.d0
      do ii=1,3
      do jj=1,3
        qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
      enddo
      enddo
      vq(ix,iy,iz)=4.d0*pi/qq2
      wwq=sqrt(omegap2+(xkf**2/3.d0)*qq2+qq2**2/4.d0)
!write(6,'(a,3f10.5)') 'qq = ',qq
!write(6,'(a,3f10.5)') 'qq2 = ',qq2
!write(6,*) 'vq = ',vq(ix,iy,iz)
!write(6,*) 'jka = ',jka
!write(6,*) 'jkk = ',jkk
!write(6,*) 'jkb = ',jkb
!write(6,*) 'igg = ',igg
!write(6,*) 'ikpt = ',ikpt
!write(6,*) 'ikptq = ',ikptq
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
!      if (.true.) then
        call mkmatelX(ibp,iband,ikpt,ikptq,igpw,igg, &
&         ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&         isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&         cmatel)
        if (ipaw.ne.0) then
!*** need to calculate pwmatel for all plane waves
!*** ipw1 indexes full plane wave set, not the local field plane wave set
!*** pwmatel is calculated for
!          call mkmatelP(pi,xred,natom,ibp,iband,ikpt,ikptq,qq,igg0,ngkpt, &
!&               pwmatel(:,:,ipw1,ix,iy,iz),tpwmatel(:,:,ipw1,ix,iy,iz), &
!&               projwf,nlmn,nkpt,nkptq,ncband,amatel)
        else
          amatel=0.d0
        endif
        if (ibp.gt.nbocc) then
          omega(ix,iy,iz)=-enrgyq(indxkbndq(ikptq)+ibp)
        else
          omega(ix,iy,iz)=enrgyq(indxkbndq(ikptq)+ibp)
        endif
      else
        call mkmatelX(ibp,iband,ikpt,ikptq,igpw,igg, &
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
          amatel=0.d0
        endif
        if (ibp.gt.nbocc) then
          omega(ix,iy,iz)=-enrgy(indxkbnd(ikptq)+ibp)
        else
          omega(ix,iy,iz)=enrgy(indxkbnd(ikptq)+ibp)
        endif
      endif
      if (iband.eq.jband) then
        cmatel2=cmatel
        amatel2=amatel
      elseif (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
!      if (.true.) then
        call mkmatelX(ibp,jband,ikpt,ikptq,igpw,igg, &
&         ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&         isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&         cmatel2)
        if (ipaw.ne.0) then
!*** need to calculate pwmatel for all plane waves
!*** ipw1 indexes full plane wave set, not the local field plane wave set
!*** pwmatel is calculated for
!          call mkmatelP(pi,xred,natom,ibp,jband,ikpt,ikptq,qq,igg0,ngkpt, &
!&               pwmatel(:,:,ipw1,ix,iy,iz),tpwmatel(:,:,ipw1,ix,iy,iz), &
!&               projwf,nlmn,nkpt,nkptq,ncband,amatel2)
        else
          amatel2=0.d0
        endif
      else
        call mkmatelX(ibp,jband,ikpt,ikptq,igpw,igg, &
&         ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&         isym,isymq,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&         cmatel2)
        if (ipaw.ne.0) then
!          call mkmatelP(pi,xred,natom,ibp,jband,ikpt,ikptq,qq,igg0,ngkpt, &
!&               pwmatel(:,:,ipw1,ix,iy,iz),tpwmatel(:,:,ipw1,ix,iy,iz), &
!&               projwf,nlmn,nkpt,nkpt,ncband,amatel2)
        else
          amatel2=0.d0
        endif
      endif
      vqmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)*vq(ix,iy,iz)
      if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) then
        cmgamma2=(cmatel+amatel)*conjg(cmatel+amatel)
      endif
!write(6,*) 'cmatel  = ',cmatel
!write(6,*) 'vqmat2 = ',vqmat2(iqq(1),iqq(2),iqq(3))
!write(56,*) '>>>>>',iqq,vqmat2(iqq(1),iqq(2),iqq(3))
!      stvec(iqq(3))=vq(ix,iy,iz)
!      stvec(iqq(3))=sqrt(qq2)
!      stvec(iqq(3))=dble(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec2(iqq(3))=dimag(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec(iqq(3))=dble(cmatel)
!      stvec2(iqq(3))=dimag(cmatel)
!      stvec(iqq(3))=dble(cmatel*conjg(cmatel))
    enddo
!    write(76,'(10f8.2)') stvec(:)
!    write(76,'(10es8.1)') stvec(:)
!    write(76,'(10es8.1)') stvec2(:)
!    write(76,*) 
    enddo
!    write(76,*) iqq(1)+1,'--------------------------------------------------',ibp
    enddo
!stop
!return
! Cut material for testing energy dependence of self energy at bottom of file

    lx=ibp.le.nbocc
    lcorr=ipwa.ne.0
    if (lx) xse1=(0.d0,0.d0)
    do ie=nbcore+1,ncband
      cse1(ie)=(0.d0,0.d0)
      if (lcorr) cse1corr(ie)=(0.d0,0.d0)
      zz1(ie)=(0.d0,0.d0)
      eps=enrgy(indxkbnd(ikpt)+ie)
      if (ibp.gt.nbocc) then
        eps1(ie)=eps
      else
        eps1(ie)=-eps
      endif
    enddo
    brd=dw


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
        iqq=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
        iqpt=iqndx(iqq(1),iqq(2),iqq(3))
!        call fespread(iqq,ngkpt,omega,esprd)
!        brd=max(esprd,dw)
        qq2=4.d0*pi/vq(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)
        wwq=sqrt(omegap2+(xkf**2/3.d0)*qq2+qq2**2/4.d0)
        ppfct=pi*omegap2/(2.d0*wwq)
        if (lcorr) then
          call locateelement(iqq,ipwa,ipwa,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
        else
          jjpw=0
        endif
        do ie=nbcore+1,ncband
          eval=omega(iqq(1),iqq(2),iqq(3))+eps1(ie)-wwq
          wppint0(ie)=dble(ppfct/(eval,brd))
          wpp2int0(ie)=-dble(ppfct/(eval,brd)**2)
          wint0(ie)=(0.d0,0.d0)
          w2int0(ie)=(0.d0,0.d0)
          if (jjpw.ne.0) then
            do iw=1,nwpt
              ww=iw*dw
              eval=omega(iqq(1),iqq(2),iqq(3))+eps1(ie)-ww
              wint0(ie)=wint0(ie)+dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)
              w2int0(ie)=w2int0(ie)-dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)**2
            enddo
          endif
        enddo
      else
        lqsing=.false.
      endif
    else
      lqsing=.false.
    endif

    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=1,1
!    do iz=1,1
      iqq=(/ix,iy,iz/)
      iqpt=iqndx(ix,iy,iz)
!      call fespread(iqq,ngkpt,omega,esprd)
!      brd=3*max(esprd,dw)
!write(6,'(3i3,es20.10)') ix,iy,iz,esprd*27.2114
      qq2=4.d0*pi/vq(ix,iy,iz)
      wwq=sqrt(omegap2+(xkf**2/3.d0)*qq2+qq2**2/4.d0)
      ppfct=pi*omegap2/(2.d0*wwq)
      if (lcorr) then
        call locateelement(iqq,ipwa,ipwa,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
      else
        jjpw=0
      endif
      do ie=nbcore+1,ncband
        eval=omega(ix,iy,iz)+eps1(ie)-wwq
        wppint(ie)=ppfct/(eval,brd)
        wpp2int(ie)=-ppfct/(eval,brd)**2
        wint(ie)=(0.d0,0.d0)
        w2int(ie)=(0.d0,0.d0)
        if (jjpw.ne.0) then
          do iw=1,nwpt
            ww=iw*dw
            eval=omega(ix,iy,iz)+eps1(ie)-ww
            wint(ie)=wint(ie)+dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)
            w2int(ie)=w2int(ie)-dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)**2
          enddo
        endif
      enddo
      if (.not.iqsing) then
        if (lx) xse1=xse1+vqmat2(ix,iy,iz)
        cse1pp=cse1pp+wppint*vqmat2(ix,iy,iz)
        zz1pp=zz1pp+wpp2int*vqmat2(ix,iy,iz)
        if (lcorr) then
          cse1=cse1+wint*vqmat2(ix,iy,iz)
          zz1=zz1+w2int*vqmat2(ix,iy,iz)
          cse1corr=cse1corr &
&                +wppint*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt))
          zz1corr=zz1corr &
&                +wpp2int*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt))
        endif
      else
        if (sqrt(qq2).lt.gfo) then
          if (ix.eq.ngkpt(1)/2.and.iy.eq.ngkpt(2)/2.and.iz.eq.ngkpt(3)/2) cycle
          gterm=cmgamma2*(1+cos(pi*sqrt(qq2)/gfo))/(2.d0)
          if (lx) xse1=xse1+vqmat2(ix,iy,iz)-gterm*vq(ix,iy,iz)
          cse1pp=cse1pp+wppint*vqmat2(ix,iy,iz)-wppint0*gterm*vq(ix,iy,iz)
          zz1pp=zz1pp+wpp2int*vqmat2(ix,iy,iz)-wpp2int0*gterm*vq(ix,iy,iz)
          if (lcorr) then
            cse1=cse1+wint*vqmat2(ix,iy,iz)-wint0*gterm*vq(ix,iy,iz)
            zz1=zz1+w2int*vqmat2(ix,iy,iz)-w2int0*gterm*vq(ix,iy,iz) 
            cse1corr=cse1corr &
&                 +wppint*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt)) &
&                 -wppint0*gterm*vq(ix,iy,iz)*(1.d0-fsum(ipwa,1))
            zz1corr=zz1corr &
&                 +wpp2int*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt)) &
&                 -wpp2int0*gterm*vq(ix,iy,iz)*(1.d0-fsum(ipwa,1))
          endif
        else
          if (lx) xse1=xse1+vqmat2(ix,iy,iz)
          cse1pp=cse1pp+wppint*vqmat2(ix,iy,iz)
          zz1pp=zz1pp+wpp2int*vqmat2(ix,iy,iz)
          if (lcorr) then
            cse1=cse1+wint*vqmat2(ix,iy,iz)
            zz1=zz1+w2int*vqmat2(ix,iy,iz)
            cse1corr=cse1corr &
&                 +wppint*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt))
            zz1corr=zz1corr &
&                 +wpp2int*vqmat2(ix,iy,iz)*(1.d0-fsum(ipwa,iqpt))
          endif
        endif
      endif
    enddo
    enddo
    enddo
    if (lx) xse1=xse1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
    cse1pp=cse1pp/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    zz1pp=zz1pp/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    if (lcorr) then
      cse1=cse1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      zz1=zz1/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      cse1corr=cse1corr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      zz1pp=zz1pp/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
    endif
    if (iqsing) then
! 1/pi Integral (d^3q/(2 pi)^3) wint0 gterm vq
      if (lx) xse1=xse1+cmgamma2*8.d0*pi*pi*gfo/((2.d0*pi)**3)
      cse1pp=cse1pp+wppint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
      zz1pp=zz1pp+wpp2int0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
      if (lcorr) then
        cse1=cse1+wint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
        zz1=zz1+w2int0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)
        cse1corr=cse1corr+ &
&         wppint0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)*(1.d0-fsum(ipwa,1))
        zz1corr=zz1corr+ &
&         wpp2int0*cmgamma2*8.d0*pi*pi*gfo/(pi*(2.d0*pi)**3)*(1.d0-fsum(ipwa,1))
      endif
    endif
    if (lx) xse1=-xse1
    cse1pp=-cse1pp*isign
    if (lcorr) then
      cse1=-cse1*isign
      zz1=-zz1*isign
      cse1corr=-cse1corr*isign
      zz1corr=-zz1corr*isign
    endif
    if (lx) xse=xse+dble(xse1)
    csepp=csepp+cse1pp
    if (lcorr) then
      cse=cse+cse1
      zz=zz+zz1
      csecorr=csecorr+cse1corr
      zzcorr=zzcorr+zz1corr
    endif

    omega=-isign*omega-ek+dw/2
    lx=.false.
    ictr=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
      iqq=(/ix,iy,iz/)
      call fval(omega,iqq,iqqp,ngkpt,enval)
      if (lqsing) call fval(vq,iqq,iqqp,ngkpt,vq2)
      do iw=1,nwpt
        call fpol(vfactor(:,:,:,iw),ngkpt,iqq,iqqp,vfv(:,iw))
      enddo
      if (lx) call fpol(vqmat2(:,:,:),ngkpt,iqq,iqqp,vfx)
      do itet=1,6
        lqcentr=.false.
        do iv=1,4
          evtx0(iv)=enval(ivndx(iv,itet))
          rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
          do iw=1,nwpt
            vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
          enddo
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
          do iv=1,4
            vqvtx0(iv)=vq2(ivndx(iv,itet))
          enddo
          bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
          do jj=1,3
            bgrad=bgrad+(evtx(jj+1)-evtx(1))*rg(:,jj)
          enddo
          xmult=4.d0*pi/sqrt(dot_product(bgrad,bgrad))
          do iw=1,nwpt
            vpyr(:)=vpyr0(:,iw)/vqvtx0
            aa0(iw)=vpyr(indxe(1))
            av(1:3,iw)=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              av(1:3,iw)=av(1:3,iw) &
&                       +(vpyr(indxe(jj+1))-vpyr(indxe(1)))*rg(1:3,jj)
            enddo
          enddo
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
                call singint2(www,xk,kvtx(:,1),evtx,sint1a,sveca)
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
            do iw=1,nwpt
              ie=iww-isign*iw
              ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
            enddo
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
          do iw=1,nwpt
            aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
            av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
            do jj=1,3
              av(1:3,iw)=av(1:3,iw) &
&                       +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
            enddo
          enddo
          if (lx) then
            xme=xpyr0(indxe(1))
            xv=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              xv=xv+(xpyr0(indxe(jj+1))-xpyr0(indxe(1)))*rg(:,jj)
            enddo
          endif
          do iww=iwwlo,iwwhi
!          do iww=2,2
!write(6,*) iww
            www=dw*iww
            if (www.lt.evtx(1)) then
              cycle
            elseif (www.lt.evtx(2)) then
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
            do iw=1,nwpt
              ie=iww-isign*iw
              ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
            enddo
            if (lx) ssx=ssx+(xme*sint1+dot_product(xv,svec))*dw
          enddo
        endif
      enddo
    enddo
    enddo
    enddo
    ssi=-isign*ssi/((2*pi)**3*pi)
    if (lx) ssx=-ssx/((2*pi)**3)




!if (lx) then 
!  write(33,'(i4,4x,i5,3x,3i3,5x,"(",f10.6,",",f10.6,")",5x,f10.6)') &
!&        ibp,ipw1,kg(ii,indxkpw(ikpt)+ipw1),cse1*27.2114,dble(xse1)*27.2114
!else
!  write(33,'(i4,4x,i5,3x,3i3,5x,"(",f10.6,",",f10.6,")")') &
!&        ibp,ipw1,kg(ii,indxkpw(ikpt)+ipw1),cse1*27.2114
!endif
  enddo
!  if (lx) then
!    write(6,'(i4,5x,"(",f10.6,",",f10.6,")",5x,f10.6)') ibp,(cse(4)-cse2(4))*27.2114,dble(xse-xse2)*27.2114
!  else
!    write(6,'(i4,5x,"(",f10.6,",",f10.6,")",f10.6)') ibp,(cse(4)-cse2(4))*27.2114
!  endif
enddo

!zz=1.d0/(1.d0-zz)
!zzpp=1.d0/(1.d0-zzpp)
!zzcorr=1.d0/(1.d0-zz-zzcorr)
!write(6,*) cse*27.2114
!write(6,*) csecorr*27.2114
!write(6,*) xse*27.2114
!write(6,*) zz


return
end subroutine mkppse

!**************************************************************************

subroutine tebound(iqq,ngkpt,omega,insd)
implicit none
integer :: iqq(3),ngkpt(3),insd
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ix,iy,iz,iqp(3)
double precision :: e1,e0

insd=0
e0=omega(iqq(1),iqq(2),iqq(3))
do ix=0,1
do iy=0,1
do iz=0,1
  iqp=mod(iqq+(/ix,iy,iz/)-1+ngkpt,ngkpt)+1
  e1=omega(iqp(1),iqp(2),iqp(3))
  if (e0*e1.lt.0.d0) insd=1
enddo
enddo
enddo

return
end subroutine tebound

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
!            vfactor(iqq(1),iqq(2),iqq(3))=vqmat2(iqq(1),iqq(2),iqq(3))*dimag(lossfn(iw,iqpt,jjpw))
!          else
!            vfactor(iqq(1),iqq(2),iqq(3))=(0.d0,0.d0)
!          endif
!!          stvec(iqq(3))=dble(vfactor(iqq(1),iqq(2),iqq(3)))
!!          stvec2(iqq(3))=dimag(vfactor(iqq(1),iqq(2),iqq(3)))
!!          stvec(iqq(3))=dimag(lossfn(iw,iqpt,jjpw))
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
!            wint0=wint0+dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)
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
!          wint=wint+dw*dimag(lossfn(iw,iqpt,jjpw))/(eval,brd)
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
