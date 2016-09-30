subroutine elloss(mass,igb,ikk,igglf,W, &
& vol,pi,nwpt,wmax,nbocc,ncband,ngkpt, &
& kg,kgq,enrgy,enrgyq,cg,cgq,npwt,npwtq,bantot,bantotq,ncg,ncgq, &
& indxkpw,indxkpwq,indxkbnd,indxkbndq,indxkcg,indxkcgq,npwarr,npwarrq, &
& kpt,kptq,nkpt,nkptq,nqpt,nsymk,nsymkq,symk,symkq,nsym,nsymq,symrel,syminv, &
& ihlf,ihlfq,lvtrans,lvtransq,bmet,ipwlf,npwlf,ipwndx,npwndx,ntpwndx, &
& npwup,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,igndxq,ikndx,ikndxq,iqndx,isymndx,isymndxq,npw,npwq, &
& nband,nbandq,nsppol,shiftk,shiftkq,cse,strate)
implicit none
integer :: igb(3),ikk(3),igglf(3),nwpt,nbocc,ncband,ngkpt(3)
integer :: igmn(3),igmx(3)
integer :: npwt,npwtq,bantot,bantotq,ncg,ncgq,nkpt,nkptq,nqpt,nsym,nsymq,nsppol
integer :: npw,npwq,ipw1,npwlf,npwndx,ntpwndx
double precision :: mass,vol,pi,wmax
double complex :: W(nwpt,nqpt,ntpwndx)
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
double precision :: bmet(3,3)
integer :: ipwlf(3,npwlf),ipwndx(2,ntpwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nband(nkpt*nsppol),nbandq(nkptq*nsppol)
double complex :: sei,seipw
double precision :: strate
double complex :: ssi(nwpt),ssc(nwpt),cse
integer :: ii,jj,kk,ix,iy,iz,iqq(3),iqqp(3),jka(3),jkk(3),iskip
integer :: ikpt,ikptq,iks(3),ikslf1(3),ikslf2(3)
integer :: igg(3),igh(3),igg0(3),icenter,isym,isymq,iqpt,iqsym,iw
double precision :: xck(3),xckq(3)
double complex :: cmatel,cmatel2,vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3)),vfv(8),lossv(8)
double precision :: vq2(8)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,ww,enval(8),dw
double precision :: abr,rlr
double complex :: tseincr,tseincrx,dse
integer :: iwh,iwl,ibmin,ibmax
double precision :: vq,qq(3),qp(3),qq2,qp2,qs(3),xk(3),xkmq(3),ek,ekmq
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),seiold,seibc
integer :: npwup,iipw,jjpw
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
integer :: ipw2,jw
double precision :: ww1,ww2
integer :: idum
double precision :: xdum

abr=1.d-6
rlr=1.d-6
sei=0.d0
strate=0.d0
cse=(0.d0,0.d0)
igg0=(/0,0,0/)
do iw=1,nwpt
  ssi(iw)=(0.d0,0.d0)
  ssc(iw)=(0.d0,0.d0)
enddo
dw=wmax/dble(nwpt)
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
isym=isymndx(ikk(1),ikk(2),ikk(3))
do ii=1,3
  xk(ii)=igb(ii)+kpt(ii,ikpt)
enddo
call fek(xk,bmet,ek)
ek=ek/mass
  seibc=sei
  iskip=1
  do iipw=1,ntpwndx
!  do iipw=1,1
!    write(6,'(a,i4)') 'iipw ',iipw
    seipw=sei
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
    if (ipw1.ne.ipw2) cycle
!    xdum=1.d26
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
      iqq=(/ix,iy,iz/)
      do ii=1,3
        iks(ii)=iqq(ii)-ngkpt(ii)/2
        ikslf1(ii)=iks(ii)+ipwlf(ii,ipw1)*ngkpt(ii)
      enddo
      if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
        call fqq(shiftk,shiftkq,ikslf1,ngkpt,qq,qs)
      else
        call fqq(shiftk,shiftk,ikslf1,ngkpt,qq,qs)
      endif
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
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))
      enddo
      ikptq=ikndxq(jkk(1),jkk(2),jkk(3))
      isymq=isymndxq(jkk(1),jkk(2),jkk(3))
      do ii=1,3
        xkmq(ii)=xk(ii)-qq(ii)
      enddo
      call fek(xkmq,bmet,ekmq)
      ekmq=ekmq/mass
      omega(iqq(1),iqq(2),iqq(3))=ek-ekmq
      if (omega(iqq(1),iqq(2),iqq(3)).gt.0.d0) iskip=0
    enddo
!    write(75,'(10f8.2)') stvec(:)*27.2113845d0
!    write(75,'(10f8.4)') stvec(:)
    enddo
!    write(75,*) iqq(1)+1,'--------------------------------------------------'
    enddo
    if (iskip.eq.1) cycle
    write(6,'(38x,i3,14x,2i3)') ipw1
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=2,2
!    do iy=2,2
!    do iz=4,5
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
      vq=4.d0*pi/sqrt(qq2*qp2)
!write(6,*) '>>>>> iqq = ',iqq
!write(6,*) 'ikslf1 = ',ikslf1
!write(6,*) 'ikslf2 = ',ikslf2
!write(6,*) 'qq = ',qq
!write(6,*) 'qp = ',qp
!write(6,*) 'qq2 = ',qq2
!write(6,*) 'qp2 = ',qp2
!write(6,*) 'vq = ',vq
      vqmat2(iqq(1),iqq(2),iqq(3))=vq
!write(6,*) 'cmatel  = ',cmatel
!write(6,*) 'cmatel2 = ',cmatel2
!write(6,*) 'vqmat2 = ',vqmat2(iqq(1),iqq(2),iqq(3))
!write(56,*) '>>>>>',iqq,vqmat2(iqq(1),iqq(2),iqq(3))
!      stvec(iqq(3))=vq
!      stvec(iqq(3))=dble(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec2(iqq(3))=dimag(vqmat2(iqq(1),iqq(2),iqq(3)))
!      stvec(iqq(3))=dble(cmatel)
!      stvec2(iqq(3))=dimag(cmatel)
!      stvec(iqq(3))=dble(cmatel2)
!      stvec2(iqq(3))=dimag(cmatel2)
!      stvec(iqq(3))=dble(cmatel*conjg(cmatel))
    enddo
!    write(76,'(10f8.2)') stvec(:)
!    write(76,'(10es8.1)') stvec(:)
!    write(76,'(10es8.1)') stvec2(:)
!    write(76,*) 
    enddo
!    write(76,*) iqq(1)+1,'--------------------------------------------------'
    enddo
!    stop
    do iw=1,nwpt
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=1,1
!      do iy=1,1
!      do iz=3,3
        iqq=(/ix,iy,iz/)
        iqpt=iqndx(iqq(1),iqq(2),iqq(3))
        call locateelement(iqq,ipw1,ipw2,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
        if (jjpw.ne.0) then
          vfactor(iqq(1),iqq(2),iqq(3))=vqmat2(iqq(1),iqq(2),iqq(3))*dimag(W(iw,iqpt,jjpw))
        else
          vfactor(iqq(1),iqq(2),iqq(3))=(0.d0,0.d0)
        endif
!        stvec(iqq(3))=dble(vfactor(iqq(1),iqq(2),iqq(3)))
!        stvec2(iqq(3))=dimag(vfactor(iqq(1),iqq(2),iqq(3)))
!        stvec(iqq(3))=dimag(W(iw,iqpt,jjpw))
!        stvec(iqq(3))=dble(iqpt)
!        stvec2(iqq(3))=dble(jjpw)
      enddo
!      write(77,'(10es8.1)') stvec(:)
!      write(77,'(10es8.1)') stvec2(:)
!      write(77,'(10i8)') nint(stvec(:))
!      write(77,'(10i8)') nint(stvec2(:))
!      write(77,*) 
      enddo
!      write(77,*) iqq(1)+1,'--------------------------------------------------'
      enddo
      seiold=sei
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=3,3
!      do iy=1,1
!      do iz=10,10
        iqq=(/ix,iy,iz/)
        call fhilo(omega,iqq,ngkpt,whi,wlo)
        iwh=min(nint(whi*dble(nwpt)/wmax),nwpt)
        iwl=nint(wlo*dble(nwpt)/wmax)
        if (iwl.gt.nwpt) cycle
        call fval(omega,iqq,iqqp,ngkpt,enval)
        call fpol(vfactor,ngkpt,iqq,iqqp,vfv)
        ww=iw*dw
!        write(42,*) '>>>>',iqq,ww
!        write(42,*)
!        write(42,'(2es12.3,2(4x,2es12.3))') enval(1:2),dble(vfv(1:2)),dimag(vfv(1:2))
!        write(42,'(2es12.3,2(4x,2es12.3))') enval(3:4),dble(vfv(3:4)),dimag(vfv(3:4))
!        write(42,*)
!        write(42,'(2es12.3,2(4x,2es12.3))') enval(5:6),dble(vfv(5:6)),dimag(vfv(5:6))
!        write(42,'(2es12.3,2(4x,2es12.3))') enval(7:8),dble(vfv(7:8)),dimag(vfv(7:8))
!        write(42,*)
        call cubeint(enval,ww,vfv,tseincr)
!        write(42,*) tseincr
        call fcenter(iqq,ngkpt,icenter)
        if (abs(tseincr).ne.0.d0) then
          rlr=1.d-2
          abr=1.d-2
!          write(6,*) ix,iy,iz
!          write(6,*) tseincr
!          call floss(W,nwpt,nqpt,iw,ngkpt,iqndx,iqq,lossv)
          call fvq2(iqq,shiftk,shiftkq,bmet,ngkpt,ipw1,ipw2,ipwlf,npwlf,vq2)
!          write(42,*)
!          write(42,'(2es12.3,2x,2(3x,2es12.3))') vq2(1:2),dble(vfv(1:2)),dimag(vfv(1:2))
!          write(42,'(2es12.3,2x,2(3x,2es12.3))') vq2(3:4),dble(vfv(3:4)),dimag(vfv(3:4))
!          write(42,*)
!          write(42,'(2es12.3,2x,2(3x,2es12.3))') vq2(5:6),dble(vfv(5:6)),dimag(vfv(5:6))
!          write(42,'(2es12.3,2x,2(3x,2es12.3))') vq2(7:8),dble(vfv(7:8)),dimag(vfv(7:8))
!          write(42,*)
          lossv=vfv*vq2/(4.d0*pi)
          call subint(mass,ikk,igb,xk,iqq,shiftk,shiftkq,lossv,ww,iw,bmet,iqndx,ngkpt,nwpt,nqpt,ek,ikndxq,vol,ipw1,ipwlf,npwlf,pi,rlr,abr,tseincr)
!          write(6,*) tseincr
        endif
!        write(42,*) tseincr
        dse=tseincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
        sei=sei+dse*dw/pi
        strate=strate+2.d0*dse*ww*dw
        ssi(iw)=ssi(iw)+dse
!        write(66,'(3i4,5x,es10.3)') iqq,dse
!        write(66,'(3i4,5x,f20.6)') iqq,dble(dse)*27.2114
!        write(66,'(3i4,a)') iqq,'-----------------------------------------------'
!        write(66,'(es10.3,5x,2es10.3)') ww,dse
!        write(66,*)
!        write(66,'(2es10.3,5x,2es10.3)') dble(vfv(1:2)),enval(1:2)-ww
!        write(66,'(2es10.3,5x,2es10.3)') dble(vfv(3:4)),enval(3:4)-ww
!        write(66,*)
!        write(66,'(2es10.3,5x,2es10.3)') dble(vfv(5:6)),enval(5:6)-ww
!        write(66,'(2es10.3,5x,2es10.3)') dble(vfv(7:8)),enval(7:8)-ww
!        write(42,*) ssi(iw)
!        stvec(iqq(3))=dble(tseincr)
!        stvec2(iqq(3))=dimag(tseincr)
!        stvec(iqq(3))=dble(icenter)
      enddo
!      write(78,'(10es8.1)') stvec(:)
!      write(78,'(10es8.1)') stvec2(:)
!      write(78,*) 
!      write(78,'(10i4)') nint(stvec(:))
      enddo
!      write(78,*) iqq(1)+1,'--------------------------------------------------'
      enddo
!      write(67,'(i4,f10.3,5x,2es10.3,5x,2es10.3)') iw,ww*27.2114,(sei-seiold)*27.2114,sei*27.2114
    enddo
    write(6,*) iipw,(sei-seipw)*pi*27.2114
  enddo
!  write(68,'(2i3,2f10.3)') (sei-seibc)*27.2114,sei*27.2114

do iw=1,nwpt
  ww1=dble(iw)*dw
  do jw=1,nwpt
    if (iw.eq.jw) cycle
    ww2=dble(jw)*dw
    ssc(iw)=ssc(iw)+2.d0*ssi(jw)*dw*(ww2/(ww2**2-ww1**2))
  enddo
  ssc(iw)=ssc(iw)-(0.d0,pi)*ssi(iw)
  cse=cse+ssc(iw)*dw/pi
!  write(84,'(3es12.4)') ww1,ssi(iw)
enddo

return
end subroutine elloss
