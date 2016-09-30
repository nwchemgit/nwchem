! for finite qq
subroutine mksecspec(qq,igglf1,igglf2,iks,vol,pi,ssnrg,ssdnrg, &
& test_bands_pol, &
& nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt, &
& ncg,bantot,npwt,nlmn, &
& nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
& xred, &
& projwf, &
& kg, &
& enrgy, &
& cg, &
& indxkpw, &
& indxkbnd, &
& indxkcg, &
& npwarr, &
& kpt, &
& nsymk,symk,symrel,syminv, &
& ihlf,lvtrans, &
& ntypepaw, &
& pwmatel1,tpwmatel1,pwmatel2,tpwmatel2, &
& igndx,ikndx,isymndx,isymg, &
& nband,secspec)
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),dq(3),vq,vol,pi,wmax,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: igglf1(3),igglf2(3)
integer :: test_bands_pol(4)
double precision :: shiftk(3)
double precision :: xred(3,natom)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,nkpt,nbcore+1:ncband)
double precision :: enrgy(bantot)
double complex :: cg(ncg)
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
integer :: nband(nkpt*nsppol)
double complex :: cmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,cmat2v(8)
double complex :: amatel,amatellf
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3)),e_cond(ngkpt(1),ngkpt(2),ngkpt(3)),e_val(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,whi2,wlo2,ww,ww1,ww2,enval(8),encond(8),wflmin,wflmax
integer :: iwh,iwl,iwh2,iwl2
double precision :: xck(3),xckq(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
integer :: ntypepaw
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double complex :: pwmatel2(ntypepaw,nlmn,nlmn),tpwmatel2(ntypepaw,nlmn,nlmn)
double precision :: ssnrg, ssdnrg, ssnrghi, ssnrglo, secspec(nwpt), specnorm
double complex :: tss


ssnrghi=ssnrg+ssdnrg/2.
ssnrglo=ssnrg-ssdnrg/2.

igg0=(/0,0,0/)
qp=qq-igglf1+igglf2
dw=wmax/dble(nwpt)
do iw=1,nwpt
  secspec(iw)=0.d0
enddo

!write(6,'(a)') 'occupied band, unoccupied band, minimum gap'
!write(6,*) nbocc, ncband-nbocc
do iocc=test_bands_pol(1),test_bands_pol(2)
!do iocc=7,7
  do iunocc=test_bands_pol(3),test_bands_pol(4)
!  do iunocc=16,16
!    write(6,*) iocc,iunocc
    gap=1.d20
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=1,1
!    do iz=1,1
      ikk=(/ix,iy,iz/)
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1
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
      call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&       ncg,nkpt,npwt,igmx,igmn,igndx, &
&       isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&       lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&       cg,indxkcg,indxkpw,npwarr,kg, &
&       cmatel)
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      e_val(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikptq)+iocc)
      e_cond(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc)
      if (ipaw.ne.0) then
        call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
& pwmatel1,tpwmatel1,nbcore,ncband,amatel)
      else
        amatel=(0.d0,0.d0)
      endif
      if (igglf1(1).eq.igglf2(1).and.igglf1(2).eq.igglf2(2).and. &
&         igglf1(3).eq.igglf2(3)) then
        cmat2(ikk(1),ikk(2),ikk(3))=dble((cmatel+amatel)*conjg(cmatel+amatel))
      else
        call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatellf)
        if (ipaw.ne.0) then
          call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qp,igg0, &
& pwmatel2,tpwmatel2,nbcore,ncband,amatellf)
        else
          amatellf=(0.d0,0.d0)
        endif
        cmat2(ikk(1),ikk(2),ikk(3))=(cmatel+amatel)*conjg(cmatellf+amatellf)
      endif
      if (omega(ikk(1),ikk(2),ikk(3)).lt.gap) gap=omega(ikk(1),ikk(2),ikk(3))
    enddo
    enddo
    enddo
!stop
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=8,8
!    do iz=2,2
      ikk=(/ix,iy,iz/)
      call fhilo(e_cond,ikk,ngkpt,whi,wlo)
      call fhilo(e_val,ikk,ngkpt,whi2,wlo2)
!      write(6,'(2i3,2x,3i3,2x,2f14.6,2x,2f14.6)') iocc,iunocc,ix,iy,iz,wlo,whi,wlo2,whi2
      whi2=whi2+ssnrg+ssdnrg
      wlo2=wlo2+ssnrg-ssdnrg
      whi=min(whi,whi2)
      wlo=max(wlo,wlo2)
      if (wlo.lt.whi.and.whi.lt.wmax) then
        iwh=min(nint(whi*dble(nwpt)/wmax),nwpt)
        iwl=max(nint(wlo*dble(nwpt)/wmax),0)
!        write(6,*)
!        write(6,'(2i3,2x,3i3,2x,3i5)') iocc,iunocc,ix,iy,iz,iwl,iwh
!        write(6,*) ssnrg
!        write(6,'(2i3,2x,3i3,2x,2f14.6,2x,2f14.6)') iocc,iunocc,ix,iy,iz,wlo,whi
        if (iwl.gt.nwpt) cycle
        call fval(e_cond,ikk,ikkp,ngkpt,encond)
        call fval(e_val,ikk,ikkp,ngkpt,enval)
        call fpol(cmat2,ngkpt,ikk,ikkp,cmat2v)
        do iw=iwl,iwh
!        do iw=16,16
          ww=iw*dw
!          write(6,*)
!          write(6,'(2i3,2x,3i3,2x,3i5)') iocc,iunocc,ix,iy,iz,iw
!          write(6,*) cmat2v
          call cubeint_d(encond,enval,ww,cmat2v,ww-ssnrghi,ww-ssnrglo,tss)
          secspec(iw)=secspec(iw)+dble(tss)/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
!          if (tss.ne.0.d0) then
!            write(6,'(2i3,2x,3i3,2x,i5,3x,e14.6)') iocc,iunocc,ix,iy,iz,iw,tss
!          endif
        enddo
      endif
    enddo
    enddo
!    write(68,*) ix+1,'------------------------------------------------------'
!    write(68,*)
    enddo
  enddo
enddo

specnorm=0.d0;
do iw=1,nwpt
  specnorm=specnorm+secspec(iw)*dw
!write(6,*) iw,secspec(iw),specnorm
enddo
!write(6,*) specnorm/sqrt(dot_product(qq,qq)*dot_product(qp,qp))
if (specnorm.eq.0.d0) then
  write(6,*) 'ERROR: spectrum normalization is zero'
  stop
endif
do iw=1,nwpt
  secspec(iw)=secspec(iw)/specnorm
!write(6,*) iw,secspec(iw)
enddo

return
end subroutine mksecspec

!*************************************************************************

! intended for use when qq->(0.d0,0.d0,0.d0).  Still use qq to get direction
subroutine mksecspecj(qq,igglf1,igglf2,iks,vol,pi,ssnrg,ssdnrg,bmet, &
& test_bands_pol, &
& nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt, &
& ncg,bantot,npwt,nlmn, &
& nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
& xred, &
& projwf, &
& kg, &
& enrgy, &
& cg, &
& indxkpw, &
& indxkbnd, &
& indxkcg, &
& npwarr, &
& kpt, &
& nsymk,symk,symrel,syminv, &
& ihlf,lvtrans, &
& ntypepaw, &
& pwjmatel,tpwjmatel,pwmatel1,tpwmatel1,pwmatel2,tpwmatel2, &
& igndx,ikndx,isymndx,isymg, &
& nband,secspec)
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),dq(3),vq,vol,pi,wmax,dw
double precision :: bmet(3,3),qqmag,qpmag,qq2,qp2,qqnorm(3),qpnorm(3)
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw
integer :: igglf1(3),igglf2(3)
integer :: test_bands_pol(4)
double precision :: shiftk(3)
double precision :: xred(3,natom)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,nkpt,nbcore+1:ncband)
double precision :: enrgy(bantot)
double complex :: cg(ncg)
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
integer :: nband(nkpt*nsppol)
double complex :: cmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,cmat2v(8)
double complex :: jmatel(3),jmatellf(3),tmat2(ngkpt(1),ngkpt(2),ngkpt(3),3,3)
double complex :: amatel,vmatel(3)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3)),e_cond(ngkpt(1),ngkpt(2),ngkpt(3)),e_val(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,whi2,wlo2,ww,ww1,ww2,enval(8),encond(8),wflmin,wflmax
integer :: iwh,iwl,iwh2,iwl2
double precision :: xck(3),xckq(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
integer :: ntypepaw
double complex :: pwjmatel(ntypepaw,nlmn,nlmn),tpwjmatel(ntypepaw,nlmn,nlmn)
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double complex :: pwmatel2(ntypepaw,nlmn,nlmn),tpwmatel2(ntypepaw,nlmn,nlmn)
double precision :: ssnrg, ssdnrg, ssnrghi, ssnrglo, secspec(nwpt), specnorm
double complex :: tss

ssnrghi=ssnrg+ssdnrg/2.
ssnrglo=ssnrg-ssdnrg/2.

igg0=(/0,0,0/)
qp=qq-igglf1+igglf2
qq2=0.d0
qp2=0.d0
do ii=1,3
do jj=1,3
  qq2 = qq2 + qq(ii)*bmet(ii,jj)*qq(jj)
  qp2 = qp2 + qp(ii)*bmet(ii,jj)*qp(jj)
enddo
enddo
qqmag=sqrt(qq2)
qpmag=sqrt(qp2)
do ii=1,3
  qqnorm(ii)=qq(ii)/qqmag
  qpnorm(ii)=qp(ii)/qpmag
enddo
dw=wmax/dble(nwpt)
do iw=1,nwpt
  secspec(iw)=0.d0
enddo

!write(6,'(a)') 'occupied band, unoccupied band, minimum gap'
!write(6,*) nbocc, ncband-nbocc
do iocc=test_bands_pol(1),test_bands_pol(2)
!do iocc=7,7
  do iunocc=test_bands_pol(3),test_bands_pol(4)
!  do iunocc=16,16
!    write(6,*) iocc,iunocc
    gap=1.d20
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=1,1
!    do iz=1,1
      ikk=(/ix,iy,iz/)
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1
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
!      write(6,'(3f8.4,4x,3f8.4)') kpt(1:3,ikpt),kpt(1:3,ikptq)
!      write(6,'(3i5)') ikpt,ikptq
!      write(6,'(3i5)') isym,isymq,nsym
!      write(6,'(2i4,4x,2i4)') mod(isym-1,nsym)+1,-(2*((isym-1)/nsym)-1),mod(isymq-1,nsym)+1,-(2*((isymq-1)/nsym)-1)
!      do ii=1,3
!        write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!      enddo
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      e_val(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikptq)+iocc)
      e_cond(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc)
      if (igglf1(1).eq.0.and.igglf1(2).eq.0.and.igglf1(3).eq.0) then
        call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         jmatel)
!write(46,'(3("(",f8.4,",",f8.4,")"))') jmatel
        if (ipaw.ne.0) then
          call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&                         pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
          jmatel(:) = jmatel(:) + vmatel(:)
        endif
        do ii=1,3
          jmatel(ii) = jmatel(ii)/omega(ikk(1),ikk(2),ikk(3))
        enddo
      else
        call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatel)
        if (ipaw.ne.0) then
          call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
&                        pwmatel1,tpwmatel1,nbcore,ncband,amatel)
          cmatel = cmatel + amatel
        endif
        jmatel(:) = cmatel*qq(:)/qq2
      endif
      if (igglf1(1).eq.igglf2(1).and.igglf1(2).eq.igglf2(2).and. &
&         igglf1(3).eq.igglf2(3)) then
        do ii=1,3
        do jj=1,3
          tmat2(ikk(1),ikk(2),ikk(3),ii,jj)=dble(jmatel(ii)*conjg(jmatel(jj)))
        enddo
        enddo
      else
        if (igglf2(1).eq.0.and.igglf2(2).eq.0.and.igglf2(3).eq.0) then
          call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&           ncg,nkpt,npwt,igmx,igmn,igndx, &
&           isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&           lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&           cg,indxkcg,indxkpw,npwarr,kg, &
&           jmatellf)
          if (ipaw.ne.0) then
            call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&                           pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
            jmatellf(:) = jmatellf(:) + vmatel(:)
          endif
          do ii=1,3
            jmatellf(ii) = jmatellf(ii)/omega(ikk(1),ikk(2),ikk(3))
          enddo
        else
          call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&           ncg,nkpt,npwt,igmx,igmn,igndx, &
&           isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&           lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&           cg,indxkcg,indxkpw,npwarr,kg, &
&           cmatel)
          if (ipaw.ne.0) then
            call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qp,igg0, &
&                          pwmatel2,tpwmatel2,nbcore,ncband,amatel)
            cmatel = cmatel + amatel
          endif
          jmatellf(:) = cmatel*qp(:)/qp2
        endif
        do ii=1,3
        do jj=1,3
          tmat2(ikk(1),ikk(2),ikk(3),ii,jj)=dble(jmatel(ii)*conjg(jmatellf(jj)))
        enddo
        enddo
      endif
      if (omega(ikk(1),ikk(2),ikk(3)).lt.gap) gap=omega(ikk(1),ikk(2),ikk(3))
!write(6,*) dble(cmat2(ikk(1),ikk(2),ikk(3)))
!write(6,*) dble(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!write(6,*) 
!      stvec(ikk(3))=enrgy(indxkbnd(ikpt)+iunocc)
!      stvec(ikk(3))=e_cond(ikk(1),ikk(2),ikk(3))-dw*58
!      stvec(ikk(3))=dble(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!      stvec2(ikk(3))=dimag(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!      stvec(ikk(3))=dble(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec2(ikk(3))=dimag(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec(ikk(3))=abs(cmatel)
!      stvec(ikk(3))=dble(cmatel)
!      stvec2(ikk(3))=dimag(cmatel)
!      stvec(ikk(3))=dble(cmatellf)
!      stvec2(ikk(3))=dimag(cmatellf)
!      stvec(ikk(3))=dble(cmatel*conjg(cmatellf))*vq
!      stvec(ikk(3))=dble(amatel*conjg(amatellf))*vq
    enddo
!    write(68,'(10f8.3)') stvec(1:ngkpt(3)),stvec2(1:ngkpt(3))
!    write(68,'(10f8.3)') stvec(1:ngkpt(3))
!    write(68,'(10f8.3)') stvec2(1:ngkpt(3))
!    write(68,*)
!    write(68,'(10f8.2)') stvec(1:ngkpt(3))*27.2113845d0
!    write(68,'(10es8.1)') stvec(1:ngkpt(3))
!    write(68,'(10es8.1)') stvec2(1:ngkpt(3))
!    write(68,*)
    enddo
!    write(68,*) ix+1,'------------------------------------------------------'
!    write(68,*) 
    enddo
!    write(6,'(2(i6,9x),f10.3)') iocc,iunocc,gap*27.2113845d0
!write(6,'(6i6))') iocc,iunocc, test_bands_pol
!stop
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=8,8
!    do iz=2,2
      ikk=(/ix,iy,iz/)
      call fhilo(e_cond,ikk,ngkpt,whi,wlo)
      call fhilo(e_val,ikk,ngkpt,whi2,wlo2)
!      write(6,'(2i3,2x,3i3,2x,2f14.6,2x,2f14.6)') iocc,iunocc,ix,iy,iz,wlo,whi,wlo2,whi2
      whi2=whi2+ssnrg+ssdnrg
      wlo2=wlo2+ssnrg-ssdnrg
      whi=min(whi,whi2)
      wlo=max(wlo,wlo2)
      if (wlo.lt.whi.and.whi.lt.wmax) then
        iwh=min(nint(whi*dble(nwpt)/wmax),nwpt)
        iwl=nint(wlo*dble(nwpt)/wmax)
!        write(6,*)
!        write(6,'(2i3,2x,3i3,2x,3i5)') iocc,iunocc,ix,iy,iz,iwl,iwh
!        write(6,*) ssnrg
!        write(6,'(2i3,2x,3i3,2x,2f14.6,2x,2f14.6)') iocc,iunocc,ix,iy,iz,wlo,whi
        if (iwl.gt.nwpt) cycle
        call fval(e_cond,ikk,ikkp,ngkpt,encond)
        call fval(e_val,ikk,ikkp,ngkpt,enval)
!if (iunocc.lt.test_bands_pol(3)) then
!  write(6,*) iunocc,test_bands_pol(3),"A"
!  stop
!endif
        do ii=1,3
!if (iunocc.lt.test_bands_pol(3)) then
!  write(6,*) iunocc,test_bands_pol(3),"B"
!  stop
!endif
        do jj=1,3
!if (iunocc.lt.test_bands_pol(3)) then
!  write(6,*) iunocc,test_bands_pol(3),"C"
!  stop
!endif
          call fpol(tmat2(:,:,:,ii,jj),ngkpt,ikk,ikkp,cmat2v)
          do iw=iwl,iwh
!          do iw=16,16
            ww=iw*dw
!            write(6,*)
!            write(6,'(2i3,2x,3i3,2x,3i5)') iocc,iunocc,ix,iy,iz,iw
!            write(6,*) cmat2v
            call cubeint_d(encond,enval,ww,cmat2v,ww-ssnrghi,ww-ssnrglo,tss)
            secspec(iw)=secspec(iw)+dble(tss)/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
!            if (tss.ne.0.d0) then
!              write(6,'(2i3,2x,3i3,2x,i5,3x,e14.6)') iocc,iunocc,ix,iy,iz,iw,tss
!            endif
          enddo
        enddo
        enddo
      endif
    enddo
    enddo
    enddo
  enddo
enddo

specnorm=0.d0;
do iw=1,nwpt
  specnorm=specnorm+secspec(iw)*dw
enddo
!write(6,*) specnorm/sqrt(dot_product(qq,qq)*dot_product(qp,qp))
if (specnorm.eq.0.d0) then
  write(6,*) 'ERROR: spectrum normalization is zero'
  stop
endif
do iw=1,nwpt
  secspec(iw)=secspec(iw)/specnorm
enddo

return
end subroutine mksecspecj

