subroutine mkexch(iband,ikpt, &
& omegap2,xkf,eigen,eigenq,wmax,nwpt,nqpt,ntpwndx,W, &
& nkpt,nkptq,npwt,npwtq,ncg,ncgq,npwarr,npwarrq,igmn,igmx,nbcore,nbocc,ncband, &
& bantot,bantotq,indxkbnd,indxkbndq,iqndx, &
& npwup,invpwndx,pwupmin,pwupmax,iqsymndx,nsym,pwsymndx,invpw2ndx, &
& kg,kgq,cg,cgq,kpt,kptq,indxkcg,indxkcgq,indxkpw,igndxq,ngkpt,ikndxq, &
& exch,ppcor,cor)
use geometry
implicit none
integer :: iband,ikpt
integer :: nkpt,nkptq,npwt,npwtq,ncg,ncgq,nbcore,nbocc,ncband,bantot,bantotq,ntpwndx,nwpt,nqpt,nsym,npwup
double precision :: omegap2,xkf,eigen(bantot),eigenq(bantotq),wmax
double complex :: W(nwpt,nqpt,ntpwndx)
integer :: igmn(3),igmx(3),ngkpt(3)
integer :: ikndxq(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: indxkpw(nkpt),indxkcg(nkpt),indxkcgq(nkptq)
integer :: indxkbnd(nkpt),indxkbndq(nkptq)
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: npwarr(nkpt),npwarrq(nkptq)
integer :: pwupmin(3),pwupmax(3)
integer :: invpwndx(pwupmin(1):pwupmax(1),pwupmin(2):pwupmax(2),pwupmin(3):pwupmax(3))
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
double precision :: kpt(3,nkpt),kptq(3,nkptq)
integer :: kg(3,npwt),kgq(3,npwtq)
double complex :: cg(ncg),cgq(ncgq)
double precision :: exch
double complex :: ppcor(nbcore+1:ncband),cor(nbcore+1:ncband)
double complex :: cse(-nwpt:nwpt)
double precision :: exch1
double complex :: ppcor1(nbcore+1:ncband),cor1(nbcore+1:ncband)
double complex :: cse1(-nwpt:nwpt)
double precision :: pi,fourpi,volelmnt,qq(3),qq2(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ikk(3),iqq(3),ikp1(3),ikp2(3),ikp3(3),jqq(3),iqqp(3),iqv(3)
integer :: ipwgg(3),ishgg(3),igg(3),ictr(3)
integer :: jband,ipw,ix,iy,iz,ii,jj,kk,ie,iw,itet,iv
double precision :: ww,dw,brd
integer :: iqsing,nqgrid,isign,ipwa,ipwb,jjpw,iqpt
integer :: ikqpt
double complex :: cme,vme2(ngkpt(1),ngkpt(2),ngkpt(3)),ame2(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: gamma(3),pln(3),dist(3),falloff,st(6),cterm
double precision :: ekgg,xgg(3)
double precision :: ebkq(ngkpt(1),ngkpt(2),ngkpt(3)),omegat
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: cedenom
double precision :: wwq,ppfct
double complex :: wppint0(nbcore+1:ncband),wpp2int0(nbcore+1:ncband)
double complex :: wint0(nbcore+1:ncband),w2int0(nbcore+1:ncband)
double complex :: wppint(nbcore+1:ncband),wpp2int(nbcore+1:ncband)
double complex :: wint(nbcore+1:ncband),w2int(nbcore+1:ncband)
logical :: lx,lc,lqsing,lqcentr
parameter (pi=3.1415926535897932384626433832795, &
& fourpi=12.566370614359172953850573533118d0)

  ictr=ngkpt/2
  dw=wmax/dble(nwpt)

  gamma=(blat(1,:)+blat(2,:)+blat(3,:))/2.d0
  do ii=1,3
    jj=mod(ii,3)+1
    kk=mod(ii+1,3)+1
    pln(1)=blat(jj,2)*blat(kk,3)-blat(jj,3)*blat(kk,2)
    pln(2)=-blat(jj,1)*blat(kk,3)+blat(jj,3)*blat(kk,1)
    pln(3)=blat(jj,1)*blat(kk,2)-blat(jj,2)*blat(kk,1)
    dist(ii)=dot_product(gamma,pln)/sqrt(dot_product(pln,pln))
  enddo
  falloff=minval(dist)

  nqgrid=ngkpt(1)*ngkpt(2)*ngkpt(3)
  volelmnt=vol*nqgrid
  exch=0.d0
  do ie=nbcore+1,ncband
    ppcor(ie)=(0.d0,0.d0)
    cor(ie)=(0.d0,0.d0)
  enddo
  ikk=nint(kpt(:,ikpt)*dble(ngkpt))
!  do jband=nbcore+1,ncband
!  do jband=nbcore+1,nbocc
  do jband=1,1
!write(6,*) jband
    if (jband.gt.nbocc) then
      isign=1
    else
      isign=-1
    endif
!    do ipw=1,npwarr(ikpt)
    do ipw=1,1
      ipwgg=kg(:,indxkpw(ikpt)+ipw)
      xgg=dble(ipwgg)+dble(ikk)/dble(ngkpt)
      ekgg=0.d0
      do ii=1,3
      do jj=1,3
        ekgg=ekgg+xgg(ii)*bmet(ii,jj)*xgg(jj)
      enddo
      enddo
      if (ekgg.gt.8.0d0) cycle
      if (ipwgg(1).ge.pwupmin(1).and.ipwgg(1).le.pwupmax(1).and.  &
&         ipwgg(2).ge.pwupmin(2).and.ipwgg(2).le.pwupmax(2).and.  &
&         ipwgg(3).ge.pwupmin(3).and.ipwgg(3).le.pwupmax(3)) then
        ipwa=invpwndx(ipwgg(1),ipwgg(2),ipwgg(3))
        ipwb=invpwndx(ipwgg(1),ipwgg(2),ipwgg(3))
      else
        ipwa=0
        ipwb=0
      endif
!write(6,'(i4,5x,3i3)') ipw,ipwgg
!write(6,*) ngkpt
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=2,2
!      do iy=2,2
!      do iz=4,4
! determine k-q point and plane wave shift to remain in B.Z.
        iqq=(/ix,iy,iz/)-ngkpt/2
!write(6,'(3i3)') iqq
        ikp1=ikk-iqq
        ikp2=mod(ikp1+ngkpt/2,ngkpt)
        do ii=1,3
          if (ikp2(ii).le.0) ikp2(ii)=ikp2(ii)+ngkpt(ii)
        enddo
        ikqpt=ikndxq(ikp2(1),ikp2(2),ikp2(3))
        ikp3=nint(kptq(:,ikqpt)*dble(ngkpt))
        ishgg=nint(dble(ikp1-ikp3)/dble(ngkpt))
        igg=ipwgg-ishgg
!write(6,'(a,3i3)') 'ikk   = ',ikk
!write(6,'(a,3i3)') 'iqq   = ',iqq
!write(6,'(a,3i3)') 'ikp1  = ',ikp1
!write(6,'(a,3i3)') 'ikp2  = ',ikp2
!write(6,'(a,3i3)') 'ikp3  = ',ikp3
!write(6,'(a,3i3)') 'ishgg = ',ishgg
!write(6,'(a,3i3)') 'ipwgg = ',ipwgg
!write(6,'(a,3i3)') 'igg   = ',igg
!write(6,'(a,i5)') 'ikpt   = ',ikpt
!write(6,'(a,i5)') 'ikptq   = ',ikqpt
!write(6,'(a,3f10.5)') 'kk (red) = ',kpt(:,ikpt)
!write(6,'(a,3f10.5)') 'kp (red) = ',kptq(:,ikqpt)
! find q^2
        qq=(/0.d0,0.d0,0.d0/)
        do ii=1,3
          qq=qq+(kpt(ii,ikpt)-kptq(ii,ikqpt)+igg(ii))*blat(ii,:)
        enddo
        qq2(ix,iy,iz)=dot_product(qq,qq)
!write(6,'(a,3f10.5)') 'qq (red) = ',kpt(:,ikpt)-kptq(:,ikqpt)+igg
!write(6,'(a,3f10.5)') 'qq = ',qq
!write(6,'(a,f10.5)') 'qq^2 = ',qq2(ix,iy,iz)
! matrix elements
        call fndme(iband,ikpt,jband,ikqpt,igg, &
& nkpt,nkptq,npwt,npwtq,ncg,ncgq,npwarr,npwarrq,igmn,igmx, &
& kg,kgq,cg,cgq,indxkcg,indxkcgq,indxkpw,igndxq, &
& cme)
        ame2(ix,iy,iz)=cme*conjg(cme)
        vme2(ix,iy,iz)=fourpi*ame2(ix,iy,iz)/qq2(ix,iy,iz)
! intermediate state energy
        ebkq(ix,iy,iz)=eigenq(indxkbndq(ikqpt)+jband)
!write(6,'(a,3f10.5)') 'cme = ',cme
!write(6,'(a,3f10.5)') 'ame2 = ',ame2(ix,iy,iz)
!write(6,'(a,3f10.5)') 'vme2 = ',vme2(ix,iy,iz)
!write(6,'(3i3,2f15.10)') ix,iy,iz,cme
      enddo
      enddo
      enddo

      lx=jband.le.nbocc
      lc=ipwa.ne.0
      brd=dw
      if (ipwgg(1).eq.0.and.ipwgg(2).eq.0.and.ipwgg(3).eq.0) then
        ix=ictr(1)
        iy=ictr(2)
        iz=ictr(3)
        st(1)=abs(vme2(ix,iy,iz)/vme2(ix+1,iy,iz))
        st(2)=abs(vme2(ix,iy,iz)/vme2(ix-1,iy,iz))
        st(3)=abs(vme2(ix,iy,iz)/vme2(ix,iy+1,iz))
        st(4)=abs(vme2(ix,iy,iz)/vme2(ix,iy-1,iz))
        st(5)=abs(vme2(ix,iy,iz)/vme2(ix,iy,iz+1))
        st(6)=abs(vme2(ix,iy,iz)/vme2(ix,iy,iz-1))
        if (st(1).gt.1.d3.and.st(2).gt.1.d3.and.st(3).gt.1.d3.and. &
&           st(4).gt.1.d3.and.st(5).gt.1.d3.and.st(6).gt.1.d3) then
          lqsing=.true.
        else
          lqsing=.false.
        endif
      else
        lqsing=.false.
      endif

      call seqptsum(ikpt,lqsing,lx,lc,ipwa,ipwb,ngkpt,iqndx,omegap2,xkf,qq2,W, &
& npwup,nsym,iqsymndx,pwsymndx,nbcore,nbocc,ncband,nwpt,nqpt,ntpwndx, &
& wmax,vme2,ame2,falloff,isign,ebkq,indxkbnd,nkpt,eigen,bantot,invpw2ndx, &
& exch,ppcor,cor)

!      if (lc) then
!        omega=ebkq-eigen(indxkbnd(ikpt)+iband)+dw/2
!        call setetint(ebkq,qq2,vme2,W,omega,isign,lqsing,wmax, &
!& ipwa,ipwb,ngkpt,iqndx,iqsymndx,npwup, &
!& nsym,pwsymndx,invpw2ndx,nwpt,nqpt,ntpwndx, &
!& cse1)
!      endif
!write(6,*) (cse1(0)+cse1(1))/2
!write(6,*) cor(1)

    enddo
!write(6,'(i4,f12.6,4x,2f12.6)') jband,exch*27.2114,cor(iband)*27.2114
  enddo

return
end subroutine mkexch

!****************************************************************************

subroutine fndme(iband,ikpt,jband,ikqpt,igg, &
& nkpt,nkptq,npwt,npwtq,ncg,ncgq,npwarr,npwarrq,igmn,igmx, &
& kg,kgq,cg,cgq,indxkcg,indxkcgq,indxkpw,igndxq, &
& cme)
implicit none
integer :: iband,ikpt,jband,ikqpt,igg(3)
integer :: nkpt,nkptq,npwt,npwtq,ncg,ncgq
integer :: igmn(3),igmx(3)
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: indxkpw(nkpt),indxkcg(nkpt),indxkcgq(nkptq)
integer :: npwarr(nkpt),npwarrq(nkptq)
integer :: kg(3,npwt),kgq(3,npwtq)
double complex :: cg(ncg),cgq(ncgq)
double complex :: cme
integer :: ipw,ipw2,iwf1,iwf2,igg2(3)

!write(6,*) iband,ikpt,jband,ikqpt
  cme=(0.d0,0.d0)
  do ipw=1,npwarr(ikpt)
    iwf1=indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw
    igg2=kg(:,indxkpw(ikpt)+ipw)-igg
    if (igg2(1).ge.igmn(1).and.igg2(1).le.igmx(1).and. &
&       igg2(2).ge.igmn(2).and.igg2(2).le.igmx(2).and. &
&       igg2(3).ge.igmn(3).and.igg2(3).le.igmx(3)) then
      ipw2=igndxq(ikqpt,igg2(1),igg2(2),igg2(3))
      if (ipw2.ne.0) then
        iwf2=indxkcgq(ikqpt)+(jband-1)*npwarrq(ikqpt)+ipw2
        cme=cme+cg(iwf1)*conjg(cgq(iwf2))
!write(6,'(i4,2x,3i3,2x,3i3,2x,i4,2es11.4,2x,2es11.4)') ipw,kg(:,indxkpw(ikpt)+ipw),igg2,ipw2,cg(iwf1),cgq(iwf2)
      endif
    endif
!write(6,'(i4,2x,3i3,2x,3i8)') ipw,kg(:,indxkpw(ikpt)+ipw),indxkcg(ikpt),(iband-1)*npwarr(ikpt),iwf1
  enddo

return
end subroutine fndme





































