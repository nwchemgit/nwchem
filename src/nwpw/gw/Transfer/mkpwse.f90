subroutine mkpwse(igg0,ikk,mode,lossfn, &
& pwse_mass, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& test_bands_se, &
& ntypepaw, &
& pwmatel,tpwmatel, pwjmatel, tpwjmatel, &
& kg,enrgy,cg,npwt,bantot,ncg, &
& indxkpw,indxkbnd,indxkcg,npwarr, &
& kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
& ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
& ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
& npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
& nband,nsppol,shiftk,zz,dcse,cse,xse)
implicit none
integer :: mode
! mode = 0 for loss rate
! mode = 1 for power loss
integer :: ipwv,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,bantot,ncg,nkpt,nqpt,nsym,nsppol
integer :: npw,npwc,npwx,ipw1,npwndx,ntpwndx,napwndx,ipaw,itetrahedron
double precision :: vol,pi,wmax,xred(3,natom)
double precision :: pwse_mass
double complex :: lossfn(nwpt,nqpt+9,ntpwndx)
integer :: test_bands_se(2)
double complex :: projwf(natom,nlmn,nkpt,ncband)
integer :: kg(3,npwt)
double precision :: enrgy(bantot)
double complex :: cg(ncg)
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt),shiftk(3)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
double precision :: bmet(3,3),blat(3,3),bmetinv(3,3),bmet_t(3,3),bmet2(3,3)
double precision :: volelmnt, volred
integer :: ipwx(3,npwx),ipwndx(2,napwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
integer :: nband(nkpt*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex, allocatable :: ssi(:)
double complex :: ssc(-nwpt:nwpt),cse,dcse,cse0,dcse0,cse1(-nwpt:nwpt),zz,ctest,cse2(-nwpt:nwpt)
double precision :: xse,xse0,ssx,ssx1,ssx2,xse2,ssx3
integer :: ii,jj,kk,ll,ix,iy,iz,iskip,ocsign,iiq
integer :: iqq(3),jka(3),jkb(3),jkk(3),iqv(3),ictr(3),iqr(3),ixr(3)
double precision :: qr(3,4),qqr(3,4)
integer :: ixx(3),ixxp(3),ixp2(3),ixp1(3),ixm1(3),ixm2(3)
integer :: ikpt,ikptq
integer :: igg(3),igh(3),igg0(3),igg1(3),igg2(3),idgg1(3),idgg2(3),isgg1(3),isgg(2),igq
integer :: isym,isymq,iqpt,iqsym,iw,iww,ie,je,ies,iqctr,ibp
double precision :: xck(3),xckq(3),qadj(6)
double complex :: cmatel,cmatel2,amatel,amatel2
double complex :: jmatel(3),jmatel2(3),vmatel(3),vmatel2(3)
double complex :: vqmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: xmat2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: jmat2(3,3), rjmat2(3,3), rjmatdum
double precision :: smat2
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3),nwpt),vjfactor(3,3,nwpt),vsfactor(3,3,nwpt)
double complex :: vfv(8,nwpt),vfx(8),vsm1(3,3,nwpt),vsm2(3,3,nwpt)
double precision :: vq2(8),vqvtx0(4),vqvtx(4),qkcvt(3,3)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,www,wwwlo,wwwhi,enval(8),dw,eshift,evtx0(4),evtx(4)
double precision :: rrpyr(3,4),kvtx(3,4),xk(3,4),ckvtx(3,4),cxk(3,4),rg(3,3),tvol,xpyr0(4),xpyr(4)
double precision :: xkc(3,4),rgc(3,3),tvolc,kvtxc(3,4)
double complex :: vpyr0(4,nwpt),vpyr(4)
double precision :: de21,de31,de32,de41,de42,de43,thresh
double precision :: fbx(4),fb(4),cmx(3,4),cm(3,4),xkt(3,3)
double complex :: aa0(nwpt),av(3,nwpt)
double precision :: xme,xv(3)
double precision :: avec(3),bgrad(3),xmult
double precision :: abr,rlr,rlr0
integer :: iwh,iwl,ibmin,ibmax,iwhi,iwlo,iehi,ielo
integer :: indxe(4),iwwlim(4)
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),kf,qq(3),qq2,qp(3),qp2,kv(3),kv2,kmq(3),kmq2,qs(3),ek,ekmq
double precision :: fratio
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),temparray(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iipw,jjpw,jjpwt
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwx,npwx)
integer :: pwsymndx(npwc,2*nsym)
integer :: ipw2,jpw1,jpw2,jw,iqcentr,jqcentr
double precision :: eps1,eps2,eps
double precision :: eval,brd,esprd
double precision :: gfo,gamma(3),pln(3),dist(3)
double precision :: sint1,sint1a,sint1b,svec(3),sveca(3),svecb(3),stens(3,3)
double complex :: cint
double complex :: wint,wint0,w2int,w2int0,gterm,vcentr
double precision :: wcentr(-nwpt:nwpt),wgrid(nwpt)
integer :: ntypepaw
double complex :: pwmatel(ntypepaw,nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(ntypepaw,nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                 pwjmatel(3,ntypepaw,nlmn,nlmn),               &
&                tpwjmatel(3,ntypepaw,nlmn,nlmn)
logical :: lqsing,lqcentr
logical :: lx,lc
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
data iter /30/  
logical centerint
integer vcenterint
integer :: i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
common /tsingint_DEBUG/ i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
integer :: idum
double precision :: xdum,xdum2,xdum3,vdum(3),mtxdum(3,3),mtxdum2(3,3)
double complex :: cdum,cdum2,cdum3,cvecdum(3),cvecdum2(3),ctensdum(3,3),ctensdum2(3,3)

!write(6,*) "AAA1"
!xdum = time()

i_prt_DEBUG = 0
j_prt_DEBUG = 0
k_prt_DEBUG = 0
i_ct_DEBUG = 0
j_ct_DEBUG = 0
k_ct_DEBUG = 0
l_ct_DEBUG = 0
i_ref_DEBUG = -1
k_ref_DEBUG = -1
abr=1.d-100
rlr0=3.d-2
rlr=1.d-3
sei=0.d0
ctest=(0.d0,0.d0)
cse=(0.d0,0.d0)
dcse=(0.d0,0.d0)
zz=(0.d0,0.d0)
do ie=-nwpt,nwpt
  cse1(ie)=(0.d0,0.d0)
enddo
xse=0.d0
dw=wmax/dble(nwpt)
!DEBUG
idum = (9.879021/27.2114)/dw+10 
!END DEBUG
do iw=1,nwpt
  wgrid(iw)=iw*dw
enddo
volelmnt=(vol*ngkpt(1)*ngkpt(2)*ngkpt(3))/((2.d0*pi)**3)
volred = 1.d0
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
isym=isymndx(ikk(1),ikk(2),ikk(3))
igsymk=isymg(:,ikpt,isym)
do ii=1,3
do jj=1,3
  qkcvt(ii,jj)=blat(ii,jj)/dble(ngkpt(ii))
enddo
enddo
ictr=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
iqctr=iqndx(ictr(1),ictr(2),ictr(3))
nsing=0
do ii=1,3
  kv(ii) = dble(ikk(ii)-ictr(ii))/dble(ngkpt(ii)) + igg0(ii)
enddo
kv2 = 0.d0
do ii=1,3
do jj=1,3
  kv2=kv2+kv(ii)*bmet(ii,jj)*kv(jj)
enddo
enddo
ek=kv2/(2.d0*pwse_mass)
kf = (3*pi*pi*(nbocc-nbcore)/vol)**(1.d0/3.d0)
ocsign=-1

do ii=1,3
do kk=1,3
  bmet2(ii,kk) = 0.d0                   ! bmet2 = bmet*bmet
  do jj=1,3
    bmet2(ii,kk) = bmet2(ii,kk) + bmet(ii,jj)*bmet(jj,kk)
  enddo
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
!write(6,'(2f10.6,3x,3f6.2)') ek,9.879021/27.2114,kv
  seibc=sei
  cse2=cse1
  xse2=xse
  do ipw1=1,npwc
!write(6,*) ipw1,npwc
!write(33,*)
!write(33,*) ipw1,npwc
    seipw=sei
    do ie=-nwpt,nwpt
      ssc(ie)=(0.d0,0.d0)
    enddo
    igg1=ipwx(:,ipw1)
    lx=.false.
    lc=ipw1.le.npwc  ! do correlation
    if ((.not.lx).and.(.not.lc)) cycle

! Step 1: Compute matrix elements
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=10,10
!    do iz=10,10
      vqmat2(ix,iy,iz)=(0.d0,0.d0)
      ixx=(/ix,iy,iz/)
      iqq=ixx-ictr
      qq = dble(iqq)/dble(ngkpt)
      jka=ikk-iqq
      jkk=mod(jka,ngkpt)
      do ii=1,3
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
      jkb=nint(kpt(:,ikptq)*dble(ngkpt))+ictr
      igg=nint(dble(jkk-jka)/dble(ngkpt))
!write(6,*) '>>>>> iqq = ',iqq
!write(6,*) '>>>>> ikk = ',ikk
!write(6,*) '>>>>> jkk = ',jkk
!write(6,*) '>>>>> jka = ',jka
!write(6,*) '>>>>> igg = ',igg
!write(6,'(a,3f10.5)') 'kpt = ',kpt(:,ikpt)
!write(6,'(a,3f10.5)') 'kptq = ',kpt(:,ikptq)
!read(*,*)
      qq=qq+ipwx(:,ipw1)
      kmq = kv - qq
      qq2=0.d0
      kmq2 = 0.d0
      do ii=1,3
      do jj=1,3
        qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
        kmq2=kmq2+kmq(ii)*bmet(ii,jj)*kmq(jj)
      enddo
      enddo
      vq(ix,iy,iz)=4.d0*pi/qq2
      ekmq = kmq2/(2.d0*pwse_mass)
!write(6,'(a,3f10.5)') 'qq = ',qq
!write(6,'(a,3f10.5)') 'qq2 = ',qq2
!write(6,*) 'vq = ',vq(ix,iy,iz)
!write(6,*) 'jka = ',jka
!write(6,*) 'jkk = ',jkk
!write(6,*) 'jkb = ',jkb
!write(6,*) 'igg = ',igg
!write(6,*) 'ikpt = ',ikpt
!write(6,*) 'ikptq = ',ikptq
      omega(ix,iy,iz)=ekmq-ek+dw/2
!write(33,'(3i3,2x,4f10.6)') ixx,ek,ekmq,ek-ekmq,9.879021/27.2114
      lqsing=.false.
      xmat2(ix,iy,iz)=(1.d0,0.d0)
      if (iqq(1).eq.0.and.iqq(2).eq.0.and.iqq(3).eq.0) then
        jmat2 = bmetinv*4.d0*pi
        if (mode.eq.1) jmat2 = jmat2*(-omega(ix,iy,iz))
        smat2=0.d0
        if (ipw1.eq.1) then
          smat2=4.d0*pi
          lqsing=.true.
          if (mode.eq.1) smat2 = smat2*(-omega(ix,iy,iz))
        endif
      else
        vqmat2(ix,iy,iz)=xmat2(ix,iy,iz)*vq(ix,iy,iz)
        if (mode.eq.1) vqmat2(ix,iy,iz) = vqmat2(ix,iy,iz)*(-omega(ix,iy,iz))
      endif
!if (ipw1.eq.59.and.(ix.eq.7.or.ix.eq.8).and.(iy.eq.9.or.iy.eq.10).and.(iz.eq.7.or.iz.eq.8)) then
!endif
    enddo
    enddo
    enddo
!read(*,*)

! DEBUG
!write(6,*) "Matrix elements computed"

    if (itetrahedron.eq.0) then
      write(6,*) "Sum over points not yet implemented in self energy calculation"
      write(6,*) "set itetrahedron=1 and run again"
    else 
! tetrahedron integration

! Don't waste time integrating over energies with no contribution
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
        ixx=(/ix,iy,iz/)
        call fhilo(omega,ixx,ngkpt,whi,wlo)
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
      if (ocsign.lt.0.d0) then
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
            ixx=(/ix,iy,iz/)
            iqpt=iqndx(ixx(1),ixx(2),ixx(3))
            call locateelement(ixx,ipw1,ipw1,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
            if (iqpt.eq.iqctr) then
              vfactor(ixx(1),ixx(2),ixx(3),iw)=(0.d0,0.d0)
            else if (jjpw.ne.0) then
              vfactor(ixx(1),ixx(2),ixx(3),iw)= &
&                 vqmat2(ixx(1),ixx(2),ixx(3))* &
&                 (lossfn(iw,iqpt,jjpw)-conjg(lossfn(iw,iqpt,jjpw)))*0.5d0
            else
              vfactor(ixx(1),ixx(2),ixx(3),iw)=(0.d0,0.d0)
            endif
!if (ipw1.eq.59.and.iw.eq.247.and.(ix.eq.7.or.ix.eq.8).and.(iy.eq.9.or.iy.eq.10).and.(iz.eq.7.or.iz.eq.8)) then
!write(34,*) ixx,vqmat2(ix,iy,iz),lossfn(iw,iqpt,jjpw)-conjg(lossfn(iw,iqpt,jjpw)),vfactor(ixx(1),ixx(2),ixx(3),iw)
!endif
!if (iw.eq.idum) write(33,'(3i3,2x,2f10.6,2x,f10.6)') ixx,lossfn(iw,iqpt,jjpw),omega(ix,iy,iz)
          enddo
          enddo
          enddo
!write(6,*) "iw = ", iw
          call locateelement(ictr,ipw1,ipw1,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
          do ii=1,3
          do jj=1,3
            if (jjpw.ne.0) then
              vjfactor(ii,jj,iw)=(0.d0,0.d0)
              do kk=1,3
              do ll=1,3
                iiq=ii + 3*(kk-1)
                vjfactor(ii,jj,iw)= &
&                     vjfactor(ii,jj,iw) + &
&                     0.5d0*(lossfn(iw,nqpt+iiq,jjpw)-conjg(lossfn(iw,nqpt+iiq,jjpw))) &
&                     *bmet(kk,ll)*jmat2(ll,jj)
              enddo
              enddo
              iiq=ii + 3*(jj-1)
              vsfactor(ii,jj,iw)=smat2*(lossfn(iw,nqpt+iiq,jjpw)-conjg(lossfn(iw,nqpt+iiq,jjpw)))*0.5d0
            else
              vjfactor(ii,jj,iw)=(0.d0,0.d0)
              vsfactor(ii,jj,iw)=(0.d0,0.d0)
            endif
          enddo
          enddo
!do ii=1,3
!write(6,*) dble(vsfactor(ii,:,iw))
!enddo
!write(6,*)
!do ii=1,3
!write(6,*) dimag(vsfactor(ii,:,iw))
!enddo
!read(*,*)
          do ii=1,3
          do jj=1,3
            vsm1(ii,jj,iw)=(0.d0,0.d0)
            do kk=1,3
              vsm1(ii,jj,iw)=vsm1(ii,jj,iw)+bmet(ii,kk)*vsfactor(kk,jj,iw)
            enddo
          enddo
          enddo
          do ii=1,3
          do jj=1,3
            vsm2(ii,jj,iw)=(0.d0,0.d0)
            do kk=1,3
              vsm2(ii,jj,iw)=vsm2(ii,jj,iw)+vsm1(ii,kk,iw)*bmet(kk,jj)
            enddo
          enddo
          enddo
        enddo
      endif

! DEBUG
!write(6,*) "Beginning integration"

! Do the integration
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!xdum = 0.d0;
!xdum2 = 0.d0;
!xdum3 = 0.d0;
        ixx=(/ix,iy,iz/)
        iqq=ixx-ictr
        ssx1=ssx
        call fval(omega,ixx,ixxp,ngkpt,enval)
        if (lqsing) call fval(vq,ixx,ixxp,ngkpt,vq2)
        if (lc) then
          do iw=1,nwpt
            call fpol(vfactor(:,:,:,iw),ngkpt,ixx,ixxp,vfv(:,iw))
          enddo
        endif
        if (lx) call fpol(vqmat2(:,:,:),ngkpt,ixx,ixxp,vfx)
        do itet=1,6
! DEBUG
!cdum = 0.d0
          centerint = .false.
          do iv=1,4
            evtx0(iv)=enval(ivndx(iv,itet))
            rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
            if (lc) then
              do iw=1,nwpt
                vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
              enddo
            endif
            if (lx) xpyr0(iv)=vfx(ivndx(iv,itet))
            if ((ixx(1)+rrpyr(1,iv).eq.ictr(1)).and.(ixx(2)+rrpyr(2,iv).eq.ictr(2)).and.(ixx(3)+rrpyr(3,iv).eq.ictr(3))) then
              centerint = .true.  ! does tetrahedron contain center vertex?
              vcenterint = iv     ! index of center vertex
            endif
          enddo
!write(6,*) ixx
!write(6,*) itet,centerint
!if (centerint) then
!write(6,*) itet
!do iv=1,4
!write(6,*) ixx+rrpyr(:,iv)
!enddo
!!read(*,*)
!endif
          call indxhpsort(4,4,evtx0,indxe)
          evtx=evtx0(indxe)
          do ii=1,4
            jj=indxe(ii)
            kvtx(:,ii) = (iqq(:) + rrpyr(:,jj))
            do kk=1,3
              kvtxc(kk,ii)=dot_product(iqq+rrpyr(:,jj),qkcvt(:,kk))
            enddo
          enddo
          xk(1:3,1)=(/0.d0,0.d0,0.d0/)
          xkc(1:3,1)=(/0.d0,0.d0,0.d0/)
          do ii=1,3
            xk(:,ii+1)=rrpyr(:,indxe(ii+1))-rrpyr(:,indxe(1))
          enddo
          do ii=2,4
            xkc(1:3,ii)=kvtxc(1:3,ii)-kvtxc(1:3,1)
          enddo
          do jj=1,3  ! make contragredient, reduced coordinates
            call cross(xk(:,mod(jj,3)+2),xk(:,mod(jj+1,3)+2),rg(:,jj))
          enddo
          do jj=1,3  ! make contragredient, Cartesian coordinates
            call cross(xkc(:,mod(jj,3)+2),xkc(:,mod(jj+1,3)+2),rgc(:,jj))
          enddo
          tvol=dot_product(xk(1:3,2),rg(1:3,1))
          rg=rg/tvol
          tvol=abs(tvol)
          tvolc=dot_product(xkc(1:3,2),rgc(1:3,1))
          rgc=rgc/tvolc
          tvolc=abs(tvolc)
          iwwlim=nint(evtx/dw)
!if (centerint) write(6,'(i3,3x,3i3)') itet,ixx
!if (centerint) write(6,*) centerint,lqsing,lx
!if (ek.gt.0.2.and.ipw1.eq.1) write(34,*) "A    ",itet,ixx,ssi(1)
          if (lqsing .and. centerint) then
! singular part
! DEBUG
!xdum = 0.d0
!write(6,*)
!write(6,*) "singular integration"
!write(6,*) 'ixx: ',ixx
!write(6,*) 'itet: ',itet,'  centerint: ',centerint
!write(6,*) 'iwwlim: ',iwwlim
            do iv=1,4
! subtraction of integrated singular function from values at vertices
              ixr=ixx+rrpyr(:,iv)
              iqr = ixr - ictr
              qr(:,iv) = dble(iqr)
              if (ixr(1).eq.ictr(1).and.ixr(2).eq.ictr(2).and.ixr(3).eq.ictr(3)) then
                if (lc) then
                  do iw=1,nwpt
                    vpyr0(iv,iw) = 0.d0
                  enddo
                endif
                if (lx) xpyr0(iv) = 0.d0
              else
                qq = qr(:,iv)/dble(ngkpt)
                qqr(:,iv)=qq
                qq2 = 0.d0
                do ii=1,3
                do jj=1,3
                  qq2 = qq2 + qq(ii)*bmet(ii,jj)*qq(jj)
                enddo
                enddo
                if (lc) then
                  do iw=1,nwpt
                    do ii=1,3
                    do jj=1,3
                      vpyr0(iv,iw)=vpyr0(iv,iw)-qq(ii)*vsm2(ii,jj,iw)*qq(jj)/(qq2**2)
                    enddo
                    enddo
                  enddo
                endif
                if (lx) then 
                  xpyr0(iv)=xpyr0(iv)-smat2/qq2
                endif
              endif
            enddo
            bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
            do jj=1,3
              bgrad=bgrad+(evtx(jj+1)-evtx(1))*rgc(:,jj)
            enddo
            xmult=1.d0/sqrt(dot_product(bgrad,bgrad))
            if (lc) then
              do iww=iwwlim(1)-1,iwwlim(4)+1
                wwwlo=dw*(iww-0.5d0)
                call singular_adaptive_tetint(qqr,evtx,wwwlo,dw,vcenterint,bmet,abr,rlr,rlr0,stens)
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  if (ie.ge.ielo.and.ie.le.iehi) then
                    do ii=1,3
                    do jj=1,3
                      ssi(ie)=ssi(ie)+vsm2(ii,jj,iw)*stens(jj,ii)
! DEBUG
!if (ie.eq.0.or.ie.eq.1) xdum3=xdum3+dimag(vsm2(ii,jj,iw)*stens(jj,ii))
                    enddo
                    enddo
                  endif
                enddo
              enddo
            endif
            if (lx) then
! DEBUG
!cdum = 0.d0
              wwwlo = min(evtx(1),evtx(2),evtx(3),evtx(4))
              dw = max(evtx(1),evtx(2),evtx(3),evtx(4))-wwwlo
              call singular_adaptive_tetint(qqr,evtx,wwwlo,dw,vcenterint,bmet,abr,rlr,rlr0,stens)
              do ii=1,3
              do jj=1,3
                ssx=ssx+smat2*bmet(ii,jj)*stens(jj,ii)
!! DEBUG
!cdum = cdum+smat2*bmet(ii,jj)*stens(jj,ii)
              enddo
              enddo
!write(6,*) ssx
!write(6,'(i3,2x,3i3,2x,2f12.8)') itet,ixx,cdum
            endif
          elseif (centerint) then
            if (lc) then
              do iww=iwwlim(1)-1,iwwlim(4)+1
                wwwlo=dw*(iww-0.5d0)
                do iw=1,nwpt
                  call vnitetint(rrpyr,evtx0,wwwlo,dw,vjfactor(:,:,iw),bmet,vcenterint,cint)
                  ie=iww-ocsign*iw
                  if (ie.ge.ielo.and.ie.le.iehi) ssi(ie)=ssi(ie)+cint/volelmnt
! DEBUG
!if (ie.eq.0.or.ie.eq.1) xdum2=xdum2+dimag(cint/volelmnt)
                enddo
              enddo
            endif
            if (lx) then
              call vnitetint(rrpyr,evtx0,evtx0(indxe(1)),evtx0(indxe(4))-evtx0(indxe(1)),jmat2,bmet,vcenterint,cint)
              ssx=ssx+cint/volelmnt
            endif
          endif
!if (ek.gt.0.2.and.ipw1.eq.1) write(34,*) "     ",itet,ixx,ssi(1)
! can use basic tetrahedron integration for everything else
! DEBUG
!cdum = 0.d0
!if ((ix.eq.1.or.ix.eq.2).and.(iy.eq.9.or.iy.eq.10).and.(iz.eq.8.or.iz.eq.9)) then
!write(6,*) "vpyr0", vpyr0(:,idum)
!write(6,*) "evtx",evtx
!write(6,*) "www",idum*dw
!write(6,*) "iw",idum
!write(6,*) "iwwlim",iwwlim
!read(*,*)
!endif
          if (lc) then
            do iw=1,nwpt
              aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
              av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
              do jj=1,3
                av(1:3,iw)=av(1:3,iw) &
&                      +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
              enddo
            enddo
            do iww=iwwlim(1)-1,iwwlim(4)+1
              wwwlo=dw*(iww-0.5d0)
!              wwwhi=dw*(iww+0.5d0)
              call mkvsint(rrpyr(1:3,indxe),xk,volred,evtx,wwwlo,dw,sint1,svec)
              do iw=1,nwpt
                ie=iww-ocsign*iw
                if (ie.ge.ielo.and.ie.le.iehi) ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
! DEBUG
!if (ipw1.eq.59.and.(ix.eq.7).and.(iy.eq.9).and.(iz.eq.7)) then
!if ((ie.eq.0.or.ie.eq.1).and.iw.eq.247) then
!write(34,*) ixx,iw,wwwlo+dw/2.d0,(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
!write(34,*) "                  ",vpyr0(indxe(1),iw)
!write(34,*) "                  ",vpyr0(indxe(2),iw)
!write(34,*) "                  ",vpyr0(indxe(3),iw)
!write(34,*) "                  ",vpyr0(indxe(4),iw)
!endif
!endif
!if (ie.eq.0.or.ie.eq.1) xdum = xdum+dimag((aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt)
!if (ie.eq.0.or.ie.eq.1) then
!if (ix.eq.1.and.iy.eq.1.and.iz.eq.4.and.itet.eq.1) then
!if ((ix.eq.1.or.ix.eq.2).and.(iy.eq.9.or.iy.eq.10).and.(iz.eq.8.or.iz.eq.9)) then
!if (iw.eq.idum) then
!write(6,*) "aa0",aa0(iw)
!write(6,*) "av",av(:,iw)
!write(6,*) "sint1",sint1
!write(6,*) "svec",svec
!write(6,*) "evtx",evtx
!write(6,*) "wwwlo",wwwlo
!write(6,*) "wwwhi",wwwlo+dw
!write(6,*) "vpyr0", vpyr0(:,iw)
!write(6,*) "vfv", vfv(:,iw)
!read(*,*)
!endif
!endif
!endif
              enddo
            enddo
          endif
! DEBUG
!write(32,'(4i3,2x,f14.10,2x,f14.10)') ixx,itet,(cdum/2.d0)
          if (lx) then
            xme=xpyr0(indxe(1))
            xv=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              xv=xv+(xpyr0(indxe(jj+1))-xpyr0(indxe(1)))*rg(:,jj)
            enddo
            call mkvsint(rrpyr(1:3,indxe),xk,volred,evtx,evtx(1),evtx(4)-evtx(1),sint1,svec)
            ssx=ssx+(xme*sint1+dot_product(svec,xv))/volelmnt
! DEBUG
!write(32,'(4i3,2x,f16.10)') ixx,itet,(xme*sint1+dot_product(svec,xv))/volelmnt
!write(32,'(4i3,2x,f14.10,2x,3f14.10)') ixx,itet, &
!& (xme*sint1+dot_product(svec,xv))/volelmnt, &
!& (xpyr0(indxe(1))+xpyr0(indxe(2))+xpyr0(indxe(3))+xpyr0(indxe(4))) &
!& /(volelmnt*6.d0*4.d0)
          endif
!write(33,'(3i3,2x,i3,2x,4f10.6,2x,2f10.6)') ixx,itet,evtx,cdum
!if (ek.gt.0.2.and.ipw1.eq.1) write(34,*) "     ",itet,ixx,ssi(1)
        enddo
! DEBUG
!write(32,'(3i3,2x,4f16.10)') ixx,ssi(0),ssi(1)
!write(33,*) ixx,xdum,xdum2,xdum3
      enddo
      enddo
      enddo
!    if (lc) ssi=-ssi/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
!    if (lx) ssx=-ssx/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
      if (lc) ssi=-ocsign*ssi/((2*pi)**3*pi)
      if (lx) ssx=-ssx/((2*pi)**3)
  
! DEBUG
!write(6,*) "Integration complete"

      if (lc) then
        do ie=-nwpt,nwpt
          eps1=dble(ie)*dw-dw/2
!        ssc(ie)=(0.d0,0.d0)
          do je=ielo,iehi
            if (ie.eq.je) then
              ssc(ie)=ssc(ie)+pi*ssi(je)
            else
              eps2=dble(je)*dw-dw/2
              ssc(ie)=ssc(ie)-(0.d0,1.d0)*ocsign*ssi(je)*dw/(eps2-eps1)
            endif
!if (ie.eq.0) write(16,'(2i4,2f10.3,es12.3,2(2x,2es11.3))') ie,je,eps1*27.2114,eps2*27.2114,1/(eps2-eps1),ssc(ie)*27.2114,dble(ssi(je)*27.2114)
!if (ie.eq.0) write(16,'(2x,i4,2x,f10.3,4x,f10.6,2x,f10.6)') je,eps2*27.2114,ssi(je)*27.2114
          enddo
!if (ipw1.eq.1) write(34,*) ie, ssi(ie)
        enddo
      endif
!do je=ielo,iehi
!write(6,*) je,ssi(je)
!enddo
      deallocate (ssi)
!write(6,'("lstv",2x,f14.10,4x,"(",f14.10,",",f14.10,")")') ssx,(ssc(0)+ssc(1))/2.

      if (lc) cse1=cse1+ssc
      if (lx) xse=xse+ssx
    endif

! DEBUG
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

!write(6,*) "AAA2", time()-xdum

return
end subroutine mkpwse

!*************************************************************************

subroutine mkxtpwse(igg0,ikk,mode,xtlossfn, &
& pwse_mass, &
& vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmn, &
& test_bands_se, &
& ntypepaw, &
& pwmatel,tpwmatel, pwjmatel, tpwjmatel, &
& kg,enrgy,cg,npwt,bantot,ncg, &
& indxkpw,indxkbnd,indxkcg,npwarr, &
& kpt,nkpt,nqpt,nxtqpt,nsymk,symk,nsym,symrel,syminv, &
& ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
& ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
& pwlo,pwhi,npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
& igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
& nxtgkpt,xtqpt,xtqpt2qpt,qpt2xtqpt,qpt2xtg,ixtqndx, &
& nband,nsppol,shiftk,zz,dcse,cse,xse)
implicit none
integer :: mode
! mode = 0 for loss rate
! mode = 1 for power loss
integer :: ipwv,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: nxtgkpt(3)
integer :: igmn(3),igmx(3)
integer :: npwt,bantot,ncg,nkpt,nqpt,nxtqpt,nsym,nsppol
integer :: npw,npwc,npwx,pwlo,pwhi,ipw1,npwndx,ntpwndx,napwndx,ipaw,itetrahedron
double precision :: vol,pi,wmax,xred(3,natom)
double precision :: pwse_mass
double complex :: xtlossfn(nwpt,nxtqpt,npwx)
integer :: test_bands_se(2)
double complex :: projwf(natom,nlmn,nkpt,ncband)
integer :: kg(3,npwt)
double precision :: enrgy(bantot)
double complex :: cg(ncg)
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt),shiftk(3)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
double precision :: bmet(3,3),blat(3,3),bmetinv(3,3),bmet_t(3,3),bmet2(3,3)
double precision :: volelmnt, volred
integer :: ipwx(3,npwx),ipwndx(2,napwndx)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ixtqndx(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
double precision :: xtqpt(3,nxtqpt)
integer :: xtqpt2qpt(nxtqpt),qpt2xtqpt(nqpt),qpt2xtg(nqpt,3)
integer :: nband(nkpt*nsppol)
double complex :: sei,seipw,seiold,seibc,csepw,cseold,csebc
double complex, allocatable :: ssi(:)
double complex :: ssc(-nwpt:nwpt),cse,dcse,cse0,dcse0,cse1(-nwpt:nwpt),zz,ctest,cse2(-nwpt:nwpt)
double precision :: xse,xse0,ssx,ssx1,ssx2,xse2,ssx3
integer :: ii,jj,kk,ll,ix,iy,iz,iskip,ocsign,iiq
integer :: iqq(3),jka(3),jkb(3),jkk(3),iqv(3),ictr(3),ixtctr(3),iqr(3),ixr(3)
double precision :: qr(3,4),qqr(3,4)
integer :: ixx(3),ixt(3),ixxp(3),ixtp(3),ixp2(3),ixp1(3),ixm1(3),ixm2(3)
integer :: ikpt,ikptq
integer :: igg(3),igh(3),igg0(3),igg1(3),igg2(3),idgg1(3),idgg2(3),isgg1(3),isgg(2),igq
integer :: isym,isymq,iqpt,ixtqpt,iqsym,iw,iww,ie,je,ies,iqctr,ibp
double precision :: xck(3),xckq(3),qadj(6)
double complex :: cmatel,cmatel2,amatel,amatel2
double complex :: jmatel(3),jmatel2(3),vmatel(3),vmatel2(3)
double complex :: vqmat2(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3))
double complex :: xmat2(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3))
double complex :: jmat2(3,3), rjmat2(3,3), rjmatdum
double precision :: smat2
double complex :: vfactor(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3),nwpt)
double complex :: vfv(8,nwpt),vfx(8)
double precision :: vq2(8),vqvtx0(4),vqvtx(4),qkcvt(3,3)
double precision :: omega(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3))
double precision :: whi,wlo,wwhi,wwlo,ww,www,wwwlo,wwwhi,enval(8),dw,eshift,evtx0(4),evtx(4)
double precision :: rrpyr(3,4),kvtx(3,4),xk(3,4),ckvtx(3,4),cxk(3,4),rg(3,3),tvol,xpyr0(4),xpyr(4)
double precision :: xkc(3,4),rgc(3,3),tvolc,kvtxc(3,4)
double complex :: vpyr0(4,nwpt),vpyr(4)
double precision :: de21,de31,de32,de41,de42,de43,thresh
double precision :: fbx(4),fb(4),cmx(3,4),cm(3,4),xkt(3,3)
double complex :: aa0(nwpt),av(3,nwpt)
double precision :: xme,xv(3)
double precision :: avec(3),bgrad(3),xmult
double precision :: abr,rlr,rlr0
integer :: iwh,iwl,ibmin,ibmax,iwhi,iwlo,iehi,ielo
integer :: indxe(4),iwwlim(4)
double precision :: vq(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3))
double precision :: kf,qq(3),qq2,qp(3),qp2,kv(3),kv2,kmq(3),kmq2,qs(3),ek,ekmq
double precision :: fratio
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),temparray(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: iipw,jjpw,jjpwt
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwx,npwx)
integer :: pwsymndx(npwc,2*nsym)
integer :: ipw2,jpw1,jpw2,jw,iqcentr,jqcentr
double precision :: eps1,eps2,eps
double precision :: eval,brd,esprd
double precision :: gfo,gamma(3),pln(3),dist(3)
double precision :: sint1,sint1a,sint1b,svec(3),sveca(3),svecb(3),stens(3,3)
double complex :: cint
double complex :: wint,wint0,w2int,w2int0,gterm,vcentr
double precision :: wcentr(-nwpt:nwpt),wgrid(nwpt)
integer :: ntypepaw
double complex :: pwmatel(ntypepaw,nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(ntypepaw,nlmn,nlmn,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                 pwjmatel(3,ntypepaw,nlmn,nlmn),               &
&                tpwjmatel(3,ntypepaw,nlmn,nlmn)
logical :: lqsing,lqcentr
logical :: lx,lc
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
data iter /30/  
logical centerint
integer vcenterint
integer :: i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
common /tsingint_DEBUG/ i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
integer :: idum
double precision :: xdum,xdum2,xdum3,vdum(3),mtxdum(3,3),mtxdum2(3,3)
double complex :: cdum,cdum2,cdum3,cvecdum(3),cvecdum2(3),ctensdum(3,3),ctensdum2(3,3)

!write(6,*) "AAA1"
!xdum = time()
!write(6,*) xse

i_prt_DEBUG = 0
j_prt_DEBUG = 0
k_prt_DEBUG = 0
i_ct_DEBUG = 0
j_ct_DEBUG = 0
k_ct_DEBUG = 0
l_ct_DEBUG = 0
i_ref_DEBUG = -1
k_ref_DEBUG = -1
abr=1.d-100
rlr0=3.d-2
rlr=1.d-3
sei=0.d0
ctest=(0.d0,0.d0)
cse=(0.d0,0.d0)
dcse=(0.d0,0.d0)
zz=(0.d0,0.d0)
do ie=-nwpt,nwpt
  cse1(ie)=(0.d0,0.d0)
enddo
xse=0.d0
dw=wmax/dble(nwpt)
!DEBUG
idum = (9.879021/27.2114)/dw+10 
!END DEBUG
do iw=1,nwpt
  wgrid(iw)=iw*dw
enddo
volelmnt=(vol*nxtgkpt(1)*nxtgkpt(2)*nxtgkpt(3))/((2.d0*pi)**3)
volred = 1.d0
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
isym=isymndx(ikk(1),ikk(2),ikk(3))
igsymk=isymg(:,ikpt,isym)
do ii=1,3
do jj=1,3
  qkcvt(ii,jj)=blat(ii,jj)/dble(nxtgkpt(ii))
enddo
enddo
ictr=(/ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2/)
ixtctr=(/nxtgkpt(1)/2,nxtgkpt(2)/2,nxtgkpt(3)/2/)
iqctr=iqndx(ictr(1),ictr(2),ictr(3))
nsing=0
do ii=1,3
  kv(ii) = dble(ikk(ii)-ictr(ii))/dble(ngkpt(ii)) + igg0(ii)
enddo
kv2 = 0.d0
do ii=1,3
do jj=1,3
  kv2=kv2+kv(ii)*bmet(ii,jj)*kv(jj)
enddo
enddo
ek=kv2/(2.d0*pwse_mass)
kf = (3*pi*pi*(nbocc-nbcore)/vol)**(1.d0/3.d0)
ocsign=-1

do ii=1,3
do kk=1,3
  bmet2(ii,kk) = 0.d0                   ! bmet2 = bmet*bmet
  do jj=1,3
    bmet2(ii,kk) = bmet2(ii,kk) + bmet(ii,jj)*bmet(jj,kk)
  enddo
enddo
enddo

!do ix=1,ngkpt(1)
!do iy=1,ngkpt(2)
!do iz=1,ngkpt(3)
!  temparray(ix,iy,iz)=(0.d0,0.d0)
!enddo
!enddo
!enddo

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
!write(6,'(2f10.6,3x,3f6.2)') ek,9.879021/27.2114,kv
  seibc=sei
  cse2=cse1
  xse2=xse
  do ipw1=pwlo,pwhi
    seipw=sei
    do ie=-nwpt,nwpt
      ssc(ie)=(0.d0,0.d0)
    enddo
    igg1=ipwx(:,ipw1)
    lx=.false.
    lc=ipw1.le.npwx  ! do correlation
    if ((.not.lx).and.(.not.lc)) cycle

! Step 1: Compute matrix elements
    do ix=1,nxtgkpt(1)
    do iy=1,nxtgkpt(2)
    do iz=1,nxtgkpt(3)
!    do ix=1,1
!    do iy=10,10
!    do iz=10,10
      vqmat2(ix,iy,iz)=(0.d0,0.d0)
      ixt=(/ix,iy,iz/)
      ixx=ixt*ngkpt/nxtgkpt
      iqq=ixx-ictr
      qq = dble(iqq)/dble(ngkpt)
      jka=ikk-iqq
      jkk=mod(jka,ngkpt)
      do ii=1,3
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
      jkb=nint(kpt(:,ikptq)*dble(ngkpt))+ictr
      igg=nint(dble(jkk-jka)/dble(ngkpt))
!write(6,*) '>>>>> ixt = ',ixt
!write(6,*) '>>>>> ixx = ',ixx
!write(6,*) '>>>>> iqq = ',iqq
!write(6,*) '>>>>> ikk = ',ikk
!write(6,*) '>>>>> jkk = ',jkk
!write(6,*) '>>>>> jka = ',jka
!write(6,*) '>>>>> igg = ',igg
!write(6,*) '>>>>> ikptq = ',ikptq,nkpt
!write(6,'(a,3f10.5)') 'kpt = ',kpt(:,ikpt)
!write(6,'(a,3f10.5)') 'kptq = ',kpt(:,ikptq)
!read(*,*)
      qq=qq+ipwx(:,ipw1)
      kmq = kv - qq
      qq2=0.d0
      kmq2 = 0.d0
      do ii=1,3
      do jj=1,3
        qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
        kmq2=kmq2+kmq(ii)*bmet(ii,jj)*kmq(jj)
      enddo
      enddo
      vq(ix,iy,iz)=4.d0*pi/qq2
      ekmq = kmq2/(2.d0*pwse_mass)
!write(6,'(a,3f10.5)') 'qq = ',qq
!write(6,'(a,3f10.5)') 'qq2 = ',qq2
!write(6,*) 'vq = ',vq(ix,iy,iz)
!write(6,*) 'jka = ',jka
!write(6,*) 'jkk = ',jkk
!write(6,*) 'jkb = ',jkb
!write(6,*) 'igg = ',igg
!write(6,*) 'ikpt = ',ikpt
!write(6,*) 'ikptq = ',ikptq
      omega(ix,iy,iz)=ekmq-ek+dw/2
!write(33,'(3i3,2x,4f10.6)') ixt,ek,ekmq,ek-ekmq,9.879021/27.2114
      lqsing=.false.
      xmat2(ix,iy,iz)=(1.d0,0.d0)
      vqmat2(ix,iy,iz)=xmat2(ix,iy,iz)*vq(ix,iy,iz)
      if (mode.eq.1) vqmat2(ix,iy,iz) = vqmat2(ix,iy,iz)*(-omega(ix,iy,iz))
    enddo
    enddo
    enddo
!read(*,*)

! DEBUG
!write(6,*) "Matrix elements computed"

    if (itetrahedron.eq.0) then
      write(6,*) "Sum over points not yet implemented in self energy calculation"
      write(6,*) "set itetrahedron=1 and run again"
    else 
! tetrahedron integration

! Don't waste time integrating over energies with no contribution
      do ix=1,nxtgkpt(1)
      do iy=1,nxtgkpt(2)
      do iz=1,nxtgkpt(3)
        ixt=(/ix,iy,iz/)
        ixx=ixt*ngkpt/nxtgkpt
        call fhilo(omega,ixt,nxtgkpt,whi,wlo)
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
      if (ocsign.lt.0.d0) then
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
          do ix=1,nxtgkpt(1)
          do iy=1,nxtgkpt(2)
          do iz=1,nxtgkpt(3)
            ixt=(/ix,iy,iz/)
            ixx=ixt*ngkpt/nxtgkpt
            ixtqpt = ixtqndx(ixt(1),ixt(2),ixt(3))
            if (ipw1.le.npwx) then
              vfactor(ixt(1),ixt(2),ixt(3),iw)= &
&                 vqmat2(ixt(1),ixt(2),ixt(3))* &
&                 (xtlossfn(iw,ixtqpt,ipw1)-conjg(xtlossfn(iw,ixtqpt,ipw1)))*0.5d0
            else
              vfactor(ixt(1),ixt(2),ixt(3),iw)=(0.d0,0.d0)
            endif
!if (ipw1.eq.59.and.iw.eq.247.and.(ix.eq.7.or.ix.eq.8).and.(iy.eq.9.or.iy.eq.10).and.(iz.eq.7.or.iz.eq.8)) then
!write(34,*) ixt,vqmat2(ix,iy,iz),xtlossfn(iw,ixtqpt,jjpw)-conjg(xtlossfn(iw,ixtqpt,jjpw)),vfactor(ixt(1),ixt(2),ixt(3),iw)
!endif
!if (iw.eq.idum) write(33,'(3i3,2x,2f10.6,2x,f10.6)') ixt,xtlossfn(iw,ixtqpt,jjpw),omega(ix,iy,iz)
          enddo
          enddo
          enddo
!write(6,*) "iw = ", iw
        enddo
      endif

! DEBUG
!write(6,*) "Beginning integration"

! Do the integration
      do ix=1,nxtgkpt(1)
      do iy=1,nxtgkpt(2)
      do iz=1,nxtgkpt(3)
!xdum = 0.d0;
!xdum2 = 0.d0;
!xdum3 = 0.d0;
        ixt=(/ix,iy,iz/)
        ixx=ixt*ngkpt/nxtgkpt
        iqq=ixt-ixtctr
        ssx1=ssx
        call fval(omega,ixt,ixtp,nxtgkpt,enval)
        if (lqsing) call fval(vq,ixt,ixtp,nxtgkpt,vq2)
        if (lc) then
          do iw=1,nwpt
            call fpol(vfactor(:,:,:,iw),nxtgkpt,ixt,ixtp,vfv(:,iw))
          enddo
        endif
        if (lx) call fpol(vqmat2(:,:,:),nxtgkpt,ixt,ixtp,vfx)
        do itet=1,6
! DEBUG
!cdum = 0.d0
          do iv=1,4
            evtx0(iv)=enval(ivndx(iv,itet))
            rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
            if (lc) then
              do iw=1,nwpt
                vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
              enddo
            endif
            if (lx) xpyr0(iv)=vfx(ivndx(iv,itet))
          enddo
          call indxhpsort(4,4,evtx0,indxe)
          evtx=evtx0(indxe)
          do ii=1,4
            jj=indxe(ii)
            kvtx(:,ii) = (iqq(:) + rrpyr(:,jj))
            do kk=1,3
              kvtxc(kk,ii)=dot_product(iqq+rrpyr(:,jj),qkcvt(:,kk))
            enddo
          enddo
          xk(1:3,1)=(/0.d0,0.d0,0.d0/)
          xkc(1:3,1)=(/0.d0,0.d0,0.d0/)
          do ii=1,3
            xk(:,ii+1)=rrpyr(:,indxe(ii+1))-rrpyr(:,indxe(1))
          enddo
          do ii=2,4
            xkc(1:3,ii)=kvtxc(1:3,ii)-kvtxc(1:3,1)
          enddo
          do jj=1,3  ! make contragredient, reduced coordinates
            call cross(xk(:,mod(jj,3)+2),xk(:,mod(jj+1,3)+2),rg(:,jj))
          enddo
          do jj=1,3  ! make contragredient, Cartesian coordinates
            call cross(xkc(:,mod(jj,3)+2),xkc(:,mod(jj+1,3)+2),rgc(:,jj))
          enddo
          tvol=dot_product(xk(1:3,2),rg(1:3,1))
          rg=rg/tvol
          tvol=abs(tvol)
          tvolc=dot_product(xkc(1:3,2),rgc(1:3,1))
          rgc=rgc/tvolc
          tvolc=abs(tvolc)
          iwwlim=nint(evtx/dw)
          if (lc) then
            do iw=1,nwpt
              aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
              av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
              do jj=1,3
                av(1:3,iw)=av(1:3,iw) &
&                      +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
              enddo
            enddo
            do iww=iwwlim(1)-1,iwwlim(4)+1
              wwwlo=dw*(iww-0.5d0)
!              wwwhi=dw*(iww+0.5d0)
              call mkvsint(rrpyr(1:3,indxe),xk,volred,evtx,wwwlo,dw,sint1,svec)
              do iw=1,nwpt
                ie=iww-ocsign*iw
                if (ie.ge.ielo.and.ie.le.iehi) ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
              enddo
            enddo
          endif
! DEBUG
!write(32,'(4i3,2x,f14.10,2x,f14.10)') ixt,itet,(cdum/2.d0)
          if (lx) then
            xme=xpyr0(indxe(1))
            xv=(/0.d0,0.d0,0.d0/)
            do jj=1,3
              xv=xv+(xpyr0(indxe(jj+1))-xpyr0(indxe(1)))*rg(:,jj)
            enddo
            call mkvsint(rrpyr(1:3,indxe),xk,volred,evtx,evtx(1),evtx(4)-evtx(1),sint1,svec)
            ssx=ssx+(xme*sint1+dot_product(svec,xv))/volelmnt
          endif
!write(33,'(3i3,2x,i3,2x,4f10.6,2x,2f10.6)') ixt,itet,evtx,cdum
        enddo
      enddo
      enddo
      enddo
!    if (lc) ssi=-ssi/(dble(nxtgkpt(1)*nxtgkpt(2)*nxtgkpt(3))*vol*pi)
!    if (lx) ssx=-ssx/(dble(nxtgkpt(1)*nxtgkpt(2)*nxtgkpt(3))*vol)
      if (lc) ssi=-ocsign*ssi/((2*pi)**3*pi)
      if (lx) ssx=-ssx/((2*pi)**3)
  
! DEBUG
!write(6,*) "Integration complete"

      if (lc) then
        do ie=-nwpt,nwpt
          eps1=dble(ie)*dw-dw/2
!        ssc(ie)=(0.d0,0.d0)
          do je=ielo,iehi
            if (ie.eq.je) then
              ssc(ie)=ssc(ie)+pi*ssi(je)
            else
              eps2=dble(je)*dw-dw/2
              ssc(ie)=ssc(ie)-(0.d0,1.d0)*ocsign*ssi(je)*dw/(eps2-eps1)
            endif
!if (ie.eq.0) write(16,'(2i4,2f10.3,es12.3,2(2x,2es11.3))') ie,je,eps1*27.2114,eps2*27.2114,1/(eps2-eps1),ssc(ie)*27.2114,dble(ssi(je)*27.2114)
!if (ie.eq.0) write(16,'(2x,i4,2x,f10.3,4x,f10.6,2x,f10.6)') je,eps2*27.2114,ssi(je)*27.2114
          enddo
        enddo
      endif
!do je=ielo,iehi
!write(6,*) je,ssi(je)
!enddo
      deallocate (ssi)
!write(6,'("lstv",2x,f14.10,4x,"(",f14.10,",",f14.10,")")') ssx,(ssc(0)+ssc(1))/2.

      if (lc) cse1=cse1+ssc
      if (lx) xse=xse+ssx
    endif

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

!write(6,*) "AAA2", time()-xdum

return
end subroutine mkxtpwse

!*************************************************************************

double precision function re_dielf_feg(omega)
implicit none
double precision :: omega
double precision :: qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi
double precision :: NN
double precision :: EE,GG,hh,pf,AA,BB,Elo,Ehi
double precision :: lnterm1,lnterm2,lnterm3,lnterm4,lnterm5
double precision :: xterm1,xterm2,xterm3,xterm4,xterm5,rterm
common /feg/ qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

  EE = omega/Ef
  GG = Eg/Ef
  hh = qq/xkf
  pf = (3.d0/4.d0)*(omegap2/(Ef*Ef))
  AA = GG + hh*hh
  BB = 0.25/(hh*hh)
  Ehi = GG + hh*abs(hh+2)
  Ehi = GG + hh*abs(hh-2)
  if (hh.le.2.d0) then
    NN = (3*GG+4*hh-0.25*GG*hh*hh)*hh*hh
    rterm = 2*(Elo-GG) - 4*hh*hh*(GG+2*hh-2*AA)*BB
  else 
    NN = 4*hh*(GG+hh*hh)
    rterm = 8*hh*AA*BB
  endif
  xterm1 = 1.d0 - (AA+EE)**2/BB
  xterm2 = 1.d0 - (AA-EE)**2/BB
  lnterm1 = log(abs((Ehi+EE)/(Elo+EE)))
  lnterm2 = log(abs((Ehi-EE)/(Elo-EE)))
  re_dielf_feg = 1.d0 + (pf/NN)*(rterm + xterm1*lnterm1 + xterm2*lnterm2)
  if (hh.le.2.d0) then
    xterm3 = -(GG+EE)
    xterm4 = -(GG-EE)
    lnterm3 = log(abs((Elo+EE)/(GG+EE)))
    lnterm4 = log(abs((Elo-EE)/(GG-EE)))
    re_dielf_feg = re_dielf_feg + (pf/NN)*(xterm3*lnterm3 + xterm4*lnterm4)
  endif

end function re_dielf_feg

!*************************************************************************

double precision function im_dielf_feg(omega)
implicit none
double precision :: omega
double precision :: qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi
double precision :: term,NN
double precision :: EE,GG,hh,pf
common /feg/ qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

!write(*,*) edens,Eg,omegap2
!read(*,*)
  EE = omega/Ef
  GG = Eg/Ef
  hh = qq/xkf
  pf = (3.d0*pi/4.d0)*(omegap2/(Ef*Ef))
!write(*,*) pf,omegap2,Ef
!read(*,*)
  if (hh.le.2.d0) then
    NN = (3*GG+4*hh-0.25*GG*hh*hh)*hh*hh
    if (EE.le.GG) then
      im_dielf_feg = 0.d0
    else if (EE.le.(GG+hh*(2.d0-hh))) then
      im_dielf_feg = pf * (EE-GG)/NN
    else if (EE.le.(GG+hh*(2.d0+hh))) then
      term = (EE-GG-hh*hh)/(2*hh)
      im_dielf_feg = pf * (1.d0 - term**2)/NN
    else
      im_dielf_feg = 0.d0
    endif
  else
    NN = 4*hh*(GG+hh*hh)
    if (EE.le.(GG+hh*(hh-2.d0))) then
      im_dielf_feg = 0.d0
    else if (EE.le.(GG+hh*(2.d0+hh))) then
      term = (EE-GG-hh*hh)/(2*hh)
      im_dielf_feg = pf * (1.d0 - term**2)/NN
    else
      im_dielf_feg = 0.d0
    endif
  endif

end function im_dielf_feg

!*************************************************************************

double complex function dielf_feg(omega)
implicit none
double precision :: omega
double precision :: qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi
double precision :: eps1, eps2
double precision :: re_dielf_feg,im_dielf_feg
external re_dielf_feg,im_dielf_feg
common /feg/ qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

  eps1 = re_dielf_feg(omega)
  eps2 = im_dielf_feg(omega)+broadening
  dielf_feg = eps1 + (0.d0,1.d0)*eps2

end function dielf_feg


!*************************************************************************

double precision function dre_dielf_feg(omega)
implicit none
double precision :: omega
double precision :: qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi
double precision :: lnterm1,lnterm2,xterm1,xterm2,xterm3,xterm4
double precision :: dlnterm1,dlnterm2
common /feg/ qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

  lnterm1 = log((Eq+qkf+omega)/(Eq-qkf+omega))
  lnterm2 = log((Eq+qkf-omega)/(Eq+qkf+omega))
  dlnterm1 = 1.d0/(Eq+qkf+omega)-1.d0/(Eq-qkf+omega)
  dlnterm2 = 1.d0/(Eq+qkf-omega)-1.d0/(Eq+qkf+omega)
  xterm1 = (4*Ef*Eq-(Eq+omega)**2)*dlnterm1/(2*qkf*(qq**2))
  xterm2 = (4*Ef*Eq-(Eq-omega)**2)*dlnterm2/(2*qkf*(qq**2))
  xterm3 = -(Eq+omega)*lnterm1/(2*qkf*(qq**2))
  xterm4 = (omega-Eq)*lnterm2/(2*qkf*(qq**2))
  dre_dielf_feg = 0.5d0*(qtf/qq)**2*(xterm1+xterm2+xterm3+xterm4)

end function dre_dielf_feg

!*************************************************************************

subroutine lossfn_free_electron_gas(omegahi,omegalo,omegap,qq0,broadening0,Ef0,qtf0,xkf0,edens0,lossfn_feg)
implicit none
double precision :: omegahi,omegalo,omegap
double precision :: qq0,Ef0,qtf0,xkf0,edens0,broadening0
double precision :: qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi
double precision :: abr,rlr,error
integer :: nsing,numcal,maxns
double precision :: xsing(20),omegaq,wtest
double complex :: lossfn_feg
double complex :: dielf_feg,cgrater
external dielf_feg,cgrater
common /feg/ qq,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

  qq = qq0
  Ef = Ef0
  qtf = qtf0
  xkf = xkf0
  edens = edens0
  Eq = 0.5d0*qq*qq
  qkf = qq*xkf
  broadening = broadening0
  omegaq = omegap + 0.3*qkf*qkf/omegap

  abr = 1.d-7
  rlr = 1.d-7

  nsing = 0
  wtest = omegaq-100.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq-10.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq-3.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq-1.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq+1.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq+3.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq+10.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif
  wtest = omegaq+100.*broadening
  if ((omegahi.gt.wtest).and.(wtest.gt.omegalo)) then
    nsing = nsing+1
    xsing(nsing) = wtest
  endif

  lossfn_feg = cgrater(dielf_feg,omegalo,omegahi,abr,rlr,nsing,xsing,error,numcal,maxns)/(omegahi-omegalo)

end subroutine lossfn_free_electron_gas


