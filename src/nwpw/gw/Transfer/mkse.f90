subroutine mkse(iband,ikk,lossfn, &
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
& nband,nsppol,shiftk,zz,cse,xse)
implicit none
integer :: iband,ikk(3),nwpt,nbcore,nbocc,ncband,ngkpt(3),natom,nlmn
integer :: igmn(3),igmx(3)
integer :: npwt,bantot,ncg,nkpt,nqpt,nsym,nsppol
integer :: npw,npwc,npwx,ipw1,npwndx,ntpwndx,napwndx,ipaw,itetrahedron
double precision :: vol,pi,wmax,xred(3,natom)
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
integer :: igg(3),igh(3),igg0(3),isym,isymq,ibp,iqpt,iqsym,iw,iww,ie,je,ies,iqctr
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
double precision :: vq(ngkpt(1),ngkpt(2),ngkpt(3)),qq(3),qp(3),qq2,qp2,qs(3),ek,ekmq
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
double precision :: xunit(3,3)

do ii=1,3
do jj=1,3
  if (ii.eq.jj) then
    xunit(ii,jj)=1.d0
  else
    xunit(ii,jj)=0.d0
  endif
enddo
enddo

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
igg0=(/0,0,0/)
do ie=-nwpt,nwpt
  cse1(ie)=(0.d0,0.d0)
enddo
xse=0.d0
dw=wmax/dble(nwpt)
do iw=1,nwpt
  wgrid(iw)=iw*dw
enddo
volelmnt=(vol*ngkpt(1)*ngkpt(2)*ngkpt(3))/((2.d0*pi)**3)
volred = 1.d0
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
nsing=0

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
do ibp=test_bands_se(1),test_bands_se(2)
!do ibp=1,1
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
!  do iipw=1,1 
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
      qp=qq+ipwx(:,ipw2)
      qq=qq+ipwx(:,ipw1)
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
      omega(ix,iy,iz)=enrgy(indxkbnd(ikptq)+ibp)-ek+dw/2
      if (iqq(1).eq.0.and.iqq(2).eq.0.and.iqq(3).eq.0) then
        if (ipw1.eq.1.and.ipw2.eq.1.and.iband.eq.ibp) then
! Find singular factor
!          call mkmatelX1(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
!&           ncg,nkpt,npwt,igmx,igmn,igndx, &
!&           isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
!&           lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
!&           cg,indxkcg,indxkpw,npwarr,kg, &
!&           cmatel)
!          if (ipaw.ne.0) then
!            call mkPAWmatelX1(pi,ibp,iband,ikk,jkk,qq,igg0, &
!&                 pwmatel,tpwmatel,nbcore,ncband,amatel)
!          else
!            amatel2=(0.d0,0.d0)
!          endif
!          smat2=4.d0*pi*dble((cmatel+amatel)*conjg(cmatel+amatel))
!          if (smat2.ne.0.d0) then
!            lqsing=.true.
!          else
!            lqsing=.false.
!          endif
          smat2=4.d0*pi
          lqsing=.true.
        else
          smat2=0.d0
          lqsing=.false.
        endif
! If not singular, find tensor response at q=0
        if (.not.lqsing) then
          if (ipwx(1,ipw1).eq.0.and.ipwx(2,ipw1).eq.0.and.ipwx(3,ipw1).eq.0) then
            call mkmatelJ1(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
&             ncg,nkpt,npwt,igmx,igmn,igndx, &
&             isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&             lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&             cg,indxkcg,indxkpw,npwarr,kg, &
&             jmatel)
            if (ipaw.ne.0) then
              call mkPAWmatelJ1(pi,ibp,iband,ikk,jkk, &
&                           pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
              jmatel(:) = jmatel(:) + vmatel(:)
            endif
            do ii=1,3
              jmatel(ii) = jmatel(ii)/(enrgy(indxkbnd(ikptq)+ibp)-ek)
            enddo
          else
            call mkmatelX1(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
&             ncg,nkpt,npwt,igmx,igmn,igndx, &
&             isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&             lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&             cg,indxkcg,indxkpw,npwarr,kg, &
&             cmatel)
            if (ipaw.ne.0) then
              call mkPAWmatelX1(pi,ibp,iband,ikk,jkk,qq,igg0, &
&                   pwmatel,tpwmatel,nbcore,ncband,amatel)
            else
              amatel=(0.d0,0.d0)
            endif
            jmatel = (cmatel+amatel)*qq/qq2
          endif
          if (ipw1.eq.ipw2) then
            jmatel2=jmatel
          else
            if (ipwx(1,ipw2).eq.0.and.ipwx(2,ipw2).eq.0.and.ipwx(3,ipw2).eq.0) then
              call mkmatelJ1(ibp,iband,ikpt,ikptq,ipwx(:,ipw2),igh, &
&               ncg,nkpt,npwt,igmx,igmn,igndx, &
&               isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&               lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&               cg,indxkcg,indxkpw,npwarr,kg, &
&               jmatel2)
              if (ipaw.ne.0) then
                call mkPAWmatelJ1(pi,ibp,iband,ikk,jkk, &
&                             pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
                jmatel2(:) = jmatel2(:) + vmatel(:)
              endif
              do ii=1,3
                jmatel2(ii) = jmatel2(ii)/(enrgy(indxkbnd(ikptq)+ibp)-ek)
              enddo
            else
              call mkmatelX1(ibp,iband,ikpt,ikptq,ipwx(:,ipw2),igh, &
&             ncg,nkpt,npwt,igmx,igmn,igndx, &
&             isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&             lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&             cg,indxkcg,indxkpw,npwarr,kg, &
&             cmatel2)
              if (ipaw.ne.0) then
                call mkPAWmatelX1(pi,ibp,iband,ikk,jkk,qp,igg0, &
&                     pwmatel,tpwmatel,nbcore,ncband,amatel2)
              else
                amatel2=(0.d0,0.d0)
              endif
              jmatel2 = (cmatel2+amatel2)*qp/qp2
            endif
          endif
!! DEBUG
!write(6,*) "jmatel"
!do ii=1,3
!write(6,*) jmatel(ii)
!enddo
!write(6,*) "jmatel2"
!do ii=1,3
!write(6,*) jmatel2(ii)
!enddo
!read(*,*)
          do ii=1,3
          do jj=1,3
            jmat2(ii,jj)=4.d0*pi*jmatel(ii)*conjg(jmatel2(jj))
          enddo
          enddo
        else
          do ii=1,3
          do jj=1,3
            jmat2(ii,jj)=(0.d0,0.d0)
          enddo
          enddo
        endif
      else
        call mkmatelX1(ibp,iband,ikpt,ikptq,ipwx(:,ipw1),igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatel)
        if (ipaw.ne.0) then
          call mkPAWmatelX1(pi,ibp,iband,ikk,jkk,qq,igg0, &
&               pwmatel,tpwmatel,nbcore,ncband,amatel)
        else
          amatel=(0.d0,0.d0)
        endif
        if (ipw1.eq.ipw2) then
          cmatel2=cmatel
          amatel2=amatel
        else
          call mkmatelX1(ibp,iband,ikpt,ikptq,ipwx(:,ipw2),igh, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatel2)
          if (ipaw.ne.0) then
            call mkPAWmatelX1(pi,ibp,iband,ikk,jkk,qp,igg0, &
&                 pwmatel,tpwmatel,nbcore,ncband,amatel2)
          else
            amatel2=(0.d0,0.d0)
          endif
        endif
        xmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)
        vqmat2(ix,iy,iz)=xmat2(ix,iy,iz)*vq(ix,iy,iz)
      endif
    enddo
    enddo
    enddo

    
! Alternate way of getting jmat2
! Interpolation of density matrix elements / q^2 to q=0
! For equally spaced samples A(-2) = -2 q, A(-1) = -1 q, A(1) = 1 q, A(2) = 2 q
! average of two quadratic interpolations gives (see Numerical Recipes 3.1)
! A(0) = (2/3) A(-1) + (2/3) A(1) - (1/6) A(-2) - (1/6) A(2)
! Expect vqmat2 to converge to well defined value at q = 0, interpolate this
! sum(ii,jj) qq(ii) * (bmet * jmat2 * bmet)(ii,jj) * qq(jj) / qq^2 = vqmat
! define rjmat2 = bmet * jmat2 * bmet
! find rjmat2 as q->0, then
! jmat2(q->0) = bmet^(-1)*rjmat2(q->0)*bmet^(-1)
!
! note qq * rjmat2 * qq = rjmat2(1,1)*qq(1)^2 + rjmat2(2,2)*qq(2)^2 + rjmat2(3,3)*qq(3)^2 + (rjmat2(1,2)+rjmat2(2,1))*qq(1)*qq(2) + (rjmat2(1,3)+rjmat2(3,1))*qq(1)*qq(3) + (rjmat2(2,3)+rjmat2(3,2))*qq(2)*qq(3)
! so longitudinal response only depends on sum of symmetric off-diagonal terms,
! independent of difference.
! Thus, we can safely set rjmat2(ii,jj) = (rjmat2(ii,jj)+rjmat2(jj,ii))/2
    if (.true.) then
! first do diagonal parts
! e.g., along direction 1
! qq = qq(1) * blat(1,:)
! qq^2 = qq(1)*qq(1) * dot_product(blat(1,:),blat(1,:))
! sum(ii,jj) qq(ii) * rjmat2(ii,jj) * qq(jj) / qq^2 
! = (1/qq^2) * rjmat2(1,1) * qq(1)*qq(1)
! = rjmat2(1,1) / dot_product(blat(1,:),blat(1,:))
! -> rjmat2(1,1) = dot_product(blat(1,:),blat(1,:)) * [response in direction 1]
      do ii=1,3
        ixx = (/0,0,0/)
        ixx(ii) = 1
        ixp2 = ictr+2*ixx
        ixp1 = ictr+ixx
        ixm1 = ictr-ixx
        ixm2 = ictr-2*ixx
        qq2 = dot_product(blat(ii,:),blat(ii,:))
        rjmatdum = 2.d0/3.d0 * vqmat2(ixp1(1),ixp1(2),ixp1(3)) &
&                + 2.d0/3.d0 * vqmat2(ixm1(1),ixm1(2),ixm1(3)) &
&                - 1.d0/6.d0 * vqmat2(ixp2(1),ixp2(2),ixp2(3)) &
&                - 1.d0/6.d0 * vqmat2(ixm2(1),ixm2(2),ixm2(3))
        rjmat2(ii,ii) =  qq2 * rjmatdum 
      enddo
! find off-diagonals using qq along directions (1,1,0), (1,0,1), and (0,1,1)
! e.g., along direction (1,1,0)
! qq = qq(1) * blat(1,:) + qq(2) * blat(2,:)
! qq^2 =    qq(1)*qq(1) * dot_product(blat(1,:),blat(1,:))
!       +   qq(2)*qq(2) * dot_product(blat(2,:),blat(2,:))
!       + 2*qq(1)*qq(2) * dot_product(blat(1,:),blat(2,:))
!      =    dot_product(blat(1,:),blat(1,:)) + dot_product(blat(2,:),blat(2,:))
!       + 2*dot_product(blat(1,:),blat(2,:))
! sum(ii,jj) qq(ii) * rjmat2(ii,jj) * qq(jj) / qq^2 
! = (1/qq^2) * ( rjmat2(1,1) * qq(1)*qq(1) + rjmat2(2,2) * qq(2)*qq(2)
!               + qq(1)*qq(2) * (rjmat2(1,2) + rjmat2(2,1)) )
! = (1/qq^2) * ( rjmat2(1,1) * qq(1)*qq(1) + rjmat2(2,2) * qq(2)*qq(2)
!               + 2 * rjmat2(1,2) qq(1)*qq(2))
! = (1/qq^2) * ( rjmat2(1,1) + rjmat2(2,2) + 2 * rjmat2(1,2) )
! -> rjmat2(1,2) = (qq^2 * [response in direction 1] - (rjmat2(1,1) + rjmat2(2,2)))/2
      do ii=1,3
        ixx = (/1,1,1/)
        ixx(ii) = 0
        ixp2 = ictr+2*ixx
        ixp1 = ictr+ixx
        ixm1 = ictr-ixx
        ixm2 = ictr-2*ixx
        qq2 = dot_product(blat(jj,:),blat(jj,:))   &
&           + dot_product(blat(kk,:),blat(kk,:))   &
&           + dot_product(blat(jj,:),blat(kk,:))*2
        rjmatdum = 2.d0/3.d0 * vqmat2(ixp1(1),ixp1(2),ixp1(3)) &
&                + 2.d0/3.d0 * vqmat2(ixm1(1),ixm1(2),ixm1(3)) &
&                - 1.d0/6.d0 * vqmat2(ixp2(1),ixp2(2),ixp2(3)) &
&                - 1.d0/6.d0 * vqmat2(ixm2(1),ixm2(2),ixm2(3))
        jj = mod(ii,3)+1   ! note ii = mod(ii-1,3)+1
        kk = mod(ii+1,3)+1
        rjmat2(jj,kk) = (qq2*rjmatdum - (rjmat2(jj,jj)+rjmat2(kk,kk)))/2.d0 
        rjmat2(kk,jj) = rjmat2(jj,kk)
!! DEBUG
!qq = ixx/dble(ngkpt)
!qq2=0.d0
!do kk=1,3
!do jj=1,3
!qq2=qq2+qq(kk)*bmet(kk,jj)*qq(jj)
!enddo
!enddo
!cvecdum = bmet(:,1)*qq(1)+bmet(:,2)*qq(2)+bmet(:,3)*qq(3)
!cvecdum2 = jmat2(:,1)*cvecdum(1)+jmat2(:,2)*cvecdum(2)+jmat2(:,3)*cvecdum(3)
!cvecdum = bmet(:,1)*cvecdum2(1)+bmet(:,2)*cvecdum2(2)+bmet(:,3)*cvecdum2(3)
!write(6,*) "qq: ",qq
!write(6,*) (cvecdum(1)*qq(1)+cvecdum(2)*qq(2)+cvecdum(3)*qq(3))/qq2
!write(6,*) ixx
!write(6,*) vqmat2(ixp2(1),ixp2(2),ixp2(3))
!write(6,*) vqmat2(ixp1(1),ixp1(2),ixp1(3))
!write(6,*) rjmatdum
!write(6,*) vqmat2(ixm1(1),ixm1(2),ixm1(3))
!write(6,*) vqmat2(ixm2(1),ixm2(2),ixm2(3))
      enddo
! DEBUG
!do ii=1,3
!do jj=1,3
!if (ii.eq.jj) then 
!rjmat2(ii,jj) = 2.9357115896d0
!else
!rjmat2(ii,jj) = 0.d0;
!endif
!enddo
!enddo
      bmet_t=bmet
      call inverse(bmet_t,bmetinv,3,3)
      do ii=1,3
      do jj=1,3
        jmat2(ii,jj) = (0.d0,0.d0)
        do kk=1,3
        do ll=1,3
          jmat2(ii,jj) = jmat2(ii,jj) + bmetinv(ii,kk)*rjmat2(kk,ll)*bmetinv(ll,jj)
        enddo
        enddo
      enddo
      enddo
    endif
!! DEBUG
!write(6,*) "rjmat2"
!do ii=1,3
!write(6,*) dble(rjmat2(:,ii))
!enddo
!write(6,*)
!do ii=1,3
!write(6,*) dimag(rjmat2(:,ii))
!enddo
!write(6,*)
!write(6,*) "jmat2"
!do ii=1,3
!write(6,*) dble(jmat2(:,ii))
!enddo
!write(6,*)
!do ii=1,3
!write(6,*) dimag(jmat2(:,ii))
!enddo
!write(6,*)
!write(6,*) "in Cartesian coordinates"
!do ii=1,3
!do jj=1,3
!ctensdum(ii,jj) = (0.d0,0.d0)
!do kk=1,3
!do ll=1,3
!ctensdum(ii,jj) = ctensdum(ii,jj) + blat(kk,ii)*jmat2(kk,ll)*blat(ll,jj)
!enddo
!enddo
!enddo
!enddo
!do ii=1,3
!write(6,*) dble(ctensdum(:,ii))
!enddo
!write(6,*)
!do ii=1,3
!write(6,*) dimag(ctensdum(:,ii))
!enddo
!stop
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
            ixx=(/ix,iy,iz/)
            iqpt=iqndx(ixx(1),ixx(2),ixx(3))
            call locateelement(ixx,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
            call locateelement(ixx,ipw2,ipw1,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpwt) ! transpose matrix coordiantes, to find anti-Hermitian part of lossfn
            if (iqpt.eq.iqctr) then
              vfactor(ixx(1),ixx(2),ixx(3),iw)=(0.d0,0.d0)
            else if (jjpw.ne.0.and.jjpwt.ne.0) then
              vfactor(ixx(1),ixx(2),ixx(3),iw)= &
&                 vqmat2(ixx(1),ixx(2),ixx(3))* &
&                 (lossfn(iw,iqpt,jjpw)-conjg(lossfn(iw,iqpt,jjpwt)))*0.5d0
            else
              vfactor(ixx(1),ixx(2),ixx(3),iw)=(0.d0,0.d0)
            endif
          enddo
          enddo
          enddo
!write(6,*) "iw = ", iw
          call locateelement(ictr,ipw1,ipw2,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpw)
          call locateelement(ictr,ipw2,ipw1,ngkpt,iqsymndx,npwc,npwx,nsym,pwsymndx,invpw2ndx,jjpwt)
          do ii=1,3
          do jj=1,3
            if (jjpw.ne.0.and.jjpwt.ne.0) then
              vjfactor(ii,jj,iw)=(0.d0,0.d0)
              do kk=1,3
              do ll=1,3
                iiq=ii + 3*(kk-1)
                vjfactor(ii,jj,iw)= &
&                     vjfactor(ii,jj,iw) + &
&                     0.5d0*(lossfn(iw,nqpt+iiq,jjpw)-conjg(lossfn(iw,nqpt+iiq,jjpwt))) &
&                     *bmet(kk,ll)*jmat2(ll,jj)
              enddo
              enddo
              iiq=ii + 3*(jj-1)
              vsfactor(ii,jj,iw)=smat2*(lossfn(iw,nqpt+iiq,jjpw)-conjg(lossfn(iw,nqpt+iiq,jjpwt)))*0.5d0
            else
              vjfactor(ii,jj,iw)=(0.d0,0.d0)
              vsfactor(ii,jj,iw)=(0.d0,0.d0)
            endif
          enddo
          enddo
!do ii=1,3
!write(6,*) dble(vjfactor(ii,:,iw))
!enddo
!write(6,*)
!do ii=1,3
!write(6,*) dimag(vjfactor(ii,:,iw))
!enddo
!write(6,*)
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
                if (lx) then
                  xpyr0(iv) = 0.d0
                endif
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
! DEBUG
!cdum = 0.d0
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
!if (ie.eq.0) cdum = cdum+vsm2(ii,jj,iw)*stens(jj,ii)
                    enddo
                    enddo
                  endif
                enddo
              enddo
!write(6,*) cdum
            endif
            if (lx) then
              wwwlo = min(evtx(1),evtx(2),evtx(3),evtx(4))
              dw = max(evtx(1),evtx(2),evtx(3),evtx(4))-wwwlo
!xdum = time()
              call singular_adaptive_tetint(qqr,evtx,wwwlo,dw,vcenterint,bmet,abr,rlr,rlr0,stens)
!! DEBUG
!cdum = 0.d0
              do ii=1,3
              do jj=1,3
                ssx=ssx+smat2*bmet(ii,jj)*stens(jj,ii)
!! DEBUG
!cdum = cdum+smat2*bmet(ii,jj)*stens(jj,ii)
              enddo
              enddo
!write(6,*) cdum
            endif
!write(6,*) 'Khst7',xdum
!write(32,'(4i3,2x,f14.10,2x,f14.10,"  sing")') ixx,itet,(cdum/2.d0)
!write(6,*) "finished singular integration"
          elseif (centerint) then
! anisotropic part
! DEBUG
!cdum = 0.d0
            if (lc) then
              do iww=iwwlim(1)-1,iwwlim(4)+1
                wwwlo=dw*(iww-0.5d0)
!                wwwhi=dw*(iww+0.5d0)
                do iw=1,nwpt
                  call vnitetint(rrpyr,evtx0,wwwlo,dw,vjfactor(:,:,iw),bmet,vcenterint,cint)
                  ie=iww-ocsign*iw
                  if (ie.ge.ielo.and.ie.le.iehi) ssi(ie)=ssi(ie)+cint/volelmnt
! DEBUG
!if (ie.eq.0.or.ie.eq.1) cdum = cdum+cint/volelmnt
                enddo
              enddo
            endif
            if (lx) then
              call vnitetint(rrpyr,evtx0,evtx0(indxe(1)),evtx0(indxe(4))-evtx0(indxe(1)),jmat2,bmet,vcenterint,cint)
              ssx=ssx+cint/volelmnt
            endif
! DEBUG
!write(32,'(4i3,2x,f14.10,2x,f14.10,"  cent")') ixx,itet,(cdum/2.d0)
          endif
! can use basic tetrahedron integration for everything else
! DEBUG
!cdum = 0.d0
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
!if (ie.eq.0.or.ie.eq.1) cdum = cdum+(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
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
        enddo
! DEBUG
!write(32,'(3i3,2x,4f16.10)') ixx,ssi(0),ssi(1)
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
        enddo
      endif
!do je=ielo,iehi
!write(6,*) je,ssi(je)
!enddo
      deallocate (ssi)
!write(6,'("lstv",2x,f14.10,4x,"(",f14.10,",",f14.10,")")') ssx,(ssc(0)+ssc(1))/2.

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

!do ie=-nwpt,nwpt
!  eps1=dble(ie)*dw-dw/2
!  write(12,'(i4,f10.3,2(3x,2es12.3))') ie,eps1*27.2114,cse1(ie)*27.2114
!enddo
if (itetrahedron.eq.1) then
  cse=(cse1(0)+cse1(1))/2
  dcse=(cse1(1)-cse1(0))/dw
endif
zz=1.d0/(1.d0-dcse)
!write(6,*) cse*27.2114
!write(6,*) xse*27.2114
!write(6,*) zz
!write(6,*) dble(zz)*dimag(cse)*27.2114

return
end subroutine mkse

!****************************************************************************
! evaluate Integral_wmin^wmax dww tsingint

subroutine tsinggrater(iv,jv,wmin,wmax,ive,xk,xkvtx,ngkpt,evtx,abr,rlr,nsing,wsing,sint)
implicit none
integer :: iv,jv,numcal,maxns,mx,nstack,ii,jj,kk,icount,ive,nsing
double precision :: wsing(20)
parameter (mx=1500)
integer :: ngkpt(3)
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
integer :: i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
common /tsingint_DEBUG/ i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater) and the region
! is removed from consideration (and a new region is the right-most)
  nstack=nsing+1
  maxns=nstack
  error=0.d0
  sint=0.d0
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
  wleft(1)=wmin
  wleft(nsing+2)=wmax
  if(nsing.gt.0) then
    do jj=1,nsing
      wleft(jj+1)=wsing(jj)
    enddo
  endif
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=wleft(ii+1)-wleft(ii)
    do jj=1,3
      ww=wleft(ii)+del*dw(jj)
if (k_ct_DEBUG.eq.k_ref_DEBUG) i_ct_DEBUG = i_ct_DEBUG+1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "%%%%%AGS  ",ww,i_ct_DEBUG
      call tsingint(iv,jv,ive,ww,xk,xkvtx,ngkpt,evtx,fval(jj,ii))
    enddo
  enddo
  numcal=nstack*3
  do
    if(nstack+3.ge.mx) then
      write(6,*) 'tsinggrater: too many regions'
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
if (k_ct_DEBUG.eq.k_ref_DEBUG) i_ct_DEBUG = i_ct_DEBUG+1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "%%%%%AGS  ",ww,i_ct_DEBUG
      call tsingint(iv,jv,ive,ww,xk,xkvtx,ngkpt,evtx,fval(1,jj))
      ww=wleft(jj)+dw(3)*del1
if (k_ct_DEBUG.eq.k_ref_DEBUG) i_ct_DEBUG = i_ct_DEBUG+1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "%%%%%AGS  ",ww,i_ct_DEBUG
      call tsingint(iv,jv,ive,ww,xk,xkvtx,ngkpt,evtx,fval(3,jj))
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
end subroutine tsinggrater

!***************************************************************************
!
! Over a surface S defined by the triangle with vertices xkvtx+xkp(1:3,1)
!                                                        xkvtx+xkp(1:3,2)
!                                                        xkvtx+xkp(1:3,3)
! where xkp defined to give S of constant energy
! (see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972))
! and vector k measured from tip of xkvtx
! find      sint = integral dS q(iv)q(jv)/q^4
! where q=k+xkvtx
! On S, define coordinates X and Y such that k=xkp(1:3,3)+X*xq1+Y*xq2
! normal=cross product(xq1,xq2)
! integral dS -> |normal|*integral^1_0 dX integral^{1-X}_0 dY
! integral^{1-X}_0 dY ((k(iv)+xkvtx(iv))*(k(jv)+xkvtx(jv)))/(dot_product(k+xkvtx,k+xkvtx))**2 = tfn1
! numerically integrate tfn1
! sint = integral^1_0 dX tfn1

subroutine tsingint(iv,jv,ive,ww,xk,xkvtx,ngkpt,evtx,sint)
use geometry
implicit none
integer :: iv,jv,ive,ivv,jvv
integer :: ngkpt(3)
double precision :: xk(3,4),xkvtx(3),evtx(4),ww
double precision :: xkp(3,3),xq1(3),xq2(3),xkv0(3),xcvt,xkvmag(3)
double precision :: xkv1(3),xkv2(3)
double precision :: vt1(3),vt2(3),vt3(3)
double precision :: aa0,aa1,aa2,bb0,bb1,cc,dd0,dd1,dd2,ee0,ee1,ff,qq,delta,brd
double precision :: xmin,xmax,abr,rlr,xsing(30),error
double precision :: tfn1,xfn2,xfn3,grater
double precision :: sint,sint1,sint2,sint3
double precision :: root1,root2,root1a,root2a,xroot1,xroot2,xroot1a,xroot2a
double precision :: xnorm(3),area
double precision :: xdum,xmx,rlxk
double precision :: vtemp(3),utemp(3),wtemp(3),vmag2,umag2,wmag2
double precision polyterm(0:4)
integer :: ii,jj,kk,ll
integer :: nsing,numcal,maxns,nroots,nrootsa
integer :: iter 
integer :: indxq(3) ! sorts vertices in ascending magnitude of q-vector
double precision :: xd01,xd02  ! xx value where delta is 0
double precision :: xdt11,xdt12,xdt21,xdt22  ! for xx = xd0 + dx, linear and quadratic term in dx for delta
double precision :: xx_dd1, xx_dd2 ! roots of dd = 0
double precision :: xx_def1, xx_def2 ! roots of dd + ee*zz + ff*zz^2 = 0
logical :: usexd1,usexd2,usexx_dd1,usexx_dd2 ! are roots in integration region
double precision :: denom,snangle,snangle1,snangle2,snangle3
double precision :: aat1,aat2,bbt1,bbt2,ddt1,ddt2,eet1,eet2,zzt1,zzt2,deltat1,deltat2,aa2ddt1,aa2ddt2
double precision :: RRA1,RRA2,RRB1,RRB2,RRC1,RRC2,d0e0f,d1e1e02f,d2e1f
double precision :: qqmag(2), q02k1, q22qdiff, qqdiff(3), qqdiffmag
double precision :: xi1,xi2,zeta12,eta0,eta2,nu2,kappa,lambda
double precision :: omxi12,omxi22,ometa22,omzeta122
double precision :: deltaprefact, d_angle
double precision :: rrootd,irootd,rrootdelta,irootdelta,rrootdezfzz,irootdezfzz
double precision :: gammapfroot1,gammapfroot2,gammapfpf
integer :: gammapfnroot
double precision :: denomroot1,denomroot2,denomroot3,denomroot4,denompf
integer :: denomnroot
common /tfn/ aa0,aa1,aa2,bb0,bb1,cc,dd0,dd1,dd2,ee0,ee1,ff,brd, &
& xkv0,xq1,ivv,jvv, &
& xd01,xd02,xdt11,xdt12,xdt21,xdt22,xx_dd1,xx_dd2,xx_def1,xx_def2, &
& usexd1,usexd2,usexx_dd1,usexx_dd2, &
& aat1,aat2,bbt1,bbt2,ddt1,ddt2,eet1,eet2,zzt1,zzt2, &
& deltat1,deltat2,aa2ddt1,aa2ddt2, qqmag, qqdiffmag, &
& q02k1, q22qdiff, deltaprefact, umag2, &
& rrootd,irootd,rrootdelta,irootdelta,rrootdezfzz,irootdezfzz, &
& gammapfroot1,gammapfroot2,denomroot1,denomroot2,denomroot3,denomroot4, &
& gammapfnroot,denomnroot, &
& gammapfpf,denompf, &
& RRA1,RRA2,RRB1,RRB2,RRC1,RRC2,d0e0f, &
& xi1,xi2,zeta12,eta0,eta2,kappa,lambda,omxi12,omxi22,ometa22,omzeta122
integer :: i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
common /tsingint_DEBUG/ i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
double precision :: sqrtm1
external tfn1,xfn2,xfn3,grater,sqrtm1
double precision :: xdummy, ydummy, zdummy,xxdummy,yydummy,zzdummy
integer :: idummy, jdummy, kdummy
double complex :: cdummy1,cdummy2

! broadening: 1/|q|^2 -> 1/(|q|^2+delta^2)
!  delta=min(dot_product(xk(:,2),xk(:,2)),dot_product(xk(:,3),xk(:,3)),dot_product(xk(:,4),xk(:,4)))*1.d-4
  ivv=iv
  jvv=jv
  if (ive.eq.1) then
    xkp(:,1)=(ww-evtx(1))*xk(:,2)/(evtx(2)-evtx(1))
    xkp(:,2)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
    xkp(:,3)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  elseif (ive.eq.2) then
    xkp(:,1)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
    xkp(:,2)=(ww-evtx(2))*(xk(:,3)-xk(:,2))/(evtx(3)-evtx(2)) + xk(:,2)
    xkp(:,3)=(ww-evtx(1))*(xk(:,4)-xk(:,1))/(evtx(4)-evtx(1))
  elseif (ive.eq.3) then
    xkp(:,1)=(ww-evtx(1))*xk(:,3)/(evtx(3)-evtx(1))
    xkp(:,2)=(ww-evtx(2))*(xk(:,3)-xk(:,2))/(evtx(3)-evtx(2)) + xk(:,2)
    xkp(:,3)=(ww-evtx(1))*(xk(:,4)-xk(:,1))/(evtx(4)-evtx(1))
  else
    xkp(:,1)=(ww-evtx(2))*(xk(:,4)-xk(:,2))/(evtx(4)-evtx(2)) + xk(:,2)
    xkp(:,2)=(ww-evtx(3))*(xk(:,4)-xk(:,3))/(evtx(4)-evtx(3)) + xk(:,3)
    xkp(:,3)=(ww-evtx(1))*xk(:,4)/(evtx(4)-evtx(1))
  endif
  do kk=1,3
    xkvmag(kk) = 0.d0
    do ii=1,3
    do jj=1,3
      xkvmag(kk) = xkvmag(kk) + ((xkp(ii,kk)+xkvtx(ii))/ngkpt(ii))*bmet(ii,jj)*((xkp(jj,kk)+xkvtx(jj))/ngkpt(jj))
    enddo
    enddo
  enddo
  call indxhpsort(3,3,xkvmag,indxq)
  xq1=(xkp(:,1)-xkp(:,3))/ngkpt
  xq2=(xkp(:,2)-xkp(:,3))/ngkpt
  xkv0=(xkp(:,3)+xkvtx)/ngkpt
  xkv1=(xkp(:,1)+xkvtx)/ngkpt
  xkv2=(xkp(:,2)+xkvtx)/ngkpt
  qqdiff=xkv1-xkv2
  vt1=xkv0+xq2
  vt2=xq1-xq2
  call crossreduced(xq1,xq2,xnorm)
  area=sqrt(dot_product(xnorm,xnorm))
  aa0=xkv0(iv)*xkv0(jv)
  aa1=xq1(iv)*xkv0(jv)+xq1(jv)*xkv0(iv)
  aa2=xq1(iv)*xq1(jv)
  bb0=xq2(iv)*xkv0(jv)+xq2(jv)*xkv0(iv)
  bb1=xq1(iv)*xq2(jv)+xq1(jv)*xq2(iv)
  cc=xq2(iv)*xq2(jv)
  dd0 = 0.d0
  dd1 = 0.d0
  dd2 = 0.d0
  ee0 = 0.d0
  ee1 = 0.d0
  ff  = 0.d0
  d0e0f = 0.d0
  d1e1e02f = 0.d0
  d2e1f = 0.d0
  eta0 = 0.d0
  eta2 = 0.d0
  nu2 = 0.d0
  qqmag(1) = 0.d0
  qqmag(2) = 0.d0
  qqdiffmag = 0.d0
  do ii=1,3
  do jj=1,3
    dd0 = dd0 + xkv0(ii)*bmet(ii,jj)*xkv0(jj)
    dd1 = dd1 +  xq1(ii)*bmet(ii,jj)*xkv0(jj)*2
    dd2 = dd2 +  xq1(ii)*bmet(ii,jj)*xq1(jj)
    ee0 = ee0 +  xq2(ii)*bmet(ii,jj)*xkv0(jj)*2
    ee1 = ee1 +  xq1(ii)*bmet(ii,jj)*xq2(jj)*2
    ff  = ff  +  xq2(ii)*bmet(ii,jj)*xq2(jj)
    d0e0f = d0e0f + vt1(ii)*bmet(ii,jj)*vt1(jj)
    d1e1e02f = d1e1e02f + vt1(ii)*bmet(ii,jj)*vt2(jj)*2
    d2e1f = d2e1f + vt2(ii)*bmet(ii,jj)*vt2(jj)
    eta0 = eta0 + qqdiff(ii)*bmet(ii,jj)*xkv0(jj)
    eta2 = eta2 + qqdiff(ii)*bmet(ii,jj)*xkv2(jj)
    nu2 = nu2 + qqdiff(ii)*bmet(ii,jj)*xq2(jj)
    qqmag(1) = qqmag(1) + xkv1(ii)*bmet(ii,jj)*xkv1(jj)
    qqmag(2) = qqmag(2) + xkv2(ii)*bmet(ii,jj)*xkv2(jj)
    qqdiffmag = qqdiffmag + qqdiff(ii)*bmet(ii,jj)*qqdiff(jj)
  enddo
  enddo
  xi1 = dd1/(2*sqrt(dd0*dd2))
  xi2 = ee0/(2*sqrt(dd0*ff))
  zeta12 = ee1/(2*sqrt(dd2*ff))
  eta0 = eta0/sqrt(qqdiffmag*dd0)
  eta2 = eta2/sqrt(qqdiffmag*qqmag(2))
  nu2 = nu2/sqrt(qqdiffmag*ff)
  call sinreducedangle(xkv0,xq1,omxi12)
  call sinreducedangle(xkv0,xq2,omxi22)
  call sinreducedangle(xq1,xq2,omzeta122)
  call sinreducedangle(xkv2,qqdiff,ometa22)
  q02k1 = sqrt(dd0/dd2)
  q22qdiff = sqrt(qqmag(2)/qqdiffmag)
  rrootd = -q02k1*xi1
  irootd = q02k1*omxi12 != q02k1*sqrt(1-xi1**2)
  call crossreduced(xq2,xkv0,utemp)
  call crossreduced(xq2,xq1,vtemp)
  utemp = 2*utemp
  vtemp = 2*vtemp
  call cross(utemp,vtemp,wtemp)
  vmag2 = dot_product(vtemp,vtemp)
  umag2 = dot_product(utemp,utemp)
  wmag2 = dot_product(wtemp,wtemp)
  deltaprefact = vmag2 ! = 4*dd2*ff-ee1*ee1
  if (vmag2.eq.0.d0) then
    rrootdelta = 0.d0
    irootdelta = 0.d0
  else
    rrootdelta = -dot_product(utemp,vtemp)/vmag2
    irootdelta = sqrt(wmag2)/vmag2
  endif
!  rrootdelta = -q02k1*kappa
!  irootdelta = q02k1*lambda
  rrootdezfzz = -q22qdiff*eta2
  irootdezfzz = q22qdiff*ometa22 != q22qdiff*sqrt(1-eta2**2)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "q1*q2:      ",dd0+0.5*(dd1+ee0+ee1)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "eta2:       ",eta2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "xi1:        ",xi1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "xi2:        ",xi2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "zeta12:     ",zeta12
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "q02k1:      ",q02k1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "dd2:        ",dd2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "ff:         ",ff
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "1-zeta12**2:",1-zeta12**2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "deltaprefact:      ",deltaprefact
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "umag2:             ",umag2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "vmag2:             ",vmag2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "wmag2:             ",wmag2
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "rroot (d):         ",rrootd
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "iroot (d):         ",irootd
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "rroot (delta):     ",rrootdelta, &
!& -q02k1*(xi1-zeta12*xi2)/(1-zeta12**2)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "iroot (delta):     ",irootdelta, &
!& q02k1*sqrt(1.d0 + 2*xi1*xi2*zeta12 - xi1*xi1 - xi2*xi2 - zeta12*zeta12)/(1-zeta12**2)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "rroot (d+ez+fz^2): ",rrootdezfzz
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "iroot (d+ez+fz^2): ",irootdezfzz
!write(77,*) xx, aa, aa2*((xx+q02k1*xi1)**2 + (1-xi1**2)*(q02k1)**2)
!  brd=max(abs(dd0),abs(dd1),abs(dd2),abs(ee0),abs(ee1),abs(ff))*1.e-10
  nsing=0
  root1 = rrootd-100*irootd
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D: ",xsing(nsing)
  endif
  root1 = rrootd-10*irootd
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D: ",xsing(nsing)
  endif
  root1 = rrootd
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D: ",xsing(nsing)
  endif
  root1 = rrootd+10*irootd
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D: ",xsing(nsing)
  endif
  root1 = rrootd+100*irootd
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D: ",xsing(nsing)
  endif
  root1 = rrootdelta-100*irootdelta
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "Delta: ",xsing(nsing)
  endif
  root1 = rrootdelta-10*irootdelta
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "Delta: ",xsing(nsing)
  endif
  root1 = rrootdelta
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "Delta: ",xsing(nsing)
  endif
  root1 = rrootdelta+10*irootdelta
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "Delta: ",xsing(nsing)
  endif
  root1 = rrootdelta+100*irootdelta
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "Delta: ",xsing(nsing)
  endif
  root1 = rrootdezfzz-100*irootdezfzz
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D+EZ+FZ^2: ",xsing(nsing)
  endif
  root1 = rrootdezfzz-10*irootdezfzz
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D+EZ+FZ^2: ",xsing(nsing)
  endif
  root1 = rrootdezfzz
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D+EZ+FZ^2: ",xsing(nsing)
  endif
  root1 = rrootdezfzz+10*irootdezfzz
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D+EZ+FZ^2: ",xsing(nsing)
  endif
  root1 = rrootdezfzz+100*irootdezfzz
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "D+EZ+FZ^2: ",xsing(nsing)
  endif
  call hpsort(nsing,nsing,xsing(1:nsing))
  xmin=0.d0
  xmax=1.d0
  ii = 1
  if (nsing.gt.0) then
    do
      root1 = xsing(nsing)*100
      if (root1.ge.xmax) then
        exit
      else
        nsing = nsing+1
        xsing(nsing)=root1
      endif
    enddo
  endif
  root1 = -ee0/ee1
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E: ",xsing(nsing)
  endif
  root1 = -(ee0+2*ff)/(ee1-2*ff)
  if (root1.gt.0.d0.and.root1.lt.1.d0) then
    nsing=nsing+1
    xsing(nsing)=root1
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E+2FZ: ",xsing(nsing)
  endif
  do
    if (ii.ge.nsing) exit
    if (xsing(ii).eq.xsing(ii+1)) then
      nsing = nsing-1
      do jj=ii+1,nsing
        xsing(jj) = xsing(jj+1)
      enddo
    endif
    ii = ii + 1
  enddo
  call hpsort(nsing,nsing,xsing(1:nsing))

  gammapfpf = 4*(cc*dd2+ff*aa2)-2*bb1*ee1
  call rquadroots(4*(cc*dd2+ff*aa2)-2*bb1*ee1, &
&                 4*(cc*dd1+ff*aa1)-2*(bb0*ee1+bb1*ee0), &
&                 4*(cc*dd0+ff*aa0)-2*bb0*ee0, &
&                 gammapfnroot,gammapfroot1,gammapfroot2)

  polyterm(4) = 2*ff*aa2*dd2 - aa2*ee1**2 + ff*aa2*ee1 + dd2*bb1*ee1 &
&             - 2*ff*dd2*bb1 + cc*dd2*ee1 - 2*cc*dd2*dd2
  polyterm(3) = 2*ff*(aa1*dd2+aa2*dd1) &
&             - (aa1*ee1**2 + 2*aa2*ee0*ee1) &
&             - ff*(aa2*(ee1-ee0)-aa1*ee1) &
&             + (dd2*(bb0*ee1+bb1*ee0)+dd1*bb1*ee1) &
&             + 2*ff*(dd2*(bb1-bb0)-dd1*bb1) &
&             - cc*(dd2*(ee1-ee0)-dd1*ee1) &
&             - 4*cc*dd1*dd2
  polyterm(2) = 2*ff*(aa2*dd0+aa1*dd1+aa0*dd2) &
&             - (aa0*ee1**2+2*aa1*ee0*ee1+aa2*ee0**2) &
&             - ff*(-aa0*ee1+aa1*(ee1-ee0)+aa2*ee0) &
&             + (dd0*bb1*ee1+dd1*(bb0*ee1+bb1*ee0)+dd2*bb0*ee0) &
&             + 2*ff*(-dd0*bb1+dd1*(bb1-bb0)+dd2*bb0) &
&             - cc*(-dd0*ee1+dd1*(ee1-ee0)+dd2*ee0) &
&             - 2*cc*(2*dd0*dd2+dd1**2)
  polyterm(1) = 2*ff*(aa0*dd1+aa1*dd0) &
&             - (2*aa0*ee0*ee1+aa1*ee0**2) &
&             - ff*(aa0*(ee1-ee0)+aa1*ee0) &
&             + (dd0*(bb0*ee1+bb1*ee0)+dd1*bb0*ee0) &
&             + 2*ff*(dd0*(bb1-bb0)+dd1*bb0) &
&             - cc*(dd0*(ee1-ee0)+dd1*ee0) &
&             - 4*cc*(dd0*dd1)
  polyterm(0) = 2*ff*aa0*dd0 - aa0*ee0**2 - ff*aa0*ee0 + dd0*bb0*ee0 &
&             + 2*ff*dd0*bb0 - cc*dd0*ee0 - 2*cc*dd0**2
  denompf = polyterm(4)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "AAA"
  call rquarticroots(polyterm(4),polyterm(3),polyterm(2),polyterm(1),polyterm(0), &
&                    denomnroot,denomroot1,denomroot2,denomroot3,denomroot4)
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "BBB"

  abr=1.d-12
  rlr=1.d-6
if (k_ct_DEBUG.eq.k_ref_DEBUG) then 
write(6,*) ww,i_ct_DEBUG
if (ww.ge.0.002125.and.ww.le.0.0021271) then
write(6,*) "nsing: ",nsing
do ii=1,nsing
write(6,*) xsing(ii)
enddo
write(6,*) "|q0|:              ",sqrt(dd0)
write(6,*) "d prefact:         ",dd2
write(6,*) "rroot (d):         ",rrootd
write(6,*) "iroot (d):         ",irootd
write(6,*) "delta prefact:     ",deltaprefact
write(6,*) "rroot (delta):     ",rrootdelta
write(6,*) "iroot (delta):     ",irootdelta
write(6,*) "(d+ez+fz^2) prfct: ",qqdiffmag
write(6,*) "rroot (d+ez+fz^2): ",rrootdezfzz
write(6,*) "iroot (d+ez+fz^2): ",irootdezfzz
write(6,*) "x^4 term: ",polyterm(4),denompf
write(6,*) "x^3 term: ",polyterm(3), &
& -denompf*(denomroot1+denomroot2+2*denomroot3)
write(6,*) "x^2 term: ",polyterm(2), &
& denompf*(denomroot1*denomroot2+2*denomroot1*denomroot3 &
&         +2*denomroot2*denomroot3+denomroot3*denomroot3+denomroot4*denomroot4)
write(6,*) "x^1 term: ",polyterm(1), &
& denompf*((denomroot1+denomroot2)*(denomroot3*denomroot3+denomroot4*denomroot4) &
&         +2*denomroot1*denomroot2*denomroot3)
write(6,*) "x^0 term: ",polyterm(0), denompf*(denomroot1*denomroot2*(denomroot3*denomroot3+denomroot4*denomroot4))
write(6,*) "gammapfpf: ",gammapfpf
write(6,*) "gammapfnroot: ",gammapfnroot
write(6,*) "roots: "
write(6,*) gammapfroot1
write(6,*) gammapfroot2
write(6,*) "denompf: ",denompf
write(6,*) "denomnroot: ",denomnroot
write(6,*) "roots: "
write(6,*) denomroot1
write(6,*) denomroot2
write(6,*) denomroot3
write(6,*) denomroot4
write(6,*)
l_ct_debug = l_ct_DEBUG + 1
i_prt_DEBUG = 1
else
i_prt_DEBUG = 0
endif
if (l_ct_DEBUG.ge.100) stop
!write(6,*)
!write(6,*) "q0^2: ",xkvmag(indxq(1))
!write(6,*) "q1^2: ",xkvmag(indxq(2))
!write(6,*) "q2^2: ",xkvmag(indxq(3))
!write(6,*) "|q0|: ",sqrt(dd0)
!write(6,*) "|k1|: ",sqrt(dd2)
!write(6,*) "|k2|: ",sqrt(ff)
!write(6,*) "q0: ",xkv0," (reduced)"
!write(6,*) "q1: ",(xkp(:,indxq(2))+xkvtx)/ngkpt," (reduced)"
!write(6,*) "q2: ",(xkp(:,indxq(3))+xkvtx)/ngkpt," (reduced)"
!write(6,*) "k1: ",xq1," (reduced)"
!write(6,*) "k2: ",xq2," (reduced)"
!write(6,*) "d: ",dd0,dd1,dd2
!write(6,*) "e: ",ee0,ee1
!write(6,*) "f: ",ff
!write(6,*) "d0+e0+f: ",dd0+ee0+ff,d0e0f
!write(6,*) "d1+e1-e0-2f: ",dd1+ee1-ee0-2*ff,d1e1e02f
!write(6,*) "d2-e1+f: ",dd2-ee1+ff,d2e1f
endif
  sint=grater(tfn1,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)*area
!!write(6,*) "%%%%%%%% ",ww,sint,i_ct_DEBUG
if (k_ct_DEBUG.eq.k_ref_DEBUG) then
write(78,*) ww,sint,i_ct_DEBUG
if (i_ct_DEBUG.gt.1000) stop
if (j_ct_DEBUG+1.gt.1) stop
endif
!if (k_ct_DEBUG.gt.k_ref_DEBUG) stop
end subroutine tsingint

!***************************************************************************

double precision function tfn1(xx)
implicit none
double precision :: xx
double precision :: aa0,aa1,aa2,bb0,bb1,cc,dd0,dd1,dd2,ee0,ee1,ff,brd
double precision :: xkv0(3),xq1(3)
double precision :: aa,bb,dd,ee
double precision :: xx2,zz,delta,sqdelta,xterm0,xterm1,xterm2,dterm1,dterm2, &
&   gammat,dgammaD,dgammaD_old
double precision :: xx0d,dxxd
double precision :: xterm1inv
double precision :: logterm1,logterm2,logterm3,logterm4,atf,darg,dargoe
double precision :: arg1,arg2,sum1,sum2,prod1,temp
double precision :: test1,test2,test3,test4,testlo1,testlo2,testlo3
double precision :: dtarg
integer :: xi
integer :: ii,jj,kk,ll
integer :: ivv,jvv
logical :: lowdarg
double precision :: xdiff1,xdiff2,dquad
double precision :: h1,h2,h3,h4,h5,h6,h7,h8,term1,term2,term3,term4
double precision :: aa2dd
double precision :: hiorderexp,quadterm1,quadterm2,quadterm3
double precision :: tt0,tt1,tt2,ta0,ta1,ta2,tx0,tx1,tx2,tc1,tc2
double precision :: xd01,xd02,xd0,dx  ! xd0 -> xx value where delta is 0
double precision :: xdt11,xdt12,xdt21,xdt22,xdt1,xdt2  ! for xx = xd0 + dx, linear and quadratic term in dx for delta
double precision :: xx_dd1, xx_dd2 ! roots of dd = 0
double precision :: xx_def1, xx_def2 ! roots of dd + ee*zz + ff*zz^2 = 0
logical :: usexd1,usexd2,usexx_dd1,usexx_dd2 ! are roots in integration region
double precision :: aat1,aat2,bbt1,bbt2,ddt1,ddt2,eet1,eet2,zzt1,zzt2,deltat1,deltat2,aa2ddt1,aa2ddt2
double precision :: RRA1,RRA2,RRB1,RRB2,RRC1,RRC2,d0e0f,d1e1e02f,d2e1f
double precision :: qqmag(2), q02k1, q22qdiff, qqdiffmag
double precision :: xi1,xi2,zeta12,eta0,eta2,nu2,kappa,lambda
double precision :: omxi12,omxi22,ometa22,omzeta122
double precision :: deltaprefact, umag2
double precision :: rrootd,irootd,rrootdelta,irootdelta,rrootdezfzz,irootdezfzz
double precision :: gammapfroot1,gammapfroot2,gammapfpf
integer :: gammapfnroot
double precision :: denomroot1,denomroot2,denomroot3,denomroot4,denompf
integer :: denomnroot
common /tfn/ aa0,aa1,aa2,bb0,bb1,cc,dd0,dd1,dd2,ee0,ee1,ff,brd, &
& xkv0,xq1,ivv,jvv, &
& xd01,xd02,xdt11,xdt12,xdt21,xdt22,xx_dd1,xx_dd2,xx_def1,xx_def2, &
& usexd1,usexd2,usexx_dd1,usexx_dd2, &
& aat1,aat2,bbt1,bbt2,ddt1,ddt2,eet1,eet2,zzt1,zzt2, &
& deltat1,deltat2,aa2ddt1,aa2ddt2, qqmag, qqdiffmag, &
& q02k1, q22qdiff, deltaprefact,umag2, &
& rrootd,irootd,rrootdelta,irootdelta,rrootdezfzz,irootdezfzz, &
& gammapfroot1,gammapfroot2,denomroot1,denomroot2,denomroot3,denomroot4, &
& gammapfnroot,denomnroot, &
& gammapfpf,denompf, &
& RRA1,RRA2,RRB1,RRB2,RRC1,RRC2,d0e0f, &
& xi1,xi2,zeta12,eta0,eta2,kappa,lambda,omxi12,omxi22,ometa22,omzeta122
double precision :: xdummy, ydummy, zdummy,xxdummy,yydummy,zzdummy
integer :: idummy,jdummy,kdummy,ldummy
integer :: i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG
common /tsingint_DEBUG/ i_prt_DEBUG, j_prt_DEBUG, k_prt_DEBUG, &
& i_ct_DEBUG, j_ct_DEBUG, k_ct_DEBUG, l_ct_DEBUG, &
& k_ref_DEBUG, i_ref_DEBUG

xdummy = 0.d0
ydummy = 0.d0
xxdummy = 0.d0
yydummy = 0.d0
  xx2 = xx*xx
  zz = 1.d0-xx
  xx0d = -q02k1*xi1
  dxxd = xx - xx0d
  aa = (xkv0(ivv)+xx*xq1(ivv))*(xkv0(jvv)+xx*xq1(jvv)) ! =aa0 + aa1*xx + aa2*xx2
  bb = bb0 + bb1*xx
  dd = dd2*((xx-rrootd)**2 + (irootd)**2) ! =dd0 + dd1*xx + dd2*xx2
  ee = ee0 + ee1*xx
  aa2dd = aa/dd
  if (dd.eq.0.d0) aa2dd = (aa1+(2*xx0d+dxxd)*aa2)/(dd1+(2*xx0d+dxxd)*dd2)
  xterm0 = zz*(ff*zz + ee)
  xterm1 = qqdiffmag*((xx-rrootdezfzz)**2 + (irootdezfzz)**2) ! = dd + xterm0
  xterm2 = xterm0/xterm1
  darg = 2*ff*zz
  dterm1 = ee+darg
  lowdarg = (abs(darg/ee).lt.1.d-4)
  if (deltaprefact.eq.0.d0) then
    delta = umag2
  else 
    delta = deltaprefact*((xx-rrootdelta)**2 + (irootdelta)**2) ! =4*dd*ff-ee*ee
  endif
xdummy = 4*dd*ff-ee*ee
ydummy = dd0 + dd1*xx + dd2*xx2
zdummy = ydummy + xterm0
xxdummy = (atan((ee+darg)/sqrt(xdummy)) - atan(ee/sqrt(xdummy)))/sqrt(xdummy)
yydummy = aa0 + aa1*xx + aa2*xx2
  if (lowdarg) then
    dargoe = darg/ee
    test1 = (dargoe+dargoe**2+(dargoe**3)/3.d0)/((ee+darg)**3)
    test2 = delta*(dargoe+2*dargoe**2+2*dargoe**3+dargoe**4+(dargoe**5)/5.d0)/((ee+darg)**5)
  else
    test1 = (1.d0/ee**3 - 1.d0/(ee+darg)**3)/3.d0
    test2 = delta*(1.d0/ee**5 - 1.d0/(ee+darg)**5)/5.d0
  endif
  if (abs(test2/test1).lt.1.d-4) then
! need series expansion to avoid numerical instability due to differences in large numbers
! dgammaD = (gamma - (1/ee - 1/(ee+darg)))/delta = (gamma - (1/ee)(darg/(ee+darg)))/delta
idummy=0
    dgammaD = test2-test1
    ii = 2
    xi = -1
    do
      jj = 2*ii+3 
      if (lowdarg) then
        ll = 1  ! binomial coefficient
        sum1 = 0.d0
        prod1 = 1.d0
        do kk=1,jj  ! binomial series expansion for ((ee+darg)**jj - ee**jj)/ee**jj
          ll = (ll*(jj-kk+1))/(kk)
          prod1 = prod1*dargoe
          sum1 = sum1 + ll*prod1
        enddo
        test3 = delta**ii * (sum1/(dble(jj)*(ee+darg)**jj))/dble(jj)
      else
        test3 = delta**ii * (1.d0/ee**jj - 1.d0/(ee+darg)**jj)/dble(jj)
      endif
      dgammaD_old = dgammaD
      dgammaD = dgammaD + xi*test3
      if (dgammaD.eq.dgammaD_old) exit
      xi = -1*xi
      ii = ii+1
    enddo
    temp = zz/(ee*dterm1*xterm1)
    tt0 =  (aa/dd)*((2*dd*ff-dterm1*(ff*zz+ee))*temp + 4*dd*ff*dgammaD)
    tt1 = -bb*(ee*temp + 2*ee*dgammaD)
    tt2 = cc*((2*dd + dterm1*zz)*temp + 4*dd*ff*dgammaD)
    tfn1 = tt0 + tt1 + tt2
  else 
    if (delta.lt.0.d0) then
idummy=1
      write(6,*) "ERROR in function tfn1: negative delta is unphysical"
      stop
    else
idummy=2
      sqdelta = sqrt(delta)
      dtarg = darg/sqdelta
      arg1 = ee/sqdelta
      arg2 = arg1 + dtarg
      testlo1 = zz/(2*dd) ! = darg/(delta+ee*ee)
      testlo2 = -zz*zz*ee/(4*dd*dd) ! = darg**2*ee/(delta+ee*ee)**2
      testlo3 = zz**3 * (ee*ee - delta/3.d0) / (8*dd**3)
      test3 = darg/(ee*(ee+darg))
      test4 = ((ee*ee*darg+ee*darg*darg+darg**3/3.d0)/(ee*(ee+darg))**3)*delta
      ta0 = 2*dd*ff-ee*ee-ee*ff*zz
      ta1 = ee+2*ff*zz
      ta2 = ee*zz+2*dd
      if ((abs(testlo1).gt.abs(1.d3*testlo2)).and.(abs(testlo1).gt.abs(1.d6*testlo3))) then
! if difference in atan arguments small, may lead to numerical error
! atan(z+d)->atan(z)+d/(1+z^2)-d^2z/(1+z^2)^2+d^3(z^2-1/3)/(1+z^2)^3+..., d->0 (taylor expansion)
        call atan_expansion(darg/sqdelta,arg1,hiorderexp,3)
        hiorderexp = hiorderexp/sqdelta
        gammat = testlo1+testlo2+testlo3 + hiorderexp
jdummy=1
        tt0 = aa2dd*(zz*ta0/xterm1 + 4*ff*dd*gammat)/delta
        tt1 = bb*zz*zz*(dd*delta-ee*zz*(2*dd*ff-ee*ee-ee*ff*zz)) &
&            /(2*dd*dd*delta*xterm1) &
&           - 2*ee*bb*hiorderexp/delta
        tt2 = cc*(zz**3*(2*ff*dd - ff*zz*ee - ee*ee))/(dd*delta*xterm1)   &
&           + 4*dd*hiorderexp*cc/delta
        tfn1 = tt0 + tt1 + tt2
!        atf = gammat*sqdelta
      elseif (abs(test3).gt.abs(1.d6*test4)) then
! if both atan arguments large, may lead to numerical error
! atan(z)->pi/2-1/z+1/(3z^3)-1/(5z^5)+..., z->infty
! atan(z)->-pi/2-1/z+1/(3z^3)-1/(5z^5)+..., z->-infty
        gammat = test3-test4
!        atf = gammat*sqdelta
        tt0 = aa2dd*(zz*ta0/xterm1 + 4*ff*dd*gammat)/delta
        tt1 = -bb*(-zz*ta1/xterm1 + 2*ee*gammat)/delta
        tt2 = cc*(-zz*ta2/xterm1 + 4*dd*gammat)/delta
        tfn1 = tt0 + tt1 + tt2
jdummy=2
      else
        atf = atan(arg2) - atan(arg1)
        gammat = atf/sqdelta
        if (gammapfnroot.gt.0) then
          quadterm3 = gammapfpf*(xx-gammapfroot1)*(xx-gammapfroot2)
        else
          quadterm3 = gammapfpf*((xx-gammapfroot1)**2 + gammapfroot2**2)
        endif
        if (denomnroot.ge.4) then
          quadterm1 = (xx-denomroot3)*(xx-denomroot4)
        else
          quadterm1 = (xx-denomroot3)**2+denomroot4**2
        endif
        if (denomnroot.ge.2) then
          quadterm2 = ((xx-denomroot1)*(xx-denomroot2))*denompf
        else
          quadterm2 = ((xx-denomroot1)**2+denomroot2**2)*denompf
        endif
        tc1 = zz*quadterm1*quadterm2/(dd*xterm1)/delta
        tc2 = quadterm3*gammat/delta
        test1 = tc1+tc2
        tt0 = aa2dd*(zz*ta0/xterm1 + 4*ff*dd*gammat)/delta
        tt1 = -bb*(-zz*ta1/xterm1 + 2*ee*gammat)/delta
        tt2 = cc*(-zz*ta2/xterm1 + 4*dd*gammat)/delta
        test2 = tt0 + tt1 + tt2
! test1 and test2 are mathematically equivalent, choose the representation that minimizes error due to difference of similar large numbers
        if (max(abs(tc1),abs(tc2)).lt.max(abs(tt0),abs(tt1),abs(tt2))) then
          tfn1 = test1
        else
          tfn1 = test2
        endif
jdummy=0
      endif
    endif
  endif

!if (i_prt_DEBUG.gt.0) then
!write(79+l_ct_DEBUG,*) xx,gammat,atf,atan(arg2),atan(arg1),sqdelta
!endif
if (i_ct_DEBUG.eq.i_ref_DEBUG) then
j_ct_DEBUG = j_ct_DEBUG+1
!write(77,*) xx,tfn1,idummy,jdummy,i_ct_DEBUG
!write(77,*) xx,tfn1,tt0,tt1,tt2
!write(77,*) xx, &
!& (aa2dd/delta)*(zz*ta0/xterm1 + 4*ff*dd*gammat), &
!& (-bb/delta)*(-zz*ta1/xterm1 + 2*ee*gammat), &
!& (cc/delta)*(-zz*ta2/xterm1 + 4*dd*gammat), &
!& (yydummy/(ydummy*xdummy))*(zz*(2*ydummy*ff-ee*ee-ee*ff*zz)/zdummy+4*ydummy*ff*xxdummy), &
!& -(bb/delta)*(-zz*(ee+2*ff*zz)/zdummy+2*ee*xxdummy), &
!& (cc/delta)*(-zz*(ee*zz+2*ydummy)/zdummy+4*ydummy*xxdummy)
!write(77,*) xx,tfn1,idummy,jdummy,(-zz*ta2/xterm1 + 4*dd*gammat),(zz**3*(2*ff*dd - ff*zz*ee - ee*ee))/(dd*xterm1)
!write(77,*) xx,tfn1,tt2,-zz*ta2/xterm1,4*dd*gammat
!write(77,*) xx,tfn1,delta,4*ff*(sqrt(dd0)*omxi22 + xx*sqrt(dd2)*omzeta122)**2
!if (jdummy.eq.1.and.idummy.eq.2) write(77,*) xx,tfn1,test1,test2,darg,ee,delta
!
!if (j_ct_DEBUG+1.gt.1000) return
if (j_ct_DEBUG+1.gt.1000) then
write(6,*) "Programmed DEBUG stop"
if (idummy.eq.0) then
  write(6,*) "Gamma from series expansion in Delta, idummy = ", idummy
else if (idummy.eq.1) then
  write(6,*) "Gamma from natural logarithm function, idummy = ", idummy
else
  write(6,*) "Gamma from arctangent function, idummy = ", idummy
endif
stop
endif
endif

if (.not.(tfn1.eq.tfn1)) then
write(6,*) "NAN error in tfn1"
write(6,*) "xx = ",xx
write(6,*) "zz = ",zz
if (idummy.eq.0) then
  write(6,*) "Gamma from series expansion in Delta"
else if (idummy.eq.1) then
  write(6,*) "Gamma from natural logarithm function"
else
  write(6,*) "Gamma from arctangent function"
endif
write(6,*) idummy,jdummy,kdummy
write(6,*) xdummy,ydummy,zdummy
write(6,*) xxdummy,yydummy,zzdummy
write(6,*) "ta0 = ",ta0
write(6,*) "ta1 = ",ta1
write(6,*) "ta2 = ",ta2
write(6,*) "tt0 = ",tt0
write(6,*) "tt1 = ",tt1
write(6,*) "tt2 = ",tt2
write(6,*) "aa = ",aa
write(6,*) "bb = ",bb
write(6,*) "cc = ",cc
write(6,*) "dd = ",dd
write(6,*) "ee = ",ee
write(6,*) "ff = ",ff
write(6,*) "dd0 = ",dd0
write(6,*) "dd1 = ",dd1
write(6,*) "dd2 = ",dd2
write(6,*) "ee0 = ",ee0
write(6,*) "ee1 = ",ee1
write(6,*) "rrootd = ",rrootd
write(6,*) "irootd = ",irootd
write(6,*) "delta = ",delta
write(6,*) "deltaprefact = ",deltaprefact
write(6,*) "q02k1 = ",q02k1
write(6,*) "kappa = ",kappa
write(6,*) "lambda = ",lambda
write(6,*) "darg = ",darg
write(6,*) "xterm1 = ",xterm1
write(6,*) "gammat = ",gammat
write(6,*) "atf = ",atf
write(6,*) "sqdelta = ",sqdelta
stop
endif

end function tfn1

!***************************************************************************
!
! expand (1+h1*x)(1+h2*x)(1+h3*x)(1+h4*x) = 1 + x*term1 + x^2*term2 + x^3*term3 + x^4*term4

subroutine expand4(h1,h2,h3,h4,term1,term2,term3,term4)
double precision :: h1,h2,h3,h4,term1,term2,term3,term4
term1 = h1+h2+h3+h4
term2 = h1*h2 + h1*h3 + h1*h4 + h2*h3 + h2*h4 + h3*h4
term3 = h2*h3*h4 + h1*h3*h4 + h1*h2*h4 + h1*h2*h3
term4 = h1*h2*h3*h4
return
end subroutine expand4

!***************************************************************************
!
! sine of angle between two vectors expressed in reduced coordinates

subroutine sinreducedangle(v1,v2,ss)
use geometry
double precision :: v1(3),v2(3),ss
double precision :: w1(3),w2(3),w3(3),w1mag,w2mag,w3mag
integer :: ii
do ii=1,3
  w1(ii) = dot_product(v1,blat(:,ii))
  w2(ii) = dot_product(v2,blat(:,ii))
enddo
call cross(w1,w2,w3)
w1mag = sqrt(dot_product(w1,w1))
w2mag = sqrt(dot_product(w2,w2))
w3mag = sqrt(dot_product(w3,w3))
ss = w3mag/(w1mag*w2mag)
return
end subroutine sinreducedangle

!***************************************************************************
!
! cross product (returned in cartesian coordinates) of two vectors expressed in reduced coordinates

subroutine crossreduced(v1,v2,w3)
use geometry
double precision :: v1(3),v2(3),ss
double precision :: w1(3),w2(3),w3(3)
integer :: ii
do ii=1,3
  w1(ii) = dot_product(v1,blat(:,ii))
  w2(ii) = dot_product(v2,blat(:,ii))
enddo
call cross(w1,w2,w3)
return
end subroutine crossreduced

!***************************************************************************
! atan(zz+dz) - atan(zz)
! using taylor expansion of atan(zz+dz) in powers of dz
! startterm > 1 ignores the first startterm-1 terms in the series
! only recommended for dz << zz

subroutine atan_expansion(dz,zz,atex,startterm)
  implicit none
  double precision :: dz,zz,atex
  integer :: startterm
  integer :: nterm
  parameter (nterm = 50)
  double precision :: termprefact(0:nterm)
  integer :: ii,jj,kk,ll,mm
  integer :: iparity,jparity
  double precision :: sum1, sum1_old, sum2, term
  double precision :: zzpow, dzpow, denompow, denomterm

! first term is 1/(1+zz*zz)
! build up subsequent terms by induction
! (d/d(zz)) zz^(ii-mm)/(1+zz*zz)^ii = (-(ii+mm)*zz^(ii-mm+1) + (ii-mm)*zz^(ii-mm-1))/(1+zz*zz)^(ii+1)
! allows power series in zz of numerator.  
! Coefficients of zz^jj stored in termprefact(jj)
! Because any term has only odd or even powers of zz, can interlace coefficients from current and next terms in one array (and can forget previous terms)
! Might as well also include factorial part of Taylor expansion in termprefact (division by ii+1)
  termprefact(0) = 1.d0;
  sum1 = 0.d0
  dzpow = 1.d0
  denomterm = 1 + zz*zz
  denompow = 1.d0
  do ii=1,nterm
    iparity = mod(ii-1,2)+1
    jparity = mod(ii,2)+1
!write(6,*) ii,iparity,jparity
    dzpow = dzpow * dz
    denompow = denompow * denomterm
    if (iparity.eq.1) then
      zzpow = 1.d0
    else
      zzpow = zz
    endif
    sum2 = 0.d0
    do jj = jparity-1,ii,2
      termprefact(jj) = 0.d0
    enddo
    do jj = iparity-1,ii-1,2
      mm = ii-jj
      termprefact(jj+1) = termprefact(jj+1) + termprefact(jj)*(-(ii+mm))/dble(ii+1)
      if (jj.gt.0) then
        termprefact(jj-1) = termprefact(jj-1) + termprefact(jj)*(ii-mm)/dble(ii+1)
      endif
!write(6,'(1x,3i3,4x,2i3,f10.3,i3,f10.3,f10.3)') ii,jj,mm,ii+1,jj+1, &
!&  termprefact(jj),-(ii+mm),(-(ii+mm))/dble(ii+1), &
!&  termprefact(jj+1)
!if (jj.gt.0) then
!write(6,'(14x,2i3,f10.3,i3,f10.3,f10.3)')                ii+1,jj-1, &
!&  termprefact(jj),(ii-mm),(ii-mm)/dble(ii+1), &
!&  termprefact(jj-1)
!else 
!write(6,*)
!endif
      sum2 = sum2 + termprefact(jj)*zzpow
      zzpow = zzpow * zz * zz
    enddo
    if (ii.ge.startterm) then
      term = dzpow * sum2 / (denompow)
      sum1_old = sum1
      sum1 = sum1 + term
!write(6,*) "**** ",ii,"      ",term,sum1
!read(*,*)
      if (sum1_old.eq.sum1) exit
    endif
  enddo
  atex = sum1

end subroutine atan_expansion










