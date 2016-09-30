subroutine mkploss(iband,ikk,lossfn, &
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
& nband,nsppol,shiftk,total_se,ploss)
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
double complex :: sei,seipw,seiold,seibc
double complex, allocatable :: ssi(:)
double precision :: ploss
double precision :: ploss1(-nwpt:nwpt)
double precision :: total_se
double complex :: ctest
integer :: ii,jj,kk,ll,ix,iy,iz,iskip,ocsign,iiq
integer :: iqq(3),jka(3),jkb(3),jkk(3),iqv(3),ictr(3),iqr(3),ixr(3)
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
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3)), omega_t
double precision :: whi,wlo,wwhi,wwlo,ww,www,wwwlo,wwwhi,enval(8),dw,rw,eshift,evtx0(4),evtx(4)
double precision :: rrpyr(3,4),kvtx(3,4),xk(3,4),ckvtx(3,4),cxk(3,4),rg(3,3),tvol,xpyr0(4),xpyr(4)
double precision :: xkc(3,4),rgc(3,3),tvolc,kvtxc(3,4)
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
logical :: lc
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

!write(6,*) "qopf ",total_se,ploss
i_prt_DEBUG = 0
j_prt_DEBUG = 0
k_prt_DEBUG = 0
i_ct_DEBUG = 0
j_ct_DEBUG = 0
k_ct_DEBUG = 0
l_ct_DEBUG = 0
i_ref_DEBUG = -1
k_ref_DEBUG = -1
abr=1.d-4
rlr=1.d-4
sei=0.d0
ctest=(0.d0,0.d0)
igg0=(/0,0,0/)
ploss=0.d0
do ie=-nwpt,nwpt
  ploss1(ie)=0.d0
enddo
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
!do ibp=2,2
!do ibp=iband,iband
  write(6,'(3(a,i4))') 'quasiparticle band ',iband,', polarized band ',ibp,' out of ',test_bands_se(2)
  seibc=sei
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
    ipw1=ipwndx(1,iipw)
    ipw2=ipwndx(2,iipw)
!    if (ipw1.ne.ipw2) cycle
!    write(6,'(1x,i3,4x,2i4)') ibp,ipw1,ipw2 
    lc=iipw.le.ntpwndx  ! do correlation
    if (.not.lc) cycle

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
      omega_t=enrgy(indxkbnd(ikptq)+ibp)-ek
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
          smat2=4.d0*pi*omega_t
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
            jmat2(ii,jj)=4.d0*pi*jmatel(ii)*conjg(jmatel2(jj))*omega_t
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
        xmat2(ix,iy,iz)=(cmatel+amatel)*conjg(cmatel2+amatel2)*omega_t
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
        call fval(omega,ixx,ixxp,ngkpt,enval)
        if (lqsing) call fval(vq,ixx,ixxp,ngkpt,vq2)
        if (lc) then
          do iw=1,nwpt
            call fpol(vfactor(:,:,:,iw),ngkpt,ixx,ixxp,vfv(:,iw))
          enddo
        endif
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
              if (ixr(1).eq.ictr(1).and.ixr(2).eq.ictr(2).and.ixr(3).eq.ictr(3)) then
                if (lc) then
                  do iw=1,nwpt
                    vpyr0(iv,iw) = 0.d0
                  enddo
                endif
              else
                qq = dble(iqr)/dble(ngkpt)
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
              endif
            enddo
            bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
            do jj=1,3
              bgrad=bgrad+(evtx(jj+1)-evtx(1))*rgc(:,jj)
            enddo
            xmult=1.d0/sqrt(dot_product(bgrad,bgrad))
! DEBUG
!cdum = 0.d0
! DEBUG
!write(6,*) "first set"
            do iww=iwwlim(1),iwwlim(2) ! energies between evtx(1) & evtx(2)
! DEBUG
!write(6,*) iww
              wwwlo=max(dw*(iww-0.5d0),evtx(1))
              wwwhi=min(dw*(iww+0.5d0),evtx(2))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              do ii=1,3
              do jj=1,3
! DEBUG
k_ct_DEBUG = k_ct_DEBUG + 1
!write(6,*) "k_ct_DEBUG: ",k_ct_DEBUG,ii,jj
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E_LIM:      ",wwwlo,wwwhi
!if (k_ct_DEBUG.le.9) write(77,*) "# ",ii,jj
                call tsinggrater(ii,jj,wwwlo,wwwhi,1,xk,kvtx(:,1),ngkpt,evtx,abr,rlr,nsing,wsing,stens(ii,jj))
              enddo
              enddo
              stens = stens * xmult
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
!! DEBUG
!cdum = 0.d0
!do ii=1,3
!do jj=1,3
!ctensdum2(ii,jj) = (0.d0,0.d0)
!do kk=1,3
!do ll=1,3
!ctensdum2(ii,jj) = ctensdum2(ii,jj) + blat(kk,ii)*vsfactor(kk,ll,iw)*blat(ll,jj)
!enddo
!enddo
!enddo
!enddo
                  do ii=1,3
                  do jj=1,3
                    ssi(ie)=ssi(ie)+vsm2(ii,jj,iw)*stens(jj,ii)
!! DEBUG
!if (ie.eq.0.or.ie.eq.1) then
!cdum = cdum+vsm2(ii,jj,iw)*stens(jj,ii)
!endif
!if (iw.eq.158) then
!write(6,*) ii,jj,dimag(ctensdum2(ii,jj)),dimag(vsfactor(ii,jj,iw)),dimag(vsm2(ii,jj,iw)),stens(jj,ii)
!endif
                  enddo
                  enddo
! DEBUG
!if (iw.eq.158) then
!write(6,'(i4,23e14.3)') iw,cdum
!stop
!endif
!if (cdum.ne.(0.d0,0.d0)) then
!write(6,'(i4,23e14.3)') iw,cdum
!read(*,*)
!endif
                enddo
              endif
! DEBUG
!write(6,'(5i3,2x,2f16.10)') ix,iy,iz,itet,iww,cdum
            enddo
! DEBUG
!write(6,*) 'Khst7',xdum
!read(*,*)
!write(6,*) "second set"
            do iww=iwwlim(2),iwwlim(3) ! energies between evtx(2) & evtx(3)
! DEBUG
!write(6,*) iww
              wwwlo=max(dw*(iww-0.5d0),evtx(2))
              wwwhi=min(dw*(iww+0.5d0),evtx(3))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              do ii=1,3
              do jj=1,3
! DEBUG
k_ct_DEBUG = k_ct_DEBUG + 1
!write(6,*) "k_ct_DEBUG: ",k_ct_DEBUG,ii,jj
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E_LIM:      ",wwwlo,wwwhi
                call tsinggrater(ii,jj,wwwlo,wwwhi,2,xk,kvtx(:,1),ngkpt,evtx,abr,rlr,nsing,wsing,sint1a)
! DEBUG
k_ct_DEBUG = k_ct_DEBUG + 1
!write(6,*) "k_ct_DEBUG: ",k_ct_DEBUG,ii,jj
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E_LIM:      ",wwwlo,wwwhi
                call tsinggrater(ii,jj,wwwlo,wwwhi,3,xk,kvtx(:,1),ngkpt,evtx,abr,rlr,nsing,wsing,sint1b)
                stens(ii,jj) = sint1a+sint1b
              enddo
              enddo
              stens = stens * xmult
! DEBUG
!cdum = 0.d0
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  do ii=1,3
                  do jj=1,3
                      ssi(ie)=ssi(ie)+vsm2(ii,jj,iw)*stens(jj,ii)
! DEBUG
!if (ie.eq.0.or.ie.eq.1) then
!cdum = cdum+vsm2(ii,jj,iw)*stens(jj,ii)
!endif
                  enddo
                  enddo
                enddo
              endif
! DEBUG
!write(6,'(5i3,2x,2f16.10)') ix,iy,iz,itet,iww,cdum
            enddo
! DEBUG
!write(6,*) 'Khst7',xdum
!write(6,*) "third set"
            do iww=iwwlim(3),iwwlim(4) ! energies between evtx(3) & evtx(4)
! DEBUG
!write(6,*) iww
              wwwlo=max(dw*(iww-0.5d0),evtx(3))
              wwwhi=min(dw*(iww+0.5d0),evtx(4))
              if (abs(wwwhi-wwwlo)*1.d14.lt.max(abs(wwwhi),abs(wwwlo),dw)) cycle
              do ii=1,3
              do jj=1,3
! DEBUG
k_ct_DEBUG = k_ct_DEBUG + 1
!write(6,*) "k_ct_DEBUG: ",k_ct_DEBUG,ii,jj
!if (k_ct_DEBUG.eq.k_ref_DEBUG) write(6,*) "E_LIM:      ",wwwlo,wwwhi
                call tsinggrater(ii,jj,wwwlo,wwwhi,4,xk,kvtx(:,1),ngkpt,evtx,abr,rlr,nsing,wsing,stens(ii,jj))
              enddo
              enddo
              stens = stens * xmult
! DEBUG
!cdum = 0.d0
              if (lc) then
                do iw=1,nwpt
                  ie=iww-ocsign*iw
                  do ii=1,3
                  do jj=1,3
                      ssi(ie)=ssi(ie)+vsm2(ii,jj,iw)*stens(jj,ii)
! DEBUG
!if (ie.eq.0.or.ie.eq.1) then
!cdum = cdum+vsm2(ii,jj,iw)*stens(jj,ii)
!endif
                  enddo
                  enddo
                enddo
              endif
! DEBUG
!write(6,'(5i3,2x,2f16.10)') ix,iy,iz,itet,iww,cdum
            enddo
! DEBUG
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
                  ssi(ie)=ssi(ie)+cint/volelmnt
! DEBUG
!if (ie.eq.0.or.ie.eq.1) cdum = cdum+cint/volelmnt
                enddo
              enddo
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
                ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
! DEBUG
!if (ie.eq.0.or.ie.eq.1) cdum = cdum+(aa0(iw)*sint1+dot_product(svec,av(:,iw)))/volelmnt
              enddo
            enddo
          endif
! DEBUG
!write(32,'(4i3,2x,f14.10,2x,f14.10)') ixx,itet,(cdum/2.d0)
        enddo
      enddo
      enddo
      enddo
!    if (lc) ssi=-ssi/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*pi)
      if (lc) ssi=-ocsign*ssi/((2*pi)**3*pi)
  
! DEBUG
!write(6,*) "Integration complete"

      if (lc) then
        do ie=max(-nwpt,ielo),min(nwpt,iehi)
          ploss1(ie)=ploss1(ie)+abs(dimag(ssi(ie)))*pi
        enddo
      endif
      deallocate (ssi)
    endif

  enddo
enddo

do ie=-nwpt,nwpt
  eps1=dble(ie)*dw-dw/2
!  write(14,'(i4,f10.3,2(3x,2es12.3))') ie,eps1*27.2114,ploss1(ie)
enddo
write(14,*)
if (itetrahedron.eq.1) then
  iw = int(total_se/dw)
  rw = total_se - iw*dw
!  write(14,*) iw,rw
  ploss=ploss1(iw)*((dw-rw)/dw)+ploss1(iw+1)*(rw/dw)
!  ploss=(ploss1(0)+ploss1(1))/2.
!  write(14,*) ploss
endif

!write(6,*) "kwls ",total_se,ploss

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
end subroutine mkploss

