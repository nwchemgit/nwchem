subroutine fpwmatel(qq,igg,pwmatel,tpwmatel)
use pawvars
use geometry
implicit none
double precision :: qq(3)
integer :: igg(3)
double complex :: pwmatel(ntypepaw,nlmnmax,nlmnmax),tpwmatel(ntypepaw,nlmnmax,nlmnmax)
integer :: iorb1,iorb2,itpaw
double precision :: qqgg(3),qqggmag,rr,dr
double precision :: phival1,phival2,tphival1,tphival2
double precision :: radintv,radintt,tripleYlmint,bessj,sphbesselj,pi
parameter (pi=3.1415926535897932384626433832795028841972d0)
double complex :: prefact,spharm,ylm
integer :: ll1,ll2,lltd,mm1,mm2,mmt,nn1,nn2,ir,ii,jj
external spharm,tripleYlmint,sphbesselj

qqgg=(/0.d0,0.d0,0.d0/)
do ii=1,3
  do jj=1,3
    qqgg(jj)=qqgg(jj)+blat(ii,jj)*(-qq(ii)-dble(igg(ii)))
  enddo
enddo
qqggmag=sqrt(dot_product(qqgg,qqgg))
!write(6,*) qqgg
!write(6,*) qqggmag
!write(6,*) rphigrid(itpaw,gridsize(itpaw,iphigrid(itpaw)))
!write(6,*) qqggmag*rphigrid(itpaw,gridsize(itpaw,iphigrid(itpaw)))

do itpaw=1,ntypepaw
  do iorb1=1,nlmn(itpaw)
!do iorb1=7,7
  do iorb2=1,nlmn(itpaw)
!do iorb2=7,7

    nn1=iorbno(itpaw,iorb1)
    nn2=iorbno(itpaw,iorb2)
    ll1=llorb(itpaw,iorb1)
    mm1=mmorb(itpaw,iorb1)
    ll2=llorb(itpaw,iorb2)
    mm2=mmorb(itpaw,iorb2)
!write(6,*) iorb1,iorb2,ll1,ll2
!write(6,*) ll1,mm1
!write(6,*) ll2,mm2

    pwmatel(itpaw,iorb1,iorb2)=0.d0
    tpwmatel(itpaw,iorb1,iorb2)=0.d0
    do lltd=abs(ll1-ll2),ll1+ll2
!write(6,*) lltd
      mmt=-mm1+mm2
      if (abs(mmt).gt.lltd) cycle
      if (qqggmag.lt.1.d-8) then
        if (lltd.ne.0) cycle
        ylm=1.d0/sqrt(4.d0*pi)
      else
        ylm=spharm(lltd,mmt,qqgg)
      endif
      prefact=4.d0*pi*tripleYlmint(ll1,ll2,lltd,-mm1,mm2)   &
&            *ylm*(-1)**(mm2)*(0.d0,1.d0)**lltd
      radintv=0.d0
      radintt=0.d0
!write(18,'(a,3i3)') 'iorb',iorb1,iorb2
!write(18,'(a,3i3)') 'nn',nn1,nn2
!write(18,'(a,3i3)') 'll',ll1,ll2,lltd
!write(18,'(a,3i3)') 'mm',mm1,mm2,mmt
!write(18,*) tripleYlmint(ll1,ll2,lltd,-mm1,mm2)
      do ir=1,gridsize(itpaw,iphigrid(itpaw))-1
        rr=(rphigrid(itpaw,ir+1)+rphigrid(itpaw,ir))/2.d0
        dr=rphigrid(itpaw,ir+1)-rphigrid(itpaw,ir)
        bessj=sphbesselj(qqggmag*rr,lltd)
        phival1=(phi(itpaw,nn1,ir+1)+phi(itpaw,nn1,ir))/2.d0
        phival2=(phi(itpaw,nn2,ir+1)+phi(itpaw,nn2,ir))/2.d0
        tphival1=(tphi(itpaw,nn1,ir+1)+tphi(itpaw,nn1,ir))/2.d0
        tphival2=(tphi(itpaw,nn2,ir+1)+tphi(itpaw,nn2,ir))/2.d0
        radintv=radintv+dr*phival1*phival2*bessj
        radintt=radintt+dr*tphival1*tphival2*bessj
!write(18,'(i5,5es14.6)') ir,rr,bessj,phival1,phival2,radintv
      enddo
!write(6,'(3i3)') iorb1,iorb2
!write(6,'(3i3)') nn1,nn2
!write(6,'(3i3)') ll1,ll2,lltd
!write(6,'(3i3)') mm1,mm2,mmt
!write(6,*) tripleYlmint(ll1,ll2,lltd,-mm1,mm2)
!write(6,*) radintv
      pwmatel(itpaw,iorb1,iorb2)=pwmatel(itpaw,iorb1,iorb2)+prefact*radintv
      tpwmatel(itpaw,iorb1,iorb2)=tpwmatel(itpaw,iorb1,iorb2)+prefact*radintt
    enddo
  enddo
  enddo
enddo

!do iorb1=1,nlmnmax
!  write(6,'(14f6.2)') dble(pwmatel(iorb1,:))
!enddo
!write(6,*)
!do iorb1=1,nlmnmax
!  write(6,'(14f6.2)') dble(tpwmatel(iorb1,:))
!enddo

return
end subroutine fpwmatel

!******************************************************************************

subroutine fpwjmatel(pwjmatel,tpwjmatel)
use pawvars
use geometry
implicit none
double complex :: pwjmatel(3,ntypepaw,nlmnmax,nlmnmax),tpwjmatel(3,ntypepaw,nlmnmax,nlmnmax)
integer :: iorb1,iorb2,itpaw
double precision :: rr,dr
double precision :: phival1,phival2,tphival1,tphival2
double precision :: dphival1,dphival2,dtphival1,dtphival2
double precision :: radintv1,radintt1,radintv2,radintt2
double precision :: tripleYlmint,pi
double complex :: coni
parameter (coni=(0.d0,1.d0))
parameter (pi=3.1415926535897932384626433832795028841972d0)
double complex :: prefact,spharm,ylm,prevec1(3),prevec2(3),vectermv(3),vectermt(3)
integer :: ll1,ll2,lltd,mm1,mm2,mmt,nn1,nn2,ir,ii,jj
external spharm,tripleYlmint

do itpaw=1,ntypepaw
  do iorb1=1,nlmn(itpaw)
!do iorb1=7,7
  do iorb2=1,nlmn(itpaw)
!do iorb2=7,7

    nn1=iorbno(itpaw,iorb1)
    nn2=iorbno(itpaw,iorb2)
    ll1=llorb(itpaw,iorb1)
    mm1=mmorb(itpaw,iorb1)
    ll2=llorb(itpaw,iorb2)
    mm2=mmorb(itpaw,iorb2)
!write(6,*) iorb1,iorb2,ll1,ll2
!write(6,*) ll1,mm1
!write(6,*) ll2,mm2

    do ii=1,3
      pwjmatel(ii,itpaw,iorb1,iorb2)=0.d0
      tpwjmatel(ii,itpaw,iorb1,iorb2)=0.d0
    enddo
    if (abs(ll1-ll2).le.1.and.(ll1+ll2).ge.1) then
      mmt=-mm1+mm2
      if (abs(mmt).gt.1) cycle
      if (mmt.eq.-1) then
        prefact = -(-1)**(mm1)*sqrt(2.d0*pi/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2)
        prevec1(1) = prefact/2.d0
        prevec1(2) = -prefact*coni/2.d0
        prevec1(3) = 0.d0
      else if (mmt.eq.1) then
        prefact = (-1)**(mm1)*sqrt(2.d0*pi/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2)
        prevec1(1) = prefact/2.d0
        prevec1(2) = prefact*coni/2.d0
        prevec1(3) = 0.d0
      else if (mmt.eq.0) then
        prevec1(1) = 0.d0
        prevec1(2) = 0.d0
        prevec1(3) = (-1)**(mm1)*sqrt(pi/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2)
      else
        prevec1(1) = 0.d0
        prevec1(2) = 0.d0
        prevec1(3) = 0.d0
      endif
      if (mmt.eq.-1) then
!        if (ll2.eq.ll1-1) then
!          prefact = -ll1*sqrt((ll1-mm1-1)*(ll1-mm1)/dble((2*ll1-1)*(2*ll1+1)))
!        else if (ll2.eq.ll1+1) then
!          prefact = -(ll1+1)*sqrt((ll1+mm1+1)*(ll1+mm1+2)/dble((2*ll1+1)*(2*ll1+3)))
!        else
!          prefact = (0.d0,0.d0)
!        endif
        prefact = ((-1)**(mm1)/2.d0)*( &
&   sqrt(8*pi/3.d0)*mm2*tripleYlmint(ll1,ll2,1,-mm1,mm2) &
& + sqrt(4*pi*(ll2-mm2)*(ll2+mm2+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2+1) &
& + sqrt(8*pi/3.d0)*mm1*tripleYlmint(ll1,ll2,1,-mm1,mm2) &
& - sqrt(4*pi*(ll1+mm1)*(ll1-mm1+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1+1,mm2) )
        prevec2(1) = prefact/2.d0
        prevec2(2) = -prefact*coni/2.d0
        prevec2(3) = 0.d0
      else if (mmt.eq.1) then
!        if (ll2.eq.ll1-1) then
!          prefact = ll1*sqrt((ll1+mm1-1)*(ll1+mm1)/dble((2*ll1-1)*(2*ll1+1)))
!        else if (ll2.eq.ll1+1) then
!          prefact = (ll1+1)*sqrt((ll1-mm1+1)*(ll1-mm1+2)/dble((2*ll1+1)*(2*ll1+3)))
!        else
!          prefact = (0.d0,0.d0)
!        endif
        prefact = ((-1)**(mm1)/2.d0)*( &
&   sqrt(8*pi/3.d0)*mm2*tripleYlmint(ll1,ll2,1,-mm1,mm2) &
& - sqrt(4*pi*(ll2+mm2)*(ll2-mm2+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2-1) &
& + sqrt(8*pi/3.d0)*mm1*tripleYlmint(ll1,ll2,1,-mm1,mm2) &
& + sqrt(4*pi*(ll1-mm1)*(ll1+mm1+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1-1,mm2) )
        prevec2(1) = prefact/2.d0
        prevec2(2) = prefact*coni/2.d0
        prevec2(3) = 0.d0
      else if (mmt.eq.0) then
!        if (ll2.eq.ll1-1) then
!          prefact = ll1*sqrt((ll1+mm1)*(ll1-mm1)/dble((2*ll1-1)*(2*ll1+1)))
!        else if (ll2.eq.ll1+1) then
!          prefact = (ll1+1)*sqrt((ll1+mm1+1)*(ll1-mm1+1)/dble((2*ll1+1)*(2*ll1+3)))
!        else
!          prefact = 0.d0
!        endif
        prefact = -((-1)**(mm1)/2.d0)*( &
&   sqrt(2*pi*(ll2-mm2)*(ll2+mm2+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2+1) &
& + sqrt(2*pi*(ll2+mm2)*(ll2-mm2+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1,mm2-1) &
& - sqrt(2*pi*(ll1-mm1)*(ll1+mm1+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1-1,mm2) &
& - sqrt(2*pi*(ll1+mm1)*(ll1-mm1+1)/3.d0)*tripleYlmint(ll1,ll2,1,-mm1+1,mm2) )
        prevec2(1) = 0.d0
        prevec2(2) = 0.d0
        prevec2(3) = prefact
      else
        prevec2(1) = 0.d0
        prevec2(2) = 0.d0
        prevec2(3) = 0.d0
      endif
!write(18,'(i2,2x,3i2,2x,3i2)') itpaw,nn1,ll1,mm1,nn2,ll2,mm2
      radintv1=0.d0
      radintv2=0.d0
      radintt1=0.d0
      radintt2=0.d0
!write(18,'(a,3i3)') 'iorb',iorb1,iorb2
!write(18,'(a,3i3)') 'nn',nn1,nn2
!write(18,'(a,3i3)') 'll',ll1,ll2,lltd
!write(18,'(a,3i3)') 'mm',mm1,mm2,mmt
!write(18,*) tripleYlmint(ll1,ll2,lltd,-mm1,mm2)
!write(18,'(i2,2x,3i2,2x,3i2,3x,4(f7.3,2x)))') &
!& itpaw,nn1,ll1,mm1,nn2,ll2,mm2, &
!& tripleYlmint(ll1,ll2,1,-mm1,mm2+1), &
!& tripleYlmint(ll1,ll2,1,-mm1,mm2-1), &
!& tripleYlmint(ll1,ll2,1,-mm1-1,mm2), &
!& tripleYlmint(ll1,ll2,1,-mm1+1,mm2)
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,prevec2
      do ir=1,gridsize(itpaw,iphigrid(itpaw))-1
        rr=(rphigrid(itpaw,ir+1)+rphigrid(itpaw,ir))/2.d0
        dr=rphigrid(itpaw,ir+1)-rphigrid(itpaw,ir)
        phival1=(phi(itpaw,nn1,ir+1)+phi(itpaw,nn1,ir))/2.d0
        phival2=(phi(itpaw,nn2,ir+1)+phi(itpaw,nn2,ir))/2.d0
        tphival1=(tphi(itpaw,nn1,ir+1)+tphi(itpaw,nn1,ir))/2.d0
        tphival2=(tphi(itpaw,nn2,ir+1)+tphi(itpaw,nn2,ir))/2.d0
        dphival1=(phi(itpaw,nn1,ir+1)-phi(itpaw,nn1,ir))/dr
        dphival2=(phi(itpaw,nn2,ir+1)-phi(itpaw,nn2,ir))/dr
        dtphival1=(tphi(itpaw,nn1,ir+1)-tphi(itpaw,nn1,ir))/dr
        dtphival2=(tphi(itpaw,nn2,ir+1)-tphi(itpaw,nn2,ir))/dr
        radintv1=radintv1+dr*(phival1*dphival2 - phival2*dphival1)
        radintv2=radintv2+dr*phival1*phival2/rr
        radintt1=radintt1+dr*(tphival1*dtphival2 - tphival2*dtphival1)
        radintt2=radintt2+dr*tphival1*tphival2/rr
!write(18,'(i5,5es14.6)') ir,rr,bessj,phival1,phival2,radintv
      enddo
!write(6,'(3i3)') iorb1,iorb2
!write(6,'(3i3)') nn1,nn2
!write(6,'(3i3)') ll1,ll2,lltd
!write(6,'(3i3)') mm1,mm2,mmt
!write(6,*) tripleYlmint(ll1,ll2,lltd,-mm1,mm2)
!write(6,*) radintv

      vectermv = radintv1*prevec1 + radintv2*prevec2
      vectermt = radintt1*prevec1 + radintt2*prevec2
! DEBUG
!      vectermv = radintv2*prevec2
!      vectermt = radintt2*prevec2
!      vectermv = radintv1*prevec1
!      vectermt = radintt1*prevec1
! END DEBUG
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,(vectermv-vectermt)
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,vectermv
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,vectermt
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,(radintv1*prevec1)
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,(radintv2*prevec2)
!write(18,'(i2,2x,3i2,2x,3i2,3x,3("(",f7.3,",",f7.3,")",2x))') itpaw,nn1,ll1,mm1,nn2,ll2,mm2,prevec1
!write(18,'(21x,3("(",f7.3,",",f7.3,")",2x))') prevec2

      do ii=1,3
        pwjmatel(ii,itpaw,iorb1,iorb2)=dot_product(rprimd(ii,:),vectermv)/(2.d0*pi)
        tpwjmatel(ii,itpaw,iorb1,iorb2)=dot_product(rprimd(ii,:),vectermt)/(2.d0*pi)
      enddo
    endif
  enddo
  enddo
enddo

!do iorb1=1,nlmnmax
!  write(6,'(14f6.2)') dble(pwmatel(iorb1,:))
!enddo
!write(6,*)
!do iorb1=1,nlmnmax
!  write(6,'(14f6.2)') dble(tpwmatel(iorb1,:))
!enddo

return
end subroutine fpwjmatel

