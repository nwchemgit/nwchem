subroutine mkbser(nbsem,nbsev,nbsec,nqpt,ncg,nkpt,npwt,igmx,igmn, &
& nsym,natom,ngkpt,ncband,nlmn,npwup,symrel,syminv,ihlf,kg, &
& igndx,lvtrans,indxkpw,indxkbnd,indxkcg,npwarr,pwupmin,pwupmax, &
& vol,pi,xred,kpt,shiftk,cg,pwmatel,tpwmatel,projwf,bser)
implicit none
integer :: nsem,nbsev,nbsec,nqpt,ncg,nkpt,npwt,igmx(3),igmn(3),nsym,natom
integer :: ngkpt(3),ncband,nlmn,npwup
integer :: symrel(3,3,nsym),syminv(3,3,nsym),ihlf(nkpt),kg(3,npwt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
integer :: pwupmin(3),pwupmax(3)
double precision :: vol,pi,xred(3,natom)
double precision :: kpt(3,nkpt),shiftk(3)
double complex :: cg(ncg)
double complex :: pwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3)), &
&                tpwmatel(nlmn,nlmn,npwup,ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: projwf(natom,nlmn,nkpt,ncband)
double precision :: bser(nbsem,nbsem)
integer :: iv,jv,ic,jc,iqpt,jqpt,ikk(3),jkk(3),ikpt,jkpt,igg(3),iqq(3),iqqpt
integer :: ivpw(3),ipwup,ipw,iipw,ii,jj,isym,jsym
double precision :: xkk(3),ykk(3),qq(3),qg2
double complex :: xidir,xix,cmatel,amatel,xmatel1,xmatel2,xmatel3,xmatel4

 igg0=(/0,0,0/)
 qq0=(/0.d0,0.d0,0.d0/)

 do ibsem=1,nbsem  ! ibsem=(iv-nbocc+nbsec-1)*nqpt*nbsec+(ic-nbocc-1)*nqpt+iqpt
   iv=(ibsem-1)/(nqpt*nbsec)+nbocc+1-nbsev
   ic=mod((ibsem-1)/nqpt,nbsec)+nbocc+1
   iqpt=mod(ibsem-1,nqpt)+1
   xkk=qpt(:,iqpt)
   ikk=nint(xkk*ngkpt)+ngkpt/2
   ikpt=ikndx(ikk(1),ikk(2),ikk(3))
   isym=isymndx(ikk(1),ikk(2),ikk(3))
   do jbsem=1,nbsem
     jv=(jbsem-1)/(nqpt*nbsec)+nbocc+1-nbsev
     jc=mod((jbsem-1)/nqpt,nbsec)+nbocc+1
     jqpt=mod(jbsem-1,nqpt)+1
     ykk=qpt(:,jqpt)
     jkk=nint(ykk*ngkpt)+ngkpt/2
     iqq=mod((ikk-jkk)-1,ngkpt)+1
     qq=dble(iqq)/dble(ngkpt)-0.5d0   ! xkk-ykk, shifted to 1st BZ
     igg=(ikk-jkk-iqq)/ngkpt  ! plane wave shift of qq to 1st BZ
     iqqpt=iqndx(iqq(1),iqq(2),iqq(3))
     jkpt=ikndx(jkk(1),jkk(2),jkk(3))
     jsym=isymndx(jkk(1),jkk(2),jkk(3))
     xidir=0.d0 ! direct matrix element
     xix=0.d0   ! exchange matrix element
     do ipw=1,npw
       ivpw=kg(:,indxkpw(ikpt)+ipw)  ! plane wave vector
       if (ivpw(1).ge.pwupmin(1).and.ivpw(1).le.pwupmax(1).and.    &
&          ivpw(2).ge.pwupmin(2).and.ivpw(2).le.pwupmax(2).and.    &
&          ivpw(3).ge.pwupmin(3).and.ivpw(3).le.pwupmax(3)) then
         ipwup=invpwndx(ivpw(1),ivpw(2),ivpw(3))
       else
         ipwup=0
       endif
       if (ipwup.ne.0) then
         qg2=0.d0
         do ii=1,3
         do jj=1,3
           qg2=qg2+(ivpw(ii)+qq(ii))*bmet(ii,jj)*(ivpw(jj)+qq(jj))
         enddo
         enddo
         iipw=invpw2ndx(ipwup,ipwup)
         call mkmatelX(iv,jv,ikpt,jkpt,igpw,igg, &
&          ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&          isym,jsym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&          lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&          cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&          cmatel)
!         call mkmatelP(pi,xred,natom,iv,jv,ikpt,jkpt,qq,igg,ngkpt, &
!&          pwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          tpwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel)
         xmatel1=cmatel+amatel
         call mkmatelX(ic,jc,ikpt,jkpt,igpw,igg, &
&          ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&          isym,jsym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&          lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&          cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&          cmatel)
!         call mkmatelP(pi,xred,natom,ic,jc,ikpt,jkpt,qq,igg,ngkpt, &
!&          pwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          tpwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel)
         xmatel2=cmatel+amatel
         xidir=xidir+lossfn(1,iqqpt,iipw)*xmatel1*xmatel2*4.d0*pi/qg2
       endif
       if (ipw.ne.1) then
         qg2=0.d0
         do ii=1,3
         do jj=1,3
           qg2=qg2+ivpw(ii)*bmet(ii,jj)*ivpw(jj)
         enddo
         enddo
         call mkmatelX(iv,ic,ikpt,ikpt,igpw,igg0, &
&          ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&          isym,jsym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&          lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&          cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&          cmatel)
!         call mkmatelP(pi,xred,natom,iv,ic,ikpt,ikpt,qq0,igg0,ngkpt, &
!&          pwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          tpwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel)
         xmatel3=cmatel+amatel
         call mkmatelX(jv,jc,jkpt,jkpt,igpw,igg0, &
&          ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&          isym,jsym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&          lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&          cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&          cmatel)
!         call mkmatelP(pi,xred,natom,jv,jc,jkpt,jkpt,qq0,igg0,ngkpt, &
!&          pwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          tpwmatel(:,:,ipwup,iqq(1),iqq(2),iqq(3)), &
!&          projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel)
         xmatel4=cmatel+amatel
         xix=xix+4.d0*pi*xmatel3*xmatel4/qg2
       endif
     enddo
     bser(ibsem,jbsem)=xix+xidir
   enddo
   bser(ibsem,ibsem)=bser(ibsem,ibsem)+ &
&         (enrgy(indxkbnd(ikpt)+ic)-enrgy(indxkbnd(ikpt)+iv))
 enddo

return
end subroutine mkbser
