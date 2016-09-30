subroutine mkhex(ncond,nval, ,hres)
implicit none

!ncond=2
!nval=4
do ix1=1,nkpt*nval*ncond
do ix2=1,ix1
  ikpt1=mod(ix1-1,nkpt)+1
  ikpt2=mod(ix2-1,nkpt)+1
  ic1=mod(ix1/nkpt-1,ncond)+1+nbocc
  ic2=mod(ix2/nkpt-1,ncond)+1+nbocc
  iv1=mod(ix1/(nkpt*ncond)-1,nval)+1+nbocc-nval
  iv2=mod(ix2/(nkpt*ncond)-1,nval)+1+nbocc-nval
  hres(ix1,ix2)=hdirlmnt(iv1,ic1,ikpt1,iv2,ic2,ikpt2, &
&                 ncg,cg,nkpt,npwarr,indxkcg,igmn(3),igmx(3), &
&                 ipwndx,ntpwndx,ipwup,npwup,igndx) &
&              +hexchlmnt(iv1,ic1,ikpt1,iv2,ic2,ikpt2, &
&                 ncg,cg,nkpt,npwarr,indxkcg,igmn(3),igmx(3),
&                 ipwup,npwup,igndx)
  if (ix1.eq.ix2) then
    hres(ix1,ix2)=hres(ix1,ix2)+enrgy(indxkbnd(ikpt1)+ic1)-enrgy(indxkbnd(ikpt1)+iv1)
  endif
enddo
enddo

return
end subroutine mkhex

!*************************************************************************

function hdirlmnt(iv1,ic1,ikpt1,iv2,ic2,ikpt2, &
& ncg,cg,nkpt,npwarr,indxkcg,igmn(3),igmx(3),ipwndx,ntpwndx,ipwup,npwup,igndx)
implicit none
integer :: iv1,ic1,ikpt1,iv2,ic2,ikpt2
integer :: ncg,nkpt,npwarr(nkpt),indxkcg(nkpt),igmn(3),igmx(3)
integer :: npwup,ipwup(3,npwup),ntpwndx,ipwndx(2,ntpwndx)
integer :: igndx(mkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
double precision :: qq(3)
integer :: iqq(3),iqp(3),iqpt,iqg(3),iipw,ipw1,ipw2,igpw,igg1(3),igg2(3),jgpw1,jgpw2,icg1,jcg1.icg2,jcg2
double complex :: hdirlmnt,amp1,amp2

 qq=kpt(:,ikpt2)-kpt(:,ikpt1)
 iqq=nint(ngkpt*(qq+0.5d0))
 iqp=mod(iqq-1,ngkpt)+1
 iqpt=iqndx(iqp(1),iqp(2),iqp(3))
 iqg=(iqq-iqp)/ngkpt
 ikk=nint(ngkpt*(kpt(ikpt1)+0.5d0))

 hdirlmnt=0.d0
 do iipw=1,ntpwndx
   ipw1=ipwndx(1,iipw)
   ipw2=ipwndx(2,iipw)
   if (ikpt1.eq.ikpt2) then
     call mkmatelX(iv1,iv2,ikpt1,ikpt2,ipw1,iqg, &
&      ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&      isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&      cmatel1)
!     call mkmatelP(pi,xred,natom,iv1,iv2,ikpt1,ikpt2,qq,ipw1,ngkpt, &
!&      pwmatel(:,:,ipw1,iqq(1),iqq(2),iqq(3)), &
!&      tpwmatel(:,:,ipw1,iqq(1),iqq(2),iqq(3)), &
!&      projwf,nlmn,nkpt,nkptq,ncband,amatel1)
     call mkmatelX(iv1,iv2,ikpt1,ikpt2,ipw2,iqg, &
&      ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
&      isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtransq(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
&      cmatel2)
!     call mkmatelP(pi,xred,natom,ic1,ic2,ikpt1,ikpt2,qq,ipw2,ngkpt, &
!&      pwmatel(:,:,ipw2,iqq(1),iqq(2),iqq(3)), &
!&      tpwmatel(:,:,ipw2,iqq(1),iqq(2),iqq(3)), &
!&      projwf,nlmn,nkpt,nkptq,ncband,amatel1)
   else
     call mkmatelX(iv1,iv2,ikpt1,ikpt2,ipw1,iqg, &
&      ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&      isym,isym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&      cmatel1)
!     call mkmatelP(pi,xred,natom,iv1,iv2,ikpt1,ikpt2,qq,ipw1,ngkpt, &
!&      pwmatel(:,:,ipw1,iqq(1),iqq(2),iqq(3)), &
!&      tpwmatel(:,:,ipw1,iqq(1),iqq(2),iqq(3)), &
!&      projwf,nlmn,nkpt,nkpt,ncband,amatel1)
     call mkmatelX(iv1,iv2,ikpt1,ikpt2,ipw2,iqg, &
&      ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&      isym,isym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&      cmatel2)
!     call mkmatelP(pi,xred,natom,ic1,ic2,ikpt1,ikpt2,qq,ipw2,ngkpt, &
!&      pwmatel(:,:,ipw2,iqq(1),iqq(2),iqq(3)), &
!&      tpwmatel(:,:,ipw2,iqq(1),iqq(2),iqq(3)), &
!&      projwf,nlmn,nkpt,nkpt,ncband,amatel1)
   endif
   hdirlmnt = hdirlmnt-(cmatel1+amatel1)*(cmatel2+amatel2)*W(1,iqpt,iipw)
 enddo

return
end function hdirlmnt

!*************************************************************************

function hexchlmnt(iv1,ic1,ikpt1,iv2,ic2,ikpt2, &
& ncg,cg,nkpt,npwarr,indxkcg,igmn(3),igmx(3),ipwup,npwup,igndx)
implicit none
integer :: iv1,ic1,ikpt1,iv2,ic2,ikpt2
integer :: ncg,nkpt,npwarr(nkpt),indxkcg(nkpt),igmn(3),igmx(3)
integer :: npwup,ipwup(3,npwup)
integer :: igndx(mkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
integer :: ipw1,jpw1,iqq(3),ipw2,jpw2,jgg2(3),icg1,jcg1,icg2,jcg2,ii,jj
double precision :: qq2,vqq
double complex :: hexchlmnt

 igg0=(/0,0,0/)
 iqctr=ngkpt/2
 qq=(/0.d0,0.d0,0.d0/)
 hexchlmnt = 0.d0
 do ipw1=1,npwarr(ikpt1)
   iqq=ipwup(:,ipw1)
   if (iqq(1).eq.0.and.iqq(2).eq.0.and.iqq(3).eq.0) cycle
   qq2=0.d0
   do ii=1,3
   do jj=1,3
     qq2=qq2+iqq(ii)*bmet(ii,jj)*iqq(jj)
   enddo
   enddo
   vqq=4.d0*pi/qq2
   call mkmatelX(iv1,ic1,ikpt1,ikpt1,ipw1,igg0, &
&      ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&      isym,isym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&      cmatel1)
!   call mkmatelP(pi,xred,natom,iv1,ic1,ikpt1,ikpt1,qq,ipw1,ngkpt, &
!&      pwmatel(:,:,ipw1,iqctr(1),iqctr(2),iqctr(3)), &
!&      tpwmatel(:,:,ipw1,iqctr(1),iqctr(2),iqctr(3)), &
!&      projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel1)
   call mkmatelX(iv2,ic2,ikpt2,ikpt2,ipw1,igg0, &
&      ncg,ncg,nkpt,nkpt,npwt,npwt,igmx,igmn,igndx,igndx, &
&      isym,isym,symrel,syminv,nsym,nsym,ihlf,ihlf,kpt,kpt, &
&      lvtrans(1:3,ikk(1),ikk(2),ikk(3)),lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&      cg,cg,indxkcg,indxkcg,indxkpw,indxkpw,npwarr,npwarr,kg,kg, &
&      cmatel2)
!   call mkmatelP(pi,xred,natom,iv2,ic2,ikpt2,ikpt2,qq,ipw1,ngkpt, &
!&      pwmatel(:,:,ipw1,iqctr(1),iqctr(2),iqctr(3)), &
!&      tpwmatel(:,:,ipw1,iqctr(1),iqctr(2),iqctr(3)), &
!&      projwf,projwf,nlmn,nkpt,nkpt,ncband,amatel2)
   hexchlmnt = hexchlmnt+(cmatel1+amatel1)*(cmatel2+amatel2)*vqq
 enddo

return
end function hexchlmnt

!*************************************************************************

subroutine felmnt(vec1,vec2,xmat,nn,elmnt)
implicit none
integer :: nn
double complex :: vec1(nn),vec2(nn),xmat(nn,nn),elmnt
double complex :: vv(nn)

 elmnt=(0.d0,0.d0)
 do ix=1,nn
   vv(ix)=0.d0
   do iy=1,nn
     vv(ix)=vv(ix)+xmat(ix,iy)*vec2(iy)
   enddo
   elmnt=elmnt+vec1(ix)*vv(ix)
 enddo

return
end subroutine felement

!*************************************************************************

subroutine pert1(nn,eigvec0,eigval0,xmat1,eigvec1)
implicit none
integer :: nn
double complex :: eigvec0(nn,nn),eigval0(nn),xmat1(nn,nn),eigvec1(nn,nn)

do ix=1,nn
  eigvec1(:,ix)=eigvec0(:,ix)
  do iy=1,nn
    if (ix.eq.iy) cycle
    call felmnt(eigvec0(:,iy),eigvec0(:,ix),xmat1,nn,elmnt)
    eigvec1(:,ix)=eigvec1(:,ix)+eigvec0(:,iy)*elmnt/(eigval(ix)-eigval(iy))
  enddo
enddo

return
end subroutine pert1








































