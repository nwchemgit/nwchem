subroutine mkxclda(iband,ikpt,ngfft,vxc,vol,ihlf,nkpt,igndx,igmx,igmn,natom,exc)
implicit none
integer :: iband,ikpt,ngfft(3),ihlf(nkpt),igmn(3),igmx(3),natom,nkpt
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double precision :: exc,vxc(ngfft(1)*ngfft(2)*ngfft(3)),vol
integer :: ir,ix,iy,iz,iatom
double precision :: rred(3),orbden
double complex :: twfval,pwfval,wfval

  exc=0.d0
  ir=0
  do iz=0,ngfft(3)-1
  do iy=0,ngfft(2)-1
  do ix=0,ngfft(1)-1
    ir=ir+1
    rred=dble((/ix,iy,iz/))/dble(ngfft)
    call pwrwf(rred,iband,ikpt,twfval,ihlf,igndx,igmx,igmn)
!    do iatom=1,natom
!      call prjrwf(rred,iatom,iband,ikpt,pwfval)
!    enddo
!    wfval=twfval+pwfval
    wfval=twfval
    orbden=dble(wfval*conjg(wfval))
    exc=exc+orbden*vxc(ir)
  enddo
  enddo
  enddo
  exc=exc*vol/(ngfft(1)*ngfft(2)*ngfft(3))

return
end subroutine mkxclda

!******************************************************************************

subroutine pwrwf(rr,iband,ikpt,wfval,ihlf,igndx,igmx,igmn)
! find plane wave expansion of wavefunction at position rr
! rr in reduced coordinates
use wfkvars
use geometry
implicit none
integer :: iband,ikpt,ipw,ihlf(nkpt),igg(3),ipw2,igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: phase,wfval,twopii
double precision :: pi,rr(3)
parameter (pi=3.1415926535897932384626433832795028841972d0)

wfval=(0.d0,0.d0)
twopii=(0.d0,2.d0)*pi
do ipw=1,npwarr(ikpt)
  phase=dot_product(rr,kpt(:,ikpt)+kg(:,indxkpw(ikpt)+ipw))*twopii
  wfval=wfval+cg(indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw)*exp(phase)
enddo
wfval=wfval/sqrt(vol)

return
end subroutine pwrwf

!******************************************************************************

subroutine prjrwf(rr,iatom,iband,ikpt,wfval)
! find PAW correction to wavefunction at position rr in reduced coordinates
use wfkvars
use pawvars
use geometry
implicit none
double precision :: rr(3),rat(3),rrel(3),phival,tphival,rrmag,linterp
integer :: iband,ikpt,nn,ll,mm,iorb,ii,jj,iatom,itpaw
double complex :: wfval,ylm,spharm
external linterp,spharm

  itpaw=typat(iatom)
  rrel=(/0.d0,0.d0,0.d0/)
  do ii=1,3
  do jj=1,3
    rrel(ii)=rrel(ii)+rprimd(ii,jj)*(rr(jj)-xred(jj,iatom))
  enddo
  enddo
  rrmag=sqrt(dot_product(rrel,rrel))
  wfval=0.d0
  if (rrmag.gt.rphigrid(itpaw,gridsize(itpaw,iphigrid(itpaw)))) return
  do iorb=1,nlmn(itpaw)
    nn=iorbno(itpaw,iorb)
    ll=llorb(itpaw,iorb)
    mm=mmorb(itpaw,iorb)
    ylm=spharm(ll,mm,rrel)
    tphival=linterp(tphi(itpaw,nn,:),rphigrid(itpaw,:),gridsizemax,rrmag)
    phival=linterp(phi(itpaw,nn,:),rphigrid(itpaw,:),gridsizemax,rrmag)
    wfval=wfval+projwf(iatom,iorb,ikpt,iband)*ylm*(phival-tphival)/rrmag
  enddo

return
end subroutine prjrwf

!******************************************************************************

subroutine prjrtwf(rr,iatom,iband,ikpt,wfval)
! find PAW correction to wavefunction at position rr in reduced coordinates
use wfkvars
use pawvars
use geometry
implicit none
double precision :: rr(3),rat(3),rrel(3),phival,tphival,rrmag,linterp
integer :: iband,ikpt,nn,ll,mm,iorb,ii,jj,iatom,itpaw
double complex :: wfval,ylm,spharm
external linterp,spharm

  itpaw=typat(iatom)
  rrel=(/0.d0,0.d0,0.d0/)
  do ii=1,3
  do jj=1,3
    rrel(ii)=rrel(ii)+rprimd(ii,jj)*(rr(jj)-xred(jj,iatom))
  enddo
  enddo
  rrmag=sqrt(dot_product(rrel,rrel))
  wfval=0.d0
  if (rrmag.gt.rphigrid(itpaw,gridsize(itpaw,iphigrid(itpaw)))) return
  do iorb=1,nlmn(itpaw)
    nn=iorbno(itpaw,iorb)
    ll=llorb(itpaw,iorb)
    mm=mmorb(itpaw,iorb)
    ylm=spharm(ll,mm,rrel)
    tphival=linterp(tphi(itpaw,nn,:),rphigrid(itpaw,:),gridsizemax,rrmag)
    wfval=wfval+projwf(iatom,iorb,ikpt,iband)*ylm*tphival/rrmag
  enddo

return
end subroutine prjrtwf

!******************************************************************************

subroutine prjrpwf(rr,iatom,iband,ikpt,wfval)
! find PAW correction to wavefunction at position rr in reduced coordinates
use wfkvars
use pawvars
use geometry
implicit none
double precision :: rr(3),rat(3),rrel(3),phival,tphival,rrmag,linterp
integer :: iband,ikpt,nn,ll,mm,iorb,ii,jj,iatom,itpaw
double complex :: wfval,ylm,spharm
external linterp,spharm

  itpaw=typat(iatom)
  rrel=(/0.d0,0.d0,0.d0/)
  do ii=1,3
  do jj=1,3
    rrel(ii)=rrel(ii)+rprimd(ii,jj)*(rr(jj)-xred(jj,iatom))
  enddo
  enddo
  rrmag=sqrt(dot_product(rrel,rrel))
  wfval=0.d0
  if (rrmag.gt.rphigrid(itpaw,gridsize(itpaw,iphigrid(itpaw)))) return
  do iorb=1,nlmn(itpaw)
    nn=iorbno(itpaw,iorb)
    ll=llorb(itpaw,iorb)
    mm=mmorb(itpaw,iorb)
    ylm=spharm(ll,mm,rrel)
    phival=linterp(phi(itpaw,nn,:),rphigrid(itpaw,:),gridsizemax,rrmag)
    wfval=wfval+projwf(iatom,iorb,ikpt,iband)*ylm*phival/rrmag
  enddo

return
end subroutine prjrpwf

