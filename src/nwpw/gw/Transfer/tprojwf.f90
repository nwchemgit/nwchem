subroutine ftprojwf(iat,ncband,nbcore,ihlf,igndx,igmn,igmx,ikndx,isymndx,syminv,isymg)
use wfkvars
use pawvars
use geometry
implicit none
integer :: ncband,nbcore,igmn(3),igmx(3),itpaw
integer :: ihlf(nkpt)
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: syminv(3,3,nsym),isymg(3,nkpt,nsym)
integer :: igndx (nkpt ,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex, allocatable :: ylm(:)
double complex :: phase,twopii,xfact
integer :: igg(3),ikk(3),nn,ll,mm,ix,iy,iz,ikv,ikpt,iband, &
&  ii,jj,ir,iorb,ipw,ipw2,nnold,lm,iat,isym
double precision :: kk(3),kkgg(3),kkggmag
double precision :: rr(gridsizemax-1),dr(gridsizemax-1)
double precision :: radint,projval,pi,sqrt4piinv,prefact
double precision, allocatable :: bessj(:,:)
parameter (pi=3.1415926535897932384626433832795028841972d0)
double precision :: sphbesselj
double complex :: spharm
external sphbesselj,spharm
double complex, allocatable :: projwfg(:,:,:)

itpaw=typat(iat)
allocate(ylm((llmax(itpaw)+1)**2))
allocate(bessj(gridsize(itpaw,iprojgrid(itpaw))-1,0:llmax(itpaw)))
sqrt4piinv=1.d0/sqrt(4.d0*pi)
prefact=4.d0*pi/sqrt(vol)
twopii=(0.d0,2.d0)*pi
do ir=1,gridsize(itpaw,iprojgrid(itpaw))-1
  rr(ir)=(rprojgrid(itpaw,ir+1)+rprojgrid(itpaw,ir))/2.d0
  dr(ir)=rprojgrid(itpaw,ir+1)-rprojgrid(itpaw,ir)
enddo

open(unit=42,file="projwf.log",status="unknown")
write(42,'(a)') "calculated orbital and band projection for: atom    k-point"

do ix=1,ngkpt(1)
do iy=1,ngkpt(2)
do iz=1,ngkpt(3)
  ikk=(/ix,iy,iz/)
  ikv=ix + (iy-1)*ngkpt(1) + (iz-1)*ngkpt(2)*ngkpt(1)
  ikpt=ikndx(ikk(1),ikk(2),ikk(3))
  isym=isymndx(ikk(1),ikk(2),ikk(3))
  do ii=1,3
    kk(ii)=dble(ikk(ii))/dble(ngkpt(ii))-0.5
  enddo
  do iorb=1,nlmn(itpaw)
  do iband=nbcore+1,ncband
    projwf(iat,iorb,ikv,iband)=(0.d0,0.d0)
  enddo
  enddo
  !========================
  allocate(projwfg(nbcore+1:ncband,nlmn(itpaw),npwarr(ikpt)))
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP& PRIVATE(ipw,ii,jj,igg,phase,kkgg,kkggmag,nnold,ll,mm,ylm,ir,bessj,iorb,nn,projval,radint,xfact,ipw2)
  !$OMP DO SCHEDULE(GUIDED)
  do ipw=1,npwarr(ikpt)
!    igg=kg(:,indxkpw(ikpt)+ipw)
    igg=-isymg(:,ikpt,isym)
    do ii=1,3
    do jj=1,3
      igg(ii) = igg(ii)+symrel(jj,ii,mod(isym-1,nsym)+1)*kg(jj,indxkpw(ikpt)+ipw)
    enddo
    enddo
    phase=exp(twopii*dot_product(kk+igg,xred(:,iat)))
    kkgg=(/0.d0,0.d0,0.d0/)
    do ii=1,3
    do jj=1,3
      kkgg(jj)=kkgg(jj)+blat(ii,jj)*(kk(ii)+dble(igg(ii)))
    enddo
    enddo
    kkggmag=sqrt(dot_product(kkgg,kkgg))
    nnold=-1
    do ll=0,llmax(itpaw)
      do mm=-ll,ll
        if (kkgg(1).eq.0.d0.and.kkgg(2).eq.0.d0.and.kkgg(3).eq.0.d0) then
          if (ll.eq.0) then
            ylm(ll**2+ll+mm+1)=sqrt4piinv*(1.d0,0.d0)
          else
            ylm(ll**2+ll+mm+1)=(0.d0,0.d0)
          endif
        else
          ylm(ll**2+ll+mm+1)=spharm(ll,-mm,kkgg)
        endif
      enddo
      do ir=1,gridsize(itpaw,iprojgrid(itpaw))-1
        bessj(ir,ll)=sphbesselj(kkggmag*rr(ir),ll)
      enddo
    enddo
    do iorb=1,nlmn(itpaw)
      nn=iorbno(itpaw,iorb)
      ll=llorb(itpaw,iorb)
      mm=mmorb(itpaw,iorb)
      if (nn.ne.nnold) then
        radint=0.d0
        do ir=1,gridsize(itpaw,iprojgrid(itpaw))-1
          projval=(tproj(itpaw,nn,ir+1)+tproj(itpaw,nn,ir))/2.d0
          radint=radint+dr(ir)*projval*bessj(ir,ll)*rr(ir)
        enddo
      endif
      nnold=nn
      xfact=prefact*radint*ylm(ll**2+ll+mm+1)*phase*(0.d0,1.d0)**ll*(-1)**mm
      do iband=nbcore+1,ncband
        !projwf(iat,iorb,ikv,iband)=projwf(iat,iorb,ikv,iband) &
!&            +xfact*cg(indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw)
        projwfg(iband,iorb,ipw)=xfact*cg(indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw)
      enddo
    enddo
    if (ihlf(ikpt).eq.1) then
      igg=-isymg(:,ikpt,isym)
      do ii=1,3
      do jj=1,3
        igg(ii) = igg(ii)+syminv(ii,jj,mod(isym-1,nsym)+1)*(-kg(jj,indxkpw(ikpt)+ipw)-nint(2.d0*kk(jj)))
      enddo
      enddo
      ipw2=igndx(ikpt,igg(1),igg(2),igg(3))
      if (ipw2.ne.0) cycle
      phase=exp(twopii*dot_product(kk+igg,xred(:,iat)))
      kkgg=(/0.d0,0.d0,0.d0/)
      do ii=1,3
      do jj=1,3
        kkgg(jj)=kkgg(jj)+blat(ii,jj)*(kk(ii)+dble(igg(ii)))
      enddo
      enddo
      kkggmag=sqrt(dot_product(kkgg,kkgg))
      nnold=-1
      do ll=0,llmax(itpaw)
        do mm=-ll,ll
          if (kkgg(1).eq.0.d0.and.kkgg(2).eq.0.d0.and.kkgg(3).eq.0.d0) then
            if (ll.eq.0) then
              ylm(ll**2+ll+mm+1)=sqrt4piinv*(1.d0,0.d0)
            else
              ylm(ll**2+ll+mm+1)=(0.d0,0.d0)
            endif
          else
            ylm(ll**2+ll+mm+1)=spharm(ll,-mm,kkgg)
          endif
        enddo
        do ir=1,gridsize(itpaw,iprojgrid(itpaw))-1
          bessj(ir,ll)=sphbesselj(kkggmag*rr(ir),ll)
        enddo
      enddo
      do iorb=1,nlmn(itpaw)
        nn=iorbno(itpaw,iorb)
        ll=llorb(itpaw,iorb)
        mm=mmorb(itpaw,iorb)
        if (nn.ne.nnold) then
          radint=0.d0
          do ir=1,gridsize(itpaw,iprojgrid(itpaw))-1
            projval=(tproj(itpaw,nn,ir+1)+tproj(itpaw,nn,ir))/2.d0
            radint=radint+dr(ir)*projval*bessj(ir,ll)*rr(ir)
          enddo
        endif
        nnold=nn
        xfact=prefact*radint*ylm(ll**2+ll+mm+1)*phase*(0.d0,1.d0)**ll*(-1)**mm
        do iband=nbcore+1,ncband
          !projwf(iat,iorb,ikv,iband)=projwf(iat,iorb,ikv,iband) &
!&             +xfact*conjg(cg(indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw))
          projwfg(iband,iorb,ipw)=projwfg(iband,iorb,ipw)&
          +xfact*cg(indxkcg(ikpt)+(iband-1)*npwarr(ikpt)+ipw)
        enddo
      enddo
    endif
  enddo
  !$OMP END DO 
  !$OMP END PARALLEL
  do iorb=1,nlmn(itpaw)
  do iband=nbcore+1,ncband
    projwf(iat,iorb,ikv,iband)=sum(projwfg(iband,iorb,:))
  enddo
  enddo
  deallocate(projwfg)
  !========================
  write(42,'(45x,i4,5x,i4)') iat,ikpt
enddo
enddo
enddo

close(unit=42,status="delete")
deallocate(bessj)
deallocate(ylm)
return
end

