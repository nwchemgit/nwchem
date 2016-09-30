subroutine setetint(ebkq,qq2,vme2,W,omega,vsign,lqsing,wmax, &
& ipwa,ipwb,ngkpt,iqndx,iqsymndx,npwup, &
& nsym,pwsymndx,invpw2ndx,nwpt,nqpt,ntpwndx, &
& cse1)
! Tetrahedron integration for self energy
use geometry
implicit none
integer :: ipwa,ipwb,ngkpt(3),npwup,nsym,nwpt,nqpt,ntpwndx
double precision :: ebkq(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: qq2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vme2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: W(nwpt,nqpt,ntpwndx)
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: cse1(-nwpt:nwpt)
double precision :: wmax
integer :: vsign
logical :: lqsing
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))

double precision :: pi,fourpi
double precision :: dw,eps1,eps2,wwlo,wwhi,whi,wlo
integer :: iwwlo,iwwhi,iwlo,iwhi,ielo,iehi
double precision :: de21,de31,de32,de41,de42,de43
double precision :: fbx(4),fb(4),cmx(3,4),cm(3,4),xkt(3,3)
double precision :: thresh
integer :: iqq(3),iqqp(3),iqpt,iqv(3),ictr(3)
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3),nwpt)
double precision :: qkcvt(3,3)
double precision :: enval(8),vqa(8),qqa(8)
double complex :: vfv(8,nwpt)
double precision :: evtx0(4),evtx(4),rrpyr(3,4),vpyr0(4,nwpt),kvtx(3,4), &
& xk(3,4),rg(3,3),tvol,vdum(3),vpyr(4),xpyr0(4),xpyr(4),vqvtx0(4),vqvtx(4)
double precision :: aa0(nwpt),av(3,nwpt)
integer :: indxe(4)
logical :: lqcentr
double precision :: bgrad(3),xmult
double precision :: sint1,sint1a,sint1b,svec(3),sveca(3),svecb(3)
double complex, allocatable :: ssi(:)
!double complex :: ssi(-20*nwpt:20*nwpt)
double complex :: ssc(-nwpt:nwpt)
double precision :: rr(3,8)
integer :: ivndx(4,6)
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
parameter (pi=3.1415926535897932384626433832795, &
& fourpi=12.566370614359172953850573533118d0)
integer :: ii,jj,kk,iv,iw,ix,iy,iz,itet,iww,ie,je
double precision :: www

  dw=wmax/dble(nwpt)
  do ii=1,3
  do jj=1,3
    qkcvt(ii,jj)=blat(ii,jj)/dble(ngkpt(ii))
  enddo
  enddo
  ictr=ngkpt/2

  do ie=-nwpt,nwpt
    ssc(ie)=(0.d0,0.d0)
  enddo
  do ix=1,ngkpt(1)
  do iy=1,ngkpt(2)
  do iz=1,ngkpt(3)
    iqq=(/ix,iy,iz/)
    call fhilo(omega,iqq,ngkpt,whi,wlo)
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
  if (vsign.eq.1) then
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
  call mkvfactor(ipwa,ipwb,ngkpt,iqndx,iqsymndx,npwup,nsym,pwsymndx, &
& invpw2ndx,nwpt,nqpt,ntpwndx,vme2,W,vfactor)

  do ix=1,ngkpt(1)
  do iy=1,ngkpt(2)
  do iz=1,ngkpt(3)
    iqq=(/ix,iy,iz/)
    call fval(omega,iqq,iqqp,ngkpt,enval)
    if (lqsing) then
      call fval(qq2,iqq,iqqp,ngkpt,qqa)
      vqa=fourpi/qqa
    endif
    do iw=1,nwpt
      call fpol(vfactor(:,:,:,iw),ngkpt,iqq,iqqp,vfv(:,iw))
    enddo
    do itet=1,6
      lqcentr=.false.
      do iv=1,4
        evtx0(iv)=enval(ivndx(iv,itet))
        rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
        do iw=1,nwpt
          vpyr0(iv,iw)=vfv(ivndx(iv,itet),iw)
        enddo
        if (lqsing.and..not.lqcentr) then
          iqv=iqq+nint(rrpyr(:,iv))
          if (iqv(1).eq.ictr(1).and.iqv(2).eq.ictr(2).and. &
& iqv(3).eq.ictr(3)) lqcentr=.true.
        endif
      enddo
      call indxhpsort(4,4,evtx0,indxe)
      evtx=evtx0(indxe)
      do ii=1,4
        jj=indxe(ii)
        do kk=1,3
          kvtx(kk,ii)=dot_product(iqq+rrpyr(:,jj)-ictr,qkcvt(:,kk))
        enddo
      enddo
      xk(1:3,1)=(/0.d0,0.d0,0.d0/)
      do ii=2,4
        xk(1:3,ii)=kvtx(1:3,ii)-kvtx(1:3,1)
      enddo
      call cross(xk(:,3),xk(:,4),vdum)
      tvol=dot_product(vdum,xk(:,2))
      do jj=1,3  ! make contragradient
        call cross(xk(1:3,mod(jj,3)+2),xk(1:3,mod(jj+1,3)+2),rg(1:3,jj))
      enddo
      rg=rg/tvol
      tvol=abs(tvol)
      iwwlo=nint(evtx(1)/dw)
      iwwhi=nint(evtx(4)/dw)
      de21=evtx(2)-evtx(1)
      de31=evtx(3)-evtx(1)
      de32=evtx(3)-evtx(2)
      de41=evtx(4)-evtx(1)
      de42=evtx(4)-evtx(2)
      de43=evtx(4)-evtx(3)
      thresh=1.d-5*de41
      if (lqcentr) then            ! integral around coulomb singularity
        do iv=1,4
          vqvtx0(iv)=vqa(ivndx(iv,itet))
        enddo
        bgrad=(/0.d0,0.d0,0.d0/)   ! energy gradient
        do jj=1,3
          bgrad=bgrad+(evtx(jj+1)-evtx(1))*rg(:,jj)
        enddo
        xmult=4.d0*pi/sqrt(dot_product(bgrad,bgrad))
        do iw=1,nwpt
          vpyr(:)=vpyr0(:,iw)/vqvtx0
          aa0(iw)=vpyr(indxe(1))
          av(1:3,iw)=(/0.d0,0.d0,0.d0/)
          do jj=1,3
            av(1:3,iw)=av(1:3,iw) &
&                     +(vpyr(indxe(jj+1))-vpyr(indxe(1)))*rg(1:3,jj)
          enddo
        enddo
        do iww=iwwlo,iwwhi
          www=dw*iww
          if (www.lt.evtx(1)) then
            cycle
          elseif (www.lt.evtx(2)) then
            call singint(1,www,xk,kvtx(:,1),evtx,sint1,svec)
          elseif (www.lt.evtx(3)) then
            if (de21.lt.thresh.and.de43.lt.thresh) then
              call singint2(www,xk,kvtx(:,1),evtx,sint1a,sveca)
            elseif (de21.ge.de43) then
              call singint(2,www,xk,kvtx(:,1),evtx,sint1a,sveca)
              call singint(1,www,xk,kvtx(:,1),evtx,sint1b,svecb)
            else
              call singint(3,www,xk,kvtx(:,1),evtx,sint1a,sveca)
              call singint(4,www,xk,kvtx(:,1),evtx,sint1b,svecb)
            endif
            sint1=sint1b-sint1a
            svec=svecb-sveca
          elseif (www.lt.evtx(4)) then
            call singint(4,www,xk,kvtx(:,1),evtx,sint1,svec)
          else
            cycle
          endif
          sint1=sint1*xmult
          svec=svec*xmult
          do iw=1,nwpt
            ie=iww-vsign*iw
            ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
          enddo
        enddo
      else
        fbx(1)=tvol/(2.d0*de21*de31*de41)
        if (de21.ge.de43) then
          fbx(2)=tvol/(2.d0*de21*de32*de42)
        else
          fbx(3)=tvol/(2.d0*de31*de32*de43)
        endif
        fbx(4)=tvol/(2.d0*de41*de42*de43)
        do ii=1,4
          cmx(1:3,ii)=(/0.d0,0.d0,0.d0/)
          do jj=1,4
            if (jj.ne.ii) cmx(1:3,ii)=cmx(1:3,ii)+(xk(1:3,jj)-xk(1:3,ii)) &
&                                                /(3*(evtx(jj)-evtx(ii)))
          enddo
        enddo
        do iw=1,nwpt
          aa0(iw)=vpyr0(indxe(1),iw)   ! base value of function
          av(1:3,iw)=(/0.d0,0.d0,0.d0/) ! gradient of function
          do jj=1,3
            av(1:3,iw)=av(1:3,iw) &
&                     +(vpyr0(indxe(jj+1),iw)-vpyr0(indxe(1),iw))*rg(1:3,jj)
          enddo
        enddo
        do iww=iwwlo,iwwhi
          www=dw*iww
          if (www.lt.evtx(1)) then
            cycle
          elseif (www.lt.evtx(2)) then
            sint1=fbx(1)*(www-evtx(1))**2
            svec=(xk(1:3,1)+(www-evtx(1))*cmx(1:3,1))*sint1
          elseif (www.lt.evtx(3)) then
            if (de21.lt.thresh.and.de43.lt.thresh) then
              xkt(:,1)=xk(:,3)*(www-evtx(1))/de31
              xkt(:,2)=xk(:,4)*(www-evtx(1))/de41
              xkt(:,3)=(xk(:,4)-xk(:,2))*(www-evtx(2))/de42+xk(:,2)
              sint1=tvol*(2*(www-evtx(1))/(de31*de41) &
&                   -(www-evtx(1))**2*(de31+de41)/(de31*de41)**2)/2
              svec=((xkt(:,1)+xkt(:,3))/2.d0)*sint1
            elseif (de21.ge.de43) then
              fb(1)=fbx(1)*(www-evtx(1))**2
              fb(2)=fbx(2)*(www-evtx(2))**2
              sint1=fb(1)-fb(2)
              cm(1:3,1)=xk(1:3,1)+(www-evtx(1))*cmx(1:3,1)
              cm(1:3,2)=xk(1:3,2)+(www-evtx(2))*cmx(1:3,2)
              svec=cm(:,1)*fb(1)-cm(:,2)*fb(2)
            else
              fb(3)=fbx(3)*(www-evtx(3))**2
              fb(4)=fbx(4)*(www-evtx(4))**2
              sint1=fb(4)-fb(3)
              cm(1:3,3)=xk(1:3,3)+(www-evtx(3))*cmx(1:3,3)
              cm(1:3,4)=xk(1:3,4)+(www-evtx(4))*cmx(1:3,4)
              svec=cm(:,4)*fb(4)-cm(:,3)*fb(3)
            endif
          elseif (www.lt.evtx(4)) then
            sint1=fbx(4)*(www-evtx(4))**2
            svec=(xk(1:3,4)+(www-evtx(4))*cmx(1:3,4))*sint1
          else
            cycle
          endif
          do iw=1,nwpt
            ie=iww-vsign*iw
            ssi(ie)=ssi(ie)+(aa0(iw)*sint1+dot_product(av(:,iw),svec))*dw
          enddo
        enddo
      endif
    enddo
  enddo
  enddo
  enddo
  ssi=-vsign*ssi/((2*pi)**3*pi)

  do ie=-nwpt,nwpt
    eps1=dble(ie)*dw-dw/2
    do je=ielo,iehi
      if (ie.eq.je) then
        ssc(ie)=ssc(ie)+(0.d0,pi)*ssi(je)
      else
        eps2=dble(je)*dw-dw/2
        ssc(ie)=ssc(ie)+vsign*ssi(je)*dw/(eps2-eps1)
      endif
    enddo
  enddo
  deallocate (ssi)

  cse1=cse1+ssc

!  cse=(cse1(0)+cse1(1))/2
!  dcse=(cse1(1)-cse1(0))/dw

return
end subroutine setetint

!*****************************************************************************

subroutine mkvfactor(ipwa,ipwb,ngkpt,iqndx,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,nwpt,nqpt,ntpwndx,vme2,W,vfactor)
implicit none
integer :: ipwa,ipwb,ngkpt(3),npwup,nsym,nwpt,nqpt,ntpwndx
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: vme2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: W(nwpt,nqpt,ntpwndx)
double complex :: vfactor(ngkpt(1),ngkpt(2),ngkpt(3),nwpt)
integer iw,ix,iy,iz,iqq(3),iqpt,jjpw

  do iw=1,nwpt
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
      iqq=(/ix,iy,iz/)
      iqpt=iqndx(iqq(1),iqq(2),iqq(3))
      call locateelement(iqq,ipwa,ipwb,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
      if (jjpw.ne.0) then
        vfactor(ix,iy,iz,iw)=vme2(ix,iy,iz)*dimag(W(iw,iqpt,jjpw))
      else
        vfactor(ix,iy,iz,iw)=(0.d0,0.d0)
      endif
    enddo
    enddo
    enddo
  enddo

return
end subroutine mkvfactor

!*****************************************************************************

subroutine seqptsum(ikpt,lqsing,lx,lc,ipwa,ipwb,ngkpt,iqndx,omegap2,xkf,qq2,W, &
& npwup,nsym,iqsymndx,pwsymndx,nbcore,nbocc,ncband,nwpt,nqpt,ntpwndx, &
& wmax,vme2,ame2,falloff,vsign,ebkq,indxkbnd,nkpt,eigen,bantot,invpw2ndx, &
& exch,ppcor,cor)
use geometry
implicit none
logical :: lqsing,lx,lc
integer :: ikpt,ipwa,ipwb,nbcore,nbocc,ncband,nwpt,nqpt,ntpwndx,vsign
integer :: nsym,npwup,nkpt,bantot
integer :: ngkpt(3)
integer :: iqndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: indxkbnd(nkpt)
double precision :: qq2(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: omegap2,xkf,wmax,falloff
integer :: iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)),invpw2ndx(npwup,npwup)
integer :: pwsymndx(npwup,2*nsym)
double complex :: W(nwpt,nqpt,ntpwndx)
double complex :: vme2(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: ame2(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: exch
double complex :: ppcor(nbcore+1:ncband),cor(nbcore+1:ncband)
double complex :: cse(-nwpt:nwpt)
double precision :: ebkq(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: eigen(bantot)

double precision :: exch1
double complex :: ppcor1(nbcore+1:ncband),cor1(nbcore+1:ncband)
double complex :: cse1(-nwpt:nwpt)

double precision :: pi,fourpi
double precision :: ww,wwq,ppfct,omegat,brd,dw,volelmnt
double complex :: cedenom,cterm
double complex :: wppint0(nbcore+1:ncband),wpp2int0(nbcore+1:ncband)
double complex :: wint0(nbcore+1:ncband),w2int0(nbcore+1:ncband)
double complex :: wppint(nbcore+1:ncband),wpp2int(nbcore+1:ncband)
double complex :: wint(nbcore+1:ncband),w2int(nbcore+1:ncband)
integer :: ictr(3),iqpt,jjpw,iqq(3)
integer :: ie,iw,ix,iy,iz,ii,jj,kk
double precision :: singterm
parameter (pi=3.1415926535897932384626433832795, &
& fourpi=12.566370614359172953850573533118d0)

  volelmnt=vol*ngkpt(1)*ngkpt(2)*ngkpt(3)
  dw=wmax/dble(nwpt)
  brd=dw
  ictr=ngkpt/2
  if (lqsing) then
    iqpt=iqndx(ictr(1),ictr(2),ictr(3))
    wwq=sqrt(omegap2+(xkf**2/3.d0)*qq2(ictr(1),ictr(2),ictr(3))  &
&                 +qq2(ictr(1),ictr(2),ictr(3))**2/4.d0)
    ppfct=pi*omegap2/(2.d0*wwq)
    if (lc) then
      call locateelement(ictr,ipwa,ipwb,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
    else
      jjpw=0
    endif
    do ie=nbcore+1,ncband
      omegat=eigen(indxkbnd(ikpt)+ie)-ebkq(ictr(1),ictr(2),ictr(3))
      cedenom=(omegat-vsign*wwq)*(1.d0,0.d0)+vsign*brd*(0.d0,1.d0)
      wppint0(ie)=(ppfct/cedenom)/pi
      wpp2int0(ie)=-(ppfct/cedenom**2)/pi
      wint0(ie)=(0.d0,0.d0)
      w2int0(ie)=(0.d0,0.d0)
      if (jjpw.ne.0) then
        do iw=1,nwpt
          ww=iw*dw
          cedenom=(omegat-vsign*ww)*(1.d0,0.d0)+vsign*brd*(0.d0,1.d0)
          wint0(ie)=wint0(ie)+dw*(dimag(W(iw,iqpt,jjpw))/cedenom)/pi
          w2int0(ie)=w2int0(ie)-dw*(dimag(W(iw,iqpt,jjpw))/cedenom**2)/pi
        enddo
      endif
    enddo
  endif

  exch1=0.d0
  do ie=nbcore+1,ncband
    ppcor1(ie)=(0.d0,0.d0)
    cor1(ie)=(0.d0,0.d0)
  enddo

  do ix=1,ngkpt(1)
  do iy=1,ngkpt(2)
  do iz=1,ngkpt(3)
    iqq=(/ix,iy,iz/)
    iqpt=iqndx(iqq(1),iqq(2),iqq(3))
    wwq=sqrt(omegap2+(xkf**2/3.d0)*qq2(ix,iy,iz)  &
&           +qq2(ix,iy,iz)**2/4.d0)
    ppfct=pi*omegap2/(2.d0*wwq)
    if (lc) then
      call locateelement(iqq,ipwa,ipwb,ngkpt,iqsymndx,npwup,nsym,pwsymndx,invpw2ndx,jjpw)
    else
      jjpw=0
    endif
    do ie=nbcore+1,ncband
      omegat=eigen(indxkbnd(ikpt)+ie)-ebkq(ix,iy,iz)
      cedenom=(omegat-vsign*wwq)*(1.d0,0.d0)+vsign*brd*(0.d0,1.d0)
      wppint(ie)=(ppfct/cedenom)/pi
      wpp2int(ie)=-(ppfct/cedenom**2)/pi
      wint(ie)=(0.d0,0.d0)
      w2int(ie)=(0.d0,0.d0)
      if (jjpw.ne.0) then
        do iw=1,nwpt
          ww=iw*dw
          cedenom=(omegat-vsign*ww)*(1.d0,0.d0)+vsign*brd*(0.d0,1.d0)
          wint(ie)=wint(ie)+dw*(dimag(W(iw,iqpt,jjpw))/cedenom)/pi
          w2int(ie)=w2int(ie)-dw*(dimag(W(iw,iqpt,jjpw))/cedenom**2)/pi
        enddo
      endif
    enddo
    if (.not.lqsing) then
      if (lx) exch1=exch1-dble(vme2(ix,iy,iz))
      ppcor1=ppcor1+vme2(ix,iy,iz)*wppint
      if (lc) cor1=cor1+vme2(ix,iy,iz)*wint
    else
!*** method 1 for dealing with singularity ***
!      if (qq2(ix,iy,iz).lt.falloff**2) then
!        if (ix.eq.ictr(1).and.iy.eq.ictr(2).and.iz.eq.ictr(3)) cycle
!        cterm=ame2(ictr(1),ictr(2),ictr(3))*(1+cos(pi*sqrt(qq2(ix,iy,iz))/falloff))/2.d0
!        if (lx) exch1=exch1-dble(vme2(ix,iy,iz)-cterm*fourpi/qq2(ix,iy,iz))
!        ppcor1=ppcor1+vme2(ix,iy,iz)*wppint-cterm*wppint0*fourpi/qq2(ix,iy,iz)
!        if (lc) cor1=cor1+vme2(ix,iy,iz)*wint-cterm*wint0*fourpi/qq2(ix,iy,iz)
!*** method 2 for dealing with singularity ***
      if (ix.eq.ictr(1).and.iy.eq.ictr(2).and.iz.eq.ictr(3)) then
        cycle
!*** end alternate singularity method
      else
        if (lx) exch1=exch1-dble(vme2(ix,iy,iz))
        ppcor1=ppcor1+vme2(ix,iy,iz)*wppint
        if (lc) cor1=cor1+vme2(ix,iy,iz)*wint
      endif
    endif
  enddo
  enddo
  enddo
  if (lx) exch1=exch1/volelmnt
  ppcor1=ppcor1/volelmnt
  if (lc) cor1=cor1/volelmnt

  if (lqsing) then
!*** method 1 for dealing with singularity ***
!    if (lx) exch1=exch1-ame2(ictr(1),ictr(2),ictr(3))*falloff/pi
!    ppcor1=ppcor1+wppint0*ame2(ictr(1),ictr(2),ictr(3))*falloff/pi
!    if (lc) cor1=cor1+wint0*ame2(ictr(1),ictr(2),ictr(3))*falloff/pi
!*** method 2 for dealing with singularity ***
    singterm=7.44d0*(vbz/(ngkpt(1)*ngkpt(2)*ngkpt(3)))**(-2.d0/3.d0)/volelmnt
    if (lx) exch1=exch1-ame2(ictr(1),ictr(2),ictr(3))*singterm
    ppcor1=ppcor1+wppint0*ame2(ictr(1),ictr(2),ictr(3))*singterm
    if (lc) cor1=cor1+wint0*ame2(ictr(1),ictr(2),ictr(3))*singterm
!*** end alternate singularity method
  endif
  exch=exch+exch1
  ppcor=ppcor+ppcor1
  cor=cor+cor1

return
end subroutine seqptsum
