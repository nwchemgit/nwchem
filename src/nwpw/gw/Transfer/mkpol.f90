subroutine mkpol1(qq,igglf1,igglf2,iks,vq,vol,pi, &
& test_bands_pol, &
& nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
& nsym,nkpt, &
& ncg,bantot,npwt,nlmn, &
& nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
& xred, &
& projwf, &
& kg, &
& enrgy, eigen, &
& cg, &
& indxkpw, &
& indxkbnd, &
& indxkcg, &
& npwarr, &
& kpt, &
& nsymk,symk,symrel,syminv, &
& ihlf,lvtrans, &
& ntypepaw, &
& pwmatel1,tpwmatel1,pwmatel2,tpwmatel2, &
& igndx,ikndx,isymndx,isymg, &
& nband,pola,polb)
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),dq(3),vq,vol,pi,wmax,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw,itetrahedron
integer :: igglf1(3),igglf2(3)
integer :: test_bands_pol(4)
double precision :: shiftk(3)
double precision :: xred(3,natom)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot), eigen(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
integer :: nband(nkpt*nsppol)
! number of bands calculated by abinit
double complex :: pola(nwpt),polb(nwpt),tpolincr,polci(npwt),dpol
! pola - hermitian part of vq * polarization matrix element for igglf1,igglf2
! polb - antihermitian part of vq * polarization matrix element 
double complex :: cmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,cmat2v(8)
! cmatel - Complex density MATrix ELement
! cmatellf - Complex density MATrix ELement for Local Fields
! cmat2 - "square" of matrix elements, cmatel*cmatellf
double complex :: amatel,amatellf
! amatel, amatellf - contribution to matrix elements from PAW
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
double precision :: volelmnt
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: testksym(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double complex :: pwmatel2(ntypepaw,nlmn,nlmn),tpwmatel2(ntypepaw,nlmn,nlmn)
! DEBUG
!double complex :: jmatel(3)
! DEBUG

volelmnt=vol*ngkpt(1)*ngkpt(2)*ngkpt(3)
igg0=(/0,0,0/)
qp=qq-igglf1+igglf2
dw=wmax/dble(nwpt)
do iw=1,nwpt
  polci(iw)=0.d0
  pola(iw)=(0.d0,0.d0)
  polb(iw)=(0.d0,0.d0)
enddo

!write(6,'(a)') 'occupied band, unoccupied band, minimum gap'
do iocc=test_bands_pol(1),test_bands_pol(2)
!do iocc=1,1
  do iunocc=test_bands_pol(3),test_bands_pol(4)
!  do iunocc=8,8
!    write(*,*) iocc,iunocc
    gap=1.d20
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=2,2
!    do iy=9,9
!    do iz=9,9
      ikk=(/ix,iy,iz/)
!write(*,*) ikk
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)  ! integer index of k-q
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
!      write(6,'(4(3i3,3x))') ikk,jkk,iks,igg
!      write(6,'(3(3f10.6,3x))') xck,xckq
!      write(6,'(3(3f10.6,3x))') qq,qq-igglf1
!write(44,*) "band (hole) = ",iunocc
!write(44,*) "band (conduction electron) = ",iocc
!write(44,'(a,3f8.4)') "k   = ",xck
!write(44,'(a,3f8.4)') "k-q = ",xckq
!write(44,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(44,'(a,3f8.4)') "k-q in irreducable Brillouin zone = ",kpt(1:3,ikptq)
!write(44,'(a,3i3,3x)') "integer index of k   = ",ikk
!write(44,'(a,3i3,3x)') "integer index of k-q = ",jkk
!write(44,'(a,3i3,3x)') "integer index of q   = ",iks
!write(44,*) "symmetry operation on k"
!write(44,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!testksym=(/0.d0,0.d0,0.d0/)
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + syminv(jj,ii,mod(isym-1,nsym)+1)*(xck(jj)-igsymk(jj))
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!testksym=igsymk
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + symrel(jj,ii,mod(isym-1,nsym)+1)*kpt(jj,ikpt)
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!write(44,*)
!write(44,*) "symmetry operation on k-q"
!write(44,*) "symmetry index ",isymq,mod(isymq-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!enddo
!testksym=(/0.d0,0.d0,0.d0/)
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + syminv(jj,ii,mod(isymq-1,nsym)+1)*(xckq(jj)-igsymq(jj))
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!testksym=igsymq
!do ii=1,3
!do jj=1,3
!testksym(ii) = testksym(ii) + symrel(jj,ii,mod(isymq-1,nsym)+1)*kpt(jj,ikptq)
!enddo
!enddo
!write(44,'(3f8.4)') testksym
!write(44,*)
!      write(6,'(a,3i5)') "ikpt  = ",ikpt
!      write(6,'(a,3i5)') "ikptq = ",ikptq
!      write(6,'(a,3i3,3x)') "integer index of k   = ",ikk
!      write(6,'(a,3i3,3x)') "integer index of k-q = ",jkk
!      write(6,*) "k-points used"
!      write(6,'(a,3f8.4)') "k   = ",xck
!      write(6,'(a,3f8.4)') "k-q = ",xckq
!      write(6,*) "k-points symmetry transform to irreducable Brillouin zone"
!      write(6,'(a,3f8.4)') "k   = ",kpt(1:3,ikpt)
!      write(6,'(a,3f8.4)') "k-q = ",kpt(1:3,ikptq)
!!      write(6,'(3f8.4,4x,3f8.4)') kpt(1:3,ikpt),kpt(1:3,ikptq)
!!      write(6,'(3i5)') isym,isymq,nsym
!      write(6,'(a,3i5)') "isym  = ",isym
!      write(6,'(a,3i5)') "isymq = ",isymq
!      write(6,'(2i4,4x,2i4)') mod(isym-1,nsym)+1,-(2*((isym-1)/nsym)-1),mod(isymq-1,nsym)+1,-(2*((isymq-1)/nsym)-1)
!      write(6,*)
!      write(6,*) "symmetry: k"
!      testksym=(/0.d0,0.d0,0.d0/)
!      do ii=1,3
!        write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!        do jj=1,3
!          testksym(ii) = testksym(ii) + symrel(ii,jj,mod(isym-1,nsym)+1)*xck(jj)
!        enddo
!      enddo
!      write(6,'(3f8.4)') testksym
!      write(6,*)
!      write(6,*) "symmetry: k-q"
!      testksym=(/0.d0,0.d0,0.d0/)
!      do ii=1,3
!        write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!        do jj=1,3
!          testksym(ii) = testksym(ii) + symrel(ii,jj,mod(isymq-1,nsym)+1)*xckq(jj)
!        enddo
!      enddo
!      write(6,'(3f8.4)') testksym
!      do ii=1,npwarr(ikptq)
!        write(91,'(i5,i10,2f12.6)') ii,indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+ii,cg(indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+ii)
!      enddo
      call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&       ncg,nkpt,npwt,igmx,igmn,igndx, &
&       isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&       lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&       cg,indxkcg,indxkpw,npwarr,kg, &
&       cmatel)
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      if (ipaw.ne.0) then
        call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
& pwmatel1,tpwmatel1,nbcore,ncband,amatel)
      else
        amatel=(0.d0,0.d0)
      endif
! DEBUG
!      call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
!&          ncg,nkpt,npwt,igmx,igmn,igndx, &
!&          isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
!&          lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
!&          cg,indxkcg,indxkpw,npwarr,kg, &
!&          jmatel)
!      if (ipaw.ne.0) then
!        call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
!&                       pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
!        jmatel(:) = jmatel(:) + vmatel(:)
!      endif
!!      if (ipaw.ne.0) then
!!        call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
!!& pwmatel1,tpwmatel1,nbcore,ncband,amatel)
!!      else
!!        amatel=(0.d0,0.d0)
!!      endif
! DEBUG
!write(6,*) jkk
!write(6,*) ihlf(ikpt),ihlf(ikptq)
!write(6,*) cmatel
!write(6,*) amatel
!write(6,*) cmatel+amatel
!write(6,*) (cmatel+amatel)*conjg(cmatel+amatel)
!write(6,*) (cmatel+amatel)*conjg(cmatel+amatel)*vq
!write(6,*) 
      if (igglf1(1).eq.igglf2(1).and.igglf1(2).eq.igglf2(2).and. &
&         igglf1(3).eq.igglf2(3)) then
        cmat2(ikk(1),ikk(2),ikk(3))=dble((cmatel+amatel)*conjg(cmatel+amatel))
      else
        call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatellf)
        if (ipaw.ne.0) then
          call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qp,igg0, &
& pwmatel2,tpwmatel2,nbcore,ncband,amatellf)
        else
          amatellf=(0.d0,0.d0)
        endif
        cmat2(ikk(1),ikk(2),ikk(3))=(cmatel+amatel)*conjg(cmatellf+amatellf)
!write(6,*) cmatellf
!write(6,*) amatellf
!write(6,*) cmatellf+amatellf
!write(6,*) 
      endif
      if (omega(ikk(1),ikk(2),ikk(3)).lt.gap) gap=omega(ikk(1),ikk(2),ikk(3))
!write(47,'(1x,3i3,f10.6,e16.6)') ikk, omega(ikk(1),ikk(2),ikk(3)), dble(cmat2(ikk(1),ikk(2),ikk(3)))
!write(6,*) cmat2(ikk(1),ikk(2),ikk(3)),omega(ikk(1),ikk(2),ikk(3))
!write(46,'(2i3,2x,3i3,2x,f7.4,1x,"(",e10.3,",",e10.3,")",2x,2i3)') &
!&  iocc,iunocc,ikk(1:3),omega(ikk(1),ikk(2),ikk(3)), &
!&  cmat2(ikk(1),ikk(2),ikk(3)) !, isym, isymq
!write(6,*) dble(cmat2(ikk(1),ikk(2),ikk(3)))
!write(6,*) dble(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!write(6,*) 
!      stvec(ikk(3))=enrgy(indxkbnd(ikpt)+iunocc)
!      stvec(ikk(3))=omega(ikk(1),ikk(2),ikk(3))-dw*58
!      stvec(ikk(3))=dble(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!      stvec2(ikk(3))=dimag(cmat2(ikk(1),ikk(2),ikk(3)))*vq
!      stvec(ikk(3))=dble(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec2(ikk(3))=dimag(cmat2(ikk(1),ikk(2),ikk(3)))
!      stvec(ikk(3))=abs(cmatel)
!      stvec(ikk(3))=dble(cmatel)
!      stvec2(ikk(3))=dimag(cmatel)
!      stvec(ikk(3))=dble(cmatellf)
!      stvec2(ikk(3))=dimag(cmatellf)
!      stvec(ikk(3))=dble(cmatel*conjg(cmatellf))*vq
!      stvec(ikk(3))=dble(amatel*conjg(amatellf))*vq
    enddo
!    write(68,'(10f8.3)') stvec(1:ngkpt(3)),stvec2(1:ngkpt(3))
!    write(68,'(10f8.3)') stvec(1:ngkpt(3))
!    write(68,'(10f8.3)') stvec2(1:ngkpt(3))
!    write(68,*)
!    write(68,'(10f8.2)') stvec(1:ngkpt(3))*27.2113845d0
!    write(68,'(10es8.1)') stvec(1:ngkpt(3))
!    write(68,'(10es8.1)') stvec2(1:ngkpt(3))
!    write(68,*)
    enddo
!    write(68,*) ix+1,'------------------------------------------------------'
!    write(68,*) 
    enddo
!    write(6,'(2(i6,9x),f10.3)') iocc,iunocc,gap*27.2113845d0
!stop

    if (itetrahedron.ne.0) then
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=10,10
!      do iy=10,10
!      do iz=10,10
        ikk=(/ix,iy,iz/)
        call fhilo(omega,ikk,ngkpt,whi,wlo)
        iwh=min(nint(whi*dble(nwpt)/wmax),nwpt)
        iwl=nint(wlo*dble(nwpt)/wmax)
        if (iwl.gt.nwpt) cycle
write(47,'(1x,2i3,2x,3i3,2f10.6,2i5)') iocc, iunocc, ikk, whi, wlo, iwh, iwl
        call fval(omega,ikk,ikkp,ngkpt,enval)
        call fpol(cmat2,ngkpt,ikk,ikkp,cmat2v)
        do iw=iwl,iwh
!        do iw=166,166
          ww=iw*dw
!          call cubeint(enval,ww,cmat2v,tpolincr)
!          dpol=tpolincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
          call vcubeint(enval,ww-dw/2.d0,dw,cmat2v,tpolincr)
          dpol=tpolincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*dw)
          polci(iw)=polci(iw)-dpol*2.d0*vq
!write(6,'(5i5,4e12.4)') ix,iy,iz,iw,nwpt,dpol*2.d0*vq
!if (iw.eq.166) then
!  if (ix.eq.4.and.iy.eq.4.and.iz.eq.4) then
!    write(27,'(2i5)') iocc,iunocc
!    write(27,'(5i5,4e12.4)') ix,iy,iz,iw,nwpt,dpol*8.d0*pi
!    write(27,'(2f12.6)') ww-dw/2.d0,ww+dw/2.d0
!    write(27,'(2f9.6,4x,2f9.6)') enval(1),enval(2),enval(5),enval(6)
!    write(27,'(2f9.6,4x,2f9.6)') enval(3),enval(4),enval(7),enval(8)
!    write(27,*)
!    write(27,'(2f9.6,4x,2f9.6)') dble(cmat2v(1)),dble(cmat2v(2)),dble(cmat2v(5)),dble(cmat2v(6))
!    write(27,'(2f9.6,4x,2f9.6)') dble(cmat2v(3)),dble(cmat2v(4)),dble(cmat2v(7)),dble(cmat2v(8))
!    write(27,*)
!    write(27,*)
!  endif
!endif
!          write(77,'(f12.6,2e20.6)') ww,dpol
        enddo
!        stvec(ikk(3))=dble(iwh)
!        stvec(ikk(3))=dble(dpol*2.d0*vq)
      enddo
!      write(68,'(10f8.3)') stvec(1:ngkpt(3))
      enddo
!      write(68,*) ix+1,'------------------------------------------------------'
!      write(68,*)
      enddo
    else
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
        iw=int(omega(ix,iy,iz)/dw)
        if (iw.le.nwpt) then
          fw=omega(ix,iy,iz)/dw-iw
          dpol=2.d0*vq*cmat2(ix,iy,iz)/(pi*vol)
!write(6,'(5i5,4e12.4)') ix,iy,iz,iw,nwpt,dpol
          polci(iw)=polci(iw)-fw*dpol
          if (iw.lt.nwpt) polci(iw+1)=polci(iw+1)-(1.d0-fw)*dpol
        endif
      enddo
      enddo
      enddo
    endif
  enddo
enddo

! polci gives "imaginary" part of polarization - i.e., the anti-hermetian part when dealing with non-diagonal components.
  do iw=1,nwpt
    ww1=dble(iw)*dw
    do jw=1,nwpt
      if (iw.eq.jw) cycle
      ww2=dble(jw)*dw
      pola(iw)=pola(iw)+2.d0*polci(jw)*dw*(ww2/(ww2**2-ww1**2))
    enddo
    polb(iw)=-(0.d0,1.d0)*pi*polci(iw)
write(48,*) iw,ww1,dble(pola(iw)),dimag(polb(iw))
  enddo

! debug
!  do iw=1,nwpt
!    write(27,'(i5,3x,4e12.4)') iw,polci(iw)
!!    write(27,'(i5,3x,4e12.4)') iw,pola(iw),polb(iw)
!  enddo
! debug

return
end subroutine mkpol1

!*************************************************************************

subroutine mkpolj1(qq,igglf1,igglf2,iks,vq,pi, &
& test_bands_pol, &
& nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
& nsym,nkpt, &
& ncg,bantot,npwt,nlmn, &
& nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
& xred, &
& projwf, &
& kg, &
& enrgy, eigen, &
& cg, &
& indxkpw, &
& indxkbnd, &
& indxkcg, &
& npwarr, &
& kpt, &
& nsymk,symk,symrel,syminv, &
& ihlf,lvtrans, &
& ntypepaw, &
& pwjmatel,tpwjmatel,pwmatel1,tpwmatel1,pwmatel2,tpwmatel2, &
& igndx,ikndx,isymndx,isymg, &
& nband,polat,polbt)
use geometry
implicit none
integer :: bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),qp(3),qq2,qp2,dq(3),vq,pi,wmax,dw
integer :: nwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
! iks: integer index of q-vector in Brillouin zone
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw,itetrahedron
integer :: igglf1(3),igglf2(3)
integer :: test_bands_pol(4)
double precision :: shiftk(3)
double precision :: xred(3,natom)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot), eigen(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
integer :: nband(nkpt*nsppol)
! number of bands calculated by abinit
double complex :: polat(nwpt,9),polbt(nwpt,9),tpolincr,polci(npwt,3,3),dpol
! polat - hermitian part of vq * polarization tensor for igglf1,igglf2
! polbt - antihermitian part of vq * polarization tensor
double complex :: cmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,cmat2v(8)
! cmatel - Complex density MATrix ELement
! cmatellf - Complex density MATrix ELement for Local Fields
! cmat2 - "square" of matrix elements, cmatel*cmatellf
double complex :: jmatel(3),jmatellf(3),tmat2(ngkpt(1),ngkpt(2),ngkpt(3),3,3)
! jmatel - current matrix element
! jmatellf - current matrix element for local fields
! tmat2 - tensor product of current matrix elements
double complex :: amatel,vmatel(3)
! amatel - contribution to density matrix elements from PAW
! vmatel - contribution to current matrix elements from PAW
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
double precision :: volelmnt
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3)),gap
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,ll,mm,nn,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq,iiq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwjmatel(ntypepaw,nlmn,nlmn),tpwjmatel(ntypepaw,nlmn,nlmn)
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double complex :: pwmatel2(ntypepaw,nlmn,nlmn),tpwmatel2(ntypepaw,nlmn,nlmn)
double precision :: matdummy1(3,3),matdummy2(3,3),matdummy3(3,3),dummy
double complex :: cmatdummy1(3,3),cmatdummy2(3,3),cmatdummy3(3,3)

volelmnt=vol*ngkpt(1)*ngkpt(2)*ngkpt(3)
igg0=(/0,0,0/)
qp=qq-igglf1+igglf2
do ii=1,3
do jj=1,3
  qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
  qp2=qp2+qp(ii)*bmet(ii,jj)*qp(jj)
enddo
enddo
dw=wmax/dble(nwpt)
do iw=1,nwpt
  do ii=1,3
  do jj=1,3
    polci(iw,ii,jj)=0.d0
  enddo
  enddo
  do ii=1,9
    polat(iw,ii)=(0.d0,0.d0)
    polbt(iw,ii)=(0.d0,0.d0)
  enddo
enddo

!write(6,'(4i5)') test_bands_pol(1),test_bands_pol(2),test_bands_pol(3),test_bands_pol(4)
!write(6,'(a)') 'occupied band, unoccupied band, minimum gap'
do iocc=test_bands_pol(1),test_bands_pol(2)
!do iocc=1,1
  do iunocc=test_bands_pol(3),test_bands_pol(4)
!  do iunocc=8,8
!write(*,*) iocc,iunocc
    gap=1.d20
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=6,6
!    do iy=5,5
!    do iz=5,5
!    do ix=6,6
!    do iy=6,6
!    do iz=6,6
!    do ix=5,5
!    do iy=6,6
!    do iz=5,5
!    do ix=5,5
!    do iy=5,5
!    do iz=4,4
!    do ix=5,6
!    do iy=5,6
!    do iz=5,6
!if (ix.eq.5.and.iy.eq.5.and.iz.eq.5) cycle
!if (ix.eq.5.and.iy.eq.6.and.iz.eq.6) cycle
!if (ix.eq.6.and.iy.eq.5.and.iz.eq.6) cycle
!if (ix.eq.6.and.iy.eq.6.and.iz.eq.5) cycle
      ikk=(/ix,iy,iz/)
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)  ! integer index of k-q
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
!      write(6,'(4(3i3,3x))') ikk,jkk,iks,igg
!      write(6,'(3(3f10.6,3x))') xck,xckq
!      write(6,'(3(3f10.6,3x))') qq,qq-igglf1
!write(6,'(a,3f8.4)') "k   = ",xck
!write(6,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(6,*) "symmetry operation on k"
!write(6,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(6,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!write(44,*) "band = ",iunocc
!write(44,'(a,3f8.4)') "k   = ",xck
!write(44,'(a,3f8.4)') "k-q = ",xckq
!write(44,'(a,3f8.4)') "k   in irreducable Brillouin zone = ",kpt(1:3,ikpt)
!write(44,'(a,3f8.4)') "k-q in irreducable Brillouin zone = ",kpt(1:3,ikptq)
!write(44,'(a,3i3,3x)') "integer index of k   = ",ikk
!write(44,'(a,3i3,3x)') "integer index of k-q = ",jkk
!write(44,'(a,3i3,3x)') "integer index of q   = ",iks
!write(44,*) "symmetry operation on k"
!write(44,*) "symmetry index ",isym,mod(isym-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isym-1,nsym)+1),syminv(1:3,ii,mod(isym-1,nsym)+1)
!enddo
!write(44,*)
!write(44,*) "symmetry operation on k-q"
!write(44,*) "symmetry index ",isymq,mod(isymq-1,nsym)+1
!do ii=1,3
!write(44,'(3i3,5x,3i3)') symrel(1:3,ii,mod(isymq-1,nsym)+1),syminv(1:3,ii,mod(isymq-1,nsym)+1)
!enddo
!write(44,*)
!write(44,'(a,3i3)') "k   folding vector: ",igsymk
!write(44,'(a,3i3)') "k-q folding vector: ",igsymq
!write(44,'(a,3i3)') "igg               : ",igg
!write(44,'(a,3i3)') "igglf             : ",igglf1
      omega(ikk(1),ikk(2),ikk(3))=enrgy(indxkbnd(ikpt)+iunocc) &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      if (igglf1(1).eq.0.and.igglf1(2).eq.0.and.igglf1(3).eq.0) then
!      if (.true.) then
        call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         jmatel)
!write(6,'(3("(",f8.4,",",f8.4,")"))') jmatel
        if (ipaw.ne.0) then
          call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&                         pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
          jmatel(:) = jmatel(:) + vmatel(:)
        endif
!write(6,'(3("(",f8.4,",",f8.4,")"))') vmatel
!write(6,'(3("(",f8.4,",",f8.4,")"))') jmatel
        do ii=1,3
!          jmatel(ii) = jmatel(ii)/omega(ikk(1),ikk(2),ikk(3))
          jmatel(ii) = jmatel(ii)/(eigen(indxkbnd(ikpt)+iunocc) &
&                                 -eigen(indxkbnd(ikptq)+iocc))
        enddo
      else
        call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatel)
        if (ipaw.ne.0) then
          call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
&                        pwmatel1,tpwmatel1,nbcore,ncband,amatel)
          cmatel = cmatel + amatel
        endif
        jmatel(:) = cmatel*qq(:)/qq2
      endif
!write(6,*) jkk
!write(6,*) ihlf(ikpt),ihlf(ikptq)
      if (igglf1(1).eq.igglf2(1).and.igglf1(2).eq.igglf2(2).and. &
&         igglf1(3).eq.igglf2(3)) then
        do ii=1,3
        do jj=1,3
          tmat2(ikk(1),ikk(2),ikk(3),ii,jj)=dble(jmatel(ii)*conjg(jmatel(jj)))
if (tmat2(ikk(1),ikk(2),ikk(3),ii,jj).ne.tmat2(ikk(1),ikk(2),ikk(3),ii,jj)) then
  write(6,*) "ii,jj: ",ii,jj
  write(6,*) "ix,iy,iz: ",ikk
  write(6,*) "tmat2: ",tmat2(ikk(1),ikk(2),ikk(3),ii,jj)
  write(6,*) "jmatel: ",jmatel(ii)
  write(6,*) "jmatellf: ",jmatellf(jj)
  stop
endif
        enddo
        enddo
      else
        if (igglf2(1).eq.0.and.igglf2(2).eq.0.and.igglf2(3).eq.0) then
          call mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&           ncg,nkpt,npwt,igmx,igmn,igndx, &
&           isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&           lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&           cg,indxkcg,indxkpw,npwarr,kg, &
&           jmatellf)
          if (ipaw.ne.0) then
            call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&                           pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
            jmatellf(:) = jmatellf(:) + vmatel(:)
          endif
          do ii=1,3
!            jmatellf(ii) = jmatellf(ii)/omega(ikk(1),ikk(2),ikk(3))
            jmatellf(ii) = jmatellf(ii)/(eigen(indxkbnd(ikpt)+iunocc) &
&                                       -eigen(indxkbnd(ikptq)+iocc))
          enddo
        else
          call mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf2,igg, &
&           ncg,nkpt,npwt,igmx,igmn,igndx, &
&           isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&           lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&           cg,indxkcg,indxkpw,npwarr,kg, &
&           cmatel)
          if (ipaw.ne.0) then
            call mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qp,igg0, &
&                          pwmatel2,tpwmatel2,nbcore,ncband,amatel)
            cmatel = cmatel + amatel
          endif
          jmatellf(:) = cmatel*qp(:)/qp2
        endif
        do ii=1,3
        do jj=1,3
          tmat2(ikk(1),ikk(2),ikk(3),ii,jj)=jmatel(ii)*conjg(jmatellf(jj))
if (tmat2(ikk(1),ikk(2),ikk(3),ii,jj).ne.tmat2(ikk(1),ikk(2),ikk(3),ii,jj)) then
  write(6,*) "ii,jj: ",ii,jj
  write(6,*) "ix,iy,iz: ",ikk
  write(6,*) "tmat2: ",tmat2(ikk(1),ikk(2),ikk(3),ii,jj)
  write(6,*) "jmatel: ",jmatel(ii)
  write(6,*) "jmatellf: ",jmatellf(jj)
  stop
endif
        enddo
        enddo
      endif
      if (omega(ikk(1),ikk(2),ikk(3)).lt.gap) gap=omega(ikk(1),ikk(2),ikk(3))
!!do ii=1,3
!!  dummy=sqrt(blat(ii,1)**2+blat(ii,2)**2+blat(ii,3)**2)  
!!  do jj=1,3
!!    if (ii.eq.jj) then 
!!      matdummy1(ii,jj)=1.d0/dummy
!!    else 
!!      matdummy1(ii,jj)=0.d0
!!    endif
!!  enddo
!!enddo
!!! T * M
!!do ii=1,3
!!do jj=1,3
!!  cmatdummy1(ii,jj)=(0.d0,0d0)
!!  do kk=1,3
!!    cmatdummy1(ii,jj)=cmatdummy1(ii,jj)+tmat2(ikk(1),ikk(2),ikk(3),ii,kk)*bmet(kk,jj)
!!  enddo
!!enddo
!!enddo
!!! M * T * M
!!do ii=1,3
!!do jj=1,3
!!  cmatdummy2(ii,jj)=(0.d0,0d0)
!!  do kk=1,3
!!    cmatdummy2(ii,jj)=cmatdummy2(ii,jj)+bmet(ii,kk)*cmatdummy1(kk,jj)
!!  enddo
!!enddo
!!enddo
!!do ii=1,3
!!do jj=1,3
!!  cmatdummy3(ii,jj)=cmatdummy3(ii,jj)+matdummy1(ii,ii)*cmatdummy2(ii,jj)*matdummy1(jj,jj)
!!enddo
!!enddo
!! M * S^t
!do ii=1,3
!do jj=1,3
!  matdummy1(ii,jj)=0.d0
!  do kk=1,3
!    matdummy1(ii,jj)=matdummy1(ii,jj)+bmet(ii,kk)*symrel(jj,kk,isym)
!  enddo
!enddo
!enddo
!! S * M * S^t
!do ii=1,3
!do jj=1,3
!  matdummy2(ii,jj)=0.d0
!  do kk=1,3
!    matdummy2(ii,jj)=matdummy2(ii,jj)+symrel(ii,kk,isym)*matdummy1(kk,jj)
!  enddo
!enddo
!enddo
!do ii=1,3
!!  write(6,'(3f10.6)') (matdummy1(ii,jj),jj=1,3)
!!  write(6,'(3f10.6)') (matdummy2(ii,jj),jj=1,3)
!!  write(6,'(3f10.6)') (dble(cmatdummy2(ii,jj)*100.d0),jj=1,3)
!!  write(6,'(3f10.6)') (dble(cmatdummy3(ii,jj)*100.d0),jj=1,3)
!!  write(6,'(3f10.6)') (dble(tmat2(ix,iy,iz,ii,jj)*100.d0),jj=1,3)
!!  write(6,'(3f10.6)') (bmet(ii,jj),jj=1,3)
!  write(6,'(3f10.6)') (blat(ii,jj),jj=1,3)
!enddo
!write(46,'(2i3,3x,3i3,3x,3("(",f6.2,",",f6.2,")"))') iocc,iunocc,ikk(1:3), &
!& tmat2(ikk(1),ikk(2),ikk(3),1,1), &
!& tmat2(ikk(1),ikk(2),ikk(3),2,2), &
!& tmat2(ikk(1),ikk(2),ikk(3),3,3)
!write(29,'(3i5,3x,f9.6,4x,2f9.6)') ix,iy,iz,omega(ix,iy,iz),tmat2(ix,iy,iz,1,1)
!write(29,'(3i5,3x,2(3("(",f9.5,",",f9.5,")"),2x))') ix,iy,iz,jmatel,jmatellf
    enddo
    enddo
    enddo

!stop

!do ii=1,3
!do jj=1,3
!  iiq = ii + 3*(jj-1)
!  do ix=1,ngkpt(1)
!    write(70+iiq,*) ix
!    do iy=1,ngkpt(1)
!      write(70+iiq,'(10f8.4)') (dble(tmat2(ix,iy,iz,ii,jj)*100.d0),iz=1,ngkpt(3))
!    enddo
!    write(70+iiq,*)
!  enddo
!enddo
!enddo

    if (itetrahedron.ne.0) then
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=1,1
!      do iy=1,1
!      do iz=1,1
        ikk=(/ix,iy,iz/)
        call fhilo(omega,ikk,ngkpt,whi,wlo)
        iwh=min(nint(whi*dble(nwpt)/wmax),nwpt)
        iwl=nint(wlo*dble(nwpt)/wmax)
!write(27,'(5i5)') ix,iy,iz,iwl,iwh
        if (iwl.gt.nwpt) cycle
        call fval(omega,ikk,ikkp,ngkpt,enval)
        do ii=1,3
        do jj=1,3
!        do ii=1,1
!        do jj=1,1
          call fpol(tmat2(:,:,:,ii,jj),ngkpt,ikk,ikkp,cmat2v)
          do iw=iwl,iwh
!          do iw=46,46
            ww=iw*dw
!            call cubeint(enval,ww,cmat2v,tpolincr)
!            dpol=tpolincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol)
            call vcubeint(enval,ww-dw/2.d0,dw,cmat2v,tpolincr)
            dpol=-8.d0*pi*tpolincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*dw)
            polci(iw,ii,jj)=polci(iw,ii,jj)+dpol
!if (polci(iw,ii,jj).ne.polci(iw,ii,jj)) then
if (dpol.ne.dpol) then
  write(6,*) "iw: ",iw
  write(6,*) "ii,jj: ",ii,jj
  write(6,*) "ix,iy,iz: ",ikk
  write(6,*) "dpol: ",dpol
  write(6,*) "tpolincr: ",tpolincr
  write(6,*) "enval: ",enval
  write(6,*) "cmat2v: ",cmat2v
  stop
endif
!if ((iw.eq.33.or.iw.eq.32).and.ii.eq.1.and.jj.eq.1) then
!  if (ix.eq.4.and.iy.eq.4.and.iz.eq.4) then
!    write(27,'(5i5,4e12.4)') ix,iy,iz,iw,nwpt,dpol*8.d0*pi
!    write(27,'(f9.6)') ww
!    write(27,'(2f9.6,4x,2f9.6)') enval(1),enval(2),enval(5),enval(6)
!    write(27,'(2f9.6,4x,2f9.6)') enval(3),enval(4),enval(7),enval(8)
!    write(27,*)
!    write(27,'(2f9.6,4x,2f9.6)') dble(cmat2v(1)),dble(cmat2v(2)),dble(cmat2v(5)),dble(cmat2v(6))
!    write(27,'(2f9.6,4x,2f9.6)') dble(cmat2v(3)),dble(cmat2v(4)),dble(cmat2v(7)),dble(cmat2v(8))
!    write(27,*)
!    write(27,*)
!  endif
!endif
!            write(77,'(f12.6,2e20.6)') ww,dpol
          enddo
        enddo
        enddo
!write(27,*)
!        stvec(ikk(3))=dble(iwh)
!        stvec(ikk(3))=dble(dpol*2.d0*vq)
      enddo
!      write(68,'(10f8.3)') stvec(1:ngkpt(3))
      enddo
!      write(68,*) ix+1,'------------------------------------------------------'
!      write(68,*)
      enddo
    else
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
        iw=int(omega(ix,iy,iz)/dw)
        if (iw.le.nwpt) then
          fw=omega(ix,iy,iz)/dw-iw
          do ii=1,3
          do jj=1,3
            dpol=8.0*pi*tmat2(ix,iy,iz,ii,jj)/(pi*vol)
!write(6,'(5i5,4e12.4)') ix,iy,iz,iw,nwpt,dpol
            polci(iw,ii,jj)=polci(iw,ii,jj)-fw*dpol
            if (iw.lt.nwpt) polci(iw+1,ii,jj)=polci(iw+1,ii,jj)-(1.d0-fw)*dpol
          enddo
          enddo
        endif
      enddo
      enddo
      enddo
    endif
  enddo
enddo

!do iw=1,nwpt
!do ii=1,3
!do jj=1,3
!  iiq = ii + 3*(jj-1)
!  write(80+iiq,'(i4,3x,2e14.6)') iw,polci(iw,ii,jj)
!enddo
!enddo
!enddo

! polci gives "imaginary" part of polarization - i.e., the anti-hermetian part when dealing with non-diagonal components.
  do iw=1,nwpt
    ww1=dble(iw)*dw
!write(28,'(i5,f9.2,9e12.4)') iw,(ww1*27.21138),((dble(polci(iw,ii,jj)),jj=1,3),ii=1,3)
!write(28,'(i5,f9.2,9e12.4)') iw,ww1,((dble(polci(iw,ii,jj)),jj=1,3),ii=1,3)
    do jw=1,nwpt
      if (iw.eq.jw) cycle
      ww2=dble(jw)*dw
      do ii=1,3
      do jj=1,3
        iiq = ii + 3*(jj-1)
        polat(iw,iiq)=polat(iw,iiq)+2.d0*polci(jw,ii,jj)*dw*(ww2/(ww2**2-ww1**2))
      enddo
      enddo
    enddo
    do ii=1,3
    do jj=1,3
      iiq = ii + 3*(jj-1)
      polbt(iw,iiq)=-(0.d0,1.d0)*pi*polci(iw,ii,jj)
    enddo
    enddo
!write(28,'(i5,f9.2,9e12.4)') iw,(ww1*27.21138),(dimag(polbt(iw,iiq)),iiq=1,9)
!write(28,'(i5,f9.2,9e12.4)') iw,(ww1*27.21138),(dble(polat(iw,iiq)),iiq=1,9)
  enddo

return
end subroutine mkpolj1

!*************************************************************************

! for approximating dielectric function by transitions of valence electrons to plane wave conduction states

subroutine mkpwpol(iqpt,qq,igglf1,iks,vq,vol,pi,bmet,blat, &
& pwse_mass, &
& test_bands_pol, &
& nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
& nsym,nkpt, &
& ncg,bantot,npwt,nlmn, &
& npwwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
& xred, &
& projwf, &
& kg, &
& enrgy, eigen, &
& cg, &
& indxkpw, &
& indxkbnd, &
& indxkcg, &
& npwarr, &
& kpt, &
& nsymk,symk,symrel,syminv, &
& ihlf,lvtrans, &
& ntypepaw, &
& pwmatel1,tpwmatel1, &
& igndx,ikndx,isymndx,isymg, &
& nband,pwpola,pwpolb)
use pwx
implicit none
integer :: iqpt,bantot,ncg,nkpt,ipw1,ipw2,ipw,nlmn,natom
double precision :: qq(3),dq(3),vq,vol,pi,wmax,dw,bmet(3,3),blat(3,3),pwse_mass
integer :: npwwpt,nbcore,nbocc,ncband,ngkpt(3),npwt,npw,iks(3)
integer :: nsym,igmx(3),igmn(3),igcut,nsppol,ipaw,itetrahedron
integer :: igglf1(3)
integer :: test_bands_pol(4)
double precision :: shiftk(3)
double precision :: xred(3,natom)
integer :: kg(3,npwt)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double precision :: enrgy(bantot), eigen(bantot)
! state energies
double complex :: cg(ncg)
! wave functions calculated by abinit
integer :: indxkpw(nkpt),indxkbnd(nkpt)
integer :: indxkcg(nkpt),npwarr(nkpt)
double precision :: kpt(3,nkpt)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: symrel(3,3,nsym),syminv(3,3,nsym)
! symmetry operations
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
! symmetry indexes
integer :: isymg(3,nkpt,nsym),igsymk(3),igsymq(3)
! translation of symmetry transformed k-point to recover tabulated k-point
integer :: nband(nkpt*nsppol)
! number of bands calculated by abinit
double complex :: pwpola(0:npwwpt),pwpolb(0:npwwpt),tpolincr,polci(0:npwwpt),dpol
! pwpola - hermitian part of vq * polarization matrix element for igglf1
! pwpolb - antihermitian part of vq * polarization matrix element 
double complex :: vmat2(ngkpt(1),ngkpt(2),ngkpt(3)),cmatel,cmatellf,vmat2v(8)
double complex :: jmatel(3)
! cmatel - Complex density MATrix ELement
! jmatel - current matrix element
! vmat2 - "square" of matrix elements times vq, cmatel*conjuagte(cmatel)
double complex :: amatel,amatellf
! amatel, amatellf - contribution to matrix elements from PAW
double complex :: tmat2(3,3),tmat2cart(3,3)
! tmat2 - tensor outer product of jmatel
! tmat2cart - as above, but in cartesian space
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
! omega - transition energy from conduction to valence state
double precision :: whi,wlo,ww,ww1,ww2,enval(8)
! energy variables
double precision :: volelmnt
integer :: iwh,iwl
double precision :: xck(3),xckq(3)
double precision :: xkkv(3),xkkv2
double precision :: testksym(3)
double precision :: stvec(ngkpt(3)),stvec2(ngkpt(3))
integer :: ivgndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3),ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,ll,iw,jw,iocc,iunocc,ikk(3),ikkp(3),jkk(3),jka(3),isym,isymq
integer :: ig,jg,igg(3),jgg(3),jgs(3),ikpt,ikptq,ix,iy,iz,igx,igy,igz,igg0(3)
double precision :: fw
integer :: ntypepaw
double complex :: pwmatel1(ntypepaw,nlmn,nlmn),tpwmatel1(ntypepaw,nlmn,nlmn)
double precision :: omegap2, val_edens
double precision :: fsum,fsum_correction
double precision :: gsum,gsum_correction
! DEBUG
! DEBUG

volelmnt=vol*ngkpt(1)*ngkpt(2)*ngkpt(3)
igg0=(/0,0,0/)
dw=wmax/dble(npwwpt)
do iw=0,npwwpt
  polci(iw)=(0.d0,0.d0)
  pwpola(iw)=(0.d0,0.d0)
  pwpolb(iw)=(0.d0,0.d0)
enddo

val_edens = dble(test_bands_pol(2) - test_bands_pol(1))/vol
omegap2 = 4.d0*pi*val_edens

!write(6,'(a)') 'occupied band, unoccupied band'
do iocc=test_bands_pol(1),test_bands_pol(2)
!do iocc=1,1
  do iunocc=1,npwx
!  do iunocc=8,8
!if (iqpt.ge.30) write(*,*) iocc,iunocc
    do ix=1,ngkpt(1)
    do iy=1,ngkpt(2)
    do iz=1,ngkpt(3)
!    do ix=1,1
!    do iy=1,1
!    do iz=5,5
      ikk=(/ix,iy,iz/)
!if (iocc .gt. 1 .and. iqpt.ge.30 ) write(*,*) ikk
      ikpt=ikndx(ikk(1),ikk(2),ikk(3))
      isym=isymndx(ikk(1),ikk(2),ikk(3))
      igsymk=isymg(:,ikpt,isym)
      do ii=1,3
        jka(ii)=ikk(ii)-iks(ii)  ! integer index of k-q
        jkk(ii)=mod(jka(ii)-1,ngkpt(ii))+1  ! shifted to 1st Brillouin zone
        if (jkk(ii).le.0) jkk(ii)=jkk(ii)+ngkpt(ii)
        xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xck(ii)=mod(xck(ii),1.d0)
        if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
        xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
        xck(ii)=xck(ii)-0.5d0
        xckq(ii)=dble(jkk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
        xckq(ii)=mod(xckq(ii),1.d0)
        if (xckq(ii).lt.0.d0) xckq(ii)=xckq(ii)+1
        xckq(ii)=1.d0-mod(1.d0-xckq(ii),1.d0)
        xckq(ii)=xckq(ii)-0.5d0
        igg(ii)=nint(xckq(ii)-xck(ii)+qq(ii))-igglf1(ii)
      enddo
      ikptq=ikndx(jkk(1),jkk(2),jkk(3))
      isymq=isymndx(jkk(1),jkk(2),jkk(3))
      igsymq=isymg(:,ikptq,isymq)
      xkkv(:) = kpt(:,ikpt)+ipwx(:,iunocc)
      xkkv2 = 0.d0
      do ii=1,3
      do jj=1,3
        xkkv2=xkkv2+xkkv(ii)*bmet(ii,jj)*xkkv(jj)
      enddo
      enddo
      omega(ikk(1),ikk(2),ikk(3))=0.5*xkkv2 &
&                                -enrgy(indxkbnd(ikptq)+iocc)
      if (igglf1(1).eq.0.and.igglf1(2).eq.0.and.igglf1(3).eq.0.and. &
&         iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
        call mkvpwmatelJ1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         jmatel)
!        if (ipaw.ne.0) then
!          call mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
!&                         pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
!          jmatel(:) = jmatel(:) + vmatel(:)
!        endif
        do ii=1,3
        do jj=1,3
          tmat2(ii,jj) = jmatel(ii)*jmatel(jj)*4.d0*pi
        enddo
        enddo
        do ii=1,3
        do jj=1,3
          tmat2cart(ii,jj) = (0.d0,0.d0)
          do kk=1,3
          do ll=1,3
            tmat2cart(ii,jj) = tmat2cart(ii,jj) + blat(kk,ii)*tmat2(kk,ll)*blat(ll,jj)
          enddo
          enddo
        enddo
        enddo
        vmat2(ikk(1),ikk(2),ikk(3))=(tmat2cart(1,1)+tmat2cart(2,2)+tmat2cart(3,3))/3.d0
      else 
        call mkvpwmatelX1(iocc,iunocc,ikpt,ikptq,igglf1,igg, &
&         ncg,nkpt,npwt,igmx,igmn,igndx, &
&         isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
&         lvtrans(1:3,ikk(1),ikk(2),ikk(3)), &
&         cg,indxkcg,indxkpw,npwarr,kg, &
&         cmatel)
!        if (ipaw.ne.0) then
!          call mkvpwPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg0, &
!& pwmatel1,tpwmatel1,nbcore,ncband,amatel)
!        else
          amatel=(0.d0,0.d0)
!        endif
        vmat2(ikk(1),ikk(2),ikk(3))=dble((cmatel+amatel)*conjg(cmatel+amatel))*vq
      endif
    enddo
    enddo
    enddo

    if (itetrahedron.ne.0) then
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
!      do ix=10,10
!      do iy=10,10
!      do iz=10,10
        ikk=(/ix,iy,iz/)
        call fhilo(omega,ikk,ngkpt,whi,wlo)
        iwh=min(nint(whi*dble(npwwpt)/wmax),npwwpt)
        iwl=max(nint(wlo*dble(npwwpt)/wmax),0)
        if (iwl.gt.npwwpt) cycle
        call fval(omega,ikk,ikkp,ngkpt,enval)
        call fpol(vmat2,ngkpt,ikk,ikkp,vmat2v)
!if (iocc .gt. 1 .and. iqpt.eq.30) write(6,*) ikk,iwl,iwh
        do iw=iwl,iwh
!        do iw=166,166
!if (iocc .gt. 1 .and. iqpt.eq.30 .and. ix.eq.1 .and. iy.eq.1 .and. iz.eq.1) write(6,*) iw
          ww=iw*dw
          call vcubeint(enval,ww-dw/2.d0,dw,vmat2v,tpolincr)
          dpol=tpolincr/(dble(ngkpt(1)*ngkpt(2)*ngkpt(3))*vol*dw)
          polci(iw)=polci(iw)-dpol*2.d0
        enddo
!if (iocc .gt. 1 .and. iqpt.eq.30) write(6,*) ikk
      enddo
      enddo
      enddo
    else
      do ix=1,ngkpt(1)
      do iy=1,ngkpt(2)
      do iz=1,ngkpt(3)
        iw=int(omega(ix,iy,iz)/dw)
        if (iw.le.npwwpt) then
          fw=omega(ix,iy,iz)/dw-iw
          dpol=2.d0*vmat2(ix,iy,iz)/(pi*vol)
          polci(iw)=polci(iw)-fw*dpol
          if (iw.lt.npwwpt) polci(iw+1)=polci(iw+1)-(1.d0-fw)*dpol
        endif
      enddo
      enddo
      enddo
    endif
  enddo
enddo

! polci gives "imaginary" part of polarization - i.e., the anti-hermetian part when dealing with non-diagonal components.
  do iw=0,npwwpt
    ww1=dble(iw)*dw
    do jw=0,npwwpt
      if (iw.eq.jw) cycle
      ww2=dble(jw)*dw
      pwpola(iw)=pwpola(iw)+2.d0*polci(jw)*dw*(ww2/(ww2**2-ww1**2))
    enddo
    pwpolb(iw)=-(0.d0,1.d0)*pi*polci(iw)
  enddo

  fsum = 0.d0
  gsum = 0.d0
  do iw=0,npwwpt
    ww1=dble(iw)*dw
    fsum = fsum + ww1*dw*dimag(1.d0/(1.d0-pwpola(iw)-pwpolb(iw)))
    gsum = gsum + ww1*dw*dimag(1.d0-pwpola(iw)-pwpolb(iw))
  enddo
  fsum_correction =  fsum*2.d0/(pi*omegap2)
  gsum_correction = -gsum*2.d0/(pi*omegap2)
write(6,*) fsum_correction, gsum_correction

! debug
!  do iw=1,npwwpt
!    write(27,'(i5,3x,4e12.4)') iw,polci(iw)
!!    write(27,'(i5,3x,4e12.4)') iw,pwpola(iw),pwpolb(iw)
!  enddo
! debug

return
end subroutine mkpwpol

!*************************************************************************

subroutine fhilo(omega,ikk,ngkpt,whi,wlo)
implicit none
integer :: ikk(3),ikkp(3),ngkpt(3),ii
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: whi,wlo

do ii=1,3
  ikkp(ii)=ikk(ii)+1
  if (ikkp(ii).gt.ngkpt(ii)) ikkp(ii)=ikkp(ii)-ngkpt(ii)
enddo

whi=max(omega(ikk (1),ikk (2),ikk (3)), &
&       omega(ikkp(1),ikk (2),ikk (3)), &
&       omega(ikk (1),ikkp(2),ikk (3)), &
&       omega(ikkp(1),ikkp(2),ikk (3)), &
&       omega(ikk (1),ikk (2),ikkp(3)), &
&       omega(ikkp(1),ikk (2),ikkp(3)), &
&       omega(ikk (1),ikkp(2),ikkp(3)), &
&       omega(ikkp(1),ikkp(2),ikkp(3)))

wlo=min(omega(ikk (1),ikk (2),ikk (3)), &
&       omega(ikkp(1),ikk (2),ikk (3)), &
&       omega(ikk (1),ikkp(2),ikk (3)), &
&       omega(ikkp(1),ikkp(2),ikk (3)), &
&       omega(ikk (1),ikk (2),ikkp(3)), &
&       omega(ikkp(1),ikk (2),ikkp(3)), &
&       omega(ikk (1),ikkp(2),ikkp(3)), &
&       omega(ikkp(1),ikkp(2),ikkp(3)))

return
end subroutine fhilo

!*************************************************************************

subroutine mkmatelX1(iocc,iunocc,ikpt,ikptq,igglf,igg, &
& ncg,nkpt,npwt,igmx,igmn,igndx, &
& isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
& lvt, &
& cg,indxkcg,indxkpw,npwarr,kg, &
& cmatel)
implicit none
integer :: iocc,iunocc,ikpt,ikptq,ncg,nkpt,npwt,igg(3)
integer :: isym,isymq,nsym,symrel(3,3,nsym),syminv(3,3,nsym)!,symproduct(3,3)
integer :: lvt(3)
integer :: igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: cmatel,xcg,xcgq
integer :: iisym,isign,iisymq,isignq,ihlf(nkpt)
double precision :: kpt(3,nkpt)
integer :: ig,igs,ii,jj,kk,igq,icg,icgq,igo
integer :: iggi(3),iggs(3),iggq(3),iggqs(3),igglf(3)
integer :: igsymk(3),igsymq(3),ishift(3)
! iggi  - initial lattice vector of k point
! iggq  - lattice vector of k-q point
! iggs  - lattice vector of kpt after symmetry transform to bring k point to calculated k point set
! iggqs - similar to iggs, but for k-q point

!write(6,*) ikpt,ikptq
!write(44,'(a,3i3)') "k   folding vector: ",igsymk
!write(44,'(a,3i3)') "k-q folding vector: ",igsymq
!write(44,'(a,3i3)') "igg               : ",igg
!write(44,'(a,3i3)') "igglf             : ",igglf
cmatel=(0.d0,0.d0)
ishift=igsymq-igg-igglf -igsymk
l1: do ig=1,npwarr(ikpt)
!write(6,*) ig
  icg=indxkcg(ikpt)+(iunocc-1)*npwarr(ikpt)+ig
  xcg=cg(icg)
  iggi(1:3)=kg(1:3,indxkpw(ikpt)+ig)
  iggs = ishift
  do ii=1,3
    do jj=1,3
      iggs(ii) = iggs(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*iggi(jj)
    enddo
  enddo
  iggq = (/0,0,0/)
  do ii=1,3
    do jj=1,3
      iggq(ii) = iggq(ii)+ syminv(jj,ii,mod(isymq-1,nsym)+1)*iggs(jj)
    enddo
    if (iggq(ii).lt.igmn(ii).or.iggq(ii).gt.igmx(ii)) cycle l1
  enddo
  igq=igndx(ikptq,iggq(1),iggq(2),iggq(3))
  if (igq.eq.0) cycle
  icgq=indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+igq
  xcgq=cg(icgq)
  cmatel=cmatel+xcg*conjg(xcgq)
!write(44,'(i4,2x,4(3i3,2x),2("(",f6.3,",",f6.3,")"))') ig, iggi, iggs-igsymq+igg+igglf, iggs-igsymq, iggq, xcg, xcgq  !, xcg*conjg(xcgq)
!write(44,'(i4,2x,4(3i3,2x),f8.5,2x,f8.5,2x,2i3)') ig, iggi, iggs-igsymq+igg+igglf, iggs-igsymq, iggq, dble(xcg), dble(xcgq) !, isym, isymq
!write(44,'(i4,3x,2(3i3,3x),2("(",f8.5,",",f8.5,")",2x))') ig, iggs, iggq, xcg, xcgq  !, xcg*conjg(xcgq)
!write(44,'(i4,3x,2(3i3,3x),2("(",f8.5,",",f8.5,")",2x))') ig, iggi, iggq, xcg, xcgq  !, xcg*conjg(xcgq)
!  write(6,'(3i3,3x,2f10.6,3x,3i3,3x,2f10.6)') &
!&               kg(1:3,indxkpw(ikpt)+ig),xcg, &
!&               iggq(1:3),xcgq
!  write(6,1000) ig,iggi(1:3),iggq(1:3),xcg,xcgq,cmatel
!1000 format(i4,4x,'a',3i3,4x,'b',3i3,3(1x,2f7.3))
enddo l1

!cmatel=(0.d0,0.d0)
!ishift=igsymq-igg-igglf
!l1: do ig=1,npwarr(ikpt)
!!write(6,*) ig
!  icg=indxkcg(ikpt)+(iunocc-1)*npwarr(ikpt)+ig
!  xcg=cg(icg)
!  iggi(1:3)=kg(1:3,indxkpw(ikpt)+ig)
!!write(44,'(i4,3x,3i3,3x,"(",f10.6,",",f10.6,")")') ig, iggi, xcg
!  iggs = -igsymk
!  do ii=1,3
!    do jj=1,3
!      iggs(ii) = iggs(ii)+ syminv(ii,jj,mod(isym-1,nsym)+1)*iggi(jj)
!    enddo
!  enddo
!  do ii=1,3
!!    iggq(ii)=iggs(ii)-igg(ii)-igglf(ii)+igsymq(ii)
!    iggq(ii)=iggs(ii)+ishift(ii)
!  enddo
!!write(6,'(3i4)') iggq
!  iggqs = (/0,0,0/)
!  do ii=1,3
!    do jj=1,3
!      iggqs(ii) = iggqs(ii)+ symrel(ii,jj,mod(isymq-1,nsym)+1)*iggq(jj)
!!      write(6,'(2i2,3x,3i4)') ii,jj, symproduct(ii,jj), iggq(jj), iggqs(ii)
!    enddo
!!write(6,'(3i4)') iggqs
!!write(6,*)
!    if (iggqs(ii).lt.igmn(ii).or.iggqs(ii).gt.igmx(ii)) cycle l1
!  enddo
!  igq=igndx(ikptq,iggqs(1),iggqs(2),iggqs(3))
!  if (igq.eq.0) cycle
!  icgq=indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+igq
!  xcgq=cg(icgq)
!  cmatel=cmatel+xcg*conjg(xcgq)
!write(44,'(i4,2x,4(3i3,2x),2("(",f6.3,",",f6.3,")"))') ig, iggi, iggs, iggq, iggqs, xcg, xcgq  !, xcg*conjg(xcgq)
!!write(44,'(i4,3x,2(3i3,3x),2("(",f8.5,",",f8.5,")",2x))') ig, iggs, iggq, xcg, xcgq  !, xcg*conjg(xcgq)
!!  write(6,'(3i3,3x,2f10.6,3x,3i3,3x,2f10.6)') &
!!&               kg(1:3,indxkpw(ikpt)+ig),xcg, &
!!&               iggq(1:3),xcgq
!!  write(6,1000) ig,iggi(1:3),iggq(1:3),xcg,xcgq,cmatel
!!1000 format(i4,4x,'a',3i3,4x,'b',3i3,3(1x,2f7.3))
!enddo l1

return
end subroutine mkmatelX1

!*************************************************************************

subroutine mkmatelJ1(iocc,iunocc,ikpt,ikptq,igglf,igg, &
& ncg,nkpt,npwt,igmx,igmn,igndx, &
& isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
& lvt, &
& cg,indxkcg,indxkpw,npwarr,kg, &
& jmatel)
implicit none
integer :: iocc,iunocc,ikpt,ikptq,ncg,nkpt,npwt,igg(3)
integer :: isym,isymq,nsym,symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvt(3)
integer :: igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: jmatel(3),cmatel,xcg,xcgq
! jmatel expressed as components of linear combination of reciprocal lattice vectors
integer :: iisym,isign,iisymq,isignq,ihlf(nkpt)
double precision :: kpt(3,nkpt),kvec(3),kvecq(3)
double precision :: dvec(3),dvecl(3)
integer :: ig,igr,ii,jj,igq,icg,icgq,igo
integer :: iggi(3),iggs(3),iggq(3),iggqs(3),igglf(3)
integer :: igsymk(3),igsymq(3),ishift(3)
! iggi  - initial lattice vector of k point
! iggq  - lattice vector of k-q point
! iggs  - lattice vector of kpt after symmetry transform to bring k point to calculated k point set
! iggqs - similar to iggs, but for k-q point

kvec = (/0.d0,0.d0,0.d0/)
kvecq = (/0.d0,0.d0,0.d0/)
do ii=1,3
  do jj=1,3
    kvec(ii) = kvec(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*kpt(jj,ikpt)
    kvecq(ii) = kvecq(ii)+ symrel(jj,ii,mod(isymq-1,nsym)+1)*kpt(jj,ikptq)
  enddo
enddo
!write(44,'(3f8.4)') kpt(:,ikpt)
!write(44,'(3f8.4)') kpt(:,ikptq)
!write(44,'(3f8.4)') kvec
!write(44,'(3f8.4)') kvecq

do ii=1,3
  jmatel(ii)=(0.d0,0.d0)
enddo
ishift=igsymq-igg-igglf -igsymk
l1: do ig=1,npwarr(ikpt)
  icg=indxkcg(ikpt)+(iunocc-1)*npwarr(ikpt)+ig
  xcg=cg(icg)
  iggi(1:3)=kg(1:3,indxkpw(ikpt)+ig)
  iggs = ishift
  do ii=1,3
    do jj=1,3
      iggs(ii) = iggs(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*iggi(jj)
    enddo
  enddo
  iggq = (/0,0,0/)
  do ii=1,3
    do jj=1,3
      iggq(ii) = iggq(ii)+ syminv(jj,ii,mod(isymq-1,nsym)+1)*iggs(jj)
    enddo
    if (iggq(ii).lt.igmn(ii).or.iggq(ii).gt.igmx(ii)) cycle l1
  enddo
  igq=igndx(ikptq,iggq(1),iggq(2),iggq(3))
  if (igq.eq.0) cycle
  icgq=indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+igq
  xcgq=cg(icgq)
  cmatel=xcg*conjg(xcgq);
  iggqs = iggs-igsymq
  iggs = iggs-igsymq+igg+igglf
  do ii=1,3
    dvecl(ii) = kvec(ii)+iggs(ii)+kvecq(ii)+iggqs(ii)
  enddo
  jmatel(1:3)=jmatel(1:3) + cmatel*dvecl(1:3)/2.d0
!write(46,'(i4,3x,4(3i3,3x),2("(",f8.5,",",f8.5,")",2x))') ig, iggi, iggs, iggq, iggqs !, xcg, xcgq  !, xcg*conjg(xcgq)
!write(46,'(i4,3x,2(3i3,3x),2("(",f8.5,",",f8.5,")",2x))') ig, iggi, iggq, xcg, xcgq  !, xcg*conjg(xcgq)
!write(46,'(i4,2x,4(3i3,2x),f8.5,2x,f8.5)') ig, iggi, iggs-igsymq+igg+igglf, iggs-igsymq, iggq, dble(xcg), dble(xcgq)  !, xcg*conjg(xcgq)
!write(46,'(i4,3x,2(3i3,3x),"(",f8.5,",",f8.5,")")') ig, iggs, iggq, xcg*conjg(xcgq)
!write(46,'(2(3i3,3x),f8.5,3f6.1,2x,3f8.5)') &
!& iggs, iggqs, dble(cmatel), dvecl, dble(jmatel)
!  write(6,'(3i3,3x,2f10.6,3x,3i3,3x,2f10.6)') &
!&               kg(1:3,indxkpw(ikpt)+ig),xcg, &
!&               iggq(1:3),xcgq
!  write(6,1000) ig,iggi(1:3),iggq(1:3),xcg,xcgq,cmatel
!1000 format(i4,4x,'a',3i3,4x,'b',3i3,3(1x,2f7.3))
enddo l1

return
end subroutine mkmatelJ1

!*************************************************************************

subroutine mkPAWmatelX1(pi,iocc,iunocc,ikk,jkk,qq,igg, &
&  pwmatel,tpwmatel,nbcore,ncband,amatel)
! PAW contribution to density matrix elements
use wfkvars
use geometry
use pawvars
implicit none
integer :: iocc,iunocc,ikk(3),jkk(3),igg(3),nbcore,ncband
double precision :: qq(3),pi
double complex :: pwmatel(ntypat,nlmnmax,nlmnmax),tpwmatel(ntypat,nlmnmax,nlmnmax)
double complex :: amatel,phase
integer :: iorb,jorb,iat,itpaw,ikv,jkv

ikv = ikk(1) + (ikk(2)-1)*ngkpt(1) + (ikk(3)-1)*ngkpt(2)*ngkpt(1)
jkv = jkk(1) + (jkk(2)-1)*ngkpt(1) + (jkk(3)-1)*ngkpt(2)*ngkpt(1)
amatel=(0.d0,0.d0)
do iat=1,natom
  itpaw=typat(iat)
  phase=exp((0.d0,2.d0)*pi*dot_product(-qq-igg,xred(:,iat)))
!write(6,'(3f10.5)') xred(:,iat)
!write(6,'(3f10.5)') -qq-igg
!write(6,'(3f10.5)') dot_product(-qq-igg,xred(:,iat))
!write(6,'(4f16.8)') phase
  do iorb=1,nlmn(itpaw)
  do jorb=1,nlmn(itpaw)
    amatel=amatel+projwf(iat,iorb,ikv,iunocc) &
&          *conjg(projwf(iat,jorb,jkv,iocc)) &
&          *(pwmatel(itpaw,iorb,jorb)-tpwmatel(itpaw,iorb,jorb))*phase
!write(6,'(2i3,4f16.8)') iorb,jorb,projwf(iat,iorb,ikv,iunocc),conjg(projwf(iat,jorb,jkv,iocc))
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb),tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb)-tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') projwf(iat,iorb,ikv,iunocc)*conjg(projwf(iat,jorb,jkv,iocc))*(pwmatel(iorb,jorb)-tpwmatel(iorb,jorb))*phase
!write(6,'(6x,4f16.8)') amatel
!write(6,*)
  enddo
  enddo
enddo

return
end subroutine mkPAWmatelX1

!*************************************************************************

subroutine mkPAWmatelJ1(pi,iocc,iunocc,ikk,jkk, &
&  pwjmatel,tpwjmatel,nbcore,ncband,vmatel)
! PAW contribution to current matrix elements
use wfkvars
use geometry
use pawvars
implicit none
integer :: iocc,iunocc,ikk(3),jkk(3),nbcore,ncband
double precision :: pi
double complex :: pwjmatel(3,ntypat,nlmnmax,nlmnmax),tpwjmatel(3,ntypat,nlmnmax,nlmnmax)
double complex :: vmatel(3),phase,factr
integer :: iorb,jorb,iat,itpaw,ikv,jkv,ii

ikv = ikk(1) + (ikk(2)-1)*ngkpt(1) + (ikk(3)-1)*ngkpt(2)*ngkpt(1)
jkv = jkk(1) + (jkk(2)-1)*ngkpt(1) + (jkk(3)-1)*ngkpt(2)*ngkpt(1)
do ii=1,3
  vmatel(ii)=(0.d0,0.d0)
enddo
do iat=1,natom
  itpaw=typat(iat)
  do iorb=1,nlmn(itpaw)
  do jorb=1,nlmn(itpaw)
    factr = projwf(iat,iorb,ikv,iunocc)*conjg(projwf(iat,jorb,jkv,iocc))
    do ii=1,3
      vmatel(ii)=vmatel(ii)+ factr &
&            *(pwjmatel(ii,itpaw,iorb,jorb)-tpwjmatel(ii,itpaw,iorb,jorb))
    enddo
!write(6,'(2i3,4f16.8)') iorb,jorb,projwf(iat,iorb,ikv,iunocc),conjg(projwf(iat,jorb,jkv,iocc))
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb),tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb)-tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') projwf(iat,iorb,ikv,iunocc)*conjg(projwf(iat,jorb,jkv,iocc))*(pwmatel(iorb,jorb)-tpwmatel(iorb,jorb))*phase
!write(6,'(6x,4f16.8)') amatel
!write(6,*)
  enddo
  enddo
enddo

return
end subroutine mkPAWmatelJ1

!*************************************************************************

subroutine mkvpwmatelX1(iocc,iunocc,ikpt,ikptq,igglf,igg, &
& ncg,nkpt,npwt,igmx,igmn,igndx, &
& isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
& lvt, &
& cg,indxkcg,indxkpw,npwarr,kg, &
& cmatel)
use pwx
implicit none
integer :: iocc,iunocc,ikpt,ikptq,ncg,nkpt,npwt,igg(3)
integer :: isym,isymq,nsym,symrel(3,3,nsym),syminv(3,3,nsym)!,symproduct(3,3)
integer :: lvt(3)
integer :: igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: cmatel,xcg,xcgq
integer :: iisym,isign,iisymq,isignq,ihlf(nkpt)
double precision :: kpt(3,nkpt)
integer :: ig,igs,ii,jj,kk,igq,icg,icgq,igo
integer :: iggi(3),iggs(3),iggq(3),iggqs(3),igglf(3)
integer :: igsymk(3),igsymq(3),ishift(3)
! iggi  - initial lattice vector of k point
! iggq  - lattice vector of k-q point
! iggs  - lattice vector of kpt after symmetry transform to bring k point to calculated k point set
! iggqs - similar to iggs, but for k-q point

!write(6,*) ikpt,ikptq
!write(44,'(a,3i3)') "k   folding vector: ",igsymk
!write(44,'(a,3i3)') "k-q folding vector: ",igsymq
!write(44,'(a,3i3)') "igg               : ",igg
!write(44,'(a,3i3)') "igglf             : ",igglf
cmatel=(0.d0,0.d0)
ishift=igsymq-igg-igglf -igsymk

  iggi(:)=ipwx(:,iunocc)
  iggs = ishift
  do ii=1,3
    do jj=1,3
      iggs(ii) = iggs(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*iggi(jj)
    enddo
  enddo
  iggq = (/0,0,0/)
  do ii=1,3
    do jj=1,3
      iggq(ii) = iggq(ii)+ syminv(jj,ii,mod(isymq-1,nsym)+1)*iggs(jj)
    enddo
    if (iggq(ii).lt.igmn(ii).or.iggq(ii).gt.igmx(ii)) return
  enddo
  igq=igndx(ikptq,iggq(1),iggq(2),iggq(3))
  if (igq.eq.0) return
  icgq=indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+igq
  xcgq=cg(icgq)
  cmatel=cmatel+conjg(xcgq)

return
end subroutine mkvpwmatelX1

!*************************************************************************

subroutine mkvpwmatelJ1(iocc,iunocc,ikpt,ikptq,igglf,igg, &
& ncg,nkpt,npwt,igmx,igmn,igndx, &
& isym,isymq,symrel,syminv,nsym,ihlf,igsymk,igsymq,kpt, &
& lvt, &
& cg,indxkcg,indxkpw,npwarr,kg, &
& jmatel)
use pwx
implicit none
integer :: iocc,iunocc,ikpt,ikptq,ncg,nkpt,npwt,igg(3)
integer :: isym,isymq,nsym,symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvt(3)
integer :: igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg)
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: jmatel(3),cmatel,xcg,xcgq
! jmatel expressed as components of linear combination of reciprocal lattice vectors
integer :: iisym,isign,iisymq,isignq,ihlf(nkpt)
double precision :: kpt(3,nkpt),kvec(3),kvecq(3)
double precision :: dvec(3),dvecl(3)
integer :: ig,igr,ii,jj,igq,icg,icgq,igo
integer :: iggi(3),iggs(3),iggq(3),iggqs(3),igglf(3)
integer :: igsymk(3),igsymq(3),ishift(3)
! iggi  - initial lattice vector of k point
! iggq  - lattice vector of k-q point
! iggs  - lattice vector of kpt after symmetry transform to bring k point to calculated k point set
! iggqs - similar to iggs, but for k-q point

kvec = (/0.d0,0.d0,0.d0/)
kvecq = (/0.d0,0.d0,0.d0/)
do ii=1,3
  do jj=1,3
    kvec(ii) = kvec(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*kpt(jj,ikpt)
    kvecq(ii) = kvecq(ii)+ symrel(jj,ii,mod(isymq-1,nsym)+1)*kpt(jj,ikptq)
  enddo
enddo

do ii=1,3
  jmatel(ii)=(0.d0,0.d0)
enddo
ishift=igsymq-igg-igglf -igsymk
  iggi(:)=ipwx(:,iunocc)
  iggs = ishift
  do ii=1,3
    do jj=1,3
      iggs(ii) = iggs(ii)+ symrel(jj,ii,mod(isym-1,nsym)+1)*iggi(jj)
    enddo
  enddo
  iggq = (/0,0,0/)
  do ii=1,3
    do jj=1,3
      iggq(ii) = iggq(ii)+ syminv(jj,ii,mod(isymq-1,nsym)+1)*iggs(jj)
    enddo
    if (iggq(ii).lt.igmn(ii).or.iggq(ii).gt.igmx(ii)) return
  enddo
  igq=igndx(ikptq,iggq(1),iggq(2),iggq(3))
  if (igq.eq.0) return
  icgq=indxkcg(ikptq)+(iocc-1)*npwarr(ikptq)+igq
  xcgq=cg(icgq)
  cmatel=conjg(xcgq);
  iggqs = iggs-igsymq
  iggs = iggs-igsymq+igg+igglf
  do ii=1,3
    dvecl(ii) = kvec(ii)+iggs(ii)+kvecq(ii)+iggqs(ii)
  enddo
  jmatel(1:3)=jmatel(1:3) + cmatel*dvecl(1:3)/2.d0

return
end subroutine mkvpwmatelJ1

!*************************************************************************

subroutine CoreCorrectPol(wmax,nwpt,nqpt,npwndx,napwndx,ipwndx, &
&                         pola,polb, &
&                         vol,pi,hart)
implicit none
integer :: nwpt,nqpt,npwndx,napwndx
integer :: ipwndx(2,napwndx)
double complex :: pola(nwpt,nqpt+9,npwndx),polb(nwpt,nqpt+9,npwndx)
double precision :: vol,hart,wmax,pi
integer :: iwpt,iqpt,iipw,ipw1,ipw2
character(80) :: line
integer :: ios,icmt
integer :: nlvlmax,nlvl,ilvl
parameter (nlvlmax=200)
integer :: nel(nlvlmax)
double precision :: elvl(nlvlmax),dw
double precision :: prefactor,ww,ww4,logarg,term
double precision :: polci,polcr

 open(unit=9,file='core.inp',status='old',iostat=ios)
 if (ios.ne.0) return
 write(6,*) "Reading core corrections"
 nlvl=0
 do 
   read(9,'(a)',iostat=ios) line
   if (ios.ne.0) exit
   icmt = index(line,"#")
   if(icmt.ne.1) then
     nlvl=nlvl+1
     read(line,*) nel(nlvl),elvl(nlvl)
     elvl(nlvl)=elvl(nlvl)/hart
     if (nlvl.gt.nlvlmax) then
       write(6,*) "Warning - number of entries in core.dat exceeds allocated number of core levels"
       write(6,*) "change parameter nlvlmax in subroutine CoreCorrectPol"
       exit
     endif
   endif
 enddo
 close(9)

 write(6,*) "Applying core corrections"
 dw=wmax/dble(nwpt)
 do ilvl=1,nlvl
   prefactor=4.d0*pi*nel(ilvl)*elvl(ilvl)**2/vol
   do iwpt=1,nwpt
     ww=iwpt*dw
!write(6,*) "energy = ",ww*hart
     ww4=ww**4
     if (ww.ge.elvl(ilvl)) then
       polci=prefactor/ww4
     else
       polci=0.d0
     endif
     logarg=((ww-elvl(ilvl))**2*(ww+elvl(ilvl))**2+dw**4)/(elvl(ilvl)**4+dw**4)
     term=elvl(ilvl)**2/(ww**2*(elvl(ilvl)**4+dw**4))
     polcr=prefactor*((0.5d0/ww4)*log(logarg)+term)
     do iipw=1,npwndx
       ipw1=ipwndx(1,iipw)
       ipw2=ipwndx(2,iipw)
       if (ipw1.ne.ipw2) cycle
       do iqpt=1,nqpt
!if (ilvl.eq.1.and.iqpt.eq.1.and.iipw.eq.1) write(77,'(7f10.5)') ww*hart, &
!& dble(pola(iwpt,iqpt,iipw)),dimag(pola(iwpt,iqpt,iipw)),polcr, &
!& dble(polb(iwpt,iqpt,iipw)),dimag(polb(iwpt,iqpt,iipw)),polci
         pola(iwpt,iqpt,iipw)=pola(iwpt,iqpt,iipw)+polcr
         polb(iwpt,iqpt,iipw)=polb(iwpt,iqpt,iipw)+(0.d0,1.d0)*polci
       enddo
       do iqpt=1,3
         pola(iwpt,nqpt+4*iqpt-3,iipw)=pola(iwpt,nqpt+4*iqpt-3,iipw)+polcr
         polb(iwpt,iqpt,iipw)=polb(iwpt,iqpt,iipw)+(0.d0,1.d0)*polci
       enddo
     enddo
   enddo
 enddo

end subroutine CoreCorrectPol

!*************************************************************************

subroutine mkmatelX(iocc,iunocc,ikpt,ikptq,igglf,igg, &
& ncg,ncgq,nkpt,nkptq,npwt,npwtq,igmx,igmn,igndx,igndxq, &
& isym,isymq,symrel,syminv,nsym,nsymq,ihlf,ihlfq,kpt,kptq, &
& lvt,lvtq, &
& cg,cgq,indxkcg,indxkcgq,indxkpw,indxkpwq,npwarr,npwarrq,kg,kgq, &
& cmatel)
implicit none
integer :: iocc,iunocc,ikpt,ikptq,ncg,ncgq,nkpt,nkptq,npwt,npwtq,igg(3)
integer :: isym,isymq,nsym,nsymq,symrel(3,3,nsym),syminv(3,3,nsym)
integer :: lvt(3),lvtq(3)
integer :: igmx(3),igmn(3)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: igndxq(nkptq,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
double complex :: cg(ncg),cgq(ncgq)
integer :: indxkcg(nkpt),indxkpw(nkpt),npwarr(nkpt),kg(3,npwt)
integer :: indxkcgq(nkptq),indxkpwq(nkptq),npwarrq(nkptq),kgq(3,npwtq)
double complex :: cmatel,xcg,xcgq
integer :: iisym,isign,iisymq,isignq,ihlf(nkpt),ihlfq(nkptq)
double precision :: kpt(3,nkpt),kptq(3,nkptq)
integer :: ig,igr,ii,jj,igq,icg,icgq,igo
integer :: iggi(3),iggr(3),iggq(3),iggqr(3),igglf(3)
! iggi  - initial lattice vector of kpt
! iggq  - lattice vector of k-q point

cmatel=(0.d0,0.d0)
l1: do ig=1,npwarr(ikpt)
  icg=indxkcg(ikpt)+(iunocc-1)*npwarr(ikpt)+ig
  xcg=cg(icg)
  iggi(1:3)=kg(1:3,indxkpw(ikpt)+ig)
  do ii=1,3
    iggq(ii)=iggi(ii)-igg(ii)-igglf(ii)
  enddo
  do ii=1,3
    if (iggq(ii).lt.igmn(ii).or.iggq(ii).gt.igmx(ii)) cycle l1
  enddo
  igq=igndxq(ikptq,iggq(1),iggq(2),iggq(3))
  if (igq.eq.0) cycle
  icgq=indxkcgq(ikptq)+(iocc-1)*npwarrq(ikptq)+igq
  xcgq=cgq(icgq)
  cmatel=cmatel+xcg*conjg(xcgq)
!  write(6,'(3i3,3x,2f10.6,3x,3i3,3x,2f10.6)') &
!&               kg(1:3,indxkpw(ikpt)+ig),xcg, &
!&               iggq(1:3),xcgq
!  write(6,1000) ig,iggi(1:3),iggq(1:3),xcg,xcgq,cmatel
!1000 format(i4,4x,'a',3i3,4x,'b',3i3,3(1x,2f7.3))
enddo l1

return
end subroutine mkmatelX

!*************************************************************************

subroutine mkmatelP(pi,xred,natom,iocc,iunocc,ikpt,ikptq,qq,igg,ngkpt, &
&  pwmatel,tpwmatel,projwf,nlmn,nkpt,nkptq,nbcore,ncband,amatel)
use geometry
implicit none
integer :: iocc,iunocc,ikpt,ikptq,igg(3),ngkpt(3),nlmn,nbcore,ncband,nkpt,nkptq
double precision :: qq(3),xred(3,natom),pi
double complex :: pwmatel(nlmn,nlmn),tpwmatel(nlmn,nlmn)
double complex :: projwf(natom,nlmn,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband)
double complex :: amatel,phase
integer :: natom,iorb,jorb,iat

!write(6,*) iocc,iunocc,ikpt,ikptq
amatel=(0.d0,0.d0)
do iat=1,natom
  phase=exp((0.d0,2.d0)*pi*dot_product(-qq-igg,xred(:,iat)))
!write(6,'(3f10.5)') xred(:,iat)
!write(6,'(3f10.5)') -qq-igg
!write(6,'(3f10.5)') dot_product(-qq-igg,xred(:,iat))
!write(6,'(4f16.8)') phase
  do iorb=1,nlmn
  do jorb=1,nlmn
    amatel=amatel+projwf(iat,iorb,ikpt,iunocc) &
&          *conjg(projwf(iat,jorb,ikptq,iocc)) &
&          *(pwmatel(iorb,jorb)-tpwmatel(iorb,jorb))*phase
!write(6,'(2i3,4f16.8)') iorb,jorb,projwf(iat,iorb,ikpt,iunocc),conjg(projwf(iat,jorb,ikptq,iocc))
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb),tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') pwmatel(iorb,jorb)-tpwmatel(iorb,jorb)
!write(6,'(6x,4f16.8)') projwf(iat,iorb,ikpt,iunocc)*conjg(projwf(iat,jorb,ikptq,iocc))*(pwmatel(iorb,jorb)-tpwmatel(iorb,jorb))*phase
!write(6,'(6x,4f16.8)') amatel
!write(6,*)
  enddo
  enddo
enddo

return
end subroutine mkmatelP
