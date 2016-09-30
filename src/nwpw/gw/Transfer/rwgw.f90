!**************************************************************************

module kssvars
  character*6 :: codvsn
  integer :: headform,fform
  character*132 :: title
  double precision :: znuclpsp,zionpsp
  integer :: pspso,pspdat,pspcod,pspxc,lmax,lloc
  integer :: nsym2,nbandksseff,npwkss
  integer,allocatable :: symrel2(:,:,:),gbig(:,:)
  double precision,allocatable :: en(:,:)
  double complex, allocatable :: wfg(:,:,:)
  double precision,allocatable :: wtk(:)
  double precision, allocatable :: kpt(:,:),occ(:)
  integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,&
& nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
  double precision :: acell(3),ecut,ecutdg,ecutsm,ecut_eff,qptn(3), &
& rprimd(3,3),stmbias,tphysel,tsmear
  integer,allocatable :: istwfk(:),nband(:),npwarr(:),so_typat(:),&
& symafm(:),symrel(:,:,:),typat(:)
  double precision, allocatable :: tnons(:,:),znucltypat(:)
  double precision :: residm,etotal,fermie
  double precision,allocatable :: xred(:,:)
  integer,allocatable :: shlim(:)
  real,allocatable :: tnons2(:,:)
  double precision,allocatable :: vkbsign(:,:),vkb(:,:,:,:),vkbd(:,:,:,:)
  integer :: ishm,mpsang
end module kssvars

!**************************************************************************

program rwgw
use kssvars
implicit none
character*80 :: line,lineold
integer :: ik,i1,ii,ngwb,ios,ird,iband,ib,ix,iy,iz
double precision :: e0,vxc,sigx,sigc,z,dsig,sige,engy
double precision, allocatable :: xk(:,:),gwc(:,:),scz(:),gw(:,:,:,:)
integer :: ncband,nkc,ist,ind,nbocc,nwpt,ngcut,nwritegw
double precision :: wmax,qq(3)
character(50) :: fname1,fname2
integer, allocatable :: iptind(:,:,:)

nkc=10
ist=1-nkc/2
ind=nkc/2
ngwb=20
nwritegw=8

fname1='../GWEo_DS2_KSS'
fname2='../GWE.out'
call rdkss(fname1,fname2)

allocate (iptind(ist:ind,ist:ind,ist:ind))
call mkiptind((/0.d0,0.d0,0.d0/),nsym,nkpt,kpt,nkc,ist,ind,symrel,iptind)

open(unit=10,file='../GWE.out',status='old')
write(6,*) 'nkpt = ',nkpt
allocate(xk(3,nkpt),gwc(ngwb,nkpt),scz(nkpt),gw(ngwb,ist:ind,ist:ind,ist:ind))
ik=0
read(10,'(a)',iostat=ios) line
do
  lineold=line
  read(10,'(a)',iostat=ios) line
  ird=index(line,'Band     E0 <VxcLDA>   SigX SigC(E0)')
  if (ird.ne.0) then
    ik=ik+1
    read(lineold(5:),*) xk(1:3,ik)
    do iband=1,ngwb
      read(10,*) i1,e0,vxc,sigx,sigc,z,dsig,sige,gwc(iband,ik),engy
    enddo
    do ii=1,3
      read(10,*)
    enddo
    read(10,'(a)',iostat=ios) line
    ird=index(line,'DeltaE^GW_gap')
    if (ird.ne.0) then
      read(line(ird+14:),*) scz(ik)
    endif
    if (ik.ge.nkpt) exit
  endif  
enddo

!do ik=1,nkpt
!!  write(6,'(i3,20f7.3)') ik,gwc(1:ngwb,ik)
!  write(6,'(i3,10f7.3)') ik,gwc(1:min(10,ngwb),ik)
!!  write(15,'(i3,10f7.3)') ik,gwc(1:min(10,ngwb),ik)
!!  write(6,'(i3,f7.3)') ik,scz(ik)
!enddo

open(unit=11,file='Scr0/gw.dat',status='unknown')
write(11,'(5i5)') nkc,ist,ind,nwritegw
do ix=ist,ind
do iy=ist,ind
do iz=ist,ind
  ik=iptind(ix,iy,iz)
!  write(6,'(3i3,3x,i4)') ix,iy,iz,ik
!  write(15,'(3i3,3x,i4)') ix,iy,iz,ik
  write(11,'(3i3,10f7.3)') ix,iy,iz,gwc(1:min(nwritegw,ngwb),ik)
enddo
enddo
enddo

end program rwgw

!**************************************************************************

subroutine rdkss(fname1,fname2)
use kssvars
implicit none
character(50) :: fname1,fname2
character*80 title1,title2,line
double precision :: hart,ryd
parameter (ryd=13.6056981,hart=2.d0*ryd)
integer :: ik,ib,ips,il,ipw,ish,is,jtypat,ieof,ikn,ikwfm,jps,kat,in,ioffmx, &
& isym,ioff,jloop,nln,npsp,iloop,iwrite,i,j,k,ig,iban,ikpt,iat,ipsp,ibwfm
double precision :: wfmax


 open(unit=10,file=fname1,status='old',form='unformatted')
 read(10) codvsn,headform,fform
 read(10) bantot,date,intxc,ixc,natom,ngfft(1:3),&
& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear
 allocate (istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),so_typat(ntypat),&
& symafm(nsym),symrel(3,3,nsym),typat(natom))
 allocate (kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat), &
& xred(3,natom))
 read(10) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
& npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym), &
& symrel(1:3,1:3,1:nsym),typat(1:natom),&
& kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat)
 nln=6

 do ipsp=1,npsp
! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, except if alchemical pseudo-atoms)
  read(10) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc
 enddo
!(final record: residm, coordinates, total energy, Fermi energy)
 read(10) residm,xred(1:3,1:natom),etotal,fermie
 read(10) title1(1:80)
 read(10) title2(1:80)
 read(10) nsym2,nbandksseff,npwkss,ishm,mpsang
 allocate (symrel2(3,3,nsym),tnons2(3,nsym),gbig(3,npwkss),shlim(ishm))
 read(10) (((symrel2(i,j,k),i=1,3),j=1,3),k=1,nsym2)
 read(10) ((tnons2(i,k),i=1,3),k=1,nsym2)
 read(10) ((gbig(i,ig),i=1,3),ig=1,npwkss)
 read(10) (shlim(in),in=1,ishm)

 allocate (vkbsign(ntypat,mpsang),vkb(npwkss,ntypat,mpsang,nkpt), &
& vkbd(npwkss,ntypat,mpsang,nkpt))
 allocate (en(nkpt,nbandksseff),wfg(npwkss,nbandksseff,nkpt))
 read(10) ((vkbsign(is,il),il=1,mpsang),is=1,ntypat)
 do ik=1,nkpt
   do is=1,ntypat
     do il=1,mpsang
       read(10) (vkb(ig,is,il,ik),ig=1,npwkss)
       read(10) (vkbd(ig,is,il,ik),ig=1,npwkss)
     end do
   end do
   read(10) (en(ik,ib),ib=1,nbandksseff)
   do ib=1,nbandksseff
     read(10) (wfg(ig,ib,ik),ig=1,npwkss)
   end do
 end do

 open(unit=12,file=fname2,status='old')
 allocate (wtk(nkpt))
 nln=6
 do
   read(12,'(a)',iostat=ieof) line
   if (ieof.ne.0) exit
   if (line(8:11).eq.'wtk4') then
     read(line(12:80),*) (wtk(ik),ik=1,min(nkpt,nln))
     do ikn=2,(nkpt+nln-1)/nln
       read(12,*) (wtk(ik),ik=nln*(ikn-1)+1,min(nkpt,nln*ikn))
     enddo
     exit
   endif
 enddo
 close(12)
 close(10)

end subroutine rdkss

!***************************************************************************

subroutine rdinput(qq,nbocc,nwpt,wmax,ncband,ngcut,nkc,ist,ind)
implicit none
character(80) :: line
integer :: ncband,nkc,ist,ind,nbocc,nwpt,ngcut,ios
integer :: iwrd,iend
double precision :: wmax,qq(3)

 open(unit=9,file='input.dat',status='old')
 do
   read(9,'(a)',iostat=ios) line
   if (ios.ne.0) exit
   iend=index(line,'!')
   if (iend.eq.0) then
     iend=len(line)
   else
     iend=iend-1
   endif
   iwrd=index(line(:iend),'qq')
   if (iwrd.ne.0) then
     read(line(iwrd+3:iend),*) qq(1),qq(2),qq(3)
   endif
   iwrd=index(line(:iend),'nbocc')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) nbocc
   endif
   iwrd=index(line(:iend),'nwpt')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) nwpt
   endif
   iwrd=index(line(:iend),'wmax')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) wmax
   endif
   iwrd=index(line(:iend),'ncband')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) ncband
   endif
   iwrd=index(line(:iend),'ngcut')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) ngcut
   endif
   iwrd=index(line(:iend),'nkc')
   if (iwrd.ne.0) then
     read(line(iwrd+4:iend),*) nkc
   endif
 enddo
 ist=1-nkc/2
 ind=nkc/2
 close(9)

return
end subroutine rdinput

!**************************************************************************

subroutine mkiptind(shiftk,nsym,nkpt,kpt,nkc,ist,ind,symrel,iptind)
implicit none
integer :: nsym,nkpt,nkc,ist,ind
double precision :: kpt(3,nkpt),shiftk(3)
integer :: iptind(ist:ind,ist:ind,ist:ind)
integer :: ipttr(ist:ind,ist:ind,ist:ind),iptsym(ist:ind,ist:ind,ist:ind)
integer :: idg(3,ist:ind,ist:ind,ist:ind)
integer :: ix,iy,iz,ik,isym1,isym,ic,jc,iwrite
integer :: symrel(3,3,nsym)
double precision :: xkc1(3),xkc(3),ykc(3),diff,tol
 tol=1.d-3
 do ix=ist,ind
! do ix=5,5
   xkc1(1)=(dble(ix)+shiftk(1))/dble(nkc)
   do iy=ist,ind
!   do iy=4,4
     xkc1(2)=(dble(iy)+shiftk(2))/dble(nkc)
     do iz=ist,ind
!     do iz=4,4
       xkc1(3)=(dble(iz)+shiftk(3))/dble(nkc)
       iptind(ix,iy,iz)=0
l1:    do ik=1,nkpt
         do isym1=1,nsym*2
           if (isym1.le.nsym) then
             isym=isym1
             do ic=1,3
               xkc(ic)=xkc1(ic)
               if (nint(xkc(ic)*nkc).gt.nint(0.5*nkc)) xkc(ic)=xkc(ic)-1.0
               if (nint(xkc(ic)*nkc).le.nint(-0.5*nkc)) xkc(ic)=xkc(ic)+1.0
             enddo
           else
             isym=isym1-nsym
             do ic=1,3
               xkc(ic)=-xkc1(ic)
               if (nint(xkc(ic)*nkc).gt.nint(0.5*nkc)) xkc(ic)=xkc(ic)-1.0
               if (nint(xkc(ic)*nkc).le.nint(-0.5*nkc)) xkc(ic)=xkc(ic)+1.0
             enddo
           endif
           do ic=1,3
             ykc(ic)=0.d0
             do jc=1,3
               ykc(ic)=ykc(ic)+symrel(jc,ic,isym)*kpt(jc,ik)
             enddo
             if (nint(ykc(ic)*nkc).gt.nint(0.5*nkc)) then
               ykc(ic)=ykc(ic)-1.0
               idg(ic,ix,iy,iz)=-1
             elseif (nint(ykc(ic)*nkc).le.nint(-0.5*nkc)) then
               ykc(ic)=ykc(ic)+1.0
               idg(ic,ix,iy,iz)=1
             else
               idg(ic,ix,iy,iz)=0
             endif
           enddo
           diff=0.d0
           do ic=1,3
             diff=diff+(ykc(ic)-xkc(ic))**2
           enddo
!           write(15,'(2i3,3f6.2,3x,3f6.2,3x,f10.6)') ik,isym,ykc,xkc,diff
           if (diff.lt.tol) then
             iptind(ix,iy,iz)=ik
             iptsym(ix,iy,iz)=isym
             if (isym1.gt.nsym) then
               ipttr(ix,iy,iz)=2
             else
               ipttr(ix,iy,iz)=1
             endif
             exit l1
           endif
         enddo
       enddo l1
     enddo
   enddo
 enddo

 iwrite=0
 if (iwrite.eq.1) then
   do ix=ist,ind
   do iy=ist,ind
   do iz=ist,ind
!     write(6,'(3i3,3f8.4,i4," e")') ix,iy,iz,(dble(ix)+shiftk(1))/dble(nkc),(dble(iy)+shiftk(2))/dble(nkc),(dble(iz)+shiftk(3))/dble(nkc),iptind(ix,iy,iz)
     write(6,'(3i3,3x,i4,1x,i2)') ix,iy,iz,iptsym(ix,iy,iz),ipttr(ix,iy,iz)
   enddo
   enddo
   enddo
 endif

end subroutine mkiptind

!***************************************************************************



