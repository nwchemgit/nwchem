subroutine rdkss(fnamekss1,fnamekss2,ipaw,npsp,nkpt,ntypat)
use kssvars
implicit none
character(80) :: fnamekss1,fnamekss2
integer :: ipaw,nkpt,ntypat
character*80 :: title1,title2
double precision :: hart,ryd
parameter (ryd=13.6056981,hart=2.d0*ryd)
integer :: ik,ib,ips,il,ipwt,ish,is,jtypat,ieof,ikn,ikwfm,jps,kat,in,ioffmx, &
& isym,ioff,jloop,npsp,iloop,i,j,k,ig,iban,ikpt,iat,ipsp,ibwfm, &
& ict,ibndt
integer :: ii,npw2,nspinor2,nbnd2,ibnd,isppol,ios,lv,lend
integer, allocatable :: occ2(:)

 open(unit=10,file=fnamekss1,status='old',form='unformatted')
 read(10) !codvsn,headform,fform

 read(10) !bantot,date,intxc,ixc,natom,ngfft(1:3),&
!& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
!& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear

! allocate (istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),so_typat(ntypat),&
!& symafm(nsym),symrel(3,3,nsym),typat(natom))
! allocate (kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat), &
!& xred(3,natom))
 read(10) !istwfk(1:nkpt),nband(1:nkpt*nsppol),&
!& npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym),symrel(1:3,1:3,1:nsym), &
!& typat(1:natom),kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym), &
!& znucltypat(1:ntypat)

 do ipsp=1,npsp
! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, except if alchemical pseudo-atoms)
  read(10) !title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc
 enddo

!(final record: residm, coordinates, total energy, Fermi energy)
 read(10) !residm,xred(1:3,1:natom),etotal,fermie

 if (ipaw.ne.0) read(10) !idummy
 if (ipaw.ne.0) read(10) !idummy

 read(10) title1(1:80)
 read(10) title2(1:80)
 read(10) nsym2,nbandkss,npwkss,ishm,mpsang
! allocate (symrel2(3,3,nsym2),tnons2(3,nsym2),gbig(3,npwkss),shlim(ishm))
 allocate (gbig(3,npwkss))
 read(10) !(((symrel2(i,j,k),i=1,3),j=1,3),k=1,nsym2)
 read(10) !((tnons2(i,k),i=1,3),k=1,nsym2)
 read(10) ((gbig(i,ig),i=1,3),ig=1,npwkss)
 read(10) !(shlim(in),in=1,ishm)

! allocate (vkbsign(ntypat,mpsang),vkb(npwkss,ntypat,mpsang,nkpt), &
!& vkbd(npwkss,ntypat,mpsang,nkpt))
 allocate (en(nkpt,nbandkss),wfg(npwkss,nbandkss,nkpt))
 read(10) !((vkbsign(is,il),il=1,mpsang),is=1,ntypat)
 do ik=1,nkpt
   do is=1,ntypat
     do il=1,mpsang
       read(10) !(vkb(ig,is,il,ik),ig=1,npwkss)
       read(10) !(vkbd(ig,is,il,ik),ig=1,npwkss)
     end do
   end do
   read(10) (en(ik,ib),ib=1,nbandkss)
   do ib=1,nbandkss
     read(10) (wfg(ig,ib,ik),ig=1,npwkss)
   end do
 end do

! deallocate (symrel2,tnons2,shlim,vkbsign,vkb,vkbd)

 close(10)

 return
end subroutine rdkss
